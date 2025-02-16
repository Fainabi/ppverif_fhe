use std::error::Error;
use std::sync::Arc;
use std::collections::BTreeMap;
use rand::{rngs::OsRng, *};
use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::commons::math::random::*;
use tfhe::core_crypto::algorithms::polynomial_algorithms::*;
use fhe::bfv::{ self, BfvParameters, BfvParametersBuilder, Ciphertext as BfvCiphertext, RGSWCiphertext};
use fhe::Result;
use fhe_math::rq::Representation::PowerBasis;
use fhe_traits::*;

use crate::params::*;
use crate::utils::*;
use crate::seeder::IdSeeder;


pub struct Client {
    // parameters
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u16>,
    bfv_param: Arc<BfvParameters>,
    bfv_lvl_params: Vec<Arc<BfvParameters>>,

    // key for encrypting templates and performing inner product
    glwe_sk_ip: GlweSecretKeyOwned<u64>,
    
    // key for blind rotation
    glwe_sk_br: GlweSecretKeyOwned<u16>,

    // key for identification
    bfv_sk: bfv::SecretKey,

    // global seed for masking templates
    id_seeder: IdSeeder,

    // for generate noise
    noise_seeder: Box<dyn Seeder>,

    // noise distribution
    distribution_ip: Gaussian<f64>,
    distribution_br: Gaussian<f64>,

    threshold: usize,  // always assume threshold to be non-negative in the current implementations
    precision: usize,

    database: BTreeMap<u128, Vec<Vec<PlaintextListOwned<u64>>>>,
}

impl Client {
    /// Instantiate a new client with a default `thread_rng`.
    /// `ip_param` is for inner product, and `br_param` is for blind rotation
    pub fn new(ip_param: GlweParameter<u64>, br_param: GlweParameter<u16>, bfv_param: Arc<BfvParameters>) -> Result<Self> {
        let mut seeder = new_seeder();
        
        let mut secret_generator =
            SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());

        let glwe_sk_ip = GlweSecretKey::generate_new_binary(
            GlweDimension(ip_param.glwe_size.0 - 1),
            ip_param.polynomial_size,
            &mut secret_generator,
        );

        let glwe_sk_br = GlweSecretKey::generate_new_binary(
            GlweDimension(br_param.glwe_size.0 - 1),
            br_param.polynomial_size,
            &mut secret_generator,
        );

        let bfv_sk = bfv::SecretKey::random(&bfv_param, &mut OsRng);

        let mut rng = thread_rng();

        let mut bfv_lvl_params = vec![];
        for i in 0..8 {
            let new_param = if i <= 1 {
                bfv_param.clone()
            } else {
                let n = bfv_param.moduli().len();
                bfv::BfvParametersBuilder::new()
                    .set_degree(bfv_param.degree())
                    .set_moduli(&bfv_param.moduli()[..n+1-i])
                    .set_plaintext_modulus(bfv_param.plaintext())
                    .build_arc()?
            };

            bfv_lvl_params.push(new_param);
        }

        Ok(Self {
            ip_param,
            br_param,
            bfv_param,
            bfv_lvl_params,
            glwe_sk_ip,
            glwe_sk_br,
            bfv_sk,
            distribution_ip: Gaussian::from_dispersion_parameter(StandardDev(ip_param.std_dev), 0.0),
            distribution_br: Gaussian::from_dispersion_parameter(StandardDev(br_param.std_dev), 0.0),
            id_seeder: IdSeeder::new(((rng.next_u64() as u128) << 64) | (rng.next_u64() as u128)),
            noise_seeder: seeder,
            threshold: 0,
            precision: 8,  // [-127, 128]
            database: BTreeMap::new(),
        })
    }

    /// generate a list of glwe ciphertext as glwe public keys for blind rotation
    pub fn new_glwe_public_keys_br(&mut self) -> Vec<GlweCiphertextOwned<u16>> {
        let plaintext_list = PlaintextList::from_container(vec![0; self.br_param.polynomial_size.0]);
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        (0..self.br_param.glwe_size.0-1)
            .map(|_| {
                let ciphertext_modulus = CiphertextModulus::new_native();
                let mut ct = GlweCiphertext::new(
                    0,
                    self.br_param.glwe_size,
                    self.br_param.polynomial_size,
                    ciphertext_modulus
                );

                encrypt_glwe_ciphertext(
                    &self.glwe_sk_br,
                    &mut ct,
                    &plaintext_list,
                    self.distribution_br,
                    &mut encryption_generator
                );

                ct
            })
            .collect()
    }

    pub fn new_bfv_relinearizatio_key(&self) -> Result<Vec<bfv::RelinearizationKey>> {
        let mut rlks = vec![];

        for i in 1..self.precision {
            let rlk = bfv::RelinearizationKey::new_leveled(&self.bfv_sk, i - 1, i - 1, &mut thread_rng())?;
            rlks.push(rlk);
        }
        
        Ok(rlks)
    }

    pub fn set_threshold(&mut self, threshold: f32) {
        self.threshold = (threshold * ((1 << self.precision) as f32)).round() as usize;
    }

    /// Encrypt a new template with a given ID. The ciphertexts are GGSW ciphertexts but only body needs transferring.
    /// The layout of GGSW ciphertexts is 
    ///     Glev { poly * s_0 }
    ///     Glev { poly * s_1 }
    ///     ...
    ///     Glev { poly }
    pub fn encrypt_new_template(&mut self, id: u128, features: &[f32], scale: f32) -> Vec<GlweBody<Vec<u64>>> {
        let ciphertext_modulus = CiphertextModulus::new_native();

        // cleartext for GGSW ciphertexts
        let cleartext: Vec<_> = features
            .iter()
            .map(|&v| {
                let v = (v * scale).round() as i64;
                if v >= 0 {
                    v as u64
                } else {
                    (self.ip_param.plaintext_modulus as i64 + v) as u64
                }
            })
            .collect();
        let clearpoly = Polynomial::from_container(cleartext);

        let mut generator = RandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
        );

        let mut glwe_ct_list = GlweCiphertextList::new(
            0, 
            self.ip_param.glwe_size, 
            self.ip_param.polynomial_size, 
            GlweCiphertextCount(self.ip_param.decomposition_level_count.0 * self.ip_param.glwe_size.0), 
            ciphertext_modulus
        );

        // id as a seed for generating masks of GGSW ciphertexts
        self.fill_with_template_masks(id, &mut glwe_ct_list);

        // computes the body of GGSW ciphertexts
        glwe_ct_list
            .chunks_mut(self.ip_param.decomposition_level_count.0)
            .enumerate()
            .flat_map(|(glev_idx, mut chunk)| {
                // poly * s
                let mut clearpoly_s = clearpoly.clone();
                if glev_idx + 1 < self.ip_param.glwe_size.0 {
                    polynomial_wrapping_mul(
                        &mut clearpoly_s, 
                        &clearpoly, 
                        &self.glwe_sk_ip.as_polynomial_list().get(glev_idx)
                    );
                }

                // each chunk is a Glev ciphertext
                chunk.iter_mut()
                    .enumerate()
                    .map(|(beta_idx, mut glwe_ct)| {
                        // poly * beta_i * s
                        let idx = self.ip_param.decomposition_level_count.0 - beta_idx;
                        let mut clearpoly_beta_s = clearpoly_s.clone();
                        clearpoly_beta_s.iter_mut()
                            .for_each(|v| {
                                // encode to the MSB
                                *v <<= self.ip_param.decomposition_base_log.0 * idx;

                                // negate if Glev for masks
                                if glev_idx + 1 < self.ip_param.glwe_size.0 {
                                    *v = u64::wrapping_add(!(*v), 1);
                                }
                            });

                        self.encrypt_with_existed_masks(&mut glwe_ct, clearpoly_beta_s, &mut generator);

                        GlweBody::from_container(
                            glwe_ct.get_body().as_polynomial().into_container().iter().copied().collect::<Vec<_>>(), 
                            ciphertext_modulus
                        )
                    })
                    .collect::<Vec<_>>()
                    .into_iter()    
            })
            .collect()
    }

    pub fn enroll_ggsw_masks(&mut self, id: u128) {
        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut glwe_ct_list = GlweCiphertextList::new(
            0,
            self.ip_param.glwe_size,
            self.ip_param.polynomial_size,
            GlweCiphertextCount(self.ip_param.decomposition_level_count.0 * self.ip_param.glwe_size.0),
            ciphertext_modulus
        );

        self.fill_with_template_masks(id, &mut glwe_ct_list);
        let decrypted_masks = glwe_ct_list
            .chunks(self.ip_param.decomposition_level_count.0)
            .map(|glev_ct| {
                glev_ct
                    .iter()
                    .map(|ct| {
                        let mut pt_list = PlaintextList::new(0, PlaintextCount(self.ip_param.polynomial_size.0));
                        decrypt_glwe_ciphertext(&self.glwe_sk_ip, &ct, &mut pt_list);
                        pt_list
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        self.database.insert(0, decrypted_masks);
    }

    /// Encrypt a new Glwe Ciphertext for inner product
    pub fn encrypt_glwe(&mut self, features: &[f32], scale: f32) -> GlweCiphertextOwned<u64> {
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut glwe = GlweCiphertext::new(
            0u64,
            self.ip_param.glwe_size,
            self.ip_param.polynomial_size,
            ciphertext_modulus,
        );

        let msgs = self.encode(features, scale);
        let plaintext_list = PlaintextList::from_container(msgs);

        encrypt_glwe_ciphertext(
            &self.glwe_sk_ip,
            &mut glwe,
            &plaintext_list,
            self.distribution_ip,
            &mut encryption_generator,
        );

        glwe
    }

    /// Given a Glwe ciphertext and id for query, transform the current masks to blind rotation masks.
    pub fn transform_mask(&mut self, id: u128, glwe_ct: GlweCiphertextView<u64>) -> GlweCiphertextListOwned<u16> {
        let masks_for_innerprod = self.transform_mask_to_body(id, glwe_ct);

        let carry = masks_for_innerprod & (1 << (63 - self.precision));
        let masks_for_innerprod = u64::wrapping_add(masks_for_innerprod >> (64 - self.precision), carry >> (63 - self.precision));

        // lookup table where `mased_innerprod + [theta,N/2)` are assigned with 1, and others with 0
        // Note that when modulo 2, -1 is regarded as one so the rotation is cyclic.
        let br_polydim = 1 << self.precision;
        let mut cleartext = vec![0u16; br_polydim];
        for idx in self.threshold..(br_polydim/2) {
            cleartext[(idx + br_polydim - masks_for_innerprod as usize) % br_polydim] = self.br_param.delta;
        }
        #[cfg(feature = "debug")]
        println!("rotation cleartext {:?}", cleartext);

        let mut cts = GlweCiphertextList::new(
            0,
            self.br_param.glwe_size,
            self.br_param.polynomial_size,
            GlweCiphertextCount((1 << self.precision) / self.br_param.polynomial_size.0),
            CiphertextModulus::new_native()
        );

        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        encrypt_glwe_ciphertext_list(
            &self.glwe_sk_br,
            &mut cts,
            &PlaintextList::from_container(cleartext),
            self.distribution_br,
            &mut encryption_generator
        );

        cts
    }

    /// Given a GLWE ciphertext and id for query, transform the current masks to blind rotation masks. The masks are dercrypted and enrolled in the database.
    pub fn transform_mask_from_database(&mut self, id: u128, glwe_ct: GlweCiphertextView<u64>) -> Option<GlweCiphertextListOwned<u16>> {
        let masks_for_innerprod = self.transform_mask_to_body_from_database(id, glwe_ct);
        if masks_for_innerprod.is_none() {
            return None;
        }
        let masks_for_innerprod = masks_for_innerprod.unwrap();
        #[cfg(feature = "debug")]
        println!("masks for innerprod: {}", masks_for_innerprod);

        let carry = masks_for_innerprod & (1 << (63 - self.precision));
        let masks_for_innerprod = u64::wrapping_add(masks_for_innerprod >> (64 - self.precision), carry >> (63 - self.precision));

        // lookup table where `mased_innerprod + [theta,N/2)` are assigned with 1, and others with 0
        // Note that when modulo 2, -1 is regarded as one so the rotation is cyclic.
        let br_polydim = 1 << self.precision;
        let mut cleartext = vec![0u16; br_polydim];
        for idx in self.threshold..(br_polydim/2) {
            cleartext[(idx + br_polydim - masks_for_innerprod as usize) % br_polydim] = self.br_param.delta;
        }
        #[cfg(feature = "debug")]
        println!("rotation cleartext {:?}", cleartext);

        let mut cts = GlweCiphertextList::new(
            0,
            self.br_param.glwe_size,
            self.br_param.polynomial_size,
            GlweCiphertextCount((1 << self.precision) / self.br_param.polynomial_size.0),
            CiphertextModulus::new_native()
        );

        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        encrypt_glwe_ciphertext_list(
            &self.glwe_sk_br,
            &mut cts,
            &PlaintextList::from_container(cleartext),
            self.distribution_br,
            &mut encryption_generator
        );

        Some(cts)
    }

    pub fn decrypt_lwe(&self, lwe_ct: &LweCiphertextOwned<u16>) -> u16 {
        let pt = decrypt_lwe_ciphertext(&self.glwe_sk_br.as_lwe_secret_key(), lwe_ct);

        // either 0 or 1
        (pt.0 as f32 / self.br_param.delta as f32).round() as u16
    }

    pub fn decyrpt_bfv_ct(&self, bfv_ct: &BfvCiphertext) -> Result<Vec<u64>> {
        let pt = self.bfv_sk.try_decrypt(&bfv_ct)?;
        Vec::<u64>::try_decode(&pt, bfv::Encoding::simd())
    }

    /// For calculated masked inner products, encrypt the bits of these values slot-wisely.
    /// Each row of the ciphertext is the grouped ciphertexts,
    /// each column is the encrypted digits for the values
    pub fn encrypt_rlwe_for_binary_decomposition(&self, vals: &[u64]) -> Result<Vec<Vec<BfvCiphertext>>> {
        let mut rng = thread_rng();

        let n = vals.len();
        let poly_deg = self.bfv_param.degree();
        let group_num = n / poly_deg + if n % poly_deg == 0 { 0 } else { 1 };
        let mut ct_group = vec![];

        for g_idx in 0..group_num {
            let idx_bg = g_idx * poly_deg;
            let idx_ed = ((g_idx + 1) * poly_deg).min(n);
            let mut cts = vec![];

            for i in 0..self.precision {
                let digits: Vec<_> = vals[idx_bg..idx_ed].iter().map(|&v| (v >> i) & 0x1).collect();
                
                let level = if i <= 1 { 0 } else { i - 1 };
                let pt = bfv::Plaintext::try_encode(&digits, bfv::Encoding::simd_at_level(0), &self.bfv_param)?;
                let ct: bfv::Ciphertext = self.bfv_sk.try_encrypt(&pt, &mut rng)?;

                cts.push(ct);
            }
            ct_group.push(cts);
        }

        Ok(ct_group)
    }

    /// Transform the template mask to decrypted body
    fn transform_mask_to_body(&mut self, id: u128, glwe_ct: GlweCiphertextView<u64>) -> u64 {
        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut glwe_ct_list = GlweCiphertextList::new(
            0,
            self.ip_param.glwe_size,
            self.ip_param.polynomial_size,
            GlweCiphertextCount(self.ip_param.decomposition_level_count.0 * self.ip_param.glwe_size.0),
            ciphertext_modulus
        );

        self.fill_with_template_masks(id, &mut glwe_ct_list);
        let masks_for_innerprod = glwe_ct_list
            .chunks(self.ip_param.decomposition_level_count.0)
            .zip(glwe_ct.as_polynomial_list().iter())
            .fold(0u64, |sum_body, (glev_ct, poly)| {
                let new_sum: u64 = glev_ct
                    .iter()
                    .enumerate()
                    .map(|(beta_idx, ct)| {
                        let idx = self.ip_param.decomposition_level_count.0 - beta_idx;
                        let poly_i: Vec<_> = poly
                            .into_container()
                            .iter()
                            .map(|&v| {
                                // decompose
                                let v = v >> (self.ip_param.decomposition_base_log.0 * idx);
                                v & ((1u64 << self.ip_param.decomposition_base_log.0) - 1)
                            })
                            .collect();

                        let mut pt_list = PlaintextList::new(0, PlaintextCount(self.ip_param.polynomial_size.0));
                        decrypt_glwe_ciphertext(&self.glwe_sk_ip, &ct, &mut pt_list);
                        let body = pt_list.as_polynomial();

                        // compute the constant term of polynomial multiplication
                        let body_container = body.into_container();
                        let mut sum: u64 = poly_i[0].wrapping_mul(body_container[0]);
                        for (&vi, &vj) in poly_i[1..].iter().zip(body_container[1..].iter().rev()) {
                            sum = sum.wrapping_sub(vi.wrapping_mul(vj));
                        }
                        sum
                    })
                    .fold(0u64, |acc, x| acc.wrapping_add(x));

                sum_body.wrapping_add(new_sum)
            });

        #[cfg(feature = "debug")]
        println!("masks for inner prod: {}", masks_for_innerprod);

        masks_for_innerprod
    }

    pub fn transform_mask_to_body_from_database(&mut self, id: u128, glwe_ct: GlweCiphertextView<u64>) -> Option<u64> {
        match self.database.get(&id) {
            None => None,
            Some(dec_masks)  => {
                let masks_for_innerprod = dec_masks
                    .iter()
                    .zip(glwe_ct.as_polynomial_list().iter())
                    .fold(0u64, |sum_body, (glev_decrypted, poly)| {
                        let poly_mul = glev_decrypted
                            .iter()
                            .enumerate()
                            .map(|(beta_idx, pt)| {
                                let idx = self.ip_param.decomposition_level_count.0 - beta_idx;
                                let poly_i: Vec<_> = poly
                                    .into_container()
                                    .iter()
                                    .map(|&v| {
                                        // decompose
                                        let v = v >> (self.ip_param.decomposition_base_log.0 * idx);
                                        v & ((1u64 << self.ip_param.decomposition_base_log.0) - 1)
                                    })
                                    .collect();

                                let body = pt.as_polynomial();
                                let body_container = body.into_container();
                                let mut sum: u64 = poly_i[0].wrapping_mul(body_container[0]);
                                for (&vi, &vj) in poly_i[1..].iter().zip(body_container[1..].iter().rev()) {
                                    sum = sum.wrapping_sub(vi.wrapping_mul(vj));
                                }
                                sum
                            })
                            .fold(0u64, |acc, x| acc.wrapping_add(x));

                        sum_body.wrapping_add(poly_mul)
                    });

                Some(masks_for_innerprod)
            }
        }
    }

    fn encode(&self, features: &[f32], scale: f32) -> Vec<u64> {
        features
            .iter()
            .map(|&v| {
                let v = (v * scale).round() as i64;
                if v >= 0 {
                    u64::wrapping_mul(v as u64, self.ip_param.delta)
                } else {
                    let t_minus_v = self.ip_param.plaintext_modulus as i64 + v as i64;
                    u64::wrapping_mul(t_minus_v as u64, self.ip_param.delta)
                }
            })
            .collect()
    }

    fn encrypt_with_existed_masks<G>(
        &mut self, 
        glwe_ct: &mut GlweCiphertextMutView<u64>, 
        pt_poly: PolynomialOwned<u64>, 
        generator: &mut RandomGenerator<G>
    ) where G : ByteRandomGenerator {
        let ciphertext_modulus = CiphertextModulus::new_native();
        let (mask, mut body) = glwe_ct.get_mut_mask_and_body();

        generator.fill_slice_with_random_from_distribution_custom_mod(
            body.as_mut(),
            self.distribution_ip,
            ciphertext_modulus,
        );

        polynomial_wrapping_add_assign(
            &mut body.as_mut_polynomial(),
            &pt_poly,
        );

        // (a, <a,s> + m)
        polynomial_wrapping_add_multisum_assign(
            &mut body.as_mut_polynomial(),
            &mask.as_polynomial_list(),
            &self.glwe_sk_ip.as_polynomial_list(),
        );
    }

    /// Compute the masks for a specific id.
    /// The masks are indexed with their GGSW construction, i.e.
    ///     Masks = GGSW Masks
    ///           = { Glev Masks of a_i }_{a_i}, { Glev Masks of b }
    /// total `GlweSize * DecompositionLevelCount` masks.
    /// The mask for the j-th decomposition of the i-th component of GlweCiphertext is at 
    ///     `i * DecompositionLevelCount + j`
    fn fill_with_template_masks(&mut self, id: u128, glwe_ct_list: &mut GlweCiphertextListOwned<u64>) {
        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut generator = RandomGenerator::<ActivatedRandomGenerator>::new(
            self.id_seeder.seed(id),
        );
        
        glwe_ct_list
            .iter_mut()
            .for_each(|mut ct| {
                generator.fill_slice_with_random_uniform_custom_mod(ct.get_mut_mask().as_mut(), ciphertext_modulus);
            });
    }
}

/*
    =====================================================
                    Malicious Client 
    =====================================================
*/

pub struct MaliciousClient {
    parameters: Arc<BfvParameters>,
    secret_key: bfv::SecretKey,
    database: BTreeMap<u128, RGSWCiphertext>,

    inv_q0_mod_q1: u128,
    inv_q1_mod_q0: u128,

    squared_indices: Vec<usize>,
}

impl MaliciousClient {
    pub fn new() -> Result<Self> {
        let (q0, q1) = (0x3fffffff000001_u64, 0x3fffffff004001_u64);
        let parameters = bfv::BfvParametersBuilder::new()
            .set_degree(4096)
            .set_moduli(&[q0, q1])
            .set_plaintext_modulus((1 << 19) + 21)
            .build_arc()?;

        let secret_key = bfv::SecretKey::random(&parameters, &mut OsRng);
        let squared_indices = find_valid_squared_indices(parameters.degree());

        Ok(Self {
            parameters,
            secret_key,
            database: std::collections::BTreeMap::new(),
            inv_q0_mod_q1: mod_inv(q0 as u128, q1 as u128),
            inv_q1_mod_q0: mod_inv(q1 as u128, q0 as u128),
            squared_indices,
        })
    }

    pub fn new_relin_keys(&self) -> Result<bfv::RelinearizationKey> {
        bfv::RelinearizationKey::new(&self.secret_key, &mut thread_rng())
    }

    pub fn new_public_key(&self) -> bfv::PublicKey {
        bfv::PublicKey::new(&self.secret_key, &mut thread_rng())
    }

    pub fn encrypt_rgsw_ciphertext(&self, features: &[f32], scale: f32) -> Result<(RGSWCiphertext, RGSWCiphertext)> {
        let cleartext = self.encode(features, scale);
        let plaintext = bfv::Plaintext::try_encode(&cleartext, bfv::Encoding::poly(), &self.parameters)?;
        let mut rng = thread_rng();
        let rgsw_ct: RGSWCiphertext = self.secret_key.try_encrypt(&plaintext, &mut rng)?;
        let mut rgsw_masks = rgsw_ct.clone();
        rgsw_masks.zeroize_body();

        Ok((rgsw_ct, rgsw_masks))
    }

    pub fn enroll_rgsw_masks(&mut self, id: u128, rgsw_masks: RGSWCiphertext) {
        self.database.insert(id, rgsw_masks);
    }

    pub fn encrypt_rlwe_ciphertext(&self, features: &[f32], scale: f32) -> Result<(BfvCiphertext, u128)> {
        let mut cleartext = self.encode(features, scale);
        let norm = cleartext.iter().map(|&v| (v * v) as u128).sum::<u128>() % self.parameters.plaintext() as u128;

        let poly_dim = features.len();
        let n = self.parameters.degree() / poly_dim;
        for i in 0..poly_dim {
            cleartext[(poly_dim - 1 - i) * n + 1] = cleartext[i * n];
        }
        let plaintext = bfv::Plaintext::try_encode(&cleartext, bfv::Encoding::poly(), &self.parameters)?;
        // let polypoly = Polynomial::from_container(cleartext.into_iter().map(|v| v as u64).collect::<Vec<_>>());
        // let mut poly2 = polypoly.clone();
        // polynomial_wrapping_mul(&mut poly2, &polypoly, &polypoly);
        // println!("{:?}, norm: {}", poly2.as_ref()[(poly_dim - 1) * n + 1], norm);

        let mut rng = thread_rng();
        let rlwe = self.secret_key.try_encrypt(&plaintext, &mut rng)?;
        Ok((rlwe, norm))
    }

    pub fn encrypt_new_lookup_tables(&self, id: u128, gamma: u128, norm: u128, d_packed: &BfvCiphertext) -> std::result::Result<(BfvCiphertext, BfvCiphertext), Box<dyn Error>> {        
        let q0 = self.parameters.moduli()[0] as u128;
        let q1 = self.parameters.moduli()[1] as u128;
        let ciphertext_modulus = q0 * q1;
        let delta = ciphertext_modulus / self.parameters.plaintext() as u128;
        let mut rng = thread_rng();

        // first calculate the mask with paritial decrypt
        let mask_for_ip = self.compute_mask_for_innerprod(id, d_packed)?;
        let pt_mask = mask_for_ip / delta % self.parameters.plaintext() as u128;

        // construct look-up table
        let precision = ((self.parameters.degree() / 4) as f32).log2().round() as usize;
        let br_polydim = 1usize << precision;
        let plaintext_precision = (self.parameters.plaintext() as f32).log2().round() as usize + 1;  // CAUTION, ceiling thus modulus should not be power-of-two
        let masks_for_innerprod_act = pt_mask >> (plaintext_precision - precision);
        let poly_dim = self.parameters.degree();

        #[cfg(feature = "debug")]
        println!("client precision: {}, rest_precision: {}, plaintext precision: {}", precision, plaintext_precision - precision, (masks_for_innerprod_act as f32).log2());

        // construct the selector
        let mut cleartext_act = vec![0u64; self.parameters.degree()];
        cleartext_act[masks_for_innerprod_act as usize] = 1;
        cleartext_act[poly_dim/2 + br_polydim - 1 - masks_for_innerprod_act as usize] = 1; // counter part
        let pt_act=  bfv::Plaintext::try_encode(&cleartext_act, bfv::Encoding::poly(), &self.parameters)?;

        let d_act = self.secret_key.try_encrypt(&pt_act, &mut rng)?;
        
        // binaries
        // println!("mask ip: {}", masks_for_innerprod & self.mal_param.ciphertext_mask);        
        let rest_precision = plaintext_precision - precision;
        let rest_mask = pt_mask & ((1 << rest_precision) - 1);
        let mut cleartext_bin = vec![0u64; poly_dim];
        for i in 0..rest_precision {
            if (rest_mask >> i) & 0x1 == 1 {
                cleartext_bin[self.squared_indices[i]] = 1;
            }
        }

        // correction flag
        let rest_noise = mask_for_ip % delta;

        #[cfg(feature = "debug")]
        println!("pt_act: {}, pt mask: {}, mask for ip: {}, delta: {}, rest noise: {}, neg noise: {}", masks_for_innerprod_act, pt_mask, mask_for_ip, delta, rest_noise, delta - rest_noise);

        // rest_noise is a_sim % delta, in the paper it shows (-a_sim) % delta
        let (u_plus, u_minus) = if rest_noise < delta / 2 {
            (0, 1)
        } else {
            (1, 0)
        };
        cleartext_bin[self.squared_indices[rest_precision]] = u_plus;
        cleartext_bin[self.squared_indices[rest_precision + 1]] = u_minus;

        // norm
        let sqrt_n = (self.parameters.degree() as f64).sqrt().round() as u128;
        let bound = gamma * sqrt_n + self.parameters.degree() as u128 / 4;
        let gamma_squared = gamma * gamma;
        let diff = (norm + bound - gamma_squared) % self.parameters.plaintext() as u128;
        let num_bits = (bound as f64).log2().ceil() as usize + 1;
        for i in 0..num_bits {
            cleartext_bin[self.squared_indices[rest_precision + 2 + i]] = ((diff >> i) & 0x1) as u64;
        }

        #[cfg(feature = "debug")]
        println!("diff: {}, norm: {}, bound: {}, val to pass: {}, 1 << numbits: {}", diff - bound, norm, bound, diff, 1 << num_bits);

        let pt_bin = bfv::Plaintext::try_encode(&cleartext_bin, bfv::Encoding::poly(), &self.parameters)?;
        let d_bin = self.secret_key.try_encrypt(&pt_bin, &mut rng)?;

        // println!("act mask: {}, rest mask: {}, plaintext: {}", masks_for_innerprod_act, rest_mask, masks_plaintext);
        
        Ok((d_act, d_bin))
    }

    pub fn compute_mask_for_innerprod(&self, id: u128, rlwe: &BfvCiphertext) -> std::result::Result<u128, Box<dyn Error>> {
        let rgsw = self.database.get(&id).ok_or(DatabaseError::KeyNotFound(id))?;
        let rlwe_ip = rgsw * rlwe;
                
        let mut poly = self.secret_key.try_phase(&rlwe_ip)?;
        poly.change_representation(PowerBasis);        
        let coeffs = poly.coefficients();

        let q0 = self.parameters.moduli()[0] as u128;
        let q1 = self.parameters.moduli()[1] as u128;
        let v_q0 = coeffs[(0, 0)] as u128;        
        let v_q1 = coeffs[(1, 0)] as u128;
        let m_0 = ((v_q0 * self.inv_q1_mod_q0) % q0) * q1;
        let m_1 = ((v_q1 * self.inv_q0_mod_q1) % q1) * q0;

        Ok((m_0 + m_1) % (q0 * q1))
    }

    pub fn compute_body_for_innerprod(&self, rgsw: &RGSWCiphertext, rlwe: &BfvCiphertext) -> std::result::Result<u128, Box<dyn Error>> {
        let rlwe_ip = rgsw * rlwe;
        
        let pt = self.secret_key.try_decrypt(&rlwe_ip)?;
        let mut poly = pt.to_poly();
        poly.change_representation(PowerBasis);
        let coeffs = poly.coefficients();

        let q0 = self.parameters.moduli()[0] as u128;
        let q1 = self.parameters.moduli()[1] as u128;
        let v_q0 = coeffs[(0, 0)] as u128;
        let v_q1 = coeffs[(1, 0)] as u128;
        let m_0 = ((v_q0 * self.inv_q1_mod_q0) % q0) * q1;
        let m_1 = ((v_q1 * self.inv_q0_mod_q1) % q1) * q0;

        Ok((m_0 + m_1) % (q0 * q1))
    }

    pub fn decrypt(&self, rlwe: &BfvCiphertext) -> Result<Vec<u64>> {
        let pt = self.secret_key.try_decrypt(rlwe)?;
        let mut poly = pt.to_poly();
        poly.change_representation(PowerBasis);

        #[cfg(feature = "debug")]
        println!("decrypted: {}", poly.coefficients());

        Vec::<u64>::try_decode(&pt, bfv::Encoding::poly())
        // let q0 = self.parameters.moduli()[0] as u128;
        // let q1 = self.parameters.moduli()[1] as u128;
        // let v_q0 = coeffs[(0, 0)] as u128;
        // let v_q1 = coeffs[(1, 0)] as u128;
        // let m_0 = ((v_q0 * self.inv_q1_mod_q0) % q0) * q1;
        // let m_1 = ((v_q1 * self.inv_q0_mod_q1) % q1) * q0;

        // Ok((m_0 + m_1) % (q0 * q1))
    }

    fn encode(&self, features: &[f32], scale: f32) -> Vec<i64> {
        let mut encoded = vec![0; self.parameters.degree()];
        let n = self.parameters.degree() / features.len();
        for (i, &fi) in features.iter().enumerate() {
            encoded[i * n] = (fi * scale).round() as i64;
        }
        
        encoded
    }
}
