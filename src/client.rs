use rand::*;
use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::commons::math::random::*;
use tfhe::core_crypto::algorithms::polynomial_algorithms::*;
use crate::extract_glwe_sample_from_rlwe_ciphertext;
use crate::params::*;
use crate::utils::*;
use crate::seeder::IdSeeder;


pub struct Client {
    // parameters
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u16>,

    // key for encrypting templates and performing inner product
    glwe_sk_ip: GlweSecretKeyOwned<u64>,
    
    // key for blind rotation
    glwe_sk_br: GlweSecretKeyOwned<u16>,

    // global seed for masking templates
    id_seeder: IdSeeder,

    // for generate noise
    noise_seeder: Box<dyn Seeder>,

    // noise distribution
    distribution_ip: Gaussian<f64>,
    distribution_br: Gaussian<f64>,

    threshold: usize,  // always assume threshold to be non-negative in the current implementations
    precision: usize,
}

impl Client {
    /// Instantiate a new client with a default `thread_rng`.
    /// `ip_param` is for inner product, and `br_param` is for blind rotation
    pub fn new(ip_param: GlweParameter<u64>, br_param: GlweParameter<u16>) -> Self {
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

        let mut rng = thread_rng();

        Self {
            ip_param,
            br_param,
            glwe_sk_ip,
            glwe_sk_br,
            distribution_ip: Gaussian::from_dispersion_parameter(StandardDev(ip_param.std_dev), 0.0),
            distribution_br: Gaussian::from_dispersion_parameter(StandardDev(br_param.std_dev), 0.0),
            id_seeder: IdSeeder::new(((rng.next_u64() as u128) << 64) | (rng.next_u64() as u128)),
            noise_seeder: seeder,
            threshold: 0,
            precision: 8,  // [-127, 128]
        }
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

    /// Create Glwe ciphertexts encrypting monomials
    /// The server verifies that the encrypted message is $X^b$ by two tests:
    ///     1. (monomial) ct acts on [1,1,...,1] would yield [1,1,...,1]
    ///     2. (faithful) ct acts on [0,1,...,n-1] would yield $b$ at the constant part
    /// Since $b$ is rounded, or decomposed, we need several Glwe Ciphertexts for wrapping these values
    pub fn transform_mask_to_glwe(&mut self, id: u128, glwe_ct: GlweCiphertextView<u64>)
        -> Vec<(GlweCiphertextOwned<u64>, GlweCiphertextOwned<u64>)>
    {
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        let masks_for_innerprod = self.transform_mask_to_body(id, glwe_ct);
        let br_polydim = 1usize << self.precision;

        let carry = masks_for_innerprod & (1 << (63 - self.precision));
        let masks_for_innerprod_act = u64::wrapping_add(masks_for_innerprod >> (64 - self.precision), carry >> (63 - self.precision));

        // construct the selector
        let mut ct_act = GlweCiphertext::new(
            0u64,
            self.ip_param.glwe_size,
            self.ip_param.polynomial_size,
            CiphertextModulus::new_native()
        );
        let mut ct_act_inv = GlweCiphertext::new(
            0u64,
            self.ip_param.glwe_size,
            self.ip_param.polynomial_size,
            CiphertextModulus::new_native()
        );


        let mut cleartext = vec![0; self.ip_param.polynomial_size.0];

        // monomial
        cleartext[(br_polydim - masks_for_innerprod_act as usize) % br_polydim] = self.ip_param.delta;

        let plaintext = PlaintextList::from_container(cleartext.clone());
        encrypt_glwe_ciphertext(
            &self.glwe_sk_ip,
            &mut ct_act,
            &plaintext,
            self.distribution_ip,
            &mut encryption_generator
        );

        // proof for monomial
        if masks_for_innerprod_act as usize != br_polydim {
            cleartext[(br_polydim - masks_for_innerprod_act as usize) % br_polydim] = 0;
            cleartext[(br_polydim - masks_for_innerprod_act as usize) % br_polydim] =
                u64::wrapping_sub(0, self.ip_param.delta);
        }
        encrypt_glwe_ciphertext(
            &self.glwe_sk_ip,
            &mut ct_act_inv,
            &plaintext,
            self.distribution_ip,
            &mut encryption_generator
        );

        vec![(ct_act, ct_act_inv)]
    }


    pub fn decrypt_lwe(&self, lwe_ct: LweCiphertextOwned<u16>) -> u16 {
        let pt = decrypt_lwe_ciphertext(&self.glwe_sk_br.as_lwe_secret_key(), &lwe_ct);

        // either 0 or 1
        (pt.0 as f32 / self.br_param.delta as f32).round() as u16
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
            .chunks_mut(self.ip_param.decomposition_level_count.0)
            .zip(glwe_ct.as_polynomial_list().iter())
            .fold(0u64, |sum_body, (mut glev_ct, poly)| {
                let new_sum: u64 = glev_ct
                    .iter_mut()
                    .enumerate()
                    .map(|(beta_idx, mut ct)| {
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
                        let (mask, mut body) = ct.get_mut_mask_and_body();
                        let mut body = body.as_mut_polynomial();

                        // TODO: can be accelerated by FFT
                        // GLWE ciphertexts are in the form of (a, <a,s>+m), so decryption is subtracting to the body
                        polynomial_wrapping_sub_multisum_assign(
                            &mut body,
                            &mask.as_polynomial_list(),
                            &self.glwe_sk_ip.as_polynomial_list(),
                        );

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
        // let mut decomposed_masked_ips = vec![];
        // let mut now_prec = self.precision;
        // let prec_mask = ((1 << self.precision) - 1) as u64;
        // while now_prec <= self.mal_precision {
        //     decomposed_masked_ips.push((masks_for_innerprod >> (64 - now_prec)) & prec_mask);
        //     now_prec += self.precision;
        // }

        // let carry = masks_for_innerprod & (1 << (63 - self.precision));
        // let masks_for_innerprod = u64::wrapping_add(masks_for_innerprod >> (64 - self.precision), carry >> (63 - self.precision));
        // #[cfg(feature = "debug")]
        // println!("truncated masked innerprod: {:?}", decomposed_masked_ips);

        // decomposed_masked_ips
        masks_for_innerprod
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

/// Malicious Client
/// The malicious client 
pub struct MalClient {
    mal_param: GlweParameter<u128>,

    // key for encrypting templates and performing inner product
    rlwe_sk: GlweSecretKeyOwned<u128>,
    glwe_sk: GlweSecretKeyOwned<u128>,  // rearranged from rlwe_sk

    // global seed for masking templates
    id_seeder: IdSeeder,

    // for generate noise
    noise_seeder: Box<dyn Seeder>,

    // noise distribution
    distribution: Gaussian<f64>,
    squared_indices: Vec<usize>,

    threshold: usize,  // always assume threshold to be non-negative in the current implementations
    precision: usize,
}

impl MalClient {
    pub fn new(mal_param: GlweParameter<u128>) -> Self {
        let mut seeder = new_seeder();
        
        let mut secret_generator =
            SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());

        let rlwe_sk = GlweSecretKey::generate_new_binary(
            GlweDimension(1),
            PolynomialSize(mal_param.polynomial_size.0 * (mal_param.glwe_size.0 - 1)),
            &mut secret_generator,
        );

        let sk_container = rlwe_sk.clone().into_container();
        let mut glwe_sk = vec![0; sk_container.len()];
        let N = mal_param.polynomial_size.0;
        let n = mal_param.glwe_size.0-1;
        for ni in 0..n {
            for Ni in 0..N {
                glwe_sk[ni * N + Ni] = sk_container[ni + Ni * n];
            }
        }

        let glwe_sk = GlweSecretKey::from_container(glwe_sk, mal_param.polynomial_size);

        let mut rng = thread_rng();

        Self {
            mal_param,
            rlwe_sk,
            glwe_sk,
            distribution: Gaussian::from_dispersion_parameter(StandardDev(mal_param.std_dev), 0.0),
            id_seeder: IdSeeder::new(((rng.next_u64() as u128) << 64) | (rng.next_u64() as u128)),
            noise_seeder: seeder,
            squared_indices: find_valid_squared_indices(N * n),
            threshold: 0,
            precision: ((N * n) as f32).log2().round() as usize - 2,  // log (Nn / 4)
        }
    }

    pub fn new_rlwe_public_key(&mut self) -> GlweCiphertextOwned<u128> {
        let plaintext_list = PlaintextList::from_container(vec![0; self.mal_param.polynomial_size.0 * (self.mal_param.glwe_size.0 - 1)]);
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut ct = GlweCiphertext::new(
            0,
            GlweSize(2),
            PolynomialSize(self.mal_param.polynomial_size.0 * (self.mal_param.glwe_size.0 - 1)),
            ciphertext_modulus
        );

        encrypt_glwe_ciphertext(
            &self.rlwe_sk,
            &mut ct,
            &plaintext_list,
            self.distribution,
            &mut encryption_generator
        );

        ct
    }

    pub fn new_rlwe_relinearizaion_keys(&mut self) -> Vec<GlweCiphertextOwned<u128>> {
        let mut s_squared = Polynomial::from_container(vec![0; self.rlwe_sk.polynomial_size().0]);
        let s_poly = Polynomial::from_container(self.rlwe_sk.as_polynomial_list().into_container());
        polynomial_wrapping_add_mul_assign(&mut s_squared, &s_poly, &s_poly);

        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        (0..self.mal_param.decomposition_level_count.0)
            .map(|beta_idx| {
                let idx = self.mal_param.decomposition_level_count.0 - beta_idx;
                let s_squared_i = s_squared
                    .iter()
                    .map(|&vi| vi * (1 << (idx * self.mal_param.decomposition_base_log.0)))
                    .collect::<Vec<_>>();
                let s_squared_i_poly = Polynomial::from_container(s_squared_i);

                let mut ct = GlweCiphertext::from_container(
                    vec![0; s_squared_i_poly.polynomial_size().0 * 2], 
                    s_squared_i_poly.polynomial_size(), 
                    CiphertextModulus::new_native()
                );
                polynomial_wrapping_add_assign(&mut ct.get_mut_body().as_mut_polynomial(), &s_squared_i_poly);

                encrypt_glwe_ciphertext_assign(
                    &self.rlwe_sk, 
                    &mut ct, 
                    self.distribution, 
                    &mut encryption_generator
                );

                ct
            })
            .collect()
    }

    pub fn encrypt_glwe(&mut self, features: &[f32], scale: f32) -> GlweCiphertextOwned<u128> {
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut glwe = GlweCiphertext::new(
            0u128,
            self.mal_param.glwe_size,
            self.mal_param.polynomial_size,
            ciphertext_modulus,
        );

        let msgs = self.encode(features, scale);
        let plaintext_list = PlaintextList::from_container(msgs);

        encrypt_glwe_ciphertext(
            &self.glwe_sk,
            &mut glwe,
            &plaintext_list,
            self.distribution,
            &mut encryption_generator,
        );

        glwe
    }

    pub fn encrypt_rlwe(&mut self, features: &[f32], scale: f32) -> GlweCiphertextOwned<u128> {
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut rlwe = GlweCiphertext::new(
            0u128,
            GlweSize(2),
            PolynomialSize(self.mal_param.polynomial_size.0 * (self.mal_param.glwe_size.0 - 1)),
            ciphertext_modulus,
        );

        let msgs = self.encode(features, scale);
        let plaintext_list = PlaintextList::from_container(msgs);

        encrypt_glwe_ciphertext(
            &self.rlwe_sk,
            &mut rlwe,
            &plaintext_list,
            self.distribution,
            &mut encryption_generator,
        );

        rlwe
    }

    pub fn decrypt_glwe(&self, ct: &GlweCiphertextOwned<u128>) -> Vec<u128> {
        let mut plaintext_list = PlaintextList::from_container(vec![0u128; self.mal_param.polynomial_size.0]);
        decrypt_glwe_ciphertext(&self.glwe_sk, ct, &mut plaintext_list);

        plaintext_list.as_mut().iter_mut().for_each(|vi| {
            // let decoded = *vi as f128 / (self.mal_param.delta as f128);
            // *vi = decoded.round() as u128;
            *vi = self.div_round(*vi & self.mal_param.ciphertext_mask, self.mal_param.delta) % self.mal_param.plaintext_modulus;
        });

        plaintext_list.into_container()
    }

    pub fn decrypt_rlwe(&self, ct: &GlweCiphertextOwned<u128>) -> Vec<u128> {
        let mut plaintext_list = PlaintextList::from_container(vec![0u128; self.mal_param.polynomial_size.0 * (self.mal_param.glwe_size.0 - 1)]);
        // decrypt_glwe_ciphertext(&self.rlwe_sk, ct, &mut plaintext_list);
        polynomial_wrapping_sub_mul_assign(&mut plaintext_list.as_mut_polynomial(), &ct.as_polynomial_list().get(0), &self.rlwe_sk.as_polynomial_list().get(0));
        // polynomial_wrapping_sub_mul_assign(&mut plaintext_list.as_mut_polynomial(), &ct.as_polynomial_list().get(1), &s_poly);
        polynomial_wrapping_add_assign(&mut plaintext_list.as_mut_polynomial(), &ct.as_polynomial_list().get(1));

        plaintext_list.as_mut().iter_mut().for_each(|vi| {
            *vi = self.div_round(*vi & self.mal_param.ciphertext_mask, self.mal_param.delta) % self.mal_param.plaintext_modulus;
        });

        plaintext_list.into_container()
    }

    pub fn decrypt_rlwe_multiplied(&self, ct: &GlweCiphertextOwned<u128>) -> Vec<u128> {
        let mut plaintext_list = PlaintextList::from_container(vec![0u128; self.mal_param.polynomial_size.0 * (self.mal_param.glwe_size.0 - 1)]);

        let mut s_squared = Polynomial::from_container(vec![0; self.rlwe_sk.polynomial_size().0]);
        let s_poly = Polynomial::from_container(self.rlwe_sk.as_polynomial_list().into_container());
        polynomial_wrapping_add_mul_assign(&mut s_squared, &s_poly, &s_poly);

        polynomial_wrapping_add_mul_assign(&mut plaintext_list.as_mut_polynomial(), &ct.as_polynomial_list().get(0), &s_squared);
        polynomial_wrapping_sub_mul_assign(&mut plaintext_list.as_mut_polynomial(), &ct.as_polynomial_list().get(1), &s_poly);
        polynomial_wrapping_add_assign(&mut plaintext_list.as_mut_polynomial(), &ct.as_polynomial_list().get(2));
        // decrypt_glwe_ciphertext(&rlwe_sk_squared, ct, &mut plaintext_list);

        plaintext_list.as_mut().iter_mut().for_each(|vi| {
            *vi = self.div_round(*vi & self.mal_param.ciphertext_mask, self.mal_param.delta) % self.mal_param.plaintext_modulus;
        });

        plaintext_list.into_container()
    }

    pub fn decrypt_lwe(&self, ct: &LweCiphertextOwned<u128>) -> u128 {
        let pt = decrypt_lwe_ciphertext(&self.rlwe_sk.as_lwe_secret_key(), ct);
        self.div_round(pt.0 & self.mal_param.ciphertext_mask, self.mal_param.delta) % self.mal_param.plaintext_modulus
    }
    
    /// Encrypt a new template with a given ID. The ciphertexts are GGSW ciphertexts but only body needs transferring.
    /// The layout of GGSW ciphertexts is 
    ///     Glev { poly * s_0 }
    ///     Glev { poly * s_1 }
    ///     ...
    ///     Glev { poly }
    pub fn encrypt_new_template_ggsw(&mut self, id: u128, features: &[f32], scale: f32) -> GlweCiphertextListOwned<u128> {
        let ciphertext_modulus = CiphertextModulus::new_native();

        // cleartext for GGSW ciphertexts
        let cleartext: Vec<_> = features
            .iter()
            .map(|&v| {
                let v = (v * scale).round() as i64;
                if v >= 0 {
                    v as u128
                } else {
                    (self.mal_param.plaintext_modulus as i64 + v) as u128
                }
            })
            .collect();
        let clearpoly = Polynomial::from_container(cleartext);

        let mut generator = RandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
        );

        let mut glwe_ct_list = GlweCiphertextList::new(
            0, 
            self.mal_param.glwe_size, 
            self.mal_param.polynomial_size, 
            GlweCiphertextCount(self.mal_param.decomposition_level_count.0 * self.mal_param.glwe_size.0), 
            ciphertext_modulus
        );

        // id as a seed for generating masks of GGSW ciphertexts
        self.fill_with_template_masks(id, &mut glwe_ct_list);

        // computes the body of GGSW ciphertexts
        glwe_ct_list
            .chunks_mut(self.mal_param.decomposition_level_count.0)
            .enumerate()
            .for_each(|(glev_idx, mut chunk)| {
                // poly * s
                let mut clearpoly_s = clearpoly.clone();
                if glev_idx + 1 < self.mal_param.glwe_size.0 {
                    polynomial_wrapping_mul(
                        &mut clearpoly_s, 
                        &clearpoly, 
                        &self.glwe_sk.as_polynomial_list().get(glev_idx)
                    );
                }

                // each chunk is a Glev ciphertext
                chunk.iter_mut()
                    .enumerate()
                    .for_each(|(beta_idx, mut glwe_ct)| {
                        // poly * beta_i * s
                        let idx = self.mal_param.decomposition_level_count.0 - beta_idx;
                        let mut clearpoly_beta_s = clearpoly_s.clone();
                        clearpoly_beta_s.iter_mut()
                            .for_each(|v| {
                                // encode to the MSB
                                *v <<= self.mal_param.decomposition_base_log.0 * idx;

                                // negate if Glev for masks
                                if glev_idx + 1 < self.mal_param.glwe_size.0 {
                                    *v = u128::wrapping_add(!(*v), 1);
                                }
                            });

                        self.encrypt_with_existed_masks(&mut glwe_ct, clearpoly_beta_s, &mut generator);
                    });
            });

        glwe_ct_list
    }

    /// Given a vector of features [f0, f1, ..., f_{N-1}]
    /// construct the rlwe ciphertext encrypting
    /// [f0, f_{N-1}, 0, ..., 0, 
    ///  f1, f_{N-2}, 0, ..., 0,
    ///  f_{N-1}, f0, 0, ..., 0]
    /// where every row contains n numbers
    pub fn encrypt_new_template_rlwe(&mut self, features: &[f32], scale: f32) -> (GlweCiphertextOwned<u128>, u128) {
        let poly_dim = self.mal_param.polynomial_size.0 * (self.mal_param.glwe_size.0 - 1);
        let n = self.mal_param.glwe_size.0 - 1;
        let N = self.mal_param.polynomial_size.0;

        let mut cleartext = vec![0u128; poly_dim];
        let norm = features.iter().map(|&v| {
            let v = (v * scale).round() as u128;
            v * v
        }).sum::<u128>() % self.mal_param.plaintext_modulus;

        for (i, vi) in self.encode(features, scale).into_iter().enumerate() {
            cleartext[n * i] = vi;
            cleartext[n * (N - 1 - i) + 1] = vi;
        }

        let mut rlwe = GlweCiphertext::new(
            0u128, 
            GlweSize(2), 
            PolynomialSize(poly_dim), 
            CiphertextModulus::new_native()
        );

        polynomial_wrapping_add_assign(
            &mut rlwe.get_mut_body().as_mut_polynomial(),
            &Polynomial::from_container(cleartext)
        );

        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        encrypt_glwe_ciphertext_assign(&self.rlwe_sk, &mut rlwe, self.distribution, &mut encryption_generator);

        (rlwe, norm)
    }

    pub fn encrypt_new_lookup_tables(&mut self, id: u128, d_packed: &GlweCiphertextOwned<u128>) -> (GlweCiphertextOwned<u128>, GlweCiphertextOwned<u128>) {
        // construct GGSW first
        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut glwe_ct_list = GlweCiphertextList::new(
            0,
            self.mal_param.glwe_size,
            self.mal_param.polynomial_size,
            GlweCiphertextCount(self.mal_param.decomposition_level_count.0 * self.mal_param.glwe_size.0),
            ciphertext_modulus
        );

        let glwe_ct = extract_glwe_sample_from_rlwe_ciphertext(d_packed, self.mal_param.polynomial_size);

        self.fill_with_template_masks(id, &mut glwe_ct_list);
        let masks_for_innerprod = glwe_ct_list
            .chunks_mut(self.mal_param.decomposition_level_count.0)
            .zip(glwe_ct.as_polynomial_list().iter())
            .fold(0u128, |sum_body, (mut glev_ct, poly)| {
                let new_sum: u128 = glev_ct
                    .iter_mut()
                    .enumerate()
                    .map(|(beta_idx, mut ct)| {
                        let idx = self.mal_param.decomposition_level_count.0 - beta_idx;
                        let poly_i: Vec<_> = poly
                            .into_container()
                            .iter()
                            .map(|&v| {
                                // decompose
                                let v = v >> (self.mal_param.decomposition_base_log.0 * idx);
                                v & ((1u128 << self.mal_param.decomposition_base_log.0) - 1)
                            })
                            .collect();
                        let (mask, mut body) = ct.get_mut_mask_and_body();
                        let mut body = body.as_mut_polynomial();

                        // TODO: can be accelerated by FFT
                        // GLWE ciphertexts are in the form of (a, <a,s>+m), so decryption is subtracting to the body
                        polynomial_wrapping_sub_multisum_assign(
                            &mut body,
                            &mask.as_polynomial_list(),
                            &self.glwe_sk.as_polynomial_list(),
                        );

                        // compute the constant term of polynomial multiplication
                        let body_container = body.into_container();
                        let mut sum: u128 = poly_i[0].wrapping_mul(body_container[0]);
                        for (&vi, &vj) in poly_i[1..].iter().zip(body_container[1..].iter().rev()) {
                            sum = sum.wrapping_sub(vi.wrapping_mul(vj));
                        }
                        sum
                    })
                    .fold(0u128, |acc, x| acc.wrapping_add(x));

                sum_body.wrapping_add(new_sum)
            });

        #[cfg(feature = "debug")]
        println!("masks for inner prod: {}", masks_for_innerprod);
        // let mut decomposed_masked_ips = vec![];
        // let mut now_prec = self.precision;
        // let prec_mask = ((1 << self.precision) - 1) as u64;
        // while now_prec <= self.mal_precision {
        //     decomposed_masked_ips.push((masks_for_innerprod >> (64 - now_prec)) & prec_mask);
        //     now_prec += self.precision;
        // }

        // let carry = masks_for_innerprod & (1 << (63 - self.precision));
        // let masks_for_innerprod = u64::wrapping_add(masks_for_innerprod >> (64 - self.precision), carry >> (63 - self.precision));
        // #[cfg(feature = "debug")]
        // println!("truncated masked innerprod: {:?}", decomposed_masked_ips);

        // decomposed_masked_ips
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        let br_polydim = 1usize << self.precision;
        let masks_plaintext = (masks_for_innerprod / self.mal_param.delta) % self.mal_param.plaintext_modulus;
        let plaintext_precision = (self.mal_param.plaintext_modulus as f32).log2().round() as usize + 1;  // CAUTION, ceiling thus modulus should not be power-of-two
        let masks_for_innerprod_act = masks_plaintext >> (plaintext_precision - self.precision);
        let poly_dim = d_packed.polynomial_size();
        println!("client precision: {}, rest_precision: {}, plaintext precision: {}", self.precision, plaintext_precision - self.precision, (masks_for_innerprod_act as f32).log2());

        // construct the selector
        let mut d_act = GlweCiphertext::new(
            0u128,
            GlweSize(2),
            poly_dim,
            CiphertextModulus::new_native()
        );

        // monomial
        let mut d_act_body = d_act.get_mut_body();
        d_act_body.as_mut()[masks_for_innerprod_act as usize] = self.mal_param.delta;
        d_act_body.as_mut()[poly_dim.0/2 + br_polydim - 1 - masks_for_innerprod_act as usize] = self.mal_param.delta; // counter part
        
        encrypt_glwe_ciphertext_assign(
            &self.rlwe_sk,
            &mut d_act,
            self.distribution,
            &mut encryption_generator
        );
        
        // binaries
        println!("mask ip: {}", masks_for_innerprod & self.mal_param.ciphertext_mask);        
        let rest_precision = plaintext_precision - self.precision;
        let rest_mask = masks_plaintext & ((1 << rest_precision) - 1);
        let mut d_bin = GlweCiphertext::new(
            0u128,
            GlweSize(2),
            poly_dim,
            CiphertextModulus::new_native()
        );
        println!("act mask: {}, rest mask: {}, plaintext: {}", masks_for_innerprod_act, rest_mask, masks_plaintext);

        let mut d_bin_body = d_bin.get_mut_body();
        for i in 0..rest_precision {
            // little endian
            if (rest_mask >> i) & 0x1 == 1 {
                d_bin_body.as_mut()[self.squared_indices[i]] = self.mal_param.delta;
            }
        }

        // correction flag
        let (u_plus, u_minus) = if rest_mask < self.mal_param.plaintext_modulus / 2 {
            (1, 0)
        } else {
            (0, 1)
        };
        d_bin_body.as_mut()[self.squared_indices[rest_precision]] = u_plus * self.mal_param.delta;
        d_bin_body.as_mut()[self.squared_indices[rest_precision + 1]] = u_minus * self.mal_param.delta;

        encrypt_glwe_ciphertext_assign(
            &self.rlwe_sk,
            &mut d_bin,
            self.distribution,
            &mut encryption_generator,
        );

        (d_act, d_bin)
    }

    fn encode(&self, features: &[f32], scale: f32) -> Vec<u128> {
        features
            .iter()
            .map(|&v| {
                let v = (v * scale).round() as i128;
                if v >= 0 {
                    u128::wrapping_mul(v as u128, self.mal_param.delta)
                } else {
                    let t_minus_v = self.mal_param.plaintext_modulus as i128 + v;
                    u128::wrapping_mul(t_minus_v as u128, self.mal_param.delta)
                }
            })
            .collect()
    }

    fn div_round(&self, v: u128, delta: u128) -> u128 {
        let q = v / delta;
        let r = v % delta;
        if r >= delta / 2 {
            q + 1
        } else {
            q
        }
    }

    fn fill_with_template_masks(&mut self, id: u128, glwe_ct_list: &mut GlweCiphertextListOwned<u128>) {
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

    fn encrypt_with_existed_masks<G>(
        &mut self, 
        glwe_ct: &mut GlweCiphertextMutView<u128>, 
        pt_poly: PolynomialOwned<u128>, 
        generator: &mut RandomGenerator<G>
    ) where G : ByteRandomGenerator {
        let ciphertext_modulus = CiphertextModulus::new_native();
        let (mask, mut body) = glwe_ct.get_mut_mask_and_body();

        generator.fill_slice_with_random_from_distribution_custom_mod(
            body.as_mut(),
            self.distribution,
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
            &self.glwe_sk.as_polynomial_list(),
        );
    }
}