use rand::*;
use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::commons::math::random::*;
use tfhe::core_crypto::algorithms::polynomial_algorithms::*;
use crate::params::*;
use crate::seeder::IdSeeder;


pub struct Client {
    // parameters
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u8>,
    mal_param: GlweParameter<u32>,

    // key for encrypting templates and performing inner product
    glwe_sk_ip: GlweSecretKeyOwned<u64>,
    
    // key for blind rotation
    glwe_sk_br: GlweSecretKeyOwned<u8>,

    // key for malicious blind rotation
    glwe_sk_mal: GlweSecretKeyOwned<u32>,

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
    pub fn new(ip_param: GlweParameter<u64>, br_param: GlweParameter<u8>, mal_param: GlweParameter<u32>) -> Self {
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

        let glwe_sk_mal = GlweSecretKey::generate_new_binary(
            GlweDimension(mal_param.glwe_size.0 - 1),
            mal_param.polynomial_size,
            &mut secret_generator,
        );


        let mut rng = thread_rng();

        Self {
            ip_param,
            br_param,
            mal_param,
            glwe_sk_ip,
            glwe_sk_br,
            glwe_sk_mal,
            distribution_ip: Gaussian::from_dispersion_parameter(StandardDev(ip_param.std_dev), 0.0),
            distribution_br: Gaussian::from_dispersion_parameter(StandardDev(br_param.std_dev), 0.0),
            id_seeder: IdSeeder::new(((rng.next_u64() as u128) << 64) | (rng.next_u64() as u128)),
            noise_seeder: seeder,
            threshold: 0,
            precision: 8,  // [-128, 127]
        }
    }

    pub fn new_lwe_public_key(&mut self) -> LwePublicKeyOwned<u8> {
        let lwe_size = LweSize((self.br_param.glwe_size.0 - 1) * self.br_param.polynomial_size.0 + 1);
        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut pk = LwePublicKey::new(0u8, lwe_size, LwePublicKeyZeroEncryptionCount(1), ciphertext_modulus);
        let mut generator =  EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            self.noise_seeder.as_mut(),
        );

        generate_lwe_public_key(
            &self.glwe_sk_br.as_lwe_secret_key(), 
            &mut pk, 
            self.distribution_br, 
            &mut generator
        );
        pk
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
    pub fn transform_mask(&mut self, id: u128, glwe_ct: GlweCiphertextView<u64>) -> GlweCiphertextListOwned<u8> {
        let masks_for_innerprod = self.transform_mask_to_body(id, glwe_ct);

        // lookup table where `mased_innerprod + [theta,N/2)` are assigned with 1, and others with 0
        // Note that when modulo 2, -1 is regarded as one so the rotation is cyclic.
        let br_polydim = 1 << self.precision;
        let mut cleartext = vec![0u8; br_polydim];
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

    /// Create a Glwe ciphertext encrypting a monomial
    /// The server verifies that the encrypted message is $X^b$ by two tests:
    ///     1. (monomial) ct acts on [1,1,...,1] would yield [1,1,...,1]
    ///     2. (faithful) ct acts on [0,1,...,n-1] would yield $b$ at the constant part
    pub fn transform_mask_to_Glwe(&mut self, id: u128, glwe_ct: GlweCiphertextView<u64>) -> GlweCiphertextOwned<u32> {
        let mut ct = GlweCiphertext::new(
            0,
            self.mal_param.glwe_size,
            self.mal_param.polynomial_size,
            CiphertextModulus::new_native()
        );

        let masks_for_innerprod = self.transform_mask_to_body(id, glwe_ct);
        let br_polydim = 1usize << self.precision;
        let mut cleartext = vec![0; self.mal_param.polynomial_size.0];

        // monomial
        cleartext[(br_polydim - masks_for_innerprod as usize) % br_polydim] = self.mal_param.delta;

        ct
    }


    pub fn decrypt_lwe(&self, lwe_ct: LweCiphertextOwned<u8>) -> u8 {
        let pt = decrypt_lwe_ciphertext(&self.glwe_sk_br.as_lwe_secret_key(), &lwe_ct);

        // either 0 or 1
        (pt.0 as f32 / self.br_param.delta as f32).round() as u8
    }


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
        // somehow modulus switch, yet noiseless
        let carry = masks_for_innerprod & (1 << (63 - self.precision));
        let masks_for_innerprod = u64::wrapping_add(masks_for_innerprod >> (64 - self.precision), carry >> (63 - self.precision));
        #[cfg(feature = "debug")]
        println!("truncated masked innerprod: {}", masks_for_innerprod);

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

