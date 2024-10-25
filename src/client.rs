use rand::*;
use tfhe::core_crypto::commons::generators::MaskRandomGenerator;
use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::commons::math::random::*;
use tfhe::core_crypto::algorithms::polynomial_algorithms::*;
use crate::params::*;
use crate::seeder::IdSeeder;


pub struct Client {
    // parameters
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u32>,

    // key for encrypting templates and performing inner product
    glwe_sk_ip: GlweSecretKeyOwned<u64>,
    
    // key for blind rotation
    glwe_sk_br: GlweSecretKeyOwned<u32>,

    // global seed for masking templates
    id_seeder: IdSeeder,

    // for generate noise
    noise_seeder: Box<dyn Seeder>,

    // noise distribution
    distribution: Gaussian<f64>,
}

impl Client {
    /// Instantiate a new client with a default `thread_rng`.
    /// `ip_param` is for inner product, and `br_param` is for blind rotation
    pub fn new(ip_param: GlweParameter<u64>, br_param: GlweParameter<u32>) -> Self {
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
            distribution: Gaussian::from_dispersion_parameter(StandardDev(ip_param.std_dev), 0.0),
            id_seeder: IdSeeder::new(((rng.next_u64() as u128) << 64) | (rng.next_u64() as u128)),
            noise_seeder: seeder,
        }
    }

    /// Encrypt a new template with a given ID. The ciphertexts are GGSW ciphertexts but only body needs transferring.
    pub fn encrypt_new_template(&mut self, id: u128, features: &[f64]) -> Vec<GlweBody<Vec<u64>>> {
        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut generator = RandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
            // seeder.as_mut(),
        );
        let masks = self.get_template_masks(id);
        masks
            .into_iter()
            .map(|mask| {
                let mut ct = GlweCiphertext::new(
                    0u64, 
                    self.ip_param.glwe_size, 
                    self.ip_param.polynomial_size, 
                    ciphertext_modulus
                );

                // copy masks to ciphertexts
                ct.get_mut_mask()
                    .as_mut()
                    .iter_mut()
                    .zip(mask.as_polynomial_list().into_container().iter())
                    .for_each(|(tgt, src)| *tgt = *src);

                // TODO: encode
                let pt = PlaintextListOwned::from_container(vec![]);

                self.encrypt_with_existed_mask(&mut ct, pt, &mut generator);

                // get the body
                GlweBody::from_container(
                    ct.get_body().as_polynomial().into_container().iter().copied().collect::<Vec<_>>(), 
                    ciphertext_modulus
                )
            })
            .collect::<Vec<_>>()
    }

    fn encrypt_with_existed_mask<G>(
        &mut self, 
        glwe_ct: &mut GlweCiphertextOwned<u64>, 
        pt: PlaintextListOwned<u64>, 
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
            &pt.as_polynomial(),
        );

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
    fn get_template_masks(&mut self, id: u128) -> Vec<GlweMask<Vec<u64>>> {
        let ciphertext_modulus = CiphertextModulus::new_native();
        let mut generator = RandomGenerator::<ActivatedRandomGenerator>::new(
            self.id_seeder.seed(id),
        );
        
        (0..self.ip_param.glwe_size.0 * self.ip_param.decomposition_level_count.0).map(|_| {
            let mut mask = GlweMask::from_container(
                vec![0u64; glwe_ciphertext_mask_size(GlweDimension(self.ip_param.glwe_size.0 - 1), self.ip_param.polynomial_size)],
                self.ip_param.polynomial_size, 
                ciphertext_modulus,
            );

            generator.fill_slice_with_random_uniform_custom_mod(mask.as_mut(), ciphertext_modulus);

            mask
        }).collect()
    }
}

