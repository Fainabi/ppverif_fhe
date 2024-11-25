use std::collections::BTreeMap;
use polynomial_algorithms::polynomial_wrapping_mul;
use tfhe::core_crypto::prelude::*;
use crate::params::*;

pub struct Server {
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u8>,
    // mal_param: GlweParameter<u32>,

    // any data structure
    database: BTreeMap<u128, Vec<GlweBody<Vec<u64>>>>,
    database_ggsw: BTreeMap<u128, FourierGgswCiphertextOwned>,

    // for masking the sample extraction to product inner product
    lwe_pk: LwePublicKeyOwned<u32>,

    seeder: Box<dyn Seeder>,

    threshold: usize,
    precision: usize,
}

impl Server {
    pub fn new(ip_param: GlweParameter<u64>, br_param: GlweParameter<u8>, mal_param: GlweParameter<u32>, lwe_pk: LwePublicKeyOwned<u8>) -> Self {
        let lwe_size = lwe_pk.lwe_size();
        let lwe_pk_u32 = LwePublicKeyOwned::from_container(
            lwe_pk.into_container().into_iter().map(|v| v as u32).collect::<Vec<_>>(), 
            lwe_size,
            CiphertextModulus::new_native()
        );

        Self {
            ip_param,
            br_param,
            // mal_param,
            database: BTreeMap::new(),
            database_ggsw: BTreeMap::new(),
            lwe_pk: lwe_pk_u32,
            seeder: new_seeder(),
            threshold: 0,
            precision: 8,
        }
    }

    pub fn set_threshold(&mut self, threshold: f32) {
        self.threshold = (threshold * ((1 << self.precision) as f32)).round() as usize;
    }

    pub fn enroll(&mut self, id: u128, template_bodies: Vec<GlweBody<Vec<u64>>>) {
        self.database.insert(id, template_bodies);
    }

    pub fn enroll_ggsw(&mut self, id: u128, template: GgswCiphertextOwned<u64>) {
        let mut template_fft = FourierGgswCiphertext::new(
            self.ip_param.glwe_size,
            self.ip_param.polynomial_size,
            self.ip_param.decomposition_base_log,
            self.ip_param.decomposition_level_count,
        );
        convert_standard_ggsw_ciphertext_to_fourier(&template, &mut template_fft);

        self.database_ggsw.insert(id, template_fft);
    }

    pub fn verify(&mut self, id: u128, query_ct: GlweCiphertextOwned<u64>, lut_ct: GlweCiphertextListOwned<u8>) -> Option<LweCiphertextOwned<u8>> {
        match self.compute_innerprod_body(id, query_ct.as_view()) {
            None => None,
            Some(innerprod_body) => {
                let ct_idx = innerprod_body as usize / self.br_param.polynomial_size.0;
                let glwe_idx = innerprod_body as usize % self.br_param.polynomial_size.0;
                let mut output_lwe = LweCiphertextOwned::new(
                    0, 
                    LweSize((self.br_param.glwe_size.0 - 1) * self.br_param.polynomial_size.0 + 1), 
                    CiphertextModulus::new_native()
                );
                extract_lwe_sample_from_glwe_ciphertext(&lut_ct.get(ct_idx), &mut output_lwe, MonomialDegree(glwe_idx));

                let zero_mask = self.new_zero_encryption();
                output_lwe.as_mut().iter_mut().zip(zero_mask.into_container().into_iter()).for_each(|(ct_v, z_v)| {
                    *ct_v = u8::wrapping_add(*ct_v, z_v);
                });
                Some(output_lwe)
            }
        }
    }

    /// `query_ct` is for calculating the inner product with masks $v + r$
    /// `act_ct` encrypts the monomial $X^(-r)$
    /// The lookup table contains binaries and `act_ct` can choose for the result on $X^v$
    /// For preventing malicious client sending some ct with backdoors, two arithemetic verification
    /// on the `monomial of act_ct` and the `faithfullness` of $r$
    pub fn verify_mal(&mut self, id: u128, query_ct: GlweCiphertextOwned<u64>, act_ct: GlweCiphertextOwned<u32>) -> Option<GlweCiphertextOwned<u32>> {
        match self.compute_innerprod_body(id, query_ct.as_view()) {
            None => None,
            Some(innerprod_body) => {
                if self.database_ggsw.get(&id).is_none() {
                    return None;
                }

                // compute the ciphertext honestly
                let ciphertext_modulus = CiphertextModulus::new_native();
                let mut out = GlweCiphertext::new(
                    0,
                    self.ip_param.glwe_size,
                    self.ip_param.polynomial_size,
                    ciphertext_modulus,
                );

                add_external_product_assign(
                    &mut out,
                    &self.database_ggsw.get(&id).unwrap(),
                    &query_ct
                );


                let br_polydim = 1 << self.precision;

                // construct inner prod lookup table
                let ip_lut = (0..self.ip_param.polynomial_size.0)
                    .map(|i| {
                        if i == 0 {
                            0 as u32
                        } else {
                            (self.ip_param.polynomial_size.0 - i) as u32
                        }
                    })
                    .collect::<Vec<_>>();
                let ip_lut_poly = PolynomialOwned::from_container(ip_lut);

                // construct comparison lookup table
                let mut lut_cleartext = vec![0; self.ip_param.polynomial_size.0];
                for idx in self.threshold..(br_polydim/2) {
                    lut_cleartext[(idx + br_polydim - innerprod_body as usize) % br_polydim] = 1;
                }

                // act on the inner prod lookup table
                // TODO: directly compute the LWE
                let mut act_on_ip_lut = act_ct.clone();
                act_on_ip_lut
                    .as_mut_polynomial_list()
                    .iter_mut()
                    .zip(act_ct.as_polynomial_list().iter())
                    .for_each(|(mut p_res, p_act)| {
                        polynomial_wrapping_mul(&mut p_res, &p_act, &ip_lut_poly);
                    });

                let lwe_size = LweSize((act_on_ip_lut.glwe_size().0 - 1) * act_on_ip_lut.polynomial_size().0 + 1);
                let mut acted_lwe_ip_lut = LweCiphertext::new(0, lwe_size, act_on_ip_lut.ciphertext_modulus());
                extract_lwe_sample_from_glwe_ciphertext(&act_on_ip_lut, &mut acted_lwe_ip_lut, MonomialDegree(0));



                // verifying monomial

                // act on comparison lookup table


                None
            }
        }

    }


    fn compute_innerprod_body(&mut self, id: u128, query_ct: GlweCiphertextView<u64>) -> Option<u64> {
         match self.database.get(&id) {
            None => None,
            Some(template_ggsw_bodies) => {
                let innerprod_body: u64 = template_ggsw_bodies
                    .chunks(self.ip_param.decomposition_level_count.0)
                    .zip(query_ct.as_polynomial_list().iter())
                    .flat_map(|(template_glev_bodies, query_poly)| {
                        let query_container = query_poly.into_container();

                        template_glev_bodies
                            .iter()
                            .enumerate()
                            .map(|(beta_idx, template_body)| {
                                let idx = self.ip_param.decomposition_level_count.0 - beta_idx;
                                let template_container = template_body.as_polynomial().into_container();

                                let mut sum = 0;
                                for (poly_i, (&vi, &wi)) in template_container
                                    .iter()
                                    .zip(query_container.iter().take(1).chain(query_container[1..].iter().rev()))
                                    .enumerate()
                                {
                                    let wi = wi >> (self.ip_param.decomposition_base_log.0 * idx);
                                    let wi = wi & ((1u64 << self.ip_param.decomposition_base_log.0) - 1);
                                    if poly_i == 0 {
                                        sum = sum.wrapping_add(vi.wrapping_mul(wi));
                                    } else {
                                        sum = sum.wrapping_sub(vi.wrapping_mul(wi));
                                    }
                                }

                                sum
                            })
                    })
                    .fold(0u64, |acc, x| acc.wrapping_add(x));

                #[cfg(feature = "debug")]
                println!("inner prod body: {}", innerprod_body);
                let carry = innerprod_body & (1 << (63 - self.precision));
                let innerprod_body = u64::wrapping_add(innerprod_body >> (64 - self.precision), carry >> (63 - self.precision));
                #[cfg(feature = "debug")]
                println!("truncated inner prod body: {}", innerprod_body);

                Some(innerprod_body)
            }
        }

    }

    fn new_zero_encryption(&mut self) -> LweCiphertextOwned<u8> {
        let mut lwe_ct = LweCiphertextOwned::new(
            0, 
            LweSize((self.br_param.glwe_size.0 - 1) * self.br_param.polynomial_size.0 + 1), 
            CiphertextModulus::new_native()
        );
        let mut generator = SecretRandomGenerator::<ActivatedRandomGenerator>::new(
            self.seeder.seed(),
        );

        // pk encryption for u8 is not currently supported in tfhe-rs
        // Use u32 for encryption first, and then take the least significant bits for equal results
        encrypt_lwe_ciphertext_with_public_key(
            &self.lwe_pk, 
            &mut lwe_ct, 
            Plaintext(0u32), 
            &mut generator
        );
        LweCiphertextOwned::from_container(
            lwe_ct.into_container().into_iter().map(|v| (v & 0xff) as u8).collect::<Vec<_>>(), 
            CiphertextModulus::new_native()
        )
    }
}
