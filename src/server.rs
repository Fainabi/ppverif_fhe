use std::collections::BTreeMap;
use tfhe::core_crypto::prelude::*;
use crate::params::*;

pub struct Server {
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u8>,
    mal_param: GlweParameter<u32>,

    // any data structure
    database: BTreeMap<u128, Vec<GlweBody<Vec<u64>>>>,

    // for masking the sample extraction to product inner product
    lwe_pk: LwePublicKeyOwned<u32>,

    seeder: Box<dyn Seeder>,

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
            mal_param,
            database: BTreeMap::new(),
            lwe_pk: lwe_pk_u32,
            seeder: new_seeder(),
            precision: 8,
        }
    }

    pub fn enroll(&mut self, id: u128, template_bodies: Vec<GlweBody<Vec<u64>>>) {
        self.database.insert(id, template_bodies);
    }

    pub fn verify(&mut self, id: u128, query_ct: GlweCiphertextOwned<u64>, lut_ct: GlweCiphertextListOwned<u8>) -> Option<LweCiphertextOwned<u8>> {
        match self.compute_innerprod_body(id, query_ct) {
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

    pub fn verify_mal(&mut self, id: u128, query_ct: GlweCiphertextOwned<u64>, act_ct: GlweCiphertextOwned<u32>) -> Option<GlweCiphertextOwned<u32>> {
        match self.compute_innerprod_body(id, query_ct) {
            None => None,
            Some(innerprod_body) => {
                // construct inner prod lookup table
                let mut ip_lut = (0..self.mal_param.polynomial_size.0)
                    .map(|i| i as u32)
                    .collect::<Vec<_>>();

                // construct comparison lookup table
                let mut lut = vec![0; self.mal_param.polynomial_size.0];

                // act on the inner prod lookup table

                // verifying monomial

                // act on comparison lookup table


                None
            }
        }

    }

    fn compute_innerprod_body(&mut self, id: u128, query_ct: GlweCiphertextOwned<u64>) -> Option<u64> {
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
