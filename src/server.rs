use std::collections::BTreeMap;
use polynomial_algorithms::{polynomial_wrapping_add_mul_assign, polynomial_wrapping_mul};
use tfhe::core_crypto::prelude::*;
use rand::{thread_rng, RngCore};
use crate::{params::*, rlwe::*};

pub struct Server {
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u16>,

    // any data structure
    database: BTreeMap<u128, Vec<GlweBody<Vec<u64>>>>,

    // for masking the sample extraction to product inner product
    glwe_pk: Vec<GlweCiphertextOwned<u16>>,

    seeder: Box<dyn Seeder>,

    precision: usize,
}

impl Server {
    pub fn new(ip_param: GlweParameter<u64>, br_param: GlweParameter<u16>, glwe_pk: Vec<GlweCiphertextOwned<u16>>) -> Self {
        Self {
            ip_param,
            br_param,
            database: BTreeMap::new(),
            glwe_pk,
            seeder: new_seeder(),
            precision: 8,
        }
    }

    pub fn enroll(&mut self, id: u128, template_bodies: Vec<GlweBody<Vec<u64>>>) {
        self.database.insert(id, template_bodies);
    }

    pub fn verify(&mut self, id: u128, query_ct: GlweCiphertextOwned<u64>, lut_ct: GlweCiphertextListOwned<u16>) -> Option<LweCiphertextOwned<u16>> {
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
                let mut zero_lwe = output_lwe.clone();
                let zero_mask = self.new_zero_encryption();

                extract_lwe_sample_from_glwe_ciphertext(&lut_ct.get(ct_idx), &mut output_lwe, MonomialDegree(glwe_idx));
                extract_lwe_sample_from_glwe_ciphertext(&zero_mask, &mut zero_lwe, MonomialDegree(glwe_idx));
                
                lwe_ciphertext_add_assign(&mut output_lwe, &mut zero_lwe);
                // output_lwe.as_mut().iter_mut().zip(zero_lwe.into_container().into_iter()).for_each(|(ct_v, z_v)| {
                //     *ct_v = u16::wrapping_add(*ct_v, z_v);
                // });
                Some(output_lwe)
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

    fn new_zero_encryption(&mut self) -> GlweCiphertextOwned<u16> {
        let mut glwe_ct = GlweCiphertextOwned::new(
            0u16, 
            self.br_param.glwe_size, 
            self.br_param.polynomial_size, 
            CiphertextModulus::new_native()
        );

        let mut rng = thread_rng();

        // glwe_ct += sum_i (r_i * pk_i)
        self.glwe_pk
            .iter()
            .for_each(|pk| {
                let mut r = vec![0; self.br_param.polynomial_size.0];
                // assert that 32 divides r.len() 
                for i in 0..(r.len() / 32) {
                    let ri = rng.next_u32();
                    for j in 0..32 {
                        r[32 * i + j] = ((ri >> j) & 0x1) as u16;
                    }
                }

                let polyr = Polynomial::from_container(r);
                pk.as_polynomial_list()
                    .iter()
                    .zip(glwe_ct.as_mut_polynomial_list().iter_mut())
                    .for_each(|(pi, mut ci)| {
                        polynomial_wrapping_add_mul_assign(&mut ci, &pi, &polyr);
                    });
            });

        // TODO: add noise

        glwe_ct
    }
}


pub struct MalServer {
    mal_param: GlweParameter<u128>,

    // storing entire GGSW ciphertexts
    database: BTreeMap<u128, GlweCiphertextListOwned<u128>>,

    // for masking the sample extraction to product inner product
    glwe_pk: GlweCiphertextOwned<u128>,

    seeder: Box<dyn Seeder>,

    precision: usize,
}

impl MalServer {
    pub fn new(mal_param: GlweParameter<u128>, glwe_pk: GlweCiphertextOwned<u128>) -> Self {
        Self {
            mal_param,
            database: BTreeMap::new(),
            glwe_pk,
            seeder: new_seeder(),
            precision: 8,
        }
    }

    pub fn enroll_ggsw(&mut self, id: u128, template: GlweCiphertextListOwned<u128>) {
        // let mut template_fft = FourierGgswCiphertext::new(
        //     self.mal_param.glwe_size,
        //     self.mal_param.polynomial_size,
        //     self.mal_param.decomposition_base_log,
        //     self.mal_param.decomposition_level_count,
        // );
        // convert_standard_ggsw_ciphertext_to_fourier(&template, &mut template_fft);

        self.database.insert(id, template);
    }

    pub fn external_product(&self, id: u128, glwe: &GlweCiphertextOwned<u128>) -> Option<GlweCiphertextOwned<u128>> {
        match self.database.get(&id) {
            None => None,
            Some(ggsw) => {
                let mut glwe_out = GlweCiphertext::from_container(
                    vec![0u128; glwe.glwe_size().0 * glwe.polynomial_size().0], 
                    self.mal_param.polynomial_size,
                    CiphertextModulus::new_native()
                );
                ggsw.chunks(self.mal_param.decomposition_level_count.0)
                    // .enumerate()
                    .zip(glwe.as_polynomial_list().iter())
                    .for_each(|(chunk, poly)| {
                        chunk.iter()
                            .enumerate()
                            .for_each(|(beta_idx, ct)| {
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
                                let poly_i = Polynomial::from_container(poly_i);

                                // TODO: can be boosted with FFT128, which needs further implementations
                                glwe_out
                                    .as_mut_polynomial_list()
                                    .iter_mut()
                                    .zip(ct.as_polynomial_list().iter())
                                    .for_each(|(mut out_poly, ct_poly)| {
                                        polynomial_wrapping_add_mul_assign(&mut out_poly, &ct_poly, &poly_i);
                                    });
                            });
                    });

                Some(glwe_out)
            }
        }
    }

    // `query_ct` is for calculating the inner product with masks $v + r$
    // `act_ct` encrypts the monomial $X^(-r)$
    // The lookup table contains binaries and `act_ct` can choose for the result on $X^v$
    // For preventing malicious client sending some ct with backdoors, two arithemetic verification
    // on the `monomial of act_ct` and the `faithfullness` of $r$
    // pub fn verify_mal(&mut self, id: u128, query_ct: GlweCiphertextOwned<u64>, act_ct: GlweCiphertextOwned<u32>) -> Option<GlweCiphertextOwned<u32>> {
    //     match self.compute_innerprod_body(id, query_ct.as_view()) {
    //         None => None,
    //         Some(innerprod_body) => {
    //             if self.database_ggsw.get(&id).is_none() {
    //                 return None;
    //             }

    //             // compute the ciphertext honestly
    //             let ciphertext_modulus = CiphertextModulus::new_native();
    //             let mut out = GlweCiphertext::new(
    //                 0,
    //                 self.ip_param.glwe_size,
    //                 self.ip_param.polynomial_size,
    //                 ciphertext_modulus,
    //             );

                // add_external_product_assign(
                //     &mut out,
                //     &self.database_ggsw.get(&id).unwrap(),
                //     &query_ct
                // );


    //             let br_polydim = 1 << self.precision;

    //             // construct inner prod lookup table
    //             let ip_lut = (0..self.ip_param.polynomial_size.0)
    //                 .map(|i| {
    //                     if i == 0 {
    //                         0 as u32
    //                     } else {
    //                         (self.ip_param.polynomial_size.0 - i) as u32
    //                     }
    //                 })
    //                 .collect::<Vec<_>>();
    //             let ip_lut_poly = PolynomialOwned::from_container(ip_lut);

    //             // construct comparison lookup table
    //             let mut lut_cleartext = vec![0; self.ip_param.polynomial_size.0];
    //             for idx in self.threshold..(br_polydim/2) {
    //                 lut_cleartext[(idx + br_polydim - innerprod_body as usize) % br_polydim] = 1;
    //             }

    //             // act on the inner prod lookup table
    //             // TODO: directly compute the LWE
    //             let mut act_on_ip_lut = act_ct.clone();
    //             act_on_ip_lut
    //                 .as_mut_polynomial_list()
    //                 .iter_mut()
    //                 .zip(act_ct.as_polynomial_list().iter())
    //                 .for_each(|(mut p_res, p_act)| {
    //                     polynomial_wrapping_mul(&mut p_res, &p_act, &ip_lut_poly);
    //                 });

    //             let lwe_size = LweSize((act_on_ip_lut.glwe_size().0 - 1) * act_on_ip_lut.polynomial_size().0 + 1);
    //             let mut acted_lwe_ip_lut = LweCiphertext::new(0, lwe_size, act_on_ip_lut.ciphertext_modulus());
    //             extract_lwe_sample_from_glwe_ciphertext(&act_on_ip_lut, &mut acted_lwe_ip_lut, MonomialDegree(0));



    //             // verifying monomial

    //             // act on comparison lookup table


    //             None
    //         }
    //     }

    // }

}

