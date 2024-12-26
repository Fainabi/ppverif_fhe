use std::collections::BTreeMap;
use polynomial_algorithms::{polynomial_wrapping_add_assign, polynomial_wrapping_add_mul_assign, polynomial_wrapping_mul};
use tfhe::core_crypto::{commons::ciphertext_modulus::CiphertextModulusKind, prelude::*};
use rand::{rngs::ThreadRng, thread_rng, Rng, RngCore};
use crate::{params::*, rlwe::{self, *}, utils::*};

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
    rlwe_pk: GlweCiphertextOwned<u128>,
    rlwe_rlk: Vec<GlweCiphertextOwned<u128>>,

    seeder: Box<dyn Seeder>,
    rng: ThreadRng,

    squared_indices: Vec<usize>,

    precision: usize,
}

impl MalServer {
    pub fn new(mal_param: GlweParameter<u128>, rlwe_pk: GlweCiphertextOwned<u128>, rlwe_rlk: Vec<GlweCiphertextOwned<u128>>) -> Self {
        let poly_dim = rlwe_pk.polynomial_size().0;
        Self {
            mal_param,
            database: BTreeMap::new(),
            rlwe_pk,
            rlwe_rlk,
            seeder: new_seeder(),
            rng: thread_rng(),
            squared_indices: find_valid_squared_indices(poly_dim),
            precision: (poly_dim as f32).log2().round() as usize - 2,
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

    pub fn relinearize(&self, ct: &GlweCiphertextOwned<u128>) -> GlweCiphertextOwned<u128> {
        let mut out = GlweCiphertext::new(0, GlweSize(2), ct.polynomial_size(), CiphertextModulus::new_native());
        let ct_poly = ct.as_polynomial_list();
        let mut out_poly = out.as_mut_polynomial_list();

        polynomial_wrapping_add_assign(&mut out_poly.get_mut(0), &ct_poly.get(1));
        polynomial_wrapping_add_assign(&mut out_poly.get_mut(1), &ct_poly.get(2));

        (0..self.mal_param.decomposition_level_count.0)
            .zip(self.rlwe_rlk.iter())
            .for_each(|(beta_idx, rlk)| {
                let idx = self.mal_param.decomposition_level_count.0 - beta_idx;
                let decomposed = ct_poly
                    .get(0)
                    .into_container()
                    .iter()
                    .map(|&vi| {
                        let v_mod = vi & self.mal_param.ciphertext_mask;
                        (v_mod >> (idx * self.mal_param.decomposition_base_log.0)) & ((1 << self.mal_param.decomposition_base_log.0) - 1)
                    })
                    .collect::<Vec<_>>();

                let dcp_poly = Polynomial::from_container(decomposed);
                out_poly.iter_mut().zip(rlk.as_polynomial_list().iter()).for_each(|(mut out_i, rlk_i)| {
                    polynomial_wrapping_add_mul_assign(&mut out_i, &rlk_i, &dcp_poly);
                });
            });

        out
    }

    pub fn construct_constraints(
        &mut self, 
        id: u128,
        norm: u128,
        d_packed: GlweCiphertextOwned<u128>, 
        d_act: GlweCiphertextOwned<u128>, 
        d_bin: GlweCiphertextOwned<u128>
    ) -> Option<GlweCiphertextOwned<u128>> {
        let glwe = extract_glwe_sample_from_rlwe_ciphertext(&d_packed, self.mal_param.polynomial_size);
        match self.external_product(id, &glwe) {
            None => None,
            Some(c_ip) => {
                // let d_packed_squared = rlwe_multiplication_u96(&d_packed, &d_packed, t)

                Some(c_ip)
            }
        }
    }

    /// Construct three constraints:
    ///     - d_packed[j] in LWE(0), at the zero coeffs
    ///     - d_packed[i] - d_packed[nN - n - i + 1] in LWE(0), loading templates
    ///     - (d_packed * d_packed)[nN - n + 1] in LWE(norm), for valid templates
    pub fn constraint_norm(&mut self, out_lwe: &mut LweCiphertextOwned<u128>, d_packed: &GlweCiphertextOwned<u128>, norm: u128) {
        let n = self.mal_param.glwe_size.0 - 1;
        let N = self.mal_param.polynomial_size.0;
        let poly_dim = d_packed.polynomial_size();

        let d_packed_squared = rlwe_multiplication_u96(d_packed, d_packed, self.mal_param.plaintext_modulus);
        let d_packed_relin = self.relinearize(&d_packed_squared);

        let mut tmp_rlwe = d_packed.clone();
        let mut tmp_lwe = LweCiphertext::new(0u128, LweSize(d_packed.polynomial_size().0 + 1), CiphertextModulus::new_native());
        let mut tmp_lwe2 = tmp_lwe.clone();

        // zeros coeffs
        let mut zero_cleartext = vec![0u128; poly_dim.0];
        for (i, vi) in zero_cleartext.iter_mut().enumerate() {
            let j = poly_dim.0 - i;
            if i != 0 && j % n != 0 && j % n != 1 {
                *vi = self.next_rand_ZZ_t_ast();
            }
        }
        let zero_poly = Polynomial::from_container(zero_cleartext);
        for (mut tmp_i, poly_i) in tmp_rlwe.as_mut_polynomial_list().iter_mut().zip(d_packed.as_polynomial_list().iter()) {
            polynomial_wrapping_mul(&mut tmp_i, &poly_i, &zero_poly);
        }
        extract_lwe_sample_from_glwe_ciphertext(&tmp_rlwe, &mut tmp_lwe, MonomialDegree(0));
        lwe_ciphertext_add_assign(out_lwe, &tmp_lwe);

        // dual coeffs
        for i in 0..N {
            extract_lwe_sample_from_glwe_ciphertext(d_packed, &mut tmp_lwe, MonomialDegree(n * i));
            extract_lwe_sample_from_glwe_ciphertext(d_packed, &mut tmp_lwe2, MonomialDegree(n * (N - 1 - i) + 1));
            lwe_ciphertext_sub_assign(&mut tmp_lwe, &tmp_lwe2);
            lwe_ciphertext_cleartext_mul_assign(&mut tmp_lwe, Cleartext(self.next_rand_ZZ_t_ast()));
            lwe_ciphertext_add_assign(out_lwe, &tmp_lwe);
        }

        // norm to be gamma
        extract_lwe_sample_from_glwe_ciphertext(&d_packed_relin, &mut tmp_lwe, MonomialDegree(n * (N - 1) + 1));
        *tmp_lwe.get_mut_body().data -= 2 * norm * self.mal_param.delta;
        lwe_ciphertext_cleartext_mul_assign(&mut tmp_lwe, Cleartext(self.next_rand_ZZ_t_ast()));
        lwe_ciphertext_add_assign(out_lwe, &tmp_lwe);
    }

    pub fn constraint_monomial(&mut self, out_lwe: &mut LweCiphertextOwned<u128>, d_act: &GlweCiphertextOwned<u128>, d_bin: &GlweCiphertextOwned<u128>) {
        let n = self.mal_param.glwe_size.0 - 1;
        let N = self.mal_param.polynomial_size.0;
        let poly_dim = d_act.polynomial_size();

        let mut tmp_rlwe = d_act.clone();
        let mut tmp_rlwe_mul = GlweCiphertext::new(0u128, GlweSize(3), poly_dim, CiphertextModulus::new_native());
        let mut tmp_lwe = LweCiphertext::new(0u128, LweSize(poly_dim.0 + 1), CiphertextModulus::new_native());
        let mut tmp_lwe2 = tmp_lwe.clone();

        // Monomial 1: zero coeffs
        // CAUTION: the total plaintext precision needs careful calculation
        let mut zero_cleartext = vec![0u128; poly_dim.0];
        for i in poly_dim.0/4..poly_dim.0/2 {
            zero_cleartext[poly_dim.0 - i] = self.next_rand_ZZ_t_ast();
            zero_cleartext[poly_dim.0 - (i + poly_dim.0/2)] = self.next_rand_ZZ_t_ast();
        }
        let zero_poly = Polynomial::from_container(zero_cleartext);

        for (mut tmp_i, poly_i) in tmp_rlwe.as_mut_polynomial_list().iter_mut().zip(d_act.as_polynomial_list().iter()) {
            polynomial_wrapping_mul(&mut tmp_i, &poly_i, &zero_poly);
        }
        extract_lwe_sample_from_glwe_ciphertext(&tmp_rlwe, &mut tmp_lwe, MonomialDegree(0));
        lwe_ciphertext_add_assign(out_lwe, &tmp_lwe);

        // Monomial 2: dual coeffs
        // TODO: optimize by polynomial multiplication
        for i in 0..poly_dim.0/4 {
            extract_lwe_sample_from_glwe_ciphertext(&d_act, &mut tmp_lwe, MonomialDegree(i));
            extract_lwe_sample_from_glwe_ciphertext(&d_act, &mut tmp_lwe2, MonomialDegree(poly_dim.0 * 3 / 4 - 1 - i));
            lwe_ciphertext_sub_assign(&mut tmp_lwe, &tmp_lwe2);
            lwe_ciphertext_cleartext_mul_assign(&mut tmp_lwe, Cleartext(self.next_rand_ZZ_t_ast()));
            lwe_ciphertext_add_assign(out_lwe, &tmp_lwe);
        }

        // Monomial 3 & 4: normalized one
        // TODO: merge relinearization for optimization
        let mut d_act_squared = rlwe_multiplication_u96(d_act, d_act, self.mal_param.plaintext_modulus);        
        let mut zero_cleartext = vec![0u128; poly_dim.0];
        for i in (poly_dim.0/2)..=(poly_dim.0 * 3 / 4 - 1) {
            zero_cleartext[poly_dim.0 - i] = self.next_rand_ZZ_t_ast();
        }
        d_act_squared.get_mut_body().as_mut()[poly_dim.0 * 3 / 4 - 1] -= 2 * self.mal_param.delta;
        let zero_poly = Polynomial::from_container(zero_cleartext);
        for (mut tmp_i, poly_i) in tmp_rlwe_mul.as_mut_polynomial_list().iter_mut().zip(d_act_squared.as_polynomial_list().iter()) {
            polynomial_wrapping_mul(&mut tmp_i, &poly_i, &zero_poly);
        }

        let d_act_relin = self.relinearize(&d_act_squared);
        extract_lwe_sample_from_glwe_ciphertext(&d_act_relin, &mut tmp_lwe, MonomialDegree(0));
        lwe_ciphertext_add_assign(out_lwe, &tmp_lwe);



        let rest_precision = (self.mal_param.plaintext_modulus as f32).log2().round() as usize + 1 - self.precision;  
    }

    fn next_rand_ZZ_t_ast(&mut self) -> u128 {
        self.rng.gen_range(1..self.mal_param.plaintext_modulus)
    }

    fn next_rand_ZZ_t_divided_2(&mut self) -> u128 {
        self.rng.gen_range(0..self.mal_param.plaintext_modulus) & 0xFFFFFFFF_FFFFFFFE
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

