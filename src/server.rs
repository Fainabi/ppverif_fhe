use std::collections::BTreeMap;
use std::ops::AddAssign;
use std::sync::Arc;
use std::vec;
use rand::{rngs::ThreadRng, thread_rng, Rng, RngCore};
use polynomial_algorithms::polynomial_wrapping_add_mul_assign;
use tfhe::core_crypto::commons::math::random::RandomGenerator;
use tfhe::core_crypto::prelude::*;
use fhe::bfv::{self, Ciphertext as BfvCiphertext, BfvParameters, RGSWCiphertext};
use fhe::Result;
use fhe_traits::*;
use fhe_math::rq::{Representation::PowerBasis, Representation, Poly, self};
use crate::{params::*, utils::*};

pub struct Server {
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u16>,
    bfv_param: Arc<BfvParameters>,

    ctx: Arc<rq::Context>,
    poly_zero: rq::Poly,
    poly_one: rq::Poly,
    bfv_rlk: bfv::RelinearizationKey,

    // any data structure
    database: BTreeMap<u128, Vec<GlweBody<Vec<u64>>>>,

    // for masking the sample extraction to product inner product
    glwe_pk: Vec<GlweCiphertextOwned<u16>>,

    noise_distribution: Gaussian<f64>,
    noise_seeder: Box<dyn Seeder>,

    precision: usize,
}

impl Server {
    pub fn new(
        ip_param: GlweParameter<u64>, 
        br_param: GlweParameter<u16>, 
        bfv_param: Arc<BfvParameters>, 
        glwe_pk: Vec<GlweCiphertextOwned<u16>>, 
        bfv_rlk: bfv::RelinearizationKey,
    ) -> Result<Self> {
        let noise_distribution = Gaussian::from_standard_dev(StandardDev(br_param.std_dev), 0.0);
        let ctx = rq::Context::new_arc(bfv_param.moduli(), bfv_param.degree())?;
        let poly_zero = Poly::zero(&ctx, Representation::Ntt);
        let mut vec_one = vec![0u64; bfv_param.degree()];
        vec_one[0] = 1;

        let poly_one = bfv::Plaintext::try_encode(&vec_one, bfv::Encoding::simd(), &bfv_param)?.to_poly();

        Ok(Self {
            ip_param,
            br_param,
            bfv_param,
            poly_zero,
            poly_one,
            ctx,
            bfv_rlk,
            database: BTreeMap::new(),
            glwe_pk,
            noise_distribution,
            noise_seeder: new_seeder(),
            precision: 8,
        })
    }

    pub fn enroll(&mut self, id: u128, template_bodies: Vec<GlweBody<Vec<u64>>>) {
        self.database.insert(id, template_bodies);
    }

    pub fn verify(&mut self, id: u128, query_ct: &GlweCiphertextOwned<u64>, lut_ct: &GlweCiphertextListOwned<u16>) -> Option<LweCiphertextOwned<u16>> {
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

    pub fn compute_innerprod_body(&mut self, id: u128, query_ct: GlweCiphertextView<u64>) -> Option<u64> {
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

    pub fn batch_fold(&self, cts: Vec<Vec<BfvCiphertext>>, ips: Vec<u64>) -> Result<Vec<Vec<BfvCiphertext>>> {
        let n = ips.len();
        let poly_deg = self.bfv_param.degree();

        let mut ct_group = vec![];
        for (g_idx, ct) in cts.iter().enumerate() {
            let idx_bg = g_idx * poly_deg;
            let idx_ed = ((g_idx + 1) * poly_deg).min(n);

            // construct digits first
            let mut digit_polys = vec![];
            for i in 0..self.precision {
                let digits = ips[idx_bg..idx_ed]
                    .iter()
                    .map(|&v| (v >> (63 - i)) & 0x1)
                    .collect::<Vec<_>>();
                
                let poly = bfv::Plaintext::try_encode(&digits, bfv::Encoding::simd(), &self.bfv_param)?;
                // poly.to_poly()
                // let poly = rq::Poly::zero(&self.ctx, Representation::Ntt);
                
                digit_polys.push(poly.to_poly());
            }

            // now fold the packed ciphertexts and digits for determine sign of (ct + ip)
            // Layout:
            //     cts: [ ct_0^{(0)} ct_0^{(1)} ... ct_0^{7};
            //            ct_1^{(0)} ct_1^{(1)} ... ct_1^{7};
            //            ct_2^{(0)} ct_2^{(1)} ... ct_2^{7};
            //                            ... ]
            //    poly: [ w_0^{(0)} w_0^{(1)} ... w_0^{7};
            //            w_1^{(0)} w_1^{(1)} ... w_1^{7};
            //            w_2^{(0)} w_2^{(1)} ... w_2^{7};
            //                            ... ]
            let mut ct_ans = vec![];
            for (i, (ct_i, w_i)) in ct.iter().zip(digit_polys.iter()).enumerate() {
                
                let mut acc_i = ct_i.clone();
                let neg_double_w = - (w_i + w_i);

                // 2 * w * c
                acc_i[0] *= &neg_double_w;
                acc_i[1] *= &neg_double_w;

                // c - 2 * w * c
                acc_i += ct_i;

                // (w, 0) + c - 2 * w * c
                acc_i[0] += w_i;

                let acc_i = if i + 1 == self.precision {
                    let mut neg_double_acci = -acc_i.clone();
                    neg_double_acci += &-acc_i;
                    neg_double_acci[0] += &self.poly_one;
                    neg_double_acci
                } else {
                    acc_i
                };

                ct_ans.push(acc_i);
            }

            self.fold_mul(&mut ct_ans)?;
            ct_group.push(ct_ans);
        }

        Ok(ct_group)
    }

    pub fn fold_mul(&self, cts: &mut [BfvCiphertext]) -> Result<()> {
        if cts.len() <= 1 {
            return Ok(());
        }

        let len = cts.len();
        let nearest_pow2 = {
            let mut cnt = 0;
            let mut len = len;
            while len > 1 {
                len >>= 1;
                cnt += 1;
            }

            1usize << cnt
        };

        if nearest_pow2 != len {
            self.fold_mul(&mut cts[..len-nearest_pow2])?;
            self.fold_mul(&mut cts[len-nearest_pow2..])?;

            for i in len-nearest_pow2..len {
                cts[i] = &cts[i] * &cts[len-nearest_pow2-1];
                self.bfv_rlk.relinearizes(&mut cts[i])?;
            }
        } else {            
            let mid = nearest_pow2 / 2;
            self.fold_mul(&mut cts[len-mid..])?;
            self.fold_mul(&mut cts[..len-mid])?;

            // then reduce
            for i in len-mid..len {
                cts[i] = &cts[i] * &cts[len-mid-1];
                self.bfv_rlk.relinearizes(&mut cts[i])?;
            }
        }

        // for ct in cts.iter_mut() {
        //     ct.mod_switch_to_next_level()?;
        // }

        Ok(())
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

        let mut generator = RandomGenerator::<ActivatedRandomGenerator>::new(
            self.noise_seeder.seed(),
        );

        glwe_ct.as_mut_polynomial_list()
            .iter_mut()
            .for_each(|mut poly| {
                let mut e = vec![0u16; self.br_param.polynomial_size.0];
                generator.fill_slice_with_random_from_distribution_custom_mod(
                    e.as_mut(),
                    self.noise_distribution,
                    CiphertextModulus::new_native(),
                );

                poly.as_mut().iter_mut().zip(e.into_iter()).for_each(|(pi, ei)| *pi += ei);
            });

        glwe_ct
    }
}

pub struct AntiMalServer {
    parameters: Arc<BfvParameters>,
    ctx: Arc<rq::Context>,

    database: BTreeMap<u128, RGSWCiphertext>,

    public_key: bfv::PublicKey,
    relin_keys: bfv::RelinearizationKey,

    inv_q0_mod_q1: u128,
    inv_q1_mod_q0: u128,
    delta_q0: u128,
    delta_q1: u128,

    threshold: usize,
    squared_indices: Vec<usize>,

    rng: ThreadRng,
    pub gamma: u128,
}

impl AntiMalServer {
    pub fn new(public_key: bfv::PublicKey, relin_keys: bfv::RelinearizationKey) -> Result<Self> {
        let (q0, q1) = (0x3fffffff000001_u64, 0x3fffffff004001_u64);
        let parameters = bfv::BfvParametersBuilder::new()
            .set_degree(4096)
            .set_moduli(&[q0, q1])
            .set_plaintext_modulus((1 << 19) + 21)
            .build_arc()?;

        // construct polynomials for multiplying

        let ctx = rq::Context::new_arc(&[q0, q1], 4096)?;
        let squared_indices = find_valid_squared_indices(parameters.degree());

        let delta = q0 as u128 * q1 as u128 / parameters.plaintext() as u128;

        Ok(Self {
            parameters,
            ctx,
            database: std::collections::BTreeMap::new(),
            public_key,
            relin_keys,
            inv_q0_mod_q1: mod_inv(q0 as u128, q1 as u128),
            inv_q1_mod_q0: mod_inv(q1 as u128, q0 as u128),
            rng: thread_rng(),
            threshold: 0,
            squared_indices,
            delta_q0: delta % q0 as u128,
            delta_q1: delta % q1 as u128,
            gamma: 512,
        })
    }

    pub fn enroll_rgsw(&mut self, id: u128, rgsw: RGSWCiphertext) {
        self.database.insert(id, rgsw);
    }

    pub fn verify_with_constraint(
        &mut self, 
        id: u128, 
        // norm: u128, 
        d_packed: &BfvCiphertext, 
        d_act: &BfvCiphertext, 
        d_bin: &BfvCiphertext
    ) -> std::result::Result<BfvCiphertext, Box<dyn std::error::Error>> {
        let rgsw = self.database.get(&id).ok_or(DatabaseError::KeyNotFound(id))?;
        let d_ip = rgsw * d_packed;
        let d_cons = self.construct_constraints(d_packed, d_act, d_bin, &d_ip)?;
        let mut d_verif = self.act_on_lookup_table(d_act, &d_ip);

        d_verif.add_assign(&d_cons);
        d_verif.add_assign(&self.new_mask_with_pk_and_pt_divided_2()?);
        d_verif.mod_switch_to_last_level()?;
        Ok(d_verif)
    }

    pub fn external_product(&self, id: u128, rlwe: &BfvCiphertext) -> std::result::Result<BfvCiphertext, DatabaseError> {
        self.database.get(&id).and_then(|rgsw| Some(rgsw * rlwe)).ok_or(DatabaseError::KeyNotFound(id))
    }

    pub fn get_body_of_innerprod(&self, rlwe: &BfvCiphertext) -> u128 {
        let mut poly = rlwe[0].clone();
        poly.change_representation(PowerBasis);
        
        let coeffs = poly.coefficients();

        let q0 = self.parameters.moduli()[0] as u128;
        let q1 = self.parameters.moduli()[1] as u128;
        let v_q0 = coeffs[(0, 0)] as u128;
        let v_q1 = coeffs[(1, 0)] as u128;
        let m_0 = ((v_q0 * self.inv_q1_mod_q0) % q0) * q1;
        let m_1 = ((v_q1 * self.inv_q0_mod_q1) % q1) * q0;

        (m_0 + m_1) % (q0 * q1)
    }

    pub fn act_on_lookup_table(&self, d_act: &BfvCiphertext, d_ip: &BfvCiphertext) -> BfvCiphertext {
        let q0 = self.parameters.moduli()[0];
        let q1 = self.parameters.moduli()[1];
        let ct_modulus = q0 as u128 * q1 as u128;
        let delta = ct_modulus / self.parameters.plaintext() as u128;
        
        let precision = ((self.parameters.degree() / 4) as f32).log2().round() as usize;
        let plaintext_precision = (self.parameters.plaintext() as f32).log2().round() as usize + 1;  // CAUTION, ceiling thus modulus should not be power-of-two
        let rest_precision = plaintext_precision - precision;
        let prec_cap = 1 << precision;

        let b = self.get_body_of_innerprod(d_ip);
        let b_tilde = ((b % ct_modulus) / delta) as usize >> rest_precision;
        #[cfg(feature = "debug")]
        println!("b: {}, b tilde: {}, prec_cap: {}", b, b_tilde, prec_cap);

        let poly_dim = self.parameters.degree();
        let mut lut = Poly::zero(&self.ctx, PowerBasis);
        let half_pt = self.parameters.plaintext() >> (rest_precision + 1);
        let pt = self.parameters.plaintext() >> rest_precision;
        for a_tilde in 0..prec_cap {
            let v = (a_tilde + b_tilde) % pt as usize;
            if v > self.threshold && v <= half_pt as usize {
                if a_tilde == 0 {
                    lut.coefficients_mut()[(0, 0)] = 1;
                    lut.coefficients_mut()[(1, 0)] = 1;
                } else {
                    lut.coefficients_mut()[(0, poly_dim - a_tilde)] = self.parameters.plaintext() - 1;
                    lut.coefficients_mut()[(1, poly_dim - a_tilde)] = self.parameters.plaintext() - 1;
                }
            }
        }
        lut.change_representation(Representation::Ntt);

        let mut d_out = d_act.clone();
        d_out[0] *= &lut;
        d_out[1] *= &lut;

        d_out
    }

    pub fn construct_constraints(
        &mut self, 
        d_packed: &BfvCiphertext, 
        d_act: &BfvCiphertext, 
        d_bin: &BfvCiphertext, 
        d_ip: &BfvCiphertext
    ) -> std::result::Result<BfvCiphertext, Box<dyn std::error::Error>> {
        
        let mut d_out = BfvCiphertext::new(vec![Poly::zero(&self.ctx, Representation::Ntt); 2], &self.parameters)?;

        self.constraint_norm(&mut d_out, d_packed, d_bin)?;
        self.constraint_monomial(&mut d_out, d_act, d_bin)?;
        self.constarint_index(&mut d_out, d_act, d_bin, &d_ip)?;

        Ok(d_out)
    }

    fn constraint_norm(&mut self, d_out: &mut BfvCiphertext, d_packed: &BfvCiphertext, d_bin: &BfvCiphertext) -> Result<()> {
        let feature_dim = 512_usize;
        let n = self.parameters.degree() / feature_dim;
        let poly_dim = self.parameters.degree();

        // constraint 1 & 2: zero coeffs and dual coeffs
        let mut poly_packed = Poly::zero(&self.ctx, PowerBasis);
        let mut poly_packed_view_mut = poly_packed.coefficients_mut();

        for i in 0..poly_dim {
            let r = self.next_rand_zz_t_ast();

            if i % n == 0 {
                // duals
                let idx_lhs = (poly_dim - i) % poly_dim;
                let idx_rhs = poly_dim - ((feature_dim - 1) * n + 1 - i);
                poly_packed_view_mut[(0, idx_lhs)] = r;
                poly_packed_view_mut[(1, idx_lhs)] = r;

                poly_packed_view_mut[(0, idx_rhs)] = self.parameters.plaintext() - r;
                poly_packed_view_mut[(1, idx_rhs)] = self.parameters.plaintext() - r;
            } else if i % n != 1 {
                // zeros
                poly_packed_view_mut[(0, poly_dim - i)] = r;
                poly_packed_view_mut[(1, poly_dim - i)] = r;
            }
        }
        poly_packed_view_mut[(0, 0)] = self.parameters.plaintext() - poly_packed_view_mut[(0, 0)];
        poly_packed_view_mut[(1, 0)] = self.parameters.plaintext() - poly_packed_view_mut[(1, 0)];
        poly_packed.change_representation(Representation::Ntt);

        let mut d_packed_with_poly = d_packed.clone();
        d_packed_with_poly[0] *= &poly_packed;
        d_packed_with_poly[1] *= &poly_packed;

        d_out.add_assign(&d_packed_with_poly);

        // constraint 3, norm to be gamma
        let mut poly_packed_squared = Poly::zero(&self.ctx, PowerBasis);
        let r = self.next_rand_zz_t_ast();
        poly_packed_squared.coefficients_mut()[(0, poly_dim - (feature_dim - 1) * n - 1)] = r;
        poly_packed_squared.coefficients_mut()[(1, poly_dim - (feature_dim - 1) * n - 1)] = r;
        poly_packed_squared.change_representation(Representation::Ntt);

        let sqrt_n = (self.parameters.degree() as f64).sqrt().round() as u128;
        let bound = self.gamma * sqrt_n + self.parameters.degree() as u128 / 4;
        let gamma_squared = self.gamma * self.gamma;
        let num_bits = (bound as f64).log2().ceil() as usize + 1;
        let norm = (self.parameters.plaintext() as u128 + bound - gamma_squared) % self.parameters.plaintext() as u128;

        let mut poly_packed_squared_offset = Poly::zero(&self.ctx, PowerBasis);
        poly_packed_squared_offset.coefficients_mut()[(0, (feature_dim - 1) * n + 1)] = ((2 * norm * self.delta_q0) % self.parameters.moduli()[0] as u128) as u64;
        poly_packed_squared_offset.coefficients_mut()[(1, (feature_dim - 1) * n + 1)] = ((2 * norm * self.delta_q1) % self.parameters.moduli()[1] as u128) as u64;
        poly_packed_squared_offset.change_representation(Representation::Ntt);

        let mut poly_bin = Poly::zero(&self.ctx, PowerBasis);
        let plaintext_precision = (self.parameters.plaintext() as f32).log2().round() as usize + 1;  // CAUTION, ceiling thus modulus should not be power-of-two
        let precision = ((self.parameters.degree() / 4) as f32).log2().round() as usize;
        let rest_precision = plaintext_precision - precision;
        for i in 0..num_bits {
            let idx = self.squared_indices[rest_precision + 2 + i];
            poly_bin.coefficients_mut()[(0, (feature_dim - 1) * n + 1 - idx)] = (2u64 << i) % self.parameters.moduli()[0];
            poly_bin.coefficients_mut()[(1, (feature_dim - 1) * n + 1 - idx)] = (2u64 << i) % self.parameters.moduli()[1];
        }
        poly_bin.change_representation(Representation::Ntt);
        let mut d_bin = d_bin.clone();
        d_bin[0] *= &poly_bin;
        d_bin[1] *= &poly_bin;

        let mut d_packed_squared = d_packed * d_packed;
        d_packed_squared[0] += &poly_packed_squared_offset;
        d_packed_squared[0] -= &d_bin[0];
        d_packed_squared[1] -= &d_bin[1];
        d_packed_squared[0] *= &poly_packed_squared;
        d_packed_squared[1] *= &poly_packed_squared;
        d_packed_squared[2] *= &poly_packed_squared;
        self.relin_keys.relinearizes(&mut d_packed_squared)?;

        d_out.add_assign(&d_packed_squared);

        Ok(())
    }

    fn constraint_monomial(&mut self, d_out: &mut BfvCiphertext, d_act: &BfvCiphertext, d_bin: &BfvCiphertext) -> Result<()> {
        let poly_dim = self.parameters.degree();
        
        // constarint 1 & 2
        let mut poly_act = Poly::zero(&self.ctx, PowerBasis);
        for i in 0..poly_dim/4 {
            // duals
            let r = self.next_rand_zz_t_ast();
            poly_act.coefficients_mut()[(0, (poly_dim - i) % poly_dim)] = r;
            poly_act.coefficients_mut()[(1, (poly_dim - i) % poly_dim)] = r;
            poly_act.coefficients_mut()[(0, poly_dim - (poly_dim * 3 / 4 - 1 - i))] = self.parameters.plaintext() - r;
            poly_act.coefficients_mut()[(1, poly_dim - (poly_dim * 3 / 4 - 1 - i))] = self.parameters.plaintext() - r;

            // zeros
            let r = self.next_rand_zz_t_ast();
            poly_act.coefficients_mut()[(0, (poly_dim - i - poly_dim/4) % poly_dim)] = r;
            poly_act.coefficients_mut()[(1, (poly_dim - i - poly_dim/4) % poly_dim)] = r;
            let r = self.next_rand_zz_t_ast();
            poly_act.coefficients_mut()[(0, (poly_dim - i - poly_dim*3/4) % poly_dim)] = r;
            poly_act.coefficients_mut()[(1, (poly_dim - i - poly_dim*3/4) % poly_dim)] = r;
        }

        poly_act.coefficients_mut()[(0, 0)] = self.parameters.plaintext() - poly_act.coefficients_mut()[(0, 0)];
        poly_act.coefficients_mut()[(1, 0)] = self.parameters.plaintext() - poly_act.coefficients_mut()[(1, 0)];
        poly_act.change_representation(Representation::Ntt);

        let mut d_act_with_poly = d_act.clone();
        d_act_with_poly[0] *= &poly_act;
        d_act_with_poly[1] *= &poly_act;
        d_out.add_assign(&d_act_with_poly);

        // constraint 3 & 4
        let mut poly_act_squared = Poly::zero(&self.ctx, PowerBasis);
        for i in (poly_dim/2)..(poly_dim*3/4) {
            let r = self.next_rand_zz_t_ast();
            poly_act_squared.coefficients_mut()[(0, poly_dim - i)] = r;
            poly_act_squared.coefficients_mut()[(1, poly_dim - i)] = r;
        }
        poly_act_squared.change_representation(Representation::Ntt);

        let mut poly_act_squared_offset = Poly::zero(&self.ctx, PowerBasis);
        poly_act_squared_offset.coefficients_mut()[(0, poly_dim * 3 / 4 - 1)] = ((2 * self.delta_q0) % self.parameters.moduli()[0] as u128) as u64;
        poly_act_squared_offset.coefficients_mut()[(1, poly_dim * 3 / 4 - 1)] = ((2 * self.delta_q1) % self.parameters.moduli()[1] as u128) as u64;
        poly_act_squared_offset.change_representation(Representation::Ntt);

        let mut d_act_squared = d_act * d_act;
        d_act_squared[0] -= &poly_act_squared_offset;
        d_act_squared[0] *= &poly_act_squared;
        d_act_squared[1] *= &poly_act_squared;
        d_act_squared[2] *= &poly_act_squared;
        self.relin_keys.relinearizes(&mut d_act_squared)?;
        
        d_out.add_assign(&d_act_squared);

        // constraint 5 & 6
        let sqrt_n = (self.parameters.degree() as f64).sqrt().round() as u128;
        let bound = self.gamma * sqrt_n + self.parameters.degree() as u128 / 4;
        let num_bits = (bound as f64).log2().ceil() as usize + 1;

        let precision = ((self.parameters.degree() / 4) as f32).log2().round() as usize;
        let plaintext_precision = (self.parameters.plaintext() as f32).log2().round() as usize + 1;  // CAUTION, ceiling thus modulus should not be power-of-two
        let rest_precision = plaintext_precision - precision;
        let mut poly_bin = Poly::zero(&self.ctx, Representation::PowerBasis);
        for i in 0..poly_dim {
            let r = self.next_rand_zz_t_ast();
            poly_bin.coefficients_mut()[(0, i)] = r;
            poly_bin.coefficients_mut()[(1, i)] = r;
        }

        let mut poly_bin_squared = Poly::zero(&self.ctx, Representation::PowerBasis);
        for i in 0..rest_precision+2+num_bits {
            let idx = self.squared_indices[i];
            poly_bin_squared.coefficients_mut()[(0, (poly_dim - 2 * idx) % poly_dim)] = self.parameters.plaintext() - poly_bin.coefficients()[(0, (poly_dim - idx) % poly_dim)];
            poly_bin_squared.coefficients_mut()[(1, (poly_dim - 2 * idx) % poly_dim)] = self.parameters.plaintext() - poly_bin.coefficients()[(1, (poly_dim - idx) % poly_dim)];
        }

        poly_bin.change_representation(Representation::Ntt);
        poly_bin_squared.change_representation(Representation::Ntt);
        let mut d_bin_squared = d_bin * d_bin;
        d_bin_squared[0] *= &poly_bin_squared;
        d_bin_squared[1] *= &poly_bin_squared;
        d_bin_squared[2] *= &poly_bin_squared;
        self.relin_keys.relinearizes(&mut d_bin_squared)?;
        d_out.add_assign(&d_bin_squared);

        let mut d_bin = d_bin.clone();
        d_bin[0] *= &poly_bin;
        d_bin[1] *= &poly_bin;
        d_out.add_assign(&d_bin);

        Ok(())
    }

    fn constarint_index(&mut self, d_out: &mut BfvCiphertext, d_act: &BfvCiphertext, d_bin: &BfvCiphertext, d_ip: &BfvCiphertext) -> Result<()> {
        let poly_dim = self.parameters.degree();
        let q0 = self.parameters.moduli()[0];
        let q1 = self.parameters.moduli()[1];

        let mut d_recons = d_act.clone();
        // lut_sim
        let plaintext_precision = (self.parameters.plaintext() as f32).log2().round() as usize + 1;  // CAUTION, ceiling thus modulus should not be power-of-two
        let precision = ((self.parameters.degree() / 4) as f32).log2().round() as usize;
        let rest_precision = plaintext_precision - precision;
        let prec_cap = 1 << precision;
        
        let mut lut_sim = Poly::zero(&self.ctx, PowerBasis);
        for i in 1..prec_cap {
            // neg coeffs
            lut_sim.coefficients_mut()[(0, poly_dim - i as usize)] = (i << rest_precision) % q0;
            lut_sim.coefficients_mut()[(1, poly_dim - i as usize)] = (i << rest_precision) % q1;
        }
        lut_sim.change_representation(Representation::Ntt);

        d_recons[0] *= &lut_sim;
        d_recons[1] *= &lut_sim;

        // lut bin
        let mut lut_bin = Poly::zero(&self.ctx, PowerBasis);
        for i in 0..rest_precision {
            let idx = self.squared_indices[i];
            // neg coeffs
            if i == 0 {
                lut_bin.coefficients_mut()[(0, 0)] = self.parameters.plaintext() - 1;
                lut_bin.coefficients_mut()[(1, 0)] = self.parameters.plaintext() - 1;
            } else {
                lut_bin.coefficients_mut()[(0, poly_dim - idx)] = 1 << i;
                lut_bin.coefficients_mut()[(1, poly_dim - idx)] = 1 << i;
            }
        }
        // correction flag
        // add b_delta % delta
        let delta = q0 as u128 * q1 as u128 / self.parameters.plaintext() as u128;
        let b_delta = self.get_body_of_innerprod(&d_ip) % delta;
        
        // correction flag, set a very loose boundary of Delta/4
        if b_delta > delta * 3 / 4 {
            // u_plus
            let idx = self.squared_indices[rest_precision];
            lut_bin.coefficients_mut()[(0, poly_dim - idx)] = 1;
            lut_bin.coefficients_mut()[(1, poly_dim - idx)] = 1;
        } else if b_delta < delta / 4 {
            // u_minus
            let idx = self.squared_indices[rest_precision + 1];
            lut_bin.coefficients_mut()[(0, poly_dim - idx)] = self.parameters.plaintext() - 1;
            lut_bin.coefficients_mut()[(1, poly_dim - idx)] = self.parameters.plaintext() - 1;
        }
        lut_bin.change_representation(Representation::Ntt);
        let mut d_bin = d_bin.clone();
        d_bin[0] *= &lut_bin;
        d_bin[1] *= &lut_bin;
        d_recons.add_assign(&d_bin);

        d_recons[1] += &d_ip[1];

        let mut b_poly = Poly::zero(&self.ctx, PowerBasis);
        // minus one for that the client encrypting `a_tilde` rather than `-a_tilde`
        b_poly.coefficients_mut()[(0, 0)] = ((b_delta + (self.parameters.plaintext() as u128 - 1) * delta) % q0 as u128) as u64;
        b_poly.coefficients_mut()[(1, 0)] = ((b_delta + (self.parameters.plaintext() as u128 - 1) * delta) % q1 as u128) as u64;
        b_poly.change_representation(Representation::Ntt);
        d_recons[0] += &b_poly;

        let mut r_poly = Poly::zero(&self.ctx, PowerBasis);
        r_poly.coefficients_mut()[(0, 0)] = self.next_rand_zz_t_ast();
        r_poly.coefficients_mut()[(1, 0)] = r_poly.coefficients()[(0, 0)];
        r_poly.change_representation(Representation::Ntt);
        d_recons[0] *= &r_poly;
        d_recons[1] *= &r_poly;

        d_out.add_assign(&d_recons);


        Ok(())
    }

    pub fn new_mask_with_pk_and_pt_divided_2(&mut self) -> Result<BfvCiphertext> {        
        let cleartext = (0..self.parameters.degree()).map(|i| {
            if i == 0 {
                self.next_rand_zz_t() & 0xFFFFFFFF_FFFFFFFE
                // 0
            } else {
                self.next_rand_zz_t()
            }
        }).collect::<Vec<_>>();
        let pt = bfv::Plaintext::try_encode(&cleartext, bfv::Encoding::poly(), &self.parameters)?;
        let ct = self.public_key.try_encrypt(&pt, &mut self.rng)?;

        Ok(ct)
    }

    fn next_rand_zz_t_ast(&mut self) -> u64 {
        self.rng.gen_range(1..self.parameters.plaintext())
    }

    fn next_rand_zz_t(&mut self) -> u64 {
        // to prevent overflow
        self.rng.gen_range(0..self.parameters.plaintext()-1)
    }
}


