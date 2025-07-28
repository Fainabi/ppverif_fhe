use std::collections::*;
use rand::*;
use polynomial_algorithms::polynomial_wrapping_add_mul_assign;
use tfhe::core_crypto::commons::math::random::RandomGenerator;
use tfhe::core_crypto::prelude::*;

use crate::{params::*, Client};

pub struct Server {
    ip_param: GlweParameter<u64>,
    br_param: GlweParameter<u16>,

    // any data structure
    database: BTreeMap<u128, Vec<GlweBody<Vec<u64>>>>,

    // for masking the sample extraction to product inner product
    glwe_pk: Vec<GlweCiphertextOwned<u16>>,

    noise_distribution: Gaussian<f64>,
    noise_seeder: Box<dyn Seeder>,

    precision: usize,
}

impl Server {
    pub fn new(ip_param: GlweParameter<u64>, br_param: GlweParameter<u16>, glwe_pk: Vec<GlweCiphertextOwned<u16>>) -> Self {
        let noise_distribution = Gaussian::from_standard_dev(StandardDev(br_param.std_dev), 0.0);
        Self {
            ip_param,
            br_param,
            database: BTreeMap::new(),
            glwe_pk,
            noise_distribution,
            noise_seeder: new_seeder(),
            precision: 8,
        }
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
                    .fold(0u64, |sum_body, (template_glev_bodies, query_poly)| {
                        // let query_container = query_poly.into_container();

                        let polly_mul = template_glev_bodies
                            .iter()
                            .enumerate()
                            .map(|(beta_idx, template_body)| {
                                let idx = self.ip_param.decomposition_level_count.0 - beta_idx;
                                let template_container = template_body.as_polynomial().into_container();
                                let query_container: Vec<_> = query_poly
                                    .into_container()
                                    .iter()
                                    .map(|&v| {
                                        let v = v >> (self.ip_param.decomposition_base_log.0 * idx);
                                        v & ((1u64 << self.ip_param.decomposition_base_log.0) - 1)
                                    })
                                    .collect();

                                let mut sum = query_container[0].wrapping_mul(template_container[0]);
                                for (&vi, &wi) in template_container[1..]
                                    .iter()
                                    .zip(query_container[1..].iter().rev())
                                {
                                    sum = sum.wrapping_sub(vi.wrapping_mul(wi));
                                }

                                sum
                            })
                            .fold(0u64, |acc, x| acc.wrapping_add(x));

                        sum_body.wrapping_add(polly_mul)
                    });

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

    pub fn collect_database_into_matrix(&self) -> Option<nalgebra::DMatrix<u64>> {
        if self.database.len() == 0 {
            return None;
        }

        let col_size = self.database
            .iter()
            .map(|(_k, v)| {
                v.iter()
                 .map(|vi| vi.as_ref().iter().count())
                 .sum()
            })
            .next()
            .unwrap();

        let mat = nalgebra::DMatrix::from_row_iterator(
            self.database.len(), 
            col_size,
            self.database
                .iter()
                .flat_map(|(_k, ggsw_dec_mask)| {
                    ggsw_dec_mask
                        .iter()
                        .flat_map(|ct_tensor_g| {
                            let cont = ct_tensor_g.as_polynomial().into_container();

                            (0..cont.len()).map(|i| if i == 0 { cont[0] } else { cont[cont.len() - i].wrapping_neg() })
                        })
                })
            );
        
        Some(mat)
    }
}


#[test]
fn test_server_db() {
    let ip_param = DEFAULT_INNER_PRODUCT_PARAMETER;
    let br_param = DEFAULT_BLIND_ROTATION_PARAMETER;
    let mut client = Client::new(ip_param, br_param);
    let glwe_pk = client.new_glwe_public_keys_br();
    let mut server = Server::new(ip_param, br_param, glwe_pk);

    for i in 0..10 {
        let mut f = vec![0.0; 512];
        f[i] = 1.0f32;
        
        let i = i as u128;
        let template = client.encrypt_new_template(i, &f, 512.0);
        server.enroll(i, template);
    }

    let mat_server = server.collect_database_into_matrix().unwrap();
    let mut query_feature = vec![0.0; 512];
    query_feature[0] = 0.3f32;

    let new_glwe = client.encrypt_glwe(&query_feature, 512.0);
    let new_glwe_cont = new_glwe.clone().into_container();
    let n = new_glwe_cont.len();
    let vec_glwe = nalgebra::DVector::from_iterator(n, new_glwe_cont.into_iter().map(|x| x >> 32));

    server.compute_innerprod_body(0, new_glwe.as_view());
    client.transform_mask_from_database(0, new_glwe.as_view());

    let mul_server_mat = mat_server * vec_glwe.clone();
    println!("{:?}", mul_server_mat);
}
