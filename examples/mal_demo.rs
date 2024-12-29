use ppverif_fhe::*;
use rand::{thread_rng, Rng};
use tfhe::{boolean::prelude::{PolynomialSize, DEFAULT_PARAMETERS_KS_PBS}, core_crypto::prelude::*};
use std::time::Instant;

pub fn main() {
    let mut mal_client = MalClient::new(DEFAULT_MALICIOUS_PARAMETER);
    let pk = mal_client.new_rlwe_public_key();
    let rlk = mal_client.new_rlwe_relinearizaion_keys();
    let mut mal_server = MalServer::new(DEFAULT_MALICIOUS_PARAMETER, pk, rlk);

    let mut features = vec![0.0f32; (DEFAULT_MALICIOUS_PARAMETER.glwe_size.0 - 1) * DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    for i in 0..DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0 {
        features[(DEFAULT_MALICIOUS_PARAMETER.glwe_size.0 - 1) * i] = i as f32;
    }

    let rlwe = mal_client.encrypt_rlwe(&features, 1.0f32);
    let glwe = extract_glwe_sample_from_rlwe_ciphertext(&rlwe, DEFAULT_MALICIOUS_PARAMETER.polynomial_size);
    
    // let decrypted_rlwe = mal_client.decrypt_rlwe(&rlwe);
    // println!("{:?}, len: {}", decrypted_rlwe, decrypted_rlwe.len());

    // let decrypted = mal_client.decrypt_glwe(&glwe);
    // println!("{:?}, len: {}", decrypted, decrypted.len());


    let mut features1 = vec![0.0f32; DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    let mut features2 = vec![0.0f32; DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    features2[0] = 1.0f32;
    features2[1] = 1.0f32;
    features1[0] = 2.0f32;
    features1[1] = 3.0f32;
    features1[511] = -5.0f32;

    let ggsw = mal_client.encrypt_new_template_ggsw(0, &features2, 1.0);
    let glwe = mal_client.encrypt_glwe(&features1, 1.0);
    mal_server.enroll_ggsw(0, ggsw);
    let ext_ct = mal_server.external_product(0, &glwe).unwrap();
    let ext_pt = mal_client.decrypt_glwe(&ext_ct);
    println!("{:?}", ext_pt);


    let mut features_rlwe_1 = vec![0.0f32; (DEFAULT_MALICIOUS_PARAMETER.glwe_size.0 - 1) * DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    let mut features_rlwe_2 = features_rlwe_1.clone();
    // let mut features_rlwe_2 = vec![0u128; (DEFAULT_MALICIOUS_PARAMETER.glwe_size.0 - 1) * DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0 * 2];
    features_rlwe_1[0] = 2.0;
    features_rlwe_1[1] = 3.0;
    features_rlwe_2[0] = 3.0;
    features_rlwe_2[2] = 5.0;
    features_rlwe_2[4095] = -1.0;
    // features_rlwe_2[(DEFAULT_MALICIOUS_PARAMETER.glwe_size.0 - 1) * DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0] = 3 * DEFAULT_MALICIOUS_PARAMETER.delta;
    let rlwe1 = mal_client.encrypt_rlwe(&features_rlwe_1, 1.0);
    let rlwe2 = mal_client.encrypt_rlwe(&features_rlwe_2, 1.0);
    // let rlwe2 = GlweCiphertext::from_container(features_rlwe_2, rlwe1.polynomial_size(), CiphertextModulus::new_native());
    // let mut rlwe_mul = GlweCiphertext::new(0, GlweSize(3), PolynomialSize(features_rlwe_1.len()), CiphertextModulus::new_native());
    // let rlwe_mul = rlwe_multiplication_u96(&rlwe1, &rlwe2, DEFAULT_MALICIOUS_PARAMETER.plaintext_modulus);
    // let rlwe_relin = mal_server.relinearize(&rlwe_mul);
    // let pt_relin = mal_client.decrypt_rlwe(&rlwe_relin);
    // let rlwe_dir_dec = mal_client.decrypt_rlwe_multiplied(&rlwe_mul);
    // println!("relin: {:?}", &pt_relin[..10]);
    // println!("dir_dec: {:?}", &rlwe_dir_dec[..10]);

    
    let mut f1 = vec![0.0f32; DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    let mut rng = thread_rng();
    for vi in f1.iter_mut().take(1) {
        *vi = 5.0; rng.gen_range(0..10) as f32;
    }
    let (d_packed, norm) = mal_client.encrypt_new_template_rlwe(&f1, 1.0);
    let d_glwe = extract_glwe_sample_from_rlwe_ciphertext(&d_packed, DEFAULT_MALICIOUS_PARAMETER.polynomial_size);
    let d_ext = mal_server.external_product(0, &d_glwe).unwrap();
    let d_ext_lwe = extract_lwe_sample_from_glwe_ciphertext_under_rlwe_secret_key(&d_ext);
    let pt_ext = mal_client.decrypt_glwe(&d_ext);    
    let pt_ext_lwe = mal_client.decrypt_lwe(&d_ext_lwe);
    println!("pt_ext: {:?}", pt_ext);
    println!("pt_ext_lwe: {:?}", pt_ext_lwe);

    
    // let mut c_cons = LweCiphertext::new(0, LweSize(d_packed.polynomial_size().0 + 1), CiphertextModulus::new_native());
    // mal_server.constraint_norm(&mut c_cons, &d_packed, norm);    

    let (d_act, d_bin) = mal_client.encrypt_new_lookup_tables(0, &d_packed);
    let time = Instant::now();
    // mal_server.constraint_monomial(&mut c_cons, &d_act, &d_bin);
    // mal_server.constarint_index(&mut c_cons, &d_act, &d_bin);
    let c_cons_64 = mal_server.construct_constraints(0, norm, d_packed, d_act, d_bin).unwrap();

    let elapsed = time.elapsed();
    println!("time: {:?} millis", elapsed.as_micros() as f32 / 1000.0);

    // let c_cons_64 = LweCiphertext::from_container(c_cons.as_ref().into_iter().map(|v| ((v >> 32) & 0xFFFFFFFF_FFFFFFFF) as u64).collect::<Vec<_>>(), CiphertextModulus::new_native());
    // let cons_pt = mal_client.decrypt_lwe(&c_cons);
    let cons_pt_64 = mal_client.decrypt_lwe(&c_cons_64);
    println!("cons pt 64: {}, norm: {}", cons_pt_64, norm);
}