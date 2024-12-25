use ppverif_fhe::*;
use std::time::Instant;

pub fn main() {
    let mut mal_client = MalClient::new(DEFAULT_MALICIOUS_PARAMETER);
    let pk = mal_client.new_rlwe_public_key();
    let mut mal_server = MalServer::new(DEFAULT_MALICIOUS_PARAMETER, pk);

    let mut features = vec![0.0f32; (DEFAULT_MALICIOUS_PARAMETER.glwe_size.0 - 1) * DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    for i in 0..DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0 {
        features[(DEFAULT_MALICIOUS_PARAMETER.glwe_size.0 - 1) * i] = i as f32;
    }

    let rlwe = mal_client.encrypt_rlwe(&features, 1.0f32);
    let glwe = extract_glwe_sample_from_rlwe_ciphertext(&rlwe, DEFAULT_MALICIOUS_PARAMETER.polynomial_size);
    
    let decrypted_rlwe = mal_client.decrypt_rlwe(&rlwe);
    println!("{:?}, len: {}", decrypted_rlwe, decrypted_rlwe.len());

    let decrypted = mal_client.decrypt_glwe(&glwe);
    println!("{:?}, len: {}", decrypted, decrypted.len());


    let mut features1 = vec![0.0f32; DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    let mut features2 = vec![0.0f32; DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    features2[0] = 1.0f32;
    features2[1] = 1.0f32;
    features1[0] = 2.0f32;
    features1[1] = 3.0f32;
    features1[511] = -5.0f32;

    let ggsw = mal_client.encrypt_new_template(0, &features2, 1.0);
    let glwe = mal_client.encrypt_glwe(&features1, 1.0);
    mal_server.enroll_ggsw(0, ggsw);
    let ext_ct = mal_server.external_product(0, &glwe).unwrap();
    let ext_pt = mal_client.decrypt_glwe(&ext_ct);
    println!("{:?}", ext_pt);
}