use ppverif_fhe::*;
use std::time::Instant;

pub fn main() {
    let mut client = Client::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER);
    let pk = client.new_glwe_public_keys_br();
    let mut server = Server::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER, pk);

    // features
    let mut f1 = vec![0f32; 512];
    let mut f2 = vec![0f32; 512];
    f1[1] = 1f32;
    f2[511] = -0.2f32;

    // enrollment
    let template = client.encrypt_new_template(0, &f1, 512.0);
    server.enroll(0, template);
    client.enroll_ggsw_masks(0);

    // verification
    // Client
    let now = Instant::now();
    let query_ct = client.encrypt_glwe(&f2, 512.0);
    let elapsed = now.elapsed();
    println!("Encrypting query ct {} micros", elapsed.as_nanos() as f32 / 1000.0);

    let now = Instant::now();
    let lut_ct = client.transform_mask_from_database(0, query_ct.as_view()).unwrap();
    let elapsed = now.elapsed();
    println!("Generate lookup table: {} micros", elapsed.as_nanos() as f32 / 1000.0);

    // Server
    let now = Instant::now();
    let verif_ct = server.verify(0, &query_ct, &lut_ct).unwrap();
    let elapsed = now.elapsed();
    println!("verification time: {} micros", elapsed.as_nanos() as f32 / 1000.0);

    // Client
    let now = Instant::now();
    let verif_res = client.decrypt_lwe(&verif_ct);
    let elapsed = now.elapsed();
    println!("decrytion time: {} micros", elapsed.as_nanos() as f32 / 1000.0);

    println!("verif res: {}", verif_res);
}
