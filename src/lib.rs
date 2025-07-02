mod params;
mod client;
mod server;
mod seeder;
mod utils;

pub use client::Client;
pub use params::{DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER};
pub use server::Server;
pub use utils::normalize;

use std::time::Instant;
use jni::objects::JClass;
use jni::sys::jfloatArray;
use jni::JNIEnv;


#[unsafe(no_mangle)]
pub extern "system" fn Java_com_example_ppverif_RustBridge_testClient(
    mut _env: JNIEnv,
    _class: JClass
) -> jfloatArray {

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
    // println!("Encrypting query ct {} micros", elapsed.as_nanos() as f32 / 1000.0);
    let encryt_query_ct_time = elapsed.as_nanos() as f32 / 1000.0;

    let now = Instant::now();
    let lut_ct = client.transform_mask_from_database(0, query_ct.as_view()).unwrap();
    let elapsed = now.elapsed();
    // println!("Generate lookup table: {} micros", elapsed.as_nanos() as f32 / 1000.0);
    let generate_lookup_table_time = elapsed.as_nanos() as f32 / 1000.0;

    // Server
    let now = Instant::now();
    let verif_ct = server.verify(0, &query_ct, &lut_ct).unwrap();
    let elapsed = now.elapsed();
    // println!("verification time: {} micros", elapsed.as_nanos() as f32 / 1000.0);
    let verification_time = elapsed.as_nanos() as f32 / 1000.0;

    // Client
    let now = Instant::now();
    let verif_res = client.decrypt_lwe(&verif_ct);
    let elapsed = now.elapsed();
    // println!("decrytion time: {} micros", elapsed.as_nanos() as f32 / 1000.0);
    let decryption_time = elapsed.as_nanos() as f32 / 1000.0;

    // println!("verif res: {}", verif_res);

    let rust_floats: Vec<f32> = vec![encryt_query_ct_time, generate_lookup_table_time, verification_time,
    decryption_time];
    let float_array: jni::objects::JPrimitiveArray<'_, f32> = _env.new_float_array(rust_floats.len() as i32).unwrap();
    _env.set_float_array_region(&float_array, 0, &rust_floats).unwrap();
    float_array.into_raw()
}
