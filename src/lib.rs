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

    let mut rust_floats: Vec<f32> = vec![];
    let repeat_time = 100;

    for ip_param in [
        INNER_PRODUCT_PARAMETER_DIMENSION_128,
        INNER_PRODUCT_PARAMETER_DIMENSION_256,
        INNER_PRODUCT_PARAMETER_DIMENSION_512,
        INNER_PRODUCT_PARAMETER_DIMENSION_1024,
        INNER_PRODUCT_PARAMETER_DIMENSION_2048,
        INNER_PRODUCT_PARAMETER_DIMENSION_4096,
    ] {
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
        let now = Instant::now();
        for i in 1..repeat_time {
            client.enroll_ggsw_masks(0);
        }
        client.enroll_ggsw_masks(0);
        let elapsed = now.elapsed();
        let preprocessing_time = elapsed.as_nanos() as f32 / 1000.0 / (repeat_time as f32);

        // Client
        let now = Instant::now();
        for i in 1..repeat_time {
            let _ = client.encrypt_glwe(&f2, 512.0);
        }
        let query_ct = client.encrypt_glwe(&f2, 512.0);
        let elapsed = now.elapsed();
        let encryt_query_ct_time = elapsed.as_nanos() as f32 / 1000.0 / (repeat_time as f32);

        let now = Instant::now();
        for i in 1..repeat_time {
            client.transform_mask_to_body_from_database(0, query_ct.as_view());
        }
        client.transform_mask_to_body_from_database(0, query_ct.as_view());
        let elapsed = now.elapsed();
        let client_innerprod_time = elapsed.as_nanos() as f32 / 1000.0 / (repeat_time as f32);

        let now = Instant::now();
        for i in 1..repeat_time {
            let _ = client.transform_mask_from_database(0, query_ct.as_view()).unwrap();
        }
        let lut_ct = client.transform_mask_from_database(0, query_ct.as_view()).unwrap();
        let elapsed = now.elapsed();
        let generate_lookup_table_time = elapsed.as_nanos() as f32 / 1000.0 / (repeat_time as f32);

        // Server
        let verif_ct = server.verify(0, &query_ct, &lut_ct).unwrap();

        // Client
        let now = Instant::now();
        for i in 1..repeat_time {
            let _ = client.decrypt_lwe(&verif_ct);
        }
        let verif_res = client.decrypt_lwe(&verif_ct);
        let elapsed = now.elapsed();
        let decryption_time = elapsed.as_nanos() as f32 / 1000.0 / (repeat_time as f32);

        vec![preprocessing_time, encryt_query_ct_time, client_innerprod_time, generate_lookup_table_time, decryption_time]
            .into_iter()
            .for_each(|t| rust_floats.push(t));
    }
    
    let float_array: jni::objects::JPrimitiveArray<'_, f32> = _env.new_float_array(rust_floats.len() as i32).unwrap();
    _env.set_float_array_region(&float_array, 0, &rust_floats).unwrap();
    float_array.into_raw()
}
