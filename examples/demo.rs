use ppverif_fhe::*;

pub fn main() {
    let mut client = Client::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER);

    client.encrypt_new_template(0, &vec![0f32; 512], 1f32);
}