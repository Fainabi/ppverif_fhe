use fhe::Result;
use ppverif_fhe::*;
use std::time::Instant;
use fhe::bfv;
use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey, RGSWCiphertext, RelinearizationKey};
use rand::{rngs::OsRng, thread_rng};
use fhe_traits::*;


fn main() -> Result<()> {
    let parameters = bfv::BfvParametersBuilder::new()
        .set_degree(4096)
        .set_moduli(&[0x3fffffff000001_u64, 0x3fffffff004001_u64, 0x3fffffff04e001, 0x3fffffff058001, 0x3fffffff07c001])
        .set_plaintext_modulus(0x88001)
        .build_arc()?;


    let mut client = Client::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER, parameters.clone());
    let pk = client.new_glwe_public_keys_br();
    let bfv_rlk = client.new_bfv_relinearizatio_key()?;
    let mut server = Server::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER, parameters.clone(), pk, bfv_rlk)?;

    let vals = (0..4096).collect::<Vec<_>>();
    let cts = client.encrypt_rlwe_for_binary_decomposition(&vals)?;

    Ok(())
}