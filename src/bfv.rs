use fhe::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_traits::*;

use rand::{rngs::OsRng, thread_rng};
use tfhe::core_crypto::prelude::IntoContainerOwned;
use std::error::Error;

// pub fn bfv_ciphertext_into_container_u128(ct: &Ciphertext) -> Vec<u128> {

// }