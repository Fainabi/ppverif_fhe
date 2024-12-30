use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey, RGSWCiphertext, RelinearizationKey};
use fhe::proto::bfv::{RgswCiphertext as RGSWCiphertextProto};
use fhe_traits::*;
use std::{borrow::BorrowMut, time::Instant};

use rand::{rngs::OsRng, thread_rng};
use std::error::Error;



fn main() -> Result<(), Box<dyn Error>> {
    let parameters = BfvParametersBuilder::new()
            .set_degree(4096)
            .set_moduli(&[0x3fffffff000001, 0x3fffffff004001])
            .set_plaintext_modulus(1 << 11)
            .build_arc()?;
    println!("{:?}", parameters);

    let mut rng = thread_rng();

    let secret_key = SecretKey::random(&parameters, &mut OsRng);
    let public_key = PublicKey::new(&secret_key, &mut rng);
    let relin_keys = RelinearizationKey::new(&secret_key, &mut rng)?;

    let plaintext_1 = Plaintext::try_encode(&[20_u64], Encoding::poly(), &parameters)?;
    let plaintext_2 = Plaintext::try_encode(&[-7_i64, -8i64], Encoding::poly(), &parameters)?;

    let mut ciphertext_1: Ciphertext = secret_key.try_encrypt(&plaintext_1, &mut rng)?;
    let ciphertext_2: Ciphertext = public_key.try_encrypt(&plaintext_2, &mut rng)?;
    println!("{:?}", ciphertext_1);

    let rgsw: RGSWCiphertext = secret_key.try_encrypt(&plaintext_2, &mut rng)?;
    let time = Instant::now();
    let mut result = &ciphertext_1 * &ciphertext_2;
    
    let muled = &rgsw * &ciphertext_1;
    // result.mod_switch_to_last_level()?;
    relin_keys.relinearizes(&mut result)?;
    
    let elapsed = time.elapsed();
    println!("hommul: {}ms", elapsed.as_micros() as f64 / 1000.0);
    println!("len: {}", rgsw.to_bytes().len());
    // rgsw.into();
    
    let res2 = &result + &muled;
    // println!("{:?}", result);
    // let decrypted_plaintext = secret_key.try_decrypt(&result)?;
    // let decrypted_vector = Vec::<i64>::try_decode(&decrypted_plaintext, Encoding::poly())?;

    println!("proto: {:?}", RGSWCiphertextProto::from(&rgsw).ksk0.unwrap().c0[0].len());
    
    // let muled = &rgsw * &ciphertext_1;
    let dec = secret_key.try_decrypt(&res2)?;
    // println!("{:?}", dec);
    // assert_eq!(decrypted_vector[0], -140);

    let x = &mut *ciphertext_1;
    let x0 = x[0].coefficients();

    Ok(())
}