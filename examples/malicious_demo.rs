use std::error::Error;
use ppverif_fhe::*;
use std::time::Instant;


fn main() -> Result<(), Box<dyn Error>> {
    let mut malclient = MaliciousClient::new()?;
    let public_key = malclient.new_public_key();
    let relin_keys = malclient.new_relin_keys()?;

    let mut malserver = AntiMalServer::new(public_key, relin_keys)?;
    
    let mut features1 = vec![0.0f32; DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    let mut features2 = vec![0.0f32; DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0];
    for i in 1..DEFAULT_MALICIOUS_PARAMETER.polynomial_size.0 {
        features1[i] = i as f32;
        features2[i] = -((512 - i) as f32);
    }
    features2[0] = 1.0;
    normalize(&mut features1);
    normalize(&mut features2);

    // enrollment
    let (rgsw, rgsw_masks) = malclient.encrypt_rgsw_ciphertext(&features1, 512.0)?;
    malclient.enroll_rgsw_masks(0, rgsw_masks);
    malserver.enroll_rgsw(0, rgsw);

    // encrypt
    let now = Instant::now();
    let (d_packed, norm) = malclient.encrypt_rlwe_ciphertext(&features2, 512.0)?;
    let elapsed = now.elapsed();
    println!("encryption d_packed time: {} millis", elapsed.as_micros() as f32 / 1000.0);
    
    let now = Instant::now();
    let (d_act, d_bin) = malclient.encrypt_new_lookup_tables(0, 512, norm, &d_packed)?;
    let elapsed = now.elapsed();
    println!("construct two cts time: {} millis", elapsed.as_micros() as f32 / 1000.0);
    
    // verification
    let now = Instant::now();
    let d_cons = malserver.verify_with_constraint(0, &d_packed, &d_act, &d_bin)?;
    // let d_cons = malserver.construct_constraints(0, norm, &d_packed, &d_act, &d_bin)?;
    let elapsed = now.elapsed();
    println!("verification time: {} millis", elapsed.as_micros() as f32 / 1000.0);

    // decryption
    let now = Instant::now();
    let pt_cons = malclient.decrypt(&d_cons)?;
    let elapsed = now.elapsed();
    println!("verification result: {:?}, dec time: {} millis", pt_cons[0], elapsed.as_micros() as f32 / 1000.0);

    Ok(())
}