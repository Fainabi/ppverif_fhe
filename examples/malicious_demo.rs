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
    features2[0] = 6100.0f32;
    features2[1] = 1.0f32;
    features1[0] = 2.0f32;
    features1[1] = 3.0f32;
    // features1[511] = -5.0f32;

    let (rgsw, rgsw_masks) = malclient.encrypt_rgsw_ciphertext(&features1, 1.0)?;
    malclient.enroll_rgsw_masks(0, rgsw_masks);
    malserver.enroll_rgsw(0, rgsw.clone());
    let (d_packed, norm) = malclient.encrypt_rlwe_ciphertext(&features2, 1.0)?;
    let (d_act, d_bin) = malclient.encrypt_new_lookup_tables(0, &d_packed)?;

    let d_ip = malserver.external_product(0, &d_packed)?;
    let pt_ip = malclient.decrypt(&d_ip)?;
    println!("pt_ip: {}", pt_ip[0]);

    let time = Instant::now();
    let d_cons = malserver.verify_with_constraint(0, norm, &d_packed, &d_act, &d_bin)?;
    // let d_cons = malserver.construct_constraints(0, norm, &d_packed, &d_act, &d_bin)?;
    let elapsed = time.elapsed();

    let pt_cons = malclient.decrypt(&d_cons)?;
    println!("pt cons: {:?}", pt_cons[0]);
    println!("cons time: {} millis", elapsed.as_micros() as f32 / 1000.0);

    // let pt_mask = malclient.compute_mask_for_innerprod(0, &rlwe)?;
    // println!("pt mask with body set to zero: {}", pt_mask);

    // let pt_bodyhalf = malclient.compute_body_for_innerprod(&rgsw_body, &rlwe)?;
    // println!("pt body with mask set to zero: {}", pt_bodyhalf);

    // let ext_rlwe = malserver.external_product(0, &rlwe)?;
    // let pt_body = malserver.get_body_of_innerprod(&ext_rlwe);
    // println!("pt body with full body and masks: {}", pt_body);

    // let pt = malclient.decrypt(&ext_rlwe)?;
    // println!("pt direct dec over external product: {}", pt);


    Ok(())
}