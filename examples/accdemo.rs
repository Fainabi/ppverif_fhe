use fhe::Result;
use ppverif_fhe::*;
use rand::RngCore;
use std::os::windows::thread;
use std::time::Instant;
use fhe::bfv;
use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey, RGSWCiphertext, RelinearizationKey};
use rand::{rngs::OsRng, thread_rng};
use fhe_traits::*;
use fhe_math::rq::{Representation::PowerBasis, Representation, Poly, self};


fn main() -> Result<()> {
    let parameters = bfv::BfvParametersBuilder::new()
        .set_degree(4096)
        .set_moduli(&[
            0x3fffffff000001_u64, 
            0x3fffffff004001_u64, 
            0x3fffffff04e001, 
            0x3fffffff058001, 
            0x3fffffff07c001, 
            0x3fffffff0c6001, 
            0x3fffffff0cc001, 
            0x3fffffff0d2001
        ])
        .set_plaintext_modulus(0x88001)
        .build_arc()?;


    let mut client = Client::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER, parameters.clone())?;
    let pk = client.new_glwe_public_keys_br();
    let bfv_rlk = client.new_bfv_relinearizatio_key()?;
    let server = Server::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER, parameters.clone(), pk, bfv_rlk)?;

    let mut rng = thread_rng();
    let vals = (0..4096).map(|_| rng.next_u64()).collect::<Vec<_>>();
    let now = Instant::now();
    let cts = client.encrypt_rlwe_for_binary_decomposition(&vals)?;
    let elapsed = now.elapsed();

    // for group_ct in cts.iter() {
    //     // for ct in group_ct {
    //         let pt = client.decyrpt_bfv_ct(group_ct.last().unwrap())?;
    //         println!("pt: {:?}", pt);
    //     // }
    // }
    println!("enc time: {} ms", elapsed.as_micros() as f32 / 1000.0);

    let vs = (0..4096).map(|_| rng.next_u64()).collect::<Vec<_>>();
    let now = Instant::now();
    let folded = server.batch_fold(cts, &vs)?;
    let elapsed = now.elapsed();
    for ct in folded.iter() {
        let pt = client.decyrpt_bfv_ct(ct)?;
        println!("pt: {:?}", pt);
    }

    let sums = vals.iter().zip(vs.iter()).map(|(&v1, &v2)| (((v1 & 0b1) + (v2 & 0b1)) >> 1) & 0x1).collect::<Vec<_>>();
    println!("sums: {:?}", sums);
    println!("fold time: {}", elapsed.as_micros() as f32 / 1000.0);
    
    // let pt = Plaintext::try_encode(&vals, Encoding::simd(), &parameters)?;
    
    // let ctx = rq::Context::new_arc(parameters.moduli(), parameters.degree())?;

    // let mut poly = unsafe {
    //     Poly::create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(&vals, &ctx)
    // };
    // poly.change_representation(Representation::PowerBasis);
    // unsafe { poly.override_representation(Representation::Ntt); } 
    // // poly.change_representation(fhe_math::rq::Representation::PowerBasis);
    // // poly *= &poly.clone();
    // println!("poly: {:?}", poly);

    // let cts = client.encrypt_rlwe_for_binary_decomposition(&vals)?;

    // let poly = ppverif_fhe::new_rq_poly_ntt_from_slice(&vals, &ctx);
    // // let mut poly = Poly::zero(&ctx, Representation::Ntt);
    // poly.coefficients().rows().into_iter().for_each(|row_i| {
    //     println!("rowi: {:?}", row_i);
    // });

    Ok(())
}