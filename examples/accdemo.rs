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
            // 0x100006001,
            // 0x100014001,
            // 0x100036001,
            // 0x100050001,
            // 0x10005c001,
            // 0x100086001,
            // 0x100098001,
            // 0x1000b0001,
            // 0x1000e6001,
            // 0x10010a001,
            // 0x100132001,
            // 0x100174001,
            // 0x100176001,
            // 0x100180001,
            // 0x1001c4001,
            // 0x1001ce001,
            // 0x1001d6001,
            // 0x1001ee001,
            // 0x100230001,
            // 0x10023a001,
            // 0x10023c001,
            // 0x100248001,
            // 0x100252001,
            // 0x10025a001,
            // 0x100270001,
            // 0x1002a8001,
            // 0x1002b2001,
            // 0x1002fa001,
            // 0x10030c001,
            // 0x100324001,
            // 0x100330001,
            // 0x100348001,
            // 0x100360001,
            // 0x100368001,
            // 0x10037e001,
            0x3fffffff000001, 
            0x3fffffff004001, 
            0x3fffffff04e001, 
            0x3fffffff058001, 
            0x3fffffff07c001, 
            0x3fffffff0c6001, 
            0x3fffffff0cc001, 
            0x3fffffff0d2001,
            0x3fffffff124001,
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

    let vs = (0..4096).map(|i| rng.next_u64()).collect::<Vec<_>>();
    let now = Instant::now();
    let folded = server.batch_fold(cts, &vs)?;
    let elapsed = now.elapsed();
    let sums = vals.iter().zip(vs.iter()).map(|(&v1, &v2)| (((v1 & 0b11111111) + (v2 & 0b11111111)) >> 7) & 0x1).collect::<Vec<_>>();
    for ct in folded.iter() {
        let pt = client.decyrpt_bfv_ct(ct)?;
        let diff = pt.iter().zip(sums.iter()).map(|(&v1, &v2)| v1 - v2).collect::<Vec<_>>();
        println!("diff: {:?}", diff);
        // println!("pt: {:?}", pt);
        // println!("sums: {:?}", sums);
    }

    
    // println!("sums: {:?}", sums);
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