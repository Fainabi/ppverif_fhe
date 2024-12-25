use tfhe::core_crypto::prelude::*;

pub fn extract_glwe_sample_from_rlwe_ciphertext<Scalar>(glwe_in: &GlweCiphertextOwned<Scalar>, polynomial_size: PolynomialSize) 
    -> GlweCiphertextOwned<Scalar> 
    where 
        Scalar: UnsignedInteger,
{
    // assume a RLWE ciphertext input
    let n = glwe_in.polynomial_size().0 / polynomial_size.0;
    let N = polynomial_size.0;
    let mut glwe_out_polys = vec![Scalar::ZERO; (n + 1) * N];

    let glwe_in_polys = glwe_in.as_polynomial_list();
    let glwe_in_mask = glwe_in_polys.get(0).into_container();
    let glwe_in_body = glwe_in_polys.get(1).into_container();

    for ni in 0..n {
        for Ni in 0..N {
            let offset = glwe_in_mask.len() - ni;

            glwe_out_polys[ni * N + Ni] = glwe_in_mask[(offset + Ni * n) % glwe_in_mask.len()];
            if Ni == 0 && ni != 0 {
                glwe_out_polys[ni * N + Ni] = Scalar::wrapping_sub(Scalar::ZERO, glwe_out_polys[ni * N + Ni]);
            }
        }
    }

    for Ni in 0..N {
        glwe_out_polys[n * N + Ni] = glwe_in_body[Ni * n];
    }

    GlweCiphertext::from_container(glwe_out_polys, polynomial_size, CiphertextModulus::new_native())
}
