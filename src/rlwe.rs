use tfhe::core_crypto::prelude::*;
use polynomial_algorithms::{polynomial_karatsuba_wrapping_mul, polynomial_wrapping_add_assign, polynomial_wrapping_add_mul_assign, polynomial_wrapping_mul};

#[allow(non_snake_case)]
pub fn extract_glwe_sample_from_rlwe_ciphertext<Scalar, Cont>(glwe_in: &GlweCiphertext<Cont>, polynomial_size: PolynomialSize) 
    -> GlweCiphertextOwned<Scalar> 
    where 
        Scalar: UnsignedInteger,
        Cont: Container<Element = Scalar>,
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

pub fn extract_lwe_sample_from_glwe_ciphertext_under_rlwe_secret_key<Scalar, Cont>(glwe_in: &GlweCiphertext<Cont>)
    -> LweCiphertextOwned<Scalar>
    where
        Scalar: UnsignedInteger,
        Cont: Container<Element = Scalar>
{
    let n = glwe_in.glwe_size().0 - 1;
    let N = glwe_in.polynomial_size().0;
    let mask_dim = N * n;
    let mut lwe_out = LweCiphertext::new(
        Scalar::ZERO, 
        LweSize(mask_dim + 1), 
        CiphertextModulus::new_native()
    );

    let mut lwe_out_mask = lwe_out.get_mut_mask();
    for (i, a_i) in glwe_in.get_mask().as_polynomial_list().iter().enumerate() {
        for (j, v) in a_i.as_ref().iter().enumerate() {
            lwe_out_mask.as_mut()[(n * ((N - j) % N) + i) % mask_dim] = 
                if j == 0 { *v } else { Scalar::wrapping_neg(*v) }; 
        }
    }
    *lwe_out.get_mut_body().data = glwe_in.get_body().as_ref()[0];

    lwe_out
}


// assume `t * lhs * rhs` does not overflow w.r.t u256
pub fn rlwe_multiplication_u96(lhs: &GlweCiphertextOwned<u128>, rhs: &GlweCiphertextOwned<u128>, t: u128) -> GlweCiphertextOwned<u128>
{
    let mut ct_out = GlweCiphertext::new(0u128, GlweSize(3), lhs.polynomial_size(), CiphertextModulus::new_native());
    let mut ct_out_poly = ct_out.as_mut_polynomial_list();

    // decompose lhs: u96 to (lhs_lower: u32, lhs_upper: u64)
    //           rhs: u96 to (rhs_lower: u32, rhs_upper: u64)
    let lhs_poly = lhs.as_polynomial_list();
    let rhs_poly = rhs.as_polynomial_list();
    let lhs_poly_lower = lhs_poly
        .iter()
        .map(|poly| {
            let containter = poly
                .into_container()
                .iter()
                .map(|&v| v & 0xFFFFFFFF)
                .collect::<Vec<_>>();
            Polynomial::from_container(containter)
        })
        .collect::<Vec<_>>();
    let lhs_poly_upper = lhs_poly
        .iter()
        .map(|poly| {
            let containter = poly
                .into_container()
                .iter()
                .map(|&v| (v >> 32) & 0xFFFFFFFF_FFFFFFFF)
                .collect::<Vec<_>>();
            Polynomial::from_container(containter)
        })
        .collect::<Vec<_>>();
    let rhs_poly_lower = rhs_poly
        .iter()
        .map(|poly| {
            let containter = poly
                .into_container()
                .iter()
                .map(|&v| v & 0xFFFFFFFF)
                .collect::<Vec<_>>();
            Polynomial::from_container(containter)
        })
        .collect::<Vec<_>>();
    let rhs_poly_upper = rhs_poly
        .iter()
        .map(|poly| {
            let containter = poly
                .into_container()
                .iter()
                .map(|&v| (v >> 32) & 0xFFFFFFFF_FFFFFFFF)
                .collect::<Vec<_>>();
            Polynomial::from_container(containter)
        })
        .collect::<Vec<_>>();
    
    // calculate:
    //     (lhs_lower * rhs_upper / 2**64) * t
    //   + (rhs_lower * lhs_upper / 2**64) * t
    //   + (lhs_upper * rhs_upper / 2**32) * t
    polynomial_wrapping_add_mul_assign_u96(&mut ct_out_poly.get_mut(2), &lhs_poly_upper[1], &lhs_poly_lower[1], &rhs_poly_upper[1], &rhs_poly_lower[1]);
    polynomial_wrapping_add_mul_assign_u96(&mut ct_out_poly.get_mut(1), &lhs_poly_upper[0], &lhs_poly_lower[0], &rhs_poly_upper[1], &rhs_poly_lower[1]);
    polynomial_wrapping_add_mul_assign_u96(&mut ct_out_poly.get_mut(1), &lhs_poly_upper[1], &lhs_poly_lower[1], &rhs_poly_upper[0], &rhs_poly_lower[0]);
    polynomial_wrapping_add_mul_assign_u96(&mut ct_out_poly.get_mut(0), &lhs_poly_upper[0], &lhs_poly_lower[0], &rhs_poly_upper[0], &rhs_poly_lower[0]);
    ct_out.as_mut().iter_mut().for_each(|v| *v = (*v * t) >> 32);

    ct_out
}

/// calculate:
///     (lhs_lower * rhs_upper / 2**32)
///   + (rhs_lower * lhs_upper / 2**32)
///   + (lhs_upper * rhs_upper)
fn polynomial_wrapping_add_mul_assign_u96<OutputCont, LhsContUpper, LhsContLower, RhsContUpper, RhsContLower>(
    output: &mut Polynomial<OutputCont>, 
    lhs_upper: &Polynomial<LhsContUpper>, 
    lhs_lower: &Polynomial<LhsContLower>, 
    rhs_upper: &Polynomial<RhsContUpper>,
    rhs_lower: &Polynomial<RhsContLower>
) where 
    OutputCont: ContainerMut<Element = u128>,
    LhsContUpper: Container<Element = u128>,
    RhsContUpper: Container<Element = u128>,
    LhsContLower: Container<Element = u128>,
    RhsContLower: Container<Element = u128>,
{
    let mut tmp = Polynomial::new(0, output.polynomial_size());
    
    // TODO: reduce to only two multiplication, and using karatsuba
    polynomial_karatsuba_wrapping_mul(&mut tmp, lhs_lower, rhs_upper);
    polynomial_wrapping_add_mul_assign(&mut tmp, rhs_lower, lhs_upper);
    tmp.as_mut().iter_mut().for_each(|v| *v = signed_shr(*v, 32));
    
    polynomial_wrapping_add_mul_assign(&mut tmp, lhs_upper, rhs_upper);
    polynomial_wrapping_add_assign(output,&tmp);
}

fn signed_shr(v: u128, n: usize) -> u128 {
    if v >> 127 == 0x1 {
        u128::wrapping_neg(u128::wrapping_neg(v) >> n)
    } else {
        v >> n
    }
}


