use tfhe::core_crypto::prelude::*;

pub(crate) fn find_valid_squared_indices(poly_dim: usize) -> Vec<usize> {
    let mut indices = vec![];
    let mut valid_indices = (0..poly_dim/2).collect::<std::collections::BTreeSet<usize>>();
    for i in 0..poly_dim/2 {
        if valid_indices.contains(&i) {
            indices.push(i);
            for &j in &indices {
                valid_indices.remove(&j);
                valid_indices.remove(&(2 * i - j));
            }
        }
    }

    indices
}

pub(crate) fn polynomial_wrapping_mul_constant_term<Scalar, LhsCont, RhsCont>(
    lhs: &Polynomial<LhsCont>, 
    rhs: &Polynomial<RhsCont>
) -> Scalar
where 
    Scalar: UnsignedInteger,
    LhsCont: Container<Element = Scalar>,
    RhsCont: Container<Element = Scalar>,
{
    let lhs = lhs.as_ref();
    let rhs = rhs.as_ref();
    let mut res = Scalar::wrapping_mul(lhs[0], rhs[0]);
    for (&vi, &vj) in lhs[1..].iter().zip(rhs[1..].iter().rev()) {
        res = Scalar::wrapping_sub(res, Scalar::wrapping_mul(vi, vj));
    }
    res
}
