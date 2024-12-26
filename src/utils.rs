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