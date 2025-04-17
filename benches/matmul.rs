use criterion::{criterion_group, criterion_main, Criterion};
use rand::{thread_rng, RngCore};
use nalgebra;


pub fn matmul(c: &mut Criterion) {
    let row = 512 * 6;
    let col = 512 * 6;
    let mut rng = thread_rng();

    let dataA = (0..row*col).map(|_| 
        rng.next_u64()
    ).collect::<Vec<_>>();

    let A = nalgebra::DMatrix::from_row_slice(row, col, &dataA);
    let B = nalgebra::DMatrix::from_vec(col, 1, (0..col).map(|_| rng.next_u64()).collect::<Vec<_>>());
    c.bench_function("MatMul", |b| {
        b.iter(|| {
            let _ = &A * &B;
        });
    });
}

criterion_group!(benches, matmul);
criterion_main!(benches);