use criterion::{criterion_group, criterion_main, Criterion};
use ppverif_fhe::*;
use rand::*;

fn prng_dec_benchmark(c: &mut Criterion) {
    let mut client = Client::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER);

    let mut rng = thread_rng();
    let mut features_to_verif = (0..DEFAULT_INNER_PRODUCT_PARAMETER.polynomial_size.0)
        .map(|_| (rng.next_u32() % 1024) as f32 / 1024.0)
        .collect::<Vec<_>>();
    let mut features_to_enroll = (0..DEFAULT_INNER_PRODUCT_PARAMETER.polynomial_size.0)
        .map(|_| (rng.next_u32() % 1024) as f32 / 1024.0)
        .collect::<Vec<_>>();

    normalize(&mut features_to_enroll);
    normalize(&mut features_to_verif);

    c.bench_function("New GGSW", |b| {
        b.iter(|| {
            client.encrypt_new_template(0, &features_to_enroll, 512.0);
            client.enroll_ggsw_masks(0);
        });
    });

    
    c.bench_function("PRNG + dec", |b| {
        b.iter(|| {
            client.enroll_ggsw_masks(0);
        });
    });

    
}

criterion_group!(benches, prng_dec_benchmark);
criterion_main!(benches);