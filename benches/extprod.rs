use tfhe::core_crypto::prelude::*;
use criterion::{criterion_group, criterion_main, Criterion};
use ppverif_fhe::DEFAULT_SELECTION_PARAMETER;


pub fn extprod(c: &mut Criterion) {
    let param = DEFAULT_SELECTION_PARAMETER;
    let mut seeder = new_seeder();
        
    let mut secret_generator =
        SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());

    let sk = GlweSecretKey::generate_new_binary(
        GlweDimension(param.glwe_size.0 - 1),
        param.polynomial_size,
        &mut secret_generator,
    );
    let distribution = Gaussian::from_dispersion_parameter(StandardDev(param.std_dev), 0.0);
    let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
        seeder.seed(),
        seeder.as_mut(),
    );


    let mut ggsw = GgswCiphertext::new(0u32, param.glwe_size, param.polynomial_size, param.decomposition_base_log, param.decomposition_level_count, CiphertextModulus::new_native());
    encrypt_constant_ggsw_ciphertext(&sk, &mut ggsw, Plaintext(1), distribution, &mut encryption_generator);
    let mut ggsw_fft = FourierGgswCiphertext::new(param.glwe_size, param.polynomial_size, param.decomposition_base_log, param.decomposition_level_count);
    convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut ggsw_fft);

    let mut glwe_ct = GlweCiphertext::new(0, param.glwe_size, param.polynomial_size, CiphertextModulus::new_native());
    encrypt_glwe_ciphertext(&sk, &mut glwe_ct, &PlaintextList::new(0, PlaintextCount(param.polynomial_size.0)), distribution, &mut encryption_generator);

    let mut glwe_out = glwe_ct.clone();
    c.bench_function("ExtProd", |b| {
        b.iter(|| {
            for _ in 0..512 {
                add_external_product_assign(&mut glwe_out, &ggsw_fft, &glwe_ct);
            }
        });
    });
}

criterion_group!(benches, extprod);
criterion_main!(benches);
