use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use ppverif_fhe::*;
use rand::{thread_rng, RngCore};
use tfhe::core_crypto::prelude::*;


fn pk_enc_benchmark(c: &mut Criterion) {
    let mut client = Client::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER);
    let glwe_pk = client.new_glwe_public_keys_br();
    let mut server= Server::new(DEFAULT_INNER_PRODUCT_PARAMETER, DEFAULT_BLIND_ROTATION_PARAMETER, glwe_pk);

    let mut rng = thread_rng();
    let mut features_to_verif = (0..DEFAULT_INNER_PRODUCT_PARAMETER.polynomial_size.0)
        .map(|_| (rng.next_u32() % 1024) as f32 / 1024.0)
        .collect::<Vec<_>>();
    let mut features_to_enroll = (0..DEFAULT_INNER_PRODUCT_PARAMETER.polynomial_size.0)
        .map(|_| (rng.next_u32() % 1024) as f32 / 1024.0)
        .collect::<Vec<_>>();

    normalize(&mut features_to_enroll);
    normalize(&mut features_to_verif);

    let ggsw = client.encrypt_new_template(0, &features_to_enroll, 512.0);
    client.enroll_ggsw_masks(0);
    server.enroll(0, ggsw);

    c.bench_function("New Template Encryption", |b| {
        b.iter(|| {
            client.encrypt_glwe(&features_to_verif, 512.0);
        });
    });

    let query_ct = client.encrypt_glwe(&features_to_verif, 512.0);
    c.bench_function("Generate LUT", |b| {
        b.iter(|| {
            client.transform_mask_from_database(0, query_ct.as_view()).unwrap();
        });
    });

    c.bench_function("Client calculate inner prod", |b| b.iter(|| {
        client.transform_mask_to_body_from_database(0, query_ct.as_view());
    }));


    let lut_ct = client.transform_mask_from_database(0, query_ct.as_view()).unwrap();

    c.bench_function("Server calculate inner prod", |b| b.iter(|| {
        server.compute_innerprod_body(0, query_ct.as_view());
    }));

    let mut tmp_lwe = LweCiphertext::new(
        0, 
        LweSize((DEFAULT_BLIND_ROTATION_PARAMETER.glwe_size.0 - 1) * DEFAULT_BLIND_ROTATION_PARAMETER.polynomial_size.0 + 1), 
        CiphertextModulus::new_native()
    );
    extract_lwe_sample_from_glwe_ciphertext(&lut_ct.get(0), &mut tmp_lwe, MonomialDegree(0));
    c.bench_function("Sample Extraction", |b| b.iter(|| {
        extract_lwe_sample_from_glwe_ciphertext(&lut_ct.get(0), &mut tmp_lwe, MonomialDegree(0));
    }));

    c.bench_function("Verification", |b| {
        b.iter(|| {
            server.verify(0, &query_ct, &lut_ct).unwrap();
        });
    });

    let ct_verif = server.verify(0, &query_ct, &lut_ct).unwrap();
    c.bench_function("Decryption", |b| {
        b.iter(|| {
            client.decrypt_lwe(&ct_verif);
        });
    });
}

criterion_group!(benches, pk_enc_benchmark);
criterion_main!(benches);