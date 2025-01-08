use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use ppverif_fhe::*;
use rand::{thread_rng, RngCore};
use std::error::Error;


fn malicious_benchmark(c: &mut Criterion) -> std::result::Result<(), Box<dyn Error>> {
    let mut client = MaliciousClient::new()?;
    let pk = client.new_public_key();
    let rlk = client.new_relin_keys()?;
    let mut server = AntiMalServer::new(pk, rlk)?;

    let mut rng = thread_rng();
    let mut features_to_verif = (0..512)
        .map(|_| (rng.next_u32() % 1024) as f32 / 1024.0)
        .collect::<Vec<_>>();
    let mut features_to_enroll = (0..512)
        .map(|_| (rng.next_u32() % 1024) as f32 / 1024.0)
        .collect::<Vec<_>>();

    normalize(&mut features_to_verif);
    normalize(&mut features_to_enroll);

    let (rgsw, rgsw_masks) = client.encrypt_rgsw_ciphertext(&features_to_enroll, 512.0)?;
    client.enroll_rgsw_masks(0, rgsw_masks);
    server.enroll_rgsw(0, rgsw);

    c.bench_function("Mal Encryption", |b| {
        b.iter(|| {
            client.encrypt_rlwe_ciphertext(&features_to_verif, 512.0).unwrap();
        });
    });

    let (d_packed, norm) = client.encrypt_rlwe_ciphertext(&features_to_verif, 512.0).unwrap();
    c.bench_function("client calculate inner prod", |b| b.iter(|| {
        client.compute_mask_for_innerprod(0, &d_packed).unwrap();
    }));

    // this includes the former one
    c.bench_function("Construct Constraints", |b| {
        b.iter(|| {
            client.encrypt_new_lookup_tables(0, 512, norm, &d_packed).unwrap();
        });
    });

    let (d_act, d_bin) = client.encrypt_new_lookup_tables(0, 512, norm, &d_packed).unwrap();

    c.bench_function("server calculate inner prod", |b| b.iter(|| {
        server.external_product(0, &d_packed).unwrap();
    }));
    c.bench_function("remask", |b| b.iter(|| {
        server.new_mask_with_pk_and_pt_divided_2().unwrap();
    }));
    c.bench_function("Mal Verification", |b| {
        b.iter(|| {
            server.verify_with_constraint(0, &d_packed, &d_act, &d_bin).unwrap();
        });
    });

    let d_cons = server.verify_with_constraint(0, &d_packed, &d_act, &d_bin)?;
    c.bench_function("Mal Decryption", |b| {
        b.iter(|| {
            client.decrypt(&d_cons).unwrap();
        });
    });

    Ok(())
}

criterion_group!(benches, malicious_benchmark);
criterion_main!(benches);