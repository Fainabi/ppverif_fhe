# Privacy-preserving Facial Verification with FHE and hints
This repository provides the implementation about privacy-preserving facial verification for semi-honest two parties in one round.

## Requirements
Current requirement is just:
- Rust version >= 1.84 
to ensure the compilation of [tfhe.rs](https://github.com/zama-ai/tfhe-rs)

## Demo and Benches
To run the demo, enter:
```sh
$ cargo run --example demo --release
```

For benchmark on the case that the client stores decrypted hints, run:
```sh
$ cargo bench --bench hint_bench
```

For benchnmark on the case that the client stores the PRNG seeds run:
```sh
$ cargo bench --bench prng_bench
```
This bench will test the time for prng generation and hint construction. Other processes are exactly same to that in `hint_bench` so one just check the `hint_bench` for the benchmark test.
