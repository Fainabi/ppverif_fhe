# Privacy-preserving biometric verification

Run semi-homest example:
```sh
$ cargo run --release --example demo
```

Run plaintext-malicious example:
```sh
$ cargo run --release --example mal_demo
```
The least significant bit of decryption is the verification result.

Run benchmark
```sh
$ cargo bench
```

