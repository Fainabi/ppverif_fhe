# Privacy-preserving Facial Verification with FHE and hints
This repository provides the implementation about privacy-preserving facial verification for semi-honest two parties in one round.

## Requirements
Current requirement is just:
- Rust version >= 1.84 
to ensure the compilation of [tfhe.rs](https://github.com/zama-ai/tfhe-rs)

## Checkout for different architectures

For x86 system, run
```sh
$ git checkout x86_64
```

For android device, run
```sh
$ git checkout android
```

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


## JNI Support
Install Rust Android tools
```sh
$ rustup target add aarch64-linux-android armv7-linux-androideabi
$ cargo install cargo-ndk
```

Install the following build tools using Android Studio in  `Settings -> Languages & Frameworks -> Android SDK -> SDK Tools`
- NDK (Side by side)
- CMake
- Android SDK Command-line Tools

Modify `Cargo.toml`
```toml
[dependencies]
tfhe = { version = "0.6.4", default-features = false, features = ["boolean", "shortint", "integer", "seeder_unix"] }
jni = "0.21"

[lib]
crate-type = ["cdylib"]
```

Compile as dynamic library for ARM platform
```sh
$ cargo ndk -t armeabi-v7a -t arm64-v8a -o ../ppverif_android/app/src/main/jniLibs build --release
```
will generate `jniLibs/arm64-v8a/libppverif_fhe.so` and `jniLibs/armeabi-v7a/libppverif_fhe.so`

Call `testClient` using Kotlin:
```kotlin
package com.example.ppverif

object RustBridge {
    init {
        System.loadLibrary("ppverif_fhe")
    }

    @JvmStatic
    external fun testClient(): FloatArray
}

val testResult: FloatArray = RustBridge.testClient() // will get the test time
```
