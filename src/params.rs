use aligned_vec::ABox;
use concrete_fft::c64;
use tfhe::core_crypto::prelude::*;

#[derive(Debug, Copy, Clone)]
pub struct GlweParameter<Scalar> {
    pub glwe_size: GlweSize,
    pub polynomial_size: PolynomialSize,
    pub std_dev: f64,
    pub plaintext_modulus: Scalar,
    pub delta: Scalar,
    pub decomposition_base_log: DecompositionBaseLog,
    pub decomposition_level_count: DecompositionLevelCount,
}

pub const DEFAULT_INNER_PRODUCT_PARAMETER: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(6),
    polynomial_size: PolynomialSize(512),
    std_dev: 1.7347234759768072e-13,
    // std_dev: 0.0,
    plaintext_modulus: 1 << 19,
    delta: 1 << 45,
    decomposition_base_log: DecompositionBaseLog(40),  // = 64 - 24
    decomposition_level_count: DecompositionLevelCount(1),
};

pub const DEFAULT_BLIND_ROTATION_PARAMETER: GlweParameter<u16> = GlweParameter {
    glwe_size: GlweSize(4),
    polynomial_size: PolynomialSize(256),
    // std_dev: 0.00000000000000029403601535432533,  // chi = 3.2
    std_dev: 0.0,
    plaintext_modulus: 2,
    delta: 1 << 15,
    decomposition_base_log: DecompositionBaseLog(16),  // unused
    decomposition_level_count: DecompositionLevelCount(1),  // unused
};

/// RLWE parameter, modulo (2 ** 96)
pub const DEFAULT_MALICIOUS_PARAMETER: GlweParameter<u128> = GlweParameter {
    glwe_size: GlweSize(2),
    polynomial_size: PolynomialSize(4096),
    std_dev: 0.00000000000000029403601535432533,
    plaintext_modulus: 1 << 19,
    delta: 1 << 45,
    decomposition_base_log: DecompositionBaseLog(16),
    decomposition_level_count: DecompositionLevelCount(1),
};

pub type FourierGgswCiphertextOwned = FourierGgswCiphertext<ABox<[c64]>>;
