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
    pub ciphertext_mask: Scalar,
    pub ciphertext_modulus: u128,
}

pub const DEFAULT_INNER_PRODUCT_PARAMETER: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(6),
    polynomial_size: PolynomialSize(512),
    std_dev: 1.7347234759768072e-19,  // 3.2 / (2 ** 64)
    plaintext_modulus: 1 << 19,
    delta: 1 << 45,
    decomposition_base_log: DecompositionBaseLog(32),
    decomposition_level_count: DecompositionLevelCount(1),
    ciphertext_mask: 0xFFFFFFFF_FFFFFFFF,
    ciphertext_modulus: 1u128 << 64,
};

pub const DEFAULT_BLIND_ROTATION_PARAMETER: GlweParameter<u16> = GlweParameter {
    glwe_size: GlweSize(4),
    polynomial_size: PolynomialSize(256),
    std_dev: 4.8828125e-05,  // 3.2 / (2 ** 16)
    plaintext_modulus: 2,
    delta: 1 << 15,
    decomposition_base_log: DecompositionBaseLog(16),  // unused
    decomposition_level_count: DecompositionLevelCount(1),  // unused
    ciphertext_mask: 0xFFFF,
    ciphertext_modulus: 1u128 << 16,
};
