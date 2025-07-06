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

pub const INNER_PRODUCT_PARAMETER_DIMENSION_128: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(1 + 20),
    polynomial_size: PolynomialSize(128),
    std_dev: 1.7347234759768072e-19,  // 3.2 / (2 ** 64)
    plaintext_modulus: 1 << 19,
    delta: 1 << 45,
    decomposition_base_log: DecompositionBaseLog(32),
    decomposition_level_count: DecompositionLevelCount(1),
    ciphertext_mask: 0xFFFFFFFF_FFFFFFFF,
    ciphertext_modulus: 1u128 << 64,
};

pub const INNER_PRODUCT_PARAMETER_DIMENSION_256: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(1 + 10),
    polynomial_size: PolynomialSize(256),
    std_dev: 1.7347234759768072e-19,  // 3.2 / (2 ** 64)
    plaintext_modulus: 1 << 19,
    delta: 1 << 45,
    decomposition_base_log: DecompositionBaseLog(32),
    decomposition_level_count: DecompositionLevelCount(1),
    ciphertext_mask: 0xFFFFFFFF_FFFFFFFF,
    ciphertext_modulus: 1u128 << 64,
};

pub const INNER_PRODUCT_PARAMETER_DIMENSION_512: GlweParameter<u64> = DEFAULT_INNER_PRODUCT_PARAMETER;

pub const INNER_PRODUCT_PARAMETER_DIMENSION_1024: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(1 + 3),
    polynomial_size: PolynomialSize(1024),
    std_dev: 1.7347234759768072e-19,  // 3.2 / (2 ** 64)
    plaintext_modulus: 1 << 19,
    delta: 1 << 45,
    decomposition_base_log: DecompositionBaseLog(32),
    decomposition_level_count: DecompositionLevelCount(1),
    ciphertext_mask: 0xFFFFFFFF_FFFFFFFF,
    ciphertext_modulus: 1u128 << 64,
};

pub const INNER_PRODUCT_PARAMETER_DIMENSION_2048: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(1 + 2),
    polynomial_size: PolynomialSize(2048),
    std_dev: 1.7347234759768072e-19,  // 3.2 / (2 ** 64)
    plaintext_modulus: 1 << 19,
    delta: 1 << 45,
    decomposition_base_log: DecompositionBaseLog(32),
    decomposition_level_count: DecompositionLevelCount(1),
    ciphertext_mask: 0xFFFFFFFF_FFFFFFFF,
    ciphertext_modulus: 1u128 << 64,
};

pub const INNER_PRODUCT_PARAMETER_DIMENSION_4096: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(1 + 1),
    polynomial_size: PolynomialSize(4096),
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
