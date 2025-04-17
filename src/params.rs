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

pub const DEFAULT_STATELESS_PARAMETER: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(6),
    polynomial_size: PolynomialSize(512),
    std_dev: 1.7347234759768072e-19,  // 3.2 / (2 ** 64)
    plaintext_modulus: (1 << 19) + 2,  // 2 * (1 << 18 + 1)
    delta: ((1u128 << 64) / ((1u128 << 19) + 2)) as u64,
    decomposition_base_log: DecompositionBaseLog(16),
    decomposition_level_count: DecompositionLevelCount(3),
    ciphertext_mask: 0xFFFFFFFF_FFFFFFFF,
    ciphertext_modulus: 1u128 << 64,
};

pub const DEFAULT_SELECTION_PARAMETER: GlweParameter<u32> = GlweParameter {
    glwe_size: GlweSize(4),
    polynomial_size: PolynomialSize(256),
    std_dev: 4.8828125e-05,  // 3.2 / (2 ** 16)
    plaintext_modulus: 2,
    delta: 1 << 31,
    decomposition_base_log: DecompositionBaseLog(8),  // unused
    decomposition_level_count: DecompositionLevelCount(3),  // unused
    ciphertext_mask: 0xFFFF_FFFF,
    ciphertext_modulus: 1u128 << 32,
};

/// RLWE parameter, modulo (2 ** 96)
pub const DEFAULT_MALICIOUS_PARAMETER: GlweParameter<u128> = GlweParameter {
    glwe_size: GlweSize(9),
    polynomial_size: PolynomialSize(512),  // 4096 dim for RLWE
    std_dev: 9.4039548065783e-39,  // =3.2 / (2 ** 128), because the Torus for CiphertextModulus::native_modulus() is in u128
    plaintext_modulus: (1u128 << 19) + 21, // 524309
    delta: 151109674856362064342866u128,  // (2**96) / ((2**19) + 21)
    decomposition_base_log: DecompositionBaseLog(32),
    decomposition_level_count: DecompositionLevelCount(2),
    ciphertext_mask: 0xFFFFFFFF_FFFFFFFF_FFFFFFFF,  // 2^96
    ciphertext_modulus: (1u128 << 96),
};
