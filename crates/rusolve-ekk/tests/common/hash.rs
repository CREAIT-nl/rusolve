//! The oracle's array hashes (FORMAT.md, "Hashes"), reimplemented exactly
//! from `tools/highs-oracle/src/emit.cpp` (`fnvF`, `fnvI`, `fnvI8`).
//!
//! FNV-1a 64 over little-endian bytes: each f64 as the 8 bytes of
//! `to_bits().to_le_bytes()`, each HighsInt (i32 in the oracle build) as 4
//! bytes, each int8 as 1 byte. Each array is hashed on its own, starting from
//! the offset basis; the empty array hashes to the offset basis.

pub const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
pub const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

/// Continue an FNV-1a 64 hash over `bytes`.
pub fn fnv1a_bytes(mut h: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

/// `fnvF`: hash of a `std::vector<double>` (bit patterns, so -0.0 != +0.0).
pub fn fnv_f64(v: &[f64]) -> u64 {
    v.iter().fold(FNV_OFFSET, |h, x| {
        fnv1a_bytes(h, &x.to_bits().to_le_bytes())
    })
}

/// `fnvI`: hash of a `std::vector<HighsInt>` (HighsInt = i32).
pub fn fnv_i32(v: &[i32]) -> u64 {
    v.iter()
        .fold(FNV_OFFSET, |h, x| fnv1a_bytes(h, &x.to_le_bytes()))
}

/// `fnvI8`: hash of a `std::vector<int8_t>` (e.g. `nonbasicMove_`).
pub fn fnv_i8(v: &[i8]) -> u64 {
    v.iter()
        .fold(FNV_OFFSET, |h, x| fnv1a_bytes(h, &x.to_le_bytes()))
}

/// Hashes are serialised as 16 lowercase hex digits.
pub fn hash_hex(h: u64) -> String {
    format!("{h:016x}")
}
