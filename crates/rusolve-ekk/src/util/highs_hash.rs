//! HiGHS source: util/HighsHash.h (HighsHashHelpers subset)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Only the helpers the simplex path needs: `pair_hash` (for `HighsRandom`),
//! `log2i`, the M61 modular arithmetic, and `sparse_combine` /
//! `sparse_inverse_combine` (basis hashing in `HEkk`). All `u64` arithmetic
//! in C++ wraps modulo 2^64: use `wrapping_add` / `wrapping_mul` /
//! `wrapping_sub` everywhere a C++ `+`, `*`, `-` acts on `u64`.
//! `HighsHashTable<uint64_t>` (the visited-basis set) is `std::collections::HashSet<u64>`.
//!
//! Presolve additions: `fibonacci_muliplier`, the small-key `hash` overload
//! (`hash_u64_bytes`, `hash_i32`, `hash_u32_pair`), `golden_ratio_reciprocal`,
//! `double_hash_code` and a bit-exact `frexp` (no libm in this crate). The
//! generic `HighsHashTable<K, V>` lives in `util/highs_hash_table.rs`.

use crate::HighsInt;

/// HiGHS: HighsHash.h::HighsHashHelpers::c (64 random constants)
pub const c: [u64; 64] = [
    0xc8497d2a400d9551,
    0x80c8963be3e4c2f3,
    0x042d8680e260ae5b,
    0x8a183895eeac1536,
    0xa94e9c75f80ad6de,
    0x7e92251dec62835e,
    0x07294165cb671455,
    0x89b0f6212b0a4292,
    0x31900011b96bf554,
    0xa44540f8eee2094f,
    0xce7ffd372e4c64fc,
    0x51c9d471bfe6a10f,
    0x758c2a674483826f,
    0xf91a20abe63f8b02,
    0xc2a069024a1fcc6f,
    0xd5bb18b70c5dbd59,
    0xd510adac6d1ae289,
    0x571d069b23050a79,
    0x60873b8872933e06,
    0x780481cc19670350,
    0x7a48551760216885,
    0xb5d68b918231e6ca,
    0xa7e5571699aa5274,
    0x7b6d309b2cfdcf01,
    0x04e77c3d474daeff,
    0x4dbf099fd7247031,
    0x5d70dca901130beb,
    0x9f8b5f0df4182499,
    0x293a74c9686092da,
    0xd09bdab6840f52b3,
    0xc05d47f3ab302263,
    0x6b79e62b884b65d6,
    0xa581106fc980c34d,
    0xf081b7145ea2293e,
    0xfb27243dd7c3f5ad,
    0x5211bf8860ea667f,
    0x9455e65cb2385e7f,
    0x0dfaf6731b449b33,
    0x4ec98b3c6f5e68c7,
    0x007bfd4a42ae936b,
    0x65c93061f8674518,
    0x640816f17127c5d1,
    0x6dd4bab17b7c3a74,
    0x34d9268c256fa1ba,
    0x0b4d0c6b5b50d7f4,
    0x30aa965bc9fadaff,
    0xc0ac1d0c2771404d,
    0xc5e64509abb76ef2,
    0xd606b11990624a36,
    0x0d3f05d242ce2fb7,
    0x469a803cb276fe32,
    0xa4a44d177a3e23f4,
    0xb9d9a120dcc1ca03,
    0x2e15af8165234a2e,
    0x10609ba2720573d4,
    0xaa4191b60368d1d5,
    0x333dd2300bc57762,
    0xdf6ec48f79fb402f,
    0x5ed20fcef1b734fa,
    0x4c94924ec8be21ee,
    0x5abe6ad9d131e631,
    0xbe10136a522e602d,
    0x53671115c340e779,
    0x9f392fe43e2144da,
];

/// HiGHS: HighsHash.h::HighsHashHelpers::M61 (Mersenne prime 2^61 - 1)
pub const fn M61() -> u64 {
    0x1fffffffffffffff
}

/// HiGHS: HighsHash.h::HighsHashHelpers::log2i(uint64_t) (builtin-clz variant, line 102)
pub fn log2i_u64(n: u64) -> i32 {
    63 - n.leading_zeros() as i32
}

/// HiGHS: HighsHash.h::HighsHashHelpers::log2i(uint32_t) (builtin-clz variant, line 104)
pub fn log2i_u32(n: u32) -> i32 {
    31 - n.leading_zeros() as i32
}

/// HiGHS: HighsHash.h::HighsHashHelpers::multiply_modM61 (lines 193-224)
pub fn multiply_mod_m61(a: u64, b: u64) -> u64 {
    let ahi = a >> 32;
    let bhi = b >> 32;
    let alo = a & 0xffff_ffffu64;
    let blo = b & 0xffff_ffffu64;

    // compute the different order terms with adicities 2^64, 2^32, 2^0
    let term_64 = ahi.wrapping_mul(bhi);
    let term_32 = ahi.wrapping_mul(blo).wrapping_add(bhi.wrapping_mul(alo));
    let mut term_0 = alo.wrapping_mul(blo);

    // Partially reduce term_0 and term_32 modulo M61() individually.
    term_0 = (term_0 & M61()).wrapping_add(term_0 >> 61);
    term_0 = term_0.wrapping_add(((term_32 >> 29).wrapping_add(term_32 << 32)) & M61());

    // The lower 61 bits of term_0 are now the lower 61 bits of the result.
    let ab61 = (term_64 << 3) | (term_0 >> 61);

    // finally take the result modulo M61
    let mut result = (term_0 & M61()).wrapping_add(ab61);
    if result >= M61() {
        result -= M61();
    }
    result
}

/// HiGHS: HighsHash.h::HighsHashHelpers::modexp_M61 (lines 226-240)
pub fn modexp_m61(a: u64, e: u64) -> u64 {
    // the exponent need to be greater than zero
    debug_assert!(e > 0);
    let mut result = a;
    let mut e = e;

    while e != 1 {
        // square
        result = multiply_mod_m61(result, result);

        // multiply with a if exponent is odd
        if (e & 1) != 0 {
            result = multiply_mod_m61(result, a);
        }

        // shift to next bit
        e >>= 1;
    }

    result
}

/// HiGHS: HighsHash.h::HighsHashHelpers::pair_hash<k> (lines 275-279).
/// The template parameter `k` becomes the first argument:
/// `pair_hash<k>(a, b)` -> `pair_hash(k, a, b)`. Arguments are `u32`; C++
/// call sites that pass a `u64` (e.g. `state >> 32`) truncate implicitly,
/// so write `(state >> 32) as u32`. The product wraps modulo 2^64.
pub fn pair_hash(k: usize, a: u32, b: u32) -> u64 {
    (u64::from(a).wrapping_add(c[2 * k])).wrapping_mul(u64::from(b).wrapping_add(c[2 * k + 1]))
}

/// HiGHS: HighsHash.h::HighsHashHelpers::sparse_combine(u64&, HighsInt, u64) (lines 281-310)
pub fn sparse_combine_value(hash: &mut u64, index: HighsInt, value: u64) {
    // make sure input value is never zero and at most 61bits are used
    let value = ((value << 1) & M61()) | 1;

    // make sure that the constant has at most 61 bits
    let a = c[(index & 63) as usize] & M61();
    let degree = ((index as u64) >> 6) + 1;

    *hash = hash.wrapping_add(multiply_mod_m61(value, modexp_m61(a, degree)));
    *hash = (*hash >> 61) + (*hash & M61());
    if *hash >= M61() {
        *hash -= M61();
    }
    debug_assert!(*hash < M61());
}

/// HiGHS: HighsHash.h::HighsHashHelpers::sparse_inverse_combine(u64&, HighsInt, u64) (lines 312-334)
pub fn sparse_inverse_combine_value(hash: &mut u64, index: HighsInt, value: u64) {
    // make sure input value is never zero and at most 61bits are used
    let value = ((value << 1) & M61()) | 1;

    let a = c[(index & 63) as usize] & M61();
    let degree = ((index as u64) >> 6) + 1;
    // add the additive inverse (M61() - hashvalue) instead of the hash value itself
    *hash = hash.wrapping_add(M61().wrapping_sub(multiply_mod_m61(value, modexp_m61(a, degree))));
    *hash = (*hash >> 61) + (*hash & M61());
    if *hash >= M61() {
        *hash -= M61();
    }
    debug_assert!(*hash < M61());
}

/// HiGHS: HighsHash.h::HighsHashHelpers::sparse_combine(u64&, HighsInt) (lines 336-346)
pub fn sparse_combine(hash: &mut u64, index: HighsInt) {
    let a = c[(index & 63) as usize] & M61();
    let degree = ((index as u64) >> 6) + 1;

    *hash = hash.wrapping_add(modexp_m61(a, degree));
    *hash = (*hash >> 61) + (*hash & M61());
    if *hash >= M61() {
        *hash -= M61();
    }
    debug_assert!(*hash < M61());
}

/// HiGHS: HighsHash.h::HighsHashHelpers::sparse_inverse_combine(u64&, HighsInt) (lines 348-365)
pub fn sparse_inverse_combine(hash: &mut u64, index: HighsInt) {
    let a = c[(index & 63) as usize] & M61();
    let degree = ((index as u64) >> 6) + 1;
    // add the additive inverse (M61() - hashvalue) instead of the hash value itself
    *hash = hash.wrapping_add(M61().wrapping_sub(modexp_m61(a, degree)));
    *hash = (*hash >> 61) + (*hash & M61());
    if *hash >= M61() {
        *hash -= M61();
    }
    debug_assert!(*hash < M61());
}

/// HiGHS: HighsHash.h::HighsHashHelpers::M31 (line 245): mersenne prime 2^31 - 1
pub const fn M31() -> u64 {
    0x7fff_ffff
}

/// HiGHS: HighsHash.h::HighsHashHelpers::multiply_modM31 (lines 248-253)
pub fn multiply_mod_m31(a: u32, b: u32) -> u32 {
    let mut result: u64 = u64::from(a) * u64::from(b);
    result = (result >> 31) + (result & M31());
    if result >= M31() {
        result -= M31();
    }
    result as u32
}

/// HiGHS: HighsHash.h::HighsHashHelpers::modexp_M31 (lines 255-272)
pub fn modexp_m31(a: u32, e: u64) -> u32 {
    // the exponent need to be greater than zero
    debug_assert!(e > 0);
    let mut result: u32 = a;
    let mut e: u64 = e;
    while e != 1 {
        // square
        result = multiply_mod_m31(result, result);
        // multiply with a if exponent is odd
        if (e & 1) != 0 {
            result = multiply_mod_m31(result, a);
        }
        // shift to next bit
        e >>= 1;
    }
    result
}

/// HiGHS: HighsHash.h::HighsHashHelpers::sparse_combine32(u32&, HighsInt, u64) (lines 367-398)
pub fn sparse_combine32(hash: &mut u32, index: HighsInt, value: u64) {
    // make sure input value is never zero and at most 31bits are used
    let value: u64 = (pair_hash(0, value as u32, (value >> 32) as u32) >> 33) | 1;

    // make sure that the constant has at most 31 bits, as otherwise the modulo
    // algorithm for multiplication mod M31 might not work properly due to
    // overflow
    let a: u32 = (c[(index & 63) as usize] & M31()) as u32;
    let degree: u64 = ((index as u64) >> 6) + 1;

    let mut result: u64 = u64::from(*hash);
    result += u64::from(multiply_mod_m31(value as u32, modexp_m31(a, degree)));
    result = (result >> 31) + (result & M31());
    if result >= M31() {
        result -= M31();
    }
    debug_assert!(result < M31());
    *hash = result as u32;
}

/// HiGHS: HighsHash.h::HighsHashHelpers::fibonacci_muliplier (line 425)
pub const fn fibonacci_muliplier() -> u64 {
    0x9e3779b97f4a7c15
}

/// HiGHS: HighsHash.h::HighsHashHelpers::hash(const T&) for `sizeof(T) <= 8`
/// (lines 565-576). `bytes[0]`, `bytes[1]` are the two little-endian `u32`
/// halves of the key (zero-padded). Note the C++ precedence:
/// `pair_hash<1>(..) ^ (pair_hash<0>(..) >> 32)`.
pub fn hash_u64_bytes(b0: u32, b1: u32) -> u64 {
    pair_hash(1, b0, b1) ^ (pair_hash(0, b0, b1) >> 32)
}

/// `HighsHashHelpers::hash(HighsInt)`: a 4-byte key, `bytes[1] = 0`.
pub fn hash_i32(v: HighsInt) -> u64 {
    hash_u64_bytes(v as u32, 0)
}

/// `HighsHashHelpers::hash(std::make_pair(uint32_t a, uint32_t b))`: an
/// 8-byte key, `bytes[0] = a`, `bytes[1] = b`.
pub fn hash_u32_pair(a: u32, b: u32) -> u64 {
    hash_u64_bytes(a, b)
}

/// HiGHS: HighsHash.h::HighsHashHelpers::golden_ratio_reciprocal (lines 732-734)
pub const fn golden_ratio_reciprocal() -> f64 {
    0.61803398874989484
}

/// C `std::frexp`: returns `(mantissa, exponent)` with `|mantissa|` in
/// `[0.5, 1)` and `x = mantissa * 2^exponent`; zero, infinities and NaN are
/// returned unchanged with exponent 0 (C99 7.12.6.4). Bit-exact.
pub fn frexp(x: f64) -> (f64, i32) {
    if x == 0.0 || !x.is_finite() {
        return (x, 0);
    }
    let mut bits = x.to_bits();
    let mut e_adj = 0;
    if (bits >> 52) & 0x7ff == 0 {
        // subnormal: scale by 2^54 (exact) and correct the exponent
        let scaled = x * f64::from_bits((1023u64 + 54) << 52);
        bits = scaled.to_bits();
        e_adj = -54;
    }
    let exp_field = ((bits >> 52) & 0x7ff) as i32;
    let e = exp_field - 1022 + e_adj;
    let mantissa = f64::from_bits((bits & !(0x7ffu64 << 52)) | (1022u64 << 52));
    (mantissa, e)
}

/// HiGHS: HighsHash.h::HighsHashHelpers::double_hash_code (lines 736-775)
///
/// `std::ldexp(hashbits, 15)` is the exact `hashbits * 32768.0`;
/// `(i16)` truncates toward zero (`as i16`, in range since `|hashbits| < 1`).
pub fn double_hash_code(val: f64) -> u32 {
    let (hashbits, exponent) = frexp(val * golden_ratio_reciprocal());
    let mut hashvalue = exponent as i16 as u16 as u32;
    hashvalue = (hashvalue << 16) | ((hashbits * 32768.0) as i16 as u16 as u32);
    hashvalue
}

/// Element types that `vector_hash` accepts. The C++ template takes any
/// trivially copyable `T` (`HighsHashable<T>`); the Rust port is instantiated
/// for `HighsInt` (`i32`) and `u32`, which are the element types the C++
/// call sites use. Both are 4 bytes, so each element is one `u32` word: its
/// native-endian bytes read back as a `u32` are its bit pattern.
pub trait HighsHashable: Copy {
    /// The 4 native-endian bytes of `self`, read as a `u32`.
    fn ne_word(self) -> u32;
}

impl HighsHashable for i32 {
    fn ne_word(self) -> u32 {
        self as u32
    }
}

impl HighsHashable for u32 {
    fn ne_word(self) -> u32 {
        self
    }
}

/// HiGHS: HighsHash.h::HighsHashHelpers::vector_hash (lines 428-557)
///
/// Hashes the raw bytes of `vals`. The C++ routine reads its input in chunks
/// of at most 256 bytes and processes each chunk in 8-byte pairs. Its `pair`
/// buffer persists across chunks, so a short final pair keeps the stale bytes
/// of the previous pair in the bytes it does not overwrite; that behavior is
/// reproduced in `vector_hash_words`.
pub fn vector_hash<T: HighsHashable>(vals: &[T]) -> u64 {
    vector_hash_words(vals.len(), vals.iter().map(|&v| v.ne_word()))
}

/// Core of `HighsHashHelpers::vector_hash` (HighsHash.h lines 428-557) over
/// `len` 4-byte elements given as `u32` words, pulled from `words` in order.
///
/// The C++ hashes the object bytes in place. Every 8-byte pair is two words, so
/// working on words gives the same pairs without first copying the bytes into a
/// buffer. The words may come from a lazy iterator, so callers can hash a
/// computed sequence without materialising it.
pub fn vector_hash_words(len: usize, mut words: impl Iterator<Item = u32>) -> u64 {
    // persistent across chunks, like the C++ `std::array<u32, 2> pair{}`
    let mut pair = [0u32; 2];
    let mut hash: u64 = 0;
    let mut k: usize = 0;
    // words still to be consumed; a chunk is at most 256 bytes = 64 words
    let mut remaining: usize = len;

    while remaining != 0 {
        let num_bytes = (remaining * 4).min(256);
        let num_pairs = num_bytes.div_ceil(8);
        let last_pair_bytes = num_bytes - (num_pairs - 1) * 8;
        let mut chunkhash = [0u64; 2];

        // The C++ switch falls through from case numPairs down to case 1, so
        // pairs are consumed in order; the pair consumed in step `j` is case
        // `n = num_pairs - j`, and only the final pair (j = num_pairs - 1)
        // copies `last_pair_bytes`.
        if num_pairs == 32 && hash != 0 {
            // make sure hash is reduced mod M61() before multiplying with the
            // next random constant
            if hash >= M61() {
                hash -= M61();
            }
            hash = multiply_mod_m61(hash, c[k & 63] & M61());
            k += 1;
        }
        for j in 0..num_pairs {
            let n = num_pairs - j;
            // memcpy(&pair[0], dataptr, B): a final 4-byte pair overwrites only
            // pair[0] and keeps the previous pair[1]; every other pair is full.
            if j == num_pairs - 1 && last_pair_bytes == 4 {
                pair[0] = words
                    .next()
                    .expect("vector_hash_words: fewer words than len");
            } else {
                pair[0] = words
                    .next()
                    .expect("vector_hash_words: fewer words than len");
                pair[1] = words
                    .next()
                    .expect("vector_hash_words: fewer words than len");
            }
            // chunkhash[N & 1] += pair_hash<32 - N>(pair[0], pair[1])
            chunkhash[n & 1] = chunkhash[n & 1].wrapping_add(pair_hash(32 - n, pair[0], pair[1]));
        }
        remaining -= num_bytes / 4;

        hash = hash.wrapping_add((chunkhash[0] >> 3) ^ (chunkhash[1] >> 32));
    }

    hash.wrapping_mul(fibonacci_muliplier())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Values produced by a clang++ program compiled against
    /// reference/HiGHS/highs/util/HighsHash.h (HighsHashHelpers::vector_hash).
    /// Element i of the i32 vector is `i*7 - 3 + (i%5)*1000003`, and of the u32
    /// vector is `i*2654435761 + 12345` (mod 2^32).
    const ORACLE: [(usize, u64, u64); 23] = [
        (0, 0, 0),
        (1, 18066512831829289727, 12738888866088011046),
        (2, 5269952368050658456, 10038416437632362649),
        (3, 14957469458717087907, 6750950801160194126),
        (5, 1660359693080669444, 6078316381795079525),
        (7, 5324066511769499736, 5891087768370553016),
        (8, 8942758374592901160, 13285590032589244479),
        (9, 14644758180582401905, 4362999017505538702),
        (15, 4953609939601876239, 14955536416479365000),
        (16, 13296875068429162256, 15243038943933789891),
        (17, 3597993975471483282, 6108392249517840999),
        (31, 13771616110269400715, 8226856800219634744),
        (32, 10434713559388625760, 1277083896547136456),
        (33, 1888785692854633354, 11427438681791274623),
        (63, 8385256035722004451, 7726427265529963742),
        (64, 7427174651085089870, 12552592342376518661),
        (65, 15248458340375526750, 15905829379904038934),
        (100, 10901125039600947850, 9530128370729217634),
        (255, 16535976398401759766, 16347655852190531916),
        (256, 5304613471544382208, 13641441469152380735),
        (257, 16926998732059236838, 6118991524853340975),
        (300, 2356879311706070306, 8441754265890285908),
        (1000, 16487454728864038364, 11586206324721406942),
    ];

    #[test]
    fn vector_hash_matches_highs_oracle() {
        for &(len, expect_i32, expect_u32) in ORACLE.iter() {
            let iv: Vec<HighsInt> = (0..len)
                .map(|i| (i as i32) * 7 - 3 + ((i % 5) as i32) * 1000003)
                .collect();
            let uv: Vec<u32> = (0..len)
                .map(|i| (i as u32).wrapping_mul(2654435761).wrapping_add(12345))
                .collect();
            assert_eq!(vector_hash(&iv), expect_i32, "i32 len {len}");
            assert_eq!(vector_hash(&uv), expect_u32, "u32 len {len}");
        }
    }
}
