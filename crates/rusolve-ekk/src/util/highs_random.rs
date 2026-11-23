//! HiGHS source: util/HighsRandom.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! All members are inline in the header. Uses [`crate::util::highs_hash`]
//! (`pair_hash`, `log2i`). `state` arithmetic wraps modulo 2^64 (xorshift +
//! hash). The 32-bit `HighsInt` build path is the one to translate (the
//! `#ifdef HIGHSINT64` branch of `integer()` is not taken).

use crate::HighsInt;
use crate::HighsUInt;
use crate::util::highs_hash::*;

/// HiGHS: HighsRandom.h::HighsRandom
///
/// `Default` is the C++ default constructor `HighsRandom(HighsUInt seed = 0)`, i.e.
/// `HighsRandom::new(0)` (`initialise(0)`), so a default-constructed member draws the same
/// sequence as its C++ counterpart. Use `HighsRandom::new(seed)` for an explicit seed.
#[derive(Clone, Debug)]
pub struct HighsRandom {
    state: u64,
}

impl Default for HighsRandom {
    fn default() -> Self {
        HighsRandom::new(0)
    }
}

impl HighsRandom {
    /// HiGHS: HighsRandom.h::HighsRandom::HighsRandom(HighsUInt seed = 0)
    pub fn new(seed: HighsUInt) -> Self {
        let mut r = HighsRandom { state: 0 };
        r.initialise(seed);
        r
    }

    /// HiGHS: HighsRandom.h::HighsRandom::drawUniform(uint32_t sup, int nbits) (lines 27-59).
    /// Tries the outputs `pair_hash(n, lo, hi)` for n = 0..=7 and 9..=31:
    /// HiGHS skips n = 8 (no `HIGHS_RAND_TRY_OUTPUT(8)`); keep the skip.
    fn draw_uniform_u32(&mut self, sup: u32, nbits: i32) -> u32 {
        debug_assert!(u64::from(sup) <= 1u64 << nbits);
        loop {
            self.advance();
            let lo = self.state as u32;
            let hi = (self.state >> 32) as u32;

            for n in 0..32usize {
                if n == 8 {
                    continue;
                }
                let val = (pair_hash(n, lo, hi) >> (64 - nbits)) as u32;
                if val < sup {
                    return val;
                }
            }
        }
    }

    /// HiGHS: HighsRandom.h::HighsRandom::drawUniform(uint64_t sup, int nbits) (lines 61-112).
    /// Tries n = 0..=7 and 9..=15 (n = 8 skipped, as in HiGHS).
    fn draw_uniform_u64(&mut self, sup: u64, nbits: i32) -> u64 {
        if nbits <= 32 {
            return u64::from(self.draw_uniform_u32(sup as u32, nbits));
        }
        debug_assert!(u128::from(sup) <= 1u128 << nbits);
        loop {
            self.advance();
            let lo = self.state as u32;
            let hi = (self.state >> 32) as u32;

            for n in 0..16usize {
                if n == 8 {
                    continue;
                }
                let val = (pair_hash(2 * n, lo, hi) >> (64 - nbits))
                    ^ (pair_hash(2 * n + 1, lo, hi) >> 32);
                if val < sup {
                    return val;
                }
            }
        }
    }

    /// HiGHS: HighsRandom.h::HighsRandom::initialise
    /// C++ default arguments: seed = 0
    pub fn initialise(&mut self, seed: HighsUInt) {
        self.state = u64::from(seed);
        loop {
            self.state = pair_hash(0, self.state as u32, (self.state >> 32) as u32);
            self.state ^= pair_hash(1, (self.state >> 32) as u32, seed) >> 32;
            if self.state != 0 {
                break;
            }
        }
    }

    /// HiGHS: HighsRandom.h::HighsRandom::advance
    pub fn advance(&mut self) {
        self.state ^= self.state >> 12;
        self.state ^= self.state << 25;
        self.state ^= self.state >> 27;
    }

    /// HiGHS: HighsRandom.h::HighsRandom::integer() (32-bit HighsInt branch)
    pub fn integer(&mut self) -> HighsInt {
        self.advance();
        (pair_hash(0, self.state as u32, (self.state >> 32) as u32) >> 33) as HighsInt
    }

    /// HiGHS: HighsRandom.h::HighsRandom::integer(HighsInt sup)
    pub fn integer_sup(&mut self, sup: HighsInt) -> HighsInt {
        if sup <= 1 {
            return 0;
        }
        let nbits = log2i_u32((sup - 1) as HighsUInt) + 1;
        self.draw_uniform_u32(sup as HighsUInt, nbits) as HighsInt
    }

    /// HiGHS: HighsRandom.h::HighsRandom::integer(HighsInt min, HighsInt sup)
    pub fn integer_range(&mut self, min: HighsInt, sup: HighsInt) -> HighsInt {
        min.wrapping_add(self.integer_sup(sup.wrapping_sub(min)))
    }

    /// HiGHS: HighsRandom.h::HighsRandom::fraction
    pub fn fraction(&mut self) -> f64 {
        self.advance();
        // 52 bit output is in interval [0,2^52-1]
        let lo = self.state as u32;
        let hi = (self.state >> 32) as u32;
        let output: u64 = (pair_hash(0, lo, hi) >> (64 - 52)) ^ (pair_hash(1, lo, hi) >> (64 - 26));
        // compute (1+output) / (2^52+1) which is strictly between 0 and 1
        (1 + output) as f64 * 2.2204460492503125e-16
    }

    /// HiGHS: HighsRandom.h::HighsRandom::closedFraction
    pub fn closed_fraction(&mut self) -> f64 {
        self.advance();
        // 53 bit result is in interval [0,2^53-1]
        let lo = self.state as u32;
        let hi = (self.state >> 32) as u32;
        let output: u64 = (pair_hash(0, lo, hi) >> (64 - 53)) ^ (pair_hash(1, lo, hi) >> 32);
        // compute output / (2^53-1) in double precision which is in the closed
        // interval [0,1]
        output as f64 * 1.1102230246251566e-16
    }

    /// HiGHS: HighsRandom.h::HighsRandom::real
    pub fn real(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.closed_fraction()
    }

    /// HiGHS: HighsRandom.h::HighsRandom::bit
    pub fn bit(&mut self) -> bool {
        self.advance();
        (pair_hash(0, self.state as u32, (self.state >> 32) as u32) >> 63) != 0
    }

    /// HiGHS: HighsRandom.h::HighsRandom::shuffle<T>(T* data, HighsInt N)
    pub fn shuffle<T>(&mut self, data: &mut [T], n: HighsInt) {
        let mut i = n;
        while i > 1 {
            let pos = self.integer_sup(i);
            data.swap(pos as usize, (i - 1) as usize);
            i -= 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::HighsRandom;

    /// C++ `HighsRandom r;` is `HighsRandom(0)`, i.e. `initialise(0)`: a default-constructed
    /// member must draw the same sequence as `HighsRandom::new(0)`, not from state 0.
    #[test]
    fn default_is_seed_zero() {
        let mut a = HighsRandom::default();
        let mut b = HighsRandom::new(0);
        for _ in 0..16 {
            assert_eq!(a.integer(), b.integer());
            assert_eq!(a.real(0.1, 0.9).to_bits(), b.real(0.1, 0.9).to_bits());
        }
    }
}
