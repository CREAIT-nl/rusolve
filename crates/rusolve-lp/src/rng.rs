//! Tiny deterministic PRNG for tests and benches (no external `rand`).
//!
//! XorShift64* (Vassilios Vrumontas's variant of Marsaglia's xorshift family,
//! Marsaglia 2003, "Xorshift RNGs", *Journal of Statistical Software* 8:1).
//! Quality is far beyond what seeded regression tests need; the point is exact
//! reproducibility across machines and platforms.

/// 64-bit xorshift generator (`xorshift64*` shape: shift trio + multiply).
#[derive(Debug, Clone)]
pub struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    /// Seeds the generator. A zero seed is remapped because xorshift has an
    /// absorbing zero state.
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 {
                0x9E37_79B9_7F4A_7C15
            } else {
                seed
            },
        }
    }

    /// Next raw `u64`.
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform float in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        // Use the top 53 bits for an exact multiple of 2^-53.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform float in `[lo, hi)`.
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rng_is_deterministic_and_bounded() {
        let mut a = XorShift64::new(42);
        let mut b = XorShift64::new(42);
        for _ in 0..1000 {
            let x = a.next_f64();
            assert_eq!(x, b.next_f64());
            assert!((0.0..1.0).contains(&x));
        }
        // Zero seed must not get stuck.
        let mut z = XorShift64::new(0);
        assert_ne!(z.next_u64(), 0);
    }
}
