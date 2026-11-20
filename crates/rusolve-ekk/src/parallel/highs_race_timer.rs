//! HiGHS source: parallel/HighsRaceTimer.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! An atomic minimum over a time limit. HiGHS 1.15.1 instantiates it nowhere outside
//! `parallel/` (it is ported for completeness of the scheduler headers). The C++ template
//! is instantiated here for `f64` (the only type a time limit takes), stored as its bit
//! pattern in an `AtomicU64`. The loop condition `current > newLimit` compares the decoded
//! values; the CAS compares bit patterns, as C++ `std::atomic<double>::compare_exchange_weak`
//! compares object representations.

use std::sync::atomic::{AtomicU64, Ordering};

/// HiGHS: HighsRaceTimer.h::HighsRaceTimer<double> (lines 15-36)
#[derive(Debug)]
pub struct HighsRaceTimer {
    /// C++ `std::atomic<T> limit` (bits of an `f64`)
    limit: AtomicU64,
}

impl Default for HighsRaceTimer {
    /// HiGHS: HighsRaceTimer.h::HighsRaceTimer::HighsRaceTimer (line 20):
    /// `limit(std::numeric_limits<T>::max())`.
    fn default() -> Self {
        HighsRaceTimer {
            limit: AtomicU64::new(f64::MAX.to_bits()),
        }
    }
}

impl HighsRaceTimer {
    /// HiGHS: HighsRaceTimer.h::HighsRaceTimer::decreaseLimit (lines 22-31)
    pub fn decrease_limit(&self, new_limit: f64) {
        let mut current_bits = self.limit.load(Ordering::Relaxed);
        let mut current = f64::from_bits(current_bits);

        while current > new_limit {
            match self.limit.compare_exchange_weak(
                current_bits,
                new_limit.to_bits(),
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual_bits) => {
                    current_bits = actual_bits;
                    current = f64::from_bits(actual_bits);
                }
            }
        }
    }

    /// HiGHS: HighsRaceTimer.h::HighsRaceTimer::limitReached (lines 33-35)
    pub fn limit_reached(&self, current_time: f64) -> bool {
        current_time > f64::from_bits(self.limit.load(Ordering::Relaxed))
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::*;

    #[test]
    fn decrease_limit_keeps_the_minimum() {
        let t = HighsRaceTimer::default();
        assert!(!t.limit_reached(1e300));
        t.decrease_limit(5.0);
        t.decrease_limit(7.0);
        assert!(!t.limit_reached(5.0));
        assert!(t.limit_reached(5.5));
        t.decrease_limit(2.0);
        assert!(t.limit_reached(2.5));
    }
}
