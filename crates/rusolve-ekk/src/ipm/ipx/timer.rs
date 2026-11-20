//! HiGHS source: ipm/ipx/timer.h, timer.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! libc++ `high_resolution_clock` is
//! `steady_clock`; `read()` is the steady time in seconds since a fixed
//! process-wide base (only differences of `read()` are ever used).

use super::*;
use std::sync::OnceLock;
use std::time::Instant;

/// HiGHS: ipx/timer.h::Timer
#[derive(Debug, Clone, Copy)]
pub struct Timer {
    t0_: f64,
    pub offset_: f64,
}

impl Timer {
    /// HiGHS: ipx/timer.cc::Timer::Timer (lines 5-8)
    ///
    /// C++ default arguments: `offset = 0`.
    pub fn new(offset: f64) -> Self {
        let mut t = Timer {
            t0_: 0.0,
            offset_: offset,
        };
        t.reset(true);
        t
    }

    /// HiGHS: ipx/timer.h::Timer::Elapsed
    pub fn elapsed(&self) -> f64 {
        Self::toc(self.t0_)
    }

    /// HiGHS: ipx/timer.cc::Timer::Reset (lines 10-14)
    ///
    /// C++ default arguments: `first = false`.
    pub fn reset(&mut self, first: bool) {
        if !first {
            self.offset_ -= self.t0_;
        }
        self.t0_ = Self::read();
        if !first {
            self.offset_ += self.t0_;
        }
    }

    /// HiGHS: ipx/timer.h::Timer::toc
    fn toc(start: f64) -> f64 {
        Self::read() - start
    }

    /// HiGHS: ipx/timer.h::Timer::read
    fn read() -> f64 {
        static BASE: OnceLock<Instant> = OnceLock::new();
        BASE.get_or_init(Instant::now).elapsed().as_secs_f64()
    }
}

impl Default for Timer {
    fn default() -> Self {
        Timer::new(0.0)
    }
}
