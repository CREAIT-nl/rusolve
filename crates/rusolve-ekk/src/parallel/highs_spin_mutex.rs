//! HiGHS source: parallel/HighsSpinMutex.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! A test-and-test-and-set spin lock. `HIGHS_HAVE_MM_PAUSE` is CMake's check that
//! `<immintrin.h>` `_mm_pause()` compiles (CMakeLists.txt lines 419-425), true on x86 only:
//! there `yieldProcessor` is `std::hint::spin_loop()` (`pause`); elsewhere (the aarch64
//! oracle build, `/* #undef HIGHS_HAVE_MM_PAUSE */`) it is `std::this_thread::yield()`,
//! i.e. `std::thread::yield_now()`. Either only changes how long a waiter spins.

use std::sync::atomic::{AtomicBool, Ordering};

/// HiGHS: HighsSpinMutex.h::HighsSpinMutex (lines 22-46)
#[derive(Debug, Default)]
pub struct HighsSpinMutex {
    /// C++ `std::atomic<bool> flag{false}`
    flag: AtomicBool,
}

impl HighsSpinMutex {
    /// C++ default constructor (`flag{false}`).
    pub const fn new() -> Self {
        HighsSpinMutex {
            flag: AtomicBool::new(false),
        }
    }

    /// HiGHS: HighsSpinMutex.h::HighsSpinMutex::yieldProcessor (lines 26-33)
    pub fn yield_processor() {
        // C++ #ifdef HIGHS_HAVE_MM_PAUSE: _mm_pause()
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        std::hint::spin_loop();
        // C++ #else: std::this_thread::yield()
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        std::thread::yield_now();
    }

    /// HiGHS: HighsSpinMutex.h::HighsSpinMutex::try_lock (line 35)
    pub fn try_lock(&self) -> bool {
        !self.flag.swap(true, Ordering::Acquire)
    }

    /// HiGHS: HighsSpinMutex.h::HighsSpinMutex::lock (lines 37-43)
    pub fn lock(&self) {
        loop {
            if !self.flag.swap(true, Ordering::Acquire) {
                return;
            }

            while self.flag.load(Ordering::Relaxed) {
                Self::yield_processor();
            }
        }
    }

    /// HiGHS: HighsSpinMutex.h::HighsSpinMutex::unlock (line 45)
    pub fn unlock(&self) {
        self.flag.store(false, Ordering::Release);
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::*;
    use std::sync::Arc;

    #[test]
    fn try_lock_fails_while_held_and_succeeds_after_unlock() {
        let m = HighsSpinMutex::new();
        assert!(m.try_lock());
        assert!(!m.try_lock());
        m.unlock();
        assert!(m.try_lock());
        m.unlock();
        HighsSpinMutex::yield_processor();
    }

    #[test]
    fn lock_gives_mutual_exclusion_across_threads() {
        struct Counter(std::cell::UnsafeCell<u64>);
        unsafe impl Sync for Counter {}
        let m = Arc::new(HighsSpinMutex::new());
        let c = Arc::new(Counter(std::cell::UnsafeCell::new(0)));
        let hs: Vec<_> = (0..4)
            .map(|_| {
                let (m, c) = (m.clone(), c.clone());
                std::thread::spawn(move || {
                    for _ in 0..10_000 {
                        m.lock();
                        // SAFETY: protected by the spin mutex m
                        unsafe { *c.0.get() += 1 };
                        m.unlock();
                    }
                })
            })
            .collect();
        for h in hs {
            h.join().unwrap();
        }
        // SAFETY: all threads have joined, no concurrent access
        assert_eq!(unsafe { *c.0.get() }, 40_000);
    }
}
