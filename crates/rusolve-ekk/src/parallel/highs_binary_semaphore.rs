//! HiGHS source: parallel/HighsBinarySemaphore.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The one-slot semaphore a sleeping worker blocks on. `count` is 1 (released), 0
//! (acquired) or -1 (a thread waits on the condvar). The spin phase doubles `spinIters`
//! from 10 until `kMicroSecsBeforeSleep` (wall clock: decides only when to block).
//!
//! `std::mutex` / `std::condition_variable` are `std::sync::Mutex<()>` / `Condvar`;
//! `std::unique_lock<std::mutex>` is `MutexGuard<'_, ()>`. Poisoning cannot happen
//! (nothing panics under the lock); recover with `PoisonError::into_inner` anyway.

use super::highs_cache_align::{CacheAligned, CacheAlignedBox};
use super::highs_scheduler_constants::kMicroSecsBeforeSleep;
use super::highs_spin_mutex::HighsSpinMutex;
use std::sync::atomic::{AtomicI32, Ordering};
use std::sync::{Condvar, Mutex, MutexGuard};
use std::time::Instant;

/// HiGHS: HighsBinarySemaphore.h::HighsBinarySemaphore::Data (lines 23-29)
#[derive(Debug)]
pub struct HighsBinarySemaphoreData {
    pub count: AtomicI32,
    /// C++ `alignas(64) std::mutex mutex`
    pub mutex: CacheAligned<Mutex<()>>,
    pub condvar: Condvar,
}

/// HiGHS: HighsBinarySemaphore.h::HighsBinarySemaphore (lines 22-106)
#[derive(Debug)]
pub struct HighsBinarySemaphore {
    /// C++ `highs::cache_aligned::unique_ptr<Data> data_`
    data_: CacheAlignedBox<HighsBinarySemaphoreData>,
}

impl HighsBinarySemaphore {
    /// HiGHS: HighsBinarySemaphore.h::HighsBinarySemaphore::HighsBinarySemaphore (lines 34-35)
    pub fn new(init: bool) -> Self {
        HighsBinarySemaphore {
            data_: Box::new(CacheAligned::new(HighsBinarySemaphoreData {
                count: AtomicI32::new(init as i32),
                mutex: CacheAligned::new(Mutex::new(())),
                condvar: Condvar::new(),
            })),
        }
    }

    /// HiGHS: HighsBinarySemaphore.h::HighsBinarySemaphore::release (lines 37-43)
    pub fn release(&self) {
        let prev = self.data_.count.swap(1, Ordering::Release);
        if prev < 0 {
            let _lg = self.data_.mutex.lock().unwrap_or_else(|e| e.into_inner());
            self.data_.condvar.notify_one();
        }
    }

    /// HiGHS: HighsBinarySemaphore.h::HighsBinarySemaphore::try_acquire (lines 45-49)
    pub fn try_acquire(&self) -> bool {
        let expected: i32 = 1;
        self.data_
            .count
            .compare_exchange_weak(expected, 0, Ordering::Acquire, Ordering::Relaxed)
            .is_ok()
    }

    /// HiGHS: HighsBinarySemaphore.h::HighsBinarySemaphore::acquire() (lines 51-87)
    pub fn acquire(&self) {
        if self.try_acquire() {
            return;
        }

        let t_start = Instant::now();
        let mut spin_iters: i32 = 10;
        loop {
            for _ in 0..spin_iters {
                if self.data_.count.load(Ordering::Acquire) == 1 {
                    if self.try_acquire() {
                        return;
                    }
                }
                HighsSpinMutex::yield_processor();
            }

            let num_micro_secs = t_start.elapsed().as_micros() as i64;

            if num_micro_secs < kMicroSecsBeforeSleep {
                spin_iters *= 2;
            } else {
                break;
            }
        }

        let lg = self.data_.mutex.lock().unwrap_or_else(|e| e.into_inner());
        let prev = self.data_.count.swap(-1, Ordering::Relaxed);
        if prev == 1 {
            self.data_.count.store(0, Ordering::Relaxed);
            return;
        }

        let mut guard = lg;
        loop {
            guard = self
                .data_
                .condvar
                .wait(guard)
                .unwrap_or_else(|e| e.into_inner());
            if self.data_.count.load(Ordering::Acquire) == 1 {
                break;
            }
        }

        self.data_.count.store(0, Ordering::Relaxed);
    }

    /// HiGHS: HighsBinarySemaphore.h::HighsBinarySemaphore::lockMutexForAcquire (lines 89-91)
    pub fn lock_mutex_for_acquire(&self) -> MutexGuard<'_, ()> {
        self.data_.mutex.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// HiGHS: HighsBinarySemaphore.h::HighsBinarySemaphore::acquire(std::unique_lock) (lines 93-105)
    pub fn acquire_locked(&self, lock_guard: MutexGuard<'_, ()>) {
        let prev = self.data_.count.swap(-1, Ordering::Relaxed);
        if prev == 1 {
            self.data_.count.store(0, Ordering::Relaxed);
            return;
        }

        let mut guard = lock_guard;
        loop {
            guard = self
                .data_
                .condvar
                .wait(guard)
                .unwrap_or_else(|e| e.into_inner());
            if self.data_.count.load(Ordering::Relaxed) == 1 {
                break;
            }
        }

        self.data_.count.store(0, Ordering::Relaxed);
    }
}

impl Default for HighsBinarySemaphore {
    /// C++ `HighsBinarySemaphore(bool init = false)`.
    fn default() -> Self {
        HighsBinarySemaphore::new(false)
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::*;
    use std::sync::Arc;

    #[test]
    fn try_acquire_takes_the_released_slot_once() {
        let s = HighsBinarySemaphore::new(true);
        assert!(s.try_acquire());
        assert!(!s.try_acquire());
        s.release();
        s.release();
        assert!(s.try_acquire());
        assert!(!s.try_acquire());
    }

    #[test]
    fn acquire_blocks_until_another_thread_releases() {
        let s = Arc::new(HighsBinarySemaphore::new(false));
        let s2 = s.clone();
        let h = std::thread::spawn(move || {
            // past kMicroSecsBeforeSleep, so the waiter reaches the condvar
            std::thread::sleep(std::time::Duration::from_millis(20));
            s2.release();
        });
        s.acquire();
        h.join().unwrap();
        assert!(!s.try_acquire());
        // the locked variant: release first, then acquire must return at once
        s.release();
        let g = s.lock_mutex_for_acquire();
        s.acquire_locked(g);
        assert!(!s.try_acquire());
    }
}
