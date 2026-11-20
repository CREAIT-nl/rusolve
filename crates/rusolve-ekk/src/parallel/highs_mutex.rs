//! HiGHS source: parallel/HighsMutex.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `HighsMutex` (alias `highs::parallel::mutex`): a lock that spins, then yields to the
//! scheduler (runs other tasks of this worker), then parks the worker on its semaphore
//! as the single "next owner". HiGHS 1.15.1 never instantiates it outside `parallel/`
//!; it is ported for completeness of the `highs::parallel` namespace.
//! The lock only decides *when* a critical section runs, never what it computes.

use std::sync::atomic::AtomicU32;
use std::sync::atomic::Ordering;

/// HiGHS: HighsMutex.h::HighsMutex::kNumSpinTries (line 17)
pub const kNumSpinTries: i32 = 10;

/// HiGHS: HighsMutex.h::HighsMutex (lines 15-122)
///
/// `state`: bit 0 = locked; bits 1.. = `(ownerId + 1)` of the worker waiting to be
/// handed the lock (0 = none).
#[derive(Debug, Default)]
pub struct HighsMutex {
    state: AtomicU32,
}

impl HighsMutex {
    pub const fn new() -> Self {
        Self {
            state: AtomicU32::new(0),
        }
    }

    /// HiGHS: HighsMutex.h::HighsMutex::try_lock (lines 20-25). `compare_exchange_weak`
    /// (may fail spuriously, as in C++).
    pub fn try_lock(&self) -> bool {
        // HiGHS: HighsMutex.h::HighsMutex::try_lock (lines 20-25)
        let uncontended_state: u32 = 0;
        self.state
            .compare_exchange_weak(uncontended_state, 1, Ordering::Acquire, Ordering::Relaxed)
            .is_ok()
    }

    /// HiGHS: HighsMutex.h::HighsMutex::lock (lines 27-109)
    pub fn lock(&self) {
        // HiGHS: HighsMutex.h::HighsMutex::lock (lines 27-109)

        // First try to acquire the lock directly to have a fast path for
        // uncontended access
        if self.try_lock() {
            return;
        }

        // Now spin a few times to check if the lock becomes available
        for _i in 0..kNumSpinTries {
            if self.state.load(Ordering::Relaxed) == 0 {
                if self.try_lock() {
                    return;
                }
            }

            // HiGHS: HighsSpinMutex.h::HighsSpinMutex::yieldProcessor (lines 26-33)
            crate::parallel::highs_spin_mutex::HighsSpinMutex::yield_processor();
        }

        // Lock is still not available, so now start a loop where we yield
        // to the scheduler after each time we observe the lock as unavailable
        // so that we can see if there are other tasks and then check
        // the lock again. We do so for some microseconds which we measure starting
        // from now.
        let this_deque =
            crate::parallel::highs_task_executor::HighsTaskExecutor::get_this_worker_deque();
        let t_start = std::time::Instant::now();

        loop {
            let mut num_tries = kNumSpinTries;

            for _i in 0..num_tries {
                if self.state.load(Ordering::Relaxed) == 0 {
                    if self.try_lock() {
                        return;
                    }
                }
                // SAFETY: this_deque is the calling thread's own worker deque (C++ thisDeque,
                // threadLocalWorkerDeque); it is valid while the executor lives and yield() is
                // called only by its owner thread, as in C++.
                unsafe { (*this_deque).yield_() };
            }

            let num_micro_secs = t_start.elapsed().as_micros() as i64;

            if num_micro_secs
                < crate::parallel::highs_scheduler_constants::kMicroSecsBeforeGlobalSync
            {
                // C++ line 65 doubles numTries, but numTries is re-initialised at the top
                // of the outer loop (line 50), so this store is dead in C++ too; kept as written.
                #[allow(unused_assignments)]
                {
                    num_tries *= 2;
                }
            } else {
                break;
            }
        }

        // The lock is still not available, now we will try to set ourselves as the
        // next worker to acquire the lock and start wait until we are
        // notified by the current worker holding the lock.
        // SAFETY: this_deque is the calling thread's own worker deque, valid while the
        // executor lives; getOwnerId() is a read of a field written only by the constructor.
        let owner_id = unsafe { (*this_deque).get_owner_id() };
        let owner_id = ((owner_id + 1) as u32) << 1;
        let mut s = self.state.load(Ordering::Relaxed);
        loop {
            // as long as we observe that the lock is available we try to lock it, so
            // that we are guaranteed to have a locked state stored within s.
            while s == 0 {
                match self
                    .state
                    .compare_exchange_weak(s, 1, Ordering::Acquire, Ordering::Relaxed)
                {
                    Ok(_) => return,
                    Err(actual) => s = actual,
                }
            }

            if s & 1u32 != 0 {
                match self.state.compare_exchange_weak(
                    s,
                    owner_id | 1,
                    Ordering::Release,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(actual) => s = actual,
                }
            } else {
                // if we observe that the semaphore is unlocked but has its state not
                // equal to zero it means another worker has just been notified to take
                // the lock, but not yet taken it. Before we can set ourselves as the
                // next worker we need to wait for its state to be in locked state.
                // HiGHS: HighsSpinMutex.h::HighsSpinMutex::yieldProcessor (lines 26-33)
                crate::parallel::highs_spin_mutex::HighsSpinMutex::yield_processor();
                s = self.state.load(Ordering::Relaxed);
            }
        }

        // once we are here we have successfully set ourselves as the next worker
        // and once the previous worker released the lock we will be notified. Hence
        // we wait until a notify and then are free to acquire the lock
        // unconditionally.
        // SAFETY: this_deque is the calling thread's own worker deque, valid while the
        // executor lives; wait() is called by its owner thread, as in C++ thisDeque->wait().
        unsafe { (*this_deque).wait() };

        // set the state to have its locked bit set and store the previous next
        // worker as the new next worker. As we only get here when the compare
        // exchange in the loop above succeeds with s having its locked bit set we
        // can simply store back s.
        self.state.store(s, Ordering::Relaxed);
    }

    /// HiGHS: HighsMutex.h::HighsMutex::unlock (lines 111-121). `fetch_add(UINT_MAX)` is
    /// a wrapping decrement.
    pub fn unlock(&self) {
        // HiGHS: HighsMutex.h::HighsMutex::unlock (lines 111-121)
        let prev_state = self.state.fetch_add(u32::MAX, Ordering::Release);

        if prev_state != 1 {
            // C++ unsigned arithmetic: (prevState >> 1) - 1 wraps; the unsigned value is
            // converted to getWorkerById's int parameter.
            let notify_worker_id = (prev_state >> 1).wrapping_sub(1) as i32;
            let this_deque =
                crate::parallel::highs_task_executor::HighsTaskExecutor::get_this_worker_deque();
            // SAFETY: this_deque is the calling thread's own worker deque, valid while the
            // executor lives; getWorkerById() returns a pointer into the executor's stable
            // workers array, as in C++ (line 117-119).
            let worker = unsafe { (*this_deque).get_worker_by_id(notify_worker_id) };
            // SAFETY: worker points into the executor's stable workers array (see above);
            // notify() is the semaphore release that the C++ calls on the notified worker.
            unsafe { (*worker).notify() };
        }
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::*;
    use crate::parallel::highs_parallel::{for_each, initialize_scheduler};
    use crate::parallel::highs_task_executor::HighsTaskExecutor;
    use std::cell::UnsafeCell;

    #[test]
    fn try_lock_then_unlock() {
        let m = HighsMutex::new();
        // weak CAS may fail spuriously: retry
        while !m.try_lock() {}
        assert!(!m.try_lock());
        m.unlock();
        while !m.try_lock() {}
        m.unlock();
    }

    struct Counter(UnsafeCell<i64>);
    // SAFETY: only touched under the mutex
    unsafe impl Sync for Counter {}
    impl Counter {
        /// # Safety
        /// The caller holds the mutex.
        unsafe fn bump(&self) {
            // SAFETY: forwarded
            unsafe { *self.0.get() += 1 };
        }
    }

    #[test]
    fn mutual_exclusion_under_contention() {
        std::thread::spawn(|| {
            initialize_scheduler(4);
            let m = HighsMutex::new();
            let c = Counter(UnsafeCell::new(0));
            for_each(
                0,
                4000,
                &|s, e| {
                    for _ in s..e {
                        m.lock();
                        // SAFETY: under the lock
                        unsafe { c.bump() };
                        m.unlock();
                    }
                },
                1,
            );
            assert_eq!(c.0.into_inner(), 4000);
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }
}
