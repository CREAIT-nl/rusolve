//! HiGHS source: parallel/HighsCombinable.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `HighsCombinable<T>`: one lazily constructed `T` per worker, indexed by the calling
//! worker's `getOwnerId()`. Its one HiGHS 1.15.1 user is
//! `HighsCliqueTable.cpp::queryNeighbourhood` (line 494), which only calls `combine_each`
//!.
//!
//! **Faithful bug.** C++ `combine` (lines 84-104) has a `break` in its second loop, so it
//! folds at most the first *two* initialised copies (in worker-id order) and drops the
//! rest. The port keeps that behaviour (no HiGHS caller reaches it); `combine_each` visits every initialised copy in
//! worker-id order. Worker-id order is fixed, but *which* workers hold a copy depends on
//! stealing, so a caller is deterministic only if its combine is order- and
//! partition-independent.

use super::highs_cache_align::CacheAligned;
use super::highs_task_executor::HighsTaskExecutor;
use std::cell::UnsafeCell;

/// HiGHS: HighsCombinable.h::HighsCombinable::PaddedData (lines 20-23).
/// `initialized_` + placement-new `data_` is an `Option<T>`.
#[derive(Debug)]
pub struct PaddedData<T> {
    pub(crate) data: UnsafeCell<Option<T>>,
}

/// HiGHS: HighsCombinable.h::HighsCombinable (lines 17-109)
pub struct HighsCombinable<T, FConstruct = fn() -> T>
where
    FConstruct: Fn() -> T,
{
    /// C++ `construct_`
    construct: FConstruct,
    /// C++ `numThreads`
    num_threads: i32,
    /// C++ `threadCopies_` (cache-aligned array)
    thread_copies: Box<[CacheAligned<PaddedData<T>>]>,
}

// SAFETY: slot `i` is only touched by the worker whose owner id is `i` (`local`), or by the
// owner of the combinable after all tasks are synced (`clear`, `combine*`), as in C++.
unsafe impl<T: Send, FConstruct: Fn() -> T + Sync> Sync for HighsCombinable<T, FConstruct> {}

impl<T: Default> HighsCombinable<T, fn() -> T> {
    /// HiGHS: HighsCombinable.h::HighsCombinable() (lines 28-34): `construct_ = T()`,
    /// `numThreads = HighsTaskExecutor::getNumWorkerThreads()`.
    pub fn new() -> Self {
        // C++ line 29: construct_ = []() { return T(); }
        let construct: fn() -> T = || T::default();
        // C++ line 30: numThreads = HighsTaskExecutor::getNumWorkerThreads()
        let num_threads = HighsTaskExecutor::get_num_worker_threads();
        // C++ lines 31-32: threadCopies_ = make_unique_array<PaddedData>(numThreads)
        // C++ line 33: for (int i = 0; i < numThreads; ++i) threadCopies_[i].initialized_ = false;
        let thread_copies: Box<[CacheAligned<PaddedData<T>>]> = (0..num_threads)
            .map(|_| {
                CacheAligned::new(PaddedData {
                    data: UnsafeCell::new(None),
                })
            })
            .collect();
        HighsCombinable {
            construct,
            num_threads,
            thread_copies,
        }
    }
}

impl<T, FConstruct: Fn() -> T> HighsCombinable<T, FConstruct> {
    /// HiGHS: HighsCombinable.h::HighsCombinable(FConstruct&&) (lines 36-43) and
    /// `makeHighsCombinable` (lines 111-114)
    pub fn with_construct(construct: FConstruct) -> Self {
        // C++ line 38: construct_ = std::forward<FConstruct>(fConstruct)
        // C++ line 39: numThreads = HighsTaskExecutor::getNumWorkerThreads()
        let num_threads = HighsTaskExecutor::get_num_worker_threads();
        // C++ lines 40-41: threadCopies_ = make_unique_array<PaddedData>(numThreads)
        // C++ line 42: for (int i = 0; i < numThreads; ++i) threadCopies_[i].initialized_ = false;
        let thread_copies: Box<[CacheAligned<PaddedData<T>>]> = (0..num_threads)
            .map(|_| {
                CacheAligned::new(PaddedData {
                    data: UnsafeCell::new(None),
                })
            })
            .collect();
        HighsCombinable {
            construct,
            num_threads,
            thread_copies,
        }
    }

    /// HiGHS: HighsCombinable.h::HighsCombinable::clear (lines 49-56)
    pub fn clear(&mut self) {
        // C++ lines 50-55: for (int i = 0; i < numThreads; ++i) {
        //   if (threadCopies_[i].initialized_) {
        //     threadCopies_[i].initialized_ = false;
        //     threadCopies_[i].data_.~T();
        //   }
        // }
        for i in 0..self.num_threads {
            // SAFETY: clear runs after all tasks are synced; no worker is concurrently
            // accessing slot i (owner-only protocol).
            let slot = unsafe { &mut *self.thread_copies[i as usize].0.data.get() };
            if slot.is_some() {
                *slot = None;
            }
        }
    }

    /// HiGHS: HighsCombinable.h::HighsCombinable::local (lines 58-76). Both C++ overloads
    /// (const and non-const) lazily construct; one `&self` method with interior
    /// mutability covers them.
    #[allow(clippy::mut_from_ref)]
    pub fn local(&self) -> &mut T {
        // C++ line 59: int threadId = HighsTaskExecutor::getThisWorkerDeque()->getOwnerId();
        // SAFETY: get_this_worker_deque() never returns null (C++ precondition).
        let thread_id = unsafe { (*HighsTaskExecutor::get_this_worker_deque()).get_owner_id() };
        // SAFETY: slot thread_id is touched only by the worker whose owner id is thread_id
        // (C++ owner-only protocol); no other thread accesses this slot concurrently.
        let slot = unsafe { &mut *self.thread_copies[thread_id as usize].0.data.get() };
        // C++ lines 60-63: if (!threadCopies_[threadId].initialized_) {
        //   threadCopies_[threadId].initialized_ = true;
        //   new (&threadCopies_[threadId].data_) T(construct_());
        // }
        if slot.is_none() {
            *slot = Some((self.construct)());
        }
        // C++ line 65: return threadCopies_[threadId].data_;
        slot.as_mut().unwrap()
    }

    /// HiGHS: HighsCombinable.h::HighsCombinable::combine_each (lines 78-82)
    pub fn combine_each(&mut self, combine: impl FnMut(&mut T)) {
        // C++ lines 80-81: for (int i = 0; i < numThreads; ++i)
        //   if (threadCopies_[i].initialized_) combine(threadCopies_[i].data_);
        let mut combine = combine;
        for i in 0..self.num_threads {
            // SAFETY: combine_each runs after all tasks are synced; no worker is
            // concurrently accessing slot i (owner-only protocol).
            let slot = unsafe { &mut *self.thread_copies[i as usize].0.data.get() };
            if let Some(data) = slot {
                combine(data);
            }
        }
    }

    /// HiGHS: HighsCombinable.h::HighsCombinable::combine (lines 84-104), including the
    /// `break` that stops after the second initialised copy. Moved-from copies stay
    /// "initialised" in C++ (moved-from `T`); here they are left as `T::default()`.
    pub fn combine(&mut self, combine: impl FnMut(T, T) -> T) -> T
    where
        T: Default,
    {
        let mut combine = combine;
        // C++ line 86: T combined;
        let mut combined = T::default();
        // C++ line 87: int i;
        let mut i: i32;
        // C++ lines 88-93: for (i = 0; i < numThreads; ++i) {
        //   if (threadCopies_[i].initialized_) {
        //     combined = std::move(threadCopies_[i].data_);
        //     break;
        //   }
        // }
        i = 0;
        while i < self.num_threads {
            // SAFETY: combine runs after all tasks are synced; no worker is concurrently
            // accessing slot i (owner-only protocol).
            let slot = unsafe { &mut *self.thread_copies[i as usize].0.data.get() };
            // C++ `std::move` leaves the copy initialised (moved-from `T`): take its
            // value and leave `T::default()` in the slot, still `Some`.
            if let Some(data) = slot.as_mut() {
                combined = std::mem::take(data);
                break;
            }
            i += 1;
        }
        // C++ lines 95-101: for (++i; i < numThreads; ++i) {
        //   if (threadCopies_[i].initialized_) {
        //     combined = combine(std::move(combined), std::move(threadCopies_[i].data_));
        //     break;
        //   }
        // }
        i += 1;
        while i < self.num_threads {
            // SAFETY: combine runs after all tasks are synced; no worker is concurrently
            // accessing slot i (owner-only protocol).
            let slot = unsafe { &mut *self.thread_copies[i as usize].0.data.get() };
            if let Some(data) = slot.as_mut() {
                combined = combine(combined, std::mem::take(data));
                break;
            }
            i += 1;
        }
        // C++ line 103: return combined;
        combined
    }
}

impl<T, FConstruct: Fn() -> T> std::fmt::Debug for HighsCombinable<T, FConstruct> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HighsCombinable")
            .field("num_threads", &self.num_threads)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::*;
    use crate::parallel::highs_parallel::{for_each, initialize_scheduler};
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    #[test]
    fn single_thread_local_is_one_slot() {
        std::thread::spawn(|| {
            initialize_scheduler(1);
            let mut c = HighsCombinable::<i64>::new();
            for_each(0, 100, &|s, e| *c.local() += (e - s) as i64, 1);
            let mut seen = Vec::new();
            c.combine_each(|v| seen.push(*v));
            assert_eq!(seen, vec![100]);
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }

    #[test]
    fn combine_folds_at_most_two_copies_like_cpp() {
        std::thread::spawn(|| {
            initialize_scheduler(4);
            let mut c = HighsCombinable::with_construct(|| 0i64);
            // slots 0, 1 and 2 initialised (as if three workers had called local())
            for (slot, v) in c.thread_copies.iter_mut().zip([1i64, 2, 4]) {
                *slot.0.data.get_mut() = Some(v);
            }
            // C++ breaks after the second copy: 1 + 2, the 4 is dropped
            assert_eq!(c.combine(|a, b| a + b), 3);
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }
}
