//! HiGHS source: parallel/HighsSplitDeque.h (lines 97-191: `HighsSplitDeque::WorkerBunk`)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The executor-wide sleeper stack and job counter. `sleeperStack` is a Treiber stack of
//! sleeping deques: the low `kAbaTagShift` bits hold `ownerId + 1` of the top sleeper (0 =
//! empty), the high bits an ABA tag bumped on every push / pop; the links are
//! `workerBunkData.nextSleeper`. `publishWork` hands tasks of the caller's shared part
//! straight to sleepers (`injectTaskAndNotify`).

use super::{HighsSplitDeque, WorkerBunkData};
use crate::parallel::highs_cache_align::CacheAligned;
use crate::parallel::highs_task::HighsTask;
use std::sync::atomic::{AtomicI32, AtomicU64, Ordering};

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunk::kAbaTagShift (line 98)
pub const kAbaTagShift: u64 = 20;
/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunk::kIndexMask (line 99)
pub const kIndexMask: u64 = (1u64 << kAbaTagShift) - 1;

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunk (lines 97-191)
#[derive(Debug)]
pub struct WorkerBunk {
    /// C++ `alignas(64) std::atomic<int> haveJobs`
    pub have_jobs: CacheAligned<AtomicI32>,
    /// C++ `alignas(64) std::atomic<uint64_t> sleeperStack`
    pub sleeper_stack: CacheAligned<AtomicU64>,
}

impl Default for WorkerBunk {
    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunk::WorkerBunk (line 103)
    fn default() -> Self {
        WorkerBunk {
            have_jobs: CacheAligned::new(AtomicI32::new(0)),
            sleeper_stack: CacheAligned::new(AtomicU64::new(0)),
        }
    }
}

impl WorkerBunk {
    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunk::pushSleeper (lines 105-124)
    #[allow(clippy::not_unsafe_ptr_arg_deref)]
    pub fn push_sleeper(&self, deque: *const HighsSplitDeque) {
        let mut stack_state = self.sleeper_stack.load(Ordering::Relaxed);
        let mut new_stack_state: u64;
        let mut head: *const HighsSplitDeque;

        loop {
            // SAFETY: reading another deque's owner-only constant fields (workers) is a
            // shared read; the workers array is fixed while the executor lives.
            head = if stack_state & kIndexMask != 0 {
                unsafe {
                    let workers = (*(*deque).owner_data.get()).workers;
                    let ptr = workers.add((stack_state & kIndexMask) as usize - 1);
                    &**ptr as *const HighsSplitDeque
                }
            } else {
                std::ptr::null()
            };

            // SAFETY: writing the nextSleeper link of the deque being pushed; the Treiber
            // stack CAS below serialises concurrent pushes/pops on the stack head.
            unsafe {
                (&*deque)
                    .worker_bunk_data
                    .next_sleeper
                    .store(head as *mut HighsSplitDeque, Ordering::Relaxed);
            }

            new_stack_state = (stack_state >> kAbaTagShift) + 1;
            // SAFETY: reading the owner-only constant ownerId of the deque being pushed.
            let owner_id = unsafe { (&*deque).worker_bunk_data.owner_id };
            new_stack_state = (new_stack_state << kAbaTagShift) | ((owner_id + 1) as u64);

            match self.sleeper_stack.compare_exchange_weak(
                stack_state,
                new_stack_state,
                Ordering::Release,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => stack_state = actual,
            }
        }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunk::popSleeper (lines 126-150)
    #[allow(clippy::not_unsafe_ptr_arg_deref)]
    pub fn pop_sleeper(&self, local_deque: *const HighsSplitDeque) -> *const HighsSplitDeque {
        let mut stack_state = self.sleeper_stack.load(Ordering::Relaxed);
        let mut head: *const HighsSplitDeque;
        let mut new_head: *const HighsSplitDeque;
        let mut new_stack_state: u64;

        loop {
            if (stack_state & kIndexMask) == 0 {
                return std::ptr::null();
            }
            // SAFETY: reading another deque's owner-only constant fields (workers) is a
            // shared read; the workers array is fixed while the executor lives.
            head = unsafe {
                let workers = (*(*local_deque).owner_data.get()).workers;
                let ptr = workers.add((stack_state & kIndexMask) as usize - 1);
                &**ptr as *const HighsSplitDeque
            };
            // SAFETY: reading the nextSleeper link of the head sleeper; the Treiber stack
            // CAS below serialises concurrent pushes/pops on the stack head.
            new_head = unsafe {
                (&*head)
                    .worker_bunk_data
                    .next_sleeper
                    .load(Ordering::Relaxed)
            };
            let new_head_id: i32 = if !new_head.is_null() {
                // SAFETY: reading the owner-only constant ownerId of the next sleeper.
                unsafe { (&*new_head).worker_bunk_data.owner_id + 1 }
            } else {
                0
            };
            new_stack_state = (stack_state >> kAbaTagShift) + 1;
            new_stack_state = (new_stack_state << kAbaTagShift) | (new_head_id as u64);

            match self.sleeper_stack.compare_exchange_weak(
                stack_state,
                new_stack_state,
                Ordering::Acquire,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => stack_state = actual,
            }
        }

        // SAFETY: clearing the nextSleeper link of the popped sleeper; the CAS above
        // established exclusive ownership of this node.
        unsafe {
            (&*head)
                .worker_bunk_data
                .next_sleeper
                .store(std::ptr::null_mut(), Ordering::Relaxed);
        }

        head
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunk::publishWork (lines 152-181).
    /// Called by the owner of `local_deque` only (it reads and writes its `ownerData`).
    #[allow(clippy::not_unsafe_ptr_arg_deref)]
    pub fn publish_work(&self, local_deque: *const HighsSplitDeque) {
        let mut sleeper = self.pop_sleeper(local_deque);
        while !sleeper.is_null() {
            // SAFETY: owner thread only; self_steal_and_get_tail reads/writes the caller's
            // own ownerData.
            let t = unsafe { (&*local_deque).self_steal_and_get_tail() };
            // SAFETY: owner thread only; reading the caller's own ownerData.
            let split_copy = unsafe { (*(*local_deque).owner_data.get()).split_copy };
            if t == split_copy {
                // SAFETY: owner thread only; reading the caller's own ownerData.
                let head = unsafe { (*(*local_deque).owner_data.get()).head };
                if head == split_copy {
                    // SAFETY: owner thread only; writing the caller's own ownerData.
                    unsafe {
                        (*(*local_deque).owner_data.get()).all_stolen_copy = true;
                    }
                    // SAFETY: owner thread only; writing the caller's own stealerData.
                    unsafe {
                        (&*local_deque)
                            .stealer_data
                            .all_stolen
                            .store(true, Ordering::Relaxed);
                    }
                    self.have_jobs.fetch_add(-1, Ordering::Release);
                }
                self.push_sleeper(sleeper);
                return;
            } else {
                // SAFETY: the task slot at index t is owned by local_deque's owner thread
                // (self_steal_and_get_tail just claimed it); injectTaskAndNotify writes the
                // stealer's injectedTask and releases its semaphore.
                let task_ptr = unsafe {
                    &(&*local_deque).task_array.0[t as usize] as *const HighsTask as *mut HighsTask
                };
                // SAFETY: sleeper is a valid deque pointer popped from the Treiber stack.
                unsafe {
                    (&*sleeper).inject_task_and_notify(task_ptr);
                }
            }

            // SAFETY: owner thread only; reading the caller's own ownerData.
            let split_copy = unsafe { (*(*local_deque).owner_data.get()).split_copy };
            if t == split_copy.wrapping_sub(1) {
                // SAFETY: owner thread only; reading the caller's own ownerData.
                let head = unsafe { (*(*local_deque).owner_data.get()).head };
                if head == split_copy {
                    // SAFETY: owner thread only; writing the caller's own ownerData.
                    unsafe {
                        (*(*local_deque).owner_data.get()).all_stolen_copy = true;
                    }
                    // SAFETY: owner thread only; writing the caller's own stealerData.
                    unsafe {
                        (&*local_deque)
                            .stealer_data
                            .all_stolen
                            .store(true, Ordering::Relaxed);
                    }
                    self.have_jobs.fetch_add(-1, Ordering::Release);
                }
                return;
            }

            sleeper = self.pop_sleeper(local_deque);
        }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunk::waitForNewTask (lines 183-190).
    /// Returns null when the executor injected its stop signal.
    #[allow(clippy::not_unsafe_ptr_arg_deref)]
    pub fn wait_for_new_task(&self, local_deque: *const HighsSplitDeque) -> *mut HighsTask {
        self.push_sleeper(local_deque);
        // SAFETY: local_deque is a valid deque pointer; the semaphore is part of its
        // stealerData which is accessed cross-thread via the semaphore protocol.
        unsafe { (&*local_deque).stealer_data.semaphore.acquire() };

        // TSAN_ANNOTATE_HAPPENS_AFTER(&localDeque->stealerData.injectedTask);
        // (line 187: TSan annotation, dropped in NDEBUG build)

        // SAFETY: the semaphore acquire above synchronises with the release in
        // injectTaskAndNotify, so injectedTask is visible here.
        unsafe {
            (&*local_deque)
                .stealer_data
                .injected_task
                .load(Ordering::Relaxed)
        }
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::*;
    use std::sync::Arc;

    /// `n` deques sharing one bunk, as `HighsTaskExecutor::HighsTaskExecutor` builds them
    /// (no threads are started).
    // Boxed deques, as `HighsTaskExecutor::worker_deques` (C++ `unique_ptr`): each deque
    // keeps a pointer into this vector.
    #[allow(clippy::vec_box)]
    fn deques(n: i32) -> (Arc<WorkerBunk>, Vec<Box<HighsSplitDeque>>) {
        let bunk = Arc::new(WorkerBunk::default());
        // the capacity is reserved up front, so `p` (the buffer) stays valid while pushing
        let mut workers: Vec<Box<HighsSplitDeque>> = Vec::with_capacity(n as usize);
        let p = workers.as_ptr();
        for i in 0..n {
            workers.push(Box::new(HighsSplitDeque::new(&bunk, p, i, n)));
        }
        (bunk, workers)
    }

    #[test]
    fn sleeper_stack_is_lifo_with_aba_tag() {
        let (bunk, w) = deques(3);
        let p = |i: usize| &*w[i] as *const HighsSplitDeque;
        assert!(bunk.pop_sleeper(p(0)).is_null());
        bunk.push_sleeper(p(1));
        bunk.push_sleeper(p(2));
        let s = bunk
            .sleeper_stack
            .load(std::sync::atomic::Ordering::Relaxed);
        assert_eq!(s & kIndexMask, 3);
        assert_eq!(s >> kAbaTagShift, 2);
        assert_eq!(bunk.pop_sleeper(p(0)), p(2));
        assert_eq!(bunk.pop_sleeper(p(0)), p(1));
        assert!(bunk.pop_sleeper(p(0)).is_null());
        let s = bunk
            .sleeper_stack
            .load(std::sync::atomic::Ordering::Relaxed);
        assert_eq!((s & kIndexMask, s >> kAbaTagShift), (0, 4));
    }

    #[test]
    fn wait_for_new_task_returns_the_injected_task() {
        let (bunk, w) = deques(2);
        let d1 = &*w[1] as *const HighsSplitDeque as usize;
        let b2 = bunk.clone();
        let h =
            std::thread::spawn(move || b2.wait_for_new_task(d1 as *const HighsSplitDeque) as usize);
        // wait until the sleeper is on the stack, then hand it the stop signal (null)
        loop {
            let s = bunk.pop_sleeper(&*w[0]);
            if !s.is_null() {
                unsafe { (&*s).inject_task_and_notify(std::ptr::null_mut()) };
                break;
            }
            std::thread::yield_now();
        }
        assert_eq!(h.join().unwrap(), 0);
    }
}
