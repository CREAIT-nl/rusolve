//! HiGHS source: parallel/HighsSplitDeque.h (lines 49-94, 193-466, 603-612)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The split deque of one worker (Dinan et al. / van Dijk: a private part
//! `[split, head)` the owner pushes and pops without atomics and a shared part
//! `[tail, split)` thieves steal from with one CAS on the packed `ts` word).
//!
//! Files: this one (layout, owner side, steal), [`worker_bunk`] (`WorkerBunk`, the sleeper
//! stack) and [`stolen`] (inject / notify / wait, `runStolenTask`, `leapfrogStolenTask`,
//! `waitForTaskToFinish`, `yield`).
//!
//! **Access protocol.** `owner_data` is touched only by the owner thread (C++ plain
//! members), hence `UnsafeCell` and [`HighsSplitDeque::owner`]. Everything a thief or a
//! notifier touches is atomic (`split_request`, `stealer_data`, `worker_bunk_data`, the
//! task metadata), except `StealerData::injectedTask`, a plain pointer published by the
//! semaphore release: it is an `AtomicPtr` here, accessed relaxed (same ordering as C++).
//!
//! **Randomness.** `randomSteal` draws the victim from `ownerData.randgen`, seeded with the
//! owner id. Which victim is robbed decides only where a task runs, never what it computes.

pub mod stolen;
pub mod worker_bunk;

pub use worker_bunk::WorkerBunk;

use super::highs_binary_semaphore::HighsBinarySemaphore;
use super::highs_cache_align::CacheAligned;
use super::highs_task::{HighsTask, TaskFn};
use crate::HighsInt;
use crate::util::highs_random::HighsRandom;
use std::cell::UnsafeCell;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicPtr, AtomicU64, Ordering};

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::kTaskArraySize (line 54)
pub const kTaskArraySize: u32 = 8192;

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::Status (lines 341-346)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HighsSplitDequeStatus {
    kEmpty,
    kStolen,
    kWork,
    kOverflown,
}

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::OwnerData (lines 59-69)
#[derive(Debug)]
pub struct OwnerData {
    /// C++ `cache_aligned::shared_ptr<WorkerBunk> workerBunk = nullptr`
    pub worker_bunk: Option<Arc<WorkerBunk>>,
    /// C++ `cache_aligned::unique_ptr<HighsSplitDeque>* workers = nullptr`: the executor's
    /// deque array (`numWorkers` boxes, addresses fixed while the executor lives).
    pub workers: *const Box<HighsSplitDeque>,
    pub randgen: HighsRandom,
    pub head: u32,
    pub split_copy: u32,
    pub num_workers: i32,
    pub owner_id: i32,
    pub root_task: *mut HighsTask,
    pub all_stolen_copy: bool,
}

impl Default for OwnerData {
    fn default() -> Self {
        OwnerData {
            worker_bunk: None,
            workers: std::ptr::null(),
            randgen: HighsRandom::default(),
            head: 0,
            split_copy: 0,
            num_workers: 0,
            owner_id: -1,
            root_task: std::ptr::null_mut(),
            all_stolen_copy: true,
        }
    }
}

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::StealerData (lines 71-76)
#[derive(Debug)]
pub struct StealerData {
    /// C++ `HighsBinarySemaphore semaphore{0}`
    pub semaphore: HighsBinarySemaphore,
    /// C++ `HighsTask* injectedTask{nullptr}`
    pub injected_task: AtomicPtr<HighsTask>,
    /// C++ `std::atomic<uint64_t> ts{0}`: `makeTailSplit(tail, split)`
    pub ts: AtomicU64,
    /// C++ `std::atomic<bool> allStolen{true}`
    pub all_stolen: AtomicBool,
}

impl Default for StealerData {
    fn default() -> Self {
        StealerData {
            semaphore: HighsBinarySemaphore::new(false),
            injected_task: AtomicPtr::new(std::ptr::null_mut()),
            ts: AtomicU64::new(0),
            all_stolen: AtomicBool::new(true),
        }
    }
}

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::WorkerBunkData (lines 91-94)
#[derive(Debug)]
pub struct WorkerBunkData {
    /// C++ `std::atomic<HighsSplitDeque*> nextSleeper{nullptr}`
    pub next_sleeper: AtomicPtr<HighsSplitDeque>,
    pub owner_id: i32,
}

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque (lines 49-613). Field layout as C++
/// (lines 201-205): each part on its own 64-byte line, the task array after them (boxed
/// here, so the deque itself stays small; the slots never move).
#[repr(C, align(64))]
pub struct HighsSplitDeque {
    pub(crate) owner_data: CacheAligned<UnsafeCell<OwnerData>>,
    pub(crate) split_request: CacheAligned<AtomicBool>,
    pub(crate) stealer_data: CacheAligned<StealerData>,
    pub(crate) worker_bunk_data: CacheAligned<WorkerBunkData>,
    pub(crate) task_array: CacheAligned<Box<[HighsTask]>>,
}

// SAFETY: the C++ access protocol (module doc): `owner_data` is owner-thread only, all
// cross-thread state is atomic.
unsafe impl Send for HighsSplitDeque {}
unsafe impl Sync for HighsSplitDeque {}

impl std::fmt::Debug for HighsSplitDeque {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HighsSplitDeque")
            .field("owner_id", &self.worker_bunk_data.owner_id)
            .finish()
    }
}

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::makeTailSplit (lines 82-84)
pub const fn make_tail_split(tail: u32, split: u32) -> u64 {
    ((tail as u64) << 32) | split as u64
}

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::tail (line 86)
pub const fn tail(tail_split: u64) -> u32 {
    (tail_split >> 32) as u32
}

/// HiGHS: HighsSplitDeque.h::HighsSplitDeque::split (lines 88-90)
pub const fn split(tail_split: u64) -> u32 {
    tail_split as u32
}

impl HighsSplitDeque {
    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::HighsSplitDeque (lines 269-289)
    ///
    /// `workers` may point at storage the executor fills after this call (the C++ array
    /// is sized first, then filled); nobody reads it before the executor is built.
    pub fn new(
        worker_bunk: &Arc<WorkerBunk>,
        workers: *const Box<HighsSplitDeque>,
        owner_id: i32,
        num_workers: i32,
    ) -> Self {
        let mut owner_data = OwnerData::default();
        owner_data.owner_id = owner_id;
        owner_data.workers = workers;
        owner_data.num_workers = num_workers;
        owner_data.randgen.initialise(owner_id as u32);
        owner_data.worker_bunk = Some(worker_bunk.clone());

        let worker_bunk_data = WorkerBunkData {
            next_sleeper: AtomicPtr::new(std::ptr::null_mut()),
            owner_id,
        };

        let split_request = AtomicBool::new(false);

        HighsSplitDeque {
            owner_data: CacheAligned::new(UnsafeCell::new(owner_data)),
            split_request: CacheAligned::new(split_request),
            stealer_data: CacheAligned::new(StealerData::default()),
            worker_bunk_data: CacheAligned::new(worker_bunk_data),
            task_array: CacheAligned::new(
                (0..kTaskArraySize)
                    .map(|_| HighsTask::default())
                    .collect::<Vec<_>>()
                    .into_boxed_slice(),
            ),
        }
    }

    /// The owner's private state (C++ `ownerData`).
    ///
    /// # Safety
    /// Only the owner thread of this deque calls it, and no two returned references live
    /// at once.
    #[allow(clippy::mut_from_ref)]
    pub(crate) unsafe fn owner(&self) -> &mut OwnerData {
        unsafe { &mut *self.owner_data.get() }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::growShared (lines 207-241)
    fn grow_shared(&self) {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        let have_jobs = owner
            .worker_bunk
            .as_ref()
            .unwrap()
            .have_jobs
            .load(Ordering::Relaxed);

        let mut split_rq = false;
        if have_jobs == owner.num_workers {
            split_rq = self.split_request.load(Ordering::Relaxed);
            if !split_rq {
                return;
            }
        }

        let new_split = crate::util::highs_utils::cpp_min(kTaskArraySize, owner.head);

        debug_assert!(new_split > owner.split_copy);

        // we want to replace the old split point with the new splitPoint
        // but not alter the upper 32 bits of tail.
        // Hence with xor we can xor or the copy of the current split point
        // to set the lower bits to zero and then xor the bits of the new split
        // point to the lower bits that are then zero. First doing the xor of the
        // old and new split point and then doing the xor with the stealerData
        // will not alter the result. Also the upper 32 bits of the xor mask are
        // zero and will therefore not alter the value of tail.
        let xor_mask = (owner.split_copy ^ new_split) as u64;

        debug_assert!((xor_mask >> 32) == 0);

        self.stealer_data.ts.fetch_xor(xor_mask, Ordering::Release);
        owner.split_copy = new_split;

        if split_rq {
            self.split_request.store(false, Ordering::Relaxed);
        } else {
            owner
                .worker_bunk
                .as_ref()
                .unwrap()
                .publish_work(self as *const HighsSplitDeque);
        }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::shrinkShared (lines 243-266)
    fn shrink_shared(&self) -> bool {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        let mut t = tail(self.stealer_data.ts.load(Ordering::Relaxed));
        let s = owner.split_copy;

        if t != s {
            owner.split_copy = (t + s) / 2;
            // C++ `uint64_t{splitCopy} - s` wraps (splitCopy < s): a modular subtraction
            // of the split bits through fetch_add.
            let delta = (owner.split_copy as u64).wrapping_sub(s as u64);
            t = tail(self.stealer_data.ts.fetch_add(delta, Ordering::AcqRel));
            if t != s {
                if t > owner.split_copy {
                    owner.split_copy = (t + s) / 2;
                    self.stealer_data
                        .ts
                        .store(make_tail_split(t, owner.split_copy), Ordering::Relaxed);
                }

                return false;
            }
        }

        self.stealer_data.all_stolen.store(true, Ordering::Relaxed);
        owner.all_stolen_copy = true;
        owner
            .worker_bunk
            .as_ref()
            .unwrap()
            .have_jobs
            .fetch_add(-1, Ordering::Relaxed);
        true
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::checkInterrupt (lines 291-296): unwinds
    /// with [`super::highs_task::Interrupt`] when the root task was cancelled.
    pub fn check_interrupt(&self) {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        // SAFETY: root_task is the owner's task pointer set by setRootTask, alive while the
        // task group that installed it runs (C++ `rootTask->isCancelled()`).
        if !owner.root_task.is_null() && unsafe { (*owner.root_task).is_cancelled() } {
            owner.root_task = std::ptr::null_mut();
            std::panic::resume_unwind(Box::new(crate::parallel::highs_task::Interrupt));
        }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::cancelTask (lines 298-302)
    pub fn cancel_task(&self, task_index: HighsInt) {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        debug_assert!(task_index < owner.head as HighsInt);
        debug_assert!(task_index >= 0);

        self.task_array[task_index as usize].cancel();
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::setRootTask (lines 304-308)
    pub fn set_root_task(&self, new_root: *mut HighsTask) -> *mut HighsTask {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        let prev_root = owner.root_task;
        owner.root_task = new_root;
        prev_root
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::push (lines 310-339). When the array is
    /// full the task runs at once (`head` still advances, so `pop` sees `kOverflown`).
    pub fn push(&self, f: TaskFn) {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        if owner.head >= kTaskArraySize {
            // task queue is full, execute task directly
            if owner.split_copy < kTaskArraySize && !owner.all_stolen_copy {
                self.grow_shared();
            }

            // SAFETY: owner thread only; re-borrowed after grow_shared took its own.
            let owner = unsafe { self.owner() };
            owner.head += 1;
            f();
            return;
        }

        // SAFETY: owner thread only, slot is not yet published.
        unsafe {
            self.task_array[owner.head as usize].set_task_data(f);
        }
        owner.head += 1;

        if owner.all_stolen_copy {
            debug_assert!(owner.head > 0);
            self.stealer_data.ts.store(
                make_tail_split(owner.head - 1, owner.head),
                Ordering::Release,
            );
            self.stealer_data.all_stolen.store(false, Ordering::Relaxed);
            owner.split_copy = owner.head;
            owner.all_stolen_copy = false;

            if self.split_request.load(Ordering::Relaxed) {
                self.split_request.store(false, Ordering::Relaxed);
            }

            let have_jobs = owner
                .worker_bunk
                .as_ref()
                .unwrap()
                .have_jobs
                .fetch_add(1, Ordering::Release);

            if have_jobs < owner.num_workers - 1 {
                owner
                    .worker_bunk
                    .as_ref()
                    .unwrap()
                    .publish_work(self as *const HighsSplitDeque);
            }
        } else {
            self.grow_shared();
        }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::pop (lines 348-378). The C++
    /// `std::pair<Status, HighsTask*>`.
    pub fn pop(&self) -> (HighsSplitDequeStatus, *mut HighsTask) {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        if owner.head == 0 {
            return (HighsSplitDequeStatus::kEmpty, std::ptr::null_mut());
        }

        if owner.head > kTaskArraySize {
            // task queue was full and the overflown tasks have
            // been directly executed
            owner.head -= 1;
            return (HighsSplitDequeStatus::kOverflown, std::ptr::null_mut());
        }

        if owner.all_stolen_copy {
            return (
                HighsSplitDequeStatus::kStolen,
                &self.task_array[(owner.head - 1) as usize] as *const HighsTask as *mut HighsTask,
            );
        }

        if owner.split_copy == owner.head {
            let shrunk = self.shrink_shared();
            // SAFETY: owner thread only; re-borrowed after shrink_shared took its own.
            let owner = unsafe { self.owner() };
            if shrunk {
                return (
                    HighsSplitDequeStatus::kStolen,
                    &self.task_array[(owner.head - 1) as usize] as *const HighsTask
                        as *mut HighsTask,
                );
            }
        }

        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };
        owner.head -= 1;

        if owner.head == 0 {
            if !owner.all_stolen_copy {
                owner.all_stolen_copy = true;
                self.stealer_data.all_stolen.store(true, Ordering::Relaxed);
                owner
                    .worker_bunk
                    .as_ref()
                    .unwrap()
                    .have_jobs
                    .fetch_add(-1, Ordering::Release);
            }
        } else if owner.head != owner.split_copy {
            self.grow_shared();
        }

        // SAFETY: owner thread only; re-borrowed after grow_shared took its own.
        let owner = unsafe { self.owner() };
        (
            HighsSplitDequeStatus::kWork,
            &self.task_array[owner.head as usize] as *const HighsTask as *mut HighsTask,
        )
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::popStolen (lines 380-387)
    pub fn pop_stolen(&self) {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        owner.head -= 1;
        if !owner.all_stolen_copy {
            owner.all_stolen_copy = true;
            self.stealer_data.all_stolen.store(true, Ordering::Relaxed);
            owner
                .worker_bunk
                .as_ref()
                .unwrap()
                .have_jobs
                .fetch_add(-1, Ordering::Release);
        }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::steal (lines 389-412)
    pub fn steal(&self) -> *mut HighsTask {
        if self.stealer_data.all_stolen.load(Ordering::Relaxed) {
            return std::ptr::null_mut();
        }

        let mut ts = self.stealer_data.ts.load(Ordering::Relaxed);
        let mut t = tail(ts);
        let mut s = split(ts);

        if t < s {
            // C++ compare_exchange_weak stores the current value into `ts` on failure.
            match self.stealer_data.ts.compare_exchange_weak(
                ts,
                make_tail_split(t + 1, s),
                Ordering::Acquire,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    return &self.task_array[t as usize] as *const HighsTask as *mut HighsTask;
                }
                Err(current) => ts = current,
            }

            t = tail(ts);
            s = split(ts);
            if t < s {
                return std::ptr::null_mut();
            }
        }

        if t < kTaskArraySize && !self.split_request.load(Ordering::Relaxed) {
            self.split_request.store(true, Ordering::Relaxed);
        }

        std::ptr::null_mut()
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::stealWithRetryLoop (lines 414-435)
    pub fn steal_with_retry_loop(&self) -> *mut HighsTask {
        if self.stealer_data.all_stolen.load(Ordering::Relaxed) {
            return std::ptr::null_mut();
        }

        let mut ts = self.stealer_data.ts.load(Ordering::Relaxed);
        let mut t = tail(ts);
        let mut s = split(ts);

        while t < s {
            // C++ compare_exchange_weak stores the current value into `ts` on failure.
            match self.stealer_data.ts.compare_exchange_weak(
                ts,
                make_tail_split(t + 1, s),
                Ordering::Acquire,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    return &self.task_array[t as usize] as *const HighsTask as *mut HighsTask;
                }
                Err(current) => ts = current,
            }

            t = tail(ts);
            s = split(ts);
        }

        if t < kTaskArraySize && !self.split_request.load(Ordering::Relaxed) {
            self.split_request.store(true, Ordering::Relaxed);
        }

        std::ptr::null_mut()
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::selfStealAndGetTail (lines 437-456)
    pub(crate) fn self_steal_and_get_tail(&self) -> u32 {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        if owner.all_stolen_copy {
            return owner.split_copy;
        }

        // when we steal from ourself we can simply do a fetch_add predictively
        // instead of a cas loop. If the tail we read like this ends up to be
        // above already equal to the splitPoint then we correct it with a simple
        // store. When tail > split instead of tail == split no wrong result can
        // occur as long as we know that the task at taskArray[split] is not
        // actually considered to be stolen and tail is corrected before the owner
        // enters shrinkShared.

        let t = tail(
            self.stealer_data
                .ts
                .fetch_add(make_tail_split(1, 0), Ordering::Relaxed),
        );

        if t == owner.split_copy {
            self.stealer_data
                .ts
                .store(make_tail_split(t, owner.split_copy), Ordering::Relaxed);
        }

        t
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::randomSteal (lines 458-466)
    pub fn random_steal(&self) -> *mut HighsTask {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };

        let mut next = owner.randgen.integer_sup(owner.num_workers - 1);
        next += (next >= owner.owner_id) as HighsInt;

        debug_assert!(next != owner.owner_id);
        debug_assert!(next >= 0);
        debug_assert!(next < owner.num_workers);

        // SAFETY: workers array is stable and valid while the executor lives.
        unsafe { (*owner.workers.offset(next as isize)).steal() }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::getOwnerId (line 603)
    pub fn get_owner_id(&self) -> i32 {
        // SAFETY: ownerId is written only by the constructor; a shared read (C++ const).
        let owner = unsafe { &*self.owner_data.get() };
        owner.owner_id
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::getNumWorkers (line 605)
    pub fn get_num_workers(&self) -> i32 {
        // SAFETY: numWorkers is written only by the constructor; a shared read (C++ const).
        let owner = unsafe { &*self.owner_data.get() };
        owner.num_workers
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::getCurrentHead (line 607)
    pub fn get_current_head(&self) -> i32 {
        // SAFETY: owner thread only (C++ protocol).
        let owner = unsafe { self.owner() };
        owner.head as i32
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::getWorkerById (lines 609-611)
    pub fn get_worker_by_id(&self, id: i32) -> *const HighsSplitDeque {
        // SAFETY: workers is written only by the constructor; a shared read (C++ const).
        let owner = unsafe { &*self.owner_data.get() };
        // SAFETY: workers array is stable and valid while the executor lives.
        unsafe { (*owner.workers.offset(id as isize)).as_ref() }
    }
}

#[cfg(test)]
mod tests {
    // Unit tests. A one-worker deque exercises the owner side without threads
    // (with one worker `push` never publishes work: `haveJobs < numWorkers - 1` is false).
    use super::*;
    use std::sync::atomic::{AtomicI32, Ordering};

    // Boxed deques, as `HighsTaskExecutor::worker_deques` (C++ `unique_ptr`): each deque
    // keeps a pointer into this vector.
    #[allow(clippy::vec_box)]
    fn single_worker() -> (Arc<WorkerBunk>, Vec<Box<HighsSplitDeque>>) {
        let bunk = Arc::new(WorkerBunk::default());
        // the capacity is reserved up front, so `p` (the buffer) stays valid while pushing
        let mut workers: Vec<Box<HighsSplitDeque>> = Vec::with_capacity(1);
        let p = workers.as_ptr();
        workers.push(Box::new(HighsSplitDeque::new(&bunk, p, 0, 1)));
        (bunk, workers)
    }

    #[test]
    fn tail_split_packing_round_trips() {
        let ts = make_tail_split(7, 9);
        assert_eq!((tail(ts), split(ts)), (7, 9));
        assert_eq!(make_tail_split(u32::MAX, 0) >> 32, u32::MAX as u64);
    }

    #[test]
    fn push_pop_is_lifo_and_empty_at_the_bottom() {
        let (_bunk, workers) = single_worker();
        let d = &*workers[0];
        assert_eq!(d.get_owner_id(), 0);
        assert_eq!(d.get_num_workers(), 1);
        let order = Arc::new(std::sync::Mutex::new(Vec::new()));
        for i in 0..3 {
            let o = order.clone();
            d.push(Box::new(move || o.lock().unwrap().push(i)));
        }
        assert_eq!(d.get_current_head(), 3);
        for _ in 0..3 {
            let (s, t) = d.pop();
            assert_eq!(s, HighsSplitDequeStatus::kWork);
            unsafe { (*t).run() };
        }
        assert_eq!(*order.lock().unwrap(), vec![2, 1, 0]);
        assert_eq!(d.get_current_head(), 0);
    }

    #[test]
    fn overflow_runs_inline_and_pops_as_overflown() {
        let (_bunk, workers) = single_worker();
        let d = &*workers[0];
        let c = Arc::new(AtomicI32::new(0));
        for _ in 0..kTaskArraySize + 2 {
            let c = c.clone();
            d.push(Box::new(move || {
                c.fetch_add(1, Ordering::Relaxed);
            }));
        }
        // the two tasks past the array ran inside push
        assert_eq!(c.load(Ordering::Relaxed), 2);
        assert_eq!(d.pop().0, HighsSplitDequeStatus::kOverflown);
        assert_eq!(d.pop().0, HighsSplitDequeStatus::kOverflown);
        assert_eq!(d.pop().0, HighsSplitDequeStatus::kWork);
    }

    #[test]
    fn cancelled_root_task_interrupts() {
        let (_bunk, workers) = single_worker();
        let d = &*workers[0];
        let mut root = HighsTask::default();
        assert!(d.set_root_task(&mut root).is_null());
        d.check_interrupt();
        root.cancel();
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| d.check_interrupt()));
        let payload = r.expect_err("checkInterrupt must unwind");
        assert!(
            payload
                .downcast_ref::<super::super::highs_task::Interrupt>()
                .is_some()
        );
        // C++ line 293: the root task is cleared before the throw
        assert!(d.set_root_task(std::ptr::null_mut()).is_null());
    }
}
