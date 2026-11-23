//! HiGHS source: parallel/HighsTask.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! One slot of a [`HighsSplitDeque`] task array. `metadata.stealer` packs the stealer's
//! (later: the waiting owner's) deque address with the two flag bits (`HighsSplitDeque` is
//! 64-byte aligned, so the low bits are free).
//!
//! The C++ inline `Callable<F>` storage (`taskData`, at most 56 bytes, `F` trivially
//! destructible) is `UnsafeCell<Option<TaskFn>>`. `Callable::operator()` moves the functor
//! out before calling it, which is `Option::take` + call. Access to the cell follows the
//! C++ protocol: the owner writes it in `setTaskData` before publishing the slot; the one
//! thread that wins `metadata.stealer` (owner `run()`: state `0`; stealer `run(stealer)`:
//! `fetch_or` returned `0`) calls it.

use super::highs_split_deque::HighsSplitDeque;
use std::cell::UnsafeCell;
use std::sync::atomic::AtomicUsize;

/// A spawned task body. `'static` here is erased: [`super::highs_parallel::TaskGroup`]
/// waits for every task it spawned before the borrowed environment can go away.
pub type TaskFn = Box<dyn FnOnce() + Send + 'static>;

/// HiGHS: HighsTask.h::HighsTask::kMaxTaskSize (line 25). Kept for reference: the boxed
/// callable has no size limit.
pub const kMaxTaskSize: usize = 64;

/// HiGHS: HighsTask.h::HighsTask::kFinishedFlag (line 31)
pub const kFinishedFlag: usize = 1;
/// HiGHS: HighsTask.h::HighsTask::kCancelFlag (line 32)
pub const kCancelFlag: usize = 2;
/// HiGHS: HighsTask.h::HighsTask::kPtrMask (line 33)
pub const kPtrMask: usize = !(kFinishedFlag | kCancelFlag);

/// HiGHS: HighsTask.h::HighsTask::Interrupt (line 28). The unwind payload of
/// `HighsSplitDeque::checkInterrupt`: raised with `std::panic::resume_unwind(Box::new(Interrupt))`,
/// caught with `catch_unwind` + `downcast_ref::<Interrupt>()`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Interrupt;

/// HiGHS: HighsTask.h::HighsTask::Metadata (lines 34-36)
#[derive(Debug, Default)]
pub struct HighsTaskMetadata {
    /// C++ `std::atomic<uintptr_t> stealer`
    pub stealer: AtomicUsize,
}

/// HiGHS: HighsTask.h::HighsTask (lines 20-168)
#[derive(Default)]
pub struct HighsTask {
    /// C++ `char taskData[kMaxTaskSize - sizeof(Metadata)]` holding a `Callable<F>`
    task_data: UnsafeCell<Option<TaskFn>>,
    metadata: HighsTaskMetadata,
}

// SAFETY: `task_data` is accessed under the C++ ownership protocol (module doc): written by
// the owner before the slot is published through the deque's release stores, called once by
// the thread that wins `metadata.stealer`.
unsafe impl Send for HighsTask {}
unsafe impl Sync for HighsTask {}

impl std::fmt::Debug for HighsTask {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HighsTask")
            .field("metadata", &self.metadata)
            .finish()
    }
}

impl HighsTask {
    /// HiGHS: HighsTask.h::HighsTask::getCallable (lines 59-67) together with
    /// `Callable::operator()` (lines 50-53): take the functor out and call it.
    ///
    /// # Safety
    /// The caller is the one thread allowed to run this task (module doc).
    unsafe fn call_callable(&self) {
        // SAFETY: The caller is the one thread allowed to run this task (module doc), so no
        // other thread touches `task_data` (C++ `taskData` is read only by the winner of
        // `metadata.stealer`).
        let cell = unsafe { &mut *self.task_data.get() };
        // C++ Callable::operator() moves the functor out before calling it (lines 51-52).
        if let Some(f) = cell.take() {
            f();
        }
    }

    /// HiGHS: HighsTask.h::HighsTask::markAsFinished (lines 69-77)
    pub(crate) fn mark_as_finished(
        &self,
        stealer: *const HighsSplitDeque,
    ) -> *const HighsSplitDeque {
        // C++ line 70-71: exchange with kFinishedFlag, memory_order_release
        let state = self
            .metadata
            .stealer
            .swap(kFinishedFlag, std::sync::atomic::Ordering::Release);
        // C++ line 72-73: reinterpret_cast<HighsSplitDeque*>(state & kPtrMask)
        let waiting_owner = (state & kPtrMask) as *const HighsSplitDeque;
        // C++ line 74: if (waitingOwner != stealer) return waitingOwner;
        if waiting_owner != stealer {
            return waiting_owner;
        }
        // C++ line 76: return nullptr;
        std::ptr::null()
    }

    /// HiGHS: HighsTask.h::HighsTask::run(HighsSplitDeque* stealer) (lines 81-87): run as
    /// stealer; returns the owner's deque if the owner waits and must be signalled.
    ///
    /// # Safety
    /// `self` was obtained by a successful steal (see `HighsSplitDeque::steal`).
    pub(crate) unsafe fn run_as_stealer(
        &self,
        stealer: *const HighsSplitDeque,
    ) -> *const HighsSplitDeque {
        // C++ line 82-83: fetch_or(stealer, memory_order_acquire)
        let state = self
            .metadata
            .stealer
            .fetch_or(stealer as usize, std::sync::atomic::Ordering::Acquire);
        // C++ line 84: if (state == 0) getCallable()();
        if state == 0 {
            // SAFETY: `fetch_or` returned 0, so this thread is the unique winner of
            // `metadata.stealer` and the only one that calls the callable (C++ line 84).
            unsafe { self.call_callable() };
        }
        // C++ line 86: return markAsFinished(stealer);
        self.mark_as_finished(stealer)
    }

    /// HiGHS: HighsTask.h::HighsTask::setTaskData (lines 92-104)
    ///
    /// # Safety
    /// Only the owning deque calls this, on a slot no other thread can reach.
    pub(crate) unsafe fn set_task_data(&self, f: TaskFn) {
        // C++ line 98: metadata.stealer.store(0, std::memory_order_relaxed);
        self.metadata
            .stealer
            .store(0, std::sync::atomic::Ordering::Relaxed);
        // C++ line 100: new (taskData) Callable<F>(std::forward<F>(f));
        // SAFETY: Only the owning deque calls this, on a slot no other thread can reach
        // (C++ `setTaskData` runs before the slot is published by a release store).
        let cell = unsafe { &mut *self.task_data.get() };
        *cell = Some(f);
        // C++ lines 102-103: assert(static_cast<CallableBase*>(reinterpret_cast<Callable<F>*>(
        //               taskData)) == reinterpret_cast<CallableBase*>(taskData));
        // This assert checks that the Callable<F> is at the start of taskData, which is
        // trivially true in our representation since we store the Box directly.
        // (no runtime check: the Box is the storage, so the layout assert has no analogue)
    }

    /// HiGHS: HighsTask.h::HighsTask::cancel (lines 106-108)
    pub fn cancel(&self) {
        // C++ line 107: metadata.stealer.fetch_or(kCancelFlag, std::memory_order_release);
        self.metadata
            .stealer
            .fetch_or(kCancelFlag, std::sync::atomic::Ordering::Release);
    }

    /// HiGHS: HighsTask.h::HighsTask::run() (lines 111-113): run as owner, if not cancelled.
    ///
    /// # Safety
    /// Only the owning deque calls this, after `pop` returned the slot as `kWork`.
    pub(crate) unsafe fn run(&self) {
        // C++ line 112: if (metadata.stealer.load(std::memory_order_relaxed) == 0) getCallable()();
        if self
            .metadata
            .stealer
            .load(std::sync::atomic::Ordering::Relaxed)
            == 0
        {
            // SAFETY: Owner-only run; `metadata.stealer` is 0, so no stealer has won the
            // task and the owner is the unique caller (C++ line 112).
            unsafe { self.call_callable() };
        }
    }

    /// HiGHS: HighsTask.h::HighsTask::requestNotifyWhenFinished (lines 120-130)
    pub(crate) fn request_notify_when_finished(
        &self,
        owner: *const HighsSplitDeque,
        stealer: *const HighsSplitDeque,
    ) -> bool {
        // C++ line 122-123: uintptr_t xormask = reinterpret_cast<uintptr_t>(owner) ^
        //                        reinterpret_cast<uintptr_t>(stealer);
        let xormask = (owner as usize) ^ (stealer as usize);
        // C++ line 124-125: uintptr_t state = metadata.stealer.fetch_xor(xormask, std::memory_order_acquire);
        let state = self
            .metadata
            .stealer
            .fetch_xor(xormask, std::sync::atomic::Ordering::Acquire);
        // C++ line 127: assert(stealer != nullptr);
        debug_assert!(!stealer.is_null());
        // C++ line 129: return (state & kFinishedFlag) == 0;
        (state & kFinishedFlag) == 0
    }

    /// HiGHS: HighsTask.h::HighsTask::isFinished (lines 133-136)
    pub fn is_finished(&self) -> bool {
        // C++ line 134: uintptr_t state = metadata.stealer.load(std::memory_order_acquire);
        let state = self
            .metadata
            .stealer
            .load(std::sync::atomic::Ordering::Acquire);
        // C++ line 135: return state & kFinishedFlag;
        (state & kFinishedFlag) != 0
    }

    /// HiGHS: HighsTask.h::HighsTask::isCancelled (lines 139-142)
    pub fn is_cancelled(&self) -> bool {
        // C++ line 140: uintptr_t state = metadata.stealer.load(std::memory_order_relaxed);
        let state = self
            .metadata
            .stealer
            .load(std::sync::atomic::Ordering::Relaxed);
        // C++ line 141: return state & kCancelFlag;
        (state & kCancelFlag) != 0
    }

    /// HiGHS: HighsTask.h::HighsTask::getStealerIfUnfinished (lines 147-167). C++
    /// `bool* cancelled = nullptr` is `Option<&mut bool>`.
    pub(crate) fn get_stealer_if_unfinished(
        &self,
        cancelled: Option<&mut bool>,
    ) -> *const HighsSplitDeque {
        // C++ line 148: uintptr_t state = metadata.stealer.load(std::memory_order_acquire);
        let mut state = self
            .metadata
            .stealer
            .load(std::sync::atomic::Ordering::Acquire);
        // C++ line 149-150: if (state & kFinishedFlag) return nullptr;
        if state & kFinishedFlag != 0 {
            return std::ptr::null();
        }
        // C++ line 151-160: else { while ((state & ~kCancelFlag) == 0) { ... } }
        // Note: ~kCancelFlag in C++ with uintptr_t is all bits except bit 1.
        // In Rust, !kCancelFlag on usize is all bits except bit 1.
        // The condition (state & ~kCancelFlag) == 0 means state is either 0 or kCancelFlag.
        while (state & !kCancelFlag) == 0 {
            // C++ line 157: HighsSpinMutex::yieldProcessor();
            std::thread::yield_now();
            // C++ line 158: state = metadata.stealer.load(std::memory_order_acquire);
            state = self
                .metadata
                .stealer
                .load(std::sync::atomic::Ordering::Acquire);
        }
        // C++ line 162: if (state & kFinishedFlag) return nullptr;
        if state & kFinishedFlag != 0 {
            return std::ptr::null();
        }
        // C++ line 164: if (cancelled) *cancelled = state & kCancelFlag;
        if let Some(c) = cancelled {
            *c = (state & kCancelFlag) != 0;
        }
        // C++ line 166: return reinterpret_cast<HighsSplitDeque*>(state & kPtrMask);
        (state & kPtrMask) as *const HighsSplitDeque
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::*;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicI32, Ordering};

    fn counting_task(counter: &Arc<AtomicI32>) -> TaskFn {
        let c = counter.clone();
        Box::new(move || {
            c.fetch_add(1, Ordering::Relaxed);
        })
    }

    /// A fake deque address: 64-byte aligned, so the flag bits are free (never dereferenced).
    fn fake_deque(n: usize) -> *const HighsSplitDeque {
        (n * 64) as *const HighsSplitDeque
    }

    #[test]
    fn owner_run_calls_once_and_skips_cancelled() {
        let c = Arc::new(AtomicI32::new(0));
        let t = HighsTask::default();
        unsafe { t.set_task_data(counting_task(&c)) };
        assert!(!t.is_finished() && !t.is_cancelled());
        unsafe { t.run() };
        assert_eq!(c.load(Ordering::Relaxed), 1);
        // owner run() does not mark the task finished (C++ line 112)
        assert!(!t.is_finished());

        unsafe { t.set_task_data(counting_task(&c)) };
        t.cancel();
        assert!(t.is_cancelled());
        unsafe { t.run() };
        assert_eq!(c.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn stealer_run_marks_finished_and_reports_waiting_owner() {
        let c = Arc::new(AtomicI32::new(0));
        let t = HighsTask::default();
        let (owner, stealer) = (fake_deque(3), fake_deque(5));

        // nobody waits: markAsFinished returns null
        unsafe { t.set_task_data(counting_task(&c)) };
        assert!(unsafe { t.run_as_stealer(stealer) }.is_null());
        assert!(t.is_finished());
        assert!(t.get_stealer_if_unfinished(None).is_null());
        assert_eq!(c.load(Ordering::Relaxed), 1);

        // the owner registers before the stealer finishes: the stealer reports the owner
        unsafe { t.set_task_data(counting_task(&c)) };
        t.metadata
            .stealer
            .store(stealer as usize, Ordering::Relaxed);
        let mut cancelled = true;
        assert_eq!(t.get_stealer_if_unfinished(Some(&mut cancelled)), stealer);
        assert!(!cancelled);
        assert!(t.request_notify_when_finished(owner, stealer));
        assert_eq!(t.mark_as_finished(stealer), owner);
        assert!(t.is_finished());
        assert!(!t.request_notify_when_finished(owner, stealer) || !t.is_finished());

        // a cancelled task is not called by the stealer but still finishes
        unsafe { t.set_task_data(counting_task(&c)) };
        t.cancel();
        assert!(unsafe { t.run_as_stealer(stealer) }.is_null());
        assert!(t.is_finished());
        assert_eq!(c.load(Ordering::Relaxed), 1);
    }
}
