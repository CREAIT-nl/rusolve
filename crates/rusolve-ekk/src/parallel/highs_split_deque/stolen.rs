//! HiGHS source: parallel/HighsSplitDeque.h (lines 468-601: the stolen-task protocol)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Waking a sleeper, running a stolen task (with `Interrupt` unwinding of its subtasks),
//! helping a thief (leapfrogging) and blocking until a stolen task finishes. The spin
//! phases double `numTries` from `kNumTryFac` until `kMicroSecsBeforeSleep` (wall clock:
//! decides only when to block).
//!
//! `try { ... } catch (const HighsTask::Interrupt&)` in `runStolenTask` is
//! `catch_unwind(AssertUnwindSafe(..))` + `downcast_ref::<Interrupt>()`. Any other payload
//! is a bug in a task: C++ would `std::terminate` (the exception leaves a worker thread),
//! so the port calls `std::process::abort()` rather than leave the owner waiting forever.

use super::HighsSplitDeque;
use crate::parallel::highs_scheduler_constants::{kMicroSecsBeforeSleep, kNumTryFac};
use crate::parallel::highs_spin_mutex::HighsSpinMutex;
use crate::parallel::highs_task::{HighsTask, Interrupt};
use std::panic::{self, AssertUnwindSafe};
use std::sync::atomic::Ordering;
use std::time::Instant;

impl HighsSplitDeque {
    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::injectTaskAndNotify (lines 468-474)
    pub fn inject_task_and_notify(&self, t: *mut HighsTask) {
        self.stealer_data.injected_task.store(t, Ordering::Relaxed);

        // TSAN_ANNOTATE_HAPPENS_BEFORE(&stealerData.injectedTask);
        // (line 471: TSan annotation, dropped in NDEBUG build)

        self.stealer_data.semaphore.release();
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::notify (line 476)
    pub fn notify(&self) {
        self.stealer_data.semaphore.release();
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::wait (line 479)
    pub fn wait(&self) {
        self.stealer_data.semaphore.acquire();
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::runStolenTask (lines 481-562)
    ///
    /// # Safety
    /// Called by the owner thread of `self`; `task` was stolen from some deque (or
    /// injected) and is run by this thread only.
    pub unsafe fn run_stolen_task(&self, task: *mut HighsTask) {
        let prev_root_task = (unsafe { self.owner() }).root_task;
        (unsafe { self.owner() }).root_task = task;
        let current_head = (unsafe { self.owner() }).head;

        // try { HighsSplitDeque* owner = task->run(this); if (owner) owner->notify(); }
        // catch (const HighsTask::Interrupt&) { ... }
        let result = panic::catch_unwind(AssertUnwindSafe(|| {
            // SAFETY: task is a valid HighsTask pointer; self is a valid deque.
            // C++ task->run(this) is the stealer overload.
            let owner = unsafe { (&*task).run_as_stealer(self as *const HighsSplitDeque) };
            if !owner.is_null() {
                // SAFETY: owner is a valid deque pointer returned by run_as_stealer.
                unsafe { (&*owner).notify() };
            }
        }));

        if let Err(payload) = result {
            if payload.downcast_ref::<Interrupt>().is_some() {
                // in case the task was interrupted we unwind and cancel all subtasks of
                // the stolen task

                // first cancel all tasks
                let head_now = (unsafe { self.owner() }).head;
                let mut i = current_head;
                while i < head_now {
                    // SAFETY: task_array slots are valid HighsTask values owned by this
                    // deque's owner thread.
                    unsafe {
                        (*(&self.task_array.0[i as usize] as *const HighsTask)).cancel();
                    }
                    i += 1;
                }

                // now remove them from our deque so that we arrive at the original state
                // before the stolen task was executed
                while (unsafe { self.owner() }).head != current_head {
                    let (status, pop_task) = self.pop();
                    debug_assert!(status != super::HighsSplitDequeStatus::kEmpty);
                    // if the task was not stolen it would be up to this worker to execute
                    // it now and we simply skip its execution as it is cancelled
                    if status != super::HighsSplitDequeStatus::kStolen {
                        continue;
                    }

                    // The task was stolen. Check if the stealer is already finished with
                    // its execution in which case we just remove it from the deque.
                    // SAFETY: pop_task is a valid HighsTask pointer from pop().
                    let stealer = unsafe { (&*pop_task).get_stealer_if_unfinished(None) };
                    if stealer.is_null() {
                        self.pop_stolen();
                        continue;
                    }

                    // The task was stolen and the stealer is still executing the task.
                    // We now wait in a spin loop until the task is marked as finished for
                    // some time. We do not proceed with stealing other tasks as when there
                    // is a cancelled task the likelihood of that task being cancelled too
                    // might be high and our priority is to finish unwinding the chain of
                    // cancelled tasks as fast as possible.
                    // When the spinning proceeds for too long we request a notification
                    // from the stealer when it is finished and yield control to the
                    // operating system until then.
                    let mut num_tries: i32 = kNumTryFac;
                    let t_start = Instant::now();

                    // SAFETY: pop_task is a valid HighsTask pointer.
                    let mut is_finished = unsafe { (&*pop_task).is_finished() };

                    while !is_finished {
                        let mut j = 0;
                        while j < num_tries {
                            HighsSpinMutex::yield_processor();
                            // SAFETY: pop_task is a valid HighsTask pointer.
                            is_finished = unsafe { (&*pop_task).is_finished() };
                            if is_finished {
                                break;
                            }
                            j += 1;
                        }

                        if !is_finished {
                            let num_micro_secs = t_start.elapsed().as_micros() as i64;

                            if num_micro_secs < kMicroSecsBeforeSleep {
                                num_tries *= 2;
                            } else {
                                self.wait_for_task_to_finish(pop_task, stealer);
                                break;
                            }
                        }
                    }

                    // the task is finished and we can safely proceed to unwind to the next
                    // task
                    self.pop_stolen();
                }

                // unwinding is finished for all subtasks and we can mark the task as
                // finished and then notify the owner if it is waiting for a signal
                // SAFETY: task is a valid HighsTask pointer; self is a valid deque.
                let owner = unsafe { (&*task).mark_as_finished(self as *const HighsSplitDeque) };
                if !owner.is_null() {
                    // SAFETY: owner is a valid deque pointer returned by mark_as_finished.
                    unsafe { (&*owner).notify() };
                }
            } else {
                // A real panic in a task: C++ would std::terminate (the exception leaves
                // a worker thread). Abort rather than leave the owner waiting forever.
                std::process::abort();
            }
        }

        (unsafe { self.owner() }).root_task = prev_root_task;
        self.check_interrupt();
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::leapfrogStolenTask (lines 567-582). C++
    /// `HighsSplitDeque*& stealer` is `&mut *const HighsSplitDeque`.
    ///
    /// # Safety
    /// Owner thread of `self`; `task` is a slot of `self` that `pop` reported `kStolen`.
    pub unsafe fn leapfrog_stolen_task(
        &self,
        task: *mut HighsTask,
        stealer: &mut *const HighsSplitDeque,
    ) -> bool {
        let mut cancelled = false;
        // SAFETY: task is a valid HighsTask pointer (owner thread, slot reported kStolen).
        *stealer = unsafe { (&*task).get_stealer_if_unfinished(Some(&mut cancelled)) };

        if stealer.is_null() {
            return true;
        }

        if !cancelled {
            loop {
                // SAFETY: *stealer is a valid deque pointer (the stealer of task).
                let t = unsafe { (**stealer).steal_with_retry_loop() };
                if t.is_null() {
                    break;
                }
                // SAFETY: t is a valid task stolen from the stealer deque; this thread
                // now owns it for execution.
                unsafe { self.run_stolen_task(t) };
                // SAFETY: task is a valid HighsTask pointer.
                if unsafe { (&*task).is_finished() } {
                    break;
                }
            }
        }

        // SAFETY: task is a valid HighsTask pointer.
        unsafe { (&*task).is_finished() }
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::waitForTaskToFinish (lines 584-596)
    #[allow(clippy::not_unsafe_ptr_arg_deref)]
    pub fn wait_for_task_to_finish(&self, t: *mut HighsTask, stealer: *const HighsSplitDeque) {
        let lg = self.stealer_data.semaphore.lock_mutex_for_acquire();
        // exchange the value stored in stealer with the pointer to owner
        // so that the stealer will see this pointer instead of nullptr
        // when it stores whether the task is finished. In that case the
        // stealer will additionally acquire the wait mutex and signal the owner
        // thread that the task is finished

        // SAFETY: t is a valid HighsTask pointer; stealer is a valid deque pointer.
        if !unsafe { (&*t).request_notify_when_finished(self as *const HighsSplitDeque, stealer) } {
            return;
        }

        self.stealer_data.semaphore.acquire_locked(lg);
    }

    /// HiGHS: HighsSplitDeque.h::HighsSplitDeque::yield (lines 598-601). `yield` is a Rust
    /// keyword; the method is `yield_`.
    pub fn yield_(&self) {
        let t = self.random_steal();
        if !t.is_null() {
            // SAFETY: yield_ is called by the owner thread; random_steal returned a task
            // from another deque that this thread now owns for execution.
            unsafe { self.run_stolen_task(t) };
        }
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::super::{HighsSplitDequeStatus, WorkerBunk};
    use super::*;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicI32, Ordering};

    #[test]
    fn stolen_task_runs_on_the_thief_and_the_owner_syncs() {
        let bunk = Arc::new(WorkerBunk::default());
        let mut w: Box<Vec<Box<HighsSplitDeque>>> = Box::new(Vec::with_capacity(2));
        let p = w.as_ptr();
        for i in 0..2 {
            w.push(Box::new(HighsSplitDeque::new(&bunk, p, i, 2)));
        }
        let c = Arc::new(AtomicI32::new(0));
        let c2 = c.clone();
        let d0 = &*w[0];
        d0.push(Box::new(move || {
            c2.fetch_add(1, Ordering::Relaxed);
        }));
        // the thief side, on another thread
        let (a0, a1) = (
            &*w[0] as *const HighsSplitDeque as usize,
            &*w[1] as *const HighsSplitDeque as usize,
        );
        let h = std::thread::spawn(move || {
            let (d0, d1) = (a0 as *const HighsSplitDeque, a1 as *const HighsSplitDeque);
            loop {
                let t = unsafe { (*d0).steal_with_retry_loop() };
                if !t.is_null() {
                    unsafe { (*d1).run_stolen_task(t) };
                    break;
                }
                std::thread::yield_now();
            }
        });
        h.join().unwrap();
        let (s, t) = d0.pop();
        assert_eq!(s, HighsSplitDequeStatus::kStolen);
        let mut stealer = std::ptr::null();
        assert!(unsafe { d0.leapfrog_stolen_task(t, &mut stealer) });
        d0.pop_stolen();
        assert_eq!(c.load(Ordering::Relaxed), 1);
        assert_eq!(d0.get_current_head(), 0);
    }
}
