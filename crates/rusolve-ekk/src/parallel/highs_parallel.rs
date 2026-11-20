//! HiGHS source: parallel/HighsParallel.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The `highs::parallel` namespace: `initialize_scheduler`, `num_threads`, `thread_num`,
//! `spawn` / `sync`, `TaskGroup` and `for_each`.
//!
//! **Determinism.** Nothing here makes a result depend on the thread count by itself:
//! `for_each` splits `[start, end)` at `(start + end) >> 1` until a range is at most
//! `grainSize` long, whatever the number of workers, and `sync` / `taskWait` return only
//! when the tasks they wait for are finished. A caller is deterministic if its tasks write
//! disjoint data and it combines their results in a fixed order.
//!
//! **Scoped groups.** C++ `TaskGroup tg; ... ~TaskGroup()` is
//! `task_group(|tg| { ... })`: [`task_group`] runs the destructor sequence (`cancel`,
//! `setRootTask(nullptr)`, `taskWait`, `setRootTask(prev)`, `checkInterrupt`) after the
//! closure returns *or unwinds*, so a spawned task may borrow anything that outlives the
//! call (the `'scope` lifetime, as in `std::thread::scope`). `TaskGroup` is neither `Send`
//! nor `Sync`: like the C++ object it belongs to the worker deque of the thread that made
//! it.
//!
//! **Free `spawn` / `sync`** (HEkkDual's `iterateTasks` / `chooseColumnSlice` use them,
//! not a `TaskGroup`) are `unsafe`: the caller promises to `sync` every task it spawned
//! before anything the task borrows goes out of scope, exactly the C++ contract.

use super::highs_split_deque::HighsSplitDeque;
use super::highs_split_deque::HighsSplitDequeStatus;
use super::highs_task::TaskFn;
use super::highs_task_executor::HighsTaskExecutor;
use crate::HighsInt;
use std::marker::PhantomData;

/// HiGHS: HighsParallel.h::highs::parallel::initialize_scheduler (lines 20-29). The
/// HiGHS C++ default `numThreads = 0`.
///
/// `0` maps to 1 (the `HIGHS_NO_DEFAULT_THREADS` branch, which the oracle build
/// selects), not `(hardware_concurrency + 1) / 2`.
pub fn initialize_scheduler(num_threads: i32) {
    // C++ lines 21-27: if (numThreads == 0) { numThreads = 1; }
    // (HIGHS_NO_DEFAULT_THREADS is defined in the oracle build)
    let mut num_threads = num_threads;
    if num_threads == 0 {
        num_threads = 1;
    }
    // C++ line 28: HighsTaskExecutor::initialize(numThreads);
    HighsTaskExecutor::initialize(num_threads);
}

/// HiGHS: HighsParallel.h::highs::parallel::num_threads (lines 31-33)
pub fn num_threads() -> i32 {
    // C++ line 32: return HighsTaskExecutor::getThisWorkerDeque()->getNumWorkers();
    let deque = HighsTaskExecutor::get_this_worker_deque();
    // SAFETY: get_this_worker_deque returns this thread's worker deque (thread-local,
    // initialized on first use as in C++).
    unsafe { (*deque).get_num_workers() }
}

/// HiGHS: HighsParallel.h::highs::parallel::thread_num (lines 35-37)
pub fn thread_num() -> i32 {
    // C++ line 36: return HighsTaskExecutor::getThisWorkerDeque()->getOwnerId();
    let deque = HighsTaskExecutor::get_this_worker_deque();
    // SAFETY: get_this_worker_deque returns this thread's worker deque (thread-local,
    // initialized on first use as in C++).
    unsafe { (*deque).get_owner_id() }
}

/// HiGHS: HighsParallel.h::highs::parallel::spawn(HighsSplitDeque*, F&&) (lines 39-42)
///
/// # Safety
/// `local_deque` is this thread's worker deque, and the caller syncs the task before
/// anything `f` borrows is dropped.
pub unsafe fn spawn_on<'a, F>(local_deque: *const HighsSplitDeque, f: F)
where
    F: FnOnce() + Send + 'a,
{
    // C++ line 41: localDeque->push(std::forward<F>(f));
    // SAFETY: local_deque is this thread's worker deque (caller contract).
    // The lifetime of f is erased to 'static; the caller guarantees the sync before
    // anything f borrows is dropped (rule (g) / (t)).
    // SAFETY: the caller guarantees the task is run or dropped before 'a ends (see
    // erase_task_lifetime); local_deque is this thread's worker deque (caller contract).
    unsafe {
        let task_fn: TaskFn = erase_task_lifetime(f);
        (*local_deque).push(task_fn);
    }
}

/// HiGHS: HighsParallel.h::highs::parallel::spawn(F&&) (lines 44-47)
///
/// # Safety
/// As [`spawn_on`].
pub unsafe fn spawn<'a, F>(f: F)
where
    F: FnOnce() + Send + 'a,
{
    // C++ line 46: spawn(HighsTaskExecutor::getThisWorkerDeque(), std::forward<F>(f));
    let deque = HighsTaskExecutor::get_this_worker_deque();
    // SAFETY: deque is this thread's worker deque (thread-local, as in C++).
    unsafe { spawn_on(deque, f) };
}

/// HiGHS: HighsParallel.h::highs::parallel::sync(HighsSplitDeque*) (lines 49-64)
///
/// # Safety
/// `local_deque` is this thread's worker deque and holds a task spawned by the caller.
pub unsafe fn sync_on(local_deque: *const HighsSplitDeque) {
    // C++ line 50: std::pair<HighsSplitDeque::Status, HighsTask*> popResult = localDeque->pop();
    // SAFETY: local_deque is this thread's worker deque (caller contract).
    let (status, task) = unsafe { (*local_deque).pop() };

    // C++ lines 51-64: switch (popResult.first)
    match status {
        HighsSplitDequeStatus::kEmpty => {
            // C++ line 53: assert(false);
            debug_assert!(false);
            // C++ line 54: // fall through
            // C++ line 55-58: kOverflown case: break;
            // fall through to kOverflown which breaks
        }
        HighsSplitDequeStatus::kOverflown => {
            // C++ lines 55-58: when the local deque is overflown the task has been
            // executed during spawn already
            // break;
        }
        HighsSplitDequeStatus::kStolen => {
            // C++ line 60: HighsTaskExecutor::sync_stolen_task(localDeque, popResult.second);
            // SAFETY: local_deque is this thread's worker deque; task is a valid stolen task.
            unsafe { HighsTaskExecutor::sync_stolen_task(local_deque, task) };
        }
        HighsSplitDequeStatus::kWork => {
            // C++ line 63: popResult.second->run();
            // SAFETY: task is a valid HighsTask pointer returned by pop().
            unsafe { (*task).run() };
        }
    }
}

/// HiGHS: HighsParallel.h::highs::parallel::sync() (line 66)
///
/// # Safety
/// As [`sync_on`].
pub unsafe fn sync() {
    // C++ line 66: sync(HighsTaskExecutor::getThisWorkerDeque());
    let deque = HighsTaskExecutor::get_this_worker_deque();
    // SAFETY: deque is this thread's worker deque (thread-local, as in C++).
    unsafe { sync_on(deque) };
}

/// Erase the lifetime of a task body (storage only; the caller guarantees the sync).
///
/// # Safety
/// The returned box must be run or dropped before `'a` ends.
pub(crate) unsafe fn erase_task_lifetime<'a, F>(f: F) -> TaskFn
where
    F: FnOnce() + Send + 'a,
{
    // C++ std::forward<F>(f) into the task storage. The lifetime is erased to 'static;
    // the scoped task group (rule (g)) or the caller's sync contract (rule (t)) makes
    // this sound: the task is run or dropped before the borrowed data goes out of scope.
    // Work counters (crate::work): the task runs on whichever worker takes it,
    // so it carries its spawner's solve counters and installs them for its duration.
    let counters = crate::work::current();
    let b: Box<dyn FnOnce() + Send + 'a> = Box::new(move || {
        let _scope = crate::work::TaskScope::enter(counters);
        f()
    });
    // SAFETY: Both Box<dyn FnOnce() + Send + 'a> and TaskFn (Box<dyn FnOnce() + Send +
    // 'static>) are fat pointers with the same size and layout. The only difference is
    // the erased lifetime. The caller guarantees the task is run or dropped before 'a
    // ends (rule (g) scoped group or rule (t) caller sync contract).
    unsafe { std::mem::transmute::<Box<dyn FnOnce() + Send + 'a>, TaskFn>(b) }
}

/// HiGHS: HighsParallel.h::highs::parallel::TaskGroup (lines 67-113)
///
/// `'scope` is the life of the group, `'env` what spawned tasks may borrow (the
/// `std::thread::Scope<'scope, 'env>` pattern).
#[derive(Debug)]
pub struct TaskGroup<'scope, 'env: 'scope> {
    /// C++ `HighsSplitDeque* workerDeque`
    worker_deque: *const HighsSplitDeque,
    /// C++ `int dequeHead`
    deque_head: i32,
    /// Invariant over both lifetimes; the raw pointer makes the group `!Send + !Sync`.
    _scope: PhantomData<&'scope mut &'scope ()>,
    _env: PhantomData<&'env mut &'env ()>,
}

impl<'scope, 'env> TaskGroup<'scope, 'env> {
    /// HiGHS: HighsParallel.h::TaskGroup::TaskGroup (lines 73-76). Private: groups are made
    /// by [`task_group`], which also runs the destructor.
    fn new() -> Self {
        // C++ line 74: workerDeque = HighsTaskExecutor::getThisWorkerDeque();
        let worker_deque = HighsTaskExecutor::get_this_worker_deque();
        // C++ line 75: dequeHead = workerDeque->getCurrentHead();
        // SAFETY: get_this_worker_deque returns this thread's worker deque (thread-local).
        let deque_head = unsafe { (*worker_deque).get_current_head() };

        TaskGroup {
            worker_deque,
            deque_head,
            _scope: PhantomData,
            _env: PhantomData,
        }
    }

    /// HiGHS: HighsParallel.h::TaskGroup::spawn (lines 78-81)
    pub fn spawn<F>(&self, f: F)
    where
        F: FnOnce() + Send + 'scope,
    {
        // C++ line 80: highs::parallel::spawn(workerDeque, std::forward<F>(f));
        // SAFETY: worker_deque is this thread's worker deque (set in new()). The scoped
        // task group guarantees the task is run or cancelled before 'scope ends (rule (g)).
        unsafe { spawn_on(self.worker_deque, f) };
    }

    /// HiGHS: HighsParallel.h::TaskGroup::sync (lines 83-86): wait for the most recently
    /// spawned unfinished task of this group (LIFO).
    pub fn sync(&self) {
        // C++ line 84: assert(workerDeque->getCurrentHead() > dequeHead);
        // SAFETY: worker_deque is a valid pointer to this thread's worker deque.
        debug_assert!(unsafe { (*self.worker_deque).get_current_head() > self.deque_head });
        // C++ line 85: highs::parallel::sync(workerDeque);
        // SAFETY: worker_deque is this thread's worker deque and holds a task spawned
        // by this group (caller contract of TaskGroup::sync).
        unsafe { sync_on(self.worker_deque) };
    }

    /// HiGHS: HighsParallel.h::TaskGroup::taskWait (lines 88-91)
    pub fn task_wait(&self) {
        // C++ lines 89-90: while (workerDeque->getCurrentHead() > dequeHead)
        //                     highs::parallel::sync(workerDeque);
        loop {
            // SAFETY: worker_deque is a valid pointer to this thread's worker deque.
            let head = unsafe { (*self.worker_deque).get_current_head() };
            if head <= self.deque_head {
                break;
            }
            // SAFETY: worker_deque is this thread's worker deque and holds tasks spawned
            // by this group (caller contract of TaskGroup::task_wait).
            unsafe { sync_on(self.worker_deque) };
        }
    }

    /// HiGHS: HighsParallel.h::TaskGroup::cancel (lines 93-96)
    pub fn cancel(&self) {
        // C++ lines 94-95: for (int i = dequeHead; i < workerDeque->getCurrentHead(); ++i)
        //                     workerDeque->cancelTask(i);
        let mut i = self.deque_head;
        // The bound is re-read every iteration, as the C++ loop condition does.
        // SAFETY: worker_deque is a valid pointer to this thread's worker deque.
        while i < unsafe { (*self.worker_deque).get_current_head() } {
            // SAFETY: worker_deque is a valid pointer to this thread's worker deque.
            unsafe { (*self.worker_deque).cancel_task(i) };
            i += 1;
        }
    }

    /// HiGHS: HighsParallel.h::TaskGroup::~TaskGroup (lines 98-112). Called by
    /// [`task_group`] on return and on unwind; may itself unwind with `Interrupt`
    /// (C++ `noexcept(false)`), but not while already unwinding.
    fn finish(&self) {
        // C++ line 99: cancel();
        self.cancel();

        // C++ line 104: HighsTask* prevRoot = workerDeque->setRootTask(nullptr);
        // SAFETY: worker_deque is a valid pointer to this thread's worker deque.
        let prev_root = unsafe { (*self.worker_deque).set_root_task(std::ptr::null_mut()) };

        // C++ line 106: taskWait();
        self.task_wait();

        // C++ line 110: workerDeque->setRootTask(prevRoot);
        // SAFETY: worker_deque is a valid pointer to this thread's worker deque.
        unsafe { (*self.worker_deque).set_root_task(prev_root) };

        // C++ line 111: workerDeque->checkInterrupt();
        // SAFETY: worker_deque is a valid pointer to this thread's worker deque.
        unsafe { (*self.worker_deque).check_interrupt() };
    }
}

/// C++ `{ TaskGroup tg; body(tg); }` (HighsParallel.h 67-113): make a group on this
/// thread's worker deque, run `body`, then run the group destructor (`finish`). The
/// destructor also runs when `body` unwinds (C++ stack unwinding runs `~TaskGroup`); the
/// panic is then resumed, as C++ rethrows after the destructor.
pub fn task_group<'env, R, F>(body: F) -> R
where
    F: for<'scope> FnOnce(&'scope TaskGroup<'scope, 'env>) -> R,
{
    // C++ lines 73-76: TaskGroup tg; (constructor)
    let tg = TaskGroup::new();

    // C++ body(tg); a panic (Interrupt or other) is caught so that ~TaskGroup runs first.
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| body(&tg)));

    // C++ ~TaskGroup(): cancel, setRootTask(nullptr), taskWait, setRootTask(prev),
    // checkInterrupt. Runs on both normal return and unwind.
    match result {
        Ok(r) => {
            // Normal scope exit: the noexcept(false) destructor may throw Interrupt.
            tg.finish();
            r
        }
        Err(p) => {
            // Stack unwinding: a destructor that throws while an exception is in flight
            // is std::terminate in C++.
            if std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| tg.finish())).is_err() {
                std::process::abort();
            }
            // Resume the body's panic after the destructor has run.
            std::panic::resume_unwind(p)
        }
    }
}

/// HiGHS: HighsParallel.h::highs::parallel::for_each (lines 115-133). C++ `grainSize = 1`
/// is passed explicitly.
pub fn for_each<F>(start: HighsInt, end: HighsInt, f: &F, grain_size: HighsInt)
where
    F: Fn(HighsInt, HighsInt) + Sync,
{
    // C++ line 117: if (end - start <= grainSize) { f(start, end); }
    if end - start <= grain_size {
        f(start, end);
    } else {
        // C++ line 120: TaskGroup tg;
        task_group(|tg| {
            // C++ lines 122-128: do { ... } while (end - start > grainSize);
            let mut start = start;
            let mut end = end;
            loop {
                // C++ line 123: HighsInt split = (start + end) >> 1;
                let split = (start + end) >> 1;

                // C++ lines 124-126: tg.spawn([split, end, grainSize, &f]() {
                //                       for_each(split, end, f, grainSize);
                //                     });
                // SAFETY: The scoped task group (rule (g)) guarantees the spawned task
                // is run or cancelled before 'scope ends. f is Sync and outlives the group.
                tg.spawn(move || {
                    for_each(split, end, f, grain_size);
                });

                // C++ line 127: end = split;
                end = split;

                // C++ line 128: } while (end - start > grainSize);
                if end - start <= grain_size {
                    break;
                }
            }

            // C++ line 130: f(start, end);
            f(start, end);

            // C++ line 131: tg.taskWait();
            tg.task_wait();
        });
    }
}

#[cfg(test)]
mod tests {
    // Unit tests. Every test runs on its own thread, so it owns its executor
    // (thread-local, as C++).
    use super::*;
    use crate::parallel::highs_task::Interrupt;
    use std::sync::Mutex;
    use std::sync::atomic::{AtomicI32, Ordering};

    fn on_fresh_thread(threads: i32, f: impl FnOnce() + Send + 'static) {
        std::thread::spawn(move || {
            initialize_scheduler(threads);
            f();
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }

    #[test]
    fn num_threads_and_thread_num_on_the_main_worker() {
        for n in [1, 2, 4, 8] {
            on_fresh_thread(n, move || {
                assert_eq!(num_threads(), n);
                assert_eq!(thread_num(), 0);
            });
        }
        // threads = 0 is the HIGHS_NO_DEFAULT_THREADS default: one worker
        on_fresh_thread(0, || assert_eq!(num_threads(), 1));
    }

    #[test]
    fn single_thread_sync_is_lifo_like_cpp() {
        on_fresh_thread(1, || {
            let order = Mutex::new(Vec::new());
            task_group(|tg| {
                tg.spawn(|| order.lock().unwrap().push('a'));
                tg.spawn(|| order.lock().unwrap().push('b'));
                tg.sync();
                order.lock().unwrap().push('m');
                tg.sync();
            });
            assert_eq!(*order.lock().unwrap(), vec!['b', 'm', 'a']);
        });
    }

    #[test]
    fn group_destructor_cancels_unsynced_tasks_at_one_thread() {
        on_fresh_thread(1, || {
            let ran = AtomicI32::new(0);
            task_group(|tg| {
                tg.spawn(|| {
                    ran.fetch_add(1, Ordering::Relaxed);
                });
            });
            // ~TaskGroup cancels first, so with one worker the task never runs
            assert_eq!(ran.load(Ordering::Relaxed), 0);
        });
    }

    #[test]
    fn for_each_covers_the_range_with_the_same_chunks_at_every_thread_count() {
        let mut reference: Option<Vec<(HighsInt, HighsInt)>> = None;
        for n in [1, 2, 4, 8] {
            let chunks = std::sync::Arc::new(Mutex::new(Vec::new()));
            let c2 = chunks.clone();
            on_fresh_thread(n, move || {
                let hits: Vec<AtomicI32> = (0..1000).map(|_| AtomicI32::new(0)).collect();
                for_each(
                    0,
                    1000,
                    &|s, e| {
                        c2.lock().unwrap().push((s, e));
                        for i in s..e {
                            hits[i as usize].fetch_add(1, Ordering::Relaxed);
                        }
                    },
                    10,
                );
                assert!(hits.iter().all(|h| h.load(Ordering::Relaxed) == 1));
            });
            let mut got = chunks.lock().unwrap().clone();
            got.sort();
            match &reference {
                None => reference = Some(got),
                Some(r) => assert_eq!(&got, r, "threads = {n}"),
            }
        }
    }

    #[test]
    fn group_destructor_runs_on_unwind_then_resumes_the_panic() {
        on_fresh_thread(1, || {
            let ran = AtomicI32::new(0);
            let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                task_group(|tg| {
                    tg.spawn(|| {
                        ran.fetch_add(1, Ordering::Relaxed);
                    });
                    std::panic::resume_unwind(Box::new(7_i32));
                })
            }));
            // the foreign payload reaches the caller unchanged, after ~TaskGroup
            assert_eq!(*r.unwrap_err().downcast::<i32>().unwrap(), 7);
            // ~TaskGroup cancelled the task (one worker) and restored the deque head
            assert_eq!(ran.load(Ordering::Relaxed), 0);
            // SAFETY: get_this_worker_deque returns this thread's worker deque.
            assert_eq!(
                unsafe { (*HighsTaskExecutor::get_this_worker_deque()).get_current_head() },
                0
            );
        });
    }

    #[test]
    fn sync_after_overflow_runs_every_task_once() {
        on_fresh_thread(1, || {
            let n = super::super::highs_split_deque::kTaskArraySize as i32 + 5;
            let ran = AtomicI32::new(0);
            task_group(|tg| {
                for _ in 0..n {
                    tg.spawn(|| {
                        ran.fetch_add(1, Ordering::Relaxed);
                    });
                }
                // the first pops see kOverflown (tasks already run in push)
                for _ in 0..n {
                    tg.sync();
                }
            });
            assert_eq!(ran.load(Ordering::Relaxed), n);
        });
    }

    #[test]
    fn for_each_within_grain_size_is_one_direct_call() {
        on_fresh_thread(2, || {
            let chunks = Mutex::new(Vec::new());
            for_each(3, 13, &|s, e| chunks.lock().unwrap().push((s, e)), 10);
            assert_eq!(*chunks.lock().unwrap(), vec![(3, 13)]);
            let chunks = Mutex::new(Vec::new());
            // default grainSize = 1: every chunk has length 1
            for_each(0, 5, &|s, e| chunks.lock().unwrap().push((s, e)), 1);
            let mut got = chunks.lock().unwrap().clone();
            got.sort();
            assert_eq!(got, vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]);
        });
    }

    #[test]
    fn cancelled_parent_interrupts_its_subtasks() {
        on_fresh_thread(4, || {
            let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                task_group(|tg| {
                    tg.spawn(|| {
                        task_group(|inner| {
                            for _ in 0..64 {
                                inner.spawn(|| {
                                    std::thread::sleep(std::time::Duration::from_millis(1))
                                });
                            }
                        })
                    });
                    tg.cancel();
                    tg.task_wait();
                })
            }));
            // either the task ran before the cancel or it was interrupted; never a foreign panic
            if let Err(p) = r {
                assert!(p.downcast_ref::<Interrupt>().is_some());
            }
        });
    }
}
