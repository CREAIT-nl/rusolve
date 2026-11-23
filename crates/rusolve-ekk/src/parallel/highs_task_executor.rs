//! HiGHS source: parallel/HighsTaskExecutor.h, parallel/HighsTaskExecutor.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The executor: `numThreads` deques, `numThreads - 1` worker threads (worker 0 is the
//! thread that called `initialize`), one shared `WorkerBunk`.
//!
//! **Thread locals.** C++ (non-Windows) has `thread_local HighsSplitDeque*
//! threadLocalWorkerDequePtr` and `thread_local ExecutorHandle globalExecutorHandle`
//! (HighsTaskExecutor.cpp 25-28): every user thread that calls `initialize` owns its own
//! executor. Here they are `thread_local!` cells; the handle's `Drop` is the C++
//! `~ExecutorHandle` (`dispose`). The executor is reference counted by hand, as C++
//! (`referenceCount` = `numThreads`, one handle per worker), and freed by the last
//! `dispose`.
//!
//! **No executor.** C++ `getThisWorkerDeque()` returns null on a thread that never called
//! `initialize` and callers dereference it. The port instead initialises a one-worker
//! executor there (as `initialize(1)`), so library code reached without a `Highs::run`
//! (unit tests, direct solver calls) behaves exactly as C++ at `threads = 1`.
//!
//! **Building the deques.** C++ sizes `workerDeques` first and passes `workerDeques.data()`
//! to every deque constructor. The port reserves the `Vec` and writes the boxes through
//! `as_mut_ptr()` before `set_len`, so the address handed to the deques never moves.

use super::highs_cache_align::CacheAligned;
use super::highs_scheduler_constants::{
    kMicroSecsBeforeGlobalSync, kMicroSecsBeforeSleep, kNumTryFac,
};
use super::highs_split_deque::{HighsSplitDeque, WorkerBunk};
use super::highs_task::HighsTask;
use std::cell::{Cell, RefCell};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicI32, Ordering};
use std::thread::JoinHandle;

/// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::ExecutorHandle (lines 27-33) with
/// HighsTaskExecutor.cpp::ExecutorHandle::dispose (lines 31-43).
#[derive(Debug)]
pub struct ExecutorHandle {
    pub ptr: *mut HighsTaskExecutor,
    pub is_main: bool,
}

impl Default for ExecutorHandle {
    fn default() -> Self {
        ExecutorHandle {
            ptr: std::ptr::null_mut(),
            is_main: false,
        }
    }
}

impl ExecutorHandle {
    /// HiGHS: HighsTaskExecutor.cpp::HighsTaskExecutor::ExecutorHandle::dispose (lines 31-43)
    pub fn dispose(&mut self) {
        // C++ line 32: if (ptr == nullptr) return;
        if self.ptr.is_null() {
            return;
        }

        // C++ line 33-35: if (isMain) { ptr->stopWorkerThreads(false); }
        if self.is_main {
            // SAFETY: self.ptr is non-null and points to a live HighsTaskExecutor
            // (the handle owns a reference to it).
            unsafe {
                (*self.ptr).stop_worker_threads_on(false, self.is_main);
            }
        }

        // C++ line 38: if (--ptr->referenceCount == 0) { cache_aligned::Deleter<HighsTaskExecutor>()(ptr); }
        // C++ --referenceCount is a SeqCst decrement (default for std::atomic<int>).
        // SAFETY: self.ptr is non-null and points to a live HighsTaskExecutor.
        let old = unsafe { (*self.ptr).reference_count.fetch_sub(1, Ordering::SeqCst) };
        if old - 1 == 0 {
            // C++ line 39: cache_aligned::Deleter<HighsTaskExecutor>()(ptr);
            // This maps to Box::from_raw on the Box<CacheAligned<HighsTaskExecutor>>.
            // SAFETY: self.ptr was obtained from Box::into_raw on a Box<CacheAligned<HighsTaskExecutor>>
            // and the reference count has reached zero, so no other handle references it.
            unsafe {
                let _ = Box::from_raw(self.ptr as *mut CacheAligned<HighsTaskExecutor>);
            }
        }

        // C++ line 42: ptr = nullptr;
        self.ptr = std::ptr::null_mut();
    }
}

impl Drop for ExecutorHandle {
    /// C++ `~ExecutorHandle() { dispose(); }`
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            self.dispose();
        }
    }
}

thread_local! {
    /// HiGHS: HighsTaskExecutor.cpp `thread_local HighsSplitDeque* threadLocalWorkerDequePtr` (lines 25-26)
    static THREAD_LOCAL_WORKER_DEQUE: Cell<*const HighsSplitDeque> = const { Cell::new(std::ptr::null()) };
    /// HiGHS: HighsTaskExecutor.cpp `thread_local ExecutorHandle globalExecutorHandle` (lines 27-28)
    static GLOBAL_EXECUTOR_HANDLE: RefCell<ExecutorHandle> = RefCell::new(ExecutorHandle::default());
}

/// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor (lines 24-215)
#[derive(Debug)]
pub struct HighsTaskExecutor {
    /// C++ `std::atomic<int> referenceCount`
    pub(crate) reference_count: AtomicI32,
    /// C++ `std::atomic<bool> hasStopped{false}`
    pub(crate) has_stopped: AtomicBool,
    /// C++ `cache_aligned::shared_ptr<HighsSplitDeque::WorkerBunk> workerBunk`
    pub(crate) worker_bunk: Arc<WorkerBunk>,
    /// C++ `std::vector<cache_aligned::unique_ptr<HighsSplitDeque>> workerDeques`. Boxed as
    /// in C++: other workers hold the deques' addresses (rule (a)).
    #[allow(clippy::vec_box)]
    pub(crate) worker_deques: Vec<Box<HighsSplitDeque>>,
    /// C++ `std::vector<std::thread> workerThreads`
    pub(crate) worker_threads: Vec<JoinHandle<()>>,
}

// SAFETY: shared by its worker threads through `ExecutorHandle::ptr`, as in C++; mutable
// state is atomic or touched only by the main thread (`worker_threads` in
// `stopWorkerThreads`).
unsafe impl Send for HighsTaskExecutor {}
unsafe impl Sync for HighsTaskExecutor {}

impl HighsTaskExecutor {
    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::random_steal_loop (lines 59-86)
    fn random_steal_loop(&self, local_deque: *const HighsSplitDeque) -> *mut HighsTask {
        // C++ line 60: const int numWorkers = static_cast<int>(workerDeques.size());
        let num_workers = self.worker_deques.len() as i32;

        // C++ line 62: int numTries = 16 * (numWorkers - 1);
        let mut num_tries = 16 * (num_workers - 1);

        // C++ line 64: auto tStart = std::chrono::high_resolution_clock::now();
        let t_start = std::time::Instant::now();

        // C++ line 66: while (true) {
        loop {
            // C++ line 67: for (int s = 0; s < numTries; ++s) {
            for _s in 0..num_tries {
                // C++ line 68: HighsTask* task = localDeque->randomSteal();
                // SAFETY: local_deque is a valid pointer to the calling worker's deque.
                let task = unsafe { (*local_deque).random_steal() };
                // C++ line 69: if (task) return task;
                if !task.is_null() {
                    return task;
                }
            }

            // C++ line 72: if (!workerBunk->haveJobs.load(std::memory_order_relaxed)) break;
            // C++ haveJobs is an int; !int means == 0.
            if self.worker_bunk.have_jobs.load(Ordering::Relaxed) == 0 {
                break;
            }

            // C++ line 74-77: auto numMicroSecs = std::chrono::duration_cast<std::chrono::microseconds>(
            //     std::chrono::high_resolution_clock::now() - tStart).count();
            let num_micro_secs = {
                let elapsed = t_start.elapsed();
                elapsed.as_micros() as i64
            };

            // C++ line 79-82: if (numMicroSecs < HighsSchedulerConstants::kMicroSecsBeforeGlobalSync)
            //     numTries *= 2;
            // else
            //     break;
            if num_micro_secs < kMicroSecsBeforeGlobalSync {
                num_tries *= 2;
            } else {
                break;
            }
        }

        // C++ line 85: return nullptr;
        std::ptr::null_mut()
    }

    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::run_worker (lines 88-108). `ptr` is
    /// passed as `usize` so the closure handed to `std::thread::spawn` is `Send`.
    fn run_worker(worker_id: i32, ptr: usize) {
        // C++ line 89: auto& executorHandle = threadLocalExecutorHandle();
        // C++ line 90: executorHandle.ptr = ptr;
        GLOBAL_EXECUTOR_HANDLE.with(|handle| {
            handle.borrow_mut().ptr = ptr as *mut HighsTaskExecutor;
        });

        // SAFETY: ptr is a valid pointer to a live HighsTaskExecutor (the main thread
        // keeps it alive until all workers have disposed their handles).
        let executor = unsafe { &*(ptr as *const HighsTaskExecutor) };

        // C++ line 92: if (!ptr->hasStopped) {
        if !executor.has_stopped.load(Ordering::SeqCst) {
            // C++ line 93: HighsSplitDeque* localDeque = ptr->workerDeques[workerId].get();
            // SAFETY: worker_id is in [1, numThreads-1] and worker_deques has numThreads entries.
            let local_deque: *const HighsSplitDeque =
                &*executor.worker_deques[worker_id as usize] as *const HighsSplitDeque;

            // C++ line 94: threadLocalWorkerDeque() = localDeque;
            THREAD_LOCAL_WORKER_DEQUE.with(|cell| {
                cell.set(local_deque);
            });

            // C++ line 96: HighsTask* currentTask = ptr->workerBunk->waitForNewTask(localDeque);
            let mut current_task = executor.worker_bunk.wait_for_new_task(local_deque);

            // C++ line 97: while (currentTask != nullptr) {
            while !current_task.is_null() {
                // C++ line 98: localDeque->runStolenTask(currentTask);
                // SAFETY: current_task is a valid task pointer returned by waitForNewTask;
                // local_deque is the calling worker's deque.
                unsafe {
                    (*local_deque).run_stolen_task(current_task);
                }

                // C++ line 100: currentTask = ptr->random_steal_loop(localDeque);
                current_task = executor.random_steal_loop(local_deque);

                // C++ line 101: if (currentTask != nullptr) continue;
                if !current_task.is_null() {
                    continue;
                }

                // C++ line 103: currentTask = ptr->workerBunk->waitForNewTask(localDeque);
                current_task = executor.worker_bunk.wait_for_new_task(local_deque);
            }
        }

        // C++ line 107: executorHandle.dispose();
        GLOBAL_EXECUTOR_HANDLE.with(|handle| {
            handle.borrow_mut().dispose();
        });
    }

    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::HighsTaskExecutor (lines 111-128).
    /// Returns the executor boxed (C++ `cache_aligned::alloc` + placement new in
    /// `initialize`): the worker threads hold its address.
    pub fn new(num_threads: i32) -> Box<CacheAligned<HighsTaskExecutor>> {
        // C++ line 112: assert(numThreads > 0);
        debug_assert!(num_threads > 0);

        // C++ line 114: workerDeques.resize(numThreads);
        // We reserve and then write through as_mut_ptr to keep addresses stable.
        let mut worker_deques: Vec<Box<HighsSplitDeque>> = Vec::with_capacity(num_threads as usize);

        // C++ line 115: workerBunk = cache_aligned::make_shared<HighsSplitDeque::WorkerBunk>();
        let worker_bunk = Arc::new(WorkerBunk::default());

        // C++ line 116-118: workerDeques[i] = make_unique<HighsSplitDeque>(workerBunk,
        // workerDeques.data(), i, numThreads). The buffer is reserved first so that its
        // address is fixed before any deque is built; it is never reallocated afterwards.
        let workers_ptr: *const Box<HighsSplitDeque> = worker_deques.as_ptr();

        // Now construct each deque and write it into the buffer.
        for i in 0..num_threads {
            let deque = Box::new(HighsSplitDeque::new(
                &worker_bunk,
                workers_ptr,
                i,
                num_threads,
            ));
            // Write into the vec's buffer at position i.
            // SAFETY: i < capacity, and we're writing to uninitialized memory
            // within the reserved capacity.
            unsafe {
                std::ptr::write(worker_deques.as_mut_ptr().add(i as usize), deque);
            }
        }
        // SAFETY: We've written num_threads elements into the buffer.
        unsafe {
            worker_deques.set_len(num_threads as usize);
        }

        // C++ line 120: threadLocalWorkerDeque() = workerDeques[0].get();
        let deque_0_ptr: *const HighsSplitDeque = &*worker_deques[0] as *const HighsSplitDeque;
        THREAD_LOCAL_WORKER_DEQUE.with(|cell| {
            cell.set(deque_0_ptr);
        });

        // C++ line 121: workerThreads.reserve(numThreads - 1);
        let worker_threads: Vec<JoinHandle<()>> = Vec::with_capacity((num_threads - 1) as usize);

        // C++ line 122: referenceCount.store(numThreads);
        let reference_count = AtomicI32::new(num_threads);

        // C++ `this` is the placement-new address from `initialize`; here the executor is
        // boxed first so the worker threads (C++ lines 124-127, below) get its final address.
        let executor = HighsTaskExecutor {
            reference_count,
            has_stopped: AtomicBool::new(false),
            worker_bunk,
            worker_deques,
            worker_threads,
        };

        let mut boxed = Box::new(CacheAligned::new(executor));
        let executor_ptr: usize = &mut boxed.0 as *mut HighsTaskExecutor as usize;

        // C++ line 124-127: for (int i = 1, numThreads = static_cast<int>(workerDeques.size());
        //      i < numThreads; ++i)
        //   workerThreads.emplace_back(&HighsTaskExecutor::run_worker, i, this);
        let num_threads_actual = boxed.0.worker_deques.len() as i32;
        for i in 1..num_threads_actual {
            let handle = std::thread::spawn(move || {
                HighsTaskExecutor::run_worker(i, executor_ptr);
            });
            boxed.0.worker_threads.push(handle);
        }

        boxed
    }

    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::stopWorkerThreads (lines 130-150)
    pub fn stop_worker_threads(&mut self, blocking: bool) {
        // C++ line 132: auto& executorHandle = threadLocalExecutorHandle();
        // C++ line 133: if (executorHandle.ptr == nullptr || hasStopped.exchange(true)) return;
        // The handle is read and the borrow released before anything else runs, so that
        // `dispose` (which calls back into this function) never sees a live borrow.
        let (handle_ptr, is_main) = GLOBAL_EXECUTOR_HANDLE.with(|handle| {
            let h = handle.borrow();
            (h.ptr, h.is_main)
        });
        if handle_ptr.is_null() {
            return;
        }
        self.stop_worker_threads_on(blocking, is_main);
    }

    /// C++ `stopWorkerThreads` body (lines 133-149) with the calling thread's handle state
    /// (`executorHandle.isMain`) passed in. Called from `ExecutorHandle::dispose`, which
    /// already holds the handle, so it must not re-enter the thread-local.
    fn stop_worker_threads_on(&mut self, blocking: bool, is_main: bool) {
        // C++ hasStopped.exchange(true) is SeqCst (default).
        if self.has_stopped.swap(true, Ordering::SeqCst) {
            return;
        }

        // C++ line 136-138: for (auto& workerDeque : workerDeques) {
        //   workerDeque->injectTaskAndNotify(nullptr);
        // }
        for worker_deque in &self.worker_deques {
            worker_deque.inject_task_and_notify(std::ptr::null_mut());
        }

        // C++ line 141-149: if (blocking && executorHandle.isMain) {
        //   for (auto& workerThread : workerThreads) { workerThread.join(); }
        // } else {
        //   for (auto& workerThread : workerThreads) { workerThread.detach(); }
        // }
        if blocking && is_main {
            for worker_thread in self.worker_threads.drain(..) {
                let _ = worker_thread.join();
            }
        } else {
            // C++ workerThread.detach() drops the std::thread without joining.
            // In Rust, dropping a JoinHandle detaches the thread.
            for worker_thread in self.worker_threads.drain(..) {
                drop(worker_thread);
            }
        }
    }

    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::getThisWorkerDeque (lines 152-154).
    /// Never null (see the module doc, "No executor").
    pub fn get_this_worker_deque() -> *const HighsSplitDeque {
        // C++ line 153: return threadLocalWorkerDeque();
        // The module doc says: if null, initialize a one-worker executor.
        let deque = THREAD_LOCAL_WORKER_DEQUE.with(|cell| cell.get());
        if deque.is_null() {
            // Initialize a one-worker executor on this thread.
            HighsTaskExecutor::initialize(1);
            THREAD_LOCAL_WORKER_DEQUE.with(|cell| cell.get())
        } else {
            deque
        }
    }

    /// rusolve: true on a thread that is the main thread of its executor (owner id 0) or has no
    /// executor yet. Debug checks only.
    pub fn this_thread_is_main() -> bool {
        let deque = THREAD_LOCAL_WORKER_DEQUE.with(|cell| cell.get());
        // SAFETY: a non-null thread-local points at a live deque (shutdown nulls it).
        deque.is_null() || unsafe { (*deque).get_owner_id() == 0 }
    }

    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::getNumWorkerThreads (lines 156-158)
    pub fn get_num_worker_threads() -> i32 {
        // C++ line 157: return threadLocalWorkerDeque()->getNumWorkers();
        let deque = HighsTaskExecutor::get_this_worker_deque();
        // SAFETY: get_this_worker_deque never returns null (it initializes if needed).
        unsafe { (*deque).get_num_workers() }
    }

    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::initialize (lines 160-167)
    pub fn initialize(num_threads: i32) {
        // C++ line 161: auto& executorHandle = threadLocalExecutorHandle();
        // C++ line 162: if (executorHandle.ptr == nullptr) {
        let already_initialized =
            GLOBAL_EXECUTOR_HANDLE.with(|handle| !handle.borrow().ptr.is_null());
        if already_initialized {
            return;
        }

        // C++ line 163: executorHandle.isMain = true;
        // C++ line 164-165: executorHandle.ptr = new (cache_aligned::alloc(sizeof(HighsTaskExecutor)))
        //     HighsTaskExecutor(numThreads);
        let boxed = HighsTaskExecutor::new(num_threads);
        let raw = Box::into_raw(boxed);

        GLOBAL_EXECUTOR_HANDLE.with(|handle| {
            let mut h = handle.borrow_mut();
            h.is_main = true;
            h.ptr = raw as *mut HighsTaskExecutor;
        });
    }

    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::shutdown (lines 171-178)
    pub fn shutdown(blocking: bool) {
        // C++ line 172: auto& executorHandle = threadLocalExecutorHandle();
        // C++ line 174: if (executorHandle.ptr != nullptr) {
        let ptr = GLOBAL_EXECUTOR_HANDLE.with(|handle| handle.borrow().ptr);
        if ptr.is_null() {
            return;
        }

        // C++ line 175: executorHandle.ptr->stopWorkerThreads(blocking);
        // SAFETY: ptr is non-null and points to a live HighsTaskExecutor.
        unsafe {
            (*ptr).stop_worker_threads(blocking);
        }

        // C++ line 176: executorHandle.dispose();
        GLOBAL_EXECUTOR_HANDLE.with(|handle| {
            handle.borrow_mut().dispose();
        });
        // rusolve: C++ leaves `threadLocalWorkerDeque()` pointing into the freed executor and
        // never reads it before `initializeMultiThreading`. Here a caller that reads it in
        // between (a later run on this thread) gets a one-worker executor instead of a
        // dangling pointer.
        THREAD_LOCAL_WORKER_DEQUE.with(|cell| cell.set(std::ptr::null()));
    }

    /// HiGHS: HighsTaskExecutor.h::HighsTaskExecutor::sync_stolen_task (lines 180-214)
    ///
    /// # Safety
    /// Called by the owner thread of `local_deque`; `stolen_task` is the slot `pop`
    /// reported `kStolen`.
    pub unsafe fn sync_stolen_task(
        local_deque: *const HighsSplitDeque,
        stolen_task: *mut HighsTask,
    ) {
        // SAFETY: sync_stolen_task is called by the owner thread of local_deque with
        // stolen_task being the slot pop reported kStolen; stealer is an out-parameter
        // as in C++ HighsSplitDeque* stealer.
        unsafe {
            // C++ line 182: HighsSplitDeque* stealer;
            let mut stealer: *const HighsSplitDeque = std::ptr::null();

            // C++ line 183: if (!localDeque->leapfrogStolenTask(stolenTask, stealer)) {
            if !(*local_deque).leapfrog_stolen_task(stolen_task, &mut stealer) {
                // C++ line 184: const int numWorkers = localDeque->getNumWorkers();
                let num_workers = (*local_deque).get_num_workers();

                // C++ line 185: int numTries = HighsSchedulerConstants::kNumTryFac * (numWorkers - 1);
                let mut num_tries = kNumTryFac * (num_workers - 1);

                // C++ line 187: auto tStart = std::chrono::high_resolution_clock::now();
                let t_start = std::time::Instant::now();

                // C++ line 189: while (true) {
                loop {
                    // C++ line 190: for (int s = 0; s < numTries; ++s) {
                    for _s in 0..num_tries {
                        // C++ line 191: if (stolenTask->isFinished()) {
                        if (*stolen_task).is_finished() {
                            // C++ line 192: localDeque->popStolen();
                            (*local_deque).pop_stolen();
                            // C++ line 193: return;
                            return;
                        }
                        // C++ line 195: localDeque->yield();
                        (*local_deque).yield_();
                    }

                    // C++ line 198-201: auto numMicroSecs = std::chrono::duration_cast<std::chrono::microseconds>(
                    //     std::chrono::high_resolution_clock::now() - tStart).count();
                    let num_micro_secs = {
                        let elapsed = t_start.elapsed();
                        elapsed.as_micros() as i64
                    };

                    // C++ line 203-206: if (numMicroSecs < HighsSchedulerConstants::kMicroSecsBeforeSleep)
                    //     numTries *= 2;
                    // else
                    //     break;
                    if num_micro_secs < kMicroSecsBeforeSleep {
                        num_tries *= 2;
                    } else {
                        break;
                    }
                }

                // C++ line 209-210: if (!stolenTask->isFinished())
                //     localDeque->waitForTaskToFinish(stolenTask, stealer);
                if !(*stolen_task).is_finished() {
                    (*local_deque).wait_for_task_to_finish(stolen_task, stealer);
                }
            }

            // C++ line 213: localDeque->popStolen();
            (*local_deque).pop_stolen();
        }
    }
}

#[cfg(test)]
mod tests {
    // Unit tests
    use super::*;

    #[test]
    fn initialize_is_per_thread_and_first_call_wins() {
        std::thread::spawn(|| {
            HighsTaskExecutor::initialize(4);
            assert_eq!(HighsTaskExecutor::get_num_worker_threads(), 4);
            assert_eq!(
                unsafe { (*HighsTaskExecutor::get_this_worker_deque()).get_owner_id() },
                0
            );
            // a second initialize is ignored while the handle is set
            HighsTaskExecutor::initialize(2);
            assert_eq!(HighsTaskExecutor::get_num_worker_threads(), 4);
            HighsTaskExecutor::shutdown(true);
            HighsTaskExecutor::initialize(2);
            assert_eq!(HighsTaskExecutor::get_num_worker_threads(), 2);
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
        // another thread is unaffected: no executor there yet, so it gets one worker
        std::thread::spawn(|| assert_eq!(HighsTaskExecutor::get_num_worker_threads(), 1))
            .join()
            .unwrap();
    }

    #[test]
    fn repeated_initialize_shutdown_cycles_join_all_workers() {
        std::thread::spawn(|| {
            for n in [8, 1, 2, 4, 8] {
                HighsTaskExecutor::initialize(n);
                assert_eq!(HighsTaskExecutor::get_num_worker_threads(), n);
                HighsTaskExecutor::shutdown(true);
            }
        })
        .join()
        .unwrap();
    }
}
