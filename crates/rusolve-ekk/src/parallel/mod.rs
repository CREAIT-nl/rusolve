//! HiGHS source: parallel/ (HighsParallel.h, HighsTaskExecutor.{h,cpp}, HighsSplitDeque.h,
//! HighsTask.h, HighsBinarySemaphore.h, HighsSpinMutex.h, HighsMutex.h, HighsRaceTimer.h,
//! HighsCombinable.h, HighsCacheAlign.h, HighsSchedulerConstants.h)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiGHS work-stealing task scheduler. Only `std` is used (no rayon, no
//! crossbeam): the scheduler decides *which thread* runs a task and never *what* a task
//! computes, and HiGHS' own scheduler has properties the port depends on (LIFO `sync`,
//! the deque head as the task-group marker, `getOwnerId` as the `HighsCombinable` slot,
//! the 8192-slot overflow that runs a task inline).
//!
//! # Module map
//!
//! | C++ | Rust |
//! |---|---|
//! | `HighsSchedulerConstants.h` | [`highs_scheduler_constants`] (complete) |
//! | `HighsCacheAlign.h` | [`highs_cache_align`] (complete) |
//! | `HighsSpinMutex.h`, `HighsBinarySemaphore.h`, `HighsRaceTimer.h` | [`highs_spin_mutex`], [`highs_binary_semaphore`], [`highs_race_timer`] |
//! | `HighsTask.h` | [`highs_task`] |
//! | `HighsSplitDeque.h` 49-94, 193-466, 603-612 (owner side, steal) | [`highs_split_deque`] |
//! | `HighsSplitDeque.h` 97-191 (`WorkerBunk`) | [`highs_split_deque::worker_bunk`] |
//! | `HighsSplitDeque.h` 468-601 (stolen-task protocol) | [`highs_split_deque::stolen`] |
//! | `HighsTaskExecutor.{h,cpp}` | [`highs_task_executor`] |
//! | `HighsParallel.h` | [`highs_parallel`] (the `highs::parallel` namespace) |
//! | `HighsMutex.h` | [`highs_mutex`] |
//! | `HighsCombinable.h` | [`highs_combinable`] |
//! | ours: the aliasing escape hatch for MIP workers | [`shared_mut`] (complete) |
//!
//! # Rust mapping
//!
//! * **Pointers.** `HighsSplitDeque*` / `HighsTask*` are `*const HighsSplitDeque` /
//!   `*mut HighsTask`. The executor owns its deques in `Box`es whose addresses never
//!   change while it lives, exactly as `cache_aligned::unique_ptr` in C++; atomics hold the
//!   addresses (`AtomicUsize` / `AtomicPtr`). Owner-only state sits in an `UnsafeCell`
//!   (C++ plain members); the types are `unsafe impl Send + Sync` with the C++ access
//!   protocol as the safety argument.
//! * **Tasks.** `HighsTask::Callable<F>` (inline, at most 56 bytes, trivially
//!   destructible) is a `Box<dyn FnOnce() + Send>` whose lifetime is erased to `'static`
//!   by [`highs_parallel::TaskGroup::spawn`]. The box is storage only: no result depends
//!   on it.
//! * **Interrupt.** `throw HighsTask::Interrupt()` is
//!   `std::panic::resume_unwind(Box::new(Interrupt))` (no panic hook, no message) and
//!   `catch (HighsTask::Interrupt&)` is `catch_unwind` plus a downcast; any other payload
//!   is re-raised.
//! * **Task groups.** `TaskGroup` is scoped: [`highs_parallel::task_group`] builds it, runs
//!   the caller's closure and then runs the C++ destructor (`cancel`, `setRootTask(null)`,
//!   `taskWait`, restore, `checkInterrupt`) also when the closure unwinds, so a spawned
//!   closure may borrow the caller's stack (like `std::thread::scope`).
//! * **Thread-local executor.** `threadLocalWorkerDequePtr` / `globalExecutorHandle` are
//!   `thread_local!`s, as in C++ (one scheduler per user thread; each test thread of the
//!   harness has its own). A thread that never called `initialize` gets a one-worker
//!   executor on first use of `getThisWorkerDeque` (C++ dereferences null there); with one
//!   worker every C++ code path is the serial one.

pub mod highs_binary_semaphore;
pub mod highs_cache_align;
pub mod highs_combinable;
pub mod highs_mutex;
pub mod highs_parallel;
pub mod highs_race_timer;
pub mod highs_scheduler_constants;
pub mod highs_spin_mutex;
pub mod highs_split_deque;
pub mod highs_task;
pub mod highs_task_executor;
pub mod shared_mut;
