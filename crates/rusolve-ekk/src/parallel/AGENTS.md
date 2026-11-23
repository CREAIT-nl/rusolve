# parallel

The HiGHS task scheduler: work-stealing deque, task executor, mutexes, semaphores, cache alignment and the parallel loop helpers the simplex (PAMI), MIP and LP code call.

## What lives here
- `highs_parallel.rs` (`HighsParallel.h`: spawn, sync, parallel_for), `highs_task_executor.rs`, `highs_task.rs`, `highs_split_deque/`.
- Synchronisation: `highs_mutex.rs`, `highs_spin_mutex.rs`, `highs_binary_semaphore.rs`, `highs_combinable.rs`, `highs_cache_align.rs`, `highs_race_timer.rs`, `highs_scheduler_constants.rs`.
- `shared_mut.rs`: the rusolve-side cell for data the C++ shares by raw pointer.

## HiGHS files mirrored
`parallel/*.h` and `HighsTaskExecutor.cpp`. Table in `docs/agent/codemap.md`.

## Invariants
- The default run is single-threaded and deterministic; threaded results must equal the serial ones where HiGHS guarantees it.
- `unsafe` sits in the lock-free structures (`highs_split_deque/`, mutexes, tasks) and in `shared_mut.rs`, where every access carries a `// SAFETY:` comment naming the C++ protocol it follows. New `unsafe` needs the same comment and a justification in the commit message.
- Thread counts come from the `threads` option.

## How to test
- `cargo test --release -p rusolve-ekk --test t_w5_threads --test t_w5_pami`; unit tests in the files (`--lib parallel`).
- CLI: `crates/rusolve-cli/tests/mip_threads.rs`, `crates/rusolve-cli/tests/parallel_clique_race.rs`.
- A/B arms must run on the same host.

## Change recipes
1. Port a scheduler piece: `docs/skills/contrib-port-function/SKILL.md`; keep the C++ memory orderings in the translation.
2. Parallelize a loop: call `highs_parallel.rs` helpers, keep the serial path as the default, add a threads test.
3. A race shows up: reproduce with `threads=2` on the build host, find the C++ protocol in the `// SAFETY:` comment of the access, and fix the access, not the test.

## Pitfalls
- Spawn, sync and `for_each` are in `highs_parallel.rs`; the executor is reference counted by hand and reached through a thread-local `ExecutorHandle` (`highs_task_executor.rs`).
- `HighsCombinable::combine_each` visits copies in worker-id order, but which workers hold a copy depends on stealing: a combine must be order- and partition-independent to stay deterministic (`highs_combinable.rs`).
