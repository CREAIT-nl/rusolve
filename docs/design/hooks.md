# MIP hook points and validity gates

Code: `crates/rusolve-ekk/src/mip/hooks.rs`. Tests: `crates/rusolve-ekk/tests/hooks_mip.rs`.

## What was already there: the HiGHS callback system (faithful part)

`lp_data/highs_callback.rs` ports `HighsCallback.{h,cpp}` of HiGHS 1.16 (`HighsCallbackOutput`,
`HighsCallbackInput`, `kCallback*`), installed with `Highs::set_callback` /
`start_callback`, and fired at the same call sites as HiGHS (`mip/mip_solver_data/part2.rs`,
`part3.rs`; HighsMipSolverData.cpp lines 920, 1310, 1564, 1895, 2413, 2603, 2731, 2753):

| callback | data out | data in |
|---|---|---|
| `kCallbackMipSolution`, `kCallbackMipImprovingSolution` | solution, objective, bounds, gap, node and LP counts | none |
| `kCallbackMipInterrupt`, `kCallbackMipLogging` | counts, bounds | `user_interrupt` stops the solve (`kInterrupt`) |
| `kCallbackMipGetCutPool` | the cut pool as CSR rows | none |
| `kCallbackMipUserSolution` | origin of the query, bounds | `user_has_solution`, `user_solution` (original space) |
| `kCallbackMipDefineLazyConstraints` | none | none: defined in HiGHS 1.16 but **never fired** (not fired here either) |

So node-LP reading, interrupt, user solutions and cut-pool reading exist. HiGHS has no way to add
cuts, lazy constraints or influence branching; the hooks below fill that gap.

## rusolve hook points

`trait MipHooks: Send` with default no-op methods, held as `Option<MipHooksHandle>` (an
`Arc<Mutex<Box<dyn MipHooks>>>`, so that the `Clone` of `HighsCallback`/`HighsMipSolver` stays
valid) in `HighsCallback::hooks`; the callback is moved into the solver and back like every other
callback state. Install with `Highs::set_mip_hooks(Box<dyn MipHooks>)`, remove with
`clear_mip_hooks`, read counters with `mip_hook_stats()`.

With no hooks the only cost at a call site is an `Option` test (`user_hooks()`); no RNG draw, no
allocation, no extra work. Hooks are never active in sub-MIPs (`submip`), like the HiGHS callbacks.
Hooks run in the presolved column space (use `presolve=off` for original indices).

| hook | call site | sits beside (HiGHS 1.16) | effect |
|---|---|---|---|
| `branching_score` | top of `HighsSearch::select_branching_candidate` (`mip/search/part1.rs`) | HighsSearch.cpp `selectBranchingCandidate` (247) | one score per fractional candidate, highest wins; skips reliability strong branching at that node |
| `select_node` | `HighsMipSolver::pop_best_node_hooked`, called where `run` pops `popBestNode` (`installNodes`, `prepareNodes`; `mip/mip_solver.rs`) | HighsNodeQueue.cpp `popBestNode` (388) | picks the open node to dive into; best-bound pops stay as HiGHS |
| `cut_score` | `HighsCutPool::separate_hooked` (`mip/cut_pool.rs`), the efficacy score of each violated pool cut; the SCIP-ported hybrid cut selector is not hooked | HighsCutPool.cpp `separate` (229), `efficacious_cuts` | replaces the score used for ordering and the min-score cut-off |
| `separate` | after the built-in separator loop in `HighsSeparation::separation_round` (`mip/separation.rs`) | HighsSeparation.cpp `separationRound` (43), the `for separator` loop | returns cuts `a x <= rhs` (assumed globally valid) for the cut pool; only violated ones reach the LP |

The first three can only make the solver slower, never wrong; their returned indices and scores
are range- and NaN-checked and ignored when invalid.

Open: primal-heuristic hook (the HiGHS user-solution callback covers solution injection).
Lazy constraints are below.

## Lazy constraints

Code: `crates/rusolve-ekk/src/mip/lazy.rs`. Tests: `crates/rusolve-ekk/tests/lazy_cuts.rs`.
Examples: `crates/rusolve-ekk/examples/tsp_lazy.rs`, `examples/user_cuts.rs` (shared code in
`examples/common/mod.rs`).

`trait LazyConstraints { fn lazy_constraints(&mut self, solution: &[f64]) -> Vec<UserCut> }`
(a closure `FnMut(&[f64]) -> Vec<UserCut>` also works) and
`Highs::run_with_lazy_constraints(&mut oracle, max_rounds) -> (HighsStatus, LazyReport)`.
Mechanism: `run`, offer the optimal solution to the oracle, add every returned row that the solution
violates as a model row (`a x <= rhs`, original column space), re-solve; stop when the oracle returns
nothing (the solution is then optimal for the model plus every row ever added), when a solve is not
optimal, or at `max_rounds` (reported as `accepted = false`, never hidden). A non-empty answer in
which every row is malformed is an oracle bug: the loop stops at once with `accepted = false`
(re-solving would give the same solution; `rows_rejected_malformed` says why). `accepted = true`
needs an empty list, or only well-formed rows the solution already satisfies. Lazy
constraints as a callback are a standard MIP-solver feature (in SCIP, the constraint-handler
`enforce` callback); the TSP use case is Dantzig, Fulkerson, Johnson, Operations Research 2 (1954) 393-410.

Design choice: an outer loop over the model, not a fire point inside the tree (HiGHS 1.16 defines
`kCallbackMipDefineLazyConstraints` but never fires it, and the search keeps no way to reject an
integral node LP and re-enter it). Consequences: independent of presolve (rows are in the original
space, unlike `separate` hooks), the search is untouched and the parity check (`tools/bench/identity.sh`) is unaffected, and every
round re-solves the MIP from scratch (no incumbent reuse), so it costs more than in-tree lazy rows
on instances that need many rounds. Revisit with a measured case.

Gates: malformed rows (length mismatch, empty, non-finite, index out of range, duplicate column)
are dropped (`rows_rejected_malformed`), rows the candidate already satisfies are ignored
(`rows_not_violated`, they could not make progress), so every added row removes the current solution
and the loop terminates. A wrong (invalid) lazy row cannot be detected here, as for user cuts.

Evidence (tests): TSP with only degree-2 rows, subtour rows from the oracle, against brute force
over all tours and against the formulation that lists every subtour row up front, on 8 generated
instances of 6 to 10 cities (equal optimum to 1e-6, accepted solution is one tour; instances
n=9 seeds 4 and 5 and n=6 seed 1 needed 2 rounds and 2 subtour rows). User cuts: the 1234 minimal
cover inequalities of a 14-item knapsack as model rows and returned by a `MipHooks::separate`
hook reach the same optimum as the plain model (180), and the hook's cuts reach the pool
(`cargo run --release -p rusolve-ekk --example user_cuts`; the example `tsp_lazy 9 4` prints 2 rounds and
the brute-force optimum 277.019016).

## Validity gates

Anything a hook adds that can cut off solutions is re-checked before it takes effect, and each
rejection is counted in `HookStats` and reported on stderr with the hook name:

- User cuts (`separate`): malformed rows (length mismatch, empty, non-finite, index out of range,
  duplicate column) are rejected (`user_cuts_rejected_malformed`). When the debug-solution harness
  is active (`--opt mip_debug_solution_file=...`, `mip/debug_sol.rs`, HighsDebugSol port), a cut that
  cuts off the debug solution is rejected and reported as
  `HOOK VALIDITY GATE: hook '<name>' separator cut ... cuts off the debug solution`
  (`user_cuts_rejected_debug_solution`); `RUSOLVE_DEBUG_SOL_ABORT=1` panics instead, as in the
  harness. Without a debug solution a cut cannot be proven valid: the gate then only protects
  against malformed rows, so validate new separators on the validity pool with a debug solution.
- User solutions (`kCallbackMipUserSolution`): checked against the original model by the ported
  `solution_feasible` (bounds, integrality, rows, HiGHS tolerances) before `add_incumbent`; a
  failure is logged and counted (`user_solutions_rejected`), a pass is counted
  (`user_solutions_checked_feasible`).

Open: proof certificates (VIPR) for hook-added cuts; hooks cannot add reductions, so there is no
gate for them.

## Python

Code: `crates/rusolve-py/src/callback.rs`. Tests: `crates/rusolve-py/tests/test_callbacks.py`.
Examples: `crates/rusolve-py/examples/{tsp_lazy,user_cuts,prototype_branching}.py` (shared code in
`hooks_common.py`, which reproduces the Rust generators). Full usage: `docs/design/python-api.md`.

`Model.solve(callback=obj)`: any Python object whose methods are hooks. A missing method costs
nothing at run time (looked up once per solve); with no `callback` the solve path is the same
code as before (target P1 of `crates/rusolve-py/bench/bench_api.py`: under 20 us).

| Python method | Rust surface it wraps | returns |
|---|---|---|
| `on_progress(ctx)` | `kCallbackMipInterrupt` (each limit check); `ctx.interrupt()` stops the solve | nothing |
| `on_improving_solution(ctx)` | `kCallbackMipImprovingSolution`: incumbent, objective, bound, gap | nothing |
| `on_cut_pool(ctx)` | `kCallbackMipGetCutPool`: `ctx.cut_pool` rows | nothing |
| `user_solution(ctx)` | `kCallbackMipUserSolution`; checked by the validity gate | full-length list or None |
| `separate(ctx)` | `MipHooks::separate`; node LP solution and bounds in `ctx` | cuts |
| `branching_score(ctx)` | `MipHooks::branching_score`; `ctx.candidates`, `ctx.solution`, bounds | one score per candidate, or None |
| `select_node(ctx)` | `MipHooks::select_node`; `ctx.open_nodes` | index into `ctx.open_nodes`, or None |
| `lazy_constraints(ctx)` | `Highs::run_with_lazy_constraints`; `ctx.solution` is the candidate | violated rows |

A cut is a `Constraint` (`x[0] + x[1] <= 1`), `(indices, values, rhs)` or `(coefficients, rhs)`
(dict or `[(col, val)]`), meaning `a x <= rhs`. Cuts pass the same validity gate as in Rust
(`Model.callback_stats["cuts_rejected"]`).

Decisions:
- **Column space.** `separate`, `branching_score` and `select_node` run in the presolved space
  (see above), so installing one switches presolve off for that solve (recorded in the
  envelope's `run.options`). `lazy_constraints`, `user_solution` and the `on_*` events use the
  model's own columns and leave presolve alone.
- **Lazy rows are per solve.** The loop adds rows to the solver's copy; they are removed after
  the solve, so the Python `Model` stays as built, `Solution.row_activity` has the model's rows
  and the next solve asks the callback again. A loop that ends at `lazy_rounds`, or interrupted,
  gives status `LazyRoundLimit` or `Interrupted` with no `x` (the candidate violates rows the
  callback would add), never a silent "Optimal".
- **Failure.** An exception (or a wrong-type or wrong-length answer, an index out of range) is
  caught inside the solver thread, never unwinds into the search, and is counted. The first one
  is kept; from then on no hook is called (each falls back to the solver's own rule) and the
  solve is interrupted at the next limit check. After the solve `Model.solve` raises
  `rusolve.CallbackError` (code RS4008, location the line that called `solve`), `__cause__` the
  original exception, `.solution` the partial result. Gates stay in force: a malformed cut or an
  infeasible user solution is rejected and counted, not an error.
- **Interrupt.** `ctx.interrupt()` works from every hook (a flag read by the `kCallbackMipInterrupt`
  call at the next limit check); the status is `Interrupted`, the best solution so far is kept.
- **Threads.** Hooks take the GIL per call (`Python::attach`), so they are also safe when the
  solver runs several threads, but they serialise. They are for prototypes: the measured Wasm
  overhead is in `wasm-hooks.md`; promote a stable component to Wasm or Rust.
- **Counters.** `Model.callback_stats`: calls per hook, errors, interrupted, cuts
  offered/accepted/duplicate/rejected, branching and node-selection overrides, user solutions
  accepted/rejected, and the lazy report (rounds, rows added, rows not violated, accepted).

With the default (recommended) profile: the 28-item knapsack of `hooks_mip.rs` is solved
without any queue pop (`select_node` is never consulted); a 34-item 5-row one pops (2 calls), so
the tests use that instance for `select_node`. The `separate` hook fires at the root only when
the heuristics are off and the gap is closed (`hooks_common.CUT_OPTIONS`): with defaults the
feasibility jump finds 180 and the root bound is 180 before the first separation round.

Evidence (`pytest crates/rusolve-py/tests`): lazy TSP equals brute force (Held-Karp)
on the 8 instances of `lazy_cuts.rs` (n=9 seed 4: 2 rounds, 2 rows, 277.019016, as in Rust); the
Python cover separator, the same 1234 covers as model rows, and the plain solve all give 180.0
(2460 cuts accepted); `ctx.interrupt()` ends with `Interrupted`; every hook kind that raises or
answers wrongly gives RS4008.

Open: in-tree lazy rows with incumbent reuse; a Python `cut_score` hook; Python primal-heuristic
hook.
