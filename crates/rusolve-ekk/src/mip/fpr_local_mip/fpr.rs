//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! precompute_var_orders, FprWorker, run (fpr.cpp:191-546).

use super::context::run_opportunistic_loop;
use super::fpr_core::{fpr_attempt_begin, fpr_attempt_finish, fpr_attempt_step};
use super::fpr_var_order::compute_var_order;
use super::rng::{S5Rng, heuristic_base_seed};
use super::sink::S5Sink;
use super::*;
use crate::mip::mip_solver::HighsMipSolver;
use crate::prelude::*;
use crate::util::highs_utils::{cpp_max, cpp_min};

// ---------------------------------------------------------------------------
// precompute_var_orders (anonymous namespace in C++)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr.cpp::precompute_var_orders (lines 211-223)
///
/// Compute variable orders for every strategy in kFprStrategies. MUST be
/// called from a sequential context: clique-based var_strategies read the
/// live HighsCliqueTable, which `addIncumbent`'s `extractObjCliques`
/// reallocates, and the live root domain.
///
/// Returns false when the solve's wall-clock deadline passed part-way
/// through, leaving `orders` incomplete (issue #117). This setup is the
/// single largest deadline-unaware unit in an FPR dispatch and it is not
/// small: eight `compute_var_order` calls, several of them a clique-cover
/// pass over the whole model. The 34.5 s on `rail02` below was measured
/// against the retired `HighsCliqueTable::cliquePartition`, not against
/// `clique_cover::build_clique_cover` (#141), which is a single linear
/// pass and plausibly much cheaper; unmeasured at that scale, so read it
/// as an upper bound of unknown tightness. Whatever it costs, the
/// dispatch's own budget — and therefore every gate derived from it — has
/// no bearing on it whatsoever. A caller that gets `false`
/// must not use the table; `fpr::run` returns instead, which is also why
/// the deadline is checked *between* strategies rather than inside the
/// loop body: one `compute_var_order` is indivisible here, and it is the
/// residual floor.
pub(crate) fn precompute_var_orders(
    ms: &HighsMipSolver,
    deadline: &Deadline,
    orders: &mut Vec<Vec<HighsInt>>,
) -> bool {
    *orders = vec![Vec::new(); kNumFprStrategies as usize];
    let base = heuristic_base_seed(ms.options_mip_.random_seed);
    for i in 0..kNumFprStrategies {
        if deadline.expired() {
            return false;
        }
        let mut rng = S5Rng::new(base.wrapping_add(i as u32) as u64);
        orders[i as usize] =
            compute_var_order(ms, kFprStrategies[i as usize].var_strategy, &mut rng, None);
    }
    true
}

// ---------------------------------------------------------------------------
// FprWorker implementation
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr.cpp::FprWorker (lines 231-488)
pub(crate) struct FprWorker<'p> {
    problem: &'p S5Problem,
    exec: ExecutionContext,
    var_orders: &'p [Vec<HighsInt>],
    worker_idx: i32,
    trace: WorkerTrace,
    rng: S5Rng,
    base: WorkerBudgetState,
    attempt_idx: i32,
    strat_idx: i32,
    mode: FrameworkMode,
    attempt_state: FprAttemptState,
    scratch: FprScratch<'p>,
}

impl<'p> FprWorker<'p> {
    /// mip-heuristics (36f7d9d): fpr.cpp::FprWorker::FprWorker (lines 231-244)
    pub(crate) fn new(
        problem: &'p S5Problem,
        exec: &ExecutionContext,
        var_orders: &'p [Vec<HighsInt>],
        worker_idx: i32,
        seed: u32,
        stale_budget: usize,
    ) -> Self {
        let mut w = Self {
            problem,
            exec: *exec,
            var_orders,
            worker_idx,
            trace: WorkerTrace {
                worker: worker_idx,
                effort_base: 0,
            },
            rng: S5Rng::new(seed as u64),
            base: WorkerBudgetState {
                total_budget: usize::MAX,
                stale_budget,
                ..Default::default()
            },
            attempt_idx: 0,
            strat_idx: 0,
            mode: FrameworkMode::kDfs,
            attempt_state: FprAttemptState::default(),
            scratch: FprScratch::default(),
        };
        w.select_config_for_current_attempt();
        w
    }

    /// mip-heuristics (36f7d9d): fpr.cpp::FprWorker::select_config_for_current_attempt (lines 246-279)
    fn select_config_for_current_attempt(&mut self) {
        // Per-worker rotation through the paper-curated `kInitialFprConfigs`
        // list (8 entries), keyed deterministically on
        // `(worker_idx + attempt_idx) % kNumInitialFprConfigs`. Each worker
        // visits every Class-1 config exactly once before wrapping. Issue #77's
        // determinism rule is satisfied because the rotation is purely a
        // function of (worker identity, attempt count) — no shared queue, no
        // rng dependency, no per-attempt randomisation.
        //
        // Why not the full 8 × 5 = 40-pair (strategy, mode) grid? An earlier
        // draft widened to it once `e_pq_mark` threading was in place, but
        // the `(kStratDomsize, kRepairSearch)` pairing exposed a second
        // latent state-restoration gap in `repair_search`'s secondary
        // backtrack: `act_mark` is not threaded through `RepairSearchNode`
        // analogously to `e_pq_mark`, so when `init_activities()` ran in
        // Phase 2 (any `kLoosedyn` value strategy) the activity vectors and
        // `vs_` diverge across the secondary backtrack. `kStratDomsize` is
        // the only entry that simultaneously uses `init_domain_pq` AND a
        // `kLoosedyn` val strategy AND was widened to a `kRepairSearch`
        // mode the curated list never exercised — so it is the smallest
        // reproducer. Fix is the same shape as `e_pq_mark` (extend
        // `RepairSearchNode` with `e_act_mark` and pass it to
        // `E.backtrack_to`); kept out of this change to bound scope. Until
        // then the curated list keeps the rotation safe. Multi-attempt
        // looping inside `run_attempt` still lets fast workers fill the slice
        // by cycling through the 8-config list, which the issue's #1
        // acceptance bullet (FPR CPU% on tbfp-network) cares about.
        let idx = (((self.worker_idx + self.attempt_idx) % kNumInitialFprConfigs)
            + kNumInitialFprConfigs)
            % kNumInitialFprConfigs;
        let cfg = kInitialFprConfigs[idx as usize];
        self.strat_idx = cfg.strat_idx;
        self.mode = cfg.mode;
    }

    /// mip-heuristics (36f7d9d): fpr.cpp::FprWorker::charge (lines 281-287)
    fn charge(&mut self, attempt: &AttemptResult) {
        if attempt.found_improvement {
            self.base.charge_improvement(attempt.effort);
        } else {
            self.base.charge_no_improvement(attempt.effort);
        }
    }

    /// mip-heuristics (36f7d9d): fpr.cpp::FprWorker::run_attempt (lines 294-488)
    ///
    /// Cognitive complexity 26 (threshold 25). Kept whole: the worker's multi-attempt loop:
    /// pause/resume of an in-flight DFS and rotation through kInitialFprConfigs. Decomposing it would
    /// move work across a worker's inner loop, and the closeout takes no unmeasured performance risk;
    /// the standards also rank fidelity to the reference algorithm above mechanical extraction.
    pub(crate) fn run_attempt(
        &mut self,
        attempt_budget: usize,
        sink: &mut S5Sink<'_>,
    ) -> AttemptResult {
        let mut attempt = AttemptResult::default();
        if self.base.finished {
            return attempt;
        }
        // The one clock FPR's search polls below this worker (issue #117).
        let deadline = self.exec.deadline;

        // Issue #77 lifecycle. Two mechanics in play:
        //
        // (1) Pause/resume *across* calls. When the DFS exhausts the per-call
        //     slice, `fpr_attempt_step` returns `kBudgetGate` and we return so
        //     peers can run their next slice. The DFS state (var_order_cursor,
        //     nodes_visited, found_complete, dfs_stack, prop_engine) lives in
        //     `attempt_state_` + `scratch_` and is preserved until the next call
        //     resumes the same attempt. Without this, a long DFS subtree gets
        //     truncated and discarded each attempt — the parallelism bottleneck
        //     this issue exists to fix.
        //
        // (2) Multi-attempt fill *within* a call. When an attempt verdicts
        //     before exhausting the slice (e.g. a fast-failing strategy on a
        //     hard instance, or a leaf found in 5% of the budget), we start
        //     the next attempt with the next rotation slot rather than idling
        //     at the runner's barrier waiting for slow peers. This is what
        //     drives FPR's CPU% up on `tbfp-network` (acceptance bullet #1) —
        //     pause/resume alone fixes truncation but leaves fast workers
        //     idle while a slow worker holds the barrier.
        //
        // The safety cap (`kMaxAttemptsPerCall`) and no-progress guard guard
        // against degenerate models where every attempt verdicts with near-zero
        // recorded effort (`infeasible-mip0` initial-propagation short-circuit):
        // without bounds we'd burn the slice on `fpr_attempt_begin`'s
        // O(ncol+nrow) setup churn alone.

        // No pool restart is pulled here (issue #122): `fpr_attempt_begin` has
        // no seed input left to feed it — the propagation engine's own
        // reset() baseline is never read before Phase 2/2.5 overwrite it — so
        // `sink_.get_restart` bought nothing but a pool-mutex acquisition, an
        // `ncol`-sized copy, and RNG draws from `rng_` on every call (a roll,
        // two parent indices, and up to `ncol` crossover coin flips whenever
        // the pool is non-empty — see `SolutionPool::get_restart`). Removing
        // it moves this worker's RNG stream from the first attempt onward, on
        // top of the shift from deleting the seeding block itself.

        // 32 attempts × 1 mutex op (`sink_.offer`, only on a feasible verdict) ×
        // N workers is a theoretical upper bound on pool-mutex acquisitions per
        // outer attempt (down from 2 before issue #122 removed the per-attempt
        // `get_restart` pull above). In practice the cap is rarely approached:
        // each attempt's begin charges O(nnz) coefficient accesses (initial
        // propagate), and the per-call slice is sized so an instance with
        // non-trivial DFS spends most of the slice inside step rather than
        // restarting attempts. HighsSpinMutex critical sections in `try_add`
        // are sub-microsecond (lower_bound over <= kPoolCapacity entries plus
        // an O(ncol) Hamming or single solution copy). Even worst-case, total
        // mutex-time per attempt is bounded ms (review R2 U-1 / Finding 3).
        const kMaxAttemptsPerCall: i32 = 32;
        let mut attempts_started: i32 = 0;
        let mut prev_loop_effort: usize = 0;

        // Issue #111: this call may not outrun what is left of the worker's
        // patience allowance. Without the clamp the *slice* bounds the call,
        // and the slice is `HeuristicBudget::attempt_cap` = `total / (10N)`,
        // which grows with the effort option — so a stalled worker would
        // still overshoot its absolute ceiling by a budget-proportional
        // amount and the gate would only half-bind. Clamping here also makes
        // the DFS pause at the ceiling rather than past it, since
        // `fpr_attempt_step` is handed `budget_remaining` derived from
        // `call_cap`.
        let patience_room = self
            .base
            .stale_budget
            .saturating_sub(self.base.effort_since_improvement);
        let call_cap = cpp_min(attempt_budget, cpp_max(patience_room, 1usize));

        while attempt.effort < call_cap {
            // `past_deadline()`, not `terminated()` (issue #114). This poll
            // exists so a 32-attempt fill cannot outrun a solver timeout, and
            // the deadline is the half it needs; `terminatorTerminated()`
            // *writes* `mipsolver.termination_status_`, so calling it from
            // this worker thread was a data race whenever a terminator was
            // attached — the one documented exception to the poller seat, now
            // retired. The terminator is still observed by the runner between
            // attempts, one attempt later than before.
            if self.exec.past_deadline() {
                break;
            }
            if attempts_started > 0 && attempt.effort == prev_loop_effort {
                // Defensive belt-and-braces guard. Today this branch is
                // unreachable: degenerate `ncol==0||nrow==0` models are filtered
                // out by `fpr::run` before workers are constructed,
                // every begin runs at least one `E.propagate(-1)` round (>0 ops
                // on any non-empty model), and finish always adds
                // `c.ar_index.size() > 0` for the LHS sum. Keep the guard so a
                // future change that relaxes any of the above (e.g., a Phase 1
                // shortcut that skips initial propagate) cannot silently turn
                // this loop into an infinite attempt-cycler.
                break;
            }
            prev_loop_effort = attempt.effort;

            // Advance the per-worker rotation BEFORE building cfg so that
            // `cfg.strategy` / `cfg.mode` / `cfg.precomputed_var_order` reflect
            // the current attempt's choice. Earlier draft built cfg from the
            // previous attempt's strat/mode and re-assigned 4 fields after the
            // rotation advance — a maintenance hazard if a future cfg field
            // is added (review R2 CF-2).
            let alive = self.attempt_state.phase != FprPhase::kIdle;
            if !alive {
                if attempts_started >= kMaxAttemptsPerCall {
                    break;
                }
                // #ifndef NDEBUG (fpr.cpp:403-407): g_multi_attempt_iters counter dropped.
                attempts_started += 1;
                self.select_config_for_current_attempt();
            }

            let strat = kFprStrategies[self.strat_idx as usize];
            let var_order = &self.var_orders[self.strat_idx as usize];
            let mut cfg = FprConfig::default();
            // Unread by the lifecycle API (issue #156): the DFS gate is the
            // per-call `effort_remaining` slice below, and Phase 3 takes no
            // effort cap from its caller at all.
            cfg.max_effort = usize::MAX;
            cfg.cont_fallback = None;
            cfg.mode = self.mode;
            cfg.strategy = Some(&strat);
            cfg.lp_ref = None;
            cfg.precomputed_var_order = Some(var_order.as_slice());
            cfg.precomputed_var_order_size = var_order.len() as HighsInt;
            cfg.scratch = Some(&mut self.scratch);

            let alive = self.attempt_state.phase != FprPhase::kIdle;
            if !alive {
                fpr_attempt_begin(
                    &mut self.attempt_state,
                    self.problem,
                    &deadline,
                    &mut cfg,
                    &mut self.rng,
                    self.attempt_idx,
                );
                // `attempt_state_.phase` is now `kDfs` (or `kReadyToFinish`
                // if Phase 1 already produced a complete fixing); either way
                // `attempt_alive()` is true on the next iteration.
                attempt.effort += self.attempt_state.effort_consumed;
            }

            if self.attempt_state.phase == FprPhase::kDfs {
                let before_step = self.attempt_state.effort_consumed;
                let budget_remaining = call_cap.saturating_sub(attempt.effort);
                let outcome = fpr_attempt_step(
                    &mut self.attempt_state,
                    self.problem,
                    &deadline,
                    &mut cfg,
                    &mut self.rng,
                    budget_remaining,
                );
                attempt.effort += self.attempt_state.effort_consumed - before_step;
                if outcome == FprStepResult::kBudgetGate {
                    // #ifndef NDEBUG (fpr.cpp:443-445): g_budget_gate_hits counter dropped.
                    // Attempt paused at the per-call slice boundary — return so
                    // peers do their next attempt's work and we resume here
                    // next call. Charge what the pause spent: the patience
                    // counter is in effort units, so slicing an attempt
                    // across calls costs it exactly what running it in one
                    // call would (see the class comment).
                    self.charge(&attempt);
                    return attempt;
                }
                // kVerdictReady — DFS ended (leaf found or stack/node-limit
                // exhausted), proceed to finish.
            }

            let before_finish = self.attempt_state.effort_consumed;
            let result = fpr_attempt_finish(
                &mut self.attempt_state,
                self.problem,
                &deadline,
                &mut cfg,
                &mut self.rng,
            );
            attempt.effort += self.attempt_state.effort_consumed - before_finish;

            // Incumbent improvement, not "this attempt reached a feasible
            // point" (#111) and not mere pool acceptance (#116): FPR reaches
            // the same neighbourhood over and over on some models, and it is
            // the heuristic the two signals differ most for — ~3.3 M
            // acceptances against 590 incumbent improvements over #113's 233
            // instances, five orders of magnitude. `offer` is still called
            // for every feasible result, so nothing that reached HiGHS before
            // stops reaching it — only the verdict read has changed.
            // `effort_at`: the worker's cumulative charge including this
            // attempt's work so far. `charge(attempt)` runs at both exits from
            // this function, so `base_.total_effort` does not yet contain it.
            if result.found_feasible {
                let effort_at = self.trace.at(self.base.total_effort + attempt.effort);
                let offer_result =
                    sink.offer(result.objective, &result.solution, &self.trace, effort_at);
                if offer_result.improved_incumbent {
                    attempt.found_improvement = true;
                }
            }

            self.attempt_idx += 1;
            // `fpr_attempt_finish` set `attempt_state_.phase = kIdle`, so
            // `attempt_alive()` is false on the next iteration.
        }

        self.charge(&attempt);
        attempt
    }
}

// ---------------------------------------------------------------------------
// Parallel FPR
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr.cpp::run (lines 494-546)
pub(crate) fn run(
    problem: &S5Problem,
    budget: &HeuristicBudget,
    exec: &ExecutionContext,
    sink: &mut S5Sink<'_>,
) -> DispatchOutcome {
    // Inlined: problem.degenerate() is ncol == 0 || nrow == 0; budget.disabled() is total == 0.
    if problem.ncol == 0 || problem.nrow == 0 || budget.total == 0 {
        return DispatchOutcome::default();
    }

    // Precompute var_orders sequentially before any parallel region. A
    // deadline that passes inside it retires the whole dispatch: the table
    // is incomplete, no worker exists yet, and the runner's first act
    // would be to stop anyway (issue #117).
    //
    // Reported as an abandoned setup rather than as a plain zero (issue
    // #119): this dispatch did not search and produce nothing, it never
    // searched, and the two are the same `effort=0 found=0` line to every
    // consumer that cannot tell them apart.
    let mut var_orders: Vec<Vec<HighsInt>> = Vec::new();
    {
        // FPR reads the clique table and domain through sink.mipsolver().
        // Must not hold that borrow across an offer.
        let ms = sink
            .mipsolver()
            .expect("fpr::run needs the IncumbentSink mipsolver (D4)");
        if !precompute_var_orders(ms, &exec.deadline, &mut var_orders) {
            // Inlined: DispatchOutcome::abandoned()
            return DispatchOutcome {
                effort: 0,
                abandoned_setup: true,
            };
        }
    }

    let mut workers: Vec<FprWorker<'_>> = Vec::with_capacity(exec.num_workers);
    for w in 0..exec.num_workers {
        let seed = exec.worker_seed(w as i32);
        // `budget.worker_stale` is this worker's share of the dispatch's
        // absolute patience ceiling. There is no attempt-wide effort hint
        // beside it any more: the `budget.total >> 2` that used to sit here
        // sized Phase 3's repair/WalkSAT sub-budgets, and issue #156
        // removed those — an attempt spans calls, so `PropEngine::effort()`
        // outgrows any such number and the cap it produced arrived as 0.
        workers.push(FprWorker::new(
            problem,
            exec,
            &var_orders,
            w as i32,
            seed,
            budget.worker_stale,
        ));
    }

    struct FprOppState {
        worker_idx: i32,
    }

    let effort = run_opportunistic_loop(
        exec,
        sink,
        budget,
        |worker_idx, _rng, _sink| FprOppState { worker_idx },
        |state: &mut FprOppState, _rng: &mut S5Rng, run_cap: usize, sink: &mut S5Sink<'_>| {
            let worker = &mut workers[state.worker_idx as usize];
            // A retired worker reports zero effort, which retires its
            // slot in `run_opportunistic_loop` (issue #111 gave FPR the
            // worker-level gate it had been doing without). Deliberately
            // no `attempt_with_rebuild`: FPR's diversity already comes
            // from the per-attempt rotation through `kInitialFprConfigs`,
            // so a rebuilt worker would resume the same rotation with a
            // fresh patience allowance and the gate would bound nothing.
            worker.run_attempt(run_cap, sink)
        },
    );

    DispatchOutcome {
        effort,
        abandoned_setup: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;

    /// mip-heuristics (36f7d9d): test_fpr.cpp::rotation matches initial configs
    #[test]
    fn s5_fpr_rotation_matches_initial_configs() {
        // Build a minimal problem and worker to test rotation logic.
        let problem = S5Problem {
            ncol: 2,
            nrow: 2,
            nnz: 4,
            col_cost: vec![1.0, 1.0],
            incumbent: vec![0.0, 0.0],
            col_lower: vec![0.0, 0.0],
            col_upper: vec![1.0, 1.0],
            row_lower: vec![0.0, 0.0],
            row_upper: vec![1.0, 1.0],
            offset: 0.0,
            integrality: vec![HighsVarType::kInteger; 2],
            ar_start: vec![0, 2, 4],
            ar_index: vec![0, 1, 0, 1],
            ar_value: vec![1.0, 1.0, 1.0, 1.0],
            csc: CscMatrix {
                col_start: vec![0, 2, 4],
                col_row: vec![0, 1, 0, 1],
                col_val: vec![1.0, 1.0, 1.0, 1.0],
            },
            uplocks: vec![0; 2],
            downlocks: vec![0; 2],
            feastol: 1e-6,
            epsilon: 1e-6,
            binary: vec![1u8; 2],
        };
        let exec = ExecutionContext {
            num_workers: 1,
            base_seed: 42,
            deadline: Deadline::default(),
            gub_swap: false,
        };
        let var_orders: Vec<Vec<HighsInt>> = vec![vec![0, 1]; kNumFprStrategies as usize];
        let mut worker = FprWorker::new(&problem, &exec, &var_orders, 0, 42, 1000);

        for attempt_idx in 0..8 {
            worker.attempt_idx = attempt_idx;
            worker.select_config_for_current_attempt();
            let expected = kInitialFprConfigs[(attempt_idx % kNumInitialFprConfigs) as usize];
            assert_eq!(worker.strat_idx, expected.strat_idx);
            assert_eq!(worker.mode, expected.mode);
        }
    }

    /// mip-heuristics (36f7d9d): test_fpr.cpp::disabled and degenerate return default
    #[test]
    fn s5_fpr_run_disabled_and_degenerate_return_default() {
        let problem = S5Problem {
            ncol: 2,
            nrow: 2,
            nnz: 4,
            col_cost: vec![1.0, 1.0],
            incumbent: vec![0.0, 0.0],
            col_lower: vec![0.0, 0.0],
            col_upper: vec![1.0, 1.0],
            row_lower: vec![0.0, 0.0],
            row_upper: vec![1.0, 1.0],
            offset: 0.0,
            integrality: vec![HighsVarType::kInteger; 2],
            ar_start: vec![0, 2, 4],
            ar_index: vec![0, 1, 0, 1],
            ar_value: vec![1.0, 1.0, 1.0, 1.0],
            csc: CscMatrix {
                col_start: vec![0, 2, 4],
                col_row: vec![0, 1, 0, 1],
                col_val: vec![1.0, 1.0, 1.0, 1.0],
            },
            uplocks: vec![0; 2],
            downlocks: vec![0; 2],
            feastol: 1e-6,
            epsilon: 1e-6,
            binary: vec![1u8; 2],
        };
        let exec = ExecutionContext {
            num_workers: 1,
            base_seed: 42,
            deadline: Deadline::default(),
            gub_swap: false,
        };
        let mut sink = S5Sink::new(&[HighsVarType::kInteger; 2], 0);

        // Disabled budget
        let budget_disabled = HeuristicBudget {
            total: 0,
            ..Default::default()
        };
        let outcome = run(&problem, &budget_disabled, &exec, &mut sink);
        assert_eq!(outcome.effort, 0);
        assert!(!outcome.abandoned_setup);

        // Degenerate problem (ncol == 0)
        let degenerate_problem = S5Problem {
            ncol: 0,
            nrow: 2,
            nnz: 0,
            col_cost: vec![],
            incumbent: vec![],
            col_lower: vec![],
            col_upper: vec![],
            row_lower: vec![0.0, 0.0],
            row_upper: vec![1.0, 1.0],
            offset: 0.0,
            integrality: vec![],
            ar_start: vec![0, 0],
            ar_index: vec![],
            ar_value: vec![],
            csc: CscMatrix {
                col_start: vec![0],
                col_row: vec![],
                col_val: vec![],
            },
            uplocks: vec![],
            downlocks: vec![],
            feastol: 1e-6,
            epsilon: 1e-6,
            binary: vec![],
        };
        let budget = HeuristicBudget {
            total: 1000,
            per_worker: 1000,
            attempt_cap: 100,
            stale: 500,
            worker_stale: 500,
            fj_stall_ends_run: false,
            restart: RestartPoint::kDefault,
            randomize_restart_weights: false,
        };
        let outcome = run(&degenerate_problem, &budget, &exec, &mut sink);
        assert_eq!(outcome.effort, 0);
        assert!(!outcome.abandoned_setup);
    }

    /// mip-heuristics (36f7d9d): test_fpr.cpp::worker run_attempt base sink records offers
    #[test]
    fn s5_fpr_worker_run_attempt_base_sink_records_offers() {
        // 2x2 binary problem, as in fpr_core.rs tests make_problem_2x2
        let problem = S5Problem {
            ncol: 2,
            nrow: 2,
            nnz: 4,
            col_cost: vec![1.0, 1.0],
            incumbent: vec![0.0, 0.0],
            col_lower: vec![0.0, 0.0],
            col_upper: vec![1.0, 1.0],
            row_lower: vec![0.0, 0.0],
            row_upper: vec![1.0, 1.0],
            offset: 0.0,
            integrality: vec![HighsVarType::kInteger; 2],
            ar_start: vec![0, 2, 4],
            ar_index: vec![0, 1, 0, 1],
            ar_value: vec![1.0, 1.0, 1.0, 1.0],
            csc: CscMatrix {
                col_start: vec![0, 2, 4],
                col_row: vec![0, 1, 0, 1],
                col_val: vec![1.0, 1.0, 1.0, 1.0],
            },
            uplocks: vec![0; 2],
            downlocks: vec![0; 2],
            feastol: 1e-6,
            epsilon: 1e-6,
            binary: vec![1u8; 2],
        };
        let exec = ExecutionContext {
            num_workers: 1,
            base_seed: 42,
            deadline: Deadline::default(),
            gub_swap: false,
        };
        let var_orders: Vec<Vec<HighsInt>> = vec![vec![0, 1]; kNumFprStrategies as usize];

        // Same seed, same slice: two workers on fresh base sinks give the
        // same effort, the same offers accepted and the same best solution.
        let run_once = || {
            let mut sink = S5Sink::new(&[HighsVarType::kInteger; 2], 0);
            let mut worker = FprWorker::new(&problem, &exec, &var_orders, 0, 42, 100000);
            let result = worker.run_attempt(100000, &mut sink);
            let mut best = Vec::new();
            let has_best = sink.copy_best(&mut best);
            (
                result.effort,
                result.found_improvement,
                worker.attempt_idx,
                sink.accepted(),
                has_best,
                best,
            )
        };
        let first = run_once();
        let second = run_once();
        assert!(first.0 > 0);
        // x = (0, 0) is feasible for 0 <= x0 + x1 <= 1, so FPR records an offer.
        assert!(first.3 > 0);
        assert!(first.4);
        assert_eq!(first, second);
    }
}
