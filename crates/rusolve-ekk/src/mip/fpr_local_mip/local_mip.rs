//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! local_mip.cpp (is_solution_feasible, resolve_worker_start, Slots, run).

use super::context::{rebuild_restart_point, run_opportunistic_loop, take_restart};
use super::lm_construction::{construct_initial_solution, construction_effort_cap};
use super::lm_worker::{LocalMipWorker, perturb_solution};
use super::rng::S5Rng;
use super::sink::S5Sink;
use super::*;
use crate::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

// kInstrumented block (local_mip.cpp:30-118, incl. g_* counters, bump_counter,
// warm_start_counters, deadline poll counters and note_* functions) is dropped.

/// mip-heuristics (36f7d9d): local_mip.cpp::is_solution_feasible (lines 123-159)
pub(crate) fn is_solution_feasible(problem: &S5Problem, solution: &[f64]) -> bool {
    let integrality = &problem.integrality;
    let ncol = problem.ncol;
    let nrow = problem.nrow;
    let feastol = problem.feastol;
    let inttol = problem.epsilon;
    let ar_start = &problem.ar_start;
    let ar_index = &problem.ar_index;
    let ar_value = &problem.ar_value;
    if solution.len() != ncol as usize {
        return false;
    }
    // Integer feasibility.
    for j in 0..ncol {
        if integrality[j as usize] == HighsVarType::kInteger
            || integrality[j as usize] == HighsVarType::kImplicitInteger
        {
            if (solution[j as usize] - solution[j as usize].round()).abs() > inttol {
                return false;
            }
        }
        if solution[j as usize] < problem.col_lower[j as usize] - feastol
            || solution[j as usize] > problem.col_upper[j as usize] + feastol
        {
            return false;
        }
    }
    // Row feasibility — walk the row-wise matrix once.
    for i in 0..nrow {
        let mut lhs = 0.0;
        let mut k = ar_start[i as usize];
        while k < ar_start[(i + 1) as usize] {
            lhs += ar_value[k as usize] * solution[ar_index[k as usize] as usize];
            k += 1;
        }
        if lhs < problem.row_lower[i as usize] - feastol
            || lhs > problem.row_upper[i as usize] + feastol
        {
            return false;
        }
    }
    true
}

/// mip-heuristics (36f7d9d): local_mip.cpp::resolve_worker_start (lines 198-261)
pub(crate) fn resolve_worker_start(
    problem: &S5Problem,
    sink: &mut S5Sink<'_>,
    trace: &WorkerTrace,
    max_effort: usize,
    seed: u32,
    mut cold_start_cache: Option<&mut Vec<f64>>,
    effort_out: Option<&mut usize>,
) -> Vec<f64> {
    let mut start = Vec::new();
    if sink.copy_best(&mut start) {
        return start;
    }
    if !problem.incumbent.is_empty() {
        return problem.incumbent.clone();
    }
    // Cold start: neither the pool nor the incumbent has a solution.
    // Re-use a cached construction if one was produced earlier in this
    // dispatch; otherwise run the paper's construction phase and cache
    // it for subsequent workers.
    if let Some(cache) = cold_start_cache.as_deref_mut() {
        if !cache.is_empty() {
            return cache.clone();
        }
    }
    let mut rng = S5Rng::new(seed as u64);
    let mut constructed = Vec::new();
    let construction_effort = construct_initial_solution(
        problem,
        &mut rng,
        construction_effort_cap(max_effort),
        &mut constructed,
    );
    if let Some(eo) = effort_out {
        *eo += construction_effort;
    }
    // If the construction happens to land on a feasible integer point,
    // publish it to the shared pool so downstream heuristics (and
    // HiGHS's own incumbent path) pick it up.
    if !constructed.is_empty() && is_solution_feasible(problem, &constructed) {
        // ProblemView::objective (heuristic_context.h:231-237)
        let mut obj = problem.offset;
        for j in 0..problem.ncol {
            obj += problem.col_cost[j as usize] * constructed[j as usize];
        }
        let _ = sink.offer(obj, &constructed, trace, trace.at(construction_effort));
    }
    if let Some(cache) = cold_start_cache {
        *cache = constructed.clone();
    }
    constructed
}

/// mip-heuristics (36f7d9d): local_mip.cpp::retire_trace (lines 269-273)
fn retire_trace(worker: Option<&LocalMipWorker<'_, '_>>, trace: &mut WorkerTrace) {
    if let Some(w) = worker {
        trace.effort_base = w.traced_effort();
    }
}

/// mip-heuristics (36f7d9d): local_mip.cpp::LmState (lines 275-282)
struct LmState<'a, 'l> {
    worker: Option<LocalMipWorker<'a, 'l>>,
    trace: WorkerTrace,
}

/// mip-heuristics (36f7d9d): local_mip.cpp::Slots (lines 289-301)
struct Slots<'a> {
    problem: &'a S5Problem,
    budget: &'a HeuristicBudget,
    exec: &'a ExecutionContext,
    construction_effort: AtomicUsize,
}

impl<'a> Slots<'a> {
    /// mip-heuristics (36f7d9d): local_mip.cpp::Slots::occupy (lines 306-321)
    fn occupy<'l>(
        &self,
        worker_idx: i32,
        seed: u32,
        mut start: Vec<f64>,
        construction: usize,
        rng: &mut S5Rng,
    ) -> LmState<'a, 'l> {
        self.construction_effort
            .fetch_add(construction, Ordering::Relaxed);
        if worker_idx != 0 {
            perturb_solution(
                &mut start,
                &self.problem.binary,
                &self.problem.integrality,
                &self.problem.col_lower,
                &self.problem.col_upper,
                self.problem.ncol,
                rng,
            );
        }
        let trace = WorkerTrace {
            worker: worker_idx,
            effort_base: construction,
        };
        let worker = LocalMipWorker::new(
            self.problem,
            *self.exec,
            self.budget.per_worker,
            self.budget.worker_stale,
            seed,
            Some(&start),
            trace,
            None,
        );
        LmState {
            worker: Some(worker),
            trace,
        }
    }

    /// mip-heuristics (36f7d9d): local_mip.cpp::Slots::rebuild (lines 327-347)
    fn rebuild<'l>(&self, state: &mut LmState<'a, 'l>, rng: &mut S5Rng, sink: &mut S5Sink<'_>) {
        retire_trace(state.worker.as_ref(), &mut state.trace);
        let mut restart_sol = Vec::new();
        if !rebuild_restart_point(
            self.budget.restart,
            None,
            sink,
            rng,
            &mut restart_sol,
            self.problem,
        ) {
            // The C++ two-branch default_restart call with an empty source is
            // one call with None (both branches are identical when the source
            // is empty).
            self.default_restart(state, rng, &mut restart_sol, sink);
        }
        let seed = rng.next_u64() as u32;
        let worker = LocalMipWorker::new(
            self.problem,
            *self.exec,
            self.budget.per_worker,
            self.budget.worker_stale,
            seed,
            Some(&restart_sol),
            state.trace,
            None,
        );
        state.worker = Some(worker);
        if self.budget.randomize_restart_weights {
            if let Some(ref mut w) = state.worker {
                w.randomize_weights(rng);
            }
        }
    }

    /// mip-heuristics (36f7d9d): local_mip.cpp::Slots::default_restart (lines 351-383)
    fn default_restart<'l>(
        &self,
        state: &mut LmState<'a, 'l>,
        rng: &mut S5Rng,
        restart_sol: &mut Vec<f64>,
        sink: &mut S5Sink<'_>,
    ) {
        if !take_restart(None, rng, restart_sol, self.problem)
            && !sink.get_restart(rng, restart_sol, RestartPolicy::kMixed)
        {
            if !self.problem.incumbent.is_empty() {
                *restart_sol = self.problem.incumbent.clone();
            } else {
                let cseed = rng.next_u64() as u32;
                let mut construct_rng = S5Rng::new(cseed as u64);
                let my_construction_effort = construct_initial_solution(
                    self.problem,
                    &mut construct_rng,
                    construction_effort_cap(self.budget.per_worker),
                    restart_sol,
                );
                self.construction_effort
                    .fetch_add(my_construction_effort, Ordering::Relaxed);
                state.trace.effort_base += my_construction_effort;
            }
        }
        perturb_solution(
            restart_sol,
            &self.problem.binary,
            &self.problem.integrality,
            &self.problem.col_lower,
            &self.problem.col_upper,
            self.problem.ncol,
            rng,
        );
    }
}

/// mip-heuristics (36f7d9d): local_mip.cpp::run (lines 388-481)
pub(crate) fn run(
    problem: &S5Problem,
    budget: &HeuristicBudget,
    exec: &ExecutionContext,
    sink: &mut S5Sink<'_>,
) -> DispatchOutcome {
    // ProblemView::degenerate (heuristic_context.h:228)
    // HeuristicBudget::disabled (heuristic_context.h:392)
    if problem.ncol == 0 || problem.nrow == 0 || budget.total == 0 {
        return DispatchOutcome::default();
    }

    let slots = Slots {
        problem,
        budget,
        exec,
        construction_effort: AtomicUsize::new(0),
    };

    // Cold-start cache shared across all workers of this dispatch.
    // Sequential path (num_workers = 1): the C++ mutex is dropped, one thread.
    let mut cold_start_cache: Vec<f64> = Vec::new();

    // Prime the cache on this thread, before any worker starts.
    if budget.total > 0 {
        let mut primed_effort = 0usize;
        resolve_worker_start(
            problem,
            sink,
            &WorkerTrace {
                worker: -1,
                effort_base: 0,
            },
            budget.per_worker,
            exec.base_seed,
            Some(&mut cold_start_cache),
            Some(&mut primed_effort),
        );
        slots
            .construction_effort
            .fetch_add(primed_effort, Ordering::Relaxed);
    }

    let total_effort = run_opportunistic_loop(
        exec,
        sink,
        budget,
        |worker_idx: i32, rng: &mut S5Rng, sink: &mut S5Sink<'_>| -> LmState<'_, '_> {
            let seed = rng.next_u64() as u32;
            let mut local_cache = cold_start_cache.clone();
            let mut my_construction_effort = 0usize;
            let start = resolve_worker_start(
                problem,
                sink,
                &WorkerTrace {
                    worker: worker_idx,
                    effort_base: 0,
                },
                budget.per_worker,
                seed,
                Some(&mut local_cache),
                Some(&mut my_construction_effort),
            );
            if !local_cache.is_empty() && cold_start_cache.is_empty() {
                cold_start_cache = local_cache;
            }
            slots.occupy(worker_idx, seed, start, my_construction_effort, rng)
        },
        |state: &mut LmState<'_, '_>,
         rng: &mut S5Rng,
         run_cap: usize,
         sink: &mut S5Sink<'_>|
         -> AttemptResult {
            if state.worker.as_ref().is_none_or(|w| w.finished()) {
                slots.rebuild(state, rng, sink);
            }
            match state.worker.as_mut() {
                Some(w) => w.run_attempt(run_cap, sink),
                None => AttemptResult::default(),
            }
        },
    );

    DispatchOutcome {
        effort: total_effort + slots.construction_effort.load(Ordering::Relaxed),
        abandoned_setup: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_problem_2x2_binary() -> S5Problem {
        let mut p = S5Problem::default();
        p.ncol = 2;
        p.nrow = 2;
        p.nnz = 4;
        p.col_lower = vec![0.0, 0.0];
        p.col_upper = vec![1.0, 1.0];
        p.row_lower = vec![0.0, 1.0];
        p.row_upper = vec![2.0, 2.0];
        p.col_cost = vec![1.0, 1.0];
        p.offset = 0.0;
        p.integrality = vec![HighsVarType::kInteger, HighsVarType::kInteger];
        p.ar_start = vec![0, 2, 4];
        p.ar_index = vec![0, 1, 0, 1];
        p.ar_value = vec![1.0, 1.0, 1.0, 1.0];
        p.csc = CscMatrix {
            col_start: vec![0, 2, 4],
            col_row: vec![0, 1, 0, 1],
            col_val: vec![1.0, 1.0, 1.0, 1.0],
        };
        p.uplocks = vec![0, 0];
        p.downlocks = vec![0, 0];
        p.feastol = 1e-7;
        p.epsilon = 1e-6;
        p.binary = vec![1, 1];
        p
    }

    #[test]
    fn s5_lm_degenerate_problem_returns_default() {
        let problem = S5Problem::default();
        let budget = HeuristicBudget {
            total: 100000,
            per_worker: 100000,
            attempt_cap: 100000,
            stale: 100000,
            worker_stale: 100000,
            fj_stall_ends_run: false,
            restart: RestartPoint::kDefault,
            randomize_restart_weights: false,
        };
        let exec = ExecutionContext::default();
        let mut sink = S5Sink::new(&problem.integrality, 0);
        let outcome = run(&problem, &budget, &exec, &mut sink);
        assert_eq!(outcome.effort, 0);
        assert!(!outcome.abandoned_setup);
        assert_eq!(sink.accepted(), 0);
    }

    #[test]
    fn s5_lm_disabled_budget_returns_default() {
        let problem = make_problem_2x2_binary();
        let budget = HeuristicBudget::default();
        let exec = ExecutionContext::default();
        let mut sink = S5Sink::new(&problem.integrality, 0);
        let outcome = run(&problem, &budget, &exec, &mut sink);
        assert_eq!(outcome.effort, 0);
        assert!(!outcome.abandoned_setup);
    }

    #[test]
    fn s5_lm_run_is_deterministic() {
        let problem = make_problem_2x2_binary();
        let budget = HeuristicBudget {
            total: 100000,
            per_worker: 100000,
            attempt_cap: 100000,
            stale: 100000,
            worker_stale: 100000,
            fj_stall_ends_run: false,
            restart: RestartPoint::kDefault,
            randomize_restart_weights: false,
        };
        let exec = ExecutionContext::default();

        let mut sink1 = S5Sink::new(&problem.integrality, 0);
        let outcome1 = run(&problem, &budget, &exec, &mut sink1);
        let accepted1 = sink1.accepted();
        let mut best1 = Vec::new();
        sink1.copy_best(&mut best1);

        let mut sink2 = S5Sink::new(&problem.integrality, 0);
        let outcome2 = run(&problem, &budget, &exec, &mut sink2);
        let accepted2 = sink2.accepted();
        let mut best2 = Vec::new();
        sink2.copy_best(&mut best2);

        assert_eq!(outcome1.effort, outcome2.effort);
        assert_eq!(outcome1.abandoned_setup, outcome2.abandoned_setup);
        assert_eq!(accepted1, accepted2);
        assert_eq!(best1, best2);
        // Non-vacuous: the run did search and found a feasible point.
        assert!(outcome1.effort > 0);
        assert!(accepted1 >= 1);
        assert!(is_solution_feasible(&problem, &best1));
    }

    #[test]
    fn s5_lm_resolve_worker_start_prefers_incumbent_snapshot() {
        // test_local_mip.cpp:513-557 (#98): with an empty pool and a
        // non-empty snapshot, the start is the incumbent and costs nothing.
        let mut problem = make_problem_2x2_binary();
        problem.incumbent = vec![1.0, 0.0];
        let mut sink = S5Sink::new(&problem.integrality, 0);
        let mut cache = Vec::new();
        let mut effort = 0usize;
        let start = resolve_worker_start(
            &problem,
            &mut sink,
            &WorkerTrace::default(),
            100000,
            0,
            Some(&mut cache),
            Some(&mut effort),
        );
        assert_eq!(start, vec![1.0, 0.0]);
        assert_eq!(effort, 0);
        assert!(cache.is_empty());
    }

    #[test]
    fn s5_lm_resolve_worker_start_cold_start_constructs_and_caches() {
        // test_local_mip.cpp:259-305 (#75): pool and incumbent empty, so the
        // construction runs once, is cached, and a second call hits the cache.
        let problem = make_problem_2x2_binary();
        let mut sink = S5Sink::new(&problem.integrality, 0);
        let mut cache = Vec::new();
        let mut effort = 0usize;
        let start = resolve_worker_start(
            &problem,
            &mut sink,
            &WorkerTrace::default(),
            100000,
            0,
            Some(&mut cache),
            Some(&mut effort),
        );
        assert_eq!(start.len(), 2);
        assert_eq!(cache, start);
        assert!(effort > 0);
        // A feasible construction is published to the pool; then the pool
        // branch wins over the cache.
        if is_solution_feasible(&problem, &start) {
            assert_eq!(sink.accepted(), 1);
        }
        let mut empty_sink = S5Sink::new(&problem.integrality, 0);
        let mut effort2 = 0usize;
        let again = resolve_worker_start(
            &problem,
            &mut empty_sink,
            &WorkerTrace::default(),
            100000,
            7,
            Some(&mut cache),
            Some(&mut effort2),
        );
        assert_eq!(again, start);
        assert_eq!(effort2, 0);
    }

    #[test]
    fn s5_lm_feasibility_checks_rows_and_integrality() {
        let problem = make_problem_2x2_binary();
        // [1,1]: row 1 lhs=2 <= 2 ok; row 2 lhs=2 >= 1 ok; integer ok
        assert!(is_solution_feasible(&problem, &[1.0, 1.0]));
        // [0.5, 0]: integrality violation
        assert!(!is_solution_feasible(&problem, &[0.5, 0.0]));
        // [0, 0]: row 2 lower=1, lhs=0 < 1 - feastol
        assert!(!is_solution_feasible(&problem, &[0.0, 0.0]));
    }
}
