//! Per-phase statistics of one MIP solve.
//!
//! Reporting only: nothing here is read by the search, so the trace cannot move
//! (`tools/bench/identity.sh`). HiGHS has most of these numbers but keeps them in locals or prints
//! them in the log. The HiGHS source for each:
//!
//! * presolve sizes and reductions: `HighsMipSolver::runPresolve` / `HighsPostsolveStack`
//!   (reductions are the postsolve stack's entries, one per reduction that needs a postsolve step;
//!   bound tightenings and coefficient changes leave no entry) and the log line
//!   `reportPresolveReductions`;
//! * root LP iterations: `HighsMipSolverData::firstrootlpiters`, `firstlpsolobj`
//!   (HighsMipSolverData.cpp, `evaluateRootNode`);
//! * separator calls and cuts found: `HighsSeparator::numCalls` / `numCutsFound`
//!   (HighsSeparator.cpp:32-45, `run`; HiGHS reads them only for the clock report);
//! * LP iterations by origin: `HighsMipSolverData::{total,heuristic,sepa,sb}_lp_iterations`, the
//!   numbers `solvingReport` prints (HighsMipSolver.cpp:1355-1370);
//! * improving incumbents by `solution_source`: the source passed to
//!   `HighsMipSolverData::addIncumbent` (HighsMipSolverData.cpp:1499), which HiGHS only prints as a
//!   letter in the log line (`solutionSourceToString`); the port tallies it;
//! * times: `HighsMipSolver::timer_` read at phase boundaries (HiGHS's `MipTimer` clocks are
//!   stubbed in the port, `mip_timer.rs`).

use crate::prelude::*;

/// Separator names in the order of `HighsSeparatorBox`.
pub const SEPARATOR_NAMES: [&str; 8] = [
    "tableau",
    "path",
    "modk",
    "machine_sched",
    "knapsack_cover",
    "ti_precedence",
    "netdesign",
    "dff",
];

/// What `HighsMipSolver` records while it runs.
/// A statistics counter that tasks of a parallel round may bump through the aliased solver
/// reference: relaxed atomic add, so the concurrent read-modify-write is not a data race
/// (rusolve). Reporting only, nothing reads it during the search.
#[derive(Default, Debug)]
pub struct RelaxedCounter(std::sync::atomic::AtomicI64);

impl RelaxedCounter {
    pub fn add(&self, n: i64) {
        self.0.fetch_add(n, std::sync::atomic::Ordering::Relaxed);
    }
    pub fn get(&self) -> i64 {
        self.0.load(std::sync::atomic::Ordering::Relaxed)
    }
}

impl Clone for RelaxedCounter {
    fn clone(&self) -> Self {
        Self(std::sync::atomic::AtomicI64::new(self.get()))
    }
}

#[derive(Clone, Default, Debug)]
pub struct MipRunStats {
    /// The work counters when `run()` started.
    pub work_start: crate::work::WorkSnapshot,
    /// The counters of the solve scope `run()` worked in (shared with the parent for a sub-MIP).
    pub work_counters: Option<std::sync::Arc<crate::work::WorkCounters>>,
    pub presolve_end_s: f64,
    pub setup_end_s: f64,
    /// Seconds spent in `evaluate_root_node` (summed over restarts).
    pub root_node_s: f64,
    /// `timer_` seconds when the last `evaluate_root_node` returned (0 if the root was never reached).
    pub root_end_s: f64,
    /// Calls, cuts found and seconds per separator (index = `SEPARATOR_NAMES`).
    pub sepa_calls: [RelaxedCounter; 8],
    pub sepa_cuts: [RelaxedCounter; 8],
    /// Nanoseconds per separator.
    pub sepa_ns: [RelaxedCounter; 8],
    /// `HighsSeparation::separation_round` calls (root and tree).
    pub sepa_rounds: RelaxedCounter,
    /// Improving incumbents by `solution_source + 1` (source `-1` = none is slot 0).
    pub improving_by_source: [i64; 28],
    /// `timer_` seconds and work at the first incumbent.
    pub first_incumbent_s: Option<f64>,
    pub first_incumbent_work: u64,
    /// mip_alns (Hendel 2022): sub-MIP calls, improving ones, nodes and work, by arm
    /// (RENS, RINS, mutation, crossover). Reporting only.
    pub alns_calls: i64,
    pub alns_found: i64,
    pub alns_nodes: i64,
    pub alns_work: u64,
    pub alns_calls_by_arm: [i64; 5],
    pub alns_found_by_arm: [i64; 5],
    /// mip_shift_and_propagate: ran, found a solution, probings, cutoffs, rows still violated
    /// at the end, seconds. Reporting only.
    pub shiftprop_calls: i64,
    pub shiftprop_found: i64,
    pub shiftprop_probings: i64,
    pub shiftprop_cutoffs: i64,
    pub shiftprop_violations_left: i64,
    pub shiftprop_s: f64,
    /// mip_node_rounding: calls, solutions found by simple rounding / activity rounding / the
    /// shifting repair, improving ones, seconds. Reporting only.
    pub noderound_calls: i64,
    pub noderound_found: [i64; 3],
    pub noderound_improving: i64,
    pub noderound_polished: i64,
    pub noderound_s: f64,
    /// mip_proximity / mip_local_branching (Fischetti-Monaci 2014, Fischetti-Lodi 2003): sub-MIP calls, improving
    /// ones, nodes and work. Reporting only.
    pub prox_calls: i64,
    pub prox_found: i64,
    pub prox_nodes: i64,
    pub prox_work: u64,
    pub lb_calls: i64,
    pub lb_found: i64,
    pub lb_nodes: i64,
    pub lb_work: u64,
}

impl MipRunStats {
    /// Records an improving incumbent found by `solution_source` (HighsMipSolverData::addIncumbent).
    pub fn record_incumbent(&mut self, solution_source: i32, time_s: f64) {
        let slot = (solution_source + 1).clamp(0, 27) as usize;
        self.improving_by_source[slot] += 1;
        if self.first_incumbent_s.is_none() {
            self.first_incumbent_s = Some(time_s);
            self.first_incumbent_work = crate::work::WorkSnapshot::now()
                .since(&self.work_start)
                .work();
        }
    }
}

/// One separator's totals.
#[derive(Clone, Debug)]
pub struct SeparatorStats {
    pub name: &'static str,
    pub calls: i64,
    pub cuts_found: i64,
    pub time_s: f64,
}

/// Improving incumbents of one `solution_source`.
#[derive(Clone, Debug)]
pub struct SourceStats {
    pub source: String,
    pub improving: i64,
}

/// Everything C5 reports, read after `HighsMipSolver::run` (see `HighsMipSolver::collect_stats`).
#[derive(Clone, Debug, Default)]
pub struct MipStats {
    /// the work of this solve (see `work.rs`).
    pub work: u64,
    /// The parts of `work`.
    pub work_parts: crate::work::WorkSnapshot,
    pub total_s: f64,
    pub presolve_s: f64,
    pub setup_s: f64,
    pub root_node_s: f64,
    pub search_s: f64,
    pub presolve_status: String,
    pub num_restarts: i64,
    pub rows_before: i64,
    pub cols_before: i64,
    pub nnz_before: i64,
    pub rows_after: i64,
    pub cols_after: i64,
    pub nnz_after: i64,
    pub reductions: Vec<(String, i64)>,
    pub root_lp_iterations: i64,
    pub root_lp_objective: f64,
    pub sepa_rounds: i64,
    pub separators: Vec<SeparatorStats>,
    pub cutpool_cuts: i64,
    pub heuristics: Vec<SourceStats>,
    pub num_nodes: i64,
    pub num_leaves: i64,
    pub lp_iterations_total: i64,
    pub lp_iterations_strong_branching: i64,
    pub lp_iterations_separation: i64,
    pub lp_iterations_heuristics: i64,
    pub first_incumbent_s: Option<f64>,
    pub first_incumbent_work: Option<u64>,
    pub primal_dual_integral: f64,
    /// Reporting-only event tallies, `crate::work::Event::NAMES` order (objective propagation and
    /// reduced-cost fixing). Include the sub-MIPs of the solve.
    pub events: [u64; crate::work::Event::COUNT],
    /// LP simplex iterations and LP solves by source, `crate::work::LP_SRC_NAMES` order
    ///. Include the sub-MIPs of the solve.
    pub lp_iters_by_source: [u64; crate::work::LP_SRC_COUNT],
    pub lp_solves_by_source: [u64; crate::work::LP_SRC_COUNT],
}

impl HighsMipSolver {
    /// the per-phase statistics of the solve that just ran. Reads recorded values only.
    pub fn collect_stats(&self) -> MipStats {
        let s = &self.stats_;
        let d = &self.mipdata_;
        let total_s = self.timer_.read();
        let parts = s
            .work_counters
            .as_ref()
            .map_or_else(crate::work::WorkSnapshot::now, |wc| wc.snapshot())
            .since(&s.work_start);
        let model = &self.orig_model_;
        let pre = &d.presolved_model;
        let mut separators = Vec::new();
        for (i, name) in SEPARATOR_NAMES.iter().enumerate() {
            if s.sepa_calls[i].get() > 0 || i < 3 {
                separators.push(SeparatorStats {
                    name,
                    calls: s.sepa_calls[i].get(),
                    cuts_found: s.sepa_cuts[i].get(),
                    time_s: s.sepa_ns[i].get() as f64 * 1e-9,
                });
            }
        }
        let mut heuristics = Vec::new();
        for (slot, &n) in s.improving_by_source.iter().enumerate() {
            if n > 0 {
                heuristics.push(SourceStats {
                    source: self.solution_source_to_string(slot as i32 - 1, false),
                    improving: n,
                });
            }
        }
        // The tree search starts when the root returns; a solve that ended in presolve or setup
        // never searched.
        let search_s = if s.root_end_s > 0.0 {
            (total_s - s.root_end_s).max(0.0)
        } else {
            0.0
        };
        MipStats {
            work: parts.work(),
            work_parts: parts,
            total_s,
            presolve_s: s.presolve_end_s,
            setup_s: (s.setup_end_s - s.presolve_end_s).max(0.0),
            root_node_s: s.root_node_s,
            search_s,
            presolve_status: format!("{:?}", d.presolve_status),
            num_restarts: d.num_restarts as i64,
            rows_before: model.num_row_ as i64,
            cols_before: model.num_col_ as i64,
            nnz_before: model.a_matrix_.value_.len() as i64,
            rows_after: pre.num_row_ as i64,
            cols_after: pre.num_col_ as i64,
            nnz_after: pre.a_matrix_.value_.len() as i64,
            reductions: d.post_solve_stack.reduction_counts(),
            root_lp_iterations: d.firstrootlpiters,
            root_lp_objective: d.firstlpsolobj,
            sepa_rounds: s.sepa_rounds.get(),
            separators,
            cutpool_cuts: d.cutpools.first().map_or(0, |cp| cp.get_num_cuts() as i64),
            heuristics,
            num_nodes: d.num_nodes,
            num_leaves: d.num_leaves,
            lp_iterations_total: d.total_lp_iterations,
            lp_iterations_strong_branching: d.sb_lp_iterations,
            lp_iterations_separation: d.sepa_lp_iterations,
            lp_iterations_heuristics: d.heuristic_lp_iterations,
            first_incumbent_s: s.first_incumbent_s,
            first_incumbent_work: s.first_incumbent_s.map(|_| s.first_incumbent_work),
            primal_dual_integral: d.primal_dual_integral.value,
            events: s
                .work_counters
                .as_ref()
                .map_or([0; crate::work::Event::COUNT], |wc| wc.events()),
            lp_iters_by_source: s
                .work_counters
                .as_ref()
                .map_or([0; crate::work::LP_SRC_COUNT], |wc| wc.lp_by_source().0),
            lp_solves_by_source: s
                .work_counters
                .as_ref()
                .map_or([0; crate::work::LP_SRC_COUNT], |wc| wc.lp_by_source().1),
        }
    }
}
