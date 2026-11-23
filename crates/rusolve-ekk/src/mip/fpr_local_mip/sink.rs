//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! solution_sink.h/.cpp, incumbent_sink.h/.cpp, effort_ledger.h/.cpp.
//!
//! `SolutionSink` is virtual in C++, with `IncumbentSink` as its one
//! HiGHS-side override of `on_accept`. uses one struct `S5Sink`:
//! `ms: None` is the base SolutionSink, `Some` is IncumbentSink.

use std::sync::atomic::{AtomicU64, Ordering};

use crate::io::highs_io::kHighsLogDevLevelVerbose;
use crate::lp_data::hconst::{HighsLogType, HighsVarType};
use crate::mip::fpr_local_mip::rng::S5Rng;
use crate::mip::fpr_local_mip::{
    RestartPolicy, S5HeurCharge, S5Offer, SolutionPool, WorkerTrace, kPoolCapacity,
};
use crate::mip::mip_solver::HighsMipSolver;
use crate::mip::mip_solver_data::{
    kSolutionSourceFPR, kSolutionSourceHeuristic, kSolutionSourceLocalMIP,
};
use crate::util::c_printf::{c_fmt_f, c_fmt_g};

/// mip-heuristics (36f7d9d): solution_sink.h::SolutionSink::OfferResult (lines 46-49)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct OfferResult {
    pub(crate) accepted: bool,
    pub(crate) improved_incumbent: bool,
}

/// Process-global dispatch counter (issue #106).
///
/// mip-heuristics (36f7d9d): incumbent_sink.cpp::<anonymous>::g_next_dispatch_id (line 24)
static S5_NEXT_DISPATCH_ID: AtomicU64 = AtomicU64::new(0);

/// The `kChain` name a solution-source tag belongs to.
///
/// mip-heuristics (36f7d9d): incumbent_sink.cpp::<anonymous>::name_for_source (lines 35-56).
/// The FJ / Scylla / FprLp arms are not ported: those heuristics are out of
/// scope (§0) and their patch source tags are not added to the crate, so
/// only the FPR and LocalMIP tags can reach this.
fn name_for_source(source: i32) -> &'static str {
    match source {
        s if s == kSolutionSourceFPR => "fpr",
        s if s == kSolutionSourceLocalMIP => "local_mip",
        _ => {
            // Reachable only for the generic kSolutionSourceHeuristic tag a
            // presolve sink is constructed with, which is replaced by
            // `set_source` before any worker runs.
            "unknown"
        }
    }
}

/// mip-heuristics (36f7d9d): solution_sink.h::SolutionSink (lines 24-143)
/// and incumbent_sink.h::IncumbentSink (lines 27-84), as one struct.
pub(crate) struct S5Sink<'a> {
    pool: SolutionPool,
    source: i32,
    accepted: usize,
    dispatch_name: &'static str,
    dispatch_id: u64,
    dispatch_start_s: f64,
    ms: Option<&'a mut HighsMipSolver>,
}

impl<'a> S5Sink<'a> {
    /// The base SolutionSink: a pool over columns of `integrality`, the mask
    /// diversity-aware admission measures Hamming distance on. `ms` is None,
    /// so accepted offers are counted but not submitted (unit tests).
    ///
    /// mip-heuristics (36f7d9d): solution_sink.cpp::SolutionSink::SolutionSink (lines 9-16)
    pub(crate) fn new(integrality: &[HighsVarType], source: i32) -> S5Sink<'static> {
        let mut pool = SolutionPool::new(kPoolCapacity);
        let mut int_mask: Vec<bool> = vec![false; integrality.len()];
        for j in 0..integrality.len() {
            int_mask[j] = integrality[j] != HighsVarType::kContinuous;
        }
        pool.set_integer_mask(int_mask);
        S5Sink {
            pool,
            source,
            accepted: 0,
            dispatch_name: "unknown",
            dispatch_id: 0,
            dispatch_start_s: 0.0,
            ms: None,
        }
    }

    /// Constructs the pool, seeds it from the current incumbent, and opens a
    /// dispatch. `source` tags everything offered until `set_source` says
    /// otherwise.
    ///
    /// mip-heuristics (36f7d9d): incumbent_sink.cpp::IncumbentSink::IncumbentSink (lines 60-86)
    pub(crate) fn new_incumbent(ms: &'a mut HighsMipSolver, source: i32) -> S5Sink<'a> {
        // SolutionSink(mipsolver.model_->integrality_, source).
        let mut pool = SolutionPool::new(kPoolCapacity);
        {
            let integrality = &crate::mip_model!(ms).integrality_;
            let mut int_mask: Vec<bool> = vec![false; integrality.len()];
            for j in 0..integrality.len() {
                int_mask[j] = integrality[j] != HighsVarType::kContinuous;
            }
            pool.set_integer_mask(int_mask);
        }

        // Seed the pool with the current incumbent, straight into the pool
        // rather than through `offer`: it came from HiGHS, so re-submitting
        // it would be pointless work, and it is the pool's initial state,
        // not an offer from a heuristic. Tagged with the generic
        // kSolutionSourceHeuristic so nothing downstream misattributes it.
        {
            let model = crate::mip_model!(ms);
            let incumbent = &ms.mipdata_.incumbent;
            if !incumbent.is_empty() {
                let mut obj = model.offset_;
                for j in 0..model.num_col_ as usize {
                    obj += model.col_cost_[j] * incumbent[j];
                }
                let _ = pool.try_add_slice(obj, incumbent, kSolutionSourceHeuristic);
            }
        }

        let mut sink = S5Sink {
            pool,
            source,
            accepted: 0,
            dispatch_name: "unknown",
            dispatch_id: 0,
            dispatch_start_s: 0.0,
            ms: Some(ms),
        };
        // Constructing a sink opens a dispatch: `fpr_lp` builds one per dive
        // and never calls `set_source`, so this is that path's only boundary.
        sink.begin_dispatch(source);
        sink
    }

    /// Take the next process-global dispatch id, remember the heuristic name
    /// the tag maps to, and stamp the dispatch's start on the solver clock.
    ///
    /// mip-heuristics (36f7d9d): incumbent_sink.cpp::IncumbentSink::begin_dispatch (lines 88-100)
    fn begin_dispatch(&mut self, source: i32) {
        self.source = source;
        self.dispatch_name = name_for_source(source);
        self.dispatch_id = S5_NEXT_DISPATCH_ID.fetch_add(1, Ordering::Relaxed);
        // Same clock as `EffortLedger::now_s()` — the solver's own timer,
        // not a raw steady_clock — so `[HeurSol] wall_ms` shares an origin
        // with `[Heur] start_s`.
        self.dispatch_start_s = match self.ms.as_ref() {
            Some(ms) => ms.timer_.read(),
            None => 0.0,
        };
    }

    /// Retarget the attribution tag for subsequent offers. This is also *the*
    /// dispatch boundary, and `[HeurSol]` uses it as one (#106).
    ///
    /// mip-heuristics (36f7d9d): incumbent_sink.h::IncumbentSink::set_source (line 46)
    pub(crate) fn set_source(&mut self, source: i32) {
        self.begin_dispatch(source);
    }

    /// mip-heuristics (36f7d9d): incumbent_sink.h::IncumbentSink::dispatch_id (line 53)
    pub(crate) fn dispatch_id(&self) -> u64 {
        self.dispatch_id
    }

    /// Number of offers the pool has accepted since construction.
    ///
    /// mip-heuristics (36f7d9d): solution_sink.h::SolutionSink::accepted (line 118)
    pub(crate) fn accepted(&self) -> usize {
        self.accepted
    }

    /// Restart material for a worker beginning a fresh attempt.
    ///
    /// mip-heuristics (36f7d9d): solution_sink.h::SolutionSink::get_restart (lines 122-125)
    pub(crate) fn get_restart(
        &mut self,
        rng: &mut S5Rng,
        out: &mut Vec<f64>,
        policy: RestartPolicy,
    ) -> bool {
        self.pool.get_restart(rng, out, policy)
    }

    /// mip-heuristics (36f7d9d): solution_sink.h::SolutionSink::copy_best (line 126)
    pub(crate) fn copy_best(&self, out: &mut Vec<f64>) -> bool {
        self.pool.copy_best(out)
    }

    /// The C++ termination lambda of highs_context.h:53, reached through the
    /// sink.
    pub(crate) fn terminator_terminated(&self) -> bool {
        match self.ms.as_ref() {
            Some(ms) => ms.terminator_terminated(),
            None => false,
        }
    }

    /// FPR reads the clique table and the domain through the sink.
    pub(crate) fn mipsolver(&self) -> Option<&HighsMipSolver> {
        self.ms.as_deref()
    }

    /// The driver charges the ledger through the sink's solver: in the
    /// C++ the EffortLedger and the IncumbentSink both hold the same
    /// `HighsMipSolver&` (mode_dispatch.cpp:184-247).
    pub(crate) fn mipsolver_mut(&mut self) -> Option<&mut HighsMipSolver> {
        self.ms.as_deref_mut()
    }

    /// Offer a candidate solution. When the pool accepts it, `on_accept` has
    /// already run, from inside this call.
    ///
    /// mip-heuristics (36f7d9d): solution_sink.cpp::SolutionSink::offer (lines 18-29)
    pub(crate) fn offer(
        &mut self,
        objective: f64,
        solution: &[f64],
        trace: &WorkerTrace,
        effort_at: usize,
    ) -> OfferResult {
        let added = self.pool.try_add_slice(objective, solution, self.source);
        // `accepted` and `on_accept` both stay on the admission verdict:
        // they feed `[Heur] found` and `[HeurSol] accepted`, whose meaning
        // external tooling depends on. Only the gates read the other flag.
        if added.accepted {
            self.accepted += 1;
            self.on_accept(objective, solution, self.source, trace, effort_at);
        }
        OfferResult {
            accepted: added.accepted,
            improved_incumbent: added.improved_best,
        }
    }

    /// Called once per accepted offer, on the offering worker's thread,
    /// after the pool lock is released and before `offer` returns.
    ///
    /// mip-heuristics (36f7d9d): solution_sink.cpp::SolutionSink::on_accept (lines 31-32)
    /// and incumbent_sink.cpp::IncumbentSink::on_accept (lines 102-109).
    /// The base sink does nothing; the HiGHS adapter submits to HiGHS, then
    /// emits the `[HeurSol]` line.
    fn on_accept(
        &mut self,
        objective: f64,
        solution: &[f64],
        source: i32,
        trace: &WorkerTrace,
        effort_at: usize,
    ) {
        if let Some(ms) = self.ms.as_deref_mut() {
            // Serialises trySolution: in the single-threaded path the
            // C++ mutex is not needed, but the call order is kept.
            ms.try_solution(solution, source);
            self.trace_offer(trace, effort_at, objective);
        }
    }

    /// Emit one `[HeurSol]` line.
    ///
    /// mip-heuristics (36f7d9d): incumbent_sink.cpp::IncumbentSink::trace_offer (lines 111-174)
    fn trace_offer(&mut self, trace: &WorkerTrace, effort_at: usize, objective: f64) {
        // Only accepted offers are traced (#113): this runs from on_accept.
        // Check the level before touching the clock: highsLogDev would drop
        // the line anyway, but the timer read and the formatting would
        // already have been paid on every offer of every run.
        let ms = self.ms.as_mut().unwrap();
        let log_options = &ms.options_mip_.log_options;
        if log_options.log_dev_level < kHighsLogDevLevelVerbose {
            // The report record is written whether or not the log line is
            // printed; the level check gates only the line.
            self.record_offer(trace, effort_at, objective);
            return;
        }
        let wall_ms = (ms.timer_.read() - self.dispatch_start_s) * 1000.0;
        // `%.17g` round-trips a double exactly: the objective is what the
        // calibration compares against a reference value, so it must not be
        // rounded on the way out.
        crate::highs_log_dev!(
            log_options,
            HighsLogType::kVerbose,
            "[HeurSol] name={} dispatch={} worker={} effort_at={} wall_ms={} obj={} accepted={}\n",
            self.dispatch_name,
            self.dispatch_id,
            trace.worker,
            effort_at,
            c_fmt_f(wall_ms, 1),
            c_fmt_g(objective, 17),
            1
        );
        self.record_offer(trace, effort_at, objective);
    }

    /// record the accepted offer in the solver's report, if the
    /// driver created one.
    fn record_offer(&mut self, trace: &WorkerTrace, effort_at: usize, objective: f64) {
        if let Some(ms) = self.ms.as_mut() {
            if let Some(report) = ms.mipdata_.s5_report.as_mut() {
                report.offers.push(S5Offer {
                    name: self.dispatch_name,
                    worker: trace.worker,
                    effort_at,
                    obj: objective,
                });
            }
        }
    }
}

/// The one place heuristic effort is booked and reported (issue #94).
///
/// The C++ holds `HighsMipSolver&`; makes it a method argument so the
/// driver can borrow the solver elsewhere.
///
/// mip-heuristics (36f7d9d): effort_ledger.h::EffortLedger (lines 39-106)
pub(crate) struct EffortLedger;

impl EffortLedger {
    /// Elapsed solve seconds, for the `t0_s` / `t1_s` arguments.
    ///
    /// mip-heuristics (36f7d9d): effort_ledger.cpp::EffortLedger::now_s (lines 9-11)
    pub(crate) fn now_s(&self, ms: &HighsMipSolver) -> f64 {
        ms.timer_.read()
    }

    /// A presolve-chain heuristic (FPR / LocalMIP) consumed `effort` units
    /// between `t0_s` and `t1_s`. `found` is whether the shared
    /// `IncumbentSink` accepted at least one of its solutions.
    ///
    /// `charge_dive` is not ported: fpr_lp is out of scope.
    ///
    /// mip-heuristics (36f7d9d): effort_ledger.cpp::EffortLedger::charge_presolve (lines 13-16)
    pub(crate) fn charge_presolve(
        &self,
        ms: &mut HighsMipSolver,
        name: &'static str,
        effort: usize,
        found: bool,
        t0_s: f64,
        t1_s: f64,
        abandoned_setup: bool,
    ) {
        self.book(
            ms,
            name,
            "presolve",
            effort,
            found,
            t0_s,
            t1_s,
            abandoned_setup,
        );
    }

    /// mip-heuristics (36f7d9d): effort_ledger.cpp::EffortLedger::book (lines 36-102)
    fn book(
        &self,
        ms: &mut HighsMipSolver,
        name: &'static str,
        phase: &str,
        effort: usize,
        found: bool,
        t0_s: f64,
        t1_s: f64,
        abandoned_setup: bool,
    ) {
        ms.mipdata_.heuristic_effort_used += effort;

        // The model's constraint-matrix nonzero count, taken from the same
        // buffer `make_problem` reads. Derived here rather than passed in so
        // no call site had to change.
        let nnz = ms.mipdata_.a_rindex_.len();
        let wall_ms = (t1_s - t0_s) * 1000.0;
        let effort_per_ms = if effort > 0 && wall_ms > 0.0 {
            (effort as f64) / wall_ms
        } else {
            0.0
        };
        let log_options = &ms.options_mip_.log_options;
        crate::highs_log_dev!(
            log_options,
            HighsLogType::kVerbose,
            "[Sequential] heur={} effort={} wall_ms={} effort_per_ms={}\n",
            name,
            effort,
            c_fmt_f(wall_ms, 1),
            c_fmt_f(effort_per_ms, 3)
        );
        // `nnz` and `abandoned_setup` are appended rather than inserted:
        // `[Heur]` was a positional pattern for two issues, and anything
        // still matching a prefix of it keeps working.
        crate::highs_log_dev!(
            log_options,
            HighsLogType::kVerbose,
            "[Heur] name={} phase={} start_s={} end_s={} effort={} wall_ms={} effort_per_ms={} found={} nnz={} abandoned_setup={}\n",
            name,
            phase,
            c_fmt_f(t0_s, 3),
            c_fmt_f(t1_s, 3),
            effort,
            c_fmt_f(wall_ms, 1),
            c_fmt_f(effort_per_ms, 3),
            if found { 1 } else { 0 },
            nnz,
            if abandoned_setup { 1 } else { 0 }
        );
        // the driver creates the report; book only appends to it.
        if let Some(report) = ms.mipdata_.s5_report.as_mut() {
            report.heur.push(S5HeurCharge {
                name,
                effort,
                found,
                abandoned_setup,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::highs_callback::HighsCallback;
    use crate::lp_data::highs_lp::HighsLp;
    use crate::lp_data::highs_options::HighsOptions;
    use crate::lp_data::hstruct::HighsSolution;
    use crate::mip::mip_solver_data::kSolutionSourceHeuristic;

    // A minimal MIP: two integer columns, one row x0 + x1 >= 1, min x0 + x1.
    fn tiny_mip() -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kInteger];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![f64::INFINITY];
        lp.a_matrix_.format_ = crate::lp_data::hconst::MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        let mut callback = HighsCallback::default();
        let options = HighsOptions::default();
        let solution = HighsSolution::default();
        HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0)
    }

    #[test]
    fn s5_pool_the_base_sink_counts_accepts() {
        // The base SolutionSink (ms None): accepted offers are counted,
        // nothing is submitted.
        let integrality = vec![HighsVarType::kInteger, HighsVarType::kInteger];
        let mut sink = S5Sink::new(&integrality, kSolutionSourceHeuristic);
        let trace = WorkerTrace::default();
        assert_eq!(sink.accepted(), 0);
        let r = sink.offer(10.0, &[1.0, 0.0], &trace, 0);
        assert!(r.accepted);
        assert!(r.improved_incumbent);
        assert_eq!(sink.accepted(), 1);
        let r = sink.offer(5.0, &[0.0, 1.0], &trace, 1);
        assert!(r.accepted);
        assert!(r.improved_incumbent);
        assert_eq!(sink.accepted(), 2);
        // Worse, but the pool (kPoolCapacity) has room: kept and counted,
        // not an improvement.
        let r = sink.offer(100.0, &[1.0, 1.0], &trace, 2);
        assert!(r.accepted);
        assert!(!r.improved_incumbent);
        assert_eq!(sink.accepted(), 3);
        // Fill the pool, then a dominated offer outside the diversity
        // tolerance is refused and not counted.
        for k in 3..kPoolCapacity {
            assert!(sink.offer(20.0 + k as f64, &[0.0, 0.0], &trace, 3).accepted);
        }
        assert_eq!(sink.accepted(), kPoolCapacity as usize);
        let r = sink.offer(1000.0, &[1.0, 1.0], &trace, 4);
        assert!(!r.accepted);
        assert!(!r.improved_incumbent);
        assert_eq!(sink.accepted(), kPoolCapacity as usize);
        // A non-finite objective is refused outright.
        assert!(!sink.offer(f64::NAN, &[0.0, 0.0], &trace, 5).accepted);
        assert_eq!(sink.accepted(), kPoolCapacity as usize);
    }

    #[test]
    fn s5_pool_the_incumbent_sink_seeds_the_pool_from_the_incumbent() {
        // The IncumbentSink ctor seeds the pool from the incumbent with
        // obj = offset_ + sum_j col_cost_[j] * incumbent[j] (incumbent_sink.
        // cpp:73-81). With an empty incumbent the pool starts empty.
        let mut ms = tiny_mip();
        assert!(ms.mipdata_.incumbent.is_empty());
        {
            let sink = S5Sink::new_incumbent(&mut ms, kSolutionSourceHeuristic);
            assert_eq!(sink.pool.size(), 0);
            assert_eq!(sink.accepted(), 0);
            assert_eq!(sink.dispatch_name, "unknown");
        }
        // With an incumbent, the pool holds it, tagged generic.
        ms.mipdata_.incumbent = vec![1.0, 0.0];
        {
            let sink = S5Sink::new_incumbent(&mut ms, kSolutionSourceHeuristic);
            assert_eq!(sink.pool.size(), 1);
            let entries = sink.pool.sorted_entries();
            assert_eq!(entries[0].objective, 1.0); // offset 0 + 1*1 + 1*0
            assert_eq!(entries[0].source, kSolutionSourceHeuristic);
            assert_eq!(entries[0].solution, vec![1.0, 0.0]);
            // Seeding is not an offer: nothing is counted as production.
            assert_eq!(sink.accepted(), 0);
        }
    }
}
