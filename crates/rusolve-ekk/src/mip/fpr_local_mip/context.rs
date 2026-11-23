//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! highs_context.h, the restart helpers of heuristic_context.cpp,
//! worker_base.h, continuous_loop.h, opportunistic_runner.h (N = 1),
//! least_violated.h and deadline.h.
//!
//! The sequential path (N = 1, Q11): the runner's parallel::for_each is
//! dropped and slots run in index order on this thread.  `Deadline`'s
//! timer is always `None`, so `expired()` is always false and every
//! deadline poll is a literal port that never fires.

use std::rc::Rc;

use crate::HighsInt;
use crate::lp_data::hconst::{HighsVarType, ObjSense, kHighsInf};
use crate::mip::fpr_local_mip::rng::{
    S5Rng, build_csc, heuristic_base_seed, s5_bernoulli, s5_uniform_int_i64, s5_uniform_real,
};
use crate::mip::fpr_local_mip::sink::S5Sink;
use crate::mip::fpr_local_mip::{
    AttemptResult, CscMatrix, Deadline, ExecutionContext, HeuristicBudget,
    K_SAFE_INT64_DOUBLE_RANGE, RestartPoint, RestartPolicy, S5Problem, WorkerBudgetState,
    WorkerTrace, kNoPoller,
};
use crate::mip::mip_solver::HighsMipSolver;
use crate::util::highs_utils::{cpp_max, cpp_min};

// ---------------------------------------------------------------------------
// heuristic_context.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): heuristic_context.h::RestartSource (line 631)
///
/// The C++ `std::function<bool(Rng&, std::vector<double>&)>`. The presolve
/// path passes `None` (the C++ empty std::function).
pub(crate) type RestartSource<'a> = dyn FnMut(&mut S5Rng, &mut Vec<f64>) -> bool + 'a;

// ---------------------------------------------------------------------------
// deadline.h
// ---------------------------------------------------------------------------

impl Deadline {
    /// mip-heuristics (36f7d9d): deadline.h::Deadline::expired (line 33)
    ///
    /// `timer` is always `None`, so this is always false; the `Some`
    /// arm is unreachable because `Infallible` is uninhabited.
    #[must_use]
    pub(crate) fn expired(&self) -> bool {
        match self.timer {
            Some(t) => match t {},
            None => false,
        }
    }

    /// mip-heuristics (36f7d9d): deadline.h::Deadline::remaining (lines 48-54)
    ///
    /// Seconds left, `kHighsInf` when there is no limit, never negative.
    /// With `timer` always `None` this is always `kHighsInf`.
    #[must_use]
    pub(crate) fn remaining(&self) -> f64 {
        match self.timer {
            Some(t) => match t {},
            None => kHighsInf,
        }
    }
}

// ---------------------------------------------------------------------------
// heuristic_context.h::ExecutionContext (methods)
// ---------------------------------------------------------------------------

impl ExecutionContext {
    /// mip-heuristics (36f7d9d): heuristic_context.h::ExecutionContext::past_deadline (lines 482-484)
    ///
    /// The write-free half of "should we stop?": the wall-clock deadline.
    /// There is no stop flag in the struct, so this is the deadline
    /// alone.
    #[must_use]
    pub(crate) fn past_deadline(&self) -> bool {
        self.deadline.expired()
    }

    /// mip-heuristics (36f7d9d): heuristic_context.h::ExecutionContext::terminated (lines 501-503)
    ///
    /// The full "should we stop?" predicate. The terminator is reached
    /// through the sink; it may write, so this one needs a single
    /// caller.
    #[must_use]
    pub(crate) fn terminated(&self, sink: &S5Sink<'_>) -> bool {
        sink.terminator_terminated() || self.past_deadline()
    }

    /// mip-heuristics (36f7d9d): heuristic_context.h::ExecutionContext::worker_seed (lines 509-511)
    ///
    /// Deterministic seed for worker `w`.
    #[must_use]
    pub(crate) fn worker_seed(&self, w: i32) -> u32 {
        crate::mip::fpr_local_mip::rng::worker_seed(self.base_seed, w)
    }
}

// ---------------------------------------------------------------------------
// highs_context.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): highs_context.h::build_binary_mask (lines 25-33)
///
/// Snapshot `HighsDomain::isBinary` for every column.
pub(crate) fn build_binary_mask(ms: &HighsMipSolver) -> Vec<u8> {
    let ncol = crate::mip_model!(ms).num_col_;
    let domain = ms.get_domain_2();
    let mut mask = vec![0u8; ncol as usize];
    for j in 0..ncol {
        mask[j as usize] = u8::from(domain.is_binary(ms, j));
    }
    mask
}

/// mip-heuristics (36f7d9d): highs_context.h::problem_view + make_problem (lines 62-109)
///
/// The owned snapshot (D3/R6): the row-major buffers plus the CSC transpose
/// built from them, with the incumbent and `isBinary` snapshots. The
/// presolved model is already in minimisation form, so the view points at
/// its costs unchanged.
pub(crate) fn make_problem(ms: &HighsMipSolver) -> S5Problem {
    let model = crate::mip_model!(ms);
    let mipdata = &*ms.mipdata_;
    let ncol = model.num_col_;
    let nrow = model.num_row_;
    debug_assert!(model.sense_ == ObjSense::kMinimize);
    let ar_start = mipdata.a_rstart_.clone();
    let ar_index = mipdata.a_rindex_.clone();
    let ar_value = mipdata.a_rvalue_.clone();
    let csc = build_csc(ncol, nrow, &ar_start, &ar_index, &ar_value);
    let binary = build_binary_mask(ms);
    S5Problem {
        col_lower: model.col_lower_.clone(),
        col_upper: model.col_upper_.clone(),
        row_lower: model.row_lower_.clone(),
        row_upper: model.row_upper_.clone(),
        col_cost: model.col_cost_.clone(),
        offset: model.offset_,
        integrality: model.integrality_.clone(),
        ar_start,
        ar_index,
        ar_value,
        csc,
        uplocks: mipdata.uplocks.clone(),
        downlocks: mipdata.downlocks.clone(),
        feastol: mipdata.feastol,
        epsilon: mipdata.epsilon,
        ncol,
        nrow,
        nnz: mipdata.a_rindex_.len(),
        incumbent: mipdata.incumbent.clone(),
        binary,
    }
}

/// mip-heuristics (36f7d9d): highs_context.h::make_exec (lines 47-54)
///
/// One dispatch's execution parameters. `num_workers` is fixed at 1;
/// the deadline has no timer; the terminator is reached through the
/// sink, so it is not a field here.
pub(crate) fn make_exec(ms: &HighsMipSolver) -> ExecutionContext {
    ExecutionContext {
        num_workers: 1,
        base_seed: heuristic_base_seed(ms.options_mip_.random_seed),
        deadline: Deadline::default(),
        // rusolve GUB swap moves (`gub_swap.rs`)
        gub_swap: ms.options_mip_.mip_ls_gub_swap,
    }
}

// ---------------------------------------------------------------------------
// heuristic_context.cpp: restart helpers
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): heuristic_context.cpp::take_restart (lines 305-323)
///
/// `source`'s point into `out`, rounded and clamped into `problem`, when
/// there is a source and its point is usable; otherwise false, with `out`
/// empty.
#[must_use]
pub(crate) fn take_restart(
    source: Option<&mut RestartSource<'_>>,
    rng: &mut S5Rng,
    out: &mut Vec<f64>,
    problem: &S5Problem,
) -> bool {
    if let Some(src) = source {
        if src(rng, out) && out.len() == problem.ncol as usize {
            let mut usable = true;
            let mut j: HighsInt = 0;
            while j < problem.ncol && usable {
                let mut v = out[j as usize];
                usable = v.is_finite();
                if problem.integrality[j as usize] != HighsVarType::kContinuous {
                    v = v.round();
                }
                out[j as usize] = cpp_max(
                    problem.col_lower[j as usize],
                    cpp_min(problem.col_upper[j as usize], v),
                );
                j += 1;
            }
            if usable {
                return true;
            }
        }
    }
    out.clear();
    false
}

/// mip-heuristics (36f7d9d): heuristic_context.cpp::hybrid_restart (lines 325-352)
///
/// Local-MIP v2's `hybrid_restart`, in place on `point`: each column keeps
/// its value or, on an even coin, takes its `sample_random_value` — uniform
/// over the column's integers when it is an integer column with both bounds
/// finite, uniform over its interval when it is continuous with both bounds
/// finite, and its value otherwise. A finite integer range wider than
/// `K_SAFE_INT64_DOUBLE_RANGE` counts as unbounded.
pub(crate) fn hybrid_restart(point: &mut [f64], problem: &S5Problem, rng: &mut S5Rng) {
    let mut j: HighsInt = 0;
    while j < problem.ncol {
        if !s5_bernoulli(rng, 0.5) {
            j += 1;
            continue;
        }
        let lb = problem.col_lower[j as usize];
        let ub = problem.col_upper[j as usize];
        if !lb.is_finite() || !ub.is_finite() {
            j += 1;
            continue;
        }
        if problem.integrality[j as usize] == HighsVarType::kContinuous {
            point[j as usize] = s5_uniform_real(rng, lb, ub);
        } else {
            // The core `make_problem` rounds integer bounds inward; the
            // adapter's view points at the model's own, which need not be
            // integral, so round them here as the reference does.
            let lo = lb.ceil();
            let hi = ub.floor();
            if hi < lo || hi - lo > K_SAFE_INT64_DOUBLE_RANGE {
                j += 1;
                continue;
            }
            let span = (hi - lo) as i64;
            point[j as usize] = lo + (s5_uniform_int_i64(rng, 0, span) as f64);
        }
        j += 1;
    }
}

/// mip-heuristics (36f7d9d): heuristic_context.cpp::rebuild_restart_point (lines 354-371)
///
/// The point a rebuilt worker restarts from under `policy`, into `out`, or
/// false — with `out` empty and nothing drawn — under `kDefault`, or when
/// there is no base point. The base is `source`'s point, else the pool's
/// best entry (drawing nothing), else `problem.incumbent`; `kHybrid` then
/// applies `hybrid_restart`.
#[must_use]
pub(crate) fn rebuild_restart_point(
    policy: RestartPoint,
    source: Option<&mut RestartSource<'_>>,
    sink: &mut S5Sink<'_>,
    rng: &mut S5Rng,
    out: &mut Vec<f64>,
    problem: &S5Problem,
) -> bool {
    out.clear();
    if policy == RestartPoint::kDefault {
        return false;
    }
    if !take_restart(source, rng, out, problem)
        && !sink.get_restart(rng, out, RestartPolicy::kBestOnly)
    {
        if problem.incumbent.is_empty() {
            return false;
        }
        *out = problem.incumbent.clone();
    }
    if policy == RestartPoint::kHybrid {
        hybrid_restart(out, problem, rng);
    }
    true
}

// ---------------------------------------------------------------------------
// worker_base.h
// ---------------------------------------------------------------------------

impl WorkerTrace {
    /// mip-heuristics (36f7d9d): worker_base.h::WorkerTrace::at (line 48)
    ///
    /// Monotone charged effort to report, given the current occupant's own
    /// running counter.
    #[must_use]
    pub(crate) fn at(&self, charged: usize) -> usize {
        self.effort_base.wrapping_add(charged)
    }
}

impl WorkerBudgetState {
    /// mip-heuristics (36f7d9d): worker_base.h::WorkerBudgetState::stale (line 88)
    ///
    /// True when this worker has exceeded its staleness budget.
    #[must_use]
    pub(crate) fn stale(&self) -> bool {
        self.effort_since_improvement > self.stale_budget
    }

    /// mip-heuristics (36f7d9d): worker_base.h::WorkerBudgetState::stale(size_t) (lines 93-95)
    ///
    /// True when already stale, or would become stale after `extra` more
    /// effort.
    #[must_use]
    pub(crate) fn stale_with(&self, extra: usize) -> bool {
        self.effort_since_improvement.wrapping_add(extra) > self.stale_budget
    }

    /// mip-heuristics (36f7d9d): worker_base.h::WorkerBudgetState::exhausted (line 98)
    ///
    /// True when this worker has consumed its total budget.
    #[must_use]
    pub(crate) fn exhausted(&self) -> bool {
        self.total_effort >= self.total_budget
    }

    /// mip-heuristics (36f7d9d): worker_base.h::WorkerBudgetState::exhausted(size_t) (lines 102-104)
    ///
    /// True when already exhausted, or would become exhausted after `extra`
    /// more effort.
    #[must_use]
    pub(crate) fn exhausted_with(&self, extra: usize) -> bool {
        self.total_effort.wrapping_add(extra) >= self.total_budget
    }

    /// mip-heuristics (36f7d9d): worker_base.h::WorkerBudgetState::reset_staleness (line 109)
    ///
    /// Clear the staleness counter; called by the worker itself when an
    /// offer moved the incumbent.
    pub(crate) fn reset_staleness(&mut self) {
        self.effort_since_improvement = 0;
    }

    /// mip-heuristics (36f7d9d): worker_base.h::WorkerBudgetState::charge_improvement (lines 113-119)
    ///
    /// Accumulate effort when the worker moved the incumbent. Resets
    /// staleness and marks finished if total budget is exceeded.
    pub(crate) fn charge_improvement(&mut self, effort: usize) {
        self.total_effort = self.total_effort.wrapping_add(effort);
        self.effort_since_improvement = 0;
        if self.exhausted() {
            self.finished = true;
        }
    }

    /// mip-heuristics (36f7d9d): worker_base.h::WorkerBudgetState::charge_no_improvement (lines 124-130)
    ///
    /// Accumulate effort when the worker did not move the incumbent.
    /// Advances the staleness counter and marks finished if either budget
    /// is exceeded.
    pub(crate) fn charge_no_improvement(&mut self, effort: usize) {
        self.total_effort = self.total_effort.wrapping_add(effort);
        self.effort_since_improvement = self.effort_since_improvement.wrapping_add(effort);
        if self.exhausted() || self.stale() {
            self.finished = true;
        }
    }
}

// ---------------------------------------------------------------------------
// continuous_loop.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState (lines 47-181)
///
/// Shared scaffold for continuous-parallel heuristic runners. The atomics
/// are plain fields here: the sequential path has one worker on one
/// thread.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ContinuousLoopState {
    pub(crate) total_effort: usize,
    pub(crate) effort_since_improvement: usize,
    pub(crate) stop: bool,
    // Index of the worker currently responsible for termination polling, or
    // kNoPoller when the seat is vacant (continuous_loop.h:52-64).
    pub(crate) poller: i32,
}

impl ContinuousLoopState {
    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::ContinuousLoopState (line 66)
    pub(crate) fn new() -> Self {
        ContinuousLoopState {
            total_effort: 0,
            effort_since_improvement: 0,
            stop: false,
            poller: 0,
        }
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::ContinuousLoopState(int) (line 70)
    ///
    /// A loop whose seat starts with worker `first_poller` rather than 0.
    pub(crate) fn with_first_poller(first_poller: i32) -> Self {
        ContinuousLoopState {
            total_effort: 0,
            effort_since_improvement: 0,
            stop: false,
            poller: first_poller,
        }
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::stopped (line 72)
    #[must_use]
    pub(crate) fn stopped(&self) -> bool {
        self.stop
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::request_stop (line 74)
    pub(crate) fn request_stop(&mut self) {
        self.stop = true;
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::claim_poller (lines 94-104)
    ///
    /// True when `w` should poll: it already holds the seat, or the seat is
    /// vacant and `w` takes it. The seat keeps the HiGHS terminator
    /// single-caller.
    pub(crate) fn claim_poller(&mut self, w: i32) -> bool {
        let cur = self.poller;
        if cur == w {
            return true;
        }
        if cur != kNoPoller {
            return false;
        }
        self.poller = w;
        true
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::release_poller (lines 108-112)
    ///
    /// Called by a worker on its way out of the loop: vacate the seat if it
    /// holds it, so a surviving peer can take over the polling duty.
    pub(crate) fn release_poller(&mut self, w: i32) {
        if self.poller == w {
            self.poller = kNoPoller;
        }
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::poll_termination (lines 117-121)
    ///
    /// Seat-holder only — the underlying HiGHS calls are not thread-safe
    /// for concurrent callers.
    pub(crate) fn poll_termination(&mut self, exec: &ExecutionContext, sink: &S5Sink<'_>) {
        if exec.terminated(sink) {
            self.request_stop();
        }
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::poll_deadline (lines 127-136)
    ///
    /// Callable from any worker thread, unlike `poll_termination`: the
    /// clock read writes nothing (issue #114). Returns whether the loop
    /// should stop.
    pub(crate) fn poll_deadline(&mut self, exec: &ExecutionContext) -> bool {
        if self.stopped() {
            return true;
        }
        if exec.past_deadline() {
            self.request_stop();
            return true;
        }
        false
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::should_stop (lines 149-157)
    ///
    /// "Has this dispatch been told to stop?", for worker `w` about to
    /// begin its `attempt_counter`-th attempt: the shared flag, the
    /// wall-clock deadline (every iteration, any thread), and the
    /// terminator (every other attempt, seat-holder only).
    pub(crate) fn should_stop(
        &mut self,
        exec: &ExecutionContext,
        sink: &S5Sink<'_>,
        w: i32,
        attempt_counter: i32,
    ) -> bool {
        if self.poll_deadline(exec) {
            return true;
        }
        if (attempt_counter & 1) == 0 && self.claim_poller(w) {
            self.poll_termination(exec, sink);
        }
        self.stopped()
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::add_effort (lines 161-167)
    ///
    /// Bump the cumulative-effort counter and set `stop` if the cumulative
    /// total crossed `budget`. Returns the new total.
    pub(crate) fn add_effort(&mut self, effort: usize, budget: usize) -> usize {
        let new_total = self.total_effort.wrapping_add(effort);
        self.total_effort = new_total;
        if new_total >= budget {
            self.request_stop();
        }
        new_total
    }

    /// mip-heuristics (36f7d9d): continuous_loop.h::ContinuousLoopState::note_staleness (lines 171-180)
    ///
    /// Update the staleness counter and set `stop` if it crossed
    /// `stale_budget`. Staleness is advisory.
    pub(crate) fn note_staleness(&mut self, effort: usize, improved: bool, stale_budget: usize) {
        if improved {
            self.effort_since_improvement = 0;
        } else {
            self.effort_since_improvement = self.effort_since_improvement.wrapping_add(effort);
        }
        if self.effort_since_improvement >= stale_budget {
            self.request_stop();
        }
    }
}

// ---------------------------------------------------------------------------
// least_violated.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): least_violated.h::LeastViolatedPoint::Point (lines 39-42)
///
/// The retained point and its violation. The Rust port is single-threaded
///, so the C++ mutex and atomic are dropped.
#[derive(Clone, Debug)]
pub(crate) struct LeastViolatedEntry {
    pub(crate) x: Rc<Vec<f64>>,
    pub(crate) violation: f64,
}

/// mip-heuristics (36f7d9d): least_violated.h::LeastViolatedPoint (lines 37-83)
///
/// The least-violated infeasible point a local search has retained (#188).
/// Caller-owned and opt-in: leave it unused and nothing is measured or
/// copied.
#[derive(Clone, Debug)]
pub(crate) struct LeastViolatedPoint {
    x: Option<Rc<Vec<f64>>>,
    violation: f64,
    // `violation`, readable without the lock; a stale read only lets an
    // offer through to the re-check (least_violated.h:80-82).
    best_violation: f64,
}

impl Default for LeastViolatedPoint {
    // least_violated.h:79, :82
    fn default() -> Self {
        LeastViolatedPoint {
            x: None,
            violation: f64::INFINITY,
            best_violation: f64::INFINITY,
        }
    }
}

impl LeastViolatedPoint {
    /// mip-heuristics (36f7d9d): least_violated.h::LeastViolatedPoint::best (lines 45-51)
    ///
    /// The retained point, or nothing before the first offer.
    #[must_use]
    pub(crate) fn best(&self) -> Option<LeastViolatedEntry> {
        self.x.as_ref().map(|x| LeastViolatedEntry {
            x: Rc::clone(x),
            violation: self.violation,
        })
    }

    /// mip-heuristics (36f7d9d): least_violated.h::LeastViolatedPoint::would_take (lines 55-57)
    ///
    /// Whether an offer at `violation` would be taken, so a worker can skip
    /// building one that would not.
    #[must_use]
    pub(crate) fn would_take(&self, violation: f64) -> bool {
        violation <= self.best_violation
    }

    /// mip-heuristics (36f7d9d): least_violated.h::LeastViolatedPoint::offer (lines 61-74)
    ///
    /// Retain `x` when `violation` is at most the retained one's; returns
    /// whether it did.
    pub(crate) fn offer(&mut self, violation: f64, x: &[f64]) -> bool {
        if !self.would_take(violation) {
            return false;
        }
        let point = Rc::new(x.to_vec());
        if violation > self.violation {
            return false;
        }
        self.x = Some(point);
        self.violation = violation;
        self.best_violation = violation;
        true
    }
}

// ---------------------------------------------------------------------------
// opportunistic_runner.h (N = 1)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): opportunistic_runner.h::attempt_with_rebuild (lines 58-69)
///
/// Run one attempt, rebuilding the worker first when it is finished and
/// again when a zero-effort attempt finds it finished.
pub(crate) fn attempt_with_rebuild<W>(
    worker: &mut W,
    run_cap: usize,
    finished: impl Fn(&W) -> bool,
    mut run_attempt: impl FnMut(&mut W, usize) -> AttemptResult,
    mut rebuild: impl FnMut(&mut W),
) -> AttemptResult {
    if finished(worker) {
        rebuild(worker);
    }
    let mut attempt = run_attempt(worker, run_cap);
    if attempt.effort == 0 && finished(worker) {
        rebuild(worker);
        attempt = run_attempt(worker, run_cap);
    }
    attempt
}

/// mip-heuristics (36f7d9d): opportunistic_runner.h::run_slot (lines 77-127)
///
/// One worker slot's loop, on the calling thread: build the slot's state,
/// then run attempts until `loop` says stop, the budget is spent or an
/// attempt charges nothing.
pub(crate) fn run_slot<State, MakeState, RunAttempt>(
    cl: &mut ContinuousLoopState,
    exec: &ExecutionContext,
    sink: &mut S5Sink<'_>,
    budget: &HeuristicBudget,
    w: i32,
    make_state: &mut MakeState,
    run_attempt: &mut RunAttempt,
) where
    MakeState: FnMut(i32, &mut S5Rng, &mut S5Sink<'_>) -> State,
    RunAttempt: FnMut(&mut State, &mut S5Rng, usize, &mut S5Sink<'_>) -> AttemptResult,
{
    let mut rng = S5Rng::new(exec.worker_seed(w) as u64);
    let mut attempt_counter: i32 = 0;

    let mut state = make_state(w, &mut rng, sink);

    while !cl.stopped() {
        // Flag, deadline and terminator in one predicate — see
        // `ContinuousLoopState::should_stop`. The deadline half runs on
        // every iteration and from this thread (issue #114).
        if cl.should_stop(exec, &*sink, w, attempt_counter) {
            break;
        }

        let current = cl.total_effort;
        let remaining = budget.total - cpp_min(budget.total, current);
        let run_cap = cpp_min(budget.attempt_cap, remaining);
        if run_cap == 0 {
            cl.request_stop();
            break;
        }

        let result = run_attempt(&mut state, &mut rng, run_cap, sink);
        attempt_counter += 1;

        // Guard against workers that make no progress: a zero-effort
        // return means this worker is done. Retire it on its own (see
        // `attempt_with_rebuild`).
        if result.effort == 0 {
            break;
        }

        cl.note_staleness(result.effort, result.found_improvement, budget.stale);
        cl.add_effort(result.effort, budget.total);
    }

    cl.release_poller(w);
}

/// mip-heuristics (36f7d9d): opportunistic_runner.h::run_opportunistic_loop (lines 130-151)
///
/// Returns total effort consumed across all workers. The parallel::for_each
/// is dropped: the slots run in index order on this thread.
#[must_use]
pub(crate) fn run_opportunistic_loop<State, MakeState, RunAttempt>(
    exec: &ExecutionContext,
    sink: &mut S5Sink<'_>,
    budget: &HeuristicBudget,
    mut make_state: MakeState,
    mut run_attempt: RunAttempt,
) -> usize
where
    MakeState: FnMut(i32, &mut S5Rng, &mut S5Sink<'_>) -> State,
    RunAttempt: FnMut(&mut State, &mut S5Rng, usize, &mut S5Sink<'_>) -> AttemptResult,
{
    let n = exec.num_workers as i64;
    if n <= 0 || budget.total == 0 {
        return 0;
    }

    let mut cl = ContinuousLoopState::new();

    let mut w: i32 = 0;
    while (w as i64) < n {
        run_slot(
            &mut cl,
            exec,
            sink,
            budget,
            w,
            &mut make_state,
            &mut run_attempt,
        );
        w += 1;
    }

    cl.total_effort
}

/// mip-heuristics (36f7d9d): opportunistic_runner.h::run_on_caller_thread (lines 159-168)
///
/// Run one slot, of any index, on the caller's thread; that slot must hold
/// the poller seat or nothing would ever poll the terminator.
#[must_use]
pub(crate) fn run_on_caller_thread<State, MakeState, RunAttempt>(
    exec: &ExecutionContext,
    sink: &mut S5Sink<'_>,
    budget: &HeuristicBudget,
    worker: i32,
    mut make_state: MakeState,
    mut run_attempt: RunAttempt,
) -> usize
where
    MakeState: FnMut(i32, &mut S5Rng, &mut S5Sink<'_>) -> State,
    RunAttempt: FnMut(&mut State, &mut S5Rng, usize, &mut S5Sink<'_>) -> AttemptResult,
{
    if budget.total == 0 {
        return 0;
    }
    let mut cl = ContinuousLoopState::with_first_poller(worker);
    run_slot(
        &mut cl,
        exec,
        sink,
        budget,
        worker,
        &mut make_state,
        &mut run_attempt,
    );
    cl.total_effort
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;
    use crate::mip::fpr_local_mip::rng::{make_budget, make_until_stopped_budget};

    fn sink() -> S5Sink<'static> {
        S5Sink::new(&[], 0)
    }

    // -----------------------------------------------------------------------
    // continuous_loop.h
    // -----------------------------------------------------------------------

    #[test]
    fn s5_ctx_continuous_loop_add_effort_accumulates() {
        // test_continuous_loop.cpp: "single-threaded add_effort accumulates correctly"
        let mut cl = ContinuousLoopState::new();
        assert_eq!(cl.add_effort(10, 100), 10);
        assert_eq!(cl.add_effort(20, 100), 30);
        assert!(!cl.stopped());
        assert_eq!(cl.add_effort(80, 100), 110);
        assert!(cl.stopped());
    }

    #[test]
    fn s5_ctx_continuous_loop_note_staleness() {
        // test_continuous_loop.cpp: the two SECTIONs
        let mut cl = ContinuousLoopState::new();
        cl.note_staleness(10, false, 50);
        assert_eq!(cl.effort_since_improvement, 10);
        assert!(!cl.stopped());
        cl.note_staleness(20, false, 50);
        assert_eq!(cl.effort_since_improvement, 30);
        assert!(!cl.stopped());
        cl.note_staleness(30, false, 50);
        assert_eq!(cl.effort_since_improvement, 60);
        assert!(cl.stopped());

        let mut cl = ContinuousLoopState::new();
        cl.note_staleness(10, false, 50);
        cl.note_staleness(20, true, 50);
        assert_eq!(cl.effort_since_improvement, 0);
        assert!(!cl.stopped());
    }

    #[test]
    fn s5_ctx_continuous_loop_request_stop_irreversible() {
        // test_continuous_loop.cpp
        let mut cl = ContinuousLoopState::new();
        assert!(!cl.stopped());
        cl.request_stop();
        assert!(cl.stopped());
        cl.request_stop();
        assert!(cl.stopped());
    }

    #[test]
    fn s5_ctx_continuous_loop_poller_seat_exclusive() {
        // test_continuous_loop.cpp
        let mut cl = ContinuousLoopState::new();
        assert!(cl.claim_poller(0));
        assert!(cl.claim_poller(0));
        assert!(!cl.claim_poller(1));
        assert!(!cl.claim_poller(2));
    }

    #[test]
    fn s5_ctx_continuous_loop_release_hands_seat_over() {
        // test_continuous_loop.cpp
        let mut cl = ContinuousLoopState::new();
        assert!(cl.claim_poller(0));
        cl.release_poller(0);
        assert!(cl.claim_poller(1));
        assert!(!cl.claim_poller(2));
        // A worker that does not hold the seat cannot vacate it.
        cl.release_poller(0);
        assert!(cl.claim_poller(1));
        assert!(!cl.claim_poller(2));
    }

    // -----------------------------------------------------------------------
    // deadline.h
    // -----------------------------------------------------------------------

    #[test]
    fn s5_ctx_deadline_never_expires() {
        // test_deadline.cpp spirit: the S5 Deadline has no timer.
        let d = Deadline::default();
        assert!(!d.expired());
        assert_eq!(d.remaining(), kHighsInf);
    }

    // -----------------------------------------------------------------------
    // opportunistic_runner.h
    // -----------------------------------------------------------------------

    #[test]
    fn s5_ctx_run_slot_stale_budget_stops() {
        let mut cl = ContinuousLoopState::new();
        let exec = ExecutionContext::default();
        let mut sink = sink();
        let budget = make_budget(1_000_000, 1, 50);
        let mut calls = 0usize;
        let mut make_state = |_w: i32, _rng: &mut S5Rng, _sink: &mut S5Sink<'_>| ();
        let mut run_attempt = |_state: &mut (),
                               _rng: &mut S5Rng,
                               _cap: usize,
                               _sink: &mut S5Sink<'_>|
         -> AttemptResult {
            calls += 1;
            AttemptResult {
                effort: 10,
                found_improvement: false,
            }
        };
        run_slot(
            &mut cl,
            &exec,
            &mut sink,
            &budget,
            0,
            &mut make_state,
            &mut run_attempt,
        );
        assert_eq!(calls, 5);
        assert_eq!(cl.total_effort, 50);
        assert!(cl.stopped());
    }

    #[test]
    fn s5_ctx_run_slot_attempt_cap_bounds_run_cap() {
        let mut cl = ContinuousLoopState::new();
        let exec = ExecutionContext::default();
        let mut sink = sink();
        let budget = make_until_stopped_budget(7, 0);
        let mut caps: Vec<usize> = Vec::new();
        let mut calls = 0usize;
        let mut make_state = |_w: i32, _rng: &mut S5Rng, _sink: &mut S5Sink<'_>| ();
        let mut run_attempt = |_state: &mut (),
                               _rng: &mut S5Rng,
                               cap: usize,
                               _sink: &mut S5Sink<'_>|
         -> AttemptResult {
            caps.push(cap);
            calls += 1;
            if calls <= 3 {
                AttemptResult {
                    effort: 7,
                    found_improvement: true,
                }
            } else {
                AttemptResult {
                    effort: 0,
                    found_improvement: false,
                }
            }
        };
        run_slot(
            &mut cl,
            &exec,
            &mut sink,
            &budget,
            0,
            &mut make_state,
            &mut run_attempt,
        );
        assert_eq!(caps, vec![7, 7, 7, 7]);
        assert_eq!(calls, 4);
    }

    #[test]
    fn s5_ctx_run_slot_total_cap_stops() {
        let mut cl = ContinuousLoopState::new();
        let exec = ExecutionContext::default();
        let mut sink = sink();
        let budget = make_budget(25, 1, usize::MAX);
        let mut calls = 0usize;
        let mut make_state = |_w: i32, _rng: &mut S5Rng, _sink: &mut S5Sink<'_>| ();
        let mut run_attempt = |_state: &mut (),
                               _rng: &mut S5Rng,
                               cap: usize,
                               _sink: &mut S5Sink<'_>|
         -> AttemptResult {
            calls += 1;
            AttemptResult {
                effort: cap,
                found_improvement: true,
            }
        };
        run_slot(
            &mut cl,
            &exec,
            &mut sink,
            &budget,
            0,
            &mut make_state,
            &mut run_attempt,
        );
        assert_eq!(calls, 13);
        assert_eq!(cl.total_effort, 25);
        assert!(cl.stopped());
    }

    #[test]
    fn s5_ctx_run_opportunistic_loop_sums_effort() {
        // N = 1 via ExecutionContext::default().
        let exec = ExecutionContext::default();
        let mut sink = sink();
        let budget = make_budget(25, 1, usize::MAX);
        let mut calls = 0usize;
        let make_state = |_w: i32, _rng: &mut S5Rng, _sink: &mut S5Sink<'_>| ();
        let mut run_attempt = |_state: &mut (),
                               _rng: &mut S5Rng,
                               cap: usize,
                               _sink: &mut S5Sink<'_>|
         -> AttemptResult {
            calls += 1;
            AttemptResult {
                effort: cap,
                found_improvement: true,
            }
        };
        let total = run_opportunistic_loop(&exec, &mut sink, &budget, make_state, &mut run_attempt);
        assert_eq!(calls, 13);
        assert_eq!(total, 25);
    }

    #[test]
    fn s5_ctx_run_on_caller_thread_zero_budget_is_noop() {
        let exec = ExecutionContext::default();
        let mut sink = sink();
        let budget = make_budget(0, 1, 0);
        let mut made = false;
        let make_state = |_w: i32, _rng: &mut S5Rng, _sink: &mut S5Sink<'_>| {
            made = true;
        };
        let mut run_attempt = |_state: &mut (),
                               _rng: &mut S5Rng,
                               _cap: usize,
                               _sink: &mut S5Sink<'_>|
         -> AttemptResult {
            AttemptResult {
                effort: 1,
                found_improvement: false,
            }
        };
        let total =
            run_on_caller_thread(&exec, &mut sink, &budget, 0, make_state, &mut run_attempt);
        assert_eq!(total, 0);
        assert!(!made);
    }

    #[test]
    fn s5_ctx_attempt_with_rebuild_retries_after_zero() {
        struct W {
            finished: bool,
            rebuilds: usize,
            runs: usize,
        }
        let mut w = W {
            finished: true,
            rebuilds: 0,
            runs: 0,
        };
        let result = attempt_with_rebuild(
            &mut w,
            10,
            |w| w.finished,
            |w, _cap| {
                w.runs += 1;
                if w.runs == 1 {
                    // The first attempt retires the worker and charges nothing.
                    w.finished = true;
                    AttemptResult {
                        effort: 0,
                        found_improvement: false,
                    }
                } else {
                    AttemptResult {
                        effort: 5,
                        found_improvement: false,
                    }
                }
            },
            |w| {
                w.rebuilds += 1;
                w.finished = false;
            },
        );
        assert_eq!(result.effort, 5);
        assert_eq!(w.rebuilds, 2);
        assert_eq!(w.runs, 2);
    }

    // -----------------------------------------------------------------------
    // least_violated.h
    // -----------------------------------------------------------------------

    #[test]
    fn s5_ctx_least_violated_retention() {
        let mut lv = LeastViolatedPoint::default();
        assert!(lv.best().is_none());
        assert!(lv.offer(5.0, &[1.0, 2.0]));
        assert!(!lv.offer(7.0, &[9.0, 9.0]));
        // Ties go to the latest.
        assert!(lv.offer(5.0, &[3.0, 4.0]));
        let best = lv.best().unwrap();
        assert_eq!(best.violation, 5.0);
        assert_eq!(&*best.x, &[3.0, 4.0]);
        assert!(!lv.would_take(6.0));
    }

    // -----------------------------------------------------------------------
    // heuristic_context.cpp: restart helpers
    // -----------------------------------------------------------------------

    #[test]
    fn s5_ctx_take_restart_rounds_and_clamps() {
        let problem = S5Problem {
            col_lower: vec![0.0, -2.0, 0.0],
            col_upper: vec![1.0, 2.0, 5.0],
            integrality: vec![
                HighsVarType::kInteger,
                HighsVarType::kContinuous,
                HighsVarType::kInteger,
            ],
            ncol: 3,
            ..Default::default()
        };
        let mut rng = S5Rng::new(1);

        let mut out = Vec::new();
        let mut src = |_rng: &mut S5Rng, out: &mut Vec<f64>| -> bool {
            out.clear();
            out.extend_from_slice(&[0.4, 5.0, -4.2]);
            true
        };
        assert!(take_restart(Some(&mut src), &mut rng, &mut out, &problem));
        assert_eq!(out, vec![0.0, 2.0, 0.0]);

        let mut out = Vec::new();
        let mut src = |_rng: &mut S5Rng, out: &mut Vec<f64>| -> bool {
            out.clear();
            out.extend_from_slice(&[0.4, f64::NAN, -4.2]);
            true
        };
        assert!(!take_restart(Some(&mut src), &mut rng, &mut out, &problem));
        assert!(out.is_empty());

        let mut out = vec![1.0, 2.0, 3.0];
        assert!(!take_restart(None, &mut rng, &mut out, &problem));
        assert!(out.is_empty());
    }

    #[test]
    fn s5_ctx_hybrid_restart_draws_within_bounds() {
        let problem = S5Problem {
            col_lower: vec![0.0, -2.0, -3.0, 0.0],
            col_upper: vec![1.0, 2.0, 3.0, kHighsInf],
            integrality: vec![
                HighsVarType::kInteger,
                HighsVarType::kContinuous,
                HighsVarType::kInteger,
                HighsVarType::kInteger,
            ],
            ncol: 4,
            ..Default::default()
        };

        let mut rng = S5Rng::new(42);
        let mut point = [0.5, 0.5, 0.5, 7.0];
        hybrid_restart(&mut point, &problem, &mut rng);
        let first = point;

        let mut rng = S5Rng::new(42);
        let mut point = [0.5, 0.5, 0.5, 7.0];
        hybrid_restart(&mut point, &problem, &mut rng);
        assert_eq!(point, first);

        // The unbounded column is never touched.
        assert_eq!(point[3], 7.0);

        let lb = [0.0, -2.0, -3.0, 0.0];
        let ub = [1.0, 2.0, 3.0, kHighsInf];
        let integ = [true, false, true, true];
        // Only the columns the coin redrew are checked: the others keep 0.5.
        let orig = [0.5, 0.5, 0.5, 7.0];
        for j in 0..4 {
            if point[j] == orig[j] {
                continue;
            }
            assert!(point[j] >= lb[j] && point[j] <= ub[j]);
            if integ[j] {
                assert_eq!(point[j].fract(), 0.0);
            }
        }
    }

    #[test]
    fn s5_ctx_rebuild_default_draws_nothing() {
        let problem = S5Problem {
            ncol: 3,
            ..Default::default()
        };
        let mut rng = S5Rng::new(7);
        let mut sink = sink();
        let mut out = vec![1.0, 2.0, 3.0];
        assert!(!rebuild_restart_point(
            RestartPoint::kDefault,
            None,
            &mut sink,
            &mut rng,
            &mut out,
            &problem
        ));
        assert!(out.is_empty());
        // Nothing was drawn: the next draw matches a fresh clone's.
        let mut clone = rng.clone();
        assert_eq!(rng.next_u64(), clone.next_u64());
    }
}
