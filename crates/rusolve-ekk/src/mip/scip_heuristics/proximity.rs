//! Translated from SCIP (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! `mip_proximity`: proximity search on the port stack.
//!
//! * M. Fischetti, M. Monaci, "Proximity search for 0-1 mixed-integer convex programming",
//!   J. Heuristics 20 (2014): replace the objective by the Hamming distance to the incumbent,
//!   keep the original objective as the row `obj <= incumbent - delta`, solve, repeat.
//! * SCIP `src/scip/heur_proximity.c` (63417ae): `heurExecProximity` (call gating, node limit
//!   `nodesquot * nnodes + nodesofs - usednodes` capped at `maxnodes`) and `SCIPapplyProximity`
//!   (incumbent and waiting-node gates, the cutoff `minimprove * lb + (1 - minimprove) * ub`,
//!   floored when the objective is integral, `limits/solutions = 1`, restart on success).
//!   Defaults: maxnodes 10000, minimprove 0.02, mingap 0.01, minnodes 1, nodesofs 50,
//!   waitingnodes 100, nodesquot 0.1, binvarquot 0.1, restart on.
//!
//! Tuning (from SCIP `heur_proximity.c`): `lpitersquot` becomes a work limit per sub-MIP, the
//! node quotient a work budget, plus a stall trigger (ours) and the cheap sub-SCIP settings
//! (separation off, no recursive LNS, no restarts). Divergences from SCIP otherwise: `maxlpiters`,
//! `minlpiters` and the node-solved event handler are not ported (a root LP in the sub-MIP is not
//! interrupted); HiGHS presolve stays on, with no restartdfs / inference branching; the
//! infinite-bound adjustment of the proximity sub-problem is not ported (only binaries carry a
//! proximity objective, so no variable with a nonzero cost can run to an infinite bound);
//! `usefinallp`, `useuct` and `uselprows` (default off) are not ported.
//! Everything runs on the main thread at the hook sites; no randomness, so threads=1 is
//! deterministic.

use super::sub_mip_s6::S6SubMipExt;
use super::*;
use crate::lp_data::hconst::HighsBasisStatus;
use crate::mip::mip_solver::HighsMipSolver;
use crate::mip::mip_solver_data::{
    kSolutionSourceProximity, kSolutionSourceTrivialL, kSolutionSourceTrivialP,
    kSolutionSourceTrivialU, kSolutionSourceTrivialZ,
};
use crate::util::sparse_matrix::HighsSparseMatrix;
use crate::work::WorkSnapshot;

/// SCIP heur_proximity.c DEFAULT_MAXNODES
pub(crate) const PROX_MAXNODES: i64 = 10000;
/// SCIP heur_proximity.c DEFAULT_MINIMPROVE
pub(crate) const PROX_MINIMPROVE: f64 = 0.02;
/// SCIP heur_proximity.c DEFAULT_MINGAP
pub(crate) const PROX_MINGAP: f64 = 0.01;
/// SCIP heur_proximity.c DEFAULT_MINNODES
pub(crate) const PROX_MINNODES: i64 = 1;
/// SCIP heur_proximity.c DEFAULT_NODESOFS
pub(crate) const PROX_NODESOFS: i64 = 50;
/// SCIP heur_proximity.c DEFAULT_WAITINGNODES
pub(crate) const PROX_WAITINGNODES: i64 = 100;
/// SCIP heur_proximity.c DEFAULT_NODESQUOT
pub(crate) const PROX_NODESQUOT: f64 = 0.1;
/// SCIP heur_proximity.c DEFAULT_BINVARQUOT
pub(crate) const PROX_BINVARQUOT: f64 = 0.1;

/// Work share of the main search that proximity may spend: SCIP's node quotient `nodesquot` 0.1,
/// in work units (`src/work.rs`) instead of nodes, as for ALNS.
pub(crate) const PROX_WORK_SHARE: f64 = 0.1;
/// Work proximity may spend before the share applies, as a fraction of the root's work (frozen at
/// the first hook): the ALNS rule.
pub(crate) const PROX_WORK_OFFSET_ROOT_FRACTION: f64 = 0.5;
/// Per sub-MIP work cap as a fraction of the root's work: SCIP's `lpitersquot` 0.2 of the root's
/// first LP iterations, `lp/rootiterlim`, in work units.
pub(crate) const PROX_LPITERSQUOT: f64 = 0.2;
/// The incumbent counts as stalled when it has not changed for this share of the main search's
/// work so far. Ours: SCIP calls proximity at every node after `waitingnodes`; models whose
/// incumbent improves on its own are not taxed.
pub(crate) const PROX_STALL_SHARE: f64 = 0.25;
/// Ceiling of the stall share after failures: an incumbent that has not changed for 90% of the
/// search still qualifies.
pub(crate) const PROX_STALL_SHARE_MAX: f64 = 0.9;

#[derive(Clone, Default, Debug)]
pub(crate) struct ProximityState {
    /// nodes already used by proximity in earlier calls (`heurdata->usednodes`)
    pub usednodes: i64,
    /// solution index of the last incumbent processed (`heurdata->lastsolidx`); None = -1
    pub lastsolidx: Option<i64>,
}

/// What proximity keeps across MIP restarts (the solver owns it, `HighsMipSolver::prox`; the
/// per-restart `ProximityState` is rebuilt with the search data): the work ledger, the frozen
/// root cost and the incumbent clock of the stall trigger.
#[derive(Clone, Default, Debug)]
pub(crate) struct ProxLedger {
    /// work spent in proximity sub-MIPs
    pub prox_work: u64,
    pub initialised: bool,
    /// the main search's work at the first hook (what a root costs)
    pub root_work: u64,
    pub offset_work: f64,
    /// index of the incumbent seen at the last hook, and the main work when it first appeared
    pub seen_inc: Option<i64>,
    pub inc_mark: u64,
    /// sub-MIPs in a row that found nothing: each doubles the stall share the next call needs
    pub fails: u32,
    /// `mip_improve_early` is running (no dual bound exists yet)
    pub early_phase: bool,
    /// `mip_primal_stall`: the dual bound when the incumbent clock was last reset
    pub lb_mark: f64,
}

impl ProxLedger {
    /// The work gate: `None` while the ledger is over the allowance, else the work one sub-MIP may
    /// spend (the smaller of the per-call cap and what is left of the allowance).
    pub(crate) fn work_limit(&self, main_work: u64) -> Option<u64> {
        let allowance = PROX_WORK_SHARE * main_work as f64 + self.offset_work;
        let left = allowance - self.prox_work as f64;
        if left < 1.0 {
            return None;
        }
        let cap = (PROX_LPITERSQUOT * self.root_work as f64).max(1.0);
        Some(left.min(cap) as u64)
    }

    /// The stall share the next call needs: `PROX_STALL_SHARE`, doubled per failed call in a row
    /// (up to `PROX_STALL_SHARE_MAX`), so a model that does not profit stops being asked.
    pub(crate) fn stall_share(&self) -> f64 {
        (PROX_STALL_SHARE * f64::from(1u32 << self.fails.min(8))).min(PROX_STALL_SHARE_MAX)
    }

    /// Records the incumbent seen at a hook; true when the incumbent has been unchanged for at
    /// least the stall share of the main work so far.
    pub(crate) fn note_incumbent(&mut self, index: i64, main_work: u64) -> bool {
        if self.seen_inc != Some(index) {
            // the first hook comes after the root: an incumbent met there is as old as the search
            self.inc_mark = if self.seen_inc.is_none() {
                0
            } else {
                main_work
            };
            self.seen_inc = Some(index);
        }
        (main_work.saturating_sub(self.inc_mark)) as f64 >= self.stall_share() * main_work as f64
    }
}

impl ProxLedger {
    /// `mip_primal_stall` trigger: the incumbent has been unchanged for
    /// `PRIMAL_STALL_SHARE` of the main work (doubled per failed call, as `stall_share`) and the dual
    /// bound has closed less than `PRIMAL_STALL_BOUND_CLOSE` of the gap that was open when the incumbent clock started.
    pub(crate) fn note_incumbent_and_bound(
        &mut self,
        index: i64,
        main_work: u64,
        lb: f64,
        ub: f64,
    ) -> bool {
        if self.seen_inc != Some(index) {
            self.inc_mark = if self.seen_inc.is_none() {
                0
            } else {
                main_work
            };
            self.lb_mark = lb;
            self.seen_inc = Some(index);
        }
        let share =
            (PRIMAL_STALL_SHARE * f64::from(1u32 << self.fails.min(8))).min(PROX_STALL_SHARE_MAX);
        let inc_stalled =
            (main_work.saturating_sub(self.inc_mark)) as f64 >= share * main_work as f64;
        let open = (ub - self.lb_mark).max(0.0);
        let bound_stalled = lb - self.lb_mark <= PRIMAL_STALL_BOUND_CLOSE * open;
        inc_stalled && bound_stalled
    }
}

/// `mip_primal_stall`: share of the main work an unchanged incumbent must age before the trigger fires
/// (`PROX_STALL_SHARE` is 0.25; ours is lower because the dual condition guards the cost).
pub(crate) const PRIMAL_STALL_SHARE: f64 = 0.1;
/// `mip_primal_stall`: the dual bound counts as stalled when it closed less than this share of the gap
/// open at the start of the incumbent's age.
pub(crate) const PRIMAL_STALL_BOUND_CLOSE: f64 = 0.01;

/// The incumbent as the pool holds it.
pub(crate) struct Incumbent {
    pub x: Vec<f64>,
    pub obj: f64,
    pub index: i64,
    pub nodenum: i64,
}

/// SCIP: the solution does not come from the trivial heuristic.
pub(crate) fn incumbent_of(ms: &HighsMipSolver) -> Option<Incumbent> {
    let n_cols = mip_model!(ms).num_col_ as usize;
    let best = ms.mipdata_.s6.as_ref()?.pool.best()?;
    if best.x.len() != n_cols {
        return None;
    }
    if let S6Origin::Highs(src) = best.origin
        && (src == kSolutionSourceTrivialL
            || src == kSolutionSourceTrivialP
            || src == kSolutionSourceTrivialU
            || src == kSolutionSourceTrivialZ)
    {
        return None;
    }
    Some(Incumbent {
        x: best.x.clone(),
        obj: best.obj,
        index: best.index,
        nodenum: best.nodenum,
    })
}

/// Unfixed binaries of the presolved model under the global domain (SCIP `SCIPgetNBinVars`).
pub(crate) fn binary_cols(ms: &HighsMipSolver) -> Vec<usize> {
    let d = &ms.mipdata_.domains[0];
    ms.mipdata_
        .integral_cols
        .iter()
        .map(|&col| col as usize)
        .filter(|&j| d.col_lower_[j] == 0.0 && d.col_upper_[j] == 1.0)
        .collect()
}

/// SCIP `SCIPgetGap`: relative gap, `|p - d| / max(|p|, |d|, 1)`, infinite when a bound is.
pub(crate) fn scip_gap(primal: f64, dual: f64) -> f64 {
    if s6_is_eq(primal, dual) {
        return 0.0;
    }
    if s6_is_inf(primal.abs()) || s6_is_inf(dual.abs()) {
        return f64::INFINITY;
    }
    let max = cpp_max(cpp_max(primal.abs(), dual.abs()), 1.0);
    ((primal - dual) / max).abs()
}

/// Lower bound, rounded up when the objective is integral (`SCIPisObjIntegral`).
pub(crate) fn rounded_lower_bound(ms: &HighsMipSolver, lb: f64) -> f64 {
    let of = &ms.mipdata_.objective_function;
    if of.is_integral() && !s6_is_inf(lb.abs()) {
        let s = of.integral_scale();
        s6_feas_ceil(lb * s, ms.mipdata_.feastol) / s
    } else {
        lb
    }
}

/// Appends the row `lo <= sum val_k x_idx_k <= hi` to `lp` (kept column-wise).
pub(crate) fn append_row(
    lp: &mut HighsLp,
    name: &str,
    idx: Vec<HighsInt>,
    val: Vec<f64>,
    lo: f64,
    hi: f64,
) {
    let nz = idx.len() as HighsInt;
    let new_rows = HighsSparseMatrix {
        format_: crate::lp_data::hconst::MatrixFormat::kRowwise,
        num_col_: lp.num_col_,
        num_row_: 1,
        start_: vec![0, nz],
        p_end_: Vec::new(),
        index_: idx,
        value_: val,
    };
    // the matrix's own dimensions drive the insertion; after a MIP restart the presolved model's
    // matrix header can lag behind the LP's
    lp.a_matrix_.num_col_ = lp.num_col_;
    lp.a_matrix_.num_row_ = lp.num_row_;
    lp.a_matrix_.add_rows(&new_rows, None);
    lp.a_matrix_.ensure_colwise();
    lp.row_lower_.push(lo);
    lp.row_upper_.push(hi);
    // a named model needs a name for every row (presolve shrinks the name vector in step)
    if lp.row_names_.len() == lp.num_row_ as usize {
        lp.row_names_.push(name.to_string());
    }
    lp.num_row_ += 1;
}

/// The parent's root basis extended by a basic slack for one appended row.
pub(crate) fn basis_with_extra_row(ms: &HighsMipSolver, num_row: usize) -> HighsBasis {
    let mut basis = ms.mipdata_.firstrootbasis.clone();
    if basis.row_status.len() == num_row {
        basis.row_status.push(HighsBasisStatus::kBasic);
    }
    basis
}

/// SCIP `heurExecProximity` + `SCIPapplyProximity` at one hook site.
pub(crate) fn proximity_step(ms: &mut HighsMipSolver) {
    let primal_stall = ms.options_mip_.mip_primal_stall && !ms.options_mip_.mip_proximity;
    if ms.submip || !(ms.options_mip_.mip_proximity || primal_stall) || ms.terminate() {
        return;
    }
    if ms.mipdata_.nodequeue.empty() || ms.mipdata_.proximity.is_none() {
        return;
    }
    let nvars = mip_model!(ms).num_col_ as usize;
    let bins = binary_cols(ms);
    // do not run heuristic when there are only few binary variables
    if (bins.len() as f64) < PROX_BINVARQUOT * nvars as f64 {
        return;
    }
    // do not run proximity if the problem does not have an objective function anyway
    if ms.mipdata_.objective_function.is_zero() {
        return;
    }
    let mut st = ms.mipdata_.proximity.take().expect("checked above");
    let total_work = WorkSnapshot::now().since(&ms.stats_.work_start).work();
    let main_work = total_work.saturating_sub(ms.prox.prox_work);
    if !ms.prox.initialised {
        ms.prox.initialised = true;
        ms.prox.root_work = main_work;
        ms.prox.offset_work = PROX_WORK_OFFSET_ROOT_FRACTION * main_work as f64;
    }
    // stall trigger: only an incumbent that has stopped improving is worth a sub-MIP
    let best_index = ms
        .mipdata_
        .s6
        .as_ref()
        .and_then(|s6| s6.pool.best())
        .map(|b| b.index);
    let lb = ms.mipdata_.lower_bound;
    let ub = ms.mipdata_.upper_bound;
    let stalled = match best_index {
        Some(index) if primal_stall => ms.prox.note_incumbent_and_bound(index, main_work, lb, ub),
        Some(index) => ms.prox.note_incumbent(index, main_work),
        None => false,
    };
    if !stalled {
        ms.mipdata_.proximity = Some(st);
        return;
    }
    let nnodes_total = ms.mipdata_.num_nodes;
    let mut nnodes = (PROX_NODESQUOT * nnodes_total as f64) as i64 + PROX_NODESOFS;
    nnodes -= st.usednodes;
    nnodes = nnodes.min(PROX_MAXNODES);
    if nnodes >= PROX_MINNODES {
        let mut last_rhs: Option<f64> = None;
        loop {
            let work_now = WorkSnapshot::now().since(&ms.stats_.work_start).work();
            let main_now = work_now.saturating_sub(ms.prox.prox_work);
            let Some(work_limit) = ms.prox.work_limit(main_now) else {
                break;
            };
            let (found, used) =
                proximity_apply(ms, &mut st, &bins, nnodes, work_limit, &mut last_rhs);
            st.usednodes += used;
            nnodes -= used;
            if !(found && !ms.terminate() && nnodes > 0) {
                break;
            }
        }
    }
    ms.mipdata_.proximity = Some(st);
}

/// Each step of the pre-root chain may spend this fraction of the chain's work budget
/// (`mip_improve_early_work`, a multiple of the work spent before the root: setup, presolve,
/// feasibility jump). Ours, measured on our benchmark set.
pub(crate) const EARLY_STEPS_PER_CHAIN: f64 = 5.0;

/// `mip_improve_early`: proximity search (Fischetti-Monaci 2014) as an improvement chain on the
/// feasibility-jump incumbent, before the root LP. Feasibility jump hands over a feasible point
/// that is far from optimal and the tree only reaches its first hook after the root; here
/// every step replaces the incumbent by a Hamming-nearby point under a
/// tightened objective cutoff, so the chain walks the objective down. It stops at the first
/// sub-MIP that does not improve, or when the work budget is spent. No stall gate (the incumbent is
/// brand new), no root basis (none exists yet), no tree separation in the sub-MIP (`cheap`).
pub(crate) fn improve_early(ms: &mut HighsMipSolver) {
    if ms.submip || ms.terminate() || ms.mipdata_.proximity.is_none() {
        return;
    }
    if ms.mipdata_.objective_function.is_zero() {
        return;
    }
    let bins = binary_cols(ms);
    if (bins.len() as f64) < PROX_BINVARQUOT * mip_model!(ms).num_col_ as f64 {
        return;
    }
    if incumbent_of(ms).is_none() {
        return;
    }
    let mut st = ms.mipdata_.proximity.take().expect("checked above");
    let pre_work = WorkSnapshot::now().since(&ms.stats_.work_start).work();
    let tm = ms.options_mip_.mip_improve_early_work;
    let call_cap = ((tm / EARLY_STEPS_PER_CHAIN * pre_work as f64) as u64).max(1);
    let total_cap = (tm * pre_work as f64) as u64;
    let mut last_rhs: Option<f64> = None;
    let mut spent: u64 = 0;
    ms.prox.early_phase = true;
    while spent < total_cap && !ms.terminate() {
        let limit = call_cap.min(total_cap - spent);
        let before = ms.prox.prox_work;
        let (found, used) =
            proximity_apply(ms, &mut st, &bins, PROX_MAXNODES, limit, &mut last_rhs);
        st.usednodes += used;
        spent += (ms.prox.prox_work - before).max(1);
        if !found {
            break;
        }
    }
    ms.prox.early_phase = false;
    ms.mipdata_.proximity = Some(st);
}

/// `SCIPapplyProximity`: one proximity sub-MIP. Returns whether the incumbent improved and the
/// sub-MIP nodes used (0 when the call did not run).
fn proximity_apply(
    ms: &mut HighsMipSolver,
    st: &mut ProximityState,
    bins: &[usize],
    nnodes: i64,
    work_limit: u64,
    last_rhs: &mut Option<f64>,
) -> (bool, i64) {
    // only call the heuristic if we have an incumbent (not a trivial one)
    let Some(inc) = incumbent_of(ms) else {
        return (false, 0);
    };
    if st.lastsolidx == Some(inc.index) {
        return (false, 0);
    }
    let num_nodes = ms.mipdata_.num_nodes;
    // waitingnodes: nodes to wait before a new incumbent is processed
    if num_nodes > 1 && num_nodes - inc.nodenum < PROX_WAITINGNODES {
        return (false, 0);
    }
    let feastol = ms.mipdata_.feastol;
    let bestobj = inc.obj;
    let lowerbound = rounded_lower_bound(ms, ms.mipdata_.lower_bound);
    // primal and dual bound already close together
    // before the root LP there is no dual bound (-inf; `SCIPrelDiff` would give NaN): the gap
    // gate only applies once a bound exists
    let no_bound = ms.prox.early_phase && s6_is_inf(lowerbound.abs());
    if !no_bound
        && (s6_feas_le(bestobj, lowerbound, feastol)
            || scip_gap(bestobj, lowerbound) <= PROX_MINGAP)
    {
        return (false, 0);
    }
    // minimum improvement in terms of the distance between incumbent objective and lower bound
    let mut objcutoff = if s6_is_inf(lowerbound.abs()) {
        if s6_feas_zero(bestobj, feastol) {
            bestobj - 1.0
        } else {
            (1.0 - PROX_MINIMPROVE) * bestobj
        }
    } else {
        PROX_MINIMPROVE * lowerbound + (1.0 - PROX_MINIMPROVE) * bestobj
    };
    // integrality of the objective strengthens the cutoff
    let of = &ms.mipdata_.objective_function;
    if of.is_integral() {
        let s = of.integral_scale();
        objcutoff = s6_feas_floor(objcutoff * s, feastol) / s;
    }
    if s6_feas_lt(objcutoff, lowerbound, feastol) {
        objcutoff = lowerbound;
    }
    // the right hand side of the objective constraint does not change: not successful before
    if let Some(r) = *last_rhs
        && s6_feas_le(r, objcutoff, feastol)
        && s6_feas_ge(r, objcutoff, feastol)
    {
        return (false, 0);
    }
    *last_rhs = Some(objcutoff);
    st.lastsolidx = Some(inc.index);

    // sub-MIP: proximity objective over the binaries, the original objective as a cutoff row
    let model = mip_model!(ms);
    let n = model.num_col_ as usize;
    let mut lp: HighsLp = model.clone();
    let mut oidx: Vec<HighsInt> = Vec::new();
    let mut oval: Vec<f64> = Vec::new();
    for j in 0..n {
        let cost = lp.col_cost_[j];
        // SCIPisFeasZero on the objective coefficient
        if cost.abs() > feastol {
            oidx.push(j as HighsInt);
            oval.push(cost);
        }
        lp.col_cost_[j] = 0.0;
    }
    for &j in bins {
        lp.col_cost_[j] = if inc.x[j] < 0.5 { 1.0 } else { -1.0 };
    }
    let num_row = lp.num_row_ as usize;
    append_row(
        &mut lp,
        "prox_objcutoff",
        oidx,
        oval,
        -f64::INFINITY,
        objcutoff,
    );
    let basis = basis_with_extra_row(ms, num_row);
    let d = &ms.mipdata_.domains[0];
    let (cl, cu) = (d.col_lower_.clone(), d.col_upper_.clone());

    let old_upper = ms.mipdata_.upper_bound;
    let w0 = WorkSnapshot::now();
    let (_ok, subnodes, _status) = detach!(ms.mipdata_.heuristics, |h| h.s6_solve_sub_mip_ext(
        ms,
        0,
        lp,
        basis,
        0.0,
        cl,
        cu,
        kHighsIInf,
        nnodes as HighsInt,
        kHighsIInf,
        f64::INFINITY,
        kSolutionSourceProximity,
        S6SubMipExt {
            own_objective: true,
            // limits/solutions = 1
            max_improving_sols: Some(1),
            work_limit: Some(work_limit),
            cheap: true,
        },
    ));
    let work = WorkSnapshot::now().since(&w0).work();
    ms.prox.prox_work += work;
    let improved = ms.mipdata_.upper_bound < old_upper;
    ms.prox.fails = if improved { 0 } else { ms.prox.fails + 1 };
    let s = &mut ms.stats_;
    s.prox_calls += 1;
    s.prox_nodes += subnodes;
    s.prox_work += work;
    if improved {
        s.prox_found += 1;
    }
    crate::io::highs_io::highs_log_dev(
        &ms.options_mip_.log_options,
        HighsLogType::kInfo,
        &format!(
            "Proximity cutoff {:.6} nodes {} work {} improved {}\n",
            objcutoff, subnodes, work, improved
        ),
    );
    (improved, subnodes)
}

/// The hook: both Fischetti heuristics, each behind its own option.
pub(crate) fn fischetti_step(ms: &mut HighsMipSolver) {
    if ms.options_mip_.mip_proximity || ms.options_mip_.mip_primal_stall {
        proximity_step(ms);
    }
    if ms.options_mip_.mip_local_branching {
        super::local_branching::local_branching_step(ms);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gap_matches_scip() {
        assert_eq!(scip_gap(5.0, 5.0), 0.0);
        assert!((scip_gap(10.0, 5.0) - 0.5).abs() < 1e-12);
        // denominator is at least 1
        assert!((scip_gap(0.5, 0.0) - 0.5).abs() < 1e-12);
        assert!(scip_gap(1.0, f64::NEG_INFINITY).is_infinite());
    }

    #[test]
    fn work_gate_share_offset_and_cap() {
        let mut led = ProxLedger {
            root_work: 1000,
            offset_work: 500.0,
            ..Default::default()
        };
        // allowance 0.1 * 2000 + 500 = 700, per-call cap 0.2 * 1000 = 200
        assert_eq!(led.work_limit(2000), Some(200));
        led.prox_work = 600;
        // 100 left: below the cap
        assert_eq!(led.work_limit(2000), Some(100));
        led.prox_work = 700;
        assert_eq!(led.work_limit(2000), None);
        // the main search keeps running, so the allowance grows with it
        assert_eq!(led.work_limit(4000), Some(200));
    }

    #[test]
    fn stall_trigger_waits_for_an_old_incumbent() {
        let mut led = ProxLedger::default();
        // first hook: the incumbent is as old as the search
        assert!(led.note_incumbent(7, 1000));
        assert!(led.note_incumbent(7, 2000));
        // a new incumbent restarts the clock
        assert!(!led.note_incumbent(8, 2000));
        assert!(!led.note_incumbent(8, 2400));
        // a quarter of the work so far (2600 * 0.25 = 650) without a change
        assert!(!led.note_incumbent(8, 2600));
        assert!(led.note_incumbent(8, 2700));
        // each failed call in a row doubles the share: 0.5 then 0.9
        led.fails = 1;
        assert!(!led.note_incumbent(8, 2700));
        assert!(led.note_incumbent(8, 4100));
        led.fails = 5;
        assert_eq!(led.stall_share(), PROX_STALL_SHARE_MAX);
    }

    #[test]
    fn primal_stall_needs_an_old_incumbent_and_a_stalled_bound() {
        let mut led = ProxLedger::default();
        // first hook: incumbent as old as the search, bound closed nothing since: fires
        assert!(led.note_incumbent_and_bound(1, 1000, 10.0, 20.0));
        // a bound that moved 20% of the open gap since the incumbent appeared: no
        assert!(!led.note_incumbent_and_bound(1, 2000, 12.0, 20.0));
        // a new incumbent resets the clock
        assert!(!led.note_incumbent_and_bound(2, 2100, 12.0, 18.0));
        // old again, bound flat: fires
        assert!(led.note_incumbent_and_bound(2, 3000, 12.0, 18.0));
    }

    #[test]
    fn appended_row_lands_in_the_matrix() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        lp.row_lower_ = vec![0.0];
        lp.row_upper_ = vec![5.0];
        append_row(&mut lp, "r", vec![0, 1], vec![2.0, -3.0], -1.0, 4.0);
        assert_eq!(lp.num_row_, 2);
        assert_eq!(lp.a_matrix_.num_row_, 2);
        assert!(lp.a_matrix_.is_colwise());
        assert_eq!(lp.a_matrix_.start_, vec![0, 2, 4]);
        assert_eq!(lp.a_matrix_.index_, vec![0, 1, 0, 1]);
        assert_eq!(lp.a_matrix_.value_, vec![1.0, 2.0, 1.0, -3.0]);
        assert_eq!(lp.row_upper_, vec![5.0, 4.0]);
    }

    /// A multi-dimensional knapsack: min -value.x, 3 rows, x binary.
    fn knapsack(n: usize) -> HighsLp {
        use crate::lp_data::hconst::{HighsVarType, MatrixFormat};
        let m = 3;
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as f64
        };
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = m as HighsInt;
        let mut weights = vec![vec![0.0; n]; m];
        for j in 0..n {
            let mut sum = 0.0;
            for row in weights.iter_mut() {
                row[j] = 10.0 + (next() % 40.0).floor();
                sum += row[j];
            }
            lp.col_cost_
                .push(-(sum / m as f64 + (next() % 10.0).floor()).floor());
        }
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf; m];
        lp.row_upper_ = weights
            .iter()
            .map(|r| (r.iter().sum::<f64>() / 3.0).floor())
            .collect();
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = m as HighsInt;
        lp.a_matrix_.start_ = (0..=n).map(|j| (j * m) as HighsInt).collect();
        for j in 0..n {
            for (i, row) in weights.iter().enumerate() {
                lp.a_matrix_.index_.push(i as HighsInt);
                lp.a_matrix_.value_.push(row[j]);
            }
        }
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp
    }

    struct Run {
        obj: f64,
        nodes: i64,
        incumbent: Vec<f64>,
        prox_state: bool,
        lb_state: bool,
        calls: (i64, i64),
    }

    fn solve(prox: bool, lb: bool, submip: bool, n: usize) -> Run {
        std::thread::spawn(move || {
            crate::parallel::highs_parallel::initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.mip_proximity = prox;
            options.mip_local_branching = lb;
            let lp = knapsack(n);
            let mut callback = HighsCallback::default();
            let solution = HighsSolution::default();
            let level = HighsInt::from(submip);
            let mut ms =
                HighsMipSolver::new(&mut callback, &options, &lp, &solution, submip, level);
            ms.run();
            let r = Run {
                obj: ms.mipdata_.upper_bound,
                nodes: ms.mipdata_.num_nodes,
                incumbent: ms.mipdata_.incumbent.clone(),
                prox_state: ms.mipdata_.proximity.is_some(),
                lb_state: ms.mipdata_.local_branching.is_some(),
                calls: (ms.stats_.prox_calls, ms.stats_.lb_calls),
            };
            drop(ms);
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            r
        })
        .join()
        .unwrap()
    }

    const N: usize = 48;

    #[test]
    fn off_creates_nothing_and_in_a_submip_too() {
        let off = solve(false, false, false, N);
        assert!(!off.prox_state && !off.lb_state);
        assert_eq!(off.calls, (0, 0));
        let sub = solve(true, true, true, N);
        assert!(!sub.prox_state && !sub.lb_state);
        assert_eq!(sub.calls, (0, 0));
    }

    #[test]
    fn on_keeps_the_optimum_and_is_deterministic() {
        let off = solve(false, false, false, N);
        let a = solve(true, true, false, N);
        let b = solve(true, true, false, N);
        assert!(a.prox_state && a.lb_state);
        assert!(
            (off.obj - a.obj).abs() < 1e-9,
            "off {} on {}",
            off.obj,
            a.obj
        );
        assert_eq!(a.nodes, b.nodes);
        assert_eq!(a.calls, b.calls);
        assert_eq!(a.obj.to_bits(), b.obj.to_bits());
        assert_eq!(a.incumbent, b.incumbent);
    }

    /// Hand the heuristic the all-zero point as incumbent on a solved model: the proximity
    /// sub-MIP (distance objective, cutoff row) runs and finishes.
    #[test]
    fn proximity_runs_from_a_weak_incumbent() {
        std::thread::spawn(|| {
            crate::parallel::highs_parallel::initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.mip_proximity = true;
            let lp = knapsack(N);
            let mut callback = HighsCallback::default();
            let solution = HighsSolution::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let (lo, up) = (
                mip_model!(ms).col_lower_.clone(),
                mip_model!(ms).col_upper_.clone(),
            );
            ms.mipdata_.domains[0].col_lower_ = lo;
            ms.mipdata_.domains[0].col_upper_ = up;
            let n = mip_model!(ms).num_col_ as usize;
            {
                let s6 = ms.mipdata_.s6.as_mut().expect("s6 state");
                s6.pool.sols.clear();
                s6.pool.sols.push(S6Sol {
                    index: 1_000_000,
                    obj: 0.0,
                    x: vec![0.0; n],
                    origin: S6Origin::Oneopt,
                    nodenum: 0,
                });
            }
            // Move the upper bound the way HiGHS does (setting upper_bound then
            // updatePrimalDualIntegral), so the P-D integral's previous bound matches.
            let prev_upper_bound = ms.mipdata_.upper_bound;
            ms.mipdata_.upper_bound = 0.0;
            let lower_bound = ms.mipdata_.lower_bound;
            ms.update_primal_dual_integral(
                lower_bound,
                lower_bound,
                prev_upper_bound,
                0.0,
                true,
                true,
            );
            let bins = binary_cols(&ms);
            assert!(!bins.is_empty() && bins.len() <= n);
            let mut st = ProximityState::default();
            let mut last = None;
            let (_found, _used) = proximity_apply(&mut ms, &mut st, &bins, 50, u64::MAX, &mut last);
            assert_eq!(ms.stats_.prox_calls, 1);
            assert_eq!(st.lastsolidx, Some(1_000_000));
            // the cutoff is 2% of the way to the lower bound, floored (integral objective)
            let lb = ms.mipdata_.lower_bound;
            assert!(lb < 0.0);
            assert_eq!(last, Some((0.02 * lb).floor()));
            // the same incumbent is not processed twice
            let again = proximity_apply(&mut ms, &mut st, &bins, 50, u64::MAX, &mut last);
            assert_eq!(again, (false, 0));
            assert_eq!(ms.stats_.prox_calls, 1);
            drop(ms);
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }
}
