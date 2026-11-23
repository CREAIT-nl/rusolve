//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! prop_engine.cpp: `PropEngine` constructor (lines 64-99), reset, fix,
//! refix, shift_domain, tighten_lb/ub, seed_worklist, propagate,
//! backtrack_to, and the marks (lines 194-708).

use super::{ActUndo, IndexedMinHeap, PqUndo, PropEngine, PropResult};
use crate::mip::fpr_local_mip::{
    Deadline, VarState, kPropagateBudgetPerNnz, kPropagateDeadlinePollWork,
};
use crate::prelude::*;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// mip-heuristics (36f7d9d): deadline.h::Deadline::expired (Q4: timer is always None, so this is always false)
fn deadline_expired(d: &Deadline) -> bool {
    match d.timer {
        Some(never) => match never {},
        None => false,
    }
}

impl<'a> PropEngine<'a> {
    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::PropEngine (lines 64-99)
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        ncol: HighsInt,
        nrow: HighsInt,
        ar_start: &'a [HighsInt],
        ar_index: &'a [HighsInt],
        ar_value: &'a [f64],
        csc_start: &'a [HighsInt],
        csc_row_p: &'a [HighsInt],
        csc_val_p: &'a [f64],
        col_lb: &'a [f64],
        col_ub: &'a [f64],
        row_lo: &'a [f64],
        row_hi: &'a [f64],
        integrality: &'a [HighsVarType],
        feastol: f64,
    ) -> Self {
        let nnz = ar_start[nrow as usize] as usize;
        let vs: Vec<VarState> = vec![VarState::default(); ncol as usize];
        let solution: Vec<f64> = vec![0.0; ncol as usize];
        let prop_in_wl: Vec<u8> = vec![0; nrow as usize];

        let mut eng = PropEngine {
            ncol_: ncol,
            nrow_: nrow,
            ar_start_: ar_start,
            ar_index_: ar_index,
            ar_value_: ar_value,
            col_start_: csc_start,
            col_row_: csc_row_p,
            col_val_: csc_val_p,
            col_lb_: col_lb,
            col_ub_: col_ub,
            row_lo_: row_lo,
            row_hi_: row_hi,
            integrality_: integrality,
            feastol_: feastol,
            nnz_: nnz,
            vs_: vs,
            solution_: solution,
            vs_undo_: Vec::new(),
            sol_undo_: Vec::new(),
            prop_worklist_: Vec::new(),
            prop_in_wl_: prop_in_wl,
            prop_work_: 0,
            deadline_: Deadline::default(),
            min_activity_: Vec::new(),
            max_activity_: Vec::new(),
            min_finite_: Vec::new(),
            max_finite_: Vec::new(),
            min_ninf_: Vec::new(),
            max_ninf_: Vec::new(),
            act_undo_: Vec::new(),
            pq_active_: false,
            domain_pq_: IndexedMinHeap::default(),
            pq_undo_: Vec::new(),
        };

        // DFS + propagation can push many more than ncol undo entries (every
        // bound tightening pushes one). Over-reserve to reduce reallocation
        // churn in the hot path; the memory cost is bounded and transient.
        eng.vs_undo_.reserve(4 * ncol as usize);
        eng.sol_undo_.reserve(4 * ncol as usize);
        eng.prop_worklist_.reserve(nrow as usize);
        // Preallocate domain-PQ storage (vector-backed indexed heap replaces
        // the prior std::set; reserve here so subsequent reset() cycles reuse
        // this allocation).
        eng.domain_pq_.reserve(ncol);
        eng.reset();
        eng
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::reset (lines 194-235)
    pub(crate) fn reset(&mut self) {
        // Reset variable states to global bounds and solution values to 0. We
        // re-derive from col_lb_/col_ub_ (the problem-data pointers, unchanged
        // across reset) rather than replaying vs_undo_ because the previous
        // attempt may have finished on a non-backtracked path (e.g. found a
        // complete assignment) and the undo stacks do not fully unwind in that
        // case.
        for j in 0..self.ncol_ {
            self.vs_[j as usize].lb = self.col_lb_[j as usize];
            self.vs_[j as usize].ub = self.col_ub_[j as usize];
            self.vs_[j as usize].val = 0.0;
            self.vs_[j as usize].fixed = false;
            self.solution_[j as usize] = 0.0;
        }
        self.vs_undo_.clear();
        self.sol_undo_.clear();
        // Clear worklist membership flags that might still be set from the
        // previous attempt's propagation. Iterating prop_worklist_ is cheapest
        // when the worklist is non-empty; if it's already empty this is a
        // no-op.
        for i in 0..self.prop_worklist_.len() {
            let wi = self.prop_worklist_[i];
            self.prop_in_wl_[wi as usize] = 0;
        }
        self.prop_worklist_.clear();
        // Activity tracking is opt-in (caller must invoke init_activities() if
        // needed). Dropping the size signals "not initialized" to downstream
        // code via activities_initialized(). Capacity is preserved so the next
        // init_activities() call does not reallocate.
        self.min_activity_.clear();
        self.max_activity_.clear();
        self.act_undo_.clear();
        self.pq_active_ = false;
        self.domain_pq_.clear();
        self.pq_undo_.clear();
        // Zero the accumulated propagation effort so the Phase 1-2 DFS gate
        // -- the per-call slice, measured as a delta from each call's own
        // starting effort -- starts each attempt fresh. Without this, a
        // second fpr_attempt call on the same PropEngine would see effort
        // already at or above the previous attempt's total and exit the DFS
        // loop immediately.
        self.prop_work_ = 0;
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::fix (lines 237-254)
    pub(crate) fn fix(&mut self, j: HighsInt, mut value: f64) -> bool {
        if value < self.vs_[j as usize].lb - self.feastol_
            || value > self.vs_[j as usize].ub + self.feastol_
        {
            return false;
        }
        value = cpp_max(
            self.vs_[j as usize].lb,
            cpp_min(self.vs_[j as usize].ub, value),
        );
        if self.is_int(j) {
            value = value.round();
        }
        let old_vs = self.vs_[j as usize];
        self.vs_undo_.push((j, old_vs));
        self.sol_undo_.push((j, self.solution_[j as usize]));
        self.vs_[j as usize].fixed = true;
        self.vs_[j as usize].val = value;
        self.solution_[j as usize] = value;
        self.update_activities(j, &old_vs);
        self.pq_notify(j, &old_vs);
        true
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::refix (lines 256-280)
    pub(crate) fn refix(&mut self, j: HighsInt, mut value: f64) -> bool {
        // Structural bounds, not the current (possibly narrower or, for a
        // decision-fixed column, possibly stale-wide) vs_[j] range -- see the
        // header comment on why refix() deliberately overrides rather than
        // refines the current domain.
        if value < self.col_lb_[j as usize] - self.feastol_
            || value > self.col_ub_[j as usize] + self.feastol_
        {
            return false;
        }
        value = cpp_max(
            self.col_lb_[j as usize],
            cpp_min(self.col_ub_[j as usize], value),
        );
        if self.is_int(j) {
            value = value.round();
        }
        let old_vs = self.vs_[j as usize];
        self.vs_undo_.push((j, old_vs));
        self.sol_undo_.push((j, self.solution_[j as usize]));
        // Narrow to a coherent singleton, unlike fix() -- see header comment.
        self.vs_[j as usize].lb = value;
        self.vs_[j as usize].ub = value;
        self.vs_[j as usize].fixed = true;
        self.vs_[j as usize].val = value;
        self.solution_[j as usize] = value;
        self.update_activities(j, &old_vs);
        self.pq_notify(j, &old_vs);
        true
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::shift_domain (lines 282-322)
    pub(crate) fn shift_domain(&mut self, j: HighsInt, delta: f64) -> bool {
        let old_vs = self.vs_[j as usize];
        // A fixed column's effective interval is the singleton at its value --
        // `fix()` deliberately leaves `[lb, ub]` wide, so reading the bounds
        // here would translate an interval the search does not actually hold.
        let cur_lb = if old_vs.fixed { old_vs.val } else { old_vs.lb };
        let cur_ub = if old_vs.fixed { old_vs.val } else { old_vs.ub };
        if cur_lb <= -kHighsInf || cur_ub >= kHighsInf {
            return false;
        }
        let mut new_lb = cur_lb + delta;
        let mut new_ub = cur_ub + delta;
        if new_lb < self.col_lb_[j as usize] - self.feastol_
            || new_ub > self.col_ub_[j as usize] + self.feastol_
        {
            return false;
        }
        // Numerical hygiene only -- the caller's clip already put the interval
        // inside the structural bounds, and clamping can only ever shrink it,
        // never widen it past what the check above accepted.
        new_lb = cpp_max(new_lb, self.col_lb_[j as usize]);
        new_ub = cpp_min(new_ub, self.col_ub_[j as usize]);
        if new_lb > new_ub + self.feastol_ {
            return false;
        }

        self.vs_undo_.push((j, old_vs));
        self.sol_undo_.push((j, self.solution_[j as usize]));
        self.vs_[j as usize].lb = new_lb;
        self.vs_[j as usize].ub = new_ub;
        if old_vs.fixed || new_ub - new_lb < self.feastol_ {
            let mut val = (new_lb + new_ub) * 0.5;
            if self.is_int(j) {
                val = val.round();
            }
            self.vs_[j as usize].fixed = true;
            self.vs_[j as usize].val = val;
            self.solution_[j as usize] = val;
        }
        self.update_activities(j, &old_vs);
        self.pq_notify(j, &old_vs);
        true
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::tighten_lb (lines 324-353)
    pub(crate) fn tighten_lb(&mut self, j: HighsInt, mut new_lb: f64) -> bool {
        if self.is_int(j) {
            new_lb = (new_lb - self.feastol_).ceil();
        }
        new_lb = cpp_max(new_lb, self.col_lb_[j as usize]);
        if new_lb > self.vs_[j as usize].ub + self.feastol_ {
            return false;
        }
        if new_lb <= self.vs_[j as usize].lb + self.feastol_ {
            return true; // no tightening
        }
        let old_vs = self.vs_[j as usize];
        self.vs_undo_.push((j, old_vs));
        self.sol_undo_.push((j, self.solution_[j as usize]));
        self.vs_[j as usize].lb = new_lb;
        // Auto-fix if domain becomes singleton
        if !self.vs_[j as usize].fixed
            && self.vs_[j as usize].ub - self.vs_[j as usize].lb < self.feastol_
        {
            let mut val = (self.vs_[j as usize].lb + self.vs_[j as usize].ub) * 0.5;
            if self.is_int(j) {
                val = val.round();
            }
            self.vs_[j as usize].fixed = true;
            self.vs_[j as usize].val = val;
            self.solution_[j as usize] = val;
        }
        self.update_activities(j, &old_vs);
        self.pq_notify(j, &old_vs);
        self.seed_worklist(j);
        true
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::tighten_ub (lines 355-384)
    pub(crate) fn tighten_ub(&mut self, j: HighsInt, mut new_ub: f64) -> bool {
        if self.is_int(j) {
            new_ub = (new_ub + self.feastol_).floor();
        }
        new_ub = cpp_min(new_ub, self.col_ub_[j as usize]);
        if new_ub < self.vs_[j as usize].lb - self.feastol_ {
            return false;
        }
        if new_ub >= self.vs_[j as usize].ub - self.feastol_ {
            return true; // no tightening
        }
        let old_vs = self.vs_[j as usize];
        self.vs_undo_.push((j, old_vs));
        self.sol_undo_.push((j, self.solution_[j as usize]));
        self.vs_[j as usize].ub = new_ub;
        // Auto-fix if domain becomes singleton
        if !self.vs_[j as usize].fixed
            && self.vs_[j as usize].ub - self.vs_[j as usize].lb < self.feastol_
        {
            let mut val = (self.vs_[j as usize].lb + self.vs_[j as usize].ub) * 0.5;
            if self.is_int(j) {
                val = val.round();
            }
            self.vs_[j as usize].fixed = true;
            self.vs_[j as usize].val = val;
            self.solution_[j as usize] = val;
        }
        self.update_activities(j, &old_vs);
        self.pq_notify(j, &old_vs);
        self.seed_worklist(j);
        true
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::seed_worklist (lines 386-398)
    pub(crate) fn seed_worklist(&mut self, j: HighsInt) {
        let kbeg = self.col_start_[j as usize];
        let kend = self.col_start_[(j + 1) as usize];
        for p in kbeg..kend {
            let i = self.col_row_[p as usize];
            if self.prop_in_wl_[i as usize] == 0 {
                self.prop_in_wl_[i as usize] = 1;
                self.prop_worklist_.push(i);
            }
        }
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::propagate (lines 405-637)
    pub(crate) fn propagate(&mut self, fixed_var: HighsInt) -> PropResult {
        // Hoist raw pointers once so the tight loop doesn't re-load member
        // fields on every iteration. All of these arrays are read-only (the
        // problem data) except vs_/solution_, which we mutate via vs_data and
        // also push undo entries for.
        let ar_start = self.ar_start_;
        let ar_index = self.ar_index_;
        let ar_value = self.ar_value_;
        let col_lb = self.col_lb_;
        let col_ub = self.col_ub_;
        let feastol = self.feastol_;

        if fixed_var >= 0 {
            // Seed only the rows containing the just-fixed variable (AC-3).
            // Clear stale prop_in_wl_ flags before clearing the worklist.
            for wi in 0..self.prop_worklist_.len() {
                let w = self.prop_worklist_[wi];
                self.prop_in_wl_[w as usize] = 0;
            }
            self.prop_worklist_.clear();
            let cbeg = self.col_start_[fixed_var as usize];
            let cend = self.col_start_[(fixed_var + 1) as usize];
            for p in cbeg..cend {
                let i = self.col_row_[p as usize];
                if self.prop_in_wl_[i as usize] == 0 {
                    self.prop_in_wl_[i as usize] = 1;
                    self.prop_worklist_.push(i);
                }
            }
        }
        // When fixed_var == -1, assume prop_worklist is already seeded by caller.

        let mut prop_work: usize = 0;
        let prop_budget = kPropagateBudgetPerNnz * self.nnz_;
        // Cadence counter for the wall-clock poll (issue #151).  SIZE_MAX when
        // no deadline is armed, so the comparison can never fire, no clock is
        // ever read, and an un-armed engine pays one integer compare per row.
        let mut next_deadline_poll: usize = if self.deadline_.timer.is_some() {
            kPropagateDeadlinePollWork
        } else {
            usize::MAX
        };

        while let Some(i) = self.prop_worklist_.pop() {
            self.prop_in_wl_[i as usize] = 0;
            let kbeg = ar_start[i as usize];
            let kend = ar_start[(i + 1) as usize];
            prop_work += (kend - kbeg) as usize;
            if prop_work > prop_budget {
                return self.truncate_fixpoint(prop_work, PropResult::kBudgetExhausted);
            }
            // The clock is polled here and nowhere else in this function: at
            // the same point as the budget check, so both truncations leave the
            // identical engine state.  Deliberately *not* on entry -- an entry
            // poll would put one clock read on every DFS node however small its
            // fixpoint, and it would make `propagate` able to return having
            // charged nothing, which is the "every begin runs at least one
            // propagate round (>0 ops)" premise `fpr.cpp`'s no-progress guard
            // reasons from.  A caller that wants to stop before propagating has
            // its own poll (`fpr_core.cpp`, `repair_search.cpp`).
            if prop_work >= next_deadline_poll {
                if deadline_expired(&self.deadline_) {
                    return self.truncate_fixpoint(prop_work, PropResult::kDeadlineExpired);
                }
                next_deadline_poll = prop_work + kPropagateDeadlinePollWork;
            }

            // Pass 1: compute row aggregates (fixed_sum, min_act, max_act) and
            // count unfixed entries. Pull the coefficient, index, and VarState
            // into locals to keep the inner loop read-dominated.
            let mut fixed_sum: f64 = 0.0;
            let mut min_act: f64 = 0.0;
            let mut max_act: f64 = 0.0;
            let mut num_unfixed: HighsInt = 0;

            for k in kbeg..kend {
                let j = ar_index[k as usize];
                let a = ar_value[k as usize];
                let vj = self.vs_[j as usize];
                if vj.fixed {
                    fixed_sum += a * vj.val;
                } else {
                    num_unfixed += 1;
                    if a > 0.0 {
                        min_act += a * vj.lb;
                        max_act += a * vj.ub;
                    } else {
                        min_act += a * vj.ub;
                        max_act += a * vj.lb;
                    }
                }
            }

            if num_unfixed == 0 {
                continue;
            }

            let has_upper = self.row_hi_[i as usize] < kHighsInf;
            let has_lower = self.row_lo_[i as usize] > -kHighsInf;
            let row_hi_i = self.row_hi_[i as usize];
            let row_lo_i = self.row_lo_[i as usize];

            // Pass 2: tighten bounds on each unfixed variable. We re-access
            // vs_data[j] here because the first pass only reads; we need a
            // mutable reference for bound tightening.
            for k in kbeg..kend {
                let j = ar_index[k as usize];
                if self.vs_[j as usize].fixed {
                    continue;
                }
                let a = ar_value[k as usize];
                if a.abs() < 1e-15 {
                    continue;
                }

                let old_lb = self.vs_[j as usize].lb;
                let old_ub = self.vs_[j as usize].ub;
                let is_integer = self.integrality_[j as usize] != HighsVarType::kContinuous;

                let min_others: f64;
                let max_others: f64;
                if a > 0.0 {
                    min_others = min_act - (a * old_lb);
                    max_others = max_act - (a * old_ub);
                } else {
                    min_others = min_act - (a * old_ub);
                    max_others = max_act - (a * old_lb);
                }

                let mut new_lb = old_lb;
                let mut new_ub = old_ub;

                if has_upper {
                    let bound = row_hi_i - fixed_sum - min_others;
                    if a > 0.0 {
                        new_ub = cpp_min(new_ub, bound / a);
                    } else {
                        new_lb = cpp_max(new_lb, bound / a);
                    }
                }

                if has_lower {
                    let bound = row_lo_i - fixed_sum - max_others;
                    if a > 0.0 {
                        new_lb = cpp_max(new_lb, bound / a);
                    } else {
                        new_ub = cpp_min(new_ub, bound / a);
                    }
                }

                if is_integer {
                    new_lb = (new_lb - feastol).ceil();
                    new_ub = (new_ub + feastol).floor();
                }

                new_lb = cpp_max(new_lb, col_lb[j as usize]);
                new_ub = cpp_min(new_ub, col_ub[j as usize]);

                if new_lb > new_ub + feastol {
                    self.prop_work_ += prop_work;
                    // Clear stale worklist markers on infeasibility too
                    for wi in 0..self.prop_worklist_.len() {
                        let w = self.prop_worklist_[wi];
                        self.prop_in_wl_[w as usize] = 0;
                    }
                    self.prop_worklist_.clear();
                    return PropResult::kInfeasible;
                }

                let tighter_lb = new_lb > old_lb + feastol;
                let tighter_ub = new_ub < old_ub - feastol;
                let mut changed = false;
                let pre_change_vs = self.vs_[j as usize];

                if tighter_lb || tighter_ub {
                    self.vs_undo_.push((j, pre_change_vs));
                    self.sol_undo_.push((j, self.solution_[j as usize]));
                    if tighter_lb {
                        self.vs_[j as usize].lb = new_lb;
                    }
                    if tighter_ub {
                        self.vs_[j as usize].ub = new_ub;
                    }
                    changed = true;
                }

                if !self.vs_[j as usize].fixed
                    && self.vs_[j as usize].ub - self.vs_[j as usize].lb < feastol
                {
                    if !changed {
                        self.vs_undo_.push((j, pre_change_vs));
                        self.sol_undo_.push((j, self.solution_[j as usize]));
                    }
                    let mut val = (self.vs_[j as usize].lb + self.vs_[j as usize].ub) * 0.5;
                    if is_integer {
                        val = val.round();
                    }
                    self.vs_[j as usize].fixed = true;
                    self.vs_[j as usize].val = val;
                    self.solution_[j as usize] = val;
                    changed = true;
                }

                if changed {
                    self.update_activities(j, &pre_change_vs);
                    self.pq_notify(j, &pre_change_vs);
                    self.seed_worklist(j);
                }
            }
        }
        self.prop_work_ += prop_work;
        PropResult::kFixpoint
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::propagate::truncate (lines 459-466)
    fn truncate_fixpoint(&mut self, prop_work: usize, why: PropResult) -> PropResult {
        self.prop_work_ += prop_work;
        for wi in 0..self.prop_worklist_.len() {
            let w = self.prop_worklist_[wi];
            self.prop_in_wl_[w as usize] = 0;
        }
        self.prop_worklist_.clear();
        why
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::backtrack_to (lines 639-686)
    pub(crate) fn backtrack_to(
        &mut self,
        vs_mark_target: HighsInt,
        sol_mark_target: HighsInt,
        act_mark_target: HighsInt,
        pq_mark_target: HighsInt,
    ) {
        let mut u = self.vs_undo_.len() as HighsInt - 1;
        while u >= vs_mark_target {
            let (var, state) = self.vs_undo_[u as usize];
            self.vs_[var as usize] = state;
            u -= 1;
        }
        self.vs_undo_.truncate(vs_mark_target as usize);
        let mut u = self.sol_undo_.len() as HighsInt - 1;
        while u >= sol_mark_target {
            let (var, val) = self.sol_undo_[u as usize];
            self.solution_[var as usize] = val;
            u -= 1;
        }
        self.sol_undo_.truncate(sol_mark_target as usize);
        // Restore row activities
        if act_mark_target >= 0 {
            let mut u = self.act_undo_.len() as HighsInt - 1;
            while u >= act_mark_target {
                let undo = self.act_undo_[u as usize];
                self.min_finite_[undo.row as usize] = undo.old_min_finite;
                self.max_finite_[undo.row as usize] = undo.old_max_finite;
                self.min_ninf_[undo.row as usize] = undo.old_min_ninf;
                self.max_ninf_[undo.row as usize] = undo.old_max_ninf;
                self.publish_activity(undo.row);
                u -= 1;
            }
            self.act_undo_.truncate(act_mark_target as usize);
        }
        // Restore domain priority queue by replaying undo entries in reverse.
        // Each entry stores the exact PQ key before and after the change, so
        // we don't depend on vs_ state during replay.  Collapse the
        // erase-then-insert pair into a single update() when both old and new
        // states had the var in the PQ (the common bound-tightening case).
        if pq_mark_target >= 0 {
            let mut u = self.pq_undo_.len() as HighsInt - 1;
            while u >= pq_mark_target {
                let undo = self.pq_undo_[u as usize];
                if undo.is_present && undo.was_present {
                    if undo.old_dom != undo.new_dom {
                        self.domain_pq_.update(undo.var, undo.old_dom);
                    }
                } else if undo.is_present {
                    self.domain_pq_.erase(undo.var);
                } else if undo.was_present {
                    self.domain_pq_.insert(undo.old_dom, undo.var);
                }
                u -= 1;
            }
            self.pq_undo_.truncate(pq_mark_target as usize);
        }
        // Clear stale worklist entries that may reference now-reverted state.
        for wi in 0..self.prop_worklist_.len() {
            let w = self.prop_worklist_[wi];
            self.prop_in_wl_[w as usize] = 0;
        }
        self.prop_worklist_.clear();
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::vs_mark (lines 688-690)
    pub(crate) fn vs_mark(&self) -> HighsInt {
        self.vs_undo_.len() as HighsInt
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::sol_mark (lines 692-694)
    pub(crate) fn sol_mark(&self) -> HighsInt {
        self.sol_undo_.len() as HighsInt
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::act_mark (lines 696-698)
    pub(crate) fn act_mark(&self) -> HighsInt {
        self.act_undo_.len() as HighsInt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::fpr_local_mip::rng::build_csc;

    /// Helper: build a small test model for PropEngine tests.
    /// 3 variables: x0 (binary), x1 (integer [0,5]), x2 (continuous [0,10])
    /// 2 constraints:
    ///   row 0: x0 + x1 >= 2   (row_lo=2, row_hi=inf)
    ///   row 1: x1 + x2 <= 8   (row_lo=-inf, row_hi=8)
    struct SmallModel {
        ar_start: Vec<HighsInt>,
        ar_index: Vec<HighsInt>,
        ar_value: Vec<f64>,
        col_lb: Vec<f64>,
        col_ub: Vec<f64>,
        row_lo: Vec<f64>,
        row_hi: Vec<f64>,
        integrality: Vec<HighsVarType>,
        csc_start: Vec<HighsInt>,
        csc_row: Vec<HighsInt>,
        csc_val: Vec<f64>,
    }

    const K_NCOL: HighsInt = 3;
    const K_NROW: HighsInt = 2;

    impl SmallModel {
        fn new() -> Self {
            let ar_start = vec![0, 2, 4];
            let ar_index = vec![0, 1, 1, 2];
            let ar_value = vec![1.0, 1.0, 1.0, 1.0];
            let col_lb = vec![0.0, 0.0, 0.0];
            let col_ub = vec![1.0, 5.0, 10.0];
            let row_lo = vec![2.0, -kHighsInf];
            let row_hi = vec![kHighsInf, 8.0];
            let integrality = vec![
                HighsVarType::kInteger,
                HighsVarType::kInteger,
                HighsVarType::kContinuous,
            ];
            let csc = build_csc(K_NCOL, K_NROW, &ar_start, &ar_index, &ar_value);
            SmallModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb,
                col_ub,
                row_lo,
                row_hi,
                integrality,
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self) -> PropEngine<'_> {
            PropEngine::new(
                K_NCOL,
                K_NROW,
                &self.ar_start,
                &self.ar_index,
                &self.ar_value,
                &self.csc_start,
                &self.csc_row,
                &self.csc_val,
                &self.col_lb,
                &self.col_ub,
                &self.row_lo,
                &self.row_hi,
                &self.integrality,
                1e-6,
            )
        }
    }

    #[test]
    fn s5_propb_fix_and_propagate() {
        let m = SmallModel::new();
        let mut eng = m.make_engine();

        // Initial state: all variables unfixed with global bounds
        assert!(!eng.var(0).fixed);
        assert!((eng.var(1).lb - 0.0).abs() < 1e-6);
        assert!((eng.var(1).ub - 5.0).abs() < 1e-6);

        // Fix x0 = 0, propagate: row 0 (x0+x1 >= 2) forces x1 >= 2
        assert!(eng.fix(0, 0.0));
        assert_eq!(eng.propagate(0), PropResult::kFixpoint);
        assert!(eng.var(0).fixed);
        assert!((eng.var(0).val - 0.0).abs() < 1e-6);
        assert!(eng.var(1).lb >= 2.0 - 1e-6);

        // Fix x1 = 5, propagate: row 1 (x1+x2 <= 8) forces x2 <= 3
        assert!(eng.fix(1, 5.0));
        assert_eq!(eng.propagate(1), PropResult::kFixpoint);
        assert!(eng.var(2).ub <= 3.0 + 1e-6);
    }

    #[test]
    fn s5_propb_backtrack_restores_state() {
        let m = SmallModel::new();
        let mut eng = m.make_engine();

        let vs_m = eng.vs_mark();
        let sol_m = eng.sol_mark();

        assert!(eng.fix(0, 1.0));
        assert_eq!(eng.propagate(0), PropResult::kFixpoint);
        assert!(eng.var(0).fixed);

        eng.backtrack_to(vs_m, sol_m, -1, -1);
        assert!(!eng.var(0).fixed);
        assert!((eng.var(0).lb - 0.0).abs() < 1e-6);
        assert!((eng.var(0).ub - 1.0).abs() < 1e-6);
        assert!((eng.var(1).lb - 0.0).abs() < 1e-6);
    }

    #[test]
    fn s5_propb_tighten_and_auto_fix() {
        let m = SmallModel::new();
        let mut eng = m.make_engine();

        assert!(eng.tighten_lb(1, 3.0));
        assert!(eng.var(1).lb >= 3.0 - 1e-6);
        assert!((eng.var(1).ub - 5.0).abs() < 1e-6);

        // Tighten ub to match lb -- should auto-fix
        assert!(eng.tighten_ub(1, 3.0));
        assert!(eng.var(1).fixed);
        assert!((eng.var(1).val - 3.0).abs() < 1e-6);
    }

    #[test]
    fn s5_propb_infeasible_propagation() {
        let m = SmallModel::new();
        let mut eng = m.make_engine();

        // Fix x0=0, tighten x1 ub to 1. Propagation from row 0 (x0+x1 >= 2)
        // tries to tighten x1 lb to 2, but ub is 1 -> lb > ub -> infeasible.
        assert!(eng.fix(0, 0.0));
        assert!(eng.tighten_ub(1, 1.0));
        assert_eq!(eng.propagate(0), PropResult::kInfeasible);
    }

    #[test]
    fn s5_propb_reset_clears_state() {
        let m = SmallModel::new();
        let mut eng = m.make_engine();

        eng.fix(0, 1.0);
        let _ = eng.propagate(0);
        assert!(eng.var(0).fixed);

        eng.reset();
        assert!(!eng.var(0).fixed);
        assert!((eng.var(0).lb - 0.0).abs() < 1e-6);
        assert!((eng.var(0).ub - 1.0).abs() < 1e-6);
        assert!((eng.var(1).lb - 0.0).abs() < 1e-6);
        assert!((eng.var(1).ub - 5.0).abs() < 1e-6);
    }

    #[test]
    fn s5_propb_effort_tracking() {
        let m = SmallModel::new();
        let mut eng = m.make_engine();

        let before = eng.effort();
        eng.fix(0, 1.0);
        let _ = eng.propagate(0);
        assert!(eng.effort() > before);

        eng.add_effort(100);
        assert!(eng.effort() >= before + 100);
    }

    /// HalvingModel (test_prop_engine.cpp:393-413): 2 variables, 2 rows.
    /// row0: x0 - 2*x1 >= 0, row1: x1 - 2*x0 >= 0, both continuous.
    /// col bounds [0, 1e300], row_lo = 0, row_hi = +infinity.
    /// Processing row0 tightens x1's ub to roughly half of x0's current ub,
    /// then seed_worklist re-adds row1; processing row1 symmetrically halves
    /// x0's ub and re-seeds row0. Convergence takes ~log2(col_ub/feastol)
    /// row visits (~1000) against nnz=4, so a single propagate() call hits
    /// the 100*4=400 budget long before reaching a fixpoint.
    struct HalvingModel {
        ar_start: Vec<HighsInt>,
        ar_index: Vec<HighsInt>,
        ar_value: Vec<f64>,
        col_lb: Vec<f64>,
        col_ub: Vec<f64>,
        row_lo: Vec<f64>,
        row_hi: Vec<f64>,
        integrality: Vec<HighsVarType>,
        csc_start: Vec<HighsInt>,
        csc_row: Vec<HighsInt>,
        csc_val: Vec<f64>,
    }

    const K_HALVING_NCOL: HighsInt = 2;
    const K_HALVING_NROW: HighsInt = 2;

    impl HalvingModel {
        fn new() -> Self {
            let ar_start = vec![0, 2, 4];
            let ar_index = vec![0, 1, 1, 0];
            let ar_value = vec![1.0, -2.0, 1.0, -2.0];
            let col_lb = vec![0.0, 0.0];
            let col_ub = vec![1e300, 1e300];
            let row_lo = vec![0.0, 0.0];
            let row_hi = vec![kHighsInf, kHighsInf];
            let integrality = vec![HighsVarType::kContinuous, HighsVarType::kContinuous];
            let csc = build_csc(
                K_HALVING_NCOL,
                K_HALVING_NROW,
                &ar_start,
                &ar_index,
                &ar_value,
            );
            HalvingModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb,
                col_ub,
                row_lo,
                row_hi,
                integrality,
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self) -> PropEngine<'_> {
            PropEngine::new(
                K_HALVING_NCOL,
                K_HALVING_NROW,
                &self.ar_start,
                &self.ar_index,
                &self.ar_value,
                &self.csc_start,
                &self.csc_row,
                &self.csc_val,
                &self.col_lb,
                &self.col_ub,
                &self.row_lo,
                &self.row_hi,
                &self.integrality,
                1e-6,
            )
        }
    }

    #[test]
    fn s5_propb_budget_exhaustion() {
        let m = HalvingModel::new();
        let mut eng = m.make_engine();

        let vs_m = eng.vs_mark();
        let sol_m = eng.sol_mark();

        // Column 0 (x) is in both rows, so seeding it starts the cascade.
        eng.seed_worklist(0);
        let result = eng.propagate(-1);

        assert_eq!(result, PropResult::kBudgetExhausted);
        // Genuinely exhausted (100 * nnz = 100 * 4 = 400), not a fluke early
        // return on the first row visited.
        assert!(eng.effort() > 400);

        // Sound: both domains are still non-empty (lb <= ub). An infeasible
        // exit -- the other early-return path in propagate() -- is exactly
        // the case lb > ub on some row; this must not be that.
        assert!(eng.var(0).lb <= eng.var(0).ub + eng.feastol());
        assert!(eng.var(1).lb <= eng.var(1).ub + eng.feastol());
        // Neither variable was spuriously auto-fixed by the partial run.
        assert!(!eng.var(0).fixed);
        assert!(!eng.var(1).fixed);
        // Genuine partial progress happened (budget was actually spent doing
        // real deductions, not just counted and discarded).
        assert!(eng.var(0).ub < 1e300);
        assert!(eng.var(1).ub < 1e300);

        // Undo-stack consistency: a caller that continues past a
        // kBudgetExhausted node (rather than pruning it) must still be able
        // to backtrack cleanly through everything this call applied.
        eng.backtrack_to(vs_m, sol_m, -1, -1);
        assert!((eng.var(0).lb - 0.0).abs() < 1e-6);
        assert!((eng.var(0).ub - 1e300).abs() < 1e294);
        assert!((eng.var(1).lb - 0.0).abs() < 1e-6);
        assert!((eng.var(1).ub - 1e300).abs() < 1e294);
    }
}
