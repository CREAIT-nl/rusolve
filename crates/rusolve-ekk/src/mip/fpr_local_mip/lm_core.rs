//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! WorkerCtx (local_mip_core.h/cpp), LiftCache (local_mip_caches.h),
//! IndexedSet, and the free helper round_tight_delta (ViolCache is not
//! ported, see lm_search.rs::compute_candidate_scores).

use super::rng::{S5Rng, clamp_round, is_integer, row_violation, s5_uniform_real};
use super::*;
use crate::prelude::*;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// mip-heuristics (36f7d9d): local_mip_core.h::round_tight_delta (lines 52-57)
pub(crate) fn round_tight_delta(delta: f64, row_violated: bool) -> f64 {
    if row_violated {
        if delta > 0.0 {
            delta.ceil()
        } else {
            delta.floor()
        }
    } else if delta > 0.0 {
        delta.floor()
    } else {
        delta.ceil()
    }
}

impl IndexedSet {
    /// mip-heuristics (36f7d9d): local_mip_caches.h::IndexedSet::IndexedSet (lines 119-119)
    pub(crate) fn new(n: HighsInt) -> Self {
        IndexedSet {
            elements: Vec::with_capacity(n as usize),
            pos: vec![-1; n as usize],
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_caches.h::IndexedSet::add (lines 121-127)
    pub(crate) fn add(&mut self, i: HighsInt) {
        if self.pos[i as usize] != -1 {
            return;
        }
        self.pos[i as usize] = self.elements.len() as HighsInt;
        self.elements.push(i);
    }

    /// mip-heuristics (36f7d9d): local_mip_caches.h::IndexedSet::remove (lines 129-139)
    pub(crate) fn remove(&mut self, i: HighsInt) {
        let p = self.pos[i as usize];
        if p == -1 {
            return;
        }
        let last = *self.elements.last().unwrap();
        self.elements[p as usize] = last;
        self.pos[last as usize] = p;
        self.elements.pop();
        self.pos[i as usize] = -1;
    }

    /// mip-heuristics (36f7d9d): local_mip_caches.h::IndexedSet::contains (lines 141-141)
    pub(crate) fn contains(&self, i: HighsInt) -> bool {
        self.pos[i as usize] != -1
    }

    /// mip-heuristics (36f7d9d): local_mip_caches.h::IndexedSet::empty (lines 142-142)
    pub(crate) fn empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// mip-heuristics (36f7d9d): local_mip_caches.h::IndexedSet::size (lines 143-143)
    pub(crate) fn size(&self) -> HighsInt {
        self.elements.len() as HighsInt
    }

    /// mip-heuristics (36f7d9d): local_mip_caches.h::IndexedSet::clear (lines 146-151)
    pub(crate) fn clear(&mut self) {
        for &e in &self.elements {
            self.pos[e as usize] = -1;
        }
        self.elements.clear();
    }
}

impl LiftCache {
    /// mip-heuristics (36f7d9d): local_mip_caches.h::LiftCache::LiftCache (lines 214-218)
    pub(crate) fn new(ncol: HighsInt) -> Self {
        LiftCache {
            lo: vec![0.0; ncol as usize],
            hi: vec![0.0; ncol as usize],
            score: vec![0.0; ncol as usize],
            dirty: vec![true; ncol as usize],
            dirty_list: Vec::with_capacity(ncol as usize),
            all_dirty: true,
            positive_list: Vec::with_capacity(ncol as usize),
            in_positive: vec![false; ncol as usize],
            costed_vars: None,
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_caches.h::LiftCache::mark_dirty (lines 220-225)
    pub(crate) fn mark_dirty(&mut self, j: HighsInt) {
        if !self.dirty[j as usize] {
            self.dirty[j as usize] = true;
            self.dirty_list.push(j);
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_caches.h::LiftCache::mark_all_dirty (lines 227-233)
    pub(crate) fn mark_all_dirty(&mut self) {
        self.all_dirty = true;
        self.dirty_list.clear();
        self.dirty.fill(true);
        self.positive_list.clear();
        self.in_positive.fill(false);
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::LiftCache::recompute_one (lines 278-349)
    ///
    /// Callers take `ctx.lift` out with `std::mem::take` for the call and put it
    /// back (lib.rs rule f); LiftCache implements Default in mod.rs.
    pub(crate) fn recompute_one(&mut self, j: HighsInt, ctx: &mut WorkerCtx) {
        let old_score = self.score[j as usize];
        if ctx.col_cost[j as usize].abs() < kEpsZero {
            self.score[j as usize] = 0.0;
            self.dirty[j as usize] = false;
            if old_score > 0.0 && self.in_positive[j as usize] {
                self.in_positive[j as usize] = false;
                // lazy removal: stale entries filtered during scan
            }
            return;
        }
        // Compute lift bounds
        ctx.ls_evals += 1;
        ctx.effort +=
            (ctx.csc.col_start[(j + 1) as usize] - ctx.csc.col_start[j as usize]) as usize;
        let mut lo_j = ctx.col_lb[j as usize];
        let mut hi_j = ctx.col_ub[j as usize];
        let mut p = ctx.csc.col_start[j as usize];
        while p < ctx.csc.col_start[(j + 1) as usize] {
            let i = ctx.csc.col_row[p as usize];
            let coeff = ctx.csc.col_val[p as usize];
            if coeff.abs() < kEpsZero {
                p += 1;
                continue;
            }
            let residual = ctx.lhs[i as usize] - (coeff * ctx.solution[j as usize]);
            if ctx.row_hi[i as usize] < kHighsInf {
                let bound = (ctx.row_hi[i as usize] - residual) / coeff;
                if coeff > 0.0 {
                    hi_j = cpp_min(hi_j, bound);
                } else {
                    lo_j = cpp_max(lo_j, bound);
                }
            }
            if ctx.row_lo[i as usize] > -kHighsInf {
                let bound = (ctx.row_lo[i as usize] - residual) / coeff;
                if coeff > 0.0 {
                    lo_j = cpp_max(lo_j, bound);
                } else {
                    hi_j = cpp_min(hi_j, bound);
                }
            }
            p += 1;
        }
        if ctx.is_int(j) {
            lo_j = (lo_j - ctx.feastol).ceil();
            hi_j = (hi_j + ctx.feastol).floor();
        }
        self.lo[j as usize] = lo_j;
        self.hi[j as usize] = hi_j;

        if lo_j > hi_j {
            self.score[j as usize] = 0.0;
        } else {
            let mut target = if ctx.col_cost[j as usize] > 0.0 {
                lo_j
            } else {
                hi_j
            };
            target = ctx.clamp_and_round(j, target);
            if (target - ctx.solution[j as usize]).abs() < kEpsZero {
                self.score[j as usize] = 0.0;
            } else {
                let obj_delta = ctx.col_cost[j as usize] * (target - ctx.solution[j as usize]);
                self.score[j as usize] = -obj_delta; // positive = improving
            }
        }
        // Maintain positive-lift list
        if self.score[j as usize] > 0.0 {
            if !self.in_positive[j as usize] {
                self.in_positive[j as usize] = true;
                self.positive_list.push(j);
            }
        } else {
            if self.in_positive[j as usize] {
                self.in_positive[j as usize] = false;
                // lazy removal: stale entries filtered during scan
            }
        }
        self.dirty[j as usize] = false;
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::LiftCache::recompute_all (lines 351-374)
    ///
    /// Callers take `ctx.lift` out with `std::mem::take` for the call and put it
    /// back (lib.rs rule f); LiftCache implements Default in mod.rs.
    pub(crate) fn recompute_all(&mut self, ctx: &mut WorkerCtx) {
        if self.all_dirty {
            // Only recompute columns with nonzero cost; zero-cost columns
            // always have score=0 and never need lift recomputation.
            let costed = self.costed_vars.clone();
            if let Some(ref cv) = costed {
                for &j in cv.as_ref() {
                    self.recompute_one(j, ctx);
                }
            } else {
                let mut j = 0 as HighsInt;
                while j < ctx.ncol {
                    self.recompute_one(j, ctx);
                    j += 1;
                }
            }
            self.dirty_list.clear();
        } else {
            let dirty_snapshot: Vec<HighsInt> = self.dirty_list.clone();
            for &j in &dirty_snapshot {
                if self.dirty[j as usize] {
                    self.recompute_one(j, ctx);
                }
            }
            self.dirty_list.clear();
        }
        self.all_dirty = false;
    }
}

impl<'a> WorkerCtx<'a> {
    /// mip-heuristics (36f7d9d): local_mip_core.cpp::WorkerCtx::WorkerCtx (lines 17-45)
    pub(crate) fn new(problem: &'a S5Problem) -> Self {
        let ncol = problem.ncol;
        let nrow = problem.nrow;
        WorkerCtx {
            ar_start: &problem.ar_start,
            ar_index: &problem.ar_index,
            ar_value: &problem.ar_value,
            col_lb: &problem.col_lower,
            col_ub: &problem.col_upper,
            col_cost: &problem.col_cost,
            offset: problem.offset,
            row_lo: &problem.row_lower,
            row_hi: &problem.row_upper,
            integrality: &problem.integrality,
            csc: &problem.csc,
            feastol: problem.feastol,
            epsilon: problem.epsilon,
            ncol,
            nrow,
            binary: &problem.binary,
            solution: vec![0.0; ncol as usize],
            lhs: vec![0.0; nrow as usize],
            weight: vec![1u64; nrow as usize],
            obj_weight: 1,
            current_obj: 0.0,
            violated: IndexedSet::new(nrow),
            satisfied: IndexedSet::new(nrow),
            lift: LiftCache::new(ncol),
            tabu_inc_until: vec![0; ncol as usize],
            tabu_dec_until: vec![0; ncol as usize],
            batch: Vec::with_capacity(kBmsBudget as usize),
            sampled: Vec::with_capacity((kBmsConstraints * 3) as usize),
            was_infeasible: true,
            feasible_recheck_counter: 0,
            effort: 0,
            track_violation: false,
            total_violation: 0.0,
            bt_live: vec![0; (ncol as usize).div_ceil(64)],
            bt_live_valid: false,
            bt_costed: 0,
            ls_evals: 0,
            ls_commits: 0,
            gub: None,
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx::is_int (lines 128-128)
    pub(crate) fn is_int(&self, j: HighsInt) -> bool {
        is_integer(self.integrality, j)
    }

    /// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx::is_binary (lines 130-130)
    pub(crate) fn is_binary(&self, j: HighsInt) -> bool {
        self.binary[j as usize] != 0
    }

    /// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx::clamp_and_round (lines 132-134)
    pub(crate) fn clamp_and_round(&self, j: HighsInt, val: f64) -> f64 {
        clamp_round(
            val,
            self.col_lb[j as usize],
            self.col_ub[j as usize],
            self.is_int(j),
        )
    }

    /// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx::compute_violation (lines 136-138)
    pub(crate) fn compute_violation(&self, i: HighsInt, l: f64) -> f64 {
        row_violation(l, self.row_lo[i as usize], self.row_hi[i as usize])
    }

    /// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx::is_violated (lines 140-142)
    pub(crate) fn is_violated(&self, i: HighsInt, l: f64) -> bool {
        l > self.row_hi[i as usize] + self.feastol || l < self.row_lo[i as usize] - self.feastol
    }

    /// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx::is_equality (lines 144-146)
    pub(crate) fn is_equality(&self, i: HighsInt) -> bool {
        self.row_lo[i as usize] == self.row_hi[i as usize] && self.row_lo[i as usize] > -kHighsInf
    }

    /// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx::is_tabu (lines 148-156)
    pub(crate) fn is_tabu(&self, j: HighsInt, delta: f64, step: HighsInt) -> bool {
        if delta > 0.0 && step < self.tabu_inc_until[j as usize] {
            return true;
        }
        if delta < 0.0 && step < self.tabu_dec_until[j as usize] {
            return true;
        }
        false
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::WorkerCtx::update_violated (lines 47-68)
    pub(crate) fn update_violated(&mut self, i: HighsInt) {
        let was = self.violated.contains(i);
        // Calls the shared predicate rather than re-deriving it (issue
        // #148 follow-up): `is_violated` uses `l > hi + feastol ||
        // l < lo - feastol`, which is not bitwise identical to
        // `row_violation(...) > feastol` at the ~1 ulp scale (`fl(hi+ft) !=
        // hi + fl(ft)` in general) — a residue of the same class of
        // mismatch #148 removed, just far smaller. Sharing the call site
        // makes set membership agree with `is_violated`,
        // `compute_tight_delta`'s branch and HiGHS's own `trySolution` up
        // to the bit, not just up to the value of `feastol`.
        let now = self.is_violated(i, self.lhs[i as usize]);
        if now && !was {
            self.violated.add(i);
            self.satisfied.remove(i);
        } else if !now && was {
            self.violated.remove(i);
            if !self.is_equality(i) {
                self.satisfied.add(i);
            }
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::WorkerCtx::apply_move (lines 70-101)
    pub(crate) fn apply_move(&mut self, j: HighsInt, new_val: f64) {
        let old_val = self.solution[j as usize];
        let delta = new_val - old_val;
        if delta.abs() < kEpsZero {
            return;
        }
        self.solution[j as usize] = new_val;
        self.ls_commits += 1;
        self.bt_update(j);
        self.current_obj += self.col_cost[j as usize] * delta;
        self.effort +=
            (self.csc.col_start[(j + 1) as usize] - self.csc.col_start[j as usize]) as usize;
        // Only maintain LiftCache during feasible mode; on the
        // infeasible→feasible transition, rebuild_state marks all dirty.
        let dirty_lift = !self.was_infeasible && !self.lift.all_dirty;
        if dirty_lift {
            self.lift.mark_dirty(j);
        }
        let mut p = self.csc.col_start[j as usize];
        while p < self.csc.col_start[(j + 1) as usize] {
            let i = self.csc.col_row[p as usize];
            if self.track_violation {
                self.total_violation -= self.compute_violation(i, self.lhs[i as usize]);
            }
            self.lhs[i as usize] += self.csc.col_val[p as usize] * delta;
            if self.track_violation {
                self.total_violation += self.compute_violation(i, self.lhs[i as usize]);
            }
            self.update_violated(i);
            if dirty_lift {
                let mut k = self.ar_start[i as usize];
                while k < self.ar_start[(i + 1) as usize] {
                    self.lift.mark_dirty(self.ar_index[k as usize]);
                    k += 1;
                }
            }
            p += 1;
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::WorkerCtx::apply_move_with_tabu (lines 103-112)
    pub(crate) fn apply_move_with_tabu(
        &mut self,
        j: HighsInt,
        new_val: f64,
        step: HighsInt,
        rng: &mut S5Rng,
    ) {
        let delta = new_val - self.solution[j as usize];
        self.apply_move(j, new_val);
        let tabu_len = kTabuBase + (rng.next_u64() % (kTabuVar as u64)) as HighsInt;
        if delta > 0.0 {
            self.tabu_dec_until[j as usize] = step + tabu_len;
        } else {
            self.tabu_inc_until[j as usize] = step + tabu_len;
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::WorkerCtx::full_recheck (lines 114-163)
    pub(crate) fn full_recheck(&mut self, update_sets: bool, early_exit: bool) -> bool {
        self.effort += self.ar_index.len();
        if update_sets {
            self.violated.clear();
            self.satisfied.clear();
        }
        let mut feasible = true;
        let mut i = 0 as HighsInt;
        while i < self.nrow {
            let mut l = 0.0;
            let mut k = self.ar_start[i as usize];
            while k < self.ar_start[(i + 1) as usize] {
                l += self.ar_value[k as usize] * self.solution[self.ar_index[k as usize] as usize];
                k += 1;
            }
            self.lhs[i as usize] = l;
            // `is_violated`, not a re-derived `compute_violation(...) >
            // feastol` (issue #148 follow-up — see `update_violated`'s
            // comment for why the two forms are not bitwise identical):
            // this is both the violated/satisfied set-membership classifier
            // and the worker's submission gate. Its *row* checks must
            // accept exactly what HiGHS's own `trySolution`
            // (mip/HighsMipSolverData.cpp) would — no looser, no stricter.
            // `trySolution` also checks column bounds and integrality
            // fractionality, which `full_recheck` does not duplicate here:
            // every candidate this worker generates is already clamped and
            // rounded to bounds/integrality by `clamp_and_round`
            // (`append_candidate` in `local_mip_search.cpp`; the starting
            // solution goes through an equivalent clamp/round in
            // `LocalMipWorker`'s constructor, `local_mip_worker.cpp`), so
            // `solution` never drifts outside them between recheck calls.
            if self.is_violated(i, l) {
                feasible = false;
                if early_exit {
                    break;
                }
                if update_sets {
                    self.violated.add(i);
                }
            } else if update_sets && !self.is_equality(i) {
                self.satisfied.add(i);
            }
            i += 1;
        }
        // Over every row, as `lhs` now stands: an early exit leaves the rows
        // past it at their previous values, which the running total also held.
        if self.track_violation {
            self.total_violation = 0.0;
            let mut i2 = 0 as HighsInt;
            while i2 < self.nrow {
                self.total_violation += self.compute_violation(i2, self.lhs[i2 as usize]);
                i2 += 1;
            }
        }
        feasible
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::WorkerCtx::rebuild_state (lines 165-174)
    pub(crate) fn rebuild_state(&mut self) {
        self.was_infeasible = true;
        self.feasible_recheck_counter = 0;
        self.full_recheck(true, false);
        self.lift.mark_all_dirty();
        self.current_obj = self.offset;
        let mut j = 0 as HighsInt;
        while j < self.ncol {
            self.current_obj += self.col_cost[j as usize] * self.solution[j as usize];
            j += 1;
        }
        // rusolve (lm-throughput): every bulk write of `solution` is
        // followed by this call, so the breakthrough live set is rebuilt here.
        let mut j = 0 as HighsInt;
        self.bt_costed = 0;
        while j < self.ncol {
            self.bt_update(j);
            self.bt_costed += usize::from(self.col_cost[j as usize].abs() >= kEpsZero);
            j += 1;
        }
        self.bt_live_valid = true;
    }

    /// rusolve (throughput): whether costed
    /// column `j` can make `compute_breakthrough_delta` + `append_candidate`
    /// push a candidate. A column whose value sits exactly on the finite
    /// bound its cost favours (`x == lb` for `c > 0`, `x == ub` for `c < 0`)
    /// and is already clamped and rounded cannot. With a non-negative
    /// `obj_gap + epsilon` the delta points out of the box (or is absorbed
    /// by `x + delta == x`), the bound clamp sends it back to
    /// `bound - x == 0`, and `append_candidate` rejects a zero move; with
    /// `cur_obj < best_obj - epsilon` the delta is 0 for every column.
    /// `infeasible_step` checks the gap precondition before it relies on
    /// this.
    fn bt_can_yield(&self, j: HighsInt) -> bool {
        let ju = j as usize;
        let cost = self.col_cost[ju];
        // The worker's `costed_vars_` predicate, NaN included.
        if !(cost.abs() >= kEpsZero) {
            return false; // not a costed column: never iterated
        }
        let x = self.solution[ju];
        let parked = if cost > 0.0 {
            x == self.col_lb[ju] && self.col_lb[ju].is_finite()
        } else if cost < 0.0 {
            x == self.col_ub[ju] && self.col_ub[ju].is_finite()
        } else {
            false // NaN cost
        };
        !(parked && (self.clamp_and_round(j, x) - x).abs() < kEpsZero)
    }

    /// Refresh bit `j` of `bt_live` after `solution[j]` changed.
    pub(crate) fn bt_update(&mut self, j: HighsInt) {
        let (word, bit) = ((j as usize) / 64, (j as usize) % 64);
        if self.bt_can_yield(j) {
            self.bt_live[word] |= 1u64 << bit;
        } else {
            self.bt_live[word] &= !(1u64 << bit);
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::WorkerCtx::compute_tight_delta (lines 181-244)
    pub(crate) fn compute_tight_delta(&self, i: HighsInt, j: HighsInt, coeff: f64) -> f64 {
        if coeff.abs() < kEpsZero {
            return 0.0;
        }
        let row = self.tight_row(i);
        self.tight_delta_in_row(&row, j, coeff)
    }

    /// The row half of `compute_tight_delta`: everything that depends on
    /// row `i` alone. rusolve (throughput):
    /// split out so the Phase 1/2/4 loops of `infeasible_step` evaluate it
    /// once per sampled row instead of once per nonzero; `lhs` is fixed
    /// while a batch is built, so the deltas are unchanged.
    pub(crate) fn tight_row(&self, i: HighsInt) -> TightRow {
        let l = self.lhs[i as usize];
        let gap: f64;
        let row_violated: bool;
        // The first two branches share an expression form but differ in which
        // bound they measure against.  (This used to carry a NOLINTBEGIN for
        // bugprone-branch-clone; verified against clang-tidy 22 that the check
        // does not fire here, so the suppression was dead and is gone.)
        if l > self.row_hi[i as usize] + self.feastol {
            gap = l - self.row_hi[i as usize]; // upper violated
            row_violated = true;
        } else if l < self.row_lo[i as usize] - self.feastol {
            gap = l - self.row_lo[i as usize]; // lower violated
            row_violated = true;
        } else {
            // Satisfied: push toward the nearest bound
            let gap_hi = if self.row_hi[i as usize] < kHighsInf {
                l - self.row_hi[i as usize]
            } else {
                kHighsInf
            };
            let gap_lo = if self.row_lo[i as usize] > -kHighsInf {
                l - self.row_lo[i as usize]
            } else {
                kHighsInf
            };
            gap = if gap_hi.abs() <= gap_lo.abs() {
                gap_hi
            } else {
                gap_lo
            };
            row_violated = false;
        }
        TightRow {
            gap,
            row_violated,
            is_equality: self.is_equality(i),
        }
    }

    /// The column half of `compute_tight_delta`, for a row from `tight_row`.
    pub(crate) fn tight_delta_in_row(&self, row: &TightRow, j: HighsInt, coeff: f64) -> f64 {
        if coeff.abs() < kEpsZero {
            return 0.0;
        }
        let gap = row.gap;
        let row_violated = row.row_violated;
        let mut delta = -gap / coeff;

        if row.is_equality {
            // Deliberately keeps the coefficient-sign rule rather than
            // `round_tight_delta` (issue #123): a violated equality row cannot be
            // made satisfied by any rounding direction — both sides of the bound
            // are infeasible — so "round away from zero" has nothing to reach.
            // The satisfied half of this branch is unreachable in practice: both
            // `satisfied.add` sites guard on `!is_equality`, so no equality row
            // enters the set the satisfied-row callers iterate.
            if self.is_int(j) {
                delta = if coeff > 0.0 {
                    delta.floor()
                } else {
                    delta.ceil()
                };
            }
            let new_val = self.solution[j as usize] + delta;
            if new_val < self.col_lb[j as usize] || new_val > self.col_ub[j as usize] {
                if (gap > 0.0 && coeff > 0.0) || (gap < 0.0 && coeff < 0.0) {
                    delta = self.col_lb[j as usize] - self.solution[j as usize];
                } else {
                    delta = self.col_ub[j as usize] - self.solution[j as usize];
                }
            }
        } else {
            // Integer rounding follows the row's state, not the coefficient's
            // sign: away from zero when the row is violated, toward zero when it
            // is satisfied.  See `round_tight_delta` for the rule and the paper
            // citation (issue #123).
            if self.is_int(j) {
                delta = round_tight_delta(delta, row_violated);
            }
            // Clamp to variable bounds (Paper Eq 5 min/max with l_j, u_j).
            let new_val = self.solution[j as usize] + delta;
            if new_val < self.col_lb[j as usize] {
                delta = self.col_lb[j as usize] - self.solution[j as usize];
            } else if new_val > self.col_ub[j as usize] {
                delta = self.col_ub[j as usize] - self.solution[j as usize];
            }
        }
        delta
    }

    /// mip-heuristics (36f7d9d): local_mip_core.cpp::WorkerCtx::update_weights (lines 246-269)
    pub(crate) fn update_weights(
        &mut self,
        rng: &mut S5Rng,
        is_feasible: bool,
        best_feasible: bool,
        best_obj: f64,
    ) {
        let coin = s5_uniform_real(rng, 0.0, 1.0);
        if coin >= kSmoothProb {
            // With probability 1 - sp: strengthen
            if is_feasible {
                self.obj_weight += 1;
            } else {
                for &ci in &self.violated.elements {
                    self.weight[ci as usize] += 1;
                }
            }
        } else {
            // With probability sp: smooth (weaken)
            let obj_better = best_feasible && self.current_obj < best_obj;
            if obj_better && self.obj_weight > 1 {
                self.obj_weight -= 1;
            }
            for &ci in &self.satisfied.elements {
                if self.weight[ci as usize] > 1 {
                    self.weight[ci as usize] -= 1;
                }
            }
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx::reset_weights (lines 190-193)
    pub(crate) fn reset_weights(&mut self) {
        self.weight.fill(1u64);
        self.obj_weight = 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;

    /// Helper: build a small S5Problem with 2 rows, 2 cols.
    /// Row 0: x0 + x1 <= 3  (row_lo=-inf, row_hi=3)
    /// Row 1: x0 + x1 >= 1  (row_lo=1, row_hi=inf)
    /// Both integer.
    fn make_problem_2x2() -> S5Problem {
        // CSR: row 0 has entries at col 0 (coeff 1) and col 1 (coeff 1)
        //       row 1 has entries at col 0 (coeff 1) and col 1 (coeff 1)
        // CSC:  col 0 has entries at row 0 (coeff 1) and row 1 (coeff 1)
        //       col 1 has entries at row 0 (coeff 1) and row 1 (coeff 1)
        S5Problem {
            col_lower: vec![0.0, 0.0],
            col_upper: vec![10.0, 10.0],
            row_lower: vec![-kHighsInf, 1.0],
            row_upper: vec![3.0, kHighsInf],
            col_cost: vec![1.0, 1.0],
            offset: 0.0,
            integrality: vec![HighsVarType::kInteger, HighsVarType::kInteger],
            ar_start: vec![0, 2, 4],
            ar_index: vec![0, 1, 0, 1],
            ar_value: vec![1.0, 1.0, 1.0, 1.0],
            csc: CscMatrix {
                col_start: vec![0, 2, 4],
                col_row: vec![0, 1, 0, 1],
                col_val: vec![1.0, 1.0, 1.0, 1.0],
            },
            uplocks: vec![0, 0],
            downlocks: vec![0, 0],
            feastol: 1e-6,
            epsilon: 1e-6,
            ncol: 2,
            nrow: 2,
            nnz: 4,
            incumbent: vec![0.0, 0.0],
            binary: vec![0, 0],
        }
    }

    #[test]
    fn s5_lmcore_apply_move_updates_lhs_and_violated() {
        let problem = make_problem_2x2();
        let mut ctx = WorkerCtx::new(&problem);
        // Start at x0=0, x1=0: lhs = [0, 0].
        // Row 0: 0 <= 3 satisfied. Row 1: 0 >= 1 violated.
        ctx.solution[0] = 0.0;
        ctx.solution[1] = 0.0;
        let feas = ctx.full_recheck(true, false);
        assert!(!feas, "should be infeasible at origin");
        assert!(ctx.violated.contains(1), "row 1 should be violated");
        assert!(!ctx.violated.contains(0), "row 0 should not be violated");
        assert!(ctx.satisfied.contains(0), "row 0 should be satisfied");
        assert!(!ctx.satisfied.contains(1), "row 1 should not be satisfied");

        // Move x0 to 2: lhs becomes [2, 2].
        // Row 0: 2 <= 3 satisfied. Row 1: 2 >= 1 satisfied.
        ctx.apply_move(0, 2.0);
        assert!((ctx.lhs[0] - 2.0).abs() < 1e-12);
        assert!((ctx.lhs[1] - 2.0).abs() < 1e-12);
        assert!(!ctx.violated.contains(0), "row 0 still satisfied");
        assert!(!ctx.violated.contains(1), "row 1 now satisfied");
        assert!(ctx.satisfied.contains(0), "row 0 still satisfied");
        assert!(ctx.satisfied.contains(1), "row 1 now satisfied");
    }

    #[test]
    fn s5_lmcore_full_recheck_counts_violations() {
        let problem = make_problem_2x2();
        let mut ctx = WorkerCtx::new(&problem);
        // x0=5, x1=5: lhs = [10, 10].
        // Row 0: 10 <= 3 violated. Row 1: 10 >= 1 satisfied.
        ctx.solution[0] = 5.0;
        ctx.solution[1] = 5.0;
        let feas = ctx.full_recheck(true, false);
        assert!(!feas, "should be infeasible");
        assert!(ctx.violated.contains(0), "row 0 violated");
        assert!(!ctx.violated.contains(1), "row 1 not violated");
        assert!(!ctx.satisfied.contains(0), "row 0 not satisfied");
        assert!(ctx.satisfied.contains(1), "row 1 satisfied");
        assert_eq!(ctx.violated.size(), 1);
    }

    #[test]
    fn s5_lmcore_compute_tight_delta_le_ge_eq() {
        // Test compute_tight_delta on <=, >=, and = rows.
        // Build a 1-row, 1-col problem for each case.

        // Case 1: <= row, violated. x0 <= 3, at x0=5 => lhs=10, violated.
        // gap = 10 - 3 = 7, coeff=1, delta = -7/1 = -7.
        // row_violated=true, is_int=true => round_tight_delta(-7, true) = floor(-7) = -7.
        // new_val = 5 + (-7) = -2, col_lb=0 => new_val < col_lb => delta = 0 - 5 = -5.
        {
            let problem = S5Problem {
                col_lower: vec![0.0],
                col_upper: vec![10.0],
                row_lower: vec![-kHighsInf],
                row_upper: vec![3.0],
                col_cost: vec![1.0],
                offset: 0.0,
                integrality: vec![HighsVarType::kInteger],
                ar_start: vec![0, 1],
                ar_index: vec![0],
                ar_value: vec![1.0],
                csc: CscMatrix {
                    col_start: vec![0, 1],
                    col_row: vec![0],
                    col_val: vec![1.0],
                },
                uplocks: vec![0],
                downlocks: vec![0],
                feastol: 1e-6,
                epsilon: 1e-6,
                ncol: 1,
                nrow: 1,
                nnz: 1,
                incumbent: vec![0.0],
                binary: vec![0],
            };
            let mut ctx = WorkerCtx::new(&problem);
            ctx.solution[0] = 5.0;
            ctx.lhs[0] = 10.0;
            let delta = ctx.compute_tight_delta(0, 0, 1.0);
            // gap = 10 - 3 = 7, delta_raw = -7, row_violated=true, int => ceil(-7)=-7
            // new_val = 5 + (-7) = -2 < col_lb=0 => delta = 0 - 5 = -5
            assert!((delta - (-5.0)).abs() < 1e-12, "delta={}", delta);
        }

        // Case 2: >= row, violated. x0 >= 5, at x0=2 => lhs=2, violated.
        // gap = 2 - 5 = -3, coeff=1, delta = -(-3)/1 = 3.
        // row_violated=true, is_int=true => round_tight_delta(3, true) = ceil(3) = 3.
        // new_val = 2 + 3 = 5, within [0,10] => delta = 3.
        {
            let problem = S5Problem {
                col_lower: vec![0.0],
                col_upper: vec![10.0],
                row_lower: vec![5.0],
                row_upper: vec![kHighsInf],
                col_cost: vec![1.0],
                offset: 0.0,
                integrality: vec![HighsVarType::kInteger],
                ar_start: vec![0, 1],
                ar_index: vec![0],
                ar_value: vec![1.0],
                csc: CscMatrix {
                    col_start: vec![0, 1],
                    col_row: vec![0],
                    col_val: vec![1.0],
                },
                uplocks: vec![0],
                downlocks: vec![0],
                feastol: 1e-6,
                epsilon: 1e-6,
                ncol: 1,
                nrow: 1,
                nnz: 1,
                incumbent: vec![0.0],
                binary: vec![0],
            };
            let mut ctx = WorkerCtx::new(&problem);
            ctx.solution[0] = 2.0;
            ctx.lhs[0] = 2.0;
            let delta = ctx.compute_tight_delta(0, 0, 1.0);
            assert!((delta - 3.0).abs() < 1e-12, "delta={}", delta);
        }

        // Case 3: = row, violated. x0 = 5, at x0=2 => lhs=2, violated.
        // gap = 2 - 5 = -3, coeff=1, delta = -(-3)/1 = 3.
        // is_equality=true, is_int=true, coeff>0 => floor(3) = 3.
        // new_val = 2 + 3 = 5, within [0,10] => delta = 3.
        {
            let problem = S5Problem {
                col_lower: vec![0.0],
                col_upper: vec![10.0],
                row_lower: vec![5.0],
                row_upper: vec![5.0],
                col_cost: vec![1.0],
                offset: 0.0,
                integrality: vec![HighsVarType::kInteger],
                ar_start: vec![0, 1],
                ar_index: vec![0],
                ar_value: vec![1.0],
                csc: CscMatrix {
                    col_start: vec![0, 1],
                    col_row: vec![0],
                    col_val: vec![1.0],
                },
                uplocks: vec![0],
                downlocks: vec![0],
                feastol: 1e-6,
                epsilon: 1e-6,
                ncol: 1,
                nrow: 1,
                nnz: 1,
                incumbent: vec![0.0],
                binary: vec![0],
            };
            let mut ctx = WorkerCtx::new(&problem);
            ctx.solution[0] = 2.0;
            ctx.lhs[0] = 2.0;
            let delta = ctx.compute_tight_delta(0, 0, 1.0);
            assert!((delta - 3.0).abs() < 1e-12, "delta={}", delta);
        }
    }

    #[test]
    fn s5_lmcore_update_weights_smoothing() {
        // Find a seed where the draw < kSmoothProb (smoothing branch).
        let mut smooth_seed: Option<u64> = None;
        let mut strengthen_seed: Option<u64> = None;
        for seed in 1u64..1_000_000 {
            let mut probe = S5Rng::new(seed);
            let draw = s5_uniform_real(&mut probe, 0.0, 1.0);
            if draw < kSmoothProb && smooth_seed.is_none() {
                smooth_seed = Some(seed);
            }
            if draw >= kSmoothProb && strengthen_seed.is_none() {
                strengthen_seed = Some(seed);
            }
            if smooth_seed.is_some() && strengthen_seed.is_some() {
                break;
            }
        }
        // Both branches must be exercised: a silent skip would hide a regression.
        assert!(
            smooth_seed.is_some(),
            "no seed below 1e6 draws < kSmoothProb"
        );
        assert!(
            strengthen_seed.is_some(),
            "no seed below 1e6 draws >= kSmoothProb"
        );

        // Test smoothing branch: draw < kSmoothProb.
        if let Some(seed) = smooth_seed {
            let problem = make_problem_2x2();
            let mut ctx = WorkerCtx::new(&problem);
            ctx.obj_weight = 5;
            ctx.weight[0] = 3;
            ctx.weight[1] = 2;
            // Put row 0 in satisfied, row 1 in violated.
            ctx.satisfied.add(0);
            ctx.violated.add(1);
            ctx.current_obj = 10.0;
            let mut rng = S5Rng::new(seed);
            // best_feasible=true, best_obj=20.0 => obj_better = true (10 < 20)
            ctx.update_weights(&mut rng, false, true, 20.0);
            // Smoothing: obj_better && obj_weight>1 => obj_weight -= 1 => 4
            assert_eq!(ctx.obj_weight, 4, "obj_weight after smoothing");
            // satisfied row 0: weight[0] > 1 => weight[0] -= 1 => 2
            assert_eq!(ctx.weight[0], 2, "weight[0] after smoothing");
            // violated row 1 is not touched in smoothing
            assert_eq!(ctx.weight[1], 2, "weight[1] unchanged in smoothing");
        }

        // Test strengthen branch: draw >= kSmoothProb, is_feasible=true.
        if let Some(seed) = strengthen_seed {
            let problem = make_problem_2x2();
            let mut ctx = WorkerCtx::new(&problem);
            ctx.obj_weight = 5;
            ctx.weight[0] = 3;
            ctx.weight[1] = 2;
            ctx.violated.add(1);
            let mut rng = S5Rng::new(seed);
            ctx.update_weights(&mut rng, true, false, 0.0);
            // Strengthen + feasible: obj_weight += 1 => 6
            assert_eq!(ctx.obj_weight, 6, "obj_weight after strengthen+feasible");
            // weights unchanged
            assert_eq!(ctx.weight[0], 3);
            assert_eq!(ctx.weight[1], 2);
        }

        // Test strengthen branch: draw >= kSmoothProb, is_feasible=false.
        if let Some(seed) = strengthen_seed {
            let problem = make_problem_2x2();
            let mut ctx = WorkerCtx::new(&problem);
            ctx.obj_weight = 5;
            ctx.weight[0] = 3;
            ctx.weight[1] = 2;
            ctx.violated.add(0);
            ctx.violated.add(1);
            let mut rng = S5Rng::new(seed);
            ctx.update_weights(&mut rng, false, false, 0.0);
            // Strengthen + infeasible: violated rows get weight += 1
            assert_eq!(
                ctx.obj_weight, 5,
                "obj_weight unchanged in strengthen+infeasible"
            );
            assert_eq!(ctx.weight[0], 4, "weight[0] incremented");
            assert_eq!(ctx.weight[1], 3, "weight[1] incremented");
        }
    }
}
