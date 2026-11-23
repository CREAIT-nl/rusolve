//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! local_mip_search.cpp — candidate scoring, aspiration,
//! breakthrough delta, select_best_from_batch, select_lift_move,
//! infeasible_step.

use super::rng::{S5Rng, row_violation, s5_uniform_real};
use super::*;
use crate::prelude::*;
use crate::util::highs_utils::cpp_min;

// ---------------------------------------------------------------------------
// libstdc++ heap operations (bits/stl_heap.h) for __partial_sort
// ---------------------------------------------------------------------------

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__parent (the `(__holeIndex - 1) / 2` of __push_heap, C++ signed division)
fn heap_parent(hole: usize) -> usize {
    if hole == 0 { 0 } else { (hole - 1) / 2 }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__push_heap (lines 135-147)
/// `_Tp __value` is held by the caller and moved in as `value`; `first` is the
/// heap base and `hole`/`top` are distances from it, as in C++.
fn push_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    top: usize,
    value: T,
    comp: &mut F,
) {
    let mut parent = heap_parent(hole);
    while hole > top && comp(&v[first + parent], &value) {
        v[first + hole] = v[first + parent].clone();
        hole = parent;
        parent = heap_parent(hole);
    }
    v[first + hole] = value;
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__adjust_heap (lines 224-249)
fn adjust_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    len: usize,
    value: T,
    comp: &mut F,
) {
    let top = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if comp(&v[first + second], &v[first + second - 1]) {
            second -= 1;
        }
        v[first + hole] = v[first + second].clone();
        hole = second;
    }
    if len.is_multiple_of(2) && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[first + hole] = v[first + second - 1].clone();
        hole = second - 1;
    }
    push_heap(v, first, hole, top, value, comp);
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__pop_heap (lines 254-265)
fn pop_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    result: usize,
    comp: &mut F,
) {
    let value = v[result].clone();
    v[result] = v[first].clone();
    adjust_heap(v, first, 0, last - first, value, comp);
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__make_heap (lines 340-362)
fn make_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    if last - first < 2 {
        return;
    }
    let len = last - first;
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[first + parent].clone();
        adjust_heap(v, first, parent, len, value, comp);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__sort_heap (lines 419-427)
fn sort_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    comp: &mut F,
) {
    while last - first > 1 {
        last -= 1;
        pop_heap(v, first, last, last, comp);
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__heap_select (lines 1631-1639)
fn heap_select<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    middle: usize,
    last: usize,
    comp: &mut F,
) {
    make_heap(v, first, middle, comp);
    for i in middle..last {
        if comp(&v[i], &v[first]) {
            pop_heap(v, first, middle, i, comp);
        }
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__partial_sort (lines 1650-1656)
/// `std::partial_sort(first, middle, last, comp)` = `__heap_select(first, middle, last, comp)`
/// then `__sort_heap(first, middle, comp)`.
fn s5_partial_sort<T: Clone, F: FnMut(&T, &T) -> bool>(v: &mut [T], middle: usize, comp: &mut F) {
    let last = v.len();
    if middle >= last {
        return;
    }
    heap_select(v, 0, middle, last, comp);
    sort_heap(v, 0, middle, comp);
}

// ---------------------------------------------------------------------------
// local_mip_search.cpp
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): local_mip_search.cpp::append_candidate (lines 22-28)
fn append_candidate(ctx: &mut WorkerCtx, batch: &mut Vec<BatchCand>, j: HighsInt, delta: f64) {
    let new_val = ctx.clamp_and_round(j, ctx.solution[j as usize] + delta);
    if (new_val - ctx.solution[j as usize]).abs() < kEpsZero {
        return;
    }
    batch.push(BatchCand {
        var_idx: j,
        new_val,
    });
}

/// mip-heuristics (36f7d9d): local_mip_search.cpp::compute_candidate_scores (lines 40-152)
pub(crate) fn compute_candidate_scores(
    ctx: &mut WorkerCtx,
    j: HighsInt,
    new_val: f64,
    best_feasible: bool,
    best_obj: f64,
) -> (f64, f64) {
    let old_val = ctx.solution[j as usize];
    let delta = new_val - old_val;
    if delta.abs() < kEpsZero {
        return (f64::NEG_INFINITY, 0.0);
    }

    let col_begin = ctx.csc.col_start[j as usize] as usize;
    let col_end = ctx.csc.col_start[(j + 1) as usize] as usize;
    ctx.effort += col_end - col_begin;
    ctx.ls_evals += 1;

    // Def 5: progress score for objective
    let obj_delta = ctx.col_cost[j as usize] * delta;
    let new_obj = ctx.current_obj + obj_delta;
    let eps = ctx.epsilon;
    let mut progress = 0.0;
    if new_obj < ctx.current_obj - eps {
        progress += ctx.obj_weight as f64; // objective improved
    } else if new_obj > ctx.current_obj + eps {
        progress -= ctx.obj_weight as f64; // objective worsened
    }

    // Def 8: breakthrough bonus (beats best-found solution)
    let mut bonus = 0.0;
    if best_feasible {
        if new_obj < best_obj - eps {
            bonus += ctx.obj_weight as f64;
        }
    }

    // Defs 6-7, 9-10: constraint progress + robustness
    //
    // rusolve (throughput): the C++ reads
    // the pre-move violation through `ViolCache` and indexes the row
    // arrays once per use; here the column is walked as slices, the row
    // bounds are read once per nonzero, and the pre-move violation is
    // recomputed with the same `row_violation`. `lhs` does not change
    // while a batch is scored, so the recomputed value is the cached one
    // bit for bit, and the scores (and the trajectory) are unchanged.
    let feastol = ctx.feastol;
    let col_rows = &ctx.csc.col_row[col_begin..col_end];
    let col_vals = &ctx.csc.col_val[col_begin..col_end];
    // Equal-length row slices let one bounds check cover all four reads.
    let nrow = ctx.lhs.len();
    let (lhs, lo, hi, weight) = (
        &ctx.lhs[..nrow],
        &ctx.row_lo[..nrow],
        &ctx.row_hi[..nrow],
        &ctx.weight[..nrow],
    );
    for (&i, &coeff) in col_rows.iter().zip(col_vals) {
        let iu = i as usize;
        let row_lo = lo[iu];
        let row_hi = hi[iu];
        let old_lhs = lhs[iu];
        let new_lhs = old_lhs + (coeff * delta);
        let old_viol = row_violation(old_lhs, row_lo, row_hi);
        let new_viol = row_violation(new_lhs, row_lo, row_hi);
        let w = weight[iu] as f64;

        // Def 6: constraint progress score.
        //
        // DELIBERATE DEVIATION FROM THE PAPER, MATCHING THE AUTHORS' CODE:
        // Def 6 awards +/- w(con_i)/2 for a constraint that stays violated
        // but improves/worsens, reserving full weight for a satisfy/violate
        // transition.  We award full weight in both cases, as the authors'
        // own implementation does (github.com/shaowei-cai-group/Local-MIP,
        // `explore_unsat.cpp`: `if (new_gap < pre_gap) score += con_weight;
        // else score -= con_weight;`).  Consequence either way: "satisfy a
        // constraint" and "reduce a violation" score identically, flattening
        // the preference Def 6 states.  Do NOT "fix" this toward the paper
        // without measuring — it would move us away from the reference
        // implementation.  (The reference also scales equality rows by 2x,
        // which we do not do; that part is a genuine gap.)
        // `feastol`, not `kScoreTol` or `kViolDeltaTol` (issue #148):
        // this is the same "is this row violated?" question
        // `WorkerCtx::is_violated` and the `violated`/`satisfied` sets
        // answer, just evaluated on a hypothetical (pre-move / post-move)
        // lhs rather than the worker's live state — it has to use the
        // one tolerance that governs violation everywhere else, or this
        // scoring pass would reintroduce the same kind of mismatch #148
        // fixed elsewhere.
        //
        // rusolve (throughput): the C++
        // branches below are written as boolean selects (`&`/`|` instead of
        // `&&`/`if`) so the loop does not mispredict on every nonzero. Each
        // row still adds exactly one of `+w`, `-w` or nothing to `progress`
        // and `w` or nothing to `bonus`, in column order; adding `0.0 * w`
        // for "nothing" is exact because both sums start at `+0.0` and only
        // ever hold integer values, so they are never `-0.0`.
        //
        // Rows violated before and after (the C++ `else if (was && now)`)
        // compare magnitudes with `kViolDeltaTol`, not `kScoreTol` and not
        // `feastol`: this asks whether the violation *magnitude* moved by
        // more than floating-point noise, not whether the row is violated
        // (already decided above) and not whether a candidate *score*
        // changed meaningfully. See `kViolDeltaTol`'s definition in
        // `local_mip_caches.h`.
        let was_viol = old_viol > feastol;
        let now_viol = new_viol > feastol;
        if was_viol && !now_viol {
            progress += w; // violated → satisfied
        } else if !was_viol && now_viol {
            progress -= w; // satisfied → violated
        } else if was_viol && now_viol {
            if new_viol < old_viol - kViolDeltaTol {
                progress += w; // still violated, improved
            } else if new_viol > old_viol + kViolDeltaTol {
                progress -= w; // still violated, worsened
            }
        }
        // Def 9: robustness bonus — only for transitions into strictly
        // satisfied (was violated or tight, now strictly interior).
        //
        // DELIBERATE DEVIATION FROM THE PAPER, MATCHING THE AUTHORS' CODE:
        // Def 9 awards w(con_i) for *every* constraint that ends strictly
        // satisfied.  Summed over all rows and cancelling the rows this
        // column does not touch, the paper's per-column form is
        // sum_i w_i * ([new strict] - [old strict]) — i.e. it also
        // *penalises* losing strict satisfaction.  We implement only the
        // positive term.  The authors' implementation is likewise
        // transition-shaped (`if (!pre_sat && now_sat) score +=
        // scaled_con_weight * 2`).  Same caution as Def 6 above.
        if !now_viol {
            let old_strict = !was_viol
                && (row_hi >= kHighsInf || old_lhs < row_hi - feastol)
                && (row_lo <= -kHighsInf || old_lhs > row_lo + feastol);
            if !old_strict {
                let new_strict = (row_hi >= kHighsInf || new_lhs < row_hi - feastol)
                    && (row_lo <= -kHighsInf || new_lhs > row_lo + feastol);
                if new_strict {
                    bonus += w;
                }
            }
        }
    }

    (progress, bonus)
}

/// mip-heuristics (36f7d9d): local_mip_search.cpp::is_aspiration (lines 154-163)
pub(crate) fn is_aspiration(
    ctx: &WorkerCtx,
    j: HighsInt,
    new_val: f64,
    best_obj: f64,
    best_feasible: bool,
) -> bool {
    if !best_feasible {
        return false;
    }
    let delta = new_val - ctx.solution[j as usize];
    let obj_delta = ctx.col_cost[j as usize] * delta;
    let new_obj = ctx.current_obj + obj_delta;
    new_obj < best_obj - ctx.epsilon
}

/// mip-heuristics (36f7d9d): local_mip_search.cpp::compute_breakthrough_delta (lines 265-292)
pub(crate) fn compute_breakthrough_delta(
    ctx: &WorkerCtx,
    j: HighsInt,
    cur_obj: f64,
    best_obj: f64,
) -> f64 {
    if cur_obj < best_obj - ctx.epsilon {
        return 0.0;
    }

    let obj_coeff = ctx.col_cost[j as usize];
    if obj_coeff.abs() < kEpsZero {
        return 0.0;
    }

    let obj_gap = cur_obj - best_obj;

    // `+ ctx.epsilon`: paper Definition 2's own eps, dropped before this
    // fix -- see the comment above, including its continuous-variable
    // and trajectory consequences.
    let mut delta = -(obj_gap + ctx.epsilon) / obj_coeff;

    if ctx.is_int(j) {
        delta = if obj_coeff > 0.0 {
            delta.floor()
        } else {
            delta.ceil()
        };
    }
    let new_val = ctx.solution[j as usize] + delta;
    if new_val < ctx.col_lb[j as usize] || new_val > ctx.col_ub[j as usize] {
        delta = if obj_coeff > 0.0 {
            ctx.col_lb[j as usize] - ctx.solution[j as usize]
        } else {
            ctx.col_ub[j as usize] - ctx.solution[j as usize]
        };
    }
    delta
}

/// mip-heuristics (36f7d9d): local_mip_search.cpp::select_best_from_batch (lines 294-327)
pub(crate) fn select_best_from_batch(
    ctx: &mut WorkerCtx,
    batch: &mut Vec<BatchCand>,
    step: HighsInt,
    aspiration: bool,
    best_obj: f64,
    best_feasible: bool,
) -> Candidate {
    let mut best = Candidate::default();
    for cc in batch.iter() {
        let delta = cc.new_val - ctx.solution[cc.var_idx as usize];
        if delta.abs() < kEpsZero {
            continue;
        }

        if ctx.is_tabu(cc.var_idx, delta, step) {
            if !(aspiration && is_aspiration(ctx, cc.var_idx, cc.new_val, best_obj, best_feasible))
            {
                continue;
            }
        }

        let (prog, bon) =
            compute_candidate_scores(ctx, cc.var_idx, cc.new_val, best_feasible, best_obj);

        // `kScoreTol`, every use from here down in this file: comparing
        // two candidate scores for a numerically meaningful difference,
        // not a violation question — see its definition in
        // `local_mip_caches.h` (issue #148).
        if prog > best.score + kScoreTol {
            best = Candidate {
                var_idx: cc.var_idx,
                new_val: cc.new_val,
                score: prog,
                bonus: bon,
            };
        } else if prog > best.score - kScoreTol {
            if bon > best.bonus {
                best = Candidate {
                    var_idx: cc.var_idx,
                    new_val: cc.new_val,
                    score: prog,
                    bonus: bon,
                };
            }
        }
    }
    best
}

/// mip-heuristics (36f7d9d): local_mip_search.cpp::select_lift_move (lines 384-420)
pub(crate) fn select_lift_move(ctx: &mut WorkerCtx) -> Candidate {
    let mut best = Candidate::default();
    best.score = 0.0;

    let mut write: HighsInt = 0;
    // Compaction, not traversal: the body writes back into the same
    // vector at `write <= read` and never resizes it (the `resize`
    // below is what shortens the list), so the bound is
    // loop-invariant.  Hoisting it drops a size() load per iteration
    // of LocalMIP's per-restart loop and keeps modernize-loop-convert
    // from proposing a range-for that would hide the rewrite.
    let n_positive = ctx.lift.positive_list.len() as HighsInt;
    let mut read: HighsInt = 0;
    while read < n_positive {
        let j = ctx.lift.positive_list[read as usize];
        if !ctx.lift.in_positive[j as usize] {
            read += 1;
            continue;
        }
        ctx.lift.positive_list[write as usize] = j;
        write += 1;
        if ctx.lift.score[j as usize] <= best.score {
            read += 1;
            continue;
        }
        let lo = ctx.lift.lo[j as usize];
        let hi = ctx.lift.hi[j as usize];
        if lo > hi {
            read += 1;
            continue;
        }
        let mut target = if ctx.col_cost[j as usize] > 0.0 {
            lo
        } else {
            hi
        };
        target = ctx.clamp_and_round(j, target);
        if (target - ctx.solution[j as usize]).abs() < kEpsZero {
            read += 1;
            continue;
        }
        best = Candidate {
            var_idx: j,
            new_val: target,
            score: ctx.lift.score[j as usize],
            bonus: 0.0,
        };
        read += 1;
    }
    ctx.lift.positive_list.truncate(write as usize);

    best
}

/// Phase 1b of `infeasible_step`: one breakthrough candidate per costed
/// column, in `costed_vars` order (the breakthrough loop of upstream
/// `local_mip_search.cpp` infeasible_step).
fn append_breakthrough_candidates(
    ctx: &mut WorkerCtx,
    batch: &mut Vec<BatchCand>,
    costed_vars: &[HighsInt],
    best_objective: f64,
) {
    // rusolve (throughput): the C++
    // visits every costed column on every step, O(ncol) work the effort
    // counter never sees. `bt_live` holds the columns that can yield a
    // candidate (see `WorkerCtx::bt_can_yield`); visiting only those, in
    // the same ascending order, appends the same candidates. The skip is
    // exact only when the Definition 2 gap is non-negative or the delta
    // is zero for every column, so other states take the C++ loop.
    let gap_ok = ctx.current_obj < best_objective - ctx.epsilon
        || (ctx.current_obj - best_objective) + ctx.epsilon >= 0.0;
    if ctx.bt_live_valid && gap_ok && costed_vars.len() == ctx.bt_costed {
        #[cfg(debug_assertions)]
        let first = batch.len();
        for w in 0..ctx.bt_live.len() {
            let mut bits = ctx.bt_live[w];
            while bits != 0 {
                let j = (w * 64 + bits.trailing_zeros() as usize) as HighsInt;
                bits &= bits - 1;
                let delta = compute_breakthrough_delta(ctx, j, ctx.current_obj, best_objective);
                append_candidate(ctx, batch, j, delta);
            }
        }
        #[cfg(debug_assertions)]
        {
            let mut full = Vec::new();
            for &j in costed_vars.iter() {
                let delta = compute_breakthrough_delta(ctx, j, ctx.current_obj, best_objective);
                append_candidate(ctx, &mut full, j, delta);
            }
            debug_assert_eq!(&batch[first..], &full[..]);
        }
    } else {
        for &j in costed_vars.iter() {
            let delta = compute_breakthrough_delta(ctx, j, ctx.current_obj, best_objective);
            append_candidate(ctx, batch, j, delta);
        }
    }
}

/// mip-heuristics (36f7d9d): local_mip_search.cpp::infeasible_step (lines 428-674)
pub(crate) fn infeasible_step(
    ctx: &mut WorkerCtx,
    rng: &mut S5Rng,
    step: HighsInt,
    best_feasible: bool,
    best_objective: f64,
    costed_vars: &[HighsInt],
    binary_vars: &[HighsInt],
) -> Candidate {
    // Issue #129: this is now also entered from a *feasible* state, on a
    // failed lift (`LocalMipWorker::run_attempt`'s fall-through) -- not
    // only from a genuinely infeasible one.  `was_infeasible` gates the
    // lift cache's incremental dirty-tracking and the caller's next
    // `need_full_recheck` (both in `local_mip_worker.cpp`), so it must
    // record whether *this* entry was genuinely infeasible, not merely
    // that this function ran: forcing it to `true` unconditionally would
    // charge a full recheck on every plateau iteration reached through
    // the fall-through path, and would suppress dirty-marking for moves
    // that never left feasible territory. Identical to the old
    // unconditional `true` for the pre-existing (genuinely infeasible)
    // caller, since `ctx.violated` is non-empty there by construction.
    ctx.was_infeasible = !ctx.violated.empty();

    let mut batch = std::mem::take(&mut ctx.batch);
    let mut sampled = std::mem::take(&mut ctx.sampled);

    let result = infeasible_step_inner(
        ctx,
        rng,
        step,
        best_feasible,
        best_objective,
        costed_vars,
        binary_vars,
        &mut batch,
        &mut sampled,
    );

    ctx.batch = batch;
    ctx.sampled = sampled;
    result
}

fn infeasible_step_inner(
    ctx: &mut WorkerCtx,
    rng: &mut S5Rng,
    step: HighsInt,
    best_feasible: bool,
    best_objective: f64,
    costed_vars: &[HighsInt],
    binary_vars: &[HighsInt],
    batch: &mut Vec<BatchCand>,
    sampled: &mut Vec<WeightedCon>,
) -> Candidate {
    // --- Phase 1: BMS tight moves from violated constraints ---
    //
    // DELIBERATE DEVIATION FROM THE PAPER, MATCHING THE AUTHORS' CODE:
    // Algorithm 2 enumerates *all* violated constraints.  BMS (Best from
    // Multiple Selections) appears nowhere in the paper, but the authors'
    // implementation samples exactly this way
    // (github.com/shaowei-cai-group/Local-MIP, `explore_unsat.cpp`:
    // `sample_idxs(..., m_bms_con, ...)` and `sample_op(m_bms_op, ...)`).
    // The same applies to the analogous caps in the satisfied-constraint
    // and boolean-flip phases below.  Reference-faithful, not
    // paper-faithful; see mip-heuristics' PARAMETERS.md `kBmsConstraints`.
    let num_to_sample = cpp_min(kBmsConstraints * 3, ctx.violated.size());
    let num_to_keep = cpp_min(kBmsConstraints, ctx.violated.size());

    sampled.clear();
    if num_to_sample == ctx.violated.size() {
        for &ci in ctx.violated.elements.iter() {
            sampled.push(WeightedCon {
                ci,
                w: ctx.weight[ci as usize],
            });
        }
    } else {
        let mut s: HighsInt = 0;
        while s < num_to_sample {
            let idx = (rng.next_u64() % (ctx.violated.size() as u64)) as HighsInt;
            let ci = ctx.violated.elements[idx as usize];
            sampled.push(WeightedCon {
                ci,
                w: ctx.weight[ci as usize],
            });
            s += 1;
        }
    }

    if sampled.len() > num_to_keep as usize {
        s5_partial_sort(
            sampled,
            num_to_keep as usize,
            &mut |a: &WeightedCon, b: &WeightedCon| a.w > b.w,
        );
        sampled.truncate(num_to_keep as usize);
    }

    batch.clear();
    let mut budget_remaining = kBmsBudget;

    for si in 0..sampled.len() {
        let ci = sampled[si].ci;
        // (void)w;
        if budget_remaining <= 0 {
            break;
        }
        let row = ctx.tight_row(ci);
        let mut k = ctx.ar_start[ci as usize];
        while k < ctx.ar_start[(ci + 1) as usize] && budget_remaining > 0 {
            let j = ctx.ar_index[k as usize];
            budget_remaining -= 1;
            let delta = ctx.tight_delta_in_row(&row, j, ctx.ar_value[k as usize]);
            append_candidate(ctx, batch, j, delta);
            k += 1;
        }
    }

    // --- Phase 1b: Breakthrough moves (only post-feasible, Alg 2 line 5-6) ---
    //
    // Definition 2's other precondition, `obj(s) >= obj(s*)`, is enforced
    // inside `compute_breakthrough_delta` rather than here (issue #150 --
    // the authors' code spells it at this loop instead; see that
    // function's comment for why the operator is the place for it).
    if best_feasible {
        append_breakthrough_candidates(ctx, batch, costed_vars, best_objective);
    }

    let mut cand = select_best_from_batch(ctx, batch, step, true, best_objective, best_feasible);

    // If positive candidate found, done (Alg 2 lines 1-6)
    if cand.var_idx != -1 && cand.score > kScoreTol {
        return cand;
    }

    // --- Phase 2: MTM in satisfied constraints (Alg 2 lines 7-8) ---
    if !ctx.satisfied.empty() {
        batch.clear();
        let num_sat_sample = cpp_min(kBmsSatCon, ctx.satisfied.size());
        let mut sat_budget = kBmsSatBudget;
        let mut s: HighsInt = 0;
        while s < num_sat_sample && sat_budget > 0 {
            let ci =
                ctx.satisfied.elements[(rng.next_u64() % (ctx.satisfied.size() as u64)) as usize];
            let row = ctx.tight_row(ci);
            let mut k = ctx.ar_start[ci as usize];
            while k < ctx.ar_start[(ci + 1) as usize] && sat_budget > 0 {
                let j = ctx.ar_index[k as usize];
                sat_budget -= 1;
                let delta = ctx.tight_delta_in_row(&row, j, ctx.ar_value[k as usize]);
                append_candidate(ctx, batch, j, delta);
                k += 1;
            }
            s += 1;
        }
        let sat_cand =
            select_best_from_batch(ctx, batch, step, false, best_objective, best_feasible);
        if sat_cand.var_idx != -1 && sat_cand.score > cand.score + kScoreTol {
            cand = sat_cand;
        }
    }

    if cand.var_idx != -1 && cand.score > kScoreTol {
        return cand;
    }

    // --- Phase 3: Boolean flip (Alg 2 lines 9-11) ---
    if !binary_vars.is_empty() {
        batch.clear();
        let nbinary = binary_vars.len() as HighsInt;
        let offset = (rng.next_u64() % (nbinary as u64)) as HighsInt;
        let mut idx: HighsInt = 0;
        while idx < nbinary && idx < kBoolFlipBudget {
            let j = binary_vars[((offset + idx) % nbinary) as usize];
            let new_val = if ctx.solution[j as usize] < 0.5 {
                1.0
            } else {
                0.0
            };
            if (new_val - ctx.solution[j as usize]).abs() < kEpsZero {
                idx += 1;
                continue;
            }
            batch.push(BatchCand {
                var_idx: j,
                new_val,
            });
            idx += 1;
        }
        if !batch.is_empty() {
            let flip_cand =
                select_best_from_batch(ctx, batch, step, true, best_objective, best_feasible);
            if flip_cand.var_idx != -1 && flip_cand.score > cand.score + kScoreTol {
                cand = flip_cand;
            }
        }
    }

    if cand.var_idx != -1 && cand.score > kScoreTol {
        return cand;
    }

    // rusolve GUB swap (`gub_swap.rs`), phase 3b: GUB swap /
    // ejection-chain compound moves at a single-move local minimum. `ctx.gub`
    // is None unless mip_ls_gub_swap is on, so this is skipped otherwise.
    if ctx.gub.is_some() {
        if let Some(compound) =
            super::gub_swap::gub_phase(ctx, rng, step, best_feasible, best_objective)
        {
            return compound;
        }
    }

    // --- Phase 4: Weight update + random constraint fallback (Alg 2 lines 12-14) ---
    //
    // `is_feasible` reads the CURRENT state (`ctx.violated`, unchanged by
    // anything above -- phases 1-3 only score candidates, never apply
    // one), not a hardcoded `false`: issue #129 also reaches this phase
    // from a feasible worker whose lift just failed, and the paper's
    // PAWS-style weight scheme strengthens `w(obj)` at a feasible local
    // optimum, `w(coni)` at an infeasible one (Algorithm 1 lines 4-7 vs.
    // Algorithm 2 lines 12-14 -- the same update, two triggers). Identical
    // to the old hardcoded `false` for the pre-existing (genuinely
    // infeasible) caller.
    ctx.update_weights(rng, ctx.violated.empty(), best_feasible, best_objective);

    if !ctx.violated.empty() {
        batch.clear();
        let ci = ctx.violated.elements[(rng.next_u64() % (ctx.violated.size() as u64)) as usize];
        let row = ctx.tight_row(ci);
        let mut k = ctx.ar_start[ci as usize];
        while k < ctx.ar_start[(ci + 1) as usize] {
            let j = ctx.ar_index[k as usize];
            let delta = ctx.tight_delta_in_row(&row, j, ctx.ar_value[k as usize]);
            append_candidate(ctx, batch, j, delta);
            k += 1;
        }
        // Breakthrough candidates already scored in Phase 1; skip re-scoring.
        let fallback =
            select_best_from_batch(ctx, batch, step, false, best_objective, best_feasible);
        if fallback.var_idx != -1
            && (cand.var_idx == -1
                || fallback.score > cand.score + kScoreTol
                || (fallback.score > cand.score - kScoreTol && fallback.bonus > cand.bonus))
        {
            cand = fallback;
        }
    }

    if cand.var_idx != -1 {
        return cand;
    }

    // --- Phase 5: Perturbation (our addition, last resort) ---
    if !ctx.violated.empty() {
        let ci = ctx.violated.elements[(rng.next_u64() % (ctx.violated.size() as u64)) as usize];
        let row_len = ctx.ar_start[(ci + 1) as usize] - ctx.ar_start[ci as usize];
        if row_len > 0 {
            let k = ctx.ar_start[ci as usize] + (rng.next_u64() % (row_len as u64)) as HighsInt;
            let j = ctx.ar_index[k as usize];
            let new_val = if ctx.is_binary(j) {
                if ctx.solution[j as usize] < 0.5 {
                    1.0
                } else {
                    0.0
                }
            } else if ctx.is_int(j) {
                // `rng() % 2 == 0` (local_mip_search.cpp:616)
                let dir: HighsInt = if rng.next_u64().is_multiple_of(2) {
                    1
                } else {
                    -1
                };
                ctx.clamp_and_round(j, ctx.solution[j as usize] + dir as f64)
            } else {
                let range = cpp_min(ctx.col_ub[j as usize], ctx.col_lb[j as usize] + 1e6)
                    - ctx.col_lb[j as usize];
                let perturbation = s5_uniform_real(rng, -0.1 * range, 0.1 * range);
                ctx.clamp_and_round(j, ctx.solution[j as usize] + perturbation)
            };
            if (new_val - ctx.solution[j as usize]).abs() > kEpsZero {
                let (prog, bon) =
                    compute_candidate_scores(ctx, j, new_val, best_feasible, best_objective);
                cand = Candidate {
                    var_idx: j,
                    new_val,
                    score: prog,
                    bonus: bon,
                };
            }
        }
    }

    if cand.var_idx != -1 {
        return cand;
    }

    // --- Phase 6: Easy moves (our addition) ---
    {
        batch.clear();
        let num_easy = cpp_min(kEasyBudget, ctx.ncol);
        let mut s: HighsInt = 0;
        while s < num_easy {
            let j = (rng.next_u64() % (ctx.ncol as u64)) as HighsInt;
            let target = if ctx.col_lb[j as usize] > 0.0 {
                ctx.col_lb[j as usize]
            } else if ctx.col_ub[j as usize] < 0.0 {
                ctx.col_ub[j as usize]
            } else {
                0.0
            };
            append_candidate(ctx, batch, j, target - ctx.solution[j as usize]);
            // Try: toward lower bound
            if ctx.col_lb[j as usize] > -1e15 && ctx.col_lb[j as usize] < 0.0 {
                append_candidate(
                    ctx,
                    batch,
                    j,
                    ctx.col_lb[j as usize] - ctx.solution[j as usize],
                );
            }
            // Try: toward upper bound
            if ctx.col_ub[j as usize] < 1e15 && ctx.col_ub[j as usize] > 0.0 {
                append_candidate(
                    ctx,
                    batch,
                    j,
                    ctx.col_ub[j as usize] - ctx.solution[j as usize],
                );
            }
            // Try: midpoint for continuous
            if !ctx.is_int(j) && ctx.col_lb[j as usize] > -1e15 && ctx.col_ub[j as usize] < 1e15 {
                append_candidate(
                    ctx,
                    batch,
                    j,
                    ((ctx.col_lb[j as usize] + ctx.col_ub[j as usize]) * 0.5)
                        - ctx.solution[j as usize],
                );
            }
            s += 1;
        }
        let easy_cand =
            select_best_from_batch(ctx, batch, step, false, best_objective, best_feasible);
        if easy_cand.var_idx != -1 {
            cand = easy_cand;
        }
    }

    cand
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

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

    fn make_ctx(problem: &S5Problem) -> WorkerCtx<'_> {
        WorkerCtx::new(problem)
    }

    #[test]
    fn s5_lmsearch_candidate_scores_progress_and_bonus() {
        // 2-col problem: x0+x1<=3, x0+x1>=1, cost [1,1], integer, bounds [0,10].
        // solution = [2, 0], lhs = [2, 2]. Row 0: 2<=3 satisfied. Row 1: 2>=1 satisfied.
        // violated empty, satisfied = {0,1}. current_obj = 1*2 + 1*0 = 2.
        // Test: j=0, new_val=1 (delta=-1). new_obj = 2 + 1*(-1) = 1 < 2 - eps => progress += obj_weight=1.
        // Row 0 (<=3): old_lhs=2, new_lhs=2+1*(-1)=1. old_viol=0 (satisfied).
        //   new_viol=0 (1<=3). was_viol=false, now_viol=false => no progress change.
        //   Def 9: now_viol=false. old_strict: !was_viol && (3>=inf? no, 2<3-eps? yes) && (-inf<=-inf? yes) => true.
        //   old_strict=true => skip bonus.
        // Row 1 (>=1): old_lhs=2, new_lhs=1. old_viol=0 (2>=1). new_viol=0 (1>=1).
        //   was_viol=false, now_viol=false => no progress change.
        //   Def 9: now_viol=false. old_strict: !was_viol && (inf>=inf? yes) && (1<=-inf? no, 2>1+eps? yes) => true.
        //   old_strict=true => skip bonus.
        // progress = 1 (obj improved). bonus = 0.
        let problem = make_problem_2x2();
        let mut ctx = make_ctx(&problem);
        ctx.solution = vec![2.0, 0.0];
        ctx.lhs = vec![2.0, 2.0];
        ctx.current_obj = 2.0;
        ctx.weight = vec![1, 1];
        ctx.obj_weight = 1;
        ctx.violated = IndexedSet::new(2);
        ctx.satisfied = IndexedSet::new(2);
        ctx.satisfied.add(0);
        ctx.satisfied.add(1);

        let (prog, bon) = compute_candidate_scores(&mut ctx, 0, 1.0, false, 0.0);
        // progress = obj_weight(1) = 1.0 (objective improved, no constraint transitions)
        assert!((prog - 1.0).abs() < 1e-12, "prog={}", prog);
        assert!((bon - 0.0).abs() < 1e-12, "bon={}", bon);
    }

    #[test]
    fn s5_lmsearch_aspiration_requires_best_feasible() {
        let problem = make_problem_2x2();
        let ctx = make_ctx(&problem);
        // best_feasible=false => always false
        assert!(!is_aspiration(&ctx, 0, 1.0, 5.0, false));
        // best_feasible=true, solution=[2,0], current_obj=2, j=0, new_val=1 => delta=-1, obj_delta=-1, new_obj=1.
        // 1 < 5 - eps => true
        let mut ctx2 = make_ctx(&problem);
        ctx2.solution = vec![2.0, 0.0];
        ctx2.current_obj = 2.0;
        assert!(is_aspiration(&ctx2, 0, 1.0, 5.0, true));
        // new_obj=1 is NOT < 0 - eps => false
        assert!(!is_aspiration(&ctx2, 0, 1.0, 0.0, true));
    }

    #[test]
    fn s5_lmsearch_breakthrough_delta_definition2() {
        // col_cost[0]=1, col_lb[0]=0, col_ub[0]=10, solution[0]=2, is_int=true.
        // epsilon = 1e-6.
        let problem = make_problem_2x2();
        let mut ctx = make_ctx(&problem);
        ctx.solution = vec![2.0, 0.0];

        // (cur_obj=6, best_obj=5): obj_gap=1, delta=-(1+eps)/1=-1-eps, floor=-2.
        // new_val=2+(-2)=0, in [0,10], so delta=-2.
        let d = compute_breakthrough_delta(&ctx, 0, 6.0, 5.0);
        assert!((d - (-2.0)).abs() < 1e-12, "d={}", d);

        // (cur_obj=5, best_obj=5): obj_gap=0, delta=-(0+eps)/1=-eps, floor=-1.
        // new_val=2+(-1)=1, in [0,10], so delta=-1.
        let d = compute_breakthrough_delta(&ctx, 0, 5.0, 5.0);
        assert!((d - (-1.0)).abs() < 1e-12, "d={}", d);

        // (cur_obj=3, best_obj=5): cur_obj < best_obj - eps => return 0.
        let d = compute_breakthrough_delta(&ctx, 0, 3.0, 5.0);
        assert!((d - 0.0).abs() < 1e-12, "d={}", d);
    }

    #[test]
    fn s5_lmsearch_select_best_bms_sampling() {
        // 2-col problem, solution=[2,0], current_obj=2, no tabu, no aspiration.
        // batch: [{0, 1.0}, {1, 3.0}]
        // j=0: delta=1-2=-1, obj_delta=-1, new_obj=1<2-eps => progress+=1.
        //   Row 0 (<=3): old_lhs=2, new_lhs=1. old_viol=0, new_viol=0. was_viol=false, now_viol=false.
        //     Def 9: now_viol=false. old_strict: !was_viol && (3>=inf? no, 2<3-eps? yes) && (-inf<=-inf? yes) => true.
        //     old_strict=true => skip.
        //   Row 1 (>=1): old_lhs=2, new_lhs=1. old_viol=0, new_viol=0. was_viol=false, now_viol=false.
        //     Def 9: now_viol=false. old_strict: !was_viol && (inf>=inf? yes) && (1<=-inf? no, 2>1+eps? yes) => true.
        //     old_strict=true => skip.
        //   progress=1, bonus=0.
        // j=1: delta=3-0=3, obj_delta=3, new_obj=5>2+eps => progress-=1.
        //   Row 0 (<=3): old_lhs=2, new_lhs=5. old_viol=0, new_viol=5-3=2>feastol.
        //     was_viol=false, now_viol=true => progress-=w=1. progress=-1-1=-2.
        //   Row 1 (>=1): old_lhs=2, new_lhs=5. old_viol=0, new_viol=0. was_viol=false, now_viol=false.
        //     Def 9: now_viol=false. old_strict: !was_viol && (inf>=inf? yes) && (1<=-inf? no, 2>1+eps? yes) => true.
        //     old_strict=true => skip.
        //   progress=-2, bonus=0.
        // best starts at score=-inf. j=0: prog=1 > -inf+tol => best={0,1,1,0}.
        // j=1: prog=-2 > 1-tol? No. So best stays j=0.
        let problem = make_problem_2x2();
        let mut ctx = make_ctx(&problem);
        ctx.solution = vec![2.0, 0.0];
        ctx.lhs = vec![2.0, 2.0];
        ctx.current_obj = 2.0;
        ctx.weight = vec![1, 1];
        ctx.obj_weight = 1;
        ctx.violated = IndexedSet::new(2);
        ctx.satisfied = IndexedSet::new(2);
        ctx.satisfied.add(0);
        ctx.satisfied.add(1);
        ctx.tabu_inc_until = vec![0; 2];
        ctx.tabu_dec_until = vec![0; 2];

        let mut batch = vec![
            BatchCand {
                var_idx: 0,
                new_val: 1.0,
            },
            BatchCand {
                var_idx: 1,
                new_val: 3.0,
            },
        ];
        let best = select_best_from_batch(&mut ctx, &mut batch, 0, false, 0.0, false);
        assert_eq!(best.var_idx, 0);
        assert!((best.score - 1.0).abs() < 1e-12);
    }

    #[test]
    fn s5_lmsearch_infeasible_step_hand_traced() {
        // 2-col problem: x0+x1<=3, x0+x1>=1, cost [1,1], integer, bounds [0,10].
        // solution=[0,0], lhs=[0,0]. Row 0: 0<=3 satisfied. Row 1: 0<1 violated.
        // violated={1}, satisfied={0}. current_obj=0, best_feasible=false, best_objective=0.
        // Phase 1: num_to_sample=min(36,1)=1, num_to_keep=min(12,1)=1.
        //   sampled: num_to_sample==violated.size()=1, push {1, weight[1]=1}.
        //   No partial_sort needed.
        //   batch: budget=2250. Row 1 (>=1): ar_start[1]=2, ar_start[2]=4.
        //     k=2: j=0, delta=compute_tight_delta(1,0,1.0). row_lo=1, row_hi=inf, lhs=0, violated.
        //       tight_delta = (1-0)/1 = 1. clamp_and_round(0, 0+1)=1. |1-0|>=eps => push {0,1}.
        //     k=3: j=1, delta=compute_tight_delta(1,1,1.0) = (1-0)/1=1. clamp_and_round(1,0+1)=1.
        //       |1-0|>=eps => push {1,1}.
        //   Phase 1b: best_feasible=false, skip.
        //   select_best_from_batch: j=0: delta=1, obj_delta=1, new_obj=1>0+eps => progress-=1.
        //     Row 0 (<=3): old_lhs=0, new_lhs=1. old_viol=0, new_viol=0. was_viol=false, now_viol=false.
        //       Def 9: now_viol=false. old_strict: !was_viol && (3>=inf? no, 0<3-eps? yes) && (-inf<=-inf? yes) => true.
        //       old_strict=true => skip.
        //     Row 1 (>=1): old_lhs=0, new_lhs=1. old_viol=1-0=1>feastol, new_viol=0 (1>=1).
        //       was_viol=true, now_viol=false => progress+=w=1. progress=-1+1=0.
        //       Def 9: now_viol=false. old_strict: was_viol=true => false.
        //       new_strict: (inf>=inf? yes) && (1<=-inf? no, 1>1+eps? no) => false. skip.
        //     progress=0, bonus=0.
        //   j=1: delta=1, obj_delta=1, new_obj=1>eps => progress-=1.
        //     Row 0 (<=3): old_lhs=0, new_lhs=1. old_viol=0, new_viol=0. was_viol=false, now_viol=false.
        //       Def 9: old_strict=true => skip.
        //     Row 1 (>=1): old_lhs=0, new_lhs=1. old_viol=1>feastol, new_viol=0.
        //       was_viol=true, now_viol=false => progress+=1. progress=0.
        //       Def 9: old_strict=false. new_strict: (inf>=inf) && (1>1+eps? no) => false. skip.
        //     progress=0, bonus=0.
        //   Both have prog=0, bon=0. First wins: best={0,1,0,0}. Second: prog=0 > 0-tol? yes. bon=0 > 0? No.
        //   cand={0,1,0,0}. score=0 > kScoreTol? No.
        // Phase 2: satisfied={0}, not empty. num_sat_sample=min(1,1)=1, sat_budget=80.
        //   s=0: rng draw for satisfied index. satisfied.size()=1, so idx=0, ci=satisfied.elements[0]=0.
        //   Row 0 (<=3): ar_start[0]=0, ar_start[1]=2.
        //     k=0: j=0, delta=compute_tight_delta(0,0,1.0). row_lo=-inf, row_hi=3, lhs=0, satisfied.
        //       tight_delta for satisfied row: need to check compute_tight_delta logic.
        //       For a satisfied row, tight_delta moves toward the bound. row_hi=3, lhs=0, coeff=1.
        //       delta = (3-0)/1 = 3? Or depends on implementation. Let's assume it gives some value.
        //       Actually compute_tight_delta for satisfied rows typically gives 0 or a small move.
        //       The exact value depends on lm_core.rs implementation. For this test, we just check
        //       that the function completes and returns a valid candidate.
        //     k=1: j=1, similar.
        //   sat_cand from select_best. If sat_cand.score > cand.score + tol, update.
        // Phase 3: binary_vars empty, skip.
        // Phase 4: update_weights. violated not empty: pick random row from violated (only row 1).
        //   rng draw: violated.size()=1, idx=0, ci=1.
        //   batch: k=2: j=0, delta=1, push {0,1}. k=3: j=1, delta=1, push {1,1}.
        //   select_best: same as Phase 1, fallback={0,1,0,0}.
        //   fallback.var_idx!=-1 && (cand.var_idx==-1 || fallback.score>cand.score+tol ||
        //     (fallback.score>cand.score-tol && fallback.bonus>cand.bonus))
        //   cand.var_idx=0!=-1. fallback.score=0 > 0+tol? No. fallback.score=0 > 0-tol? Yes.
        //   fallback.bonus=0 > cand.bonus=0? No. So cand unchanged.
        //   cand.var_idx!=-1 => return cand={0,1,0,0}.
        let problem = make_problem_2x2();
        let mut ctx = make_ctx(&problem);
        ctx.solution = vec![0.0, 0.0];
        ctx.lhs = vec![0.0, 0.0];
        ctx.current_obj = 0.0;
        ctx.weight = vec![1, 1];
        ctx.obj_weight = 1;
        ctx.violated = IndexedSet::new(2);
        ctx.violated.add(1);
        ctx.satisfied = IndexedSet::new(2);
        ctx.satisfied.add(0);
        ctx.tabu_inc_until = vec![0; 2];
        ctx.tabu_dec_until = vec![0; 2];
        ctx.lift = LiftCache::default();

        let mut rng = S5Rng::new(42);
        let binary_vars: Vec<HighsInt> = Vec::new();
        let costed_vars: Vec<HighsInt> = vec![0, 1];

        let cand = infeasible_step(
            &mut ctx,
            &mut rng,
            0,
            false,
            0.0,
            &costed_vars,
            &binary_vars,
        );
        assert_eq!(cand.var_idx, 0);
        assert!((cand.new_val - 1.0).abs() < 1e-12);
    }

    #[test]
    fn s5_lmsearch_partial_sort_ties_libstdcxx() {
        // Oracle: g++ 13.3 libstdc++ `std::partial_sort(v.begin(), v.begin() + 7, v.end(),
        // [](a, b) { return a.w > b.w; })` then `v.resize(7)` on these 20 (ci, w) pairs
        // prints the ci order 7,4,11,14,2,10,0 (tie-dependent heap order).
        let ws: [u64; 20] = [3, 1, 3, 2, 3, 1, 2, 3, 1, 2, 3, 3, 2, 1, 3, 2, 1, 3, 2, 3];
        let mut v: Vec<WeightedCon> = ws
            .iter()
            .enumerate()
            .map(|(i, &w)| WeightedCon {
                ci: i as HighsInt,
                w,
            })
            .collect();
        s5_partial_sort(&mut v, 7, &mut |a: &WeightedCon, b: &WeightedCon| a.w > b.w);
        v.truncate(7);
        let order: Vec<HighsInt> = v.iter().map(|x| x.ci).collect();
        assert_eq!(order, vec![7, 4, 11, 14, 2, 10, 0]);
    }

    /// n rows `x_i >= 1`, one column each (cost 0, integer, [0,10]); all rows violated at x = 0.
    fn make_problem_diag(n: HighsInt) -> S5Problem {
        let nu = n as usize;
        let idx: Vec<HighsInt> = (0..n).collect();
        let start: Vec<HighsInt> = (0..=n).collect();
        S5Problem {
            col_lower: vec![0.0; nu],
            col_upper: vec![10.0; nu],
            row_lower: vec![1.0; nu],
            row_upper: vec![kHighsInf; nu],
            col_cost: vec![0.0; nu],
            offset: 0.0,
            integrality: vec![HighsVarType::kInteger; nu],
            ar_start: start.clone(),
            ar_index: idx.clone(),
            ar_value: vec![1.0; nu],
            csc: CscMatrix {
                col_start: start,
                col_row: idx,
                col_val: vec![1.0; nu],
            },
            uplocks: vec![0; nu],
            downlocks: vec![0; nu],
            feastol: 1e-6,
            epsilon: 1e-9,
            ncol: n,
            nrow: n,
            nnz: nu,
            incumbent: vec![0.0; nu],
            binary: vec![0; nu],
        }
    }

    #[test]
    fn s5_lmsearch_infeasible_step_bms_sampling_draws() {
        // 40 violated rows > kBmsConstraints * 3 = 36, so Phase 1 draws 36 times
        // `rng() % 40` (local_mip_search.cpp:469); row i has distinct weight i + 1.
        // Every sampled row's single column moves 0 -> 1 (tight delta 1), satisfying
        // the row: progress = w > kScoreTol, so Phase 1 returns (no further draws).
        // Expected: the kept 12 are the 12 largest sampled weights (descending; equal
        // weights are the same row, so ties are indistinguishable), the candidate is
        // the heaviest sampled row's column, and exactly 36 draws were consumed.
        let n: HighsInt = 40;
        let problem = make_problem_diag(n);
        let mut ctx = make_ctx(&problem);
        ctx.weight = (1..=n as u64).collect();
        for i in 0..n {
            ctx.violated.add(i);
        }
        let seed = 1039u64;
        let mut replay = S5Rng::new(seed);
        let mut draws: Vec<HighsInt> = (0..36)
            .map(|_| (replay.next_u64() % (n as u64)) as HighsInt)
            .collect();
        draws.sort_by(|a, b| b.cmp(a));
        let expected_kept: Vec<HighsInt> = draws[..12].to_vec();

        let mut rng = S5Rng::new(seed);
        let cand = infeasible_step(&mut ctx, &mut rng, 0, false, 0.0, &[], &[]);
        assert!(ctx.was_infeasible);
        let kept: Vec<HighsInt> = ctx.sampled.iter().map(|s| s.ci).collect();
        assert_eq!(kept, expected_kept);
        assert_eq!(cand.var_idx, expected_kept[0]);
        assert_eq!(cand.new_val, 1.0);
        assert_eq!(cand.score, (expected_kept[0] + 1) as f64);
        assert_eq!(rng.next_u64(), replay.next_u64());
    }

    /// Random problem with every column kind the breakthrough skip must
    /// handle: binary, integer, integer with fractional bounds, continuous
    /// with infinite bounds; costs of both signs, zero, and below kEpsZero.
    fn make_problem_random(rng: &mut S5Rng, ncol: HighsInt, nrow: HighsInt) -> S5Problem {
        let (n, m) = (ncol as usize, nrow as usize);
        let mut p = S5Problem {
            feastol: 1e-6,
            epsilon: 1e-6,
            ncol,
            nrow,
            ..S5Problem::default()
        };
        for _ in 0..n {
            let (lb, ub, ty) = match rng.next_u64() % 5 {
                0 => (0.0, 1.0, HighsVarType::kInteger),
                1 => (-3.0, 7.0, HighsVarType::kInteger),
                2 => (0.5, 4.5, HighsVarType::kInteger),
                3 => (-kHighsInf, kHighsInf, HighsVarType::kContinuous),
                _ => (-2.0, kHighsInf, HighsVarType::kContinuous),
            };
            let cost = match rng.next_u64() % 6 {
                0 => 0.0,
                1 => 1e-16,
                2 => s5_uniform_real(rng, -5.0, -0.5),
                3 => -1e-14,
                _ => s5_uniform_real(rng, 0.5, 5.0),
            };
            p.col_lower.push(lb);
            p.col_upper.push(ub);
            p.integrality.push(ty);
            p.col_cost.push(cost);
            p.binary.push(u8::from(lb == 0.0 && ub == 1.0));
        }
        // CSR first, then its transpose as the CSC.
        let mut csc: Vec<Vec<(HighsInt, f64)>> = vec![Vec::new(); n];
        p.ar_start.push(0);
        for i in 0..m {
            for j in 0..n {
                if rng.next_u64().is_multiple_of(8) {
                    let a = s5_uniform_real(rng, -3.0, 3.0);
                    p.ar_index.push(j as HighsInt);
                    p.ar_value.push(a);
                    csc[j].push((i as HighsInt, a));
                }
            }
            p.ar_start.push(p.ar_index.len() as HighsInt);
            p.row_lower.push(s5_uniform_real(rng, -4.0, 0.0));
            p.row_upper.push(s5_uniform_real(rng, 0.0, 4.0));
        }
        p.csc.col_start.push(0);
        for col in &csc {
            for &(i, a) in col {
                p.csc.col_row.push(i);
                p.csc.col_val.push(a);
            }
            p.csc.col_start.push(p.csc.col_row.len() as HighsInt);
        }
        p.nnz = p.ar_index.len();
        p.uplocks = vec![0; n];
        p.downlocks = vec![0; n];
        p.incumbent = vec![0.0; n];
        p
    }

    /// A random value for column `j`: a bound, a non-canonical interior
    /// point, or a rounded interior point.
    fn random_value(rng: &mut S5Rng, p: &S5Problem, j: usize) -> f64 {
        let (lb, ub) = (p.col_lower[j], p.col_upper[j]);
        let lo = if lb.is_finite() { lb } else { -6.0 };
        let hi = if ub.is_finite() { ub } else { 6.0 };
        match rng.next_u64() % 4 {
            0 => lb.max(-1e3),
            1 => ub.min(1e3),
            2 => s5_uniform_real(rng, lo, hi),
            _ => s5_uniform_real(rng, lo, hi).round().clamp(lo, hi),
        }
    }

    #[test]
    fn s5_lmsearch_breakthrough_live_set_matches_full_scan() {
        // rusolve (throughput): the
        // `bt_live` fast path of Phase 1b must append exactly the
        // candidates of the full costed-column loop, in the same order,
        // across bulk writes (rebuild_state), single moves (apply_move)
        // and best objectives above, at and below the current one.
        let mut rng = S5Rng::new(20261115);
        for trial in 0..40 {
            let problem = make_problem_random(&mut rng, 150, 12);
            let mut ctx = make_ctx(&problem);
            let costed: Vec<HighsInt> = (0..problem.ncol)
                .filter(|&j| problem.col_cost[j as usize].abs() >= kEpsZero)
                .collect();
            for j in 0..problem.ncol as usize {
                ctx.solution[j] = random_value(&mut rng, &problem, j);
            }
            ctx.rebuild_state();
            for step in 0..60 {
                if step % 20 == 19 {
                    for j in 0..problem.ncol as usize {
                        ctx.solution[j] = random_value(&mut rng, &problem, j);
                    }
                    ctx.rebuild_state();
                } else {
                    let j = (rng.next_u64() % problem.ncol as u64) as usize;
                    let v = random_value(&mut rng, &problem, j);
                    ctx.apply_move(j as HighsInt, v);
                }
                let cur = ctx.current_obj;
                let eps = ctx.epsilon;
                for best in [
                    cur,
                    cur - 0.5 * eps,
                    cur + 0.5 * eps,
                    cur - eps,
                    cur + 2.0 * eps,
                    cur - 1.0,
                    cur + 1.0,
                    cur - 1e3,
                    f64::INFINITY,
                    f64::NEG_INFINITY,
                ] {
                    let mut fast = Vec::new();
                    assert!(ctx.bt_live_valid);
                    append_breakthrough_candidates(&mut ctx, &mut fast, &costed, best);
                    ctx.bt_live_valid = false;
                    let mut full = Vec::new();
                    append_breakthrough_candidates(&mut ctx, &mut full, &costed, best);
                    ctx.bt_live_valid = true;
                    assert_eq!(fast, full, "trial {trial} step {step} best {best}");
                }
            }
        }
    }
}
