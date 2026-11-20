//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! mip-heuristics (36f7d9d): src/walksat.h + src/walksat.cpp (lines 1-407)
//!
//! The leaf-time (complete-assignment) WalkSAT repair of paper Fig. 4, plus
//! `walksat_select_move` (the move selector shared with `repair_search`) and
//! the greedy 1-opt improvement of Sect. 6.  Unlike `repair_walk` (which
//! shifts *domains* at an infeasible DFS node), these operate on a complete
//! `solution` vector and its `lhs_cache`, and the effort valve is *internal*
//! (`kWalkSatBudgetPerNnz`) rather than a caller argument (issue #156).

use crate::HighsInt;
use crate::lp_data::hconst::HighsVarType;
use crate::mip::fpr_local_mip::PropEngine;
use crate::mip::fpr_local_mip::WalkSatCandidate;
use crate::mip::fpr_local_mip::WalkSatMove;
use crate::mip::fpr_local_mip::WalkSatScratch;
use crate::mip::fpr_local_mip::WalkSatUndoEntry;
use crate::mip::fpr_local_mip::kWalkSatBudgetPerNnz;
use crate::mip::fpr_local_mip::rng::S5Rng;
use crate::mip::fpr_local_mip::rng::is_row_violated;
use crate::mip::fpr_local_mip::rng::row_violation;
use crate::mip::fpr_local_mip::rng::s5_uniform_int_i32;
use crate::mip::fpr_local_mip::rng::s5_uniform_real;
use crate::util::highs_utils::cpp_min;

/// mip-heuristics (36f7d9d): src/walksat.cpp::walksat_select_move (lines 39-172)
///
/// Cognitive complexity 31 (threshold 25).  Kept whole: WalkSAT move
/// selection: the noise / greedy / tabu case analysis.
///
/// Select a WalkSAT repair move for the given violated row.  Uses PropEngine
/// for read-only matrix data, bounds, integrality, feastol.
/// `solution`/`lhs_cache` are the current (possibly infeasible) assignment;
/// `col_lb`/`col_ub` are clamping bounds for candidate shifts.  Returns
/// `{-1, 0.0}` if no valid candidate found.  Increments `effort` by the
/// coefficient accesses consumed.
pub(crate) fn walksat_select_move(
    row: HighsInt,
    solution: &[f64],
    lhs_cache: &[f64],
    col_lb: &[f64],
    col_ub: &[f64],
    data: &PropEngine,
    noise: f64,
    rng: &mut S5Rng,
    effort: &mut usize,
    scratch: &mut WalkSatScratch,
) -> WalkSatMove {
    // The C++ `__restrict` locals are dropped (no semantic change).
    let ar_start = data.ar_start();
    let ar_index = data.ar_index();
    let ar_value = data.ar_value();
    let csc_start = data.csc_start();
    let csc_row = data.csc_row();
    let csc_val = data.csc_val();
    let row_lo = data.row_lo();
    let row_hi = data.row_hi();
    let integrality = data.integrality();
    let feastol = data.feastol();

    let row_begin = ar_start[row as usize];
    let row_end = ar_start[(row + 1) as usize];
    let row_len = row_end - row_begin;
    if row_len == 0 {
        return WalkSatMove::default();
    }

    // Loop-invariant row metadata.
    let row_lo_r = row_lo[row as usize];
    let row_hi_r = row_hi[row as usize];
    let ci_lhs = lhs_cache[row as usize];
    let ci_viol = row_violation(ci_lhs, row_lo_r, row_hi_r);
    let overshoot_above = ci_lhs > row_hi_r + feastol;
    let target_rhs = if overshoot_above { row_hi_r } else { row_lo_r };
    let ci_viol_threshold = ci_viol - feastol; // early-reject bound

    let cand = &mut scratch.cand;
    cand.clear();
    cand.reserve(row_len as usize);
    let mut best_damage = f64::INFINITY;

    *effort += row_len as usize;
    let mut k = row_begin;
    while k < row_end {
        let j = ar_index[k as usize];
        let a = ar_value[k as usize];
        if a.abs() < 1e-15 {
            k += 1;
            continue;
        }

        let old_val = solution[j as usize];
        let mut new_val = old_val + ((target_rhs - ci_lhs) / a);

        if integrality[j as usize] != HighsVarType::kContinuous {
            // Rounding direction depends on (overshoot_above) XOR (a < 0):
            //   overshoot_above && a > 0: floor to decrease lhs
            //   overshoot_above && a < 0: ceil  to decrease lhs
            //   undershoot    && a > 0: ceil  to increase lhs
            //   undershoot    && a < 0: floor to increase lhs
            let floor_it = overshoot_above == (a > 0.0);
            new_val = if floor_it {
                (new_val + feastol).floor()
            } else {
                (new_val - feastol).ceil()
            };
        }
        let lb = col_lb[j as usize];
        let ub = col_ub[j as usize];
        if new_val < lb {
            new_val = lb;
        } else if new_val > ub {
            new_val = ub;
        }

        let delta_change = new_val - old_val;
        if delta_change.abs() < 1e-15 {
            k += 1;
            continue;
        }

        let new_ci_lhs = ci_lhs + (a * delta_change);
        let new_ci_viol = row_violation(new_ci_lhs, row_lo_r, row_hi_r);
        if new_ci_viol >= ci_viol_threshold {
            k += 1;
            continue;
        }

        let col_begin = csc_start[j as usize];
        let col_end = csc_start[(j + 1) as usize];
        *effort += (col_end - col_begin) as usize;

        let mut damage = 0.0f64;
        let col_deg = col_end - col_begin;
        let mut p: HighsInt = 0;
        while p < col_deg {
            let i2 = csc_row[(col_begin + p) as usize];
            if i2 == row {
                p += 1;
                continue;
            }
            let coeff = csc_val[(col_begin + p) as usize];
            let old_lhs = lhs_cache[i2 as usize];
            let new_lhs = old_lhs + (coeff * delta_change);
            let lo_i = row_lo[i2 as usize];
            let hi_i = row_hi[i2 as usize];
            let dv = row_violation(new_lhs, lo_i, hi_i) - row_violation(old_lhs, lo_i, hi_i);
            if dv > 0.0 {
                damage += dv;
            }
            p += 1;
        }

        best_damage = cpp_min(damage, best_damage);
        cand.push(WalkSatCandidate {
            var: j,
            new_val,
            damage,
        });
        k += 1;
    }

    let ncand = cand.len();
    if ncand == 0 {
        return WalkSatMove::default();
    }

    // WalkSAT selection (paper Fig. 4, lines 17-21)
    if best_damage > feastol && s5_uniform_real(rng, 0.0, 1.0) < noise {
        let idx = s5_uniform_int_i32(rng, 0, ncand as HighsInt - 1);
        return WalkSatMove {
            var: cand[idx as usize].var,
            val: cand[idx as usize].new_val,
        };
    }
    let best_indices = &mut scratch.best_indices;
    best_indices.clear();
    best_indices.reserve(ncand);
    let damage_threshold = best_damage + feastol;
    let mut ci: usize = 0;
    while ci < ncand {
        if cand[ci].damage <= damage_threshold {
            best_indices.push(ci as HighsInt);
        }
        ci += 1;
    }
    let idx = best_indices[s5_uniform_int_i32(rng, 0, best_indices.len() as HighsInt - 1) as usize];
    WalkSatMove {
        var: cand[idx as usize].var,
        val: cand[idx as usize].new_val,
    }
}

/// mip-heuristics (36f7d9d): src/walksat.cpp::walksat_repair (lines 183-335)
///
/// Cognitive complexity 46 (threshold 25).  Kept whole: the WalkSAT repair
/// loop with its violated-set bookkeeping and undo trail, which must stay in
/// step with each other.
///
/// Run flat WalkSAT repair walk (paper Fig. 4) on `solution`/`lhs_cache`.
/// Returns true if feasible.  Modifies `solution`/`lhs_cache` in place.
/// `total_viol` is the current total violation (updated incrementally).
/// `scratch` holds reusable buffers -- cleared at entry.  `max_iterations` is
/// the paper's per-call step limit; the effort valve is *internal*
/// (`kWalkSatBudgetPerNnz`, issue #156).  `effort` is in-out: the caller's
/// running total is added to, so the valve gates the *delta* this call
/// charges.
pub(crate) fn walksat_repair(
    data: &PropEngine,
    solution: &mut Vec<f64>,
    lhs_cache: &mut Vec<f64>,
    col_lb: &[f64],
    col_ub: &[f64],
    max_iterations: HighsInt,
    noise: f64,
    track_best: bool,
    rng: &mut S5Rng,
    effort: &mut usize,
    scratch: &mut WalkSatScratch,
) -> bool {
    let nrow = data.nrow();
    let feastol = data.feastol();
    let row_lo = data.row_lo();
    let row_hi = data.row_hi();
    let csc_start = data.csc_start();
    let csc_row = data.csc_row();
    let csc_val = data.csc_val();

    // The internal effort valve (issue #156).  `effort` is the caller's
    // running total, so the budget is measured against the delta this call
    // charges -- which is what keeps the valve per call whatever a caller
    // hands in.
    let effort_at_entry = *effort;
    let effort_budget = kWalkSatBudgetPerNnz * (data.ar_start()[nrow as usize] as usize);

    // Reuse scratch's violated / violated_pos / sol_undo / lhs_undo to avoid
    // per-call heap allocations.  clear() retains capacity; violated_pos must
    // be sized to nrow and reset to -1 before use.
    scratch.violated.clear();
    scratch.sol_undo.clear();
    scratch.lhs_undo.clear();
    if scratch.violated_pos.len() < nrow as usize {
        scratch.violated_pos.clear();
        scratch.violated_pos.resize(nrow as usize, -1);
    } else {
        for slot in scratch.violated_pos[..nrow as usize].iter_mut() {
            *slot = -1;
        }
    }
    if scratch.violated.capacity() < nrow as usize {
        scratch.violated.reserve(nrow as usize);
    }

    // Initialize violated set and total violation
    let mut total_viol = 0.0f64;
    let mut i: HighsInt = 0;
    while i < nrow {
        let v = row_violation(
            lhs_cache[i as usize],
            row_lo[i as usize],
            row_hi[i as usize],
        );
        total_viol += v;
        if is_row_violated(
            lhs_cache[i as usize],
            row_lo[i as usize],
            row_hi[i as usize],
            feastol,
        ) {
            // add_violated(i)
            if scratch.violated_pos[i as usize] == -1 {
                scratch.violated_pos[i as usize] = scratch.violated.len() as HighsInt;
                scratch.violated.push(i);
            }
        }
        i += 1;
    }

    if scratch.violated.is_empty() {
        return true;
    }

    // Delta-based best-state tracking: instead of copying full solution/
    // lhs_cache vectors on each improvement (O(ncol + nrow)), we keep an undo
    // log of changed entries and record a mark at each improvement.  Restore
    // replays only the entries after the best mark -- O(changes_since_best).
    let mut best_viol = total_viol;
    let mut best_sol_mark: HighsInt = 0;
    let mut best_lhs_mark: HighsInt = 0;

    let mut step: HighsInt = 0;
    // stats.rs throughput counters (never read by the walk).
    let (mut ls_evals, mut ls_commits) = (0u64, 0u64);
    while step < max_iterations && !scratch.violated.is_empty() {
        if *effort - effort_at_entry >= effort_budget {
            break;
        }

        let pick = s5_uniform_int_i32(rng, 0, scratch.violated.len() as HighsInt - 1);
        let row = scratch.violated[pick as usize];

        let move_ = walksat_select_move(
            row, solution, lhs_cache, col_lb, col_ub, data, noise, rng, effort, scratch,
        );
        ls_evals += scratch.cand.len() as u64;
        if move_.var < 0 {
            step += 1;
            continue;
        }

        let old_val = solution[move_.var as usize];
        let delta = move_.val - old_val;
        if track_best {
            scratch.sol_undo.push(WalkSatUndoEntry {
                idx: move_.var,
                old_val,
            });
        }
        solution[move_.var as usize] = move_.val;
        ls_commits += 1;

        // Apply move: update lhs_cache, violated set, and total_viol
        // incrementally.
        *effort += (csc_start[(move_.var + 1) as usize] - csc_start[move_.var as usize]) as usize;
        let mut p = csc_start[move_.var as usize];
        while p < csc_start[(move_.var + 1) as usize] {
            let i2 = csc_row[p as usize];
            let old_v = row_violation(
                lhs_cache[i2 as usize],
                row_lo[i2 as usize],
                row_hi[i2 as usize],
            );
            if track_best {
                scratch.lhs_undo.push(WalkSatUndoEntry {
                    idx: i2,
                    old_val: lhs_cache[i2 as usize],
                });
            }
            lhs_cache[i2 as usize] += csc_val[p as usize] * delta;
            let new_v = row_violation(
                lhs_cache[i2 as usize],
                row_lo[i2 as usize],
                row_hi[i2 as usize],
            );
            total_viol += new_v - old_v;

            let was = scratch.violated_pos[i2 as usize] != -1;
            let now = is_row_violated(
                lhs_cache[i2 as usize],
                row_lo[i2 as usize],
                row_hi[i2 as usize],
                feastol,
            );
            if was && !now {
                // remove_violated(i2)
                let pos = scratch.violated_pos[i2 as usize];
                if pos != -1 {
                    let last = *scratch.violated.last().unwrap();
                    scratch.violated[pos as usize] = last;
                    scratch.violated_pos[last as usize] = pos;
                    scratch.violated.pop();
                    scratch.violated_pos[i2 as usize] = -1;
                }
            } else if !was && now {
                // add_violated(i2)
                if scratch.violated_pos[i2 as usize] == -1 {
                    scratch.violated_pos[i2 as usize] = scratch.violated.len() as HighsInt;
                    scratch.violated.push(i2);
                }
            }
            p += 1;
        }

        // Update best-state tracking (paper Fig. 4, lines 23-26)
        if track_best && total_viol < best_viol {
            best_viol = total_viol;
            best_sol_mark = scratch.sol_undo.len() as HighsInt;
            best_lhs_mark = scratch.lhs_undo.len() as HighsInt;
        }

        step += 1;
    }
    super::stats::FPR.add_moves(ls_evals, ls_commits);

    // Restore best state if tracking enabled (paper Fig. 4, line 27)
    if track_best && !scratch.violated.is_empty() {
        // Undo changes back to the best mark
        let mut u: HighsInt = scratch.sol_undo.len() as HighsInt - 1;
        while u >= best_sol_mark {
            solution[scratch.sol_undo[u as usize].idx as usize] =
                scratch.sol_undo[u as usize].old_val;
            u -= 1;
        }
        let mut u: HighsInt = scratch.lhs_undo.len() as HighsInt - 1;
        while u >= best_lhs_mark {
            lhs_cache[scratch.lhs_undo[u as usize].idx as usize] =
                scratch.lhs_undo[u as usize].old_val;
            u -= 1;
        }
        // Rebuild violated set from restored state
        for vi in scratch.violated.iter() {
            scratch.violated_pos[*vi as usize] = -1;
        }
        scratch.violated.clear();
        let mut i: HighsInt = 0;
        while i < nrow {
            if is_row_violated(
                lhs_cache[i as usize],
                row_lo[i as usize],
                row_hi[i as usize],
                feastol,
            ) {
                if scratch.violated_pos[i as usize] == -1 {
                    scratch.violated_pos[i as usize] = scratch.violated.len() as HighsInt;
                    scratch.violated.push(i);
                }
            }
            i += 1;
        }
    }

    scratch.violated.is_empty()
}

/// mip-heuristics (36f7d9d): src/walksat.cpp::greedy_1opt (lines 341-407)
///
/// Greedy 1-opt: shift each integer variable by +/-1 toward better objective
/// if the shift maintains feasibility.  Modifies `solution`/`lhs_cache` in
/// place.  `col_cost` is in minimisation form (`ProblemView::col_cost`).
pub(crate) fn greedy_1opt(
    data: &PropEngine,
    solution: &mut Vec<f64>,
    lhs_cache: &mut Vec<f64>,
    col_cost: &[f64],
    effort: &mut usize,
) {
    let ncol = data.ncol();
    let csc_start = data.csc_start();
    let csc_row = data.csc_row();
    let csc_val = data.csc_val();
    let col_lb = data.col_lb();
    let col_ub = data.col_ub();
    let row_lo = data.row_lo();
    let row_hi = data.row_hi();
    let integrality = data.integrality();
    let feastol = data.feastol();

    let mut j: HighsInt = 0;
    while j < ncol {
        if integrality[j as usize] == HighsVarType::kContinuous {
            j += 1;
            continue;
        }
        let cj = col_cost[j as usize];
        if cj.abs() < 1e-15 {
            j += 1;
            continue;
        }

        // Against the cost: `col_cost` is the view's minimisation form.
        let direction = if cj > 0.0 { -1.0 } else { 1.0 };
        let old_val = solution[j as usize];
        let mut new_val = old_val + direction;
        let lb = col_lb[j as usize];
        let ub = col_ub[j as usize];
        if new_val < lb {
            new_val = lb;
        } else if new_val > ub {
            new_val = ub;
        }
        let delta = new_val - old_val;
        if delta.abs() < 1e-15 {
            j += 1;
            continue;
        }

        let col_begin = csc_start[j as usize];
        let col_end = csc_start[(j + 1) as usize];
        let col_deg = col_end - col_begin;
        *effort += col_deg as usize;

        let mut shift_feasible = true;
        let mut p: HighsInt = 0;
        while p < col_deg {
            let row = csc_row[(col_begin + p) as usize];
            let new_lhs = lhs_cache[row as usize] + (csc_val[(col_begin + p) as usize] * delta);
            if is_row_violated(new_lhs, row_lo[row as usize], row_hi[row as usize], feastol) {
                shift_feasible = false;
                break;
            }
            p += 1;
        }

        if shift_feasible {
            let mut p: HighsInt = 0;
            while p < col_deg {
                lhs_cache[csc_row[(col_begin + p) as usize] as usize] +=
                    csc_val[(col_begin + p) as usize] * delta;
                p += 1;
            }
            solution[j as usize] = new_val;
        }

        j += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;
    use crate::lp_data::hconst::kHighsInf;
    use crate::mip::fpr_local_mip::PropEngine;
    use crate::mip::fpr_local_mip::WalkSatScratch;
    use crate::mip::fpr_local_mip::rng::S5Rng;
    use crate::mip::fpr_local_mip::rng::build_csc;

    // Two binary rows, each violated by the all-zero assignment and each
    // repaired by flipping one column to 1: `x0 >= 1` and `x1 >= 1`.  The
    // repair is deterministic (one candidate per row, zero damage), so the
    // walk converges regardless of the seed.
    struct TwoRowModel {
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

    impl TwoRowModel {
        fn new() -> Self {
            let ar_start = vec![0, 1, 2];
            let ar_index = vec![0, 1];
            let ar_value = vec![1.0, 1.0];
            let csc = build_csc(2, 2, &ar_start, &ar_index, &ar_value);
            TwoRowModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0, 0.0],
                col_ub: vec![1.0, 1.0],
                row_lo: vec![1.0, 1.0],
                row_hi: vec![kHighsInf, kHighsInf],
                integrality: vec![HighsVarType::kInteger, HighsVarType::kInteger],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                2,
                2,
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
                feastol,
            )
        }
    }

    #[test]
    fn s5_walk_walksat_repairs_two_row_violation() {
        let m = TwoRowModel::new();
        let data = m.make_engine(1e-6);
        let mut solution = vec![0.0, 0.0];
        let mut lhs_cache = vec![0.0, 0.0];
        let mut scratch = WalkSatScratch::default();
        let mut rng = S5Rng::new(42);
        let mut effort: usize = 0;

        let feasible = walksat_repair(
            &data,
            &mut solution,
            &mut lhs_cache,
            &m.col_lb,
            &m.col_ub,
            200,
            0.75,
            true,
            &mut rng,
            &mut effort,
            &mut scratch,
        );
        assert!(feasible);
        assert!(effort > 0);
        // Both rows satisfied: lhs >= 1 for each.
        assert!(!is_row_violated(
            lhs_cache[0],
            m.row_lo[0],
            m.row_hi[0],
            1e-6
        ));
        assert!(!is_row_violated(
            lhs_cache[1],
            m.row_lo[1],
            m.row_hi[1],
            1e-6
        ));
    }

    // One integer column y in [0, 10], one row `-inf <= y <= inf` (always
    // satisfied), cost +1.  The all-zero assignment is feasible, so greedy
    // 1-opt can shift y down by 1 (toward the cost) without violating the
    // row, and the objective cost*y decreases.
    struct GreedyModel {
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

    impl GreedyModel {
        fn new() -> Self {
            let ar_start = vec![0, 1];
            let ar_index = vec![0];
            let ar_value = vec![1.0];
            let csc = build_csc(1, 1, &ar_start, &ar_index, &ar_value);
            GreedyModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0],
                col_ub: vec![10.0],
                row_lo: vec![-kHighsInf],
                row_hi: vec![kHighsInf],
                integrality: vec![HighsVarType::kInteger],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                1,
                1,
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
                feastol,
            )
        }
    }

    #[test]
    fn s5_walk_greedy_1opt_improves_objective() {
        let m = GreedyModel::new();
        let data = m.make_engine(1e-6);
        let mut solution = vec![5.0];
        let mut lhs_cache = vec![5.0];
        let col_cost = vec![1.0];
        let mut effort: usize = 0;

        greedy_1opt(&data, &mut solution, &mut lhs_cache, &col_cost, &mut effort);

        // cost*x decreased (5 -> 4).
        assert!(col_cost[0] * solution[0] < col_cost[0] * 5.0);
        // No row became violated.
        assert!(!is_row_violated(
            lhs_cache[0],
            m.row_lo[0],
            m.row_hi[0],
            1e-6
        ));
    }

    // `y0 + y1 >= 5` with both columns binary: unsatisfiable, and the walk
    // cannot even stall quietly on it.  The first two steps push both columns
    // to their upper bound; from then on every candidate clips to the bound it
    // already holds, so `walksat_select_move` returns no move -- after
    // charging the row's length -- and the loop `continue`s.  With
    // `max_iterations` at 10^9 nothing but the internal valve can stop it.
    // One binary row `sum y_j >= ncol + 4`, which no assignment satisfies, so
    // every step finds a candidate, charges a row length and makes no
    // progress.  `ncol` is a parameter because the valve is
    // `kWalkSatBudgetPerNnz * nnz` and a single width cannot tell that apart
    // from a flat constant.
    struct UnrepairableRowModel {
        ncol: HighsInt,
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

    impl UnrepairableRowModel {
        fn new(n: HighsInt) -> Self {
            let ar_start = vec![0, n];
            let ar_value = vec![1.0; n as usize];
            let mut ar_index = Vec::new();
            for j in 0..n {
                ar_index.push(j);
            }
            let csc = build_csc(n, 1, &ar_start, &ar_index, &ar_value);
            UnrepairableRowModel {
                ncol: n,
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0; n as usize],
                col_ub: vec![1.0; n as usize],
                row_lo: vec![(n as f64) + 4.0],
                row_hi: vec![kHighsInf],
                integrality: vec![HighsVarType::kInteger; n as usize],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                self.ncol,
                1,
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
                feastol,
            )
        }
    }

    #[test]
    fn s5_walk_walksat_valve_bounds_walk() {
        // Port of test_repair_walk.cpp "the internal per-nnz valve bounds a
        // walk that cannot converge (#156)" (:1340).  Two widths, because the
        // valve scales with `nnz`: a single width would be satisfied by a flat
        // constant just as well.
        for width in [2, 8] {
            let m = UnrepairableRowModel::new(width);
            let data = m.make_engine(1e-6);
            let mut solution = vec![0.0; width as usize];
            let mut lhs_cache = vec![0.0];
            let mut scratch = WalkSatScratch::default();
            let mut rng = S5Rng::new(42);
            let mut effort: usize = 0;

            let feasible = walksat_repair(
                &data,
                &mut solution,
                &mut lhs_cache,
                &m.col_lb,
                &m.col_ub,
                1000000000,
                0.75,
                true,
                &mut rng,
                &mut effort,
                &mut scratch,
            );
            assert!(!feasible);

            // `kWalkSatBudgetPerNnz` is file-local to `src/walksat.cpp`, so the
            // value is spelled out here -- the same way `tests/test_prop_engine.cpp`
            // spells `100 * nnz` for `kPropagateBudgetPerNnz`.  The gate is
            // checked at the top of a step, and one step charges at most a row
            // length plus the columns it scores, so the overrun past the budget
            // is under one step -- well inside `2 * nnz`.  Without the valve
            // this call runs 10^9 steps.
            let nnz = m.ar_start[1] as usize;
            assert!(effort >= 100 * nnz);
            assert!(effort < (100 * nnz) + (2 * nnz));
        }
    }
}
