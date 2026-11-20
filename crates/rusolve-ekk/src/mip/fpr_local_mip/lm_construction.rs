//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! Construction phase for Local-MIP cold start
//! (local_mip_construction.h/cpp). Produces a (possibly infeasible)
//! starting assignment for the search loop.

use super::lm_core::round_tight_delta;
use super::rng::{S5Rng, clamp_round, is_integer, is_row_violated, row_violation, s5_std_sort};
use super::*;
use crate::prelude::*;

/// mip-heuristics (36f7d9d): local_mip_construction.h::construction_effort_cap (lines 59-61)
pub(crate) fn construction_effort_cap(max_effort: usize) -> usize {
    (max_effort as f64 * kConstructionEffortFraction) as usize
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::clamp_round_var (lines 54-56)
fn clamp_round_var(v: f64, lb: f64, ub: f64, integer: bool) -> f64 {
    clamp_round(v, lb, ub, integer)
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::zero_start_value (lines 60-70)
fn zero_start_value(lb: f64, ub: f64, integer: bool) -> f64 {
    let v: f64 = if lb > 0.0 {
        lb
    } else if ub < 0.0 {
        ub
    } else {
        0.0
    };
    clamp_round_var(v, lb, ub, integer)
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::compute_all_lhs (lines 74-86)
fn compute_all_lhs(
    nrow: HighsInt,
    ar_start: &[HighsInt],
    ar_index: &[HighsInt],
    ar_value: &[f64],
    solution: &[f64],
    lhs: &mut [f64],
    effort: &mut usize,
) {
    for i in 0..nrow {
        let mut l = 0.0;
        let mut k = ar_start[i as usize];
        while k < ar_start[(i + 1) as usize] {
            l += ar_value[k as usize] * solution[ar_index[k as usize] as usize];
            k += 1;
        }
        lhs[i as usize] = l;
    }
    *effort += ar_index.len();
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::CandidateEffect (lines 92-96)
struct CandidateEffect {
    delta_viol: f64,
    rows_newly_violated: HighsInt,
    rows_newly_satisfied: HighsInt,
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::evaluate_move (lines 98-123)
fn evaluate_move(
    j: HighsInt,
    old_val: f64,
    new_val: f64,
    csc: &CscMatrix,
    lhs: &[f64],
    row_lo: &[f64],
    row_hi: &[f64],
    feastol: f64,
    effort: &mut usize,
) -> CandidateEffect {
    let mut eff = CandidateEffect {
        delta_viol: 0.0,
        rows_newly_violated: 0,
        rows_newly_satisfied: 0,
    };
    let delta = new_val - old_val;
    if delta.abs() < 1e-15 {
        return eff;
    }
    *effort += (csc.col_start[(j + 1) as usize] - csc.col_start[j as usize]) as usize;
    let mut p = csc.col_start[j as usize];
    while p < csc.col_start[(j + 1) as usize] {
        let i = csc.col_row[p as usize];
        let old_lhs = lhs[i as usize];
        let new_lhs = old_lhs + (csc.col_val[p as usize] * delta);
        let old_v = row_violation(old_lhs, row_lo[i as usize], row_hi[i as usize]);
        let new_v = row_violation(new_lhs, row_lo[i as usize], row_hi[i as usize]);
        eff.delta_viol += new_v - old_v;
        let was_viol = is_row_violated(old_lhs, row_lo[i as usize], row_hi[i as usize], feastol);
        let now_viol = is_row_violated(new_lhs, row_lo[i as usize], row_hi[i as usize], feastol);
        if was_viol && !now_viol {
            eff.rows_newly_satisfied += 1;
        } else if !was_viol && now_viol {
            eff.rows_newly_violated += 1;
        }
        p += 1;
    }
    eff
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::apply_move_inplace (lines 126-139)
fn apply_move_inplace(
    j: HighsInt,
    new_val: f64,
    solution: &mut [f64],
    lhs: &mut [f64],
    csc: &CscMatrix,
    effort: &mut usize,
) {
    let old_val = solution[j as usize];
    let delta = new_val - old_val;
    if delta.abs() < 1e-15 {
        return;
    }
    *effort += (csc.col_start[(j + 1) as usize] - csc.col_start[j as usize]) as usize;
    solution[j as usize] = new_val;
    let mut p = csc.col_start[j as usize];
    while p < csc.col_start[(j + 1) as usize] {
        let i = csc.col_row[p as usize];
        lhs[i as usize] += csc.col_val[p as usize] * delta;
        p += 1;
    }
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::weighted_order (lines 183-208)
fn weighted_order(csc: &CscMatrix, ncol: HighsInt, rng: &mut S5Rng, out_order: &mut Vec<HighsInt>) {
    let mut tiebreak: Vec<u32> = vec![0u32; ncol as usize];
    out_order.resize(ncol as usize, 0);
    for j in 0..ncol {
        out_order[j as usize] = j;
        tiebreak[j as usize] = rng.next_u64() as u32;
    }
    // std::ranges::sort is libstdc++ introsort (R9); the comparator can tie
    // on (nnz, tiebreak) only when two variables share both keys, which is
    // astronomically unlikely but not impossible, so we use s5_std_sort.
    s5_std_sort(out_order, |a, b| {
        let nnz_a = csc.col_start[(*a + 1) as usize] - csc.col_start[*a as usize];
        let nnz_b = csc.col_start[(*b + 1) as usize] - csc.col_start[*b as usize];
        if nnz_a != nnz_b {
            return nnz_a > nnz_b;
        }
        tiebreak[*a as usize] < tiebreak[*b as usize]
    });
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::tight_delta_for_row (lines 212-244)
pub(crate) fn tight_delta_for_row(
    i: HighsInt,
    j: HighsInt,
    coeff: f64,
    lhs: &[f64],
    row_lo: &[f64],
    row_hi: &[f64],
    col_lb: &[f64],
    col_ub: &[f64],
    solution: &[f64],
    feastol: f64,
    integer: bool,
) -> f64 {
    if coeff.abs() < 1e-15 {
        return 0.0;
    }
    let l = lhs[i as usize];
    let gap: f64 = if l > row_hi[i as usize] + feastol {
        l - row_hi[i as usize]
    } else if l < row_lo[i as usize] - feastol {
        l - row_lo[i as usize]
    } else {
        return 0.0;
    };
    let mut delta = -gap / coeff;
    if integer {
        delta = round_tight_delta(delta, true);
    }
    let new_val = solution[j as usize] + delta;
    if new_val < col_lb[j as usize] {
        delta = col_lb[j as usize] - solution[j as usize];
    } else if new_val > col_ub[j as usize] {
        delta = col_ub[j as usize] - solution[j as usize];
    }
    delta
}

/// mip-heuristics (36f7d9d): local_mip_construction.cpp::construct_initial_solution (lines 252-426)
pub(crate) fn construct_initial_solution(
    problem: &S5Problem,
    rng: &mut S5Rng,
    max_effort: usize,
    out_solution: &mut Vec<f64>,
) -> usize {
    let ncol = problem.ncol;
    let nrow = problem.nrow;

    out_solution.clear();
    out_solution.resize(ncol as usize, 0.0);
    if ncol == 0 {
        return 0;
    }

    let col_lb = &problem.col_lower;
    let col_ub = &problem.col_upper;
    let row_lo = &problem.row_lower;
    let row_hi = &problem.row_upper;
    let integrality = &problem.integrality;
    let ar_start = &problem.ar_start;
    let ar_index = &problem.ar_index;
    let ar_value = &problem.ar_value;
    let csc = &problem.csc;
    let feastol = problem.feastol;

    // --- Phase A: zero-start (paper Alg 1 Line 1) -----------------------
    for j in 0..ncol {
        out_solution[j as usize] = zero_start_value(
            col_lb[j as usize],
            col_ub[j as usize],
            is_integer(integrality, j),
        );
    }

    if nrow == 0 || max_effort == 0 {
        return ncol as usize;
    }

    // Compute initial lhs[] vector; this is the anchor the greedy sweep
    // mutates incrementally.  Also acts as the effort "upfront" cost
    // (one full pass over nnz) against the construction cap.
    let mut lhs: Vec<f64> = vec![0.0; nrow as usize];
    let mut effort: usize = 0;
    compute_all_lhs(
        nrow,
        ar_start,
        ar_index,
        ar_value,
        out_solution,
        &mut lhs,
        &mut effort,
    );

    // --- Phase B: greedy variable sweep --------------------------------
    let mut order: Vec<HighsInt> = Vec::new();
    weighted_order(csc, ncol, rng, &mut order);

    // Small cap on candidates per variable — bounded so one variable
    // sweep stays O(col_nnz) coefficient accesses.
    const kMaxTightPerVar: HighsInt = 4;

    // Per-variable candidate buffer.
    let mut candidates: Vec<f64> = Vec::with_capacity(8);

    for &j in &order {
        if effort >= max_effort {
            break;
        }
        let integer = is_integer(integrality, j);
        let old_val = out_solution[j as usize];

        candidates.clear();
        // Always consider the current value as a "no-op" baseline.
        candidates.push(old_val);

        // Bounds (finite only; clamp_round_var handles integer rounding).
        if col_lb[j as usize] > -kHighsInf {
            candidates.push(clamp_round_var(
                col_lb[j as usize],
                col_lb[j as usize],
                col_ub[j as usize],
                integer,
            ));
        }
        if col_ub[j as usize] < kHighsInf {
            candidates.push(clamp_round_var(
                col_ub[j as usize],
                col_lb[j as usize],
                col_ub[j as usize],
                integer,
            ));
        }
        // Zero-clamp baseline.
        candidates.push(zero_start_value(
            col_lb[j as usize],
            col_ub[j as usize],
            integer,
        ));

        // Tight deltas from up to kMaxTightPerVar rows containing j that
        // are currently violated: each row gives the mtm operator's
        // "make this row tight" candidate (paper Def 4 / §3.2).
        let mut tight_added: HighsInt = 0;
        let mut p = csc.col_start[j as usize];
        while p < csc.col_start[(j + 1) as usize] && tight_added < kMaxTightPerVar {
            let i = csc.col_row[p as usize];
            if !is_row_violated(
                lhs[i as usize],
                row_lo[i as usize],
                row_hi[i as usize],
                feastol,
            ) {
                p += 1;
                continue;
            }
            let delta = tight_delta_for_row(
                i,
                j,
                csc.col_val[p as usize],
                &lhs,
                row_lo,
                row_hi,
                col_lb,
                col_ub,
                out_solution,
                feastol,
                integer,
            );
            if delta.abs() < 1e-15 {
                p += 1;
                continue;
            }
            let cand = clamp_round_var(
                old_val + delta,
                col_lb[j as usize],
                col_ub[j as usize],
                integer,
            );
            candidates.push(cand);
            tight_added += 1;
            p += 1;
        }

        // Evaluate each candidate.  Best = lowest (delta_viol), with a
        // feasibility-first rule: prefer candidates with strictly more
        // rows_newly_satisfied than rows_newly_violated over candidates
        // that merely reduce delta_viol.
        let mut best_new_val = old_val;
        let mut best_delta_viol = 0.0;
        let mut best_net_feas: HighsInt = 0;

        for &cand in &candidates {
            if (cand - old_val).abs() < 1e-15 {
                continue;
            }
            let eff = evaluate_move(
                j,
                old_val,
                cand,
                csc,
                &lhs,
                row_lo,
                row_hi,
                feastol,
                &mut effort,
            );
            let net_feas = eff.rows_newly_satisfied - eff.rows_newly_violated;

            // Feasibility-first: a strictly larger net_feas always wins,
            // regardless of delta_viol sign.  Within the same net_feas,
            // pick the smaller delta_viol.
            if net_feas > best_net_feas
                || (net_feas == best_net_feas && eff.delta_viol < best_delta_viol - 1e-12)
            {
                best_net_feas = net_feas;
                best_delta_viol = eff.delta_viol;
                best_new_val = cand;
            }
            if effort >= max_effort {
                break;
            }
        }

        if (best_new_val - old_val).abs() > 1e-15 {
            apply_move_inplace(j, best_new_val, out_solution, &mut lhs, csc, &mut effort);
        }
    }

    // Defensive final clamp — clamp_round_var + evaluate_move respect
    // bounds but rounding drift on pathological coefficients could push
    // a continuous variable ~feastol outside.  Cheap to enforce.
    for j in 0..ncol {
        let mut v = out_solution[j as usize];
        if v < col_lb[j as usize] {
            v = col_lb[j as usize];
        } else if v > col_ub[j as usize] {
            v = col_ub[j as usize];
        }
        if is_integer(integrality, j) {
            v = v.round();
            if v < col_lb[j as usize] {
                v = col_lb[j as usize].ceil();
            }
            if v > col_ub[j as usize] {
                v = col_ub[j as usize].floor();
            }
        }
        out_solution[j as usize] = v;
    }

    // Total effort = Phase B greedy sweep only.  Phase A's column-write
    // loop (~ncol cheap writes) is dropped from the normal-path total
    // for unit consistency.
    effort
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;

    /// Helper: build a 0/1 knapsack with 3 items, one <= row.
    /// Row 0: 2*x0 + 3*x1 + 4*x2 <= 5  (row_lo=-inf, row_hi=5)
    /// All integer, col bounds [0,1].
    fn make_knapsack_3() -> S5Problem {
        // CSR: row 0 has entries at col 0 (2), col 1 (3), col 2 (4)
        // CSC: col 0 -> row 0 (2), col 1 -> row 0 (3), col 2 -> row 0 (4)
        S5Problem {
            col_lower: vec![0.0, 0.0, 0.0],
            col_upper: vec![1.0, 1.0, 1.0],
            row_lower: vec![-kHighsInf],
            row_upper: vec![5.0],
            col_cost: vec![3.0, 4.0, 5.0],
            offset: 0.0,
            integrality: vec![
                HighsVarType::kInteger,
                HighsVarType::kInteger,
                HighsVarType::kInteger,
            ],
            ar_start: vec![0, 3],
            ar_index: vec![0, 1, 2],
            ar_value: vec![2.0, 3.0, 4.0],
            csc: CscMatrix {
                col_start: vec![0, 1, 2, 3],
                col_row: vec![0, 0, 0],
                col_val: vec![2.0, 3.0, 4.0],
            },
            uplocks: vec![0, 0, 0],
            downlocks: vec![0, 0, 0],
            feastol: 1e-6,
            epsilon: 1e-6,
            ncol: 3,
            nrow: 1,
            nnz: 3,
            incumbent: vec![0.0, 0.0, 0.0],
            binary: vec![0, 0, 0],
        }
    }

    /// Helper: build a 3-col problem with mixed bounds for zero-start test.
    fn make_mixed_bounds() -> S5Problem {
        // col 0: [2, 5] integer
        // col 1: [-7, -3] continuous
        // col 2: [-1, 4] integer
        // 1 row: x0 + x1 + x2 <= 10 (row_lo=-inf, row_hi=10)
        S5Problem {
            col_lower: vec![2.0, -7.0, -1.0],
            col_upper: vec![5.0, -3.0, 4.0],
            row_lower: vec![-kHighsInf],
            row_upper: vec![10.0],
            col_cost: vec![1.0, 1.0, 1.0],
            offset: 0.0,
            integrality: vec![
                HighsVarType::kInteger,
                HighsVarType::kContinuous,
                HighsVarType::kInteger,
            ],
            ar_start: vec![0, 3],
            ar_index: vec![0, 1, 2],
            ar_value: vec![1.0, 1.0, 1.0],
            csc: CscMatrix {
                col_start: vec![0, 1, 2, 3],
                col_row: vec![0, 0, 0],
                col_val: vec![1.0, 1.0, 1.0],
            },
            uplocks: vec![0, 0, 0],
            downlocks: vec![0, 0, 0],
            feastol: 1e-6,
            epsilon: 1e-6,
            ncol: 3,
            nrow: 1,
            nnz: 3,
            incumbent: vec![0.0, 0.0, 0.0],
            binary: vec![0, 0, 0],
        }
    }

    #[test]
    fn s5_lmcons_construct_knapsack_is_feasible_within_cap() {
        // The C++ checks the effort cap before each variable and after each
        // evaluate_move, then may still call apply_move_inplace, so the
        // returned effort can reach cap - 1 + 2 * (max col nnz).
        // Here max col nnz = 1, so the bound is cap + 1.
        let problem = make_knapsack_3();
        let mut rng = S5Rng::new(42);
        let cap = construction_effort_cap(1000);
        let mut sol: Vec<f64> = Vec::new();
        let effort = construct_initial_solution(&problem, &mut rng, cap, &mut sol);

        // Effort bound: cap + max col nnz (1 here).
        assert!(
            effort <= cap + 1,
            "effort {} exceeds cap {} + max col nnz 1",
            effort,
            cap
        );

        // Every value is 0.0 or 1.0.
        for &v in &sol {
            assert!(
                (v - 0.0).abs() < 1e-12 || (v - 1.0).abs() < 1e-12,
                "value {} is not 0 or 1",
                v
            );
        }

        // Row activity: 2*x0 + 3*x1 + 4*x2 <= 5 (feasibility).
        let activity = 2.0 * sol[0] + 3.0 * sol[1] + 4.0 * sol[2];
        assert!(
            activity <= 5.0 + 1e-6,
            "knapsack violated: activity {} > 5",
            activity
        );
    }

    #[test]
    fn s5_lmcons_zero_start_respects_bounds_lb_pos_ub_neg() {
        // max_effort 0 triggers the early return (Phase A only).
        let problem = make_mixed_bounds();
        let mut rng = S5Rng::new(42);
        let mut sol: Vec<f64> = Vec::new();
        let ret = construct_initial_solution(&problem, &mut rng, 0, &mut sol);

        // Return value is ncol (3) on the early-exit path.
        assert_eq!(ret, 3);

        // Phase A zero-start:
        // col 0: lb=2 > 0 -> v=2, clamp_round(2, 2, 5, integer) = 2.0
        // col 1: ub=-3 < 0 -> v=-3, clamp_round(-3, -7, -3, continuous) = -3.0
        // col 2: lb=-1 <= 0, ub=4 >= 0 -> v=0, clamp_round(0, -1, 4, integer) = 0.0
        assert_eq!(sol[0], 2.0);
        assert_eq!(sol[1], -3.0);
        assert_eq!(sol[2], 0.0);
    }

    #[test]
    fn s5_lmcons_clamp_round_var_edges() {
        // zero_start_value(0.4, 3.0, true):
        //   lb=0.4 > 0 -> v = 0.4
        //   clamp_round(0.4, 0.4, 3.0, true):
        //     round(0.4) = 0.0 (C std::round rounds half away from zero, 0.4 -> 0)
        //     cpp_min(3.0, 0.0) = 0.0
        //     cpp_max(0.4, 0.0) = 0.4
        //   result: 0.4
        let v = zero_start_value(0.4, 3.0, true);
        assert_eq!(v, 0.4);

        // zero_start_value(-2.5, -0.5, true):
        //   ub=-0.5 < 0 -> v = -0.5
        //   clamp_round(-0.5, -2.5, -0.5, true):
        //     round(-0.5) = -1.0 (C std::round rounds half away from zero)
        //     cpp_min(-0.5, -1.0) = -1.0
        //     cpp_max(-2.5, -1.0) = -1.0
        //   result: -1.0
        let v = zero_start_value(-2.5, -0.5, true);
        assert_eq!(v, -1.0);

        // zero_start_value(0.0, 0.0, false):
        //   lb=0 not > 0, ub=0 not < 0 -> v = 0.0
        //   clamp_round(0.0, 0.0, 0.0, false) = 0.0
        let v = zero_start_value(0.0, 0.0, false);
        assert_eq!(v, 0.0);
    }
}
