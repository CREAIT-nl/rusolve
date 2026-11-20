//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! value selection strategies for FPR (fpr_val_select.cpp/h).

use super::{CscMatrix, ValStrategy};
use crate::HighsInt;
use crate::lp_data::hconst::kHighsInf;
use crate::mip::fpr_local_mip::rng::{S5Rng, s5_uniform_real};
use crate::util::highs_utils::{cpp_max, cpp_min};

/// mip-heuristics (36f7d9d): fpr_val_select.cpp::val_up (lines 15-17)
fn val_up(ub: f64) -> f64 {
    ub
}

/// mip-heuristics (36f7d9d): fpr_val_select.cpp::val_random (lines 19-25)
fn val_random(lb: f64, ub: f64, is_int: bool, rng: &mut S5Rng) -> f64 {
    let mut v = s5_uniform_real(rng, lb, ub);
    if is_int {
        v = v.round();
    }
    cpp_max(lb, cpp_min(ub, v))
}

/// mip-heuristics (36f7d9d): fpr_val_select.cpp::val_goodobj (lines 29-34)
fn val_goodobj(lb: f64, ub: f64, cost: f64) -> f64 {
    if cost.abs() < 1e-15 {
        return lb; // paper Section 4.2: "always pick the lower bound"
    }
    if cost > 0.0 { lb } else { ub }
}

/// mip-heuristics (36f7d9d): fpr_val_select.cpp::val_badobj (lines 36-42)
fn val_badobj(lb: f64, ub: f64, cost: f64) -> f64 {
    if cost.abs() < 1e-15 {
        return lb; // paper Section 4.2: "always pick the lower bound"
    }
    // Opposite of goodobj
    if cost > 0.0 { ub } else { lb }
}

/// mip-heuristics (36f7d9d): fpr_val_select.cpp::val_lp_based (lines 46-61)
fn val_lp_based(lb: f64, ub: f64, is_int: bool, lp_val: f64, rng: &mut S5Rng) -> f64 {
    let clamped = cpp_max(lb, cpp_min(ub, lp_val));
    if !is_int {
        return clamped;
    }

    let f = clamped - clamped.floor();
    let v;
    if f < 1e-10 {
        v = clamped.round();
    } else {
        let u = s5_uniform_real(rng, 0.0, 1.0);
        if u < f {
            v = clamped.ceil();
        } else {
            v = clamped.floor();
        }
    }
    cpp_max(lb, cpp_min(ub, v))
}

/// mip-heuristics (36f7d9d): fpr_val_select.cpp::val_loosedyn (lines 72-134)
#[allow(clippy::too_many_arguments)]
fn val_loosedyn(
    j: HighsInt,
    lb: f64,
    ub: f64,
    _is_int: bool,
    cost: f64,
    row_lo: &[f64],
    row_hi: &[f64],
    min_act: &[f64],
    max_act: &[f64],
    csc: &CscMatrix,
) -> f64 {
    // Count dynamic up-locks and down-locks for variable j
    let mut up_locks: HighsInt = 0;
    let mut down_locks: HighsInt = 0;
    let total_rows = csc.col_start[(j + 1) as usize] - csc.col_start[j as usize];

    for r in 0..total_rows {
        let p = csc.col_start[j as usize] + r;
        let i = csc.col_row[p as usize];
        let a = csc.col_val[p as usize];

        // Check if constraint is infeasible or redundant
        let has_upper = row_hi[i as usize] < kHighsInf;
        let has_lower = row_lo[i as usize] > -kHighsInf;

        if has_upper && min_act[i as usize] > row_hi[i as usize] {
            // Upper bound infeasible — skip
            continue;
        }
        if has_lower && max_act[i as usize] < row_lo[i as usize] {
            // Lower bound infeasible — skip
            continue;
        }
        // Paper Section 4.2: "constraints that have already become redundant
        // do not contribute to the lock counters."
        if (!has_lower || min_act[i as usize] >= row_lo[i as usize])
            && (!has_upper || max_act[i as usize] <= row_hi[i as usize])
        {
            continue;
        }
        // Active constraint: count locks
        if has_upper && a > 0.0 {
            up_locks += 1;
        }
        if has_lower && a < 0.0 {
            up_locks += 1;
        }
        if has_upper && a < 0.0 {
            down_locks += 1;
        }
        if has_lower && a > 0.0 {
            down_locks += 1;
        }

        // Early exit: if one direction already wins by more than remaining rows
        let remaining = total_rows - r - 1;
        if up_locks > down_locks + remaining {
            return lb;
        }
        if down_locks > up_locks + remaining {
            return ub;
        }
    }

    // Pick direction with fewer locks
    if up_locks < down_locks {
        return ub;
    }
    if down_locks < up_locks {
        return lb;
    }
    // Tie: fall back to objective direction
    val_goodobj(lb, ub, cost)
}

/// mip-heuristics (36f7d9d): fpr_val_select.cpp::choose_value (lines 138-182)
#[allow(clippy::too_many_arguments)]
pub(crate) fn choose_value(
    j: HighsInt,
    lb: f64,
    ub: f64,
    is_int: bool,
    cost: f64,
    strategy: ValStrategy,
    rng: &mut S5Rng,
    lp_ref: Option<&[f64]>,
    row_lo: Option<&[f64]>,
    row_hi: Option<&[f64]>,
    min_act: Option<&[f64]>,
    max_act: Option<&[f64]>,
    csc: Option<&CscMatrix>,
) -> f64 {
    let mut v;
    match strategy {
        ValStrategy::kUp => {
            v = val_up(ub);
        }
        ValStrategy::kRandom => {
            v = val_random(lb, ub, is_int, rng);
        }
        ValStrategy::kGoodobj => {
            v = val_goodobj(lb, ub, cost);
        }
        ValStrategy::kBadobj => {
            v = val_badobj(lb, ub, cost);
        }
        ValStrategy::kLoosedyn => {
            if let (Some(row_lo), Some(row_hi), Some(min_act), Some(max_act), Some(csc)) =
                (row_lo, row_hi, min_act, max_act, csc)
            {
                v = val_loosedyn(
                    j, lb, ub, is_int, cost, row_lo, row_hi, min_act, max_act, csc,
                );
            } else {
                v = val_goodobj(lb, ub, cost);
            }
        }
        ValStrategy::kZerocore | ValStrategy::kZerolp | ValStrategy::kCore | ValStrategy::kLp => {
            if let Some(lp_ref) = lp_ref {
                v = val_lp_based(lb, ub, is_int, lp_ref[j as usize], rng);
            } else {
                v = val_goodobj(lb, ub, cost);
            }
        }
    }
    if is_int {
        v = v.round();
    }
    cpp_max(lb, cpp_min(ub, v))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::fpr_local_mip::rng::S5Rng;

    #[test]
    fn s5_varsel_val_up() {
        let mut rng = S5Rng::new(42);
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            1.0,
            ValStrategy::kUp,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 10.0);
    }

    #[test]
    fn s5_varsel_val_goodobj_pos_cost() {
        let mut rng = S5Rng::new(42);
        // cost > 0 -> lb
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            1.0,
            ValStrategy::kGoodobj,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 0.0);
    }

    #[test]
    fn s5_varsel_val_goodobj_neg_cost() {
        let mut rng = S5Rng::new(42);
        // cost < 0 -> ub
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            -1.0,
            ValStrategy::kGoodobj,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 10.0);
    }

    #[test]
    fn s5_varsel_val_goodobj_zero_cost() {
        let mut rng = S5Rng::new(42);
        // |cost| < 1e-15 -> lb
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            0.0,
            ValStrategy::kGoodobj,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 0.0);
    }

    #[test]
    fn s5_varsel_val_badobj_pos_cost() {
        let mut rng = S5Rng::new(42);
        // cost > 0 -> ub (opposite of goodobj)
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            1.0,
            ValStrategy::kBadobj,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 10.0);
    }

    #[test]
    fn s5_varsel_val_badobj_neg_cost() {
        let mut rng = S5Rng::new(42);
        // cost < 0 -> lb (opposite of goodobj)
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            -1.0,
            ValStrategy::kBadobj,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 0.0);
    }

    #[test]
    fn s5_varsel_val_badobj_zero_cost() {
        let mut rng = S5Rng::new(42);
        // |cost| < 1e-15 -> lb
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            0.0,
            ValStrategy::kBadobj,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 0.0);
    }

    #[test]
    fn s5_varsel_val_random_int() {
        let mut rng = S5Rng::new(42);
        let v = choose_value(
            0,
            0.0,
            10.0,
            true,
            1.0,
            ValStrategy::kRandom,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert!(v >= 0.0 && v <= 10.0);
        assert_eq!(v, v.round());
    }

    #[test]
    fn s5_varsel_val_lp_based_int() {
        let mut rng = S5Rng::new(42);
        let lp_ref = vec![2.5];
        let v = choose_value(
            0,
            0.0,
            10.0,
            true,
            1.0,
            ValStrategy::kLp,
            &mut rng,
            Some(lp_ref.as_slice()),
            None,
            None,
            None,
            None,
            None,
        );
        assert!(v >= 0.0 && v <= 10.0);
        assert_eq!(v, v.round());
    }

    #[test]
    fn s5_varsel_val_loosedyn_fallback() {
        let mut rng = S5Rng::new(42);
        // csc is None, so fallback to goodobj
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            1.0,
            ValStrategy::kLoosedyn,
            &mut rng,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 0.0); // goodobj with cost > 0 is lb
    }

    #[test]
    fn s5_varsel_val_clamp() {
        let mut rng = S5Rng::new(42);
        // lp_ref is out of bounds
        let lp_ref = vec![20.0];
        let v = choose_value(
            0,
            0.0,
            10.0,
            false,
            1.0,
            ValStrategy::kLp,
            &mut rng,
            Some(lp_ref.as_slice()),
            None,
            None,
            None,
            None,
            None,
        );
        assert_eq!(v, 10.0);
    }
}
