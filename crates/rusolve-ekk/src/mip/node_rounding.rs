//! rusolve `mip_node_rounding`: LP-free rounding of a tree node's LP solution, run inside dives.
//!
//! Not part of HiGHS, whose tree search rounds LP solutions only while it has no incumbent
//! (`randomizedRounding` at the start of a plunge). This is our own implementation of the three
//! LP-free rounding heuristics described in Achterberg, "Constraint Integer Programming" (PhD
//! thesis, TU Berlin, 2007), section 9.1, which SCIP (Apache-2.0) implements as
//! `heur_simplerounding.c`, `heur_rounding.c` and `heur_shifting.c`; see also Berthold, "Primal
//! Heuristics for Mixed Integer Programs" (diploma thesis, TU Berlin, 2006), chapter 3. SCIP runs
//! the first two at every node and shifting every tenth.
//!
//! The three stages, each run only if the previous one left work:
//! 1. Simple rounding: a fractional integer column with no down-locks (no row that a decrease
//!    can violate) is rounded down, one with no up-locks up; with neither lock, the direction
//!    that does not raise the objective. Rows stay feasible by construction.
//! 2. Rounding: the remaining fractional columns are rounded, in LP order, in the direction
//!    that violates the fewest rows given the current row activities; ties go to the smaller
//!    relative violation, then to the objective.
//! 3. Shifting repair: while a row is violated, shift one column of that row (an integer one
//!    by an integral step, a continuous one by the exact amount) toward feasibility, within
//!    the global bounds, choosing the shift with the smallest net change in the number of
//!    violated rows (then the smallest relative violation, then the objective). A column is
//!    never shifted back against its previous shift. The repair fails when the number of
//!    violated rows has not reached a new minimum for `MAX_NON_IMPROVING` shifts, or after a
//!    shift budget.
//!
//! With `mip_node_rounding_polish`, a found point is then improved by a 1-opt pass
//! ([`one_opt`]). Continuous columns keep their LP values unless stage 3 or the polish moves them. No LP is solved, so
//! one call costs a few passes over the fractional columns and the rows.

use crate::prelude::*;

/// Shifts without a new minimum of violated rows before the repair gives up.
const MAX_NON_IMPROVING: usize = 25;

/// What a call produced: the point, its objective and the stage (0 simple rounding,
/// 1 rounding, 2 shifting repair) that completed it.
pub(crate) struct NodeRoundingPoint {
    pub sol: Vec<f64>,
    pub objective: f64,
    pub stage: usize,
    /// whether the 1-opt polish lowered the objective
    pub polished: bool,
}

/// Row violation of activity `act` against `[lo, hi]`, beyond the tolerance.
#[inline]
fn excess(act: f64, lo: f64, hi: f64, tol: f64) -> f64 {
    if act > hi + tol {
        act - hi
    } else if act < lo - tol {
        lo - act
    } else {
        0.0
    }
}

/// Scale of row `r` for comparing violations across rows: the magnitude of its finite sides.
#[inline]
fn row_scale(model: &HighsLp, r: usize) -> f64 {
    let lo = model.row_lower_[r];
    let hi = model.row_upper_[r];
    let mut s: f64 = 1.0;
    if lo.is_finite() {
        s = s.max(lo.abs());
    }
    if hi.is_finite() {
        s = s.max(hi.abs());
    }
    s
}

/// Effect of moving column `j` by `delta`: (net change in the number of violated rows,
/// change of the scaled total violation).
#[inline]
fn move_effect(model: &HighsLp, act: &[f64], j: usize, delta: f64, tol: f64) -> (i64, f64) {
    let start = model.a_matrix_.start_[j] as usize;
    let end = model.a_matrix_.start_[j + 1] as usize;
    let mut net = 0i64;
    let mut rel = 0.0;
    for k in start..end {
        let r = model.a_matrix_.index_[k] as usize;
        let a = model.a_matrix_.value_[k];
        let lo = model.row_lower_[r];
        let hi = model.row_upper_[r];
        let before = excess(act[r], lo, hi, tol);
        let after = excess(act[r] + a * delta, lo, hi, tol);
        net += (after > 0.0) as i64 - (before > 0.0) as i64;
        rel += (after - before) / row_scale(model, r);
    }
    (net, rel)
}

#[inline]
fn apply_move(model: &HighsLp, act: &mut [f64], x: &mut [f64], j: usize, delta: f64) {
    x[j] += delta;
    let start = model.a_matrix_.start_[j] as usize;
    let end = model.a_matrix_.start_[j + 1] as usize;
    for k in start..end {
        act[model.a_matrix_.index_[k] as usize] += model.a_matrix_.value_[k] * delta;
    }
}

/// Lexicographic "a is better than b" on (net violated rows, scaled violation, objective).
#[inline]
fn better(a: (i64, f64, f64), b: (i64, f64, f64)) -> bool {
    a.0 < b.0 || (a.0 == b.0 && (a.1 < b.1 - 1e-12 || (a.1 <= b.1 + 1e-12 && a.2 < b.2)))
}

/// Rounds `lpsol` (the node LP solution; `frac` its fractional integer columns) to a point
/// that is integral and satisfies every model row, or returns `None`. `glb`/`gub` are the
/// global column bounds, which a shift must respect. With `polish`, a found point is improved
/// by [`one_opt`].
pub(crate) fn round_lp_point(
    ms: &HighsMipSolver,
    glb: &[f64],
    gub: &[f64],
    lpsol: &[f64],
    frac: &[(HighsInt, f64)],
    polish: bool,
) -> Option<NodeRoundingPoint> {
    let ncol = ms.num_col() as usize;
    if lpsol.len() < ncol || frac.is_empty() {
        return None;
    }
    let model: &HighsLp = mip_model!(ms);
    let tol = ms.mipdata_.feastol;
    let mut x: Vec<f64> = lpsol[..ncol].to_vec();

    // Stage 1: simple rounding along the locks.
    let mut hard: Vec<usize> = Vec::new();
    for &(col, val) in frac {
        let j = col as usize;
        let free_down = ms.mipdata_.downlocks[j] == 0;
        let free_up = ms.mipdata_.uplocks[j] == 0;
        x[j] = match (free_down, free_up) {
            (true, true) => {
                if model.col_cost_[j] >= 0.0 {
                    val.floor()
                } else {
                    val.ceil()
                }
            }
            (true, false) => val.floor(),
            (false, true) => val.ceil(),
            (false, false) => {
                hard.push(j);
                continue;
            }
        };
    }
    let mut stage = 0;
    // Row activities of the point so far (the hard columns still at their LP values).
    let nrow = ms.num_row() as usize;
    let mut act = vec![0.0; nrow];
    if !hard.is_empty() || polish {
        for (r, a) in act.iter_mut().enumerate() {
            let start = ms.mipdata_.a_rstart_[r] as usize;
            let end = ms.mipdata_.a_rstart_[r + 1] as usize;
            let mut s = 0.0;
            for k in start..end {
                s += x[ms.mipdata_.a_rindex_[k] as usize] * ms.mipdata_.a_rvalue_[k];
            }
            *a = s;
        }
    }

    if !hard.is_empty() {
        stage = 1;
        // Stage 2: activity-aware rounding of the locked columns.
        for &j in &hard {
            let v = x[j];
            let dn = v.floor() - v;
            let up = v.ceil() - v;
            let cost = model.col_cost_[j];
            let (nd, rd) = move_effect(model, &act, j, dn, tol);
            let (nu, ru) = move_effect(model, &act, j, up, tol);
            let delta = if better((nu, ru, cost * up), (nd, rd, cost * dn)) {
                up
            } else {
                dn
            };
            apply_move(model, &mut act, &mut x, j, delta);
        }

        // Stage 3: shifting repair of violated rows.
        let is_viol = |act: &[f64], r: usize| {
            excess(act[r], model.row_lower_[r], model.row_upper_[r], tol) > 0.0
        };
        let mut nviol = (0..nrow).filter(|&r| is_viol(&act, r)).count() as i64;
        if nviol > 0 {
            stage = 2;
            let max_shifts = 4 * hard.len() + 2 * nrow.min(1000) + 20;
            let mut last_dir: Vec<i8> = vec![0; ncol];
            let mut best_nviol = nviol;
            let mut non_improving = 0usize;
            let mut shifts = 0usize;
            let mut r = 0usize;
            while nviol > 0 {
                if shifts >= max_shifts || non_improving > MAX_NON_IMPROVING {
                    return None;
                }
                // next violated row, cyclically
                let mut found = false;
                for _ in 0..nrow {
                    if is_viol(&act, r) {
                        found = true;
                        break;
                    }
                    r = (r + 1) % nrow;
                }
                if !found {
                    break;
                }
                let lo = model.row_lower_[r];
                let hi = model.row_upper_[r];
                // sign of the activity change the row needs, and its size
                let (sign, need) = if act[r] > hi + tol {
                    (-1.0, act[r] - hi)
                } else {
                    (1.0, lo - act[r])
                };
                let start = ms.mipdata_.a_rstart_[r] as usize;
                let end = ms.mipdata_.a_rstart_[r + 1] as usize;
                let mut best: Option<((i64, f64, f64), usize, f64)> = None;
                for k in start..end {
                    let j = ms.mipdata_.a_rindex_[k] as usize;
                    let a = ms.mipdata_.a_rvalue_[k];
                    if a == 0.0 || glb[j] == gub[j] {
                        continue;
                    }
                    let dir: i8 = if sign * a > 0.0 { 1 } else { -1 };
                    if last_dir[j] == -dir {
                        continue;
                    }
                    let room = if dir > 0 {
                        gub[j] - x[j]
                    } else {
                        x[j] - glb[j]
                    };
                    if room <= tol {
                        continue;
                    }
                    let mut step = need / a.abs();
                    if ms.is_col_integer(j as HighsInt) {
                        step = (step - tol).ceil().max(1.0);
                        step = step.min((room + tol).floor());
                        if step < 1.0 {
                            continue;
                        }
                    } else {
                        step = step.min(room);
                    }
                    let delta = dir as f64 * step;
                    let (net, rel) = move_effect(model, &act, j, delta, tol);
                    if rel >= 0.0 && net >= 0 {
                        continue;
                    }
                    let key = (net, rel, model.col_cost_[j] * delta);
                    if best.as_ref().is_none_or(|(bk, _, _)| better(key, *bk)) {
                        best = Some((key, j, delta));
                    }
                }
                let ((net, _, _), j, delta) = best?;
                apply_move(model, &mut act, &mut x, j, delta);
                last_dir[j] = if delta > 0.0 { 1 } else { -1 };
                nviol += net;
                shifts += 1;
                if nviol < best_nviol {
                    best_nviol = nviol;
                    non_improving = 0;
                } else {
                    non_improving += 1;
                }
                if nviol <= 0 {
                    // recount: the incremental count is exact, but be safe against drift
                    nviol = (0..nrow).filter(|&i| is_viol(&act, i)).count() as i64;
                }
            }
        }
    }

    let mut objective = 0.0;
    for (j, &v) in x.iter().enumerate() {
        objective += model.col_cost_[j] * v;
    }
    let mut polished = false;
    if polish {
        let gain = one_opt(ms, model, glb, gub, &mut x, &mut act, tol);
        if gain > 0.0 {
            polished = true;
            objective = 0.0;
            for (j, &v) in x.iter().enumerate() {
                objective += model.col_cost_[j] * v;
            }
        }
    }
    Some(NodeRoundingPoint {
        sol: x,
        objective,
        stage,
        polished,
    })
}

/// Largest move of column `j` in direction `dir` (+1 up, -1 down) that keeps every row
/// feasible and the column within its global bounds; integral for an integer column.
fn max_shift(
    ms: &HighsMipSolver,
    model: &HighsLp,
    glb: &[f64],
    gub: &[f64],
    x: &[f64],
    act: &[f64],
    j: usize,
    dir: f64,
    tol: f64,
) -> f64 {
    let mut s = if dir > 0.0 {
        gub[j] - x[j]
    } else {
        x[j] - glb[j]
    };
    if s <= tol {
        return 0.0;
    }
    let start = model.a_matrix_.start_[j] as usize;
    let end = model.a_matrix_.start_[j + 1] as usize;
    for k in start..end {
        let r = model.a_matrix_.index_[k] as usize;
        let rate = model.a_matrix_.value_[k] * dir;
        if rate > 0.0 {
            let hi = model.row_upper_[r];
            if hi.is_finite() {
                s = s.min(((hi - act[r]).max(0.0) + tol) / rate);
            }
        } else if rate < 0.0 {
            let lo = model.row_lower_[r];
            if lo.is_finite() {
                s = s.min(((act[r] - lo).max(0.0) + tol) / -rate);
            }
        }
        if s <= tol {
            return 0.0;
        }
    }
    if ms.is_col_integer(j as HighsInt) {
        s = (s + 1e-9).floor();
    }
    if s.is_finite() { s } else { 0.0 }
}

/// 1-opt polish of a row-feasible point: move single columns toward a lower objective as far
/// as the rows and global bounds allow, best estimated gain first, in up to two passes. The
/// idea is SCIP's `heur_oneopt.c` (Apache-2.0), there applied to new incumbents with an LP
/// for the continuous columns; here every column moves on its own and no LP is solved.
/// Returns the objective decrease.
fn one_opt(
    ms: &HighsMipSolver,
    model: &HighsLp,
    glb: &[f64],
    gub: &[f64],
    x: &mut [f64],
    act: &mut [f64],
    tol: f64,
) -> f64 {
    let mut total = 0.0;
    for _pass in 0..2 {
        let mut cand: Vec<(f64, usize)> = Vec::new();
        for (j, &cost) in model.col_cost_.iter().enumerate() {
            if cost == 0.0 || glb[j] == gub[j] {
                continue;
            }
            let dir = if cost > 0.0 { -1.0 } else { 1.0 };
            let s = max_shift(ms, model, glb, gub, x, act, j, dir, tol);
            if s > 0.0 {
                cand.push((cost.abs() * s, j));
            }
        }
        if cand.is_empty() {
            break;
        }
        cand.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        let mut pass_gain = 0.0;
        for &(_, j) in &cand {
            let cost = model.col_cost_[j];
            let dir = if cost > 0.0 { -1.0 } else { 1.0 };
            let s = max_shift(ms, model, glb, gub, x, act, j, dir, tol);
            if s > 0.0 {
                apply_move(model, act, x, j, dir * s);
                pass_gain += cost.abs() * s;
            }
        }
        total += pass_gain;
        if pass_gain <= 0.0 {
            break;
        }
    }
    total
}

/// `RUSOLVE_NODE_ROUNDING_DRY=1`: count the points the heuristic finds without adding them,
/// so that the search is the one without the heuristic (a ceiling measurement).
pub(crate) fn node_rounding_dry() -> bool {
    static DRY: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *DRY.get_or_init(|| std::env::var("RUSOLVE_NODE_ROUNDING_DRY").is_ok_and(|v| v == "1"))
}

#[cfg(test)]
mod tests {
    use super::{better, excess};

    #[test]
    fn excess_respects_tolerance_and_infinite_bounds() {
        assert_eq!(excess(1.0, f64::NEG_INFINITY, 1.0, 1e-6), 0.0);
        assert_eq!(excess(1.5, f64::NEG_INFINITY, 1.0, 1e-6), 0.5);
        assert_eq!(excess(-2.0, -1.0, f64::INFINITY, 1e-6), 1.0);
        assert_eq!(excess(1.0 + 1e-7, 0.0, 1.0, 1e-6), 0.0);
    }

    #[test]
    fn better_is_lexicographic() {
        assert!(better((0, 5.0, 9.0), (1, 0.0, 0.0)));
        assert!(better((1, 0.1, 9.0), (1, 0.2, 0.0)));
        assert!(better((1, 0.1, 1.0), (1, 0.1, 2.0)));
        assert!(!better((1, 0.1, 2.0), (1, 0.1, 2.0)));
    }
}
