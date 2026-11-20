//! Activity-based bound propagation.
//!
//! For a row `lo ≤ aᵀx ≤ hi` the minimum activity over the current box
//! bounds each variable: `a_j x_j ≤ hi − (minact − mincontrib_j)`. Rows with
//! two or more infinite contributions cannot propagate; with exactly one, only
//! that variable can be tightened. With an incumbent, the objective row
//! `costᵀx ≤ cutoff` is propagated the same way (objective propagation).
//!
//! Provenance: bound strengthening as in Savelsbergh, "Preprocessing and
//! probing techniques for mixed integer programming problems", ORSA J.
//! Computing 6(4):445–454 (1994); acceptance rules and constants from HiGHS
//! `HighsDomain::adjustedUb/adjustedLb/propagateRowUpper/propagateRowLower`
//! (reference/HiGHS/highs/mip/HighsDomain.cpp:1349–1498, MIT) — integer
//! bounds rounded with feastol slack and accepted only if they move by more
//! than `1000·feastol·|bound|`; continuous bounds only if they cut ≥ 30% of
//! the domain; candidate bounds with `|b|·1e-14 > feastol` ignored.
// TODO(perf): activities are recomputed per visited row; HiGHS maintains
// them incrementally on every bound change (HighsDomain.cpp:1530–1870).

use crate::problem::MipProblem;
use rusolve_core::{CscMatrix, CsrMatrix, Triplet};

/// `|b|·TINY > feastol` ⇒ the candidate bound is numerically meaningless.
const TINY: f64 = 1e-14;
/// Continuous bounds snap to the opposite bound within this distance.
const EPSILON: f64 = 1e-9;
/// Continuous tightenings must remove at least this share of the domain.
const MIN_RELATIVE_IMPROVEMENT: f64 = 0.3;

pub struct Propagator {
    /// Constraint rows, then the objective row at index `num_row`.
    rows: CsrMatrix,
    /// Column-wise view of `rows`, to find the rows a bound change touches.
    cols: CscMatrix,
    row_lower: Vec<f64>,
    row_upper: Vec<f64>,
    integer: Vec<bool>,
    feastol: f64,
}

impl Propagator {
    pub fn new(p: &MipProblem, feastol: f64) -> Self {
        let mut trip: Vec<Triplet> = Vec::with_capacity(p.matrix.nnz() + p.num_col);
        let (ptr, idx, val) = (p.matrix.col_ptr(), p.matrix.row_idx(), p.matrix.values());
        for j in 0..p.num_col {
            for k in ptr[j]..ptr[j + 1] {
                if val[k] != 0.0 {
                    trip.push((idx[k] as usize, j, val[k]));
                }
            }
            if p.cost[j] != 0.0 {
                trip.push((p.num_row, j, p.cost[j]));
            }
        }
        let cols = CscMatrix::from_triplets(p.num_row + 1, p.num_col, &trip);
        let mut row_lower = p.row_lower.clone();
        let mut row_upper = p.row_upper.clone();
        row_lower.push(f64::NEG_INFINITY);
        row_upper.push(f64::INFINITY);
        Self {
            rows: cols.to_csr(),
            cols,
            row_lower,
            row_upper,
            integer: p.integer.clone(),
            feastol,
        }
    }

    /// Tightens `lower`/`upper` to a propagation fixpoint (or a work limit).
    /// `changed_col = None` seeds every row; `Some(j)` only the rows of `j`.
    /// `cutoff` bounds the objective row. Returns `false` iff infeasible.
    pub fn propagate(
        &mut self,
        lower: &mut [f64],
        upper: &mut [f64],
        changed_col: Option<usize>,
        cutoff: f64,
    ) -> bool {
        match changed_col {
            None => self.propagate_rows(lower, upper, None, cutoff),
            Some(j) => self.propagate_rows(lower, upper, Some(&[j]), cutoff),
        }
    }

    /// [`Self::propagate`] seeded with the rows of every column in
    /// `changed_cols` (bounds moved outside propagation, e.g. by reduced-cost
    /// fixing).
    pub fn propagate_cols(
        &mut self,
        lower: &mut [f64],
        upper: &mut [f64],
        changed_cols: &[usize],
        cutoff: f64,
    ) -> bool {
        self.propagate_rows(lower, upper, Some(changed_cols), cutoff)
    }

    fn propagate_rows(
        &mut self,
        lower: &mut [f64],
        upper: &mut [f64],
        changed_cols: Option<&[usize]>,
        cutoff: f64,
    ) -> bool {
        let num_rows = self.row_upper.len();
        self.row_upper[num_rows - 1] = cutoff;
        let mut queued = vec![false; num_rows];
        let mut queue: Vec<usize> = match changed_cols {
            None => (0..num_rows).collect(),
            Some(cols) => {
                let mut q = Vec::new();
                for &j in cols {
                    for r in self.rows_of(j) {
                        if !queued[r] {
                            queued[r] = true;
                            q.push(r);
                        }
                    }
                }
                q
            }
        };
        for &r in &queue {
            queued[r] = true;
        }
        if cutoff.is_finite() && !queued[num_rows - 1] {
            queued[num_rows - 1] = true;
            queue.push(num_rows - 1);
        }
        let mut changes: Vec<usize> = Vec::new();
        let mut budget = 20 * num_rows + 1000;
        while let Some(r) = queue.pop() {
            queued[r] = false;
            if budget == 0 {
                break;
            }
            budget -= 1;
            changes.clear();
            if !self.propagate_row(r, lower, upper, &mut changes) {
                return false;
            }
            for &j in &changes {
                for rr in self.rows_of(j) {
                    if !queued[rr] {
                        queued[rr] = true;
                        queue.push(rr);
                    }
                }
            }
        }
        true
    }

    fn rows_of(&self, j: usize) -> impl Iterator<Item = usize> + '_ {
        let ptr = self.cols.col_ptr();
        self.cols.row_idx()[ptr[j]..ptr[j + 1]]
            .iter()
            .map(|&r| r as usize)
    }

    /// Propagates one row; pushes changed columns. `false` iff infeasible.
    fn propagate_row(
        &self,
        r: usize,
        lower: &mut [f64],
        upper: &mut [f64],
        changes: &mut Vec<usize>,
    ) -> bool {
        let ptr = self.rows.row_ptr();
        let (cols, vals) = (
            &self.rows.col_idx()[ptr[r]..ptr[r + 1]],
            &self.rows.values()[ptr[r]..ptr[r + 1]],
        );
        let (lo, hi) = (self.row_lower[r], self.row_upper[r]);
        let (mut min_act, mut min_inf) = (Sum::default(), 0usize);
        let (mut max_act, mut max_inf) = (Sum::default(), 0usize);
        for (&j, &a) in cols.iter().zip(vals) {
            let (cmin, cmax) = contributions(a, lower[j as usize], upper[j as usize]);
            if cmin == f64::NEG_INFINITY {
                min_inf += 1;
            } else {
                min_act.add(cmin);
            }
            if cmax == f64::INFINITY {
                max_inf += 1;
            } else {
                max_act.add(cmax);
            }
        }
        let (min_act, max_act) = (min_act.value(), max_act.value());
        if (min_inf == 0 && min_act > hi + self.feastol)
            || (max_inf == 0 && max_act < lo - self.feastol)
        {
            return false;
        }
        // Upper side: a_j x_j ≤ hi − (residual minimum activity).
        if hi.is_finite() && min_inf <= 1 && (max_inf != 0 || max_act > hi + self.feastol) {
            for (&j, &a) in cols.iter().zip(vals) {
                let j = j as usize;
                let (cmin, _) = contributions(a, lower[j], upper[j]);
                let residual = match (min_inf, cmin == f64::NEG_INFINITY) {
                    (1, true) => min_act,
                    (1, false) => continue,
                    _ => min_act - cmin,
                };
                let bound = (hi - residual) / a;
                if !self.tighten(j, bound, a > 0.0, lower, upper, changes) {
                    return false;
                }
            }
        }
        // Lower side: a_j x_j ≥ lo − (residual maximum activity). The
        // activities above may be stale after upper-side changes; they are
        // still valid (looser) bounds, and the row is requeued on change.
        if lo.is_finite() && max_inf <= 1 && (min_inf != 0 || min_act < lo - self.feastol) {
            for (&j, &a) in cols.iter().zip(vals) {
                let j = j as usize;
                let (_, cmax) = contributions(a, lower[j], upper[j]);
                let residual = match (max_inf, cmax == f64::INFINITY) {
                    (1, true) => max_act,
                    (1, false) => continue,
                    _ => max_act - cmax,
                };
                let bound = (lo - residual) / a;
                if !self.tighten(j, bound, a < 0.0, lower, upper, changes) {
                    return false;
                }
            }
        }
        true
    }

    /// Applies candidate bound `b` as an upper (`is_upper`) or lower bound on
    /// column `j` under the acceptance rules. `false` iff the domain empties.
    fn tighten(
        &self,
        j: usize,
        b: f64,
        is_upper: bool,
        lower: &mut [f64],
        upper: &mut [f64],
        changes: &mut Vec<usize>,
    ) -> bool {
        if !b.is_finite() || (b * TINY).abs() > self.feastol {
            return true;
        }
        let ft = self.feastol;
        let (lo, hi) = (lower[j], upper[j]);
        if is_upper {
            let accept_bound = if self.integer[j] {
                let b = (b + ft).floor();
                (b < hi && hi - b > 1000.0 * ft * b.abs()).then_some(b)
            } else {
                let b = if (b - lo).abs() <= EPSILON { lo } else { b };
                let accept = if hi == f64::INFINITY {
                    true
                } else if b + 1000.0 * ft < hi {
                    let span = if lo.is_finite() {
                        hi - lo
                    } else {
                        hi.abs().max(b.abs())
                    };
                    (hi - b) / span >= MIN_RELATIVE_IMPROVEMENT
                } else {
                    false
                };
                accept.then_some(b)
            };
            let Some(b) = accept_bound else { return true };
            if b < lo - ft {
                return false;
            }
            upper[j] = b.max(lo);
        } else {
            let accept_bound = if self.integer[j] {
                let b = (b - ft).ceil();
                (b > lo && b - lo > 1000.0 * ft * b.abs()).then_some(b)
            } else {
                let b = if (hi - b).abs() <= EPSILON { hi } else { b };
                let accept = if lo == f64::NEG_INFINITY {
                    true
                } else if b - 1000.0 * ft > lo {
                    let span = if hi.is_finite() {
                        hi - lo
                    } else {
                        lo.abs().max(b.abs())
                    };
                    (b - lo) / span >= MIN_RELATIVE_IMPROVEMENT
                } else {
                    false
                };
                accept.then_some(b)
            };
            let Some(b) = accept_bound else { return true };
            if b > hi + ft {
                return false;
            }
            lower[j] = b.min(hi);
        }
        changes.push(j);
        true
    }
}

/// (min, max) of `a·x` over `x ∈ [lo, hi]`, with ±∞ for unbounded sides.
/// Zero coefficients never reach here (filtered when the rows are built).
fn contributions(a: f64, lo: f64, hi: f64) -> (f64, f64) {
    if a > 0.0 {
        (a * lo, a * hi)
    } else {
        (a * hi, a * lo)
    }
}

/// Neumaier-compensated sum: activities of big-M rows cancel badly in plain
/// f64, and a false "min activity > rhs" prunes feasible nodes. HiGHS uses
/// double-double (`HighsCDouble`) for the same reason.
#[derive(Default)]
struct Sum {
    s: f64,
    c: f64,
}

impl Sum {
    fn add(&mut self, v: f64) {
        let t = self.s + v;
        self.c += if self.s.abs() >= v.abs() {
            (self.s - t) + v
        } else {
            (v - t) + self.s
        };
        self.s = t;
    }

    fn value(&self) -> f64 {
        self.s + self.c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn problem(
        rows: &[(&[f64], f64, f64)],
        lower: &[f64],
        upper: &[f64],
        integer: bool,
    ) -> MipProblem {
        let n = lower.len();
        let mut trip = Vec::new();
        for (i, (coefs, _, _)) in rows.iter().enumerate() {
            for (j, &a) in coefs.iter().enumerate() {
                if a != 0.0 {
                    trip.push((i, j, a));
                }
            }
        }
        MipProblem {
            num_col: n,
            num_row: rows.len(),
            matrix: CscMatrix::from_triplets(rows.len(), n, &trip),
            row_lower: rows.iter().map(|r| r.1).collect(),
            row_upper: rows.iter().map(|r| r.2).collect(),
            cost: vec![1.0; n],
            maximize: false,
            col_lower: lower.to_vec(),
            col_upper: upper.to_vec(),
            integer: vec![integer; n],
        }
    }

    #[test]
    fn knapsack_row_fixes_heavy_items() {
        // 5a + 3b + 2c ≤ 6 with b = 1 fixed ⇒ a ≤ 0.6 ⇒ a = 0 (integer).
        let inf = f64::INFINITY;
        let p = problem(
            &[(&[5.0, 3.0, 2.0], -inf, 6.0)],
            &[0.0, 1.0, 0.0],
            &[1.0, 1.0, 1.0],
            true,
        );
        let mut prop = Propagator::new(&p, 1e-6);
        let (mut lo, mut hi) = (p.col_lower.clone(), p.col_upper.clone());
        assert!(prop.propagate(&mut lo, &mut hi, None, inf));
        assert_eq!(hi, vec![0.0, 1.0, 1.0]);
    }

    #[test]
    fn chained_equalities_propagate_and_detect_infeasibility() {
        // x0 = x1, x1 = x2, x0 ≥ 2, x2 ≤ 1 ⇒ infeasible through the chain.
        let p = problem(
            &[(&[1.0, -1.0, 0.0], 0.0, 0.0), (&[0.0, 1.0, -1.0], 0.0, 0.0)],
            &[2.0, 0.0, 0.0],
            &[5.0, 5.0, 1.0],
            true,
        );
        let mut prop = Propagator::new(&p, 1e-6);
        let (mut lo, mut hi) = (p.col_lower.clone(), p.col_upper.clone());
        assert!(!prop.propagate(&mut lo, &mut hi, None, f64::INFINITY));
    }

    #[test]
    fn objective_cutoff_tightens_bounds() {
        // min x0 + x1, cutoff 1.5 ⇒ integer x0, x1 ≤ 1.
        let p = problem(&[], &[0.0, 0.0], &[9.0, 9.0], true);
        let mut prop = Propagator::new(&p, 1e-6);
        let (mut lo, mut hi) = (p.col_lower.clone(), p.col_upper.clone());
        assert!(prop.propagate(&mut lo, &mut hi, None, 1.5));
        assert_eq!(hi, vec![1.0, 1.0]);
    }

    #[test]
    fn continuous_tightening_needs_real_improvement() {
        // x0 + x1 ≤ 9.9 on [0,10]²: cuts 1% of x's domain ⇒ rejected.
        let inf = f64::INFINITY;
        let p = problem(
            &[(&[1.0, 1.0], -inf, 9.9)],
            &[0.0, 0.0],
            &[10.0, 10.0],
            false,
        );
        let mut prop = Propagator::new(&p, 1e-6);
        let (mut lo, mut hi) = (p.col_lower.clone(), p.col_upper.clone());
        assert!(prop.propagate(&mut lo, &mut hi, None, inf));
        assert_eq!(hi, vec![10.0, 10.0]);
    }
}
