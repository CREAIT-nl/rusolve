//! rusolve (not a HiGHS rule): `mip_presolve_eq_subst`.
//!
//! Substitution of a unit-coefficient integer column out of a short all-integer equation when
//! the column is not implied free: the bounds the rest of the equation does not already
//! imply are moved into one new row on the remaining columns, then the column is relaxed to
//! free and substituted by the existing HiGHS free-column substitution
//! (`substitute_free_col`, the same postsolve record the aggregator writes).
//!
//! Sources. Implied-free column substitution: Andersen and Andersen, "Presolving in linear
//! programming", Math. Programming 71, 1995; Achterberg, Bixby, Gu, Rothberg, Weninger,
//! "Presolve reductions in mixed integer programming", INFORMS J. Computing 32(2), 2020
//! (aggregation). The extension here, keeping the non-implied bound sides as one row on the
//! remaining columns so the column becomes free, is plain algebra and is not in SCIP (its
//! `convertLongEquality` in cons_linear.c, Apache-2.0, multi-aggregates only when the bounds are
//! implied). The integrality test is the one SCIP uses there: the right-hand side divided by
//! the coefficient and every residual coefficient integral over integer columns. With every
//! other column integer, every coefficient integer, the right-hand side integer and
//! |a_k| = 1, the substituted column is integer valued, so dropping its integrality loses
//! nothing. No SCIP code is translated.
//!
//! Needs the row-append infrastructure (`mip_presolve_row_append`).

use super::*;
use crate::prelude::*;

/// Longest equation considered (columns).
const EQ_SUBST_MAX_ROW: HighsInt = 4;
/// Longest column considered (rows): bounds the fill-in of the substitution.
const EQ_SUBST_MAX_COL: HighsInt = 3;

impl HPresolve {
    /// Returns kOk; the number of columns removed is visible through `num_deleted_cols`.
    pub fn eq_subst_bound_rows(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        if !self.options.mip_presolve_eq_subst
            || !self.options.mip_presolve_row_append
            || self.mipsolver.is_none()
        {
            return HPresolveResult::kOk;
        }
        let frac = |x: f64| (x - x.round()).abs();
        let small = self.options.small_matrix_value;
        let tol = self.primal_feastol;
        let initial_rows = self.model.num_row_;

        for row in 0..initial_rows {
            let ru = row as usize;
            if self.row_deleted[ru] != 0 || !self.is_equation(row) {
                continue;
            }
            let rs = self.rowsize[ru];
            if !(3..=EQ_SUBST_MAX_ROW).contains(&rs) || self.rowsize_integer[ru] != rs {
                continue;
            }
            let rhs = self.model.row_upper_[ru];
            if frac(rhs) > tol {
                continue;
            }
            let entries: Vec<(HighsInt, f64)> = self
                .get_row_vector(row)
                .iter(self)
                .map(|nz| (nz.index(), nz.value()))
                .collect();
            if entries.iter().any(|&(_, vv)| frac(vv) > small) {
                continue;
            }
            // pick the unit-coefficient column with the shortest column
            let mut best: Option<(HighsInt, f64)> = None;
            for &(cc, vv) in &entries {
                if (vv.abs() - 1.0).abs() > small || self.colsize[cc as usize] > EQ_SUBST_MAX_COL {
                    continue;
                }
                if best.is_none_or(|(b, _)| self.colsize[cc as usize] < self.colsize[b as usize]) {
                    best = Some((cc, vv));
                }
            }
            let Some((k, ak)) = best else { continue };
            if self.col_deleted[k as usize] != 0 || self.is_implied_free(k) {
                continue;
            }
            let (lk, uk) = (
                self.model.col_lower_[k as usize],
                self.model.col_upper_[k as usize],
            );

            // activity range of S = sum of the other terms
            let (mut smin, mut smax) = (0.0_f64, 0.0_f64);
            for &(cc, vv) in &entries {
                if cc == k {
                    continue;
                }
                let (l, u) = (
                    self.model.col_lower_[cc as usize],
                    self.model.col_upper_[cc as usize],
                );
                let (a, b) = if vv > 0.0 {
                    (vv * l, vv * u)
                } else {
                    (vv * u, vv * l)
                };
                smin += a;
                smax += b;
            }
            // x_k = ak * (rhs - S); the bounds of x_k as bounds on S
            let (mut lo, mut hi) = if ak > 0.0 {
                (rhs - uk, rhs - lk)
            } else {
                (rhs + lk, rhs + uk)
            };
            if !(lo > smin + tol) {
                lo = -kHighsInf;
            }
            if !(hi < smax - tol) {
                hi = kHighsInf;
            }
            if lo == -kHighsInf && hi == kHighsInf {
                continue;
            }

            let mut idx: Vec<HighsInt> = Vec::new();
            let mut val: Vec<f64> = Vec::new();
            for &(cc, vv) in &entries {
                if cc != k {
                    idx.push(cc);
                    val.push(vv);
                }
            }
            if !self.add_to_matrix_row(postsolve_stack, lo, hi, &idx, &val) {
                return HPresolveResult::kStopped;
            }
            hpresolve_checked_call!(self.change_col_bounds(k, -kHighsInf, kHighsInf));
            self.store_row(row);
            self.substitute_free_col(postsolve_stack, row, k, false);
            hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));
            hpresolve_checked_call!(self.check_limits(postsolve_stack));
        }
        HPresolveResult::kOk
    }
}

#[cfg(test)]
mod tests {
    use crate::lp_data::hconst::{HighsModelStatus, HighsVarType, ObjSense};
    use crate::lp_data::highs_run::Highs;
    use crate::prelude::*;

    const G: usize = 3; // gadgets
    const M: f64 = 3.0;

    /// Per gadget i: integers x_i, y_i in [0,2]; binaries b0, b1, b2 with
    /// x - y - M b1 <= -1, y - x - M b2 <= -1, -b0 + b1 + b2 = 1 (b0 = [x == y] style
    /// reification pair, the gfd gadget); one coupling row sum x >= 3; non-integral costs.
    /// Columns: x,y,b0,b1,b2 per gadget.
    fn build(on: bool) -> (Highs, Vec<Vec<(usize, f64)>>, Vec<(f64, f64)>, Vec<f64>) {
        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        let n = 5 * G;
        let mut cost = vec![0.0; n];
        for g in 0..G {
            let o = 5 * g;
            cost[o] = 1.3 + g as f64; // x
            cost[o + 1] = 0.7; // y
            cost[o + 2] = -2.1 + 0.4 * g as f64; // b0
        }
        for j in 0..n {
            let (ub, int) = (if j % 5 < 2 { 2.0 } else { 1.0 }, HighsVarType::kInteger);
            h.add_vars(1, &[0.0], &[ub]);
            h.change_col_cost(j as HighsInt, cost[j]);
            h.change_col_integrality(j as HighsInt, int);
        }
        let mut rows: Vec<Vec<(usize, f64)>> = Vec::new();
        let mut bounds: Vec<(f64, f64)> = Vec::new();
        for g in 0..G {
            let o = 5 * g;
            rows.push(vec![(o, 1.0), (o + 1, -1.0), (o + 3, -M)]);
            bounds.push((-kHighsInf, -1.0));
            rows.push(vec![(o + 1, 1.0), (o, -1.0), (o + 4, -M)]);
            bounds.push((-kHighsInf, -1.0));
            rows.push(vec![(o + 2, -1.0), (o + 3, 1.0), (o + 4, 1.0)]);
            bounds.push((1.0, 1.0));
        }
        rows.push((0..G).map(|g| (5 * g, 1.0)).collect());
        bounds.push((3.0, kHighsInf));
        for (r, row) in rows.iter().enumerate() {
            let idx: Vec<HighsInt> = row.iter().map(|&(j, _)| j as HighsInt).collect();
            let val: Vec<f64> = row.iter().map(|&(_, v)| v).collect();
            h.add_row(bounds[r].0, bounds[r].1, idx.len() as HighsInt, &idx, &val);
        }
        h.change_objective_sense(ObjSense::kMinimize);
        h.set_option_value_bool("mip_presolve_row_append", on);
        h.set_option_value_bool("mip_presolve_eq_subst", on);
        (h, rows, bounds, cost)
    }

    #[test]
    fn eq_subst_postsolve_is_feasible_and_optimal() {
        let mut objs = Vec::new();
        for on in [false, true] {
            let (mut h, rows, bounds, cost) = build(on);
            let _ = h.run();
            assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
            let x = h.get_solution().col_value.clone();
            assert_eq!(x.len(), cost.len());
            for (j, v) in x.iter().enumerate() {
                assert!((v - v.round()).abs() < 1e-6, "x[{j}] = {v} not integral");
                let ub = if j % 5 < 2 { 2.0 } else { 1.0 };
                assert!(*v >= -1e-6 && *v <= ub + 1e-6, "x[{j}] = {v} out of bounds");
            }
            for (r, row) in rows.iter().enumerate() {
                let a: f64 = row.iter().map(|&(j, v)| v * x[j]).sum();
                assert!(
                    a >= bounds[r].0 - 1e-6 && a <= bounds[r].1 + 1e-6,
                    "on={on} row {r} violated: {a}"
                );
            }
            let obj: f64 = cost.iter().zip(&x).map(|(cj, v)| cj * v).sum();
            assert!((obj - h.get_info().objective_function_value).abs() < 1e-6);
            objs.push(obj);
        }
        assert!((objs[0] - objs[1]).abs() < 1e-6, "{objs:?}");
    }
}
