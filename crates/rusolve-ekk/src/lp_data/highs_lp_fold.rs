//! LP folding by colour refinement (rusolve lever `lp_fold`).
//!
//! Sources: Grohe, Kersting, Mladenov, Selman, "Dimension reduction via colour refinement",
//! ESA 2014; Mladenov, Ahmadi, Kersting, "Lifted linear programming", AISTATS 2012; Bödi, Herr,
//! Joswig, "Algorithms for highly symmetric linear and integer programs", Math. Prog. 2013
//! (orbit averaging).
//!
//! An equitable partition of rows and columns (1-WL on the coefficient graph, started from
//! cost/bounds colours) lets the LP be replaced by a quotient LP whose variable `y_P` stands for
//! every `x_j`, `j in P`: row class `Q` becomes `sum_P (sum_{j in P} a_ij) y_P` (any `i in Q`),
//! cost `|P| c_j`. An optimal primal/dual pair of the quotient lifts to one of the original:
//! `x_j = y_P`, row dual `y_i = y_Q / |Q|`, reduced cost `d_j = d_P / |P|`.

use std::collections::HashMap;

use crate::prelude::*;

fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Hash key of a value, equal for values within about 1e-9 relative (cont1: orbit-symmetric
/// right-hand sides differ in the last digits).
fn fbits(x: f64) -> u64 {
    if x == 0.0 {
        0
    } else if !x.is_finite() {
        x.to_bits()
    } else {
        let e = x.abs().log10().floor();
        let q = (x / 10f64.powf(e - 3.0)).round() as i64;
        mix(q as u64 ^ ((e as i64 as u64) << 48))
    }
}

/// Result of a fold: class of every row / column and the quotient LP.
pub struct LpFold {
    pub row_class: Vec<u32>,
    pub col_class: Vec<u32>,
    pub num_row_class: usize,
    pub num_col_class: usize,
    pub row_class_size: Vec<u32>,
    pub col_class_size: Vec<u32>,
    pub rounds: usize,
    /// False when classes merged values that differ (tolerance merge): the quotient optimum is
    /// then only a starting point for the original.
    pub exact: bool,
    pub quotient: HighsLp,
}

/// Renumber signatures to dense class ids; returns the class count.
fn renumber(sig: &[u64], out: &mut [u32]) -> usize {
    let mut map: HashMap<u64, u32> = HashMap::with_capacity(sig.len() / 2 + 1);
    for (o, s) in out.iter_mut().zip(sig) {
        let next = map.len() as u32;
        *o = *map.entry(*s).or_insert(next);
    }
    map.len()
}

/// Colour refinement. Returns `None` when the final quotient would exceed `max_ratio` of
/// (rows + cols) (checked every round: the class count only grows) or on a round cap.
pub fn fold_lp(lp: &HighsLp, max_ratio: f64, budget_s: f64) -> Option<LpFold> {
    let t_fold = std::time::Instant::now();
    let m = lp.num_row_ as usize;
    let n = lp.num_col_ as usize;
    if m == 0 || n == 0 || !lp.a_matrix_.is_colwise() {
        return None;
    }
    let a = &lp.a_matrix_;
    let nnz = a.start_[n] as usize;
    // Row-wise copy (CSR).
    let mut rstart = vec![0usize; m + 1];
    for k in 0..nnz {
        rstart[a.index_[k] as usize + 1] += 1;
    }
    for i in 0..m {
        rstart[i + 1] += rstart[i];
    }
    let mut fill = rstart.clone();
    let mut ridx = vec![0u32; nnz];
    let mut rval = vec![0f64; nnz];
    for j in 0..n {
        for k in a.start_[j] as usize..a.start_[j + 1] as usize {
            let i = a.index_[k] as usize;
            ridx[fill[i]] = j as u32;
            rval[fill[i]] = a.value_[k];
            fill[i] += 1;
        }
    }
    // Initial colours.
    let mut rc = vec![0u32; m];
    let mut cc = vec![0u32; n];
    let rsig: Vec<u64> = (0..m)
        .map(|i| mix(fbits(lp.row_lower_[i]) ^ mix(fbits(lp.row_upper_[i]) ^ 0x5151)))
        .collect();
    let csig: Vec<u64> = (0..n)
        .map(|j| {
            mix(fbits(lp.col_cost_[j])
                ^ mix(fbits(lp.col_lower_[j]) ^ mix(fbits(lp.col_upper_[j]) ^ 0x7777)))
        })
        .collect();
    let mut nr = renumber(&rsig, &mut rc);
    let mut nc = renumber(&csig, &mut cc);
    let limit = max_ratio * (m + n) as f64;
    let mut rounds = 0usize;
    let mut sig_r = vec![0u64; m];
    let mut sig_c = vec![0u64; n];
    loop {
        // rusolve: fhnw-binschedule1 spent 239 s folding at TL 30 s
        if (nr + nc) as f64 > limit || rounds > 1000 || t_fold.elapsed().as_secs_f64() > budget_s {
            return None;
        }
        rounds += 1;
        for i in 0..m {
            let mut h = 0u64;
            for k in rstart[i]..rstart[i + 1] {
                h = h.wrapping_add(mix(
                    ((cc[ridx[k] as usize] as u64) << 32 | 1) ^ mix(fbits(rval[k]))
                ));
            }
            sig_r[i] = mix(h ^ mix(rc[i] as u64 + 0x1234_5678));
        }
        let nr2 = renumber(&sig_r, &mut rc);
        for j in 0..n {
            let mut h = 0u64;
            for k in a.start_[j] as usize..a.start_[j + 1] as usize {
                h = h.wrapping_add(mix(
                    ((rc[a.index_[k] as usize] as u64) << 32 | 2) ^ mix(fbits(a.value_[k]))
                ));
            }
            sig_c[j] = mix(h ^ mix(cc[j] as u64 + 0x8765_4321));
        }
        let nc2 = renumber(&sig_c, &mut cc);
        let stable = nr2 == nr && nc2 == nc;
        nr = nr2;
        nc = nc2;
        if stable {
            break;
        }
    }
    if (nr + nc) as f64 > limit {
        return None;
    }
    // Class sizes and representatives.
    let mut rsize = vec![0u32; nr];
    let mut rrep = vec![usize::MAX; nr];
    for i in 0..m {
        let q = rc[i] as usize;
        rsize[q] += 1;
        if rrep[q] == usize::MAX {
            rrep[q] = i;
        }
    }
    let mut csize = vec![0u32; nc];
    let mut crep = vec![usize::MAX; nc];
    for j in 0..n {
        let p = cc[j] as usize;
        csize[p] += 1;
        if crep[p] == usize::MAX {
            crep[p] = j;
        }
    }
    // Quotient LP (column-wise): entry (Q,P) = sum_{j in P} a_{rep(Q), j}.
    let mut rows_of_col: Vec<Vec<(u32, f64)>> = vec![Vec::new(); nc];
    let mut acc: HashMap<u32, f64> = HashMap::new();
    for q in 0..nr {
        let i = rrep[q];
        acc.clear();
        for k in rstart[i]..rstart[i + 1] {
            *acc.entry(cc[ridx[k] as usize]).or_insert(0.0) += rval[k];
        }
        for (&p, &v) in acc.iter() {
            if v != 0.0 {
                rows_of_col[p as usize].push((q as u32, v));
            }
        }
    }
    let mut start = vec![0 as HighsInt; nc + 1];
    let mut index = Vec::new();
    let mut value = Vec::new();
    for p in 0..nc {
        rows_of_col[p].sort_by_key(|e| e.0);
        for &(q, v) in &rows_of_col[p] {
            index.push(q as HighsInt);
            value.push(v);
        }
        start[p + 1] = index.len() as HighsInt;
    }
    let mut q = HighsLp::default();
    q.num_col_ = nc as HighsInt;
    q.num_row_ = nr as HighsInt;
    q.sense_ = lp.sense_;
    q.offset_ = lp.offset_;
    q.model_name_ = lp.model_name_.clone();
    q.col_cost_ = (0..nc)
        .map(|p| lp.col_cost_[crep[p]] * csize[p] as f64)
        .collect();
    q.col_lower_ = (0..nc).map(|p| lp.col_lower_[crep[p]]).collect();
    q.col_upper_ = (0..nc).map(|p| lp.col_upper_[crep[p]]).collect();
    q.row_lower_ = (0..nr).map(|r| lp.row_lower_[rrep[r]]).collect();
    q.row_upper_ = (0..nr).map(|r| lp.row_upper_[rrep[r]]).collect();
    q.a_matrix_.format_ = MatrixFormat::kColwise;
    q.a_matrix_.num_col_ = nc as HighsInt;
    q.a_matrix_.num_row_ = nr as HighsInt;
    q.a_matrix_.start_ = start;
    q.a_matrix_.index_ = index;
    q.a_matrix_.value_ = value;
    // Exactness: every member must carry its representative's data.
    let close = |a: f64, b: f64| a == b || (a - b).abs() <= 1e-13 * a.abs().max(b.abs());
    let mut exact = (0..n).all(|j| {
        let r = crep[cc[j] as usize];
        close(lp.col_cost_[j], lp.col_cost_[r])
            && close(lp.col_lower_[j], lp.col_lower_[r])
            && close(lp.col_upper_[j], lp.col_upper_[r])
    }) && (0..m).all(|i| {
        let r = rrep[rc[i] as usize];
        close(lp.row_lower_[i], lp.row_lower_[r]) && close(lp.row_upper_[i], lp.row_upper_[r])
    });
    if exact {
        let mut qrow: Vec<Vec<(u32, f64)>> = vec![Vec::new(); nr];
        for (p, col) in rows_of_col.iter().enumerate() {
            for &(r, v) in col {
                qrow[r as usize].push((p as u32, v));
            }
        }
        for row in qrow.iter_mut() {
            row.sort_by_key(|e| e.0);
        }
        let mut cur: Vec<(u32, f64)> = Vec::new();
        for i in 0..m {
            acc.clear();
            for k in rstart[i]..rstart[i + 1] {
                *acc.entry(cc[ridx[k] as usize]).or_insert(0.0) += rval[k];
            }
            cur.clear();
            cur.extend(acc.iter().filter(|e| *e.1 != 0.0).map(|(&p, &v)| (p, v)));
            cur.sort_by_key(|e| e.0);
            let rep = &qrow[rc[i] as usize];
            if cur.len() != rep.len()
                || cur
                    .iter()
                    .zip(rep)
                    .any(|(a, b)| a.0 != b.0 || !close(a.1, b.1))
            {
                exact = false;
                break;
            }
        }
    }
    Some(LpFold {
        exact,
        row_class: rc,
        col_class: cc,
        num_row_class: nr,
        num_col_class: nc,
        row_class_size: rsize,
        col_class_size: csize,
        rounds,
        quotient: q,
    })
}

impl LpFold {
    /// Lift a quotient solution to the original LP (see module doc).
    pub fn lift(&self, sol: &HighsSolution) -> HighsSolution {
        let m = self.row_class.len();
        let n = self.col_class.len();
        let mut out = HighsSolution::default();
        out.col_value = (0..n)
            .map(|j| sol.col_value[self.col_class[j] as usize])
            .collect();
        out.row_value = (0..m)
            .map(|i| sol.row_value[self.row_class[i] as usize])
            .collect();
        out.col_dual = (0..n)
            .map(|j| {
                let p = self.col_class[j] as usize;
                sol.col_dual[p] / self.col_class_size[p] as f64
            })
            .collect();
        out.row_dual = (0..m)
            .map(|i| {
                let q = self.row_class[i] as usize;
                sol.row_dual[q] / self.row_class_size[q] as f64
            })
            .collect();
        out.value_valid = true;
        out.dual_valid = true;
        out
    }
}

/// Largest violation of the row/column bounds by `col_value` (absolute).
pub fn max_primal_violation(lp: &HighsLp, x: &[f64]) -> f64 {
    let m = lp.num_row_ as usize;
    let n = lp.num_col_ as usize;
    let a = &lp.a_matrix_;
    let mut act = vec![0.0; m];
    let mut worst: f64 = 0.0;
    for j in 0..n {
        worst = worst
            .max(lp.col_lower_[j] - x[j])
            .max(x[j] - lp.col_upper_[j]);
        for k in a.start_[j] as usize..a.start_[j + 1] as usize {
            act[a.index_[k] as usize] += a.value_[k] * x[j];
        }
    }
    for i in 0..m {
        worst = worst
            .max(lp.row_lower_[i] - act[i])
            .max(act[i] - lp.row_upper_[i]);
    }
    worst
}

impl Highs {
    /// Fold the LP, solve the quotient, lift and cross over on the original (lever `lp_fold`).
    /// Returns `None` to fall back to the normal path (not an LP, nothing to fold, quotient not
    /// optimal, or the lift failed the primal check).
    pub(crate) fn try_lp_fold(&mut self) -> Option<HighsStatus> {
        let status = self.try_lp_fold_inner()?;
        // Crossover / the lifted point leave the iteration counts at -1; reports expect >= 0.
        self.info_.simplex_iteration_count = self.info_.simplex_iteration_count.max(0);
        self.info_.ipm_iteration_count = self.info_.ipm_iteration_count.max(0);
        self.info_.crossover_iteration_count = self.info_.crossover_iteration_count.max(0);
        Some(status)
    }

    fn try_lp_fold_inner(&mut self) -> Option<HighsStatus> {
        if !self.options_.lp_fold
            || self.model_lp_.is_mip()
            || self.model_hessian_.dim_ != 0
            || self.model_lp_.has_infinite_cost_
            || self.basis_.valid
            || self.solution_.value_valid
        {
            return None;
        }
        self.model_lp_.ensure_colwise();
        let t0 = std::time::Instant::now();
        let t_start = t0;
        let (m, n) = (self.model_lp_.num_row_, self.model_lp_.num_col_);
        // ratio 0.6 (was 0.8): irish-electricity folds to 0.73 and its quotient solve (44 s, status
        // Unknown, primal infeasibility 3e-5) was pure waste; every fold that paid off is <= 0.51.
        // A fold may use at most a fifth of the user's time limit
        let fold = fold_lp(&self.model_lp_, 0.6, 0.2 * self.options_.time_limit);
        let dt = t0.elapsed().as_secs_f64();
        let fold = match fold {
            Some(f) => f,
            None => {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "LP fold: no useful fold of {} rows, {} cols ({:.2}s)\n",
                    m,
                    n,
                    dt
                );
                return None;
            }
        };
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "LP fold: {} rows, {} cols -> {} rows, {} cols in {} rounds ({:.2}s){}\n",
            m,
            n,
            fold.num_row_class,
            fold.num_col_class,
            fold.rounds,
            dt,
            if fold.exact {
                ""
            } else {
                " (tolerance merge, not used)"
            }
        );
        if !fold.exact {
            // Tolerance merge (cont1/cont11 reach a 20500-row fold only this way): the
            // lifted point is not complementary on the original and IPX crossover rejects it
            // so report the size and solve normally.
            return None;
        }
        let mut opts = self.options_.clone();
        opts.lp_fold = false;
        // ex10 quotient: dualised simplex 15 s vs 5 s
        opts.lp_dualize_auto = false;
        opts.lp_solver_auto = false;
        // the quotient solve gets what is left of the user's limit, not a fresh one
        if opts.time_limit.is_finite() {
            opts.time_limit = (opts.time_limit - t_start.elapsed().as_secs_f64()).max(0.0);
        }
        let mut sub = Highs::with_options(opts);
        sub.pass_model(fold.quotient.clone());
        let sub_run = sub.run();
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "LP fold: quotient solve ended with status {} (run {:?}, value {}, dual {}, simplex its {}, obj {}, pinf {:.2e}, dinf {:.2e}, pstat {}, dstat {}) at {:.2}s\n",
            self.model_status_to_string(sub.get_model_status()),
            sub_run,
            sub.get_solution().value_valid,
            sub.get_solution().dual_valid,
            sub.info_.simplex_iteration_count,
            sub.info_.objective_function_value,
            sub.info_.max_primal_infeasibility,
            sub.info_.max_dual_infeasibility,
            sub.info_.primal_solution_status,
            sub.info_.dual_solution_status,
            t_start.elapsed().as_secs_f64()
        );
        if sub.get_model_status() != HighsModelStatus::kOptimal {
            return None;
        }
        let qs = sub.get_solution();
        if !qs.value_valid || !qs.dual_valid {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "LP fold: quotient solution invalid (value {}, dual {})\n",
                qs.value_valid,
                qs.dual_valid
            );
            return None;
        }
        let lifted = fold.lift(qs);
        let viol = max_primal_violation(&self.model_lp_, &lifted.col_value);
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "LP fold: quotient optimal ({} simplex, {} ipm iterations), lifted primal violation {:.2e}\n",
            sub.info_.simplex_iteration_count,
            sub.info_.ipm_iteration_count,
            viol
        );
        if viol > 1e-6 && fold.exact {
            return None;
        }
        if self.options_.lp_fold_basis && sub.get_basis().valid {
            let qb = sub.get_basis();
            let mut basis = HighsBasis::default();
            basis.col_status = fold
                .col_class
                .iter()
                .map(|&p| qb.col_status[p as usize])
                .collect();
            basis.row_status = fold
                .row_class
                .iter()
                .map(|&q| qb.row_status[q as usize])
                .collect();
            basis.alien = true;
            self.clear_solver();
            self.set_basis(&basis, "lp_fold");
            let status = self.optimize_highs();
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "LP fold: simplex from the lifted quotient basis, {} iterations\n",
                self.info_.simplex_iteration_count
            );
            return Some(status);
        }
        if self.options_.run_crossover != kHighsOffString || !fold.exact {
            // Bounded effort: a crossover that costs far more than the fold saved is not worth it.
            let cap = 10.0 + 10.0 * t_start.elapsed().as_secs_f64();
            let user_limit = self.options_.time_limit;
            self.options_.time_limit = user_limit.min(cap);
            let status = self.crossover(&lifted);
            self.options_.time_limit = user_limit;
            let ms = self.model_status_;
            if ms == HighsModelStatus::kOptimal {
                return Some(status);
            }
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "LP fold: crossover ended with status {} (cap {:.0}s)\n",
                self.model_status_to_string(ms),
                cap
            );
            if !fold.exact {
                self.clear_solver();
                return None;
            }
        }
        // No basis wanted (or crossover failed): the lifted optimal primal/dual pair is the answer.
        self.clear_solver();
        self.solution_ = lifted;
        self.basis_.clear();
        self.model_status_ = HighsModelStatus::kOptimal;
        self.info_.objective_function_value = self
            .model_lp_
            .objective_value(&self.solution_.col_value, self.options_.h16_obj_cdouble);
        crate::lp_data::presolve_support::get_lp_kkt_failures(
            &self.options_,
            &self.model_lp_,
            &self.solution_,
            &mut self.info_,
            false,
        );
        Some(self.return_from_highs(HighsStatus::kOk))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cycle cover: min sum x s.t. x_j + x_{j+1} >= 1 on a 6-cycle, plus a distinguished
    /// column (cost 2, bound 3) to check classes stay apart. Folds to 2 cols + 2 rows.
    fn cycle_lp() -> HighsLp {
        let n = 7;
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = 7;
        lp.col_cost_ = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0];
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0];
        lp.row_lower_ = vec![1.0; 7];
        lp.row_upper_ = vec![kHighsInf; 7];
        let mut cols: Vec<Vec<(HighsInt, f64)>> = vec![Vec::new(); n];
        for r in 0..6 {
            cols[r].push((r as HighsInt, 1.0));
            cols[(r + 1) % 6].push((r as HighsInt, 1.0));
        }
        cols[6].push((6 as HighsInt, 1.0));
        let mut start: Vec<HighsInt> = vec![0];
        let mut index: Vec<HighsInt> = Vec::new();
        let mut value: Vec<f64> = Vec::new();
        for col in cols.iter_mut() {
            col.sort_by_key(|e| e.0);
            for &(i, v) in col.iter() {
                index.push(i);
                value.push(v);
            }
            start.push(index.len() as HighsInt);
        }
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = 7;
        lp.a_matrix_.start_ = start;
        lp.a_matrix_.index_ = index;
        lp.a_matrix_.value_ = value;
        lp
    }

    #[test]
    fn detects_cycle_symmetry() {
        let f = fold_lp(&cycle_lp(), 0.8, f64::INFINITY).expect("fold");
        assert_eq!((f.num_col_class, f.num_row_class), (2, 2));
        assert_eq!(f.quotient.col_cost_, vec![6.0, 2.0]);
    }

    #[test]
    fn fold_matches_unfolded_optimum() {
        let solve = |fold: bool, crossover: &str| {
            let mut h = Highs::default();
            h.options_.output_flag = false;
            h.options_.lp_fold = fold;
            h.options_.run_crossover = crossover.to_string();
            h.pass_model(cycle_lp());
            h.run();
            assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
            (h.info_.objective_function_value, h.get_solution().clone())
        };
        let (base, _) = solve(false, "choose");
        for crossover in ["off", "on"] {
            let (obj, sol) = solve(true, crossover);
            assert!((obj - base).abs() < 1e-8, "{obj} vs {base}");
            assert!(max_primal_violation(&cycle_lp(), &sol.col_value) < 1e-8);
        }
    }
}
