//! rusolve (not a HiGHS rule): `mip_presolve_xor_gf2`.
//!
//! Parity (XOR) reasoning over GF(2). A linear equation `sum_j s_j x_j + 2 s y = r` with binary
//! `x_j`, unit coefficients `s_j = +-1`, one integer column `y` with `|coef| = 2` that appears
//! in this row only, and integral `r`, holds only if `sum_j x_j = r (mod 2)`. All such rows form a
//! linear system over GF(2). Gaussian elimination decides infeasibility (a row `0 = 1`), and
//! fixes every binary that the reduced system determines (a reduced row with one variable).
//!
//! Sources. SCIP 10 (Apache-2.0) `cons_xor.c`: `linconsUpgdXor` (the row pattern and the
//! "parity column occurs in this row only" test) and `checkSystemGF2` /
//! `computeRowEchelonGF2` (Gaussian elimination over GF(2); infeasible system is a cutoff, a
//! full-rank system fixes every variable; size caps `MAXXORCONSSSYSTEM` and `MAXXORVARSSYSTEM`
//! are both 1000). The row-reduced form here additionally fixes variables of rank-deficient
//! systems (a reduced row with a single variable), which follows from the same equations and is
//! not in SCIP. Motivation: SCIP solves enlight_hard at the root through this propagation and
//! needs 800+ nodes with it disabled. Re-implemented from
//! the algorithm description, no SCIP code translated verbatim.

use super::*;
use crate::prelude::*;

/// Largest system (rows and distinct binaries), the SCIP caps.
const XOR_MAX_ROWS: usize = 1000;
const XOR_MAX_VARS: usize = 1000;

impl HPresolve {
    pub fn xor_gf2_presolve(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        if !self.options.mip_presolve_xor_gf2 || self.mipsolver.is_none() {
            return HPresolveResult::kOk;
        }
        let tol = self.primal_feastol;
        let small = self.options.small_matrix_value;
        let frac = |x: f64| (x - x.round()).abs();

        // collect XOR rows: (binary columns, parity of the right-hand side)
        let mut xor_rows: Vec<(Vec<HighsInt>, bool)> = Vec::new();
        'rows: for row in 0..self.model.num_row_ {
            let ru = row as usize;
            if self.row_deleted[ru] != 0 || !self.is_equation(row) {
                continue;
            }
            let rhs = self.model.row_upper_[ru];
            if !rhs.is_finite() || frac(rhs) > tol {
                continue;
            }
            let mut cols: Vec<HighsInt> = Vec::new();
            let mut has_parity_col = false;
            let mut fixed_ones = 0i64;
            let entries: Vec<(HighsInt, f64)> = self
                .get_row_vector(row)
                .iter(self)
                .map(|nz| (nz.index(), nz.value()))
                .collect();
            for &(cj, coef) in &entries {
                let cu = cj as usize;
                if self.model.integrality_[cu] == HighsVarType::kContinuous {
                    continue 'rows;
                }
                if (coef.abs() - 2.0).abs() <= small {
                    if has_parity_col || self.colsize[cu] != 1 {
                        continue 'rows;
                    }
                    has_parity_col = true;
                } else if (coef.abs() - 1.0).abs() <= small {
                    let (lb, ub) = (self.model.col_lower_[cu], self.model.col_upper_[cu]);
                    if lb >= 0.5 && ub <= 1.5 {
                        fixed_ones += 1; // fixed at one
                    } else if lb <= -0.5 || ub >= 1.5 {
                        continue 'rows; // not binary
                    } else if ub >= 0.5 {
                        cols.push(cj); // free binary
                    }
                } else {
                    continue 'rows;
                }
            }
            if !has_parity_col {
                continue;
            }
            let parity = ((rhs.round() as i64) + fixed_ones).rem_euclid(2) == 1;
            xor_rows.push((cols, parity));
        }
        if xor_rows.is_empty() {
            return HPresolveResult::kOk;
        }

        // columns of the system
        let mut col_index = std::collections::HashMap::<HighsInt, usize>::new();
        let mut col_of: Vec<HighsInt> = Vec::new();
        for (cols, _) in &xor_rows {
            for &cj in cols {
                col_index.entry(cj).or_insert_with(|| {
                    col_of.push(cj);
                    col_of.len() - 1
                });
            }
        }
        let (m, n) = (xor_rows.len(), col_of.len());
        if m > XOR_MAX_ROWS || n > XOR_MAX_VARS {
            return HPresolveResult::kOk;
        }
        let words = n.div_ceil(64) + 1; // last bit of the row holds the right-hand side
        let rhs_bit = n;
        let mut a = vec![vec![0u64; words]; m];
        for (i, (cols, parity)) in xor_rows.iter().enumerate() {
            for &cj in cols {
                let j = col_index[&cj];
                a[i][j / 64] ^= 1u64 << (j % 64);
            }
            if *parity {
                a[i][rhs_bit / 64] ^= 1u64 << (rhs_bit % 64);
            }
        }
        let bit = |r: &Vec<u64>, j: usize| (r[j / 64] >> (j % 64)) & 1 == 1;

        // reduced row echelon form over GF(2)
        let mut rank = 0usize;
        for j in 0..n {
            let Some(p) = (rank..m).find(|&i| bit(&a[i], j)) else {
                continue;
            };
            a.swap(rank, p);
            let pivot = a[rank].clone();
            for (i, row) in a.iter_mut().enumerate() {
                if i != rank && (row[j / 64] >> (j % 64)) & 1 == 1 {
                    for (w, pw) in row.iter_mut().zip(&pivot) {
                        *w ^= *pw;
                    }
                }
            }
            rank += 1;
            if rank == m {
                break;
            }
        }
        // a zero row with right-hand side one: infeasible
        for row in a.iter().skip(rank) {
            if bit(row, rhs_bit) {
                return HPresolveResult::kPrimalInfeasible;
            }
        }
        // reduced rows with a single variable fix it
        let mut fixings: Vec<(HighsInt, f64)> = Vec::new();
        for row in a.iter().take(rank) {
            let mut only: Option<usize> = None;
            let mut count = 0;
            for j in 0..n {
                if bit(row, j) {
                    count += 1;
                    only = Some(j);
                    if count > 1 {
                        break;
                    }
                }
            }
            if count == 1 {
                let j = only.unwrap();
                fixings.push((col_of[j], if bit(row, rhs_bit) { 1.0 } else { 0.0 }));
            }
        }
        crate::highs_log_user!(
            &self.options.log_options,
            HighsLogType::kInfo,
            "XOR parity rows: {} rows, {} columns, rank {}, {} fixings\n",
            m,
            n,
            rank,
            fixings.len()
        );
        for (cj, val) in fixings {
            if self.col_deleted[cj as usize] != 0 {
                continue;
            }
            let (lb, ub) = (
                self.model.col_lower_[cj as usize],
                self.model.col_upper_[cj as usize],
            );
            if val < lb - tol || val > ub + tol {
                return HPresolveResult::kPrimalInfeasible;
            }
            hpresolve_checked_call!(self.change_col_bounds(cj, val, val));
            hpresolve_checked_call!(self.fix_col_to_lower(postsolve_stack, cj));
        }
        hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));
        HPresolveResult::kOk
    }
}

impl HPresolve {
    /// rusolve `mip_ti_precedence_sepa` / `mip_ti_list_sched`: the
    /// levers keep presolve from sparsifying, because sparsification rewrites the machine-slot
    /// cliques they separate over. Only models that look time-indexed are protected: at least ten
    /// disjoint assignment rows (`sum x = 1`, unit coefficients, binaries) covering 90 % of the
    /// integer columns, and, over those columns, a machine-slot row (`sum x <= 1` across two or more
    /// sets), a precedence row (one-sided, two sets) and a makespan link (one column outside the sets).
    pub fn ti_structure_protect(&self) -> bool {
        // the separator reads the structure from the main problem only (sub-MIPs shift row
        // sides), so a sub-MIP presolve keeps sparsify
        if !(self.options.mip_ti_precedence_sepa || self.options.mip_ti_list_sched)
            || self.mipsolver.is_none()
            || self.mip().submip
        {
            return false;
        }
        let ncol = self.model.num_col_ as usize;
        let n_int = (0..ncol)
            .filter(|&cj| self.model.integrality_[cj] != HighsVarType::kContinuous)
            .count();
        let is_unit_bin = |cu: usize, v: f64| {
            v == 1.0
                && self.model.integrality_[cu] != HighsVarType::kContinuous
                && self.model.col_lower_[cu] == 0.0
                && self.model.col_upper_[cu] == 1.0
        };
        // assignment sets
        let mut set_of = vec![-1i32; ncol];
        let (mut sets, mut members) = (0i32, 0usize);
        for row in 0..self.model.num_row_ {
            let ru = row as usize;
            if self.row_deleted[ru] != 0
                || self.model.row_lower_[ru] != 1.0
                || self.model.row_upper_[ru] != 1.0
            {
                continue;
            }
            let ent: Vec<(usize, f64)> = self
                .get_row_vector(row)
                .iter(self)
                .map(|nz| (nz.index() as usize, nz.value()))
                .collect();
            if ent.len() < 2
                || ent
                    .iter()
                    .any(|&(cu, v)| !is_unit_bin(cu, v) || set_of[cu] >= 0)
            {
                continue;
            }
            for &(cu, _) in &ent {
                set_of[cu] = sets;
            }
            sets += 1;
            members += ent.len();
        }
        if sets < 10 || (members + 4) * 10 < n_int * 9 {
            return false;
        }
        let (mut slots, mut precs, mut links) = (0usize, 0usize, 0usize);
        for row in 0..self.model.num_row_ {
            let ru = row as usize;
            if self.row_deleted[ru] != 0 {
                continue;
            }
            let (rl, rup) = (self.model.row_lower_[ru], self.model.row_upper_[ru]);
            let ent: Vec<(usize, f64)> = self
                .get_row_vector(row)
                .iter(self)
                .map(|nz| (nz.index() as usize, nz.value()))
                .collect();
            if ent.len() < 2 {
                continue;
            }
            // columns outside every assignment set (the makespan, continuous or integer typed)
            let ncont = ent.iter().filter(|&&(cu, _)| set_of[cu] < 0).count();
            if ncont == 1 && (rl.is_finite() != rup.is_finite()) {
                links += 1;
            } else if ncont == 0 && ent.iter().all(|&(cu, _)| set_of[cu] >= 0) {
                let mut ids: Vec<i32> = ent.iter().map(|&(cu, _)| set_of[cu]).collect();
                ids.sort_unstable();
                ids.dedup();
                if rl <= -kHighsInf && rup == 1.0 && ent.iter().all(|&(_, v)| v == 1.0) {
                    if ids.len() >= 2 {
                        slots += 1;
                    }
                } else if (rl.is_finite() != rup.is_finite()) && ids.len() == 2 {
                    precs += 1;
                }
            }
            if slots > 0 && precs > 0 && links > 0 {
                return true;
            }
        }
        false
    }
}
