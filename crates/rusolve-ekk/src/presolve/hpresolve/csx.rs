//! HiGHS source: presolve/HPresolve.cpp (lines 2559-3034)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! CSC/CSR conversion, fill-in, column transform and substitution. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::okFromCSC (lines 2559-2620, 62 lines)
    ///
    /// Callers pass vectors taken out of `self.model.a_matrix_` (module doc).
    /// `okAssign(v, n, x)` -> `ok_assign(&mut self.v, n, x)`; `setBoundArrays` -> `set_bound_arrays()` (no-op, keep the call).
    /// `eqiters.assign(n, equations.end())` -> `self.eqiters = vec![None; n as usize]` (the try/catch is dropped).
    /// `Acol.insert(Acol.end(), collen, i)` -> `extend(std::iter::repeat(i).take(collen as usize))`.
    pub fn ok_from_csc(&mut self, aval: &[f64], aindex: &[HighsInt], astart: &[HighsInt]) -> bool {
        self.a_value.clear();
        self.a_col.clear();
        self.a_row.clear();

        self.freeslots.clear();
        if !ok_assign(&mut self.colhead, self.model.num_col_, -1) {
            return false;
        }
        if !ok_assign(&mut self.rowroot, self.model.num_row_, -1) {
            return false;
        }
        if !ok_assign(&mut self.colsize, self.model.num_col_, 0) {
            return false;
        }
        if !ok_assign(&mut self.rowsize, self.model.num_row_, 0) {
            return false;
        }
        if !ok_assign(&mut self.rowsize_integer, self.model.num_row_, 0) {
            return false;
        }
        if !ok_assign(&mut self.rowsize_impl_int, self.model.num_row_, 0) {
            return false;
        }
        if !ok_assign(&mut self.rowsize_singleton, self.model.num_row_, 0) {
            return false;
        }

        self.implied_row_bounds.set_num_sums(0);
        self.implied_dual_row_bounds.set_num_sums(0);
        self.implied_row_bounds.set_bound_arrays();
        self.implied_row_bounds.set_num_sums(self.model.num_row_);
        self.implied_dual_row_bounds.set_bound_arrays();
        self.implied_dual_row_bounds
            .set_num_sums(self.model.num_col_);

        let ncol = astart.len() as HighsInt - 1;
        debug_assert!(ncol as usize == self.colhead.len());
        let nnz = aval.len() as HighsInt;

        self.a_value = aval.to_vec();
        if !ok_reserve(&mut self.a_col, nnz) {
            return false;
        }
        if !ok_reserve(&mut self.a_row, nnz) {
            return false;
        }

        for i in 0..ncol {
            let collen = astart[(i + 1) as usize] - astart[i as usize];
            self.a_col.extend(std::iter::repeat_n(i, collen as usize));
            self.a_row.extend_from_slice(
                &aindex[astart[i as usize] as usize..astart[(i + 1) as usize] as usize],
            );
        }

        if !ok_resize(&mut self.a_next, nnz, 0) {
            return false;
        }
        if !ok_resize(&mut self.a_prev, nnz, 0) {
            return false;
        }
        if !ok_resize(&mut self.a_r_left, nnz, 0) {
            return false;
        }
        if !ok_resize(&mut self.a_r_right, nnz, 0) {
            return false;
        }
        for pos in 0..nnz {
            self.link(pos);
            // HiGHS #3318 (mip_coef_change_impl_bounds): HPresolve.cpp lines 1883-1887 at v1.16.0
            // (`link` no longer adds the nonzero; matrixNonZeroChanged would reset the
            // implied bounds preserved by shrinkProblem)
            if self.options.mip_coef_change_impl_bounds {
                let p = pos as usize;
                self.implied_row_bounds.add(
                    &col_bound_arrays!(self),
                    self.a_row[p],
                    self.a_col[p],
                    self.a_value[p],
                );
                self.implied_dual_row_bounds.add(
                    &row_dual_bound_arrays!(self),
                    self.a_col[p],
                    self.a_row[p],
                    self.a_value[p],
                );
            }
        }

        if self.equations.is_empty() {
            self.eqiters = vec![None; self.model.num_row_ as usize];
            for i in 0..self.model.num_row_ {
                // register equation
                if self.is_equation(i) {
                    let key = (self.rowsize[i as usize], i);
                    self.equations.insert(key);
                    self.eqiters[i as usize] = Some(key);
                }
            }
        }
        true
    }

    /// HiGHS: HPresolve.cpp::HPresolve::okFromCSR (lines 2622-2683, 62 lines)
    ///
    /// As `ok_from_csc`.
    pub fn ok_from_csr(
        &mut self,
        arval: &[f64],
        arindex: &[HighsInt],
        arstart: &[HighsInt],
    ) -> bool {
        self.a_value.clear();
        self.a_col.clear();
        self.a_row.clear();

        self.freeslots.clear();
        if !ok_assign(&mut self.colhead, self.model.num_col_, -1) {
            return false;
        }
        if !ok_assign(&mut self.rowroot, self.model.num_row_, -1) {
            return false;
        }
        if !ok_assign(&mut self.colsize, self.model.num_col_, 0) {
            return false;
        }
        if !ok_assign(&mut self.rowsize, self.model.num_row_, 0) {
            return false;
        }
        if !ok_assign(&mut self.rowsize_integer, self.model.num_row_, 0) {
            return false;
        }
        if !ok_assign(&mut self.rowsize_impl_int, self.model.num_row_, 0) {
            return false;
        }
        if !ok_assign(&mut self.rowsize_singleton, self.model.num_row_, 0) {
            return false;
        }

        self.implied_row_bounds.set_num_sums(0);
        self.implied_dual_row_bounds.set_num_sums(0);
        self.implied_row_bounds.set_bound_arrays();
        self.implied_row_bounds.set_num_sums(self.model.num_row_);
        self.implied_dual_row_bounds.set_bound_arrays();
        self.implied_dual_row_bounds
            .set_num_sums(self.model.num_col_);

        let nrow = arstart.len() as HighsInt - 1;
        debug_assert!(nrow as usize == self.rowroot.len());
        let nnz = arval.len() as HighsInt;

        self.a_value = arval.to_vec();
        if !ok_reserve(&mut self.a_col, nnz) {
            return false;
        }
        if !ok_reserve(&mut self.a_row, nnz) {
            return false;
        }
        //  entries.reserve(nnz);

        for i in 0..nrow {
            let rowlen = arstart[(i + 1) as usize] - arstart[i as usize];
            self.a_row.extend(std::iter::repeat_n(i, rowlen as usize));
            self.a_col.extend_from_slice(
                &arindex[arstart[i as usize] as usize..arstart[(i + 1) as usize] as usize],
            );
        }

        if !ok_resize(&mut self.a_next, nnz, 0) {
            return false;
        }
        if !ok_resize(&mut self.a_prev, nnz, 0) {
            return false;
        }
        if !ok_resize(&mut self.a_r_left, nnz, 0) {
            return false;
        }
        if !ok_resize(&mut self.a_r_right, nnz, 0) {
            return false;
        }
        for pos in 0..nnz {
            self.link(pos);
            // HiGHS #3318 (mip_coef_change_impl_bounds): HPresolve.cpp lines 1954-1958 at v1.16.0
            // (`link` no longer adds the nonzero; matrixNonZeroChanged would reset the
            // implied bounds preserved by shrinkProblem)
            if self.options.mip_coef_change_impl_bounds {
                let p = pos as usize;
                self.implied_row_bounds.add(
                    &col_bound_arrays!(self),
                    self.a_row[p],
                    self.a_col[p],
                    self.a_value[p],
                );
                self.implied_dual_row_bounds.add(
                    &row_dual_bound_arrays!(self),
                    self.a_col[p],
                    self.a_row[p],
                    self.a_value[p],
                );
            }
        }

        if self.equations.is_empty() {
            self.eqiters = vec![None; nrow as usize];
            for i in 0..nrow {
                // register equation
                if self.is_equation(i) {
                    let key = (self.rowsize[i as usize], i);
                    self.equations.insert(key);
                    self.eqiters[i as usize] = Some(key);
                }
            }
        }
        true
    }

    /// HiGHS: HPresolve.cpp::HPresolve::countFillin (lines 2685-2692, 8 lines)
    pub fn count_fillin(&mut self, row: HighsInt) -> HighsInt {
        let mut fillin: HighsInt = 0;
        for k in 0..self.rowpositions.len() {
            let rowiter = self.rowpositions[k];
            let col = self.a_col[rowiter as usize];
            if self.find_nonzero(row, col) == -1 {
                fillin += 1;
            }
        }

        fillin
    }

    /// HiGHS: HPresolve.cpp::HPresolve::checkFillin (lines 2694-2752, 59 lines)
    ///
    /// `fillinCache.find(k)` -> `fillin_cache.find(&k)` (`Option<&HighsInt>`);
    /// `HighsInt& x = fillinCache[k]` -> `fillin_cache.index_or_insert(&k)` (`&mut HighsInt`).
    pub fn check_fillin(
        &mut self,
        fillin_cache: &mut HighsHashTable<HighsInt, HighsInt>,
        row: HighsInt,
        col: HighsInt,
    ) -> bool {
        // check numerics against markowitz tolerance
        debug_assert!(self.rowpositions.len() == self.rowsize[row as usize] as usize);

        // check fillin against max fillin
        let mut fillin: HighsInt = -(self.rowsize[row as usize] + self.colsize[col as usize] - 1);

        // C++ `#if 1` branch (the `#else` branch is not compiled).
        // first use fillin for rows where it is already computed
        let mut coliter = self.colhead[col as usize];
        while coliter != -1 {
            let cur = coliter;
            coliter = self.a_next[cur as usize];
            if self.a_row[cur as usize] == row {
                continue;
            }

            let key = self.a_row[cur as usize];
            let cached_fillin = match fillin_cache.find(&key) {
                Some(&v) => v,
                None => continue,
            };

            fillin += cached_fillin - 1;
            if fillin > self.options.presolve_substitution_maxfillin {
                return false;
            }
        }

        // iterate over rows of substituted column again to count the fillin for the
        // remaining rows
        let mut coliter = self.colhead[col as usize];
        while coliter != -1 {
            let cur = coliter;
            coliter = self.a_next[cur as usize];
            debug_assert!(self.a_col[cur as usize] == col);

            if self.a_row[cur as usize] == row {
                continue;
            }

            let key = self.a_row[cur as usize];
            if *fillin_cache.index_or_insert(&key) != 0 {
                continue;
            }

            let rowfillin = self.count_fillin(key);
            *fillin_cache.index_or_insert(&key) = rowfillin + 1;
            fillin += rowfillin;

            if fillin > self.options.presolve_substitution_maxfillin {
                return false;
            }
            // we count a fillin of 1 if the column is not present in the row and
            // a fillin of zero otherwise. the fillin for the substituted column
            // itself was already counted before the loop so we skip that entry.
        }

        true
    }

    /// HiGHS: HPresolve.cpp::HPresolve::reinsertEquation (lines 2754-2763, 10 lines)
    pub fn reinsert_equation(&mut self, row: HighsInt) {
        let r = row as usize;
        // check if this is an equation row and it now has a different size
        if self.is_equation(row) {
            if let Some(old_key) = self.eqiters[r] {
                if old_key.0 != self.rowsize[r] {
                    // if that is the case reinsert it into the equation set that is ordered
                    // by sparsity
                    self.equations.remove(&old_key);
                    let new_key = (self.rowsize[r], row);
                    self.equations.insert(new_key);
                    self.eqiters[r] = Some(new_key);
                }
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::transformColumn (lines 2765-2855, 91 lines)
    pub fn transform_column(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
        scale: f64,
        constant: f64,
    ) -> HPresolveResult {
        // replace column x by x = scale * x' + constant
        if self.mipsolver.is_some() {
            self.mip_mut()
                .mipdata_
                .implications
                .column_transformed(col, scale, constant);
        }

        postsolve_stack.linear_transform(col, scale, constant);

        // new variable x' has the following bounds:
        // scale > 0 --> (lb - constant) / scale <= x' <= (ub - constant) / scale
        // scale < 0 --> (ub - constant) / scale <= x' <= (lb - constant) / scale
        // each matrix coefficient a will be replaced by scale * a and, therefore, the
        // contributions to the constraint activities will change as follows:
        // a * lb --> (a * scale) * (lb - constant) / scale = a * (lb - constant)
        // a * ub --> (a * scale) * (ub - constant) / scale = a * (ub - constant).
        // therefore, for now the scaling can be neglected and the bounds on
        // constraint activities can be updated using the constant term.
        let j = col as usize;
        if constant != 0.0 {
            let old_lower = self.model.col_lower_[j];
            let old_upper = self.model.col_upper_[j];
            self.model.col_upper_[j] -= constant;
            self.model.col_lower_[j] -= constant;

            let mut it = self.get_column_vector(col);
            while let Some(nonzero) = it.next_nz(self) {
                self.implied_row_bounds.updated_var_lower(
                    &col_bound_arrays!(self),
                    nonzero.index(),
                    col,
                    nonzero.value(),
                    old_lower,
                );
                self.implied_row_bounds.updated_var_upper(
                    &col_bound_arrays!(self),
                    nonzero.index(),
                    col,
                    nonzero.value(),
                    old_upper,
                );
            }

            let old_impl_lower = self.impl_col_lower[j];
            let old_impl_upper = self.impl_col_upper[j];
            self.impl_col_lower[j] -= constant;
            self.impl_col_upper[j] -= constant;

            let mut it = self.get_column_vector(col);
            while let Some(nonzero) = it.next_nz(self) {
                let lower_source = self.col_lower_source[j];
                self.implied_row_bounds.updated_impl_var_lower(
                    &col_bound_arrays!(self),
                    nonzero.index(),
                    col,
                    nonzero.value(),
                    old_impl_lower,
                    lower_source,
                );
                let upper_source = self.col_upper_source[j];
                self.implied_row_bounds.updated_impl_var_upper(
                    &col_bound_arrays!(self),
                    nonzero.index(),
                    col,
                    nonzero.value(),
                    old_impl_upper,
                    upper_source,
                );
            }
        }

        // now apply the scaling, which does not change the contributions to the
        // implied row bounds, but requires adjusting the implied bounds of the
        // columns dual constraint
        self.implied_dual_row_bounds.sum_scaled(col, scale);

        let bound_scale = 1.0 / scale;
        self.model.col_lower_[j] *= bound_scale;
        self.model.col_upper_[j] *= bound_scale;
        self.impl_col_lower[j] *= bound_scale;
        self.impl_col_upper[j] *= bound_scale;

        if scale < 0.0 {
            std::mem::swap(&mut self.model.col_lower_[j], &mut self.model.col_upper_[j]);
            std::mem::swap(&mut self.impl_col_lower[j], &mut self.impl_col_upper[j]);
            std::mem::swap(&mut self.col_lower_source[j], &mut self.col_upper_source[j]);
        }

        self.model.offset_ += self.model.col_cost_[j] * constant;
        self.model.col_cost_[j] *= scale;

        let mut coliter = self.colhead[j];
        while coliter != -1 {
            let ci = coliter as usize;
            let val = self.a_value[ci];
            self.a_value[ci] *= scale;
            let row = self.a_row[ci] as usize;
            let row_constant = val * constant;
            if self.model.row_lower_[row] != -kHighsInf {
                self.model.row_lower_[row] -= row_constant;
            }
            if self.model.row_upper_[row] != kHighsInf {
                self.model.row_upper_[row] -= row_constant;
            }
            coliter = self.a_next[ci];
        }

        // finally, use utility methods for rounding scaled bounds of
        // integer-constrained variables and updating bounds on constraint activities
        // accordingly. this should not be done before the preceding bound updates
        // (scaling and swaps) and matrix updates. we rely on the integrality status
        // being already updated to the newly scaled column by the caller, if
        // necessary.
        if self.model.integrality_[j] != HighsVarType::kContinuous {
            let new_lower = self.model.col_lower_[j];
            let new_upper = self.model.col_upper_[j];
            hpresolve_checked_call!(self.change_col_bounds(col, new_lower, new_upper));
        }

        self.mark_changed_col(col);
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::scaleRow (lines 2857-2861, 5 lines)
    ///
    /// C++ default arguments: `integral = false`.
    pub fn scale_row(&mut self, row: HighsInt, scale: f64, integral: bool) {
        self.store_row(row);

        self.scale_stored_row(row, scale, integral);
    }

    /// HiGHS: HPresolve.cpp::HPresolve::scaleStoredRow (lines 2863-2889, 27 lines)
    ///
    /// C++ default arguments: `integral = false`.
    pub fn scale_stored_row(&mut self, row: HighsInt, scale: f64, integral: bool) {
        let x1_clique_normalise = self.options.mip_clique_normalise;
        let r = row as usize;
        self.model.row_upper_[r] *= scale;
        self.model.row_lower_[r] *= scale;
        self.impl_row_dual_lower[r] /= scale;
        self.impl_row_dual_upper[r] /= scale;

        // HiGHS #3286 (mip_clique_normalise): HPresolve.cpp lines 3269-3270 at 843e9fdb6
        if x1_clique_normalise {
            self.row_dual_lower[r] /= 1.0_f64.copysign(scale);
            self.row_dual_upper[r] /= 1.0_f64.copysign(scale);
        } else if integral {
            if self.model.row_upper_[r] != kHighsInf {
                self.model.row_upper_[r] = self.model.row_upper_[r].round();
            }
            if self.model.row_lower_[r] != -kHighsInf {
                self.model.row_lower_[r] = self.model.row_lower_[r].round();
            }
        }

        let mut j = 0usize;
        while j < self.rowpositions.len() {
            let pos = self.rowpositions[j];
            let p = pos as usize;
            self.a_value[p] *= scale;
            if self.a_value[p].abs() <= self.options.small_matrix_value {
                self.unlink(pos);
            }
            j += 1;
        }

        self.implied_row_bounds.sum_scaled(row, scale);
        if scale < 0.0 {
            std::mem::swap(&mut self.row_dual_lower[r], &mut self.row_dual_upper[r]);
            std::mem::swap(
                &mut self.impl_row_dual_lower[r],
                &mut self.impl_row_dual_upper[r],
            );
            std::mem::swap(
                &mut self.row_dual_lower_source[r],
                &mut self.row_dual_upper_source[r],
            );
            std::mem::swap(&mut self.model.row_lower_[r], &mut self.model.row_upper_[r]);
        }

        // HiGHS #3286 (mip_clique_normalise): HPresolve.cpp lines 3285-3291 at 843e9fdb6
        if x1_clique_normalise && integral {
            if self.model.row_upper_[r] != kHighsInf {
                self.model.row_upper_[r] = (self.model.row_upper_[r] + self.primal_feastol).floor();
            }
            if self.model.row_lower_[r] != -kHighsInf {
                self.model.row_lower_[r] = (self.model.row_lower_[r] - self.primal_feastol).ceil();
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::substitute (lines 2891-2976, 86 lines)
    ///
    /// C++ overload `substitute(row, col, rhs)`.
    pub fn substitute(&mut self, row: HighsInt, col: HighsInt, rhs: f64) {
        debug_assert!(self.row_deleted[row as usize] == 0);
        debug_assert!(self.col_deleted[col as usize] == 0);
        let pos = self.find_nonzero(row, col);
        debug_assert!(pos != -1);

        debug_assert!(self.a_row[pos as usize] == row);
        debug_assert!(self.a_col[pos as usize] == col);
        let substrowscale = -1.0 / self.a_value[pos as usize];
        debug_assert!(self.is_implied_free(col));

        self.mark_row_deleted(row);
        self.mark_col_deleted(col);

        let h16_row_duals = self.options.h16_row_duals;
        // substitute the column in each row where it occurs
        let mut coliter = self.colhead[col as usize];
        while coliter != -1 {
            let colrow = self.a_row[coliter as usize];
            let colval = self.a_value[coliter as usize];

            // walk to the next position before doing any modifications, because
            // the current position will be deleted in the loop below
            debug_assert!(self.a_col[coliter as usize] == col);
            let colpos = coliter;
            coliter = self.a_next[coliter as usize];

            // skip the row that is used for substitution
            if row == colrow {
                continue;
            }

            debug_assert!(self.find_nonzero(colrow, col) != -1);

            // cancels out and bounds of dual row for this column do not need to be
            // updated
            self.unlink(colpos);

            // determine the scale for the substitution row for addition to this row
            let scale = colval * substrowscale;

            // adjust the sides
            let cr = colrow as usize;
            if h16_row_duals {
                // HiGHS 1.16 #3270 (h16_row_duals): HPresolveAggregator.cpp lines 204-205 (1.16.0)
                self.h16_add_to_row_lower(colrow, HighsCDouble::from(scale) * rhs);
                self.h16_add_to_row_upper(colrow, HighsCDouble::from(scale) * rhs);
            } else {
                if self.model.row_lower_[cr] != -kHighsInf {
                    self.model.row_lower_[cr] += scale * rhs;
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 3455 at f2ae2b320 (1.15.1 2933)
                    self.x1_row_modified(colrow);
                }

                if self.model.row_upper_[cr] != kHighsInf {
                    self.model.row_upper_[cr] += scale * rhs;
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 3456 at f2ae2b320 (1.15.1 2936)
                    self.x1_row_modified(colrow);
                }
            }

            for k in 0..self.rowpositions.len() {
                let rowiter = self.rowpositions[k];
                debug_assert!(self.a_row[rowiter as usize] == row);

                if self.a_col[rowiter as usize] != col {
                    let add_col = self.a_col[rowiter as usize];
                    let add_val = scale * self.a_value[rowiter as usize];
                    self.add_to_matrix(colrow, add_col, add_val);
                }
            }

            // check if this is an equation row and it now has a different size
            self.reinsert_equation(colrow);
        }

        debug_assert!(self.colsize[col as usize] == 1);

        // substitute column in the objective function
        if self.model.col_cost_[col as usize] != 0.0 {
            let objscale = HighsCDouble::from(self.model.col_cost_[col as usize] * substrowscale);
            self.model.offset_ = f64::from(self.model.offset_ - objscale * rhs);
            debug_assert!(self.model.offset_.is_finite());
            for k in 0..self.rowpositions.len() {
                let rowiter = self.rowpositions[k];
                let cc = self.a_col[rowiter as usize] as usize;
                let add = objscale * self.a_value[rowiter as usize];
                self.model.col_cost_[cc] = f64::from(self.model.col_cost_[cc] + add);
                if self.model.col_cost_[cc].abs() <= self.options.small_matrix_value {
                    self.model.col_cost_[cc] = 0.0;
                }
            }
            debug_assert!(
                self.model.col_cost_[col as usize].abs()
                    <= cpp_max(
                        self.options.dual_feasibility_tolerance,
                        kHighsTiny * f64::from(objscale).abs()
                    )
            );
            self.model.col_cost_[col as usize] = 0.0;
        }

        // finally remove the entries of the row that was used for substitution
        for k in 0..self.rowpositions.len() {
            let rowiter = self.rowpositions[k];
            self.unlink(rowiter);
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::toCSC (lines 2978-3006, 29 lines)
    ///
    /// Callers pass vectors taken out of `self.model.a_matrix_` (module doc).
    pub fn to_csc(
        &mut self,
        aval: &mut Vec<f64>,
        aindex: &mut Vec<HighsInt>,
        astart: &mut Vec<HighsInt>,
    ) {
        // set up the column starts using the column size array
        let numcol = self.colsize.len();
        astart.resize(numcol + 1, 0);
        let mut nnz: HighsInt = 0;
        for i in 0..numcol {
            astart[i] = nnz;
            nnz += self.colsize[i];
        }
        astart[numcol] = nnz;

        // now setup the entries of the CSC matrix
        // we reuse the colsize array to count down to zero
        // for determining the position of each nonzero
        aval.resize(nnz as usize, 0.0);
        aindex.resize(nnz as usize, 0);
        let numslots = self.a_value.len();
        debug_assert!(numslots - self.freeslots.len() == nnz as usize);
        for i in 0..numslots {
            if self.a_value[i] == 0.0 {
                continue;
            }
            let cc = self.a_col[i] as usize;
            debug_assert!(self.a_col[i] >= 0 && self.a_col[i] < self.model.num_col_);
            let pos = astart[cc + 1] - self.colsize[cc];
            self.colsize[cc] -= 1;
            debug_assert!(self.colsize[cc] >= 0);
            aval[pos as usize] = self.a_value[i];
            aindex[pos as usize] = self.a_row[i];
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::toCSR (lines 3008-3034, 27 lines)
    pub fn to_csr(
        &mut self,
        arval: &mut Vec<f64>,
        arindex: &mut Vec<HighsInt>,
        arstart: &mut Vec<HighsInt>,
    ) {
        // set up the row starts using the row size array
        let numrow = self.rowsize.len();
        arstart.resize(numrow + 1, 0);
        let mut nnz: HighsInt = 0;
        for i in 0..numrow {
            arstart[i] = nnz;
            nnz += self.rowsize[i];
        }
        arstart[numrow] = nnz;

        // now setup the entries of the CSC matrix
        // we reuse the colsize array to count down to zero
        // for determining the position of each nonzero
        arval.resize(nnz as usize, 0.0);
        arindex.resize(nnz as usize, 0);
        for i in 0..nnz {
            let iu = i as usize;
            if self.a_value[iu] == 0.0 {
                continue;
            }
            let r = self.a_row[iu] as usize;
            let pos = arstart[r + 1] - self.rowsize[r];
            self.rowsize[r] -= 1;
            debug_assert!(self.rowsize[r] >= 0);
            arval[pos as usize] = self.a_value[iu];
            arindex[pos as usize] = self.a_col[iu];
        }
    }
}
