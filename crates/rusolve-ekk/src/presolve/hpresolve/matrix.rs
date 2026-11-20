//! HiGHS source: presolve/HPresolve.cpp (lines 2069-2550)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Matrix modification and bound changes. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;
use crate::util::highs_matrix_slice::SliceStack;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::addToMatrix (lines 2069-2132, 64 lines)
    ///
    /// `Avalue.push_back(...)` etc. grow the triplet arrays; `freeslots.back()`/`pop_back()`.
    pub fn add_to_matrix(&mut self, row: HighsInt, col: HighsInt, val: f64) {
        let mut pos = self.find_nonzero(row, col);

        self.mark_changed_row(row);
        self.mark_changed_col(col);

        if pos == -1 {
            if self.freeslots.is_empty() {
                pos = self.a_value.len() as HighsInt;
                self.a_value.push(val);
                self.a_row.push(row);
                self.a_col.push(col);
                self.a_next.push(-1);
                self.a_prev.push(-1);
                self.a_r_left.push(-1);
                self.a_r_right.push(-1);
            } else {
                pos = *self.freeslots.last().unwrap();
                self.freeslots.pop();
                self.a_value[pos as usize] = val;
                self.a_row[pos as usize] = row;
                self.a_col[pos as usize] = col;
                self.a_prev[pos as usize] = -1;
            }

            self.link(pos);

            // HiGHS #3318 (mip_coef_change_impl_bounds): HPresolve.cpp line 1305 at v1.16.0
            if self.options.mip_coef_change_impl_bounds {
                self.matrix_non_zero_changed(row, col, 0.0, val, false, false);
                return;
            }

            // remove implied bounds on row duals that where implied by this column's
            // dual constraint
            self.reset_row_dual_implied_bounds_derived_from_col(col);

            // remove implied bounds on columns that where implied by this row
            self.reset_col_implied_bounds_derived_from_row(row);

            // modifications to row invalidate lifting opportunities
            self.clear_lifting_opportunities(row);

            // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 961-962 at f2ae2b320
            self.x1_row_modified(row);
        } else {
            let sum = self.a_value[pos as usize] + val;
            if sum.abs() <= self.options.small_matrix_value {
                self.unlink(pos);
            } else if self.options.mip_coef_change_impl_bounds {
                // HiGHS #3318 (mip_coef_change_impl_bounds): HPresolve.cpp lines 1312-1313 at v1.16.0
                self.matrix_non_zero_changed(
                    row,
                    col,
                    self.a_value[pos as usize],
                    sum,
                    false,
                    false,
                );
                self.a_value[pos as usize] = sum;
            } else {
                // remove implied bounds on row duals that where implied by this column's
                // dual constraint
                self.reset_row_dual_implied_bounds_derived_from_col(col);

                // remove implied bounds on columns that where implied by this row
                self.reset_col_implied_bounds_derived_from_row(row);

                // modifications to row invalidate lifting opportunities
                self.clear_lifting_opportunities(row);

                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 961-962 at f2ae2b320
                self.x1_row_modified(row);

                // remove the locks and contribution to implied (dual) row bounds, then
                // add then again
                self.implied_row_bounds.remove(
                    &col_bound_arrays!(self),
                    row,
                    col,
                    self.a_value[pos as usize],
                );
                self.implied_dual_row_bounds.remove(
                    &row_dual_bound_arrays!(self),
                    col,
                    row,
                    self.a_value[pos as usize],
                );
                self.a_value[pos as usize] = sum;
                // value not zero, add new contributions and locks with opposite sign
                self.implied_row_bounds.add(
                    &col_bound_arrays!(self),
                    row,
                    col,
                    self.a_value[pos as usize],
                );
                self.implied_dual_row_bounds.add(
                    &row_dual_bound_arrays!(self),
                    col,
                    row,
                    self.a_value[pos as usize],
                );
            }
        }
    }

    /// HiGHS latest (#3264, f2ae2b320): HPresolve.cpp::HPresolve::matrixNonZeroChanged (lines 961-962), changeRowLower (line 990), changeRowUpper (line 1019)
    ///
    /// Rust-only gate (mip_aggregate_implications): the C++ `if (mipsolver != nullptr)
    /// mipsolver->mipdata_->implications.rowModified(row)`, called at the 1.15.1 row writes.
    pub fn x1_row_modified(&mut self, row: HighsInt) {
        if !self.options.mip_aggregate_implications {
            return;
        }
        if let Some(ms) = self.mipsolver.as_mut() {
            ms.mipdata_.implications.row_modified(row);
        }
    }

    /// HiGHS 1.16 (#3218): HPresolve.cpp::HPresolve::addToMatrix(postsolve_stack, row_lower,
    /// row_upper, row_indices, row_values) (1.16.0 lines 1321-1407): append rows to the
    /// presolve model. Reached only from Fourier-Motzkin presolve (`mip_fm_presolve`,
    /// which requires `mip_presolve_row_append`).
    ///
    /// `impliedDualRowBounds.setBoundArrays` is a no-op in the port (no raw pointers).
    /// `mipsolver->numCol()/numRow()` read the presolved model, which the port keeps in
    /// `self.model` while HPresolve runs.
    pub fn add_to_matrix_rows(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row_lower: &[f64],
        row_upper: &[f64],
        row_indices: &[Vec<HighsInt>],
        row_values: &[Vec<f64>],
    ) -> bool {
        // update number of rows
        let num_rows = row_indices.len() as HighsInt;
        if num_rows == 0 {
            return true;
        }
        let old_num_rows = self.model.num_row_;
        self.model.num_row_ += num_rows;
        self.model.a_matrix_.num_row_ += num_rows;
        let nr = self.model.num_row_;

        // resize postsolve vectors
        postsolve_stack.append_rows_to_model(num_rows);

        // add row bounds
        self.model.row_lower_.extend_from_slice(row_lower);
        self.model.row_upper_.extend_from_slice(row_upper);

        // initialise row sizes
        if !ok_resize(&mut self.rowroot, nr, -1) {
            return false;
        }
        if !ok_resize(&mut self.rowsize, nr, 0) {
            return false;
        }
        if !ok_resize(&mut self.rowsize_integer, nr, 0) {
            return false;
        }
        if !ok_resize(&mut self.rowsize_impl_int, nr, 0) {
            return false;
        }
        if !ok_resize(&mut self.rowsize_singleton, nr, 0) {
            return false;
        }

        // initialise row duals
        if !ok_resize(&mut self.row_dual_lower, nr, -kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_upper, nr, kHighsInf) {
            return false;
        }
        for i in old_num_rows as usize..nr as usize {
            if self.model.row_lower_[i] == -kHighsInf {
                self.row_dual_upper[i] = 0.0;
            }
            if self.model.row_upper_[i] == kHighsInf {
                self.row_dual_lower[i] = 0.0;
            }
        }

        // initialise implied row duals
        if !ok_resize(&mut self.impl_row_dual_lower, nr, -kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.impl_row_dual_upper, nr, kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_lower_source, nr, -1) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_upper_source, nr, -1) {
            return false;
        }
        if !ok_resize(&mut self.col_impl_source_by_row, nr, BTreeSet::new()) {
            return false;
        }

        // initialise flags
        if !ok_resize(&mut self.changed_row_flag, nr, 0) {
            return false;
        }
        if !ok_resize(&mut self.row_deleted, nr, 0) {
            return false;
        }
        if !ok_resize(&mut self.single_equation_checked, nr, 0) {
            return false;
        }

        // initialise row names
        if !ok_resize(&mut self.model.row_names_, nr, String::new()) {
            return false;
        }

        // resize vector for equations
        if !ok_resize(&mut self.eqiters, nr, None) {
            return false;
        }

        // resize vectors for implied row bounds
        self.implied_row_bounds.set_num_sums(nr);

        // set bound arrays again (no-op in the port)
        self.implied_dual_row_bounds.set_bound_arrays();

        let ncols = self.model.num_col_;
        if let Some(ms) = self.mipsolver.as_mut() {
            let dual_fix = ms.options_mip_.mip_dual_fix_probing;
            ms.mipdata_.implications.resize_x1b(ncols, nr, dual_fix);
        }

        for i in 0..num_rows as usize {
            // new row index
            let row = old_num_rows + i as HighsInt;

            // add non-zeros
            for j in 0..row_indices[i].len() {
                self.add_to_matrix(row, row_indices[i][j], row_values[i][j]);
            }

            // add row singleton
            if self.rowsize[row as usize] == 1 {
                self.singleton_rows.push(row);
            }

            // add equation
            if self.is_equation(row) {
                let key = (self.rowsize[row as usize], row);
                self.equations.insert(key);
                self.eqiters[row as usize] = Some(key);
            }
        }

        true
    }

    /// HiGHS 1.16 (#3218): HPresolve.cpp::HPresolve::addToMatrix single-row overload
    /// (1.16.0 lines 1409-1417)
    pub fn add_to_matrix_row(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row_lower: f64,
        row_upper: f64,
        row_indices: &[HighsInt],
        row_values: &[f64],
    ) -> bool {
        self.add_to_matrix_rows(
            postsolve_stack,
            &[row_lower],
            &[row_upper],
            &[row_indices.to_vec()],
            &[row_values.to_vec()],
        )
    }

    /// HiGHS: HPresolve.cpp::HPresolve::markRowDeleted (lines 2149-2162, 14 lines)
    pub fn mark_row_deleted(&mut self, row: HighsInt) {
        debug_assert!(self.row_deleted[row as usize] == 0);

        // remove equations from set of equations
        if self.is_equation(row) {
            if let Some(key) = self.eqiters[row as usize].take() {
                self.equations.remove(&key);
            }
        }

        // prevents row from being added to change vector
        self.changed_row_flag[row as usize] = 1;
        self.row_deleted[row as usize] = 1;
        self.num_deleted_rows += 1;
    }

    /// HiGHS: HPresolve.cpp::HPresolve::markColDeleted (lines 2164-2171, 8 lines)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::markColDeleted (lines 2679-2691)
    pub fn mark_col_deleted(&mut self, col: HighsInt) {
        debug_assert!(self.col_deleted[col as usize] == 0);

        // prevents col from being added to change vector
        self.changed_col_flag[col as usize] = 1;
        self.col_deleted[col as usize] = 1;
        self.num_deleted_cols += 1;
        // HiGHS 1.16 (#3255): HPresolve.cpp line 1462; -1 unless FM reformulated.
        if col == self.model.fme_obj_col_ {
            self.model.fme_obj_col_ = -1;
        }
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 2688-2690 at 7b890c255
        if self.options.mip_presolve_clique_table {
            self.x1b_eliminate_clique_col(col);
        }
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp lines 2688-2690 (also 4210-4212, 9407-9409):
    /// `if (mipsolver != nullptr && mipsolver->mipdata_->cliquesExtracted)
    /// presolveCliqueTable.eliminateCol(col);`. Called only from gated sites.
    pub fn x1b_eliminate_clique_col(&mut self, col: HighsInt) {
        if let Some(ms) = self.mipsolver.as_mut() {
            if ms.mipdata_.cliques_extracted {
                self.presolve_clique_table
                    .eliminate_col(&mut ms.mipdata_.cliquetable, col);
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeColUpper (lines 2173-2190, 18 lines)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::changeColUpper (lines 2693-2715)
    pub fn change_col_upper(&mut self, col: HighsInt, new_upper: f64) -> HPresolveResult {
        let mut new_upper = new_upper;
        if self.model.integrality_[col as usize] != HighsVarType::kContinuous {
            new_upper = (new_upper + self.primal_feastol).floor();
            if new_upper == self.model.col_upper_[col as usize] {
                return HPresolveResult::kOk;
            }
        }

        let old_upper = self.model.col_upper_[col as usize];
        self.model.col_upper_[col as usize] = new_upper;

        // HiGHS #3332 (mip_dual_side_lost): HPresolve.cpp lines 1477-1479 at v1.16.0
        self.x1_reset_row_dual_implied_bounds_if_dual_side_lost(col, false);

        hpresolve_checked_call!(self.check_col_bounds(col, None));

        let mut it = self.get_column_vector(col);
        while let Some(nonzero) = it.next_nz(self) {
            self.implied_row_bounds.updated_var_upper(
                &col_bound_arrays!(self),
                nonzero.index(),
                col,
                nonzero.value(),
                old_upper,
            );
            self.mark_changed_row(nonzero.index());
        }
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 2709-2713 at 7b890c255
        if self.options.mip_presolve_clique_table
            && self.mipsolver.is_some()
            && self.model.col_lower_[col as usize] == self.model.col_upper_[col as usize]
        {
            hpresolve_checked_call!(
                self.update_clique_table_fixed_col(col, self.model.col_lower_[col as usize])
            );
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeColLower (lines 2192-2209, 18 lines)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::changeColLower (lines 2717-2739)
    pub fn change_col_lower(&mut self, col: HighsInt, new_lower: f64) -> HPresolveResult {
        let mut new_lower = new_lower;
        if self.model.integrality_[col as usize] != HighsVarType::kContinuous {
            new_lower = (new_lower - self.primal_feastol).ceil();
            if new_lower == self.model.col_lower_[col as usize] {
                return HPresolveResult::kOk;
            }
        }

        let old_lower = self.model.col_lower_[col as usize];
        self.model.col_lower_[col as usize] = new_lower;

        // HiGHS #3332 (mip_dual_side_lost): HPresolve.cpp lines 1505-1507 at v1.16.0
        self.x1_reset_row_dual_implied_bounds_if_dual_side_lost(col, true);

        hpresolve_checked_call!(self.check_col_bounds(col, None));

        let mut it = self.get_column_vector(col);
        while let Some(nonzero) = it.next_nz(self) {
            self.implied_row_bounds.updated_var_lower(
                &col_bound_arrays!(self),
                nonzero.index(),
                col,
                nonzero.value(),
                old_lower,
            );
            self.mark_changed_row(nonzero.index());
        }
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 2733-2737 at 7b890c255
        if self.options.mip_presolve_clique_table
            && self.mipsolver.is_some()
            && self.model.col_lower_[col as usize] == self.model.col_upper_[col as usize]
        {
            hpresolve_checked_call!(
                self.update_clique_table_fixed_col(col, self.model.col_lower_[col as usize])
            );
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeColBounds (lines 2211-2223, 13 lines)
    pub fn change_col_bounds(
        &mut self,
        col: HighsInt,
        new_lower: f64,
        new_upper: f64,
    ) -> HPresolveResult {
        if new_lower > self.model.col_upper_[col as usize] {
            // change upper bound first
            hpresolve_checked_call!(self.change_col_upper(col, new_upper));
            hpresolve_checked_call!(self.change_col_lower(col, new_lower));
        } else {
            // change lower bound first
            hpresolve_checked_call!(self.change_col_lower(col, new_lower));
            hpresolve_checked_call!(self.change_col_upper(col, new_upper));
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::checkColBounds (lines 2225-2241, 17 lines)
    ///
    /// C++ default arguments: `isFixed = nullptr` (`None`).
    pub fn check_col_bounds(
        &mut self,
        col: HighsInt,
        is_fixed: Option<&mut bool>,
    ) -> HPresolveResult {
        debug_assert!(self.col_deleted[col as usize] == 0);
        let bound_diff = self.model.col_upper_[col as usize] - self.model.col_lower_[col as usize];
        let mut is_fixed = is_fixed;
        if let Some(f) = is_fixed.as_deref_mut() {
            *f = false;
        }
        if bound_diff <= self.primal_feastol
            && (bound_diff <= self.options.small_matrix_value
                || self.get_max_abs_col_val(col) * bound_diff <= self.primal_feastol)
        {
            // check for primal infeasibility
            if bound_diff < -self.primal_feastol {
                return HPresolveResult::kPrimalInfeasible;
            }
            // check for unboundedness
            if self.model.col_lower_[col as usize].abs() == kHighsInf {
                return HPresolveResult::kDualInfeasible;
            }
            // column is fixed
            if let Some(f) = is_fixed {
                *f = true;
            }
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeRowDualUpper (lines 2243-2255, 13 lines)
    pub fn change_row_dual_upper(&mut self, row: HighsInt, new_upper: f64) {
        let old_upper = self.row_dual_upper[row as usize];
        self.row_dual_upper[row as usize] = new_upper;

        let mut it = self.get_row_vector(row);
        while let Some(nonzero) = it.next_nz(self) {
            self.implied_dual_row_bounds.updated_var_upper(
                &row_dual_bound_arrays!(self),
                nonzero.index(),
                row,
                nonzero.value(),
                old_upper,
            );
            self.mark_changed_col(nonzero.index());
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeRowDualLower (lines 2257-2269, 13 lines)
    pub fn change_row_dual_lower(&mut self, row: HighsInt, new_lower: f64) {
        let old_lower = self.row_dual_lower[row as usize];
        self.row_dual_lower[row as usize] = new_lower;

        let mut it = self.get_row_vector(row);
        while let Some(nonzero) = it.next_nz(self) {
            self.implied_dual_row_bounds.updated_var_lower(
                &row_dual_bound_arrays!(self),
                nonzero.index(),
                row,
                nonzero.value(),
                old_lower,
            );
            self.mark_changed_col(nonzero.index());
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeImplColUpper (lines 2271-2313, 43 lines)
    pub fn change_impl_col_upper(&mut self, col: HighsInt, new_upper: f64, origin_row: HighsInt) {
        let cl = col as usize;
        let old_impl_upper = self.impl_col_upper[cl];
        let old_upper_source = self.col_upper_source[cl];
        if old_impl_upper >= self.model.col_upper_[cl] - self.primal_feastol
            && new_upper < self.model.col_upper_[cl] - self.primal_feastol
        {
            // the dual constraint can be considered a >= constraint and was free, or a
            // <= constraint before
            self.mark_changed_col(col);
        }
        let new_implied_free = self.is_lower_implied(col)
            && old_impl_upper > self.model.col_upper_[cl] + self.primal_feastol
            && new_upper <= self.model.col_upper_[cl] + self.primal_feastol;

        // remember the source of this upper bound, so that we can correctly identify
        // weak domination
        if old_upper_source != origin_row {
            if old_upper_source != -1 && old_upper_source != self.col_lower_source[cl] {
                self.col_impl_source_by_row[old_upper_source as usize].remove(&col);
            }
            if origin_row != -1 {
                self.col_impl_source_by_row[origin_row as usize].insert(col);
            }

            self.col_upper_source[cl] = origin_row;
        }

        // update implied bound
        self.impl_col_upper[cl] = new_upper;

        // HiGHS #3332 (mip_dual_side_lost): HPresolve.cpp lines 1612-1614 at v1.16.0
        self.x1_reset_row_dual_implied_bounds_if_dual_side_lost(col, false);

        // if the old and the new implied bound are not better than the upper bound,
        // nothing needs to be updated
        if !new_implied_free && cpp_min(old_impl_upper, new_upper) >= self.model.col_upper_[cl] {
            return;
        }

        let mut it = self.get_column_vector(col);
        while let Some(nonzero) = it.next_nz(self) {
            self.implied_row_bounds.updated_impl_var_upper(
                &col_bound_arrays!(self),
                nonzero.index(),
                col,
                nonzero.value(),
                old_impl_upper,
                old_upper_source,
            );
            if new_implied_free && self.is_dual_implied_free(nonzero.index()) {
                self.substitution_opportunities.push((nonzero.index(), col));
            }

            self.mark_changed_row(nonzero.index());
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeImplColLower (lines 2315-2357, 43 lines)
    pub fn change_impl_col_lower(&mut self, col: HighsInt, new_lower: f64, origin_row: HighsInt) {
        let cl = col as usize;
        let old_impl_lower = self.impl_col_lower[cl];
        let old_lower_source = self.col_lower_source[cl];
        if old_impl_lower <= self.model.col_lower_[cl] + self.primal_feastol
            && new_lower > self.model.col_lower_[cl] + self.primal_feastol
        {
            // the dual constraint can additionally be considered a <= constraint and
            // was free, or a >= constraint before
            self.mark_changed_col(col);
        }
        let new_implied_free = self.is_upper_implied(col)
            && old_impl_lower < self.model.col_lower_[cl] - self.primal_feastol
            && new_lower >= self.model.col_lower_[cl] - self.primal_feastol;

        // remember the source of this lower bound, so that we can correctly identify
        // weak domination
        if old_lower_source != origin_row {
            if old_lower_source != -1 && old_lower_source != self.col_upper_source[cl] {
                self.col_impl_source_by_row[old_lower_source as usize].remove(&col);
            }
            if origin_row != -1 {
                self.col_impl_source_by_row[origin_row as usize].insert(col);
            }

            self.col_lower_source[cl] = origin_row;
        }

        // update implied bound
        self.impl_col_lower[cl] = new_lower;

        // HiGHS #3332 (mip_dual_side_lost): HPresolve.cpp lines 1660-1662 at v1.16.0
        self.x1_reset_row_dual_implied_bounds_if_dual_side_lost(col, true);

        // if the old and the new implied bound are not better than the lower bound,
        // nothing needs to be updated
        if !new_implied_free && cpp_max(old_impl_lower, new_lower) <= self.model.col_lower_[cl] {
            return;
        }

        let mut it = self.get_column_vector(col);
        while let Some(nonzero) = it.next_nz(self) {
            self.implied_row_bounds.updated_impl_var_lower(
                &col_bound_arrays!(self),
                nonzero.index(),
                col,
                nonzero.value(),
                old_impl_lower,
                old_lower_source,
            );
            if new_implied_free && self.is_dual_implied_free(nonzero.index()) {
                self.substitution_opportunities.push((nonzero.index(), col));
            }

            self.mark_changed_row(nonzero.index());
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeImplRowDualUpper (lines 2359-2398, 40 lines)
    /// HiGHS latest (#3437, 6d361ef85): HPresolve.cpp::HPresolve::changeImplRowDualUpper (lines 1634-1680)
    pub fn change_impl_row_dual_upper(
        &mut self,
        row: HighsInt,
        new_upper: f64,
        origin_col: HighsInt,
    ) {
        let r = row as usize;
        let ranged_dual_implied = self.options.mip_ranged_dual_implied;
        let old_impl_upper = self.impl_row_dual_upper[r];
        let old_upper_source = self.row_dual_upper_source[r];

        if old_impl_upper >= -self.options.dual_feasibility_tolerance
            && new_upper < -self.options.dual_feasibility_tolerance
        {
            self.mark_changed_row(row);
        }

        // HiGHS #3437 (mip_ranged_dual_implied): HPresolve.cpp lines 1655-1656 at 6d361ef85
        let new_dual_implied = if ranged_dual_implied {
            !self.is_dual_implied_free(row)
                && self.model.row_upper_[r] != kHighsInf
                && new_upper <= self.options.dual_feasibility_tolerance
        } else {
            !self.is_dual_implied_free(row)
                && old_impl_upper > self.row_dual_upper[r] + self.options.dual_feasibility_tolerance
                && new_upper <= self.row_dual_upper[r] + self.options.dual_feasibility_tolerance
        };

        // remember the source of this upper bound, so that we can correctly identify
        // weak domination
        if old_upper_source != origin_col {
            if old_upper_source != -1 && old_upper_source != self.row_dual_lower_source[r] {
                self.impl_row_dual_source_by_col[old_upper_source as usize].remove(&row);
            }
            if origin_col != -1 {
                self.impl_row_dual_source_by_col[origin_col as usize].insert(row);
            }

            self.row_dual_upper_source[r] = origin_col;
        }

        // update implied bound
        self.impl_row_dual_upper[r] = new_upper;

        // nothing needs to be updated
        if !new_dual_implied && cpp_min(old_impl_upper, new_upper) >= self.row_dual_upper[r] {
            return;
        }

        let mut it = self.get_row_vector(row);
        while let Some(nonzero) = it.next_nz(self) {
            self.implied_dual_row_bounds.updated_impl_var_upper(
                &row_dual_bound_arrays!(self),
                nonzero.index(),
                row,
                nonzero.value(),
                old_impl_upper,
                old_upper_source,
            );
            self.mark_changed_col(nonzero.index());

            if new_dual_implied && self.is_implied_free(nonzero.index()) {
                self.substitution_opportunities.push((row, nonzero.index()));
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::changeImplRowDualLower (lines 2400-2439, 40 lines)
    /// HiGHS latest (#3437, 6d361ef85): HPresolve.cpp::HPresolve::changeImplRowDualLower (lines 1682-1728)
    pub fn change_impl_row_dual_lower(
        &mut self,
        row: HighsInt,
        new_lower: f64,
        origin_col: HighsInt,
    ) {
        let r = row as usize;
        let ranged_dual_implied = self.options.mip_ranged_dual_implied;
        let old_impl_lower = self.impl_row_dual_lower[r];
        let old_lower_source = self.row_dual_lower_source[r];

        if old_impl_lower <= self.options.dual_feasibility_tolerance
            && new_lower > self.options.dual_feasibility_tolerance
        {
            self.mark_changed_row(row);
        }

        // HiGHS #3437 (mip_ranged_dual_implied): HPresolve.cpp lines 1702-1703 at 6d361ef85
        let new_dual_implied = if ranged_dual_implied {
            !self.is_dual_implied_free(row)
                && self.model.row_lower_[r] != -kHighsInf
                && new_lower >= -self.options.dual_feasibility_tolerance
        } else {
            !self.is_dual_implied_free(row)
                && old_impl_lower < self.row_dual_lower[r] - self.options.dual_feasibility_tolerance
                && new_lower >= self.row_dual_lower[r] - self.options.dual_feasibility_tolerance
        };

        // remember the source of this lower bound, so that we can correctly identify
        // weak domination
        if old_lower_source != origin_col {
            if old_lower_source != -1 && old_lower_source != self.row_dual_upper_source[r] {
                self.impl_row_dual_source_by_col[old_lower_source as usize].remove(&row);
            }
            if origin_col != -1 {
                self.impl_row_dual_source_by_col[origin_col as usize].insert(row);
            }

            self.row_dual_lower_source[r] = origin_col;
        }

        // update implied bound
        self.impl_row_dual_lower[r] = new_lower;

        // nothing needs to be updated
        if !new_dual_implied && cpp_max(old_impl_lower, new_lower) <= self.row_dual_lower[r] {
            return;
        }

        let mut it = self.get_row_vector(row);
        while let Some(nonzero) = it.next_nz(self) {
            self.implied_dual_row_bounds.updated_impl_var_lower(
                &row_dual_bound_arrays!(self),
                nonzero.index(),
                row,
                nonzero.value(),
                old_impl_lower,
                old_lower_source,
            );
            self.mark_changed_col(nonzero.index());

            if new_dual_implied && self.is_implied_free(nonzero.index()) {
                self.substitution_opportunities.push((row, nonzero.index()));
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::getRowPositions (lines 2543-2550, 8 lines)
    ///
    /// `for (auto rowiter = rowVector.begin(); ...) myrowpositions.push_back(rowiter.position())` ->
    /// iterate `self.get_row_vector(row).iter(self)` and push `nz.position()`.
    pub fn get_row_positions(&self, row: HighsInt, myrowpositions: &mut Vec<HighsInt>) {
        myrowpositions.clear();

        // PORT NOTE (perf): the in-order walk of
        // HighsMatrixSlice<HighsTripletTreeSliceInOrder>::iterator (begin():
        // push -1 then descend left; operator++: right child then descend left,
        // else pop) written out on local state, so the cursor lives in
        // registers and only the two link arrays are read. Same push/pop
        // sequence, so the same positions in the same order.
        let left = &self.a_r_left[..];
        let right = &self.a_r_right[..];
        let mut cur = self.rowroot[row as usize];
        if cur == -1 {
            return;
        }
        // reporting-only. One unit per nonzero of the row walked.
        crate::work::add_presolve_units(self.rowsize[row as usize] as u64);
        let mut stack = SliceStack::new();
        stack.push(-1);
        loop {
            // descend to the leftmost node of the subtree rooted at `cur`
            let mut l = left[cur as usize];
            while l != -1 {
                stack.push(cur);
                cur = l;
                l = left[cur as usize];
            }
            // visit, then climb until a node with a right subtree is visited
            loop {
                myrowpositions.push(cur);
                let r = right[cur as usize];
                if r != -1 {
                    cur = r;
                    break;
                }
                cur = stack.pop().unwrap_or(-1);
                if cur == -1 {
                    return;
                }
            }
        }
    }
}

#[cfg(test)]
mod x1_rdif_tests {
    use super::*;
    use crate::prelude::*;
    use crate::presolve::highs_postsolve_stack::HighsPostsolveStack;
    use crate::util::highs_timer::HighsTimer;

    fn ranged_row_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![-1.0, 0.0];
        lp.col_lower_ = vec![-kHighsInf, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![2.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        lp
    }

    /// Latest HPresolveTest.cpp::presolveRuleTestAggregator (#3437, diff 202-230; the
    /// entry point is not ported to `col_stuffing.rs`, so the test drives the same
    /// steps): compute the implied row dual bounds of every column, then run the
    /// aggregator. Then shrink the problem, as `HPresolve::run` does after presolve.
    fn rule_test_aggregator(flag: bool) -> (HPresolve, HighsPostsolveStack, HighsOptions) {
        let lp = ranged_row_lp();
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_ranged_dual_implied = flag;
        let mut postsolve_stack = HighsPostsolveStack::default();
        postsolve_stack.initialize_index_maps(lp.num_row_, lp.num_col_);
        let mut timer = HighsTimer::default();
        timer.start();
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, Some(timer)));
        let (model, opts) = (presolve.model.clone(), presolve.options.clone());
        presolve.analysis_.setup(&model, &opts, true);
        for col in 0..presolve.model.num_col_ {
            if presolve.col_deleted[col as usize] == 0 {
                presolve.update_row_dual_implied_bounds_col(col);
            }
        }
        assert_eq!(
            presolve.aggregator(&mut postsolve_stack),
            HPresolveResult::kOk
        );
        presolve.shrink_problem(&mut postsolve_stack);
        (presolve, postsolve_stack, options)
    }

    /// HiGHS #3437: port of TestPresolveRules `test-aggregator-ranged-row`.
    #[test]
    fn x1_rdif_aggregator_ranged_row() {
        // flag off (1.15.1): the ranged row is not a substitution candidate
        let (presolve, _, _) = rule_test_aggregator(false);
        assert_eq!(presolve.model.num_col_, 2);
        assert_eq!(presolve.model.num_row_, 1);

        // flag on: x and the row have been removed
        let (presolve, mut postsolve_stack, options) = rule_test_aggregator(true);
        assert_eq!(presolve.model.num_col_, 1);
        assert_eq!(presolve.model.num_row_, 0);

        // the presolved problem is min -2 + z, z >= 0, with optimal solution z = 0
        // and reduced cost 1
        let mut solution = HighsSolution::default();
        solution.col_value = vec![0.0];
        solution.col_dual = vec![1.0];
        solution.value_valid = true;
        solution.dual_valid = true;
        let mut basis = HighsBasis::default();
        postsolve_stack.undo(&options, &mut solution, &mut basis, -1, false);

        // postsolved solution: x = 2, z = 0, row dual -1
        assert!((solution.col_value[0] - 2.0).abs() < 1e-8);
        assert!(solution.col_value[1].abs() < 1e-8);
        assert!((solution.row_dual[0] + 1.0).abs() < 1e-8);
    }

    /// Full presolve + solve with the flag off and on give the same optimum.
    #[test]
    fn x1_rdif_full_solve_ranged_row() {
        let mut solutions = Vec::new();
        for flag in [false, true] {
            let mut h = Highs::default();
            h.set_option_value_bool("output_flag", false);
            h.set_option_value_bool("mip_ranged_dual_implied", flag);
            assert_eq!(h.pass_model(ranged_row_lp()), HighsStatus::kOk);
            h.run();
            assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
            let sol = h.get_solution().clone();
            assert!((sol.col_value[0] - 2.0).abs() < 1e-8);
            assert!(sol.col_value[1].abs() < 1e-8);
            assert!((sol.row_dual[0] + 1.0).abs() < 1e-8);
            solutions.push(sol);
        }
        assert_eq!(solutions[0].col_value, solutions[1].col_value);
        assert_eq!(solutions[0].row_dual, solutions[1].row_dual);
    }
}

#[cfg(test)]
mod x1_agg_c_tests {
    use super::*;
    use crate::prelude::*;

    fn x1_agg_c_mip(flag: bool) -> Vec<f64> {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![-1.0, -1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![1.0, 10.0];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kInteger];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![0.0, 7.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![-5.0, 1.0, 1.0, 1.0];

        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_bool("mip_aggregate_implications", flag);
        h.pass_model(lp);
        h.run();
        assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
        h.get_solution().clone().col_value
    }

    #[test]
    fn x1_agg_c_off_solve() {
        let sol = x1_agg_c_mip(false);
        assert!((sol[0] - 1.0).abs() < 1e-6);
        assert!((sol[1] - 5.0).abs() < 1e-6);
    }

    #[test]
    fn x1_agg_c_on_matches_off() {
        let sol_off = x1_agg_c_mip(false);
        let sol_on = x1_agg_c_mip(true);
        assert_eq!(sol_off.len(), sol_on.len());
        for i in 0..sol_off.len() {
            assert!((sol_off[i] - sol_on[i]).abs() < 1e-6);
        }
    }

    /// Row 0: `x - 10 y <= 0`, row 1: `x + y <= 7`; `x` continuous in [0, 10], `y` binary.
    /// Sets the MIP solver on HPresolve with the option `flag`, extracts the variable bounds
    /// of row 0, and returns (origin of the VUB x <= 10 y, number of row-0 entries).
    fn vub_origin(presolve: &HPresolve) -> (HighsInt, usize) {
        let imp = &presolve.mip().mipdata_.implications;
        let origin = imp.get_vubs(0).find(&1).map_or(-2, |vb| vb.origin);
        let mut n = 0usize;
        imp.get_row_var_bounds(0).for_each(&mut |_c, _b| n += 1);
        (origin, n)
    }

    fn vub_presolve(flag: bool) -> HPresolve {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, -10.0, 1.0];
        lp.col_cost_ = vec![-1.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 1.0];
        lp.integrality_ = vec![HighsVarType::kContinuous, HighsVarType::kInteger];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![0.0, 7.0];

        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_aggregate_implications = flag;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        // the MIP solver data sizes the implication table before presolve
        let imp = crate::mip::implications::HighsImplications::new(&ms);
        *ms.mipdata_.implications = imp;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input_mip(Box::new(ms), options.presolve_reduction_limit));
        presolve.extract_var_bounds(0);
        presolve
    }

    /// Option on: the size-2 one-sided row is recorded as the VUB origin, and a coefficient
    /// change of that row (`add_to_matrix`, HPresolve.cpp 2523-2530 / 961-962 at f2ae2b320)
    /// clears it.
    #[test]
    fn x1_agg_c_on_records_and_clears_origin() {
        let mut presolve = vub_presolve(true);
        assert_eq!(vub_origin(&presolve), (0, 1));
        presolve.add_to_matrix(0, 1, 1.0);
        assert_eq!(vub_origin(&presolve), (-1, 0));
    }

    /// Option off (1.15.1): the VUB carries origin -1 and no row is recorded.
    #[test]
    fn x1_agg_c_off_records_no_origin() {
        let mut presolve = vub_presolve(false);
        assert_eq!(vub_origin(&presolve), (-1, 0));
        presolve.add_to_matrix(0, 1, 1.0);
        assert_eq!(vub_origin(&presolve), (-1, 0));
    }
}

#[cfg(test)]
mod x1b_psr_c_tests {
    use super::*;
    use crate::prelude::*;
    use crate::presolve::highs_postsolve_stack::HighsPostsolveStack;
    use crate::util::highs_timer::HighsTimer;

    /// HiGHS #3218: port of TestPresolve.cpp "add-to-matrix" (1.16.0 lines 1096-1155).
    #[test]
    fn x1b_psr_c_add_to_matrix() {
        // min x0 + x1 s.t. x0 + x1 >= 2, 0 <= x0,x1 <= 3
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![3.0, 3.0];
        lp.row_lower_ = vec![2.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];

        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_presolve_row_append = true;
        let mut postsolve_stack = HighsPostsolveStack::default();
        postsolve_stack.initialize_index_maps(lp.num_row_, lp.num_col_);
        let mut timer = HighsTimer::default();
        timer.start();
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, &options, 0, Some(timer)));

        // a single redundant row: x0 + x1 <= 10
        assert!(presolve.add_to_matrix_row(
            &mut postsolve_stack,
            -kHighsInf,
            10.0,
            &[0, 1],
            &[1.0, 1.0]
        ));
        assert_eq!(presolve.model.num_row_, 2);

        // multiple redundant rows: x0 <= 5, x1 <= 5
        assert!(presolve.add_to_matrix_rows(
            &mut postsolve_stack,
            &[-kHighsInf, -kHighsInf],
            &[5.0, 5.0],
            &[vec![0], vec![1]],
            &[vec![1.0], vec![1.0]],
        ));
        assert_eq!(presolve.model.num_row_, 4);
        assert_eq!(presolve.singleton_rows, vec![2, 3]);
        assert!(postsolve_stack.is_appended_row(1) && postsolve_stack.is_appended_row(3));
        assert_eq!(postsolve_stack.get_next_row_index(), 4);

        // rebuild a_matrix_ from the internal triplet structure
        presolve.shrink_problem(&mut postsolve_stack);
        let lp = presolve.model.clone();
        assert_eq!(lp.num_row_, 4);

        for presolve_opt in ["off", "on"] {
            let mut h = Highs::default();
            h.set_option_value_bool("output_flag", false);
            h.set_option_value_str("presolve", presolve_opt);
            assert_eq!(h.pass_model(lp.clone()), HighsStatus::kOk);
            h.run();
            assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
            assert_eq!(h.get_info().objective_function_value, 2.0);
        }
    }
}
