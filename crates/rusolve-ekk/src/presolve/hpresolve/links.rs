//! HiGHS source: presolve/HPresolve.cpp (lines 452-825)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Nonzero linking, change tracking and implied-bound propagation. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::link (lines 452-475, 24 lines)
    ///
    /// `highs_splay_link(pos, rowroot[Arow[pos]], get_row_left, get_row_right, get_row_key)` ->
    /// `highs_splay_link(pos, &mut self.rowroot[r], &mut self.a_r_left, &mut self.a_r_right, &self.a_col)`
    /// with `r = self.a_row[pos as usize] as usize` (disjoint fields).
    /// `impliedRowBounds.add(...)` -> `self.implied_row_bounds.add(&col_bound_arrays!(self), ...)`.
    pub fn link(&mut self, pos: HighsInt) {
        let pi = pos as usize;
        let ci = self.a_col[pi] as usize;
        let ri = self.a_row[pi] as usize;

        self.a_next[pi] = self.colhead[ci];
        self.a_prev[pi] = -1;
        self.colhead[ci] = pos;
        if self.a_next[pi] != -1 {
            let next = self.a_next[pi] as usize;
            self.a_prev[next] = pos;
        }

        self.colsize[ci] += 1;

        self.a_r_left[pi] = -1;
        self.a_r_right[pi] = -1;
        highs_splay_link(
            pos,
            &mut self.rowroot[ri],
            &mut self.a_r_left,
            &mut self.a_r_right,
            &self.a_col,
        );

        // HiGHS #3318 (mip_coef_change_impl_bounds): HPresolve.cpp link (lines 536-562 at
        // v1.16.0) no longer adds the nonzero to the implied (dual) row bounds; its callers do.
        if !self.options.mip_coef_change_impl_bounds {
            self.implied_row_bounds.add(
                &col_bound_arrays!(self),
                self.a_row[pi],
                self.a_col[pi],
                self.a_value[pi],
            );
            self.implied_dual_row_bounds.add(
                &row_dual_bound_arrays!(self),
                self.a_col[pi],
                self.a_row[pi],
                self.a_value[pi],
            );
        }
        self.rowsize[ri] += 1;
        if self.model.integrality_[ci] == HighsVarType::kInteger {
            self.rowsize_integer[ri] += 1;
        } else if self.model.integrality_[ci] == HighsVarType::kImplicitInteger {
            self.rowsize_impl_int[ri] += 1;
        }

        // HiGHS 1.16 #3367: HPresolve.cpp lines 558-561 (1.16.0)
        if self.colsize[ci] == 1 {
            self.rowsize_singleton[ri] += 1;
        } else if self.colsize[ci] == 2 {
            let next_row = self.a_row[self.a_next[pi] as usize] as usize;
            self.rowsize_singleton[next_row] -= 1;
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::unlink (lines 477-531, 55 lines)
    ///
    /// Splay call as in `link` (`highs_splay_unlink`).
    pub fn unlink(&mut self, pos: HighsInt) {
        let pi = pos as usize;
        let ci = self.a_col[pi] as usize;
        let ri = self.a_row[pi] as usize;
        let x1_coef_change = self.options.mip_coef_change_impl_bounds;

        let next = self.a_next[pi];
        let prev = self.a_prev[pi];

        if next != -1 {
            self.a_prev[next as usize] = prev;
        }

        if prev != -1 {
            self.a_next[prev as usize] = next;
        } else {
            self.colhead[ci] = next;
        }
        self.colsize[ci] -= 1;

        if self.col_deleted[ci] == 0 {
            if self.colsize[ci] == 1 {
                self.singleton_columns.push(self.a_col[pi]);
            } else {
                self.mark_changed_col(self.a_col[pi]);
            }

            if !x1_coef_change {
                self.implied_dual_row_bounds.remove(
                    &row_dual_bound_arrays!(self),
                    self.a_col[pi],
                    self.a_row[pi],
                    self.a_value[pi],
                );
            }
        }

        highs_splay_unlink(
            pos,
            &mut self.rowroot[ri],
            &mut self.a_r_left,
            &mut self.a_r_right,
            &self.a_col,
        );
        self.rowsize[ri] -= 1;
        if self.model.integrality_[ci] == HighsVarType::kInteger {
            self.rowsize_integer[ri] -= 1;
        } else if self.model.integrality_[ci] == HighsVarType::kImplicitInteger {
            self.rowsize_impl_int[ri] -= 1;
        }

        // HiGHS 1.16 #3367: HPresolve.cpp lines 594-597 (1.16.0)
        if self.colsize[ci] == 0 {
            self.rowsize_singleton[ri] -= 1;
        } else if self.colsize[ci] == 1 {
            let head_row = self.a_row[self.colhead[ci] as usize] as usize;
            self.rowsize_singleton[head_row] += 1;
        }

        if self.row_deleted[ri] == 0 {
            if self.rowsize[ri] == 1 {
                self.singleton_rows.push(self.a_row[pi]);
            } else {
                self.mark_changed_row(self.a_row[pi]);
            }
            if !x1_coef_change {
                self.implied_row_bounds.remove(
                    &col_bound_arrays!(self),
                    self.a_row[pi],
                    self.a_col[pi],
                    self.a_value[pi],
                );
            }
        }

        if x1_coef_change {
            // HiGHS #3318 (mip_coef_change_impl_bounds): HPresolve.cpp lines 606-607 at v1.16.0
            self.matrix_non_zero_changed(
                self.a_row[pi],
                self.a_col[pi],
                self.a_value[pi],
                0.0,
                self.row_deleted[ri] != 0,
                self.col_deleted[ci] != 0,
            );
        } else {
            // remove implied bounds on row duals that where implied by this column's dual
            // constraint
            self.reset_row_dual_implied_bounds_derived_from_col(self.a_col[pi]);

            // remove implied bounds on columns that where implied by this row
            self.reset_col_implied_bounds_derived_from_row(self.a_row[pi]);

            // modifications to row invalidate lifting opportunities
            self.clear_lifting_opportunities(self.a_row[pi]);

            // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 639-640, 961-962 at f2ae2b320
            self.x1_row_modified(self.a_row[pi]);
        }

        // remove non-zero
        self.a_value[pi] = 0.0;

        self.freeslots.push(pos);
    }

    /// HiGHS: HPresolve.cpp::HPresolve::markChangedRow (lines 533-539, 7 lines)
    pub fn mark_changed_row(&mut self, row: HighsInt) {
        let ri = row as usize;
        if self.changed_row_flag[ri] == 0 {
            self.changed_row_indices.push(row);
            self.changed_row_flag[ri] = 1;
        }
        self.single_equation_checked[ri] = 0;
    }

    /// HiGHS: HPresolve.cpp::HPresolve::markChangedCol (lines 541-546, 6 lines)
    pub fn mark_changed_col(&mut self, col: HighsInt) {
        let ci = col as usize;
        if self.changed_col_flag[ci] == 0 {
            self.changed_col_indices.push(col);
            self.changed_col_flag[ci] = 1;
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::getMaxAbsColVal (lines 548-555, 8 lines)
    pub fn get_max_abs_col_val(&self, col: HighsInt) -> f64 {
        let mut max_val = 0.0;

        for nz in self.get_column_vector(col).iter(self) {
            max_val = cpp_max(nz.value().abs(), max_val);
        }

        max_val
    }

    /// HiGHS: HPresolve.cpp::HPresolve::getMaxAbsRowVal (lines 557-564, 8 lines)
    pub fn get_max_abs_row_val(&self, row: HighsInt) -> f64 {
        let mut max_val = 0.0;

        for nz in self.get_row_vector(row).iter(self) {
            max_val = cpp_max(nz.value().abs(), max_val);
        }

        max_val
    }

    /// HiGHS: HPresolve.cpp::HPresolve::checkUpdateRowDualImpliedBounds (lines 566-595, 30 lines)
    ///
    /// C++ default arguments: `dualRowLower = nullptr`, `dualRowUpper = nullptr`.
    /// `if (p != nullptr) *p = x` -> `if let Some(p) = p { *p = x }` (use `as_deref_mut()` to test twice).
    pub fn check_update_row_dual_implied_bounds(
        &self,
        col: HighsInt,
        dual_row_lower: Option<&mut f64>,
        dual_row_upper: Option<&mut f64>,
    ) -> bool {
        // check if implied bounds of row duals in given column can be updated (i.e.
        // dual row has finite bounds and number of infinite contributions to
        // corresponding activity bounds is at most one)

        // if the column has an infinite lower bound the reduced cost cannot be
        // positive, i.e. the column corresponds to a <= constraint in the dual with
        // right hand side -cost which becomes a >= constraint with side +cost.
        // Furthermore, we can ignore strictly redundant primal
        // column bounds and treat them as if they are infinite
        let ci = col as usize;
        let implied_margin = if self.colsize[ci] != 1 {
            self.primal_feastol
        } else {
            -self.primal_feastol
        };

        let my_dual_row_lower = if self.is_lower_strictly_implied(col, Some(&implied_margin)) {
            self.model.col_cost_[ci]
        } else {
            -kHighsInf
        };

        let my_dual_row_upper = if self.is_upper_strictly_implied(col, Some(&implied_margin)) {
            self.model.col_cost_[ci]
        } else {
            kHighsInf
        };

        if let Some(p) = dual_row_lower {
            *p = my_dual_row_lower;
        }
        if let Some(p) = dual_row_upper {
            *p = my_dual_row_upper;
        }

        (my_dual_row_lower != -kHighsInf
            && self.implied_dual_row_bounds.get_num_inf_sum_upper_orig(col) <= 1)
            || (my_dual_row_upper != kHighsInf
                && self.implied_dual_row_bounds.get_num_inf_sum_lower_orig(col) <= 1)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::updateRowDualImpliedBounds (lines 597-638, 42 lines)
    ///
    /// C++ overload `updateRowDualImpliedBounds(row, col, val)`.
    pub fn update_row_dual_implied_bounds(&mut self, row: HighsInt, col: HighsInt, val: f64) {
        /// `checkImplied` lambda (HPresolve.cpp lines 603-621). `this` is the
        /// lambda's captured `this`; `threshold` is its captured `threshold`.
        #[allow(clippy::too_many_arguments)]
        fn check_implied_bound_row_dual(
            this: &mut HPresolve,
            row: HighsInt,
            col: HighsInt,
            val: f64,
            dual_row_bnd: f64,
            residual_act: f64,
            direction: HighsInt,
            threshold: f64,
        ) {
            let ri = row as usize;
            if f64::from(direction) * residual_act <= -kHighsInf {
                return;
            }
            let implied_bound = f64::from((HighsCDouble::from(dual_row_bnd) - residual_act) / val);

            if implied_bound.abs() * kHighsTiny > this.options.dual_feasibility_tolerance {
                return;
            }

            if f64::from(direction) * val > 0.0 {
                // only tighten bound if it is tighter by a wide enough margin
                if implied_bound < this.impl_row_dual_upper[ri] - threshold {
                    this.change_impl_row_dual_upper(row, implied_bound, col);
                }
            } else if implied_bound > this.impl_row_dual_lower[ri] + threshold {
                this.change_impl_row_dual_lower(row, implied_bound, col);
            }
        }

        // propagate implied row dual bound
        let mut dual_row_lower = 0.0;
        let mut dual_row_upper = 0.0;
        if !self.check_update_row_dual_implied_bounds(
            col,
            Some(&mut dual_row_lower),
            Some(&mut dual_row_upper),
        ) {
            return;
        }

        let threshold = 1000.0 * self.options.dual_feasibility_tolerance;

        if dual_row_upper != kHighsInf {
            let residual_act = self.implied_dual_row_bounds.get_residual_sum_lower_orig(
                &row_dual_bound_arrays!(self),
                col,
                row,
                val,
                -1,
                -kHighsInf,
                -kHighsInf,
            );
            check_implied_bound_row_dual(
                self,
                row,
                col,
                val,
                dual_row_upper,
                residual_act,
                1,
                threshold,
            );
        }

        if dual_row_lower != -kHighsInf {
            let residual_act = self.implied_dual_row_bounds.get_residual_sum_upper_orig(
                &row_dual_bound_arrays!(self),
                col,
                row,
                val,
                -1,
                -kHighsInf,
                -kHighsInf,
            );
            check_implied_bound_row_dual(
                self,
                row,
                col,
                val,
                dual_row_lower,
                residual_act,
                -1,
                threshold,
            );
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::updateRowDualImpliedBounds (lines 640-646, 7 lines)
    ///
    /// C++ overload `updateRowDualImpliedBounds(col)`.
    pub fn update_row_dual_implied_bounds_col(&mut self, col: HighsInt) {
        // update dual implied bounds of all rows in given column
        debug_assert!(col >= 0 && col < self.model.num_col_);
        if !self.check_update_row_dual_implied_bounds(col, None, None) {
            return;
        }
        let mut it = self.get_column_vector(col);
        while let Some(nonzero) = it.next_nz(self) {
            self.update_row_dual_implied_bounds(nonzero.index(), col, nonzero.value());
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::checkUpdateColImpliedBounds (lines 648-667, 20 lines)
    ///
    /// C++ default arguments: `rowLower = nullptr`, `rowUpper = nullptr`.
    pub fn check_update_col_implied_bounds(
        &self,
        row: HighsInt,
        row_lower: Option<&mut f64>,
        row_upper: Option<&mut f64>,
    ) -> bool {
        // check if implied bounds of columns in given row can be updated (i.e. if
        // row's left-hand or right-hand side is finite and number of infinite
        // contributions to corresponding activity bounds is at most one)
        let ri = row as usize;
        let my_row_lower = if self.is_implied_equation_at_upper(row) {
            self.model.row_upper_[ri]
        } else {
            self.model.row_lower_[ri]
        };
        let my_row_upper = if self.is_implied_equation_at_lower(row) {
            self.model.row_lower_[ri]
        } else {
            self.model.row_upper_[ri]
        };
        debug_assert!(my_row_lower != kHighsInf);
        debug_assert!(my_row_upper != -kHighsInf);

        if let Some(p) = row_lower {
            *p = my_row_lower;
        }
        if let Some(p) = row_upper {
            *p = my_row_upper;
        }

        (my_row_lower != -kHighsInf && self.implied_row_bounds.get_num_inf_sum_upper_orig(row) <= 1)
            || (my_row_upper != kHighsInf
                && self.implied_row_bounds.get_num_inf_sum_lower_orig(row) <= 1)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::updateColImpliedBounds (lines 669-742, 74 lines)
    ///
    /// C++ overload `updateColImpliedBounds(row, col, val)`.
    ///
    /// PORT NOTE: signature changed because the C++ `useImplBound` test reads
    /// `mipsolver->mipdata_->postSolveStack`, which is the `postsolve_stack` HPresolve
    /// runs on (HighsMipSolverData.cpp:819 passes it). In the port that stack is moved
    /// out of `mipdata_` while presolve runs (`run_mip_presolve_data`), so the live stack
    /// is passed down instead.
    pub fn update_col_implied_bounds(
        &mut self,
        postsolve_stack: &HighsPostsolveStack,
        row: HighsInt,
        col: HighsInt,
        val: f64,
    ) -> HPresolveResult {
        /// `checkImpliedBound` lambda (HPresolve.cpp lines 680-717). `this` is the
        /// lambda's captured `this`; `threshold` is its captured `threshold`.
        #[allow(clippy::too_many_arguments)]
        fn check_implied_bound_col(
            this: &mut HPresolve,
            postsolve_stack: &HighsPostsolveStack,
            row: HighsInt,
            col: HighsInt,
            val: f64,
            row_bnd: f64,
            residual_act: f64,
            direction: HighsInt,
            threshold: f64,
        ) -> HPresolveResult {
            let ci = col as usize;
            if f64::from(direction) * residual_act <= -kHighsInf {
                return HPresolveResult::kOk;
            }
            let implied_bound = f64::from((HighsCDouble::from(row_bnd) - residual_act) / val);

            if implied_bound.abs() * kHighsTiny > this.primal_feastol {
                return HPresolveResult::kOk;
            }

            // do not use the implied bound if this a not a model row, since the
            // row can be removed and should not be used, e.g., to identify a
            // column as implied free
            // HiGHS: `useImplBound` (HPresolve.cpp 691-693)
            // HiGHS #3218 (mip_presolve_row_append): `!isCutRow(row)` (1.16.0 HPresolve.cpp
            // lines 778-779), so rows appended by FM presolve count as model rows.
            let use_impl_bound = this.mipsolver.is_none()
                || if this.options.mip_presolve_row_append {
                    !postsolve_stack.is_cut_row(row)
                } else {
                    postsolve_stack.get_orig_row_index(row) < this.mip().orig_model_.num_row_
                };

            if f64::from(direction) * val > 0.0 {
                // upper bound
                // when solving a MIP, keep (1) tighter bounds on integer variables or (2)
                // implied bound (as column bound)
                // HPresolve.cpp 699-703
                let update_col_bound = this.mipsolver.is_some()
                    && ((this.model.integrality_[ci] != HighsVarType::kContinuous
                        && implied_bound < this.model.col_upper_[ci] - this.primal_feastol)
                        || (!use_impl_bound
                            && implied_bound < this.model.col_upper_[ci] - threshold));

                // only tighten bound if it is tighter by a wide enough margin
                if use_impl_bound && implied_bound < this.impl_col_upper[ci] - threshold {
                    this.change_impl_col_upper(col, implied_bound, row);
                }
                if update_col_bound {
                    hpresolve_checked_call!(this.change_col_upper(col, implied_bound));
                }
            } else {
                // lower bound
                // when solving a MIP, keep (1) tighter bounds on integer variables or (2)
                // implied bound (as column bound)
                // HPresolve.cpp 715-719
                let update_col_bound = this.mipsolver.is_some()
                    && ((this.model.integrality_[ci] != HighsVarType::kContinuous
                        && implied_bound > this.model.col_lower_[ci] + this.primal_feastol)
                        || (!use_impl_bound
                            && implied_bound > this.model.col_lower_[ci] + threshold));

                // only tighten bound if it is tighter by a wide enough margin
                if use_impl_bound && implied_bound > this.impl_col_lower[ci] + threshold {
                    this.change_impl_col_lower(col, implied_bound, row);
                }
                if update_col_bound {
                    hpresolve_checked_call!(this.change_col_lower(col, implied_bound));
                }
            }
            HPresolveResult::kOk
        }

        // propagate implied column bound upper bound if row has an upper bound
        let mut row_lower = 0.0;
        let mut row_upper = 0.0;
        if !self.check_update_col_implied_bounds(row, Some(&mut row_lower), Some(&mut row_upper)) {
            return HPresolveResult::kOk;
        }

        let threshold = 1000.0 * self.primal_feastol;

        if row_upper != kHighsInf {
            let residual_act = self.implied_row_bounds.get_residual_sum_lower_orig(
                &col_bound_arrays!(self),
                row,
                col,
                val,
                -1,
                -kHighsInf,
                -kHighsInf,
            );
            hpresolve_checked_call!(check_implied_bound_col(
                self,
                postsolve_stack,
                row,
                col,
                val,
                row_upper,
                residual_act,
                1,
                threshold,
            ));
        }

        if row_lower != -kHighsInf {
            let residual_act = self.implied_row_bounds.get_residual_sum_upper_orig(
                &col_bound_arrays!(self),
                row,
                col,
                val,
                -1,
                -kHighsInf,
                -kHighsInf,
            );
            hpresolve_checked_call!(check_implied_bound_col(
                self,
                postsolve_stack,
                row,
                col,
                val,
                row_lower,
                residual_act,
                -1,
                threshold,
            ));
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::updateColImpliedBounds (lines 744-752, 9 lines)
    ///
    /// C++ overload `updateColImpliedBounds(row)`.
    ///
    /// PORT NOTE: takes `postsolve_stack` for `update_col_implied_bounds` (see there).
    pub fn update_col_implied_bounds_row(
        &mut self,
        postsolve_stack: &HighsPostsolveStack,
        row: HighsInt,
    ) -> HPresolveResult {
        // update implied bounds of all columns in given row
        debug_assert!(row >= 0 && row < self.model.num_row_);
        if !self.check_update_col_implied_bounds(row, None, None) {
            return HPresolveResult::kOk;
        }
        let mut it = self.get_row_vector(row);
        while let Some(nonzero) = it.next_nz(self) {
            hpresolve_checked_call!(self.update_col_implied_bounds(
                postsolve_stack,
                row,
                nonzero.index(),
                nonzero.value()
            ));
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::resetColImpliedBounds (lines 754-767, 14 lines)
    ///
    /// C++ default arguments: `row = -1`.
    pub fn reset_col_implied_bounds(&mut self, col: HighsInt, row: HighsInt) {
        let ci = col as usize;
        debug_assert!(
            row == -1 || self.col_lower_source[ci] == row || self.col_upper_source[ci] == row
        );
        if self.col_deleted[ci] == 0 {
            // set implied bounds to infinite values if (1) they were deduced from the
            // given row or (2) no row was given
            if self.col_lower_source[ci] != -1 && (row == -1 || self.col_lower_source[ci] == row) {
                self.change_impl_col_lower(col, -kHighsInf, -1);
            }
            if self.col_upper_source[ci] != -1 && (row == -1 || self.col_upper_source[ci] == row) {
                self.change_impl_col_upper(col, kHighsInf, -1);
            }
        } else if row != -1 && self.row_deleted[row as usize] == 0 {
            // remove column from row-wise implied bound storage
            self.col_impl_source_by_row[row as usize].remove(&col);
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::resetRowDualImpliedBounds (lines 769-785, 17 lines)
    ///
    /// C++ default arguments: `col = -1`.
    pub fn reset_row_dual_implied_bounds(&mut self, row: HighsInt, col: HighsInt) {
        let ri = row as usize;
        debug_assert!(
            col == -1
                || self.row_dual_lower_source[ri] == col
                || self.row_dual_upper_source[ri] == col
        );
        if self.row_deleted[ri] == 0 {
            // set implied bounds to infinite values if (1) they were deduced from the
            // given column or (2) no column was given
            if self.row_dual_lower_source[ri] != -1
                && (col == -1 || self.row_dual_lower_source[ri] == col)
            {
                self.change_impl_row_dual_lower(row, -kHighsInf, -1);
            }
            if self.row_dual_upper_source[ri] != -1
                && (col == -1 || self.row_dual_upper_source[ri] == col)
            {
                self.change_impl_row_dual_upper(row, kHighsInf, -1);
            }
        } else if col != -1 && self.col_deleted[col as usize] == 0 {
            // remove row from column-wise implied bound storage
            self.impl_row_dual_source_by_col[col as usize].remove(&row);
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::resetColImpliedBoundsDerivedFromRow (lines 787-798, 12 lines)
    ///
    /// `std::set<HighsInt> affectedCols` -> `BTreeSet<HighsInt>`; iterate in set order.
    pub fn reset_col_implied_bounds_derived_from_row(&mut self, row: HighsInt) {
        // reset implied column bounds affected by a modification in a row
        // (removed / added non-zeros, etc.)
        let ri = row as usize;
        if self.col_impl_source_by_row[ri].is_empty() {
            return;
        }
        let affected_cols = std::mem::take(&mut self.col_impl_source_by_row[ri]);
        for col in affected_cols {
            // set implied bounds to infinite values if they were deduced from the
            // given row
            self.reset_col_implied_bounds(col, row);
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::resetRowDualImpliedBoundsDerivedFromCol (lines 800-811, 12 lines)
    pub fn reset_row_dual_implied_bounds_derived_from_col(&mut self, col: HighsInt) {
        // reset implied row dual bounds affected by a modification in a column
        // (removed / added non-zeros, etc.)
        let ci = col as usize;
        if self.impl_row_dual_source_by_col[ci].is_empty() {
            return;
        }
        let affected_rows = std::mem::take(&mut self.impl_row_dual_source_by_col[ci]);
        for row in affected_rows {
            // set implied bounds to infinite values if they were deduced from the
            // given column
            self.reset_row_dual_implied_bounds(row, col);
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::findNonzero (lines 813-825, 13 lines)
    ///
    /// Non-`const` in C++: the splay restructures the row tree.
    /// `self.rowroot[row] = highs_splay(col, self.rowroot[row], &mut self.a_r_left, &mut self.a_r_right, &self.a_col)`.
    pub fn find_nonzero(&mut self, row: HighsInt, col: HighsInt) -> HighsInt {
        let ri = row as usize;
        if self.rowroot[ri] == -1 {
            return -1;
        }

        let root = self.rowroot[ri];
        self.rowroot[ri] = highs_splay(
            col,
            root,
            &mut self.a_r_left,
            &mut self.a_r_right,
            &self.a_col,
        );

        if self.a_col[self.rowroot[ri] as usize] == col {
            return self.rowroot[ri];
        }

        -1
    }
}
