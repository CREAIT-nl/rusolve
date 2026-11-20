//! HiGHS source: presolve/HPresolve.cpp (lines 47-450)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Input setup and the implied-bound predicates. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::debugPrintRow (lines 47-67, 21 lines)
    ///
    /// C++ compiles this only `#ifndef NDEBUG`; its call sites are commented out.
    /// `printf` -> `print!`. The output is diagnostic only: any rendering of the same
    /// values is acceptable for `%.15g`.
    pub fn debug_print_row(&mut self, postsolve_stack: &mut HighsPostsolveStack, row: HighsInt) {
        print!(
            "(row {}) {} (impl: {}) <= ",
            postsolve_stack.get_orig_row_index(row),
            self.model.row_lower_[row as usize],
            self.implied_row_bounds.get_sum_lower(row, 0.0)
        );
        for nonzero in self.get_sorted_row_vector(row).iter(self) {
            let colchar =
                if self.model.integrality_[nonzero.index() as usize] == HighsVarType::kInteger {
                    'y'
                } else {
                    'x'
                };
            let signchar = if nonzero.value() < 0.0 { '-' } else { '+' };
            print!(
                "{}{} {}{} ",
                signchar,
                nonzero.value().abs(),
                colchar,
                postsolve_stack.get_orig_col_index(nonzero.index())
            );
        }
        println!(
            "<= {} (impl: {})",
            self.model.row_upper_[row as usize],
            self.implied_row_bounds.get_sum_upper(row, 0.0)
        );
    }

    /// HiGHS: HPresolve.cpp::HPresolve::okSetInput (lines 70-148, 79 lines)
    ///
    /// C++ default arguments: `timer = nullptr` (`None`).
    /// `model = &model_` -> `self.model = model_` (owned; see the module doc),
    /// `options = &options_` -> `self.options = options_.clone()`,
    /// `this->timer = timer` -> `self.timer = timer`. After that every `model_.x` /
    /// `model->x` is `self.model.x`.
    /// `okResize(v, n, x)` -> `ok_resize(&mut self.v, n, x)`, `okReserve(v, n)` ->
    /// `ok_reserve(&mut self.v, n)`, `okReserve(liftingOpportunities, n)` ->
    /// `self.lifting_opportunities.reserve(n as usize)` (all return `true`).
    /// `okFromCSC(model->a_matrix_.value_, ...)`: take the three vectors out
    /// (`std::mem::take`), call `self.ok_from_csc(&value, &index, &start)`, put them back.
    /// `kHighsSize_tInf` is `usize::MAX`; `presolve_reduction_limit as usize`.
    pub fn ok_set_input(
        &mut self,
        model_: HighsLp,
        options_: &HighsOptions,
        presolve_reduction_limit: HighsInt,
        timer: Option<HighsTimer>,
    ) -> bool {
        self.model = model_;
        self.options = options_.clone();
        self.timer = timer;

        // HiGHS 1.16 #3167 (h16_presolve_rules): HPresolve.cpp::setInput lines 65-97 (1.16.0): no
        // data structures here; presolve_h16.rs builds them inside presolve.
        let h16_presolve_rules = self.options.h16_presolve_rules;
        if h16_presolve_rules {
            return self.h16_set_input(presolve_reduction_limit);
        }

        if !ok_resize(&mut self.col_lower_source, self.model.num_col_, -1) {
            return false;
        }
        if !ok_resize(&mut self.col_upper_source, self.model.num_col_, -1) {
            return false;
        }
        if !ok_resize(&mut self.impl_col_lower, self.model.num_col_, -kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.impl_col_upper, self.model.num_col_, kHighsInf) {
            return false;
        }
        if !ok_resize(
            &mut self.col_impl_source_by_row,
            self.model.num_row_,
            BTreeSet::new(),
        ) {
            return false;
        }
        if !ok_resize(
            &mut self.impl_row_dual_source_by_col,
            self.model.num_col_,
            BTreeSet::new(),
        ) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_lower, self.model.num_row_, -kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_upper, self.model.num_row_, kHighsInf) {
            return false;
        }
        if !ok_resize(
            &mut self.impl_row_dual_lower,
            self.model.num_row_,
            -kHighsInf,
        ) {
            return false;
        }
        if !ok_resize(
            &mut self.impl_row_dual_upper,
            self.model.num_row_,
            kHighsInf,
        ) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_upper_source, self.model.num_row_, -1) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_lower_source, self.model.num_row_, -1) {
            return false;
        }

        for i in 0..self.model.num_row_ {
            if self.model.row_lower_[i as usize] == -kHighsInf {
                self.row_dual_upper[i as usize] = 0.0;
            }
            if self.model.row_upper_[i as usize] == kHighsInf {
                self.row_dual_lower[i as usize] = 0.0;
            }
        }

        if self.mipsolver.is_none() {
            self.primal_feastol = self.options.primal_feasibility_tolerance;
            ok_assign(
                &mut self.model.integrality_,
                self.model.num_col_,
                HighsVarType::kContinuous,
            );
        } else {
            self.primal_feastol = self.options.mip_feasibility_tolerance;
        }

        if self.model.a_matrix_.is_rowwise() {
            // Does this even happen?
            debug_assert!(self.model.a_matrix_.is_colwise());
            let value = std::mem::take(&mut self.model.a_matrix_.value_);
            let index = std::mem::take(&mut self.model.a_matrix_.index_);
            let start = std::mem::take(&mut self.model.a_matrix_.start_);
            let ok = self.ok_from_csr(&value, &index, &start);
            self.model.a_matrix_.value_ = value;
            self.model.a_matrix_.index_ = index;
            self.model.a_matrix_.start_ = start;
            if !ok {
                return false;
            }
        } else {
            let value = std::mem::take(&mut self.model.a_matrix_.value_);
            let index = std::mem::take(&mut self.model.a_matrix_.index_);
            let start = std::mem::take(&mut self.model.a_matrix_.start_);
            let ok = self.ok_from_csc(&value, &index, &start);
            self.model.a_matrix_.value_ = value;
            self.model.a_matrix_.index_ = index;
            self.model.a_matrix_.start_ = start;
            if !ok {
                return false;
            }
        }

        // initialize everything as changed, but do not add all indices
        // since the first thing presolve will do is a scan for easy reductions
        // of each row and column and set the flag of processed columns to false
        // from then on they are added to the vector whenever there are changes
        if !ok_resize(&mut self.changed_row_flag, self.model.num_row_, 1u8) {
            return false;
        }
        if !ok_resize(&mut self.row_deleted, self.model.num_row_, 0u8) {
            return false;
        }
        if !ok_reserve(&mut self.changed_row_indices, self.model.num_row_) {
            return false;
        }
        if !ok_resize(&mut self.changed_col_flag, self.model.num_col_, 1u8) {
            return false;
        }
        if !ok_resize(&mut self.col_deleted, self.model.num_col_, 0u8) {
            return false;
        }
        if !ok_reserve(&mut self.changed_col_indices, self.model.num_col_) {
            return false;
        }
        // C++ `okReserve(liftingOpportunities, num_row_)`: a `BTreeMap` has no capacity to
        // reserve, and the reserve is an allocation hint with no observable effect.
        if !ok_resize(&mut self.single_equation_checked, self.model.num_row_, 0u8) {
            return false;
        }
        self.num_deleted_cols = 0;
        self.num_deleted_rows = 0;
        // initialize substitution opportunities
        for row in 0..self.model.num_row_ {
            if !self.is_dual_implied_free(row) {
                continue;
            }
            let mut it = self.get_row_vector(row);
            while let Some(nonzero) = it.next_nz(self) {
                if self.is_implied_free(nonzero.index()) {
                    self.substitution_opportunities.push((row, nonzero.index()));
                }
            }
        }
        // Take value passed in as reduction limit, allowing different
        // values to be used for initial presolve, and after restart
        self.reduction_limit = if presolve_reduction_limit < 0 {
            kHighsSize_tInf
        } else {
            presolve_reduction_limit as usize
        };
        if self.options.presolve != kHighsOffString && self.reduction_limit < kHighsSize_tInf {
            crate::highs_log_dev!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "HPresolve::okSetInput reductionLimit = {}\n",
                self.reduction_limit as i32
            );
        }
        true
    }

    /// HiGHS: HPresolve.cpp::HPresolve::rowCoefficientsIntegral (lines 174-181, 8 lines)
    pub fn row_coefficients_integral(&self, row: HighsInt, scale: f64) -> bool {
        let fractionality = |x: f64| (x - x.round()).abs();
        for nz in self.get_row_vector(row).iter(self) {
            if fractionality(nz.value() * scale) > self.options.small_matrix_value {
                return false;
            }
        }

        true
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isLowerImplied (lines 183-186, 4 lines)
    pub fn is_lower_implied(&self, col: HighsInt) -> bool {
        self.model.col_lower_[col as usize] == -kHighsInf
            || self.impl_col_lower[col as usize]
                >= self.model.col_lower_[col as usize] - self.primal_feastol
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isLowerStrictlyImplied (lines 188-193, 6 lines)
    ///
    /// C++ default arguments: `tolerance = nullptr` (`None`).
    pub fn is_lower_strictly_implied(&self, col: HighsInt, tolerance: Option<&f64>) -> bool {
        let tol = match tolerance {
            Some(t) => *t,
            None => self.primal_feastol,
        };
        self.model.col_lower_[col as usize] == -kHighsInf
            || self.impl_col_lower[col as usize] > self.model.col_lower_[col as usize] + tol
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isUpperImplied (lines 195-198, 4 lines)
    pub fn is_upper_implied(&self, col: HighsInt) -> bool {
        self.model.col_upper_[col as usize] == kHighsInf
            || self.impl_col_upper[col as usize]
                <= self.model.col_upper_[col as usize] + self.primal_feastol
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isUpperStrictlyImplied (lines 200-205, 6 lines)
    ///
    /// C++ default arguments: `tolerance = nullptr` (`None`).
    pub fn is_upper_strictly_implied(&self, col: HighsInt, tolerance: Option<&f64>) -> bool {
        let tol = match tolerance {
            Some(t) => *t,
            None => self.primal_feastol,
        };
        self.model.col_upper_[col as usize] == kHighsInf
            || self.impl_col_upper[col as usize] < self.model.col_upper_[col as usize] - tol
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isImpliedFree (lines 207-209, 3 lines)
    pub fn is_implied_free(&self, col: HighsInt) -> bool {
        self.is_lower_implied(col) && self.is_upper_implied(col)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isDualImpliedFree (lines 211-217, 7 lines)
    pub fn is_dual_implied_free(&self, row: HighsInt) -> bool {
        self.is_equation(row)
            || (self.model.row_upper_[row as usize] != kHighsInf
                && self.impl_row_dual_upper[row as usize]
                    <= self.options.dual_feasibility_tolerance)
            || (self.model.row_lower_[row as usize] != -kHighsInf
                && self.impl_row_dual_lower[row as usize]
                    >= -self.options.dual_feasibility_tolerance)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::dualImpliedFreeGetRhsAndRowType (lines 219-236, 18 lines)
    ///
    /// C++ default arguments: `relaxRowDualBounds = false`.
    /// `HighsPostsolveStack::RowType::kEq` -> `RowType::kEq`.
    pub fn dual_implied_free_get_rhs_and_row_type(
        &mut self,
        row: HighsInt,
        rhs: &mut f64,
        row_type: &mut RowType,
        relax_row_dual_bounds: bool,
    ) {
        debug_assert!(self.is_dual_implied_free(row));
        if self.is_equation(row) {
            *row_type = RowType::kEq;
            *rhs = self.model.row_upper_[row as usize];
        } else if self.model.row_upper_[row as usize] != kHighsInf
            && self.impl_row_dual_upper[row as usize] <= self.options.dual_feasibility_tolerance
        {
            *row_type = RowType::kLeq;
            *rhs = self.model.row_upper_[row as usize];
            if relax_row_dual_bounds {
                self.change_row_dual_upper(row, kHighsInf);
            }
        } else {
            *row_type = RowType::kGeq;
            *rhs = self.model.row_lower_[row as usize];
            if relax_row_dual_bounds {
                self.change_row_dual_lower(row, -kHighsInf);
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isEquation (lines 238-240, 3 lines)
    pub fn is_equation(&self, row: HighsInt) -> bool {
        self.model.row_lower_[row as usize] == self.model.row_upper_[row as usize]
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isRanged (lines 242-245, 4 lines)
    pub fn is_ranged(&self, row: HighsInt) -> bool {
        self.model.row_lower_[row as usize] != -kHighsInf
            && self.model.row_upper_[row as usize] != kHighsInf
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isRedundant (lines 247-252, 6 lines)
    pub fn is_redundant(&self, row: HighsInt) -> bool {
        self.implied_row_bounds.get_sum_lower(row, 0.0)
            >= self.model.row_lower_[row as usize] - self.primal_feastol
            && self.implied_row_bounds.get_sum_upper(row, 0.0)
                <= self.model.row_upper_[row as usize] + self.primal_feastol
    }

    /// HiGHS: HPresolve.cpp::HPresolve::yieldsImpliedLowerBound (lines 254-257, 4 lines)
    pub fn yields_implied_lower_bound(&self, row: HighsInt, val: f64) -> bool {
        (val < 0.0 && self.model.row_upper_[row as usize] != kHighsInf)
            || (val > 0.0 && self.model.row_lower_[row as usize] != -kHighsInf)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::yieldsImpliedUpperBound (lines 259-261, 3 lines)
    pub fn yields_implied_upper_bound(&self, row: HighsInt, val: f64) -> bool {
        self.yields_implied_lower_bound(row, -val)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isImpliedEquationAtLower (lines 263-271, 9 lines)
    /// HiGHS 1.16 (#3385): HPresolve.cpp::isImpliedEquationAtLower (lines 351-357)
    pub fn is_implied_equation_at_lower(&self, row: HighsInt) -> bool {
        // HiGHS 1.16 #3385 (h16_impl_eq): HPresolve.cpp lines 351-357 (1.16.0)
        if self.options.h16_impl_eq {
            return self.model.row_lower_[row as usize] != -kHighsInf
                && self.impl_row_dual_lower[row as usize]
                    > self.options.dual_feasibility_tolerance;
        }
        // if the implied lower bound on a row dual is strictly positive then the row
        // is an implied equation (using its lower bound) due to complementary
        // slackness
        let is_lbnd_positive =
            self.impl_row_dual_lower[row as usize] > self.options.dual_feasibility_tolerance;
        debug_assert!(!is_lbnd_positive || self.model.row_lower_[row as usize] != -kHighsInf);
        is_lbnd_positive
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isImpliedEquationAtUpper (lines 273-281, 9 lines)
    /// HiGHS 1.16 (#3385): HPresolve.cpp::isImpliedEquationAtUpper (lines 359-365)
    pub fn is_implied_equation_at_upper(&self, row: HighsInt) -> bool {
        // HiGHS 1.16 #3385 (h16_impl_eq): HPresolve.cpp lines 359-365 (1.16.0)
        if self.options.h16_impl_eq {
            return self.model.row_upper_[row as usize] != kHighsInf
                && self.impl_row_dual_upper[row as usize]
                    < -self.options.dual_feasibility_tolerance;
        }
        // if the implied upper bound on a row dual is strictly negative then the row
        // is an implied equation (using its upper bound) due to complementary
        // slackness
        let is_ubnd_negative =
            self.impl_row_dual_upper[row as usize] < -self.options.dual_feasibility_tolerance;
        debug_assert!(!is_ubnd_negative || self.model.row_upper_[row as usize] != kHighsInf);
        is_ubnd_negative
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isImpliedIntegral (lines 283-357, 75 lines)
    ///
    /// Reached on LP only for `kInteger` columns, i.e. never (integrality is all
    /// `kContinuous` after `ok_set_input`); translate it anyway.
    pub fn is_implied_integral(&mut self, col: HighsInt) -> StatusResult {
        // check if the integer constraint on a variable is implied by the model
        debug_assert!(self.model.integrality_[col as usize] == HighsVarType::kInteger);
        let fractionality = |x: f64| (x - x.round()).abs();

        let mut run_dual_detection = true;

        for nz in self.get_column_vector(col).iter(self) {
            let row_idx = nz.index() as usize;
            // if not all other columns are integer, skip row and also do not try the
            // dual detection in the second loop as it must hold for all rows
            if self.rowsize[row_idx] < 2 || self.rowsize_integer[row_idx] < self.rowsize[row_idx] {
                run_dual_detection = false;
                continue;
            }

            let row_lower = if self.is_implied_equation_at_upper(nz.index()) {
                self.model.row_upper_[row_idx]
            } else {
                self.model.row_lower_[row_idx]
            };

            let row_upper = if self.is_implied_equation_at_lower(nz.index()) {
                self.model.row_lower_[row_idx]
            } else {
                self.model.row_upper_[row_idx]
            };

            if row_upper == row_lower {
                // if there is an equation the dual detection does not need to be tried
                run_dual_detection = false;
                let scale = 1.0 / nz.value();

                if !self.row_coefficients_integral(nz.index(), scale) {
                    continue;
                }

                if fractionality(row_lower * scale) > self.primal_feastol {
                    return StatusResult::from(HPresolveResult::kPrimalInfeasible);
                }

                return StatusResult::from(true);
            }
        }

        if !run_dual_detection {
            return StatusResult::from(false);
        }

        let h16_row_duals = self.options.h16_row_duals;
        let mut it = self.get_column_vector(col);
        while let Some(nz) = it.next_nz(self) {
            let row_idx = nz.index() as usize;
            let scale = 1.0 / nz.value();
            // if row coefficients are not integral, variable is not (implied) integral
            if !self.row_coefficients_integral(nz.index(), scale) {
                return StatusResult::from(false);
            }
            if self.model.row_upper_[row_idx] != kHighsInf {
                // right-hand side: scale, round down and unscale again
                let r_upper = nz.value().abs()
                    * (self.model.row_upper_[row_idx] * scale.abs() + self.primal_feastol).floor();
                // check if modification is large enough
                if (self.model.row_upper_[row_idx] - r_upper).abs()
                    > self.options.small_matrix_value
                {
                    if h16_row_duals {
                        // HiGHS 1.16 #3270 (h16_row_duals): HPresolveImpliedInts.cpp line 65 (1.16.0)
                        self.h16_change_row_upper(nz.index(), r_upper, false);
                    } else {
                        // update right-hand side and mark row as changed
                        self.model.row_upper_[row_idx] = r_upper;
                        self.mark_changed_row(nz.index());
                        // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 381 at f2ae2b320 (1.15.1 336)
                        self.x1_row_modified(nz.index());
                    }
                }
            }
            if self.model.row_lower_[row_idx] != -kHighsInf {
                // left-hand side: scale, round up and unscale again
                let r_lower = nz.value().abs()
                    * (self.model.row_lower_[row_idx] * scale.abs() - self.primal_feastol).ceil();
                // check if modification is large enough
                if (self.model.row_lower_[row_idx] - r_lower).abs()
                    > self.options.small_matrix_value
                {
                    if h16_row_duals {
                        // HiGHS 1.16 #3270 (h16_row_duals): HPresolveImpliedInts.cpp line 76 (1.16.0)
                        self.h16_change_row_lower(nz.index(), r_lower, false);
                    } else {
                        // update left-hand side and mark row as changed
                        self.model.row_lower_[row_idx] = r_lower;
                        self.mark_changed_row(nz.index());
                        // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 392 at f2ae2b320 (1.15.1 350)
                        self.x1_row_modified(nz.index());
                    }
                }
            }
        }

        StatusResult::from(true)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::isImpliedInteger (lines 359-418, 60 lines)
    ///
    /// LP-reachable: `convertImpliedInteger` calls it for every continuous column.
    pub fn is_implied_integer(&self, col: HighsInt) -> StatusResult {
        // check if a continuous variable is implied integer
        debug_assert!(self.model.integrality_[col as usize] == HighsVarType::kContinuous);
        let fractionality = |x: f64| (x - x.round()).abs();

        let mut run_dual_detection = true;

        for nz in self.get_column_vector(col).iter(self) {
            let row_idx = nz.index() as usize;
            // if not all other columns are integer, skip row and also do not try the
            // dual detection in the second loop as it must hold for all rows
            if self.rowsize[row_idx] < 2
                || self.rowsize_integer[row_idx] + self.rowsize_impl_int[row_idx]
                    < self.rowsize[row_idx] - 1
            {
                run_dual_detection = false;
                continue;
            }

            let row_lower = if self.is_implied_equation_at_upper(nz.index()) {
                self.model.row_upper_[row_idx]
            } else {
                self.model.row_lower_[row_idx]
            };

            let row_upper = if self.is_implied_equation_at_lower(nz.index()) {
                self.model.row_lower_[row_idx]
            } else {
                self.model.row_upper_[row_idx]
            };

            if row_upper == row_lower {
                // if there is an equation the dual detection does not need to be tried
                run_dual_detection = false;
                let scale = 1.0 / nz.value();

                if fractionality(row_lower * scale) > self.primal_feastol {
                    continue;
                }

                if !self.row_coefficients_integral(nz.index(), scale) {
                    continue;
                }

                return StatusResult::from(true);
            }
        }

        if !run_dual_detection {
            return StatusResult::from(false);
        }

        if (self.model.col_lower_[col as usize] != -kHighsInf
            && fractionality(self.model.col_lower_[col as usize]) > self.options.small_matrix_value)
            || (self.model.col_upper_[col as usize] != kHighsInf
                && fractionality(self.model.col_upper_[col as usize])
                    > self.options.small_matrix_value)
        {
            return StatusResult::from(false);
        }

        for nz in self.get_column_vector(col).iter(self) {
            let row_idx = nz.index() as usize;
            let scale = 1.0 / nz.value();
            if self.model.row_upper_[row_idx] != kHighsInf
                && fractionality(self.model.row_upper_[row_idx] * scale) > self.primal_feastol
            {
                return StatusResult::from(false);
            }

            if self.model.row_lower_[row_idx] != -kHighsInf
                && fractionality(self.model.row_lower_[row_idx] * scale) > self.primal_feastol
            {
                return StatusResult::from(false);
            }

            if !self.row_coefficients_integral(nz.index(), scale) {
                return StatusResult::from(false);
            }
        }

        StatusResult::from(true)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::convertImpliedInteger (lines 420-450, 31 lines)
    /// HiGHS 1.16 (#3276): HPresolveImpliedInts.cpp::convertImpliedInteger (lines 144-184)
    ///
    /// C++ default arguments: `row = -1`, `skipInputChecks = false`.
    /// Flag `h16_implint_bounds` (off = 1.15.1 bounds, unchanged): strengthen the column
    /// bounds from the implied bounds before `changeColBounds` rounds them.
    pub fn convert_implied_integer(
        &mut self,
        col: HighsInt,
        row: HighsInt,
        skip_input_checks: bool,
    ) -> StatusResult {
        // return if column was deleted
        if self.col_deleted[col as usize] != 0 {
            return StatusResult::from(false);
        }

        // return if column is not continuous or cannot be converted to an implied
        // integer
        if !skip_input_checks {
            if self.model.integrality_[col as usize] != HighsVarType::kContinuous {
                return StatusResult::from(false);
            }
            let implied_integer = self.is_implied_integer(col);
            if !implied_integer.flag() {
                return implied_integer;
            }
        }

        // convert to implied integer
        self.model.integrality_[col as usize] = HighsVarType::kImplicitInteger;

        if row != -1 {
            // use row index supplied by caller (e.g. singleton)
            self.rowsize_impl_int[row as usize] += 1;
        } else {
            // iterate over rows
            let mut it = self.get_column_vector(col);
            while let Some(nonzero) = it.next_nz(self) {
                self.rowsize_impl_int[nonzero.index() as usize] += 1;
            }
        }

        // HiGHS 1.16 #3276 (h16_implint_bounds): HPresolveImpliedInts.cpp::convertImpliedInteger lines 171-181 (1.16.0)
        // round and update bounds
        let (col_lower, col_upper) = if self.options.h16_implint_bounds {
            // Potentially strengthen existing bound using implied bounds
            // changeColBounds will perform rounding
            let mut new_lower = self.model.col_lower_[col as usize];
            let mut new_upper = self.model.col_upper_[col as usize];
            if self.impl_col_lower[col as usize] > new_lower + self.primal_feastol {
                new_lower = self.impl_col_lower[col as usize];
            }
            if self.impl_col_upper[col as usize] < new_upper - self.primal_feastol {
                new_upper = self.impl_col_upper[col as usize];
            }
            (new_lower, new_upper)
        } else {
            (
                self.model.col_lower_[col as usize],
                self.model.col_upper_[col as usize],
            )
        };
        StatusResult::from(self.change_col_bounds(col, col_lower, col_upper))
    }
}

#[cfg(test)]
mod h16_impl_eq_tests {
    use super::*;

    /// One-row presolve with the given row bounds and implied row-dual bounds.
    fn presolve_with_row(
        row_lower: f64,
        row_upper: f64,
        dual_lower: f64,
        dual_upper: f64,
        flag: bool,
    ) -> HPresolve {
        let mut presolve = HPresolve::default();
        presolve.options.h16_impl_eq = flag;
        presolve.options.dual_feasibility_tolerance = 1e-7;
        presolve.model.num_row_ = 1;
        presolve.model.num_col_ = 1;
        presolve.model.row_lower_ = vec![row_lower];
        presolve.model.row_upper_ = vec![row_upper];
        presolve.impl_row_dual_lower = vec![dual_lower];
        presolve.impl_row_dual_upper = vec![dual_upper];
        presolve
    }

    #[test]
    fn h16_impl_eq_on_unbounded_row_lower_is_not_equation() {
        // Row lower is -inf and the implied row-dual lower bound is positive: the 1.16.0
        // predicate returns false, so no implied equation at the (infinite) lower bound.
        let presolve = presolve_with_row(-kHighsInf, kHighsInf, 1.0, 0.0, true);
        assert!(!presolve.is_implied_equation_at_lower(0));
    }

    #[test]
    fn h16_impl_eq_on_finite_row_lower_is_equation() {
        let presolve = presolve_with_row(1.0, kHighsInf, 1.0, 0.0, true);
        assert!(presolve.is_implied_equation_at_lower(0));
    }

    #[test]
    fn h16_impl_eq_on_nonpositive_dual_is_not_equation() {
        let presolve = presolve_with_row(1.0, kHighsInf, 0.0, 0.0, true);
        assert!(!presolve.is_implied_equation_at_lower(0));
    }

    #[test]
    fn h16_impl_eq_on_unbounded_row_upper_is_not_equation() {
        // Row upper is +inf and the implied row-dual upper bound is negative.
        let presolve = presolve_with_row(-kHighsInf, kHighsInf, 0.0, -1.0, true);
        assert!(!presolve.is_implied_equation_at_upper(0));
    }

    #[test]
    fn h16_impl_eq_on_finite_row_upper_is_equation() {
        let presolve = presolve_with_row(-kHighsInf, 1.0, 0.0, -1.0, true);
        assert!(presolve.is_implied_equation_at_upper(0));
    }

    /// The 1.15.1 body has a debug_assert that fires on the unbounded case; it is
    /// compiled out in release, where the lib tests run, so the off path returns true.
    #[cfg(not(debug_assertions))]
    #[test]
    fn h16_impl_eq_off_unbounded_row_lower_matches_1_15_1() {
        let presolve = presolve_with_row(-kHighsInf, kHighsInf, 1.0, 0.0, false);
        assert!(presolve.is_implied_equation_at_lower(0));
        let presolve = presolve_with_row(-kHighsInf, kHighsInf, 0.0, -1.0, false);
        assert!(presolve.is_implied_equation_at_upper(0));
    }

    #[test]
    fn h16_impl_eq_off_finite_row_matches_1_15_1() {
        let presolve = presolve_with_row(1.0, kHighsInf, 1.0, 0.0, false);
        assert!(presolve.is_implied_equation_at_lower(0));
        let presolve = presolve_with_row(-kHighsInf, 1.0, 0.0, -1.0, false);
        assert!(presolve.is_implied_equation_at_upper(0));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// One continuous column `x` in `[0, 10]` with the row `x <= 10`, built by hand.
    fn one_col_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![10.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        lp.integrality_ = vec![HighsVarType::kContinuous];
        lp
    }

    /// Port of the 1.16.0 `issue-3274-1` case (TestMipSolver.cpp, #3276):
    /// three integer columns and one continuous column, three rows.
    fn issue_3274_1_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 3;
        lp.col_cost_ = vec![0.0, 1.0, 0.0, 0.0];
        lp.col_lower_ = vec![-kHighsInf, -kHighsInf, 0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf, kHighsInf, kHighsInf];
        lp.row_lower_ = vec![0.0, 0.0, -12.0];
        lp.row_upper_ = vec![kHighsInf, 0.0, kHighsInf];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 1, 2, 5, 7];
        lp.a_matrix_.index_ = vec![1, 0, 0, 1, 2, 1, 2];
        lp.a_matrix_.value_ = vec![1.0, 1.0, -1.0, 1.0, -68.0, -1.0, -1.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
        ];
        lp
    }

    /// Solves `lp` (default presolve on) and returns the objective.
    fn solve_objective(lp: HighsLp, on: bool) -> f64 {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_implint_bounds", on),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        highs.get_info().objective_function_value
    }

    /// Port of the 1.16.0 `issue-3274-2` case (TestMipSolver.cpp, #3276):
    /// a semicontinuous column, two integer columns, two rows.
    fn issue_3274_2_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![1.0, 0.0, 0.0];
        lp.col_lower_ = vec![1.0, -kHighsInf, -kHighsInf];
        lp.col_upper_ = vec![3.0, 1.0, kHighsInf];
        lp.row_lower_ = vec![13.0, 0.0];
        lp.row_upper_ = vec![kHighsInf, 0.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 2, 4, 5];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 1];
        lp.a_matrix_.value_ = vec![9.0, 1.0, 7.0, -1.0, -1.0];
        lp.integrality_ = vec![
            HighsVarType::kSemiContinuous,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp
    }

    /// Runs `convert_implied_integer` on column 0 whose implied bounds are `[2.5, 7.5]`
    /// (tighter than `[0, 10]`) and returns the column bounds afterwards.
    fn implied_integer_bounds(on: bool) -> (f64, f64) {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.h16_implint_bounds = on;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(
            one_col_lp(),
            &options,
            options.presolve_reduction_limit,
            None
        ));
        presolve.impl_col_lower[0] = 2.5;
        presolve.impl_col_upper[0] = 7.5;
        assert!(presolve.convert_implied_integer(0, -1, true).flag());
        assert_eq!(
            presolve.model.integrality_[0],
            HighsVarType::kImplicitInteger
        );
        (presolve.model.col_lower_[0], presolve.model.col_upper_[0])
    }

    #[test]
    fn h16_implint_bounds_on_uses_implied_bounds() {
        // on: the implied bounds [2.5, 7.5] strengthen [0, 10]; changeColBounds rounds them
        // to the integer bounds [3, 7].
        assert_eq!(implied_integer_bounds(true), (3.0, 7.0));
    }

    #[test]
    fn h16_implint_bounds_off_keeps_col_bounds() {
        // off: 1.15.1 ignores the implied bounds; the column keeps [0, 10].
        assert_eq!(implied_integer_bounds(false), (0.0, 10.0));
    }

    #[test]
    fn h16_implint_bounds_solve_issue_3274_1() {
        // 1.16.0 (TestMipSolver.cpp issue-3274-1): optimum 0.0. Off: the 1.15.1 C++ result
        // with presolve on is 1.0 (checked against the 1.15.1 library).
        for (on, expected) in [(true, 0.0), (false, 1.0)] {
            let obj = solve_objective(issue_3274_1_lp(), on);
            assert!(
                (obj - expected).abs() <= 1e-6,
                "on={on} objective {obj}, expected {expected}"
            );
        }
    }

    #[test]
    fn h16_implint_bounds_solve_issue_3274_2() {
        // 1.16.0 (TestMipSolver.cpp issue-3274-2): optimum 1.0. Off: the 1.15.1 C++ result
        // with presolve on is 2.0 (checked against the 1.15.1 library), which #3276 fixes.
        for (on, expected) in [(true, 1.0), (false, 2.0)] {
            let obj = solve_objective(issue_3274_2_lp(), on);
            assert!(
                (obj - expected).abs() <= 1e-6,
                "on={on} objective {obj}, expected {expected}"
            );
        }
    }
}
