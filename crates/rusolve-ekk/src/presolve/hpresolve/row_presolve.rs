//! HiGHS source: presolve/HPresolve.cpp (lines 3440-4378)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! rowPresolve. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::rowPresolve (lines 3440-4378, 939 lines)
    ///
    /// The largest function (939 lines). `std::set<double> scalars = {1.0}` (C++ 4064) ->
    /// `BTreeSet<OrderedF64>` initialised with `OrderedF64(1.0)`; iterate in set order.
    /// `computeDynamism(getStoredRow())` -> `self.compute_dynamism(self.get_stored_row())`.
    /// Lambdas capturing `this`: see the module doc.
    pub fn row_presolve(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row: HighsInt,
    ) -> HPresolveResult {
        debug_assert!(self.row_deleted[row as usize] == 0);

        let logging_on = self.analysis_.logging_on_;
        // HiGHS 1.16 #3392, #3260 (h16_row_presolve): read the option once, at the top of rowPresolve
        let h16_row_presolve = self.options.h16_row_presolve;

        // C++ lambdas that capture `this` are nested fns taking `this` explicitly
        // (see the module doc of `hpresolve/mod.rs`).

        /// C++ `computeOffset` lambda.
        fn compute_offset(model: &HighsLp, col: HighsInt, val: f64) -> HighsCDouble {
            val.abs()
                * (HighsCDouble::from(model.col_upper_[col as usize])
                    - HighsCDouble::from(model.col_lower_[col as usize]))
        }

        /// C++ `checkRowInfeasible` lambda.
        fn check_row_infeasible(this: &HPresolve, row: HighsInt) -> HPresolveResult {
            if this.implied_row_bounds.get_sum_lower(row, 0.0)
                > this.model.row_upper_[row as usize] + this.primal_feastol
                || this.implied_row_bounds.get_sum_upper(row, 0.0)
                    < this.model.row_lower_[row as usize] - this.primal_feastol
            {
                return HPresolveResult::kPrimalInfeasible;
            }
            HPresolveResult::kOk
        }

        /// HiGHS 1.16 (#3392): HPresolve.cpp::HPresolve::rowPresolve, lambda
        /// checkRowDualInfeasible (lines 2523-2530, 1.16.0). Called only when
        /// `h16_row_presolve` is set.
        fn check_row_dual_infeasible(this: &HPresolve, row: HighsInt) -> HPresolveResult {
            let ri = row as usize;
            if (this.impl_row_dual_lower[ri] > this.options.dual_feasibility_tolerance
                && this.model.row_lower_[ri] == -kHighsInf)
                || (this.impl_row_dual_upper[ri] < -this.options.dual_feasibility_tolerance
                    && this.model.row_upper_[ri] == kHighsInf)
            {
                return HPresolveResult::kDualInfeasible;
            }
            HPresolveResult::kOk
        }

        /// C++ `checkRowRedundant` lambda.
        fn check_row_redundant(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
            row: HighsInt,
            logging_on: bool,
        ) -> HPresolveResult {
            if this.is_redundant(row) {
                // row is redundant
                let presolve_rule: HighsInt = if this.rowsize[row as usize] != 0 {
                    kPresolveRuleRedundantRow
                } else {
                    kPresolveRuleEmptyRow
                };
                if logging_on {
                    this.analysis_start_presolve_rule_log(presolve_rule);
                }
                postsolve_stack.redundant_row(row);
                this.remove_row(row);
                this.analysis_.logging_on_ = logging_on;
                if logging_on {
                    this.analysis_stop_presolve_rule_log(presolve_rule);
                }
                return this.check_limits(postsolve_stack);
            }
            HPresolveResult::kOk
        }

        /// C++ `checkRedundantBounds` lambda.
        fn check_redundant_bounds(
            this: &mut HPresolve,
            col: HighsInt,
            row: HighsInt,
        ) -> HPresolveResult {
            // check if column singleton has redundant bounds
            let ci = col as usize;
            debug_assert!(this.model.col_cost_[ci] != 0.0);
            if this.colsize[ci] != 1 {
                return HPresolveResult::kOk;
            }
            if this.model.col_cost_[ci] > 0.0 {
                debug_assert!(
                    this.model.col_lower_[ci] == -kHighsInf
                        || !this.is_lower_implied(col)
                        || this.col_lower_source[ci] == row
                );
                if this.model.col_lower_[ci] > this.impl_col_lower[ci] - this.primal_feastol {
                    hpresolve_checked_call!(this.change_col_lower(col, -kHighsInf));
                }
            } else {
                debug_assert!(
                    this.model.col_upper_[ci] == kHighsInf
                        || !this.is_upper_implied(col)
                        || this.col_upper_source[ci] == row
                );
                if this.model.col_upper_[ci] < this.impl_col_upper[ci] + this.primal_feastol {
                    hpresolve_checked_call!(this.change_col_upper(col, kHighsInf));
                }
            }
            HPresolveResult::kOk
        }

        /// C++ `degree1Tests` lambda.
        fn degree1_tests(
            this: &mut HPresolve,
            col: HighsInt,
            val: f64,
            direction: HighsInt,
            row_activity_bound: f64,
            row_bound: f64,
        ) -> HPresolveResult {
            // return if row is still feasible
            if (direction as f64) * row_activity_bound
                >= (direction as f64) * row_bound - this.primal_feastol
            {
                return HPresolveResult::kOk;
            }
            // tighten bound
            if (direction as f64) * val > 0.0 {
                let new_lower = this.model.col_lower_[col as usize] + 1.0;
                hpresolve_checked_call!(this.change_col_lower(col, new_lower));
            } else {
                let new_upper = this.model.col_upper_[col as usize] - 1.0;
                hpresolve_checked_call!(this.change_col_upper(col, new_upper));
            }
            HPresolveResult::kOk
        }

        /// C++ `remDoubletonEq` lambda (inside the simple probing loop).
        /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::rowPresolve (lines 4229-4239)
        /// `x1b` = `mip_presolve_clique_table`; off, the result is always `kOk`.
        #[allow(clippy::too_many_arguments)]
        fn rem_doubleton_eq(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
            col: HighsInt,
            bin_col: HighsInt,
            direction: HighsInt,
            lower: f64,
            upper: f64,
            x1b: bool,
        ) -> HPresolveResult {
            let bound = if direction >= 0 { upper } else { lower };
            let scale = (direction as f64) * (lower - upper);
            let offset = bound - this.model.col_lower_[bin_col as usize] * scale;
            postsolve_stack.doubleton_equation(
                -1,
                col,
                bin_col,
                1.0,
                -scale,
                offset,
                lower,
                upper,
                0.0,
                false,
                false,
                RowType::kEq,
                HighsEmptySlice,
            );
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 4238 at 7b890c255
            if x1b {
                return this.substitute_x1b(col, bin_col, offset, scale, -1);
            }
            this.substitute_col(col, bin_col, offset, scale);
            HPresolveResult::kOk
        }

        // check for infeasibility
        hpresolve_checked_call!(check_row_infeasible(self, row));
        // HiGHS 1.16 #3392 (h16_row_presolve): HPresolve.cpp lines 2553-2554 (1.16.0)
        if h16_row_presolve {
            hpresolve_checked_call!(check_row_dual_infeasible(self, row));
        }

        // handle special cases directly via a call to the specialized procedure
        if self.rowsize[row as usize] == 1 {
            return self.singleton_row(postsolve_stack, row);
        }

        // check for redundancy
        hpresolve_checked_call!(check_row_redundant(self, postsolve_stack, row, logging_on));
        if self.row_deleted[row as usize] != 0 {
            return HPresolveResult::kOk;
        }

        if self.rowsize_integer[row as usize] != 0 || self.rowsize_impl_int[row as usize] != 0 {
            // check if setting variable to opposite bound (with respect to bound used
            // in calculation of bounds on row activity) makes row infeasible.
            let mut it = self.get_row_vector(row);
            while let Some(nonzero) = it.next_nz(self) {
                // get column index and coefficient
                let col = nonzero.index();
                let val = nonzero.value();

                // skip continuous variables and non-binary variables
                if self.model.integrality_[col as usize] == HighsVarType::kContinuous
                    || self.model.col_upper_[col as usize]
                        != self.model.col_lower_[col as usize] + 1.0
                {
                    continue;
                }

                // perform tests
                let activity_upper = self
                    .implied_row_bounds
                    .get_sum_upper_orig(row, -compute_offset(&self.model, col, val));
                let row_lower_v = self.model.row_lower_[row as usize];
                hpresolve_checked_call!(degree1_tests(
                    self,
                    col,
                    val,
                    1,
                    activity_upper,
                    row_lower_v
                ));
                let activity_lower = self
                    .implied_row_bounds
                    .get_sum_lower_orig(row, compute_offset(&self.model, col, val));
                let row_upper_v = self.model.row_upper_[row as usize];
                hpresolve_checked_call!(degree1_tests(
                    self,
                    col,
                    val,
                    -1,
                    activity_lower,
                    row_upper_v
                ));
            }
        }

        // check for redundancy again
        hpresolve_checked_call!(check_row_redundant(self, postsolve_stack, row, logging_on));
        if self.row_deleted[row as usize] != 0 {
            return HPresolveResult::kOk;
        }

        // Store original bounds
        let orig_row_upper = self.model.row_upper_[row as usize];
        let orig_row_lower = self.model.row_lower_[row as usize];

        if !self.is_equation(row) && self.options.h16_row_duals {
            // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp lines 2646-2663 (1.16.0)
            // Convert to equality constraint and record for dual postsolve
            if self.is_implied_equation_at_lower(row) {
                // save source before changeRowUpper resets it
                let dual_lower_source_col = self.row_dual_lower_source[row as usize];
                let new_upper = self.model.row_lower_[row as usize];
                self.h16_change_row_upper(row, new_upper, false);
                self.h16_record_implied_equation(postsolve_stack, row, true);
                if self.mipsolver.is_none() {
                    hpresolve_checked_call!(check_redundant_bounds(
                        self,
                        dual_lower_source_col,
                        row
                    ));
                }
            } else if self.is_implied_equation_at_upper(row) {
                // save source before changeRowLower resets it
                let dual_upper_source_col = self.row_dual_upper_source[row as usize];
                let new_lower = self.model.row_upper_[row as usize];
                self.h16_change_row_lower(row, new_lower, false);
                self.h16_record_implied_equation(postsolve_stack, row, false);
                if self.mipsolver.is_none() {
                    hpresolve_checked_call!(check_redundant_bounds(
                        self,
                        dual_upper_source_col,
                        row
                    ));
                }
            }
        } else if !self.is_equation(row) {
            if self.is_implied_equation_at_lower(row) {
                // Convert to equality constraint (note that currently postsolve will not
                // know about this conversion)
                self.model.row_upper_[row as usize] = self.model.row_lower_[row as usize];
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 4086 at f2ae2b320 (1.15.1 3571)
                self.x1_row_modified(row);
                // Since row upper bound is now finite, lower bound on row dual is
                // -kHighsInf
                self.change_row_dual_lower(row, -kHighsInf);
                if self.mipsolver.is_none() {
                    let src = self.row_dual_lower_source[row as usize];
                    hpresolve_checked_call!(check_redundant_bounds(self, src, row));
                }
            } else if self.is_implied_equation_at_upper(row) {
                // Convert to equality constraint (note that currently postsolve will not
                // know about this conversion)
                self.model.row_lower_[row as usize] = self.model.row_upper_[row as usize];
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 4093 at f2ae2b320 (1.15.1 3581)
                self.x1_row_modified(row);
                // Since row lower bound is now finite, upper bound on row dual is
                // kHighsInf
                self.change_row_dual_upper(row, kHighsInf);
                if self.mipsolver.is_none() {
                    let src = self.row_dual_upper_source[row as usize];
                    hpresolve_checked_call!(check_redundant_bounds(self, src, row));
                }
            }
        }

        // Get row bounds
        let row_upper = self.model.row_upper_[row as usize];
        let row_lower = self.model.row_lower_[row as usize];

        // Handle doubleton equations
        if self.rowsize[row as usize] == 2
            && row_lower == row_upper
            && self.analysis_.allow_rule_[kPresolveRuleDoubletonEquation as usize]
        {
            let row_type = if orig_row_lower == orig_row_upper {
                RowType::kEq
            } else if orig_row_upper != kHighsInf {
                RowType::kLeq
            } else {
                debug_assert!(orig_row_lower != -kHighsInf);
                RowType::kGeq
            };
            return self.doubleton_eq(postsolve_stack, row, row_type);
        }

        /// C++ `roundRhs` lambda (inside `checkScaleRow`).
        fn round_rhs(
            this: &HPresolve,
            rhs: HighsCDouble,
            rounded_rhs: &mut HighsCDouble,
            fraction_rhs: &mut HighsCDouble,
            min_rhs_tightening: f64,
            rhs_tightened: &mut bool,
            direction: HighsInt,
        ) {
            // round rhs
            *rounded_rhs =
                (direction as f64) * ((direction as f64) * rhs + this.primal_feastol).floor();
            // compute fractional part
            *fraction_rhs = (direction as f64) * (rhs - *rounded_rhs);
            // check if tightened
            *rhs_tightened = *fraction_rhs >= min_rhs_tightening - this.options.small_matrix_value;
        }

        /// C++ `checkScaleRow` lambda.
        #[allow(clippy::too_many_arguments)]
        fn check_scale_row(
            this: &HPresolve,
            row: HighsInt,
            row_coefs: &mut [f64],
            row_index: &[HighsInt],
            rounded_lhs: &mut HighsCDouble,
            rounded_rhs: &mut HighsCDouble,
            fraction_lhs: &mut HighsCDouble,
            fraction_rhs: &mut HighsCDouble,
            min_lhs_tightening: &mut f64,
            min_rhs_tightening: &mut f64,
            max_val: &mut f64,
            lhs_tightened: &mut bool,
            rhs_tightened: &mut bool,
            is_infeasible: &mut bool,
            int_scale: f64,
        ) -> bool {
            let small = this.options.small_matrix_value;
            let mut lhs = HighsCDouble::from(this.model.row_lower_[row as usize]);
            let mut rhs = HighsCDouble::from(this.model.row_upper_[row as usize]);
            let lhs_finite = lhs != -kHighsInf;
            let rhs_finite = rhs != kHighsInf;
            if lhs_finite {
                lhs = lhs * int_scale;
            }
            if rhs_finite {
                rhs = rhs * int_scale;
            }
            *rounded_lhs = HighsCDouble::from(-kHighsInf);
            *rounded_rhs = HighsCDouble::from(kHighsInf);
            *fraction_lhs = HighsCDouble::from(0.0);
            *fraction_rhs = HighsCDouble::from(0.0);
            *min_rhs_tightening = 0.0;
            *min_lhs_tightening = 0.0;
            *max_val = 0.0;
            *lhs_tightened = false;
            *rhs_tightened = false;
            *is_infeasible = false;
            for i in 0..row_coefs.len() {
                // computed scaled coefficient
                let scale_coef = HighsCDouble::from(row_coefs[i]) * int_scale;
                // round to the nearest integer
                let int_coef = (scale_coef + 0.5).floor();
                // compute difference
                let coef_delta = int_coef - scale_coef;
                // store integral coefficient and maximum absolute value
                row_coefs[i] = f64::from(int_coef);
                *max_val = cpp_max(row_coefs[i].abs(), *max_val);
                // get column upper bound
                let ub = this.model.col_upper_[row_index[i] as usize];
                if coef_delta < -small {
                    // for the >= side of the constraint a smaller coefficient is
                    // stronger: Therefore we relax the left hand side using the
                    // bound constraint, if the bound is infinite, abort
                    if lhs_finite {
                        if ub == kHighsInf {
                            return false;
                        }
                        lhs += ub * coef_delta;
                    }
                    *min_rhs_tightening = cpp_max(-f64::from(coef_delta), *min_rhs_tightening);
                } else if coef_delta > small {
                    if rhs_finite {
                        if ub == kHighsInf {
                            return false;
                        }
                        rhs += ub * coef_delta;
                    }
                    // the coefficient was relaxed regarding the rows lower bound.
                    // Therefore the lower bound should be tightened by at least
                    // this amount for the scaled constraint to dominate the
                    // unscaled constraint be rounded by at least this value
                    *min_lhs_tightening = cpp_max(f64::from(coef_delta), *min_lhs_tightening);
                }
            }
            // round left-hand and right-hand sides
            if lhs_finite {
                round_rhs(
                    this,
                    lhs,
                    rounded_lhs,
                    fraction_lhs,
                    *min_lhs_tightening,
                    lhs_tightened,
                    -1,
                );
            }
            if rhs_finite {
                round_rhs(
                    this,
                    rhs,
                    rounded_rhs,
                    fraction_rhs,
                    *min_rhs_tightening,
                    rhs_tightened,
                    1,
                );
            }
            // check for infeasibility
            *is_infeasible = lhs_finite && rhs_finite && *rounded_rhs < *rounded_lhs - 0.5;
            true
        }

        /// C++ `scaleRow` lambda (inside the inequality branch).
        #[allow(clippy::too_many_arguments)]
        fn lam_scale_row(
            this: &mut HPresolve,
            row: HighsInt,
            rounded_lhs: HighsCDouble,
            rounded_rhs: HighsCDouble,
            scalar: f64,
            check_delta: bool,
            row_coefs: &[f64],
            row_index: &[HighsInt],
        ) {
            // scale the row
            if this.options.h16_row_duals {
                // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp lines 3077-3080 (1.16.0)
                if rounded_lhs != -kHighsInf {
                    this.h16_change_row_lower(row, f64::from(rounded_lhs / scalar), false);
                }
                if rounded_rhs != kHighsInf {
                    this.h16_change_row_upper(row, f64::from(rounded_rhs / scalar), false);
                }
            } else {
                if rounded_lhs != -kHighsInf {
                    this.model.row_lower_[row as usize] = f64::from(rounded_lhs / scalar);
                }
                if rounded_rhs != kHighsInf {
                    this.model.row_upper_[row as usize] = f64::from(rounded_rhs / scalar);
                }
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 4506-4508 at f2ae2b320 (1.15.1 3985-3987)
                this.x1_row_modified(row);
            }
            for i in 0..row_coefs.len() {
                let delta = f64::from(
                    HighsCDouble::from(row_coefs[i]) / scalar
                        - this.a_value[this.rowpositions[i] as usize],
                );
                if !check_delta || delta.abs() > this.options.small_matrix_value {
                    this.add_to_matrix(row, row_index[i], delta);
                }
            }
        }

        /// C++ `scaleRowIntVals` lambda (inside the inequality branch).
        #[allow(clippy::too_many_arguments)]
        fn lam_scale_row_int_vals(
            this: &mut HPresolve,
            row: HighsInt,
            rounded_lhs: HighsCDouble,
            rounded_rhs: HighsCDouble,
            int_scale: f64,
            max_val: f64,
            row_coefs: &[f64],
            row_index: &[HighsInt],
        ) -> bool {
            // return if maximum value in scaled row or scalar are too large
            if max_val > 1000.0 && int_scale > 100.0 {
                return false;
            }
            // the scale value is reasonably small, change the row values to be
            // integral
            lam_scale_row(
                this,
                row,
                rounded_lhs,
                rounded_rhs,
                1.0,
                false,
                row_coefs,
                row_index,
            );
            true
        }

        /// C++ `complementOrShift` lambda.
        fn complement_or_shift(
            this: &HPresolve,
            row_coefs: &[f64],
            row_index: &[HighsInt],
            direction: HighsInt,
            rhs: &mut HighsCDouble,
            min_abs_coef: &mut f64,
            max_abs_coef: &mut f64,
            dynamism: &mut f64,
        ) -> bool {
            *min_abs_coef = kHighsInf;
            *max_abs_coef = 0.0;
            for i in 0..row_coefs.len() {
                // get column index and (absolute) coefficient
                let col = row_index[i] as usize;
                let val = (direction as f64) * row_coefs[i];
                let absval = val.abs();
                // compute minimum and maximum absolute coefficients along the way
                *min_abs_coef = cpp_min(*min_abs_coef, absval);
                *max_abs_coef = cpp_max(*max_abs_coef, absval);
                if val < 0.0 && this.model.col_upper_[col] != kHighsInf {
                    // complement
                    *rhs -= val * HighsCDouble::from(this.model.col_upper_[col]);
                } else if val > 0.0 && this.model.col_lower_[col] != -kHighsInf {
                    // shift
                    *rhs -= val * HighsCDouble::from(this.model.col_lower_[col]);
                } else {
                    // unbounded variable; cannot shift or complement!
                    return false;
                }
            }
            *dynamism = *max_abs_coef / *min_abs_coef;
            true
        }

        /// C++ `undoComplementOrShift` lambda.
        fn undo_complement_or_shift(
            this: &HPresolve,
            row_coefs: &[f64],
            row_index: &[HighsInt],
            rounded_row_coefs: &mut [f64],
            direction: HighsInt,
            rounded_rhs: &mut HighsCDouble,
        ) {
            for i in 0..row_coefs.len() {
                let col = row_index[i] as usize;
                let val = (direction as f64) * row_coefs[i];
                if val < 0.0 && this.model.col_upper_[col] != kHighsInf {
                    // uncomplement
                    rounded_row_coefs[i] = -rounded_row_coefs[i];
                    *rounded_rhs +=
                        rounded_row_coefs[i] * HighsCDouble::from(this.model.col_upper_[col]);
                } else if val > 0.0 && this.model.col_lower_[col] != -kHighsInf {
                    // unshift
                    *rounded_rhs +=
                        rounded_row_coefs[i] * HighsCDouble::from(this.model.col_lower_[col]);
                }
                // flip coefficient sign for <= inequality
                rounded_row_coefs[i] *= direction as f64;
            }
            // flip rhs sign for <= inequality
            *rounded_rhs *= direction as f64;
        }

        /// C++ `roundRow` lambda.
        fn round_row(
            this: &HPresolve,
            s: f64,
            rhs: &HighsCDouble,
            rounded_rhs: &mut HighsCDouble,
            row_coefs: &[f64],
            rounded_row_coefs: &mut [f64],
        ) -> bool {
            let small = this.options.small_matrix_value;
            let mut accept = false;
            // round rhs (using feasibility tolerance)
            let scalar = HighsCDouble::from(s);
            *rounded_rhs = (*rhs * scalar - this.primal_feastol).ceil();
            if *rounded_rhs <= this.primal_feastol {
                return false;
            }
            let rhs_ratio = *rhs / *rounded_rhs;

            for i in 0..row_coefs.len() {
                // coefficient sign has not been flipped for complemented
                // variables; take absolute value of coefficient.
                let abs_coef = row_coefs[i].abs();
                // round coefficient
                rounded_row_coefs[i] = f64::from((abs_coef * scalar - kHighsTiny).ceil());
                // compare "normalised" coefficients, i.e. coefficients divided by
                // corresponding rhs.
                let threshold = f64::from(rounded_row_coefs[i] * rhs_ratio);
                // return if coefficient is weaker
                if abs_coef < threshold - small {
                    return false;
                }
                // accept rounding if at least one coefficient is improved
                accept = accept || (abs_coef > threshold + small);
            }
            accept
        }

        /// C++ `setScalars` lambda.
        fn set_scalars(scalars: &mut BTreeSet<OrderedF64>, min_abs_coef: f64, max_abs_coef: f64) {
            for t in 1..=5i32 {
                let tf = t as f64;
                scalars.insert(OrderedF64(tf / max_abs_coef));
                scalars.insert(OrderedF64(tf / min_abs_coef));
                scalars.insert(OrderedF64(((2 * t - 1) as f64) / (2.0 * min_abs_coef)));
            }
        }

        /// C++ `rowCanBeTightened` lambda.
        fn row_can_be_tightened(
            this: &HPresolve,
            row_coefs: &[f64],
            rounded_row_coefs: &mut [f64],
            scalars: &BTreeSet<OrderedF64>,
            rhs: &HighsCDouble,
            rounded_rhs: &mut HighsCDouble,
        ) -> bool {
            for s in scalars.iter() {
                if round_row(this, s.0, rhs, rounded_rhs, row_coefs, rounded_row_coefs) {
                    return true;
                }
            }
            false
        }

        /// C++ `updateRow` lambda: replace the model row by the rounded one.
        fn update_row(
            this: &mut HPresolve,
            row: HighsInt,
            row_coefs: &[f64],
            row_index: &[HighsInt],
            rounded_row_coefs: &[f64],
            direction: HighsInt,
            rounded_rhs: &HighsCDouble,
        ) {
            if this.options.h16_row_duals {
                // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp lines 3264-3267 (1.16.0)
                if direction < 0 {
                    this.h16_change_row_upper(row, f64::from(*rounded_rhs), false);
                } else {
                    this.h16_change_row_lower(row, f64::from(*rounded_rhs), false);
                }
            } else {
                if direction < 0 {
                    this.model.row_upper_[row as usize] = f64::from(*rounded_rhs);
                } else {
                    this.model.row_lower_[row as usize] = f64::from(*rounded_rhs);
                }
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 4693-4695 at f2ae2b320 (1.15.1 4172-4174)
                this.x1_row_modified(row);
            }
            for i in 0..row_coefs.len() {
                let delta = f64::from(HighsCDouble::from(rounded_row_coefs[i]) - row_coefs[i]);
                if delta.abs() > this.options.small_matrix_value {
                    this.add_to_matrix(row, row_index[i], delta);
                }
            }
        }

        /// rusolve (mip_presolve_exact_strengthen): make the strengthened coefficient equal to
        /// `target` (the nearest double to the intended `±maxAbsCoefValue`). `strengthenCoefs`
        /// adds `double(delta)` to the old coefficient, and when the old coefficient is large
        /// (a big-M of 1e8) `delta` is close to it in magnitude, so rounding `delta` to a
        /// double loses the low digits: -1e8 + double(1e8 - 5.923) is -5.9229999929666519, not
        /// -5.923. The row is then tighter than intended and can cut off the optimum (with
        /// coefficients 0.784 and 5.139 the point with all three binaries at 1 is lost). When
        /// the coefficient left by the first add is within a factor of two of `target` (the
        /// usual case), `target - current` is exact (Sterbenz lemma) and the second add lands
        /// exactly on `target`. Going through `add_to_matrix` keeps the implied-bound bookkeeping.
        fn set_strengthened_coef(this: &mut HPresolve, row: HighsInt, col: HighsInt, target: f64) {
            let pos = this.find_nonzero(row, col);
            if pos == -1 {
                return;
            }
            let correction = target - this.a_value[pos as usize];
            // Errors at round-off level (at most 7e-14 relative on 20 MIPLIB 2017 instances)
            // are left alone, so the search only changes where the cancellation lost
            // real digits (about 1e-9 on a 1e8 big-M).
            if correction.abs() > 1e-12 * target.abs().max(1.0) {
                this.add_to_matrix(row, col, correction);
            }
        }

        /// C++ `strengthenCoefs` lambda (inside the inequality branch).
        fn strengthen_coefs(
            this: &mut HPresolve,
            row: HighsInt,
            rhs: &mut HighsCDouble,
            direction: HighsInt,
            max_abs_coef_value: HighsCDouble,
        ) {
            // iterate over non-zero positions (see the C++ comment: the stored
            // positions are read by index, the slices are not used here)
            for i in 0..this.rowpositions.len() {
                let rowiter = this.rowpositions[i] as usize;
                // max. absolute coefficient should not be negative
                debug_assert!(max_abs_coef_value >= 0.0);

                // get column index and coefficient
                let col = this.a_col[rowiter];
                let val = (direction as f64) * this.a_value[rowiter];

                // get lower and upper bounds
                let col_lower =
                    this.implied_row_bounds
                        .get_impl_var_lower(&col_bound_arrays!(this), row, col);
                let col_upper =
                    this.implied_row_bounds
                        .get_impl_var_upper(&col_bound_arrays!(this), row, col);

                // skip continuous variables
                if this.model.integrality_[col as usize] == HighsVarType::kContinuous {
                    continue;
                }

                if val > max_abs_coef_value + this.primal_feastol {
                    debug_assert!(col_upper != kHighsInf);
                    // new matrix coefficient is direction * maxAbsCoefValue; subtract
                    // existing matrix coefficient to get delta
                    let delta = (direction as f64) * (max_abs_coef_value - val);
                    this.add_to_matrix(row, col, f64::from(delta));
                    if this.options.mip_presolve_exact_strengthen {
                        let target = (direction as f64) * f64::from(max_abs_coef_value);
                        set_strengthened_coef(this, row, col, target);
                    }
                    *rhs += delta * col_upper;
                } else if val < -max_abs_coef_value - this.primal_feastol {
                    debug_assert!(col_lower != -kHighsInf);
                    // new matrix coefficient is (-direction) * maxAbsCoefValue;
                    // subtract existing matrix coefficient to get delta
                    let delta = ((-direction) as f64) * (max_abs_coef_value + val);
                    this.add_to_matrix(row, col, f64::from(delta));
                    if this.options.mip_presolve_exact_strengthen {
                        let target = ((-direction) as f64) * f64::from(max_abs_coef_value);
                        set_strengthened_coef(this, row, col, target);
                    }
                    *rhs += delta * col_lower;
                }
            }
        }

        // todo: do additional single row presolve for mip here. It may assume a
        // non-redundant and non-infeasible row when considering variable and implied
        // bounds
        if self.rowsize_integer[row as usize] != 0 || self.rowsize_impl_int[row as usize] != 0 {
            if row_lower == row_upper {
                // equation
                let implied_row_lower = self.implied_row_bounds.get_sum_lower(row, 0.0);
                let implied_row_upper = self.implied_row_bounds.get_sum_upper(row, 0.0);
                // see section 3.6 "Simple probing on a single equation", Achterberg et
                // al., Presolve Reductions in Mixed Integer Programming, INFORMS Journal
                // on Computing 32(2):473-506.
                if implied_row_lower != -kHighsInf
                    && implied_row_upper != kHighsInf
                    && (implied_row_lower + implied_row_upper - 2.0 * row_upper).abs()
                        <= self.options.small_matrix_value
                {
                    // simple probing on equation case
                    let mut bin_col: HighsInt = -1;
                    let mut bin_coef = (implied_row_upper - row_upper).abs();
                    // store row
                    self.store_row(row);
                    for nonz in self.get_stored_row().iter(self) {
                        let ci = nonz.index() as usize;
                        if (nonz.value().abs() - bin_coef).abs() <= self.options.small_matrix_value
                            && self.model.integrality_[ci] == HighsVarType::kInteger
                            && (self.model.col_upper_[ci] - self.model.col_lower_[ci] - 1.0).abs()
                                <= self.primal_feastol
                        {
                            // found a binary variable that implies all other variables to be
                            // fixed when it sits at one of its bounds therefore we can
                            // substitute all other variables in the row
                            bin_col = nonz.index();
                            // store the binary coefficient with its actual sign
                            bin_coef = nonz.value();
                            break;
                        }
                    }

                    if bin_col != -1 {
                        // found binary column for substituting all other columns
                        // (the stored positions are read by index, as in the C++ range-for)
                        // HiGHS 1.16 #3260 (h16_row_presolve): HPresolve.cpp lines 2726-2735 (1.16.0):
                        // snapshot the implied variable bounds before the substitution loop
                        let mut impl_var_snapshot: Vec<(f64, f64)> = Vec::new();
                        if h16_row_presolve {
                            impl_var_snapshot.reserve(self.rowpositions.len());
                            for i in 0..self.rowpositions.len() {
                                let rowiter = self.rowpositions[i] as usize;
                                let col = self.a_col[rowiter];
                                impl_var_snapshot.push((
                                    self.implied_row_bounds.get_impl_var_lower(
                                        &col_bound_arrays!(self),
                                        row,
                                        col,
                                    ),
                                    self.implied_row_bounds.get_impl_var_upper(
                                        &col_bound_arrays!(self),
                                        row,
                                        col,
                                    ),
                                ));
                            }
                        }
                        for i in 0..self.rowpositions.len() {
                            let rowiter = self.rowpositions[i] as usize;
                            let col = self.a_col[rowiter];
                            debug_assert!(self.a_row[rowiter] == row);

                            // skip binary column
                            if col == bin_col {
                                continue;
                            }

                            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 4208-4213 at 7b890c255
                            if self.options.mip_presolve_clique_table {
                                // Use eliminateCol rather than substituteCol to avoid
                                // cascading clique fixings that could corrupt this loop
                                self.x1b_eliminate_clique_col(col);
                            }

                            // get column lower and upper bounds used to compute bounds on row
                            // activities
                            let (col_lower, col_upper) = if h16_row_presolve {
                                // HiGHS 1.16 #3260 (h16_row_presolve): HPresolve.cpp lines 2758-2759 (1.16.0)
                                (impl_var_snapshot[i].0, impl_var_snapshot[i].1)
                            } else {
                                (
                                    self.implied_row_bounds.get_impl_var_lower(
                                        &col_bound_arrays!(self),
                                        row,
                                        col,
                                    ),
                                    self.implied_row_bounds.get_impl_var_upper(
                                        &col_bound_arrays!(self),
                                        row,
                                        col,
                                    ),
                                )
                            };
                            debug_assert!(col_lower != -kHighsInf);
                            debug_assert!(col_upper != kHighsInf);

                            if self.model.col_lower_[col as usize]
                                == self.model.col_upper_[col as usize]
                            {
                                let fix_value = self.model.col_lower_[col as usize];
                                postsolve_stack.removed_fixed_col(
                                    col,
                                    fix_value,
                                    0.0,
                                    HighsEmptySlice,
                                );
                                // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 4224 at 7b890c255
                                if self.options.mip_presolve_clique_table {
                                    hpresolve_checked_call!(self.remove_fixed_col_x1b(col));
                                } else {
                                    self.remove_fixed_col(col);
                                }
                                continue;
                            }

                            // 1. binary coefficient is positive: (see C++ HPresolve.cpp)
                            // 2. the opposite case
                            let direction: HighsInt = if bin_coef.is_sign_negative()
                                == self.a_value[rowiter].is_sign_negative()
                            {
                                1
                            } else {
                                -1
                            };
                            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 4261 at 7b890c255
                            let x1b = self.options.mip_presolve_clique_table;
                            let rem_result = rem_doubleton_eq(
                                self,
                                postsolve_stack,
                                col,
                                bin_col,
                                direction,
                                col_lower,
                                col_upper,
                                x1b,
                            );
                            if x1b {
                                hpresolve_checked_call!(rem_result);
                            }
                        }

                        self.remove_row(row);
                        hpresolve_checked_call!(self.check_limits(postsolve_stack));
                        return self.remove_row_singletons(postsolve_stack);
                    }
                }

                if self.rowsize_integer[row as usize] + self.rowsize_impl_int[row as usize]
                    >= self.rowsize[row as usize] - 1
                {
                    let mut continuous_col: HighsInt = -1;
                    let mut continuous_coef: f64 = 0.0;
                    let mut row_coefs_int: Vec<f64> = Vec::new();
                    self.store_row(row);

                    for nonz in self.get_stored_row().iter(self) {
                        if self.model.integrality_[nonz.index() as usize]
                            == HighsVarType::kContinuous
                        {
                            debug_assert!(continuous_coef == 0.0);
                            continuous_coef = nonz.value();
                            continuous_col = nonz.index();
                            continue;
                        }

                        row_coefs_int.push(nonz.value());
                    }

                    if continuous_coef != 0.0 {
                        row_coefs_int.push(row_upper);

                        let int_scale: f64 = HighsIntegers::integral_scale_vec(
                            &row_coefs_int,
                            self.options.small_matrix_value,
                            self.options.small_matrix_value,
                        );

                        if int_scale != 0.0 && int_scale <= 1e3 {
                            let scale = 1.0 / (continuous_coef * int_scale).abs();
                            if scale != 1.0 {
                                // transform continuous column x to implicit integer z with x = scale * z
                                hpresolve_checked_call!(self.transform_column(
                                    postsolve_stack,
                                    continuous_col,
                                    scale,
                                    0.0
                                ));

                                hpresolve_checked_call!(
                                    self.convert_implied_integer(continuous_col, -1, true)
                                        .result()
                                );

                                if int_scale != 1.0 {
                                    self.scale_stored_row(row, int_scale, true);
                                }
                            }
                        }
                    } else {
                        let mut int_scale: f64 = HighsIntegers::integral_scale_vec(
                            &row_coefs_int,
                            self.options.small_matrix_value,
                            self.options.small_matrix_value,
                        );

                        if int_scale != 0.0 && int_scale <= 1e3 {
                            let mut rhs = row_upper * int_scale;
                            let frac_rhs: f64 = fractionality(rhs, None);
                            if frac_rhs > self.primal_feastol {
                                return HPresolveResult::kPrimalInfeasible;
                            }

                            rhs = rhs.round();

                            let mut x1_cand: HighsInt = -1;
                            let mut d: i64 = 0;

                            for i in 0..self.rowpositions.len() {
                                let scaled =
                                    int_scale * self.a_value[self.rowpositions[i] as usize];
                                let newgcd: i64 = if d == 0 {
                                    scaled.round().abs() as i64
                                } else {
                                    // C++ passes the double `std::abs(std::round(..))`
                                    // to the int64_t parameter (truncating conversion).
                                    HighsIntegers::gcd(scaled.round().abs() as i64, d)
                                };
                                if newgcd == 1 {
                                    // adding this variable would set the gcd to 1, therefore it
                                    // must be our candidate x1 for substitution. If another
                                    // candidate already exists no reduction is possible except for
                                    // scaling the equation
                                    if x1_cand != -1 {
                                        x1_cand = -1;
                                        break;
                                    }
                                    x1_cand = i as HighsInt;
                                } else {
                                    d = newgcd;
                                }
                            }

                            if x1_cand != -1 {
                                let x1_pos = self.rowpositions[x1_cand as usize];
                                let x1 = self.a_col[x1_pos as usize];
                                let rhs2 = rhs / (d as f64);
                                let frac_rhs2: f64 = fractionality(rhs2, None);
                                let mip_epsilon: f64 = self.mip().mipdata_.epsilon;
                                if frac_rhs2 <= mip_epsilon {
                                    // the right hand side is integral, so we can substitute
                                    // x1 = d * z
                                    hpresolve_checked_call!(self.transform_column(
                                        postsolve_stack,
                                        x1,
                                        d as f64,
                                        0.0
                                    ));
                                } else {
                                    // we can substitute x1 = d * z + b, with b = a1^-1 rhs (mod d)

                                    // first compute the modular multiplicative inverse of a1^-1
                                    // (mod d) of a1
                                    let mut a1: i64 =
                                        (int_scale * self.a_value[x1_pos as usize]).round() as i64;
                                    a1 = HighsIntegers::mod_i64(a1, d);
                                    let a1_inverse: i64 = HighsIntegers::modular_inverse(a1, d);

                                    // now compute b = a1^-1 rhs (mod d)
                                    let b: f64 =
                                        HighsIntegers::mod_f64(a1_inverse as f64 * rhs, d as f64);

                                    // before we substitute, we check whether the resulting variable
                                    // z is fixed after rounding its new bounds. If that is the case
                                    // we directly fix x1 instead of first substituting with d * z + b.
                                    let z_lower = ((self.model.col_lower_[x1 as usize] - b)
                                        / (d as f64)
                                        - self.primal_feastol)
                                        .ceil();
                                    let z_upper = ((self.model.col_upper_[x1 as usize] - b)
                                        / (d as f64)
                                        + self.primal_feastol)
                                        .floor();

                                    if z_lower == z_upper {
                                        // Rounded bounds are equal
                                        //
                                        // Adjust bounds if variable is fixed to a value in between
                                        // its bounds
                                        let fix_val = z_lower * (d as f64) + b;
                                        debug_assert!(
                                            fix_val
                                                > self.model.col_lower_[x1 as usize]
                                                    - self.primal_feastol
                                        );
                                        debug_assert!(
                                            fix_val
                                                < self.model.col_upper_[x1 as usize]
                                                    + self.primal_feastol
                                        );
                                        // change bounds
                                        hpresolve_checked_call!(
                                            self.change_col_bounds(x1, fix_val, fix_val)
                                        );
                                        // remove column
                                        let cost = self.model.col_cost_[x1 as usize];
                                        postsolve_stack.removed_fixed_col(
                                            x1,
                                            fix_val,
                                            cost,
                                            self.get_column_vector(x1).iter(self),
                                        );
                                        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 4410 at 7b890c255
                                        if self.options.mip_presolve_clique_table {
                                            hpresolve_checked_call!(self.remove_fixed_col_x1b(x1));
                                        } else {
                                            self.remove_fixed_col(x1);
                                        }
                                        self.rowpositions.remove(x1_cand as usize);
                                        self.sync_rowpositions_mem();
                                    } else {
                                        hpresolve_checked_call!(self.transform_column(
                                            postsolve_stack,
                                            x1,
                                            d as f64,
                                            b
                                        ));
                                    }
                                }

                                int_scale /= d as f64;
                            }

                            if int_scale != 1.0 {
                                self.scale_stored_row(row, int_scale, true);
                            }
                        }
                    }
                }
            } else {
                // inequality or ranged row, first store row positions
                self.store_row(row);

                if self.rowsize[row as usize]
                    == self.rowsize_integer[row as usize] + self.rowsize_impl_int[row as usize]
                {
                    let mut row_coefs: Vec<f64> = Vec::new();
                    let mut row_index: Vec<HighsInt> = Vec::new();

                    for nonz in self.get_stored_row().iter(self) {
                        debug_assert!(nonz.value() != 0.0);
                        row_coefs.push(nonz.value());
                        row_index.push(nonz.index());
                    }

                    let int_scale: f64 = HighsIntegers::integral_scale_vec(
                        &row_coefs,
                        if self.model.row_lower_[row as usize] == -kHighsInf {
                            self.primal_feastol
                        } else {
                            self.options.small_matrix_value
                        },
                        if self.model.row_upper_[row as usize] == kHighsInf {
                            self.primal_feastol
                        } else {
                            self.options.small_matrix_value
                        },
                    );

                    if int_scale != 0.0 {
                        let mut rounded_lhs = HighsCDouble::from(0.0);
                        let mut rounded_rhs = HighsCDouble::from(0.0);
                        let mut fraction_lhs = HighsCDouble::from(0.0);
                        let mut fraction_rhs = HighsCDouble::from(0.0);
                        let mut min_lhs_tightening = 0.0;
                        let mut min_rhs_tightening = 0.0;
                        let mut max_val = 0.0;
                        let mut lhs_tightened = false;
                        let mut rhs_tightened = false;
                        let mut is_infeasible = false;
                        if check_scale_row(
                            self,
                            row,
                            &mut row_coefs,
                            &row_index,
                            &mut rounded_lhs,
                            &mut rounded_rhs,
                            &mut fraction_lhs,
                            &mut fraction_rhs,
                            &mut min_lhs_tightening,
                            &mut min_rhs_tightening,
                            &mut max_val,
                            &mut lhs_tightened,
                            &mut rhs_tightened,
                            &mut is_infeasible,
                            int_scale,
                        ) {
                            // check for infeasibility
                            if is_infeasible {
                                return HPresolveResult::kPrimalInfeasible;
                            }
                            // only accept row whose sides were tightened
                            let ranged_or_equation_row = lhs_tightened && rhs_tightened;
                            if ranged_or_equation_row
                                || (lhs_tightened
                                    && self.model.row_upper_[row as usize] == kHighsInf)
                                || (rhs_tightened
                                    && self.model.row_lower_[row as usize] == -kHighsInf)
                            {
                                // check if constraint can be scaled to integral values
                                if !lam_scale_row_int_vals(
                                    self,
                                    row,
                                    rounded_lhs,
                                    rounded_rhs,
                                    int_scale,
                                    max_val,
                                    &row_coefs,
                                    &row_index,
                                ) {
                                    if ranged_or_equation_row {
                                        // ranged or equation row
                                        // scale value is large, just tighten the sides
                                        rounded_lhs /= int_scale;
                                        rounded_rhs /= int_scale;
                                        if rounded_rhs
                                            < self.model.row_upper_[row as usize]
                                                - self.primal_feastol
                                        {
                                            if self.options.h16_row_duals {
                                                // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp line 3132 (1.16.0)
                                                self.h16_change_row_upper(
                                                    row,
                                                    f64::from(rounded_rhs),
                                                    false,
                                                );
                                            } else {
                                                self.model.row_upper_[row as usize] =
                                                    f64::from(rounded_rhs);
                                                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 4560 at f2ae2b320 (1.15.1 4039)
                                                self.x1_row_modified(row);
                                            }
                                        }
                                        if rounded_lhs
                                            > self.model.row_lower_[row as usize]
                                                + self.primal_feastol
                                        {
                                            if self.options.h16_row_duals {
                                                // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp line 3134 (1.16.0)
                                                self.h16_change_row_lower(
                                                    row,
                                                    f64::from(rounded_lhs),
                                                    false,
                                                );
                                            } else {
                                                self.model.row_lower_[row as usize] =
                                                    f64::from(rounded_lhs);
                                                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 4562 at f2ae2b320 (1.15.1 4041)
                                                self.x1_row_modified(row);
                                            }
                                        }
                                    } else if (rhs_tightened
                                        && fraction_rhs < min_rhs_tightening - self.primal_feastol)
                                        || (lhs_tightened
                                            && fraction_lhs
                                                < min_lhs_tightening - self.primal_feastol)
                                    {
                                        // <= or >= inequality
                                        // scale value is large, so we scale back the altered
                                        // constraint the scaled back constraint must be stronger
                                        // than the original constraint for this to make sense with
                                        // is checked with the condition above
                                        lam_scale_row(
                                            self,
                                            row,
                                            rounded_lhs,
                                            rounded_rhs,
                                            int_scale,
                                            true,
                                            &row_coefs,
                                            &row_index,
                                        );
                                    }
                                }
                            }
                        }
                    } else if !self.is_ranged(row) {
                        // Chvatal-Gomory strengthening
                        // See section 3.4 "Chvatal-Gomory strengthening of inequalities",
                        // Achterberg et al., Presolve Reductions in Mixed Integer
                        // Programming, INFORMS Journal on Computing 32(2):473-506.
                        let max_dynamism: f64 = 1e5;
                        let mut rounded_row_coefs: Vec<f64> =
                            vec![0.0; self.rowsize[row as usize] as usize];
                        let mut scalars: BTreeSet<OrderedF64> = BTreeSet::new();
                        scalars.insert(OrderedF64(1.0));

                        // convert to >= inequality
                        // direction = -1: <= constraint, direction = 1: >= constraint
                        let direction: HighsInt =
                            if self.model.row_upper_[row as usize] != kHighsInf {
                                -1
                            } else {
                                1
                            };
                        // get rhs
                        let rhs_value = if direction < 0 {
                            -self.model.row_upper_[row as usize]
                        } else {
                            self.model.row_lower_[row as usize]
                        };
                        let mut rhs = HighsCDouble::from(rhs_value);

                        // initialise
                        let mut rounded_rhs = HighsCDouble::from(0.0);
                        let mut min_abs_coef = 0.0;
                        let mut max_abs_coef = 0.0;
                        let mut dynamism = 0.0;

                        // complement or shift variables to have lower bounds of zero
                        if complement_or_shift(
                            self,
                            &row_coefs,
                            &row_index,
                            direction,
                            &mut rhs,
                            &mut min_abs_coef,
                            &mut max_abs_coef,
                            &mut dynamism,
                        ) {
                            // identify scalars for row
                            if dynamism <= max_dynamism {
                                set_scalars(&mut scalars, min_abs_coef, max_abs_coef);
                            }
                            // find a scalar that produces improved coefficients
                            if row_can_be_tightened(
                                self,
                                &row_coefs,
                                &mut rounded_row_coefs,
                                &scalars,
                                &rhs,
                                &mut rounded_rhs,
                            ) {
                                // undo complementation and shifting
                                undo_complement_or_shift(
                                    self,
                                    &row_coefs,
                                    &row_index,
                                    &mut rounded_row_coefs,
                                    direction,
                                    &mut rounded_rhs,
                                );
                                // replace row by rounded one
                                update_row(
                                    self,
                                    row,
                                    &row_coefs,
                                    &row_index,
                                    &rounded_row_coefs,
                                    direction,
                                    &rounded_rhs,
                                );
                            }
                        }
                    }
                }

                // check for redundancy again
                hpresolve_checked_call!(check_row_redundant(
                    self,
                    postsolve_stack,
                    row,
                    logging_on
                ));
                if self.row_deleted[row as usize] != 0 {
                    return HPresolveResult::kOk;
                }

                if self.model.row_lower_[row as usize] == -kHighsInf {
                    // <= constraint: try to strengthen coefficients
                    let implied_row_upper = self.implied_row_bounds.get_sum_upper(row, 0.0);
                    if implied_row_upper != kHighsInf {
                        let mut rhs = HighsCDouble::from(self.model.row_upper_[row as usize]);
                        let max_abs = HighsCDouble::from(implied_row_upper)
                            - self.model.row_upper_[row as usize];
                        strengthen_coefs(self, row, &mut rhs, 1, max_abs);
                        if self.options.h16_row_duals {
                            // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp line 3358 (1.16.0)
                            self.h16_change_row_upper(row, f64::from(rhs), false);
                        } else {
                            self.model.row_upper_[row as usize] = f64::from(rhs);
                            // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 4786 at f2ae2b320 (1.15.1 4265)
                            self.x1_row_modified(row);
                        }
                    }
                }

                if self.model.row_upper_[row as usize] == kHighsInf {
                    // >= constraint: try to strengthen coefficients
                    let implied_row_lower = self.implied_row_bounds.get_sum_lower(row, 0.0);
                    if implied_row_lower != -kHighsInf {
                        let mut rhs = HighsCDouble::from(self.model.row_lower_[row as usize]);
                        let max_abs = self.model.row_lower_[row as usize]
                            - HighsCDouble::from(implied_row_lower);
                        strengthen_coefs(self, row, &mut rhs, -1, max_abs);
                        if self.options.h16_row_duals {
                            // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp line 3370 (1.16.0)
                            self.h16_change_row_lower(row, f64::from(rhs), false);
                        } else {
                            self.model.row_lower_[row as usize] = f64::from(rhs);
                            // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 4798 at f2ae2b320 (1.15.1 4277)
                            self.x1_row_modified(row);
                        }
                    }
                }
            }
        }

        // C++ `checkForcingRow` lambda.
        #[allow(clippy::too_many_arguments)]
        fn check_forcing_row(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
            row: HighsInt,
            direction: HighsInt,
            row_side: f64,
            implied_row_bound: f64,
            dynamism: HighsCDouble,
            row_type: RowType,
            logging_on: bool,
        ) -> HPresolveResult {
            // 1. direction =  1 (>=): forcing row if upper bound on constraint activity
            //                         is equal to row's lower bound
            // 2. direction = -1 (<=): forcing row if lower bound on constraint activity
            //                         is equal to row's upper bound
            // scale tolerance using dynamism (ratio between absolute largest and
            // smallest coefficients)
            // HiGHS 1.16 #3224 (h16_row_presolve): HPresolve.cpp lines 3389-3394 (1.16.0): the
            // tolerance is floored at std::numeric_limits<double>::epsilon()
            if (direction as f64) * row_side == -kHighsInf
                || (direction as f64) * implied_row_bound == kHighsInf
                || (HighsCDouble::from(row_side) - HighsCDouble::from(implied_row_bound)).abs()
                    > if this.options.h16_row_presolve {
                        HighsCDouble::from(cpp_max(
                            f64::from(this.primal_feastol / dynamism),
                            f64::EPSILON,
                        ))
                    } else {
                        this.primal_feastol / dynamism
                    }
            {
                return HPresolveResult::kOk;
            }

            // get stored row
            let row_vector = this.get_stored_row();

            let mut nfixings: HighsInt = 0;
            for nonzero in row_vector.clone().iter(&*this) {
                let idx = nonzero.index() as usize;
                if (direction as f64) * nonzero.value() > 0.0 {
                    if this.model.col_upper_[idx] <= this.impl_col_upper[idx] {
                        nfixings += 1;
                    }
                } else if this.model.col_lower_[idx] >= this.impl_col_lower[idx] {
                    nfixings += 1;
                }
            }

            if nfixings != this.rowsize[row as usize] {
                return HPresolveResult::kOk;
            }

            if logging_on {
                this.analysis_start_presolve_rule_log(kPresolveRuleForcingRow);
            }
            postsolve_stack.forcing_row(row, row_vector.clone().iter(&*this), row_side, row_type);
            // already mark the row as deleted, since otherwise it would be
            // registered as changed/singleton in the process of fixing and
            // removing the contained columns
            this.mark_row_deleted(row);
            let mut it = row_vector.clone();
            while let Some(nonzero) = it.next_nz(this) {
                if (direction as f64) * nonzero.value() > 0.0 {
                    // the upper bound of the column is as tight as the implied upper
                    // bound or comes from this row, which means it is not used in the
                    // rows implied bounds. Therefore we can fix the variable at its
                    // upper bound.
                    hpresolve_checked_call!(
                        this.fix_col_to_upper(postsolve_stack, nonzero.index())
                    );
                } else {
                    hpresolve_checked_call!(
                        this.fix_col_to_lower(postsolve_stack, nonzero.index())
                    );
                }
            }
            // now the row might be empty, but not necessarily because the implied
            // column bounds might be implied by other rows in which case we
            // cannot fix the column
            postsolve_stack.redundant_row(row);
            // Row removal accounted for above

            // if there are any new row singletons, also remove them immediately
            hpresolve_checked_call!(this.remove_row_singletons(postsolve_stack));
            this.analysis_.logging_on_ = logging_on;
            if logging_on {
                this.analysis_stop_presolve_rule_log(kPresolveRuleForcingRow);
            }
            this.check_limits(postsolve_stack)
        }

        if self.analysis_.allow_rule_[kPresolveRuleForcingRow as usize] {
            // Allow rule to consider forcing rows

            // store row and compute dynamism
            self.store_row(row);
            let dynamism = self.compute_dynamism(self.get_stored_row());

            // >= inequality
            let row_lower_v = self.model.row_lower_[row as usize];
            let activity_upper = self.implied_row_bounds.get_sum_upper_orig(row, 0.0);
            hpresolve_checked_call!(check_forcing_row(
                self,
                postsolve_stack,
                row,
                1,
                row_lower_v,
                activity_upper,
                dynamism,
                RowType::kGeq,
                logging_on
            ));
            if self.row_deleted[row as usize] != 0 {
                return HPresolveResult::kOk;
            }

            // <= inequality
            let row_upper_v = self.model.row_upper_[row as usize];
            let activity_lower = self.implied_row_bounds.get_sum_lower_orig(row, 0.0);
            hpresolve_checked_call!(check_forcing_row(
                self,
                postsolve_stack,
                row,
                -1,
                row_upper_v,
                activity_lower,
                dynamism,
                RowType::kLeq,
                logging_on
            ));
            if self.row_deleted[row as usize] != 0 {
                return HPresolveResult::kOk;
            }
        }

        // update implied bounds of all columns in given row
        hpresolve_checked_call!(self.update_col_implied_bounds_row(postsolve_stack, row));

        // extract variable bound constraints
        self.extract_var_bounds(row);

        self.check_limits(postsolve_stack)
    }
}

#[cfg(test)]
mod h16_row_presolve_tests {
    use super::*;
    use crate::lp_data::highs_run::Highs;
    use crate::prelude::*;

    /// Builds an `HPresolve` on `lp` with `h16_row_presolve` set to `on`.
    fn presolve_with(lp: HighsLp, on: bool, primal_feastol: f64) -> HPresolve {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.h16_row_presolve = on;
        options.primal_feasibility_tolerance = primal_feastol;
        let mut presolve = HPresolve::default();
        // `checkTimeLimit` asserts the timer (HPresolve.cpp:4258), and `rowPresolve` reaches it
        let timer = Some(crate::util::highs_timer::HighsTimer::new());
        assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, timer));
        presolve
            .analysis_
            .setup(&presolve.model, &presolve.options, true);
        presolve
    }

    /// The model of HiGHS 1.16 #3392 (check/TestPresolve.cpp, issue-3391): three free
    /// columns, rows `-x0 - x2 <= 0` and `-x0 + x1 - 2 x2 <= 0`.
    fn issue_3391_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.model_name_ = "issue-3391".to_string();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![0.0, 1.0, 0.0];
        lp.col_lower_ = vec![-kHighsInf; 3];
        lp.col_upper_ = vec![kHighsInf; 3];
        lp.row_lower_ = vec![-kHighsInf; 2];
        lp.row_upper_ = vec![0.0; 2];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 3, 5];
        lp.a_matrix_.index_ = vec![0, 1, 1, 0, 1];
        lp.a_matrix_.value_ = vec![-1.0, -1.0, 1.0, -1.0, -2.0];
        lp
    }

    fn issue_3391_row_presolve(on: bool) -> HPresolveResult {
        let mut presolve = presolve_with(issue_3391_lp(), on, 1e-7);
        // hand-set implied row dual bound (the initial sweep is not run here); its
        // source is column 0, which is not a singleton, so no column bound changes
        presolve.impl_row_dual_lower[1] = 1.0;
        presolve.row_dual_lower_source[1] = 0;
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(2, 3);
        presolve.row_presolve(&mut stack, 1)
    }

    /// HiGHS 1.16 #3392: a row whose implied row dual is positive while it has no lower
    /// bound is dual infeasible. Only the option-on path reports it.
    #[test]
    fn h16_row_presolve_dual_infeasible_row() {
        assert!(matches!(
            issue_3391_row_presolve(true),
            HPresolveResult::kDualInfeasible
        ));
    }

    /// Off: 1.15.1 `isImpliedEquationAtLower` asserts that this state cannot occur
    /// (1.15.1 HPresolve.cpp:268), so its outcome is only observable without debug asserts.
    #[test]
    #[cfg(not(debug_assertions))]
    fn h16_row_presolve_dual_infeasible_row_off() {
        assert!(!matches!(
            issue_3391_row_presolve(false),
            HPresolveResult::kDualInfeasible
        ));
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "is_lbnd_positive")]
    fn h16_row_presolve_dual_infeasible_row_off_trips_cpp_assert() {
        issue_3391_row_presolve(false);
    }

    /// The model of HiGHS 1.16 #3260 (check/TestMipSolver.cpp, TEST_CASE "pr-3260"):
    /// `x0 + x1 + x2 = -1` twice, `x0` integer fixed at 0, `x1` binary, `x2` in [-3, 0].
    fn pr_3260_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![0.0; 3];
        lp.col_lower_ = vec![0.0, 0.0, -3.0];
        lp.col_upper_ = vec![0.0, 1.0, 0.0];
        lp.row_lower_ = vec![-1.0, -1.0];
        lp.row_upper_ = vec![-1.0, -1.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 4, 6];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![1.0; 6];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
        ];
        lp
    }

    /// HiGHS 1.16 #3260: in the simple probing case on row 1 (binary `x1`), removing the
    /// fixed column `x0` changes the live implied bounds of `x2` from [-2, -1] to [-3, 0]
    /// before `x2` is substituted. With the option on the substitution uses the snapshot
    /// taken before the loop; the C++ test requires the MIP to solve to optimal with
    /// objective 0. With the option off the stale-bound substitution of 1.15.1 makes
    /// presolve declare the feasible MIP infeasible.
    #[test]
    fn h16_row_presolve_simple_probing_snapshot() {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        highs.options_.h16_row_presolve = true;
        assert_eq!(highs.pass_model(pr_3260_lp()), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        assert_eq!(highs.get_info().objective_function_value, 0.0);

        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(pr_3260_lp()), HighsStatus::kOk);
        assert_eq!(highs.presolve(), HighsStatus::kOk);
        assert!(matches!(
            highs.model_presolve_status_,
            HighsPresolveStatus::kInfeasible
        ));
    }

    /// HiGHS 1.16 #3392, TEST_CASE "issue-3391" end to end: IPM without crossover, the
    /// model status must be unbounded.
    #[test]
    fn h16_row_presolve_issue_3391_run() {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        highs.options_.h16_row_presolve = true;
        assert_eq!(
            highs.set_option_value_str("solver", "ipm"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("run_crossover", "off"),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(issue_3391_lp()), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kUnbounded);
    }

    /// HiGHS 1.16 #3224: a forcing-row tolerance floor of `std::numeric_limits<double>::epsilon()`.
    /// Row `x + 1e-300 y >= 1 - 2^-53` has dynamism 1e300, so `primal_feastol / dynamism` is
    /// about 1e-307 and the gap 2^-53 is above it. The row is forcing only when the floor
    /// (2.2e-16, above the gap) applies, which the option-on path uses.
    #[test]
    fn h16_row_presolve_forcing_row_tolerance_floor() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0; 2];
        lp.col_lower_ = vec![0.0; 2];
        lp.col_upper_ = vec![1.0; 2];
        lp.row_lower_ = vec![1.0 - f64::EPSILON / 2.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1e-300];

        let mut deleted = Vec::new();
        for on in [true, false] {
            let mut presolve = presolve_with(lp.clone(), on, 1e-7);
            let mut stack = HighsPostsolveStack::default();
            stack.initialize_index_maps(1, 2);
            let _ = presolve.row_presolve(&mut stack, 0);
            deleted.push(presolve.row_deleted[0] != 0);
        }
        assert!(deleted[0], "option on: forcing row must be removed");
        assert!(!deleted[1], "option off: 1.15.1 tolerance keeps the row");
    }

    /// HiGHS 1.16 #3224, singleton row: the bound tolerance is floored at epsilon. With
    /// `primal_feastol` 1e-18 and `x >= 1e-17`, the bound 0 -> 1e-17 is tightened only
    /// when the tolerance is the 1.15.1 value 1e-18.
    #[test]
    fn h16_row_presolve_singleton_tolerance_floor() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![1.0];
        lp.row_lower_ = vec![1e-17];
        lp.row_upper_ = vec![kHighsInf];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];

        let mut lower = Vec::new();
        for on in [true, false] {
            let mut presolve = presolve_with(lp.clone(), on, 1e-18);
            let mut stack = HighsPostsolveStack::default();
            stack.initialize_index_maps(1, 1);
            let _ = presolve.row_presolve(&mut stack, 0);
            lower.push(presolve.model.col_lower_[0]);
        }
        assert_eq!(
            lower[0], 0.0,
            "option on: epsilon floor rejects the tightening"
        );
        assert_eq!(lower[1], 1e-17, "option off: 1.15.1 tolerance accepts it");
    }

    /// rusolve `mip_presolve_exact_strengthen`: `min -10 x1 - 10 x2 + x3` subject to
    /// `0.784 x1 + 5.139 x2 - 1e8 x3 <= 0`, all binary. The optimum is x = (1, 1, 1),
    /// objective -19, exactly feasible. Coefficient strengthening replaces -1e8 by
    /// -(0.784 + 5.139); the 1.16.0 code adds a delta rounded to a double and gets
    /// -5.9229999929666519, which cuts off the optimum (HiGHS 1.15.1 and 1.16.0 return -9).
    fn bigm_strengthen_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![-10.0, -10.0, 1.0];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![1.0; 3];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![0.0];
        lp.integrality_ = vec![HighsVarType::kInteger; 3];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2, 3];
        lp.a_matrix_.index_ = vec![0, 0, 0];
        lp.a_matrix_.value_ = vec![0.784, 5.139, -1e8];
        lp
    }

    /// The x3 coefficient after one `row_presolve` of the big-M row.
    fn bigm_strengthened_coef(exact: bool) -> f64 {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_presolve_exact_strengthen = exact;
        let mut presolve = HPresolve::default();
        let timer = Some(crate::util::highs_timer::HighsTimer::new());
        assert!(presolve.ok_set_input(
            bigm_strengthen_lp(),
            &options,
            options.presolve_reduction_limit,
            timer
        ));
        presolve
            .analysis_
            .setup(&presolve.model, &presolve.options, true);
        // LP mode (no mipsolver) makes every column continuous; strengthening needs integers
        presolve.model.integrality_ = vec![HighsVarType::kInteger; 3];
        presolve.rowsize_integer[0] = 3;
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(1, 3);
        let _ = presolve.row_presolve(&mut stack, 0);
        let pos = presolve.find_nonzero(0, 2);
        assert!(pos != -1);
        presolve.a_value[pos as usize]
    }

    #[test]
    fn exact_strengthen_coef() {
        let intended = -(0.784f64 + 5.139);
        assert_eq!(bigm_strengthened_coef(true), intended);
        // the 1.16.0 path: the rounded delta cancels to a tighter coefficient
        let faithful = bigm_strengthened_coef(false);
        assert!(faithful > intended + 1e-9, "{faithful}");
    }

    fn bigm_strengthen_objective(exact: bool) -> f64 {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        highs.options_.mip_presolve_exact_strengthen = exact;
        assert_eq!(highs.pass_model(bigm_strengthen_lp()), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        highs.get_info().objective_function_value
    }

    #[test]
    fn exact_strengthen_keeps_optimum() {
        assert_eq!(bigm_strengthen_objective(true), -19.0);
        // HiGHS 1.15.1 and 1.16.0 behaviour with the lever off
        assert_eq!(bigm_strengthen_objective(false), -9.0);
    }

    #[test]
    fn exact_strengthen_recommended_on() {
        use crate::lp_data::highs_options::RECOMMENDED_MIP_ON;
        assert!(!HighsOptions::default().mip_presolve_exact_strengthen);
        assert!(HighsOptions::recommended().mip_presolve_exact_strengthen);
        assert!(RECOMMENDED_MIP_ON.contains(&"mip_presolve_exact_strengthen"));
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("mip_presolve_exact_strengthen", true),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_presolve_exact_strengthen"),
            Some(true)
        );
    }
}
