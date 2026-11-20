//! HiGHS source: presolve/HPresolve.cpp (lines 4670-4985)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! dualFixing. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::dualFixing (lines 4670-4985, 316 lines)
    ///
    /// Several `[&]` lambdas (`fixCols`, `handleSingleEquation`, `hasTighterBound`, `lockCallback`):
    /// see the module doc for capturing `this`.
    pub fn dual_fixing(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        // fix variables or tighten bounds using dual arguments
        // see section 4.4 "Dual fixing, substitution and bound strengthening",
        // Achterberg et al., Presolve Reductions in Mixed Integer Programming,
        // INFORMS Journal on Computing 32(2):473-506.
        debug_assert!(self.col_deleted[col as usize] == 0);

        // return if variable is already fixed
        if self.model.col_lower_[col as usize] == self.model.col_upper_[col as usize] {
            return HPresolveResult::kOk;
        }

        // struct for single equation handling
        struct equationNonZero {
            col: HighsInt,
            val: f64,
            mark: HighsInt,
        }

        // lambda that checks column for single equation handling
        fn check_column(
            this: &HPresolve,
            col: HighsInt,
            col_direction: HighsInt,
            row: HighsInt,
        ) -> bool {
            if (col_direction as f64) * this.model.col_cost_[col as usize] < 0.0 {
                return false;
            }
            for col_nz in this.get_column_vector(col).iter(this) {
                // skip equation we are inspecting; coefficient is positive (see sign
                // adjustment performed by caller)
                if col_nz.index() == row {
                    continue;
                }
                // skip redundant rows
                if this.is_redundant(col_nz.index()) {
                    continue;
                }
                // another ranged row blocks variable in both directions -> cannot be a
                // candidate for fixing
                if this.is_ranged(col_nz.index()) {
                    return false;
                }
                // compute row direction
                let row_direction: HighsInt = if this.model.row_lower_[col_nz.index() as usize]
                    == -kHighsInf
                    && this.model.row_upper_[col_nz.index() as usize] != kHighsInf
                {
                    1
                } else {
                    -1
                };
                // coefficient must have the correct sign
                if ((col_direction * row_direction) as f64) * col_nz.value() < 0.0 {
                    return false;
                }
            }
            true
        }

        // lambda that computes activities for single equation handling
        fn compute_activity(
            model: &HighsLp,
            col: HighsInt,
            val: f64,
            activity: &mut HighsCDouble,
            activity_finite: &mut bool,
            direction: HighsInt,
        ) {
            // direction =  1 -> update upper bound on activity (supremum)
            // direction = -1 -> update lower bound on activity (infimum)
            let bound = if (direction as f64) * val > 0.0 {
                model.col_upper_[col as usize]
            } else {
                model.col_lower_[col as usize]
            };
            *activity_finite = *activity_finite && bound.abs() != kHighsInf;
            if *activity_finite {
                *activity += HighsCDouble::from(val) * bound;
            }
        }

        // lambda for fixing variables
        fn fix_cols(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
            equation_row_vector: &[equationNonZero],
            direction: HighsInt,
        ) -> HPresolveResult {
            for row_nz in equation_row_vector {
                if direction * row_nz.mark >= 0 {
                    continue;
                }
                if (direction as f64) * row_nz.val > 0.0 {
                    hpresolve_checked_call!(this.fix_col_to_upper(postsolve_stack, row_nz.col));
                } else {
                    hpresolve_checked_call!(this.fix_col_to_lower(postsolve_stack, row_nz.col));
                }
            }
            HPresolveResult::kOk
        }

        // lambda for handling single equations
        fn handle_single_equation(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
            row: HighsInt,
        ) -> HPresolveResult {
            debug_assert!(this.is_equation(row));
            let mut equation_row_vector: Vec<equationNonZero> =
                Vec::with_capacity(this.rowsize[row as usize] as usize);
            let mut num_s_plus: HighsInt = 0;
            let mut num_s_minus: HighsInt = 0;
            for row_nz in this.get_row_vector(row).iter(&*this) {
                // flip column direction if coefficient is negative
                // (equivalent to negating the column's coefficients)
                let col_direction: HighsInt = (1.0f64).copysign(row_nz.value()) as HighsInt;
                let mut mark: HighsInt = 0;
                if check_column(&*this, row_nz.index(), col_direction, row) {
                    // store in S_+
                    num_s_plus += 1;
                    mark = 1;
                } else if check_column(&*this, row_nz.index(), -col_direction, row) {
                    // store in S_-
                    num_s_minus += 1;
                    mark = -1;
                }
                // store non-zero
                equation_row_vector.push(equationNonZero {
                    col: row_nz.index(),
                    val: row_nz.value(),
                    mark,
                });
            }
            // return if both sets are empty
            if num_s_plus == 0 && num_s_minus == 0 {
                return HPresolveResult::kOk;
            }
            // compute activities
            let mut activity_t_plus = HighsCDouble::from(0.0);
            let mut activity_t_minus = HighsCDouble::from(0.0);
            let mut activity_sc_plus = HighsCDouble::from(0.0);
            let mut activity_sc_minus = HighsCDouble::from(0.0);
            let mut activity_t_plus_finite = true;
            let mut activity_t_minus_finite = true;
            let mut activity_sc_plus_finite = true;
            let mut activity_sc_minus_finite = true;
            for row_nz in equation_row_vector.iter() {
                if row_nz.mark <= 0
                    || this.model.integrality_[row_nz.col as usize] != HighsVarType::kContinuous
                {
                    // T_+
                    compute_activity(
                        &this.model,
                        row_nz.col,
                        row_nz.val,
                        &mut activity_t_plus,
                        &mut activity_t_plus_finite,
                        1,
                    );
                } else {
                    // S^C_+
                    compute_activity(
                        &this.model,
                        row_nz.col,
                        row_nz.val,
                        &mut activity_sc_plus,
                        &mut activity_sc_plus_finite,
                        -1,
                    );
                }
                if row_nz.mark >= 0
                    || this.model.integrality_[row_nz.col as usize] != HighsVarType::kContinuous
                {
                    // T_-
                    compute_activity(
                        &this.model,
                        row_nz.col,
                        row_nz.val,
                        &mut activity_t_minus,
                        &mut activity_t_minus_finite,
                        -1,
                    );
                } else {
                    // S^C_-
                    compute_activity(
                        &this.model,
                        row_nz.col,
                        row_nz.val,
                        &mut activity_sc_minus,
                        &mut activity_sc_minus_finite,
                        1,
                    );
                }
                // break if activities are not finite
                if (num_s_minus == 0 || !activity_t_plus_finite || !activity_sc_plus_finite)
                    && (num_s_plus == 0 || !activity_t_minus_finite || !activity_sc_minus_finite)
                {
                    break;
                }
            }
            // fix variables
            if num_s_minus > 0
                && activity_t_plus_finite
                && activity_sc_plus_finite
                && activity_t_plus + activity_sc_plus
                    <= this.model.row_lower_[row as usize] + this.primal_feastol
            {
                // fix all variables in S_-
                hpresolve_checked_call!(fix_cols(this, postsolve_stack, &equation_row_vector, 1));
            } else if num_s_plus > 0
                && activity_t_minus_finite
                && activity_sc_minus_finite
                && activity_t_minus + activity_sc_minus
                    >= this.model.row_lower_[row as usize] - this.primal_feastol
            {
                // fix all variables in S_+
                hpresolve_checked_call!(fix_cols(this, postsolve_stack, &equation_row_vector, -1));
            }
            HPresolveResult::kOk
        }

        // lambda for variable substitution
        /// HiGHS latest (#3404, 3fa5161b1): HPresolveDualFixing.cpp::HPresolve::dualFixing substituteCol lambda (lines 31-58)
        fn substitute_col(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
            col: HighsInt,
            row: HighsInt,
            direction: HighsInt,
            col_bound: f64,
            other_col_bound: f64,
        ) -> HPresolveResult {
            let x1_dsub = this.options.mip_dual_subst_limits;
            if x1_dsub {
                debug_assert!(
                    this.model.col_lower_[col as usize] != -kHighsInf
                        && this.model.col_upper_[col as usize] != kHighsInf
                );
                const MAX_SUBST_ROW_SIZE: HighsInt = 1000;
                if this.rowsize_integer[row as usize] == 0
                    || this.rowsize[row as usize] > MAX_SUBST_ROW_SIZE
                {
                    return HPresolveResult::kOk;
                }
                let lhs_finite = this.model.row_lower_[row as usize] != -kHighsInf;
                let rhs_finite = this.model.row_upper_[row as usize] != kHighsInf;
                if lhs_finite && rhs_finite {
                    // Row cannot be redundant if the column to be substituted can
                    // account for all row activity
                    let pos = this.find_nonzero(row, col);
                    debug_assert!(pos != -1);
                    let col_coef = this.a_value[pos as usize];
                    // double * HighsCDouble as in C++ (friend operator*, i.e. range * |coef|)
                    let activity_range = col_coef.abs()
                        * (HighsCDouble::from(this.model.col_upper_[col as usize])
                            - this.model.col_lower_[col as usize]);
                    if activity_range
                        > HighsCDouble::from(this.model.row_upper_[row as usize])
                            - this.model.row_lower_[row as usize]
                            + 2.0 * this.primal_feastol
                    {
                        return HPresolveResult::kOk;
                    }
                }
            }
            // check lhs and rhs for finiteness
            let lhs_finite = this.model.row_lower_[row as usize] != -kHighsInf;
            let rhs_finite = this.model.row_upper_[row as usize] != kHighsInf;

            // use storeRow and getStoredRow since getRowVector's rowroot[row] would be
            // overwritten by subsequent findNonZero calls, which would produce
            // undefined behavior
            this.store_row(row);
            let mut it = this.get_stored_row();
            while let Some(nz) = it.next_nz(&*this) {
                let idx = nz.index();
                // skip column index that was passed to this lambda
                if idx == col {
                    continue;
                }

                // only consider non-fixed binary variables
                if this.model.integrality_[idx as usize] != HighsVarType::kInteger
                    || this.model.col_lower_[idx as usize] != 0.0
                    || this.model.col_upper_[idx as usize] != 1.0
                {
                    continue;
                }

                let val = nz.value();

                // skip binary variable if setting it to its lower bound does not make the
                // row redundant
                let row_upper = this.model.row_upper_[row as usize];
                let row_lower = this.model.row_lower_[row as usize];
                if (rhs_finite
                    && this.implied_row_bounds.get_residual_sum_upper_orig(
                        &col_bound_arrays!(this),
                        row,
                        idx,
                        val,
                        -1,
                        -kHighsInf,
                        -kHighsInf,
                    ) > row_upper + this.primal_feastol)
                    || (lhs_finite
                        && this.implied_row_bounds.get_residual_sum_lower_orig(
                            &col_bound_arrays!(this),
                            row,
                            idx,
                            val,
                            -1,
                            -kHighsInf,
                            -kHighsInf,
                        ) < row_lower - this.primal_feastol)
                {
                    continue;
                }

                // now compute the implied lower bound (direction = 1) or implied upper
                // bound (direction = -1) provided that the binary variable is set to its
                // upper bound.
                let col_upper_idx = this.model.col_upper_[idx as usize];
                let best_bound = if direction > 0 {
                    this.compute_implied_lower_bound(col, idx, col_upper_idx, 0)
                } else {
                    -this.compute_implied_upper_bound(col, idx, col_upper_idx, 0)
                };

                // check if lower / upper bound is implied
                if best_bound >= (direction as f64) * col_bound - this.primal_feastol {
                    // substitute variable
                    let offset = other_col_bound;
                    let scale = col_bound - other_col_bound;
                    let lower = this.model.col_lower_[col as usize];
                    let upper = this.model.col_upper_[col as usize];
                    postsolve_stack.doubleton_equation(
                        -1,
                        col,
                        idx,
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
                    // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 5333 at 7b890c255
                    if this.options.mip_presolve_clique_table {
                        hpresolve_checked_call!(this.substitute_x1b(col, idx, offset, scale, -1));
                    } else {
                        this.mark_col_deleted(col);
                        this.substitute_col(col, idx, offset, scale);
                    }
                    hpresolve_checked_call!(this.check_limits(postsolve_stack));
                    break;
                }
            }
            HPresolveResult::kOk
        }

        // lambda for computing tighter bounds
        fn has_tighter_bound(
            this: &mut HPresolve,
            col: HighsInt,
            direction: HighsInt,
            current_bound: f64,
            new_bound: &mut f64,
        ) -> bool {
            // return if objective coefficient has wrong sign
            if (direction as f64) * this.model.col_cost_[col as usize] < 0.0 {
                return false;
            }

            // do not accept huge bounds
            let huge_bound = this.primal_feastol / kHighsTiny;

            // compute worst-case bounds (direction = 1: lower bound,
            // direction = -1: upper bound)
            *new_bound = if direction > 0 {
                this.compute_worst_case_lower_bound(col, -1, kHighsInf, 0)
            } else {
                -this.compute_worst_case_upper_bound(col, -1, kHighsInf, 0)
            };

            // return if bound is too large
            if *new_bound != -kHighsInf
                && (*new_bound >= (direction as f64) * current_bound - this.primal_feastol
                    || (*new_bound).abs() > huge_bound)
            {
                return false;
            }

            // round up to make sure that all rows are redundant
            // HiGHS 1.16 #3411 (h16_int_bound_rounding): HPresolveDualFixing.cpp::dualFixing
            // (3411.diff 195-209): 1.16.0 removes this rounding
            if !this.options.h16_int_bound_rounding
                && this.model.integrality_[col as usize] != HighsVarType::kContinuous
            {
                *new_bound = (*new_bound - this.primal_feastol).ceil();
            }

            // flip sign
            *new_bound *= direction as f64;
            true
        }

        // lambda callback for lock computation
        let mut num_down_locks: HighsInt = 0;
        let mut num_up_locks: HighsInt = 0;
        let mut down_lock_row: HighsInt = -1;
        let mut up_lock_row: HighsInt = -1;
        let mut lock_callback = |row: HighsInt, has_down_lock: bool, has_up_lock: bool| -> bool {
            // count locks and remember row index
            if has_up_lock {
                num_up_locks += 1;
                up_lock_row = row;
            }
            if has_down_lock {
                num_down_locks += 1;
                down_lock_row = row;
            }
            // stop early if there are locks in both directions, since the variable
            // cannot be fixed in this case.
            num_down_locks > 1 && num_up_locks > 1
        };
        // compute locks
        self.compute_locks(col, true, &mut lock_callback);

        // check if variable can be fixed
        if num_down_locks == 0 || num_up_locks == 0 {
            if self.options.h16_zero_cost_singleton {
                // HiGHS 1.16 #2962 (h16_zero_cost_singleton): HPresolveDualFixing.cpp::dualFixing
                // (1.16.0 lines 288-298, 2962.diff hunk at 5334-5360): if |cost| <= dual
                // tolerance the column is dual feasible at any value, so only fix when the
                // bound it is fixed to is finite
                if num_down_locks == 0 {
                    if self.model.col_cost_[col as usize] > self.options.dual_feasibility_tolerance
                        || self.model.col_lower_[col as usize] > -kHighsInf
                    {
                        hpresolve_checked_call!(self.fix_col_to_lower(postsolve_stack, col));
                    }
                } else {
                    debug_assert!(num_up_locks == 0);
                    if self.model.col_cost_[col as usize] < -self.options.dual_feasibility_tolerance
                        || self.model.col_upper_[col as usize] < kHighsInf
                    {
                        hpresolve_checked_call!(self.fix_col_to_upper(postsolve_stack, col));
                    }
                }
            } else {
                // fix variable
                if num_down_locks == 0 {
                    hpresolve_checked_call!(self.fix_col_to_lower(postsolve_stack, col));
                } else {
                    hpresolve_checked_call!(self.fix_col_to_upper(postsolve_stack, col));
                }
            }
        } else {
            let has_single_down_lock = num_down_locks == 1 && down_lock_row != -1;
            let has_single_up_lock = num_up_locks == 1 && up_lock_row != -1;
            if has_single_down_lock || has_single_up_lock {
                let equation_row: HighsInt =
                    if has_single_down_lock && self.is_equation(down_lock_row) {
                        down_lock_row
                    } else if has_single_up_lock && self.is_equation(up_lock_row) {
                        up_lock_row
                    } else {
                        -1
                    };
                if equation_row != -1 && self.single_equation_checked[equation_row as usize] == 0 {
                    // see section 6.1 "Extension of dual fixing for single equations",
                    // Achterberg et al., Presolve Reductions in Mixed Integer
                    // Programming, INFORMS Journal on Computing 32(2):473-506.
                    hpresolve_checked_call!(handle_single_equation(
                        self,
                        postsolve_stack,
                        equation_row
                    ));
                    self.single_equation_checked[equation_row as usize] = 1;
                    if self.col_deleted[col as usize] != 0 {
                        return HPresolveResult::kOk;
                    }
                } else if self.mipsolver.is_some()
                    && self.model.col_lower_[col as usize] != -kHighsInf
                    && self.model.col_upper_[col as usize] != kHighsInf
                {
                    // try substitution
                    if has_single_down_lock {
                        let col_upper = self.model.col_upper_[col as usize];
                        let col_lower = self.model.col_lower_[col as usize];
                        hpresolve_checked_call!(substitute_col(
                            self,
                            postsolve_stack,
                            col,
                            down_lock_row,
                            1,
                            col_upper,
                            col_lower
                        ));
                        if self.col_deleted[col as usize] != 0 {
                            return HPresolveResult::kOk;
                        }
                    }
                    if has_single_up_lock {
                        let col_upper = self.model.col_upper_[col as usize];
                        let col_lower = self.model.col_lower_[col as usize];
                        hpresolve_checked_call!(substitute_col(
                            self,
                            postsolve_stack,
                            col,
                            up_lock_row,
                            -1,
                            col_lower,
                            col_upper
                        ));
                        if self.col_deleted[col as usize] != 0 {
                            return HPresolveResult::kOk;
                        }
                    }
                }
            }
            // try to strengthen bounds
            let mut new_bound: f64 = 0.0;
            let col_upper = self.model.col_upper_[col as usize];
            if has_tighter_bound(self, col, 1, col_upper, &mut new_bound) {
                // do not make bounds inconsistent
                new_bound = cpp_max(new_bound, self.model.col_lower_[col as usize]);
                // update upper bound
                if new_bound < self.model.col_upper_[col as usize] - self.primal_feastol {
                    if new_bound == self.model.col_lower_[col as usize] {
                        hpresolve_checked_call!(self.fix_col_to_lower(postsolve_stack, col));
                    } else if self.model.integrality_[col as usize] != HighsVarType::kContinuous {
                        hpresolve_checked_call!(self.change_col_upper(col, new_bound));
                    }
                }
            } else {
                let col_lower = self.model.col_lower_[col as usize];
                if has_tighter_bound(self, col, -1, col_lower, &mut new_bound) {
                    // do not make bounds inconsistent
                    new_bound = cpp_min(new_bound, self.model.col_upper_[col as usize]);
                    if new_bound > self.model.col_lower_[col as usize] + self.primal_feastol {
                        // update lower bound
                        if new_bound == self.model.col_upper_[col as usize] {
                            hpresolve_checked_call!(self.fix_col_to_upper(postsolve_stack, col));
                        } else if self.model.integrality_[col as usize] != HighsVarType::kContinuous
                        {
                            hpresolve_checked_call!(self.change_col_lower(col, new_bound));
                        }
                    }
                }
            }
        }
        HPresolveResult::kOk
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Runs dual fixing on col 0 (`x`) of one ranged row `20 <= x + 10 y <= 25`, `x` continuous
    /// in [0, 10], `y` binary, with the option set to `on`; returns whether col 0 was deleted.
    fn ranged_lock_mip(on: bool) -> bool {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 10.0];
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 1.0];
        lp.integrality_ = vec![HighsVarType::kContinuous, HighsVarType::kInteger];
        lp.row_lower_ = vec![20.0];
        lp.row_upper_ = vec![25.0];

        let mut options = HighsOptions::default();
        options.mip_dual_subst_limits = on;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input_mip(Box::new(ms), options.presolve_reduction_limit));
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(1, 2);
        let _ = presolve.dual_fixing(&mut stack, 0);
        presolve.col_deleted[0] != 0
    }

    #[test]
    fn x1_dsub_range_early_stop() {
        // Activity range of x (1 * (10 - 0) = 10) exceeds the row range (25 - 20) plus
        // 2 * feastol, so with the option on substitute_col returns before its binary loop.
        // On the 1.15.1 path that loop substitutes nothing here either (the residual
        // activity check skips y), so the option must not change the outcome.
        assert_eq!(ranged_lock_mip(true), ranged_lock_mip(false));
    }

    /// Row `x - 10 y + sum_i z_i >= 0` with `x` continuous in [0, 10] and cost 1 (its only
    /// down lock is the row), `y` binary and `num_fixed` integer columns `z_i` fixed at 0.
    /// Runs dual fixing on `x` with the option set to `on`; returns whether `x` was deleted.
    fn large_row_subst_mip(num_fixed: HighsInt, on: bool) -> bool {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2 + num_fixed;
        lp.num_row_ = 1;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, -10.0];
        lp.col_cost_ = vec![1.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 1.0];
        lp.integrality_ = vec![HighsVarType::kContinuous, HighsVarType::kInteger];
        for _ in 0..num_fixed {
            lp.col_cost_.push(0.0);
            lp.col_lower_.push(0.0);
            lp.col_upper_.push(0.0);
            lp.integrality_.push(HighsVarType::kInteger);
            lp.a_matrix_.index_.push(0);
            lp.a_matrix_.value_.push(1.0);
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        lp.row_lower_ = vec![0.0];
        lp.row_upper_ = vec![kHighsInf];

        let mut options = HighsOptions::default();
        options.mip_dual_subst_limits = on;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input_mip(Box::new(ms), options.presolve_reduction_limit));
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(1, 2 + num_fixed);
        let _ = presolve.dual_fixing(&mut stack, 0);
        presolve.col_deleted[0] != 0
    }

    #[test]
    fn x1_dsub_subst_row_too_large_skips() {
        // x = 10 y is substituted when the row is small, with the option on or off.
        assert!(large_row_subst_mip(10, false));
        assert!(large_row_subst_mip(10, true));
        // Row of 1001 entries (> 1000): substitution is skipped only with the option on.
        assert!(large_row_subst_mip(999, false));
        assert!(!large_row_subst_mip(999, true));
    }
}
