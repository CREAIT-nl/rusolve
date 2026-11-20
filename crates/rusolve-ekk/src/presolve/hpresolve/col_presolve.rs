//! HiGHS source: presolve/HPresolve.cpp (lines 4380-4668)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Empty columns, column presolve, dominated columns, locks. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::emptyCol (lines 4380-4407, 28 lines)
    pub fn empty_col(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleEmptyCol);
        }
        if (self.model.col_cost_[col as usize] > 0.0
            && self.model.col_lower_[col as usize] == -kHighsInf)
            || (self.model.col_cost_[col as usize] < 0.0
                && self.model.col_upper_[col as usize] == kHighsInf)
        {
            if self.model.col_cost_[col as usize].abs() <= self.options.dual_feasibility_tolerance {
                self.model.col_cost_[col as usize] = 0.0;
            } else {
                return HPresolveResult::kDualInfeasible;
            }
        }

        if self.model.col_cost_[col as usize] > 0.0 {
            hpresolve_checked_call!(self.fix_col_to_lower(postsolve_stack, col));
        } else if self.model.col_cost_[col as usize] < 0.0
            || self.model.col_upper_[col as usize].abs() < self.model.col_lower_[col as usize].abs()
        {
            hpresolve_checked_call!(self.fix_col_to_upper(postsolve_stack, col));
        } else if self.model.col_lower_[col as usize] != -kHighsInf {
            hpresolve_checked_call!(self.fix_col_to_lower(postsolve_stack, col));
        } else {
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 4953 at 7b890c255
            if self.options.mip_presolve_clique_table {
                hpresolve_checked_call!(self.fix_col_to_zero_x1b(postsolve_stack, col));
            } else {
                self.fix_col_to_zero(postsolve_stack, col);
            }
        }

        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleEmptyCol);
        }
        self.check_limits(postsolve_stack)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::colPresolve (lines 4409-4512, 104 lines)
    ///
    /// The `modifyImpliedRowDualBound` lambda sits in a `mipsolver != nullptr` branch.
    pub fn col_presolve(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        debug_assert!(self.col_deleted[col as usize] == 0);

        // check bounds
        let mut is_fixed = false;
        hpresolve_checked_call!(self.check_col_bounds(col, Some(&mut is_fixed)));
        if is_fixed {
            // remove fixed column
            postsolve_stack.removed_fixed_col(
                col,
                self.model.col_lower_[col as usize],
                self.model.col_cost_[col as usize],
                self.get_column_vector(col).iter(&*self),
            );
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 4975 at 7b890c255
            if self.options.mip_presolve_clique_table {
                hpresolve_checked_call!(self.remove_fixed_col_x1b(col));
            } else {
                self.remove_fixed_col(col);
            }
            return self.check_limits(postsolve_stack);
        }

        match self.colsize[col as usize] {
            0 => return self.empty_col(postsolve_stack, col),
            1 => return self.singleton_col(postsolve_stack, col),
            _ => {}
        }

        // detect strong / weak domination
        // HiGHS 1.16 #3348 (h16_dual_postsolve_rules): HPresolve.cpp::colPresolve (1.16.0
        // lines 3503-3634, the dominated column call-site guard)
        if !self.options.h16_dual_postsolve_rules
            || self.analysis_.allow_rule_[kPresolveRuleDominatedCol as usize]
        {
            hpresolve_checked_call!(self.detect_dominated_col(postsolve_stack, col, true));
        }
        if self.col_deleted[col as usize] != 0 {
            return HPresolveResult::kOk;
        }

        // column is not (weakly) dominated

        // integer columns cannot be used to tighten bounds on dual multipliers
        if self.mipsolver.is_some() {
            // lambda for changing implied row dual
            fn modify_implied_row_dual_bound(
                this: &mut HPresolve,
                col: HighsInt,
                row: HighsInt,
                direction: HighsInt,
                is_bound_implied: bool,
                num_inf: HighsInt,
            ) {
                if is_bound_implied
                    && row != -1
                    && num_inf == 1
                    && (direction as f64) * this.model.col_cost_[col as usize] >= 0.0
                    && !this.is_ranged(row)
                {
                    let nz_pos = this.find_nonzero(row, col);
                    let a_nz = this.a_value[nz_pos as usize];

                    if this.model.integrality_[col as usize] != HighsVarType::kInteger
                        || (this.rowsize_integer[row as usize] == this.rowsize[row as usize]
                            && this.row_coefficients_integral(row, 1.0 / a_nz))
                    {
                        if (direction as f64) * a_nz > 0.0 {
                            this.change_impl_row_dual_lower(row, 0.0, col);
                        } else {
                            this.change_impl_row_dual_upper(row, 0.0, col);
                        }
                    }
                }
            }

            // if there is an infinite or redundant column lower bound, the reduced
            // cost of the column must not be positive (i.e. <= 0; the upper bound on
            // the reduced cost is zero).
            let col_lower_row = self.col_lower_source[col as usize];
            let lower_implied = self.is_lower_implied(col);
            let num_inf_upper = self.implied_dual_row_bounds.get_num_inf_sum_upper_orig(col);
            modify_implied_row_dual_bound(
                self,
                col,
                col_lower_row,
                1,
                lower_implied,
                num_inf_upper,
            );

            // if there is an infinite or redundant column upper bound, the reduced
            // cost of the column must not be negative (i.e. >= 0; the lower bound on
            // the reduced cost is zero).
            let col_upper_row = self.col_upper_source[col as usize];
            let upper_implied = self.is_upper_implied(col);
            let num_inf_lower = self.implied_dual_row_bounds.get_num_inf_sum_lower_orig(col);
            modify_implied_row_dual_bound(
                self,
                col,
                col_upper_row,
                -1,
                upper_implied,
                num_inf_lower,
            );

            // check if variable is implied integer
            hpresolve_checked_call!(self.convert_implied_integer(col, -1, false).result());

            // HiGHS #2990 (mip_machine_sched_sepa): HPresolve.cpp lines 4480-4486 at dc83470a9
            // shift "binary" variables to have a lower bound of zero
            let col_lower = self.model.col_lower_[col as usize];
            let col_upper = self.model.col_upper_[col as usize];
            if self.model.integrality_[col as usize] != HighsVarType::kContinuous
                && col_lower != 0.0
                && (col_lower != -kHighsInf || col_upper != kHighsInf)
                && col_upper - col_lower > 0.5
                && (!self.options.mip_machine_sched_sepa || col_upper - col_lower < 1.5)
            {
                // substitute with the bound that is smaller in magnitude and only
                // substitute if bound is not large for an integer
                if col_upper.abs() > col_lower.abs() {
                    if col_lower.abs() < 1000.5 {
                        hpresolve_checked_call!(self.transform_column(
                            postsolve_stack,
                            col,
                            1.0,
                            col_lower
                        ));
                    }
                } else if col_upper.abs() < 1000.5 {
                    hpresolve_checked_call!(self.transform_column(
                        postsolve_stack,
                        col,
                        -1.0,
                        col_upper
                    ));
                }
            }
        }

        // dual fixing
        // HiGHS 1.16 #3348 (h16_dual_postsolve_rules): HPresolve.cpp::colPresolve (1.16.0
        // lines 3503-3634, the dual fixing call-site guard)
        if !self.options.h16_dual_postsolve_rules
            || self.analysis_.allow_rule_[kPresolveRuleDualFixing as usize]
        {
            hpresolve_checked_call!(self.dual_fixing(postsolve_stack, col));
        }
        if self.col_deleted[col as usize] != 0 {
            return HPresolveResult::kOk;
        }

        // singleton column stuffing
        hpresolve_checked_call!(self.singleton_col_stuffing(postsolve_stack, col));
        if self.col_deleted[col as usize] != 0 {
            return HPresolveResult::kOk;
        }

        // update dual implied bounds of all rows in given column
        if self.model.integrality_[col as usize] != HighsVarType::kInteger {
            self.update_row_dual_implied_bounds_col(col);
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::detectDominatedCol (lines 4514-4640, 127 lines)
    ///
    /// C++ default arguments: `handleSingletonRows = true`.
    pub fn detect_dominated_col(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
        handle_singleton_rows: bool,
    ) -> HPresolveResult {
        // HiGHS 1.16 #3348 (h16_dual_postsolve_rules): HPresolve.cpp::detectDominatedCol (1.16.0
        // lines 3636-3771, the allow_rule_ assert)
        if self.options.h16_dual_postsolve_rules {
            debug_assert!(self.analysis_.allow_rule_[kPresolveRuleDominatedCol as usize]);
        }
        debug_assert!(self.col_deleted[col as usize] == 0);

        // get bounds on column dual
        let col_dual_upper = -self
            .implied_dual_row_bounds
            .get_sum_lower(col, -self.model.col_cost_[col as usize]);
        let col_dual_lower = -self
            .implied_dual_row_bounds
            .get_sum_upper(col, -self.model.col_cost_[col as usize]);

        // compute dynamism
        let dynamism: HighsCDouble = self.compute_dynamism(self.get_column_vector(col));

        let logging_on = self.analysis_.logging_on_;

        // C++ lambda `dominatedCol`: captures `this`, `logging_on`, `handleSingletonRows`.
        let dominated_col = |this: &mut HPresolve,
                             postsolve_stack: &mut HighsPostsolveStack,
                             col: HighsInt,
                             dual_bound: f64,
                             bound: f64,
                             direction: HighsInt|
         -> HPresolveResult {
            // column is (strongly) dominated if the bounds on the column dual satisfy:
            // 1. lower bound >  dual feasibility tolerance (direction =  1) or
            // 2. upper bound < -dual feasibility tolerance (direction = -1).
            if (direction as f64) * dual_bound <= this.options.dual_feasibility_tolerance {
                return HPresolveResult::kOk;
            }
            // cannot fix to +-infinity -> infeasible
            if (direction as f64) * bound == -kHighsInf {
                return HPresolveResult::kDualInfeasible;
            }
            if logging_on {
                this.analysis_start_presolve_rule_log(kPresolveRuleDominatedCol);
            }
            // fix variable
            if direction > 0 {
                hpresolve_checked_call!(this.fix_col_to_lower(postsolve_stack, col));
            } else {
                hpresolve_checked_call!(this.fix_col_to_upper(postsolve_stack, col));
            }
            this.analysis_.logging_on_ = logging_on;
            if logging_on {
                this.analysis_stop_presolve_rule_log(kPresolveRuleDominatedCol);
            }
            // handle row singletons (if requested)
            if handle_singleton_rows {
                hpresolve_checked_call!(this.remove_row_singletons(postsolve_stack));
            }
            this.check_limits(postsolve_stack)
        };

        // C++ lambda `weaklyDominatedCol`: captures `this`, `logging_on`, `handleSingletonRows`.
        let weakly_dominated_col = |this: &mut HPresolve,
                                    postsolve_stack: &mut HighsPostsolveStack,
                                    col: HighsInt,
                                    dual_bound: f64,
                                    bound: f64,
                                    other_bound: f64,
                                    dynamism: &HighsCDouble,
                                    direction: HighsInt|
         -> HPresolveResult {
            // column is weakly dominated if the bounds on the column dual satisfy:
            // 1. lower bound >= -dual feasibility tolerance (direction =  1) or
            // 2. upper bound <=  dual feasibility tolerance (direction = -1).
            if (direction as f64) * dual_bound < -this.options.dual_feasibility_tolerance {
                return HPresolveResult::kOk;
            }
            if (direction as f64) * bound != -kHighsInf {
                if logging_on {
                    this.analysis_start_presolve_rule_log(kPresolveRuleDominatedCol);
                }
                // fix variable
                if direction > 0 {
                    hpresolve_checked_call!(this.fix_col_to_lower(postsolve_stack, col));
                } else {
                    hpresolve_checked_call!(this.fix_col_to_upper(postsolve_stack, col));
                }
                this.analysis_.logging_on_ = logging_on;
                if logging_on {
                    this.analysis_stop_presolve_rule_log(kPresolveRuleDominatedCol);
                }
                // handle row singletons (if requested)
                if handle_singleton_rows {
                    hpresolve_checked_call!(this.remove_row_singletons(postsolve_stack));
                }
                return this.check_limits(postsolve_stack);
            } else if this.analysis_.allow_rule_[kPresolveRuleForcingCol as usize] {
                // check for forcing column (see Andersen and Andersen, Presolving in
                // linear programming. Math. Program. 71, 221-245, 1995).
                // the column's lower bound is infinite (direction = 1) or its upper
                // bound is infinite (direction = -1).
                // now get lower bound (direction = 1) or upper bound (direction = -1) on
                // column dual using the original bounds on the row duals.
                let bound_on_col_dual = if direction > 0 {
                    -this
                        .implied_dual_row_bounds
                        .get_sum_upper_orig(col, -this.model.col_cost_[col as usize])
                } else {
                    -this
                        .implied_dual_row_bounds
                        .get_sum_lower_orig(col, -this.model.col_cost_[col as usize])
                };
                // HiGHS 1.16 #3224 (h16_row_presolve): HPresolve.cpp::detectDominatedCol (1.16.0
                // lines 3636-3771; hunk reference/h16-prs/3224.diff 26-36): the dual
                // tolerance is floored at numeric_limits<double>::epsilon()
                let dual_tol_scaled = this.options.dual_feasibility_tolerance / *dynamism;
                let within_dual_tol = if this.options.h16_row_presolve {
                    bound_on_col_dual.abs() <= cpp_max(f64::from(dual_tol_scaled), f64::EPSILON)
                } else {
                    bound_on_col_dual.abs() <= dual_tol_scaled
                };
                if within_dual_tol {
                    // 1. column dual's upper bound is zero (since the column's lower bound
                    // is infinite) and column dual's lower bound is zero as well
                    // (direction = 1) or
                    // 2. column dual's lower bound is zero (since the column's upper bound
                    // is infinite) and column dual's upper bound is zero as well
                    // (direction = -1).
                    // thus, the column dual is zero, and we can remove the column
                    // and all its rows
                    if logging_on {
                        this.analysis_start_presolve_rule_log(kPresolveRuleForcingCol);
                    }
                    let integral = this.model.integrality_[col as usize] == HighsVarType::kInteger;
                    let cost = this.model.col_cost_[col as usize];
                    postsolve_stack.forcing_column(
                        col,
                        this.get_column_vector(col).iter(&*this),
                        cost,
                        other_bound,
                        direction < 0,
                        integral,
                    );
                    this.mark_col_deleted(col);
                    let mut coliter = this.colhead[col as usize];
                    while coliter != -1 {
                        let row = this.a_row[coliter as usize];
                        let rhs = if (direction as f64) * this.a_value[coliter as usize] > 0.0 {
                            this.model.row_upper_[row as usize]
                        } else {
                            this.model.row_lower_[row as usize]
                        };
                        coliter = this.a_next[coliter as usize];

                        postsolve_stack.forcing_column_removed_row(
                            col,
                            row,
                            rhs,
                            this.get_row_vector(row).iter(&*this),
                        );
                        this.remove_row(row);
                    }
                    this.analysis_.logging_on_ = logging_on;
                    if logging_on {
                        this.analysis_stop_presolve_rule_log(kPresolveRuleForcingCol);
                    }
                    return this.check_limits(postsolve_stack);
                }
            }
            HPresolveResult::kOk
        };

        // check for dominated column
        let bound = self.model.col_lower_[col as usize];
        hpresolve_checked_call!(dominated_col(
            self,
            postsolve_stack,
            col,
            col_dual_lower,
            bound,
            1
        ));
        if self.col_deleted[col as usize] != 0 {
            return HPresolveResult::kOk;
        }

        let bound = self.model.col_upper_[col as usize];
        hpresolve_checked_call!(dominated_col(
            self,
            postsolve_stack,
            col,
            col_dual_upper,
            bound,
            -1
        ));
        if self.col_deleted[col as usize] != 0 {
            return HPresolveResult::kOk;
        }

        // HiGHS 1.16 #3329 (h16_presolve_rules): HPresolve.cpp lines 3759-3769 (1.16.0)
        // PORT NOTE: the 1.16 rule log of a weakly dominated column is
        // kPresolveRuleWeaklyDominatedCol; HPresolveAnalysis keeps the 1.15.1 rule arrays, so
        // the log stays on kPresolveRuleDominatedCol (MIP presolve never logs rules).
        if self.options.h16_presolve_rules
            && !self.analysis_.allow_rule_[kH16PresolveRuleWeaklyDominatedCol as usize]
        {
            return HPresolveResult::kOk;
        }

        // check for weakly dominated column
        let bound = self.model.col_lower_[col as usize];
        let other_bound = self.model.col_upper_[col as usize];
        hpresolve_checked_call!(weakly_dominated_col(
            self,
            postsolve_stack,
            col,
            col_dual_lower,
            bound,
            other_bound,
            &dynamism,
            1
        ));
        if self.col_deleted[col as usize] != 0 {
            return HPresolveResult::kOk;
        }

        let bound = self.model.col_upper_[col as usize];
        let other_bound = self.model.col_lower_[col as usize];
        hpresolve_checked_call!(weakly_dominated_col(
            self,
            postsolve_stack,
            col,
            col_dual_upper,
            bound,
            other_bound,
            &dynamism,
            -1
        ));
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::computeLocks (lines 4642-4668, 27 lines)
    ///
    /// `std::function` by value -> `&mut dyn FnMut`; callers pass `&mut |row, down, up| { ... }`.
    pub fn compute_locks(
        &self,
        col: HighsInt,
        consider_objective: bool,
        lock_callback: &mut dyn FnMut(HighsInt, bool, bool) -> bool,
    ) {
        // if the callback returns true, we stop examining the locks
        if consider_objective {
            // consider objective function
            if self.model.col_cost_[col as usize] < 0.0 {
                // downlock
                if lock_callback(-1, true, false) {
                    return;
                }
            } else if self.model.col_cost_[col as usize] > 0.0 {
                // uplock
                if lock_callback(-1, false, true) {
                    return;
                }
            }
        }

        // check coefficients
        for nz in self.get_column_vector(col).iter(self) {
            // implied lower bound -> downlock
            let has_down_lock = self.yields_implied_lower_bound(nz.index(), nz.value());
            // implied upper bound -> uplock
            let has_up_lock = self.yields_implied_upper_bound(nz.index(), nz.value());
            // callback
            if (has_down_lock || has_up_lock)
                && lock_callback(nz.index(), has_down_lock, has_up_lock)
            {
                break;
            }
        }
    }
}
