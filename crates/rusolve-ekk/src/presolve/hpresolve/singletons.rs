//! HiGHS source: presolve/HPresolve.cpp (lines 3036-3438)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Doubleton equations, singleton rows/columns, free-column substitution. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::doubletonEq (lines 3036-3203, 168 lines)
    pub fn doubleton_eq(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row: HighsInt,
        row_type: RowType,
    ) -> HPresolveResult {
        debug_assert!(self.analysis_.allow_rule_[kPresolveRuleDoubletonEquation as usize]);
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleDoubletonEquation);
        }
        debug_assert!(self.row_deleted[row as usize] == 0);
        debug_assert!(self.rowsize[row as usize] == 2);
        debug_assert!(self.is_equation(row));

        // printf("doubleton equation: ");
        // debugPrintRow(row);
        let nz_pos1: HighsInt = self.rowroot[row as usize];
        let nz_pos2: HighsInt = if self.a_r_right[nz_pos1 as usize] != -1 {
            self.a_r_right[nz_pos1 as usize]
        } else {
            self.a_r_left[nz_pos1 as usize]
        };

        let col_at_pos1_better = || -> bool {
            if self.model.integrality_[self.a_col[nz_pos1 as usize] as usize]
                == HighsVarType::kInteger
            {
                if self.model.integrality_[self.a_col[nz_pos2 as usize] as usize]
                    == HighsVarType::kInteger
                {
                    // both columns integer. For substitution choose smaller absolute
                    // coefficient value, or sparser column if values are equal
                    if self.a_value[nz_pos1 as usize].abs()
                        < self.a_value[nz_pos2 as usize].abs() - self.options.small_matrix_value
                    {
                        true
                    } else if self.a_value[nz_pos2 as usize].abs()
                        < self.a_value[nz_pos1 as usize].abs() - self.options.small_matrix_value
                    {
                        false
                    } else {
                        self.colsize[self.a_col[nz_pos1 as usize] as usize]
                            < self.colsize[self.a_col[nz_pos2 as usize] as usize]
                    }
                } else {
                    // one col is integral, substitute the continuous one
                    false
                }
            } else if self.model.integrality_[self.a_col[nz_pos2 as usize] as usize]
                == HighsVarType::kInteger
            {
                // one col is integral, substitute the continuous one
                true
            } else {
                // both columns continuous the one with a larger absolute coefficient
                // value if the difference is more than factor 2, and otherwise the one
                // with fewer nonzeros if those are equal
                let col1_size = self.colsize[self.a_col[nz_pos1 as usize] as usize];
                if col1_size == 1 {
                    true
                } else {
                    let col2_size = self.colsize[self.a_col[nz_pos2 as usize] as usize];
                    if col2_size == 1 {
                        false
                    } else {
                        let abs1_val = self.a_value[nz_pos1 as usize].abs();
                        let abs2_val = self.a_value[nz_pos2 as usize].abs();
                        if col1_size != col2_size
                            && cpp_max(abs1_val, abs2_val) <= 2.0 * cpp_min(abs1_val, abs2_val)
                        {
                            col1_size < col2_size
                        } else {
                            abs1_val > abs2_val
                        }
                    }
                }
            }
        };

        let (mut substcol, mut staycol, mut substcoef, mut staycoef): (
            HighsInt,
            HighsInt,
            f64,
            f64,
        ) = if col_at_pos1_better() {
            (
                self.a_col[nz_pos1 as usize],
                self.a_col[nz_pos2 as usize],
                self.a_value[nz_pos1 as usize],
                self.a_value[nz_pos2 as usize],
            )
        } else {
            (
                self.a_col[nz_pos2 as usize],
                self.a_col[nz_pos1 as usize],
                self.a_value[nz_pos2 as usize],
                self.a_value[nz_pos1 as usize],
            )
        };

        let mut rhs: f64 = self.model.row_upper_[row as usize];
        if self.model.integrality_[substcol as usize] == HighsVarType::kInteger
            && self.model.integrality_[staycol as usize] == HighsVarType::kInteger
        {
            // check integrality conditions
            let round_coef: f64 = (staycoef / substcoef).round() * substcoef;
            if (round_coef - staycoef).abs() > self.options.small_matrix_value {
                return HPresolveResult::kOk;
            }
            staycoef = round_coef;
            let round_rhs: f64 = (rhs / substcoef).round() * substcoef;
            if (rhs - round_rhs).abs() > self.primal_feastol {
                return HPresolveResult::kPrimalInfeasible;
            }
            rhs = round_rhs;
        }

        let old_stay_lower: f64 = self.model.col_lower_[staycol as usize];
        let old_stay_upper: f64 = self.model.col_upper_[staycol as usize];
        let subst_lower: f64 = self.model.col_lower_[substcol as usize];
        let subst_upper: f64 = self.model.col_upper_[substcol as usize];

        let (stay_impl_lower, stay_impl_upper): (f64, f64) =
            if substcoef.is_sign_negative() != staycoef.is_sign_negative() {
                // coefficients have the opposite sign, therefore the implied lower bound of
                // the stay column is computed from the lower bound of the substituted
                // column:
                // staycol * staycoef + substcol * substcoef = rhs
                // staycol = (rhs - substcol * substcoef) / staycoef
                // staycol >= rhs / staycoef + lower(-substcoef/staycoef * substcol)
                // lower(-substcoef/staycoef * substcol) is (-substcoef/staycoef) *
                // substLower if (-substcoef/staycoef) is positive, i.e. if the coefficients
                // have opposite sign
                (
                    if subst_lower == -kHighsInf {
                        -kHighsInf
                    } else {
                        f64::from((HighsCDouble::from(rhs) - substcoef * subst_lower) / staycoef)
                    },
                    if subst_upper == kHighsInf {
                        kHighsInf
                    } else {
                        f64::from((HighsCDouble::from(rhs) - substcoef * subst_upper) / staycoef)
                    },
                )
            } else {
                (
                    if subst_upper == kHighsInf {
                        -kHighsInf
                    } else {
                        f64::from((HighsCDouble::from(rhs) - substcoef * subst_upper) / staycoef)
                    },
                    if subst_lower == -kHighsInf {
                        kHighsInf
                    } else {
                        f64::from((HighsCDouble::from(rhs) - substcoef * subst_lower) / staycoef)
                    },
                )
            };

        // possibly tighten bounds of the column that stays
        let lower_tightened: bool = stay_impl_lower > old_stay_lower + self.primal_feastol;
        if lower_tightened {
            hpresolve_checked_call!(self.change_col_lower(staycol, stay_impl_lower));
        }

        let upper_tightened: bool = stay_impl_upper < old_stay_upper - self.primal_feastol;
        if upper_tightened {
            hpresolve_checked_call!(self.change_col_upper(staycol, stay_impl_upper));
        }

        postsolve_stack.doubleton_equation(
            row,
            substcol,
            staycol,
            substcoef,
            staycoef,
            rhs,
            subst_lower,
            subst_upper,
            self.model.col_cost_[substcol as usize],
            lower_tightened,
            upper_tightened,
            row_type,
            self.get_column_vector(substcol).iter(self),
        );

        // finally modify matrix
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 3736-3737 at 7b890c255
        if self.options.mip_presolve_clique_table {
            hpresolve_checked_call!(self.substitute_x1b(
                substcol,
                staycol,
                rhs / substcoef,
                -staycoef / substcoef,
                row
            ));
        } else {
            self.mark_col_deleted(substcol);
            self.remove_row(row);
            self.substitute_col(substcol, staycol, rhs / substcoef, -staycoef / substcoef);
        }

        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleDoubletonEquation);
        }

        // since a column was deleted we might have new row singletons which we
        // immediately remove
        hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));

        self.check_limits(postsolve_stack)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::singletonRow (lines 3205-3346, 142 lines)
    pub fn singleton_row(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row: HighsInt,
    ) -> HPresolveResult {
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleSingletonRow);
        }
        debug_assert!(self.row_deleted[row as usize] == 0);
        debug_assert!(self.rowsize[row as usize] == 1);

        // the tree of nonzeros of this row should just contain the single nonzero
        let nz_pos: HighsInt = self.rowroot[row as usize];
        debug_assert!(nz_pos != -1);
        // nonzero should have the row in the row array
        debug_assert!(self.a_row[nz_pos as usize] == row);
        // tree with one element should not have children
        debug_assert!(self.a_r_left[nz_pos as usize] == -1);
        debug_assert!(self.a_r_right[nz_pos as usize] == -1);

        let col: HighsInt = self.a_col[nz_pos as usize];
        let val: f64 = self.a_value[nz_pos as usize];

        // printf("singleton row\n");
        // debugPrintRow(row);
        // delete row singleton nonzero directly, we have all information that we need
        // in local variables
        self.mark_row_deleted(row);
        self.unlink(nz_pos);

        // check for simple
        if val > 0.0 {
            if self.model.col_upper_[col as usize] * val
                <= self.model.row_upper_[row as usize] + self.primal_feastol
                && self.model.col_lower_[col as usize] * val
                    >= self.model.row_lower_[row as usize] - self.primal_feastol
            {
                postsolve_stack.redundant_row(row);
                self.analysis_.logging_on_ = logging_on;
                if logging_on {
                    self.analysis_stop_presolve_rule_log(kPresolveRuleSingletonRow);
                }
                return self.check_limits(postsolve_stack);
            }
        } else if self.model.col_lower_[col as usize] * val
            <= self.model.row_upper_[row as usize] + self.primal_feastol
            && self.model.col_upper_[col as usize] * val
                >= self.model.row_lower_[row as usize] - self.primal_feastol
        {
            postsolve_stack.redundant_row(row);
            self.analysis_.logging_on_ = logging_on;
            if logging_on {
                self.analysis_stop_presolve_rule_log(kPresolveRuleSingletonRow);
            }
            return self.check_limits(postsolve_stack);
        }

        // zeros should not be linked in the matrix
        debug_assert!(val.abs() > self.options.small_matrix_value);

        let mut new_col_upper: f64 = kHighsInf;
        let mut new_col_lower: f64 = -kHighsInf;
        if val > 0.0 {
            if self.model.row_upper_[row as usize] != kHighsInf {
                new_col_upper = self.model.row_upper_[row as usize] / val;
            }
            if self.model.row_lower_[row as usize] != -kHighsInf {
                new_col_lower = self.model.row_lower_[row as usize] / val;
            }
        } else {
            if self.model.row_upper_[row as usize] != kHighsInf {
                new_col_lower = self.model.row_upper_[row as usize] / val;
            }
            if self.model.row_lower_[row as usize] != -kHighsInf {
                new_col_upper = self.model.row_lower_[row as usize] / val;
            }
        }

        // use either the primal feasibility tolerance for the bound constraint or
        // for the singleton row including scaling, whichever is tighter.
        let bound_tol: f64 = if self.options.h16_row_presolve {
            // HiGHS 1.16 #3224 (h16_row_presolve): HPresolveUtils.cpp lines 49-50 (1.16.0,
            // computeSingletonRowBounds; this port keeps the bounds inline as in 1.15.1)
            cpp_max(self.primal_feastol / cpp_max(1.0, val.abs()), f64::EPSILON)
        } else {
            self.primal_feastol / cpp_max(1.0, val.abs())
        };
        let is_integral: bool = self.model.integrality_[col as usize] != HighsVarType::kContinuous;

        let mut lower_tightened: bool =
            new_col_lower > self.model.col_lower_[col as usize] + bound_tol;
        let mut upper_tightened: bool =
            new_col_upper < self.model.col_upper_[col as usize] - bound_tol;

        let mut lb: f64 = if lower_tightened {
            if is_integral {
                new_col_lower = (new_col_lower - bound_tol).ceil();
            }
            new_col_lower
        } else {
            self.model.col_lower_[col as usize]
        };

        let mut ub: f64 = if upper_tightened {
            if is_integral {
                new_col_upper = (new_col_upper + bound_tol).floor();
            }
            new_col_upper
        } else {
            self.model.col_upper_[col as usize]
        };

        // printf("old bounds [%.15g,%.15g], new bounds [%.15g,%.15g] ... ",
        //        model->col_lower_[col], model->col_upper_[col], lb, ub);
        // check whether the bounds are equal in tolerances
        if ub <= lb + self.primal_feastol {
            // bounds could be infeasible or equal in tolerances, first check infeasible
            if ub < lb - self.primal_feastol {
                return HPresolveResult::kPrimalInfeasible;
            }

            // bounds are equal in tolerances, if they have a slight infeasibility below
            // those tolerances or they have a slight numerical distance which changes
            // the largest contribution below feasibility tolerance then we can safely
            // set the bound to one of the values. To heuristically get rid of numerical
            // errors we choose the bound that was not tightened, or the midpoint if
            // both where tightened.
            if ub < lb
                || (ub > lb
                    && (ub - lb) * cpp_max(val.abs(), self.get_max_abs_col_val(col))
                        <= self.primal_feastol)
            {
                if lower_tightened && upper_tightened {
                    ub = 0.5 * (ub + lb);
                    lb = ub;
                    lower_tightened = lb > self.model.col_lower_[col as usize];
                    upper_tightened = ub < self.model.col_upper_[col as usize];
                } else if lower_tightened {
                    lb = ub;
                    lower_tightened = lb > self.model.col_lower_[col as usize];
                } else {
                    ub = lb;
                    upper_tightened = ub < self.model.col_upper_[col as usize];
                }
            }
        }

        // printf("final bounds: [%.15g,%.15g]\n", lb, ub);

        postsolve_stack.singleton_row(row, col, val, lower_tightened, upper_tightened);

        // just update bounds (and row activities)
        if lower_tightened {
            hpresolve_checked_call!(self.change_col_lower(col, lb));
        }
        // update bounds, or remove as fixed column directly
        if ub == lb {
            postsolve_stack.removed_fixed_col(
                col,
                lb,
                self.model.col_cost_[col as usize],
                self.get_column_vector(col).iter(self),
            );
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 3806 at 7b890c255
            if self.options.mip_presolve_clique_table {
                hpresolve_checked_call!(self.remove_fixed_col_x1b(col));
            } else {
                self.remove_fixed_col(col);
            }
        } else if upper_tightened {
            hpresolve_checked_call!(self.change_col_upper(col, ub));
        }

        if self.col_deleted[col as usize] == 0 && self.colsize[col as usize] == 0 {
            let result: HPresolveResult = self.empty_col(postsolve_stack, col);
            self.analysis_.logging_on_ = logging_on;
            if logging_on {
                self.analysis_stop_presolve_rule_log(kPresolveRuleSingletonRow);
            }
            return result;
        }
        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleSingletonRow);
        }
        self.check_limits(postsolve_stack)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::singletonCol (lines 3348-3418, 71 lines)
    pub fn singleton_col(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        debug_assert!(self.colsize[col as usize] == 1);
        debug_assert!(self.col_deleted[col as usize] == 0);
        let nz_pos: HighsInt = self.colhead[col as usize];
        let row: HighsInt = self.a_row[nz_pos as usize];
        let col_coef: f64 = self.a_value[nz_pos as usize];

        if self.rowsize[row as usize] == 1 {
            hpresolve_checked_call!(self.singleton_row(postsolve_stack, row));

            if self.col_deleted[col as usize] == 0 {
                debug_assert!(self.colsize[col as usize] == 0);
                return self.empty_col(postsolve_stack, col);
            }
            return HPresolveResult::kOk;
        }

        // detect strong / weak domination
        // HiGHS 1.16 #3342/#3348 (h16_dual_postsolve_rules): HPresolve.cpp::singletonCol (1.16.0
        // lines 2389-2494, the dominated column call-site guard)
        if !self.options.h16_dual_postsolve_rules
            || self.analysis_.allow_rule_[kPresolveRuleDominatedCol as usize]
        {
            hpresolve_checked_call!(self.detect_dominated_col(postsolve_stack, col, false));
        }
        if self.col_deleted[col as usize] != 0 {
            return HPresolveResult::kOk;
        }

        // check if variable is implied integer
        if self.mipsolver.is_some() {
            hpresolve_checked_call!(self.convert_implied_integer(col, row, false).result());
        }

        // dual fixing
        // HiGHS 1.16 #3342/#3348 (h16_dual_postsolve_rules): HPresolve.cpp::singletonCol (1.16.0
        // lines 2389-2494, the dual fixing call-site guard)
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

        // update column implied bounds
        hpresolve_checked_call!(self.update_col_implied_bounds(
            postsolve_stack,
            row,
            col,
            col_coef
        ));

        // update row dual implied bounds
        if self.model.integrality_[col as usize] != HighsVarType::kInteger {
            self.update_row_dual_implied_bounds(row, col, col_coef);
        }

        // now check if column is implied free within an equation and substitute the
        // column if that is the case
        if self.is_dual_implied_free(row)
            && self.is_implied_free(col)
            && self.analysis_.allow_rule_[kPresolveRuleFreeColSubstitution as usize]
        {
            if self.model.integrality_[col as usize] == HighsVarType::kInteger {
                let implied_integral: StatusResult = self.is_implied_integral(col);
                hpresolve_checked_call!(implied_integral.result());
                if !implied_integral.flag() {
                    return HPresolveResult::kOk;
                }
            }
            let logging_on = self.analysis_.logging_on_;

            if logging_on {
                self.analysis_start_presolve_rule_log(kPresolveRuleFreeColSubstitution);
            }

            // todo, store which side of an implied free dual variable needs to be used
            // for substitution
            self.store_row(row);

            self.substitute_free_col(postsolve_stack, row, col, false);

            self.analysis_.logging_on_ = logging_on;
            if logging_on {
                self.analysis_stop_presolve_rule_log(kPresolveRuleFreeColSubstitution);
            }
            return self.check_limits(postsolve_stack);
        }

        // HiGHS 1.16 #2962/#3306/#3455 (h16_zero_cost_singleton): HPresolve.cpp::singletonCol
        // (1.16.0 lines 2389-2494, the zeroCostSingleton call). The C++ gate is
        // allow_rule_[ZeroCostSingleton], which 1.16 switches off under mayRequireDualPostsolve
        // (HPresolve.cpp::chooseRules, mipsolver == nullptr); there is no rule enum here, so the
        // same condition is tested directly.
        // HiGHS 1.16 #3167 (h16_presolve_rules): with the 1.16 rule slots the gate is
        // allow_rule_[ZeroCostSingleton] itself (chooseRules switches it off).
        let h16_zcs_allowed = if self.options.h16_presolve_rules {
            self.analysis_.allow_rule_[kH16PresolveRuleZeroCostSingleton as usize]
        } else {
            !(self.mipsolver.is_none()
                && super::presolve_loop::h16_may_require_dual_postsolve(&self.options))
        };
        if self.options.h16_zero_cost_singleton && h16_zcs_allowed {
            hpresolve_checked_call!(self.zero_cost_singleton(postsolve_stack, col));
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::substituteFreeCol (lines 3420-3438, 19 lines)
    ///
    /// C++ default arguments: `relaxRowDualBounds = false`.
    pub fn substitute_free_col(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row: HighsInt,
        col: HighsInt,
        relax_row_dual_bounds: bool,
    ) {
        debug_assert!(self.row_deleted[row as usize] == 0);
        debug_assert!(self.col_deleted[col as usize] == 0);
        debug_assert!(self.is_dual_implied_free(row));
        debug_assert!(self.is_implied_free(col));

        let mut rhs: f64 = 0.0;
        let mut row_type: RowType = RowType::default();
        self.dual_implied_free_get_rhs_and_row_type(
            row,
            &mut rhs,
            &mut row_type,
            relax_row_dual_bounds,
        );

        postsolve_stack.free_col_substitution(
            row,
            col,
            rhs,
            self.model.col_cost_[col as usize],
            row_type,
            self.get_stored_row().iter(self),
            self.get_column_vector(col).iter(self),
        );

        // todo, check integrality of coefficients and allow this
        self.substitute(row, col, rhs);
    }
}
