//! HiGHS source: lp_data/HighsInterface.cpp (handleInfCost, restoreInfCost)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Columns whose cost is infinite (|c| >= `options_.infinite_cost`)
//! are fixed at the bound the objective sense drives them to, with cost 0,
//! before `calledOptimizeModel` solves (`handle_inf_cost`); the costs, bounds
//! and basis statuses are put back by `returnFromOptimizeModel` when
//! `undo_mods` is true (`restore_inf_cost`). The records live in
//! `HighsLp::mods_` (`HighsLpMods::save_inf_cost_variable_*`, hstruct.rs).
//!
//! Both call sites are wired in highs_run.rs: `handle_inf_cost` only when
//! `model_lp_.has_infinite_cost_` (as in C++); `restore_inf_cost` only when
//! `undo_mods`, and its opening guard (no saved infinite-cost columns) is in
//! place so that semi-variable-only modifications never reach this code.

use crate::prelude::*;

impl Highs {
    /// HiGHS: HighsInterface.cpp::Highs::handleInfCost (lines 2960-3045).
    /// `std::ceil`/`std::floor` -> `f64::ceil`/`f64::floor`; the two passes
    /// (`k = 0` checks, `k = 1` modifies) are kept so the model is unchanged
    /// on `kError`.
    pub fn handle_inf_cost(&mut self) -> HighsStatus {
        // HiGHS: HighsInterface.cpp::Highs::handleInfCost (lines 2960-3045)
        let lp = &mut self.model_lp_;
        if !lp.has_infinite_cost_ {
            return HighsStatus::kOk;
        }
        let inf_cost = self.options_.infinite_cost;
        for k in 0..2 {
            // Pass twice: first checking that infinite costs can be handled,
            // then handling them, so that model is unmodified if infinite
            // costs cannot be handled
            for i_col in 0..lp.num_col_ {
                let cost = lp.col_cost_[i_col as usize];
                if cost > -inf_cost && cost < inf_cost {
                    continue;
                }
                let mut lower = lp.col_lower_[i_col as usize];
                let mut upper = lp.col_upper_[i_col as usize];
                if lp.is_mip() {
                    if lp.integrality_[i_col as usize] == HighsVarType::kInteger {
                        lower = lower.ceil();
                        upper = upper.floor();
                    }
                }
                if cost <= -inf_cost {
                    if lp.sense_ == ObjSense::kMinimize {
                        // Minimizing with -inf cost so try to fix at upper bound
                        if upper < kHighsInf {
                            if k != 0 {
                                lp.col_lower_[i_col as usize] = upper;
                            }
                        } else {
                            highs_log_user(
                                &self.options_.log_options,
                                HighsLogType::kError,
                                &format!(
                                    "Cannot minimize with a cost on variable {} of {} and upper bound of {}\n",
                                    i_col,
                                    crate::util::c_printf::c_fmt_g(cost, 6),
                                    crate::util::c_printf::c_fmt_g(upper, 6),
                                ),
                            );
                            return HighsStatus::kError;
                        }
                    } else {
                        // Maximizing with -inf cost so try to fix at lower bound
                        if lower > -kHighsInf {
                            if k != 0 {
                                lp.col_upper_[i_col as usize] = lower;
                            }
                        } else {
                            highs_log_user(
                                &self.options_.log_options,
                                HighsLogType::kError,
                                &format!(
                                    "Cannot maximize with a cost on variable {} of {} and lower bound of {}\n",
                                    i_col,
                                    crate::util::c_printf::c_fmt_g(cost, 6),
                                    crate::util::c_printf::c_fmt_g(lower, 6),
                                ),
                            );
                            return HighsStatus::kError;
                        }
                    }
                } else {
                    if lp.sense_ == ObjSense::kMinimize {
                        // Minimizing with inf cost so try to fix at lower bound
                        if lower > -kHighsInf {
                            if k != 0 {
                                lp.col_upper_[i_col as usize] = lower;
                            }
                        } else {
                            highs_log_user(
                                &self.options_.log_options,
                                HighsLogType::kError,
                                &format!(
                                    "Cannot minimize with a cost on variable {} of {} and lower bound of {}\n",
                                    i_col,
                                    crate::util::c_printf::c_fmt_g(cost, 6),
                                    crate::util::c_printf::c_fmt_g(lower, 6),
                                ),
                            );
                            return HighsStatus::kError;
                        }
                    } else {
                        // Maximizing with inf cost so try to fix at upper bound
                        if upper < kHighsInf {
                            if k != 0 {
                                lp.col_lower_[i_col as usize] = upper;
                            }
                        } else {
                            highs_log_user(
                                &self.options_.log_options,
                                HighsLogType::kError,
                                &format!(
                                    "Cannot maximize with a cost on variable {} of {} and upper bound of {}\n",
                                    i_col,
                                    crate::util::c_printf::c_fmt_g(cost, 6),
                                    crate::util::c_printf::c_fmt_g(upper, 6),
                                ),
                            );
                            return HighsStatus::kError;
                        }
                    }
                }
                if k != 0 {
                    lp.mods_.save_inf_cost_variable_index.push(i_col);
                    lp.mods_.save_inf_cost_variable_cost.push(cost);
                    lp.mods_.save_inf_cost_variable_lower.push(lower);
                    lp.mods_.save_inf_cost_variable_upper.push(upper);
                    lp.col_cost_[i_col as usize] = 0.0;
                }
            }
        }
        // Infinite costs have been removed, but their presence in the
        // original model is known from lp.mods_.save_inf_cost_variable_*, so
        // set lp.has_infinite_cost_ to be false to avoid assert when run()
        // is called using copy of model in MIP solver (See #1446)
        lp.has_infinite_cost_ = false;

        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::restoreInfCost (lines 3053-3091),
    /// `void restoreInfCost(HighsStatus& return_status)`.
    ///
    /// The opening guard (lines 3057-3058) is in place.
    pub fn restore_inf_cost(&mut self, return_status: &mut HighsStatus) {
        // HiGHS: HighsInterface.cpp::Highs::restoreInfCost (lines 3053-3091)
        let num_inf_cost = self.model_lp_.mods_.save_inf_cost_variable_index.len() as HighsInt;
        if num_inf_cost <= 0 {
            return;
        }
        debug_assert!(num_inf_cost != 0);
        for ix in 0..num_inf_cost {
            let i_col = self.model_lp_.mods_.save_inf_cost_variable_index[ix as usize];
            let cost = self.model_lp_.mods_.save_inf_cost_variable_cost[ix as usize];
            let lower = self.model_lp_.mods_.save_inf_cost_variable_lower[ix as usize];
            let upper = self.model_lp_.mods_.save_inf_cost_variable_upper[ix as usize];
            let value = if self.solution_.value_valid {
                self.solution_.col_value[i_col as usize]
            } else {
                0.0
            };
            if self.basis_.valid {
                debug_assert!(self.basis_.col_status[i_col as usize] != HighsBasisStatus::kBasic);
                if self.model_lp_.col_lower_[i_col as usize] == lower {
                    self.basis_.col_status[i_col as usize] = HighsBasisStatus::kLower;
                } else {
                    self.basis_.col_status[i_col as usize] = HighsBasisStatus::kUpper;
                }
            }
            debug_assert!(self.model_lp_.col_cost_[i_col as usize] == 0.0);
            if value != 0.0 {
                self.info_.objective_function_value += value * cost;
            }
            self.model_lp_.col_cost_[i_col as usize] = cost;
            self.model_lp_.col_lower_[i_col as usize] = lower;
            self.model_lp_.col_upper_[i_col as usize] = upper;
        }
        // Infinite costs have been reintroduced, so reset to true the flag
        // that was set false in Highs::handleInfCost() (See #1446)
        self.model_lp_.has_infinite_cost_ = true;

        if self.model_status_ == HighsModelStatus::kInfeasible {
            // Model is infeasible with the infinite cost variables fixed at
            // appropriate values, so model status cannot be determined
            self.model_status_ = HighsModelStatus::kUnknown;
            self.set_highs_model_status_and_clear_solution_and_basis(self.model_status_);
            *return_status = crate::lp_data::highs_solve::highs_status_from_highs_model_status(
                self.model_status_,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// min x0 + inf*x1, x0 in [0, 4], x1 in [1, 3], no rows.
    fn inf_cost_highs(sense: ObjSense) -> Highs {
        let mut highs = Highs::default();
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.col_cost_ = vec![1.0, kHighsInf];
        lp.col_lower_ = vec![0.0, 1.0];
        lp.col_upper_ = vec![4.0, 3.0];
        lp.a_matrix_.start_ = vec![0, 0, 0];
        lp.sense_ = sense;
        lp.has_infinite_cost_ = true;
        highs.model_lp_ = lp;
        highs
    }

    #[test]
    fn handle_then_restore_round_trips() {
        let mut highs = inf_cost_highs(ObjSense::kMinimize);
        let before = highs.model_lp_.clone();
        assert_eq!(highs.handle_inf_cost(), HighsStatus::kOk);
        // Minimizing with +inf cost: fixed at its lower bound, cost zeroed
        assert_eq!(highs.model_lp_.col_upper_[1], 1.0);
        assert_eq!(highs.model_lp_.col_cost_[1], 0.0);
        assert!(!highs.model_lp_.has_infinite_cost_);
        assert_eq!(highs.model_lp_.mods_.save_inf_cost_variable_index, vec![1]);
        let mut status = HighsStatus::kOk;
        highs.restore_inf_cost(&mut status);
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(highs.model_lp_.col_cost_, before.col_cost_);
        assert_eq!(highs.model_lp_.col_upper_, before.col_upper_);
        assert!(highs.model_lp_.has_infinite_cost_);
    }

    #[test]
    fn unbounded_fix_is_an_error_and_leaves_model_unchanged() {
        // Maximizing with +inf cost needs a finite upper bound
        let mut highs = inf_cost_highs(ObjSense::kMaximize);
        highs.model_lp_.col_upper_[1] = kHighsInf;
        let before = highs.model_lp_.clone();
        assert_eq!(highs.handle_inf_cost(), HighsStatus::kError);
        assert_eq!(highs.model_lp_.col_lower_, before.col_lower_);
        assert_eq!(highs.model_lp_.col_cost_, before.col_cost_);
        assert!(
            highs
                .model_lp_
                .mods_
                .save_inf_cost_variable_index
                .is_empty()
        );
    }

    #[test]
    fn minus_inf_cost_maximize_fixes_at_lower() {
        // HighsInterface.cpp:2993-2995: maximizing with -inf cost fixes at lower
        let mut highs = inf_cost_highs(ObjSense::kMaximize);
        highs.model_lp_.col_cost_[1] = -kHighsInf;
        assert_eq!(highs.handle_inf_cost(), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.col_lower_[1], 1.0);
        assert_eq!(highs.model_lp_.col_upper_[1], 1.0);
        assert_eq!(
            highs.model_lp_.mods_.save_inf_cost_variable_cost,
            vec![-kHighsInf]
        );
    }

    #[test]
    fn integer_bounds_are_rounded_and_saved_rounded() {
        // HighsInterface.cpp:2974-2979 and 3029-3035: the rounded bounds are
        // used for the fixing and are what is saved (and later restored)
        let mut highs = inf_cost_highs(ObjSense::kMinimize);
        highs.model_lp_.col_lower_[1] = 0.5;
        highs.model_lp_.col_upper_[1] = 2.5;
        highs.model_lp_.integrality_ = vec![HighsVarType::kContinuous, HighsVarType::kInteger];
        assert_eq!(highs.handle_inf_cost(), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.col_lower_[1], 0.5);
        assert_eq!(highs.model_lp_.col_upper_[1], 1.0);
        assert_eq!(
            highs.model_lp_.mods_.save_inf_cost_variable_lower,
            vec![1.0]
        );
        assert_eq!(
            highs.model_lp_.mods_.save_inf_cost_variable_upper,
            vec![2.0]
        );
        let mut status = HighsStatus::kOk;
        highs.restore_inf_cost(&mut status);
        assert_eq!(highs.model_lp_.col_lower_[1], 1.0);
        assert_eq!(highs.model_lp_.col_upper_[1], 2.0);
    }

    #[test]
    fn restore_sets_basis_objective_and_infeasible_status() {
        let mut highs = inf_cost_highs(ObjSense::kMinimize);
        // -inf cost minimizing: fixed at upper (3), so col_lower_ != lower
        highs.model_lp_.col_cost_[1] = -kHighsInf;
        assert_eq!(highs.handle_inf_cost(), HighsStatus::kOk);
        highs.solution_.value_valid = true;
        highs.solution_.col_value = vec![0.0, 3.0];
        highs.basis_.valid = true;
        highs.basis_.col_status = vec![HighsBasisStatus::kBasic, HighsBasisStatus::kLower];
        highs.info_.objective_function_value = 2.0;
        let mut status = HighsStatus::kOk;
        highs.restore_inf_cost(&mut status);
        assert_eq!(highs.basis_.col_status[1], HighsBasisStatus::kUpper);
        assert_eq!(highs.info_.objective_function_value, -kHighsInf);
        assert_eq!(status, HighsStatus::kOk);

        // HighsInterface.cpp:3084-3090: kInfeasible becomes kUnknown
        let mut highs = inf_cost_highs(ObjSense::kMinimize);
        assert_eq!(highs.handle_inf_cost(), HighsStatus::kOk);
        highs.model_status_ = HighsModelStatus::kInfeasible;
        highs.solution_.value_valid = true;
        highs.solution_.col_value = vec![0.0, 1.0];
        highs.basis_.valid = true;
        highs.basis_.col_status = vec![HighsBasisStatus::kBasic, HighsBasisStatus::kUpper];
        let mut status = HighsStatus::kOk;
        highs.restore_inf_cost(&mut status);
        assert_eq!(highs.model_status_, HighsModelStatus::kUnknown);
        assert_eq!(status, HighsStatus::kWarning);
        assert!(!highs.solution_.value_valid);
        assert!(!highs.basis_.valid);
    }

    #[test]
    fn restore_without_records_is_a_no_op() {
        let mut highs = Highs::default();
        let mut status = HighsStatus::kWarning;
        highs.restore_inf_cost(&mut status);
        assert_eq!(status, HighsStatus::kWarning);
    }
}
