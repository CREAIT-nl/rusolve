//! HiGHS source: lp_data/Highs.cpp (completeSolutionFromDiscreteAssignment,
//! assessPrimalSolution, getSavedMipSolutions), lp_data/HighsLpUtils.cpp
//! (assessLpPrimalSolution, assessColPrimalSolution)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Completing a user's (partial) assignment of a
//! MIP's discrete variables. `calledOptimizeModel` calls
//! `complete_solution_from_discrete_assignment` when `model_.isMip() &&
//! solution_.value_valid` (wired in highs_run.rs): a feasible complete
//! solution is kept as is; otherwise the integer-valued discrete variables
//! are fixed and the LP (or a sub-MIP, when some discrete values are
//! missing or fractional) is solved by a recursive `optimize_model` call,
//! after which the bounds and integrality are restored.
//!
//! The recursion relies on `called_return_from_optimize_model` being true
//! on entry (C++ comment at Highs.cpp:1227-1233); this port does not track
//! that flag, so nothing extra is needed.

use crate::lp_data::presolve_support::calculate_row_values_quad_vec;
use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::completeSolutionFromDiscreteAssignment
    /// (lines 3976-4119).
    /// HiGHS 1.16 (#3556): Highs.cpp::completeSolutionFromDiscreteAssignment
    /// (lines 4228-4374): only the #3556 bound-infeasibility test is gated.
    pub fn complete_solution_from_discrete_assignment(&mut self) -> HighsStatus {
        // Determine whether the current solution of a MIP is feasible and,
        // if not, try to assign values to continuous variables and discrete
        // variables not at integer values to achieve a feasible
        // solution. Valuable in the case where users make a heuristic
        // (partial) assignment of discrete variables
        debug_assert!(self.model_lp_.is_mip() && self.solution_.value_valid);
        // Determine whether the solution contains undefined values, in
        // order to decide whether to check its feasibility
        let contains_undefined_values = self.solution_.has_undefined();
        if !contains_undefined_values {
            let mut valid = false;
            let mut integral = false;
            let mut feasible = false;
            // Determine whether this solution is integer feasible
            let return_status = assess_lp_primal_solution(
                "",
                &self.options_,
                &self.model_lp_,
                &self.solution_,
                &mut valid,
                &mut integral,
                &mut feasible,
            );
            debug_assert!(return_status != HighsStatus::kError);
            debug_assert!(valid);
            // If the current solution is integer feasible, then it can be
            // used by MIP solver to get a primal bound
            if feasible {
                return HighsStatus::kOk;
            }
        }
        // Save the column bounds and integrality in preparation for fixing
        // the discrete variables when user-supplied values are
        // integer
        let save_col_lower = self.model_lp_.col_lower_.clone();
        let save_col_upper = self.model_lp_.col_upper_.clone();
        let save_integrality = self.model_lp_.integrality_.clone();
        let have_integrality = self.model_lp_.integrality_.len() != 0;
        debug_assert!(have_integrality);
        // Count the number of fixed and unfixed discrete variables
        let mut num_fixed_discrete_variable: HighsInt = 0;
        let mut num_unfixed_discrete_variable: HighsInt = 0;
        // HiGHS 1.16 #3556 (h16_user_solution_assess): Highs.cpp lines 4262, 4276-4282 (1.16.0)
        let h16_user_solution_assess = self.options_.h16_user_solution_assess;
        for i_col in 0..self.model_lp_.num_col_ {
            let primal = self.solution_.col_value[i_col as usize];
            // Default value is lower bound, unless primal is integer for a
            // discrete variable
            self.solution_.col_value[i_col as usize] = self.model_lp_.col_lower_[i_col as usize];
            if self.model_lp_.integrality_[i_col as usize] == HighsVarType::kContinuous {
                continue;
            }
            // Fix discrete variable if its value is defined and integer
            if primal == kHighsUndefined {
                num_unfixed_discrete_variable += 1;
            } else {
                let lower = self.model_lp_.col_lower_[i_col as usize];
                let upper = self.model_lp_.col_upper_[i_col as usize];
                let type_ = if have_integrality {
                    self.model_lp_.integrality_[i_col as usize]
                } else {
                    HighsVarType::kContinuous
                };
                let mut col_infeasibility: f64 = 0.0;
                let mut integer_infeasibility: f64 = 0.0;
                if h16_user_solution_assess {
                    // HiGHS 1.16 #3556 (h16_user_solution_assess): Highs.cpp lines 4262,
                    // 4278-4279 (1.16.0); col_infeasibility is bound_infeasibility
                    let assess_as_mip = true;
                    assess_col_primal_solution_h16(
                        assess_as_mip,
                        &self.options_,
                        primal,
                        lower,
                        upper,
                        type_,
                        &mut col_infeasibility,
                        &mut integer_infeasibility,
                    );
                } else {
                    assess_col_primal_solution(
                        &self.options_,
                        primal,
                        lower,
                        upper,
                        type_,
                        &mut col_infeasibility,
                        &mut integer_infeasibility,
                    );
                }
                let unfixed = if h16_user_solution_assess {
                    // HiGHS 1.16 #3556 (h16_user_solution_assess): Highs.cpp lines 4280-4281 (1.16.0)
                    integer_infeasibility > self.options_.mip_feasibility_tolerance
                        || col_infeasibility > self.options_.mip_feasibility_tolerance
                } else {
                    integer_infeasibility > self.options_.mip_feasibility_tolerance
                };
                if unfixed {
                    num_unfixed_discrete_variable += 1;
                } else {
                    // Variable is integer feasible, so fix it at this value and
                    // remove its integrality
                    num_fixed_discrete_variable += 1;
                    self.model_lp_.col_lower_[i_col as usize] = primal;
                    self.model_lp_.col_upper_[i_col as usize] = primal;
                    self.model_lp_.integrality_[i_col as usize] = HighsVarType::kContinuous;
                }
            }
        }
        debug_assert!(!self.solution_.has_undefined());
        let num_discrete_variable = num_unfixed_discrete_variable + num_fixed_discrete_variable;
        let num_continuous_variable = self.model_lp_.num_col_ - num_discrete_variable;
        debug_assert!(num_continuous_variable >= 0);
        let mut call_run = true;
        let few_fixed_discrete_variables = 10 * num_fixed_discrete_variable < num_discrete_variable;
        if num_unfixed_discrete_variable == 0 {
            // Solution is integer valued
            if num_continuous_variable == 0 {
                // There are no continuous variables, so no feasible solution can be
                // deduced
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "User-supplied values of discrete variables cannot yield \
                     feasible solution\n"
                );
                call_run = false;
            } else {
                // Solve an LP, so clear all integrality
                self.model_lp_.integrality_.clear();
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Attempting to find feasible solution \
                     by solving LP for user-supplied values of discrete variables\n"
                );
            }
        } else {
            // There are unfixed discrete variables
            if few_fixed_discrete_variables {
                // Too few discrete variables are fixed so warn, but still
                // attempt to complete a feasible solution
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kWarning,
                    "User-supplied values fix only {} / {} discrete variables, \
                     so attempt to complete a feasible solution may be expensive\n",
                    num_fixed_discrete_variable,
                    num_discrete_variable
                );
            } else {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Attempting to find feasible solution \
                     by solving MIP for user-supplied values of {} / {} \
                     discrete variables\n",
                    num_fixed_discrete_variable,
                    num_discrete_variable
                );
            }
        }
        let mut return_status = HighsStatus::kOk;
        // Clear the current solution since either the user solution has
        // been used to fix (a subset of) discrete variables - so a valid
        // solution will be obtained from optimizeModel() if the local model is
        // feasible - or it's not worth using the user solution
        self.solution_.clear();
        if call_run {
            // Solve the model, using mip_max_start_nodes for
            // mip_max_nodes...
            let mip_max_nodes = self.options_.mip_max_nodes;
            self.options_.mip_max_nodes = self.options_.mip_max_start_nodes;
            // Solve the model
            self.basis_.clear();
            // profiling_ branch dropped (lines 4096-4100): profiling out of scope
            return_status = self.optimize_model();
            // resetProfiling dropped (line 4103): profiling out of scope
            // ... remembering to recover the original value of mip_max_nodes
            self.options_.mip_max_nodes = mip_max_nodes;
        }
        // Recover the column bounds and integrality
        self.model_lp_.col_lower_ = save_col_lower;
        self.model_lp_.col_upper_ = save_col_upper;
        self.model_lp_.integrality_ = save_integrality;
        // Handle the error return
        if return_status == HighsStatus::kError {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Highs::optimizeModel() error trying to find feasible solution\n"
            );
            return HighsStatus::kError;
        }
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::assessPrimalSolution (lines 3681-3685)
    pub fn assess_primal_solution(
        &self,
        valid: &mut bool,
        integral: &mut bool,
        feasible: &mut bool,
    ) -> HighsStatus {
        assess_lp_primal_solution(
            "",
            &self.options_,
            &self.model_lp_,
            &self.solution_,
            valid,
            integral,
            feasible,
        )
    }

    /// HiGHS: Highs.h::Highs::getSavedMipSolutions (lines 570-572)
    pub fn get_saved_mip_solutions(&self) -> &[HighsObjectiveSolution] {
        &self.saved_objective_and_solution_
    }
}

/// HiGHS: HighsLpUtils.cpp::assessColPrimalSolution (lines 2591-2619)
pub fn assess_col_primal_solution(
    options: &HighsOptions,
    primal: f64,
    lower: f64,
    upper: f64,
    type_: HighsVarType,
    col_infeasibility: &mut f64,
    integer_infeasibility: &mut f64,
) {
    // @primal_infeasibility calculation
    *col_infeasibility = 0.0;
    if primal < lower - options.primal_feasibility_tolerance {
        *col_infeasibility = lower - primal;
    } else if primal > upper + options.primal_feasibility_tolerance {
        *col_infeasibility = primal - upper;
    }
    *integer_infeasibility = 0.0;
    if type_ == HighsVarType::kInteger || type_ == HighsVarType::kSemiInteger {
        *integer_infeasibility = fractionality(primal, None);
    }
    if *col_infeasibility > 0.0
        && (type_ == HighsVarType::kSemiContinuous || type_ == HighsVarType::kSemiInteger)
    {
        // Semi-variables at zero will have positive col
        // infeasibility, so possibly zero this
        if primal.abs() <= options.mip_feasibility_tolerance {
            *col_infeasibility = 0.0;
        }
        // If there is (still) column infeasibility due to value being
        // off zero or below lower bound, then this counts as an integer
        // infeasibility
        // C++ tests `col_infeasibility` as a bool: nonzero
        if *col_infeasibility != 0.0 && primal < upper {
            *integer_infeasibility = cpp_max(*col_infeasibility, *integer_infeasibility);
        }
    }
}

/// HiGHS 1.16 (#3525, the callee #3556 calls): HighsLpUtils.cpp::assessColPrimalSolution
/// (lines 2606-2643), with the leading `assess_as_mip` parameter. Called only
/// when `h16_user_solution_assess` is on (from
/// `complete_solution_from_discrete_assignment` and `assess_user_col_values`);
/// the 1.15.1 `assess_col_primal_solution` above is unchanged.
#[allow(clippy::too_many_arguments)]
pub fn assess_col_primal_solution_h16(
    assess_as_mip: bool,
    options: &HighsOptions,
    primal: f64,
    lower: f64,
    upper: f64,
    type_: HighsVarType,
    col_infeasibility: &mut f64,
    integer_infeasibility: &mut f64,
) {
    let k_primal_feasibility_tolerance = if assess_as_mip {
        options.mip_feasibility_tolerance
    } else {
        options.primal_feasibility_tolerance
    };
    // @primal_infeasibility calculation
    *col_infeasibility = 0.0;
    if primal < lower - k_primal_feasibility_tolerance {
        *col_infeasibility = lower - primal;
    } else if primal > upper + k_primal_feasibility_tolerance {
        *col_infeasibility = primal - upper;
    }
    *integer_infeasibility = 0.0;
    if assess_as_mip {
        if type_ == HighsVarType::kInteger || type_ == HighsVarType::kSemiInteger {
            *integer_infeasibility = fractionality(primal, None);
            debug_assert!(*integer_infeasibility >= 0.0);
            if *integer_infeasibility <= options.mip_feasibility_tolerance {
                *integer_infeasibility = 0.0;
            }
        }
        if *col_infeasibility > 0.0
            && (type_ == HighsVarType::kSemiContinuous || type_ == HighsVarType::kSemiInteger)
        {
            // Semi-variables at zero will have positive col infeasibility,
            // so possibly zero this
            if primal.abs() <= options.mip_feasibility_tolerance {
                *col_infeasibility = 0.0;
            }
            // If there is (still) column infeasibility due to value being
            // off zero or below lower bound, then this counts as an integer
            // infeasibility
            // C++ tests `col_infeasibility` as a bool: nonzero
            if *col_infeasibility != 0.0 && primal < upper {
                *integer_infeasibility = cpp_max(*col_infeasibility, *integer_infeasibility);
            }
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::assessLpPrimalSolution (lines 2623-2769).
/// The `highsLogUser` reports are kept (they are user-visible output of
/// `completeSolutionFromDiscreteAssignment`); `message` is printed first.
pub fn assess_lp_primal_solution(
    message: &str,
    options: &HighsOptions,
    lp: &HighsLp,
    solution: &HighsSolution,
    valid: &mut bool,
    integral: &mut bool,
    feasible: &mut bool,
) -> HighsStatus {
    *valid = false;
    *integral = false;
    *feasible = false;
    let mut num_col_infeasibilities: HighsInt = 0;
    let mut max_col_infeasibility: f64 = 0.0;
    let mut sum_col_infeasibilities: f64 = 0.0;
    let mut num_integer_infeasibilities: HighsInt = 0;
    let mut max_integer_infeasibility: f64 = 0.0;
    let mut sum_integer_infeasibilities: f64 = 0.0;
    let mut num_row_infeasibilities: HighsInt = 0;
    let mut max_row_infeasibility: f64 = 0.0;
    let mut sum_row_infeasibilities: f64 = 0.0;
    let mut num_row_residuals: HighsInt = 0;
    let mut max_row_residual: f64 = 0.0;
    let mut sum_row_residuals: f64 = 0.0;
    let k_row_residual_tolerance = options.primal_feasibility_tolerance; // 1e-12;
    let k_primal_feasibility_tolerance = if lp.is_mip() {
        options.mip_feasibility_tolerance
    } else {
        options.primal_feasibility_tolerance
    };
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kInfo,
        "{}Assessing feasibility of {} tolerance of {}\n",
        message,
        if lp.is_mip() {
            "MIP using primal feasibility and integrality"
        } else {
            "LP using primal feasibility"
        },
        crate::util::c_printf::c_fmt_g_width(k_primal_feasibility_tolerance, 4, 11)
    );
    let mut row_value: Vec<f64> = vec![0.0; lp.num_row_ as usize];
    let have_integrality = lp.integrality_.len() != 0;
    if !solution.value_valid {
        return HighsStatus::kError;
    }
    let have_col_names = lp.col_names_.len() == lp.num_col_ as usize;
    for i_col in 0..lp.num_col_ {
        let name_string = if have_col_names {
            format!(" ({})", lp.col_names_[i_col as usize])
        } else {
            String::new()
        };
        let primal = solution.col_value[i_col as usize];
        let lower = lp.col_lower_[i_col as usize];
        let upper = lp.col_upper_[i_col as usize];
        let type_ = if have_integrality {
            lp.integrality_[i_col as usize]
        } else {
            HighsVarType::kContinuous
        };

        let mut col_infeasibility: f64 = 0.0;
        let mut integer_infeasibility: f64 = 0.0;

        assess_col_primal_solution(
            options,
            primal,
            lower,
            upper,
            type_,
            &mut col_infeasibility,
            &mut integer_infeasibility,
        );
        if col_infeasibility > 0.0 {
            if col_infeasibility > k_primal_feasibility_tolerance {
                if col_infeasibility > 2.0 * max_col_infeasibility {
                    crate::highs_log_user!(
                        &options.log_options,
                        HighsLogType::kWarning,
                        "Col {:6}{} has         infeasibility of {} from \
                         [lower, value, upper] = [{}; {}; {}]\n",
                        i_col,
                        name_string,
                        crate::util::c_printf::c_fmt_g_width(col_infeasibility, 4, 11),
                        crate::util::c_printf::c_fmt_g_width(lower, 8, 15),
                        crate::util::c_printf::c_fmt_g_width(primal, 8, 15),
                        crate::util::c_printf::c_fmt_g_width(upper, 8, 15)
                    );
                }
                num_col_infeasibilities += 1;
            }
            max_col_infeasibility = cpp_max(col_infeasibility, max_col_infeasibility);
            sum_col_infeasibilities += col_infeasibility;
        }
        if integer_infeasibility > 0.0 {
            if integer_infeasibility > options.mip_feasibility_tolerance {
                if integer_infeasibility > 2.0 * max_integer_infeasibility {
                    crate::highs_log_user!(
                        &options.log_options,
                        HighsLogType::kWarning,
                        "Col {:6}{} has integer infeasibility of {}\n",
                        i_col,
                        name_string,
                        crate::util::c_printf::c_fmt_g_width(integer_infeasibility, 4, 11)
                    );
                }
                num_integer_infeasibilities += 1;
            }
            max_integer_infeasibility = cpp_max(integer_infeasibility, max_integer_infeasibility);
            sum_integer_infeasibilities += integer_infeasibility;
        }
    }
    let return_status = calculate_row_values_quad_vec(lp, &solution.col_value, &mut row_value, -1);
    if return_status != HighsStatus::kOk {
        return return_status;
    }
    let have_row_names = lp.row_names_.len() >= lp.num_row_ as usize;
    for i_row in 0..lp.num_row_ {
        let primal = solution.row_value[i_row as usize];
        let lower = lp.row_lower_[i_row as usize];
        let upper = lp.row_upper_[i_row as usize];
        let name_string = if have_row_names {
            format!(" ({})", lp.row_names_[i_row as usize])
        } else {
            String::new()
        };
        // @primal_infeasibility calculation
        let mut row_infeasibility: f64 = 0.0;
        if primal < lower - k_primal_feasibility_tolerance {
            row_infeasibility = lower - primal;
        } else if primal > upper + k_primal_feasibility_tolerance {
            row_infeasibility = primal - upper;
        }
        if row_infeasibility > 0.0 {
            if row_infeasibility > k_primal_feasibility_tolerance {
                if row_infeasibility > 2.0 * max_row_infeasibility {
                    crate::highs_log_user!(
                        &options.log_options,
                        HighsLogType::kWarning,
                        "Row {:6}{} has         infeasibility of {} from \
                         [lower, value, upper] = [{}; {}; {}]\n",
                        i_row,
                        name_string,
                        crate::util::c_printf::c_fmt_g_width(row_infeasibility, 4, 11),
                        crate::util::c_printf::c_fmt_g_width(lower, 8, 15),
                        crate::util::c_printf::c_fmt_g_width(primal, 8, 15),
                        crate::util::c_printf::c_fmt_g_width(upper, 8, 15)
                    );
                }
                num_row_infeasibilities += 1;
            }
            max_row_infeasibility = cpp_max(row_infeasibility, max_row_infeasibility);
            sum_row_infeasibilities += row_infeasibility;
        }
        let row_residual = (primal - row_value[i_row as usize]).abs();
        if row_residual > k_row_residual_tolerance {
            if row_residual > 2.0 * max_row_residual {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kWarning,
                    "Row {:6}{} has         residual      of {}\n",
                    i_row,
                    name_string,
                    crate::util::c_printf::c_fmt_g_width(row_residual, 4, 11)
                );
            }
            num_row_residuals += 1;
        }
        max_row_residual = cpp_max(row_residual, max_row_residual);
        sum_row_residuals += row_residual;
    }
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kInfo,
        "Solution has               num          max          sum\n"
    );
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kInfo,
        "Col     infeasibilities {:6}  {}  {}\n",
        num_col_infeasibilities,
        crate::util::c_printf::c_fmt_g_width(max_col_infeasibility, 4, 11),
        crate::util::c_printf::c_fmt_g_width(sum_col_infeasibilities, 4, 11)
    );
    if lp.is_mip() {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "Integer infeasibilities {:6}  {}  {}\n",
            num_integer_infeasibilities,
            crate::util::c_printf::c_fmt_g_width(max_integer_infeasibility, 4, 11),
            crate::util::c_printf::c_fmt_g_width(sum_integer_infeasibilities, 4, 11)
        );
    }
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kInfo,
        "Row     infeasibilities {:6}  {}  {}\n",
        num_row_infeasibilities,
        crate::util::c_printf::c_fmt_g_width(max_row_infeasibility, 4, 11),
        crate::util::c_printf::c_fmt_g_width(sum_row_infeasibilities, 4, 11)
    );
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kInfo,
        "Row     residuals       {:6}  {}  {}\n",
        num_row_residuals,
        crate::util::c_printf::c_fmt_g_width(max_row_residual, 4, 11),
        crate::util::c_printf::c_fmt_g_width(sum_row_residuals, 4, 11)
    );
    *valid = num_row_residuals == 0;
    *integral = *valid && num_integer_infeasibilities == 0;
    *feasible = *valid
        && num_col_infeasibilities == 0
        && num_integer_infeasibilities == 0
        && num_row_infeasibilities == 0;
    if !(*integral && *feasible) {
        return HighsStatus::kWarning;
    }
    HighsStatus::kOk
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// x0 integer in [0, 3], x1 continuous in [0, 5]; row x0 + x1 in [2, 4].
    fn mip_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 2.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![3.0, 5.0];
        lp.row_lower_ = vec![2.0];
        lp.row_upper_ = vec![4.0];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        lp
    }

    fn solution(col_value: Vec<f64>, row_value: Vec<f64>) -> HighsSolution {
        let mut s = HighsSolution::default();
        s.col_value = col_value;
        s.row_value = row_value;
        s.value_valid = true;
        s
    }

    #[test]
    fn col_assessment_measures_bound_and_integer_infeasibility() {
        let options = HighsOptions::default();
        let (mut col_inf, mut int_inf) = (0.0, 0.0);
        assess_col_primal_solution(
            &options,
            3.5,
            0.0,
            3.0,
            HighsVarType::kInteger,
            &mut col_inf,
            &mut int_inf,
        );
        assert_eq!(col_inf, 0.5);
        assert_eq!(int_inf, 0.5);
        assess_col_primal_solution(
            &options,
            1.0,
            0.0,
            3.0,
            HighsVarType::kContinuous,
            &mut col_inf,
            &mut int_inf,
        );
        assert_eq!((col_inf, int_inf), (0.0, 0.0));
    }

    #[test]
    fn lp_assessment_flags_fractional_and_infeasible() {
        let options = HighsOptions::default();
        let lp = mip_lp();
        let (mut valid, mut integral, mut feasible) = (false, false, false);
        let ok = solution(vec![1.0, 2.0], vec![3.0]);
        assert_eq!(
            assess_lp_primal_solution(
                "",
                &options,
                &lp,
                &ok,
                &mut valid,
                &mut integral,
                &mut feasible
            ),
            HighsStatus::kOk
        );
        assert!(valid && integral && feasible);
        let fractional = solution(vec![1.5, 1.0], vec![2.5]);
        assess_lp_primal_solution(
            "",
            &options,
            &lp,
            &fractional,
            &mut valid,
            &mut integral,
            &mut feasible,
        );
        assert!(valid && !integral && !feasible);
    }

    #[test]
    fn feasible_solution_is_kept() {
        let mut highs = Highs::default();
        highs.model_lp_ = mip_lp();
        highs.solution_ = solution(vec![1.0, 2.0], vec![3.0]);
        assert_eq!(
            highs.complete_solution_from_discrete_assignment(),
            HighsStatus::kOk
        );
        assert_eq!(highs.solution_.col_value, vec![1.0, 2.0]);
        assert!(highs.solution_.value_valid);
    }

    /// Semi-variable branch of assessColPrimalSolution (HighsLpUtils.cpp 2607-2618).
    #[test]
    fn col_assessment_semi_variables() {
        let options = HighsOptions::default();
        let (mut col_inf, mut int_inf) = (0.0, 0.0);
        // Semi-continuous at zero: column infeasibility zeroed
        assess_col_primal_solution(
            &options,
            0.0,
            2.0,
            5.0,
            HighsVarType::kSemiContinuous,
            &mut col_inf,
            &mut int_inf,
        );
        assert_eq!((col_inf, int_inf), (0.0, 0.0));
        // Semi-continuous off zero below lower: counts as integer infeasibility
        assess_col_primal_solution(
            &options,
            1.0,
            2.0,
            5.0,
            HighsVarType::kSemiContinuous,
            &mut col_inf,
            &mut int_inf,
        );
        assert_eq!((col_inf, int_inf), (1.0, 1.0));
        // Semi-integer: max of column and fractionality
        assess_col_primal_solution(
            &options,
            1.25,
            2.0,
            5.0,
            HighsVarType::kSemiInteger,
            &mut col_inf,
            &mut int_inf,
        );
        assert_eq!((col_inf, int_inf), (0.75, 0.75));
    }

    /// No continuous variables and all discrete values integer: no run
    /// (Highs.cpp 4044-4051), solution cleared, bounds and integrality restored.
    #[test]
    fn integer_assignment_without_continuous_variables_is_not_run() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![3.0];
        lp.row_lower_ = vec![2.0];
        lp.row_upper_ = vec![4.0];
        lp.integrality_ = vec![HighsVarType::kInteger];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let mut highs = Highs::default();
        highs.model_lp_ = lp;
        // x0 = 1 is integer but violates the row lower bound
        highs.solution_ = solution(vec![1.0], vec![1.0]);
        assert_eq!(
            highs.complete_solution_from_discrete_assignment(),
            HighsStatus::kOk
        );
        assert!(!highs.solution_.value_valid);
        assert!(highs.solution_.col_value.is_empty());
        assert_eq!(highs.model_lp_.col_lower_, vec![0.0]);
        assert_eq!(highs.model_lp_.col_upper_, vec![3.0]);
        assert_eq!(highs.model_lp_.integrality_, vec![HighsVarType::kInteger]);
    }

    // ---- HiGHS 1.16 (#3556, h16_user_solution_assess) ----

    /// x0 integer in [57, 75], no continuous variable; row x0 in [0, 100].
    fn bound_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![57.0];
        lp.col_upper_ = vec![75.0];
        lp.row_lower_ = vec![0.0];
        lp.row_upper_ = vec![100.0];
        lp.integrality_ = vec![HighsVarType::kInteger];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        lp
    }

    /// check-set-mip-solution test4 (3556.diff 176-205): an integer but
    /// bound-infeasible value (50 for [57, 75]) is no longer fixed. 1.16
    /// Highs.cpp 4280-4282: counted unfixed, so the MIP is solved (x0 = 57).
    #[test]
    fn h16_usersol_on_bound_infeasible_integer_is_not_fixed() {
        let mut highs = Highs::default();
        highs.options_.output_flag = false;
        highs.options_.h16_user_solution_assess = true;
        highs.model_lp_ = bound_lp();
        highs.solution_ = solution(vec![50.0], vec![50.0]);
        assert_eq!(
            highs.complete_solution_from_discrete_assignment(),
            HighsStatus::kOk
        );
        assert!(highs.solution_.value_valid);
        assert_eq!(highs.solution_.col_value, vec![57.0]);
        assert_eq!(highs.model_lp_.col_lower_, vec![57.0]);
        assert_eq!(highs.model_lp_.col_upper_, vec![75.0]);
        assert_eq!(highs.model_lp_.integrality_, vec![HighsVarType::kInteger]);
    }

    /// 1.15.1 Highs.cpp 4031-4041: the same value is fixed, every discrete
    /// variable is fixed and there is no continuous variable, so no run.
    #[test]
    fn h16_usersol_off_bound_infeasible_integer_is_fixed() {
        let mut highs = Highs::default();
        highs.options_.output_flag = false;
        highs.model_lp_ = bound_lp();
        highs.solution_ = solution(vec![50.0], vec![50.0]);
        assert_eq!(
            highs.complete_solution_from_discrete_assignment(),
            HighsStatus::kOk
        );
        assert!(!highs.solution_.value_valid);
        assert!(highs.solution_.col_value.is_empty());
    }

    /// 1.16 assessColPrimalSolution (HighsLpUtils.cpp 2606-2643).
    #[test]
    fn h16_usersol_col_assessment_as_mip_or_lp() {
        let mut options = HighsOptions::default();
        let (mut col_inf, mut int_inf) = (0.0, 0.0);
        // Not as MIP: no integrality test
        assess_col_primal_solution_h16(
            false,
            &options,
            3.5,
            0.0,
            3.0,
            HighsVarType::kInteger,
            &mut col_inf,
            &mut int_inf,
        );
        assert_eq!((col_inf, int_inf), (0.5, 0.0));
        // As MIP: the bound test uses mip_feasibility_tolerance
        options.mip_feasibility_tolerance = 1e-3;
        assess_col_primal_solution_h16(
            true,
            &options,
            3.0005,
            0.0,
            3.0,
            HighsVarType::kInteger,
            &mut col_inf,
            &mut int_inf,
        );
        assert_eq!((col_inf, int_inf), (0.0, 0.0));
        assess_col_primal_solution_h16(
            true,
            &options,
            1.25,
            2.0,
            5.0,
            HighsVarType::kSemiInteger,
            &mut col_inf,
            &mut int_inf,
        );
        assert_eq!((col_inf, int_inf), (0.75, 0.75));
    }
}
