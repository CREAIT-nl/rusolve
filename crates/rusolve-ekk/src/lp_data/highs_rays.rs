//! HiGHS source: lp_data/Highs.cpp (getFixedLp, getDualRay,
//! getDualUnboundednessDirection, getPrimalRay), lp_data/HighsInterface.cpp
//! (getDualRayInterface, getPrimalRayInterface)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Certificates of infeasibility / unboundedness and the fixed
//! LP of a MIP. The interfaces use the simplex dual / primal ray records
//! (`ekk_instance_.dual_ray_record_`, `primal_ray_record_`) when there is an
//! INVERT; otherwise they solve an LP (the original with zero costs for a
//! dual ray, or a feasibility problem for a primal ray) with presolve off,
//! then restore costs, Hessian, `solve_relaxation` and `presolve`.
//! `getDualRaySparse` is already ported (highs_api.rs).

use crate::highs_log_user;
use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::getFixedLp (lines 2096-2139).
    /// HiGHS quirk: a semi-continuous column is fixed (at its value) when
    /// `value < lp.col_lower_[iCol] - value`, as written.
    pub fn get_fixed_lp(&self, lp: &mut HighsLp) -> HighsStatus {
        if !self.model_lp_.is_mip() {
            highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Incumbent model is not a MIP, so cannot form fixed LP\n"
            );
            return HighsStatus::kError;
        }
        if !self.solution_.value_valid {
            highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Incumbent model does not have a valid solution, so cannot form fixed LP\n"
            );
            return HighsStatus::kError;
        }
        *lp = self.model_lp_.clone();
        let integrality = self.model_lp_.integrality_.clone();
        lp.integrality_.clear();
        let mut num_non_conts_fractional: HighsInt = 0;
        let mut max_fractional: f64 = 0.0;
        for i_col in 0..self.model_lp_.num_col_ {
            let value = self.solution_.col_value[i_col as usize];
            // Fix integer and semi-integer variables at their
            // value. Semi-continuous variables are fixed at zero if they are
            // closer to zero than their lower bound
            if integrality[i_col as usize] == HighsVarType::kInteger
                || integrality[i_col as usize] == HighsVarType::kSemiInteger
                || (integrality[i_col as usize] == HighsVarType::kSemiContinuous
                    && value < lp.col_lower_[i_col as usize] - value)
            {
                let fractional = fractionality(value, None);
                if fractional > self.options_.mip_feasibility_tolerance {
                    num_non_conts_fractional += 1;
                    max_fractional = cpp_max(fractional, max_fractional);
                }
                lp.col_lower_[i_col as usize] = value;
                lp.col_upper_[i_col as usize] = value;
            }
        }
        if num_non_conts_fractional != 0 {
            highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "Fixed LP has {} variables fixed at max fractional value of {}\n",
                num_non_conts_fractional as i32,
                crate::util::c_printf::c_fmt_g(max_fractional, 6)
            );
            return HighsStatus::kWarning;
        }
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::getDualRay (lines 2141-2144).
    /// `double* dual_ray_value = nullptr` -> `Option<&mut [f64]>`.
    pub fn get_dual_ray(
        &mut self,
        has_dual_ray: &mut bool,
        dual_ray_value: Option<&mut [f64]>,
    ) -> HighsStatus {
        *has_dual_ray = false;
        self.get_dual_ray_interface(has_dual_ray, dual_ray_value)
    }

    /// HiGHS: Highs.cpp::Highs::getDualUnboundednessDirection (lines 2164-2183).
    /// `productTransposeQuad` is util::sparse_matrix.
    pub fn get_dual_unboundedness_direction(
        &mut self,
        has_dual_unboundedness_direction: &mut bool,
        dual_unboundedness_direction_value: Option<&mut [f64]>,
    ) -> HighsStatus {
        if let Some(dual_unboundedness_direction_value) = dual_unboundedness_direction_value {
            let mut dual_ray_value = vec![0.0f64; self.model_lp_.num_row_ as usize];
            let status =
                self.get_dual_ray(has_dual_unboundedness_direction, Some(&mut dual_ray_value));
            if status != HighsStatus::kOk || !*has_dual_unboundedness_direction {
                return HighsStatus::kError;
            }
            let mut dual_unboundedness_direction: Vec<f64> = Vec::new();
            self.model_lp_.a_matrix_.product_transpose_quad(
                &mut dual_unboundedness_direction,
                &dual_ray_value,
                kDebugReportOff,
            );
            for i_col in 0..self.model_lp_.num_col_ {
                dual_unboundedness_direction_value[i_col as usize] =
                    dual_unboundedness_direction[i_col as usize];
            }
        } else {
            return self.get_dual_ray(has_dual_unboundedness_direction, None);
        }
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::getPrimalRay (lines 2185-2189)
    pub fn get_primal_ray(
        &mut self,
        has_primal_ray: &mut bool,
        primal_ray_value: Option<&mut [f64]>,
    ) -> HighsStatus {
        *has_primal_ray = false;
        self.get_primal_ray_interface(has_primal_ray, primal_ray_value)
    }

    /// HiGHS: HighsInterface.cpp::Highs::getDualRayInterface (lines 1571-1700).
    /// HiGHS 1.16 (#3435): HighsInterface.cpp::Highs::getDualRayInterface (lines 1580-1726)
    /// `basisSolveInterface` is highs_api_misc.rs; the nested `run()` keeps
    /// its C++ position.
    pub(crate) fn get_dual_ray_interface(
        &mut self,
        has_dual_ray: &mut bool,
        dual_ray_value: Option<&mut [f64]>,
    ) -> HighsStatus {
        // HiGHS 1.16 #3435 (h16_ray_simplex): read once at entry
        let h16_ray_simplex = self.options_.h16_ray_simplex;
        let mut return_status = HighsStatus::kOk;
        let num_row = self.model_lp_.num_row_;
        // For an LP with no rows the dual ray is vacuous
        if num_row == 0 {
            return return_status;
        }
        let mut has_invert = self.ekk_instance_.status_.has_invert;
        debug_assert!(!self.model_lp_.is_moved_);
        *has_dual_ray = self.ekk_instance_.dual_ray_record_.index != kNoRayIndex;

        // Declare identifiers to save column costs, integrality, any
        // Hessian and the presolve setting, and a flag to know when they
        // should be recovered
        let mut col_cost: Vec<f64> = Vec::new();
        let mut hessian = HighsHessian::default();
        let mut solve_relaxation = false;
        let mut presolve = String::new();
        // HiGHS 1.16 #3435 (h16_ray_simplex): HighsInterface.cpp lines 1593-1599 (1.16.0)
        let mut solver = String::new();
        let mut allow_unbounded_or_infeasible = false;
        let mut solve_feasibility_problem = false;
        let is_qp = self.model_hessian_.dim_ != 0;

        if let Some(dual_ray_value) = dual_ray_value {
            // User wants a dual ray whatever
            if !*has_dual_ray || !has_invert {
                // No dual ray is known, or no INVERT to compute it
                //
                // No point in trying to get a dual ray if the model status is
                // optimal
                if self.model_status_ == HighsModelStatus::kOptimal {
                    highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "Model status is optimal, so no dual ray is available\n"
                    );
                    return return_status;
                }
                highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Solving LP to try to compute dual ray\n"
                );
                // Save the column costs, integrality, any Hessian and the
                // presolve setting
                col_cost = self.model_lp_.col_cost_.clone();
                if is_qp {
                    hessian = self.model_hessian_.clone();
                }
                if h16_ray_simplex {
                    // HiGHS 1.16 #3435 (h16_ray_simplex): HighsInterface.cpp lines 1617-1626 (1.16.0)
                    solve_relaxation = self.options_.solve_relaxation;
                    presolve = self.options_.presolve.clone();
                    solver = self.options_.solver.clone();
                    allow_unbounded_or_infeasible = self.options_.allow_unbounded_or_infeasible;
                } else {
                    presolve = self.get_option_value_str("presolve").unwrap_or_default();
                    solve_relaxation = self
                        .get_option_value_bool("solve_relaxation")
                        .unwrap_or(false);
                }
                solve_feasibility_problem = true;
                // Zero the costs, integrality and Hessian
                let zero_costs = vec![0.0f64; self.model_lp_.num_col_ as usize];
                // Take a copy of the primal ray record, since this will be
                // cleared by calling changeColsCost
                let primal_ray_record = self.ekk_instance_.primal_ray_record_.get_ray_record();
                let status =
                    self.change_cols_cost_interval(0, self.model_lp_.num_col_ - 1, &zero_costs);
                debug_assert!(status == HighsStatus::kOk);
                // Reinstate the primal ray record
                self.ekk_instance_
                    .primal_ray_record_
                    .set_ray_record(&primal_ray_record);
                if is_qp {
                    let zero_hessian = HighsHessian::default();
                    self.pass_hessian(zero_hessian);
                }
                let call_status = if h16_ray_simplex {
                    // HiGHS 1.16 #3435 (h16_ray_simplex): HighsInterface.cpp lines 1644-1659 (1.16.0)
                    self.options_.solve_relaxation = true;
                    self.options_.presolve = kHighsOffString.to_string();
                    self.options_.solver = kSimplexString.to_string();
                    self.options_.allow_unbounded_or_infeasible = false;
                    let call_status = self.run();
                    // Restore the column costs, any Hessian, the relaxation,
                    // presolve and solver settings (1.16.0 restores them here)
                    self.options_.solve_relaxation = solve_relaxation;
                    self.options_.presolve = presolve.clone();
                    self.options_.solver = solver.clone();
                    self.options_.allow_unbounded_or_infeasible = allow_unbounded_or_infeasible;
                    self.model_lp_.col_cost_ = col_cost.clone();
                    if is_qp {
                        self.model_hessian_ = hessian.clone();
                    }
                    call_status
                } else {
                    self.set_option_value_str("presolve", kHighsOffString);
                    self.set_option_value_bool("solve_relaxation", true);
                    self.run()
                };
                if call_status != HighsStatus::kOk {
                    return_status = call_status;
                }
                *has_dual_ray = self.ekk_instance_.dual_ray_record_.index != kNoRayIndex;
                has_invert = self.ekk_instance_.status_.has_invert;
                debug_assert!(has_invert);
            }
            if *has_dual_ray {
                if !self.ekk_instance_.dual_ray_record_.value.is_empty() {
                    // Dual ray is already computed
                    highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "Copying known dual ray\n"
                    );
                    for i_row in 0..num_row {
                        dual_ray_value[i_row as usize] =
                            self.ekk_instance_.dual_ray_record_.value[i_row as usize];
                    }
                } else if has_invert {
                    // Dual ray is known and can be calculated
                    highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "Solving linear system to compute dual ray\n"
                    );
                    let mut rhs = vec![0.0f64; num_row as usize];
                    let i_row = self.ekk_instance_.dual_ray_record_.index;
                    rhs[i_row as usize] = self.ekk_instance_.dual_ray_record_.sign as f64;
                    // C++: HighsInt* dual_ray_num_nz = 0 (a null pointer)
                    self.basis_solve_interface(&rhs, dual_ray_value, None, None, true);
                    // Now save the dual ray itself
                    self.ekk_instance_
                        .dual_ray_record_
                        .value
                        .resize(num_row as usize, 0.0);
                    for i_row in 0..num_row {
                        self.ekk_instance_.dual_ray_record_.value[i_row as usize] =
                            dual_ray_value[i_row as usize];
                    }
                } else {
                    debug_assert!(!has_invert);
                    // Dual ray is known but cannot be calculated
                    highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "No LP invertible representation to compute dual ray\n"
                    );
                    return_status = HighsStatus::kError;
                }
            } else {
                highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "No dual ray found\n"
                );
                return_status = HighsStatus::kOk;
            }
        }
        if solve_feasibility_problem {
            // Feasibility problem has been solved, so any objective-related
            // information has been lost. Reverting the objective function via
            // Highs calls clears info_, so better to just copy the data
            // directly and set the info_ entries that are no longer valid
            // HiGHS 1.16 #3435 (h16_ray_simplex): with the flag on, 1.16.0 has already
            // restored the costs, Hessian and options after run() (lines 1658-1659)
            if !h16_ray_simplex {
                self.model_lp_.col_cost_ = col_cost;
                if is_qp {
                    self.model_hessian_ = hessian;
                }
                self.set_option_value_str("presolve", &presolve);
                self.set_option_value_bool("solve_relaxation", solve_relaxation);
            }
            // The relaxation for an infeasible MIP may be feasible - so no
            // ray is generated - so make sure (#2415) that the primal
            // solution status is reset
            self.info_.primal_solution_status = kSolutionStatusNone;
            // Modify the objective-related information
            self.info_.dual_solution_status = kSolutionStatusNone;
            self.info_.objective_function_value = 0.0;
            self.info_.invalidate_dual_kkt();
            if *has_dual_ray {
                // HiGHS 1.16 #3435 (h16_ray_simplex): HighsInterface.cpp lines 1713-1716 (1.16.0)
                if h16_ray_simplex {
                    debug_assert!(
                        self.info_.num_primal_infeasibilities > 0
                            || self.info_.num_semi_infeasibilities > 0
                    );
                } else {
                    debug_assert!(self.info_.num_primal_infeasibilities > 0);
                }
                debug_assert!(self.model_status_ == HighsModelStatus::kInfeasible);
            } else {
                // If someone has tried to get a dual ray for a feasible problem
                // - or if the relaxation is feasible - then any model and
                // primal KKT status of the original model has been lost
                self.info_.invalidate_primal_kkt();
                self.model_status_ = HighsModelStatus::kNotset;
            }
        }
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::getPrimalRayInterface (lines 1702-1823)
    /// HiGHS 1.16 (#3435): HighsInterface.cpp::Highs::getPrimalRayInterface (lines 1728-1854)
    pub(crate) fn get_primal_ray_interface(
        &mut self,
        has_primal_ray: &mut bool,
        primal_ray_value: Option<&mut [f64]>,
    ) -> HighsStatus {
        // HiGHS 1.16 #3435 (h16_ray_simplex): read once at entry
        let h16_ray_simplex = self.options_.h16_ray_simplex;
        let mut return_status = HighsStatus::kOk;
        let num_row = self.model_lp_.num_row_;
        let num_col = self.model_lp_.num_col_;
        // For an LP with no rows the primal ray is vacuous
        if num_row == 0 {
            return return_status;
        }
        if self.model_hessian_.dim_ != 0 {
            highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Cannot find primal ray for unbounded QP\n"
            );
            return HighsStatus::kError;
        }
        let mut has_invert = self.ekk_instance_.status_.has_invert;
        debug_assert!(!self.model_lp_.is_moved_);
        *has_primal_ray = self.ekk_instance_.primal_ray_record_.index != kNoRayIndex;

        let mut presolve = String::new();
        let mut solve_relaxation = false;
        // HiGHS 1.16 #3435 (h16_ray_simplex): HighsInterface.cpp lines 1745-1748 (1.16.0)
        let mut solver = String::new();
        let mut allow_unbounded_or_infeasible = false;
        let mut solve_unboundedness_problem = false;

        if let Some(primal_ray_value) = primal_ray_value {
            // User wants a primal ray whatever
            if !*has_primal_ray || !has_invert {
                // No primal ray is known, or no INVERT to compute it
                //
                // No point in trying to get a primal ray if the model status is
                // optimal
                if self.model_status_ == HighsModelStatus::kOptimal {
                    highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "Model status is optimal, so no primal ray is available\n"
                    );
                    return return_status;
                }
                highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Solving LP to try to compute primal ray\n"
                );
                if h16_ray_simplex {
                    // HiGHS 1.16 #3435 (h16_ray_simplex): HighsInterface.cpp lines 1761-1772 (1.16.0)
                    solve_relaxation = self.options_.solve_relaxation;
                    presolve = self.options_.presolve.clone();
                    solver = self.options_.solver.clone();
                    allow_unbounded_or_infeasible = self.options_.allow_unbounded_or_infeasible;
                } else {
                    presolve = self.get_option_value_str("presolve").unwrap_or_default();
                    solve_relaxation = self
                        .get_option_value_bool("solve_relaxation")
                        .unwrap_or(false);
                    allow_unbounded_or_infeasible = self
                        .get_option_value_bool("allow_unbounded_or_infeasible")
                        .unwrap_or(false);
                }
                solve_unboundedness_problem = true;
                let call_status = if h16_ray_simplex {
                    // HiGHS 1.16 #3435 (h16_ray_simplex): HighsInterface.cpp lines 1773-1783 (1.16.0)
                    self.options_.solve_relaxation = true;
                    self.options_.presolve = kHighsOffString.to_string();
                    self.options_.solver = kSimplexString.to_string();
                    self.options_.allow_unbounded_or_infeasible = false;
                    let call_status = self.run();
                    // Restore the relaxation, presolve and solver settings
                    self.options_.solve_relaxation = solve_relaxation;
                    self.options_.presolve = presolve.clone();
                    self.options_.solver = solver.clone();
                    self.options_.allow_unbounded_or_infeasible = allow_unbounded_or_infeasible;
                    call_status
                } else {
                    self.set_option_value_str("presolve", kHighsOffString);
                    self.set_option_value_bool("solve_relaxation", true);
                    self.set_option_value_bool("allow_unbounded_or_infeasible", false);
                    self.run()
                };
                if call_status != HighsStatus::kOk {
                    return_status = call_status;
                }
                *has_primal_ray = self.ekk_instance_.primal_ray_record_.index != kNoRayIndex;
                has_invert = self.ekk_instance_.status_.has_invert;
                debug_assert!(has_invert);
            }
            if *has_primal_ray {
                if !self.ekk_instance_.primal_ray_record_.value.is_empty() {
                    // Primal ray is already computed
                    highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "Copying known primal ray\n"
                    );
                    for i_col in 0..num_col {
                        primal_ray_value[i_col as usize] =
                            self.ekk_instance_.primal_ray_record_.value[i_col as usize];
                    }
                    return return_status;
                } else if has_invert {
                    // Primal ray is known and can be calculated
                    highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "Solving linear system to compute primal ray\n"
                    );
                    let col = self.ekk_instance_.primal_ray_record_.index;
                    debug_assert!(
                        self.ekk_instance_.basis_.nonbasic_flag_[col as usize] == kNonbasicFlagTrue
                    );
                    // Get this pivotal column
                    let mut rhs = vec![0.0f64; num_row as usize];
                    let mut column = vec![0.0f64; num_row as usize];
                    self.model_lp_.ensure_colwise();
                    let primal_ray_sign = self.ekk_instance_.primal_ray_record_.sign;
                    if col < num_col {
                        let mut i_el = self.model_lp_.a_matrix_.start_[col as usize];
                        while i_el < self.model_lp_.a_matrix_.start_[(col + 1) as usize] {
                            rhs[self.model_lp_.a_matrix_.index_[i_el as usize] as usize] =
                                primal_ray_sign as f64
                                    * self.model_lp_.a_matrix_.value_[i_el as usize];
                            i_el += 1;
                        }
                    } else {
                        rhs[(col - num_col) as usize] = primal_ray_sign as f64;
                    }
                    // C++: HighsInt* column_num_nz = 0 (a null pointer)
                    self.basis_solve_interface(&rhs, &mut column, None, None, false);
                    // Now zero primal_ray_value and scatter the column according to
                    // the basic variables.
                    for i_col in 0..num_col {
                        primal_ray_value[i_col as usize] = 0.0;
                    }
                    for i_row in 0..num_row {
                        let i_col = self.ekk_instance_.basis_.basic_index_[i_row as usize];
                        if i_col < num_col {
                            primal_ray_value[i_col as usize] = column[i_row as usize];
                        }
                    }
                    if col < num_col {
                        primal_ray_value[col as usize] = -(primal_ray_sign as f64);
                    }
                    // Now save the primal ray itself
                    self.ekk_instance_
                        .primal_ray_record_
                        .value
                        .resize(num_col as usize, 0.0);
                    for i_col in 0..num_col {
                        self.ekk_instance_.primal_ray_record_.value[i_col as usize] =
                            primal_ray_value[i_col as usize];
                    }
                }
            } else {
                highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "No primal ray found\n"
                );
                return_status = HighsStatus::kOk;
            }
        }
        let is_mip = self.model_lp_.is_mip();
        if solve_unboundedness_problem {
            if is_mip {
                // Unboundedness LP has been solved, but that will give dual
                // solution status kInfeasible which, for a MIP is not correct
                self.info_.dual_solution_status = kSolutionStatusNone;
                self.info_.invalidate_dual_kkt();
            }
            // HiGHS 1.16 #3435 (h16_ray_simplex): with the flag on, 1.16.0 has already
            // restored the options after run() (lines 1778-1783)
            if !h16_ray_simplex {
                // Restore the option values
                self.set_option_value_str("presolve", &presolve);
                self.set_option_value_bool("solve_relaxation", solve_relaxation);
                self.set_option_value_bool(
                    "allow_unbounded_or_infeasible",
                    allow_unbounded_or_infeasible,
                );
            }
            if *has_primal_ray {
                debug_assert!(is_mip || self.info_.num_dual_infeasibilities > 0);
                debug_assert!(self.model_status_ == HighsModelStatus::kUnbounded);
            }
        }
        return_status
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// x0 integer in [0, 4], x1 continuous in [0, 4], x0 + x1 >= 1.
    fn mip_highs() -> Highs {
        let mut highs = Highs::default();
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![4.0, 4.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        highs.model_lp_ = lp;
        highs
    }

    #[test]
    fn fixed_lp_fixes_discrete_columns() {
        let mut highs = mip_highs();
        let mut fixed = HighsLp::default();
        // No solution yet
        assert_eq!(highs.get_fixed_lp(&mut fixed), HighsStatus::kError);
        highs.solution_.col_value = vec![2.0, 0.5];
        highs.solution_.value_valid = true;
        assert_eq!(highs.get_fixed_lp(&mut fixed), HighsStatus::kOk);
        assert!(fixed.integrality_.is_empty());
        assert_eq!((fixed.col_lower_[0], fixed.col_upper_[0]), (2.0, 2.0));
        assert_eq!((fixed.col_lower_[1], fixed.col_upper_[1]), (0.0, 4.0));
        // A fractional discrete value is fixed with a warning
        highs.solution_.col_value[0] = 2.5;
        assert_eq!(highs.get_fixed_lp(&mut fixed), HighsStatus::kWarning);
    }

    #[test]
    fn rays_of_an_lp_without_rows_are_vacuous() {
        let mut highs = Highs::default();
        highs.model_lp_.num_col_ = 1;
        highs.model_lp_.col_cost_ = vec![1.0];
        highs.model_lp_.col_lower_ = vec![0.0];
        highs.model_lp_.col_upper_ = vec![1.0];
        highs.model_lp_.a_matrix_.start_ = vec![0, 0];
        let mut has = true;
        assert_eq!(highs.get_dual_ray(&mut has, None), HighsStatus::kOk);
        assert!(!has);
    }

    /// A quiet Highs with the HiGHS 1.16 #3435 option set as given and `solver` set.
    fn ray_highs(h16_ray_simplex: bool, solver: &str) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_ray_simplex", h16_ray_simplex),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", solver),
            HighsStatus::kOk
        );
        highs
    }

    /// One column and one row, column-wise (TestRays issue-3413 models).
    fn one_by_one_lp(cost: f64, col: (f64, f64), row: (f64, f64), value: f64) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![cost];
        lp.col_lower_ = vec![col.0];
        lp.col_upper_ = vec![col.1];
        lp.row_lower_ = vec![row.0];
        lp.row_upper_ = vec![row.1];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![value];
        lp
    }

    /// TestRays issue-3413, infeasible LP: x0 >= 0 and x0 <= -1, with IPM.
    /// The 1.16.0 nested solve forces simplex, so a valid dual ray results.
    #[test]
    fn h16_ray_simplex_dual_ray_ipm_infeasible() {
        let mut highs = ray_highs(true, kIpmString);
        let lp = one_by_one_lp(0.0, (0.0, kHighsInf), (-kHighsInf, -1.0), 1.0);
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs.run();
        let mut has_dual_ray = false;
        let mut dual_ray_value = vec![f64::NAN];
        assert_eq!(
            highs.get_dual_ray(&mut has_dual_ray, Some(&mut dual_ray_value)),
            HighsStatus::kOk
        );
        assert!(has_dual_ray);
        assert_eq!(dual_ray_value[0], -1.0);
    }

    /// TestRays issue-3413, unbounded LP: min -x0, x0 >= 0, x0 >= 1, with IPM.
    #[test]
    fn h16_ray_simplex_primal_ray_ipm_unbounded() {
        let mut highs = ray_highs(true, kIpmString);
        let lp = one_by_one_lp(-1.0, (0.0, kHighsInf), (1.0, kHighsInf), 1.0);
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs.run();
        let mut has_primal_ray = false;
        let mut primal_ray_value = vec![f64::NAN];
        assert_eq!(
            highs.get_primal_ray(&mut has_primal_ray, Some(&mut primal_ray_value)),
            HighsStatus::kOk
        );
        assert!(has_primal_ray);
        assert_eq!(primal_ray_value[0], 1.0);
    }

    /// Flag off keeps 1.15.1: the same IPM dual-ray request yields no ray.
    /// In a debug build the 1.15.1 `debug_assert!(has_invert)` after the
    /// nested IPM run fires instead, so this runs only without debug asserts.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "1.15.1 debug_assert!(has_invert) fires")]
    fn h16_ray_simplex_dual_ray_ipm_flag_off_keeps_1_15_1() {
        let mut highs = ray_highs(false, kIpmString);
        let lp = one_by_one_lp(0.0, (0.0, kHighsInf), (-kHighsInf, -1.0), 1.0);
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs.run();
        let mut has_dual_ray = true;
        let mut dual_ray_value = vec![f64::NAN];
        highs.get_dual_ray(&mut has_dual_ray, Some(&mut dual_ray_value));
        assert!(!has_dual_ray);
    }

    /// Rays-2415: an infeasible MIP (integer x in [0, 1] with x = 0.5). The
    /// primal solution status must stay kSolutionStatusNone after the ray
    /// request, with the flag both off and on.
    #[test]
    fn h16_ray_simplex_rays_2415_both_settings() {
        for h16 in [false, true] {
            let mut highs = ray_highs(h16, kSimplexString);
            let mut lp = one_by_one_lp(0.0, (0.0, 1.0), (0.5, 0.5), 1.0);
            lp.integrality_ = vec![HighsVarType::kInteger];
            assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
            assert_eq!(highs.run(), HighsStatus::kOk);
            assert_ne!(
                highs.get_info().primal_solution_status,
                kSolutionStatusFeasible
            );
            assert_eq!(highs.get_info().primal_solution_status, kSolutionStatusNone);
            let mut has_dual_ray = false;
            let mut dual_ray_value = vec![f64::NAN];
            assert_eq!(
                highs.get_dual_ray(&mut has_dual_ray, Some(&mut dual_ray_value)),
                HighsStatus::kOk
            );
            assert_ne!(
                highs.get_info().primal_solution_status,
                kSolutionStatusFeasible
            );
            assert_eq!(highs.get_info().primal_solution_status, kSolutionStatusNone);
        }
    }

    /// TestRays issue-3413, infeasible problem with a semi-integer column:
    /// x0 in {0} U [2, 3], x1 = 4, x0 - x1 >= 0. The dual ray is (1).
    #[test]
    fn h16_ray_simplex_semi_integer_dual_ray() {
        let mut highs = ray_highs(true, "choose");
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![2.0, 4.0];
        lp.col_upper_ = vec![3.0, 4.0];
        lp.row_lower_ = vec![0.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.integrality_ = vec![HighsVarType::kSemiInteger, HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, -1.0];
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs.run();
        let mut has_dual_ray = false;
        let mut dual_ray_value = vec![f64::NAN, f64::NAN];
        assert_eq!(
            highs.get_dual_ray(&mut has_dual_ray, Some(&mut dual_ray_value)),
            HighsStatus::kOk
        );
        assert!(has_dual_ray);
        assert_eq!(dual_ray_value[0], 1.0);
    }

    /// TestRays issue-3413, unbounded problem with a semi-integer column:
    /// min -x0 - x1, x0 in {0} U [2, 3], x1 >= 0, x1 >= 1. The primal ray is (0, 1).
    #[test]
    fn h16_ray_simplex_semi_integer_primal_ray() {
        let mut highs = ray_highs(true, "choose");
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![-1.0, -1.0];
        lp.col_lower_ = vec![2.0, 0.0];
        lp.col_upper_ = vec![3.0, kHighsInf];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.integrality_ = vec![HighsVarType::kSemiInteger, HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs.run();
        let mut has_primal_ray = false;
        let mut primal_ray_value = vec![f64::NAN, f64::NAN];
        assert_eq!(
            highs.get_primal_ray(&mut has_primal_ray, Some(&mut primal_ray_value)),
            HighsStatus::kOk
        );
        assert!(has_primal_ray);
        assert_eq!(primal_ray_value[0], 0.0);
        assert_eq!(primal_ray_value[1], 1.0);
    }
}
