//! HiGHS source: lp_data/HighsSolution.cpp (subset; getEffectiveCosts from 1.16.0),
//! lp_data/HighsInfo.cpp (invalidate*)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The solution/info bookkeeping that `solveLp`/`solveLpSimplex` and
//! `Highs::calledOptimizeModel` call. The KKT-failure machinery
//! (`getLpKktFailures`, `reportKktFailures`) is out of scope.

use crate::lp_data::highs_complete_solution::assess_col_primal_solution_h16;
use crate::prelude::*;

/// HiGHS: HighsSolution.cpp::resetModelStatusAndHighsInfo(HighsLpSolverObject&) (lines 1975-1978, 4 lines)
pub fn reset_model_status_and_highs_info_so(solver_object: &mut HighsLpSolverObject) {
    reset_model_status_and_highs_info(&mut solver_object.model_status_, solver_object.highs_info_);
}

/// HiGHS: HighsSolution.cpp::resetModelStatusAndHighsInfo(HighsModelStatus&, HighsInfo&) (lines 1980-1987, 8 lines)
pub fn reset_model_status_and_highs_info(
    model_status: &mut HighsModelStatus,
    highs_info: &mut HighsInfo,
) {
    *model_status = HighsModelStatus::kNotset;
    highs_info.objective_function_value = 0.0;
    highs_info.primal_solution_status = kSolutionStatusNone;
    highs_info.dual_solution_status = kSolutionStatusNone;
    highs_info.invalidate_kkt();
}

/// HiGHS: HighsSolution.cpp::isColPrimalSolutionRightSize (lines 2004-2007, 4 lines)
pub fn is_col_primal_solution_right_size(lp: &HighsLp, solution: &HighsSolution) -> bool {
    solution.col_value.len() == lp.num_col_ as usize
}

/// HiGHS: HighsSolution.cpp::isRowPrimalSolutionRightSize (lines 2009-2012, 4 lines)
pub fn is_row_primal_solution_right_size(lp: &HighsLp, solution: &HighsSolution) -> bool {
    solution.row_value.len() == lp.num_row_ as usize
}

/// HiGHS: HighsSolution.cpp::isPrimalSolutionRightSize (lines 2014-2018, 5 lines)
pub fn is_primal_solution_right_size(lp: &HighsLp, solution: &HighsSolution) -> bool {
    is_col_primal_solution_right_size(lp, solution)
        && is_row_primal_solution_right_size(lp, solution)
}

/// HiGHS: HighsSolution.cpp::isColDualSolutionRightSize (lines 2020-2023, 4 lines)
pub fn is_col_dual_solution_right_size(lp: &HighsLp, solution: &HighsSolution) -> bool {
    solution.col_dual.len() == lp.num_col_ as usize
}

/// HiGHS: HighsSolution.cpp::isRowDualSolutionRightSize (lines 2025-2028, 4 lines)
pub fn is_row_dual_solution_right_size(lp: &HighsLp, solution: &HighsSolution) -> bool {
    solution.row_dual.len() == lp.num_row_ as usize
}

/// HiGHS: HighsSolution.cpp::isDualSolutionRightSize (lines 2030-2033, 4 lines)
pub fn is_dual_solution_right_size(lp: &HighsLp, solution: &HighsSolution) -> bool {
    is_col_dual_solution_right_size(lp, solution) && is_row_dual_solution_right_size(lp, solution)
}

/// HiGHS: HighsSolution.cpp::isSolutionRightSize (lines 2035-2038, 4 lines)
pub fn is_solution_right_size(lp: &HighsLp, solution: &HighsSolution) -> bool {
    is_primal_solution_right_size(lp, solution) && is_dual_solution_right_size(lp, solution)
}

/// HiGHS: HighsSolution.cpp::isBasisRightSize (lines 2040-2043, 4 lines)
pub fn is_basis_right_size(lp: &HighsLp, basis: &HighsBasis) -> bool {
    basis.col_status.len() == lp.num_col_ as usize && basis.row_status.len() == lp.num_row_ as usize
}

impl HighsInfo {
    /// HiGHS: HighsInfo.cpp::HighsInfo::invalidate (lines 17-32, 16 lines)
    pub fn invalidate(&mut self) {
        self.valid = false;
        self.mip_node_count = -1;
        self.simplex_iteration_count = -1;
        self.ipm_iteration_count = -1;
        self.crossover_iteration_count = -1;
        self.pdlp_iteration_count = -1;
        self.qp_iteration_count = -1;
        self.basis_validity = kBasisValidityInvalid;
        self.objective_function_value = 0.0;
        self.mip_dual_bound = 0.0;
        self.mip_gap = kHighsInf;
        self.max_integrality_violation = kHighsIllegalInfeasibilityMeasure;
        self.invalidate_kkt();
        self.primal_dual_integral = -kHighsInf;
    }

    /// HiGHS: HighsInfo.cpp::HighsInfo::invalidateKkt (lines 34-37, 4 lines)
    pub fn invalidate_kkt(&mut self) {
        self.invalidate_primal_kkt();
        self.invalidate_dual_kkt();
    }

    /// HiGHS: HighsInfo.cpp::HighsInfo::invalidatePrimalKkt (lines 39-56, 18 lines)
    pub fn invalidate_primal_kkt(&mut self) {
        self.primal_solution_status = kSolutionStatusNone;
        self.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
        self.num_semi_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.max_semi_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.sum_semi_infeasibilities = kHighsIllegalInfeasibilityMeasure;
        self.num_relative_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.max_relative_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.num_primal_residual_errors = kHighsIllegalResidualCount;
        self.max_primal_residual_error = kHighsIllegalResidualMeasure;
        self.num_relative_primal_residual_errors = kHighsIllegalResidualCount;
        self.max_relative_primal_residual_error = kHighsIllegalResidualMeasure;
        self.num_complementarity_violations = kHighsIllegalComplementarityCount;
        self.max_complementarity_violation = kHighsIllegalComplementarityViolation;
        self.primal_dual_objective_error = kHighsIllegalComplementarityViolation;
    }

    /// HiGHS: HighsInfo.cpp::HighsInfo::invalidateDualKkt (lines 58-72, 15 lines)
    pub fn invalidate_dual_kkt(&mut self) {
        self.dual_solution_status = kSolutionStatusNone;
        self.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
        self.num_relative_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.max_relative_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.num_dual_residual_errors = kHighsIllegalResidualCount;
        self.max_dual_residual_error = kHighsIllegalResidualMeasure;
        self.num_relative_dual_residual_errors = kHighsIllegalResidualCount;
        self.max_relative_dual_residual_error = kHighsIllegalResidualMeasure;
        self.num_complementarity_violations = kHighsIllegalComplementarityCount;
        self.max_complementarity_violation = kHighsIllegalComplementarityViolation;
        self.primal_dual_objective_error = kHighsIllegalComplementarityViolation;
    }
}

/// HiGHS 1.16 (#3556): HighsSolution.cpp::assessUserColValues (lines 2277-2378).
/// Called only when `h16_user_solution_assess` is on. `index == nullptr` ->
/// `None`; `std::vector<HighsBool> is_set` -> `Vec<bool>`;
/// `highsVarTypeToString` (HighsLpUtils.cpp lines 3745-3760) is inlined as a
/// match, as in highs_pass_model_raw.rs. The C++ default argument
/// `perform_feasibility_check = true` is explicit at every call site.
pub fn assess_user_col_values(
    lp: &HighsLp,
    options: &HighsOptions,
    num_entries: HighsInt,
    index: Option<&[HighsInt]>,
    value: &[f64],
    perform_feasibility_check: bool,
) -> HighsStatus {
    let packed = index.is_some();
    if packed {
        if num_entries < 0 {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "assessUserColValues: packed solution has num_entries = {} < 0\n",
                num_entries
            );
            return HighsStatus::kError;
        }
        if num_entries == 0 {
            // No solution provided, which suggests a user error so, rather
            // than going on to define a completely undefined start
            // solution, issue a warning amd return
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kWarning,
                "assessUserColValues: packed solution has num_entries = 0, \
                 so no user solution is defined\n"
            );
            return HighsStatus::kWarning;
        }
    } else {
        // lp.num_col_ entries of value will be used, so num_entries
        // should reflect this
        if num_entries != lp.num_col_ {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "assessUserColValues: full solution has num_entries ({}) \
                 not equal to num_col ({})\n",
                num_entries,
                lp.num_col_
            );
            return HighsStatus::kError;
        }
    }
    let mut return_status = HighsStatus::kOk;
    let to_ix = if packed { num_entries } else { lp.num_col_ };
    let mut num_duplicates: HighsInt = 0;
    let mut is_set: Vec<bool> = vec![false; lp.num_col_ as usize];
    let mut num_infeasible_value: HighsInt = 0;
    let max_warning: HighsInt = 10;
    let has_integrality = !lp.integrality_.is_empty();
    let assess_as_mip = lp.is_mip() && !options.solve_relaxation;
    let k_primal_feasibility_tolerance = if assess_as_mip {
        options.mip_feasibility_tolerance
    } else {
        options.primal_feasibility_tolerance
    };
    for i_x in 0..to_ix {
        let i_col = match index {
            Some(index) => index[i_x as usize],
            None => i_x,
        };
        if i_col < 0 || i_col >= lp.num_col_ {
            debug_assert!(packed);
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "assessUserColValues: Entry {} of user solution has index {} out of \
                 range [0, {})\n",
                i_x,
                i_col,
                lp.num_col_
            );
            return HighsStatus::kError;
        }
        // Count duplicates, and record that this value has been set, so
        // future duplicates can be spotted
        if is_set[i_col as usize] {
            num_duplicates += 1;
        }
        is_set[i_col as usize] = true;
        // Don't test undefined values for feasibility!
        if !perform_feasibility_check || value[i_x as usize] == kHighsUndefined {
            continue;
        }
        let mut col_infeasibility: f64 = 0.0;
        let mut integer_infeasibility: f64 = 0.0;
        let type_ = if has_integrality {
            lp.integrality_[i_col as usize]
        } else {
            HighsVarType::kContinuous
        };
        assess_col_primal_solution_h16(
            assess_as_mip,
            options,
            value[i_x as usize],
            lp.col_lower_[i_col as usize],
            lp.col_upper_[i_col as usize],
            type_,
            &mut col_infeasibility,
            &mut integer_infeasibility,
        );
        if col_infeasibility > k_primal_feasibility_tolerance
            || integer_infeasibility > options.mip_feasibility_tolerance
        {
            num_infeasible_value += 1;
            if num_infeasible_value <= max_warning {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kWarning,
                    "assessUserColValues: Entry {} of user solution has value \
                     {}: infeasible \
                     for bounds [{}, {}]{}{}\n",
                    i_x,
                    crate::util::c_printf::c_fmt_g(value[i_x as usize], 6),
                    crate::util::c_printf::c_fmt_g(lp.col_lower_[i_col as usize], 6),
                    crate::util::c_printf::c_fmt_g(lp.col_upper_[i_col as usize], 6),
                    if has_integrality { " and type " } else { "" },
                    if has_integrality {
                        match type_ {
                            HighsVarType::kContinuous => "continuous",
                            HighsVarType::kInteger => "integer",
                            HighsVarType::kSemiContinuous => "semi continuous",
                            HighsVarType::kSemiInteger => "semi integer",
                            HighsVarType::kImplicitInteger => "implicit integer",
                        }
                    } else {
                        ""
                    }
                );
            }
            return_status = HighsStatus::kWarning;
        }
    }
    if num_duplicates > 0 {
        debug_assert!(packed);
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "assessUserColValues: Set of user solution indices has {} \
             duplicate{}: last \
             value used\n",
            num_duplicates,
            if num_duplicates == 1 { "" } else { "s" }
        );
        return_status = HighsStatus::kWarning;
    }
    if num_infeasible_value > max_warning {
        debug_assert!(return_status == HighsStatus::kWarning);
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "assessUserColValues: {} entries of user solution are infeasible\n",
            num_infeasible_value
        );
        return_status = HighsStatus::kWarning;
    }
    return_status
}

/// HiGHS 1.16 (#3556): HighsSolution.cpp::assignUserColValues (lines 2380-2392).
/// Called only when `h16_user_solution_assess` is on. `index == nullptr` ->
/// `None`.
pub fn assign_user_col_values(
    num_col: HighsInt,
    col_value: &mut Vec<f64>,
    num_entries: HighsInt,
    index: Option<&[HighsInt]>,
    value: &[f64],
) {
    debug_assert!(num_col >= 0);
    let packed = index.is_some();
    if !packed {
        debug_assert!(num_entries == num_col);
    }
    let to_ix = if packed { num_entries } else { num_col };
    col_value.clear();
    col_value.resize(num_col as usize, kHighsUndefined);
    for i_x in 0..to_ix {
        let i_col = match index {
            Some(index) => index[i_x as usize],
            None => i_x,
        };
        col_value[i_col as usize] = value[i_x as usize];
    }
}

/// HiGHS 1.16 (h16_scaling_report): HighsSolution.cpp::getEffectiveCosts (lines 685-781, 1.16.0)
/// C++: `std::vector<double> getEffectiveCosts(const HighsLp& lp, const HighsOptions& options)`
///
/// Only called from the 1.16.0 coefficient-range report. `options` is unused, as in C++.
/// The 1.16.0 `getKktFailures` call (HighsSolution.cpp line 169) is not ported.
pub fn get_effective_costs(lp: &HighsLp, _options: &HighsOptions) -> Vec<f64> {
    // Look for free column singletons that allow equations to be added
    // into the objective and, hence, compute the effective costs of an
    // LP: handles the case where the objective is min z, with a
    // constraint z - c^Tx so the "true" costs are c.
    //
    // Constants to distinguish row status: kFree (row may contain a
    // free column singleton); >0 (position in "ID" vectors to indicate
    // that its singleton column is to be added into costs); kUsed (row
    // contains a free column singleton already added into costs);
    const K_USED: HighsInt = -2;
    const K_FREE: HighsInt = -1;
    let mut effective_costs: Vec<f64> = lp.col_cost_.clone();
    let mut col_count: Vec<HighsInt> = Vec::new();
    let mut row_id: Vec<HighsInt> = vec![K_FREE; lp.num_row_ as usize];
    let mut row_of_id: Vec<HighsInt> = Vec::new();
    let mut col_of_id: Vec<HighsInt> = Vec::new();
    let mut row_mu_of_id: Vec<f64> = Vec::new();

    debug_assert!(lp.a_matrix_.is_colwise());
    let start = &lp.a_matrix_.start_;
    let index = &lp.a_matrix_.index_;
    let value = &lp.a_matrix_.value_;
    col_count.reserve(lp.num_col_ as usize);
    for i_col in 0..lp.num_col_ as usize {
        col_count.push(start[i_col + 1] - start[i_col]);
    }
    loop {
        col_of_id.clear();
        row_of_id.clear();
        row_mu_of_id.clear();
        let mut num_id: HighsInt = 0;
        for i_col in 0..lp.num_col_ as usize {
            let free_column_singleton = col_count[i_col] == 1
                && lp.col_lower_[i_col] == -kHighsInf
                && lp.col_upper_[i_col] == kHighsInf
                && effective_costs[i_col] != 0.0;
            if !free_column_singleton {
                continue;
            }
            // Potential free column singleton - so long as the row is an
            // equation and has not already been selected for substitution
            // into the objective
            //
            // Find the row in the column (which may originally have had
            // more than one entry)
            let mut i_el = start[i_col];
            let mut i_row: HighsInt = -1;
            while i_el < start[i_col + 1] {
                let this_row = index[i_el as usize];
                if row_id[this_row as usize] == K_FREE {
                    i_row = this_row;
                    break;
                }
                i_el += 1;
            }
            if i_row < 0 {
                continue;
            }
            debug_assert!(i_row == index[i_el as usize]);
            // Row must be an equation
            if lp.row_lower_[i_row as usize] != lp.row_upper_[i_row as usize] {
                continue;
            }
            // Found a free column singleton, so can zero its column cost,
            // as its effective cost is zero after substitution, so the
            // elimination need not be done
            col_count[i_col] = 0;
            // Retain the index of the free column singleton, its
            // corresponding row, the multiplier for the elimination, and
            // row_id[iRow]: this is both a marker to indicate that the row
            // cannot be chosen for a subsequent free column singleton, and
            // the position in col_of_id/row_of_id/row_mu_of_id so that
            // these values can be extracted when passing through the matrix
            col_of_id.push(i_col as HighsInt);
            row_of_id.push(i_row);
            row_mu_of_id.push(effective_costs[i_col] / value[i_el as usize]);
            row_id[i_row as usize] = num_id;
            num_id += 1;
        }
        if num_id == 0 {
            break;
        }
        // Now pass through the matrix, adding multiples of selected rows
        // into the objective
        for i_col in 0..lp.num_col_ as usize {
            if col_count[i_col] == 0 {
                continue;
            }
            for i_el in start[i_col] as usize..start[i_col + 1] as usize {
                let id = row_id[index[i_el] as usize];
                if id < 0 {
                    continue;
                }
                // Entry in a row to be added into effective_costs: get the
                // multiplier, reduce the column's count and update its
                // effective cost
                let row_mu = row_mu_of_id[id as usize];
                col_count[i_col] -= 1;
                effective_costs[i_col] -= row_mu * value[i_el];
            }
        }
        // Update row_id so that these rows cannot be chosen again, and
        // zero the effective cost for the free column singletons
        for id in 0..num_id as usize {
            let i_col = col_of_id[id];
            let i_row = row_of_id[id];
            row_id[i_row as usize] = K_USED;
            effective_costs[i_col as usize] = 0.0;
        }
    }
    effective_costs
}
