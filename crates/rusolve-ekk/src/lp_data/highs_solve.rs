//! HiGHS source: lp_data/HighsSolve.cpp (solveLp, solveUnconstrainedLp,
//! assessExcessiveObjectiveBoundScaling; its 1.16.0 version behind `h16_scaling_report`),
//! lp_data/HighsModelUtils.cpp (highsStatusFromHighsModelStatus)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `solveLp` is ported for its simplex, IPX, HiPO (`crate::ipm::hipo`)
//! and PDLP (cuPDLP-C and HiPDLP, `crate::pdlp`) paths. `useHipo` returns
//! `HighsExternalApi::isAvailable<HighsExtras::hipo>()` for `"ipm"`/`"hipo"`,
//! modelled by the non-upstream switch `HighsOptions::hipo_available`
//! (default `false`). With the switch off `"ipm"`/`"hipo"` route to IPX
//! exactly as in the IPX oracle build (`ipx_header.hipo_available = 0` in
//! the `highs-oracle ipx` fixtures); the HiPO oracle build has HiPO and the
//! HiPO tests set the switch. `debugHighsLpSolution` returns
//! `kNotChecked` at the default `highs_debug_level`, so it is omitted.

#![allow(
    clippy::int_plus_one,
    clippy::needless_late_init,
    clippy::collapsible_match,
    clippy::assertions_on_constants
)]

use crate::ipm::hipo_wrapper::solve_lp_hipo_so;
use crate::ipm::ipx_wrapper::solve_lp_ipx_so;
use crate::lp_data::highs_solution::{
    is_solution_right_size, reset_model_status_and_highs_info, reset_model_status_and_highs_info_so,
};
use crate::lp_data::hstruct::HighsUserScaleData;
use crate::lp_data::lp_assess::assess_lp;
use crate::pdlp::cupdlp_wrapper::solve_lp_cupdlp_so;
use crate::pdlp::hipdlp_wrapper::solve_lp_hipdlp_so;
use crate::prelude::*;
use crate::simplex::happ::solve_lp_simplex;
use crate::util::c_printf::{c_fmt_e, c_pad_right};

/// HiGHS: HighsSolve.cpp::solveLp (lines 20-176, 157 lines); simplex, IPX and PDLP paths.
pub fn solve_lp(
    ekk: &mut HEkk,
    solver_object: &mut HighsLpSolverObject,
    message: &str,
) -> HighsStatus {
    let mut return_status = HighsStatus::kOk;
    let mut call_status: HighsStatus;
    // Reset unscaled model status and solution params - except for
    // iteration counts
    reset_model_status_and_highs_info_so(solver_object);
    crate::highs_log_user!(
        &solver_object.options_.log_options,
        HighsLogType::kInfo,
        "{}\n",
        message
    );
    if solver_object.options_.highs_debug_level > kHighsDebugLevelMin {
        // Shouldn't have to check validity of the LP since this is done when it is
        // loaded or modified
        let options = solver_object.options_.clone();
        call_status = assess_lp(solver_object.lp_, &options);
        // If any errors have been found or normalisation carried out,
        // call_status will be ERROR or WARNING, so only valid return is OK.
        debug_assert!(call_status == HighsStatus::kOk);
        return_status = interpret_call_status(
            &solver_object.options_.log_options,
            call_status,
            return_status,
            "assessLp",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
    }
    let use_only_ipm =
        use_ipm(&solver_object.options_.solver) || solver_object.options_.run_centring;
    let use_hipo = use_hipo(
        &*solver_object.options_,
        kSolverString,
        &*solver_object.lp_,
        false,
    );
    let use_ipx = use_only_ipm && !use_hipo;
    if solver_object.lp_.num_row_ == 0 || solver_object.lp_.a_matrix_.num_nz() == 0 {
        // LP is unconstrained due to having no rows or a zero constraint
        // matrix, so solve directly
        call_status = solve_unconstrained_lp_so(solver_object);
        return_status = interpret_call_status(
            &solver_object.options_.log_options,
            call_status,
            return_status,
            "solveUnconstrainedLp",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
    } else if use_only_ipm || use_pdlp(&solver_object.options_.solver) {
        // Use IPM or PDLP
        if use_only_ipm {
            // Use IPM to solve the LP
            if use_hipo {
                // Use HIPO to solve the LP. The `try` / `catch` is as for
                // IPX below.
                call_status = if crate::ipm::hipo_wrapper::solve_lp_hipo_dualized(solver_object) {
                    HighsStatus::kOk
                } else {
                    solve_lp_hipo_so(solver_object)
                };
                return_status = interpret_call_status(
                    &solver_object.options_.log_options,
                    call_status,
                    return_status,
                    "solveLpHipo",
                );
            } else if use_ipx {
                // C++ `try { ... } catch (const std::exception&)`: IPX throws
                // become panics.
                call_status = solve_lp_ipx_so(solver_object);
                return_status = interpret_call_status(
                    &solver_object.options_.log_options,
                    call_status,
                    return_status,
                    "solveLpIpx",
                );
            }
        } else {
            // Use cuPDLP-C or HiPDLP to solve the LP. The
            // `profiling->start/stop(kSubSolverPdlp)` calls are dropped
            // (`profiling_` is out of scope); the `try` / `catch` is as for
            // IPX.
            if solver_object.options_.solver == kPdlpString {
                call_status = solve_lp_cupdlp_so(solver_object);
            } else {
                call_status = solve_lp_hipdlp_so(solver_object);
            }
            return_status = interpret_call_status(
                &solver_object.options_.log_options,
                call_status,
                return_status,
                "solveLp-Pdlp",
            );
        }
        // rusolve: an IPM that ends in an error (HiPO internal
        // error, IPX failure) hands the LP to the simplex instead of giving up
        if return_status == HighsStatus::kError
            && use_only_ipm
            && solver_object.options_.lp_ipm_failure_policy
            && solver_object.model_status_ != HighsModelStatus::kTimeLimit
        {
            crate::highs_log_user!(
                &solver_object.options_.log_options,
                HighsLogType::kWarning,
                "IPM failed, so solve with simplex (lp_ipm_failure_policy)\n"
            );
            solver_object.basis_.valid = false;
            return simplex_solve(ekk, solver_object);
        }
        // Check for error return
        if return_status == HighsStatus::kError {
            return return_status;
        }

        // Non-error return requires a primal solution
        debug_assert!(solver_object.solution_.value_valid);

        if use_ipm(&solver_object.options_.solver) || solver_object.options_.run_centring {
            // Setting the IPM-specific values of (highs_)info_ has been done in
            // solveLpHipo/Ipx
            let unwelcome_ipx_status = solver_object.model_status_ == HighsModelStatus::kUnknown
                || (solver_object.model_status_ == HighsModelStatus::kUnboundedOrInfeasible
                    && !solver_object.options_.allow_unbounded_or_infeasible);
            if unwelcome_ipx_status {
                // When performing an analytic centre calculation, the setting
                // of options.run_crossover is ignored, so simplex clean-up is
                // not possible - or desirable, anyway!
                crate::highs_log_user!(
                    &solver_object.options_.log_options,
                    HighsLogType::kWarning,
                    "Unwelcome IPX status of {}: basis is {}valid; solution is {}valid; run_crossover is \"{}\"\n",
                    util_model_status_to_string(solver_object.model_status_),
                    if solver_object.basis_.valid {
                        ""
                    } else {
                        "not "
                    },
                    if solver_object.solution_.value_valid {
                        ""
                    } else {
                        "not "
                    },
                    if solver_object.options_.run_centring {
                        kHighsOffString
                    } else {
                        solver_object.options_.run_crossover.as_str()
                    }
                );
                let allow_simplex_cleanup = solver_object.options_.run_crossover != kHighsOffString
                    && !solver_object.options_.run_centring;
                if allow_simplex_cleanup {
                    // IPX has returned a model status that HiGHS would rather
                    // avoid, so perform simplex clean-up if crossover was allowed.
                    crate::highs_log_user!(
                        &solver_object.options_.log_options,
                        HighsLogType::kWarning,
                        "IPM solution is imprecise, so clean up with simplex\n"
                    );
                    return_status = simplex_solve(ekk, solver_object);
                    if return_status == HighsStatus::kError {
                        return return_status;
                    }
                }
            }
        }
    } else {
        // Use Simplex
        return_status = simplex_solve(ekk, solver_object);
        if return_status == HighsStatus::kError {
            return return_status;
        }
    }
    // Analyse the HiGHS (basic) solution: debugHighsLpSolution returns
    // kNotChecked below kHighsDebugLevelCheap
    return_status
}

/// HiGHS: HighsSolve.cpp::solveLp, the `simplexSolve` lambda (lines 47-59).
/// It assigns the enclosing `return_status` and returns it; every C++ call
/// is `return_status = simplexSolve();`, so returning the value is the same.
fn simplex_solve(ekk: &mut HEkk, solver_object: &mut HighsLpSolverObject) -> HighsStatus {
    let mut return_status = HighsStatus::kOk;
    let call_status = solve_lp_simplex(ekk, solver_object);
    return_status = interpret_call_status(
        &solver_object.options_.log_options,
        call_status,
        return_status,
        "solveLpSimplex",
    );
    if return_status == HighsStatus::kError {
        return return_status;
    }
    if !is_solution_right_size(solver_object.lp_, solver_object.solution_) {
        crate::highs_log_user!(
            &solver_object.options_.log_options,
            HighsLogType::kError,
            "Inconsistent solution returned from solver\n"
        );
        return_status = HighsStatus::kError;
    }
    return_status
}

/// HiGHS: HighsSolve.cpp::solveUnconstrainedLp(HighsLpSolverObject&) (lines 180-185, 6 lines)
pub fn solve_unconstrained_lp_so(solver_object: &mut HighsLpSolverObject) -> HighsStatus {
    solve_unconstrained_lp(
        solver_object.options_,
        solver_object.lp_,
        &mut solver_object.model_status_,
        solver_object.highs_info_,
        solver_object.solution_,
        solver_object.basis_,
    )
}

/// HiGHS: HighsSolve.cpp::solveUnconstrainedLp(const HighsOptions&, ...) (lines 189-371, 183 lines)
pub fn solve_unconstrained_lp(
    options: &HighsOptions,
    lp: &HighsLp,
    model_status: &mut HighsModelStatus,
    highs_info: &mut HighsInfo,
    solution: &mut HighsSolution,
    basis: &mut HighsBasis,
) -> HighsStatus {
    // Aliase to model status and solution parameters
    reset_model_status_and_highs_info(model_status, highs_info);

    // Check that the LP really is unconstrained!
    debug_assert!(lp.num_row_ == 0 || lp.a_matrix_.num_nz() == 0);
    if lp.num_row_ > 0 {
        // LP has rows, but should only be here if the constraint matrix
        // is zero
        if lp.a_matrix_.num_nz() > 0 {
            return HighsStatus::kError;
        }
    }

    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kInfo,
        "Solving an unconstrained LP with {} columns\n",
        lp.num_col_
    );
    let nc = lp.num_col_ as usize;
    solution.col_value.clear();
    solution.col_value.resize(nc, 0.0);
    solution.col_dual.clear();
    solution.col_dual.resize(nc, 0.0);
    basis.col_status.clear();
    basis.col_status.resize(nc, HighsBasisStatus::kNonbasic);
    // No rows for primal solution, dual solution or basis
    solution.row_value.clear();
    solution.row_dual.clear();
    basis.row_status.clear();

    let primal_feasibility_tolerance = options.primal_feasibility_tolerance;
    let dual_feasibility_tolerance = options.dual_feasibility_tolerance;

    // Initialise the objective value calculation. Done using
    // HighsSolution so offset is vanilla
    let mut objective = lp.offset_;

    highs_info.num_primal_infeasibilities = 0;
    highs_info.max_primal_infeasibility = 0.0;
    highs_info.sum_primal_infeasibilities = 0.0;
    highs_info.num_dual_infeasibilities = 0;
    highs_info.max_dual_infeasibility = 0.0;
    highs_info.sum_dual_infeasibilities = 0.0;

    if lp.num_row_ > 0 {
        // Assign primal, dual and basis status for rows, checking for
        // infeasibility
        for i_row in 0..lp.num_row_ as usize {
            let mut primal_infeasibility = 0.0;
            let lower = lp.row_lower_[i_row];
            let upper = lp.row_upper_[i_row];
            if lower > primal_feasibility_tolerance {
                // Lower bound too large for zero activity
                primal_infeasibility = lower;
            } else if upper < -primal_feasibility_tolerance {
                // Upper bound too small for zero activity
                primal_infeasibility = -upper;
            }
            solution.row_value.push(0.0);
            solution.row_dual.push(0.0);
            basis.row_status.push(HighsBasisStatus::kBasic);
            if primal_infeasibility > primal_feasibility_tolerance {
                highs_info.num_primal_infeasibilities += 1;
            }
            highs_info.sum_primal_infeasibilities += primal_infeasibility;
            highs_info.max_primal_infeasibility =
                cpp_max(primal_infeasibility, highs_info.max_primal_infeasibility);
        }
    }

    let sense = lp.sense_ as HighsInt;
    for i_col in 0..nc {
        let cost = lp.col_cost_[i_col];
        let dual = (sense as f64) * cost;
        let lower = lp.col_lower_[i_col];
        let upper = lp.col_upper_[i_col];
        let value: f64;
        let mut primal_infeasibility = 0.0;
        let dual_infeasibility: f64;
        let status: HighsBasisStatus;
        if lower > upper {
            // Inconsistent bounds, so set the variable to lower bound,
            // unless it's infinite. Otherwise set the variable to upper
            // bound, unless it's infinite. Otherwise set the variable to
            // zero.
            if highs_is_infinity(lower) {
                // Lower bound of +inf
                if highs_is_infinity(-upper) {
                    // Upper bound of -inf
                    value = 0.0;
                    status = HighsBasisStatus::kZero;
                    primal_infeasibility = kHighsInf;
                    dual_infeasibility = dual.abs();
                } else {
                    // Finite upper bound - since lower exceeds it
                    value = upper;
                    status = HighsBasisStatus::kUpper;
                    primal_infeasibility = lower - value;
                    dual_infeasibility = cpp_max(dual, 0.0);
                }
            } else {
                // Finite lower bound
                value = lower;
                status = HighsBasisStatus::kLower;
                primal_infeasibility = value - upper;
                dual_infeasibility = cpp_max(-dual, 0.0);
            }
        } else if highs_is_infinity(-lower) && highs_is_infinity(upper) {
            // Free column: set to zero and record dual infeasibility
            value = 0.0;
            status = HighsBasisStatus::kZero;
            dual_infeasibility = dual.abs();
        } else if dual >= dual_feasibility_tolerance {
            // Column with sufficiently positive dual
            if !highs_is_infinity(-lower) {
                // Set to this finite lower bound
                value = lower;
                status = HighsBasisStatus::kLower;
                dual_infeasibility = 0.0;
            } else {
                // Infinite lower bound so set to upper bound and record dual
                // infeasibility
                value = upper;
                status = HighsBasisStatus::kUpper;
                dual_infeasibility = dual;
            }
        } else if dual <= -dual_feasibility_tolerance {
            // Column with sufficiently negative dual
            if !highs_is_infinity(upper) {
                // Set to this finite upper bound
                value = upper;
                status = HighsBasisStatus::kUpper;
                dual_infeasibility = 0.0;
            } else {
                // Infinite upper bound so set to lower bound and record dual
                // infeasibility
                value = lower;
                status = HighsBasisStatus::kLower;
                dual_infeasibility = -dual;
            }
        } else {
            // Column with sufficiently small dual: set to lower bound (if
            // finite) otherwise upper bound
            if highs_is_infinity(-lower) {
                value = upper;
                status = HighsBasisStatus::kUpper;
            } else {
                value = lower;
                status = HighsBasisStatus::kLower;
            }
            dual_infeasibility = dual.abs();
        }
        debug_assert!(status != HighsBasisStatus::kNonbasic);
        debug_assert!(dual_infeasibility >= 0.0);
        solution.col_value[i_col] = value;
        solution.col_dual[i_col] = (sense as f64) * dual;
        basis.col_status[i_col] = status;
        objective += value * cost;
        if primal_infeasibility > primal_feasibility_tolerance {
            highs_info.num_primal_infeasibilities += 1;
        }
        highs_info.sum_primal_infeasibilities += primal_infeasibility;
        highs_info.max_primal_infeasibility =
            cpp_max(primal_infeasibility, highs_info.max_primal_infeasibility);
        if dual_infeasibility > dual_feasibility_tolerance {
            highs_info.num_dual_infeasibilities += 1;
        }
        highs_info.sum_dual_infeasibilities += dual_infeasibility;
        highs_info.max_dual_infeasibility =
            cpp_max(dual_infeasibility, highs_info.max_dual_infeasibility);
    }
    highs_info.objective_function_value = objective;
    solution.value_valid = true;
    solution.dual_valid = true;
    basis.valid = true;
    basis.useful = true;
    highs_info.basis_validity = kBasisValidityValid;
    set_solution_status(highs_info);
    if highs_info.num_primal_infeasibilities != 0 {
        // Primal infeasible
        *model_status = HighsModelStatus::kInfeasible;
    } else if highs_info.num_dual_infeasibilities != 0 {
        // Dual infeasible => primal unbounded for unconstrained LP
        *model_status = HighsModelStatus::kUnbounded;
    } else {
        *model_status = HighsModelStatus::kOptimal;
    }

    HighsStatus::kOk
}

/// HiGHS: HighsModelUtils.cpp::highsStatusFromHighsModelStatus (lines 1530-1575, 46 lines)
pub fn highs_status_from_highs_model_status(model_status: HighsModelStatus) -> HighsStatus {
    match model_status {
        HighsModelStatus::kNotset => HighsStatus::kError,
        HighsModelStatus::kLoadError => HighsStatus::kError,
        HighsModelStatus::kModelError => HighsStatus::kError,
        HighsModelStatus::kPresolveError => HighsStatus::kError,
        HighsModelStatus::kSolveError => HighsStatus::kError,
        HighsModelStatus::kPostsolveError => HighsStatus::kError,
        HighsModelStatus::kMemoryLimit => HighsStatus::kError,
        HighsModelStatus::kModelEmpty => HighsStatus::kOk,
        HighsModelStatus::kOptimal => HighsStatus::kOk,
        HighsModelStatus::kInfeasible => HighsStatus::kOk,
        HighsModelStatus::kUnboundedOrInfeasible => HighsStatus::kOk,
        HighsModelStatus::kUnbounded => HighsStatus::kOk,
        HighsModelStatus::kObjectiveBound => HighsStatus::kOk,
        HighsModelStatus::kObjectiveTarget => HighsStatus::kOk,
        HighsModelStatus::kTimeLimit => HighsStatus::kWarning,
        HighsModelStatus::kIterationLimit => HighsStatus::kWarning,
        HighsModelStatus::kSolutionLimit => HighsStatus::kWarning,
        HighsModelStatus::kInterrupt => HighsStatus::kWarning,
        HighsModelStatus::kHighsInterrupt => HighsStatus::kWarning,
        HighsModelStatus::kUnknown => HighsStatus::kWarning,
        #[allow(unreachable_patterns)]
        _ => HighsStatus::kError,
    }
}

/// HiGHS: HighsSolve.cpp::assessExcessiveObjectiveBoundScaling (lines 377-676)
/// C++: `void assessExcessiveObjectiveBoundScaling(const HighsLogOptions log_options, const HighsModel& model, HighsUserScaleData& user_scale_data)`
///
/// PORT NOTE: the port's model is the LP alone (no `HighsHessian`), so `model.hessian_`
/// has no nonzeros: the Hessian extremes keep their initial `kHighsInf` / `-kHighsInf`.
/// The port passes the whole `HighsOptions` (the 1.16.0 signature) so that
/// `h16_scaling_report` can select the 1.16.0 report
/// (`assess_excessive_objective_bound_scaling_h16`); with it off only
/// `options.log_options` is read, as in 1.15.1.
pub fn assess_excessive_objective_bound_scaling(
    options: &HighsOptions,
    lp: &HighsLp,
    user_scale_data: &mut HighsUserScaleData,
) {
    if options.h16_scaling_report {
        assess_excessive_objective_bound_scaling_h16(options, lp, user_scale_data);
        return;
    }
    let log_options = &options.log_options;
    if lp.num_col_ == 0 || lp.num_row_ == 0 {
        return;
    }
    let user_cost_or_bound_scale =
        user_scale_data.user_objective_scale != 0 || user_scale_data.user_bound_scale != 0;
    let small_objective_coefficient = kExcessivelySmallObjectiveCoefficient;
    let large_objective_coefficient = kExcessivelyLargeObjectiveCoefficient;
    let small_bound = kExcessivelySmallBoundValue;
    let large_bound = kExcessivelyLargeBoundValue;
    let mut message = String::new();
    if user_cost_or_bound_scale {
        if user_scale_data.user_objective_scale != 0 {
            message += &format!(
                " user_objective_scale option value of {}",
                user_scale_data.user_objective_scale
            );
        }
        if user_scale_data.user_bound_scale != 0 {
            if user_scale_data.user_objective_scale != 0 {
                message += " and";
            }
            message += &format!(
                " user_bound_scale option value of {}",
                user_scale_data.user_bound_scale
            );
        }
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "Assessing costs and bounds after applying{}\n",
            message
        );
    }
    // Lambda for assessing a finite nonzero
    let assess_finite_nonzero = |value: f64, min_value: &mut f64, max_value: &mut f64| {
        let abs_value = value.abs();
        if abs_value > 0.0 && abs_value < kHighsInf {
            *min_value = cpp_min(abs_value, *min_value);
            *max_value = cpp_max(abs_value, *max_value);
        }
    };
    let mut min_continuous_col_cost = kHighsInf;
    let mut min_noncontinuous_col_cost = kHighsInf;
    let mut max_continuous_col_cost = -kHighsInf;
    let mut max_noncontinuous_col_cost = -kHighsInf;
    let mut min_continuous_col_bound = kHighsInf;
    let mut min_noncontinuous_col_bound = kHighsInf;
    let mut max_continuous_col_bound = -kHighsInf;
    let mut max_noncontinuous_col_bound = -kHighsInf;
    let is_mip = !lp.integrality_.is_empty();
    for i_col in 0..lp.num_col_ as usize {
        if is_mip && lp.integrality_[i_col] != HighsVarType::kContinuous {
            assess_finite_nonzero(
                lp.col_cost_[i_col],
                &mut min_noncontinuous_col_cost,
                &mut max_noncontinuous_col_cost,
            );
            assess_finite_nonzero(
                lp.col_lower_[i_col],
                &mut min_noncontinuous_col_bound,
                &mut max_noncontinuous_col_bound,
            );
            assess_finite_nonzero(
                lp.col_upper_[i_col],
                &mut min_noncontinuous_col_bound,
                &mut max_noncontinuous_col_bound,
            );
        } else {
            assess_finite_nonzero(
                lp.col_cost_[i_col],
                &mut min_continuous_col_cost,
                &mut max_continuous_col_cost,
            );
            assess_finite_nonzero(
                lp.col_lower_[i_col],
                &mut min_continuous_col_bound,
                &mut max_continuous_col_bound,
            );
            assess_finite_nonzero(
                lp.col_upper_[i_col],
                &mut min_continuous_col_bound,
                &mut max_continuous_col_bound,
            );
        }
    }
    let mut min_col_cost = cpp_min(min_continuous_col_cost, min_noncontinuous_col_cost);
    let mut max_col_cost = cpp_max(max_continuous_col_cost, max_noncontinuous_col_cost);
    let mut min_col_bound = cpp_min(min_continuous_col_bound, min_noncontinuous_col_bound);
    let mut max_col_bound = cpp_max(max_continuous_col_bound, max_noncontinuous_col_bound);

    let mut min_matrix_value = kHighsInf;
    let mut max_matrix_value = -kHighsInf;
    let num_matrix_nz = lp.a_matrix_.num_nz();
    for i_el in 0..num_matrix_nz as usize {
        assess_finite_nonzero(
            lp.a_matrix_.value_[i_el],
            &mut min_matrix_value,
            &mut max_matrix_value,
        );
    }

    let mut min_row_bound = kHighsInf;
    let mut max_row_bound = -kHighsInf;
    for i_row in 0..lp.num_row_ as usize {
        assess_finite_nonzero(lp.row_lower_[i_row], &mut min_row_bound, &mut max_row_bound);
        assess_finite_nonzero(lp.row_upper_[i_row], &mut min_row_bound, &mut max_row_bound);
    }

    let min_continuous_hessian_value = kHighsInf;
    let max_continuous_hessian_value = -kHighsInf;
    let num_hessian_nz: HighsInt = 0;

    let mut min_scalable_bound = cpp_min(min_continuous_col_bound, min_row_bound);
    let mut max_scalable_bound = cpp_max(max_continuous_col_bound, max_row_bound);
    if min_scalable_bound == kHighsInf {
        min_scalable_bound = 0.0;
    }
    if max_scalable_bound == -kHighsInf {
        max_scalable_bound = 0.0;
    }
    let _ = min_scalable_bound;

    if min_col_cost == kHighsInf {
        min_col_cost = 0.0;
    }
    if max_col_cost == -kHighsInf {
        max_col_cost = 0.0;
    }
    if min_col_bound == kHighsInf {
        min_col_bound = 0.0;
    }
    if max_col_bound == -kHighsInf {
        max_col_bound = 0.0;
    }
    if min_row_bound == kHighsInf {
        min_row_bound = 0.0;
    }
    if max_row_bound == -kHighsInf {
        max_row_bound = 0.0;
    }

    let mut min_hessian_value = min_continuous_hessian_value;
    let mut max_hessian_value = max_continuous_hessian_value;
    if min_hessian_value == kHighsInf {
        min_hessian_value = 0.0;
    }
    if max_hessian_value == -kHighsInf {
        max_hessian_value = 0.0;
    }

    // C `%5.0e`
    let e50 = |x: f64| c_pad_right(&c_fmt_e(x, 0), 5);
    // Report on the coefficient ranges
    crate::highs_log_user!(log_options, HighsLogType::kInfo, "Coefficient ranges:\n");
    if num_matrix_nz != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "  Matrix  [{}, {}]\n",
            e50(min_matrix_value),
            e50(max_matrix_value)
        );
    }
    if lp.num_col_ != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "  Cost    [{}, {}]\n",
            e50(min_col_cost),
            e50(max_col_cost)
        );
        if num_hessian_nz != 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "  Hessian [{}, {}]\n",
                e50(min_hessian_value),
                e50(max_hessian_value)
            );
        }
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "  Bound   [{}, {}]\n",
            e50(min_col_bound),
            e50(max_col_bound)
        );
    }
    if lp.num_row_ != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "  RHS     [{}, {}]\n",
            e50(min_row_bound),
            e50(max_row_bound)
        );
    }

    debug_assert!(max_col_cost >= 0.0);
    debug_assert!(max_col_bound >= 0.0);
    debug_assert!(max_row_bound >= 0.0);

    let problem = if user_cost_or_bound_scale {
        "User-scaled problem"
    } else {
        "Problem"
    };

    let warn = |what: &str| {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kWarning,
            "{} has some excessively {}\n",
            problem,
            what
        );
    };
    if 0.0 < min_col_cost && min_col_cost < small_objective_coefficient {
        warn("small costs");
    }
    if max_col_cost > large_objective_coefficient {
        warn("large costs");
    }
    if 0.0 < min_hessian_value && min_hessian_value < small_objective_coefficient {
        warn("small Hessian values");
    }
    if max_hessian_value > large_objective_coefficient {
        warn("large Hessian values");
    }
    if 0.0 < min_col_bound && min_col_bound < small_bound {
        warn("small column bounds");
    }
    if max_col_bound > large_bound {
        warn("large column bounds");
    }
    if 0.0 < min_row_bound && min_row_bound < small_bound {
        warn("small row bounds");
    }
    if max_row_bound > large_bound {
        warn("large row bounds");
    }

    // Lambda to determine recommended user scaling values
    let suggest_scaling = |_min_value: f64, max_value: f64, small_value: f64, large_value: f64| {
        let mut ratio = 1.0;
        if max_value > large_value {
            ratio = large_value / max_value;
        } else if 0.0 < max_value && max_value < small_value {
            ratio = small_value / max_value;
        }
        debug_assert!(ratio != 0.0);
        ratio
    };

    // Lambda to for outward rounding of log
    let outer_rounded_log = |value: f64, base: HighsInt| -> HighsInt {
        debug_assert!(value > 0.0);
        debug_assert!(base == 2 || base == 10);
        // C++ assigns the double floor/ceil to a HighsInt
        if base == 2 {
            (if value < 1.0 {
                value.log2().floor()
            } else {
                value.log2().ceil()
            }) as HighsInt
        } else {
            (if value < 1.0 {
                value.log10().floor()
            } else {
                value.log10().ceil()
            }) as HighsInt
        }
    };
    let suggested_bound_scaling = suggest_scaling(
        min_scalable_bound,
        max_scalable_bound,
        small_bound,
        large_bound,
    );
    let dl_user_bound_scale = outer_rounded_log(suggested_bound_scaling, 2);
    user_scale_data.suggested_user_bound_scale =
        user_scale_data.user_bound_scale + dl_user_bound_scale;
    let suggested_bound_scale_order_of_magnitude = outer_rounded_log(suggested_bound_scaling, 10);
    let suggested_user_bound_scale_value =
        2.0f64.powf(user_scale_data.suggested_user_bound_scale as f64);
    min_noncontinuous_col_cost *= suggested_user_bound_scale_value;
    max_noncontinuous_col_cost *= suggested_user_bound_scale_value;
    min_hessian_value /= suggested_user_bound_scale_value;
    max_hessian_value /= suggested_user_bound_scale_value;
    let _ = (min_hessian_value, max_hessian_value);

    min_col_cost = cpp_min(min_continuous_col_cost, min_noncontinuous_col_cost);
    max_col_cost = cpp_max(max_continuous_col_cost, max_noncontinuous_col_cost);
    let mut min_objective_coefficient = cpp_min(min_col_cost, min_continuous_hessian_value);
    let mut max_objective_coefficient = cpp_max(max_col_cost, max_continuous_hessian_value);
    if min_objective_coefficient == kHighsInf {
        min_objective_coefficient = 0.0;
    }
    if max_objective_coefficient == -kHighsInf {
        max_objective_coefficient = 0.0;
    }

    let suggested_objective_scaling = suggest_scaling(
        min_objective_coefficient,
        max_objective_coefficient,
        small_objective_coefficient,
        large_objective_coefficient,
    );
    let dl_user_objective_scale = outer_rounded_log(suggested_objective_scaling, 2);
    user_scale_data.suggested_user_objective_scale =
        user_scale_data.user_objective_scale + dl_user_objective_scale;
    let suggested_objective_scale_order_of_magnitude =
        outer_rounded_log(suggested_objective_scaling, 10);

    let mut order_of_magnitude_message = suggested_objective_scale_order_of_magnitude != 0
        && user_scale_data.user_objective_scale == 0;
    message.clear();
    if order_of_magnitude_message {
        message += &format!(
            "   Consider scaling the objective by 1e{:+}",
            suggested_objective_scale_order_of_magnitude
        );
    }
    if dl_user_objective_scale != 0 {
        if !order_of_magnitude_message {
            message += "   Consider";
        } else {
            message += ", or";
        }
        message += &format!(
            " setting the user_objective_scale option to {}",
            user_scale_data.suggested_user_objective_scale
        );
    }
    if order_of_magnitude_message || dl_user_objective_scale != 0 {
        crate::highs_log_user!(log_options, HighsLogType::kWarning, "{}\n", message);
    }

    message.clear();
    order_of_magnitude_message =
        suggested_bound_scale_order_of_magnitude != 0 && user_scale_data.user_bound_scale == 0;
    if order_of_magnitude_message {
        message += &format!(
            "   Consider scaling the    bounds by 1e{:+}",
            suggested_bound_scale_order_of_magnitude
        );
    }
    if dl_user_bound_scale != 0 {
        if !order_of_magnitude_message {
            message += "   Consider";
        } else {
            message += ", or";
        }
        message += &format!(
            " setting the user_bound_scale option to {}",
            user_scale_data.suggested_user_bound_scale
        );
    }
    if order_of_magnitude_message || dl_user_bound_scale != 0 {
        crate::highs_log_user!(log_options, HighsLogType::kWarning, "{}\n", message);
    }
}

/// HiGHS 1.16 (h16_scaling_report): HighsSolve.cpp::assessExcessiveObjectiveBoundScaling
/// (lines 376-805, 1.16.0)
/// C++: `void assessExcessiveObjectiveBoundScaling(const HighsOptions& options, const HighsModel& model, HighsUserScaleData& user_scale_data)`
///
/// The 1.16.0 report: effective cost ranges ([`get_effective_costs`]), bound ranges split
/// into continuous and non-continuous variables, the reworded warnings, objective scaling
/// suggested from the effective costs, and the tolerance-ratio warnings. The caller has
/// made the matrix column-wise (`ensureColwise` at Highs.cpp lines 1238 and 2443).
///
/// PORT NOTE: as in the 1.15.1 port, the model is the LP alone, so the Hessian extremes
/// keep their initial `kHighsInf` / `-kHighsInf` and `num_hessian_nz` is 0.
fn assess_excessive_objective_bound_scaling_h16(
    options: &HighsOptions,
    lp: &HighsLp,
    user_scale_data: &mut HighsUserScaleData,
) {
    use crate::lp_data::highs_solution::get_effective_costs;
    use crate::util::c_printf::c_fmt_g;
    let log_options = &options.log_options;
    if lp.num_col_ == 0 || lp.num_row_ == 0 {
        return;
    }
    let effective_cost = get_effective_costs(lp, options);
    let user_cost_or_bound_scale =
        user_scale_data.user_objective_scale != 0 || user_scale_data.user_bound_scale != 0;
    let small_objective_coefficient = kExcessivelySmallObjectiveCoefficient;
    let large_objective_coefficient = kExcessivelyLargeObjectiveCoefficient;
    let small_bound = kExcessivelySmallBoundValue;
    let large_bound = kExcessivelyLargeBoundValue;
    let large_integer_bound = kExcessivelyLargeIntegerBoundValue;
    let mut message = String::new();
    if user_cost_or_bound_scale {
        if user_scale_data.user_objective_scale != 0 {
            message += &format!(
                " user_objective_scale option value of {}",
                user_scale_data.user_objective_scale
            );
        }
        if user_scale_data.user_bound_scale != 0 {
            if user_scale_data.user_objective_scale != 0 {
                message += " and";
            }
            message += &format!(
                " user_bound_scale option value of {}",
                user_scale_data.user_bound_scale
            );
        }
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "Assessing costs and bounds after applying{}\n",
            message
        );
    }
    // Lambda for assessing a finite nonzero
    let assess_finite_nonzero = |value: f64, min_value: &mut f64, max_value: &mut f64| {
        let abs_value = value.abs();
        if abs_value > 0.0 && abs_value < kHighsInf {
            *min_value = cpp_min(abs_value, *min_value);
            *max_value = cpp_max(abs_value, *max_value);
        }
    };
    let mut min_continuous_col_cost = kHighsInf;
    let mut min_noncontinuous_col_cost = kHighsInf;
    let mut max_continuous_col_cost = -kHighsInf;
    let mut max_noncontinuous_col_cost = -kHighsInf;
    let mut min_effective_continuous_col_cost = kHighsInf;
    let mut min_effective_noncontinuous_col_cost = kHighsInf;
    let mut max_effective_continuous_col_cost = -kHighsInf;
    let mut max_effective_noncontinuous_col_cost = -kHighsInf;
    let mut min_continuous_col_bound = kHighsInf;
    let mut min_noncontinuous_col_bound = kHighsInf;
    let mut max_continuous_col_bound = -kHighsInf;
    let mut max_noncontinuous_col_bound = -kHighsInf;
    let is_mip = !lp.integrality_.is_empty();
    let mut num_continuous_variable: HighsInt = 0;
    let mut num_noncontinuous_variable: HighsInt = 0;
    for i_col in 0..lp.num_col_ as usize {
        if is_mip && lp.integrality_[i_col] != HighsVarType::kContinuous {
            assess_finite_nonzero(
                lp.col_cost_[i_col],
                &mut min_noncontinuous_col_cost,
                &mut max_noncontinuous_col_cost,
            );
            assess_finite_nonzero(
                effective_cost[i_col],
                &mut min_effective_noncontinuous_col_cost,
                &mut max_effective_noncontinuous_col_cost,
            );
            assess_finite_nonzero(
                lp.col_lower_[i_col],
                &mut min_noncontinuous_col_bound,
                &mut max_noncontinuous_col_bound,
            );
            assess_finite_nonzero(
                lp.col_upper_[i_col],
                &mut min_noncontinuous_col_bound,
                &mut max_noncontinuous_col_bound,
            );
            num_noncontinuous_variable += 1;
        } else {
            assess_finite_nonzero(
                lp.col_cost_[i_col],
                &mut min_continuous_col_cost,
                &mut max_continuous_col_cost,
            );
            assess_finite_nonzero(
                effective_cost[i_col],
                &mut min_effective_continuous_col_cost,
                &mut max_effective_continuous_col_cost,
            );
            assess_finite_nonzero(
                lp.col_lower_[i_col],
                &mut min_continuous_col_bound,
                &mut max_continuous_col_bound,
            );
            assess_finite_nonzero(
                lp.col_upper_[i_col],
                &mut min_continuous_col_bound,
                &mut max_continuous_col_bound,
            );
            num_continuous_variable += 1;
        }
    }
    let mut min_col_cost = cpp_min(min_continuous_col_cost, min_noncontinuous_col_cost);
    let mut max_col_cost = cpp_max(max_continuous_col_cost, max_noncontinuous_col_cost);

    let mut min_effective_col_cost = cpp_min(
        min_effective_continuous_col_cost,
        min_effective_noncontinuous_col_cost,
    );
    let mut max_effective_col_cost = cpp_max(
        max_effective_continuous_col_cost,
        max_effective_noncontinuous_col_cost,
    );

    let mut min_matrix_value = kHighsInf;
    let mut max_matrix_value = -kHighsInf;
    let num_matrix_nz = lp.a_matrix_.num_nz();
    for i_el in 0..num_matrix_nz as usize {
        assess_finite_nonzero(
            lp.a_matrix_.value_[i_el],
            &mut min_matrix_value,
            &mut max_matrix_value,
        );
    }

    let mut min_row_bound = kHighsInf;
    let mut max_row_bound = -kHighsInf;
    for i_row in 0..lp.num_row_ as usize {
        assess_finite_nonzero(lp.row_lower_[i_row], &mut min_row_bound, &mut max_row_bound);
        assess_finite_nonzero(lp.row_upper_[i_row], &mut min_row_bound, &mut max_row_bound);
    }

    let min_continuous_hessian_value = kHighsInf;
    let max_continuous_hessian_value = -kHighsInf;
    let num_hessian_nz: HighsInt = 0;

    // Determine the minimum and maximum overall bounds that can be
    // scaled with user_bound_scale before zeroing extrema due to
    // absence of finite nonzero bounds

    let mut min_scalable_bound = cpp_min(min_continuous_col_bound, min_row_bound);
    let mut max_scalable_bound = cpp_max(max_continuous_col_bound, max_row_bound);
    if min_scalable_bound == kHighsInf {
        min_scalable_bound = 0.0;
    }
    if max_scalable_bound == -kHighsInf {
        max_scalable_bound = 0.0;
    }
    let _ = min_scalable_bound;

    if min_col_cost == kHighsInf {
        min_col_cost = 0.0;
    }
    if max_col_cost == -kHighsInf {
        max_col_cost = 0.0;
    }
    if min_effective_col_cost == kHighsInf {
        min_effective_col_cost = 0.0;
    }
    if max_effective_col_cost == -kHighsInf {
        max_effective_col_cost = 0.0;
    }
    if min_continuous_col_bound == kHighsInf {
        min_continuous_col_bound = 0.0;
    }
    if max_continuous_col_bound == -kHighsInf {
        max_continuous_col_bound = 0.0;
    }
    if min_noncontinuous_col_bound == kHighsInf {
        min_noncontinuous_col_bound = 0.0;
    }
    if max_noncontinuous_col_bound == -kHighsInf {
        max_noncontinuous_col_bound = 0.0;
    }
    if min_row_bound == kHighsInf {
        min_row_bound = 0.0;
    }
    if max_row_bound == -kHighsInf {
        max_row_bound = 0.0;
    }

    let mut min_hessian_value = min_continuous_hessian_value;
    let mut max_hessian_value = max_continuous_hessian_value;
    if min_hessian_value == kHighsInf {
        min_hessian_value = 0.0;
    }
    if max_hessian_value == -kHighsInf {
        max_hessian_value = 0.0;
    }

    // C `%5.0e`
    let e50 = |x: f64| c_pad_right(&c_fmt_e(x, 0), 5);
    // Report on the coefficient ranges
    crate::highs_log_user!(log_options, HighsLogType::kInfo, "Coefficient ranges:\n");
    if num_matrix_nz != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "  Matrix  [{}, {}]\n",
            e50(min_matrix_value),
            e50(max_matrix_value)
        );
    }
    if lp.num_col_ != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "  Cost    [{}, {}]\n",
            e50(min_col_cost),
            e50(max_col_cost)
        );
        if min_effective_col_cost != min_col_cost || max_effective_col_cost != max_col_cost {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "  Cost    [{}, {}] (effective)\n",
                e50(min_effective_col_cost),
                e50(max_effective_col_cost)
            );
        }
        if num_hessian_nz != 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "  Hessian [{}, {}]\n",
                e50(min_hessian_value),
                e50(max_hessian_value)
            );
        }
        if num_continuous_variable != 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "  Bound   [{}, {}]{}\n",
                e50(min_continuous_col_bound),
                e50(max_continuous_col_bound),
                if num_noncontinuous_variable > 0 {
                    " (continuous)"
                } else {
                    ""
                }
            );
        }
        if num_noncontinuous_variable != 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "  Bound   [{}, {}]{}\n",
                e50(min_noncontinuous_col_bound),
                e50(max_noncontinuous_col_bound),
                if num_continuous_variable > 0 {
                    " (non-continuous)"
                } else {
                    ""
                }
            );
        }
    }
    if lp.num_row_ != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "  RHS     [{}, {}]\n",
            e50(min_row_bound),
            e50(max_row_bound)
        );
    }

    // LPs with no columns or no finite nonzero costs will have
    // max_col_cost = 0
    debug_assert!(max_col_cost >= 0.0);
    // LPs with no columns or no finite nonzero bounds will have
    // max_continuous_col_bound = 0 and max_noncontinuous_col_bound = 0
    debug_assert!(max_continuous_col_bound >= 0.0);
    debug_assert!(max_noncontinuous_col_bound >= 0.0);
    // LPs with no rows or no finite nonzero bounds will have
    // max_row_bound = 0
    debug_assert!(max_row_bound >= 0.0);

    let problem = if user_cost_or_bound_scale {
        "User-scaled problem"
    } else {
        "Problem"
    };

    let warn = |what: &str| {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kWarning,
            "{} has some excessively {}\n",
            problem,
            what
        );
    };
    // If the effective costs aren't small or large, don't report costs
    // that are small or large
    if 0.0 < min_col_cost
        && min_col_cost < small_objective_coefficient
        && 0.0 < min_effective_col_cost
        && min_effective_col_cost < small_objective_coefficient
    {
        warn("small costs");
    }
    if max_col_cost > large_objective_coefficient
        && max_effective_col_cost > large_objective_coefficient
    {
        warn("large costs");
    }

    // If extreme effective costs aren't different from extreme costs,
    // don't report extreme costs
    if 0.0 < min_effective_col_cost
        && min_effective_col_cost != min_col_cost
        && min_effective_col_cost < small_objective_coefficient
    {
        warn("small effective costs");
    }
    if max_effective_col_cost != max_col_cost
        && max_effective_col_cost > large_objective_coefficient
    {
        warn("large effective costs");
    }

    if 0.0 < min_hessian_value && min_hessian_value < small_objective_coefficient {
        warn("small Hessian values");
    }
    if max_hessian_value > large_objective_coefficient {
        warn("large Hessian values");
    }
    let continuous = if is_mip { " continuous" } else { "" };
    if 0.0 < min_continuous_col_bound && min_continuous_col_bound < small_bound {
        warn(&format!("small bounds on{} variables", continuous));
    }
    if max_continuous_col_bound > large_bound {
        warn(&format!("large bounds on{} variables", continuous));
    }
    if 0.0 < min_noncontinuous_col_bound && min_noncontinuous_col_bound < small_bound {
        warn("small bounds on non-continuous variables");
    }
    if max_noncontinuous_col_bound > large_integer_bound {
        warn("large bounds on non-continuous variables");
    }
    if 0.0 < min_row_bound && min_row_bound < small_bound {
        warn("small bounds on constraints");
    }
    if max_row_bound > large_bound {
        warn("large bounds on constraints");
    }

    // Lambda to determine recommended user scaling values
    let suggest_scaling = |_min_value: f64, max_value: f64, small_value: f64, large_value: f64| {
        let mut ratio = 1.0;
        if max_value > large_value {
            // Max scalable value is large, so suggest scaling values down
            // so that the max value is large_value
            ratio = large_value / max_value;
        } else if 0.0 < max_value && max_value < small_value {
            // All scalable values are small, so suggest scaling them up so
            // the max value is small_value
            ratio = small_value / max_value;
        }
        debug_assert!(ratio != 0.0);
        ratio
    };

    // Lambda to for outward rounding of log
    let outer_rounded_log = |value: f64, base: HighsInt| -> HighsInt {
        debug_assert!(value > 0.0);
        debug_assert!(base == 2 || base == 10);
        // C++ assigns the double floor/ceil to a HighsInt
        if base == 2 {
            (if value < 1.0 {
                value.log2().floor()
            } else {
                value.log2().ceil()
            }) as HighsInt
        } else {
            (if value < 1.0 {
                value.log10().floor()
            } else {
                value.log10().ceil()
            }) as HighsInt
        }
    };
    let suggested_bound_scaling = suggest_scaling(
        min_scalable_bound,
        max_scalable_bound,
        small_bound,
        large_bound,
    );
    // Determine the suggested (new) value for user_bound_scale,
    // allowing for the fact that the current value has been applied
    let dl_user_bound_scale = outer_rounded_log(suggested_bound_scaling, 2);
    user_scale_data.suggested_user_bound_scale =
        user_scale_data.user_bound_scale + dl_user_bound_scale;
    // Determine the order of magnitude of the suggested bound scaling -
    // just for logging
    let suggested_bound_scale_order_of_magnitude = outer_rounded_log(suggested_bound_scaling, 10);
    // Applying the suggested bound scaling requires the costs and
    // matrix columns of non-continuous variables to be scaled, and any
    // Hessian entries are also scaled
    //
    // Determine the corresponding extreme non-continuous costs and
    // update the extreme costs so that objective scaling can be
    // suggested
    let suggested_user_bound_scale_value =
        2.0f64.powf(user_scale_data.suggested_user_bound_scale as f64);
    min_effective_noncontinuous_col_cost *= suggested_user_bound_scale_value;
    max_effective_noncontinuous_col_cost *= suggested_user_bound_scale_value;
    min_hessian_value /= suggested_user_bound_scale_value;
    max_hessian_value /= suggested_user_bound_scale_value;
    let _ = min_hessian_value;

    min_effective_col_cost = cpp_min(
        min_effective_continuous_col_cost,
        min_effective_noncontinuous_col_cost,
    );
    max_effective_col_cost = cpp_max(
        max_effective_continuous_col_cost,
        max_effective_noncontinuous_col_cost,
    );
    let mut min_objective_coefficient =
        cpp_min(min_effective_col_cost, min_continuous_hessian_value);
    let mut max_objective_coefficient =
        cpp_max(max_effective_col_cost, max_continuous_hessian_value);
    if min_objective_coefficient == kHighsInf {
        min_objective_coefficient = 0.0;
    }
    if max_objective_coefficient == -kHighsInf {
        max_objective_coefficient = 0.0;
    }

    let suggested_objective_scaling = suggest_scaling(
        min_objective_coefficient,
        max_objective_coefficient,
        small_objective_coefficient,
        large_objective_coefficient,
    );
    // Determine the suggested (new) value for user_objective_scale,
    // allowing for the fact that the current value has been applied
    let dl_user_objective_scale = outer_rounded_log(suggested_objective_scaling, 2);
    user_scale_data.suggested_user_objective_scale =
        user_scale_data.user_objective_scale + dl_user_objective_scale;
    // Determine the order of magnitude of the suggested objective scaling -
    // just for logging
    let suggested_objective_scale_order_of_magnitude =
        outer_rounded_log(suggested_objective_scaling, 10);

    let mut warning_issued = false;
    // Only report the order of magnitude scaling if there is no user
    // scaling
    let mut order_of_magnitude_message = suggested_objective_scale_order_of_magnitude != 0
        && user_scale_data.user_objective_scale == 0;
    message.clear();
    if order_of_magnitude_message {
        message += &format!(
            "   Consider scaling the objective by 1e{:+}",
            suggested_objective_scale_order_of_magnitude
        );
    }
    if dl_user_objective_scale != 0 {
        if !order_of_magnitude_message {
            message += "   Consider";
        } else {
            message += ", or";
        }
        message += &format!(
            " setting the user_objective_scale option to {}",
            user_scale_data.suggested_user_objective_scale
        );
    }
    if order_of_magnitude_message || dl_user_objective_scale != 0 {
        crate::highs_log_user!(log_options, HighsLogType::kWarning, "{}\n", message);
        warning_issued = true;
    }

    message.clear();
    order_of_magnitude_message =
        suggested_bound_scale_order_of_magnitude != 0 && user_scale_data.user_bound_scale == 0;
    if order_of_magnitude_message {
        message += &format!(
            "   Consider scaling the    bounds by 1e{:+}",
            suggested_bound_scale_order_of_magnitude
        );
    }
    if dl_user_bound_scale != 0 {
        if !order_of_magnitude_message {
            message += "   Consider";
        } else {
            message += ", or";
        }
        message += &format!(
            " setting the user_bound_scale option to {}",
            user_scale_data.suggested_user_bound_scale
        );
    }
    if order_of_magnitude_message || dl_user_bound_scale != 0 {
        crate::highs_log_user!(log_options, HighsLogType::kWarning, "{}\n", message);
        warning_issued = true;
    }

    // Assess the reationship between extreme objective/bound values and
    // feasibility tolerances
    let use_mip_feasibility_tolerances = lp.is_mip() && !options.solve_relaxation;
    // HighsOptions::getMipPrimalFeasibilityTolerance / getMipDualFeasibilityTolerance
    // (HighsOptions.h lines 1795-1800, 1.16.0)
    let use_primal_feasibility_tolerance = if use_mip_feasibility_tolerances {
        options.mip_feasibility_tolerance
    } else {
        options.primal_feasibility_tolerance
    };
    let use_dual_feasibility_tolerance = if use_mip_feasibility_tolerances {
        options.mip_feasibility_tolerance * 0.1
    } else {
        options.dual_feasibility_tolerance
    };
    let max_objective = cpp_max(max_effective_col_cost, max_hessian_value);
    let max_bound = cpp_max(
        max_continuous_col_bound,
        cpp_max(max_noncontinuous_col_bound, max_row_bound),
    );
    let primal_ratio = max_bound / use_primal_feasibility_tolerance;
    let dual_ratio = max_objective / use_dual_feasibility_tolerance;
    let primal_fragile = primal_ratio > 1e-1 / kHighsMacheps;
    let dual_fragile = dual_ratio > 1e-1 / kHighsMacheps;
    let numerically_fragile = primal_fragile || dual_fragile;
    // C `%.1g`
    let g1 = |x: f64| c_fmt_g(x, 1);
    if numerically_fragile {
        let blanks = " ".repeat(problem.len());
        if primal_fragile {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kWarning,
                "{} has excessively large ratio ({}) between max variable or constraint bound ({}) and primal feasibility tolerance ({})\n",
                problem,
                g1(primal_ratio),
                g1(max_bound),
                g1(use_primal_feasibility_tolerance)
            );
        }
        if dual_fragile {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kWarning,
                "{} has excessively large ratio ({}) between max objective coefficient ({}) and dual feasibility tolerance ({})\n",
                if primal_fragile {
                    blanks.as_str()
                } else {
                    problem
                },
                g1(dual_ratio),
                g1(max_objective),
                g1(use_dual_feasibility_tolerance)
            );
        }
        crate::highs_log_user!(
            log_options,
            HighsLogType::kWarning,
            "   Improve scaling and/or increase tolerance\n"
        );
        warning_issued = true;
    }
    // Give a scary warning to users asking too much of HiGHS!
    if warning_issued {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kWarning,
            "{} is badly scaled, {}which may compromise the speed, accuracy and reliability of solvers in HiGHS\n",
            problem,
            if numerically_fragile {
                "particularly relative to tolerances, "
            } else {
                ""
            }
        );
    }
}

#[cfg(test)]
mod h16_scaling_report_tests {
    use super::*;
    use crate::lp_data::highs_callback::{
        HighsCallbackFunctionType, HighsCallbackInput, HighsCallbackOutput,
    };
    use std::sync::{Arc, Mutex};

    /// min z s.t. z - 3 x1 - 200 x2 = 0, x1 + x2 >= 1, z free, x1, x2 in [0, x_upper].
    /// z is a free column singleton in an equation, so the effective costs of x1 and x2
    /// are 3 and 200 (HighsSolution.cpp::getEffectiveCosts, 1.16.0).
    fn objective_row_lp(x_upper: f64, integer: bool) -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 3,
            num_row_: 2,
            col_cost_: vec![1.0, 0.0, 0.0],
            col_lower_: vec![-kHighsInf, 0.0, 0.0],
            col_upper_: vec![kHighsInf, 10.0, x_upper],
            row_lower_: vec![0.0, 1.0],
            row_upper_: vec![0.0, kHighsInf],
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 1, 3, 5];
        lp.a_matrix_.index_ = vec![0, 0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![1.0, -3.0, 1.0, -200.0, 1.0];
        if integer {
            lp.integrality_ = vec![
                HighsVarType::kContinuous,
                HighsVarType::kInteger,
                HighsVarType::kInteger,
            ];
        }
        lp
    }

    /// The user log of assess_excessive_objective_bound_scaling, captured via the callback.
    fn report(lp: &HighsLp, h16_scaling_report: bool) -> String {
        let out = Arc::new(Mutex::new(String::new()));
        let sink = out.clone();
        let cb = HighsCallbackFunctionType::new(
            move |_t: HighsInt,
                  msg: &str,
                  _o: &HighsCallbackOutput,
                  _i: &mut HighsCallbackInput| {
                sink.lock().unwrap().push_str(msg);
            },
        );
        let mut options = HighsOptions {
            h16_scaling_report,
            ..Default::default()
        };
        options.log_options = HighsLogOptions {
            output_flag: true,
            log_to_console: false,
            user_callback: Some(cb),
            user_callback_active: true,
            ..Default::default()
        };
        let mut data = HighsUserScaleData::default();
        assess_excessive_objective_bound_scaling(&options, lp, &mut data);
        out.lock().unwrap().clone()
    }

    #[test]
    fn h16_scaling_report_effective_costs() {
        let lp = objective_row_lp(10.0, false);
        let effective =
            crate::lp_data::highs_solution::get_effective_costs(&lp, &HighsOptions::default());
        assert_eq!(effective, vec![0.0, 3.0, 200.0]);
    }

    #[test]
    fn h16_scaling_report_lp_on() {
        assert_eq!(
            report(&objective_row_lp(10.0, false), true),
            "Coefficient ranges:\n\
             \x20 Matrix  [1e+00, 2e+02]\n\
             \x20 Cost    [1e+00, 1e+00]\n\
             \x20 Cost    [3e+00, 2e+02] (effective)\n\
             \x20 Bound   [1e+01, 1e+01]\n\
             \x20 RHS     [1e+00, 1e+00]\n"
        );
    }

    #[test]
    fn h16_scaling_report_lp_off() {
        assert_eq!(
            report(&objective_row_lp(10.0, false), false),
            "Coefficient ranges:\n\
             \x20 Matrix  [1e+00, 2e+02]\n\
             \x20 Cost    [1e+00, 1e+00]\n\
             \x20 Bound   [1e+01, 1e+01]\n\
             \x20 RHS     [1e+00, 1e+00]\n"
        );
    }

    #[test]
    fn h16_scaling_report_mip_on() {
        assert_eq!(
            report(&objective_row_lp(2e6, true), true),
            "Coefficient ranges:\n\
             \x20 Matrix  [1e+00, 2e+02]\n\
             \x20 Cost    [1e+00, 1e+00]\n\
             \x20 Cost    [3e+00, 2e+02] (effective)\n\
             \x20 Bound   [0e+00, 0e+00] (continuous)\n\
             \x20 Bound   [1e+01, 2e+06] (non-continuous)\n\
             \x20 RHS     [1e+00, 1e+00]\n\
             WARNING: Problem has some excessively large bounds on non-continuous variables\n"
        );
    }

    #[test]
    fn h16_scaling_report_mip_off() {
        assert_eq!(
            report(&objective_row_lp(2e6, true), false),
            "Coefficient ranges:\n\
             \x20 Matrix  [1e+00, 2e+02]\n\
             \x20 Cost    [1e+00, 1e+00]\n\
             \x20 Bound   [1e+01, 2e+06]\n\
             \x20 RHS     [1e+00, 1e+00]\n\
             WARNING: Problem has some excessively large column bounds\n"
        );
    }

    /// A large constraint bound: the 1.16.0 wording, the bound suggestion, and the
    /// tolerance-ratio and badly-scaled warnings (neos-4413714-turia prints the same
    /// lines with its own values).
    #[test]
    fn h16_scaling_report_fragile_on() {
        let mut lp = objective_row_lp(10.0, false);
        lp.row_lower_[1] = 1e9;
        assert_eq!(
            report(&lp, true),
            "Coefficient ranges:\n\
             \x20 Matrix  [1e+00, 2e+02]\n\
             \x20 Cost    [1e+00, 1e+00]\n\
             \x20 Cost    [3e+00, 2e+02] (effective)\n\
             \x20 Bound   [1e+01, 1e+01]\n\
             \x20 RHS     [1e+09, 1e+09]\n\
             WARNING: Problem has some excessively large bounds on constraints\n\
             WARNING:    Consider scaling the    bounds by 1e-3, or setting the user_bound_scale option to -10\n\
             WARNING: Problem has excessively large ratio (1e+16) between max variable or constraint bound (1e+09) and primal feasibility tolerance (1e-07)\n\
             WARNING:    Improve scaling and/or increase tolerance\n\
             WARNING: Problem is badly scaled, particularly relative to tolerances, which may compromise the speed, accuracy and reliability of solvers in HiGHS\n"
        );
    }
}
