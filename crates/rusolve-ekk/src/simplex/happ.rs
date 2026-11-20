//! HiGHS source: simplex/HApp.h (returnFromSolveLpSimplex, solveLpSimplex)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The driver around `HEkk::solve`: scale, move the LP into EKK, solve the
//! scaled LP, unscale, and (if the unscaled solution has infeasibilities)
//! refine by solving the unscaled LP with scaled NLA.
//!
//! Ownership (crate rules (a), (d), (e)): `HighsLpSolverObject` has no `HEkk`,
//! so the `HEkk` is passed as `ekk`. `HEkk::options_` is an owned clone made
//! by `moveLp`, whereas in C++ it points at `solver_object.options_`; every
//! write that `solveLpSimplex` makes to `options` after a `moveLp` is
//! therefore mirrored into `ekk.options_`, so that the solver sees it as the
//! C++ one does.
//!
//! A useful but unvalidated HiGHS basis goes through
//! `formSimplexLpBasisAndFactor` + `refineBasis` (lp_data/highs_basis_form.rs). Profiling clocks are dropped.

#![allow(
    clippy::int_plus_one,
    clippy::needless_late_init,
    clippy::collapsible_match,
    clippy::assertions_on_constants
)]

use crate::lp_data::highs_solve::highs_status_from_highs_model_status;
use crate::lp_data::lp_scale::{consider_scaling, unscale_solution};
use crate::lp_data::presolve_support::refine_basis;
use crate::prelude::*;

/// `HEkk::getHighsBasis(ekk_lp)` where `ekk_lp` is `ekk.lp_` itself: the
/// borrow checker cannot pass `&mut ekk.lp_` to a `&self` method, so the
/// fields `getHighsBasis` reads from `use_lp` (dimensions and bounds) are
/// moved into a shallow `HighsLp` for the call and moved back afterwards.
fn get_highs_basis_of_own_lp(ekk: &mut HEkk) -> HighsBasis {
    let mut use_lp = HighsLp {
        num_col_: ekk.lp_.num_col_,
        num_row_: ekk.lp_.num_row_,
        col_lower_: std::mem::take(&mut ekk.lp_.col_lower_),
        col_upper_: std::mem::take(&mut ekk.lp_.col_upper_),
        row_lower_: std::mem::take(&mut ekk.lp_.row_lower_),
        row_upper_: std::mem::take(&mut ekk.lp_.row_upper_),
        ..HighsLp::default()
    };
    let basis = ekk.get_highs_basis(&mut use_lp);
    ekk.lp_.col_lower_ = use_lp.col_lower_;
    ekk.lp_.col_upper_ = use_lp.col_upper_;
    ekk.lp_.row_lower_ = use_lp.row_lower_;
    ekk.lp_.row_upper_ = use_lp.row_upper_;
    basis
}

/// HiGHS: HApp.h::returnFromSolveLpSimplex (lines 37-92, 56 lines)
pub fn return_from_solve_lp_simplex(
    ekk: &mut HEkk,
    solver_object: &mut HighsLpSolverObject,
    return_status: HighsStatus,
) -> HighsStatus {
    let mut return_status = return_status;
    // Copy the simplex iteration count to highs_info_ from ekk_instance
    solver_object.highs_info_.simplex_iteration_count = ekk.iteration_count_;
    // Stop whichever clock was running: profiling out of scope
    //
    // Ensure that the incumbent LP is neither moved, nor scaled
    debug_assert!(!solver_object.lp_.is_moved_);
    debug_assert!(!solver_object.lp_.is_scaled_);
    // Cannot expect any more with an error return, and safer to clear
    // HEkk than try to retain any data
    if return_status == HighsStatus::kError {
        ekk.clear();
        return return_status;
    }
    //
    // Ensure that there is an invert for the current LP
    debug_assert!(ekk.status_.has_invert);
    // Ensure that simplex NLA is set up
    debug_assert!(ekk.status_.has_nla);
    // Set the simplex NLA scaling
    ekk.set_nla_pointers_for_lp_and_scale(solver_object.lp_);
    debug_assert!(ekk.debug_nla_scaling_ok(solver_object.lp_));
    let alt_debug_level: HighsInt = -1;
    if ekk.debug_nla_check_invert("HApp: returnFromSolveLpSimplex", alt_debug_level)
        == HighsDebugStatus::kError
    {
        crate::highs_log_user!(
            &solver_object.options_.log_options,
            HighsLogType::kError,
            "Error in basis matrix inverse after solving the LP\n"
        );
        return_status = HighsStatus::kError;
    }
    if solver_object.model_status_ == HighsModelStatus::kOptimal {
        solver_object.highs_info_.num_complementarity_violations = 0;
        solver_object.highs_info_.max_complementarity_violation = 0.0;
    }
    return_status
}

/// HiGHS: HApp.h::solveLpSimplex (lines 94-532, 439 lines)
#[allow(unused_assignments)]
pub fn solve_lp_simplex(ekk: &mut HEkk, solver_object: &mut HighsLpSolverObject) -> HighsStatus {
    let mut return_status = HighsStatus::kOk;
    let mut call_status: HighsStatus;
    let mut scaled_model_status = HighsModelStatus::kUnknown;

    // Check that any retained Ekk data - basis and NLA - are OK on entry
    let retained_ekk_data_ok =
        ekk.debug_retained_data_ok(solver_object.lp_) != HighsDebugStatus::kLogicalError;
    if !retained_ekk_data_ok {
        crate::highs_log_user!(
            &solver_object.options_.log_options,
            HighsLogType::kError,
            "solveLpSimplex: Retained Ekk data not OK on entry\n"
        );
        debug_assert!(retained_ekk_data_ok);
        return_status = HighsStatus::kError;
    }
    // Profiling: out of scope
    //
    // Copy the simplex iteration count from highs_info_ to ekk_instance, just for
    // convenience
    ekk.iteration_count_ = solver_object.highs_info_.simplex_iteration_count;

    // Reset the model status and HighsInfo values in case of premature
    // return
    crate::lp_data::highs_solution::reset_model_status_and_highs_info_so(solver_object);

    // Initialise the simplex stats
    ekk.initialise_simplex_stats();

    // Assumes that the LP has a positive number of rows, since
    // unconstrained LPs should be solved in solveLp
    let positive_num_row = solver_object.lp_.num_row_ > 0;
    debug_assert!(positive_num_row);
    if !positive_num_row {
        crate::highs_log_user!(
            &solver_object.options_.log_options,
            HighsLogType::kError,
            "solveLpSimplex called for LP with non-positive ({}) number of constraints\n",
            solver_object.lp_.num_row_
        );
        return return_from_solve_lp_simplex(ekk, solver_object, HighsStatus::kError);
    }
    // On entry to solveLpSimplex, the incumbent LP is assumed to be
    // unscaled and not moved
    debug_assert!(!solver_object.lp_.is_scaled_);
    debug_assert!(!solver_object.lp_.is_moved_);
    // Consider scaling the LP - either with any existing scaling, or by
    // considering computing scaling factors if there are none - and
    // then move to EKK
    consider_scaling(solver_object.options_, solver_object.lp_);
    //
    let was_scaled = solver_object.lp_.is_scaled_;
    if !ekk.status_.has_basis && !solver_object.basis_.valid && solver_object.basis_.useful {
        // There is no simplex basis, but there is a useful HiGHS basis
        // that is not validated (HApp.h:179-190; reached after deleteRows
        // removed a nonbasic row)
        debug_assert!(solver_object.basis_.col_status.len() == solver_object.lp_.num_col_ as usize);
        debug_assert!(solver_object.basis_.row_status.len() == solver_object.lp_.num_row_ as usize);
        let return_status = form_simplex_lp_basis_and_factor(ekk, solver_object, false);
        if return_status != HighsStatus::kOk {
            // Known unfaithful spot: C++
            // returns with the incumbent LP still scaled, which trips the
            // `!is_scaled_` assert of returnFromSolveLpSimplex; unscale first.
            if solver_object.lp_.is_scaled_ {
                solver_object.lp_.unapply_scale();
            }
            return return_from_solve_lp_simplex(ekk, solver_object, HighsStatus::kError);
        }
        // formSimplexLpBasisAndFactor may introduce variables with
        // HighsBasisStatus::kNonbasic, so refine it
        refine_basis(
            solver_object.lp_,
            solver_object.solution_,
            solver_object.basis_,
        );
        solver_object.basis_.valid = true;
    }
    // Check that scaling has not been removed
    debug_assert!(solver_object.lp_.is_scaled_ == was_scaled);
    // Move the LP to EKK, updating other EKK pointers and any simplex
    // NLA pointers, since they may have moved if the LP has been
    // modified
    ekk.move_lp(solver_object);
    if !ekk.status_.has_basis {
        // There is no simplex basis, so use any HiGHS basis
        if solver_object.basis_.valid {
            call_status = ekk.set_basis_from(solver_object.basis_);
            if call_status == HighsStatus::kError {
                solver_object
                    .lp_
                    .move_back_lp_and_unapply_scaling(&mut ekk.lp_);
                return return_from_solve_lp_simplex(ekk, solver_object, call_status);
            }
        } else {
            // Starting from a logical basis, so consider dualizing and/or
            // permuting the LP
            let options = &*solver_object.options_;
            // lp_dualize_auto: rusolve lever, also
            // dualize when the LP has at least 3x more rows than columns
            let auto_dualize =
                options.lp_dualize_auto && ekk.lp_.num_row_ >= 3 * ekk.lp_.num_col_.max(1);
            if options.simplex_dualize_strategy == kHighsOptionChoose
                || options.simplex_dualize_strategy == kHighsOptionOn
                || auto_dualize
            {
                // Dualize unless we choose not to
                let mut dualize_lp = true;
                if options.simplex_dualize_strategy == kHighsOptionChoose && !auto_dualize {
                    // C++ reads the moved-from incumbent_lp, whose scalars
                    // num_row_/num_col_ are left as they were: read ekk.lp_
                    if ekk.lp_.num_row_ < 10 * ekk.lp_.num_col_ {
                        dualize_lp = false;
                    }
                }
                if dualize_lp {
                    ekk.dualize();
                }
            }
            if options.simplex_permute_strategy == kHighsOptionChoose
                || options.simplex_permute_strategy == kHighsOptionOn
            {
                // Permute the LP
                ekk.permute();
            }
        }
    }
    // These local illegal values are over-written with correct values
    // if the scaled LP is solved in order to take correct algorithmic
    // decisions if the unscaled LP is solved later
    let mut num_unscaled_primal_infeasibilities: HighsInt = kHighsIllegalInfeasibilityCount;
    let mut num_unscaled_dual_infeasibilities: HighsInt = kHighsIllegalInfeasibilityCount;
    // Record whether the unscaled LP is to be solved, and has been
    // solved. It's not solved if it's proved to be infeasible using the
    // ray from the scaled LP, in which case the solution (FWIW) and
    // basis are taken from the unscaled solution of the scaled LP.
    let mut solve_unscaled_lp = false;
    let mut solved_unscaled_lp = false;
    // C++ reads incumbent_lp.scale_.has_scaling after moveLp. HighsLp's
    // implicit move leaves the moved-from scalars as they were, so it
    // still reads the original value; `mem::take` reset it here, so the
    // same value is read from ekk.lp_.
    if !ekk.lp_.scale_.has_scaling {
        //
        // Solve the unscaled LP with unscaled NLA
        //
        solve_unscaled_lp = true;
        return_status = ekk.solve(false);
        solved_unscaled_lp = true;
        ekk.unpermute();
        ekk.undualize();
        debug_assert!(!ekk.status_.is_permuted && !ekk.status_.is_dualized);
        let options = &*solver_object.options_;
        if options.cost_scale_factor != 0 {
            let cost_scale_factor = 2.0_f64.powf(-(options.cost_scale_factor as f64));
            crate::highs_log_dev!(
                &options.log_options,
                HighsLogType::kInfo,
                "Objective = {:11.4}\n",
                cost_scale_factor * ekk.info_.dual_objective_value
            );
            ekk.model_status_ = HighsModelStatus::kNotset;
            return_status = HighsStatus::kError;
        }
        //
    } else {
        // Indicate that there is no (current) need to refine the solution
        // by solving the unscaled LP with scaled NLA
        let mut refine_solution = false;
        let simplex_unscaled_solution_strategy =
            solver_object.options_.simplex_unscaled_solution_strategy;
        if simplex_unscaled_solution_strategy == kSimplexUnscaledSolutionStrategyNone
            || simplex_unscaled_solution_strategy == kSimplexUnscaledSolutionStrategyRefine
        {
            // Check that the incumbent LP has scaling and is scaled
            debug_assert!(ekk.lp_.scale_.has_scaling);
            debug_assert!(ekk.lp_.is_scaled_);
            //
            // Solve the scaled LP!
            //
            return_status = ekk.solve(false);
            ekk.unpermute();
            ekk.undualize();
            debug_assert!(!ekk.status_.is_permuted && !ekk.status_.is_dualized);
            //
            if solver_object.options_.cost_scale_factor != 0 {
                let cost_scale_factor =
                    2.0_f64.powf(-(solver_object.options_.cost_scale_factor as f64));
                crate::highs_log_dev!(
                    &solver_object.options_.log_options,
                    HighsLogType::kInfo,
                    "Objective = {:11.4}\n",
                    cost_scale_factor * ekk.info_.dual_objective_value
                );
                ekk.model_status_ = HighsModelStatus::kNotset;
                return_status = HighsStatus::kError;
            }
            if return_status == HighsStatus::kError {
                solver_object
                    .lp_
                    .move_back_lp_and_unapply_scaling(&mut ekk.lp_);
                return return_from_solve_lp_simplex(ekk, solver_object, return_status);
            }
            // Copy solution data from the EKK instance
            scaled_model_status = ekk.model_status_;
            solver_object.highs_info_.objective_function_value = ekk.info_.primal_objective_value;
            solver_object.highs_info_.simplex_iteration_count = ekk.iteration_count_;
            *solver_object.solution_ = ekk.get_solution();
            *solver_object.basis_ = get_highs_basis_of_own_lp(ekk);
            debug_assert!(solver_object.basis_.valid);
            solver_object.highs_info_.basis_validity = kBasisValidityValid;
            solver_object
                .lp_
                .move_back_lp_and_unapply_scaling(&mut ekk.lp_);
            // Now that the incumbent LP is unscaled, to use the simplex NLA
            // requires scaling to be applied
            ekk.set_nla_pointers_for_lp_and_scale(solver_object.lp_);
            unscale_solution(solver_object.solution_, &solver_object.lp_.scale_);
            // Determine whether the unscaled LP has been solved
            get_unscaled_infeasibilities(
                solver_object.options_,
                &solver_object.lp_.scale_,
                &ekk.basis_,
                &ekk.info_,
                solver_object.highs_info_,
            );
            let highs_info = &*solver_object.highs_info_;
            num_unscaled_primal_infeasibilities = highs_info.num_primal_infeasibilities;
            num_unscaled_dual_infeasibilities = highs_info.num_dual_infeasibilities;
            // Determine whether the unscaled solution has infeasibilities
            // after the scaled LP has been solved to optimality
            let scaled_optimality_but_unscaled_infeasibilities = scaled_model_status
                == HighsModelStatus::kOptimal
                && (num_unscaled_primal_infeasibilities != 0
                    || num_unscaled_dual_infeasibilities != 0);
            // Determine whether the unscaled solution has primal
            // infeasibilities after the scaled LP has been solved to the
            // objective target
            let scaled_objective_target_but_unscaled_primal_infeasibilities = scaled_model_status
                == HighsModelStatus::kObjectiveTarget
                && highs_info.num_primal_infeasibilities > 0;
            // Determine whether the unscaled solution has dual
            // infeasibilities after the scaled LP has been solved to the
            // objective bound
            let scaled_objective_bound_but_unscaled_dual_infeasibilities = scaled_model_status
                == HighsModelStatus::kObjectiveBound
                && highs_info.num_dual_infeasibilities > 0;
            if scaled_optimality_but_unscaled_infeasibilities
                || scaled_objective_target_but_unscaled_primal_infeasibilities
                || scaled_objective_bound_but_unscaled_dual_infeasibilities
            {
                crate::highs_log_dev!(
                    &solver_object.options_.log_options,
                    HighsLogType::kInfo,
                    "After unscaling with status {}, have num/max/sum primal ({}/{}/{}) and dual ({}/{}/{}) unscaled infeasibilities\n",
                    util_model_status_to_string(scaled_model_status),
                    highs_info.num_primal_infeasibilities,
                    highs_info.max_primal_infeasibility,
                    highs_info.sum_primal_infeasibilities,
                    highs_info.num_dual_infeasibilities,
                    highs_info.max_dual_infeasibility,
                    highs_info.sum_dual_infeasibilities
                );
            }
            // Determine whether refinement will take place
            refine_solution = simplex_unscaled_solution_strategy
                == kSimplexUnscaledSolutionStrategyRefine
                && (scaled_optimality_but_unscaled_infeasibilities
                    || scaled_model_status == HighsModelStatus::kInfeasible
                    || scaled_model_status == HighsModelStatus::kUnboundedOrInfeasible
                    || scaled_model_status == HighsModelStatus::kUnbounded
                    || scaled_objective_bound_but_unscaled_dual_infeasibilities
                    || scaled_objective_target_but_unscaled_primal_infeasibilities
                    || scaled_model_status == HighsModelStatus::kUnknown);
            // Handle the case when refinement will not take place
            if !refine_solution {
                solver_object.model_status_ = scaled_model_status;
                return_status = highs_status_from_highs_model_status(solver_object.model_status_);
                return return_from_solve_lp_simplex(ekk, solver_object, return_status);
            }
        } else {
            // LP is scaled, but simplex_unscaled_solution_strategy is
            // kSimplexUnscaledSolutionStrategyDirect, so have to move back
            // the LP and unscale it
            debug_assert!(
                simplex_unscaled_solution_strategy == kSimplexUnscaledSolutionStrategyDirect
            );
            solver_object
                .lp_
                .move_back_lp_and_unapply_scaling(&mut ekk.lp_);
        }
        debug_assert!(
            simplex_unscaled_solution_strategy == kSimplexUnscaledSolutionStrategyDirect
                || refine_solution
        );
        // Solve the unscaled LP using scaled NLA. This requires pointers of
        // a scaled matrix to be passed to the HFactor instance. Use the
        // incumbent LP for this.
        //
        // Check that the incumbent LP has been moved back and is unscaled
        debug_assert!(!solver_object.lp_.is_moved_);
        debug_assert!(!solver_object.lp_.is_scaled_);
        // Move the incumbent LP
        ekk.move_lp(solver_object);
        // If refining after proving primal infeasibility of the scaled
        // LP, see whether the proof still holds for the unscaled LP. If
        // it does, then there's no need to solve the unscaled LP
        solve_unscaled_lp = true;
        if scaled_model_status == HighsModelStatus::kInfeasible
            && ekk.exit_algorithm_ == SimplexAlgorithm::kDual
        {
            debug_assert!(ekk.dual_ray_record_.index != kNoRayIndex);
        }
        if scaled_model_status == HighsModelStatus::kInfeasible
            && ekk.dual_ray_record_.index != kNoRayIndex
        {
            // setNlaPointersForLpAndScale: no-op (crate rule (b))
            if ekk.proof_of_primal_infeasibility() {
                solve_unscaled_lp = false;
            }
        }
        if solve_unscaled_lp {
            // Save options/strategies that may be changed
            let simplex_strategy = solver_object.options_.simplex_strategy;
            let dual_simplex_cost_perturbation_multiplier = solver_object
                .options_
                .dual_simplex_cost_perturbation_multiplier;
            let simplex_dual_edge_weight_strategy = ekk.info_.dual_edge_weight_strategy;
            // #1865: see HApp.h
            if num_unscaled_primal_infeasibilities == 0
                || scaled_model_status == HighsModelStatus::kObjectiveTarget
            {
                // Only dual infeasibilities, or objective target reached (in
                // primal phase 2) - so primal infeasibilities should be small
                // - so use primal simplex phase 2
                solver_object.options_.simplex_strategy = kSimplexStrategyPrimal;
                ekk.options_.simplex_strategy = kSimplexStrategyPrimal;
                if scaled_model_status == HighsModelStatus::kObjectiveTarget {
                    crate::highs_log_dev!(
                        &solver_object.options_.log_options,
                        HighsLogType::kInfo,
                        "solveLpSimplex: Calling primal simplex after scaled_model_status == HighsModelStatus::kObjectiveTarget: solve = {}; tick = {}; iter = {}\n",
                        ekk.debug_solve_call_num_,
                        ekk.debug_initial_build_synthetic_tick_,
                        ekk.iteration_count_
                    );
                }
            } else {
                // There are unscaled primal infeasibilities, so use dual
                // simplex
                debug_assert!(num_unscaled_primal_infeasibilities > 0);
                if solver_object.options_.simplex_strategy != kSimplexStrategyDual {
                    crate::highs_log_dev!(
                        &solver_object.options_.log_options,
                        HighsLogType::kInfo,
                        "Forcing change from {} to {}\n",
                        ekk.simplex_strategy_to_string(solver_object.options_.simplex_strategy),
                        ekk.simplex_strategy_to_string(kSimplexStrategyDual)
                    );
                }

                solver_object.options_.simplex_strategy = kSimplexStrategyDual;
                ekk.options_.simplex_strategy = kSimplexStrategyDual;
                // Using dual simplex, so force Devex if starting from an advanced
                // basis with no steepest edge weights
                if (ekk.status_.has_basis || solver_object.basis_.valid)
                    && !ekk.status_.has_dual_steepest_edge_weights
                {
                    ekk.info_.dual_edge_weight_strategy = kSimplexEdgeWeightStrategyDevex;
                }
            }
            //
            // Solve the unscaled LP with scaled NLA
            //
            // Force the simplex solver to start in phase 1 if solving the
            // LP directly as unscaled, or using primal simplex to clean up
            // small dual infeasibilities after the scaled LP yielded model
            // status HighsModelStatus::kObjectiveTarget. Otherwise force
            // the simplex solver to start in phase 2
            //
            let force_phase1 = (simplex_unscaled_solution_strategy
                == kSimplexUnscaledSolutionStrategyDirect)
                || (scaled_model_status == HighsModelStatus::kObjectiveTarget);
            let force_phase2 = (simplex_unscaled_solution_strategy
                != kSimplexUnscaledSolutionStrategyDirect)
                && (scaled_model_status != HighsModelStatus::kObjectiveTarget);
            debug_assert!(force_phase2 == !force_phase1);
            return_status = ekk.solve(force_phase2);
            solved_unscaled_lp = true;
            if scaled_model_status != HighsModelStatus::kObjectiveBound
                && ekk.model_status_ == HighsModelStatus::kObjectiveBound
            {
                // it may happen that the unscaled LP detected status kObjectiveBound
                // for the first time in which case we again call solve with primal
                // simplex if not dual feasible
                let objective_bound_refinement = ekk.info_.num_dual_infeasibilities > 0;
                if objective_bound_refinement {
                    solver_object.options_.simplex_strategy = kSimplexStrategyPrimal;
                    ekk.options_.simplex_strategy = kSimplexStrategyPrimal;
                    return_status = ekk.solve(force_phase2);
                }
            }
            // Restore the options/strategies that may have been changed
            solver_object.options_.simplex_strategy = simplex_strategy;
            solver_object
                .options_
                .dual_simplex_cost_perturbation_multiplier =
                dual_simplex_cost_perturbation_multiplier;
            ekk.options_.simplex_strategy = simplex_strategy;
            ekk.options_.dual_simplex_cost_perturbation_multiplier =
                dual_simplex_cost_perturbation_multiplier;
            ekk.info_.dual_edge_weight_strategy = simplex_dual_edge_weight_strategy;
        }
    }
    if solved_unscaled_lp {
        // Copy solution data from the EKK istance
        scaled_model_status = ekk.model_status_;
        solver_object.highs_info_.objective_function_value = ekk.info_.primal_objective_value;
        solver_object.highs_info_.simplex_iteration_count = ekk.iteration_count_;
        *solver_object.solution_ = ekk.get_solution();
        *solver_object.basis_ = get_highs_basis_of_own_lp(ekk);
        debug_assert!(solver_object.basis_.valid);
        solver_object.highs_info_.basis_validity = kBasisValidityValid;
    }
    // Move the incumbent LP back from Ekk
    *solver_object.lp_ = std::mem::take(&mut ekk.lp_);
    solver_object.lp_.is_moved_ = false;
    ekk.set_nla_pointers_for_lp_and_scale(solver_object.lp_);
    if return_status == HighsStatus::kError {
        // Error return, so make sure that the unscaled and scaled model
        // status values are the same, and that they correspond to an
        // error return
        solver_object.model_status_ = scaled_model_status;
        return_status = highs_status_from_highs_model_status(solver_object.model_status_);
        debug_assert!(return_status == HighsStatus::kError);
        return return_from_solve_lp_simplex(ekk, solver_object, HighsStatus::kError);
    }
    let highs_info = &mut *solver_object.highs_info_;
    if solved_unscaled_lp {
        // The unscaled LP has been solved - either directly, or because
        // there was no scaling. Copy values into HighsInfo from solving
        // the scaled LP.
        debug_assert!(solve_unscaled_lp);
        highs_info.num_primal_infeasibilities = ekk.info_.num_primal_infeasibilities;
        highs_info.max_primal_infeasibility = ekk.info_.max_primal_infeasibility;
        highs_info.sum_primal_infeasibilities = ekk.info_.sum_primal_infeasibilities;
        highs_info.num_dual_infeasibilities = ekk.info_.num_dual_infeasibilities;
        highs_info.max_dual_infeasibility = ekk.info_.max_dual_infeasibility;
        highs_info.sum_dual_infeasibilities = ekk.info_.sum_dual_infeasibilities;
    } else {
        // The unscaled LP has not been solved because the scaled LP was
        // infeasible and the proof of infeasibility held for the unscaled
        // LP. Hence the values in HighsInfo that are set (above) by the
        // call to getUnscaledInfeasibilities are correct
        debug_assert!(!solve_unscaled_lp);
        debug_assert!(scaled_model_status == HighsModelStatus::kInfeasible);
    }
    set_solution_status(highs_info);
    solver_object.model_status_ = scaled_model_status;
    return_status = highs_status_from_highs_model_status(solver_object.model_status_);
    return_from_solve_lp_simplex(ekk, solver_object, return_status)
}
