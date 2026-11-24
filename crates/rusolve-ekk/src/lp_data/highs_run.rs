//! HiGHS source: lp_data/Highs.cpp (passModel, calledOptimizeModel, callSolveLp,
//! returnFromOptimizeModel, invalidate*), lp_data/HighsInterface.cpp
//! (infeasibleBoundsOk, zeroIterationCounts, callLpKktCheck)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The public entry point of the crate: a `Highs` object holding one LP,
//! with `pass_model` + `run` reproducing `Highs::passModel(HighsLp)` +
//! `Highs::run()` for an LP with `solver = "simplex"`/"choose", `threads = 1`.
//! `run` keeps the C++ gate of `calledOptimizeModel`: with `presolve = "off"`
//! (or an unconstrained LP) it takes the "Solving LP without presolve or
//! useful basis" path, unchanged and bit-exact; otherwise (default
//! `presolve = "choose"`) it calls `called_optimize_model_presolve`
//! (highs_presolve.rs, the LP presolve path).
//!
//! A MIP is routed to `call_solve_mip` (the MIP port). The semi-variable,
//! scaling and modification call sites are wired where C++ calls them:
//! `user_scale` / `user_unscale` and `optimize_highs` in `run`;
//! `clear_solver` (when `!use_warm_start`), `handle_inf_cost`,
//! `complete_solution_from_discrete_assignment`, `assess_semi_variables`,
//! `relax_semi_variables` (`solve_relaxation`) and `basis_for_solution` in
//! `called_optimize_model`; `restore_inf_cost` / `unapply_mods` (when
//! `undo_mods`) in `return_from_optimize_model`; `without_semi_variables` and
//! `active_modified_upper_bounds` in `call_solve_mip`.
//!
//! Out of scope: names, ranging, IIS, the run_data time fields, logging header, timers.

#![allow(
    clippy::int_plus_one,
    clippy::needless_late_init,
    clippy::collapsible_match,
    clippy::assertions_on_constants
)]

use crate::lp_data::highs_solve::assess_excessive_objective_bound_scaling;
use crate::lp_data::highs_solve::{highs_status_from_highs_model_status, solve_lp};
use crate::lp_data::hstruct::{HighsLinearObjective, HighsUserScaleData};
use crate::lp_data::lp_assess::{assess_lp, lp_dimensions_ok};
use crate::lp_data::lp_utils_semi::{assess_semi_variables, relax_semi_variables};
use crate::lp_data::lp_utils_semi_mip::{active_modified_upper_bounds, without_semi_variables};
use crate::lp_data::presolve_support::{get_kkt_failures, refine_basis};
use crate::prelude::*;
use crate::util::sparse_matrix::HighsSparseMatrix;

thread_local! {
    /// rusolve: nesting depth of `Highs::run` on this thread (debug check in
    /// `reset_scheduler_for_threads`).
    static RUN_DEPTH: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
}

struct RunDepth;

impl RunDepth {
    fn enter() -> Self {
        RUN_DEPTH.with(|d| d.set(d.get() + 1));
        RunDepth
    }
}

impl Drop for RunDepth {
    fn drop(&mut self) {
        RUN_DEPTH.with(|d| d.set(d.get() - 1));
    }
}

/// rusolve: whether a `Highs::run` is active on this thread.
pub(crate) fn in_run() -> bool {
    RUN_DEPTH.with(|d| d.get() > 0)
}

/// HiGHS: Highs.h::Highs (the members the LP simplex path uses)
pub struct Highs {
    pub options_: HighsOptions,
    pub timer_: HighsTimer,
    pub callback_: HighsCallback,
    /// HiGHS: `model_.lp_`
    pub model_lp_: HighsLp,
    pub solution_: HighsSolution,
    pub basis_: HighsBasis,
    pub info_: HighsInfo,
    /// HiGHS: `run_data_` (Highs.h), the HighsRunDataStruct counts (highs_run_data.rs)
    pub run_data_: HighsRunData,
    pub model_status_: HighsModelStatus,
    pub ekk_instance_: HEkk,
    /// HiGHS: `presolve_` (Highs.h)
    pub presolve_: PresolveComponent,
    /// HiGHS: `model_presolve_status_` (Highs.h initialiser:
    /// `HighsPresolveStatus::kNotPresolved`)
    pub model_presolve_status_: HighsPresolveStatus,
    /// HiGHS: `presolve_log_` (Highs.h)
    pub presolve_log_: HighsPresolveLog,
    /// HiGHS: `model_.hessian_` (Highs.h `HighsModel model_`; kept apart
    /// from `model_lp_`)
    pub model_hessian_: HighsHessian,
    /// HiGHS: `icrash_info_` (Highs.h)
    pub icrash_info_: ICrashInfo,
    /// HiGHS: `ranging_` (Highs.h)
    pub ranging_: HighsRanging,
    /// HiGHS: `iis_` (Highs.h)
    pub iis_: HighsIis,
    /// HiGHS: `multi_linear_objective_` (Highs.h:1604)
    pub multi_linear_objective_: Vec<HighsLinearObjective>,
    /// HiGHS: `presolved_model_.lp_` (Highs.h:1606 `HighsModel presolved_model_`), set
    /// only by the public `Highs::presolve`. The presolved model has no
    /// Hessian here: `Highs::presolve` rejects a QP.
    pub presolved_model_lp_: HighsLp,
    /// HiGHS: `saved_objective_and_solution_` (Highs.h:1615)
    pub saved_objective_and_solution_: Vec<HighsObjectiveSolution>,
    /// C1/C5 (no HiGHS counterpart): per-phase statistics of the last MIP solve, reporting only
    /// (`mip/mip_stats.rs`).
    pub mip_stats_: Option<crate::mip::mip_stats::MipStats>,
    /// HiGHS: `standard_form_valid_` (Highs.h:1621)
    pub standard_form_valid_: bool,
    /// HiGHS: `standard_form_offset_` (Highs.h:1622)
    pub standard_form_offset_: f64,
    /// HiGHS: `standard_form_cost_` (Highs.h:1623)
    pub standard_form_cost_: Vec<f64>,
    /// HiGHS: `standard_form_rhs_` (Highs.h:1624)
    pub standard_form_rhs_: Vec<f64>,
    /// HiGHS: `standard_form_matrix_` (Highs.h:1625)
    pub standard_form_matrix_: HighsSparseMatrix,
    /// HiGHS: `max_threads_` (Highs.h:1633, initialiser `0`); set by
    /// `initializeMultiThreading`.
    pub max_threads_: HighsInt,
    /// C4 (ours): when set, `call_solve_mip` copies the MIP solver counters into `mip_probe_` after the solve (read-only; default off).
    pub capture_mip_probe_: bool,
    /// C4 (ours): see `capture_mip_probe_`.
    pub mip_probe_: Option<crate::mip::mip_probe::MipProbeSnapshot>,
}

impl Default for Highs {
    fn default() -> Self {
        let mut highs = Self {
            options_: HighsOptions::default(),
            timer_: HighsTimer::new(),
            callback_: HighsCallback::default(),
            model_lp_: HighsLp::default(),
            solution_: HighsSolution::default(),
            basis_: HighsBasis::default(),
            info_: HighsInfo::default(),
            run_data_: HighsRunData::default(),
            model_status_: HighsModelStatus::kNotset,
            ekk_instance_: HEkk::default(),
            presolve_: PresolveComponent::default(),
            model_presolve_status_: HighsPresolveStatus::kNotPresolved,
            presolve_log_: HighsPresolveLog::default(),
            model_hessian_: HighsHessian::default(),
            icrash_info_: ICrashInfo::default(),
            ranging_: HighsRanging::default(),
            iis_: HighsIis::default(),
            multi_linear_objective_: Vec::new(),
            presolved_model_lp_: HighsLp::default(),
            saved_objective_and_solution_: Vec::new(),
            mip_stats_: None,
            standard_form_valid_: false,
            standard_form_offset_: 0.0,
            standard_form_cost_: Vec::new(),
            standard_form_rhs_: Vec::new(),
            standard_form_matrix_: HighsSparseMatrix::default(),
            max_threads_: 0,
            capture_mip_probe_: false,
            mip_probe_: None,
        };
        // Highs::Highs() calls initialiseInfo-type setup; info starts invalid
        highs.info_.invalidate();
        highs
    }
}

impl Highs {
    /// A `Highs` with the given options (as after a series of `setOptionValue`).
    pub fn with_options(options: HighsOptions) -> Self {
        let mut highs = Self::default();
        highs.options_ = options;
        highs
    }

    /// HiGHS: Highs.cpp::Highs::passModel(HighsModel) (lines 443-523, 81 lines),
    /// reached through `Highs::passModel(HighsLp)` (lines 525-529). No Hessian.
    pub fn pass_model(&mut self, model_lp: HighsLp) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;
        // Clear the incumbent model and any associated data
        self.clear_model();
        // Move the model's LP to the internal LP
        self.model_lp_ = model_lp;
        let lp = &mut self.model_lp_;
        lp.origin_name_ = "Original".to_string();
        debug_assert!(lp.a_matrix_.format_ok());
        if lp.num_col_ == 0 || lp.num_row_ == 0 {
            // Model constraint matrix has either no columns or no
            // rows. Clearly the matrix is empty, so may have no orientation
            // or starts assigned. HiGHS assumes that such a model will have
            // null starts, so make it column-wise
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Model has either no columns or no rows, so ignoring user constraint matrix data and initialising empty matrix\n"
            );
            lp.a_matrix_.format_ = MatrixFormat::kColwise;
            lp.a_matrix_.start_.clear();
            lp.a_matrix_.start_.resize((lp.num_col_ + 1) as usize, 0);
            lp.a_matrix_.index_.clear();
            lp.a_matrix_.value_.clear();
        } else {
            // Matrix has rows and columns, so a_matrix format must be valid
            if !lp.a_matrix_.format_ok() {
                return HighsStatus::kError;
            }
        }
        // Dimensions in a_matrix_ may not be set, so take them from lp.
        lp.set_matrix_dimensions();
        // Residual scale factors may be present. ToDo Allow user-defined
        // scale factors to be passed
        debug_assert!(!lp.is_scaled_);
        debug_assert!(!lp.is_moved_);
        lp.reset_scale();
        // Check that the LP array dimensions are valid
        if !lp_dimensions_ok("passModel", lp, &self.options_.log_options) {
            return HighsStatus::kError;
        }
        // Check validity of the LP, normalising its values
        return_status = interpret_call_status(
            &self.options_.log_options,
            assess_lp(lp, &self.options_),
            return_status,
            "assessLp",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        // Now legality of matrix is established, ensure that it is
        // column-wise
        lp.ensure_colwise();
        // Clear solver status, solution, basis and info associated with any
        // previous model
        let clear_solver_status = self.clear_solver();
        return_status = interpret_call_status(
            &self.options_.log_options,
            clear_solver_status,
            return_status,
            "clearSolver",
        );
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::clearModel (lines 63-68, 6 lines). Public for
    /// HighsLpRelaxation::loadModel (HighsLpRelaxation.cpp:252).
    pub fn clear_model(&mut self) {
        // model_.clear(): HighsModel::clear (lp_, then hessian_)
        self.model_lp_.clear();
        self.model_hessian_.clear();
        self.multi_linear_objective_.clear();
        self.saved_objective_and_solution_.clear();
        self.clear_solver();
    }

    /// HiGHS: Highs.cpp::Highs::clearSolver (lines 70-76, 7 lines)
    pub fn clear_solver(&mut self) -> HighsStatus {
        let return_status = HighsStatus::kOk;
        self.clear_derived_model_properties();
        self.invalidate_solver_data();
        self.ekk_instance_.clear();
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::invalidateSolverData (lines 3917-3923, 7 lines)
    pub(crate) fn invalidate_solver_data(&mut self) {
        // invalidateSolverDualData
        self.model_status_ = HighsModelStatus::kNotset;
        self.invalidate_ranging();
        self.info_.invalidate();
        self.run_data_.invalidate();
        self.invalidate_solution();
        self.invalidate_basis();
        // invalidateEkk
        self.ekk_instance_.invalidate();
        self.clear_iis();
    }

    /// HiGHS: Highs.cpp::Highs::invalidateSolution (lines 3949-3959, 11 lines)
    pub(crate) fn invalidate_solution(&mut self) {
        self.info_.primal_solution_status = kSolutionStatusNone;
        self.info_.dual_solution_status = kSolutionStatusNone;
        self.info_.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.info_.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
        self.info_.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.info_.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
        self.solution_.invalidate();
    }

    /// HiGHS: Highs.cpp::Highs::invalidateBasis (lines 3961-3964, 4 lines)
    pub(crate) fn invalidate_basis(&mut self) {
        self.info_.basis_validity = kBasisValidityInvalid;
        self.basis_.invalidate();
    }

    /// HiGHS: Highs.cpp::Highs::setHighsModelStatusAndClearSolutionAndBasis (lines 4739-4745, 7 lines)
    pub(crate) fn set_highs_model_status_and_clear_solution_and_basis(
        &mut self,
        model_status: HighsModelStatus,
    ) {
        self.model_status_ = model_status;
        self.invalidate_solution();
        self.invalidate_basis();
        self.info_.valid = true;
    }

    /// HiGHS: HighsInterface.cpp::Highs::zeroIterationCounts (lines 1563-1569, 7 lines)
    fn zero_iteration_counts(&mut self) {
        self.info_.simplex_iteration_count = 0;
        self.info_.ipm_iteration_count = 0;
        self.info_.crossover_iteration_count = 0;
        self.info_.pdlp_iteration_count = 0;
        self.info_.qp_iteration_count = 0;
    }

    /// HiGHS: HighsInterface.cpp::Highs::infeasibleBoundsOk (lines 3747-3860, 114 lines).
    /// HiGHS 1.16 (#3576): HighsInterface.cpp::processedBoundsFindsInfeasibility (lines 3756-3877)
    fn infeasible_bounds_ok(&mut self) -> bool {
        let log_options = &self.options_.log_options;
        let primal_feasibility_tolerance = self.options_.primal_feasibility_tolerance;
        let mip_feasibility_tolerance = self.options_.mip_feasibility_tolerance;
        let perform_inward_integer_rounding = !self.options_.solve_relaxation;
        // HiGHS 1.16 #3576 (h16_int_round_assert_tol): read the flag once, before `lp` borrows model_lp_
        let h16_int_round_assert_tol = self.options_.h16_int_round_assert_tol;
        let lp = &mut self.model_lp_;

        let mut num_true_infeasible_bound: HighsInt = 0;
        let mut num_ok_infeasible_bound: HighsInt = 0;
        // C++ `has_integrality = lp.integrality_.size() > 0`
        let has_integrality = !lp.integrality_.is_empty();
        // Lambda for assessing infeasible bounds. C++ captures `performed_inward_integer_rounding`
        // by reference and resets it per column; the port passes it as the last argument.
        let mut infeasible_bound_ok = |type_: &str,
                                       i_x: HighsInt,
                                       lower: &mut f64,
                                       upper: &mut f64,
                                       performed_inward_integer_rounding: bool|
         -> bool {
            let range = *upper - *lower;
            // Should only be called if lower > upper, so range < 0
            debug_assert!(range < 0.0);
            if range >= 0.0 {
                return true;
            }
            if range > -primal_feasibility_tolerance {
                // Infeasibility is less than feasibility tolerance, so fix
                // bounds at lower (upper) if lower (upper) is an integer - and
                // both can't be integer, otherwise the range <= -1 - otherwise
                // fix at 0.5 * (lower + upper)
                num_ok_infeasible_bound += 1;
                let report = num_ok_infeasible_bound <= 10;
                let integer_lower = *lower == (*lower + 0.5).floor();
                let integer_upper = *upper == (*upper + 0.5).floor();
                debug_assert!(!integer_lower || !integer_upper);
                if integer_lower {
                    if report {
                        crate::highs_log_user!(
                            log_options,
                            HighsLogType::kInfo,
                            "{} {} bounds [{}, {}] have infeasibility = {} so set upper bound to {}\n",
                            type_,
                            i_x,
                            *lower,
                            *upper,
                            range,
                            *lower
                        );
                    }
                    *upper = *lower;
                } else if integer_upper {
                    if report {
                        crate::highs_log_user!(
                            log_options,
                            HighsLogType::kInfo,
                            "{} {} bounds [{}, {}] have infeasibility = {} so set lower bound to {}\n",
                            type_,
                            i_x,
                            *lower,
                            *upper,
                            range,
                            *upper
                        );
                    }
                    *lower = *upper;
                } else {
                    let mid = 0.5 * (*lower + *upper);
                    if report {
                        crate::highs_log_user!(
                            log_options,
                            HighsLogType::kInfo,
                            "{} {} bounds [{}, {}] have infeasibility = {} so set both bounds to {}\n",
                            type_,
                            i_x,
                            *lower,
                            *upper,
                            range,
                            mid
                        );
                    }
                    *lower = mid;
                    *upper = mid;
                }
                return true;
            }
            // Infeasibility is greater than feasibility tolerance, so report
            // this (up to 10 times)
            num_true_infeasible_bound += 1;
            if num_true_infeasible_bound <= 10 {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kInfo,
                    "{} {} bounds [{}, {}] have excessive infeasibility = {}{}\n",
                    type_,
                    i_x,
                    *lower,
                    *upper,
                    range,
                    if performed_inward_integer_rounding {
                        // HiGHS 1.16 #3576 (h16_int_round_assert_tol): HighsInterface.cpp lines 3813-3816 (1.16.0)
                        if h16_int_round_assert_tol {
                            " after inward integer rounding"
                        } else {
                            " due to inward integer rounding"
                        }
                    } else {
                        ""
                    }
                );
            }
            false
        };

        for i_col in 0..lp.num_col_ {
            let ic = i_col as usize;
            let mut performed_inward_integer_rounding = false;
            let mut lower = lp.col_lower_[ic];
            let mut upper = lp.col_upper_[ic];
            if has_integrality {
                // Semi-variables cannot have inconsistent bounds
                if lp.integrality_[ic] == HighsVarType::kSemiContinuous
                    || lp.integrality_[ic] == HighsVarType::kSemiInteger
                {
                    continue;
                }
                if perform_inward_integer_rounding && lp.integrality_[ic] == HighsVarType::kInteger
                {
                    // Assess bounds after inward integer rounding
                    let feastol = mip_feasibility_tolerance;
                    let integer_lower = (lower - feastol).ceil();
                    let integer_upper = (upper + feastol).floor();
                    // HiGHS 1.16 #3576 (h16_int_round_assert_tol): HighsInterface.cpp lines 3836-3837 (1.16.0)
                    if h16_int_round_assert_tol {
                        debug_assert!(integer_lower >= lower - feastol);
                        debug_assert!(integer_upper <= upper + feastol);
                    } else {
                        debug_assert!(integer_lower >= lower);
                        debug_assert!(integer_upper <= upper);
                    }
                    performed_inward_integer_rounding =
                        integer_lower > lower || integer_upper < upper;
                    lower = integer_lower;
                    upper = integer_upper;
                }
            }
            //
            if lower > upper {
                if infeasible_bound_ok(
                    "Column",
                    i_col,
                    &mut lower,
                    &mut upper,
                    performed_inward_integer_rounding,
                ) {
                    // Bound infeasibility is OK (less than the tolerance), so can
                    // change the model data. Inward integer rounding is not written
                    // back: it may be a significant change (C++ comment, cpp 3852-3857).
                    lp.col_lower_[ic] = lower;
                    lp.col_upper_[ic] = upper;
                }
            }
        }
        for i_row in 0..lp.num_row_ {
            let ir = i_row as usize;
            if lp.row_lower_[ir] > lp.row_upper_[ir] {
                let mut lower = lp.row_lower_[ir];
                let mut upper = lp.row_upper_[ir];
                infeasible_bound_ok("Row", i_row, &mut lower, &mut upper, false);
                lp.row_lower_[ir] = lower;
                lp.row_upper_[ir] = upper;
            }
        }
        drop(infeasible_bound_ok);
        if num_ok_infeasible_bound > 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "Model has {} small inconsistent bound(s): rectified\n",
                num_ok_infeasible_bound
            );
        }
        if num_true_infeasible_bound > 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "Model has {} significant inconsistent bound(s): infeasible\n",
                num_true_infeasible_bound
            );
        }
        num_true_infeasible_bound == 0
    }

    /// HiGHS: Highs.cpp::Highs::callSolveLp (lines 4122-4137, 16 lines)
    pub(crate) fn call_solve_lp(&mut self, message: &str) -> HighsStatus {
        // Check that the model is column-wise
        debug_assert!(self.model_lp_.a_matrix_.is_colwise());
        let mut solver_object = HighsLpSolverObject::new(
            &mut self.model_lp_,
            &mut self.basis_,
            &mut self.solution_,
            &mut self.info_,
            &mut self.callback_,
            &mut self.options_,
            &mut self.timer_,
        );
        // Solve the LP
        let return_status = solve_lp(&mut self.ekk_instance_, &mut solver_object, message);
        // Extract the model status
        self.model_status_ = solver_object.model_status_;
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::callLpKktCheck(model_.lp_) (lines
    /// 2948-2951) -> HighsSolution.cpp::lpKktCheck (presolve_support.rs
    /// `lp_kkt_check`).
    pub(crate) fn call_lp_kkt_check(&mut self) {
        self.call_lp_kkt_check_target(
            crate::lp_data::highs_presolve::HighsSolveTarget::Incumbent,
            "",
        );
    }

    /// HiGHS: HighsInterface.cpp::Highs::reportModelStats (lines 25-125)
    ///
    /// PORT NOTE: the port's model has no Hessian, so `q_num_nz` is 0 and the problem
    /// type is "MIP" or "LP".
    pub fn report_model_stats(&self) {
        let lp = &self.model_lp_;
        let log_options = &self.options_.log_options;
        if !log_options.output_flag {
            return;
        }
        let mut num_integer: HighsInt = 0;
        let mut num_binary: HighsInt = 0;
        let mut num_semi_continuous: HighsInt = 0;
        let mut num_semi_integer: HighsInt = 0;
        for i_col in 0..lp.integrality_.len() {
            match lp.integrality_[i_col] {
                HighsVarType::kInteger => {
                    num_integer += 1;
                    if lp.col_lower_[i_col] == 0.0 && lp.col_upper_[i_col] == 1.0 {
                        num_binary += 1;
                    }
                }
                HighsVarType::kSemiContinuous => num_semi_continuous += 1,
                HighsVarType::kSemiInteger => num_semi_integer += 1,
                _ => {}
            }
        }
        let non_continuous = num_integer + num_semi_continuous + num_semi_integer != 0;
        let problem_type = if non_continuous { "MIP" } else { "LP" };
        let a_num_nz = lp.a_matrix_.num_nz();
        let q_num_nz: HighsInt = 0;
        let s = |n: HighsInt| if n == 1 { "" } else { "s" };
        if log_options.log_dev_level != 0 {
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kInfo,
                "{:>4}      : {}\n",
                problem_type,
                lp.model_name_
            );
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kInfo,
                "Row{}      : {}\n",
                s(lp.num_row_),
                lp.num_row_
            );
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kInfo,
                "Col{}      : {}\n",
                s(lp.num_col_),
                lp.num_col_
            );
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kInfo,
                "Nonzero{}  : {}\n",
                s(a_num_nz),
                a_num_nz
            );
            if num_integer != 0 {
                crate::highs_log_dev!(
                    log_options,
                    HighsLogType::kInfo,
                    "Integer   : {} ({} binary)\n",
                    num_integer,
                    num_binary
                );
            }
            if num_semi_continuous != 0 {
                crate::highs_log_dev!(
                    log_options,
                    HighsLogType::kInfo,
                    "SemiConts : {}\n",
                    num_semi_continuous
                );
            }
            if num_semi_integer != 0 {
                crate::highs_log_dev!(
                    log_options,
                    HighsLogType::kInfo,
                    "SemiInt   : {}\n",
                    num_semi_integer
                );
            }
        } else {
            let mut stats_line = problem_type.to_string();
            if !lp.model_name_.is_empty() {
                stats_line += &format!(" {}", lp.model_name_);
            }
            stats_line += &format!(
                " has {} row{}; {} col{}",
                lp.num_row_,
                s(lp.num_row_),
                lp.num_col_,
                s(lp.num_col_)
            );
            let _ = q_num_nz;
            stats_line += &format!("; {} nonzero{}", a_num_nz, s(a_num_nz));
            if num_integer != 0 {
                // HiGHS quirk: the plural of "integer variable" tests a_num_nz
                stats_line += &format!(
                    "; {} integer variable{} ({} binary)",
                    num_integer,
                    s(a_num_nz),
                    num_binary
                );
            }
            if num_semi_continuous != 0 {
                stats_line += &format!("; {} semi-continuous variables", num_semi_continuous);
            }
            if num_semi_integer != 0 {
                stats_line += &format!("; {} semi-integer variables", num_semi_integer);
            }
            crate::highs_log_user!(log_options, HighsLogType::kInfo, "{}\n", stats_line);
        }
    }

    /// HiGHS: Highs.cpp::Highs::run (lines 1088-1151) -> optimizeModel (lines
    /// 1153-1171) -> calledOptimizeModel (lines 1173-2059), the LP path with
    /// presolve off and no useful basis. See the module doc for what is out
    /// of scope.
    pub fn run(&mut self) -> HighsStatus {
        let _run_depth = RunDepth::enter();
        // HiGHS: Highs.cpp::Highs::run (lines 1088-1114), level 0: the read/write file
        // options are not ported (always ""), so only the model report and the user
        // scaling assessment run before optimizeModel.
        self.report_model_stats();
        // Apply any user objective and bound scaling
        let mut user_scale_data = HighsUserScaleData::default();
        if self.user_scale(&mut user_scale_data) == HighsStatus::kError {
            return HighsStatus::kError;
        }
        // Determine coefficient ranges and possibly warn the user about
        // excessive values, obtaining suggested values for user_objective_scale
        // and user_bound_scale
        // HiGHS 1.16 (h16_scaling_report): Highs.cpp line 1238 (1.16.0) makes the matrix
        // column-wise first
        if self.options_.h16_scaling_report {
            self.model_lp_.a_matrix_.ensure_colwise();
        }
        assess_excessive_objective_bound_scaling(
            &self.options_,
            &self.model_lp_,
            &mut user_scale_data,
        );
        // Optimize the model in the Highs instance
        // rusolve: the fold attempt and the fallback solve used
        // to each get the full time limit (bdry2: TL 60 s -> 130.6 s, killed by the harness at
        // TL 300 s). The fallback now gets what the fold attempt left.
        let fold_start = std::time::Instant::now();
        let user_limit = self.options_.time_limit;
        let status = match self.try_lp_concurrent().or_else(|| self.try_lp_fold()) {
            Some(s) => s,
            None => {
                if user_limit.is_finite() {
                    let spent = fold_start.elapsed().as_secs_f64();
                    self.options_.time_limit = (user_limit - spent).max(0.0);
                }
                // The fold attempt used the whole limit (dlr1: the quotient solve ran to the
                // limit, then presolve of 1.7M rows and the setup added 6 s): nothing is left
                // for a fallback solve.
                if user_limit.is_finite()
                    && user_limit > 0.0
                    && self.options_.time_limit <= 0.0
                    && self.options_.lp_fold
                {
                    self.options_.time_limit = user_limit;
                    self.model_status_ = HighsModelStatus::kTimeLimit;
                    return HighsStatus::kWarning;
                }
                // rusolve lever lp_solver_auto: after the fold has
                // declined, so the fold's quotient keeps the default solver
                let auto_start = std::time::Instant::now();
                let auto_checked = self.options_.lp_solver_auto
                    && self.options_.time_limit > 0.0
                    && self.options_.threads < 2
                    && self.options_.solver == kHighsChooseString
                    && !self.model_lp_.is_mip();
                let auto_solver = auto_checked && self.lp_auto_picks_ipm();
                // the rule's symbolic analysis is part of the user's time budget, whichever
                // way it decided
                if auto_checked && self.options_.time_limit.is_finite() {
                    let t = self.options_.time_limit - auto_start.elapsed().as_secs_f64();
                    self.options_.time_limit = t.max(0.0);
                }
                let (saved_solver, saved_hipo) =
                    (self.options_.solver.clone(), self.options_.hipo_available);
                if auto_solver {
                    self.options_.solver = kIpmString.to_string();
                    self.options_.hipo_available = true;
                }
                let s = self.optimize_highs();
                self.options_.solver = saved_solver;
                self.options_.hipo_available = saved_hipo;
                self.options_.time_limit = user_limit;
                s
            }
        };
        if status == HighsStatus::kError {
            return status;
        }
        // Assume that run data are valid
        self.run_data_.valid = true;
        //
        // Undo any user objective and bound scaling
        if self.user_unscale(&mut user_scale_data) == HighsStatus::kError {
            return HighsStatus::kError;
        }
        // The write_iis_model_file / solution_file / write_basis_file writes are not
        // ported (the file options are always "")
        status
    }

    /// HiGHS: Highs.cpp::Highs::calledOptimizeModel (lines 1173-2059), the part of
    /// [`Highs::run`] below level 0. C++ `Highs::optimizeLp` (what `HighsLpRelaxation`
    /// calls) enters here directly, skipping `reportModelStats` and
    /// `assessExcessiveObjectiveBoundScaling`.
    pub fn called_optimize_model(&mut self) -> HighsStatus {
        let min_highs_debug_level = kHighsDebugLevelMin;
        if self.options_.highs_debug_level < min_highs_debug_level {
            self.options_.highs_debug_level = min_highs_debug_level;
        }
        if !self.options_.use_warm_start {
            self.clear_solver();
        }
        // HiGHS solvers require models with no infinite costs, and no semi-variables
        //
        // Set undo_mods = false so that returnFromOptimizeModel() doesn't undo any
        // mods that must be preserved - such as when solving a MIP node
        let mut undo_mods = false;
        if self.model_lp_.has_infinite_cost_ {
            // If the model has infinite costs, then try to remove them. The
            // return_status indicates the success of this operation and, if
            // it's unsuccessful, the model will not have been modified and
            // optimizeModel() can simply return an error with model status
            // HighsModelStatus::kUnknown
            let return_status = self.handle_inf_cost();
            if return_status != HighsStatus::kOk {
                debug_assert!(return_status == HighsStatus::kError);
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kUnknown,
                );
                return return_status;
            }
            // Modifications have been performed, so must be undone before
            // this call to optimizeModel() returns
            debug_assert!(!self.model_lp_.has_infinite_cost_);
            undo_mods = true;
        }
        // Ensure that all vectors in the model have exactly the right size
        self.model_lp_.exact_resize();

        if self.model_lp_.is_mip() && self.solution_.value_valid {
            // Determine whether the current (partial) solution of a MIP is
            // feasible and, if not, try to complete the assignment with
            // integer values (if necessary) and continuous values (if
            // necessary) to achieve a feasible solution. Valuable in the case
            // where users make a heuristic (partial) assignment of discrete variables
            let call_status = self.complete_solution_from_discrete_assignment();
            if call_status != HighsStatus::kOk {
                return HighsStatus::kError;
            }
        }

        let mut return_status = HighsStatus::kOk;
        // Initialise the HiGHS model status
        self.model_status_ = HighsModelStatus::kNotset;
        // Clear the run info and data
        self.info_.invalidate();
        self.run_data_.invalidate();
        // Zero the iteration counts
        self.zero_iteration_counts();
        // Start the HiGHS run clock
        self.timer_.start();
        // Return immediately if the model has no columns
        if self.model_lp_.num_col_ == 0 {
            self.set_highs_model_status_and_clear_solution_and_basis(HighsModelStatus::kModelEmpty);
            return self.return_from_optimize_model(HighsStatus::kOk, undo_mods);
        }
        // Return immediately if the model is infeasible due to inconsistent
        // bounds, modifying any bounds with tiny infeasibilities
        if !self.infeasible_bounds_ok() {
            self.set_highs_model_status_and_clear_solution_and_basis(HighsModelStatus::kInfeasible);
            return self.return_from_optimize_model(return_status, undo_mods);
        }
        // Ensure that the LP (and any simplex LP) has the matrix column-wise
        self.model_lp_.ensure_colwise();
        // Ensure that the matrix has no large values
        if self
            .model_lp_
            .a_matrix_
            .has_large_value(self.options_.large_matrix_value)
        {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot solve a model with a |value| exceeding {} in constraint matrix\n",
                self.options_.large_matrix_value
            );
            return self.return_from_optimize_model(HighsStatus::kError, undo_mods);
        }
        if self.options_.highs_debug_level > min_highs_debug_level {
            let options = self.options_.clone();
            let call_status = assess_lp(&mut self.model_lp_, &options);
            debug_assert!(call_status == HighsStatus::kOk);
            return_status = interpret_call_status(
                &self.options_.log_options,
                call_status,
                return_status,
                "assessLp",
            );
            if return_status == HighsStatus::kError {
                return self.return_from_optimize_model(return_status, undo_mods);
            }
        }
        if !self.options_.solve_relaxation {
            // Not solving the relaxation, so check validity of any
            // integrality, keeping a record of any bound and type
            // modifications for semi-variables
            //
            // C++ calls assessSemiVariables whenever !solve_relaxation; without
            // semi-variables it only counts and may warn "No semi-integer/integer
            // variables ...".
            let mut made_semi_variable_mods = false;
            // PORT NOTE (perf): disjoint field borrows; no options copy per call.
            let call_status = assess_semi_variables(
                &mut self.model_lp_,
                &self.options_,
                &mut made_semi_variable_mods,
            );
            undo_mods = undo_mods || made_semi_variable_mods;
            if call_status == HighsStatus::kError {
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kSolveError,
                );
                return self.return_from_optimize_model(HighsStatus::kError, undo_mods);
            }
        }
        // Choose method according to model class (lines 1347-1390)
        //
        // (lines 1347-1372) the QP branch, Reached only with a
        // Hessian (`model_.isQp()`), so the LP and MIP paths never enter it.
        if self.model_hessian_.dim_ != 0 {
            if self.model_lp_.is_mip() {
                if self.options_.solve_relaxation {
                    // Relax any semi-variables (without semi-variables
                    // relaxSemiVariables changes nothing)
                    let mut made_semi_variable_mods = false;
                    relax_semi_variables(&mut self.model_lp_, &mut made_semi_variable_mods);
                    undo_mods = undo_mods || made_semi_variable_mods;
                } else {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Cannot solve MIQP problems with HiGHS\n"
                    );
                    return self.return_from_optimize_model(HighsStatus::kError, undo_mods);
                }
            }
            // Ensure that its diagonal entries are OK in the context of the
            // objective sense. It's OK to be semi-definite
            if !ok_hessian_diagonal(
                &self.options_,
                &mut self.model_hessian_,
                self.model_lp_.sense_,
            ) {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Cannot solve non-convex QP problems with HiGHS\n"
                );
                return self.return_from_optimize_model(HighsStatus::kError, undo_mods);
            }
            // solverValidForQp / warnSolverInvalid: not ported (as for LP/MIP)
            let call_status = self.call_solve_qp();
            return_status = interpret_call_status(
                &self.options_.log_options,
                call_status,
                return_status,
                "callSolveQp",
            );
            return self.return_from_optimize_model(return_status, undo_mods);
        }
        // (lines 1373-1390) the MIP branch
        if self.model_lp_.is_mip() {
            if self.options_.solve_relaxation {
                // Solving just the relaxation, so relax any semi-variables
                let mut made_semi_variable_mods = false;
                relax_semi_variables(&mut self.model_lp_, &mut made_semi_variable_mods);
                undo_mods = undo_mods || made_semi_variable_mods;
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Solving LP relaxation since solve_relaxation is true\n"
                );
            } else {
                // Solve model as a MIP (solverValidForMip/warnSolverInvalid: the
                // CLI passes the default solver)
                let call_status = self.call_solve_mip();
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "callSolveMip",
                );
                return self.return_from_optimize_model(return_status, undo_mods);
            }
        }
        // Solve the model as an LP
        //
        // (lines 1408-1411)
        let ipm_no_crossover =
            use_ipm(&self.options_.solver) && self.options_.run_crossover == kHighsOffString;
        let lp_no_solution_basis = ipm_no_crossover || self.options_.solver == kPdlpString;
        // (lines 1414-1462) iCrash, Gated on options_.icrash
        // (default false), so the LP path never enters it.
        if self.options_.icrash {
            let mut strategy = ICrashStrategy::kICA;
            let strategy_ok = crate::presolve::icrash::parse_icrash_strategy(
                &self.options_.icrash_strategy,
                &mut strategy,
            );
            if !strategy_ok {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "ICrash error: unknown strategy.\n"
                );
                return HighsStatus::kError;
            }
            let icrash_options = ICrashOptions {
                dualize: self.options_.icrash_dualize,
                strategy,
                starting_weight: self.options_.icrash_starting_weight,
                iterations: self.options_.icrash_iterations,
                approximate_minimization_iterations: self.options_.icrash_approx_iter,
                exact: self.options_.icrash_exact,
                breakpoints: self.options_.icrash_breakpoints,
                log_options: self.options_.log_options.clone(),
            };
            let icrash_status = crate::presolve::icrash::call_icrash(
                &self.model_lp_,
                &icrash_options,
                &mut self.icrash_info_,
            );
            if icrash_status != HighsStatus::kOk {
                return self.return_from_optimize_model(icrash_status, undo_mods);
            }
            // for now set the solution_.col_value
            self.solution_.col_value = self.icrash_info_.x_values.clone();
            // Better not to use Highs::crossover (use_highs_crossover = false)
            let crossover_status = crate::presolve::icrash_x::call_crossover(
                &self.options_,
                &self.model_lp_,
                &mut self.basis_,
                &mut self.solution_,
                &mut self.model_status_,
                &mut self.info_,
                &mut self.callback_,
            );
            // callCrossover can return HighsStatus::kWarning due to
            // imprecise dual values. Ignore this since primal simplex will
            // be called to clean up duals
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Crossover following iCrash has return status of {}, and problem status is {}\n",
                highs_status_to_string(crossover_status),
                self.model_status_to_string(self.model_status_)
            );
            if crossover_status == HighsStatus::kError {
                return self.return_from_optimize_model(crossover_status, undo_mods);
            }
            debug_assert!(self.options_.simplex_strategy == kSimplexStrategyPrimal);
        }
        let mut no_incumbent_lp_solution_or_basis = false;
        // Even if options_.solver == kHighsChooseString in isolation will,
        // ultimately lead to a choice between simplex and IPM, if a basis
        // is available, simplex should surely be chosen.
        let solver_will_use_basis =
            self.options_.solver == kSimplexString || self.options_.solver == kHighsChooseString;
        if solver_will_use_basis {
            if !self.basis_.valid && self.solution_.value_valid {
                // There is no valid basis, but there is a valid solution, so use
                // it to construct a basis
                let call_status = self.basis_for_solution();
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "basisForSolution",
                );
                if return_status == HighsStatus::kError {
                    return self.return_from_optimize_model(return_status, undo_mods);
                }
                debug_assert!(self.basis_.valid);
            }
        } else {
            // The basis won't be used, so clear it to ensure that, after any
            // presolve the solver choice won't be over-ruled by choosing
            // simplex due to the existence of a basis - which mustn't be used
            // after a strict reduction due to presolve!
            self.basis_.clear();
        }
        let unconstrained_lp = self.model_lp_.a_matrix_.num_nz() == 0;
        debug_assert!(self.model_lp_.num_row_ != 0 || unconstrained_lp);
        let has_basis = self.basis_.useful;
        if has_basis {
            debug_assert!(self.basis_.col_status.len() == self.model_lp_.num_col_ as usize);
            debug_assert!(self.basis_.row_status.len() == self.model_lp_.num_row_ as usize);
        }
        if self.basis_.valid {
            debug_assert!(self.basis_.useful);
        }
        let without_presolve = self.options_.presolve == kHighsOffString;
        if (unconstrained_lp || has_basis || without_presolve) && solver_will_use_basis {
            // There is a valid basis for the problem, presolve is off, or LP
            // has no constraint matrix, and the solver will use the basis
            //
            // Determine a coherent message about how the LP is being solved
            let lp_solve = if unconstrained_lp {
                "Solving unconstrained LP"
            } else if has_basis {
                if without_presolve {
                    "Solving LP with useful basis"
                } else {
                    "Solving LP with useful basis so presolve not used"
                }
            } else {
                // One of unconstrained_lp, has_basis and without_presolve
                // must be true, and the first two aren't
                debug_assert!(without_presolve);
                "Solving LP without presolve or useful basis"
            };
            self.ekk_instance_.lp_name_ = lp_solve.to_string();
            // If there is a valid HiGHS basis, refine any status values that
            // are simply HighsBasisStatus::kNonbasic
            if self.basis_.useful {
                refine_basis(&self.model_lp_, &self.solution_, &mut self.basis_);
            }
            let call_status = self.call_solve_lp(lp_solve);
            return_status = interpret_call_status(
                &self.options_.log_options,
                call_status,
                return_status,
                "callSolveLp",
            );
            if return_status == HighsStatus::kError {
                return self.return_from_optimize_model(return_status, undo_mods);
            }
        } else {
            // Otherwise, consider presolve (highs_presolve.rs)
            match self.called_optimize_model_presolve(
                return_status,
                lp_no_solution_basis,
                undo_mods,
            ) {
                Ok(no_incumbent) => no_incumbent_lp_solution_or_basis = no_incumbent,
                Err(early_return_status) => return early_return_status,
            }
        }
        // Unless the model status was determined using the strictly reduced
        // LP, the HiGHS info is valid
        if !no_incumbent_lp_solution_or_basis {
            self.call_lp_kkt_check();
            self.info_.valid = true;
        }
        // Assess success according to the model status, regardless of
        // whether anything worse has happened earlier
        return_status = highs_status_from_highs_model_status(self.model_status_);
        self.return_from_optimize_model(return_status, undo_mods)
    }

    /// HiGHS: Highs.cpp::Highs::callSolveMip (lines 4327-4448): copy the model and options into a `HighsMipSolver`,
    /// run it, and copy back status, solution, info and the MIP bounds.
    ///
    /// Ownership: `HighsMipSolver solver(callback_,
    /// options_, lp, solution_)` -> `let mut solver = HighsMipSolver::new(&mut
    /// self.callback_, &self.options_, lp, &self.solution_, false, 0)` (moves the callback
    /// in), `solver.set_profiling()`, `solver.run()`, then **first**
    /// `self.callback_ = std::mem::take(&mut solver.callback_)`, then the C++ copies in
    /// order. `lp` is `&self.model_lp_`, or the `without_semi_variables` copy when
    /// the model has semi-variables.
    pub(crate) fn call_solve_mip(&mut self) -> HighsStatus {
        // Record whether there is a valid primal solution on entry
        let user_solution = self.solution_.value_valid;
        let saved_user_solution = if user_solution {
            // Save the col and row values
            Some((
                std::mem::take(&mut self.solution_.col_value),
                std::mem::take(&mut self.solution_.row_value),
            ))
        } else {
            None
        };
        // Ensure that any solver data for users in Highs class members are cleared
        self.invalidate_solver_data();
        if let Some((col_value, row_value)) = saved_user_solution {
            // Recover the col and row values
            self.solution_.col_value = col_value;
            self.solution_.row_value = row_value;
            self.solution_.value_valid = true;
        }
        // rusolve: independent components solved separately
        if let Some(st) = self.call_solve_mip_decomposed() {
            return st;
        }
        // rusolve: small pure-integer models in reduced lattice coordinates
        if let Some(st) = self.call_solve_mip_lattice() {
            return st;
        }
        // Run the MIP solver
        let log_dev_level = self.options_.log_dev_level;
        //  options_.log_dev_level = kHighsLogDevLevelInfo;
        // Check that the model isn't row-wise
        debug_assert!(self.model_lp_.a_matrix_.format_ != MatrixFormat::kRowwise);
        let has_semi_variables = self.model_lp_.has_semi_variables();
        // HighsLp& lp = has_semi_variables ? use_lp : model_.lp_, with
        // use_lp = withoutSemiVariables(model_.lp_, solution_,
        // options_.primal_feasibility_tolerance), which extends solution_
        let use_lp = if has_semi_variables {
            Some(without_semi_variables(
                &self.model_lp_,
                &mut self.solution_,
                self.options_.primal_feasibility_tolerance,
            ))
        } else {
            None
        };
        let lp = use_lp.as_ref().unwrap_or(&self.model_lp_);
        // rusolve mip_activation_dff: the solver sees the model
        // plus implied rows; the solution's row values are recomputed from the user model below.
        let dff_lp = if self.options_.mip_activation_dff {
            let mut aug = lp.clone();
            let k = crate::mip::activation_dff::add_activation_dff_rows(&mut aug);
            if k > 0 {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Activation DFF rows added before presolve: {}\n",
                    k
                );
            }
            (k > 0).then_some(aug)
        } else {
            None
        };
        let lp = dff_lp.as_ref().unwrap_or(lp);
        // HighsMipSolver solver(callback_, options_, lp, solution_); the callback is
        // moved in and moved back below before any other use of `solver`.
        let racing = self.options_.mip_racing_threads > 1;
        let mut solver = if racing {
            // rusolve: diversified copies race and share incumbents
            race_solve_mip(
                &mut self.callback_,
                &self.options_,
                lp,
                &self.solution_,
                self.options_.mip_racing_threads as usize,
            )
        } else {
            let mut solver = HighsMipSolver::new(
                &mut self.callback_,
                &self.options_,
                lp,
                &self.solution_,
                false,
                0,
            );
            // PORT NOTE: profiling_->start/stop(kSubSolverMip) are no-op stubs.
            solver.set_profiling();
            solver.run();
            solver
        };
        if self.capture_mip_probe_ {
            self.mip_probe_ = Some(solver.probe_snapshot());
        }
        self.mip_stats_ = Some(solver.collect_stats()); // read-only
        self.callback_ = std::mem::take(&mut solver.callback_);
        self.options_.log_dev_level = log_dev_level;
        // Set the return_status and model status
        let mut return_status = highs_status_from_highs_model_status(solver.modelstatus_);
        self.model_status_ = solver.modelstatus_;
        // Extract the solution
        if solver.solution_objective_ != kHighsInf {
            // There is a primal solution
            //
            // If the original model has semi-variables, its solution is
            // (still) given by the first model_.lp_.num_col_ entries of the
            // solution from the MIP solver
            //
            // HiGHS quirk (#2547): the resize to model_.lp_.num_col_ is commented
            // out in C++, so with semi-variables col_value keeps the extra entries.
            self.solution_.col_value = std::mem::take(&mut solver.solution_);
            self.saved_objective_and_solution_ =
                std::mem::take(&mut solver.saved_objective_and_solution_);
            self.model_lp_.a_matrix_.product_quad(
                &mut self.solution_.row_value,
                &self.solution_.col_value,
                kDebugReportOff,
            );
            self.solution_.value_valid = true;
        }
        // HiGHS quirk: the C++ `assert(!solution_.value_valid)` here is not reproduced:
        // a user solution restored above is still valid when the MIP finds none.
        //
        // Check that no modified upper bounds for semi-variables are active
        if self.solution_.value_valid
            && active_modified_upper_bounds(
                &self.options_,
                &self.model_lp_,
                &self.solution_.col_value,
            )
        {
            self.solution_.value_valid = false;
            self.model_status_ = HighsModelStatus::kSolveError;
            return_status = HighsStatus::kError;
        }
        // There is no dual solution or basis: invalidateSolverData cleared them.
        debug_assert!(!self.solution_.dual_valid);
        debug_assert!(!self.basis_.valid);
        // Get the objective and any KKT failures
        self.info_.objective_function_value = solver.solution_objective_;
        // Remember to judge primal feasibility according to
        // mip_feasibility_tolerance, so take a copy of the original value...
        let primal_feasibility_tolerance = self.options_.primal_feasibility_tolerance;
        self.options_.primal_feasibility_tolerance = self.options_.mip_feasibility_tolerance;
        // NB getKktFailures sets the primal and dual solution status. The model overload
        // of C++ getKktFailures also runs getPrimalDualBasisErrors and
        // getPrimalDualGlpsolErrors, which write only the discarded primal_dual_errors.
        get_kkt_failures(
            &self.options_,
            false,
            &self.model_lp_,
            &self.model_lp_.col_cost_,
            &self.solution_,
            &mut self.info_,
            false,
        );
        // Set the MIP-specific values of info_
        self.info_.mip_node_count = solver.node_count_;
        self.info_.mip_dual_bound = solver.dual_bound_;
        self.info_.mip_gap = solver.gap_;
        self.info_.primal_dual_integral = solver.primal_dual_integral_;
        // Get the number of LP iterations, avoiding overflow if the int64_t value is too large
        let mip_total_lp_iterations = solver.total_lp_iterations_;
        self.info_.simplex_iteration_count = if mip_total_lp_iterations > i64::from(kHighsIInf) {
            -1
        } else {
            mip_total_lp_iterations as HighsInt
        };
        self.info_.valid = true;
        if self.model_status_ == HighsModelStatus::kOptimal {
            return_status = self.check_optimality("MIP");
        }
        // Overwrite max infeasibility to include integrality if there is a solution
        if solver.solution_objective_ != kHighsInf {
            let mip_max_bound_violation = cpp_max(solver.row_violation_, solver.bound_violation_);
            let delta_max_bound_violation =
                (mip_max_bound_violation - self.info_.max_primal_infeasibility).abs();
            // Possibly report a mis-match between the max bound violation returned by the
            // MIP solver, and the value obtained from the solution
            if delta_max_bound_violation > 1e-12 {
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kWarning,
                    "Inconsistent max bound violation: MIP solver ({:10.4}); LP ({:10.4}); Difference of {:10.4}\n",
                    mip_max_bound_violation,
                    self.info_.max_primal_infeasibility,
                    delta_max_bound_violation
                );
            }
            self.info_.max_integrality_violation = solver.integrality_violation_;
            if self.info_.max_integrality_violation > self.options_.mip_feasibility_tolerance {
                self.info_.primal_solution_status = kSolutionStatusInfeasible;
                debug_assert!(self.model_status_ == HighsModelStatus::kInfeasible);
            }
        }
        // ... and remember to recover the primal feasibility tolerance
        self.options_.primal_feasibility_tolerance = primal_feasibility_tolerance;
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::checkOptimality (lines 2902-2952, 51 lines).
    /// Private to this module: only `call_solve_mip` reaches it in this port.
    pub(crate) fn check_optimality(&mut self, solver_type: &str) -> HighsStatus {
        debug_assert!(self.model_status_ == HighsModelStatus::kOptimal);
        // Cannot expect to have no dual_infeasibilities since the QP solver
        // (and, of course, the MIP solver) give no dual information
        if self.info_.num_primal_infeasibilities == 0 && self.info_.num_dual_infeasibilities <= 0 {
            // Consider semi-continuous infeasibilities
            if self.info_.num_semi_infeasibilities > 0 {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "{} solver claims optimality, but with num/max/sum {}/{}/{} semi-variable infeasibilities: consider solving with smaller mip_feasibility_tolerance\n",
                    solver_type,
                    self.info_.num_semi_infeasibilities,
                    self.info_.max_semi_infeasibility,
                    self.info_.sum_semi_infeasibilities
                );
                self.model_status_ = HighsModelStatus::kSolveError;
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Setting model status to {}\n",
                    util_model_status_to_string(self.model_status_)
                );
                return HighsStatus::kError;
            }
            return HighsStatus::kOk;
        }
        self.model_status_ = HighsModelStatus::kSolveError;
        let mut report_string = format!(
            "{} solver claims optimality, but with num/max/sum primal({}/{}/{})",
            solver_type,
            self.info_.num_primal_infeasibilities,
            self.info_.max_primal_infeasibility,
            self.info_.sum_primal_infeasibilities
        );
        // HiGHS quirk: the "and dual" clause has no leading space in C++, so the report
        // reads "primal(a/b/c)and dual(d/e/f)".
        if self.info_.num_dual_infeasibilities > 0 {
            report_string.push_str(&format!(
                "and dual({}/{}/{})",
                self.info_.num_dual_infeasibilities,
                self.info_.max_dual_infeasibility,
                self.info_.sum_dual_infeasibilities
            ));
        }
        report_string.push_str(" infeasibilities\n");
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kError,
            "{}",
            report_string
        );
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kError,
            "Setting model status to {}\n",
            util_model_status_to_string(self.model_status_)
        );
        HighsStatus::kError
    }

    /// HiGHS: Highs.cpp::Highs::returnFromOptimizeModel (lines 4784-4960),
    /// ending in `returnFromHighs` (highs_api_support.rs), which resizes an
    /// invalid basis to `kNonbasic` columns / `kBasic` rows
    /// (`forceHighsSolutionBasisSize`); the debug checks and reporting are out
    /// of scope.
    pub(crate) fn return_from_optimize_model(
        &mut self,
        run_return_status: HighsStatus,
        undo_mods: bool,
    ) -> HighsStatus {
        let mut return_status = highs_status_from_highs_model_status(self.model_status_);
        if return_status != run_return_status {
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Highs::returnFromOptimizeModel: run_return_status = {} != {} = return_status = highsStatusFromHighsModelStatus(model_status_ = {})\n",
                run_return_status as i32,
                return_status as i32,
                util_model_status_to_string(self.model_status_)
            );
        }
        match self.model_status_ {
            // First consider the error returns
            HighsModelStatus::kNotset
            | HighsModelStatus::kLoadError
            | HighsModelStatus::kModelError
            | HighsModelStatus::kPresolveError
            | HighsModelStatus::kSolveError
            | HighsModelStatus::kPostsolveError
            | HighsModelStatus::kMemoryLimit
            | HighsModelStatus::kModelEmpty => {
                // Don't clear the model status!
                self.info_.invalidate();
                self.run_data_.invalidate();
                self.invalidate_solution();
                self.invalidate_basis();
            }
            HighsModelStatus::kUnboundedOrInfeasible => {
                if self.options_.allow_unbounded_or_infeasible
                    || (use_ipm(&self.options_.solver)
                        && self.options_.run_crossover == kHighsOnString)
                    || (self.options_.solver == kPdlpString)
                    || self.model_lp_.is_mip()
                {
                    debug_assert!(return_status == HighsStatus::kOk);
                } else {
                    // This model status is not permitted unless IPM is run without
                    // crossover, or if PDLP is used
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "returnFromHighs: HighsModelStatus::kUnboundedOrInfeasible is not permitted\n"
                    );
                    debug_assert!(self.options_.allow_unbounded_or_infeasible);
                    return_status = HighsStatus::kError;
                }
            }
            _ => {}
        }
        // The debug checks and reportSolvedLpQpStats (logging) are out of scope.
        //
        // (lines 4937-4941)
        if undo_mods {
            // Restore any infinite costs
            self.restore_inf_cost(&mut return_status);
            // Unapply any modifications that have not yet been unapplied
            self.model_lp_.unapply_mods();
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.h::Highs::getModelStatus
    pub fn get_model_status(&self) -> HighsModelStatus {
        self.model_status_
    }

    /// HiGHS: Highs.h::Highs::getInfo
    pub fn get_info(&self) -> &HighsInfo {
        &self.info_
    }

    /// HiGHS: Highs.h::Highs::getRunData (line 430)
    pub fn get_run_data(&self) -> &HighsRunData {
        &self.run_data_
    }

    /// HiGHS: Highs.h::Highs::getSolution
    pub fn get_solution(&self) -> &HighsSolution {
        &self.solution_
    }

    /// HiGHS: Highs.h::Highs::getBasis
    pub fn get_basis(&self) -> &HighsBasis {
        &self.basis_
    }
}

/// rusolve (racing, Shinano et al. 2018; Carvajal et al. 2014): run `n` diversified copies of the MIP search on `n`
/// threads sharing incumbents and a stop flag, and return the solver whose result stands.
///
/// The winner is the first racer to reach a terminal status (ties: best objective, lowest
/// index). With no terminal racer (time limit), it is the racer with the best incumbent, and
/// the dual bound is the best any racer proved. Node and LP-iteration counts are summed.
pub(crate) fn race_solve_mip(
    callback: &mut HighsCallback,
    options: &HighsOptions,
    lp: &HighsLp,
    solution: &HighsSolution,
    n: usize,
) -> HighsMipSolver {
    use crate::mip::race::{RaceShared, racer_options, run_race};
    let shared = std::sync::Arc::new(RaceShared::default());
    let cb0 = std::sync::Mutex::new(Some(std::mem::take(callback)));
    let terminal = |m: HighsModelStatus| {
        !matches!(
            m,
            HighsModelStatus::kTimeLimit | HighsModelStatus::kInterrupt | HighsModelStatus::kNotset
        )
    };
    let mut racers = run_race(n, &shared, |k, shared| {
        let mut cb = if k == 0 {
            cb0.lock().unwrap().take().unwrap_or_default()
        } else {
            HighsCallback::default()
        };
        let opts = racer_options(options, k);
        let mut solver = HighsMipSolver::new(&mut cb, &opts, lp, solution, false, 0);
        solver.race = Some(std::sync::Arc::clone(shared));
        solver.set_profiling();
        solver.run();
        solver.callback_ = cb;
        let t = terminal(solver.modelstatus_);
        (solver, t)
    });
    let sense = lp.sense_ as HighsInt as f64; // 1 min, -1 max
    let nodes: i64 = racers.iter().map(|r| r.node_count_).sum();
    let iters: i64 = racers.iter().map(|r| r.total_lp_iterations_).sum();
    let best_dual = racers
        .iter()
        .map(|r| sense * r.dual_bound_)
        .fold(-kHighsInf, f64::max);
    let key = |r: &HighsMipSolver| sense * r.solution_objective_;
    let mut win = 0;
    for i in 1..racers.len() {
        let (t, w_t) = (
            terminal(racers[i].modelstatus_),
            terminal(racers[win].modelstatus_),
        );
        if (t && !w_t) || (t == w_t && key(&racers[i]) < key(&racers[win])) {
            win = i;
        }
    }
    let any_time_limit = racers
        .iter()
        .any(|r| r.modelstatus_ == HighsModelStatus::kTimeLimit);
    let mut solver = racers.swap_remove(win);
    // racer 0 owns the user callback: hand it back whichever racer won
    if win != 0 {
        solver.callback_ = std::mem::take(&mut racers[0].callback_);
    }
    solver.node_count_ = nodes;
    solver.total_lp_iterations_ = iters;
    if !terminal(solver.modelstatus_) {
        // limit stop: report the time limit (a racer stopped by the flag reads as an interrupt)
        if any_time_limit && solver.modelstatus_ == HighsModelStatus::kInterrupt {
            solver.modelstatus_ = HighsModelStatus::kTimeLimit;
        }
        // only the bound is shared across racers; each racer's bound is valid for the problem
        let d = sense * best_dual;
        solver.dual_bound_ = if sense > 0.0 {
            d.min(solver.primal_bound_)
        } else {
            d.max(solver.primal_bound_)
        };
        let (p, d) = (solver.primal_bound_, solver.dual_bound_);
        solver.gap_ = if p == 0.0 {
            if d == 0.0 { 0.0 } else { kHighsInf }
        } else if p.abs() != kHighsInf {
            (p - d).abs() / p.abs()
        } else {
            kHighsInf
        };
    }
    solver
}

#[cfg(test)]
mod h16_intassert_tests {
    use super::*;

    /// HiGHS 1.16 (#3576) check/TestMipSolver.cpp issue-3573 (3576.diff lines 5-45): two
    /// integer columns, presolve off, eps = 0.5 * mip_feasibility_tolerance. The three
    /// column-0 bounds are run in turn. Tests run in --release, so the debug asserts that
    /// differ between the flag values are compiled out; the statuses are the same either way.
    fn issue_3573_statuses(flag: bool) -> [HighsModelStatus; 3] {
        let eps = 0.5 * HighsOptions::default().mip_feasibility_tolerance;
        let lp_with = |lower0: f64, upper0: f64| HighsLp {
            num_col_: 2,
            num_row_: 1,
            col_cost_: vec![1.0, 1.0],
            col_lower_: vec![lower0, 0.0],
            col_upper_: vec![upper0, 3.0],
            integrality_: vec![HighsVarType::kInteger, HighsVarType::kInteger],
            row_lower_: vec![3.0],
            row_upper_: vec![kHighsInf],
            a_matrix_: HighsSparseMatrix {
                num_col_: 2,
                num_row_: 1,
                start_: vec![0, 1, 2],
                index_: vec![0, 0],
                value_: vec![1.0, 1.0],
                ..Default::default()
            },
            ..Default::default()
        };
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("h16_int_round_assert_tol", flag),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("presolve", "off"),
            HighsStatus::kOk
        );
        let cases = [
            (2.0 + eps, 3.0),
            (2.0, 3.0 - eps),
            (2.0 + 3.0 * eps, 3.0 - 3.0 * eps),
        ];
        let mut statuses = [HighsModelStatus::kNotset; 3];
        for (k, (lower0, upper0)) in cases.iter().enumerate() {
            assert_eq!(
                highs.pass_model(lp_with(*lower0, *upper0)),
                HighsStatus::kOk
            );
            assert_eq!(highs.run(), HighsStatus::kOk);
            statuses[k] = highs.get_model_status();
        }
        statuses
    }

    #[test]
    fn h16_intround_on_issue_3573_statuses() {
        assert_eq!(
            issue_3573_statuses(true),
            [
                HighsModelStatus::kOptimal,
                HighsModelStatus::kOptimal,
                HighsModelStatus::kInfeasible
            ]
        );
    }

    /// Flag off keeps the 1.15.1 asserts `integer_lower >= lower`, which fail for
    /// lower = 2 + eps (ceil(2 - eps) = 2): debug builds panic here, release builds
    /// (asserts compiled out) return the same statuses as the flag-on path.
    #[test]
    #[cfg_attr(debug_assertions, should_panic)]
    fn h16_intround_off_issue_3573_statuses() {
        assert_eq!(
            issue_3573_statuses(false),
            [
                HighsModelStatus::kOptimal,
                HighsModelStatus::kOptimal,
                HighsModelStatus::kInfeasible
            ]
        );
    }
}
