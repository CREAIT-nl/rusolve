//! HiGHS source: simplex/HEkk.cpp (lines 31-450)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Clearing/invalidation, status updates, pointer plumbing, LP move-in.

use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkk.cpp::HEkk::clear (lines 31-43, 13 lines)
    pub fn clear(&mut self) {
        self.clear_ekk_lp();
        self.clear_ekk_dualize();
        self.clear_ekk_data();
        self.clear_ekk_dual_edge_weight_data();
        self.clear_ekk_pointers();
        self.basis_.clear();
        self.simplex_nla_.clear();
        self.clear_ekk_all_status();
        self.clear_ray_records();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkAllStatus (lines 45-54, 10 lines)
    pub fn clear_ekk_all_status(&mut self) {
        self.status_.initialised_for_new_lp = false;
        self.status_.initialised_for_solve = false;
        self.clear_nla_status();
        self.clear_ekk_data_status();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkDataStatus (lines 56-66, 11 lines)
    pub fn clear_ekk_data_status(&mut self) {
        self.status_.has_ar_matrix = false;
        self.status_.has_dual_steepest_edge_weights = false;
        self.status_.has_fresh_rebuild = false;
        self.status_.has_dual_objective_value = false;
        self.status_.has_primal_objective_value = false;
    }

    /// HiGHS: HEkk.cpp::HEkk::clearNlaStatus (lines 68-75, 8 lines)
    pub fn clear_nla_status(&mut self) {
        self.status_.has_basis = false;
        self.status_.has_nla = false;
        self.clear_nla_invert_status();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearNlaInvertStatus (lines 77-80, 4 lines)
    pub fn clear_nla_invert_status(&mut self) {
        self.status_.has_invert = false;
        self.status_.has_fresh_invert = false;
    }

    /// HiGHS: HEkk.cpp::HEkk::clearRayRecords (lines 82-85, 4 lines)
    pub fn clear_ray_records(&mut self) {
        self.dual_ray_record_.clear();
        self.primal_ray_record_.clear();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkPointers (lines 87-91, 5 lines)
    pub fn clear_ekk_pointers(&mut self) {
        // Pattern (d): C++ nulls the pointers. The owned stand-ins are reset to their defaults;
        // `move_lp` / `set_pointers` set them again before any use.
        self.callback_ = HighsCallback::default();
        self.options_ = HighsOptions::default();
        self.timer_ = HighsTimer::default();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkLp (lines 93-96, 4 lines)
    pub fn clear_ekk_lp(&mut self) {
        self.lp_.clear();
        self.lp_name_ = String::new();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkDualize (lines 98-106, 9 lines)
    pub fn clear_ekk_dualize(&mut self) {
        self.original_col_cost_.clear();
        self.original_col_lower_.clear();
        self.original_col_upper_.clear();
        self.original_row_lower_.clear();
        self.original_row_upper_.clear();
        self.upper_bound_col_.clear();
        self.upper_bound_row_.clear();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkDualEdgeWeightData (lines 108-111, 4 lines)
    pub fn clear_ekk_dual_edge_weight_data(&mut self) {
        self.dual_edge_weight_.clear();
        self.scattered_dual_edge_weight_.clear();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkData (lines 113-166, 54 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion): resets `build_synthetic_tick_`/`total_synthetic_tick_`.
    pub fn clear_ekk_data(&mut self) {
        self.clear_ekk_data_info();
        self.model_status_ = HighsModelStatus::kNotset;
        // random_; Has no data

        self.simplex_in_scaled_space_ = false;
        self.ar_matrix_.clear();
        self.scaled_a_matrix_.clear();

        self.cost_scale_ = 1.0;
        self.iteration_count_ = 0;
        self.dual_simplex_cleanup_level_ = 0;
        self.dual_simplex_phase1_cleanup_level_ = 0;

        self.previous_iteration_cycling_detected = -kHighsIInf;

        self.solve_bailout_ = false;
        self.called_return_from_solve_ = false;
        self.exit_algorithm_ = SimplexAlgorithm::kNone;
        self.return_primal_solution_status_ = 0;
        self.return_dual_solution_status_ = 0;

        self.proof_index_.clear();
        self.proof_value_.clear();

        self.clear_ray_records();

        self.build_synthetic_tick_ = 0.0;
        self.total_synthetic_tick_ = 0.0;

        // Clear values used for debugging
        self.debug_solve_call_num_ = 0;
        self.debug_basis_id_ = 0;
        self.time_report_ = false;
        self.debug_initial_build_synthetic_tick_ = 0;
        self.debug_solve_report_ = false;
        self.debug_iteration_report_ = false;
        self.debug_basis_report_ = false;
        self.debug_dual_feasible = false;
        self.debug_max_relative_dual_steepest_edge_weight_error = 0.0;

        self.clear_bad_basis_change(BadBasisChangeReason::kAll);

        self.primal_phase1_dual_.clear();
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkDataInfo (lines 168-249, 82 lines)
    /// MUST-PORT (determinism, random stream): clears the random vectors (`num_tot_random_value_` etc.).
    pub fn clear_ekk_data_info(&mut self) {
        let info = &mut self.info_;
        info.work_cost_.clear();
        info.work_dual_.clear();
        info.work_shift_.clear();
        info.work_lower_.clear();
        info.work_upper_.clear();
        info.work_range_.clear();
        info.work_value_.clear();
        info.work_lower_shift_.clear();
        info.work_upper_shift_.clear();
        info.base_lower_.clear();
        info.base_upper_.clear();
        info.base_value_.clear();
        info.num_tot_random_value_.clear();
        info.num_tot_permutation_.clear();
        info.num_col_permutation_.clear();
        info.devex_index_.clear();
        info.pivot_.clear();
        info.index_chosen_.clear();
        info.phase1_backtracking_test_done = 0;
        info.phase2_backtracking_test_done = 0;
        info.backtracking_ = false;
        info.valid_backtracking_basis_ = false;
        info.backtracking_basis_.clear();
        info.backtracking_basis_costs_shifted_ = 0;
        info.backtracking_basis_costs_perturbed_ = 0;
        info.backtracking_basis_bounds_shifted_ = 0;
        info.backtracking_basis_bounds_perturbed_ = 0;
        info.backtracking_basis_work_shift_.clear();
        info.backtracking_basis_work_lower_shift_.clear();
        info.backtracking_basis_work_upper_shift_.clear();
        info.backtracking_basis_edge_weight_.clear();
        info.simplex_strategy = 0;
        info.dual_edge_weight_strategy = 0;
        info.primal_edge_weight_strategy = 0;
        info.price_strategy = 0;
        info.dual_simplex_cost_perturbation_multiplier = 1.0;
        info.primal_simplex_phase1_cost_perturbation_multiplier = 1.0;
        info.primal_simplex_bound_perturbation_multiplier = 1.0;
        info.allow_dual_steepest_edge_to_devex_switch = false;
        info.dual_steepest_edge_weight_log_error_threshold = 0.0;
        info.run_quiet = false;
        info.store_squared_primal_infeasibility = false;
        info.report_simplex_inner_clock = false;
        info.report_simplex_outer_clock = false;
        info.report_simplex_phases_clock = false;
        info.report_h_factor_clock = false;
        info.analyse_lp = false;
        info.analyse_iterations = false;
        info.analyse_invert_form = false;

        info.allow_cost_shifting = true;
        info.allow_cost_perturbation = true;
        info.allow_bound_perturbation = true;
        info.costs_shifted = false;
        info.costs_perturbed = false;
        info.bounds_shifted = false;
        info.bounds_perturbed = false;

        info.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
        info.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
        info.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
        info.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
        info.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
        info.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
        info.dual_phase1_iteration_count = 0;
        info.dual_phase2_iteration_count = 0;
        info.primal_phase1_iteration_count = 0;
        info.primal_phase2_iteration_count = 0;
        info.primal_bound_swap = 0;
        info.min_concurrency = 1;
        info.num_concurrency = 1;
        info.max_concurrency = kSimplexConcurrencyLimit;
        info.multi_iteration = 0;
        info.update_count = 0;
        info.dual_objective_value = 0.0;
        info.primal_objective_value = 0.0;
        info.updated_dual_objective_value = 0.0;
        info.updated_primal_objective_value = 0.0;
        info.num_basic_logicals = 0;
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkControlInfo (lines 251-269, 19 lines)
    pub fn clear_ekk_control_info(&mut self) {
        let info = &mut self.info_;
        info.control_iteration_count0 = 0;
        info.col_aq_density = 0.0;
        info.row_ep_density = 0.0;
        info.row_ap_density = 0.0;
        info.row_dse_density = 0.0;
        info.col_steepest_edge_density = 0.0;
        info.col_basic_feasibility_change_density = 0.0;
        info.row_basic_feasibility_change_density = 0.0;
        info.col_bfrt_density = 0.0;
        info.primal_col_density = 0.0;
        info.dual_col_density = 0.0;
        info.costly_dse_frequency = 0.0;
        info.num_costly_dse_iteration = 0;
        info.costly_dse_measure = 0.0;
        info.average_log_low_dse_weight_error = 0.0;
        info.average_log_high_dse_weight_error = 0.0;
    }

    /// HiGHS: HEkk.cpp::HEkk::clearEkkNlaInfo (lines 271-275, 5 lines)
    pub fn clear_ekk_nla_info(&mut self) {
        self.info_.factor_pivot_threshold = 0.0;
        self.info_.update_limit = 0;
    }

    /// HiGHS: HEkk.cpp::HEkk::invalidate (lines 277-284, 8 lines)
    pub fn invalidate(&mut self) {
        self.status_.initialised_for_new_lp = false;
        debug_assert!(!self.status_.is_dualized);
        debug_assert!(!self.status_.is_permuted);
        self.status_.initialised_for_solve = false;
        self.invalidate_basis_matrix();
        self.simplex_stats_.initialise(0);
    }

    /// HiGHS: HEkk.cpp::HEkk::invalidateBasisMatrix (lines 286-292, 7 lines)
    pub fn invalidate_basis_matrix(&mut self) {
        self.status_.has_nla = false;
        self.invalidate_basis();
    }

    /// HiGHS: HEkk.cpp::HEkk::invalidateBasis (lines 294-299, 6 lines)
    pub fn invalidate_basis(&mut self) {
        self.status_.has_basis = false;
        self.invalidate_basis_artifacts();
    }

    /// HiGHS: HEkk.cpp::HEkk::invalidateBasisArtifacts (lines 301-311, 11 lines)
    pub fn invalidate_basis_artifacts(&mut self) {
        self.status_.has_ar_matrix = false;
        self.status_.has_dual_steepest_edge_weights = false;
        self.status_.has_invert = false;
        self.status_.has_fresh_invert = false;
        self.status_.has_fresh_rebuild = false;
        self.status_.has_dual_objective_value = false;
        self.status_.has_primal_objective_value = false;
        self.clear_ray_records();
    }

    /// HiGHS: HEkk.cpp::HEkk::updateStatus (lines 313-378, 66 lines)
    pub fn update_status(&mut self, action: LpAction) {
        debug_assert!(!self.status_.is_dualized);
        debug_assert!(!self.status_.is_permuted);
        match action {
            LpAction::kScale => {
                self.invalidate_basis_matrix();
            }
            LpAction::kNewCosts => {
                self.status_.has_fresh_rebuild = false;
                self.status_.has_dual_objective_value = false;
                self.status_.has_primal_objective_value = false;
            }
            LpAction::kNewBounds => {
                self.status_.has_fresh_rebuild = false;
                self.status_.has_dual_objective_value = false;
                self.status_.has_primal_objective_value = false;
            }
            LpAction::kNewBasis => {
                self.invalidate_basis();
            }
            LpAction::kNewCols => {
                self.clear();
                //    self.invalidate_basis_artifacts();
            }
            LpAction::kNewRows => {
                if kExtendInvertWhenAddingRows {
                    // Just clear Ekk data
                    self.clear_ekk_data();
                } else {
                    // Clear everything
                    self.clear();
                }
                //    self.invalidate_basis_artifacts();
            }
            LpAction::kDelCols => {
                self.clear();
                //    self.invalidate_basis();
            }
            LpAction::kDelNonbasicCols => {
                self.clear();
                //    self.invalidate_basis();
            }
            LpAction::kDelRows => {
                self.clear();
                //   self.invalidate_basis();
            }
            LpAction::kDelRowsBasisOk => {
                // C++ `assert(1 == 0)`: unreachable, debug builds only.
                if cfg!(debug_assertions) {
                    panic!("kDelRowsBasisOk is not a valid LpAction in updateStatus");
                }
                //      info.lp_ = true;
            }
            LpAction::kScaledCol => {
                self.invalidate_basis_matrix();
            }
            LpAction::kScaledRow => {
                self.invalidate_basis_matrix();
            }
            LpAction::kBacktracking => {
                self.status_.has_ar_matrix = false;
                self.status_.has_fresh_rebuild = false;
                self.status_.has_dual_objective_value = false;
                self.status_.has_primal_objective_value = false;
            }
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::setNlaPointersForLpAndScale (lines 380-383, 4 lines)
    /// No-op (pattern (b)): HSimplexNla receives `NlaRefs` per call. The C++ `assert(status_.has_nla)`
    /// may be kept as a `debug_assert!`.
    pub fn set_nla_pointers_for_lp_and_scale(&mut self, lp: &HighsLp) {}

    /// HiGHS: HEkk.cpp::HEkk::setNlaPointersForTrans (lines 385-390, 6 lines)
    /// No-op (pattern (b)); see `set_nla_pointers_for_lp_and_scale`.
    pub fn set_nla_pointers_for_trans(&mut self, lp: &HighsLp) {}

    /// HiGHS: HEkk.cpp::HEkk::setNlaRefactorInfo (lines 392-395, 4 lines)
    pub fn set_nla_refactor_info(&mut self) {
        self.simplex_nla_.factor_.refactor_info_ = self.hot_start_.refactor_info.clone();
        self.simplex_nla_.factor_.refactor_info_.r#use = true;
    }

    /// HiGHS: HEkk.cpp::HEkk::btran (lines 397-400, 4 lines)
    /// Caveat: in C++ `simplex_nla_.lp_` is whatever LP was last set by
    /// `setNlaPointersForLpAndScale`. Throughout `HEkk::solve` that is `HEkk::lp_`; only callers in
    /// the Highs interface layer (out of scope) use it after the LP has been moved back to the
    /// incumbent. Such a caller must build `NlaRefs` from the incumbent LP instead.
    /// `simplex_nla_.btran(rhs, d)` -> `let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
    /// self.simplex_nla_.btran(&refs, rhs, d, None)`.
    pub fn btran(&mut self, rhs: &mut HVector, expected_density: f64) {
        debug_assert!(self.status_.has_nla);
        let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
        self.simplex_nla_.btran(&refs, rhs, expected_density, None);
    }

    /// HiGHS: HEkk.cpp::HEkk::ftran (lines 402-405, 4 lines)
    /// As `btran`, with `ftran`.
    pub fn ftran(&mut self, rhs: &mut HVector, expected_density: f64) {
        debug_assert!(self.status_.has_nla);
        let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
        self.simplex_nla_.ftran(&refs, rhs, expected_density, None);
    }

    /// HiGHS: HEkk.cpp::HEkk::moveLp (lines 407-429, 23 lines)
    /// `this->lp_ = std::move(incumbent_lp)` -> `self.lp_ = std::mem::take(&mut *solver_object.lp_)`;
    /// `setPointers(&solver_object.callback_, ...)` passes the referents (`&*solver_object.callback_`, ...).
    pub fn move_lp(&mut self, solver_object: &mut HighsLpSolverObject) {
        // Move the incumbent LP to EKK
        self.lp_ = std::mem::take(&mut *solver_object.lp_);
        solver_object.lp_.is_moved_ = true;
        //
        // Invalidate the row-wise matrix
        self.status_.has_ar_matrix = false;
        //
        // The simplex algorithm runs in the same space as the LP that has
        // just been moved in. This is a scaled space if the LP is scaled.
        self.simplex_in_scaled_space_ = self.lp_.is_scaled_;
        //
        // Update other EKK pointers. Currently just pointers to the
        // HighsOptions and HighsTimer members of the Highs class that are
        // communicated by reference via the HighsLpSolverObject instance.
        self.set_pointers(
            &*solver_object.callback_,
            &*solver_object.options_,
            &*solver_object.timer_,
        );
        // Initialise Ekk if this has not been done. Ekk isn't initialised
        // if moveLp hasn't been called for this instance of HiGHS, or if
        // the Ekk instance is junked due to removing rows from the LP
        self.initialise_ekk();
    }

    /// HiGHS: HEkk.cpp::HEkk::setPointers (lines 431-437, 7 lines)
    /// Pattern (d): HEkk owns these. Translate as `self.callback_ = callback.clone(); self.options_ =
    /// options.clone(); self.timer_ = timer.clone();` (`analysis_.timer_` has no counterpart).
    pub fn set_pointers(
        &mut self,
        callback: &HighsCallback,
        options: &HighsOptions,
        timer: &HighsTimer,
    ) {
        self.callback_ = callback.clone();
        self.options_ = options.clone();
        self.timer_ = timer.clone();
    }

    /// HiGHS: HEkk.cpp::HEkk::getScaledAMatrixPointer (lines 439-450, 12 lines)
    /// Call site (`initialiseSimplexLpBasisAndFactor`): the returned borrow of `self` cannot be
    /// held across `self.simplex_nla_.set_pointers(..)`. Call this for its side effect (refreshing
    /// `scaled_a_matrix_`), drop the result, then pass the same matrix by disjoint field borrow:
    /// `if self.lp_.scale_.has_scaling && !self.lp_.is_scaled_ { &self.scaled_a_matrix_ } else { &self.lp_.a_matrix_ }`.
    pub fn get_scaled_a_matrix_pointer(&mut self) -> &mut HighsSparseMatrix {
        if self.lp_.scale_.has_scaling && !self.lp_.is_scaled_ {
            self.scaled_a_matrix_ = self.lp_.a_matrix_.clone();
            self.scaled_a_matrix_.apply_scale(&self.lp_.scale_);
            &mut self.scaled_a_matrix_
        } else {
            &mut self.lp_.a_matrix_
        }
    }
}
