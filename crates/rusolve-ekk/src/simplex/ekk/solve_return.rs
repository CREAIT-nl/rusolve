//! HiGHS source: simplex/HEkk.cpp (lines 3444-4052)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Bailout/return-from-solve, basis condition, iterate freeze, solve error, bad basis
//! change bookkeeping and taboo application.

use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkk.cpp::HEkk::bailout (lines 3444-3472, 29 lines)
    pub fn bailout(&mut self) -> bool {
        if self.solve_bailout_ {
            // Bailout has already been decided: check that it's for one of these
            // reasons
            debug_assert!(
                self.model_status_ == HighsModelStatus::kTimeLimit
                    || self.model_status_ == HighsModelStatus::kIterationLimit
                    || self.model_status_ == HighsModelStatus::kObjectiveBound
                    || self.model_status_ == HighsModelStatus::kObjectiveTarget
            );
        } else if self.options_.time_limit < kHighsInf
            && self.timer_.read() > self.options_.time_limit
        {
            self.solve_bailout_ = true;
            self.model_status_ = HighsModelStatus::kTimeLimit;
        } else if self.iteration_count_ >= self.options_.simplex_iteration_limit {
            self.solve_bailout_ = true;
            self.model_status_ = HighsModelStatus::kIterationLimit;
        } else if self.callback_.user_callback.is_some()
            && self.callback_.active[kCallbackSimplexInterrupt as usize]
        {
            self.callback_.clear_highs_callback_output();
            self.callback_.data_out.simplex_iteration_count = self.iteration_count_;
            if self
                .callback_
                .callback_action(kCallbackSimplexInterrupt, "Simplex interrupt")
            {
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "User interrupt\n"
                );
                self.solve_bailout_ = true;
                self.model_status_ = HighsModelStatus::kInterrupt;
            }
        }
        self.solve_bailout_
    }

    /// HiGHS: HEkk.cpp::HEkk::returnFromEkkSolve (lines 3474-3498, 25 lines)
    pub fn return_from_ekk_solve(&mut self, return_status: HighsStatus) -> HighsStatus {
        if self.analysis_.analyse_simplex_time {
            self.analysis_.simplex_timer_stop(SimplexTotalClock, 0);
        }
        // Restore any modified development or timing settings and analyse
        // solver timing
        if self.debug_solve_report_ {
            self.debug_reporting(1, kHighsLogDevLevelDetailed);
        }
        if self.time_report_ {
            self.time_reporting(1);
        }
        // Note that in timeReporting(1), analysis_.analyse_simplex_time
        // reverts to its value given by options_
        if self.analysis_.analyse_simplex_time {
            self.analysis_.report_simplex_timer();
        }
        self.simplex_stats_.valid = true;
        // Since HEkk::iteration_count_ includes iteration on presolved LP,
        // simplex_stats_.iteration_count is initialised to -
        // HEkk::iteration_count_
        self.simplex_stats_.iteration_count += self.iteration_count_;
        //  simplex_stats_.num_invert is incremented internally
        self.simplex_stats_.last_invert_num_el = self.simplex_nla_.factor_.invert_num_el;
        self.simplex_stats_.last_factored_basis_num_el =
            self.simplex_nla_.factor_.basis_matrix_num_el;
        self.simplex_stats_.col_aq_density = self.analysis_.col_aq_density;
        self.simplex_stats_.row_ep_density = self.analysis_.row_ep_density;
        self.simplex_stats_.row_ap_density = self.analysis_.row_ap_density;
        self.simplex_stats_.row_dse_density = self.analysis_.row_DSE_density;
        return_status
    }

    /// HiGHS: HEkk.cpp::HEkk::returnFromSolve (lines 3500-3647, 148 lines)
    pub fn return_from_solve(&mut self, return_status: HighsStatus) -> HighsStatus {
        // Always called before returning from HEkkPrimal/Dual::solve()
        if self.solve_bailout_ {
            // If bailout has already been decided: check that it's for one of
            // these reasons
            debug_assert!(
                self.model_status_ == HighsModelStatus::kTimeLimit
                    || self.model_status_ == HighsModelStatus::kIterationLimit
                    || self.model_status_ == HighsModelStatus::kObjectiveBound
                    || self.model_status_ == HighsModelStatus::kObjectiveTarget
                    || self.model_status_ == HighsModelStatus::kInterrupt
            );
        }
        // Check that returnFromSolve has not already been called: it should
        // be called exactly once per solve
        debug_assert!(!self.called_return_from_solve_);
        self.called_return_from_solve_ = true;
        self.info_.valid_backtracking_basis_ = false;

        // Initialise the status of the primal and dual solutions
        self.return_primal_solution_status_ = kSolutionStatusNone;
        self.return_dual_solution_status_ = kSolutionStatusNone;
        // Nothing more is known about the solve after an error return
        if return_status == HighsStatus::kError {
            return return_status;
        }

        // Check that an invert exists
        debug_assert!(self.status_.has_invert);

        // Determine a primal and dual solution, removing the effects of
        // perturbations and shifts
        //
        // Unless the solution is optimal, invalidate the infeasibility data
        if self.model_status_ != HighsModelStatus::kOptimal {
            self.invalidate_primal_infeasibility_record();
            self.invalidate_dual_infeasibility_record();
        }
        // The simplex algorithm used on exit should be set
        debug_assert!(self.exit_algorithm_ != SimplexAlgorithm::kNone);
        match self.model_status_ {
            HighsModelStatus::kOptimal => {
                if self.info_.num_primal_infeasibilities != 0 {
                    // Optimal - but not to desired primal feasibility tolerance
                    self.return_primal_solution_status_ = kSolutionStatusInfeasible;
                } else {
                    self.return_primal_solution_status_ = kSolutionStatusFeasible;
                }
                if self.info_.num_dual_infeasibilities != 0 {
                    // Optimal - but not to desired dual feasibility tolerance
                    self.return_dual_solution_status_ = kSolutionStatusInfeasible;
                } else {
                    self.return_dual_solution_status_ = kSolutionStatusFeasible;
                }
            }
            HighsModelStatus::kInfeasible => {
                // Primal infeasibility has been identified in primal phase 1,
                // or proved in dual phase 2. There should be no primal
                // perturbations or shifts
                debug_assert!(!self.info_.bounds_shifted);
                debug_assert!(!self.info_.bounds_perturbed);
                if self.exit_algorithm_ == SimplexAlgorithm::kPrimal {
                    // Reset the simplex costs and recompute duals after primal
                    // phase 1
                    self.initialise_cost(SimplexAlgorithm::kDual, kSolvePhase2, false);
                    self.compute_dual();
                }
                self.compute_simplex_infeasible();
                // Primal solution shouldn't be feasible
                debug_assert!(self.info_.num_primal_infeasibilities > 0);
            }
            HighsModelStatus::kUnboundedOrInfeasible => {
                // Dual simplex has identified dual infeasibility in phase
                // 1. There should be no dual perturbations
                debug_assert!(self.exit_algorithm_ == SimplexAlgorithm::kDual);
                debug_assert!(!self.info_.costs_shifted);
                debug_assert!(!self.info_.costs_perturbed);
                // Reset the simplex bounds and recompute primals
                self.initialise_bound(SimplexAlgorithm::kDual, kSolvePhase2, false);
                self.compute_primal();
                self.compute_simplex_infeasible();
                // Dual solution shouldn't be feasible
                debug_assert!(self.info_.num_dual_infeasibilities > 0);
            }
            HighsModelStatus::kUnbounded => {
                // Primal simplex has identified unboundedness in phase 2. There
                // should be no primal or dual perturbations
                debug_assert!(self.exit_algorithm_ == SimplexAlgorithm::kPrimal);
                debug_assert!(!self.info_.costs_shifted);
                debug_assert!(!self.info_.costs_perturbed);
                debug_assert!(!self.info_.bounds_shifted);
                debug_assert!(!self.info_.bounds_perturbed);
                self.compute_simplex_infeasible();
                // Primal solution should be feasible
                debug_assert!(self.info_.num_primal_infeasibilities == 0);
            }
            HighsModelStatus::kObjectiveBound
            | HighsModelStatus::kObjectiveTarget
            | HighsModelStatus::kTimeLimit
            | HighsModelStatus::kIterationLimit
            | HighsModelStatus::kInterrupt
            | HighsModelStatus::kUnknown => {
                // Simplex has failed to conclude a model property. Either it
                // has bailed out due to reaching the objective bound, target,
                // time, iteration limit or user interrupt, or it has not been
                // set (cycling is the only reason). Could happen anywhere.
                //
                // Reset the simplex bounds and recompute primals
                self.initialise_bound(SimplexAlgorithm::kDual, kSolvePhase2, false);
                self.initialise_nonbasic_value_and_move();
                self.compute_primal();
                // Reset the simplex costs and recompute duals
                self.initialise_cost(SimplexAlgorithm::kDual, kSolvePhase2, false);
                self.compute_dual();
                self.compute_simplex_infeasible();
            }
            _ => {
                let algorithm_name = if self.exit_algorithm_ == SimplexAlgorithm::kPrimal {
                    "primal"
                } else {
                    "dual"
                };
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "{} simplex solver returns status {}\n",
                    algorithm_name,
                    util_model_status_to_string(self.model_status_)
                );
                return HighsStatus::kError;
            }
        }
        self.zero_basic_duals();
        debug_assert!(self.info_.num_primal_infeasibilities >= 0);
        debug_assert!(self.info_.num_dual_infeasibilities >= 0);
        if self.info_.num_primal_infeasibilities == 0 {
            self.return_primal_solution_status_ = kSolutionStatusFeasible;
        } else {
            self.return_primal_solution_status_ = kSolutionStatusInfeasible;
        }
        if self.info_.num_dual_infeasibilities == 0 {
            self.return_dual_solution_status_ = kSolutionStatusFeasible;
        } else {
            self.return_dual_solution_status_ = kSolutionStatusInfeasible;
        }
        debug_assert!(self.debug_no_shifts_or_perturbations());
        self.compute_primal_objective_value();
        if self.options_.log_dev_level == 0 {
            let force = true;
            self.analysis_.user_invert_report(force);
        }
        return_status
    }

    /// HiGHS: HEkk.cpp::HEkk::computeBasisCondition (lines 3649-3774, 126 lines)
    /// C++ default arguments: exact = false; report = false
    /// The header-inline no-argument overload is `compute_basis_condition` (in `ekk/mod.rs`).
    pub fn compute_basis_condition_lp(&self, lp: &HighsLp, exact: bool, report: bool) -> f64 {
        let solver_num_row: HighsInt = lp.num_row_;
        let solver_num_col: HighsInt = lp.num_col_;
        let mut bs_cond_x: Vec<f64> = Vec::new();
        let mut bs_cond_y: Vec<f64> = Vec::new();
        let mut bs_cond_z: Vec<f64> = Vec::new();
        let mut bs_cond_w: Vec<f64> = Vec::new();
        let mut row_ep = HVector::default();
        row_ep.setup(solver_num_row);

        let astart: &[HighsInt] = &lp.a_matrix_.start_;
        let avalue: &[f64] = &lp.a_matrix_.value_;
        // C++ simplex_nla_.ftran/btran read HEkk's own lp_ and basis through
        // the HSimplexNla pointers (setLpAndScalePointers), so the refs are
        // built from self, not from the `lp` argument.
        let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
        let mut exact_norm_binv: f64 = 0.0;
        if exact {
            // Compute the exact norm of B^{-1}
            for r_n in 0..solver_num_row {
                row_ep.clear();
                row_ep.index[row_ep.count as usize] = r_n;
                row_ep.array[r_n as usize] = 1.0;
                row_ep.count += 1;
                row_ep.pack_flag = false;
                self.simplex_nla_.ftran(&refs, &mut row_ep, 0.1, None);
                debug_assert!(row_ep.count <= solver_num_row);
                let mut c_norm: f64 = 0.0;
                for i_x in 0..row_ep.count {
                    c_norm += row_ep.array[row_ep.index[i_x as usize] as usize].abs();
                }
                exact_norm_binv = cpp_max(c_norm, exact_norm_binv);
            }
        }
        // Compute the Hager condition number estimate for the basis matrix
        let expected_density: f64 = 1.0;
        bs_cond_x.resize(solver_num_row as usize, 0.0);
        bs_cond_y.resize(solver_num_row as usize, 0.0);
        bs_cond_z.resize(solver_num_row as usize, 0.0);
        bs_cond_w.resize(solver_num_row as usize, 0.0);
        // x = ones(n,1)/n;
        // y = A\x;
        let mu: f64 = 1.0 / (solver_num_row as f64);
        let mut norm_binv: f64 = 0.0;
        for r_n in 0..solver_num_row {
            bs_cond_x[r_n as usize] = mu;
        }
        row_ep.clear();
        for r_n in 0..solver_num_row {
            let value = bs_cond_x[r_n as usize];
            if value != 0.0 {
                row_ep.index[row_ep.count as usize] = r_n;
                row_ep.array[r_n as usize] = value;
                row_ep.count += 1;
            }
        }
        for _ps_n in 1..=5 {
            row_ep.pack_flag = false;
            self.simplex_nla_
                .ftran(&refs, &mut row_ep, expected_density, None);

            // zeta = sign(y);
            for r_n in 0..solver_num_row {
                bs_cond_y[r_n as usize] = row_ep.array[r_n as usize];
                if bs_cond_y[r_n as usize] > 0.0 {
                    bs_cond_w[r_n as usize] = 1.0;
                } else if bs_cond_y[r_n as usize] < 0.0 {
                    bs_cond_w[r_n as usize] = -1.0;
                } else {
                    bs_cond_w[r_n as usize] = 0.0;
                }
            }
            // z=A'\zeta;
            row_ep.clear();
            for r_n in 0..solver_num_row {
                let value = bs_cond_w[r_n as usize];
                if value != 0.0 {
                    row_ep.index[row_ep.count as usize] = r_n;
                    row_ep.array[r_n as usize] = value;
                    row_ep.count += 1;
                }
            }
            row_ep.pack_flag = false;
            self.simplex_nla_
                .btran(&refs, &mut row_ep, expected_density, None);
            let mut norm_z: f64 = 0.0;
            let mut ztx: f64 = 0.0;
            norm_binv = 0.0;
            let mut argmax_z: HighsInt = -1;
            for r_n in 0..solver_num_row {
                bs_cond_z[r_n as usize] = row_ep.array[r_n as usize];
                let abs_z_v = bs_cond_z[r_n as usize].abs();
                if abs_z_v > norm_z {
                    norm_z = abs_z_v;
                    argmax_z = r_n;
                }
                ztx += bs_cond_z[r_n as usize] * bs_cond_x[r_n as usize];
                norm_binv += bs_cond_y[r_n as usize].abs();
            }
            if norm_z <= ztx {
                break;
            }
            // x = zeros(n,1);
            // x(fd_i) = 1;
            for r_n in 0..solver_num_row {
                bs_cond_x[r_n as usize] = 0.0;
            }
            row_ep.clear();
            row_ep.count = 1;
            row_ep.index[0] = argmax_z;
            row_ep.array[argmax_z as usize] = 1.0;
            bs_cond_x[argmax_z as usize] = 1.0;
        }
        let mut norm_b: f64 = 0.0;
        for r_n in 0..solver_num_row {
            let vr_n = self.basis_.basic_index_[r_n as usize];
            let mut c_norm: f64 = 0.0;
            if vr_n < solver_num_col {
                for el_n in astart[vr_n as usize]..astart[(vr_n + 1) as usize] {
                    c_norm += avalue[el_n as usize].abs();
                }
            } else {
                c_norm += 1.0;
            }
            norm_b = cpp_max(c_norm, norm_b);
        }
        let cond_b = norm_binv * norm_b;
        let exact_cond_b = exact_norm_binv * norm_b;
        if exact {
            debug_assert!(exact_norm_binv > 0.0);
            if report {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "HEkk::computeBasisCondition: grep_kappa model,||B||_1,approx ||B^{{-1}}||_1,approx_kappa,||B^{{-1}}||_1,kappa = ,{},{},{},{},{},{}\n",
                    lp.model_name_,
                    norm_b,
                    norm_binv,
                    cond_b,
                    exact_norm_binv,
                    exact_cond_b
                );
            }
            return exact_cond_b;
        }
        cond_b
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseAnalysis (lines 3776-3778, 3 lines)
    pub fn initialise_analysis(&mut self) {
        self.analysis_.setup(
            &self.lp_name_,
            &self.lp_,
            &self.options_,
            self.iteration_count_,
        );
    }

    /// HiGHS: HEkk.cpp::HEkk::rebuildReason (lines 3780-3809, 30 lines)
    /// HiGHS 1.16 (#3351): HEkk.cpp::HEkk::rebuildReason (lines 3788-3819)
    pub fn rebuild_reason(&self, rebuild_reason: HighsInt) -> String {
        // HiGHS 1.16 #3351 (h16_rebuild_reason): read the flag once, before the chain
        let h16_rebuild_reason: bool = self.options_.h16_rebuild_reason;
        let rebuild_reason_string: &str = if rebuild_reason == kRebuildReasonCleanup {
            "Perturbation cleanup"
        } else if rebuild_reason == kRebuildReasonNo {
            "No reason"
        } else if rebuild_reason == kRebuildReasonUpdateLimitReached {
            "Update limit reached"
        } else if rebuild_reason == kRebuildReasonSyntheticClockSaysInvert {
            "Synthetic clock"
        } else if rebuild_reason == kRebuildReasonPossiblyOptimal {
            "Possibly optimal"
        } else if rebuild_reason == kRebuildReasonPossiblyPhase1Feasible {
            "Possibly phase 1 feasible"
        } else if rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded {
            "Possibly primal unbounded"
        } else if rebuild_reason == kRebuildReasonPossiblyDualUnbounded {
            "Possibly dual unbounded"
        } else if rebuild_reason == kRebuildReasonPossiblySingularBasis {
            "Possibly singular basis"
        } else if rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex {
            "Primal infeasible in primal simplex"
        } else if rebuild_reason == kRebuildReasonChooseColumnFail {
            "Choose column failure"
        } else if h16_rebuild_reason && rebuild_reason == kRebuildReasonExcessivePrimalValue {
            // HiGHS 1.16 #3351 (h16_rebuild_reason): HEkk.cpp lines 3812-3813 (1.16.0)
            "Excessive primal value"
        } else {
            #[allow(clippy::assertions_on_constants)]
            {
                debug_assert!(1 == 0);
            }
            "Unidentified"
        };
        rebuild_reason_string.to_string()
    }

    /// HiGHS: HEkk.cpp::HEkk::putIterate (lines 3811-3823, 13 lines)
    pub fn put_iterate(&mut self) {
        debug_assert!(self.status_.has_invert);
        self.simplex_nla_.put_invert();
        let iterate = &mut self.simplex_nla_.simplex_iterate_;
        iterate.basis_ = self.basis_.clone();
        if self.status_.has_dual_steepest_edge_weights {
            // Copy the dual edge weights
            iterate.dual_edge_weight_ = self.dual_edge_weight_.clone();
        } else {
            // Clear to indicate no weights
            iterate.dual_edge_weight_.clear();
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::getIterate (lines 3825-3837, 13 lines)
    pub fn get_iterate(&mut self) -> HighsStatus {
        if !self.simplex_nla_.simplex_iterate_.valid_ {
            return HighsStatus::kError;
        }
        self.simplex_nla_.get_invert();
        self.basis_ = self.simplex_nla_.simplex_iterate_.basis_.clone();
        if !self
            .simplex_nla_
            .simplex_iterate_
            .dual_edge_weight_
            .is_empty()
        {
            self.dual_edge_weight_ = self.simplex_nla_.simplex_iterate_.dual_edge_weight_.clone();
        } else {
            self.status_.has_dual_steepest_edge_weights = false;
        }
        self.status_.has_invert = true;
        HighsStatus::kOk
    }

    /// HiGHS: HEkk.cpp::HEkk::factorSolveError (lines 3839-3919, 81 lines)
    pub fn factor_solve_error(&mut self) -> f64 {
        // Cheap assessment of factor accuracy.
        //
        // Forms a random solution with at most 50 nonzeros, solves for
        // the corresponding RHS, and then checks the 50 solution values.
        let num_col: HighsInt = self.lp_.num_col_;
        let num_row: HighsInt = self.lp_.num_row_;
        // PORT NOTE (perf): recycled storage, same values as fresh vectors (HVectorPool).
        let mut btran_rhs = self.hvector_pool_.take();
        let mut ftran_rhs = self.hvector_pool_.take();
        btran_rhs.setup(num_row);
        ftran_rhs.setup(num_row);

        // Solve for a random solution
        let mut random = HighsRandom::new(1);

        ftran_rhs.clear();
        let ideal_solution_num_nz: HighsInt = 50;
        let solution_num_nz: HighsInt = cpp_min(ideal_solution_num_nz, (num_row + 1) / 2);
        debug_assert!(solution_num_nz > 0);
        // PORT NOTE (perf): sized up front; C++ grows them by push_back.
        let mut solution_value: Vec<f64> = Vec::with_capacity(solution_num_nz as usize);
        let mut solution_index: Vec<HighsInt> = Vec::with_capacity(solution_num_nz as usize);
        let mut solution_nonzero: Vec<i8> = vec![0; num_row as usize];
        loop {
            let i_row = random.integer_sup(num_row);
            debug_assert!(i_row < num_row);
            if solution_nonzero[i_row as usize] != 0 {
                continue;
            }
            let value = random.fraction();
            solution_value.push(value);
            solution_index.push(i_row);
            solution_nonzero[i_row as usize] = 1;
            let i_col = self.basis_.basic_index_[i_row as usize];
            self.lp_.a_matrix_.collect_aj(&mut ftran_rhs, i_col, value);
            if solution_value.len() as HighsInt == solution_num_nz {
                break;
            }
        }

        btran_rhs.clear();
        let mut btran_solution: Vec<f64> = vec![0.0; num_row as usize];
        for i_x in 0..solution_value.len() {
            btran_solution[solution_index[i_x] as usize] = solution_value[i_x];
        }
        let mut btran_scattered_rhs: Vec<f64> = vec![0.0; (num_col + num_row) as usize];
        let ar_matrix = &self.ar_matrix_;
        for i_x in 0..solution_value.len() {
            let i_row = solution_index[i_x];
            for i_el in ar_matrix.p_end_[i_row as usize]..ar_matrix.start_[(i_row + 1) as usize] {
                let i_col = ar_matrix.index_[i_el as usize];
                btran_scattered_rhs[i_col as usize] +=
                    ar_matrix.value_[i_el as usize] * solution_value[i_x];
            }
            let i_col = num_col + i_row;
            if self.basis_.nonbasic_flag_[i_col as usize] == 0 {
                btran_scattered_rhs[i_col as usize] = solution_value[i_x];
            }
        }
        for i_row in 0..num_row {
            let i_col = self.basis_.basic_index_[i_row as usize];
            if btran_scattered_rhs[i_col as usize] == 0.0 {
                continue;
            }
            btran_rhs.array[i_row as usize] = btran_scattered_rhs[i_col as usize];
            btran_rhs.index[btran_rhs.count as usize] = i_row;
            btran_rhs.count += 1;
        }

        let expected_density: f64 = (solution_num_nz as f64) * self.info_.col_aq_density;
        self.ftran(&mut ftran_rhs, expected_density);
        self.btran(&mut btran_rhs, expected_density);

        let mut ftran_solution_error: f64 = 0.0;
        for i_x in 0..solution_value.len() {
            ftran_solution_error = cpp_max(
                (ftran_rhs.array[solution_index[i_x] as usize] - solution_value[i_x]).abs(),
                ftran_solution_error,
            );
        }
        let mut btran_solution_error: f64 = 0.0;
        for i_x in 0..solution_value.len() {
            btran_solution_error = cpp_max(
                (btran_rhs.array[solution_index[i_x] as usize] - solution_value[i_x]).abs(),
                btran_solution_error,
            );
        }
        self.hvector_pool_.give(btran_rhs);
        self.hvector_pool_.give(ftran_rhs);
        cpp_max(ftran_solution_error, btran_solution_error)
    }

    /// HiGHS: HEkk.cpp::HEkk::clearBadBasisChange (lines 3921-3933, 13 lines)
    /// C++ default arguments: reason = BadBasisChangeReason::kAll
    pub fn clear_bad_basis_change(&mut self, reason: BadBasisChangeReason) {
        if reason == BadBasisChangeReason::kAll {
            self.bad_basis_change_.clear();
        } else {
            self.bad_basis_change_
                .retain(|record| record.reason != reason);
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::updateBadBasisChange (lines 3935-3946, 12 lines)
    pub fn update_bad_basis_change(&mut self, col_aq: &HVector, theta_primal: f64) {
        if !self.bad_basis_change_.is_empty() {
            let primal_feasibility_tolerance = self.options_.primal_feasibility_tolerance;
            // remove_if removes a record when the predicate holds, so keep
            // the records for which it does not (`!(... >= tol)`, not `< tol`,
            // to match the C++ result on NaN).
            self.bad_basis_change_.retain(|record| {
                !((col_aq.array[record.row_out as usize] * theta_primal).abs()
                    >= primal_feasibility_tolerance)
            });
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::addBadBasisChange (lines 3948-3986, 39 lines)
    /// C++ default arguments: taboo = false
    pub fn add_bad_basis_change(
        &mut self,
        row_out: HighsInt,
        variable_out: HighsInt,
        variable_in: HighsInt,
        reason: BadBasisChangeReason,
        taboo: bool,
    ) -> HighsInt {
        debug_assert!(0 <= row_out && row_out <= self.lp_.num_row_);
        debug_assert!(0 <= variable_out && variable_out <= self.lp_.num_col_ + self.lp_.num_row_);
        if variable_in >= 0 {
            debug_assert!(0 <= variable_in && variable_in <= self.lp_.num_col_ + self.lp_.num_row_);
        } else {
            debug_assert!(variable_in == -1);
        }
        // Check that this is not already on the list
        let num_bad_basis_change: HighsInt = self.bad_basis_change_.len() as HighsInt;
        let mut bad_basis_change_num: HighsInt = -1;
        for ix in 0..num_bad_basis_change {
            let record = &self.bad_basis_change_[ix as usize];
            if record.row_out == row_out
                && record.variable_out == variable_out
                && record.variable_in == variable_in
                && record.reason == reason
            {
                bad_basis_change_num = ix;
                break;
            }
        }
        if bad_basis_change_num < 0 {
            // Not on the list so create record
            let record = HighsSimplexBadBasisChangeRecord {
                taboo,
                row_out,
                variable_out,
                variable_in,
                reason,
                ..Default::default()
            };
            self.bad_basis_change_.push(record);
            bad_basis_change_num = self.bad_basis_change_.len() as HighsInt - 1;
        } else {
            // On the list so just update whether it is taboo
            self.bad_basis_change_[bad_basis_change_num as usize].taboo = taboo;
        }
        bad_basis_change_num
    }

    /// HiGHS: HEkk.cpp::HEkk::clearBadBasisChangeTabooFlag (lines 3988-3990, 3 lines)
    pub fn clear_bad_basis_change_taboo_flag(&mut self) {
        for change in self.bad_basis_change_.iter_mut() {
            change.taboo = false;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::tabooBadBasisChange (lines 3992-3997, 6 lines)
    pub fn taboo_bad_basis_change(&self) -> bool {
        for change in self.bad_basis_change_.iter() {
            if change.taboo {
                return true;
            }
        }
        false
    }

    /// HiGHS: HEkk.cpp::HEkk::applyTabooRowOut (lines 3999-4009, 11 lines)
    pub fn apply_taboo_row_out(&mut self, values: &mut Vec<f64>, overwrite_with: f64) {
        debug_assert!(values.len() >= self.lp_.num_row_ as usize);
        for change in self.bad_basis_change_.iter_mut() {
            if change.taboo {
                let i_row = change.row_out as usize;
                change.save_value = values[i_row];
                values[i_row] = overwrite_with;
            }
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::unapplyTabooRowOut (lines 4011-4020, 10 lines)
    pub fn unapply_taboo_row_out(&mut self, values: &mut Vec<f64>) {
        debug_assert!(values.len() as HighsInt >= self.lp_.num_row_);
        // Unapply taboo rows in opposite order in case the row appears
        // twice in the list. This way the first saved value for the row is
        // what remains, not overwrite_with
        for change in self.bad_basis_change_.iter().rev() {
            if change.taboo {
                values[change.row_out as usize] = change.save_value;
            }
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::applyTabooVariableIn (lines 4022-4033, 12 lines)
    pub fn apply_taboo_variable_in(&mut self, values: &mut Vec<f64>, overwrite_with: f64) {
        debug_assert!(values.len() >= (self.lp_.num_col_ as usize) + (self.lp_.num_row_ as usize));
        for change in self.bad_basis_change_.iter_mut() {
            if change.taboo {
                let i_col = change.variable_in as usize;
                change.save_value = values[i_col];
                values[i_col] = overwrite_with;
            }
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::unapplyTabooVariableIn (lines 4035-4045, 11 lines)
    pub fn unapply_taboo_variable_in(&mut self, values: &mut Vec<f64>) {
        debug_assert!(values.len() as HighsInt >= self.lp_.num_col_ + self.lp_.num_row_);
        // Unapply taboo variables in opposite order in case the row appears
        // twice in the list. This way the first saved value for the
        // variable is what remains, not overwrite_with
        for change in self.bad_basis_change_.iter().rev() {
            if change.taboo {
                values[change.variable_in as usize] = change.save_value;
            }
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::logicalBasis (lines 4047-4052, 6 lines)
    pub fn logical_basis(&self) -> bool {
        for i_row in 0..self.lp_.num_row_ {
            if self.basis_.basic_index_[i_row as usize] < self.lp_.num_col_ {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod h16_rebuild_reason_tests {
    //! HiGHS 1.16 #3351 (h16_rebuild_reason) tests: the message mapping, plus the end-to-end
    //! TestEkk.cpp "Ekk-excessive-primal-value" case (3351.diff lines 1-43).

    use super::*;
    use crate::io::highs_io::kHighsLogDevLevelInfo;
    use crate::lp_data::highs_callback::{
        HighsCallbackFunctionType, HighsCallbackInput, HighsCallbackOutput,
    };
    use crate::lp_data::highs_run::Highs;
    use std::sync::{Arc, Mutex};

    /// TestEkk.cpp "Ekk-excessive-primal-value" (3351.diff lines 5-37): returns the
    /// run status, the captured kCallbackLogging text and the last rebuild reason
    /// string set by HEkkDual::reportRebuild.
    fn excessive_primal_value_run(flag: bool) -> (HighsStatus, String, String) {
        let log = Arc::new(Mutex::new(String::new()));
        let sink = log.clone();
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("log_to_console", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_int("log_dev_level", kHighsLogDevLevelInfo),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("presolve", "off"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", "simplex"),
            HighsStatus::kOk
        );
        assert_eq!(highs.set_option_value_int("threads", 1), HighsStatus::kOk);
        assert_eq!(
            highs.set_option_value_bool("h16_rebuild_reason", flag),
            HighsStatus::kOk
        );
        let cb = HighsCallbackFunctionType::new(
            move |callback_type: HighsInt,
                  message: &str,
                  _o: &HighsCallbackOutput,
                  _i: &mut HighsCallbackInput| {
                if callback_type == kCallbackLogging {
                    sink.lock().unwrap().push_str(message);
                }
            },
        );
        assert_eq!(highs.set_callback(cb), HighsStatus::kOk);
        assert_eq!(highs.start_callback(kCallbackLogging), HighsStatus::kOk);

        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 3;
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.col_upper_ = vec![kHighsInf, kHighsInf];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![0.0, 0.0, -1e16];
        lp.a_matrix_.start_ = vec![0, 2, 5];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 2];
        lp.a_matrix_.value_ = vec![1e-6, -1e6, -2.0, -1.0, 1.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        let status = highs.run();
        let text = log.lock().unwrap().clone();
        let reason = highs.ekk_instance_.analysis_.rebuild_reason_string.clone();
        (status, text, reason)
    }

    #[test]
    fn h16_rebuild_reason_excessive_primal_value_run_on() {
        // TestEkk.cpp: run() returns kError and the log names "Excessive primal value".
        // The C++ log line comes from HighsSimplexAnalysis::invertReport (dev level),
        // which reportRebuild calls with analysis->rebuild_reason_string. The Rust
        // invert_report is a no-op stub (simplex/analysis.rs), so the test checks the
        // string reportRebuild hands to it, plus the user error that ends phase 2.
        let (status, log, reason) = excessive_primal_value_run(true);
        assert_eq!(status, HighsStatus::kError);
        assert!(
            log.contains("Dual simplex detected excessive primal values"),
            "{log}"
        );
        assert_eq!(reason, "Excessive primal value");
    }

    #[test]
    #[cfg(not(debug_assertions))]
    fn h16_rebuild_reason_excessive_primal_value_run_off() {
        // Same run with the flag off: 1.15.1 names the reason "Unidentified".
        let (status, log, reason) = excessive_primal_value_run(false);
        assert_eq!(status, HighsStatus::kError);
        assert!(
            log.contains("Dual simplex detected excessive primal values"),
            "{log}"
        );
        assert_eq!(reason, "Unidentified");
    }

    #[test]
    fn h16_rebuild_reason_excessive_primal_value_on() {
        let mut ekk = HEkk::default();
        ekk.options_.h16_rebuild_reason = true;
        assert_eq!(
            ekk.rebuild_reason(kRebuildReasonExcessivePrimalValue),
            "Excessive primal value"
        );
        // Existing reasons are unchanged with the flag on
        assert_eq!(
            ekk.rebuild_reason(kRebuildReasonCleanup),
            "Perturbation cleanup"
        );
    }

    #[test]
    #[cfg(not(debug_assertions))]
    fn h16_rebuild_reason_excessive_primal_value_off() {
        // With the flag off the 1.15.1 chain is unchanged: the reason is not named.
        // The 1.15.1 debug_assert only fires in debug builds, so this runs in release.
        let ekk = HEkk::default();
        assert!(!ekk.options_.h16_rebuild_reason);
        assert_eq!(
            ekk.rebuild_reason(kRebuildReasonExcessivePrimalValue),
            "Unidentified"
        );
    }
}
