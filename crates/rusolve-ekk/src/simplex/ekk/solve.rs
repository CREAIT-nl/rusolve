//! HiGHS source: simplex/HEkk.cpp (lines 1002-1693)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Solve driver, basis set/get, add/delete cols/rows, unscaling, solution, basis+factor init.

use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkk.cpp::HEkk::solve (lines 1002-1120, 119 lines)
    /// C++ default arguments: force_phase2 = false
    /// Dual: `HEkkDual dual_solver(*this); call_status = dual_solver.solve(force_phase2);` ->
    /// `let mut dual_solver = HEkkDual::new(self); call_status = dual_solver.solve(self, force_phase2);`
    /// (primal likewise with `HEkkPrimal::new(self)`).
    pub fn solve(&mut self, force_phase2: bool) -> HighsStatus {
        self.debug_initialise();

        self.initialise_analysis();
        self.initialise_control();

        self.dual_simplex_cleanup_level_ = 0;
        self.dual_simplex_phase1_cleanup_level_ = 0;

        self.previous_iteration_cycling_detected = -kHighsIInf;

        self.initialise_for_solve();

        let simplex_nla_status = {
            let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
            self.simplex_nla_
                .debug_check_data(&refs, "Before HEkk::solve()")
        };
        let simplex_nla_ok = simplex_nla_status == HighsDebugStatus::kOk;
        if !simplex_nla_ok {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Error in simplex NLA data\n"
            );
            debug_assert!(simplex_nla_ok);
            return self.return_from_ekk_solve(HighsStatus::kError);
        }

        let report_initial_basis = false;
        if report_initial_basis {
            self.debug_report_initial_basis();
        }

        debug_assert!(self.status_.has_basis);
        debug_assert!(self.status_.has_invert);
        debug_assert!(self.status_.initialised_for_solve);
        if self.model_status_ == HighsModelStatus::kOptimal {
            return self.return_from_ekk_solve(HighsStatus::kOk);
        }

        let mut return_status = HighsStatus::kOk;
        let mut call_status: HighsStatus;
        let algorithm_name: &str;

        // Indicate that dual and primal rays are not known
        self.clear_ray_records();

        // Allow primal and dual perturbations in case a block on them is
        // hanging over from a previous call
        self.info_.allow_cost_shifting = true;
        self.info_.allow_cost_perturbation = true;
        self.info_.allow_bound_perturbation = true;

        // C++ passes (*options_, info_), both members of this
        Self::choose_simplex_strategy_threads(&self.options_, &mut self.info_);
        let simplex_strategy = self.info_.simplex_strategy;

        // Initial solve according to strategy
        if simplex_strategy == kSimplexStrategyPrimal {
            algorithm_name = "primal";
            report_simplex_phase_iterations(
                &self.options_.log_options,
                self.iteration_count_,
                &mut self.info_,
                true,
            );
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Using {}\n",
                self.simplex_strategy_to_string(kSimplexStrategyPrimal)
            );
            let mut primal_solver = HEkkPrimal::new(self);
            call_status = primal_solver.solve(self, force_phase2);
            debug_assert!(self.called_return_from_solve_);
            return_status = interpret_call_status(
                &self.options_.log_options,
                call_status,
                return_status,
                "HEkkPrimal::solve",
            );
        } else {
            algorithm_name = "dual";
            report_simplex_phase_iterations(
                &self.options_.log_options,
                self.iteration_count_,
                &mut self.info_,
                true,
            );
            // Solve, depending on the particular strategy
            if simplex_strategy == kSimplexStrategyDualTasks {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Using {} with concurrency of {}\n",
                    self.simplex_strategy_to_string(kSimplexStrategyDualTasks),
                    self.info_.num_concurrency
                );
            } else if simplex_strategy == kSimplexStrategyDualMulti {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Using {} with concurrency of {}\n",
                    self.simplex_strategy_to_string(kSimplexStrategyDualMulti),
                    self.info_.num_concurrency
                );
            } else {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Using {}\n",
                    self.simplex_strategy_to_string(kSimplexStrategyDual)
                );
            }
            let mut dual_solver = HEkkDual::new(self);
            call_status = dual_solver.solve(self, force_phase2);
            dual_solver.recycle(self);
            debug_assert!(self.called_return_from_solve_);
            return_status = interpret_call_status(
                &self.options_.log_options,
                call_status,
                return_status,
                "HEkkDual::solve",
            );

            // Dual simplex solver may set model_status to be
            // kUnboundedOrInfeasible, and Highs::run() may not allow that to
            // be returned, so use primal simplex to distinguish
            if self.model_status_ == HighsModelStatus::kUnboundedOrInfeasible
                && !self.options_.allow_unbounded_or_infeasible
            {
                let mut primal_solver = HEkkPrimal::new(self);
                call_status = primal_solver.solve(self, false);
                debug_assert!(self.called_return_from_solve_);
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "HEkkPrimal::solve",
                );
            }
        }

        report_simplex_phase_iterations(
            &self.options_.log_options,
            self.iteration_count_,
            &mut self.info_,
            false,
        );
        if return_status == HighsStatus::kError {
            return self.return_from_ekk_solve(return_status);
        }
        crate::highs_log_dev!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "{} simplex solver returns {} primal and {} dual infeasibilities: Status {}\n",
            algorithm_name,
            self.info_.num_primal_infeasibilities,
            self.info_.num_dual_infeasibilities,
            util_model_status_to_string(self.model_status_)
        );
        // Can model_status_ = HighsModelStatus::kNotset be returned?
        debug_assert!(self.model_status_ != HighsModelStatus::kNotset);

        self.return_from_ekk_solve(return_status)
    }

    /// HiGHS: HEkk.cpp::HEkk::setBasis (lines 1122-1173, 52 lines)
    pub fn set_basis(&mut self) -> HighsStatus {
        // Set up nonbasicFlag and basicIndex for a logical basis
        let num_col = self.lp_.num_col_;
        let num_row = self.lp_.num_row_;

        self.basis_.setup(num_col, num_row);
        self.basis_.debug_origin_name = "HEkk::setBasis - logical".to_string();

        for i_col in 0..num_col {
            self.basis_.nonbasic_flag_[i_col as usize] = kNonbasicFlagTrue;
            let lower = self.lp_.col_lower_[i_col as usize];
            let upper = self.lp_.col_upper_[i_col as usize];
            let mut move_: i8 = kIllegalMoveValue;
            if lower == upper {
                // Fixed
                move_ = kNonbasicMoveZe;
            } else if !highs_is_infinity(-lower) {
                // Finite lower bound so boxed or lower
                if !highs_is_infinity(upper) {
                    // Finite upper bound so boxed. Set to bound of LP that is closer to
                    // zero
                    if move_ == kIllegalMoveValue {
                        if lower.abs() < upper.abs() {
                            move_ = kNonbasicMoveUp;
                        } else {
                            move_ = kNonbasicMoveDn;
                        }
                    }
                } else {
                    // Lower (since upper bound is infinite)
                    move_ = kNonbasicMoveUp;
                }
            } else if !highs_is_infinity(upper) {
                // Upper
                move_ = kNonbasicMoveDn;
            } else {
                // FREE
                move_ = kNonbasicMoveZe;
            }
            debug_assert!(move_ != kIllegalMoveValue);
            self.basis_.nonbasic_move_[i_col as usize] = move_;
        }
        for i_row in 0..num_row {
            let i_var = num_col + i_row;
            self.basis_.nonbasic_flag_[i_var as usize] = kNonbasicFlagFalse;
            sparse_combine(&mut self.basis_.hash, i_var);
            self.basis_.basic_index_[i_row as usize] = i_var;
        }
        self.info_.num_basic_logicals = num_row;
        self.status_.has_basis = true;
        HighsStatus::kOk
    }

    /// HiGHS: HEkk.cpp::HEkk::setBasis (lines 1175-1250, 76 lines)
    pub fn set_basis_from(&mut self, highs_basis: &HighsBasis) -> HighsStatus {
        // Shouldn't have to check the incoming basis since this is an
        // internal call, but it may be a basis that's set up internally
        // with errors :-) ...
        //
        if kDebugMipNodeDualFeasible {
            // The basis should be dual feasible unless it was alien
            self.debug_dual_feasible = !highs_basis.was_alien;
        } else {
            debug_assert!(!self.debug_dual_feasible);
        }
        if debug_highs_basis_consistent(&self.options_, &self.lp_, highs_basis)
            == HighsDebugStatus::kLogicalError
        {
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Supposed to be a Highs basis, but not valid\n"
            );
            return HighsStatus::kError;
        }
        let num_col = self.lp_.num_col_;
        let num_row = self.lp_.num_row_;
        // Set up the basis in case it has not yet been done for this LP
        self.basis_.setup(num_col, num_row);
        self.basis_.debug_id = highs_basis.debug_id;
        self.basis_.debug_update_count = highs_basis.debug_update_count;
        self.basis_.debug_origin_name = highs_basis.debug_origin_name.clone();
        debug_assert!(!self.basis_.debug_origin_name.is_empty());
        let mut num_basic_variables: HighsInt = 0;
        for i_col in 0..num_col {
            let i_var = i_col;
            let lower = self.lp_.col_lower_[i_col as usize];
            let upper = self.lp_.col_upper_[i_col as usize];
            if highs_basis.col_status[i_col as usize] == HighsBasisStatus::kBasic {
                self.basis_.nonbasic_flag_[i_var as usize] = kNonbasicFlagFalse;
                self.basis_.nonbasic_move_[i_var as usize] = 0;
                self.basis_.basic_index_[num_basic_variables as usize] = i_var;
                num_basic_variables += 1;
                sparse_combine(&mut self.basis_.hash, i_var);
            } else {
                self.basis_.nonbasic_flag_[i_var as usize] = kNonbasicFlagTrue;
                if lower == upper {
                    self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveZe;
                } else if highs_basis.col_status[i_col as usize] == HighsBasisStatus::kLower {
                    self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveUp;
                } else if highs_basis.col_status[i_col as usize] == HighsBasisStatus::kUpper {
                    self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveDn;
                } else {
                    debug_assert!(
                        highs_basis.col_status[i_col as usize] == HighsBasisStatus::kZero
                    );
                    self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveZe;
                }
            }
        }
        for i_row in 0..num_row {
            let i_var = num_col + i_row;
            let lower = self.lp_.row_lower_[i_row as usize];
            let upper = self.lp_.row_upper_[i_row as usize];
            if highs_basis.row_status[i_row as usize] == HighsBasisStatus::kBasic {
                self.basis_.nonbasic_flag_[i_var as usize] = kNonbasicFlagFalse;
                self.basis_.nonbasic_move_[i_var as usize] = 0;
                self.basis_.basic_index_[num_basic_variables as usize] = i_var;
                num_basic_variables += 1;
                sparse_combine(&mut self.basis_.hash, i_var);
            } else {
                self.basis_.nonbasic_flag_[i_var as usize] = kNonbasicFlagTrue;
                if lower == upper {
                    self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveZe;
                } else if highs_basis.row_status[i_row as usize] == HighsBasisStatus::kLower {
                    self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveDn;
                } else if highs_basis.row_status[i_row as usize] == HighsBasisStatus::kUpper {
                    self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveUp;
                } else {
                    debug_assert!(
                        highs_basis.row_status[i_row as usize] == HighsBasisStatus::kZero
                    );
                    self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveZe;
                }
            }
        }
        self.status_.has_basis = true;
        HighsStatus::kOk
    }

    /// HiGHS: HEkk.cpp::HEkk::addCols (lines 1252-1266, 15 lines)
    pub fn add_cols(&mut self, lp: &HighsLp, scaled_a_matrix: &HighsSparseMatrix) {
        if self.status_.has_nla {
            self.simplex_nla_.add_cols(lp);
        }
        self.update_status(LpAction::kNewCols);
    }

    /// HiGHS: HEkk.cpp::HEkk::addRows (lines 1268-1291, 24 lines)
    pub fn add_rows(&mut self, lp: &HighsLp, scaled_ar_matrix: &HighsSparseMatrix) {
        if kExtendInvertWhenAddingRows && self.status_.has_nla {
            self.simplex_nla_
                .add_rows(lp, &mut self.basis_.basic_index_, scaled_ar_matrix);
            self.set_nla_pointers_for_trans(lp);
            self.debug_nla_check_invert("HEkk::addRows - on entry", kHighsDebugLevelExpensive + 1);
        }
        // Update the number of rows in the simplex LP so that it's
        // consistent with simplex basis information
        self.lp_.num_row_ = lp.num_row_;
        self.update_status(LpAction::kNewRows);
    }

    /// HiGHS: HEkk.cpp::HEkk::deleteCols (lines 1293-1295, 3 lines)
    pub fn delete_cols(&mut self, index_collection: &HighsIndexCollection) {
        self.update_status(LpAction::kDelCols);
    }

    /// HiGHS: HEkk.cpp::HEkk::deleteRows (lines 1296-1298, 3 lines)
    pub fn delete_rows(&mut self, index_collection: &HighsIndexCollection) {
        // rusolve lever `mip_lp_keep_dse_weights` (not in HiGHS, dse_keep.rs).
        if self.options_.mip_lp_keep_dse_weights > 0 {
            self.dse_keep_save();
            self.dse_keep_delete_rows(index_collection);
        }
        self.update_status(LpAction::kDelRows);
    }

    /// HiGHS: HEkk.cpp::HEkk::unscaleSimplex (lines 1300-1346, 47 lines)
    pub fn unscale_simplex(&mut self, incumbent_lp: &HighsLp) {
        if !self.simplex_in_scaled_space_ {
            return;
        }
        debug_assert!(incumbent_lp.scale_.has_scaling);
        let num_col = incumbent_lp.num_col_;
        let num_row = incumbent_lp.num_row_;
        let col_scale = &incumbent_lp.scale_.col;
        let row_scale = &incumbent_lp.scale_.row;
        for i_col in 0..num_col {
            let i_var = i_col as usize;
            let factor = col_scale[i_col as usize];
            self.info_.work_cost_[i_var] /= factor;
            self.info_.work_dual_[i_var] /= factor;
            self.info_.work_shift_[i_var] /= factor;
            self.info_.work_lower_[i_var] *= factor;
            self.info_.work_upper_[i_var] *= factor;
            self.info_.work_range_[i_var] *= factor;
            self.info_.work_value_[i_var] *= factor;
            self.info_.work_lower_shift_[i_var] *= factor;
            self.info_.work_upper_shift_[i_var] *= factor;
        }
        for i_row in 0..num_row {
            let i_var = (num_col + i_row) as usize;
            let factor = row_scale[i_row as usize];
            self.info_.work_cost_[i_var] *= factor;
            self.info_.work_dual_[i_var] *= factor;
            self.info_.work_shift_[i_var] *= factor;
            self.info_.work_lower_[i_var] /= factor;
            self.info_.work_upper_[i_var] /= factor;
            self.info_.work_range_[i_var] /= factor;
            self.info_.work_value_[i_var] /= factor;
            self.info_.work_lower_shift_[i_var] /= factor;
            self.info_.work_upper_shift_[i_var] /= factor;
        }
        for i_row in 0..num_row {
            let i_var = self.basis_.basic_index_[i_row as usize];
            let factor: f64 = if i_var < num_col {
                col_scale[i_var as usize]
            } else {
                1.0 / row_scale[(i_var - num_col) as usize]
            };
            self.info_.base_lower_[i_row as usize] *= factor;
            self.info_.base_upper_[i_row as usize] *= factor;
            self.info_.base_value_[i_row as usize] *= factor;
        }
        self.simplex_in_scaled_space_ = false;
    }

    /// HiGHS: HEkk.cpp::HEkk::getSolution (lines 1348-1376, 29 lines)
    pub fn get_solution(&mut self) -> HighsSolution {
        let mut solution = HighsSolution::default();
        // Scatter the basic primal values
        for i_row in 0..self.lp_.num_row_ {
            self.info_.work_value_[self.basis_.basic_index_[i_row as usize] as usize] =
                self.info_.base_value_[i_row as usize];
        }
        // Zero the basic dual values
        for i_row in 0..self.lp_.num_row_ {
            self.info_.work_dual_[self.basis_.basic_index_[i_row as usize] as usize] = 0.0;
        }

        // Now we can get the solution
        solution.col_value.resize(self.lp_.num_col_ as usize, 0.0);
        solution.col_dual.resize(self.lp_.num_col_ as usize, 0.0);
        solution.row_value.resize(self.lp_.num_row_ as usize, 0.0);
        solution.row_dual.resize(self.lp_.num_row_ as usize, 0.0);

        let sense = self.lp_.sense_ as HighsInt;
        for i_col in 0..self.lp_.num_col_ {
            solution.col_value[i_col as usize] = self.info_.work_value_[i_col as usize];
            solution.col_dual[i_col as usize] =
                (sense as f64) * self.info_.work_dual_[i_col as usize];
        }
        for i_row in 0..self.lp_.num_row_ {
            solution.row_value[i_row as usize] =
                -self.info_.work_value_[(self.lp_.num_col_ + i_row) as usize];
            // @FlipRowDual negate RHS
            solution.row_dual[i_row as usize] =
                ((-sense) as f64) * self.info_.work_dual_[(self.lp_.num_col_ + i_row) as usize];
        }
        solution.value_valid = true;
        solution.dual_valid = true;
        solution
    }

    /// HiGHS: HEkk.cpp::HEkk::getHighsBasis (lines 1378-1439, 62 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn get_highs_basis(&self, use_lp: &mut HighsLp) -> HighsBasis {
        let num_col = use_lp.num_col_;
        let num_row = use_lp.num_row_;
        let mut highs_basis = HighsBasis::default();
        highs_basis
            .col_status
            .resize(num_col as usize, HighsBasisStatus::kLower);
        highs_basis
            .row_status
            .resize(num_row as usize, HighsBasisStatus::kLower);
        debug_assert!(self.status_.has_basis);
        highs_basis.valid = false;
        let sense = self.lp_.sense_ as HighsInt;
        for i_col in 0..num_col {
            let i_var = i_col;
            let lower = use_lp.col_lower_[i_col as usize];
            let upper = use_lp.col_upper_[i_col as usize];
            let mut basis_status = HighsBasisStatus::kNonbasic;
            if self.basis_.nonbasic_flag_[i_var as usize] == 0 {
                basis_status = HighsBasisStatus::kBasic;
            } else if self.basis_.nonbasic_move_[i_var as usize] == kNonbasicMoveUp {
                basis_status = HighsBasisStatus::kLower;
            } else if self.basis_.nonbasic_move_[i_var as usize] == kNonbasicMoveDn {
                basis_status = HighsBasisStatus::kUpper;
            } else if self.basis_.nonbasic_move_[i_var as usize] == kNonbasicMoveZe {
                if lower == upper {
                    let dual = (sense as f64) * self.info_.work_dual_[i_col as usize];
                    basis_status = if dual >= 0.0 {
                        HighsBasisStatus::kLower
                    } else {
                        HighsBasisStatus::kUpper
                    };
                } else {
                    basis_status = HighsBasisStatus::kZero;
                }
            }
            highs_basis.col_status[i_col as usize] = basis_status;
        }
        for i_row in 0..num_row {
            let i_var = num_col + i_row;
            let lower = use_lp.row_lower_[i_row as usize];
            let upper = use_lp.row_upper_[i_row as usize];
            let mut basis_status = HighsBasisStatus::kNonbasic;
            if self.basis_.nonbasic_flag_[i_var as usize] == 0 {
                basis_status = HighsBasisStatus::kBasic;
            } else if self.basis_.nonbasic_move_[i_var as usize] == kNonbasicMoveUp {
                basis_status = HighsBasisStatus::kUpper;
            } else if self.basis_.nonbasic_move_[i_var as usize] == kNonbasicMoveDn {
                basis_status = HighsBasisStatus::kLower;
            } else if self.basis_.nonbasic_move_[i_var as usize] == kNonbasicMoveZe {
                if lower == upper {
                    let dual = (sense as f64) * self.info_.work_dual_[i_var as usize];
                    basis_status = if dual >= 0.0 {
                        HighsBasisStatus::kLower
                    } else {
                        HighsBasisStatus::kUpper
                    };
                } else {
                    basis_status = HighsBasisStatus::kZero;
                }
            }
            highs_basis.row_status[i_row as usize] = basis_status;
        }
        highs_basis.valid = true;
        highs_basis.alien = false;
        highs_basis.useful = true;
        highs_basis.was_alien = false;
        highs_basis.debug_id =
            (self.build_synthetic_tick_ + self.total_synthetic_tick_) as HighsInt;
        highs_basis.debug_update_count = self.info_.update_count;
        highs_basis.debug_origin_name = self.basis_.debug_origin_name.clone();
        highs_basis
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseSimplexLpBasisAndFactor (lines 1441-1537, 97 lines)
    /// C++ default arguments: only_from_known_basis = false
    /// `simplex_nla_.setup(&lp_, basis_.basicIndex_.data(), options_, timer_, &analysis_, factor_a_matrix, thr)`
    /// -> `self.simplex_nla_.setup(&self.lp_, &self.basis_.basic_index_, &self.options_, factor_a_matrix, thr)`;
    /// `getScaledAMatrixPointer()` borrows `self` mutably, so clone or restructure the borrow
    /// (e.g. compute the scaled matrix into `self.scaled_a_matrix_` first, then borrow fields disjointly).
    pub fn initialise_simplex_lp_basis_and_factor(
        &mut self,
        only_from_known_basis: bool,
    ) -> HighsStatus {
        // If only_from_known_basis is true, then there should be a simplex
        // basis to use
        if only_from_known_basis {
            debug_assert!(self.status_.has_basis);
        }
        // If there is no simplex basis, set up a logical basis
        if !self.status_.has_basis {
            self.set_basis();
        }
        // getScaledAMatrixPointer() returns a pointer to either the
        // constraint matrix or a scaled copy (that is a member of the HEkk
        // class). Call for its side effect, then borrow disjointly.
        let _ = self.get_scaled_a_matrix_pointer();
        let local_scaled_a_matrix: &HighsSparseMatrix =
            if self.lp_.scale_.has_scaling && !self.lp_.is_scaled_ {
                &self.scaled_a_matrix_
            } else {
                &self.lp_.a_matrix_
            };
        //
        // If simplex NLA is set up, pass the pointers that it uses.
        if self.status_.has_nla {
            debug_assert!(self.lp_factor_row_compatible());
            self.simplex_nla_
                .set_pointers(&self.lp_, Some(local_scaled_a_matrix));
        } else {
            // todo @ Julian: this fails on glass4
            debug_assert!(
                self.info_.factor_pivot_threshold >= self.options_.factor_pivot_threshold
            );
            self.simplex_nla_.setup(
                &self.lp_,
                &self.basis_.basic_index_,
                &self.options_,
                local_scaled_a_matrix,
                self.info_.factor_pivot_threshold,
            );
            self.status_.has_nla = true;
        }

        if !self.status_.has_invert {
            let rank_deficiency = self.compute_factor();
            if rank_deficiency != 0 {
                // Basis is rank deficient
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "HEkk::initialiseSimplexLpBasisAndFactor ({}) Rank_deficiency {}: Id = {}; UpdateCount = {}\n",
                    self.basis_.debug_origin_name,
                    rank_deficiency,
                    self.basis_.debug_id,
                    self.basis_.debug_update_count
                );
                if only_from_known_basis {
                    // If only this basis should be used, then return error
                    crate::highs_log_dev!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Supposed to be a full-rank basis, but incorrect\n"
                    );
                    return HighsStatus::kError;
                }
                // Account for rank deficiency by correcting nonbasicFlag
                self.handle_rank_deficiency();
                self.update_status(LpAction::kNewBasis);
                self.set_nonbasic_move();
                self.status_.has_basis = true;
                self.status_.has_invert = true;
                self.status_.has_fresh_invert = true;
            }
            // Record the synthetic clock for INVERT, and zero it for UPDATE
            self.reset_synthetic_clock();
        }
        debug_assert!(self.status_.has_invert);
        HighsStatus::kOk
    }

    /// HiGHS: HEkk.cpp::HEkk::handleRankDeficiency (lines 1539-1569, 31 lines)
    pub fn handle_rank_deficiency(&mut self) {
        let rank_deficiency = self.simplex_nla_.factor_.rank_deficiency;
        for k in 0..rank_deficiency {
            let row_in = self.simplex_nla_.factor_.row_with_no_pivot[k as usize];
            let variable_in = self.lp_.num_col_ + row_in;
            let variable_out = self.simplex_nla_.factor_.var_with_no_pivot[k as usize];
            self.basis_.nonbasic_flag_[variable_in as usize] = kNonbasicFlagFalse;
            self.basis_.nonbasic_flag_[variable_out as usize] = kNonbasicFlagTrue;
            let row_out = self.simplex_nla_.factor_.row_with_no_pivot[k as usize];
            debug_assert!(self.basis_.basic_index_[row_out as usize] == variable_in);
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "HEkk::handleRankDeficiency: {:4}: Basic row of leaving variable ({:4} is {} {:4}) is {:4}; Entering logical = {:4} is variable {})\n",
                k,
                variable_out,
                if variable_out < self.lp_.num_col_ {
                    " column"
                } else {
                    "logical"
                },
                if variable_out < self.lp_.num_col_ {
                    variable_out
                } else {
                    variable_out - self.lp_.num_col_
                },
                row_out,
                row_in,
                variable_in
            );
            // NB Parameters row_out, variable_in, variable_out, since
            // variable_in is the logical that must not come out to be
            // replaced by the structural variable_out
            self.add_bad_basis_change(
                row_out,
                variable_in,
                variable_out,
                BadBasisChangeReason::kSingular,
                true,
            );
        }
        self.status_.has_ar_matrix = false;
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseEkk (lines 1573-1581, 9 lines)
    /// MUST-PORT (determinism, random stream): first `initialise_simplex_lp_random_vectors` draw; its values are overwritten but the draws advance `random_`.
    pub fn initialise_ekk(&mut self) {
        if self.status_.initialised_for_new_lp {
            return;
        }
        self.set_simplex_options();
        self.initialise_control();
        self.initialise_simplex_lp_random_vectors();
        self.simplex_nla_.clear();
        self.clear_bad_basis_change(BadBasisChangeReason::kAll);
        self.status_.initialised_for_new_lp = true;
    }

    /// HiGHS: HEkk.cpp::HEkk::isUnconstrainedLp (lines 1583-1593, 11 lines)
    pub fn is_unconstrained_lp(&self) -> bool {
        let is_unconstrained_lp = self.lp_.num_row_ <= 0;
        if is_unconstrained_lp {
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kError,
                "HEkkDual::solve called for LP with non-positive ({}) number of constraints\n",
                self.lp_.num_row_
            );
        }
        debug_assert!(!is_unconstrained_lp);
        is_unconstrained_lp
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseForSolve (lines 1595-1621, 27 lines)
    /// MUST-PORT (determinism, random stream): second `initialise_simplex_lp_random_vectors` draw on the same stream; this is the one the solve uses.
    pub fn initialise_for_solve(&mut self) {
        let return_status = self.initialise_simplex_lp_basis_and_factor(false);
        debug_assert!(return_status == HighsStatus::kOk);
        debug_assert!(self.status_.has_basis);

        self.update_simplex_options();
        self.initialise_simplex_lp_random_vectors();
        self.initialise_partitioned_rowwise_matrix(); // Timed
        self.allocate_work_and_base_arrays();
        self.initialise_cost(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown, false);
        self.initialise_bound(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown, false);
        self.initialise_nonbasic_value_and_move();
        self.compute_primal(); // Timed
        self.compute_dual(); // Timed
        self.compute_simplex_infeasible(); // Timed
        self.compute_dual_objective_value(2); // Timed
        self.compute_primal_objective_value(); // Timed
        self.status_.initialised_for_solve = true;

        let primal_feasible = self.info_.num_primal_infeasibilities == 0;
        let dual_feasible = self.info_.num_dual_infeasibilities == 0;
        self.visited_basis_.clear();
        self.visited_basis_.insert(self.basis_.hash);
        self.model_status_ = HighsModelStatus::kNotset;
        if primal_feasible && dual_feasible {
            self.model_status_ = HighsModelStatus::kOptimal;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::setSimplexOptions (lines 1623-1643, 21 lines)
    /// MUST-PORT (determinism, random stream): seeds `self.random_.initialise(options_.random_seed as HighsUInt)`; every later draw depends on this call happening exactly here.
    pub fn set_simplex_options(&mut self) {
        // Copy values of HighsOptions for the simplex solver
        //
        // NB simplex_strategy is set by chooseSimplexStrategyThreads in each call
        //
        self.info_.dual_edge_weight_strategy = self.options_.simplex_dual_edge_weight_strategy;
        self.info_.price_strategy = self.options_.simplex_price_strategy;
        self.info_.dual_simplex_cost_perturbation_multiplier =
            self.options_.dual_simplex_cost_perturbation_multiplier;
        self.info_.primal_simplex_bound_perturbation_multiplier =
            self.options_.primal_simplex_bound_perturbation_multiplier;
        self.info_.factor_pivot_threshold = self.options_.factor_pivot_threshold;
        self.info_.update_limit = self.options_.simplex_update_limit;
        self.random_
            .initialise(self.options_.random_seed as HighsUInt);

        // Set values of internal options
        self.info_.store_squared_primal_infeasibility = true;
    }

    /// HiGHS: HEkk.cpp::HEkk::updateSimplexOptions (lines 1645-1657, 13 lines)
    pub fn update_simplex_options(&mut self) {
        // Update some simplex option values from HighsOptions when
        // (re-)solving an LP.
        //
        // NB simplex_strategy is set by chooseSimplexStrategyThreads in each call
        //
        self.info_.dual_simplex_cost_perturbation_multiplier =
            self.options_.dual_simplex_cost_perturbation_multiplier;
        self.info_.primal_simplex_bound_perturbation_multiplier =
            self.options_.primal_simplex_bound_perturbation_multiplier;
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseSimplexLpRandomVectors (lines 1659-1693, 35 lines)
    /// MUST-PORT (determinism, random stream): draw order is: shuffle `num_col_permutation_` (only if num_col > 0), shuffle `num_tot_permutation_`, then `num_tot` fractions; all on `self.random_` (not a fresh generator).
    pub fn initialise_simplex_lp_random_vectors(&mut self) {
        let num_col = self.lp_.num_col_;
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        if num_tot == 0 {
            return;
        }
        // The random number generator is the member random_ (not re-initialised)
        let random = &mut self.random_;

        if num_col != 0 {
            // Generate a random permutation of the column indices
            let num_col_permutation = &mut self.info_.num_col_permutation_;
            num_col_permutation.resize(num_col as usize, 0);
            for i in 0..num_col {
                num_col_permutation[i as usize] = i;
            }
            random.shuffle(num_col_permutation, num_col);
        }

        // Generate a random permutation of all the indices
        let num_tot_permutation = &mut self.info_.num_tot_permutation_;
        num_tot_permutation.resize(num_tot as usize, 0);
        for i in 0..num_tot {
            num_tot_permutation[i as usize] = i;
        }
        random.shuffle(num_tot_permutation, num_tot);

        // Generate a vector of random reals
        self.info_
            .num_tot_random_value_
            .resize(num_tot as usize, 0.0);
        let num_tot_random_value = &mut self.info_.num_tot_random_value_;
        for i in 0..num_tot {
            num_tot_random_value[i as usize] = random.fraction();
        }
    }
}
