//! HiGHS source: simplex/HEkkDual.cpp (lines 1185-1940)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Iteration: CHUZR (`chooseRow`), CHUZC with BFRT (`chooseColumn`), DSE weight acceptance,
//! Devex framework.
//! See `HEkkDual` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::highs_log_dev;
use crate::parallel::highs_parallel::{for_each, spawn, sync};
use crate::prelude::*;

/// Raw pointer captured by a task closure (rule (s)): the C++ tasks share `this` and `ekk`
/// and touch disjoint members under the sync protocol of `iterateTasks`.
struct TaskPtr<T>(*mut T);
impl<T> Clone for TaskPtr<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for TaskPtr<T> {}
// SAFETY: the pointee is only accessed under the C++ task protocol (disjoint fields, joins
// after sync in slice order).
unsafe impl<T> Send for TaskPtr<T> {}
// SAFETY: as above.
unsafe impl<T> Sync for TaskPtr<T> {}
impl<T> TaskPtr<T> {
    fn get(self) -> *mut T {
        self.0
    }
}

impl HEkkDual {
    /// HiGHS: HEkkDual.cpp::HEkkDual::iterate (lines 1185-1268, 84 lines)
    pub fn iterate(&mut self, ekk: &mut HEkk) {
        // This is the main iteration loop for dual revised simplex. All the
        // methods have as their first line if (rebuild_reason) return;, where
        // rebuild_reason is, for example, set to 1 when CHUZR finds no
        // candidate. This causes a break from the inner loop of
        // solvePhase% and, hence, a call to rebuild()
        let from_check_iter: HighsInt = 0;
        let to_check_iter: HighsInt = from_check_iter + 100;
        if ekk.debug_solve_report_ {
            ekk.debug_iteration_report_ =
                ekk.iteration_count_ >= from_check_iter && ekk.iteration_count_ <= to_check_iter;
            if ekk.debug_iteration_report_ {
                println!("HEkkDual::iterate Debug iteration {}", ekk.iteration_count_);
            }
        }

        ekk.analysis_.simplex_timer_start(IterateChuzrClock, 0);
        self.choose_row(ekk);
        ekk.analysis_.simplex_timer_stop(IterateChuzrClock, 0);

        ekk.analysis_.simplex_timer_start(IterateChuzcClock, 0);
        let mut row_ep = std::mem::take(&mut self.row_ep);
        self.choose_column(ekk, &mut row_ep);
        self.row_ep = row_ep;
        ekk.analysis_.simplex_timer_stop(IterateChuzcClock, 0);

        if self.is_bad_basis_change(ekk) {
            return;
        }

        ekk.analysis_.simplex_timer_start(IterateFtranClock, 0);
        self.update_ftran_bfrt(ekk);

        // updateFtran(); computes the pivotal column in the data structure "column"
        self.update_ftran(ekk);

        // updateFtranDSE performs the DSE FTRAN on pi_p
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            let mut row_ep = std::mem::take(&mut self.row_ep);
            self.update_ftran_dse(ekk, &mut row_ep);
            self.row_ep = row_ep;
        }
        ekk.analysis_.simplex_timer_stop(IterateFtranClock, 0);

        // updateVerify() Checks row-wise pivot against column-wise pivot for
        // numerical trouble
        ekk.analysis_.simplex_timer_start(IterateVerifyClock, 0);
        self.update_verify(ekk);
        ekk.analysis_.simplex_timer_stop(IterateVerifyClock, 0);

        // updateDual() Updates the dual values
        ekk.analysis_.simplex_timer_start(IterateDualClock, 0);
        self.update_dual(ekk);
        ekk.analysis_.simplex_timer_stop(IterateDualClock, 0);

        // updatePrimal(&row_ep); Updates the primal values and the edge weights
        ekk.analysis_.simplex_timer_start(IteratePrimalClock, 0);
        let mut row_ep = std::mem::take(&mut self.row_ep);
        self.update_primal(ekk, &mut row_ep);
        self.row_ep = row_ep;
        ekk.analysis_.simplex_timer_stop(IteratePrimalClock, 0);
        // After primal update in dual simplex the primal objective value is not known
        ekk.status_.has_primal_objective_value = false;

        // Update the basis representation
        ekk.analysis_.simplex_timer_start(IteratePivotsClock, 0);
        self.update_pivots(ekk);
        ekk.analysis_.simplex_timer_stop(IteratePivotsClock, 0);

        if self.new_devex_framework {
            // Initialise new Devex framework
            ekk.analysis_.simplex_timer_start(IterateDevexIzClock, 0);
            self.initialise_devex_framework(ekk);
            ekk.analysis_.simplex_timer_stop(IterateDevexIzClock, 0);
        }

        // Analyse the iteration: possibly report; possibly switch strategy
        self.iteration_analysis(ekk);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::iterateTasks (lines 1270-1310, 41 lines)
    /// SIP (`simplex_strategy = kSimplexStrategyDualTasks`): free `spawn` / `sync` at
    /// 1284-1302 (`highs_parallel::{spawn, sync}`).
    pub fn iterate_tasks(&mut self, ekk: &mut HEkk) {
        self.slice_price = 1;

        // Group 1
        self.choose_row(ekk);

        // Disable slice when too sparse
        if 1.0 * self.row_ep.count as f64 * self.inv_solver_num_row < 0.01 {
            self.slice_price = 0;
        }

        ekk.analysis_.simplex_timer_start(Group1Clock, 0);
        // #pragma omp parallel
        // #pragma omp single
        {
            // The C++ lambdas capture `this` (and through it `ekk_instance_`) by reference:
            // raw pointers taken before the spawns (rule (s)).
            let this = TaskPtr(self as *mut HEkkDual);
            let ekk_ptr = TaskPtr(ekk as *mut HEkk);
            // #pragma omp task
            let dse_task = move || {
                // SAFETY: C++ 1284-1287 task protocol: until the sync at 1302 this task
                // alone writes col_DSE and info_.row_dse_density; it reads row_ep (not
                // written by chooseColumn / chooseColumnSlice, which only read it), the
                // factor and the LP. `rebuild_reason` is read here while chooseColumn(Slice)
                // may write it: the C++ race, kept as is (rule (f)).
                let this = unsafe { &mut *this.get() };
                // SAFETY: as above.
                let ekk = unsafe { &mut *ekk_ptr.get() };
                // C++ passes &col_DSE while `this` stays live: alias via a raw pointer.
                let col_dse: *mut HVector = std::ptr::addr_of_mut!(this.col_dse);
                // SAFETY: col_DSE is owned by this task until the sync at 1302 (above).
                unsafe { (*col_dse).copy(&this.row_ep) };
                // SAFETY: as above; updateFtranDSE does not touch col_DSE through `this`.
                this.update_ftran_dse(ekk, unsafe { &mut *col_dse });
            };
            // SAFETY: spawned on this thread's worker deque and joined by the sync at
            // C++ 1302 below, before anything it borrows goes out of scope.
            unsafe { spawn(dse_task) };
            // #pragma omp task
            {
                // C++ passes &row_ep while `this` stays live (and the DSE task reads it):
                // alias via a raw pointer.
                let row_ep: *mut HVector = std::ptr::addr_of_mut!(self.row_ep);
                if self.slice_price != 0 {
                    // SAFETY: chooseColumnSlice reads row_ep through this pointer only and
                    // does not touch `self.row_ep` otherwise.
                    self.choose_column_slice(ekk, unsafe { &mut *row_ep });
                } else {
                    // SAFETY: as above, for chooseColumn.
                    self.choose_column(ekk, unsafe { &mut *row_ep });
                }
                // #pragma omp task
                let bfrt_task = move || {
                    // SAFETY: C++ 1295 task protocol: until the sync at 1299 this task alone
                    // writes col_BFRT, info_.col_bfrt_density and (updateFlip) the flipped
                    // bounds; updateFtran on the spawning thread writes col_aq,
                    // info_.col_aq_density and alpha_col, and the DSE task writes col_DSE.
                    let this = unsafe { &mut *this.get() };
                    // SAFETY: as above.
                    let ekk = unsafe { &mut *ekk_ptr.get() };
                    this.update_ftran_bfrt(ekk);
                };
                // SAFETY: spawned on this thread's worker deque and joined by the sync at
                // C++ 1299 below.
                unsafe { spawn(bfrt_task) };
                // #pragma omp task
                self.update_ftran(ekk);
                // #pragma omp taskwait
                // SAFETY: pairs with the updateFtranBFRT spawn (C++ 1295 / 1299, LIFO).
                unsafe { sync() };
            }

            // SAFETY: pairs with the col_DSE spawn (C++ 1284 / 1302, LIFO).
            unsafe { sync() };
        }
        ekk.analysis_.simplex_timer_stop(Group1Clock, 0);

        self.update_verify(ekk);
        self.update_dual(ekk);
        let mut col_dse = std::mem::take(&mut self.col_dse);
        self.update_primal(ekk, &mut col_dse);
        self.col_dse = col_dse;
        self.update_pivots(ekk);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::iterationAnalysisData (lines 1312-1370, 59 lines)
    pub fn iteration_analysis_data(&mut self, ekk: &mut HEkk) {
        let cost_scale_factor = 2.0f64.powf(-(ekk.options_.cost_scale_factor as f64));
        let info = &ekk.info_;
        let analysis = &mut ekk.analysis_;
        analysis.simplex_strategy = info.simplex_strategy;
        analysis.edge_weight_mode = self.edge_weight_mode;
        analysis.solve_phase = self.solve_phase;
        analysis.simplex_iteration_count = ekk.iteration_count_;
        analysis.devex_iteration_count = self.num_devex_iterations;
        analysis.pivotal_row_index = self.row_out;
        analysis.leaving_variable = self.variable_out;
        analysis.entering_variable = self.variable_in;
        analysis.rebuild_reason = self.rebuild_reason;
        analysis.reduced_rhs_value = 0.0;
        analysis.reduced_cost_value = 0.0;
        analysis.edge_weight = 0.0;
        analysis.primal_delta = self.delta_primal;
        analysis.primal_step = self.theta_primal;
        analysis.dual_step = self.theta_dual * cost_scale_factor;
        analysis.pivot_value_from_column = self.alpha_col;
        analysis.pivot_value_from_row = self.alpha_row;
        analysis.factor_pivot_threshold = info.factor_pivot_threshold;
        analysis.numerical_trouble = self.numerical_trouble;
        analysis.edge_weight_error = ekk.edge_weight_error_;
        analysis.objective_value = info.updated_dual_objective_value;
        // Since maximization is achieved by minimizing the LP with negated
        // costs, in phase 2 the dual objective value is negated, so flip
        // its sign according to the LP sense
        if self.solve_phase == kSolvePhase2 {
            analysis.objective_value *= (ekk.lp_.sense_ as HighsInt) as f64;
        }
        analysis.num_primal_infeasibility = info.num_primal_infeasibilities;
        analysis.sum_primal_infeasibility = info.sum_primal_infeasibilities;
        if self.solve_phase == kSolvePhase1 {
            analysis.num_dual_infeasibility = analysis.num_dual_phase_1_lp_dual_infeasibility;
            analysis.sum_dual_infeasibility = analysis.sum_dual_phase_1_lp_dual_infeasibility;
        } else {
            analysis.num_dual_infeasibility = info.num_dual_infeasibilities;
            analysis.sum_dual_infeasibility = info.sum_dual_infeasibilities;
        }
        if (self.edge_weight_mode == EdgeWeightMode::kDevex) && (self.num_devex_iterations == 0) {
            analysis.num_devex_framework += 1;
        }
        analysis.col_aq_density = info.col_aq_density;
        analysis.row_ep_density = info.row_ep_density;
        analysis.row_ap_density = info.row_ap_density;
        analysis.row_DSE_density = info.row_dse_density;
        analysis.col_basic_feasibility_change_density = info.col_basic_feasibility_change_density;
        analysis.row_basic_feasibility_change_density = info.row_basic_feasibility_change_density;
        analysis.col_BFRT_density = info.col_bfrt_density;
        analysis.primal_col_density = info.primal_col_density;
        analysis.dual_col_density = info.dual_col_density;
        analysis.num_costly_DSE_iteration = info.num_costly_dse_iteration;
        analysis.costly_DSE_measure = info.costly_dse_measure;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::iterationAnalysis (lines 1372-1395, 24 lines)
    pub fn iteration_analysis(&mut self, ekk: &mut HEkk) {
        // Compute the infeasibility data (expensive) if analysing run-time
        // data and the log level is at least kIterationReportLogType
        // (Verbose)
        let make_iteration_report = ekk.analysis_.analyse_simplex_runtime_data
            && ekk.options_.log_dev_level >= (HighsLogType::kVerbose as HighsInt);
        if make_iteration_report {
            ekk.compute_infeasibilities_for_reporting(SimplexAlgorithm::kDual, self.solve_phase);
        }
        // Possibly report on the iteration
        self.iteration_analysis_data(ekk);
        ekk.analysis_.iteration_report();

        // Possibly switch from DSE to Devex
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            let switch_to_devex = ekk.switch_to_devex();
            if switch_to_devex {
                self.edge_weight_mode = EdgeWeightMode::kDevex;
                self.initialise_devex_framework(ekk);
            }
        }
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.iteration_record();
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::reportRebuild (lines 1397-1405, 9 lines)
    pub fn report_rebuild(&mut self, ekk: &mut HEkk, reason_for_rebuild: HighsInt) {
        ekk.analysis_.simplex_timer_start(ReportRebuildClock, 0);
        self.iteration_analysis_data(ekk);
        ekk.analysis_.rebuild_reason = reason_for_rebuild;
        ekk.analysis_.rebuild_reason_string = ekk.rebuild_reason(reason_for_rebuild);
        if ekk.options_.output_flag {
            ekk.analysis_.invert_report();
        }
        ekk.analysis_.simplex_timer_stop(ReportRebuildClock, 0);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::chooseRow (lines 1407-1502, 96 lines)
    pub fn choose_row(&mut self, ekk: &mut HEkk) {
        // Choose the index of a row to leave the basis (CHUZR)
        //
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }
        // Zero the infeasibility of any taboo rows
        ekk.apply_taboo_row_out(&mut self.dual_rhs.work_infeasibility, 0.0);
        // Choose candidates repeatedly until candidate is OK or optimality is
        // detected
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            let return_status = ekk.dev_debug_dual_steepest_edge_weights("chooseRow");
            debug_assert!(
                return_status == HighsDebugStatus::kNotChecked
                    || return_status == HighsDebugStatus::kOk
            );
        }
        loop {
            // Choose the index of a good row to leave the basis
            self.dual_rhs.choose_normal(ekk, &mut self.row_out);
            if self.row_out == kNoRowChosen {
                // No index found so may be dual optimal.
                self.rebuild_reason = kRebuildReasonPossiblyOptimal;
                return;
            }
            // Compute pi_p = B^{-T}e_p in row_ep
            ekk.analysis_.simplex_timer_start(BtranClock, 0);
            // Set up RHS for BTRAN
            self.row_ep.clear();
            self.row_ep.count = 1;
            self.row_ep.index[0] = self.row_out;
            self.row_ep.array[self.row_out as usize] = 1.0;
            self.row_ep.pack_flag = true;
            if ekk.analysis_.analyse_simplex_summary_data {
                ekk.analysis_.operation_record_before(
                    kSimplexNlaBtranEp,
                    &self.row_ep,
                    ekk.info_.row_ep_density,
                );
            }
            // Perform BTRAN
            {
                let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
                ekk.simplex_nla_
                    .btran(&refs, &mut self.row_ep, ekk.info_.row_ep_density, None);
            }
            if ekk.analysis_.analyse_simplex_summary_data {
                ekk.analysis_
                    .operation_record_after(kSimplexNlaBtranEp, &self.row_ep);
            }
            ekk.analysis_.simplex_timer_stop(BtranClock, 0);
            // Verify DSE weight
            if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                // For DSE, see how accurate the updated weight is
                // Save the updated weight
                let updated_edge_weight = ekk.dual_edge_weight_[self.row_out as usize];
                // Compute the weight from row_ep and over-write the updated weight
                if ekk.simplex_in_scaled_space_ {
                    let w = self.row_ep.norm2();
                    ekk.dual_edge_weight_[self.row_out as usize] = w;
                    self.computed_edge_weight = w;
                } else {
                    let w = {
                        let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
                        ekk.simplex_nla_.row_ep2_norm_in_scaled_space(
                            &refs,
                            self.row_out,
                            &self.row_ep,
                        )
                    };
                    ekk.dual_edge_weight_[self.row_out as usize] = w;
                    self.computed_edge_weight = w;
                }
                // If the weight error is acceptable then break out of the
                // loop. All we worry about is accepting rows with weights
                // which are not too small, since this can make the row look
                // unreasonably attractive
                if self.accept_dual_steepest_edge_weight(ekk, updated_edge_weight) {
                    break;
                }
                // Weight error is unacceptable so look for another
                // candidate. Of course, it's possible that the same
                // candidate is chosen, but the weight will be correct (so
                // no infinite loop).
            } else {
                // If not using DSE then accept the row by breaking out of
                // the loop
                break;
            }
        }
        // Recover the infeasibility of any taboo rows
        ekk.unapply_taboo_row_out(&mut self.dual_rhs.work_infeasibility);

        // Index of row to leave the basis has been found
        //
        // Assign basic info:
        //
        // Record the column (variable) associated with the leaving row
        self.variable_out = ekk.basis_.basic_index_[self.row_out as usize];
        // Record the change in primal variable associated with the move to the bound
        // being violated
        let row_out = self.row_out as usize;
        if ekk.info_.base_value_[row_out] < ekk.info_.base_lower_[row_out] {
            // Below the lower bound so set delta_primal = value - LB < 0
            self.delta_primal = ekk.info_.base_value_[row_out] - ekk.info_.base_lower_[row_out];
        } else {
            // Above the upper bound so set delta_primal = value - UB > 0
            self.delta_primal = ekk.info_.base_value_[row_out] - ekk.info_.base_upper_[row_out];
        }
        // Set move_out to be -1 if delta_primal<0, otherwise +1 (since
        // delta_primal>0)
        self.move_out = if self.delta_primal < 0.0 { -1 } else { 1 };
        // Update the record of average row_ep (pi_p) density. This ignores
        // any BTRANs done for skipped candidates
        let local_row_ep_density = (self.row_ep.count as f64) * self.inv_solver_num_row;
        let mut density = ekk.info_.row_ep_density;
        ekk.update_operation_result_density(local_row_ep_density, &mut density);
        ekk.info_.row_ep_density = density;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::acceptDualSteepestEdgeWeight (lines 1504-1514, 11 lines)
    pub fn accept_dual_steepest_edge_weight(
        &mut self,
        ekk: &mut HEkk,
        updated_edge_weight: f64,
    ) -> bool {
        // Accept the updated weight if it is at least a quarter of the
        // computed weight. Excessively large updated weights don't matter!
        let accept_weight =
            updated_edge_weight >= kAcceptDseWeightThreshold * self.computed_edge_weight;
        ekk.assess_dse_weight_error(self.computed_edge_weight, updated_edge_weight);
        ekk.analysis_
            .dual_steepest_edge_weight_error(self.computed_edge_weight, updated_edge_weight);
        accept_weight
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::newDevexFramework (lines 1516-1544, 29 lines)
    pub fn new_devex_framework(&mut self, ekk: &mut HEkk, updated_edge_weight: f64) -> bool {
        // Analyse the Devex weight to determine whether a new framework
        // should be set up
        //
        // There is a new Devex framework if either
        //
        // 1) The weight inaccuracy ratio exceeds kMaxAllowedDevexWeightRatio
        //
        // 2) There have been max(kMinAbsNumberDevexIterations,
        // numRow/kMinRlvNumberDevexIterations) Devex iterations
        //
        let kMinAbsNumberDevexIterations: HighsInt = 25;
        let kMinRlvNumberDevexIterations: f64 = 1e-2;
        let kMaxAllowedDevexWeightRatio: f64 = 3.0;

        let devex_ratio = cpp_max(
            updated_edge_weight / self.computed_edge_weight,
            self.computed_edge_weight / updated_edge_weight,
        );
        let mut i_te: HighsInt =
            (self.solver_num_row as f64 / kMinRlvNumberDevexIterations) as HighsInt;
        i_te = cpp_max(kMinAbsNumberDevexIterations, i_te);
        // Square kMaxAllowedDevexWeightRatio due to keeping squared
        // weights
        let accept_ratio_threshold = kMaxAllowedDevexWeightRatio * kMaxAllowedDevexWeightRatio;
        let accept_ratio = devex_ratio <= accept_ratio_threshold;
        let accept_it = self.num_devex_iterations <= i_te;
        let return_new_devex_framework: bool = !accept_ratio || !accept_it;
        return_new_devex_framework
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::chooseColumn (lines 1546-1735, 190 lines)
    pub fn choose_column(&mut self, ekk: &mut HEkk, row_ep: &mut HVector) {
        // Compute pivot row (PRICE) and choose the index of a column to enter the
        // basis (CHUZC)
        //
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }

        let mut debug_price_report: HighsInt = kDebugReportOff;
        let debug_price_report_on = false;
        let debug_small_pivot_issue_report_on = false;
        let mut debug_small_pivot_issue_report = false;
        let mut debug_rows_report = false;
        if ekk.debug_iteration_report_ {
            if debug_price_report_on {
                debug_price_report = kDebugReportAll;
            }
            debug_rows_report = debug_price_report_on;
            debug_small_pivot_issue_report = debug_small_pivot_issue_report_on;
            if debug_price_report != kDebugReportOff
                || debug_rows_report
                || debug_small_pivot_issue_report
            {
                println!(
                    "HEkkDual::chooseColumn Check iter = {}",
                    ekk.iteration_count_
                );
            }
        }
        //
        // PRICE
        //
        let quad_precision = false;
        ekk.tableau_row_price(quad_precision, row_ep, &mut self.row_ap, debug_price_report);
        if debug_rows_report {
            ekk.simplex_nla_
                .report_array_offset("Row a_p", 0, &self.row_ap, true);
            ekk.simplex_nla_
                .report_array_offset("Row e_p", ekk.lp_.num_col_, row_ep, true);
        }
        //
        // CHUZC
        //
        // Section 0: Clear data and call createFreemove to set a value of
        // nonbasicMove for all free columns to prevent their dual values
        // from being changed.
        ekk.analysis_.simplex_timer_start(Chuzc0Clock, 0);
        self.dual_row.clear(ekk);
        self.dual_row.work_delta = self.delta_primal;
        self.dual_row.create_freemove(ekk, row_ep);
        ekk.analysis_.simplex_timer_stop(Chuzc0Clock, 0);
        //
        // Section 1: Pack row_ap and row_ep
        ekk.analysis_.simplex_timer_start(Chuzc1Clock, 0);
        // Pack row_ap into the packIndex/Value of HEkkDualRow
        self.dual_row.choose_makepack(ekk, &self.row_ap, 0);
        // Pack row_ep into the packIndex/Value of HEkkDualRow
        self.dual_row
            .choose_makepack(ekk, row_ep, self.solver_num_col);
        let row_ep_scale = ekk.get_value_scale(self.dual_row.pack_count, &self.dual_row.pack_value);
        ekk.analysis_.simplex_timer_stop(Chuzc1Clock, 0);
        // Loop until an acceptable pivot is found. Each pass either finds a
        // pivot, identifies possible unboundedness, or reduced the number
        // of nonzeros in dualRow.pack_value
        let mut chuzc_pass: HighsInt = 0;
        loop {
            //
            // Section 2: Determine the possible variables - candidates for CHUZC
            ekk.analysis_.simplex_timer_start(Chuzc2Clock, 0);
            self.dual_row.choose_possible(ekk);
            ekk.analysis_.simplex_timer_stop(Chuzc2Clock, 0);
            //
            // Take action if the step to an expanded bound is not positive, or
            // there are no candidates for CHUZC
            self.variable_in = -1;
            if self.dual_row.work_theta <= 0.0 || self.dual_row.work_count == 0 {
                if chuzc_pass > 0 && debug_small_pivot_issue_report {
                    println!(
                        "                                                       Negative step or no candidates after {:2} CHUZC passes",
                        chuzc_pass
                    );
                }
                if debug_rows_report {
                    ekk.simplex_nla_.report_vector(
                        "dualRow.packValue/Index",
                        self.dual_row.pack_count,
                        &self.dual_row.pack_value,
                        &self.dual_row.pack_index,
                        true,
                    );
                }
                self.rebuild_reason = kRebuildReasonPossiblyDualUnbounded;
                return;
            }
            //
            // Sections 3 and 4: Perform (bound-flipping) ratio test. This can
            // fail if the dual values are excessively large
            let choose_column_fail = self.dual_row.choose_final(ekk) != 0;
            if choose_column_fail {
                self.rebuild_reason = kRebuildReasonChooseColumnFail;
                return;
            }
            if self.dual_row.work_pivot >= 0 {
                let growth_tolerance = ekk.options_.dual_simplex_pivot_growth_tolerance; // kHighsMacheps; //
                // A pivot has been chosen
                let alpha_row = self.dual_row.work_alpha;
                debug_assert!(alpha_row != 0.0);
                let scaled_value = row_ep_scale * alpha_row;
                if scaled_value.abs() <= growth_tolerance {
                    if chuzc_pass == 0 && debug_small_pivot_issue_report {
                        let wp = self.dual_row.work_pivot as usize;
                        print!(
                            "CHUZC: Solve {:6}; Iter {:4}; ||e_p|| = {:11.4e}: Variable {:6} Pivot {:11.4e} (dual {:11.4e}; ratio = {:11.4e}) is small",
                            ekk.debug_solve_call_num_,
                            ekk.iteration_count_,
                            1.0 / row_ep_scale,
                            self.dual_row.work_pivot,
                            self.dual_row.work_alpha,
                            ekk.info_.work_dual_[wp],
                            ekk.info_.work_dual_[wp] / self.dual_row.work_alpha
                        );
                    }
                    // On the first pass, try to make the pivotal row more accurate
                    if chuzc_pass == 0 {
                        if debug_small_pivot_issue_report {
                            println!(": improve row");
                        }
                        ekk.analysis_.num_improve_choose_column_row_call += 1;
                        self.improve_choose_column_row(ekk, row_ep);
                    } else {
                        // Remove the pivot
                        ekk.analysis_.num_remove_pivot_from_pack += 1;
                        for i in 0..self.dual_row.pack_count {
                            let iu = i as usize;
                            if self.dual_row.pack_index[iu] == self.dual_row.work_pivot {
                                let last = (self.dual_row.pack_count - 1) as usize;
                                self.dual_row.pack_index[iu] = self.dual_row.pack_index[last];
                                self.dual_row.pack_value[iu] = self.dual_row.pack_value[last];
                                self.dual_row.pack_count -= 1;
                                if chuzc_pass == 0 && debug_small_pivot_issue_report {
                                    print!(
                                        ": removing pivot gives pack count = {:6}",
                                        self.dual_row.pack_count
                                    );
                                }
                                break;
                            }
                        }
                        if chuzc_pass == 0 && debug_small_pivot_issue_report {
                            println!(
                                " so {}",
                                if self.dual_row.pack_count > 0 {
                                    "repeat CHUZC"
                                } else {
                                    "consider unbounded"
                                }
                            );
                        }
                    }
                    // Indicate that no pivot has been chosen
                    self.dual_row.work_pivot = -1;
                } else if chuzc_pass > 0 && debug_small_pivot_issue_report {
                    let wp = self.dual_row.work_pivot as usize;
                    println!(
                        "                                                       Variable {:6} Pivot {:11.4e} (dual {:11.4e}; ratio = {:11.4e}) is OK after {:2} CHUZC passes",
                        self.dual_row.work_pivot,
                        self.dual_row.work_alpha,
                        ekk.info_.work_dual_[wp],
                        ekk.info_.work_dual_[wp] / self.dual_row.work_alpha,
                        chuzc_pass
                    );
                }
            } else {
                // No pivot has been chosen
                debug_assert!(self.dual_row.work_pivot == -1);
                if chuzc_pass > 0 && debug_small_pivot_issue_report {
                    println!(
                        "                                                       No pivot after {:2} CHUZC passes",
                        chuzc_pass
                    );
                }
                break;
            }
            // If a pivot has been chosen, or there are no more packed values
            // then end CHUZC
            if self.dual_row.work_pivot >= 0 || self.dual_row.pack_count <= 0 {
                break;
            }
            chuzc_pass += 1;
        }
        //
        // Section 5: Reset the nonbasicMove values for free columns
        ekk.analysis_.simplex_timer_start(Chuzc5Clock, 0);
        self.dual_row.delete_freemove(ekk);
        ekk.analysis_.simplex_timer_stop(Chuzc5Clock, 0);
        // Record values for basis change, checking for numerical problems and update
        // of dual variables
        self.variable_in = self.dual_row.work_pivot; // Index of the column entering the basis
        self.alpha_row = self.dual_row.work_alpha; // Pivot value computed row-wise - used for
        // numerical checking
        self.theta_dual = self.dual_row.work_theta; // Dual step length

        if self.edge_weight_mode == EdgeWeightMode::kDevex && !self.new_devex_framework {
            // When using Devex, unless a new framework is to be used, get the
            // exact weight for the pivotal row and, based on its accuracy,
            // determine that a new framework is to be used. In serial
            // new_devex_framework should only ever be false at this point in
            // this method, but in PAMI, this method may be called multiple
            // times in minor iterations and the new framework is set up in
            // majorUpdate.
            ekk.analysis_.simplex_timer_start(DevexWtClock, 0);
            // Determine the exact Devex weight
            self.dual_row.compute_devex_weight(ekk, -1);
            self.computed_edge_weight = self.dual_row.computed_edge_weight;
            self.computed_edge_weight = cpp_max(1.0, self.computed_edge_weight);
            ekk.analysis_.simplex_timer_stop(DevexWtClock, 0);
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::improveChooseColumnRow (lines 1737-1783, 47 lines)
    pub fn improve_choose_column_row(&mut self, ekk: &mut HEkk, row_ep: &mut HVector) {
        let debug_price_report_on = false;
        let mut debug_price_report: HighsInt = kDebugReportOff;
        let mut debug_rows_report = false;
        if ekk.debug_iteration_report_ {
            if debug_price_report_on {
                debug_price_report = kDebugReportAll;
            }
            debug_rows_report = debug_price_report_on;
            if debug_price_report != kDebugReportOff {
                println!(
                    "HEkkDual::chooseColumn Check iter = {}",
                    ekk.iteration_count_
                );
            }
        }

        ekk.analysis_.simplex_timer_start(Chuzc5Clock, 0);
        self.dual_row.delete_freemove(ekk);
        ekk.analysis_.simplex_timer_stop(Chuzc5Clock, 0);

        if debug_rows_report {
            ekk.simplex_nla_
                .report_array_offset("Row e_p.0", ekk.lp_.num_col_, row_ep, true);
        }
        ekk.unit_btran_iterative_refinement(self.row_out, row_ep);
        if debug_rows_report {
            ekk.simplex_nla_
                .report_array_offset("Row e_p.1", ekk.lp_.num_col_, row_ep, true);
        }

        let quad_precision = true;
        ekk.tableau_row_price(quad_precision, row_ep, &mut self.row_ap, debug_price_report);
        if debug_rows_report {
            ekk.simplex_nla_
                .report_array_offset("Row a_p", 0, &self.row_ap, true);
            ekk.simplex_nla_
                .report_array_offset("Row e_p", ekk.lp_.num_col_, row_ep, true);
        }
        ekk.analysis_.simplex_timer_start(Chuzc0Clock, 0);
        self.dual_row.clear(ekk);
        self.dual_row.work_delta = self.delta_primal;
        self.dual_row.create_freemove(ekk, row_ep);
        ekk.analysis_.simplex_timer_stop(Chuzc0Clock, 0);
        //
        // Section 1: Pack row_ap and row_ep
        ekk.analysis_.simplex_timer_start(Chuzc1Clock, 0);
        // Pack row_ap into the packIndex/Value of HEkkDualRow
        self.dual_row.choose_makepack(ekk, &self.row_ap, 0);
        // Pack row_ep into the packIndex/Value of HEkkDualRow
        self.dual_row
            .choose_makepack(ekk, row_ep, self.solver_num_col);
        ekk.analysis_.simplex_timer_stop(Chuzc1Clock, 0);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::chooseColumnSlice (lines 1785-1940, 156 lines)
    /// Sliced PRICE + CHUZC1 (`spawn` per slice, 1824-1857), used by SIP and PAMI.
    pub fn choose_column_slice(&mut self, ekk: &mut HEkk, row_ep: &mut HVector) {
        // Choose the index of a column to enter the basis (CHUZC) by
        // exploiting slices of the pivotal row - for SIP and PAMI
        //
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }

        ekk.analysis_.simplex_timer_start(Chuzc0Clock, 0);
        self.dual_row.clear(ekk);
        self.dual_row.work_delta = self.delta_primal;
        self.dual_row.create_freemove(ekk, row_ep);
        ekk.analysis_.simplex_timer_stop(Chuzc0Clock, 0);

        //  const HighsInt solver_num_row = ekk_instance_.lp_.num_row_;
        let local_density = 1.0 * row_ep.count as f64 * self.inv_solver_num_row;
        let mut use_col_price: bool = false;
        let mut use_row_price_w_switch: bool = false;
        let info = &ekk.info_;
        ekk.choose_price_technique(
            info.price_strategy,
            local_density,
            &mut use_col_price,
            &mut use_row_price_w_switch,
        );

        if ekk.analysis_.analyse_simplex_summary_data {
            let row_ep_count = row_ep.count;
            if use_col_price {
                ekk.analysis_
                    .operation_record_before_count(kSimplexNlaPriceAp, row_ep_count, 0.0);
                ekk.analysis_.num_col_price += 1;
            } else if use_row_price_w_switch {
                ekk.analysis_.operation_record_before_count(
                    kSimplexNlaPriceAp,
                    row_ep_count,
                    ekk.info_.row_ep_density,
                );
                ekk.analysis_.num_row_price_with_switch += 1;
            } else {
                ekk.analysis_.operation_record_before_count(
                    kSimplexNlaPriceAp,
                    row_ep_count,
                    ekk.info_.row_ep_density,
                );
                ekk.analysis_.num_row_price += 1;
            }
        }
        ekk.analysis_.simplex_timer_start(PriceChuzc1Clock, 0);
        // Row_ep:         PACK + CC1

        // The C++ lambdas capture `this`, `row_ep` and the locals by reference: field-level
        // raw pointers taken before the spawn (rule (s)).
        let dual_row_ptr = TaskPtr(std::ptr::addr_of_mut!(self.dual_row));
        let solver_num_col_ptr = TaskPtr(std::ptr::addr_of_mut!(self.solver_num_col));
        let row_ep_ptr = TaskPtr(row_ep as *mut HVector);
        let ekk_ptr = TaskPtr(ekk as *mut HEkk);
        let pack_task = move || {
            // SAFETY: C++ 1824-1827 task protocol: until the sync at 1857 this task alone
            // touches dualRow; row_ep, solver_num_col and the ekk data choosePossible reads
            // (info_.update_count, options_, basis_.nonbasic_move_, info_.workDual_) are not
            // written by the for_each slices.
            let dual_row = unsafe { &mut *dual_row_ptr.get() };
            // SAFETY: as above (read only).
            let row_ep = unsafe { &*row_ep_ptr.get() };
            // SAFETY: as above (read only).
            let solver_num_col = unsafe { *solver_num_col_ptr.get() };
            // SAFETY: as above.
            let ekk = unsafe { &mut *ekk_ptr.get() };
            dual_row.choose_makepack(ekk, row_ep, solver_num_col);
            dual_row.choose_possible(ekk);
        };
        // SAFETY: spawned on this thread's worker deque and joined by the sync at C++ 1857
        // below, before anything it borrows goes out of scope.
        unsafe { spawn(pack_task) };

        // Row_ap: PRICE + PACK + CC1
        let slice_row_ap_ptr = TaskPtr(self.slice_row_ap.as_mut_ptr());
        let slice_dual_row_ptr = TaskPtr(self.slice_dual_row.as_mut_ptr());
        let slice_a_matrix_ptr = TaskPtr(self.slice_a_matrix.as_mut_ptr());
        let slice_ar_matrix_ptr = TaskPtr(self.slice_ar_matrix.as_mut_ptr());
        let slice_start_ptr = TaskPtr(self.slice_start.as_mut_ptr());
        let delta_primal_ptr = TaskPtr(std::ptr::addr_of_mut!(self.delta_primal));
        let use_col_price_ptr = TaskPtr(std::ptr::addr_of_mut!(use_col_price));
        let use_row_price_w_switch_ptr = TaskPtr(std::ptr::addr_of_mut!(use_row_price_w_switch));
        for_each(
            0,
            self.slice_num,
            &move |start: HighsInt, end: HighsInt| {
                // SAFETY (whole body): C++ 1830-1855 for_each protocol: index i alone writes
                // slice_row_ap[i] and slice_dualRow[i]; slice_a_matrix, slice_ar_matrix,
                // slice_start, delta_primal, the two PRICE flags, row_ep and the ekk data the
                // kernels read (info_.row_ap_density, info_.update_count, options_,
                // basis_.nonbasic_move_, info_.workDual_) are not written until the sync.
                let use_col_price = unsafe { *use_col_price_ptr.get() };
                let use_row_price_w_switch = unsafe { *use_row_price_w_switch_ptr.get() };
                let row_ep = unsafe { &*row_ep_ptr.get() };
                let quad_precision = false;
                for i in start..end {
                    let iu = i as usize;
                    let slice_row_ap = unsafe { &mut *slice_row_ap_ptr.get().add(iu) };
                    slice_row_ap.clear();

                    if use_col_price {
                        // Perform column-wise PRICE
                        let slice_a_matrix = unsafe { &*slice_a_matrix_ptr.get().add(iu) };
                        slice_a_matrix.price_by_column(
                            quad_precision,
                            slice_row_ap,
                            row_ep,
                            kDebugReportOff,
                        );
                    } else if use_row_price_w_switch {
                        // Perform hyper-sparse row-wise PRICE, but switch if the density of
                        // row_ap becomes extreme
                        let slice_ar_matrix = unsafe { &*slice_ar_matrix_ptr.get().add(iu) };
                        let row_ap_density = unsafe { (*ekk_ptr.get()).info_.row_ap_density };
                        slice_ar_matrix.price_by_row_with_switch(
                            quad_precision,
                            slice_row_ap,
                            row_ep,
                            row_ap_density,
                            0,
                            kHyperPriceDensity,
                            kDebugReportOff,
                        );
                    } else {
                        // Perform hyper-sparse row-wise PRICE
                        let slice_ar_matrix = unsafe { &*slice_ar_matrix_ptr.get().add(iu) };
                        slice_ar_matrix.price_by_row(
                            quad_precision,
                            slice_row_ap,
                            row_ep,
                            kDebugReportOff,
                        );
                    }

                    let slice_dual_row = unsafe { &mut *slice_dual_row_ptr.get().add(iu) };
                    let ekk = unsafe { &mut *ekk_ptr.get() };
                    slice_dual_row.clear(ekk);
                    slice_dual_row.work_delta = unsafe { *delta_primal_ptr.get() };
                    slice_dual_row.choose_makepack(ekk, slice_row_ap, unsafe {
                        *slice_start_ptr.get().add(iu)
                    });
                    slice_dual_row.choose_possible(ekk);
                }
            },
            1,
        );

        // SAFETY: pairs with the chooseMakepack / choosePossible spawn (C++ 1824 / 1857).
        unsafe { sync() };

        if ekk.analysis_.analyse_simplex_summary_data {
            // Determine the nonzero count of the whole row
            let mut row_ap_count: HighsInt = 0;
            for i in 0..self.slice_num {
                row_ap_count += self.slice_row_ap[i as usize].count;
            }
            ekk.analysis_
                .operation_record_after_count(kSimplexNlaPriceAp, row_ap_count);
        }

        // Join CC1 results here
        for i in 0..self.slice_num {
            self.dual_row
                .choose_joinpack(ekk, &self.slice_dual_row[i as usize]);
        }

        ekk.analysis_.simplex_timer_stop(PriceChuzc1Clock, 0);

        // Infeasible we created before
        self.variable_in = -1;
        if self.dual_row.work_theta <= 0.0 || self.dual_row.work_count == 0 {
            self.rebuild_reason = kRebuildReasonPossiblyDualUnbounded;
            return;
        }

        // Choose column 2, This only happens if didn't go out
        let return_code = self.dual_row.choose_final(ekk);
        if return_code != 0 {
            // Only returns -1, if not zero
            debug_assert!(return_code == -1);
            if return_code < 0 {
                self.rebuild_reason = kRebuildReasonChooseColumnFail;
            } else {
                self.rebuild_reason = kRebuildReasonPossiblyDualUnbounded;
            }
            return;
        }

        if self.slice_num == 0 {
            // This check is only done for serial code - since packIndex/Value
            // is distributed
            let num_infeasibility = self.dual_row.debug_choose_column_infeasibilities(ekk);
            if num_infeasibility != 0 {
                highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kError,
                    "chooseFinal would create {} dual infeasibilities\n",
                    num_infeasibility
                );
                ekk.analysis_.simplex_timer_stop(Chuzc4dClock, 0);
                self.rebuild_reason = kRebuildReasonChooseColumnFail;
                return;
            }
        }

        ekk.analysis_.simplex_timer_start(Chuzc5Clock, 0);
        self.dual_row.delete_freemove(ekk);
        ekk.analysis_.simplex_timer_stop(Chuzc5Clock, 0);

        self.variable_in = self.dual_row.work_pivot;
        self.alpha_row = self.dual_row.work_alpha;
        self.theta_dual = self.dual_row.work_theta;

        if self.edge_weight_mode == EdgeWeightMode::kDevex && !self.new_devex_framework {
            // When using Devex, unless a new framework is to be used, get the
            // exact weight for the pivotal row and, based on its accuracy,
            // determine that a new framework is to be used. In serial
            // new_devex_framework should only ever be false at this point in
            // this method, but in PAMI, this method may be called multiple
            // times in minor iterations and the new framework is set up in
            // majorUpdate.
            ekk.analysis_.simplex_timer_start(DevexWtClock, 0);
            // Determine the partial sums of the exact Devex weight
            // First the partial sum for row_ep
            self.dual_row.compute_devex_weight(ekk, -1);
            // Second the partial sums for the slices of row_ap
            for i in 0..self.slice_num {
                self.slice_dual_row[i as usize].compute_devex_weight(ekk, i);
            }
            // Accumulate the partial sums
            // Initialise with the partial sum for row_ep
            self.computed_edge_weight = self.dual_row.computed_edge_weight;
            // Update with the partial sum for row_ep
            for i in 0..self.slice_num {
                self.computed_edge_weight += self.slice_dual_row[i as usize].computed_edge_weight;
            }
            self.computed_edge_weight = cpp_max(1.0, self.computed_edge_weight);
            ekk.analysis_.simplex_timer_stop(DevexWtClock, 0);
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests
    use crate::lp_data::highs_run::Highs;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;
    use crate::prelude::*;

    /// Deterministic LP `min c'x s.t. A x >= b, 0 <= x <= u`: the slack basis is dual
    /// feasible and primal infeasible, so dual phase 2 runs `iterate` / `iterateTasks`.
    /// `dense` = false gives a two-band matrix whose `row_ep` stays sparse.
    fn test_lp(num_row: HighsInt, num_col: HighsInt, dense: bool) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = num_col;
        lp.num_row_ = num_row;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = num_col;
        lp.a_matrix_.num_row_ = num_row;
        lp.a_matrix_.start_ = vec![0];
        let mut seed: u64 = 12345;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as HighsInt
        };
        for j in 0..num_col {
            lp.col_cost_.push((1 + next() % 9) as f64);
            lp.col_lower_.push(0.0);
            lp.col_upper_.push(if j % 3 == 0 {
                2.0 + (next() % 3) as f64
            } else {
                kHighsInf
            });
            for i in 0..num_row {
                let nz = if dense {
                    next() % 2 == 0
                } else {
                    i == j % num_row || (i + 1) % num_row == j % num_row
                };
                if nz {
                    lp.a_matrix_.index_.push(i);
                    lp.a_matrix_.value_.push((1 + next() % 9) as f64);
                }
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        for _ in 0..num_row {
            lp.row_lower_.push((10 + next() % 11) as f64);
            lp.row_upper_.push(kHighsInf);
        }
        lp
    }

    /// The bits of one solve: status, iteration count, objective, primal and dual values.
    #[derive(Debug, PartialEq)]
    struct Run {
        status: HighsModelStatus,
        iterations: HighsInt,
        strategy: HighsInt,
        objective: u64,
        col_value: Vec<u64>,
        row_dual: Vec<u64>,
    }

    fn solve(lp: &HighsLp, strategy: HighsInt, concurrency: HighsInt) -> Run {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", "simplex"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("presolve", "off"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_int("simplex_strategy", strategy),
            HighsStatus::kOk
        );
        // `simplex_max_concurrency` has no entry in the port's option setter.
        highs.options_.simplex_max_concurrency = concurrency;
        assert_eq!(highs.pass_model(lp.clone()), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        Run {
            status: highs.model_status_,
            iterations: highs.info_.simplex_iteration_count,
            strategy: highs.ekk_instance_.info_.simplex_strategy,
            objective: highs.info_.objective_function_value.to_bits(),
            col_value: highs
                .solution_
                .col_value
                .iter()
                .map(|v| v.to_bits())
                .collect(),
            row_dual: highs
                .solution_
                .row_dual
                .iter()
                .map(|v| v.to_bits())
                .collect(),
        }
    }

    /// Solve with SIP on a fresh thread owning a `threads`-worker scheduler.
    fn solve_sip(lp: &HighsLp, threads: i32, concurrency: HighsInt) -> Run {
        let lp = lp.clone();
        std::thread::spawn(move || {
            initialize_scheduler(threads);
            let run = solve(&lp, kSimplexStrategyDualTasks, concurrency);
            HighsTaskExecutor::shutdown(true);
            run
        })
        .join()
        .unwrap()
    }

    fn check_sip_matches_serial(lp: &HighsLp) {
        let serial = solve(lp, kSimplexStrategyDual, 1);
        assert_eq!(serial.status, HighsModelStatus::kOptimal);
        assert!(serial.iterations > 0);
        let serial_objective = f64::from_bits(serial.objective);
        for concurrency in [3, 4, 8] {
            // Threads = 1 is the reference: the tasks run in the C++ LIFO order on the
            // main worker.
            let reference = solve_sip(lp, 1, concurrency);
            assert_eq!(reference.status, HighsModelStatus::kOptimal);
            assert_eq!(reference.strategy, kSimplexStrategyDualTasks);
            let objective = f64::from_bits(reference.objective);
            assert!(
                (objective - serial_objective).abs() <= 1e-9 * (1.0 + serial_objective.abs()),
                "SIP objective {objective} != serial {serial_objective}"
            );
            for threads in [1, 2, 4, 8] {
                for _ in 0..2 {
                    assert_eq!(
                        solve_sip(lp, threads, concurrency),
                        reference,
                        "threads = {threads}, concurrency = {concurrency}"
                    );
                }
            }
        }
    }

    #[test]
    fn w5_sip_tasks_iterate_tasks_matches_serial() {
        // `iterateTasks` / `chooseColumnSlice` (HEkkDual.cpp 1270-1310, 1785-1940): a dense
        // row_ep keeps slice_PRICE on (1277), so PRICE / CHUZC1 run on the slices and are
        // joined in slice order (1867-1869). Same optimum as the serial `iterate`, and the
        // same bits at threads 1, 2, 4, 8, run twice each.
        check_sip_matches_serial(&test_lp(30, 48, true));
    }

    #[test]
    fn w5_sip_tasks_sparse_row_ep_uses_choose_column() {
        // 120 rows and a banded matrix: row_ep.count * inv_solver_num_row < 0.01 turns
        // slice_PRICE off (HEkkDual.cpp 1277), so `iterateTasks` calls `chooseColumn`.
        check_sip_matches_serial(&test_lp(120, 130, false));
    }
}
