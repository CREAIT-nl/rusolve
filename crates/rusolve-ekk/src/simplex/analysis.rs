//! HiGHS source: simplex/HighsSimplexAnalysis.h (stub)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Out of scope: simplex analysis/reporting. This stub keeps every field and method that the
//! in-scope simplex code touches, with HiGHS defaults; all `analyse_*` flags are `false` so
//! the analysis branches in HEkk* are dead, and every method is a no-op (value-returning
//! methods return a neutral value). Filler agents: translate `analysis->x = y` as
//! `ekk.analysis_.x = y` and calls as calls on this stub; do not add logic here.

use crate::lp_data::highs_lp::HighsLp;
use crate::lp_data::highs_options::HighsOptions;
use crate::prelude::*;
use crate::util::hfactor::HFactor;
use crate::util::highs_timer::HighsTimerClock;
use crate::util::highs_utils::HighsValueDistribution;
use crate::util::hvector::HVector;

/// HiGHS: HighsSimplexAnalysis.h::HighsSimplexAnalysis (stub; field subset)
#[derive(Clone, Debug)]
pub struct HighsSimplexAnalysis {
    pub analyse_factor_data: bool,
    pub analyse_factor_time: bool,
    pub analyse_simplex_runtime_data: bool,
    pub analyse_simplex_summary_data: bool,
    pub analyse_simplex_time: bool,
    pub analyse_simplex_data: bool,
    pub analyse_lp_data: bool,
    pub col_aq_density: f64,
    pub col_basic_feasibility_change_density: f64,
    pub col_BFRT_density: f64,
    pub col_steepest_edge_density: f64,
    pub costly_DSE_measure: f64,
    pub devex_iteration_count: HighsInt,
    pub dual_col_density: f64,
    pub dual_step: f64,
    pub edge_weight: f64,
    pub edge_weight_error: f64,
    pub edge_weight_mode: EdgeWeightMode,
    pub entering_variable: HighsInt,
    pub factor_pivot_threshold: f64,
    pub leaving_variable: HighsInt,
    pub max_correct_dual_cost_shift: f64,
    pub max_correct_dual_cost_shift_dual_infeasibility: f64,
    pub max_correct_dual_primal_flip: f64,
    pub max_dual_phase_1_lp_dual_infeasibility: f64,
    pub max_heap_chuzc_size: HighsInt,
    pub max_quad_chuzc_size: HighsInt,
    pub max_single_cost_shift: f64,
    pub min_correct_dual_primal_flip_dual_infeasibility: f64,
    pub net_num_single_cost_shift: HighsInt,
    pub num_col_price: HighsInt,
    pub num_correct_dual_cost_shift: HighsInt,
    pub num_correct_dual_primal_flip: HighsInt,
    pub num_costly_DSE_iteration: HighsInt,
    pub num_devex_framework: HighsInt,
    pub num_dual_cycling_detections: HighsInt,
    pub num_dual_infeasibility: HighsInt,
    pub num_dual_phase_1_lp_dual_infeasibility: HighsInt,
    pub num_heap_chuzc: HighsInt,
    pub num_improve_choose_column_row_call: HighsInt,
    pub num_primal_cycling_detections: HighsInt,
    pub num_primal_infeasibility: HighsInt,
    pub num_quad_chuzc: HighsInt,
    pub num_remove_pivot_from_pack: HighsInt,
    pub num_row_price: HighsInt,
    pub num_row_price_with_switch: HighsInt,
    pub num_single_cost_shift: HighsInt,
    pub numerical_trouble: f64,
    pub objective_value: f64,
    pub pivot_value_from_column: f64,
    pub pivot_value_from_row: f64,
    pub pivotal_row_index: HighsInt,
    pub primal_col_density: f64,
    pub primal_delta: f64,
    pub primal_step: f64,
    pub rebuild_reason: HighsInt,
    pub rebuild_reason_string: String,
    pub reduced_cost_value: f64,
    pub reduced_rhs_value: f64,
    pub row_ap_density: f64,
    pub row_basic_feasibility_change_density: f64,
    pub row_DSE_density: f64,
    pub row_ep_density: f64,
    pub simplex_iteration_count: HighsInt,
    pub simplex_strategy: HighsInt,
    pub solve_phase: HighsInt,
    pub sum_dual_infeasibility: f64,
    pub sum_dual_phase_1_lp_dual_infeasibility: f64,
    pub sum_heap_chuzc_size: f64,
    pub sum_primal_infeasibility: f64,
    pub sum_quad_chuzc_size: f64,
    pub sum_single_cost_shift: f64,
    pub cost_perturbation1_distribution: HighsValueDistribution,
    pub cost_perturbation2_distribution: HighsValueDistribution,
}

impl Default for HighsSimplexAnalysis {
    fn default() -> Self {
        Self {
            analyse_factor_data: false,
            analyse_factor_time: false,
            analyse_simplex_runtime_data: false,
            analyse_simplex_summary_data: false,
            analyse_simplex_time: false,
            analyse_simplex_data: false,
            analyse_lp_data: false,
            col_aq_density: 0.0,
            col_basic_feasibility_change_density: 0.0,
            col_BFRT_density: 0.0,
            col_steepest_edge_density: 0.0,
            costly_DSE_measure: 0.0,
            devex_iteration_count: 0,
            dual_col_density: 0.0,
            dual_step: 0.0,
            edge_weight: 0.0,
            edge_weight_error: 0.0,
            edge_weight_mode: EdgeWeightMode::kSteepestEdge,
            entering_variable: 0,
            factor_pivot_threshold: 0.0,
            leaving_variable: 0,
            max_correct_dual_cost_shift: 0.0,
            max_correct_dual_cost_shift_dual_infeasibility: 0.0,
            max_correct_dual_primal_flip: 0.0,
            max_dual_phase_1_lp_dual_infeasibility: 0.0,
            max_heap_chuzc_size: 0,
            max_quad_chuzc_size: 0,
            max_single_cost_shift: 0.0,
            min_correct_dual_primal_flip_dual_infeasibility: kHighsInf,
            net_num_single_cost_shift: 0,
            num_col_price: 0,
            num_correct_dual_cost_shift: 0,
            num_correct_dual_primal_flip: 0,
            num_costly_DSE_iteration: 0,
            num_devex_framework: 0,
            num_dual_cycling_detections: 0,
            num_dual_infeasibility: 0,
            num_dual_phase_1_lp_dual_infeasibility: 0,
            num_heap_chuzc: 0,
            num_improve_choose_column_row_call: 0,
            num_primal_cycling_detections: 0,
            num_primal_infeasibility: 0,
            num_quad_chuzc: 0,
            num_remove_pivot_from_pack: 0,
            num_row_price: 0,
            num_row_price_with_switch: 0,
            num_single_cost_shift: 0,
            numerical_trouble: 0.0,
            objective_value: 0.0,
            pivot_value_from_column: 0.0,
            pivot_value_from_row: 0.0,
            pivotal_row_index: 0,
            primal_col_density: 0.0,
            primal_delta: 0.0,
            primal_step: 0.0,
            rebuild_reason: 0,
            rebuild_reason_string: Default::default(),
            reduced_cost_value: 0.0,
            reduced_rhs_value: 0.0,
            row_ap_density: 0.0,
            row_basic_feasibility_change_density: 0.0,
            row_DSE_density: 0.0,
            row_ep_density: 0.0,
            simplex_iteration_count: 0,
            simplex_strategy: 0,
            solve_phase: 0,
            sum_dual_infeasibility: 0.0,
            sum_dual_phase_1_lp_dual_infeasibility: 0.0,
            sum_heap_chuzc_size: 0.0,
            sum_primal_infeasibility: 0.0,
            sum_quad_chuzc_size: 0.0,
            sum_single_cost_shift: 0.0,
            cost_perturbation1_distribution: Default::default(),
            cost_perturbation2_distribution: Default::default(),
        }
    }
}

impl HighsSimplexAnalysis {
    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::setup (stub: no-op)
    pub fn setup(
        &mut self,
        lp_name: &str,
        lp: &HighsLp,
        options: &HighsOptions,
        simplex_iteration_count: HighsInt,
    ) {
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::iterationReport (stub: no-op)
    pub fn iteration_report(&mut self) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::invertReport (stub: no-op)
    pub fn invert_report(&mut self) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::userInvertReport (stub: no-op)
    pub fn user_invert_report(&mut self, force: bool) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::simplexTimerStart (stub: no-op)
    pub fn simplex_timer_start(&mut self, simplex_clock: HighsInt, thread_id: HighsInt) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::simplexTimerStop (stub: no-op)
    pub fn simplex_timer_stop(&mut self, simplex_clock: HighsInt, thread_id: HighsInt) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::simplexTimerRunning (stub: no-op)
    pub fn simplex_timer_running(&self, simplex_clock: HighsInt, thread_id: HighsInt) -> bool {
        false
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::simplexTimerNumCall (stub: no-op)
    pub fn simplex_timer_num_call(&self, simplex_clock: HighsInt, thread_id: HighsInt) -> HighsInt {
        0
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::simplexTimerRead (stub: no-op)
    pub fn simplex_timer_read(&self, simplex_clock: HighsInt, thread_id: HighsInt) -> f64 {
        0.0
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::getThreadFactorTimerClockPtr (stub: no-op)
    pub fn get_thread_factor_timer_clock_ptr(&self, i: HighsInt) -> Option<&HighsTimerClock> {
        None
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::iterationRecord (stub: no-op)
    pub fn iteration_record(&mut self) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::iterationRecordMajor (stub: no-op)
    pub fn iteration_record_major(&mut self) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::operationRecordBefore (stub: no-op)
    pub fn operation_record_before(
        &mut self,
        operation_type: HighsInt,
        vector: &HVector,
        historical_density: f64,
    ) {
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::operationRecordBefore (stub: no-op)
    pub fn operation_record_before_count(
        &mut self,
        operation_type: HighsInt,
        current_count: HighsInt,
        historical_density: f64,
    ) {
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::operationRecordAfter (stub: no-op)
    pub fn operation_record_after(&mut self, operation_type: HighsInt, vector: &HVector) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::operationRecordAfter (stub: no-op)
    pub fn operation_record_after_count(
        &mut self,
        operation_type: HighsInt,
        result_count: HighsInt,
    ) {
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::summaryReport (stub: no-op)
    pub fn summary_report(&mut self) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::reportSimplexTimer (stub: no-op)
    pub fn report_simplex_timer(&mut self) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::reportFactorTimer (stub: no-op)
    pub fn report_factor_timer(&mut self) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::updateInvertFormData (stub: no-op)
    pub fn update_invert_form_data(&mut self, factor: &HFactor) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::reportInvertFormData (stub: no-op)
    pub fn report_invert_form_data(&mut self) {}

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::numInvert (stub: no-op)
    pub fn num_invert(&self) -> HighsInt {
        0
    }

    /// HiGHS: HighsSimplexAnalysis.cpp::HighsSimplexAnalysis::dualSteepestEdgeWeightError (stub: no-op)
    pub fn dual_steepest_edge_weight_error(
        &mut self,
        computed_edge_weight: f64,
        updated_edge_weight: f64,
    ) {
    }
}
