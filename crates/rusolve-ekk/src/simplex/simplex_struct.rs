//! HiGHS source: simplex/SimplexStruct.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Data only. `SimplexBasis::clear/setup` are in `simplex/hsimplex.rs` (HSimplex.cpp);
//! `HighsRayRecord` members are in `simplex/ekk/proof.rs` (HEkk.cpp).

use crate::prelude::*;

/// HiGHS: SimplexStruct.h::SimplexBasis
#[derive(Clone, Debug)]
pub struct SimplexBasis {
    pub basic_index_: Vec<HighsInt>,
    pub nonbasic_flag_: Vec<i8>,
    pub nonbasic_move_: Vec<i8>,
    pub hash: u64,
    pub debug_id: HighsInt,
    pub debug_update_count: HighsInt,
    pub debug_origin_name: String,
}

impl Default for SimplexBasis {
    fn default() -> Self {
        Self {
            basic_index_: Default::default(),
            nonbasic_flag_: Default::default(),
            nonbasic_move_: Default::default(),
            hash: Default::default(),
            debug_id: -1,
            debug_update_count: -1,
            debug_origin_name: String::new(),
        }
    }
}

/// HiGHS: SimplexStruct.h::HighsSimplexStatus
#[derive(Clone, Debug)]
pub struct HighsSimplexStatus {
    pub initialised_for_new_lp: bool,
    pub is_dualized: bool,
    pub is_permuted: bool,
    pub initialised_for_solve: bool,
    pub has_basis: bool,
    pub has_ar_matrix: bool,
    pub has_nla: bool,
    pub has_dual_steepest_edge_weights: bool,
    pub has_invert: bool,
    pub has_fresh_invert: bool,
    pub has_fresh_rebuild: bool,
    pub has_dual_objective_value: bool,
    pub has_primal_objective_value: bool,
}

impl Default for HighsSimplexStatus {
    fn default() -> Self {
        Self {
            initialised_for_new_lp: false,
            is_dualized: false,
            is_permuted: false,
            initialised_for_solve: false,
            has_basis: false,
            has_ar_matrix: false,
            has_nla: false,
            has_dual_steepest_edge_weights: false,
            has_invert: false,
            has_fresh_invert: false,
            has_fresh_rebuild: false,
            has_dual_objective_value: false,
            has_primal_objective_value: false,
        }
    }
}

/// HiGHS: SimplexStruct.h::HighsSimplexInfo
#[derive(Clone, Debug)]
pub struct HighsSimplexInfo {
    pub work_cost_: Vec<f64>,
    pub work_dual_: Vec<f64>,
    pub work_shift_: Vec<f64>,
    pub work_lower_: Vec<f64>,
    pub work_upper_: Vec<f64>,
    pub work_range_: Vec<f64>,
    pub work_value_: Vec<f64>,
    pub work_lower_shift_: Vec<f64>,
    pub work_upper_shift_: Vec<f64>,
    pub base_lower_: Vec<f64>,
    pub base_upper_: Vec<f64>,
    pub base_value_: Vec<f64>,
    pub num_tot_random_value_: Vec<f64>,
    pub num_tot_permutation_: Vec<HighsInt>,
    pub num_col_permutation_: Vec<HighsInt>,
    pub devex_index_: Vec<HighsInt>,
    pub index_chosen_: Vec<HighsInt>,
    pub pivot_: Vec<f64>,
    pub phase1_backtracking_test_done: HighsInt,
    pub phase2_backtracking_test_done: HighsInt,
    pub backtracking_: bool,
    pub valid_backtracking_basis_: bool,
    pub backtracking_basis_: SimplexBasis,
    pub backtracking_basis_costs_shifted_: HighsInt,
    pub backtracking_basis_costs_perturbed_: HighsInt,
    pub backtracking_basis_bounds_shifted_: HighsInt,
    pub backtracking_basis_bounds_perturbed_: HighsInt,
    pub backtracking_basis_work_shift_: Vec<f64>,
    pub backtracking_basis_work_lower_shift_: Vec<f64>,
    pub backtracking_basis_work_upper_shift_: Vec<f64>,
    pub backtracking_basis_edge_weight_: Vec<f64>,
    pub simplex_strategy: HighsInt,
    pub dual_edge_weight_strategy: HighsInt,
    pub primal_edge_weight_strategy: HighsInt,
    pub price_strategy: HighsInt,
    pub dual_simplex_cost_perturbation_multiplier: f64,
    pub primal_simplex_phase1_cost_perturbation_multiplier: f64,
    pub primal_simplex_bound_perturbation_multiplier: f64,
    pub factor_pivot_threshold: f64,
    pub update_limit: HighsInt,
    pub control_iteration_count0: HighsInt,
    pub col_aq_density: f64,
    pub row_ep_density: f64,
    pub row_ap_density: f64,
    pub row_dse_density: f64,
    pub col_steepest_edge_density: f64,
    pub col_basic_feasibility_change_density: f64,
    pub row_basic_feasibility_change_density: f64,
    pub col_bfrt_density: f64,
    pub primal_col_density: f64,
    pub dual_col_density: f64,
    pub allow_dual_steepest_edge_to_devex_switch: bool,
    pub dual_steepest_edge_weight_log_error_threshold: f64,
    pub costly_dse_frequency: f64,
    pub num_costly_dse_iteration: HighsInt,
    pub costly_dse_measure: f64,
    pub average_log_low_dse_weight_error: f64,
    pub average_log_high_dse_weight_error: f64,
    pub run_quiet: bool,
    pub store_squared_primal_infeasibility: bool,
    pub report_simplex_inner_clock: bool,
    pub report_simplex_outer_clock: bool,
    pub report_simplex_phases_clock: bool,
    pub report_h_factor_clock: bool,
    pub analyse_lp: bool,
    pub analyse_iterations: bool,
    pub analyse_invert_form: bool,
    pub allow_cost_shifting: bool,
    pub allow_cost_perturbation: bool,
    pub allow_bound_perturbation: bool,
    pub costs_shifted: bool,
    pub costs_perturbed: bool,
    pub bounds_shifted: bool,
    pub bounds_perturbed: bool,
    pub num_primal_infeasibilities: HighsInt,
    pub max_primal_infeasibility: f64,
    pub sum_primal_infeasibilities: f64,
    pub num_dual_infeasibilities: HighsInt,
    pub max_dual_infeasibility: f64,
    pub sum_dual_infeasibilities: f64,
    pub dual_phase1_iteration_count: HighsInt,
    pub dual_phase2_iteration_count: HighsInt,
    pub primal_phase1_iteration_count: HighsInt,
    pub primal_phase2_iteration_count: HighsInt,
    pub primal_bound_swap: HighsInt,
    pub iteration_count0: HighsInt,
    pub dual_phase1_iteration_count0: HighsInt,
    pub dual_phase2_iteration_count0: HighsInt,
    pub primal_phase1_iteration_count0: HighsInt,
    pub primal_phase2_iteration_count0: HighsInt,
    pub primal_bound_swap0: HighsInt,
    pub min_concurrency: HighsInt,
    pub num_concurrency: HighsInt,
    pub max_concurrency: HighsInt,
    pub multi_iteration: HighsInt,
    pub update_count: HighsInt,
    pub dual_objective_value: f64,
    pub primal_objective_value: f64,
    pub updated_dual_objective_value: f64,
    pub updated_primal_objective_value: f64,
    pub num_basic_logicals: HighsInt,
}

impl Default for HighsSimplexInfo {
    fn default() -> Self {
        Self {
            work_cost_: Default::default(),
            work_dual_: Default::default(),
            work_shift_: Default::default(),
            work_lower_: Default::default(),
            work_upper_: Default::default(),
            work_range_: Default::default(),
            work_value_: Default::default(),
            work_lower_shift_: Default::default(),
            work_upper_shift_: Default::default(),
            base_lower_: Default::default(),
            base_upper_: Default::default(),
            base_value_: Default::default(),
            num_tot_random_value_: Default::default(),
            num_tot_permutation_: Default::default(),
            num_col_permutation_: Default::default(),
            devex_index_: Default::default(),
            index_chosen_: Default::default(),
            pivot_: Default::default(),
            phase1_backtracking_test_done: 0,
            phase2_backtracking_test_done: 0,
            backtracking_: false,
            valid_backtracking_basis_: false,
            backtracking_basis_: Default::default(),
            backtracking_basis_costs_shifted_: Default::default(),
            backtracking_basis_costs_perturbed_: Default::default(),
            backtracking_basis_bounds_shifted_: Default::default(),
            backtracking_basis_bounds_perturbed_: Default::default(),
            backtracking_basis_work_shift_: Default::default(),
            backtracking_basis_work_lower_shift_: Default::default(),
            backtracking_basis_work_upper_shift_: Default::default(),
            backtracking_basis_edge_weight_: Default::default(),
            simplex_strategy: Default::default(),
            dual_edge_weight_strategy: Default::default(),
            primal_edge_weight_strategy: Default::default(),
            price_strategy: Default::default(),
            dual_simplex_cost_perturbation_multiplier: Default::default(),
            primal_simplex_phase1_cost_perturbation_multiplier: 1.0,
            primal_simplex_bound_perturbation_multiplier: Default::default(),
            factor_pivot_threshold: Default::default(),
            update_limit: Default::default(),
            control_iteration_count0: Default::default(),
            col_aq_density: Default::default(),
            row_ep_density: Default::default(),
            row_ap_density: Default::default(),
            row_dse_density: Default::default(),
            col_steepest_edge_density: Default::default(),
            col_basic_feasibility_change_density: Default::default(),
            row_basic_feasibility_change_density: Default::default(),
            col_bfrt_density: Default::default(),
            primal_col_density: Default::default(),
            dual_col_density: Default::default(),
            allow_dual_steepest_edge_to_devex_switch: Default::default(),
            dual_steepest_edge_weight_log_error_threshold: Default::default(),
            costly_dse_frequency: Default::default(),
            num_costly_dse_iteration: Default::default(),
            costly_dse_measure: Default::default(),
            average_log_low_dse_weight_error: Default::default(),
            average_log_high_dse_weight_error: Default::default(),
            run_quiet: false,
            store_squared_primal_infeasibility: false,
            report_simplex_inner_clock: false,
            report_simplex_outer_clock: false,
            report_simplex_phases_clock: false,
            report_h_factor_clock: false,
            analyse_lp: false,
            analyse_iterations: false,
            analyse_invert_form: false,
            allow_cost_shifting: true,
            allow_cost_perturbation: true,
            allow_bound_perturbation: true,
            costs_shifted: false,
            costs_perturbed: false,
            bounds_shifted: false,
            bounds_perturbed: false,
            num_primal_infeasibilities: -1,
            max_primal_infeasibility: Default::default(),
            sum_primal_infeasibilities: Default::default(),
            num_dual_infeasibilities: -1,
            max_dual_infeasibility: Default::default(),
            sum_dual_infeasibilities: Default::default(),
            dual_phase1_iteration_count: 0,
            dual_phase2_iteration_count: 0,
            primal_phase1_iteration_count: 0,
            primal_phase2_iteration_count: 0,
            primal_bound_swap: 0,
            iteration_count0: 0,
            dual_phase1_iteration_count0: 0,
            dual_phase2_iteration_count0: 0,
            primal_phase1_iteration_count0: 0,
            primal_phase2_iteration_count0: 0,
            primal_bound_swap0: 0,
            min_concurrency: 1,
            num_concurrency: 1,
            max_concurrency: kSimplexConcurrencyLimit,
            multi_iteration: 0,
            update_count: Default::default(),
            dual_objective_value: Default::default(),
            primal_objective_value: Default::default(),
            updated_dual_objective_value: Default::default(),
            updated_primal_objective_value: Default::default(),
            num_basic_logicals: Default::default(),
        }
    }
}

/// HiGHS: SimplexStruct.h::HighsSimplexBadBasisChangeRecord
#[derive(Clone, Debug)]
pub struct HighsSimplexBadBasisChangeRecord {
    pub taboo: bool,
    pub row_out: HighsInt,
    pub variable_out: HighsInt,
    pub variable_in: HighsInt,
    pub reason: BadBasisChangeReason,
    pub save_value: f64,
}

impl Default for HighsSimplexBadBasisChangeRecord {
    fn default() -> Self {
        Self {
            taboo: Default::default(),
            row_out: Default::default(),
            variable_out: Default::default(),
            variable_in: Default::default(),
            reason: Default::default(),
            save_value: 0.0,
        }
    }
}

/// HiGHS: SimplexStruct.h::HighsRayRecord
#[derive(Clone, Debug)]
pub struct HighsRayRecord {
    pub index: HighsInt,
    pub sign: HighsInt,
    pub value: Vec<f64>,
}

impl Default for HighsRayRecord {
    fn default() -> Self {
        Self {
            index: Default::default(),
            sign: Default::default(),
            value: Default::default(),
        }
    }
}
