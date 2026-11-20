//! HiGHS source: simplex/HEkk.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The `HEkk` struct (EKK simplex instance). Method bodies are split by HEkk.cpp line
//! ranges: `clear.rs`, `dualize.rs`, `solve.rs`, `factor.rs`, `init.rs`, `iterate.rs`,
//! `solve_return.rs`, `proof.rs`; plus `control.rs` (HEkkControl.cpp, HEkkInterface.cpp) and
//! `debug.rs` (HEkkDebug.cpp, stubs).
//!
//! Ownership (crate docs, pattern (d)): `callback_`, `options_` and `timer_` are owned values
//! (C++ pointers): translate `options_->x` as `self.options_.x`, `timer_->read()` as
//! `self.timer_.read()`. `HEkkDual`/`HEkkPrimal` are created locally in `solve` and receive
//! `self` as a parameter (pattern (a)); `simplex_nla_` methods receive `NlaRefs` (pattern (b)).

use crate::prelude::*;

mod clear;
mod control;
mod debug;
mod dualize;
mod factor;
mod init;
mod iterate;
mod proof;
mod solve;
mod solve_return;

/// HiGHS: HEkk.h::HEkk
/// All members are `pub` (C++ private members included) so the split impl files and
/// the dual/primal solvers can reach them.
/// `Default` = the C++ constructor's initialiser list.
#[derive(Clone, Debug)]
pub struct HEkk {
    pub callback_: HighsCallback,
    pub options_: HighsOptions,
    pub timer_: HighsTimer,
    pub analysis_: HighsSimplexAnalysis,
    pub lp_: HighsLp,
    pub lp_name_: String,
    pub status_: HighsSimplexStatus,
    pub info_: HighsSimplexInfo,
    pub model_status_: HighsModelStatus,
    pub basis_: SimplexBasis,
    pub visited_basis_: HashSet<u64>,
    pub random_: HighsRandom,
    pub dual_edge_weight_: Vec<f64>,
    pub scattered_dual_edge_weight_: Vec<f64>,
    pub simplex_in_scaled_space_: bool,
    pub ar_matrix_: HighsSparseMatrix,
    pub scaled_a_matrix_: HighsSparseMatrix,
    pub simplex_nla_: HSimplexNla,
    pub hot_start_: HotStart,
    pub cost_scale_: f64,
    pub cost_perturbation_base_: f64,
    pub cost_perturbation_max_abs_cost_: f64,
    pub iteration_count_: HighsInt,
    pub dual_simplex_cleanup_level_: HighsInt,
    pub dual_simplex_phase1_cleanup_level_: HighsInt,
    pub previous_iteration_cycling_detected: HighsInt,
    pub solve_bailout_: bool,
    pub called_return_from_solve_: bool,
    pub exit_algorithm_: SimplexAlgorithm,
    pub return_primal_solution_status_: HighsInt,
    pub return_dual_solution_status_: HighsInt,
    pub proof_index_: Vec<HighsInt>,
    pub proof_value_: Vec<f64>,
    pub dual_ray_record_: HighsRayRecord,
    pub primal_ray_record_: HighsRayRecord,
    pub original_num_col_: HighsInt,
    pub original_num_row_: HighsInt,
    pub original_num_nz_: HighsInt,
    pub original_offset_: f64,
    pub original_col_cost_: Vec<f64>,
    pub original_col_lower_: Vec<f64>,
    pub original_col_upper_: Vec<f64>,
    pub original_row_lower_: Vec<f64>,
    pub original_row_upper_: Vec<f64>,
    pub upper_bound_col_: Vec<HighsInt>,
    pub upper_bound_row_: Vec<HighsInt>,
    pub edge_weight_error_: f64,
    pub build_synthetic_tick_: f64,
    pub total_synthetic_tick_: f64,
    pub debug_solve_call_num_: HighsInt,
    pub debug_basis_id_: HighsInt,
    pub time_report_: bool,
    pub debug_initial_build_synthetic_tick_: HighsInt,
    pub debug_solve_report_: bool,
    pub debug_iteration_report_: bool,
    pub debug_basis_report_: bool,
    pub debug_dual_feasible: bool,
    pub debug_max_relative_dual_steepest_edge_weight_error: f64,
    pub bad_basis_change_: Vec<HighsSimplexBadBasisChangeRecord>,
    pub primal_phase1_dual_: Vec<f64>,
    pub simplex_stats_: HighsSimplexStats,
}

impl Default for HEkk {
    fn default() -> Self {
        Self {
            callback_: Default::default(),
            options_: Default::default(),
            timer_: Default::default(),
            analysis_: Default::default(),
            lp_: Default::default(),
            lp_name_: String::new(),
            status_: Default::default(),
            info_: Default::default(),
            model_status_: Default::default(),
            basis_: Default::default(),
            visited_basis_: Default::default(),
            random_: Default::default(),
            dual_edge_weight_: Default::default(),
            scattered_dual_edge_weight_: Default::default(),
            simplex_in_scaled_space_: Default::default(),
            ar_matrix_: Default::default(),
            scaled_a_matrix_: Default::default(),
            simplex_nla_: Default::default(),
            hot_start_: Default::default(),
            cost_scale_: 1.0,
            cost_perturbation_base_: Default::default(),
            cost_perturbation_max_abs_cost_: Default::default(),
            iteration_count_: Default::default(),
            dual_simplex_cleanup_level_: Default::default(),
            dual_simplex_phase1_cleanup_level_: Default::default(),
            previous_iteration_cycling_detected: -kHighsIInf,
            solve_bailout_: Default::default(),
            called_return_from_solve_: Default::default(),
            exit_algorithm_: Default::default(),
            return_primal_solution_status_: Default::default(),
            return_dual_solution_status_: Default::default(),
            proof_index_: Default::default(),
            proof_value_: Default::default(),
            dual_ray_record_: Default::default(),
            primal_ray_record_: Default::default(),
            original_num_col_: Default::default(),
            original_num_row_: Default::default(),
            original_num_nz_: Default::default(),
            original_offset_: Default::default(),
            original_col_cost_: Default::default(),
            original_col_lower_: Default::default(),
            original_col_upper_: Default::default(),
            original_row_lower_: Default::default(),
            original_row_upper_: Default::default(),
            upper_bound_col_: Default::default(),
            upper_bound_row_: Default::default(),
            edge_weight_error_: Default::default(),
            build_synthetic_tick_: Default::default(),
            total_synthetic_tick_: Default::default(),
            debug_solve_call_num_: Default::default(),
            debug_basis_id_: Default::default(),
            time_report_: Default::default(),
            debug_initial_build_synthetic_tick_: Default::default(),
            debug_solve_report_: Default::default(),
            debug_iteration_report_: Default::default(),
            debug_basis_report_: Default::default(),
            debug_dual_feasible: Default::default(),
            debug_max_relative_dual_steepest_edge_weight_error: Default::default(),
            bad_basis_change_: Default::default(),
            primal_phase1_dual_: Default::default(),
            simplex_stats_: Default::default(),
        }
    }
}

impl HEkk {
    /// HiGHS: HEkk.h::HEkk::getSimplexBasis (inline getter)
    pub fn get_simplex_basis(&self) -> &SimplexBasis {
        &self.basis_
    }

    /// HiGHS: HEkk.h::HEkk::computeBasisCondition() (inline no-argument overload)
    pub fn compute_basis_condition(&self) -> f64 {
        self.compute_basis_condition_lp(&self.lp_, false, false)
    }

    /// HiGHS: HEkk.h::HEkk::getSimplexStats (inline getter)
    pub fn get_simplex_stats(&self) -> &HighsSimplexStats {
        &self.simplex_stats_
    }

    /// HiGHS: HEkk.h::HEkk::initialiseSimplexStats (inline)
    pub fn initialise_simplex_stats(&mut self) {
        let iteration_count = self.iteration_count_;
        self.simplex_stats_.initialise(iteration_count);
    }

    /// HiGHS: HEkk.h::HEkk::reportSimplexStats (inline)
    /// C++ default arguments: message = ""
    pub fn report_simplex_stats(&self, file: Option<HighsLogStream>, message: &str) {
        self.simplex_stats_.report(file, message);
    }
}
