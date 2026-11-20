//! HiGHS source: simplex/HEkkPrimal.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The primal simplex solver struct. Bodies: `solve.rs`, `iterate.rs`, `row.rs`, `hyper.rs`,
//! `weights.rs`, `report.rs` (HEkkPrimal.cpp by line range).

use crate::prelude::*;

mod hyper;
mod iterate;
mod report;
mod row;
mod solve;
mod weights;

/// HiGHS: HEkkPrimal.h::HEkkPrimal
/// # Removed C++ aliases (pattern (a))
///
/// | C++ member      | Rust expression          |
/// |-----------------|--------------------------|
/// | `ekk_instance_` | `ekk` (method parameter) |
/// | `analysis`      | `ekk.analysis_`          |
#[derive(Clone, Debug)]
pub struct HEkkPrimal {
    pub num_col: HighsInt,
    pub num_row: HighsInt,
    pub num_tot: HighsInt,
    pub solve_phase: HighsInt,
    pub edge_weight_mode: EdgeWeightMode,
    pub primal_feasibility_tolerance: f64,
    pub dual_feasibility_tolerance: f64,
    pub objective_target: f64,
    pub rebuild_reason: HighsInt,
    pub variable_in: HighsInt,
    pub move_in: HighsInt,
    pub row_out: HighsInt,
    pub variable_out: HighsInt,
    pub move_out: HighsInt,
    pub theta_dual: f64,
    pub theta_primal: f64,
    pub value_in: f64,
    pub alpha_col: f64,
    pub alpha_row: f64,
    pub numerical_trouble: f64,
    pub num_flip_since_rebuild: HighsInt,
    pub ph1_sorter_r: Vec<(f64, HighsInt)>,
    pub ph1_sorter_t: Vec<(f64, HighsInt)>,
    pub edge_weight_: Vec<f64>,
    pub num_devex_iterations_: HighsInt,
    pub num_bad_devex_weight_: HighsInt,
    pub devex_index_: Vec<HighsInt>,
    pub num_free_col: HighsInt,
    pub nonbasic_free_col_set: HSet,
    pub use_hyper_chuzc: bool,
    pub initialise_hyper_chuzc: bool,
    pub done_next_chuzc: bool,
    pub max_num_hyper_chuzc_candidates: HighsInt,
    pub num_hyper_chuzc_candidates: HighsInt,
    pub hyper_chuzc_candidate: Vec<HighsInt>,
    pub hyper_chuzc_measure: Vec<f64>,
    pub hyper_chuzc_candidate_set: HSet,
    pub max_hyper_chuzc_non_candidate_measure: f64,
    pub max_changed_measure_value: f64,
    pub max_changed_measure_column: HighsInt,
    pub report_hyper_chuzc: bool,
    pub row_ep: HVector,
    pub row_ap: HVector,
    pub col_aq: HVector,
    pub col_basic_feasibility_change: HVector,
    pub row_basic_feasibility_change: HVector,
    pub col_steepest_edge: HVector,
    pub random_: HighsRandom,
    pub max_max_local_primal_infeasibility_: f64,
    pub max_max_ignored_violation_: f64,
    pub max_max_primal_correction_: f64,
    pub last_header_iteration_count_: HighsInt,
    pub primal_correction_strategy: HighsInt,
    pub debug_max_relative_primal_steepest_edge_weight_error: f64,
    pub check_iter: HighsInt,
    pub check_column: HighsInt,
}

impl Default for HEkkPrimal {
    fn default() -> Self {
        Self {
            num_col: Default::default(),
            num_row: Default::default(),
            num_tot: Default::default(),
            solve_phase: Default::default(),
            edge_weight_mode: Default::default(),
            primal_feasibility_tolerance: Default::default(),
            dual_feasibility_tolerance: Default::default(),
            objective_target: Default::default(),
            rebuild_reason: Default::default(),
            variable_in: Default::default(),
            move_in: Default::default(),
            row_out: Default::default(),
            variable_out: Default::default(),
            move_out: Default::default(),
            theta_dual: Default::default(),
            theta_primal: Default::default(),
            value_in: Default::default(),
            alpha_col: Default::default(),
            alpha_row: Default::default(),
            numerical_trouble: Default::default(),
            num_flip_since_rebuild: Default::default(),
            ph1_sorter_r: Default::default(),
            ph1_sorter_t: Default::default(),
            edge_weight_: Default::default(),
            num_devex_iterations_: Default::default(),
            num_bad_devex_weight_: Default::default(),
            devex_index_: Default::default(),
            num_free_col: Default::default(),
            nonbasic_free_col_set: Default::default(),
            use_hyper_chuzc: false,
            initialise_hyper_chuzc: Default::default(),
            done_next_chuzc: Default::default(),
            max_num_hyper_chuzc_candidates: 50,
            num_hyper_chuzc_candidates: Default::default(),
            hyper_chuzc_candidate: Default::default(),
            hyper_chuzc_measure: Default::default(),
            hyper_chuzc_candidate_set: Default::default(),
            max_hyper_chuzc_non_candidate_measure: Default::default(),
            max_changed_measure_value: Default::default(),
            max_changed_measure_column: Default::default(),
            report_hyper_chuzc: false,
            row_ep: Default::default(),
            row_ap: Default::default(),
            col_aq: Default::default(),
            col_basic_feasibility_change: Default::default(),
            row_basic_feasibility_change: Default::default(),
            col_steepest_edge: Default::default(),
            random_: Default::default(),
            max_max_local_primal_infeasibility_: Default::default(),
            max_max_ignored_violation_: Default::default(),
            max_max_primal_correction_: Default::default(),
            last_header_iteration_count_: Default::default(),
            primal_correction_strategy: kSimplexPrimalCorrectionStrategyAlways,
            debug_max_relative_primal_steepest_edge_weight_error: 0.0,
            check_iter: 9999999,
            check_column: -2133,
        }
    }
}

impl HEkkPrimal {
    /// HiGHS: HEkkPrimal.h::HEkkPrimal::HEkkPrimal(HEkk& simplex) (header constructor):
    /// `ekk_instance_(simplex) { initialiseInstance(); }`
    pub fn new(ekk: &mut HEkk) -> Self {
        let mut primal = HEkkPrimal::default();
        primal.initialise_instance(ekk);
        primal
    }
}
