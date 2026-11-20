//! HiGHS source: lp_data/HighsInfo.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HighsInfoStruct fields only; the records machinery is dropped.
//! Default is all zero/false (C++ fills values via `invalidate`, out of scope).

use crate::prelude::*;

/// HiGHS: HighsInfo.h::HighsInfoStruct
#[derive(Clone, Debug, Default)]
pub struct HighsInfo {
    pub valid: bool,
    pub mip_node_count: i64,
    pub simplex_iteration_count: HighsInt,
    pub ipm_iteration_count: HighsInt,
    pub crossover_iteration_count: HighsInt,
    pub pdlp_iteration_count: HighsInt,
    pub qp_iteration_count: HighsInt,
    pub primal_solution_status: HighsInt,
    pub dual_solution_status: HighsInt,
    pub basis_validity: HighsInt,
    /// HiGHS 1.16 #3456 (h16_ac_centring): HighsInfo.h::centring_status (1.16.0). Inert when the
    /// option is off: written only by the IPX record under `h16_ac_centring`, read only by
    /// the analytic-centre task under `h16_ac_centring`. Default is `kCentringStatusNone`.
    pub centring_status: HighsInt,
    pub objective_function_value: f64,
    pub mip_dual_bound: f64,
    pub mip_gap: f64,
    pub max_integrality_violation: f64,
    pub num_primal_infeasibilities: HighsInt,
    pub max_primal_infeasibility: f64,
    pub sum_primal_infeasibilities: f64,
    pub num_dual_infeasibilities: HighsInt,
    pub max_dual_infeasibility: f64,
    pub sum_dual_infeasibilities: f64,
    pub num_semi_infeasibilities: HighsInt,
    pub max_semi_infeasibility: f64,
    pub sum_semi_infeasibilities: f64,
    pub num_relative_primal_infeasibilities: HighsInt,
    pub max_relative_primal_infeasibility: f64,
    pub num_relative_dual_infeasibilities: HighsInt,
    pub max_relative_dual_infeasibility: f64,
    pub num_primal_residual_errors: HighsInt,
    pub max_primal_residual_error: f64,
    pub num_dual_residual_errors: HighsInt,
    pub max_dual_residual_error: f64,
    pub num_relative_primal_residual_errors: HighsInt,
    pub max_relative_primal_residual_error: f64,
    pub num_relative_dual_residual_errors: HighsInt,
    pub max_relative_dual_residual_error: f64,
    pub num_complementarity_violations: HighsInt,
    pub max_complementarity_violation: f64,
    pub primal_dual_objective_error: f64,
    pub primal_dual_integral: f64,
}
