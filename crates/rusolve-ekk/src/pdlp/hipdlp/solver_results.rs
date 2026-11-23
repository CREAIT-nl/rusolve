//! HiGHS source: pdlp/hipdlp/solver_results.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Termination codes and the per-check convergence record of HiPDLP.

/// HiGHS: solver_results.hpp `enum class TerminationStatus`
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TerminationStatus {
    #[default]
    NOTSET,
    OPTIMAL,
    INFEASIBLE,
    UNBOUNDED,
    MAXITER,
    TIMEOUT,
    ERROR,
}

/// HiGHS: solver_results.hpp `enum class TerminationIterate`
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TerminationIterate {
    #[default]
    LAST_ITERATE,
    AVERAGE_ITERATE,
}

/// HiGHS: solver_results.hpp `struct SolverResults`
#[derive(Clone, Debug)]
pub struct SolverResults {
    pub term_code: TerminationStatus,
    pub term_iterate: TerminationIterate,
    pub primal_obj: f64,
    pub dual_obj: f64,
    pub duality_gap: f64,
    pub complementarity: f64,
    pub primal_feasibility: f64,
    pub dual_feasibility: f64,
    pub relative_obj_gap: f64,
    // --- Averaged Metrics ---
    pub primal_obj_average: f64,
    pub dual_obj_average: f64,
    pub duality_gap_average: f64,
    pub primal_feasibility_average: f64,
    pub dual_feasibility_average: f64,
    // --- Residual Vectors ---
    pub primal_residual: Vec<f64>,
    pub dual_residual: Vec<f64>,
    pub primal_residual_average: Vec<f64>,
    pub dual_residual_average: Vec<f64>,
    // --- Infeasibility Detection Metrics ---
    pub primal_infeasibility_obj: f64,
    pub dual_infeasibility_obj: f64,
    pub primal_infeasibility_res: f64,
    pub dual_infeasibility_res: f64,
    // --- Restart Statistics ---
    pub primal_feasibility_last_restart: f64,
    pub dual_feasibility_last_restart: f64,
    pub duality_gap_last_restart: f64,
}

impl Default for SolverResults {
    /// Note: `term_code` defaults to `TIMEOUT`, not `NOTSET` (header value).
    fn default() -> Self {
        Self {
            term_code: TerminationStatus::TIMEOUT,
            term_iterate: TerminationIterate::LAST_ITERATE,
            primal_obj: 0.0,
            dual_obj: 0.0,
            duality_gap: 0.0,
            complementarity: 0.0,
            primal_feasibility: 0.0,
            dual_feasibility: 0.0,
            relative_obj_gap: 0.0,
            primal_obj_average: 0.0,
            dual_obj_average: 0.0,
            duality_gap_average: 0.0,
            primal_feasibility_average: 0.0,
            dual_feasibility_average: 0.0,
            primal_residual: Vec::new(),
            dual_residual: Vec::new(),
            primal_residual_average: Vec::new(),
            dual_residual_average: Vec::new(),
            primal_infeasibility_obj: 0.0,
            dual_infeasibility_obj: 0.0,
            primal_infeasibility_res: 0.0,
            dual_infeasibility_res: 0.0,
            primal_feasibility_last_restart: 0.0,
            dual_feasibility_last_restart: 0.0,
            duality_gap_last_restart: 0.0,
        }
    }
}
