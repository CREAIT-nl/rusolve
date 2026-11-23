//! HiGHS source: pdlp/HiPdlpWrapper.h, pdlp/HiPdlpWrapper.cpp (lines 1-141)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiGHS -> HiPDLP bridge (`solver=hipdlp`). HiGHS C++ functions, so
//! snake_case names (repo convention). The `PDLP_DEBUG_LOG` and
//! `PDLP_PROFILE` blocks are not compiled.

use super::hipdlp::*;
use crate::lp_data::highs_solution::reset_model_status_and_highs_info;
use crate::prelude::*;

/// HiGHS: HiPdlpWrapper.cpp::solveLpHiPdlp(HighsLpSolverObject&) (lines 19-24)
pub fn solve_lp_hipdlp_so(solver_object: &mut HighsLpSolverObject) -> HighsStatus {
    solve_lp_hipdlp(
        solver_object.options_,
        solver_object.timer_,
        solver_object.lp_,
        solver_object.basis_,
        solver_object.solution_,
        &mut solver_object.model_status_,
        solver_object.highs_info_,
        solver_object.callback_,
    )
}

/// HiGHS: HiPdlpWrapper.cpp::solveLpHiPdlp (lines 26-141)
///
/// `resetModelStatusAndHighsInfo`; `PDLPSolver pdlp`; `setup(options,
/// &*timer)` (the solver keeps a shared borrow of the timer);
/// `passLp(lp)`; `preprocessLp`; `scaleProblem`; `x` / `y` zero vectors
/// of `getNumCol()` / `getNumRow()`; `solve`; `unscaleSolution`;
/// `postprocess` into a fresh `HighsSolution` (not OK -> `assert(111 ==
/// 444)` as `debug_assert!`, return `kError`); `logSummary` (print only);
/// `highs_info.pdlp_iteration_count = getIterationCount()`; `model_status
/// = kUnknown`; `highs_solution = pdlp_solution` (the `clear()` before the
/// assignment has no observable effect); `highs_basis.valid = false`.
/// Status: OPTIMAL -> `kOptimal`; INFEASIBLE -> `debug_assert!` (111 ==
/// 222), `kInfeasible`, return `kOk` without setting the valid flags;
/// UNBOUNDED -> likewise (111 == 333) `kUnbounded`; MAXITER ->
/// `kIterationLimit`; TIMEOUT -> `kTimeLimit`; anything else (ERROR) ->
/// return `kError`. Then `value_valid = dual_valid = true`, `kOk`.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::assertions_on_constants)]
pub fn solve_lp_hipdlp(
    options: &HighsOptions,
    timer: &mut HighsTimer,
    lp: &HighsLp,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
    model_status: &mut HighsModelStatus,
    highs_info: &mut HighsInfo,
    callback: &mut HighsCallback,
) -> HighsStatus {
    // Indicate that no imprecise solution has (yet) been found
    reset_model_status_and_highs_info(model_status, highs_info);

    /*** Order of operations
     * Preprocess with HiPdlp
     * Scale with HiPdlp
     * Solve with HiPdlp
     * Unscale with HiPdlp
     * Postprocess with HiPDLP
     * ***/

    // 0. Set up logger and params
    // PDLPSolver<'a> borrows the timer for the life of the solver.
    let mut pdlp = PDLPSolver::default();
    pdlp.setup(options, &*timer);

    // 1. Pass the LP to be solved
    pdlp.passLp(lp);

    // 2. Preprocess with HiPdlp
    pdlp.preprocessLp();

    // 3. Scale with HiPdlp
    pdlp.scaleProblem();

    // 4. Solve with HiPdlp
    let mut x: Vec<f64> = vec![0.0; pdlp.getNumCol() as usize];
    let mut y: Vec<f64> = vec![0.0; pdlp.getNumRow() as usize];

    pdlp.solve(&mut x, &mut y);

    // 5. Unscale with HiPdlp
    pdlp.unscaleSolution(&mut x, &mut y);

    // 6. Postprocess with HiPDLP
    let mut pdlp_solution = HighsSolution::default();
    let retcode = pdlp.postprocess(&mut pdlp_solution);
    if retcode != PostSolveRetcode::OK {
        debug_assert!(111 == 444);
        return HighsStatus::kError;
    }

    // #if PDLP_DEBUG_LOG (lines 73-79): not compiled.
    // #if PDLP_PROFILE (lines 82-84): not compiled.

    // --- Print Summary ---
    pdlp.logSummary();

    highs_info.pdlp_iteration_count = pdlp.getIterationCount();

    *model_status = HighsModelStatus::kUnknown;
    // highs_solution.clear(); has no observable effect before assignment.
    *highs_solution = pdlp_solution;
    highs_basis.valid = false;

    let termination_status = pdlp.getTerminationCode();
    match termination_status {
        TerminationStatus::OPTIMAL => {
            *model_status = HighsModelStatus::kOptimal;
        }
        TerminationStatus::INFEASIBLE => {
            debug_assert!(111 == 222);
            *model_status = HighsModelStatus::kInfeasible;
            return HighsStatus::kOk;
        }
        TerminationStatus::UNBOUNDED => {
            debug_assert!(111 == 333);
            *model_status = HighsModelStatus::kUnbounded;
            return HighsStatus::kOk;
        }
        TerminationStatus::MAXITER => {
            *model_status = HighsModelStatus::kIterationLimit;
        }
        TerminationStatus::TIMEOUT => {
            *model_status = HighsModelStatus::kTimeLimit;
        }
        _ => {
            debug_assert!(termination_status == TerminationStatus::ERROR);
            return HighsStatus::kError;
        }
    }
    debug_assert!(
        termination_status == TerminationStatus::OPTIMAL
            || termination_status == TerminationStatus::MAXITER
            || termination_status == TerminationStatus::TIMEOUT
    );
    // Commented-out code (lines 129-137): dropped.
    highs_solution.value_valid = true;
    highs_solution.dual_valid = true;
    HighsStatus::kOk
}

#[cfg(test)]
mod tests {
    //! Unit tests. End-to-end parity is the oracle
    //! test `tests/t_w3_hipdlp.rs`.
    use super::*;

    /// min x s.t. x >= 1, 0 <= x <= 10: HiPDLP reaches kOptimal near x = 1.
    #[test]
    fn one_variable_lp_is_optimal() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let mut options = HighsOptions::default();
        options.output_flag = false;
        let mut timer = HighsTimer::default();
        let mut basis = HighsBasis::default();
        let mut solution = HighsSolution::default();
        let mut status = HighsModelStatus::kNotset;
        let mut info = HighsInfo::default();
        let mut callback = HighsCallback::default();
        let rc = solve_lp_hipdlp(
            &options,
            &mut timer,
            &lp,
            &mut basis,
            &mut solution,
            &mut status,
            &mut info,
            &mut callback,
        );
        assert_eq!(rc, HighsStatus::kOk);
        assert_eq!(status, HighsModelStatus::kOptimal);
        assert!((solution.col_value[0] - 1.0).abs() < 1e-4);
        assert!(solution.value_valid && solution.dual_valid);
        assert!(!basis.valid);
        assert!(info.pdlp_iteration_count > 0);
    }
}
