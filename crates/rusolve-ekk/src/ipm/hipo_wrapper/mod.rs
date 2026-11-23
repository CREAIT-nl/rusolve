//! HiGHS source: ipm/IpxWrapper.h, ipm/IpxWrapper.cpp (HiPO part: lines 408-657, 1131-1142, 1170-1188, 1393-1511)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiGHS -> HiPO bridge (`solver=hipo`, or `solver=ipm` when HiPO is
//! available). In HiGHS these functions sit in `IpxWrapper.cpp` next to the
//! IPX ones ([`crate::ipm::ipx_wrapper`], which this module reuses:
//! `fill_in_rhs_and_constraints`, `ipx_status_error`,
//! `illegal_ipx_stopped_crossover_status`, and
//! `lp_data::highs_solution_ipx::{ipx_solution_to_highs_solution,
//! ipx_basic_solution_to_highs_basic_solution}`). HiGHS-side functions, so
//! snake_case names (crate convention).
//!
//! * `HighsExtras::blas::openblas_set_num_threads(1)` is a no-op: the BLAS is
//!   [`crate::ipm::hipo::blas`] (serial by construction).
//! * `hipo.setCallback(callback)` stores a pointer in C++; here the
//!   `HighsCallback` is moved into the solver with
//!   `hipo.setCallback(std::mem::take(callback))` and moved back with
//!   `*callback = hipo.take_callback().unwrap_or_default()` on every return
//!   path after it (same pattern as `solve_lp_ipx`).
//! * `hipo.setTimer(timer)` keeps a clone of the timer state (module doc of
//!   `crate::ipm::hipo`, "Ownership").
//! * Oracle trace hook: `oracle_trace::hipo_end(&hipo_info)` right after
//!   `const hipo::Info hipo_info = hipo.getInfo()` (`tools/highs-oracle/FORMAT.md`, HiPO records).
//! * C++ `try { ... } catch (const std::exception&)` around `solveLpHipo` in
//!   `solveLp` is not reproduced: a throw is a panic (as for IPX).

pub mod report;

pub use report::*;

use crate::prelude::*;

use crate::ipm::hipo::Int;
use crate::ipm::hipo::ipm::solver::Solver;
use crate::ipm::hipo::ipm::status::*;
use crate::ipm::hipo::oracle_trace;
use crate::ipm::ipx;
use crate::ipm::ipx_wrapper::{
    IpxSolution, fill_in_rhs_and_constraints, illegal_ipx_stopped_crossover_status,
    ipx_status_error,
};
use crate::lp_data::highs_solution::reset_model_status_and_highs_info;
use crate::lp_data::highs_solution_ipx::ipx_basic_solution_to_highs_basic_solution;
use crate::util::c_printf::c_fmt_g;

/// HiGHS: IpxWrapper.cpp::solveLpHipo(HighsLpSolverObject&) (lines 408-413)
///
/// The C++ passes a temporary `HighsHessian{}` (LP: no Hessian).
pub fn solve_lp_hipo_so(solver_object: &mut HighsLpSolverObject) -> HighsStatus {
    solve_hipo(
        &*solver_object.options_,
        &mut *solver_object.timer_,
        &*solver_object.lp_,
        &HighsHessian::default(),
        &mut *solver_object.basis_,
        &mut *solver_object.solution_,
        &mut solver_object.model_status_,
        &mut *solver_object.highs_info_,
        &mut *solver_object.callback_,
    )
}

/// rusolve lever `ipm_hipo_dualize` (Andersen-Andersen 1995,
/// Gondzio 1997). An LP with at least 3x more rows than columns is dualized with the ported
/// `HEkk::dualize`, the dual is solved by HiPO with crossover, and the dual basis is mapped
/// back by `HEkk::undualize`, which also re-solves the primal from that basis. Returns true
/// when `basis_`, `solution_` and `model_status_` of the solver object hold an optimal
/// primal answer; false leaves the solver object untouched so the caller solves the LP as is.
/// rusolve (dualize-guard): predicted factor flops and nnz of the
/// Newton system HiPO would pick for `lp`, from its symbolic analysis alone (no iterations,
/// no log). `None` when the model or the setup fails.
pub(crate) fn hipo_analysed_flops(
    options: &HighsOptions,
    timer: &mut HighsTimer,
    lp: &HighsLp,
) -> Option<(f64, f64)> {
    let mut quiet = options.clone();
    quiet.output_flag = false;
    quiet.sync_log_options();
    let mut hipo = Solver::default();
    hipo.setIpxTimerOffset(timer.read());
    hipo.setOptions(&quiet);
    hipo.setTimer(timer);
    hipo.analyse_only_ = true;
    if hipo.load(lp, &HighsHessian::default()) != 0 {
        return None;
    }
    hipo.solve();
    (hipo.analysed_flops_ > 0.0).then_some((hipo.analysed_flops_, hipo.analysed_nz_))
}

/// Dualize only when the dual's predicted factor is not clearly costlier than the primal's
/// (rule in `Solver::initialise`: dual flops above 1e9 and above 1.5x the primal's, unless
/// HiPO refuses the dual as too expensive). Below 1e9 flops a factorisation takes well under a
/// second, so the comparison is moot. Source: the fill-in comparison behind Andersen-Andersen 1995 (dualize
/// when the dual normal equations factor cheaper) with margins measured on our benchmark set.
pub fn solve_lp_hipo_dualized(solver_object: &mut HighsLpSolverObject) -> bool {
    use crate::simplex::ekk::HEkk;
    {
        let lp = &*solver_object.lp_;
        if !solver_object.options_.ipm_hipo_dualize
            || lp.num_row_ < 3 * lp.num_col_.max(1)
            || lp.a_matrix_.num_nz() == 0
            || lp.is_scaled_
            || lp.scale_.has_scaling
            || !lp.a_matrix_.is_colwise()
        {
            return false;
        }
    }
    let mut ekk = HEkk::default();
    ekk.set_pointers(
        &*solver_object.callback_,
        &*solver_object.options_,
        &*solver_object.timer_,
    );
    ekk.lp_ = solver_object.lp_.clone();
    ekk.initialise_ekk();
    if ekk.dualize() != HighsStatus::kOk {
        return false;
    }
    let dual_lp = ekk.lp_.clone();
    // The primal is the side that is analysed as extra work; the dual's own analysis is the
    // one of the real solve, which stops after it when the dual is clearly costlier.
    let pf = hipo_analysed_flops(
        solver_object.options_,
        &mut *solver_object.timer_,
        &*solver_object.lp_,
    );
    if let Some((f, n)) = pf {
        crate::highs_log_user!(
            &solver_object.options_.log_options,
            HighsLogType::kInfo,
            "HiPO dualize guard: primal flops {:.2e} nnz {:.2e}\n",
            f,
            n
        );
    }
    crate::ipm::hipo::blas::GUARD_PRIMAL_FLOPS.store(
        pf.map_or(0, |(f, _)| f.to_bits()),
        std::sync::atomic::Ordering::Relaxed,
    );
    let mut dual_basis = HighsBasis::default();
    let mut dual_solution = HighsSolution::default();
    let mut dual_info = HighsInfo::default();
    let mut dual_callback = HighsCallback::default();
    let mut dual_model_status = HighsModelStatus::kNotset;
    let mut options = solver_object.options_.clone();
    options.ipm_hipo_dualize = false;
    let status = solve_hipo(
        &options,
        &mut *solver_object.timer_,
        &dual_lp,
        &HighsHessian::default(),
        &mut dual_basis,
        &mut dual_solution,
        &mut dual_model_status,
        &mut dual_info,
        &mut dual_callback,
    );
    crate::ipm::hipo::blas::GUARD_PRIMAL_FLOPS.store(0, std::sync::atomic::Ordering::Relaxed);
    if status == HighsStatus::kError
        || !(dual_model_status == HighsModelStatus::kOptimal
            // An imprecise dual crossover (kUnknown) still leaves a basis;
            // `undualize` re-solves the primal from it by simplex.
            || dual_model_status == HighsModelStatus::kUnknown)
        || !dual_basis.valid
    {
        return false;
    }
    solver_object.highs_info_.ipm_iteration_count += dual_info.ipm_iteration_count;
    solver_object.highs_info_.crossover_iteration_count += dual_info.crossover_iteration_count;
    if ekk.set_basis_from(&dual_basis) != HighsStatus::kOk {
        return false;
    }
    if ekk.undualize() != HighsStatus::kOk || ekk.model_status_ != HighsModelStatus::kOptimal {
        return false;
    }
    let mut lp_copy = solver_object.lp_.clone();
    *solver_object.basis_ = ekk.get_highs_basis(&mut lp_copy);
    *solver_object.solution_ = ekk.get_solution();
    solver_object.model_status_ = HighsModelStatus::kOptimal;
    solver_object.highs_info_.simplex_iteration_count += ekk.iteration_count_;
    true
}

/// HiGHS: IpxWrapper.cpp::solveHipo (lines 415-657)
/// Also called by `Highs::callSolveQp` (Highs.cpp 4161-4168, `lp_data/highs_qp.rs`)
/// with the model Hessian. `assert(1 == 0)` in the stopped-status chain is
/// `debug_assert!(false)`. `Q.empty()` is `Q.dim_ == 0` (`HighsHessian::empty`).
#[allow(clippy::too_many_arguments)]
pub fn solve_hipo(
    options: &HighsOptions,
    timer: &mut HighsTimer,
    lp: &HighsLp,
    Q: &HighsHessian,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
    model_status: &mut HighsModelStatus,
    highs_info: &mut HighsInfo,
    callback: &mut HighsCallback,
) -> HighsStatus {
    // HiGHS: IpxWrapper.cpp::solveHipo line 454 (no-op: BLAS is serial)

    // Indicate that there is no valid primal solution, dual solution or basis
    highs_basis.valid = false;
    highs_solution.value_valid = false;
    highs_solution.dual_valid = false;
    // Indicate that no imprecise solution has (yet) been found
    reset_model_status_and_highs_info(model_status, highs_info);

    // Create solver instance
    let mut hipo = Solver::default();
    // This creates ipx::LpSolver ipx_lps_, in case HiPO has to switch
    // to IPX, so use the current HiGHS time as an offset for the
    // ipx_lps.control_ elapsed time
    hipo.setIpxTimerOffset(timer.read());

    if options.kkt_tolerance != kDefaultKktTolerance {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "IpxWrapper: feasibility_tol = {}; optimality_tol = {}; crossover_tol = {}\n",
            c_fmt_g(options.kkt_tolerance, 6),
            c_fmt_g(1e-1 * options.kkt_tolerance, 6),
            c_fmt_g(1e-1 * options.kkt_tolerance, 6)
        );
    }

    hipo.setOptions(options);
    hipo.setTimer(timer);
    hipo.setCallback(std::mem::take(callback));

    // Load the problem
    let load_status: Int = hipo.load(lp, Q);
    if load_status != 0 {
        *model_status = HighsModelStatus::kSolveError;
        *callback = hipo.take_callback().unwrap_or_default();
        return HighsStatus::kError;
    }

    // This information about the problem loaded into HiPO is needed for later
    let mut num_row: HighsInt = 0;
    let mut num_col: HighsInt = 0;
    hipo.getOriginalDims(&mut num_row, &mut num_col);
    let mut rhs: Vec<f64> = Vec::new();
    let mut constraints: Vec<u8> = Vec::new();
    fill_in_rhs_and_constraints(lp, &mut rhs, &mut constraints);

    hipo.solve();
    if hipo.guard_aborted_ {
        *callback = hipo.take_callback().unwrap_or_default();
        return HighsStatus::kError;
    }

    // const bool report_solve_data =
    //    kHighsAnalysisLevelSolverSummaryData & options.highs_analysis_level;

    // Differently from IPX, HiPO returns a single status. So, dealing with
    // statuses is a bit different.
    // hipo.solved(), hipo.stopped(), hipo.failed() can be used to query if the
    // status belongs to the solved, stopped or failed group.
    // If primal-dual feasible solution is found (non-vertex solution), then the
    // status is kStatusPDfeas.
    // If crossover is successful, then the status is kStatusBasic.
    // Otherwise, the specific crossover status can be accessed through the
    // ipx_info stored in hipo_info.

    // Get solver and solution information.
    // HiGHS: IpxWrapper.cpp::solveHipo line 505
    let hipo_info = *hipo.getInfo();
    oracle_trace::hipo_end(&hipo_info);
    let solve_status: Int = hipo_info.status;
    highs_info.ipm_iteration_count += hipo_info.ipm_iter + hipo_info.ipx_info.iter;
    highs_info.crossover_iteration_count += hipo_info.ipx_info.updates_crossover;

    let status: HighsStatus = 'body: {
        // Report hipo status
        let solve_return_status = report_hipo_status(options, solve_status, &hipo);
        if solve_return_status == HighsStatus::kError {
            *model_status = HighsModelStatus::kSolveError;
            break 'body HighsStatus::kError;
        }

        // Report crossover status
        let crossover_return_status =
            report_hipo_crossover_status(options, hipo_info.ipx_info.status_crossover, !Q.empty());
        if crossover_return_status == HighsStatus::kError {
            *model_status = HighsModelStatus::kSolveError;
            break 'body HighsStatus::kError;
        }

        // Failures should have been handled. Status should be stopped or solved.
        if ipx_status_error(
            !hipo.solved() && !hipo.stopped(),
            options,
            "Hipo",
            "status should be solved or stopped but value is",
            solve_status,
        ) {
            break 'body HighsStatus::kError;
        }

        if hipo.stopped() {
            let local_model_status = HighsModelStatus::kUnknown;

            get_hipo_non_vertex_solution(
                options,
                lp,
                num_col,
                num_row,
                &rhs,
                &constraints,
                &hipo,
                local_model_status,
                highs_solution,
            );

            // For crossover
            if illegal_ipx_stopped_crossover_status(&hipo_info.ipx_info, options) {
                break 'body HighsStatus::kError;
            }
            // Can stop and reach time limit
            if hipo_info.ipx_info.status_crossover == ipx::IPX_STATUS_time_limit {
                *model_status = HighsModelStatus::kTimeLimit;
                break 'body HighsStatus::kWarning;
            }

            // if crossover didn't time out, then solver can only stop as follows
            if solve_status == kStatusUserInterrupt {
                *model_status = HighsModelStatus::kInterrupt;
                break 'body HighsStatus::kWarning;
            } else if solve_status == kStatusTimeLimit {
                *model_status = HighsModelStatus::kTimeLimit;
                break 'body HighsStatus::kWarning;
            } else if solve_status == kStatusMaxIter {
                *model_status = HighsModelStatus::kIterationLimit;
                break 'body HighsStatus::kWarning;
            } else if solve_status == kStatusNoProgress {
                report_hipo_no_progress(options, &hipo_info);
                *model_status = HighsModelStatus::kUnknown;
                break 'body HighsStatus::kWarning;
            } else {
                debug_assert!(false);
            }
        }

        // Stopped status should have been handled. Status should be solved.
        if ipx_status_error(
            !hipo.solved(),
            options,
            "Hipo",
            "status should be solved but value is",
            solve_status,
        ) {
            break 'body HighsStatus::kError;
        }

        // primal/dual infeasible
        if solve_status == kStatusPrimalInfeasible || solve_status == kStatusDualInfeasible {
            if solve_status == kStatusPrimalInfeasible {
                *model_status = HighsModelStatus::kInfeasible;
            } else {
                *model_status = HighsModelStatus::kUnboundedOrInfeasible;
            }

            get_hipo_non_vertex_solution(
                options,
                lp,
                num_col,
                num_row,
                &rhs,
                &constraints,
                &hipo,
                *model_status,
                highs_solution,
            );

            break 'body HighsStatus::kOk;
        }

        // Status should be optimal or imprecise
        if ipx_status_error(
            solve_status != kStatusPDFeas
                && solve_status != kStatusBasic
                && solve_status != kStatusImprecise,
            options,
            "Hipo",
            "status should be optimal or imprecise but value is",
            solve_status,
        ) {
            break 'body HighsStatus::kError;
        }

        let have_basic_solution =
            hipo_info.ipx_used && hipo_info.ipx_info.status_crossover != ipx::IPX_STATUS_not_run;

        let imprecise_solution = hipo_info.status == kStatusImprecise
            || hipo_info.ipx_info.status_crossover == ipx::IPX_STATUS_imprecise;

        if have_basic_solution {
            let mut hipo_solution = IpxSolution::default();
            hipo_solution.num_col = num_col;
            hipo_solution.num_row = num_row;
            hipo_solution.ipx_col_value.resize(num_col as usize, 0.0);
            hipo_solution.ipx_row_value.resize(num_row as usize, 0.0);
            hipo_solution.ipx_col_dual.resize(num_col as usize, 0.0);
            hipo_solution.ipx_row_dual.resize(num_row as usize, 0.0);
            hipo_solution.ipx_row_status.resize(num_row as usize, 0);
            hipo_solution.ipx_col_status.resize(num_col as usize, 0);

            let errflag: Int = hipo.getBasicSolution(
                &mut hipo_solution.ipx_col_value,
                &mut hipo_solution.ipx_row_value,
                &mut hipo_solution.ipx_row_dual,
                &mut hipo_solution.ipx_col_dual,
                &mut hipo_solution.ipx_row_status,
                &mut hipo_solution.ipx_col_status,
            );

            if errflag != 0 {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "IPX crossover getting basic solution: flag = {}\n",
                    errflag
                );
                break 'body HighsStatus::kError;
            }
            // Convert the IPX basic solution to a HiGHS basic solution
            let status = ipx_basic_solution_to_highs_basic_solution(
                &options.log_options,
                lp,
                &rhs,
                &constraints,
                &hipo_solution,
                highs_basis,
                highs_solution,
            );
            if status != HighsStatus::kOk {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "Failed to convert IPX basic solution to Highs basic solution\n"
                );
                break 'body HighsStatus::kError;
            }
        } else {
            let local_model_status = if imprecise_solution {
                HighsModelStatus::kUnknown
            } else {
                HighsModelStatus::kOptimal
            };
            get_hipo_non_vertex_solution(
                options,
                lp,
                num_col,
                num_row,
                &rhs,
                &constraints,
                &hipo,
                local_model_status,
                highs_solution,
            );
            debug_assert!(!highs_basis.valid);
        }

        highs_info.basis_validity = if highs_basis.valid {
            kBasisValidityValid
        } else {
            kBasisValidityInvalid
        };

        if imprecise_solution {
            *model_status = HighsModelStatus::kUnknown;
            HighsStatus::kWarning
        } else {
            *model_status = HighsModelStatus::kOptimal;
            HighsStatus::kOk
        }
    };

    *callback = hipo.take_callback().unwrap_or_default();
    status
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn hipo_on_afiro_like_lp_is_optimal() {
        // min -x0 - x1 s.t. x0 + x1 <= 1, 0 <= x <= 1
        let lp = HighsLp {
            num_col_: 2,
            num_row_: 1,
            col_cost_: vec![-1.0, -1.0],
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![1.0, 1.0],
            row_lower_: vec![-kHighsInf],
            row_upper_: vec![1.0],
            a_matrix_: HighsSparseMatrix {
                format_: MatrixFormat::kColwise,
                num_col_: 2,
                num_row_: 1,
                start_: vec![0, 1, 2],
                index_: vec![0, 0],
                value_: vec![1.0, 1.0],
                ..HighsSparseMatrix::default()
            },
            ..HighsLp::default()
        };
        let options = HighsOptions {
            output_flag: false,
            hipo_available: true,
            ..HighsOptions::default()
        };
        let mut timer = HighsTimer::default();
        let mut basis = HighsBasis::default();
        let mut solution = HighsSolution::default();
        let mut model_status = HighsModelStatus::kNotset;
        let mut info = HighsInfo::default();
        let mut callback = HighsCallback::default();
        let status = solve_hipo(
            &options,
            &mut timer,
            &lp,
            &HighsHessian::default(),
            &mut basis,
            &mut solution,
            &mut model_status,
            &mut info,
            &mut callback,
        );
        assert_ne!(status, HighsStatus::kError);
        assert!(solution.value_valid);
        assert!((solution.col_value[0] + solution.col_value[1] - 1.0).abs() < 1e-6);
    }
}
