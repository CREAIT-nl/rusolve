//! HiGHS source: ipm/IpxWrapper.cpp (IPX part), ipm/IpxWrapper.h, ipm/IpxSolution.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiGHS LP -> IPX model translation, the IPX `Info` -> HiGHS status
//! mapping and the IPX solution -> HiGHS solution copy. The HiPO functions
//! of the same file (`solveLpHipo`, `solveHipo`, `reportHipoStatus`,
//! `reportHipoCrossoverStatus`, `reportHipoNoProgress`,
//! `getHipoNonVertexSolution`, lines 408-657, 1131-1142, 1170-1188,
//! 1393-1511) are in [`crate::ipm::hipo_wrapper`].
//!
//! # Call surface
//!
//! * `std::vector::data()` of an empty libc++ vector is `nullptr`, and
//!   `ipx::Model::CopyInput` rejects null pointers. [`solve_lp_ipx`] passes
//!   [`vec_data`] of each vector to `LpSolver::load_model`, so an empty
//!   vector reaches IPX as `None` exactly as in C++.
//! * `lps.SetCallback(&callback)` moves the `HighsCallback` into the
//!   `LpSolver` ([`LpSolver::set_callback`]); move it back with
//!   [`LpSolver::take_callback`] on every return path after it was set.
//! * C++ `throw` inside IPX is a Rust `panic!`; `solveLp`'s
//!   `catch (const std::exception&)` is not reproduced .

use crate::ipm::ipx::lp_solver::LpSolver;
use crate::ipm::ipx::{
    self, IPX_ERROR_argument_null, IPX_ERROR_invalid_basis, IPX_ERROR_invalid_dimension,
    IPX_ERROR_invalid_matrix, IPX_ERROR_invalid_vector, IPX_STATUS_debug, IPX_STATUS_dual_infeas,
    IPX_STATUS_failed, IPX_STATUS_imprecise, IPX_STATUS_internal_error, IPX_STATUS_iter_limit,
    IPX_STATUS_no_model, IPX_STATUS_no_progress, IPX_STATUS_not_run, IPX_STATUS_optimal,
    IPX_STATUS_out_of_memory, IPX_STATUS_primal_infeas, IPX_STATUS_solved, IPX_STATUS_stopped,
    IPX_STATUS_time_limit, IPX_STATUS_user_interrupt, Int,
};
use crate::lp_data::highs_solution::reset_model_status_and_highs_info;
use crate::lp_data::highs_solution_ipx::ipx_basic_solution_to_highs_basic_solution;
use crate::prelude::*;
use crate::util::c_printf::{c_fmt_f_width, c_fmt_g_width};
use std::io::Write;

/// HiGHS: ipm/IpxSolution.h::IpxSolution
#[derive(Debug, Clone, Default)]
pub struct IpxSolution {
    pub num_col: Int,
    pub num_row: Int,
    pub ipx_col_value: Vec<f64>,
    pub ipx_row_value: Vec<f64>,
    pub ipx_col_dual: Vec<f64>,
    pub ipx_row_dual: Vec<f64>,
    pub ipx_col_status: Vec<Int>,
    pub ipx_row_status: Vec<Int>,
}

/// HiGHS 1.16 #3456 (h16_ac_centring): the `centring_status` value of IpxWrapper.cpp::solveLpIpx
/// lines 205-208 (1.16.0): `!centring_tried ? kCentringStatusNone :
/// centring_success ? kCentringStatusSuccess : kCentringStatusFailed`.
/// Only called from the `h16_ac_centring` branch of `solve_lp_ipx`.
pub fn h16_centring_status(centring_tried: Int, centring_success: Int) -> HighsInt {
    if centring_tried == 0 {
        kCentringStatusNone
    } else if centring_success != 0 {
        kCentringStatusSuccess
    } else {
        kCentringStatusFailed
    }
}

/// libc++ `std::vector<T>::data()` as an IPX input pointer: `None` (null)
/// for an empty vector, the whole vector otherwise.
#[inline]
pub fn vec_data<T>(v: &[T]) -> Option<&[T]> {
    if v.is_empty() { None } else { Some(v) }
}

/// libc++ `std::vector<T>::data()` as an IPX output pointer: `None` (null)
/// for an empty vector, the whole vector otherwise.
#[inline]
pub fn vec_data_mut<T>(v: &mut [T]) -> Option<&mut [T]> {
    if v.is_empty() { None } else { Some(v) }
}

/// HiGHS: IpxWrapper.cpp::solveLpIpx(HighsLpSolverObject&) (lines 22-27)
pub fn solve_lp_ipx_so(solver_object: &mut HighsLpSolverObject) -> HighsStatus {
    solve_lp_ipx(
        &*solver_object.options_,
        &mut *solver_object.timer_,
        &*solver_object.lp_,
        &mut *solver_object.basis_,
        &mut *solver_object.solution_,
        &mut solver_object.model_status_,
        &mut *solver_object.highs_info_,
        &mut *solver_object.callback_,
    )
}

/// HiGHS: IpxWrapper.cpp::solveLpIpx (lines 29-406)
///
/// PORT NOTE: the C++ `return`s after `lps.SetCallback(&callback)` are kept
/// as `return`s inside a closure, so that `take_callback` runs on every path.
/// The closure is only a control-flow frame:
/// it is called once, in place.
#[allow(clippy::redundant_closure_call)]
pub fn solve_lp_ipx(
    options: &HighsOptions,
    timer: &mut HighsTimer,
    lp: &HighsLp,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
    model_status: &mut HighsModelStatus,
    highs_info: &mut HighsInfo,
    callback: &mut HighsCallback,
) -> HighsStatus {
    // Use IPX to try to solve the LP
    //
    // Can return HighsModelStatus (HighsStatus) values:
    //
    // 1. kSolveError (kError) if various unlikely solution errors occur
    //
    // 2. kTimeLimit (kWarning) if time limit is reached
    //
    // 3. kIterationLimit (kWarning) if iteration limit is reached
    //
    // 4. kUnknown (kWarning) if IPM makes no progress or if
    // IPM/crossover are imprecise
    //
    // 5. kInfeasible (kOk) if IPM identifies primal infeasibility
    //
    // 6. kUnboundedOrInfeasible (kOk) if IPM identifies dual
    // infeasibility
    //
    // kOptimal (kOk) if IPM/crossover identify optimality
    //
    // With a non-error return, if just IPM has been run then a
    // non-vertex primal solution is obtained; if crossover has been run
    // then a basis and primal+dual solution are obtained.
    //
    //
    // Indicate that there is no valid primal solution, dual solution or basis
    highs_basis.valid = false;
    highs_solution.value_valid = false;
    highs_solution.dual_valid = false;
    // Indicate that no imprecise solution has (yet) been found
    reset_model_status_and_highs_info(model_status, highs_info);
    // Create the LpSolver instance
    let mut lps = LpSolver::new();
    // Use the current HiGHS time as an offset for the lps.control_
    // elapsed time
    lps.set_timer_offset(timer.read());

    // Set IPX parameters
    //
    // Cannot set internal IPX parameters directly since they are
    // private, so create instance of parameters
    let mut parameters = ipx::Parameters::default();
    // Set IPX parameters from options
    //
    // Set display according to output
    parameters.display = 1;
    if !options.output_flag || !options.log_to_console {
        parameters.display = 0;
    }
    // Modify parameters.debug according to log_dev_level
    parameters.debug = 0;
    if options.log_dev_level == kHighsLogDevLevelDetailed {
        // Default options.log_dev_level setting is kHighsLogDevLevelNone, yielding
        // default setting debug = 0
        parameters.debug = 0;
    } else if options.log_dev_level == kHighsLogDevLevelInfo {
        parameters.debug = 2;
    } else if options.log_dev_level == kHighsLogDevLevelVerbose {
        parameters.debug = 4;
    }
    parameters.highs_logging = true;
    parameters.timeless_log = options.timeless_log;
    // Use this to change print_interval from default value (5s) to
    // other value for debugging
    //
    //  parameters.print_interval = 1e-4;
    parameters.log_options = Some(options.log_options.clone());
    // Just test feasibility and optimality tolerances for now
    // ToDo Set more parameters
    //
    // Translate dualization option
    //
    // parameters.dualize = -2 => Possibly dualize - Filippo style
    // parameters.dualize = -1 => Possibly dualize - Lukas style
    // parameters.dualize = 0 => No dualization
    // parameters.dualize = 1 => Perform dualization
    if options.ipx_dualize_strategy == kIpxDualizeStrategyOn {
        parameters.dualize = 1;
    } else if options.ipx_dualize_strategy == kIpxDualizeStrategyOff {
        parameters.dualize = 0;
    } else if options.ipx_dualize_strategy == kIpxDualizeStrategyLukas {
        parameters.dualize = -1;
    } else if options.ipx_dualize_strategy == kIpxDualizeStrategyFilippo {
        parameters.dualize = -2;
    } else {
        // C++ `assert(111 == 222)`. Reached also for kIpxDualizeStrategyChoose
        // (0), which keeps the Parameters default dualize = -1 in release builds.
        debug_assert!(
            false,
            "ipx_dualize_strategy {} not handled",
            options.ipx_dualize_strategy
        );
    }

    parameters.ipm_feasibility_tol = cpp_min(
        options.primal_feasibility_tolerance,
        options.dual_feasibility_tolerance,
    );
    parameters.ipm_optimality_tol = options.ipm_optimality_tolerance;
    parameters.start_crossover_tol = options.start_crossover_tolerance;

    if options.kkt_tolerance != kDefaultKktTolerance {
        parameters.ipm_feasibility_tol = options.kkt_tolerance;
        parameters.ipm_optimality_tol = 1e-1 * options.kkt_tolerance;
        parameters.start_crossover_tol = 1e-1 * options.kkt_tolerance;
    }

    parameters.analyse_basis_data =
        (kHighsAnalysisLevelNlaData & options.highs_analysis_level) != 0;
    // Now that the lps.control_ elapsed time includes the HiGHS time,
    // can use the HiGHS time limit
    parameters.time_limit = options.time_limit;
    parameters.ipm_maxiter = options.ipm_iteration_limit - highs_info.ipm_iteration_count;
    // Determine if crossover is to be run or not
    //
    // When doing analytic centring calculations, crossover must not be
    // run
    if options.run_centring {
        parameters.run_crossover = 0;
    } else if options.run_crossover == kHighsOnString {
        parameters.run_crossover = 1;
    } else if options.run_crossover == kHighsOffString {
        parameters.run_crossover = 0;
    } else {
        debug_assert!(options.run_crossover == kHighsChooseString);
        parameters.run_crossover = -1;
    }
    if parameters.run_crossover == 0 {
        // If crossover is sure not to be run, then set crossover_start_ to
        // -1 so that IPX can terminate according to its feasibility and
        // optimality tolerances
        parameters.start_crossover_tol = -1.0;
    }

    parameters.run_centring = if options.run_centring { 1 } else { 0 };
    parameters.max_centring_steps = options.max_centring_steps;
    parameters.centring_ratio_tolerance = options.centring_ratio_tolerance;
    // rusolve (HiGHS 1.16 backport, default off): copied into the IPX parameters
    parameters.h16_ipx_stall = options.h16_ipx_stall;
    parameters.h16_ipx_scaled = options.h16_ipx_scaled;
    parameters.h16_ipm_gap_offset = options.h16_ipm_gap_offset;
    parameters.stall_to_crossover = options.lp_ipm_failure_policy;
    parameters.crossover_ipm_crash = options.lp_crossover_ipm_crash;
    parameters.crossover_batch_push = options.lp_crossover_batch_push;
    if options.lp_ipx_no_drop {
        parameters.ipm_drop_primal = 0.0;
        parameters.ipm_drop_dual = 0.0;
    }
    // LP path only: the MIP analytic-centre solves (run_centring) keep HiGHS's drop rule so the
    // MIP search and its oracle replay are unchanged
    parameters.ipm_drop_dual_soft = options.lp_ipx_soft_drop && !options.run_centring;

    // Set the internal IPX parameters
    lps.set_parameters(parameters);

    // Set pointer to any callback. The callback is moved into `lps` here and
    // moved back below, after the closure has returned on any path.
    lps.set_callback(Some(std::mem::take(callback)));

    let status = (|| -> HighsStatus {
        let mut num_col: Int = 0;
        let mut num_row: Int = 0;
        let mut offset: f64 = 0.0;
        let mut objective: Vec<f64> = Vec::new();
        let mut col_lb: Vec<f64> = Vec::new();
        let mut col_ub: Vec<f64> = Vec::new();
        let mut Ap: Vec<Int> = Vec::new();
        let mut Ai: Vec<Int> = Vec::new();
        let mut Av: Vec<f64> = Vec::new();
        let mut rhs: Vec<f64> = Vec::new();
        let mut constraint_type: Vec<u8> = Vec::new();
        fill_in_ipx_data(
            lp,
            &mut num_col,
            &mut num_row,
            &mut offset,
            &mut objective,
            &mut col_lb,
            &mut col_ub,
            &mut Ap,
            &mut Ai,
            &mut Av,
            &mut rhs,
            &mut constraint_type,
            options.h16_ipx_free_rows,
        );
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "IPX model has {} rows, {} columns and {} nonzeros\n",
            num_row,
            num_col,
            Ap[num_col as usize]
        );

        let load_status = lps.load_model(
            num_col,
            offset,
            vec_data(&objective),
            vec_data(&col_lb),
            vec_data(&col_ub),
            num_row,
            vec_data(&Ap),
            vec_data(&Ai),
            vec_data(&Av),
            vec_data(&rhs),
            vec_data(&constraint_type),
        );

        if load_status != 0 {
            *model_status = HighsModelStatus::kSolveError;
            return HighsStatus::kError;
        }

        // Use IPX to solve the LP!
        let solve_status = lps.solve();

        let want_solve_data =
            (kHighsAnalysisLevelSolverSummaryData & options.highs_analysis_level) != 0;
        // Get solver and solution information.
        // Struct ipx_info defined in ipx/ipx_info.h
        let ipx_info = lps.get_info();
        if want_solve_data {
            report_solve_data(&options.log_options, &ipx_info);
        }
        highs_info.ipm_iteration_count += ipx_info.iter;
        highs_info.crossover_iteration_count += ipx_info.updates_crossover;
        // HiGHS 1.16 #3456 (h16_ac_centring): IpxWrapper.cpp::solveLpIpx lines 201-206 (1.16.0).
        // Records the outcome of any centring, since the model status only indicates whether
        // the LP is solved. Centring is not tried if the IPM stops before reaching optimality.
        if options.h16_ac_centring {
            highs_info.centring_status =
                h16_centring_status(ipx_info.centring_tried, ipx_info.centring_success);
        }

        // If not solved...
        if solve_status != IPX_STATUS_solved {
            let solve_return_status =
                report_ipx_solve_status(options, solve_status, ipx_info.errflag);
            // Return error if IPX solve error has occurred
            if solve_return_status == HighsStatus::kError {
                *model_status = HighsModelStatus::kSolveError;
                return HighsStatus::kError;
            }
        }
        let mut ipm_status = true;
        let ipm_return_status =
            report_ipx_ipm_crossover_status(options, ipx_info.status_ipm, ipm_status);
        ipm_status = false;
        let crossover_return_status =
            report_ipx_ipm_crossover_status(options, ipx_info.status_crossover, ipm_status);
        // Return error if IPX IPM or crossover error has occurred
        if ipm_return_status == HighsStatus::kError
            || crossover_return_status == HighsStatus::kError
        {
            *model_status = HighsModelStatus::kSolveError;
            return HighsStatus::kError;
        }
        // Should only reach here if Solve() returned IPX_STATUS_solved or
        // IPX_STATUS_stopped
        if ipx_status_error(
            solve_status != IPX_STATUS_solved && solve_status != IPX_STATUS_stopped,
            options,
            "Ipx",
            "solve_status should be solved or stopped here but value is",
            solve_status,
        ) {
            return HighsStatus::kError;
        }

        // Only error returns so far
        //

        if solve_status == IPX_STATUS_stopped {
            // IPX stopped, so there's certainly no basic solution. Get the
            // non-vertex solution, though. This needs the model status to
            // know whether to worry about dual infeasibilities.
            let local_model_status = HighsModelStatus::kUnknown;
            get_highs_non_vertex_solution(
                options,
                lp,
                num_col,
                num_row,
                &rhs,
                &constraint_type,
                &lps,
                local_model_status,
                highs_solution,
            );
            //
            // Look at the reason why IPX stopped
            //
            // Return error if stopped status settings occur that JAJH doesn't
            // think should happen
            //
            //==============
            // For crossover
            //==============
            if illegal_ipx_stopped_crossover_status(&ipx_info, options) {
                return HighsStatus::kError;
            }
            // Can stop and reach time limit
            if ipx_info.status_crossover == IPX_STATUS_time_limit {
                *model_status = HighsModelStatus::kTimeLimit;
                return HighsStatus::kWarning;
            }
            //========
            // For IPM
            //========
            //
            // Note that IPX can stop with IPM optimal, imprecise,
            // primal_infeas or dual_infeas, due to crossover stopping with
            // time limit, and this is why crossover returns are tested first
            if illegal_ipx_stopped_ipm_status(&ipx_info, options) {
                return HighsStatus::kError;
            }
            // Can stop with user interrupt
            // Can stop with time limit
            // Can stop with iter limit
            // Can stop with no progress
            if ipx_info.status_ipm == IPX_STATUS_user_interrupt {
                *model_status = HighsModelStatus::kInterrupt;
                return HighsStatus::kWarning;
            } else if ipx_info.status_ipm == IPX_STATUS_time_limit {
                *model_status = HighsModelStatus::kTimeLimit;
                return HighsStatus::kWarning;
            } else if ipx_info.status_ipm == IPX_STATUS_iter_limit {
                *model_status = HighsModelStatus::kIterationLimit;
                return HighsStatus::kWarning;
            } else {
                debug_assert!(ipx_info.status_ipm == IPX_STATUS_no_progress);
                report_ipm_no_progress(options, &ipx_info);
                *model_status = HighsModelStatus::kUnknown;
                return HighsStatus::kWarning;
            }
        }
        // Should only reach here if Solve() returned IPX_STATUS_solved
        if ipx_status_error(
            solve_status != IPX_STATUS_solved,
            options,
            "Ipx",
            "solve_status should be solved here but value is",
            solve_status,
        ) {
            return HighsStatus::kError;
        }
        // Return error if solved status settings occur that JAJH doesn't
        // think should happen
        if illegal_ipx_solved_status(&ipx_info, options) {
            return HighsStatus::kError;
        }
        //==============
        // For crossover
        //==============
        // Can be not run
        // Can solve and be optimal
        // Can solve and be imprecise
        //========
        // For IPM
        //========
        // Can solve and be optimal
        // Can solve and be imprecise
        // Can solve and be primal_infeas
        // Can solve and be dual_infeas
        if ipx_info.status_ipm == IPX_STATUS_primal_infeas
            || ipx_info.status_ipm == IPX_STATUS_dual_infeas
        {
            // IPM identified primal or dual infeasibility: crossover will not
            // have run, so get the non-vertex solution and return
            if ipx_info.status_ipm == IPX_STATUS_primal_infeas {
                *model_status = HighsModelStatus::kInfeasible;
            } else if ipx_info.status_ipm == IPX_STATUS_dual_infeas {
                *model_status = HighsModelStatus::kUnboundedOrInfeasible;
            }
            get_highs_non_vertex_solution(
                options,
                lp,
                num_col,
                num_row,
                &rhs,
                &constraint_type,
                &lps,
                *model_status,
                highs_solution,
            );
            return HighsStatus::kOk;
        }

        // Should only reach here if IPM is optimal or imprecise
        if ipx_status_error(
            ipx_info.status_ipm != IPX_STATUS_optimal
                && ipx_info.status_ipm != IPX_STATUS_imprecise,
            options,
            "Ipx",
            "ipm status should be not run, optimal or imprecise but value is",
            ipx_info.status_ipm,
        ) {
            return HighsStatus::kError;
        }

        // Should only reach here if crossover is not run, optimal or imprecise
        if ipx_status_error(
            ipx_info.status_crossover != IPX_STATUS_not_run
                && ipx_info.status_crossover != IPX_STATUS_optimal
                && ipx_info.status_crossover != IPX_STATUS_imprecise,
            options,
            "Ipx",
            "crossover status should be not run, optimal or imprecise but value is",
            ipx_info.status_crossover,
        ) {
            return HighsStatus::kError;
        }

        // Basic solution depends on crossover being run
        let have_basic_solution = ipx_info.status_crossover != IPX_STATUS_not_run;
        // Both crossover and IPM can be imprecise
        let imprecise_solution = ipx_info.status_crossover == IPX_STATUS_imprecise
            || ipx_info.status_ipm == IPX_STATUS_imprecise;
        if have_basic_solution {
            let mut ipx_solution = IpxSolution {
                num_col,
                num_row,
                ..Default::default()
            };
            ipx_solution.ipx_col_value.resize(num_col as usize, 0.0);
            ipx_solution.ipx_row_value.resize(num_row as usize, 0.0);
            ipx_solution.ipx_col_dual.resize(num_col as usize, 0.0);
            ipx_solution.ipx_row_dual.resize(num_row as usize, 0.0);
            ipx_solution.ipx_row_status.resize(num_row as usize, 0);
            ipx_solution.ipx_col_status.resize(num_col as usize, 0);
            let errflag = lps.get_basic_solution(
                vec_data_mut(&mut ipx_solution.ipx_col_value),
                vec_data_mut(&mut ipx_solution.ipx_row_value),
                vec_data_mut(&mut ipx_solution.ipx_row_dual),
                vec_data_mut(&mut ipx_solution.ipx_col_dual),
                vec_data_mut(&mut ipx_solution.ipx_row_status),
                vec_data_mut(&mut ipx_solution.ipx_col_status),
            );
            if errflag != 0 {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "IPX crossover getting basic solution: flag = {}\n",
                    errflag
                );
                return HighsStatus::kError;
            }
            // Convert the IPX basic solution to a HiGHS basic solution
            let status = ipx_basic_solution_to_highs_basic_solution(
                &options.log_options,
                lp,
                &rhs,
                &constraint_type,
                &ipx_solution,
                highs_basis,
                highs_solution,
            );
            if status != HighsStatus::kOk {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "Failed to convert IPX basic solution to Highs basic solution\n"
                );
                return HighsStatus::kError;
            }
        } else {
            // No basic solution, so get a non-vertex HiGHS solution. This
            // needs the model status to know whether to worry about dual
            // infeasibilities.
            let local_model_status = if imprecise_solution {
                HighsModelStatus::kUnknown
            } else {
                HighsModelStatus::kOptimal
            };
            get_highs_non_vertex_solution(
                options,
                lp,
                num_col,
                num_row,
                &rhs,
                &constraint_type,
                &lps,
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
    })();

    // The callback goes back to the caller on every path (see above).
    *callback = lps
        .take_callback()
        .expect("IPX callback was moved in by SetCallback");
    status
}

/// HiGHS: IpxWrapper.cpp::fillInIpxData (lines 659-773)
/// HiGHS 1.16 (#3319): IpxWrapper.cpp::fillInIpxData (lines 684-799)
///
/// `constraint_type` is C++ `std::vector<char>`.
/// `h16_ipx_free_rows` is the HiGHS 1.16 #3319 free-row flag (`options.h16_ipx_free_rows`):
/// off gives the 1.15.1 slack row index, on maps it through `reduced_rowmap`.
pub fn fill_in_ipx_data(
    lp: &HighsLp,
    num_col: &mut Int,
    num_row: &mut Int,
    offset: &mut f64,
    obj: &mut Vec<f64>,
    col_lb: &mut Vec<f64>,
    col_ub: &mut Vec<f64>,
    Ap: &mut Vec<Int>,
    Ai: &mut Vec<Int>,
    Ax: &mut Vec<f64>,
    rhs: &mut Vec<f64>,
    constraint_type: &mut Vec<u8>,
    // rusolve (HiGHS 1.16 backport, default off): the #3319 free-row flag (gates below)
    h16_ipx_free_rows: bool,
) {
    *num_col = lp.num_col_;
    *num_row = lp.num_row_;

    // For each row with both a lower and an upper bound introduce one new column
    // so num_col may increase. Ignore each free row so num_row may decrease.
    // lba <= a'x <= uba becomes
    // a'x-s = 0 and lba <= s <= uba.

    // For each row with bounds on both sides introduce explicit slack and
    // transfer bounds.
    debug_assert!(lp.row_lower_.len() as Int == *num_row);
    debug_assert!(lp.row_upper_.len() as Int == *num_row);

    let mut general_bounded_rows: Vec<Int> = Vec::new();
    let mut free_rows: Vec<Int> = Vec::new();

    for row in 0..*num_row {
        let r = row as usize;
        if lp.row_lower_[r] < lp.row_upper_[r]
            && lp.row_lower_[r] > -kHighsInf
            && lp.row_upper_[r] < kHighsInf
        {
            general_bounded_rows.push(row);
        } else if lp.row_lower_[r] <= -kHighsInf && lp.row_upper_[r] >= kHighsInf {
            free_rows.push(row);
        }
    }

    let num_slack = general_bounded_rows.len() as Int;

    fill_in_rhs_and_constraints(lp, rhs, constraint_type);

    let mut reduced_rowmap: Vec<Int> = vec![-1; lp.num_row_ as usize];
    if !free_rows.is_empty() {
        let mut counter: Int = 0;
        let mut findex: usize = 0;
        for row in 0..lp.num_row_ {
            // IPX quirk: C++ `free_rows[findex]` reads one element past the end
            // once every free row has been passed (undefined behaviour in C++).
            // `get` treats that read as "no match": free_rows is only ever
            // pushed to, so the slot past its end is never written.
            // HiGHS 1.16 #3319 (h16_ipx_free_rows): IpxWrapper.cpp lines 722-724 (1.16.0)
            // add the bound `findex < free_rows_size`, which `get` already gives.
            if free_rows.get(findex).copied() == Some(row) {
                findex += 1;
                continue;
            } else {
                reduced_rowmap[row as usize] = counter;
                counter += 1;
            }
        }
    } else {
        for (k, entry) in reduced_rowmap.iter_mut().enumerate() {
            *entry = k as Int;
        }
    }
    *num_row -= free_rows.len() as Int;
    *num_col += num_slack;

    let mut sizes: Vec<Int> = vec![0; *num_col as usize];

    for col in 0..lp.num_col_ as usize {
        for k in lp.a_matrix_.start_[col]..lp.a_matrix_.start_[col + 1] {
            let row = lp.a_matrix_.index_[k as usize] as usize;
            if lp.row_lower_[row] > -kHighsInf || lp.row_upper_[row] < kHighsInf {
                sizes[col] += 1;
            }
        }
    }
    // Copy Astart and Aindex to ipx::Int array.
    let nnz = lp.a_matrix_.index_.len() as Int;
    Ap.resize(*num_col as usize + 1, 0);
    Ai.reserve((nnz + num_slack) as usize);
    Ax.reserve((nnz + num_slack) as usize);

    // Set starting points of original and newly introduced columns.
    Ap[0] = 0;
    for col in 0..lp.num_col_ as usize {
        Ap[col + 1] = Ap[col] + sizes[col];
    }
    for col in lp.num_col_ as usize..*num_col as usize {
        Ap[col + 1] = Ap[col] + 1;
    }
    for k in 0..nnz as usize {
        let row = lp.a_matrix_.index_[k] as usize;
        if lp.row_lower_[row] > -kHighsInf || lp.row_upper_[row] < kHighsInf {
            Ai.push(reduced_rowmap[lp.a_matrix_.index_[k] as usize]);
            Ax.push(lp.a_matrix_.value_[k]);
        }
    }

    for k in 0..num_slack as usize {
        // HiGHS 1.16 #3319 (h16_ipx_free_rows): IpxWrapper.cpp line 770 (1.16.0)
        if h16_ipx_free_rows {
            Ai.push(reduced_rowmap[general_bounded_rows[k] as usize]);
        } else {
            Ai.push(general_bounded_rows[k] as Int);
        }
        Ax.push(-1.0);
    }

    // Column bound vectors.
    col_lb.resize(*num_col as usize, 0.0);
    col_ub.resize(*num_col as usize, 0.0);
    for col in 0..lp.num_col_ as usize {
        if lp.col_lower_[col] <= -kHighsInf {
            col_lb[col] = f64::NEG_INFINITY;
        } else {
            col_lb[col] = lp.col_lower_[col];
        }

        if lp.col_upper_[col] >= kHighsInf {
            col_ub[col] = f64::INFINITY;
        } else {
            col_ub[col] = lp.col_upper_[col];
        }
    }
    for slack in 0..num_slack as usize {
        let row = general_bounded_rows[slack] as usize;
        col_lb[lp.num_col_ as usize + slack] = lp.row_lower_[row];
        col_ub[lp.num_col_ as usize + slack] = lp.row_upper_[row];
    }

    *offset = (lp.sense_ as Int) as f64 * lp.offset_;
    obj.resize(*num_col as usize, 0.0);
    for col in 0..lp.num_col_ as usize {
        obj[col] = (lp.sense_ as Int) as f64 * lp.col_cost_[col];
    }
}

/// HiGHS: IpxWrapper.cpp::fillInRhsAndConstraints (lines 775-803)
pub fn fill_in_rhs_and_constraints(
    lp: &HighsLp,
    rhs: &mut Vec<f64>,
    constraint_type: &mut Vec<u8>,
) {
    // For each row except free rows add entry to char array and set up rhs
    // vector

    let num_row = lp.num_row_;

    rhs.reserve(num_row as usize);
    constraint_type.reserve(num_row as usize);

    for row in 0..num_row as usize {
        let row_lower = lp.row_lower_[row];
        let row_upper = lp.row_upper_[row];
        if row_lower > -kHighsInf && row_upper >= kHighsInf {
            rhs.push(row_lower);
            constraint_type.push(b'>');
        } else if row_lower <= -kHighsInf && row_upper < kHighsInf {
            rhs.push(row_upper);
            constraint_type.push(b'<');
        } else if row_lower == row_upper {
            rhs.push(row_upper);
            constraint_type.push(b'=');
        } else if row_lower > -kHighsInf && row_upper < kHighsInf {
            // general bounded
            rhs.push(0.0);
            constraint_type.push(b'=');
        }
    }
}

/// HiGHS: IpxWrapper.cpp::reportIpxSolveStatus (lines 805-847)
pub fn report_ipx_solve_status(
    options: &HighsOptions,
    solve_status: Int,
    error_flag: Int,
) -> HighsStatus {
    if solve_status == IPX_STATUS_solved {
        crate::highs_log_user!(&options.log_options, HighsLogType::kInfo, "Ipx: Solved\n");
        HighsStatus::kOk
    } else if solve_status == IPX_STATUS_stopped {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Ipx: Stopped\n"
        );
        HighsStatus::kWarning
        // Remaining cases are errors so drop through to return HighsStatus::kError;
    } else {
        if solve_status == IPX_STATUS_no_model {
            if error_flag == IPX_ERROR_argument_null {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "Ipx: Invalid input - argument_null\n"
                );
            } else if error_flag == IPX_ERROR_invalid_dimension {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "Ipx: Invalid input - invalid dimension\n"
                );
            } else if error_flag == IPX_ERROR_invalid_matrix {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "Ipx: Invalid input - invalid matrix\n"
                );
            } else if error_flag == IPX_ERROR_invalid_vector {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "Ipx: Invalid input - invalid vector\n"
                );
            } else if error_flag == IPX_ERROR_invalid_basis {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "Ipx: Invalid input - invalid basis\n"
                );
            } else {
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kError,
                    "Ipx: Invalid input - unrecognised error\n"
                );
            }
        } else if solve_status == IPX_STATUS_out_of_memory {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "Ipx: Out of memory\n"
            );
        } else if solve_status == IPX_STATUS_internal_error {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "Ipx: Internal error {}\n",
                error_flag
            );
        } else {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "Ipx: unrecognised solve status = {}\n",
                solve_status
            );
        }
        HighsStatus::kError
    }
}

/// HiGHS: IpxWrapper.cpp::reportIpxIpmCrossoverStatus (lines 849-914)
pub fn report_ipx_ipm_crossover_status(
    options: &HighsOptions,
    status: Int,
    ipm_status: bool,
) -> HighsStatus {
    let method_name = if ipm_status { "IPM      " } else { "Crossover" };
    if status == IPX_STATUS_not_run {
        if ipm_status || options.run_crossover == kHighsOnString {
            // Warn if method not run is IPM or method not run is crossover
            // and run_crossover option is "on"
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kWarning,
                "Ipx: {} not run\n",
                method_name
            );
            HighsStatus::kWarning
        } else {
            // OK if method not run is crossover and run_crossover option is
            // not "on"
            HighsStatus::kOk
        }
    } else if status == IPX_STATUS_optimal {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "Ipx: {} optimal\n",
            method_name
        );
        HighsStatus::kOk
    } else if status == IPX_STATUS_imprecise {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Ipx: {} imprecise\n",
            method_name
        );
        HighsStatus::kWarning
    } else if status == IPX_STATUS_primal_infeas {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Ipx: {} primal infeasible\n",
            method_name
        );
        HighsStatus::kWarning
    } else if status == IPX_STATUS_dual_infeas {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Ipx: {} dual infeasible\n",
            method_name
        );
        HighsStatus::kWarning
    } else if status == IPX_STATUS_user_interrupt {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Ipx: {} user interrupt\n",
            method_name
        );
        HighsStatus::kOk
    } else if status == IPX_STATUS_time_limit {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Ipx: {} reached time limit\n",
            method_name
        );
        HighsStatus::kWarning
    } else if status == IPX_STATUS_iter_limit {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Ipx: {} reached iteration limit\n",
            method_name
        );
        HighsStatus::kWarning
    } else if status == IPX_STATUS_no_progress {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Ipx: {} no progress\n",
            method_name
        );
        HighsStatus::kWarning
    } else if status == IPX_STATUS_failed {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Ipx: {} failed\n",
            method_name
        );
        HighsStatus::kError
    } else if status == IPX_STATUS_debug {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Ipx: {} debug\n",
            method_name
        );
        HighsStatus::kError
    } else {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Ipx: {} unrecognised status\n",
            method_name
        );
        HighsStatus::kError
    }
}

/// HiGHS: IpxWrapper.cpp::ipxStatusError (lines 916-931)
///
/// C++ default arguments: `value = -1` (IpxWrapper.h line 76). `value` is
/// C++ `int`; `message` is printed with `%s`, not used as a format.
pub fn ipx_status_error(
    status_error: bool,
    options: &HighsOptions,
    solver: &str,
    message: &str,
    value: i32,
) -> bool {
    if status_error {
        if value < 0 {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "{}: {}\n",
                solver,
                message
            );
        } else {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "{}: {} {}\n",
                solver,
                message,
                value
            );
        }
        // C++ fflush(NULL): flush the output stream.
        let _ = std::io::stdout().flush();
    }
    debug_assert!(!status_error);
    status_error
}

/// HiGHS: IpxWrapper.cpp::illegalIpxSolvedStatus (lines 933-1019)
///
/// PORT NOTE: each check is `found || ipxStatusError(..)` as in C++, so the
/// short-circuit decides which messages are logged.
pub fn illegal_ipx_solved_status(ipx_info: &ipx::Info, options: &HighsOptions) -> bool {
    let mut found_illegal_status = false;
    //========
    // For IPX
    //========
    // Can solve and be optimal
    // Can solve and be imprecise
    // Can solve and be primal_infeas
    // Can solve and be dual_infeas
    // Cannot solve and reach time limit
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_time_limit,
            options,
            "Ipx",
            "solved  status_ipm should not be IPX_STATUS_time_limit",
            -1,
        );
    // Cannot solve and reach iteration limit
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_iter_limit,
            options,
            "Ipx",
            "solved  status_ipm should not be IPX_STATUS_iter_limit",
            -1,
        );
    // Cannot solve and make no progress
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_no_progress,
            options,
            "Ipx",
            "solved  status_ipm should not be IPX_STATUS_no_progress",
            -1,
        );
    // Cannot solve and failed
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_failed,
            options,
            "Ipx",
            "solved  status_ipm should not be IPX_STATUS_failed",
            -1,
        );
    // Cannot solve and debug
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_debug,
            options,
            "Ipx",
            "solved  status_ipm should not be IPX_STATUS_debug",
            -1,
        );
    //==============
    // For crossover
    //==============
    // Can solve and be optimal
    // Can solve and be imprecise
    // Cannot solve with primal infeasibility
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_primal_infeas,
            options,
            "Ipx",
            "solved  status_crossover should not be IPX_STATUS_primal_infeas",
            -1,
        );
    // Cannot solve with dual infeasibility
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_dual_infeas,
            options,
            "Ipx",
            "solved  status_crossover should not be IPX_STATUS_dual_infeas",
            -1,
        );
    // Cannot solve and reach time limit
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_time_limit,
            options,
            "Ipx",
            "solved  status_crossover should not be IPX_STATUS_time_limit",
            -1,
        );
    // Cannot solve and reach time limit
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_iter_limit,
            options,
            "Ipx",
            "solved  status_crossover should not be IPX_STATUS_iter_limit",
            -1,
        );
    // Cannot solve and make no progress
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_no_progress,
            options,
            "Ipx",
            "solved  status_crossover should not be IPX_STATUS_no_progress",
            -1,
        );
    // Cannot solve and failed
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_failed,
            options,
            "Ipx",
            "solved  status_crossover should not be IPX_STATUS_failed",
            -1,
        );
    // Cannot solve and debug
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_debug,
            options,
            "Ipx",
            "solved  status_crossover should not be IPX_STATUS_debug",
            -1,
        );
    found_illegal_status
}

/// HiGHS: IpxWrapper.cpp::illegalIpxStoppedIpmStatus (lines 1021-1061)
pub fn illegal_ipx_stopped_ipm_status(ipx_info: &ipx::Info, options: &HighsOptions) -> bool {
    let mut found_illegal_status = false;
    // Cannot stop and be optimal
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_optimal,
            options,
            "Ipx",
            "stopped status_ipm should not be IPX_STATUS_optimal",
            -1,
        );
    // Cannot stop and be imprecise
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_imprecise,
            options,
            "Ipx",
            "stopped status_ipm should not be IPX_STATUS_imprecise",
            -1,
        );
    // Cannot stop with primal infeasibility
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_primal_infeas,
            options,
            "Ipx",
            "stopped status_ipm should not be IPX_STATUS_primal_infeas",
            -1,
        );
    // Cannot stop with dual infeasibility
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_dual_infeas,
            options,
            "Ipx",
            "stopped status_ipm should not be IPX_STATUS_dual_infeas",
            -1,
        );
    // Can stop with time limit
    // Can stop with iter limit
    // Can stop with no progress
    // Cannot stop and failed - should be error return earlier
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_failed,
            options,
            "Ipx",
            "stopped status_ipm should not be IPX_STATUS_failed",
            -1,
        );
    // Cannot stop and debug - should be error return earlier
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_ipm == IPX_STATUS_debug,
            options,
            "Ipx",
            "stopped status_ipm should not be IPX_STATUS_debug",
            -1,
        );
    found_illegal_status
}

/// HiGHS: IpxWrapper.cpp::illegalIpxStoppedCrossoverStatus (lines 1063-1116)
pub fn illegal_ipx_stopped_crossover_status(ipx_info: &ipx::Info, options: &HighsOptions) -> bool {
    let mut found_illegal_status = false;
    // Cannot stop and be optimal
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_optimal,
            options,
            "Ipx",
            "stopped status_crossover should not be IPX_STATUS_optimal",
            -1,
        );
    // Cannot stop and be imprecise
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_imprecise,
            options,
            "Ipx",
            "stopped status_crossover should not be IPX_STATUS_imprecise",
            -1,
        );
    // Cannot stop with primal infeasibility
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_primal_infeas,
            options,
            "Ipx",
            "stopped status_crossover should not be IPX_STATUS_primal_infeas",
            -1,
        );
    // Cannot stop with dual infeasibility
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_dual_infeas,
            options,
            "Ipx",
            "stopped status_crossover should not be IPX_STATUS_dual_infeas",
            -1,
        );
    // Cannot stop and reach iteration limit
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_iter_limit,
            options,
            "Ipx",
            "stopped status_crossover should not be IPX_STATUS_iter_limit",
            -1,
        );
    // Can stop and reach time limit
    // Cannot stop with no_progress
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_no_progress,
            options,
            "Ipx",
            "stopped status_crossover should not be IPX_STATUS_no_progress",
            -1,
        );
    // Cannot stop and failed - should be error return earlier
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_failed,
            options,
            "Ipx",
            "stopped status_crossover should not be IPX_STATUS_failed",
            -1,
        );
    // Cannot stop and debug - should be error return earlier
    found_illegal_status = found_illegal_status
        || ipx_status_error(
            ipx_info.status_crossover == IPX_STATUS_debug,
            options,
            "Ipx",
            "stopped status_crossover should not be IPX_STATUS_debug",
            -1,
        );
    found_illegal_status
}

/// HiGHS: IpxWrapper.cpp::reportIpmNoProgress (lines 1118-1129)
pub fn report_ipm_no_progress(options: &HighsOptions, ipx_info: &ipx::Info) {
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kWarning,
        "No progress: primal objective value       = {}\n",
        c_fmt_g_width(ipx_info.pobjval, 4, 11)
    );
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kWarning,
        "No progress: max absolute primal residual = {}\n",
        c_fmt_g_width(ipx_info.abs_presidual, 4, 11)
    );
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kWarning,
        "No progress: max absolute   dual residual = {}\n",
        c_fmt_g_width(ipx_info.abs_dresidual, 4, 11)
    );
}

/// HiGHS: IpxWrapper.cpp::getHighsNonVertexSolution (lines 1144-1168)
pub fn get_highs_non_vertex_solution(
    options: &HighsOptions,
    lp: &HighsLp,
    num_col: Int,
    num_row: Int,
    rhs: &[f64],
    constraint_type: &[u8],
    lps: &LpSolver,
    model_status: HighsModelStatus,
    highs_solution: &mut HighsSolution,
) {
    // Get the interior solution (available if IPM was started).
    // GetInteriorSolution() returns the final IPM iterate, regardless if the
    // IPM terminated successfully or not. (Only in case of out-of-memory no
    // solution exists.)
    let mut x: Vec<f64> = vec![0.0; num_col as usize];
    let mut xl: Vec<f64> = vec![0.0; num_col as usize];
    let mut xu: Vec<f64> = vec![0.0; num_col as usize];
    let mut zl: Vec<f64> = vec![0.0; num_col as usize];
    let mut zu: Vec<f64> = vec![0.0; num_col as usize];
    let mut slack: Vec<f64> = vec![0.0; num_row as usize];
    let mut y: Vec<f64> = vec![0.0; num_row as usize];

    lps.get_interior_solution(
        vec_data_mut(&mut x),
        vec_data_mut(&mut xl),
        vec_data_mut(&mut xu),
        vec_data_mut(&mut slack),
        vec_data_mut(&mut y),
        vec_data_mut(&mut zl),
        vec_data_mut(&mut zu),
    );

    // The C++ return value (HighsStatus) is not used here.
    let _ = crate::lp_data::highs_solution_ipx::ipx_solution_to_highs_solution(
        options,
        lp,
        rhs,
        constraint_type,
        num_col,
        num_row,
        &x,
        &slack,
        &y,
        &zl,
        &zu,
        highs_solution,
    );
}

/// HiGHS: IpxWrapper.cpp::reportSolveData (lines 1190-1391)
pub fn report_solve_data(log_options: &HighsLogOptions, ipx_info: &ipx::Info) {
    crate::highs_log_dev!(log_options, HighsLogType::kInfo, "\nIPX Solve data\n");
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    IPX       status = {:>4}\n",
        ipx_info.status
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    IPM       status = {:>4}\n",
        ipx_info.status_ipm
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Crossover status = {:>4}\n",
        ipx_info.status_crossover
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    IPX errflag      = {:>4}\n\n",
        ipx_info.errflag
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    LP variables   = {:>8}\n",
        ipx_info.num_var
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    LP constraints = {:>8}\n",
        ipx_info.num_constr
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    LP entries     = {:>8}\n\n",
        ipx_info.num_entries
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Solver columns = {:>8}\n",
        ipx_info.num_cols_solver
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Solver rows    = {:>8}\n",
        ipx_info.num_rows_solver
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Solver entries = {:>8}\n\n",
        ipx_info.num_entries_solver
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Dualized = {}\n",
        ipx_info.dualized
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Number of dense columns detected = {}\n\n",
        ipx_info.dense_cols
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Dependent rows    = {}\n",
        ipx_info.dependent_rows
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Dependent cols    = {}\n",
        ipx_info.dependent_cols
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Inconsistent rows = {}\n",
        ipx_info.rows_inconsistent
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Inconsistent cols = {}\n",
        ipx_info.cols_inconsistent
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Primal dropped    = {}\n",
        ipx_info.primal_dropped
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Dual   dropped    = {}\n\n",
        ipx_info.dual_dropped
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    |Absolute primal residual| = {}\n",
        c_fmt_g_width(ipx_info.abs_presidual, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    |Absolute   dual residual| = {}\n",
        c_fmt_g_width(ipx_info.abs_dresidual, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    |Relative primal residual| = {}\n",
        c_fmt_g_width(ipx_info.rel_presidual, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    |Relative   dual residual| = {}\n\n",
        c_fmt_g_width(ipx_info.rel_dresidual, 4, 11)
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Primal objective value     = {}\n",
        c_fmt_g_width(ipx_info.pobjval, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Dual   objective value     = {}\n",
        c_fmt_g_width(ipx_info.dobjval, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Relative objective gap     = {}\n",
        c_fmt_g_width(ipx_info.rel_objgap, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Complementarity            = {}\n\n",
        c_fmt_g_width(ipx_info.complementarity, 4, 11)
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    |x| = {}\n",
        c_fmt_g_width(ipx_info.normx, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    |y| = {}\n",
        c_fmt_g_width(ipx_info.normy, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    |z| = {}\n\n",
        c_fmt_g_width(ipx_info.normz, 4, 11)
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Objective value       = {}\n",
        c_fmt_g_width(ipx_info.objval, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Primal infeasibility = {}\n",
        c_fmt_g_width(ipx_info.primal_infeas, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Dual infeasibility   = {}\n\n",
        c_fmt_g_width(ipx_info.dual_infeas, 4, 11)
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    IPM iter   = {}\n",
        ipx_info.iter
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    KKT iter 1 = {}\n",
        ipx_info.kktiter1
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    KKT iter 2 = {}\n",
        ipx_info.kktiter2
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Basis repairs = {}\n",
        ipx_info.basis_repairs
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Updates start     = {}\n",
        ipx_info.updates_start
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Updates ipm       = {}\n",
        ipx_info.updates_ipm
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Updates crossover = {}\n\n",
        ipx_info.updates_crossover
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time total          = {}\n\n",
        c_fmt_f_width(ipx_info.time_total, 2, 8)
    );
    let mut sum_time: f64 = 0.0;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time IPM 1          = {}\n",
        c_fmt_f_width(ipx_info.time_ipm1, 2, 8)
    );
    sum_time += ipx_info.time_ipm1;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time IPM 2          = {}\n",
        c_fmt_f_width(ipx_info.time_ipm2, 2, 8)
    );
    sum_time += ipx_info.time_ipm2;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time starting basis = {}\n",
        c_fmt_f_width(ipx_info.time_starting_basis, 2, 8)
    );
    sum_time += ipx_info.time_starting_basis;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time crossover      = {}\n",
        c_fmt_f_width(ipx_info.time_crossover, 2, 8)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Sum                 = {}\n\n",
        c_fmt_f_width(sum_time, 2, 8)
    );

    sum_time = 0.0;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time kkt_factorize  = {}\n",
        c_fmt_f_width(ipx_info.time_kkt_factorize, 2, 8)
    );
    sum_time += ipx_info.time_kkt_factorize;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time kkt_solve      = {}\n",
        c_fmt_f_width(ipx_info.time_kkt_solve, 2, 8)
    );
    sum_time += ipx_info.time_kkt_solve;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Sum                 = {}\n\n",
        c_fmt_f_width(sum_time, 2, 8)
    );

    sum_time = 0.0;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time maxvol         = {}\n",
        c_fmt_f_width(ipx_info.time_maxvol, 2, 8)
    );
    sum_time += ipx_info.time_maxvol;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time cr1            = {}\n",
        c_fmt_f_width(ipx_info.time_cr1, 2, 8)
    );
    sum_time += ipx_info.time_cr1;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time cr2            = {}\n",
        c_fmt_f_width(ipx_info.time_cr2, 2, 8)
    );
    sum_time += ipx_info.time_cr2;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Sum                 = {}\n\n",
        c_fmt_f_width(sum_time, 2, 8)
    );

    sum_time = 0.0;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time cr1_AAt        = {}\n",
        c_fmt_f_width(ipx_info.time_cr1_AAt, 2, 8)
    );
    sum_time += ipx_info.time_cr1_AAt;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time cr1_pre        = {}\n",
        c_fmt_f_width(ipx_info.time_cr1_pre, 2, 8)
    );
    sum_time += ipx_info.time_cr1_pre;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Sum  cr1            = {}\n\n",
        c_fmt_f_width(sum_time, 2, 8)
    );

    sum_time = 0.0;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time cr2_NNt        = {}\n",
        c_fmt_f_width(ipx_info.time_cr2_NNt, 2, 8)
    );
    sum_time += ipx_info.time_cr2_NNt;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time cr2_B          = {}\n",
        c_fmt_f_width(ipx_info.time_cr2_B, 2, 8)
    );
    sum_time += ipx_info.time_cr2_B;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time cr2_Bt         = {}\n",
        c_fmt_f_width(ipx_info.time_cr2_Bt, 2, 8)
    );
    sum_time += ipx_info.time_cr2_Bt;
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Sum  cr2            = {}\n\n",
        c_fmt_f_width(sum_time, 2, 8)
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Proportion of sparse FTRAN = {}\n",
        c_fmt_g_width(ipx_info.ftran_sparse, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Proportion of sparse BTRAN = {}\n\n",
        c_fmt_g_width(ipx_info.btran_sparse, 4, 11)
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time FTRAN       = {}\n",
        c_fmt_f_width(ipx_info.time_ftran, 2, 8)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time BTRAN       = {}\n",
        c_fmt_f_width(ipx_info.time_btran, 2, 8)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time LU INVERT   = {}\n",
        c_fmt_f_width(ipx_info.time_lu_invert, 2, 8)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time LU UPDATE   = {}\n",
        c_fmt_f_width(ipx_info.time_lu_update, 2, 8)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Mean fill-in     = {}\n",
        c_fmt_g_width(ipx_info.mean_fill, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Max fill-in      = {}\n",
        c_fmt_g_width(ipx_info.max_fill, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Time symb INVERT = {}\n\n",
        c_fmt_g_width(ipx_info.time_symb_invert, 4, 11)
    );

    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Maxvol updates       = {}\n",
        ipx_info.maxvol_updates
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Maxvol skipped       = {}\n",
        ipx_info.maxvol_skipped
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Maxvol passes        = {}\n",
        ipx_info.maxvol_passes
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Tableau num nonzeros = {}\n",
        ipx_info.tbl_nnz
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Tbl max?             = {}\n",
        c_fmt_g_width(ipx_info.tbl_max, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Frobnorm squared     = {}\n",
        c_fmt_g_width(ipx_info.frobnorm_squared, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Lambda max           = {}\n",
        c_fmt_g_width(ipx_info.lambdamax, 4, 11)
    );
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "    Volume increase      = {}\n\n",
        c_fmt_g_width(ipx_info.volume_increase, 4, 11)
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tiny LP with every row kind: row 0 general-bounded (slack column),
    /// row 1 free (dropped), row 2 equality. Hand-computed from fillInIpxData.
    fn tiny_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 3,
            col_cost_: vec![3.0, -4.0],
            col_lower_: vec![0.0, -kHighsInf],
            col_upper_: vec![10.0, kHighsInf],
            row_lower_: vec![1.0, -kHighsInf, 2.0],
            row_upper_: vec![4.0, kHighsInf, 2.0],
            sense_: ObjSense::kMaximize,
            offset_: 0.5,
            ..Default::default()
        };
        // Column-wise: col 0 = rows {0, 1} with values {2, 3},
        //              col 1 = rows {1, 2} with values {5, 7}.
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 1, 2];
        lp.a_matrix_.value_ = vec![2.0, 3.0, 5.0, 7.0];
        lp
    }

    #[test]
    fn fill_in_ipx_data_tiny_lp() {
        let lp = tiny_lp();
        let (mut num_col, mut num_row, mut offset) = (0, 0, 0.0);
        let (mut obj, mut col_lb, mut col_ub) = (Vec::new(), Vec::new(), Vec::new());
        let (mut ap, mut ai, mut ax) = (Vec::new(), Vec::new(), Vec::new());
        let (mut rhs, mut ctype) = (Vec::new(), Vec::new());
        fill_in_ipx_data(
            &lp,
            &mut num_col,
            &mut num_row,
            &mut offset,
            &mut obj,
            &mut col_lb,
            &mut col_ub,
            &mut ap,
            &mut ai,
            &mut ax,
            &mut rhs,
            &mut ctype,
            false,
        );
        assert_eq!(num_col, 3);
        assert_eq!(num_row, 2);
        assert_eq!(offset, -0.5);
        assert_eq!(obj, vec![-3.0, 4.0, 0.0]);
        assert_eq!(col_lb, vec![0.0, f64::NEG_INFINITY, 1.0]);
        assert_eq!(col_ub, vec![10.0, f64::INFINITY, 4.0]);
        assert_eq!(ap, vec![0, 1, 2, 3]);
        assert_eq!(ai, vec![0, 1, 0]);
        assert_eq!(ax, vec![2.0, 7.0, -1.0]);
        assert_eq!(rhs, vec![0.0, 2.0]);
        assert_eq!(ctype, vec![b'=', b'=']);
    }

    /// LP with a free row (row 0) before a boxed row (row 1), so reduced_rowmap is not
    /// the identity. Max 3x0 + 4x1 s.t. 1 <= x0 + 2x1 <= 4, x0 - x1 <= 2, x >= 0 and
    /// x0 + x1 free. Optimum 32/3 at (8/3, 2/3).
    fn free_before_boxed_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 3,
            col_cost_: vec![3.0, 4.0],
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![10.0, 10.0],
            row_lower_: vec![-kHighsInf, 1.0, -kHighsInf],
            row_upper_: vec![kHighsInf, 4.0, 2.0],
            sense_: ObjSense::kMaximize,
            offset_: 0.0,
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 3, 5];
        lp.a_matrix_.index_ = vec![0, 1, 2, 1, 2];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 2.0, -1.0];
        lp
    }

    fn fill_free_before_boxed(h16_ipx_free_rows: bool) -> (Vec<Int>, Vec<Int>, Vec<f64>, Int, Int) {
        let lp = free_before_boxed_lp();
        let (mut num_col, mut num_row, mut offset) = (0, 0, 0.0);
        let (mut obj, mut col_lb, mut col_ub) = (Vec::new(), Vec::new(), Vec::new());
        let (mut ap, mut ai, mut ax) = (Vec::new(), Vec::new(), Vec::new());
        let (mut rhs, mut ctype) = (Vec::new(), Vec::new());
        fill_in_ipx_data(
            &lp,
            &mut num_col,
            &mut num_row,
            &mut offset,
            &mut obj,
            &mut col_lb,
            &mut col_ub,
            &mut ap,
            &mut ai,
            &mut ax,
            &mut rhs,
            &mut ctype,
            h16_ipx_free_rows,
        );
        (ai, ap, ax, num_col, num_row)
    }

    /// `h16_ipx_free_rows` on: the slack of boxed row 1 points at reduced row 0.
    #[test]
    fn h16_ipx_free_rows_fill_on_reduced_slack() {
        let (ai, ap, ax, num_col, num_row) = fill_free_before_boxed(true);
        assert_eq!(num_col, 3);
        assert_eq!(num_row, 2);
        assert_eq!(ap, vec![0, 2, 4, 5]);
        assert_eq!(ai, vec![0, 1, 0, 1, 0]);
        assert_eq!(ax, vec![1.0, 1.0, 2.0, -1.0, -1.0]);
    }

    /// `h16_ipx_free_rows` off: 1.15.1 keeps the unreduced slack index (row 1), which
    /// points at reduced row 1 (the upper-only row). Recorded as the 1.15.1 result.
    #[test]
    fn h16_ipx_free_rows_fill_off_unreduced_slack() {
        let (ai, ap, _ax, num_col, num_row) = fill_free_before_boxed(false);
        assert_eq!(num_col, 3);
        assert_eq!(num_row, 2);
        assert_eq!(ap, vec![0, 2, 4, 5]);
        assert_eq!(ai, vec![0, 1, 0, 1, 1]);
    }

    /// `h16_ipx_free_rows` on, end to end: IPX (solver=ipx, crossover on) reaches the
    /// optimum 32/3 on the free-before-boxed LP.
    #[test]
    fn h16_ipx_free_rows_solve_on_optimum() {
        let mut highs = crate::lp_data::highs_run::Highs::default();
        assert_eq!(
            highs.set_option_value_bool("log_to_console", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("presolve", "off"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", "ipx"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_ipx_free_rows", true),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(free_before_boxed_lp()), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        let objective = highs.get_info().objective_function_value;
        assert!((objective - 32.0 / 3.0).abs() < 1e-6, "{objective}");
        let col = &highs.get_solution().col_value;
        assert!((col[0] - 8.0 / 3.0).abs() < 1e-6, "{col:?}");
        assert!((col[1] - 2.0 / 3.0).abs() < 1e-6, "{col:?}");
    }

    #[test]
    fn fill_in_rhs_and_constraints_kinds() {
        let lp = HighsLp {
            num_row_: 5,
            // > , < , = , general bounded, free
            row_lower_: vec![1.0, -kHighsInf, 2.0, 1.0, -kHighsInf],
            row_upper_: vec![kHighsInf, 5.0, 2.0, 4.0, kHighsInf],
            ..Default::default()
        };
        let (mut rhs, mut ctype) = (Vec::new(), Vec::new());
        fill_in_rhs_and_constraints(&lp, &mut rhs, &mut ctype);
        // The free row is not pushed: the vectors index the reduced rows.
        assert_eq!(rhs, vec![1.0, 5.0, 2.0, 0.0]);
        assert_eq!(ctype, vec![b'>', b'<', b'=', b'=']);
    }
}
