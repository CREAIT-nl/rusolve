//! HiGHS source: ipm/IpxWrapper.cpp (lines 1131-1142, 1170-1188, 1393-1511)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiPO status reporting and the interior solution copy. The messages are
//! `highs_log_user!` calls with the C++ text verbatim.

use crate::ipm::hipo::Int;
use crate::ipm::hipo::ipm::info::Info;
use crate::ipm::hipo::ipm::solver::Solver;
use crate::ipm::hipo::ipm::status::*;
use crate::ipm::ipx;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g_width;

/// HiGHS: IpxWrapper.cpp::reportHipoNoProgress (lines 1131-1142)
pub fn report_hipo_no_progress(options: &HighsOptions, hipo_info: &Info) {
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kWarning,
        "No progress: primal objective value       = {}\n",
        c_fmt_g_width(hipo_info.p_obj, 4, 11)
    );
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kWarning,
        "No progress: max absolute primal residual = {}\n",
        c_fmt_g_width(hipo_info.p_res_abs, 4, 11)
    );
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kWarning,
        "No progress: max absolute   dual residual = {}\n",
        c_fmt_g_width(hipo_info.d_res_abs, 4, 11)
    );
}

/// HiGHS: IpxWrapper.cpp::getHipoNonVertexSolution (lines 1170-1188)
/// `constraint_type` is the `std::vector<char>` of `fill_in_rhs_and_constraints` (bytes).
#[allow(clippy::too_many_arguments)]
pub fn get_hipo_non_vertex_solution(
    options: &HighsOptions,
    lp: &HighsLp,
    num_col: Int,
    num_row: Int,
    rhs: &[f64],
    constraint_type: &[u8],
    hipo: &Solver,
    model_status: HighsModelStatus,
    highs_solution: &mut HighsSolution,
) {
    // Called from `Solver::checkTerminationKkt`.
    let _ = model_status;
    let mut x = vec![0.0; num_col as usize];
    let mut xl = vec![0.0; num_col as usize];
    let mut xu = vec![0.0; num_col as usize];
    let mut zl = vec![0.0; num_col as usize];
    let mut zu = vec![0.0; num_col as usize];
    let mut slack = vec![0.0; num_row as usize];
    let mut y = vec![0.0; num_row as usize];

    hipo.getInteriorSolution(
        &mut x, &mut xl, &mut xu, &mut slack, &mut y, &mut zl, &mut zu,
    );
    // The C++ discards the returned status.
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

/// HiGHS: IpxWrapper.cpp::reportHipoStatus (lines 1393-1447)
pub fn report_hipo_status(options: &HighsOptions, status: Int, hipo: &Solver) -> HighsStatus {
    if hipo.solved() {
        crate::highs_log_user!(&options.log_options, HighsLogType::kInfo, "Hipo: Solved\n");
        return HighsStatus::kOk;
    }
    // these are warnings
    else if status == kStatusTimeLimit {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Time limit\n"
        );
        return HighsStatus::kWarning;
    } else if status == kStatusUserInterrupt {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: User interrupt\n"
        );
        return HighsStatus::kWarning;
    } else if status == kStatusMaxIter {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Reached maximum iterations\n"
        );
        return HighsStatus::kWarning;
    } else if status == kStatusNoProgress {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: No progress\n"
        );
        return HighsStatus::kWarning;
    } else if status == kStatusImprecise {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Imprecise solution\n"
        );
        return HighsStatus::kWarning;
    }
    // these are errors
    else if status == kStatusError {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Internal error\n"
        );
    } else if status == kStatusOverflow {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Integer overflow\n"
        );
    } else if status == kStatusErrorAnalyse {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Error in analyse phase\n"
        );
    } else if status == kStatusErrorFactorise {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Error in factorise phase\n"
        );
    } else if status == kStatusErrorSolve {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Error in solve phase\n"
        );
    } else if status == kStatusBadModel {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Invalid model\n"
        );
    } else {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Unrecognized status\n"
        );
    }
    HighsStatus::kError
}

/// HiGHS: IpxWrapper.cpp::reportHipoCrossoverStatus (lines 1449-1511)
pub fn report_hipo_crossover_status(
    options: &HighsOptions,
    status: ipx::Int,
    is_qp: bool,
) -> HighsStatus {
    if is_qp && options.run_crossover == kHighsOnString {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Crossover not available for QP\n"
        );
        return HighsStatus::kWarning;
    }

    if status == ipx::IPX_STATUS_not_run {
        if options.run_crossover == kHighsOnString {
            // Warn if crossover not run and run_crossover option is "on"
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kWarning,
                "Hipo: Crossover not run\n"
            );
            return HighsStatus::kWarning;
        }
        return HighsStatus::kOk;
    } else if status == ipx::IPX_STATUS_optimal {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "Hipo: Crossover optimal\n"
        );
        return HighsStatus::kOk;
    } else if status == ipx::IPX_STATUS_imprecise {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Crossover imprecise\n"
        );
        return HighsStatus::kWarning;
    } else if status == ipx::IPX_STATUS_primal_infeas {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Crossover primal infeasible\n"
        );
        return HighsStatus::kWarning;
    } else if status == ipx::IPX_STATUS_dual_infeas {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Crossover dual infeasible\n"
        );
        return HighsStatus::kWarning;
    } else if status == ipx::IPX_STATUS_user_interrupt {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Crossover user interrupt\n"
        );
        return HighsStatus::kOk;
    } else if status == ipx::IPX_STATUS_time_limit {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Crossover reached time limit\n"
        );
        return HighsStatus::kWarning;
    } else if status == ipx::IPX_STATUS_iter_limit {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Crossover reached iteration limit\n"
        );
        return HighsStatus::kWarning;
    } else if status == ipx::IPX_STATUS_no_progress {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Hipo: Crossover no progress\n"
        );
        return HighsStatus::kWarning;
    } else if status == ipx::IPX_STATUS_failed {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Crossover failed\n"
        );
        return HighsStatus::kError;
    } else if status == ipx::IPX_STATUS_debug {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Crossover debug\n"
        );
        return HighsStatus::kError;
    } else {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "Hipo: Crossover unrecognised status\n"
        );
        return HighsStatus::kError;
    }
    HighsStatus::kError
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn crossover_not_run_is_ok_unless_forced() {
        let mut options = HighsOptions {
            output_flag: false,
            // The default run_crossover is "on" (C++ HighsOptions.h:773 record
            // default; highs_options.rs:214 matches), under which not_run warns.
            // The first assertion needs a non-"on" value, so set one explicitly.
            run_crossover: String::new(),
            ..HighsOptions::default()
        };
        assert_eq!(
            report_hipo_crossover_status(&options, ipx::IPX_STATUS_not_run, false),
            HighsStatus::kOk
        );
        options.run_crossover = kHighsOnString.to_string();
        assert_eq!(
            report_hipo_crossover_status(&options, ipx::IPX_STATUS_not_run, false),
            HighsStatus::kWarning
        );
        assert_eq!(
            report_hipo_crossover_status(&options, ipx::IPX_STATUS_optimal, true),
            HighsStatus::kWarning
        );
    }
}
