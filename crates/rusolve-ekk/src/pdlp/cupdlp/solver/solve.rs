//! HiGHS source: pdlp/cupdlp/cupdlp_solver.c (lines 899-1215)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The PDHG main loop.

use super::super::defs::*;
use super::super::proj::*;
use super::super::step::*;
use super::super::utils::*;
use super::*;

/// HiGHS: cupdlp_solver.c::PDHG_Solve (lines 899-1215)
///
/// Loop, for `nIter` in `0..nIterLim` (`timers->nIter` is the loop
/// variable):
/// 1. `PDHG_Compute_SolvingTime`;
/// 2. `bool_checking = nIter < 10 || nIter == nIterLim - 1 ||
///    dSolvingTime > dTimeLim || !(nIter % CUPDLP_RELEASE_INTERVAL)` with
///    `CUPDLP_RELEASE_INTERVAL == 40`; the `bool_print` logic is print-only;
/// 3. when checking: `Compute_Average_Iterate`, `Compute_Residuals`,
///    `Compute_Infeas_Residuals`, then the oracle hook `pdlp_check`
///    ([`crate::pdlp::oracle_trace`]), then in order
///    `Check_Termination` (LAST_ITERATE, OPTIMAL, break),
///    `Check_Termination_Average` (copy the average iterate and
///    `dSlackPos/NegAverage` into the current ones; AVERAGE_ITERATE,
///    OPTIMAL, break), `Check_Infeasibility(pdhg, 0) ==
///    INFEASIBLE_OR_UNBOUNDED` (break), time limit and `nIter >= nIterLim -
///    1` (TIMELIMIT_OR_ITERLIMIT, break), `PDHG_Restart_Iterate`;
/// 4. `PDHG_Update_Iterate == RETCODE_FAILED` -> TIMELIMIT_OR_ITERLIMIT,
///    break.
///
/// Right after the loop (on every exit path, before the post-loop block)
/// the oracle hook `pdlp_end(pdhg)` fires. The post-loop block is
/// logging only. `CUPDLP_CALL(PDHG_Init_Step_Sizes)`
/// is written out.
pub fn PDHG_Solve(has_variables: &cupdlp_int, pdhg: &mut CUPDLPwork) -> cupdlp_retcode {
    let mut retcode: cupdlp_retcode = RETCODE_OK;

    // CUPDLPproblem *problem = pdhg->problem;
    // CUPDLPstepsize *stepsize = pdhg->stepsize;
    // CUPDLPsettings *settings = pdhg->settings;
    // CUPDLPresobj *resobj = pdhg->resobj;
    // CUPDLPiterates *iterates = pdhg->iterates;
    // CUPDLPtimers *timers = pdhg->timers;
    // CUPDLPscaling *scaling = pdhg->scaling;

    pdhg.timers.nIter = 0;
    pdhg.timers.dSolvingBeg = getTimeStamp();

    // #if PDLP_DEBUG_LOG (lines 913-917): dropped (PDLP_DEBUG_LOG == 0).

    // PDHG_Init_Data does nothing!
    PDHG_Init_Data(pdhg);

    'exit_cleanup: {
        // CUPDLP_CALL(PDHG_Init_Step_Sizes(pdhg));
        if PDHG_Init_Step_Sizes(pdhg) != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        PDHG_Init_Variables(has_variables, pdhg);

        // todo: translate check_data into cuda or do it on cpu
        // PDHG_Check_Data(pdhg);

        // PDHG_Print_Header(pdhg);

        // Repeat the iteration logging header periodically if logging style
        // is minimal (pdhg->settings->nLogLevel=1), initialising
        // iter_log_since_header so that an initial header is printed
        const iter_log_between_header: i32 = 50;
        let mut iter_log_since_header: i32 = iter_log_between_header;

        // #if PDLP_DEBUG_LOG (line 936): dropped.

        // Main loop!
        // for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter) {
        pdhg.timers.nIter = 0;
        while pdhg.timers.nIter < pdhg.settings.nIterLim {
            // #if PDLP_DEBUG_LOG (lines 940-947): dropped.

            PDHG_Compute_SolvingTime(pdhg);

            // #if CUPDLP_DUMP_ITERATES_STATS && CUPDLP_DEBUG (lines 950-952): dropped (CUPDLP_DEBUG not defined).

            // int bool_checking = (timers->nIter < 10) ||
            //     (timers->nIter == (settings->nIterLim - 1)) ||
            //     (timers->dSolvingTime > settings->dTimeLim);
            let mut bool_checking: bool = (pdhg.timers.nIter < 10)
                || (pdhg.timers.nIter == (pdhg.settings.nIterLim - 1))
                || (pdhg.timers.dSolvingTime > pdhg.settings.dTimeLim);

            // int bool_print = 0;
            // #if CUPDLP_DEBUG (lines 957-959): dropped (CUPDLP_DEBUG not defined).
            // #else (lines 960-971):
            // bool_checking =
            //     (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
            bool_checking = bool_checking || (pdhg.timers.nIter % CUPDLP_RELEASE_INTERVAL == 0);

            // bool_print =
            //     (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
            //     settings->nLogInterval))) ||
            //     (timers->nIter == (settings->nIterLim - 1)) ||
            //     (timers->dSolvingTime > settings->dTimeLim);
            let bool_print: bool = (bool_checking
                && (pdhg.timers.nIter % (CUPDLP_RELEASE_INTERVAL * pdhg.settings.nLogInterval)
                    == 0))
                || (pdhg.timers.nIter == (pdhg.settings.nIterLim - 1))
                || (pdhg.timers.dSolvingTime > pdhg.settings.dTimeLim);

            // // Ensure that bool_print is false if the logging level has been
            // // set to 0 via the HiGHS option
            // bool_print = pdhg->settings->nLogLevel>0 && bool_print;
            let bool_print: bool = pdhg.settings.nLogLevel > 0 && bool_print;

            // // Full printing is false only if the logging level has been set
            // // to 0 or 1 via the HiGHS option
            // int full_print = pdhg->settings->nLogLevel >= 2;
            let full_print: bool = pdhg.settings.nLogLevel >= 2;

            if bool_checking {
                PDHG_Compute_Average_Iterate(pdhg);
                PDHG_Compute_Residuals(pdhg);
                PDHG_Compute_Infeas_Residuals(pdhg);

                // Oracle hook: pdlp_check (after PDHG_Compute_Infeas_Residuals).
                crate::pdlp::oracle_trace::pdlp_check(pdhg);

                if bool_print {
                    // With reduced printing, the header is only needed for the
                    // first iteration since only average iteration printing is
                    // carried out
                    if full_print || iter_log_since_header == iter_log_between_header {
                        // PDHG_Print_Header(pdhg);
                        iter_log_since_header = 0;
                    }
                    // int print_iter = timers->nIter == 0 || full_print;
                    let _print_iter: bool = pdhg.timers.nIter == 0 || full_print;
                    // int print_iter_average = timers->nIter > 0 || full_print;
                    let _print_iter_average: bool = pdhg.timers.nIter > 0 || full_print;
                    // if (print_iter) PDHG_Print_Iter(pdhg);
                    // if (print_iter_average) PDHG_Print_Iter_Average(pdhg);
                    iter_log_since_header += 1;
                }

                // #if PDLP_DEBUG_LOG (lines 996-1012): dropped.

                // // Termination check printing is only done when printing is full
                // int termination_print = bool_print && full_print;
                let termination_print: i32 = (bool_print && full_print) as i32;

                if PDHG_Check_Termination(pdhg, termination_print) {
                    // resobj->termIterate = LAST_ITERATE;
                    pdhg.resobj.termIterate = LAST_ITERATE;
                    // resobj->termCode = OPTIMAL;
                    pdhg.resobj.termCode = OPTIMAL;
                    break;
                }

                if PDHG_Check_Termination_Average(pdhg, termination_print) {
                    // cupdlp_int iter = pdhg->timers->nIter;
                    let iter: cupdlp_int = pdhg.timers.nIter;

                    // CUPDLP_COPY_VEC(x->data, iterates->xAverage->data, cupdlp_float, problem->nCols);
                    // CUPDLP_COPY_VEC(y->data, iterates->yAverage->data, cupdlp_float, problem->nRows);
                    // CUPDLP_COPY_VEC(ax->data, iterates->axAverage->data, cupdlp_float, problem->nRows);
                    // CUPDLP_COPY_VEC(aty->data, iterates->atyAverage->data, cupdlp_float, problem->nCols);
                    // CUPDLP_COPY_VEC(resobj->dSlackPos, resobj->dSlackPosAverage, cupdlp_float, problem->nCols);
                    // CUPDLP_COPY_VEC(resobj->dSlackNeg, resobj->dSlackNegAverage, cupdlp_float, problem->nCols);
                    {
                        let idx = (iter % 2) as usize;
                        let n_cols = pdhg.problem.nCols;
                        let n_rows = pdhg.problem.nRows;
                        cupdlp_copy(
                            &mut pdhg.iterates.x[idx].data,
                            &pdhg.iterates.xAverage.data,
                            n_cols,
                        );
                        cupdlp_copy(
                            &mut pdhg.iterates.y[idx].data,
                            &pdhg.iterates.yAverage.data,
                            n_rows,
                        );
                        cupdlp_copy(
                            &mut pdhg.iterates.ax[idx].data,
                            &pdhg.iterates.axAverage.data,
                            n_rows,
                        );
                        cupdlp_copy(
                            &mut pdhg.iterates.aty[idx].data,
                            &pdhg.iterates.atyAverage.data,
                            n_cols,
                        );
                        cupdlp_copy(
                            &mut pdhg.resobj.dSlackPos,
                            &pdhg.resobj.dSlackPosAverage,
                            n_cols,
                        );
                        cupdlp_copy(
                            &mut pdhg.resobj.dSlackNeg,
                            &pdhg.resobj.dSlackNegAverage,
                            n_cols,
                        );
                    }

                    // resobj->termIterate = AVERAGE_ITERATE;
                    pdhg.resobj.termIterate = AVERAGE_ITERATE;
                    // resobj->termCode = OPTIMAL;
                    pdhg.resobj.termCode = OPTIMAL;
                    break;
                }

                if PDHG_Check_Infeasibility(pdhg, 0) == INFEASIBLE_OR_UNBOUNDED {
                    // resobj->termCode = INFEASIBLE_OR_UNBOUNDED;
                    pdhg.resobj.termCode = INFEASIBLE_OR_UNBOUNDED;
                    break;
                }

                if pdhg.timers.dSolvingTime > pdhg.settings.dTimeLim {
                    // resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
                    pdhg.resobj.termCode = TIMELIMIT_OR_ITERLIMIT;
                    break;
                }

                if pdhg.timers.nIter >= (pdhg.settings.nIterLim - 1) {
                    // resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
                    pdhg.resobj.termCode = TIMELIMIT_OR_ITERLIMIT;
                    break;
                }

                PDHG_Restart_Iterate(pdhg);
            }

            // // CUPDLP_CALL(PDHG_Update_Iterate(pdhg));
            // if (PDHG_Update_Iterate(pdhg) == RETCODE_FAILED) {
            //     resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
            //     break;
            // }
            if PDHG_Update_Iterate(pdhg) == RETCODE_FAILED {
                pdhg.resobj.termCode = TIMELIMIT_OR_ITERLIMIT;
                break;
            }

            pdhg.timers.nIter += 1;
        }

        // Oracle hook: pdlp_end (right after the for loop, on every exit path).
        crate::pdlp::oracle_trace::pdlp_end(pdhg);

        // // print at last
        // if (pdhg->settings->nLogLevel>0 && timers->nIter > 0) {
        //     int full_print = pdhg->settings->nLogLevel >= 2;
        //     if (full_print) PDHG_Print_Header(pdhg);
        //     int print_iter = resobj->termIterate == LAST_ITERATE || full_print;
        //     int print_iter_average = resobj->termIterate == AVERAGE_ITERATE || full_print;
        //     if (print_iter) PDHG_Print_Iter(pdhg);
        //     if (print_iter_average) PDHG_Print_Iter_Average(pdhg);
        // }
        // All print calls dropped.

        // if (pdhg->settings->nLogLevel>0) {
        //     ... switch on resobj->termCode, cupdlp_printf calls ...
        // }
        // All print calls dropped.

        // #if PDHG_USE_TIMERS (lines 1173-1197): print-only, dropped.

        // #ifndef CUPDLP_CPU (lines 1199-1210): GPU code, dropped.
    }

    // exit_cleanup:
    // //fclose(pdhg->debug_pdlp_log_file_);
    // return retcode;
    retcode
}

#[cfg(test)]
mod tests {
    //! Unit tests. The loop itself is covered by
    //! the oracle test `tests/t_w3_cupdlp.rs`.
    use super::*;

    /// nIterLim = 1: the first iteration is a check iteration and hits the
    /// iteration limit before any update, on an empty problem.
    #[test]
    fn iteration_limit_one_stops_at_first_check() {
        let mut w = CUPDLPwork::default();
        PDHG_Alloc(&mut w);
        // PDHG_Init_Step_Sizes reads problem->data->csc_matrix->MatElemNormInf
        // (C dereferences it unconditionally on the fixed line-search path);
        // the empty problem still needs that matrix object present.
        w.problem.data.csc_matrix = Some(Box::default());
        w.settings.nIterLim = 1;
        w.settings.dTimeLim = 1e9;
        w.settings.dPrimalTol = -1.0;
        w.settings.dDualTol = -1.0;
        w.settings.dGapTol = -1.0;
        let rc = PDHG_Solve(&0, &mut w);
        assert_eq!(rc, RETCODE_OK);
        assert_eq!(w.resobj.termCode, TIMELIMIT_OR_ITERLIMIT);
        assert_eq!(w.timers.nIter, 0);
    }
}
