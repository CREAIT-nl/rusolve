//! HiGHS source: pdlp/cupdlp/cupdlp_solver.c (lines 598-897)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Termination and infeasibility tests, the printers and the solve clock.
//!
//! `CUPDLP_CHECK_TIMEOUT(pdhg)` (cupdlp_solver.h) expands to
//! `PDHG_Compute_SolvingTime(pdhg); if (pdhg->timers->dSolvingTime >
//! pdhg->settings->dTimeLim) { retcode = RETCODE_FAILED; goto exit_cleanup; }`
//! and is written out at each use; HiGHS uses it
//! nowhere in the CPU build.

use super::super::defs::*;
use super::super::utils::*;

/// HiGHS: cupdlp_solver.c::PDHG_Check_Data (lines 598-708)
///
/// Never called in HiGHS (no call site in `pdlp/`); ported for completeness.
/// The counting loops and the `CUPDLP_ASSERT`s (`assert`, active in debug
/// builds only, so `debug_assert!`) stay; the printf lines go. The row loop
/// `for (iSeq = lp->nCols; iSeq < lp->nCols; ...)` runs zero times upstream
/// and is kept as written.
pub fn PDHG_Check_Data(work: &mut CUPDLPwork) {
    let problem = &work.problem;
    let lp = &problem.data;
    // stepsize unused in this function (lines 601); iterates used below.
    let iterates = &work.iterates;

    let mut nFreeCol: cupdlp_int = 0;
    let mut nFixedCol: cupdlp_int = 0;
    let mut nUpperCol: cupdlp_int = 0;
    let mut nLowerCol: cupdlp_int = 0;
    let mut nRangedCol: cupdlp_int = 0;
    let mut nFreeRow: cupdlp_int = 0;
    let mut nFixedRow: cupdlp_int = 0;
    let mut nUpperRow: cupdlp_int = 0;
    let mut nLowerRow: cupdlp_int = 0;
    let mut nRangedRow: cupdlp_int = 0;

    let iter = work.timers.nIter;
    let x = &iterates.x[(iter % 2) as usize];
    let y = &iterates.y[(iter % 2) as usize];

    for iSeq in 0..lp.nCols {
        let hasLower = problem.lower[iSeq as usize] > -f64::INFINITY;
        let hasUpper = problem.upper[iSeq as usize] < f64::INFINITY;

        if !hasLower && !hasUpper {
            nFreeCol += 1;
            // printf dropped (lines 625-626)
        }

        if hasLower && hasUpper {
            if problem.lower[iSeq as usize] == problem.upper[iSeq as usize] {
                nFixedCol += 1;
                // printf commented out (line 632)
            } else {
                nRangedCol += 1;
            }
        }

        if hasLower {
            // CUPDLP_ASSERT commented out (line 639)
            nLowerCol += (!hasUpper) as cupdlp_int;
        }

        if hasUpper {
            // CUPDLP_ASSERT commented out (line 645)
            nUpperCol += (!hasLower) as cupdlp_int;
        }
    }

    // Loop runs zero times upstream (iSeq starts at lp->nCols, condition is
    // iSeq < lp->nCols). Kept as written.
    for iSeq in lp.nCols..lp.nCols {
        let hasLower = problem.lower[iSeq as usize] > -f64::INFINITY;
        let hasUpper = problem.upper[iSeq as usize] < f64::INFINITY;

        if !hasLower && !hasUpper {
            nFreeRow += 1;
            // printf dropped (lines 656-657)
        }

        if hasLower && hasUpper {
            if problem.lower[iSeq as usize] == problem.upper[iSeq as usize] {
                nFixedRow += 1;
            } else {
                nRangedRow += 1;
            }
        }

        if hasLower {
            // CUPDLP_ASSERT commented out (line 668)
            nLowerRow += (!hasUpper) as cupdlp_int;
        }

        if hasUpper {
            // CUPDLP_ASSERT commented out (line 673)
            nUpperRow += (!hasLower) as cupdlp_int;
        }
    }

    for iRow in 0..lp.nRows {
        debug_assert!(y.data[iRow as usize] < f64::INFINITY);
        debug_assert!(y.data[iRow as usize] > -f64::INFINITY);
    }

    for iRow in 0..lp.nRows {
        let csr = lp.csr_matrix.as_ref().unwrap();
        if csr.rowMatBeg[(iRow + 1) as usize] - csr.rowMatBeg[iRow as usize] == 1 {
            // printf dropped (lines 687-688)
        }
    }

    debug_assert!(nRangedRow == 0);
    // printf block dropped (lines 693-704); the counters are read only there.
    let _ = (
        nFreeCol, nFixedCol, nRangedCol, nLowerCol, nUpperCol, nFreeRow, nFixedRow, nRangedRow,
        nLowerRow, nUpperRow,
    );

    // We need to test problems ranged row-bounds more carefully.
    debug_assert!(nRangedRow == 0);
}

/// HiGHS: cupdlp_solver.c::PDHG_Check_Primal_Infeasibility (lines 710-723)
pub fn PDHG_Check_Primal_Infeasibility(
    pdhg: &CUPDLPwork,
    dPrimalInfeasObj: cupdlp_float,
    dPrimalInfeasRes: cupdlp_float,
) -> termination_code {
    let resobj = &pdhg.resobj;

    let mut primalCode: termination_code = FEASIBLE;

    if dPrimalInfeasObj > 0.0 {
        if dPrimalInfeasRes < resobj.dFeasTol * dPrimalInfeasObj {
            primalCode = INFEASIBLE;
        }
    }

    primalCode
}

/// HiGHS: cupdlp_solver.c::PDHG_Check_Dual_Infeasibility (lines 725-738)
pub fn PDHG_Check_Dual_Infeasibility(
    pdhg: &CUPDLPwork,
    dDualInfeasObj: cupdlp_float,
    dDualInfeasRes: cupdlp_float,
) -> termination_code {
    let resobj = &pdhg.resobj;

    let mut dualCode: termination_code = FEASIBLE;

    if dDualInfeasObj < 0.0 {
        if dDualInfeasRes < -resobj.dFeasTol * dDualInfeasObj {
            dualCode = INFEASIBLE;
        }
    }

    dualCode
}

/// HiGHS: cupdlp_solver.c::PDHG_Check_Infeasibility (lines 740-795)
///
/// Current primal, current dual, average primal, average dual, in that
/// order; a later hit overwrites `termInfeasIterate`. The `bool_print`
/// block is printf only.
pub fn PDHG_Check_Infeasibility(pdhg: &mut CUPDLPwork, bool_print: i32) -> termination_code {
    let mut t_code: termination_code = FEASIBLE;

    // current solution

    // primal infeasibility
    {
        let obj = pdhg.resobj.dPrimalInfeasObj;
        let res = pdhg.resobj.dPrimalInfeasRes;
        if PDHG_Check_Primal_Infeasibility(pdhg, obj, res) == INFEASIBLE {
            pdhg.resobj.primalCode = INFEASIBLE;
            pdhg.resobj.termInfeasIterate = LAST_ITERATE;
            t_code = INFEASIBLE_OR_UNBOUNDED;
        }
    }

    // dual infeasibility
    {
        let obj = pdhg.resobj.dDualInfeasObj;
        let res = pdhg.resobj.dDualInfeasRes;
        if PDHG_Check_Dual_Infeasibility(pdhg, obj, res) == INFEASIBLE {
            pdhg.resobj.dualCode = INFEASIBLE;
            pdhg.resobj.termInfeasIterate = LAST_ITERATE;
            t_code = INFEASIBLE_OR_UNBOUNDED;
        }
    }

    // average solution
    // primal infeasibility
    {
        let obj = pdhg.resobj.dPrimalInfeasObjAverage;
        let res = pdhg.resobj.dPrimalInfeasResAverage;
        if PDHG_Check_Primal_Infeasibility(pdhg, obj, res) == INFEASIBLE {
            pdhg.resobj.primalCode = INFEASIBLE;
            pdhg.resobj.termInfeasIterate = AVERAGE_ITERATE;
            t_code = INFEASIBLE_OR_UNBOUNDED;
        }
    }

    // dual infeasibility
    {
        let obj = pdhg.resobj.dDualInfeasObjAverage;
        let res = pdhg.resobj.dDualInfeasResAverage;
        if PDHG_Check_Dual_Infeasibility(pdhg, obj, res) == INFEASIBLE {
            pdhg.resobj.dualCode = INFEASIBLE;
            pdhg.resobj.termInfeasIterate = AVERAGE_ITERATE;
            t_code = INFEASIBLE_OR_UNBOUNDED;
        }
    }

    // printf block dropped (lines 781-792)

    t_code
}

/// HiGHS: cupdlp_solver.c::PDHG_Check_Termination (lines 797-818)
///
/// The `#if PDHG_DISPLAY_TERMINATION_CHECK` block (defined as 1) is printf
/// only and is dropped.
pub fn PDHG_Check_Termination(pdhg: &CUPDLPwork, bool_print: i32) -> cupdlp_bool {
    // problem unused (line 798)
    let settings = &pdhg.settings;
    let resobj = &pdhg.resobj;
    let scaling = &pdhg.scaling;

    // #if PDHG_DISPLAY_TERMINATION_CHECK block dropped (lines 802-812)

    (resobj.dPrimalFeas < settings.dPrimalTol * (1.0 + scaling.dNormRhs))
        && (resobj.dDualFeas < settings.dDualTol * (1.0 + scaling.dNormCost))
        && (resobj.dRelObjGap < settings.dGapTol)
}

/// HiGHS: cupdlp_solver.c::PDHG_Check_Termination_Average (lines 820-841)
pub fn PDHG_Check_Termination_Average(pdhg: &CUPDLPwork, bool_print: i32) -> cupdlp_bool {
    // problem unused (line 821)
    let settings = &pdhg.settings;
    let resobj = &pdhg.resobj;
    let scaling = &pdhg.scaling;

    // #if PDHG_DISPLAY_TERMINATION_CHECK block dropped (lines 825-834)

    (resobj.dPrimalFeasAverage < settings.dPrimalTol * (1.0 + scaling.dNormRhs))
        && (resobj.dDualFeasAverage < settings.dDualTol * (1.0 + scaling.dNormCost))
        && (resobj.dRelObjGapAverage < settings.dGapTol)
}

/// HiGHS: cupdlp_solver.c::PDHG_Print_Header (lines 843-849)
///
/// Print-only: the faithful translation is an empty body.
pub fn PDHG_Print_Header(pdhg: &CUPDLPwork) {
    let _ = pdhg;
}

/// HiGHS: cupdlp_solver.c::PDHG_Print_Iter (lines 851-870)
///
/// Print-only.
pub fn PDHG_Print_Iter(pdhg: &CUPDLPwork) {
    let _ = pdhg;
}

/// HiGHS: cupdlp_solver.c::PDHG_Print_Iter_Average (lines 872-892)
///
/// Print-only.
pub fn PDHG_Print_Iter_Average(pdhg: &CUPDLPwork) {
    let _ = pdhg;
}

/// HiGHS: cupdlp_solver.c::PDHG_Compute_SolvingTime (lines 894-897)
///
/// `getTimeStamp()` is `my_clock()` = `(double)time(NULL)` (whole seconds,
/// `CUPDLP_TIMER` defined).
pub fn PDHG_Compute_SolvingTime(pdhg: &mut CUPDLPwork) {
    let timers = &mut pdhg.timers;
    timers.dSolvingTime = getTimeStamp() - timers.dSolvingBeg;
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn infeasibility_tests_use_feas_tol() {
        let mut w = CUPDLPwork::default();
        w.resobj.dFeasTol = 1e-8;
        assert_eq!(PDHG_Check_Primal_Infeasibility(&w, 1.0, 1e-9), INFEASIBLE);
        assert_eq!(PDHG_Check_Primal_Infeasibility(&w, 0.0, 0.0), FEASIBLE);
        assert_eq!(PDHG_Check_Dual_Infeasibility(&w, -1.0, 1e-9), INFEASIBLE);
        w.resobj.dDualInfeasObjAverage = -1.0;
        w.resobj.dDualInfeasResAverage = 0.0;
        assert_eq!(PDHG_Check_Infeasibility(&mut w, 0), INFEASIBLE_OR_UNBOUNDED);
        assert_eq!(w.resobj.dualCode, INFEASIBLE);
        assert_eq!(w.resobj.termInfeasIterate, AVERAGE_ITERATE);
    }

    #[test]
    fn termination_is_relative_to_norms() {
        let mut w = CUPDLPwork::default();
        w.settings.dPrimalTol = 1e-4;
        w.settings.dDualTol = 1e-4;
        w.settings.dGapTol = 1e-4;
        w.scaling.dNormRhs = 9.0;
        w.scaling.dNormCost = 0.0;
        w.resobj.dPrimalFeas = 9e-4;
        w.resobj.dDualFeas = 9e-5;
        w.resobj.dRelObjGap = 0.0;
        assert!(PDHG_Check_Termination(&w, 0));
        w.resobj.dDualFeas = 1e-4;
        assert!(!PDHG_Check_Termination(&w, 0));
        w.resobj.dPrimalFeasAverage = 9e-4;
        w.resobj.dDualFeasAverage = 9e-5;
        w.resobj.dRelObjGapAverage = 0.0;
        assert!(PDHG_Check_Termination_Average(&w, 0));
        w.resobj.dDualFeasAverage = 1e-4;
        assert!(!PDHG_Check_Termination_Average(&w, 0));
    }
}
