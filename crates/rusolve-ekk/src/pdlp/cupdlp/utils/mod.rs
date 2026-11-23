//! HiGHS source: pdlp/cupdlp/cupdlp_utils.h, pdlp/cupdlp/cupdlp_utils.c (lines 651-1030)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `*_SetUserParam` (651-811) and the `*_Alloc` initialisers (813-1030).
//! The rest of `cupdlp_utils.c` is in [`work`] (1032-1158: `vec_Alloc`,
//! `PDHG_Alloc`, `PDHG_Create`, clocks) and [`matrix`] (1160-1525: matrix
//! copies, conversions and allocation).
//!
//! Not ported: the `*_clear` functions 22-417 (`Drop`), the printers
//! 419-558 and 1527-1628 (`PDHG_PrintPDHGParam`, called at the end of
//! `PDHG_SetUserParam`, prints only and is dropped), `getUserParam`
//! 560-649 (command-line parsing), `writeJson` / `writeSol` 1630-1783
//! (file writers, deferred), the `PDLP_DEBUG_LOG` functions 1785-1850.
//!
//! The `cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, ...` arrays
//! are slices of length `N_INT_USER_PARAM` / `N_FLOAT_USER_PARAM`.
//! `CUPDLP_INIT_*` (`malloc` of the struct) is the field already present
//! in [`CUPDLPwork`]; every field the C later reads is assigned by the
//! matching `*_Alloc`.

pub mod matrix;
pub mod work;
pub use matrix::*;
pub use work::*;

use super::defs::*;
use super::linalg::*;

/// HiGHS: cupdlp_utils.c::settings_SetUserParam (lines 651-704)
pub fn settings_SetUserParam(
    settings: &mut CUPDLPsettings,
    ifChangeIntParam: &[cupdlp_bool],
    intParam: &[cupdlp_int],
    ifChangeFloatParam: &[cupdlp_bool],
    floatParam: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    if ifChangeIntParam[N_ITER_LIM] {
        settings.nIterLim = intParam[N_ITER_LIM];
    }

    if ifChangeIntParam[N_LOG_LEVEL] {
        settings.nLogLevel = intParam[N_LOG_LEVEL];
    }

    if ifChangeIntParam[N_LOG_INTERVAL] {
        settings.nLogInterval = intParam[N_LOG_INTERVAL];
    }

    if ifChangeIntParam[IF_SCALING] {
        settings.ifScaling = intParam[IF_SCALING];
    }

    if ifChangeIntParam[I_SCALING_METHOD] {
        settings.iScalingMethod = intParam[I_SCALING_METHOD];
    }

    if ifChangeFloatParam[D_SCALING_LIMIT] {
        settings.dScalingLimit = floatParam[D_SCALING_LIMIT];
    }

    if ifChangeFloatParam[D_PRIMAL_TOL] {
        settings.dPrimalTol = floatParam[D_PRIMAL_TOL];
    }

    if ifChangeFloatParam[D_DUAL_TOL] {
        settings.dDualTol = floatParam[D_DUAL_TOL];
    }

    if ifChangeFloatParam[D_GAP_TOL] {
        settings.dGapTol = floatParam[D_GAP_TOL];
    }

    if ifChangeFloatParam[D_TIME_LIM] {
        settings.dTimeLim = floatParam[D_TIME_LIM];
    }

    if ifChangeIntParam[E_RESTART_METHOD] {
        settings.eRestartMethod = intParam[E_RESTART_METHOD];
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::resobj_SetUserParam (lines 706-719)
pub fn resobj_SetUserParam(
    resobj: &mut CUPDLPresobj,
    ifChangeIntParam: &[cupdlp_bool],
    intParam: &[cupdlp_int],
    ifChangeFloatParam: &[cupdlp_bool],
    floatParam: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    if ifChangeFloatParam[D_FEAS_TOL] {
        resobj.dFeasTol = floatParam[D_FEAS_TOL];
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::iterates_SetUserParam (lines 721-730) (empty upstream)
pub fn iterates_SetUserParam(
    iterates: &mut CUPDLPiterates,
    ifChangeIntParam: &[cupdlp_bool],
    intParam: &[cupdlp_int],
    ifChangeFloatParam: &[cupdlp_bool],
    floatParam: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::stepsize_SetUserParam (lines 732-745)
pub fn stepsize_SetUserParam(
    stepsize: &mut CUPDLPstepsize,
    ifChangeIntParam: &[cupdlp_bool],
    intParam: &[cupdlp_int],
    ifChangeFloatParam: &[cupdlp_bool],
    floatParam: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    if ifChangeIntParam[E_LINE_SEARCH_METHOD] {
        stepsize.eLineSearchMethod = intParam[E_LINE_SEARCH_METHOD];
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::scaling_SetUserParam (lines 747-768)
pub fn scaling_SetUserParam(
    scaling: &mut CUPDLPscaling,
    ifChangeIntParam: &[cupdlp_bool],
    intParam: &[cupdlp_int],
    ifChangeFloatParam: &[cupdlp_bool],
    floatParam: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    if ifChangeIntParam[IF_RUIZ_SCALING] {
        scaling.ifRuizScaling = intParam[IF_RUIZ_SCALING];
    }

    if ifChangeIntParam[IF_L2_SCALING] {
        scaling.ifL2Scaling = intParam[IF_L2_SCALING];
    }

    if ifChangeIntParam[IF_PC_SCALING] {
        scaling.ifPcScaling = intParam[IF_PC_SCALING];
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::timers_SetUserParam (lines 770-779) (empty upstream)
pub fn timers_SetUserParam(
    timers: &mut CUPDLPtimers,
    ifChangeIntParam: &[cupdlp_bool],
    intParam: &[cupdlp_int],
    ifChangeFloatParam: &[cupdlp_bool],
    floatParam: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::PDHG_SetUserParam (lines 781-811)
///
/// Calls the six setters in C order through `CUPDLP_CALL`; the trailing
/// `PDHG_PrintPDHGParam(w)` prints only and is dropped.
pub fn PDHG_SetUserParam(
    w: &mut CUPDLPwork,
    ifChangeIntParam: &[cupdlp_bool],
    intParam: &[cupdlp_int],
    ifChangeFloatParam: &[cupdlp_bool],
    floatParam: &[cupdlp_float],
) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;

    // CUPDLPsettings *settings = w->settings;
    // CUPDLPstepsize *stepsize = w->stepsize;
    // CUPDLPresobj *resobj = w->resobj;
    // CUPDLPiterates *iterates = w->iterates;
    // CUPDLPscaling *scaling = w->scaling;
    // CUPDLPtimers *timers = w->timers;

    'exit_cleanup: {
        // CUPDLP_CALL(settings_SetUserParam(...));
        let rc = settings_SetUserParam(
            &mut w.settings,
            ifChangeIntParam,
            intParam,
            ifChangeFloatParam,
            floatParam,
        );
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(stepsize_SetUserParam(...));
        let rc = stepsize_SetUserParam(
            &mut w.stepsize,
            ifChangeIntParam,
            intParam,
            ifChangeFloatParam,
            floatParam,
        );
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(resobj_SetUserParam(...));
        let rc = resobj_SetUserParam(
            &mut w.resobj,
            ifChangeIntParam,
            intParam,
            ifChangeFloatParam,
            floatParam,
        );
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(iterates_SetUserParam(...));
        let rc = iterates_SetUserParam(
            &mut w.iterates,
            ifChangeIntParam,
            intParam,
            ifChangeFloatParam,
            floatParam,
        );
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(scaling_SetUserParam(...));
        let rc = scaling_SetUserParam(
            &mut w.scaling,
            ifChangeIntParam,
            intParam,
            ifChangeFloatParam,
            floatParam,
        );
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(timers_SetUserParam(...));
        let rc = timers_SetUserParam(
            &mut w.timers,
            ifChangeIntParam,
            intParam,
            ifChangeFloatParam,
            floatParam,
        );
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // PDHG_PrintPDHGParam(w); // prints only, dropped
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::settings_Alloc (lines 813-832)
pub fn settings_Alloc(settings: &mut CUPDLPsettings) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    // settings->nIterLim = INFINITY;
    settings.nIterLim = i32::MAX; // INFINITY cause bug on MacOS
    settings.nLogLevel = 2; // Ensures that, by default, cuPDLP-C printing is unchanged
    settings.nLogInterval = 100;
    // settings->dTimeLim = INFINITY;
    settings.dTimeLim = 3600.0;
    settings.ifScaling = true as cupdlp_int;
    settings.iScalingMethod = 3; // no use
    settings.dScalingLimit = 5.0; // no use
    settings.eRestartMethod = PDHG_GPU_RESTART;

    // termination criteria
    settings.dPrimalTol = 1e-4;
    settings.dDualTol = 1e-4;
    settings.dGapTol = 1e-4;

    retcode
}

/// HiGHS: cupdlp_utils.c::resobj_Alloc (lines 834-932)
///
/// CPU build: every residual / ray vector is allocated (the `NULL`
/// branches are `!defined(CUPDLP_CPU) && USE_KERNELS`). Ends with
/// `cupdlp_filterlb(dLowerFiltered, lower, -INFINITY, ncols)` and
/// `cupdlp_filterub(dUpperFiltered, upper, +INFINITY, ncols)`.
pub fn resobj_Alloc(
    resobj: &mut CUPDLPresobj,
    problem: &CUPDLPproblem,
    ncols: cupdlp_int,
    nrows: cupdlp_int,
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    // #if !defined(CUPDLP_CPU) && USE_KERNELS — GPU branch dropped; CPU branch kept:
    resobj.primalResidual = vec![0.0; nrows as usize];
    resobj.dualResidual = vec![0.0; ncols as usize];
    resobj.primalResidualAverage = vec![0.0; nrows as usize];
    resobj.dualResidualAverage = vec![0.0; ncols as usize];

    resobj.dSlackPos = vec![0.0; ncols as usize];
    resobj.dSlackNeg = vec![0.0; ncols as usize];
    resobj.dSlackPosAverage = vec![0.0; ncols as usize];
    resobj.dSlackNegAverage = vec![0.0; ncols as usize];
    resobj.dLowerFiltered = vec![0.0; ncols as usize];
    resobj.dUpperFiltered = vec![0.0; ncols as usize];

    // #if !defined(CUPDLP_CPU) && USE_KERNELS — GPU branch dropped; CPU branch kept:
    resobj.primalInfeasRay = vec![0.0; ncols as usize];
    resobj.primalInfeasConstr = vec![0.0; nrows as usize];
    resobj.primalInfeasBound = vec![0.0; ncols as usize];
    resobj.dualInfeasRay = vec![0.0; nrows as usize];
    resobj.dualInfeasLbRay = vec![0.0; ncols as usize];
    resobj.dualInfeasUbRay = vec![0.0; ncols as usize];

    // #if !defined(CUPDLP_CPU) && USE_KERNELS — GPU branch dropped; CPU branch kept:
    resobj.dualInfeasConstr = vec![0.0; ncols as usize];

    // CUPDLP_INIT_DOUBLE_ZERO_VEC(resobj->dualInfeasBound, nrows); // commented out upstream

    cupdlp_filterlb(
        &mut resobj.dLowerFiltered,
        &problem.lower,
        -f64::INFINITY,
        ncols,
    );
    cupdlp_filterub(
        &mut resobj.dUpperFiltered,
        &problem.upper,
        f64::INFINITY,
        ncols,
    );

    // initialization
    resobj.dFeasTol = 1e-8;
    resobj.dPrimalObj = 0.0;
    resobj.dDualObj = 0.0;
    resobj.dDualityGap = 0.0;
    resobj.dComplementarity = 0.0;
    resobj.dPrimalFeas = 0.0;
    resobj.dDualFeas = 0.0;
    resobj.dRelObjGap = 0.0;
    resobj.dPrimalObjAverage = 0.0;
    resobj.dDualObjAverage = 0.0;
    resobj.dDualityGapAverage = 0.0;
    resobj.dComplementarityAverage = 0.0;
    resobj.dPrimalFeasAverage = 0.0;
    resobj.dDualFeasAverage = 0.0;
    resobj.dRelObjGapAverage = 0.0;
    resobj.dPrimalFeasLastRestart = 0.0;
    resobj.dDualFeasLastRestart = 0.0;
    resobj.dDualityGapLastRestart = 0.0;
    resobj.dPrimalFeasLastCandidate = 0.0;
    resobj.dDualFeasLastCandidate = 0.0;
    resobj.dDualityGapLastCandidate = 0.0;

    resobj.primalCode = FEASIBLE;
    resobj.dualCode = FEASIBLE;
    resobj.termInfeasIterate = LAST_ITERATE;
    resobj.dPrimalInfeasObj = 0.0;
    resobj.dDualInfeasObj = 0.0;
    resobj.dPrimalInfeasRes = 1.0;
    resobj.dDualInfeasRes = 1.0;
    resobj.dPrimalInfeasObjAverage = 0.0;
    resobj.dDualInfeasObjAverage = 0.0;
    resobj.dPrimalInfeasResAverage = 1.0;
    resobj.dDualInfeasResAverage = 1.0;

    resobj.termCode = TIMELIMIT_OR_ITERLIMIT;
    resobj.termIterate = LAST_ITERATE;

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::iterates_Alloc (lines 934-979)
pub fn iterates_Alloc(
    iterates: &mut CUPDLPiterates,
    ncols: cupdlp_int,
    nrows: cupdlp_int,
) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;

    iterates.nCols = ncols;
    iterates.nRows = nrows;

    iterates.xSum = vec![0.0; ncols as usize];
    iterates.ySum = vec![0.0; nrows as usize];
    iterates.xLastRestart = vec![0.0; ncols as usize];
    iterates.yLastRestart = vec![0.0; nrows as usize];

    // CUPDLP_INIT_CUPDLP_VEC(iterates->x[0], 1); // already default
    // CUPDLP_INIT_CUPDLP_VEC(iterates->x[1], 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->xAverage, 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->y[0], 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->y[1], 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->yAverage, 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->ax[0], 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->ax[1], 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->axAverage, 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->aty[0], 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->aty[1], 1);
    // CUPDLP_INIT_CUPDLP_VEC(iterates->atyAverage, 1);

    'exit_cleanup: {
        // CUPDLP_CALL(vec_Alloc(iterates->x[0], ncols));
        let rc = vec_Alloc(&mut iterates.x[0], ncols);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->x[1], ncols));
        let rc = vec_Alloc(&mut iterates.x[1], ncols);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->xAverage, ncols));
        let rc = vec_Alloc(&mut iterates.xAverage, ncols);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->y[0], nrows));
        let rc = vec_Alloc(&mut iterates.y[0], nrows);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->y[1], nrows));
        let rc = vec_Alloc(&mut iterates.y[1], nrows);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->yAverage, nrows));
        let rc = vec_Alloc(&mut iterates.yAverage, nrows);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->ax[0], nrows));
        let rc = vec_Alloc(&mut iterates.ax[0], nrows);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->ax[1], nrows));
        let rc = vec_Alloc(&mut iterates.ax[1], nrows);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->axAverage, nrows));
        let rc = vec_Alloc(&mut iterates.axAverage, nrows);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->aty[0], ncols));
        let rc = vec_Alloc(&mut iterates.aty[0], ncols);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->aty[1], ncols));
        let rc = vec_Alloc(&mut iterates.aty[1], ncols);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(vec_Alloc(iterates->atyAverage, ncols));
        let rc = vec_Alloc(&mut iterates.atyAverage, ncols);
        if rc != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // initialization
        iterates.iLastRestartIter = 0;
        iterates.dLastRestartDualityGap = 0.0;
        iterates.dLastRestartBeta = 0.0;
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::stepsize_Alloc (lines 981-997)
pub fn stepsize_Alloc(stepsize: &mut CUPDLPstepsize) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    stepsize.eLineSearchMethod = PDHG_ADAPTIVE_LINESEARCH;

    // initialization
    stepsize.nStepSizeIter = 0;
    stepsize.dPrimalStep = 0.0;
    stepsize.dDualStep = 0.0;
    stepsize.dSumPrimalStep = 0.0;
    stepsize.dSumDualStep = 0.0;
    stepsize.dBeta = 0.0;
    stepsize.dTheta = 0.0;

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::timers_Alloc (lines 999-1030)
pub fn timers_Alloc(timers: &mut CUPDLPtimers) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    timers.nIter = 0;
    timers.dSolvingTime = 0.0;
    timers.dSolvingBeg = 0.0;
    timers.dScalingTime = 0.0;
    timers.dPresolveTime = 0.0;

    // #if PDHG_USE_TIMERS (== 1 in CPU build)
    timers.dAtyTime = 0.0;
    timers.dAxTime = 0.0;
    timers.dComputeResidualsTime = 0.0;
    timers.dUpdateIterateTime = 0.0;
    timers.nAtyCalls = 0;
    timers.nAxCalls = 0;
    timers.nComputeResidualsCalls = 0;
    timers.nUpdateIterateCalls = 0;
    // #endif

    // #ifndef CUPDLP_CPU — GPU timers dropped (lines 1018-1026)

    // exit_cleanup:
    retcode
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn settings_alloc_defaults_then_user_params() {
        let mut s = CUPDLPsettings::default();
        assert_eq!(settings_Alloc(&mut s), RETCODE_OK);
        assert_eq!(s.nIterLim, i32::MAX);
        assert_eq!(s.nLogLevel, 2);
        assert_eq!(s.nLogInterval, 100);
        assert_eq!(s.dTimeLim, 3600.0);
        assert_eq!(s.ifScaling, 1);
        assert_eq!(s.eRestartMethod, PDHG_GPU_RESTART);
        assert_eq!(s.dPrimalTol, 1e-4);
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        ici[N_ITER_LIM] = true;
        ip[N_ITER_LIM] = 77;
        icf[D_GAP_TOL] = true;
        fp[D_GAP_TOL] = 1e-7;
        assert_eq!(
            settings_SetUserParam(&mut s, &ici, &ip, &icf, &fp),
            RETCODE_OK
        );
        assert_eq!(s.nIterLim, 77);
        assert_eq!(s.dGapTol, 1e-7);
        assert_eq!(s.nLogLevel, 2);
    }

    #[test]
    fn resobj_and_iterates_alloc() {
        let mut p = CUPDLPproblem::default();
        p.lower = vec![f64::NEG_INFINITY, 1.0];
        p.upper = vec![2.0, f64::INFINITY];
        let mut r = CUPDLPresobj::default();
        assert_eq!(resobj_Alloc(&mut r, &p, 2, 3), RETCODE_OK);
        assert_eq!(r.dFeasTol, 1e-8);
        assert_eq!(r.dPrimalInfeasRes, 1.0);
        assert_eq!(r.termCode, TIMELIMIT_OR_ITERLIMIT);
        assert_eq!(r.primalCode, FEASIBLE);
        assert_eq!(r.dLowerFiltered, vec![0.0, 1.0]);
        assert_eq!(r.dUpperFiltered, vec![2.0, 0.0]);
        assert_eq!(r.primalResidual.len(), 3);
        let mut it = CUPDLPiterates::default();
        assert_eq!(iterates_Alloc(&mut it, 2, 3), RETCODE_OK);
        assert_eq!(it.x[1].len, 2);
        assert_eq!(it.ax[0].data.len(), 3);
        assert_eq!(it.atyAverage.len, 2);
        assert_eq!(it.xSum.len(), 2);
        assert_eq!(it.ySum.len(), 3);
        assert_eq!(it.nCols, 2);
        assert_eq!(it.nRows, 3);
    }

    #[test]
    fn pdhg_set_user_param_and_other_allocs() {
        let mut w = CUPDLPwork::default();
        assert_eq!(stepsize_Alloc(&mut w.stepsize), RETCODE_OK);
        assert_eq!(w.stepsize.eLineSearchMethod, PDHG_ADAPTIVE_LINESEARCH);
        assert_eq!(timers_Alloc(&mut w.timers), RETCODE_OK);
        assert_eq!(w.timers.nIter, 0);
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        ici[E_LINE_SEARCH_METHOD] = true;
        ip[E_LINE_SEARCH_METHOD] = 0;
        ici[IF_PC_SCALING] = true;
        ip[IF_PC_SCALING] = 1;
        icf[D_FEAS_TOL] = true;
        fp[D_FEAS_TOL] = 1e-6;
        icf[D_TIME_LIM] = true;
        fp[D_TIME_LIM] = 12.5;
        assert_eq!(PDHG_SetUserParam(&mut w, &ici, &ip, &icf, &fp), RETCODE_OK);
        assert_eq!(w.stepsize.eLineSearchMethod, 0);
        assert_eq!(w.scaling.ifPcScaling, 1);
        assert_eq!(w.resobj.dFeasTol, 1e-6);
        assert_eq!(w.settings.dTimeLim, 12.5);
    }
}
