//! HiGHS source: pdlp/cupdlp/cupdlp_proj.h, pdlp/cupdlp/cupdlp_proj.c
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Bound projection and the restart of the iterates. Lines 26-73 are
//! commented-out code upstream.

use super::defs::*;
use super::linalg::*;
use super::restart::*;
use super::solver::PDHG_Compute_Residuals;
use super::step::PDHG_Compute_Step_Size_Ratio;

/// HiGHS: cupdlp_proj.c::PDHG_Project_Bounds (lines 17-25)
///
/// `work` replaced by `problem`: the callers
/// (`PDHG_Init_Variables`, `PDHG_Compute_Average_Iterate`) pass `iterates`
/// vectors as `r`. It only calls `cupdlp_projub` /
/// `cupdlp_projlb`; its unit test
/// is `linalg::part2::tests::project_bounds_clamps_to_box`.
pub fn PDHG_Project_Bounds(problem: &CUPDLPproblem, r: &mut [cupdlp_float]) {
    cupdlp_projub(r, &problem.upper, problem.nCols);
    cupdlp_projlb(r, &problem.lower, problem.nCols);
}

/// HiGHS: cupdlp_proj.c::PDHG_Restart_Iterate (lines 75-86)
///
/// `PDHG_CPU_RESTART` is a `TODO` upstream and does nothing.
pub fn PDHG_Restart_Iterate(pdhg: &mut CUPDLPwork) {
    match pdhg.settings.eRestartMethod {
        PDHG_WITHOUT_RESTART => {
            // break
        }
        PDHG_GPU_RESTART => {
            PDHG_Restart_Iterate_GPU(pdhg);
        }
        PDHG_CPU_RESTART => {
            // TODO: implement PDHG_Restart_Iterate_CPU(pdhg);
        }
        _ => {}
    }
}

/// HiGHS: cupdlp_proj.c::PDHG_Restart_Iterate_GPU (lines 88-148)
///
/// The `PDLP_DEBUG_LOG` block is not compiled. The oracle hook
/// `pdlp_restart(pdhg, restart_choice)` fires as the last statement (after the final
/// `PDHG_Compute_Residuals`), i.e. only when a restart happened.
pub fn PDHG_Restart_Iterate_GPU(pdhg: &mut CUPDLPwork) {
    // CUPDLPproblem *problem = pdhg->problem;
    // CUPDLPiterates *iterates = pdhg->iterates;
    // CUPDLPstepsize *stepsize = pdhg->stepsize;
    // CUPDLPresobj *resobj = pdhg->resobj;
    // CUPDLPtimers *timers = pdhg->timers;

    // PDHG_Compute_Average_Iterate(pdhg);
    let restart_choice = PDHG_Check_Restart_GPU(pdhg);

    if restart_choice == PDHG_NO_RESTART {
        return;
    }

    let iter = pdhg.timers.nIter;
    let idx = (iter % 2) as usize;

    // CUPDLPvec *x = iterates->x[iter % 2];  (and y, ax, aty likewise):
    // accessed below as `pdhg.iterates.{x,y,ax,aty}[idx]`.

    pdhg.stepsize.dSumPrimalStep = 0.0;
    pdhg.stepsize.dSumDualStep = 0.0;

    // CUPDLP_ZERO_VEC(iterates->xSum, cupdlp_float, problem->nCols);
    cupdlp_zero(&mut pdhg.iterates.xSum, pdhg.problem.nCols);
    // CUPDLP_ZERO_VEC(iterates->ySum, cupdlp_float, problem->nRows);
    cupdlp_zero(&mut pdhg.iterates.ySum, pdhg.problem.nRows);

    if restart_choice == PDHG_RESTART_TO_AVERAGE {
        pdhg.resobj.dPrimalFeasLastRestart = pdhg.resobj.dPrimalFeasAverage;
        pdhg.resobj.dDualFeasLastRestart = pdhg.resobj.dDualFeasAverage;
        pdhg.resobj.dDualityGapLastRestart = pdhg.resobj.dDualityGapAverage;

        // CUPDLP_COPY_VEC(x->data, iterates->xAverage->data, cupdlp_float, problem->nCols);
        cupdlp_copy(
            &mut pdhg.iterates.x[idx].data,
            &pdhg.iterates.xAverage.data,
            pdhg.problem.nCols,
        );
        // CUPDLP_COPY_VEC(y->data, iterates->yAverage->data, cupdlp_float, problem->nRows);
        cupdlp_copy(
            &mut pdhg.iterates.y[idx].data,
            &pdhg.iterates.yAverage.data,
            pdhg.problem.nRows,
        );
        // CUPDLP_COPY_VEC(ax->data, iterates->axAverage->data, cupdlp_float, problem->nRows);
        cupdlp_copy(
            &mut pdhg.iterates.ax[idx].data,
            &pdhg.iterates.axAverage.data,
            pdhg.problem.nRows,
        );
        // CUPDLP_COPY_VEC(aty->data, iterates->atyAverage->data, cupdlp_float, problem->nCols);
        cupdlp_copy(
            &mut pdhg.iterates.aty[idx].data,
            &pdhg.iterates.atyAverage.data,
            pdhg.problem.nCols,
        );
    } else {
        pdhg.resobj.dPrimalFeasLastRestart = pdhg.resobj.dPrimalFeas;
        pdhg.resobj.dDualFeasLastRestart = pdhg.resobj.dDualFeas;
        pdhg.resobj.dDualityGapLastRestart = pdhg.resobj.dDualityGap;
    }

    PDHG_Compute_Step_Size_Ratio(pdhg);

    // CUPDLP_COPY_VEC(iterates->xLastRestart, x->data, cupdlp_float, problem->nCols);
    cupdlp_copy(
        &mut pdhg.iterates.xLastRestart,
        &pdhg.iterates.x[idx].data,
        pdhg.problem.nCols,
    );
    // CUPDLP_COPY_VEC(iterates->yLastRestart, y->data, cupdlp_float, problem->nRows);
    cupdlp_copy(
        &mut pdhg.iterates.yLastRestart,
        &pdhg.iterates.y[idx].data,
        pdhg.problem.nRows,
    );

    pdhg.iterates.iLastRestartIter = pdhg.timers.nIter;

    PDHG_Compute_Residuals(pdhg);

    // Oracle hook: fires only when a restart happened (after the early return).
    super::super::oracle_trace::pdlp_restart(pdhg, restart_choice);
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn restart_without_restart_is_a_no_op() {
        let mut w = CUPDLPwork::default();
        w.settings.eRestartMethod = PDHG_WITHOUT_RESTART;
        w.iterates.iLastRestartIter = 7;
        PDHG_Restart_Iterate(&mut w);
        assert_eq!(w.iterates.iLastRestartIter, 7);
    }
}
