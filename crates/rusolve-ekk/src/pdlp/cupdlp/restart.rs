//! HiGHS source: pdlp/cupdlp/cupdlp_restart.h, pdlp/cupdlp/cupdlp_restart.c
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The adaptive restart test of PDLP (the "GPU" restart scheme, which is
//! the one the CPU build runs too). `PDHG_restart_choice` is in
//! [`super::defs`]. The `PDLP_DEBUG_LOG` blocks are not compiled.

use super::defs::*;

/// HiGHS: cupdlp_restart.c::PDHG_Check_Restart_GPU (lines 3-99)
pub fn PDHG_Check_Restart_GPU(work: &mut CUPDLPwork) -> PDHG_restart_choice {
    // CUPDLPproblem *problem = work->problem;
    // CUPDLPdata *lp = problem->data;
    // CUPDLPstepsize *stepsize = work->stepsize;
    // CUPDLPiterates *iterates = work->iterates;
    // CUPDLPresobj *resobj = work->resobj;
    // CUPDLPtimers *timers = work->timers;

    // Is it first time called?
    if work.timers.nIter == work.iterates.iLastRestartIter {
        work.resobj.dPrimalFeasLastRestart = work.resobj.dPrimalFeas;
        work.resobj.dDualFeasLastRestart = work.resobj.dDualFeas;
        work.resobj.dDualityGapLastRestart = work.resobj.dDualityGap;

        work.resobj.dPrimalFeasLastCandidate = work.resobj.dPrimalFeas;
        work.resobj.dDualFeasLastCandidate = work.resobj.dDualFeas;
        work.resobj.dDualityGapLastCandidate = work.resobj.dDualityGap;

        return PDHG_NO_RESTART;
    }

    let muCurrent = PDHG_Restart_Score_GPU(
        work.stepsize.dBeta,
        work.resobj.dPrimalFeas,
        work.resobj.dDualFeas,
        work.resobj.dDualityGap,
    );
    let muAverage = PDHG_Restart_Score_GPU(
        work.stepsize.dBeta,
        work.resobj.dPrimalFeasAverage,
        work.resobj.dDualFeasAverage,
        work.resobj.dDualityGapAverage,
    );
    // #if PDLP_DEBUG_LOG ... #endif (not compiled)

    let mut muCandidate: cupdlp_float;
    let mut restart_choice = PDHG_RESTART_TO_AVERAGE;
    if muCurrent < muAverage {
        restart_choice = PDHG_RESTART_TO_CURRENT;
        muCandidate = muCurrent;
    } else {
        // restart_choice = PDHG_RESTART_TO_AVERAGE;
        muCandidate = muAverage;
    }

    // Or should we do artificial restart based on iteration count?
    // if (2 * (timers->nIter - iterates->iLastRestartIter) >= timers->nIter) {
    if ((work.timers.nIter - work.iterates.iLastRestartIter) as f64)
        >= 0.36 * (work.timers.nIter as f64)
    {
        // #if CUPDLP_DEBUG ... #endif (not compiled)
    } else {
        let muLastRestart = PDHG_Restart_Score_GPU(
            work.stepsize.dBeta,
            work.resobj.dPrimalFeasLastRestart,
            work.resobj.dDualFeasLastRestart,
            work.resobj.dDualityGapLastRestart,
        );

        // Sufficient decay
        if muCandidate < 0.2 * muLastRestart {
            // #if CUPDLP_DEBUG ... #endif (not compiled)
        } else {
            let muLastCandidate = PDHG_Restart_Score_GPU(
                work.stepsize.dBeta,
                work.resobj.dPrimalFeasLastCandidate,
                work.resobj.dDualFeasLastCandidate,
                work.resobj.dDualityGapLastCandidate,
            );

            // Necessary decay
            if muCandidate < 0.8 * muLastRestart && muCandidate > muLastCandidate {
                // #if CUPDLP_DEBUG ... #endif (not compiled)
            } else {
                restart_choice = PDHG_NO_RESTART;
            }
        }
    }

    // record candidate
    if muCurrent < muAverage {
        work.resobj.dPrimalFeasLastCandidate = work.resobj.dPrimalFeas;
        work.resobj.dDualFeasLastCandidate = work.resobj.dDualFeas;
        work.resobj.dDualityGapLastCandidate = work.resobj.dDualityGap;
    } else {
        work.resobj.dPrimalFeasLastCandidate = work.resobj.dPrimalFeasAverage;
        work.resobj.dDualFeasLastCandidate = work.resobj.dDualFeasAverage;
        work.resobj.dDualityGapLastCandidate = work.resobj.dDualityGapAverage;
    }

    // if (restart_choice != PDHG_NO_RESTART) { ... logging ... }
    // Logging dropped.

    restart_choice
}

/// HiGHS: cupdlp_restart.c::PDHG_Check_Restart_Merit_Function (lines 101-111)
pub fn PDHG_Check_Restart_Merit_Function(work: &CUPDLPwork) -> cupdlp_bool {
    // CUPDLPproblem *problem = work->problem;
    // CUPDLPdata *lp = problem->data;
    // CUPDLPstepsize *stepsize = work->stepsize;
    // CUPDLPiterates *iterates = work->iterates;
    // CUPDLPresobj *resobj = work->resobj;

    (work.resobj.dDualityGap.abs() > 2.0 * work.resobj.dDualityGapAverage.abs())
        && (work.resobj.dPrimalFeas > 2.0 * work.resobj.dPrimalFeasAverage)
}

/// HiGHS: cupdlp_restart.c::PDHG_Restart_Score_GPU (lines 113-124)
pub fn PDHG_Restart_Score_GPU(
    weightSquared: cupdlp_float,
    dPrimalFeas: cupdlp_float,
    dDualFeas: cupdlp_float,
    dDualityGap: cupdlp_float,
) -> cupdlp_float {
    (weightSquared * dPrimalFeas * dPrimalFeas
        + dDualFeas * dDualFeas / weightSquared
        + dDualityGap * dDualityGap)
        .sqrt()
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn score_formula() {
        // sqrt(4*1*1 + 2*2/4 + 0) = sqrt(5)
        assert_eq!(PDHG_Restart_Score_GPU(4.0, 1.0, 2.0, 0.0), 5.0f64.sqrt());
    }

    #[test]
    fn first_call_records_last_restart_and_does_not_restart() {
        let mut w = CUPDLPwork::default();
        w.timers.nIter = 3;
        w.iterates.iLastRestartIter = 3;
        w.resobj.dPrimalFeas = 1.5;
        w.resobj.dDualityGap = 2.5;
        assert_eq!(PDHG_Check_Restart_GPU(&mut w), PDHG_NO_RESTART);
        assert_eq!(w.resobj.dPrimalFeasLastRestart, 1.5);
        assert_eq!(w.resobj.dDualityGapLastCandidate, 2.5);
    }
}
