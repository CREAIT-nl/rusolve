//! HiGHS source: pdlp/cupdlp/cupdlp_solver.c (lines 313-591)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Dual infeasibility, the residual drivers and the iterate initialisation.
//! The drivers call the functions of [`super`] with vectors moved
//! out of `work` and back (see the module doc of [`super`]).

use super::super::defs::*;
use super::super::linalg::*;
use super::super::proj::PDHG_Project_Bounds;
use super::*;

/// HiGHS: cupdlp_solver.c::PDHG_Compute_Dual_Infeasibility (lines 313-429)
///
/// CPU branch. Writes `resobj->primalInfeasRay` / `primalInfeasConstr` /
/// `primalInfeasBound` through `work`.
pub fn PDHG_Compute_Dual_Infeasibility(
    work: &mut CUPDLPwork,
    x: &[cupdlp_float],
    ax: &[cupdlp_float],
    primalObj: cupdlp_float,
    dDualInfeasObj: &mut cupdlp_float,
    dDualInfeasRes: &mut cupdlp_float,
) {
    let nCols = work.problem.data.nCols;
    let nRows = work.problem.data.nRows;
    let nEqs = work.problem.nEqs;
    let ifScaled = work.scaling.ifScaled;

    let mut pScale: cupdlp_float = 1.0;
    let mut pConstrResSq: cupdlp_float = 0.0;
    let mut pBoundLbResSq: cupdlp_float = 0.0;
    let mut pBoundUbResSq: cupdlp_float = 0.0;

    // GPU branch dropped (lines 326-367: #if !defined(CUPDLP_CPU) && USE_KERNELS)

    // CPU branch (lines 368-428):

    // x ray
    // CUPDLP_COPY_VEC(resobj->primalInfeasRay, x, cupdlp_float, problem->data->nCols);
    // cupdlp_twoNorm(work, problem->data->nCols, resobj->primalInfeasRay, &pScale);
    let mut primalInfeasRay = std::mem::take(&mut work.resobj.primalInfeasRay);
    cupdlp_copy(&mut primalInfeasRay, x, nCols);
    cupdlp_twoNorm(nCols, &primalInfeasRay, &mut pScale);
    // pScale /= sqrt(problem->data->nCols);
    if pScale < 1e-8 {
        pScale = 1.0;
    }
    // cupdlp_scaleVector(work, 1.0 / pScale, resobj->primalInfeasRay, problem->data->nCols);
    cupdlp_scaleVector(1.0 / pScale, &mut primalInfeasRay, nCols);

    // primal obj
    *dDualInfeasObj = (primalObj - work.problem.offset) / work.problem.sense_origin / pScale;

    // primal constraints [Ax, min(Gx, 0)]
    // CUPDLP_COPY_VEC(resobj->primalInfeasConstr, ax, cupdlp_float, problem->data->nRows);
    let mut primalInfeasConstr = std::mem::take(&mut work.resobj.primalInfeasConstr);
    cupdlp_copy(&mut primalInfeasConstr, ax, nRows);
    // cupdlp_scaleVector(work, 1.0 / pScale, resobj->primalInfeasConstr, problem->data->nRows);
    cupdlp_scaleVector(1.0 / pScale, &mut primalInfeasConstr, nRows);
    // cupdlp_projNeg(resobj->primalInfeasConstr + problem->nEqs, problem->data->nRows - problem->nEqs);
    {
        let nEqs_usize = nEqs as usize;
        let len = (nRows - nEqs) as usize;
        cupdlp_projNeg(
            &mut primalInfeasConstr[nEqs_usize..nEqs_usize + len],
            nRows - nEqs,
        );
    }
    if ifScaled != 0 {
        // cupdlp_edot(resobj->primalInfeasConstr, work->rowScale, problem->data->nRows);
        cupdlp_edot(&mut primalInfeasConstr, &work.rowScale, nRows);
    }
    // cupdlp_twoNormSquared(work, problem->data->nRows, resobj->primalInfeasConstr, &pConstrResSq);
    cupdlp_twoNormSquared(nRows, &primalInfeasConstr, &mut pConstrResSq);

    // primal bound
    // lb
    // CUPDLP_COPY_VEC(resobj->primalInfeasBound, resobj->primalInfeasRay, cupdlp_float, problem->data->nCols);
    let mut primalInfeasBound = std::mem::take(&mut work.resobj.primalInfeasBound);
    cupdlp_copy(&mut primalInfeasBound, &primalInfeasRay, nCols);
    // cupdlp_projNeg(resobj->primalInfeasBound, problem->data->nCols);
    cupdlp_projNeg(&mut primalInfeasBound, nCols);
    // cupdlp_edot(resobj->primalInfeasBound, problem->hasLower, problem->data->nCols);
    cupdlp_edot(&mut primalInfeasBound, &work.problem.hasLower, nCols);
    if ifScaled != 0 {
        // cupdlp_ediv(resobj->primalInfeasBound, work->colScale, problem->data->nCols);
        cupdlp_ediv(&mut primalInfeasBound, &work.colScale, nCols);
    }
    // cupdlp_twoNormSquared(work, problem->data->nCols, resobj->primalInfeasBound, &pBoundLbResSq);
    cupdlp_twoNormSquared(nCols, &primalInfeasBound, &mut pBoundLbResSq);

    // ub
    // CUPDLP_COPY_VEC(resobj->primalInfeasBound, resobj->primalInfeasRay, cupdlp_float, problem->data->nCols);
    cupdlp_copy(&mut primalInfeasBound, &primalInfeasRay, nCols);
    // cupdlp_projPos(resobj->primalInfeasBound, problem->data->nCols);
    cupdlp_projPos(&mut primalInfeasBound, nCols);
    // cupdlp_edot(resobj->primalInfeasBound, problem->hasUpper, problem->data->nCols);
    cupdlp_edot(&mut primalInfeasBound, &work.problem.hasUpper, nCols);
    if ifScaled != 0 {
        // cupdlp_ediv(resobj->primalInfeasBound, work->colScale, problem->data->nCols);
        cupdlp_ediv(&mut primalInfeasBound, &work.colScale, nCols);
    }
    // cupdlp_twoNormSquared(work, problem->data->nCols, resobj->primalInfeasBound, &pBoundUbResSq);
    cupdlp_twoNormSquared(nCols, &primalInfeasBound, &mut pBoundUbResSq);

    // sum up
    *dDualInfeasRes = (pConstrResSq + pBoundLbResSq + pBoundUbResSq).sqrt();

    // Put vectors back.
    work.resobj.primalInfeasRay = primalInfeasRay;
    work.resobj.primalInfeasConstr = primalInfeasConstr;
    work.resobj.primalInfeasBound = primalInfeasBound;
}

/// HiGHS: cupdlp_solver.c::PDHG_Compute_Infeas_Residuals (lines 433-471)
///
/// Current iterate first (`x[nIter % 2]` etc.), then the average. The
/// `#if problem_USE_TIMERS` blocks are dropped: that macro is never
/// defined (only `PDHG_USE_TIMERS` is).
pub fn PDHG_Compute_Infeas_Residuals(work: &mut CUPDLPwork) {
    // #if problem_USE_TIMERS blocks dropped (lines 434-437, 468-470):
    // problem_USE_TIMERS is never defined; only PDHG_USE_TIMERS is.

    let iter = work.timers.nIter;
    let idx = (iter % 2) as usize;

    // current solution
    // PDHG_Compute_Primal_Infeasibility(work, y->data, resobj->dSlackPos,
    //                                    resobj->dSlackNeg, aty->data,
    //                                    resobj->dDualObj, &resobj->dPrimalInfeasObj,
    //                                    &resobj->dPrimalInfeasRes);
    {
        let y_vec = std::mem::take(&mut work.iterates.y[idx]);
        let aty_vec = std::mem::take(&mut work.iterates.aty[idx]);
        let dSlackPos = std::mem::take(&mut work.resobj.dSlackPos);
        let dSlackNeg = std::mem::take(&mut work.resobj.dSlackNeg);
        let dualObj = work.resobj.dDualObj;
        let mut dPrimalInfeasObj: cupdlp_float = work.resobj.dPrimalInfeasObj;
        let mut dPrimalInfeasRes: cupdlp_float = work.resobj.dPrimalInfeasRes;
        PDHG_Compute_Primal_Infeasibility(
            work,
            &y_vec.data,
            &dSlackPos,
            &dSlackNeg,
            &aty_vec.data,
            dualObj,
            &mut dPrimalInfeasObj,
            &mut dPrimalInfeasRes,
        );
        work.iterates.y[idx] = y_vec;
        work.iterates.aty[idx] = aty_vec;
        work.resobj.dSlackPos = dSlackPos;
        work.resobj.dSlackNeg = dSlackNeg;
        work.resobj.dPrimalInfeasObj = dPrimalInfeasObj;
        work.resobj.dPrimalInfeasRes = dPrimalInfeasRes;
    }

    // PDHG_Compute_Dual_Infeasibility(work, x->data, ax->data,
    //                                  resobj->dPrimalObj, &resobj->dDualInfeasObj,
    //                                  &resobj->dDualInfeasRes);
    {
        let x_vec = std::mem::take(&mut work.iterates.x[idx]);
        let ax_vec = std::mem::take(&mut work.iterates.ax[idx]);
        let primalObj = work.resobj.dPrimalObj;
        let mut dDualInfeasObj: cupdlp_float = work.resobj.dDualInfeasObj;
        let mut dDualInfeasRes: cupdlp_float = work.resobj.dDualInfeasRes;
        PDHG_Compute_Dual_Infeasibility(
            work,
            &x_vec.data,
            &ax_vec.data,
            primalObj,
            &mut dDualInfeasObj,
            &mut dDualInfeasRes,
        );
        work.iterates.x[idx] = x_vec;
        work.iterates.ax[idx] = ax_vec;
        work.resobj.dDualInfeasObj = dDualInfeasObj;
        work.resobj.dDualInfeasRes = dDualInfeasRes;
    }

    // average solution
    // PDHG_Compute_Primal_Infeasibility(
    //     work, iterates->yAverage->data, resobj->dSlackPosAverage,
    //     resobj->dSlackNegAverage, iterates->atyAverage->data,
    //     resobj->dDualObjAverage, &resobj->dPrimalInfeasObjAverage,
    //     &resobj->dPrimalInfeasResAverage);
    {
        let yAverage_vec = std::mem::take(&mut work.iterates.yAverage);
        let atyAverage_vec = std::mem::take(&mut work.iterates.atyAverage);
        let dSlackPosAverage = std::mem::take(&mut work.resobj.dSlackPosAverage);
        let dSlackNegAverage = std::mem::take(&mut work.resobj.dSlackNegAverage);
        let dualObjAverage = work.resobj.dDualObjAverage;
        let mut dPrimalInfeasObjAverage: cupdlp_float = work.resobj.dPrimalInfeasObjAverage;
        let mut dPrimalInfeasResAverage: cupdlp_float = work.resobj.dPrimalInfeasResAverage;
        PDHG_Compute_Primal_Infeasibility(
            work,
            &yAverage_vec.data,
            &dSlackPosAverage,
            &dSlackNegAverage,
            &atyAverage_vec.data,
            dualObjAverage,
            &mut dPrimalInfeasObjAverage,
            &mut dPrimalInfeasResAverage,
        );
        work.iterates.yAverage = yAverage_vec;
        work.iterates.atyAverage = atyAverage_vec;
        work.resobj.dSlackPosAverage = dSlackPosAverage;
        work.resobj.dSlackNegAverage = dSlackNegAverage;
        work.resobj.dPrimalInfeasObjAverage = dPrimalInfeasObjAverage;
        work.resobj.dPrimalInfeasResAverage = dPrimalInfeasResAverage;
    }

    // PDHG_Compute_Dual_Infeasibility(
    //     work, iterates->xAverage->data, iterates->axAverage->data,
    //     resobj->dPrimalObjAverage, &resobj->dDualInfeasObjAverage,
    //     &resobj->dDualInfeasResAverage);
    {
        let xAverage_vec = std::mem::take(&mut work.iterates.xAverage);
        let axAverage_vec = std::mem::take(&mut work.iterates.axAverage);
        let primalObjAverage = work.resobj.dPrimalObjAverage;
        let mut dDualInfeasObjAverage: cupdlp_float = work.resobj.dDualInfeasObjAverage;
        let mut dDualInfeasResAverage: cupdlp_float = work.resobj.dDualInfeasResAverage;
        PDHG_Compute_Dual_Infeasibility(
            work,
            &xAverage_vec.data,
            &axAverage_vec.data,
            primalObjAverage,
            &mut dDualInfeasObjAverage,
            &mut dDualInfeasResAverage,
        );
        work.iterates.xAverage = xAverage_vec;
        work.iterates.axAverage = axAverage_vec;
        work.resobj.dDualInfeasObjAverage = dDualInfeasObjAverage;
        work.resobj.dDualInfeasResAverage = dDualInfeasResAverage;
    }
}

/// HiGHS: cupdlp_solver.c::PDHG_Compute_Residuals (lines 473-529)
///
/// Order: primal / dual feasibility of the current iterate, primal / dual
/// feasibility of the average, then the gaps:
/// `dDualityGap = dPrimalObj - dDualObj`,
/// `dRelObjGap = fabs(dPrimalObj - dDualObj) / (1.0 + fabs(dPrimalObj) + fabs(dDualObj))`,
/// and the same for the averages. The `nComputeResidualsCalls` /
/// `dComputeResidualsTime` updates sit under `#if problem_USE_TIMERS`
/// (never defined, not `PDHG_USE_TIMERS`) and are dropped.
///
/// The output locals start from the current `resobj` fields because the C
/// passes `&resobj->...` and the CPU `PDHG_Compute_Dual_Feasibility` never
/// writes `*dComplementarity`: the field keeps its previous value.
pub fn PDHG_Compute_Residuals(work: &mut CUPDLPwork) {
    // #if problem_USE_TIMERS blocks dropped (lines 474-477, 526-528):
    // problem_USE_TIMERS is never defined; only PDHG_USE_TIMERS is.

    let iter = work.timers.nIter;
    let idx = (iter % 2) as usize;

    // PDHG_Compute_Primal_Feasibility(work, resobj->primalResidual,
    //                                  ax->data, x->data,
    //                                  &resobj->dPrimalFeas, &resobj->dPrimalObj);
    {
        let ax_vec = std::mem::take(&mut work.iterates.ax[idx]);
        let x_vec = std::mem::take(&mut work.iterates.x[idx]);
        let mut primalResidual = std::mem::take(&mut work.resobj.primalResidual);
        let mut dPrimalFeas: cupdlp_float = work.resobj.dPrimalFeas;
        let mut dPrimalObj: cupdlp_float = work.resobj.dPrimalObj;
        PDHG_Compute_Primal_Feasibility(
            work,
            &mut primalResidual,
            &ax_vec.data,
            &x_vec.data,
            &mut dPrimalFeas,
            &mut dPrimalObj,
        );
        work.iterates.ax[idx] = ax_vec;
        work.iterates.x[idx] = x_vec;
        work.resobj.primalResidual = primalResidual;
        work.resobj.dPrimalFeas = dPrimalFeas;
        work.resobj.dPrimalObj = dPrimalObj;
    }

    // PDHG_Compute_Dual_Feasibility(
    //     work, resobj->dualResidual, aty->data, x->data,
    //     y->data, &resobj->dDualFeas, &resobj->dDualObj,
    //     &resobj->dComplementarity, resobj->dSlackPos, resobj->dSlackNeg);
    {
        let aty_vec = std::mem::take(&mut work.iterates.aty[idx]);
        let x_vec = std::mem::take(&mut work.iterates.x[idx]);
        let y_vec = std::mem::take(&mut work.iterates.y[idx]);
        let mut dualResidual = std::mem::take(&mut work.resobj.dualResidual);
        let mut dSlackPos = std::mem::take(&mut work.resobj.dSlackPos);
        let mut dSlackNeg = std::mem::take(&mut work.resobj.dSlackNeg);
        let mut dDualFeas: cupdlp_float = work.resobj.dDualFeas;
        let mut dDualObj: cupdlp_float = work.resobj.dDualObj;
        let mut dComplementarity: cupdlp_float = work.resobj.dComplementarity;
        PDHG_Compute_Dual_Feasibility(
            work,
            &mut dualResidual,
            &aty_vec.data,
            &x_vec.data,
            &y_vec.data,
            &mut dDualFeas,
            &mut dDualObj,
            &mut dComplementarity,
            &mut dSlackPos,
            &mut dSlackNeg,
        );
        work.iterates.aty[idx] = aty_vec;
        work.iterates.x[idx] = x_vec;
        work.iterates.y[idx] = y_vec;
        work.resobj.dualResidual = dualResidual;
        work.resobj.dSlackPos = dSlackPos;
        work.resobj.dSlackNeg = dSlackNeg;
        work.resobj.dDualFeas = dDualFeas;
        work.resobj.dDualObj = dDualObj;
        work.resobj.dComplementarity = dComplementarity;
    }

    // PDHG_Compute_Primal_Feasibility(
    //     work, resobj->primalResidualAverage, iterates->axAverage->data,
    //     iterates->xAverage->data, &resobj->dPrimalFeasAverage,
    //     &resobj->dPrimalObjAverage);
    {
        let axAverage_vec = std::mem::take(&mut work.iterates.axAverage);
        let xAverage_vec = std::mem::take(&mut work.iterates.xAverage);
        let mut primalResidualAverage = std::mem::take(&mut work.resobj.primalResidualAverage);
        let mut dPrimalFeasAverage: cupdlp_float = work.resobj.dPrimalFeasAverage;
        let mut dPrimalObjAverage: cupdlp_float = work.resobj.dPrimalObjAverage;
        PDHG_Compute_Primal_Feasibility(
            work,
            &mut primalResidualAverage,
            &axAverage_vec.data,
            &xAverage_vec.data,
            &mut dPrimalFeasAverage,
            &mut dPrimalObjAverage,
        );
        work.iterates.axAverage = axAverage_vec;
        work.iterates.xAverage = xAverage_vec;
        work.resobj.primalResidualAverage = primalResidualAverage;
        work.resobj.dPrimalFeasAverage = dPrimalFeasAverage;
        work.resobj.dPrimalObjAverage = dPrimalObjAverage;
    }

    // PDHG_Compute_Dual_Feasibility(
    //     work, resobj->dualResidualAverage, iterates->atyAverage->data,
    //     iterates->xAverage->data, iterates->yAverage->data,
    //     &resobj->dDualFeasAverage, &resobj->dDualObjAverage,
    //     &resobj->dComplementarityAverage, resobj->dSlackPosAverage,
    //     resobj->dSlackNegAverage);
    {
        let atyAverage_vec = std::mem::take(&mut work.iterates.atyAverage);
        let xAverage_vec = std::mem::take(&mut work.iterates.xAverage);
        let yAverage_vec = std::mem::take(&mut work.iterates.yAverage);
        let mut dualResidualAverage = std::mem::take(&mut work.resobj.dualResidualAverage);
        let mut dSlackPosAverage = std::mem::take(&mut work.resobj.dSlackPosAverage);
        let mut dSlackNegAverage = std::mem::take(&mut work.resobj.dSlackNegAverage);
        let mut dDualFeasAverage: cupdlp_float = work.resobj.dDualFeasAverage;
        let mut dDualObjAverage: cupdlp_float = work.resobj.dDualObjAverage;
        let mut dComplementarityAverage: cupdlp_float = work.resobj.dComplementarityAverage;
        PDHG_Compute_Dual_Feasibility(
            work,
            &mut dualResidualAverage,
            &atyAverage_vec.data,
            &xAverage_vec.data,
            &yAverage_vec.data,
            &mut dDualFeasAverage,
            &mut dDualObjAverage,
            &mut dComplementarityAverage,
            &mut dSlackPosAverage,
            &mut dSlackNegAverage,
        );
        work.iterates.atyAverage = atyAverage_vec;
        work.iterates.xAverage = xAverage_vec;
        work.iterates.yAverage = yAverage_vec;
        work.resobj.dualResidualAverage = dualResidualAverage;
        work.resobj.dSlackPosAverage = dSlackPosAverage;
        work.resobj.dSlackNegAverage = dSlackNegAverage;
        work.resobj.dDualFeasAverage = dDualFeasAverage;
        work.resobj.dDualObjAverage = dDualObjAverage;
        work.resobj.dComplementarityAverage = dComplementarityAverage;
    }

    // resobj->dPrimalObj /= (scaling->dObjScale * scaling->dObjScale);
    // resobj->dDualObj /= (scaling->dObjScale * scaling->dObjScale);
    // (commented out in C)
    work.resobj.dDualityGap = work.resobj.dPrimalObj - work.resobj.dDualObj;
    work.resobj.dRelObjGap = (work.resobj.dPrimalObj - work.resobj.dDualObj).abs()
        / (1.0 + work.resobj.dPrimalObj.abs() + work.resobj.dDualObj.abs());

    // resobj->dPrimalObjAverage /= scaling->dObjScale * scaling->dObjScale;
    // resobj->dDualObjAverage /= scaling->dObjScale * scaling->dObjScale;
    // (commented out in C)
    work.resobj.dDualityGapAverage = work.resobj.dPrimalObjAverage - work.resobj.dDualObjAverage;
    work.resobj.dRelObjGapAverage = (work.resobj.dPrimalObjAverage - work.resobj.dDualObjAverage)
        .abs()
        / (1.0 + work.resobj.dPrimalObjAverage.abs() + work.resobj.dDualObjAverage.abs());
}

/// HiGHS: cupdlp_solver.c::PDHG_Init_Variables (lines 531-591)
///
/// The `x_norm` `cupdlp_twoNorm` feeds only a printf and is not computed.
pub fn PDHG_Init_Variables(has_variables: &cupdlp_int, work: &mut CUPDLPwork) {
    let nCols = work.problem.data.nCols;
    let nRows = work.problem.data.nRows;

    let iter = work.timers.nIter;
    let idx = (iter % 2) as usize;

    // if (!*has_variables)
    //   CUPDLP_ZERO_VEC(x->data, cupdlp_float, lp->nCols);
    if *has_variables == 0 {
        cupdlp_zero(&mut work.iterates.x[idx].data, nCols);
    }

    // XXX: PDLP Does not project x0,  so we uncomment for 1-1 comparison
    // PDHG_Project_Bounds(work, x->data);
    PDHG_Project_Bounds(&work.problem, &mut work.iterates.x[idx].data);

    // double x_norm;
    // cupdlp_twoNorm(work, lp->nCols, x->data, &x_norm);
    // if (work->settings->nLogLevel>0)
    //   cupdlp_printf("||x0||_2 = %e\n", x_norm);
    // x_norm feeds only a printf; not computed.

    // cupdlp_zero(iterates->y, cupdlp_float, lp->nRows);
    // if (!*has_variables)
    //   CUPDLP_ZERO_VEC(y->data, cupdlp_float, lp->nRows);
    if *has_variables == 0 {
        cupdlp_zero(&mut work.iterates.y[idx].data, nRows);
    }

    // Ax(work, iterates->ax, iterates->x);
    // ATyCPU(work, iterates->aty, iterates->y);
    // Ax(work, ax, x);
    // ATy(work, aty, y);
    {
        let mut ax_vec = std::mem::take(&mut work.iterates.ax[idx]);
        let x_vec = std::mem::take(&mut work.iterates.x[idx]);
        Ax(&work.problem, &mut work.timers, &mut ax_vec, &x_vec);
        work.iterates.ax[idx] = ax_vec;
        work.iterates.x[idx] = x_vec;
    }
    {
        let mut aty_vec = std::mem::take(&mut work.iterates.aty[idx]);
        let y_vec = std::mem::take(&mut work.iterates.y[idx]);
        ATy(&work.problem, &mut work.timers, &mut aty_vec, &y_vec);
        work.iterates.aty[idx] = aty_vec;
        work.iterates.y[idx] = y_vec;
    }

    // #if PDLP_DEBUG_LOG block dropped (lines 563-572): PDLP_DEBUG_LOG == 0.

    // cupdlp_zero(iterates->xSum, cupdlp_float, lp->nCols);
    // cupdlp_zero(iterates->ySum, cupdlp_float, lp->nRows);
    // cupdlp_zero(iterates->xAverage, cupdlp_float, lp->nCols);
    // cupdlp_zero(iterates->yAverage, cupdlp_float, lp->nRows);
    // CUPDLP_ZERO_VEC(iterates->xSum, cupdlp_float, lp->nCols);
    // CUPDLP_ZERO_VEC(iterates->ySum, cupdlp_float, lp->nRows);
    // CUPDLP_ZERO_VEC(iterates->xAverage->data, cupdlp_float, lp->nCols);
    // CUPDLP_ZERO_VEC(iterates->yAverage->data, cupdlp_float, lp->nRows);
    cupdlp_zero(&mut work.iterates.xSum, nCols);
    cupdlp_zero(&mut work.iterates.ySum, nRows);
    cupdlp_zero(&mut work.iterates.xAverage.data, nCols);
    cupdlp_zero(&mut work.iterates.yAverage.data, nRows);

    // PDHG_Project_Bounds(work, iterates->xSum);
    // PDHG_Project_Bounds(work, iterates->xAverage->data);
    PDHG_Project_Bounds(&work.problem, &mut work.iterates.xSum);
    PDHG_Project_Bounds(&work.problem, &mut work.iterates.xAverage.data);

    work.stepsize.dSumPrimalStep = 0.0;
    work.stepsize.dSumDualStep = 0.0;

    // CUPDLP_ZERO_VEC(iterates->xLastRestart, cupdlp_float, lp->nCols);
    // CUPDLP_ZERO_VEC(iterates->yLastRestart, cupdlp_float, lp->nRows);
    cupdlp_zero(&mut work.iterates.xLastRestart, nCols);
    cupdlp_zero(&mut work.iterates.yLastRestart, nRows);
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn init_variables_cold_start_zeroes_and_projects() {
        let mut w = CUPDLPwork::default();
        w.problem.nRows = 1;
        w.problem.nCols = 1;
        w.problem.lower = vec![1.0];
        w.problem.upper = vec![2.0];
        w.problem.hasLower = vec![1.0];
        w.problem.hasUpper = vec![1.0];
        let v1 = |v: f64| CUPDLPvec {
            len: 1,
            data: vec![v],
        };
        w.iterates.x = [v1(5.0), v1(5.0)];
        w.iterates.y = [v1(5.0), v1(5.0)];
        w.iterates.ax = [v1(0.0), v1(0.0)];
        w.iterates.aty = [v1(0.0), v1(0.0)];
        w.iterates.xAverage = v1(5.0);
        w.iterates.yAverage = v1(5.0);
        w.iterates.xSum = vec![5.0];
        w.iterates.ySum = vec![5.0];
        w.iterates.xLastRestart = vec![5.0];
        w.iterates.yLastRestart = vec![5.0];
        w.problem.data.nRows = 1;
        w.problem.data.nCols = 1;
        w.problem.data.matrix_format = CSR_CSC;
        w.problem.data.csc_matrix = Some(Box::new(CUPDLPcsc {
            nRows: 1,
            nCols: 1,
            nMatElem: 1,
            colMatBeg: vec![0, 1],
            colMatIdx: vec![0],
            colMatElem: vec![1.0],
            MatElemNormInf: 1.0,
        }));
        w.problem.data.csr_matrix = Some(Box::new(CUPDLPcsr {
            nRows: 1,
            nCols: 1,
            nMatElem: 1,
            rowMatBeg: vec![0, 1],
            rowMatIdx: vec![0],
            rowMatElem: vec![1.0],
        }));
        PDHG_Init_Variables(&0, &mut w);
        // zero then project onto [1, 2]
        assert_eq!(w.iterates.x[0].data, vec![1.0]);
        assert_eq!(w.iterates.y[0].data, vec![0.0]);
        assert_eq!(w.iterates.xSum, vec![1.0]);
        assert_eq!(w.iterates.xLastRestart, vec![0.0]);
        assert_eq!(w.stepsize.dSumPrimalStep, 0.0);
    }
}
