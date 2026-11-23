//! HiGHS source: pdlp/cupdlp/cupdlp_step.h, pdlp/cupdlp/cupdlp_step.c
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! PDHG iterate updates and step sizes. Borrowing: `x[iter % 2]` and
//! `x[(iter + 1) % 2]` are disjoint elements of one `[CUPDLPvec; 2]`; the
//! bodies take them with `let [x0, x1] = &mut iterates.x;` and select by
//! parity.
//!
//! `PDHG_Init_Step_Sizes` calls `PDHG_Power_Method(pdhg,
//! &stepsize->dPrimalStep)`: the lambda out-parameter aliases a work field
//! the method never reads, so the port passes a local and stores it back.

use super::defs::*;
use super::linalg::*;
use super::solver::check::PDHG_Compute_SolvingTime;
use super::utils::*;

/// HiGHS: cupdlp_step.c::PDHG_primalGradientStep (lines 16-40)
///
/// `work` replaced by `problem`; CPU branch.
pub fn PDHG_primalGradientStep(
    problem: &CUPDLPproblem,
    xUpdate: &mut CUPDLPvec,
    x: &CUPDLPvec,
    ATy: &CUPDLPvec,
    dPrimalStepSize: cupdlp_float,
) {
    // CUPDLP_COPY_VEC(xUpdate->data, x->data, cupdlp_float, problem->nCols);
    xUpdate.data[..problem.nCols as usize].copy_from_slice(&x.data[..problem.nCols as usize]);

    // cupdlp_float alpha = -dPrimalStepSize;
    let alpha = -dPrimalStepSize;
    // cupdlp_axpy(work, problem->nCols, &alpha, problem->cost, xUpdate->data);
    cupdlp_axpy(problem.nCols, alpha, &problem.cost, &mut xUpdate.data);

    // alpha = dPrimalStepSize;
    let alpha = dPrimalStepSize;
    // cupdlp_axpy(work, problem->nCols, &alpha, ATy->data, xUpdate->data);
    cupdlp_axpy(problem.nCols, alpha, &ATy.data, &mut xUpdate.data);

    // cupdlp_projub(xUpdate->data, problem->upper, problem->nCols);
    cupdlp_projub(&mut xUpdate.data, &problem.upper, problem.nCols);
    // cupdlp_projlb(xUpdate->data, problem->lower, problem->nCols);
    cupdlp_projlb(&mut xUpdate.data, &problem.lower, problem.nCols);
}

/// HiGHS: cupdlp_step.c::PDHG_dualGradientStep (lines 43-69)
///
/// `work` replaced by `problem`; CPU branch. Only the
/// inequality rows `[nEqs, nRows)` are projected onto `y >= 0`.
pub fn PDHG_dualGradientStep(
    problem: &CUPDLPproblem,
    yUpdate: &mut CUPDLPvec,
    y: &CUPDLPvec,
    Ax: &CUPDLPvec,
    AxUpdate: &CUPDLPvec,
    dDualStepSize: cupdlp_float,
) {
    // CUPDLP_COPY_VEC(yUpdate->data, y->data, cupdlp_float, problem->nRows);
    yUpdate.data[..problem.nRows as usize].copy_from_slice(&y.data[..problem.nRows as usize]);

    // cupdlp_float alpha = dDualStepSize;
    let alpha = dDualStepSize;
    // cupdlp_axpy(work, problem->nRows, &alpha, problem->rhs, yUpdate->data);
    cupdlp_axpy(problem.nRows, alpha, &problem.rhs, &mut yUpdate.data);

    // alpha = -2.0 * dDualStepSize;
    let alpha = -2.0 * dDualStepSize;
    // cupdlp_axpy(work, problem->nRows, &alpha, AxUpdate->data, yUpdate->data);
    cupdlp_axpy(problem.nRows, alpha, &AxUpdate.data, &mut yUpdate.data);

    // alpha = dDualStepSize;
    let alpha = dDualStepSize;
    // cupdlp_axpy(work, problem->nRows, &alpha, Ax->data, yUpdate->data);
    cupdlp_axpy(problem.nRows, alpha, &Ax.data, &mut yUpdate.data);

    // cupdlp_projPos(yUpdate->data + problem->nEqs, problem->nRows - problem->nEqs);
    let n_eqs = problem.nEqs as usize;
    let n_rows = problem.nRows as usize;
    cupdlp_projPos(
        &mut yUpdate.data[n_eqs..n_rows],
        problem.nRows - problem.nEqs,
    );
}

/// HiGHS: cupdlp_step.c::PDHG_Power_Method (lines 71-145)
///
/// The residual `res` feeds only the `log_iters` printf, and upstream
/// computes it as `cupdlp_twoNormSquared(work, lp->nCols, ax->data, &res)`
/// over an `nRows` vector (out-of-bounds read when `nCols > nRows`): it is
/// not computed. `dl_lambda` / `previous_lambda` feed only the printf too.
pub fn PDHG_Power_Method(work: &mut CUPDLPwork, lambda: &mut f64) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    let problem = &work.problem;
    let lp = &problem.data;
    let iterates = &work.iterates;

    // cupdlp_float *q = work->buffer->data;
    // cupdlp_initvec(q, 1.0, lp->nRows);
    cupdlp_initvec(&mut work.buffer.data, 1.0, lp.nRows);

    // cupdlp_int iter = work->timers->nIter;
    let iter = work.timers.nIter;
    // CUPDLPvec *ax = iterates->ax[iter % 2];
    // CUPDLPvec *aty = iterates->aty[iter % 2];
    let ax_idx = (iter % 2) as usize;
    let aty_idx = (iter % 2) as usize;

    // double res = 0.0;
    // double previous_lambda = 0.0;
    // int log_iters = work->settings->nLogLevel > 1;

    // for (cupdlp_int iter = 0; iter < 20; ++iter) {
    for _iter in 0..20 {
        // ATy(work, aty, work->buffer);
        {
            let aty = &mut work.iterates.aty[aty_idx];
            let buffer = &work.buffer;
            ATy(&work.problem, &mut work.timers, aty, buffer);
        }
        // Ax(work, ax, aty);
        {
            let ax = &mut work.iterates.ax[ax_idx];
            let aty = &work.iterates.aty[aty_idx];
            Ax(&work.problem, &mut work.timers, ax, aty);
        }

        // CUPDLP_COPY_VEC(q, ax->data, cupdlp_float, lp->nRows);
        {
            let ax_data = &work.iterates.ax[ax_idx].data;
            let n_rows = lp.nRows as usize;
            work.buffer.data[..n_rows].copy_from_slice(&ax_data[..n_rows]);
        }

        // cupdlp_float qNorm = 0.0;
        let mut q_norm: cupdlp_float = 0.0;
        // cupdlp_twoNorm(work, lp->nRows, q, &qNorm);
        {
            let q = &work.buffer.data;
            cupdlp_twoNorm(lp.nRows, q, &mut q_norm);
        }
        // cupdlp_scaleVector(work, 1.0 / qNorm, q, lp->nRows);
        {
            let q = &mut work.buffer.data;
            cupdlp_scaleVector(1.0 / q_norm, q, lp.nRows);
        }

        // ATy(work, aty, work->buffer);
        {
            let aty = &mut work.iterates.aty[aty_idx];
            let buffer = &work.buffer;
            ATy(&work.problem, &mut work.timers, aty, buffer);
        }

        // cupdlp_twoNormSquared(work, lp->nCols, aty->data, lambda);
        {
            let aty_data = &work.iterates.aty[aty_idx].data;
            cupdlp_twoNormSquared(lp.nCols, aty_data, lambda);
        }

        // cupdlp_float alpha = -(*lambda);
        let alpha = -(*lambda);
        // cupdlp_axpy(work, lp->nRows, &alpha, q, ax->data);
        {
            let q = &work.buffer.data;
            let ax_data = &mut work.iterates.ax[ax_idx].data;
            cupdlp_axpy(lp.nRows, alpha, q, ax_data);
        }

        // cupdlp_twoNormSquared(work, lp->nCols, ax->data, &res);
        // (res feeds only printf; not computed)

        // double dl_lambda = fabs(*lambda - previous_lambda);
        // previous_lambda = *lambda;
        // (dl_lambda / previous_lambda feed only printf)
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_step.c::PDHG_Compute_Step_Size_Ratio (lines 147-176)
pub fn PDHG_Compute_Step_Size_Ratio(pdhg: &mut CUPDLPwork) {
    let problem = &pdhg.problem;
    let iterates = &pdhg.iterates;
    let stepsize = &pdhg.stepsize;

    // cupdlp_float dMeanStepSize = sqrt(stepsize->dPrimalStep * stepsize->dDualStep);
    let d_mean_step_size = (stepsize.dPrimalStep * stepsize.dDualStep).sqrt();

    // cupdlp_int iter = pdhg->timers->nIter;
    let iter = pdhg.timers.nIter;
    // CUPDLPvec *x = iterates->x[iter % 2];
    // CUPDLPvec *y = iterates->y[iter % 2];
    let x_idx = (iter % 2) as usize;
    let y_idx = (iter % 2) as usize;

    // cupdlp_float dDiffPrimal = 0.0;
    let mut d_diff_primal: cupdlp_float = 0.0;
    // cupdlp_diffTwoNorm(pdhg, x->data, iterates->xLastRestart, problem->nCols, &dDiffPrimal);
    {
        let x_data = &pdhg.iterates.x[x_idx].data;
        let x_last_restart = &pdhg.iterates.xLastRestart;
        let buffer2 = &mut pdhg.buffer2;
        cupdlp_diffTwoNorm(
            buffer2,
            x_data,
            x_last_restart,
            problem.nCols,
            &mut d_diff_primal,
        );
    }

    // cupdlp_float dDiffDual = 0.0;
    let mut d_diff_dual: cupdlp_float = 0.0;
    // cupdlp_diffTwoNorm(pdhg, y->data, iterates->yLastRestart, problem->nRows, &dDiffDual);
    {
        let y_data = &pdhg.iterates.y[y_idx].data;
        let y_last_restart = &pdhg.iterates.yLastRestart;
        let buffer2 = &mut pdhg.buffer2;
        cupdlp_diffTwoNorm(
            buffer2,
            y_data,
            y_last_restart,
            problem.nRows,
            &mut d_diff_dual,
        );
    }

    // if (fmin(dDiffPrimal, dDiffDual) > 1e-10) {
    // PORT NOTE: C fmin
    if d_diff_primal.min(d_diff_dual) > 1e-10 {
        // cupdlp_float dBetaUpdate = dDiffDual / dDiffPrimal;
        let d_beta_update = d_diff_dual / d_diff_primal;
        // cupdlp_float dLogBetaUpdate = 0.5 * log(dBetaUpdate) + 0.5 * log(sqrt(stepsize->dBeta));
        let d_log_beta_update = 0.5 * d_beta_update.ln() + 0.5 * (stepsize.dBeta.sqrt()).ln();
        // stepsize->dBeta = exp(dLogBetaUpdate) * exp(dLogBetaUpdate);
        pdhg.stepsize.dBeta = d_log_beta_update.exp() * d_log_beta_update.exp();
    }

    // stepsize->dPrimalStep = dMeanStepSize / sqrt(stepsize->dBeta);
    pdhg.stepsize.dPrimalStep = d_mean_step_size / pdhg.stepsize.dBeta.sqrt();
    // stepsize->dDualStep = stepsize->dPrimalStep * stepsize->dBeta;
    pdhg.stepsize.dDualStep = pdhg.stepsize.dPrimalStep * pdhg.stepsize.dBeta;
    // stepsize->dTheta = 1.0;
    pdhg.stepsize.dTheta = 1.0;
}

/// HiGHS: cupdlp_step.c::PDHG_Update_Iterate_Constant_Step_Size (lines 178-206)
pub fn PDHG_Update_Iterate_Constant_Step_Size(pdhg: &mut CUPDLPwork) {
    // CUPDLP_ASSERT(0); dropped (commented out in C)

    let problem = &pdhg.problem;
    let iterates = &pdhg.iterates;
    let stepsize = &pdhg.stepsize;

    // cupdlp_int iter = pdhg->timers->nIter;
    let iter = pdhg.timers.nIter;
    let idx = (iter % 2) as usize;
    let idx_update = ((iter + 1) % 2) as usize;

    // Ax(pdhg, ax, x);
    {
        let ax = &mut pdhg.iterates.ax[idx];
        let x = &pdhg.iterates.x[idx];
        Ax(&pdhg.problem, &mut pdhg.timers, ax, x);
    }
    // ATy(pdhg, aty, y);
    {
        let aty = &mut pdhg.iterates.aty[idx];
        let y = &pdhg.iterates.y[idx];
        ATy(&pdhg.problem, &mut pdhg.timers, aty, y);
    }

    // PDHG_primalGradientStep(pdhg, xUpdate, x, aty, stepsize->dPrimalStep);
    {
        let [x0, x1] = &mut pdhg.iterates.x;
        let (x, x_update): (&CUPDLPvec, &mut CUPDLPvec) =
            if idx == 0 { (&*x0, x1) } else { (&*x1, x0) };
        let aty = &pdhg.iterates.aty[idx];
        let d_primal_step = stepsize.dPrimalStep;
        PDHG_primalGradientStep(&pdhg.problem, x_update, x, aty, d_primal_step);
    }

    // Ax(pdhg, axUpdate, xUpdate);
    {
        let ax_update = &mut pdhg.iterates.ax[idx_update];
        let x_update = &pdhg.iterates.x[idx_update];
        Ax(&pdhg.problem, &mut pdhg.timers, ax_update, x_update);
    }

    // PDHG_dualGradientStep(pdhg, yUpdate, y, ax, axUpdate, stepsize->dDualStep);
    {
        let [y0, y1] = &mut pdhg.iterates.y;
        let (y, y_update): (&CUPDLPvec, &mut CUPDLPvec) =
            if idx == 0 { (&*y0, y1) } else { (&*y1, y0) };
        let ax = &pdhg.iterates.ax[idx];
        let ax_update = &pdhg.iterates.ax[idx_update];
        let d_dual_step = stepsize.dDualStep;
        PDHG_dualGradientStep(&pdhg.problem, y_update, y, ax, ax_update, d_dual_step);
    }

    // ATy(pdhg, atyUpdate, yUpdate);
    {
        let aty_update = &mut pdhg.iterates.aty[idx_update];
        let y_update = &pdhg.iterates.y[idx_update];
        ATy(&pdhg.problem, &mut pdhg.timers, aty_update, y_update);
    }
}

/// HiGHS: cupdlp_step.c::PDHG_Update_Iterate_Malitsky_Pock (lines 208-213)
///
/// Upstream prints and calls `exit(-1)`: the port panics with the same
/// message. Unreachable from HiGHS (`getUserParamsFromOptions` sets only
/// `PDHG_FIXED_LINESEARCH` or `PDHG_ADAPTIVE_LINESEARCH`).
pub fn PDHG_Update_Iterate_Malitsky_Pock(_pdhg: &mut CUPDLPwork) {
    // cupdlp_printf("Malitsky-Pock is not implemented\n");
    // cupdlp_printf(" - use %d and %d instead", PDHG_FIXED_LINESEARCH, PDHG_ADAPTIVE_LINESEARCH);
    // exit(-1);
    panic!(
        "Malitsky-Pock is not implemented\n - use {} and {} instead",
        PDHG_FIXED_LINESEARCH, PDHG_ADAPTIVE_LINESEARCH
    );
}

/// HiGHS: cupdlp_step.c::PDHG_Update_Iterate_Adaptive_Step_Size (lines 215-310)
pub fn PDHG_Update_Iterate_Adaptive_Step_Size(pdhg: &mut CUPDLPwork) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;
    let problem = &pdhg.problem;
    let iterates = &pdhg.iterates;
    let stepsize = &pdhg.stepsize;

    // cupdlp_int iter = pdhg->timers->nIter;
    let iter = pdhg.timers.nIter;
    let idx = (iter % 2) as usize;
    let idx_update = ((iter + 1) % 2) as usize;

    // cupdlp_float dStepSizeUpdate = sqrt(stepsize->dPrimalStep * stepsize->dDualStep);
    let mut d_step_size_update = (stepsize.dPrimalStep * stepsize.dDualStep).sqrt();

    // cupdlp_bool isDone = false;
    let mut is_done = false;
    // int stepIterThis = 0;
    let mut _step_iter_this = 0;

    'exit_cleanup: {
        // while (!isDone) {
        while !is_done {
            // ++stepsize->nStepSizeIter;
            pdhg.stepsize.nStepSizeIter += 1;
            // ++stepIterThis;
            _step_iter_this += 1;

            // cupdlp_float dPrimalStepUpdate = dStepSizeUpdate / sqrt(stepsize->dBeta);
            let d_primal_step_update = d_step_size_update / pdhg.stepsize.dBeta.sqrt();
            // cupdlp_float dDualStepUpdate = dStepSizeUpdate * sqrt(stepsize->dBeta);
            let d_dual_step_update = d_step_size_update * pdhg.stepsize.dBeta.sqrt();

            // PDHG_primalGradientStep(pdhg, xUpdate, x, aty, dPrimalStepUpdate);
            {
                let [x0, x1] = &mut pdhg.iterates.x;
                let (x, x_update): (&CUPDLPvec, &mut CUPDLPvec) =
                    if idx == 0 { (&*x0, x1) } else { (&*x1, x0) };
                let aty = &pdhg.iterates.aty[idx];
                PDHG_primalGradientStep(&pdhg.problem, x_update, x, aty, d_primal_step_update);
            }

            // Ax(pdhg, axUpdate, xUpdate);
            {
                let ax_update = &mut pdhg.iterates.ax[idx_update];
                let x_update = &pdhg.iterates.x[idx_update];
                Ax(&pdhg.problem, &mut pdhg.timers, ax_update, x_update);
            }

            // PDHG_dualGradientStep(pdhg, yUpdate, y, ax, axUpdate, dDualStepUpdate);
            {
                let [y0, y1] = &mut pdhg.iterates.y;
                let (y, y_update): (&CUPDLPvec, &mut CUPDLPvec) =
                    if idx == 0 { (&*y0, y1) } else { (&*y1, y0) };
                let ax = &pdhg.iterates.ax[idx];
                let ax_update = &pdhg.iterates.ax[idx_update];
                PDHG_dualGradientStep(
                    &pdhg.problem,
                    y_update,
                    y,
                    ax,
                    ax_update,
                    d_dual_step_update,
                );
            }

            // ATy(pdhg, atyUpdate, yUpdate);
            {
                let aty_update = &mut pdhg.iterates.aty[idx_update];
                let y_update = &pdhg.iterates.y[idx_update];
                ATy(&pdhg.problem, &mut pdhg.timers, aty_update, y_update);
            }

            // cupdlp_float dMovement = 0.0;
            let mut d_movement: cupdlp_float = 0.0;
            // cupdlp_float dInteraction = 0.0;
            let mut d_interaction: cupdlp_float = 0.0;

            // cupdlp_compute_interaction_and_movement(pdhg, &dMovement, &dInteraction);
            cupdlp_compute_interaction_and_movement(pdhg, &mut d_movement, &mut d_interaction);

            // #if CUPDLP_DUMP_LINESEARCH_STATS && CUPDLP_DEBUG ... dropped

            // cupdlp_float dStepSizeLimit;
            // if (dInteraction != 0.0) {
            let d_step_size_limit = if d_interaction != 0.0 {
                // dStepSizeLimit = dMovement / fabs(dInteraction);
                d_movement / d_interaction.abs()
            } else {
                // dStepSizeLimit = INFINITY;
                f64::INFINITY
            };

            // if (dStepSizeUpdate <= dStepSizeLimit) {
            if d_step_size_update <= d_step_size_limit {
                // isDone = true;
                is_done = true;
            } else {
                // CUPDLP_CHECK_TIMEOUT(pdhg);
                // {
                //   PDHG_Compute_SolvingTime(pdhg);
                //   if (pdhg->timers->dSolvingTime > pdhg->settings->dTimeLim) {
                //     retcode = RETCODE_FAILED;
                //     goto exit_cleanup;
                //   }
                // }
                PDHG_Compute_SolvingTime(pdhg);
                if pdhg.timers.dSolvingTime > pdhg.settings.dTimeLim {
                    retcode = RETCODE_FAILED;
                    // goto exit_cleanup;
                    break 'exit_cleanup;
                }
            }

            // cupdlp_float dFirstTerm = (1.0 - pow(stepsize->nStepSizeIter + 1.0, -PDHG_STEPSIZE_REDUCTION_EXP)) * dStepSizeLimit;
            let d_first_term = (1.0
                - ((pdhg.stepsize.nStepSizeIter as f64) + 1.0).powf(-PDHG_STEPSIZE_REDUCTION_EXP))
                * d_step_size_limit;

            // cupdlp_float dSecondTerm = (1.0 + pow(stepsize->nStepSizeIter + 1.0, -PDHG_STEPSIZE_GROWTH_EXP)) * dStepSizeUpdate;
            let d_second_term = (1.0
                + ((pdhg.stepsize.nStepSizeIter as f64) + 1.0).powf(-PDHG_STEPSIZE_GROWTH_EXP))
                * d_step_size_update;

            // dStepSizeUpdate = fmin(dFirstTerm, dSecondTerm);
            // PORT NOTE: C fmin
            d_step_size_update = d_first_term.min(d_second_term);

            // #if CUPDLP_DUMP_LINESEARCH_STATS && CUPDLP_DEBUG ... dropped
        }

        // stepsize->dPrimalStep = dStepSizeUpdate / sqrt(stepsize->dBeta);
        pdhg.stepsize.dPrimalStep = d_step_size_update / pdhg.stepsize.dBeta.sqrt();
        // stepsize->dDualStep = dStepSizeUpdate * sqrt(stepsize->dBeta);
        pdhg.stepsize.dDualStep = d_step_size_update * pdhg.stepsize.dBeta.sqrt();
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_step.c::PDHG_Init_Step_Sizes (lines 312-375)
pub fn PDHG_Init_Step_Sizes(pdhg: &mut CUPDLPwork) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;

    'exit_cleanup: {
        // if (stepsize->eLineSearchMethod == PDHG_FIXED_LINESEARCH) {
        if pdhg.stepsize.eLineSearchMethod == PDHG_FIXED_LINESEARCH {
            // CUPDLP_CALL(PDHG_Power_Method(pdhg, &stepsize->dPrimalStep));
            // The lambda out-param aliases stepsize->dPrimalStep which Power_Method never reads.
            // Pass a local and store back.
            let mut power_method_lambda: f64 = 0.0;
            retcode = PDHG_Power_Method(pdhg, &mut power_method_lambda);
            if retcode != RETCODE_OK {
                break 'exit_cleanup;
            }
            pdhg.stepsize.dPrimalStep = power_method_lambda;

            // cupdlp_float power_method_lambda = stepsize->dPrimalStep;
            // (already stored above)

            // cupdlp_float a = 0.0;
            let mut a: cupdlp_float = 0.0;
            // cupdlp_float b = 0.0;
            let mut b: cupdlp_float = 0.0;
            // cupdlp_twoNormSquared(pdhg, problem->nCols, problem->cost, &a);
            cupdlp_twoNormSquared(pdhg.problem.nCols, &pdhg.problem.cost, &mut a);
            // cupdlp_twoNormSquared(pdhg, problem->nRows, problem->rhs, &b);
            cupdlp_twoNormSquared(pdhg.problem.nRows, &pdhg.problem.rhs, &mut b);

            // if (fmin(a, b) > 1e-6) {
            // PORT NOTE: C fmin
            if a.min(b) > 1e-6 {
                // stepsize->dBeta = a / b;
                pdhg.stepsize.dBeta = a / b;
            } else {
                // stepsize->dBeta = 1.0;
                pdhg.stepsize.dBeta = 1.0;
            }

            // stepsize->dPrimalStep = 0.8 / sqrt(stepsize->dPrimalStep);
            pdhg.stepsize.dPrimalStep = 0.8 / pdhg.stepsize.dPrimalStep.sqrt();
            // stepsize->dDualStep = stepsize->dPrimalStep;
            pdhg.stepsize.dDualStep = pdhg.stepsize.dPrimalStep;
            // stepsize->dPrimalStep /= sqrt(stepsize->dBeta);
            pdhg.stepsize.dPrimalStep /= pdhg.stepsize.dBeta.sqrt();
            // stepsize->dDualStep *= sqrt(stepsize->dBeta);
            pdhg.stepsize.dDualStep *= pdhg.stepsize.dBeta.sqrt();

            // if (pdhg->settings->nLogLevel > 1) cupdlp_printf(...); dropped
        } else {
            // stepsize->dTheta = 1.0;
            pdhg.stepsize.dTheta = 1.0;

            // cupdlp_float a = 0.0;
            let mut a: cupdlp_float = 0.0;
            // cupdlp_float b = 0.0;
            let mut b: cupdlp_float = 0.0;
            // cupdlp_twoNormSquared(pdhg, problem->nCols, problem->cost, &a);
            cupdlp_twoNormSquared(pdhg.problem.nCols, &pdhg.problem.cost, &mut a);
            // cupdlp_twoNormSquared(pdhg, problem->nRows, problem->rhs, &b);
            cupdlp_twoNormSquared(pdhg.problem.nRows, &pdhg.problem.rhs, &mut b);

            // if (fmin(a, b) > 1e-6) {
            // PORT NOTE: C fmin
            if a.min(b) > 1e-6 {
                // stepsize->dBeta = a / b;
                pdhg.stepsize.dBeta = a / b;
            } else {
                // stepsize->dBeta = 1.0;
                pdhg.stepsize.dBeta = 1.0;
            }

            // stepsize->dPrimalStep = (1.0 / problem->data->csc_matrix->MatElemNormInf) / sqrt(stepsize->dBeta);
            pdhg.stepsize.dPrimalStep = (1.0
                / pdhg
                    .problem
                    .data
                    .csc_matrix
                    .as_ref()
                    .unwrap()
                    .MatElemNormInf)
                / pdhg.stepsize.dBeta.sqrt();
            // stepsize->dDualStep = stepsize->dPrimalStep * stepsize->dBeta;
            pdhg.stepsize.dDualStep = pdhg.stepsize.dPrimalStep * pdhg.stepsize.dBeta;
            // iterates->dLastRestartBeta = stepsize->dBeta;
            pdhg.iterates.dLastRestartBeta = pdhg.stepsize.dBeta;
        }

        // iterates->iLastRestartIter = 0;
        pdhg.iterates.iLastRestartIter = 0;
        // stepsize->dSumPrimalStep = 0;
        pdhg.stepsize.dSumPrimalStep = 0.0;
        // stepsize->dSumDualStep = 0;
        pdhg.stepsize.dSumDualStep = 0.0;
    }

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_step.c::PDHG_Compute_Average_Iterate (lines 377-420)
pub fn PDHG_Compute_Average_Iterate(work: &mut CUPDLPwork) {
    let problem = &work.problem;
    let lp = &problem.data;
    let stepsize = &work.stepsize;
    let iterates = &work.iterates;

    // cupdlp_float dPrimalScale = stepsize->dSumPrimalStep > 0.0 ? 1.0 / stepsize->dSumPrimalStep : 1.0;
    let d_primal_scale: cupdlp_float = if stepsize.dSumPrimalStep > 0.0 {
        1.0 / stepsize.dSumPrimalStep
    } else {
        1.0
    };
    // cupdlp_float dDualScale = stepsize->dSumDualStep > 0.0 ? 1.0 / stepsize->dSumDualStep : 1.0;
    let d_dual_scale: cupdlp_float = if stepsize.dSumDualStep > 0.0 {
        1.0 / stepsize.dSumDualStep
    } else {
        1.0
    };

    // CUPDLP_COPY_VEC(iterates->xAverage->data, iterates->xSum, cupdlp_float, lp->nCols);
    {
        let n_cols = lp.nCols as usize;
        let x_sum = &work.iterates.xSum;
        let x_average = &mut work.iterates.xAverage;
        x_average.data[..n_cols].copy_from_slice(&x_sum[..n_cols]);
    }
    // CUPDLP_COPY_VEC(iterates->yAverage->data, iterates->ySum, cupdlp_float, lp->nRows);
    {
        let n_rows = lp.nRows as usize;
        let y_sum = &work.iterates.ySum;
        let y_average = &mut work.iterates.yAverage;
        y_average.data[..n_rows].copy_from_slice(&y_sum[..n_rows]);
    }

    // cupdlp_scaleVector(work, dPrimalScale, iterates->xAverage->data, lp->nCols);
    {
        let x_average = &mut work.iterates.xAverage;
        cupdlp_scaleVector(d_primal_scale, &mut x_average.data, lp.nCols);
    }
    // cupdlp_scaleVector(work, dDualScale, iterates->yAverage->data, lp->nRows);
    {
        let y_average = &mut work.iterates.yAverage;
        cupdlp_scaleVector(d_dual_scale, &mut y_average.data, lp.nRows);
    }

    // Ax(work, iterates->axAverage, iterates->xAverage);
    {
        let ax_average = &mut work.iterates.axAverage;
        let x_average = &work.iterates.xAverage;
        Ax(&work.problem, &mut work.timers, ax_average, x_average);
    }
    // ATy(work, iterates->atyAverage, iterates->yAverage);
    {
        let aty_average = &mut work.iterates.atyAverage;
        let y_average = &work.iterates.yAverage;
        ATy(&work.problem, &mut work.timers, aty_average, y_average);
    }

    // #if PDLP_DEBUG_LOG ... dropped (lines 404-419)
}

/// HiGHS: cupdlp_step.c::PDHG_Update_Average (lines 422-442)
pub fn PDHG_Update_Average(work: &mut CUPDLPwork) {
    let problem = &work.problem;
    let lp = &problem.data;
    let stepsize = &work.stepsize;
    let iterates = &work.iterates;

    // cupdlp_int iter = work->timers->nIter;
    let iter = work.timers.nIter;
    let idx_update = ((iter + 1) % 2) as usize;

    // cupdlp_float dMeanStepSize = sqrt(stepsize->dPrimalStep * stepsize->dDualStep);
    let d_mean_step_size = (stepsize.dPrimalStep * stepsize.dDualStep).sqrt();

    // cupdlp_axpy(work, lp->nCols, &dMeanStepSize, xUpdate->data, iterates->xSum);
    {
        let x_update_data = &work.iterates.x[idx_update].data;
        let x_sum = &mut work.iterates.xSum;
        cupdlp_axpy(lp.nCols, d_mean_step_size, x_update_data, x_sum);
    }
    // cupdlp_axpy(work, lp->nRows, &dMeanStepSize, yUpdate->data, iterates->ySum);
    {
        let y_update_data = &work.iterates.y[idx_update].data;
        let y_sum = &mut work.iterates.ySum;
        cupdlp_axpy(lp.nRows, d_mean_step_size, y_update_data, y_sum);
    }

    // stepsize->dSumPrimalStep += dMeanStepSize;
    work.stepsize.dSumPrimalStep += d_mean_step_size;
    // stepsize->dSumDualStep += dMeanStepSize;
    work.stepsize.dSumDualStep += d_mean_step_size;
}

/// HiGHS: cupdlp_step.c::PDHG_Update_Iterate (lines 444-478)
///
/// Returns `RETCODE_OK` even when the adaptive step fails (upstream
/// `exit_cleanup: return RETCODE_OK;`); that path skips
/// `PDHG_Update_Average` and the `dUpdateIterateTime` update, as in C.
pub fn PDHG_Update_Iterate(pdhg: &mut CUPDLPwork) -> cupdlp_retcode {
    // #if PDHG_USE_TIMERS (defined to 1 in cupdlp_defs.h)
    // ++timers->nUpdateIterateCalls;
    pdhg.timers.nUpdateIterateCalls += 1;
    // cupdlp_float dStartTime = getTimeStamp();
    let d_start_time: cupdlp_float = getTimeStamp();
    // #endif

    'exit_cleanup: {
        // switch (stepsize->eLineSearchMethod) {
        match pdhg.stepsize.eLineSearchMethod {
            // case PDHG_FIXED_LINESEARCH:
            PDHG_FIXED_LINESEARCH => {
                // PDHG_Update_Iterate_Constant_Step_Size(pdhg);
                PDHG_Update_Iterate_Constant_Step_Size(pdhg);
            }
            // case PDHG_MALITSKY_POCK_LINESEARCH:
            PDHG_MALITSKY_POCK_LINESEARCH => {
                // PDHG_Update_Iterate_Malitsky_Pock(pdhg);
                PDHG_Update_Iterate_Malitsky_Pock(pdhg);
            }
            // case PDHG_ADAPTIVE_LINESEARCH:
            PDHG_ADAPTIVE_LINESEARCH => {
                // CUPDLP_CALL(PDHG_Update_Iterate_Adaptive_Step_Size(pdhg));
                let retcode = PDHG_Update_Iterate_Adaptive_Step_Size(pdhg);
                if retcode != RETCODE_OK {
                    // goto exit_cleanup;
                    break 'exit_cleanup;
                }
            }
            _ => {}
        }

        // PDHG_Update_Average(pdhg);
        PDHG_Update_Average(pdhg);

        // #if PDHG_USE_TIMERS
        // timers->dUpdateIterateTime += getTimeStamp() - dStartTime;
        pdhg.timers.dUpdateIterateTime += getTimeStamp() - d_start_time;
        // #endif
    }

    // exit_cleanup:
    RETCODE_OK
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn primal_gradient_step_projects() {
        let mut p = CUPDLPproblem::default();
        p.nCols = 2;
        p.cost = vec![1.0, -1.0];
        p.lower = vec![0.0, 0.0];
        p.upper = vec![10.0, 1.0];
        let x = CUPDLPvec {
            len: 2,
            data: vec![0.5, 0.5],
        };
        let aty = CUPDLPvec {
            len: 2,
            data: vec![0.0, 0.0],
        };
        let mut xu = CUPDLPvec {
            len: 2,
            data: vec![0.0; 2],
        };
        // x - 1 * (c - A'y) = [-0.5, 1.5] -> projected [0, 1]
        PDHG_primalGradientStep(&p, &mut xu, &x, &aty, 1.0);
        assert_eq!(xu.data, vec![0.0, 1.0]);
    }

    #[test]
    fn dual_gradient_step_projects_inequalities_only() {
        let mut p = CUPDLPproblem::default();
        p.nRows = 2;
        p.nEqs = 1;
        p.rhs = vec![-1.0, -1.0];
        let y = CUPDLPvec {
            len: 2,
            data: vec![0.0, 0.0],
        };
        let ax = CUPDLPvec {
            len: 2,
            data: vec![0.0, 0.0],
        };
        let axu = CUPDLPvec {
            len: 2,
            data: vec![0.0, 0.0],
        };
        let mut yu = CUPDLPvec {
            len: 2,
            data: vec![0.0; 2],
        };
        PDHG_dualGradientStep(&p, &mut yu, &y, &ax, &axu, 1.0);
        assert_eq!(yu.data, vec![-1.0, 0.0]);
    }

    #[test]
    fn step_size_ratio_keeps_mean_step() {
        let mut w = CUPDLPwork::default();
        w.problem.nCols = 1;
        w.problem.nRows = 1;
        w.buffer2 = vec![0.0; 2048];
        w.buffer3 = vec![0.0; 2048];
        w.stepsize.dPrimalStep = 0.5;
        w.stepsize.dDualStep = 2.0;
        w.stepsize.dBeta = 4.0;
        w.iterates.x = [
            CUPDLPvec {
                len: 1,
                data: vec![0.0],
            },
            CUPDLPvec {
                len: 1,
                data: vec![0.0],
            },
        ];
        w.iterates.y = [
            CUPDLPvec {
                len: 1,
                data: vec![0.0],
            },
            CUPDLPvec {
                len: 1,
                data: vec![0.0],
            },
        ];
        w.iterates.xLastRestart = vec![0.0];
        w.iterates.yLastRestart = vec![0.0];
        // movement below 1e-10: beta unchanged, steps = 1/2 and 2.
        PDHG_Compute_Step_Size_Ratio(&mut w);
        assert_eq!(w.stepsize.dBeta, 4.0);
        assert_eq!(w.stepsize.dPrimalStep, 0.5);
        assert_eq!(w.stepsize.dDualStep, 2.0);
        assert_eq!(w.stepsize.dTheta, 1.0);
    }

    #[test]
    fn init_step_sizes_adaptive_branch() {
        let mut w = CUPDLPwork::default();
        w.stepsize.eLineSearchMethod = PDHG_ADAPTIVE_LINESEARCH;
        w.problem.nCols = 1;
        w.problem.nRows = 1;
        w.problem.cost = vec![2.0];
        w.problem.rhs = vec![1.0];
        w.problem.data.csc_matrix = Some(Box::new(CUPDLPcsc {
            MatElemNormInf: 2.0,
            ..Default::default()
        }));
        w.iterates.iLastRestartIter = 7;
        w.stepsize.dSumPrimalStep = 3.0;
        assert_eq!(PDHG_Init_Step_Sizes(&mut w), RETCODE_OK);
        // a = 4, b = 1 -> beta = 4; primal = (1 / 2) / 2; dual = primal * beta.
        assert_eq!(w.stepsize.dBeta, 4.0);
        assert_eq!(w.stepsize.dPrimalStep, 0.25);
        assert_eq!(w.stepsize.dDualStep, 1.0);
        assert_eq!(w.stepsize.dTheta, 1.0);
        assert_eq!(w.iterates.dLastRestartBeta, 4.0);
        assert_eq!(w.iterates.iLastRestartIter, 0);
        assert_eq!(w.stepsize.dSumPrimalStep, 0.0);
    }

    #[test]
    fn update_average_accumulates_mean_step() {
        let mut w = CUPDLPwork::default();
        w.problem.data.nCols = 1;
        w.problem.data.nRows = 1;
        w.stepsize.dPrimalStep = 0.5;
        w.stepsize.dDualStep = 2.0;
        w.timers.nIter = 0;
        w.iterates.x[1] = CUPDLPvec {
            len: 1,
            data: vec![3.0],
        };
        w.iterates.y[1] = CUPDLPvec {
            len: 1,
            data: vec![-2.0],
        };
        w.iterates.xSum = vec![1.0];
        w.iterates.ySum = vec![1.0];
        PDHG_Update_Average(&mut w);
        assert_eq!(w.iterates.xSum, vec![4.0]);
        assert_eq!(w.iterates.ySum, vec![-1.0]);
        assert_eq!(w.stepsize.dSumPrimalStep, 1.0);
        assert_eq!(w.stepsize.dSumDualStep, 1.0);
    }
}
