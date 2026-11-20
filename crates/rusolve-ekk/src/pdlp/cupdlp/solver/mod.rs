//! HiGHS source: pdlp/cupdlp/cupdlp_solver.h, pdlp/cupdlp/cupdlp_solver.c (lines 12-311)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The PDHG driver. `cupdlp_solver.c` is 1498 lines, split by contiguous
//! line range:
//!
//! | Rust module  | lines     | contents |
//! |--------------|-----------|----------|
//! | this file    | 12-311    | primal / dual feasibility, primal infeasibility |
//! | [`part2`]    | 313-591   | dual infeasibility, residuals, `Init_Variables` |
//! | [`check`]    | 598-897   | `Check_Data`, termination and infeasibility tests, printers |
//! | [`solve`]    | 899-1215  | `PDHG_Solve` (the main loop) |
//! | [`post`]     | 1217-1498 | `PDHG_PreSolve`, `PDHG_PostSolve`, `LP_SolvePDHG` |
//!
//! Borrowing: the four `Compute_*_Feasibility` /
//! `*_Infeasibility` functions receive `work` together with pointers into
//! `*work` (`resobj->primalResidual`, `iterates->ax[k]->data`, `&resobj->...`).
//! The caller `std::mem::take`s every vector argument out of `work`, passes
//! `&mut CUPDLPwork` plus the slices, and moves them back; scalar
//! out-parameters are locals assigned back after the call. The callee never
//! reads those fields through `work` (verified per function).

pub mod check;
pub mod part2;
pub mod post;
pub mod solve;

pub use check::*;
pub use part2::*;
pub use post::*;
pub use solve::*;

use super::defs::*;
use super::linalg::*;

/// HiGHS: cupdlp_solver.c::PDHG_Compute_Primal_Feasibility (lines 12-67)
///
/// CPU branch; `CUPDLP_DEBUG` blocks dropped.
pub fn PDHG_Compute_Primal_Feasibility(
    work: &mut CUPDLPwork,
    primalResidual: &mut [cupdlp_float],
    ax: &[cupdlp_float],
    x: &[cupdlp_float],
    dPrimalFeasibility: &mut cupdlp_float,
    dPrimalObj: &mut cupdlp_float,
) {
    let nRows = work.problem.data.nRows;
    let nCols = work.problem.data.nCols;
    let nEqs = work.problem.nEqs;
    let ifScaled = work.scaling.ifScaled;

    // primal variable violation

    // todo, add this
    //    *dPrimalObj = Dotprod_Neumaier(problem->cost, x, lp->nCols);
    cupdlp_dot(nCols, x, &work.problem.cost, dPrimalObj);
    *dPrimalObj = *dPrimalObj * work.problem.sense_origin + work.problem.offset;

    // GPU branch dropped (lines 27-35: #if !defined(CUPDLP_CPU) && USE_KERNELS)

    // CPU branch (lines 36-66):

    // cupdlp_copy(primalResidual, ax, cupdlp_float, lp->nRows);
    cupdlp_copy(primalResidual, ax, nRows);

    // AddToVector(primalResidual, -1.0, problem->rhs, lp->nRows);
    let alpha = -1.0;
    cupdlp_axpy(nRows, alpha, &work.problem.rhs, primalResidual);

    //  double dPrimalFeas = 0.0; // Redundant

    // todo, check
    //  cupdlp_projNegative(primalResidual + problem->nEqs, primalResidual +
    //  problem->nEqs, lp->nRows - problem->nEqs);
    //

    cupdlp_projNeg(&mut primalResidual[nEqs as usize..], nRows - nEqs);

    if ifScaled != 0 {
        // cupdlp_edot(primalResidual, scaling->rowScale, lp->nRows);
        // cupdlp_edot(primalResidual, scaling->rowScale_gpu, lp->nRows);

        // cupdlp_copy_vec(work->buffer3, scaling->rowScale, cupdlp_float,
        // lp->nRows); cupdlp_edot(primalResidual, work->buffer3, lp->nRows);

        cupdlp_edot(primalResidual, &work.rowScale, nRows);
    }

    cupdlp_twoNorm(nRows, primalResidual, dPrimalFeasibility);
}

/// HiGHS: cupdlp_solver.c::PDHG_Compute_Dual_Feasibility (lines 69-204)
///
/// CPU branch; `CUPDLP_DEBUG` blocks dropped. `work->buffer2` is used only by
/// the dropped GPU branch; `x` and `dComplementarity` are unused in the CPU
/// branch (as upstream).
#[allow(clippy::too_many_arguments)]
pub fn PDHG_Compute_Dual_Feasibility(
    work: &mut CUPDLPwork,
    dualResidual: &mut [cupdlp_float],
    aty: &[cupdlp_float],
    x: &[cupdlp_float],
    y: &[cupdlp_float],
    dDualFeasibility: &mut cupdlp_float,
    dDualObj: &mut cupdlp_float,
    dComplementarity: &mut cupdlp_float,
    dSlackPos: &mut [cupdlp_float],
    dSlackNeg: &mut [cupdlp_float],
) {
    let nRows = work.problem.data.nRows;
    let nCols = work.problem.data.nCols;
    let ifScaled = work.scaling.ifScaled;

    // todo, compute Neumaier
    //    *dDualObj = Dotprod_Neumaier(problem->rhs, y, lp->nRows);
    cupdlp_dot(nRows, y, &work.problem.rhs, dDualObj);

    // GPU branch dropped (lines 82-113: #if !defined(CUPDLP_CPU) && USE_KERNELS)

    // CPU branch (lines 114-203):

    // *dComplementarity = 0.0;

    // @note:
    // original dual residual in pdlp:
    // they compute:
    //    violation   +  reduced cost
    //  |max(-y, 0)|  + |(I-\Pi)(c-Α'\nu)|
    // compute c - A'y

    cupdlp_copy(dualResidual, aty, nCols);
    let alpha = -1.0;
    cupdlp_scaleVector(alpha, dualResidual, nCols);

    let alpha = 1.0;
    cupdlp_axpy(nCols, alpha, &work.problem.cost, dualResidual);

    //    julia version
    //        function compute_reduced_costs_from_primal_gradient_kernel!(
    //            primal_gradient::CuDeviceVector{Float64},
    //            isfinite_variable_lower_bound::CuDeviceVector{Bool},
    //            isfinite_variable_upper_bound::CuDeviceVector{Bool},
    //            num_variables::Int64,
    //            reduced_costs::CuDeviceVector{Float64},
    //            reduced_costs_violation::CuDeviceVector{Float64},
    //        )
    //            tx = threadIdx().x + (blockDim().x * (blockIdx().x - 0x1))
    //            if tx <= num_variables
    //                @inbounds begin
    //                    reduced_costs[tx] = max(primal_gradient[tx], 0.0) *
    //                    isfinite_variable_lower_bound[tx] +
    //                    min(primal_gradient[tx], 0.0) *
    //                    isfinite_variable_upper_bound[tx]
    //
    //                    reduced_costs_violation[tx] = primal_gradient[tx] -
    //                    reduced_costs[tx]
    //                end
    //            end
    //            return
    //        end

    cupdlp_copy(dSlackPos, dualResidual, nCols);

    cupdlp_projPos(dSlackPos, nCols);

    cupdlp_edot(dSlackPos, &work.problem.hasLower, nCols);

    let mut temp = 0.0;
    // cupdlp_dot(work, lp->nCols, x, dSlackPos, &temp);
    // *dComplementarity += temp;
    // cupdlp_dot(work, lp->nCols, dSlackPos, resobj->dLowerFiltered, &temp);
    // *dComplementarity -= temp;
    cupdlp_dot(nCols, dSlackPos, &work.resobj.dLowerFiltered, &mut temp);
    *dDualObj += temp;

    cupdlp_copy(dSlackNeg, dualResidual, nCols);

    cupdlp_projNeg(dSlackNeg, nCols);

    cupdlp_scaleVector(-1.0, dSlackNeg, nCols);

    cupdlp_edot(dSlackNeg, &work.problem.hasUpper, nCols);

    // cupdlp_dot(work, lp->nCols, x, dSlackNeg, &temp);
    // *dComplementarity -= temp;
    // cupdlp_dot(work, lp->nCols, dSlackNeg, resobj->dUpperFiltered, &temp);
    // *dComplementarity += temp;
    cupdlp_dot(nCols, dSlackNeg, &work.resobj.dUpperFiltered, &mut temp);
    *dDualObj -= temp;

    *dDualObj = *dDualObj * work.problem.sense_origin + work.problem.offset;

    let alpha = -1.0;
    cupdlp_axpy(nCols, alpha, dSlackPos, dualResidual);
    let alpha = 1.0;
    cupdlp_axpy(nCols, alpha, dSlackNeg, dualResidual);

    if ifScaled != 0 {
        // cupdlp_edot(dualResidual, scaling->colScale, lp->nCols);
        // cupdlp_edot(dualResidual, scaling->colScale_gpu, lp->nCols);

        // cupdlp_copy_vec(work->buffer3, scaling->colScale, cupdlp_float,
        // lp->nCols); cupdlp_edot(dualResidual, work->buffer3, lp->nCols);

        cupdlp_edot(dualResidual, &work.colScale, nCols);
    }

    cupdlp_twoNorm(nCols, dualResidual, dDualFeasibility);
}

/// HiGHS: cupdlp_solver.c::PDHG_Compute_Primal_Infeasibility (lines 206-311)
///
/// CPU branch. Writes the ray vectors `resobj->dualInfeasRay` /
/// `dualInfeasLbRay` / `dualInfeasUbRay` / `dualInfeasConstr` through `work`.
#[allow(clippy::too_many_arguments)]
pub fn PDHG_Compute_Primal_Infeasibility(
    work: &mut CUPDLPwork,
    y: &[cupdlp_float],
    dSlackPos: &[cupdlp_float],
    dSlackNeg: &[cupdlp_float],
    aty: &[cupdlp_float],
    dualObj: cupdlp_float,
    dPrimalInfeasObj: &mut cupdlp_float,
    dPrimalInfeasRes: &mut cupdlp_float,
) {
    let nRows = work.problem.data.nRows;
    let nCols = work.problem.data.nCols;
    let ifScaled = work.scaling.ifScaled;

    let mut yNrmSq = 1.0;
    let mut slackPosNrmSq = 1.0;
    let mut slackNegSq = 1.0;
    // C `double dScale = 1.0;` is always overwritten by the sqrt below before it is read.
    let mut dScale;

    // cupdlp_float dConstrResSq = 0.0;
    // y and lambda must be feasible, no need to check bound
    // cupdlp_float dBoundLbResSq = 0.0;
    // cupdlp_float dBoundUbResSq = 0.0;

    // GPU branch dropped (lines 229-256: #if !defined(CUPDLP_CPU) && USE_KERNELS)

    // CPU branch (lines 257-309):

    // y, ldb ray
    cupdlp_copy(&mut work.resobj.dualInfeasRay, y, nRows);
    cupdlp_copy(&mut work.resobj.dualInfeasLbRay, dSlackPos, nCols);
    cupdlp_copy(&mut work.resobj.dualInfeasUbRay, dSlackNeg, nCols);
    cupdlp_twoNormSquared(nRows, &work.resobj.dualInfeasRay, &mut yNrmSq);
    cupdlp_twoNormSquared(nCols, &work.resobj.dualInfeasLbRay, &mut slackPosNrmSq);
    cupdlp_twoNormSquared(nCols, &work.resobj.dualInfeasUbRay, &mut slackNegSq);
    dScale = (yNrmSq + slackPosNrmSq + slackNegSq).sqrt();
    // dScale /= sqrt(problem->data->nRows + 2 * problem->data->nCols);
    if dScale < 1e-8 {
        dScale = 1.0;
    }
    cupdlp_scaleVector(1.0 / dScale, &mut work.resobj.dualInfeasRay, nRows);
    cupdlp_scaleVector(1.0 / dScale, &mut work.resobj.dualInfeasLbRay, nCols);
    cupdlp_scaleVector(1.0 / dScale, &mut work.resobj.dualInfeasUbRay, nCols);

    // dual obj
    *dPrimalInfeasObj = (dualObj - work.problem.offset) / work.problem.sense_origin / dScale;

    // dual constraints [ATy1 + GTy2 + lambda]
    cupdlp_copy(&mut work.resobj.dualInfeasConstr, aty, nCols);
    cupdlp_scaleVector(1.0 / dScale, &mut work.resobj.dualInfeasConstr, nCols);
    let alpha = 1.0;
    cupdlp_axpy(
        nCols,
        alpha,
        &work.resobj.dualInfeasLbRay,
        &mut work.resobj.dualInfeasConstr,
    );
    let alpha = -1.0;
    cupdlp_axpy(
        nCols,
        alpha,
        &work.resobj.dualInfeasUbRay,
        &mut work.resobj.dualInfeasConstr,
    );
    if ifScaled != 0 {
        cupdlp_edot(&mut work.resobj.dualInfeasConstr, &work.colScale, nCols);
    }
    // cupdlp_twoNormSquared(work, problem->data->nCols, resobj->dualInfeasConstr,
    //                       &dConstrResSq);
    cupdlp_twoNorm(nCols, &work.resobj.dualInfeasConstr, dPrimalInfeasRes);

    // dual bound
    // always satisfied, no need to check
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// min x0 s.t. x0 >= 1 (one inequality row, unscaled); x = 2.
    fn one_row_work() -> CUPDLPwork {
        let mut w = CUPDLPwork::default();
        w.problem.nRows = 1;
        w.problem.nCols = 1;
        // C reads `lp = problem->data; lp->nRows / lp->nCols`, so the sizes
        // must also be set on `problem.data` (test fix, see commit message).
        w.problem.data.nRows = 1;
        w.problem.data.nCols = 1;
        w.problem.nEqs = 0;
        w.problem.cost = vec![1.0];
        w.problem.rhs = vec![1.0];
        w.problem.lower = vec![0.0];
        w.problem.upper = vec![f64::INFINITY];
        w.problem.hasLower = vec![1.0];
        w.problem.hasUpper = vec![0.0];
        w.problem.sense_origin = 1.0;
        w.problem.offset = 0.0;
        w.buffer2 = vec![0.0; 2048];
        w.buffer3 = vec![0.0; 2048];
        w
    }

    #[test]
    fn primal_feasibility_projects_inequality_residual() {
        let mut w = one_row_work();
        let mut r = vec![0.0];
        let (mut feas, mut obj) = (0.0, 0.0);
        // Ax - b = 2 - 1 = 1 on an inequality row -> projNeg -> 0.
        PDHG_Compute_Primal_Feasibility(&mut w, &mut r, &[2.0], &[2.0], &mut feas, &mut obj);
        assert_eq!(obj, 2.0);
        assert_eq!(r, vec![0.0]);
        assert_eq!(feas, 0.0);
    }

    #[test]
    fn primal_infeasibility_zero_ray_has_zero_objective() {
        let mut w = one_row_work();
        w.resobj.primalInfeasRay = vec![0.0];
        w.resobj.primalInfeasConstr = vec![0.0];
        w.resobj.primalInfeasBound = vec![0.0];
        // The C writes `resobj->dualInfeas*` unconditionally (CPU branch); they
        // are allocated by PDHG_Alloc, so the test allocates them too.
        w.resobj.dualInfeasRay = vec![0.0];
        w.resobj.dualInfeasLbRay = vec![0.0];
        w.resobj.dualInfeasUbRay = vec![0.0];
        w.resobj.dualInfeasConstr = vec![0.0];
        let (mut obj, mut res) = (1.0, 1.0);
        PDHG_Compute_Primal_Infeasibility(
            &mut w,
            &[0.0],
            &[0.0],
            &[0.0],
            &[0.0],
            0.0,
            &mut obj,
            &mut res,
        );
        assert_eq!(obj, 0.0);
    }
}
