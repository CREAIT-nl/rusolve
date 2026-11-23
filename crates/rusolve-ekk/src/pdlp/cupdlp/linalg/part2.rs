//! HiGHS source: pdlp/cupdlp/cupdlp_linalg.c (lines 457-802)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The "Apis compatible with both CPU and GPU", CPU branch only. Every
//! `#ifndef CUPDLP_CPU` branch is dropped; in `Ax` / `ATy` the
//! `SINGLE_GPU` / `MULTI_GPU` / default cases of the CPU build
//! (`printf(...); exit(1);`) become `panic!`.
//!
//! Ownership: the CPU bodies of `cupdlp_axpy`,
//! `cupdlp_dot`, `cupdlp_twoNorm`, `cupdlp_infNormIndex`,
//! `cupdlp_scaleVector` and `cupdlp_twoNormSquared` never read `w`, so the
//! parameter is dropped; `cupdlp_diffTwoNormSquared` / `cupdlp_diffTwoNorm`
//! take `w->buffer2` and `cupdlp_diffDotDiff` takes `w->buffer2` and
//! `w->buffer3`. `cupdlp_infNorm` is declared in the header but never
//! defined, so it has no port. `cupdlp_sub` (760-769) is commented out
//! upstream.

use super::super::defs::*;
use super::super::utils::getTimeStamp;
use super::*;

/// HiGHS: cupdlp_linalg.c::Ax (lines 460-494)
///
/// `w` replaced by `problem` and `timers`. The
/// `PDHG_USE_TIMERS` block (`dAxTime += getTimeStamp() - begin;
/// nAxCalls++`) is kept.
pub fn Ax(problem: &CUPDLPproblem, timers: &mut CUPDLPtimers, ax: &mut CUPDLPvec, x: &CUPDLPvec) {
    let begin = getTimeStamp();

    let d = &problem.data;
    match d.device {
        CPU => {
            AxCPU(problem, &mut ax.data, &x.data);
        }
        SINGLE_GPU => {
            panic!("GPU not supported in CPU build");
        }
        MULTI_GPU => {
            panic!("GPU not supported in CPU build");
        }
        _ => {
            panic!("Error: Unknown device type in Ax");
        }
    }

    // PDHG_USE_TIMERS
    timers.dAxTime += getTimeStamp() - begin;
    timers.nAxCalls += 1;
}

/// HiGHS: cupdlp_linalg.c::ATy (lines 496-530)
///
/// `w` replaced by `problem` and `timers`.
pub fn ATy(problem: &CUPDLPproblem, timers: &mut CUPDLPtimers, aty: &mut CUPDLPvec, y: &CUPDLPvec) {
    let begin = getTimeStamp();

    let d = &problem.data;
    match d.device {
        CPU => {
            ATyCPU(problem, &mut aty.data, &y.data);
        }
        SINGLE_GPU => {
            panic!("GPU not supported in CPU build");
        }
        MULTI_GPU => {
            panic!("GPU not supported in CPU build");
        }
        _ => {
            panic!("Error: Unknown device type in ATy");
        }
    }

    // PDHG_USE_TIMERS
    timers.dAtyTime += getTimeStamp() - begin;
    timers.nAtyCalls += 1;
}

/// HiGHS: cupdlp_linalg.c::cupdlp_axpy (lines 537-551)
///
/// `w` dropped; `const cupdlp_float *alpha` passed by value.
pub fn cupdlp_axpy(
    n: cupdlp_int,
    alpha: cupdlp_float,
    x: &[cupdlp_float],
    y: &mut [cupdlp_float],
) -> cupdlp_int {
    // CPU branch: AddToVector(y, *alpha, x, n);
    AddToVector(y, alpha, x, n);
    0
}

/// HiGHS: cupdlp_linalg.c::cupdlp_dot (lines 553-565)
///
/// `w` dropped.
pub fn cupdlp_dot(
    n: cupdlp_int,
    x: &[cupdlp_float],
    y: &[cupdlp_float],
    res: &mut cupdlp_float,
) -> cupdlp_int {
    // CPU branch: *res = dot(n, x, 1, y, 1);
    *res = dot(n, x, 1, y, 1);
    0
}

/// HiGHS: cupdlp_linalg.c::cupdlp_twoNorm (lines 567-579)
///
/// `w` dropped.
pub fn cupdlp_twoNorm(n: cupdlp_int, x: &[cupdlp_float], res: &mut cupdlp_float) -> cupdlp_int {
    // CPU branch: *res = nrm2(n, x, 1);
    *res = nrm2(n, x, 1);
    0
}

/// HiGHS: cupdlp_linalg.c::cupdlp_infNormIndex (lines 581-593)
///
/// `w` dropped.
pub fn cupdlp_infNormIndex(n: cupdlp_int, x: &[cupdlp_float], res: &mut cupdlp_int) -> cupdlp_int {
    // CPU branch: *res = nrminfindex(n, x, 1);
    *res = nrminfindex(n, x, 1);
    0
}

/// HiGHS: cupdlp_linalg.c::cupdlp_scaleVector (lines 595-607)
///
/// `w` dropped.
pub fn cupdlp_scaleVector(
    weight: cupdlp_float,
    x: &mut [cupdlp_float],
    n: cupdlp_int,
) -> cupdlp_int {
    // CPU branch: ScaleVector(weight, x, n);
    ScaleVector(weight, x, n);
    0
}

/// HiGHS: cupdlp_linalg.c::cupdlp_twoNormSquared (lines 609-612)
///
/// `w` dropped. Upstream is `cupdlp_dot(w, n, x, x, res)`, not
/// `pow(nrm2, 2)` as `twoNormSquared`.
pub fn cupdlp_twoNormSquared(n: cupdlp_int, x: &[cupdlp_float], res: &mut cupdlp_float) {
    cupdlp_dot(n, x, x, res);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_diffTwoNormSquared (lines 615-622)
///
/// `w` replaced by `w->buffer2`.
pub fn cupdlp_diffTwoNormSquared(
    buffer2: &mut [cupdlp_float],
    x: &[cupdlp_float],
    y: &[cupdlp_float],
    len: cupdlp_int,
    res: &mut cupdlp_float,
) {
    cupdlp_copy(buffer2, x, len);
    let alpha: cupdlp_float = -1.0;
    cupdlp_axpy(len, alpha, y, buffer2);
    cupdlp_twoNormSquared(len, buffer2, res);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_diffTwoNorm (lines 625-632)
///
/// `w` replaced by `w->buffer2`.
pub fn cupdlp_diffTwoNorm(
    buffer2: &mut [cupdlp_float],
    x: &[cupdlp_float],
    y: &[cupdlp_float],
    len: cupdlp_int,
    res: &mut cupdlp_float,
) {
    cupdlp_copy(buffer2, x, len);
    let alpha: cupdlp_float = -1.0;
    cupdlp_axpy(len, alpha, y, buffer2);
    cupdlp_twoNorm(len, buffer2, res);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_diffDotDiff (lines 635-646)
///
/// `w` replaced by `w->buffer2` and `w->buffer3`.
pub fn cupdlp_diffDotDiff(
    buffer2: &mut [cupdlp_float],
    buffer3: &mut [cupdlp_float],
    x1: &[cupdlp_float],
    x2: &[cupdlp_float],
    y1: &[cupdlp_float],
    y2: &[cupdlp_float],
    len: cupdlp_int,
    res: &mut cupdlp_float,
) {
    cupdlp_copy(buffer2, x1, len);
    let alpha: cupdlp_float = -1.0;
    cupdlp_axpy(len, alpha, x2, buffer2);
    cupdlp_copy(buffer3, y1, len);
    cupdlp_axpy(len, alpha, y2, buffer3);
    // reduce step
    cupdlp_dot(len, buffer2, buffer3, res);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_edot (lines 651-657)
pub fn cupdlp_edot(x: &mut [cupdlp_float], y: &[cupdlp_float], len: cupdlp_int) {
    // CPU branch: cupdlp_cdot(x, y, len);
    cupdlp_cdot(x, y, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_ediv (lines 660-666)
pub fn cupdlp_ediv(x: &mut [cupdlp_float], y: &[cupdlp_float], len: cupdlp_int) {
    // CPU branch: cupdlp_cdiv(x, y, len);
    cupdlp_cdiv(x, y, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projlb (lines 668-675)
pub fn cupdlp_projlb(x: &mut [cupdlp_float], lb: &[cupdlp_float], len: cupdlp_int) {
    // CPU branch: cupdlp_projLowerBound(x, lb, len);
    cupdlp_projLowerBound(x, lb, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projub (lines 677-684)
pub fn cupdlp_projub(x: &mut [cupdlp_float], ub: &[cupdlp_float], len: cupdlp_int) {
    // CPU branch: cupdlp_projUpperBound(x, ub, len);
    cupdlp_projUpperBound(x, ub, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projSamelb (lines 686-693)
pub fn cupdlp_projSamelb(x: &mut [cupdlp_float], lb: cupdlp_float, len: cupdlp_int) {
    // CPU branch: cupdlp_projSameLowerBound(x, lb, len);
    cupdlp_projSameLowerBound(x, lb, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projSameub (lines 695-702)
pub fn cupdlp_projSameub(x: &mut [cupdlp_float], ub: cupdlp_float, len: cupdlp_int) {
    // CPU branch: cupdlp_projSameUpperBound(x, ub, len);
    cupdlp_projSameUpperBound(x, ub, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projPos (lines 705-707)
pub fn cupdlp_projPos(x: &mut [cupdlp_float], len: cupdlp_int) {
    cupdlp_projSamelb(x, 0.0, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projNeg (lines 710-712)
pub fn cupdlp_projNeg(x: &mut [cupdlp_float], len: cupdlp_int) {
    cupdlp_projSameub(x, 0.0, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_filterlb (lines 732-739)
pub fn cupdlp_filterlb(
    x: &mut [cupdlp_float],
    lb: &[cupdlp_float],
    bound: cupdlp_float,
    len: cupdlp_int,
) {
    // CPU branch: cupdlp_filter_lower_bound(x, lb, bound, len);
    cupdlp_filter_lower_bound(x, lb, bound, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_filterub (lines 741-748)
pub fn cupdlp_filterub(
    x: &mut [cupdlp_float],
    ub: &[cupdlp_float],
    bound: cupdlp_float,
    len: cupdlp_int,
) {
    // CPU branch: cupdlp_filter_upper_bound(x, ub, bound, len);
    cupdlp_filter_upper_bound(x, ub, bound, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_initvec (lines 750-757)
pub fn cupdlp_initvec(x: &mut [cupdlp_float], val: cupdlp_float, len: cupdlp_int) {
    // CPU branch: cupdlp_init_vector(x, val, len);
    cupdlp_init_vector(x, val, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_compute_interaction_and_movement (lines 772-801)
///
/// Takes the whole work: it reads `iterates`, `problem`, `stepsize`,
/// `timers` and writes `buffer2` / `buffer3`; the body splits the borrows
/// field by field.
pub fn cupdlp_compute_interaction_and_movement(
    w: &mut CUPDLPwork,
    dMovement: &mut cupdlp_float,
    dInteraction: &mut cupdlp_float,
) {
    let nCols = w.problem.nCols;
    let nRows = w.problem.nRows;
    let beta = w.stepsize.dBeta.sqrt();
    let mut dX: cupdlp_float = 0.0;
    let mut dY: cupdlp_float = 0.0;

    let iter = w.timers.nIter;
    let iterates = &w.iterates;
    let x = &iterates.x[(iter % 2) as usize];
    let y = &iterates.y[(iter % 2) as usize];
    let aty = &iterates.aty[(iter % 2) as usize];
    let xUpdate = &iterates.x[((iter + 1) % 2) as usize];
    let yUpdate = &iterates.y[((iter + 1) % 2) as usize];
    let atyUpdate = &iterates.aty[((iter + 1) % 2) as usize];

    // CPU branch (the `USE_KERNELS` path is GPU-only).
    cupdlp_diffTwoNormSquared(&mut w.buffer2, &x.data, &xUpdate.data, nCols, &mut dX);
    cupdlp_diffTwoNormSquared(&mut w.buffer2, &y.data, &yUpdate.data, nRows, &mut dY);
    //      Δx' (AΔy)
    cupdlp_diffDotDiff(
        &mut w.buffer2,
        &mut w.buffer3,
        &x.data,
        &xUpdate.data,
        &aty.data,
        &atyUpdate.data,
        nCols,
        dInteraction,
    );
    *dMovement = dX * 0.5 * beta + dY / (2.0 * beta);
}

#[cfg(test)]
mod tests {
    //! Unit tests, including
    //! `proj.rs::PDHG_Project_Bounds`).
    use super::*;
    use crate::pdlp::cupdlp::proj::PDHG_Project_Bounds;

    #[test]
    fn project_bounds_clamps_to_box() {
        let mut p = CUPDLPproblem::default();
        p.nCols = 3;
        p.lower = vec![0.0, f64::NEG_INFINITY, -1.0];
        p.upper = vec![1.0, 2.0, f64::INFINITY];
        let mut r = vec![-5.0, 5.0, -0.5];
        PDHG_Project_Bounds(&p, &mut r);
        assert_eq!(r, vec![0.0, 2.0, -0.5]);
    }

    #[test]
    fn diff_norms_use_buffer() {
        let mut b2 = vec![0.0; 4];
        let mut b3 = vec![0.0; 4];
        let mut r = 0.0;
        cupdlp_diffTwoNormSquared(&mut b2, &[4.0, 1.0], &[1.0, 5.0], 2, &mut r);
        assert_eq!(r, 25.0);
        cupdlp_diffTwoNorm(&mut b2, &[4.0, 1.0], &[1.0, 5.0], 2, &mut r);
        assert_eq!(r, 5.0);
        cupdlp_diffDotDiff(
            &mut b2,
            &mut b3,
            &[2.0, 0.0],
            &[1.0, 0.0],
            &[0.0, 3.0],
            &[0.0, 1.0],
            2,
            &mut r,
        );
        assert_eq!(r, 0.0);
        let mut y = vec![1.0, 1.0];
        assert_eq!(cupdlp_axpy(2, -2.0, &[1.0, 0.5], &mut y), 0);
        assert_eq!(y, vec![-1.0, 0.0]);
    }

    #[test]
    fn ax_counts_calls() {
        let mut p = CUPDLPproblem::default();
        p.nRows = 1;
        p.nCols = 1;
        p.data.nRows = 1;
        p.data.nCols = 1;
        p.data.device = CPU;
        p.data.csc_matrix = Some(Box::new(CUPDLPcsc {
            nRows: 1,
            nCols: 1,
            nMatElem: 1,
            colMatBeg: vec![0, 1],
            colMatIdx: vec![0],
            colMatElem: vec![2.0],
            MatElemNormInf: 2.0,
        }));
        let mut t = CUPDLPtimers::default();
        let mut ax = CUPDLPvec {
            len: 1,
            data: vec![0.0],
        };
        Ax(
            &p,
            &mut t,
            &mut ax,
            &CUPDLPvec {
                len: 1,
                data: vec![3.0],
            },
        );
        assert_eq!(ax.data, vec![6.0]);
        assert_eq!(t.nAxCalls, 1);
    }

    #[test]
    fn compute_interaction_uses_odd_iterate_parity() {
        // iter = 1: x[1] is the current iterate, x[0] the update (C: x[iter % 2]).
        let mut w = CUPDLPwork::default();
        w.problem.nRows = 1;
        w.problem.nCols = 2;
        w.timers.nIter = 1;
        w.stepsize.dBeta = 1.0;
        w.buffer2 = vec![0.0; 2];
        w.buffer3 = vec![0.0; 2];
        w.iterates.x[1].data = vec![1.0, 0.0];
        w.iterates.x[0].data = vec![1.0, 2.0];
        w.iterates.y[1].data = vec![0.0];
        w.iterates.y[0].data = vec![0.0];
        w.iterates.aty[1].data = vec![0.0, 0.0];
        w.iterates.aty[0].data = vec![0.0, 0.0];
        let mut movement = 0.0;
        let mut interaction = 0.0;
        cupdlp_compute_interaction_and_movement(&mut w, &mut movement, &mut interaction);
        // dX = ||(1,0) - (1,2)||^2 = 4, dY = 0, dMovement = 4 * 0.5 * 1 + 0 / 2 = 2.
        assert_eq!(movement, 2.0);
        assert_eq!(interaction, 0.0);
    }
}
