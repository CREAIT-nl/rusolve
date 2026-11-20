//! HiGHS source: pdlp/cupdlp/cupdlp_linalg.h, pdlp/cupdlp/cupdlp_linalg.c (lines 17-390)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The CPU primitives (`USE_MY_BLAS` branch). The CPU / GPU dispatching
//! wrappers (lines 457-802) are in [`part2`]; the GPU functions 394-456 are
//! not compiled in the CPU build.
//!
//! Ownership: a function that takes
//! `CUPDLPwork *w` only to reach one or two of its members takes those
//! members instead, so that call sites can pass `w->iterates` vectors and
//! `w->problem` at the same time. `ScatterCol` / `ScatterRow` / `AxCPU` /
//! `ATyCPU` read only `w->problem`.
//!
//! The `incx` / `incy` parameters are kept and `assert`ed to be 1 as in C.

pub mod part2;
pub use part2::*;

use super::defs::*;

/// HiGHS: cupdlp_linalg.c::ScatterCol (lines 17-24)
///
/// `w` replaced by `problem`.
pub fn ScatterCol(
    problem: &CUPDLPproblem,
    iCol: cupdlp_int,
    multiplier: cupdlp_float,
    target: &mut [cupdlp_float],
) {
    let matrix = problem.data.csc_matrix.as_ref().unwrap();

    for p in matrix.colMatBeg[iCol as usize]..matrix.colMatBeg[(iCol + 1) as usize] {
        target[matrix.colMatIdx[p as usize] as usize] += matrix.colMatElem[p as usize] * multiplier;
    }
}

/// HiGHS: cupdlp_linalg.c::ScatterRow (lines 26-33)
///
/// `w` replaced by `problem`.
pub fn ScatterRow(
    problem: &CUPDLPproblem,
    iRow: cupdlp_int,
    multiplier: cupdlp_float,
    target: &mut [cupdlp_float],
) {
    let matrix = problem.data.csr_matrix.as_ref().unwrap();

    for p in matrix.rowMatBeg[iRow as usize]..matrix.rowMatBeg[(iRow + 1) as usize] {
        target[matrix.rowMatIdx[p as usize] as usize] += matrix.rowMatElem[p as usize] * multiplier;
    }
}

/// HiGHS: cupdlp_linalg.c::AxCPU (lines 35-71)
///
/// `w` replaced by `problem`.
pub fn AxCPU(problem: &CUPDLPproblem, ax: &mut [cupdlp_float], x: &[cupdlp_float]) {
    // #if PDHG_USE_TIMERS (lines 36-39) dropped.

    // memset(ax, 0, sizeof(cupdlp_float) * lp->data->nRows); (line 62)
    for i in 0..problem.data.nRows as usize {
        ax[i] = 0.0;
    }

    for iCol in 0..problem.data.nCols {
        ScatterCol(problem, iCol, x[iCol as usize], ax);
    }

    // #if PDHG_USE_TIMERS (lines 68-70) dropped.
}

/// HiGHS: cupdlp_linalg.c::ATyCPU (lines 73-109)
///
/// `w` replaced by `problem`.
pub fn ATyCPU(problem: &CUPDLPproblem, aty: &mut [cupdlp_float], y: &[cupdlp_float]) {
    // #if PDHG_USE_TIMERS (lines 74-77) dropped.

    // memset(aty, 0, sizeof(cupdlp_float) * lp->data->nCols); (line 101)
    for i in 0..problem.data.nCols as usize {
        aty[i] = 0.0;
    }

    for iRow in 0..problem.data.nRows {
        ScatterRow(problem, iRow, y[iRow as usize], aty);
    }

    // #if PDHG_USE_TIMERS (lines 106-108) dropped.
}

/// HiGHS: cupdlp_linalg.c::nrm2 (lines 111-125)
pub fn nrm2(n: cupdlp_int, x: &[f64], incx: cupdlp_int) -> f64 {
    debug_assert!(incx == 1);

    let mut nrm = 0.0;

    for i in 0..n {
        nrm += x[i as usize] * x[i as usize];
    }

    nrm.sqrt()
}

/// HiGHS: cupdlp_linalg.c::nrminf (lines 127-142)
pub fn nrminf(n: cupdlp_int, x: &[f64], incx: cupdlp_int) -> f64 {
    debug_assert!(incx == 1);

    let mut nrm = 0.0;

    for i in 0..n {
        let tmp = x[i as usize].abs();
        if tmp > nrm {
            nrm = tmp;
        }
    }

    nrm
}

/// HiGHS: cupdlp_linalg.c::nrminfindex (lines 144-163) (`static`)
pub fn nrminfindex(n: cupdlp_int, x: &[f64], incx: cupdlp_int) -> cupdlp_int {
    debug_assert!(incx == 1);

    let mut nrm = 0.0;
    let mut index = 0;

    for i in 0..n {
        let tmp = x[i as usize].abs();
        if tmp > nrm {
            nrm = tmp;
            index = i;
        }
    }

    index
}

/// HiGHS: cupdlp_linalg.c::twoNorm (line 165)
pub fn twoNorm(x: &[f64], n: cupdlp_int) -> f64 {
    nrm2(n, x, 1)
}

/// HiGHS: cupdlp_linalg.c::twoNormSquared (line 167)
///
/// `pow(twoNorm(x, n), 2)`: `powf(2.0)`, never `x * x` by hand.
pub fn twoNormSquared(x: &[f64], n: cupdlp_int) -> f64 {
    twoNorm(x, n).powf(2.0)
}

/// HiGHS: cupdlp_linalg.c::infNorm (line 169)
pub fn infNorm(x: &[f64], n: cupdlp_int) -> f64 {
    nrminf(n, x, 1)
}

/// HiGHS: cupdlp_linalg.c::GenNorm (lines 172-186)
pub fn GenNorm(x: &[f64], n: cupdlp_int, p: cupdlp_float) -> f64 {
    if p == 2.0 {
        twoNorm(x, n)
    } else if p == f64::INFINITY {
        infNorm(x, n)
    } else {
        let mut nrm = 0.0;

        for i in 0..n {
            nrm += x[i as usize].abs().powf(p);
        }

        nrm.powf(1.0 / p)
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_cdot (lines 189-193)
pub fn cupdlp_cdot(x: &mut [cupdlp_float], y: &[cupdlp_float], len: cupdlp_int) {
    for i in 0..len {
        x[i as usize] *= y[i as usize];
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_cdiv (lines 196-200)
pub fn cupdlp_cdiv(x: &mut [cupdlp_float], y: &[cupdlp_float], len: cupdlp_int) {
    for i in 0..len {
        x[i as usize] /= y[i as usize];
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projLowerBound (lines 213-218)
pub fn cupdlp_projLowerBound(x: &mut [cupdlp_float], lb: &[cupdlp_float], len: cupdlp_int) {
    for i in 0..len {
        x[i as usize] = if x[i as usize] > lb[i as usize] {
            x[i as usize]
        } else {
            lb[i as usize]
        };
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projUpperBound (lines 221-226)
pub fn cupdlp_projUpperBound(x: &mut [cupdlp_float], ub: &[cupdlp_float], len: cupdlp_int) {
    for i in 0..len {
        x[i as usize] = if x[i as usize] < ub[i as usize] {
            x[i as usize]
        } else {
            ub[i as usize]
        };
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projSameLowerBound (lines 229-234)
pub fn cupdlp_projSameLowerBound(x: &mut [cupdlp_float], lb: cupdlp_float, len: cupdlp_int) {
    for i in 0..len {
        x[i as usize] = if x[i as usize] > lb {
            x[i as usize]
        } else {
            lb
        };
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projSameUpperBound (lines 237-242)
pub fn cupdlp_projSameUpperBound(x: &mut [cupdlp_float], ub: cupdlp_float, len: cupdlp_int) {
    for i in 0..len {
        x[i as usize] = if x[i as usize] < ub {
            x[i as usize]
        } else {
            ub
        };
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projPositive (lines 245-247)
pub fn cupdlp_projPositive(x: &mut [cupdlp_float], len: cupdlp_int) {
    cupdlp_projSameLowerBound(x, 0.0, len);
}

/// HiGHS: cupdlp_linalg.c::cupdlp_projNegative (lines 250-252)
pub fn cupdlp_projNegative(x: &mut [cupdlp_float], len: cupdlp_int) {
    cupdlp_projSameUpperBound(x, 0.0, len);
}

/// HiGHS: cupdlp_linalg.c::dot (lines 318-333)
pub fn dot(
    n: cupdlp_int,
    x: &[cupdlp_float],
    incx: cupdlp_int,
    y: &[cupdlp_float],
    incy: cupdlp_int,
) -> f64 {
    debug_assert!(incx == 1 && incy == 1);

    let mut dres = 0.0;

    for i in 0..n {
        dres += x[i as usize] * y[i as usize];
    }

    dres
}

/// HiGHS: cupdlp_linalg.c::Dotprod (lines 335-337)
pub fn Dotprod(x: &[cupdlp_float], y: &[cupdlp_float], n: cupdlp_int) -> f64 {
    dot(n, x, 1, y, 1)
}

/// HiGHS: cupdlp_linalg.c::Dotprod_Neumaier (lines 339-341) (plain `dot` upstream)
pub fn Dotprod_Neumaier(x: &[cupdlp_float], y: &[cupdlp_float], n: cupdlp_int) -> f64 {
    dot(n, x, 1, y, 1)
}

/// HiGHS: cupdlp_linalg.c::AddToVector (lines 344-355)
pub fn AddToVector(
    x: &mut [cupdlp_float],
    weight: cupdlp_float,
    y: &[cupdlp_float],
    n: cupdlp_int,
) {
    for i in 0..n {
        x[i as usize] += weight * y[i as usize];
    }
}

/// HiGHS: cupdlp_linalg.c::ScaleVector (lines 358-368)
pub fn ScaleVector(weight: cupdlp_float, x: &mut [cupdlp_float], n: cupdlp_int) {
    for i in 0..n {
        x[i as usize] *= weight;
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_filter_lower_bound (lines 371-376) (`static`)
pub fn cupdlp_filter_lower_bound(
    x: &mut [cupdlp_float],
    lb: &[cupdlp_float],
    bound: cupdlp_float,
    len: cupdlp_int,
) {
    for i in 0..len {
        x[i as usize] = if lb[i as usize] > bound {
            lb[i as usize]
        } else {
            0.0
        };
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_filter_upper_bound (lines 378-383) (`static`)
pub fn cupdlp_filter_upper_bound(
    x: &mut [cupdlp_float],
    ub: &[cupdlp_float],
    bound: cupdlp_float,
    len: cupdlp_int,
) {
    for i in 0..len {
        x[i as usize] = if ub[i as usize] < bound {
            ub[i as usize]
        } else {
            0.0
        };
    }
}

/// HiGHS: cupdlp_linalg.c::cupdlp_init_vector (lines 385-390) (`static`)
pub fn cupdlp_init_vector(x: &mut [cupdlp_float], val: cupdlp_float, len: cupdlp_int) {
    for i in 0..len {
        x[i as usize] = val;
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// 2x2 CSC matrix [[1 2] [0 3]].
    fn problem_2x2() -> CUPDLPproblem {
        let mut p = CUPDLPproblem::default();
        p.nRows = 2;
        p.nCols = 2;
        p.data.nRows = 2;
        p.data.nCols = 2;
        p.data.csc_matrix = Some(Box::new(CUPDLPcsc {
            nRows: 2,
            nCols: 2,
            nMatElem: 3,
            colMatBeg: vec![0, 1, 3],
            colMatIdx: vec![0, 0, 1],
            colMatElem: vec![1.0, 2.0, 3.0],
            MatElemNormInf: 3.0,
        }));
        p.data.csr_matrix = Some(Box::new(CUPDLPcsr {
            nRows: 2,
            nCols: 2,
            nMatElem: 3,
            rowMatBeg: vec![0, 2, 3],
            rowMatIdx: vec![0, 1, 1],
            rowMatElem: vec![1.0, 2.0, 3.0],
        }));
        p
    }

    #[test]
    fn ax_and_aty_on_small_matrix() {
        let p = problem_2x2();
        let mut ax = vec![7.0; 2];
        AxCPU(&p, &mut ax, &[1.0, 1.0]);
        assert_eq!(ax, vec![3.0, 3.0]);
        let mut aty = vec![7.0; 2];
        ATyCPU(&p, &mut aty, &[1.0, 1.0]);
        assert_eq!(aty, vec![1.0, 5.0]);
    }

    #[test]
    fn norms_and_projections() {
        let x = [3.0, -4.0];
        assert_eq!(twoNorm(&x, 2), 5.0);
        assert_eq!(twoNormSquared(&x, 2), 25.0);
        assert_eq!(infNorm(&x, 2), 4.0);
        assert_eq!(nrminfindex(2, &x, 1), 1);
        assert_eq!(GenNorm(&x, 2, f64::INFINITY), 4.0);
        assert_eq!(dot(2, &x, 1, &x, 1), 25.0);
        let mut v = vec![-1.0, 2.0];
        cupdlp_projPositive(&mut v, 2);
        assert_eq!(v, vec![0.0, 2.0]);
        let mut f = vec![9.0; 2];
        cupdlp_filter_lower_bound(&mut f, &[-1e30, 1.0], -1e20, 2);
        assert_eq!(f, vec![0.0, 1.0]);
    }
}
