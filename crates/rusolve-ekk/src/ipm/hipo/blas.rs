//! HiGHS source: extern/blas/mycblas.h (interface), tools/highs-oracle/refblas/refblas.c (reference implementation)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The BLAS kernels FactorHiGHS calls (`HighsExtras::blas::*`). HiGHS links
//! whatever BLAS the build finds; the oracle is pinned to
//! `tools/highs-oracle/refblas/refblas.c` (MIT, written for rusolve, netlib
//! reference loop order, compiled `-O2 -ffp-contract=off`), and this module is
//! its statement-by-statement mirror, so both sides perform the same floating
//! point operations in the same order. Column-major
//! only (`CblasColMajor` is the only order HiPO passes); arguments are the
//! CBLAS ones with the enum values below. Array arguments are slices starting
//! at the C pointer; the kernels index them with the netlib `(i-1)*inc`
//! convention for negative increments.

use super::Int;

/// HiGHS: mycblas.h `enum CBLAS_ORDER`
pub const CblasRowMajor: i32 = 101;
pub const CblasColMajor: i32 = 102;
/// HiGHS: mycblas.h `enum CBLAS_TRANSPOSE`
pub const CblasNoTrans: i32 = 111;
pub const CblasTrans: i32 = 112;
pub const CblasConjTrans: i32 = 113;
/// HiGHS: mycblas.h `enum CBLAS_UPLO`
pub const CblasUpper: i32 = 121;
pub const CblasLower: i32 = 122;
/// HiGHS: mycblas.h `enum CBLAS_DIAG`
pub const CblasNonUnit: i32 = 131;
pub const CblasUnit: i32 = 132;
/// HiGHS: mycblas.h `enum CBLAS_SIDE`
pub const CblasLeft: i32 = 141;
pub const CblasRight: i32 = 142;

/// HiGHS: mycblas.h `typedef int32_t blasint`
pub type blasint = i32;

/// HiGHS: refblas.c::refblas_col_major_only (lines 40-46)
fn refblas_col_major_only(name: &str, order: i32) {
    if order != CblasColMajor {
        eprintln!(
            "refblas: {} called with order {} (only ColMajor)",
            name, order
        );
        panic!(
            "refblas: {} called with order {} (only ColMajor)",
            name, order
        );
    }
}

/// HiGHS: refblas.c::refblas_start (lines 48-51)
fn refblas_start(n: blasint, inc: blasint) -> blasint {
    if inc > 0 { 0 } else { -(n - 1) * inc }
}

/// HiGHS: refblas.c `A_(i, j)` / `B_(i, j)` / `C_(i, j)` macros (lines 36-38):
/// `(size_t)(i) + (size_t)(j) * (size_t)ld`, computed in `usize` as C does.
#[inline]
fn ij(i: blasint, j: blasint, ld: blasint) -> usize {
    i as usize + j as usize * ld as usize
}

/// rusolve hot-locals: block of `MI x NJ` dot products of the netlib
/// `C := alpha*A**T*B + beta*C` loop, `temp(i,j) = sum_l a(l,i) * b(l,j)` accumulated in
/// increasing `l` from 0.0. The netlib loop runs one such chain at a time and is bound by the
/// floating-point add latency; here `MI * NJ` independent chains advance together. Every
/// chain performs exactly the same additions in the same order, so each `temp` and hence
/// every `c` entry is bit-identical (no reassociation, no FMA).
#[inline(always)]
fn dgemm_tn_block<const MI: usize, const NJ: usize>(
    alpha: f64,
    beta: f64,
    a: &[f64],
    lda: blasint,
    b: &[f64],
    ldb: blasint,
    c: &mut [f64],
    ldc: blasint,
    k: usize,
    i0: blasint,
    j0: blasint,
) {
    let aa: [&[f64]; MI] = std::array::from_fn(|x| &a[ij(0, i0 + x as blasint, lda)..][..k]);
    let bb: [&[f64]; NJ] = std::array::from_fn(|y| &b[ij(0, j0 + y as blasint, ldb)..][..k]);
    let mut t = [[0.0f64; MI]; NJ];
    for l in 0..k {
        for y in 0..NJ {
            let bv = bb[y][l];
            for x in 0..MI {
                t[y][x] = t[y][x] + aa[x][l] * bv;
            }
        }
    }
    for y in 0..NJ {
        for x in 0..MI {
            let idx = ij(i0 + x as blasint, j0 + y as blasint, ldc);
            if beta == 0.0 {
                c[idx] = alpha * t[y][x];
            } else {
                c[idx] = alpha * t[y][x] + beta * c[idx];
            }
        }
    }
}

/// rusolve hot-locals: `NJ` columns of the netlib dtrsm
/// Left / Upper / Trans loop together. Column `j` of `B := alpha*inv(A**T)*B` is solved by its
/// own sequential chain, `temp = temp - a(k,i) * b(k,j)` for increasing `k`; the columns do not
/// depend on each other, so advancing `NJ` chains in lockstep hides the add latency and leaves
/// every operation and its order, hence every result, bit-identical.
#[inline(always)]
fn dtrsm_lut_block<const NJ: usize>(
    nounit: bool,
    alpha: f64,
    a: &[f64],
    lda: blasint,
    b: &mut [f64],
    ldb: blasint,
    m: blasint,
    j0: blasint,
) {
    for i in 0..m {
        let acol = &a[ij(0, i, lda)..][..i as usize];
        let mut t = [0.0f64; NJ];
        for y in 0..NJ {
            t[y] = alpha * b[ij(i, j0 + y as blasint, ldb)];
        }
        for (k, &av) in acol.iter().enumerate() {
            for y in 0..NJ {
                t[y] = t[y] - av * b[ij(k as blasint, j0 + y as blasint, ldb)];
            }
        }
        for y in 0..NJ {
            if nounit {
                t[y] = t[y] / a[ij(i, i, lda)];
            }
            b[ij(i, j0 + y as blasint, ldb)] = t[y];
        }
    }
}

/// HiGHS: refblas.c::cblas_daxpy (netlib daxpy order; mycblas.h lines 22-23)
pub fn cblas_daxpy(n: blasint, alpha: f64, x: &[f64], incx: blasint, y: &mut [f64], incy: blasint) {
    // HiGHS: refblas.c lines 55-67
    if n <= 0 {
        return;
    }
    if alpha == 0.0 {
        return;
    }
    let mut ix = refblas_start(n, incx);
    let mut iy = refblas_start(n, incy);
    for _i in 0..n {
        y[iy as usize] = y[iy as usize] + alpha * x[ix as usize];
        ix += incx;
        iy += incy;
    }
}

/// HiGHS: refblas.c::cblas_dcopy (netlib dcopy order; mycblas.h lines 24-25)
pub fn cblas_dcopy(n: blasint, x: &[f64], incx: blasint, y: &mut [f64], incy: blasint) {
    // HiGHS: refblas.c lines 69-80
    if n <= 0 {
        return;
    }
    let mut ix = refblas_start(n, incx);
    let mut iy = refblas_start(n, incy);
    for _i in 0..n {
        y[iy as usize] = x[ix as usize];
        ix += incx;
        iy += incy;
    }
}

/// HiGHS: refblas.c::cblas_dscal (netlib dscal order; mycblas.h lines 26-27)
pub fn cblas_dscal(n: blasint, alpha: f64, x: &mut [f64], incx: blasint) {
    // HiGHS: refblas.c lines 82-91
    if n <= 0 || incx <= 0 || alpha == 1.0 {
        return;
    }
    let mut ix: blasint = 0;
    for _i in 0..n {
        x[ix as usize] = alpha * x[ix as usize];
        ix += incx;
    }
}

/// HiGHS: refblas.c::cblas_dswap (netlib dswap order; mycblas.h lines 28-29)
/// `x` and `y` never alias in HiPO (`swapCols` swaps two distinct columns
/// of one array: the caller splits it with `split_at_mut`).
pub fn cblas_dswap(n: blasint, x: &mut [f64], incx: blasint, y: &mut [f64], incy: blasint) {
    // HiGHS: refblas.c lines 93-107
    if n <= 0 {
        return;
    }
    let mut ix = refblas_start(n, incx);
    let mut iy = refblas_start(n, incy);
    for _i in 0..n {
        std::mem::swap(&mut x[ix as usize], &mut y[iy as usize]);
        ix += incx;
        iy += incy;
    }
}

/// HiGHS: refblas.c::cblas_dgemv (netlib dgemv order; mycblas.h lines 30-34)
pub fn cblas_dgemv(
    order: i32,
    transa: i32,
    m: blasint,
    n: blasint,
    alpha: f64,
    a: &[f64],
    lda: blasint,
    x: &[f64],
    incx: blasint,
    beta: f64,
    y: &mut [f64],
    incy: blasint,
) {
    // HiGHS: refblas.c lines 111-172
    // A_(i, j) = a[(size_t)(i) + (size_t)(j) * (size_t)lda]
    refblas_col_major_only("cblas_dgemv", order);
    if m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0) {
        return;
    }
    let lenx: blasint;
    let leny: blasint;
    if transa == CblasNoTrans {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }
    let kx = refblas_start(lenx, incx);
    let ky = refblas_start(leny, incy);
    // y := beta*y
    if beta != 1.0 {
        let mut iy = ky;
        if beta == 0.0 {
            for _i in 0..leny {
                y[iy as usize] = 0.0;
                iy += incy;
            }
        } else {
            for _i in 0..leny {
                y[iy as usize] = beta * y[iy as usize];
                iy += incy;
            }
        }
    }
    if alpha == 0.0 {
        return;
    }
    if transa == CblasNoTrans {
        // y := alpha*A*x + y
        let mut jx = kx;
        for j in 0..n {
            let temp = alpha * x[jx as usize];
            let mut iy = ky;
            for i in 0..m {
                y[iy as usize] = y[iy as usize] + temp * a[ij(i, j, lda)];
                iy += incy;
            }
            jx += incx;
        }
    } else {
        // y := alpha*A**T*x + y
        let mut jy = ky;
        for j in 0..n {
            let mut temp = 0.0;
            let mut ix = kx;
            for i in 0..m {
                temp = temp + a[ij(i, j, lda)] * x[ix as usize];
                ix += incx;
            }
            y[jy as usize] = y[jy as usize] + alpha * temp;
            jy += incy;
        }
    }
}

/// HiGHS: refblas.c::cblas_dtpsv (netlib dtpsv order; mycblas.h lines 35-38)
pub fn cblas_dtpsv(
    order: i32,
    uplo: i32,
    transa: i32,
    diag: i32,
    n: blasint,
    ap: &[f64],
    x: &mut [f64],
    incx: blasint,
) {
    // HiGHS: refblas.c lines 174-256
    let nounit = diag == CblasNonUnit;
    refblas_col_major_only("cblas_dtpsv", order);
    if n == 0 {
        return;
    }
    let kx = refblas_start(n, incx);
    // kk is the 0-based packed index of the current diagonal element
    if transa == CblasNoTrans {
        // x := inv(A)*x
        if uplo == CblasUpper {
            let mut kk: blasint = (n * (n + 1)) / 2 - 1;
            let mut jx = kx + (n - 1) * incx;
            let mut j = n - 1;
            while j >= 0 {
                if x[jx as usize] != 0.0 {
                    if nounit {
                        x[jx as usize] = x[jx as usize] / ap[kk as usize];
                    }
                    let temp = x[jx as usize];
                    let mut ix = jx;
                    let mut k = kk - 1;
                    while k >= kk - j {
                        ix -= incx;
                        x[ix as usize] = x[ix as usize] - temp * ap[k as usize];
                        k -= 1;
                    }
                }
                jx -= incx;
                kk -= j + 1;
                j -= 1;
            }
        } else {
            let mut kk: blasint = 0;
            let mut jx = kx;
            for j in 0..n {
                if x[jx as usize] != 0.0 {
                    if nounit {
                        x[jx as usize] = x[jx as usize] / ap[kk as usize];
                    }
                    let temp = x[jx as usize];
                    let mut ix = jx;
                    let mut k = kk + 1;
                    while k <= kk + n - 1 - j {
                        ix += incx;
                        x[ix as usize] = x[ix as usize] - temp * ap[k as usize];
                        k += 1;
                    }
                }
                jx += incx;
                kk += n - j;
            }
        }
    } else {
        // x := inv(A**T)*x
        if uplo == CblasUpper {
            let mut kk: blasint = 0; // packed index of the first element of column j
            let mut jx = kx;
            for j in 0..n {
                let mut temp = x[jx as usize];
                let mut ix = kx;
                let mut k = kk;
                while k < kk + j {
                    temp = temp - ap[k as usize] * x[ix as usize];
                    ix += incx;
                    k += 1;
                }
                if nounit {
                    temp = temp / ap[(kk + j) as usize];
                }
                x[jx as usize] = temp;
                jx += incx;
                kk += j + 1;
            }
        } else {
            let mut kk: blasint = (n * (n + 1)) / 2 - 1; // packed index of the last element of col j
            let mut kx = kx + (n - 1) * incx;
            let mut jx = kx;
            let mut j = n - 1;
            while j >= 0 {
                let mut temp = x[jx as usize];
                let mut ix = kx;
                let mut k = kk;
                while k > kk - (n - 1 - j) {
                    temp = temp - ap[k as usize] * x[ix as usize];
                    ix -= incx;
                    k -= 1;
                }
                if nounit {
                    temp = temp / ap[(kk - (n - 1 - j)) as usize];
                }
                x[jx as usize] = temp;
                jx -= incx;
                kk -= n - j;
                j -= 1;
            }
        }
    }
}

/// HiGHS: refblas.c::cblas_dtrsv (netlib dtrsv order; mycblas.h lines 39-42)
pub fn cblas_dtrsv(
    order: i32,
    uplo: i32,
    transa: i32,
    diag: i32,
    n: blasint,
    a: &[f64],
    lda: blasint,
    x: &mut [f64],
    incx: blasint,
) {
    // HiGHS: refblas.c lines 258-330
    // A_(i, j) = a[(size_t)(i) + (size_t)(j) * (size_t)lda]
    let nounit = diag == CblasNonUnit;
    refblas_col_major_only("cblas_dtrsv", order);
    if n == 0 {
        return;
    }
    let kx = refblas_start(n, incx);
    if transa == CblasNoTrans {
        // x := inv(A)*x
        if uplo == CblasUpper {
            let mut jx = kx + (n - 1) * incx;
            let mut j = n - 1;
            while j >= 0 {
                if x[jx as usize] != 0.0 {
                    if nounit {
                        x[jx as usize] = x[jx as usize] / a[ij(j, j, lda)];
                    }
                    let temp = x[jx as usize];
                    let mut ix = jx;
                    let mut i = j - 1;
                    while i >= 0 {
                        ix -= incx;
                        x[ix as usize] = x[ix as usize] - temp * a[ij(i, j, lda)];
                        i -= 1;
                    }
                }
                jx -= incx;
                j -= 1;
            }
        } else {
            let mut jx = kx;
            for j in 0..n {
                if x[jx as usize] != 0.0 {
                    if nounit {
                        x[jx as usize] = x[jx as usize] / a[ij(j, j, lda)];
                    }
                    let temp = x[jx as usize];
                    let mut ix = jx;
                    for i in (j + 1)..n {
                        ix += incx;
                        x[ix as usize] = x[ix as usize] - temp * a[ij(i, j, lda)];
                    }
                }
                jx += incx;
            }
        }
    } else {
        // x := inv(A**T)*x
        if uplo == CblasUpper {
            let mut jx = kx;
            for j in 0..n {
                let mut temp = x[jx as usize];
                let mut ix = kx;
                for i in 0..j {
                    temp = temp - a[ij(i, j, lda)] * x[ix as usize];
                    ix += incx;
                }
                if nounit {
                    temp = temp / a[ij(j, j, lda)];
                }
                x[jx as usize] = temp;
                jx += incx;
            }
        } else {
            let mut kx = kx + (n - 1) * incx;
            let mut jx = kx;
            let mut j = n - 1;
            while j >= 0 {
                let mut temp = x[jx as usize];
                let mut ix = kx;
                let mut i = n - 1;
                while i > j {
                    temp = temp - a[ij(i, j, lda)] * x[ix as usize];
                    ix -= incx;
                    i -= 1;
                }
                if nounit {
                    temp = temp / a[ij(j, j, lda)];
                }
                x[jx as usize] = temp;
                jx -= incx;
                j -= 1;
            }
        }
    }
}

/// HiGHS: refblas.c::cblas_dger (netlib dger order; mycblas.h lines 43-46)
pub fn cblas_dger(
    order: i32,
    m: blasint,
    n: blasint,
    alpha: f64,
    x: &[f64],
    incx: blasint,
    y: &[f64],
    incy: blasint,
    a: &mut [f64],
    lda: blasint,
) {
    // HiGHS: refblas.c lines 332-353
    // A_(i, j) = a[(size_t)(i) + (size_t)(j) * (size_t)lda]
    refblas_col_major_only("cblas_dger", order);
    if m == 0 || n == 0 || alpha == 0.0 {
        return;
    }
    let mut jy = refblas_start(n, incy);
    let kx = refblas_start(m, incx);
    for j in 0..n {
        if y[jy as usize] != 0.0 {
            let temp = alpha * y[jy as usize];
            let mut ix = kx;
            for i in 0..m {
                a[ij(i, j, lda)] = a[ij(i, j, lda)] + x[ix as usize] * temp;
                ix += incx;
            }
        }
        jy += incy;
    }
}

/// HiGHS: refblas.c::cblas_dgemm (netlib dgemm order; mycblas.h lines 47-53)
pub fn cblas_dgemm(
    order: i32,
    transa: i32,
    transb: i32,
    m: blasint,
    n: blasint,
    k: blasint,
    alpha: f64,
    a: &[f64],
    lda: blasint,
    b: &[f64],
    ldb: blasint,
    beta: f64,
    c: &mut [f64],
    ldc: blasint,
) {
    // HiGHS: refblas.c lines 357-435
    // A_(i, j) = a[(size_t)(i) + (size_t)(j) * (size_t)lda]
    // B_(i, j) = b[(size_t)(i) + (size_t)(j) * (size_t)ldb]
    // C_(i, j) = c[(size_t)(i) + (size_t)(j) * (size_t)ldc]
    let nota = transa == CblasNoTrans;
    let notb = transb == CblasNoTrans;
    refblas_col_major_only("cblas_dgemm", order);
    if m == 0 || n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0) {
        return;
    }
    if alpha == 0.0 {
        if beta == 0.0 {
            for j in 0..n {
                for i in 0..m {
                    c[ij(i, j, ldc)] = 0.0;
                }
            }
        } else {
            for j in 0..n {
                for i in 0..m {
                    c[ij(i, j, ldc)] = beta * c[ij(i, j, ldc)];
                }
            }
        }
        return;
    }
    if notb {
        if nota {
            // C := alpha*A*B + beta*C
            for j in 0..n {
                if beta == 0.0 {
                    for i in 0..m {
                        c[ij(i, j, ldc)] = 0.0;
                    }
                } else if beta != 1.0 {
                    for i in 0..m {
                        c[ij(i, j, ldc)] = beta * c[ij(i, j, ldc)];
                    }
                }
                for l in 0..k {
                    let temp = alpha * b[ij(l, j, ldb)];
                    for i in 0..m {
                        c[ij(i, j, ldc)] = c[ij(i, j, ldc)] + temp * a[ij(i, l, lda)];
                    }
                }
            }
        } else {
            // C := alpha*A**T*B + beta*C
            let ku = k as usize;
            let mut j = 0;
            while j < n {
                let nj2 = j + 2 <= n;
                let mut i = 0;
                while i < m {
                    let mi4 = i + 4 <= m;
                    match (mi4, nj2) {
                        (true, true) => {
                            dgemm_tn_block::<4, 2>(alpha, beta, a, lda, b, ldb, c, ldc, ku, i, j)
                        }
                        (true, false) => {
                            dgemm_tn_block::<4, 1>(alpha, beta, a, lda, b, ldb, c, ldc, ku, i, j)
                        }
                        (false, true) => {
                            dgemm_tn_block::<1, 2>(alpha, beta, a, lda, b, ldb, c, ldc, ku, i, j)
                        }
                        (false, false) => {
                            dgemm_tn_block::<1, 1>(alpha, beta, a, lda, b, ldb, c, ldc, ku, i, j)
                        }
                    }
                    i += if mi4 { 4 } else { 1 };
                }
                j += if nj2 { 2 } else { 1 };
            }
        }
    } else {
        if nota {
            // C := alpha*A*B**T + beta*C
            for j in 0..n {
                if beta == 0.0 {
                    for i in 0..m {
                        c[ij(i, j, ldc)] = 0.0;
                    }
                } else if beta != 1.0 {
                    for i in 0..m {
                        c[ij(i, j, ldc)] = beta * c[ij(i, j, ldc)];
                    }
                }
                for l in 0..k {
                    let temp = alpha * b[ij(j, l, ldb)];
                    for i in 0..m {
                        c[ij(i, j, ldc)] = c[ij(i, j, ldc)] + temp * a[ij(i, l, lda)];
                    }
                }
            }
        } else {
            // C := alpha*A**T*B**T + beta*C
            for j in 0..n {
                for i in 0..m {
                    let mut temp = 0.0;
                    for l in 0..k {
                        temp = temp + a[ij(l, i, lda)] * b[ij(j, l, ldb)];
                    }
                    if beta == 0.0 {
                        c[ij(i, j, ldc)] = alpha * temp;
                    } else {
                        c[ij(i, j, ldc)] = alpha * temp + beta * c[ij(i, j, ldc)];
                    }
                }
            }
        }
    }
}

/// HiGHS: refblas.c::cblas_dsyrk (netlib dsyrk order; mycblas.h lines 54-58)
pub fn cblas_dsyrk(
    order: i32,
    uplo: i32,
    trans: i32,
    n: blasint,
    k: blasint,
    alpha: f64,
    a: &[f64],
    lda: blasint,
    beta: f64,
    c: &mut [f64],
    ldc: blasint,
) {
    // HiGHS: refblas.c lines 437-500
    // A_(i, j) = a[(size_t)(i) + (size_t)(j) * (size_t)lda]
    // C_(i, j) = c[(size_t)(i) + (size_t)(j) * (size_t)ldc]
    let upper = uplo == CblasUpper;
    refblas_col_major_only("cblas_dsyrk", order);
    if n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0) {
        return;
    }
    // the stored triangle of column j is rows [i0, i1)
    // REFBLAS_TRI(j): if upper { i0 = 0; i1 = j + 1; } else { i0 = j; i1 = n; }
    if alpha == 0.0 {
        for j in 0..n {
            let (i0, i1) = if upper { (0, j + 1) } else { (j, n) };
            if beta == 0.0 {
                for i in i0..i1 {
                    c[ij(i, j, ldc)] = 0.0;
                }
            } else {
                for i in i0..i1 {
                    c[ij(i, j, ldc)] = beta * c[ij(i, j, ldc)];
                }
            }
        }
        return;
    }
    if trans == CblasNoTrans {
        // C := alpha*A*A**T + beta*C
        for j in 0..n {
            let (i0, i1) = if upper { (0, j + 1) } else { (j, n) };
            if beta == 0.0 {
                for i in i0..i1 {
                    c[ij(i, j, ldc)] = 0.0;
                }
            } else if beta != 1.0 {
                for i in i0..i1 {
                    c[ij(i, j, ldc)] = beta * c[ij(i, j, ldc)];
                }
            }
            for l in 0..k {
                if a[ij(j, l, lda)] != 0.0 {
                    let temp = alpha * a[ij(j, l, lda)];
                    for i in i0..i1 {
                        c[ij(i, j, ldc)] = c[ij(i, j, ldc)] + temp * a[ij(i, l, lda)];
                    }
                }
            }
        }
    } else {
        // C := alpha*A**T*A + beta*C
        for j in 0..n {
            let (i0, i1) = if upper { (0, j + 1) } else { (j, n) };
            for i in i0..i1 {
                let mut temp = 0.0;
                for l in 0..k {
                    temp = temp + a[ij(l, i, lda)] * a[ij(l, j, lda)];
                }
                if beta == 0.0 {
                    c[ij(i, j, ldc)] = alpha * temp;
                } else {
                    c[ij(i, j, ldc)] = alpha * temp + beta * c[ij(i, j, ldc)];
                }
            }
        }
    }
}

/// HiGHS: refblas.c::cblas_dtrsm (netlib dtrsm order; mycblas.h lines 59-63)
pub fn cblas_dtrsm(
    order: i32,
    side: i32,
    uplo: i32,
    transa: i32,
    diag: i32,
    m: blasint,
    n: blasint,
    alpha: f64,
    a: &[f64],
    lda: blasint,
    b: &mut [f64],
    ldb: blasint,
) {
    // HiGHS: refblas.c lines 502-635
    // A_(i, j) = a[(size_t)(i) + (size_t)(j) * (size_t)lda]
    // B_(i, j) = b[(size_t)(i) + (size_t)(j) * (size_t)ldb]
    let lside = side == CblasLeft;
    let upper = uplo == CblasUpper;
    let nounit = diag == CblasNonUnit;
    refblas_col_major_only("cblas_dtrsm", order);
    if m == 0 || n == 0 {
        return;
    }
    if alpha == 0.0 {
        for j in 0..n {
            for i in 0..m {
                b[ij(i, j, ldb)] = 0.0;
            }
        }
        return;
    }
    if lside {
        if transa == CblasNoTrans {
            // B := alpha*inv(A)*B
            if upper {
                for j in 0..n {
                    if alpha != 1.0 {
                        for i in 0..m {
                            b[ij(i, j, ldb)] = alpha * b[ij(i, j, ldb)];
                        }
                    }
                    let mut k = m - 1;
                    while k >= 0 {
                        if b[ij(k, j, ldb)] != 0.0 {
                            if nounit {
                                b[ij(k, j, ldb)] = b[ij(k, j, ldb)] / a[ij(k, k, lda)];
                            }
                            for i in 0..k {
                                b[ij(i, j, ldb)] =
                                    b[ij(i, j, ldb)] - b[ij(k, j, ldb)] * a[ij(i, k, lda)];
                            }
                        }
                        k -= 1;
                    }
                }
            } else {
                for j in 0..n {
                    if alpha != 1.0 {
                        for i in 0..m {
                            b[ij(i, j, ldb)] = alpha * b[ij(i, j, ldb)];
                        }
                    }
                    for k in 0..m {
                        if b[ij(k, j, ldb)] != 0.0 {
                            if nounit {
                                b[ij(k, j, ldb)] = b[ij(k, j, ldb)] / a[ij(k, k, lda)];
                            }
                            for i in (k + 1)..m {
                                b[ij(i, j, ldb)] =
                                    b[ij(i, j, ldb)] - b[ij(k, j, ldb)] * a[ij(i, k, lda)];
                            }
                        }
                    }
                }
            }
        } else {
            // B := alpha*inv(A**T)*B
            if upper {
                let mut j = 0;
                while j < n {
                    if j + 4 <= n {
                        dtrsm_lut_block::<4>(nounit, alpha, a, lda, b, ldb, m, j);
                        j += 4;
                    } else {
                        dtrsm_lut_block::<1>(nounit, alpha, a, lda, b, ldb, m, j);
                        j += 1;
                    }
                }
            } else {
                for j in 0..n {
                    let mut i = m - 1;
                    while i >= 0 {
                        let mut temp = alpha * b[ij(i, j, ldb)];
                        for k in (i + 1)..m {
                            temp = temp - a[ij(k, i, lda)] * b[ij(k, j, ldb)];
                        }
                        if nounit {
                            temp = temp / a[ij(i, i, lda)];
                        }
                        b[ij(i, j, ldb)] = temp;
                        i -= 1;
                    }
                }
            }
        }
    } else {
        if transa == CblasNoTrans {
            // B := alpha*B*inv(A)
            if upper {
                for j in 0..n {
                    if alpha != 1.0 {
                        for i in 0..m {
                            b[ij(i, j, ldb)] = alpha * b[ij(i, j, ldb)];
                        }
                    }
                    for k in 0..j {
                        if a[ij(k, j, lda)] != 0.0 {
                            for i in 0..m {
                                b[ij(i, j, ldb)] =
                                    b[ij(i, j, ldb)] - a[ij(k, j, lda)] * b[ij(i, k, ldb)];
                            }
                        }
                    }
                    if nounit {
                        let temp = 1.0 / a[ij(j, j, lda)];
                        for i in 0..m {
                            b[ij(i, j, ldb)] = temp * b[ij(i, j, ldb)];
                        }
                    }
                }
            } else {
                let mut j = n - 1;
                while j >= 0 {
                    if alpha != 1.0 {
                        for i in 0..m {
                            b[ij(i, j, ldb)] = alpha * b[ij(i, j, ldb)];
                        }
                    }
                    for k in (j + 1)..n {
                        if a[ij(k, j, lda)] != 0.0 {
                            for i in 0..m {
                                b[ij(i, j, ldb)] =
                                    b[ij(i, j, ldb)] - a[ij(k, j, lda)] * b[ij(i, k, ldb)];
                            }
                        }
                    }
                    if nounit {
                        let temp = 1.0 / a[ij(j, j, lda)];
                        for i in 0..m {
                            b[ij(i, j, ldb)] = temp * b[ij(i, j, ldb)];
                        }
                    }
                    j -= 1;
                }
            }
        } else {
            // B := alpha*B*inv(A**T)
            if upper {
                let mut k = n - 1;
                while k >= 0 {
                    if nounit {
                        let temp = 1.0 / a[ij(k, k, lda)];
                        for i in 0..m {
                            b[ij(i, k, ldb)] = temp * b[ij(i, k, ldb)];
                        }
                    }
                    for j in 0..k {
                        if a[ij(j, k, lda)] != 0.0 {
                            let temp = a[ij(j, k, lda)];
                            for i in 0..m {
                                b[ij(i, j, ldb)] = b[ij(i, j, ldb)] - temp * b[ij(i, k, ldb)];
                            }
                        }
                    }
                    if alpha != 1.0 {
                        for i in 0..m {
                            b[ij(i, k, ldb)] = alpha * b[ij(i, k, ldb)];
                        }
                    }
                    k -= 1;
                }
            } else {
                for k in 0..n {
                    if nounit {
                        let temp = 1.0 / a[ij(k, k, lda)];
                        for i in 0..m {
                            b[ij(i, k, ldb)] = temp * b[ij(i, k, ldb)];
                        }
                    }
                    for j in (k + 1)..n {
                        if a[ij(j, k, lda)] != 0.0 {
                            let temp = a[ij(j, k, lda)];
                            for i in 0..m {
                                b[ij(i, j, ldb)] = b[ij(i, j, ldb)] - temp * b[ij(i, k, ldb)];
                            }
                        }
                    }
                    if alpha != 1.0 {
                        for i in 0..m {
                            b[ij(i, k, ldb)] = alpha * b[ij(i, k, ldb)];
                        }
                    }
                }
            }
        }
    }
}

/// HiGHS: myblas.cpp::highs_openblas_set_num_threads (lines 3-7). A no-op:
/// the reference BLAS is single threaded.
pub fn highs_openblas_set_num_threads(num_threads: i32) {}

// ---------------------------------------------------------------------------
// rusolve lever `ipm_hipo_fast_blas` (default off): a packed, register-tiled
// dgemm for the one shape the hybrid dense factorisation spends its flops on,
// C := alpha*A^T*B + beta*C (column-major). The reference loop above is a
// strictly ordered dot product that the compiler cannot vectorise (about 3
// GFlop/s); this follows the packing + micro-kernel structure of Goto and van
// de Geijn, "Anatomy of high-performance matrix multiplication", ACM TOMS 34(3),
// 2008. The summation order differs from the reference, so the result is not
// bit-identical to the oracle; the lever is therefore off by default.
// ---------------------------------------------------------------------------

use std::sync::atomic::{AtomicBool, Ordering};

/// Set from `HighsOptions::ipm_hipo_fast_blas` when a HiPO solve starts.
pub static FAST_BLAS: AtomicBool = AtomicBool::new(false);

/// Set from `HighsOptions::ipm_hipo_cache_perm` when a HiPO solve starts.
pub static CACHE_PERM: AtomicBool = AtomicBool::new(false);

/// Set from `HighsOptions::ipm_hipo_par_ne` when a HiPO solve starts.
pub static PAR_NE: AtomicBool = AtomicBool::new(false);
/// rusolve `ipm_hipo_dense_cols`: process-wide for the solve
/// rusolve (dualize guard): predicted factor flops of the primal while a dualized LP is solved;
/// f64 bits, 0 = no guard. The dual's solve stops after its own analysis when it is clearly
/// costlier, so the chosen side is analysed once.
pub static GUARD_PRIMAL_FLOPS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static DENSE_COLS: AtomicBool = AtomicBool::new(false);
/// rusolve `ipm_hipo_dense_schur`: Woodbury preconditioner for the dense-column CG.
pub static DENSE_SCHUR: AtomicBool = AtomicBool::new(false);

const MR: usize = 8;
const NR: usize = 6;

#[inline(always)]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    #[cfg(target_arch = "aarch64")]
    {
        a.mul_add(b, c)
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        a * b + c
    }
}

type Acc = [[f64; MR]; NR];

/// Portable micro-kernel: `acc[jj][ii] += sum_l ap[l*MR+ii] * bp[l*NR+jj]`.
#[inline(always)]
fn micro_generic(ap: &[f64], bp: &[f64], k: usize, acc: &mut Acc) {
    for l in 0..k {
        let av = &ap[l * MR..l * MR + MR];
        let bv = &bp[l * NR..l * NR + NR];
        for jj in 0..NR {
            let bj = bv[jj];
            for ii in 0..MR {
                acc[jj][ii] = fma(av[ii], bj, acc[jj][ii]);
            }
        }
    }
}

/// AVX2 + FMA micro-kernel, 8 x 6 tile in 12 ymm accumulators (Goto and van de Geijn 2008,
/// section 4; the BLIS `haswell` dgemm kernel has this shape).
///
/// SAFETY: callers must have checked `is_x86_feature_detected!("avx2")` and `"fma"`. The
/// asserts below guarantee every unaligned load stays inside `ap` (k * MR) and `bp`
/// (k * NR); `acc` is a plain `[[f64; 8]; 6]` and is written with unaligned stores.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn micro_avx2(ap: &[f64], bp: &[f64], k: usize, acc: &mut Acc) {
    use std::arch::x86_64::*;
    assert!(ap.len() >= k * MR && bp.len() >= k * NR);
    let mut c = [[_mm256_setzero_pd(); 2]; NR];
    let ap = ap.as_ptr();
    let bp = bp.as_ptr();
    for l in 0..k {
        let a0 = unsafe { _mm256_loadu_pd(ap.add(l * MR)) };
        let a1 = unsafe { _mm256_loadu_pd(ap.add(l * MR + 4)) };
        for jj in 0..NR {
            let b = unsafe { _mm256_broadcast_sd(&*bp.add(l * NR + jj)) };
            c[jj][0] = _mm256_fmadd_pd(a0, b, c[jj][0]);
            c[jj][1] = _mm256_fmadd_pd(a1, b, c[jj][1]);
        }
    }
    for jj in 0..NR {
        unsafe {
            _mm256_storeu_pd(acc[jj].as_mut_ptr(), c[jj][0]);
            _mm256_storeu_pd(acc[jj].as_mut_ptr().add(4), c[jj][1]);
        }
    }
}

/// NEON 8 x 6 tile: 24 q-register accumulators (4 per column, 6 columns) plus 4 A vectors and
/// 3 B vectors, 31 of the 32 vector registers (Goto and van de Geijn 2008, section 4; BLIS
/// `armv8a` dgemm kernels use the same scheme).
///
/// SAFETY: NEON is mandatory on aarch64, so no feature detection is needed. The assert
/// guarantees every load stays inside `ap` (k * MR) and `bp` (k * NR); the stores target the
/// plain `[[f64; 8]; 6]` accumulator with in-bounds offsets.
#[cfg(target_arch = "aarch64")]
unsafe fn micro_neon(ap: &[f64], bp: &[f64], k: usize, acc: &mut Acc) {
    use std::arch::aarch64::*;
    assert!(ap.len() >= k * MR && bp.len() >= k * NR);
    unsafe {
        let mut c = [[vdupq_n_f64(0.0); 4]; NR];
        let ap = ap.as_ptr();
        let bp = bp.as_ptr();
        for l in 0..k {
            let a = [
                vld1q_f64(ap.add(l * MR)),
                vld1q_f64(ap.add(l * MR + 2)),
                vld1q_f64(ap.add(l * MR + 4)),
                vld1q_f64(ap.add(l * MR + 6)),
            ];
            let b0 = vld1q_f64(bp.add(l * NR));
            let b1 = vld1q_f64(bp.add(l * NR + 2));
            let b2 = vld1q_f64(bp.add(l * NR + 4));
            for q in 0..4 {
                c[0][q] = vfmaq_laneq_f64::<0>(c[0][q], a[q], b0);
                c[1][q] = vfmaq_laneq_f64::<1>(c[1][q], a[q], b0);
                c[2][q] = vfmaq_laneq_f64::<0>(c[2][q], a[q], b1);
                c[3][q] = vfmaq_laneq_f64::<1>(c[3][q], a[q], b1);
                c[4][q] = vfmaq_laneq_f64::<0>(c[4][q], a[q], b2);
                c[5][q] = vfmaq_laneq_f64::<1>(c[5][q], a[q], b2);
            }
        }
        for jj in 0..NR {
            for q in 0..4 {
                vst1q_f64(acc[jj].as_mut_ptr().add(2 * q), c[jj][q]);
            }
        }
    }
}

#[inline]
fn have_avx2_fma() -> bool {
    #[cfg(target_arch = "aarch64")]
    {
        true
    }
    #[cfg(target_arch = "x86_64")]
    {
        is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma")
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        false
    }
}

#[inline(always)]
fn micro(avx: bool, ap: &[f64], bp: &[f64], k: usize, acc: &mut Acc) {
    #[cfg(target_arch = "x86_64")]
    if avx {
        // SAFETY: `avx` is only true when `have_avx2_fma()` returned true.
        unsafe { micro_avx2(ap, bp, k, acc) };
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if avx {
        // SAFETY: NEON is always present on aarch64; bounds are asserted inside.
        unsafe { micro_neon(ap, bp, k, acc) };
        return;
    }
    let _ = avx;
    micro_generic(ap, bp, k, acc);
}

thread_local! {
    static SCRATCH: std::cell::RefCell<Vec<f64>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Pack `A` (k x m, lda) into MR-row panels, layout `[panel][l][ii]`, zero-padded.
fn pack_a(ap: &mut Vec<f64>, m: usize, k: usize, a: &[f64], lda: usize) {
    let mp = m.div_ceil(MR) * MR;
    ap.clear();
    ap.resize(mp * k, 0.0);
    for p in 0..mp / MR {
        let base = p * MR * k;
        let rows = (m - p * MR).min(MR);
        for ii in 0..rows {
            let col = &a[(p * MR + ii) * lda..(p * MR + ii) * lda + k];
            for (l, &v) in col.iter().enumerate() {
                ap[base + l * MR + ii] = v;
            }
        }
    }
}

/// `C(m x n) := alpha * A^T B + beta * C` where `B` (k x n, `ldb`) and `C` (m x n, `ldc`) live
/// in the same buffer `bc` at element offsets `b_off` and `c_off` (disjoint row ranges). Each B
/// panel is packed before any C entry of that panel is written, so the overlap is harmless.
fn gemm_tn_core(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    bc: &mut [f64],
    b_off: usize,
    ldb: usize,
    beta: f64,
    c_off: usize,
    ldc: usize,
) {
    if m == 0 || n == 0 {
        return;
    }
    let avx = have_avx2_fma();
    let mp = m.div_ceil(MR) * MR;
    let mut ap = SCRATCH.with(|s| std::mem::take(&mut *s.borrow_mut()));
    pack_a(&mut ap, m, k, a, lda);
    let mut bp = vec![0.0f64; NR * k];
    let mut j = 0;
    while j < n {
        let cols = (n - j).min(NR);
        for l in 0..k {
            for jj in 0..NR {
                bp[l * NR + jj] = if jj < cols {
                    bc[b_off + (j + jj) * ldb + l]
                } else {
                    0.0
                };
            }
        }
        for p in 0..mp / MR {
            let apan = &ap[p * MR * k..(p + 1) * MR * k];
            let mut acc = [[0.0f64; MR]; NR];
            micro(avx, apan, &bp, k, &mut acc);
            let rows = (m - p * MR).min(MR);
            for jj in 0..cols {
                for ii in 0..rows {
                    let idx = c_off + (j + jj) * ldc + p * MR + ii;
                    let v = alpha * acc[jj][ii];
                    bc[idx] = if beta == 0.0 { v } else { v + beta * bc[idx] };
                }
            }
        }
        j += NR;
    }
    SCRATCH.with(|s| *s.borrow_mut() = ap);
}

/// C(m x n) := alpha * A^T B + beta * C with A (k x m, lda) and B (k x n, ldb).
pub fn fast_dgemm_tn(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    beta: f64,
    c: &mut [f64],
    ldc: usize,
) {
    if m == 0 || n == 0 {
        return;
    }
    let avx = have_avx2_fma();
    let mp = m.div_ceil(MR) * MR;
    let mut ap = SCRATCH.with(|s| std::mem::take(&mut *s.borrow_mut()));
    pack_a(&mut ap, m, k, a, lda);
    let mut bp = vec![0.0f64; NR * k];
    let mut j = 0;
    while j < n {
        let cols = (n - j).min(NR);
        for l in 0..k {
            for jj in 0..NR {
                bp[l * NR + jj] = if jj < cols {
                    b[(j + jj) * ldb + l]
                } else {
                    0.0
                };
            }
        }
        for p in 0..mp / MR {
            let apan = &ap[p * MR * k..(p + 1) * MR * k];
            let mut acc = [[0.0f64; MR]; NR];
            micro(avx, apan, &bp, k, &mut acc);
            let rows = (m - p * MR).min(MR);
            for jj in 0..cols {
                for ii in 0..rows {
                    let idx = (j + jj) * ldc + p * MR + ii;
                    let v = alpha * acc[jj][ii];
                    c[idx] = if beta == 0.0 { v } else { v + beta * c[idx] };
                }
            }
        }
        j += NR;
    }
    SCRATCH.with(|s| *s.borrow_mut() = ap);
}

/// B := inv(U^T) * B for upper unit-diagonal U (m x m, lda), B (m x n, ldb): the
/// `dtrsm('L','U','T','U')` call of the hybrid dense factorisation. Blocked: for each block of
/// `TB` rows, the contribution of all earlier rows is one packed `fast_dgemm_tn` update
/// (`B_I -= U_{L,I}^T B_L`), then the small diagonal block is solved with eight independent
/// column chains. Most flops therefore run in the micro-kernel (Goto and van de Geijn 2008;
/// the blocked triangular solve of the BLIS `trsm` is built the same way).
pub fn fast_dtrsm_lutu(m: usize, n: usize, a: &[f64], lda: usize, b: &mut [f64], ldb: usize) {
    const TB: usize = 32;
    const W: usize = 8;
    let mut i0 = 0;
    while i0 < m {
        let ib = (m - i0).min(TB);
        if i0 > 0 {
            // B[i0..i0+ib, :] -= U[0..i0, i0..i0+ib]^T * B[0..i0, :]
            gemm_tn_core(
                ib,
                n,
                i0,
                -1.0,
                &a[i0 * lda..],
                lda,
                b,
                0,
                ldb,
                1.0,
                i0,
                ldb,
            );
        }
        let mut j = 0;
        while j < n {
            let w = (n - j).min(W);
            for i in i0..i0 + ib {
                let ucol = &a[i * lda + i0..i * lda + i];
                let mut acc = [0.0f64; W];
                for (l, &u) in ucol.iter().enumerate() {
                    for c in 0..w {
                        acc[c] = fma(u, b[(j + c) * ldb + i0 + l], acc[c]);
                    }
                }
                for c in 0..w {
                    b[(j + c) * ldb + i] -= acc[c];
                }
            }
            j += W;
        }
        i0 += ib;
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn daxpy_unit_and_strided() {
        let x = [1.0, 2.0, 3.0];
        let mut y = [1.0, 1.0, 1.0];
        cblas_daxpy(3, 2.0, &x, 1, &mut y, 1);
        assert_eq!(y, [3.0, 5.0, 7.0]);
        let mut z = [0.0; 4];
        cblas_daxpy(2, 1.0, &x, 2, &mut z, 3);
        assert_eq!(z, [1.0, 0.0, 0.0, 3.0]);
    }

    #[test]
    fn dgemm_small() {
        // C = A * B^T, A 2x2 col-major [1 3; 2 4], B = I
        let a = [1.0, 2.0, 3.0, 4.0];
        let b = [1.0, 0.0, 0.0, 1.0];
        let mut c = [0.0; 4];
        cblas_dgemm(
            CblasColMajor,
            CblasNoTrans,
            CblasTrans,
            2,
            2,
            2,
            1.0,
            &a,
            2,
            &b,
            2,
            0.0,
            &mut c,
            2,
        );
        assert_eq!(c, a);
    }

    /// rusolve hot-locals: the blocked `A**T*B` kernel equals the netlib loop bit for bit,
    /// for every block remainder (m, n not multiples of the block) and both beta paths.
    #[test]
    fn dgemm_tn_blocked_is_bit_identical_to_netlib_loop() {
        let mut seed = 12345u64;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 7.0
        };
        for &(m, n, k) in &[
            (1, 1, 1),
            (3, 5, 7),
            (4, 2, 9),
            (9, 7, 33),
            (13, 6, 5),
            (8, 8, 64),
        ] {
            for &beta in &[0.0, 1.0, -0.75] {
                let (lda, ldb, ldc) = (k + 2, k + 1, m + 3);
                let a: Vec<f64> = (0..lda as usize * m as usize).map(|_| rnd()).collect();
                let b: Vec<f64> = (0..ldb as usize * n as usize).map(|_| rnd()).collect();
                let c0: Vec<f64> = (0..ldc as usize * n as usize).map(|_| rnd()).collect();
                let mut c = c0.clone();
                cblas_dgemm(
                    CblasColMajor,
                    CblasTrans,
                    CblasNoTrans,
                    m,
                    n,
                    k,
                    -1.3,
                    &a,
                    lda,
                    &b,
                    ldb,
                    beta,
                    &mut c,
                    ldc,
                );
                let mut want = c0.clone();
                for j in 0..n {
                    for i in 0..m {
                        let mut temp = 0.0;
                        for l in 0..k {
                            temp = temp + a[ij(l, i, lda)] * b[ij(l, j, ldb)];
                        }
                        want[ij(i, j, ldc)] = if beta == 0.0 {
                            -1.3 * temp
                        } else {
                            -1.3 * temp + beta * want[ij(i, j, ldc)]
                        };
                    }
                }
                for (x, y) in c.iter().zip(&want) {
                    assert_eq!(x.to_bits(), y.to_bits());
                }
            }
        }
    }

    /// rusolve hot-locals: the column-blocked Left/Upper/Trans dtrsm equals the netlib loop
    /// bit for bit, for unit and non-unit diagonals and every column remainder.
    #[test]
    fn dtrsm_lut_blocked_is_bit_identical_to_netlib_loop() {
        let mut seed = 99u64;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 3.0 + 0.1
        };
        for &(m, n) in &[(1, 1), (5, 3), (7, 4), (12, 9), (30, 6)] {
            for &diag in &[CblasUnit, CblasNonUnit] {
                let (lda, ldb) = (m + 2, m + 1);
                let a: Vec<f64> = (0..lda as usize * m as usize).map(|_| rnd()).collect();
                let b0: Vec<f64> = (0..ldb as usize * n as usize).map(|_| rnd()).collect();
                let mut b = b0.clone();
                cblas_dtrsm(
                    CblasColMajor,
                    CblasLeft,
                    CblasUpper,
                    CblasTrans,
                    diag,
                    m,
                    n,
                    0.7,
                    &a,
                    lda,
                    &mut b,
                    ldb,
                );
                let mut want = b0.clone();
                for j in 0..n {
                    for i in 0..m {
                        let mut temp = 0.7 * want[ij(i, j, ldb)];
                        for k in 0..i {
                            temp = temp - a[ij(k, i, lda)] * want[ij(k, j, ldb)];
                        }
                        if diag == CblasNonUnit {
                            temp = temp / a[ij(i, i, lda)];
                        }
                        want[ij(i, j, ldb)] = temp;
                    }
                }
                for (x, y) in b.iter().zip(&want) {
                    assert_eq!(x.to_bits(), y.to_bits());
                }
            }
        }
    }

    #[test]
    fn dtrsv_lower_unit() {
        // L = [1 0; 2 1], solve L x = [1, 4] -> x = [1, 2]
        let a = [1.0, 2.0, 0.0, 1.0];
        let mut x = [1.0, 4.0];
        cblas_dtrsv(
            CblasColMajor,
            CblasLower,
            CblasNoTrans,
            CblasUnit,
            2,
            &a,
            2,
            &mut x,
            1,
        );
        assert_eq!(x, [1.0, 2.0]);
    }
}

#[cfg(test)]
mod fast_tests {
    use super::*;

    fn fill(len: usize, seed: u64) -> Vec<f64> {
        let mut s = seed;
        (0..len)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((s >> 33) as f64 / (1u64 << 31) as f64) - 0.5
            })
            .collect()
    }

    #[test]
    fn fast_dgemm_tn_matches_reference() {
        for &(m, n, k, beta) in &[
            (1, 1, 1, 0.0),
            (7, 5, 3, 1.0),
            (13, 9, 17, 0.0),
            (64, 70, 33, 1.0),
        ] {
            let (lda, ldb, ldc) = (k + 2, k + 1, m + 3);
            let a = fill(lda * m, 1);
            let b = fill(ldb * n, 2);
            let c0 = fill(ldc * n, 3);
            let mut c_ref = c0.clone();
            let mut c_fast = c0.clone();
            cblas_dgemm(
                CblasColMajor,
                CblasTrans,
                CblasNoTrans,
                m as i32,
                n as i32,
                k as i32,
                -1.0,
                &a,
                lda as i32,
                &b,
                ldb as i32,
                beta,
                &mut c_ref,
                ldc as i32,
            );
            fast_dgemm_tn(m, n, k, -1.0, &a, lda, &b, ldb, beta, &mut c_fast, ldc);
            for (x, y) in c_ref.iter().zip(&c_fast) {
                assert!((x - y).abs() <= 1e-12 * (1.0 + x.abs()), "{x} vs {y}");
            }
        }
    }

    #[test]
    fn micro_kernels_agree() {
        let k = 37;
        let ap = fill(k * MR, 11);
        let bp = fill(k * NR, 12);
        let mut g: Acc = [[0.0; MR]; NR];
        micro_generic(&ap, &bp, k, &mut g);
        let mut v: Acc = [[0.0; MR]; NR];
        micro(have_avx2_fma(), &ap, &bp, k, &mut v);
        for jj in 0..NR {
            for ii in 0..MR {
                assert!((g[jj][ii] - v[jj][ii]).abs() <= 1e-13 * (1.0 + g[jj][ii].abs()));
            }
        }
    }

    #[test]
    fn fast_dtrsm_lutu_matches_reference() {
        for &(m, n) in &[
            (1usize, 1usize),
            (5, 3),
            (17, 19),
            (64, 40),
            (100, 13),
            (257, 9),
        ] {
            let (lda, ldb) = (m + 2, m + 1);
            let a = fill(lda * m, 7);
            let b0 = fill(ldb * n, 8);
            let mut b_ref = b0.clone();
            let mut b_fast = b0.clone();
            cblas_dtrsm(
                CblasColMajor,
                CblasLeft,
                CblasUpper,
                CblasTrans,
                CblasUnit,
                m as i32,
                n as i32,
                1.0,
                &a,
                lda as i32,
                &mut b_ref,
                ldb as i32,
            );
            fast_dtrsm_lutu(m, n, &a, lda, &mut b_fast, ldb);
            for (x, y) in b_ref.iter().zip(&b_fast) {
                assert!((x - y).abs() <= 1e-9 * (1.0 + x.abs()), "{x} vs {y}");
            }
        }
    }

    #[test]
    #[ignore]
    fn fast_dgemm_tn_speed() {
        let (m, n, k) = (128usize, 4000usize, 128usize);
        let a = fill(k * m, 1);
        let b = fill(k * n, 2);
        let mut c = vec![0.0; m * n];
        let flops = 2.0 * (m * n * k) as f64 * 20.0;
        let t = std::time::Instant::now();
        for _ in 0..20 {
            cblas_dgemm(
                CblasColMajor,
                CblasTrans,
                CblasNoTrans,
                m as i32,
                n as i32,
                k as i32,
                -1.0,
                &a,
                k as i32,
                &b,
                k as i32,
                0.0,
                &mut c,
                m as i32,
            );
        }
        println!(
            "ref  {:.2} GFlop/s",
            flops / t.elapsed().as_secs_f64() / 1e9
        );
        let t = std::time::Instant::now();
        for _ in 0..20 {
            fast_dgemm_tn(m, n, k, -1.0, &a, k, &b, k, 0.0, &mut c, m);
        }
        println!(
            "fast {:.2} GFlop/s",
            flops / t.elapsed().as_secs_f64() / 1e9
        );
    }
}
