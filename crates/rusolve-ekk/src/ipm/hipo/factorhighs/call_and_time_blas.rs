//! HiGHS source: ipm/hipo/factorhighs/CallAndTimeBlas.h, ipm/hipo/factorhighs/CallAndTimeBlas.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The only callers of [`crate::ipm::hipo::blas`]: char arguments (`'N'`,
//! `'T'`, `'U'`, `'L'`, ...) map to the CBLAS enums with the `TRANS`, `UPLO`,
//! `DIAG`, `SIDE` macros (`x == 'N' ? CblasNoTrans : CblasTrans`, ...), order
//! is always `CblasColMajor`; the `HIPO_CLOCK_*` timing is empty.

use super::super::Int;
use super::super::blas::*;
use super::data_collector::DataCollector;

// rusolve: phase profile behind env `RUSOLVE_HIPO_PROF` (summed thread time per BLAS routine
// and per dense-factor phase).
use std::sync::atomic::{AtomicU64, Ordering as AtomOrd};
pub const PROF_NAMES: [&str; 29] = [
    "daxpy",
    "dcopy",
    "dscal",
    "dgemv",
    "dtpsv",
    "dtrsv",
    "dger",
    "dgemm",
    "dsyrk",
    "dtrsm",
    "ph:diag",
    "ph:solvecols",
    "ph:frontal",
    "ph:schur",
    "ph:schurcopy",
    "fa:new",
    "fa:fhnew",
    "fa:assemble",
    "fa:dense",
    "fa:extreme",
    "fa:terminate",
    "it:prepare",
    "it:predictor",
    "it:correctors",
    "it:makestep",
    "it:solve6x6",
    "it:residual",
    "it:omega",
    "it:refinecopy",
];
pub static PROF_NS: [AtomicU64; 29] = [const { AtomicU64::new(0) }; 29];
pub static PROF_CALLS: [AtomicU64; 29] = [const { AtomicU64::new(0) }; 29];
fn prof_on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("RUSOLVE_HIPO_PROF").is_some())
}
/// Front-size buckets (front dimension k + clique: <8, <16, <32, <64, <128, <256, <512, more):
/// supernode count, assemble ns, dense-factor ns and dense flops (k^3/3 + k^2 c + k c^2).
pub static BK_CNT: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];
pub static BK_ASM_NS: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];
pub static BK_DENSE_NS: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];
pub static BK_FLOPS: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];
pub fn bucket_on() -> bool {
    prof_on()
}
/// Record one supernode (k pivots, clique size c) with its assemble and dense-factor time.
pub fn bucket_add(k: i64, c: i64, asm_ns: u64, dense_ns: u64) {
    let n = k + c;
    let b = match n {
        0..=7 => 0,
        8..=15 => 1,
        16..=31 => 2,
        32..=63 => 3,
        64..=127 => 4,
        128..=255 => 5,
        256..=511 => 6,
        _ => 7,
    };
    let (kf, cf) = (k as f64, c as f64);
    let fl = kf * kf * kf / 3.0 + kf * kf * cf + kf * cf * cf;
    BK_CNT[b].fetch_add(1, AtomOrd::Relaxed);
    BK_ASM_NS[b].fetch_add(asm_ns, AtomOrd::Relaxed);
    BK_DENSE_NS[b].fetch_add(dense_ns, AtomOrd::Relaxed);
    BK_FLOPS[b].fetch_add(fl as u64, AtomOrd::Relaxed);
}
pub struct ProfGuard(usize, Option<std::time::Instant>);
impl Drop for ProfGuard {
    fn drop(&mut self) {
        if let Some(t) = self.1 {
            PROF_NS[self.0].fetch_add(t.elapsed().as_nanos() as u64, AtomOrd::Relaxed);
            PROF_CALLS[self.0].fetch_add(1, AtomOrd::Relaxed);
        }
    }
}
/// `RUSOLVE_HIPO_PROF=2`: per-BLAS-call guards (indices 0-9, hundreds of millions of calls) are off,
/// so the phase and iteration timers are not inflated by the clock reads.
fn prof_blas_on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("RUSOLVE_HIPO_PROF").is_ok_and(|v| v != "2"))
}
pub fn prof(i: usize) -> ProfGuard {
    let on = if i < 10 { prof_blas_on() } else { prof_on() };
    ProfGuard(i, on.then(std::time::Instant::now))
}
/// One line per routine with time and calls; empty when profiling is off.
pub fn prof_report() -> String {
    if !prof_on() {
        return String::new();
    }
    let mut s = String::from("HiPO profile (summed thread seconds, calls):\n");
    for i in 0..PROF_NAMES.len() {
        s.push_str(&format!(
            "  {:13} {:9.3} {:10}\n",
            PROF_NAMES[i],
            PROF_NS[i].load(AtomOrd::Relaxed) as f64 * 1e-9,
            PROF_CALLS[i].load(AtomOrd::Relaxed)
        ));
    }
    s.push_str("front size   supernodes  assemble_s  dense_s  dense_GF/s\n");
    let lab = ["<8", "<16", "<32", "<64", "<128", "<256", "<512", ">=512"];
    for b in 0..8 {
        let d = BK_DENSE_NS[b].load(AtomOrd::Relaxed) as f64 * 1e-9;
        let f = BK_FLOPS[b].load(AtomOrd::Relaxed) as f64;
        s.push_str(&format!(
            "  {:8} {:10} {:10.3} {:8.3} {:8.2}\n",
            lab[b],
            BK_CNT[b].load(AtomOrd::Relaxed),
            BK_ASM_NS[b].load(AtomOrd::Relaxed) as f64 * 1e-9,
            d,
            if d > 0.0 { f / d * 1e-9 } else { 0.0 }
        ));
    }
    s
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_daxpy (lines 20-25)
pub fn callAndTime_daxpy(
    n: Int,
    da: f64,
    dx: &[f64],
    incx: Int,
    dy: &mut [f64],
    incy: Int,
    data: &DataCollector,
) {
    let _prof = prof(0);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_axpy) dropped (lines 22, 24).
    let _ = data;
    cblas_daxpy(n, da, dx, incx, dy, incy);
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dcopy (lines 27-32)
pub fn callAndTime_dcopy(
    n: Int,
    dx: &[f64],
    incx: Int,
    dy: &mut [f64],
    incy: Int,
    data: &DataCollector,
) {
    let _prof = prof(1);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_copy) dropped (lines 29, 31).
    let _ = data;
    cblas_dcopy(n, dx, incx, dy, incy);
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dscal (lines 34-39)
pub fn callAndTime_dscal(n: Int, da: f64, dx: &mut [f64], incx: Int, data: &DataCollector) {
    let _prof = prof(2);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_scal) dropped (lines 36, 38).
    let _ = data;
    cblas_dscal(n, da, dx, incx);
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dswap (lines 41-46)
pub fn callAndTime_dswap(
    n: Int,
    dx: &mut [f64],
    incx: Int,
    dy: &mut [f64],
    incy: Int,
    data: &DataCollector,
) {
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_swap) dropped (lines 43, 45).
    let _ = data;
    cblas_dswap(n, dx, incx, dy, incy);
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dgemv (lines 50-57)
pub fn callAndTime_dgemv(
    trans: u8,
    m: Int,
    n: Int,
    alpha: f64,
    A: &[f64],
    lda: Int,
    x: &[f64],
    incx: Int,
    beta: f64,
    y: &mut [f64],
    incy: Int,
    data: &DataCollector,
) {
    let _prof = prof(3);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_gemv) dropped (lines 53, 56).
    let _ = data;
    // TRANS(trans): (x) == 'N' ? CblasNoTrans : CblasTrans
    let trans_cblas = if trans == b'N' {
        CblasNoTrans
    } else {
        CblasTrans
    };
    cblas_dgemv(
        CblasColMajor,
        trans_cblas,
        m,
        n,
        alpha,
        A,
        lda,
        x,
        incx,
        beta,
        y,
        incy,
    );
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dtpsv (lines 59-66)
pub fn callAndTime_dtpsv(
    uplo: u8,
    trans: u8,
    diag: u8,
    n: Int,
    ap: &[f64],
    x: &mut [f64],
    incx: Int,
    data: &DataCollector,
) {
    let _prof = prof(4);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_tpsv) dropped (lines 62, 65).
    let _ = data;
    // UPLO(uplo): (x) == 'U' ? CblasUpper : CblasLower
    let uplo_cblas = if uplo == b'U' { CblasUpper } else { CblasLower };
    // TRANS(trans): (x) == 'N' ? CblasNoTrans : CblasTrans
    let trans_cblas = if trans == b'N' {
        CblasNoTrans
    } else {
        CblasTrans
    };
    // DIAG(diag): (x) == 'N' ? CblasNonUnit : CblasUnit
    let diag_cblas = if diag == b'N' {
        CblasNonUnit
    } else {
        CblasUnit
    };
    cblas_dtpsv(
        CblasColMajor,
        uplo_cblas,
        trans_cblas,
        diag_cblas,
        n,
        ap,
        x,
        incx,
    );
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dtrsv (lines 68-74)
pub fn callAndTime_dtrsv(
    uplo: u8,
    trans: u8,
    diag: u8,
    n: Int,
    A: &[f64],
    lda: Int,
    x: &mut [f64],
    incx: Int,
    data: &DataCollector,
) {
    let _prof = prof(5);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_trsv) dropped (lines 70, 73).
    let _ = data;
    // UPLO(uplo): (x) == 'U' ? CblasUpper : CblasLower
    let uplo_cblas = if uplo == b'U' { CblasUpper } else { CblasLower };
    // TRANS(trans): (x) == 'N' ? CblasNoTrans : CblasTrans
    let trans_cblas = if trans == b'N' {
        CblasNoTrans
    } else {
        CblasTrans
    };
    // DIAG(diag): (x) == 'N' ? CblasNonUnit : CblasUnit
    let diag_cblas = if diag == b'N' {
        CblasNonUnit
    } else {
        CblasUnit
    };
    cblas_dtrsv(
        CblasColMajor,
        uplo_cblas,
        trans_cblas,
        diag_cblas,
        n,
        A,
        lda,
        x,
        incx,
    );
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dger (lines 76-82)
pub fn callAndTime_dger(
    m: Int,
    n: Int,
    alpha: f64,
    x: &[f64],
    incx: Int,
    y: &[f64],
    incy: Int,
    A: &mut [f64],
    lda: Int,
    data: &DataCollector,
) {
    let _prof = prof(6);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_ger) dropped (lines 79, 81).
    let _ = data;
    cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dgemm (lines 86-94)
pub fn callAndTime_dgemm(
    transa: u8,
    transb: u8,
    m: Int,
    n: Int,
    k: Int,
    alpha: f64,
    A: &[f64],
    lda: Int,
    B: &[f64],
    ldb: Int,
    beta: f64,
    C: &mut [f64],
    ldc: Int,
    data: &DataCollector,
) {
    let _prof = prof(7);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_gemm) dropped (lines 90, 93).
    let _ = data;
    // rusolve lever `ipm_hipo_fast_blas`: packed micro-kernel for the T,N shape.
    if transa == b'T'
        && transb == b'N'
        && m > 0
        && n > 0
        && k > 0
        && crate::ipm::hipo::blas::FAST_BLAS.load(std::sync::atomic::Ordering::Relaxed)
    {
        crate::ipm::hipo::blas::fast_dgemm_tn(
            m as usize,
            n as usize,
            k as usize,
            alpha,
            A,
            lda as usize,
            B,
            ldb as usize,
            beta,
            C,
            ldc as usize,
        );
        return;
    }
    // TRANS(transa): (x) == 'N' ? CblasNoTrans : CblasTrans
    let transa_cblas = if transa == b'N' {
        CblasNoTrans
    } else {
        CblasTrans
    };
    // TRANS(transb): (x) == 'N' ? CblasNoTrans : CblasTrans
    let transb_cblas = if transb == b'N' {
        CblasNoTrans
    } else {
        CblasTrans
    };
    cblas_dgemm(
        CblasColMajor,
        transa_cblas,
        transb_cblas,
        m,
        n,
        k,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc,
    );
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dsyrk (lines 96-103)
pub fn callAndTime_dsyrk(
    uplo: u8,
    trans: u8,
    n: Int,
    k: Int,
    alpha: f64,
    A: &[f64],
    lda: Int,
    beta: f64,
    C: &mut [f64],
    ldc: Int,
    data: &DataCollector,
) {
    let _prof = prof(8);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_syrk) dropped (lines 99, 102).
    let _ = data;
    // UPLO(uplo): (x) == 'U' ? CblasUpper : CblasLower
    let uplo_cblas = if uplo == b'U' { CblasUpper } else { CblasLower };
    // TRANS(trans): (x) == 'N' ? CblasNoTrans : CblasTrans
    let trans_cblas = if trans == b'N' {
        CblasNoTrans
    } else {
        CblasTrans
    };
    cblas_dsyrk(
        CblasColMajor,
        uplo_cblas,
        trans_cblas,
        n,
        k,
        alpha,
        A,
        lda,
        beta,
        C,
        ldc,
    );
}

/// HiGHS: CallAndTimeBlas.cpp::callAndTime_dtrsm (lines 105-112)
pub fn callAndTime_dtrsm(
    side: u8,
    uplo: u8,
    trans: u8,
    diag: u8,
    m: Int,
    n: Int,
    alpha: f64,
    A: &[f64],
    lda: Int,
    B: &mut [f64],
    ldb: Int,
    data: &DataCollector,
) {
    let _prof = prof(9);
    // HIPO_CLOCK_CREATE / HIPO_CLOCK_STOP(3, data, kTimeBlas_trsm) dropped (lines 108, 111).
    let _ = data;
    // rusolve lever `ipm_hipo_fast_blas`: grouped-column solve for the L,U,T,U shape.
    if side == b'L'
        && uplo == b'U'
        && trans == b'T'
        && diag == b'U'
        && alpha == 1.0
        && m > 0
        && n > 0
        && crate::ipm::hipo::blas::FAST_BLAS.load(std::sync::atomic::Ordering::Relaxed)
    {
        crate::ipm::hipo::blas::fast_dtrsm_lutu(
            m as usize,
            n as usize,
            A,
            lda as usize,
            B,
            ldb as usize,
        );
        return;
    }
    // SIDE(side): (x) == 'L' ? CblasLeft : CblasRight
    let side_cblas = if side == b'L' { CblasLeft } else { CblasRight };
    // UPLO(uplo): (x) == 'U' ? CblasUpper : CblasLower
    let uplo_cblas = if uplo == b'U' { CblasUpper } else { CblasLower };
    // TRANS(trans): (x) == 'N' ? CblasNoTrans : CblasTrans
    let trans_cblas = if trans == b'N' {
        CblasNoTrans
    } else {
        CblasTrans
    };
    // DIAG(diag): (x) == 'N' ? CblasNonUnit : CblasUnit
    let diag_cblas = if diag == b'N' {
        CblasNonUnit
    } else {
        CblasUnit
    };
    cblas_dtrsm(
        CblasColMajor,
        side_cblas,
        uplo_cblas,
        trans_cblas,
        diag_cblas,
        m,
        n,
        alpha,
        A,
        lda,
        B,
        ldb,
    );
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn dcopy_passthrough() {
        let d = DataCollector::new();
        let mut y = [0.0; 2];
        callAndTime_dcopy(2, &[4.0, 5.0], 1, &mut y, 1, &d);
        assert_eq!(y, [4.0, 5.0]);
    }
}
