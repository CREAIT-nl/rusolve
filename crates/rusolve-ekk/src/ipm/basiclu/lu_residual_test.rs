//! HiGHS source: ipm/basiclu/lu_residual_test.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::lu_matrix_norm::lu_matrix_norm;
use super::*;

/// HiGHS: basiclu/lu_residual_test.c::lu_onenorm (lines 12-21)
fn lu_onenorm(m: LuInt, x: &[f64]) -> f64 {
    let mut d: f64 = 0.0;
    for i in 0..m {
        d += x[i as usize].abs();
    }
    d
}

/// HiGHS: basiclu/lu_residual_test.c::lu_residual_test (lines 23-155)
///
/// PORT NOTE: `rhs` is `work0`, `lhs` is `work1`. `p` is `wblink[m+1..]`,
/// `Ltbegin_p` is `wflink[m+1..]`, `pivotcol`/`pivotrow` are the aliases
/// `colcount_flink`/`colcount_blink`. C resets `rhs` after the final
/// `lu_matrix_norm` call, which does not touch `work0`; the reset is written
/// through `this.work0` once `lu_matrix_norm` has taken `this`.
pub fn lu_residual_test(
    this: &mut Lu<'_>,
    bbegin: &[LuInt],
    bend: &[LuInt],
    bi: &[LuInt],
    bx: &[f64],
) {
    let m = this.m;
    let rank = this.rank;
    let p: &[LuInt] = &this.wblink[(m + 1) as usize..];
    let pivotcol: &[LuInt] = &*this.colcount_flink;
    let pivotrow: &[LuInt] = &*this.colcount_blink;
    let lbegin_p: &[LuInt] = &*this.lbegin_p;
    let ltbegin_p: &[LuInt] = &this.wflink[(m + 1) as usize..];
    let ubegin: &[LuInt] = &*this.ubegin;
    let row_pivot: &[f64] = &*this.row_pivot;
    let lindex: &[LuInt] = &*this.lindex;
    let lvalue: &[f64] = &*this.lvalue;
    let uindex: &[LuInt] = &*this.uindex;
    let uvalue: &[f64] = &*this.uvalue;
    let rhs: &mut [f64] = &mut *this.work0;
    let lhs: &mut [f64] = &mut *this.work1;

    let mut d: f64;

    /* --------------------------------- */
    /* Residual Test with Forward System */
    /* --------------------------------- */

    /* Compute lhs = L\rhs and build rhs on-the-fly. */
    for k in 0..m {
        d = 0.0;
        let mut pos = ltbegin_p[k as usize];
        while lindex[pos as usize] >= 0 {
            let i = lindex[pos as usize];
            d += lhs[i as usize] * lvalue[pos as usize];
            pos += 1;
        }
        let ipivot = p[k as usize] as usize;
        rhs[ipivot] = if d <= 0.0 { 1.0 } else { -1.0 };
        lhs[ipivot] = rhs[ipivot] - d;
    }

    /* Overwrite lhs by U\lhs. */
    for k in (0..m).rev() {
        let ipivot = pivotrow[k as usize] as usize;
        lhs[ipivot] /= row_pivot[ipivot];
        d = lhs[ipivot];
        let mut pos = ubegin[ipivot];
        while uindex[pos as usize] >= 0 {
            let i = uindex[pos as usize];
            lhs[i as usize] -= d * uvalue[pos as usize];
            pos += 1;
        }
    }

    /* Overwrite rhs by the residual rhs-B*lhs. */
    for k in 0..rank {
        let ipivot = pivotrow[k as usize] as usize;
        let jpivot = pivotcol[k as usize] as usize;
        d = lhs[ipivot];
        let mut pos = bbegin[jpivot];
        while pos < bend[jpivot] {
            rhs[bi[pos as usize] as usize] -= d * bx[pos as usize];
            pos += 1;
        }
    }
    for k in rank..m {
        let ipivot = pivotrow[k as usize] as usize;
        rhs[ipivot] -= lhs[ipivot];
    }
    let norm_ftran = lu_onenorm(m, lhs);
    let norm_ftran_res = lu_onenorm(m, rhs);

    /* ---------------------------------- */
    /* Residual Test with Backward System */
    /* ---------------------------------- */

    /* Compute lhs = U'\rhs and build rhs on-the-fly. */
    for k in 0..m {
        let ipivot = pivotrow[k as usize] as usize;
        d = 0.0;
        let mut pos = ubegin[ipivot];
        while uindex[pos as usize] >= 0 {
            let i = uindex[pos as usize];
            d += lhs[i as usize] * uvalue[pos as usize];
            pos += 1;
        }
        rhs[ipivot] = if d <= 0.0 { 1.0 } else { -1.0 };
        lhs[ipivot] = (rhs[ipivot] - d) / row_pivot[ipivot];
    }

    /* Overwrite lhs by L'\lhs. */
    for k in (0..m).rev() {
        d = 0.0;
        let mut pos = lbegin_p[k as usize];
        while lindex[pos as usize] >= 0 {
            let i = lindex[pos as usize];
            d += lhs[i as usize] * lvalue[pos as usize];
            pos += 1;
        }
        lhs[p[k as usize] as usize] -= d;
    }

    /* Overwrite rhs by the residual rhs-B'*lhs. */
    for k in 0..rank {
        let ipivot = pivotrow[k as usize] as usize;
        let jpivot = pivotcol[k as usize] as usize;
        d = 0.0;
        let mut pos = bbegin[jpivot];
        while pos < bend[jpivot] {
            d += lhs[bi[pos as usize] as usize] * bx[pos as usize];
            pos += 1;
        }
        rhs[ipivot] -= d;
    }
    for k in rank..m {
        let ipivot = pivotrow[k as usize] as usize;
        rhs[ipivot] -= lhs[ipivot];
    }
    let norm_btran = lu_onenorm(m, lhs);
    let norm_btran_res = lu_onenorm(m, rhs);

    /* -------- */
    /* Finalize */
    /* -------- */

    lu_matrix_norm(this, bbegin, bend, bi, bx);
    debug_assert!(this.onenorm > 0.0);
    debug_assert!(this.infnorm > 0.0);
    this.residual_test = (norm_ftran_res / (m as f64 + this.onenorm * norm_ftran))
        .max(norm_btran_res / (m as f64 + this.infnorm * norm_btran));

    /* reset workspace */
    for i in 0..m {
        this.work0[i as usize] = 0.0;
    }
}
