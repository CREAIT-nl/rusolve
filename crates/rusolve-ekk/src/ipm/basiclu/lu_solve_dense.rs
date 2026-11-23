//! HiGHS source: ipm/basiclu/lu_solve_dense.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_garbage_perm::lu_garbage_perm;
use super::lu_internal::Lu;
use super::*;

/// HiGHS: basiclu/lu_solve_dense.c::lu_solve_dense (lines 10-129)
///
/// C aliases (see the table on [`Lu`]): `p` = `wblink[m+1..]`, `eta_row` =
/// `rowcount_flink[m+1..]`, `Rbegin` = `rowcount_flink`, `pivotcol` =
/// `colcount_flink`, `pivotrow` = `colcount_blink`, `Ltbegin_p` =
/// `wflink[m+1..]`. `work1` holds the dense work vector (`memcpy` of `rhs`).
pub fn lu_solve_dense(this: &mut Lu<'_>, rhs: &[f64], lhs: &mut [f64], trans: u8) {
    let m = this.m;
    let nforrest = this.nforrest;
    let mu = m as usize;

    lu_garbage_perm(this);
    debug_assert!(this.pivotlen == m);

    let p = &this.wblink[mu + 1..];
    let eta_row = &this.rowcount_flink[mu + 1..];
    let pivotcol = &this.colcount_flink[..];
    let pivotrow = &this.colcount_blink[..];
    let lbegin_p = &this.lbegin_p[..];
    let ltbegin_p = &this.wflink[mu + 1..];
    let ubegin = &this.ubegin[..];
    let rbegin = &this.rowcount_flink[..];
    let wbegin = &this.wbegin[..];
    let wend = &this.wend[..];
    let col_pivot = &this.col_pivot[..];
    let row_pivot = &this.row_pivot[..];
    let lindex = &this.lindex[..];
    let lvalue = &this.lvalue[..];
    let uindex = &this.uindex[..];
    let uvalue = &this.uvalue[..];
    let windex = &this.windex[..];
    let wvalue = &this.wvalue[..];
    let work1 = &mut this.work1[..];

    if trans == b't' || trans == b'T' {
        // ----------------------- //
        // Solve transposed system //
        // ----------------------- //

        work1[..mu].copy_from_slice(&rhs[..mu]);

        // Solve with U'.
        for k in 0..mu {
            let jpivot = pivotcol[k] as usize;
            let ipivot = pivotrow[k] as usize;
            let x: f64 = work1[jpivot] / col_pivot[jpivot];
            for pos in wbegin[jpivot]..wend[jpivot] {
                let pos = pos as usize;
                work1[windex[pos] as usize] -= x * wvalue[pos];
            }
            lhs[ipivot] = x;
        }

        // Solve with update ETAs backwards.
        for t in (0..nforrest).rev() {
            let ipivot = eta_row[t as usize] as usize;
            let x: f64 = lhs[ipivot];
            for pos in rbegin[t as usize]..rbegin[t as usize + 1] {
                let pos = pos as usize;
                let i = lindex[pos] as usize;
                lhs[i] -= x * lvalue[pos];
            }
        }

        // Solve with L'.
        for k in (0..mu).rev() {
            let mut x: f64 = 0.0;
            let mut pos = lbegin_p[k];
            loop {
                let i = lindex[pos as usize];
                if i < 0 {
                    break;
                }
                x += lhs[i as usize] * lvalue[pos as usize];
                pos += 1;
            }
            lhs[p[k] as usize] -= x;
        }
    } else {
        // -------------------- //
        // Solve forward system //
        // -------------------- //

        work1[..mu].copy_from_slice(&rhs[..mu]);

        // Solve with L.
        for k in 0..mu {
            let mut x: f64 = 0.0;
            let mut pos = ltbegin_p[k];
            loop {
                let i = lindex[pos as usize];
                if i < 0 {
                    break;
                }
                x += work1[i as usize] * lvalue[pos as usize];
                pos += 1;
            }
            work1[p[k] as usize] -= x;
        }

        // Solve with update ETAs.
        let mut pos = rbegin[0];
        for t in 0..nforrest {
            let ipivot = eta_row[t as usize] as usize;
            let mut x: f64 = 0.0;
            while pos < rbegin[t as usize + 1] {
                x += work1[lindex[pos as usize] as usize] * lvalue[pos as usize];
                pos += 1;
            }
            work1[ipivot] -= x;
        }

        // Solve with U.
        for k in (0..mu).rev() {
            let jpivot = pivotcol[k] as usize;
            let ipivot = pivotrow[k] as usize;
            let x: f64 = work1[ipivot] / row_pivot[ipivot];
            let mut pos = ubegin[ipivot];
            loop {
                let i = uindex[pos as usize];
                if i < 0 {
                    break;
                }
                work1[i as usize] -= x * uvalue[pos as usize];
                pos += 1;
            }
            lhs[jpivot] = x;
        }
    }
}
