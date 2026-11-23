//! HiGHS source: ipm/basiclu/basiclu_get_factors.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::{Lu, lu_load, lu_save};
use super::*;

/// HiGHS: basiclu/basiclu_get_factors.c::basiclu_get_factors (lines 10-148)
///
/// `rowperm`, `colperm` and each of the L/U triples may be `NULL` in C.
pub fn basiclu_get_factors(
    istore: &mut [LuInt],
    xstore: &mut [f64],
    li: &mut [LuInt],
    lx: &mut [f64],
    ui: &mut [LuInt],
    ux: &mut [f64],
    wi: &mut [LuInt],
    wx: &mut [f64],
    rowperm: Option<&mut [LuInt]>,
    colperm: Option<&mut [LuInt]>,
    lcolptr: Option<&mut [LuInt]>,
    lrowidx: Option<&mut [LuInt]>,
    lvalue_: Option<&mut [f64]>,
    ucolptr: Option<&mut [LuInt]>,
    urowidx: Option<&mut [LuInt]>,
    uvalue_: Option<&mut [f64]>,
) -> LuInt {
    let mut this: Lu<'_> = match lu_load(istore, xstore, li, lx, ui, ux, wi, wx) {
        Ok(this) => this,
        Err(status) => return status,
    };
    if this.nupdate != 0 {
        let status = BASICLU_ERROR_invalid_call;
        return lu_save(&mut this, status);
    }
    let m = this.m;
    let mu = m as usize;

    // C pivotrow / pivotcol are the second-list aliases colcount_blink /
    // colcount_flink (see the table on `Lu`).
    if let Some(rowperm) = rowperm {
        rowperm[..mu].copy_from_slice(&this.colcount_blink[..mu]);
    }
    if let Some(colperm) = colperm {
        colperm[..mu].copy_from_slice(&this.colcount_flink[..mu]);
    }

    if let (Some(lcolptr), Some(lrowidx), Some(lvalue_)) = (lcolptr, lrowidx, lvalue_) {
        let lbegin_p = &this.lbegin_p[..];
        let ltbegin_p = &this.wflink[mu + 1..]; // Ltbegin_p = Wflink + m + 1
        let lindex = &this.lindex[..];
        let lvalue = &this.lvalue[..];
        let p = &this.wblink[mu + 1..]; // p = Wblink + m + 1
        let colptr = &mut this.rowcount_blink[..]; // iwork1: size m workspace

        // L[:,k] will hold the elimination factors from the k-th pivot step.
        // First set the column pointers and store the unit diagonal elements
        // at the front of each column. Then scatter each row of L' into the
        // columnwise L so that the row indices become sorted.
        let mut put: LuInt = 0;
        for k in 0..mu {
            lcolptr[k] = put;
            lrowidx[put as usize] = k as LuInt;
            lvalue_[put as usize] = 1.0;
            put += 1;
            colptr[p[k] as usize] = put; // next free position in column
            put += lbegin_p[k + 1] - lbegin_p[k] - 1;
            // subtract 1 because internal storage uses (-1) terminators
        }
        lcolptr[mu] = put;
        debug_assert!(put == this.lnz + m);

        for k in 0..mu {
            let mut pos = ltbegin_p[k];
            loop {
                let i = lindex[pos as usize];
                if i < 0 {
                    break;
                }
                put = colptr[i as usize];
                colptr[i as usize] += 1;
                lrowidx[put as usize] = k as LuInt;
                lvalue_[put as usize] = lvalue[pos as usize];
                pos += 1;
            }
        }

        for k in 0..mu {
            debug_assert!(colptr[p[k] as usize] == lcolptr[k + 1]);
        }
    }

    if let (Some(ucolptr), Some(urowidx), Some(uvalue_)) = (ucolptr, urowidx, uvalue_) {
        let wbegin = &this.wbegin[..]; // Wbegin
        let wend = &this.wend[..]; // Wend
        let windex = &this.windex[..];
        let wvalue = &this.wvalue[..];
        let col_pivot = &this.col_pivot[..];
        let pivotcol = &this.colcount_flink[..];
        let colptr = &mut this.rowcount_blink[..]; // iwork1: size m workspace

        // U[:,k] will hold the column of B from the k-th pivot step.
        // First set the column pointers and store the pivot element at the end
        // of each column. Then scatter each row of U' into the columnwise U so
        // that the row indices become sorted.
        colptr[..mu].fill(0); // column counts
        for j in 0..mu {
            for pos in wbegin[j]..wend[j] {
                colptr[windex[pos as usize] as usize] += 1;
            }
        }
        let mut put: LuInt = 0;
        for k in 0..mu {
            // set column pointers
            let j = pivotcol[k] as usize;
            ucolptr[k] = put;
            put += colptr[j];
            colptr[j] = ucolptr[k]; // next free position in column
            urowidx[put as usize] = k as LuInt;
            uvalue_[put as usize] = col_pivot[j];
            put += 1;
        }
        ucolptr[mu] = put;
        debug_assert!(put == this.unz + m);
        for k in 0..mu {
            // scatter row k
            let j = pivotcol[k] as usize;
            for pos in wbegin[j]..wend[j] {
                let idx = windex[pos as usize] as usize;
                put = colptr[idx];
                colptr[idx] += 1;
                urowidx[put as usize] = k as LuInt;
                uvalue_[put as usize] = wvalue[pos as usize];
            }
        }

        for k in 0..mu {
            debug_assert!(colptr[pivotcol[k] as usize] == ucolptr[k + 1] - 1);
        }
        for k in 0..mu {
            debug_assert!(urowidx[(ucolptr[k + 1] - 1) as usize] == k as LuInt);
        }
    }

    // BASICLU quirk: the C success path returns BASICLU_OK without lu_save
    // (only the error path above saves). This routine changes no scalar of
    // `this`, so skipping the write-back leaves istore/xstore as lu_save
    // would have written them.
    BASICLU_OK
}
