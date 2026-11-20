//! HiGHS source: ipm/basiclu/basiclu_solve_sparse.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::{Lu, lu_load, lu_save};
use super::lu_solve_sparse::lu_solve_sparse;
use super::*;

/// HiGHS: basiclu/basiclu_solve_sparse.c::basiclu_solve_sparse (lines 10-63)
///
/// PORT NOTE: the C `NULL` test of the arrays and `p_nzlhs` (`argument_missing`)
/// is dropped: slices and references are never `NULL` (see basiclu_factorize.rs).
pub fn basiclu_solve_sparse(
    istore: &mut [LuInt],
    xstore: &mut [f64],
    li: &mut [LuInt],
    lx: &mut [f64],
    ui: &mut [LuInt],
    ux: &mut [f64],
    wi: &mut [LuInt],
    wx: &mut [f64],
    nzrhs: LuInt,
    irhs: &[LuInt],
    xrhs: &[f64],
    p_nzlhs: &mut LuInt,
    ilhs: &mut [LuInt],
    lhs: &mut [f64],
    trans: u8,
) -> LuInt {
    let mut this = match lu_load(istore, xstore, li, lx, ui, ux, wi, wx) {
        Ok(this) => this,
        Err(status) => return status,
    };

    let mut status = BASICLU_OK;
    if this.nupdate < 0 {
        status = BASICLU_ERROR_invalid_call;
    } else {
        /* check RHS indices */
        let mut ok = nzrhs >= 0 && nzrhs <= this.m;
        let mut n: LuInt = 0;
        while n < nzrhs && ok {
            ok = ok && irhs[n as usize] >= 0 && irhs[n as usize] < this.m;
            n += 1;
        }
        if !ok {
            status = BASICLU_ERROR_invalid_argument;
        }
    }

    if status == BASICLU_OK {
        lu_solve_sparse(&mut this, nzrhs, irhs, xrhs, p_nzlhs, ilhs, lhs, trans);
    }

    lu_save(&mut this, status)
}
