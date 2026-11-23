//! HiGHS source: ipm/basiclu/basiclu_solve_dense.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::{Lu, lu_load, lu_save};
use super::lu_solve_dense::lu_solve_dense;
use super::*;

/// HiGHS: basiclu/basiclu_solve_dense.c::basiclu_solve_dense (lines 10-46)
///
/// PORT NOTE: the C `NULL` test of `Li..Lx..rhs, lhs` (`argument_missing`)
/// is dropped: slices are never `NULL` (see basiclu_factorize.rs).
pub fn basiclu_solve_dense(
    istore: &mut [LuInt],
    xstore: &mut [f64],
    li: &mut [LuInt],
    lx: &mut [f64],
    ui: &mut [LuInt],
    ux: &mut [f64],
    wi: &mut [LuInt],
    wx: &mut [f64],
    rhs: &[f64],
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
        lu_solve_dense(&mut this, rhs, lhs, trans);
    }
    lu_save(&mut this, status)
}
