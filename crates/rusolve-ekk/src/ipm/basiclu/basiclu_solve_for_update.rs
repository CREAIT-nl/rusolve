//! HiGHS source: ipm/basiclu/basiclu_solve_for_update.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::{Lu, lu_load, lu_save};
use super::lu_solve_for_update::lu_solve_for_update;
use super::*;

/// HiGHS: basiclu/basiclu_solve_for_update.c::basiclu_solve_for_update (lines 10-79)
///
/// `xrhs` (unused for `trans == 'T'`), `p_nzlhs`, `ilhs`, `lhs` are `NULL`
/// in the "without lhs" calls of `BasicLu`.
///
/// PORT NOTE: the C `NULL` test of the arrays (`argument_missing`) is dropped:
/// slices are never `NULL`. The `xrhs` test is kept, because `xrhs` is an
/// `Option`.
///
/// BASICLU quirk: for `trans == 't' / 'T'` the C checks `irhs[0]` without
/// looking at `nzrhs`, so it reads `irhs[0]` even when `nzrhs == 0`. Kept.
pub fn basiclu_solve_for_update(
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
    xrhs: Option<&[f64]>,
    p_nzlhs: Option<&mut LuInt>,
    ilhs: Option<&mut [LuInt]>,
    lhs: Option<&mut [f64]>,
    trans: u8,
) -> LuInt {
    let mut this = match lu_load(istore, xstore, li, lx, ui, ux, wi, wx) {
        Ok(this) => this,
        Err(status) => return status,
    };

    let mut status = BASICLU_OK;
    if trans != b't' && trans != b'T' && xrhs.is_none() {
        status = BASICLU_ERROR_argument_missing;
    } else if this.nupdate < 0 {
        status = BASICLU_ERROR_invalid_call;
    } else if this.nforrest == this.m {
        status = BASICLU_ERROR_maximum_updates;
    } else {
        /* check RHS indices */
        let mut ok: bool;
        if trans == b't' || trans == b'T' {
            ok = irhs[0] >= 0 && irhs[0] < this.m;
        } else {
            ok = nzrhs >= 0 && nzrhs <= this.m;
            let mut n: LuInt = 0;
            while n < nzrhs && ok {
                ok = ok && irhs[n as usize] >= 0 && irhs[n as usize] < this.m;
                n += 1;
            }
        }
        if !ok {
            status = BASICLU_ERROR_invalid_argument;
        }
    }

    if status == BASICLU_OK {
        /* may request reallocation */
        status = lu_solve_for_update(&mut this, nzrhs, irhs, xrhs, p_nzlhs, ilhs, lhs, trans);
    }

    lu_save(&mut this, status)
}
