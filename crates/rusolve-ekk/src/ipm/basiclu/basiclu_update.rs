//! HiGHS source: ipm/basiclu/basiclu_update.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::{Lu, lu_load, lu_save};
use super::lu_update::lu_update;
use super::*;

/// HiGHS: basiclu/basiclu_update.c::basiclu_update (lines 10-44)
///
/// PORT NOTE: the C `NULL` test of the factor arrays is dropped: slices are
/// never `NULL` (see basiclu_factorize.rs).
pub fn basiclu_update(
    istore: &mut [LuInt],
    xstore: &mut [f64],
    li: &mut [LuInt],
    lx: &mut [f64],
    ui: &mut [LuInt],
    ux: &mut [f64],
    wi: &mut [LuInt],
    wx: &mut [f64],
    xtbl: f64,
) -> LuInt {
    let mut this = match lu_load(istore, xstore, li, lx, ui, ux, wi, wx) {
        Ok(this) => this,
        Err(status) => return status,
    };

    let status = if this.nupdate < 0 || this.ftran_for_update < 0 || this.btran_for_update < 0 {
        BASICLU_ERROR_invalid_call
    } else {
        lu_update(&mut this, xtbl)
    };
    lu_save(&mut this, status)
}
