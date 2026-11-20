//! HiGHS source: ipm/basiclu/lu_initialize.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU source)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::{lu_load, lu_reset, lu_save};
use super::*;

/// HiGHS: basiclu/lu_initialize.c::lu_initialize (lines 18-56)
///
/// PORT NOTE: C calls `lu_load(&this, istore, xstore, NULL, ...)`. The factor
/// arrays are passed as empty slices (the C NULL): they are only stored, never
/// read, by `lu_reset` and `lu_save`.
pub fn lu_initialize(m: LuInt, istore: &mut [LuInt], xstore: &mut [f64]) {
    // set constant entries
    istore[0] = BASICLU_HASH;
    xstore[0] = BASICLU_HASH as f64;
    xstore[BASICLU_DIM] = m as f64;

    // set default parameters
    xstore[BASICLU_MEMORYL] = 0.0;
    xstore[BASICLU_MEMORYU] = 0.0;
    xstore[BASICLU_MEMORYW] = 0.0;
    xstore[BASICLU_DROP_TOLERANCE] = 1e-20;
    xstore[BASICLU_ABS_PIVOT_TOLERANCE] = 1e-14;
    xstore[BASICLU_REL_PIVOT_TOLERANCE] = 0.1;
    xstore[BASICLU_BIAS_NONZEROS] = 1.0;
    xstore[BASICLU_MAXN_SEARCH_PIVOT] = 3.0;
    xstore[BASICLU_PAD] = 4.0;
    xstore[BASICLU_STRETCH] = 0.3;
    xstore[BASICLU_COMPRESSION_THRESHOLD] = 0.5;
    xstore[BASICLU_SPARSE_THRESHOLD] = 0.05;
    xstore[BASICLU_REMOVE_COLUMNS] = 0.0;
    xstore[BASICLU_SEARCH_ROWS] = 1.0;
    xstore[BASICLU_TIME_START] = 0.0;
    xstore[BASICLU_TIME_LIMIT] = f64::INFINITY;

    // initialize global counters
    xstore[BASICLU_NFACTORIZE] = 0.0;
    xstore[BASICLU_NUPDATE_TOTAL] = 0.0;
    xstore[BASICLU_NFORREST_TOTAL] = 0.0;
    xstore[BASICLU_NSYMPERM_TOTAL] = 0.0;
    xstore[BASICLU_TIME_FACTORIZE_TOTAL] = 0.0;
    xstore[BASICLU_TIME_SOLVE_TOTAL] = 0.0;
    xstore[BASICLU_TIME_UPDATE_TOTAL] = 0.0;

    // lu_reset() and lu_save() initializes the remaining slots
    // The hash was set above, so lu_load cannot reject the store.
    let mut this = lu_load(
        istore,
        xstore,
        &mut [],
        &mut [],
        &mut [],
        &mut [],
        &mut [],
        &mut [],
    )
    .expect("lu_initialize: store shorter than the scalar block");
    lu_reset(&mut this);
    lu_save(&mut this, BASICLU_OK);
}
