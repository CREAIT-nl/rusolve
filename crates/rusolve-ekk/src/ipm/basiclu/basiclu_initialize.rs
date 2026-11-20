//! HiGHS source: ipm/basiclu/basiclu_initialize.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_initialize::lu_initialize;
use super::lu_internal::Lu;
use super::*;

/// HiGHS: basiclu/basiclu_initialize.c::basiclu_initialize (lines 10-24)
///
/// PORT NOTE: the C `!istore || !xstore` test is dropped: slices are never
/// `NULL`.
pub fn basiclu_initialize(m: LuInt, istore: &mut [LuInt], xstore: &mut [f64]) -> LuInt {
    if m <= 0 {
        return BASICLU_ERROR_invalid_argument;
    }

    lu_initialize(m, istore, xstore);
    BASICLU_OK
}
