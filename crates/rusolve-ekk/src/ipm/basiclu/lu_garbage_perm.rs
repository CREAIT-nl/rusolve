//! HiGHS source: ipm/basiclu/lu_garbage_perm.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::*;

/// HiGHS: basiclu/lu_garbage_perm.c::lu_garbage_perm (lines 25-53)
///
/// `pivotcol` / `pivotrow` are the aliases `colcount_flink` / `colcount_blink`,
/// `marked` is `iwork0`.
pub fn lu_garbage_perm(this: &mut Lu<'_>) {
    let m = this.m;
    let pivotlen = this.pivotlen;

    if pivotlen > m {
        this.marker += 1;
        let marker = this.marker;
        let mut put = pivotlen;
        let mut get = pivotlen - 1;
        while get >= 0 {
            let j = this.colcount_flink[get as usize];
            if this.iwork0[j as usize] != marker {
                this.iwork0[j as usize] = marker;
                put -= 1;
                this.colcount_flink[put as usize] = j;
                this.colcount_blink[put as usize] = this.colcount_blink[get as usize];
            }
            get -= 1;
        }
        debug_assert_eq!(put + m, pivotlen);
        // memmove(pivotcol, pivotcol+put, m*sizeof(lu_int)), same for pivotrow
        this.colcount_flink
            .copy_within(put as usize..(put + m) as usize, 0);
        this.colcount_blink
            .copy_within(put as usize..(put + m) as usize, 0);
        this.pivotlen = m;
    }
}
