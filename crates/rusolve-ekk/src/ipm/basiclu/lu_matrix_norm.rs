//! HiGHS source: ipm/basiclu/lu_matrix_norm.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::*;

/// HiGHS: basiclu/lu_matrix_norm.c::lu_matrix_norm (lines 13-51)
///
/// `pivotcol` / `pivotrow` are the second-list aliases `colcount_flink` /
/// `colcount_blink`; `rowsum` is `work1`.
pub fn lu_matrix_norm(
    this: &mut Lu<'_>,
    bbegin: &[LuInt],
    bend: &[LuInt],
    bi: &[LuInt],
    bx: &[f64],
) {
    debug_assert_eq!(this.nupdate, 0);

    let m = this.m;
    let rank = this.rank;
    let pivotcol: &[LuInt] = &*this.colcount_flink;
    let pivotrow: &[LuInt] = &*this.colcount_blink;
    let rowsum: &mut [f64] = &mut *this.work1;

    for i in 0..m {
        rowsum[i as usize] = 0.0;
    }
    let mut onenorm: f64 = 0.0;
    let mut infnorm: f64 = 0.0;
    for k in 0..rank {
        let jpivot = pivotcol[k as usize];
        let mut colsum: f64 = 0.0;
        let mut pos = bbegin[jpivot as usize];
        while pos < bend[jpivot as usize] {
            colsum += bx[pos as usize].abs();
            rowsum[bi[pos as usize] as usize] += bx[pos as usize].abs();
            pos += 1;
        }
        onenorm = onenorm.max(colsum);
    }
    for k in rank..m {
        let ipivot = pivotrow[k as usize];
        rowsum[ipivot as usize] += 1.0;
        onenorm = onenorm.max(1.0);
    }
    for i in 0..m {
        infnorm = infnorm.max(rowsum[i as usize]);
    }

    this.onenorm = onenorm;
    this.infnorm = infnorm;
}
