//! HiGHS source: ipm/basiclu/lu_factorize_bump.c
//!
//! Copyright (C) 2016-2019  ERGO-Code (BASICLU source)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::{Lu, lu_time};
use super::lu_list::lu_list_remove;
use super::lu_markowitz::lu_markowitz;
use super::lu_pivot::lu_pivot;
use super::*;

/// HiGHS: basiclu/lu_factorize_bump.c::lu_factorize_bump (lines 13-56)
pub fn lu_factorize_bump(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let mut status = BASICLU_OK;

    while this.rank + this.rankdef < m {
        // Find pivot element. Markowitz search need not be called if the
        // previous call to lu_pivot() returned for reallocation. In this case
        // this->pivot_col is valid.
        if this.pivot_col < 0 {
            lu_markowitz(this);
        }
        debug_assert!(this.pivot_col >= 0);

        if this.pivot_row < 0 {
            // Eliminate empty column without choosing a pivot.
            lu_list_remove(
                &mut this.colcount_flink[..],
                &mut this.colcount_blink[..],
                this.pivot_col,
            );
            this.pivot_col = -1;
            this.rankdef += 1;
        } else {
            // Eliminate pivot. This may require reallocation.
            debug_assert!(this.pinv[this.pivot_row as usize] == -1);
            debug_assert!(this.qinv[this.pivot_col as usize] == -1);
            status = lu_pivot(this);
            if status != BASICLU_OK {
                break;
            }
            this.pinv[this.pivot_row as usize] = this.rank;
            this.qinv[this.pivot_col as usize] = this.rank;
            this.pivot_col = -1;
            this.pivot_row = -1;
            this.rank += 1;
        }
        // rusolve BASICLU time limit (not in the HiGHS reference): check for time-out every 20 passes
        if (this.rank + this.rankdef) % 20 == 0
            && lu_time() - this.basiclu_time_start > this.basiclu_time_limit
        {
            return BASICLU_WARNING_timeout;
        }
    }
    status
}
