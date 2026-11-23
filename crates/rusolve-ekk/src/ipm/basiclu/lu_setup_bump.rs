//! HiGHS source: ipm/basiclu/lu_setup_bump.c
//!
//! Copyright (C) 2016-2019  ERGO-Code (BASICLU source)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_file::{lu_file_diff, lu_file_empty};
use super::lu_internal::{Lu, lu_time};
use super::lu_list::{lu_list_add, lu_list_init, lu_list_move};
use super::*;

/// HiGHS: basiclu/lu_setup_bump.c::lu_setup_bump (lines 53-198)
///
/// PORT NOTE: `Wbegin2 = Wbegin + m` and `Wend2 = Wend + m` are the row half of
/// the W file, written as `wbegin[m + i]` and `wend[m + i]`. `colmax` is
/// `col_pivot`, `iwork0` is `iwork0`.
///
/// BASICLU quirk: `put += stretch*cnz + pad` adds a double to the lu_int `put`
/// and truncates the sum to an integer (C `put = (lu_int)((double)put + x)`).
#[allow(clippy::needless_range_loop)]
pub fn lu_setup_bump(
    this: &mut Lu<'_>,
    bbegin: &[LuInt],
    bend: &[LuInt],
    bi: &[LuInt],
    bx: &[f64],
) -> LuInt {
    let m = this.m;
    let mu = m as usize;
    let rank = this.rank;
    let wmem = this.wmem;
    let bnz = this.matrix_nz;
    let lnz = this.lbegin_p[rank as usize] - rank;
    let unz = this.ubegin[rank as usize];
    let abstol = this.abstol;
    let pad = this.pad;
    let stretch = this.stretch;

    let mut bump_nz: LuInt = bnz - lnz - unz - rank; // will change if columns are dropped
    let mut min_colnz: LuInt = 0;
    let mut min_rownz: LuInt = 0;

    debug_assert!(lnz >= 0);
    debug_assert!(unz >= 0);
    debug_assert!(bump_nz >= 0);
    for i in 0..mu {
        debug_assert!(this.iwork0[i] == 0);
    }

    // Calculate memory and reallocate. For each row/column with nz nonzeros
    // add stretch*nz+pad elements extra space for fill-in.
    let b = bump_nz as f64;
    let need_f: f64 = b + stretch * b + ((m - rank) * pad) as f64;
    let mut need: LuInt = need_f as LuInt;
    need = 2 * need; // rowwise + columnwise
    if need > wmem {
        this.addmem_w = need - wmem;
        return BASICLU_REALLOCATE;
    }

    lu_file_empty(
        2 * m,
        &mut this.wbegin[..],
        &mut this.wend[..],
        &mut this.wflink[..],
        &mut this.wblink[..],
        wmem,
    );

    // Build columnwise storage. Build row counts in iwork0.
    lu_list_init(
        &mut this.colcount_flink[..],
        &mut this.colcount_blink[..],
        m,
        m + 2,
        Some(&mut min_colnz),
    );
    let mut put: LuInt = 0;
    for j in 0..mu {
        if this.qinv[j] >= 0 {
            continue;
        }
        let mut cnz: LuInt = 0; // count nz per column
        let mut cmx: f64 = 0.0; // find column maximum
        for pos in bbegin[j]..bend[j] {
            let i = bi[pos as usize];
            if this.pinv[i as usize] >= 0 {
                continue;
            }
            cmx = cmx.max(bx[pos as usize].abs());
            cnz += 1;
        }
        if cmx == 0.0 || cmx < abstol {
            // Leave column of active submatrix empty.
            this.col_pivot[j] = 0.0;
            lu_list_add(
                j as LuInt,
                0,
                &mut this.colcount_flink[..],
                &mut this.colcount_blink[..],
                m,
                Some(&mut min_colnz),
            );
            bump_nz -= cnz;
        } else {
            // Copy column into active submatrix.
            this.col_pivot[j] = cmx;
            lu_list_add(
                j as LuInt,
                cnz,
                &mut this.colcount_flink[..],
                &mut this.colcount_blink[..],
                m,
                Some(&mut min_colnz),
            );
            this.wbegin[j] = put;
            for pos in bbegin[j]..bend[j] {
                let i = bi[pos as usize];
                if this.pinv[i as usize] >= 0 {
                    continue;
                }
                this.windex[put as usize] = i;
                this.wvalue[put as usize] = bx[pos as usize];
                put += 1;
                this.iwork0[i as usize] += 1;
            }
            this.wend[j] = put;
            // BASICLU quirk: put += stretch*cnz + pad truncates to lu_int.
            put = ((put as f64) + (stretch * (cnz as f64) + (pad as f64))) as LuInt;
            // reappend line to list end
            lu_list_move(
                j as LuInt,
                0,
                &mut this.wflink[..],
                &mut this.wblink[..],
                2 * m,
                None,
            );
        }
    }

    // Build rowwise storage (pattern only).
    lu_list_init(
        &mut this.rowcount_flink[..],
        &mut this.rowcount_blink[..],
        m,
        m + 2,
        Some(&mut min_rownz),
    );
    for i in 0..mu {
        // set row pointers
        if this.pinv[i] >= 0 {
            continue;
        }
        let rnz = this.iwork0[i];
        this.iwork0[i] = 0;
        lu_list_add(
            i as LuInt,
            rnz,
            &mut this.rowcount_flink[..],
            &mut this.rowcount_blink[..],
            m,
            Some(&mut min_rownz),
        );
        this.wbegin[mu + i] = put;
        this.wend[mu + i] = put;
        put += rnz;
        // reappend line to list end
        lu_list_move(
            m + i as LuInt,
            0,
            &mut this.wflink[..],
            &mut this.wblink[..],
            2 * m,
            None,
        );
        // BASICLU quirk: put += stretch*rnz + pad truncates to lu_int.
        put = ((put as f64) + (stretch * (rnz as f64) + (pad as f64))) as LuInt;
    }
    for j in 0..mu {
        // fill rows
        for pos in this.wbegin[j]..this.wend[j] {
            let i = this.windex[pos as usize] as usize;
            let p = this.wend[mu + i];
            this.wend[mu + i] += 1; // Windex[Wend2[i]++] = j
            this.windex[p as usize] = j as LuInt;
        }
    }
    this.wbegin[2 * mu] = put; // set beginning of free space
    debug_assert!(this.wbegin[2 * mu] <= this.wend[2 * mu]);

    debug_assert_eq!(
        lu_file_diff(
            m,
            &this.wbegin[..],
            &this.wend[..],
            &this.wbegin[mu..],
            &this.wend[mu..],
            &this.windex[..],
            None,
        ),
        0
    );
    debug_assert_eq!(
        lu_file_diff(
            m,
            &this.wbegin[mu..],
            &this.wend[mu..],
            &this.wbegin[..],
            &this.wend[..],
            &this.windex[..],
            None,
        ),
        0
    );

    for i in 0..mu {
        debug_assert!(this.iwork0[i] == 0);
    }

    this.bump_nz = bump_nz;
    this.bump_size = m - rank;
    this.min_colnz = min_colnz;
    this.min_rownz = min_rownz;
    if lu_time() - this.basiclu_time_start > this.basiclu_time_limit {
        return BASICLU_WARNING_timeout;
    }
    BASICLU_OK
}
