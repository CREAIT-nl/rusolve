//! HiGHS source: ipm/basiclu/lu_markowitz.c
//!
//! Copyright (C) 2016-2019  ERGO-Code (BASICLU source)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::lu_list::lu_list_move;
use super::*;

/// HiGHS: basiclu/lu_markowitz.c::lu_markowitz (lines 37-188)
///
/// PORT NOTE: the C `goto done` is `break 'done` on the outer labeled block.
/// The C `continue` of the column loop (which still executes the list step
/// `j = colcount_flink[j]`) is `break 'col` on a labeled block inside the loop.
/// The C `int_least64_t` Markowitz integers are `i64`.
#[allow(clippy::too_many_lines)]
pub fn lu_markowitz(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let mu = m as usize;
    let wbegin = &this.wbegin[..];
    let wend = &this.wend[..];
    let windex = &this.windex[..];
    let wvalue = &this.wvalue[..];
    let colcount_flink = &this.colcount_flink[..];
    let colmax = &this.col_pivot[..];
    let abstol = this.abstol;
    let reltol = this.reltol;
    let maxsearch = this.maxsearch;
    let search_rows = this.search_rows;
    let nz_start: LuInt = if search_rows != 0 {
        this.min_colnz.min(this.min_rownz)
    } else {
        this.min_colnz
    };
    let rowcount_flink = &mut this.rowcount_flink[..];
    let rowcount_blink = &mut this.rowcount_blink[..];

    // integers for Markowitz cost must be 64 bit to prevent overflow
    let m64 = m as i64;

    let mut pivot_row: LuInt = -1; // row of best pivot so far
    let mut pivot_col: LuInt = -1; // col of best pivot so far
    let mut mc_best: i64 = m64 * m64; // Markowitz cost of best pivot so far
    let mut nsearch: LuInt = 0; // count rows/columns searched
    let mut min_colnz: LuInt = -1; // minimum col count in active submatrix
    let mut min_rownz: LuInt = -1; // minimum row count in active submatrix
    debug_assert!(nz_start >= 1);

    'done: {
        // If the active submatrix contains empty columns, choose one and return
        // with pivot_row = -1.
        if colcount_flink[mu] != m {
            pivot_col = colcount_flink[mu];
            debug_assert!(pivot_col >= 0 && pivot_col < m);
            debug_assert!(wend[pivot_col as usize] == wbegin[pivot_col as usize]);
            break 'done;
        }

        for nz in nz_start..=m {
            // Search columns with nz nonzeros.
            let mut j = colcount_flink[(m + nz) as usize];
            while j < m {
                'col: {
                    if min_colnz == -1 {
                        min_colnz = nz;
                    }
                    let ju = j as usize;
                    debug_assert!(wend[ju] - wbegin[ju] == nz);
                    let cmx = colmax[ju];
                    debug_assert!(cmx >= 0.0);
                    if cmx == 0.0 || cmx < abstol {
                        break 'col;
                    }
                    let tol = abstol.max(reltol * cmx);
                    for pos in wbegin[ju]..wend[ju] {
                        let pu = pos as usize;
                        let x = wvalue[pu].abs();
                        if x == 0.0 || x < tol {
                            continue;
                        }
                        let i = windex[pu];
                        debug_assert!(i >= 0 && i < m);
                        let nz1 = nz as i64;
                        let nz2 = (wend[mu + i as usize] - wbegin[mu + i as usize]) as i64;
                        debug_assert!(nz2 >= 1);
                        let mc = (nz1 - 1) * (nz2 - 1);
                        if mc < mc_best {
                            mc_best = mc;
                            pivot_row = i;
                            pivot_col = j;
                            if search_rows != 0 && mc_best <= (nz1 - 1) * (nz1 - 1) {
                                break 'done;
                            }
                        }
                    }
                    // We have seen at least one eligible pivot in column j.
                    debug_assert!(mc_best < m64 * m64);
                    nsearch += 1;
                    if nsearch >= maxsearch {
                        break 'done;
                    }
                }
                j = colcount_flink[j as usize];
            }
            debug_assert!(j == m + nz);

            if search_rows == 0 {
                continue;
            }

            // Search rows with nz nonzeros.
            let mut i = rowcount_flink[(m + nz) as usize];
            while i < m {
                if min_rownz == -1 {
                    min_rownz = nz;
                }
                // rowcount_flink[i] might be changed below, so keep a copy
                let inext = rowcount_flink[i as usize];
                let iu = i as usize;
                debug_assert!(wend[mu + iu] - wbegin[mu + iu] == nz);
                let mut cheap = false; // row has entries with Markowitz cost < MC?
                let mut found = false; // eligible pivot found?
                for pos in wbegin[mu + iu]..wend[mu + iu] {
                    let pu = pos as usize;
                    let jj = windex[pu];
                    let jju = jj as usize;
                    debug_assert!(jj >= 0 && jj < m);
                    let nz1 = nz as i64;
                    let nz2 = (wend[jju] - wbegin[jju]) as i64;
                    debug_assert!(nz2 >= 1);
                    let mc = (nz1 - 1) * (nz2 - 1);
                    if mc >= mc_best {
                        continue;
                    }
                    cheap = true;
                    let cmx = colmax[jju];
                    debug_assert!(cmx >= 0.0);
                    if cmx == 0.0 || cmx < abstol {
                        continue;
                    }
                    // find position of pivot in column file
                    let mut wh = wbegin[jju];
                    while windex[wh as usize] != i {
                        debug_assert!(wh < wend[jju] - 1);
                        wh += 1;
                    }
                    let x = wvalue[wh as usize].abs();
                    if x >= abstol && x >= reltol * cmx {
                        found = true;
                        mc_best = mc;
                        pivot_row = i;
                        pivot_col = jj;
                        if mc_best <= nz1 * (nz1 - 1) {
                            break 'done;
                        }
                    }
                }
                // If row i has cheap entries but none of them is numerically
                // acceptable, then don't search the row again until updated.
                if cheap && !found {
                    lu_list_move(
                        i,
                        m + 1,
                        &mut *rowcount_flink,
                        &mut *rowcount_blink,
                        m,
                        None,
                    );
                } else {
                    debug_assert!(mc_best < m64 * m64);
                    nsearch += 1;
                    if nsearch >= maxsearch {
                        break 'done;
                    }
                }
                i = inext;
            }
            debug_assert!(i == m + nz);
        }
    }

    // done:
    this.pivot_row = pivot_row;
    this.pivot_col = pivot_col;
    this.nsearch_pivot += nsearch;
    if min_colnz >= 0 {
        this.min_colnz = min_colnz;
    }
    if min_rownz >= 0 {
        this.min_rownz = min_rownz;
    }
    BASICLU_OK
}
