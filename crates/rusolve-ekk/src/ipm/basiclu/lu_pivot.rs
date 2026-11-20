//! HiGHS source: ipm/basiclu/lu_pivot.c
//!
//! Copyright (C) 2016-2019 ERGO-Code (BASICLU, lu_pivot.c)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `DEBUG_EXTRA` consistency checks (`lu_file_diff` after each pivot) are not
//! compiled: the HiGHS build does not define `DEBUG_EXTRA`.

use super::lu_file::{lu_file_compress, lu_file_reappend};
use super::lu_internal::Lu;
use super::lu_list::{lu_list_move, lu_list_remove};
use super::*;

/// HiGHS: basiclu/lu_pivot.c `#define MAXROW_SMALL 64`
const MAXROW_SMALL: LuInt = 64;

/// The double expression `nz+r + stretch*(nz+r) + pad` of lu_pivot.c, with
/// `s = nz+r` already formed as an integer (C evaluates the `nz+r` terms
/// first, as int, then converts). Evaluated left to right in double.
#[inline]
fn fill_room(s: LuInt, stretch: f64, pad: LuInt) -> f64 {
    (s as f64) + stretch * (s as f64) + (pad as f64)
}

/// HiGHS: basiclu/lu_pivot.c::lu_pivot (lines 63-141)
pub fn lu_pivot(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let rank = this.rank;
    let pivot_col = this.pivot_col;
    let pivot_row = this.pivot_row;
    let nz_col = this.wend[pivot_col as usize] - this.wbegin[pivot_col as usize];
    let nz_row = this.wend[(m + pivot_row) as usize] - this.wbegin[(m + pivot_row) as usize];

    debug_assert!(nz_row >= 1);
    debug_assert!(nz_col >= 1);

    // Check if room is available in L and U.
    let mut status = BASICLU_OK;
    let mut room = this.lmem - this.lbegin_p[rank as usize];
    let mut need = nz_col; // # off-diagonals in pivot col + end marker (-1)
    if room < need {
        this.addmem_l = need - room;
        status = BASICLU_REALLOCATE;
    }
    room = this.umem - this.ubegin[rank as usize];
    need = nz_row - 1; // # off-diagonals in pivot row
    if room < need {
        this.addmem_u = need - room;
        status = BASICLU_REALLOCATE;
    }
    if status != BASICLU_OK {
        return status;
    }

    // Branch out implementation of pivot operation.
    status = if nz_row == 1 {
        lu_pivot_singleton_row(this)
    } else if nz_col == 1 {
        lu_pivot_singleton_col(this)
    } else if nz_col == 2 {
        lu_pivot_doubleton_col(this)
    } else if nz_col - 1 <= MAXROW_SMALL {
        lu_pivot_small(this)
    } else {
        lu_pivot_any(this)
    };

    // Remove all entries in columns whose maximum entry has dropped below
    // absolute pivot tolerance.
    if status == BASICLU_OK {
        for pos in this.ubegin[rank as usize]..this.ubegin[rank as usize + 1] {
            let j = this.uindex[pos as usize];
            debug_assert!(j != pivot_col);
            let colmax = this.col_pivot[j as usize];
            if colmax == 0.0 || colmax < this.abstol {
                lu_remove_col(this, j);
            }
        }
    }

    // BASICLU quirk: C `factor_flops` is a signed int and the product may
    // overflow; C leaves that undefined, the arm64 build wraps.
    this.factor_flops = this
        .factor_flops
        .wrapping_add((nz_col - 1).wrapping_mul(nz_row - 1));
    status
}

/// HiGHS: basiclu/lu_pivot.c::lu_pivot_any (lines 148-432)
fn lu_pivot_any(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let rank = this.rank;
    let droptol = this.droptol;
    let pad = this.pad;
    let stretch = this.stretch;
    let pivot_col = this.pivot_col;
    let pivot_row = this.pivot_row;
    let pc = pivot_col as usize;
    let prow_m = (m + pivot_row) as usize; // row head index in Wbegin/Wend
    let m2 = (2 * m) as usize; // head of the free-space list

    let mut cbeg = this.wbegin[pc]; // changed by file compression
    let mut cend = this.wend[pc];
    let mut rbeg = this.wbegin[prow_m];
    let mut rend = this.wend[prow_m];
    let cnz1 = cend - cbeg - 1; // nz in pivot column except pivot
    let rnz1 = rend - rbeg - 1; // nz in pivot row except pivot

    // Check if room is available in W. At most each updated row and each
    // updated column will be reappended and filled-in with rnz1 respectively
    // cnz1 elements. Move pivot to the front of pivot row and pivot column.
    let mut grow: LuInt = 0;
    let mut where_: LuInt = -1;
    for pos in cbeg..cend {
        let i = this.windex[pos as usize];
        if i == pivot_row {
            where_ = pos;
        } else {
            let nz = this.wend[(m + i) as usize] - this.wbegin[(m + i) as usize];
            // BASICLU quirk: `grow += <double>` truncates the double sum to lu_int.
            grow = ((grow as f64) + fill_room(nz + rnz1, stretch, pad)) as LuInt;
        }
    }
    debug_assert!(where_ >= 0);
    lu_iswap(this.windex, cbeg as usize, where_ as usize);
    lu_fswap(this.wvalue, cbeg as usize, where_ as usize);
    let pivot = this.wvalue[cbeg as usize];
    debug_assert!(pivot != 0.0);
    where_ = -1;
    for rpos in rbeg..rend {
        let j = this.windex[rpos as usize];
        if j == pivot_col {
            where_ = rpos;
        } else {
            let nz = this.wend[j as usize] - this.wbegin[j as usize];
            grow = ((grow as f64) + fill_room(nz + cnz1, stretch, pad)) as LuInt;
        }
    }
    debug_assert!(where_ >= 0);
    lu_iswap(this.windex, rbeg as usize, where_ as usize);
    let mut room = this.wend[m2] - this.wbegin[m2];
    if grow > room {
        lu_file_compress(
            2 * m,
            this.wbegin,
            this.wend,
            this.wflink,
            this.windex,
            this.wvalue,
            stretch,
            pad,
        );
        cbeg = this.wbegin[pc];
        cend = this.wend[pc];
        rbeg = this.wbegin[prow_m];
        rend = this.wend[prow_m];
        room = this.wend[m2] - this.wbegin[m2];
        this.ngarbage += 1;
    }
    if grow > room {
        this.addmem_w = grow - room;
        return BASICLU_REALLOCATE;
    }

    // get pointer to U
    let mut uput = this.ubegin[rank as usize];
    debug_assert!(uput >= 0);
    debug_assert!(uput < this.umem);

    // ----------------------------------------------------------------------
    // Column file update
    // ----------------------------------------------------------------------

    // For each row i to be updated set marked[i] > 0 to its position
    // in the (packed) pivot column.
    // position of pos in the packed column is pos - cbeg (position starts at 1)
    for pos in (cbeg + 1)..cend {
        let i = this.windex[pos as usize];
        this.iwork0[i as usize] = pos - cbeg;
    }

    // C: wi = Windex + cbeg, wx = Wvalue + cbeg (read live, absolute indices)
    for rpos in (rbeg + 1)..rend {
        let j = this.windex[rpos as usize];
        debug_assert!(j != pivot_col);
        let mut cmx: f64 = 0.0; // column maximum

        // Compress unmodified column entries. Store entries to be updated
        // in workspace. Move pivot row entry to the front of column.
        let mut where_: LuInt = -1;
        let pos1 = this.wbegin[j as usize];
        let mut put = pos1;
        let wend_j = this.wend[j as usize];
        for pos in pos1..wend_j {
            let i = this.windex[pos as usize];
            let position = this.iwork0[i as usize];
            if position > 0 {
                debug_assert!(i != pivot_row);
                this.work0[position as usize] = this.wvalue[pos as usize];
            } else {
                debug_assert!(position == 0);
                if i == pivot_row {
                    where_ = put;
                } else {
                    let x = this.wvalue[pos as usize].abs();
                    if x > cmx {
                        cmx = x;
                    }
                }
                this.windex[put as usize] = this.windex[pos as usize];
                this.wvalue[put as usize] = this.wvalue[pos as usize];
                put += 1;
            }
        }
        debug_assert!(where_ >= 0);
        this.wend[j as usize] = put;
        lu_iswap(this.windex, pos1 as usize, where_ as usize);
        lu_fswap(this.wvalue, pos1 as usize, where_ as usize);
        let xrj = this.wvalue[pos1 as usize]; // pivot row entry

        // Reappend column if no room for update.
        let mut room = this.wbegin[this.wflink[j as usize] as usize] - put;
        if room < cnz1 {
            let nz = this.wend[j as usize] - this.wbegin[j as usize];
            room = (cnz1 as f64 + stretch * ((nz + cnz1) as f64) + (pad as f64)) as LuInt;
            lu_file_reappend(
                j,
                2 * m,
                this.wbegin,
                this.wend,
                this.wflink,
                this.wblink,
                this.windex,
                this.wvalue,
                room,
            );
            put = this.wend[j as usize];
            debug_assert!(this.wbegin[this.wflink[j as usize] as usize] - put == room);
            this.nexpand += 1;
        }

        // Compute update in workspace and append to column.
        let a = xrj / pivot;
        for pos in 1..=cnz1 {
            let wx_pos = this.wvalue[(cbeg + pos) as usize];
            this.work0[pos as usize] -= a * wx_pos;
        }
        for pos in 1..=cnz1 {
            this.windex[put as usize] = this.windex[(cbeg + pos) as usize];
            this.wvalue[put as usize] = this.work0[pos as usize];
            put += 1;
            let x = this.work0[pos as usize].abs();
            if x > cmx {
                cmx = x;
            }
            this.work0[pos as usize] = 0.0;
        }
        this.wend[j as usize] = put;

        // Write pivot row entry to U and remove from file.
        if xrj.abs() > droptol {
            debug_assert!(uput < this.umem);
            this.uindex[uput as usize] = j;
            this.uvalue[uput as usize] = xrj;
            uput += 1;
        }
        debug_assert!(this.windex[this.wbegin[j as usize] as usize] == pivot_row);
        this.wbegin[j as usize] += 1;

        // Move column to new list and update min_colnz.
        let nz = this.wend[j as usize] - this.wbegin[j as usize];
        lu_list_move(
            j,
            nz,
            this.colcount_flink,
            this.colcount_blink,
            m,
            Some(&mut this.min_colnz),
        );

        this.col_pivot[j as usize] = cmx;
    }
    for pos in (cbeg + 1)..cend {
        let i = this.windex[pos as usize];
        this.iwork0[i as usize] = 0;
    }

    // ----------------------------------------------------------------------
    // Row file update
    // ----------------------------------------------------------------------

    for rpos in rbeg..rend {
        let i = this.windex[rpos as usize];
        this.iwork0[i as usize] = 1;
    }
    debug_assert!(this.iwork0[pc] == 1);

    for pos in (cbeg + 1)..cend {
        let i = this.windex[pos as usize];
        let mi = (m + i) as usize;
        debug_assert!(i != pivot_row);

        // Compress unmodified row entries (not marked). Remove
        // overlap with pivot row, including pivot column entry.
        let mut found = false;
        let mut put = this.wbegin[mi];
        for rp in this.wbegin[mi]..this.wend[mi] {
            let j = this.windex[rp as usize];
            if j == pivot_col {
                found = true;
            }
            if this.iwork0[j as usize] == 0 {
                this.windex[put as usize] = j;
                put += 1;
            }
        }
        debug_assert!(found);
        this.wend[mi] = put;

        // Reappend row if no room for update. Append pattern of pivot row.
        let mut room = this.wbegin[this.wflink[mi] as usize] - put;
        if room < rnz1 {
            let nz = this.wend[mi] - this.wbegin[mi];
            room = (rnz1 as f64 + stretch * ((nz + rnz1) as f64) + (pad as f64)) as LuInt;
            lu_file_reappend(
                m + i,
                2 * m,
                this.wbegin,
                this.wend,
                this.wflink,
                this.wblink,
                this.windex,
                this.wvalue,
                room,
            );
            put = this.wend[mi];
            debug_assert!(this.wbegin[this.wflink[mi] as usize] - put == room);
            this.nexpand += 1;
        }
        for rp in (rbeg + 1)..rend {
            this.windex[put as usize] = this.windex[rp as usize];
            put += 1;
        }
        this.wend[mi] = put;

        // Move to new list. The row must be reinserted even if nz are
        // unchanged since it might have been taken out in Markowitz search.
        let nz = this.wend[mi] - this.wbegin[mi];
        lu_list_move(
            i,
            nz,
            this.rowcount_flink,
            this.rowcount_blink,
            m,
            Some(&mut this.min_rownz),
        );
    }
    for rpos in rbeg..rend {
        let i = this.windex[rpos as usize];
        this.iwork0[i as usize] = 0;
    }

    // Store column in L.
    let mut put = this.lbegin_p[rank as usize];
    for pos in (cbeg + 1)..cend {
        let x = this.wvalue[pos as usize] / pivot;
        if x.abs() > droptol {
            this.lindex[put as usize] = this.windex[pos as usize];
            this.lvalue[put as usize] = x;
            put += 1;
        }
    }
    this.lindex[put as usize] = -1; // terminate column
    put += 1;
    this.lbegin_p[rank as usize + 1] = put;
    this.ubegin[rank as usize + 1] = uput;

    // Cleanup:
    // store pivot element;
    // remove pivot column from column file, pivot row from row file;
    // remove pivot column/row from column/row counts
    this.col_pivot[pc] = pivot;
    this.wend[pc] = cbeg;
    this.wend[prow_m] = rbeg;
    lu_list_remove(this.colcount_flink, this.colcount_blink, pivot_col);
    lu_list_remove(this.rowcount_flink, this.rowcount_blink, pivot_row);

    BASICLU_OK
}

/// HiGHS: basiclu/lu_pivot.c::lu_pivot_small (lines 439-748)
///
/// `int64_t *cancelled = (void *) this->row_pivot`: one `i64` per
/// `row_pivot` slot through [`pun_set_int`]/[`pun_get_int`].
///
/// PORT NOTE: the C `int64_t` array overlays `row_pivot`, which is `m`
/// doubles, and only `rnz1 <= m` entries are written. Nothing reads
/// `row_pivot` as doubles inside this routine, and `row_pivot` is only
/// meaningful again after `lu_pivot` returns (the stale bits are never read as
/// `f64` here). The Rust pun keeps one `i64` per slot, so the layout differs
/// from C only in the unused halves, which C never touches.
fn lu_pivot_small(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let rank = this.rank;
    let droptol = this.droptol;
    let pad = this.pad;
    let stretch = this.stretch;
    let pivot_col = this.pivot_col;
    let pivot_row = this.pivot_row;
    let pc = pivot_col as usize;
    let prow_m = (m + pivot_row) as usize;
    let m2 = (2 * m) as usize;

    let mut cbeg = this.wbegin[pc]; // changed by file compression
    let mut cend = this.wend[pc];
    let mut rbeg = this.wbegin[prow_m];
    let mut rend = this.wend[prow_m];
    let cnz1 = cend - cbeg - 1; // nz in pivot column except pivot
    let rnz1 = rend - rbeg - 1; // nz in pivot row except pivot

    debug_assert!(cnz1 <= MAXROW_SMALL);

    // Check if room is available in W (see lu_pivot_any).
    let mut grow: LuInt = 0;
    let mut where_: LuInt = -1;
    for pos in cbeg..cend {
        let i = this.windex[pos as usize];
        if i == pivot_row {
            where_ = pos;
        } else {
            let nz = this.wend[(m + i) as usize] - this.wbegin[(m + i) as usize];
            grow = ((grow as f64) + fill_room(nz + rnz1, stretch, pad)) as LuInt;
        }
    }
    debug_assert!(where_ >= 0);
    lu_iswap(this.windex, cbeg as usize, where_ as usize);
    lu_fswap(this.wvalue, cbeg as usize, where_ as usize);
    let pivot = this.wvalue[cbeg as usize];
    debug_assert!(pivot != 0.0);
    where_ = -1;
    for rpos in rbeg..rend {
        let j = this.windex[rpos as usize];
        if j == pivot_col {
            where_ = rpos;
        } else {
            let nz = this.wend[j as usize] - this.wbegin[j as usize];
            grow = ((grow as f64) + fill_room(nz + cnz1, stretch, pad)) as LuInt;
        }
    }
    debug_assert!(where_ >= 0);
    lu_iswap(this.windex, rbeg as usize, where_ as usize);
    let mut room = this.wend[m2] - this.wbegin[m2];
    if grow > room {
        lu_file_compress(
            2 * m,
            this.wbegin,
            this.wend,
            this.wflink,
            this.windex,
            this.wvalue,
            stretch,
            pad,
        );
        cbeg = this.wbegin[pc];
        cend = this.wend[pc];
        rbeg = this.wbegin[prow_m];
        rend = this.wend[prow_m];
        room = this.wend[m2] - this.wbegin[m2];
        this.ngarbage += 1;
    }
    if grow > room {
        this.addmem_w = grow - room;
        return BASICLU_REALLOCATE;
    }

    // get pointer to U
    let mut uput = this.ubegin[rank as usize];
    debug_assert!(uput >= 0);
    debug_assert!(uput < this.umem);

    // ----------------------------------------------------------------------
    // Column file update
    // ----------------------------------------------------------------------

    // position of pos in the packed column is pos - cbeg (position starts at 1)
    for pos in (cbeg + 1)..cend {
        let i = this.windex[pos as usize];
        this.iwork0[i as usize] = pos - cbeg;
    }

    // C: col_number indexes `cancelled`, one bit mask per pivot-row entry
    for rpos in (rbeg + 1)..rend {
        let col_number = rpos - (rbeg + 1); // C: col_number++ alongside rpos
        let j = this.windex[rpos as usize];
        debug_assert!(j != pivot_col);
        let mut cmx: f64 = 0.0; // column maximum

        // Compress unmodified column entries. Store entries to be updated
        // in workspace. Move pivot row entry to the front of column.
        let mut where_: LuInt = -1;
        let pos1 = this.wbegin[j as usize];
        let mut put = pos1;
        let wend_j = this.wend[j as usize];
        for pos in pos1..wend_j {
            let i = this.windex[pos as usize];
            let position = this.iwork0[i as usize];
            if position > 0 {
                debug_assert!(i != pivot_row);
                this.work0[position as usize] = this.wvalue[pos as usize];
            } else {
                debug_assert!(position == 0);
                if i == pivot_row {
                    where_ = put;
                } else {
                    let x = this.wvalue[pos as usize].abs();
                    if x > cmx {
                        cmx = x;
                    }
                }
                this.windex[put as usize] = this.windex[pos as usize];
                this.wvalue[put as usize] = this.wvalue[pos as usize];
                put += 1;
            }
        }
        debug_assert!(where_ >= 0);
        this.wend[j as usize] = put;
        lu_iswap(this.windex, pos1 as usize, where_ as usize);
        lu_fswap(this.wvalue, pos1 as usize, where_ as usize);
        let xrj = this.wvalue[pos1 as usize]; // pivot row entry

        // Reappend column if no room for update.
        let mut room = this.wbegin[this.wflink[j as usize] as usize] - put;
        if room < cnz1 {
            let nz = this.wend[j as usize] - this.wbegin[j as usize];
            room = (cnz1 as f64 + stretch * ((nz + cnz1) as f64) + (pad as f64)) as LuInt;
            lu_file_reappend(
                j,
                2 * m,
                this.wbegin,
                this.wend,
                this.wflink,
                this.wblink,
                this.windex,
                this.wvalue,
                room,
            );
            put = this.wend[j as usize];
            debug_assert!(this.wbegin[this.wflink[j as usize] as usize] - put == room);
            this.nexpand += 1;
        }

        // Compute update in workspace and append to column.
        let a = xrj / pivot;
        for pos in 1..=cnz1 {
            let wx_pos = this.wvalue[(cbeg + pos) as usize];
            this.work0[pos as usize] -= a * wx_pos;
        }
        let mut mask: i64 = 0;
        for pos in 1..=cnz1 {
            let x = this.work0[pos as usize].abs();
            if x > droptol {
                this.windex[put as usize] = this.windex[(cbeg + pos) as usize];
                this.wvalue[put as usize] = this.work0[pos as usize];
                put += 1;
                if x > cmx {
                    cmx = x;
                }
            } else {
                // cancellation in row wi[pos]
                mask |= 1i64 << (pos - 1);
            }
            this.work0[pos as usize] = 0.0;
        }
        this.wend[j as usize] = put;
        pun_set_int(this.row_pivot, col_number as usize, mask);

        // Write pivot row entry to U and remove from file.
        if xrj.abs() > droptol {
            debug_assert!(uput < this.umem);
            this.uindex[uput as usize] = j;
            this.uvalue[uput as usize] = xrj;
            uput += 1;
        }
        debug_assert!(this.windex[this.wbegin[j as usize] as usize] == pivot_row);
        this.wbegin[j as usize] += 1;

        // Move column to new list and update min_colnz.
        let nz = this.wend[j as usize] - this.wbegin[j as usize];
        lu_list_move(
            j,
            nz,
            this.colcount_flink,
            this.colcount_blink,
            m,
            Some(&mut this.min_colnz),
        );

        this.col_pivot[j as usize] = cmx;
    }
    for pos in (cbeg + 1)..cend {
        let i = this.windex[pos as usize];
        this.iwork0[i as usize] = 0;
    }

    // ----------------------------------------------------------------------
    // Row file update
    // ----------------------------------------------------------------------

    for rpos in rbeg..rend {
        let i = this.windex[rpos as usize];
        this.iwork0[i as usize] = 1;
    }
    debug_assert!(this.iwork0[pc] == 1);

    let mut mask: i64 = 1;
    for pos in (cbeg + 1)..cend {
        debug_assert!(mask != 0);
        let i = this.windex[pos as usize];
        let mi = (m + i) as usize;
        debug_assert!(i != pivot_row);

        // Compress unmodified row entries (not marked). Remove
        // overlap with pivot row, including pivot column entry.
        let mut found = false;
        let mut put = this.wbegin[mi];
        for rp in this.wbegin[mi]..this.wend[mi] {
            let j = this.windex[rp as usize];
            if j == pivot_col {
                found = true;
            }
            if this.iwork0[j as usize] == 0 {
                this.windex[put as usize] = j;
                put += 1;
            }
        }
        debug_assert!(found);
        this.wend[mi] = put;

        // Reappend row if no room for update. Append pattern of pivot row.
        let mut room = this.wbegin[this.wflink[mi] as usize] - put;
        if room < rnz1 {
            let nz = this.wend[mi] - this.wbegin[mi];
            room = (rnz1 as f64 + stretch * ((nz + rnz1) as f64) + (pad as f64)) as LuInt;
            lu_file_reappend(
                m + i,
                2 * m,
                this.wbegin,
                this.wend,
                this.wflink,
                this.wblink,
                this.windex,
                this.wvalue,
                room,
            );
            put = this.wend[mi];
            debug_assert!(this.wbegin[this.wflink[mi] as usize] - put == room);
            this.nexpand += 1;
        }
        for rp in (rbeg + 1)..rend {
            let col_number = rp - (rbeg + 1); // C: col_number++ alongside rp
            if (pun_get_int(this.row_pivot, col_number as usize) & mask) == 0 {
                this.windex[put as usize] = this.windex[rp as usize];
                put += 1;
            }
        }
        this.wend[mi] = put;

        // Move to new list. The row must be reinserted even if nz are
        // unchanged since it might have been taken out in Markowitz search.
        let nz = this.wend[mi] - this.wbegin[mi];
        lu_list_move(
            i,
            nz,
            this.rowcount_flink,
            this.rowcount_blink,
            m,
            Some(&mut this.min_rownz),
        );
        mask <<= 1;
    }
    for rpos in rbeg..rend {
        let i = this.windex[rpos as usize];
        this.iwork0[i as usize] = 0;
    }

    // Store column in L.
    let mut put = this.lbegin_p[rank as usize];
    for pos in (cbeg + 1)..cend {
        let x = this.wvalue[pos as usize] / pivot;
        if x.abs() > droptol {
            this.lindex[put as usize] = this.windex[pos as usize];
            this.lvalue[put as usize] = x;
            put += 1;
        }
    }
    this.lindex[put as usize] = -1; // terminate column
    put += 1;
    this.lbegin_p[rank as usize + 1] = put;
    this.ubegin[rank as usize + 1] = uput;

    // Cleanup (see lu_pivot_any).
    this.col_pivot[pc] = pivot;
    this.wend[pc] = cbeg;
    this.wend[prow_m] = rbeg;
    lu_list_remove(this.colcount_flink, this.colcount_blink, pivot_col);
    lu_list_remove(this.rowcount_flink, this.rowcount_blink, pivot_row);

    BASICLU_OK
}

/// HiGHS: basiclu/lu_pivot.c::lu_pivot_singleton_row (lines 755-835)
fn lu_pivot_singleton_row(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let rank = this.rank;
    let droptol = this.droptol;
    let pivot_col = this.pivot_col;
    let pivot_row = this.pivot_row;
    let pc = pivot_col as usize;
    let prow_m = (m + pivot_row) as usize;

    let cbeg = this.wbegin[pc];
    let cend = this.wend[pc];
    let rbeg = this.wbegin[prow_m];
    let rend = this.wend[prow_m];
    let rnz1 = rend - rbeg - 1; // nz in pivot row except pivot

    debug_assert!(rnz1 == 0);

    // Find pivot.
    let mut where_ = cbeg;
    while this.windex[where_ as usize] != pivot_row {
        debug_assert!(where_ < cend - 1);
        where_ += 1;
    }
    let pivot = this.wvalue[where_ as usize];
    debug_assert!(pivot != 0.0);

    // Store column in L.
    let mut put = this.lbegin_p[rank as usize];
    for pos in cbeg..cend {
        let x = this.wvalue[pos as usize] / pivot;
        if pos != where_ && x.abs() > droptol {
            this.lindex[put as usize] = this.windex[pos as usize];
            this.lvalue[put as usize] = x;
            put += 1;
        }
    }
    this.lindex[put as usize] = -1; // terminate column
    put += 1;
    this.lbegin_p[rank as usize + 1] = put;
    this.ubegin[rank as usize + 1] = this.ubegin[rank as usize];

    // Remove pivot column from row file. Update row lists.
    for pos in cbeg..cend {
        let i = this.windex[pos as usize];
        if i == pivot_row {
            continue;
        }
        let mi = (m + i) as usize;
        let mut w = this.wbegin[mi];
        while this.windex[w as usize] != pivot_col {
            debug_assert!(w < this.wend[mi] - 1);
            w += 1;
        }
        where_ = w;
        this.wend[mi] -= 1;
        this.windex[where_ as usize] = this.windex[this.wend[mi] as usize];
        let nz = this.wend[mi] - this.wbegin[mi];
        lu_list_move(
            i,
            nz,
            this.rowcount_flink,
            this.rowcount_blink,
            m,
            Some(&mut this.min_rownz),
        );
    }

    // Cleanup:
    // store pivot element;
    // remove pivot column from column file, pivot row from row file;
    // remove pivot column/row from column/row counts
    this.col_pivot[pc] = pivot;
    this.wend[pc] = cbeg;
    this.wend[prow_m] = rbeg;
    lu_list_remove(this.colcount_flink, this.colcount_blink, pivot_col);
    lu_list_remove(this.rowcount_flink, this.rowcount_blink, pivot_row);

    BASICLU_OK
}

/// HiGHS: basiclu/lu_pivot.c::lu_pivot_singleton_col (lines 842-935)
fn lu_pivot_singleton_col(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let rank = this.rank;
    let droptol = this.droptol;
    let pivot_col = this.pivot_col;
    let pivot_row = this.pivot_row;
    let pc = pivot_col as usize;
    let prow_m = (m + pivot_row) as usize;

    let cbeg = this.wbegin[pc];
    let cend = this.wend[pc];
    let rbeg = this.wbegin[prow_m];
    let rend = this.wend[prow_m];
    let cnz1 = cend - cbeg - 1; // nz in pivot column except pivot

    debug_assert!(cnz1 == 0);

    // Remove pivot row from column file and store in U. Update column lists.
    let mut put = this.ubegin[rank as usize];
    let pivot = this.wvalue[cbeg as usize];
    debug_assert!(pivot != 0.0);
    let mut found = false;
    let mut xrj: f64 = 0.0; // initialize to make gcc happy
    for rpos in rbeg..rend {
        let j = this.windex[rpos as usize];
        if j == pivot_col {
            found = true;
            continue;
        }
        let mut where_: LuInt = -1;
        let mut cmx: f64 = 0.0; // column maximum
        for pos in this.wbegin[j as usize]..this.wend[j as usize] {
            if this.windex[pos as usize] == pivot_row {
                where_ = pos;
                xrj = this.wvalue[pos as usize];
            } else {
                let x = this.wvalue[pos as usize].abs();
                if x > cmx {
                    cmx = x;
                }
            }
        }
        debug_assert!(where_ >= 0);
        if xrj.abs() > droptol {
            this.uindex[put as usize] = j;
            this.uvalue[put as usize] = xrj;
            put += 1;
        }
        let ju = j as usize;
        this.wend[ju] -= 1;
        let last = this.wend[ju] as usize;
        this.windex[where_ as usize] = this.windex[last];
        this.wvalue[where_ as usize] = this.wvalue[last];
        let nz = this.wend[ju] - this.wbegin[ju];
        lu_list_move(
            j,
            nz,
            this.colcount_flink,
            this.colcount_blink,
            m,
            Some(&mut this.min_colnz),
        );
        this.col_pivot[ju] = cmx;
    }
    debug_assert!(found);
    this.ubegin[rank as usize + 1] = put;

    // Store empty column in L.
    let mut put = this.lbegin_p[rank as usize];
    this.lindex[put as usize] = -1; // terminate column
    put += 1;
    this.lbegin_p[rank as usize + 1] = put;

    // Cleanup:
    // store pivot element;
    // remove pivot column from column file, pivot row from row file;
    // remove pivot column/row from column/row counts
    this.col_pivot[pc] = pivot;
    this.wend[pc] = cbeg;
    this.wend[prow_m] = rbeg;
    lu_list_remove(this.colcount_flink, this.colcount_blink, pivot_col);
    lu_list_remove(this.rowcount_flink, this.rowcount_blink, pivot_row);

    BASICLU_OK
}

/// HiGHS: basiclu/lu_pivot.c::lu_pivot_doubleton_col (lines 942-1208)
fn lu_pivot_doubleton_col(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let rank = this.rank;
    let droptol = this.droptol;
    let pad = this.pad;
    let stretch = this.stretch;
    let pivot_col = this.pivot_col;
    let pivot_row = this.pivot_row;
    let pc = pivot_col as usize;
    let prow_m = (m + pivot_row) as usize;
    let m2 = (2 * m) as usize;

    let mut cbeg = this.wbegin[pc]; // changed by file compression
    let mut rbeg = this.wbegin[prow_m];
    let mut rend = this.wend[prow_m];
    let cnz1 = this.wend[pc] - cbeg - 1; // nz in pivot column except pivot
    let rnz1 = rend - rbeg - 1; // nz in pivot row except pivot

    debug_assert!(cnz1 == 1);

    // Move pivot element to front of pivot column and pivot row.
    if this.windex[cbeg as usize] != pivot_row {
        lu_iswap(this.windex, cbeg as usize, (cbeg + 1) as usize);
        lu_fswap(this.wvalue, cbeg as usize, (cbeg + 1) as usize);
    }
    debug_assert!(this.windex[cbeg as usize] == pivot_row);
    let pivot = this.wvalue[cbeg as usize];
    debug_assert!(pivot != 0.0);
    let other_row = this.windex[(cbeg + 1) as usize];
    let other_value = this.wvalue[(cbeg + 1) as usize];
    let mut where_: LuInt = rbeg;
    while this.windex[where_ as usize] != pivot_col {
        debug_assert!(where_ < rend - 1);
        where_ += 1;
    }
    lu_iswap(this.windex, rbeg as usize, where_ as usize);

    // Check if room is available in W.
    // Columns can be updated in place but the updated row may need to be
    // expanded.
    let om = (m + other_row) as usize;
    let nz = this.wend[om] - this.wbegin[om];
    // BASICLU quirk: plain assignment here (not `+=`), so the double is
    // truncated to lu_int directly.
    let grow = fill_room(nz + rnz1, stretch, pad) as LuInt;
    let mut room = this.wend[m2] - this.wbegin[m2];
    if grow > room {
        lu_file_compress(
            2 * m,
            this.wbegin,
            this.wend,
            this.wflink,
            this.windex,
            this.wvalue,
            stretch,
            pad,
        );
        cbeg = this.wbegin[pc];
        rbeg = this.wbegin[prow_m];
        rend = this.wend[prow_m];
        room = this.wend[m2] - this.wbegin[m2];
        this.ngarbage += 1;
    }
    if grow > room {
        this.addmem_w = grow - room;
        return BASICLU_REALLOCATE;
    }

    // ----------------------------------------------------------------------
    // Column file update
    // ----------------------------------------------------------------------

    let mut uput = this.ubegin[rank as usize];
    let mut put = rbeg + 1;
    let mut ncancelled: LuInt = 0;
    for rpos in (rbeg + 1)..rend {
        let j = this.windex[rpos as usize];
        debug_assert!(j != pivot_col);
        let mut cmx: f64 = 0.0; // column maximum

        // Find position of pivot row entry and possibly other row entry in
        // column j.
        let mut where_pivot: LuInt = -1;
        let mut where_other: LuInt = -1;
        let mut end = this.wend[j as usize];
        for pos in this.wbegin[j as usize]..end {
            if this.windex[pos as usize] == pivot_row {
                where_pivot = pos;
            } else if this.windex[pos as usize] == other_row {
                where_other = pos;
            } else {
                let x = this.wvalue[pos as usize].abs();
                if x > cmx {
                    cmx = x;
                }
            }
        }
        debug_assert!(where_pivot >= 0);
        let xrj = this.wvalue[where_pivot as usize];

        // Store pivot row entry in U.
        if this.wvalue[where_pivot as usize].abs() > droptol {
            this.uindex[uput as usize] = j;
            this.uvalue[uput as usize] = this.wvalue[where_pivot as usize];
            uput += 1;
        }

        let ju = j as usize;
        if where_other == -1 {
            // Compute fill-in element.
            let x = -xrj * (other_value / pivot);
            let xabs = x.abs();
            if xabs > droptol {
                // Store fill-in where pivot row entry was.
                this.windex[where_pivot as usize] = other_row;
                this.wvalue[where_pivot as usize] = x;
                this.windex[put as usize] = j;
                put += 1;
                if xabs > cmx {
                    cmx = xabs;
                }
            } else {
                // Remove pivot row entry.
                this.wend[ju] -= 1;
                end = this.wend[ju];
                this.windex[where_pivot as usize] = this.windex[end as usize];
                this.wvalue[where_pivot as usize] = this.wvalue[end as usize];

                // Decrease column count.
                let nz = end - this.wbegin[ju];
                lu_list_move(
                    j,
                    nz,
                    this.colcount_flink,
                    this.colcount_blink,
                    m,
                    Some(&mut this.min_colnz),
                );
            }
        } else {
            // Remove pivot row entry and update other row entry.
            this.wend[ju] -= 1;
            end = this.wend[ju];
            this.windex[where_pivot as usize] = this.windex[end as usize];
            this.wvalue[where_pivot as usize] = this.wvalue[end as usize];
            if where_other == end {
                where_other = where_pivot;
            }
            let wo = where_other as usize;
            this.wvalue[wo] -= xrj * (other_value / pivot);

            // If we have numerical cancellation, then remove the entry and mark
            // the column.
            let x = this.wvalue[wo].abs();
            if x <= droptol {
                this.wend[ju] -= 1;
                end = this.wend[ju];
                this.windex[wo] = this.windex[end as usize];
                this.wvalue[wo] = this.wvalue[end as usize];
                this.iwork0[ju] = 1;
                ncancelled += 1;
            } else if x > cmx {
                cmx = x;
            }

            // Decrease column count.
            let nz = this.wend[ju] - this.wbegin[ju];
            lu_list_move(
                j,
                nz,
                this.colcount_flink,
                this.colcount_blink,
                m,
                Some(&mut this.min_colnz),
            );
        }
        this.col_pivot[ju] = cmx;
    }
    rend = put;
    this.ubegin[rank as usize + 1] = uput;

    // ----------------------------------------------------------------------
    // Row file update
    // ----------------------------------------------------------------------

    // If we have numerical cancellation, then we have to remove these entries
    // (marked) from the row pattern. In any case remove pivot column entry.
    if ncancelled != 0 {
        debug_assert!(this.iwork0[pc] == 0);
        this.iwork0[pc] = 1; // treat as cancelled
        let mut put = this.wbegin[om]; // compress remaining entries
        let end = this.wend[om];
        let mut pos = put;
        while pos < end {
            let j = this.windex[pos as usize];
            if this.iwork0[j as usize] != 0 {
                this.iwork0[j as usize] = 0;
            } else {
                this.windex[put as usize] = j;
                put += 1;
            }
            pos += 1;
        }
        debug_assert!(end - put == ncancelled + 1);
        this.wend[om] = put;
    } else {
        let mut w = this.wbegin[om];
        while this.windex[w as usize] != pivot_col {
            debug_assert!(w < this.wend[om] - 1);
            w += 1;
        }
        this.wend[om] -= 1;
        let end = this.wend[om];
        this.windex[w as usize] = this.windex[end as usize];
    }

    // Reappend row if no room for update.
    let nfill = rend - (rbeg + 1);
    let room = this.wbegin[this.wflink[om] as usize] - this.wend[om];
    if nfill > room {
        let nz = this.wend[om] - this.wbegin[om];
        let space = (nfill as f64 + stretch * ((nz + nfill) as f64) + (pad as f64)) as LuInt;
        lu_file_reappend(
            m + other_row,
            2 * m,
            this.wbegin,
            this.wend,
            this.wflink,
            this.wblink,
            this.windex,
            this.wvalue,
            space,
        );
        this.nexpand += 1;
    }

    // Append fill-in to row pattern.
    let mut put = this.wend[om];
    for pos in (rbeg + 1)..rend {
        this.windex[put as usize] = this.windex[pos as usize];
        put += 1;
    }
    this.wend[om] = put;

    // Reinsert other row into row counts.
    let nz = this.wend[om] - this.wbegin[om];
    lu_list_move(
        other_row,
        nz,
        this.rowcount_flink,
        this.rowcount_blink,
        m,
        Some(&mut this.min_rownz),
    );

    // Store column in L.
    let mut put = this.lbegin_p[rank as usize];
    let x = other_value / pivot;
    if x.abs() > droptol {
        this.lindex[put as usize] = other_row;
        this.lvalue[put as usize] = x;
        put += 1;
    }
    this.lindex[put as usize] = -1; // terminate column
    put += 1;
    this.lbegin_p[rank as usize + 1] = put;

    // Cleanup (see lu_pivot_any).
    this.col_pivot[pc] = pivot;
    this.wend[pc] = cbeg;
    this.wend[prow_m] = rbeg;
    lu_list_remove(this.colcount_flink, this.colcount_blink, pivot_col);
    lu_list_remove(this.rowcount_flink, this.rowcount_blink, pivot_row);

    BASICLU_OK
}

/// HiGHS: basiclu/lu_pivot.c::lu_remove_col (lines 1215-1247)
fn lu_remove_col(this: &mut Lu<'_>, j: LuInt) {
    let m = this.m;
    let ju = j as usize;
    let cbeg = this.wbegin[ju];
    let cend = this.wend[ju];

    // Remove column j from row file.
    for pos in cbeg..cend {
        let i = this.windex[pos as usize];
        let mi = (m + i) as usize;
        let mut w = this.wbegin[mi];
        while this.windex[w as usize] != j {
            debug_assert!(w < this.wend[mi] - 1);
            w += 1;
        }
        this.wend[mi] -= 1;
        this.windex[w as usize] = this.windex[this.wend[mi] as usize];
        let nz = this.wend[mi] - this.wbegin[mi];
        lu_list_move(
            i,
            nz,
            this.rowcount_flink,
            this.rowcount_blink,
            m,
            Some(&mut this.min_rownz),
        );
    }

    // Remove column j from column file.
    this.col_pivot[ju] = 0.0;
    this.wend[ju] = cbeg;
    lu_list_move(
        j,
        0,
        this.colcount_flink,
        this.colcount_blink,
        m,
        Some(&mut this.min_colnz),
    );
}
