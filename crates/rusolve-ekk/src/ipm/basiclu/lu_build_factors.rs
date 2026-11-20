//! HiGHS source: ipm/basiclu/lu_build_factors.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_file::lu_file_empty;
use super::lu_internal::{Lu, lu_time};
use super::lu_list::lu_list_move;
use super::*;

/// HiGHS: basiclu/lu_build_factors.c::lu_build_factors (lines 120-441)
///
/// PORT NOTE: the C aliases are carved out of the first-list regions at the
/// C offsets (see the table on [`Lu`]), once, after the memory check:
/// `Wbegin`/`Lbegin` split `wbegin` at `m+1`, `Wend`/`Ltbegin` split `wend`,
/// `Wflink`/`Ltbegin_p` split `wflink`, `Wblink`/`p` split `wblink`. The
/// `_w` halves hold `Wbegin[m]` / `Wend[m]` / head `m` used by `lu_file_empty`
/// and `lu_list_move`. `pmap`/`qmap` are `pinv`/`qinv` (same memory, as in C).
pub fn lu_build_factors(this: &mut Lu<'_>) -> LuInt {
    let m = this.m;
    let mu = m as usize;
    let rank = this.rank;
    let lmem = this.lmem;
    let umem = this.umem;
    let wmem = this.wmem;
    let pad = this.pad;
    let stretch = this.stretch;

    let mut lnz: LuInt;
    let mut unz: LuInt;
    let mut status = BASICLU_OK;

    // So far L is stored columnwise in Lindex, Lvalue and U stored rowwise
    // in Uindex, Uvalue. The factorization has computed rank columns of L
    // and rank rows of U. If rank < m, then the columns which have not been
    // pivotal will be removed from U.
    lnz = this.lbegin_p[rank as usize];
    lnz -= rank; // because each column is terminated by -1
    unz = this.ubegin[rank as usize]; // might be decreased when rank < m

    // Calculate memory and reallocate. The rowwise and columnwise storage of
    // L both need space for Lnz nonzeros + m terminators. The same for the
    // columnwise storage of U except that Uindex[0] = -1 is reserved to
    // accommodate pointers to empty columns. In the rowwise storage of U each
    // row with nz nonzeros is padded by stretch*nz + pad elements.
    let mut need = 2 * (lnz + m);
    if lmem < need {
        this.addmem_l = need - lmem;
        status = BASICLU_REALLOCATE;
    }
    need = unz + m + 1;
    if umem < need {
        this.addmem_u = need - umem;
        status = BASICLU_REALLOCATE;
    }
    // C: need = Unz + stretch*Unz + m*pad;  (double, then truncated to lu_int)
    need = (unz as f64 + stretch * unz as f64 + (m * pad) as f64) as LuInt;
    if wmem < need {
        this.addmem_w = need - wmem;
        status = BASICLU_REALLOCATE;
    }
    if status != BASICLU_OK {
        return status;
    }

    // C aliases. Each split is at the C offset m+1 (see the PORT NOTE).
    let (wbegin_w, lbegin) = this.wbegin.split_at_mut(mu + 1);
    let (wend_w, ltbegin) = this.wend.split_at_mut(mu + 1);
    let (wflink_w, ltbegin_p) = this.wflink.split_at_mut(mu + 1);
    let (wblink_w, p) = this.wblink.split_at_mut(mu + 1);
    let pinv = &mut this.pinv[..]; // pinv, pmap (same memory)
    let qinv = &mut this.qinv[..]; // qinv, qmap (same memory)
    let pivotcol = &mut this.colcount_flink[..];
    let pivotrow = &mut this.colcount_blink[..];
    let lbegin_p = &mut this.lbegin_p[..];
    let ubegin = &mut this.ubegin[..];
    let rbegin = &mut this.rowcount_flink[..]; // Rbegin = rowcount_flink
    let iwork1 = &mut this.rowcount_blink[..];
    let col_pivot = &mut this.col_pivot[..];
    let row_pivot = &mut this.row_pivot[..];
    let lindex = &mut this.lindex[..];
    let lvalue = &mut this.lvalue[..];
    let uindex = &mut this.uindex[..];
    let uvalue = &mut this.uvalue[..];
    let windex = &mut this.windex[..];
    let wvalue = &mut this.wvalue[..];

    // ------------------ //
    // Build permutations //
    // ------------------ //

    // Append columns/rows which have not been pivotal to the end of the
    // pivot sequence. Build pivotrow, pivotcol as inverse of pinv, qinv.
    if cfg!(debug_assertions) {
        for k in 0..mu {
            pivotrow[k] = -1;
        }
        for k in 0..mu {
            pivotcol[k] = -1;
        }
    }

    let mut lrank = rank;
    for i in 0..mu {
        if pinv[i] < 0 {
            pinv[i] = lrank;
            lrank += 1;
        }
        pivotrow[pinv[i] as usize] = i as LuInt;
    }
    debug_assert_eq!(lrank, m);
    lrank = rank;
    for j in 0..mu {
        if qinv[j] < 0 {
            qinv[j] = lrank;
            lrank += 1;
        }
        pivotcol[qinv[j] as usize] = j as LuInt;
    }
    debug_assert_eq!(lrank, m);

    for k in 0..mu {
        debug_assert!(pivotrow[k] >= 0);
    }
    for k in 0..mu {
        debug_assert!(pivotcol[k] >= 0);
    }

    // Dependent columns get unit pivot elements.
    for k in rank..m {
        col_pivot[pivotcol[k as usize] as usize] = 1.0;
    }

    // ----------------------- //
    // Lower triangular factor //
    // ----------------------- //

    // L columnwise. If rank < m, then complete with unit columns (no
    // off-diagonals, so nothing to store here).
    let mut put = lbegin_p[rank as usize];
    for k in rank..m {
        lindex[put as usize] = -1;
        put += 1;
        lbegin_p[k as usize + 1] = put;
    }
    debug_assert!(lbegin_p[mu] == lnz + m);
    for i in 0..mu {
        lbegin[i] = lbegin_p[pinv[i] as usize];
    }

    // L rowwise.
    iwork1[..mu].fill(0); // row counts
    for get in 0..lnz + m {
        let i = lindex[get as usize];
        if i >= 0 {
            iwork1[i as usize] += 1;
        }
    }
    put = lnz + m; // L rowwise starts here
    for k in 0..mu {
        let i = pivotrow[k] as usize;
        ltbegin_p[k] = put;
        ltbegin[i] = put;
        put += iwork1[i];
        lindex[put as usize] = -1; // terminate row
        put += 1;
        iwork1[i] = ltbegin_p[k];
    }
    debug_assert!(put == 2 * (lnz + m));
    for k in 0..mu {
        // fill rows
        let ipivot = pivotrow[k];
        let mut get = lbegin_p[k];
        loop {
            let i = lindex[get as usize];
            if i < 0 {
                break;
            }
            let slot = iwork1[i as usize]; // put into row i
            iwork1[i as usize] += 1;
            lindex[slot as usize] = ipivot;
            lvalue[slot as usize] = lvalue[get as usize];
            get += 1;
        }
    }

    for i in 0..mu {
        debug_assert!(lindex[iwork1[i] as usize] == -1);
    }
    rbegin[0] = 2 * (lnz + m); // beginning of update etas

    // ----------------------- //
    // Upper triangular factor //
    // ----------------------- //

    // U rowwise.
    lu_file_empty(m, wbegin_w, wend_w, wflink_w, wblink_w, wmem);
    iwork1[..mu].fill(0); // column counts
    put = 0;

    // Use separate loops for full rank and rank deficient factorizations. In
    // the first case no elements are removed from U, so skip the test.
    if rank == m {
        for k in 0..mu {
            let jpivot = pivotcol[k];
            wbegin_w[jpivot as usize] = put;
            let mut nz: LuInt = 0;
            for pos in ubegin[k]..ubegin[k + 1] {
                let j = uindex[pos as usize];
                windex[put as usize] = j;
                wvalue[put as usize] = uvalue[pos as usize];
                put += 1;
                iwork1[j as usize] += 1;
                nz += 1;
            }
            wend_w[jpivot as usize] = put;
            // put += stretch*nz + pad;  (double, then truncated to lu_int)
            put = (put as f64 + (stretch * nz as f64 + pad as f64)) as LuInt;
            lu_list_move(jpivot, 0, wflink_w, wblink_w, m, None);
        }
    } else {
        unz = 0; // actual number of nonzeros
        for k in 0..rank as usize {
            let jpivot = pivotcol[k];
            wbegin_w[jpivot as usize] = put;
            let mut nz: LuInt = 0;
            for pos in ubegin[k]..ubegin[k + 1] {
                let j = uindex[pos as usize];
                if qinv[j as usize] < rank {
                    windex[put as usize] = j;
                    wvalue[put as usize] = uvalue[pos as usize];
                    put += 1;
                    iwork1[j as usize] += 1;
                    nz += 1;
                }
            }
            wend_w[jpivot as usize] = put;
            put = (put as f64 + (stretch * nz as f64 + pad as f64)) as LuInt;
            lu_list_move(jpivot, 0, wflink_w, wblink_w, m, None);
            unz += nz;
        }
        for k in rank..m {
            let jpivot = pivotcol[k as usize];
            wbegin_w[jpivot as usize] = put;
            wend_w[jpivot as usize] = put;
            put += pad;
            lu_list_move(jpivot, 0, wflink_w, wblink_w, m, None);
        }
    }
    debug_assert!(put <= wend_w[mu]);
    wbegin_w[mu] = put; // beginning of free space

    // U columnwise.
    uindex[0] = -1;
    put = 1;
    for k in 0..mu {
        // set column pointers
        let j = pivotcol[k] as usize;
        let i = pivotrow[k] as usize;
        let nz = iwork1[j];
        if nz == 0 {
            ubegin[i] = 0; // empty columns all in position 0
        } else {
            ubegin[i] = put;
            put += nz;
            uindex[put as usize] = -1; // terminate column
            put += 1;
        }
        iwork1[j] = ubegin[i];
    }
    ubegin[mu] = put;
    for k in 0..mu {
        // fill columns
        let jpivot = pivotcol[k] as usize;
        let i = pivotrow[k];
        for pos in wbegin_w[jpivot]..wend_w[jpivot] {
            let j = windex[pos as usize] as usize;
            let slot = iwork1[j];
            iwork1[j] += 1;
            debug_assert!(slot >= 1);
            uindex[slot as usize] = i;
            uvalue[slot as usize] = wvalue[pos as usize];
        }
    }

    for j in 0..mu {
        debug_assert!(uindex[iwork1[j] as usize] == -1);
    }

    // -------------------- //
    // Build pivot sequence //
    // -------------------- //

    // Build row-column mappings, overwriting pinv, qinv.
    for k in 0..mu {
        let i = pivotrow[k];
        let j = pivotcol[k];
        pinv[j as usize] = i; // pmap[j] = i
        qinv[i as usize] = j; // qmap[i] = j
    }

    // Build pivots by row index.
    let mut max_pivot: f64 = 0.0;
    let mut min_pivot: f64 = f64::INFINITY;
    for i in 0..mu {
        row_pivot[i] = col_pivot[qinv[i] as usize];
        let pivot = row_pivot[i].abs();
        max_pivot = pivot.max(max_pivot); // fmax
        min_pivot = pivot.min(min_pivot); // fmin
    }

    // memcpy(this->p, pivotrow, m*sizeof(lu_int));
    p[..mu].copy_from_slice(&pivotrow[..mu]);

    this.min_pivot = min_pivot;
    this.max_pivot = max_pivot;
    this.pivotlen = m;
    this.lnz = lnz;
    this.unz = unz;
    this.rnz = 0;
    if lu_time() - this.basiclu_time_start > this.basiclu_time_limit {
        return BASICLU_WARNING_timeout;
    }
    status
}
