//! HiGHS source: ipm/basiclu/lu_singletons.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU source)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::{Lu, lu_time};
use super::*;

/// HiGHS: basiclu/lu_singletons.c::lu_singletons (lines 140-286)
///
/// B rowwise is built in W (`Btp = Wbegin`, `Bti = Windex`, `Btx = Wvalue`).
///
/// PORT NOTE: `iwork1` (row counts, then `iset`) and `iwork2` (`queue`) are the
/// two halves of `rowcount_blink` (`iwork2 = iwork1 + m`), split at `m`.
pub fn lu_singletons(
    this: &mut Lu<'_>,
    bbegin: &[LuInt],
    bend: &[LuInt],
    bi: &[LuInt],
    bx: &[f64],
) -> LuInt {
    let m = this.m;
    let mu = m as usize;
    let lmem = this.lmem;
    let umem = this.umem;
    let wmem = this.wmem;
    let abstol = this.abstol;
    let nzbias = this.nzbias;

    let (iwork1, iwork2) = this.rowcount_blink.split_at_mut(mu);

    // -------------------------------- */
    // Check matrix and build transpose */
    // -------------------------------- */
    // Check pointers and count nnz(B).
    let mut bnz: LuInt = 0;
    let mut ok = true;
    for j in 0..mu {
        if bend[j] < bbegin[j] {
            ok = false;
            break;
        }
        bnz += bend[j] - bbegin[j];
    }
    if !ok {
        return BASICLU_ERROR_invalid_argument;
    }

    // Check if sufficient memory in L, U, W.
    ok = true;
    if lmem < bnz {
        this.addmem_l = bnz - lmem;
        ok = false;
    }
    if umem < bnz {
        this.addmem_u = bnz - umem;
        ok = false;
    }
    if wmem < bnz {
        this.addmem_w = bnz - wmem;
        ok = false;
    }
    if !ok {
        return BASICLU_REALLOCATE;
    }

    // Count nz per row, check indices.
    iwork1[..mu].fill(0); // row counts
    ok = true;
    for j in 0..mu {
        for pos in bbegin[j]..bend[j] {
            let i = bi[pos as usize];
            if i < 0 || i >= m {
                ok = false;
                break;
            }
            iwork1[i as usize] += 1;
        }
        if !ok {
            break;
        }
    }
    if !ok {
        return BASICLU_ERROR_invalid_argument;
    }

    // Pack matrix rowwise, check for duplicates.
    let mut put: LuInt = 0;
    for i in 0..mu {
        // set row pointers
        this.wbegin[i] = put;
        put += iwork1[i];
        iwork1[i] = this.wbegin[i];
    }
    this.wbegin[mu] = put;
    debug_assert_eq!(put, bnz);
    ok = true;
    for j in 0..mu {
        // fill rows
        for pos in bbegin[j]..bend[j] {
            let i = bi[pos as usize] as usize;
            let p = iwork1[i];
            iwork1[i] += 1; // put = iwork1[i]++
            this.windex[p as usize] = j as LuInt;
            this.wvalue[p as usize] = bx[pos as usize];
            if p > this.wbegin[i] && this.windex[p as usize - 1] == j as LuInt {
                ok = false;
            }
        }
    }
    if !ok {
        return BASICLU_ERROR_invalid_argument;
    }

    // ---------------- */
    // Pivot singletons */
    // ---------------- */
    // No pivot rows or pivot columns so far.
    for i in 0..mu {
        this.pinv[i] = -1;
    }
    for j in 0..mu {
        this.qinv[j] = -1;
    }

    // Lbegin_p[0] = Ubegin[0] = rank = 0;
    this.lbegin_p[0] = 0;
    this.ubegin[0] = 0;

    let rank = if nzbias >= 0 {
        // put more in U
        let rank = singleton_cols(
            m,
            bbegin,
            bend,
            bi,
            &this.wbegin[..],
            &this.windex[..],
            &this.wvalue[..],
            &mut this.ubegin[..],
            &mut this.uindex[..],
            &mut this.uvalue[..],
            &mut this.lbegin_p[..],
            &mut this.lindex[..],
            &mut this.col_pivot[..],
            &mut this.pinv[..],
            &mut this.qinv[..],
            iwork1,
            iwork2,
            0,
            abstol,
        );
        singleton_rows(
            m,
            bbegin,
            bend,
            bi,
            bx,
            &this.wbegin[..],
            &this.windex[..],
            &mut this.ubegin[..],
            &mut this.lbegin_p[..],
            &mut this.lindex[..],
            &mut this.lvalue[..],
            &mut this.col_pivot[..],
            &mut this.pinv[..],
            &mut this.qinv[..],
            iwork1,
            iwork2,
            rank,
            abstol,
        )
    } else {
        // put more in L
        let rank = singleton_rows(
            m,
            bbegin,
            bend,
            bi,
            bx,
            &this.wbegin[..],
            &this.windex[..],
            &mut this.ubegin[..],
            &mut this.lbegin_p[..],
            &mut this.lindex[..],
            &mut this.lvalue[..],
            &mut this.col_pivot[..],
            &mut this.pinv[..],
            &mut this.qinv[..],
            iwork1,
            iwork2,
            0,
            abstol,
        );
        singleton_cols(
            m,
            bbegin,
            bend,
            bi,
            &this.wbegin[..],
            &this.windex[..],
            &this.wvalue[..],
            &mut this.ubegin[..],
            &mut this.uindex[..],
            &mut this.uvalue[..],
            &mut this.lbegin_p[..],
            &mut this.lindex[..],
            &mut this.col_pivot[..],
            &mut this.pinv[..],
            &mut this.qinv[..],
            iwork1,
            iwork2,
            rank,
            abstol,
        )
    };

    // pinv, qinv were used as nonzero counters. Reset to -1 if not pivoted.
    for i in 0..mu {
        if this.pinv[i] < 0 {
            this.pinv[i] = -1;
        }
    }
    for j in 0..mu {
        if this.qinv[j] < 0 {
            this.qinv[j] = -1;
        }
    }

    if lu_time() - this.basiclu_time_start > this.basiclu_time_limit {
        return BASICLU_WARNING_timeout;
    }
    this.matrix_nz = bnz;
    this.rank = rank;
    BASICLU_OK
}

/// HiGHS: basiclu/lu_singletons.c::singleton_cols (lines 315-408)
///
/// PORT NOTE: the C `Bx` and `Lx` parameters are unused by this routine and are
/// dropped from the Rust helper (they would be unused-variable warnings).
#[allow(clippy::too_many_arguments)]
fn singleton_cols(
    m: LuInt,
    bbegin: &[LuInt],
    bend: &[LuInt],
    bi: &[LuInt],
    btp: &[LuInt],
    bti: &[LuInt],
    btx: &[f64],
    up: &mut [LuInt],
    ui: &mut [LuInt],
    ux: &mut [f64],
    lp: &mut [LuInt],
    li: &mut [LuInt],
    col_pivot: &mut [f64],
    pinv: &mut [LuInt],
    qinv: &mut [LuInt],
    iset: &mut [LuInt],
    queue: &mut [LuInt],
    rank: LuInt,
    abstol: f64,
) -> LuInt {
    let mu = m as usize;
    let mut rank = rank;
    let rk_start = rank;
    let mut tail: LuInt = 0;

    // Build index sets and initialize queue.
    for j in 0..mu {
        if qinv[j] < 0 {
            let nz = bend[j] - bbegin[j];
            let mut i: LuInt = 0;
            for pos in bbegin[j]..bend[j] {
                i ^= bi[pos as usize]; // put row into set j
            }
            iset[j] = i;
            qinv[j] = -nz - 1; // use as nonzero counter
            if nz == 1 {
                queue[tail as usize] = j as LuInt;
                tail += 1;
            }
        }
    }

    // Eliminate singleton columns.
    let mut put = up[rank as usize];
    let mut front: LuInt = 0;
    while front < tail {
        let j = queue[front as usize] as usize;
        front += 1;
        debug_assert!(qinv[j] == -2 || qinv[j] == -1);
        if qinv[j] == -1 {
            continue; // empty column in active submatrix
        }
        let i = iset[j];
        debug_assert!(i >= 0 && i < m);
        debug_assert!(pinv[i as usize] < 0);
        let end = btp[i as usize + 1];
        let mut pos = btp[i as usize];
        while bti[pos as usize] != j as LuInt {
            // find pivot
            debug_assert!(pos < end - 1);
            pos += 1;
        }
        let piv = btx[pos as usize];
        if piv == 0.0 || piv.abs() < abstol {
            continue; // skip singularity
        }

        // Eliminate pivot.
        qinv[j] = rank;
        pinv[i as usize] = rank;
        for pos2 in btp[i as usize]..end {
            let j2 = bti[pos2 as usize] as usize;
            if qinv[j2] < 0 {
                // test is mandatory because the initial active submatrix may
                // not be the entire matrix (rows eliminated before)
                ui[put as usize] = j2 as LuInt;
                ux[put as usize] = btx[pos2 as usize];
                put += 1;
                iset[j2] ^= i; // remove i from set j2
                qinv[j2] += 1;
                if qinv[j2] == -2 {
                    queue[tail as usize] = j2 as LuInt; // new singleton
                    tail += 1;
                }
            }
        }
        up[rank as usize + 1] = put;
        col_pivot[j] = piv;
        rank += 1;
    }

    // Put empty columns into L.
    let mut pos = lp[rk_start as usize];
    let mut rk = rk_start;
    while rk < rank {
        li[pos as usize] = -1;
        pos += 1;
        lp[rk as usize + 1] = pos;
        rk += 1;
    }
    rank
}

/// HiGHS: basiclu/lu_singletons.c::singleton_rows (lines 419-511)
///
/// PORT NOTE: the C `Btx`, `Ui` and `Ux` parameters are unused by this routine
/// and are dropped from the Rust helper.
#[allow(clippy::too_many_arguments)]
fn singleton_rows(
    m: LuInt,
    bbegin: &[LuInt],
    bend: &[LuInt],
    bi: &[LuInt],
    bx: &[f64],
    btp: &[LuInt],
    bti: &[LuInt],
    up: &mut [LuInt],
    lp: &mut [LuInt],
    li: &mut [LuInt],
    lx: &mut [f64],
    col_pivot: &mut [f64],
    pinv: &mut [LuInt],
    qinv: &mut [LuInt],
    iset: &mut [LuInt],
    queue: &mut [LuInt],
    rank: LuInt,
    abstol: f64,
) -> LuInt {
    let mu = m as usize;
    let mut rank = rank;
    let rk_start = rank;
    let mut tail: LuInt = 0;

    // Build index sets and initialize queue.
    for i in 0..mu {
        if pinv[i] < 0 {
            let nz = btp[i + 1] - btp[i];
            let mut j: LuInt = 0;
            for pos in btp[i]..btp[i + 1] {
                j ^= bti[pos as usize]; // put column into set i
            }
            iset[i] = j;
            pinv[i] = -nz - 1; // use as nonzero counter
            if nz == 1 {
                queue[tail as usize] = i as LuInt;
                tail += 1;
            }
        }
    }

    // Eliminate singleton rows.
    let mut put = lp[rank as usize];
    let mut front: LuInt = 0;
    while front < tail {
        let i = queue[front as usize] as usize;
        front += 1;
        debug_assert!(pinv[i] == -2 || pinv[i] == -1);
        if pinv[i] == -1 {
            continue; // empty column in active submatrix
        }
        let j = iset[i];
        debug_assert!(j >= 0 && j < m);
        debug_assert!(qinv[j as usize] < 0);
        let end = bend[j as usize];
        let mut pos = bbegin[j as usize];
        while bi[pos as usize] != i as LuInt {
            // find pivot
            debug_assert!(pos < end - 1);
            pos += 1;
        }
        let piv = bx[pos as usize];
        if piv == 0.0 || piv.abs() < abstol {
            continue; // skip singularity
        }

        // Eliminate pivot.
        qinv[j as usize] = rank;
        pinv[i] = rank;
        for pos2 in bbegin[j as usize]..end {
            let i2 = bi[pos2 as usize] as usize;
            if pinv[i2] < 0 {
                // test is mandatory because the initial active submatrix may
                // not be the entire matrix (columns eliminated before)
                li[put as usize] = i2 as LuInt;
                lx[put as usize] = bx[pos2 as usize] / piv;
                put += 1;
                iset[i2] ^= j; // remove j from set i2
                pinv[i2] += 1;
                if pinv[i2] == -2 {
                    queue[tail as usize] = i2 as LuInt; // new singleton
                    tail += 1;
                }
            }
        }
        li[put as usize] = -1; // terminate column
        put += 1;
        lp[rank as usize + 1] = put;
        col_pivot[j as usize] = piv;
        rank += 1;
    }

    // Put empty rows into U.
    let mut pos = up[rk_start as usize];
    let mut rk = rk_start;
    while rk < rank {
        up[rk as usize + 1] = pos;
        rk += 1;
    }
    rank
}

#[cfg(test)]
mod tests {
    use super::super::lu_initialize::lu_initialize;
    use super::super::lu_internal::lu_load;
    use super::*;

    /// Hand-computed on B = [[2, 0], [1, 4]] (col0 = {r0: 2, r1: 1}, col1 = {r1: 4}).
    /// BIAS >= 0, so singleton columns go first: col1 (pivot 4, row1) at stage 0,
    /// which leaves col0 as a singleton with pivot 2 (row0) at stage 1. The
    /// off-diagonal 1 of row1 in col0 is the only U entry.
    #[test]
    fn singletons_pivot_order_on_2x2() {
        let m: LuInt = 2;
        let mut istore = vec![0; 1024];
        let mut xstore = vec![0.0; 1024];
        lu_initialize(m, &mut istore, &mut xstore);
        xstore[BASICLU_MEMORYL] = 100.0;
        xstore[BASICLU_MEMORYU] = 100.0;
        xstore[BASICLU_MEMORYW] = 100.0;
        let (mut li, mut lx) = (vec![0; 100], vec![0.0; 100]);
        let (mut ui, mut ux) = (vec![0; 100], vec![0.0; 100]);
        let (mut wi, mut wx) = (vec![0; 100], vec![0.0; 100]);
        let mut lu = lu_load(
            &mut istore,
            &mut xstore,
            &mut li,
            &mut lx,
            &mut ui,
            &mut ux,
            &mut wi,
            &mut wx,
        )
        .expect("valid store");

        let bbegin = [0, 2];
        let bend = [2, 3];
        let bi = [0, 1, 1];
        let bx = [2.0, 1.0, 4.0];
        let status = lu_singletons(&mut lu, &bbegin, &bend, &bi, &bx);

        assert_eq!(status, BASICLU_OK);
        assert_eq!(lu.rank, 2);
        assert_eq!(lu.matrix_nz, 3);
        assert_eq!(&lu.pinv[..], &[1, 0]);
        assert_eq!(&lu.qinv[..], &[1, 0]);
        assert_eq!(&lu.col_pivot[..], &[2.0, 4.0]);
        assert_eq!(&lu.ubegin[..3], &[0, 1, 1]);
        assert_eq!(lu.uindex[0], 0);
        assert_eq!(lu.uvalue[0], 1.0);
        assert_eq!(&lu.lbegin_p[..3], &[0, 1, 2]);
    }
}
