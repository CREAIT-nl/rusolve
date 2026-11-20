//! HiGHS source: ipm/basiclu/lu_solve_sparse.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::lu_solve_symbolic::lu_solve_symbolic;
use super::lu_solve_triangular::lu_solve_triangular;
use super::*;

/// HiGHS: basiclu/lu_solve_sparse.c::lu_solve_sparse (lines 10-284)
///
/// `pattern_symb = iwork1`, `pattern = iwork1 + m`: split `rowcount_blink`
/// at `m`. `pstack = (void *) work1`: pun ([`pun_set_int`]).
pub fn lu_solve_sparse(
    this: &mut Lu<'_>,
    nrhs: LuInt,
    irhs: &[LuInt],
    xrhs: &[f64],
    p_nlhs: &mut LuInt,
    ilhs: &mut [LuInt],
    xlhs: &mut [f64],
    trans: u8,
) {
    let m = this.m;
    let mu = m as usize;
    let nforrest = this.nforrest;
    let pivotlen = this.pivotlen;
    // C: const lu_int nz_sparse = this->sparse_thres * m;  (double -> int)
    let nz_sparse: LuInt = (this.sparse_thres * (m as f64)) as LuInt;
    let droptol = this.droptol;

    // C aliases (see the table on `Lu`)
    let p = &this.wblink[mu + 1..]; // p = Wblink + m + 1
    let pmap = &this.pinv[..];
    let qmap = &this.qinv[..];
    let eta_row = &this.rowcount_flink[mu + 1..];
    let pivotcol = &this.colcount_flink[..];
    let pivotrow = &this.colcount_blink[..];
    let lbegin = &this.wbegin[mu + 1..]; // Lbegin = Wbegin + m + 1
    let ltbegin = &this.wend[mu + 1..]; // Ltbegin = Wend + m + 1
    let ltbegin_p = &this.wflink[mu + 1..]; // Ltbegin_p = Wflink + m + 1
    let ubegin = &this.ubegin[..];
    let rbegin = &this.rowcount_flink[..]; // Rbegin = rowcount_flink
    let wbegin = &this.wbegin[..];
    let wend = &this.wend[..];
    let col_pivot = &this.col_pivot[..];
    let row_pivot = &this.row_pivot[..];
    let lindex = &this.lindex[..];
    let lvalue = &this.lvalue[..];
    let uindex = &this.uindex[..];
    let uvalue = &this.uvalue[..];
    let windex = &this.windex[..];
    let wvalue = &this.wvalue[..];
    let marked = &mut this.iwork0[..]; // marked = iwork0

    // pattern_symb = iwork1, pattern = iwork1 + m
    let (pattern_symb_all, pattern_all) = this.rowcount_blink.split_at_mut(mu);
    let pattern_symb = &mut pattern_symb_all[..mu];
    let pattern = &mut pattern_all[..mu];
    // work = work0; pstack = (void *) work1
    let work = &mut this.work0[..];
    let pstack = &mut this.work1[..];

    let mut lflops: LuInt = 0;
    let mut uflops: LuInt = 0;
    let mut rflops: LuInt = 0;

    // Declared without initialisers: every path assigns before the first read.
    let mut i: LuInt;
    let mut j: LuInt;
    let mut t: LuInt;
    let mut n: LuInt;
    let mut top: LuInt;
    let mut pos: LuInt;
    let mut ipivot: LuInt;
    let mut jpivot: LuInt;
    let mut nz: LuInt;
    let mut nz_symb: LuInt;
    let mut big_m: LuInt; // C: M

    if trans == b't' || trans == b'T' {
        // ----------------------- //
        // Solve transposed system //
        // ----------------------- //

        // Sparse triangular solve with U'.
        // Solution scattered into work, indices in pattern[0..nz-1].
        this.marker += 1;
        big_m = this.marker;
        top = lu_solve_symbolic(
            m,
            wbegin,
            Some(wend),
            windex,
            nrhs,
            irhs,
            pattern_symb,
            pstack,
            marked,
            big_m,
        );
        nz_symb = m - top;

        for n in 0..nrhs {
            work[irhs[n as usize] as usize] = xrhs[n as usize];
        }
        nz = lu_solve_triangular(
            nz_symb,
            &pattern_symb[top as usize..],
            wbegin,
            Some(wend),
            windex,
            wvalue,
            Some(col_pivot),
            droptol,
            work,
            pattern,
            &mut uflops,
        );

        // Permute solution into xlhs.
        // Map pattern from column indices to row indices.
        this.marker += 1;
        big_m = this.marker;
        for n in 0..nz {
            let nu = n as usize;
            j = pattern[nu];
            i = pmap[j as usize];
            pattern[nu] = i;
            xlhs[i as usize] = work[j as usize];
            work[j as usize] = 0.0;
            marked[i as usize] = big_m;
        }

        // Solve with update etas.
        // Append fill-in to pattern.
        for t in (0..nforrest).rev() {
            ipivot = eta_row[t as usize];
            if xlhs[ipivot as usize] != 0.0 {
                let x = xlhs[ipivot as usize];
                for pos in rbegin[t as usize]..rbegin[t as usize + 1] {
                    let pu = pos as usize;
                    i = lindex[pu];
                    if marked[i as usize] != big_m {
                        marked[i as usize] = big_m;
                        pattern[nz as usize] = i;
                        nz += 1;
                    }
                    xlhs[i as usize] -= x * lvalue[pu];
                    rflops += 1;
                }
            }
        }

        if nz <= nz_sparse {
            // Sparse triangular solve with L'.
            // Solution scattered into xlhs, indices in ilhs[0..nz-1].
            this.marker += 1;
            big_m = this.marker;
            top = lu_solve_symbolic(
                m,
                ltbegin,
                None,
                lindex,
                nz,
                &pattern[..nz as usize],
                pattern_symb,
                pstack,
                marked,
                big_m,
            );
            nz_symb = m - top;

            nz = lu_solve_triangular(
                nz_symb,
                &pattern_symb[top as usize..],
                ltbegin,
                None,
                lindex,
                lvalue,
                None,
                droptol,
                xlhs,
                ilhs,
                &mut lflops,
            );
            *p_nlhs = nz;
        } else {
            // Sequential triangular solve with L'.
            // Solution scattered into xlhs, indices in ilhs[0..nz-1].
            nz = 0;
            for k in (0..mu).rev() {
                ipivot = p[k];
                if xlhs[ipivot as usize] != 0.0 {
                    let x = xlhs[ipivot as usize];
                    // for (pos = Ltbegin_p[k]; (i = Lindex[pos]) >= 0; pos++)
                    pos = ltbegin_p[k];
                    loop {
                        i = lindex[pos as usize];
                        if i < 0 {
                            break;
                        }
                        xlhs[i as usize] -= x * lvalue[pos as usize];
                        lflops += 1;
                        pos += 1;
                    }
                    if x.abs() > droptol {
                        ilhs[nz as usize] = ipivot;
                        nz += 1;
                    } else {
                        xlhs[ipivot as usize] = 0.0;
                    }
                }
            }
            *p_nlhs = nz;
        }
    } else {
        // -------------------- //
        // Solve forward system //
        // -------------------- //

        // Sparse triangular solve with L.
        // Solution scattered into work, indices in pattern[0..nz-1].
        this.marker += 1;
        big_m = this.marker;
        top = lu_solve_symbolic(
            m,
            lbegin,
            None,
            lindex,
            nrhs,
            irhs,
            pattern_symb,
            pstack,
            marked,
            big_m,
        );
        nz_symb = m - top;

        for n in 0..nrhs {
            work[irhs[n as usize] as usize] = xrhs[n as usize];
        }
        nz = lu_solve_triangular(
            nz_symb,
            &pattern_symb[top as usize..],
            lbegin,
            None,
            lindex,
            lvalue,
            None,
            droptol,
            work,
            pattern,
            &mut lflops,
        );

        // unmark cancellation
        if nz < nz_symb {
            t = top;
            n = 0;
            while n < nz {
                i = pattern_symb[t as usize];
                if i == pattern[n as usize] {
                    n += 1;
                } else {
                    marked[i as usize] -= 1;
                }
                t += 1;
            }
            while t < m {
                marked[pattern_symb[t as usize] as usize] -= 1;
                t += 1;
            }
        }

        // Solve with update etas.
        // Append fill-in to pattern.
        pos = rbegin[0];
        for t in 0..nforrest {
            ipivot = eta_row[t as usize];
            let mut x: f64 = 0.0;
            while pos < rbegin[t as usize + 1] {
                x += work[lindex[pos as usize] as usize] * lvalue[pos as usize];
                pos += 1;
            }
            work[ipivot as usize] -= x;
            if x != 0.0 && marked[ipivot as usize] != big_m {
                marked[ipivot as usize] = big_m;
                pattern[nz as usize] = ipivot;
                nz += 1;
            }
        }
        rflops += rbegin[nforrest as usize] - rbegin[0];

        if nz <= nz_sparse {
            // Sparse triangular solve with U.
            // Solution scattered into work, indices in ilhs[0..nz-1].
            this.marker += 1;
            big_m = this.marker;
            top = lu_solve_symbolic(
                m,
                ubegin,
                None,
                uindex,
                nz,
                &pattern[..nz as usize],
                pattern_symb,
                pstack,
                marked,
                big_m,
            );
            nz_symb = m - top;

            nz = lu_solve_triangular(
                nz_symb,
                &pattern_symb[top as usize..],
                ubegin,
                None,
                uindex,
                uvalue,
                Some(row_pivot),
                droptol,
                work,
                ilhs,
                &mut uflops,
            );

            // Permute solution into xlhs.
            // Map pattern from row indices to column indices.
            for n in 0..nz {
                let nu = n as usize;
                i = ilhs[nu];
                j = qmap[i as usize];
                ilhs[nu] = j;
                xlhs[j as usize] = work[i as usize];
                work[i as usize] = 0.0;
            }
        } else {
            // Sequential triangular solve with U.
            // Solution computed in work and permuted into xlhs.
            // Pattern (in column indices) stored in ilhs[0..nz-1].
            nz = 0;
            for k in (0..pivotlen).rev() {
                ipivot = pivotrow[k as usize];
                jpivot = pivotcol[k as usize];
                if work[ipivot as usize] != 0.0 {
                    let x = work[ipivot as usize] / row_pivot[ipivot as usize];
                    work[ipivot as usize] = 0.0;
                    // for (pos = Ubegin[ipivot]; (i = Uindex[pos]) >= 0; pos++)
                    pos = ubegin[ipivot as usize];
                    loop {
                        i = uindex[pos as usize];
                        if i < 0 {
                            break;
                        }
                        work[i as usize] -= x * uvalue[pos as usize];
                        uflops += 1;
                        pos += 1;
                    }
                    if x.abs() > droptol {
                        ilhs[nz as usize] = jpivot;
                        nz += 1;
                        xlhs[jpivot as usize] = x;
                    }
                }
            }
        }
        *p_nlhs = nz;
    }

    this.lflops += lflops;
    this.uflops += uflops;
    this.rflops += rflops;
    this.update_cost_numer += rflops as f64;
}
