//! HiGHS source: ipm/basiclu/lu_solve_for_update.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::lu_solve_symbolic::lu_solve_symbolic;
use super::lu_solve_triangular::lu_solve_triangular;
use super::*;

/// HiGHS: basiclu/lu_solve_for_update.c::lu_solve_for_update (lines 10-360)
///
/// `pattern_symb = iwork1`, `pattern = iwork1 + m`: split `rowcount_blink`
/// at `m`. `pstack = (void *) work1`: passed as `work1` directly, since
/// [`lu_solve_symbolic`] takes the `f64` workspace.
///
/// PORT NOTE: `goto done` is a labeled block `'done` with `break 'done`. The
/// C `want_solution` test (`p_nlhs && ilhs && xlhs`) is `lhs.is_some()`; the
/// lhs outputs are only touched after the `want_solution` test, as in C.
///
/// PORT NOTE: C dereferences `xrhs` in the forward branch unconditionally
/// (it is only read when `nrhs > 0`); a `None` with `nrhs > 0` panics here,
/// where C would read through NULL.
///
/// PORT NOTE: the C second-list aliases (`Rbegin`/`eta_row`, `Lbegin`,
/// `Ltbegin`, `Ltbegin_p`, `p`, `Wbegin`/`Wend`) are carved out of the same
/// `rowcount_flink`/`wbegin`/`wend`/`wflink`/`wblink` arrays with
/// `split_at_mut` at offset `m+1`.
#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::needless_range_loop
)]
pub fn lu_solve_for_update(
    this: &mut Lu<'_>,
    nrhs: LuInt,
    irhs: &[LuInt],
    xrhs: Option<&[f64]>,
    p_nlhs: Option<&mut LuInt>,
    ilhs: Option<&mut [LuInt]>,
    xlhs: Option<&mut [f64]>,
    trans: u8,
) -> LuInt {
    let m = this.m;
    let mu = m as usize;
    let nforrest = this.nforrest;
    let pivotlen = this.pivotlen;
    let nz_sparse = (this.sparse_thres * m as f64) as LuInt;
    let droptol = this.droptol;

    // C: want_solution = p_nlhs && ilhs && xlhs
    let mut lhs = match (p_nlhs, ilhs, xlhs) {
        (Some(p), Some(il), Some(xl)) => Some((p, il, xl)),
        _ => None,
    };
    let want_solution = lhs.is_some();
    let mut lflops: LuInt = 0;
    let mut uflops: LuInt = 0;
    let mut rflops: LuInt = 0;

    'done: {
        // Carve the second-list aliases out of the first lists (see lu_internal.rs).
        let (rbegin, eta_row) = this.rowcount_flink.split_at_mut(mu + 1);
        let (wbegin_head, lbegin) = this.wbegin.split_at_mut(mu + 1);
        let (wend_head, ltbegin) = this.wend.split_at_mut(mu + 1);
        let (_, ltbegin_p) = this.wflink.split_at_mut(mu + 1);
        let (_, p) = this.wblink.split_at_mut(mu + 1);
        let (pattern_symb, pattern) = this.rowcount_blink.split_at_mut(mu);

        if trans == b't' || trans == b'T' {
            /* ----------------------- */
            /* Solve transposed system */
            /* ----------------------- */

            let jpivot = irhs[0];
            let ipivot = this.pinv[jpivot as usize];
            let jbegin = wbegin_head[jpivot as usize];
            let jend = wend_head[jpivot as usize];

            /*
             * Compute row eta vector.
             * Symbolic pattern in pattern_symb[top..m-1], indices of (actual)
             * nonzeros in pattern[0..nz-1], values scattered into work.
             * We do not drop small elements to zero, but the symbolic and the
             * numeric pattern will still be different when we have exact
             * cancellation.
             */
            this.marker += 1;
            let mm = this.marker;
            let top = lu_solve_symbolic(
                m,
                &*wbegin_head,
                Some(&*wend_head),
                &*this.windex,
                jend - jbegin,
                &this.windex[jbegin as usize..jend as usize],
                pattern_symb,
                this.work1,
                this.iwork0,
                mm,
            );
            let nz_symb = m - top;

            /* reallocate if not enough memory in Li, Lx (where we store R) */
            let room = this.lmem - rbegin[nforrest as usize];
            if room < nz_symb {
                this.addmem_l = nz_symb - room;
                return BASICLU_REALLOCATE;
            }

            for pos in jbegin..jend {
                this.work0[this.windex[pos as usize] as usize] = this.wvalue[pos as usize];
            }
            let _ = lu_solve_triangular(
                nz_symb,
                &pattern_symb[top as usize..],
                &*wbegin_head,
                Some(&*wend_head),
                &*this.windex,
                &*this.wvalue,
                Some(&*this.col_pivot),
                0.0,
                &mut *this.work0,
                &mut *pattern,
                &mut uflops,
            );

            /*
             * Compress row eta into L, pattern mapped from column to row indices.
             * The triangularity test in lu_update requires the symbolic pattern.
             */
            let mut put = rbegin[nforrest as usize];
            for t in top..m {
                let j = pattern_symb[t as usize];
                let i = this.pinv[j as usize];
                this.lindex[put as usize] = i;
                this.lvalue[put as usize] = this.work0[j as usize];
                put += 1;
                this.work0[j as usize] = 0.0;
            }
            rbegin[(nforrest + 1) as usize] = put;
            eta_row[nforrest as usize] = ipivot;
            this.btran_for_update = jpivot;

            // C: if (!want_solution) goto done;
            let Some((p_nlhs, ilhs, xlhs)) = lhs.take() else {
                break 'done;
            };

            /*
             * Scatter the row eta into xlhs and scale it to become the solution
             * to U^{-1}*[unit vector]. Now we can drop small entries to zero and
             * recompute the numerical pattern.
             */
            this.marker += 1;
            let mm = this.marker;
            pattern[0] = ipivot;
            this.iwork0[ipivot as usize] = mm;
            let pivot = this.col_pivot[jpivot as usize];
            xlhs[ipivot as usize] = 1.0 / pivot;

            let xdrop = droptol * pivot.abs();
            let mut nz: LuInt = 1;
            for pos in rbegin[nforrest as usize]..rbegin[(nforrest + 1) as usize] {
                if this.lvalue[pos as usize].abs() > xdrop {
                    let i = this.lindex[pos as usize];
                    pattern[nz as usize] = i;
                    nz += 1;
                    this.iwork0[i as usize] = mm;
                    xlhs[i as usize] = -this.lvalue[pos as usize] / pivot;
                }
            }

            /*
             * Solve with update etas.
             * Append fill-in to pattern.
             */
            for t in (0..nforrest).rev() {
                let ip = eta_row[t as usize];
                if xlhs[ip as usize] != 0.0 {
                    let x = xlhs[ip as usize];
                    for pos in rbegin[t as usize]..rbegin[(t + 1) as usize] {
                        let i = this.lindex[pos as usize];
                        if this.iwork0[i as usize] != mm {
                            this.iwork0[i as usize] = mm;
                            pattern[nz as usize] = i;
                            nz += 1;
                        }
                        xlhs[i as usize] -= x * this.lvalue[pos as usize];
                        rflops += 1;
                    }
                }
            }

            if nz <= nz_sparse {
                /*
                 * Sparse triangular solve with L'.
                 * Solution scattered into xlhs, indices in ilhs[0..nz-1].
                 */
                this.marker += 1;
                let mm = this.marker;
                let top = lu_solve_symbolic(
                    m,
                    &*ltbegin,
                    None,
                    &*this.lindex,
                    nz,
                    &pattern[..nz as usize],
                    pattern_symb,
                    this.work1,
                    this.iwork0,
                    mm,
                );
                let nz_symb = m - top;
                nz = lu_solve_triangular(
                    nz_symb,
                    &pattern_symb[top as usize..],
                    &*ltbegin,
                    None,
                    &*this.lindex,
                    &*this.lvalue,
                    None,
                    droptol,
                    xlhs,
                    ilhs,
                    &mut lflops,
                );
                *p_nlhs = nz;
            } else {
                /*
                 * Sequential triangular solve with L'.
                 * Solution scattered into xlhs, indices in ilhs[0..nz-1].
                 */
                nz = 0;
                for k in (0..m).rev() {
                    let ip = p[k as usize];
                    if xlhs[ip as usize] != 0.0 {
                        let x = xlhs[ip as usize];
                        let mut pos = ltbegin_p[k as usize];
                        loop {
                            let i = this.lindex[pos as usize];
                            if i < 0 {
                                break;
                            }
                            xlhs[i as usize] -= x * this.lvalue[pos as usize];
                            lflops += 1;
                            pos += 1;
                        }
                        if x.abs() > droptol {
                            ilhs[nz as usize] = ip;
                            nz += 1;
                        } else {
                            xlhs[ip as usize] = 0.0;
                        }
                    }
                }
                *p_nlhs = nz;
            }
        } else {
            /* -------------------- */
            /* Solve forward system */
            /* -------------------- */

            /*
             * Sparse triangular solve with L.
             * Solution scattered into work, indices in pattern[0..nz-1].
             */
            this.marker += 1;
            let mm = this.marker;
            let top = lu_solve_symbolic(
                m,
                &*lbegin,
                None,
                &*this.lindex,
                nrhs,
                irhs,
                pattern_symb,
                this.work1,
                this.iwork0,
                mm,
            );
            let nz_symb = m - top;

            let xr = xrhs.unwrap_or(&[]);
            for n in 0..nrhs as usize {
                this.work0[irhs[n] as usize] = xr[n];
            }
            let mut nz = lu_solve_triangular(
                nz_symb,
                &pattern_symb[top as usize..],
                &*lbegin,
                None,
                &*this.lindex,
                &*this.lvalue,
                None,
                droptol,
                &mut *this.work0,
                &mut *pattern,
                &mut lflops,
            );

            /* unmark cancellation */
            if nz < nz_symb {
                let mut t = top;
                let mut n: LuInt = 0;
                while n < nz {
                    let i = pattern_symb[t as usize];
                    if i == pattern[n as usize] {
                        n += 1;
                    } else {
                        this.iwork0[i as usize] -= 1;
                    }
                    t += 1;
                }
                while t < m {
                    this.iwork0[pattern_symb[t as usize] as usize] -= 1;
                    t += 1;
                }
            }

            /*
             * Solve with update etas.
             * Append fill-in to pattern.
             */
            let mut pos = rbegin[0];
            for t in 0..nforrest {
                let ip = eta_row[t as usize];
                let mut x = 0.0_f64;
                while pos < rbegin[(t + 1) as usize] {
                    x += this.work0[this.lindex[pos as usize] as usize] * this.lvalue[pos as usize];
                    pos += 1;
                }
                this.work0[ip as usize] -= x;
                if x != 0.0 && this.iwork0[ip as usize] != mm {
                    this.iwork0[ip as usize] = mm;
                    pattern[nz as usize] = ip;
                    nz += 1;
                }
            }
            rflops += rbegin[nforrest as usize] - rbegin[0];

            /* reallocate if not enough memory in U */
            let room = this.umem - this.ubegin[mu];
            let need = nz + 1;
            if room < need {
                for n in 0..nz as usize {
                    this.work0[pattern[n] as usize] = 0.0;
                }
                this.addmem_u = need - room;
                return BASICLU_REALLOCATE;
            }

            /* Compress spike into U. */
            let mut put = this.ubegin[mu];
            for n in 0..nz as usize {
                let i = pattern[n];
                this.uindex[put as usize] = i;
                this.uvalue[put as usize] = this.work0[i as usize];
                put += 1;
                if !want_solution {
                    this.work0[i as usize] = 0.0;
                }
            }
            this.uindex[put as usize] = -1; /* terminate column */
            this.ftran_for_update = 0;

            // C: if (!want_solution) goto done;
            let Some((p_nlhs, ilhs, xlhs)) = lhs.take() else {
                break 'done;
            };

            if nz <= nz_sparse {
                /*
                 * Sparse triangular solve with U.
                 * Solution scattered into work, indices in ilhs[0..nz-1].
                 */
                this.marker += 1;
                let mm = this.marker;
                let top = lu_solve_symbolic(
                    m,
                    &*this.ubegin,
                    None,
                    &*this.uindex,
                    nz,
                    &pattern[..nz as usize],
                    pattern_symb,
                    this.work1,
                    this.iwork0,
                    mm,
                );
                let nz_symb = m - top;

                nz = lu_solve_triangular(
                    nz_symb,
                    &pattern_symb[top as usize..],
                    &*this.ubegin,
                    None,
                    &*this.uindex,
                    &*this.uvalue,
                    Some(&*this.row_pivot),
                    droptol,
                    &mut *this.work0,
                    &mut *ilhs,
                    &mut uflops,
                );

                /*
                 * Permute solution into xlhs.
                 * Map pattern from row indices to column indices.
                 */
                for n in 0..nz as usize {
                    let i = ilhs[n];
                    let j = this.qinv[i as usize];
                    ilhs[n] = j;
                    xlhs[j as usize] = this.work0[i as usize];
                    this.work0[i as usize] = 0.0;
                }
            } else {
                /*
                 * Sequential triangular solve with U.
                 * Solution computed in work and permuted into xlhs.
                 * Pattern (in column indices) stored in ilhs[0..nz-1].
                 */
                nz = 0;
                for k in (0..pivotlen).rev() {
                    let ip = this.colcount_blink[k as usize]; // pivotrow
                    let jp = this.colcount_flink[k as usize]; // pivotcol
                    if this.work0[ip as usize] != 0.0 {
                        let x = this.work0[ip as usize] / this.row_pivot[ip as usize];
                        this.work0[ip as usize] = 0.0;
                        let mut pos = this.ubegin[ip as usize];
                        loop {
                            let i = this.uindex[pos as usize];
                            if i < 0 {
                                break;
                            }
                            this.work0[i as usize] -= x * this.uvalue[pos as usize];
                            uflops += 1;
                            pos += 1;
                        }
                        if x.abs() > droptol {
                            ilhs[nz as usize] = jp;
                            nz += 1;
                            xlhs[jp as usize] = x;
                        }
                    }
                }
            }
            *p_nlhs = nz;
        }
    }

    // done:
    this.lflops += lflops;
    this.uflops += uflops;
    this.rflops += rflops;
    this.update_cost_numer += rflops as f64;
    BASICLU_OK
}
