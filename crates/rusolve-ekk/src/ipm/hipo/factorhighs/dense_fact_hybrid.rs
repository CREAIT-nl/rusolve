//! HiGHS source: ipm/hipo/factorhighs/DenseFact.h, ipm/hipo/factorhighs/DenseFactHybrid.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The blocked hybrid dense factorisation of a frontal matrix (`denseFactFH`)
//! and the full-to-hybrid format conversion (`denseFactFP2FH`). `B` is the
//! clique (Schur complement) storage; an empty slice stands for `nullptr`.
//! `parnode` is always false in a serial run.

use super::super::{Int, Int64};
use super::data_collector::DataCollector;
use super::*;
use crate::ipm::hipo::auxiliary::getDiagStart;
use crate::ipm::hipo::factorhighs::call_and_time_blas::*;
use crate::ipm::hipo::factorhighs::dense_fact_kernel::denseFactK;
use crate::ipm::hipo::factorhighs::dgemm_parallel::dgemmParallel;
use crate::ipm::hipo::factorhighs::swaps::{applySwaps, permuteWithSwaps};
use crate::util::highs_utils::{cpp_max, cpp_min};

/// HiGHS: DenseFactHybrid.cpp::denseFactFH (lines 13-235)
pub fn denseFactFH(
    format: u8,
    n: Int,
    k: Int,
    nb: Int,
    A: &mut [f64],
    B: &mut [f64],
    pivot_sign: &[Int],
    thresh: f64,
    regval: &Regul,
    totalreg: &mut [f64],
    swaps: &mut [Int],
    pivot_2x2: &mut [f64],
    parnode: bool,
    data: &DataCollector,
) -> Int {
    // Partial blocked factorisation
    // Matrix A is in format FH
    // Matrix B is in format FP, if format == 'P'
    //                       FH, if format == 'H'
    // BLAS calls: dcopy, dscal, daxpy, dgemm, dtrsm

    // HIPO_CLOCK_CREATE: timing dropped (HIPO_TIMING_LEVEL == 0)

    // check input
    if n < 0 || k < 0 || A.is_empty() || (k < n && B.is_empty()) {
        return kRetInvalidInput;
    }

    // quick return
    if n == 0 {
        return kRetOk;
    }

    // number of blocks of columns
    let n_blocks: Int = (k - 1) / nb + 1;

    // start of diagonal blocks
    let mut diag_start: Vec<Int64> = vec![0; n_blocks as usize];
    getDiagStart(n, k, nb, n_blocks, &mut diag_start, false);

    // buffer for copy of block column
    let mut T: Vec<f64> = vec![0.0; (n * nb) as usize];

    // number of rows/columns in the Schur complement
    let ns: Int = n - k;

    // number of blocks in Schur complement
    let s_blocks: Int = (ns - 1) / nb + 1;

    // buffer for full-format of block of columns of Schur complement
    let mut schur_buf: Vec<f64> = Vec::new();
    if format == b'P' {
        schur_buf.resize((ns * nb) as usize, 0.0);
    }

    // ===========================================================================
    // LOOP OVER BLOCKS
    // ===========================================================================
    for j in 0..n_blocks {
        // j is the index of the block column

        // jb is the number of columns
        let jb: Int = cpp_min(nb, k - nb * j);

        // size of current block could be smaller than nb*nb
        let this_diag_size: Int = jb * jb;

        // diagonal block j is A[d_start..]
        let d_start: usize = diag_start[j as usize] as usize;

        // number of rows left below block j
        let M: Int = n - nb * j - jb;

        // block of columns below diagonal block j is A[r_start..]
        let r_start: usize = (diag_start[j as usize] + this_diag_size as Int64) as usize;

        // ===========================================================================
        // FACTORISE DIAGONAL BLOCK
        // ===========================================================================
        let mut max_in_R: f64 = -1.0;
        if jb == 1 {
            for i in 0..(M * jb) as Int64 {
                max_in_R = cpp_max(max_in_R, A[r_start + i as usize].abs());
            }
        }

        let reg_start: usize = (j * nb) as usize;
        let mut pivot_sign_current: Vec<Int> =
            pivot_sign[reg_start..reg_start + jb as usize].to_vec();
        let sw_start: usize = (j * nb) as usize;
        let p2_start: usize = (j * nb) as usize;

        // D = A[d_start..] (lo), R = A[r_start..] (hi): disjoint halves of A
        let (lo, hi) = A.split_at_mut(r_start);
        let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(10);
        let info = denseFactK(
            b'U',
            jb,
            &mut lo[d_start..],
            jb,
            &mut pivot_sign_current,
            thresh,
            regval,
            &mut totalreg[reg_start..],
            &mut swaps[sw_start..],
            &mut pivot_2x2[p2_start..],
            data,
        );
        drop(_pf);
        if info != 0 {
            return info;
        }

        // #ifdef HIPO_PIVOTING (defined in FactorHiGHSSettings.h)
        // swap columns in R
        applySwaps(&swaps[sw_start..], M, jb, hi, data);

        // unswap regularisation, to keep it with original ordering
        permuteWithSwaps(&mut totalreg[reg_start..], &swaps[sw_start..], jb, true);

        if M > 0 {
            // ===========================================================================
            // SOLVE COLUMNS
            // ===========================================================================
            let d: &[f64] = &lo[d_start..];

            let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(11);
            // solve block R with D
            callAndTime_dtrsm(b'L', b'U', b'T', b'U', jb, M, 1.0, d, jb, hi, jb, data);

            // make copy of partially solved columns
            callAndTime_dcopy(jb * M, &hi[..], 1, &mut T[..], 1, data);

            // solve block R with pivots
            let mut col: Int = 0;
            while col < jb {
                let piv_off: f64 = pivot_2x2[p2_start + col as usize];
                if piv_off == 0.0 {
                    // 1x1 pivots
                    let step: Int = 1;
                    callAndTime_dscal(
                        M,
                        d[(col + jb * col) as usize],
                        &mut hi[col as usize..],
                        jb,
                        data,
                    );
                    col += step;
                } else {
                    // 2x2 pivots
                    let step: Int = 2;

                    // inverse of 2x2 pivot
                    let i_d1: f64 = d[(col + jb * col) as usize];
                    let i_d2: f64 = d[(col + 1 + jb * (col + 1)) as usize];
                    let i_off: f64 = piv_off;

                    // copy of original col1
                    let mut c1_temp: Vec<f64> = vec![0.0; M as usize];
                    callAndTime_dcopy(M, &hi[col as usize..], jb, &mut c1_temp[..], 1, data);

                    // copy of col2 (read before col1 is updated, C++ aliasing of R)
                    let c2_copy: Vec<f64> =
                        (0..M).map(|i| hi[(col + 1 + i * jb) as usize]).collect();

                    // solve col and col+1
                    callAndTime_dscal(M, i_d1, &mut hi[col as usize..], jb, data);
                    callAndTime_daxpy(M, i_off, &c2_copy[..], 1, &mut hi[col as usize..], jb, data);
                    callAndTime_dscal(M, i_d2, &mut hi[(col + 1) as usize..], jb, data);
                    callAndTime_daxpy(
                        M,
                        i_off,
                        &c1_temp[..],
                        1,
                        &mut hi[(col + 1) as usize..],
                        jb,
                        data,
                    );
                    col += step;
                }
            }

            // ===========================================================================
            // UPDATE FRONTAL
            // ===========================================================================
            drop(_pf);
            let mut offset: Int64 = 0;
            let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(12);

            if j + 1 < n_blocks {
                // Q regions (diagonal blocks jj > j) start at or after ds_next;
                // R lies entirely before ds_next.
                let ds_next: usize = diag_start[(j + 1) as usize] as usize;
                let (head, tail) = A.split_at_mut(ds_next);

                // go through remaining blocks of columns
                for jj in (j + 1)..n_blocks {
                    // number of columns in block jj
                    let col_jj: Int = cpp_min(nb, k - nb * jj);

                    // number of rows in block jj
                    let row_jj: Int = n - nb * jj;

                    let P: &[f64] = &T[offset as usize..];
                    let q_off: usize = diag_start[jj as usize] as usize - ds_next;
                    let Q: &mut [f64] = &mut tail[q_off..];
                    let Rjj: &[f64] = &head[r_start + offset as usize..];

                    // perform gemm (potentially) in parallel
                    if parnode {
                        dgemmParallel(P, Rjj, Q, col_jj, jb, row_jj, nb, 1.0, data);
                    } else {
                        callAndTime_dgemm(
                            b'T', b'N', col_jj, row_jj, jb, -1.0, P, jb, Rjj, jb, 1.0, Q, col_jj,
                            data,
                        );
                    }

                    offset += (jb * col_jj) as Int64;
                }
            }

            // ===========================================================================
            // UPDATE SCHUR COMPLEMENT
            // ===========================================================================
            drop(_pf);
            let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(13);
            if k < n {
                let mut B_offset: Int64 = 0;

                // go through blocks of columns of the Schur complement
                for sb in 0..s_blocks {
                    // number of rows of the block
                    let nrow: Int = ns - nb * sb;

                    // number of columns of the block
                    let ncol: Int = cpp_min(nb, nrow);

                    let P: &[f64] = &T[offset as usize..];
                    let Rjj: &[f64] = &A[r_start + offset as usize..];

                    // beta is 0 to avoid initialising schur_buf if format=='P'
                    let beta: f64 = if format == b'P' { 0.0 } else { 1.0 };

                    {
                        let Q: &mut [f64] = if format == b'P' {
                            &mut schur_buf[..]
                        } else {
                            &mut B[B_offset as usize..]
                        };

                        // perform gemm (potentially) in parallel
                        if parnode {
                            dgemmParallel(P, Rjj, Q, ncol, jb, nrow, nb, beta, data);
                        } else {
                            callAndTime_dgemm(
                                b'T', b'N', ncol, nrow, jb, -1.0, P, jb, Rjj, jb, beta, Q, ncol,
                                data,
                            );
                        }
                    }

                    if format == b'P' {
                        // schur_buf contains Schur complement in hybrid format (with full
                        // diagonal blocks). Store it by columns in B (with full diagonal
                        // blocks).
                        for buf_row in 0..nrow {
                            let N: Int = ncol;
                            callAndTime_daxpy(
                                N,
                                1.0,
                                &schur_buf[(buf_row * ncol) as usize..],
                                1,
                                &mut B[(B_offset + buf_row as Int64) as usize..],
                                nrow,
                                data,
                            );
                        }
                    }

                    B_offset += (nrow * ncol) as Int64;
                    offset += (jb * ncol) as Int64;
                }
            }
        }
    }

    kRetOk
}

/// HiGHS: DenseFactHybrid.cpp::denseFactFP2FH (lines 237-276)
pub fn denseFactFP2FH(A: &mut [f64], nrow: Int, ncol: Int, nb: Int, data: &DataCollector) -> Int {
    // Packed to Hybrid conversion
    // Matrix A on  input is in format FP
    // Matrix A on output is in format FH
    // BLAS calls: dcopy

    // HIPO_CLOCK_CREATE: timing dropped (HIPO_TIMING_LEVEL == 0)

    let mut buf: Vec<f64> = vec![0.0; (nrow * nb) as usize];

    let mut startAtoBuf: Int64 = 0;
    let mut startBuftoA: Int64 = 0;

    for k in 0..=((ncol - 1) / nb) {
        // Number of columns in the block. Can be smaller than nb for last block.
        let block_size: Int = cpp_min(nb, ncol - k * nb);

        // Number of rows in the block
        let row_size: Int = nrow - k * nb;

        // Copy block into buf
        callAndTime_dcopy(
            row_size * block_size,
            &A[startAtoBuf as usize..],
            1,
            &mut buf[..],
            1,
            data,
        );
        startAtoBuf += (row_size * block_size) as Int64;

        // Copy columns back into A, row by row.
        // One call of dcopy_ for each row of the block of columns.
        for i in 0..row_size {
            let N: Int = block_size;
            callAndTime_dcopy(
                N,
                &buf[i as usize..],
                row_size,
                &mut A[startBuftoA as usize..],
                1,
                data,
            );
            startBuftoA += N as Int64;
        }
    }

    // HIPO_CLOCK_STOP: timing dropped (HIPO_TIMING_LEVEL == 0)

    kRetOk
}

#[cfg(test)]
mod tests {
    //! Unit tests. Hand-derived from DenseFactHybrid.cpp;
    //! all values are exact in binary floating point.
    use super::*;

    fn fh(
        format: u8,
        n: Int,
        k: Int,
        nb: Int,
        a: &mut [f64],
        b: &mut [f64],
        sign: &[Int],
    ) -> (Int, Vec<f64>, Vec<Int>, Vec<f64>) {
        let len = k.max(0) as usize;
        let mut reg = vec![0.0; len];
        let mut swaps = vec![0 as Int; len];
        let mut p2 = vec![0.0; len];
        let d = DataCollector::new();
        let r = denseFactFH(
            format,
            n,
            k,
            nb,
            a,
            b,
            sign,
            1e-10,
            &Regul::default(),
            &mut reg,
            &mut swaps,
            &mut p2,
            false,
            &d,
        );
        (r, reg, swaps, p2)
    }

    #[test]
    fn invalid_input_and_quick_return() {
        // line 27: n < 0, k < 0, !A, (k < n && !B); line 30: n == 0.
        let mut a = [1.0];
        let mut empty: [f64; 0] = [];
        assert_eq!(fh(b'H', -1, 0, 1, &mut a, &mut [], &[]).0, kRetInvalidInput);
        assert_eq!(fh(b'H', 1, -1, 1, &mut a, &mut [], &[]).0, kRetInvalidInput);
        assert_eq!(
            fh(b'H', 1, 1, 1, &mut empty, &mut [], &[1]).0,
            kRetInvalidInput
        );
        assert_eq!(fh(b'H', 2, 1, 1, &mut a, &mut [], &[1]).0, kRetInvalidInput);
        assert_eq!(fh(b'H', 0, 0, 1, &mut a, &mut [], &[]).0, kRetOk);
        assert_eq!(a, [1.0]);
    }

    #[test]
    fn two_blocks_frontal_update() {
        // [4 2; 2 3], nb = 1: two block columns, FH storage [a00 a10 | a11].
        // Block 0: D = 1/4, R = 2 * 1/4 = 0.5, Q = 3 - 2 * 0.5 = 2.
        // Block 1: D = 1/2.
        let mut a = [4.0, 2.0, 3.0];
        let (r, reg, swaps, p2) = fh(b'H', 2, 2, 1, &mut a, &mut [], &[1, 1]);
        assert_eq!(r, kRetOk);
        assert_eq!(a, [0.25, 0.5, 0.5]);
        assert_eq!(reg, [0.0, 0.0]);
        assert_eq!(swaps, [0, 0]);
        assert_eq!(p2, [0.0, 0.0]);
    }

    #[test]
    fn schur_hybrid_format() {
        // n = 2, k = 1: A = [4 2], clique B = [3] in FH (beta = 1, line 192).
        let mut a = [4.0, 2.0];
        let mut b = [3.0];
        let (r, ..) = fh(b'H', 2, 1, 1, &mut a, &mut b, &[1]);
        assert_eq!(r, kRetOk);
        assert_eq!(a, [0.25, 0.5]);
        assert_eq!(b, [2.0]);
    }

    #[test]
    fn schur_packed_format_accumulates_into_b() {
        // format 'P': gemm into schur_buf with beta = 0, then daxpy adds the
        // buffer into B (lines 204-208), so B keeps its previous contents:
        // 3 + (-2 * 0.5) = 2.
        let mut a = [4.0, 2.0];
        let mut b = [3.0];
        let (r, ..) = fh(b'P', 2, 1, 1, &mut a, &mut b, &[1]);
        assert_eq!(r, kRetOk);
        assert_eq!(b, [2.0]);
    }

    #[test]
    fn swap_applied_to_r() {
        // jb = 2, D = diag(1, 4), R row = [2 8], Schur entry 30.
        // Bunch-Kaufman picks the larger diagonal first: swaps = [1, 1], and
        // applySwaps (line 106) turns R into [8 2]. Then R = [8/4, 2/1] = [2 2]
        // and B = 30 - (8 * 2 + 2 * 2) = 10.
        let mut a = [1.0, 0.0, 0.0, 4.0, 2.0, 8.0];
        let mut b = [30.0];
        let (r, reg, swaps, p2) = fh(b'H', 3, 2, 2, &mut a, &mut b, &[1, 1]);
        assert_eq!(r, kRetOk);
        assert_eq!(swaps, [1, 1]);
        assert_eq!(&a[4..6], &[2.0, 2.0]);
        assert_eq!(a[0], 0.25);
        assert_eq!(a[3], 1.0);
        assert_eq!(reg, [0.0, 0.0]);
        assert_eq!(p2, [0.0, 0.0]);
        assert_eq!(b, [10.0]);
    }

    #[test]
    fn two_by_two_pivot_columns() {
        // D = [0 1; 1 0] forces a 2x2 pivot with inverse [0 1; 1 0]
        // (i_off = 1). R row = [2 3] becomes [3 2] (lines 136-148) and
        // B = 20 - (2 * 3 + 3 * 2) = 8.
        let mut a = [0.0, 1.0, 1.0, 0.0, 2.0, 3.0];
        let mut b = [20.0];
        let (r, _reg, swaps, p2) = fh(b'H', 3, 2, 2, &mut a, &mut b, &[1, -1]);
        assert_eq!(r, kRetOk);
        assert_eq!(swaps, [0, 1]);
        assert_eq!(p2[0], 1.0);
        assert_eq!(&a[4..6], &[3.0, 2.0]);
        assert_eq!(b, [8.0]);
    }

    #[test]
    fn fp_to_fh_conversion() {
        // nrow = 3, ncol = 3, nb = 2. Block 0 is 3x2 column-major
        // [0 1 2 | 3 4 5] and is written back row by row: [0 3 1 4 2 5].
        // Block 1 is the single entry 6.
        let mut a = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let d = DataCollector::new();
        assert_eq!(denseFactFP2FH(&mut a, 3, 3, 2, &d), kRetOk);
        assert_eq!(a, [0.0, 3.0, 1.0, 4.0, 2.0, 5.0, 6.0]);
    }
}
