//! HiGHS source: ipm/hipo/factorhighs/DenseFact.h, ipm/hipo/factorhighs/DenseFactKernel.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The dense LDL^T kernel with block Bunch-Kaufman pivoting and static /
//! dynamic regularisation (`denseFactK`).

use super::super::{Int, Int64};
use super::call_and_time_blas::{
    callAndTime_daxpy, callAndTime_dcopy, callAndTime_dger, callAndTime_dscal,
};
use super::data_collector::DataCollector;
use super::swaps::swapCols;
use super::*;
use crate::util::highs_utils::cpp_max;

/// HiGHS: DenseFactKernel.cpp::maxInCol (lines 14-39, file static)
pub fn maxInCol(j: Int, n: Int, m: Int, A: &[f64], lda: Int) -> (Int, f64) {
    // Given the symemtric matrix A, of size nxn, accessed with leading dimension
    // lda, in upper triangular format, ignoring rows 0:j-1, find the maximum in
    // row/col m.

    let mut maxval: f64 = -1.0;
    let mut r: Int = -1;
    // scan column m, rows j:m-1
    let mut row = j;
    while row < m {
        let val = A[row as usize + lda as usize * m as usize].abs();
        if val > maxval {
            maxval = val;
            r = row;
        }
        row += 1;
    }
    // scan row m, columns m+1:n
    let mut col = m + 1;
    while col < n {
        let val = A[m as usize + lda as usize * col as usize].abs();
        if val > maxval {
            maxval = val;
            r = col;
        }
        col += 1;
    }
    (r, maxval)
}

/// HiGHS: DenseFactKernel.cpp::staticReg (lines 41-51, file static)
pub fn staticReg(pivot: &mut f64, sign: Int, regval: &Regul, totalreg: &mut f64) {
    // apply static regularisation

    let old_pivot = *pivot;
    if sign > 0 {
        *pivot += regval.primal;
    } else {
        *pivot -= regval.dual;
    }
    *totalreg = *pivot - old_pivot;
}

/// HiGHS: DenseFactKernel.cpp::blockBunchKaufman (lines 53-147, file static)
pub fn blockBunchKaufman(
    j: Int,
    n: Int,
    A: &mut [f64],
    lda: Int,
    swaps: &mut [Int],
    sign: &mut [Int],
    thresh: f64,
    regval: &Regul,
    totalreg: &mut [f64],
    data: &DataCollector,
) -> bool {
    // Perform Bunch-Kaufman pivoting within a block of the supernode (see Schenk,
    // Gartner, ETNA 2006).
    // It works only for upper triangular A.
    // Return true if 2x2 pivot should be used, false otherwise.
    // Swap of columns may be performed.
    // Regularisation of pivot may be performed.

    // HIPO_CLOCK_CREATE dropped (line 63).

    let mut flag_2x2 = false;

    // Find largest diagonal entry in the residual part of the block
    let mut ind_max_diag: Int = -1;
    let mut max_diag: f64 = -1.0;
    let mut i = j;
    while i < n {
        let val = A[i as usize + lda as usize * i as usize].abs();
        if val > max_diag {
            max_diag = val;
            ind_max_diag = i;
        }
        i += 1;
    }
    debug_assert!(ind_max_diag >= j);

    // put column with max pivot as first column
    swapCols(b'U', n, A, lda, j, ind_max_diag, swaps, sign, data);

    // Max in column j of diagonal block
    let res = maxInCol(j, n, j, A, lda);
    let gamma_j = res.1;
    let r = res.0;
    let Ajj = if sign[j as usize] > 0 {
        A[j as usize + lda as usize * j as usize] + regval.dual
    } else {
        A[j as usize + lda as usize * j as usize] - regval.primal
    };

    if cpp_max(Ajj.abs(), gamma_j) <= thresh || (sign[j as usize] as f64) * Ajj < 0.0 || j == n - 1
    {
        // Must accept current pivot
        let old_pivot = A[j as usize + lda as usize * j as usize];
        staticReg(
            &mut A[j as usize + lda as usize * j as usize],
            sign[j as usize],
            regval,
            &mut totalreg[j as usize],
        );

        if (sign[j as usize] as f64) * A[j as usize + lda as usize * j as usize] < 0.0 {
            data.setWrongSign(A[j as usize + lda as usize * j as usize]);
            // A[j + lda * j] = sign[j] * thresh; (commented out, line 97)
        }

        if cpp_max(Ajj.abs(), gamma_j) < thresh {
            // perturbe pivot
            A[j as usize + lda as usize * j as usize] = sign[j as usize] as f64 * thresh;
            data.countRegPiv();
        }
        totalreg[j as usize] = A[j as usize + lda as usize * j as usize] - old_pivot;
        data.setMaxReg(totalreg[j as usize].abs());
    } else {
        // Max in column r of diagonal block
        debug_assert!(r >= 0);
        let res = maxInCol(j, n, r, A, lda);
        let gamma_r = res.1;
        let Arr = if sign[r as usize] > 0 {
            A[r as usize + lda as usize * r as usize] + regval.primal
        } else {
            A[r as usize + lda as usize * r as usize] - regval.dual
        };

        if Ajj.abs() >= kAlphaBK * gamma_j || Ajj.abs() * gamma_r >= kAlphaBK * gamma_j * gamma_j {
            // Accept current pivot
            staticReg(
                &mut A[j as usize + lda as usize * j as usize],
                sign[j as usize],
                regval,
                &mut totalreg[j as usize],
            );

            if (sign[j as usize] as f64) * A[j as usize + lda as usize * j as usize] < 0.0 {
                data.setWrongSign(A[j as usize + lda as usize * j as usize]);
            }
        } else if Arr.abs() >= kAlphaBK * gamma_r {
            // Use pivot r

            swapCols(b'U', n, A, lda, j, r, swaps, sign, data);
            staticReg(
                &mut A[j as usize + lda as usize * j as usize],
                sign[j as usize],
                regval,
                &mut totalreg[j as usize],
            );

            if (sign[j as usize] as f64) * A[j as usize + lda as usize * j as usize] < 0.0 {
                data.setWrongSign(A[j as usize + lda as usize * j as usize]);
            }
        } else {
            // Use 2x2 pivot (j,r)

            swapCols(b'U', n, A, lda, j + 1, r, swaps, sign, data);
            flag_2x2 = true;
            staticReg(
                &mut A[j as usize + lda as usize * j as usize],
                sign[j as usize],
                regval,
                &mut totalreg[j as usize],
            );
            staticReg(
                &mut A[(j + 1) as usize + lda as usize * (j + 1) as usize],
                sign[(j + 1) as usize],
                regval,
                &mut totalreg[(j + 1) as usize],
            );

            if (sign[j as usize] as f64) * A[j as usize + lda as usize * j as usize] < 0.0 {
                data.setWrongSign(A[j as usize + lda as usize * j as usize]);
            }
            if (sign[(j + 1) as usize] as f64)
                * A[(j + 1) as usize + lda as usize * (j + 1) as usize]
                < 0.0
            {
                data.setWrongSign(A[(j + 1) as usize + lda as usize * (j + 1) as usize]);
            }
        }
    }

    // HIPO_CLOCK_STOP(2, data, kTimeDenseFact_pivoting) dropped (line 145).
    flag_2x2
}

/// Rust-only helper (no HiGHS counterpart): copies `m` entries of `a`
/// starting at `off` with stride `inc` into `out[..m]`. Pure loads/stores, no
/// arithmetic; used only to hand dger an unaliased copy of its `y` row.
fn gatherStrided(a: &[f64], off: usize, inc: Int, m: Int, out: &mut [f64]) {
    for k in 0..m as usize {
        out[k] = a[off + k * inc as usize];
    }
}

/// HiGHS: DenseFactKernel.cpp::denseFactK (lines 149-281)
pub fn denseFactK(
    uplo: u8,
    n: Int,
    A: &mut [f64],
    lda: Int,
    pivot_sign: &mut [Int],
    thresh: f64,
    regval: &Regul,
    totalreg: &mut [f64],
    swaps: &mut [Int],
    pivot_2x2: &mut [f64],
    data: &DataCollector,
) -> Int {
    // ===========================================================================
    // Factorisation kernel
    // Matrix A is in format F
    // BLAS calls: dscal, dcopy, dger
    // ===========================================================================

    // check input
    if n < 0 || A.is_empty() || lda < n {
        return kRetInvalidInput;
    }

    // quick return
    if n == 0 {
        return kRetOk;
    }

    // HIPO_CLOCK_CREATE dropped (line 164).

    if uplo == b'L' {
        debug_assert!(false); // assert(1 == 0) at line 167
    } else {
        if swaps.is_empty() || pivot_2x2.is_empty() {
            return kRetInvalidInput;
        }

        // initialise order of pivots
        for i in 0..n {
            swaps[i as usize] = i;
        }

        // allocate space for copy of col(s)
        let mut temp = vec![0.0f64; (n - 1) as usize];
        let mut temp2 = vec![0.0f64; (n - 1) as usize];
        // Rust-only scratch for the dger `y` row (see the PORT NOTEs below).
        let mut ybuf = vec![0.0f64; (n - 1) as usize];

        let mut step: Int;

        let mut j: Int = 0;
        while j < n {
            // #ifdef HIPO_PIVOTING (active)
            let flag_2x2 = blockBunchKaufman(
                j, n, A, lda, swaps, pivot_sign, thresh, regval, totalreg, data,
            );

            // cannot do 2x2 pivoting on last column
            debug_assert!(j < n - 1 || !flag_2x2);

            if !flag_2x2 {
                // 1x1 pivots
                step = 1;

                // diagonal element
                let mut Ajj = A[j as usize + lda as usize * j as usize];

                if Ajj.is_nan() {
                    return kRetInvalidPivot;
                }

                // add regularisation
                staticReg(
                    &mut Ajj,
                    pivot_sign[j as usize],
                    regval,
                    &mut totalreg[j as usize],
                );

                // #ifndef HIPO_PIVOTING branch dropped (lines 203-207)

                // save reciprocal of pivot
                A[j as usize + lda as usize * j as usize] = 1.0 / Ajj;

                let M = n - j - 1;
                if M > 0 {
                    // make copy of row
                    let off = j as usize + (j + 1) as usize * lda as usize;
                    callAndTime_dcopy(M, &A[off..], lda, &mut temp[..M as usize], 1, data);

                    // scale row j
                    let off = j as usize + (j + 1) as usize * lda as usize;
                    callAndTime_dscal(M, 1.0 / Ajj, &mut A[off..], lda, data);

                    // update rest of the matrix
                    // PORT NOTE: C++ passes y = &A[j + (j + 1) * lda] (stride
                    // lda) and the output &A[j + 1 + (j + 1) * lda] into the
                    // same array; row j interleaves with the output columns,
                    // so no split_at_mut separates them. dger only writes rows
                    // j+1..n-1, never row j, so y is gathered unchanged into
                    // ybuf and passed with stride 1: same values, same order.
                    let y_off = j as usize + (j + 1) as usize * lda as usize;
                    gatherStrided(&A[..], y_off, lda, M, &mut ybuf[..]);
                    let out_off = (j + 1) as usize + (j + 1) as usize * lda as usize;
                    callAndTime_dger(
                        M,
                        M,
                        -1.0,
                        &temp[..M as usize],
                        1,
                        &ybuf[..M as usize],
                        1,
                        &mut A[out_off..],
                        lda,
                        data,
                    );
                }
            } else {
                // 2x2 pivots
                step = 2;

                // diagonal pivot elements
                let d1 = A[j as usize + lda as usize * j as usize];
                let d2 = A[(j + 1) as usize + lda as usize * (j + 1) as usize];

                if d1.is_nan() || d2.is_nan() {
                    return kRetInvalidPivot;
                }

                // off-diagonal pivot element
                let offd = A[j as usize + lda as usize * (j + 1) as usize];
                A[j as usize + lda as usize * (j + 1) as usize] = 0.0;

                // compute coefficients of 2x2 inverse
                let denom = d1 * d2 - offd * offd;
                let i_d1 = d2 / denom;
                let i_d2 = d1 / denom;
                let i_off = -offd / denom;

                // save them in place of pivots
                A[j as usize + lda as usize * j as usize] = i_d1;
                A[(j + 1) as usize + lda as usize * (j + 1) as usize] = i_d2;
                pivot_2x2[j as usize] = i_off;

                let M = n - j - 2;
                if M > 0 {
                    let r1_off = j as usize + (j + 2) as usize * lda as usize;
                    let r2_off = (j + 1) as usize + (j + 2) as usize * lda as usize;

                    // make a copy of first row
                    callAndTime_dcopy(M, &A[r1_off..], lda, &mut temp[..M as usize], 1, data);

                    // make a copy of second row
                    callAndTime_dcopy(M, &A[r2_off..], lda, &mut temp2[..M as usize], 1, data);

                    // solve rows j,j+1
                    callAndTime_dscal(M, i_d1, &mut A[r1_off..], lda, data);
                    callAndTime_daxpy(
                        M,
                        i_off,
                        &temp2[..M as usize],
                        1,
                        &mut A[r1_off..],
                        lda,
                        data,
                    );
                    callAndTime_dscal(M, i_d2, &mut A[r2_off..], lda, data);
                    callAndTime_daxpy(
                        M,
                        i_off,
                        &temp[..M as usize],
                        1,
                        &mut A[r2_off..],
                        lda,
                        data,
                    );

                    // update rest of the matrix
                    // PORT NOTE: as in the 1x1 branch, r1 / r2 (rows j, j+1,
                    // stride lda) interleave with the output block at
                    // (j+2, j+2); dger writes only rows j+2..n-1, so each row
                    // is gathered unchanged into ybuf right before its call.
                    let out_off = (j + 2) as usize + (j + 2) as usize * lda as usize;
                    gatherStrided(&A[..], r1_off, lda, M, &mut ybuf[..]);
                    callAndTime_dger(
                        M,
                        M,
                        -1.0,
                        &temp[..M as usize],
                        1,
                        &ybuf[..M as usize],
                        1,
                        &mut A[out_off..],
                        lda,
                        data,
                    );
                    gatherStrided(&A[..], r2_off, lda, M, &mut ybuf[..]);
                    callAndTime_dger(
                        M,
                        M,
                        -1.0,
                        &temp2[..M as usize],
                        1,
                        &ybuf[..M as usize],
                        1,
                        &mut A[out_off..],
                        lda,
                        data,
                    );
                }

                data.count2x2();
            }

            j += step;
        }
    }

    // HIPO_CLOCK_STOP(2, data, kTimeDenseFact_kernel) dropped (line 279).
    kRetOk
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn factor_spd_2x2() {
        // Upper format (C++ `uplo == 'U'`, DenseFactKernel.cpp 169-281). The
        // 'L' branch is `assert(1 == 0)` with no body (lines 166-167), so under
        // NDEBUG it leaves A untouched: the skeleton's 'L' test contradicted the
        // C++. Column-major a = [4 0 2 3] is the upper matrix [4 2; 0 3]:
        // D = [4, 2], U12 = 2 / 4 = 0.5, and the Schur complement is 3 - 0.5 * 2 = 2.
        // The diagonal keeps the reciprocal of each 1x1 pivot (line 210): 1/4 and 1/2.
        let mut a = [4.0, 0.0, 2.0, 3.0];
        let mut sign = [1, 1];
        let mut reg = [0.0; 2];
        let mut swaps = [0, 1];
        let mut p2 = [0.0; 2];
        let d = DataCollector::new();
        let r = denseFactK(
            b'U',
            2,
            &mut a,
            2,
            &mut sign,
            1e-10,
            &Regul::default(),
            &mut reg,
            &mut swaps,
            &mut p2,
            &d,
        );
        assert_eq!(r, kRetOk);
        assert_eq!(a[0], 0.25);
        assert_eq!(a[2], 0.5);
        assert_eq!(a[3], 0.5);
    }

    #[test]
    fn factor_1x1_n3_reaches_dger_with_m2() {
        // Upper, lda = 3, symmetric [[4,2,1],[2,3,0.5],[1,0.5,2]] (lower entries 0).
        // j=0: no swap, gamma_j=2, |Ajj|=4 >= kAlphaBK*2 -> accept; 1/4 on the
        // diagonal, row scaled by 1/4, dger with temp=[2,1], y=[0.5,0.25] updates
        // rows 1..2 of columns 1..2 (lower (2,1) becomes -0.5).
        // j=1: gamma_j=|A(1,2)|=0 -> accept; 1/2, A(1,2)=0, dger skipped (y=0).
        // j=2: last column, must accept; 1/1.75.
        let mut a = [4.0, 0.0, 0.0, 2.0, 3.0, 0.0, 1.0, 0.5, 2.0];
        let mut sign = [1, 1, 1];
        let mut reg = [0.0; 3];
        let mut swaps = [0; 3];
        let mut p2 = [0.0; 3];
        let d = DataCollector::new();
        let r = denseFactK(
            b'U',
            3,
            &mut a,
            3,
            &mut sign,
            1e-10,
            &Regul::default(),
            &mut reg,
            &mut swaps,
            &mut p2,
            &d,
        );
        assert_eq!(r, kRetOk);
        assert_eq!(a, [0.25, 0.0, 0.0, 0.5, 0.5, -0.5, 0.25, 0.0, 1.0 / 1.75]);
        assert_eq!(swaps, [0, 1, 2]);
        assert_eq!(reg, [0.0; 3]);
    }

    #[test]
    fn factor_2x2_pivot_n4() {
        // Upper, lda = 4, zero diagonal and ones above it, signs [1,-1,-1,-1].
        // j=0: all diagonals 0 -> no swap; gamma_j=1 (r=1), gamma_r=1, Ajj=Arr=0
        // fail both BK tests -> 2x2 pivot (0,1) (lines 129-140). denom=-1, so
        // i_d1=i_d2=-0.0, i_off=1. Rows 0,1 of columns 2..3 become [1,1]; the two
        // dger calls (M=2) subtract 1 twice from the block rows 2..3, cols 2..3.
        // j=2: diag -2,-2 -> no swap; gamma_j=|A(2,3)|=1, |Ajj|=2 -> accept;
        // A(2,2)=-0.5, A(2,3)=-1*(1/-2)=0.5, A(3,3)=-2+(-1)*(-0.5)=-1.5.
        // j=3: last column -> 1/-1.5.
        let mut a = [
            0.0, 0.0, 0.0, 0.0, // col 0
            1.0, 0.0, 0.0, 0.0, // col 1
            1.0, 1.0, 0.0, 0.0, // col 2
            1.0, 1.0, 1.0, 0.0, // col 3
        ];
        let mut sign = [1, -1, -1, -1];
        let mut reg = [0.0; 4];
        let mut swaps = [0; 4];
        let mut p2 = [0.0; 4];
        let d = DataCollector::new();
        let r = denseFactK(
            b'U',
            4,
            &mut a,
            4,
            &mut sign,
            1e-10,
            &Regul::default(),
            &mut reg,
            &mut swaps,
            &mut p2,
            &d,
        );
        assert_eq!(r, kRetOk);
        let expected: [f64; 16] = [
            -0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            -0.5,
            -2.0,
            1.0,
            1.0,
            0.5,
            1.0 / -1.5,
        ];
        for k in 0..16 {
            assert_eq!(a[k].to_bits(), expected[k].to_bits(), "a[{k}]");
        }
        assert_eq!(p2[0], 1.0);
        assert_eq!(swaps, [0, 1, 2, 3]);
    }

    #[test]
    fn factor_bk_uses_pivot_r() {
        // Ajj adds regval.dual for sign > 0 but Arr adds regval.primal (lines
        // 89-90, 111-112), so with primal=1, dual=0 and a zero diagonal the
        // "use pivot r" branch (lines 120-127) is taken: swapCols(0,1) sets
        // swaps[0]=1, staticReg makes A(0,0)=1 (totalreg 1); the 1x1 step adds
        // primal again (line 200): 1/2. A(0,1)=0.5, A(1,1)=0-0.5. j=1 is last:
        // -0.5 -> -0.5+1=0.5, totalreg[1]=0.5-(-0.5)=1; 1x1 step: 1/1.5, totalreg 1.
        let mut a = [0.0, 0.0, 1.0, 0.0];
        let mut sign = [1, 1];
        let mut reg = [0.0; 2];
        let mut swaps = [0; 2];
        let mut p2 = [0.0; 2];
        let d = DataCollector::new();
        let regval = Regul {
            primal: 1.0,
            dual: 0.0,
        };
        let r = denseFactK(
            b'U', 2, &mut a, 2, &mut sign, 1e-10, &regval, &mut reg, &mut swaps, &mut p2, &d,
        );
        assert_eq!(r, kRetOk);
        assert_eq!(a, [0.5, 0.0, 0.5, 1.0 / 1.5]);
        assert_eq!(swaps, [1, 1]);
        assert_eq!(reg, [1.0, 1.0]);
    }

    #[test]
    fn factor_small_pivot_perturbed() {
        // Must-accept branch with max(|Ajj|, gamma_j) < thresh (lines 100-104):
        // the pivot is set to sign*thresh, then the 1x1 step regularises a local
        // copy again (line 200) and overwrites totalreg[0].
        let thresh = 1e-10;
        let regval = Regul {
            primal: 1e-14,
            dual: 0.0,
        };
        let mut a = [1e-12];
        let mut sign = [1];
        let mut reg = [0.0];
        let mut swaps = [0];
        let mut p2 = [0.0];
        let d = DataCollector::new();
        let r = denseFactK(
            b'U', 1, &mut a, 1, &mut sign, thresh, &regval, &mut reg, &mut swaps, &mut p2, &d,
        );
        assert_eq!(r, kRetOk);
        let ajj = thresh + regval.primal;
        assert_eq!(a[0], 1.0 / ajj);
        assert_eq!(reg[0], ajj - thresh);
    }

    #[test]
    fn factor_wrong_sign_last_column() {
        // sign = -1 on a positive pivot: sign*Ajj < 0 forces the must-accept
        // branch (line 92); staticReg subtracts dual: 2 - 0.5 = 1.5.
        // 1x1 step subtracts dual again: 1.0, so A = 1/1.
        let regval = Regul {
            primal: 0.0,
            dual: 0.5,
        };
        let mut a = [2.0];
        let mut sign = [-1];
        let mut reg = [0.0];
        let mut swaps = [0];
        let mut p2 = [0.0];
        let d = DataCollector::new();
        let r = denseFactK(
            b'U', 1, &mut a, 1, &mut sign, 1e-10, &regval, &mut reg, &mut swaps, &mut p2, &d,
        );
        assert_eq!(r, kRetOk);
        assert_eq!(a[0], 1.0);
        assert_eq!(reg[0], -0.5);
    }

    #[test]
    fn factor_invalid_input() {
        // lda < n (line 158). (A NaN-only residual diagonal is not tested: the
        // C++ then calls swapCols(j, -1), UB under NDEBUG, lines 66-79.)
        let d = DataCollector::new();
        let r = denseFactK(
            b'U',
            2,
            &mut [0.0; 4],
            1,
            &mut [1, 1],
            1e-10,
            &Regul::default(),
            &mut [0.0; 2],
            &mut [0; 2],
            &mut [0.0; 2],
            &d,
        );
        assert_eq!(r, kRetInvalidInput);
    }
}
