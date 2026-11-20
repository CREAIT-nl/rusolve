//! HiGHS source: ipm/hipo/factorhighs/Swaps.h, ipm/hipo/factorhighs/Swaps.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Row/column interchanges of the Bunch-Kaufman pivoting.

use super::super::Int;
use super::data_collector::DataCollector;

/// HiGHS: Swaps.cpp::permuteWithSwaps (lines 8-22)
/// C++ default arguments: `reverse = false`.
pub fn permuteWithSwaps(x: &mut [f64], swaps: &[Int], n: Int, reverse: bool) {
    // Apply swaps to vector x of length n

    if !reverse {
        // apply the swaps in forward order
        let mut i: Int = 0;
        while i < n {
            if swaps[i as usize] != i {
                x.swap(i as usize, swaps[i as usize] as usize);
            }
            i += 1;
        }
    } else {
        // apply the swaps in backward order
        let mut i: Int = n - 1;
        while i >= 0 {
            if swaps[i as usize] != i {
                x.swap(i as usize, swaps[i as usize] as usize);
            }
            i -= 1;
        }
    }
}

/// HiGHS: Swaps.cpp::swapCols (lines 24-57)
pub fn swapCols(
    uplo: u8,
    n: Int,
    A: &mut [f64],
    lda: Int,
    i: Int,
    j: Int,
    swaps: &mut [Int],
    sign: &mut [Int],
    data: &DataCollector,
) {
    // Exchange rows/cols i and j of symmetric matrix A

    // make sure that i < j
    if i == j {
        return;
    }
    let (mut i, mut j) = (i, j);
    if i > j {
        std::mem::swap(&mut i, &mut j);
    }

    // swap diagonal elements
    A.swap((i + i * lda) as usize, (j + j * lda) as usize);

    // swap rest of rows/cols
    if uplo == b'L' {
        // PORT NOTE: C++ passes two pointers into the same array A to
        // callAndTime_dswap (strided rows/columns); Rust cannot express the
        // aliasing as two &mut slices, so the swap is done in place on A with
        // the same element order as refblas dswap. C lines 37-41.
        dswap_same_array(A, i, (i) as usize, lda, (j) as usize, lda);
        dswap_same_array(
            A,
            n - j - 1,
            (j + 1 + i * lda) as usize,
            1,
            (j + 1 + j * lda) as usize,
            1,
        );
        dswap_same_array(
            A,
            j - i - 1,
            (i + 1 + i * lda) as usize,
            1,
            (j + (i + 1) * lda) as usize,
            lda,
        );
    } else {
        // PORT NOTE: same aliasing as above. C lines 43-47.
        dswap_same_array(A, i, (i * lda) as usize, 1, (j * lda) as usize, 1);
        dswap_same_array(
            A,
            n - j - 1,
            (i + (j + 1) * lda) as usize,
            lda,
            (j + (j + 1) * lda) as usize,
            lda,
        );
        dswap_same_array(
            A,
            j - i - 1,
            (i + (i + 1) * lda) as usize,
            lda,
            (i + 1 + j * lda) as usize,
            1,
        );
    }

    // swap pivot sign
    sign.swap(i as usize, j as usize);

    // keep track of order of swaps
    swaps[i as usize] = j;

    data.countSwap();
}

/// HiGHS: Swaps.cpp::applySwaps (lines 59-68)
pub fn applySwaps(swaps: &[Int], nrow: Int, ncol: Int, R: &mut [f64], data: &DataCollector) {
    // apply the column swaps to block R
    let mut i: Int = 0;
    while i < ncol {
        if swaps[i as usize] != i {
            // swap col i and col swaps[i]
            // PORT NOTE: C++ passes two pointers into the same array R to
            // callAndTime_dswap (column stride ncol); Rust cannot express the
            // aliasing as two &mut slices, so the swap is done in place on R
            // with the same element order as refblas dswap. C line 65.
            dswap_same_array(
                R,
                nrow,
                (i) as usize,
                ncol,
                (swaps[i as usize]) as usize,
                ncol,
            );
        }
        i += 1;
    }
}

/// PORT NOTE: private helper replacing callAndTime_dswap when both pointers
/// alias the same array (C++ Swaps.cpp lines 37-41, 43-47, 65). Performs
/// exactly the refblas dswap element swaps (refblas.c lines 93-107) in index
/// order. All callers pass positive increments (`1`, `lda`, `ncol`), for which
/// refblas starts at offset 0, and the two strided sets never overlap, so the
/// result is the same permutation of `a` as the C++. `kTimeBlas_swap` timing
/// is a no-op (timing level 0).
fn dswap_same_array(a: &mut [f64], n: Int, ox: usize, incx: Int, oy: usize, incy: Int) {
    if n <= 0 {
        return;
    }
    let mut k: Int = 0;
    while k < n {
        let ix = (ox as i64 + (k as i64) * (incx as i64)) as usize;
        let iy = (oy as i64 + (k as i64) * (incy as i64)) as usize;
        a.swap(ix, iy);
        k += 1;
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn permute_with_swaps_roundtrip() {
        let mut x = [1.0, 2.0, 3.0];
        let swaps = [2, 1, 2];
        permuteWithSwaps(&mut x, &swaps, 3, false);
        permuteWithSwaps(&mut x, &swaps, 3, true);
        assert_eq!(x, [1.0, 2.0, 3.0]);
    }

    /// Full symmetric 4x4 matrix, column-major: swapping rows/cols i and j of
    /// the stored triangle must equal the triangle of P A P.
    fn sym4() -> [[f64; 4]; 4] {
        let mut m = [[0.0; 4]; 4];
        for r in 0..4 {
            for c in 0..4 {
                let (a, b) = if r >= c { (r, c) } else { (c, r) };
                m[r][c] = (10 * a + b) as f64 + 1.0;
            }
        }
        m
    }

    fn check_swap_cols(uplo: u8, i: Int, j: Int) {
        let n = 4;
        let m = sym4();
        let mut A = vec![0.0; 16];
        for c in 0..4 {
            for r in 0..4 {
                if (uplo == b'L' && r >= c) || (uplo == b'U' && r <= c) {
                    A[r + 4 * c] = m[r][c];
                }
            }
        }
        let mut swaps = [0, 1, 2, 3];
        let mut sign = [1, -1, 1, -1];
        let d = DataCollector::new();
        swapCols(uplo, n, &mut A, n, i, j, &mut swaps, &mut sign, &d);
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        let p = |k: usize| {
            if k == lo as usize {
                hi as usize
            } else if k == hi as usize {
                lo as usize
            } else {
                k
            }
        };
        for c in 0..4 {
            for r in 0..4 {
                let stored = (uplo == b'L' && r >= c) || (uplo == b'U' && r <= c);
                let expect = if stored { m[p(r)][p(c)] } else { 0.0 };
                assert_eq!(
                    A[r + 4 * c],
                    expect,
                    "uplo {} r {} c {}",
                    uplo as char,
                    r,
                    c
                );
            }
        }
        let mut sign_expect = [1, -1, 1, -1];
        sign_expect.swap(lo as usize, hi as usize);
        assert_eq!(sign, sign_expect);
        assert_eq!(swaps[lo as usize], hi);
    }

    #[test]
    fn swap_cols_lower_and_upper() {
        check_swap_cols(b'L', 1, 3);
        check_swap_cols(b'L', 3, 0);
        check_swap_cols(b'U', 1, 3);
        check_swap_cols(b'U', 2, 0);
        // i == j: early return, nothing changes
        let mut A = vec![1.0, 2.0, 3.0, 4.0];
        let mut swaps = [0, 1];
        let mut sign = [1, -1];
        swapCols(
            b'L',
            2,
            &mut A,
            2,
            1,
            1,
            &mut swaps,
            &mut sign,
            &DataCollector::new(),
        );
        assert_eq!(A, [1.0, 2.0, 3.0, 4.0]);
        assert_eq!(swaps, [0, 1]);
        assert_eq!(sign, [1, -1]);
    }

    #[test]
    fn apply_swaps_row_major_block() {
        // R is nrow x ncol with stride ncol: entries R[i + r * ncol].
        let (nrow, ncol) = (2, 3);
        let mut R = vec![0.0, 1.0, 2.0, 10.0, 11.0, 12.0];
        let swaps = [2, 1, 2];
        applySwaps(&swaps, nrow, ncol, &mut R, &DataCollector::new());
        assert_eq!(R, [2.0, 1.0, 0.0, 12.0, 11.0, 10.0]);
    }
}
