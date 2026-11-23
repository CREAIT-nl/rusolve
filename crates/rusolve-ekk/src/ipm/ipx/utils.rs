//! HiGHS source: ipm/ipx/utils.h, utils.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `Sortperm` uses HiGHS pdqsort with a strict (value, index) pair order:
//! translate with `util::highs_pdqsort::pdqsort_by` (not `sort_by`, which
//! differs on NaN values).

use super::*;
use crate::util::highs_pdqsort::pdqsort_by;
use crate::util::highs_utils::cpp_max;

/// libc++ `std::pair` `operator<`: `a.first < b.first || (!(b.first < a.first) && a.second < b.second)`.
///
/// PORT NOTE: C++ `std::pair` comparison helper shared by `Sortperm` and
/// `SparseMatrix::SortIndices`. Written out so that NaN values compare exactly
/// as the C++ does (the second field is consulted when neither first field is
/// less than the other).
pub fn pair_less<A: PartialOrd, B: PartialOrd>(a: (A, B), b: (A, B)) -> bool {
    a.0 < b.0 || (!(b.0 < a.0) && a.1 < b.1)
}

/// HiGHS: ipx/utils.cc::AllFinite (lines 10-15)
pub fn all_finite(x: &[f64]) -> bool {
    for &xi in x {
        if !xi.is_finite() {
            return false;
        }
    }
    true
}

/// HiGHS: ipx/utils.cc::Onenorm(const Vector&) (lines 17-22)
pub fn onenorm(x: &[f64]) -> f64 {
    let mut norm = 0.0;
    for &xi in x {
        norm += xi.abs();
    }
    norm
}

/// HiGHS: ipx/utils.cc::Twonorm (lines 24-29)
pub fn twonorm(x: &[f64]) -> f64 {
    let mut norm = 0.0;
    for &xi in x {
        norm += xi * xi;
    }
    norm.sqrt()
}

/// HiGHS: ipx/utils.cc::Infnorm(const Vector&) (lines 31-36)
pub fn infnorm(x: &[f64]) -> f64 {
    let mut norm = 0.0;
    for &xi in x {
        norm = cpp_max(norm, xi.abs());
    }
    norm
}

/// HiGHS: ipx/utils.cc::Dot(const Vector&, const Vector&) (lines 38-44)
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    debug_assert_eq!(x.len(), y.len());
    let mut d = 0.0;
    for (xi, yi) in x.iter().zip(y.iter()) {
        d += xi * yi;
    }
    d
}

/// HiGHS: ipx/utils.cc::FindMaxAbs (lines 46-56)
pub fn find_max_abs(x: &[f64]) -> Int {
    let mut xmax = 0.0;
    let mut imax: usize = 0;
    for (i, &xi) in x.iter().enumerate() {
        if xi.abs() > xmax {
            xmax = xi.abs();
            imax = i;
        }
    }
    imax as Int
}

/// HiGHS: ipx/utils.cc::Permute (lines 58-62)
pub fn permute(permuted_index: &[Int], rhs: &[f64], lhs: &mut [f64]) {
    for (i, &pi) in permuted_index.iter().enumerate() {
        lhs[pi as usize] = rhs[i];
    }
}

/// HiGHS: ipx/utils.cc::PermuteBack (lines 64-68)
pub fn permute_back(permuted_index: &[Int], rhs: &[f64], lhs: &mut [f64]) {
    for (i, &pi) in permuted_index.iter().enumerate() {
        lhs[i] = rhs[pi as usize];
    }
}

/// HiGHS: ipx/utils.cc::InversePerm (lines 70-76)
///
/// C++ `invperm.at(perm[i])` throws `std::out_of_range` on a bad index; the
/// Rust indexing panics in the same situation.
pub fn inverse_perm(perm: &[Int]) -> Vec<Int> {
    let mut invperm: Vec<Int> = vec![0; perm.len()];
    for (i, &p) in perm.iter().enumerate() {
        invperm[p as usize] = i as Int;
    }
    invperm
}

/// HiGHS: ipx/utils.cc::Sortperm (lines 78-93)
///
/// `values` may be `nullptr` (identity permutation).
pub fn sortperm(m: Int, values: Option<&[f64]>, reverse: bool) -> Vec<Int> {
    let mut perm: Vec<Int> = (0..m).collect();
    let values = match values {
        None => return perm,
        Some(v) => v,
    };
    if reverse {
        // make_pair(values[i], i) > make_pair(values[j], j)  <=>  pair(j) < pair(i)
        pdqsort_by(&mut perm, |&i, &j| {
            pair_less((values[j as usize], j), (values[i as usize], i))
        });
    } else {
        pdqsort_by(&mut perm, |&i, &j| {
            pair_less((values[i as usize], i), (values[j as usize], j))
        });
    }
    perm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sortperm_ties_and_reverse() {
        let values = [2.0, 1.0, 2.0, 0.0];
        // ascending by (value, index): 3(0.0), 1(1.0), 0(2.0), 2(2.0)
        assert_eq!(sortperm(4, Some(&values[..]), false), vec![3, 1, 0, 2]);
        // descending by (value, index): 2(2.0,idx 2) before 0(2.0,idx 0)
        assert_eq!(sortperm(4, Some(&values[..]), true), vec![2, 0, 1, 3]);
        assert_eq!(sortperm(3, None, false), vec![0, 1, 2]);
    }

    #[test]
    fn inverse_perm_roundtrip() {
        assert_eq!(inverse_perm(&[2, 0, 1]), vec![1, 2, 0]);
        let mut lhs = [0.0; 3];
        permute(&[2, 0, 1], &[10.0, 20.0, 30.0], &mut lhs);
        assert_eq!(lhs, [20.0, 30.0, 10.0]);
        let mut back = [0.0; 3];
        permute_back(&[2, 0, 1], &lhs, &mut back);
        assert_eq!(back, [10.0, 20.0, 30.0]);
    }
}
