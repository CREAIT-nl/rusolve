//! HiGHS source: pdlp/hipdlp/linalg.hpp, pdlp/hipdlp/linalg.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `namespace linalg` of HiPDLP. Reached by path (`hipdlp::linalg::ax`):
//! the names collide with the cuPDLP and simplex helpers.
//!
//! `std::max` / `std::min` are `cpp_max` / `cpp_min`; `std::pow` is `powf`;
//! sums are left to right as written (no FMA).

use crate::prelude::*;

/// HiGHS: linalg.cc::linalg::projectBox (lines 17-19)
pub fn projectBox(x: f64, l: f64, u: f64) -> f64 {
    cpp_max(l, cpp_min(x, u))
}

/// HiGHS: linalg.cc::linalg::projectNonNegative (line 21)
pub fn projectNonNegative(y: f64) -> f64 {
    cpp_max(0.0, y)
}

/// HiGHS: linalg.cc::linalg::projectBounds (lines 23-34)
pub fn projectBounds(lp: &HighsLp, x: &mut [f64]) {
    for i in 0..lp.num_col_ {
        // Project to upper bound
        if x[i as usize] > lp.col_upper_[i as usize] {
            x[i as usize] = lp.col_upper_[i as usize];
        }
        // Project to lower bound
        if x[i as usize] < lp.col_lower_[i as usize] {
            x[i as usize] = lp.col_lower_[i as usize];
        }
    }
}

/// HiGHS: linalg.cc::linalg::ax (lines 36-47)
///
/// `result` is filled with zero (not resized) then accumulated column by
/// column (CSC).
pub fn ax(lp: &HighsLp, x: &[f64], result: &mut [f64]) {
    result.fill(0.0);
    // Assumes column-wise matrix format
    for col in 0..lp.num_col_ {
        for i in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
            let row = lp.a_matrix_.index_[i as usize];
            result[row as usize] += lp.a_matrix_.value_[i as usize] * x[col as usize];
        }
    }
}

/// HiGHS: linalg.cc::linalg::aTy (lines 49-61)
pub fn aTy(lp: &HighsLp, y: &[f64], result: &mut [f64]) {
    result.fill(0.0);
    // Assumes column-wise matrix format. For each column `col` of A,
    // this loop calculates Dot(column_col, y) and adds it to result[col].
    for col in 0..lp.num_col_ {
        for i in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
            let row = lp.a_matrix_.index_[i as usize];
            result[col as usize] += lp.a_matrix_.value_[i as usize] * y[row as usize];
        }
    }
}

/// HiGHS: linalg.cc::linalg::dot (lines 63-73)
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    assert!(
        a.len() == b.len(),
        "Vectors must be of the same size for Dot product."
    );
    let mut result = 0.0;
    for i in 0..a.len() {
        result += a[i] * b[i];
    }
    result
}

/// HiGHS: linalg.cc::linalg::norm2 (lines 79-81)
pub fn norm2(vec: &[f64]) -> f64 {
    dot(vec, vec).sqrt()
}

/// HiGHS: linalg.cc::linalg::vectorNormSquared (lines 83-85)
pub fn vectorNormSquared(vec: &[f64]) -> f64 {
    dot(vec, vec)
}

/// HiGHS: linalg.cc::linalg::vectorNorm(const double*, size_t, double) (lines 87-108)
///
/// C++ default arguments: `p = 2.0`.
pub fn vectorNorm_ptr(values: &[f64], size: usize, p: f64) -> f64 {
    if p.is_infinite() {
        // Infinity norm
        let mut max_val = 0.0;
        for i in 0..size {
            max_val = cpp_max(max_val, values[i].abs());
        }
        return max_val;
    }
    if p == 1.0 {
        // L1 norm
        let mut sum = 0.0;
        for i in 0..size {
            sum += values[i].abs();
        }
        return sum;
    }
    // General p-norm (including L2)
    let mut sum = 0.0;
    for i in 0..size {
        sum += values[i].abs().powf(p);
    }
    sum.powf(1.0 / p)
}

/// HiGHS: linalg.cc::linalg::vectorNorm(const std::vector<double>&, double) (lines 110-116)
///
/// C++ default arguments: `p = 2.0`.
pub fn vectorNorm(vec: &[f64], p: f64) -> f64 {
    // For L2 norm, call the optimized Norm2 to avoid pow()
    if p == 2.0 {
        return norm2(vec);
    }
    vectorNorm_ptr(vec, vec.len(), p)
}

/// HiGHS: linalg.cc::linalg::scale (lines 122-126)
pub fn scale(vec: &mut [f64], factor: f64) {
    for val in vec.iter_mut() {
        *val *= factor;
    }
}

/// HiGHS: linalg.cc::linalg::normalize (lines 128-133)
pub fn normalize(vec: &mut [f64]) {
    let norm = norm2(vec);
    if norm > 1e-9 {
        // Use a small tolerance
        scale(vec, 1.0 / norm);
    }
}

/// HiGHS: linalg.cc::linalg::diffTwoNorm (lines 135-146)
pub fn diffTwoNorm(v1: &[f64], v2: &[f64]) -> f64 {
    assert!(
        v1.len() == v2.len(),
        "Vectors must have the same dimension."
    );
    let mut norm_sq = 0.0;
    for i in 0..v1.len() {
        let diff = v1[i] - v2[i];
        norm_sq += diff * diff;
    }
    norm_sq.sqrt()
}

/// HiGHS: linalg.cc::linalg::computeCostNorm (lines 152-154)
///
/// C++ default arguments: `p = 2.0`.
pub fn computeCostNorm(lp: &HighsLp, p: f64) -> f64 {
    vectorNorm(&lp.col_cost_, p)
}

/// HiGHS: linalg.cc::linalg::computeRhsNorm (lines 156-178)
///
/// C++ default arguments: `p = 2.0`.
pub fn computeRhsNorm(lp: &HighsLp, p: f64) -> f64 {
    if p.is_infinite() {
        let mut max_val = 0.0;
        for val in lp.row_lower_.iter() {
            if val.is_finite() {
                max_val = cpp_max(max_val, val.abs());
            }
        }
        return max_val;
    }

    let mut sum = 0.0;
    for val in lp.row_lower_.iter() {
        if val.is_finite() {
            if p == 1.0 {
                sum += val.abs();
            } else {
                sum += val.abs().powf(p);
            }
        }
    }

    if p == 2.0 {
        return sum.sqrt();
    }
    if p == 1.0 {
        return sum;
    }
    sum.powf(1.0 / p)
}

/// HiGHS: linalg.cc::linalg::computeColumnNorms (lines 180-190)
///
/// C++ default arguments: `p = std::numeric_limits<double>::infinity()`.
pub fn computeColumnNorms(lp: &HighsLp, p: f64) -> Vec<f64> {
    let mut col_norms = vec![0.0; lp.num_col_ as usize];
    for col in 0..lp.num_col_ {
        let start = lp.a_matrix_.start_[col as usize];
        let end = lp.a_matrix_.start_[(col + 1) as usize];
        if start < end {
            col_norms[col as usize] = vectorNorm_ptr(
                &lp.a_matrix_.value_[start as usize..end as usize],
                (end - start) as usize,
                p,
            );
        }
    }
    col_norms
}

/// HiGHS: linalg.cc::linalg::computeRowNorms (lines 192-216)
///
/// C++ default arguments: `p = std::numeric_limits<double>::infinity()`.
pub fn computeRowNorms(lp: &HighsLp, p: f64) -> Vec<f64> {
    let mut row_norms = vec![0.0; lp.num_row_ as usize];
    for col in 0..lp.num_col_ {
        for i in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
            let row = lp.a_matrix_.index_[i as usize];
            let abs_val = lp.a_matrix_.value_[i as usize].abs();
            if p.is_infinite() {
                row_norms[row as usize] = cpp_max(row_norms[row as usize], abs_val);
            } else if p == 1.0 {
                row_norms[row as usize] += abs_val;
            } else {
                row_norms[row as usize] += abs_val.powf(p);
            }
        }
    }

    if p != 1.0 && !p.is_infinite() {
        let p_inv = 1.0 / p;
        for row in 0..lp.num_row_ {
            row_norms[row as usize] = row_norms[row as usize].powf(p_inv);
        }
    }
    row_norms
}

/// HiGHS: linalg.cc::linalg::vectorSubtract (lines 218-229)
pub fn vectorSubtract(a: &[f64], b: &[f64]) -> Vec<f64> {
    assert!(
        a.len() == b.len(),
        "Vectors must be of the same size for subtraction."
    );
    let mut result = vec![0.0; a.len()];
    for i in 0..a.len() {
        result[i] = a[i] - b[i];
    }
    result
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn projections_and_norms() {
        assert_eq!(projectBox(5.0, 0.0, 2.0), 2.0);
        assert_eq!(projectNonNegative(-1.0), 0.0);
        assert_eq!(norm2(&[3.0, 4.0]), 5.0);
        assert_eq!(vectorNorm(&[3.0, -4.0], f64::INFINITY), 4.0);
        assert_eq!(vectorNorm(&[3.0, -4.0], 1.0), 7.0);
        assert_eq!(diffTwoNorm(&[1.0, 1.0], &[4.0, 5.0]), 5.0);
        assert_eq!(vectorSubtract(&[3.0, 2.0], &[1.0, 1.0]), vec![2.0, 1.0]);
    }
}
