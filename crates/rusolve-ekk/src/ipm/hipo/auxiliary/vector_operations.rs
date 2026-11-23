//! HiGHS source: ipm/hipo/auxiliary/VectorOperations.h, ipm/hipo/auxiliary/VectorOperations.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Dense vector kernels of the IPM (plain loops, no BLAS).

use super::super::Int;

/// HiGHS: VectorOperations.cpp::vectorAdd (lines 8-13)
/// `v1[i] = alpha * v1[i] + beta * v2[i]`. C++ default arguments: `beta = 1.0`, `alpha = 1.0`.
pub fn vectorAdd(v1: &mut Vec<f64>, v2: &[f64], beta: f64, alpha: f64) {
    for i in 0..v1.len() as Int {
        v1[i as usize] = alpha * v1[i as usize] + beta * v2[i as usize];
    }
}

/// HiGHS: VectorOperations.cpp::vectorAdd (lines 15-19)
/// Scalar overload: `v1 += alpha`.
pub fn vectorAdd_scalar(v1: &mut Vec<f64>, alpha: f64) {
    for i in 0..v1.len() as Int {
        v1[i as usize] += alpha;
    }
}

/// HiGHS: VectorOperations.cpp::vectorMultiply (lines 21-26)
/// C++ default arguments (positional; the header names them `beta = 1.0`,
/// `alpha = 0.0`, the .cpp `alpha`, `beta`): first `1.0`, second `0.0`.
pub fn vectorMultiply(v1: &mut Vec<f64>, v2: &[f64], alpha: f64, beta: f64) {
    for i in 0..v1.len() as Int {
        v1[i as usize] = alpha * v1[i as usize] * v2[i as usize] + beta;
    }
}

/// HiGHS: VectorOperations.cpp::vectorAddMult (lines 28-33)
/// C++ default arguments: `beta = 1.0`, `alpha = 1.0`.
pub fn vectorAddMult(v1: &mut Vec<f64>, v2: &[f64], v3: &[f64], beta: f64, alpha: f64) {
    for i in 0..v1.len() as Int {
        v1[i as usize] = alpha * v1[i as usize] + beta * v2[i as usize] * v3[i as usize];
    }
}

/// HiGHS: VectorOperations.cpp::vectorDivide (lines 35-39)
pub fn vectorDivide(v1: &mut Vec<f64>, v2: &[f64]) {
    for i in 0..v1.len() as Int {
        v1[i as usize] /= v2[i as usize];
    }
}

/// HiGHS: VectorOperations.cpp::vectorScale (lines 41-45)
pub fn vectorScale(v1: &mut Vec<f64>, alpha: f64) {
    for i in 0..v1.len() as Int {
        v1[i as usize] *= alpha;
    }
}

/// HiGHS: VectorOperations.cpp::dotProd (lines 47-53)
pub fn dotProd(v1: &[f64], v2: &[f64]) -> f64 {
    let mut result: f64 = 0.0;
    for i in 0..v1.len() as Int {
        result += v1[i as usize] * v2[i as usize];
    }
    result
}

/// HiGHS: VectorOperations.cpp::norm2 (lines 55-61)
pub fn norm2(x: &[f64]) -> f64 {
    let mut norm: f64 = 0.0;
    for i in 0..x.len() as Int {
        norm += x[i as usize] * x[i as usize];
    }
    norm.sqrt()
}

/// HiGHS: VectorOperations.cpp::norm2 (lines 63-72)
/// Two-vector overload: norm of the concatenation.
pub fn norm2_2(x0: &[f64], x1: &[f64]) -> f64 {
    let mut norm: f64 = 0.0;
    for i in 0..x0.len() as Int {
        norm += x0[i as usize] * x0[i as usize];
    }
    for i in 0..x1.len() as Int {
        norm += x1[i as usize] * x1[i as usize];
    }
    norm.sqrt()
}

/// HiGHS: VectorOperations.cpp::infNorm (lines 74-80)
pub fn infNorm(x: &[f64]) -> f64 {
    let mut norm: f64 = 0.0;
    for i in 0..x.len() as Int {
        norm = crate::util::highs_utils::cpp_max(norm, x[i as usize].abs());
    }
    norm
}

/// HiGHS: VectorOperations.cpp::infNorm (lines 82-91)
/// Two-vector overload.
pub fn infNorm_2(x0: &[f64], x1: &[f64]) -> f64 {
    let mut norm: f64 = 0.0;
    for i in 0..x0.len() as Int {
        norm = crate::util::highs_utils::cpp_max(norm, x0[i as usize].abs());
    }
    for i in 0..x1.len() as Int {
        norm = crate::util::highs_utils::cpp_max(norm, x1[i as usize].abs());
    }
    norm
}

/// HiGHS: VectorOperations.cpp::infNormDiff (lines 93-101)
pub fn infNormDiff(x: &[f64], y: &[f64]) -> f64 {
    debug_assert!(x.len() == y.len());
    let mut inf_norm_diff: f64 = 0.0;
    for i in 0..x.len() as Int {
        let diff = (x[i as usize] - y[i as usize]).abs();
        inf_norm_diff = crate::util::highs_utils::cpp_max(diff, inf_norm_diff);
    }
    inf_norm_diff
}

/// HiGHS: VectorOperations.cpp::isNanVector (lines 103-108)
pub fn isNanVector(x: &[f64]) -> bool {
    for i in 0..x.len() as Int {
        if x[i as usize].is_nan() {
            return true;
        }
    }
    false
}

/// HiGHS: VectorOperations.cpp::isInfVector (lines 110-115)
pub fn isInfVector(x: &[f64]) -> bool {
    for i in 0..x.len() as Int {
        if x[i as usize].is_infinite() {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn kernels() {
        let mut v = vec![1.0, 2.0];
        vectorAdd(&mut v, &[1.0, 1.0], 2.0, 1.0);
        assert_eq!(v, vec![3.0, 4.0]);
        assert_eq!(norm2(&[3.0, 4.0]), 5.0);
        assert_eq!(infNorm_2(&[1.0, -7.0], &[3.0]), 7.0);
        assert_eq!(dotProd(&[1.0, 2.0], &[3.0, 4.0]), 11.0);
        assert!(isNanVector(&[0.0, f64::NAN]));
    }
}
