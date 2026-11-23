//! HiGHS source: ipm/ipx/power_method.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::utils::twonorm;
use super::*;

/// HiGHS: ipx/power_method.h::PowerMethod (template, lines 13-44)
///
/// `func(v, fv)` computes `fv = A*v`. `v /= Twonorm(v)` and
/// `v = fv/lambda` are elementwise valarray expressions
/// (one division per element, in index order).
pub fn power_method<F: FnMut(&Vector, &mut Vector)>(func: F, v: &mut Vector) -> f64 {
    let mut func = func;
    let maxiter: Int = 100;
    let tol = 1e-3;
    let dim = v.size();
    let mut fv = Vector::new(dim);

    // Construct starting vector and normalize.
    for i in 0..dim {
        v[i as usize] = 1.0 + 1.0 / ((i + 1) as f64);
    }
    let nrm = twonorm(v);
    for i in 0..dim as usize {
        v[i] = v[i] / nrm;
    }

    // Run power method
    let mut lambda = 0.0;
    let mut iter: Int = 0;
    // C++ `while (iter++ < maxiter)`
    while iter < maxiter {
        iter += 1;
        func(&*v, &mut fv);
        let lambda_old = lambda;
        lambda = twonorm(&fv);
        for i in 0..dim as usize {
            v[i] = fv[i] / lambda;
        }
        if (lambda - lambda_old).abs() <= tol * lambda {
            break;
        }
    }
    lambda
}
