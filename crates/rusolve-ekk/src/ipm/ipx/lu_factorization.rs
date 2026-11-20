//! HiGHS source: ipm/ipx/lu_factorization.h, lu_factorization.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Abstract base class -> trait; base member `stability_` through
//! `stability_mut`. Implementor: `basiclu_kernel::BasicLuKernel`.

use super::sparse_matrix::{
    SparseMatrix, dot_column, infnorm_matrix, multiply_add, onenorm_matrix, scatter_column,
    triangular_solve,
};
use super::utils::{inverse_perm, onenorm};
use super::*;
use crate::util::highs_utils::cpp_max;

/// HiGHS: ipx/lu_factorization.h::LuFactorization (abstract)
pub trait LuFactorization {
    /// base member `stability_` (read)
    fn stability_ref(&self) -> f64;
    /// base member `stability_` (write)
    fn stability_mut(&mut self) -> &mut f64;

    /// HiGHS: ipx/lu_factorization.cc::LuFactorization::Factorize (lines 114-125)
    fn factorize(
        &mut self,
        dim: Int,
        Bbegin: &[Int],
        Bend: &[Int],
        Bi: &[Int],
        Bx: &[f64],
        pivottol: f64,
        strict_abs_pivottol: bool,
        L: &mut SparseMatrix,
        U: &mut SparseMatrix,
        rowperm: &mut Vec<Int>,
        colperm: &mut Vec<Int>,
        dependent_cols: &mut Vec<Int>,
    ) {
        self._factorize(
            dim,
            Bbegin,
            Bend,
            Bi,
            Bx,
            pivottol,
            strict_abs_pivottol,
            L,
            U,
            rowperm,
            colperm,
            dependent_cols,
        );
        let stability = stability_estimate(
            Bbegin,
            Bend,
            Bi,
            Bx,
            L,
            U,
            rowperm.as_slice(),
            colperm.as_slice(),
            dependent_cols.as_slice(),
        );
        *self.stability_mut() = stability;
    }

    /// HiGHS: ipx/lu_factorization.cc::LuFactorization::stability (lines 127-129)
    fn stability(&self) -> f64 {
        self.stability_ref()
    }

    /// HiGHS: ipx/lu_factorization.h::LuFactorization::_Factorize (pure virtual)
    #[allow(clippy::too_many_arguments)]
    fn _factorize(
        &mut self,
        dim: Int,
        Bbegin: &[Int],
        Bend: &[Int],
        Bi: &[Int],
        Bx: &[f64],
        pivottol: f64,
        strict_abs_pivottol: bool,
        L: &mut SparseMatrix,
        U: &mut SparseMatrix,
        rowperm: &mut Vec<Int>,
        colperm: &mut Vec<Int>,
        dependent_cols: &mut Vec<Int>,
    );
}

/// HiGHS: ipx/lu_factorization.cc::PermutedMatrix (static, lines 11-34)
fn permuted_matrix(
    Bbegin: &[Int],
    Bend: &[Int],
    Bi: &[Int],
    Bx: &[f64],
    rowperm: &[Int],
    colperm: &[Int],
    dependent_cols: &[Int],
) -> SparseMatrix {
    let dim = rowperm.len() as Int;
    let permuted_row = inverse_perm(rowperm);
    let mut dependent = vec![false; dim as usize];
    for &k in dependent_cols {
        dependent[k as usize] = true;
    }

    let mut b = SparseMatrix::with_dims(dim, 0);
    for k in 0..dim {
        if dependent[k as usize] {
            b.push_back(k, 1.0);
        } else {
            let j = colperm[k as usize];
            for p in Bbegin[j as usize]..Bend[j as usize] {
                b.push_back(permuted_row[Bi[p as usize] as usize], Bx[p as usize]);
            }
        }
        b.add_column();
    }
    b
}

/// HiGHS: ipx/lu_factorization.cc::SolveForward (static, lines 39-49)
fn solve_forward(L: &SparseMatrix, U: &SparseMatrix, rhs: &mut Vector, lhs: &mut Vector) {
    let dim = rhs.size();
    lhs.fill(0.0);
    for i in 0..dim {
        let iu = i as usize;
        rhs[iu] = if lhs[iu] >= 0.0 { 1.0 } else { -1.0 };
        lhs[iu] += rhs[iu];
        let alpha = -lhs[iu];
        scatter_column(L, i, alpha, lhs);
    }
    triangular_solve(U, lhs, b'n', "upper", 0);
}

/// HiGHS: ipx/lu_factorization.cc::SolveBackward (static, lines 54-67)
fn solve_backward(L: &SparseMatrix, U: &SparseMatrix, rhs: &mut Vector, lhs: &mut Vector) {
    let dim = rhs.size();
    lhs.fill(0.0);
    for j in 0..dim {
        let ju = j as usize;
        // C: `lhs[j] -= DotColumn(U, j, lhs);` (right operand sequenced first)
        let d = dot_column(U, j, lhs);
        lhs[ju] -= d;
        rhs[ju] = if lhs[ju] >= 0.0 { 1.0 } else { -1.0 };
        lhs[ju] += rhs[ju];
        let p = U.end(j) - 1;
        debug_assert_eq!(U.index(p), j);
        lhs[ju] /= U.value(p);
    }
    triangular_solve(L, lhs, b't', "lower", 1);
}

/// HiGHS: ipx/lu_factorization.cc::StabilityEstimate (static, lines 85-112)
#[allow(clippy::too_many_arguments)]
fn stability_estimate(
    Bbegin: &[Int],
    Bend: &[Int],
    Bi: &[Int],
    Bx: &[f64],
    L: &SparseMatrix,
    U: &SparseMatrix,
    rowperm: &[Int],
    colperm: &[Int],
    dependent_cols: &[Int],
) -> f64 {
    let dim = rowperm.len() as Int;
    let mut rhs = Vector::new(dim);
    let mut lhs = Vector::new(dim);

    // Compute 1-norm and infinity-norm of B.
    let b = permuted_matrix(Bbegin, Bend, Bi, Bx, rowperm, colperm, dependent_cols);
    let b_onenorm = onenorm_matrix(&b);
    let b_infnorm = infnorm_matrix(&b);

    solve_forward(L, U, &mut rhs, &mut lhs); // builds some rhs and computes lhs = B\rhs
    let norm_ftran = onenorm(&lhs);
    multiply_add(&b, &lhs, -1.0, &mut rhs, b'N'); // overwrites rhs by residual
    let norm_ftran_res = onenorm(&rhs);

    solve_backward(L, U, &mut rhs, &mut lhs); // builds some rhs and computes lhs = B'\rhs
    let norm_btran = onenorm(&lhs);
    multiply_add(&b, &lhs, -1.0, &mut rhs, b'T'); // overwrites rhs by residual
    let norm_btran_res = onenorm(&rhs);

    cpp_max(
        norm_ftran_res / (dim as f64 + b_onenorm * norm_ftran),
        norm_btran_res / (dim as f64 + b_infnorm * norm_btran),
    )
}
