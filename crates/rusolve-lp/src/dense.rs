//! Dense LU factorization with partial pivoting — the reference-grade
//! baseline that product-form updates are validated against.
//!
//! Provenance: textbook right-looking Gaussian elimination with row
//! pivoting (Golub & Van Loan, *Matrix Computations*, 4th ed., §3.1.4;
//! also described in HiGHS' `HFactor` refactorization path, MIT — read for
//! API shape, no code copied). The pivot-tolerance guard follows the
//! singular-detection convention of Forrest & Gold's BASP/TRAM
//! (ACM TOMS 18(3), 1992): a pivot below tolerance after max-magnitude
//! selection means the basis is (numerically) singular.

use crate::FactorError;

/// Square dense matrix, column-major (`data[i + j * n]`).
#[derive(Debug, Clone, PartialEq)]
pub struct DenseMatrix {
    n: usize,
    data: Vec<f64>,
}

impl DenseMatrix {
    /// All-zero `n x n` matrix.
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            data: vec![0.0; n * n],
        }
    }

    /// Identity matrix.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n);
        for i in 0..n {
            m.data[i + i * n] = 1.0;
        }
        m
    }

    /// Builds from basis columns: `columns[j]` becomes column `j`.
    ///
    /// # Panics
    ///
    /// Panics if any column length differs from the number of columns.
    pub fn from_columns(columns: &[Vec<f64>]) -> Self {
        let n = columns.len();
        for col in columns {
            assert_eq!(col.len(), n, "from_columns: column length != n");
        }
        let mut data = Vec::with_capacity(n * n);
        for col in columns {
            data.extend_from_slice(col);
        }
        Self { n, data }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i + j * self.n]
    }

    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i + j * self.n] = v;
    }

    /// Column `j` as a slice.
    pub fn col(&self, j: usize) -> &[f64] {
        &self.data[j * self.n..(j + 1) * self.n]
    }
}

/// LU factorization `P * A = L * U` of a square matrix with partial
/// (max-magnitude row) pivoting. `L` is unit lower triangular and stored
/// implicitly below the diagonal of the combined factor; `pivot[k]` records
/// the row index swapped into position `k` at elimination step `k`.
#[derive(Debug, Clone)]
pub struct DenseLu {
    n: usize,
    lu: Vec<f64>,
    pivot: Vec<u32>,
}

impl DenseLu {
    /// Factors `a` in place (the input is consumed). Returns
    /// [`FactorError::Singular`] if some column's largest remaining entry
    /// is below `pivot_tol` — never panics on rank deficiency.
    pub fn factor(a: &DenseMatrix, pivot_tol: f64) -> Result<Self, FactorError> {
        let n = a.n();
        let mut lu = a.data.clone();
        let mut pivot = vec![0u32; n];

        for k in 0..n {
            // Partial pivoting: largest |value| in column k, rows k..n.
            let mut best = k;
            let mut best_val = lu[k + k * n].abs();
            for i in k + 1..n {
                let v = lu[i + k * n].abs();
                if v > best_val {
                    best_val = v;
                    best = i;
                }
            }
            let piv = lu[best + k * n];
            if piv.abs() < pivot_tol {
                return Err(FactorError::Singular {
                    column: k,
                    pivot: piv,
                });
            }
            pivot[k] = best as u32;

            if best != k {
                // Swap entire rows (including already-eliminated L part,
                // which keeps the permutation representable as sequential
                // swaps of the working vector).
                for j in 0..n {
                    lu.swap(k + j * n, best + j * n);
                }
            }

            let inv = 1.0 / piv;
            for i in k + 1..n {
                let l_ik = lu[i + k * n] * inv;
                lu[i + k * n] = l_ik;
                for jj in k + 1..n {
                    let s = l_ik * lu[k + jj * n];
                    lu[i + jj * n] -= s;
                }
            }
        }

        Ok(Self { n, lu, pivot })
    }

    pub fn n(&self) -> usize {
        self.n
    }

    /// Solves `A x = b` via `P A = L U`: permute, forward-substitute L,
    /// back-substitute U.
    ///
    /// # Panics
    ///
    /// Panics if `b.len()` or `x.len()` differ from `n`.
    pub fn solve(&self, b: &[f64], x: &mut [f64]) {
        assert_eq!(b.len(), self.n, "solve: wrong rhs length");
        assert_eq!(x.len(), self.n, "solve: wrong solution length");
        x.copy_from_slice(b);

        // Apply row permutation P (same swap order used during factor).
        for k in 0..self.n {
            x.swap(k, self.pivot[k] as usize);
        }
        // Forward substitution with unit lower triangular L.
        for j in 0..self.n {
            let xj = x[j];
            if xj != 0.0 {
                for (i, xi) in x.iter_mut().enumerate().skip(j + 1) {
                    *xi -= self.lu[i + j * self.n] * xj;
                }
            }
        }
        // Back substitution with U.
        for j in (0..self.n).rev() {
            x[j] /= self.lu[j + j * self.n];
            let xj = x[j];
            for (i, xi) in x.iter_mut().enumerate().take(j) {
                *xi -= self.lu[i + j * self.n] * xj;
            }
        }
    }

    /// Solves `A^T x = b`: `U^T w = b`, `L^T v = w`, then undo permutation
    /// (`x = P^{-1} v`, i.e. replay the pivot swaps in reverse).
    ///
    /// # Panics
    ///
    /// Panics if `b.len()` or `x.len()` differ from `n`.
    pub fn solve_transpose(&self, b: &[f64], x: &mut [f64]) {
        assert_eq!(b.len(), self.n, "solve_t: wrong rhs length");
        assert_eq!(x.len(), self.n, "solve_t: wrong solution length");
        x.copy_from_slice(b);

        // Forward substitution with U^T.
        for j in 0..self.n {
            let mut s = x[j];
            for (i, &xi) in x.iter().take(j).enumerate() {
                s -= self.lu[i + j * self.n] * xi;
            }
            x[j] = s / self.lu[j + j * self.n];
        }
        // Back substitution with unit upper triangular L^T:
        // row i of L^T v = w reads v_i = w_i − Σ_{k>i} L[k][i]·v_k.
        for i in (0..self.n).rev() {
            let mut s = x[i];
            for (k, &xk) in x.iter().enumerate().take(self.n).skip(i + 1) {
                s -= self.lu[k + i * self.n] * xk;
            }
            x[i] = s;
        }
        // Undo permutation: swaps applied in reverse order.
        for k in (0..self.n).rev() {
            x.swap(k, self.pivot[k] as usize);
        }
    }

    /// Expands the combined factor into explicit `L` (unit diagonal) and `U`
    /// of `P A = L U`, for testing and diagnostics.
    pub fn expand(&self) -> (DenseMatrix, DenseMatrix) {
        let mut l = DenseMatrix::identity(self.n);
        let mut u = DenseMatrix::zeros(self.n);
        for j in 0..self.n {
            for i in 0..=j {
                u.set(i, j, self.lu[i + j * self.n]);
            }
            for i in j + 1..self.n {
                l.set(i, j, self.lu[i + j * self.n]);
            }
        }
        (l, u)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::XorShift64;
    use crate::test_support::{max_abs_diff, naive_solve, random_well_conditioned};

    #[test]
    fn lu_solve_matches_naive_gaussian_reference() {
        for seed in [1, 2, 3, 12345] {
            let a = random_well_conditioned(40, seed);
            let lu = DenseLu::factor(&a, 1e-12).expect("factorizable");
            let mut rng = XorShift64::new(seed ^ 0xABCD);
            let b: Vec<f64> = (0..40).map(|_| rng.uniform(-1.0, 1.0)).collect();

            let mut x = vec![0.0; 40];
            lu.solve(&b, &mut x);
            let ref_x = naive_solve(&a, &b);
            assert!(
                max_abs_diff(&x, &ref_x) < 1e-9,
                "seed {seed}: LU solve diverged from reference"
            );
        }
    }

    #[test]
    fn lu_solve_transpose_matches_naive_reference() {
        for seed in [7, 8] {
            let a = random_well_conditioned(35, seed);
            let lu = DenseLu::factor(&a, 1e-12).expect("factorizable");
            let mut rng = XorShift64::new(seed ^ 0x1234);
            let b: Vec<f64> = (0..35).map(|_| rng.uniform(-1.0, 1.0)).collect();

            let mut x = vec![0.0; 35];
            lu.solve_transpose(&b, &mut x);

            // Reference: solve A^T x = b via the naive solver on A^T.
            let n = a.n();
            let mut at = DenseMatrix::zeros(n);
            for j in 0..n {
                for i in 0..n {
                    at.set(i, j, a.get(j, i));
                }
            }
            let ref_x = naive_solve(&at, &b);
            assert!(
                max_abs_diff(&x, &ref_x) < 1e-9,
                "seed {seed}: transpose solve diverged from reference"
            );
        }
    }

    #[test]
    fn lu_expand_satisfies_pa_equals_lu() {
        let a = random_well_conditioned(24, 99);
        let lu = DenseLu::factor(&a, 1e-12).expect("factorizable");
        let (l, u) = lu.expand();

        // Rebuild P by replaying pivot swaps on identity rows.
        let n = a.n();
        let mut pa = a.clone();
        for k in 0..n {
            let p = lu.pivot[k] as usize;
            if p != k {
                for j in 0..n {
                    pa.data.swap(k + j * n, p + j * n);
                }
            }
        }
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for r in 0..n {
                    s += l.get(i, r) * u.get(r, j);
                }
                assert!((s - pa.get(i, j)).abs() < 1e-10, "L*U != P*A at {i},{j}");
            }
        }
    }

    #[test]
    fn singular_matrix_returns_error_not_panic() {
        // Rank-deficient: column 2 = column 0 + column 1.
        let mut a = random_well_conditioned(12, 5);
        for i in 0..12 {
            a.set(i, 2, a.get(i, 0) + a.get(i, 1));
        }
        let err = DenseLu::factor(&a, 1e-9).expect_err("must not factor");
        assert!(matches!(err, FactorError::Singular { .. }));
    }

    #[test]
    fn near_zero_pivot_below_tolerance_is_singular() {
        let mut a = DenseMatrix::identity(4);
        a.set(3, 3, 1e-14);
        assert!(DenseLu::factor(&a, 1e-12).is_err());
        // Same matrix above tolerance factors fine.
        assert!(DenseLu::factor(&a, 1e-16).is_ok());
    }
}
