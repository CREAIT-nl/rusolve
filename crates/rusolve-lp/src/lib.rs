//! Linear-programming basis numerics: the factorization layer that the dual
//! simplex engine will sit on.
//!
//! Contents:
//!
//! * [`dense::DenseLu`] — dense LU with partial pivoting: the
//!   reference-grade correctness baseline and the refactorization engine for
//!   small bases.
//! * [`pfi::ProductFormFactorization`] — product-form basis factorization:
//!   LU of the initial basis plus a file of eta vectors recording successive
//!   single-column basis updates, with `solve_b` / `solve_bt` and
//!   threshold-triggered refactorization.
//!
//! Provenance (every algorithmic idea cites a source):
//!
//! * Product-form / eta techniques: Forrest & Gold, "BASP — the Bell ABS
//!   package: revised simplex programs", ACM Trans. Math. Softw. 18(3):349–375
//!   (1992), which defines the eta-file layout, update, and refactorize-on-age
//!   policy; and F. Huangfu & J. A. J. Hall, "Parallelizing the dual revised
//!   simplex method", arXiv:1503.01889 (Math. Prog. Comp. 2018), whose
//!   published description of the product-form solve/update is followed
//!   (equations in the `pfi` module docs). Module boundaries follow HiGHS'
//!   `HFactor` (MIT).
//! * Dense LU: Golub & Van Loan, *Matrix Computations*, 4th ed., §3.1.
//! * Pivot-tolerance singular detection: BASP convention (above).
//!
//! No `unsafe`. f64 throughout; sparse indices are `i32` to match
//! `rusolve-core`.

// Numeric kernels index several parallel arrays per loop; clippy 1.99 now flags
// those (rusolve-ekk allows the same lint).
#![allow(clippy::needless_range_loop)]
#![allow(clippy::collapsible_match)]

/// The LP input type IS the core CSC `Lp` — the crate-local dense
/// duplicate is gone.
pub use rusolve_core::Lp;
/// Sparse LP front door: the type
/// `SimplexSolver::solve_sparse` consumes. Zero-cost re-export of the
/// `rusolve-core` definition (the conversions live there to avoid a cycle).
pub use rusolve_core::SparseLp;
pub use simplex::{
    Basis, BasisStatus, SimplexConfig, SimplexError, SimplexSolution, SimplexSolver, SimplexStatus,
};

pub mod dense;
pub mod dual;
mod hvec;
pub mod ipm;
pub mod lu;
// FT update scaffolding (see module doc).
// TODO: thread FtConfig into the factorization call sites.
pub mod lu_ft;
pub mod pfi;
pub mod presolve;
pub mod rng;
pub mod scaling;
pub mod simplex;

pub use dense::{DenseLu, DenseMatrix};
pub use dual::{
    CutRow, DseMode, DseWeights, ReoptResult, ReoptSession, ReoptStatus, StdCache, StdLayout,
    TableauBasic, TableauRow, reoptimize, reoptimize_dse, reoptimize_dse_cached,
    reoptimize_dse_cached_until,
};
pub use pfi::ProductFormFactorization;

/// Failure modes of basis factorization and update.
///
/// The distinction matters to the (future) simplex loop: [`FactorError::Singular`]
/// from an initial factorization means the basis itself is unusable, while
/// [`FactorError::RefactorizeNeeded`] from an update means accumulated eta-file
/// error has made a pivot unreliable and the standard recovery is to rebuild
/// the factorization from the explicit basis columns before retrying.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FactorError {
    /// A pivot of magnitude below the tolerance was found after max-magnitude
    /// row selection: the matrix is (numerically) singular.
    Singular { column: usize, pivot: f64 },
    /// An update pivot fell below tolerance against the *updated* factorization;
    /// results computed through the eta file can no longer be trusted.
    RefactorizeNeeded { pivot: f64 },
}

impl std::fmt::Display for FactorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FactorError::Singular { column, pivot } => {
                write!(
                    f,
                    "singular basis: pivot {pivot:e} below tol at column {column}"
                )
            }
            FactorError::RefactorizeNeeded { pivot } => {
                write!(
                    f,
                    "eta-file pivot {pivot:e} below tol: refactorize required"
                )
            }
        }
    }
}

impl std::error::Error for FactorError {}

/// Tunables for [`ProductFormFactorization`].
#[derive(Debug, Clone, Copy)]
pub struct PfiConfig {
    /// Pivots below this magnitude (after max-magnitude row selection) are
    /// treated as singularity / refactorize signals rather than divided by.
    pub pivot_tol: f64,
    /// Maximum eta-file length before [`ProductFormFactorization::update`]
    /// auto-refactorizes. Mirrors BASP's refactorization interval.
    pub eta_limit: usize,
    /// Minimum updates before the cost-based refactorization trigger may
    /// fire (sparse engine; see `ProductFormFactorization::update`).
    pub min_updates: usize,
    /// Refactor-clock scale `K`: the cost trigger fires once the eta solve
    /// ticks reach `K × build_work` (sparse engine). `1.0` is the
    /// calibrated PFI clock; the dual simplex runs 4.0 (measured on our
    /// benchmark set). `RUSOLVE_RFCAD`
    /// overrides both for experiments.
    pub clock_scale: f64,
    /// Basis factorization engine. `Sparse` swaps the dense LU for
    /// [`crate::lu::SparseLu`] (sparse solves + 10× cheaper refactorization
    /// on basis-like matrices). Netlib A/B: geo-mean time ratio 50.8×→31.3×,
    /// 27/28 solved, tuff no longer times out ⇒ DEFAULT. A sparse factorization that hits Singular falls
    /// back to dense STICKILY (nasty bases keep the safety net).
    pub engine: PfiEngine,
}

/// Factorization engine selector ([`PfiConfig::engine`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PfiEngine {
    Dense,
    Sparse,
}

impl Default for PfiConfig {
    fn default() -> Self {
        Self {
            pivot_tol: 1e-10,
            eta_limit: 1000,
            min_updates: 2,
            clock_scale: 1.0,
            engine: PfiEngine::Sparse,
        }
    }
}

#[cfg(test)]
mod test_support {
    use crate::dense::DenseMatrix;
    use crate::rng::XorShift64;

    /// Well-conditioned random `n x n`: diagonally dominant by construction
    /// (diag in [1,2], dense perturbation in [-0.1, 0.1]).
    pub(crate) fn random_well_conditioned(n: usize, seed: u64) -> DenseMatrix {
        let mut rng = XorShift64::new(seed);
        let mut a = DenseMatrix::zeros(n);
        for j in 0..n {
            for i in 0..n {
                a.set(i, j, rng.uniform(-0.1, 0.1));
            }
            a.set(j, j, a.get(j, j) + rng.uniform(1.0, 2.0));
        }
        a
    }

    /// Naive dense Gaussian elimination with partial pivoting on
    /// `Vec<Vec<f64>>` (row-major) — an *independent* reference for the
    /// column-major in-place LU and the product-form solves. Not fast, not
    /// clever: that is the point.
    pub(crate) fn naive_solve(a: &DenseMatrix, b: &[f64]) -> Vec<f64> {
        let n = a.n();
        let mut m: Vec<Vec<f64>> = Vec::with_capacity(n);
        for (i, &bi) in b.iter().enumerate() {
            let mut row: Vec<f64> = (0..n).map(|j| a.get(i, j)).collect();
            row.push(bi);
            m.push(row);
        }
        for k in 0..n {
            let mut p = k;
            for i in k + 1..n {
                if m[i][k].abs() > m[p][k].abs() {
                    p = i;
                }
            }
            m.swap(k, p);
            let piv = m[k][k];
            assert!(piv.abs() > 1e-300, "naive reference: singular");
            let tail: Vec<f64> = m[k][k..].to_vec();
            for row in m.iter_mut().skip(k + 1) {
                let f = row[k] / piv;
                for (off, &mk) in tail.iter().enumerate() {
                    row[k + off] -= f * mk;
                }
            }
        }
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let s: f64 = m[i][i + 1..n]
                .iter()
                .zip(&x[i + 1..n])
                .map(|(&a, &x)| a * x)
                .sum();
            x[i] = (m[i][n] - s) / m[i][i];
        }
        x
    }

    /// `||a - b||_inf`, shared by the numerics tests.
    pub(crate) fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0_f64, f64::max)
    }
}
