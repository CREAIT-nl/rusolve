//! Core data structures for the rusolve solver.
//!
//! The sparse linear-algebra core is deliberately hand-rolled (no `sprs`
//! dependency): the solver will eventually want a bespoke memory layout for
//! factorization- and SIMD-friendliness, so we keep full control here.
//!
//! Two complementary sparse formats are provided:
//!
//! * [`CscMatrix`] — compressed sparse column (the natural layout for simplex
//!   column pricing / B-column extraction).
//! * [`CsrMatrix`] — compressed sparse row, with conversions both ways
//!   ([`CscMatrix::to_csr`], [`CsrMatrix::to_csc`]) and a cheap
//!   [`CscMatrix::transpose`].
//!
//! On top of the formats:
//!
//! * [`CscBuilder`] — append-only builder with the growth lease
//!   (appending rows/columns never re-indexes what exists).
//! * [`SparseLp`] / [`Lp`] — the caller-space and solver-space LP input
//!   plumbing: CSC matrices + row
//!   bounds / senses, with bit-identical conversion between the two forms.
//!
mod csc_builder;
pub mod diagnostic;
pub mod qp_convexity;
mod sparse_lp;

pub use csc_builder::CscBuilder;
pub use sparse_lp::{Lp, NO_ROW, SparseLp};

/// A coordinate (triplet) entry `(row, col, value)` of a sparse matrix.
pub type Triplet = (usize, usize, f64);

/// Compressed sparse column matrix: `A[i, j]` with `f64` values and `i32`
/// row indices. Column `j` occupies entries `col_ptr[j]..col_ptr[j + 1]`,
/// with row indices sorted ascending within each column.
#[derive(Debug, Clone, PartialEq)]
pub struct CscMatrix {
    nrow: usize,
    ncol: usize,
    col_ptr: Vec<usize>,
    row_idx: Vec<i32>,
    val: Vec<f64>,
}

impl CscMatrix {
    /// Builds a CSC matrix from unsorted triplets. Duplicate positions are
    /// combined by summing their values; entries with value `0.0` are kept
    /// (structure-preserving — presolve and factorization may care).
    ///
    /// # Panics
    ///
    /// Panics if any triplet index is out of bounds.
    pub fn from_triplets(nrow: usize, ncol: usize, triplets: &[Triplet]) -> Self {
        let mut sorted = triplets.to_vec();
        sorted.sort_unstable_by_key(|&(r, c, _)| (c, r));

        let mut col_ptr = vec![0usize; ncol + 1];
        let mut row_idx = Vec::with_capacity(sorted.len());
        let mut val = Vec::with_capacity(sorted.len());
        let mut prev: Option<(usize, usize)> = None;

        for (r, c, v) in sorted {
            assert!(r < nrow && c < ncol, "triplet index out of bounds");
            if prev == Some((r, c)) {
                // Duplicate position: adjacent after sorting, combine by summing.
                *val.last_mut().expect("non-empty val") += v;
            } else {
                col_ptr[c + 1] += 1;
                row_idx.push(r as i32);
                val.push(v);
                prev = Some((r, c));
            }
        }
        for c in 0..ncol {
            col_ptr[c + 1] += col_ptr[c];
        }

        Self {
            nrow,
            ncol,
            col_ptr,
            row_idx,
            val,
        }
    }

    pub fn nrow(&self) -> usize {
        self.nrow
    }

    pub fn ncol(&self) -> usize {
        self.ncol
    }

    pub fn nnz(&self) -> usize {
        self.val.len()
    }

    /// Column boundaries, length `ncol + 1`.
    pub fn col_ptr(&self) -> &[usize] {
        &self.col_ptr
    }

    /// Row index of each entry (sorted within each column).
    pub fn row_idx(&self) -> &[i32] {
        &self.row_idx
    }

    /// Value of each entry, parallel to [`CscMatrix::row_idx`].
    pub fn values(&self) -> &[f64] {
        &self.val
    }

    /// Mutable values, parallel to [`CscMatrix::row_idx`] (structure fixed).
    /// The one in-place mutation the pipeline needs: value scaling
    /// rewrites entries at their stored positions.
    pub fn values_mut(&mut self) -> &mut [f64] {
        &mut self.val
    }

    /// Dense access `A[i, j]`, returning `0.0` for structural zeros.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        assert!(i < self.nrow && j < self.ncol, "index out of bounds");
        let col = self.row_idx[self.col_ptr[j]..self.col_ptr[j + 1]].binary_search(&(i as i32));
        match col {
            Ok(k) => self.val[self.col_ptr[j] + k],
            Err(_) => 0.0,
        }
    }

    /// `y = A * x`.
    ///
    /// # Panics
    ///
    /// Panics if `x.len() != ncol` or `y.len() != nrow`.
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.ncol, "matvec: wrong x length");
        assert_eq!(y.len(), self.nrow, "matvec: wrong y length");
        y.fill(0.0);
        for (xc, rows) in x.iter().zip(self.col_ptr.windows(2)) {
            for k in rows[0]..rows[1] {
                y[self.row_idx[k] as usize] += self.val[k] * xc;
            }
        }
    }

    /// `y = A^T * x` without materializing the transpose.
    ///
    /// # Panics
    ///
    /// Panics if `x.len() != nrow` or `y.len() != ncol`.
    pub fn matvec_transpose(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.nrow, "matvec_t: wrong x length");
        assert_eq!(y.len(), self.ncol, "matvec_t: wrong y length");
        for (yc, rows) in y.iter_mut().zip(self.col_ptr.windows(2)) {
            let mut s = 0.0;
            for k in rows[0]..rows[1] {
                s += self.val[k] * x[self.row_idx[k] as usize];
            }
            *yc = s;
        }
    }

    /// Converts to row-major storage (counting sort over row indices).
    pub fn to_csr(&self) -> CsrMatrix {
        let mut row_ptr = vec![0usize; self.nrow + 1];
        for &r in &self.row_idx {
            row_ptr[r as usize + 1] += 1;
        }
        for r in 0..self.nrow {
            row_ptr[r + 1] += row_ptr[r];
        }

        let nnz = self.nnz();
        let mut col_idx = vec![0i32; nnz];
        let mut val = vec![0.0; nnz];
        let mut pos = row_ptr.clone();
        for c in 0..self.ncol {
            for k in self.col_ptr[c]..self.col_ptr[c + 1] {
                let r = self.row_idx[k] as usize;
                let dst = pos[r];
                col_idx[dst] = c as i32;
                val[dst] = self.val[k];
                pos[r] += 1;
            }
        }

        CsrMatrix {
            nrow: self.nrow,
            ncol: self.ncol,
            row_ptr,
            col_idx,
            val,
        }
    }

    /// Returns `A^T` in CSC form. Column indices of `A` become sorted row
    /// indices of `A^T` (and vice versa), so the transpose is exactly the
    /// CSR view reinterpreted as CSC.
    pub fn transpose(&self) -> CscMatrix {
        let csr = self.to_csr();
        CscMatrix {
            nrow: csr.ncol,
            ncol: csr.nrow,
            col_ptr: csr.row_ptr,
            row_idx: csr.col_idx,
            val: csr.val,
        }
    }
}

/// Compressed sparse row matrix: `A[i, j]` with `f64` values and `i32`
/// column indices, sorted ascending within each row.
#[derive(Debug, Clone, PartialEq)]
pub struct CsrMatrix {
    nrow: usize,
    ncol: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<i32>,
    val: Vec<f64>,
}

impl CsrMatrix {
    pub fn nrow(&self) -> usize {
        self.nrow
    }

    pub fn ncol(&self) -> usize {
        self.ncol
    }

    pub fn nnz(&self) -> usize {
        self.val.len()
    }

    /// Row boundaries, length `nrow + 1`.
    pub fn row_ptr(&self) -> &[usize] {
        &self.row_ptr
    }

    /// Column index of each entry (sorted within each row).
    pub fn col_idx(&self) -> &[i32] {
        &self.col_idx
    }

    /// Value of each entry, parallel to [`CsrMatrix::col_idx`].
    pub fn values(&self) -> &[f64] {
        &self.val
    }

    /// Dense access `A[i, j]`, returning `0.0` for structural zeros.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        assert!(i < self.nrow && j < self.ncol, "index out of bounds");
        let hit = self.col_idx[self.row_ptr[i]..self.row_ptr[i + 1]].binary_search(&(j as i32));
        match hit {
            Ok(k) => self.val[self.row_ptr[i] + k],
            Err(_) => 0.0,
        }
    }

    /// `y = A * x`.
    ///
    /// # Panics
    ///
    /// Panics if `x.len() != ncol` or `y.len() != nrow`.
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.ncol, "matvec: wrong x length");
        assert_eq!(y.len(), self.nrow, "matvec: wrong y length");
        for (yr, cols) in y.iter_mut().zip(self.row_ptr.windows(2)) {
            let mut s = 0.0;
            for k in cols[0]..cols[1] {
                s += self.val[k] * x[self.col_idx[k] as usize];
            }
            *yr = s;
        }
    }

    /// Converts to column-major storage (counting sort over column indices).
    pub fn to_csc(&self) -> CscMatrix {
        let mut col_ptr = vec![0usize; self.ncol + 1];
        for &c in &self.col_idx {
            col_ptr[c as usize + 1] += 1;
        }
        for c in 0..self.ncol {
            col_ptr[c + 1] += col_ptr[c];
        }

        let nnz = self.nnz();
        let mut row_idx = vec![0i32; nnz];
        let mut val = vec![0.0; nnz];
        let mut pos = col_ptr.clone();
        for r in 0..self.nrow {
            for k in self.row_ptr[r]..self.row_ptr[r + 1] {
                let c = self.col_idx[k] as usize;
                let dst = pos[c];
                row_idx[dst] = r as i32;
                val[dst] = self.val[k];
                pos[c] += 1;
            }
        }

        CscMatrix {
            nrow: self.nrow,
            ncol: self.ncol,
            col_ptr,
            row_idx,
            val,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A 3x4 matrix:
    //     [ 1  0  2  0 ]
    //     [ 0  3  0  4 ]
    //     [ 5  0  0  0 ]
    fn sample() -> CscMatrix {
        CscMatrix::from_triplets(
            3,
            4,
            &[
                (0, 0, 1.0),
                (2, 0, 5.0),
                (1, 1, 3.0),
                (0, 2, 2.0),
                (1, 3, 4.0),
            ],
        )
    }

    #[test]
    fn from_triplets_sorts_and_indexes_correctly() {
        let a = CscMatrix::from_triplets(2, 2, &[(1, 1, 7.0), (0, 1, 2.0), (1, 0, 9.0)]);
        assert_eq!(a.col_ptr(), &[0, 1, 3]);
        // Column 1 rows must be sorted ascending.
        assert_eq!(&a.row_idx()[a.col_ptr()[1]..a.col_ptr()[2]], &[0, 1]);
        assert_eq!(a.get(1, 1), 7.0);
        assert_eq!(a.get(0, 1), 2.0);
        assert_eq!(a.get(1, 0), 9.0);
        assert_eq!(a.get(0, 0), 0.0);
    }

    #[test]
    fn from_triplets_combines_duplicates_by_summing() {
        let a =
            CscMatrix::from_triplets(2, 2, &[(0, 0, 1.5), (0, 0, -0.5), (1, 1, 2.0), (0, 0, 1.0)]);
        assert_eq!(a.nnz(), 2);
        assert_eq!(a.get(0, 0), 2.0);
        assert_eq!(a.get(1, 1), 2.0);
    }

    #[test]
    #[should_panic]
    fn from_triplets_rejects_out_of_bounds() {
        CscMatrix::from_triplets(2, 2, &[(2, 0, 1.0)]);
    }

    #[test]
    fn matvec_matches_dense_expectation() {
        let a = sample();
        let x = [1.0, 2.0, 3.0, 4.0];
        let mut y = [0.0; 3];
        a.matvec(&x, &mut y);
        assert_eq!(y, [1.0 * 1.0 + 2.0 * 3.0, 3.0 * 2.0 + 4.0 * 4.0, 5.0 * 1.0]);

        let w = [1.0, 1.0, 1.0];
        let mut z = [0.0; 4];
        a.matvec_transpose(&w, &mut z);
        assert_eq!(z, [6.0, 3.0, 2.0, 4.0]);
    }

    #[test]
    fn csr_round_trip_preserves_matrix() {
        let a = sample();
        let csr = a.to_csr();
        assert_eq!(csr.nnz(), a.nnz());
        let back = csr.to_csc();
        assert_eq!(back, a);
    }

    #[test]
    fn transpose_matches_matvec_transpose() {
        let a = sample();
        let at = a.transpose();
        assert_eq!(
            (at.nrow(), at.ncol(), at.nnz()),
            (a.ncol(), a.nrow(), a.nnz())
        );

        let w = [1.5, -2.0, 0.25];
        let mut y = [0.0; 4];
        let mut z = [0.0; 4];
        at.matvec(&w, &mut y);
        a.matvec_transpose(&w, &mut z);
        assert_eq!(y, z);
        // (A^T)^T == A
        assert_eq!(at.transpose(), a);
    }
}
