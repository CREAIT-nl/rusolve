//! Append-only CSC builder.
//!
//! [`CscMatrix`] is immutable; `SparseLp` load and a future `add_rows`
//! need growth with the **I2 lease**: appending rows or columns never
//! re-indexes existing rows or columns. This type pins that mutable-side
//! layout contract (only tests consume it today).
//!
//! Layout: one `(row, value)` staging `Vec` per column, rows ascending within
//! each column — the same column layout [`CscMatrix`] guarantees, so
//! [`CscBuilder::into_matrix`] is a flat copy, and [`CscBuilder::push_row`]
//! stays O(nnz of the row): a new row index is larger than every stored index,
//! so pushing it at each touched column's tail preserves ascending order.
//!
//! Conventions (mirroring the crate's existing formats):
//!
//! * Out-of-range indices are always asserted (like [`CscMatrix::from_triplets`]).
//! * Strictly-ascending (sorted, duplicate-free) entry lists are a caller
//!   contract, debug-asserted only — release builds trust the caller, matching
//!   the "sorted, dedup-checked (debug)" rule of [`crate::SparseLp`].
//! * Entries with `v == 0.0` (which catches `-0.0` too) are dropped at push
//!   time (invariant I3 of [`crate::SparseLp`]); [`CscBuilder::from_matrix`] is the
//!   one exception: absorbing an existing matrix copies its stored structure
//!   verbatim, zeros included.
//!
//! Provenance: layout/boundary level only, HiGHS `HighsSparseMatrix`
//! (MIT). No code copied.

use crate::CscMatrix;

/// Mutable, append-only builder producing [`CscMatrix`] values.
///
/// Column `j` keeps its index forever (invariant I1); row `r` keeps its index
/// forever and new rows land at the tail (invariant I2).
#[derive(Debug, Clone, PartialEq)]
pub struct CscBuilder {
    nrow: usize,
    /// Per-column staging `(row_idx, value)` lists, rows ascending.
    cols: Vec<Vec<(i32, f64)>>,
}

impl CscBuilder {
    /// Empty builder for matrices with `nrow` rows (rows still arrive via
    /// [`CscBuilder::push_row`]; `nrow` is the *current* row count).
    pub fn new(nrow: usize) -> Self {
        Self {
            nrow,
            cols: Vec::new(),
        }
    }

    /// Absorb an existing matrix for later extension (the `add_rows`
    /// entry point). Entries are copied verbatim — stored zeros included,
    /// since dropping them here would rewrite structure the caller gave us.
    pub fn from_matrix(a: &CscMatrix) -> Self {
        let ptr = a.col_ptr();
        let cols = (0..a.ncol())
            .map(|j| {
                a.row_idx()[ptr[j]..ptr[j + 1]]
                    .iter()
                    .copied()
                    .zip(a.values()[ptr[j]..ptr[j + 1]].iter().copied())
                    .collect()
            })
            .collect();
        Self {
            nrow: a.nrow(),
            cols,
        }
    }

    pub fn nrow(&self) -> usize {
        self.nrow
    }

    pub fn ncol(&self) -> usize {
        self.cols.len()
    }

    pub fn nnz(&self) -> usize {
        self.cols.iter().map(|c| c.len()).sum()
    }

    /// Append a column of `(row, value)` pairs, returning its index.
    ///
    /// Rows must be strictly ascending (debug-asserted) and `< nrow`
    /// (always asserted). Zero-valued entries are dropped (I3).
    ///
    /// # Panics
    ///
    /// Panics (release too) on an out-of-range row index; debug builds also
    /// panic on unsorted or duplicated rows.
    pub fn push_col(&mut self, entries: &[(usize, f64)]) -> usize {
        let mut col: Vec<(i32, f64)> = Vec::with_capacity(entries.len());
        let mut prev: Option<i32> = None;
        for &(r, v) in entries {
            assert!(
                r < self.nrow,
                "push_col: row {r} out of bounds (nrow {})",
                self.nrow
            );
            let r = r as i32;
            debug_assert!(
                prev.is_none_or(|p| p < r),
                "push_col: rows must be strictly ascending (got {r} after {p:?})",
                p = prev
            );
            prev = Some(r);
            if v == 0.0 {
                continue; // I3: zeros are structural, not stored
            }
            col.push((r, v));
        }
        let j = self.cols.len();
        self.cols.push(col);
        j
    }

    /// Append a row of `(col, value)` pairs at index `nrow`, growing the row
    /// count by one, returning the new row index.
    ///
    /// Columns must be strictly ascending (debug-asserted) and `< ncol`
    /// (always asserted). Zero-valued entries are dropped (I3). Existing
    /// columns keep their indices and entries; the new index exceeds every
    /// stored row index, so per-column ascending order is preserved by
    /// construction (I2).
    ///
    /// # Panics
    ///
    /// Panics (release too) on an out-of-range column index; debug builds
    /// also panic on unsorted or duplicated columns.
    pub fn push_row(&mut self, entries: &[(usize, f64)]) -> usize {
        let r = self.nrow;
        let mut prev: Option<usize> = None;
        for &(c, v) in entries {
            assert!(
                c < self.cols.len(),
                "push_row: column {c} out of bounds (ncol {})",
                self.cols.len()
            );
            debug_assert!(
                prev.is_none_or(|p| p < c),
                "push_row: columns must be strictly ascending (got {c} after {p:?})",
                p = prev
            );
            prev = Some(c);
            if v == 0.0 {
                continue; // I3
            }
            // r exceeds every stored index of cols[c] ⇒ tail push keeps it ascending.
            self.cols[c].push((r as i32, v));
        }
        self.nrow = r + 1;
        r
    }

    /// Finalize into an immutable [`CscMatrix`] (flat copy, no re-sort).
    pub fn into_matrix(self) -> CscMatrix {
        let nnz = self.nnz();
        let mut col_ptr = vec![0usize; self.cols.len() + 1];
        let mut row_idx = Vec::with_capacity(nnz);
        let mut val = Vec::with_capacity(nnz);
        for (j, col) in self.cols.iter().enumerate() {
            row_idx.extend(col.iter().map(|&(r, _)| r));
            val.extend(col.iter().map(|&(_, v)| v));
            col_ptr[j + 1] = row_idx.len();
        }
        CscMatrix {
            nrow: self.nrow,
            ncol: self.cols.len(),
            col_ptr,
            row_idx,
            val,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_col_builds_sorted_csc_layout() {
        let mut b = CscBuilder::new(3);
        assert_eq!(b.push_col(&[(0, 1.0), (2, 5.0)]), 0);
        assert_eq!(b.push_col(&[(1, 3.0)]), 1);
        assert_eq!(b.push_col(&[]), 2); // empty column is legal
        assert_eq!(b.push_col(&[(0, 2.0), (1, 4.0)]), 3);
        let a = b.into_matrix();
        assert_eq!((a.nrow(), a.ncol(), a.nnz()), (3, 4, 5));
        assert_eq!(a.col_ptr(), &[0, 2, 3, 3, 5]);
        assert_eq!(&a.row_idx()[a.col_ptr()[1]..a.col_ptr()[2]], &[1]);
        for (i, j, v) in [
            (0, 0, 1.0),
            (2, 0, 5.0),
            (1, 1, 3.0),
            (0, 2, 0.0),
            (0, 3, 2.0),
            (1, 3, 4.0),
            (2, 3, 0.0),
        ] {
            assert_eq!(a.get(i, j), v, "get({i}, {j})");
        }
    }

    #[test]
    fn push_col_and_push_row_drop_zeros() {
        let mut b = CscBuilder::new(2);
        b.push_col(&[(0, 1.0), (1, 0.0)]);
        b.push_col(&[(0, -0.0), (1, 2.0)]); // -0.0 == 0.0 ⇒ dropped
        b.push_row(&[(0, 0.0), (1, -0.0)]); // whole row lands zero-free
        assert_eq!(b.nnz(), 2);
        let a = b.into_matrix();
        assert_eq!((a.nrow(), a.ncol(), a.nnz()), (3, 2, 2));
        assert_eq!(a.get(0, 0), 1.0);
        assert_eq!(a.get(1, 1), 2.0);
        assert_eq!(a.get(1, 0), 0.0);
        assert_eq!(a.get(2, 0), 0.0);
        assert_eq!(a.get(2, 1), 0.0);
    }

    #[test]
    fn into_matrix_matches_from_triplets_bit_exact() {
        let triplets = [(0, 1, 2.5), (2, 0, -1.25), (1, 1, 7.0), (2, 2, 0.5)];
        let mut b = CscBuilder::new(3);
        b.push_col(&[(2, -1.25)]);
        b.push_col(&[(0, 2.5), (1, 7.0)]);
        b.push_col(&[(2, 0.5)]);
        assert_eq!(
            b.clone().into_matrix(),
            CscMatrix::from_triplets(3, 3, &triplets)
        );
    }

    /// The growth-lease test: absorbing a matrix and
    /// appending rows and columns leaves every pre-existing entry bit-identical
    /// at its pre-existing index (I1 + I2), with growth only at the tail.
    #[test]
    fn append_preserves_prefix() {
        let base = CscMatrix::from_triplets(
            3,
            4,
            &[
                (0, 0, 1.5),
                (2, 0, -3.0),
                (1, 1, 2.25),
                (0, 2, 0.125),
                (2, 3, 9.0),
            ],
        );
        let mut b = CscBuilder::from_matrix(&base);
        // Append one cut row referencing existing columns…
        assert_eq!(b.push_row(&[(0, -4.5), (3, 6.75)]), 3);
        // …and its new slack column (plus one more structural column).
        b.push_col(&[(3, 1.0)]);
        b.push_col(&[(0, -0.5), (3, 2.0)]);
        let grown = b.into_matrix();

        assert_eq!(grown.nrow(), 4);
        assert_eq!(grown.ncol(), 6);
        // I1/I2: every original entry unchanged, bit for bit.
        for j in 0..base.ncol() {
            for i in 0..base.nrow() {
                assert_eq!(
                    grown.get(i, j).to_bits(),
                    base.get(i, j).to_bits(),
                    "entry ({i}, {j}) moved or changed"
                );
            }
        }
        // New content landed where appended, nowhere else.
        assert_eq!(grown.get(3, 0), -4.5);
        assert_eq!(grown.get(3, 3), 6.75);
        assert_eq!(grown.get(3, 4), 1.0);
        assert_eq!(grown.get(0, 5), -0.5);
        assert_eq!(grown.get(2, 4), 0.0);
        // Absorbed stored zeros survive verbatim (from_matrix exception to I3).
        let with_zeros = CscMatrix::from_triplets(2, 1, &[(0, 0, 1.0), (0, 0, -1.0)]);
        assert_eq!(with_zeros.get(0, 0), 0.0);
        assert_eq!(with_zeros.nnz(), 1); // duplicates summed to a stored zero
        let kept = CscBuilder::from_matrix(&with_zeros).into_matrix();
        assert_eq!(kept, with_zeros);
        assert_eq!(kept.nnz(), 1);
    }

    #[test]
    fn row_only_then_column_growth_keeps_ascending_columns() {
        // Interleaved growth: rows and columns added alternately.
        let mut b = CscBuilder::new(0);
        for k in 0..4 {
            let row_entries: Vec<(usize, f64)> =
                (0..b.ncol()).map(|c| (c, c as f64 + k as f64)).collect();
            b.push_row(&row_entries);
            b.push_col(
                &(0..b.nrow())
                    .map(|r| (r, r as f64 - k as f64))
                    .collect::<Vec<_>>(),
            );
        }
        let a = b.into_matrix();
        assert_eq!((a.nrow(), a.ncol()), (4, 4));
        // Every column strictly ascending in row index (no debug-assert trips).
        for j in 0..a.ncol() {
            let rows = &a.row_idx()[a.col_ptr()[j]..a.col_ptr()[j + 1]];
            assert!(rows.windows(2).all(|w| w[0] < w[1]), "column {j} unsorted");
        }
        assert_eq!(a.get(2, 3), 2.0 - 3.0);
        assert_eq!(a.get(3, 3), 0.0); // (3,3) was pushed as zero ⇒ dropped
    }

    #[test]
    #[should_panic(expected = "out of bounds")]
    fn push_col_rejects_out_of_bounds_row() {
        let mut b = CscBuilder::new(2);
        b.push_col(&[(2, 1.0)]);
    }

    #[test]
    #[should_panic(expected = "out of bounds")]
    fn push_row_rejects_out_of_bounds_col() {
        let mut b = CscBuilder::new(2);
        b.push_col(&[(0, 1.0)]);
        b.push_row(&[(1, 1.0)]);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "strictly ascending")]
    fn push_col_rejects_unsorted() {
        let mut b = CscBuilder::new(4);
        b.push_col(&[(2, 1.0), (1, 2.0)]);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "strictly ascending")]
    fn push_col_rejects_duplicates() {
        let mut b = CscBuilder::new(4);
        b.push_col(&[(1, 1.0), (1, 2.0)]);
    }
}
