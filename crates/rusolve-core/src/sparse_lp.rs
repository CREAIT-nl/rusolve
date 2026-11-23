//! Sparse LP input plumbing. [`SparseLp`] is the caller-space LP: CSC constraint matrix +
//! row bounds (`row_lower ≤ Ax ≤ row_upper`) + column bounds + costs + sense,
//! the exact `LpEngine::load` type of the native MIP search.
//!
//! [`Lp`] is the solver-space row form (one sense per row, `rhs`, column
//! bounds): its matrix is the SAME [`CscMatrix`] storage — the
//! dense `col: Vec<Vec<f64>>` layout is gone from the pipeline (it cost
//! n·m f64 per copy; krka-class instances measured 12 GB peak RSS through
//! it). `to_lp` materializes the ranged-row
//! expansion into CSC with replicated entries; `from_lp` is a bounds
//! mapping plus a matrix copy.
//!
//! The acceptance spine is bit-identity: every conversion here preserves the
//! fp operation sequence — in fact it performs *no* fp arithmetic at all
//! (entry copies and sense/bound routing), and the ranged-row expansion
//! reproduces the historical `RusolveEngine::new` row order (equality rows
//! once; ranged rows hi-L first, then lo-G) that the dense differential
//! oracle's iteration counts are load-bearing on.
//!
//! Invariants (asserted in [`SparseLp::new`], the only
//! constructor used across the MIP seam):
//!
//! * **I1** — column `j` of `a` IS caller variable `j`; nothing here ever
//!   reorders columns, and growth (via [`crate::CscBuilder`]) appends only.
//! * **I2** — ranged-row expansion is a function of row order only and is
//!   append-only at the tail; [`SparseLp::row_map`] exposes the expansion map.
//! * **I3** — stored zeros are dropped once at load (debug-asserting sorted,
//!   duplicate-free columns), so consumers may treat every stored entry as
//!   nonzero.
//!
//! Provenance: HiGHS `HighsLp` (MIT) —
//! "CSC matrix + row bounds is the LP's public shape"; layout independently
//! expressed in Rust, no code copied.

use crate::CscBuilder;
use crate::CscMatrix;

/// The solver-space LP in row form: one sense per row (`rhs` + `senses`),
/// CSC matrix, column bounds and costs. [`SparseLp`] is the caller-space
/// ranged form; [`SparseLp::to_lp`] / [`SparseLp::from_lp`] convert.
#[derive(Debug, Clone, PartialEq)]
pub struct Lp {
    pub num_col: usize,
    pub num_row: usize,
    /// Constraint matrix, `num_row × num_col`, CSC. Stored zeros are legal
    /// here (the io reader keeps them); every consumer skips `v == 0.0`,
    /// matching the dense-era scans that dropped zeros on the floor.
    pub a: CscMatrix,
    pub rhs: Vec<f64>,
    pub objective: Vec<f64>,
    /// true ⇒ maximize (internally minimizes −cᵀx).
    pub maximize: bool,
    /// Row senses: 'L' (≤), 'G' (≥), 'E' (=).
    pub senses: Vec<char>,
    /// Column lower bounds (per variable, not per row).
    pub lower: Vec<f64>,
    /// Column upper bounds.
    pub upper: Vec<f64>,
}

/// Sentinel in [`SparseLp::row_map`] slots: this slot produced no dense row
/// (e.g. the hi slot of a `≥`-only row).
pub const NO_ROW: usize = usize::MAX;

/// Caller-space LP in CSC form: `min cᵀx` (or max, see `maximize`) s.t.
/// `row_lower ≤ Ax ≤ row_upper`, `col_lower ≤ x ≤ col_upper`.
///
/// Bounds may be ±inf (`lo == hi` ⇒ equality row; `lo == -inf, hi == +inf`
/// ⇒ vacuous row, expanded away). Semantics beyond shape (empty bounds,
/// NaN-ridden data) remain the solver's `to_standard` checks — this type
/// asserts dimensions and `lo ≤ hi` only.
#[derive(Debug, Clone, PartialEq)]
pub struct SparseLp {
    /// Constraint matrix, `ncol == col_lower.len()`, `nrow == row_lower.len()`.
    pub a: CscMatrix,
    /// Row lower bounds, `-inf` allowed, length `nrow`.
    pub row_lower: Vec<f64>,
    /// Row upper bounds, `+inf` allowed; `lo == hi` ⇒ equality row.
    pub row_upper: Vec<f64>,
    /// Column lower bounds, length `ncol`.
    pub col_lower: Vec<f64>,
    /// Column upper bounds.
    pub col_upper: Vec<f64>,
    /// Objective coefficients, length `ncol`.
    pub objective: Vec<f64>,
    /// true ⇒ maximize.
    pub maximize: bool,
}

impl SparseLp {
    /// The only constructor across the MIP seam: asserts shape invariants
    /// (dimensions, `lo ≤ hi` — which also rejects NaN bounds) and applies
    /// I3 (stored zeros dropped once at load; columns debug-asserted sorted
    /// and duplicate-free).
    ///
    /// # Panics
    ///
    /// Panics on inconsistent dimensions or inverted/NaN bounds.
    pub fn new(
        a: CscMatrix,
        row_lower: Vec<f64>,
        row_upper: Vec<f64>,
        col_lower: Vec<f64>,
        col_upper: Vec<f64>,
        objective: Vec<f64>,
        maximize: bool,
    ) -> Self {
        let (m, n) = (a.nrow(), a.ncol());
        assert_eq!(row_lower.len(), m, "row_lower length != nrow");
        assert_eq!(row_upper.len(), m, "row_upper length != nrow");
        assert_eq!(col_lower.len(), n, "col_lower length != ncol");
        assert_eq!(col_upper.len(), n, "col_upper length != ncol");
        assert_eq!(objective.len(), n, "objective length != ncol");
        for (i, (&lo, &hi)) in row_lower.iter().zip(&row_upper).enumerate() {
            assert!(lo <= hi, "row {i}: row_lower {lo} > row_upper {hi}");
        }
        for (j, (&lo, &hi)) in col_lower.iter().zip(&col_upper).enumerate() {
            assert!(lo <= hi, "column {j}: col_lower {lo} > col_upper {hi}");
        }
        Self {
            a: Self::drop_stored_zeros(a),
            row_lower,
            row_upper,
            col_lower,
            col_upper,
            objective,
            maximize,
        }
    }

    /// I3: rebuild the matrix without `0.0` entries (catches `-0.0` too),
    /// debug-checking the column layout via [`CscBuilder`] when a rebuild is
    /// needed. A zero-free matrix is passed through with only a debug check —
    /// load is once-per-solve, but copying nnz needlessly is still waste.
    fn drop_stored_zeros(a: CscMatrix) -> CscMatrix {
        let ptr = a.col_ptr();
        if !a.values().contains(&0.0) {
            debug_assert!(
                (0..a.ncol()).all(|j| a.row_idx()[ptr[j]..ptr[j + 1]]
                    .windows(2)
                    .all(|w| w[0] < w[1])),
                "SparseLp::new: columns must be strictly ascending"
            );
            return a;
        }
        let mut b = CscBuilder::new(a.nrow());
        for j in 0..a.ncol() {
            let col: Vec<(usize, f64)> = a.row_idx()[ptr[j]..ptr[j + 1]]
                .iter()
                .zip(a.values()[ptr[j]..ptr[j + 1]].iter())
                .map(|(&r, &v)| (r as usize, v))
                .collect();
            b.push_col(&col);
        }
        b.into_matrix()
    }

    /// Load from the row-form [`Lp`]: senses/rhs route to row bounds
    /// (`'L' → (-inf, rhs]`, `'G' → [rhs, +inf)`, `'E' → [rhs, rhs]`),
    /// column bounds/costs are copied verbatim, and the matrix is carried
    /// over as-is (same row space) with stored zeros dropped (I3) — explicit
    /// `0.0` entries are not stored, and `-0.0` entries normalize to absent
    /// (no consumer's fp sequence distinguishes them — design §3.1).
    ///
    /// # Panics
    ///
    /// Panics on inconsistent dimensions or a row sense outside `L/G/E`
    /// (plus everything [`SparseLp::new`] rejects).
    pub fn from_lp(lp: &Lp) -> Self {
        let (n, m) = (lp.num_col, lp.num_row);
        assert_eq!(lp.a.ncol(), n, "matrix ncol != num_col");
        assert_eq!(lp.a.nrow(), m, "matrix nrow != num_row");
        assert_eq!(lp.rhs.len(), m, "rhs length != num_row");
        assert_eq!(lp.senses.len(), m, "senses length != num_row");
        assert_eq!(lp.lower.len(), n, "lower length != num_col");
        assert_eq!(lp.upper.len(), n, "upper length != num_col");
        assert_eq!(lp.objective.len(), n, "objective length != num_col");

        let mut row_lower = Vec::with_capacity(m);
        let mut row_upper = Vec::with_capacity(m);
        for (i, &s) in lp.senses.iter().enumerate() {
            let (lo, hi) = match s {
                'L' => (f64::NEG_INFINITY, lp.rhs[i]),
                'G' => (lp.rhs[i], f64::INFINITY),
                'E' => (lp.rhs[i], lp.rhs[i]),
                other => panic!("row {i} sense '{other}'"),
            };
            row_lower.push(lo);
            row_upper.push(hi);
        }
        Self::new(
            lp.a.clone(),
            row_lower,
            row_upper,
            lp.lower.clone(),
            lp.upper.clone(),
            lp.objective.clone(),
            lp.maximize,
        )
    }

    /// The ranged-row expansion map (I2): slot layout
    /// `[hi_slot, lo_slot]` —
    ///
    /// * equality row (`lo == hi`): `[idx, NO_ROW]`, the dense row is `'E'`;
    /// * `hi` finite: `[idx, _]` gets the `'L'` row `aᵀx ≤ hi` **first**;
    /// * `lo` finite and `lo != hi`: `[_ , idx]` gets the `'G'` row `aᵀx ≥ lo`
    ///   **second**;
    /// * unbounded row (`-inf ≤ row ≤ +inf`): `[NO_ROW, NO_ROW]`, expanded away.
    ///
    /// Dense row indices rise with caller row order, so appending caller rows
    /// appends dense rows at the tail — the property `add_rows` leases.
    pub fn row_map(&self) -> Vec<[usize; 2]> {
        self.expansion().0
    }

    /// `(row_map, senses, rhs)` of the dense expansion, one ascending pass.
    fn expansion(&self) -> (Vec<[usize; 2]>, Vec<char>, Vec<f64>) {
        let mut map = vec![[NO_ROW, NO_ROW]; self.a.nrow()];
        let mut senses = Vec::new();
        let mut rhs = Vec::new();
        for (slot, (&lo, &hi)) in map
            .iter_mut()
            .zip(self.row_lower.iter().zip(&self.row_upper))
        {
            if lo == hi {
                slot[0] = senses.len();
                senses.push('E');
                rhs.push(lo);
            } else {
                if hi.is_finite() {
                    slot[0] = senses.len();
                    senses.push('L');
                    rhs.push(hi);
                }
                if lo.is_finite() {
                    slot[1] = senses.len();
                    senses.push('G');
                    rhs.push(lo);
                }
            }
        }
        (map, senses, rhs)
    }

    /// Materialize the row-form [`Lp`], bit-identical to the MIP engine's
    /// historical template build (`rusolve_mip::engine::RusolveEngine::new`):
    /// same row order (equality once; ranged hi-L then lo-G — I2), and CSC
    /// entries produced by plain copies of the CSC values replicated into
    /// both mapped rows (no fp arithmetic).
    pub fn to_lp(&self) -> Lp {
        let (map, senses, rhs) = self.expansion();
        let num_row = senses.len();
        let n = self.a.ncol();
        let (ptr, idx, val) = (self.a.col_ptr(), self.a.row_idx(), self.a.values());
        let mut b = CscBuilder::new(num_row);
        let mut col: Vec<(usize, f64)> = Vec::new();
        for j in 0..n {
            col.clear();
            for k in ptr[j]..ptr[j + 1] {
                let (r, v) = (idx[k] as usize, val[k]);
                for &t in map[r].iter() {
                    if t != NO_ROW {
                        col.push((t, v));
                    }
                }
            }
            // map[r] = [hi, lo] is ascending across and within rows, so the
            // replicated entries land row-ascending like the dense scatter.
            b.push_col(&col);
        }
        Lp {
            num_col: n,
            num_row,
            a: b.into_matrix(),
            rhs,
            objective: self.objective.clone(),
            maximize: self.maximize,
            senses,
            lower: self.col_lower.clone(),
            upper: self.col_upper.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test-only dense materialization (dense columns are allowed in tests,
    /// design §6): one `Vec<f64>` per column, length `nrow`, stored zeros
    /// landing at their positions.
    fn dense_cols(a: &CscMatrix) -> Vec<Vec<f64>> {
        let (m, n) = (a.nrow(), a.ncol());
        let mut cols = vec![vec![0.0; m]; n];
        let ptr = a.col_ptr();
        for j in 0..n {
            for k in ptr[j]..ptr[j + 1] {
                cols[j][a.row_idx()[k] as usize] = a.values()[k];
            }
        }
        cols
    }

    /// Test-only constructor standing in for the dense-era `Lp { col, .. }`
    /// literal: dense columns (zeros included, like the io reader's stored
    /// zeros) become a CSC matrix with every position stored, bit-copied.
    fn lp_from_dense_cols(
        col: &[Vec<f64>],
        rhs: Vec<f64>,
        objective: Vec<f64>,
        maximize: bool,
        senses: Vec<char>,
        lower: Vec<f64>,
        upper: Vec<f64>,
    ) -> Lp {
        let n = col.len();
        let m = col.first().map_or(0, |c| c.len());
        let triplets: Vec<(usize, usize, f64)> = col
            .iter()
            .enumerate()
            .flat_map(|(j, c)| c.iter().enumerate().map(move |(i, &v)| (i, j, v)))
            .collect();
        Lp {
            num_col: n,
            num_row: m,
            a: CscMatrix::from_triplets(m, n, &triplets),
            rhs,
            objective,
            maximize,
            senses,
            lower,
            upper,
        }
    }

    /// Overwrite the stored entry at `(i, j)` (must exist — crafted matrices
    /// store every position). Test seam for `-0.0` normalization.
    fn set_stored(a: &mut CscMatrix, i: usize, j: usize, v: f64) {
        let ptr = a.col_ptr().to_vec();
        let mut found = false;
        for k in ptr[j]..ptr[j + 1] {
            if a.row_idx()[k] as usize == i {
                a.values_mut()[k] = v;
                found = true;
            }
        }
        assert!(found, "no stored entry at ({i}, {j})");
    }

    /// Mixed-bound, mixed-sense LP: L/G/E rows, negative rhs, a flipped
    /// variable (-inf, ub], a free variable, a fixed variable, and structural
    /// zeros — the crafted shape the dense→sparse→dense bit-identity spine
    /// gates on.
    fn crafted_lp() -> Lp {
        lp_from_dense_cols(
            &[
                vec![1.0, 0.0, -2.5, 4.0, 0.0],
                vec![0.5, 3.0, 0.0, 0.0, -1.25],
                vec![-0.75, 0.0, 6.0, 2.0, 0.0],
                vec![0.0, 0.0, 0.0, -8.5, 9.0],
            ],
            vec![7.0, -3.0, 0.0, 12.5, -0.5],
            vec![2.0, -1.0, 0.0, 3.5],
            true,
            vec!['L', 'G', 'E', 'L', 'G'],
            vec![0.0, f64::NEG_INFINITY, -2.0, 5.0],
            vec![10.0, 4.0, f64::INFINITY, 5.0],
        )
    }

    fn bits_eq(a: &[Vec<f64>], b: &[Vec<f64>]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| bits_eq_flat(x, y))
    }

    fn bits_eq_flat(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(p, q)| p.to_bits() == q.to_bits())
    }

    #[test]
    fn round_trip_dense_sparse_dense_bit_identical() {
        let lp = crafted_lp();
        let back = SparseLp::from_lp(&lp).to_lp();
        assert_eq!(back.num_col, lp.num_col);
        assert_eq!(back.num_row, lp.num_row);
        assert_eq!(back.senses, lp.senses);
        assert_eq!(back.maximize, lp.maximize);
        assert!(
            bits_eq(&dense_cols(&back.a), &dense_cols(&lp.a)),
            "dense entries not bit-identical"
        );
        assert!(bits_eq_flat(&back.rhs, &lp.rhs));
        assert!(bits_eq_flat(&back.objective, &lp.objective));
        assert!(bits_eq_flat(&back.lower, &lp.lower));
        assert!(bits_eq_flat(&back.upper, &lp.upper));
        // And the sparse form itself is canonical: I3 dropped every zero.
        let sp = SparseLp::from_lp(&lp);
        assert_eq!(
            sp.a.nnz(),
            dense_cols(&lp.a)
                .iter()
                .map(|c| c.iter().filter(|&&v| v != 0.0).count())
                .sum::<usize>()
        );
        assert!(sp.a.values().iter().all(|&v| v != 0.0));
    }

    #[test]
    fn from_lp_normalizes_negative_zero_entries() {
        // Documented I3 edge: -0.0 entries are dropped like zeros; the dense
        // rebuild carries +0.0 there. Value-equal, bit-normalized — no
        // consumer's fp sequence distinguishes them (design §3.1).
        let mut lp = crafted_lp();
        set_stored(&mut lp.a, 1, 0, -0.0); // was a structural 0.0 ⇒ stays 0 entries
        set_stored(&mut lp.a, 0, 1, -0.0); // was 0.5 ⇒ the value is dropped as a zero
        let sp = SparseLp::from_lp(&lp);
        assert_eq!(sp.a.get(1, 0), 0.0);
        assert_eq!(sp.a.get(0, 1), 0.0);
        assert!(sp.a.values().iter().all(|&v| v != 0.0));
        let back = sp.to_lp();
        let cols = dense_cols(&back.a);
        assert_eq!(cols[1][0].to_bits(), 0.0f64.to_bits()); // +0.0, not -0.0
        assert_eq!(cols[0][1], 0.0);
    }

    // Name mirrors the design's "hi-L then lo-G" vocabulary (grep target).
    #[allow(non_snake_case)]
    #[test]
    fn ranged_row_expands_hi_L_then_lo_G() {
        // I2 order, matching RusolveEngine::new: hi-L first, lo-G second.
        let mut b = CscBuilder::new(2);
        b.push_col(&[(0, 1.0), (1, -2.0)]);
        b.push_col(&[(0, 3.0), (1, 0.5)]);
        let sp = SparseLp::new(
            b.into_matrix(),
            vec![2.0, -3.0],
            vec![7.0, f64::INFINITY],
            vec![0.0, 0.0],
            vec![f64::INFINITY, f64::INFINITY],
            vec![1.0, 1.0],
            false,
        );
        let dense = sp.to_lp();
        assert_eq!(dense.senses, vec!['L', 'G', 'G']);
        assert_eq!(dense.rhs, vec![7.0, 2.0, -3.0]);
        // The hi-L row carries the original column entries, bit-copied…
        assert_eq!(dense_cols(&dense.a)[0], vec![1.0, 1.0, -2.0]);
        // …and the lo-G row is the same entries (no negation at this layer).
        assert_eq!(dense_cols(&dense.a)[1], vec![3.0, 3.0, 0.5]);
        assert_eq!(sp.row_map(), vec![[0, 1], [NO_ROW, 2]]);
    }

    #[test]
    fn row_map_covers_every_row_shape() {
        let mut b = CscBuilder::new(5);
        b.push_col(&[(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0)]);
        let sp = SparseLp::new(
            b.into_matrix(),
            // row shapes: equality, hi-only, lo-only, ranged, vacuous
            vec![3.0, f64::NEG_INFINITY, -1.0, 2.0, f64::NEG_INFINITY],
            vec![3.0, 9.0, f64::INFINITY, 7.0, f64::INFINITY],
            vec![0.0],
            vec![f64::INFINITY],
            vec![1.0],
            false,
        );
        let map = sp.row_map();
        assert_eq!(
            map,
            vec![
                [0, NO_ROW],
                [1, NO_ROW],
                [NO_ROW, 2],
                [3, 4],
                [NO_ROW, NO_ROW]
            ]
        );
        let dense = sp.to_lp();
        assert_eq!(dense.num_row, 5); // 1 + 1 + 1 + 2 + 0
        assert_eq!(dense.senses, vec!['E', 'L', 'G', 'L', 'G']);
        assert_eq!(dense.rhs, vec![3.0, 9.0, -1.0, 7.0, 2.0]);
        // Each mapped slot points at the dense row carrying the entries.
        for (i, slots) in map.iter().enumerate() {
            for &t in slots.iter() {
                if t != NO_ROW {
                    assert_eq!(
                        dense_cols(&dense.a)[0][t],
                        1.0,
                        "row {i} slot lost its entry"
                    );
                }
            }
        }
    }

    #[test]
    fn column_append_preserves_dense_prefix() {
        // I1 end-to-end: extending the matrix with new columns only appends —
        // the original dense block stays bit-identical at the same indices.
        let lp = crafted_lp();
        let sp = SparseLp::from_lp(&lp);
        let mut b = CscBuilder::from_matrix(&sp.a);
        b.push_col(&[(0, 2.5), (4, -1.5)]);
        b.push_col(&[]);
        let mut col_lower = sp.col_lower.clone();
        let mut col_upper = sp.col_upper.clone();
        let mut objective = sp.objective.clone();
        col_lower.extend([0.0, f64::NEG_INFINITY]);
        col_upper.extend([f64::INFINITY, 1.0]);
        objective.extend([0.25, -0.5]);
        let grown = SparseLp::new(
            b.into_matrix(),
            sp.row_lower.clone(),
            sp.row_upper.clone(),
            col_lower,
            col_upper,
            objective,
            sp.maximize,
        );
        let dense = grown.to_lp();
        assert_eq!(dense.num_col, lp.num_col + 2);
        assert_eq!(dense.num_row, lp.num_row);
        assert_eq!(dense.senses, lp.senses);
        let (grown_cols, lp_cols) = (dense_cols(&dense.a), dense_cols(&lp.a));
        for j in 0..lp.num_col {
            for i in 0..lp.num_row {
                assert_eq!(
                    grown_cols[j][i].to_bits(),
                    lp_cols[j][i].to_bits(),
                    "column {j} row {i} moved or changed under appends"
                );
            }
        }
    }

    #[test]
    fn new_accepts_stored_zeros_by_dropping_them() {
        // CscMatrix::from_triplets keeps stored zeros; SparseLp::new must
        // normalize them away (I3) without touching anything else.
        let a = CscMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (0, 0, -1.0), (1, 1, 2.0)]);
        assert_eq!(a.nnz(), 2); // summed duplicate ⇒ stored zero at (0,0)
        let sp = SparseLp::new(
            a,
            vec![0.0, 0.0],
            vec![1.0, 1.0],
            vec![0.0, 0.0],
            vec![1.0, 1.0],
            vec![0.0, 0.0],
            false,
        );
        assert_eq!(sp.a.nnz(), 1);
        assert_eq!(sp.a.get(0, 0), 0.0);
        assert_eq!(sp.a.get(1, 1), 2.0);
    }

    #[test]
    #[should_panic(expected = "row_upper")]
    fn new_rejects_inverted_row_bounds() {
        let mut b = CscBuilder::new(1);
        b.push_col(&[(0, 1.0)]);
        SparseLp::new(
            b.into_matrix(),
            vec![5.0],
            vec![3.0],
            vec![0.0],
            vec![1.0],
            vec![1.0],
            false,
        );
    }

    #[test]
    #[should_panic(expected = "col_lower")]
    fn new_rejects_nan_col_bounds() {
        let mut b = CscBuilder::new(1);
        b.push_col(&[(0, 1.0)]);
        SparseLp::new(
            b.into_matrix(),
            vec![f64::NEG_INFINITY],
            vec![f64::INFINITY],
            vec![f64::NAN],
            vec![1.0],
            vec![1.0],
            false,
        );
    }

    #[test]
    #[should_panic(expected = "ncol")]
    fn new_rejects_dimension_mismatch() {
        let mut b = CscBuilder::new(1);
        b.push_col(&[(0, 1.0)]);
        b.push_col(&[]);
        SparseLp::new(
            b.into_matrix(),
            vec![f64::NEG_INFINITY],
            vec![f64::INFINITY],
            vec![0.0], // one bound for two columns
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            false,
        );
    }

    #[test]
    #[should_panic(expected = "sense")]
    fn from_lp_rejects_unknown_sense() {
        let mut lp = crafted_lp();
        lp.senses[1] = 'N';
        SparseLp::from_lp(&lp);
    }
}
