//! In-memory optimization model shared by all file formats.
//!
//! [`Model::fill_hessian`] is translated from HiGHS `io/HMpsFF.cpp`
//! (`HMpsFF::fillHessian`).
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use rusolve_core::{CscBuilder, CscMatrix, SparseLp, Triplet};

use crate::gc::GcModel;

/// Constraint sense for one row of the constraint matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowSense {
    /// `L` row: `aᵀx ≤ rhs`
    LessEqual,
    /// `G` row: `aᵀx ≥ rhs`
    GreaterEqual,
    /// `E` row: `aᵀx = rhs`
    Equal,
}

/// Objective direction (`OBJSENSE`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ObjSense {
    #[default]
    Minimize,
    Maximize,
}

/// Bounds and integrality for one variable.
///
/// `lower`/`upper` may be infinite. This is the *resolved* view: MPS default
/// bounds (`0 ≤ x < ∞`) plus any `BOUNDS` directives have been applied.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VariableBound {
    pub lower: f64,
    pub upper: f64,
    pub integer: bool,
}

impl Default for VariableBound {
    /// MPS default: `0 ≤ x < ∞`, continuous.
    fn default() -> Self {
        Self {
            lower: 0.0,
            upper: f64::INFINITY,
            integer: false,
        }
    }
}

/// A linear (later MIP) model as loaded from or to be written to a file.
///
/// `N` (free) rows are not constraints: the first `N` row is the objective
/// (`obj_name`, its COLUMNS entries land in [`Model::objective`]); further
/// `N` rows have their entries ignored, matching common reader conventions.
/// `rhs` and `ranges` hold the *raw* MPS values (defaults 0.0); no sense
/// normalization of ranges is applied — that is the solver's interpretation.
#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    pub name: String,
    pub num_variables: usize,
    pub num_constraints: usize,
    pub obj_sense: ObjSense,
    /// Name of the `N` row used as objective, if any.
    pub obj_name: Option<String>,
    /// Objective constant: an `RHS` entry `v` on the objective
    /// (`N`) row, applied HiGHS-style as `-v` (HiGHS `io/HMpsFF.cpp:1065–1084`
    /// MIT: `obj_offset = -val`, carried as `lp.offset_` at :90 and added to
    /// the *reported* objective at `HighsInterface.cpp:3195`). First entry on
    /// the objective row wins, mirroring HiGHS's `has_obj_entry_` guard;
    /// `0.0` when the file carries none. **Reporting-only:** the solver-side
    /// views (`to_sparse_lp`, `MipProblem::from_model`) never read it, so
    /// simplex objectives and MIP bound/incumbent comparisons stay on raw
    /// `c·x` and no numeric behavior can shift.
    pub obj_offset: f64,
    /// Objective coefficients, length `num_variables`.
    pub objective: Vec<f64>,
    /// One sense per constraint row, length `num_constraints`.
    pub senses: Vec<RowSense>,
    pub rhs: Vec<f64>,
    pub ranges: Vec<f64>,
    pub bounds: Vec<VariableBound>,
    /// Constraint matrix only (objective row excluded), CSC.
    pub matrix: CscMatrix,
    pub row_names: Vec<String>,
    pub col_names: Vec<String>,
    /// Objective Hessian entries `(row, col, value)` in file order: HiGHS
    /// `HMpsFF::q_entries`, filled by the MPS `QUADOBJ` / `QMATRIX` sections
    /// and a `QSECTION` on the objective row (`HMpsFF.cpp` 1735-1889). Every
    /// off-diagonal entry of a triangular section (`QUADOBJ`, `QSECTION`) is
    /// also recorded mirrored, so this is always the full square matrix.
    /// Zero values are skipped and duplicates kept. Empty for an LP and for
    /// every LP-format model. [`Model::fill_hessian`] builds the CSC.
    pub hessian: Vec<Triplet>,
    // ---- general constraints ----
    /// Semi-continuous columns and general constraints, as the user wrote
    /// them (never their linearization). Empty for a plain LP/MIP, which then
    /// takes exactly the plain LP/MIP code path.
    pub gc: GcModel,
    // --- quadratic rows and cones, as the user wrote them.
    /// Empty for every LP/MIP/QP model. The legacy readers reject a model where it is not empty
    /// (the LP/MIP code never sees a dropped constraint); `read_mps_conic` / `read_lp_conic` keep it.
    pub conic: crate::conic::ConicBlock,
}

/// `HessianFormat::kSquare` (`model/HighsHessian.h`), the format of every
/// non-empty [`MpsHessian`].
pub const HESSIAN_FORMAT_SQUARE: i32 = 2;

/// The column-wise Hessian HiGHS's MPS reader builds from [`Model::hessian`]:
/// `HMpsFF::q_dim`, `q_start`, `q_index` and `q_value`. When `dim > 0` the
/// format is [`HESSIAN_FORMAT_SQUARE`]. When no entry was read, `dim` is 0 and
/// the vectors are empty, as `fillHessian` leaves them.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MpsHessian {
    /// `q_dim`: `num_variables` when any entry was read, else 0.
    pub dim: i32,
    /// `q_start`, `dim + 1` entries (with the end sentinel).
    pub start: Vec<i32>,
    /// `q_index`: row indices, each column's entries in file order.
    pub index: Vec<i32>,
    /// `q_value`.
    pub value: Vec<f64>,
}

impl Model {
    /// HiGHS: HMpsFF.cpp::HMpsFF::fillHessian (lines 177-216)
    ///
    /// `num_col` is `num_variables`. The C++ status return is always 0, so
    /// it is dropped.
    pub fn fill_hessian(&self) -> MpsHessian {
        let q_entries = &self.hessian;
        let num_col = self.num_variables as i32;
        let num_entries: usize = q_entries.len();
        if num_entries == 0 {
            // q_dim = 0, and q_start / q_index / q_value are left empty.
            return MpsHessian {
                dim: 0,
                ..MpsHessian::default()
            };
        }
        let q_dim: i32 = num_col;

        let mut q_start: Vec<i32> = vec![0; (q_dim + 1) as usize];
        let mut q_index: Vec<i32> = vec![0; num_entries];
        let mut q_value: Vec<f64> = vec![0.0; num_entries];

        // Use q_length to determine the number of entries in each column,
        // and then as workspace to point to the next entry to be filled in
        // each column
        let mut q_length: Vec<i32> = vec![0; q_dim as usize];

        for &(_, i_col, _) in q_entries.iter() {
            q_length[i_col] += 1;
        }
        q_start[0] = 0;
        for i_col in 0..num_col {
            let c = i_col as usize;
            q_start[c + 1] = q_start[c] + q_length[c];
            q_length[c] = q_start[c];
        }

        for &(i_row, i_col, value) in q_entries.iter() {
            q_index[q_length[i_col] as usize] = i_row as i32;
            q_value[q_length[i_col] as usize] = value;
            q_length[i_col] += 1;
        }
        // q_format = HessianFormat::kSquare (HESSIAN_FORMAT_SQUARE)
        MpsHessian {
            dim: q_dim,
            start: q_start,
            index: q_index,
            value: q_value,
        }
    }

    /// The solver-space LP of this model:
    /// [`SparseLp`] with every row single-sided, ranged MPS rows expanded as
    /// TWO rows — the base row in place, the complementary row appended at
    /// the tail — and matrix entries replicated into both. This is row-for-
    /// row, entry-for-entry the expansion the former dense CLI `model_to_lp`
    /// built, so `to_sparse_lp().to_dense()` reproduces that `Lp`'s senses,
    /// rhs, bounds and entries bit-identically (the historical row order is
    /// load-bearing for iteration-count differentials).
    ///
    /// Ranged-row semantics (matching HiGHS and `MipProblem::from_model`):
    /// `L b|r|` ⇒ also `G b−|r|`; `G b|r|` ⇒ also `L b+|r|`; `E b|r|` ⇒
    /// two-sided with the sign of `r` picking the sides.
    ///
    /// `Err` on integer variables (the LP pipeline is continuous-only).
    pub fn to_sparse_lp(&self) -> Result<SparseLp, String> {
        if self.bounds.iter().any(|b| b.integer) {
            return Err(
                "integer variables not supported yet (the LP pipeline is continuous-only)".into(),
            );
        }
        if let Some(what) = self.gc.first_general() {
            return Err(format!("{what} is not supported by the LP pipeline"));
        }
        if !self.gc.semi.is_empty() {
            return Err("semi-continuous columns are not supported by the LP pipeline".into());
        }
        let inf = f64::INFINITY;
        let n = self.num_variables;
        let m = self.num_constraints;
        let mut senses: Vec<char> =
            Vec::with_capacity(m + self.ranges.iter().filter(|&&r| r != 0.0).count());
        let mut rhs: Vec<f64> = Vec::with_capacity(senses.capacity());
        // (source row, complementary sense, rhs) — appended as tail rows.
        let mut extra: Vec<(usize, char, f64)> = Vec::new();
        for i in 0..m {
            let base = match self.senses[i] {
                RowSense::LessEqual => 'L',
                RowSense::GreaterEqual => 'G',
                RowSense::Equal => 'E',
            };
            senses.push(base);
            rhs.push(self.rhs[i]);
            let r = self.ranges[i];
            if r != 0.0 {
                match base {
                    'L' => extra.push((i, 'G', self.rhs[i] - r.abs())),
                    'G' => extra.push((i, 'L', self.rhs[i] + r.abs())),
                    'E' if r > 0.0 => {
                        senses[i] = 'G';
                        extra.push((i, 'L', self.rhs[i] + r));
                    }
                    'E' => {
                        senses[i] = 'L';
                        extra.push((i, 'G', self.rhs[i] + r));
                    }
                    _ => unreachable!(),
                }
            }
        }
        for &(row, sense, v) in &extra {
            let _ = row;
            senses.push(sense);
            rhs.push(v);
        }
        let nr = senses.len();
        let mut row_lower = Vec::with_capacity(nr);
        let mut row_upper = Vec::with_capacity(nr);
        for (i, &s) in senses.iter().enumerate() {
            let (lo, hi) = match s {
                'L' => (-inf, rhs[i]),
                'G' => (rhs[i], inf),
                'E' => (rhs[i], rhs[i]),
                other => return Err(format!("row {i} sense '{other}'")),
            };
            row_lower.push(lo);
            row_upper.push(hi);
        }
        // CSC → CSC with the ranged entries replicated into their extra row.
        // Column entries stay row-ascending: base entries first (all < m),
        // then each entry's replicated extra row — extras of ascending
        // source rows are ascending (m + k tracks row order).
        let (ptr, idx, val) = (
            self.matrix.col_ptr(),
            self.matrix.row_idx(),
            self.matrix.values(),
        );
        // Source row → its extra row (dense lookup; at most one extra per row).
        let mut extra_of = vec![usize::MAX; m];
        for (kx, &(row, _, _)) in extra.iter().enumerate() {
            extra_of[row] = m + kx;
        }
        let mut b = CscBuilder::new(nr);
        let mut col: Vec<(usize, f64)> = Vec::new();
        for j in 0..n {
            col.clear();
            // Base entries (all rows < m, ascending by the CSC layout) …
            for k in ptr[j]..ptr[j + 1] {
                col.push((idx[k] as usize, val[k]));
            }
            // … then each entry's replicated extra row, in the same
            // source-row order (≥ m, ascending), keeping the column sorted.
            for k in ptr[j]..ptr[j + 1] {
                let i = idx[k] as usize;
                if extra_of[i] != usize::MAX {
                    col.push((extra_of[i], val[k]));
                }
            }
            b.push_col(&col);
        }
        Ok(SparseLp::new(
            b.into_matrix(),
            row_lower,
            row_upper,
            self.bounds.iter().map(|b| b.lower).collect(),
            self.bounds.iter().map(|b| b.upper).collect(),
            self.objective.clone(),
            self.obj_sense == ObjSense::Maximize,
        ))
    }
}
