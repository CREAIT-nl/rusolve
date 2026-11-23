//! Append-only model builder: the Rust core of the Python `rusolve.Model` API
//! (docs/design/python-api.md).
//!
//! Adding a column or a row only appends to buffers. Values are checked when they arrive
//! (NaN, infinite coefficients, column indices out of range). Duplicates are merged and exact
//! zeros dropped when a row is finished, in O(row length), with a reusable marker array. The
//! whole model becomes a [`HighsLp`] in one O(nnz) pass ([`ModelBuilder::take`]). Nothing is
//! ever O(num_col) per row.
//!
//! Persistence: [`ModelBuilder::take`] moves everything added since the previous `take` out
//! of the builder (the solver then owns it, so the model is never stored twice). The builder
//! keeps only the counts, so new columns and rows keep their global indices. The first `take`
//! is the whole model ([`Batch::into_lp`] for `Highs::pass_model`); later ones are appended
//! with `Highs::add_cols` / `Highs::add_rows`. Entities already taken are "synced": changing
//! them is the caller's job (through the `Highs` change API). The `set_*` methods here change
//! only unsynced entities and return [`BuildError::Synced`] otherwise.

//! Provenance: the duplicate merge with a marker array and the counting-sort transpose are
//! textbook sparse techniques (T. Davis, "Direct Methods for Sparse Linear Systems", SIAM 2006,
//! ch. 2 on transpose and ch. 2.4 on summing duplicates). The `c`/`r` name prefixes are HiGHS's
//! (HConst.h `kHighsMinimalColNamePrefix`, `kHighsMinimalRowNamePrefix`).

use std::fmt;
use std::ops::Range;

use rusolve_ekk::HighsInt;
use rusolve_ekk::lp_data::hconst::{HighsVarType, MatrixFormat, ObjSense};
use rusolve_ekk::lp_data::highs_lp::HighsLp;
use rusolve_ekk::util::sparse_matrix::HighsSparseMatrix;

/// A per-column (or per-row) value: one scalar for all, or one value each.
#[derive(Clone, Copy, Debug)]
pub enum Fill<'a> {
    Scalar(f64),
    Slice(&'a [f64]),
}

/// Integrality for a block of columns: one flag for all, or one flag each.
#[derive(Clone, Copy, Debug)]
pub enum FillBool<'a> {
    Scalar(bool),
    Slice(&'a [bool]),
}

/// Names for a block of columns or rows.
#[derive(Clone, Debug)]
pub enum Names {
    /// One name each (length must equal the block length).
    List(Vec<String>),
    /// `"<prefix>[k]"` for k = 0..len (k counts within the block).
    Prefix(String),
}

/// Objective sense.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Sense {
    #[default]
    Minimize,
    Maximize,
}

/// Index type of an incoming sparse matrix (scipy gives int32 or int64).
pub trait Idx: Copy {
    fn to_i64(self) -> i64;
}
impl Idx for i32 {
    #[inline]
    fn to_i64(self) -> i64 {
        self as i64
    }
}
impl Idx for i64 {
    #[inline]
    fn to_i64(self) -> i64 {
        self
    }
}

/// What can go wrong; every variant names the row and/or column.
#[derive(Clone, Debug, PartialEq)]
pub enum BuildError {
    /// A value is NaN. `what` is e.g. "lower bound", "cost", "coefficient", "row lower bound".
    Nan {
        what: &'static str,
        row: Option<usize>,
        col: Option<usize>,
    },
    /// A matrix coefficient or a cost is ±inf.
    Infinite {
        what: &'static str,
        row: Option<usize>,
        col: Option<usize>,
    },
    /// A column index outside 0..num_cols (or negative).
    ColumnOutOfRange {
        row: Option<usize>,
        col: i64,
        num_cols: usize,
    },
    /// A row index outside 0..num_rows.
    RowOutOfRange { row: usize, num_rows: usize },
    /// A slice whose length is not the block length.
    Length {
        what: &'static str,
        expected: usize,
        got: usize,
    },
    /// The entity was already taken by the solver; change it through the solver instead.
    Synced { what: &'static str, index: usize },
    /// Anything else, e.g. a malformed sparse matrix (decreasing indptr).
    Invalid(String),
}

fn location(f: &mut fmt::Formatter<'_>, row: Option<usize>, col: Option<usize>) -> fmt::Result {
    match (col, row) {
        (Some(c), Some(r)) => write!(f, " of column {c} in row {r}"),
        (Some(c), None) => write!(f, " of column {c}"),
        (None, Some(r)) => write!(f, " of row {r}"),
        (None, None) => Ok(()),
    }
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BuildError::Nan { what, row, col } => {
                write!(f, "{what}")?;
                location(f, *row, *col)?;
                write!(f, " is NaN")
            }
            BuildError::Infinite { what, row, col } => {
                write!(f, "{what}")?;
                location(f, *row, *col)?;
                write!(f, " is infinite")
            }
            BuildError::ColumnOutOfRange { row, col, num_cols } => {
                write!(f, "column index {col}")?;
                if let Some(r) = row {
                    write!(f, " in row {r}")?;
                }
                write!(f, " is out of range: the model has {num_cols} columns")
            }
            BuildError::RowOutOfRange { row, num_rows } => {
                write!(
                    f,
                    "row index {row} is out of range: the model has {num_rows} rows"
                )
            }
            BuildError::Length {
                what,
                expected,
                got,
            } => write!(f, "{what} has length {got}, expected {expected}"),
            BuildError::Synced { what, index } => write!(
                f,
                "{what} {index} was already passed to the solver and cannot be changed here"
            ),
            BuildError::Invalid(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for BuildError {}

/// New columns since the previous [`ModelBuilder::take`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ColBatch {
    /// Global index of the first column in this batch.
    pub first: usize,
    pub cost: Vec<f64>,
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
    /// Empty if every column of the whole model so far is continuous.
    pub integrality: Vec<HighsVarType>,
    /// Empty if no column of the model has ever been named; otherwise one per column, with
    /// HiGHS's default `c<j>` (global j; HiGHS HConst.h kHighsMinimalColNamePrefix) for unnamed ones.
    pub names: Vec<String>,
}

/// New rows since the previous [`ModelBuilder::take`], row-wise, duplicates merged, zeros
/// dropped, column indices < the model's column count at take time.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RowBatch {
    /// Global index of the first row in this batch.
    pub first: usize,
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
    /// Row starts, length = rows + 1 (so `start[r]..start[r+1]` are row r's entries).
    pub start: Vec<HighsInt>,
    pub index: Vec<HighsInt>,
    pub value: Vec<f64>,
    /// As [`ColBatch::names`], with `r<i>` defaults.
    pub names: Vec<String>,
}

/// Everything added since the previous take.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Batch {
    pub cols: ColBatch,
    pub rows: RowBatch,
    pub sense: Sense,
    pub offset: f64,
    pub name: String,
}

impl Batch {
    /// The whole model as a `HighsLp` (only meaningful for the first take, when
    /// `cols.first == 0 && rows.first == 0`). The matrix is column-wise (`kColwise`), built
    /// by a counting-sort transpose in O(nnz), with row indices ascending within each column,
    /// which is the layout the MPS reader produces.
    pub fn into_lp(self) -> HighsLp {
        let Batch {
            cols,
            rows,
            sense,
            offset,
            name,
        } = self;
        debug_assert!(
            cols.first == 0 && rows.first == 0,
            "into_lp needs the first take"
        );
        let ncols = cols.cost.len();
        let nrows = rows.lower.len();
        let nnz = rows.index.len();
        // counting sort: column starts, then scatter rows in ascending order
        let mut start: Vec<HighsInt> = vec![0; ncols + 1];
        for &j in &rows.index {
            start[j as usize + 1] += 1;
        }
        for j in 0..ncols {
            start[j + 1] += start[j];
        }
        let mut next: Vec<HighsInt> = start[..ncols].to_vec();
        let mut index: Vec<HighsInt> = vec![0; nnz];
        let mut value: Vec<f64> = vec![0.0; nnz];
        for r in 0..nrows {
            let (a, b) = (rows.start[r] as usize, rows.start[r + 1] as usize);
            for p in a..b {
                let j = rows.index[p] as usize;
                let q = next[j] as usize;
                next[j] += 1;
                index[q] = r as HighsInt;
                value[q] = rows.value[p];
            }
        }
        let a_matrix = HighsSparseMatrix {
            format_: MatrixFormat::kColwise,
            num_col_: ncols as HighsInt,
            num_row_: nrows as HighsInt,
            start_: start,
            p_end_: Vec::new(),
            index_: index,
            value_: value,
        };
        HighsLp {
            num_col_: ncols as HighsInt,
            num_row_: nrows as HighsInt,
            col_cost_: cols.cost,
            col_lower_: cols.lower,
            col_upper_: cols.upper,
            row_lower_: rows.lower,
            row_upper_: rows.upper,
            a_matrix_: a_matrix,
            sense_: match sense {
                Sense::Minimize => ObjSense::kMinimize,
                Sense::Maximize => ObjSense::kMaximize,
            },
            offset_: offset,
            model_name_: name,
            col_names_: cols.names,
            row_names_: rows.names,
            integrality_: cols.integrality,
            ..Default::default()
        }
    }
}

/// The builder. See the module docs.
#[derive(Debug)]
pub struct ModelBuilder {
    // columns not yet taken
    col_cost: Vec<f64>,
    col_lower: Vec<f64>,
    col_upper: Vec<f64>,
    col_integer: Vec<bool>,
    /// Unsynced columns (local index) that are semi-continuous (`false`) or
    /// semi-integer (`true`); they override `col_integer` in `take`.
    col_semi: Vec<(usize, bool)>,
    col_names: Vec<String>, // "" = unnamed; empty Vec while no column was ever named
    // rows not yet taken (row-wise; the open row's entries are at the tail of index/value)
    row_lower: Vec<f64>,
    row_upper: Vec<f64>,
    row_start: Vec<usize>, // len = untaken rows + 1
    row_index: Vec<i32>,
    row_value: Vec<f64>,
    row_names: Vec<String>,
    // counts already taken
    synced_cols: usize,
    synced_rows: usize,
    any_integer: bool,
    any_col_name: bool,
    any_row_name: bool,
    // marker for duplicate merging: mark[j] = position+1 of column j in the open row, else 0
    mark: Vec<u32>,
    sense: Sense,
    offset: f64,
    name: String,
}

impl Default for ModelBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cold]
fn bad_value(what: &'static str, v: f64, row: Option<usize>, col: Option<usize>) -> BuildError {
    if v.is_nan() {
        BuildError::Nan { what, row, col }
    } else {
        BuildError::Infinite { what, row, col }
    }
}

#[inline]
fn fill_at(f: &Fill<'_>, i: usize) -> f64 {
    match f {
        Fill::Scalar(v) => *v,
        Fill::Slice(s) => s[i],
    }
}

fn check_fill_len(what: &'static str, f: &Fill<'_>, n: usize) -> Result<(), BuildError> {
    if let Fill::Slice(s) = f
        && s.len() != n
    {
        return Err(BuildError::Length {
            what,
            expected: n,
            got: s.len(),
        });
    }
    Ok(())
}

/// Map matrix column `k` to a model column and range-check it.
#[inline]
fn map_col(
    col_map: Option<&[i32]>,
    k: i64,
    num_cols: usize,
    row: Option<usize>,
) -> Result<usize, BuildError> {
    match col_map {
        Some(m) => {
            if (k as u64) >= m.len() as u64 {
                return Err(BuildError::ColumnOutOfRange {
                    row,
                    col: k,
                    num_cols: m.len(),
                });
            }
            let j = m[k as usize] as i64;
            if (j as u64) >= num_cols as u64 {
                return Err(BuildError::ColumnOutOfRange {
                    row,
                    col: j,
                    num_cols,
                });
            }
            Ok(j as usize)
        }
        None => {
            if (k as u64) >= num_cols as u64 {
                return Err(BuildError::ColumnOutOfRange {
                    row,
                    col: k,
                    num_cols,
                });
            }
            Ok(k as usize)
        }
    }
}

/// A source of (row, matrix column, value) triplets, visited twice by the transpose (count,
/// then scatter). Monomorphised, so the visit is a plain loop.
trait Triplets {
    fn visit<F: FnMut(i64, i64, f64) -> Result<(), BuildError>>(
        &self,
        f: &mut F,
    ) -> Result<(), BuildError>;
}

struct CscSrc<'a, P: Idx, I: Idx> {
    indptr: &'a [P],
    indices: &'a [I],
    data: &'a [f64],
}

impl<P: Idx, I: Idx> Triplets for CscSrc<'_, P, I> {
    #[inline]
    fn visit<F: FnMut(i64, i64, f64) -> Result<(), BuildError>>(
        &self,
        f: &mut F,
    ) -> Result<(), BuildError> {
        for k in 0..self.indptr.len() - 1 {
            let (a, b) = (
                self.indptr[k].to_i64() as usize,
                self.indptr[k + 1].to_i64() as usize,
            );
            for p in a..b {
                f(self.indices[p].to_i64(), k as i64, self.data[p])?;
            }
        }
        Ok(())
    }
}

struct CooSrc<'a, I: Idx, J: Idx> {
    row: &'a [I],
    col: &'a [J],
    data: &'a [f64],
}

impl<I: Idx, J: Idx> Triplets for CooSrc<'_, I, J> {
    #[inline]
    fn visit<F: FnMut(i64, i64, f64) -> Result<(), BuildError>>(
        &self,
        f: &mut F,
    ) -> Result<(), BuildError> {
        for p in 0..self.data.len() {
            f(self.row[p].to_i64(), self.col[p].to_i64(), self.data[p])?;
        }
        Ok(())
    }
}

impl ModelBuilder {
    pub fn new() -> Self {
        ModelBuilder {
            col_cost: Vec::new(),
            col_lower: Vec::new(),
            col_upper: Vec::new(),
            col_integer: Vec::new(),
            col_semi: Vec::new(),
            col_names: Vec::new(),
            row_lower: Vec::new(),
            row_upper: Vec::new(),
            row_start: vec![0],
            row_index: Vec::new(),
            row_value: Vec::new(),
            row_names: Vec::new(),
            synced_cols: 0,
            synced_rows: 0,
            any_integer: false,
            any_col_name: false,
            any_row_name: false,
            mark: Vec::new(),
            sense: Sense::Minimize,
            offset: 0.0,
            name: String::new(),
        }
    }

    /// Total columns (taken + not taken).
    #[inline]
    pub fn num_cols(&self) -> usize {
        self.synced_cols + self.col_cost.len()
    }
    /// Total rows (taken + not taken), not counting an open row.
    #[inline]
    pub fn num_rows(&self) -> usize {
        self.synced_rows + self.row_lower.len()
    }
    /// Nonzeros in rows not yet taken (after merging for finished rows).
    pub fn num_pending_nonzeros(&self) -> usize {
        *self.row_start.last().unwrap()
    }
    pub fn synced_cols(&self) -> usize {
        self.synced_cols
    }
    pub fn synced_rows(&self) -> usize {
        self.synced_rows
    }

    // ---------------------------------------------------------------- columns

    /// Append `n` columns. Checks: no NaN in any input, cost finite, slice lengths == n.
    /// Bounds may be ±inf; lower > upper is allowed (the model is then infeasible).
    /// Returns the global index range of the new columns.
    pub fn add_cols(
        &mut self,
        n: usize,
        cost: Fill<'_>,
        lower: Fill<'_>,
        upper: Fill<'_>,
        integer: FillBool<'_>,
    ) -> Result<Range<usize>, BuildError> {
        check_fill_len("cost", &cost, n)?;
        check_fill_len("lower bound", &lower, n)?;
        check_fill_len("upper bound", &upper, n)?;
        if let FillBool::Slice(s) = integer
            && s.len() != n
        {
            return Err(BuildError::Length {
                what: "integrality",
                expected: n,
                got: s.len(),
            });
        }
        let first = self.num_cols();
        if first + n > i32::MAX as usize {
            return Err(BuildError::Invalid(format!(
                "adding {n} columns to {first} exceeds the limit of {} columns",
                i32::MAX
            )));
        }
        for k in 0..n {
            let c = fill_at(&cost, k);
            if !c.is_finite() {
                return Err(bad_value("cost", c, None, Some(first + k)));
            }
            if fill_at(&lower, k).is_nan() {
                return Err(BuildError::Nan {
                    what: "lower bound",
                    row: None,
                    col: Some(first + k),
                });
            }
            if fill_at(&upper, k).is_nan() {
                return Err(BuildError::Nan {
                    what: "upper bound",
                    row: None,
                    col: Some(first + k),
                });
            }
            if matches!(&cost, Fill::Scalar(_))
                && matches!(&lower, Fill::Scalar(_))
                && matches!(&upper, Fill::Scalar(_))
            {
                break; // scalars: one check is enough
            }
        }
        match &cost {
            Fill::Scalar(v) => self.col_cost.resize(self.col_cost.len() + n, *v),
            Fill::Slice(s) => self.col_cost.extend_from_slice(s),
        }
        match &lower {
            Fill::Scalar(v) => self.col_lower.resize(self.col_lower.len() + n, *v),
            Fill::Slice(s) => self.col_lower.extend_from_slice(s),
        }
        match &upper {
            Fill::Scalar(v) => self.col_upper.resize(self.col_upper.len() + n, *v),
            Fill::Slice(s) => self.col_upper.extend_from_slice(s),
        }
        match integer {
            FillBool::Scalar(b) => {
                self.any_integer |= b && n > 0;
                self.col_integer.resize(self.col_integer.len() + n, b);
            }
            FillBool::Slice(s) => {
                self.any_integer |= s.iter().any(|&b| b);
                self.col_integer.extend_from_slice(s);
            }
        }
        if self.any_col_name {
            self.col_names.resize(self.col_cost.len(), String::new());
        }
        Ok(first..first + n)
    }

    fn check_unsynced_cols(&self, cols: &Range<usize>) -> Result<(), BuildError> {
        if cols.end > self.num_cols() {
            return Err(BuildError::ColumnOutOfRange {
                row: None,
                col: cols.end as i64 - 1,
                num_cols: self.num_cols(),
            });
        }
        if cols.start < self.synced_cols && cols.start < cols.end {
            return Err(BuildError::Synced {
                what: "column",
                index: cols.start,
            });
        }
        Ok(())
    }

    fn check_unsynced_col(&self, col: usize) -> Result<usize, BuildError> {
        if col >= self.num_cols() {
            return Err(BuildError::ColumnOutOfRange {
                row: None,
                col: col as i64,
                num_cols: self.num_cols(),
            });
        }
        if col < self.synced_cols {
            return Err(BuildError::Synced {
                what: "column",
                index: col,
            });
        }
        Ok(col - self.synced_cols)
    }

    fn check_unsynced_row(&self, row: usize) -> Result<usize, BuildError> {
        if row >= self.num_rows() {
            return Err(BuildError::RowOutOfRange {
                row,
                num_rows: self.num_rows(),
            });
        }
        if row < self.synced_rows {
            return Err(BuildError::Synced {
                what: "row",
                index: row,
            });
        }
        Ok(row - self.synced_rows)
    }

    fn make_names(len: usize, names: Names) -> Result<Vec<String>, BuildError> {
        match names {
            Names::List(v) => {
                if v.len() != len {
                    return Err(BuildError::Length {
                        what: "names",
                        expected: len,
                        got: v.len(),
                    });
                }
                Ok(v)
            }
            Names::Prefix(p) => Ok((0..len).map(|k| format!("{p}[{k}]")).collect()),
        }
    }

    /// Name a block of columns (global indices; must be unsynced).
    pub fn set_col_names(&mut self, cols: Range<usize>, names: Names) -> Result<(), BuildError> {
        self.check_unsynced_cols(&cols)?;
        let names = Self::make_names(cols.len(), names)?;
        if cols.is_empty() {
            return Ok(());
        }
        self.any_col_name = true;
        self.col_names.resize(self.col_cost.len(), String::new());
        let base = cols.start - self.synced_cols;
        for (k, nm) in names.into_iter().enumerate() {
            self.col_names[base + k] = nm;
        }
        Ok(())
    }

    /// The name of an unsynced column if it has one.
    pub fn col_name(&self, col: usize) -> Option<&str> {
        let i = col.checked_sub(self.synced_cols)?;
        self.col_names
            .get(i)
            .map(String::as_str)
            .filter(|s| !s.is_empty())
    }

    pub fn set_col_bounds(&mut self, col: usize, lower: f64, upper: f64) -> Result<(), BuildError> {
        let i = self.check_unsynced_col(col)?;
        if lower.is_nan() {
            return Err(BuildError::Nan {
                what: "lower bound",
                row: None,
                col: Some(col),
            });
        }
        if upper.is_nan() {
            return Err(BuildError::Nan {
                what: "upper bound",
                row: None,
                col: Some(col),
            });
        }
        self.col_lower[i] = lower;
        self.col_upper[i] = upper;
        Ok(())
    }
    pub fn set_col_cost(&mut self, col: usize, cost: f64) -> Result<(), BuildError> {
        let i = self.check_unsynced_col(col)?;
        if !cost.is_finite() {
            return Err(bad_value("cost", cost, None, Some(col)));
        }
        self.col_cost[i] = cost;
        Ok(())
    }
    pub fn set_col_integer(&mut self, col: usize, integer: bool) -> Result<(), BuildError> {
        let i = self.check_unsynced_col(col)?;
        self.col_integer[i] = integer;
        self.any_integer |= integer;
        self.col_semi.retain(|&(k, _)| k != i);
        Ok(())
    }
    /// Make an unsynced column semi-continuous (`integer = false`) or semi-integer.
    pub fn set_col_semi(&mut self, col: usize, integer: bool) -> Result<(), BuildError> {
        let i = self.check_unsynced_col(col)?;
        self.col_semi.retain(|&(k, _)| k != i);
        self.col_semi.push((i, integer));
        self.any_integer = true;
        Ok(())
    }
    /// Costs of all unsynced columns become 0 (for "replace the objective").
    pub fn clear_unsynced_costs(&mut self) {
        self.col_cost.fill(0.0);
    }

    // ---------------------------------------------------------------- rows, entry by entry

    /// Append one entry to the open row (opening it if none is open). Checks the value for
    /// NaN/inf and the column for range against `num_cols()`; on error the open row is left as
    /// is (call [`Self::abort_row`]). This is the hot path of dict ingest: keep it two pushes
    /// and two compares.
    #[inline]
    pub fn push(&mut self, col: i64, value: f64) -> Result<(), BuildError> {
        if (col as u64) >= self.num_cols() as u64 {
            return Err(self.push_col_error(col));
        }
        if !value.is_finite() {
            return Err(bad_value(
                "coefficient",
                value,
                Some(self.num_rows()),
                Some(col as usize),
            ));
        }
        self.row_index.push(col as i32);
        self.row_value.push(value);
        Ok(())
    }

    #[cold]
    fn push_col_error(&self, col: i64) -> BuildError {
        BuildError::ColumnOutOfRange {
            row: Some(self.num_rows()),
            col,
            num_cols: self.num_cols(),
        }
    }

    fn grow_mark(&mut self) {
        let n = self.num_cols();
        if self.mark.len() < n {
            self.mark.resize(n, 0);
        }
    }

    /// Merge duplicates (sum) and drop exact zeros of the raw entries `a..b` of the row buffers,
    /// writing the result from position `w` (`w <= a`). Returns the new end, and the first
    /// column whose merged value is not finite: `push` rejects NaN and inf, so this is a sum of
    /// finite duplicates that overflowed to ±inf, which the spec forbids like any inf coefficient. Marks are clear again on return. Davis
    /// 2006, the marker-array form of summing duplicates. Needs `mark.len() >= num_cols`.
    fn merge_range(&mut self, a: usize, b: usize, w: usize) -> (usize, Option<usize>) {
        let mut q = w;
        for p in a..b {
            let j = self.row_index[p];
            let v = self.row_value[p];
            let m = self.mark[j as usize];
            if m == 0 {
                self.row_index[q] = j;
                self.row_value[q] = v;
                self.mark[j as usize] = (q - w) as u32 + 1;
                q += 1;
            } else {
                self.row_value[w + m as usize - 1] += v;
            }
        }
        let mut z = w;
        let mut bad = None;
        for p in w..q {
            let j = self.row_index[p];
            let v = self.row_value[p];
            self.mark[j as usize] = 0;
            if v != 0.0 {
                if !v.is_finite() && bad.is_none() {
                    bad = Some(j as usize);
                }
                self.row_index[z] = j;
                self.row_value[z] = v;
                z += 1;
            }
        }
        (z, bad)
    }

    /// Close the open row (an empty row is fine) with bounds `lower..=upper`; merge duplicate
    /// columns (sum), drop exact zeros. Returns the global row index. NaN bounds are an error
    /// (and abort the row).
    pub fn finish_row(&mut self, lower: f64, upper: f64) -> Result<usize, BuildError> {
        let row = self.num_rows();
        if lower.is_nan() || upper.is_nan() {
            self.abort_row();
            return Err(BuildError::Nan {
                what: if lower.is_nan() {
                    "row lower bound"
                } else {
                    "row upper bound"
                },
                row: Some(row),
                col: None,
            });
        }
        let s = *self.row_start.last().unwrap();
        let len = self.row_index.len();
        self.grow_mark();
        let (end, bad) = self.merge_range(s, len, s);
        self.row_index.truncate(end);
        self.row_value.truncate(end);
        if let Some(col) = bad {
            self.row_index.truncate(s);
            self.row_value.truncate(s);
            return Err(BuildError::Infinite {
                what: "summed duplicate coefficient",
                row: Some(row),
                col: Some(col),
            });
        }
        self.row_start.push(end);
        self.row_lower.push(lower);
        self.row_upper.push(upper);
        if self.any_row_name {
            self.row_names.push(String::new());
        }
        Ok(row)
    }

    /// Discard the open row's entries.
    pub fn abort_row(&mut self) {
        let s = *self.row_start.last().unwrap();
        self.row_index.truncate(s);
        self.row_value.truncate(s);
    }

    /// `push` every pair, then `finish_row`; aborts the row on error.
    pub fn add_row(
        &mut self,
        index: &[i32],
        value: &[f64],
        lower: f64,
        upper: f64,
    ) -> Result<usize, BuildError> {
        if index.len() != value.len() {
            return Err(BuildError::Length {
                what: "row values",
                expected: index.len(),
                got: value.len(),
            });
        }
        for (&j, &v) in index.iter().zip(value) {
            if let Err(e) = self.push(j as i64, v) {
                self.abort_row();
                return Err(e);
            }
        }
        self.finish_row(lower, upper)
    }

    // ---------------------------------------------------------------- rows, whole matrices

    fn check_no_open_row(&self) -> Result<(), BuildError> {
        if self.row_index.len() != *self.row_start.last().unwrap() {
            return Err(BuildError::Invalid(
                "a row is open: finish_row or abort_row it before adding a matrix".to_string(),
            ));
        }
        Ok(())
    }

    fn check_row_fills(
        &self,
        nrows: usize,
        lower: &Fill<'_>,
        upper: &Fill<'_>,
    ) -> Result<(), BuildError> {
        check_fill_len("row lower bounds", lower, nrows)?;
        check_fill_len("row upper bounds", upper, nrows)?;
        let first = self.num_rows();
        for r in 0..nrows {
            if fill_at(lower, r).is_nan() {
                return Err(BuildError::Nan {
                    what: "row lower bound",
                    row: Some(first + r),
                    col: None,
                });
            }
            if fill_at(upper, r).is_nan() {
                return Err(BuildError::Nan {
                    what: "row upper bound",
                    row: Some(first + r),
                    col: None,
                });
            }
        }
        Ok(())
    }

    /// Row bounds and names of rows just appended to the entry buffers.
    fn push_row_meta(&mut self, nrows: usize, lower: &Fill<'_>, upper: &Fill<'_>) {
        match lower {
            Fill::Scalar(v) => self.row_lower.resize(self.row_lower.len() + nrows, *v),
            Fill::Slice(s) => self.row_lower.extend_from_slice(s),
        }
        match upper {
            Fill::Scalar(v) => self.row_upper.resize(self.row_upper.len() + nrows, *v),
            Fill::Slice(s) => self.row_upper.extend_from_slice(s),
        }
        if self.any_row_name {
            self.row_names.resize(self.row_lower.len(), String::new());
        }
    }

    /// Undo a failed bulk add.
    fn rollback(&mut self, entries: usize, rows_start_len: usize) {
        self.row_index.truncate(entries);
        self.row_value.truncate(entries);
        self.row_start.truncate(rows_start_len);
    }

    /// Append the rows of a CSR matrix (`indptr.len() == nrows + 1`). Column `k` of the matrix
    /// is model column `col_map[k]` if given, else `k`. `lower`/`upper` are per row. On error
    /// no row is added. Must not be slower than a few ns per nonzero.
    pub fn add_rows_csr<P: Idx, I: Idx>(
        &mut self,
        indptr: &[P],
        indices: &[I],
        data: &[f64],
        col_map: Option<&[i32]>,
        lower: Fill<'_>,
        upper: Fill<'_>,
    ) -> Result<Range<usize>, BuildError> {
        self.check_no_open_row()?;
        if indptr.is_empty() {
            return Err(BuildError::Invalid("indptr is empty".to_string()));
        }
        let nrows = indptr.len() - 1;
        if indices.len() != data.len() {
            return Err(BuildError::Length {
                what: "data",
                expected: indices.len(),
                got: data.len(),
            });
        }
        let mut prev = indptr[0].to_i64();
        if prev != 0 {
            return Err(BuildError::Invalid(format!(
                "indptr[0] is {prev}, expected 0"
            )));
        }
        for (r, p) in indptr.iter().enumerate().skip(1) {
            let cur = p.to_i64();
            if cur < prev || cur as u64 > indices.len() as u64 {
                return Err(BuildError::Invalid(format!(
                    "indptr is not non-decreasing within 0..={} at position {r} (value {cur})",
                    indices.len()
                )));
            }
            prev = cur;
        }
        if prev as usize != indices.len() {
            return Err(BuildError::Invalid(format!(
                "indptr ends at {prev} but there are {} entries",
                indices.len()
            )));
        }
        self.check_row_fills(nrows, &lower, &upper)?;
        let first = self.num_rows();
        if first + nrows > i32::MAX as usize {
            return Err(BuildError::Invalid("too many rows".to_string()));
        }
        let nc = self.num_cols();
        self.grow_mark();
        let entries0 = self.row_index.len();
        let starts0 = self.row_start.len();
        self.row_index.reserve(indices.len());
        self.row_value.reserve(indices.len());
        self.row_start.reserve(nrows);
        for r in 0..nrows {
            let (a, b) = (indptr[r].to_i64() as usize, indptr[r + 1].to_i64() as usize);
            let w = self.row_index.len();
            let mut dirty = false;
            let mut err = None;
            for p in a..b {
                let v = data[p];
                let j = match map_col(col_map, indices[p].to_i64(), nc, Some(first + r)) {
                    Ok(j) => j,
                    Err(e) => {
                        err = Some(e);
                        break;
                    }
                };
                if !v.is_finite() {
                    err = Some(bad_value("coefficient", v, Some(first + r), Some(j)));
                    break;
                }
                let m = &mut self.mark[j];
                dirty |= (*m != 0) | (v == 0.0);
                *m = 1;
                self.row_index.push(j as i32);
                self.row_value.push(v);
            }
            let len = self.row_index.len();
            for p in w..len {
                self.mark[self.row_index[p] as usize] = 0;
            }
            if let Some(e) = err {
                self.rollback(entries0, starts0);
                return Err(e);
            }
            if dirty {
                let (end, bad) = self.merge_range(w, len, w);
                self.row_index.truncate(end);
                self.row_value.truncate(end);
                if let Some(col) = bad {
                    self.rollback(entries0, starts0);
                    return Err(BuildError::Infinite {
                        what: "summed duplicate coefficient",
                        row: Some(first + r),
                        col: Some(col),
                    });
                }
            }
            self.row_start.push(self.row_index.len());
        }
        self.push_row_meta(nrows, &lower, &upper);
        Ok(first..first + nrows)
    }

    /// Transpose-ingest shared by CSC and COO: validate everything and count (pass 1), scatter
    /// row-wise (pass 2), merge each row in place.
    fn add_rows_triplets<T: Triplets>(
        &mut self,
        nrows: usize,
        src: &T,
        nnz: usize,
        col_map: Option<&[i32]>,
        lower: Fill<'_>,
        upper: Fill<'_>,
    ) -> Result<Range<usize>, BuildError> {
        self.check_no_open_row()?;
        self.check_row_fills(nrows, &lower, &upper)?;
        let first = self.num_rows();
        if first + nrows > i32::MAX as usize {
            return Err(BuildError::Invalid("too many rows".to_string()));
        }
        let nc = self.num_cols();
        let mut pos: Vec<usize> = vec![0; nrows + 1];
        src.visit(&mut |r, k, v| {
            if (r as u64) >= nrows as u64 {
                return Err(BuildError::RowOutOfRange {
                    row: r as usize,
                    num_rows: nrows,
                });
            }
            let row = Some(first + r as usize);
            let j = map_col(col_map, k, nc, row)?;
            if !v.is_finite() {
                return Err(bad_value("coefficient", v, row, Some(j)));
            }
            pos[r as usize + 1] += 1;
            Ok(())
        })?;
        for r in 0..nrows {
            pos[r + 1] += pos[r];
        }
        let raw_starts = pos.clone();
        let base = self.row_index.len();
        let starts0 = self.row_start.len();
        self.row_index.resize(base + nnz, 0);
        self.row_value.resize(base + nnz, 0.0);
        {
            let (index, value) = (&mut self.row_index, &mut self.row_value);
            src.visit(&mut |r, k, v| {
                let j = map_col(col_map, k, nc, None)?;
                let q = base + pos[r as usize];
                pos[r as usize] += 1;
                index[q] = j as i32;
                value[q] = v;
                Ok(())
            })?;
        }
        self.grow_mark();
        let mut w = base;
        for r in 0..nrows {
            let (a, b) = (base + raw_starts[r], base + raw_starts[r + 1]);
            let (end, bad) = self.merge_range(a, b, w);
            if let Some(col) = bad {
                self.rollback(base, starts0);
                return Err(BuildError::Infinite {
                    what: "summed duplicate coefficient",
                    row: Some(first + r),
                    col: Some(col),
                });
            }
            w = end;
            self.row_start.push(w);
        }
        self.row_index.truncate(w);
        self.row_value.truncate(w);
        self.push_row_meta(nrows, &lower, &upper);
        Ok(first..first + nrows)
    }

    /// The same for a CSC matrix with `nrows` rows (`indptr.len() == ncols + 1`), via a
    /// counting-sort transpose.
    #[allow(clippy::too_many_arguments)]
    pub fn add_rows_csc<P: Idx, I: Idx>(
        &mut self,
        nrows: usize,
        indptr: &[P],
        indices: &[I],
        data: &[f64],
        col_map: Option<&[i32]>,
        lower: Fill<'_>,
        upper: Fill<'_>,
    ) -> Result<Range<usize>, BuildError> {
        if indptr.is_empty() {
            return Err(BuildError::Invalid("indptr is empty".to_string()));
        }
        let ncols = indptr.len() - 1;
        if indices.len() != data.len() {
            return Err(BuildError::Length {
                what: "data",
                expected: indices.len(),
                got: data.len(),
            });
        }
        if let Some(m) = col_map
            && m.len() != ncols
        {
            return Err(BuildError::Length {
                what: "col_map",
                expected: ncols,
                got: m.len(),
            });
        }
        let mut prev = indptr[0].to_i64();
        if prev != 0 {
            return Err(BuildError::Invalid(format!(
                "indptr[0] is {prev}, expected 0"
            )));
        }
        for (k, p) in indptr.iter().enumerate().skip(1) {
            let cur = p.to_i64();
            if cur < prev || cur as u64 > indices.len() as u64 {
                return Err(BuildError::Invalid(format!(
                    "indptr is not non-decreasing within 0..={} at position {k} (value {cur})",
                    indices.len()
                )));
            }
            prev = cur;
        }
        if prev as usize != indices.len() {
            return Err(BuildError::Invalid(format!(
                "indptr ends at {prev} but there are {} entries",
                indices.len()
            )));
        }
        let src = CscSrc {
            indptr,
            indices,
            data,
        };
        self.add_rows_triplets(nrows, &src, indices.len(), col_map, lower, upper)
    }

    /// The same for COO triplets with `nrows` rows.
    #[allow(clippy::too_many_arguments)]
    pub fn add_rows_coo<I: Idx, J: Idx>(
        &mut self,
        nrows: usize,
        row: &[I],
        col: &[J],
        data: &[f64],
        col_map: Option<&[i32]>,
        lower: Fill<'_>,
        upper: Fill<'_>,
    ) -> Result<Range<usize>, BuildError> {
        if row.len() != data.len() || col.len() != data.len() {
            return Err(BuildError::Length {
                what: "row/col/data",
                expected: data.len(),
                got: if row.len() != data.len() {
                    row.len()
                } else {
                    col.len()
                },
            });
        }
        let src = CooSrc { row, col, data };
        self.add_rows_triplets(nrows, &src, data.len(), col_map, lower, upper)
    }

    /// The same for a dense row-major `nrows × ncols` array (nonzeros only are kept).
    pub fn add_rows_dense(
        &mut self,
        nrows: usize,
        ncols: usize,
        data: &[f64],
        col_map: Option<&[i32]>,
        lower: Fill<'_>,
        upper: Fill<'_>,
    ) -> Result<Range<usize>, BuildError> {
        self.check_no_open_row()?;
        if data.len() != nrows * ncols {
            return Err(BuildError::Length {
                what: "dense data",
                expected: nrows * ncols,
                got: data.len(),
            });
        }
        if let Some(m) = col_map
            && m.len() != ncols
        {
            return Err(BuildError::Length {
                what: "col_map",
                expected: ncols,
                got: m.len(),
            });
        }
        self.check_row_fills(nrows, &lower, &upper)?;
        let first = self.num_rows();
        if first + nrows > i32::MAX as usize {
            return Err(BuildError::Invalid("too many rows".to_string()));
        }
        let nc = self.num_cols();
        // validate the column map once, not per row
        let mapped: Vec<usize> = (0..ncols)
            .map(|k| map_col(col_map, k as i64, nc, None))
            .collect::<Result<_, _>>()?;
        self.grow_mark();
        let entries0 = self.row_index.len();
        let starts0 = self.row_start.len();
        for r in 0..nrows {
            let w = self.row_index.len();
            let mut err = None;
            for (k, &v) in data[r * ncols..(r + 1) * ncols].iter().enumerate() {
                if v == 0.0 {
                    continue;
                }
                if !v.is_finite() {
                    err = Some(bad_value(
                        "coefficient",
                        v,
                        Some(first + r),
                        Some(mapped[k]),
                    ));
                    break;
                }
                self.row_index.push(mapped[k] as i32);
                self.row_value.push(v);
            }
            if let Some(e) = err {
                self.rollback(entries0, starts0);
                return Err(e);
            }
            if col_map.is_some() {
                let len = self.row_index.len();
                let (end, bad) = self.merge_range(w, len, w);
                self.row_index.truncate(end);
                self.row_value.truncate(end);
                if let Some(col) = bad {
                    self.rollback(entries0, starts0);
                    return Err(BuildError::Infinite {
                        what: "summed duplicate coefficient",
                        row: Some(first + r),
                        col: Some(col),
                    });
                }
            }
            self.row_start.push(self.row_index.len());
        }
        self.push_row_meta(nrows, &lower, &upper);
        Ok(first..first + nrows)
    }

    pub fn set_row_names(&mut self, rows: Range<usize>, names: Names) -> Result<(), BuildError> {
        if rows.end > self.num_rows() {
            return Err(BuildError::RowOutOfRange {
                row: rows.end - 1,
                num_rows: self.num_rows(),
            });
        }
        if rows.start < self.synced_rows && rows.start < rows.end {
            return Err(BuildError::Synced {
                what: "row",
                index: rows.start,
            });
        }
        let names = Self::make_names(rows.len(), names)?;
        if rows.is_empty() {
            return Ok(());
        }
        self.any_row_name = true;
        self.row_names.resize(self.row_lower.len(), String::new());
        let base = rows.start - self.synced_rows;
        for (k, nm) in names.into_iter().enumerate() {
            self.row_names[base + k] = nm;
        }
        Ok(())
    }
    pub fn row_name(&self, row: usize) -> Option<&str> {
        let i = row.checked_sub(self.synced_rows)?;
        self.row_names
            .get(i)
            .map(String::as_str)
            .filter(|s| !s.is_empty())
    }
    pub fn set_row_bounds(&mut self, row: usize, lower: f64, upper: f64) -> Result<(), BuildError> {
        let i = self.check_unsynced_row(row)?;
        if lower.is_nan() {
            return Err(BuildError::Nan {
                what: "row lower bound",
                row: Some(row),
                col: None,
            });
        }
        if upper.is_nan() {
            return Err(BuildError::Nan {
                what: "row upper bound",
                row: Some(row),
                col: None,
            });
        }
        self.row_lower[i] = lower;
        self.row_upper[i] = upper;
        Ok(())
    }
    /// Set (or insert, or with 0.0 remove) the coefficient of an unsynced row. Finding the
    /// entry is O(row length); inserting or removing also shifts the entries after it (a
    /// `Vec` insert), so those two are O(pending nonzeros) in the worst case.
    pub fn set_coef(&mut self, row: usize, col: usize, value: f64) -> Result<(), BuildError> {
        let i = self.check_unsynced_row(row)?;
        if col >= self.num_cols() {
            return Err(BuildError::ColumnOutOfRange {
                row: Some(row),
                col: col as i64,
                num_cols: self.num_cols(),
            });
        }
        if !value.is_finite() {
            return Err(bad_value("coefficient", value, Some(row), Some(col)));
        }
        let (a, b) = (self.row_start[i], self.row_start[i + 1]);
        let pos = self.row_index[a..b].iter().position(|&j| j as usize == col);
        match (pos, value == 0.0) {
            (Some(p), false) => self.row_value[a + p] = value,
            (Some(p), true) => {
                self.row_index.remove(a + p);
                self.row_value.remove(a + p);
                for s in &mut self.row_start[i + 1..] {
                    *s -= 1;
                }
            }
            (None, true) => {}
            (None, false) => {
                self.row_index.insert(b, col as i32);
                self.row_value.insert(b, value);
                for s in &mut self.row_start[i + 1..] {
                    *s += 1;
                }
            }
        }
        Ok(())
    }

    // ---------------------------------------------------------------- objective, model

    pub fn set_sense(&mut self, sense: Sense) {
        self.sense = sense;
    }
    pub fn sense(&self) -> Sense {
        self.sense
    }
    pub fn set_offset(&mut self, offset: f64) {
        self.offset = offset;
    }
    pub fn offset(&self) -> f64 {
        self.offset
    }
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    // ---------------------------------------------------------------- hand-over

    /// Move everything not yet taken out of the builder (see the module docs). Panics if a row
    /// is open. Names: if any column (row) of the model was ever named, every column (row) in
    /// the batch gets a name, with `c<j>` (`r<i>`) defaults. Integrality is empty unless some
    /// column of the model was ever integer.
    ///
    /// A column (row) taken before the first name was set stays unnamed in the solver: names
    /// only reach entities taken after the naming.
    pub fn take(&mut self) -> Batch {
        assert!(
            self.row_index.len() == *self.row_start.last().unwrap(),
            "take() with an open row"
        );
        let first_col = self.synced_cols;
        let first_row = self.synced_rows;
        let ncols = self.col_cost.len();
        let nrows = self.row_lower.len();
        let names = if self.any_col_name {
            let mut v = std::mem::take(&mut self.col_names);
            v.resize(ncols, String::new());
            for (k, n) in v.iter_mut().enumerate() {
                if n.is_empty() {
                    *n = format!("c{}", first_col + k);
                }
            }
            v
        } else {
            Vec::new()
        };
        let integrality = if self.any_integer {
            let mut t: Vec<HighsVarType> = std::mem::take(&mut self.col_integer)
                .into_iter()
                .map(|b| {
                    if b {
                        HighsVarType::kInteger
                    } else {
                        HighsVarType::kContinuous
                    }
                })
                .collect();
            for (i, integer) in self.col_semi.drain(..) {
                t[i] = if integer {
                    HighsVarType::kSemiInteger
                } else {
                    HighsVarType::kSemiContinuous
                };
            }
            t
        } else {
            self.col_integer.clear();
            Vec::new()
        };
        let cols = ColBatch {
            first: first_col,
            cost: std::mem::take(&mut self.col_cost),
            lower: std::mem::take(&mut self.col_lower),
            upper: std::mem::take(&mut self.col_upper),
            integrality,
            names,
        };
        let rnames = if self.any_row_name {
            let mut v = std::mem::take(&mut self.row_names);
            v.resize(nrows, String::new());
            for (k, n) in v.iter_mut().enumerate() {
                if n.is_empty() {
                    *n = format!("r{}", first_row + k);
                }
            }
            v
        } else {
            Vec::new()
        };
        // HighsInt is i32: more than i32::MAX nonzeros cannot be passed to HiGHS anyway.
        let start: Vec<HighsInt> = self.row_start.iter().map(|&s| s as HighsInt).collect();
        self.row_start.clear();
        self.row_start.push(0);
        let rows = RowBatch {
            first: first_row,
            lower: std::mem::take(&mut self.row_lower),
            upper: std::mem::take(&mut self.row_upper),
            start,
            index: std::mem::take(&mut self.row_index),
            value: std::mem::take(&mut self.row_value),
            names: rnames,
        };
        self.synced_cols += ncols;
        self.synced_rows += nrows;
        Batch {
            cols,
            rows,
            sense: self.sense,
            offset: self.offset,
            name: self.name.clone(),
        }
    }
}

#[cfg(test)]
mod tests;
mod truncate;
