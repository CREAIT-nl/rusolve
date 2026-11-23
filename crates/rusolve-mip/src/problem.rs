//! The solver's MIP form, built from a file-level [`rusolve_io::Model`].

use rusolve_core::{CscBuilder, CscMatrix};
use rusolve_io::{Model, ObjSense, RowSense};

/// `min costᵀx  s.t.  row_lower ≤ Ax ≤ row_upper,  col_lower ≤ x ≤ col_upper,
/// x_j ∈ ℤ for integer[j]`.
///
/// Maximization models are stored negated (`cost = −c`); [`Self::caller_objective`]
/// maps internal values back.
#[derive(Debug, Clone)]
pub struct MipProblem {
    pub num_col: usize,
    pub num_row: usize,
    /// Constraint matrix, `num_row × num_col`.
    pub matrix: CscMatrix,
    pub row_lower: Vec<f64>,
    pub row_upper: Vec<f64>,
    /// Minimization-form objective.
    pub cost: Vec<f64>,
    /// The caller's model maximizes (internal costs are negated).
    pub maximize: bool,
    pub col_lower: Vec<f64>,
    pub col_upper: Vec<f64>,
    pub integer: Vec<bool>,
}

/// Worst violations of a candidate point, for acceptance checks and tests.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Violation {
    pub bound: f64,
    pub row: f64,
    pub integrality: f64,
}

impl MipProblem {
    /// Row bounds follow MPS semantics (matching HiGHS and the LP CLI):
    /// `L b |r|` ⇒ `[b−|r|, b]`, `G b |r|` ⇒ `[b, b+|r|]`, `E b r` ⇒ `[b, b+r]`
    /// for r > 0 and `[b+r, b]` for r < 0.
    pub fn from_model(m: &Model) -> Self {
        let inf = f64::INFINITY;
        let mut row_lower = Vec::with_capacity(m.num_constraints);
        let mut row_upper = Vec::with_capacity(m.num_constraints);
        for i in 0..m.num_constraints {
            let (b, r) = (m.rhs[i], m.ranges[i]);
            let (lo, hi) = match m.senses[i] {
                RowSense::LessEqual if r != 0.0 => (b - r.abs(), b),
                RowSense::LessEqual => (-inf, b),
                RowSense::GreaterEqual if r != 0.0 => (b, b + r.abs()),
                RowSense::GreaterEqual => (b, inf),
                RowSense::Equal if r > 0.0 => (b, b + r),
                RowSense::Equal => (b + r.min(0.0), b),
            };
            row_lower.push(lo);
            row_upper.push(hi);
        }
        let maximize = m.obj_sense == ObjSense::Maximize;
        let cost = m
            .objective
            .iter()
            .map(|&c| if maximize { -c } else { c })
            .collect();
        Self {
            num_col: m.num_variables,
            num_row: m.num_constraints,
            matrix: m.matrix.clone(),
            row_lower,
            row_upper,
            cost,
            maximize,
            col_lower: m.bounds.iter().map(|b| b.lower).collect(),
            col_upper: m.bounds.iter().map(|b| b.upper).collect(),
            integer: m.bounds.iter().map(|b| b.integer).collect(),
        }
    }

    /// Internal (minimization) objective value of `x`.
    pub fn objective(&self, x: &[f64]) -> f64 {
        self.cost.iter().zip(x).map(|(c, v)| c * v).sum()
    }

    /// Maps an internal objective (or bound) value to the caller's sense.
    pub fn caller_objective(&self, internal: f64) -> f64 {
        if self.maximize { -internal } else { internal }
    }

    /// Row activities `Ax`.
    pub fn activities(&self, x: &[f64]) -> Vec<f64> {
        let mut ax = vec![0.0; self.num_row];
        self.matrix.matvec(x, &mut ax);
        ax
    }

    /// Worst absolute bound, row and integrality violation of `x`.
    pub fn violation(&self, x: &[f64]) -> Violation {
        let mut v = Violation {
            bound: 0.0,
            row: 0.0,
            integrality: 0.0,
        };
        for (j, &xj) in x.iter().enumerate() {
            v.bound = v
                .bound
                .max(self.col_lower[j] - xj)
                .max(xj - self.col_upper[j]);
            if self.integer[j] {
                v.integrality = v.integrality.max((xj - xj.round()).abs());
            }
        }
        for (i, &a) in self.activities(x).iter().enumerate() {
            v.row = v.row.max(self.row_lower[i] - a).max(a - self.row_upper[i]);
        }
        v
    }

    /// Same constraints, zero objective: the feasibility question asked when
    /// the LP relaxation is unbounded (see [`crate::bnb`]).
    pub fn feasibility_version(&self) -> Self {
        Self {
            cost: vec![0.0; self.num_col],
            ..self.clone()
        }
    }

    /// The same problem with columns and rows reordered: column `k` of the
    /// result is column `cols[k]` of `self`, row `i` is row `rows[i]`. Both
    /// slices must be permutations. Racing searches permuted copies
    /// ([`crate::race`]).
    pub fn permuted(&self, cols: &[usize], rows: &[usize]) -> Self {
        debug_assert_eq!(cols.len(), self.num_col);
        debug_assert_eq!(rows.len(), self.num_row);
        let mut new_row = vec![0; self.num_row];
        for (i, &r) in rows.iter().enumerate() {
            new_row[r] = i;
        }
        let a = &self.matrix;
        let mut b = CscBuilder::new(self.num_row);
        let mut entries = Vec::new();
        for &j in cols {
            entries.clear();
            let range = a.col_ptr()[j]..a.col_ptr()[j + 1];
            for (&r, &v) in a.row_idx()[range.clone()].iter().zip(&a.values()[range]) {
                entries.push((new_row[r as usize], v));
            }
            entries.sort_unstable_by_key(|e| e.0);
            b.push_col(&entries);
        }
        let by_col = |v: &[f64]| cols.iter().map(|&j| v[j]).collect::<Vec<_>>();
        let by_row = |v: &[f64]| rows.iter().map(|&r| v[r]).collect::<Vec<_>>();
        Self {
            num_col: self.num_col,
            num_row: self.num_row,
            matrix: b.into_matrix(),
            row_lower: by_row(&self.row_lower),
            row_upper: by_row(&self.row_upper),
            cost: by_col(&self.cost),
            maximize: self.maximize,
            col_lower: by_col(&self.col_lower),
            col_upper: by_col(&self.col_upper),
            integer: cols.iter().map(|&j| self.integer[j]).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_io::read_mps;

    /// A permuted copy has the same objective and violations at the
    /// correspondingly permuted point.
    #[test]
    fn permuted_preserves_objective_and_violation() {
        let trip = [
            (0, 0, 1.0),
            (0, 2, -2.0),
            (1, 1, 3.0),
            (2, 0, 0.5),
            (2, 2, 1.0),
        ];
        let p = MipProblem {
            num_col: 3,
            num_row: 3,
            matrix: CscMatrix::from_triplets(3, 3, &trip),
            row_lower: vec![-1.0, 0.0, f64::NEG_INFINITY],
            row_upper: vec![1.0, 2.0, 1.5],
            cost: vec![1.0, -2.0, 3.0],
            maximize: false,
            col_lower: vec![0.0, -1.0, 0.0],
            col_upper: vec![2.0, 1.0, 1.0],
            integer: vec![true, false, true],
        };
        let (cols, rows) = ([2, 0, 1], [1, 2, 0]);
        let q = p.permuted(&cols, &rows);
        for x in [[1.0, 0.5, 1.0], [2.0, -1.0, 0.0], [0.3, 0.0, 1.0]] {
            let y: Vec<f64> = cols.iter().map(|&j| x[j]).collect();
            assert_eq!(p.objective(&x), q.objective(&y));
            assert_eq!(p.violation(&x), q.violation(&y));
        }
        assert_eq!(q.integer, [true, true, false]);
        assert_eq!(q.row_upper, [2.0, 1.5, 1.0]);
    }

    #[test]
    fn ranged_rows_and_max_sense() {
        let mps = "\
NAME T
OBJSENSE
 MAX
ROWS
 N  OBJ
 L  R0
 G  R1
 E  R2
 E  R3
COLUMNS
    MARKER  'MARKER'  'INTORG'
    X0  OBJ  2.0  R0  1.0
    X0  R1  1.0  R2  1.0
    MARKER  'MARKER'  'INTEND'
    X1  R3  1.0
RHS
    RHS  R0  4.0  R1  1.0
    RHS  R2  2.0  R3  3.0
RANGES
    RNG  R0  1.5  R1  -2.0
    RNG  R2  -1.0  R3  0.5
ENDATA
";
        let p = MipProblem::from_model(&read_mps(mps).unwrap());
        assert_eq!(p.row_lower, vec![2.5, 1.0, 1.0, 3.0]);
        assert_eq!(p.row_upper, vec![4.0, 3.0, 2.0, 3.5]);
        assert_eq!(p.cost, vec![-2.0, 0.0]);
        assert_eq!(p.integer, vec![true, false]);
        assert_eq!(p.caller_objective(-6.0), 6.0);
    }
}
