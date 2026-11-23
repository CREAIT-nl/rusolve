//! Quadratic constraints and cones, stored as the user wrote them
//! (the original row is kept for the checker; no reformulation lives here).
//!
//! ## Conventions
//!
//! * A [`QuadRow`] adds `½ xᵀQx` to the activity of constraint row `row`, the
//!   same convention as [`Model::hessian`](crate::Model::hessian) for the
//!   objective. `Q` is symmetric and stored as its lower triangle
//!   (`i >= j`, each off-diagonal pair once): the row activity is
//!   `aᵀx + ½ Σ_i Q_ii x_i² + Σ_{i>j} Q_ij x_i x_j`.
//! * File formats write the *monomial* coefficients without the half (MPS
//!   `QCMATRIX`, LP `[ ... ]` in a constraint). [`QuadRow::from_monomials`]
//!   converts: a monomial `m x_i²` is `Q_ii = 2m`, and `m x_i x_j` (i != j,
//!   from either triangle, summed) is `Q_ij = m`.
//! * A [`Cone`] constrains the listed *variables* (in order) to the cone.
//!   `Soc`: `x0 >= ‖x1..‖`. `RotatedSoc`: `2 x0 x1 >= ‖x2..‖²`, `x0, x1 >= 0`.
//!   `Exp` and `Pow` are reserved in the type and rejected by the readers and
//!   the router until a solver path implements them.
//!
//! Provenance: MPS `QCMATRIX`/`CSECTION` follow the public MOSEK, CPLEX and
//! Gurobi file-format documentation; the `QCMATRIX` scale (monomial
//! coefficients, no `½`) matches the file SCIP 10 writes for
//! `2x² + 3xy + y² + x <= 5`.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};

use rusolve_core::Triplet;

/// The quadratic part of one constraint row: `½ xᵀQx`, `Q` as its lower
/// triangle `(row, col, value)` with `row >= col`, sorted, no zeros.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadRow {
    /// Index into the model's constraint rows.
    pub row: usize,
    pub q: Vec<Triplet>,
}

impl QuadRow {
    /// Builds the row from monomial coefficients `m · x_i x_j`, in any
    /// triangle, duplicates summed (see the module doc for the scale).
    pub fn from_monomials(row: usize, monomials: impl IntoIterator<Item = Triplet>) -> QuadRow {
        let mut acc: BTreeMap<(usize, usize), f64> = BTreeMap::new();
        for (i, j, v) in monomials {
            *acc.entry((i.max(j), i.min(j))).or_insert(0.0) += v;
        }
        let q = acc
            .into_iter()
            .filter(|&(_, v)| v != 0.0)
            .map(|((i, j), v)| (i, j, if i == j { 2.0 * v } else { v }))
            .collect();
        QuadRow { row, q }
    }

    /// The monomial coefficients `(i, j, m)`, `i >= j`: the inverse of
    /// [`QuadRow::from_monomials`], what file writers emit.
    pub fn monomials(&self) -> impl Iterator<Item = Triplet> + '_ {
        self.q
            .iter()
            .map(|&(i, j, v)| (i, j, if i == j { v / 2.0 } else { v }))
    }

    /// `½ xᵀQx` at `x`.
    pub fn value(&self, x: &[f64]) -> f64 {
        self.monomials().map(|(i, j, m)| m * x[i] * x[j]).sum()
    }
}

/// The kind of a cone. Only `Soc` and `RotatedSoc` are implemented.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConeKind {
    /// `x0 >= ‖x1..‖₂`.
    Soc,
    /// `2 x0 x1 >= ‖x2..‖₂²`, `x0, x1 >= 0`.
    RotatedSoc,
    /// Reserved: primal exponential cone (3 variables). Unimplemented.
    Exp,
    /// Reserved: power cone with exponent `alpha`. Unimplemented.
    Pow(f64),
}

impl ConeKind {
    pub fn label(self) -> &'static str {
        match self {
            ConeKind::Soc => "SOC",
            ConeKind::RotatedSoc => "rotated SOC",
            ConeKind::Exp => "exponential",
            ConeKind::Pow(_) => "power",
        }
    }
}

/// One cone over model columns.
#[derive(Debug, Clone, PartialEq)]
pub struct Cone {
    pub name: String,
    pub kind: ConeKind,
    /// Column indices, in cone order.
    pub vars: Vec<usize>,
}

/// Everything a [`Model`](crate::Model) carries beyond LP, MIP and the
/// objective Hessian. Empty for every model that takes today's code path.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ConicBlock {
    pub quad_rows: Vec<QuadRow>,
    pub cones: Vec<Cone>,
}

impl ConicBlock {
    pub fn is_empty(&self) -> bool {
        self.quad_rows.is_empty() && self.cones.is_empty()
    }
}

/// The message the LP/MIP readers give when a file has quadratic rows or
/// cones, so the caller can tell it from a parse error.
pub const CONIC_CONTENT_MSG: &str =
    "model has quadratic constraints or cones; it takes the conic route (rusolve-cli), not LP/MIP";

static CONIC_REJECTED: AtomicBool = AtomicBool::new(false);

/// Records that a legacy (LP/MIP) reader rejected a conic model. The CLI
/// reads it after a failed solve to re-dispatch to the conic route, which
/// keeps the path of every non-conic model exactly as it was.
pub fn note_conic_rejection() {
    CONIC_REJECTED.store(true, Ordering::Relaxed);
}

/// Whether a legacy reader has rejected a conic model in this process.
pub fn conic_rejected() -> bool {
    CONIC_REJECTED.load(Ordering::Relaxed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monomials_round_trip_and_value() {
        // 2x² + 3xy + y²  (x = 0, y = 1)
        let r = QuadRow::from_monomials(0, [(0, 0, 2.0), (0, 1, 1.5), (1, 0, 1.5), (1, 1, 1.0)]);
        assert_eq!(r.q, vec![(0, 0, 4.0), (1, 0, 3.0), (1, 1, 2.0)]);
        let m: Vec<_> = r.monomials().collect();
        assert_eq!(m, vec![(0, 0, 2.0), (1, 0, 3.0), (1, 1, 1.0)]);
        assert_eq!(r.value(&[1.0, 2.0]), 2.0 + 6.0 + 4.0);
    }
}
