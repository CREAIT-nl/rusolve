//! MPS → `SparseLp` load round-trip.
//!
//! Integration test (not a `src/` unit test) because rusolve-io depends on
//! rusolve-core: a dev-dependency cycle typechecks only from the integration
//! harness, where the lib is linked once. This is the load seam
//! `MipProblem::from_model` will call once rusolve-mip moves to `SparseLp`
//! (sense/rhs → row bounds; the matrix is already CSC) — rusolve-io and
//! rusolve-mip themselves are not changed here.

use rusolve_core::SparseLp;
use rusolve_io::{ObjSense, RowSense, read_mps};

fn bits_eq(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(p, q)| p.to_bits() == q.to_bits())
}

// Mixed senses, negative rhs/values, free / flipped / fixed columns.
// No RANGES: the seam maps r = 0 exactly as MipProblem::from_model does.
const MPS: &str = "\
NAME          B1S1
OBJSENSE
  MAX
ROWS
  N  cost
  L  cap
  G  demand
  E  balance
  L  negrhs
COLUMNS
  xa        cost        2.5
  xa        cap         1.0
  xa        demand      -3.0
  xa        balance     1.0
  xb        cost        -1.5
  xb        cap         0.5
  xb        negrhs      -2.0
  xc        balance     1.0
  xc        cap         -0.25
  xd        demand      4.0
  xd        cost        0.75
RHS
  r0        cap         10.0
  r0        demand      -2.0
  r0        balance     3.5
  r0        negrhs      -7.0
BOUNDS
  FR bnd xa
  MI bnd xb
  UP bnd xb 6.0
  LO bnd xc 1.0
  FX bnd xd 2.0
ENDATA
";

#[test]
fn mps_load_round_trip_bit_identical() {
    let model = read_mps(MPS).expect("reader regression?");
    assert_eq!((model.num_variables, model.num_constraints), (4, 4));

    let inf = f64::INFINITY;
    let (row_lower, row_upper): (Vec<f64>, Vec<f64>) = model
        .senses
        .iter()
        .zip(&model.rhs)
        .map(|(s, &b)| match s {
            RowSense::LessEqual => (-inf, b),
            RowSense::GreaterEqual => (b, inf),
            RowSense::Equal => (b, b),
        })
        .unzip();
    let col_lower: Vec<f64> = model.bounds.iter().map(|b| b.lower).collect();
    let col_upper: Vec<f64> = model.bounds.iter().map(|b| b.upper).collect();
    let sp = SparseLp::new(
        model.matrix.clone(),
        row_lower.clone(),
        row_upper.clone(),
        col_lower.clone(),
        col_upper.clone(),
        model.objective.clone(),
        model.obj_sense == ObjSense::Maximize,
    );

    // Matrix arrives CSC from the reader and is kept verbatim (MPS entries
    // are nonzero ⇒ I3 drops nothing).
    assert_eq!(sp.a, model.matrix);
    assert_eq!(sp.a.nnz(), model.matrix.nnz());
    assert!(sp.maximize);

    let dense = sp.to_lp();
    assert_eq!((dense.num_col, dense.num_row), (4, 4));
    assert_eq!(dense.senses, vec!['L', 'G', 'E', 'L']);
    assert!(bits_eq(&dense.rhs, &model.rhs));
    assert!(bits_eq(&dense.objective, &model.objective));
    assert!(bits_eq(&dense.lower, &col_lower));
    assert!(bits_eq(&dense.upper, &col_upper));
    // The solver-space CSC agrees with the reader's CSC entry-for-entry,
    // bit for bit (structural zeros included, +0.0 on both sides).
    for j in 0..4 {
        for i in 0..4 {
            assert_eq!(
                dense.a.get(i, j).to_bits(),
                model.matrix.get(i, j).to_bits(),
                "entry ({i}, {j})"
            );
        }
    }
    // And the solver-space form survives a second round-trip unchanged.
    let twice = SparseLp::from_lp(&dense).to_lp();
    assert!(bits_eq(dense.a.values(), twice.a.values()));
    assert!(bits_eq(&twice.rhs, &dense.rhs));
    assert_eq!(twice.senses, dense.senses);
    assert_eq!(twice.maximize, dense.maximize);
}
