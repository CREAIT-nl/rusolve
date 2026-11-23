//! Semi-continuous columns and general
//! constraints survive LP -> Model -> LP and MPS -> Model -> MPS.

use rusolve_io::lp::read_lp;
use rusolve_io::lp_write::write_lp;
use rusolve_io::mps::{read_mps, write_mps};
use rusolve_io::{FuncKind, GcModel, Model, SosKind};

const LP: &str = "\
Minimize
 obj: x + 2 y + z - w
Subject To
 c1: x + y + z >= 2
 ind1: b = 1 -> x + y <= 4
 ind2: b2 = 0 -> y - z >= -1
Bounds
 0 <= x <= 10
 2 <= y <= 8
 0 <= z <= 5
 0 <= w <= 9
 0 <= r <= 20
 0 <= q <= 20
Binaries
 b b2 a1 a2
Generals
 k
Semi-continuous
 y k
SOS
 s1: S1 :: x : 1 z : 2 w : 3
 s2: S2 5 :: x : 1 y : 2 z : 3
General Constraints
 m1: r = MAX ( x , z , 0.5 )
 m2: q = ABS ( w )
 n1: a1 = AND ( b , b2 )
 pw: z = PWL ( x ) : (0, 0) (2, 1) (4, 4)
End
";

fn gc(m: &Model) -> &GcModel {
    &m.gc
}

#[test]
fn lp_reads_every_object() {
    let m = read_lp(LP).unwrap();
    let g = gc(&m);
    assert_eq!(g.indicators.len(), 2);
    assert!(g.indicators[0].bin_val && !g.indicators[1].bin_val);
    assert_eq!(g.semi.len(), 2);
    assert!(g.semi.iter().any(|s| s.integer), "k is semi-integer");
    assert_eq!(g.sos.len(), 2);
    assert_eq!(g.sos[1].kind, SosKind::Sos2);
    assert_eq!(g.sos[1].priority, 5.0);
    assert_eq!(g.funcs.len(), 3);
    assert_eq!(g.funcs[0].kind, FuncKind::Max);
    assert_eq!(g.funcs[0].constant, Some(0.5));
    assert_eq!(g.pwl.len(), 1);
    assert_eq!(g.pwl[0].xs, [0.0, 2.0, 4.0]);
    // indicator rows are not regular rows
    assert_eq!(m.num_constraints, 1);
}

fn same(a: &Model, b: &Model) {
    assert_eq!(a.num_constraints, b.num_constraints);
    assert_eq!(a.num_variables, b.num_variables);
    let (x, y) = (&a.gc, &b.gc);
    assert_eq!(x.semi, y.semi);
    assert_eq!(x.indicators, y.indicators);
    assert_eq!(x.sos, y.sos);
    assert_eq!(x.funcs, y.funcs);
    assert_eq!(x.pwl, y.pwl);
}

#[test]
fn lp_round_trip() {
    let m = read_lp(LP).unwrap();
    let text = write_lp(&m);
    let m2 = read_lp(&text).unwrap_or_else(|e| panic!("{e:?}\n{text}"));
    same(&m, &m2);
    assert_eq!(write_lp(&m2), text, "writer is a fixed point");
}

#[test]
fn mps_round_trip_without_pwl_and_funcs() {
    let mut m = read_lp(LP).unwrap();
    m.gc.pwl.clear();
    m.gc.funcs.clear();
    let text = write_mps(&m);
    let m2 = read_mps(&text).unwrap_or_else(|e| panic!("{e:?}\n{text}"));
    same(&m, &m2);
}

#[test]
fn mps_sc_si_bounds() {
    let text = "\
NAME T
ROWS
 N obj
 G c1
COLUMNS
 x obj 1 c1 1
 k obj 1 c1 1
RHS
 rhs c1 3
BOUNDS
 LO bnd x 2
 SC bnd x 10
 LO bnd k 1
 SI bnd k 7
ENDATA
";
    let m = read_mps(text).unwrap();
    assert_eq!(m.gc.semi.len(), 2);
    assert!(!m.gc.semi[0].integer && m.gc.semi[1].integer);
    assert_eq!(m.bounds[0].upper, 10.0);
    assert_eq!(m.bounds[0].lower, 2.0);
    assert!(!m.bounds[1].integer);
}

#[test]
fn bad_objects_are_typed_errors() {
    let e = read_lp("Minimize\n obj: x\nSubject To\n c: x >= 0\nSOS\n s: S1 :: x\nEnd\n");
    assert!(e.is_ok(), "single-member SOS1 is legal");
    let e = read_lp(
        "Minimize\n obj: x\nSubject To\n c: x >= 0\nGeneral Constraints\n g: y = ABS ( x , x )\nEnd\n",
    );
    assert!(e.unwrap_err().message.contains("abs"));
}

/// Removing an indicator row renumbers the rows after it; a QCMATRIX on a later row must follow.
#[test]
fn mps_qcmatrix_row_is_renumbered_after_an_indicator_row_is_removed() {
    let text = "NAME          T
ROWS
 N  obj
 L  ind1
 L  c1
 L  q1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    z  obj  -1  ind1  0
    MARKER  'MARKER'  'INTEND'
    x  obj  -1  ind1  1  c1  1
    y  obj  -1  q1  1
RHS
    rhs  ind1  1  c1  5  q1  4
BOUNDS
 UP bnd z 1
QCMATRIX   q1
    x  x  1
    y  y  2
INDICATORS
 IF ind1  z  1
ENDATA
";
    let m = rusolve_io::mps::read_mps_conic(text).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(m.row_names, ["c1", "q1"]);
    assert_eq!(m.gc.indicators.len(), 1);
    assert_eq!(m.gc.indicators[0].name, "ind1");
    assert_eq!(m.conic.quad_rows.len(), 1);
    assert_eq!(m.row_names[m.conic.quad_rows[0].row], "q1");
}
