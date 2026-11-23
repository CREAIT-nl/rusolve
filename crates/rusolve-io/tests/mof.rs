//! MathOptFormat reader against hand-written files that follow the schema
//! (`schemas/mof.1.9.schema.json`, github.com/jump-dev/MathOptFormat; all
//! three validate against it with the `jsonschema` Python package):
//! * `mof/milp.mof.json`: the specification's own MILP example, verbatim
//!   (MIT, Copyright (c) 2020 Oscar Dowson; see THIRD_PARTY_NOTICES.md);
//! * `mof/features.mof.json`: the `tests/nl/features` Pyomo model written by
//!   hand, which must read to the same problem as Pyomo's .nl, LP and MPS
//!   exports of it;
//! * `dialects/handwritten.mof.json`: the dialect-corpus model, checked in
//!   `tests/dialects.rs`.

mod common;

use common::canon;
use rusolve_io::{ObjSense, RowSense, read_model_file};

fn path(f: &str) -> String {
    format!("{}/tests/{f}", env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn specification_milp_example() {
    // min{x | x + y >= 1, x in [0, 1], y in {0, 1}}
    let m = read_model_file(&path("mof/milp.mof.json")).unwrap();
    assert_eq!(m.col_names, ["x", "y"]);
    assert_eq!((m.obj_sense, m.obj_offset), (ObjSense::Minimize, 0.0));
    assert_eq!(m.objective, [1.0, 0.0]);
    assert_eq!(m.num_constraints, 1);
    assert_eq!(m.row_names, ["x + y >= 1"]);
    assert_eq!((m.senses[0], m.rhs[0]), (RowSense::GreaterEqual, 1.0));
    assert_eq!(m.matrix.nnz(), 2);
    let b: Vec<_> = m
        .bounds
        .iter()
        .map(|b| (b.lower, b.upper, b.integer))
        .collect();
    assert_eq!(b, [(0.0, 1.0, false), (0.0, 1.0, true)]);
}

#[test]
fn mof_is_the_same_problem_as_pyomos_exports() {
    let mof = canon(&read_model_file(&path("mof/features.mof.json")).unwrap());
    for other in ["nl/features.nl", "nl/features.lp", "nl/features.mps"] {
        let m = read_model_file(&path(other)).unwrap();
        assert_eq!(canon(&m), mof, "{other}");
    }
}
