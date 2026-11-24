//! The translator and the C API, end to end, on tiny models. The tests that call CP-SAT need the
//! OR-Tools library at run time (`RUSOLVE_ORTOOLS_LIB=<path of libortools.so.9>`, or any other
//! place `rusolve_cpsat::runtime` searches) and skip, saying so, when it is not found.

use rusolve_cpsat::{CpError, CpParams, CpStatus, ObjectiveKind, Refusal, solve, translate};
use rusolve_ekk::prelude::{HighsInt, HighsLp, HighsVarType, MatrixFormat, ObjSense};
use rusolve_io::cp::*;

/// Skips (returns from) a test that calls CP-SAT when the library is not installed.
macro_rules! need_lib {
    () => {
        if rusolve_cpsat::runtime::library().is_err() {
            eprintln!("skipped: the OR-Tools library was not found (set RUSOLVE_ORTOOLS_LIB)");
            return;
        }
    };
}

type Row = (Vec<(usize, f64)>, f64, f64);

fn lp(cols: &[(f64, f64)], cost: &[f64], maximize: bool, rows: &[Row]) -> HighsLp {
    let n = cols.len();
    let mut l = HighsLp {
        num_col_: n as HighsInt,
        num_row_: rows.len() as HighsInt,
        col_cost_: cost.to_vec(),
        col_lower_: cols.iter().map(|c| c.0).collect(),
        col_upper_: cols.iter().map(|c| c.1).collect(),
        row_lower_: rows.iter().map(|r| r.1).collect(),
        row_upper_: rows.iter().map(|r| r.2).collect(),
        sense_: if maximize {
            ObjSense::kMaximize
        } else {
            ObjSense::kMinimize
        },
        integrality_: vec![HighsVarType::kInteger; n],
        ..Default::default()
    };
    let (mut start, mut index, mut value) = (vec![0 as HighsInt], vec![], vec![]);
    for j in 0..n {
        for (i, r) in rows.iter().enumerate() {
            for &(c, v) in &r.0 {
                if c == j {
                    index.push(i as HighsInt);
                    value.push(v);
                }
            }
        }
        start.push(index.len() as HighsInt);
    }
    l.a_matrix_.format_ = MatrixFormat::kColwise;
    l.a_matrix_.num_col_ = n as HighsInt;
    l.a_matrix_.num_row_ = rows.len() as HighsInt;
    l.a_matrix_.start_ = start;
    l.a_matrix_.index_ = index;
    l.a_matrix_.value_ = value;
    l
}

fn params() -> CpParams {
    CpParams {
        time_limit: Some(30.0),
        ..CpParams::default()
    }
}

fn optimal(lp: &HighsLp, cp: Option<&CpModel>) -> (f64, Vec<f64>) {
    let o = solve(lp, cp, &params()).expect("solves");
    assert_eq!(o.status, CpStatus::Optimal);
    (o.objective.expect("objective"), o.x.expect("point"))
}

fn iv(name: &str, start: usize, size: CpTerm, end: usize, presence: Option<usize>) -> Interval {
    Interval {
        name: name.into(),
        start,
        size,
        end,
        presence,
    }
}

#[test]
fn small_integer_program() {
    need_lib!();
    // max 5x + 4y  s.t.  6x + 4y <= 24,  x + 2y <= 6
    let m = lp(
        &[(0.0, 10.0), (0.0, 10.0)],
        &[5.0, 4.0],
        true,
        &[
            (vec![(0, 6.0), (1, 4.0)], f64::NEG_INFINITY, 24.0),
            (vec![(0, 1.0), (1, 2.0)], f64::NEG_INFINITY, 6.0),
        ],
    );
    let (obj, x) = optimal(&m, None);
    assert_eq!(obj, 20.0);
    assert_eq!(x, vec![4.0, 0.0]);
}

#[test]
fn objective_offset_and_sense() {
    need_lib!();
    let mut m = lp(&[(0.0, 4.0)], &[1.0], true, &[]);
    m.offset_ = 7.0;
    assert_eq!(optimal(&m, None).0, 11.0);
    m.sense_ = ObjSense::kMinimize;
    assert_eq!(optimal(&m, None).0, 7.0);
}

#[test]
fn fractional_cost_goes_floating_point() {
    need_lib!();
    let m = lp(&[(3.0, 9.0)], &[0.5], false, &[]);
    let t = translate(&m, None).unwrap();
    assert_eq!(t.objective_kind, ObjectiveKind::FloatingPoint);
    let o = solve(&m, None, &params()).unwrap();
    assert_eq!(o.objective_kind, ObjectiveKind::FloatingPoint);
    assert!((o.objective.unwrap() - 1.5).abs() < 1e-9);
}

#[test]
fn infeasible_is_proved() {
    need_lib!();
    let m = lp(&[(0.0, 5.0)], &[1.0], false, &[(vec![(0, 1.0)], 7.0, 9.0)]);
    let o = solve(&m, None, &params()).unwrap();
    assert_eq!(o.status, CpStatus::Infeasible);
    assert!(o.x.is_none());
}

#[test]
fn refusals_name_the_first_offender() {
    let mut m = lp(&[(0.0, 5.0), (0.0, 5.0)], &[1.0, 1.0], false, &[]);
    m.integrality_[1] = HighsVarType::kContinuous;
    m.col_names_ = vec!["a".into(), "b".into()];
    assert_eq!(
        translate(&m, None).unwrap_err(),
        Refusal::ContinuousColumn {
            col: 1,
            name: "b".into()
        }
    );
    let m = lp(&[(0.0, f64::INFINITY)], &[1.0], false, &[]);
    assert!(matches!(
        translate(&m, None).unwrap_err(),
        Refusal::InfiniteColumnBound {
            col: 0,
            side: "upper",
            ..
        }
    ));
    let m = lp(&[(0.0, 5.0)], &[1.0], false, &[(vec![(0, 0.5)], 0.0, 1.0)]);
    assert!(matches!(
        translate(&m, None).unwrap_err(),
        Refusal::NonIntegralCoefficient { row: 0, .. }
    ));
    let m = lp(&[(0.0, 5.0)], &[1.0], false, &[(vec![(0, 1.0)], 0.5, 1.0)]);
    assert!(matches!(
        translate(&m, None).unwrap_err(),
        Refusal::NonIntegralRowBound {
            row: 0,
            side: "lower",
            ..
        }
    ));
    // through `solve` the refusal is a CpError
    assert!(matches!(
        solve(&m, None, &params()),
        Err(CpError::Refused(_))
    ));
}

#[test]
fn no_overlap_two_jobs() {
    need_lib!();
    // columns: s1 e1 s2 e2 makespan; sizes 3 and 2
    let cols = vec![(0.0, 20.0); 5];
    let m = lp(
        &cols,
        &[0.0, 0.0, 0.0, 0.0, 1.0],
        false,
        &[
            (vec![(4, 1.0), (1, -1.0)], 0.0, f64::INFINITY),
            (vec![(4, 1.0), (3, -1.0)], 0.0, f64::INFINITY),
        ],
    );
    let cp = CpModel {
        intervals: vec![
            iv("a", 0, CpTerm::Const(3), 1, None),
            iv("b", 2, CpTerm::Const(2), 3, None),
        ],
        no_overlap: vec![NoOverlap {
            name: "m".into(),
            intervals: vec![0, 1],
        }],
        ..Default::default()
    };
    assert_eq!(optimal(&m, Some(&cp)).0, 5.0);
    // without the no-overlap both run at once
    let free = CpModel {
        no_overlap: vec![],
        ..cp
    };
    assert_eq!(optimal(&m, Some(&free)).0, 3.0);
}

#[test]
fn no_overlap_zero_size_overlaps_nothing() {
    need_lib!();
    // B = [0, 10] fixed; A starts at 5 with a size in [0, 5], maximize the size: our definition
    // lets a size-0 interval sit inside B, so the answer is 0 (plain CP-SAT would be infeasible)
    let m = lp(
        &[
            (0.0, 0.0),
            (10.0, 10.0),
            (5.0, 5.0),
            (0.0, 30.0),
            (0.0, 5.0),
        ],
        &[0.0, 0.0, 0.0, 0.0, 1.0],
        true,
        &[],
    );
    let cp = CpModel {
        intervals: vec![
            iv("B", 0, CpTerm::Const(10), 1, None),
            iv("A", 2, CpTerm::Col(4), 3, None),
        ],
        no_overlap: vec![NoOverlap {
            name: "m".into(),
            intervals: vec![0, 1],
        }],
        ..Default::default()
    };
    let o = solve(&m, Some(&cp), &params()).unwrap();
    assert_eq!(o.status, CpStatus::Optimal);
    assert_eq!(o.objective, Some(0.0));
    assert!(o.aux_variables > 0);
}

#[test]
fn optional_interval_absent() {
    need_lib!();
    // A = [0, 10] fixed; B (size 5, start in [0, 3]) is optional with presence column p;
    // maximize p: B would overlap A, so p = 0
    let m = lp(
        &[
            (0.0, 0.0),
            (10.0, 10.0),
            (0.0, 3.0),
            (0.0, 30.0),
            (0.0, 1.0),
        ],
        &[0.0, 0.0, 0.0, 0.0, 1.0],
        true,
        &[],
    );
    let cp = CpModel {
        intervals: vec![
            iv("A", 0, CpTerm::Const(10), 1, None),
            iv("B", 2, CpTerm::Const(5), 3, Some(4)),
        ],
        no_overlap: vec![NoOverlap {
            name: "m".into(),
            intervals: vec![0, 1],
        }],
        ..Default::default()
    };
    let (obj, x) = optimal(&m, Some(&cp));
    assert_eq!(obj, 0.0);
    assert_eq!(x[4], 0.0);
    // a presence column that is not binary is refused
    let bad = lp(
        &[
            (0.0, 0.0),
            (10.0, 10.0),
            (0.0, 3.0),
            (0.0, 30.0),
            (0.0, 4.0),
        ],
        &[0.0; 5],
        false,
        &[],
    );
    assert!(matches!(
        translate(&bad, Some(&cp)).unwrap_err(),
        Refusal::Cp { .. }
    ));
}

#[test]
fn cumulative_capacity_two() {
    need_lib!();
    // three jobs of size 2, demand 1, capacity 2: two in parallel, then one: makespan 4.
    // columns: s1 e1 s2 e2 s3 e3 makespan
    let m = lp(
        &[(0.0, 20.0); 7],
        &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        false,
        &[
            (vec![(6, 1.0), (1, -1.0)], 0.0, f64::INFINITY),
            (vec![(6, 1.0), (3, -1.0)], 0.0, f64::INFINITY),
            (vec![(6, 1.0), (5, -1.0)], 0.0, f64::INFINITY),
        ],
    );
    let cp = CpModel {
        intervals: (0..3)
            .map(|k| iv(&format!("j{k}"), 2 * k, CpTerm::Const(2), 2 * k + 1, None))
            .collect(),
        cumulative: vec![Cumulative {
            name: "r".into(),
            intervals: vec![0, 1, 2],
            demands: vec![CpTerm::Const(1); 3],
            capacity: CpTerm::Const(2),
        }],
        ..Default::default()
    };
    assert_eq!(optimal(&m, Some(&cp)).0, 4.0);
}

#[test]
fn all_different_permutation() {
    need_lib!();
    // min x + 2y + 3z over distinct values in [1, 3]: z = 1, y = 2, x = 3: 3 + 4 + 3 = 10
    let m = lp(&[(1.0, 3.0); 3], &[1.0, 2.0, 3.0], false, &[]);
    let cp = CpModel {
        all_different: vec![AllDifferent {
            name: "d".into(),
            cols: vec![0, 1, 2],
        }],
        ..Default::default()
    };
    let (obj, x) = optimal(&m, Some(&cp));
    assert_eq!(obj, 10.0);
    assert_eq!(x, vec![3.0, 2.0, 1.0]);
}

#[test]
fn element_picks_by_index() {
    need_lib!();
    // target = values[index], the index is 0-based and limited to the list: max target = 30
    let m = lp(&[(0.0, 5.0), (0.0, 100.0)], &[0.0, 1.0], true, &[]);
    let cp = CpModel {
        element: vec![Element {
            name: "e".into(),
            index: 0,
            values: vec![CpTerm::Const(10), CpTerm::Const(20), CpTerm::Const(30)],
            target: 1,
        }],
        ..Default::default()
    };
    let (obj, x) = optimal(&m, Some(&cp));
    assert_eq!(obj, 30.0);
    assert_eq!(x[0], 2.0);
}

#[test]
fn table_allowed_and_forbidden_tuples() {
    need_lib!();
    let m = lp(&[(0.0, 3.0), (0.0, 3.0)], &[1.0, 1.0], true, &[]);
    let allowed = CpModel {
        table: vec![Table {
            name: "t".into(),
            cols: vec![0, 1],
            tuples: vec![vec![0, 1], vec![2, 3], vec![1, 1]],
            negated: false,
        }],
        ..Default::default()
    };
    assert_eq!(optimal(&m, Some(&allowed)).0, 5.0);
    let forbidden = CpModel {
        table: vec![Table {
            name: "t".into(),
            cols: vec![0, 1],
            tuples: vec![vec![3, 3], vec![2, 3], vec![3, 2]],
            negated: true,
        }],
        ..Default::default()
    };
    assert_eq!(optimal(&m, Some(&forbidden)).0, 4.0);
    // a ragged tuple is a typed refusal
    let ragged = CpModel {
        table: vec![Table {
            name: "t".into(),
            cols: vec![0, 1],
            tuples: vec![vec![0]],
            negated: false,
        }],
        ..Default::default()
    };
    assert!(matches!(
        translate(&m, Some(&ragged)).unwrap_err(),
        Refusal::Cp { .. }
    ));
}

#[test]
fn bad_indices_are_refused() {
    let m = lp(&[(0.0, 3.0)], &[1.0], false, &[]);
    let cp = CpModel {
        all_different: vec![AllDifferent {
            name: "d".into(),
            cols: vec![0, 7],
        }],
        ..Default::default()
    };
    assert!(matches!(
        translate(&m, Some(&cp)).unwrap_err(),
        Refusal::Cp { .. }
    ));
    let cp = CpModel {
        no_overlap: vec![NoOverlap {
            name: "n".into(),
            intervals: vec![0],
        }],
        ..Default::default()
    };
    assert!(matches!(
        translate(&m, Some(&cp)).unwrap_err(),
        Refusal::Cp { .. }
    ));
}

#[test]
fn a_missing_library_is_a_typed_error_with_the_search_list() {
    // a library variable that names nothing: the engine says so and lists where it looked. This
    // test needs a process that has not loaded the library, so it checks the pure parts.
    let c = rusolve_cpsat::runtime::candidates_for(
        Some("/nonexistent/libortools.so.9".into()),
        None,
        None,
    );
    assert_eq!(c.len(), 1);
    let e = rusolve_cpsat::runtime::LoadError {
        searched: c.iter().map(|c| c.to_string()).collect(),
        problems: vec![],
    };
    let text = e.to_string();
    assert!(text.contains("/nonexistent/libortools.so.9"), "{text}");
    assert!(text.contains("install.sh --with-cpsat"), "{text}");
    assert!(text.contains("pip install 'rusolve[cpsat]'"), "{text}");
}
