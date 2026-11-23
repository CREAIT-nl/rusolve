//! `read_lp(write_lp(m)) == m`, bit for bit, on the hand-written edge-case
//! models of `tests/cplex_lp/models/` and the small netlib instances.
//!
//! The LP text writes a ranged row as `name: expr - Rg<name> = rhs` with a
//! bounded range column (CPLEX LP has no ranged-row syntax; see
//! `src/lp_write.rs`), so the read-back model has one extra column per ranged
//! row. The check folds them: each must sit after the model columns, carry
//! the single entry -1 in its own row and the bounds that give the row its
//! MPS interval.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use rusolve_io::{Model, RowSense, lp_names, read_lp, read_mps, write_lp_report};

fn rg_bounds(sense: RowSense, r: f64) -> (f64, f64) {
    match sense {
        RowSense::LessEqual => (-r.abs(), 0.0),
        RowSense::GreaterEqual => (0.0, r.abs()),
        RowSense::Equal if r > 0.0 => (0.0, r),
        RowSense::Equal => (r, 0.0),
    }
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
}

/// Bit-equal, or within `tol` relative when the writer reported values it
/// wrote in CPLEX's (not strtod's) exact form.
fn close(a: f64, b: f64, tol: f64) -> bool {
    same(a, b) || (a - b).abs() <= tol * a.abs().max(b.abs())
}

fn columns_close(a: &[(usize, f64)], b: &[(usize, f64)], tol: f64) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(x, y)| x.0 == y.0 && close(x.1, y.1, tol))
}

fn column(m: &Model, j: usize) -> Vec<(usize, f64)> {
    let a = &m.matrix;
    (a.col_ptr()[j]..a.col_ptr()[j + 1])
        .map(|k| (a.row_idx()[k] as usize, a.values()[k]))
        .collect()
}

fn hessian_upper(m: &Model) -> BTreeMap<(usize, usize), f64> {
    let mut q = BTreeMap::new();
    for &(r, c, v) in &m.hessian {
        let v = if r == c { v } else { v / 2.0 };
        *q.entry((r.min(c), r.max(c))).or_insert(0.0) += v;
    }
    q.retain(|_, v| *v != 0.0);
    q
}

/// Every difference between `m` and `back = read_lp(write_lp(m))`.
fn diff(m: &Model, back: &Model, tol: f64) -> Vec<String> {
    let names = lp_names(m);
    let mut d = Vec::new();
    let n = m.num_variables;
    let ranged: Vec<usize> = (0..m.num_constraints)
        .filter(|&i| m.ranges[i] != 0.0)
        .collect();
    if back.num_variables != n + ranged.len() {
        d.push(format!(
            "columns: {} read back, expected {} + {} range columns",
            back.num_variables,
            n,
            ranged.len()
        ));
        return d;
    }
    if back.num_constraints != m.num_constraints {
        d.push(format!(
            "rows: {} vs {}",
            back.num_constraints, m.num_constraints
        ));
        return d;
    }
    if back.col_names[..n] != names.columns[..] {
        d.push("column names differ".into());
    }
    if back.row_names != names.rows {
        d.push("row names differ".into());
    }
    if back.obj_name != names.objective {
        d.push(format!(
            "objective name {:?} vs {:?}",
            back.obj_name, names.objective
        ));
    }
    if back.obj_sense != m.obj_sense {
        d.push("objective sense".into());
    }
    if back.obj_offset != m.obj_offset {
        d.push(format!("offset {} vs {}", back.obj_offset, m.obj_offset));
    }
    for j in 0..n {
        if !close(back.objective[j], m.objective[j], tol) && back.objective[j] != m.objective[j] {
            d.push(format!(
                "cost {j}: {:?} vs {:?}",
                back.objective[j], m.objective[j]
            ));
        }
        let (b, o) = (back.bounds[j], m.bounds[j]);
        if !same(b.lower, o.lower) || !same(b.upper, o.upper) || b.integer != o.integer {
            d.push(format!("bounds {j}: {b:?} vs {o:?}"));
        }
        if !columns_close(&column(back, j), &column(m, j), tol) {
            d.push(format!("column {j} entries differ"));
        }
    }
    for (k, &i) in ranged.iter().enumerate() {
        let j = n + k;
        if Some(&back.col_names[j]) != names.range_columns[i].as_ref() {
            d.push(format!("range column {j} name {}", back.col_names[j]));
        }
        if column(back, j) != [(i, -1.0)] || back.objective[j] != 0.0 {
            d.push(format!("range column {j} is not -1 in row {i} only"));
        }
        let (lo, hi) = rg_bounds(m.senses[i], m.ranges[i]);
        let b = back.bounds[j];
        if !same(b.lower, lo) || !same(b.upper, hi) || b.integer {
            d.push(format!(
                "range column {j} bounds {b:?}, expected [{lo}, {hi}]"
            ));
        }
    }
    for i in 0..m.num_constraints {
        let sense = if m.ranges[i] != 0.0 {
            RowSense::Equal
        } else {
            m.senses[i]
        };
        if back.senses[i] != sense || !same(back.rhs[i], m.rhs[i]) || back.ranges[i] != 0.0 {
            d.push(format!(
                "row {i}: {:?} {:?} {:?} vs {:?} {:?} {:?}",
                back.senses[i], back.rhs[i], back.ranges[i], m.senses[i], m.rhs[i], m.ranges[i]
            ));
        }
    }
    if hessian_upper(back) != hessian_upper(m) {
        d.push("hessian differs".into());
    }
    d
}

fn check(path: &Path) -> Result<(), String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let m = read_mps(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let (lp, report) = write_lp_report(&m);
    let back = read_lp(&lp).map_err(|e| format!("{}: read back: {e}", path.display()))?;
    // Values written for CPLEX's reader may sit a few ulps off for a
    // correctly rounded one (only beyond about 1e+-90).
    let tol = if report.cplex_only_numbers > 0 {
        1e-15
    } else {
        0.0
    };
    let d = diff(&m, &back, tol);
    if d.is_empty() {
        Ok(())
    } else {
        Err(format!("{}:\n  {}", path.display(), d.join("\n  ")))
    }
}

fn mps_in(dir: &Path, max_bytes: u64) -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .map(|r| r.filter_map(|e| e.ok().map(|e| e.path())).collect())
        .unwrap_or_default();
    v.retain(|p| {
        p.extension().is_some_and(|e| e == "mps")
            && p.metadata().is_ok_and(|md| md.len() <= max_bytes)
    });
    v.sort();
    v
}

#[test]
fn edge_case_models_round_trip_bit_exactly() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/cplex_lp/models");
    let files = mps_in(&dir, u64::MAX);
    assert!(
        files.len() >= 6,
        "edge-case models missing in {}",
        dir.display()
    );
    let errs: Vec<String> = files.iter().filter_map(|p| check(p).err()).collect();
    assert!(errs.is_empty(), "{}", errs.join("\n"));
}

#[test]
fn netlib_models_round_trip_bit_exactly() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../benches/data/netlib");
    let files = mps_in(&dir, 400_000);
    if files.is_empty() {
        eprintln!("no netlib instances in {}; skipped", dir.display());
        return;
    }
    let errs: Vec<String> = files.iter().filter_map(|p| check(p).err()).collect();
    assert!(errs.is_empty(), "{}", errs.join("\n"));
}

#[test]
fn renamed_and_truncated_names_stay_unique() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/cplex_lp/models");
    let m = read_mps(&std::fs::read_to_string(dir.join("edge_long.mps")).unwrap()).unwrap();
    let names = lp_names(&m);
    assert!(names.columns.iter().all(|c| c.len() <= 254));
    assert!(names.rows.iter().all(|r| r.len() <= 255));
    let mut cols = names.columns.clone();
    cols.sort();
    cols.dedup();
    assert_eq!(cols.len(), m.num_variables);
    let m = read_mps(&std::fs::read_to_string(dir.join("edge_names.mps")).unwrap()).unwrap();
    let names = lp_names(&m);
    // Legal names are kept as they are, including ones that look odd.
    for keep in [
        "x(1,2)",
        "B0(998)",
        "eabc",
        "e_",
        "a_b",
        "a!\"#$%&'?@_`{}~.",
    ] {
        assert!(names.columns.iter().any(|c| c == keep), "{keep} renamed");
    }
    for keep in ["e1", "end", "x(1,2)"] {
        assert!(names.rows.iter().any(|r| r == keep), "row {keep} renamed");
    }
    // `Rge1` is a model column, so the range column of row `e1` moves aside.
    assert_eq!(names.range_columns[0].as_deref(), Some("Rge1~1"));
}

/// A quadratic *objective* typed in LP syntax, with a
/// cross term, survives write and read; a quadratic *constraint* stays an error.
#[test]
fn lp_quadratic_objective_round_trips_and_quadratic_rows_are_refused() {
    let text = "\
Minimize
 obj: - x - 2 y + [ 2 x ^ 2 + 2 x * y + 4 y ^ 2 ] / 2
Subject To
 c1: x + y <= 1
Bounds
 0 <= x <= 10
 0 <= y <= 10
End
";
    let m = read_lp(text).expect("quadratic objective reads");
    // x'Hx/2 = (2x^2 + 2xy + 4y^2)/2  =>  H = [[2, 1], [1, 4]]
    let h = hessian_upper(&m);
    assert_eq!(h.get(&(0, 0)), Some(&2.0));
    assert_eq!(h.get(&(0, 1)), Some(&1.0));
    assert_eq!(h.get(&(1, 1)), Some(&4.0));
    assert_eq!(m.hessian.len(), 4, "both triangles stored");
    let (lp, _) = write_lp_report(&m);
    let back = read_lp(&lp).expect("written LP reads back");
    assert!(
        diff(&m, &back, 0.0).is_empty(),
        "{:?}",
        diff(&m, &back, 0.0)
    );

    let qc = "\
Minimize
 obj: x
Subject To
 q1: [ x ^ 2 + y ^ 2 ] <= 4
End
";
    assert!(
        read_lp(qc).is_err(),
        "quadratic constraints need read_lp_conic"
    );
}
