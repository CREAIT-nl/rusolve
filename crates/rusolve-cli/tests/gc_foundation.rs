//! Semi-continuous columns and general
//! constraints through the real binary: `rusolve check` planted violations,
//! one per object, and solves with a known optimum.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-check-1.schema.json");

fn tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-gc-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn sol(cols: &[(&str, f64)]) -> String {
    let mut s = format!(
        "Model status\nOptimal\n\n# Primal solution values\nFeasible\nObjective 0\n# Columns {}\n",
        cols.len()
    );
    for (n, v) in cols {
        s.push_str(&format!("{n} {v}\n"));
    }
    s.push_str("# Rows 0\n\n# Dual solution values\nNone\n\n# Basis\nHiGHS_basis_file 2\nNone\n");
    s
}

fn check(tag: &str, model: &str, cols: &[(&str, f64)]) -> Value {
    let m = tmp(&format!("{tag}.lp"), model);
    let s = tmp(&format!("{tag}.sol"), &sol(cols));
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args([
            "check",
            m.to_str().unwrap(),
            s.to_str().unwrap(),
            "--format",
            "json",
        ])
        .output()
        .unwrap();
    assert_eq!(
        out.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    // every kind the report uses is in the schema's enum
    for f in v["findings"].as_array().unwrap() {
        assert!(SCHEMA.contains(f["kind"].as_str().unwrap()));
    }
    v
}

fn kinds(v: &Value) -> Vec<String> {
    v["findings"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| {
            format!(
                "{}:{}",
                f["kind"].as_str().unwrap(),
                f["name"].as_str().unwrap()
            )
        })
        .collect()
}

const SEMI: &str = "Minimize\n obj: 3 x + 2 y\nSubject To\n c1: x + y >= 1\nBounds\n 2 <= x <= 5\n 3 <= y <= 8\nSemi-continuous\n x y\nEnd\n";

#[test]
fn semi_continuous_membership_is_checked() {
    assert_eq!(
        check("semi_ok0", SEMI, &[("x", 2.0), ("y", 0.0)])["verdict"],
        "feasible"
    );
    assert_eq!(
        check("semi_ok1", SEMI, &[("x", 0.0), ("y", 4.5)])["verdict"],
        "feasible"
    );
    let v = check("semi_bad", SEMI, &[("x", 1.0), ("y", 0.0)]);
    assert_eq!(v["verdict"], "infeasible");
    assert!(kinds(&v).contains(&"semi_continuous:x".to_string()), "{v}");
    assert_eq!(v["summary"]["general_constraints"]["count"], 1);
}

#[test]
fn semi_continuous_solve_has_the_known_optimum() {
    // x in {0} U [2,5], y in {0} U [3,8], x + y >= 1: x = 2, y = 0 costs 6
    let m = tmp("semi_solve.lp", SEMI);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["mip", m.to_str().unwrap(), "--format", "json"])
        .output()
        .unwrap();
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["status"], "Optimal");
    assert!((v["objective"].as_f64().unwrap() - 6.0).abs() < 1e-9, "{v}");
}

#[test]
fn semi_integer_solve_and_check() {
    // k semi-integer in {0} U {2,..,6}: k >= 1 forces k = 2
    let text = "Minimize\n obj: k\nSubject To\n c: k >= 1\nBounds\n 2 <= k <= 6\nGenerals\n k\nSemi-continuous\n k\nEnd\n";
    let m = tmp("semi_int.lp", text);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["mip", m.to_str().unwrap(), "--format", "json"])
        .output()
        .unwrap();
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert!((v["objective"].as_f64().unwrap() - 2.0).abs() < 1e-9, "{v}");
    let bad = check("si_bad", text, &[("k", 2.5)]);
    assert!(
        kinds(&bad).contains(&"semi_continuous:k".to_string()),
        "{bad}"
    );
}

const IND: &str = "Minimize\n obj: x\nSubject To\n c: x >= 0\n i1: b = 1 -> x <= 1\nBounds\n 0 <= x <= 9\nBinaries\n b\nEnd\n";

#[test]
fn indicator_is_checked_on_the_rounded_binary() {
    assert_eq!(
        check("ind_ok", IND, &[("x", 5.0), ("b", 0.0)])["verdict"],
        "feasible"
    );
    assert_eq!(
        check("ind_ok1", IND, &[("x", 1.0), ("b", 1.0)])["verdict"],
        "feasible"
    );
    let v = check("ind_bad", IND, &[("x", 5.0), ("b", 1.0)]);
    assert!(kinds(&v).contains(&"indicator:i1".to_string()), "{v}");
    // trickle flow: b = 1e-9 rounds to 0, so the indicator is inactive
    assert_eq!(
        check("ind_trickle", IND, &[("x", 5.0), ("b", 1e-9)])["verdict"],
        "feasible"
    );
    // a fractional binary is itself a violation
    let f = check("ind_frac", IND, &[("x", 5.0), ("b", 0.5)]);
    assert!(kinds(&f).contains(&"indicator:i1".to_string()), "{f}");
}

const SOS: &str = "Minimize\n obj: x + y + z\nSubject To\n c: x + y + z >= 0\nBounds\n 0 <= x <= 5\n 0 <= y <= 5\n 0 <= z <= 5\nSOS\n s1: S1 :: x : 1 y : 2 z : 3\n s2: S2 :: x : 1 y : 2 z : 3\nEnd\n";

#[test]
fn sos_sets_are_checked() {
    // adjacent pair: fine for S2, violates S1
    let v = check("sos_adj", SOS, &[("x", 1.0), ("y", 1.0), ("z", 0.0)]);
    assert_eq!(kinds(&v), ["sos:s1"], "{v}");
    // x and z are not adjacent: violates both
    let v = check("sos_gap", SOS, &[("x", 1.0), ("y", 0.0), ("z", 1.0)]);
    assert_eq!(v["summary"]["general_constraints"]["count"], 2, "{v}");
    assert_eq!(
        check("sos_ok", SOS, &[("x", 0.0), ("y", 2.0), ("z", 0.0)])["verdict"],
        "feasible"
    );
}

const PWL: &str = "Minimize\n obj: y\nSubject To\n c: y >= 0\nBounds\n 0 <= x <= 4\n 0 <= y <= 9\nGeneral Constraints\n p: y = PWL ( x ) : (0, 0) (2, 1) (4, 4)\nEnd\n";

#[test]
fn pwl_values_are_checked() {
    assert_eq!(
        check("pwl_ok", PWL, &[("x", 3.0), ("y", 2.5)])["verdict"],
        "feasible"
    );
    let v = check("pwl_bad", PWL, &[("x", 3.0), ("y", 3.0)]);
    assert!(kinds(&v).contains(&"pwl:p".to_string()), "{v}");
}

const FUNCS: &str = "Minimize\n obj: r + q\nSubject To\n c: r >= -100\nBounds\n -5 <= x <= 5\n -5 <= w <= 5\n r free\n q free\n 0 <= u <= 1\n 0 <= v <= 1\nBinaries\n b1 b2 a o\nGeneral Constraints\n m: r = MAX ( x , w , 0.5 )\n q1: q = ABS ( x )\n a1: a = AND ( b1 , b2 )\n o1: o = OR ( b1 , b2 )\nEnd\n";

#[test]
fn functions_are_evaluated() {
    let ok = [
        ("x", -3.0),
        ("w", 1.0),
        ("r", 1.0),
        ("q", 3.0),
        ("b1", 1.0),
        ("b2", 0.0),
        ("a", 0.0),
        ("o", 1.0),
        ("u", 0.0),
        ("v", 0.0),
    ];
    assert_eq!(check("fn_ok", FUNCS, &ok)["verdict"], "feasible");
    let mut bad = ok;
    bad[2].1 = 0.5; // r is the max of (-3, 1, 0.5) = 1
    bad[6].1 = 1.0; // a = AND(1, 0) = 0
    let v = check("fn_bad", FUNCS, &bad);
    let k = kinds(&v);
    assert!(
        k.contains(&"function:m".to_string()) && k.contains(&"function:a1".to_string()),
        "{v}"
    );
}

#[test]
fn plain_models_have_no_general_constraints_entry() {
    let m = "Minimize\n obj: x\nSubject To\n c: x >= 1\nEnd\n";
    let v = check("plain", m, &[("x", 1.0)]);
    assert!(v["summary"].get("general_constraints").is_none());
}

// ------------------------------------------------------------------ general-constraint solves

/// Solves `model` with `rusolve mip`, checks the written solution with
/// `rusolve check`, and returns (objective, the `general_constraints` report).
fn solve(tag: &str, model: &str, extra: &[&str]) -> (f64, Value) {
    let m = tmp(&format!("{tag}.lp"), model);
    let sf = tmp(&format!("{tag}.out.sol"), "");
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args([
            "mip",
            m.to_str().unwrap(),
            "--format",
            "json",
            "--solution-file",
            sf.to_str().unwrap(),
        ])
        .args(extra)
        .output()
        .unwrap();
    let v: Value = serde_json::from_slice(&out.stdout)
        .unwrap_or_else(|e| panic!("{e}: {}", String::from_utf8_lossy(&out.stderr)));
    assert_eq!(v["status"], "Optimal", "{tag}: {v}");
    let chk = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args([
            "check",
            m.to_str().unwrap(),
            sf.to_str().unwrap(),
            "--format",
            "json",
        ])
        .output()
        .unwrap();
    let c: Value = serde_json::from_slice(&chk.stdout).unwrap();
    assert_eq!(
        c["verdict"], "feasible",
        "{tag}: the solve's own solution fails the checker: {c}"
    );
    (
        v["objective"].as_f64().unwrap(),
        v["general_constraints"].clone(),
    )
}

fn route(g: &Value, object: &str) -> String {
    g["routes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["object"].as_str().unwrap().starts_with(object))
        .unwrap_or_else(|| panic!("no route for {object}: {g}"))["route"]
        .as_str()
        .unwrap()
        .to_string()
}

#[test]
fn indicator_solve_known_optimum() {
    // b = 0 -> x >= 5: either pay x = 5 (cost 5) or b = 1 (cost 3)
    let m = "Minimize\n obj: x + 3 b\nSubject To\n c: x + b >= 0\n i: b = 0 -> x >= 5\nBounds\n 0 <= x <= 10\nBinaries\n b\nEnd\n";
    let (obj, g) = solve("s_ind", m, &[]);
    assert!((obj - 3.0).abs() < 1e-6, "{obj}");
    assert_eq!(route(&g, "indicator"), "big_m");
}

#[test]
fn indicator_without_finite_bounds_is_a_typed_error() {
    let m = "Minimize\n obj: x + 3 b\nSubject To\n c: x + b >= 0\n i: b = 0 -> x <= 5\nBinaries\n b\nEnd\n";
    let p = tmp("s_ind_inf.lp", m);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["mip", p.to_str().unwrap()])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(3));
    let e = String::from_utf8_lossy(&out.stderr);
    assert!(
        e.contains("indicator 'i'") && e.contains("unbounded"),
        "{e}"
    );
}

#[test]
fn sos_solves() {
    let s1 = "Maximize\n obj: x + y\nSubject To\n c: x + y <= 100\nBounds\n 0 <= x <= 4\n 0 <= y <= 3\nSOS\n s: S1 :: x : 1 y : 2\nEnd\n";
    let (o, g) = solve("s_sos1", s1, &[]);
    assert!((o - 4.0).abs() < 1e-6);
    assert_eq!(route(&g, "sos"), "sos1_binary");
    let s2 = "Maximize\n obj: x + y + z\nSubject To\n c: x + y + z <= 100\nBounds\n 0 <= x <= 4\n 0 <= y <= 3\n 0 <= z <= 2\nSOS\n s: S2 :: x : 1 y : 2 z : 3\nEnd\n";
    let (o, _) = solve("s_sos2", s2, &[]);
    assert!((o - 7.0).abs() < 1e-6, "{o}");
}

#[test]
fn pwl_convex_is_an_lp_and_nonconvex_is_incremental() {
    let cvx = "Minimize\n obj: y - 2 x\nSubject To\n c: x >= 0\nBounds\n 0 <= x <= 4\n 0 <= y <= 9\nGeneral Constraints\n p: y = PWL ( x ) : (0, 0) (2, 1) (4, 4)\nEnd\n";
    let (o, g) = solve("s_pwl_cvx", cvx, &[]);
    assert!((o + 4.0).abs() < 1e-6, "{o}");
    assert_eq!(route(&g, "pwl"), "epigraph_lp");
    assert_eq!(g["added_columns"], 0);
    let non = "Minimize\n obj: y\nSubject To\n c: x >= 1.5\nBounds\n 0 <= x <= 3\n 0 <= y <= 9\nGeneral Constraints\n p: y = PWL ( x ) : (0, 0) (1, 2) (2, 1) (3, 3)\nEnd\n";
    let (o, g) = solve("s_pwl_non", non, &[]);
    assert!((o - 1.0).abs() < 1e-6, "{o}");
    assert_eq!(route(&g, "pwl"), "incremental");
}

#[test]
fn functions_solve() {
    let abs = "Minimize\n obj: x\nSubject To\n c: r >= 3\nBounds\n -5 <= x <= 5\n 0 <= r <= 9\nGeneral Constraints\n a: r = ABS ( x )\nEnd\n";
    let (o, g) = solve("s_abs", abs, &[]);
    assert!((o + 5.0).abs() < 1e-6, "{o}");
    assert_eq!(route(&g, "abs"), "abs_binary");
    let abs_lp = "Minimize\n obj: r\nSubject To\n c: x >= -2\nBounds\n -5 <= x <= 5\n 0 <= r <= 9\nGeneral Constraints\n a: r = ABS ( x )\nEnd\n";
    let (o, g) = solve("s_abs_lp", abs_lp, &[]);
    assert!(o.abs() < 1e-6, "{o}");
    assert_eq!(route(&g, "abs"), "abs_lp");
    let mx = "Minimize\n obj: - r\nSubject To\n c: x + y >= 0\nBounds\n 0 <= x <= 3\n 0 <= y <= 2\n 0 <= r <= 9\nGeneral Constraints\n m: r = MAX ( x , y , 1 )\nEnd\n";
    let (o, g) = solve("s_max", mx, &[]);
    assert!((o + 3.0).abs() < 1e-6, "{o}");
    assert_eq!(route(&g, "max"), "max_binary");
    let mn = "Minimize\n obj: r\nSubject To\n c: x + y >= 0\nBounds\n 0 <= x <= 3\n 0 <= y <= 2\n -9 <= r <= 9\nGeneral Constraints\n m: r = MIN ( x , y , 1 )\nEnd\n";
    // min is concave: minimizing r is the non-convex direction -> binary route
    let (o, g) = solve("s_min", mn, &[]);
    assert!(o.abs() < 1e-6, "{o}");
    assert_eq!(route(&g, "min"), "min_binary");
    let and = "Minimize\n obj: - a\nSubject To\n c: b1 + b2 >= 0\nBinaries\n b1 b2 a\nGeneral Constraints\n g: a = AND ( b1 , b2 )\nEnd\n";
    let (o, g) = solve("s_and", and, &[]);
    assert!((o + 1.0).abs() < 1e-6);
    assert_eq!(route(&g, "and"), "logic_rows");
}

#[test]
fn plain_models_take_the_unchanged_path() {
    let m = "Minimize\n obj: x\nSubject To\n c: x >= 1\nEnd\n";
    let p = tmp("plain_solve.lp", m);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["mip", p.to_str().unwrap(), "--format", "json"])
        .output()
        .unwrap();
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert!(v.get("general_constraints").is_none());
}
