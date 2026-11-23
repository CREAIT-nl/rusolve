//! `rusolve ranging`, sensitivity by name, run as the real binary.
//!
//! The Netlib comparison fixture `fixtures/ranging_afiro_highspy.json` was
//! written by highspy 1.15.1 (`Highs.getRanging` on afiro, `presolve=off`); the
//! port is HiGHS 1.16.0, whose simplex and ranging code is unchanged for this
//! case. Regenerate with a script that dumps `getRanging()[1]`, `getSolution()`
//! and the names (in_var/out_var are simplex indices: columns, then rows).

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-ranging-1.schema.json");
const FIXTURE: &str = include_str!("fixtures/ranging_afiro_highspy.json");

/// min -3x - 4y + 5z + w, x + 2y + z <= 8, z + 3w <= 10, x<=4, y<=5, z<=3 integer.
/// Optimum x=4 y=2 z=0 w=0, objective -20.
const TINY_MIP: &str = "\
NAME RNG1
ROWS
 N  cost
 L  cap1
 L  cap2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  cost  -3  cap1  1
    y  cost  -4  cap1  2
    z  cost   5  cap1  1  cap2  1
    MARKER  'MARKER'  'INTEND'
    w  cost   1  cap2  3
RHS
    rhs  cap1  8  cap2  10
BOUNDS
 UP bnd x 4
 UP bnd y 5
 UP bnd z 3
ENDATA
";

fn afiro() -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../benches/data/netlib/afiro.mps")
        .to_str()
        .unwrap()
        .to_string()
}

fn write_tmp(name: &str, text: &str) -> String {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-ranging-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p.to_str().unwrap().to_string()
}

fn ranging(model: &str, extra: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("ranging")
        .arg(model)
        .args(extra)
        .output()
        .unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

fn json(model: &str, extra: &[&str]) -> Value {
    let (code, stdout, stderr) = ranging(model, &[&["--format", "json"][..], extra].concat());
    assert_eq!(code, 0, "stderr: {stderr}");
    serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not one JSON object: {e}\n{stdout}"))
}

fn close(got: &Value, want: &Value, what: &str) {
    match (got.as_f64(), want.as_f64()) {
        (None, None) => assert!(got.is_null() && want.is_null(), "{what}: {got} vs {want}"),
        (Some(g), Some(w)) => {
            let tol = 1e-7 * w.abs().max(1.0);
            assert!((g - w).abs() <= tol, "{what}: got {g}, highspy {w}");
        }
        _ => panic!("{what}: got {got}, highspy {want}"),
    }
}

/// The name a simplex variable index stands for in the fixture.
fn fixture_var(fx: &Value, v: i64) -> Value {
    let cols = fx["col_names"].as_array().unwrap();
    let rows = fx["row_names"].as_array().unwrap();
    let v = v as usize;
    if v < cols.len() {
        serde_json::json!({"kind": "col", "name": cols[v]})
    } else if v - cols.len() < rows.len() {
        serde_json::json!({"kind": "row", "name": rows[v - cols.len()]})
    } else {
        Value::Null
    }
}

/// in_var/out_var are exact except for ties. afiro's ratio tests tie to 1e-15
/// (the limit values agree to 1e-7 and are checked above), and which of the tied
/// variables is reported depends on last-bit rounding, which differs between
/// highspy 1.15.1 and the 1.16.0 port: the leaving variable is one of X06, X14,
/// X16 and the entering one X15 (a column) or X27 (a row) in exactly those cases.
fn same_var(role: &str, got: &Value, want: &Value) -> bool {
    if got == want {
        return true;
    }
    let tied: &[&str] = if role == "out_var" {
        &["X06", "X14", "X16"]
    } else {
        &["X15", "X27"]
    };
    let is_tied = |v: &Value| tied.contains(&v["name"].as_str().unwrap_or(""));
    is_tied(got) && is_tied(want)
}

fn check_side(side: &Value, rec: &Value, i: usize, fx: &Value, what: &str) {
    close(&side["value"], &rec["value"][i], &format!("{what} value"));
    close(
        &side["objective"],
        &rec["objective"][i],
        &format!("{what} objective"),
    );
    for (k, f) in [("in_var", "in_var"), ("out_var", "out_var")] {
        let idx = rec[f][i].as_i64().unwrap();
        let want = if idx < 0 {
            Value::Null
        } else {
            fixture_var(fx, idx)
        };
        assert!(
            same_var(k, &side[k], &want),
            "{what} {k}: {} vs highspy {want}",
            side[k]
        );
    }
}

#[test]
fn ranging_afiro_matches_highspy() {
    let fx: Value = serde_json::from_str(FIXTURE).unwrap();
    let out = json(&afiro(), &[]);
    assert_eq!(out["basis_of"], "lp");
    assert_eq!(out["status"], "Optimal");
    close(&out["objective"], &fx["objective"], "objective");
    let cols = out["columns"].as_array().unwrap();
    let rows = out["rows"].as_array().unwrap();
    assert_eq!(cols.len(), fx["col_names"].as_array().unwrap().len());
    assert_eq!(rows.len(), fx["row_names"].as_array().unwrap().len());
    for (j, c) in cols.iter().enumerate() {
        let n = c["name"].as_str().unwrap();
        assert_eq!(c["name"], fx["col_names"][j]);
        close(&c["value"], &fx["col_value"][j], &format!("{n} value"));
        close(&c["reduced_cost"], &fx["col_dual"][j], &format!("{n} dual"));
        for (dir, key) in [("down", "dn"), ("up", "up")] {
            check_side(
                &c["cost_ranging"][dir],
                &fx[format!("col_cost_{key}")],
                j,
                &fx,
                &format!("{n} cost {dir}"),
            );
            check_side(
                &c["bound_ranging"][dir],
                &fx[format!("col_bound_{key}")],
                j,
                &fx,
                &format!("{n} bound {dir}"),
            );
        }
    }
    for (i, r) in rows.iter().enumerate() {
        let n = r["name"].as_str().unwrap();
        assert_eq!(r["name"], fx["row_names"][i]);
        close(
            &r["activity"],
            &fx["row_value"][i],
            &format!("{n} activity"),
        );
        close(&r["dual"], &fx["row_dual"][i], &format!("{n} dual"));
        for (dir, key) in [("down", "dn"), ("up", "up")] {
            check_side(
                &r["bound_ranging"][dir],
                &fx[format!("row_bound_{key}")],
                i,
                &fx,
                &format!("{n} row bound {dir}"),
            );
        }
    }
}

#[test]
fn ranging_names_present_and_filters_apply() {
    let all = json(&afiro(), &[]);
    for c in all["columns"].as_array().unwrap() {
        assert!(!c["name"].as_str().unwrap().is_empty());
    }
    for r in all["rows"].as_array().unwrap() {
        assert!(!r["name"].as_str().unwrap().is_empty());
    }
    let some = json(&afiro(), &["--cols", "X01,X02", "--rows", "R09"]);
    let cols: Vec<&str> = some["columns"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c["name"].as_str().unwrap())
        .collect();
    assert_eq!(cols, ["X01", "X02"]);
    assert_eq!(some["rows"].as_array().unwrap().len(), 1);
    assert_eq!(some["rows"][0]["name"], "R09");
    // a glob and an unknown name
    let glob = json(&afiro(), &["--cols", "X0*", "--rows", "R0*"]);
    assert!(glob["columns"].as_array().unwrap().len() >= 2);
    let (code, _, err) = ranging(&afiro(), &["--cols", "NOPE"]);
    assert_eq!(code, 2, "{err}");
    // text mode is not JSON and names the rows
    let (code, text, _) = ranging(&afiro(), &["--rows", "R09"]);
    assert_eq!(code, 0);
    assert!(text.contains("R09") && !text.trim_start().starts_with('{'));
    let (code, _, _) = ranging("/nonexistent.mps", &[]);
    assert_eq!(code, 3);
}

#[test]
fn ranging_mip_is_labelled_fixed_integer_lp() {
    let model = write_tmp("tiny.mps", TINY_MIP);
    let out = json(&model, &[]);
    assert_eq!(out["basis_of"], "fixed_integer_lp");
    assert_eq!(out["mip_status"], "Optimal");
    assert_eq!(out["status"], "Optimal");
    assert_eq!(out["input"]["integer_cols"], 3);
    close(&out["objective"], &serde_json::json!(-20.0), "objective");
    let by_name = |n: &str| {
        out["columns"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["name"] == n)
            .unwrap()
            .clone()
    };
    // the integer columns are fixed at the incumbent, so their bounds in the
    // fixed LP are the incumbent values
    close(&by_name("x")["value"], &serde_json::json!(4.0), "x");
    close(&by_name("x")["lower"], &serde_json::json!(4.0), "x lower");
    close(&by_name("y")["upper"], &serde_json::json!(2.0), "y upper");
    // an LP is labelled 'lp', not fixed_integer_lp
    assert_eq!(json(&afiro(), &["--rows", "R09"])["basis_of"], "lp");
}

#[test]
fn ranging_json_has_every_required_field_of_the_schema() {
    let schema: Value = serde_json::from_str(SCHEMA).unwrap();
    let req = |node: &Value| -> Vec<String> {
        node["required"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect()
    };
    let check = |node: &Value, def: &Value, what: &str| {
        for f in req(def) {
            assert!(node.get(&f).is_some(), "{what}.{f} missing from {node}");
        }
    };
    for model in [afiro(), write_tmp("tiny2.mps", TINY_MIP)] {
        let out = json(&model, &[]);
        check(&out, &schema, "top");
        check(&out["solver"], &schema["properties"]["solver"], "solver");
        check(&out["input"], &schema["properties"]["input"], "input");
        assert_eq!(out["schema_version"], "rusolve.ranging/1");
        let defs = &schema["$defs"];
        for c in out["columns"].as_array().unwrap() {
            check(c, &defs["column"], "column");
            for r in ["cost_ranging", "bound_ranging"] {
                for d in ["down", "up"] {
                    check(&c[r][d], &defs["side"], r);
                }
            }
        }
        for r in out["rows"].as_array().unwrap() {
            check(r, &defs["row"], "row");
            for d in ["down", "up"] {
                check(&r["bound_ranging"][d], &defs["side"], "row bound");
            }
        }
    }
}

#[test]
fn ranging_infeasible_model_reports_status_without_ranging() {
    let model = write_tmp(
        "infeasible.mps",
        "NAME INF\nROWS\n N  obj\n G  r1\nCOLUMNS\n    x  obj  1  r1  1\nRHS\n    rhs  r1  5\nBOUNDS\n UP bnd x 1\nENDATA\n",
    );
    let out = json(&model, &[]);
    assert_eq!(out["ranging_valid"], false);
    assert_eq!(out["status"], "Infeasible");
    assert!(out["columns"].as_array().unwrap().is_empty());
}
