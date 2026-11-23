//! Continuous SOCP / convex QCQP through Clarabel, end to end on the real binary:
//! `solve` and `mip` route a model with cones or quadratic rows, `check` verifies the answer,
//! a nonconvex row or an integer column is a typed error (exit 3), and an LP is untouched.

use serde_json::Value;
use std::path::PathBuf;
use std::process::{Command, Output};

const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-result-1.schema.json");
const CHECK_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-check-conic-1.schema.json");

/// min t s.t. x + y = 2, t >= ‖(x, y)‖ and the disc row x² + y² <= 4 (inactive at the optimum
/// x = y = 1, t = sqrt 2: x² + y² = 2).
const MODEL: &str = "\
NAME conic1
ROWS
 N  obj
 E  e1
 L  disc
COLUMNS
    t  obj  1
    x  e1  1
    y  e1  1
RHS
    rhs  e1  2  disc  4
BOUNDS
 FR bnd t
 FR bnd x
 FR bnd y
QCMATRIX disc
    x  x  1
    y  y  1
CSECTION k1 0.0 QUAD
    t
    x
    y
ENDATA
";

fn tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-conic-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(args)
        .output()
        .unwrap()
}

fn json_of(o: &Output) -> Value {
    serde_json::from_slice(&o.stdout).unwrap_or_else(|e| {
        panic!(
            "stdout is not JSON ({e}): {}",
            String::from_utf8_lossy(&o.stdout)
        )
    })
}

/// Every `required` key of the committed schema is present (and no key is unknown), as the
/// other envelope tests do.
fn against(schema: &Value, node: &Value, at: &str) {
    if let Some(req) = schema["required"].as_array() {
        for k in req {
            assert!(node.get(k.as_str().unwrap()).is_some(), "{at} lacks {k}");
        }
    }
    if let (Some(props), Some(obj)) = (schema["properties"].as_object(), node.as_object()) {
        if schema["additionalProperties"] == Value::Bool(false) {
            for k in obj.keys() {
                assert!(props.contains_key(k), "{at}.{k} is not in the schema");
            }
        }
        for (k, sub) in props {
            if let Some(v) = obj.get(k) {
                against(sub, v, &format!("{at}.{k}"));
            }
        }
    }
}

#[test]
fn solve_and_mip_route_a_conic_model_to_clarabel() {
    let m = tmp("m.mps", MODEL);
    for cmd in ["solve", "mip"] {
        let o = run(&[
            cmd,
            m.to_str().unwrap(),
            "--format",
            "json",
            "--include-solution",
        ]);
        assert!(
            o.status.success(),
            "{cmd}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        let v = json_of(&o);
        against(&serde_json::from_str(SCHEMA).unwrap(), &v, "$");
        assert_eq!(v["status"], "Optimal");
        assert_eq!(v["solver"]["engine"], "clarabel");
        assert!((v["objective"].as_f64().unwrap() - 2f64.sqrt()).abs() < 1e-6);
        assert_eq!(v["conic"]["checker"]["verdict"], "feasible");
        assert_eq!(v["conic"]["cone_duals"][0]["name"], "k1");
        assert_eq!(v["solution"]["columns"].as_array().unwrap().len(), 3);
    }
}

#[test]
fn the_solution_checks_through_rusolve_check_and_a_planted_violation_is_caught() {
    let m = tmp("m2.mps", MODEL);
    let v = json_of(&run(&[
        "solve",
        m.to_str().unwrap(),
        "--format",
        "json",
        "--include-solution",
    ]));
    let sol = tmp("m2.json", &serde_json::to_string(&v["solution"]).unwrap());
    let c = run(&[
        "check",
        m.to_str().unwrap(),
        sol.to_str().unwrap(),
        "--format",
        "json",
    ]);
    assert!(c.status.success(), "{}", String::from_utf8_lossy(&c.stderr));
    assert!(c.stderr.is_empty(), "check prints no stray read error");
    let r = json_of(&c);
    against(&serde_json::from_str(CHECK_SCHEMA).unwrap(), &r, "$");
    assert_eq!(r["verdict"], "feasible", "{r}");
    assert!((r["objective"]["value"].as_f64().unwrap() - 2f64.sqrt()).abs() < 1e-6);

    // planted: t below the cone, x + y off the equality
    let bad = tmp(
        "m2bad.json",
        r#"{"columns":[{"name":"t","value":0.5},{"name":"x","value":1.0},{"name":"y","value":1.0}]}"#,
    );
    let r = json_of(&run(&[
        "check",
        m.to_str().unwrap(),
        bad.to_str().unwrap(),
        "--format",
        "json",
    ]));
    assert_eq!(r["verdict"], "infeasible");
    assert_eq!(r["findings"][0]["kind"], "cone");
    assert_eq!(r["findings"][0]["name"], "k1");
}

#[test]
fn a_nonconvex_row_and_an_integer_column_are_typed_errors() {
    let nc = tmp("nc.mps", &MODEL.replace("    y  y  1", "    y  y  -1"));
    let o = run(&["solve", nc.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(3));
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(
        e.contains("quadratic row disc") && e.contains("column y"),
        "{e}"
    );

    let int = tmp("int.mps", &MODEL.replace(" FR bnd x\n", " BV bnd x\n"));
    let o = run(&["mip", int.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(3));
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(
        e.contains("mi_conic") && e.contains("column x is integer"),
        "{e}"
    );
}

#[test]
fn a_model_without_cones_or_quadratic_rows_takes_the_old_path() {
    let lp = tmp(
        "lp.mps",
        "NAME lp\nROWS\n N obj\n L c1\nCOLUMNS\n    x obj -1 c1 1\nRHS\n    rhs c1 4\nENDATA\n",
    );
    let o = run(&["solve", lp.to_str().unwrap(), "--format", "json"]);
    assert!(o.status.success());
    let v = json_of(&o);
    assert_ne!(v["solver"]["engine"], "clarabel");
    assert!(v.get("conic").is_none());
    assert_eq!(v["status"], "Optimal");
}

#[test]
fn a_conic_solve_prints_nothing_on_stderr_and_a_refusal_exits_nonzero() {
    let lp = "Minimize\n obj: x + y\nSubject To\n q1: x + [ x ^ 2 + y ^ 2 ] <= 4\nBounds\n x free\n y free\nEnd\n";
    let ok = tmp("qc.lp", lp);
    let bad = tmp("qcn.lp", &lp.replace("+ y ^ 2", "- y ^ 2"));
    for cmd in ["solve", "mip"] {
        let o = run(&[cmd, ok.to_str().unwrap(), "--format", "json"]);
        assert!(o.status.success(), "{cmd}");
        assert_eq!(String::from_utf8_lossy(&o.stderr), "", "{cmd} stderr");
        assert_eq!(json_of(&o)["status"], "Optimal");
        // the MPS spelling too
        let o = run(&[
            cmd,
            tmp("m3.mps", MODEL).to_str().unwrap(),
            "--format",
            "json",
        ]);
        assert_eq!(String::from_utf8_lossy(&o.stderr), "", "{cmd} mps stderr");

        let o = run(&[cmd, bad.to_str().unwrap(), "--format", "json"]);
        assert_eq!(o.status.code(), Some(3), "{cmd}");
        // The refusal is the typed rusolve.error/1 object, naming the class
        let err = json_of(&o);
        assert_eq!(err["schema_version"], "rusolve.error/1", "{cmd}");
        assert_eq!(err["error"]["class"], "nonconvex_qcqp", "{cmd}");
        let e = String::from_utf8_lossy(&o.stderr);
        assert!(e.contains("quadratic row q1 is not convex"), "{e}");
        assert!(!e.contains("LP parse error"), "{e}");
    }
}

#[test]
fn an_lp_with_only_an_objective_quadratic_keeps_the_legacy_path() {
    // `[` makes the pre-scan look closer, but there is no quadratic row or cone
    let qp = tmp(
        "qpobj.lp",
        "Minimize\n obj: x + [ x ^ 2 ] / 2\nSubject To\n c1: x >= 1\nEnd\n",
    );
    let o = run(&["solve", qp.to_str().unwrap(), "--format", "json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert_ne!(json_of(&o)["solver"]["engine"], "clarabel");
}
