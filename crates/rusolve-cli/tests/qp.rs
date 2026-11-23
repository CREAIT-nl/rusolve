//! Convex QP front: the convexity check, the `qp`
//! block of the envelope and the typed refusals. Every test runs the real binary.

use serde_json::Value;
use std::path::{Path, PathBuf};
use std::process::Command;

/// min -x - 2y + (x^2 + y^2)/2, x + y <= 1, 0 <= x, y <= 10. Optimum (0, 1), objective -1.5.
const CONVEX: &str = "\
Minimize
 obj: - x - 2 y + [ x ^ 2 + y ^ 2 ] / 2
Subject To
 c1: x + y <= 1
Bounds
 0 <= x <= 10
 0 <= y <= 10
End
";

/// Q = [[1, 2], [2, 1]] has eigenvalues 3 and -1.
const NONCONVEX: &str = "\
Minimize
 obj: - x - 2 y + [ x ^ 2 + 4 x * y + y ^ 2 ] / 2
Subject To
 c1: x + y <= 1
Bounds
 0 <= x <= 10
 0 <= y <= 10
End
";

const MIQP: &str = "\
Minimize
 obj: - x - 2 y + [ x ^ 2 + y ^ 2 ] / 2
Subject To
 c1: x + y <= 1.5
Bounds
 0 <= x <= 10
 0 <= y <= 10
General
 x
End
";

fn write(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-qp-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn run(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(args)
        .output()
        .expect("run rusolve-cli");
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[test]
fn convex_qp_solves_and_the_envelope_carries_the_class() {
    let p = write("convex.lp", CONVEX);
    let (code, out, _) = run(&["solve", p.to_str().unwrap(), "--format", "json"]);
    assert_eq!(code, 0);
    let v: Value = serde_json::from_str(&out).unwrap();
    assert_eq!(v["status"], "Optimal");
    assert!((v["objective"].as_f64().unwrap() + 1.5).abs() < 1e-6);
    let qp = &v["qp"];
    assert_eq!(qp["class"], "DCL");
    assert_eq!(qp["convexity"], "convex");
    // Engine choice: HiPO first; a non-Optimal HiPO would fall back to QUASS and say so.
    assert_eq!(qp["engine_selection"], "default");
    assert_eq!(
        qp["engine"],
        if qp["fallback"] == true {
            "quass"
        } else {
            "hipo"
        }
    );
    assert_eq!(qp["shift"], 0.0);
}

#[test]
fn a_successful_qp_solve_leaves_stderr_empty_in_both_formats() {
    // Agents read stderr as a failure signal: class, engine and status go to stdout.
    let p = write("quiet.lp", CONVEX);
    let (code, out, err) = run(&["solve", p.to_str().unwrap(), "--format", "json"]);
    assert_eq!(code, 0);
    assert_eq!(err, "", "json mode");
    serde_json::from_str::<Value>(&out).expect("stdout is one JSON document");
    let (code, out, err) = run(&["solve", p.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert_eq!(err, "", "text mode");
    assert!(out.contains("qp: class DCL"), "{out}");
    assert!(out.contains("engine "), "{out}");
}

#[test]
fn an_explicit_solver_is_honoured_and_never_falls_back() {
    let p = write("explicit.lp", CONVEX);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["solve", p.to_str().unwrap(), "--format", "json"])
        .env("RUSOLVE_EKK_OPTS", "solver=simplex")
        .output()
        .expect("run rusolve-cli");
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["qp"]["engine"], "quass");
    assert_eq!(v["qp"]["engine_selection"], "user");
    assert_eq!(v["qp"]["fallback"], false);
    assert!(out.stderr.is_empty());
}

#[test]
fn an_lp_envelope_has_no_qp_block() {
    let p = write(
        "lp.lp",
        "Minimize\n obj: x + y\nSubject To\n c1: x + y >= 1\nEnd\n",
    );
    let (code, out, _) = run(&["solve", p.to_str().unwrap(), "--format", "json"]);
    assert_eq!(code, 0);
    let v: Value = serde_json::from_str(&out).unwrap();
    assert!(v.get("qp").is_none());
}

#[test]
fn nonconvex_q_is_a_typed_input_error_naming_variables_and_eigenvalue() {
    let p = write("nonconvex.lp", NONCONVEX);
    let (code, out, err) = run(&["solve", p.to_str().unwrap(), "--format", "json"]);
    assert_eq!(code, 3);
    assert!(err.contains("nonconvex"), "{err}");
    let v: Value = serde_json::from_str(&out).unwrap();
    assert_eq!(v["schema_version"], "rusolve.error/1");
    let e = &v["error"];
    assert_eq!(e["kind"], "nonconvex_quadratic_objective");
    assert!((e["eigenvalue"].as_f64().unwrap() + 1.0).abs() < 1e-6);
    let vars: Vec<&str> = e["variables"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_str().unwrap())
        .collect();
    assert!(vars.contains(&"x") && vars.contains(&"y"), "{vars:?}");
}

#[test]
fn miqp_is_refused_not_relaxed() {
    let p = write("miqp.lp", MIQP);
    let (code, out, _) = run(&["solve", p.to_str().unwrap(), "--format", "json"]);
    assert_eq!(code, 3);
    let v: Value = serde_json::from_str(&out).unwrap();
    assert_eq!(v["error"]["kind"], "miqp_unsupported");
}

#[test]
fn mip_solves_a_convex_qp_like_solve() {
    // Router: `mip` used to refuse a Hessian ("use rusolve solve"); both commands now take
    // the QP route, so the answers agree.
    let p = write("convex_mip.lp", CONVEX);
    let mut objs = Vec::new();
    for cmd in ["solve", "mip"] {
        let (code, out, err) = run(&[cmd, p.to_str().unwrap(), "--format", "json"]);
        assert_eq!(code, 0, "{cmd}: {err}");
        let v: Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["class"], "convex_qp");
        objs.push(v["objective"].as_f64().unwrap());
    }
    assert!(
        (objs[0] - objs[1]).abs() < 1e-9 && (objs[0] + 1.5).abs() < 1e-6,
        "{objs:?}"
    );
}

// ---- the checker on QP answers -----------------------------------------------------------

/// Solves `text`, writes the raw solution file with duals and returns `(model path, solution
/// path)`.
fn solve_to_file(name: &str, text: &str) -> (PathBuf, PathBuf) {
    let model = write(name, text);
    let sol = model.with_extension("sol");
    // The planted-error tests need an exact vertex answer, so they name the active-set engine;
    // the default engine's answer is checked in `the_default_engines_answer_passes_the_checker`.
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["solve", model.to_str().unwrap(), "--solution-file"])
        .arg(&sol)
        .env("RUSOLVE_EKK_OPTS", "solver=simplex")
        .output()
        .expect("run rusolve-cli");
    assert_eq!(out.status.code(), Some(0));
    (model, sol)
}

#[test]
fn the_default_engines_answer_passes_the_checker() {
    let model = write("default_engine.lp", CONVEX);
    let sol = model.with_extension("sol");
    let (code, _, err) = run(&[
        "solve",
        model.to_str().unwrap(),
        "--solution-file",
        sol.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "{err}");
    assert_eq!(check(&model, &sol)["qp"]["kkt_verdict"], "optimal");
}

fn check(model: &Path, sol: &Path) -> Value {
    let (code, out, err) = run(&[
        "check",
        model.to_str().unwrap(),
        sol.to_str().unwrap(),
        "--format",
        "json",
    ]);
    assert_eq!(code, 0, "{err}");
    serde_json::from_str(&out).unwrap()
}

fn assert_kkt_optimal(rep: &Value, objective: f64) {
    assert_eq!(rep["verdict"], "feasible", "{rep}");
    let q = &rep["qp"];
    assert_eq!(q["kkt_verdict"], "optimal", "{q}");
    assert!((rep["objective"]["value"].as_f64().unwrap() - objective).abs() < 1e-6);
    assert!(q["duality_gap_rel"].as_f64().unwrap() < 1e-6, "{q}");
}

#[test]
fn checker_confirms_a_convex_qp_answer_and_recomputes_the_objective_with_q() {
    let (m, s) = solve_to_file("chk_convex.lp", CONVEX);
    let rep = check(&m, &s);
    assert_kkt_optimal(&rep, -1.5);
    // the linear part alone is -2: the quadratic term is what the checker adds
    assert!((rep["qp"]["objective_quadratic"].as_f64().unwrap() - 0.5).abs() < 1e-9);
}

#[test]
fn checker_confirms_a_maximisation_with_equality_and_two_sided_bounds() {
    // max 3x + 2y + z - (x^2 + y^2 + 2 z^2 + x y)/2, x + y + z = 4, 0 <= x <= 1, y free-ish, z >= 0
    let text = "\
Maximize
 obj: 3 x + 2 y + z - [ x ^ 2 + y ^ 2 + 2 z ^ 2 + 2 x * y ] / 2
Subject To
 e: x + y + z = 4
 r: y - z <= 2
Bounds
 0 <= x <= 1
 -5 <= y <= 5
 z >= 0
End
";
    let (m, s) = solve_to_file("chk_max.lp", text);
    let rep = check(&m, &s);
    assert_eq!(rep["verdict"], "feasible", "{rep}");
    assert_eq!(rep["objective"]["sense"], "maximize");
    assert_eq!(rep["qp"]["kkt_verdict"], "optimal", "{}", rep["qp"]);
}

#[test]
fn checker_flags_a_planted_wrong_primal_point() {
    let (m, s) = solve_to_file("chk_bad_primal.lp", CONVEX);
    let text = std::fs::read_to_string(&s).unwrap();
    // move y off the optimum (still feasible, wrong objective): stationarity must fail
    let bad = text.replacen("y 1\n", "y 0.5\n", 1);
    assert_ne!(bad, text);
    let s2 = s.with_file_name("chk_bad_primal_moved.sol");
    std::fs::write(&s2, bad).unwrap();
    let rep = check(&m, &s2);
    assert_eq!(rep["verdict"], "feasible");
    assert_eq!(rep["qp"]["kkt_verdict"], "not_optimal", "{}", rep["qp"]);
    assert!(rep["qp"]["stationarity"].as_f64().unwrap() > 1e-3);
}

#[test]
fn checker_flags_a_planted_wrong_dual_sign_and_a_wrong_gap() {
    let (m, s) = solve_to_file("chk_bad_dual.lp", CONVEX);
    let text = std::fs::read_to_string(&s).unwrap();
    let dual_at = text.find("# Dual solution values").unwrap();
    let (head, tail) = text.split_at(dual_at);
    // row dual -1 becomes +1: wrong sign for an active <= row of a minimisation
    let tail = tail.replacen("c1 -0.9999999", "c1 1", 1);
    let s2 = s.with_file_name("chk_bad_dual_flipped.sol");
    std::fs::write(&s2, format!("{head}{tail}")).unwrap();
    let rep = check(&m, &s2);
    assert_eq!(rep["qp"]["kkt_verdict"], "not_optimal", "{}", rep["qp"]);
}

#[test]
fn checker_flags_an_infeasible_qp_point_and_works_without_duals() {
    let (m, s) = solve_to_file("chk_infeasible.lp", CONVEX);
    let text = std::fs::read_to_string(&s).unwrap();
    let s2 = s.with_file_name("chk_infeasible_x.sol");
    std::fs::write(&s2, text.replacen("x 0\n", "x 3\n", 1)).unwrap();
    assert_eq!(check(&m, &s2)["verdict"], "infeasible");
    // no dual section: the objective is still recomputed, the KKT verdict is unknown
    let cut = &text[..text.find("# Dual solution values").unwrap()];
    let s3 = s.with_file_name("chk_nodual.sol");
    std::fs::write(&s3, format!("{cut}\n# Dual solution values\nNone\n")).unwrap();
    let rep = check(&m, &s3);
    assert_eq!(rep["qp"]["kkt_verdict"], "unknown");
    assert_eq!(rep["qp"]["duals_present"], false);
    assert!((rep["objective"]["value"].as_f64().unwrap() + 1.5).abs() < 1e-6);
}

#[test]
fn an_lp_check_report_has_no_qp_block() {
    let (m, s) = solve_to_file(
        "chk_lp.lp",
        "Minimize\n obj: x + y\nSubject To\n c1: x + y >= 1\nEnd\n",
    );
    assert!(check(&m, &s).get("qp").is_none());
}

#[test]
fn qp_block_has_exactly_the_keys_of_the_check_schema() {
    let schema: Value = serde_json::from_str(include_str!(
        "../../../docs/schema/rusolve-check-1.schema.json"
    ))
    .unwrap();
    let props = schema["properties"]["qp"]["properties"]
        .as_object()
        .unwrap();
    let (m, s) = solve_to_file("chk_schema.lp", CONVEX);
    let rep = check(&m, &s);
    let got = rep["qp"].as_object().unwrap();
    let mut a: Vec<_> = props.keys().collect();
    let mut b: Vec<_> = got.keys().collect();
    a.sort();
    b.sort();
    assert_eq!(a, b);
}
