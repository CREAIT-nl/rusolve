//! Problem-class router: `rusolve solve` and `rusolve mip` classify a
//! model the same way and hand it to the same engine. Every test runs the real binary, both
//! commands, and compares status, objective, stderr, exit code and the new envelope fields.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

const LP: &str =
    "Minimize\n obj: x + y\nSubject To\n c1: x + y >= 2\nBounds\n 0 <= x <= 5\n 0 <= y <= 5\nEnd\n";
const MIP: &str = "Maximize\n obj: 3 x + 2 y\nSubject To\n c1: x + y <= 4.5\nBounds\n 0 <= x <= 3\n 0 <= y <= 3\nGeneral\n x y\nEnd\n";
const QP: &str =
    "Minimize\n obj: - x - y + [ 2 x ^ 2 + 2 y ^ 2 ] / 2\nSubject To\n c1: x + y <= 1\nEnd\n";
const QPN: &str =
    "Minimize\n obj: - x - y + [ 2 x ^ 2 + - 2 y ^ 2 ] / 2\nSubject To\n c1: x + y <= 1\nEnd\n";
const MIQP: &str = "Minimize\n obj: - x - y + [ 2 x ^ 2 + 2 y ^ 2 ] / 2\nSubject To\n c1: x + y <= 1.5\nBounds\n 0 <= x <= 10\n 0 <= y <= 10\nGeneral\n x\nEnd\n";
const SOCP: &str = "Minimize\n obj: - x - y\nSubject To\n q1: [ x ^ 2 + y ^ 2 ] <= 1\nBounds\n x free\n y free\nEnd\n";
const QCQP: &str = "Minimize\n obj: - x - y + [ 2 x ^ 2 ] / 2\nSubject To\n q1: [ x ^ 2 + y ^ 2 ] <= 1\nBounds\n x free\n y free\nEnd\n";
const QCN: &str = "Minimize\n obj: - x - y\nSubject To\n q1: [ x ^ 2 + y ^ 2 ] >= 1\nBounds\n x free\n y free\nEnd\n";
const MICONIC: &str = "Minimize\n obj: - x - y\nSubject To\n q1: [ x ^ 2 + y ^ 2 ] <= 1\nBounds\n x free\n y free\nGeneral\n x\nEnd\n";
/// An indicator row, a MAX and a binary: a general-constraint MIP, optimum 20.
const GC_MIP: &str = "Maximize\n obj: 3 x + 2 y - 10 b\nSubject To\n c1: x + y <= 10\n ind1: b = 0 -> x <= 2\nGeneral Constraints\n g1: m = MAX ( x , y )\nBounds\n 0 <= x <= 8\n 0 <= y <= 8\n m <= 7\nBinary\n b\nEnd\n";
/// Indicator on a column that is binary by bounds only is still a MIP; here the only object is a
/// MAX over continuous columns: general constraints without integers. max m = 7 at x = 7 or y = 7.
const GC: &str = "Maximize\n obj: m\nSubject To\n c1: x + y <= 10\nGeneral Constraints\n g1: m = MAX ( x , y )\nBounds\n 0 <= x <= 8\n 0 <= y <= 8\n m <= 7\nEnd\n";
const GCQ: &str = "Maximize\n obj: 3 x + 2 y - 10 b - [ x ^ 2 ] / 2\nSubject To\n c1: x + y <= 10\n ind1: b = 0 -> x <= 2\nBounds\n 0 <= x <= 8\n 0 <= y <= 8\nBinary\n b\nEnd\n";

struct Run {
    code: i32,
    stdout: String,
    stderr: String,
}

fn run(cmd: &str, model: &str, name: &str, extra: &[&str]) -> Run {
    let dir: PathBuf = std::env::temp_dir().join(format!("rusolve-routing-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let ext = if model.starts_with("NAME") {
        "mps"
    } else {
        "lp"
    };
    let path = dir.join(format!("{name}.{ext}"));
    std::fs::write(&path, model).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg(cmd)
        .arg(&path)
        .args(extra)
        .output()
        .unwrap();
    Run {
        code: out.status.code().unwrap_or(-1),
        stdout: String::from_utf8_lossy(&out.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
    }
}

fn json(r: &Run) -> Value {
    serde_json::from_str(&r.stdout)
        .unwrap_or_else(|e| panic!("stdout is not json ({e}): {}", r.stdout))
}

/// Both commands solve `model` to `objective` with an empty stderr and the same `class`, `engine`
/// and `convexity` in the JSON envelope.
fn solves_alike(
    name: &str,
    model: &str,
    objective: f64,
    class: &str,
    engine: &str,
    convexity: Option<&str>,
) {
    let mut seen = Vec::new();
    for cmd in ["solve", "mip"] {
        let r = run(cmd, model, name, &["--format", "json"]);
        assert_eq!(r.code, 0, "{cmd} {name}: {}\n{}", r.stdout, r.stderr);
        // the LP path keeps its `rusolve ekk:` trace line on stderr (identity gate): only the
        // routed classes promise an empty stderr
        if !matches!(class, "lp" | "mip") {
            assert_eq!(
                r.stderr, "",
                "{cmd} {name}: stderr must be empty on success"
            );
        }
        let v = json(&r);
        assert_eq!(v["status"], "Optimal", "{cmd} {name}: {v}");
        let obj = v["objective"].as_f64().unwrap();
        assert!(
            (obj - objective).abs() < 1e-5,
            "{cmd} {name}: {obj} vs {objective}"
        );
        assert_eq!(v["class"], class, "{cmd} {name}");
        if class == "convex_qp" {
            // the engine that ran, named as the `qp` block names it (hipo, or quass on a fallback)
            assert_eq!(v["engine"], v["qp"]["engine"], "{cmd} {name}");
        } else {
            assert_eq!(v["engine"], engine, "{cmd} {name}");
        }
        match convexity {
            Some(c) => assert_eq!(v["convexity"], c, "{cmd} {name}"),
            None => assert!(v.get("convexity").is_none(), "{cmd} {name}: {v}"),
        }
        seen.push(obj);
    }
    assert!(
        (seen[0] - seen[1]).abs() < 1e-9,
        "{name}: solve and mip disagree: {seen:?}"
    );
}

/// Both commands refuse with exit 3, a stderr that names the class, and a typed JSON error.
fn refuses_alike(name: &str, model: &str, class: &str, convexity: Option<&str>) {
    for cmd in ["solve", "mip"] {
        let r = run(cmd, model, name, &["--format", "json"]);
        assert_eq!(r.code, 3, "{cmd} {name}: {}\n{}", r.stdout, r.stderr);
        assert!(r.stderr.contains(class), "{cmd} {name}: {}", r.stderr);
        let v = json(&r);
        assert_eq!(v["schema_version"], "rusolve.error/1", "{cmd} {name}: {v}");
        assert_eq!(v["error"]["class"], class, "{cmd} {name}: {v}");
        // a QP refusal keeps the QP front's kind, eigenvalue and variables
        assert!(
            [
                "unsupported_class",
                "nonconvex_quadratic_objective",
                "miqp_unsupported"
            ]
            .contains(&v["error"]["kind"].as_str().unwrap()),
            "{v}"
        );
        assert_eq!(
            v["error"].get("convexity").and_then(Value::as_str),
            convexity
        );
        // never points to the other command, which refuses too
        for other in ["`rusolve solve`", "`rusolve mip`", "use rusolve"] {
            assert!(!r.stderr.contains(other), "{cmd} {name}: {}", r.stderr);
        }
        let text = run(cmd, model, name, &[]);
        assert_eq!(text.code, 3, "{cmd} {name} text mode");
        assert!(text.stdout.is_empty(), "{cmd} {name}: {}", text.stdout);
    }
}

#[test]
fn lp_keeps_its_path_and_reports_the_class() {
    solves_alike("lp", LP, 2.0, "lp", "highs-port", None);
}

#[test]
fn mip_keeps_its_path_and_reports_the_class() {
    // x + y <= 4.5, integers: 3x + 2y, x = 3, y = 1 -> 11. `solve` on a MIP is the LP
    // relaxation, as before the router (12), so only `mip` is compared with the integer optimum.
    let r = run("mip", MIP, "mip", &["--format", "json"]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    let v = json(&r);
    assert_eq!(v["class"], "mip");
    assert_eq!(v["engine"], "highs-port");
    assert!(
        (v["objective"].as_f64().unwrap() - 11.0).abs() < 1e-6,
        "{v}"
    );
    let r = run("solve", MIP, "mip", &["--format", "json"]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    let v = json(&r);
    // the class is what ran: the solve path hands the port the relaxation
    assert_eq!(v["class"], "lp");
    assert!(
        (v["objective"].as_f64().unwrap() - 12.0).abs() < 1e-6,
        "{v}"
    );
}

#[test]
fn convex_qp_solves_from_both_commands() {
    solves_alike("qp", QP, -0.5, "convex_qp", "hipo", Some("convex"));
}

#[test]
fn conic_solves_from_both_commands() {
    solves_alike(
        "socp",
        SOCP,
        -(2f64.sqrt()),
        "conic",
        "clarabel",
        Some("convex"),
    );
    solves_alike("qcqp", QCQP, -1.16509, "conic", "clarabel", Some("convex"));
}

#[test]
fn general_constraints_solve_from_both_commands() {
    solves_alike("gc_mip", GC_MIP, 20.0, "gc_mip", "gc-linearization", None);
    solves_alike("gc", GC, 7.0, "gc", "gc-linearization", None);
}

#[test]
fn nonconvex_qp_is_refused_by_both_commands() {
    refuses_alike("qpn", QPN, "nonconvex_qp", Some("nonconvex"));
}

#[test]
fn nonconvex_qcqp_is_refused_by_both_commands() {
    refuses_alike("qcn", QCN, "nonconvex_qcqp", Some("nonconvex"));
}

#[test]
fn miqp_is_refused_by_both_commands() {
    refuses_alike("miqp", MIQP, "miqp", None);
    // the refusal names the experimental engine and the exact option, and why it is off
    let r = run("mip", MIQP, "miqp_hint", &["--format", "json"]);
    for needle in [
        "--opt miqp_oa=true",
        "experimental",
        "off by default",
        "slower",
    ] {
        assert!(r.stderr.contains(needle), "{needle}: {}", r.stderr);
    }
    assert!(
        json(&r)["error"]["hint"]
            .as_str()
            .unwrap()
            .contains("miqp_oa=true")
    );
    // a nonconvex Q has no lever to point to
    let nonconvex = MIQP.replace("+ 2 y ^ 2", "+ - 2 y ^ 2");
    let r = run("mip", &nonconvex, "miqp_nohint", &[]);
    assert!(!r.stderr.contains("miqp_oa"), "{}", r.stderr);
}

#[test]
fn mixed_integer_conic_is_refused_by_both_commands() {
    refuses_alike("miconic", MICONIC, "mi_conic", None);
    let r = run("solve", MICONIC, "miconic_hint", &["--format", "json"]);
    for needle in [
        "--opt misocp_oa=true",
        "experimental",
        "off by default",
        "slower",
    ] {
        assert!(r.stderr.contains(needle), "{needle}: {}", r.stderr);
    }
    assert!(
        json(&r)["error"]["hint"]
            .as_str()
            .unwrap()
            .contains("misocp_oa=true")
    );
}

#[test]
fn general_constraints_with_quadratic_terms_are_refused_by_both_commands() {
    refuses_alike("gcq", GCQ, "gc_quadratic", None);
}

/// MPS `INDICATORS` with a `QCMATRIX` on another row: the indicator row `ind1` comes first, so
/// the reader must renumber `q1` when it removes `ind1`. The model then reaches the router like
/// its LP twin and gets the same typed refusal, not an MPS parse error.
const GCQ_MPS: &str = "NAME          GCQMPS
ROWS
 N  obj
 L  ind1
 L  q1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    z  obj  -1  ind1  0
    MARKER  'MARKER'  'INTEND'
    x  obj  -1  ind1  1
    y  obj  -1  q1  0
RHS
    rhs  ind1  1  q1  4
BOUNDS
 UP bnd z 1
 UP bnd x 10
 UP bnd y 10
QCMATRIX   q1
    x  x  1
    y  y  1
INDICATORS
 IF ind1  z  1
ENDATA
";

#[test]
fn mps_indicators_with_qcmatrix_are_refused_like_the_lp_form() {
    refuses_alike("gcq_mps", GCQ_MPS, "gc_quadratic", None);
}

#[test]
fn mps_qcmatrix_on_an_indicator_row_is_a_reader_error_naming_the_row() {
    let model = GCQ_MPS.replace("QCMATRIX   q1", "QCMATRIX   ind1");
    for cmd in ["solve", "mip"] {
        let r = run(cmd, &model, "gcq_mps_indrow", &[]);
        assert_eq!(r.code, 3, "{cmd}: {}\n{}", r.stdout, r.stderr);
        assert!(
            r.stderr
                .contains("indicator row 'ind1' has a QCMATRIX section"),
            "{cmd}: {}",
            r.stderr
        );
    }
}

#[test]
fn classify_reports_the_router_class() {
    for (name, model, class) in [
        ("lp", LP, "lp"),
        ("mip", MIP, "mip"),
        ("qp", QP, "convex_qp"),
        ("socp", SOCP, "conic"),
        ("gc", GC, "gc"),
    ] {
        let dir = std::env::temp_dir().join(format!("rusolve-routing-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("{name}-classify.lp"));
        std::fs::write(&path, model).unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
            .args(["classify"])
            .arg(&path)
            .args(["--format", "json"])
            .output()
            .unwrap();
        if out.status.success() {
            let v: Value = serde_json::from_slice(&out.stdout).unwrap();
            assert_eq!(v["class"], class, "{name}");
        } else {
            // classify has no general-constraint or conic reader yet: it refuses, it does not guess
            assert!(
                matches!(class, "conic" | "gc"),
                "{name}: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        }
    }
}

#[test]
fn flags_a_class_cannot_honour_are_usage_errors() {
    let r = run("mip", QP, "qp_flags", &["--node-limit", "5"]);
    assert_eq!(r.code, 2, "{}", r.stderr);
    let r = run("solve", GC_MIP, "gc_flags", &["--duals"]);
    assert_eq!(r.code, 2, "{}", r.stderr);
}

/// Convex MIQP with the `miqp_oa` lever: solved by outer approximation, class `miqp`, engine
/// `oa`, and the answer passes `rusolve check`.
#[test]
fn convex_miqp_with_the_lever_solves_and_checks() {
    let sol = std::env::temp_dir().join(format!("rusolve-routing-{}-oa.sol", std::process::id()));
    let sol_arg = sol.to_string_lossy().into_owned();
    let r = run(
        "mip",
        MIQP,
        "miqp_oa_on",
        &[
            "--opt",
            "miqp_oa=true",
            "--format",
            "json",
            "--solution-file",
            &sol_arg,
        ],
    );
    assert_eq!(r.code, 0, "{}\n{}", r.stdout, r.stderr);
    let v = json(&r);
    assert_eq!(v["class"], "miqp");
    assert_eq!(v["engine"], "oa");
    assert_eq!(v["convexity"], "convex");
    assert_eq!(v["miqp_oa"]["engine"], "oa");
    // x integer, y continuous, x + y <= 1.5: both x = 0 and x = 1 give -0.25
    let obj = v["objective"].as_f64().unwrap();
    assert!((obj + 0.25).abs() < 1e-4, "objective {obj}");
    let model = std::env::temp_dir()
        .join(format!("rusolve-routing-{}", std::process::id()))
        .join("miqp_oa_on.lp");
    let c = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("check")
        .arg(&model)
        .arg(&sol)
        .args(["--format", "json"])
        .output()
        .unwrap();
    let cv: Value = serde_json::from_slice(&c.stdout).unwrap();
    assert_eq!(cv["verdict"], "feasible", "{cv}");
    assert!((cv["objective"]["value"].as_f64().unwrap() - obj).abs() < 1e-6);
}

/// Without the lever the convex MIQP is refused as before; a nonconvex Q with an integer column
/// stays refused even with the lever on.
#[test]
fn miqp_oa_lever_off_or_nonconvex_is_still_refused() {
    let r = run("mip", MIQP, "miqp_oa_off", &["--opt", "miqp_oa=false"]);
    assert_eq!(r.code, 3, "{}\n{}", r.stdout, r.stderr);
    let nonconvex = MIQP.replace("+ 2 y ^ 2", "+ - 2 y ^ 2");
    let r = run(
        "mip",
        &nonconvex,
        "miqp_oa_nonconvex",
        &["--opt", "miqp_oa=true"],
    );
    assert_eq!(r.code, 3, "{}\n{}", r.stdout, r.stderr);
}

/// Mixed-integer conic model with the `misocp_oa` lever: solved by outer approximation, class
/// `mi_conic`, engine `oa`, and the answer passes `rusolve check`.
#[test]
fn mixed_integer_conic_with_the_lever_solves_and_checks() {
    let sol = std::env::temp_dir().join(format!("rusolve-routing-{}-soc.json", std::process::id()));
    let sol_arg = sol.to_string_lossy().into_owned();
    let r = run(
        "mip",
        MICONIC,
        "misocp_oa_on",
        &[
            "--opt",
            "misocp_oa=true",
            "--format",
            "json",
            "--solution-file",
            &sol_arg,
        ],
    );
    assert_eq!(r.code, 0, "{}\n{}", r.stdout, r.stderr);
    let v = json(&r);
    assert_eq!(v["class"], "mi_conic");
    assert_eq!(v["engine"], "oa");
    assert_eq!(v["convexity"], "convex");
    assert_eq!(v["status"], "Optimal", "{v}");
    // x integer in the unit disc: x = 0, y = 1 or x = 1, y = 0 give -1
    let obj = v["objective"].as_f64().unwrap();
    assert!((obj + 1.0).abs() < 1e-4, "objective {obj}");
    let model = std::env::temp_dir()
        .join(format!("rusolve-routing-{}", std::process::id()))
        .join("misocp_oa_on.lp");
    let c = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("check")
        .arg(&model)
        .arg(&sol)
        .args(["--format", "json"])
        .output()
        .unwrap();
    let cv: Value = serde_json::from_slice(&c.stdout).unwrap();
    assert_eq!(cv["verdict"], "feasible", "{cv}");
    assert!((cv["objective"]["value"].as_f64().unwrap() - obj).abs() < 1e-6);
}

/// Without the lever the model is refused as before; a nonconvex quadratic row with an integer
/// column stays refused even with the lever on.
#[test]
fn misocp_oa_lever_off_or_nonconvex_is_still_refused() {
    let r = run(
        "mip",
        MICONIC,
        "misocp_oa_off",
        &["--opt", "misocp_oa=false"],
    );
    assert_eq!(r.code, 3, "{}\n{}", r.stdout, r.stderr);
    let nonconvex = MICONIC.replace("<= 1", ">= 1");
    let r = run(
        "mip",
        &nonconvex,
        "misocp_oa_nonconvex",
        &["--opt", "misocp_oa=true"],
    );
    assert_eq!(r.code, 3, "{}\n{}", r.stdout, r.stderr);
}
