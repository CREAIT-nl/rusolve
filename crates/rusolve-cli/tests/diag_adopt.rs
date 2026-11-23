//! Diagnostics on a run: every reader error at once, lint warnings
//! on solve and mip, the infeasible set in the envelope, `--no-warnings`.

use serde_json::Value;
use std::process::Command;

fn write(name: &str, text: &str) -> String {
    let dir = std::env::temp_dir().join("rusolve-diag-adopt-test");
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p.to_str().unwrap().to_string()
}

fn run(args: &[&str]) -> (i32, String, String) {
    let o = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(args)
        .output()
        .unwrap();
    (
        o.status.code().unwrap(),
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

const BAD: &str = "NAME BAD\nROWS\n N OBJ\n L R1\n G R2\nCOLUMNS\n X OBJ 1 R1 1\n X R7 2\n Y OBJ 1 R2 1\nRHS\n RHS R1 abc\n RHS R2 4\nENDATA\n";

#[test]
fn malformed_mps_reports_both_errors_with_position_source_and_help() {
    let p = write("bad.mps", BAD);
    let (code, out, err) = run(&["solve", &p, "--format", "json"]);
    assert_eq!(code, 3);
    let v: Value = serde_json::from_str(&out).unwrap();
    let ds = v["diagnostics"].as_array().unwrap();
    assert_eq!(ds.len(), 2, "{out}");
    assert_eq!(ds[0]["code"], "RS1004");
    assert_eq!(ds[0]["location"]["line"], 8);
    assert_eq!(ds[0]["location"]["column"], 4);
    assert_eq!(ds[0]["location"]["source_line"], " X R7 2");
    let help = ds[0]["help"].to_string();
    assert!(help.contains("declared rows: OBJ, R1, R2"), "{help}");
    assert!(help.contains("'R1'"), "did you mean: {help}");
    assert_eq!(ds[1]["code"], "RS1005");
    assert_eq!(ds[1]["location"]["line"], 11);
    assert!(
        err.contains("error[RS1004]") && err.contains("error[RS1005]"),
        "{err}"
    );
}

const INFEASIBLE: &str = "NAME INF\nROWS\n N OBJ\n L C1\n G C2\n L C3\n L C4\nCOLUMNS\n X OBJ 1 C1 1\n X C2 1 C4 1\n Y OBJ 1 C3 1\n Y C4 1\nRHS\n RHS C1 2\n RHS C2 5\n RHS C3 10\n RHS C4 20\nENDATA\n";

#[test]
fn infeasible_envelope_names_the_conflicting_rows() {
    let p = write("infeasible.mps", INFEASIBLE);
    let (code, out, err) = run(&["solve", &p, "--format", "json"]);
    assert_eq!(code, 0);
    let v: Value = serde_json::from_str(&out).unwrap();
    assert_eq!(v["status"], "Infeasible");
    let d = &v["diagnostics"][0];
    assert_eq!(d["code"], "RS2010");
    let msg = d["message"].as_str().unwrap();
    assert!(msg.contains("C1") && msg.contains("C2"), "{msg}");
    assert!(!msg.contains("C3") && !msg.contains("C4"), "{msg}");
    let help = d["help"].to_string();
    assert!(
        help.contains("rusolve explain") && help.contains("rusolve feasrelax"),
        "{help}"
    );
    assert!(err.contains("warning[RS2010]"), "{err}");
}

const BIGM: &str = "NAME BIGM\nROWS\n N OBJ\n L LINK\nCOLUMNS\n MARKER 'MARKER' 'INTORG'\n Z OBJ 1 LINK -1000000\n MARKER 'MARKER' 'INTEND'\n X OBJ -1 LINK 1\nRHS\n RHS LINK 0\nBOUNDS\n UP BND Z 1\n UP BND X 10\nENDATA\n";

#[test]
fn big_m_warns_on_solve_and_mip_and_no_warnings_silences_it() {
    let p = write("bigm.mps", BIGM);
    for cmd in ["solve", "mip"] {
        let (code, out, err) = run(&[cmd, &p, "--format", "json"]);
        assert_eq!(code, 0, "{cmd}: {err}");
        assert!(err.contains("warning[RS2001]"), "{cmd}: {err}");
        let v: Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["diagnostics"][0]["code"], "RS2001", "{cmd}: {out}");
        let (code, out, err) = run(&[cmd, &p, "--format", "json", "--no-warnings"]);
        assert_eq!(code, 0);
        assert!(!err.contains("RS2001"), "{cmd}: {err}");
        let v: Value = serde_json::from_str(&out).unwrap();
        assert!(v.get("diagnostics").is_none(), "{cmd}: {out}");
        // the answer does not depend on the warnings
        let (_, plain, _) = run(&[cmd, &p, "--format", "json"]);
        let plain: Value = serde_json::from_str(&plain).unwrap();
        assert_eq!(v["objective"], plain["objective"]);
    }
}

#[test]
fn big_m_via_singleton_row_names_row_column_and_value() {
    let p = write(
        "bigm_sing.mps",
        "NAME SING\nROWS\n N  cost\n L  c1\n L  c2\nCOLUMNS\n    x  cost  -1  c1  1\n    x  c2  1\n    MARKER  'MARKER'  'INTORG'\n    b  c1  -1000000\n    MARKER  'MARKER'  'INTEND'\nRHS\n    rhs  c1  0\n    rhs  c2  5\nBOUNDS\n UP bnd b 1\nENDATA\n",
    );
    let (code, out, err) = run(&["mip", &p, "--format", "json"]);
    assert_eq!(code, 0, "{err}");
    assert!(
        err.contains("row 'c1'")
            && err.contains("binary 'b'")
            && err.contains("tightest valid value -5"),
        "{err}"
    );
    // the location line shows the row as written, not the finding's raw evidence
    assert!(err.contains("= row: c1: x - 1e6 b <= 0"), "{err}");
    assert!(!err.contains("= row: {"), "{err}");
    let v: Value = serde_json::from_str(&out).unwrap();
    let msg = v["diagnostics"][0]["message"].as_str().unwrap();
    assert!(
        msg.contains("row 'c1'") && msg.contains("binary 'b'"),
        "{msg}"
    );
    let (_, _, err) = run(&["mip", &p, "--format", "json", "--no-warnings"]);
    assert!(!err.contains("RS2001"), "{err}");
}

#[test]
fn clean_model_has_no_diagnostics_key() {
    let p = write(
        "clean.mps",
        "NAME T\nROWS\n N obj\n L c1\nCOLUMNS\n x obj -1 c1 1\nRHS\n RHS c1 4\nBOUNDS\n UP BND x 10\nENDATA\n",
    );
    let (code, out, _) = run(&["mip", &p, "--format", "json"]);
    assert_eq!(code, 0);
    let v: Value = serde_json::from_str(&out).unwrap();
    assert!(v.get("diagnostics").is_none(), "{out}");
}

const QP_LP: &str =
    "Minimize\n obj: - x - y + [ 2 x ^ 2 + 2 y ^ 2 ] / 2\nSubject To\n c1: x + y <= 1\nEnd\n";
const SOCP_LP: &str = "Minimize\n obj: - x - y\nSubject To\n q1: [ x ^ 2 + y ^ 2 ] <= 1\nBounds\n x free\n y free\nEnd\n";
const GC_LP: &str = "Maximize\n obj: m\nSubject To\n c1: x + y <= 10\nGeneral Constraints\n g1: m = MAX ( x , y )\nBounds\n 0 <= x <= 8\n 0 <= y <= 8\n m <= 7\nEnd\n";

/// `--no-warnings` is a no-op on every route (QP, conic, general constraints, native).
#[test]
fn no_warnings_is_accepted_on_every_route() {
    let qp = write("nw_qp.lp", QP_LP);
    let socp = write("nw_socp.lp", SOCP_LP);
    let gc = write("nw_gc.lp", GC_LP);
    for (cmd, p) in [
        ("mip", &qp),
        ("solve", &qp),
        ("mip", &socp),
        ("solve", &socp),
        ("mip", &gc),
        ("solve", &gc),
    ] {
        let (code, out, err) = run(&[cmd, p, "--no-warnings", "--format", "json"]);
        assert_eq!(code, 0, "{cmd} {p}: {err}");
        let v: Value = serde_json::from_str(&out).unwrap_or_else(|e| panic!("{cmd}: {e}: {out}"));
        assert!(v.get("status").is_some(), "{cmd}: {out}");
    }
    let lp = write(
        "nw_clean.mps",
        "NAME T\nROWS\n N obj\n L c1\nCOLUMNS\n x obj -1 c1 1\nRHS\n RHS c1 4\nBOUNDS\n UP BND x 10\nENDATA\n",
    );
    // the native search (RUSOLVE_MIP_HIGHS=0) takes it as a no-op
    let o = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["mip", &lp, "--no-warnings"])
        .env("RUSOLVE_MIP_HIGHS", "0")
        .output()
        .unwrap();
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
}

/// n x n integer transport model: supply rows `<= 10`, demand rows `>= 10.1`. Infeasible, and
/// the infeasible set takes many reduction steps.
fn transport_lp(n: usize) -> String {
    let mut s = String::from("Minimize\n obj:");
    for i in 0..n * n {
        s.push_str(&format!(" + x{i}"));
    }
    s.push_str("\nSubject To\n");
    for a in 0..n {
        s.push_str(&format!(" s{a}:"));
        for d in 0..n {
            s.push_str(&format!(" + x{}", a * n + d));
        }
        s.push_str(" <= 10\n");
    }
    for d in 0..n {
        s.push_str(&format!(" d{d}:"));
        for a in 0..n {
            s.push_str(&format!(" + x{}", a * n + d));
        }
        s.push_str(" >= 10.1\n");
    }
    s.push_str("General\n");
    for i in 0..n * n {
        s.push_str(&format!(" x{i}"));
    }
    s.push_str("\nEnd\n");
    s
}

/// A proved status stays proved under `--watchdog` even when the post-solve diagnostics
/// run past the hard deadline; the status is claimed before they start.
#[test]
fn infeasible_stays_infeasible_under_the_watchdog() {
    let p = write("tr250i.lp", &transport_lp(250));
    let t = std::time::Instant::now();
    let (code, out, err) = run(&[
        "mip",
        &p,
        "--format",
        "json",
        "--watchdog",
        "--time-limit",
        "8",
        "--watchdog-grace",
        "2",
    ]);
    let wall = t.elapsed().as_secs_f64();
    assert_eq!(code, 0, "{err}");
    let v: Value = serde_json::from_str(&out).unwrap_or_else(|e| panic!("{e}: {out}"));
    // the solve takes about 1 s on a typical machine; the 8 s limit leaves room for a loaded test machine
    // (with 3 s, a 3x slowdown made the solve itself time out), while the old bug, a post-solve
    // that ignored the limit, overran by about 20 s and still fails both checks
    assert_eq!(v["status"], "Infeasible", "{out}");
    // time limit + grace + slack: the post-solve never runs past the user's limit
    assert!(wall < 12.0, "wall {wall}");
}

/// An infeasible general-constraint model must not be blamed on integrality (it has no
/// integer columns); the conflict is the MAX constraint against the bounds.
#[test]
fn infeasible_gc_model_is_not_blamed_on_integrality() {
    let p = write(
        "gcinf.lp",
        "Maximize\n obj: m\nSubject To\n c1: x + y <= 20\nGeneral Constraints\n g1: m = MAX ( x , y )\nBounds\n 8 <= x <= 9\n 0 <= y <= 8\n m <= 7\nEnd\n",
    );
    for cmd in ["mip", "solve"] {
        let (_, out, _) = run(&[cmd, &p, "--format", "json"]);
        let v: Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["status"], "Infeasible", "{cmd}");
        let text = v["diagnostics"].to_string();
        assert!(!text.contains("integrality"), "{cmd}: {text}");
    }
}
