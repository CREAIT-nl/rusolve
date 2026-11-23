//! Warm-start inputs: `--start`, `--basis`, `--add-rows` on `rusolve mip` (the port path).
//! Runs the real binary, as an agent would.

use std::path::PathBuf;
use std::process::Command;

const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-result-1.schema.json");

/// min -3x - 4y + 5z + w, x + 2y + z <= 8, z + 3w <= 10, x in [0,4], y in [0,5], z in [0,3]
/// integer, w continuous: optimum x=4, y=2, z=0, w=0, objective -20.
const TINY_MIP: &str = "\
NAME WARM1
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

/// min -3x - 2y, c1: x + y <= 4, c2: x + 3y <= 6: optimum x=4, y=0, objective -12.
const TINY_LP: &str = "\
NAME WARMLP
ROWS
 N  cost
 L  c1
 L  c2
COLUMNS
    x  cost  -3  c1  1  c2  1
    y  cost  -2  c1  1  c2  3
RHS
    rhs  c1  4  c2  6
ENDATA
";

fn dir() -> PathBuf {
    let d = std::env::temp_dir().join(format!("rusolve-cli-warm-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write(name: &str, text: &str) -> PathBuf {
    let p = dir().join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn mip(model: &PathBuf, extra: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
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

fn envelope(model: &PathBuf, extra: &[&str]) -> serde_json::Value {
    let (code, stdout, stderr) = mip(model, &[&["--format", "json"][..], extra].concat());
    assert_eq!(code, 0, "exit code, stderr: {stderr}");
    serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not JSON: {e}\n{stdout}"))
}

fn required(section: &str) -> Vec<String> {
    let schema: serde_json::Value = serde_json::from_str(SCHEMA).unwrap();
    schema["properties"][section]["required"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect()
}

fn has_required(env: &serde_json::Value, section: &str) {
    for f in required(section) {
        assert!(
            env[section].get(&f).is_some(),
            "{section}.{f} missing: {env}"
        );
    }
}

#[test]
fn warm_start_roundtrip_keeps_names_and_values_in_both_formats() {
    let model = write("rt.mps", TINY_MIP);
    for ext in ["sol", "json"] {
        let sol = dir().join(format!("rt-out.{ext}"));
        let sol_s = sol.to_str().unwrap();
        let first = envelope(&model, &["--solution-file", sol_s]);
        assert_eq!(first["objective"], -20.0);
        // the file carries the model's names
        let text = std::fs::read_to_string(&sol).unwrap();
        for n in ["x", "y", "z", "w"] {
            assert!(text.contains(&format!("\"{n}\"")) || text.contains(&format!("\n{n} ")));
        }
        let env = envelope(&model, &["--start", sol_s]);
        has_required(&env, "start");
        assert_eq!(env["start"]["accepted"], true);
        assert_eq!(env["start"]["complete"], true);
        assert_eq!(env["start"]["columns_given"], 4);
        assert_eq!(env["start"]["objective"], -20.0);
        assert_eq!(env["start"]["became_incumbent"], true);
        assert_eq!(env["start"]["improved_on_start"], false);
        assert_eq!(env["start"]["final_incumbent_is_start"], true);
        assert_eq!(env["objective"], -20.0);
    }
}

#[test]
fn warm_start_feasible_but_poor_start_is_improved() {
    let model = write("poor.mps", TINY_MIP);
    let start = write(
        "poor.json",
        r#"{"columns":[{"name":"x","value":0},{"name":"y","value":0},{"name":"z","value":0},{"name":"w","value":0}]}"#,
    );
    let env = envelope(&model, &["--start", start.to_str().unwrap()]);
    assert_eq!(env["start"]["accepted"], true);
    assert_eq!(env["start"]["objective"], 0.0);
    assert_eq!(env["start"]["became_incumbent"], true);
    assert_eq!(env["start"]["improved_on_start"], true);
    assert_eq!(env["start"]["final_incumbent_is_start"], false);
    assert_eq!(env["objective"], -20.0);
}

#[test]
fn warm_start_partial_start_is_accepted_and_completed() {
    let model = write("part.mps", TINY_MIP);
    let start = write("part.json", r#"{"columns":[{"name":"x","value":4}]}"#);
    let env = envelope(&model, &["--start", start.to_str().unwrap()]);
    has_required(&env, "start");
    assert_eq!(env["start"]["complete"], false);
    assert_eq!(env["start"]["columns_given"], 1);
    assert!(env["start"]["objective"].is_null());
    assert_eq!(env["objective"], -20.0);
}

#[test]
fn warm_start_infeasible_start_is_an_actionable_input_error() {
    let model = write("inf.mps", TINY_MIP);
    // x + 2y + z = 4 + 10 = 14 > 8 breaks cap1
    let start = write(
        "inf.json",
        r#"{"columns":[{"name":"x","value":4},{"name":"y","value":5},{"name":"z","value":0},{"name":"w","value":0}]}"#,
    );
    let (code, stdout, stderr) = mip(
        &model,
        &["--format", "json", "--start", start.to_str().unwrap()],
    );
    assert_eq!(code, 3, "{stderr}");
    // No result envelope, one rusolve.diagnostics/1 document
    let doc: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(doc["schema_version"], "rusolve.diagnostics/1", "{stdout}");
    assert!(
        doc.get("status").is_none(),
        "no result on an input error: {stdout}"
    );
    assert!(stderr.contains("infeasible"), "{stderr}");
    assert!(stderr.contains("cap1"), "names the violated row: {stderr}");
    assert!(
        stderr.contains("rusolve check"),
        "says what to do: {stderr}"
    );
}

#[test]
fn warm_start_misnamed_column_gets_the_closest_names() {
    let model = write("mis.mps", TINY_MIP);
    let start = write("mis.json", r#"{"columns":[{"name":"xx","value":4}]}"#);
    let (code, _o, stderr) = mip(&model, &["--start", start.to_str().unwrap()]);
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("xx"), "{stderr}");
    assert!(stderr.contains("closest model names"), "{stderr}");
    let (code, _o, _e) = mip(
        &model,
        &["--start", dir().join("nope.sol").to_str().unwrap()],
    );
    assert_eq!(code, 3);
    let (code, _o, _e) = mip(&model, &["--start"]);
    assert_eq!(code, 2);
}

fn basis_file(c1: i32, c2: i32) -> String {
    format!("HiGHS_basis_file v2\nValid\n# Columns 2\nx 1\ny 0\n# Rows 2\nc1 {c1}\nc2 {c2}\n")
}

#[test]
fn warm_basis_lp_warm_start_is_read_by_name_and_used() {
    let model = write("lp.mps", TINY_LP);
    let cold = envelope(&model, &[]);
    assert_eq!(cold["objective"], -12.0);
    // the optimal basis: x basic, c2 basic (slack), c1 at its upper bound
    let f = write("lp.bas", &basis_file(2, 1));
    let env = envelope(&model, &["--basis", f.to_str().unwrap()]);
    has_required(&env, "basis");
    assert_eq!(env["basis"]["installed"], true);
    assert_eq!(env["objective"], -12.0);
    assert_eq!(
        env["stats"]["lp_iterations"], 0,
        "the optimal basis needs no simplex iteration: {env}"
    );
    // a wrong size is an input error naming the section
    let bad = write(
        "lp-bad.bas",
        "HiGHS_basis_file v2\nValid\n# Columns 3\nx 1\ny 0\nq 0\n# Rows 2\nc1 2\nc2 1\n",
    );
    let (code, _o, stderr) = mip(&model, &["--basis", bad.to_str().unwrap()]);
    assert_eq!(code, 3);
    assert!(stderr.contains("3 columns, not 2"), "{stderr}");
    // an unknown row name
    let bad = write(
        "lp-bad2.bas",
        "HiGHS_basis_file v2\nValid\n# Columns 2\nx 1\ny 0\n# Rows 2\nc1 2\nc9 1\n",
    );
    let (code, _o, stderr) = mip(&model, &["--basis", bad.to_str().unwrap()]);
    assert_eq!(code, 3);
    assert!(stderr.contains("\"c9\" is not found"), "{stderr}");
}

#[test]
fn warm_basis_on_a_mip_is_read_but_reported_unused() {
    let model = write("bm.mps", TINY_MIP);
    let f = write(
        "bm.bas",
        "HiGHS_basis_file v2\nValid\n# Columns 4\nx 1\ny 1\nz 0\nw 0\n# Rows 2\ncap1 0\ncap2 1\n",
    );
    let env = envelope(&model, &["--basis", f.to_str().unwrap()]);
    assert_eq!(env["basis"]["installed"], false);
    assert_eq!(env["objective"], -20.0);
}

#[test]
fn add_rows_valid_rows_pass_the_gate_and_malformed_ones_are_rejected_by_name() {
    let model = write("ar.mps", TINY_MIP);
    let frag = write(
        "ar.json",
        r#"{"rows":[
            {"name":"cutA","upper":6,"coefficients":{"x":1,"y":1}},
            {"name":"bad_col","upper":6,"coefficients":{"xx":1}},
            {"name":"bad_empty","upper":6,"coefficients":{}},
            {"name":"cap1","upper":6,"coefficients":{"x":1}}
        ]}"#,
    );
    let env = envelope(&model, &["--add-rows", frag.to_str().unwrap()]);
    has_required(&env, "added_rows");
    let a = &env["added_rows"];
    assert_eq!(a["requested"], 4);
    assert_eq!(a["added"], 1);
    assert_eq!(a["rejected"], 3);
    assert_eq!(a["accepted_names"][0], "cutA");
    let names: Vec<&str> = a["rejected_rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| r["name"].as_str().unwrap())
        .collect();
    assert_eq!(names, ["bad_col", "bad_empty", "cap1"]);
    assert!(
        a["rejected_rows"][0]["reason"]
            .as_str()
            .unwrap()
            .contains("closest model names")
    );
    assert_eq!(env["input"]["rows"], 3, "the model plus the added row");
    assert_eq!(env["objective"], -20.0, "x + y <= 6 holds at the optimum");
}

#[test]
fn add_rows_mps_fragment_is_read_with_ranges() {
    let model = write("arm.mps", TINY_MIP);
    let frag = write(
        "arm.mpsfrag",
        "ROWS\n G  lo\n E  eq\nCOLUMNS\n    x  lo  1  eq  1\n    y  lo  1  eq  -1\nRHS\n    rhs  lo  1  eq  2\nRANGES\n    rng  lo  10\nENDATA\n",
    );
    let env = envelope(
        &model,
        &["--add-rows", frag.to_str().unwrap(), "--solution"],
    );
    let a = &env["added_rows"];
    assert_eq!(a["format"], "mps-fragment");
    assert_eq!(a["added"], 2);
    // x - y = 2 and 1 <= x + y <= 11: x=4, y=2 satisfies both
    assert_eq!(env["objective"], -20.0);
}

#[test]
fn add_rows_gate_rejects_a_row_that_cuts_off_the_debug_solution() {
    let model = write("gate.mps", TINY_MIP);
    let sol = dir().join("gate-opt.sol");
    envelope(&model, &["--solution-file", sol.to_str().unwrap()]);
    let frag = write(
        "gate.json",
        r#"[{"name":"too_tight","upper":3,"coefficients":{"x":1}},
            {"name":"fine","upper":6,"coefficients":{"x":1,"y":1}}]"#,
    );
    let opt = format!("mip_debug_solution_file={}", sol.to_str().unwrap());
    let env = envelope(
        &model,
        &["--add-rows", frag.to_str().unwrap(), "--opt", &opt],
    );
    let a = &env["added_rows"];
    assert_eq!(a["debug_solution_gate"], true);
    assert_eq!(a["accepted_names"][0], "fine");
    assert_eq!(a["rejected_rows"][0]["name"], "too_tight");
    assert!(
        a["rejected_rows"][0]["reason"]
            .as_str()
            .unwrap()
            .contains("cuts off the debug solution")
    );
    assert_eq!(env["objective"], -20.0, "the optimum survived the gate");
    // without a debug solution the gate cannot prove it invalid: it is added and changes the optimum
    let env = envelope(&model, &["--add-rows", frag.to_str().unwrap()]);
    assert_eq!(env["added_rows"]["debug_solution_gate"], false);
    assert_eq!(env["added_rows"]["added"], 2);
    assert_ne!(env["objective"], -20.0);
}

#[test]
fn warm_start_can_follow_added_rows_and_text_output_reports_them() {
    let model = write("both.mps", TINY_MIP);
    let frag = write(
        "both.json",
        r#"[{"name":"cutB","upper":6,"coefficients":{"x":1,"y":1}}]"#,
    );
    let start = write(
        "both.json.start.json",
        r#"{"columns":[{"name":"x","value":4},{"name":"y","value":2},{"name":"z","value":0},{"name":"w","value":0}]}"#,
    );
    let (code, stdout, stderr) = mip(
        &model,
        &[
            "--add-rows",
            frag.to_str().unwrap(),
            "--start",
            start.to_str().unwrap(),
        ],
    );
    assert_eq!(code, 0, "{stderr}");
    assert!(
        stdout.contains("added_rows:") && stdout.contains("start:"),
        "{stdout}"
    );
    // a start that breaks an added row is refused, naming the added row
    let tight = write(
        "both2.json",
        r#"[{"name":"cutC","upper":5,"coefficients":{"x":1,"y":1}}]"#,
    );
    let (code, _o, stderr) = mip(
        &model,
        &[
            "--add-rows",
            tight.to_str().unwrap(),
            "--start",
            start.to_str().unwrap(),
        ],
    );
    assert_eq!(code, 3);
    assert!(stderr.contains("cutC"), "{stderr}");
}
