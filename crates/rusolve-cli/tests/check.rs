//! `rusolve check <model> <solution>`, the named solution checker. Every test
//! runs the real binary.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

/// `docs/schema/rusolve-check-1.schema.json`, the normative report schema.
const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-check-1.schema.json");

/// min -3x - 4y + 5z + w - 7 (objective offset: RHS of the N row is -offset)
/// with x + 2y + z <= 8 and z + 3w <= 10, x in [0,4], y in [0,5], z in [0,3]
/// integer and w continuous. Optimum x=4, y=2, z=0, w=0: objective -27.
const MODEL: &str = "\
NAME CHK1
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
    rhs  cost  7  cap1  8  cap2  10
BOUNDS
 UP bnd x 4
 UP bnd y 5
 UP bnd z 3
ENDATA
";

fn tmp_dir() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-check-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Writes `text` once: tests run in parallel and share file names, so a rewrite
/// would truncate a file another test is reading. A name must keep one content.
fn write_tmp(name: &str, text: &str) -> PathBuf {
    let p = tmp_dir().join(name);
    if std::fs::read_to_string(&p).ok().as_deref() != Some(text) {
        let tmp = tmp_dir().join(format!("{name}.{:?}.part", std::thread::current().id()));
        std::fs::write(&tmp, text).unwrap();
        std::fs::rename(&tmp, &p).unwrap();
    }
    p
}

/// A HiGHS raw solution file for the four columns, as `--solution-file` writes it.
fn sol_text(x: f64, y: f64, z: f64, w: f64) -> String {
    format!(
        "Model status\nOptimal\n\n# Primal solution values\nFeasible\nObjective 0\n\
         # Columns 4\nx {x}\ny {y}\nz {z}\nw {w}\n# Rows 2\ncap1 0\ncap2 0\n\n\
         # Dual solution values\nNone\n\n# Basis\nHiGHS_basis_file 2\nNone\n"
    )
}

fn run(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("check")
        .args(args)
        .output()
        .unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// `rusolve check --format json` on `MODEL` and a solution text, parsed.
fn check_json(sol: &str) -> Value {
    let model = write_tmp("m.mps", MODEL);
    // one file per call: the tests run in parallel
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    let sol = write_tmp(&format!("s{n}.sol"), sol);
    let (code, stdout, stderr) = run(&[
        model.to_str().unwrap(),
        sol.to_str().unwrap(),
        "--format",
        "json",
    ]);
    assert_eq!(code, 0, "stderr: {stderr}");
    serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not one JSON object: {e}\n{stdout}"))
}

fn finding<'a>(report: &'a Value, kind: &str, name: &str) -> &'a Value {
    report["findings"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["kind"] == kind && f["name"] == name)
        .unwrap_or_else(|| panic!("no {kind} finding for {name} in {report}"))
}

#[test]
fn check_planted_column_bound_is_found_by_name() {
    // x = -1 breaks its lower bound 0 by 1 and nothing else
    let r = check_json(&sol_text(-1.0, 2.0, 0.0, 0.0));
    assert_eq!(r["verdict"], "infeasible");
    assert_eq!(r["findings_total"], 1);
    let f = finding(&r, "column_bound", "x");
    assert_eq!(f["amount"], 1.0);
    assert_eq!(f["value"], -1.0);
    assert_eq!(f["lower"], 0.0);
    assert_eq!(f["upper"], 4.0);
    assert_eq!(f["violated_bound"], 0.0);
    assert_eq!(r["summary"]["column_bounds"]["count"], 1);
    assert_eq!(r["summary"]["rows"]["count"], 0);
    assert_eq!(r["summary"]["integrality"]["count"], 0);
    // objective recomputed with the offset: HiGHS reads an N-row RHS of 7 as the
    // offset -7, so -3*-1 - 4*2 - 7 = -12
    assert_eq!(r["objective"]["value"], -12.0);
    assert_eq!(r["objective"]["offset"], -7.0);
}

#[test]
fn check_planted_row_violation_is_found_by_name() {
    // x=4, y=3: cap1 activity 10 > 8 by 2, all bounds met
    let r = check_json(&sol_text(4.0, 3.0, 0.0, 0.0));
    assert_eq!(r["verdict"], "infeasible");
    assert_eq!(r["findings_total"], 1);
    let f = finding(&r, "row", "cap1");
    assert_eq!(f["amount"], 2.0);
    assert_eq!(f["value"], 10.0);
    assert_eq!(f["upper"], 8.0);
    assert_eq!(f["lower"], Value::Null);
    assert_eq!(f["violated_bound"], 8.0);
    assert_eq!(r["summary"]["rows"]["max"], 2.0);
}

#[test]
fn check_planted_integrality_violation_is_found_by_name() {
    // z = 0.5 is integer-infeasible by 0.5; rows and bounds are fine
    let r = check_json(&sol_text(3.0, 2.0, 0.5, 0.0));
    assert_eq!(r["verdict"], "infeasible");
    assert_eq!(r["findings_total"], 1);
    let f = finding(&r, "integrality", "z");
    assert_eq!(f["amount"], 0.5);
    assert_eq!(f["value"], 0.5);
    assert!(f["nearest_integer"] == 0.0 || f["nearest_integer"] == 1.0);
    assert_eq!(r["summary"]["integrality"]["count"], 1);
    assert_eq!(r["highs"]["integral"], false);
    assert_eq!(r["agrees_with_highs"], true);
}

#[test]
fn check_findings_are_sorted_by_amount_and_cut_at_max_findings() {
    // x = 6 (bound by 2, cap1 by 2 + ...): several findings of different size
    let model = write_tmp("m.mps", MODEL);
    let sol = write_tmp("sorted.sol", &sol_text(6.0, 5.0, 0.25, 0.0));
    let args = |extra: &[&'static str]| {
        let mut a = vec![model.to_str().unwrap(), sol.to_str().unwrap()];
        a.extend_from_slice(extra);
        a.extend_from_slice(&["--format", "json"]);
        run(&a)
    };
    let (code, stdout, _) = args(&[]);
    assert_eq!(code, 0);
    let all: Value = serde_json::from_str(&stdout).unwrap();
    let amounts: Vec<f64> = all["findings"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f["amount"].as_f64().unwrap())
        .collect();
    assert!(amounts.len() >= 3, "{all}");
    assert!(amounts.windows(2).all(|w| w[0] >= w[1]), "{amounts:?}");
    let total = all["findings_total"].as_i64().unwrap();
    let (code, stdout, _) = args(&["--max-findings", "1"]);
    assert_eq!(code, 0);
    let cut: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(cut["findings"].as_array().unwrap().len(), 1);
    assert_eq!(cut["findings_shown"], 1);
    assert_eq!(cut["findings_total"], total, "totals survive the cut");
    assert_eq!(cut["summary"], all["summary"]);
    assert_eq!(cut["findings"][0]["amount"].as_f64().unwrap(), amounts[0]);
}

#[test]
fn check_optimum_from_mip_solution_file_is_feasible_with_the_solves_objective() {
    let model = write_tmp("opt.mps", MODEL);
    let sol = tmp_dir().join("opt.sol");
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(&model)
        .args(["--format", "json", "--solution-file", sol.to_str().unwrap()])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0));
    let solve: Value = serde_json::from_slice(&out.stdout).unwrap();
    let (code, stdout, stderr) = run(&[
        model.to_str().unwrap(),
        sol.to_str().unwrap(),
        "--format",
        "json",
    ]);
    assert_eq!(code, 0, "{stderr}");
    let r: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(r["verdict"], "feasible");
    assert_eq!(r["findings_total"], 0);
    assert_eq!(r["findings"].as_array().unwrap().len(), 0);
    assert_eq!(r["objective"]["value"], solve["objective"]);
    assert_eq!(r["objective"]["value"], -27.0);
    assert_eq!(r["agrees_with_highs"], true);
    assert_eq!(r["highs"]["feasible"], true);
    assert_eq!(r["highs"]["kkt_num_primal_infeasibilities"], 0);

    // the JSON solution file form checks the same
    let json_sol = tmp_dir().join("opt.json");
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(&model)
        .args(["--solution-file", json_sol.to_str().unwrap()])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0));
    let (code, stdout, stderr) = run(&[
        model.to_str().unwrap(),
        json_sol.to_str().unwrap(),
        "--format",
        "json",
    ]);
    assert_eq!(code, 0, "{stderr}");
    let r2: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(r2["verdict"], "feasible");
    assert_eq!(r2["objective"]["value"], -27.0);
}

/// The `required` lists of the committed schema, so the test checks what the
/// file demands (the approach of `tests/envelope.rs`).
fn check_against_schema(report: &Value) {
    let schema: Value = serde_json::from_str(SCHEMA).unwrap();
    fn walk(schema: &Value, root: &Value, node: &Value, at: &str) {
        let Some(obj) = node.as_object() else { return };
        if let Some(req) = schema["required"].as_array() {
            for field in req {
                let field = field.as_str().unwrap();
                assert!(obj.contains_key(field), "{at}.{field} is required: {node}");
            }
        }
        if let Some(props) = schema["properties"].as_object() {
            for (key, value) in obj {
                let sub = props
                    .get(key)
                    .unwrap_or_else(|| panic!("{at}.{key} is not in the schema"));
                let sub = match sub["$ref"].as_str() {
                    Some(r) => &root["$defs"][r.rsplit('/').next().unwrap()],
                    None => sub,
                };
                walk(sub, root, value, &format!("{at}.{key}"));
                if let (Some(items), Some(arr)) = (sub.get("items"), value.as_array()) {
                    for item in arr {
                        walk(items, root, item, &format!("{at}.{key}[]"));
                    }
                }
            }
        }
    }
    walk(&schema, &schema, report, "$");
    assert_eq!(report["schema_version"], "rusolve.check/1");
}

#[test]
fn check_json_matches_its_schema() {
    // a report with every kind of finding, and a clean one
    check_against_schema(&check_json(&sol_text(6.0, 5.0, 0.25, 0.0)));
    check_against_schema(&check_json(&sol_text(4.0, 2.0, 0.0, 0.0)));
}

#[test]
fn check_text_report_names_the_violations() {
    let model = write_tmp("m.mps", MODEL);
    let sol = write_tmp("text.sol", &sol_text(3.0, 3.0, 0.5, 0.0));
    let (code, stdout, _) = run(&[model.to_str().unwrap(), sol.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(stdout.starts_with("verdict: infeasible"), "{stdout}");
    assert!(stdout.contains("cap1"), "{stdout}");
    assert!(stdout.contains("integrality  z"), "{stdout}");
}

#[test]
fn check_missing_column_is_an_input_error_naming_it() {
    let model = write_tmp("m.mps", MODEL);
    let sol = write_tmp(
        "missing.sol",
        "Model status\nOptimal\n\n# Primal solution values\nFeasible\nObjective 0\n\
         # Columns 3\nx 1\ny 1\nz 1\n",
    );
    let (code, stdout, stderr) = run(&[model.to_str().unwrap(), sol.to_str().unwrap()]);
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("column w"), "{stderr}");
    assert!(stderr.contains("missing"), "{stderr}");
    assert!(stdout.is_empty(), "{stdout}");
}

#[test]
fn check_unknown_column_is_an_input_error_naming_it() {
    let model = write_tmp("m.mps", MODEL);
    let sol = write_tmp(
        "unknown.sol",
        "Model status\nOptimal\n\n# Primal solution values\nFeasible\nObjective 0\n\
         # Columns 5\nx 1\ny 1\nz 1\nw 0\nghost 3\n",
    );
    let (code, _, stderr) = run(&[model.to_str().unwrap(), sol.to_str().unwrap()]);
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("ghost"), "{stderr}");
}

#[test]
fn check_json_solution_with_missing_and_unknown_columns_is_an_input_error() {
    let model = write_tmp("m.mps", MODEL);
    let missing = write_tmp(
        "missing.json",
        r#"{"columns":[{"name":"x","value":1},{"name":"y","value":1},{"name":"z","value":1}]}"#,
    );
    let (code, _, stderr) = run(&[model.to_str().unwrap(), missing.to_str().unwrap()]);
    assert_eq!(code, 3);
    assert!(stderr.contains("column w"), "{stderr}");
    let unknown = write_tmp(
        "unknown.json",
        r#"{"columns":[{"name":"x","value":1},{"name":"y","value":1},{"name":"z","value":1},
            {"name":"w","value":1},{"name":"nope","value":1}]}"#,
    );
    let (code, _, stderr) = run(&[model.to_str().unwrap(), unknown.to_str().unwrap()]);
    assert_eq!(code, 3);
    assert!(stderr.contains("nope"), "{stderr}");
}

#[test]
fn check_bad_inputs_and_usage_exit_codes() {
    let model = write_tmp("m.mps", MODEL);
    let sol = write_tmp("ok.sol", &sol_text(4.0, 2.0, 0.0, 0.0));
    // unreadable solution file, unreadable model: input errors naming the path
    let (code, _, stderr) = run(&[model.to_str().unwrap(), "/nonexistent/none.sol"]);
    assert_eq!(code, 3);
    assert!(stderr.contains("/nonexistent/none.sol"), "{stderr}");
    let (code, _, stderr) = run(&["/nonexistent/none.mps", sol.to_str().unwrap()]);
    assert_eq!(code, 3);
    assert!(stderr.contains("/nonexistent/none.mps"), "{stderr}");
    // a value that is not a number names its column
    let bad = write_tmp(
        "nan.sol",
        "Model status\nOptimal\n\n# Primal solution values\nFeasible\nObjective 0\n\
         # Columns 4\nx 1\ny abc\nz 1\nw 0\n",
    );
    let (code, _, stderr) = run(&[model.to_str().unwrap(), bad.to_str().unwrap()]);
    assert_eq!(code, 3);
    assert!(stderr.contains("column y"), "{stderr}");
    // usage errors are 2
    let (code, _, _) = run(&[model.to_str().unwrap()]);
    assert_eq!(code, 2);
    let (code, _, _) = run(&[
        model.to_str().unwrap(),
        sol.to_str().unwrap(),
        "--max-findings",
        "x",
    ]);
    assert_eq!(code, 2);
}

#[test]
fn check_tolerance_flags_move_the_verdict() {
    // z = 1e-4 off an integer: infeasible at the default 1e-6, feasible at 1e-3
    let model = write_tmp("m.mps", MODEL);
    let sol = write_tmp("tol.sol", &sol_text(4.0, 2.0, 0.0001, 0.0));
    let (_, stdout, _) = run(&[
        model.to_str().unwrap(),
        sol.to_str().unwrap(),
        "--format",
        "json",
    ]);
    let r: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(r["verdict"], "infeasible");
    let (_, stdout, _) = run(&[
        model.to_str().unwrap(),
        sol.to_str().unwrap(),
        "--format",
        "json",
        "--mip-feasibility-tolerance",
        "1e-3",
    ]);
    let r: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(r["verdict"], "feasible");
    assert_eq!(r["agrees_with_highs"], true);
}
