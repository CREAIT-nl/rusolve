//! `rusolve feasrelax` on the real binary.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-feasrelax-1.schema.json");

/// min x + 2y with x + y >= 10, 10x <= 30, 10y <= 40: infeasible (x + y <= 7).
/// Relaxing `demand` by 3 costs 3; relaxing a cap row costs 30 (its row is scaled by 10).
const ONE_ROW: &str = "\
NAME ONE
ROWS
 N  cost
 G  demand
 L  capx
 L  capy
COLUMNS
    x  cost  1  demand  1
    x  capx  10
    y  cost  2  demand  1
    y  capy  10
RHS
    rhs  demand  10  capx  30
    rhs  capy  40
ENDATA
";

/// x >= 10 against x <= 9 three times: one big violation (row `big`, by 1)
/// or three small ones (s1..s3, by 1 each).
const SEVERAL: &str = "\
NAME SEVERAL
ROWS
 N  cost
 G  big
 L  s1
 L  s2
 L  s3
COLUMNS
    x  cost  1  big  1
    x  s1  1  s2  1
    x  s3  1
RHS
    rhs  big  10  s1  9
    rhs  s2  9  s3  9
ENDATA
";

const FEASIBLE: &str = "\
NAME FEAS
ROWS
 N  cost
 G  r1
COLUMNS
    x  cost  1  r1  1
RHS
    rhs  r1  2
ENDATA
";

/// Integer x in [0, 3] with the row x >= 4 (and a continuous y that does nothing).
const COLBOUND: &str = "\
NAME COLB
ROWS
 N  cost
 G  need
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  cost  1  need  1
    MARKER  'MARKER'  'INTEND'
    y  cost  1
RHS
    rhs  need  4
BOUNDS
 UP bnd x 3
 UP bnd y 5
ENDATA
";

fn write_tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-feasrelax-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn run(model: &PathBuf, extra: &[&str]) -> (i32, Value, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("feasrelax")
        .arg(model)
        .args(["--format", "json"])
        .args(extra)
        .output()
        .unwrap();
    let stdout = String::from_utf8(out.stdout).unwrap();
    let json = serde_json::from_str(&stdout).unwrap_or(Value::Null);
    (
        out.status.code().unwrap_or(-1),
        json,
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

fn names(doc: &Value) -> Vec<String> {
    doc["relaxed"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| r["name"].as_str().unwrap().to_string())
        .collect()
}

fn close(v: &Value, want: f64) -> bool {
    (v.as_f64().unwrap() - want).abs() < 1e-6
}

#[test]
fn feasrelax_one_row_gives_way_in_both_modes() {
    let model = write_tmp("one.mps", ONE_ROW);
    for mode in ["weighted", "min-count"] {
        let (code, doc, err) = run(&model, &["--mode", mode]);
        assert_eq!(code, 0, "{err}");
        assert_eq!(doc["status"], "Relaxed", "{mode}");
        assert_eq!(names(&doc), ["demand"], "{mode}");
        let r = &doc["relaxed"][0];
        assert_eq!(r["bound"], "lower");
        assert_eq!(r["original_bound"], 10.0);
        assert!(close(&r["violation"], 3.0));
        assert!(close(&r["relaxed_to"], 7.0));
        assert!(close(&r["penalty"], 3.0));
        assert_eq!(
            doc["edits"][0]["text"],
            "change lower bound of row demand from 10 to 7"
        );
        // x + 2y with x + y = 7, 10x <= 30, 10y <= 40: x = 3, y = 4
        assert!(close(&doc["summary"]["objective_original"], 11.0));
        assert_eq!(doc["verified"], true, "{mode}");
    }
}

#[test]
fn feasrelax_weighted_and_min_count_differ() {
    let model = write_tmp("several.mps", SEVERAL);
    let pen = write_tmp("several.json", r#"{"rows":{"big":5}}"#);
    let pen = pen.to_str().unwrap();
    // weighted: three small violations (total 3) beat the weight-5 big one
    let (_, w, _) = run(&model, &["--penalties", pen]);
    assert_eq!(names(&w), ["s1", "s2", "s3"]);
    assert!(close(&w["summary"]["total_penalty"], 3.0));
    // min-count: one relaxed constraint, however heavy
    let (_, c, _) = run(&model, &["--penalties", pen, "--mode", "min-count"]);
    assert_eq!(names(&c), ["big"]);
    assert_eq!(c["summary"]["relaxed_count"], 1);
    assert!(close(&c["summary"]["total_violation"], 1.0));
    assert_eq!(w["verified"], true);
    assert_eq!(c["verified"], true);
}

#[test]
fn feasrelax_verification_applies_the_edits() {
    let model = write_tmp("several_v.mps", SEVERAL);
    let (_, doc, _) = run(&model, &["--mode", "min-count"]);
    assert_eq!(doc["verified"], true);
    // the original model with the edits applied: x = 9, objective 9
    assert!(close(&doc["summary"]["objective_after_edits"], 9.0));
    // the model itself is infeasible as given
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .args(["mip", model.to_str().unwrap(), "--format", "json"])
        .output()
        .unwrap();
    let env: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(env["status"], "Infeasible");
}

#[test]
fn feasrelax_feasible_model_reports_zero_relaxation() {
    let model = write_tmp("feas.mps", FEASIBLE);
    for mode in ["weighted", "min-count"] {
        let (code, doc, _) = run(&model, &["--mode", mode]);
        assert_eq!(code, 0);
        assert_eq!(doc["status"], "Feasible");
        assert_eq!(doc["summary"]["relaxed_count"], 0);
        assert_eq!(doc["relaxed"].as_array().unwrap().len(), 0);
        assert_eq!(doc["edits"].as_array().unwrap().len(), 0);
        assert_eq!(doc["verified"], true);
    }
}

#[test]
fn feasrelax_column_bound_by_name_with_negative_row_weight() {
    let model = write_tmp("colb.mps", COLBOUND);
    for mode in ["weighted", "min-count"] {
        let (code, doc, err) = run(&model, &["--mode", mode, "--global-rhs-penalty", "-1"]);
        assert_eq!(code, 0, "{err}");
        assert_eq!(names(&doc), ["x"], "{mode}");
        assert_eq!(doc["relaxed"][0]["kind"], "column");
        assert_eq!(doc["relaxed"][0]["bound"], "upper");
        assert!(close(&doc["relaxed"][0]["relaxed_to"], 4.0));
        assert_eq!(doc["verified"], true, "{mode}");
    }
}

#[test]
fn feasrelax_not_relaxable_when_every_weight_is_negative() {
    let model = write_tmp("several_n.mps", SEVERAL);
    let all_negative = [
        "--global-lower-penalty",
        "-1",
        "--global-upper-penalty",
        "-1",
        "--global-rhs-penalty",
        "-1",
    ];
    for mode in ["weighted", "min-count"] {
        let (code, doc, _) = run(&model, &[&all_negative[..], &["--mode", mode]].concat());
        assert_eq!(code, 0);
        assert_eq!(doc["status"], "NotRelaxable", "{mode}");
        assert_eq!(doc["verified"], false);
    }
}

#[test]
fn feasrelax_usage_and_input_errors_use_the_a5_codes() {
    let model = write_tmp("one_e.mps", ONE_ROW);
    let (code, _, _) = run(&model, &["--mode", "bogus"]);
    assert_eq!(code, 2);
    let bad_name = write_tmp("badname.json", r#"{"rows":{"nope":1}}"#);
    let (code, _, err) = run(&model, &["--penalties", bad_name.to_str().unwrap()]);
    assert_eq!(code, 2, "{err}");
    let (code, _, _) = run(&PathBuf::from("/nonexistent/x.mps"), &[]);
    assert_eq!(code, 3);
    let junk = write_tmp("junk.json", "{not json");
    let (code, _, _) = run(&model, &["--penalties", junk.to_str().unwrap()]);
    assert_eq!(code, 3);
}

#[test]
fn feasrelax_json_carries_every_required_field_of_the_schema() {
    let schema: Value = serde_json::from_str(SCHEMA).unwrap();
    let required = |node: &Value| -> Vec<String> {
        node["required"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect()
    };
    let model = write_tmp("one_s.mps", ONE_ROW);
    let (_, doc, _) = run(&model, &["--mode", "min-count"]);
    for f in required(&schema) {
        assert!(doc.get(&f).is_some(), "top-level {f} missing");
    }
    for obj in ["input", "weights", "summary", "verification"] {
        for f in required(&schema["properties"][obj]) {
            assert!(doc[obj].get(&f).is_some(), "{obj}.{f} missing");
        }
    }
    for arr in ["relaxed", "edits"] {
        let items = &schema["properties"][arr]["items"];
        assert!(!doc[arr].as_array().unwrap().is_empty());
        for entry in doc[arr].as_array().unwrap() {
            for f in required(items) {
                assert!(entry.get(&f).is_some(), "{arr}[].{f} missing");
            }
        }
    }
    assert_eq!(doc["schema_version"], "rusolve.feasrelax/1");
}

#[test]
fn feasrelax_text_report_names_the_edit() {
    let model = write_tmp("one_t.mps", ONE_ROW);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["feasrelax", model.to_str().unwrap()])
        .output()
        .unwrap();
    let text = String::from_utf8(out.stdout).unwrap();
    assert!(
        text.contains("change lower bound of row demand from 10 to 7"),
        "{text}"
    );
    assert!(text.contains("verified true"), "{text}");
}
