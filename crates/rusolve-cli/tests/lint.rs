//! `rusolve lint`. One crafted model per rule, each of which must produce that rule and no
//! other; clean fixtures must produce no errors; the JSON must carry every field the committed
//! schema requires. Every test runs the real binary.

use std::path::PathBuf;
use std::process::Command;

const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-lint-1.schema.json");

fn write_tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-lint-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

/// `(exit, stdout, stderr)` of `rusolve lint <model> <extra>`.
fn lint(model: &PathBuf, extra: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("lint")
        .arg(model)
        .args(extra)
        .output()
        .unwrap();
    (
        out.status.code().unwrap(),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8(out.stderr).unwrap(),
    )
}

fn lint_json(model: &PathBuf, extra: &[&str]) -> serde_json::Value {
    let mut args = vec!["--format", "json"];
    args.extend_from_slice(extra);
    let (code, out, err) = lint(model, &args);
    assert_eq!(code, 0, "stderr: {err}");
    serde_json::from_str(out.trim()).unwrap_or_else(|e| panic!("not one JSON object ({e}): {out}"))
}

/// The distinct rule ids of the findings, in report order.
fn rules(v: &serde_json::Value) -> Vec<String> {
    let mut r: Vec<String> = Vec::new();
    for f in v["findings"].as_array().unwrap() {
        let id = f["rule"].as_str().unwrap().to_string();
        if !r.contains(&id) {
            r.push(id);
        }
    }
    r
}

/// Asserts the model yields exactly the one rule, and returns that rule's first finding.
fn only(name: &str, mps: &str, rule: &str, extra: &[&str]) -> serde_json::Value {
    let p = write_tmp(&format!("{name}.mps"), mps);
    let v = lint_json(&p, extra);
    assert_eq!(rules(&v), vec![rule.to_string()], "findings: {v}");
    v["findings"][0].clone()
}

/// A plain model with no problem a rule would report: a base for the crafted ones.
const CLEAN: &str = "\
NAME CLEAN
ROWS
 N  cost
 L  r1
 L  r2
COLUMNS
    x  cost  -1  r1  1
    x  r2  1
    y  cost  -2  r1  1
    y  r2  -1
RHS
    rhs  r1  5  r2  2
ENDATA
";

#[test]
fn lint_clean_model_has_no_findings() {
    let p = write_tmp("clean.mps", CLEAN);
    let v = lint_json(&p, &[]);
    assert!(v["findings"].as_array().unwrap().is_empty(), "{v}");
    assert_eq!(v["summary"]["errors"], 0);
}

#[test]
fn lint_big_m() {
    // y <= 10, link: y - 100 z <= 0 with z binary. The row's maximal activity is 10 (z at its
    // lower bound), so 100 can be 10: the tightest valid big-M.
    let f = only(
        "bigm",
        "\
NAME BIGM
ROWS
 N  cost
 L  link
 L  cap
COLUMNS
    MARKER  'MARKER'  'INTORG'
    z  cost  3  link  -100
    z  cap  2
    MARKER  'MARKER'  'INTEND'
    y  cost  -1  link  1
    y  cap  1
RHS
    rhs  cap  8
BOUNDS
 UP bnd z 1
 UP bnd y 10
ENDATA
",
        "big-m",
        &[],
    );
    assert_eq!(f["severity"], "warning");
    let it = &f["items"][0];
    assert_eq!(it["rows"][0], "link");
    assert_eq!(it["columns"][0], "z");
    assert_eq!(it["evidence"]["coefficient"], -100.0);
    assert_eq!(it["evidence"]["tightest_coefficient"], -10.0);
}

#[test]
fn lint_big_m_no_false_positive() {
    // y - 5 z >= 0 with y <= 10: -5 is within the slack 10, not far above it.
    let p = write_tmp(
        "bigm_ok.mps",
        "\
NAME BIGMOK
ROWS
 N  cost
 G  link
COLUMNS
    MARKER  'MARKER'  'INTORG'
    z  cost  3  link  -5
    MARKER  'MARKER'  'INTEND'
    y  cost  1  link  1
RHS
    rhs  link  0
BOUNDS
 UP bnd z 1
 UP bnd y 10
ENDATA
",
    );
    let v = lint_json(&p, &["--rule", "big-m"]);
    assert!(rules(&v).is_empty(), "{v}");
}

const SINGLETON_BIGM: &str = "NAME SING\nROWS\n N  cost\n L  c1\n L  c2\nCOLUMNS\n    x  cost  -1  c1  1\n    x  c2  1\n    MARKER  'MARKER'  'INTORG'\n    b  c1  -1000000\n    MARKER  'MARKER'  'INTEND'\nRHS\n    rhs  c1  0\n    rhs  c2  5\nBOUNDS\n UP bnd b 1\nENDATA\n";

#[test]
fn lint_big_m_uses_singleton_rows() {
    // x - 1e6 b <= 0 with x <= 5 only stated as a row: the tightest valid coefficient is -5.
    let p = write_tmp("bigm_sing.mps", SINGLETON_BIGM);
    let v = lint_json(&p, &["--rule", "big-m"]);
    let it = &v["findings"][0]["items"][0];
    assert_eq!(it["rows"][0], "c1");
    assert_eq!(it["columns"][0], "b");
    assert_eq!(it["evidence"]["coefficient"], -1000000.0);
    assert_eq!(it["evidence"]["tightest_coefficient"], -5.0);
}

#[test]
fn lint_big_m_non_binding_singleton_row_is_no_finding() {
    // x <= 1e9 is no tighter than the coefficient: no finding.
    let p = write_tmp(
        "bigm_sing_loose.mps",
        &SINGLETON_BIGM.replace(" 5\n", " 1000000000\n"),
    );
    let v = lint_json(&p, &["--rule", "big-m"]);
    assert!(rules(&v).is_empty(), "{v}");
}

#[test]
fn lint_coef_range() {
    let f = only(
        "coef",
        "\
NAME COEF
ROWS
 N  cost
 L  r1
 L  r2
COLUMNS
    x  cost  -1  r1  0.0001
    x  r2  1
    y  cost  -1  r1  1000
    y  r2  1
RHS
    rhs  r1  5  r2  5
ENDATA
",
        "coef-range",
        &[],
    );
    assert_eq!(f["severity"], "info");
    assert_eq!(f["items"][0]["evidence"]["ratio"], 1e7);
}

#[test]
fn lint_rhs_range() {
    let f = only(
        "rhs",
        "\
NAME RHS
ROWS
 N  cost
 L  r1
 L  r2
COLUMNS
    x  cost  -1  r1  1
    x  r2  1
    y  cost  -1  r1  1
    y  r2  -1
RHS
    rhs  r1  2000000  r2  3
ENDATA
",
        "rhs-range",
        &[],
    );
    assert_eq!(f["severity"], "warning");
    assert_eq!(f["items"][0]["evidence"]["max_abs"], 2e6);
}

#[test]
fn lint_obj_range() {
    let f = only(
        "obj",
        &CLEAN.replace("x  cost  -1", "x  cost  -0.00001"),
        "obj-range",
        &[],
    );
    assert_eq!(f["items"][0]["columns"][0], "x");
}

#[test]
fn lint_bound_range() {
    let f = only(
        "bound",
        &CLEAN.replace("ENDATA", "BOUNDS\n UP bnd x 5000000\nENDATA"),
        "bound-range",
        &[],
    );
    assert_eq!(f["severity"], "warning");
    assert_eq!(f["items"][0]["columns"][0], "x");
}

#[test]
fn lint_scaling() {
    // r3 has coefficients about 1e4 times those of the other rows; columns stay level.
    let f = only(
        "scale",
        "\
NAME SCALE
ROWS
 N  cost
 L  r1
 L  r2
 L  r3
COLUMNS
    x  cost  -1  r1  1
    x  r2  1  r3  10000
    y  cost  -1  r1  1
    y  r2  -1  r3  20000
RHS
    rhs  r1  5  r2  2  r3  100000
ENDATA
",
        "scaling",
        &[],
    );
    assert_eq!(f["items"][0]["rows"][0], "r3");
}

#[test]
fn lint_near_parallel_rows() {
    let f = only(
        "par",
        "\
NAME PAR
ROWS
 N  cost
 L  r1
 L  r2
 L  r3
COLUMNS
    x  cost  -1  r1  1
    x  r2  3  r3  1
    y  cost  -1  r1  2
    y  r2  6  r3  -1
RHS
    rhs  r1  4  r2  9  r3  1
ENDATA
",
        "near-parallel-rows",
        &[],
    );
    assert_eq!(f["severity"], "warning");
    let mut names: Vec<&str> = f["items"][0]["rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|n| n.as_str().unwrap())
        .collect();
    names.sort();
    assert_eq!(names, ["r1", "r2"]);
}

#[test]
fn lint_near_parallel_rows_contradiction_is_an_error() {
    // x + 2y <= 4 and 3x + 6y >= 15 (that is x + 2y >= 5): cannot both hold.
    let p = write_tmp(
        "par_infeasible.mps",
        "\
NAME PARINF
ROWS
 N  cost
 L  r1
 G  r2
 L  r3
COLUMNS
    x  cost  -1  r1  1
    x  r2  3  r3  1
    y  cost  -1  r1  2
    y  r2  6  r3  -1
RHS
    rhs  r1  4  r2  15  r3  1
ENDATA
",
    );
    let v = lint_json(&p, &[]);
    assert_eq!(v["summary"]["errors"], 1, "{v}");
    assert_eq!(v["findings"][0]["rule"], "near-parallel-rows");
    assert_eq!(v["findings"][0]["severity"], "error");
}

#[test]
fn lint_empty_row() {
    let f = only(
        "emptyrow",
        &CLEAN.replace(" L  r2\n", " L  r2\n L  e1\n"),
        "empty-row",
        &[],
    );
    assert_eq!(f["items"][0]["rows"][0], "e1");
    assert_eq!(f["severity"], "warning");
}

#[test]
fn lint_empty_row_infeasible_is_an_error() {
    let p = write_tmp(
        "emptyrow_inf.mps",
        &CLEAN
            .replace(" L  r2\n", " L  r2\n G  e1\n")
            .replace("rhs  r1  5  r2  2", "rhs  r1  5  r2  2  e1  3"),
    );
    let v = lint_json(&p, &[]);
    assert_eq!(v["findings"][0]["severity"], "error", "{v}");
}

#[test]
fn lint_empty_col() {
    let f = only(
        "emptycol",
        &CLEAN.replace("RHS", "    z  cost  0  r1  0\nRHS"),
        "empty-col",
        &[],
    );
    assert_eq!(f["items"][0]["columns"][0], "z");
}

#[test]
fn lint_free_row() {
    // a >= row with a -1e30 (infinite) bound has no finite side
    let f = only(
        "freerow",
        &CLEAN
            .replace(" L  r2\n", " L  r2\n G  f1\n")
            .replace("    y  r2  -1", "    y  r2  -1\n    x  f1  1\n    y  f1  2")
            .replace("rhs  r1  5  r2  2", "rhs  r1  5  r2  2  f1  -1e30"),
        "free-row",
        &[],
    );
    assert_eq!(f["items"][0]["rows"][0], "f1");
}

#[test]
fn lint_free_col_no_cost() {
    // w is free, costless, and appears only in r3: r3 can always be satisfied through w
    let f = only(
        "freecol",
        &CLEAN
            .replace(" L  r2\n", " L  r2\n L  r3\n")
            .replace("RHS", "    w  r3  1\n    x  r3  1\n    y  r3  2\nRHS")
            .replace("rhs  r1  5  r2  2", "rhs  r1  5  r2  2  r3  7")
            .replace("ENDATA", "BOUNDS\n FR bnd w\nENDATA"),
        "free-col-no-cost",
        &[],
    );
    assert_eq!(f["items"][0]["columns"][0], "w");
}

#[test]
fn lint_fixed_col() {
    let f = only(
        "fixed",
        &CLEAN.replace("ENDATA", "BOUNDS\n FX bnd y 2\nENDATA"),
        "fixed-col",
        &[],
    );
    assert_eq!(f["severity"], "info");
    assert_eq!(f["items"][0]["evidence"]["value"], 2.0);
}

#[test]
fn lint_fractional_int_bound() {
    let f = only(
        "fracint",
        &CLEAN
            .replace(
                "    x  cost  -1  r1  1",
                "    MARKER  'MARKER'  'INTORG'\n    x  cost  -1  r1  1",
            )
            .replace(
                "    y  cost  -2",
                "    MARKER  'MARKER'  'INTEND'\n    y  cost  -2",
            )
            .replace("ENDATA", "BOUNDS\n UP bnd x 4.5\nENDATA"),
        "fractional-int-bound",
        &[],
    );
    assert_eq!(f["items"][0]["evidence"]["rounded_upper"], 4.0);
}

#[test]
fn lint_condition_estimate_is_opt_in() {
    // x + y = 2 and x + (1+1e-6) y = 2+1e-6: the solution (1, 1) needs both columns basic (the
    // point (0, 2) violates r2 by 1e-6, above the feasibility tolerance), and the basis
    // B = [1 1; 1 1+1e-6] has condition about 4e6. The rows differ by more than the parallel
    // tolerance, so no other rule fires.
    let mps = "\
NAME COND
ROWS
 N  cost
 E  r1
 E  r2
COLUMNS
    x  cost  1  r1  1
    x  r2  1
    y  cost  1  r1  1
    y  r2  1.000001
RHS
    rhs  r1  2  r2  2.000001
ENDATA
";
    let p = write_tmp("cond.mps", mps);
    let without = lint_json(&p, &[]);
    assert!(rules(&without).is_empty(), "{without}");
    let f = only("cond2", mps, "condition-estimate", &["--with-solve"]);
    let kappa = f["items"][0]["evidence"]["kappa_estimate"]
        .as_f64()
        .unwrap();
    assert!(kappa > 1e6, "{kappa}");
}

#[test]
fn lint_rule_filter_and_min_severity() {
    let p = write_tmp(
        "filter.mps",
        &CLEAN.replace("ENDATA", "BOUNDS\n FX bnd y 2\n UP bnd x 5000000\nENDATA"),
    );
    let all = lint_json(&p, &[]);
    assert_eq!(rules(&all), ["bound-range", "fixed-col"]);
    let one = lint_json(&p, &["--rule", "fixed-col"]);
    assert_eq!(rules(&one), ["fixed-col"]);
    let warn = lint_json(&p, &["--min-severity", "warning"]);
    assert_eq!(rules(&warn), ["bound-range"]);
    let (code, _, err) = lint(&p, &["--rule", "nope"]);
    assert_eq!(code, 2);
    assert!(err.contains("unknown rule nope"), "{err}");
}

#[test]
fn lint_text_output_names_rows_and_columns() {
    let p = write_tmp(
        "text.mps",
        &CLEAN.replace("ENDATA", "BOUNDS\n FX bnd y 2\nENDATA"),
    );
    let (code, out, _) = lint(&p, &[]);
    assert_eq!(code, 0);
    assert!(out.contains("[info] fixed-col (1)"), "{out}");
    assert!(out.contains("cols y"), "{out}");
    assert!(out.contains("fix:"), "{out}");
}

#[test]
fn lint_exit_codes() {
    let missing = PathBuf::from("/nonexistent/model.mps");
    assert_eq!(lint(&missing, &[]).0, 3);
    let p = write_tmp("usage.mps", CLEAN);
    assert_eq!(lint(&p, &["--min-severity", "loud"]).0, 2);
    assert_eq!(lint(&p, &["--bogus"]).0, 2);
}

/// Every fixture MPS under rusolve-ekk/tests/fixtures lints without an error.
#[test]
fn lint_fixtures_have_no_errors() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../rusolve-ekk/tests/fixtures");
    let mut n = 0;
    for dir in std::fs::read_dir(&root).unwrap() {
        let dir = dir.unwrap().path();
        if !dir.is_dir() {
            continue;
        }
        for f in std::fs::read_dir(&dir).unwrap() {
            let f = f.unwrap().path();
            if f.extension().is_some_and(|e| e == "mps") {
                let v = lint_json(&f, &[]);
                assert_eq!(v["summary"]["errors"], 0, "{}: {v}", f.display());
                n += 1;
            }
        }
    }
    assert!(n >= 6, "only {n} fixtures found");
}

#[test]
fn lint_json_carries_every_required_field_of_the_schema() {
    let schema: serde_json::Value = serde_json::from_str(SCHEMA).unwrap();
    let required = |node: &serde_json::Value| -> Vec<String> {
        node["required"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect()
    };
    let p = write_tmp(
        "schema.mps",
        &CLEAN.replace("ENDATA", "BOUNDS\n FX bnd y 2\nENDATA"),
    );
    let v = lint_json(&p, &[]);
    for k in required(&schema) {
        assert!(v.get(&k).is_some(), "top-level {k} missing: {v}");
    }
    for k in required(&schema["properties"]["input"]) {
        assert!(v["input"].get(&k).is_some(), "input.{k} missing");
    }
    for k in required(&schema["properties"]["summary"]) {
        assert!(v["summary"].get(&k).is_some(), "summary.{k} missing");
    }
    let fschema = &schema["properties"]["findings"]["items"];
    let finding = &v["findings"][0];
    for k in required(fschema) {
        assert!(finding.get(&k).is_some(), "finding.{k} missing: {finding}");
    }
    for k in required(&fschema["properties"]["items"]["items"]) {
        assert!(finding["items"][0].get(&k).is_some(), "item.{k} missing");
    }
    assert_eq!(v["schema_version"], "rusolve.lint/1");
}
