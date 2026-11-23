//! Farkas certificates, primal rays and the verified LP IIS:
//! `rusolve explain <model.mps> [--format json]`, through the real binary.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

/// `docs/schema/rusolve-explain-1.schema.json`, the normative output schema.
const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-explain-1.schema.json");

/// min x + 2y + z with cap: x + y <= 3, needx: x >= 2, needy: y >= 2 and a
/// decoy row x + y + z <= 100. The IIS is exactly {cap, needx, needy}.
const INF3: &str = "\
NAME INF3
ROWS
 N  obj
 L  cap
 G  needx
 G  needy
 L  decoy
COLUMNS
    x  obj  1  cap  1
    x  needx  1  decoy  1
    y  obj  2  cap  1
    y  needy  1  decoy  1
    z  obj  1  decoy  1
RHS
    rhs  cap  3  needx  2
    rhs  needy  2  decoy  100
ENDATA
";

/// The same system with x and y declared integer (bounds explicit, since MPS
/// integer columns default to [0, 1]): a MIP whose relaxation is infeasible.
const INF3_MIP: &str = "\
NAME INF3M
ROWS
 N  obj
 L  cap
 G  needx
 G  needy
 L  decoy
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  cap  1
    x  needx  1  decoy  1
    y  obj  2  cap  1
    y  needy  1  decoy  1
    MARKER  'MARKER'  'INTEND'
    z  obj  1  decoy  1
RHS
    rhs  cap  3  needx  2
    rhs  needy  2  decoy  100
BOUNDS
 UP bnd x 10
 UP bnd y 10
ENDATA
";

/// x + y >= 10 with x, y <= 4: infeasible only through the column upper bounds.
const COLB: &str = "\
NAME COLB
ROWS
 N  obj
 G  r1
 L  r2
 L  other
COLUMNS
    x  obj  1  r1  1
    x  other  1
    y  obj  1  r1  1
    y  r2  1
    z  other  1
RHS
    rhs  r1  10  r2  3
    rhs  other  50
BOUNDS
 UP bnd x 4
 UP bnd y 4
ENDATA
";

/// min -x - y with x - y >= -5 and y <= 7: unbounded along x.
const UNB: &str = "\
NAME UNB
ROWS
 N  obj
 G  c1
 L  c2
COLUMNS
    x  obj  -1  c1  1
    y  obj  -1  c1  -1
    y  c2  1
RHS
    rhs  c1  -5  c2  7
ENDATA
";

/// min x + y with x + y >= 2, x, y <= 4: feasible and bounded.
const FEAS: &str = "\
NAME FEAS
ROWS
 N  obj
 G  r1
COLUMNS
    x  obj  1  r1  1
    y  obj  1  r1  1
RHS
    rhs  r1  2
BOUNDS
 UP bnd x 4
 UP bnd y 4
ENDATA
";

/// Integer-infeasible but LP-feasible: 2x = 1 with x integer.
const INT_ONLY: &str = "\
NAME INTONLY
ROWS
 N  obj
 E  e1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  e1  2
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  e1  1
ENDATA
";

fn write_tmp(name: &str, text: &str) -> PathBuf {
    // one directory per test thread: tests run in parallel threads of one process and several
    // write the same file name, and a rewrite truncates the file another thread may be reading;
    // within a test the path stays the same (some tests compare outputs that include it)
    let thread: String = format!("{:?}", std::thread::current().id())
        .chars()
        .filter(char::is_ascii_digit)
        .collect();
    let dir = std::env::temp_dir().join(format!(
        "rusolve-cli-explain-{}-{thread}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn explain_raw(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("explain")
        .args(args)
        .output()
        .unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

fn explain(name: &str, mps: &str, extra: &[&str]) -> Value {
    let p = write_tmp(name, mps);
    let mut args = vec![p.to_str().unwrap(), "--format", "json"];
    args.extend_from_slice(extra);
    let (code, stdout, stderr) = explain_raw(&args);
    assert_eq!(code, 0, "exit code, stderr: {stderr}");
    serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not one JSON object: {e}\n{stdout}"))
}

fn names(v: &Value) -> Vec<String> {
    let mut n: Vec<String> = v
        .as_array()
        .unwrap()
        .iter()
        .map(|m| m["name"].as_str().unwrap().to_string())
        .collect();
    n.sort();
    n
}

#[test]
fn explain_iis_three_rows_by_name_and_verified() {
    let e = explain("inf3.mps", INF3, &[]);
    assert_eq!(e["status"], "Infeasible");
    assert_eq!(e["kind"], "lp");
    assert_eq!(e["verified"], true, "{e}");
    let x = &e["explanation"];
    assert_eq!(x["type"], "infeasible");
    let iis = &x["iis"];
    assert_eq!(iis["kind"], "lp");
    assert_eq!(iis["highs_claim"], "irreducible");
    assert_eq!(
        names(&iis["rows"]),
        ["cap", "needx", "needy"],
        "the decoy row is not in it"
    );
    for r in iis["rows"].as_array().unwrap() {
        let want = if r["name"] == "cap" { "upper" } else { "lower" };
        assert_eq!(r["bound"], want, "{r}");
    }
    // the verification ran: the set is infeasible and every side is needed
    let v = &x["verification"];
    assert_eq!(v["performed"], true);
    assert_eq!(v["set_infeasible"], true);
    assert_eq!(v["irreducible"], true);
    assert_eq!(v["failures"].as_array().unwrap().len(), 0);
    assert_eq!(
        v["checks"].as_array().unwrap().len(),
        3,
        "one drop per member bound: {v}"
    );
    assert!(
        v["checks"]
            .as_array()
            .unwrap()
            .iter()
            .all(|c| c["rest_feasible"] == true)
    );
}

#[test]
fn explain_iis_column_bounds_are_named_with_their_side() {
    let e = explain("colb.mps", COLB, &[]);
    assert_eq!(e["verified"], true, "{e}");
    let iis = &e["explanation"]["iis"];
    assert_eq!(names(&iis["rows"]), ["r1"]);
    assert_eq!(names(&iis["columns"]), ["x", "y"]);
    for c in iis["columns"].as_array().unwrap() {
        assert_eq!(c["bound"], "upper", "{c}");
        assert_eq!(c["upper"], 4.0);
    }
}

#[test]
fn explain_farkas_proof_checks_numerically() {
    let e = explain("inf3.mps", INF3, &[]);
    let f = &e["explanation"]["farkas"];
    assert_eq!(f["available"], true);
    // multipliers by row name, nonzeros only: the decoy row is absent
    let mult: std::collections::BTreeMap<String, f64> = f["multipliers"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| {
            (
                m["row"].as_str().unwrap().to_string(),
                m["multiplier"].as_f64().unwrap(),
            )
        })
        .collect();
    assert_eq!(
        mult.keys().cloned().collect::<Vec<_>>(),
        ["cap", "needx", "needy"]
    );
    assert!(mult.values().all(|v| *v != 0.0));
    // recompute the combination from the model's own rows: cap = x + y (<= 3),
    // needx = x (>= 2), needy = y (>= 2). The combined form on x and y is
    // m_cap*(x+y) + m_needx*x + m_needy*y.
    let cx = mult["cap"] + mult["needx"];
    let cy = mult["cap"] + mult["needy"];
    let (cap_lo, cap_hi) = (f64::NEG_INFINITY, 3.0);
    let (need_lo, need_hi) = (2.0, f64::INFINITY);
    let interval = |m: f64, lo: f64, hi: f64| {
        if m > 0.0 {
            (m * lo, m * hi)
        } else {
            (m * hi, m * lo)
        }
    };
    let parts = [
        interval(mult["cap"], cap_lo, cap_hi),
        interval(mult["needx"], need_lo, need_hi),
        interval(mult["needy"], need_lo, need_hi),
    ];
    let row_lo: f64 = parts.iter().map(|p| p.0).sum();
    let row_hi: f64 = parts.iter().map(|p| p.1).sum();
    // x, y >= 0 with no upper bounds: c'x over the column bounds
    let col_range = |c: f64| {
        if c.abs() < 1e-9 {
            (0.0, 0.0)
        } else if c > 0.0 {
            (0.0, f64::INFINITY)
        } else {
            (f64::NEG_INFINITY, 0.0)
        }
    };
    let (a, b) = (col_range(cx), col_range(cy));
    let (col_lo, col_hi) = (a.0 + b.0, a.1 + b.1);
    let gap = (row_lo - col_hi).max(col_lo - row_hi);
    assert!(
        gap.is_finite() && gap > 1e-9,
        "the multipliers must prove 0 >= gap > 0, got {gap}"
    );
    let p = &f["proof"];
    assert_eq!(p["verified"], true);
    assert!(
        (p["gap"].as_f64().unwrap() - gap).abs() < 1e-9,
        "{p} vs {gap}"
    );
    assert_eq!(
        p["combined_form_nonzero"], 0,
        "x and y cancel in the combination"
    );
}

#[test]
fn explain_unbounded_ray_is_verified_by_name() {
    let e = explain("unb.mps", UNB, &[]);
    assert_eq!(e["status"], "Unbounded");
    assert_eq!(e["verified"], true, "{e}");
    let ray = &e["explanation"]["ray"];
    assert_eq!(ray["available"], true);
    let d: std::collections::BTreeMap<String, f64> = ray["direction"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| {
            (
                m["column"].as_str().unwrap().to_string(),
                m["value"].as_f64().unwrap(),
            )
        })
        .collect();
    // recompute: c1: x - y >= -5 needs d_x - d_y >= 0, c2: y <= 7 needs d_y <= 0,
    // y >= 0 needs d_y >= 0, min -x - y must decrease.
    let (dx, dy) = (*d.get("x").unwrap_or(&0.0), *d.get("y").unwrap_or(&0.0));
    assert!(dx - dy >= -1e-9 && (-1e-9..=1e-9).contains(&dy));
    assert!(
        -dx - dy < -1e-9,
        "the objective must decrease along the ray"
    );
    assert!(ray["objective_rate"].as_f64().unwrap() < 0.0);
    assert_eq!(ray["failures"].as_array().unwrap().len(), 0);
}

#[test]
fn explain_feasible_model_reports_nothing_and_exits_zero() {
    let e = explain("feas.mps", FEAS, &[]);
    assert_eq!(e["status"], "Optimal");
    assert!(e["explanation"].is_null());
    assert!(e["verified"].is_null());
    let p = write_tmp("feas.mps", FEAS);
    let (code, stdout, _) = explain_raw(&[p.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(stdout.contains("nothing to explain"), "{stdout}");
}

#[test]
fn explain_mip_is_labelled_relaxation() {
    let e = explain("inf3m.mps", INF3_MIP, &[]);
    assert_eq!(e["kind"], "relaxation");
    assert!(e["note"].as_str().unwrap().contains("relaxation"));
    assert_eq!(e["explanation"]["iis"]["kind"], "relaxation");
    assert_eq!(
        names(&e["explanation"]["iis"]["rows"]),
        ["cap", "needx", "needy"]
    );
    assert_eq!(e["verified"], true);
    // an integer-only infeasibility has nothing to say about the relaxation
    let i = explain("intonly.mps", INT_ONLY, &[]);
    assert_eq!(i["kind"], "relaxation");
    assert_eq!(i["status"], "Optimal");
    assert!(i["explanation"].is_null());
    assert!(i["note"].as_str().unwrap().contains("true MIP IIS"));
}

#[test]
fn explain_no_verify_and_max_members_options() {
    let e = explain("inf3.mps", INF3, &["--no-verify"]);
    assert_eq!(e["explanation"]["verification"]["performed"], false);
    assert_eq!(
        e["verified"], false,
        "an unverified set must not claim verified"
    );
    let t = explain("inf3.mps", INF3, &["--max-members", "1"]);
    let iis = &t["explanation"]["iis"];
    assert_eq!(iis["rows"].as_array().unwrap().len(), 1);
    assert_eq!(iis["num_rows"], 3);
    assert_eq!(iis["truncated"], true);
    assert_eq!(
        t["explanation"]["farkas"]["multipliers"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        t["verified"], true,
        "truncating the listing does not truncate the check"
    );
}

#[test]
fn explain_exit_codes() {
    assert_eq!(explain_raw(&[]).0, 2, "missing path is a usage error");
    let p = write_tmp("feas.mps", FEAS);
    assert_eq!(explain_raw(&[p.to_str().unwrap(), "--bogus"]).0, 2);
    assert_eq!(
        explain_raw(&["/nonexistent/model.mps"]).0,
        3,
        "unreadable input"
    );
}

#[test]
fn explain_text_names_the_members() {
    let p = write_tmp("inf3.mps", INF3);
    let (code, out, _) = explain_raw(&[p.to_str().unwrap()]);
    assert_eq!(code, 0);
    for want in [
        "row cap",
        "row needx",
        "row needy",
        "verified: true",
        "VERIFIED",
    ] {
        assert!(out.contains(want), "{want} missing from:\n{out}");
    }
    assert!(!out.contains("row decoy"), "{out}");
}

// ---- JSON matches the committed schema -------------------------------------

/// A small draft-2020-12 subset checker: `$ref` into `$defs`, `oneOf`, `type`
/// (including type lists), `const`, `enum`, `required`, `properties`, `items`,
/// `additionalProperties:false`. Enough to hold the output to the schema's
/// shape rather than to a list written next to the test.
fn conforms(v: &Value, s: &Value, root: &Value, path: &str, errs: &mut Vec<String>) -> bool {
    let mut local = Vec::new();
    if let Some(r) = s.get("$ref").and_then(Value::as_str) {
        let name = r.trim_start_matches("#/$defs/");
        return conforms(v, &root["$defs"][name], root, path, errs);
    }
    if let Some(alts) = s.get("oneOf").and_then(Value::as_array) {
        let ok = alts
            .iter()
            .filter(|a| conforms(v, a, root, path, &mut Vec::new()))
            .count();
        if ok != 1 {
            errs.push(format!("{path}: matches {ok} of the oneOf alternatives"));
            return false;
        }
    }
    if let Some(t) = s.get("type") {
        let types: Vec<&str> = match t {
            Value::String(x) => vec![x.as_str()],
            Value::Array(a) => a.iter().filter_map(Value::as_str).collect(),
            _ => vec![],
        };
        let hit = types.iter().any(|t| match *t {
            "null" => v.is_null(),
            "object" => v.is_object(),
            "array" => v.is_array(),
            "string" => v.is_string(),
            "boolean" => v.is_boolean(),
            "integer" => v.is_i64() || v.is_u64(),
            "number" => v.is_number(),
            _ => false,
        });
        if !hit {
            local.push(format!("{path}: {v} is not of type {t}"));
        }
    }
    if let Some(c) = s.get("const").filter(|c| *c != v) {
        local.push(format!("{path}: {v} != const {c}"));
    }
    if let Some(e) = s
        .get("enum")
        .and_then(Value::as_array)
        .filter(|e| !e.contains(v))
    {
        local.push(format!("{path}: {v} not in enum {e:?}"));
    }
    if let Some(obj) = v.as_object() {
        for f in s
            .get("required")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
        {
            let f = f.as_str().unwrap();
            if !obj.contains_key(f) {
                local.push(format!("{path}.{f} is required and missing"));
            }
        }
        if let Some(props) = s.get("properties").and_then(Value::as_object) {
            for (k, sub) in props {
                if let Some(val) = obj.get(k) {
                    conforms(val, sub, root, &format!("{path}.{k}"), &mut local);
                }
            }
            if s.get("additionalProperties") == Some(&Value::Bool(false)) {
                for k in obj.keys() {
                    if !props.contains_key(k) {
                        local.push(format!("{path}.{k} is not in the schema"));
                    }
                }
            }
        }
    }
    if let (Some(arr), Some(items)) = (v.as_array(), s.get("items")) {
        for (i, x) in arr.iter().enumerate() {
            conforms(x, items, root, &format!("{path}[{i}]"), &mut local);
        }
    }
    let ok = local.is_empty();
    errs.extend(local);
    ok
}

#[test]
fn explain_json_matches_its_schema() {
    let schema: Value = serde_json::from_str(SCHEMA).unwrap();
    assert_eq!(
        schema["properties"]["schema_version"]["const"],
        "rusolve.explain/1"
    );
    for (name, mps) in [
        ("inf3.mps", INF3),
        ("colb.mps", COLB),
        ("unb.mps", UNB),
        ("feas.mps", FEAS),
        ("inf3m.mps", INF3_MIP),
    ] {
        let e = explain(name, mps, &[]);
        let mut errs = Vec::new();
        conforms(&e, &schema, &schema, "$", &mut errs);
        assert!(
            errs.is_empty(),
            "{name} does not match the schema: {errs:?}\n{e}"
        );
        assert_eq!(e["schema_version"], "rusolve.explain/1");
    }
    let skipped = explain("inf3.mps", INF3, &["--no-verify"]);
    let mut errs = Vec::new();
    conforms(&skipped, &schema, &schema, "$", &mut errs);
    assert!(errs.is_empty(), "{errs:?}");
}

// ---- IIS of the MIP itself (--mip-iis) ----------------------------------

/// 2x + 2y = 1 (x, y integer) with two decoy rows: the LP relaxation is feasible,
/// the MIP is not, and the equality row alone is the IIS.
const MIP_PARITY: &str = "\
NAME PAR
ROWS
 N  obj
 E  par
 L  decoy1
 G  decoy2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  par  2
    x  decoy1  1  decoy2  1
    y  obj  1  par  2
    y  decoy1  1
    MARKER  'MARKER'  'INTEND'
    z  obj  1  decoy1  1
    z  decoy2  1
RHS
    rhs  par  1  decoy1  50
    rhs  decoy2  -5
BOUNDS
 UP bnd x 10
 UP bnd y 10
 UP bnd z 10
ENDATA
";

/// r1: x - y >= 0.5 and r2: x - y <= 0.5 (x, y integer) hidden among decoys d1, d2:
/// each alone is feasible, both together force x - y = 0.5.
const MIP_TWO: &str = "\
NAME TWO
ROWS
 N  obj
 G  r1
 L  r2
 L  d1
 G  d2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  r1  1
    x  r2  1  d1  1
    y  obj  1  r1  -1
    y  r2  -1  d2  1
    MARKER  'MARKER'  'INTEND'
    z  obj  1  d1  1
    z  d2  1
RHS
    rhs  r1  0.5  r2  0.5
    rhs  d1  100  d2  -3
BOUNDS
 UP bnd x 10
 UP bnd y 10
 UP bnd z 10
ENDATA
";

/// 15 is not a subset sum of 12, 9, 13, 8, 10 (binaries) while the LP is
/// feasible; k2 and decoy are extra. Some tests need more than the root.
const MIP_HARD: &str = "\
NAME HARD
ROWS
 N  obj
 E  k1
 E  k2
 L  decoy
COLUMNS
    MARKER  'MARKER'  'INTORG'
    a  obj  1  k1  12
    a  k2  7  decoy  1
    b  obj  1  k1  9
    b  k2  11  decoy  1
    c  obj  1  k1  13
    c  k2  5  decoy  1
    d  obj  1  k1  8
    d  k2  14  decoy  1
    e  obj  1  k1  10
    e  k2  9  decoy  1
    f  obj  1  k1  11
    f  k2  6  decoy  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  k1  15  k2  20
    rhs  decoy  500
BOUNDS
 UP bnd a 1
 UP bnd b 1
 UP bnd c 1
 UP bnd d 1
 UP bnd e 1
 UP bnd f 1
ENDATA
";

#[test]
fn explain_mip_integer_only_infeasibility_is_the_equality_row() {
    let e = explain("par.mps", MIP_PARITY, &["--mip-iis"]);
    assert_eq!(e["kind"], "mip", "{e}");
    assert_eq!(e["status"], "Infeasible");
    assert_eq!(e["verified"], true, "{e}");
    let x = &e["explanation"];
    let iis = &x["iis"];
    assert_eq!(iis["kind"], "mip");
    assert_eq!(names(&iis["rows"]), ["par"], "decoys are not in it: {e}");
    assert_eq!(iis["rows"][0]["bound"], "both");
    assert_eq!(names(&iis["columns"]), Vec::<String>::new());
    assert_eq!(iis["irreducible"], true);
    assert_eq!(iis["undetermined"].as_array().unwrap().len(), 0);
    assert!(x["mip_search"]["num_solves"].as_u64().unwrap() > 0);
    // the relaxation is feasible, so the seed is every bound
    assert_eq!(x["mip_search"]["seed"], "all_bounds");
    // verification: set infeasible, both sides needed, witnesses checked
    let v = &x["verification"];
    assert_eq!(v["set_infeasible"], true);
    assert_eq!(v["irreducible"], true);
    let checks = v["checks"].as_array().unwrap();
    assert_eq!(checks.len(), 2, "lower and upper side of par");
    for c in checks {
        assert_eq!(c["rest_feasible"], true);
        assert_eq!(c["witness_checked"], true);
    }
}

#[test]
fn explain_mip_two_member_iis_hidden_among_decoys() {
    let e = explain("two.mps", MIP_TWO, &["--mip-iis"]);
    assert_eq!(e["verified"], true, "{e}");
    let iis = &e["explanation"]["iis"];
    assert_eq!(names(&iis["rows"]), ["r1", "r2"], "{e}");
    for r in iis["rows"].as_array().unwrap() {
        let want = if r["name"] == "r1" { "lower" } else { "upper" };
        assert_eq!(r["bound"], want, "{r}");
    }
    assert_eq!(iis["num_cols"], 0);
    assert_eq!(
        e["explanation"]["verification"]["checks"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    // the default path explains the relaxation, which is feasible: nothing to explain
    let lp = explain("two.mps", MIP_TWO, &[]);
    assert_eq!(lp["kind"], "relaxation");
    assert!(lp["explanation"].is_null());
}

#[test]
fn explain_mip_relaxation_infeasible_model_starts_from_the_relaxation_iis() {
    let e = explain("inf3m.mps", INF3_MIP, &["--mip-iis"]);
    assert_eq!(e["kind"], "mip");
    assert_eq!(e["verified"], true, "{e}");
    assert_eq!(e["explanation"]["mip_search"]["seed"], "relaxation_iis");
    assert_eq!(
        names(&e["explanation"]["iis"]["rows"]),
        ["cap", "needx", "needy"]
    );
}

#[test]
fn explain_mip_limit_hit_reports_undetermined_members() {
    // node limit 0: tests that need branching come back without a verdict
    let e = explain(
        "hard.mps",
        MIP_HARD,
        &["--mip-iis", "--per-test-node-limit", "0"],
    );
    assert_eq!(e["kind"], "mip");
    let x = &e["explanation"];
    assert_eq!(x["type"], "infeasible", "{e}");
    let und = x["mip_search"]["undetermined"].as_array().unwrap();
    assert!(!und.is_empty(), "limit 0 leaves members undetermined: {e}");
    assert_eq!(x["iis"]["irreducible"], false);
    assert_eq!(x["iis"]["undetermined"], x["mip_search"]["undetermined"]);
    assert_eq!(e["verified"], false, "never claim more than was proven");
    // with the default limit the same model gets the exact set
    let full = explain("hard.mps", MIP_HARD, &["--mip-iis"]);
    assert_eq!(full["verified"], true, "{full}");
    assert_eq!(names(&full["explanation"]["iis"]["rows"]), ["k1"]);
    assert_eq!(full["explanation"]["iis"]["irreducible"], true);
}

#[test]
fn explain_mip_time_limit_zero_is_unknown_not_a_claim() {
    let e = explain("par.mps", MIP_PARITY, &["--mip-iis", "--time-limit", "0"]);
    assert_eq!(e["kind"], "mip");
    assert_eq!(e["status"], "Unknown", "{e}");
    assert!(e["explanation"].is_null());
    assert!(e["verified"].is_null());
}

#[test]
fn explain_mip_feasible_model_and_lp_model_unchanged() {
    let feasible = "\
NAME F
ROWS
 N  obj
 G  r1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  r1  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  r1  2
BOUNDS
 UP bnd x 4
ENDATA
";
    let e = explain("f.mps", feasible, &["--mip-iis"]);
    assert_eq!(e["kind"], "mip");
    assert_eq!(e["status"], "Feasible");
    assert!(e["explanation"].is_null());
    // a pure LP ignores --mip-iis: same output as without it
    let a = explain("inf3.mps", INF3, &["--mip-iis"]);
    let b = explain("inf3.mps", INF3, &[]);
    assert_eq!(a, b);
}

#[test]
fn explain_mip_no_verify_and_text_output() {
    let e = explain("par.mps", MIP_PARITY, &["--mip-iis", "--no-verify"]);
    assert_eq!(e["verified"], false);
    assert_eq!(e["explanation"]["verification"]["performed"], false);
    assert_eq!(e["explanation"]["iis"]["irreducible"], true);
    let p = write_tmp("par.mps", MIP_PARITY);
    let (code, out, _) = explain_raw(&[p.to_str().unwrap(), "--mip-iis"]);
    assert_eq!(code, 0);
    for want in [
        "IIS (mip)",
        "row par: both bound",
        "irreducible: yes",
        "MIP solves:",
        "verified: true",
    ] {
        assert!(out.contains(want), "{want} missing from:\n{out}");
    }
    assert!(!out.contains("row decoy"), "{out}");
    let (code, _, _) = explain_raw(&[p.to_str().unwrap(), "--time-limit", "abc"]);
    assert_eq!(code, 2);
}

#[test]
fn explain_mip_json_matches_its_schema() {
    let schema: Value = serde_json::from_str(SCHEMA).unwrap();
    for (name, mps, extra) in [
        ("par.mps", MIP_PARITY, vec!["--mip-iis"]),
        ("two.mps", MIP_TWO, vec!["--mip-iis", "--no-verify"]),
        ("inf3m.mps", INF3_MIP, vec!["--mip-iis"]),
        (
            "hard.mps",
            MIP_HARD,
            vec!["--mip-iis", "--per-test-node-limit", "0"],
        ),
        (
            "par.mps",
            MIP_PARITY,
            vec!["--mip-iis", "--time-limit", "0"],
        ),
    ] {
        let e = explain(name, mps, &extra);
        let mut errs = Vec::new();
        conforms(&e, &schema, &schema, "$", &mut errs);
        assert!(errs.is_empty(), "{name} {extra:?}: {errs:?}\n{e}");
    }
}
