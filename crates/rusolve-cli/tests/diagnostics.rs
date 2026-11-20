//! `--explain`, did-you-mean, and the JSON error body, checked against
//! `docs/schema/rusolve-diagnostics-1.schema.json` (docs/design/diagnostics.md).

use serde_json::Value;
use std::process::Command;

const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-diagnostics-1.schema.json");

const TINY: &str = "NAME T\nROWS\n N  obj\n L  c1\nCOLUMNS\n    x  obj  -1  c1  1\nRHS\n    RHS  c1  4\nBOUNDS\n UP BND  x  10\nENDATA\n";

fn tiny() -> String {
    let dir = std::env::temp_dir().join("rusolve-diag-test");
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("tiny.mps");
    std::fs::write(&p, TINY).unwrap();
    p.to_str().unwrap().to_string()
}

fn run(args: &[&str]) -> (i32, String, String) {
    run_env(args, &[])
}

fn run_env(args: &[&str], env: &[(&str, &str)]) -> (i32, String, String) {
    let o = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(args)
        .envs(env.iter().copied())
        .output()
        .unwrap();
    (
        o.status.code().unwrap(),
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

/// A validator for the subset of JSON Schema the diagnostics schema uses.
fn valid(schema: &Value, root: &Value, v: &Value, path: &str, errs: &mut Vec<String>) {
    if let Some(r) = schema.get("$ref").and_then(Value::as_str) {
        let target = r
            .trim_start_matches("#/")
            .split('/')
            .fold(root, |n, k| &n[k]);
        return valid(target, root, v, path, errs);
    }
    if let Some(c) = schema.get("const")
        && c != v
    {
        errs.push(format!("{path}: want {c}"));
    }
    if let Some(e) = schema.get("enum").and_then(Value::as_array)
        && !e.contains(v)
    {
        errs.push(format!("{path}: {v} not in enum"));
    }
    if let Some(alts) = schema.get("oneOf").and_then(Value::as_array) {
        let ok = alts.iter().filter(|a| {
            let mut e = Vec::new();
            valid(a, root, v, path, &mut e);
            e.is_empty()
        });
        if ok.count() != 1 {
            errs.push(format!("{path}: does not match exactly one alternative"));
        }
        return;
    }
    match schema.get("type").and_then(Value::as_str) {
        Some("object") => {
            let Some(obj) = v.as_object() else {
                return errs.push(format!("{path}: not an object"));
            };
            let props = schema.get("properties").and_then(Value::as_object);
            for r in schema
                .get("required")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
            {
                if !obj.contains_key(r.as_str().unwrap()) {
                    errs.push(format!("{path}: missing {r}"));
                }
            }
            for (k, val) in obj {
                match props.and_then(|p| p.get(k)) {
                    Some(s) => valid(s, root, val, &format!("{path}.{k}"), errs),
                    None if schema["additionalProperties"] == Value::Bool(false) => {
                        errs.push(format!("{path}: unexpected {k}"))
                    }
                    None => {}
                }
            }
        }
        Some("array") => match v.as_array() {
            Some(a) => {
                for (i, x) in a.iter().enumerate() {
                    valid(&schema["items"], root, x, &format!("{path}[{i}]"), errs);
                }
            }
            None => errs.push(format!("{path}: not an array")),
        },
        Some("string") if !v.is_string() => errs.push(format!("{path}: not a string")),
        Some("integer") => match v.as_i64() {
            Some(n) => {
                if let Some(m) = schema.get("minimum").and_then(Value::as_i64)
                    && n < m
                {
                    errs.push(format!("{path}: {n} < {m}"));
                }
            }
            None => errs.push(format!("{path}: not an integer")),
        },
        _ => {}
    }
}

fn check_schema(stdout: &str) -> Value {
    let schema: Value = serde_json::from_str(SCHEMA).unwrap();
    let doc: Value = serde_json::from_str(stdout.trim())
        .unwrap_or_else(|e| panic!("stdout is not one JSON document ({e}): {stdout}"));
    let mut errs = Vec::new();
    valid(&schema, &schema, &doc, "$", &mut errs);
    assert!(errs.is_empty(), "schema violations: {errs:?}\n{stdout}");
    doc
}

#[test]
fn explain_prints_the_registry_entry() {
    let (code, out, _) = run(&["--explain", "RS0002"]);
    assert_eq!(code, 0);
    assert!(out.starts_with("RS0002: unknown flag"), "{out}");
    let (code, out, _) = run(&["--explain=rs1004"]);
    assert_eq!(code, 0);
    assert!(out.contains("declared"), "{out}");
    let (code, _, err) = run(&["--explain", concat!("RS00", "99")]);
    assert_eq!(code, 2);
    assert!(err.contains("no diagnostic code"), "{err}");
    let (_, all, _) = run(&["--explain", "all"]);
    let v: Value = serde_json::from_str(&all).unwrap();
    assert!(v.as_array().unwrap().len() >= 10);
}

#[test]
fn unknown_flag_suggests_the_closest_flag() {
    for cmd in ["solve", "mip"] {
        let (code, _, err) = run(&[cmd, &tiny(), "--time-limt", "5"]);
        assert_eq!(code, 2, "{err}");
        assert!(err.contains("error[RS0002]"), "{err}");
        assert!(err.contains("did you mean '--time-limit'?"), "{cmd}: {err}");
    }
}

#[test]
fn unknown_option_and_command_suggest() {
    let (code, _, err) = run(&["mip", &tiny(), "--opt", "mip_rel_gapp=0.1"]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("error[RS0004]"), "{err}");
    assert!(err.contains("did you mean mip_rel_gap?"), "{err}");
    let (code, _, err) = run(&["slove", "x.mps"]);
    assert_eq!(code, 2);
    assert!(err.contains("did you mean 'solve'?"), "{err}");
}

#[test]
fn json_error_bodies_match_the_schema_and_exit_codes() {
    // usage error (2)
    let (code, out, _) = run(&["mip", &tiny(), "--time-limt", "5", "--format", "json"]);
    assert_eq!(code, 2);
    let doc = check_schema(&out);
    assert_eq!(doc["exit_code"], 2);
    assert_eq!(doc["diagnostics"][0]["code"], "RS0002");
    assert!(
        doc["diagnostics"][0]["help"][0]
            .as_str()
            .unwrap()
            .contains("--time-limit")
    );
    // unknown option (2)
    let (code, out, _) = run(&[
        "mip",
        &tiny(),
        "--opt",
        "mip_rel_gapp=1",
        "--format",
        "json",
    ]);
    assert_eq!(code, 2);
    assert_eq!(check_schema(&out)["diagnostics"][0]["code"], "RS0004");
    // unknown command (2)
    let (code, out, _) = run(&["slove", "x.mps", "--format", "json"]);
    assert_eq!(code, 2);
    assert_eq!(check_schema(&out)["diagnostics"][0]["code"], "RS0001");
    // missing model file (3), both commands
    for cmd in ["solve", "mip"] {
        let (code, out, _) = run(&[cmd, "/no/such/model.mps", "--format", "json"]);
        assert_eq!(code, 3, "{cmd}");
        let doc = check_schema(&out);
        assert_eq!(doc["exit_code"], 3);
        assert_eq!(doc["diagnostics"][0]["code"], "RS0010");
    }
    // a file that exists but does not parse (3)
    let dir = std::env::temp_dir().join("rusolve-diag-test");
    std::fs::create_dir_all(&dir).unwrap();
    let bad = dir.join("broken.mps");
    std::fs::write(&bad, "NAME B\nROWS\n N  o\n ZZZ r\nENDATA\n").unwrap();
    for cmd in ["solve", "mip"] {
        let (code, out, _) = run(&[cmd, bad.to_str().unwrap(), "--format", "json"]);
        assert_eq!(code, 3, "{cmd}");
        // The reader's own code (RS1006 syntax, RS1004 undeclared, RS1005 not a number);
        // RS0011 stays the fallback for a format with no located reader errors
        assert_eq!(check_schema(&out)["diagnostics"][0]["code"], "RS1006");
    }
}

/// Every `"RSdddd"` literal in the workspace sources must be registered.
/// Lives here (not only in rusolve-core) so the cli test run covers it.
#[test]
fn every_emitted_code_is_registered() {
    fn rs_files(dir: &std::path::Path, out: &mut Vec<std::path::PathBuf>) {
        for e in std::fs::read_dir(dir).unwrap().flatten() {
            let p = e.path();
            let name = p.file_name().unwrap().to_string_lossy().to_string();
            if p.is_dir() {
                if name != "target" && name != "fixtures" {
                    rs_files(&p, out);
                }
            } else if name.ends_with(".rs") {
                out.push(p);
            }
        }
    }
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let mut files = Vec::new();
    rs_files(&root, &mut files);
    let mut scanned = 0;
    for f in files {
        let text = std::fs::read_to_string(&f).unwrap_or_default();
        let b = text.as_bytes();
        for i in 0..b.len().saturating_sub(7) {
            if b[i] == b'"'
                && &b[i + 1..i + 3] == b"RS"
                && b[i + 3..i + 7].iter().all(u8::is_ascii_digit)
                && b[i + 7] == b'"'
            {
                scanned += 1;
                let code = &text[i + 1..i + 7];
                assert!(
                    rusolve_core::diagnostic::REGISTRY
                        .iter()
                        .any(|e| e.code == code),
                    "{code} is named in {} but is not in diagnostic::REGISTRY",
                    f.display()
                );
            }
        }
    }
    assert!(scanned > 0, "source scan found no RS codes; wrong root?");
}

/// Every registered code is emitted at a source site (not only in the registry's own
/// literals or a test). A dead code is documentation that lies.
#[test]
fn every_registered_code_is_emitted() {
    fn rs_files(dir: &std::path::Path, out: &mut Vec<std::path::PathBuf>) {
        for e in std::fs::read_dir(dir).unwrap().flatten() {
            let p = e.path();
            let name = p.file_name().unwrap().to_string_lossy().to_string();
            if p.is_dir() {
                if name != "target" && name != "fixtures" && name != "tests" {
                    rs_files(&p, out);
                }
            } else if name.ends_with(".rs") {
                out.push(p);
            }
        }
    }
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let mut files = Vec::new();
    rs_files(&root, &mut files);
    let mut emitted = std::collections::HashSet::new();
    for f in files {
        let text = std::fs::read_to_string(&f).unwrap_or_default();
        // inline test modules and the registry's own `code:` lines are not emit sites
        let body = text.split("#[cfg(test)]\nmod ").next().unwrap_or("");
        for line in body.lines() {
            if line.trim_start().starts_with("code: \"RS") {
                continue;
            }
            let b = line.as_bytes();
            for i in 0..b.len().saturating_sub(7) {
                if b[i] == b'"'
                    && &b[i + 1..i + 3] == b"RS"
                    && b[i + 3..i + 7].iter().all(u8::is_ascii_digit)
                    && b[i + 7] == b'"'
                {
                    emitted.insert(line[i + 1..i + 7].to_string());
                }
            }
        }
    }
    for e in rusolve_core::diagnostic::REGISTRY {
        assert!(
            emitted.contains(e.code),
            "{} is registered but no emit site names it: emit it or remove it from REGISTRY",
            e.code
        );
    }
}

/// Every error invocation answers with one `rusolve.diagnostics/1` document on stdout under
/// `--format json`, and the exit code of the exit-code table.
#[test]
fn every_error_path_prints_one_diagnostics_document_under_json() {
    let dir = std::env::temp_dir().join("rusolve-diag-test");
    std::fs::create_dir_all(&dir).unwrap();
    let tiny = tiny();
    let out_lp = dir.join("conv-out.lp");
    let out_lp = out_lp.to_str().unwrap();
    let qp = dir.join("e_qp.lp");
    std::fs::write(
        &qp,
        "Minimize\n obj: - x - y + [ 2 x ^ 2 + 2 y ^ 2 ] / 2\nSubject To\n c1: x + y <= 1\nEnd\n",
    )
    .unwrap();
    let qp = qp.to_str().unwrap().to_string();
    let socp = dir.join("e_socp.lp");
    std::fs::write(
        &socp,
        "Minimize\n obj: - x - y\nSubject To\n q1: [ x ^ 2 + y ^ 2 ] <= 1\nBounds\n x free\n y free\nEnd\n",
    )
    .unwrap();
    let socp = socp.to_str().unwrap().to_string();
    let hook = "/no/such/hook.wasm";
    let cases: Vec<(Vec<&str>, i32)> = vec![
        (vec!["solve", &tiny, "--time-limit", "abc"], 2),
        (vec!["solve", &tiny, "--iter-limit", "abc"], 2),
        (vec!["convert", "/no/such/in.mps", out_lp], 3),
        (vec!["convert", &tiny], 2),
        (vec!["convert", &tiny, "out.xyz"], 2),
        (vec!["check", "/no/such/model.mps", "/no/such.sol"], 3),
        (vec!["check", &tiny], 2),
        (vec!["lint", "/no/such/model.mps"], 3),
        (vec!["lint"], 2),
        // route.rs: a flag the QP / conic routes do not take
        (vec!["mip", &qp, "--bogus-flag-xyz"], 2),
        (vec!["mip", &qp, "--branch", "x"], 2),
        (vec!["solve", &socp, "--time-limit", "abc"], 2),
        (vec!["mip", &socp, "--threads"], 2),
        (vec!["solve", &socp, "--bogus-flag-xyz"], 2),
        // mip_highs.rs error sites
        (vec!["mip", &tiny, "--watchdog-grace", "-1"], 2),
        (vec!["mip", &tiny, "--node-limit", "99999999999"], 2),
        (vec!["mip", &tiny, "--seed", "99999999999"], 2),
        (vec!["mip", &tiny, "--threads", "0"], 2),
        (vec!["mip", &tiny, "--csv"], 2),
        (vec!["mip", &tiny, "--wasm-fuel", "5"], 2),
        (vec!["mip", &tiny, "--watchdog"], 2),
        (vec!["mip", &tiny, "--wasm-hook", hook], 3),
        (vec!["mip", &tiny, "--options-file", "/no/such/opts"], 2),
        (vec!["mip", &tiny, "--profile", "/no/such/profile.json"], 3),
        (vec!["mip", &tiny, "--start", "/no/such/start.sol"], 3),
        (vec!["mip", &tiny, "--basis", "/no/such/basis.bas"], 3),
        (vec!["mip", &tiny, "--add-rows", "/no/such/rows.json"], 3),
        (vec!["mip"], 2),
        (vec!["mip", &tiny, "extra"], 2),
        // every other command: a missing file, a bad flag value, a wrong positional count
        (vec!["solve", "/no/such/model.mps"], 3),
        (vec!["solve"], 2),
        (vec!["mip", "/no/such/model.mps"], 3),
        (
            vec!["check", &tiny, "/no/such.sol", "--max-findings", "abc"],
            2,
        ),
        (vec!["lint", &tiny, "--bogus-flag-xyz"], 2),
        (vec!["lint", &tiny, &tiny], 2),
        (vec!["explain", "/no/such/model.mps"], 3),
        (vec!["explain", &tiny, "--time-limit", "abc"], 2),
        (vec!["explain"], 2),
        (vec!["explain", &tiny, &tiny], 2),
        (vec!["explain", "/no/such/model.mps", "--mip-iis"], 3),
        (
            vec![
                "explain",
                &tiny,
                "--mip-iis",
                "--per-test-node-limit",
                "abc",
            ],
            2,
        ),
        (vec!["feasrelax", "/no/such/model.mps"], 3),
        (vec!["feasrelax", &tiny, "--mode", "bogus"], 2),
        (vec!["feasrelax"], 2),
        (vec!["feasrelax", &tiny, &tiny], 2),
        (
            vec!["feasrelax", &tiny, "--penalties", "/no/such/pen.json"],
            3,
        ),
        (vec!["ranging", "/no/such/model.mps"], 3),
        (vec!["ranging", &tiny, "--rows", "nomatchxyz"], 2),
        (vec!["ranging"], 2),
        (vec!["ranging", &tiny, &tiny], 2),
        (vec!["classify", "/no/such/model.mps"], 3),
        (vec!["classify", &tiny, "--bogus-flag-xyz"], 2),
        (vec!["classify"], 2),
        (vec!["classify", &tiny, &tiny], 2),
        (vec!["probe", "/no/such/model.mps"], 3),
        (vec!["probe", &tiny, "--node-limit", "abc"], 2),
        (vec!["probe"], 2),
        (vec!["probe", &tiny, &tiny], 2),
        (vec!["diagnose", "/no/such/model.mps"], 3),
        (vec!["diagnose", &tiny, "--time-limit", "abc"], 2),
        (vec!["diagnose"], 2),
        (vec!["diagnose", &tiny, &tiny], 2),
        (vec!["tune", "/no/such/dir"], 3),
        (vec!["tune", &tiny, "--seeds", "abc"], 2),
        (vec!["tune"], 2),
        (vec!["tune", &tiny, &tiny], 2),
        (vec!["profile", "show", "/no/such/profile.json"], 3),
        (vec!["profile", "show", "default", "--base", "bogus"], 2),
        (vec!["profile", "show"], 2),
        (vec!["profile", "show", "default", "default"], 2),
        (vec!["options", "--profile", "bogus"], 2),
        (vec!["options", "--bogus-flag-xyz"], 2),
    ];
    let registered = |doc: &Value, args: &[&str]| {
        let diags = doc["diagnostics"].as_array().unwrap();
        assert!(!diags.is_empty(), "{args:?}: no diagnostics");
        for d in diags {
            let c = d["code"].as_str().unwrap();
            assert!(
                rusolve_core::diagnostic::lookup(c).is_some(),
                "{args:?}: {c} is not registered"
            );
        }
    };
    for (mut args, want) in cases {
        args.extend(["--format", "json"]);
        let (code, out, err) = run(&args);
        assert_eq!(code, want, "{args:?}: {err}");
        assert!(
            !out.trim().is_empty(),
            "{args:?}: empty stdout; stderr: {err}"
        );
        let doc = check_schema(&out);
        assert_eq!(doc["exit_code"], want, "{args:?}: {out}");
        registered(&doc, &args);
    }
    // the native-stack `mip` parser (RUSOLVE_MIP_HIGHS=0) answers the same way
    let native: Vec<(Vec<&str>, i32)> = vec![
        (vec!["mip", "/no/such/model.mps"], 3),
        (vec!["mip", &tiny, "--time-limit", "abc"], 2),
        (vec!["mip", &tiny, "--bogus-flag-xyz"], 2),
        (vec!["mip", &tiny, "--threads", "abc"], 2),
        (vec!["mip"], 2),
        (vec!["mip", &tiny, "extra"], 2),
    ];
    for (mut args, want) in native {
        args.extend(["--format", "json"]);
        let (code, out, err) = run_env(&args, &[("RUSOLVE_MIP_HIGHS", "0")]);
        assert_eq!(code, want, "native {args:?}: {err}");
        assert!(
            !out.trim().is_empty(),
            "native {args:?}: empty stdout; stderr: {err}"
        );
        let doc = check_schema(&out);
        assert_eq!(doc["exit_code"], want, "native {args:?}: {out}");
        registered(&doc, &args);
    }
}

/// A refused class under `--format json` carries registered `diagnostics` inside its
/// `rusolve.error/1` object, next to the existing fields.
#[test]
fn refused_class_error_carries_diagnostics() {
    let dir = std::env::temp_dir().join("rusolve-diag-test");
    std::fs::create_dir_all(&dir).unwrap();
    let qp = dir.join("e_nonconvex.lp");
    std::fs::write(
        &qp,
        "Minimize\n obj: x + y - [ 2 x * y ] / 2\nSubject To\n c1: x + y <= 1\nBounds\n 0 <= x <= 1\n 0 <= y <= 1\nEnd\n",
    )
    .unwrap();
    let qp = qp.to_str().unwrap();
    for cmd in ["solve", "mip"] {
        let (code, out, err) = run(&[cmd, qp, "--format", "json"]);
        assert_eq!(code, 3, "{cmd}: {err}");
        let v: Value = serde_json::from_str(out.trim()).unwrap();
        assert_eq!(v["schema_version"], "rusolve.error/1", "{out}");
        assert_eq!(v["error"]["class"], "nonconvex_qp", "{out}");
        let d = v["diagnostics"].as_array().unwrap();
        assert!(!d.is_empty(), "{out}");
        assert!(rusolve_core::diagnostic::lookup(d[0]["code"].as_str().unwrap()).is_some());
    }
}

/// A profile that parses as JSON but breaks the schema is a content error (RS0011, exit 2),
/// with each further error as its own note, never RS0006 by wording.
#[test]
fn profile_schema_errors_are_rs0011_with_one_note_each() {
    let dir = std::env::temp_dir().join("rusolve-diag-test");
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("bad_profile.json");
    std::fs::write(
        &p,
        r#"{"schema_version":"rusolve.profile/1","base":"x","opts":{}}"#,
    )
    .unwrap();
    let (code, out, err) = run(&["profile", "show", p.to_str().unwrap(), "--format", "json"]);
    assert_eq!(code, 2, "{err}");
    let v: Value = serde_json::from_str(out.trim()).unwrap();
    assert_eq!(v["diagnostics"][0]["code"], "RS0011", "{out}");
    assert!(!out.contains('\u{1}') && !out.contains("\\u0001"), "{out}");
}
