//! `rusolve profile` and `rusolve mip --profile`. Every test runs the real binary.

use std::path::{Path, PathBuf};
use std::process::Command;

const PROFILE_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-profile-1.schema.json");
const DIFF_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-profile-diff-1.schema.json");

const TINY_MIP: &str = "\
NAME PROF1
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

fn tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-profile-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn cli(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(args)
        .env_remove("RUSOLVE_PROFILE")
        .env_remove("RUSOLVE_MIP_HIGHS")
        .output()
        .unwrap();
    (
        out.status.code().unwrap(),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8(out.stderr).unwrap(),
    )
}

fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}

fn assert_required(schema: &str, v: &serde_json::Value) {
    let schema: serde_json::Value = serde_json::from_str(schema).unwrap();
    for k in schema["required"].as_array().unwrap() {
        let k = k.as_str().unwrap();
        assert!(v.get(k).is_some(), "{k} missing: {v}");
    }
}

#[test]
fn profile_write_then_read_round_trips() {
    let file = tmp("rt.json", "");
    let (code, _, err) = cli(&[
        "profile",
        "new",
        "--name",
        "rt",
        "--opt",
        "mip_heuristic_effort=0.15",
        "--opt",
        "mip_pscost_minreliable=16",
        "--opt",
        "mip_detect_symmetry=false",
        "-o",
        p(&file),
    ]);
    assert_eq!(code, 0, "{err}");
    let written: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&file).unwrap()).unwrap();
    assert_eq!(written["schema_version"], "rusolve.profile/1");
    assert_required(PROFILE_SCHEMA, &written);
    assert_eq!(written["options"]["mip_pscost_minreliable"], 16);
    // read it back: validate, show (JSON round trip) and diff
    let (code, out, err) = cli(&["profile", "validate", p(&file)]);
    assert_eq!(code, 0, "{err}");
    assert!(out.contains("3 option(s)"), "{out}");
    let (code, shown, _) = cli(&["profile", "show", p(&file), "--format", "json"]);
    assert_eq!(code, 0);
    let shown: serde_json::Value = serde_json::from_str(shown.trim()).unwrap();
    assert_eq!(shown, written);
    let (code, out, _) = cli(&["profile", "diff", p(&file), "recommended"]);
    assert_eq!(code, 0);
    assert!(out.contains("3 option(s) differ"), "{out}");
    assert!(out.contains("mip_heuristic_effort : 0.15 -> 0.05"), "{out}");
}

#[test]
fn profile_diff_json_against_default_recommended_and_files() {
    let (code, out, _) = cli(&[
        "profile",
        "diff",
        "default",
        "recommended",
        "--format",
        "json",
    ]);
    assert_eq!(code, 0);
    let v: serde_json::Value = serde_json::from_str(out.trim()).unwrap();
    assert_eq!(v["schema_version"], "rusolve.profile-diff/1");
    assert_required(DIFF_SCHEMA, &v);
    assert_eq!(v["identical"], false);
    assert!(!v["differences"].as_array().unwrap().is_empty());
    // a profile on base default with no options equals `default`
    let f = tmp(
        "d.json",
        r#"{"schema_version":"rusolve.profile/1","base":"default"}"#,
    );
    let (code, out, _) = cli(&["profile", "diff", p(&f), "default", "--format", "json"]);
    assert_eq!(code, 0);
    let v: serde_json::Value = serde_json::from_str(out.trim()).unwrap();
    assert_eq!(v["identical"], true);
}

#[test]
fn profile_invalid_file_fails_with_a7_errors() {
    let f = tmp(
        "bad.json",
        r#"{"schema_version":"rusolve.profile/1","options":{
            "mip_pscost_minrelaible": 4, "mip_heuristic_effort": 9.5, "time_limit": 3}}"#,
    );
    let (code, _, err) = cli(&["profile", "validate", p(&f)]);
    assert_eq!(code, 2, "{err}");
    assert!(
        err.contains("did you mean mip_pscost_minreliable?"),
        "{err}"
    );
    assert!(err.contains("outside the range"), "{err}");
    assert!(err.contains("may not set this option"), "{err}");
    // the same file refuses to start a solve
    let model = tmp("tiny.mps", TINY_MIP);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["mip", p(&model), "--profile", p(&f)])
        .env("RUSOLVE_MIP_HIGHS", "1")
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("did you mean"));
    // wrong schema version, not JSON, missing file
    let f2 = tmp("v2.json", r#"{"schema_version":"rusolve.profile/2"}"#);
    assert_eq!(cli(&["profile", "validate", p(&f2)]).0, 2);
    let f3 = tmp("nojson.json", "not json");
    assert_eq!(cli(&["profile", "validate", p(&f3)]).0, 3);
    assert_eq!(cli(&["profile", "validate", "/nonexistent/x.json"]).0, 3);
}

#[test]
fn profile_applies_in_mip_and_the_command_line_still_wins() {
    let model = tmp("tiny2.mps", TINY_MIP);
    let f = tmp(
        "ok.json",
        r#"{"schema_version":"rusolve.profile/1","base":"default",
            "options":{"mip_pscost_minreliable": 3}}"#,
    );
    let envs = |extra: &[&str]| {
        let mut args = vec!["mip", p(&model), "--format", "json", "--profile", p(&f)];
        args.extend_from_slice(extra);
        let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
            .args(&args)
            .env("RUSOLVE_MIP_HIGHS", "1")
            .env_remove("RUSOLVE_PROFILE")
            .output()
            .unwrap();
        assert_eq!(
            out.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        serde_json::from_slice::<serde_json::Value>(&out.stdout).unwrap()
    };
    let v = envs(&[]);
    assert_eq!(v["status"], "Optimal");
    // --opt after the profile is applied last: an invalid value there is still exit 2
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args([
            "mip",
            p(&model),
            "--profile",
            p(&f),
            "--opt",
            "mip_pscost_minreliable=oops",
        ])
        .env("RUSOLVE_MIP_HIGHS", "1")
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(2));
}
