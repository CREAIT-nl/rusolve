//! The environment-variable registry, run against the real binary: a misspelled
//! `RUSOLVE_*` name is a usage error with a suggestion, a registered name is
//! not, and `rusolve env` is how an agent discovers them.
//!
//! The registry itself (and the tests that keep it complete) lives in
//! `src/env_registry.rs`; these tests are the process-level contract.

use std::path::PathBuf;
use std::process::Command;

/// A knapsack MPS, 20 columns and 3 rows: the LP path solves its relaxation
/// to Optimal in milliseconds, which is all these tests ask of a solve.
fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../rusolve-ekk/tests/fixtures/s5/e2e_knap.mps")
}

/// Runs the CLI with a clean-ish environment plus `extra`, so a variable the
/// developer happens to have exported cannot change the answer.
fn cli(args: &[&str], extra: &[(&str, &str)]) -> (i32, String, String) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    for key in [
        "RUSOLVE_PROFLE",
        "RUSOLVE_PROFILE",
        "RUSOLVE_DUAL_STATS",
        "RUSOLVE_MIP_HIGHS",
    ] {
        cmd.env_remove(key);
    }
    for (k, v) in extra {
        cmd.env(k, v);
    }
    let out = cmd.args(args).output().expect("rusolve-cli runs");
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// A typo is a usage error (2), not a silent engine switch: the message names
/// the variable, suggests the nearest registered name and points at
/// `rusolve env`.
#[test]
fn env_registry_misspelled_variable_exits_2_with_a_suggestion() {
    let (code, _out, err) = cli(
        &["solve", fixture().to_str().unwrap()],
        &[("RUSOLVE_PROFLE", "highs-1.15.1")],
    );
    assert_eq!(code, 2, "stderr: {err}");
    assert!(
        err.contains(
            "rusolve: unknown environment variable RUSOLVE_PROFLE (rusolve env lists them all)"
        ) && err.contains("= help: did you mean RUSOLVE_PROFILE?"),
        "unexpected stderr: {err}"
    );
    // It never got as far as solving, so nothing was silently switched.
    assert!(!err.contains("legacy-engine lever"), "stderr: {err}");
}

/// A name no registry entry resembles is still an error, without a suggestion.
#[test]
fn env_registry_unrelated_variable_exits_2_without_a_suggestion() {
    let (code, _out, err) = cli(
        &["solve", fixture().to_str().unwrap()],
        &[("RUSOLVE_NOT_A_REAL_LEVER_AT_ALL", "1")],
    );
    assert_eq!(code, 2, "stderr: {err}");
    assert!(
        err.contains("unknown environment variable RUSOLVE_NOT_A_REAL_LEVER_AT_ALL"),
        "stderr: {err}"
    );
    assert!(!err.contains("did you mean"), "stderr: {err}");
}

/// A registered variable that is not a lever runs the solve normally: the
/// profile is honoured and the engine stays on the default (ekk) path.
#[test]
fn env_registry_registered_non_lever_runs_normally() {
    let (code, out, err) = cli(
        &["solve", fixture().to_str().unwrap()],
        &[("RUSOLVE_PROFILE", "recommended")],
    );
    assert_eq!(code, 0, "stderr: {err}\nstdout: {out}");
    assert!(out.contains("status:  Optimal"), "stdout: {out}");
    assert!(!err.contains("legacy-engine lever"), "stderr: {err}");
}

/// The same solve with no environment at all, as the reference for the test
/// above: the profile must not move the answer.
#[test]
fn env_registry_profile_and_default_agree() {
    let path = fixture();
    let default = cli(&["solve", path.to_str().unwrap()], &[]);
    let profile = cli(
        &["solve", path.to_str().unwrap()],
        &[("RUSOLVE_PROFILE", "recommended")],
    );
    assert_eq!(default.0, profile.0);
    assert_eq!(default.1, profile.1, "the profile changed the solve output");
}

/// `rusolve env --json` parses, and lists every name the plain table lists.
#[test]
fn env_registry_env_json_lists_every_name() {
    let (code, out, err) = cli(&["env", "--json"], &[]);
    assert_eq!(code, 0, "stderr: {err}");
    let json: serde_json::Value = serde_json::from_str(&out).expect("`env --json` is valid JSON");
    let arr = json.as_array().expect("`env --json` is a JSON array");
    assert!(arr.len() > 100, "only {} entries listed", arr.len());

    let mut names = Vec::new();
    for e in arr {
        for field in ["name", "read_by", "kind", "scope", "values", "doc"] {
            assert!(e.get(field).is_some(), "entry {e} has no `{field}`");
        }
        assert!(
            matches!(
                e["kind"].as_str(),
                Some("lever") | Some("profile") | Some("debug") | Some("harness")
            ),
            "entry {e} has an unknown kind"
        );
        assert!(
            matches!(
                e["scope"].as_str(),
                Some("legacy_lp")
                    | Some("ekk")
                    | Some("port_mip")
                    | Some("native_mip")
                    | Some("any")
            ),
            "entry {e} has an unknown scope"
        );
        assert!(
            e["read_by"]
                .as_array()
                .is_some_and(|r| !r.is_empty() && r.iter().all(|f| f.is_string())),
            "entry {e} has no reader file list"
        );
        names.push(e["name"].as_str().expect("name is a string").to_string());
    }
    // Names an agent must be able to find, across every kind and scope.
    for want in [
        "RUSOLVE_PROFILE",
        "RUSOLVE_MIP_HIGHS",
        "RUSOLVE_EKK",
        "RUSOLVE_DUAL_STATS",
        "RUSOLVE_PDLP_TRACE",
        "RUSOLVE_MIP_LOG",
        "RUSOLVE_HIGHS_FIXTURES",
        "PD_DETAIL",
    ] {
        assert!(names.iter().any(|n| n == want), "{want} is not listed");
    }

    let (_tcode, text, _terr) = cli(&["env"], &[]);
    let text_names: Vec<&str> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().next().unwrap_or_default())
        .collect();
    assert_eq!(
        text_names.len(),
        names.len(),
        "the plain table and --json list a different number of variables"
    );
    for n in &names {
        assert!(
            text_names.contains(&n.as_str()),
            "{n} missing from the plain table"
        );
    }
}
