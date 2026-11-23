//! (`--log-format jsonl`) and (`--watchdog`) on the port MIP path
//! (`RUSOLVE_MIP_HIGHS=1 rusolve mip`): the NDJSON event stream, its schema, and the
//! wall-clock watchdog, run through the real binary.
//!
//! Note: `progress_watchdog_kills_the_fiball_cover_loop` is a timing test calibrated for an
//! optimised binary and is ignored in debug builds. Run it with
//! `cargo test --release -p rusolve-cli --test progress`.

use std::path::PathBuf;
use std::process::Command;
use std::time::Instant;

/// `docs/schema/rusolve-event-1.schema.json`, the normative event schema.
const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-event-1.schema.json");

const TINY_MIP: &str = "\
NAME PRG1
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

fn tiny_model() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-progress-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("tiny.mps");
    std::fs::write(&p, TINY_MIP).unwrap();
    p
}

/// fiball from the benchmark corpus (`tools/bench/fetch_corpus.py`), `None` when absent.
fn fiball() -> Option<PathBuf> {
    // the shell expands `~`; reading the home variable here would need an env registry entry
    let out = Command::new("sh")
        .arg("-c")
        .arg("echo ~/.cache/rusolve-bench/miplib-benchmark/fiball.mps.gz")
        .output()
        .ok()?;
    let p = PathBuf::from(String::from_utf8(out.stdout).ok()?.trim());
    p.exists().then_some(p)
}

/// `(exit, stdout, stderr)` of `RUSOLVE_MIP_HIGHS=1 rusolve mip <model> <extra>`.
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

/// The `required` list the schema demands for `event` (top level plus its `then` branch).
fn required_for(event: &str) -> Vec<String> {
    let schema: serde_json::Value = serde_json::from_str(SCHEMA).unwrap();
    let list = |v: &serde_json::Value| -> Vec<String> {
        v.as_array()
            .map(|a| a.iter().map(|s| s.as_str().unwrap().to_string()).collect())
            .unwrap_or_default()
    };
    let mut out = list(&schema["required"]);
    for branch in schema["allOf"].as_array().unwrap() {
        if branch["if"]["properties"]["event"]["const"] == event {
            out.extend(list(&branch["then"]["required"]));
        }
    }
    out
}

/// Parses every line of `stderr` as an event and checks it against the schema's `required`
/// lists, `seq` and `schema_version`.
fn events(stderr: &str) -> Vec<serde_json::Value> {
    let events: Vec<serde_json::Value> = stderr
        .lines()
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("not JSON: {e}: {l}")))
        .collect();
    assert!(!events.is_empty(), "no events on stderr");
    for (i, ev) in events.iter().enumerate() {
        assert_eq!(ev["schema_version"], "rusolve.event/1", "{ev}");
        assert_eq!(ev["seq"], (i + 1) as u64, "seq counts lines: {ev}");
        let name = ev["event"].as_str().unwrap();
        for field in required_for(name) {
            assert!(ev.get(&field).is_some(), "{name}.{field} missing: {ev}");
        }
    }
    events
}

/// Every line parses, the stream opens with the presolve phase, and the last event is
/// the `final` event whose status, objective, bound and gap are the envelope's.
#[test]
fn progress_jsonl_last_event_matches_the_envelope() {
    let model = tiny_model();
    let (code, stdout, stderr) = mip(&model, &["--format", "json", "--log-format", "jsonl"]);
    assert_eq!(code, 0, "stderr: {stderr}");
    let env: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    let evs = events(&stderr);
    assert_eq!(evs[0]["event"], "phase");
    assert_eq!(evs[0]["phase"], "presolve");
    let last = evs.last().unwrap();
    assert_eq!(last["event"], "final");
    for field in ["status", "objective", "dual_bound", "gap", "has_solution"] {
        assert_eq!(last[field], env[field], "final.{field} vs the envelope");
    }
    assert_eq!(last["nodes"], env["stats"]["nodes"]);
    assert_eq!(last["status"], "Optimal");
    assert_eq!(last["objective"], -20.0);
    assert!(
        evs.iter().any(|e| e["event"] == "incumbent"),
        "an incumbent event is expected: {stderr}"
    );
    for ev in evs.iter().filter(|e| e["event"] == "incumbent") {
        assert!(
            ev["source"].is_string(),
            "an incumbent names its source: {ev}"
        );
    }
}

/// The event stream is orthogonal to the stdout format: the CSV line is unchanged by it.
#[test]
fn progress_jsonl_keeps_stdout_as_it_was() {
    let model = tiny_model();
    let (_, plain, _) = mip(&model, &["--csv"]);
    let (code, with, stderr) = mip(&model, &["--csv", "--log-format", "jsonl"]);
    assert_eq!(code, 0);
    // time_s differs, the first five fields do not
    let key = |s: &str| {
        s.lines()
            .last()
            .unwrap()
            .split(',')
            .take(5)
            .collect::<Vec<_>>()
            .join(",")
    };
    assert_eq!(key(&plain), key(&with));
    events(&stderr);
}

/// `--watchdog` needs a time limit; a bad `--log-format` is a usage error.
#[test]
fn progress_usage_errors_exit_2() {
    let model = tiny_model();
    assert_eq!(mip(&model, &["--watchdog"]).0, 2);
    assert_eq!(mip(&model, &["--log-format", "xml"]).0, 2);
    assert_eq!(
        mip(
            &model,
            &["--watchdog", "--time-limit", "5", "--watchdog-grace", "x"]
        )
        .0,
        2
    );
}

/// A solve that finishes inside its limit is untouched by the watchdog: no `limit`.
#[test]
fn progress_watchdog_is_silent_when_the_solve_finishes() {
    let model = tiny_model();
    let (code, stdout, stderr) = mip(
        &model,
        &["--format", "json", "--watchdog", "--time-limit", "30"],
    );
    assert_eq!(code, 0, "stderr: {stderr}");
    let env: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    assert_eq!(env["status"], "Optimal");
    assert!(env.get("limit").is_none(), "no limit hit: {env}");
}

/// Kill test: fiball, seed 1, hybrid cut selection on, the cover guard
/// off - the 1.16.0 cover loop never reaches a limit check, so only the hard deadline can
/// end the solve. It must return within limit plus grace, with the incumbent, bound, gap
/// and the improvement information in the envelope and a matching `final` event.
#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "timing test: a debug build reaches no limit check within the watchdog window; run with --release"
)]
fn progress_watchdog_kills_the_fiball_cover_loop() {
    let Some(model) = fiball() else {
        eprintln!("fiball.mps.gz is not in the bench cache: kill test skipped");
        return;
    };
    // The hang sits in the RENS sub-MIP, so the run has to get that far before the limit.
    // On a starved machine it may still be in code that checks the limit when it fires:
    // that is the legitimate `watchdog_interrupt` outcome (see the next test), not a kill.
    // Retry with a doubled limit until the hang is reached.
    let grace = 2.0;
    let (mut limit, mut attempt) = (4.0, 0);
    let (stderr, env) = loop {
        let start = Instant::now();
        let (code, stdout, stderr) = mip(
            &model,
            &[
                "--format",
                "json",
                "--log-format",
                "jsonl",
                "--watchdog",
                "--watchdog-grace",
                "2",
                "--time-limit",
                &limit.to_string(),
                "--seed",
                "1",
                "--opt",
                "mip_cut_selection_hybrid=true",
                "--opt",
                "mip_strengthen_cover_guard=false",
            ],
        );
        let wall = start.elapsed().as_secs_f64();
        assert_eq!(code, 0, "stderr: {stderr}");
        assert!(
            wall < limit + grace + 5.0,
            "returned after {wall:.1}s, limit {limit} + grace {grace}"
        );
        let env: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
        assert_eq!(env["status"], "TimeLimit");
        assert_eq!(env["limit"]["kind"], "time");
        let by = env["limit"]["enforced_by"].as_str().unwrap().to_string();
        if by == "watchdog_kill" {
            break (stderr, env);
        }
        assert_eq!(by, "watchdog_interrupt", "{env}");
        attempt += 1;
        assert!(
            attempt < 4,
            "the hang was not reached within {limit}s: {env}"
        );
        limit *= 2.0;
    };
    assert_eq!(env["limit"]["enforced_by"], "watchdog_kill");
    assert_eq!(
        env["has_solution"], true,
        "an incumbent was found before the hang"
    );
    assert!(env["objective"].is_number() && env["dual_bound"].is_number());
    assert!(env["gap"].is_number());
    let improving = &env["limit"]["improving"];
    assert_eq!(improving["tracked"], true);
    assert!(improving["count"].as_u64().unwrap() >= 1);
    assert_eq!(
        improving["still_improving"], false,
        "stuck, not improving: {improving}"
    );
    let evs = events(&stderr);
    let last = evs.last().unwrap();
    assert_eq!(last["event"], "final");
    for field in ["status", "objective", "dual_bound", "gap"] {
        assert_eq!(last[field], env[field], "final.{field} vs the envelope");
    }
    let actions: Vec<_> = evs
        .iter()
        .filter(|e| e["event"] == "watchdog")
        .map(|e| e["action"].as_str().unwrap().to_string())
        .collect();
    assert_eq!(actions, ["interrupt_requested", "hard_stop"]);
}

/// With the cover guard on (the default) fiball reaches limit checks, so the interrupt
/// path ends the solve: `watchdog_interrupt`, status TimeLimit, within limit plus grace.
#[test]
fn progress_watchdog_interrupts_a_solve_that_reaches_limit_checks() {
    let Some(model) = fiball() else {
        eprintln!("fiball.mps.gz is not in the bench cache: test skipped");
        return;
    };
    let start = Instant::now();
    let (code, stdout, stderr) = mip(
        &model,
        &[
            "--format",
            "json",
            "--watchdog",
            "--time-limit",
            "3",
            "--seed",
            "1",
        ],
    );
    let wall = start.elapsed().as_secs_f64();
    assert_eq!(code, 0, "stderr: {stderr}");
    assert!(wall < 3.0 + 5.0 + 5.0, "returned after {wall:.1}s");
    let env: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    assert_eq!(env["status"], "TimeLimit", "{env}");
    assert!(
        ["solver", "watchdog_interrupt"].contains(&env["limit"]["enforced_by"].as_str().unwrap()),
        "{env}"
    );
    assert_eq!(env["limit"]["improving"]["tracked"], true);
}
