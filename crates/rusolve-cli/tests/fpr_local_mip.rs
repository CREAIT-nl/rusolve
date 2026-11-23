//! Integration tests for the mip-heuristics FPR / LocalMIP chain
//!
//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! This file is not edited for the other tests (highs_backport_options.rs, mip_trace.rs, mip_trace_threads.rs).

use std::path::{Path, PathBuf};
use std::process::Command;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

// ---------------------------------------------------------------------------
// Helpers copied from conflict_resolution.rs
// ---------------------------------------------------------------------------

/// The main checkout's root (from a git worktree) or this checkout's root, whichever is
/// first found to hold `rel` (mip_trace_threads.rs `find_in_checkouts`).
fn find_in_checkouts(rel: &str) -> Option<PathBuf> {
    let local = Path::new(CRATE_DIR).join("../..").join(rel);
    if local.exists() {
        return Some(local);
    }
    let out = Command::new("git")
        .args(["rev-parse", "--git-common-dir"])
        .current_dir(CRATE_DIR)
        .output()
        .ok()?;
    let common = PathBuf::from(String::from_utf8(out.stdout).ok()?.trim());
    let common = if common.is_absolute() {
        common
    } else {
        Path::new(CRATE_DIR).join(common)
    };
    let main = common.parent()?.join(rel);
    main.exists().then_some(main)
}

fn miplib_instance(name: &str) -> Option<PathBuf> {
    find_in_checkouts(&format!("benches/data/miplib/{name}.mps"))
}

/// (status, objective, bound, nodes, lp_iters, time_s, x) of one `mip --csv` line.
type MipCsv = (String, f64, f64, i64, i64, f64, Vec<f64>);

/// Parse the CSV output of `rusolve-cli mip --csv`.
/// Returns (status, objective, bound, nodes, lp_iters, time_s, x_values).
fn parse_mip_csv(output: &str) -> MipCsv {
    // `rusolve-cli mip --csv` prints the HiGHS log (including the [HeurSol] and
    // [Sequential] lines) before the CSV header, so find the header line first.
    let lines: Vec<&str> = output.lines().collect();
    let header = lines
        .iter()
        .position(|l| l.starts_with("status,"))
        .unwrap_or_else(|| panic!("no CSV header in output: {output}"));
    assert!(lines.len() > header + 1, "CSV output too short: {output}");
    let fields: Vec<&str> = lines[header + 1].split(',').collect();
    assert!(
        fields.len() >= 7,
        "CSV line has too few fields: {}",
        lines[header + 1]
    );
    let status = fields[0].to_string();
    let objective: f64 = fields[1].parse().unwrap();
    let bound: f64 = fields[2].parse().unwrap();
    let nodes: i64 = fields[3].parse().unwrap();
    let lp_iters: i64 = fields[4].parse().unwrap();
    let time_s: f64 = fields[5].parse().unwrap();
    let x_joined = fields[6..].join(",");
    // the x field is `[x0;x1;...]`
    let x_str = x_joined.trim_start_matches('[').trim_end_matches(']');
    let x: Vec<f64> = x_str
        .split(';')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<f64>().unwrap())
        .collect();
    (status, objective, bound, nodes, lp_iters, time_s, x)
}

/// `rusolve-cli mip --csv` on the port; exit 0 iff the run ended Optimal.
fn run_mip_csv_unchecked(path: &Path, extra_args: &[&str]) -> (bool, String, MipCsv) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    // RUSOLVE_MIP_HIGHS routes `mip` to the HiGHS port, the only path that takes --opt
    cmd.env("RUSOLVE_MIP_HIGHS", "1");
    cmd.arg("mip").arg(path).arg("--csv").args(extra_args);
    let out = cmd.output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(
        stdout.lines().any(|l| l.starts_with("status,")),
        "mip produced no CSV: {stderr}"
    );
    (out.status.success(), stderr, parse_mip_csv(&stdout))
}

/// [`run_mip_csv_unchecked`], requiring status Optimal (TestMipSolver `solve()` requires
/// `kOptimal` before it checks the objective).
fn run_mip_csv(path: &Path, extra_args: &[&str]) -> MipCsv {
    let (ok, stderr, parsed) = run_mip_csv_unchecked(path, extra_args);
    assert!(
        ok && parsed.0 == "Optimal",
        "mip not Optimal ({}): {stderr}",
        parsed.0
    );
    parsed
}

/// TestMipSolver.cpp `objectiveOk`: relative error to max(1, |required|) below 1e-10.
fn objective_ok(objective: f64, required: f64) -> bool {
    (objective - required).abs() / required.abs().max(1.0) < 1e-10
}

/// Trace-exempt instances: the objective must agree with the fixture to the
/// MIP gap, HiGHS's default mip_rel_gap (1e-4), not bit for bit.
fn objective_within_mip_gap(objective: f64, required: f64) -> bool {
    (objective - required).abs() / required.abs().max(1.0) < 1e-4
}

/// Check feasibility of solution x against the original model read by rusolve_io.
fn check_feasibility(model_path: &Path, x: &[f64]) {
    let model = rusolve_io::read_mps_file(model_path.to_str().unwrap())
        .expect("failed to read MPS for feasibility check");

    let num_vars = model.num_variables;
    assert_eq!(x.len(), num_vars, "x length mismatch with model variables");

    // Check bounds and integrality
    for (i, &xi) in x.iter().enumerate() {
        let bound = &model.bounds[i];
        let tol = 1e-6;
        assert!(
            xi >= bound.lower - tol,
            "var {i}: x={xi} < lower={}",
            bound.lower
        );
        assert!(
            xi <= bound.upper + tol,
            "var {i}: x={xi} > upper={}",
            bound.upper
        );
        if bound.integer {
            let rounded = xi.round();
            assert!(
                (xi - rounded).abs() <= tol,
                "var {i}: x={xi} not integral (rounded={rounded})"
            );
        }
    }

    // Check row activities against rhs
    let matrix = &model.matrix;
    let nrow = matrix.nrow();
    let ncol = matrix.ncol();
    assert_eq!(ncol, num_vars, "matrix ncol != num_variables");

    let col_ptr = matrix.col_ptr();
    let row_idx = matrix.row_idx();
    let values = matrix.values();

    // Compute row norms (max |a_ij| per row)
    let mut row_norm = vec![0.0f64; nrow];
    for j in 0..ncol {
        let start = col_ptr[j];
        let end = col_ptr[j + 1];
        for k in start..end {
            let i = row_idx[k] as usize;
            let a = values[k].abs();
            if a > row_norm[i] {
                row_norm[i] = a;
            }
        }
    }

    // Compute row activities
    let mut activity = vec![0.0f64; nrow];
    for j in 0..ncol {
        let start = col_ptr[j];
        let end = col_ptr[j + 1];
        let xj = x[j];
        for k in start..end {
            let i = row_idx[k] as usize;
            activity[i] += values[k] * xj;
        }
    }

    // Check each row against [lo, hi], the MPS RANGES rules (rhs and ranges are raw)
    for i in 0..nrow {
        let rhs = model.rhs[i];
        let r = model.ranges[i];
        let (lo, hi) = match model.senses[i] {
            rusolve_io::RowSense::LessEqual if r != 0.0 => (rhs - r.abs(), rhs),
            rusolve_io::RowSense::LessEqual => (f64::NEG_INFINITY, rhs),
            rusolve_io::RowSense::GreaterEqual if r != 0.0 => (rhs, rhs + r.abs()),
            rusolve_io::RowSense::GreaterEqual => (rhs, f64::INFINITY),
            rusolve_io::RowSense::Equal if r > 0.0 => (rhs, rhs + r),
            rusolve_io::RowSense::Equal if r < 0.0 => (rhs + r, rhs),
            rusolve_io::RowSense::Equal => (rhs, rhs),
        };
        let tol = 1e-6 * row_norm[i].max(1.0);
        assert!(
            activity[i] >= lo - tol && activity[i] <= hi + tol,
            "row {i}: activity={} outside [{lo}, {hi}] (tol={tol})",
            activity[i]
        );
    }
}

// ---------------------------------------------------------------------------
// Helpers copied from mip_trace.rs
// ---------------------------------------------------------------------------

fn check_trace(instance: &str, fixture: &str, args: &[&str]) {
    run_trace(instance, fixture, args, false);
}

/// `check_trace`; with `crosscheck_ipm`, `RUSOLVE_MIP_AC_FIXTURE=<fixture>.ac` makes the port
/// compare every IPM result bit for bit with the oracle's (`ipm_replay_check`).
fn run_trace(instance: &str, fixture: &str, args: &[&str], crosscheck_ipm: bool) {
    let Some(mps) = miplib_instance(instance) else {
        eprintln!("skipping {fixture}: benches/data/miplib/{instance}.mps not found");
        return;
    };
    let dir = Path::new(CRATE_DIR).join("../rusolve-ekk/tests/fixtures/mip");
    let expected = std::fs::read_to_string(dir.join(format!("{fixture}.trace"))).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    cmd.arg("mip-trace").arg(&mps).args(args);
    if crosscheck_ipm {
        cmd.env("RUSOLVE_MIP_AC_FIXTURE", dir.join(format!("{fixture}.ac")));
    } else {
        cmd.env_remove("RUSOLVE_MIP_AC_FIXTURE");
    }
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "mip-trace failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let actual = String::from_utf8(out.stdout).unwrap();
    if let Some((i, (e, a))) = expected
        .lines()
        .zip(actual.lines())
        .enumerate()
        .find(|(_, (e, a))| e != a)
    {
        panic!(
            "{fixture}: first divergence at line {}\n  oracle: {e}\n  port:   {a}",
            i + 1
        );
    }
    assert_eq!(
        expected.lines().count(),
        actual.lines().count(),
        "{fixture}: trace lengths differ"
    );
}

// ---------------------------------------------------------------------------
// Helpers copied from mip_trace_threads.rs
// ---------------------------------------------------------------------------

/// `target/highs-fixtures/w5/mip/<file>` from this checkout or the main checkout.
fn w5_mip_trace(file: &str) -> Option<PathBuf> {
    find_in_checkouts(&format!("target/highs-fixtures/w5/mip/{file}"))
}

/// Skip banner, or a failure when `RUSOLVE_REQUIRE_TARGET_FIXTURES=1`.
fn skip(test: &str, reason: &str) {
    if std::env::var("RUSOLVE_REQUIRE_TARGET_FIXTURES").as_deref() == Ok("1") {
        panic!("{test}: {reason} (RUSOLVE_REQUIRE_TARGET_FIXTURES=1)");
    }
    eprintln!("SKIPPED, NOT PASSED: {test}: {reason}");
}

/// Diff `mip-trace <instance> <args>` against `expected`, the oracle trace text.
fn check_trace_w5(test: &str, instance: &str, expected: &str, args: &[&str]) {
    let Some(mps) = miplib_instance(instance) else {
        skip(
            test,
            &format!("benches/data/miplib/{instance}.mps not found"),
        );
        return;
    };
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    cmd.arg("mip-trace").arg(&mps).args(args);
    cmd.env_remove("RUSOLVE_MIP_AC_FIXTURE");
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "{test}: mip-trace failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let actual = String::from_utf8(out.stdout).unwrap();
    if let Some((i, (e, a))) = expected
        .lines()
        .zip(actual.lines())
        .enumerate()
        .find(|(_, (e, a))| e != a)
    {
        panic!(
            "{test}: first divergence at line {}\n  oracle: {e}\n  port:   {a}",
            i + 1
        );
    }
    assert_eq!(
        expected.lines().count(),
        actual.lines().count(),
        "{test}: trace lengths differ"
    );
}

/// The oracle trace `w5/mip/<instance>.t<threads>.trace`, or a skip.
fn oracle_w5(test: &str, instance: &str, threads: u32) -> Option<String> {
    let file = format!("{instance}.t{threads}.trace");
    let Some(path) = w5_mip_trace(&file) else {
        skip(
            test,
            &format!("target/highs-fixtures/w5/mip/{file} not found"),
        );
        return None;
    };
    Some(std::fs::read_to_string(path).unwrap())
}

// ---------------------------------------------------------------------------
// Local helpers
// ---------------------------------------------------------------------------

/// Mask helper: keep only lines starting with "[HeurSol] " or "[Sequential] ";
/// remove the fields " wall_ms=<token>" and " effort_per_ms=<token>" from each line.
fn mask_heur_lines(output: &str) -> Vec<String> {
    let mut result = Vec::new();
    for line in output.lines() {
        if !line.starts_with("[HeurSol] ") && !line.starts_with("[Sequential] ") {
            continue;
        }
        let mut masked = String::new();
        let mut chars = line.chars().peekable();
        while let Some(c) = chars.next() {
            if c == ' ' {
                // Check for " wall_ms="
                let rest: String = chars.clone().collect();
                if rest.starts_with("wall_ms=") {
                    // skip the token (non-whitespace chars after "wall_ms=")
                    for _ in "wall_ms=".chars() {
                        chars.next();
                    }
                    while chars.peek().is_some_and(|&c| c != ' ') {
                        chars.next();
                    }
                    continue;
                }
                if rest.starts_with("effort_per_ms=") {
                    for _ in "effort_per_ms=".chars() {
                        chars.next();
                    }
                    while chars.peek().is_some_and(|&c| c != ' ') {
                        chars.next();
                    }
                    continue;
                }
            }
            masked.push(c);
        }
        result.push(masked);
    }
    result
}

/// Fixture parser: read crates/rusolve-ekk/tests/fixtures/s5/<inst>.<cfg>.heur;
/// the lines before the last line are the expected masked trace; the last line
/// "final obj=<f64> status=<status>" gives the expected objective and status.
fn parse_heur_fixture(inst: &str, cfg: &str) -> Option<(Vec<String>, f64, String)> {
    let path = Path::new(CRATE_DIR)
        .join("../rusolve-ekk/tests/fixtures/s5")
        .join(format!("{inst}.{cfg}.heur"));
    if !path.exists() {
        return None;
    }
    let content = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert!(!lines.is_empty(), "fixture {inst}.{cfg}.heur is empty");
    let last_line = lines.last().unwrap();
    assert!(
        last_line.starts_with("final obj="),
        "fixture {inst}.{cfg}.heur last line does not start with 'final obj=': {last_line}"
    );
    // Parse "final obj=<f64> status=<status>"
    let parts: Vec<&str> = last_line.split_whitespace().collect();
    assert!(
        parts.len() >= 3,
        "fixture {inst}.{cfg}.heur last line has too few fields: {last_line}"
    );
    let obj_str = parts[1]
        .strip_prefix("obj=")
        .expect("no obj= in final line");
    let obj: f64 = obj_str.parse().unwrap();
    let status_str = parts[2]
        .strip_prefix("status=")
        .expect("no status= in final line");
    let status = status_str.to_string();
    let expected_trace: Vec<String> = lines[..lines.len() - 1]
        .iter()
        .map(|s| s.to_string())
        .collect();
    Some((expected_trace, obj, status))
}

/// Resolve the MPS path for an instance: first try s5 fixtures, then miplib.
fn s5_instance_path(name: &str) -> Option<PathBuf> {
    let s5_path = Path::new(CRATE_DIR)
        .join("../rusolve-ekk/tests/fixtures/s5")
        .join(format!("{name}.mps"));
    if s5_path.exists() {
        return Some(s5_path);
    }
    miplib_instance(name)
}

/// Run `rusolve-cli mip <mps> --csv` with heuristic options and return
/// (masked_trace_lines, MipCsv).
fn run_s5_mip(mps: &Path, run_fpr: bool, run_lm: bool) -> (Vec<String>, MipCsv) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    cmd.env("RUSOLVE_MIP_HIGHS", "1");
    cmd.arg("mip")
        .arg(mps)
        .arg("--csv")
        .arg("--opt")
        .arg("output_flag=true")
        .arg("--opt")
        .arg("log_dev_level=3")
        .arg("--opt")
        .arg("log_to_console=true")
        .arg("--opt")
        .arg("mip_heuristic_run_feasibility_jump=false")
        .arg("--opt")
        .arg("mip_first_inc=false")
        .arg("--opt")
        .arg("mip_first_inc_root=false")
        .arg("--opt")
        .arg(format!("mip_heuristic_run_fpr={run_fpr}"))
        .arg("--opt")
        .arg(format!("mip_heuristic_run_local_mip={run_lm}"));
    let out = cmd.output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(
        stdout.lines().any(|l| l.starts_with("status,")),
        "mip produced no CSV: {stderr}"
    );
    let csv = parse_mip_csv(&stdout);
    let masked = mask_heur_lines(&stdout);
    (masked, csv)
}

/// Core test logic for a flag-on integration test.
fn s5_it_core(inst: &str, cfg: &str, trace_exempt: bool) {
    let Some(mps) = s5_instance_path(inst) else {
        skip(
            &format!("s5_it_{inst}_{cfg}"),
            &format!("MPS for {inst} not found"),
        );
        return;
    };
    let Some((expected_trace, expected_obj, expected_status)) = parse_heur_fixture(inst, cfg)
    else {
        skip(
            &format!("s5_it_{inst}_{cfg}"),
            &format!("fixture {inst}.{cfg}.heur not found"),
        );
        return;
    };
    let (run_fpr, run_lm) = match cfg {
        "fpr" => (true, false),
        "lm" => (false, true),
        "both" => (true, true),
        _ => panic!("unknown cfg: {cfg}"),
    };
    let (masked, csv) = run_s5_mip(&mps, run_fpr, run_lm);

    if !trace_exempt {
        assert_eq!(
            masked, expected_trace,
            "s5_it_{inst}_{cfg}: masked trace mismatch"
        );
    }

    assert_eq!(
        csv.0, expected_status,
        "s5_it_{inst}_{cfg}: status mismatch (got {}, expected {expected_status})",
        csv.0
    );
    assert_eq!(
        csv.1.to_bits(),
        expected_obj.to_bits(),
        "s5_it_{inst}_{cfg}: objective bits mismatch (got {}, expected {expected_obj})",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

/// Core test logic for a trace-exempt instance (markshare_4_0, neos5, pk1):
/// check only status, objective (objective_ok) and feasibility.
fn s5_it_core_exempt(inst: &str, cfg: &str) {
    let Some(mps) = s5_instance_path(inst) else {
        skip(
            &format!("s5_it_{inst}_{cfg}"),
            &format!("MPS for {inst} not found"),
        );
        return;
    };
    let Some((_expected_trace, expected_obj, expected_status)) = parse_heur_fixture(inst, cfg)
    else {
        skip(
            &format!("s5_it_{inst}_{cfg}"),
            &format!("fixture {inst}.{cfg}.heur not found"),
        );
        return;
    };
    let (run_fpr, run_lm) = match cfg {
        "fpr" => (true, false),
        "lm" => (false, true),
        "both" => (true, true),
        _ => panic!("unknown cfg: {cfg}"),
    };
    let (_masked, csv) = run_s5_mip(&mps, run_fpr, run_lm);

    assert_eq!(
        csv.0, expected_status,
        "s5_it_{inst}_{cfg}: status mismatch (got {}, expected {expected_status})",
        csv.0
    );
    assert!(
        objective_within_mip_gap(csv.1, expected_obj),
        "s5_it_{inst}_{cfg}: objective mismatch (got {}, expected {expected_obj})",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

// ---------------------------------------------------------------------------
// Flag-on tests: e2e_facility (s5 fixture)
// ---------------------------------------------------------------------------

#[test]
fn s5_it_e2e_facility_fpr() {
    s5_it_core("e2e_facility", "fpr", false);
}

#[test]
fn s5_it_e2e_facility_lm() {
    s5_it_core("e2e_facility", "lm", false);
}

#[test]
fn s5_it_e2e_facility_both() {
    s5_it_core("e2e_facility", "both", false);
}

// ---------------------------------------------------------------------------
// Flag-on tests: e2e_knap (s5 fixture)
// ---------------------------------------------------------------------------

#[test]
fn s5_it_e2e_knap_fpr() {
    s5_it_core("e2e_knap", "fpr", false);
}

#[test]
fn s5_it_e2e_knap_lm() {
    s5_it_core("e2e_knap", "lm", false);
}

#[test]
fn s5_it_e2e_knap_both() {
    s5_it_core("e2e_knap", "both", false);
}

// ---------------------------------------------------------------------------
// Flag-on tests: e2e_setcover (s5 fixture)
// ---------------------------------------------------------------------------

#[test]
fn s5_it_e2e_setcover_fpr() {
    s5_it_core("e2e_setcover", "fpr", false);
}

#[test]
fn s5_it_e2e_setcover_lm() {
    s5_it_core("e2e_setcover", "lm", false);
}

#[test]
fn s5_it_e2e_setcover_both() {
    s5_it_core("e2e_setcover", "both", false);
}

// ---------------------------------------------------------------------------
// Flag-on tests: gen-ip002 (miplib, slow)
// ---------------------------------------------------------------------------

#[test]
#[ignore = "gen-ip002 slow"]
fn s5_it_gen_ip002_fpr() {
    s5_it_core("gen-ip002", "fpr", false);
}

#[test]
#[ignore = "gen-ip002 slow"]
fn s5_it_gen_ip002_lm() {
    s5_it_core("gen-ip002", "lm", false);
}

#[test]
#[ignore = "gen-ip002 slow"]
fn s5_it_gen_ip002_both() {
    s5_it_core("gen-ip002", "both", false);
}

// ---------------------------------------------------------------------------
// Flag-on tests: mas74 (miplib, slow)
// ---------------------------------------------------------------------------

#[test]
#[ignore = "mas74 slow"]
fn s5_it_mas74_fpr() {
    s5_it_core("mas74", "fpr", false);
}

#[test]
#[ignore = "mas74 slow"]
fn s5_it_mas74_lm() {
    s5_it_core("mas74", "lm", false);
}

#[test]
#[ignore = "mas74 slow"]
fn s5_it_mas74_both() {
    s5_it_core("mas74", "both", false);
}

// ---------------------------------------------------------------------------
// Flag-on tests: pg (miplib, slow)
// ---------------------------------------------------------------------------

#[test]
#[ignore = "pg slow"]
fn s5_it_pg_fpr() {
    s5_it_core("pg", "fpr", false);
}

#[test]
#[ignore = "pg slow"]
fn s5_it_pg_lm() {
    s5_it_core("pg", "lm", false);
}

#[test]
#[ignore = "pg slow"]
fn s5_it_pg_both() {
    s5_it_core("pg", "both", false);
}

// ---------------------------------------------------------------------------
// Flag-on tests: markshare_4_0 (miplib, slow, trace-exempt)
// ---------------------------------------------------------------------------

#[test]
#[ignore = "markshare_4_0 slow"]
fn s5_it_markshare_4_0_fpr() {
    s5_it_core_exempt("markshare_4_0", "fpr");
}

#[test]
#[ignore = "markshare_4_0 slow"]
fn s5_it_markshare_4_0_lm() {
    s5_it_core_exempt("markshare_4_0", "lm");
}

#[test]
#[ignore = "markshare_4_0 slow"]
fn s5_it_markshare_4_0_both() {
    s5_it_core_exempt("markshare_4_0", "both");
}

// ---------------------------------------------------------------------------
// Flag-on tests: neos5 (miplib, slow, trace-exempt)
// ---------------------------------------------------------------------------

#[test]
#[ignore = "neos5 slow"]
fn s5_it_neos5_fpr() {
    s5_it_core_exempt("neos5", "fpr");
}

#[test]
#[ignore = "neos5 slow"]
fn s5_it_neos5_lm() {
    s5_it_core_exempt("neos5", "lm");
}

#[test]
#[ignore = "neos5 slow"]
fn s5_it_neos5_both() {
    s5_it_core_exempt("neos5", "both");
}

// ---------------------------------------------------------------------------
// Flag-on tests: pk1 (miplib, slow, trace-exempt)
// ---------------------------------------------------------------------------

#[test]
#[ignore = "pk1 slow"]
fn s5_it_pk1_fpr() {
    s5_it_core_exempt("pk1", "fpr");
}

#[test]
#[ignore = "pk1 slow"]
fn s5_it_pk1_lm() {
    s5_it_core_exempt("pk1", "lm");
}

#[test]
#[ignore = "pk1 slow"]
fn s5_it_pk1_both() {
    s5_it_core_exempt("pk1", "both");
}

// ---------------------------------------------------------------------------
// FJ-on tests: default FJ (do not pass mip_heuristic_run_feasibility_jump),
// both flags true, check only the optimum (objective_ok) and feasibility.
// ---------------------------------------------------------------------------

const S5_BOTH_ON: [&str; 4] = [
    "--opt",
    "mip_heuristic_run_fpr=true",
    "--opt",
    "mip_heuristic_run_local_mip=true",
];

fn s5_it_fj_on(test: &str, instance: &str, optimum: f64) {
    let Some(mps) = miplib_instance(instance) else {
        skip(
            test,
            &format!("benches/data/miplib/{instance}.mps not found"),
        );
        return;
    };
    let csv = run_mip_csv(&mps, &S5_BOTH_ON);
    assert!(
        objective_ok(csv.1, optimum),
        "{test}: objective mismatch (got {}, expected {optimum})",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
#[ignore = "pg slow"]
fn s5_it_fj_on_pg() {
    s5_it_fj_on("s5_it_fj_on_pg", "pg", -8674.34260712);
}

#[test]
#[ignore = "sp150x300d slow"]
fn s5_it_fj_on_sp150x300d() {
    s5_it_fj_on("s5_it_fj_on_sp150x300d", "sp150x300d", 69.0);
}

#[test]
#[ignore = "air05 slow"]
fn s5_it_fj_on_air05() {
    s5_it_fj_on("s5_it_fj_on_air05", "air05", 26374.0);
}

#[test]
#[ignore = "enlight_hard slow"]
fn s5_it_fj_on_enlight_hard() {
    s5_it_fj_on("s5_it_fj_on_enlight_hard", "enlight_hard", 37.0);
}

// ---------------------------------------------------------------------------
// Flag-off tests: both flags false, FJ off, mip-trace compared to existing
// committed trace fixtures.
// ---------------------------------------------------------------------------

#[test]
fn s5_it_off_sp150x300d() {
    check_trace(
        "sp150x300d",
        "sp150x300d",
        &[
            "mip_heuristic_run_fpr=false",
            "mip_heuristic_run_local_mip=false",
        ],
    );
}

#[test]
#[ignore = "s5_it_off_pg slow"]
fn s5_it_off_pg() {
    check_trace(
        "pg",
        "pg",
        &[
            "mip_heuristic_run_fpr=false",
            "mip_heuristic_run_local_mip=false",
        ],
    );
}

#[test]
#[ignore = "s5_it_off_air05 slow"]
fn s5_it_off_air05() {
    check_trace(
        "air05",
        "air05",
        &[
            "mip_heuristic_run_fpr=false",
            "mip_heuristic_run_local_mip=false",
        ],
    );
}

#[test]
#[ignore = "s5_it_off_enlight_hard slow"]
fn s5_it_off_enlight_hard() {
    check_trace(
        "enlight_hard",
        "enlight_hard",
        &[
            "mip_heuristic_run_fpr=false",
            "mip_heuristic_run_local_mip=false",
        ],
    );
}

// Both flags false, threads = 4, parallel = on, against the multi-worker oracle trace
// (mip_trace_threads.rs gates sp150x300d and pg; the air05 / enlight_hard traces are
// ignored there).

fn s5_it_off_w5(test: &str, instance: &str) {
    let Some(expected) = oracle_w5(test, instance, 4) else {
        return;
    };
    check_trace_w5(
        test,
        instance,
        &expected,
        &[
            "threads=4",
            "parallel=on",
            "mip_heuristic_run_fpr=false",
            "mip_heuristic_run_local_mip=false",
        ],
    );
}

#[test]
#[ignore = "s5_it_off_w5_sp150x300d_t4 slow"]
fn s5_it_off_w5_sp150x300d_t4() {
    s5_it_off_w5("s5_it_off_w5_sp150x300d_t4", "sp150x300d");
}

#[test]
#[ignore = "s5_it_off_w5_pg_t4 slow"]
fn s5_it_off_w5_pg_t4() {
    s5_it_off_w5("s5_it_off_w5_pg_t4", "pg");
}
