//! Integration tests for mip_cut_selection_hybrid (SCIP 11.0.0 hybrid cut selection, Apache-2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB))
//!
//! This file is not edited for the other tests (highs_backport_options.rs, mip_trace.rs, mip_trace_threads.rs).

use std::path::{Path, PathBuf};
use std::process::Command;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

// ---------------------------------------------------------------------------
// Helpers copied from highs_backport_options.rs
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
    let lines: Vec<&str> = output.lines().collect();
    assert!(lines.len() >= 2, "CSV output too short: {output}");
    let fields: Vec<&str> = lines[1].split(',').collect();
    assert!(
        fields.len() >= 7,
        "CSV line has too few fields: {}",
        lines[1]
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
        stdout.starts_with("status,"),
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

/// Parse the final line of `rusolve-cli mip-trace` output.
/// Requires `status=Optimal`; returns the objective value from `obj=<value>`.
fn parse_mip_trace_final_obj(output: &str) -> f64 {
    let last_line = output.lines().last().expect("mip-trace produced no output");
    assert!(
        last_line.starts_with("final status=Optimal "),
        "last line is not an Optimal final status line: {last_line}"
    );
    let obj_str = last_line
        .split_whitespace()
        .find(|s| s.starts_with("obj="))
        .expect("no obj= in final line");
    obj_str[4..].parse::<f64>().unwrap()
}

fn run_mip_trace(path: &Path, time_limit: &str, extra_args: &[&str]) -> f64 {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    cmd.arg("mip-trace")
        .arg(path)
        .arg(time_limit)
        .args(extra_args);
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "mip-trace failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    parse_mip_trace_final_obj(&stdout)
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
// Flag off: mip-trace byte-identical to fixtures with mip_cut_selection_hybrid=false
// ---------------------------------------------------------------------------

#[test]
fn s4_it_off_sp150x300d() {
    check_trace(
        "sp150x300d",
        "sp150x300d",
        &["mip_cut_selection_hybrid=false"],
    );
}

#[test]
fn s4_it_off_sp150x300d_presolve_off() {
    check_trace(
        "sp150x300d",
        "sp150x300d.presolve_off",
        &["presolve=off", "mip_cut_selection_hybrid=false"],
    );
}

#[test]
#[ignore = "s4_it_off_pg slow"]
fn s4_it_off_pg() {
    check_trace("pg", "pg", &["mip_cut_selection_hybrid=false"]);
}

#[test]
#[ignore = "s4_it_off_air05 slow"]
fn s4_it_off_air05() {
    check_trace("air05", "air05", &["mip_cut_selection_hybrid=false"]);
}

#[test]
#[ignore = "s4_it_off_enlight_hard slow"]
fn s4_it_off_enlight_hard() {
    check_trace(
        "enlight_hard",
        "enlight_hard",
        &["mip_cut_selection_hybrid=false"],
    );
}

#[test]
#[ignore = "s4_it_off_w5_sp150x300d_t4 slow"]
fn s4_it_off_w5_sp150x300d_t4() {
    let expected = match oracle_w5("s4_it_off_w5_sp150x300d_t4", "sp150x300d", 4) {
        Some(e) => e,
        None => return,
    };
    check_trace_w5(
        "s4_it_off_w5_sp150x300d_t4",
        "sp150x300d",
        &expected,
        &["threads=4", "parallel=on", "mip_cut_selection_hybrid=false"],
    );
}

// ---------------------------------------------------------------------------
// Flag on: objective and feasibility
// ---------------------------------------------------------------------------

#[test]
fn s4_it_on_sp150x300d() {
    let mps = miplib_instance("sp150x300d").expect("sp150x300d.mps not found");

    // mip-trace with mip_cut_selection_hybrid=true, objective 69
    let obj_trace = run_mip_trace(&mps, "3600", &["mip_cut_selection_hybrid=true"]);
    assert!(
        objective_ok(obj_trace, 69.0),
        "mip-trace objective {obj_trace} != 69"
    );

    // csv with --opt mip_cut_selection_hybrid=true, objective 69, plus feasibility
    let csv = run_mip_csv(&mps, &["--opt", "mip_cut_selection_hybrid=true"]);
    assert!(objective_ok(csv.1, 69.0), "csv objective {} != 69", csv.1);
    check_feasibility(&mps, &csv.6);
}

#[test]
#[ignore = "s4_it_on_pg slow"]
fn s4_it_on_pg() {
    let mps = miplib_instance("pg").expect("pg.mps not found");
    let obj = run_mip_trace(&mps, "3600", &["mip_cut_selection_hybrid=true"]);
    assert!(
        objective_ok(obj, -8674.34260712),
        "mip-trace objective {obj} != -8674.34260712"
    );
    let csv = run_mip_csv(&mps, &["--opt", "mip_cut_selection_hybrid=true"]);
    assert!(
        objective_ok(csv.1, -8674.34260712),
        "csv objective {} != -8674.34260712",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
#[ignore = "s4_it_on_air05 slow"]
fn s4_it_on_air05() {
    let mps = miplib_instance("air05").expect("air05.mps not found");
    let obj = run_mip_trace(&mps, "3600", &["mip_cut_selection_hybrid=true"]);
    assert!(
        objective_ok(obj, 26374.0),
        "mip-trace objective {obj} != 26374"
    );
    let csv = run_mip_csv(&mps, &["--opt", "mip_cut_selection_hybrid=true"]);
    assert!(
        objective_ok(csv.1, 26374.0),
        "csv objective {} != 26374",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
#[ignore = "s4_it_on_enlight_hard slow"]
fn s4_it_on_enlight_hard() {
    let mps = miplib_instance("enlight_hard").expect("enlight_hard.mps not found");
    let obj = run_mip_trace(&mps, "3600", &["mip_cut_selection_hybrid=true"]);
    assert!(objective_ok(obj, 37.0), "mip-trace objective {obj} != 37");
    let csv = run_mip_csv(&mps, &["--opt", "mip_cut_selection_hybrid=true"]);
    assert!(objective_ok(csv.1, 37.0), "csv objective {} != 37", csv.1);
    check_feasibility(&mps, &csv.6);
}
