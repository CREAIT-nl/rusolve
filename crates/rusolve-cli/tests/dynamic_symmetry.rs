//! Integration tests for mip_symmetry_dynamic (SCIP 11.0.0 dynamic symmetry handling, Apache-2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB))
//!
//! This file is not edited for the other tests (highs_backport_options.rs, mip_trace.rs, mip_trace_threads.rs, conflict_resolution.rs).

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

fn mip50_instance(name: &str) -> Option<PathBuf> {
    find_in_checkouts(&format!("benches/data/mip50/{name}.mps"))
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

fn rel_close(a: f64, b: f64, tol: f64) -> bool {
    if a == b {
        return true;
    }
    let diff = (a - b).abs();
    let scale = if b.abs() > 1.0 { b.abs() } else { 1.0 };
    diff <= tol * scale
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
// Node limit for node-limited tests
// ---------------------------------------------------------------------------

const NODE_LIMIT: u32 = 500;

// ---------------------------------------------------------------------------
// Flag on: optimum kept and x feasible (MIPLIB instances)
// ---------------------------------------------------------------------------

#[test]
#[ignore = "s3_it_on_pg slow"]
fn s3_it_on_pg() {
    let mps = miplib_instance("pg").expect("pg.mps not found");
    let csv = run_mip_csv(&mps, &["--opt", "mip_symmetry_dynamic=true"]);
    assert!(
        rel_close(csv.1, -8674.34260712, 1e-6),
        "pg objective {} not close to -8674.34260712",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
fn s3_it_on_sp150x300d() {
    let mps = miplib_instance("sp150x300d").expect("sp150x300d.mps not found");
    let csv = run_mip_csv(&mps, &["--opt", "mip_symmetry_dynamic=true"]);
    assert!(
        rel_close(csv.1, 69.0, 1e-6),
        "sp150x300d objective {} not close to 69",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
#[ignore = "s3_it_on_air05 slow"]
fn s3_it_on_air05() {
    let mps = miplib_instance("air05").expect("air05.mps not found");
    let csv = run_mip_csv(&mps, &["--opt", "mip_symmetry_dynamic=true"]);
    assert!(
        rel_close(csv.1, 26374.0, 1e-6),
        "air05 objective {} not close to 26374",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
#[ignore = "s3_it_on_enlight_hard slow"]
fn s3_it_on_enlight_hard() {
    let mps = miplib_instance("enlight_hard").expect("enlight_hard.mps not found");
    let csv = run_mip_csv(&mps, &["--opt", "mip_symmetry_dynamic=true"]);
    assert!(
        rel_close(csv.1, 37.0, 1e-6),
        "enlight_hard objective {} not close to 37",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

// ---------------------------------------------------------------------------
// Flag on: neos instances
// ---------------------------------------------------------------------------

#[test]
#[ignore = "s3_it_on_neos911970 slow"]
fn s3_it_on_neos911970() {
    let mps = miplib_instance("neos-911970").expect("neos-911970.mps not found");
    let csv = run_mip_csv(&mps, &["--opt", "mip_symmetry_dynamic=true"]);
    assert!(
        rel_close(csv.1, 54.76, 1e-6),
        "neos-911970 objective {} not close to 54.76",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
#[ignore = "s3_it_on_neos5 slow"]
fn s3_it_on_neos5() {
    let mps = miplib_instance("neos5").expect("neos5.mps not found");
    let csv = run_mip_csv(&mps, &["--opt", "mip_symmetry_dynamic=true"]);
    assert!(
        rel_close(csv.1, 15.0, 1e-6),
        "neos5 objective {} not close to 15",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

// ---------------------------------------------------------------------------
// Both mip_conflict_resolution=true and mip_symmetry_dynamic=true
// ---------------------------------------------------------------------------

#[test]
#[ignore = "s3_it_both_pg slow"]
fn s3_it_both_pg() {
    let mps = miplib_instance("pg").expect("pg.mps not found");
    let csv = run_mip_csv(
        &mps,
        &[
            "--opt",
            "mip_conflict_resolution=true",
            "--opt",
            "mip_symmetry_dynamic=true",
        ],
    );
    assert!(
        rel_close(csv.1, -8674.34260712, 1e-6),
        "pg objective {} not close to -8674.34260712",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
fn s3_it_both_sp150x300d() {
    let mps = miplib_instance("sp150x300d").expect("sp150x300d.mps not found");
    let csv = run_mip_csv(
        &mps,
        &[
            "--opt",
            "mip_conflict_resolution=true",
            "--opt",
            "mip_symmetry_dynamic=true",
        ],
    );
    assert!(
        rel_close(csv.1, 69.0, 1e-6),
        "sp150x300d objective {} not close to 69",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

#[test]
#[ignore = "s3_it_both_air05 slow"]
fn s3_it_both_air05() {
    let mps = miplib_instance("air05").expect("air05.mps not found");
    let csv = run_mip_csv(
        &mps,
        &[
            "--opt",
            "mip_conflict_resolution=true",
            "--opt",
            "mip_symmetry_dynamic=true",
        ],
    );
    assert!(
        rel_close(csv.1, 26374.0, 1e-6),
        "air05 objective {} not close to 26374",
        csv.1
    );
    check_feasibility(&mps, &csv.6);
}

// ---------------------------------------------------------------------------
// Node-limited bounds (mip50 instances, flag on, minimization)
// ---------------------------------------------------------------------------

/// Flag on, node limit [`NODE_LIMIT`]: dual bound <= opt + tol and, when an incumbent
/// exists (the CSV objective is NaN without one), primal >= opt - tol.
fn check_node_limited(name: &str, opt: f64) {
    let mps = mip50_instance(name).unwrap_or_else(|| panic!("{name}.mps not found"));
    let tol = 1e-6 * opt.abs().max(1.0);
    let (_, stderr, csv) = run_mip_csv_unchecked(
        &mps,
        &[
            "--opt",
            "mip_symmetry_dynamic=true",
            "--node-limit",
            &NODE_LIMIT.to_string(),
        ],
    );
    assert!(
        csv.2 <= opt + tol,
        "{name} bound {} > opt+tol ({}): {stderr}",
        csv.2,
        opt + tol
    );
    if !csv.1.is_nan() {
        assert!(
            csv.1 >= opt - tol,
            "{name} objective {} < opt-tol ({}): {stderr}",
            csv.1,
            opt - tol
        );
    }
}

#[test]
#[ignore = "s3_it_nodelim_cvs16r128_89 slow"]
fn s3_it_nodelim_cvs16r128_89() {
    check_node_limited("cvs16r128-89", -97.0);
}

#[test]
#[ignore = "s3_it_nodelim_wachplan slow"]
fn s3_it_nodelim_wachplan() {
    check_node_limited("wachplan", -8.0);
}

#[test]
#[ignore = "s3_it_nodelim_fastxgemm_n2r6s0t2 slow"]
fn s3_it_nodelim_fastxgemm_n2r6s0t2() {
    check_node_limited("fastxgemm-n2r6s0t2", 230.0);
}

#[test]
#[ignore = "s3_it_nodelim_brazil3 slow"]
fn s3_it_nodelim_brazil3() {
    check_node_limited("brazil3", 24.0);
}

#[test]
#[ignore = "s3_it_nodelim_comp07_2idx slow"]
fn s3_it_nodelim_comp07_2idx() {
    check_node_limited("comp07-2idx", 6.0);
}

#[test]
#[ignore = "s3_it_nodelim_fast0507 slow"]
fn s3_it_nodelim_fast0507() {
    check_node_limited("fast0507", 174.0);
}

#[test]
#[ignore = "s3_it_nodelim_seymour slow"]
fn s3_it_nodelim_seymour() {
    check_node_limited("seymour", 423.0);
}

// ---------------------------------------------------------------------------
// Report: off vs on objective identity and node counts
// ---------------------------------------------------------------------------

#[test]
#[ignore = "s3_it_report slow"]
fn s3_it_report() {
    let instances: &[(&str, f64)] = &[
        ("cvs16r128-89", -97.0),
        ("wachplan", -8.0),
        ("fastxgemm-n2r6s0t2", 230.0),
        ("brazil3", 24.0),
        ("comp07-2idx", 6.0),
        ("fast0507", 174.0),
        ("seymour", 423.0),
    ];

    let node_limit_str = NODE_LIMIT.to_string();

    eprintln!(
        "{:<22} {:>12} {:>12} {:>12} {:>9} {:>10} {:>10}",
        "instance", "opt", "off obj", "on obj", "obj ident", "off nodes", "on nodes"
    );
    eprintln!("{}", "-".repeat(93));

    for (name, opt) in instances {
        let mps = mip50_instance(name).unwrap_or_else(|| panic!("{name}.mps not found"));

        let (_, _, off_csv) = run_mip_csv_unchecked(
            &mps,
            &[
                "--opt",
                "mip_symmetry_dynamic=false",
                "--node-limit",
                &node_limit_str,
            ],
        );

        let (_, _, on_csv) = run_mip_csv_unchecked(
            &mps,
            &[
                "--opt",
                "mip_symmetry_dynamic=true",
                "--node-limit",
                &node_limit_str,
            ],
        );

        // bit identity (NaN == NaN when neither run has an incumbent)
        let ident = off_csv.1.to_bits() == on_csv.1.to_bits();
        eprintln!(
            "{:<22} {:>12.6} {:>12.6} {:>12.6} {:>9} {:>10} {:>10}",
            name,
            opt,
            off_csv.1,
            on_csv.1,
            if ident { "yes" } else { "no" },
            off_csv.3,
            on_csv.3
        );
    }
}
