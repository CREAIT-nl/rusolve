//! HiGHS backport integration tests: PR tests #3220 (merge b7d4e98c5), #3223 (merge 7db1c8180),
//! #3209 (merge f7b87ae01), #3244 (merge 8a5f9c1d8).
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! (1) each option of this set on alone and all 3 together on pg,
//! sp150x300d, air05 and enlight_hard (optimum and feasibility on the original model);
//! (2) keep_implications together with dual_fix_probing and aggregate_implications;
//! (3) traininstance6 with cdouble_proofs (skip if absent); (4) the eligibility report
//! (flag-on mip-trace byte-compared to the 1.15.1 fixtures, printed, never a failure).

use std::path::{Path, PathBuf};
use std::process::Command;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

const X1B_OPTIONS: &[&str] = &[
    "mip_cdouble_proofs",
    "mip_fix_last_active",
    "mip_keep_implications",
];

fn miplib_instance(name: &str) -> Option<PathBuf> {
    repo_file(&format!("benches/data/miplib/{name}.mps"))
}

/// `rel` in this worktree, else in the main checkout (gitignored data lives there).
fn repo_file(rel: &str) -> Option<PathBuf> {
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

/// Run mip-trace with a single option on, check objective and feasibility via mip --csv.
fn check_option_on_instance(inst_name: &str, expected_obj: f64, option: &str) {
    let path =
        miplib_instance(inst_name).unwrap_or_else(|| panic!("instance {inst_name} not found"));
    let opt_arg = format!("{option}=true");
    let obj = run_mip_trace(&path, "3600", &[&opt_arg]);
    assert!(
        rel_close(obj, expected_obj, 1e-6),
        "{inst_name} with {option}=true: obj={obj}, expected={expected_obj}"
    );
    let csv = run_mip_csv(&path, &["--opt", &opt_arg]);
    assert!(
        rel_close(csv.1, expected_obj, 1e-6),
        "{inst_name} mip --csv with {option}=true: obj={}, expected={expected_obj}",
        csv.1
    );
    check_feasibility(&path, &csv.6);
}

/// Run mip-trace with all given options on, check objective and feasibility via mip --csv.
fn check_all_options_on_instance(inst_name: &str, expected_obj: f64, options: &[&str]) {
    let path =
        miplib_instance(inst_name).unwrap_or_else(|| panic!("instance {inst_name} not found"));
    let opt_args: Vec<String> = options.iter().map(|o| format!("{o}=true")).collect();
    let opt_refs: Vec<&str> = opt_args.iter().map(|s| s.as_str()).collect();
    let obj = run_mip_trace(&path, "3600", &opt_refs);
    assert!(
        rel_close(obj, expected_obj, 1e-6),
        "{inst_name} with all options on: obj={obj}, expected={expected_obj}"
    );
    let mut csv_args: Vec<&str> = Vec::new();
    for s in &opt_args {
        csv_args.push("--opt");
        csv_args.push(s.as_str());
    }
    let csv = run_mip_csv(&path, &csv_args);
    assert!(
        rel_close(csv.1, expected_obj, 1e-6),
        "{inst_name} mip --csv with all options on: obj={}, expected={expected_obj}",
        csv.1
    );
    check_feasibility(&path, &csv.6);
}

/// Helper: run mip-trace with keep_implications plus extra options, check objective
/// and feasibility.
fn check_b7_instance(inst_name: &str, expected_obj: f64, extra_options: &[&str]) {
    let path =
        miplib_instance(inst_name).unwrap_or_else(|| panic!("instance {inst_name} not found"));
    let mut opt_args: Vec<String> = vec!["mip_keep_implications=true".to_string()];
    for o in extra_options {
        opt_args.push(format!("{o}=true"));
    }
    let opt_refs: Vec<&str> = opt_args.iter().map(|s| s.as_str()).collect();
    let obj = run_mip_trace(&path, "3600", &opt_refs);
    assert!(
        rel_close(obj, expected_obj, 1e-6),
        "{inst_name} keep_implications: obj={obj}, expected={expected_obj}"
    );
    let mut csv_args: Vec<&str> = Vec::new();
    for s in &opt_args {
        csv_args.push("--opt");
        csv_args.push(s.as_str());
    }
    let csv = run_mip_csv(&path, &csv_args);
    assert!(
        rel_close(csv.1, expected_obj, 1e-6),
        "{inst_name} keep_implications mip --csv: obj={}, expected={expected_obj}",
        csv.1
    );
    check_feasibility(&path, &csv.6);
}

// ---------------------------------------------------------------------------
// 1. Each HiGHS backport option on alone
// ---------------------------------------------------------------------------

#[test]
#[ignore = "slow"]
fn x1b_it_each_option_on_pg() {
    for &opt in X1B_OPTIONS {
        check_option_on_instance("pg", -8674.34260712, opt);
    }
}

#[test]
#[ignore = "slow"]
fn x1b_it_each_option_on_sp150x300d() {
    for &opt in X1B_OPTIONS {
        check_option_on_instance("sp150x300d", 69.0, opt);
    }
}

#[test]
#[ignore = "slow"]
fn x1b_it_each_option_on_air05() {
    for &opt in X1B_OPTIONS {
        check_option_on_instance("air05", 26374.0, opt);
    }
}

#[test]
#[ignore = "slow"]
fn x1b_it_each_option_on_enlight_hard() {
    for &opt in X1B_OPTIONS {
        check_option_on_instance("enlight_hard", 37.0, opt);
    }
}

// ---------------------------------------------------------------------------
// 2. All HiGHS backport options on together
// ---------------------------------------------------------------------------

#[test]
#[ignore = "slow"]
fn x1b_it_all_options_on_pg() {
    check_all_options_on_instance("pg", -8674.34260712, X1B_OPTIONS);
}

#[test]
#[ignore = "slow"]
fn x1b_it_all_options_on_sp150x300d() {
    check_all_options_on_instance("sp150x300d", 69.0, X1B_OPTIONS);
}

#[test]
#[ignore = "slow"]
fn x1b_it_all_options_on_air05() {
    check_all_options_on_instance("air05", 26374.0, X1B_OPTIONS);
}

#[test]
#[ignore = "slow"]
fn x1b_it_all_options_on_enlight_hard() {
    check_all_options_on_instance("enlight_hard", 37.0, X1B_OPTIONS);
}

// ---------------------------------------------------------------------------
// 3. keep_implications with dual_fix_probing and aggregate_implications
// ---------------------------------------------------------------------------

#[test]
#[ignore = "slow"]
fn x1b_it_keep_implications_b7_pg() {
    check_b7_instance(
        "pg",
        -8674.34260712,
        &["mip_dual_fix_probing", "mip_aggregate_implications"],
    );
}

#[test]
#[ignore = "slow"]
fn x1b_it_keep_implications_b7_sp150x300d() {
    check_b7_instance(
        "sp150x300d",
        69.0,
        &["mip_dual_fix_probing", "mip_aggregate_implications"],
    );
}

#[test]
#[ignore = "slow"]
fn x1b_it_keep_implications_b7_air05() {
    check_b7_instance(
        "air05",
        26374.0,
        &["mip_dual_fix_probing", "mip_aggregate_implications"],
    );
}

#[test]
#[ignore = "slow"]
fn x1b_it_keep_implications_b7_enlight_hard() {
    check_b7_instance(
        "enlight_hard",
        37.0,
        &["mip_dual_fix_probing", "mip_aggregate_implications"],
    );
}

// ---------------------------------------------------------------------------
// 4. traininstance6 with cdouble_proofs (skip if absent)
// ---------------------------------------------------------------------------

#[test]
fn x1b_it_traininstance6() {
    // HiGHS backport-IT: the instance is looked up under reference/HiGHS-latest/check/instances/
    // (where the job places it), then under benches/data/miplib/.
    let path = match repo_file("reference/HiGHS-latest/check/instances/traininstance6.mps")
        .or_else(|| miplib_instance("traininstance6"))
    {
        Some(p) => p,
        None => {
            eprintln!("x1b_it_traininstance6: traininstance6.mps not found, skipping");
            return;
        }
    };
    let obj = run_mip_trace(&path, "3600", &["mip_cdouble_proofs=true"]);
    eprintln!("x1b_it_traininstance6: mip_cdouble_proofs=true obj={obj}");
}

// ---------------------------------------------------------------------------
// 5. Eligibility report: flag-on mip-trace vs 1.15.1 fixtures (report only)
// ---------------------------------------------------------------------------

#[test]
#[ignore = "slow"]
fn x1b_it_eligibility_report() {
    let cases: &[(&str, &str, &[&str])] = &[
        ("pg", "pg", &[]),
        ("sp150x300d", "sp150x300d", &[]),
        ("sp150x300d", "sp150x300d.presolve_off", &["presolve=off"]),
    ];
    let flags = &["mip_cdouble_proofs", "mip_fix_last_active"];
    for &flag in flags {
        for &(inst, fixture_name, extra) in cases {
            let path = match miplib_instance(inst) {
                Some(p) => p,
                None => {
                    eprintln!("x1b_it_eligibility_report: instance {inst} not found, skipping");
                    continue;
                }
            };
            let fixture_path = Path::new(CRATE_DIR)
                .join("../rusolve-ekk/tests/fixtures/mip")
                .join(format!("{fixture_name}.trace"));
            if !fixture_path.exists() {
                eprintln!(
                    "x1b_it_eligibility_report: fixture {} not found, skipping",
                    fixture_path.display()
                );
                continue;
            }
            let opt_arg = format!("{flag}=true");
            let mut args: Vec<&str> = vec![opt_arg.as_str()];
            args.extend_from_slice(extra);

            let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
            cmd.arg("mip-trace").arg(&path).arg("3600").args(&args);
            let out = cmd.output().unwrap();
            assert!(
                out.status.success(),
                "mip-trace failed for {inst} with {flag}: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            let stdout = String::from_utf8(out.stdout).unwrap();

            let fixture_content = std::fs::read_to_string(&fixture_path)
                .unwrap_or_else(|e| panic!("cannot read fixture {}: {e}", fixture_path.display()));

            let stdout_lines: Vec<&str> = stdout.lines().collect();
            let fixture_lines: Vec<&str> = fixture_content.lines().collect();

            let identical = stdout_lines == fixture_lines;
            if identical {
                eprintln!(
                    "x1b_it_eligibility_report: {flag} on {inst} ({fixture_name}): IDENTICAL"
                );
            } else {
                let min_len = stdout_lines.len().min(fixture_lines.len());
                let mut first_diff = None;
                for i in 0..min_len {
                    if stdout_lines[i] != fixture_lines[i] {
                        first_diff = Some(i);
                        break;
                    }
                }
                match first_diff {
                    Some(i) => {
                        eprintln!(
                            "x1b_it_eligibility_report: {flag} on {inst} ({fixture_name}): DIFFERS at line {i}\n  stdout:  {}\n  fixture: {}",
                            stdout_lines[i], fixture_lines[i]
                        );
                    }
                    None => {
                        eprintln!(
                            "x1b_it_eligibility_report: {flag} on {inst} ({fixture_name}): DIFFERS (length {} vs {})",
                            stdout_lines.len(),
                            fixture_lines.len()
                        );
                    }
                }
            }
        }
    }
}
