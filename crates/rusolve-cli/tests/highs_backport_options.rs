//! Integration tests: PR tests TestMipSolver issue-3170 (#3177, merge 2100494d8),
//! issue-3171 (#3179, merge 28339ce4b), issue-3461 (#3468, merge 3af3ff719).
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! (2) each option on alone and all 15 together on pg,
//! sp150x300d, air05 and enlight_hard (optimum and feasibility on the original model);
//! (3) the issue instances with the fix flag on (asserted) and off (printed); (4) the
//! default-on eligibility report (flag-on `mip-trace` byte-compared to the 1.15.1
//! fixtures, printed, never a failure).

use std::path::{Path, PathBuf};
use std::process::Command;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

const ALL_OPTIONS: &[&str] = &[
    "mip_sb_budget_pre_restart",
    "mip_dfs_node_limit",
    "mip_reduced_stack_fast",
    "mip_cdouble_act_contrib",
    "mip_keep_branch_reason",
    "mip_redundant_vbd",
    "mip_path_same_bds",
    "mip_strengthen_small_contrib",
    "mip_lurking_bounds_robust",
    "mip_dual_subst_limits",
    "mip_ranged_dual_implied",
    "mip_machine_sched_sepa",
    "mip_dual_fix_probing",
    "mip_aggregate_implications",
    "mip_clique_normalise",
    "mip_coef_change_impl_bounds",
    "mip_dual_side_lost",
];

fn miplib_instance(name: &str) -> Option<PathBuf> {
    let rel = format!("benches/data/miplib/{name}.mps");
    let local = Path::new(CRATE_DIR).join("../..").join(&rel);
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
    let main = common.parent()?.join(&rel);
    main.exists().then_some(main)
}

fn fixture_instance(name: &str) -> PathBuf {
    Path::new(CRATE_DIR)
        .join("../rusolve-ekk/tests/fixtures/x1")
        .join(format!("{name}.mps"))
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

/// A flag-off run: printed, never a failure (the 1.15.1 wrong value
/// documents the bug).
fn report_flag_off(label: &str, path: &Path, extra_args: &[&str]) {
    let (_ok, _stderr, (status, obj, ..)) = run_mip_csv_unchecked(path, extra_args);
    eprintln!("{label}: status={status} obj={obj}");
}

/// TestMipSolver.cpp `objectiveOk`: relative error to max(1, |required|) below 1e-10.
fn objective_ok(objective: f64, required: f64) -> bool {
    (objective - required).abs() / required.abs().max(1.0) < 1e-10
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

fn check_option_on_instance(inst_name: &str, opt: &str, expected_obj: f64) {
    let Some(mps) = miplib_instance(inst_name) else {
        eprintln!("skipping {inst_name}: benches/data/miplib/{inst_name}.mps not found");
        return;
    };

    // mip-trace check
    let obj = run_mip_trace(&mps, "3600", &[&format!("{opt}=true")]);
    assert!(
        rel_close(obj, expected_obj, 1e-6),
        "{inst_name} with {opt}=true: obj={obj}, expected={expected_obj}"
    );

    // mip --csv for feasibility
    let (_status, obj_csv, _bound, _nodes, _iters, _time, x) =
        run_mip_csv(&mps, &["--opt", &format!("{opt}=true")]);
    assert!(
        rel_close(obj_csv, expected_obj, 1e-6),
        "{inst_name} with {opt}=true (mip --csv): obj={obj_csv}, expected={expected_obj}"
    );
    check_feasibility(&mps, &x);
}

fn check_all_options_on_instance(inst_name: &str, expected_obj: f64) {
    let Some(mps) = miplib_instance(inst_name) else {
        eprintln!("skipping {inst_name}: benches/data/miplib/{inst_name}.mps not found");
        return;
    };

    let mut args: Vec<String> = Vec::new();
    for opt in ALL_OPTIONS {
        args.push(format!("{opt}=true"));
    }
    let opt_args: Vec<&str> = args.iter().map(|s| s.as_str()).collect();

    // mip-trace check
    let obj = run_mip_trace(&mps, "3600", &opt_args);
    assert!(
        rel_close(obj, expected_obj, 1e-6),
        "{inst_name} with all options on: obj={obj}, expected={expected_obj}"
    );

    // mip --csv for feasibility
    let mut csv_args: Vec<String> = Vec::new();
    for opt in ALL_OPTIONS {
        csv_args.push("--opt".to_string());
        csv_args.push(format!("{opt}=true"));
    }
    let csv_opt_args: Vec<&str> = csv_args.iter().map(|s| s.as_str()).collect();
    let (_status, obj_csv, _bound, _nodes, _iters, _time, x) = run_mip_csv(&mps, &csv_opt_args);
    assert!(
        rel_close(obj_csv, expected_obj, 1e-6),
        "{inst_name} with all options on (mip --csv): obj={obj_csv}, expected={expected_obj}"
    );
    check_feasibility(&mps, &x);
}

/// TestMipSolver `issue-3170` (#3177): `mip_redundant_vbd`. The test runs both
/// presolve settings with gaps 0/0 (the C++ `solve()` calls `resetOptions()`, so its
/// second, presolve-on run has default gaps).
#[test]
fn x1_it_issue_3170_1() {
    let mps = fixture_instance("issue-3170-1");
    assert!(
        mps.exists(),
        "fixture issue-3170-1.mps not found at {mps:?}"
    );
    let expected = -21446579.3647;
    for presolve in ["off", "on"] {
        let presolve_opt = format!("presolve={presolve}");
        let args = |flag: &'static str| {
            [
                "--mip-rel-gap",
                "0",
                "--mip-abs-gap",
                "0",
                "--opt",
                presolve_opt.as_str(),
                "--opt",
                flag,
            ]
        };
        let (_status, obj, ..) = run_mip_csv(&mps, &args("mip_redundant_vbd=true"));
        assert!(
            objective_ok(obj, expected),
            "issue-3170-1 flag=true presolve={presolve}: obj={obj}, expected={expected}"
        );
        report_flag_off(
            &format!("issue-3170-1 flag=false presolve={presolve}"),
            &mps,
            &args("mip_redundant_vbd=false"),
        );
    }
}

/// TestMipSolver `issue-3171` (#3179): `mip_path_same_bds`, default options.
#[test]
fn x1_it_issue_3171() {
    let mps = fixture_instance("issue-3171");
    assert!(mps.exists(), "fixture issue-3171.mps not found at {mps:?}");
    let expected = 42215.5250005;
    let (_status, obj, ..) = run_mip_csv(&mps, &["--opt", "mip_path_same_bds=true"]);
    assert!(
        objective_ok(obj, expected),
        "issue-3171 flag=true: obj={obj}, expected={expected}"
    );
    report_flag_off(
        "issue-3171 flag=false",
        &mps,
        &["--opt", "mip_path_same_bds=false"],
    );
}

/// TestMipSolver `issue-3461` (#3468): `mip_strengthen_small_contrib`, gaps 0/0.
/// Flag off gives the wrong 23138 the test comment names.
#[test]
fn x1_it_issue_3461() {
    let mps = fixture_instance("issue-3461");
    assert!(mps.exists(), "fixture issue-3461.mps not found at {mps:?}");
    let expected = 21166.0;
    let args = |flag: &'static str| ["--mip-rel-gap", "0", "--mip-abs-gap", "0", "--opt", flag];
    let (_status, obj, ..) = run_mip_csv(&mps, &args("mip_strengthen_small_contrib=true"));
    assert!(
        objective_ok(obj, expected),
        "issue-3461 flag=true: obj={obj}, expected={expected}"
    );
    report_flag_off(
        "issue-3461 flag=false",
        &mps,
        &args("mip_strengthen_small_contrib=false"),
    );
}

#[test]
#[ignore = "slow"]
fn x1_it_each_option_on_pg() {
    for opt in ALL_OPTIONS {
        check_option_on_instance("pg", opt, -8674.34260712);
    }
}

#[test]
#[ignore = "slow"]
fn x1_it_each_option_on_sp150x300d() {
    for opt in ALL_OPTIONS {
        check_option_on_instance("sp150x300d", opt, 69.0);
    }
}

#[test]
#[ignore = "slow"]
fn x1_it_each_option_on_air05() {
    for opt in ALL_OPTIONS {
        check_option_on_instance("air05", opt, 26374.0);
    }
}

#[test]
#[ignore = "slow"]
fn x1_it_each_option_on_enlight_hard() {
    for opt in ALL_OPTIONS {
        check_option_on_instance("enlight_hard", opt, 37.0);
    }
}

#[test]
#[ignore = "slow"]
fn x1_it_all_options_on_pg() {
    check_all_options_on_instance("pg", -8674.34260712);
}

#[test]
#[ignore = "slow"]
fn x1_it_all_options_on_sp150x300d() {
    check_all_options_on_instance("sp150x300d", 69.0);
}

#[test]
#[ignore = "slow"]
fn x1_it_all_options_on_air05() {
    check_all_options_on_instance("air05", 26374.0);
}

#[test]
#[ignore = "slow"]
fn x1_it_all_options_on_enlight_hard() {
    check_all_options_on_instance("enlight_hard", 37.0);
}

/// The flags of the default-on eligibility table: the wrong-answer
/// fixes and #3336b.
const ELIGIBILITY_FLAGS: &[(&str, &str)] = &[
    ("#3184", "mip_cdouble_act_contrib"),
    ("#3357", "mip_keep_branch_reason"),
    ("#3177", "mip_redundant_vbd"),
    ("#3179", "mip_path_same_bds"),
    ("#3468", "mip_strengthen_small_contrib"),
    ("#3474/#3481", "mip_lurking_bounds_robust"),
    ("#3336b", "mip_reduced_stack_fast"),
];

/// A report and never a failure: whether `mip-trace` with the flag on
/// is byte-identical to the 1.15.1 fixture (`mip_trace.rs` runs the same instances and
/// arguments with every flag off).
#[test]
#[ignore = "slow"]
fn x1_it_default_on_eligibility_report() {
    let fixtures = Path::new(CRATE_DIR).join("../rusolve-ekk/tests/fixtures/mip");
    let cases: &[(&str, &str, &[&str])] = &[
        ("pg", "pg", &[]),
        ("sp150x300d", "sp150x300d", &[]),
        ("sp150x300d", "sp150x300d.presolve_off", &["presolve=off"]),
    ];
    let mut table = Vec::new();
    for &(pr, flag) in ELIGIBILITY_FLAGS {
        for &(instance, fixture, args) in cases {
            let Some(mps) = miplib_instance(instance) else {
                eprintln!("skipping {fixture}: benches/data/miplib/{instance}.mps not found");
                continue;
            };
            let expected = std::fs::read(fixtures.join(format!("{fixture}.trace"))).unwrap();
            let on = format!("{flag}=true");
            let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
                .arg("mip-trace")
                .arg(&mps)
                .args(args)
                .arg(&on)
                .env_remove("RUSOLVE_MIP_AC_FIXTURE")
                .output()
                .unwrap();
            assert!(
                out.status.success(),
                "{fixture} {on}: mip-trace failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            let verdict = if out.stdout == expected {
                "IDENTICAL".to_string()
            } else {
                let e = String::from_utf8_lossy(&expected);
                let a = String::from_utf8_lossy(&out.stdout);
                match e.lines().zip(a.lines()).position(|(x, y)| x != y) {
                    Some(i) => format!("differs at line {}", i + 1),
                    None => "differs in length".to_string(),
                }
            };
            table.push(format!("{pr:<12} {flag:<34} {fixture:<24} {verdict}"));
        }
    }
    eprintln!("X1-IT default-on eligibility (flag on vs 1.15.1 mip-trace fixture):");
    for row in &table {
        eprintln!("  {row}");
    }
}
