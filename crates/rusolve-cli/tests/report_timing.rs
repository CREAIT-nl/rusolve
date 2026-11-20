//! Regression: the MIP "Solving report" timing block (Presolve / Solve / Postsolve) once
//! printed zeros because the profiling clocks were no-op stubs. HiGHS reports them from
//! `HighsMipSolver::profiling_`; they are display-only and never read by the search.

use std::path::{Path, PathBuf};
use std::process::Command;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

/// The main checkout's root (from a git worktree) or this checkout's root, whichever is
/// first found to hold `rel` (benches/data is gitignored, so worktrees do not have it).
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

/// The leading number of the first line containing `tag`.
fn timing_of(text: &str, tag: &str) -> f64 {
    let line = text
        .lines()
        .find(|l| l.contains(tag))
        .unwrap_or_else(|| panic!("no `{tag}` line in:\n{text}"));
    line.split_whitespace().next().unwrap().parse().unwrap()
}

#[test]
fn solving_report_timing_breakdown_is_nonzero_and_bounded() {
    let Some(mps) = find_in_checkouts("benches/data/mip50/app1-1.mps") else {
        eprintln!("skip: benches/data/mip50/app1-1.mps not present");
        return;
    };
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(&mps)
        .args(["--opt", "time_limit=3"])
        .args(["--opt", "output_flag=true", "--opt", "log_to_console=true"])
        .output()
        .unwrap();
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    let total = text
        .lines()
        .find(|l| l.trim_start().starts_with("Timing"))
        .and_then(|l| l.split_whitespace().nth(1))
        .unwrap_or_else(|| panic!("no Timing line in:\n{text}"))
        .parse::<f64>()
        .unwrap();
    let sum = timing_of(&text, "(Presolve)")
        + timing_of(&text, "(Solve)")
        + timing_of(&text, "(Postsolve)");
    // two-decimal rounding of each of the three terms: allow 0.015 slack
    assert!(sum > 0.0, "timing block is all zeros:\n{text}");
    assert!(
        sum <= total + 0.015,
        "parts {sum} exceed total {total}:\n{text}"
    );
}
