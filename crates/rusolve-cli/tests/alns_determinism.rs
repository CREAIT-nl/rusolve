//! `mip_alns` at threads=1 is deterministic per process: two separate solves give the same
//! search. The in-process unit test cannot show this while other tests run beside it,
//! because the work counters that budget ALNS are process-wide (crates/rusolve-ekk/src/work.rs);
//! hosts that run several solves at once in one process (Python,
//! MCP) therefore get a budget that depends on the other solves.

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

fn solve(mps: &Path) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(mps)
        .args(["--opt", "mip_alns=true", "--opt", "threads=1"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    // Everything except the wall-clock part of the nodes line.
    stdout
        .lines()
        .map(|l| match l.find(", ") {
            Some(i) if l.starts_with("nodes:") => l[..i].to_string(),
            _ => l.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn alns_on_is_deterministic_across_processes() {
    let Some(mps) = find_in_checkouts("benches/data/mip50/app1-1.mps") else {
        eprintln!("skip: benches/data/mip50/app1-1.mps not present");
        return;
    };
    let a = solve(&mps);
    let b = solve(&mps);
    assert!(a.contains("objective: -3"), "{a}");
    assert_eq!(a, b);
}
