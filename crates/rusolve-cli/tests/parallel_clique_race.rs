//! Regression: parallel MIP workers (threads=8, parallel=on) call
//! `HighsCliqueTable::add_implications` concurrently on the one shared clique table. Its
//! id-list scratch once lived inside the table, so concurrent calls raced on one Vec and
//! panicked with "range end index N out of range" (part2.rs). The scratch is per thread now.
//!
//! Second race (same file, `parallel_threads8_nidda_has_no_stale_pool_activities`): a heuristic's
//! local domain copied under the parallel lock dropped its pool-0 registration but kept a worker
//! domain's stale cut-pool arrays, so its first bound change indexed past them (propagation.rs
//! "index out of bounds"); and releasing such a copy gc'd the shared conflict-pool log.

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

// A debug build needs about 460 s for this solve, so it runs in release
// only; `locked_node_separator_solution_is_buffered_then_synced` (rusolve-ekk) covers the
// parallel-lock incumbent path in debug.
#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "about 460 s in a debug build, over the 120 s limit; run with --release"
)]
fn parallel_threads8_beasley_c3_reaches_optimum() {
    let Some(mps) = find_in_checkouts("benches/data/miplib/beasleyC3.mps") else {
        eprintln!("skip: benches/data/miplib/beasleyC3.mps not present");
        return;
    };
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(&mps)
        .args(["--opt", "threads=8", "--opt", "parallel=on"])
        // release solves in about 35 s
        .args(["--opt", "time_limit=120"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        out.status.success() && !stderr.contains("panicked"),
        "parallel solve failed:\n{stderr}\n{stdout}"
    );
    assert!(
        stdout.contains("status:    Optimal"),
        "not optimal:\n{stdout}"
    );
    assert!(
        stdout.contains("objective: 754"),
        "wrong objective:\n{stdout}"
    );
}

/// Before the fix this panicked in every one of 36 runs within 8 s (RINS on a worker's local
/// domain: `activitycutsinf_` shorter than the shared cut pool's row count).
#[test]
fn parallel_threads8_nidda_has_no_stale_pool_activities() {
    let Some(mps) = find_in_checkouts("benches/data/miplib/neos-3754480-nidda.mps") else {
        eprintln!("skip: benches/data/miplib/neos-3754480-nidda.mps not present");
        return;
    };
    for h16 in ["false", "true"] {
        let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
            .env("RUSOLVE_MIP_HIGHS", "1")
            .arg("mip")
            .arg(&mps)
            .args(["--opt", "threads=8", "--opt", "parallel=on"])
            .args(["--opt", &format!("h16_mip_workers={h16}")])
            .args(["--opt", "time_limit=8"])
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            out.status.success() && !stderr.contains("panicked"),
            "parallel solve failed (h16_mip_workers={h16}):\n{stderr}"
        );
    }
}
