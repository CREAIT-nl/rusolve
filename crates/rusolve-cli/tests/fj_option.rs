//! `mip_heuristic_run_feasibility_jump` (HiGHS 1.15.1, HighsOptions.h:1229) is
//! settable by name through `rusolve-cli mip-trace key=value`. The heuristics self-check and
//! the oracle comparisons run with FJ off on both sides.

use std::path::{Path, PathBuf};
use std::process::Command;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

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

/// The `final ...` line of `rusolve-cli mip-trace <path> <args>`.
fn mip_trace_final(path: &Path, args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("mip-trace")
        .arg(path)
        .args(args)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "mip-trace failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    stdout
        .lines()
        .find(|l| l.starts_with("final "))
        .unwrap_or_else(|| panic!("no final line in:\n{stdout}"))
        .to_string()
}

#[test]
fn s5_data_mip_trace_accepts_fj_option() {
    let knap = Path::new(CRATE_DIR).join("../rusolve-ekk/tests/fixtures/s5/e2e_knap.mps");
    let fin = mip_trace_final(&knap, &["mip_heuristic_run_feasibility_jump=false"]);
    // selfcheck.txt: e2e_knap off oracle run
    assert_eq!(
        fin,
        "final status=Optimal obj=-417.00000000000074 bound=-417.00000000000074 nodes=1 lp_iters=24"
    );
}

#[test]
fn s5_data_mip_trace_rejects_bad_fj_value() {
    let knap = Path::new(CRATE_DIR).join("../rusolve-ekk/tests/fixtures/s5/e2e_knap.mps");
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("mip-trace")
        .arg(&knap)
        .arg("mip_heuristic_run_feasibility_jump=maybe")
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1));
}

/// pg with FJ off equals the stock and patched HiGHS FJ-off trace (selfcheck.txt); with
/// FJ on (the default) the LP iteration count differs, so the option takes effect.
#[test]
fn s5_data_pg_fj_off_matches_oracle() {
    let Some(pg) = find_in_checkouts("benches/data/miplib/pg.mps") else {
        eprintln!("skipping: benches/data/miplib/pg.mps not found");
        return;
    };
    let base = ["threads=1", "log_dev_level=3"];
    let off = mip_trace_final(
        &pg,
        &[base[0], base[1], "mip_heuristic_run_feasibility_jump=false"],
    );
    assert_eq!(
        off,
        "final status=Optimal obj=-8674.342607117027 bound=-8675.1721976191238 nodes=27 lp_iters=34191"
    );
    let on = mip_trace_final(&pg, &base);
    assert_ne!(on, off);
}
