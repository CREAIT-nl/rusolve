//! `rusolve-cli mip-trace` must reproduce the HiGHS
//! oracle trace (`tools/highs-oracle/mip_trace.cpp`, threads=1, seed 0) line for line.
//! The traces are in `crates/rusolve-ekk/tests/fixtures/mip/`. The port runs
//! its own IPX for the analytic centre and the LP IPM fallback:
//! `RUSOLVE_MIP_AC_FIXTURE` is not set, and the `.ac` files there are kept only as a
//! debugging reference (`mip/ipm_replay.rs`).
//!
//! The instances live in `benches/data/miplib/`, which is gitignored: the test looks in
//! this checkout first, then in the main checkout (`git rev-parse --git-common-dir`), and
//! skips with a note when the file is in neither.

use std::path::{Path, PathBuf};
use std::process::Command;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

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

/// Diff `mip-trace <instance> <args>` against `<fixture>.trace`.
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
        "{fixture}: mip-trace failed: {}",
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

#[test]
fn l1_sp150x300d_presolve_off() {
    check_trace("sp150x300d", "sp150x300d.presolve_off", &["presolve=off"]);
}

#[test]
fn l2_sp150x300d() {
    check_trace("sp150x300d", "sp150x300d", &[]);
}

/// The analytic-centre IPX results themselves equal the oracle's `.ac` records.
#[test]
fn ipm_crosscheck_sp150x300d() {
    run_trace(
        "sp150x300d",
        "sp150x300d.presolve_off",
        &["presolve=off"],
        true,
    );
    run_trace("sp150x300d", "sp150x300d", &[], true);
}

#[test]
#[ignore = "slow; cargo test --release -p rusolve-cli --test mip_trace -- --ignored"]
fn l2_pg() {
    check_trace("pg", "pg", &[]);
}

#[test]
#[ignore = "slow; cargo test --release -p rusolve-cli --test mip_trace -- --ignored"]
fn l2_air05() {
    check_trace("air05", "air05", &[]);
}

#[test]
#[ignore = "slow; cargo test --release -p rusolve-cli --test mip_trace -- --ignored"]
fn l2_enlight_hard() {
    check_trace("enlight_hard", "enlight_hard", &[]);
}
