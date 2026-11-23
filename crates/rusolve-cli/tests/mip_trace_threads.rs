//! `rusolve-cli mip-trace <instance> threads=N parallel=on` must reproduce
//! the oracle trace (`tools/highs-oracle/mip_trace.cpp` with the same overrides)
//! line for line, as `mip_trace.rs` does at threads=1.
//!
//! MIP is bit-identical to HiGHS *per thread count*, not across thread counts:
//! `getMaxNumWorkers` is ceil(1.7 * threads) with `parallel=on` (HighsMipSolver.h
//! 161-168), and the number of workers changes the search. The traces are target-only (`target/highs-fixtures/w5/mip/`); with none
//! generated the test prints a SKIP banner (`RUSOLVE_REQUIRE_TARGET_FIXTURES=1` turns a
//! skip into a failure).
//!
//! Gates: `sp150x300d` and `pg` at threads 2, 4, 8 (`parallel=on`), and `parallel=off` at
//! threads = 4 against the committed threads = 1 trace. `air05` and `enlight_hard` are
//! written but ignored: their `w5/mip` traces come from a libstdc++ oracle and need regenerating on macOS first.

use std::path::{Path, PathBuf};
use std::process::Command;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

/// The main checkout's root (from a git worktree) or this checkout's root, whichever is
/// first found to hold `rel`.
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
fn check_trace(test: &str, instance: &str, expected: &str, args: &[&str]) {
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

/// The committed threads = 1 trace `tests/fixtures/mip/<instance>.trace` (rusolve-ekk).
fn oracle_threads_1(instance: &str) -> String {
    let path = Path::new(CRATE_DIR)
        .join("../rusolve-ekk/tests/fixtures/mip")
        .join(format!("{instance}.trace"));
    std::fs::read_to_string(path).unwrap()
}

fn check_parallel_on(test: &str, instance: &str) {
    for threads in [2, 4, 8] {
        let Some(expected) = oracle_w5(test, instance, threads) else {
            return;
        };
        let t = threads.to_string();
        let args = [format!("threads={t}"), "parallel=on".to_string()];
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        check_trace(&format!("{test} t{threads}"), instance, &expected, &args);
    }
}

#[test]
#[ignore = "slow; cargo test --release -p rusolve-cli --test mip_trace_threads -- --ignored"]
fn w5_mip_it_trace_matches_highs_per_thread_count() {
    check_parallel_on("sp150x300d", "sp150x300d");
    check_parallel_on("pg", "pg");
}

#[test]
#[ignore = "slow; cargo test --release -p rusolve-cli --test mip_trace_threads -- --ignored"]
fn w5_mip_it_parallel_off_matches_threads_1() {
    for instance in ["sp150x300d", "pg"] {
        // threads = 4, parallel = off: one worker, so the committed threads = 1 trace.
        check_trace(
            &format!("{instance} parallel=off"),
            instance,
            &oracle_threads_1(instance),
            &["threads=4", "parallel=off"],
        );
    }
}

#[test]
#[ignore = "libstdc++ oracle: regenerate w5/mip on macOS"]
fn w5_mip_it_air05_matches_highs_per_thread_count() {
    check_parallel_on("air05", "air05");
}

#[test]
#[ignore = "libstdc++ oracle: regenerate w5/mip on macOS"]
fn w5_mip_it_enlight_hard_matches_highs_per_thread_count() {
    check_parallel_on("enlight_hard", "enlight_hard");
}
