//! `rusolve mip --wasm-hook` on the port MIP path (`RUSOLVE_MIP_HIGHS=1`). A component
//! can only reorder the search, so the optimum must not move; `stats.wasm` reports the
//! counters; a bad component is an input error (exit 3) and the budgets without `--wasm-hook`
//! are a usage error (exit 2). Every test runs the real binary.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

/// The small knapsack MIP of the ekk fixtures (root pool separation closes it; the
/// `most_fractional` component still gets its `cut_score` calls on it).
fn fixture() -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../crates/rusolve-ekk/tests/fixtures/s5/e2e_knap.mps")
        .to_string_lossy()
        .into_owned()
}

fn component(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../crates/rusolve-wasm/components")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn lcg(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

/// The 28-item, 4-row knapsack of the `rusolve-wasm` hook tests as MPS text (minimise
/// `-value x`): with presolve and heuristics off it builds a real tree, which is what a
/// `select_node` component needs. The same generator as `wasm_hooks.rs` (same seed).
fn knapsack_mps() -> String {
    const N: usize = 28;
    const M: usize = 4;
    let mut s = 12345u64;
    let values: Vec<f64> = (0..N).map(|_| (20 + lcg(&mut s) % 60) as f64).collect();
    let mut weights = vec![vec![0.0; N]; M];
    for row in weights.iter_mut() {
        for w in row.iter_mut() {
            *w = (10 + lcg(&mut s) % 50) as f64;
        }
    }
    let caps: Vec<f64> = weights
        .iter()
        .map(|r| (r.iter().sum::<f64>() * 0.35).floor())
        .collect();
    let mut t = String::from("NAME  WASMK\nROWS\n N  OBJ\n");
    for i in 0..M {
        t.push_str(&format!(" L  k{i}\n"));
    }
    t.push_str("COLUMNS\n    MARKER  'MARKER'  'INTORG'\n");
    for j in 0..N {
        t.push_str(&format!("    x{j}  OBJ  -{}\n", values[j]));
        for (i, row) in weights.iter().enumerate() {
            t.push_str(&format!("    x{j}  k{i}  {}\n", row[j]));
        }
    }
    t.push_str("    MARKER  'MARKER'  'INTEND'\nRHS\n    rhs");
    for (i, cap) in caps.iter().enumerate() {
        t.push_str(&format!("  k{i}  {cap}"));
    }
    t.push_str("\nBOUNDS\n");
    for j in 0..N {
        t.push_str(&format!(" UP bnd x{j} 1\n"));
    }
    t.push_str("ENDATA\n");
    t
}

/// Writes `text` once per name (tests run in parallel and share file names).
fn write_tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-wasm-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    if std::fs::read_to_string(&p).ok().as_deref() != Some(text) {
        let staged = dir.join(format!(".{name}.{}", std::process::id()));
        std::fs::write(&staged, text).unwrap();
        std::fs::rename(&staged, &p).unwrap();
    }
    p
}

/// `RUSOLVE_MIP_HIGHS=1 rusolve mip <args>`: `(exit, stdout, stderr)`.
fn mip(args: &[&str]) -> (i32, String, String) {
    mip_env(&[], args)
}

/// The same with extra environment (e.g. `RUSOLVE_PROFILE`).
fn mip_env(env: &[(&str, &str)], args: &[&str]) -> (i32, String, String) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    cmd.env("RUSOLVE_MIP_HIGHS", "1").arg("mip").args(args);
    for (k, v) in env {
        cmd.env(k, v);
    }
    let out = cmd.output().unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

fn json(stdout: &str) -> Value {
    serde_json::from_str(stdout).unwrap_or_else(|e| panic!("not one JSON object: {e}\n{stdout}"))
}

#[test]
fn wasm_cli_hook_reaches_the_same_optimum_and_reports_stats() {
    let model = fixture();
    let (code, stdout, _) = mip(&[&model, "--no-presolve", "--format", "json"]);
    assert_eq!(code, 0);
    let plain = json(&stdout);
    assert_eq!(plain["status"], "Optimal");
    assert!(
        plain["stats"].get("wasm").is_none(),
        "no hook: no stats.wasm"
    );

    let hook = component("most_fractional.wat");
    let (code, stdout, stderr) = mip(&[
        &model,
        "--no-presolve",
        "--wasm-hook",
        &hook,
        "--format",
        "json",
    ]);
    assert_eq!(code, 0, "stderr: {stderr}");
    let env = json(&stdout);
    assert_eq!(env["status"], "Optimal");
    let (obj, plain_obj) = (
        env["objective"].as_f64().unwrap(),
        plain["objective"].as_f64().unwrap(),
    );
    assert!(
        (obj - plain_obj).abs() < 1e-6,
        "a component can only reorder, never change the answer: {obj} vs {plain_obj}"
    );
    let w = &env["stats"]["wasm"];
    assert_eq!(w["component"], "most_fractional");
    assert!(w["calls"].as_u64().unwrap() > 0, "{w}");
    assert_eq!(
        w["ok"], w["calls"],
        "the component answered every call: {w}"
    );
    assert_eq!(w["traps"].as_u64().unwrap(), 0);
    assert_eq!(w["out_of_fuel"].as_u64().unwrap(), 0);
    assert_eq!(w["invalid"].as_u64().unwrap(), 0);
    assert!(w["fuel_used"].as_u64().unwrap() > 0);
}

/// The `select_node` component through the CLI: on the 1.15.1 tree loop (where HiGHS pops
/// `popBestNode`, the hook site) it dives, and the optimum is unchanged. The budgets
/// `--wasm-fuel` and `--wasm-memory` ride along.
#[test]
fn wasm_cli_select_node_component_reaches_the_same_optimum() {
    let model = write_tmp("wasmk.mps", &knapsack_mps());
    let model = model.to_str().unwrap();
    let off = [
        "--no-presolve",
        "--opt",
        "mip_heuristic_effort=0",
        "--opt",
        "mip_heuristic_run_feasibility_jump=false",
    ];
    let profile = [("RUSOLVE_PROFILE", "highs-1.15.1")];
    let mut plain_args: Vec<&str> = vec![model, "--format", "json"];
    plain_args.extend_from_slice(&off);
    let (code, stdout, stderr) = mip_env(&profile, &plain_args);
    assert_eq!(code, 0, "stderr: {stderr}");
    let plain = json(&stdout);
    assert_eq!(plain["status"], "Optimal");

    let hook = component("deepest_first.wat");
    let mut args: Vec<&str> = vec![
        model,
        "--wasm-hook",
        &hook,
        "--wasm-fuel",
        "100000",
        "--wasm-memory",
        "1048576",
        "--format",
        "json",
    ];
    args.extend_from_slice(&off);
    let (code, stdout, stderr) = mip_env(&profile, &args);
    assert_eq!(code, 0, "stderr: {stderr}");
    let env = json(&stdout);
    assert_eq!(env["status"], "Optimal");
    let (obj, plain_obj) = (
        env["objective"].as_f64().unwrap(),
        plain["objective"].as_f64().unwrap(),
    );
    assert!(
        (obj - plain_obj).abs() < 1e-6,
        "{obj} vs {plain_obj}: the component only reorders the tree"
    );
    let w = &env["stats"]["wasm"];
    assert_eq!(w["component"], "deepest_first");
    assert!(w["calls"].as_u64().unwrap() > 0, "{w}");
    assert_eq!(w["ok"], w["calls"], "{w}");
}

#[test]
fn wasm_cli_importing_component_is_an_input_error_naming_the_file() {
    let bad = write_tmp(
        "bad_import.wat",
        r#"(module (import "env" "now" (func (result i64))))"#,
    );
    let model = fixture();
    let (code, _stdout, stderr) = mip(&[
        &model,
        "--wasm-hook",
        bad.to_str().unwrap(),
        "--format",
        "json",
    ]);
    assert_eq!(code, 3, "stderr: {stderr}");
    assert!(
        stderr.contains("bad_import.wat") && stderr.contains("no host access"),
        "stderr: {stderr}"
    );
}

#[test]
fn wasm_cli_budgets_without_a_hook_are_a_usage_error() {
    let model = fixture();
    let (code, _stdout, stderr) = mip(&[&model, "--wasm-fuel", "10", "--format", "json"]);
    assert_eq!(code, 2, "stderr: {stderr}");
    assert!(stderr.contains("--wasm-hook"), "stderr: {stderr}");
    let (code, _stdout, stderr) = mip(&[&model, "--wasm-memory", "4096"]);
    assert_eq!(code, 2, "stderr: {stderr}");
    assert!(stderr.contains("--wasm-hook"), "stderr: {stderr}");
}
