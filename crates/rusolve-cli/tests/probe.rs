//! `rusolve probe` and `rusolve classify --format json`. Every test runs the real binary.

use std::path::{Path, PathBuf};
use std::process::Command;

const PROBE_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-probe-1.schema.json");
const CLASSIFY_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-classify-1.schema.json");

fn write_tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-probe-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

/// A multi-row binary knapsack from a fixed LCG: the LP relaxation is fractional and the tree
/// needs thousands of nodes at 40 columns and 4 rows (a node limit of 50 stops it).
fn knapsack_mps_sized(n: usize, m: usize) -> String {
    let mut x: u64 = 12345;
    let mut r = move || {
        x = (x * 1103515245 + 12345) % (1u64 << 31);
        (x >> 8) as i64
    };
    let w: Vec<Vec<i64>> = (0..m)
        .map(|_| (0..n).map(|_| r() % 47 + 10).collect())
        .collect();
    let mut s = String::from("NAME KNAP\nROWS\n N  cost\n");
    for i in 0..m {
        s += &format!(" L  k{i}\n");
    }
    s += "COLUMNS\n    MARKER  'MARKER'  'INTORG'\n";
    for j in 0..n {
        let c = (0..m).map(|i| w[i][j]).sum::<i64>() / m as i64 + r() % 7;
        s += &format!("    x{j}  cost  -{c}");
        for (i, row) in w.iter().enumerate() {
            s += &format!("  k{i}  {}", row[j]);
        }
        s += "\n";
    }
    s += "    MARKER  'MARKER'  'INTEND'\nRHS\n";
    for (i, row) in w.iter().enumerate() {
        s += &format!("    rhs  k{i}  {}\n", row.iter().sum::<i64>() / 3);
    }
    s += "BOUNDS\n";
    for j in 0..n {
        s += &format!(" UP bnd x{j} 1\n");
    }
    s + "ENDATA\n"
}

fn knapsack_mps() -> String {
    knapsack_mps_sized(40, 4)
}

const PURE_LP: &str = "\
NAME PLP
ROWS
 N  cost
 L  r1
 L  r2
COLUMNS
    x  cost  -1  r1  1
    x  r2  1
    y  cost  -2  r1  1
    y  r2  -1
RHS
    rhs  r1  5  r2  2
ENDATA
";

fn cli(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(args)
        .env_remove("RUSOLVE_PROFILE")
        .output()
        .unwrap();
    (
        out.status.code().unwrap(),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8(out.stderr).unwrap(),
    )
}

fn json_of(args: &[&str]) -> serde_json::Value {
    let (code, out, err) = cli(args);
    assert_eq!(code, 0, "stderr: {err}");
    serde_json::from_str(out.trim()).unwrap_or_else(|e| panic!("not one JSON object ({e}): {out}"))
}

fn probe(model: &Path, extra: &[&str]) -> serde_json::Value {
    let mut args = vec!["probe", model.to_str().unwrap(), "--format", "json"];
    args.extend_from_slice(extra);
    json_of(&args)
}

fn required(node: &serde_json::Value) -> Vec<String> {
    node["required"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_str().unwrap().to_string())
        .collect()
}

fn assert_required(schema: &serde_json::Value, v: &serde_json::Value, at: &str) {
    for k in required(schema) {
        assert!(v.get(&k).is_some(), "{at}.{k} missing: {v}");
    }
}

fn feature(v: &serde_json::Value, name: &str) -> serde_json::Value {
    v["features"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["name"] == name)
        .unwrap_or_else(|| panic!("no feature {name}"))["value"]
        .clone()
}

#[test]
fn probe_json_carries_every_required_field_of_the_schema() {
    let schema: serde_json::Value = serde_json::from_str(PROBE_SCHEMA).unwrap();
    let m = write_tmp("knap.mps", &knapsack_mps());
    let v = probe(&m, &[]);
    assert_eq!(v["schema_version"], "rusolve.probe/1");
    assert_required(&schema, &v, "");
    let p = &schema["properties"];
    for (sec, ps) in [
        ("input", &p["input"]),
        ("budget", &p["budget"]),
        ("presolve", &p["presolve"]),
        ("root_lp", &p["root_lp"]),
        ("root", &p["root"]),
        ("tree", &p["tree"]),
        ("bounds", &p["bounds"]),
        ("incumbent", &p["incumbent"]),
    ] {
        assert_required(ps, &v[sec], sec);
    }
    for f in v["features"].as_array().unwrap() {
        assert_required(&p["features"]["items"], f, "features[]");
    }
    for u in v["unavailable"].as_array().unwrap() {
        assert_required(&p["unavailable"]["items"], u, "unavailable[]");
    }
    assert!(v["features_version"].as_u64().unwrap() >= 1);
}

#[test]
fn probe_reports_root_and_tree_features_on_a_knapsack() {
    let m = write_tmp("knap2.mps", &knapsack_mps());
    let v = probe(&m, &[]);
    assert_eq!(v["is_mip"], true);
    assert_eq!(v["input"]["integer_columns"], 40);
    // The relaxation of a knapsack with 3 rows is fractional.
    assert!(feature(&v, "root_lp.fractional_frac").as_f64().unwrap() > 0.0);
    assert!(v["root_lp"]["iterations"].as_i64().unwrap() > 0);
    // The tree explores nodes after the root and reports per-node work.
    assert!(v["tree"]["nodes"].as_i64().unwrap() > 1, "{v}");
    assert!(feature(&v, "tree.lp_iterations_per_node").as_f64().unwrap() > 0.0);
    assert_eq!(v["incumbent"]["found"], true);
    // The root bound is no worse than the root LP (minimisation).
    let lp = v["bounds"]["root_lp_objective"].as_f64().unwrap();
    let rb = v["bounds"]["root_bound_after_cuts"].as_f64().unwrap();
    assert!(rb >= lp - 1e-6, "root bound {rb} below root LP {lp}");
    assert!(
        v["bounds"]["dual_bound"].as_f64().unwrap() >= rb - 1e-6,
        "final dual bound below root bound: {v}"
    );
}

#[test]
fn probe_node_limit_stops_the_tree() {
    let m = write_tmp("knap3.mps", &knapsack_mps());
    let v = probe(&m, &["--node-limit", "3"]);
    assert_eq!(v["budget"]["node_limit"], 3);
    assert!(v["tree"]["nodes"].as_i64().unwrap() <= 4, "{v}");
    assert_eq!(v["status"], "SolutionLimit");
    // Open nodes are left, so the depth feature is present.
    assert!(
        v["tree"]["open_nodes_max_depth"].as_i64().unwrap() >= 1,
        "{v}"
    );
}

#[test]
fn probe_is_deterministic_apart_from_times() {
    let m = write_tmp("knap4.mps", &knapsack_mps());
    let a = probe(&m, &[]);
    let b = probe(&m, &[]);
    for k in ["presolve", "bounds", "incumbent", "input"] {
        assert_eq!(a[k], b[k], "{k}");
    }
    for name in [
        "root_lp.iterations",
        "root_lp.fractional_frac",
        "root.gap_closed_by_cuts",
        "root.cuts_in_pool",
        "tree.lp_iterations_per_node",
        "tree.open_nodes_max_depth",
    ] {
        assert_eq!(feature(&a, name), feature(&b, name), "{name}");
    }
}

#[test]
fn probe_counters_equal_the_default_solve() {
    // The probe adds nothing to the search: its node-limited run ends on the same node count and
    // LP iteration total as `rusolve mip` on a model the limit does not stop.
    let m = write_tmp("knap5.mps", &knapsack_mps_sized(24, 3));
    let v = probe(&m, &["--node-limit", "100000"]);
    assert_eq!(v["status"], "Optimal");
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["mip", m.to_str().unwrap(), "--csv", "--seed", "0"])
        .env("RUSOLVE_MIP_HIGHS", "1")
        .env_remove("RUSOLVE_PROFILE")
        .output()
        .unwrap();
    let text = String::from_utf8(out.stdout).unwrap();
    let f: Vec<&str> = text.lines().last().unwrap().split(',').collect();
    assert_eq!(f[3], v["tree"]["nodes"].to_string(), "nodes: {text}");
    assert_eq!(f[4], v["tree"]["lp_iterations_total"].to_string(), "{text}");
}

#[test]
fn probe_pure_lp_has_no_mip_sections() {
    let m = write_tmp("plp.mps", PURE_LP);
    let v = probe(&m, &[]);
    assert_eq!(v["is_mip"], false);
    for k in ["presolve", "root", "tree", "bounds"] {
        assert!(v[k].is_null(), "{k}: {v}");
    }
    assert_eq!(v["root_lp"]["status"], "Optimal");
    assert!(v["root_lp"]["fractional_fraction"].is_null());
    assert!(feature(&v, "tree.nodes_per_second").is_null());
}

#[test]
fn probe_unrecorded_features_are_null_and_listed() {
    let m = write_tmp("knap6.mps", &knapsack_mps());
    let v = probe(&m, &[]);
    assert!(v["incumbent"]["first_time_s"].is_null());
    assert!(v["incumbent"]["first_source"].is_null());
    let un: Vec<&str> = v["unavailable"]
        .as_array()
        .unwrap()
        .iter()
        .map(|u| u["feature"].as_str().unwrap())
        .collect();
    for k in [
        "incumbent.first_time_s",
        "incumbent.first_source",
        "tree.depth_max",
    ] {
        assert!(un.contains(&k), "{k} not listed: {un:?}");
    }
}

#[test]
fn probe_text_output_lists_the_features() {
    let m = write_tmp("knap7.mps", &knapsack_mps());
    let (code, out, _) = cli(&["probe", m.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(out.contains("root_lp.fractional_frac"), "{out}");
    assert!(out.contains("features (version 1)"), "{out}");
}

#[test]
fn probe_exit_codes_follow_a5() {
    let m = write_tmp("knap8.mps", &knapsack_mps());
    assert_eq!(cli(&["probe"]).0, 2);
    assert_eq!(
        cli(&["probe", m.to_str().unwrap(), "--node-limit", "0"]).0,
        2
    );
    assert_eq!(cli(&["probe", m.to_str().unwrap(), "--bogus"]).0, 2);
    assert_eq!(cli(&["probe", "/nonexistent/model.mps"]).0, 3);
}

#[test]
fn classify_json_carries_every_required_field_of_the_schema() {
    let schema: serde_json::Value = serde_json::from_str(CLASSIFY_SCHEMA).unwrap();
    let m = write_tmp("c1.mps", &knapsack_mps());
    let v = json_of(&["classify", m.to_str().unwrap(), "--format", "json"]);
    assert_eq!(v["schema_version"], "rusolve.classify/1");
    assert_required(&schema, &v, "");
    assert_required(&schema["properties"]["input"], &v["input"], "input");
    assert_required(&schema["properties"]["options"], &v["options"], "options");
    for f in v["features"].as_array().unwrap() {
        assert_required(&schema["properties"]["features"]["items"], f, "features[]");
        assert!(!f["description"].as_str().unwrap().is_empty());
    }
}

#[test]
fn classify_json_lists_the_44_static_features_with_their_values() {
    let m = write_tmp("c2.mps", &knapsack_mps());
    let v = json_of(&["classify", m.to_str().unwrap(), "--format", "json"]);
    assert_eq!(v["features"].as_array().unwrap().len(), 44);
    assert_eq!(v["features_version"], 1);
    // The same values as the pre-existing `--json` object.
    let old = json_of(&["classify", m.to_str().unwrap(), "--json"]);
    assert_eq!(v["detail"], old);
    for f in v["features"].as_array().unwrap() {
        assert_eq!(
            f["value"].as_f64(),
            old["features"][f["name"].as_str().unwrap()].as_f64(),
            "{f}"
        );
    }
    // A 40-column binary model with knapsack rows.
    assert!(feature(&v, "vars.binary_frac").as_f64().unwrap() > 0.99);
    assert!(feature(&v, "rows.knapsack").as_f64().unwrap() > 0.0);
}

#[test]
fn classify_json_text_output_is_unchanged_by_the_format_flag() {
    let m = write_tmp("c3.mps", &knapsack_mps());
    let plain = cli(&["classify", m.to_str().unwrap()]);
    let text = cli(&["classify", m.to_str().unwrap(), "--format", "text"]);
    assert_eq!(plain.0, 0);
    assert_eq!(plain, text);
    assert!(plain.1.contains("rows"), "{}", plain.1);
}

#[test]
fn classify_json_rejects_a_bad_format_and_a_missing_model() {
    let m = write_tmp("c4.mps", &knapsack_mps());
    assert_eq!(
        cli(&["classify", m.to_str().unwrap(), "--format", "xml"]).0,
        2
    );
    assert_eq!(
        cli(&["classify", "/nonexistent/model.mps", "--format", "json"]).0,
        3
    );
}
