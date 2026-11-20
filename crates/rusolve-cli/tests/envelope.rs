//! The result envelope, exit codes, the solution file and model names on the port MIP path (`RUSOLVE_MIP_HIGHS=1 rusolve mip`):
//! the JSON result envelope, the named solution files, and the exit codes.
//!
//! Every test runs the real binary, so the envelope is checked as an agent reads
//! it - on the process's stdout and stderr - through the same entry point the
//! CLI uses (`mip_highs::run_mip_highs`).

use std::path::PathBuf;
use std::process::Command;

/// `docs/schema/rusolve-result-1.schema.json`, the normative envelope schema.
const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-result-1.schema.json");

/// min -3x - 4y + 5z + w with x + 2y + z <= 8 and z + 3w <= 10, x in [0,4],
/// y in [0,5], z in [0,3] integer and w continuous: the optimum is x=4, y=2,
/// z=0, w=0, objective -20. Four columns and two rows, all named, so every name
/// has to come back out.
const TINY_MIP: &str = "\
NAME ENV1
ROWS
 N  cost
 L  cap1
 L  cap2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  cost  -3  cap1  1
    y  cost  -4  cap1  2
    z  cost   5  cap1  1  cap2  1
    MARKER  'MARKER'  'INTEND'
    w  cost   1  cap2  3
RHS
    rhs  cap1  8  cap2  10
BOUNDS
 UP bnd x 4
 UP bnd y 5
 UP bnd z 3
ENDATA
";

fn tmp_dir() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-envelope-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_tmp(name: &str, text: &str) -> PathBuf {
    // Tests run in parallel threads and several rewrite the same file while a solver started by
    // another thread may be reading it; `fs::write` truncates first, so that reader could see an
    // empty or partial model. Write a unique file and rename it into place (atomic).
    static SEQ: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let p = tmp_dir().join(name);
    let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let staged = tmp_dir().join(format!(".{name}.{}.{n}", std::process::id()));
    std::fs::write(&staged, text).unwrap();
    std::fs::rename(&staged, &p).unwrap();
    p
}

/// Runs `rusolve mip <model> <extra>` on the port: `(exit, stdout, stderr)`.
fn mip(model: &PathBuf, extra: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(model)
        .args(extra)
        .output()
        .unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// `rusolve mip <TINY_MIP> --format json <extra>` as a parsed envelope.
fn envelope(extra: &[&str]) -> serde_json::Value {
    let model = write_tmp("tiny.mps", TINY_MIP);
    let (code, stdout, stderr) = mip(&model, &[&["--format", "json"][..], extra].concat());
    assert_eq!(code, 0, "exit code, stderr: {stderr}");
    serde_json::from_str(&stdout)
        .unwrap_or_else(|e| panic!("stdout is not one JSON object: {e}\n{stdout}"))
}

/// The `required` list of the schema and of each object it declares, so the
/// checks below test what the committed file demands rather than a list written
/// next to the test.
fn schema_required() -> Vec<(String, Vec<String>)> {
    let schema: serde_json::Value = serde_json::from_str(SCHEMA).unwrap();
    let required_of = |node: &serde_json::Value, name: &str| {
        node["required"]
            .as_array()
            .unwrap_or_else(|| panic!("{name} has no required list in the schema"))
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect()
    };
    let mut out = vec![("top".to_string(), required_of(&schema, "top"))];
    for name in ["solver", "input", "stats", "run", "solution"] {
        out.push((
            name.to_string(),
            required_of(&schema["properties"][name], name),
        ));
    }
    out
}

/// `--format json` writes one object carrying every field the committed
/// schema requires, with the values of the model that was solved, and the
/// incumbent by name when asked for.
#[test]
fn envelope_json_carries_every_required_field_of_the_schema() {
    let model = write_tmp("tiny.mps", TINY_MIP);
    let env = envelope(&["--solution"]);
    for (name, required) in schema_required() {
        let node = if name == "top" { &env } else { &env[&name] };
        assert!(
            node.is_object(),
            "{name} is missing from the envelope: {env}"
        );
        for field in &required {
            assert!(
                node.get(field).is_some(),
                "{name}.{field} is required by the schema and missing from {node}"
            );
        }
    }
    assert_eq!(env["schema_version"], "rusolve.result/1");
    assert_eq!(env["solver"]["name"], "rusolve");
    assert_eq!(env["solver"]["engine"], "highs-port");
    assert_eq!(
        env["solver"]["version"],
        env!("CARGO_PKG_VERSION"),
        "solver.version is the crate version"
    );
    assert_eq!(env["status"], "Optimal");
    assert_eq!(env["objective"], -20.0);
    assert_eq!(env["dual_bound"], -20.0);
    assert_eq!(env["gap"], 0.0);
    assert_eq!(env["has_solution"], true);
    assert_eq!(env["input"]["path"], model.to_str().unwrap());
    assert_eq!(env["input"]["rows"], 2, "the N row is not a constraint row");
    assert_eq!(env["input"]["cols"], 4);
    assert_eq!(env["input"]["integer_cols"], 3, "x, y and z are integer");
    assert_eq!(
        env["input"]["nonzeros"], 5,
        "the structural nonzeros: the objective row is not in the matrix"
    );
    assert_eq!(env["run"]["threads"], 1);
    assert_eq!(env["run"]["seed"], 0);
    assert_eq!(env["run"]["time_limit"], serde_json::Value::Null);
    assert_eq!(env["run"]["node_limit"], serde_json::Value::Null);
    assert_eq!(env["run"]["options"], serde_json::json!({}));
    assert!(env["stats"]["time_s"].as_f64().unwrap() >= 0.0);
    assert!(env["stats"]["nodes"].as_i64().unwrap() >= 1);

    // The incumbent is keyed by the model's own names, in original column
    // order (postsolve returns the solution in original space).
    let named: Vec<(&str, f64)> = env["solution"]["columns"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| (c["name"].as_str().unwrap(), c["value"].as_f64().unwrap()))
        .collect();
    assert_eq!(
        named,
        vec![("x", 4.0), ("y", 2.0), ("z", 0.0), ("w", 0.0)],
        "the solution must come back by name, in column order"
    );
}

/// The envelope stays small - no solution unless asked for - and echoes the
/// `--opt` settings as given.
#[test]
fn envelope_omits_the_solution_unless_asked_and_echoes_options() {
    let plain = envelope(&[]);
    assert!(
        plain.get("solution").is_none(),
        "the default envelope must not carry a solution: {plain}"
    );
    let with = envelope(&["--include-solution"]);
    assert_eq!(with["solution"]["columns"].as_array().unwrap().len(), 4);
    let opted = envelope(&["--opt", "mip_rel_gap=0.01", "--opt", "output_flag=false"]);
    assert_eq!(
        opted["run"]["options"],
        serde_json::json!({"mip_rel_gap": "0.01", "output_flag": "false"}),
        "--opt pairs are echoed as given, values as strings"
    );
}

/// With the HiGHS log switched on, stdout is still exactly one JSON object
/// and the whole log is on stderr.
#[test]
fn envelope_stdout_is_one_json_object_with_logging_on() {
    let model = write_tmp("tiny.mps", TINY_MIP);
    let (code, stdout, stderr) = mip(
        &model,
        &[
            "--format",
            "json",
            "--opt",
            "output_flag=true",
            "--solution",
        ],
    );
    assert_eq!(code, 0, "stderr: {stderr}");
    let env: serde_json::Value = serde_json::from_str(stdout.trim())
        .unwrap_or_else(|e| panic!("stdout is not one JSON object: {e}\n{stdout}"));
    assert_eq!(env["status"], "Optimal");
    assert_eq!(env["solution"]["columns"].as_array().unwrap().len(), 4);
    assert!(
        stderr.contains("Solving report"),
        "the HiGHS solving report did not reach stderr: {stderr}"
    );
    assert!(
        !stdout.contains("Solving report"),
        "the log leaked into stdout: {stdout}"
    );
}

/// The raw HiGHS solution file lists the MPS names, not `NoName`, and
/// carries the primal solution, the MIP bound and gap, and the row activities.
#[test]
fn envelope_solution_file_lists_the_mps_names() {
    let model = write_tmp("tiny.mps", TINY_MIP);
    let sol = tmp_dir().join("names.sol");
    let _ = std::fs::remove_file(&sol);
    let (code, _, stderr) = mip(&model, &["--solution-file", sol.to_str().unwrap()]);
    assert_eq!(code, 0, "stderr: {stderr}");
    let text = std::fs::read_to_string(&sol).unwrap();
    assert!(
        !text.contains("NoName"),
        "names did not reach the file:\n{text}"
    );
    for line in [
        "Model status",
        "Optimal",
        "# Primal solution values",
        "Feasible",
        "Objective -20",
        "Objective_bound -20",
        "MIP_rel_gap 0",
        "# Columns 4",
        "x 4",
        "y 2",
        "z 0",
        "w 0",
        "# Rows 2",
        "cap1 8",
        "cap2 0",
        "# Dual solution values",
        "# Basis",
        "HiGHS_basis_file v2",
    ] {
        assert!(
            text.contains(line),
            "solution file has no {line:?}:\n{text}"
        );
    }
}

/// A `.json` solution file is the `solution` object of the envelope.
#[test]
fn envelope_solution_file_json_is_the_solution_object() {
    let model = write_tmp("tiny.mps", TINY_MIP);
    let sol = tmp_dir().join("tiny.sol.json");
    let _ = std::fs::remove_file(&sol);
    let (code, _, stderr) = mip(
        &model,
        &["--format", "json", "--solution-file", sol.to_str().unwrap()],
    );
    assert_eq!(code, 0, "stderr: {stderr}");
    let written: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sol).unwrap()).expect("valid JSON");
    assert_eq!(
        written,
        envelope(&["--solution"])["solution"],
        "the .json solution file is the envelope's solution object"
    );
}

/// Exit 0 whenever the solve ran (whatever the model status), 2 for usage
/// errors, 3 for an input error.
#[test]
fn envelope_exit_codes_separate_ran_usage_and_input() {
    let model = write_tmp("tiny.mps", TINY_MIP);
    for mode in [&["--csv"][..], &[], &["--format", "json"][..]] {
        let (code, stdout, stderr) = mip(&model, mode);
        assert_eq!(code, 0, "`rusolve mip tiny.mps {mode:?}`: {stdout}{stderr}");
    }
    // A node limit that stops the search is a result, not a failure: the status
    // is in the output, not in the exit code.
    let (code, stdout, _) = mip(&model, &["--node-limit", "1", "--format", "json"]);
    assert_eq!(code, 0, "stdout: {stdout}");
    let status: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    println!("node-limit 1 gives status {}", status["status"]);

    for (args, why) in [
        (vec!["--threads", "0"], "thread count below 1"),
        (vec!["--no-such-flag"], "unsupported flag"),
        (vec!["--time-limit", "soon"], "value that will not parse"),
        (vec!["--opt", "mip_rel_gap"], "option with no value"),
        (vec!["--opt", "no_such_option=1"], "unknown option"),
        (vec!["--format", "yaml"], "unknown format"),
        (vec!["--solution-file"], "missing path"),
        (
            vec!["--csv", "--format", "json"],
            "two output modes at once",
        ),
    ] {
        let (code, _, _) = mip(&model, &args);
        assert_eq!(code, 2, "{why} {args:?} exit code");
    }

    let missing = tmp_dir().join("definitely-not-here.mps");
    let (code, _, _) = mip(&missing, &[]);
    assert_eq!(code, 3, "a missing model file");
    let broken = write_tmp("broken.mps", "NAME B\nROWS\n N  o\n ZZZ r\nENDATA\n");
    let (code, _, _) = mip(&broken, &[]);
    assert_eq!(code, 3, "a model that does not parse");
}

/// `--csv` is untouched - the header line and the seven fields the harness
/// parsers read, and nothing else on stdout.
#[test]
fn envelope_csv_contract_is_unchanged() {
    let model = write_tmp("tiny.mps", TINY_MIP);
    let (code, stdout, _) = mip(&model, &["--csv"]);
    assert_eq!(code, 0);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2, "stdout gained a line: {stdout:?}");
    assert_eq!(lines[0], "status,objective,bound,nodes,lp_iters,time_s,x");
    let f: Vec<&str> = lines[1].split(',').collect();
    assert_eq!(f.len(), 7, "csv field count: {}", lines[1]);
    assert_eq!(f[0], "Optimal");
    assert_eq!(f[1], "-20");
    assert_eq!(f[2], "-20");
    assert_eq!(
        f[6], "[4;2;-0;0]",
        "x is still the positional value blob, signed zeros and all"
    );
}

/// The status enum of the committed schema is exactly the set the port can
/// report, so a consumer switching on it sees no surprises. The per-variant side
/// of this check is `envelope_status_mapping_covers_every_model_status` in
/// `src/envelope.rs`.
#[test]
fn envelope_reports_a_status_the_schema_allows() {
    let schema: serde_json::Value = serde_json::from_str(SCHEMA).unwrap();
    let allowed: Vec<String> = schema["properties"]["status"]["enum"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect();
    for extra in [&[][..], &["--node-limit", "1"][..]] {
        let status = envelope(extra)["status"].as_str().unwrap().to_string();
        assert!(
            allowed.contains(&status),
            "{status} is not in the schema's status enum {allowed:?}"
        );
    }
}

/// HiGHS's `kSolutionLimit` covers the node, leaf, improving-solution and stall
/// limits. The envelope says `NodeLimit` only when the node limit is what
/// stopped the solve, so an agent that raises `--node-limit` in reply is not sent
/// after the wrong knob.
#[test]
fn envelope_separates_the_node_limit_from_the_other_solution_limits() {
    let status = |extra: &[&str]| envelope(extra)["status"].as_str().unwrap().to_string();
    assert_eq!(status(&["--node-limit", "0"]), "NodeLimit");
    assert_eq!(status(&["--opt", "mip_max_leaves=0"]), "SolutionLimit");
}

/// `input.sha256` is the SHA-256 of the model file's bytes as read, so for a
/// `.gz` model it is the hash of the compressed file.
#[test]
fn envelope_input_sha256_is_the_hash_of_the_file_bytes() {
    use sha2::{Digest, Sha256};
    let hex = |bytes: &[u8]| -> String {
        Sha256::digest(bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect()
    };
    // the published SHA-256 test vector for "abc" pins the helper itself
    assert_eq!(
        hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    let env = envelope(&[]);
    assert_eq!(env["input"]["sha256"], hex(TINY_MIP.as_bytes()));
    // a gzip file holding the model: the hash is of the .gz bytes on disk
    let plain = write_tmp("tiny_for_gz.mps", TINY_MIP);
    let gz = tmp_dir().join("tiny_for_gz.mps.gz");
    let _ = std::fs::remove_file(&gz);
    let status = Command::new("gzip")
        .arg("-kf")
        .arg(&plain)
        .status()
        .expect("gzip available");
    assert!(status.success());
    let (code, stdout, stderr) = mip(&gz, &["--format", "json"]);
    assert_eq!(code, 0, "stderr: {stderr}");
    let env: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(env["input"]["sha256"], hex(&std::fs::read(&gz).unwrap()));
    assert_ne!(env["input"]["sha256"], hex(TINY_MIP.as_bytes()));
}

/// `solver.profile` is the `RUSOLVE_PROFILE` in force and `solver.git` is the
/// short hash `build.rs` captured (absent only when the build had no repository).
#[test]
fn envelope_names_the_profile_and_the_commit() {
    let env = envelope(&[]);
    assert_eq!(
        env["solver"]["profile"], "recommended",
        "the default profile"
    );
    let model = write_tmp("tiny.mps", TINY_MIP);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .env("RUSOLVE_PROFILE", "highs-1.15.1")
        .args(["mip", model.to_str().unwrap(), "--format", "json"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0));
    let other: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(other["solver"]["profile"], "highs-1.15.1");
    let commit = option_env!("RUSOLVE_GIT");
    assert_eq!(
        env["solver"].get("git").and_then(|g| g.as_str()),
        commit,
        "solver.git is what the build script captured"
    );
    if let Some(commit) = commit {
        let hash = commit.strip_suffix("-dirty").unwrap_or(commit);
        assert_eq!(hash.len(), 12, "short=12: {commit}");
        assert!(hash.chars().all(|c| c.is_ascii_hexdigit()), "{commit}");
    }
}

/// The same model with the same options gives the same envelope, wall time
/// aside.
#[test]
fn envelope_rerun_is_identical() {
    let strip = |mut env: serde_json::Value| {
        strip_times(&mut env);
        env
    };
    for extra in [
        &["--solution"][..],
        &[
            "--seed",
            "7",
            "--opt",
            "mip_rel_gap=0.01",
            "--node-limit",
            "5",
        ][..],
    ] {
        let first = strip(envelope(extra));
        let second = strip(envelope(extra));
        assert_eq!(first, second, "two runs of {extra:?} differ");
    }
}

/// A multi-knapsack with 30 binaries and 4 rows from a fixed LCG: small enough to solve in
/// well under a second, big enough that the root LP is fractional so the LP, the cuts and the
/// tree all run (the work and phase tests need a solve that does work).
fn knapsack_mps() -> String {
    let (n, m) = (30usize, 4usize);
    let mut state: u64 = 12345;
    let mut next = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as i64
    };
    let a: Vec<Vec<i64>> = (0..m)
        .map(|_| (0..n).map(|_| 5 + next() % 40).collect())
        .collect();
    let profit: Vec<i64> = (0..n).map(|_| 10 + next() % 50).collect();
    let mut s = String::from("NAME KNAP\nROWS\n N  obj\n");
    for r in 0..m {
        s += &format!(" L  c{r}\n");
    }
    s += "COLUMNS\n    MARKER  'MARKER'  'INTORG'\n";
    for j in 0..n {
        s += &format!("    x{j}  obj  {}", -profit[j]);
        for (r, row) in a.iter().enumerate() {
            s += &format!("  c{r}  {}", row[j]);
        }
        s += "\n";
    }
    s += "    MARKER  'MARKER'  'INTEND'\nRHS\n";
    for (r, row) in a.iter().enumerate() {
        let cap: i64 = row.iter().sum::<i64>() / 3;
        s += &format!("    rhs  c{r}  {cap}\n");
    }
    s += "BOUNDS\n";
    for j in 0..n {
        s += &format!(" UP bnd x{j} 1\n");
    }
    s + "ENDATA\n"
}

fn knapsack_envelope(extra: &[&str]) -> serde_json::Value {
    let model = write_tmp("knap.mps", &knapsack_mps());
    let (code, stdout, stderr) = mip(&model, &[&["--format", "json"][..], extra].concat());
    assert_eq!(code, 0, "exit code, stderr: {stderr}");
    serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not one JSON object: {e}\n{stdout}"))
}

/// Removes every `*_s` (wall seconds) key at any depth: what is left must repeat exactly.
fn strip_times(value: &mut serde_json::Value) {
    match value {
        serde_json::Value::Object(map) => {
            map.retain(|k, _| !k.ends_with("_s") && k != "primal_dual_integral"); // the integral is over wall time
            map.values_mut().for_each(strip_times);
        }
        serde_json::Value::Array(items) => items.iter_mut().for_each(strip_times),
        _ => {}
    }
}

/// Checks `node` against the `required` lists of the schema node `schema`, recursively through
/// every property that is itself an object schema with a `required` list.
fn assert_schema_required(node: &serde_json::Value, schema: &serde_json::Value, path: &str) {
    if node.is_null() {
        return;
    }
    if let Some(required) = schema["required"].as_array() {
        for field in required {
            let field = field.as_str().unwrap();
            assert!(
                node.get(field).is_some(),
                "{path}.{field} is required by the schema and missing from {node}"
            );
        }
    }
    if let Some(props) = schema["properties"].as_object() {
        for (key, sub) in props {
            if let Some(child) = node.get(key).filter(|c| c.is_object()) {
                assert_schema_required(child, sub, &format!("{path}.{key}"));
            }
        }
    }
}

/// `stats.work` is a positive integer for a solve that did LP work, and the same run gives
/// the same number every time (it is a count of operations, not a clock).
#[test]
fn envelope_work_is_a_deterministic_integer() {
    let first = knapsack_envelope(&[]);
    let work = first["stats"]["work"]
        .as_u64()
        .expect("stats.work is an integer");
    assert!(work > 0, "a solve with LP iterations has work: {first}");
    for _ in 0..2 {
        let again = knapsack_envelope(&[]);
        assert_eq!(
            again["stats"]["work"].as_u64(),
            Some(work),
            "work moved between runs"
        );
    }
    // Another seed is another search, so it may differ; the same seed must not.
    let a = knapsack_envelope(&["--seed", "3"]);
    let b = knapsack_envelope(&["--seed", "3"]);
    assert_eq!(a["stats"]["work"], b["stats"]["work"]);
}

/// `stats.phases` has every field the schema requires, with the right types, and all of it
/// except the `*_s` times repeats exactly between runs.
#[test]
fn envelope_phases_have_the_schema_fields_and_repeat() {
    let schema: serde_json::Value = serde_json::from_str(SCHEMA).unwrap();
    let env = knapsack_envelope(&[]);
    let p = &env["stats"]["phases"];
    assert!(p.is_object(), "stats.phases is missing: {env}");
    assert_schema_required(
        p,
        &schema["properties"]["stats"]["properties"]["phases"],
        "stats.phases",
    );
    assert!(p["presolve"]["time_s"].as_f64().unwrap() >= 0.0);
    assert!(p["presolve"]["status"].is_string());
    assert!(p["presolve"]["reductions"].is_object());
    assert_eq!(p["presolve"]["rows_before"], 4);
    assert_eq!(p["presolve"]["cols_before"], 30);
    assert!(p["root"]["lp_iterations"].as_i64().unwrap() > 0);
    assert!(p["cuts"]["rounds"].as_i64().unwrap() >= 0);
    let separators = p["cuts"]["separators"].as_array().unwrap();
    assert!(!separators.is_empty());
    for s in separators {
        assert!(s["name"].is_string() && s["calls"].is_i64() && s["cuts_found"].is_i64());
    }
    assert!(p["heuristics"]["improving_incumbents_by_source"].is_array());
    assert_eq!(p["search"]["nodes"], env["stats"]["nodes"]);
    assert_eq!(p["lp_iterations"]["total"], env["stats"]["lp_iterations"]);
    let lp = &p["lp_iterations"];
    let parts: i64 = ["strong_branching", "separation", "heuristics"]
        .iter()
        .map(|k| lp[k].as_i64().unwrap())
        .sum();
    assert!(parts <= lp["total"].as_i64().unwrap(), "{lp}");
    let first = &p["first_incumbent"];
    assert!(first.is_object(), "the knapsack has an incumbent: {p}");
    assert!(first["work"].as_u64().unwrap() <= env["stats"]["work"].as_u64().unwrap());
    assert!(p["total_time_s"].as_f64().unwrap() >= p["presolve"]["time_s"].as_f64().unwrap());

    let (mut a, mut b) = (env.clone(), knapsack_envelope(&[]));
    strip_times(&mut a);
    strip_times(&mut b);
    assert_eq!(a, b, "phases (apart from times) differ between two runs");
}
