//! `rusolve tune`. Every test runs the real binary on a small generated family.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const TUNE_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-tune-1.schema.json");
const PROFILE_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-profile-1.schema.json");

/// A multi-row 0/1 knapsack from an LCG seeded by `seed`: a family of related instances.
fn knapsack(seed: u64, n: usize, m: usize) -> String {
    let mut x: u64 = 12345 + seed * 7919;
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
        s += &format!("    x{j}  cost  {}\n", -(r() % 40 + 10));
        for (i, row) in w.iter().enumerate() {
            s += &format!("    x{j}  k{i}  {}\n", row[j]);
        }
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

fn family(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-tune-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for seed in 0..6 {
        std::fs::write(dir.join(format!("k{seed}.mps")), knapsack(seed, 18, 3)).unwrap();
    }
    dir
}

fn tune(args: &[&str], stdin: Option<&str>) -> (i32, String, String) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("tune")
        .args(args)
        .env_remove("RUSOLVE_PROFILE")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    if let Some(text) = stdin {
        std::io::Write::write_all(child.stdin.as_mut().unwrap(), text.as_bytes()).unwrap();
    }
    drop(child.stdin.take());
    let out = child.wait_with_output().unwrap();
    (
        out.status.code().unwrap(),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8(out.stderr).unwrap(),
    )
}

fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}

fn assert_required(schema: &str, v: &serde_json::Value) {
    let schema: serde_json::Value = serde_json::from_str(schema).unwrap();
    for k in schema["required"].as_array().unwrap() {
        let k = k.as_str().unwrap();
        assert!(v.get(k).is_some(), "{k} missing: {v}");
    }
}

#[test]
fn tune_json_has_the_schema_fields_and_writes_a_valid_profile() {
    let dir = family("json");
    let out = dir.join("profile.json");
    let (code, stdout, stderr) = tune(
        &[
            p(&dir),
            "--seeds",
            "2",
            "--max-candidates",
            "3",
            "--budget-work",
            "3e5",
            "--format",
            "json",
            "-o",
            p(&out),
        ],
        None,
    );
    assert_eq!(code, 0, "{stderr}");
    let v: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    assert_eq!(v["schema_version"], "rusolve.tune/1");
    assert_required(TUNE_SCHEMA, &v);
    assert_eq!(
        v["split"]["train"].as_u64().unwrap() + v["split"]["holdout"].as_u64().unwrap(),
        6
    );
    assert_required(PROFILE_SCHEMA, &v["profile"]);
    assert!(
        v["profile"]["evidence"]["race"]["blocks_run"]
            .as_u64()
            .unwrap()
            >= 1
    );
    // the evidence an agent proposes candidates from: features and cost per candidate, per
    // training instance
    let rows = v["profile"]["evidence"]["train_instances"]
        .as_array()
        .unwrap();
    assert_eq!(rows.len() as u64, v["split"]["train"].as_u64().unwrap());
    for r in rows {
        assert!(r["features"].as_object().unwrap().len() >= 40, "{r}");
        let costs = r["mean_cost_by_candidate"].as_object().unwrap();
        assert!(costs.contains_key("base:recommended") || costs.contains_key("base:default"));
    }
    // the tiny budget may not reach every instance, but it reaches some
    assert!(rows.iter().any(|r| {
        r["mean_cost_by_candidate"]
            .as_object()
            .unwrap()
            .values()
            .any(|c| c.is_number())
    }));
    // the file it wrote passes the profile validator, and is the same profile
    let file: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    assert_eq!(file, v["profile"]);
    let check = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(["profile", "validate", p(&out)])
        .output()
        .unwrap();
    assert_eq!(check.status.code(), Some(0));
}

#[test]
fn tune_wrong_answer_check_flags_a_false_known_optimum() {
    let dir = family("wrong");
    let known = dir.join("known.txt");
    // no instance has optimum -1, so every optimal answer is "wrong"
    std::fs::write(&known, "=opt= k0 -1\nk1,-1\nk2 -1\nk3 -1\nk4 -1\nk5 -1\n").unwrap();
    let (code, stdout, stderr) = tune(
        &[
            p(&dir),
            "--seeds",
            "1",
            "--max-candidates",
            "2",
            "--budget-work",
            "2e5",
            "--known",
            p(&known),
            "--format",
            "json",
        ],
        None,
    );
    assert_eq!(code, 0, "{stderr}");
    let v: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    assert!(v["wrong_answers"].as_u64().unwrap() >= 1, "{v}");
    let none = dir.join("none.txt");
    std::fs::write(&none, "# nothing known\n").unwrap();
    let (code, stdout, _) = tune(
        &[
            p(&dir),
            "--seeds",
            "1",
            "--max-candidates",
            "2",
            "--budget-work",
            "2e5",
            "--known",
            p(&none),
            "--format",
            "json",
        ],
        None,
    );
    assert_eq!(code, 0);
    let v: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    assert_eq!(v["wrong_answers"], 0, "the solver agrees with itself");
}

#[test]
fn tune_proposals_from_stdin_are_validated_and_raced() {
    let dir = family("prop");
    let good = r#"[{"name": "agent-guess", "options": {"mip_heuristic_effort": 0.3, "mip_detect_symmetry": false}}]"#;
    let (code, stdout, stderr) = tune(
        &[
            p(&dir),
            "--seeds",
            "1",
            "--budget-work",
            "3e5",
            "--candidates",
            "-",
            "--candidates-only",
            "--format",
            "json",
        ],
        Some(good),
    );
    assert_eq!(code, 0, "{stderr}");
    let v: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    let cands = v["profile"]["evidence"]["candidates"].as_array().unwrap();
    assert_eq!(cands.len(), 2, "baseline plus the one proposal");
    assert_eq!(cands[1]["name"], "agent-guess");
    assert_eq!(cands[1]["origin"], "proposal");
    // an invalid proposal fails with errors before any run
    let bad = r#"{"candidates": [{"options": {"mip_heuristic_efort": 0.3}}]}"#;
    let (code, _, stderr) = tune(
        &[p(&dir), "--candidates", "-", "--candidates-only"],
        Some(bad),
    );
    assert_eq!(code, 2, "{stderr}");
    assert!(
        stderr.contains("did you mean mip_heuristic_effort?"),
        "{stderr}"
    );
}

#[test]
fn tune_usage_and_input_errors_use_the_a5_codes() {
    assert_eq!(tune(&[], None).0, 2);
    assert_eq!(tune(&["--bogus"], None).0, 2);
    assert_eq!(tune(&["/nonexistent-dir"], None).0, 3);
}

// ------------------------------------------------------ iterated racing (irace)

fn json_of(text: &str) -> serde_json::Value {
    serde_json::from_str(text.trim()).expect("stdout is one JSON object")
}

/// The race's round list, `evidence.race.rounds`.
fn race_rounds(v: &serde_json::Value) -> Vec<serde_json::Value> {
    v["profile"]["evidence"]["race"]["rounds"]
        .as_array()
        .unwrap_or_else(|| panic!("no rounds in {v}"))
        .clone()
}

fn strip_timestamps(v: &mut serde_json::Value) {
    match v {
        serde_json::Value::Object(map) => {
            map.remove("created");
            for child in map.values_mut() {
                strip_timestamps(child);
            }
        }
        serde_json::Value::Array(items) => items.iter_mut().for_each(strip_timestamps),
        _ => {}
    }
}

/// `--iterations 3`: round 2 must race candidates that change two options (the pair moves the
/// single race could never reach), the JSON must list every round, and a budget that runs out
/// part-way must say so.
#[test]
fn tune_iter_three_rounds_races_pair_moves_and_lists_its_rounds() {
    let dir = family("iter3");
    let (code, stdout, stderr) = tune(
        &[
            p(&dir),
            "--seeds",
            "2",
            "--max-candidates",
            "4",
            "--budget-work",
            "2e12",
            "--iterations",
            "3",
            "--format",
            "json",
        ],
        None,
    );
    assert_eq!(code, 0, "{stderr}");
    let v = json_of(&stdout);
    assert_required(TUNE_SCHEMA, &v);
    let rounds = race_rounds(&v);
    let asked = 3;
    let pairs = v["profile"]["evidence"]["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|c| {
            c["changed_options"].as_u64().unwrap() >= 2
                && c["blocks_evaluated"].as_u64().unwrap_or(0) > 0
        })
        .count();
    assert!(
        pairs >= 1,
        "no candidate changing two options was raced: {v}"
    );
    assert!(
        v["profile"]["evidence"]["candidates"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c["origin"] == "elite_mutation"),
        "the later rounds must add mutated elites: {v}"
    );
    assert_eq!(v["iterations"].as_u64().unwrap(), rounds.len() as u64);
    if rounds.len() < asked {
        // fewer rounds only ever happens because the budget or the blocks ran out, and the JSON
        // has to name the reason, both for the run and for the round that was cut short
        let race = &v["profile"]["evidence"]["race"];
        let stopped = race["stopped"].as_str().unwrap();
        assert!(
            stopped == "work_budget_spent" || stopped == "blocks_exhausted",
            "{v}"
        );
        assert_eq!(rounds.last().unwrap()["stopped"].as_str().unwrap(), stopped);
    } else {
        assert_eq!(rounds.len(), asked);
        assert_eq!(
            v["profile"]["evidence"]["race"]["stopped"]
                .as_str()
                .unwrap(),
            "iterations_done"
        );
        // round 2 and 3 each own their own blocks, every round keeps the baseline, and every
        // round after the first is seeded with new candidates
        for (i, r) in rounds.iter().enumerate() {
            assert_eq!(r["round"].as_u64().unwrap(), (i + 1) as u64);
            assert!(r["candidates_raced"].as_u64().unwrap() >= 2, "{r}");
            assert!(r["blocks_run"].as_u64().unwrap() >= 1, "{r}");
            assert!(
                r["candidates"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|c| c.as_str() == Some("base:recommended")),
                "the baseline must race in every round: {r}"
            );
            if i + 1 < rounds.len() {
                // every round that is not the last one seeded a following round
                assert!(r["next_round_candidates"].as_u64().unwrap() >= 1, "{r}");
            } else {
                assert_eq!(r["next_round_candidates"].as_u64().unwrap(), 0, "{r}");
            }
        }
        assert!(rounds[0]["elites"].as_array().unwrap().len() >= 2);
    }
    // a budget that is too small stops the rounds, and says which round and why
    let (code, stdout, stderr) = tune(
        &[
            p(&dir),
            "--seeds",
            "2",
            "--max-candidates",
            "4",
            "--budget-work",
            "1",
            "--iterations",
            "3",
            "--format",
            "json",
        ],
        None,
    );
    assert_eq!(code, 0, "{stderr}");
    let v = json_of(&stdout);
    assert!(race_rounds(&v).len() < asked, "one block of work: {v}");
    assert_eq!(
        v["iterations"].as_u64().unwrap(),
        race_rounds(&v).len() as u64
    );
    assert_eq!(
        v["profile"]["evidence"]["race"]["stopped"]
            .as_str()
            .unwrap(),
        "work_budget_spent"
    );
    assert_eq!(
        race_rounds(&v).last().unwrap()["stopped"].as_str().unwrap(),
        "work_budget_spent"
    );
}

/// The same command twice gives the same JSON: the race, the round list and the winner are all
/// driven by `--race-seed`, not by process state.
#[test]
fn tune_iter_same_command_twice_gives_identical_json() {
    let dir = family("repeat");
    let mut runs = Vec::new();
    for _ in 0..2 {
        let (code, stdout, stderr) = tune(
            &[
                p(&dir),
                "--seeds",
                "1",
                "--max-candidates",
                "3",
                "--budget-work",
                "2e12",
                "--iterations",
                "2",
                "--elite",
                "2",
                "--race-seed",
                "7",
                "--format",
                "json",
            ],
            None,
        );
        assert_eq!(code, 0, "{stderr}");
        let mut v = json_of(&stdout);
        strip_timestamps(&mut v);
        runs.push(v);
    }
    assert_eq!(runs[0], runs[1], "the same command must repeat itself");
    assert_eq!(runs[0]["iterations"].as_u64().unwrap(), 2);
}

/// The cost metric follows `--jobs`: the deterministic work count only when solves do not overlap.
#[test]
fn tune_iter_cost_metric_follows_jobs() {
    let dir = family("metric");
    let common = [
        p(&dir),
        "--seeds",
        "1",
        "--max-candidates",
        "2",
        "--budget-work",
        "2e12",
        "--format",
        "json",
    ];
    let (code, stdout, stderr) = tune(&common, None);
    assert_eq!(code, 0, "{stderr}");
    let v1 = json_of(&stdout);
    assert_eq!(v1["cost_metric"], "c1_work", "{v1}");
    assert_eq!(
        v1["profile"]["provenance"]["metric"]["kind"]
            .as_str()
            .unwrap(),
        "c1_work"
    );
    // the work-counter run counted real work, not the LP-iteration placeholder
    assert!(
        v1["profile"]["evidence"]["race"]["work_spent"]
            .as_f64()
            .unwrap()
            > 1e5,
        "{v1}"
    );
    let mut jobs2: Vec<&str> = common.to_vec();
    jobs2.extend(["--jobs", "2"]);
    let (code, stdout, stderr) = tune(&jobs2, None);
    assert_eq!(code, 0, "{stderr}");
    let v = json_of(&stdout);
    assert_eq!(v["cost_metric"], "iters_nodes", "{v}");
    assert_eq!(
        v["profile"]["provenance"]["metric"]["kind"]
            .as_str()
            .unwrap(),
        "iters_nodes"
    );
}

#[test]
fn tune_iter_zero_iterations_or_elite_is_a_usage_error() {
    let dir = family("zeros");
    assert_eq!(tune(&[p(&dir), "--iterations", "0"], None).0, 2);
    assert_eq!(tune(&[p(&dir), "--elite", "0"], None).0, 2);
    let (_, _, stderr) = tune(&[p(&dir), "--iterations", "0"], None);
    assert!(stderr.contains("--iterations"), "{stderr}");
    assert_eq!(tune(&[p(&dir), "--iterations", "1"], None).0, 0);
}

#[test]
fn tune_space_lists_the_searched_options_without_instances() {
    let (code, stdout, stderr) = tune(&["--space"], None);
    assert_eq!(code, 0, "{stderr}");
    let v: serde_json::Value = serde_json::from_str(stdout.trim()).unwrap();
    assert_eq!(v["schema_version"], "rusolve.tune-space/1");
    let opts = v["options"].as_array().unwrap();
    assert_eq!(opts.len(), 16);
    assert!(
        opts.iter()
            .all(|o| !o["values"].as_array().unwrap().is_empty())
    );
}
