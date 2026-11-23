//! One command/option registry ( the value kinds behind every flag, the `--help --json`
//! types, the MCP `solve` options against the option registry), `--help` and
//! `--help --json`, and one CLI across the HiGHS port and the native stack. Every test
//! runs the real binary.

use serde_json::{Value, json};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const HELP_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-help-1.schema.json");
const RESULT_SCHEMA: &str = include_str!("../../../docs/schema/rusolve-result-1.schema.json");

/// min -3x - 4y + 5z + w, x + 2y + z <= 8, z + 3w <= 10, x in [0,4], y in [0,5], z in [0,3]
/// integer, w continuous: optimum -20 at x=4 y=2.
const TINY_MPS: &str = "\
NAME A1
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

/// The same model in LP format.
const TINY_LP: &str = "\
\\ LP format for model a1
Minimize
 obj: - 3 x - 4 y + 5 z + w
Subject To
 cap1: x + 2 y + z <= 8
 cap2: z + 3 w <= 10
Bounds
 0 <= x <= 4
 0 <= y <= 5
 0 <= z <= 3
General
 x y z
End
";

const INFEASIBLE_LP: &str = "\
\\ LP format for model inf
Minimize
 obj: x
Subject To
 lo: x >= 2
 hi: x <= 1
End
";

fn tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-a1-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn gzip(path: &PathBuf) -> PathBuf {
    let gz = PathBuf::from(format!("{}.gz", path.display()));
    let out = Command::new("gzip")
        .arg("-kf")
        .arg(path)
        .output()
        .expect("gzip");
    assert!(out.status.success());
    gz
}

/// `(exit code, stdout, stderr)` of `rusolve <args>`, on the port when `port`, else on the native search (`RUSOLVE_MIP_HIGHS=0`).
fn run(port: bool, args: &[&str]) -> (i32, String, String) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    if port {
        cmd.env("RUSOLVE_MIP_HIGHS", "1");
    } else {
        cmd.env("RUSOLVE_MIP_HIGHS", "0");
    }
    let out = cmd.args(args).output().unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

fn json_of(stdout: &str) -> Value {
    serde_json::from_str(stdout.trim())
        .unwrap_or_else(|e| panic!("not one JSON object ({e}): {stdout}"))
}

// ------------------------------------------------------------------------- help

#[test]
fn cli_help_top_level_lists_every_command_and_exits_zero() {
    for flag in ["--help", "-h"] {
        let (code, out, _) = run(false, &[flag]);
        assert_eq!(code, 0, "{flag}");
        assert!(out.starts_with("usage: rusolve"), "{out}");
        let help = json_of(&run(false, &["--help", "--json"]).1);
        for c in help["commands"].as_array().unwrap() {
            let name = c["name"].as_str().unwrap();
            assert!(out.contains(name), "{flag} does not list {name}");
        }
    }
}

#[test]
fn cli_help_every_subcommand_prints_usage_and_never_reads_the_flag_as_a_path() {
    let help = json_of(&run(false, &["--help", "--json"]).1);
    let commands = help["commands"].as_array().unwrap();
    assert!(commands.len() >= 15, "{} commands", commands.len());
    for c in commands {
        let name = c["name"].as_str().unwrap();
        let (code, out, err) = run(false, &[name, "--help"]);
        assert_eq!(code, 0, "{name} --help: {err}");
        assert!(
            out.starts_with(&format!("usage: rusolve {name}")),
            "{name} --help: {out}"
        );
        assert!(
            err.is_empty(),
            "{name} --help wrote to stderr (read --help as a path?): {err}"
        );
        // every registered flag is documented in the usage text
        for f in c["flags"].as_array().unwrap() {
            assert!(
                out.contains(f["name"].as_str().unwrap()),
                "{name} --help misses {}",
                f["name"]
            );
        }
        // and `-h`, and the JSON form of the one command
        assert_eq!(run(false, &[name, "-h"]).0, 0);
        let one = json_of(&run(false, &[name, "--help", "--json"]).1);
        assert_eq!(one["commands"].as_array().unwrap().len(), 1);
        assert_eq!(one["commands"][0]["name"], name);
    }
}

fn required(node: &Value) -> Vec<String> {
    node["required"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect()
}

#[test]
fn cli_help_json_carries_every_required_field_of_the_committed_schema() {
    let schema: Value = serde_json::from_str(HELP_SCHEMA).unwrap();
    let help = json_of(&run(false, &["--help", "--json"]).1);
    assert_eq!(help["schema_version"], "rusolve.help/1");
    assert_eq!(
        schema["properties"]["schema_version"]["const"],
        "rusolve.help/1"
    );
    let has = |node: &Value, what: &str, fields: Vec<String>| {
        for f in fields {
            assert!(
                node.get(&f).is_some(),
                "{what}.{f} is required by the schema: {node}"
            );
        }
    };
    has(&help, "top", required(&schema));
    let cmd_schema = &schema["properties"]["commands"]["items"];
    let pos_schema = &cmd_schema["properties"]["positionals"]["items"];
    let flag_schema = &cmd_schema["properties"]["flags"]["items"];
    let enum_stack: Vec<&str> = flag_schema["properties"]["stack"]["enum"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    for c in help["commands"].as_array().unwrap() {
        has(c, "command", required(cmd_schema));
        for p in c["positionals"].as_array().unwrap() {
            has(p, "positional", required(pos_schema));
        }
        for f in c["flags"].as_array().unwrap() {
            has(f, "flag", required(flag_schema));
            assert!(f["name"].as_str().unwrap().starts_with("--"));
            assert!(enum_stack.contains(&f["stack"].as_str().unwrap()), "{f}");
        }
    }
    for code in required(&schema["properties"]["exit_codes"]) {
        assert!(help["exit_codes"].get(&code).is_some(), "exit code {code}");
    }
    // the `mip` entry is the stack table: shared flags on both stacks, the rest per stack
    let mip = help["commands"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["name"] == "mip")
        .unwrap();
    let stack_of = |flag: &str| {
        mip["flags"]
            .as_array()
            .unwrap()
            .iter()
            .find(|f| f["name"] == flag)
            .unwrap_or_else(|| panic!("mip has no flag {flag}"))["stack"]
            .clone()
    };
    assert_eq!(stack_of("--threads"), "both");
    assert_eq!(stack_of("--format"), "both");
    assert_eq!(stack_of("--opt"), "port");
    assert_eq!(stack_of("--no-cuts"), "native");
}

#[test]
fn cli_help_unknown_command_and_bare_invocation_are_usage_errors() {
    let (code, _, err) = run(false, &["frobnicate"]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("frobnicate"));
    let (code, out, err) = run(false, &[]);
    assert_eq!(code, 2);
    assert!(out.is_empty());
    assert!(err.contains("commands:"));
}

// ------------------------------------------------------------------------- one CLI

/// The flag set every test below hands both stacks: the `both` flags of `commands::MIP_FLAGS`.
const SHARED_FLAGS: &[&str] = &[
    "--format",
    "json",
    "--solution",
    "--seed",
    "3",
    "--threads",
    "1",
    "--time-limit",
    "60",
    "--node-limit",
    "100000",
    "--no-presolve",
];

fn envelope(port: bool, model: &Path) -> Value {
    let mut args = vec!["mip", model.to_str().unwrap()];
    args.extend_from_slice(SHARED_FLAGS);
    let (code, out, err) = run(port, &args);
    assert_eq!(code, 0, "port={port}: {err}");
    json_of(&out)
}

fn assert_same_envelope(a: &Value, b: &Value, what: &str) {
    // the schema's required fields are all present on both
    let schema: Value = serde_json::from_str(RESULT_SCHEMA).unwrap();
    for env in [a, b] {
        for f in required(&schema) {
            assert!(env.get(&f).is_some(), "{what}: {f} missing from {env}");
        }
    }
    assert_eq!(a["schema_version"], b["schema_version"], "{what}");
    assert_eq!(a["status"], b["status"], "{what}");
    assert_eq!(a["has_solution"], b["has_solution"], "{what}");
    assert_eq!(a["input"], b["input"], "{what}: input");
    for k in ["seed", "threads", "time_limit", "node_limit"] {
        assert_eq!(a["run"][k], b["run"][k], "{what}: run.{k}");
    }
    let close = |x: &Value, y: &Value, what2: &str| match (x.as_f64(), y.as_f64()) {
        (Some(x), Some(y)) => assert!(
            (x - y).abs() <= 1e-6 * (1.0 + x.abs()),
            "{what}: {what2} {x} vs {y}"
        ),
        (None, None) => {}
        _ => panic!("{what}: {what2} {x} vs {y}"),
    };
    close(&a["objective"], &b["objective"], "objective");
    close(&a["dual_bound"], &b["dual_bound"], "dual_bound");
    close(&a["gap"], &b["gap"], "gap");
    let (ca, cb) = (&a["solution"]["columns"], &b["solution"]["columns"]);
    assert_eq!(
        ca.as_array().map(Vec::len),
        cb.as_array().map(Vec::len),
        "{what}"
    );
    for (x, y) in ca.as_array().unwrap().iter().zip(cb.as_array().unwrap()) {
        assert_eq!(x["name"], y["name"], "{what}");
        close(&x["value"], &y["value"], "solution value");
    }
}

#[test]
fn a1_same_flags_same_envelope_on_both_stacks_for_mps_lp_and_gz() {
    let mps = tmp("tiny.mps", TINY_MPS);
    let lp = tmp("tiny.lp", TINY_LP);
    let mps_gz = gzip(&mps);
    let lp_gz = gzip(&lp);
    for model in [&mps, &lp, &mps_gz, &lp_gz] {
        let what = model.file_name().unwrap().to_str().unwrap();
        let native = envelope(false, model);
        let port = envelope(true, model);
        assert_eq!(native["solver"]["engine"], "native", "{what}");
        assert_eq!(port["solver"]["engine"], "highs-port", "{what}");
        assert_eq!(native["status"], "Optimal", "{what}");
        assert!(
            (native["objective"].as_f64().unwrap() + 20.0).abs() < 1e-6,
            "{what}"
        );
        assert_eq!(native["input"]["rows"], 2);
        assert_eq!(native["input"]["integer_cols"], 3);
        assert_eq!(native["solution"]["columns"][0]["name"], "x");
        assert_same_envelope(&native, &port, what);
    }
}

#[test]
fn a1_infeasible_model_is_status_infeasible_and_exit_zero_on_both_stacks() {
    let lp = tmp("inf.lp", INFEASIBLE_LP);
    for port in [false, true] {
        let (code, out, err) = run(port, &["mip", lp.to_str().unwrap(), "--format", "json"]);
        assert_eq!(code, 0, "port={port}: {err}");
        let env = json_of(&out);
        assert_eq!(env["status"], "Infeasible", "port={port}");
        assert_eq!(env["has_solution"], false);
        assert!(env["objective"].is_null());
    }
}

#[test]
fn a1_exit_codes_are_the_a5_table_on_both_stacks() {
    let mps = tmp("codes.mps", TINY_MPS);
    let m = mps.to_str().unwrap();
    for port in [false, true] {
        // 3: input error
        assert_eq!(
            run(port, &["mip", "/nonexistent/model.mps"]).0,
            3,
            "port={port}"
        );
        let junk = tmp("junk.mps", "this is not a model\n");
        assert_eq!(
            run(port, &["mip", junk.to_str().unwrap()]).0,
            3,
            "port={port}"
        );
        // 2: usage errors
        assert_eq!(run(port, &["mip"]).0, 2, "port={port}: missing path");
        assert_eq!(
            run(port, &["mip", m, "--bogus"]).0,
            2,
            "port={port}: unknown flag"
        );
        assert_eq!(
            run(port, &["mip", m, "--seed"]).0,
            2,
            "port={port}: no value"
        );
        assert_eq!(run(port, &["mip", m, "--threads", "0"]).0, 2, "port={port}");
        assert_eq!(
            run(port, &["mip", m, "--format", "xml"]).0,
            2,
            "port={port}"
        );
        assert_eq!(
            run(port, &["mip", m, "--csv", "--format", "json"]).0,
            2,
            "port={port}"
        );
        // 0: a limit stopped the solve; the status is in the output
        let (code, out, _) = run(port, &["mip", m, "--node-limit", "0", "--csv"]);
        assert_eq!(code, 0, "port={port}: {out}");
    }
}

#[test]
fn a1_a_flag_of_the_other_stack_is_a_usage_error_that_names_the_stack() {
    let mps = tmp("stacks.mps", TINY_MPS);
    let m = mps.to_str().unwrap();
    let (code, _, err) = run(false, &["mip", m, "--mip-rel-gap", "0.1"]);
    assert_eq!(code, 2);
    assert!(
        err.contains("HiGHS port") && err.contains("unset RUSOLVE_MIP_HIGHS"),
        "{err}"
    );
    let (code, _, err) = run(true, &["mip", m, "--no-cuts"]);
    assert_eq!(code, 2);
    assert!(
        err.contains("native search") && err.contains("set RUSOLVE_MIP_HIGHS=0"),
        "{err}"
    );
}

#[test]
fn a1_solution_file_json_is_written_by_both_stacks() {
    let mps = tmp("solfile.mps", TINY_MPS);
    for port in [false, true] {
        let path = tmp(&format!("sol-{port}.json"), "");
        let (code, _, err) = run(
            port,
            &[
                "mip",
                mps.to_str().unwrap(),
                "--solution-file",
                path.to_str().unwrap(),
            ],
        );
        assert_eq!(code, 0, "{err}");
        let sol: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        assert_eq!(sol["columns"].as_array().unwrap().len(), 4, "port={port}");
        assert_eq!(sol["columns"][0]["name"], "x");
    }
    // the raw HiGHS format is port-only, said before the solve
    let (code, _, err) = run(
        false,
        &[
            "mip",
            mps.to_str().unwrap(),
            "--solution-file",
            "/tmp/x.sol",
        ],
    );
    assert_eq!(code, 2, "{err}");
}

#[test]
fn a1_solve_json_envelopes_agree_between_the_port_and_the_legacy_engine() {
    // the legacy engine has no integers: the continuous relaxation of the tiny model
    let lp = tmp("solve.lp", &TINY_LP.replace("General\n x y z\n", ""));
    let mps = tmp(
        "solve.mps",
        &TINY_MPS
            .replace("    MARKER  'MARKER'  'INTORG'\n", "")
            .replace("    MARKER  'MARKER'  'INTEND'\n", ""),
    );
    for model in [&lp, &mps] {
        let m = model.to_str().unwrap();
        let (c1, o1, e1) = run(false, &["solve", m, "--format", "json", "--solution"]);
        let (c2, o2, e2) = run(
            false,
            &["solve", m, "--format", "json", "--solution", "--legacy-lp"],
        );
        assert_eq!((c1, c2), (0, 0), "{e1} {e2}");
        let (port, legacy) = (json_of(&o1), json_of(&o2));
        assert_eq!(port["solver"]["engine"], "highs-port");
        assert_eq!(legacy["solver"]["engine"], "native");
        assert_eq!(port["status"], "Optimal");
        assert_eq!(legacy["status"], "Optimal");
        assert_eq!(port["input"], legacy["input"]);
        let (a, b) = (
            port["objective"].as_f64().unwrap(),
            legacy["objective"].as_f64().unwrap(),
        );
        assert!((a - b).abs() < 1e-6, "{a} vs {b}");
    }
}

#[test]
fn a1_solve_is_strict_about_flags_and_exits_by_the_a5_table() {
    let lp = tmp("strict.lp", INFEASIBLE_LP);
    let l = lp.to_str().unwrap();
    assert_eq!(run(false, &["solve", l, "--bogus"]).0, 2);
    assert_eq!(run(false, &["solve", l, "--time-limit"]).0, 2);
    assert_eq!(run(false, &["solve", "/nonexistent.mps"]).0, 3);
    assert_eq!(run(false, &["solve"]).0, 2);
    // infeasible: the solve ran
    for extra in [&[][..], &["--legacy-lp"][..], &["--ipm"][..]] {
        let mut args = vec!["solve", l];
        args.extend_from_slice(extra);
        assert_eq!(run(false, &args).0, 0, "{extra:?}");
    }
}

// ------------------------------------------------------------------------- registry

/// The `solve` tool's `options` property as an MCP client sees it, from `tools/list`.
fn mcp_solve_option_properties() -> Value {
    fn ask(
        stdin: &mut std::process::ChildStdin,
        out: &mut BufReader<std::process::ChildStdout>,
        line: &str,
    ) -> Value {
        writeln!(stdin, "{line}").unwrap();
        stdin.flush().unwrap();
        let mut reply = String::new();
        out.read_line(&mut reply).unwrap();
        serde_json::from_str(&reply).unwrap_or_else(|e| panic!("not JSON ({e}): {reply}"))
    }
    let mut child = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("mcp")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let mut out = BufReader::new(child.stdout.take().unwrap());
    ask(
        &mut stdin,
        &mut out,
        &json!({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
            "protocolVersion": "2025-11-25", "capabilities": {},
            "clientInfo": {"name": "a0", "version": "0"}}})
        .to_string(),
    );
    // a notification, so no line to read back
    writeln!(
        stdin,
        "{}",
        json!({"jsonrpc": "2.0", "method": "notifications/initialized"})
    )
    .unwrap();
    stdin.flush().unwrap();
    let list = ask(
        &mut stdin,
        &mut out,
        &json!({"jsonrpc": "2.0", "id": 2, "method": "tools/list", "params": {}}).to_string(),
    );
    let solve = list["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .find(|t| t["name"] == "solve")
        .expect("the solve tool");
    let properties = solve["inputSchema"]["properties"]["options"]["properties"].clone();
    assert!(
        properties.is_object(),
        "solve.options has no properties: {solve}"
    );
    drop(stdin);
    let _ = child.wait();
    properties
}

/// The registry contract on the options side, end to end: what `rusolve options --json` lists is
/// exactly what the MCP `solve` tool's `options` property offers, name for name and type for
/// type, because both read the option registry (`commands::option_properties`).
#[test]
fn a0_the_mcp_solve_options_are_the_option_registry_names_and_types() {
    let (code, stdout, stderr) = run(false, &["options", "--json"]);
    assert_eq!(code, 0, "{stderr}");
    let records = json_of(&stdout);
    let properties = mcp_solve_option_properties();
    let listed = records.as_array().expect("an array of option records");
    let object = properties.as_object().unwrap();
    assert!(listed.len() > 200, "{} options listed", listed.len());
    assert_eq!(
        object.len(),
        listed.len(),
        "the MCP schema and `options --json` list a different number of options"
    );
    // the registry type as JSON Schema: a string is accepted for every type, as the HiGHS
    // setters are, so the permissive pair is part of the one source too
    let schema_type = |type_name: &str| match type_name {
        "bool" => json!(["boolean", "string"]),
        "HighsInt" => json!(["integer", "string"]),
        "double" => json!(["number", "string"]),
        "string" => json!("string"),
        other => panic!("option type {other} is not in the registry's vocabulary"),
    };
    let mut missing = Vec::new();
    for record in listed {
        let name = record["name"].as_str().unwrap();
        let want = schema_type(record["type_name"].as_str().unwrap());
        match object.get(name) {
            // the property is `{"type": ...}`: the same type the registry states
            Some(p) => assert_eq!(p["type"], want, "{name}"),
            None => missing.push(name.to_string()),
        }
    }
    assert!(missing.is_empty(), "{missing:?} are not in solve.options");
    let listed_names: Vec<&str> = listed.iter().map(|r| r["name"].as_str().unwrap()).collect();
    let extra: Vec<&String> = object
        .keys()
        .filter(|k| !listed_names.contains(&k.as_str()))
        .collect();
    assert!(
        extra.is_empty(),
        "solve.options offers {extra:?} that the registry has not"
    );
}

/// On the CLI side, through the real binary: the port `mip` parser is the registry, so
/// each way a flag can be wrong is exit 2 naming the flag and what it expects.
#[test]
fn a0_the_port_mip_parser_refuses_by_the_registry_and_exits_2() {
    let mps = tmp("a0.mps", TINY_MPS);
    let m = mps.to_str().unwrap();
    for (args, flag, expects) in [
        (vec!["--node-limit", "abc"], "--node-limit", "an integer"),
        (
            vec!["--format", "yaml"],
            "--format",
            "one of text, csv or json",
        ),
        (vec!["--no-such-flag"], "--no-such-flag", "unknown flag"),
        (vec!["--time-limit"], "--time-limit", "a number"),
        (vec!["--seed", "1", "--seed", "2"], "--seed", "only once"),
        (vec!["--opt", "mip_rel_gap"], "--opt", "key=value"),
        (vec!["--wasm-fuel", "-1"], "--wasm-fuel", "non-negative"),
    ] {
        let mut line = vec!["mip", m];
        line.extend(args.iter().copied());
        let (code, _, err) = run(true, &line);
        assert_eq!(code, 2, "{args:?}: {err}");
        assert!(
            err.contains(flag) && err.contains(expects),
            "{args:?}: {err}"
        );
    }
    // the checks *between* flags stay the command's own, and still exit 2
    let (code, _, err) = run(true, &["mip", m, "--csv", "--format", "json"]);
    assert_eq!(code, 2, "{err}");
    let (code, _, err) = run(true, &["mip", m, "--wasm-fuel", "10"]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("--wasm-hook"), "{err}");
    // a flag of the other stack still names it
    let (code, _, err) = run(true, &["mip", m, "--no-cuts"]);
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("native search"), "{err}");
}

/// In the machine-readable surface: `--help --json` types every flag, and a `Choice`
/// lists its values, so a schema or binding can be generated from it (one
/// source).
#[test]
fn a0_help_json_types_every_flag() {
    let help = json_of(&run(false, &["--help", "--json"]).1);
    let mut typed = 0;
    let mut choices = 0;
    for c in help["commands"].as_array().unwrap() {
        for f in c["flags"].as_array().unwrap() {
            let name = format!("{} {}", c["name"], f["name"]);
            let kind = f["type"]
                .as_str()
                .unwrap_or_else(|| panic!("{name} has no type"));
            assert!(
                [
                    "switch",
                    "int",
                    "uint",
                    "float",
                    "path",
                    "text",
                    "choice",
                    "key=value"
                ]
                .contains(&kind),
                "{name} has type {kind}"
            );
            assert_eq!(
                (kind == "switch"),
                f["value"].is_null(),
                "{name}: type {kind} vs placeholder {}",
                f["value"]
            );
            typed += 1;
            match f.get("choices") {
                Some(list) => {
                    assert_eq!(kind, "choice", "{name} lists choices but is {kind}");
                    assert_eq!(
                        list.as_array()
                            .unwrap()
                            .iter()
                            .map(|v| v.as_str().unwrap())
                            .collect::<Vec<_>>()
                            .join("|"),
                        f["value"].as_str().unwrap(),
                        "{name}"
                    );
                    choices += 1;
                }
                None => assert_ne!(kind, "choice", "{name} is a choice with no list"),
            }
        }
    }
    assert!(typed > 100, "{typed} flags typed");
    assert!(choices >= 10, "{choices} choice flags");
}
