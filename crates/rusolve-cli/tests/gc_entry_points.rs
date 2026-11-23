//! General-constraint regression: no entry point that reads a model may drop
//! semi-continuous columns or general constraints and report an answer for the wrong model.
//! Each command either refuses with a typed error naming the object (exit != 0, nothing that
//! reads "Optimal" on stdout) or, for `mip`/`check`/`convert`, honours the objects. Since the
//! router (route.rs) `solve` honours them too, by the same route as `mip`.

use serde_json::{Value, json};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};

const GC: &str = "Maximize\n obj: 3 x + 2 y - 10 b\nSubject To\n c1: x + y <= 10\n ind1: b = 0 -> x <= 2\nGeneral Constraints\n g1: m = MAX ( x , y )\nBounds\n 0 <= x <= 8\n 0 <= y <= 8\n m <= 7\nBinary\n b\nEnd\n";
const SEMI: &str = "Maximize\n obj: x\nSubject To\n c1: x + y <= 10\nBounds\n 0 <= y <= 3\n x <= 8\nSemi-continuous\n x\nEnd\n";

fn tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-gcep-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn run(args: &[&str]) -> (Option<i32>, String, String) {
    let o = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .args(args)
        .output()
        .unwrap();
    (
        o.status.code(),
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

/// `cmd` (with its flags) on the model must refuse and name the offending object.
fn refuses(cmd: &[&str], model: &std::path::Path, needle: &str) {
    let mut a: Vec<&str> = vec![cmd[0], model.to_str().unwrap()];
    a.extend_from_slice(&cmd[1..]);
    let (code, out, err) = run(&a);
    assert!(
        matches!(code, Some(c) if c != 0),
        "{cmd:?} must refuse, exit {code:?}\n{out}\n{err}"
    );
    assert!(
        !out.to_lowercase().contains("optimal"),
        "{cmd:?} printed an answer: {out}"
    );
    assert!(
        err.contains(needle),
        "{cmd:?}: stderr must name `{needle}`: {err}"
    );
}

#[test]
fn solve_and_friends_refuse_general_constraints() {
    let m = tmp("gc.lp", GC);
    for cmd in [
        // Plain `solve` now routes general constraints to the `mip` engine (see
        // `solve_routes_general_constraints_like_mip`); only an engine choice refuses
        &["solve", "--legacy-lp"][..],
        &["classify"],
        &["classify", "--format", "json"],
        &["dump"],
        &["explain"],
        &["lint"],
        &["ranging"],
        &["diagnose"],
        &["feasrelax"],
        &["probe"],
        &["presolve-dump"],
        &["mip-trace"],
    ] {
        let needle = if cmd.contains(&"--legacy-lp") {
            "--legacy-lp"
        } else {
            "ind1"
        };
        refuses(cmd, &m, needle);
    }
}

#[test]
fn solve_routes_general_constraints_like_mip() {
    // The same model, the same answer, from both commands
    for (name, model, objective) in [("gc_route.lp", GC, 20.0), ("semi_route.lp", SEMI, 8.0)] {
        let m = tmp(name, model);
        let mut seen = Vec::new();
        for cmd in ["solve", "mip"] {
            let (code, out, err) = run(&[cmd, m.to_str().unwrap(), "--format", "json"]);
            assert_eq!(code, Some(0), "{cmd} {name}: {err}");
            let v: Value = serde_json::from_str(&out).unwrap();
            assert_eq!(v["status"], "Optimal", "{cmd} {name}");
            seen.push(v["objective"].as_f64().unwrap());
        }
        assert!(
            (seen[0] - objective).abs() < 1e-6 && (seen[1] - objective).abs() < 1e-6,
            "{name}: {seen:?}"
        );
    }
}

#[test]
fn solve_and_friends_refuse_semi_continuous_columns() {
    let m = tmp("semi.lp", SEMI);
    for cmd in [&["solve", "--legacy-lp"][..], &["dump"], &["ranging"]] {
        let needle = if cmd.contains(&"--legacy-lp") {
            "--legacy-lp"
        } else {
            "semi-continuous"
        };
        refuses(cmd, &m, needle);
    }
}

#[test]
fn profile_features_refuse_general_constraints() {
    // `profile` has no model-reading subcommand; its feature vector feeds `mip --profile`
    // and `tune`, which must not featurise a model they cannot solve natively.
    let m = tmp("gc_prof.lp", GC);
    let (code, out, _) = run(&["tune", m.to_str().unwrap()]);
    assert!(!matches!(code, Some(0)) || !out.to_lowercase().contains("optimal"));
}

#[test]
fn mip_is_the_honest_route_for_general_constraints() {
    let m = tmp("gc_mip.lp", GC);
    let (code, out, err) = run(&["mip", m.to_str().unwrap(), "--format", "json"]);
    assert_eq!(code, Some(0), "{err}");
    let v: Value = serde_json::from_str(&out).unwrap();
    assert_eq!(v["status"], "Optimal");
    assert!(
        (v["objective"].as_f64().unwrap() - 20.0).abs() < 1e-6,
        "{out}"
    );
    assert_eq!(v["general_constraints"]["verified"], true);
}

#[test]
fn convert_keeps_the_objects_or_names_what_it_cannot() {
    let m = tmp("gc_conv.lp", GC);
    let out = tmp("gc_conv_out.lp", "");
    let (code, _, err) = run(&["convert", m.to_str().unwrap(), out.to_str().unwrap()]);
    assert_eq!(code, Some(0), "{err}");
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(
        text.contains("->") && text.to_lowercase().contains("g1: m = max"),
        "{text}"
    );
    // MPS has no function constraint: a typed error, not a file without g1
    let mps = tmp("gc_conv_out.mps", "");
    let (code, _, err) = run(&["convert", m.to_str().unwrap(), mps.to_str().unwrap()]);
    assert!(matches!(code, Some(c) if c != 0), "{err}");
    assert!(err.contains("g1"), "{err}");
    // semi columns survive a conversion
    let s = tmp("semi_conv.lp", SEMI);
    let so = tmp("semi_conv_out.lp", "");
    let (code, _, err) = run(&["convert", s.to_str().unwrap(), so.to_str().unwrap()]);
    assert_eq!(code, Some(0), "{err}");
    assert!(
        std::fs::read_to_string(&so)
            .unwrap()
            .to_lowercase()
            .contains("semi")
    );
}

struct Mcp {
    child: std::process::Child,
    stdin: std::process::ChildStdin,
    out: BufReader<std::process::ChildStdout>,
}

impl Mcp {
    fn call(&mut self, id: i64, method: &str, params: Value) -> Value {
        writeln!(
            self.stdin,
            "{}",
            json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params})
        )
        .unwrap();
        self.stdin.flush().unwrap();
        let mut l = String::new();
        self.out.read_line(&mut l).unwrap();
        serde_json::from_str(&l).unwrap()
    }
}

#[test]
fn mcp_load_model_routes_general_constraints() {
    let m = tmp("gc_mcp.lp", GC);
    let mut child = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("mcp")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let stdin = child.stdin.take().unwrap();
    let out = BufReader::new(child.stdout.take().unwrap());
    let mut c = Mcp { child, stdin, out };
    c.call(
        1,
        "initialize",
        json!({"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "t", "version": "0"}}),
    );
    let r = c.call(
        2,
        "tools/call",
        json!({"name": "load_model", "arguments": {"path": m.to_str().unwrap()}}),
    );
    // The MCP load routes a general-constraint model like `solve` does, so the handle
    // loads with its class and solves; the LP-only tools refuse it, naming the object's class.
    let res = &r["result"];
    assert_eq!(res["isError"], json!(false), "{r}");
    let class = res["structuredContent"]["class"].as_str().unwrap();
    assert!(["gc", "gc_mip"].contains(&class), "{r}");
    let id = res["structuredContent"]["model_id"].as_str().unwrap();
    let r = c.call(
        3,
        "tools/call",
        json!({"name": "lint", "arguments": {"model_id": id}}),
    );
    assert_eq!(r["result"]["isError"], json!(true), "{r}");
    assert!(r["result"].to_string().contains(class), "{r}");
    let _ = c.child.kill();
}

// Mixed classes (general constraints with a quadratic objective or quadratic rows): the combination is
// refused with a typed error, never answered for a model with a part dropped.
const IND_QROW: &str = "Minimize\n obj: x + 3 b\nSubject To\n q: x + [ y ^ 2 ] <= 9\n i: b = 0 -> x >= 5\nBounds\n 0 <= x <= 10\n 0 <= y <= 3\nBinaries\n b\nEnd\n";
const IND_QOBJ: &str = "Minimize\n obj: x + 3 b + [ x ^ 2 ] / 2\nSubject To\n c: x + b >= 1\n i: b = 0 -> x >= 5\nBounds\n 0 <= x <= 10\nBinaries\n b\nEnd\n";
const IND_QIND: &str = "Minimize\n obj: x + 3 b\nSubject To\n c: x + b >= 1\n i: b = 0 -> x + [ x ^ 2 ] >= 5\nBounds\n 0 <= x <= 10\nBinaries\n b\nEnd\n";
const SEMI_QOBJ: &str = "Minimize\n obj: x + [ x ^ 2 ] / 2\nSubject To\n c: x + y >= 1\nBounds\n 0 <= y <= 3\n 2 <= x <= 8\nSemi-continuous\n x\nEnd\n";

#[test]
fn mixed_classes_are_typed_errors() {
    for (tag, text) in [
        ("ind_qrow", IND_QROW),
        ("ind_qobj", IND_QOBJ),
        ("ind_qind", IND_QIND),
        ("semi_qobj", SEMI_QOBJ),
    ] {
        let m = tmp(&format!("{tag}.lp"), text);
        for cmd in [&["mip"][..], &["solve"]] {
            let mut a = vec![cmd[0], m.to_str().unwrap()];
            a.extend_from_slice(&cmd[1..]);
            let (code, out, err) = run(&a);
            assert!(
                matches!(code, Some(c) if c != 0),
                "{tag} {cmd:?} must refuse, exit {code:?}\n{out}\n{err}"
            );
            assert!(
                !out.to_lowercase().contains("optimal"),
                "{tag} {cmd:?} printed an answer: {out}"
            );
            assert!(!err.trim().is_empty(), "{tag} {cmd:?}: no message");
        }
    }
}
