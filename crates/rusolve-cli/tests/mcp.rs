//! The MCP server over stdio, and cheap modify and re-solve on model handles.
//!
//! Every test spawns the real `rusolve-cli mcp`, speaks newline-delimited
//! JSON-RPC to it, and checks what an MCP client would see. `Schema` below is a
//! small JSON Schema checker (type, required, properties, items, enum, const,
//! minimum, oneOf/anyOf, $ref, additionalProperties) used to validate every
//! tool's `inputSchema`/`outputSchema` and every `structuredContent` against
//! its `outputSchema`.

use serde_json::{Value, json};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

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

/// x + y >= 10 and x + y <= 5: infeasible.
const INFEASIBLE: &str = "\
NAME INF3
ROWS
 N  obj
 G  r1
 L  r2
COLUMNS
    x  obj  1  r1  1
    x  r2  1
    y  obj  1  r1  1
    y  r2  1
RHS
    rhs  r1  10  r2  5
ENDATA
";

fn afiro_path() -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../benches/data/netlib/afiro.mps")
        .to_string_lossy()
        .into_owned()
}

/// afiro with extra bound lines (afiro has no BOUNDS section).
fn afiro_with_bounds(bounds: &str) -> String {
    let text = std::fs::read_to_string(afiro_path()).unwrap();
    assert!(!text.contains("BOUNDS"));
    text.replace("ENDATA", &format!("BOUNDS\n{bounds}ENDATA"))
}

// ------------------------------------------------------------------ client

struct Client {
    child: Child,
    stdin: ChildStdin,
    out: BufReader<ChildStdout>,
    next: i64,
}

impl Client {
    fn spawn() -> Self {
        let mut child = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
            .arg("mcp")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .unwrap();
        let stdin = child.stdin.take().unwrap();
        let out = BufReader::new(child.stdout.take().unwrap());
        Client {
            child,
            stdin,
            out,
            next: 1,
        }
    }

    fn send_raw(&mut self, line: &str) -> Value {
        writeln!(self.stdin, "{line}").unwrap();
        self.stdin.flush().unwrap();
        let mut reply = String::new();
        self.out.read_line(&mut reply).unwrap();
        assert!(reply.ends_with('\n'), "no reply to {line}");
        serde_json::from_str(&reply)
            .unwrap_or_else(|e| panic!("stdout line is not JSON ({e}): {reply}"))
    }

    /// Writes one request line without reading the reply.
    fn send_request(&mut self, method: &str, params: Value) -> i64 {
        let id = self.next;
        self.next += 1;
        writeln!(
            self.stdin,
            "{}",
            json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params})
        )
        .unwrap();
        self.stdin.flush().unwrap();
        id
    }

    fn send_notification(&mut self, method: &str, params: Value) {
        writeln!(
            self.stdin,
            "{}",
            json!({"jsonrpc": "2.0", "method": method, "params": params})
        )
        .unwrap();
        self.stdin.flush().unwrap();
    }

    fn read_line(&mut self) -> Value {
        let mut line = String::new();
        self.out.read_line(&mut line).unwrap();
        assert!(line.ends_with('\n'), "server closed stdout");
        serde_json::from_str(&line).unwrap_or_else(|e| panic!("not JSON ({e}): {line}"))
    }

    fn request(&mut self, method: &str, params: Value) -> Value {
        let id = self.next;
        self.next += 1;
        let reply = self.send_raw(
            &json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params}).to_string(),
        );
        assert_eq!(reply["id"], json!(id));
        assert_eq!(reply["jsonrpc"], json!("2.0"));
        reply
    }

    fn handshake(&mut self) {
        let r = self.request(
            "initialize",
            json!({"protocolVersion": "2025-11-25", "capabilities": {}, "clientInfo": {"name": "test", "version": "0"}}),
        );
        assert_eq!(r["result"]["protocolVersion"], json!("2025-11-25"));
        // a notification gets no reply: the next request's reply is the next line
        writeln!(
            self.stdin,
            "{}",
            json!({"jsonrpc": "2.0", "method": "notifications/initialized"})
        )
        .unwrap();
    }

    /// `tools/call`: the raw result object (`isError` may be true).
    fn call(&mut self, tool: &str, args: Value) -> Value {
        let r = self.request("tools/call", json!({"name": tool, "arguments": args}));
        assert!(r.get("error").is_none(), "protocol error for {tool}: {r}");
        r["result"].clone()
    }

    /// A successful call's `structuredContent`, after checking the text twin and
    /// the output schema.
    fn ok(&mut self, tool: &str, args: Value) -> Value {
        let res = self.call(tool, args);
        assert_eq!(res["isError"], json!(false), "{tool} failed: {res}");
        let structured = res["structuredContent"].clone();
        let text: Value =
            serde_json::from_str(res["content"][0]["text"].as_str().unwrap()).unwrap();
        assert_eq!(
            text, structured,
            "text content must repeat structuredContent"
        );
        let schema = self.output_schema(tool);
        let errs = Schema::new(&schema).check(&structured);
        assert!(
            errs.is_empty(),
            "{tool} output violates its outputSchema: {errs:?}\n{structured}"
        );
        structured
    }

    fn output_schema(&mut self, tool: &str) -> Value {
        let list = self.request("tools/list", json!({}));
        list["result"]["tools"]
            .as_array()
            .unwrap()
            .iter()
            .find(|t| t["name"] == json!(tool))
            .unwrap()["outputSchema"]
            .clone()
    }

    /// The raw text of a failed call: the `rusolve.diagnostics/1` document.
    fn err_raw(&mut self, tool: &str, args: Value) -> String {
        let res = self.call(tool, args);
        assert_eq!(res["isError"], json!(true), "{tool} should fail: {res}");
        assert!(res.get("structuredContent").is_none());
        res["content"][0]["text"].as_str().unwrap().to_string()
    }

    /// The message of the first diagnostic of a failed call.
    fn err_text(&mut self, tool: &str, args: Value) -> String {
        let text = self.err_raw(tool, args);
        let d: Value = serde_json::from_str(&text).expect("error text is JSON");
        d["diagnostics"][0]["message"].as_str().unwrap().to_string()
    }

    fn load_text(&mut self, mps: &str) -> String {
        self.ok("load_model", json!({"mps_text": mps}))["model_id"]
            .as_str()
            .unwrap()
            .to_string()
    }
}

impl Drop for Client {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

// ------------------------------------------------------- JSON Schema checker

struct Schema<'a> {
    root: &'a Value,
}

impl<'a> Schema<'a> {
    fn new(root: &'a Value) -> Self {
        Schema { root }
    }

    fn check(&self, v: &Value) -> Vec<String> {
        let mut errs = Vec::new();
        self.walk(self.root, v, "$", &mut errs);
        errs
    }

    fn resolve(&self, r: &str) -> &'a Value {
        let mut cur = self.root;
        for part in r.trim_start_matches("#/").split('/') {
            cur = &cur[part];
        }
        assert!(!cur.is_null(), "unresolved $ref {r}");
        cur
    }

    fn type_ok(t: &str, v: &Value) -> bool {
        match t {
            "object" => v.is_object(),
            "array" => v.is_array(),
            "string" => v.is_string(),
            "boolean" => v.is_boolean(),
            "null" => v.is_null(),
            "number" => v.is_number(),
            "integer" => v.is_i64() || v.is_u64() || v.as_f64().is_some_and(|x| x.fract() == 0.0),
            other => panic!("unknown JSON Schema type {other}"),
        }
    }

    fn walk(&self, s: &Value, v: &Value, at: &str, errs: &mut Vec<String>) {
        if let Some(r) = s.get("$ref").and_then(Value::as_str) {
            return self.walk(self.resolve(r), v, at, errs);
        }
        if let Some(t) = s.get("type") {
            let ok = match t {
                Value::String(t) => Self::type_ok(t, v),
                Value::Array(ts) => ts.iter().any(|t| Self::type_ok(t.as_str().unwrap(), v)),
                _ => panic!("bad type keyword"),
            };
            if !ok {
                errs.push(format!("{at}: expected {t}, got {v}"));
                return;
            }
        }
        if let Some(c) = s.get("const")
            && c != v
        {
            errs.push(format!("{at}: expected const {c}, got {v}"));
        }
        if let Some(e) = s.get("enum").and_then(Value::as_array)
            && !e.contains(v)
        {
            errs.push(format!("{at}: {v} not in enum {e:?}"));
        }
        if let (Some(m), Some(x)) = (s.get("minimum").and_then(Value::as_f64), v.as_f64())
            && x < m
        {
            errs.push(format!("{at}: {x} < minimum {m}"));
        }
        for key in ["oneOf", "anyOf"] {
            if let Some(alts) = s.get(key).and_then(Value::as_array) {
                let passing = alts
                    .iter()
                    .filter(|a| {
                        let mut e = Vec::new();
                        self.walk(a, v, at, &mut e);
                        e.is_empty()
                    })
                    .count();
                if passing == 0 || (key == "oneOf" && passing > 1) {
                    errs.push(format!(
                        "{at}: {passing} of the {key} alternatives match {v}"
                    ));
                }
            }
        }
        if let Some(obj) = v.as_object() {
            for r in s
                .get("required")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
            {
                if !obj.contains_key(r.as_str().unwrap()) {
                    errs.push(format!("{at}: missing required {r}"));
                }
            }
            let props = s.get("properties").and_then(Value::as_object);
            for (k, x) in obj {
                let sub = format!("{at}.{k}");
                match (props.and_then(|p| p.get(k)), s.get("additionalProperties")) {
                    (Some(ps), _) => self.walk(ps, x, &sub, errs),
                    (None, Some(Value::Bool(false))) => {
                        errs.push(format!("{sub}: unexpected property"))
                    }
                    (None, Some(a)) if a.is_object() => self.walk(a, x, &sub, errs),
                    _ => {}
                }
            }
        }
        if let (Some(arr), Some(item)) = (v.as_array(), s.get("items")) {
            for (i, x) in arr.iter().enumerate() {
                self.walk(item, x, &format!("{at}[{i}]"), errs);
            }
        }
    }
}

/// A schema document is itself well formed: every keyword we rely on has the
/// right shape, recursively.
fn assert_schema_wellformed(s: &Value, at: &str) {
    let obj = s
        .as_object()
        .unwrap_or_else(|| panic!("{at}: schema is not an object"));
    if let Some(t) = obj.get("type") {
        let known = [
            "object", "array", "string", "boolean", "null", "number", "integer",
        ];
        let names: Vec<&str> = match t {
            Value::String(t) => vec![t.as_str()],
            Value::Array(ts) => ts.iter().map(|t| t.as_str().unwrap()).collect(),
            _ => panic!("{at}: bad type"),
        };
        for n in names {
            assert!(known.contains(&n), "{at}: unknown type {n}");
        }
    }
    if let Some(req) = obj.get("required") {
        let props = obj.get("properties").and_then(Value::as_object);
        for r in req
            .as_array()
            .unwrap_or_else(|| panic!("{at}: required is not an array"))
        {
            let r = r.as_str().unwrap();
            // $ref/oneOf schemas may declare required names elsewhere
            if let Some(p) = props
                && !obj.contains_key("oneOf")
            {
                assert!(p.contains_key(r), "{at}: required '{r}' is not a property");
            }
        }
    }
    if let Some(p) = obj.get("properties") {
        for (k, v) in p
            .as_object()
            .unwrap_or_else(|| panic!("{at}: properties is not an object"))
        {
            assert_schema_wellformed(v, &format!("{at}.{k}"));
        }
    }
    if let Some(i) = obj.get("items") {
        assert_schema_wellformed(i, &format!("{at}[]"));
    }
    if let Some(a) = obj.get("additionalProperties")
        && a.is_object()
    {
        assert_schema_wellformed(a, &format!("{at}.*"));
    }
    for key in ["oneOf", "anyOf", "allOf"] {
        if let Some(alts) = obj.get(key) {
            for (i, a) in alts.as_array().unwrap().iter().enumerate() {
                assert_schema_wellformed(a, &format!("{at}.{key}[{i}]"));
            }
        }
    }
    for key in ["$defs", "definitions"] {
        if let Some(d) = obj.get(key) {
            for (k, v) in d.as_object().unwrap() {
                assert_schema_wellformed(v, &format!("{at}.{key}.{k}"));
            }
        }
    }
}

// -------------------------------------------------------------------- tests

#[test]
fn mcp_handshake_lists_tools_with_valid_schemas() {
    let mut c = Client::spawn();
    let r = c.request(
        "initialize",
        json!({"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "t", "version": "0"}}),
    );
    assert_eq!(r["result"]["protocolVersion"], json!("2025-06-18"));
    assert_eq!(r["result"]["serverInfo"]["name"], json!("rusolve"));
    assert!(r["result"]["capabilities"]["tools"].is_object());
    // an unknown requested version is answered with one the server supports
    let r = c.request(
        "initialize",
        json!({"protocolVersion": "1999-01-01", "capabilities": {}}),
    );
    assert_eq!(r["result"]["protocolVersion"], json!("2025-11-25"));
    c.handshake();
    assert_eq!(c.request("ping", json!({}))["result"], json!({}));

    let list = c.request("tools/list", json!({}));
    let tools = list["result"]["tools"].as_array().unwrap();
    let names: Vec<&str> = tools.iter().map(|t| t["name"].as_str().unwrap()).collect();
    assert_eq!(
        names,
        [
            "load_model",
            "build_model",
            "solve",
            "get_solution",
            "modify",
            "check_solution",
            "explain",
            "lint",
            "ranging",
            "feasrelax",
            "free_model"
        ]
    );
    for t in tools {
        let name = t["name"].as_str().unwrap();
        assert!(
            t["description"].as_str().is_some_and(|d| d.len() > 20),
            "{name} description"
        );
        for key in ["inputSchema", "outputSchema"] {
            assert_eq!(t[key]["type"], json!("object"), "{name}.{key}");
            assert_schema_wellformed(&t[key], &format!("{name}.{key}"));
        }
        assert_eq!(
            t["inputSchema"]["additionalProperties"],
            json!(false),
            "{name}"
        );
    }
    // deterministic order
    assert_eq!(c.request("tools/list", json!({}))["result"], list["result"]);

    // the 2026-07-28 era: per-request _meta, server/discover, resultType
    let meta = json!({"io.modelcontextprotocol/protocolVersion": "2026-07-28",
                      "io.modelcontextprotocol/clientInfo": {"name": "t", "version": "0"},
                      "io.modelcontextprotocol/clientCapabilities": {}});
    let d = c.request("server/discover", json!({"_meta": meta}));
    assert_eq!(d["result"]["resultType"], json!("complete"));
    assert_eq!(d["result"]["supportedVersions"][0], json!("2026-07-28"));
    assert_eq!(
        d["result"]["_meta"]["io.modelcontextprotocol/serverInfo"]["name"],
        json!("rusolve")
    );
    let l = c.request("tools/list", json!({"_meta": meta}));
    assert_eq!(l["result"]["resultType"], json!("complete"));
    assert_eq!(l["result"]["tools"].as_array().unwrap().len(), 11);
}

#[test]
fn mcp_solve_roundtrip() {
    let mut c = Client::spawn();
    c.handshake();
    let loaded = c.ok("load_model", json!({"mps_text": TINY_MIP}));
    assert_eq!(loaded["rows"], json!(2));
    assert_eq!(loaded["cols"], json!(4));
    assert_eq!(loaded["integer_cols"], json!(3));
    assert_eq!(
        loaded["names"]["cols"]["first"],
        json!(["x", "y", "z", "w"])
    );
    let id = loaded["model_id"].as_str().unwrap();

    let env = c.ok("solve", json!({"model_id": id, "include_solution": true}));
    assert_eq!(env["schema_version"], json!("rusolve.result/1"));
    assert_eq!(env["status"], json!("Optimal"));
    assert!((env["objective"].as_f64().unwrap() + 20.0).abs() < 1e-9);
    assert_eq!(env["input"]["integer_cols"], json!(3));
    assert_eq!(
        env["solution"]["columns"][0],
        json!({"name": "x", "value": 4.0})
    );

    let sol = c.ok("get_solution", json!({"model_id": id, "names": ["y", "x"]}));
    assert_eq!(
        sol["columns"],
        json!([{"name": "y", "value": 2.0}, {"name": "x", "value": 4.0}])
    );

    // an option and a time limit go through the ported setters
    let env = c.ok(
        "solve",
        json!({"model_id": id, "options": {"mip_rel_gap": 0.5, "presolve": "off"}, "time_limit": 30}),
    );
    assert_eq!(env["run"]["options"]["presolve"], json!("off"));
    assert_eq!(env["status"], json!("Optimal"));

    // a handle is independent of the others; free_model releases it
    let other = c.load_text(INFEASIBLE);
    assert_ne!(other, id);
    let freed = c.ok("free_model", json!({"model_id": id}));
    assert_eq!(freed["remaining"], json!(1));
    assert!(
        c.err_text("solve", json!({"model_id": id}))
            .contains("model_id")
    );
}

/// An LP re-solve after a bound edit starts from the kept basis and takes
/// fewer simplex iterations than a cold solve; the result equals a fresh solve
/// of the edited model (the same model written out as MPS).
#[test]
fn mcp_edit_warm_resolve_takes_fewer_iterations() {
    let mut c = Client::spawn();
    c.handshake();
    let id = c.ok("load_model", json!({"path": afiro_path()}))["model_id"]
        .as_str()
        .unwrap()
        .to_string();
    let first = c.ok("solve", json!({"model_id": id}));
    assert_eq!(first["status"], json!("Optimal"));
    let first_iters = first["stats"]["lp_iterations"].as_i64().unwrap();
    assert!(first_iters > 0);

    // X02 is 25.5 at the optimum; cap it at 20
    let edit = json!([{"kind": "col_bounds", "name": "X02", "upper": 20.0}]);
    let m = c.ok("modify", json!({"model_id": id, "edits": edit}));
    assert_eq!(m["applied"], json!(1));
    assert_eq!(m["basis_kept"], json!(true));
    // no solution until the next solve
    assert!(
        c.err_text("get_solution", json!({"model_id": id}))
            .contains("solve")
    );

    let warm = c.ok("solve", json!({"model_id": id}));
    let warm_iters = warm["stats"]["lp_iterations"].as_i64().unwrap();
    assert_eq!(warm["input"]["edits_applied"], json!(1));
    assert!(warm["input"]["path"].as_str().unwrap().contains("(edited)"));

    // the same edited model, solved cold on the same handle...
    let cold = c.ok("solve", json!({"model_id": id, "warm_start": false}));
    let cold_iters = cold["stats"]["lp_iterations"].as_i64().unwrap();
    // ...and as a fresh model written to MPS
    let fresh_id = c.load_text(&afiro_with_bounds(" UP BND X02 20\n"));
    let fresh = c.ok("solve", json!({"model_id": fresh_id}));
    let fresh_iters = fresh["stats"]["lp_iterations"].as_i64().unwrap();

    eprintln!(
        "afiro: first solve {first_iters} iterations; after the X02 bound edit warm {warm_iters}, cold {cold_iters}, fresh-from-MPS {fresh_iters}"
    );
    assert!(
        warm_iters < cold_iters,
        "warm {warm_iters} >= cold {cold_iters}"
    );
    let (w, f, k) = (
        warm["objective"].as_f64().unwrap(),
        fresh["objective"].as_f64().unwrap(),
        cold["objective"].as_f64().unwrap(),
    );
    assert!(
        (w - f).abs() <= 1e-9 * f.abs().max(1.0),
        "warm {w} vs fresh {f}"
    );
    assert!(
        (k - f).abs() <= 1e-9 * f.abs().max(1.0),
        "cold {k} vs fresh {f}"
    );
    assert!(
        (w - first["objective"].as_f64().unwrap()).abs() > 1e-6,
        "the edit must bite"
    );
    let all_w = c.ok("get_solution", json!({"model_id": id}));
    let all_f = c.ok("get_solution", json!({"model_id": fresh_id}));
    assert_eq!(all_w["columns"].as_array().unwrap().len(), 32);
    for (a, b) in all_w["columns"]
        .as_array()
        .unwrap()
        .iter()
        .zip(all_f["columns"].as_array().unwrap())
    {
        assert_eq!(a["name"], b["name"]);
    }
}

/// Every edit kind, each checked against a fresh solve of the equivalent MPS.
#[test]
fn mcp_modify_every_kind_matches_fresh_models() {
    let mut c = Client::spawn();
    c.handshake();
    let objective = |c: &mut Client, id: &str| {
        let e = c.ok("solve", json!({"model_id": id}));
        assert_eq!(e["status"], json!("Optimal"), "{e}");
        e["objective"].as_f64().unwrap()
    };
    let base = "\
NAME BASE
ROWS
 N  obj
 L  cap
 L  mix
COLUMNS
    a  obj  -3  cap  1
    a  mix  1
    b  obj  -2  cap  1
    b  mix  -1
RHS
    rhs  cap  10  mix  4
BOUNDS
 UP bnd a 8
 UP bnd b 8
ENDATA
";
    let id = c.load_text(base);
    assert!((objective(&mut c, &id) + 27.0).abs() < 1e-9);
    // cost, row bound
    c.ok(
        "modify",
        json!({"model_id": id, "edits": [
        {"kind": "cost", "name": "b", "cost": -5.0},
        {"kind": "row_bounds", "name": "cap", "upper": 12.0}]}),
    );
    let fresh = c.load_text(
        &base
            .replace("b  obj  -2", "b  obj  -5")
            .replace("cap  10", "cap  12"),
    );
    let (x, y) = (objective(&mut c, &id), objective(&mut c, &fresh));
    assert!((x - y).abs() < 1e-9, "{x} vs {y}");
    // coef + add_row + add_col + del_col... one at a time against fresh models
    c.ok("modify", json!({"model_id": id, "edits": [{"kind": "coef", "row": "mix", "col": "b", "value": -2.0}]}));
    let fresh = c.load_text(
        &base
            .replace("b  obj  -2", "b  obj  -5")
            .replace("cap  10", "cap  12")
            .replace("b  mix  -1", "b  mix  -2"),
    );
    let (x, y) = (objective(&mut c, &id), objective(&mut c, &fresh));
    assert!((x - y).abs() < 1e-9, "{x} vs {y}");

    c.ok(
        "modify",
        json!({"model_id": id, "edits": [
        {"kind": "add_col", "name": "d", "cost": -1.0, "upper": 3.0, "coefs": {"cap": 1.0}},
        {"kind": "add_row", "name": "lim", "upper": 7.0, "coefs": {"a": 1.0, "d": 1.0}}]}),
    );
    let grown = "\
NAME GROWN
ROWS
 N  obj
 L  cap
 L  mix
 L  lim
COLUMNS
    a  obj  -3  cap  1
    a  mix  1  lim  1
    b  obj  -5  cap  1
    b  mix  -2
    d  obj  -1  cap  1
    d  lim  1
RHS
    rhs  cap  12  mix  4
    rhs  lim  7
BOUNDS
 UP bnd a 8
 UP bnd b 8
 UP bnd d 3
ENDATA
";
    let fresh = c.load_text(grown);
    let (x, y) = (objective(&mut c, &id), objective(&mut c, &fresh));
    assert!((x - y).abs() < 1e-9, "{x} vs {y}");
    let names = c.ok("get_solution", json!({"model_id": id, "names": ["d", "a"]}));
    assert_eq!(names["columns"][0]["name"], json!("d"));

    c.ok(
        "modify",
        json!({"model_id": id, "edits": [
        {"kind": "del_row", "name": "mix"}, {"kind": "del_col", "name": "d"}]}),
    );
    let shrunk = "\
NAME SHRUNK
ROWS
 N  obj
 L  cap
 L  lim
COLUMNS
    a  obj  -3  cap  1
    a  lim  1
    b  obj  -5  cap  1
RHS
    rhs  cap  12  lim  7
BOUNDS
 UP bnd a 8
 UP bnd b 8
ENDATA
";
    let fresh = c.load_text(shrunk);
    let (x, y) = (objective(&mut c, &id), objective(&mut c, &fresh));
    assert!((x - y).abs() < 1e-9, "{x} vs {y}");

    // errors name the edit and the field, and say what was applied
    let e = c.err_text("modify", json!({"model_id": id, "edits": [
        {"kind": "cost", "name": "a", "cost": 1.0}, {"kind": "cost", "name": "nope", "cost": 1.0}]}));
    assert!(
        e.contains("edits[1].name") && e.contains("nope") && e.contains("1 edit(s)"),
        "{e}"
    );
    let e = c.err_text(
        "modify",
        json!({"model_id": id, "edits": [{"kind": "col_bounds", "name": "a", "lower": "x"}]}),
    );
    assert!(e.contains("edits[0].lower"), "{e}");
    let e = c.err_text(
        "modify",
        json!({"model_id": id, "edits": [{"kind": "wobble"}]}),
    );
    assert!(e.contains("edits[0].kind"), "{e}");
}

#[test]
fn mcp_infeasible_model_through_explain() {
    let mut c = Client::spawn();
    c.handshake();
    let id = c.load_text(INFEASIBLE);
    let env = c.ok("solve", json!({"model_id": id}));
    assert_eq!(env["status"], json!("Infeasible"));
    let ex = c.ok("explain", json!({"model_id": id}));
    assert_eq!(ex["schema_version"], json!("rusolve.explain/1"));
    assert_eq!(ex["status"], json!("Infeasible"));
    assert_eq!(ex["explanation"]["type"], json!("infeasible"));
    assert_eq!(ex["verified"], json!(true), "{ex}");
    let text = ex.to_string();
    assert!(
        text.contains("r1") && text.contains("r2"),
        "names in the certificate: {text}"
    );

    // the cheapest fix, then apply it through modify and re-solve
    let fr = c.ok("feasrelax", json!({"model_id": id}));
    assert_eq!(fr["status"], json!("Relaxed"));
    assert_eq!(fr["verified"], json!(true));
    let edit = &fr["edits"][0];
    assert_eq!(edit["kind"], json!("row"));
    let lint = c.ok("lint", json!({"model_id": id}));
    assert_eq!(lint["schema_version"], json!("rusolve.lint/1"));
    let to = edit["to"].as_f64().unwrap();
    let name = edit["name"].as_str().unwrap();
    let field = if edit["bound"] == json!("lower") {
        "lower"
    } else {
        "upper"
    };
    c.ok(
        "modify",
        json!({"model_id": id, "edits": [{"kind": "row_bounds", "name": name, field: to}]}),
    );
    assert_eq!(
        c.ok("solve", json!({"model_id": id}))["status"],
        json!("Optimal")
    );
}

#[test]
fn mcp_check_lint_ranging_on_an_optimal_lp() {
    let mut c = Client::spawn();
    c.handshake();
    let id = c.load_text(TINY_MIP);
    let env = c.ok("solve", json!({"model_id": id, "include_solution": true}));
    // the solve result's solution is accepted as is
    let ok = c.ok(
        "check_solution",
        json!({"model_id": id, "solution": env["solution"]}),
    );
    assert_eq!(ok["schema_version"], json!("rusolve.check/1"));
    assert_eq!(ok["verdict"], json!("feasible"), "{ok}");
    let bad = c.ok(
        "check_solution",
        json!({"model_id": id, "solution": {"values": {"x": 4.0, "y": 5.0, "z": 0.5, "w": 0.0}}}),
    );
    assert_ne!(bad["verdict"], json!("feasible"), "{bad}");
    assert!(bad["findings"].to_string().contains("cap1"));
    let e = c.err_text(
        "check_solution",
        json!({"model_id": id, "solution": {"values": {"x": 1.0}}}),
    );
    assert!(e.contains("solution") && e.contains("missing"), "{e}");

    let lint = c.ok(
        "lint",
        json!({"model_id": id, "rules": ["big-m"], "min_severity": "info"}),
    );
    assert_eq!(lint["rules_run"], json!(["big-m"]));
    let e = c.err_text("lint", json!({"model_id": id, "rules": ["nope"]}));
    assert!(e.starts_with("rules:"), "{e}");

    // ranging by name on the fixed-integer LP of the MIP
    let rg = c.ok(
        "ranging",
        json!({"model_id": id, "cols": ["w"], "rows": ["cap*"]}),
    );
    assert_eq!(rg["columns"][0]["name"], json!("w"));
    assert_eq!(rg["basis_of"], json!("fixed_integer_lp"));
    let e = c.err_text("ranging", json!({"model_id": id, "cols": ["zz*"]}));
    assert!(e.contains("cols"), "{e}");
}

/// A failed tool call carries the `rusolve.diagnostics/1` document
/// (docs/design/diagnostics.md) as its text, with a registered code.
#[test]
fn mcp_tool_errors_carry_diagnostics_json() {
    let mut c = Client::spawn();
    c.handshake();
    let code_of = |text: &str| -> Value {
        let d: Value = serde_json::from_str(text).expect("error text is JSON");
        assert_eq!(d["schema_version"], json!("rusolve.diagnostics/1"), "{d}");
        d["diagnostics"][0].clone()
    };
    let d = code_of(&c.err_raw("load_model", json!({"mps_text": "not an mps file"})));
    assert_eq!(d["code"], json!("RS0011"), "{d}");
    assert_eq!(d["severity"], json!("error"), "{d}");
    let d = code_of(&c.err_raw("solve", json!({"model_id": "m99"})));
    assert_eq!(d["code"], json!("RS3002"), "{d}");
    let d = code_of(&c.err_raw("solve", json!({})));
    assert_eq!(d["code"], json!("RS3001"), "{d}");
}

#[test]
fn mcp_bad_tool_argument_is_an_error_result_naming_the_field() {
    let mut c = Client::spawn();
    c.handshake();
    let id = c.load_text(TINY_MIP);
    let e = c.err_text("solve", json!({"model_id": id, "time_limit": "soon"}));
    assert!(e.contains("time_limit"), "{e}");
    let e = c.err_text(
        "solve",
        json!({"model_id": id, "options": {"no_such_option": 1}}),
    );
    assert!(e.contains("options.no_such_option"), "{e}");
    let e = c.err_text("solve", json!({"model_id": id, "threads": 4}));
    assert!(
        e.contains("threads") && e.contains("unknown argument"),
        "{e}"
    );
    let e = c.err_text("solve", json!({"model_id": id, "wasm_fuel": 10}));
    assert!(e.contains("wasm_fuel") && e.contains("wasm_hook"), "{e}");
    let e = c.err_text("solve", json!({"model_id": id, "wasm_memory": 4096}));
    assert!(e.contains("wasm_memory") && e.contains("wasm_hook"), "{e}");
    let e = c.err_text("solve", json!({}));
    assert!(e.contains("model_id"), "{e}");
    let e = c.err_text("solve", json!({"model_id": "m99"}));
    assert!(e.contains("model_id") && e.contains("m99"), "{e}");
    let e = c.err_text("load_model", json!({"path": "/no/such/model.mps"}));
    assert!(e.contains("path"), "{e}");
    let e = c.err_text("load_model", json!({"path": "a", "mps_text": "b"}));
    assert!(e.contains("path / mps_text"), "{e}");
    let e = c.err_text("load_model", json!({"mps_text": "not an mps file"}));
    assert!(e.contains("mps_text"), "{e}");
    let e = c.err_text("get_solution", json!({"model_id": id}));
    assert!(e.contains("solve first") || e.contains("call solve"), "{e}");
    c.ok("solve", json!({"model_id": id}));
    let e = c.err_text("get_solution", json!({"model_id": id, "names": ["x", "q"]}));
    assert!(e.contains("names[1]") && e.contains("'q'"), "{e}");
    let e = c.err_text("feasrelax", json!({"model_id": id, "mode": "fast"}));
    assert!(e.contains("mode"), "{e}");
    // the server is still alive and the handle still works
    assert_eq!(
        c.ok("solve", json!({"model_id": id}))["status"],
        json!("Optimal")
    );
}

#[test]
fn mcp_protocol_errors_use_json_rpc_codes() {
    let mut c = Client::spawn();
    c.handshake();
    let r = c.request("no/such/method", json!({}));
    assert_eq!(r["error"]["code"], json!(-32601));
    let r = c.request(
        "tools/call",
        json!({"name": "no_such_tool", "arguments": {}}),
    );
    assert_eq!(r["error"]["code"], json!(-32602));
    assert!(
        r["error"]["message"]
            .as_str()
            .unwrap()
            .contains("no_such_tool")
    );
    let r = c.request("tools/call", json!({"arguments": {}}));
    assert_eq!(r["error"]["code"], json!(-32602));
    let r = c.request("tools/call", json!({"name": "solve", "arguments": [1]}));
    assert_eq!(r["error"]["code"], json!(-32602));
    // malformed JSON: parse error with a null id, and the server keeps going
    let r = c.send_raw("{not json");
    assert_eq!(r["error"]["code"], json!(-32700));
    assert_eq!(r["id"], Value::Null);
    let r = c.send_raw("[1,2]");
    assert!(r.is_array() || r.get("error").is_some());
    let r = c.send_raw("\"just a string\"");
    assert_eq!(r["error"]["code"], json!(-32600));
    // an unsupported per-request protocol version
    let r = c.request(
        "tools/list",
        json!({"_meta": {"io.modelcontextprotocol/protocolVersion": "1900-01-01"}}),
    );
    assert_eq!(r["error"]["code"], json!(-32022));
    assert_eq!(r["error"]["data"]["requested"], json!("1900-01-01"));
    assert!(
        r["error"]["data"]["supported"]
            .as_array()
            .unwrap()
            .contains(&json!("2026-07-28"))
    );
    assert_eq!(c.request("ping", json!({}))["result"], json!({}));
}

/// One server process solves the same handle with 1, 4, 2 and 1 threads.
/// Without `reset_scheduler_for_threads` the second solve failed with "global scheduler has
/// already been initialized to use 1 threads".
#[test]
fn mcp_solve_can_change_threads_between_solves() {
    let mut c = Client::spawn();
    c.handshake();
    let id = c.load_text(TINY_MIP);
    let mut objective = None;
    for threads in [1, 4, 2, 1] {
        let env = c.ok(
            "solve",
            json!({"model_id": id, "options": {"threads": threads}}),
        );
        assert_eq!(env["status"], json!("Optimal"), "threads {threads}: {env}");
        assert_eq!(env["run"]["threads"], json!(threads), "{env}");
        let obj = env["objective"].as_f64().unwrap();
        assert_eq!(*objective.get_or_insert(obj), obj, "threads {threads}");
    }
}

/// An LP that takes the automatic HiPO pick (at least 3 nonzeros per column, fewer than 3 rows
/// per column), solved with threads 4 and then 1 in one server: `Highs::run` reads the worker
/// deque for that pick before it initialises multi-threading, so the deque of the scheduler
/// that the thread change shut down must not be read.
#[test]
fn mcp_lp_threads_4_then_1_reads_a_live_scheduler() {
    let (rows, cols) = (30usize, 40usize);
    let mut mps = String::from("NAME LPCOV\nROWS\n N  cost\n");
    for r in 0..rows {
        mps += &format!(" G  r{r}\n");
    }
    mps += "COLUMNS\n";
    for j in 0..cols {
        mps += &format!("    c{j}  cost  {}", 1 + j % 7);
        for k in 0..4 {
            mps += &format!("  r{}  1", (j + k) % rows);
        }
        mps += "\n";
    }
    mps += "RHS\n";
    for r in 0..rows {
        mps += &format!("    rhs  r{r}  1\n");
    }
    mps += "ENDATA\n";
    let mut c = Client::spawn();
    c.handshake();
    let id = c.load_text(&mps);
    let mut objective = None;
    for threads in [4, 1, 4, 1, 2, 1] {
        let env = c.ok(
            "solve",
            json!({"model_id": id, "options": {"threads": threads}}),
        );
        assert_eq!(env["status"], json!("Optimal"), "threads {threads}: {env}");
        let obj = env["objective"].as_f64().unwrap();
        assert!(
            (*objective.get_or_insert(obj) - obj).abs() < 1e-6,
            "threads {threads}"
        );
    }
}

/// `solve` with `wasm_hook` loads the sandboxed Wasm component, its counters come back as
/// `stats.wasm`, and the objective equals the hook-less solve (a component can only reorder).
/// A bad component is a tool error naming the file; the server stays alive.
#[test]
fn mcp_solve_with_wasm_hook_reports_stats() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let model = dir
        .join("crates/rusolve-ekk/tests/fixtures/s5/e2e_knap.mps")
        .to_string_lossy()
        .into_owned();
    let hook = dir
        .join("crates/rusolve-wasm/components/most_fractional.wat")
        .to_string_lossy()
        .into_owned();
    let mut c = Client::spawn();
    c.handshake();
    let id = c.ok("load_model", json!({"path": model}))["model_id"]
        .as_str()
        .unwrap()
        .to_string();
    let opts = json!({"presolve": "off"});
    let with = c.ok(
        "solve",
        json!({"model_id": id, "options": opts, "wasm_hook": hook, "wasm_fuel": 1000000}),
    );
    assert_eq!(with["status"], json!("Optimal"));
    let w = &with["stats"]["wasm"];
    assert_eq!(w["component"], json!("most_fractional"));
    assert!(w["calls"].as_u64().unwrap() > 0, "{w}");
    assert_eq!(
        w["ok"], w["calls"],
        "the component answered every call: {w}"
    );
    assert_eq!(w["traps"].as_u64().unwrap(), 0);
    assert_eq!(w["out_of_fuel"].as_u64().unwrap(), 0);
    // the same solve without the hook: same objective, no stats.wasm
    let plain = c.ok("solve", json!({"model_id": id, "options": opts}));
    assert_eq!(plain["objective"], with["objective"]);
    assert!(plain["stats"].get("wasm").is_none());
    // `wasm_memory` is the linear-memory cap, the same knob as `--wasm-memory`. A generous
    // cap changes nothing; a cap below the component's declared memory is a load failure that
    // names the file, exactly like the CLI's exit 3.
    let roomy = c.ok(
        "solve",
        json!({"model_id": id, "options": opts, "wasm_hook": hook, "wasm_memory": 1048576}),
    );
    assert_eq!(roomy["objective"], plain["objective"]);
    assert_eq!(
        roomy["stats"]["wasm"]["component"],
        json!("most_fractional")
    );
    assert!(roomy["stats"]["wasm"]["calls"].as_u64().unwrap() > 0);
    let e = c.err_text(
        "solve",
        json!({"model_id": id, "wasm_hook": hook, "wasm_memory": 4096}),
    );
    assert!(e.contains("most_fractional.wat"), "{e}");
    // a missing component is a tool error naming the file, not a crash
    let e = c.err_text(
        "solve",
        json!({"model_id": id, "wasm_hook": "/no/such/component.wat"}),
    );
    assert!(e.contains("component.wat"), "{e}");
    assert_eq!(
        c.ok("solve", json!({"model_id": id}))["status"],
        json!("Optimal")
    );
}

// ------------------------------------------------------------------ MCP reach

const LP_TEXT: &str = "Maximize\n obj: 3 x + 2 y\nSubject To\n c1: x + y <= 4\n c2: x + 3 y <= 6\nBounds\n 0 <= x <= 3\nEnd\n";

/// min -x - 2y + (x^2 + y^2)/2, x + y <= 1, 0 <= x, y <= 10. Optimum (0, 1), objective -1.5.
const QP_TEXT: &str = "Minimize\n obj: - x - 2 y + [ x ^ 2 + y ^ 2 ] / 2\nSubject To\n c1: x + y <= 1\nBounds\n 0 <= x <= 10\n 0 <= y <= 10\nEnd\n";

const NONCONVEX_TEXT: &str = "Minimize\n obj: - x - 2 y + [ x ^ 2 + 4 x * y + y ^ 2 ] / 2\nSubject To\n c1: x + y <= 1\nBounds\n 0 <= x <= 10\n 0 <= y <= 10\nEnd\n";

/// Maximize 3x + 2y - 10b, x + y <= 10, b = 0 -> x <= 2, m = max(x, y) <= 7.
const GC_TEXT: &str = "Maximize\n obj: 3 x + 2 y - 10 b\nSubject To\n c1: x + y <= 10\n ind1: b = 0 -> x <= 2\nGeneral Constraints\n g1: m = MAX ( x , y )\nBounds\n 0 <= x <= 8\n 0 <= y <= 8\n m <= 7\nBinary\n b\nEnd\n";

const MOF_TEXT: &str = r#"{"version": {"major": 1, "minor": 7},
 "variables": [{"name": "x"}, {"name": "y"}],
 "objective": {"sense": "max", "function": {"type": "ScalarAffineFunction", "terms": [
   {"coefficient": 3, "variable": "x"}, {"coefficient": 2, "variable": "y"}], "constant": 0}},
 "constraints": [
  {"name": "c1", "function": {"type": "ScalarAffineFunction", "terms": [
    {"coefficient": 1, "variable": "x"}, {"coefficient": 1, "variable": "y"}], "constant": 0},
   "set": {"type": "LessThan", "upper": 4}},
  {"name": "c2", "function": {"type": "ScalarAffineFunction", "terms": [
    {"coefficient": 1, "variable": "x"}, {"coefficient": 3, "variable": "y"}], "constant": 0},
   "set": {"type": "LessThan", "upper": 6}},
  {"name": "bx", "function": {"type": "ScalarAffineFunction", "terms": [
    {"coefficient": 1, "variable": "x"}], "constant": 0},
   "set": {"type": "Interval", "lower": 0, "upper": 3}}]}"#;

fn value_of(sol: &Value, name: &str) -> f64 {
    sol["columns"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["name"] == json!(name))
        .unwrap_or_else(|| panic!("no column {name}"))["value"]
        .as_f64()
        .unwrap()
}

#[test]
fn mcp_lp_and_mof_text_load_and_solve() {
    let mut c = Client::spawn();
    // optimum: x = 3, y = 1, objective 11
    for (field, text) in [("lp_text", LP_TEXT), ("mof_text", MOF_TEXT)] {
        let loaded = c.ok("load_model", json!({field: text}));
        assert_eq!(
            loaded["rows"],
            json!(if field == "lp_text" { 2 } else { 3 }),
            "{field}"
        );
        assert_eq!(loaded["cols"], json!(2), "{field}");
        assert_eq!(loaded["sense"], json!("maximize"), "{field}");
        assert_eq!(loaded["class"], json!("lp"), "{field}");
        let id = loaded["model_id"].as_str().unwrap().to_string();
        let r = c.ok("solve", json!({"model_id": id}));
        assert_eq!(r["status"], json!("Optimal"), "{field}");
        assert!(
            (r["objective"].as_f64().unwrap() - 11.0).abs() < 1e-7,
            "{field}: {r}"
        );
        let sol = c.ok("get_solution", json!({"model_id": id}));
        assert!((value_of(&sol, "x") - 3.0).abs() < 1e-7, "{field}");
    }
    // two texts, or a text and a path, is an error naming the fields
    let e = c.err_text(
        "load_model",
        json!({"lp_text": LP_TEXT, "mps_text": TINY_MIP}),
    );
    assert!(e.contains("lp_text") && e.contains("mps_text"), "{e}");
    // a broken LP text is an error that names the field, not the temp file
    let e = c.err_text(
        "load_model",
        json!({"lp_text": "Maximize\n obj: 3 x +\nSubject To\n c1: x <= \nEnd\n"}),
    );
    assert!(
        e.starts_with("lp_text:") && !e.contains("rusolve-mcp"),
        "{e}"
    );
}

#[test]
fn mcp_build_model_from_json_matches_the_lp_text() {
    let mut c = Client::spawn();
    let built = c.ok(
        "build_model",
        json!({
            "sense": "maximize",
            "variables": [
                {"name": "x", "upper": 3, "objective": 3},
                {"name": "y"}
            ],
            "objective": {"terms": {"y": 2}, "offset": 5},
            "rows": [
                {"name": "c1", "terms": {"x": 1, "y": 1}, "sense": "<=", "rhs": 4},
                {"name": "c2", "terms": [{"var": "x", "coef": 1}, {"var": "y", "coef": 3}], "upper": 6}
            ]
        }),
    );
    assert_eq!(built["rows"], json!(2));
    assert_eq!(built["cols"], json!(2));
    let id = built["model_id"].as_str().unwrap().to_string();
    let r = c.ok("solve", json!({"model_id": id}));
    assert_eq!(r["status"], json!("Optimal"));
    // the LP text's optimum 11 plus the offset 5
    assert!(
        (r["objective"].as_f64().unwrap() - 16.0).abs() < 1e-7,
        "{r}"
    );

    // integers, binaries, a range row and a free variable: x binary, n integer in [0, 10], f free
    // max 5x + 4n + f, f <= 3, 1 <= x + n + f <= 6 (range), n <= 2.5 -> n = 2, x = 1, f = 3 gives 5 + 8 + 3 with x + n + f = 6
    let m = c.ok(
        "build_model",
        json!({
            "sense": "maximize",
            "variables": [
                {"name": "x", "type": "binary", "objective": 5},
                {"name": "n", "type": "integer", "upper": 2.5, "objective": 4},
                {"name": "f", "lower": null, "upper": 3, "objective": 1}
            ],
            "rows": [{"name": "rg", "terms": {"x": 1, "n": 1, "f": 1}, "lower": 1, "upper": 6}]
        }),
    );
    assert_eq!(m["integer_cols"], json!(2));
    let id = m["model_id"].as_str().unwrap().to_string();
    let r = c.ok("solve", json!({"model_id": id}));
    assert_eq!(r["status"], json!("Optimal"), "{r}");
    assert!(
        (r["objective"].as_f64().unwrap() - 16.0).abs() < 1e-7,
        "{r}"
    );

    // errors carry the JSON path of the offending entry
    let e = c.err_text(
        "build_model",
        json!({"variables": [{"name": "x"}], "rows": [{"name": "r", "terms": {"z": 1}, "upper": 1}]}),
    );
    assert!(e.starts_with("rows[0].terms.z:"), "{e}");
    let e = c.err_text(
        "build_model",
        json!({"variables": [{"name": "x", "type": "complex"}]}),
    );
    assert!(e.starts_with("variables[0].type:"), "{e}");
    let e = c.err_text(
        "build_model",
        json!({"variables": [{"name": "x", "lower": 2, "upper": 1}]}),
    );
    assert!(e.contains("lower 2 > upper 1"), "{e}");
}

#[test]
fn mcp_convex_qp_and_general_constraints_solve_through_the_router() {
    let mut c = Client::spawn();
    let q = c.ok("load_model", json!({"lp_text": QP_TEXT}));
    assert_eq!(q["class"], json!("convex_qp"));
    let qid = q["model_id"].as_str().unwrap().to_string();
    let r = c.ok("solve", json!({"model_id": qid}));
    assert_eq!(r["status"], json!("Optimal"), "{r}");
    assert_eq!(r["class"], json!("convex_qp"));
    assert!((r["objective"].as_f64().unwrap() + 1.5).abs() < 1e-6, "{r}");
    assert!(r.get("solution").is_none(), "no solution unless asked");
    let sol = c.ok("get_solution", json!({"model_id": qid, "names": ["y"]}));
    assert!((value_of(&sol, "y") - 1.0).abs() < 1e-3, "{sol}");
    let with = c.ok("solve", json!({"model_id": qid, "include_solution": true}));
    assert!(with["solution"]["columns"].is_array());
    // the LP/MIP-only tools say so and name the tools that work
    let e = c.err_text("lint", json!({"model_id": qid}));
    assert!(e.contains("convex_qp") && e.contains("solve"), "{e}");
    // an engine option is refused rather than ignored
    let e = c.err_text(
        "solve",
        json!({"model_id": qid, "options": {"presolve": "off"}}),
    );
    assert!(e.starts_with("options:"), "{e}");

    let g = c.ok("load_model", json!({"lp_text": GC_TEXT}));
    assert!(
        ["gc", "gc_mip"].contains(&g["class"].as_str().unwrap()),
        "{g}"
    );
    let gid = g["model_id"].as_str().unwrap().to_string();
    let r = c.ok("solve", json!({"model_id": gid}));
    assert_eq!(r["status"], json!("Optimal"), "{r}");
    // b = 0 forces x <= 2: x = 2, y = 7 (m = max(x, y) <= 7) gives 6 + 14 = 20; b = 1 costs 10
    // and allows x = 3, y = 7: 9 + 14 - 10 = 13
    assert!(
        (r["objective"].as_f64().unwrap() - 20.0).abs() < 1e-6,
        "{r}"
    );
    c.ok("free_model", json!({"model_id": gid}));

    // a class nothing solves is refused at load, with the class and the reason
    let e = c.err_text("load_model", json!({"lp_text": NONCONVEX_TEXT}));
    assert!(
        e.contains("nonconvex_qp") && e.contains("not supported"),
        "{e}"
    );
}

#[test]
fn mcp_prompts_are_the_skills_and_resources_are_the_schemas() {
    let mut c = Client::spawn();
    let skills_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/skills");
    let mut skills: Vec<String> = std::fs::read_dir(&skills_dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().join("SKILL.md").is_file())
        // `audience: contributor` skills are not agent-facing, so not prompts
        .filter(|e| {
            !std::fs::read_to_string(e.path().join("SKILL.md"))
                .unwrap()
                .contains("\naudience: contributor\n")
        })
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    skills.sort();
    assert!(skills.len() >= 4, "{skills:?}");
    assert!(
        skills.iter().all(|s| !s.starts_with("contrib-")),
        "{skills:?}"
    );
    let list = c.request("prompts/list", json!({}));
    let mut names: Vec<String> = list["result"]["prompts"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p["name"].as_str().unwrap().to_string())
        .collect();
    names.sort();
    assert_eq!(names, skills);
    for p in list["result"]["prompts"].as_array().unwrap() {
        assert!(
            p["description"].as_str().is_some_and(|d| d.len() > 10),
            "{p}"
        );
    }
    for s in &skills {
        let g = c.request("prompts/get", json!({"name": s}));
        let m = &g["result"]["messages"][0];
        assert_eq!(m["role"], json!("user"));
        let text = m["content"]["text"].as_str().unwrap();
        let on_disk = std::fs::read_to_string(skills_dir.join(s).join("SKILL.md")).unwrap();
        assert!(!text.starts_with("---"), "front matter stripped");
        assert!(
            on_disk.contains(text.trim()),
            "{s}: body is the file's body"
        );
    }
    let e = c.request("prompts/get", json!({"name": "nope"}));
    assert_eq!(e["error"]["code"], json!(-32602));

    // resources: every schema file, plus the agent reference
    let schema_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/schema");
    let mut stems: Vec<String> = std::fs::read_dir(&schema_dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|f| f.ends_with(".schema.json"))
        .map(|f| f.trim_end_matches(".schema.json").to_string())
        .collect();
    stems.sort();
    let list = c.request("resources/list", json!({}));
    let uris: Vec<String> = list["result"]["resources"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| r["uri"].as_str().unwrap().to_string())
        .collect();
    for stem in &stems {
        assert!(
            uris.contains(&format!("rusolve://schema/{stem}")),
            "{stem} in {uris:?}"
        );
    }
    assert!(uris.contains(&"rusolve://docs/agent/reference".to_string()));
    let r = c.request(
        "resources/read",
        json!({"uri": "rusolve://schema/rusolve-result-1"}),
    );
    let doc: Value =
        serde_json::from_str(r["result"]["contents"][0]["text"].as_str().unwrap()).unwrap();
    assert!(doc.is_object());
    assert_eq!(
        r["result"]["contents"][0]["mimeType"],
        json!("application/schema+json")
    );
    let r = c.request(
        "resources/read",
        json!({"uri": "rusolve://docs/agent/reference"}),
    );
    assert!(r["result"]["contents"][0]["text"].as_str().unwrap().len() > 1000);
    let e = c.request("resources/read", json!({"uri": "rusolve://nothing"}));
    assert_eq!(e["error"]["code"], json!(-32002));
    // the capabilities say so
    let init = c.request(
        "initialize",
        json!({"protocolVersion": "2025-11-25", "capabilities": {}}),
    );
    assert!(init["result"]["capabilities"]["prompts"].is_object());
    assert!(init["result"]["capabilities"]["resources"].is_object());
}

/// An equality-knapsack MIP (market split) that no solver closes in seconds.
fn hard_mip_text() -> String {
    let (m, n) = (5, 45);
    let mut seed: u64 = 12345;
    let mut next = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 33) % 9000 + 1000
    };
    let mut rows = String::new();
    for r in 0..m {
        let coefs: Vec<u64> = (0..n).map(|_| next()).collect();
        let rhs: u64 = coefs.iter().sum::<u64>() / 2;
        let terms: Vec<String> = coefs
            .iter()
            .enumerate()
            .map(|(j, a)| format!("{a} x{j}"))
            .collect();
        rows += &format!(" r{r}: {} = {rhs}\n", terms.join(" + "));
    }
    let obj: Vec<String> = (0..n).map(|j| format!("x{j}")).collect();
    format!(
        "Minimize\n obj: {}\nSubject To\n{rows}Binary\n {}\nEnd\n",
        obj.join(" + "),
        obj.join(" ")
    )
}

#[test]
fn mcp_solve_reports_progress_with_a_token() {
    let mut c = Client::spawn();
    let id = c.ok("load_model", json!({"lp_text": hard_mip_text()}))["model_id"]
        .as_str()
        .unwrap()
        .to_string();
    let rid = c.send_request(
        "tools/call",
        json!({"name": "solve", "arguments": {"model_id": id, "time_limit": 1.5},
               "_meta": {"progressToken": "tok-1"}}),
    );
    let mut progress = Vec::new();
    let reply = loop {
        let m = c.read_line();
        if m.get("id") == Some(&json!(rid)) {
            break m;
        }
        assert_eq!(m["method"], json!("notifications/progress"), "{m}");
        assert!(m.get("id").is_none());
        assert_eq!(m["params"]["progressToken"], json!("tok-1"));
        progress.push(m["params"]["progress"].as_f64().unwrap());
    };
    assert!(!progress.is_empty(), "no progress notification");
    assert!(
        progress.windows(2).all(|w| w[1] > w[0]),
        "progress must increase: {progress:?}"
    );
    assert_eq!(reply["result"]["isError"], json!(false), "{reply}");
    // without a token nothing is sent: the next line is the reply itself
    let rid = c.send_request(
        "tools/call",
        json!({"name": "solve", "arguments": {"model_id": id, "time_limit": 0.5}}),
    );
    assert_eq!(c.read_line()["id"], json!(rid));
}

#[test]
fn mcp_cancelled_solve_stops_and_gets_no_response() {
    let mut c = Client::spawn();
    let id = c.ok("load_model", json!({"lp_text": hard_mip_text()}))["model_id"]
        .as_str()
        .unwrap()
        .to_string();
    let start = std::time::Instant::now();
    let rid = c.send_request(
        "tools/call",
        json!({"name": "solve", "arguments": {"model_id": id, "time_limit": 120}}),
    );
    std::thread::sleep(std::time::Duration::from_millis(700));
    c.send_notification(
        "notifications/cancelled",
        json!({"requestId": rid, "reason": "test"}),
    );
    // the server is free again: a ping is answered, and it is the first thing answered
    let pid = c.send_request("ping", json!({}));
    let m = c.read_line();
    assert_eq!(
        m["id"],
        json!(pid),
        "a cancelled request gets no response: {m}"
    );
    assert!(
        start.elapsed().as_secs_f64() < 30.0,
        "the solve ran on after the cancel ({:?})",
        start.elapsed()
    );
    // the handle survives: it solves again, and a cancel that arrives before the call starts
    // (queued behind another one) also drops it
    let r = c.ok("solve", json!({"model_id": id, "time_limit": 0.3}));
    assert!(r["status"].is_string());
    // cancel a request that is already queued behind a running solve
    let first = c.send_request(
        "tools/call",
        json!({"name": "solve", "arguments": {"model_id": id, "time_limit": 1.0}}),
    );
    let queued = c.send_request(
        "tools/call",
        json!({"name": "solve", "arguments": {"model_id": id, "time_limit": 60}}),
    );
    c.send_notification("notifications/cancelled", json!({"requestId": queued}));
    let pid = c.send_request("ping", json!({}));
    assert_eq!(c.read_line()["id"], json!(first));
    assert_eq!(c.read_line()["id"], json!(pid));
}

/// `hard_mip_text` plus a general constraint, so the router sends it to the gc-linearization
/// engine, a child process that runs the whole branch and bound.
fn hard_routed_text() -> String {
    hard_mip_text().replace(
        "Binary\n",
        "General Constraints\n g1: m = MAX ( x0 , x1 )\nBounds\n 0 <= m <= 1\nBinary\n",
    )
}

#[test]
fn mcp_cancel_kills_a_routed_engine() {
    // a routed handle's solve is a child process that would run for its whole 60 s limit; the
    // cancel must kill it (the call waits for the child after the kill, so a child left running
    // delays the next reply by the rest of the solve)
    let mut c = Client::spawn();
    let m = c.ok("load_model", json!({"lp_text": hard_routed_text()}));
    assert_eq!(m["engine"], json!("gc-linearization"), "{m}");
    let id = m["model_id"].as_str().unwrap().to_string();
    let rid = c.send_request(
        "tools/call",
        json!({"name": "solve", "arguments": {"model_id": id, "time_limit": 60}}),
    );
    std::thread::sleep(std::time::Duration::from_millis(1500));
    let cancelled_at = std::time::Instant::now();
    c.send_notification("notifications/cancelled", json!({"requestId": rid}));
    let pid = c.send_request("ping", json!({}));
    let m = c.read_line();
    assert_eq!(m["id"], json!(pid), "no reply for the cancelled solve: {m}");
    assert!(
        cancelled_at.elapsed().as_secs_f64() < 5.0,
        "the ping waited {:?}: the routed child was not killed",
        cancelled_at.elapsed()
    );
    // nothing for the solve arrives later either: the next reply is the next ping's
    let pid2 = c.send_request("ping", json!({}));
    assert_eq!(c.read_line()["id"], json!(pid2));
}

#[test]
fn mcp_load_model_errors_carry_the_reader_line() {
    // RS0011 with the file or inline text, the line the reader named, and that line's text
    let mut c = Client::spawn();
    c.handshake();
    let diag = |c: &mut Client, args: Value| -> Value {
        let text = c.err_raw("load_model", args);
        let d: Value = serde_json::from_str(&text).expect("error text is JSON");
        d["diagnostics"][0].clone()
    };
    // LP text: the objective is fine, line 4 is not a constraint
    let lp = "Minimize\n obj: x + y\nSubject To\n c1: x + y >= <> 2\nEnd\n";
    let d = diag(&mut c, json!({"lp_text": lp}));
    assert_eq!(d["code"], json!("RS0011"), "{d}");
    let at = &d["location"];
    assert_eq!(at["kind"], json!("file"), "{d}");
    assert_eq!(at["path"], json!("<lp_text>"), "{d}");
    assert_eq!(at["line"], json!(4), "{d}");
    assert_eq!(at["source_line"], json!(" c1: x + y >= <> 2"), "{d}");
    assert!(at["column"].as_u64().unwrap() >= 1);
    // a bad MPS file: line 6 holds a column entry naming a row that is not declared
    let bad = "NAME BAD\nROWS\n N  obj\n L  r1\nCOLUMNS\n    x  obj  1  r9  1\nRHS\n    rhs  r1  4\nENDATA\n";
    let dir = std::env::temp_dir().join(format!("rusolve_mcp_loc_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("bad.mps");
    std::fs::write(&file, bad).unwrap();
    let path = file.to_string_lossy().into_owned();
    let d = diag(&mut c, json!({"path": path}));
    assert_eq!(d["code"], json!("RS0011"), "{d}");
    assert_eq!(d["location"]["path"], json!(path), "{d}");
    assert_eq!(d["location"]["line"], json!(6), "{d}");
    assert_eq!(
        d["location"]["source_line"],
        json!("    x  obj  1  r9  1"),
        "{d}"
    );
    // the same text inline
    let d = diag(&mut c, json!({"mps_text": bad}));
    assert_eq!(d["location"]["path"], json!("<mps_text>"), "{d}");
    assert_eq!(d["location"]["line"], json!(6), "{d}");
    // a missing file is RS0010, with no location
    let d = diag(
        &mut c,
        json!({"path": dir.join("nope.mps").to_string_lossy()}),
    );
    assert_eq!(d["code"], json!("RS0010"), "{d}");
    assert!(d.get("location").is_none(), "{d}");
    // MathOptFormat: invalid JSON names a line and a column
    let d = diag(&mut c, json!({"mof_text": "{\n \"version\": ,\n}"}));
    assert_eq!(d["code"], json!("RS0011"), "{d}");
    assert_eq!(d["location"]["line"], json!(2), "{d}");
    assert!(d["location"]["column"].as_u64().unwrap() > 1, "{d}");
    let _ = std::fs::remove_dir_all(&dir);
}
