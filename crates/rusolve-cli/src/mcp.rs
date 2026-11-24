//! `rusolve mcp`, a Model Context
//! Protocol server over stdio, so an agent can drive the solver as tools.
//!
//! Transport: the stdio binding, one JSON-RPC 2.0 message per line on stdin and
//! stdout, no embedded newlines, logs on stderr only
//! (<https://modelcontextprotocol.io/specification/2026-07-28/basic/transports/stdio>).
//! HTTP transport is not implemented yet.
//!
//! ## Why hand-rolled, not `rmcp`
//!
//! The protocol surface used here is five methods and a JSON-RPC envelope; the
//! tools are thin wrappers over crate code that already speaks `serde_json`. The
//! current revision (2026-07-28) dropped the `initialize` handshake for
//! per-request `_meta`, and a dual-era server has to serve both; whether `rmcp`
//! tracks that revision could not be verified, while the
//! hand-rolled loop is ~300 lines and adds no dependency (`serde_json` is
//! already in the tree).
//!
//! ## Eras (spec: basic/versioning)
//!
//! <https://modelcontextprotocol.io/specification/2026-07-28/basic/versioning>:
//! the current protocol version is `2026-07-28`; it has no handshake, every
//! request carries `params._meta["io.modelcontextprotocol/protocolVersion"]`,
//! and a server that does not implement the requested version answers
//! `UnsupportedProtocolVersionError` (code -32022, `data.supported`,
//! `data.requested`). Servers MUST implement `server/discover`
//! (<https://modelcontextprotocol.io/specification/2026-07-28/server/discover>).
//! This server is dual-era (the spec's term): a request with modern `_meta` is
//! served statelessly and its result carries `resultType: "complete"`; an
//! `initialize` request selects the legacy semantics
//! (<https://modelcontextprotocol.io/specification/2025-11-25/basic/lifecycle#version-negotiation>:
//! the server echoes the requested version when it supports it, else answers
//! with another version it supports) for `2025-11-25`, `2025-06-18`,
//! `2025-03-26` and `2024-11-05`. A request with neither (a client that skips
//! `initialize`) is served leniently as legacy.
//!
//! Methods: `initialize`, `notifications/initialized` (and any other
//! notification: no response), `ping`, `server/discover`, `tools/list`,
//! `tools/call`, `prompts/list`, `prompts/get`, `resources/list`, `resources/read`
//! (`mcp_resources.rs`).
//!
//! ## Progress and cancellation
//!
//! A `tools/call` whose `params._meta.progressToken` is set gets `notifications/progress`
//! while it runs (<https://modelcontextprotocol.io/specification/2026-07-28/basic/utilities/progress>:
//! `progress` increases with every notification). `notifications/cancelled {requestId}` stops a
//! running or queued call (<https://modelcontextprotocol.io/specification/2026-07-28/basic/utilities/cancellation>):
//! the receiver "SHOULD NOT send a response" for a cancelled request, so the call ends and
//! nothing is answered. A solve blocks the thread that runs it, so a reader thread owns stdin
//! and handles the cancel notification itself (it sets the flag the solve's interrupt callback
//! polls); every other line is queued for the request loop. Both threads write stdout through
//! one lock, so lines never interleave. Anything else is `-32601`. Malformed JSON is `-32700` with a
//! null id, a message that is not a JSON-RPC request is `-32600`, bad `params`
//! (unknown tool, missing name) are `-32602`
//! (<https://modelcontextprotocol.io/specification/2026-07-28/server/tools#error-handling>).
//! Failures inside a tool (bad argument, unknown model id, solver failure) are
//! tool results with `isError: true` and an actionable message instead.
//!
//! ## Tool results
//!
//! Each tool declares `inputSchema` and `outputSchema` and returns the object in
//! `structuredContent` plus the same JSON serialised in a `text` content item
//! ("a tool that returns structured content SHOULD also return the serialized
//! JSON in a TextContent block",
//! <https://modelcontextprotocol.io/specification/2026-07-28/server/tools#structured-content>).
//! Tools are listed in a fixed order
//! (<https://modelcontextprotocol.io/specification/2026-07-28/server/tools#capabilities>).

use crate::mcp_resources;
use crate::mcp_tools::{Ctx, Registry};
use rusolve_core::diagnostic::{Diagnostic, Location};
use serde_json::{Map, Value, json};
use std::collections::HashSet;
use std::io::{BufRead, Write};
use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, mpsc};

/// The stdout writer shared by the request loop and the progress notifications of a running
/// tool: one JSON message per line, written and flushed under one lock.
pub struct Out(Mutex<std::io::Stdout>);

impl Out {
    pub fn stdout() -> Self {
        Out(Mutex::new(std::io::stdout()))
    }

    /// `false` when stdout is closed.
    pub fn send(&self, msg: &Value) -> bool {
        let mut w = self.0.lock().unwrap_or_else(|e| e.into_inner());
        // serde_json's compact form has no newline, as the stdio binding requires.
        writeln!(w, "{msg}").and_then(|_| w.flush()).is_ok()
    }
}

/// Which requests were cancelled, written by the reader thread and read by the request loop.
/// A cancel can arrive before the request starts (it is queued behind a running solve) or while
/// it runs; both orders end in the running call's flag being set.
pub struct Control {
    cancelled: Mutex<HashSet<String>>,
    active: Mutex<Option<String>>,
    flag: Arc<AtomicBool>,
}

impl Control {
    pub fn new() -> Arc<Self> {
        Arc::new(Self {
            cancelled: Mutex::new(HashSet::new()),
            active: Mutex::new(None),
            flag: Arc::new(AtomicBool::new(false)),
        })
    }

    fn cancel(&self, key: String) {
        self.cancelled
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .insert(key.clone());
        if self
            .active
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .as_deref()
            == Some(key.as_str())
        {
            self.flag.store(true, Ordering::SeqCst);
        }
    }

    /// Starts request `key`: returns its flag, already set when a cancel came first.
    fn begin(&self, key: &str) -> Arc<AtomicBool> {
        self.flag.store(false, Ordering::SeqCst);
        *self.active.lock().unwrap_or_else(|e| e.into_inner()) = Some(key.to_string());
        if self
            .cancelled
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .contains(key)
        {
            self.flag.store(true, Ordering::SeqCst);
        }
        self.flag.clone()
    }

    fn end(&self, key: &str) {
        *self.active.lock().unwrap_or_else(|e| e.into_inner()) = None;
        self.cancelled
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .remove(key);
    }
}

/// The key a request id is tracked under (an id is a string or a number).
fn id_key(id: &Value) -> String {
    id.to_string()
}

/// `Some(key)` of a `notifications/cancelled` line.
fn cancel_target(line: &str) -> Option<String> {
    if !line.contains("notifications/cancelled") {
        return None;
    }
    let v: Value = serde_json::from_str(line).ok()?;
    if v.get("method")?.as_str()? != "notifications/cancelled" || v.get("id").is_some() {
        return None;
    }
    Some(id_key(v.get("params")?.get("requestId")?))
}

/// The revision whose `_meta` protocol version the server serves statelessly.
pub const MODERN_VERSION: &str = "2026-07-28";
/// The handshake-based revisions `initialize` negotiates, newest first.
pub const LEGACY_VERSIONS: [&str; 4] = ["2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05"];
const META_VERSION_KEY: &str = "io.modelcontextprotocol/protocolVersion";
const META_SERVER_INFO_KEY: &str = "io.modelcontextprotocol/serverInfo";

const PARSE_ERROR: i64 = -32700;
const INVALID_REQUEST: i64 = -32600;
const METHOD_NOT_FOUND: i64 = -32601;
const INVALID_PARAMS: i64 = -32602;
/// `UnsupportedProtocolVersionError` (versioning page).
const UNSUPPORTED_VERSION: i64 = -32022;

const INSTRUCTIONS: &str = "rusolve solves LPs and MIPs (a faithful HiGHS port), convex QPs, SOCPs and convex QCQPs, \
and general-constraint models; a convex MIQP and a mixed-integer conic model load too, and solve runs their experimental \
engine only with the option miqp_oa or misocp_oa set to true (off by default, slower). Call load_model first (path, or mps_text / lp_text / mof_text), or build_model \
with JSON variables and rows; it returns a model_id that every other tool takes. solve returns a rusolve.result/1 \
envelope. After modify (edits by row/column name) an LP re-solves from the kept basis, so batteries of \
what-if solves are cheap. explain gives an infeasibility or unboundedness certificate by name, feasrelax \
the cheapest bound relaxation, lint static modelling problems, ranging sensitivity. free_model releases a handle. \
The prompts are workflows and problem-class playbooks (solve, debug-infeasible, tune, extend, conic, convex-qp, ...); the resources are the JSON schemas and the \
agent reference.";

fn capabilities() -> Value {
    json!({
        "tools": {"listChanged": false},
        "prompts": {"listChanged": false},
        "resources": {"subscribe": false, "listChanged": false},
    })
}

fn supported_versions() -> Vec<&'static str> {
    std::iter::once(MODERN_VERSION)
        .chain(LEGACY_VERSIONS)
        .collect()
}

fn error_response(id: Value, code: i64, message: &str, data: Option<Value>) -> Value {
    let mut err = json!({"code": code, "message": message});
    if let Some(d) = data {
        err["data"] = d;
    }
    json!({"jsonrpc": "2.0", "id": id, "error": err})
}

fn ok_response(id: Value, result: Value) -> Value {
    json!({"jsonrpc": "2.0", "id": id, "result": result})
}

/// The protocol state: the model handles, and nothing else. There is no
/// session-global "current model".
pub struct Server {
    registry: Registry,
    /// The version `initialize` negotiated (legacy era); informational.
    negotiated: Option<String>,
    out: Arc<Out>,
    control: Arc<Control>,
}

/// An error code of the request loop that means "answer nothing" (a cancelled request).
const NO_RESPONSE: i64 = 0;

impl Default for Server {
    fn default() -> Self {
        Self::new()
    }
}

impl Server {
    pub fn new() -> Self {
        Self::with_control(Arc::new(Out::stdout()), Control::new())
    }

    pub fn with_control(out: Arc<Out>, control: Arc<Control>) -> Self {
        Self {
            registry: Registry::new(),
            negotiated: None,
            out,
            control,
        }
    }

    /// Handles one decoded JSON-RPC message (a request, a notification, or a
    /// batch). `None` means nothing is sent back.
    pub fn handle_message(&mut self, msg: Value) -> Option<Value> {
        if let Value::Array(items) = msg {
            // JSON-RPC batch (2025-03-26 allowed it; later revisions dropped it).
            if items.is_empty() {
                return Some(error_response(
                    Value::Null,
                    INVALID_REQUEST,
                    "empty batch",
                    None,
                ));
            }
            let out: Vec<Value> = items
                .into_iter()
                .filter_map(|m| self.handle_one(m))
                .collect();
            return (!out.is_empty()).then_some(Value::Array(out));
        }
        self.handle_one(msg)
    }

    fn handle_one(&mut self, msg: Value) -> Option<Value> {
        let Some(obj) = msg.as_object() else {
            return Some(error_response(
                Value::Null,
                INVALID_REQUEST,
                "a JSON-RPC message must be an object",
                None,
            ));
        };
        let id = obj.get("id").cloned();
        let Some(method) = obj.get("method").and_then(Value::as_str) else {
            // A response from the client (we send no requests) or garbage.
            return if obj.contains_key("result") || obj.contains_key("error") {
                None
            } else {
                Some(error_response(
                    id.unwrap_or(Value::Null),
                    INVALID_REQUEST,
                    "missing string field \"method\"",
                    None,
                ))
            };
        };
        let params = obj.get("params").cloned().unwrap_or(Value::Null);
        let Some(id) = id else {
            // Notification (notifications/initialized, notifications/cancelled,
            // ...): never answered.
            return None;
        };
        if !params.is_object() && !params.is_null() {
            return Some(error_response(
                id,
                INVALID_PARAMS,
                "params must be an object",
                None,
            ));
        }
        // Per-request protocol version (modern era).
        let requested = params
            .get("_meta")
            .and_then(|m| m.get(META_VERSION_KEY))
            .and_then(Value::as_str)
            .map(str::to_string);
        let mut modern = false;
        if let Some(v) = &requested {
            if !supported_versions().contains(&v.as_str()) {
                return Some(error_response(
                    id,
                    UNSUPPORTED_VERSION,
                    "Unsupported protocol version",
                    Some(json!({"supported": supported_versions(), "requested": v})),
                ));
            }
            modern = v == MODERN_VERSION;
        }
        let result = match method {
            "initialize" => self.initialize(&params),
            "ping" => Ok(json!({})),
            "server/discover" => Ok(discover_result()),
            "tools/list" => Ok(json!({"tools": self.registry.tool_list()})),
            "tools/call" => self.tools_call(&id, &params),
            "prompts/list" => Ok(mcp_resources::prompts_list()),
            "prompts/get" => match params.get("name").and_then(Value::as_str) {
                Some(n) => mcp_resources::prompts_get(n).map_err(|e| (INVALID_PARAMS, e)),
                None => Err((
                    INVALID_PARAMS,
                    "prompts/get: params.name is required (a name from prompts/list)".into(),
                )),
            },
            "resources/list" => Ok(mcp_resources::resources_list()),
            "resources/read" => match params.get("uri").and_then(Value::as_str) {
                Some(u) => mcp_resources::resources_read(u)
                    .map_err(|e| (mcp_resources::RESOURCE_NOT_FOUND, e)),
                None => Err((
                    INVALID_PARAMS,
                    "resources/read: params.uri is required (a uri from resources/list)".into(),
                )),
            },
            other => Err((METHOD_NOT_FOUND, format!("method not found: {other}"))),
        };
        Some(match result {
            Ok(mut r) => {
                if modern && let Some(o) = r.as_object_mut() {
                    o.insert("resultType".into(), json!("complete"));
                }
                ok_response(id, r)
            }
            Err((NO_RESPONSE, _)) => return None,
            Err((code, message)) => error_response(id, code, &message, None),
        })
    }

    fn initialize(&mut self, params: &Value) -> Result<Value, (i64, String)> {
        let Some(requested) = params.get("protocolVersion").and_then(Value::as_str) else {
            return Err((
                INVALID_PARAMS,
                "initialize: params.protocolVersion is required (a string such as \"2025-11-25\")"
                    .into(),
            ));
        };
        // Lifecycle, version negotiation: echo the requested version if
        // supported, otherwise answer with the latest one we support.
        let chosen = if LEGACY_VERSIONS.contains(&requested) {
            requested
        } else {
            LEGACY_VERSIONS[0]
        };
        self.negotiated = Some(chosen.to_string());
        Ok(json!({
            "protocolVersion": chosen,
            "capabilities": capabilities(),
            "serverInfo": server_info(),
            "instructions": INSTRUCTIONS,
        }))
    }

    fn tools_call(&mut self, id: &Value, params: &Value) -> Result<Value, (i64, String)> {
        let Some(name) = params.get("name").and_then(Value::as_str) else {
            return Err((
                INVALID_PARAMS,
                "tools/call: params.name is required (a tool name from tools/list)".into(),
            ));
        };
        if !self.registry.has_tool(name) {
            return Err((INVALID_PARAMS, format!("Unknown tool: {name}")));
        }
        let args = match params.get("arguments") {
            None | Some(Value::Null) => Map::new(),
            Some(Value::Object(m)) => m.clone(),
            Some(_) => {
                return Err((
                    INVALID_PARAMS,
                    "tools/call: params.arguments must be an object".into(),
                ));
            }
        };
        // A panic inside a tool (a port assertion, a vanished model file) must
        // not take the server and every other handle down with it.
        let key = id_key(id);
        let token = params
            .get("_meta")
            .and_then(|m| m.get("progressToken"))
            .cloned();
        let cancel = self.control.begin(&key);
        let ctx = Arc::new(Ctx::new(token, cancel.clone(), self.out.clone()));
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if cancel.load(Ordering::SeqCst) {
                // cancelled while queued: never started; the Err is dropped below, the
                // request gets no reply
                return Err("cancelled".to_string());
            }
            self.registry.call(name, &args, &ctx)
        }))
        .unwrap_or_else(|p| {
            let what = p
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "unknown panic".into());
            Err(format!(
                "internal error in tool {name}: {what}; the handle may be unusable, free_model it and load again"
            ))
        });
        let cancelled = cancel.load(Ordering::SeqCst);
        self.control.end(&key);
        if cancelled {
            // cancellation spec: the receiver SHOULD NOT answer a cancelled request
            return Err((NO_RESPONSE, String::new()));
        }
        Ok(tool_result(outcome, &args))
    }
}

fn server_info() -> Value {
    json!({
        "name": "rusolve",
        "title": "rusolve LP/MIP solver",
        "version": env!("CARGO_PKG_VERSION"),
    })
}

fn discover_result() -> Value {
    json!({
        "supportedVersions": supported_versions(),
        "capabilities": capabilities(),
        "_meta": {META_SERVER_INFO_KEY: server_info()},
        "instructions": INSTRUCTIONS,
    })
}

/// A `CallToolResult`: structured content plus the same JSON as text, or
/// `isError: true` with the message as text
/// (<https://modelcontextprotocol.io/specification/2026-07-28/server/tools#tool-result>).
fn tool_result(outcome: Result<Value, String>, args: &Map<String, Value>) -> Value {
    match outcome {
        Ok(v) => json!({
            "content": [{"type": "text", "text": v.to_string()}],
            "structuredContent": v,
            "isError": false,
        }),
        Err(msg) => json!({
            "content": [{"type": "text", "text": error_diagnostic(&msg, args).mcp_error_content(1)}],
            "isError": true,
        }),
    }
}

/// Wraps a tool failure message in a registered diagnostic (docs/design/diagnostics.md) so an
/// agent parses the same `rusolve.diagnostics/1` document as on the CLI. The code follows the
/// failure: a model file that is missing RS0010, a model that does not parse RS0011 (with the
/// reader's line as a location), a class nothing solves RS0012, an unknown handle RS3002, a
/// bad argument RS3001, anything else RS0090. (A cancelled call is never answered: MCP's
/// cancellation spec says the receiver SHOULD NOT reply, so no code exists for it.)
fn error_diagnostic(msg: &str, args: &Map<String, Value>) -> Diagnostic {
    let missing_file = || {
        msg.strip_prefix("path: read ")
            .and_then(|r| r.split_once(": "))
            .is_some_and(|(p, _)| std::fs::File::open(p).is_err())
    };
    if msg.contains("models are not supported") {
        Diagnostic::error("RS0012", msg)
            .with_help("rusolve classify names the class and the engine that would take it")
    } else if msg.contains("call load_model first") {
        Diagnostic::error("RS3002", msg)
            .with_help("call load_model or build_model and use the model_id it returns")
    } else if missing_file() {
        Diagnostic::error("RS0010", msg)
            .with_help("check the path; `rusolve --explain RS0010` describes the failure")
    } else if msg.contains(" parse error")
        || msg.contains(" read error")
        || msg.starts_with("path: read ")
        || msg.starts_with("build_model:")
    {
        let mut d = Diagnostic::error("RS0011", msg).with_help(
            "fix the model text and load it again; `rusolve --explain RS0011` describes the failure",
        );
        d.location = reader_location(msg, args);
        d
    } else if msg.starts_with("internal error") {
        Diagnostic::error("RS0090", msg).with_help("report the model and the call")
    } else {
        Diagnostic::error("RS3001", msg).with_help("tools/list gives the inputSchema of each tool")
    }
}

/// Where a model reader failed. LP and MPS readers say `at line N` (rusolve-io `lp.rs`,
/// `mps.rs`), MathOptFormat says `at line N column M` for bad JSON; the source line is read back
/// from the inline text (`<lp_text>`, `<mps_text>`, `<mof_text>`) or the file named by `path`.
/// A reader that names no line (a bad MOF field, a missing file) gives no location. A column the
/// reader does not give is the first non-blank one, a one-character caret.
fn reader_location(msg: &str, args: &Map<String, Value>) -> Option<Location> {
    let rest = &msg[msg.find(" at line ")? + " at line ".len()..];
    let digits = |s: &str| -> Option<usize> {
        let n: String = s.chars().take_while(char::is_ascii_digit).collect();
        n.parse().ok()
    };
    let line = digits(rest)?;
    let reader_col = rest
        .split_once(" column ")
        .and_then(|(l, c)| l.chars().all(|ch| ch.is_ascii_digit()).then(|| digits(c))?);
    let what = msg.split_once(": ")?.0;
    let (path, text) = match what {
        "lp_text" | "mps_text" | "mof_text" => {
            (format!("<{what}>"), args.get(what)?.as_str()?.to_string())
        }
        "path" => {
            let p = args.get("path")?.as_str()?;
            (p.to_string(), std::fs::read_to_string(p).ok()?)
        }
        _ => return None,
    };
    let source_line = text.lines().nth(line.checked_sub(1)?)?.to_string();
    let column = reader_col.unwrap_or_else(|| {
        source_line
            .chars()
            .position(|c| !c.is_whitespace())
            .map_or(1, |i| i + 1)
    });
    Some(Location::File {
        path,
        line,
        column,
        len: 1,
        source_line,
    })
}

/// `rusolve mcp`: the stdio loop. Exits 0 at end of input.
pub fn run(args: &[String]) -> ExitCode {
    if let Some(a) = args.first() {
        eprintln!("usage: rusolve mcp   (stdio transport only; unexpected argument '{a}')");
        return ExitCode::from(2);
    }
    let out = Arc::new(Out::stdout());
    let control = Control::new();
    let mut server = Server::with_control(out.clone(), control.clone());
    // The reader thread owns stdin so a cancel is seen while a solve runs.
    let (tx, rx) = mpsc::channel::<String>();
    std::thread::spawn(move || {
        let stdin = std::io::stdin();
        for line in stdin.lock().lines() {
            let Ok(line) = line else { break };
            if let Some(key) = cancel_target(&line) {
                control.cancel(key);
                continue;
            }
            if tx.send(line).is_err() {
                break;
            }
        }
    });
    for line in rx {
        if line.trim().is_empty() {
            continue;
        }
        let reply = match serde_json::from_str::<Value>(&line) {
            Ok(msg) => server.handle_message(msg),
            Err(e) => Some(error_response(
                Value::Null,
                PARSE_ERROR,
                &format!("parse error: {e}"),
                None,
            )),
        };
        if let Some(r) = reply
            && !out.send(&r)
        {
            break;
        }
    }
    server.registry.shutdown();
    ExitCode::SUCCESS
}
