//! CLI side of the diagnostics (`docs/design/diagnostics.md`): the one funnel every error
//! path goes through, the `--explain` entry, and the did-you-mean helpers over the registries.
//!
//! The type, registry and renderers live in `rusolve_core::diagnostic`.

use rusolve_core::diagnostic::{self, Diagnostic};
use std::process::ExitCode;

/// `true` when the command line asks for `--format json` (or `--format=json`). Read before the
/// command's own parser runs, so even a parse error can be answered in JSON.
pub fn wants_json(args: &[String]) -> bool {
    args.windows(2)
        .any(|w| w[0] == "--format" && w[1] == "json")
        || args.iter().any(|a| a == "--format=json")
}

static DEFER_CONIC_JSON: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// `main` arms this around a `solve`/`mip` run whose read error may hand over to the conic route
/// (the safety net in `main.rs`): a JSON read error of the first attempt is then not printed,
/// so the command still ends with exactly one JSON document on stdout.
pub fn defer_conic_json(on: bool) {
    DEFER_CONIC_JSON.store(on, std::sync::atomic::Ordering::SeqCst);
}

/// `true` when a JSON error document must not be printed now: the conic safety net will rerun
/// the command and answer in its place.
pub fn json_deferred() -> bool {
    DEFER_CONIC_JSON.load(std::sync::atomic::Ordering::SeqCst)
        && rusolve_io::conic::conic_rejected()
}

/// Report `d` and return the exit code: the rustc-style text on stderr always, and under
/// `--format json` also the `rusolve.diagnostics/1` document on stdout. `exit` is the code
/// of the exit-code table (2 usage, 3 input, 4 internal).
pub fn fail(json: bool, d: Diagnostic, exit: u8) -> ExitCode {
    eprint!("{}", d.render_text());
    if json && !json_deferred() {
        println!("{}", diagnostic::diagnostics_json(&[d], i32::from(exit)));
    }
    ExitCode::from(exit)
}

/// An error with a code and a message, no location.
pub fn fail_msg(json: bool, code: &'static str, exit: u8, message: impl Into<String>) -> ExitCode {
    fail(json, Diagnostic::error(code, message), exit)
}

/// A model or input file that could not be read (exit 3): RS0010 when the file is missing or
/// unreadable, RS0011 when it was read but does not parse. `detail` is the reader's text.
pub fn fail_read(json: bool, path: &str, detail: &str) -> ExitCode {
    let missing = std::fs::File::open(path).is_err();
    // A model that was read but does not parse is reported in full: every independent
    // error of the file (at most 20) with line, column, source line and help
    if !missing {
        let all = rusolve_io::diag::read_diagnostics(path);
        if !all.is_empty() {
            for d in &all {
                eprint!("{}", d.render_text());
            }
            if json && !json_deferred() {
                println!("{}", diagnostic::diagnostics_json(&all, 3));
            }
            return ExitCode::from(3);
        }
    }
    let d = if missing {
        Diagnostic::error("RS0010", format!("read {path}: {detail}"))
            .with_help("check the path; `rusolve --explain RS0010`")
    } else {
        Diagnostic::error("RS0011", format!("read {path}: {detail}"))
            .with_help("`rusolve --explain RS0011`")
    };
    fail(json, d, 3)
}

/// [`fail_read`] for the `Err(String)` of `mip_highs::read_mip*`, which already carries the
/// `read <path>: ` prefix. A model with general constraints that the command cannot take is
/// RS0012 (not supported), not a parse error.
pub fn fail_read_str(json: bool, path: &str, err: &str) -> ExitCode {
    let detail = err.strip_prefix(&format!("read {path}: ")).unwrap_or(err);
    if detail.contains("general constraints are supported by") {
        return fail_msg(json, "RS0012", 3, format!("read {path}: {detail}"));
    }
    fail_read(json, path, detail)
}

/// The registered usage code (RS0002 unknown flag, RS0003 bad or missing flag value, RS0006
/// missing or unexpected argument) that fits a usage message, by its wording.
pub fn usage_code(msg: &str) -> &'static str {
    let m = msg.to_ascii_lowercase();
    if m.contains("unknown option")
        || m.contains("unknown flag")
        || m.contains("unsupported flag")
        || m.contains("is not supported")
        || m.contains("unknown subcommand")
    {
        "RS0002"
    } else if m.contains("missing")
        || m.contains("more than one")
        || m.starts_with("usage")
        || m.contains("unexpected")
    {
        "RS0006"
    } else {
        "RS0003"
    }
}

/// A usage error (exit 2) of `command` from a free-text message: its first line is the message,
/// the following lines (the usage text) the help. The code comes from [`usage_code`].
pub fn fail_usage_text(json: bool, command: &str, msg: &str) -> ExitCode {
    let mut lines = msg.lines();
    let first = lines.next().unwrap_or("");
    let mut d = Diagnostic::error(usage_code(msg), format!("rusolve {command}: {first}"));
    let rest: Vec<&str> = lines.collect();
    if !rest.is_empty() {
        d = d.with_help(rest.join("\n"));
    }
    fail(json, d, 2)
}

/// A [`crate::cli_core::CoreError`] of a shared core (explain, ranging, feasrelax): usage text
/// for exit 2, RS0012 for an input the command cannot honour (3), RS0090 for an internal failure.
pub fn fail_core(json: bool, command: &str, e: &crate::cli_core::CoreError) -> ExitCode {
    match e.code {
        2 => fail_usage_text(json, command, &e.msg),
        3 => fail_msg(json, "RS0012", 3, format!("rusolve {command}: {}", e.msg)),
        c => fail_msg(json, "RS0090", c, format!("rusolve {command}: {}", e.msg)),
    }
}

/// Marks an exit-2 message of [`fail_coded`] as a content error (RS0011), not a usage error.
/// Stripped before printing; the lines after the first become notes.
pub const INVALID_CONTENT: &str = "\u{1}invalid: ";

/// The `(exit, message)` error pair of the profile and tune code: 2 is a usage error, 3 an input
/// that is missing (RS0010) or does not parse (RS0011, by its wording), anything else internal.
pub fn fail_coded(json: bool, command: &str, exit: u8, msg: &str) -> ExitCode {
    match exit {
        2 if msg.starts_with(INVALID_CONTENT) => {
            // a file that parses but breaks its schema: RS0011, one note per further error
            let mut lines = msg[INVALID_CONTENT.len()..].lines();
            let mut d = Diagnostic::error(
                "RS0011",
                format!("rusolve {command}: {}", lines.next().unwrap_or("")),
            );
            for l in lines {
                d = d.with_note(l);
            }
            fail(json, d, 2)
        }
        2 => fail_usage_text(json, command, msg),
        3 => {
            let code = if msg.contains("not JSON") || msg.contains("parse") {
                "RS0011"
            } else {
                "RS0010"
            };
            fail_msg(json, code, 3, format!("rusolve {command}: {msg}"))
        }
        c => fail_msg(json, "RS0090", c, format!("rusolve {command}: {msg}")),
    }
}

/// A [`crate::commands::UsageError`] as exit 2.
pub fn fail_usage(json: bool, e: &crate::commands::UsageError) -> ExitCode {
    fail(json, e.diagnostic(), 2)
}

/// An unknown command (RS0001) with the closest command names as help.
pub fn unknown_command(json: bool, cmd: &str) -> ExitCode {
    let mut d = Diagnostic::error(
        "RS0001",
        format!("rusolve: unknown command '{cmd}' (see `rusolve --help`)"),
    );
    d.help.extend(diagnostic::did_you_mean(
        cmd,
        crate::commands::COMMANDS.iter().map(|c| c.name),
    ));
    fail(json, d, 2)
}

/// An unknown option name (`--opt`, options file): RS0004 with the closest registry names.
pub fn unknown_option(json: bool, command: &str, name: &str) -> ExitCode {
    let mut d = Diagnostic::error(
        "RS0004",
        format!("rusolve {command}: unknown option '{name}'"),
    );
    // the registry's own hint: `did you mean mip_rel_gap?` or `closest options are ...`
    d.help
        .extend(rusolve_ekk::lp_data::highs_option_records::unknown_option_hint(name));
    fail(json, d, 2)
}

/// `rusolve --explain RSnnnn | all [--format json]`: `Some(code)` when args are an explain call.
pub fn explain(args: &[String]) -> Option<ExitCode> {
    let first = args.first()?;
    let (flag, rest) = match first.split_once('=') {
        Some(("--explain", v)) => ("--explain", Some(v.to_string())),
        _ if first == "--explain" => ("--explain", args.get(1).cloned()),
        _ => return None,
    };
    debug_assert_eq!(flag, "--explain");
    let json = wants_json(args);
    let Some(what) = rest else {
        return Some(fail_msg(
            json,
            "RS0003",
            2,
            "rusolve --explain needs a code, for example `rusolve --explain RS0002`, or `all`",
        ));
    };
    if what == "all" {
        println!("{}", diagnostic::registry_json());
        return Some(ExitCode::SUCCESS);
    }
    match diagnostic::lookup(&what) {
        Some(e) => {
            print!("{}", diagnostic::explain_text(e));
            Some(ExitCode::SUCCESS)
        }
        None => {
            let mut d = Diagnostic::error("RS0003", format!("no diagnostic code '{what}'"));
            d.help.extend(diagnostic::did_you_mean(
                &what,
                diagnostic::REGISTRY.iter().map(|e| e.code),
            ));
            d.help
                .push("`rusolve --explain all` lists every code".to_string());
            Some(fail(json, d, 2))
        }
    }
}
