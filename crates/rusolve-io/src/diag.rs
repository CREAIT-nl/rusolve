//! Parse diagnostics (`docs/design/diagnostics.md`): every parse error
//! of an MPS, LP or MathOptFormat file as a [`Diagnostic`] with line, column, source line and help,
//! and all of them in one pass (at most [`MAX_ERRORS`]) instead of the first.
//!
//! The readers stay strict and fast: they return their first error. This module runs on the error
//! path only. For MPS and LP it re-reads the file with the line of the reported error turned into a
//! comment, so the next error surfaces, and repeats; recovery continues only where skipping one
//! data line cannot change how the lines after it are read (a bad number, an undeclared row or
//! column on a COLUMNS/RHS/RANGES/BOUNDS line). An error in ROWS, a section header or a
//! structural error ends the list, since what follows would be a cascade, not a finding.
//!
//! Layout and recovery follow rustc's guidance that a front end should report as many independent
//! errors as it can and never invent dependent ones
//! (<https://rustc-dev-guide.rust-lang.org/diagnostics.html>); the did-you-mean rule is
//! [`rusolve_core::diagnostic::did_you_mean`].

use rusolve_core::diagnostic::{Diagnostic, Location, did_you_mean};

use crate::{Format, format_from_extension, read_file_bytes, sniff_format};

/// Most errors reported for one file.
pub const MAX_ERRORS: usize = 20;
/// Above this many bytes of text only the first error is reported (each recovery step re-reads).
const RECOVERY_MAX_BYTES: usize = 32 << 20;
/// Declared names listed in a help line.
const LISTED: usize = 10;

/// The diagnostics of the model file at `path`: empty when it parses (or the format has no
/// recovery to offer, `.nl`), otherwise one error per independent problem, in file order.
pub fn read_diagnostics(path: &str) -> Vec<Diagnostic> {
    let Ok(bytes) = read_file_bytes(path) else {
        return Vec::new();
    };
    let fmt = format_from_extension(path).unwrap_or_else(|| sniff_format(&bytes));
    match fmt {
        Format::Mof => mof_diagnostics(path, &bytes),
        Format::Mps | Format::Lp => match String::from_utf8(bytes) {
            Ok(text) => text_diagnostics(path, &text, fmt),
            Err(_) => Vec::new(),
        },
        Format::Nl => Vec::new(),
    }
}

/// [`read_diagnostics`] on text already in memory (`name` is shown as the path).
pub fn text_diagnostics(name: &str, text: &str, fmt: Format) -> Vec<Diagnostic> {
    let mut lines: Vec<String> = text.lines().map(str::to_string).collect();
    let mut out: Vec<Diagnostic> = Vec::new();
    let declared = Declared::scan(text);
    let recover = text.len() <= RECOVERY_MAX_BYTES;
    while out.len() < MAX_ERRORS {
        let joined = lines.join("\n");
        let (line, message) = match fmt {
            Format::Mps => match crate::mps::read_mps_conic(&joined) {
                Ok(_) => break,
                Err(e) => (e.line, e.message),
            },
            Format::Lp => match crate::lp::read_lp_conic(&joined) {
                Ok(_) => break,
                Err(e) => (e.line, e.message),
            },
            _ => break,
        };
        // an "undeclared" name the file does declare only lost its declaration to an earlier
        // skipped line: a dependent error, never reported (rule: do not invent errors)
        if fmt == Format::Mps
            && line >= 1
            && line <= lines.len()
            && is_cascade(&message, &lines[line - 1], &declared)
        {
            lines[line - 1] = "*".into();
            continue;
        }
        let d = if line >= 1 && line <= lines.len() {
            located(name, fmt, line, &lines[line - 1], &message, &declared)
        } else {
            Diagnostic::error("RS1006", format!("{}: {message}", fmt_name(fmt)))
        };
        out.push(d);
        let safe = recover && line >= 1 && line <= lines.len() && skippable(fmt, &message);
        if !safe {
            break;
        }
        // a comment line in either format (`*` in column 1 for MPS, `\` for LP)
        lines[line - 1] = if fmt == Format::Mps {
            "*".into()
        } else {
            "\\".into()
        };
    }
    out
}

fn is_cascade(message: &str, source: &str, declared: &Declared) -> bool {
    if !message.contains("undeclared") {
        return false;
    }
    let is_row = message.contains("undeclared row");
    let mut tok = quoted(message);
    if message.starts_with("bounds reference undeclared column") {
        let t: Vec<&str> = source.split_whitespace().collect();
        if t.len() >= 4 {
            tok = Some(t[2]);
        }
    }
    let names = if is_row {
        &declared.rows
    } else {
        &declared.cols
    };
    tok.is_some_and(|n| names.iter().any(|x| x == n))
}

fn fmt_name(f: Format) -> &'static str {
    match f {
        Format::Mps => "MPS",
        Format::Lp => "LP",
        Format::Nl => "NL",
        Format::Mof => "MathOptFormat",
    }
}

/// Whether the line of this error can be dropped without changing how the following lines read.
fn skippable(fmt: Format, message: &str) -> bool {
    match fmt {
        Format::Mps => {
            message.contains("undeclared")
                || message.starts_with("BOUNDS entry needs")
                || (message.starts_with("bad ")
                    && !message.starts_with("bad row sense")
                    && !message.starts_with("bad OBJSENSE")
                    && !message.starts_with("bad marker"))
        }
        Format::Lp => message.starts_with("bad number"),
        _ => false,
    }
}

/// Names declared in the file, for help lines.
struct Declared {
    rows: Vec<String>,
    cols: Vec<String>,
}

impl Declared {
    fn scan(text: &str) -> Declared {
        let mut rows = Vec::new();
        let mut cols: Vec<String> = Vec::new();
        let mut section = "";
        for raw in text.lines() {
            if raw.starts_with('*') {
                continue;
            }
            let t: Vec<&str> = raw.split_whitespace().collect();
            let Some(&first) = t.first() else { continue };
            if !raw.starts_with(char::is_whitespace)
                && matches!(
                    first,
                    "NAME" | "ROWS" | "COLUMNS" | "RHS" | "RANGES" | "BOUNDS" | "ENDATA"
                )
                && (t.len() == 1 || first == "NAME")
            {
                section = match first {
                    "ROWS" => "ROWS",
                    "COLUMNS" => "COLUMNS",
                    _ => "",
                };
                continue;
            }
            match section {
                "ROWS" if t.len() >= 2 => rows.push(t[1].to_string()),
                "COLUMNS"
                    if t.len() >= 3
                        && t[1] != "'MARKER'"
                        && cols.last().map(String::as_str) != Some(first) =>
                {
                    cols.push(first.to_string());
                }
                _ => {}
            }
        }
        Declared { rows, cols }
    }
}

fn list_help(kind: &str, names: &[String]) -> String {
    let shown: Vec<&str> = names.iter().take(LISTED).map(String::as_str).collect();
    let more = if names.len() > LISTED {
        format!(", ... ({} declared)", names.len())
    } else {
        format!(" ({} declared)", names.len())
    };
    format!("declared {kind}: {}{more}", shown.join(", "))
}

/// The quoted token of a reader message (`'R7'`), if any.
fn quoted(message: &str) -> Option<&str> {
    let b = message.rfind('\'')?;
    let a = message[..b].rfind('\'')?;
    Some(&message[a + 1..b])
}

/// Column (1-based, in chars) and width of `token` in `source`, else of its first token.
fn span_of(source: &str, token: Option<&str>) -> (usize, usize) {
    let find = |tok: &str| {
        source.match_indices(tok).find(|(i, _)| {
            let before = source[..*i].chars().next_back();
            let after = source[i + tok.len()..].chars().next();
            before.is_none_or(char::is_whitespace) && after.is_none_or(char::is_whitespace)
        })
    };
    if let Some(tok) = token.filter(|t| !t.is_empty())
        && let Some((i, _)) = find(tok)
    {
        return (source[..i].chars().count() + 1, tok.chars().count().max(1));
    }
    let lead = source.len() - source.trim_start().len();
    let width = source
        .split_whitespace()
        .next()
        .map_or(1, |w| w.chars().count());
    (source[..lead].chars().count() + 1, width.max(1))
}

fn located(
    path: &str,
    fmt: Format,
    line: usize,
    source: &str,
    message: &str,
    declared: &Declared,
) -> Diagnostic {
    let mut tok = quoted(message);
    // `UP set col value`: the reader names the set when the column is the undeclared one
    if message.starts_with("bounds reference undeclared column") {
        let t: Vec<&str> = source.split_whitespace().collect();
        if t.len() >= 4 {
            tok = Some(t[2]);
        }
    }
    let (column, len) = span_of(source, tok);
    let loc = Location::File {
        path: path.to_string(),
        line,
        column,
        len,
        source_line: source.to_string(),
    };
    let fmt_n = fmt_name(fmt);
    if fmt == Format::Mps && message.contains("undeclared") {
        let is_row = message.contains("undeclared row");
        let names = if is_row {
            &declared.rows
        } else {
            &declared.cols
        };
        let kind = if is_row { "row" } else { "column" };
        let name = tok.unwrap_or("");
        let message = &message.replace(quoted(message).unwrap_or(""), name);
        let mut d = Diagnostic::error("RS1004", format!("{kind} '{name}' is not declared"))
            .with_location(loc)
            .with_note(format!("{fmt_n} reader: {message}"));
        if let Some(h) = did_you_mean(name, names.iter().map(String::as_str)) {
            d = d.with_help(h);
        }
        return d.with_help(list_help(if is_row { "rows" } else { "columns" }, names));
    }
    if message.starts_with("bad ") && message.contains("value") || message.starts_with("bad number")
    {
        return Diagnostic::error(
            "RS1005",
            format!("'{}' is not a number", tok.unwrap_or("this value")),
        )
        .with_location(loc)
        .with_note(format!("{fmt_n} reader: {message}"))
        .with_help("write the value as a decimal or exponent number such as 12, -0.5 or 1e6; `inf` and `-inf` are accepted for bounds");
    }
    Diagnostic::error("RS1006", format!("{fmt_n} syntax: {message}"))
        .with_location(loc)
        .with_help("`rusolve --explain RS1006`")
}

fn mof_diagnostics(path: &str, bytes: &[u8]) -> Vec<Diagnostic> {
    let Err(e) = crate::mof::read_mof(bytes) else {
        return Vec::new();
    };
    let mut d = Diagnostic::error("RS1006", format!("MathOptFormat: {}", e.message));
    // invalid JSON carries `line N column M`; point at it
    let mut it = e.location.split_whitespace();
    if let (Some("line"), Some(l), Some("column"), Some(c)) =
        (it.next(), it.next(), it.next(), it.next())
        && let (Ok(l), Ok(c)) = (l.parse::<usize>(), c.parse::<usize>())
        && l >= 1
    {
        let src = String::from_utf8_lossy(bytes)
            .lines()
            .nth(l - 1)
            .unwrap_or("")
            .to_string();
        d = d.with_location(Location::File {
            path: path.to_string(),
            line: l,
            column: c.max(1),
            len: 1,
            source_line: src,
        });
    } else {
        d = d.with_note(format!("at {}", e.location));
    }
    vec![d.with_help("`rusolve --explain RS1006`")]
}

#[cfg(test)]
mod tests {
    use super::*;

    const BAD: &str = "\
NAME BAD
ROWS
 N OBJ
 L R1
 G R2
COLUMNS
 X OBJ 1 R1 1
 X R7 2
 Y OBJ 1 R2 1
RHS
 RHS R1 abc
 RHS R2 4
BOUNDS
 UP BND GHOST 3
ENDATA
";

    #[test]
    fn skipped_column_line_does_not_cascade() {
        let t = "NAME C\nROWS\n N OBJ\n L R1\nCOLUMNS\n X OBJ 1 R1 1\n Y R9 1\nRHS\n RHS R1 4\nBOUNDS\n UP BND Y 3\nENDATA\n";
        let ds = text_diagnostics("c.mps", t, Format::Mps);
        assert_eq!(ds.len(), 1, "{ds:?}");
        assert_eq!(ds[0].code, "RS1004");
    }

    #[test]
    fn reports_every_independent_mps_error() {
        let ds = text_diagnostics("bad.mps", BAD, Format::Mps);
        assert_eq!(ds.len(), 3, "{ds:#?}");
        assert_eq!(ds[0].code, "RS1004");
        assert!(
            ds[0].help.iter().any(|h| h.contains("R1, R2")),
            "{:?}",
            ds[0].help
        );
        match &ds[0].location {
            Some(Location::File {
                line,
                column,
                source_line,
                ..
            }) => {
                assert_eq!((*line, *column), (8, 4));
                assert_eq!(source_line, " X R7 2");
            }
            l => panic!("{l:?}"),
        }
        assert_eq!(ds[1].code, "RS1005");
        assert_eq!(ds[2].code, "RS1004");
    }

    #[test]
    fn a_valid_file_has_no_diagnostics() {
        let ok = "NAME T\nROWS\n N OBJ\n L R1\nCOLUMNS\n X OBJ 1 R1 1\nRHS\n RHS R1 4\nENDATA\n";
        assert!(text_diagnostics("ok.mps", ok, Format::Mps).is_empty());
    }

    #[test]
    fn rows_error_ends_the_list() {
        let t = "NAME T\nROWS\n Q OBJ\nCOLUMNS\n X OBJ 1\nENDATA\n";
        assert_eq!(text_diagnostics("t.mps", t, Format::Mps).len(), 1);
    }

    #[test]
    fn lp_bad_number_is_located() {
        let t = "Minimize\n obj: x + y\nSubject To\n c1: x + y >= 1e\nEnd\n";
        let ds = text_diagnostics("t.lp", t, Format::Lp);
        assert!(!ds.is_empty());
        assert!(
            matches!(ds[0].location, Some(Location::File { line: 4, .. })),
            "{ds:?}"
        );
    }
}
