//! CPLEX LP format reader (the plain-text `.lp` format, as produced by
//! CPLEX/Gurobi/HiGHS writers and used by many public benchmark instances).
//!
//! Reader covers: `Minimize`/`Maximize` (with optional objective label and
//! objective constant), `Subject To`/`st` with named and unnamed rows,
//! `L`/`G`/`E` rows and two-sided rows (`lo <= expr <= hi` and the mirrored
//! `hi >= expr >= lo`, stored as base sense + RHS + RANGE on [`Model`]),
//! `Bounds` (per-variable `free` / `<= b` / `>= b` / `= b` / `in [l,u]` /
//! two-sided chains, comma lists, and a GLPK-style `default` directive),
//! `Binary`, `General` and `Free` sections, `\` comments anywhere,
//! case-insensitive keywords, scientific notation, implicit coefficients
//! (`x`, `-x`, `3x`, `3 * x`), variables first mentioned in `Bounds` or the
//! integer sections, and line-wrapped statements (a statement continues
//! while it has no relational operator, or its line ends with an operator).
//!
//! Names are single tokens, as in the HiGHS grammar (identifiers cannot
//! contain whitespace): two names in an expression with no `+`/`-` between
//! them are an error, as are space-separated labels. Relational operators
//! accept the HiGHS spellings `<`, `<=`, `=<` and `>`, `>=`, `=>`, and `inf` /
//! `infinity` (any case, with a sign) are infinite bound and RHS values. List
//! sections (`Binary`/`General`/`Free`) take one name per token,
//! comma-separated. Names may not contain the operator characters
//! `+ - = < > * : [ ] ^ ; \`; parentheses are name characters (`x(1)`), and
//! so is a comma between two name characters (`x(1,2)`; `a, b` and `a ,b`
//! are two names), as in CPLEX and HiGHS, so MPS names such as `B0(998)`
//! (supportcase10) survive an LP round trip.
//!
//! The objective may carry a quadratic part `[ q x ^ 2 + q x * y ] / 2`
//! (CPLEX LP format), read into [`Model::hessian`] as the full square
//! matrix of `c'x + 1/2 x'Hx`: `q x ^ 2` is `H_xx = q`, `q x * y` is
//! `H_xy = H_yx = q / 2`.
//!
//! A constant on the left-hand side of a constraint (`c: x + 2 >= 5`) is
//! ignored, as HiGHS does (the row is `x >= 5`; the reader has no warning
//! channel, so it is ignored silently). An objective constant goes to
//! [`Model::obj_offset`]. An empty `semi-continuous` / `semi` section is
//! accepted. A sign must be followed by a term, two-token keywords
//! (`subject to`, `such that`, `semi-continuous`) must be on one line, and
//! nothing may follow `end`.
//!
//! Provenance: grammar follows the published IBM ILOG CPLEX LP format
//! documentation, with the `default`
//! bounds directive and the bare-expression `= 0` row convention as also
//! documented for GLPK's LP reader.
//!
//! General-constraint additions, grammar from the public
//! Gurobi "LP format" page: a `semi-continuous` section, indicator rows
//! `name: z = 1 -> row`, an `SOS` section (`name: S1 :: x : 1 y : 2`) and a
//! `General Constraints` section (`r = MAX|MIN|AND|OR|ABS|NOT ( ... )`,
//! `y = PWL ( x ) : (x0, y0) ...`). They land in [`Model::gc`].
//!
//! A constraint may carry quadratic terms `c: x + [ a x ^ 2 + b x * y ] <= r`
//! (CPLEX LP format; coefficients are monomial coefficients, no `/ 2`), read
//! into [`Model::conic`]. [`read_lp`] rejects such a model, as the LP/MIP
//! route would otherwise drop the terms; [`read_lp_conic`] keeps them. A
//! two-sided row with quadratic terms is an error.
//!
//! Not supported (explicit errors): multi-objective extensions, Gurobi
//! function constraints (POLY, EXP, ...).

use std::collections::HashMap;

use rusolve_core::Triplet;

use crate::conic::{ConicBlock, QuadRow};
use crate::gc::{FuncCon, FuncKind, GcModel, Indicator, Pwl, SemiVar, SosKind, SosSet};
use crate::model::{Model, ObjSense, RowSense, VariableBound};

/// A parse failure with the offending 1-based line number.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LpError {
    pub line: usize,
    pub message: String,
}

impl std::fmt::Display for LpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "LP parse error at line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for LpError {}

fn err<T>(line: usize, message: impl Into<String>) -> Result<T, LpError> {
    Err(LpError {
        line,
        message: message.into(),
    })
}

// ---------------------------------------------------------------- tokenizer

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Num(f64),
    Name(String),
    Le,
    Ge,
    Eq,
    Plus,
    Minus,
    Star,
    Colon,
    Comma,
    Lsq,
    Rsq,
    Caret,
    /// `->` of an indicator constraint (`z = 1 -> row`).
    Arrow,
}

fn is_rel(t: &Tok) -> bool {
    matches!(t, Tok::Le | Tok::Ge | Tok::Eq)
}

fn is_op(t: &Tok) -> bool {
    matches!(
        t,
        Tok::Plus | Tok::Minus | Tok::Star | Tok::Le | Tok::Ge | Tok::Eq | Tok::Comma
    )
}

/// Characters that may NOT appear in a variable name (operators and
/// delimiters; everything else printable is name material, so `x.1`,
/// `x'`, `p/q`, `x{1}`, `x(1)` all work — names with `-`/`+` do not, the
/// common convention shared with GLPK). A comma is a delimiter except
/// between two name characters (see [`tokenize`]).
fn is_delim(c: char) -> bool {
    c.is_whitespace()
        || matches!(
            c,
            '<' | '>' | '=' | '+' | '-' | '*' | ',' | ':' | '[' | ']' | '^' | ';' | '\\'
        )
}

/// `inf` / `infinity`, any case (HiGHS `<number>` grammar).
fn is_inf_word(s: &str) -> bool {
    s.eq_ignore_ascii_case("inf") || s.eq_ignore_ascii_case("infinity")
}

fn tokenize(text: &str) -> Result<Vec<(Tok, usize)>, LpError> {
    let b: Vec<char> = text.chars().collect();
    let mut out: Vec<(Tok, usize)> = Vec::new();
    let mut i = 0usize;
    let mut line = 1usize;
    // Inside `in [lo, hi]` a comma always separates.
    let mut in_bracket = false;
    while i < b.len() {
        let c = b[i];
        match c {
            '\n' => {
                line += 1;
                i += 1;
            }
            c if c.is_whitespace() || c == ';' => i += 1,
            // `\` comments run to end of line, anywhere in the line.
            '\\' => {
                while i < b.len() && b[i] != '\n' {
                    i += 1;
                }
            }
            // `<`, `<=`, `=<` are all Le; `>`, `>=`, `=>` are all Ge.
            '<' => {
                out.push((Tok::Le, line));
                i += if i + 1 < b.len() && b[i + 1] == '=' {
                    2
                } else {
                    1
                };
            }
            '>' => {
                out.push((Tok::Ge, line));
                i += if i + 1 < b.len() && b[i + 1] == '=' {
                    2
                } else {
                    1
                };
            }
            '=' => {
                match b.get(i + 1) {
                    Some('<') => {
                        out.push((Tok::Le, line));
                        i += 2;
                    }
                    Some('>') => {
                        out.push((Tok::Ge, line));
                        i += 2;
                    }
                    Some('=') => {
                        // `==` same as `=`
                        out.push((Tok::Eq, line));
                        i += 2;
                    }
                    _ => {
                        out.push((Tok::Eq, line));
                        i += 1;
                    }
                }
            }
            '+' => {
                out.push((Tok::Plus, line));
                i += 1;
            }
            '-' if b.get(i + 1) == Some(&'>') => {
                out.push((Tok::Arrow, line));
                i += 2;
            }
            '-' => {
                out.push((Tok::Minus, line));
                i += 1;
            }
            '*' => {
                out.push((Tok::Star, line));
                i += 1;
            }
            ':' => {
                out.push((Tok::Colon, line));
                i += 1;
            }
            ',' => {
                out.push((Tok::Comma, line));
                i += 1;
            }
            '[' => {
                in_bracket =
                    matches!(out.last(), Some((Tok::Name(n), _)) if n.eq_ignore_ascii_case("in"));
                out.push((Tok::Lsq, line));
                i += 1;
            }
            ']' => {
                in_bracket = false;
                out.push((Tok::Rsq, line));
                i += 1;
            }
            '^' => {
                out.push((Tok::Caret, line));
                i += 1;
            }
            c if c.is_ascii_digit()
                || (c == '.' && b.get(i + 1).is_some_and(|d| d.is_ascii_digit())) =>
            {
                let start = i;
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
                if i < b.len() && b[i] == '.' {
                    i += 1;
                    while i < b.len() && b[i].is_ascii_digit() {
                        i += 1;
                    }
                }
                if i < b.len() && (b[i] == 'e' || b[i] == 'E') {
                    let mut j = i + 1;
                    if j < b.len() && (b[j] == '+' || b[j] == '-') {
                        j += 1;
                    }
                    if j < b.len() && b[j].is_ascii_digit() {
                        i = j;
                        while i < b.len() && b[i].is_ascii_digit() {
                            i += 1;
                        }
                    }
                }
                let s: String = b[start..i].iter().collect();
                let v = s.parse::<f64>().map_err(|_| LpError {
                    line,
                    message: format!("bad number '{s}'"),
                })?;
                out.push((Tok::Num(v), line));
            }
            c if !is_delim(c) => {
                let start = i;
                loop {
                    while i < b.len() && !is_delim(b[i]) {
                        i += 1;
                    }
                    // `x(1,2)`: a comma between two name characters belongs
                    // to the name, as in CPLEX and HiGHS (neither treats `,`
                    // as a delimiter); `a, b` and `a ,b` stay two names.
                    if !in_bracket && i + 1 < b.len() && b[i] == ',' && !is_delim(b[i + 1]) {
                        i += 1;
                        continue;
                    }
                    break;
                }
                out.push((Tok::Name(b[start..i].iter().collect()), line));
            }
            c => return err(line, format!("unexpected character '{c}'")),
        }
    }
    Ok(out)
}

/// Model name from the `\ LP format for model <name>` marker comment
/// (written by CPLEX/HiGHS-family LP writers); empty when absent.
fn model_name(text: &str) -> String {
    for raw in text.lines() {
        let t = raw.trim_start();
        if !t.starts_with('\\') {
            if !t.is_empty() {
                break; // marker only counts while still in the comment header
            }
            continue;
        }
        let lower = t.to_ascii_lowercase();
        if let Some(k) = lower.find("lp format for model") {
            let rest = &t[k + "lp format for model".len()..];
            if let Some(w) = rest.split_whitespace().next() {
                return w.trim_matches(|c| c == '"' || c == '\'').to_string();
            }
        }
    }
    String::new()
}

// ------------------------------------------------------------------ parser

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kw {
    Min,
    Max,
    SubjectTo,
    Bounds,
    General,
    Binary,
    Free,
    End,
    Semi,
}

#[derive(Debug, Clone, Copy)]
enum BoundOp {
    Lower(f64),
    Upper(f64),
    Fixed(f64),
    Free,
    Range(f64, f64),
}

#[derive(Debug, Clone, Copy)]
enum DefaultBound {
    Free,
    Lower(f64),
    Upper(f64),
    Fixed(f64),
    Range(f64, f64),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ListKind {
    Binary,
    General,
    Free,
    /// `semi-continuous` section (a column also listed under `General` is
    /// semi-integer).
    Semi,
}

/// Quadratic monomials `(x, y, m)` by variable name.
type QuadTerms = Vec<(String, String, f64)>;

#[derive(Debug)]
struct RowAcc {
    name: String,
    terms: Vec<(String, f64)>,
    sense: RowSense,
    rhs: f64,
    range: f64,
    /// Quadratic monomials `(x, y, m)` from `[ ... ]`.
    quad: QuadTerms,
}

struct P {
    toks: Vec<(Tok, usize)>,
    pos: usize,
    sense: ObjSense,
    obj_name: Option<String>,
    obj_terms: Vec<(String, f64)>,
    obj_const: f64,
    /// Quadratic objective terms `(x, y, q)` from `[ ... ] / 2`.
    obj_quad: Vec<(String, String, f64)>,
    obj_seen: bool,
    rows: Vec<RowAcc>,
    bound_ops: Vec<(String, BoundOp)>,
    default: Option<DefaultBound>,
    lists: Vec<(String, ListKind)>,
    /// Indicator rows `(binary name, active value, row)`, in file order.
    inds: Vec<(String, bool, RowAcc)>,
    /// The `SOS` and `General Constraints` sections, parsed by
    /// [`extract_raw_sections`] (they have their own line grammar).
    raw: RawSections,
}

impl P {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos).map(|t| &t.0)
    }

    /// Line of the token at `i` (or the last token's line at EOF).
    fn line_at(&self, i: usize) -> usize {
        let i = i.min(self.toks.len().saturating_sub(1));
        self.toks.get(i).map_or(1, |t| t.1)
    }

    fn cur_line(&self) -> usize {
        self.line_at(self.pos)
    }

    /// First index past the line that token `i` lives on.
    fn line_end_from(&self, i: usize) -> usize {
        match self.toks.get(i) {
            Some(&(_, l)) => {
                let mut j = i;
                while j < self.toks.len() && self.toks[j].1 == l {
                    j += 1;
                }
                j
            }
            None => self.toks.len(),
        }
    }

    /// If a `Name :` label starts at `i`, its name and the index after the
    /// colon. A label is a single identifier (HiGHS grammar) and never
    /// crosses a line boundary.
    fn label_at(&self, i: usize) -> Option<(String, usize)> {
        let line = self.toks.get(i)?.1;
        let Tok::Name(s) = self.peek_tok(i)? else {
            return None;
        };
        if self.toks[i].1 == line && matches!(self.peek_tok(i + 1), Some(Tok::Colon)) {
            Some((s.clone(), i + 2))
        } else {
            None
        }
    }

    fn peek_tok(&self, i: usize) -> Option<&Tok> {
        self.toks.get(i).map(|t| &t.0)
    }

    /// Section keyword starting at `i` (case-insensitive), with the number
    /// of tokens it consumes. A keyword immediately followed by a (possibly
    /// space-joined) `label:` is a label, not a keyword.
    fn kw_at(&self, i: usize) -> Option<(Kw, usize)> {
        let Tok::Name(s) = self.peek_tok(i)? else {
            return None;
        };
        if self.label_at(i).is_some() {
            return None;
        }
        let lower = s.to_ascii_lowercase();
        // `semi-continuous`: three tokens on one line (HiGHS writes this).
        if lower == "semi"
            && matches!(self.peek_tok(i + 1), Some(Tok::Minus))
            && matches!(self.peek_tok(i + 2), Some(Tok::Name(c)) if c.eq_ignore_ascii_case("continuous"))
            && self.line_at(i + 1) == self.line_at(i)
            && self.line_at(i + 2) == self.line_at(i)
        {
            return Some((Kw::Semi, 3));
        }
        let kw = match lower.as_str() {
            "min" | "minimize" | "minimise" | "minimum" => Kw::Min,
            "max" | "maximize" | "maximise" | "maximum" => Kw::Max,
            "st" | "s.t." | "st." => Kw::SubjectTo,
            // `subject to` / `such that`: both tokens on one line (1.16.0
            // LpReader.cpp keywordAt).
            "subject" | "such" => match self.peek_tok(i + 1) {
                Some(Tok::Name(n))
                    if n.eq_ignore_ascii_case(if lower == "subject" { "to" } else { "that" })
                        && self.line_at(i + 1) == self.line_at(i) =>
                {
                    Kw::SubjectTo
                }
                _ => return None,
            },
            "bound" | "bounds" => Kw::Bounds,
            "gen" | "general" | "generals" | "integer" | "integers" => Kw::General,
            "bin" | "binary" | "binaries" => Kw::Binary,
            "free" | "fre" => Kw::Free,
            "end" => Kw::End,
            "semi" | "semis" | "semicontinuous" | "semiinteger" => Kw::Semi,
            _ => return None,
        };
        let n = if kw == Kw::SubjectTo && (lower == "subject" || lower == "such") {
            2
        } else {
            1
        };
        Some((kw, n))
    }

    /// Collect the statement starting at `pos` (at least to end of line,
    /// extended over wrapped continuation lines) and return its exclusive
    /// end index. A statement is complete when it ends on a plain token and
    /// either carries >= 2 relational tokens (`lo <= expr <= hi`) or exactly
    /// one relation whose RHS number terminates it. Otherwise it extends to
    /// the next line, unless that line starts a keyword or a new label
    /// (or EOF).
    fn collect_stmt(&mut self) -> usize {
        let mut end = self.line_end_from(self.pos);
        // Relation positions (relative to `self.pos`), found incrementally: a statement that
        // continues over many lines (a 1M-term objective) is scanned once, not once per line.
        let mut rels: Vec<usize> = Vec::new();
        let mut scanned = self.pos;
        loop {
            let slice = &self.toks[self.pos..end];
            rels.extend(
                self.toks[scanned..end]
                    .iter()
                    .enumerate()
                    .filter(|(_, (t, _))| is_rel(t))
                    .map(|(k, _)| scanned - self.pos + k),
            );
            scanned = end;
            let last_op = slice.last().is_some_and(|(t, _)| is_op(t));
            // `checked_sub`: a statement of fewer than two tokens has no sign before its last
            // token (a bare `len - 2` underflows there, a debug-build panic).
            let ends_in_num = matches!(slice.last(), Some((Tok::Num(_), _)))
                || matches!(
                    (
                        slice.len().checked_sub(2).and_then(|k| slice.get(k)),
                        slice.last()
                    ),
                    (Some((Tok::Plus | Tok::Minus, _)), Some((Tok::Num(_), _)))
                );
            let complete = !last_op
                && (rels.len() >= 2
                    || (rels.len() == 1 && ends_in_num && rels[0] + 1 < slice.len()));
            if complete || end >= self.toks.len() {
                return end;
            }
            if self.kw_at(end).is_some() || self.label_at(end).is_some() {
                return end;
            }
            end = self.line_end_from(end);
        }
    }

    fn parse_model(&mut self) -> Result<Model, LpError> {
        loop {
            if self.pos >= self.toks.len() {
                break;
            }
            let Some((kw, n)) = self.kw_at(self.pos) else {
                return err(
                    self.cur_line(),
                    format!(
                        "expected a section keyword, found '{}'",
                        self.describe(self.pos)
                    ),
                );
            };
            self.pos += n;
            match kw {
                Kw::Min | Kw::Max => {
                    if self.obj_seen {
                        return err(self.line_at(self.pos - n), "second objective section");
                    }
                    self.obj_seen = true;
                    self.sense = if kw == Kw::Min {
                        ObjSense::Minimize
                    } else {
                        ObjSense::Maximize
                    };
                    self.parse_obj()?;
                }
                Kw::SubjectTo => self.parse_constraints()?,
                Kw::Bounds => self.parse_bounds()?,
                Kw::General => self.parse_list(ListKind::General)?,
                Kw::Binary => self.parse_list(ListKind::Binary)?,
                Kw::Free => self.parse_list(ListKind::Free)?,
                // Names listed here are semi-continuous (semi-integer when
                // also listed under `General`); an empty section is accepted.
                Kw::Semi => self.parse_list(ListKind::Semi)?,
                Kw::End => {
                    // Nothing may follow `end` (1.16.0 LpReader.cpp parse).
                    if self.pos < self.toks.len() {
                        return err(
                            self.cur_line(),
                            format!(
                                "expected end of file, found '{}' after `end`",
                                self.describe(self.pos)
                            ),
                        );
                    }
                    break;
                }
            }
        }
        self.build()
    }

    fn describe(&self, i: usize) -> String {
        match self.peek_tok(i) {
            Some(Tok::Name(s)) => s.clone(),
            Some(Tok::Num(v)) => format!("{v}"),
            Some(_) => "operator".to_string(),
            None => "end of file".to_string(),
        }
    }

    /// One `Name` token, never past `end`; a `stop` word (case-insensitive)
    /// is not a name. Names are single tokens (HiGHS grammar), so a second
    /// name in the same statement is left for the caller to reject.
    fn name_token(&mut self, stop: &[&str], end: usize) -> Result<String, LpError> {
        let line = self.cur_line();
        if self.pos < end
            && let Some(Tok::Name(s)) = self.peek().cloned()
            && !stop.iter().any(|w| s.eq_ignore_ascii_case(w))
        {
            self.pos += 1;
            return Ok(s);
        }
        err(
            line,
            format!(
                "expected a variable name, found '{}'",
                self.describe(self.pos)
            ),
        )
    }

    /// Linear expression: signed terms `[+-] [num] [*] name`, implicit
    /// coefficients (`3x`, `3 x`, `3 * x`, `x`, `-x`) and plain constants.
    /// Returns the terms and the accumulated constant.
    /// C++ `Parser::parseTerm`/`LpReader.cpp` (1.16.0): every term after the
    /// first must be preceded by a `+` or `-`, so `x y` and `3 x 4 y` are
    /// errors rather than sums.
    fn parse_expr(&mut self, end: usize) -> Result<(Vec<(String, f64)>, f64), LpError> {
        let mut terms: Vec<(String, f64)> = Vec::new();
        let mut konst = 0.0;
        let mut first = true;
        loop {
            let mut sign = 1.0f64;
            let mut signed = false;
            let sign_start = self.pos;
            while self.pos < end && matches!(self.peek(), Some(Tok::Plus | Tok::Minus)) {
                if matches!(self.peek(), Some(Tok::Minus)) {
                    sign = -sign;
                }
                signed = true;
                self.pos += 1;
            }
            // `+ [ ... ] / 2`: the quadratic part belongs to the caller
            // (objective only).
            if self.pos < end && matches!(self.peek(), Some(Tok::Lsq)) {
                self.pos = sign_start;
                break;
            }
            // A sign must be followed by a term (1.16.0 parseTerm: "expected
            // a term"), e.g. `x + >= 1` or a dangling `x +` before a keyword.
            if signed
                && (self.pos >= end || !matches!(self.peek(), Some(Tok::Num(_) | Tok::Name(_))))
            {
                return err(
                    self.cur_line(),
                    format!("expected a term, found '{}'", self.describe(self.pos)),
                );
            }
            if self.pos >= end {
                break;
            }
            if !first && !signed && matches!(self.peek(), Some(Tok::Num(_) | Tok::Name(_))) {
                return err(
                    self.cur_line(),
                    format!("expected `+` or `-`, found '{}'", self.describe(self.pos)),
                );
            }
            first = false;
            match self.peek() {
                Some(Tok::Num(_)) => {
                    let Some(Tok::Num(v)) = self.peek().cloned() else {
                        unreachable!()
                    };
                    self.pos += 1;
                    let star = matches!(self.peek(), Some(Tok::Star));
                    if star {
                        self.pos += 1;
                    }
                    let is_name = self.pos < end && matches!(self.peek(), Some(Tok::Name(_)));
                    if is_name || star {
                        let name = self.name_token(&[], end)?;
                        terms.push((name, sign * v));
                    } else {
                        konst += sign * v;
                    }
                }
                Some(Tok::Name(_)) => {
                    let name = self.name_token(&[], end)?;
                    terms.push((name, sign));
                }
                other => {
                    let t = other.cloned();
                    // Trailing relation/bracket/colon tokens belong to the caller
                    // (relation + RHS, or the `label:` of a following row on
                    // the same line); anything else is junk in the expression.
                    match t {
                        Some(Tok::Le) | Some(Tok::Ge) | Some(Tok::Eq) | Some(Tok::Colon) => break,
                        _ => {
                            return err(
                                self.cur_line(),
                                format!("unexpected '{}' in expression", self.describe(self.pos)),
                            );
                        }
                    }
                }
            }
        }
        Ok((terms, konst))
    }

    /// Value of the number token at `i`: a constant, or `inf` / `infinity`
    /// (any case) as `+inf`. Sign handling is the caller's.
    fn num_value_at(&self, i: usize) -> Option<f64> {
        match self.peek_tok(i)? {
            Tok::Num(v) => Some(*v),
            Tok::Name(s) if is_inf_word(s) => Some(f64::INFINITY),
            _ => None,
        }
    }

    /// RHS number with optional sign; leaves `pos` after it.
    /// `<number> := [<signs>] (<constant> | "inf" | "infinity")`.
    fn parse_num(&mut self, end: usize) -> Result<f64, LpError> {
        let mut sign = 1.0;
        while self.pos < end && matches!(self.peek(), Some(Tok::Plus | Tok::Minus)) {
            if matches!(self.peek(), Some(Tok::Minus)) {
                sign = -sign;
            }
            self.pos += 1;
        }
        match self.num_value_at(self.pos) {
            Some(v) => {
                self.pos += 1;
                Ok(sign * v)
            }
            _ => err(
                self.cur_line(),
                format!("expected a number, found '{}'", self.describe(self.pos)),
            ),
        }
    }

    /// Consume the next token if it is a relational operator.
    fn take_rel(&mut self) -> Option<Tok> {
        let t = match self.peek() {
            Some(Tok::Le) => Tok::Le,
            Some(Tok::Ge) => Tok::Ge,
            Some(Tok::Eq) => Tok::Eq,
            _ => return None,
        };
        self.pos += 1;
        Some(t)
    }

    fn parse_obj(&mut self) -> Result<(), LpError> {
        let end = self.collect_stmt();
        if let Some((label, after)) = self.label_at(self.pos) {
            self.obj_name = Some(label);
            self.pos = after;
        }
        let (mut terms, mut konst) = self.parse_expr(end)?;
        if self.pos < end {
            self.parse_quadratic(end, &mut terms)?;
            let (more, k2) = self.parse_expr(end)?;
            terms.extend(more);
            konst += k2;
        }
        // An objective relation (`obj: expr = 10`) is not LP format; it used
        // to be skipped silently, now it is an error (1.16.0 LpReader.cpp).
        if matches!(self.peek(), Some(Tok::Le | Tok::Ge | Tok::Eq)) && self.pos < end {
            return err(self.cur_line(), "relation in objective");
        }
        if self.pos != end {
            return err(
                self.cur_line(),
                format!("trailing '{}' in objective", self.describe(self.pos)),
            );
        }
        self.obj_terms = terms;
        self.obj_const = konst;
        Ok(())
    }

    /// `[+-] [ q x ^ 2 + q x * y ... ] / 2` in the objective, if one starts
    /// at `pos`: `q x ^ 2` is the Hessian entry `H_xx = q` and `q x * y` is
    /// `H_xy = H_yx = q / 2` (objective `c'x + 1/2 x'Hx`, as the CPLEX LP
    /// format defines it). Each variable is also pushed as
    /// a zero linear term so columns keep their order of first appearance.
    fn parse_quadratic(
        &mut self,
        end: usize,
        terms: &mut Vec<(String, f64)>,
    ) -> Result<(), LpError> {
        let save = self.pos;
        let mut outer = 1.0f64;
        while self.pos < end && matches!(self.peek(), Some(Tok::Plus | Tok::Minus)) {
            if matches!(self.peek(), Some(Tok::Minus)) {
                outer = -outer;
            }
            self.pos += 1;
        }
        if !matches!(self.peek(), Some(Tok::Lsq)) || self.pos >= end {
            self.pos = save;
            return Ok(());
        }
        self.pos += 1;
        let mut first = true;
        loop {
            if self.pos >= end {
                return err(self.cur_line(), "unterminated `[` in objective");
            }
            if matches!(self.peek(), Some(Tok::Rsq)) {
                self.pos += 1;
                break;
            }
            let mut sign = outer;
            let mut signed = false;
            while self.pos < end && matches!(self.peek(), Some(Tok::Plus | Tok::Minus)) {
                if matches!(self.peek(), Some(Tok::Minus)) {
                    sign = -sign;
                }
                signed = true;
                self.pos += 1;
            }
            if !first && !signed {
                return err(
                    self.cur_line(),
                    format!("expected `+` or `-`, found '{}'", self.describe(self.pos)),
                );
            }
            first = false;
            let mut q = 1.0;
            if let Some(Tok::Num(v)) = self.peek().cloned() {
                q = v;
                self.pos += 1;
            }
            let a = self.name_token(&[], end)?;
            let b = match self.peek() {
                Some(Tok::Caret) => {
                    self.pos += 1;
                    if !matches!(self.peek(), Some(Tok::Num(v)) if *v == 2.0) {
                        return err(self.cur_line(), "expected `^ 2` in quadratic term");
                    }
                    self.pos += 1;
                    a.clone()
                }
                Some(Tok::Star) => {
                    self.pos += 1;
                    self.name_token(&[], end)?
                }
                _ => {
                    return err(
                        self.cur_line(),
                        format!(
                            "expected `^ 2` or `* name` in quadratic term, found '{}'",
                            self.describe(self.pos)
                        ),
                    );
                }
            };
            terms.push((a.clone(), 0.0));
            if b != a {
                terms.push((b.clone(), 0.0));
            }
            self.obj_quad.push((a, b, sign * q));
        }
        // `/ 2`: `/` is a name character, so it is `/` + `2` or `/2`.
        let halved = match self.peek().cloned() {
            Some(Tok::Name(n)) if n == "/" && self.pos + 1 < end => {
                self.pos += 1;
                matches!(self.peek(), Some(Tok::Num(v)) if *v == 2.0) && {
                    self.pos += 1;
                    true
                }
            }
            Some(Tok::Name(n)) if n == "/2" => {
                self.pos += 1;
                true
            }
            _ => false,
        };
        if !halved {
            return err(
                self.cur_line(),
                "expected `/ 2` after the quadratic objective",
            );
        }
        Ok(())
    }

    /// A `[` (after optional signs) starts at `pos`.
    fn at_bracket(&self, end: usize) -> bool {
        let mut k = self.pos;
        while k < end && matches!(self.toks.get(k).map(|t| &t.0), Some(Tok::Plus | Tok::Minus)) {
            k += 1;
        }
        k < end && matches!(self.toks.get(k).map(|t| &t.0), Some(Tok::Lsq))
    }

    /// `[+-] [ m x ^ 2 + m x * y ... ]` in a constraint: monomial
    /// coefficients, no `/ 2`.
    fn parse_row_quadratic(
        &mut self,
        end: usize,
        quad: &mut Vec<(String, String, f64)>,
    ) -> Result<(), LpError> {
        let mut outer = 1.0f64;
        while self.pos < end && matches!(self.peek(), Some(Tok::Plus | Tok::Minus)) {
            if matches!(self.peek(), Some(Tok::Minus)) {
                outer = -outer;
            }
            self.pos += 1;
        }
        self.pos += 1; // `[`
        let mut first = true;
        loop {
            if self.pos >= end {
                return err(self.cur_line(), "unterminated `[` in constraint");
            }
            if matches!(self.peek(), Some(Tok::Rsq)) {
                self.pos += 1;
                return Ok(());
            }
            let mut sign = outer;
            let mut signed = false;
            while self.pos < end && matches!(self.peek(), Some(Tok::Plus | Tok::Minus)) {
                if matches!(self.peek(), Some(Tok::Minus)) {
                    sign = -sign;
                }
                signed = true;
                self.pos += 1;
            }
            if !first && !signed {
                return err(
                    self.cur_line(),
                    format!("expected `+` or `-`, found '{}'", self.describe(self.pos)),
                );
            }
            first = false;
            let mut q = 1.0;
            if let Some(Tok::Num(v)) = self.peek().cloned() {
                q = v;
                self.pos += 1;
                if matches!(self.peek(), Some(Tok::Star)) && self.pos + 1 < end {
                    // `3 * x ^ 2`
                    self.pos += 1;
                }
            }
            let a = self.name_token(&[], end)?;
            let b = match self.peek() {
                Some(Tok::Caret) => {
                    self.pos += 1;
                    if !matches!(self.peek(), Some(Tok::Num(v)) if *v == 2.0) {
                        return err(self.cur_line(), "expected `^ 2` in quadratic term");
                    }
                    self.pos += 1;
                    a.clone()
                }
                Some(Tok::Star) => {
                    self.pos += 1;
                    self.name_token(&[], end)?
                }
                _ => {
                    return err(
                        self.cur_line(),
                        format!(
                            "expected `^ 2` or `* name` in quadratic term, found '{}'",
                            self.describe(self.pos)
                        ),
                    );
                }
            };
            quad.push((a, b, sign * q));
        }
    }

    fn parse_constraints(&mut self) -> Result<(), LpError> {
        loop {
            if self.pos >= self.toks.len() {
                return Ok(());
            }
            if self.kw_at(self.pos).is_some() {
                return Ok(());
            }
            let end = self.collect_stmt();
            self.parse_row(end)?;
            if self.pos != end {
                return err(
                    self.cur_line(),
                    format!("trailing '{}' in constraint row", self.describe(self.pos)),
                );
            }
        }
    }

    fn parse_row(&mut self, end: usize) -> Result<(), LpError> {
        let name = match self.label_at(self.pos) {
            Some((label, after)) => {
                self.pos = after;
                label
            }
            None => String::new(),
        };
        // Indicator: `[label:] z = 1 -> row` (CPLEX and Gurobi LP formats).
        if let (Some(Tok::Name(z)), Some(Tok::Eq), Some(Tok::Num(v)), Some(Tok::Arrow)) = (
            self.peek_tok(self.pos).cloned(),
            self.peek_tok(self.pos + 1),
            self.peek_tok(self.pos + 2),
            self.peek_tok(self.pos + 3),
        ) {
            let line = self.cur_line();
            let v = *v;
            if v != 0.0 && v != 1.0 {
                return err(line, format!("indicator value must be 0 or 1, got {v}"));
            }
            self.pos += 4;
            let before = self.rows.len();
            self.parse_row_body(end, name)?;
            if self.rows.len() != before + 1 {
                return err(line, "indicator constraint needs a row after `->`");
            }
            let row = self.rows.pop().expect("one row was just pushed");
            if !row.quad.is_empty() {
                return err(line, "an indicator row cannot have quadratic terms");
            }
            if row.range != 0.0 {
                return err(line, "an indicator row cannot be two-sided");
            }
            self.inds.push((z, v == 1.0, row));
            return Ok(());
        }
        self.parse_row_body(end, name)
    }

    /// The row after an optional label and indicator prefix.
    fn parse_row_body(&mut self, end: usize, name: String) -> Result<(), LpError> {
        // Two-sided form: `[sign] num (<=|>=) expr [(<=|>=) num]`.
        let save = self.pos;
        let mut sign = 1.0;
        while self.pos < end && matches!(self.peek(), Some(Tok::Plus | Tok::Minus)) {
            if matches!(self.peek(), Some(Tok::Minus)) {
                sign = -sign;
            }
            self.pos += 1;
        }
        let mut two_sided = None;
        if self.pos < end
            && let Some(Tok::Num(v)) = self.peek().cloned()
        {
            let after_num = self.pos + 1;
            let rel = if after_num < end {
                self.toks.get(after_num).map(|t| &t.0)
            } else {
                None
            };
            if matches!(rel, Some(Tok::Le) | Some(Tok::Ge)) {
                let left = sign * v;
                let r1 = if matches!(rel, Some(Tok::Le)) {
                    Tok::Le
                } else {
                    Tok::Ge
                };
                self.pos = after_num + 1;
                two_sided = Some((left, r1));
            }
        }
        // A constant inside the expression is ignored (not moved to the RHS),
        // as HiGHS does (LpReader.cpp parseTerm, 1.16.0).
        if let Some((left, r1)) = two_sided {
            let (terms, _ignored_konst) = self.parse_expr(end)?;
            if self.at_bracket(end) {
                return err(
                    self.cur_line(),
                    "quadratic terms in a two-sided row are not supported",
                );
            }
            let r2 = if self.pos < end {
                self.take_rel()
            } else {
                None
            };
            match r2 {
                None => {
                    // `lo <= expr` (G) / `hi >= expr` (L), one-sided.
                    let sense = if r1 == Tok::Le {
                        RowSense::GreaterEqual
                    } else {
                        RowSense::LessEqual
                    };
                    self.rows.push(RowAcc {
                        name,
                        terms,
                        sense,
                        rhs: left,
                        range: 0.0,
                        quad: Vec::new(),
                    });
                }
                Some(r2) if r2 == r1 && !matches!(r1, Tok::Eq) => {
                    let right = self.parse_num(end)?;
                    // Two-sided row as base sense + RHS + RANGE (the same
                    // representation MPS RANGES produces):
                    //   `lo <= e <= hi` ⇒ L hi | hi−lo ; `hi >= e >= lo` same.
                    let (rhs, span) = if r1 == Tok::Le {
                        (right, right - left)
                    } else {
                        (left, left - right)
                    };
                    if span < 0.0 {
                        return err(
                            self.line_at(self.pos),
                            format!("inverted two-sided row bounds ({left} / {right})"),
                        );
                    }
                    self.rows.push(RowAcc {
                        name,
                        terms,
                        sense: RowSense::LessEqual,
                        rhs,
                        range: span,
                        quad: Vec::new(),
                    });
                }
                Some(other) => {
                    let _ = other;
                    return err(
                        self.line_at(self.pos),
                        "mixed relational operators in a row (e.g. `a <= expr >= b`)",
                    );
                }
            }
            return Ok(());
        }
        self.pos = save;
        let (mut terms, mut konst) = self.parse_expr(end)?;
        let mut quad: Vec<(String, String, f64)> = Vec::new();
        while self.at_bracket(end) {
            self.parse_row_quadratic(end, &mut quad)?;
            let (more, k2) = self.parse_expr(end)?;
            terms.extend(more);
            konst += k2;
        }
        let (sense, rhs) = match if self.pos < end {
            self.take_rel()
        } else {
            None
        } {
            Some(r) => {
                let rhs = self.parse_num(end)?;
                let s = match r {
                    Tok::Le => RowSense::LessEqual,
                    Tok::Ge => RowSense::GreaterEqual,
                    _ => RowSense::Equal,
                };
                (s, rhs)
            }
            None => {
                // Bare expression ⇒ `= 0` (GLPK LP-reader convention).
                if konst != 0.0 {
                    return err(
                        self.line_at(self.pos),
                        "a constant-only row needs a relational operator",
                    );
                }
                (RowSense::Equal, 0.0)
            }
        };
        self.rows.push(RowAcc {
            name,
            terms,
            sense,
            rhs,
            range: 0.0,
            quad,
        });
        Ok(())
    }

    fn parse_bounds(&mut self) -> Result<(), LpError> {
        while self.pos < self.toks.len() && self.kw_at(self.pos).is_none() {
            let end = self.line_end_from(self.pos);
            self.parse_bound_line(end)?;
            if self.pos != end {
                return err(
                    self.cur_line(),
                    format!("trailing '{}' in Bounds line", self.describe(self.pos)),
                );
            }
        }
        Ok(())
    }

    /// Comma-separated variable list (each entry space-joined per the
    /// name convention), bounded by `end`.
    fn name_list(&mut self, end: usize) -> Result<Vec<String>, LpError> {
        let mut names = vec![self.name_token(&["free", "in", "default"], end)?];
        while self.pos < end && matches!(self.peek(), Some(Tok::Comma)) {
            self.pos += 1;
            names.push(self.name_token(&["free", "in", "default"], end)?);
        }
        Ok(names)
    }

    fn push_ops(&mut self, names: &[String], op: BoundOp) {
        for n in names {
            self.bound_ops.push((n.clone(), op));
        }
    }

    fn parse_bound_line(&mut self, end: usize) -> Result<(), LpError> {
        if let Some(Tok::Name(s)) = self.peek().cloned()
            && s.eq_ignore_ascii_case("default")
        {
            self.pos += 1;
            let dflt = self.parse_default_spec(end)?;
            self.default = Some(dflt);
            return Ok(());
        }
        // Left-constant form: `[sign] num (<=|>=) list [(<=|>=) num]`.
        let save = self.pos;
        let mut sign = 1.0;
        while self.pos < end && matches!(self.peek(), Some(Tok::Plus | Tok::Minus)) {
            if matches!(self.peek(), Some(Tok::Minus)) {
                sign = -sign;
            }
            self.pos += 1;
        }
        if self.pos < end
            && let Some(v) = self.num_value_at(self.pos)
        {
            let rel = if self.pos + 1 < end {
                self.toks.get(self.pos + 1).map(|t| &t.0)
            } else {
                None
            };
            // `<number> <comparison> <identifier> [<comparison> <number>]`
            // (1.16.0 LpReader.cpp parseBound): `5 = x` fixes x.
            if matches!(rel, Some(Tok::Le) | Some(Tok::Ge) | Some(Tok::Eq)) {
                let left = sign * v;
                let r1 = match rel {
                    Some(Tok::Le) => Tok::Le,
                    Some(Tok::Ge) => Tok::Ge,
                    _ => Tok::Eq,
                };
                self.pos += 2;
                let names = self.name_list(end)?;
                let r2 = if self.pos < end {
                    self.take_rel()
                } else {
                    None
                };
                if let Some(r2) = r2 {
                    if r1 == Tok::Eq {
                        return err(
                            self.line_at(self.pos),
                            "a bound on both sides of a variable cannot use `=`",
                        );
                    }
                    let right = self.parse_num(end)?;
                    if r2 != r1 {
                        return err(
                            self.line_at(self.pos),
                            "mixed relational operators in a bound chain",
                        );
                    }
                    let (lo, hi) = if r1 == Tok::Le {
                        (left, right)
                    } else {
                        (right, left)
                    };
                    if lo > hi {
                        return err(
                            self.line_at(self.pos),
                            format!("inverted bound range ({lo} / {hi})"),
                        );
                    }
                    self.push_ops(&names, BoundOp::Range(lo, hi));
                } else if r1 == Tok::Le {
                    self.push_ops(&names, BoundOp::Lower(left));
                } else if r1 == Tok::Ge {
                    self.push_ops(&names, BoundOp::Upper(left));
                } else {
                    self.push_ops(&names, BoundOp::Fixed(left));
                }
                return Ok(());
            }
        }
        self.pos = save;
        let names = self.name_list(end)?;
        match self.peek().cloned() {
            Some(Tok::Name(s)) if s.eq_ignore_ascii_case("free") && self.pos < end => {
                self.pos += 1;
                self.push_ops(&names, BoundOp::Free);
            }
            Some(Tok::Name(s)) if s.eq_ignore_ascii_case("in") && self.pos < end => {
                self.pos += 1;
                let (lo, hi) = self.parse_in_brackets(end)?;
                self.push_ops(&names, BoundOp::Range(lo, hi));
            }
            Some(Tok::Le) if self.pos < end => {
                self.pos += 1;
                let v = self.parse_num(end)?;
                self.push_ops(&names, BoundOp::Upper(v));
            }
            Some(Tok::Ge) if self.pos < end => {
                self.pos += 1;
                let v = self.parse_num(end)?;
                self.push_ops(&names, BoundOp::Lower(v));
            }
            Some(Tok::Eq) if self.pos < end => {
                self.pos += 1;
                let v = self.parse_num(end)?;
                self.push_ops(&names, BoundOp::Fixed(v));
            }
            _ => {
                return err(
                    self.cur_line(),
                    "bound needs a directive: free, <= b, >= b, = b, in [l, u] or a range chain",
                );
            }
        }
        Ok(())
    }

    fn parse_in_brackets(&mut self, end: usize) -> Result<(f64, f64), LpError> {
        if !matches!(self.peek(), Some(Tok::Lsq)) {
            return err(self.cur_line(), "expected '[' after 'in'");
        }
        self.pos += 1;
        let lo = self.parse_num(end)?;
        if !matches!(self.peek(), Some(Tok::Comma)) {
            return err(self.cur_line(), "expected ',' in `in [lo, hi]`");
        }
        self.pos += 1;
        let hi = self.parse_num(end)?;
        if !matches!(self.peek(), Some(Tok::Rsq)) {
            return err(self.cur_line(), "expected ']' in `in [lo, hi]`");
        }
        self.pos += 1;
        if lo > hi {
            return err(
                self.cur_line(),
                format!("inverted bound range ({lo} / {hi})"),
            );
        }
        Ok((lo, hi))
    }

    fn parse_default_spec(&mut self, end: usize) -> Result<DefaultBound, LpError> {
        match self.peek().cloned() {
            Some(Tok::Name(s)) if s.eq_ignore_ascii_case("free") && self.pos < end => {
                self.pos += 1;
                Ok(DefaultBound::Free)
            }
            Some(Tok::Name(s)) if s.eq_ignore_ascii_case("in") && self.pos < end => {
                self.pos += 1;
                let (lo, hi) = self.parse_in_brackets(end)?;
                Ok(DefaultBound::Range(lo, hi))
            }
            Some(Tok::Le) if self.pos < end => {
                self.pos += 1;
                Ok(DefaultBound::Upper(self.parse_num(end)?))
            }
            Some(Tok::Ge) if self.pos < end => {
                self.pos += 1;
                Ok(DefaultBound::Lower(self.parse_num(end)?))
            }
            Some(Tok::Eq) if self.pos < end => {
                self.pos += 1;
                Ok(DefaultBound::Fixed(self.parse_num(end)?))
            }
            _ => err(
                self.cur_line(),
                "default bound needs: free, <= b, >= b, = b or in [l, u]",
            ),
        }
    }

    fn parse_list(&mut self, kind: ListKind) -> Result<(), LpError> {
        while self.pos < self.toks.len() && self.kw_at(self.pos).is_none() {
            let end = self.line_end_from(self.pos);
            while self.pos < end {
                match self.peek().cloned() {
                    Some(Tok::Name(s)) => {
                        self.lists.push((s, kind));
                        self.pos += 1;
                    }
                    Some(Tok::Comma) => self.pos += 1,
                    _ => {
                        return err(
                            self.cur_line(),
                            format!(
                                "expected a variable name, found '{}'",
                                self.describe(self.pos)
                            ),
                        );
                    }
                }
            }
        }
        Ok(())
    }

    fn build(&mut self) -> Result<Model, LpError> {
        // Column order: first appearance in the objective, then rows, then
        // Bounds, then the Binary/General/Free lists.
        let mut col_index: HashMap<String, usize> = HashMap::new();
        let mut col_names: Vec<String> = Vec::new();
        let col = |names: &mut Vec<String>, index: &mut HashMap<String, usize>, n: &str| -> usize {
            if let Some(&j) = index.get(n) {
                return j;
            }
            let j = names.len();
            index.insert(n.to_string(), j);
            names.push(n.to_string());
            j
        };

        let mut objective;
        {
            let mut obj = vec![0.0f64; col_names.len()];
            for (n, c) in &self.obj_terms {
                let j = col(&mut col_names, &mut col_index, n);
                if obj.len() < col_names.len() {
                    obj.resize(col_names.len(), 0.0);
                }
                obj[j] += c;
            }
            obj.resize(col_names.len(), 0.0);
            objective = obj;
        }

        let mut triplets: Vec<Triplet> = Vec::new();
        let mut senses = Vec::new();
        let mut rhs = Vec::new();
        let mut ranges = Vec::new();
        let mut row_names = Vec::new();
        let mut row_quads: Vec<(usize, &QuadTerms)> = Vec::new();
        for r in &self.rows {
            let i = senses.len();
            for (n, c) in &r.terms {
                let j = col(&mut col_names, &mut col_index, n);
                triplets.push((i, j, *c));
            }
            senses.push(r.sense);
            rhs.push(r.rhs);
            ranges.push(r.range);
            row_names.push(r.name.clone());
            if !r.quad.is_empty() {
                // quadratic variables become columns without a matrix entry
                for (a, b, _) in &r.quad {
                    col(&mut col_names, &mut col_index, a);
                    col(&mut col_names, &mut col_index, b);
                }
                row_quads.push((i, &r.quad));
            }
        }
        let nc = senses.len();
        objective.resize(col_names.len(), 0.0);
        // Indicator rows. Their columns come after the regular rows'.
        let mut gc = GcModel::default();
        for (k, (z, val, r)) in self.inds.iter().enumerate() {
            let bin = col(&mut col_names, &mut col_index, z);
            let terms = r
                .terms
                .iter()
                .map(|(n, c)| (col(&mut col_names, &mut col_index, n), *c))
                .collect();
            gc.indicators.push(Indicator {
                name: if r.name.is_empty() {
                    format!("ind{}", k + 1)
                } else {
                    r.name.clone()
                },
                bin,
                bin_val: *val,
                terms,
                sense: r.sense,
                rhs: r.rhs,
            });
        }
        objective.resize(col_names.len(), 0.0);

        let mut bounds = vec![VariableBound::default(); col_names.len()];
        let mut explicit = vec![false; col_names.len()];
        for (n, op) in self.bound_ops.clone() {
            let j = col(&mut col_names, &mut col_index, &n);
            bounds.resize(col_names.len(), VariableBound::default());
            explicit.resize(col_names.len(), false);
            let b = &mut bounds[j];
            explicit[j] = true;
            match op {
                BoundOp::Lower(v) => b.lower = v,
                BoundOp::Upper(v) => b.upper = v,
                BoundOp::Fixed(v) => {
                    b.lower = v;
                    b.upper = v;
                }
                BoundOp::Free => {
                    b.lower = f64::NEG_INFINITY;
                    b.upper = f64::INFINITY;
                }
                BoundOp::Range(lo, hi) => {
                    b.lower = lo;
                    b.upper = hi;
                }
            }
        }
        if let Some(d) = self.default {
            for (b, &exp) in bounds.iter_mut().zip(&explicit) {
                if exp {
                    continue;
                }
                match d {
                    DefaultBound::Free => {
                        b.lower = f64::NEG_INFINITY;
                        b.upper = f64::INFINITY;
                    }
                    DefaultBound::Lower(v) => b.lower = v,
                    DefaultBound::Upper(v) => b.upper = v,
                    DefaultBound::Fixed(v) => {
                        b.lower = v;
                        b.upper = v;
                    }
                    DefaultBound::Range(lo, hi) => {
                        b.lower = lo;
                        b.upper = hi;
                    }
                }
            }
        }
        let mut semi_cols: Vec<usize> = Vec::new();
        for (n, kind) in self.lists.clone() {
            let j = col(&mut col_names, &mut col_index, &n);
            bounds.resize(col_names.len(), VariableBound::default());
            let b = &mut bounds[j];
            match kind {
                ListKind::Semi => semi_cols.push(j),
                ListKind::Binary => {
                    b.lower = 0.0;
                    b.upper = 1.0;
                    b.integer = true;
                }
                ListKind::General => b.integer = true,
                ListKind::Free => {
                    b.lower = f64::NEG_INFINITY;
                    b.upper = f64::INFINITY;
                }
            }
        }
        // A column in the semi-continuous section is semi-integer when it
        // is also integer (General, or Binary with its [0, 1] box).
        for &j in &semi_cols {
            if gc.semi.iter().any(|s| s.col == j) {
                continue;
            }
            let integer = bounds[j].integer;
            bounds[j].integer = false;
            gc.semi.push(SemiVar { col: j, integer });
        }
        // SOS sets and general constraints.
        for s in &self.raw.sos {
            let mut members = Vec::with_capacity(s.members.len());
            for (i, (n, w)) in s.members.iter().enumerate() {
                members.push((
                    col(&mut col_names, &mut col_index, n),
                    w.unwrap_or(i as f64 + 1.0),
                ));
            }
            gc.sos.push(SosSet {
                name: s.name.clone(),
                kind: s.kind,
                priority: s.priority,
                members,
            });
        }
        for (k, g) in self.raw.gens.iter().enumerate() {
            let res = col(&mut col_names, &mut col_index, &g.res);
            let name = if g.name.is_empty() {
                format!("gc{}", k + 1)
            } else {
                g.name.clone()
            };
            match &g.what {
                RawGenKind::Pwl(x, pts) => {
                    let x = col(&mut col_names, &mut col_index, x);
                    gc.pwl.push(Pwl {
                        name,
                        x,
                        y: res,
                        xs: pts.iter().map(|p| p.0).collect(),
                        ys: pts.iter().map(|p| p.1).collect(),
                    });
                }
                RawGenKind::Func(kind, args) => {
                    let mut cols = Vec::new();
                    let mut constant: Option<f64> = None;
                    for a in args {
                        match a {
                            RawArg::Var(n) => cols.push(col(&mut col_names, &mut col_index, n)),
                            RawArg::Const(v) => {
                                // several constants fold into one
                                constant = Some(match (constant, kind) {
                                    (None, _) => *v,
                                    (Some(c), FuncKind::Min) => c.min(*v),
                                    (Some(c), _) => c.max(*v),
                                });
                            }
                        }
                    }
                    if constant.is_some() && !matches!(kind, FuncKind::Min | FuncKind::Max) {
                        return err(g.line, format!("{} takes variables only", kind.word()));
                    }
                    gc.funcs.push(FuncCon {
                        name,
                        kind: *kind,
                        res,
                        args: cols,
                        constant,
                    });
                }
            }
        }
        bounds.resize(col_names.len(), VariableBound::default());
        let nv = col_names.len();
        objective.resize(nv, 0.0);
        gc.validate(&bounds, &col_names)
            .map_err(|message| LpError { line: 0, message })?;
        let mut hessian: Vec<Triplet> = Vec::new();
        for (a, b, q) in &self.obj_quad {
            if *q == 0.0 {
                continue; // zero entries are skipped, as in the MPS reader
            }
            let (i, j) = (col_index[a], col_index[b]);
            if i == j {
                hessian.push((i, i, *q));
            } else {
                hessian.push((i, j, q / 2.0));
                hessian.push((j, i, q / 2.0));
            }
        }
        let matrix = rusolve_core::CscMatrix::from_triplets(nc, nv, &triplets);
        let quad_rows: Vec<QuadRow> = row_quads
            .into_iter()
            .map(|(i, qs)| {
                QuadRow::from_monomials(
                    i,
                    qs.iter().map(|(a, b, m)| (col_index[a], col_index[b], *m)),
                )
            })
            .filter(|q| !q.q.is_empty())
            .collect();

        Ok(Model {
            name: String::new(), // filled by read_lp from the marker comment
            num_variables: nv,
            num_constraints: nc,
            obj_sense: self.sense,
            obj_name: self.obj_name.clone(),
            obj_offset: self.obj_const,
            objective,
            senses,
            rhs,
            ranges,
            bounds,
            matrix,
            row_names,
            col_names,
            hessian,
            gc,
            conic: ConicBlock {
                quad_rows,
                cones: Vec::new(),
            },
        })
    }
}

// ----------------------------------------------- SOS and General Constraints

/// One `SOS` line: `name: S1|S2 [priority] :: var : weight ...`.
#[derive(Debug, Clone)]
struct RawSos {
    name: String,
    kind: SosKind,
    priority: f64,
    members: Vec<(String, Option<f64>)>,
}

/// One `General Constraints` line.
#[derive(Debug, Clone)]
struct RawGen {
    line: usize,
    name: String,
    res: String,
    what: RawGenKind,
}

#[derive(Debug, Clone)]
enum RawGenKind {
    /// `MAX`/`MIN`/`AND`/`OR`/`ABS`/`NOT` with variable and numeric arguments.
    Func(FuncKind, Vec<RawArg>),
    /// `PWL ( x ) : (x0, y0) ...`
    Pwl(String, Vec<(f64, f64)>),
}

#[derive(Debug, Clone)]
enum RawArg {
    Var(String),
    Const(f64),
}

/// The `SOS` and `General Constraints` sections, pulled out of the text
/// before tokenizing (their grammar uses `::`, `(`, `)` and `,` which the
/// shared tokenizer treats as name characters).
#[derive(Debug, Clone, Default)]
struct RawSections {
    sos: Vec<RawSos>,
    gens: Vec<RawGen>,
}

/// First words that start a section of the shared grammar (they end a raw
/// section).
fn is_section_word(w: &str) -> bool {
    matches!(
        w,
        "min"
            | "max"
            | "minimize"
            | "maximize"
            | "minimise"
            | "maximise"
            | "minimum"
            | "maximum"
            | "subject"
            | "such"
            | "st"
            | "s.t."
            | "st."
            | "bound"
            | "bounds"
            | "gen"
            | "general"
            | "generals"
            | "integer"
            | "integers"
            | "bin"
            | "binary"
            | "binaries"
            | "free"
            | "fre"
            | "semi"
            | "semis"
            | "semicontinuous"
            | "semi-continuous"
            | "semiinteger"
            | "end"
            | "sos"
    )
}

/// Splits the `SOS` and `General Constraints` sections out of `text`. The
/// returned text has those lines blanked (line numbers are kept).
///
/// Grammars: Gurobi "LP format" documentation: `SOS` header, then
/// `name: S1 :: x1 : 1 x2 : 2` (S1 or S2, pair of colons, `var : weight`
/// members); `general constraints` / `general constraint` / `gencons` /
/// `g.c.` header, then `[label:] r = MAX|MIN|AND|OR|ABS ( a , b , 0.7 )` and
/// `[label:] y = PWL ( x ) : (x0, y0) (x1, y1) ...`. CPLEX's LP format uses
/// the same `SOS` grammar (`name: S1 :: x1:1 x2:2`). `NOT ( a )` is a rusolve
/// extension (Gurobi has no NOT general constraint).
fn extract_raw_sections(text: &str) -> Result<(std::borrow::Cow<'_, str>, RawSections), LpError> {
    // Allocation-free prefilter: a model without these headers (every model
    // without general constraints) is tokenized from the original text, untouched.
    let has_header = text.lines().any(|l| {
        let code = l.split('\\').next().unwrap_or("").trim();
        code.len() <= 20
            && [
                "sos",
                "general constraints",
                "general constraint",
                "gencons",
                "g.c.",
            ]
            .iter()
            .any(|h| code.eq_ignore_ascii_case(h))
    });
    if !has_header {
        return Ok((std::borrow::Cow::Borrowed(text), RawSections::default()));
    }
    #[derive(PartialEq, Clone, Copy)]
    enum Mode {
        None,
        Sos,
        Gen,
    }
    let mut mode = Mode::None;
    let mut main = String::with_capacity(text.len());
    let mut raw = RawSections::default();
    for (i, line_text) in text.lines().enumerate() {
        let line = i + 1;
        let code = line_text.split('\\').next().unwrap_or("").trim();
        let lower = code.to_ascii_lowercase();
        let first = lower.split_whitespace().next().unwrap_or("");
        let gen_header = matches!(
            lower.as_str(),
            "general constraints" | "general constraint" | "gencons" | "g.c."
        );
        let is_header = gen_header || lower == "sos";
        if is_header || (mode != Mode::None && !code.is_empty() && is_section_word(first)) {
            mode = if gen_header {
                Mode::Gen
            } else if lower == "sos" {
                Mode::Sos
            } else {
                Mode::None
            };
            if is_header {
                main.push('\n');
                continue;
            }
        }
        if mode == Mode::None || code.is_empty() {
            main.push_str(line_text);
            main.push('\n');
            continue;
        }
        main.push('\n');
        match mode {
            Mode::Sos => parse_sos_line(code, line, &mut raw)?,
            Mode::Gen => raw.gens.push(parse_gen_line(code, line)?),
            Mode::None => unreachable!(),
        }
    }
    Ok((std::borrow::Cow::Owned(main), raw))
}

/// Spaces out the punctuation of the raw sections.
fn spaced(code: &str) -> Vec<String> {
    let mut s = String::with_capacity(code.len() * 2);
    for c in code.chars() {
        if matches!(c, ':' | '(' | ')' | ',' | '=') {
            s.push(' ');
            s.push(c);
            s.push(' ');
        } else {
            s.push(c);
        }
    }
    s.split_whitespace().map(str::to_string).collect()
}

fn parse_sos_line(code: &str, line: usize, raw: &mut RawSections) -> Result<(), LpError> {
    let t = spaced(code);
    let is_type = |s: &str| s.eq_ignore_ascii_case("s1") || s.eq_ignore_ascii_case("s2");
    let mut k = 0;
    let header = if t.len() >= 3 && t[1] == ":" && is_type(&t[2]) {
        k = 3;
        Some(t[0].clone())
    } else if !t.is_empty() && is_type(&t[0]) {
        k = 1;
        Some(format!("sos{}", raw.sos.len() + 1))
    } else {
        None
    };
    let set = match header {
        Some(name) => {
            let kind = if t[k - 1].eq_ignore_ascii_case("s1") {
                SosKind::Sos1
            } else {
                SosKind::Sos2
            };
            let mut priority = 0.0;
            if let Some(v) = t.get(k).and_then(|s| s.parse::<f64>().ok()) {
                priority = v;
                k += 1;
            }
            if t.get(k).map(String::as_str) != Some(":")
                || t.get(k + 1).map(String::as_str) != Some(":")
            {
                return err(line, "SOS constraint needs `::` after its type (S1 or S2)");
            }
            k += 2;
            raw.sos.push(RawSos {
                name,
                kind,
                priority,
                members: Vec::new(),
            });
            raw.sos.last_mut().expect("just pushed")
        }
        // a wrapped member list continues the previous set
        None => match raw.sos.last_mut() {
            Some(s) => s,
            None => {
                return err(
                    line,
                    "SOS section line must start with `name: S1 ::` or `name: S2 ::`",
                );
            }
        },
    };
    while k < t.len() {
        let name = t[k].clone();
        k += 1;
        let mut w = None;
        if t.get(k).map(String::as_str) == Some(":") {
            let Some(v) = t.get(k + 1).and_then(|s| s.parse::<f64>().ok()) else {
                return err(
                    line,
                    format!("SOS member '{name}' has no numeric weight after ':'"),
                );
            };
            w = Some(v);
            k += 2;
        }
        set.members.push((name, w));
    }
    Ok(())
}

fn parse_gen_line(code: &str, line: usize) -> Result<RawGen, LpError> {
    let t = spaced(code);
    let mut k = 0;
    let mut name = String::new();
    if t.len() >= 2 && t[1] == ":" {
        name = t[0].clone();
        k = 2;
    }
    let (Some(res), Some("=")) = (t.get(k), t.get(k + 1).map(String::as_str)) else {
        return err(
            line,
            "general constraint must be `[label:] result = KIND ( ... )`",
        );
    };
    let res = res.clone();
    let Some(kind) = t.get(k + 2).map(|s| s.to_ascii_uppercase()) else {
        return err(
            line,
            "general constraint has no type (MAX, MIN, AND, OR, ABS, PWL)",
        );
    };
    k += 3;
    if t.get(k).map(String::as_str) != Some("(") {
        return err(line, format!("expected '(' after {kind}"));
    }
    k += 1;
    let what = if kind == "PWL" {
        let Some(x) = t.get(k) else {
            return err(line, "PWL needs an operand variable");
        };
        let x = x.clone();
        if t.get(k + 1).map(String::as_str) != Some(")")
            || t.get(k + 2).map(String::as_str) != Some(":")
        {
            return err(
                line,
                "PWL syntax is `y = PWL ( x ) : (x0, y0) (x1, y1) ...`",
            );
        }
        k += 3;
        let mut pts = Vec::new();
        while k < t.len() {
            // ( x , y )
            let ok = t[k] == "("
                && t.get(k + 2).map(String::as_str) == Some(",")
                && t.get(k + 4).map(String::as_str) == Some(")");
            let (Some(a), Some(b)) = (
                t.get(k + 1).and_then(|s| s.parse::<f64>().ok()),
                t.get(k + 3).and_then(|s| s.parse::<f64>().ok()),
            ) else {
                return err(line, "PWL breakpoints are `(x, y)` pairs of numbers");
            };
            if !ok {
                return err(line, "PWL breakpoints are `(x, y)` pairs of numbers");
            }
            pts.push((a, b));
            k += 5;
        }
        RawGenKind::Pwl(x, pts)
    } else {
        let fk = match kind.as_str() {
            "MAX" => FuncKind::Max,
            "MIN" => FuncKind::Min,
            "ABS" => FuncKind::Abs,
            "AND" => FuncKind::And,
            "OR" => FuncKind::Or,
            "NOT" => FuncKind::Not,
            other => {
                return err(
                    line,
                    format!(
                        "general constraint type {other} is not supported (MAX, MIN, AND, OR, ABS, PWL)"
                    ),
                );
            }
        };
        let mut args = Vec::new();
        let mut closed = false;
        while k < t.len() {
            match t[k].as_str() {
                ")" => {
                    closed = true;
                    k += 1;
                    break;
                }
                "," => k += 1,
                s => {
                    args.push(match s.parse::<f64>() {
                        Ok(v) => RawArg::Const(v),
                        Err(_) => RawArg::Var(s.to_string()),
                    });
                    k += 1;
                }
            }
        }
        if !closed || k != t.len() {
            return err(
                line,
                format!("{kind} constraint is not closed by ')' at the end of the line"),
            );
        }
        RawGenKind::Func(fk, args)
    };
    Ok(RawGen {
        line,
        name,
        res,
        what,
    })
}

/// Parse a CPLEX LP-format document into the shared [`Model`]. A model with
/// quadratic constraint terms is an error here (see [`read_lp_conic`]).
pub fn read_lp(text: &str) -> Result<Model, LpError> {
    let m = read_lp_conic(text)?;
    if !m.conic.is_empty() {
        crate::conic::note_conic_rejection();
        return Err(LpError {
            line: 0,
            message: crate::conic::CONIC_CONTENT_MSG.to_string(),
        });
    }
    Ok(m)
}

/// [`read_lp`] keeping quadratic constraint terms in [`Model::conic`].
pub fn read_lp_conic(text: &str) -> Result<Model, LpError> {
    let (main_text, raw) = extract_raw_sections(text)?;
    let mut p = P {
        toks: tokenize(&main_text)?,
        pos: 0,
        sense: ObjSense::default(),
        obj_name: None,
        obj_terms: Vec::new(),
        obj_const: 0.0,
        obj_quad: Vec::new(),
        obj_seen: false,
        rows: Vec::new(),
        bound_ops: Vec::new(),
        default: None,
        lists: Vec::new(),
        inds: Vec::new(),
        raw,
    };
    let mut m = p.parse_model()?;
    m.name = model_name(text);
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A statement wrapped over many lines (HiGHS and CPLEX wrap a long objective) is scanned
    /// once: before the fix each continuation line rescanned the statement from its start,
    /// so a 1M-term objective took 16 s instead of 2 (rail4284.lp).
    const QLP: &str = "\
Minimize
 obj: x + y
Subject To
 q1: x + [ 2 x ^ 2 + 3 x * y + y ^ 2 ] <= 5
 l1: x + y >= 1
Bounds
 x <= 10
 y <= 10
End
";

    #[test]
    fn constraint_quadratic_terms_fill_the_conic_block() {
        let m = read_lp_conic(QLP).unwrap();
        assert_eq!(m.conic.quad_rows.len(), 1);
        let q = &m.conic.quad_rows[0];
        assert_eq!(q.row, 0);
        assert_eq!(
            q.monomials().collect::<Vec<_>>(),
            [(0, 0, 2.0), (1, 0, 3.0), (1, 1, 1.0)]
        );
        // the linear part and the rhs are untouched
        assert_eq!(m.rhs, [5.0, 1.0]);
        assert_eq!(m.matrix.values().len(), 3);
        let e = read_lp(QLP).unwrap_err();
        assert_eq!(e.message, crate::conic::CONIC_CONTENT_MSG);
    }

    #[test]
    fn constraint_quadratic_terms_round_trip_and_match_mps() {
        let m = read_lp_conic(QLP).unwrap();
        let again = read_lp_conic(&crate::lp_write::write_lp(&m)).unwrap();
        assert_eq!(again.conic, m.conic);
        let via_mps = crate::mps::read_mps_conic(&crate::mps::write_mps(&m)).unwrap();
        assert_eq!(via_mps.conic, m.conic);
    }

    #[test]
    fn quadratic_terms_in_two_sided_rows_are_errors() {
        let t = "Minimize\n obj: x\nSubject To\n q: 1 <= x + [ x ^ 2 ] <= 4\nEnd\n";
        let e = read_lp_conic(t).unwrap_err();
        assert!(e.message.contains("two-sided"), "{e}");
    }

    #[test]
    fn long_wrapped_statement_is_linear() {
        let n = 200_000;
        let mut text = String::from("Minimize\n obj:");
        for j in 0..n {
            text.push_str(&format!(" + {} x{j}\n", j % 7 + 1));
        }
        text.push_str("Subject To\n c: x0 + x1 >= 1\nEnd\n");
        let t = std::time::Instant::now();
        let m = read_lp(&text).unwrap();
        assert_eq!(m.num_variables, n);
        assert_eq!(m.objective[n - 1], ((n - 1) % 7 + 1) as f64);
        // quadratic rescanning takes minutes here; linear takes well under a second
        assert!(t.elapsed().as_secs_f64() < 5.0, "took {:?}", t.elapsed());
    }

    #[test]
    fn parses_basic_minimize_with_comments_and_semicolons() {
        let text = "\
\\ LP format for model diet
\\ a header comment
Minimize
 obj: 2 x + 3 y ;
Subject To
 c1: x + y >= 1
 c2: 2 x <= 4
Bounds
 y <= 3
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.name, "diet", "model name from the marker comment");
        assert_eq!(m.obj_sense, ObjSense::Minimize);
        assert_eq!(m.obj_name.as_deref(), Some("obj"));
        assert_eq!(m.col_names, ["x", "y"]);
        assert_eq!(m.row_names, ["c1", "c2"]);
        assert_eq!(m.objective, [2.0, 3.0]);
        assert_eq!(m.senses, [RowSense::GreaterEqual, RowSense::LessEqual]);
        assert_eq!(m.rhs, [1.0, 4.0]);
        assert_eq!(m.ranges, [0.0, 0.0]);
        assert_eq!(m.matrix.get(1, 0), 2.0);
        // LP default bounds 0 ≤ x < ∞; explicit UP on y only.
        assert_eq!((m.bounds[0].lower, m.bounds[0].upper), (0.0, f64::INFINITY));
        assert_eq!(m.bounds[1].upper, 3.0);
    }

    #[test]
    fn parses_implicit_coefficients_and_scientific_notation() {
        let text = "\
Maximize
 obj: x - y + 3z + 2 w + 1.5e3 a + 2.5E-4 b + .5 c
Subject To
 c1: x + y + z + w + a + b + c <= 10
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.obj_sense, ObjSense::Maximize);
        assert_eq!(m.objective, [1.0, -1.0, 3.0, 2.0, 1500.0, 2.5e-4, 0.5]);
        assert_eq!(m.col_names, ["x", "y", "z", "w", "a", "b", "c"]);
    }

    #[test]
    fn parses_two_sided_and_left_led_rows_as_sense_rhs_range() {
        let text = "\
Min
 obj: x + y
s.t.
 r1: 1 <= x + y <= 5
 r2: 5 >= x >= 1
 r3: 1 <= x + y
 r4: 9 >= x
 r5: -2 <= x <= 2
End
";
        let m = read_lp(text).unwrap();
        // Two-sided rows land on the same base-sense + RANGE representation
        // the MPS RANGES reader produces: L hi | hi−lo.
        assert_eq!(m.senses[0], RowSense::LessEqual);
        assert_eq!(m.rhs[0], 5.0);
        assert_eq!(m.ranges[0], 4.0);
        assert_eq!(m.rhs[1], 5.0);
        assert_eq!(m.ranges[1], 4.0);
        // Left-led one-sided rows keep a single side, no range.
        assert_eq!(m.senses[2], RowSense::GreaterEqual);
        assert_eq!(m.rhs[2], 1.0);
        assert_eq!(m.senses[3], RowSense::LessEqual);
        assert_eq!(m.rhs[3], 9.0);
        assert_eq!(m.ranges[2..], [0.0, 0.0, 4.0]);
        assert_eq!(m.rhs[4], 2.0);
        // A left-hand constant is ignored, as HiGHS does: x + 2 <= 10 ⇒ rhs 10.
        let m2 = read_lp("Min\n o: x\nSt\n c: x + 2 <= 10\nEnd\n").unwrap();
        assert_eq!(m2.rhs, [10.0]);
    }

    #[test]
    fn parses_equality_and_bare_expression_rows() {
        let text = "\
Minimize
 obj: x
Subject To
 e1: x + y = 4
 bare: x - y
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.senses, [RowSense::Equal, RowSense::Equal]);
        assert_eq!(m.rhs, [4.0, 0.0], "bare expression row is `= 0`");
    }

    #[test]
    fn parses_all_bound_forms() {
        let text = "\
Minimize
 obj: x + y + z + w + v + u + a + b
Subject To
 c: x + y + z + w + v + u + a + b <= 100
Bounds
 x free
 1 <= y <= 4
 z in [2, 3]
 w <= 10
 v >= -2
 u = 7
 a, b <= 5
End
";
        let m = read_lp(text).unwrap();
        let lu = |n: &str| {
            let j = m.col_names.iter().position(|c| c == n).unwrap();
            (m.bounds[j].lower, m.bounds[j].upper)
        };
        assert_eq!(lu("x"), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(lu("y"), (1.0, 4.0));
        assert_eq!(lu("z"), (2.0, 3.0));
        assert_eq!(lu("w"), (0.0, 10.0));
        assert_eq!(lu("v"), (-2.0, f64::INFINITY));
        assert_eq!(lu("u"), (7.0, 7.0));
        assert_eq!(lu("a"), (0.0, 5.0));
        assert_eq!(lu("b"), (0.0, 5.0));
    }

    #[test]
    fn bounds_default_applies_only_to_variables_without_explicit_bounds() {
        let text = "\
Minimize
 obj: x + y
Subject To
 c: x + y <= 10
Bounds
 default free
 y <= 3
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(
            (m.bounds[0].lower, m.bounds[0].upper),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
        assert_eq!((m.bounds[1].lower, m.bounds[1].upper), (0.0, 3.0));

        let m2 = read_lp("Min\n o: x\nSt\n c: x <= 1\nBounds\n default in [1, 2]\nEnd\n").unwrap();
        assert_eq!((m2.bounds[0].lower, m2.bounds[0].upper), (1.0, 2.0));
    }

    #[test]
    fn parses_binary_general_free_sections() {
        let text = "\
Minimize
 obj: x + y + z + w
Subject To
 c: x + y + z + w <= 4
Binaries
 x, y
General z
Free
 w
End
";
        let m = read_lp(text).unwrap();
        let b = |n: &str| m.bounds[m.col_names.iter().position(|c| c == n).unwrap()];
        assert_eq!(
            b("x"),
            VariableBound {
                lower: 0.0,
                upper: 1.0,
                integer: true,
            }
        );
        assert!(b("y").integer && b("y").upper == 1.0);
        // General: integer, bounds untouched (LP default 0..inf).
        assert!(b("z").integer);
        assert_eq!((b("z").lower, b("z").upper), (0.0, f64::INFINITY));
        // Free section: continuous free variable.
        assert!(!b("w").integer);
        assert_eq!(
            (b("w").lower, b("w").upper),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
    }

    #[test]
    fn keywords_are_case_insensitive() {
        let text = "\
MAXIMUM
 obj: x
St
 c: x >= 1
bOuNdS
 x <= 5
EnD
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.obj_sense, ObjSense::Maximize);
        assert_eq!(m.rhs, [1.0]);
        assert_eq!(m.bounds[0].upper, 5.0);
    }

    #[test]
    fn space_separated_names_are_rejected() {
        // Names are single tokens (HiGHS grammar): no space-joined names.
        let e = read_lp("Min\n obj: my var\nEnd\n").unwrap_err();
        assert!(e.message.contains("expected `+` or `-`"), "{e}");
        let e = read_lp("Min\n total cost: x\nEnd\n").unwrap_err();
        assert!(e.message.contains("expected `+` or `-`"), "{e}");
        let e = read_lp("Min\n o: x\nSt\n c: my var + x <= 10\nEnd\n").unwrap_err();
        assert!(e.message.contains("expected `+` or `-`"), "{e}");
    }

    #[test]
    fn parses_punctuated_variable_names() {
        let text = "\
Minimize
 obj: x_1 + a.b + x' + p/q + x{1}
Subject To
 c: x_1 + a.b + x' + p/q + x{1} <= 5
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.col_names, ["x_1", "a.b", "x'", "p/q", "x{1}"]);
        assert!(m.objective.iter().all(|&c| c == 1.0));
    }

    #[test]
    fn objective_constant_becomes_obj_offset() {
        let text = "\
Maximize
 obj: 2 x + 5
Subject To
 c: x <= 3
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.obj_offset, 5.0);
        assert_eq!(m.objective, [2.0]);
    }

    #[test]
    fn parses_wrapped_statements() {
        let text = "\
Minimize
 obj: 2 x +
   3 y
Subject To
 c1: x +
   2 y <=
   10
 c2: 1 <=
   x + y <=
   4
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.objective, [2.0, 3.0]);
        assert_eq!(m.senses, [RowSense::LessEqual, RowSense::LessEqual]);
        assert_eq!(m.rhs, [10.0, 4.0]);
        assert_eq!(m.ranges, [0.0, 3.0]);
        assert_eq!(m.matrix.get(1, 1), 1.0);
    }

    #[test]
    fn variables_first_seen_in_bounds_or_lists_are_created() {
        let text = "\
Minimize
 obj: x
Subject To
 c: x <= 5
Bounds
 ghost free
Binaries
 hidden
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.num_variables, 3);
        assert_eq!(m.col_names, ["x", "ghost", "hidden"]);
        assert_eq!(m.bounds[1].lower, f64::NEG_INFINITY);
        assert!(m.bounds[2].integer);
    }

    #[test]
    fn unnamed_rows_get_empty_names() {
        let text = "\
Minimize
 2 x + 3 y
Subject To
 x + y >= 2
 x <= 5
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.obj_name, None);
        assert_eq!(m.row_names, ["", ""]);
        assert_eq!(m.num_constraints, 2);
        assert_eq!(m.senses, [RowSense::GreaterEqual, RowSense::LessEqual]);
    }

    #[test]
    fn two_sided_rows_survive_to_sparse_lp() {
        let m = read_lp("Min\n o: x\nSt\n r: 1 <= x <= 5\nEnd\n").unwrap();
        let sp = m.to_sparse_lp().unwrap();
        // L 5 | 4 ⇒ base row x ≤ 5 plus appended complement x ≥ 1.
        let dense = sp.to_lp();
        assert_eq!(dense.senses, ['L', 'G']);
        assert_eq!(dense.rhs, [5.0, 1.0]);
    }

    #[test]
    fn errors_are_reported_with_line_numbers() {
        // `<>` is two comparisons, not a comparison followed by a number
        let e = read_lp("Min\n o: x\nSt\n c: x <> 5\nEnd\n").unwrap_err();
        assert!(e.message.contains("expected a number"), "{e}");
        assert_eq!(e.line, 4);
        // parentheses are name characters (CPLEX and HiGHS), not grouping
        let m = read_lp("Min\n o: x\nSt\n c: (x + y) <= 5\nEnd\n").unwrap();
        assert_eq!(m.col_names, ["x", "(x", "y)"]);
        // stray token outside a section
        let e = read_lp("foo\nMin\n o: x\nEnd\n").unwrap_err();
        assert!(e.message.contains("section keyword"), "{e}");
        // mixed operators in a two-sided chain
        let e = read_lp("Min\n o: x\nSt\n c: 1 <= x >= 5\nEnd\n").unwrap_err();
        assert!(e.message.contains("mixed relational"), "{e}");
        // inverted range
        let e = read_lp("Min\n o: x\nSt\n c: 5 <= x <= 1\nEnd\n").unwrap_err();
        assert!(e.message.contains("inverted"), "{e}");
        // second objective
        let e = read_lp("Min\n o: x\nMax\n o2: x\nEnd\n").unwrap_err();
        assert!(e.message.contains("second objective"), "{e}");
        // semi-continuous is a model object
        let m = read_lp("Min\n o: x\nSt\n c: x <= 1\nSemi-continuous\n x\nEnd\n").unwrap();
        assert_eq!(m.gc.semi.len(), 1);
    }

    #[test]
    fn h16_lpio_accepts_bare_and_equals_relations() {
        // `<` / `>` and `=<` / `=>` are the HiGHS spellings of `<=` / `>=`.
        let text = "\
Min
 obj: x + y
St
 c1: x < 5
 c2: x > 1
 c3: x =< 6
 c4: y => 2
Bounds
 x > 1
 y =< 3
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(
            m.senses,
            [
                RowSense::LessEqual,
                RowSense::GreaterEqual,
                RowSense::LessEqual,
                RowSense::GreaterEqual
            ]
        );
        assert_eq!(m.rhs, [5.0, 1.0, 6.0, 2.0]);
        assert_eq!((m.bounds[0].lower, m.bounds[0].upper), (1.0, f64::INFINITY));
        assert_eq!((m.bounds[1].lower, m.bounds[1].upper), (0.0, 3.0));
    }

    #[test]
    fn h16_lpio_accepts_infinity_bounds() {
        let text = "\
Min
 obj: x + y + z
St
 c: x + y + z <= 10
Bounds
 -inf <= x <= 5
 y >= -infinity
 z <= +Inf
End
";
        let m = read_lp(text).unwrap();
        let lu = |n: &str| {
            let j = m.col_names.iter().position(|c| c == n).unwrap();
            (m.bounds[j].lower, m.bounds[j].upper)
        };
        assert_eq!(lu("x"), (f64::NEG_INFINITY, 5.0));
        assert_eq!(lu("y"), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(lu("z"), (0.0, f64::INFINITY));
    }

    #[test]
    fn h16_lpio_accepts_empty_semi_sections() {
        // HiGHS writes an empty semi section for a model without semi variables.
        for text in [
            "Min\n obj: x\nSt\n c: x >= 1\nSemi-continuous\nEnd\n",
            "Min\n obj: x\nSt\n c: x >= 1\nsemi\nEnd\n",
            "Min\n obj: x\nSt\n c: x >= 1\nSemi-continuous\nBinaries\n x\nEnd\n",
        ] {
            let m = read_lp(text).unwrap_or_else(|e| panic!("{text:?}: {e}"));
            assert_eq!(m.num_variables, 1, "{text:?}");
        }
    }

    #[test]
    fn h16_lpio_rejects_adjacent_terms_without_sign() {
        for text in [
            "Min\n obj: x\nSt\n c: 3 x 4 y <= 5\nEnd\n",
            "Min\n obj: x\nSt\n c: x 3 <= 5\nEnd\n",
            "Min\n obj: 2 3 x\nEnd\n",
            "Min\n obj: x\nSt\n c: x y <= 1\nEnd\n",
        ] {
            let e = read_lp(text).unwrap_err();
            assert!(e.message.contains("expected `+` or `-`"), "{text:?}: {e}");
        }
    }

    #[test]
    fn h16_lpio_rejects_relation_in_objective_and_keyword_mid_line() {
        let e = read_lp("Min\n obj: x = 10\nSt\n c: x >= 1\nEnd\n").unwrap_err();
        assert!(e.message.contains("relation in objective"), "{e}");
        assert_eq!(e.line, 2);
        assert!(read_lp("Minimize a subject to a >= 1\nEnd\n").is_err());
    }

    #[test]
    fn h16_lpio_ignores_lhs_constant_but_keeps_objective_constant() {
        let m = read_lp("Min\n o: x\nSt\n c: x + 2 >= 5\nEnd\n").unwrap();
        assert_eq!(m.rhs, [5.0]);
        let m = read_lp("Max\n obj: 2 x + 5\nSt\n c: x <= 3\nEnd\n").unwrap();
        assert_eq!(m.obj_offset, 5.0);
        assert_eq!(m.rhs, [3.0]);
    }

    #[test]
    fn h16_lpio_reads_highs_written_mip() {
        // Byte-for-byte the 1.16.0 FilereaderLp::writeModelToFile layout
        // (%+.17g values, `lo`/`up` boxed rows, `-inf <= x <= 5`, empty
        // `bin` and `semi` sections) for min x + 2 y + 3, x integer.
        let text = "\\ File written by HiGHS .lp file handler
min
 obj: +1 x +2 y  +3
st
 c1: +1 x +1 y +1 w >= +1
 c2lo: +1 x -1 y >= +0
 c2up: +1 x -1 y <= +4
bounds
 -inf <= x <= 5
 w free
bin
gen
 x
semi
end
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.col_names, ["x", "y", "w"]);
        assert_eq!(m.objective, [1.0, 2.0, 0.0]);
        assert_eq!(m.obj_offset, 3.0);
        assert_eq!(m.row_names, ["c1", "c2lo", "c2up"]);
        assert_eq!(
            m.senses,
            [
                RowSense::GreaterEqual,
                RowSense::GreaterEqual,
                RowSense::LessEqual
            ]
        );
        assert_eq!(m.rhs, [1.0, 0.0, 4.0]);
        assert_eq!(
            (m.bounds[0].lower, m.bounds[0].upper),
            (f64::NEG_INFINITY, 5.0)
        );
        assert!(m.bounds[0].integer);
        assert_eq!((m.bounds[1].lower, m.bounds[1].upper), (0.0, f64::INFINITY));
        assert!(!m.bounds[1].integer);
        assert_eq!(
            (m.bounds[2].lower, m.bounds[2].upper),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
    }

    #[test]
    fn h16_lpio_reads_issue_1451_constant_lhs() {
        // check/instances/1451.lp as changed by #3283: the LHS constant
        // `- 1` is ignored in every row, so x >= 9 and HiGHS's optimum is 9
        // (TestFilereader.cpp); a row named `end` is a label, not a keyword.
        let text = "\
Minimize
  x
Subject to
\\ Name this constraint end to test for keywords as constraint names
 end: x - 1 >= -1
 cs0: x - 1 >= 0
 cs1: x - 1 >= 1
 cs2: x - 1 >= 2
 cs3: x - 1 >= 3
 cs4: x - 1 >= 4
 cs5: x - 1 >= 5
 cs6: x - 1 >= 6
 cs7: x - 1 >= 7
 cs8: x - 1 >= 8
 cs9: x - 1 >= 9
Integer
  x
End
";
        let m = read_lp(text).unwrap();
        assert_eq!(m.row_names.len(), 11);
        assert_eq!(m.row_names[0], "end");
        assert_eq!(
            m.rhs,
            [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        );
        assert!(m.senses.iter().all(|s| *s == RowSense::GreaterEqual));
        assert!(m.bounds[0].integer);
    }

    #[test]
    fn h16_lpio_corrupt_and_good_files() {
        // #3343 TestFilereader "lp-file-format-corrupt-files" (linear cases)
        // and #3355 TestLPFileFormat error cases the linear reader can see.
        for text in [
            "min\n obj: x y\nst\n c: x + y >= 1\nend\n",
            "min\n obj: x + y\nst\n c: x + y >= 1 + 1\nend\n",
            "min\n obj: x + y\nst\n c: x y >= 1\nend\n",
            "min\n obj: 2 3 x\nst\n c: x >= 1\nend\n",
            "min\n obj: x 3\nst\n c: x >= 1\nend\n",
            "x + y\nmin\n obj: x + y\nst\n c: x + y >= 1\nend\n",
            "comment\nmin\nend",
            "subject\nto\n c: x >= 1\nend\n",
            "min\n x\nend\nx\n",
            "min\n obj: x semi\nend\n",
            "min\n obj: x such that y\nend\n",
            "min\n obj: x + 0.2 5C0ST\nend\n",
            "min\n obj: x +\nst\nend\n",
            "min\n obj: x\n y\nend\n",
            "min\n x\nst\n c: x + >= 1\nend\n",
            "min\n x\nst\n c: 3 x >= 2 y\nend\n",
            "min\n x\nsemi\n- continuous\nend\n",
            "min\n x\nsemi -\ncontinuous\nend\n",
            "min\n x\nbounds\n x 1\nend\n",
            "min\n x\nbounds\n 1 = x <= 2\nend\n",
            "min\n x\nbounds\n 1 <= x >= 2\nend\n",
            "min\n x\nbounds\n x <= 1 y <= 2\nend\n",
        ] {
            assert!(read_lp(text).is_err(), "{text:?} should be rejected");
        }
        for text in [
            "min\n obj: x + y\nst\n c: x + y >= 1\nend\n",
            "min\n obj: - x - 2 y + 3\nst\n c: x\n + y >= 1\nend\n",
            "min\n obj: 2 + 3 x\nst\n c: x >= 1\nend\n",
            "\\ comment\nmin\n obj: x + y\nst\n c: x + y >= 1\nend\n",
            "\\ comment",
            "minimise\n obj: x\nsuch that\n c: x >= 1\nintegers\n x\nend\n",
        ] {
            read_lp(text).unwrap_or_else(|e| panic!("{text:?}: {e}"));
        }
        // `5 = x` fixes x (parseBound `<number> <comparison> <identifier>`).
        let m = read_lp("min\n obj: x\nbounds\n 5 = x\nend\n").unwrap();
        assert_eq!((m.bounds[0].lower, m.bounds[0].upper), (5.0, 5.0));
    }
}
