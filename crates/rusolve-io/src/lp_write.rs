//! CPLEX LP format writer: a [`Model`] as an `.lp` text that IBM ILOG CPLEX
//! reads back to the same columns, rows, names, values and types.
//!
//! Provenance: the grammar is the public IBM ILOG CPLEX "LP file format"
//! documentation (sections, `name:` labels, `[ ... ] / 2` quadratic objective,
//! `Bounds` / `Binaries` / `Generals`, `inf`, `free`). Beyond the grammar,
//! the writer follows these reading rules of the format so that every
//! value, name and index survives the round trip:
//!
//! - **Column order** is first appearance, objective first; a `+ 0 x` term
//!   in the objective creates `x`. So every column is written in the
//!   objective, in model order, zero cost included, and every column keeps
//!   its index.
//! - **Ranged rows** have no LP syntax (`lo <= expr <= hi` is rejected,
//!   `expr >= lo <= hi` silently drops the second side). A ranged row is
//!   written `name: expr - Rg<name> = rhs` with the range column `Rg<name>`
//!   bounded so that the row interval is exactly the MPS one (no arithmetic):
//!   `L`: `[-|r|, 0]`, `G`: `[0, |r|]`, `E`: `[0, r]` or `[r, 0]`. The range
//!   columns come after every model column. CPLEX keeps them as columns.
//! - **Numbers**: CPLEX's decimal reader is not correctly rounded (a
//!   17-digit decimal can land a few ulps away). It accumulates integer and
//!   fraction digits in doubles, adds `int + frac / 10^m`, and scales by a
//!   power of ten built by repeated multiplication; an exponent below -307
//!   reads as zero. [`cplex_read`] models that reader. Each value is written as the
//!   shortest decimal that both CPLEX ([`cplex_read`]) and a correctly
//!   rounded reader (HiGHS, rusolve) map to the same bits. Where no
//!   such string exists (only for magnitudes beyond about 1e±90), the
//!   string CPLEX reads exactly wins and [`LpWriteReport::cplex_only_numbers`]
//!   counts it; values CPLEX cannot read at all (subnormals, exponents below
//!   -307) are written shortest-exact and counted in
//!   [`LpWriteReport::unreadable_numbers`]. `-0` is written as `-0` (CPLEX
//!   keeps the sign in bounds and right-hand sides).
//! - **Finite values >= 1e20** are written as they are. CPLEX treats bounds
//!   at or above 1e20 as infinite, both when it reads this file and when it
//!   reads the MPS original, so the two readings agree.
//! - **Names**: CPLEX accepts letters, digits, `. ! " # $ % & ' ( ) , ; ? @ _
//!   \` { } ~` and non-ASCII bytes; it rejects `* + - / : < = > [ ] ^ | \`; a
//!   leading digit or `.` is an error; a column name that starts with `e`/`E`
//!   followed by nothing, a digit, `.` or `e` is read as an exponent; section
//!   keywords (`end`, `gen`, `bin`, ...) at the start of a `Bounds` or list
//!   line start a section; column names are cut at 254 bytes, row names at
//!   255. A name breaking any of those, or one that the HiGHS and rusolve
//!   readers split (`;`, a comma at either end or doubled, whitespace), is
//!   renamed deterministically: bad characters become `_`, a `_` is put in
//!   front when the start is the problem, empty names become `x<j+1>` /
//!   `c<i+1>` (CPLEX's own default names), and a clash with a name already
//!   taken gets `~1`, `~2`, ... Legal, unique names are never touched.
//! - **Line length** has no documented limit; lines
//!   are wrapped near 255 characters anyway, each continuation line starting
//!   with a sign so a keyword-like name can never open it.

use std::collections::{BTreeMap, HashSet};
use std::fmt::Write as _;
use std::sync::OnceLock;

use crate::gc::SosKind;
use crate::model::{Model, ObjSense, RowSense};

/// Wrap lines once they would pass this many bytes.
const LINE_WRAP: usize = 255;
/// The CPLEX LP reader keeps at most this many bytes of a column name ...
const MAX_COL_NAME: usize = 254;
/// ... and of a row or objective name.
const MAX_ROW_NAME: usize = 255;

/// What [`write_lp_report`] had to change or could not keep exact.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LpWriteReport {
    /// Columns written under a new name (illegal, empty, too long or a
    /// duplicate).
    pub renamed_columns: usize,
    /// Rows written under a new name.
    pub renamed_rows: usize,
    /// The objective label was renamed.
    pub renamed_objective: bool,
    /// Ranged rows, each written with an extra `Rg<row>` column.
    pub ranged_rows: usize,
    /// Values written in the form CPLEX reads exactly although a correctly
    /// rounded reader lands on a neighbouring double.
    pub cplex_only_numbers: usize,
    /// Values CPLEX reads as a different number whatever the string
    /// (subnormals, magnitudes below about 1e-307).
    pub unreadable_numbers: usize,
}

/// The names [`write_lp`] uses: one per model column and row, the objective
/// label, and the range column of each ranged row (`None` when not ranged).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LpNames {
    pub columns: Vec<String>,
    pub rows: Vec<String>,
    pub objective: Option<String>,
    pub range_columns: Vec<Option<String>>,
}

/// The model as CPLEX LP text. Quadratic rows are written; cones have no LP
/// form and are not ([`write_model_file`] refuses a model that has any).
pub fn write_lp(m: &Model) -> String {
    write_lp_report(m).0
}

/// [`write_lp`] plus a report of renames and inexact values.
pub fn write_lp_report(m: &Model) -> (String, LpWriteReport) {
    let mut report = LpWriteReport::default();
    let names = lp_names_counted(m, &mut report);
    let mut w = Writer {
        out: String::with_capacity(estimate_size(m)),
        line_len: 0,
        num: NumFmt::default(),
    };
    w.model(m, &names);
    report.cplex_only_numbers = w.num.cplex_only;
    report.unreadable_numbers = w.num.unreadable;
    (w.out, report)
}

/// The names [`write_lp`] gives every column, row, the objective and the
/// range columns.
pub fn lp_names(m: &Model) -> LpNames {
    lp_names_counted(m, &mut LpWriteReport::default())
}

/// Write `m` to `path`, LP or MPS by extension (`.lp`, `.mps`, either with
/// `.gz` for gzip).
pub fn write_model_file(path: &str, m: &Model) -> std::io::Result<()> {
    write_model_file_report(path, m).map(|_| ())
}

/// [`write_model_file`], returning the [`LpWriteReport`] when the format is
/// LP. An unknown extension is an `InvalidInput` error and writes nothing.
pub fn write_model_file_report(path: &str, m: &Model) -> std::io::Result<Option<LpWriteReport>> {
    let lower = path.to_ascii_lowercase();
    let gz = lower.ends_with(".gz");
    let stem = lower.strip_suffix(".gz").unwrap_or(&lower);
    let (text, report) = if stem.ends_with(".lp") {
        if let Some(c) = m.conic.cones.first() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "{path}: the LP format has no cones (cone '{}'); write .mps instead",
                    c.name
                ),
            ));
        }
        let (t, r) = write_lp_report(m);
        (t, Some(r))
    } else if stem.ends_with(".mps") {
        if !m.gc.pwl.is_empty() || !m.gc.funcs.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "{path}: MPS output cannot hold {} (write .lp instead)",
                    m.gc.pwl.first().map_or_else(
                        || {
                            let f = &m.gc.funcs[0];
                            format!("{} '{}'", f.kind.word(), f.name)
                        },
                        |p| format!("pwl '{}'", p.name)
                    )
                ),
            ));
        }
        (crate::mps::write_mps(m), None)
    } else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("{path}: unknown output format (expected .lp, .mps, .lp.gz or .mps.gz)"),
        ));
    };
    if gz {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::io::Write as _;
        let f = std::fs::File::create(path)?;
        let mut e = GzEncoder::new(std::io::BufWriter::new(f), Compression::default());
        e.write_all(text.as_bytes())?;
        e.finish()?.flush()?;
    } else {
        std::fs::write(path, text)?;
    }
    Ok(report)
}

fn estimate_size(m: &Model) -> usize {
    64 + 16 * (m.num_variables + m.num_constraints) + 24 * (m.matrix.values().len())
}

// ------------------------------------------------------------------ numbers

/// Powers of ten as the CPLEX LP reader builds them: repeated multiplication by 10
/// (exact up to 1e22, drifting after; `inf` past 1e308).
fn pow10_rep(n: usize) -> f64 {
    static TABLE: OnceLock<Vec<f64>> = OnceLock::new();
    let t = TABLE.get_or_init(|| {
        let mut t = Vec::with_capacity(400);
        let mut p = 1.0f64;
        for _ in 0..400 {
            t.push(p);
            p *= 10.0;
        }
        t
    });
    t.get(n).copied().unwrap_or(f64::INFINITY)
}

/// The value the CPLEX LP reader gives the decimal `s` (`[+-]digits[.digits][e[+-]digits]`):
/// integer and fraction digits accumulate in doubles
/// (`acc * 10 + d`), the value is `int + frac / 10^m`, then it is multiplied
/// or divided by `10^|e|` from [`pow10_rep`]; an exponent below -307 gives 0.
/// `None` when `s` is not of that shape.
pub fn cplex_read(s: &str) -> Option<f64> {
    let b = s.as_bytes();
    let mut i = 0;
    let neg = match b.first() {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let mut int = 0.0f64;
    let mut nd = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        int = int * 10.0 + f64::from(b[i] - b'0');
        i += 1;
        nd += 1;
    }
    let mut frac = 0.0f64;
    let mut nf = 0usize;
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            frac = frac * 10.0 + f64::from(b[i] - b'0');
            i += 1;
            nf += 1;
        }
    }
    if nd + nf == 0 {
        return None;
    }
    let mut v = if nf > 0 {
        int + frac / pow10_rep(nf)
    } else {
        int
    };
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        i += 1;
        let eneg = match b.get(i) {
            Some(b'-') => {
                i += 1;
                true
            }
            Some(b'+') => {
                i += 1;
                false
            }
            _ => false,
        };
        let start = i;
        let mut e: i64 = 0;
        while i < b.len() && b[i].is_ascii_digit() {
            e = (e * 10 + i64::from(b[i] - b'0')).min(100_000);
            i += 1;
        }
        if i == start {
            return None;
        }
        if eneg && e > 307 {
            v = 0.0;
        } else if eneg {
            v /= pow10_rep(e as usize);
        } else if e > 0 {
            v *= pow10_rep(e as usize);
        }
    }
    if i != b.len() {
        return None;
    }
    Some(if neg { -v } else { v })
}

#[derive(Default)]
struct NumFmt {
    cplex_only: usize,
    unreadable: usize,
}

/// Significant digits and decimal exponent of `v > 0` from Rust's `{:e}` /
/// `{:.*e}`: `v ~ d.ddd * 10^e10`.
fn sci_digits(v: f64, sig: Option<usize>) -> (String, i32) {
    let s = match sig {
        None => format!("{v:e}"),
        Some(n) => format!("{:.*e}", n - 1, v),
    };
    let (mant, exp) = s.split_once('e').expect("{:e} has an exponent");
    (mant.replace('.', ""), exp.parse().expect("{:e} exponent"))
}

/// Every placement of the decimal point in the digit string `d` (value
/// `d * 10^(e10 - len + 1)`), with a compensating exponent: `0.0ddd`,
/// `d.dd`, `dd.d`, `ddd0`, each `e<q>`; plus the placements that make the
/// exponent +-22, the largest power of ten CPLEX builds exactly (`1e25` is
/// read inexactly, `1000e22` exactly).
fn layouts(d: &str, e10: i32, out: &mut Vec<String>) {
    let n = d.len() as i32;
    let near = -2..=n + 2;
    let far = [e10 + 1 - 22, e10 + 1 + 22]
        .into_iter()
        .filter(|p| !near.contains(p) && (-40..=n + 40).contains(p));
    for p in near.clone().chain(far) {
        let q = e10 - p + 1;
        let mut s = String::with_capacity(d.len() + 8);
        if p <= 0 {
            s.push_str("0.");
            for _ in 0..-p {
                s.push('0');
            }
            s.push_str(d);
        } else if p >= n {
            s.push_str(d);
            for _ in 0..p - n {
                s.push('0');
            }
        } else {
            s.push_str(&d[..p as usize]);
            s.push('.');
            s.push_str(&d[p as usize..]);
        }
        if q != 0 {
            let _ = write!(s, "e{q}");
        }
        out.push(s);
    }
}

impl NumFmt {
    /// Append `v` (finite or infinite, any sign) to `out`.
    fn push(&mut self, out: &mut String, v: f64) {
        if v.is_nan() {
            out.push_str("nan");
            return;
        }
        if v.is_sign_negative() {
            out.push('-');
        }
        self.push_mag(out, v.abs());
    }

    /// Append `|v|` (non-negative).
    fn push_mag(&mut self, out: &mut String, a: f64) {
        if a == f64::INFINITY {
            out.push_str("inf");
            return;
        }
        if a == 0.0 {
            out.push('0');
            return;
        }
        // Fast path: the conventional shortest form.
        let conv = if (1e-4..1e16).contains(&a) {
            format!("{a}")
        } else {
            format!("{a:e}")
        };
        if cplex_read(&conv) == Some(a) {
            out.push_str(&conv);
            return;
        }
        out.push_str(&self.search(a, conv));
    }

    /// A string both CPLEX and a correctly rounded reader map to `a`, else
    /// one CPLEX maps to `a`, else the shortest exact one (`conv`).
    fn search(&mut self, a: f64, conv: String) -> String {
        let mut cands = Vec::new();
        // Same digits, other point placements: every one of them is the same
        // decimal, so a correctly rounded reader gives `a` for all.
        let (d, e10) = sci_digits(a, None);
        layouts(&d, e10, &mut cands);
        if let Some(s) = cands
            .iter()
            .filter(|s| cplex_read(s) == Some(a))
            .min_by_key(|s| s.len())
        {
            return s.clone();
        }
        // Longer digit strings near `a`.
        let mut cplex_only: Option<String> = None;
        for nd in 16..=19usize {
            let (d, e10) = sci_digits(a, Some(nd));
            let base: i128 = d.parse().expect("digits");
            let mut deltas: Vec<i128> = (-12..=12).collect();
            deltas.sort_by_key(|x| x.abs());
            for delta in deltas {
                let dd = (base + delta).to_string();
                if dd.len() != nd {
                    continue;
                }
                let exact = format!("{dd}e{}", e10 - nd as i32 + 1).parse::<f64>() == Ok(a);
                cands.clear();
                layouts(&dd, e10, &mut cands);
                for s in &cands {
                    if cplex_read(s) == Some(a) {
                        if exact {
                            return s.clone();
                        }
                        if cplex_only.is_none() {
                            cplex_only = Some(s.clone());
                        }
                    }
                }
            }
        }
        if let Some(s) = cplex_only {
            self.cplex_only += 1;
            return s;
        }
        self.unreadable += 1;
        conv
    }
}

// ------------------------------------------------------------------ names

/// Lower-case words that are section keywords or bound words in the CPLEX,
/// HiGHS or rusolve LP readers.
const RESERVED: &[&str] = &[
    "min",
    "minimize",
    "minimise",
    "minimum",
    "max",
    "maximize",
    "maximise",
    "maximum",
    "st",
    "s.t.",
    "st.",
    "subject",
    "such",
    "bound",
    "bounds",
    "gen",
    "general",
    "generals",
    "integer",
    "integers",
    "bin",
    "binary",
    "binaries",
    "free",
    "fre",
    "end",
    "semi",
    "semis",
    "semicontinuous",
    "semiinteger",
    "sos",
    "inf",
    "infinity",
    "in",
    "default",
    "nan",
];

/// A character every reader keeps inside a name.
fn name_char_ok(c: char) -> bool {
    if c.is_ascii() {
        c.is_ascii_alphanumeric() || "!\"#$%&'(),.?@_`{}~".contains(c)
    } else {
        !c.is_whitespace() && !c.is_control()
    }
}

/// The start of the name breaks a reader: a digit or `.` (a number), or for
/// a column `e`/`E` followed by nothing, a digit, `.` or `e`/`E` (CPLEX reads
/// an exponent), or a reserved word.
fn bad_start(s: &str, column: bool) -> bool {
    let b = s.as_bytes();
    if b[0].is_ascii_digit() || b[0] == b'.' {
        return true;
    }
    if column
        && (b[0] == b'e' || b[0] == b'E')
        && (b.len() == 1 || b[1].is_ascii_digit() || matches!(b[1], b'.' | b'e' | b'E'))
    {
        return true;
    }
    column && RESERVED.iter().any(|w| s.eq_ignore_ascii_case(w))
}

/// A comma is kept inside a name only between two other name characters
/// (the HiGHS-compatible rusolve tokenizer splits on any other comma).
fn bad_commas(s: &str) -> bool {
    s.starts_with(',') || s.ends_with(',') || s.contains(",,")
}

fn name_ok(s: &str, column: bool) -> bool {
    let max = if column { MAX_COL_NAME } else { MAX_ROW_NAME };
    !s.is_empty()
        && s.len() <= max
        && s.chars().all(name_char_ok)
        && !bad_commas(s)
        && !bad_start(s, column)
}

/// Longest prefix of `s` of at most `max` bytes, on a char boundary.
fn truncate(s: &str, max: usize) -> &str {
    if s.len() <= max {
        return s;
    }
    let mut k = max;
    while !s.is_char_boundary(k) {
        k -= 1;
    }
    &s[..k]
}

/// The legal form of `s` (`fallback` when empty), before de-duplication.
fn sanitize(s: &str, column: bool, fallback: impl FnOnce() -> String) -> String {
    if s.is_empty() {
        return fallback();
    }
    let mut t: String = s
        .chars()
        .map(|c| if name_char_ok(c) { c } else { '_' })
        .collect();
    while bad_commas(&t) {
        if t.starts_with(',') {
            t.replace_range(0..1, "_");
        }
        if t.ends_with(',') {
            let n = t.len();
            t.replace_range(n - 1..n, "_");
        }
        t = t.replace(",,", ",_");
    }
    if bad_start(&t, column) {
        t.insert(0, '_');
    }
    let max = if column { MAX_COL_NAME } else { MAX_ROW_NAME };
    truncate(&t, max).to_string()
}

/// Final names for one namespace: legal unique names stay; the others are
/// sanitized and suffixed `~k` until unique. `used` is shared so range
/// columns avoid model column names.
fn assign(
    raw: &[&str],
    column: bool,
    used: &mut HashSet<String>,
    fallback: impl Fn(usize) -> String,
    renamed: &mut usize,
) -> Vec<String> {
    let max = if column { MAX_COL_NAME } else { MAX_ROW_NAME };
    let mut out: Vec<Option<String>> = vec![None; raw.len()];
    for (k, s) in raw.iter().enumerate() {
        if name_ok(s, column) && !used.contains(*s) {
            used.insert(s.to_string());
            out[k] = Some(s.to_string());
        }
    }
    for (k, s) in raw.iter().enumerate() {
        if out[k].is_some() {
            continue;
        }
        let base = sanitize(s, column, || fallback(k));
        let mut cand = base.clone();
        let mut n = 0usize;
        while used.contains(&cand) || !name_ok(&cand, column) {
            n += 1;
            let suffix = format!("~{n}");
            cand = format!("{}{suffix}", truncate(&base, max - suffix.len()));
        }
        used.insert(cand.clone());
        out[k] = Some(cand);
        *renamed += 1;
    }
    out.into_iter().map(|s| s.expect("assigned")).collect()
}

fn lp_names_counted(m: &Model, report: &mut LpWriteReport) -> LpNames {
    let mut used_cols = HashSet::with_capacity(m.num_variables);
    let raw_cols: Vec<&str> = m.col_names.iter().map(String::as_str).collect();
    let columns = assign(
        &raw_cols,
        true,
        &mut used_cols,
        |j| format!("x{}", j + 1),
        &mut report.renamed_columns,
    );
    let mut used_rows = HashSet::with_capacity(m.num_constraints);
    let raw_rows: Vec<&str> = m.row_names.iter().map(String::as_str).collect();
    let rows = assign(
        &raw_rows,
        false,
        &mut used_rows,
        |i| format!("c{}", i + 1),
        &mut report.renamed_rows,
    );
    // The objective label may equal a row name (CPLEX accepts it).
    let objective = match m.obj_name.as_deref() {
        None | Some("") => None,
        Some(o) if name_ok(o, false) => Some(o.to_string()),
        Some(o) => {
            report.renamed_objective = true;
            Some(sanitize(o, false, || "obj".to_string()))
        }
    };
    let ranged: Vec<usize> = (0..m.num_constraints)
        .filter(|&i| m.ranges[i] != 0.0)
        .collect();
    report.ranged_rows = ranged.len();
    let rg_raw: Vec<String> = ranged.iter().map(|&i| format!("Rg{}", rows[i])).collect();
    let rg_refs: Vec<&str> = rg_raw.iter().map(String::as_str).collect();
    let mut not_counted = 0;
    let rg = assign(
        &rg_refs,
        true,
        &mut used_cols,
        |k| format!("Rg{}", k + 1),
        &mut not_counted,
    );
    let mut range_columns = vec![None; m.num_constraints];
    for (k, &i) in ranged.iter().enumerate() {
        range_columns[i] = Some(rg[k].clone());
    }
    LpNames {
        columns,
        rows,
        objective,
        range_columns,
    }
}

// ------------------------------------------------------------------ writer

struct Writer {
    out: String,
    line_len: usize,
    num: NumFmt,
}

impl Writer {
    fn line(&mut self, s: &str) {
        self.out.push_str(s);
        self.out.push('\n');
        self.line_len = 0;
    }

    /// Start a statement line (` label:` or empty).
    fn start(&mut self, label: Option<&str>) {
        let before = self.out.len();
        self.out.push(' ');
        if let Some(l) = label {
            self.out.push_str(l);
            self.out.push(':');
        }
        self.line_len = self.out.len() - before;
    }

    /// Append a piece that starts with a space, wrapping first when it would
    /// pass [`LINE_WRAP`].
    fn piece(&mut self, s: &str) {
        if self.line_len > 1 && self.line_len + s.len() > LINE_WRAP {
            self.out.push('\n');
            self.line_len = 0;
        }
        self.out.push_str(s);
        self.line_len += s.len();
    }

    fn end_line(&mut self) {
        self.out.push('\n');
        self.line_len = 0;
    }

    /// ` + c name` / ` - c name`, the magnitude CPLEX-exact.
    fn term(&mut self, c: f64, name: &str, buf: &mut String) {
        buf.clear();
        buf.push_str(if c.is_sign_negative() { " - " } else { " + " });
        self.num.push_mag(buf, c.abs());
        buf.push(' ');
        buf.push_str(name);
        let piece = std::mem::take(buf);
        self.piece(&piece);
        *buf = piece;
    }

    /// ` <rel> value` / a signed constant ` + v`.
    fn signed_piece(&mut self, prefix: &str, v: f64, buf: &mut String) {
        buf.clear();
        buf.push_str(prefix);
        self.num.push(buf, v);
        let piece = std::mem::take(buf);
        self.piece(&piece);
        *buf = piece;
    }

    fn num_string(&mut self, v: f64) -> String {
        let mut s = String::new();
        self.num.push(&mut s, v);
        s
    }

    fn model(&mut self, m: &Model, names: &LpNames) {
        let mut buf = String::new();
        let _ = writeln!(self.out, "\\ LP format for model {}", m.name);
        self.line(match m.obj_sense {
            ObjSense::Minimize => "Minimize",
            ObjSense::Maximize => "Maximize",
        });

        // Objective: every column in model order (zero cost included) so
        // CPLEX creates the columns in model order.
        self.start(names.objective.as_deref());
        for j in 0..m.num_variables {
            self.term(m.objective[j], &names.columns[j], &mut buf);
        }
        self.quadratic(m, names, &mut buf);
        if m.obj_offset != 0.0 {
            buf.clear();
            buf.push_str(if m.obj_offset < 0.0 { " - " } else { " + " });
            self.num.push_mag(&mut buf, m.obj_offset.abs());
            let piece = std::mem::take(&mut buf);
            self.piece(&piece);
            buf = piece;
        }
        self.end_line();

        // Rows, from the CSC by a counting transpose.
        self.line("Subject To");
        let a = &m.matrix;
        let (cp, ri, va) = (a.col_ptr(), a.row_idx(), a.values());
        let nr = m.num_constraints;
        let mut rp = vec![0usize; nr + 1];
        for &i in ri {
            rp[i as usize + 1] += 1;
        }
        for i in 0..nr {
            rp[i + 1] += rp[i];
        }
        let mut next = rp.clone();
        let mut cols = vec![0usize; ri.len()];
        let mut vals = vec![0.0f64; ri.len()];
        for j in 0..m.num_variables {
            for k in cp[j]..cp[j + 1] {
                let i = ri[k] as usize;
                cols[next[i]] = j;
                vals[next[i]] = va[k];
                next[i] += 1;
            }
        }
        let mut quad_of: Vec<Option<&crate::conic::QuadRow>> = vec![None; nr];
        for qr in &m.conic.quad_rows {
            quad_of[qr.row] = Some(qr);
        }
        for i in 0..nr {
            self.start(Some(&names.rows[i]));
            for k in rp[i]..rp[i + 1] {
                self.term(vals[k], &names.columns[cols[k]], &mut buf);
            }
            // `+ [ m x ^ 2 + m x * y ]`, monomial coefficients (no `/ 2`).
            if let Some(qr) = quad_of[i] {
                self.piece(" + [");
                for (a, b, mono) in qr.monomials() {
                    let name = if a == b {
                        format!("{} ^ 2", names.columns[a])
                    } else {
                        format!("{} * {}", names.columns[a], names.columns[b])
                    };
                    self.term(mono, &name, &mut buf);
                }
                self.piece(" ]");
            }
            let rel = if let Some(rg) = &names.range_columns[i] {
                self.term(-1.0, rg, &mut buf);
                " = "
            } else {
                match m.senses[i] {
                    RowSense::LessEqual => " <= ",
                    RowSense::GreaterEqual => " >= ",
                    RowSense::Equal => " = ",
                }
            };
            self.signed_piece(rel, m.rhs[i], &mut buf);
            self.end_line();
        }
        // Indicators: `name: z = v -> terms rel rhs` (Gurobi/CPLEX LP).
        for ind in &m.gc.indicators {
            let label = if name_ok(&ind.name, false) {
                Some(ind.name.as_str())
            } else {
                None
            };
            self.start(label);
            let head = format!(" {} = {} ->", names.columns[ind.bin], u8::from(ind.bin_val));
            self.piece(&head);
            for &(j, c) in &ind.terms {
                self.term(c, &names.columns[j], &mut buf);
            }
            let rel = match ind.sense {
                RowSense::LessEqual => " <= ",
                RowSense::GreaterEqual => " >= ",
                RowSense::Equal => " = ",
            };
            self.signed_piece(rel, ind.rhs, &mut buf);
            self.end_line();
        }

        // Bounds.
        let mut bounds = String::new();
        let mut binaries: Vec<usize> = Vec::new();
        let mut generals: Vec<usize> = Vec::new();
        for (j, b) in m.bounds.iter().enumerate() {
            let semi = m.gc.is_semi(j);
            let binary = !semi && b.integer && b.lower.to_bits() == 0 && b.upper == 1.0;
            if binary {
                binaries.push(j);
                continue;
            }
            if b.integer || m.gc.semi.iter().any(|s| s.col == j && s.integer) {
                generals.push(j);
            }
            self.bound(&mut bounds, &names.columns[j], b.lower, b.upper);
        }
        for i in 0..nr {
            if let Some(rg) = &names.range_columns[i] {
                let r = m.ranges[i];
                let (lo, hi) = match m.senses[i] {
                    RowSense::LessEqual => (-r.abs(), 0.0),
                    RowSense::GreaterEqual => (0.0, r.abs()),
                    RowSense::Equal if r > 0.0 => (0.0, r),
                    RowSense::Equal => (r, 0.0),
                };
                self.bound(&mut bounds, rg, lo, hi);
            }
        }
        if !bounds.is_empty() {
            self.line("Bounds");
            self.out.push_str(&bounds);
        }
        if !binaries.is_empty() {
            self.line("Binaries");
            for j in binaries {
                let _ = writeln!(self.out, " {}", names.columns[j]);
            }
        }
        if !generals.is_empty() {
            self.line("Generals");
            for j in generals {
                let _ = writeln!(self.out, " {}", names.columns[j]);
            }
        }
        if !m.gc.semi.is_empty() {
            self.line("Semi-continuous");
            for s in &m.gc.semi {
                let _ = writeln!(self.out, " {}", names.columns[s.col]);
            }
        }
        if !m.gc.sos.is_empty() {
            self.line("SOS");
            for s in &m.gc.sos {
                let pr = if s.priority != 0.0 {
                    format!(" {}", self.num_string(s.priority))
                } else {
                    String::new()
                };
                let mut l = format!(
                    " {}: {}{pr} ::",
                    s.name,
                    if s.kind == SosKind::Sos1 { "S1" } else { "S2" }
                );
                for &(j, w) in &s.members {
                    let ws = self.num_string(w);
                    let _ = write!(l, " {} : {}", names.columns[j], ws);
                }
                self.line(&l);
            }
        }
        if !m.gc.funcs.is_empty() || !m.gc.pwl.is_empty() {
            self.line("General Constraints");
            for f in &m.gc.funcs {
                let mut args: Vec<String> =
                    f.args.iter().map(|&j| names.columns[j].clone()).collect();
                if let Some(c) = f.constant {
                    args.push(self.num_string(c));
                }
                let l = format!(
                    " {}: {} = {} ( {} )",
                    f.name,
                    names.columns[f.res],
                    f.kind.word(),
                    args.join(" , ")
                );
                self.line(&l);
            }
            for p in &m.gc.pwl {
                let mut l = format!(
                    " {}: {} = PWL ( {} ) :",
                    p.name, names.columns[p.y], names.columns[p.x]
                );
                for (x, y) in p.xs.iter().zip(&p.ys) {
                    let (xs, ys) = (self.num_string(*x), self.num_string(*y));
                    let _ = write!(l, " ({xs}, {ys})");
                }
                self.line(&l);
            }
        }
        self.line("End");
    }

    /// `+ [ a x ^ 2 + b x * y ... ] / 2`: the Hessian `H` of
    /// `c'x + 1/2 x'Hx`. CPLEX reads `a x ^ 2` as `Q_xx = a` and `b x * y` as
    /// `Q_xy = Q_yx = b / 2`, so the diagonal is written as is and each
    /// off-diagonal pair as `H_xy + H_yx` (exactly `2 H_xy` when symmetric).
    fn quadratic(&mut self, m: &Model, names: &LpNames, buf: &mut String) {
        if m.hessian.is_empty() {
            return;
        }
        let mut q: BTreeMap<(usize, usize), f64> = BTreeMap::new();
        for &(r, c, v) in &m.hessian {
            *q.entry((r.min(c), r.max(c))).or_insert(0.0) += v;
        }
        q.retain(|_, v| *v != 0.0);
        if q.is_empty() {
            return;
        }
        self.piece(" + [");
        for (&(i, j), &v) in &q {
            let name = if i == j {
                format!("{} ^ 2", names.columns[i])
            } else {
                format!("{} * {}", names.columns[i], names.columns[j])
            };
            self.term(v, &name, buf);
        }
        self.piece(" ] / 2");
    }

    /// One variable's `Bounds` lines (nothing for the default `[0, +inf)`).
    fn bound(&mut self, out: &mut String, name: &str, lo: f64, hi: f64) {
        let inf = f64::INFINITY;
        let default_lo = lo.to_bits() == 0;
        if default_lo && hi == inf {
            return;
        }
        if lo == -inf && hi == inf {
            let _ = writeln!(out, " {name} free");
        } else if lo.to_bits() == hi.to_bits() {
            let v = self.num_string(lo);
            let _ = writeln!(out, " {name} = {v}");
        } else if default_lo && hi >= 0.0 {
            let h = self.num_string(hi);
            let _ = writeln!(out, " {name} <= {h}");
        } else if hi == inf {
            let l = self.num_string(lo);
            let _ = writeln!(out, " {name} >= {l}");
        } else if lo <= hi && !default_lo {
            let (l, h) = (self.num_string(lo), self.num_string(hi));
            let _ = writeln!(out, " {l} <= {name} <= {h}");
        } else {
            // An inverted pair or a negative upper bound over the default
            // lower: two lines, the lower one explicit (CPLEX keeps both
            // sides as written; a chain `lo <= x <= hi` with lo > hi is a
            // rusolve/HiGHS reader error).
            let (l, h) = (self.num_string(lo), self.num_string(hi));
            let _ = writeln!(out, " {name} >= {l}");
            let _ = writeln!(out, " {name} <= {h}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cplex_reader_known_values() {
        // Known CPLEX LP reader values for edge-case decimals.
        assert_eq!(cplex_read("1e30"), Some(9.999999999999999e29));
        assert_eq!(cplex_read("2.5e-308"), Some(0.0));
        assert_eq!(cplex_read("-0"), Some(-0.0));
        assert_eq!(cplex_read("0.1"), Some(0.1));
        assert_eq!(cplex_read("3"), Some(3.0));
        assert_eq!(cplex_read("x"), None);
        assert_eq!(cplex_read("1e"), None);
    }

    fn fmt(v: f64) -> (String, NumFmt) {
        let mut f = NumFmt::default();
        let mut s = String::new();
        f.push(&mut s, v);
        (s, f)
    }

    #[test]
    fn numbers_read_back_exactly_in_cplex_and_strtod() {
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut checked = 0;
        for k in 0..20_000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let mant = (x >> 11) as f64 / (1u64 << 53) as f64;
            let e = k % 41 - 20;
            let v = mant * 10f64.powi(e);
            if v == 0.0 {
                continue;
            }
            let (s, f) = fmt(v);
            assert_eq!(f.cplex_only + f.unreadable, 0, "{v:?} -> {s}");
            assert_eq!(cplex_read(&s), Some(v), "{v:?} -> {s}");
            assert_eq!(s.parse::<f64>(), Ok(v), "{v:?} -> {s}");
            checked += 1;
        }
        assert!(checked > 19_000);
        assert_eq!(fmt(-0.0).0, "-0");
        assert_eq!(fmt(f64::NEG_INFINITY).0, "-inf");
        assert_eq!(fmt(2.5).0, "2.5");
        assert_eq!(fmt(1e30).0.parse::<f64>(), Ok(1e30));
        assert_eq!(cplex_read(&fmt(1e30).0), Some(1e30));
    }

    #[test]
    fn names_are_renamed_only_when_needed() {
        assert!(name_ok("x(1,2)", true));
        assert!(name_ok("e_1", true));
        assert!(name_ok("e1", false), "row labels may look like exponents");
        assert!(!name_ok("e1", true));
        assert!(!name_ok("E.5", true));
        assert!(!name_ok("1x", true));
        assert!(!name_ok("a+b", true));
        assert!(!name_ok("x;y", true));
        assert!(!name_ok("x,", true));
        assert!(!name_ok("End", true));
        assert!(name_ok("End", false));
        assert!(!name_ok(&"y".repeat(255), true));
        assert!(name_ok(&"y".repeat(255), false));
        let mut used = HashSet::new();
        let mut n = 0;
        let got = assign(
            &["a", "a+b", "a_b", "", "1x", "a"],
            true,
            &mut used,
            |j| format!("x{}", j + 1),
            &mut n,
        );
        assert_eq!(got, ["a", "a_b~1", "a_b", "x4", "_1x", "a~1"]);
        assert_eq!(n, 4);
    }
}
