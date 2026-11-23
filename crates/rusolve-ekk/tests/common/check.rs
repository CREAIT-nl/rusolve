//! Bit-exact comparison helpers. No tolerances anywhere: the oracle is built
//! with `-ffp-contract=off`, so a faithful port reproduces every f64 bit
//! pattern (including the sign of zero and `kHighsZero`).
//!
//! Every failure message names the location (`fixture:line [type].field`),
//! the first differing index, and expected vs actual as bits + decimal.

use std::fmt::Write as _;

use serde_json::Value;

use super::emit::{HVectorView, MatrixView};
use super::fixture::short;
use super::hash::hash_hex;
use super::records::{HVectorRec, MatrixRec};
use super::{K_HIGHS_ZERO_BITS, hex64};

/// How many differing elements a slice report lists after the first one.
const MORE_DIFFS: usize = 4;

/// `0x3ff0000000000000 (1.0)`; marks `-0.0`, `kHighsZero`, NaN payloads.
pub fn fmt_f64(x: f64) -> String {
    let bits = x.to_bits();
    let note = if bits == K_HIGHS_ZERO_BITS {
        " = kHighsZero"
    } else if bits == 0x8000_0000_0000_0000 {
        " = -0.0"
    } else {
        ""
    };
    format!("{bits:016x} ({x:e}{note})")
}

/// Distance in units in the last place between two finite f64s of any sign
/// (as ordered integers); `None` if either is NaN. Diagnostic only: it says
/// whether a mismatch is a rounding-order slip (1-2 ulp) or a logic bug.
pub fn ulp_distance(a: f64, b: f64) -> Option<u64> {
    if a.is_nan() || b.is_nan() {
        return None;
    }
    fn ordered(x: f64) -> i128 {
        let b = x.to_bits() as i64;
        if b < 0 {
            i128::from(i64::MIN) - i128::from(b)
        } else {
            i128::from(b)
        }
    }
    Some((ordered(a) - ordered(b)).unsigned_abs() as u64)
}

fn f64_diff_note(e: f64, a: f64) -> String {
    match ulp_distance(e, a) {
        Some(u) => format!("{u} ulp apart"),
        None => "NaN involved".into(),
    }
}

/// Panic unless `actual` has exactly the bit pattern of `expected`.
#[track_caller]
pub fn assert_f64_bits(ctx: &str, expected: f64, actual: f64) {
    if expected.to_bits() != actual.to_bits() {
        panic!(
            "{ctx}: f64 mismatch\n  expected {}\n  actual   {}\n  ({})",
            fmt_f64(expected),
            fmt_f64(actual),
            f64_diff_note(expected, actual)
        );
    }
}

/// Describe the differences of two slices, `None` if equal.
fn slice_report<T: PartialEq>(
    expected: &[T],
    actual: &[T],
    eq: impl Fn(&T, &T) -> bool,
    fmt: impl Fn(&T) -> String,
    extra: impl Fn(&T, &T) -> String,
) -> Option<String> {
    let n = expected.len().min(actual.len());
    let diffs: Vec<usize> = (0..n).filter(|&i| !eq(&expected[i], &actual[i])).collect();
    if diffs.is_empty() && expected.len() == actual.len() {
        return None;
    }
    let mut s = String::new();
    if expected.len() != actual.len() {
        let _ = write!(
            s,
            "length mismatch: expected {}, actual {}",
            expected.len(),
            actual.len()
        );
        if diffs.is_empty() {
            let _ = write!(s, " (common prefix of {n} equal)");
        }
        s.push('\n');
    }
    if let Some(&first) = diffs.first() {
        let _ = writeln!(
            s,
            "{} of {n} compared elements differ; first at [{first}]",
            diffs.len()
        );
        for &i in diffs.iter().take(1 + MORE_DIFFS) {
            let _ = writeln!(
                s,
                "  [{i}] expected {}  actual {}{}",
                fmt(&expected[i]),
                fmt(&actual[i]),
                extra(&expected[i], &actual[i])
            );
        }
        if diffs.len() > 1 + MORE_DIFFS {
            let _ = writeln!(s, "  ... {} more", diffs.len() - 1 - MORE_DIFFS);
        }
    }
    Some(s.trim_end().to_string())
}

fn f64_slice_report(expected: &[f64], actual: &[f64]) -> Option<String> {
    slice_report(
        expected,
        actual,
        |a, b| a.to_bits() == b.to_bits(),
        |x| fmt_f64(*x),
        |e, a| format!("  ({})", f64_diff_note(*e, *a)),
    )
}

fn plain_slice_report<T: PartialEq + std::fmt::Display>(
    expected: &[T],
    actual: &[T],
) -> Option<String> {
    slice_report(
        expected,
        actual,
        |a, b| a == b,
        |x| x.to_string(),
        |_, _| String::new(),
    )
}

/// Bit-exact f64 slice comparison.
#[track_caller]
pub fn assert_f64_slice(ctx: &str, expected: &[f64], actual: &[f64]) {
    if let Some(r) = f64_slice_report(expected, actual) {
        panic!("{ctx}: f64 slice mismatch\n{r}");
    }
}

/// HighsInt slice comparison.
#[track_caller]
pub fn assert_int_slice(ctx: &str, expected: &[i32], actual: &[i32]) {
    if let Some(r) = plain_slice_report(expected, actual) {
        panic!("{ctx}: HighsInt slice mismatch\n{r}");
    }
}

/// int8 slice comparison (`nonbasicFlag_`, `nonbasicMove_`).
#[track_caller]
pub fn assert_i8_slice(ctx: &str, expected: &[i8], actual: &[i8]) {
    if let Some(r) = plain_slice_report(expected, actual) {
        panic!("{ctx}: int8 slice mismatch\n{r}");
    }
}

/// One difference found by [`Checker`] or [`diff_values`].
#[derive(Clone, Debug)]
pub struct Mismatch {
    /// `fixture:line [type].field[.sub]`
    pub path: String,
    pub expected: String,
    pub actual: String,
    pub note: String,
}

impl std::fmt::Display for Mismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\n    expected {}\n    actual   {}",
            self.path, self.expected, self.actual
        )?;
        if !self.note.is_empty() {
            write!(f, "\n    {}", self.note.replace('\n', "\n    "))?;
        }
        Ok(())
    }
}

/// Collects every field mismatch of one record (or one test step) and panics
/// once, in [`Checker::finish`], listing all of them. Use it so a single run
/// shows the whole picture rather than only the first wrong field.
///
/// ```ignore
/// let mut c = Checker::new(r.ctx.clone());
/// c.int("count", r.int("count"), actual.count);
/// c.f64s("array", &expected, &actual.array);
/// c.finish();
/// ```
#[derive(Debug)]
pub struct Checker {
    pub ctx: String,
    pub failures: Vec<Mismatch>,
}

impl Checker {
    pub fn new(ctx: impl Into<String>) -> Self {
        Checker {
            ctx: ctx.into(),
            failures: Vec::new(),
        }
    }

    fn path(&self, field: &str) -> String {
        format!("{}.{}", self.ctx, field)
    }

    pub fn push(&mut self, field: &str, expected: String, actual: String, note: String) {
        let path = self.path(field);
        self.failures.push(Mismatch {
            path,
            expected,
            actual,
            note,
        });
    }

    pub fn ok(&self) -> bool {
        self.failures.is_empty()
    }

    pub fn f64(&mut self, field: &str, expected: f64, actual: f64) -> &mut Self {
        if expected.to_bits() != actual.to_bits() {
            self.push(
                field,
                fmt_f64(expected),
                fmt_f64(actual),
                f64_diff_note(expected, actual),
            );
        }
        self
    }

    pub fn int(
        &mut self,
        field: &str,
        expected: impl Into<i64>,
        actual: impl Into<i64>,
    ) -> &mut Self {
        let (e, a) = (expected.into(), actual.into());
        if e != a {
            self.push(field, e.to_string(), a.to_string(), String::new());
        }
        self
    }

    pub fn bool(&mut self, field: &str, expected: bool, actual: bool) -> &mut Self {
        if expected != actual {
            self.push(
                field,
                expected.to_string(),
                actual.to_string(),
                String::new(),
            );
        }
        self
    }

    pub fn str(&mut self, field: &str, expected: &str, actual: &str) -> &mut Self {
        if expected != actual {
            self.push(
                field,
                format!("{expected:?}"),
                format!("{actual:?}"),
                String::new(),
            );
        }
        self
    }

    /// FNV-1a hash comparison (see [`super::hash`]).
    pub fn hash(&mut self, field: &str, expected: u64, actual: u64) -> &mut Self {
        if expected != actual {
            self.push(
                field,
                hash_hex(expected),
                hash_hex(actual),
                "hash mismatch: compare the full array at the nearest record that dumps it".into(),
            );
        }
        self
    }

    pub fn f64s(&mut self, field: &str, expected: &[f64], actual: &[f64]) -> &mut Self {
        if let Some(r) = f64_slice_report(expected, actual) {
            self.push(
                field,
                format!("{} f64", expected.len()),
                format!("{} f64", actual.len()),
                r,
            );
        }
        self
    }

    pub fn ints(&mut self, field: &str, expected: &[i32], actual: &[i32]) -> &mut Self {
        if let Some(r) = plain_slice_report(expected, actual) {
            self.push(
                field,
                format!("{} int", expected.len()),
                format!("{} int", actual.len()),
                r,
            );
        }
        self
    }

    pub fn i8s(&mut self, field: &str, expected: &[i8], actual: &[i8]) -> &mut Self {
        if let Some(r) = plain_slice_report(expected, actual) {
            self.push(
                field,
                format!("{} int8", expected.len()),
                format!("{} int8", actual.len()),
                r,
            );
        }
        self
    }

    /// Compare an HVector against its fixture dump exactly as the oracle sees
    /// it: `size`, `count`, `index[0..count)` (order matters), `synthetic_tick`,
    /// the full dense `array` (bit patterns, so a stray `-0.0` or
    /// `kHighsZero` fails), and `pack*` when the fixture has them.
    pub fn hv(&mut self, field: &str, expected: &HVectorRec, actual: &HVectorView) -> &mut Self {
        self.int(&format!("{field}.size"), expected.size, actual.size);
        self.int(&format!("{field}.count"), expected.count, actual.count);
        if let Some(ei) = &expected.index {
            let n = actual.count.max(0) as usize;
            let ai = actual.index.get(..n).unwrap_or(actual.index);
            self.ints(&format!("{field}.index"), ei, ai);
        }
        self.f64(
            &format!("{field}.synthetic_tick"),
            expected.synthetic_tick,
            actual.synthetic_tick,
        );
        self.f64s(&format!("{field}.array"), &expected.array(), actual.array);
        if let Some(ep) = &expected.pack {
            match &actual.pack {
                None => self.push(
                    &format!("{field}.pack"),
                    "pack fields".into(),
                    "none".into(),
                    String::new(),
                ),
                Some(ap) => {
                    let n = ap.pack_count.max(0) as usize;
                    self.bool(&format!("{field}.packFlag"), ep.pack_flag, ap.pack_flag);
                    self.int(&format!("{field}.packCount"), ep.pack_count, ap.pack_count);
                    self.ints(
                        &format!("{field}.packIndex"),
                        &ep.pack_index,
                        ap.pack_index.get(..n).unwrap_or(ap.pack_index),
                    );
                    self.f64s(
                        &format!("{field}.packValue"),
                        &ep.pack_value,
                        ap.pack_value.get(..n).unwrap_or(ap.pack_value),
                    );
                }
            }
        }
        self
    }

    /// Compare a HighsSparseMatrix: format, dimensions, start, p_end (format
    /// 3), index and value (bit-exact, in stored order).
    pub fn mat(&mut self, field: &str, expected: &MatrixRec, actual: &MatrixView) -> &mut Self {
        self.int(&format!("{field}.format"), expected.format, actual.format);
        self.int(
            &format!("{field}.num_col"),
            expected.num_col,
            actual.num_col,
        );
        self.int(
            &format!("{field}.num_row"),
            expected.num_row,
            actual.num_row,
        );
        self.ints(&format!("{field}.start"), &expected.start, actual.start);
        if let Some(ep) = &expected.p_end {
            self.ints(&format!("{field}.p_end"), ep, actual.p_end.unwrap_or(&[]));
        }
        self.ints(&format!("{field}.index"), &expected.index, actual.index);
        self.f64s(&format!("{field}.value"), &expected.value, actual.value);
        self
    }

    /// Encoding-aware diff of a whole JSON value (e.g. a fixture record vs
    /// the `Obj` the port built), see [`diff_values`].
    pub fn value(&mut self, field: &str, expected: &Value, actual: &Value) -> &mut Self {
        let path = if field.is_empty() {
            self.ctx.clone()
        } else {
            self.path(field)
        };
        self.failures.extend(diff_values(&path, expected, actual));
        self
    }

    /// Panic listing every mismatch, if any.
    #[track_caller]
    pub fn finish(&self) {
        if self.failures.is_empty() {
            return;
        }
        let mut s = format!("{}: {} mismatch(es)\n", self.ctx, self.failures.len());
        for m in &self.failures {
            let _ = writeln!(s, "  {m}");
        }
        panic!("{}", s.trim_end());
    }
}

fn is_hash_key(key: &str) -> bool {
    key.starts_with("h_")
}

/// Render a JSON scalar for a message, decoding hex16 strings as f64 (or as a
/// hash when `as_hash`).
fn render(v: &Value, as_hash: bool) -> String {
    match v.as_str().map(hex64) {
        Some(Ok(bits)) if as_hash => format!("hash {bits:016x}"),
        Some(Ok(bits)) => fmt_f64(f64::from_bits(bits)),
        _ => short(v),
    }
}

/// Diff two JSON values in the oracle encoding, returning every mismatch.
///
/// - 16-hex-digit strings are compared as bit patterns and shown as f64
///   (bits + decimal), or as hashes for `h_*` keys and inside `hash` objects.
/// - Missing and extra object keys are reported (the oracle's conditional
///   fields, e.g. `edge_weight_row_out_after`, must be present exactly when
///   HiGHS emits them).
/// - Arrays report their length mismatch and the first differing elements.
pub fn diff_values(path: &str, expected: &Value, actual: &Value) -> Vec<Mismatch> {
    let mut out = Vec::new();
    diff_rec(path, expected, actual, false, &mut out);
    out
}

fn diff_rec(path: &str, e: &Value, a: &Value, as_hash: bool, out: &mut Vec<Mismatch>) {
    match (e, a) {
        (Value::Object(em), Value::Object(am)) => {
            for (k, ev) in em {
                let p = format!("{path}.{k}");
                let h = as_hash || is_hash_key(k) || k == "hash";
                match am.get(k) {
                    Some(av) => diff_rec(&p, ev, av, h, out),
                    None => out.push(Mismatch {
                        path: p,
                        expected: render(ev, h),
                        actual: "<missing>".into(),
                        note: String::new(),
                    }),
                }
            }
            for (k, av) in am {
                if !em.contains_key(k) {
                    out.push(Mismatch {
                        path: format!("{path}.{k}"),
                        expected: "<absent>".into(),
                        actual: render(av, as_hash || is_hash_key(k)),
                        note: "the oracle does not emit this field here".into(),
                    });
                }
            }
        }
        (Value::Array(ev), Value::Array(av)) => {
            let scalar = ev.iter().chain(av).all(|x| !x.is_object() && !x.is_array());
            if scalar {
                let diffs: Vec<usize> = (0..ev.len().min(av.len()))
                    .filter(|&i| !scalar_eq(&ev[i], &av[i]))
                    .collect();
                if diffs.is_empty() && ev.len() == av.len() {
                    return;
                }
                let mut note = String::new();
                if ev.len() != av.len() {
                    let _ = writeln!(
                        note,
                        "length mismatch: expected {}, actual {}",
                        ev.len(),
                        av.len()
                    );
                }
                if let Some(&first) = diffs.first() {
                    let _ = writeln!(
                        note,
                        "{} element(s) differ; first at [{first}]",
                        diffs.len()
                    );
                    for &i in diffs.iter().take(1 + MORE_DIFFS) {
                        let _ = writeln!(
                            note,
                            "[{i}] expected {}  actual {}",
                            render(&ev[i], as_hash),
                            render(&av[i], as_hash)
                        );
                    }
                    if diffs.len() > 1 + MORE_DIFFS {
                        let _ = writeln!(note, "... {} more", diffs.len() - 1 - MORE_DIFFS);
                    }
                }
                out.push(Mismatch {
                    path: path.to_string(),
                    expected: format!("array of {}", ev.len()),
                    actual: format!("array of {}", av.len()),
                    note: note.trim_end().to_string(),
                });
            } else {
                if ev.len() != av.len() {
                    out.push(Mismatch {
                        path: path.to_string(),
                        expected: format!("array of {}", ev.len()),
                        actual: format!("array of {}", av.len()),
                        note: String::new(),
                    });
                }
                for (i, (x, y)) in ev.iter().zip(av).enumerate() {
                    diff_rec(&format!("{path}[{i}]"), x, y, as_hash, out);
                }
            }
        }
        _ => {
            if !scalar_eq(e, a) {
                let note = match (hexval(e), hexval(a)) {
                    (Some(x), Some(y)) if !as_hash => {
                        f64_diff_note(f64::from_bits(x), f64::from_bits(y))
                    }
                    _ => String::new(),
                };
                out.push(Mismatch {
                    path: path.to_string(),
                    expected: render(e, as_hash),
                    actual: render(a, as_hash),
                    note,
                });
            }
        }
    }
}

fn hexval(v: &Value) -> Option<u64> {
    v.as_str().and_then(|s| hex64(s).ok())
}

fn scalar_eq(e: &Value, a: &Value) -> bool {
    match (hexval(e), hexval(a)) {
        (Some(x), Some(y)) => x == y,
        _ => e == a,
    }
}
