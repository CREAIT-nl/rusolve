//! AMPL `.nl` reader (linear and mixed-integer linear models).
//!
//! The format is what AMPL and Pyomo write for a solver. Sources:
//! * D. M. Gay, "Writing .nl Files", Sandia report SAND2005-7907P (2005):
//!   the 10-line text header, the segment kinds (F S V C L O d x r b k J G),
//!   the `r`/`b` line types (0 range, 1 upper, 2 lower, 3 free, 4 equal,
//!   5 complementarity), variable ordering (linear binaries then linear
//!   integers last), suffix kinds (`k & 3` entity, `k & 4` real-valued) and
//!   the binary variant ("binary int (4 byte) and double (8 byte) values ...
//!   a string consists of an int length followed by that many characters").
//! * The AMPL MP library's public reader (`include/mp/nl-reader.h`,
//!   `src/nl-reader.cc`, github.com/ampl/mp) for the binary details that
//!   paper leaves implicit: the segment key and the `r`/`b` line type are
//!   one ASCII byte each, `s`/`l` numeric nodes are a 2-/4-byte integer,
//!   `K` carries non-cumulative column counts, and header field `arith`
//!   names the byte order (1 IEEE little-endian, 2 IEEE big-endian; MP
//!   rejects any other value in a binary file, and so do we).
//!
//! Supported: text (`g`) and binary (`b`) files; one linear objective
//! (sense, coefficients and constant); linear constraints with any `r`
//! line type (ranges, equalities, one-sided, a constant expression body
//! moved into the bounds); variable bounds; binary and integer variables;
//! names from the `.col` / `.row` sidecar files AMPL writes with
//! `option auxfiles rc` (Pyomo with `symbolic_solver_labels`). Ignored on
//! purpose: initial primal/dual guesses (`x`, `d`), Jacobian column counts
//! (`k`, `K`) and suffixes that only hint (branching priority, basis
//! status, ...).
//!
//! Rejected with an error naming the feature (never read into a different
//! model): nonlinear constraints or objectives (header counts or any
//! operator/variable node in an expression), more than one objective,
//! complementarity constraints, logical constraints, defined variables
//! (common expressions), imported functions, SOS constraints (the
//! `sosno`/`ref` suffixes AMPL and Pyomo use for them), and binary files
//! whose arithmetic kind is not IEEE.
//!
//! Model conventions: constraints with no bound at all (`r` type 3) are
//! dropped, as the MPS reader drops extra `N` rows; without sidecar files
//! variables are named `x<j>` and constraints `c<i>` with the file's
//! 0-based indices.

use crate::model::{Model, ObjSense, RowSense, VariableBound};
use rusolve_core::{CscMatrix, Triplet};

/// A `.nl` read failure: where (text line or binary byte offset) and why.
#[derive(Debug, Clone, PartialEq)]
pub struct NlError {
    pub location: String,
    pub message: String,
}

impl std::fmt::Display for NlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, ".nl read error at {}: {}", self.location, self.message)
    }
}

impl std::error::Error for NlError {}

/// Contents of the `.col` / `.row` sidecar files (one name per line), when
/// they exist. `.row` lists the constraints, then the objectives.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NlNames {
    pub col: Option<Vec<String>>,
    pub row: Option<Vec<String>>,
}

fn unsupported(feature: &str) -> String {
    format!("{feature} not supported (rusolve reads linear and mixed-integer linear .nl models)")
}

/// Cursor over the segment part of the file, text or binary.
struct Src<'a> {
    b: &'a [u8],
    pos: usize,
    binary: bool,
    /// Binary only: the file's byte order differs from little-endian.
    big_endian: bool,
}

impl<'a> Src<'a> {
    fn err<T>(&self, message: impl Into<String>) -> Result<T, NlError> {
        let location = if self.binary {
            format!("byte {}", self.pos)
        } else {
            let line = 1 + self.b[..self.pos.min(self.b.len())]
                .iter()
                .filter(|&&c| c == b'\n')
                .count();
            format!("line {line}")
        };
        Err(NlError {
            location,
            message: message.into(),
        })
    }

    /// Text: skip blanks, newlines and `#` comments.
    fn skip_ws(&mut self) {
        while let Some(&c) = self.b.get(self.pos) {
            if c == b'#' {
                while self.pos < self.b.len() && self.b[self.pos] != b'\n' {
                    self.pos += 1;
                }
            } else if c.is_ascii_whitespace() {
                self.pos += 1;
            } else {
                break;
            }
        }
    }

    /// Text: the next whitespace-delimited token.
    fn token(&mut self) -> Result<&'a [u8], NlError> {
        self.skip_ws();
        let start = self.pos;
        while let Some(&c) = self.b.get(self.pos) {
            if c.is_ascii_whitespace() || c == b'#' {
                break;
            }
            self.pos += 1;
        }
        if start == self.pos {
            return self.err("unexpected end of file");
        }
        Ok(&self.b[start..self.pos])
    }

    fn bytes(&mut self, n: usize) -> Result<&'a [u8], NlError> {
        if self.b.len() - self.pos < n {
            return self.err("unexpected end of file");
        }
        let s = &self.b[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    /// The next single byte: a segment key, a bound line type, an
    /// expression node. `None` at end of file.
    fn byte(&mut self) -> Option<u8> {
        if !self.binary {
            self.skip_ws();
        }
        let c = *self.b.get(self.pos)?;
        self.pos += 1;
        Some(c)
    }

    fn int(&mut self) -> Result<i64, NlError> {
        if self.binary {
            let s: [u8; 4] = self.bytes(4)?.try_into().expect("4 bytes");
            Ok(if self.big_endian {
                i32::from_be_bytes(s)
            } else {
                i32::from_le_bytes(s)
            } as i64)
        } else {
            // Text integers may be glued to their key letter (`C12`), so
            // read digits only.
            self.skip_ws();
            let start = self.pos;
            if matches!(self.b.get(self.pos), Some(b'-' | b'+')) {
                self.pos += 1;
            }
            while self.b.get(self.pos).is_some_and(u8::is_ascii_digit) {
                self.pos += 1;
            }
            let s = std::str::from_utf8(&self.b[start..self.pos]).unwrap_or("");
            match s.parse::<i64>() {
                Ok(v) => Ok(v),
                Err(_) => self.err("expected an integer"),
            }
        }
    }

    fn short(&mut self) -> Result<i64, NlError> {
        if self.binary {
            let s: [u8; 2] = self.bytes(2)?.try_into().expect("2 bytes");
            Ok(if self.big_endian {
                i16::from_be_bytes(s)
            } else {
                i16::from_le_bytes(s)
            } as i64)
        } else {
            self.int()
        }
    }

    fn real(&mut self) -> Result<f64, NlError> {
        if self.binary {
            let s: [u8; 8] = self.bytes(8)?.try_into().expect("8 bytes");
            Ok(if self.big_endian {
                f64::from_be_bytes(s)
            } else {
                f64::from_le_bytes(s)
            })
        } else {
            let t = self.token()?;
            match std::str::from_utf8(t)
                .ok()
                .and_then(|s| s.parse::<f64>().ok())
            {
                Some(v) => Ok(v),
                None => self.err(format!(
                    "expected a number, found '{}'",
                    String::from_utf8_lossy(t)
                )),
            }
        }
    }

    /// A count or index in `0..limit`.
    fn index(&mut self, limit: usize, what: &str) -> Result<usize, NlError> {
        let v = self.int()?;
        if v < 0 || v as usize >= limit {
            return self.err(format!("{what} {v} out of range 0..{limit}"));
        }
        Ok(v as usize)
    }

    fn count(&mut self) -> Result<usize, NlError> {
        let v = self.int()?;
        if v < 0 {
            return self.err(format!("negative count {v}"));
        }
        Ok(v as usize)
    }

    fn name(&mut self) -> Result<String, NlError> {
        if self.binary {
            let n = self.count()?;
            Ok(String::from_utf8_lossy(self.bytes(n)?).into_owned())
        } else {
            Ok(String::from_utf8_lossy(self.token()?).into_owned())
        }
    }

    /// An `r`/`b` line type: one ASCII digit in both encodings.
    fn line_type(&mut self) -> Result<u8, NlError> {
        match self.byte() {
            Some(c @ b'0'..=b'5') => Ok(c - b'0'),
            Some(c) => self.err(format!("expected a bound type 0-5, found '{}'", c as char)),
            None => self.err("unexpected end of file"),
        }
    }

    /// An expression graph that must be a numeric constant (the nonlinear
    /// part of a linear constraint or objective): its value.
    fn constant_expr(&mut self, what: &str) -> Result<f64, NlError> {
        match self.byte() {
            Some(b'n') => self.real(),
            Some(b's') => Ok(self.short()? as f64),
            Some(b'l') => Ok(self.int()? as f64),
            Some(b'o') => {
                let op = self.int()?;
                self.err(unsupported(&format!(
                    "nonlinear expression in {what} (operator o{op})"
                )))
            }
            Some(b'v') => self.err(unsupported(&format!(
                "nonlinear expression in {what} (variable node)"
            ))),
            Some(b'f') => self.err(unsupported(&format!("imported function call in {what}"))),
            Some(b'h') => self.err(unsupported(&format!("string expression in {what}"))),
            Some(c) => self.err(format!("unexpected expression node '{}'", c as char)),
            None => self.err("unexpected end of file"),
        }
    }
}

/// The 10 header lines (always text).
struct Header {
    binary: bool,
    name: String,
    n_var: usize,
    n_con: usize,
    n_obj: usize,
    arith: i64,
    n_bin: usize,
    n_int: usize,
}

fn header_err<T>(line: usize, message: impl Into<String>) -> Result<T, NlError> {
    Err(NlError {
        location: format!("line {line}"),
        message: message.into(),
    })
}

/// Parse the header; returns it and the byte offset where segments start.
fn read_header(b: &[u8]) -> Result<(Header, usize), NlError> {
    let mut pos = 0;
    let mut lines: Vec<&[u8]> = Vec::with_capacity(10);
    for _ in 0..10 {
        let Some(len) = b[pos..].iter().position(|&c| c == b'\n') else {
            return header_err(lines.len() + 1, "file ends inside the 10-line header");
        };
        lines.push(&b[pos..pos + len]);
        pos += len + 1;
    }
    let binary = match lines[0].first() {
        Some(b'g') => false,
        Some(b'b') => true,
        _ => {
            return header_err(
                1,
                "not an .nl file: first byte must be 'g' (text) or 'b' (binary)",
            );
        }
    };
    // Line 1: `g<n> <options...> # problem <name>`; the options matter to
    // AMPL only.
    let line1 = String::from_utf8_lossy(lines[0]);
    let name = line1
        .split_once('#')
        .map(|(_, c)| c.trim())
        .map(|c| c.strip_prefix("problem").map(str::trim).unwrap_or(c))
        .unwrap_or("")
        .to_string();
    let ints = |k: usize| -> Result<Vec<i64>, NlError> {
        let s = String::from_utf8_lossy(lines[k]);
        let s = s.split('#').next().unwrap_or("");
        s.split_ascii_whitespace()
            .map(|t| {
                t.parse::<i64>().or_else(|_| {
                    header_err(k + 1, format!("header: expected an integer, found '{t}'"))
                })
            })
            .collect()
    };
    let at = |v: &[i64], i: usize| v.get(i).copied().unwrap_or(0);
    let need = |v: &[i64], n: usize, k: usize| -> Result<(), NlError> {
        if v.len() < n {
            return header_err(k + 1, format!("header line {} needs {n} integers", k + 1));
        }
        if v.iter().any(|&x| x < 0) {
            return header_err(k + 1, "header: negative count");
        }
        Ok(())
    };
    let l2 = ints(1)?;
    need(&l2, 3, 1)?;
    let l3 = ints(2)?;
    need(&l3, 2, 2)?;
    let l4 = ints(3)?;
    need(&l4, 2, 3)?;
    let l5 = ints(4)?;
    need(&l5, 2, 4)?;
    let l6 = ints(5)?;
    need(&l6, 2, 5)?;
    let l7 = ints(6)?;
    need(&l7, 2, 6)?;
    let l10 = ints(9)?;
    need(&l10, 0, 9)?;
    let fail = |k: usize, feature: &str| header_err(k + 1, unsupported(feature));
    let (n_var, n_con, n_obj) = (l2[0] as usize, l2[1] as usize, l2[2] as usize);
    if n_obj > 1 {
        return fail(1, &format!("multiple objectives ({n_obj})"));
    }
    if at(&l2, 5) > 0 {
        return fail(1, &format!("logical constraints ({})", at(&l2, 5)));
    }
    if l3[0] > 0 {
        return fail(2, &format!("nonlinear constraints ({})", l3[0]));
    }
    if l3[1] > 0 {
        return fail(2, "nonlinear objective");
    }
    if at(&l3, 2) + at(&l3, 3) > 0 {
        return fail(
            2,
            &format!("complementarity constraints ({})", at(&l3, 2) + at(&l3, 3)),
        );
    }
    if l4[0] > 0 {
        return fail(3, &format!("nonlinear network constraints ({})", l4[0]));
    }
    if l5[0] + l5[1] + at(&l5, 2) > 0 {
        return fail(4, "nonlinear variables");
    }
    if l6[1] > 0 {
        return fail(5, &format!("imported functions ({})", l6[1]));
    }
    if at(&l7, 2) + at(&l7, 3) + at(&l7, 4) > 0 {
        return fail(6, "nonlinear integer variables");
    }
    if l10.iter().sum::<i64>() > 0 {
        return fail(9, "defined variables (common expressions)");
    }
    let (n_bin, n_int) = (l7[0] as usize, l7[1] as usize);
    if n_bin + n_int > n_var {
        return header_err(
            7,
            format!("{n_bin} binary + {n_int} integer variables > {n_var}"),
        );
    }
    let arith = at(&l6, 2);
    if binary && arith != 1 && arith != 2 {
        return header_err(
            6,
            unsupported(&format!(
                "binary .nl with arithmetic kind {arith} (only 1 = IEEE little-endian, 2 = IEEE big-endian)"
            )),
        );
    }
    Ok((
        Header {
            binary,
            name,
            n_var,
            n_con,
            n_obj,
            arith,
            n_bin,
            n_int,
        },
        pos,
    ))
}

/// Suffix names that define SOS constraints (AMPL `.sosno`/`.ref`, also
/// what Pyomo writes for `SOSConstraint`).
const SOS_SUFFIXES: [&str; 4] = ["sosno", "ref", "sos", "sostype"];

/// Read an `.nl` file (text or binary) into a [`Model`]. `names` are the
/// sidecar `.col`/`.row` contents, if any (see [`NlNames`]).
pub fn read_nl(bytes: &[u8], names: &NlNames) -> Result<Model, NlError> {
    let (h, start) = read_header(bytes)?;
    let mut s = Src {
        b: bytes,
        pos: start,
        binary: h.binary,
        big_endian: h.arith == 2,
    };
    let (n, m) = (h.n_var, h.n_con);
    let mut con_const = vec![0.0f64; m];
    let mut con_lo = vec![f64::NEG_INFINITY; m];
    let mut con_hi = vec![f64::INFINITY; m];
    let mut have_r = false;
    let mut var_lo = vec![0.0f64; n];
    let mut var_hi = vec![0.0f64; n];
    let mut have_b = false;
    let mut obj_sense = ObjSense::Minimize;
    let mut obj_const = 0.0;
    let mut objective = vec![0.0f64; n];
    let mut triplets: Vec<Triplet> = Vec::new();

    while let Some(key) = s.byte() {
        match key {
            b'F' => return s.err(unsupported("imported functions (F segment)")),
            b'V' => return s.err(unsupported("defined variables (V segment)")),
            b'L' => return s.err(unsupported("logical constraints (L segment)")),
            b'S' => {
                let kind = s.int()?;
                let count = s.count()?;
                let name = s.name()?;
                if SOS_SUFFIXES.contains(&name.as_str()) {
                    return s.err(unsupported(&format!("SOS constraints (suffix .{name})")));
                }
                for _ in 0..count {
                    s.int()?;
                    if kind & 4 != 0 {
                        s.real()?;
                    } else {
                        s.int()?;
                    }
                }
            }
            b'C' => {
                let i = s.index(m, "constraint")?;
                con_const[i] = s.constant_expr(&format!("constraint {i}"))?;
            }
            b'O' => {
                let i = s.index(h.n_obj, "objective")?;
                let sigma = s.int()?;
                let c = s.constant_expr(&format!("objective {i}"))?;
                obj_sense = if sigma != 0 {
                    ObjSense::Maximize
                } else {
                    ObjSense::Minimize
                };
                obj_const = c;
            }
            b'd' | b'x' => {
                let (limit, what) = if key == b'd' {
                    (m, "constraint")
                } else {
                    (n, "variable")
                };
                for _ in 0..s.count()? {
                    s.index(limit, what)?;
                    s.real()?;
                }
            }
            b'r' => {
                for i in 0..m {
                    let (lo, hi) = match s.line_type()? {
                        0 => (s.real()?, s.real()?),
                        1 => (f64::NEG_INFINITY, s.real()?),
                        2 => (s.real()?, f64::INFINITY),
                        3 => (f64::NEG_INFINITY, f64::INFINITY),
                        4 => {
                            let v = s.real()?;
                            (v, v)
                        }
                        _ => {
                            return s.err(unsupported(&format!(
                                "complementarity constraints (constraint {i})"
                            )));
                        }
                    };
                    con_lo[i] = lo;
                    con_hi[i] = hi;
                }
                have_r = true;
            }
            b'b' => {
                for j in 0..n {
                    let (lo, hi) = match s.line_type()? {
                        0 => (s.real()?, s.real()?),
                        1 => (f64::NEG_INFINITY, s.real()?),
                        2 => (s.real()?, f64::INFINITY),
                        3 => (f64::NEG_INFINITY, f64::INFINITY),
                        4 => {
                            let v = s.real()?;
                            (v, v)
                        }
                        _ => return s.err(format!("bound type 5 is invalid for variable {j}")),
                    };
                    var_lo[j] = lo;
                    var_hi[j] = hi;
                }
                have_b = true;
            }
            b'k' | b'K' => {
                let c = s.count()?;
                if c + 1 != n.max(1) {
                    return s.err(format!(
                        "k segment has {c} entries, expected {}",
                        n.max(1) - 1
                    ));
                }
                for _ in 0..c {
                    s.int()?;
                }
            }
            b'J' => {
                let i = s.index(m, "constraint")?;
                for _ in 0..s.count()? {
                    let j = s.index(n, "variable")?;
                    let v = s.real()?;
                    triplets.push((i, j, v));
                }
            }
            b'G' => {
                let i = s.index(h.n_obj, "objective")?;
                for _ in 0..s.count()? {
                    let j = s.index(n, "variable")?;
                    let v = s.real()?;
                    if i == 0 {
                        objective[j] += v;
                    }
                }
            }
            c => {
                s.pos -= 1;
                return s.err(format!("unknown segment '{}'", c.escape_ascii()));
            }
        }
    }
    if !have_b && n > 0 {
        return s.err("missing variable bounds (b segment)");
    }
    if !have_r && m > 0 {
        return s.err("missing constraint bounds (r segment)");
    }

    // Names: sidecar files or generated, checked against the counts.
    let col_names = match &names.col {
        Some(v) => {
            let v: Vec<String> = trim_trailing_blank(v);
            if v.len() != n {
                return s.err(format!(
                    ".col file lists {} names for {n} variables (stale sidecar? delete it or rewrite both)",
                    v.len()
                ));
            }
            v
        }
        None => (0..n).map(|j| format!("x{j}")).collect(),
    };
    let (mut all_row_names, obj_name) = match &names.row {
        Some(v) => {
            let v: Vec<String> = trim_trailing_blank(v);
            if v.len() != m && v.len() != m + h.n_obj {
                return s.err(format!(
                    ".row file lists {} names for {m} constraints and {} objectives (stale sidecar? delete it or rewrite both)",
                    v.len(),
                    h.n_obj
                ));
            }
            let obj = v.get(m).cloned();
            (v.into_iter().take(m).collect::<Vec<_>>(), obj)
        }
        None => ((0..m).map(|i| format!("c{i}")).collect(), None),
    };

    // Rows: bounds shifted by the constant body, free rows dropped.
    let mut row_of = vec![usize::MAX; m];
    let mut senses = Vec::with_capacity(m);
    let mut rhs = Vec::with_capacity(m);
    let mut ranges = Vec::with_capacity(m);
    let mut row_names = Vec::with_capacity(m);
    for i in 0..m {
        let (lo, hi) = (con_lo[i] - con_const[i], con_hi[i] - con_const[i]);
        let (sense, b, r) = match (lo.is_finite(), hi.is_finite()) {
            (false, false) => continue,
            (true, true) if lo == hi => (RowSense::Equal, lo, 0.0),
            (true, true) if lo < hi => (RowSense::LessEqual, hi, hi - lo),
            (true, true) => {
                return s.err(format!(
                    "constraint {} ({}) has lower bound {lo} > upper bound {hi}",
                    i, all_row_names[i]
                ));
            }
            (false, true) => (RowSense::LessEqual, hi, 0.0),
            (true, false) => (RowSense::GreaterEqual, lo, 0.0),
        };
        row_of[i] = senses.len();
        senses.push(sense);
        rhs.push(b);
        ranges.push(r);
        row_names.push(std::mem::take(&mut all_row_names[i]));
    }
    let nr = senses.len();
    let triplets: Vec<Triplet> = triplets
        .into_iter()
        .filter(|&(i, _, _)| row_of[i] != usize::MAX)
        .map(|(i, j, v)| (row_of[i], j, v))
        .collect();
    let matrix = CscMatrix::from_triplets(nr, n, &triplets);

    // Variables: linear binaries then linear integers come last.
    let first_bin = n - h.n_bin - h.n_int;
    let first_int = n - h.n_int;
    let bounds = (0..n)
        .map(|j| {
            let (mut lower, mut upper) = (var_lo[j], var_hi[j]);
            if (first_bin..first_int).contains(&j) {
                lower = lower.max(0.0);
                upper = upper.min(1.0);
            }
            VariableBound {
                lower,
                upper,
                integer: j >= first_bin,
            }
        })
        .collect();

    Ok(Model {
        name: h.name,
        num_variables: n,
        num_constraints: nr,
        obj_sense,
        obj_name,
        obj_offset: obj_const,
        objective,
        senses,
        rhs,
        ranges,
        bounds,
        matrix,
        row_names,
        col_names,
        hessian: Vec::new(),
        gc: Default::default(),
        conic: Default::default(),
        cp: Default::default(),
    })
}

fn trim_trailing_blank(v: &[String]) -> Vec<String> {
    let mut v = v.to_vec();
    while v.last().is_some_and(|l| l.trim().is_empty()) {
        v.pop();
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pyomo 6.10's text writer output for the dialect-corpus model
    /// (`tests/dialects/pyomo.nl`).
    const PYOMO: &str = include_str!("../tests/dialects/pyomo.nl");

    fn header(body_lines: &[&str]) -> String {
        let mut s = String::new();
        for l in body_lines {
            s.push_str(l);
            s.push('\n');
        }
        s
    }

    fn read(text: &str) -> Result<Model, NlError> {
        read_nl(text.as_bytes(), &NlNames::default())
    }

    #[test]
    fn corpus_file_without_sidecars_gets_generated_names() {
        let m = read(PYOMO).unwrap();
        assert_eq!(m.name, "unknown");
        assert_eq!((m.num_variables, m.num_constraints), (6, 5));
        assert_eq!(m.col_names[0], "x0");
        assert_eq!(m.row_names[4], "c4");
        assert_eq!(m.obj_sense, ObjSense::Maximize);
        assert_eq!(m.obj_offset, 7.0);
        // b (index 4) binary, x (index 5) integer: the last nbv + niv.
        let ints: Vec<bool> = m.bounds.iter().map(|b| b.integer).collect();
        assert_eq!(ints, [false, false, false, false, true, true]);
    }

    const HEADER: [&str; 10] = [
        "g3 1 1 0",
        " 2 1 1 0 0",
        " 0 0",
        " 0 0",
        " 0 0 0",
        " 0 0 0 1",
        " 0 0 0 0 0",
        " 2 0",
        " 0 0",
        " 0 0 0 0 0",
    ];

    /// min 0 s.t. (r line) on x + y, x and y free; header line `k`
    /// (0-based) replaced when given.
    fn tiny_with(r_lines: &str, k: usize, line: &str) -> String {
        let mut h = HEADER;
        if !line.is_empty() {
            h[k] = line;
        }
        format!(
            "{}C0\nn0\nO0 0\nn0\nr\n{r_lines}b\n3\n3\nJ0 2\n0 1\n1 1\n",
            header(&h)
        )
    }

    fn tiny(r_lines: &str) -> String {
        tiny_with(r_lines, 0, "")
    }

    #[test]
    fn every_constraint_line_type_and_a_constant_body() {
        let row = |r: &str, c: &str| {
            let m = read(&tiny(r).replace("C0\nn0\n", &format!("C0\n{c}\n"))).unwrap();
            if m.num_constraints == 0 {
                return None;
            }
            Some((m.senses[0], m.rhs[0], m.ranges[0]))
        };
        assert_eq!(row("0 1 4\n", "n0"), Some((RowSense::LessEqual, 4.0, 3.0)));
        assert_eq!(row("1 4\n", "n0"), Some((RowSense::LessEqual, 4.0, 0.0)));
        assert_eq!(row("2 1\n", "n0"), Some((RowSense::GreaterEqual, 1.0, 0.0)));
        assert_eq!(row("3\n", "n0"), None, "free row dropped");
        assert_eq!(row("4 2\n", "n0"), Some((RowSense::Equal, 2.0, 0.0)));
        // body = x + y + 1.5, so 1 <= body <= 4 is -0.5 <= x + y <= 2.5
        assert_eq!(
            row("0 1 4\n", "n1.5"),
            Some((RowSense::LessEqual, 2.5, 3.0))
        );
        assert_eq!(
            row("2 1\n", "s2"),
            Some((RowSense::GreaterEqual, -1.0, 0.0))
        );
        assert_eq!(
            row("2 1\n", "l3"),
            Some((RowSense::GreaterEqual, -2.0, 0.0))
        );
    }

    #[test]
    fn unsupported_features_are_named() {
        let msg = |t: String| read(&t).unwrap_err().message;
        let t = tiny("2 1\n");
        assert!(
            msg(t.replace("C0\nn0\n", "C0\no2\nv0\nv1\n"))
                .contains("nonlinear expression in constraint 0")
        );
        assert!(
            msg(t.replace("O0 0\nn0\n", "O0 0\nv0\n"))
                .contains("nonlinear expression in objective 0")
        );
        assert!(msg(t.replace("r\n2 1\n", "r\n5 1 1\n")).contains("complementarity"));
        assert!(msg(t.replace(" 2 1 1 0 0\n", " 2 1 2 0 0\n")).contains("multiple objectives (2)"));
        assert!(msg(t.replace(" 2 1 1 0 0\n", " 2 1 1 0 0 1\n")).contains("logical constraints"));
        assert!(msg(tiny_with("2 1\n", 2, " 1 0")).contains("nonlinear constraints (1)"));
        assert!(msg(tiny_with("2 1\n", 5, " 0 1 0 1")).contains("imported functions"));
        assert!(msg(tiny_with("2 1\n", 6, " 0 0 1 0 0")).contains("nonlinear integer"));
        assert!(msg(tiny_with("2 1\n", 9, " 1 0 0 0 0")).contains("defined variables"));
        assert!(msg(tiny_with("2 1\n", 2, " 0 0 1 0 0 0")).contains("complementarity"));
        assert!(
            msg(t.replace("C0\nn0\n", "S0 1 sosno\n0 1\nC0\nn0\n"))
                .contains("SOS constraints (suffix .sosno)")
        );
        assert!(msg(t.replace("C0\nn0\n", "V2 0 0\nn0\nC0\nn0\n")).contains("defined variables"));
        assert!(msg(t.replacen("g3", "b3", 1)).contains("arithmetic kind 0"));
        // hint suffixes are skipped
        let m = read(&t.replace(
            "C0\nn0\n",
            "S0 2 priority\n0 5\n1 7\nS4 1 zork\n0 1.5\nC0\nn0\n",
        ))
        .unwrap();
        assert_eq!(m.num_constraints, 1);
    }

    #[test]
    fn sidecar_names_are_used_and_checked() {
        let t = tiny("2 1\n");
        let names = NlNames {
            col: Some(vec!["a".into(), "b".into(), String::new()]),
            row: Some(vec!["r0".into(), "obj".into()]),
        };
        let m = read_nl(t.as_bytes(), &names).unwrap();
        assert_eq!(m.col_names, ["a", "b"]);
        assert_eq!(m.row_names, ["r0"]);
        assert_eq!(m.obj_name.as_deref(), Some("obj"));
        let stale = NlNames {
            col: Some(vec!["a".into()]),
            row: None,
        };
        assert!(
            read_nl(t.as_bytes(), &stale)
                .unwrap_err()
                .message
                .contains(".col file lists 1 names")
        );
    }

    #[test]
    fn truncated_file_is_an_error() {
        let t = tiny("2 1\n");
        let cut = &t[..t.len() - 4];
        assert!(read(cut).is_err());
        assert!(read(&t[..40]).unwrap_err().message.contains("header"));
    }
}
