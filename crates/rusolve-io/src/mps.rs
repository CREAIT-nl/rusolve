//! Free-format MPS reader and writer.
//!
//! Reader covers: `NAME`, `OBJSENSE`, `ROWS` (`N`/`L`/`G`/`E`), `COLUMNS`
//! (including `INTORG`/`INTEND` markers), `RHS` (first-set-only entries; an
//! entry on the objective `N` row becomes [`Model::obj_offset`]), `RANGES`,
//! `BOUNDS` (`UP`, `LO`, `FX`, `FR`, `MI`, `PL`, `BV`, `LI`, `UI`;
//! first-set-only; `SC`/`SI` give semi-continuous/semi-integer columns an
//! upper bound, stored in [`Model::gc`]), the sections `INDICATORS`
//! (`IF row binvar 0|1`) and `SOS` (`S1`/`name`, then `col weight`), the quadratic objective sections `QUADOBJ`,
//! `QMATRIX` and `QSECTION` (into [`Model::hessian`]) and `ENDATA`.
//! `QCMATRIX` (quadratic constraints, monomial coefficients without a half)
//! and `CSECTION` (`QUAD` / `RQUAD` cones; other cone types are an explicit
//! error) fill [`Model::conic`]; the
//! LP/MIP entry points ([`read_mps`], [`read_mps_file_order`]) reject a file
//! that has either, [`read_mps_conic`] keeps them.
//!
//! Bound-type semantics follow HiGHS's free-format reader (MIT,
//! `reference/HiGHS/highs/io/HMpsFF.cpp` `parseBounds`, :1340–1537):
//! `LI v` = integer with lower bound `v`, `UI v` = integer with upper bound
//! `v` (the other bound is left as it stands: MPS default `0` / `+∞`, or an
//! earlier directive), `PL` = upper bound `+∞`, `MI` = lower bound `−∞`
//! (upper untouched), `BV` = integer `[0, 1]`, `FR` = `(−∞, +∞)`. Like HiGHS
//! (and unlike some older readers) a negative `UP`/`UI` value does NOT move
//! the default lower bound `0` to `−∞`: `UP x −5` with no `LO`/`MI` is an
//! infeasible column, exactly as HiGHS reads it.
//! Writer emits free-format MPS for any [`Model`], so that
//! `read_mps(&write_mps(&m))` round-trips semantically.
//!
//! Provenance: MPS format as described in the public IBM COIN-OR MPS
//! documentation and HiGHS user documentation. `SC`/`SI`,
//! `INDICATORS` and `SOS` follow the public Gurobi "MPS format" page
//! and HiGHS `HMpsFF.cpp` (`parseSos`).
//!
//! The quadratic sections are translated from HiGHS
//! `io/HMpsFF.cpp`: section keywords (lines 418-447), dispatch (296-308),
//! `parseHessian`, `parseQuadRows` and `parseQuadMatrix` (1735-1889).
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use crate::conic::{Cone, ConeKind, ConicBlock, QuadRow};
use std::collections::HashMap as StdHashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// FxHash (rustc's multiply-rotate hasher). The name tables are only probed
/// and inserted into, never iterated, so the hash function cannot change the
/// parse; SipHash cost 0.1+ s on a 1.2M-line file (measured).
#[derive(Default, Clone, Copy)]
struct FxHasher(u64);

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        const K: u64 = 0x517c_c1b7_2722_0a95;
        let mut h = self.0;
        let mut b = bytes;
        while b.len() >= 8 {
            let w = u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
            h = (h.rotate_left(5) ^ w).wrapping_mul(K);
            b = &b[8..];
        }
        if b.len() >= 4 {
            let w = u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as u64;
            h = (h.rotate_left(5) ^ w).wrapping_mul(K);
            b = &b[4..];
        }
        for &x in b {
            h = (h.rotate_left(5) ^ x as u64).wrapping_mul(K);
        }
        self.0 = h;
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        const K: u64 = 0x517c_c1b7_2722_0a95;
        self.0 = (self.0.rotate_left(5) ^ i as u64).wrapping_mul(K);
    }
    #[inline]
    fn finish(&self) -> u64 {
        // The low bits of a multiply depend only on the low bits of its inputs, and hashbrown
        // indexes with the low bits: numbered names (x1234567) clustered so badly that reading
        // dlr1 (9M columns) took 155 s instead of about 20. rustc-hash 2 rotates the same way.
        self.0.rotate_left(26)
    }
}

type HashMap<K, V> = StdHashMap<K, V, BuildHasherDefault<FxHasher>>;
use std::fmt::Write as _;

use rusolve_core::Triplet;

use crate::gc::{GcModel, Indicator, SemiVar, SosKind, SosSet};
use crate::model::{Model, ObjSense, RowSense, VariableBound};

/// A parse failure with the offending 1-based line number.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MpsError {
    pub line: usize,
    pub message: String,
}

impl std::fmt::Display for MpsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MPS parse error at line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for MpsError {}

fn err<T>(line: usize, message: impl Into<String>) -> Result<T, MpsError> {
    Err(MpsError {
        line,
        message: message.into(),
    })
}

fn parse_f64(s: &str, line: usize, what: &str) -> Result<f64, MpsError> {
    s.parse::<f64>().map_err(|_| MpsError {
        line,
        message: format!("bad {what} value '{s}'"),
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Section {
    None,
    Name,
    ObjSense,
    Rows,
    Columns,
    Rhs,
    Ranges,
    Bounds,
    /// `QUADOBJ`, `QMATRIX` or a `QSECTION` on a known row
    /// (`parseQuadMatrix`).
    Quad(QuadSection),
    /// A `QSECTION` on an undefined or free row: lines are read until the
    /// next section (`parseQuadRows`, HMpsFF.cpp:1764-1789).
    QuadSkip,
    /// `INDICATORS`: `IF row binvar value` (CPLEX/Gurobi MPS extension).
    Indicators,
    /// `SOS`: `S1|S2 [SOS] name [priority]`, then member lines.
    Sos,
    /// `QCMATRIX <row>`: entries of the quadratic part of constraint row `.0`.
    QcMatrix(usize),
    /// `CSECTION`: one column name per line, for the last cone opened.
    Cone,
}

/// The `parseQuadMatrix` arguments of one quadratic section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct QuadSection {
    /// `section_name`.
    name: &'static str,
    /// `triangular`: `QUADOBJ` / `QSECTION` give one triangle.
    triangular: bool,
    /// Which `q_entries` the section fills.
    dest: QuadDest,
}

/// `parseQuadRows`: `rowidx == -1 ? q_entries : qrows_entries[rowidx]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QuadDest {
    /// The objective Hessian (`q_entries`, [`Model::hessian`]).
    Objective,
    /// `qrows_entries[row]` of a constraint row (parsed, not used).
    Row(usize),
}

/// HiGHS: HMpsFF.cpp::HMpsFF::getColIdx (lines 491-506), with the default
/// `add_if_new = true` that `parseQuadMatrix` uses.
///
/// A new column is continuous with the default bounds `[0, +inf)`
/// (`VariableBound::default()`, which the bounds pass gives every column).
fn get_col_idx<'a>(
    colname: &'a str,
    col_index: &mut HashMap<&'a str, usize>,
    col_names: &mut Vec<String>,
    objective: &mut Vec<f64>,
    marker_binary: &mut Vec<bool>,
) -> usize {
    // look up column name
    if let Some(&j) = col_index.get(colname) {
        return j;
    }
    // add new continuous column with default bounds
    col_index.insert(colname, col_names.len());
    col_names.push(colname.to_string());
    objective.push(0.0);
    marker_binary.push(false);
    col_names.len() - 1
}

/// A BOUNDS directive, applied in order after all columns are known.
#[derive(Debug, Clone, Copy)]
enum BoundOp {
    Upper(f64),
    Lower(f64),
    Fixed(f64),
    Free,
    MinusInf,
    /// `PL`: upper bound `+∞` (redundant with the default unless an earlier
    /// directive set one).
    PlusInf,
    Binary,
    /// Integrality only (an `INTORG`/`INTEND` marker column).
    Integer,
    /// `LI v`: integer with lower bound `v`.
    IntLower(f64),
    /// `UI v`: integer with upper bound `v`.
    IntUpper(f64),
    /// `SC v`: semi-continuous, upper bound `v` (HMpsFF.cpp:1392).
    Semi(f64),
    /// `SI v`: semi-integer, upper bound `v` (HMpsFF.cpp:1386).
    SemiInt(f64),
}

/// Parses a free-format MPS document.
pub fn read_mps(text: &str) -> Result<Model, MpsError> {
    read_mps_file_order(text).map(|(m, _)| m)
}

/// The COLUMNS-section entries exactly as they appear in the file, before
/// [`CscMatrix::from_triplets`](rusolve_core::CscMatrix::from_triplets)
/// sorts and sums them (and before objective entries are summed).
///
/// A reader that must reproduce another solver's LP bit for bit needs this
/// raw order: HiGHS's free-format reader (MIT, `io/HMpsFF.cpp`
/// `parseColumns`, :890-960 and `fillMatrix`, :136-170) keeps each column's
/// nonzeros in file order, skips zero values and keeps the FIRST of
/// duplicate entries, and the simplex pivot sequence depends on that order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MpsFileOrder {
    /// `(row, col, value)` constraint entries in file order (zeros and
    /// duplicates included).
    pub entries: Vec<Triplet>,
    /// `(col, value)` objective-row entries in file order.
    pub cost_entries: Vec<(usize, f64)>,
}

impl MpsFileOrder {
    /// The order a reader without a COLUMNS section keeps (LP, `.nl`,
    /// MathOptFormat input): the CSC matrix column by column, costs in
    /// column order.
    pub fn of_csc(m: &Model) -> MpsFileOrder {
        let (ptr, idx, val) = (m.matrix.col_ptr(), m.matrix.row_idx(), m.matrix.values());
        let mut entries = Vec::with_capacity(val.len());
        for j in 0..m.num_variables {
            for k in ptr[j]..ptr[j + 1] {
                entries.push((idx[k] as usize, j, val[k]));
            }
        }
        MpsFileOrder {
            entries,
            cost_entries: m.objective.iter().copied().enumerate().collect(),
        }
    }
}

/// [`read_mps_file_order`], but a file with `QCMATRIX` or `CSECTION` is kept
/// (in [`Model::conic`]) instead of rejected. Only the conic route calls it.
pub fn read_mps_conic(text: &str) -> Result<Model, MpsError> {
    read_mps_conic_order(text).map(|(m, _)| m)
}

/// [`read_mps`], also returning the COLUMNS entries in file order. A model
/// with quadratic rows or cones is an error here: the LP/MIP code must never
/// see a model with a constraint silently dropped.
pub fn read_mps_file_order(text: &str) -> Result<(Model, MpsFileOrder), MpsError> {
    let (m, order) = read_mps_conic_order(text)?;
    if !m.conic.is_empty() {
        crate::conic::note_conic_rejection();
        return err(0, crate::conic::CONIC_CONTENT_MSG);
    }
    Ok((m, order))
}

/// [`read_mps_conic`], also returning the COLUMNS entries in file order.
pub fn read_mps_conic_order(text: &str) -> Result<(Model, MpsFileOrder), MpsError> {
    let mut section = Section::None;
    let mut name = String::new();
    let mut obj_sense = ObjSense::default();

    // Rows, in file order. Constraint rows get dense indices; the first `N`
    // row becomes the objective, further `N` rows are ignored.
    let mut senses: Vec<RowSense> = Vec::new();
    let mut row_names: Vec<String> = Vec::new();
    let mut con_index: HashMap<&str, usize> = HashMap::default();
    let mut n_row_names: Vec<&str> = Vec::new();
    let mut obj_name: Option<&str> = None;

    let mut col_names: Vec<String> = Vec::new();
    let mut col_index: HashMap<&str, usize> = HashMap::default();
    let mut triplets: Vec<Triplet> = Vec::new();
    let mut objective: Vec<f64> = Vec::new();
    let mut cost_entries: Vec<(usize, f64)> = Vec::new();

    let mut rhs_entries: Vec<(usize, f64)> = Vec::new();
    let mut range_entries: Vec<(usize, f64)> = Vec::new();
    let mut bound_ops: Vec<(usize, BoundOp)> = Vec::new();
    // RHS entry on the objective (N) row → Model::obj_offset (first
    // entry wins, HiGHS has_obj_entry_ guard, HMpsFF.cpp:1065–1084).
    let mut obj_rhs: Option<f64> = None;
    // The MPS spec keeps only the FIRST RHS/RANGES/BOUNDS set
    // ("only the first RHS vector is used", IBM/CPLEX MPS format
    // reference). Corpus-wide a no-op (no Netlib or MIPLIB instance we test has
    // >1 set name) but spec-correct; e226's single set is ZZZZZZ01.
    let mut rhs_set: Option<&str> = None;
    let mut range_set: Option<&str> = None;
    let mut bound_set: Option<&str> = None;

    let mut in_integer = false; // INTORG/INTEND marker state
    // Marker-integer columns with no BOUNDS entry are binary (MPS convention;
    // HiGHS HMpsFF.cpp:327-331 + :1536, kintegerVarsInColumnsAreBinary, MIT —
    // MIPLIB instances rely on it). Cleared by any BOUNDS line for the column.
    let mut marker_binary: Vec<bool> = Vec::new();

    // HiGHS HMpsFF::q_entries (the objective Hessian) and qrows_entries
    // (one list per constraint row from QSECTION; parsed, not used).
    let mut q_entries: Vec<Triplet> = Vec::new();
    let mut qrows_entries: Vec<Vec<Triplet>> = Vec::new();
    // `INDICATORS` entries (line, row, binary column, active value) and
    // `SOS` sets, resolved after the model is complete.
    let mut ind_entries: Vec<(usize, usize, usize, bool)> = Vec::new();
    let mut sos_raw: Vec<(usize, SosSet)> = Vec::new();
    // QCMATRIX monomials `(row, [(i, j, m)])` per section, and the cones.
    let mut qc_entries: Vec<(usize, Vec<Triplet>)> = Vec::new();
    let mut cones: Vec<Cone> = Vec::new();

    // One token buffer for the whole file (a Vec per line cost an allocation each).
    let mut t: Vec<&str> = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        let line = i + 1;
        t.clear();
        if raw.is_ascii() {
            // Same tokens as `split_whitespace` on ASCII text (Unicode
            // White_Space on ASCII is space and 0x09..=0x0D), without the
            // char decoding.
            let bytes = raw.as_bytes();
            let mut k = 0;
            let n = bytes.len();
            while k < n {
                while k < n && (bytes[k] == b' ' || (9..=13).contains(&bytes[k])) {
                    k += 1;
                }
                let st = k;
                while k < n && !(bytes[k] == b' ' || (9..=13).contains(&bytes[k])) {
                    k += 1;
                }
                if k > st {
                    t.push(&raw[st..k]);
                }
            }
        } else {
            t.extend(raw.split_whitespace());
        }
        // A comment is a `*` in column 1 only (HiGHS HMpsFF.cpp:226, MIT:
        // `strline[0] == '*'` before trimming). An indented `*` starts a name:
        // irish-electricity's RHS set is `*RHS*`, and treating those lines as
        // comments dropped every right-hand side (wrong Optimal 0).
        if t.is_empty() || raw.starts_with('*') {
            continue;
        }

        // Section headers take priority over data lines — EXCEPT inside
        // RHS/RANGES, where `RHS R0 5` (set-name-led entry) is data: a
        // keyword-led line with ≥3 tokens there is never a header.
        let header = match t[0] {
            "NAME" | "NAMES" => Some(Section::Name),
            "OBJSENSE" => Some(Section::ObjSense),
            "ROWS" => Some(Section::Rows),
            "COLUMNS" => Some(Section::Columns),
            "RHS" if !matches!(section, Section::Rhs | Section::Ranges) || t.len() < 3 => {
                Some(Section::Rhs)
            }
            "RANGES" if !(section == Section::Ranges && t.len() >= 3) => Some(Section::Ranges),
            "BOUNDS" => Some(Section::Bounds),
            // HMpsFF.cpp:418-447 (checkFirstWord): QUADOBJ / QMATRIX are a
            // section only when nothing follows them on the line; QSECTION,
            // QCMATRIX and CSECTION take the rest of the line as argument.
            // HiGHS: HMpsFF.cpp::HMpsFF::parseHessian (lines 1738-1746)
            "QUADOBJ" if t.len() == 1 => Some(Section::Quad(QuadSection {
                name: "QUADOBJ",
                triangular: true,
                dest: QuadDest::Objective,
            })),
            "QMATRIX" if t.len() == 1 => Some(Section::Quad(QuadSection {
                name: "QMATRIX",
                triangular: false,
                dest: QuadDest::Objective,
            })),
            // HiGHS: HMpsFF.cpp::HMpsFF::parseQuadRows (lines 1748-1801)
            "QSECTION" => {
                let section_name = "QSECTION";
                // Get row name from section argument
                let rowname = t.get(1).copied().unwrap_or("");
                if rowname.is_empty() {
                    return err(
                        line,
                        format!("No row name given in argument of {section_name}"),
                    );
                }
                let rowidx: Option<QuadDest> = if Some(rowname) == obj_name {
                    Some(QuadDest::Objective)
                } else {
                    con_index.get(rowname).map(|&r| QuadDest::Row(r))
                };
                // if row of section does not exist or is free (index -2),
                // then skip
                match rowidx {
                    None => {
                        // An undefined row (not a free N row) is ignored with
                        // a kWarning in HiGHS (HMpsFF.cpp:1767-1772:
                        // "Row name \"%s\" in %s section is not defined:
                        // ignored"); rusolve-io has no log channel.
                        Some(Section::QuadSkip)
                    }
                    Some(dest) => {
                        if let QuadDest::Row(_) = dest {
                            qrows_entries.resize(senses.len(), Vec::new());
                        }
                        Some(Section::Quad(QuadSection {
                            name: section_name,
                            triangular: true,
                            dest,
                        }))
                    }
                }
            }
            // QCMATRIX / CSECTION: public MOSEK / CPLEX / Gurobi MPS docs
            // (the HiGHS reader skips both). One name on the QCMATRIX line;
            // `CSECTION <name> <param> <type>`.
            "QCMATRIX" => {
                let Some(&rname) = t.get(1) else {
                    return err(line, "QCMATRIX needs a row name");
                };
                match con_index.get(rname) {
                    Some(&r) => Some(Section::QcMatrix(r)),
                    None if Some(rname) == obj_name => {
                        return err(
                            line,
                            "QCMATRIX on the objective row is not supported (use QMATRIX)",
                        );
                    }
                    None => {
                        return err(
                            line,
                            format!("QCMATRIX references undeclared row '{rname}'"),
                        );
                    }
                }
            }
            "CSECTION" => {
                if t.len() < 4 {
                    return err(line, "CSECTION needs `<name> <param> <type>`");
                }
                let kind = match t[3] {
                    "QUAD" => ConeKind::Soc,
                    "RQUAD" => ConeKind::RotatedSoc,
                    "PEXP" | "DEXP" => {
                        return err(
                            line,
                            format!(
                                "cone '{}' has type {}: exponential cones are not implemented",
                                t[1], t[3]
                            ),
                        );
                    }
                    "PPOW" | "DPOW" => {
                        return err(
                            line,
                            format!(
                                "cone '{}' has type {}: power cones are not implemented",
                                t[1], t[3]
                            ),
                        );
                    }
                    other => return err(line, format!("unknown CSECTION type '{other}'")),
                };
                cones.push(Cone {
                    name: t[1].to_string(),
                    kind,
                    vars: Vec::new(),
                });
                Some(Section::Cone)
            }
            // CPLEX / Gurobi MPS
            // extensions for general constraints.
            "INDICATORS" if t.len() == 1 => Some(Section::Indicators),
            "SOS" if t.len() == 1 => Some(Section::Sos),
            "ENDATA" => break,
            _ => None,
        };
        if let Some(h) = header {
            section = h;
            if let Section::QcMatrix(r) = h {
                qc_entries.push((r, Vec::new()));
            }
            if h == Section::Name
                && let Some(n) = t.get(1)
            {
                name = (*n).to_string();
                section = Section::None;
            }
            // Same-line sense, as Gurobi (`OBJSENSE MAX`) and Xpress (`OBJSENSE  MAXIMIZE`)
            // write it: HiGHS HMpsFF.cpp:528-545 (MIT) takes `MAX`/`MIN` there and keeps the
            // section open for a later "proper" line; we also take the long spellings
            // (`MAXIMIZE`, `MINIMIZE`). Ignoring the same-line sense would silently
            // solve the minimisation (corpus: gurobi.mps optimum 1 instead of 25).
            if h == Section::ObjSense
                && let Some(s) = t.get(1)
            {
                match s.to_ascii_uppercase().as_str() {
                    "MAX" | "MAXIMIZE" => obj_sense = ObjSense::Maximize,
                    "MIN" | "MINIMIZE" => obj_sense = ObjSense::Minimize,
                    other => return err(line, format!("bad OBJSENSE '{other}'")),
                }
            }
            continue;
        }

        match section {
            Section::None => {
                return err(
                    line,
                    format!("unknown section or stray data '{t0}'", t0 = t[0]),
                );
            }
            Section::Name => {
                name = t[0].to_string();
                section = Section::None;
            }
            Section::ObjSense => {
                obj_sense = match t[0].to_ascii_uppercase().as_str() {
                    "MIN" | "MINIMIZE" => ObjSense::Minimize,
                    "MAX" | "MAXIMIZE" => ObjSense::Maximize,
                    other => return err(line, format!("bad OBJSENSE '{other}'")),
                };
            }
            Section::Rows => {
                if t.len() < 2 {
                    return err(line, "ROWS entry needs a sense and a name");
                }
                let (sense, rname) = (t[0], t[1]);
                if con_index.contains_key(rname) || n_row_names.contains(&rname) {
                    return err(line, format!("duplicate row '{rname}'"));
                }
                match sense {
                    "N" => {
                        if obj_name.is_none() {
                            obj_name = Some(rname);
                        }
                        n_row_names.push(rname);
                    }
                    "L" | "G" | "E" => {
                        let s = match sense {
                            "L" => RowSense::LessEqual,
                            "G" => RowSense::GreaterEqual,
                            _ => RowSense::Equal,
                        };
                        con_index.insert(rname, senses.len());
                        senses.push(s);
                    }
                    other => return err(line, format!("bad row sense '{other}'")),
                }
                // Only constraint rows go into row_names: the objective row is
                // carried by obj_name, extra N rows are dropped.
                if sense != "N" {
                    row_names.push(rname.to_string());
                }
            }
            Section::Columns => {
                if t.len() >= 3 && t[1].trim_matches('\'') == "MARKER" {
                    in_integer = match t[2].trim_matches('\'') {
                        "INTORG" => true,
                        "INTEND" => false,
                        other => return err(line, format!("bad marker '{other}'")),
                    };
                    continue;
                }
                if t.len() < 3 || t.len() % 2 != 1 {
                    return err(line, "COLUMNS entry needs a column and (row, value) pairs");
                }
                let cname = t[0];
                let j = match col_index.get(cname) {
                    Some(&j) => j,
                    None => {
                        col_index.insert(cname, col_names.len());
                        col_names.push(cname.to_string());
                        objective.push(0.0);
                        marker_binary.push(in_integer);
                        if in_integer {
                            // Recorded as an Integer op so bound application stays uniform.
                            bound_ops.push((col_names.len() - 1, BoundOp::Integer));
                        }
                        col_names.len() - 1
                    }
                };
                let mut k = 1;
                while k < t.len() {
                    let rname = t[k];
                    let v = parse_f64(t[k + 1], line, "column")?;
                    if Some(rname) == obj_name {
                        objective[j] += v;
                        cost_entries.push((j, v));
                    } else if let Some(&r) = con_index.get(rname) {
                        triplets.push((r, j, v));
                    } else if !n_row_names.contains(&rname) {
                        return err(
                            line,
                            format!("column '{cname}' references undeclared row '{rname}'"),
                        );
                    }
                    k += 2;
                }
            }
            Section::Rhs | Section::Ranges => {
                // Optional trailing set-name field: with it the token count is
                // odd (name + pairs), without it even.
                let has_set = t.len() % 2 == 1;
                // First-set-only rule: a line led by a set name other
                // than the section's first is dropped; unnamed lines join the
                // current set (tolerant, as elsewhere in this reader).
                if has_set {
                    let seen = if section == Section::Rhs {
                        &mut rhs_set
                    } else {
                        &mut range_set
                    };
                    if seen.is_none() {
                        *seen = Some(t[0]);
                    } else if *seen != Some(t[0]) {
                        continue;
                    }
                }
                let start = if has_set { 1 } else { 0 };
                let mut k = start;
                while k + 1 < t.len() {
                    let rname = t[k];
                    let v = parse_f64(t[k + 1], line, "rhs/ranges")?;
                    if section == Section::Rhs && Some(rname) == obj_name {
                        // Objective-row RHS = objective constant, first
                        // entry wins; HiGHS maps it to obj_offset = -val
                        // (HMpsFF.cpp:1065–1084). Further N rows stay ignored.
                        if obj_rhs.is_none() {
                            obj_rhs = Some(v);
                        }
                    } else if let Some(&r) = con_index.get(rname) {
                        if section == Section::Rhs {
                            rhs_entries.push((r, v));
                        } else {
                            range_entries.push((r, v));
                        }
                    } else if !n_row_names.contains(&rname) {
                        let what = if section == Section::Rhs {
                            "RHS"
                        } else {
                            "RANGES"
                        };
                        return err(line, format!("{what} references undeclared row '{rname}'"));
                    }
                    k += 2;
                }
            }
            Section::Bounds => {
                const TYPES: [&str; 11] = [
                    "UP", "LO", "FX", "FR", "MI", "PL", "BV", "LI", "UI", "SC", "SI",
                ];
                // Optional bound-set-name field before the type keyword.
                let ty_idx = if TYPES.contains(&t[0]) { 0 } else { 1 };
                // First-set-only rule: the set name is t[0] when it
                // precedes the type keyword, else the standard `TYPE <set>
                // <col>` slot — whichever token after the type is not a
                // declared column (same tolerance as the column lookup).
                let set_name: Option<&str> = if ty_idx == 1 {
                    Some(t[0])
                } else if !t.get(1).is_some_and(|c| col_index.contains_key(*c)) {
                    t.get(1).copied()
                } else {
                    None
                };
                if let Some(s) = set_name {
                    if bound_set.is_none() {
                        bound_set = Some(s);
                    } else if bound_set != Some(s) {
                        continue;
                    }
                }
                let ty = match t.get(ty_idx) {
                    Some(&s) if TYPES.contains(&s) => s,
                    Some(s) => return err(line, format!("bad bound type '{s}'")),
                    None => return err(line, "BOUNDS entry needs a type"),
                };
                // Accept both `UP <col> [val]` and standard `UP <set> <col> [val]`:
                // the column is whichever of the next tokens is declared.
                let (cname, val_tok) = match (t.get(ty_idx + 1), t.get(ty_idx + 2)) {
                    (Some(a), _) if col_index.contains_key(*a) => (*a, t.get(ty_idx + 2)),
                    (Some(_), Some(b)) if col_index.contains_key(*b) => (*b, t.get(ty_idx + 3)),
                    (Some(a), _) => {
                        return err(line, format!("bounds reference undeclared column '{a}'"));
                    }
                    (None, _) => return err(line, "BOUNDS entry needs a column name"),
                };
                let &j = col_index.get(cname).ok_or_else(|| MpsError {
                    line,
                    message: format!("bounds reference undeclared column '{cname}'"),
                })?;
                let val = match val_tok {
                    Some(s) => Some(parse_f64(s, line, "bound")?),
                    None => None,
                };
                let op = match ty {
                    "UP" => BoundOp::Upper(need_val(val, line, "UP")?),
                    "LO" => BoundOp::Lower(need_val(val, line, "LO")?),
                    "FX" => BoundOp::Fixed(need_val(val, line, "FX")?),
                    "FR" => BoundOp::Free,
                    "MI" => BoundOp::MinusInf,
                    "PL" => BoundOp::PlusInf,
                    "BV" => BoundOp::Binary,
                    "LI" => BoundOp::IntLower(need_val(val, line, "LI")?),
                    "UI" => BoundOp::IntUpper(need_val(val, line, "UI")?),
                    "SC" => BoundOp::Semi(need_val(val, line, "SC")?),
                    "SI" => BoundOp::SemiInt(need_val(val, line, "SI")?),
                    _ => unreachable!("checked above"),
                };
                marker_binary[j] = false;
                bound_ops.push((j, op));
            }
            // HiGHS: HMpsFF.cpp::HMpsFF::parseQuadMatrix (lines 1803-1889)
            //
            // QMATRIX should have all entries, whereas every off-diagonal
            // QUADOBJ/QSECTION entry also defines its entry in the opposite
            // triangle, so add this explicitly to unify the record
            // regardless of section.
            Section::Quad(q) => {
                let section_name = q.name;
                let triangular = q.triangular;
                let col_name = t[0];
                // Get the column index from the name
                let colidx = get_col_idx(
                    col_name,
                    &mut col_index,
                    &mut col_names,
                    &mut objective,
                    &mut marker_binary,
                );
                // Loop over the maximum of two entries per row of the file
                let mut end = 1;
                for _entry in 0..2 {
                    // Get the row name
                    let Some(&row_name) = t.get(end) else {
                        break;
                    };
                    let Some(&coeff_name) = t.get(end + 1) else {
                        return err(
                            line,
                            format!(
                                "{section_name} has no coefficient for entry \"{row_name}\" in column \"{col_name}\""
                            ),
                        );
                    };
                    let rowidx = get_col_idx(
                        row_name,
                        &mut col_index,
                        &mut col_names,
                        &mut objective,
                        &mut marker_binary,
                    );
                    // HiGHS getValue is atof after a D -> E swap and never
                    // sets is_nan (HMpsFF.cpp:2088-2112), so the NaN error
                    // (:1859-1865) is unreachable; the coefficient is read
                    // with this reader's parse_f64.
                    let coeff = parse_f64(coeff_name, line, "Hessian coefficient")?;
                    if coeff != 0.0 {
                        let entries = match q.dest {
                            QuadDest::Objective => &mut q_entries,
                            QuadDest::Row(r) => &mut qrows_entries[r],
                        };
                        entries.push((rowidx, colidx, coeff));
                        // For a triangular section (QUADOBJ/QSECTION), make a
                        // separate record of the entry in the opposite
                        // triangle
                        if triangular && rowidx != colidx {
                            entries.push((colidx, rowidx, coeff));
                        }
                    }
                    end += 2;
                    // Don't read more if end of line reached
                    if end == t.len() {
                        break;
                    }
                }
            }
            // HMpsFF.cpp:1776-1788: lines of a skipped QSECTION are read
            // and dropped.
            Section::QuadSkip => {}
            Section::Indicators => {
                // `IF row binvar value` (the leading `IF` is the only
                // indicator type the public CPLEX/Gurobi formats define).
                if t.len() != 4 || t[0] != "IF" {
                    return err(line, "INDICATORS entry must be `IF row binvar value`");
                }
                let Some(&r) = con_index.get(t[1]) else {
                    return err(
                        line,
                        format!("INDICATORS references undeclared constraint row '{}'", t[1]),
                    );
                };
                let Some(&j) = col_index.get(t[2]) else {
                    return err(
                        line,
                        format!("INDICATORS references undeclared column '{}'", t[2]),
                    );
                };
                let v = parse_f64(t[3], line, "indicator value")?;
                if v != 0.0 && v != 1.0 {
                    return err(line, format!("indicator value must be 0 or 1, got {v}"));
                }
                ind_entries.push((line, r, j, v == 1.0));
            }
            Section::Sos => {
                if t[0] == "S1" || t[0] == "S2" {
                    let kind = if t[0] == "S1" {
                        SosKind::Sos1
                    } else {
                        SosKind::Sos2
                    };
                    // ` S1 SOS name [priority]` or ` S1 name [priority]`
                    let rest = if t.get(1) == Some(&"SOS") && t.len() >= 3 {
                        &t[2..]
                    } else {
                        &t[1..]
                    };
                    let Some(set_name) = rest.first() else {
                        return err(line, "SOS set has no name");
                    };
                    let priority = match rest.get(1) {
                        Some(s) => parse_f64(s, line, "SOS priority")?,
                        None => 0.0,
                    };
                    sos_raw.push((
                        line,
                        SosSet {
                            name: (*set_name).to_string(),
                            kind,
                            priority,
                            members: Vec::new(),
                        },
                    ));
                    continue;
                }
                let Some((_, cur)) = sos_raw.last_mut() else {
                    return err(line, "SOS member line before any `S1`/`S2` set line");
                };
                // `col [weight]` (HiGHS kSos) or `set col [weight]`
                let (cname, wtok) =
                    if t.len() >= 2 && t[0] == cur.name && col_index.contains_key(t[1]) {
                        (t[1], t.get(2))
                    } else {
                        (t[0], t.get(1))
                    };
                let Some(&j) = col_index.get(cname) else {
                    return err(line, format!("SOS references undeclared column '{cname}'"));
                };
                let w = match wtok {
                    Some(s) => parse_f64(s, line, "SOS weight")?,
                    None => cur.members.len() as f64 + 1.0,
                };
                cur.members.push((j, w));
            }
            // `<col> <col> <value>` per line, one monomial `value * x_i x_j`.
            Section::QcMatrix(_) => {
                if t.len() != 3 {
                    return err(line, "QCMATRIX entry needs `<col> <col> <value>`");
                }
                let mut ix = [0usize; 2];
                for (k, slot) in ix.iter_mut().enumerate() {
                    *slot = *col_index.get(t[k]).ok_or_else(|| MpsError {
                        line,
                        message: format!("QCMATRIX references undeclared column '{}'", t[k]),
                    })?;
                }
                let v = parse_f64(t[2], line, "QCMATRIX")?;
                if let Some((_, ent)) = qc_entries.last_mut() {
                    ent.push((ix[0], ix[1], v));
                }
            }
            // One column name per line.
            Section::Cone => {
                if t.len() != 1 {
                    return err(line, "a CSECTION line is one column name");
                }
                let &j = col_index.get(t[0]).ok_or_else(|| MpsError {
                    line,
                    message: format!("CSECTION references undeclared column '{}'", t[0]),
                })?;
                if let Some(c) = cones.last_mut() {
                    c.vars.push(j);
                }
            }
        }
    }
    // qrows_entries (QSECTION on a constraint row) are parsed but not used.
    drop(qrows_entries);

    let num_variables = col_names.len();
    let num_constraints = senses.len();
    let mut objective_final = objective;
    objective_final.resize(num_variables, 0.0);

    let mut rhs = vec![0.0; num_constraints];
    for (r, v) in rhs_entries {
        rhs[r] = v;
    }
    let mut ranges = vec![0.0; num_constraints];
    for (r, v) in range_entries {
        ranges[r] = v;
    }

    let mut bounds = vec![VariableBound::default(); num_variables];
    // SC / SI columns (last directive on a column wins, as in HiGHS,
    // where the integrality of the column is overwritten).
    let mut semi_kind: Vec<Option<bool>> = vec![None; num_variables];
    for (j, op) in bound_ops {
        let b = &mut bounds[j];
        match op {
            BoundOp::Semi(v) | BoundOp::SemiInt(v) => {
                b.upper = v;
                b.integer = false;
                semi_kind[j] = Some(matches!(op, BoundOp::SemiInt(_)));
                continue;
            }
            BoundOp::Binary | BoundOp::Integer | BoundOp::IntLower(_) | BoundOp::IntUpper(_) => {
                semi_kind[j] = None;
            }
            _ => {}
        }
        match op {
            BoundOp::Semi(_) | BoundOp::SemiInt(_) => unreachable!("handled above"),
            BoundOp::Upper(v) => b.upper = v,
            BoundOp::Lower(v) => b.lower = v,
            BoundOp::Fixed(v) => {
                b.lower = v;
                b.upper = v;
            }
            BoundOp::Free => {
                b.lower = f64::NEG_INFINITY;
                b.upper = f64::INFINITY;
            }
            BoundOp::MinusInf => b.lower = f64::NEG_INFINITY,
            BoundOp::PlusInf => b.upper = f64::INFINITY,
            BoundOp::Binary => {
                b.lower = 0.0;
                b.upper = 1.0;
                b.integer = true;
            }
            BoundOp::Integer => b.integer = true,
            BoundOp::IntLower(v) => {
                b.integer = true;
                b.lower = v;
            }
            BoundOp::IntUpper(v) => {
                b.integer = true;
                b.upper = v;
            }
        }
    }

    for (b, _) in bounds
        .iter_mut()
        .zip(&marker_binary)
        .filter(|(_, bin)| **bin)
    {
        b.upper = 1.0;
    }

    let mut gc = GcModel::default();
    for (j, k) in semi_kind.iter().enumerate() {
        if let Some(integer) = *k {
            gc.semi.push(SemiVar { col: j, integer });
        }
    }
    // An indicator row is declared in ROWS/COLUMNS/RHS like any row and
    // then taken out of the constraint set; it only holds when its binary
    // takes the active value.
    // Indicator rows leave the row list, so the QCMATRIX row numbers are remapped
    // below; a QCMATRIX on an indicator row itself is refused (an indicator holds
    // a linear row only).
    let (senses, row_names, rhs, ranges, triplets, num_constraints) = if ind_entries.is_empty() {
        (senses, row_names, rhs, ranges, triplets, num_constraints)
    } else {
        let mut is_ind = vec![false; num_constraints];
        for &(line, r, _, _) in &ind_entries {
            if is_ind[r] {
                return err(
                    line,
                    format!("row '{}' is named by two INDICATORS entries", row_names[r]),
                );
            }
            if ranges[r] != 0.0 {
                return err(
                    line,
                    format!("indicator row '{}' has a RANGES entry", row_names[r]),
                );
            }
            is_ind[r] = true;
        }
        let mut new_index = vec![usize::MAX; num_constraints];
        let mut kept = 0usize;
        for r in 0..num_constraints {
            if !is_ind[r] {
                new_index[r] = kept;
                kept += 1;
            }
        }
        for (r, _) in qc_entries.iter_mut() {
            if is_ind[*r] {
                let line = ind_entries.iter().find(|e| e.1 == *r).map_or(0, |e| e.0);
                return err(
                    line,
                    format!(
                        "indicator row '{}' has a QCMATRIX section; an indicator row must be linear",
                        row_names[*r]
                    ),
                );
            }
            *r = new_index[*r];
        }
        let mut terms_of: Vec<Vec<(usize, f64)>> = vec![Vec::new(); num_constraints];
        let mut kept_triplets: Vec<Triplet> = Vec::with_capacity(triplets.len());
        for &(r, j, v) in &triplets {
            if is_ind[r] {
                if v != 0.0 {
                    terms_of[r].push((j, v));
                }
            } else {
                kept_triplets.push((new_index[r], j, v));
            }
        }
        for &(_, r, j, bin_val) in &ind_entries {
            gc.indicators.push(Indicator {
                name: row_names[r].clone(),
                bin: j,
                bin_val,
                terms: std::mem::take(&mut terms_of[r]),
                sense: senses[r],
                rhs: rhs[r],
            });
        }
        fn pick<T: Clone>(v: &[T], is_ind: &[bool]) -> Vec<T> {
            v.iter()
                .zip(is_ind)
                .filter(|(_, i)| !**i)
                .map(|(x, _)| x.clone())
                .collect()
        }
        (
            pick(&senses, &is_ind),
            pick(&row_names, &is_ind),
            pick(&rhs, &is_ind),
            pick(&ranges, &is_ind),
            kept_triplets,
            kept,
        )
    };
    for (line, s) in sos_raw {
        if s.members.is_empty() {
            return err(line, format!("SOS set '{}' has no members", s.name));
        }
        gc.sos.push(s);
    }
    gc.validate(&bounds, &col_names)
        .map_err(|message| MpsError { line: 0, message })?;

    let matrix = rusolve_core::CscMatrix::from_triplets(num_constraints, num_variables, &triplets);

    // Sections for one row merge, rows ascend.
    let mut by_row: std::collections::BTreeMap<usize, Vec<Triplet>> = Default::default();
    for (r, ent) in qc_entries {
        by_row.entry(r).or_default().extend(ent);
    }
    let quad_rows: Vec<QuadRow> = by_row
        .into_iter()
        .map(|(r, ent)| QuadRow::from_monomials(r, ent))
        .filter(|q| !q.q.is_empty())
        .collect();
    for c in &cones {
        let need = match c.kind {
            ConeKind::RotatedSoc => 3,
            _ => 2,
        };
        if c.vars.len() < need {
            return err(
                0,
                format!(
                    "cone '{}' ({}) needs at least {need} columns, has {}",
                    c.name,
                    c.kind.label(),
                    c.vars.len()
                ),
            );
        }
    }

    let file_order = MpsFileOrder {
        entries: triplets,
        cost_entries,
    };
    let model = Model {
        name,
        num_variables,
        num_constraints,
        obj_sense,
        obj_name: obj_name.map(str::to_string),
        // HiGHS semantics: the constant added to the reported objective is
        // the NEGATED N-row RHS entry (HMpsFF.cpp:1081 `obj_offset = -val`).
        obj_offset: obj_rhs.map_or(0.0, |v| -v),
        objective: objective_final,
        senses,
        rhs,
        ranges,
        bounds,
        matrix,
        row_names,
        col_names,
        hessian: q_entries,
        gc,
        conic: ConicBlock { quad_rows, cones },
        cp: Default::default(),
    };
    Ok((model, file_order))
}

fn need_val(v: Option<f64>, line: usize, ty: &str) -> Result<f64, MpsError> {
    v.ok_or_else(|| MpsError {
        line,
        message: format!("bound type {ty} requires a value"),
    })
}

fn fmt_val(v: f64) -> String {
    // Rust's Display for f64 is shortest-round-trip: read-back is exact.
    format!("{v}")
}

fn sense_char(s: RowSense) -> char {
    match s {
        RowSense::LessEqual => 'L',
        RowSense::GreaterEqual => 'G',
        RowSense::Equal => 'E',
    }
}

/// BOUND directives reproducing `b` (relative to the MPS default `0 ≤ x < ∞`).
/// Integrality is carried by the `INTORG`/`INTEND` markers [`write_mps`]
/// puts around integer columns; a marker column with no BOUNDS line reads
/// back as binary, so an integer column always gets at least one line here
/// (`BV` for `[0, 1]`, `PL` for the default box `[0, +∞)`).
fn bound_lines(col: &str, b: &VariableBound) -> Vec<String> {
    let mut out = Vec::new();
    if b.integer && b.lower == 0.0 && b.upper == 1.0 {
        out.push(format!("BV bnd {col}"));
        return out;
    }
    let fin_lo = b.lower.is_finite();
    let fin_up = b.upper.is_finite();
    // Standard field order: TYPE <set-name> <col> [value]
    if !fin_lo && !fin_up {
        out.push(format!("FR bnd {col}"));
    } else if fin_lo && fin_up && b.lower == b.upper {
        out.push(format!("FX bnd {col} {}", fmt_val(b.lower)));
    } else {
        if !fin_lo {
            out.push(format!("MI bnd {col}"));
        }
        if fin_lo && b.lower != 0.0 {
            out.push(format!("LO bnd {col} {}", fmt_val(b.lower)));
        }
        if fin_up {
            out.push(format!("UP bnd {col} {}", fmt_val(b.upper)));
        }
    }
    if b.integer && out.is_empty() {
        out.push(format!("PL bnd {col}"));
    }
    out
}

/// Writes `m` as a free-format MPS document. `read_mps` of the output
/// reproduces `m` semantically (matrix/bounds/rhs/senses equal).
pub fn write_mps(m: &Model) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "NAME {}", m.name);
    let _ = writeln!(
        s,
        "OBJSENSE\n {}",
        match m.obj_sense {
            ObjSense::Minimize => "MIN",
            ObjSense::Maximize => "MAX",
        }
    );

    let obj = m.obj_name.as_deref().unwrap_or("OBJ");
    let _ = writeln!(s, "ROWS");
    let _ = writeln!(s, " N  {obj}");
    for (r, sense) in m.senses.iter().enumerate() {
        let _ = writeln!(s, " {}  {}", sense_char(*sense), m.row_names[r]);
    }
    // Indicator rows are declared like rows and listed in INDICATORS.
    for o in &m.gc.indicators {
        let _ = writeln!(s, " {}  {}", sense_char(o.sense), o.name);
    }
    let mut ind_terms: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m.num_variables];
    for (k, o) in m.gc.indicators.iter().enumerate() {
        for &(j, a) in &o.terms {
            ind_terms[j].push((k, a));
        }
    }

    let _ = writeln!(s, "COLUMNS");
    let mut in_marker = false;
    #[allow(clippy::needless_range_loop)] // `j` also indexes the model's columns
    for j in 0..m.num_variables {
        let col = &m.col_names[j];
        if m.bounds[j].integer != in_marker {
            in_marker = m.bounds[j].integer;
            let tag = if in_marker { "INTORG" } else { "INTEND" };
            let _ = writeln!(s, "  MARKER 'MARKER' '{tag}'");
        }
        let empty = m.matrix.col_ptr()[j] == m.matrix.col_ptr()[j + 1];
        if m.objective[j] != 0.0 || empty {
            // An entry-less column still needs one line to be declared.
            let _ = writeln!(s, "  {col} {obj} {}", fmt_val(m.objective[j]));
        }
        for k in m.matrix.col_ptr()[j]..m.matrix.col_ptr()[j + 1] {
            let r = m.matrix.row_idx()[k] as usize;
            let _ = writeln!(
                s,
                "  {col} {} {}",
                m.row_names[r],
                fmt_val(m.matrix.values()[k])
            );
        }
        for &(k, a) in &ind_terms[j] {
            let _ = writeln!(s, "  {col} {} {}", m.gc.indicators[k].name, fmt_val(a));
        }
    }
    if in_marker {
        let _ = writeln!(s, "  MARKER 'MARKER' 'INTEND'");
    }

    // The objective constant round-trips as an N-row RHS entry: the
    // reader maps an N-row value v to obj_offset = -v (HiGHS semantics,
    // HMpsFF.cpp:1065–1084), so emit the negation under the set name the
    // constraint entries already use.
    if m.rhs.iter().any(|&v| v != 0.0)
        || m.obj_offset != 0.0
        || m.gc.indicators.iter().any(|o| o.rhs != 0.0)
    {
        let _ = writeln!(s, "RHS");
        for o in &m.gc.indicators {
            if o.rhs != 0.0 {
                let _ = writeln!(s, "  rhs {} {}", o.name, fmt_val(o.rhs));
            }
        }
        if m.obj_offset != 0.0 {
            let _ = writeln!(s, "  rhs {obj} {}", fmt_val(-m.obj_offset));
        }
        for (r, &v) in m.rhs.iter().enumerate() {
            if v != 0.0 {
                let _ = writeln!(s, "  rhs {} {}", m.row_names[r], fmt_val(v));
            }
        }
    }
    if m.ranges.iter().any(|&v| v != 0.0) {
        let _ = writeln!(s, "RANGES");
        for (r, &v) in m.ranges.iter().enumerate() {
            if v != 0.0 {
                let _ = writeln!(s, "  rng {} {}", m.row_names[r], fmt_val(v));
            }
        }
    }

    let mut bounds_written = false;
    for j in 0..m.num_variables {
        let lines = match m.gc.semi.iter().find(|sv| sv.col == j) {
            // `LO` when the lower bound is not the default, then `SC`/`SI`
            // carrying the upper bound (HMpsFF.cpp:1386-1396).
            Some(sv) => {
                let b = &m.bounds[j];
                let mut l = Vec::new();
                if b.lower == f64::NEG_INFINITY {
                    l.push(format!("MI bnd {}", m.col_names[j]));
                } else if b.lower != 0.0 {
                    l.push(format!("LO bnd {} {}", m.col_names[j], fmt_val(b.lower)));
                }
                let ty = if sv.integer { "SI" } else { "SC" };
                l.push(format!("{ty} bnd {} {}", m.col_names[j], fmt_val(b.upper)));
                l
            }
            None => bound_lines(&m.col_names[j], &m.bounds[j]),
        };
        if !lines.is_empty() {
            if !bounds_written {
                let _ = writeln!(s, "BOUNDS");
                bounds_written = true;
            }
            for l in lines {
                let _ = writeln!(s, " {l}");
            }
        }
    }

    if !m.gc.indicators.is_empty() {
        let _ = writeln!(s, "INDICATORS");
        for o in &m.gc.indicators {
            let _ = writeln!(
                s,
                " IF {} {} {}",
                o.name,
                m.col_names[o.bin],
                u8::from(o.bin_val)
            );
        }
    }
    if !m.gc.sos.is_empty() {
        let _ = writeln!(s, "SOS");
        for set in &m.gc.sos {
            let t = if set.kind == SosKind::Sos1 {
                "S1"
            } else {
                "S2"
            };
            let _ = writeln!(s, " {t} SOS {} {}", set.name, fmt_val(set.priority));
            for &(j, w) in &set.members {
                let _ = writeln!(s, " {} {} {}", set.name, m.col_names[j], fmt_val(w));
            }
        }
    }
    if !m.gc.pwl.is_empty() || !m.gc.funcs.is_empty() {
        // The MPS format has no section for these (CPLEX/Gurobi write them
        // to LP files only); say so rather than drop them silently.
        let _ = writeln!(
            s,
            "* NOTE: piecewise-linear and abs/min/max/and/or objects have no MPS section; use the LP format"
        );
    }
    // Quadratic rows as a symmetric QCMATRIX (monomial coefficients, an
    // off-diagonal monomial split over both triangles), then the cones.
    for qr in &m.conic.quad_rows {
        let _ = writeln!(s, "QCMATRIX {}", m.row_names[qr.row]);
        for (i, j, mono) in qr.monomials() {
            let (ci, cj) = (&m.col_names[i], &m.col_names[j]);
            if i == j {
                let _ = writeln!(s, "  {ci} {cj} {}", fmt_val(mono));
            } else {
                let _ = writeln!(s, "  {ci} {cj} {}", fmt_val(mono / 2.0));
                let _ = writeln!(s, "  {cj} {ci} {}", fmt_val(mono / 2.0));
            }
        }
    }
    for c in &m.conic.cones {
        let ty = match c.kind {
            ConeKind::Soc => "QUAD",
            ConeKind::RotatedSoc => "RQUAD",
            ConeKind::Exp => "PEXP",
            ConeKind::Pow(_) => "PPOW",
        };
        let _ = writeln!(s, "CSECTION {} 0 {ty}", c.name);
        for &j in &c.vars {
            let _ = writeln!(s, "  {}", m.col_names[j]);
        }
    }

    let _ = writeln!(s, "ENDATA");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn objsense_on_the_same_line() {
        for (head, want) in [
            ("OBJSENSE MAX", ObjSense::Maximize),
            ("OBJSENSE  MAXIMIZE", ObjSense::Maximize),
            ("OBJSENSE max", ObjSense::Maximize),
            ("OBJSENSE MIN", ObjSense::Minimize),
            ("OBJSENSE\n    MAX", ObjSense::Maximize),
        ] {
            let text = format!(
                "NAME t\n{head}\nROWS\n N obj\n L c\nCOLUMNS\n x obj 1 c 1\nRHS\n rhs c 4\nENDATA\n"
            );
            let m = read_mps(&text).unwrap();
            assert_eq!(m.obj_sense, want, "{head}");
        }
        assert!(read_mps("NAME t\nOBJSENSE SIDEWAYS\nROWS\n N obj\nENDATA\n").is_err());
    }

    const TOL: f64 = 1e-9;

    /// Mixed L/G/E rows with objective, RHS and a dangling-free column set.
    const SMALL_LP: &str = "\
* tiny diet-ish LP
NAME     SMALL
OBJSENSE
  MAX
ROWS
  N  COST
  L  FAT
  G  PROT
  E  EXACT
COLUMNS
  X1  COST  -2.0  FAT  1.0
  X1  PROT  1.0   EXACT 1.0
  X2  COST  -3.0  FAT  2.0
  X2  PROT  1.0
RHS
  rhs  FAT  10.0  PROT  4.0  EXACT  2.0
BOUNDS
  UP X1 6.0
  LO X2 0.5
ENDATA
";

    #[test]
    fn parses_small_lp_with_mixed_row_senses() {
        let m = read_mps(SMALL_LP).unwrap();
        assert_eq!(m.name, "SMALL");
        assert_eq!(m.obj_sense, ObjSense::Maximize);
        assert_eq!(m.obj_name.as_deref(), Some("COST"));
        assert_eq!(m.num_variables, 2);
        assert_eq!(m.num_constraints, 3);
        assert_eq!(
            m.senses,
            [RowSense::LessEqual, RowSense::GreaterEqual, RowSense::Equal]
        );
        assert_eq!(m.rhs, [10.0, 4.0, 2.0]);
        assert_eq!(m.objective, [-2.0, -3.0]);
        // Constraint matrix (rows FAT/PROT/EXACT x cols HiGHS backport/X2), N row excluded.
        assert_eq!(m.matrix.get(0, 0), 1.0);
        assert_eq!(m.matrix.get(0, 1), 2.0);
        assert_eq!(m.matrix.get(1, 1), 1.0);
        assert_eq!(m.matrix.get(2, 1), 0.0);
        assert_eq!(m.bounds[0].upper, 6.0);
        assert_eq!(m.bounds[1].lower, 0.5);
        assert_eq!(m.bounds[0].lower, 0.0); // MPS default
    }

    #[test]
    fn parses_all_bound_types() {
        let text = "\
NAME BNDTEST
ROWS
  N  OF
  L  C1
COLUMNS
  A  OF  1  C1  1
  B  OF  1  C1  1
  C  OF  1  C1  1
  D  OF  1  C1  1
  E  OF  1  C1  1
  F  OF  1  C1  1
RHS
  rhs C1 5
BOUNDS
  UP A 4
  LO B 2
  FR C
  BV D
  LI E 3
  FX F 7
ENDATA
";
        let m = read_mps(text).unwrap();
        let b = |n: &str| m.bounds[m.col_names.iter().position(|c| c == n).unwrap()];
        assert_eq!(
            b("A"),
            VariableBound {
                lower: 0.0,
                upper: 4.0,
                integer: false,
            }
        );
        assert_eq!(b("B").lower, 2.0);
        assert_eq!(b("C").lower, f64::NEG_INFINITY);
        assert_eq!(b("C").upper, f64::INFINITY);
        assert_eq!(
            b("D"),
            VariableBound {
                lower: 0.0,
                upper: 1.0,
                integer: true,
            }
        );
        // LI is an integer LOWER bound,
        // not an initial value (HMpsFF.cpp:1371–1375, MIT).
        assert_eq!(
            b("E"),
            VariableBound {
                lower: 3.0,
                upper: f64::INFINITY,
                integer: true,
            }
        );
        assert_eq!(b("F").lower, 7.0);
        assert_eq!(b("F").upper, 7.0);
    }

    /// One column per bound type, in both field orders (`TYPE set col v`
    /// and `TYPE col v`), each on its own default box. Semantics are
    /// HiGHS's `parseBounds` (HMpsFF.cpp:1340–1537, MIT).
    #[test]
    fn every_bound_type_matches_highs_semantics() {
        let text = "\
NAME ALLBND
ROWS
  N  OF
  L  C1
COLUMNS
  UPc  OF  1  C1  1
  LOc  OF  1  C1  1
  FXc  OF  1  C1  1
  FRc  OF  1  C1  1
  MIc  OF  1  C1  1
  PLc  OF  1  C1  1
  BVc  OF  1  C1  1
  LIc  OF  1  C1  1
  UIc  OF  1  C1  1
  UIn  OF  1  C1  1
  UPn  OF  1  C1  1
  LIUI OF  1  C1  1
  MIUP OF  1  C1  1
RHS
  rhs C1 5
BOUNDS
 UP bnd UPc 4
 LO bnd LOc -2
 FX bnd FXc 7
 FR bnd FRc
 MI bnd MIc
 UP bnd PLc 3
 PL bnd PLc
 BV bnd BVc
 LI bnd LIc 3
 UI bnd UIc 9
 UI bnd UIn -4
 UP bnd UPn -4
 LI bnd LIUI -5
 UI LIUI 5
 MI MIUP
 UP MIUP 6
ENDATA
";
        let m = read_mps(text).unwrap();
        let b = |n: &str| {
            let v = m.bounds[m.col_names.iter().position(|c| c == n).unwrap()];
            (v.lower, v.upper, v.integer)
        };
        let inf = f64::INFINITY;
        assert_eq!(b("UPc"), (0.0, 4.0, false));
        assert_eq!(b("LOc"), (-2.0, inf, false));
        assert_eq!(b("FXc"), (7.0, 7.0, false));
        assert_eq!(b("FRc"), (-inf, inf, false));
        assert_eq!(b("MIc"), (-inf, inf, false), "MI leaves the upper bound");
        assert_eq!(b("PLc"), (0.0, inf, false), "PL resets the upper bound");
        assert_eq!(b("BVc"), (0.0, 1.0, true));
        assert_eq!(b("LIc"), (3.0, inf, true));
        assert_eq!(b("UIc"), (0.0, 9.0, true));
        // HiGHS applies no "negative upper bound ⇒ lower −∞" rule, for UI
        // nor UP: the default lower bound 0 stays.
        assert_eq!(b("UIn"), (0.0, -4.0, true));
        assert_eq!(b("UPn"), (0.0, -4.0, false));
        assert_eq!(b("LIUI"), (-5.0, 5.0, true));
        assert_eq!(b("MIUP"), (-inf, 6.0, false));
    }

    /// Regression: a `UI`
    /// line FIRST in BOUNDS used to error, later ones were swallowed as a
    /// bound-set name. HiGHS's writer emits general integers as marker
    /// columns with `LI`/`UI` bounds.
    #[test]
    fn ui_and_li_on_marker_columns_as_highs_writes_them() {
        let text = "\
NAME HWRITE
ROWS
 N  obj
 L  r0
COLUMNS
    MARKER                 'MARKER'                 'INTORG'
    c0        obj       1   r0        1
    c1        obj       1   r0        1
    c2        obj       1   r0        1
    MARKER                 'MARKER'                 'INTEND'
RHS
    RHS_V     r0        10
BOUNDS
 UI BOUND     c0        7
 LI BOUND     c1        -3
 UI BOUND     c1        12
 LI BOUND     c2        2
ENDATA
";
        let m = read_mps(text).unwrap();
        let lu: Vec<(f64, f64, bool)> = m
            .bounds
            .iter()
            .map(|b| (b.lower, b.upper, b.integer))
            .collect();
        assert_eq!(
            lu,
            [
                (0.0, 7.0, true),
                (-3.0, 12.0, true),
                (2.0, f64::INFINITY, true)
            ]
        );
    }

    #[test]
    fn li_ui_sc_need_a_value() {
        let head = "\
NAME E
ROWS
  N  OF
  L  C1
COLUMNS
  X  OF  1  C1  1
BOUNDS
";
        for (line, what) in [
            (" LI bnd X", "requires a value"),
            (" UI bnd X", "requires a value"),
            (" SC bnd X", "requires a value"),
        ] {
            let e = read_mps(&format!("{head}{line}\nENDATA\n")).unwrap_err();
            assert!(e.message.contains(what), "{line}: {e}");
            assert_eq!(e.line, 8, "{line}");
        }
    }

    /// PL (like any BOUNDS line) clears the marker-binary default: a marker
    /// integer with `PL` is a general integer on `[0, +∞)`.
    #[test]
    fn pl_on_marker_column_is_a_general_integer() {
        let text = "\
NAME PLM
ROWS
  N  OF
  L  R1
COLUMNS
  M  'MARKER' 'INTORG'
  A  OF  1  R1  1
  B  OF  1  R1  1
  M  'MARKER' 'INTEND'
BOUNDS
 PL bnd A
ENDATA
";
        let m = read_mps(text).unwrap();
        assert_eq!((m.bounds[0].lower, m.bounds[0].upper), (0.0, f64::INFINITY));
        assert_eq!((m.bounds[1].lower, m.bounds[1].upper), (0.0, 1.0));
        assert!(m.bounds.iter().all(|b| b.integer));
    }

    #[test]
    fn parses_intorg_markers_and_ranges() {
        let text = "\
NAME MARK
ROWS
  N  OF
  L  R1
COLUMNS
  X1  OF  1  R1  1
  X2  'MARKER' 'INTORG'
  X2  OF  1  R1  1
  X2  'MARKER' 'INTEND'
RANGES
  rhs R1 2.5
RHS
  rhs R1 8
ENDATA
";
        let m = read_mps(text).unwrap();
        assert!(m.bounds[1].integer);
        assert!(!m.bounds[0].integer);
        // marker integer without BOUNDS entry defaults to binary
        assert_eq!((m.bounds[1].lower, m.bounds[1].upper), (0.0, 1.0));
        assert_eq!(m.bounds[0].upper, f64::INFINITY);
        assert_eq!(m.ranges, [2.5]);
        assert_eq!(m.rhs, [8.0]);
    }

    /// `Model::to_sparse_lp`: ranged MPS rows expand base-in-place +
    /// complement-at-tail (the historical `model_to_lp` order), with the
    /// column entries replicated into both rows. Covers every ranged shape:
    /// L, G, E r>0, E r<0.
    #[test]
    fn to_sparse_lp_expands_ranges_in_historical_row_order() {
        let text = "\
NAME RANGED
ROWS
  N  OF
  L  RL
  G  RG
  E  REP
  E  REM
COLUMNS
  X1  OF  1  RL  1
  X1  RG  2  REP  3
  X1  REM  4
  X2  RL  5  RG  6
  X2  REP  7  REM  8
RHS
  rhs RL 10  RG 1
  rhs REP 2  REM 3
RANGES
  rng RL 4   RG 2
  rng REP 1.5  REM -0.5
ENDATA
";
        let m = read_mps(text).unwrap();
        let sp = m.to_sparse_lp().unwrap();
        let dense = sp.to_lp();
        // Base rows keep their slots with the ranged-E sense flip
        // (E r>0 ⇒ base G; E r<0 ⇒ base L); complements append in row order.
        assert_eq!(dense.senses, vec!['L', 'G', 'G', 'L', 'G', 'L', 'L', 'G']);
        assert_eq!(dense.rhs, vec![10.0, 1.0, 2.0, 3.0, 6.0, 3.0, 3.5, 2.5]);
        // Entries replicate into the complement row of their source row:
        // HiGHS backport = [RL, RG, REP, REM, RL', RG', REP', REM'] with source values.
        let col0: Vec<f64> = (0..8).map(|i| dense.a.get(i, 0)).collect();
        assert_eq!(col0, vec![1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0]);
        let col1: Vec<f64> = (0..8).map(|i| dense.a.get(i, 1)).collect();
        assert_eq!(col1, vec![5.0, 6.0, 7.0, 8.0, 5.0, 6.0, 7.0, 8.0]);
        // And the SparseLp rows themselves are single-sided (or equality),
        // so a solver-side expansion is the identity on them.
        assert_eq!(sp.a.nrow(), 8);
        for (r, (&lo, &hi)) in sp.row_lower.iter().zip(&sp.row_upper).enumerate() {
            assert!(
                lo == hi || lo.is_infinite() || hi.is_infinite(),
                "row {r} must be single-sided, got [{lo}, {hi}]"
            );
        }
    }

    #[test]
    fn to_sparse_lp_rejects_integer_models() {
        let text = "\
NAME INT
ROWS
  N  OF
  L  R1
COLUMNS
  X1  OF  1  R1  1
  X2  'MARKER' 'INTORG'
  X2  OF  1  R1  1
  X2  'MARKER' 'INTEND'
RHS
  rhs R1 8
ENDATA
";
        let m = read_mps(text).unwrap();
        assert_eq!(
            m.to_sparse_lp().unwrap_err(),
            "integer variables not supported yet (the LP pipeline is continuous-only)"
        );
    }

    #[test]
    fn marker_binary_default_cleared_by_any_bound() {
        let text = "\
NAME MB
ROWS
  N  OF
  L  R1
COLUMNS
  M  'MARKER' 'INTORG'
  A  OF  1  R1  1
  B  OF  1  R1  1
  C  OF  1  R1  1
  M  'MARKER' 'INTEND'
RHS
  rhs R1 8
BOUNDS
 UP bnd B 5
 LO bnd C 2
ENDATA
";
        let m = read_mps(text).unwrap();
        let lu: Vec<(f64, f64)> = m.bounds.iter().map(|b| (b.lower, b.upper)).collect();
        assert_eq!(lu, [(0.0, 1.0), (0.0, 5.0), (2.0, f64::INFINITY)]);
        assert!(m.bounds.iter().all(|b| b.integer));
    }

    /// An RHS entry on the objective (N) row is the objective
    /// constant, HiGHS-style obj_offset = −val, first entry wins
    /// (HMpsFF.cpp:1065–1084). Reporting-only — rhs itself stays constraint-only.
    #[test]
    fn n_row_rhs_becomes_objective_offset() {
        let text = "\
NAME OBJCONST
ROWS
  N  OF
  L  R1
COLUMNS
  X  OF  1  R1  1
RHS
  s1 OF 5.0  R1 8.0
  s1 OF 100.0
ENDATA
";
        let m = read_mps(text).unwrap();
        assert_eq!(m.obj_offset, -5.0, "first N-row entry wins, negated");
        assert_eq!(m.rhs, [8.0]);
    }

    /// The spec keeps only the first RHS/RANGES/BOUNDS set; entries
    /// under a later set name (including a would-be N-row constant) are
    /// ignored. Corpus-wide a no-op (verified over Netlib and MIPLIB: no
    /// instance carries >1 set name) but spec-correct.
    #[test]
    fn only_the_first_rhs_ranges_bounds_set_is_applied() {
        let text = "\
NAME TWOST
ROWS
  N  OF
  L  R1
  L  R2
COLUMNS
  X  OF  1  R1  1  R2  1
RHS
  A R1 8
  B R1 9  R2 7  OF 50
RANGES
  A R1 2
  B R1 3
BOUNDS
  UP s1 X 4
  UP s2 X 9
ENDATA
";
        let m = read_mps(text).unwrap();
        assert_eq!(m.rhs, [8.0, 0.0], "set B rhs ignored");
        assert_eq!(m.ranges, [2.0, 0.0], "set B ranges ignored");
        assert_eq!(m.bounds[0].upper, 4.0, "set s2 bounds ignored");
        assert_eq!(
            m.obj_offset, 0.0,
            "set B N-row entry must not set a constant"
        );
    }

    /// Corpus pin: e226's single RHS set ZZZZZZ01 carries −7.113 on the
    /// N row ...000 → obj_offset +7.113; with our raw c·x = −18.75192906637053
    /// the reported objective becomes HiGHS's −11.638929.
    #[test]
    fn e226_carries_the_objective_constant() {
        let path = "../../benches/data/netlib/e226.mps";
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let m = read_mps(&text).unwrap();
        assert_eq!(m.obj_name.as_deref(), Some("...000"));
        assert_eq!(m.obj_offset, 7.113);
    }

    #[test]
    fn round_trip_preserves_model_semantics() {
        let m = read_mps(SMALL_LP).unwrap();
        let text = write_mps(&m);
        let m2 = read_mps(&text).unwrap();

        assert_eq!(m2.name, m.name);
        assert_eq!(m2.obj_sense, m.obj_sense);
        assert_eq!(m2.senses, m.senses);
        assert_eq!(m2.rhs, m.rhs);
        assert_eq!(m2.ranges, m.ranges);
        assert_eq!(m2.bounds, m.bounds);
        for j in 0..m.num_variables {
            assert!((m2.objective[j] - m.objective[j]).abs() < TOL);
        }
        assert_eq!(m2.matrix.nrow(), m.matrix.nrow());
        assert_eq!(m2.matrix.ncol(), m.matrix.ncol());
        for j in 0..m.matrix.ncol() {
            for i in 0..m.matrix.nrow() {
                assert!(
                    (m2.matrix.get(i, j) - m.matrix.get(i, j)).abs() < TOL,
                    "({i},{j}): {} vs {}",
                    m2.matrix.get(i, j),
                    m.matrix.get(i, j)
                );
            }
        }
    }

    /// Round-trip: `write_mps` re-emits the objective constant as an
    /// N-row RHS entry, so read→write→read preserves `obj_offset` bit-exactly.
    #[test]
    fn round_trip_preserves_objective_offset() {
        let text = "\
NAME OBJRT
ROWS
  N  OF
  L  R1
COLUMNS
  X  OF  1  R1  1
RHS
  s1 OF -7.113  R1 8.0
ENDATA
";
        let m = read_mps(text).unwrap();
        assert_eq!(m.obj_offset, 7.113);
        let out = write_mps(&m);
        assert!(out.contains("rhs OF -7.113"), "{out}");
        let m2 = read_mps(&out).unwrap();
        assert_eq!(m2.obj_offset.to_bits(), m.obj_offset.to_bits());
        assert_eq!(m2.rhs, m.rhs);
    }

    #[test]
    fn round_trip_covers_all_bound_kinds() {
        let text = "\
NAME RT
ROWS
  N  OF
  G  R1
  E  R2
COLUMNS
  A  OF  1  R1  1
  B  OF  1  R1  1  R2  1
  C  OF  1  R2  1
  D  OF  1  R1  1
RHS
  rhs R1 3 R2 1
RANGES
  rng R2 0.5
BOUNDS
  FR A
  BV B
  LI C 2
  UP D 9
ENDATA
";
        let m = read_mps(text).unwrap();
        let m2 = read_mps(&write_mps(&m)).unwrap();
        assert_eq!(m2.bounds, m.bounds);
        assert_eq!((m.bounds[2].lower, m.bounds[2].integer), (2.0, true));
        assert_eq!(m2.ranges, m.ranges);
        assert_eq!(m2.senses, m.senses);
        for j in 0..m.matrix.ncol() {
            for i in 0..m.matrix.nrow() {
                assert!((m2.matrix.get(i, j) - m.matrix.get(i, j)).abs() < TOL);
            }
        }
    }

    /// Every integer box shape (binary, default `[0, ∞)`, fixed, free,
    /// half-open, general) plus an entry-less column survives write → read.
    #[test]
    fn round_trip_covers_integer_boxes() {
        let text = "\
NAME RTI
ROWS
  N  OF
  L  R1
COLUMNS
  Y  OF  1  R1  1
  M  'MARKER' 'INTORG'
  A  OF  1  R1  1
  B  OF  1  R1  1
  C  OF  1  R1  1
  D  OF  1  R1  1
  E  OF  1  R1  1
  F  OF  1  R1  1
  G  OF  0
  M  'MARKER' 'INTEND'
  Z  OF  1
BOUNDS
 PL bnd B
 FX bnd C 4
 FR bnd D
 MI bnd E
 UI bnd E 3
 LI bnd F -2
 UI bnd F 8
 UI bnd G 5
 UP bnd Z 2.5
ENDATA
";
        let m = read_mps(text).unwrap();
        let out = write_mps(&m);
        let m2 = read_mps(&out).unwrap();
        assert_eq!(m2.col_names, m.col_names, "{out}");
        assert_eq!(m2.bounds, m.bounds, "{out}");
        assert_eq!(m2.objective, m.objective);
        let inf = f64::INFINITY;
        let lu: Vec<(f64, f64, bool)> = m
            .bounds
            .iter()
            .map(|b| (b.lower, b.upper, b.integer))
            .collect();
        assert_eq!(
            lu,
            [
                (0.0, inf, false),
                (0.0, 1.0, true),
                (0.0, inf, true),
                (4.0, 4.0, true),
                (-inf, inf, true),
                (-inf, 3.0, true),
                (-2.0, 8.0, true),
                (0.0, 5.0, true),
                (0.0, 2.5, false),
            ]
        );
    }

    #[test]
    fn unknown_section_is_an_error() {
        let text = "\
NAME BAD
FOOBAR
  X Y 1
ROWS
  N OF
ENDATA
";
        let e = read_mps(text).unwrap_err();
        assert!(e.message.contains("unknown section"), "{e}");
        assert_eq!(e.line, 2);
    }

    #[test]
    fn dangling_row_reference_is_an_error() {
        let text = "\
NAME BAD
ROWS
  N OF
  L R1
COLUMNS
  X  OF  1  NOSUCH  2
ENDATA
";
        let e = read_mps(text).unwrap_err();
        assert!(e.message.contains("undeclared row"), "{e}");
    }

    #[test]
    fn dangling_column_in_bounds_is_an_error() {
        let text = "\
NAME BAD
ROWS
  N OF
  L R1
COLUMNS
  X  OF  1  R1  1
BOUNDS
  UP GHOST 5
ENDATA
";
        let e = read_mps(text).unwrap_err();
        assert!(e.message.contains("undeclared column"), "{e}");
    }

    // Quadratic sections, hand-derived from HiGHS HMpsFF.cpp
    // parseQuadRows / parseQuadMatrix (1748-1889) and fillHessian (177-216).
    // The instances are reference/HiGHS/check/instances/qjh*.mps.

    /// The qjh linear part, followed by `quad` and ENDATA.
    fn qjh(quad: &str) -> String {
        format!(
            "\
NAME          qjh
ROWS
 N  obj
 L  c1
COLUMNS
    x1        obj        0.0
    x1        c1        1.0
    x2        obj       -1.0
    x3        obj       -3.0
    x3        c1        1.0
RHS
    rhs       c1        2.0
{quad}ENDATA
"
        )
    }

    /// The square Hessian of qjh in file order: `x1 x3 -1` gives (x3, x1)
    /// then its mirror (x1, x3).
    const QJH_TRIPLETS: [Triplet; 5] = [
        (0, 0, 2.0),
        (2, 0, -1.0),
        (0, 2, -1.0),
        (1, 1, 0.2),
        (2, 2, 2.0),
    ];

    fn assert_qjh_csc(m: &Model) {
        let h = m.fill_hessian();
        assert_eq!(h.dim, 3);
        assert_eq!(h.start, [0, 2, 3, 5]);
        assert_eq!(h.index, [0, 2, 1, 0, 2]);
        assert_eq!(h.value, [2.0, -1.0, 0.2, -1.0, 2.0]);
    }

    #[test]
    fn qp_quadobj_mirrors_off_diagonal() {
        let text = qjh("\
QUADOBJ
    x1        x1        2.0
    x1        x3        -1.0
    x2        x2        0.2
    x3        x3        2.0
");
        let m = read_mps(&text).unwrap();
        assert_eq!(m.num_variables, 3);
        assert_eq!(m.num_constraints, 1);
        assert_eq!(m.hessian, QJH_TRIPLETS);
        assert_qjh_csc(&m);
    }

    #[test]
    fn qp_qmatrix_is_the_full_matrix() {
        let text = qjh("\
QMATRIX
    x1        x1        2.0
    x1        x3        -1.0
    x3        x1        -1.0
    x2        x2        0.2
    x3        x3        2.0
");
        let m = read_mps(&text).unwrap();
        // No mirroring: (x3, x1) and (x1, x3) come from the file.
        assert_eq!(m.hessian, QJH_TRIPLETS);
        assert_qjh_csc(&m);
    }

    #[test]
    fn qp_qsection_on_objective_fills_hessian() {
        let text = qjh("\
QSECTION      obj
    x1        x1        2.0
    x1        x3        -1.0
    x2        x2        0.2
    x3        x3        2.0
");
        let m = read_mps(&text).unwrap();
        assert_eq!(m.hessian, QJH_TRIPLETS);
        assert_qjh_csc(&m);
    }

    #[test]
    fn qp_qsection_on_constraint_row_leaves_hessian_empty() {
        // A constraint row fills qrows_entries, which is parsed but not used.
        let text = qjh("\
QSECTION      c1
    x1        x1        2.0
    x1        x3        -1.0
");
        let m = read_mps(&text).unwrap();
        assert!(m.hessian.is_empty());
        assert_eq!(m.fill_hessian(), crate::model::MpsHessian::default());
    }

    #[test]
    fn qp_qsection_on_undefined_or_free_row_is_skipped() {
        for row in ["nosuch", "free"] {
            let text = format!(
                "\
NAME q
ROWS
 N  obj
 N  free
 L  c1
COLUMNS
    x1  obj  1.0  c1  1.0
QSECTION {row}
    x1  x1  2.0
    ghost  ghost  1.0
QUADOBJ
    x1  x1  3.0
ENDATA
"
            );
            let m = read_mps(&text).unwrap();
            // The skipped lines add no column; the next section is parsed.
            assert_eq!(m.num_variables, 1, "{row}");
            assert_eq!(m.hessian, [(0, 0, 3.0)], "{row}");
        }
    }

    #[test]
    fn qp_qsection_without_row_name_is_an_error() {
        let e = read_mps(&qjh("QSECTION\n    x1  x1  2.0\n")).unwrap_err();
        assert_eq!(e.message, "No row name given in argument of QSECTION");
        assert_eq!(e.line, 13);
    }

    #[test]
    fn qp_missing_coefficient_is_an_error() {
        let e = read_mps(&qjh("QUADOBJ\n    x1  x1  2.0  x2\n")).unwrap_err();
        assert_eq!(
            e.message,
            "QUADOBJ has no coefficient for entry \"x2\" in column \"x1\""
        );
        assert_eq!(e.line, 14);
        let e = read_mps(&qjh("QMATRIX\n    x1  x2\n")).unwrap_err();
        assert!(e.message.starts_with("QMATRIX has no coefficient"), "{e}");
    }

    #[test]
    fn qp_two_entries_per_line_zeros_skipped_duplicates_kept() {
        // Up to two (row, value) pairs per line, the rest is ignored; a zero
        // value is not stored; a repeated entry is kept, not summed.
        let text = qjh("\
QUADOBJ
    x1  x1  2.0  x2  0.5  x3  9.0
    x2  x2  0.0  x3  1.5
    x1  x1  1.0
    x3
");
        let m = read_mps(&text).unwrap();
        assert_eq!(
            m.hessian,
            [
                (0, 0, 2.0),
                (1, 0, 0.5),
                (0, 1, 0.5),
                (2, 1, 1.5),
                (1, 2, 1.5),
                (0, 0, 1.0),
            ]
        );
        let h = m.fill_hessian();
        assert_eq!(h.dim, 3);
        assert_eq!(h.start, [0, 3, 5, 6]);
        assert_eq!(h.index, [0, 1, 0, 0, 2, 1]);
        assert_eq!(h.value, [2.0, 0.5, 1.0, 0.5, 1.5, 1.5]);
    }

    #[test]
    fn qp_unknown_column_is_added() {
        // getColIdx(name) with add_if_new = true: a continuous column with
        // bounds [0, inf) and no cost.
        let text = qjh("QUADOBJ\n    x4  x1  1.0\n");
        let m = read_mps(&text).unwrap();
        assert_eq!(m.num_variables, 4);
        assert_eq!(m.col_names[3], "x4");
        assert_eq!(m.objective, [0.0, -1.0, -3.0, 0.0]);
        assert_eq!(m.bounds[3], VariableBound::default());
        assert_eq!(m.matrix.ncol(), 4);
        assert_eq!(m.hessian, [(0, 3, 1.0), (3, 0, 1.0)]);
        let h = m.fill_hessian();
        assert_eq!(h.dim, 4);
        assert_eq!(h.start, [0, 1, 1, 1, 2]);
        assert_eq!(h.index, [3, 0]);
    }

    #[test]
    fn qp_keyword_followed_by_text_is_data() {
        // checkFirstWord: QUADOBJ / QMATRIX with text after them are data.
        let text = "\
NAME q
ROWS
 N  obj
COLUMNS
    QUADOBJ  obj  1.0
    QMATRIX  obj  2.0
QUADOBJ
    QUADOBJ  QMATRIX  4.0
ENDATA
";
        let m = read_mps(text).unwrap();
        assert_eq!(m.col_names, ["QUADOBJ", "QMATRIX"]);
        assert_eq!(m.objective, [1.0, 2.0]);
        assert_eq!(m.hessian, [(1, 0, 4.0), (0, 1, 4.0)]);
    }

    const CONIC: &str = "\
NAME conic
ROWS
 N  obj
 L  q1
 E  lin
COLUMNS
    x  obj  1  q1  1
    x  lin  1
    y  obj  1  q1  0
    y  lin  1
    t  obj  0  lin  1
RHS
    rhs  q1  5  lin  2
BOUNDS
 UP bnd x 10
 UP bnd y 10
QCMATRIX q1
    x  x  1
    x  x  1
    x  y  1.5
    y  x  1.5
    y  y  1
CSECTION k1 0.0 QUAD
    t
    x
    y
CSECTION k2 0.0 RQUAD
    x
    y
    t
ENDATA
";

    #[test]
    fn qcmatrix_and_csection_fill_the_conic_block() {
        let m = read_mps_conic(CONIC).unwrap();
        // 2x^2 + 3xy + y^2 monomials, SCIP's reading of the same file.
        assert_eq!(m.conic.quad_rows.len(), 1);
        let q = &m.conic.quad_rows[0];
        assert_eq!(q.row, 0);
        assert_eq!(
            q.monomials().collect::<Vec<_>>(),
            [(0, 0, 2.0), (1, 0, 3.0), (1, 1, 1.0)]
        );
        assert_eq!(m.conic.cones.len(), 2);
        assert_eq!(m.conic.cones[0].kind, ConeKind::Soc);
        assert_eq!(m.conic.cones[0].vars, [2, 0, 1]);
        assert_eq!(m.conic.cones[1].kind, ConeKind::RotatedSoc);
    }

    #[test]
    fn legacy_readers_reject_conic_models_by_name() {
        let e = read_mps(CONIC).unwrap_err();
        assert_eq!(e.message, crate::conic::CONIC_CONTENT_MSG);
        assert!(crate::conic::conic_rejected());
    }

    #[test]
    fn conic_write_read_round_trip() {
        let m = read_mps_conic(CONIC).unwrap();
        let m2 = read_mps_conic(&write_mps(&m)).unwrap();
        assert_eq!(m2.conic, m.conic);
    }

    #[test]
    fn conic_sections_report_typed_errors() {
        let bad = [
            (
                CONIC.replace("QCMATRIX q1", "QCMATRIX nope"),
                "undeclared row",
            ),
            (
                CONIC.replace("QCMATRIX q1", "QCMATRIX obj"),
                "objective row",
            ),
            (CONIC.replace("RQUAD", "PEXP"), "exponential"),
            (CONIC.replace("RQUAD", "DPOW"), "power"),
            (
                CONIC.replace("    y  y  1\n", "    y  zz  1\n"),
                "undeclared column 'zz'",
            ),
            (
                CONIC.replace("    t\n    x\n    y\nCSECTION k2", "    t\nCSECTION k2"),
                "needs at least 2",
            ),
        ];
        for (text, want) in bad {
            let e = read_mps_conic(&text).unwrap_err();
            assert!(e.message.contains(want), "{want}: {e}");
        }
    }

    #[test]
    fn star_led_names_are_data_not_comments() {
        // irish-electricity (MIPLIB 2017) names its sets `*RHS*` and `*BND*`;
        // only a `*` in column 1 is a comment (HMpsFF.cpp:226).
        let text = "\
NAME     STAR
ROWS
 N  *OBJ*
 G  DEMAND
COLUMNS
    X  *OBJ*  1.0  DEMAND  1.0
RHS
    *RHS*  DEMAND  5.0
BOUNDS
 UP *BND*  X  9.0
* a real comment
ENDATA
";
        let m = read_mps(text).unwrap();
        assert_eq!(m.rhs, [5.0]);
        assert_eq!(m.objective, [1.0]);
        assert_eq!(m.obj_name.as_deref(), Some("*OBJ*"));
    }

    #[test]
    fn qp_lp_has_empty_hessian() {
        let m = read_mps(SMALL_LP).unwrap();
        assert!(m.hessian.is_empty());
        let h = m.fill_hessian();
        assert_eq!(h.dim, 0);
        assert!(h.start.is_empty() && h.index.is_empty() && h.value.is_empty());
    }
}
