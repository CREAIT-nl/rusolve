//! Dialect regression: one small MILP written by each tool's own LP / MPS
//! (and AMPL .nl) writer must read back to the model the tool meant.
//!
//! The files in `tests/dialects/` are the same model as written by the LP/MPS
//! writers of common modelling tools and solvers, so the readers stay compatible
//! with the dialects found in the wild; they are kept verbatim apart from
//! timestamp comment lines.
//!
//! The model: max 3x + 2y − z + 4b + 0.5w + v + 7 over
//! c1: x + y + z ≤ 10, c2: 2 ≤ x − y ≤ 8, c3: x + b = 3, c4: y + v ≥ −5,
//! c5: z − y ≥ −4; x integer in [0, 6], y in [−2, 5], z free, b binary,
//! w fixed at 2, v in (−∞, 3].
//!
//! Writers legitimately differ in how they carry the range on c2 (one ranged
//! row, two rows, or an explicit slack column), the objective constant (an
//! offset, or a column fixed at 1 with cost 7), and some lose information:
//! CPLEX's MPS writer stores no OBJSENSE (the file is the negated min
//! problem), PuLP drops the constant (and its MPS negates a max objective,
//! the sense only in a comment), and the Xpress file carries only the upper
//! side of c2. Each case states what the file means; the test checks every
//! row as an interval over named terms, every bound, every integrality flag
//! and every objective coefficient.

use rusolve_io::{Model, ObjSense, RowSense, read_model_file};

const INF: f64 = f64::INFINITY;

/// How the file carries `2 ≤ x − y ≤ 8`.
#[derive(Clone, Copy)]
enum C2 {
    /// One ranged row.
    Ranged,
    /// `x − y ≥ 2` and `x − y ≤ 8`.
    Split,
    /// `x − y − s = 2`, `s ∈ [0, 6]` (CPLEX, COPT LP writers).
    SlackLow(&'static str),
    /// `x − y + s = 8`, `s ∈ [0, 6]` (Gurobi).
    SlackHigh(&'static str),
    /// Only `x − y ≤ 8` (the Xpress file has no lower side).
    UpperOnly,
}

/// How the file carries the objective constant 7.
#[derive(Clone, Copy)]
enum Constant {
    /// `obj_offset` (with this value).
    Offset(f64),
    /// A column with cost 7 and bounds [1, 1].
    FixedColumn(&'static str),
    /// A column with cost 7, default bounds [0, ∞), and a row `col = 1`
    /// (Pyomo's MPS writer).
    RowFixedColumn(&'static str),
}

struct Case {
    file: &'static str,
    sense: ObjSense,
    /// +1: the objective as written in the model; −1: negated (min of −f).
    obj_sign: f64,
    c2: C2,
    constant: Constant,
}

fn case(file: &'static str, c2: C2, constant: Constant) -> Case {
    Case {
        file,
        sense: ObjSense::Maximize,
        obj_sign: 1.0,
        c2,
        constant,
    }
}

type Row = (Vec<(String, f64)>, f64, f64);

fn row(terms: &[(&str, f64)], lo: f64, hi: f64) -> Row {
    let mut t: Vec<(String, f64)> = terms.iter().map(|&(n, c)| (n.to_string(), c)).collect();
    t.sort_by(|a, b| a.0.cmp(&b.0));
    (t, lo, hi)
}

/// Every constraint row as (terms sorted by column name, lower, upper),
/// using the MPS range semantics `Model` documents.
fn rows_of(m: &Model) -> Vec<Row> {
    let mut terms: Vec<Vec<(String, f64)>> = vec![Vec::new(); m.num_constraints];
    let (ptr, idx, val) = (m.matrix.col_ptr(), m.matrix.row_idx(), m.matrix.values());
    for j in 0..m.num_variables {
        for k in ptr[j]..ptr[j + 1] {
            terms[idx[k] as usize].push((m.col_names[j].clone(), val[k]));
        }
    }
    let mut out = Vec::new();
    for (i, mut t) in terms.into_iter().enumerate() {
        t.sort_by(|a, b| a.0.cmp(&b.0));
        let (b, r) = (m.rhs[i], m.ranges[i]);
        let (lo, hi) = match m.senses[i] {
            RowSense::LessEqual if r != 0.0 => (b - r.abs(), b),
            RowSense::LessEqual => (-INF, b),
            RowSense::GreaterEqual if r != 0.0 => (b, b + r.abs()),
            RowSense::GreaterEqual => (b, INF),
            RowSense::Equal if r > 0.0 => (b, b + r),
            RowSense::Equal if r < 0.0 => (b + r, b),
            RowSense::Equal => (b, b),
        };
        out.push((t, lo, hi));
    }
    out
}

fn sorted(mut v: Vec<Row>) -> Vec<Row> {
    v.sort_by(|a, b| format!("{a:?}").cmp(&format!("{b:?}")));
    v
}

fn col(m: &Model, name: &str) -> usize {
    m.col_names
        .iter()
        .position(|n| n == name)
        .unwrap_or_else(|| panic!("no column {name} in {:?}", m.col_names))
}

fn check(c: &Case) {
    let path = format!("{}/tests/dialects/{}", env!("CARGO_MANIFEST_DIR"), c.file);
    let m = read_model_file(&path).unwrap_or_else(|e| panic!("{}: {e}", c.file));
    let f = c.file;

    // Rows.
    let mut want = vec![
        row(&[("x", 1.0), ("y", 1.0), ("z", 1.0)], -INF, 10.0),
        row(&[("x", 1.0), ("b", 1.0)], 3.0, 3.0),
        row(&[("y", 1.0), ("v", 1.0)], -5.0, INF),
        row(&[("y", -1.0), ("z", 1.0)], -4.0, INF),
    ];
    let mut extra_cols: Vec<(&str, f64, f64, f64)> = Vec::new(); // name, cost, lo, hi
    match c.c2 {
        C2::Ranged => want.push(row(&[("x", 1.0), ("y", -1.0)], 2.0, 8.0)),
        C2::Split => {
            want.push(row(&[("x", 1.0), ("y", -1.0)], 2.0, INF));
            want.push(row(&[("x", 1.0), ("y", -1.0)], -INF, 8.0));
        }
        C2::SlackLow(s) => {
            want.push(row(&[("x", 1.0), ("y", -1.0), (s, -1.0)], 2.0, 2.0));
            extra_cols.push((s, 0.0, 0.0, 6.0));
        }
        C2::SlackHigh(s) => {
            want.push(row(&[("x", 1.0), ("y", -1.0), (s, 1.0)], 8.0, 8.0));
            extra_cols.push((s, 0.0, 0.0, 6.0));
        }
        C2::UpperOnly => want.push(row(&[("x", 1.0), ("y", -1.0)], -INF, 8.0)),
    }
    let mut offset = 0.0;
    match c.constant {
        Constant::Offset(v) => offset = v,
        Constant::FixedColumn(n) => extra_cols.push((n, 7.0, 1.0, 1.0)),
        Constant::RowFixedColumn(n) => {
            extra_cols.push((n, 7.0, 0.0, INF));
            want.push(row(&[(n, 1.0)], 1.0, 1.0));
        }
    }
    assert_eq!(m.num_constraints, want.len(), "{f}: row count");
    assert_eq!(sorted(rows_of(&m)), sorted(want), "{f}: rows");

    // Columns: the six model variables plus the file's extra columns.
    assert_eq!(
        m.num_variables,
        6 + extra_cols.len(),
        "{f}: column count {:?}",
        m.col_names
    );
    assert_eq!(m.obj_sense, c.sense, "{f}: sense");
    assert_eq!(m.obj_offset, offset, "{f}: objective offset");
    let vars: [(&str, f64, f64, f64, bool); 6] = [
        ("x", 3.0, 0.0, 6.0, true),
        ("y", 2.0, -2.0, 5.0, false),
        ("z", -1.0, -INF, INF, false),
        ("b", 4.0, 0.0, 1.0, true),
        ("w", 0.5, 2.0, 2.0, false),
        ("v", 1.0, -INF, 3.0, false),
    ];
    for (n, cost, lo, hi, int) in vars {
        let j = col(&m, n);
        assert_eq!(m.objective[j], c.obj_sign * cost, "{f}: cost of {n}");
        let b = m.bounds[j];
        assert_eq!(
            (b.lower, b.upper, b.integer),
            (lo, hi, int),
            "{f}: bounds of {n}"
        );
    }
    for (n, cost, lo, hi) in extra_cols {
        let j = col(&m, n);
        assert_eq!(m.objective[j], cost, "{f}: cost of {n}");
        let b = m.bounds[j];
        assert_eq!(
            (b.lower, b.upper, b.integer),
            (lo, hi, false),
            "{f}: bounds of {n}"
        );
    }
}

#[test]
fn every_writer_dialect_reads_to_the_model_it_means() {
    use C2::*;
    use Constant::*;
    let cases = [
        case("highs.lp", Split, Offset(7.0)),
        case("highs.mps", Ranged, Offset(7.0)),
        case("scip.lp", Split, Offset(7.0)),
        case("scip.mps", Split, Offset(7.0)),
        case("copt.lp", SlackLow("Rg_c2"), FixedColumn("ObjConstant")),
        case("copt.mps", Ranged, Offset(7.0)),
        case("cplex.lp", SlackLow("Rgc2"), Offset(7.0)),
        // No OBJSENSE in the file: the negated min problem, constant kept.
        Case {
            sense: ObjSense::Minimize,
            obj_sign: -1.0,
            ..case("cplex.mps", Ranged, Offset(7.0))
        },
        case("gurobi.lp", SlackHigh("Rgc2"), FixedColumn("Constant")),
        case("gurobi.mps", SlackHigh("Rgc2"), Offset(7.0)),
        case("xpress.lp", UpperOnly, Offset(7.0)),
        case("xpress.mps", UpperOnly, Offset(7.0)),
        // PuLP drops the constant; its MPS negates the max objective and
        // keeps the sense only in a `*SENSE:Maximize` comment.
        case("pulp.lp", Split, Offset(0.0)),
        Case {
            sense: ObjSense::Minimize,
            obj_sign: -1.0,
            ..case("pulp.mps", Split, Offset(0.0))
        },
        case("pyomo.lp", Split, FixedColumn("ONE_VAR_CONSTANT")),
        case("pyomo.mps", Split, RowFixedColumn("ONE_VAR_CONSTANT")),
        // AMPL .nl (text), names from the pyomo.col / pyomo.row sidecars.
        case("pyomo.nl", Ranged, Offset(7.0)),
        // MathOptFormat, hand-written against the 1.9 schema.
        case("handwritten.mof.json", Ranged, Offset(7.0)),
    ];
    for c in &cases {
        check(c);
    }
}
