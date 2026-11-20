//! Instance classification: MIPLIB 2017 row types, column / objective
//! statistics, cheap structure detectors and a fixed-length feature vector.
//!
//! Purpose: recognise problem structure from the formulation alone (no LP
//! solve), so that (a) `rusolve-cli classify` can describe an instance and
//! (b) a tuning profile can match a new
//! instance to the training set it was tuned on by feature distance.
//!
//! Everything here is a pure function of the [`MipProblem`]; no state is kept
//! between calls.
//!
//! # Provenance
//!
//! * **Row types** — the 17 linear-constraint classes of A. Gleixner, G.
//!   Hendel, G. Gamrath, T. Achterberg, M. Bastubbe, T. Berthold, P.
//!   Christophel, K. Jarck, T. Koch, J. Linderoth, M. Lübbecke, H. D.
//!   Mittelmann, D. Ozyurt, T. K. Ralphs, D. Salvagnin, Y. Shinano, "MIPLIB
//!   2017: data-driven compilation of the 6th mixed-integer programming
//!   library", Math. Prog. Comp. 13:443–490 (2021), §4.4 (the classification
//!   table; it is Table 5 in the 2019 preprint, Optimization Online 7285).
//!   The table is sorted most-specific first and "a constraint is always
//!   assigned the first (topmost) type that applies". Where the table is
//!   informal (how ranged rows, `≥` rows and negative coefficients are
//!   treated) we follow the reference implementation MIPLIB used,
//!   `SCIPclassifyConstraintTypesLinear` in SCIP's `src/scip/cons_linear.c`
//!   (Apache-2.0, © Zuse Institute Berlin; read as a reference text and
//!   reimplemented, not translated). Consequences of following it:
//!   - each finite side of a row is classified, so a ranged row can count
//!     twice (e.g. `set packing` on its `≤` side + `set covering` on its `≥`
//!     side), and singleton / precedence / variable-bound / mixed-binary /
//!     general rows that are ranged count 2;
//!   - negative coefficients on binaries are complemented (`x → 1 − x̄`)
//!     before the set-type and knapsack tests;
//!   - "binary" means an integer column with bounds inside `[0, 1]`.
//!
//!   Two deliberate deviations, both switched off by
//!   [`ClassifyOptions::miplib_strict`] (default `false`):
//!   - SCIP's knapsack and integer-knapsack tests only look at the `≤`
//!     side, so a pure `Σ aᵢxᵢ ≥ b` binary knapsack (e.g. a knapsack-cover
//!     row) falls through to *mixed binary*. By default we also test the
//!     negated `≥` side (`Σ −aᵢxᵢ ≤ −b`, complemented), which the table's
//!     definitions admit after complementation.
//!   - The table defines a variable bound as `ax + by ≤ c, x ∈ {0, 1}`
//!     (current SCIP checks for the binary too), but the *published*
//!     MIPLIB 2017 per-instance counts file every two-term non-precedence
//!     inequality as a variable bound, binary or not (verified on dano3_3:
//!     684 continuous–continuous two-term rows, all published as variable
//!     bounds). Default: require the binary, as the table says.
//!
//!   `miplib_strict = true` reproduces the published MIPLIB 2017 "Original"
//!   counts (all 47 benchmark instances we checked match exactly).
//! * **Embedded network rows** — greedy *row-scanning addition* heuristic of
//!   R. E. Bixby & R. Fourer, "Finding embedded network rows in linear
//!   programs I. Extraction heuristics", Management Science 34(3):342–376
//!   (1988): rows are added in scan order, each possibly reflected (×−1), as
//!   long as every column keeps at most one `+1` and one `−1` among the
//!   accepted rows. Restricted to `±1` entries (pure, unscaled networks).
//! * **Cliques from knapsack rows** — a binary `≤` row (after complementing)
//!   whose two smallest weights among its `k` largest exceed the capacity
//!   makes those `k` literals pairwise conflicting: A. Atamtürk, G. L.
//!   Nemhauser, M. W. P. Savelsbergh, "Conflict graphs in solving integer
//!   programming problems", Eur. J. Oper. Res. 121:40–55 (2000), §3.
//! * **Feature normalisation** — MIPLIB 2017 §4.3 (log scaling for sizes and
//!   dynamism, proportions for variable types and constraint classes, row
//!   normalisation by `‖aᵢ‖∞`). Our features are computed on the *original*
//!   formulation, not after trivial presolve as in MIPLIB.
//!
//! No `unsafe`.

use crate::MipProblem;
use std::collections::HashSet;
use std::fmt::Write as _;

/// Version of [`FEATURE_NAMES`] / [`Classification::features`]. Bump on any
/// change to the list, order or normalisation; profiles store it and refuse
/// to compare across versions.
pub const FEATURE_VERSION: u32 = 1;

/// Absolute magnitude at or above which a bound counts as infinite (SCIP's
/// convention; MPS files often write `1e30`).
const INF: f64 = 1e20;
/// Comparison tolerance, SCIP's default `numerics/epsilon`.
const EPS: f64 = 1e-9;

fn is_inf(v: f64) -> bool {
    v.abs() >= INF
}

fn eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= EPS * a.abs().max(b.abs()).max(1.0)
}

fn integral(v: f64) -> bool {
    (v - v.round()).abs() <= EPS * v.abs().max(1.0)
}

/// MIPLIB 2017 constraint classes, most specific first (the table order).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum RowType {
    Empty,
    Free,
    Singleton,
    Aggregation,
    Precedence,
    VariableBound,
    SetPartitioning,
    SetPacking,
    SetCovering,
    Cardinality,
    InvariantKnapsack,
    EquationKnapsack,
    BinPacking,
    Knapsack,
    IntegerKnapsack,
    MixedBinary,
    GeneralLinear,
}

impl RowType {
    /// All 17 types in table order (the histogram index order).
    pub const ALL: [RowType; 17] = [
        RowType::Empty,
        RowType::Free,
        RowType::Singleton,
        RowType::Aggregation,
        RowType::Precedence,
        RowType::VariableBound,
        RowType::SetPartitioning,
        RowType::SetPacking,
        RowType::SetCovering,
        RowType::Cardinality,
        RowType::InvariantKnapsack,
        RowType::EquationKnapsack,
        RowType::BinPacking,
        RowType::Knapsack,
        RowType::IntegerKnapsack,
        RowType::MixedBinary,
        RowType::GeneralLinear,
    ];

    /// Stable snake_case name (used in JSON and feature names).
    pub fn name(self) -> &'static str {
        match self {
            RowType::Empty => "empty",
            RowType::Free => "free",
            RowType::Singleton => "singleton",
            RowType::Aggregation => "aggregation",
            RowType::Precedence => "precedence",
            RowType::VariableBound => "variable_bound",
            RowType::SetPartitioning => "set_partitioning",
            RowType::SetPacking => "set_packing",
            RowType::SetCovering => "set_covering",
            RowType::Cardinality => "cardinality",
            RowType::InvariantKnapsack => "invariant_knapsack",
            RowType::EquationKnapsack => "equation_knapsack",
            RowType::BinPacking => "bin_packing",
            RowType::Knapsack => "knapsack",
            RowType::IntegerKnapsack => "integer_knapsack",
            RowType::MixedBinary => "mixed_binary",
            RowType::GeneralLinear => "general_linear",
        }
    }

    fn index(self) -> usize {
        self as usize
    }
}

/// Column kind as the MIPLIB classification sees it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarKind {
    /// Integer with bounds inside `[0, 1]`.
    Binary,
    Integer,
    Continuous,
}

/// Knobs for [`classify`].
#[derive(Debug, Clone, Copy, Default)]
pub struct ClassifyOptions {
    /// Reproduce SCIP/MIPLIB 2017 counts exactly: knapsack and integer
    /// knapsack only on the `≤` side (see the module docs).
    pub miplib_strict: bool,
}

/// Class(es) of one row: one entry per classified side. A ranged row whose
/// single class covers both sides has the same type twice (SCIP counting).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RowClass {
    pub sides: [Option<RowType>; 2],
}

impl RowClass {
    /// The topmost type assigned to the row.
    pub fn primary(&self) -> RowType {
        self.sides[0].expect("every row gets at least one class")
    }
}

/// Sign pattern of the minimization-form objective.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjSign {
    Zero,
    NonNegative,
    NonPositive,
    Mixed,
}

impl ObjSign {
    pub fn name(self) -> &'static str {
        match self {
            ObjSign::Zero => "zero",
            ObjSign::NonNegative => "nonnegative",
            ObjSign::NonPositive => "nonpositive",
            ObjSign::Mixed => "mixed",
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ColumnStats {
    pub binary: usize,
    pub integer: usize,
    pub continuous: usize,
    /// Bound types (all columns): both infinite / only lower finite / only
    /// upper finite / both finite and different / fixed.
    pub free: usize,
    pub lower_only: usize,
    pub upper_only: usize,
    pub boxed: usize,
    pub fixed: usize,
    /// Columns without a matrix nonzero.
    pub empty: usize,
    /// Non-binary integers whose domain is unbounded or spans more than
    /// 10⁴ values (where general-integer branching gets hard).
    pub wide_integer: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ObjectiveStats {
    /// Nonzero cost coefficients.
    pub nonzeros: usize,
    pub zero_objective: bool,
    /// Every nonzero cost sits on a binary column (vacuously false if zero).
    pub only_binaries: bool,
    /// Distinct nonzero cost values (relative tolerance 1e-9).
    pub distinct_values: usize,
    /// Sign pattern of the *minimization-form* costs (maximization negated).
    pub sign: ObjSign,
    /// All nonzero costs are integral and sit on integer columns: the optimal
    /// value is integral (objective-integrality cutoff applies).
    pub integral: bool,
    /// `max|c| / min|c|` over nonzeros (1 when ≤ 1 nonzero).
    pub dynamism: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StructureStats {
    /// Equality rows (`lo = hi`, finite).
    pub equality_rows: usize,
    /// Ranged rows (both sides finite, different).
    pub ranged_rows: usize,
    /// Distinct `(y, x)` pairs with `y` non-binary, `x` binary, from 2-term
    /// rows that bound `y` from above / below by an affine function of `x`.
    pub vub_pairs: usize,
    pub vlb_pairs: usize,
    /// Largest `|coef(x)| / |coef(y)|` over those rows: the big-M size.
    pub max_vb_ratio: f64,
    /// Rows accepted by the greedy embedded-network scan.
    pub network_rows: usize,
    /// Set-packing / set-partitioning shaped rows (binary, unit coefficients
    /// after scaling and complementing, capacity 1) with ≥ 2 literals —
    /// includes 2-term edge rows the MIPLIB table files as variable bounds.
    pub clique_rows: usize,
    /// `clique_rows` after removing rows with identical literal sets.
    pub distinct_cliques: usize,
    /// Largest clique row.
    pub max_clique: usize,
    /// Further cliques (≥ 2 literals) extracted from binary knapsack rows.
    pub knapsack_cliques: usize,
    /// Connected components of the row–column graph over non-empty columns.
    pub components: usize,
    /// Columns in the largest component.
    pub largest_component: usize,
    /// `nnz / (rows · cols)`.
    pub density: f64,
    /// Smallest / largest absolute nonzero matrix coefficient (0 if none).
    pub coef_min: f64,
    pub coef_max: f64,
    /// Largest and mean `log10(max|aᵢⱼ| / min|aᵢⱼ|)` over non-empty rows.
    pub row_dynamism_max_log10: f64,
    pub row_dynamism_mean_log10: f64,
}

/// Everything [`classify`] finds out about an instance.
#[derive(Debug, Clone, PartialEq)]
pub struct Classification {
    pub num_rows: usize,
    pub num_cols: usize,
    pub nnz: usize,
    pub row_class: Vec<RowClass>,
    /// Side counts per [`RowType`] in [`RowType::ALL`] order (SCIP counting:
    /// sums to ≥ `num_rows`).
    pub row_histogram: [usize; 17],
    pub columns: ColumnStats,
    pub objective: ObjectiveStats,
    pub structure: StructureStats,
}

/// Names of [`Classification::features`], index-aligned. Normalisation
/// (all values land in `[0, 1]`; `clamp01` caps outliers):
///
/// | features | value |
/// |---|---|
/// | `size.rows`, `size.cols` | `log10(1 + count) / 7` |
/// | `size.nnz` | `log10(1 + nnz) / 8` |
/// | `size.density` | `(log10(density) + 6) / 6`, 0 when no nonzeros |
/// | `size.row_nnz_mean` | `log10(1 + nnz / rows) / 4` |
/// | `vars.*_frac` | share of columns |
/// | `bounds.*_frac` | share of columns |
/// | `obj.zero`, `obj.only_binaries`, `obj.integral` | 0 / 1 |
/// | `obj.density` | nonzero costs / columns |
/// | `obj.distinct` | `log10(1 + distinct) / 6` |
/// | `obj.pos_frac` | share of nonzero min-form costs that are > 0 |
/// | `obj.dynamism` | `log10(dynamism) / 12` |
/// | `rows.<type>` | side count / total side count (sums to 1) |
/// | `rows.equality_frac` | equality rows / rows |
/// | `struct.network_frac` | network rows / non-empty rows |
/// | `struct.clique_frac` | clique rows / rows |
/// | `struct.vb_per_col` | `(vub + vlb pairs) / cols` |
/// | `struct.components` | `log10(components) / 7` |
/// | `struct.largest_component` | largest component columns / cols |
/// | `coef.dynamism` | `log10(coef_max / coef_min) / 12` |
/// | `coef.row_dynamism_max` | `row_dynamism_max_log10 / 12` |
/// | `coef.big_m` | `log10(1 + max_vb_ratio) / 12` |
pub const FEATURE_NAMES: [&str; 44] = [
    "size.rows",
    "size.cols",
    "size.nnz",
    "size.density",
    "size.row_nnz_mean",
    "vars.binary_frac",
    "vars.integer_frac",
    "vars.continuous_frac",
    "bounds.free_frac",
    "bounds.fixed_frac",
    "bounds.boxed_frac",
    "obj.zero",
    "obj.density",
    "obj.only_binaries",
    "obj.distinct",
    "obj.pos_frac",
    "obj.integral",
    "obj.dynamism",
    "rows.empty",
    "rows.free",
    "rows.singleton",
    "rows.aggregation",
    "rows.precedence",
    "rows.variable_bound",
    "rows.set_partitioning",
    "rows.set_packing",
    "rows.set_covering",
    "rows.cardinality",
    "rows.invariant_knapsack",
    "rows.equation_knapsack",
    "rows.bin_packing",
    "rows.knapsack",
    "rows.integer_knapsack",
    "rows.mixed_binary",
    "rows.general_linear",
    "rows.equality_frac",
    "struct.network_frac",
    "struct.clique_frac",
    "struct.vb_per_col",
    "struct.components",
    "struct.largest_component",
    "coef.dynamism",
    "coef.row_dynamism_max",
    "coef.big_m",
];

/// One-line meaning of every entry of [`FEATURE_NAMES`], index-aligned (`rusolve classify
/// --format json` ships them next to the values so an agent need not read this file).
pub const FEATURE_DESCRIPTIONS: [&str; 44] = [
    "number of rows, log10(1 + rows) / 7",
    "number of columns, log10(1 + cols) / 7",
    "number of nonzeros, log10(1 + nnz) / 8",
    "matrix density, (log10(nnz / (rows * cols)) + 6) / 6; 0 without nonzeros",
    "mean nonzeros per row, log10(1 + nnz / rows) / 4",
    "share of columns that are binary (integer within [0, 1])",
    "share of columns that are general integer (not binary)",
    "share of continuous columns",
    "share of columns with no finite bound",
    "share of columns with lower bound equal to upper bound",
    "share of columns with both bounds finite and different",
    "1 when the objective is all zero (feasibility problem), else 0",
    "share of columns with a nonzero cost",
    "1 when every nonzero cost sits on a binary column, else 0",
    "number of distinct nonzero cost values, log10(1 + distinct) / 6",
    "share of nonzero costs that are positive in minimisation form",
    "1 when all nonzero costs are integral on integer columns, else 0",
    "cost magnitude spread, log10(max |c| / min |c|) / 12",
    "share of row sides classified as empty rows (MIPLIB 2017 row type)",
    "share of row sides classified as free rows",
    "share of row sides classified as singleton rows",
    "share of row sides classified as aggregation rows",
    "share of row sides classified as precedence rows",
    "share of row sides classified as variable-bound rows",
    "share of row sides classified as set-partitioning rows",
    "share of row sides classified as set-packing rows",
    "share of row sides classified as set-covering rows",
    "share of row sides classified as cardinality rows",
    "share of row sides classified as invariant-knapsack rows",
    "share of row sides classified as equation-knapsack rows",
    "share of row sides classified as bin-packing rows",
    "share of row sides classified as knapsack rows",
    "share of row sides classified as integer-knapsack rows",
    "share of row sides classified as mixed-binary rows",
    "share of row sides classified as general linear rows",
    "equality rows / rows",
    "embedded network rows / non-empty rows (Bixby-Fourer row scanning, +-1 entries)",
    "rows that are cliques (also found in knapsack rows) / rows",
    "variable upper plus lower bound pairs per column",
    "connected components of the column-row graph, log10(components) / 7",
    "columns of the largest connected component / columns",
    "matrix coefficient spread, log10(coef_max / coef_min) / 12",
    "largest within-row coefficient spread, log10 scale / 12",
    "big-M size: log10(1 + largest variable-bound ratio) / 12",
];

/// Column kind of every column of `p`.
pub fn var_kinds(p: &MipProblem) -> Vec<VarKind> {
    (0..p.num_col)
        .map(|j| {
            if !p.integer[j] {
                VarKind::Continuous
            } else if (p.col_lower[j] - EPS).ceil() >= 0.0 && (p.col_upper[j] + EPS).floor() <= 1.0
            {
                VarKind::Binary
            } else {
                VarKind::Integer
            }
        })
        .collect()
}

/// One row with duplicates merged and |a| ≤ EPS entries dropped (SCIP
/// `mergeMultiples` + `isZero`), sorted by column.
struct Row {
    cols: Vec<usize>,
    vals: Vec<f64>,
}

fn rows_of(p: &MipProblem) -> Vec<Row> {
    let csr = p.matrix.to_csr();
    (0..p.num_row)
        .map(|i| {
            let (s, e) = (csr.row_ptr()[i], csr.row_ptr()[i + 1]);
            let mut cols: Vec<usize> = Vec::with_capacity(e - s);
            let mut vals: Vec<f64> = Vec::with_capacity(e - s);
            for k in s..e {
                let j = csr.col_idx()[k] as usize;
                let v = csr.values()[k];
                if cols.last() == Some(&j) {
                    *vals.last_mut().expect("paired with cols") += v;
                } else {
                    cols.push(j);
                    vals.push(v);
                }
            }
            let (cols, vals) = cols
                .into_iter()
                .zip(vals)
                .filter(|&(_, v)| v.abs() > EPS)
                .unzip();
            Row { cols, vals }
        })
        .collect()
}

/// Knapsack-stage test on one `≤` side `Σ sign·aᵢxᵢ ≤ bound` (SCIP's
/// EQKNAPSACK/BINPACKING/KNAPSACK block): all columns binary, coefficients
/// integral; capacity after complementing must be finite, ≥ 0, integral.
/// Returns the complemented capacity on success.
fn knapsack_capacity(row: &Row, kinds: &[VarKind], sign: f64, bound: f64) -> Option<f64> {
    if is_inf(bound) {
        return None;
    }
    let mut b = bound;
    for (&j, &a) in row.cols.iter().zip(&row.vals) {
        let a = sign * a;
        if kinds[j] != VarKind::Binary || !integral(a) {
            return None;
        }
        if a < -EPS {
            b -= a;
        }
    }
    (b > -EPS && integral(b)).then_some(b)
}

fn bin_packing_or_knapsack(row: &Row, b: f64) -> RowType {
    if row.vals.iter().any(|&a| eq(b, a.abs())) {
        RowType::BinPacking
    } else {
        RowType::Knapsack
    }
}

/// Integer-knapsack test on one `≤` side: integral columns with lb ≥ 0,
/// nonnegative integral coefficients, finite nonnegative integral bound.
fn is_int_knapsack(row: &Row, p: &MipProblem, sign: f64, bound: f64) -> bool {
    !is_inf(bound)
        && bound > -EPS
        && integral(bound)
        && row.cols.iter().zip(&row.vals).all(|(&j, &a)| {
            let a = sign * a;
            p.integer[j] && p.col_lower[j] > -EPS && integral(a) && a > -EPS
        })
}

/// Classifies one row following the MIPLIB 2017 table / SCIP reference
/// implementation (module docs). `lo`/`hi` are the row bounds.
fn classify_row(
    row: &Row,
    lo: f64,
    hi: f64,
    p: &MipProblem,
    kinds: &[VarKind],
    opts: ClassifyOptions,
) -> RowClass {
    let mut lhs = if is_inf(lo) && lo < 0.0 {
        f64::NEG_INFINITY
    } else {
        lo
    };
    let mut rhs = if is_inf(hi) && hi > 0.0 {
        f64::INFINITY
    } else {
        hi
    };
    let n = row.cols.len();
    let one = |t| RowClass {
        sides: [Some(t), None],
    };
    let ranged = |lhs: f64, rhs: f64| !(eq(lhs, rhs) || lhs.is_infinite() || rhs.is_infinite());
    let counted = |t, lhs: f64, rhs: f64| RowClass {
        sides: [Some(t), ranged(lhs, rhs).then_some(t)],
    };
    if n == 0 {
        return one(RowType::Empty);
    }
    if lhs.is_infinite() && rhs.is_infinite() {
        return one(RowType::Free);
    }
    if n == 1 {
        return counted(RowType::Singleton, lhs, rhs);
    }
    let is_eq = |lhs: f64, rhs: f64| lhs.is_finite() && rhs.is_finite() && eq(lhs, rhs);
    if n == 2 && is_eq(lhs, rhs) {
        return one(RowType::Aggregation);
    }
    if n == 2 {
        let (k0, k1) = (kinds[row.cols[0]], kinds[row.cols[1]]);
        if k0 == k1 && eq(row.vals[0], -row.vals[1]) {
            return counted(RowType::Precedence, lhs, rhs);
        }
        if k0 == VarKind::Binary || k1 == VarKind::Binary || opts.miplib_strict {
            return counted(RowType::VariableBound, lhs, rhs);
        }
    }

    let mut out: Vec<RowType> = Vec::with_capacity(2);
    let finish = |out: Vec<RowType>| RowClass {
        sides: [out.first().copied(), out.get(1).copied()],
    };

    // Set partitioning / packing / covering, cardinality, invariant knapsack:
    // binary columns with equal |coefficient| (complementing negatives).
    let scale = row.vals[0].abs();
    let unit = row
        .cols
        .iter()
        .zip(&row.vals)
        .all(|(&j, &a)| kinds[j] == VarKind::Binary && eq(a.abs(), scale));
    if unit {
        let nneg = row.vals.iter().filter(|&&a| a < 0.0).count() as f64;
        if is_eq(lhs, rhs) {
            let b = rhs / scale + nneg;
            if eq(1.0, b) {
                return one(RowType::SetPartitioning);
            }
            if integral(b) && b > -EPS {
                return one(RowType::Cardinality);
            }
        }
        if rhs.is_finite() {
            let b = rhs / scale + nneg;
            if eq(1.0, b) {
                out.push(RowType::SetPacking);
                rhs = f64::INFINITY;
            } else if integral(b) && b > -EPS {
                out.push(RowType::InvariantKnapsack);
                rhs = f64::INFINITY;
            }
        }
        if lhs.is_finite() && eq(1.0, lhs / scale + nneg) {
            out.push(RowType::SetCovering);
            lhs = f64::NEG_INFINITY;
        }
        if lhs.is_infinite() && rhs.is_infinite() {
            return finish(out);
        }
    }

    // Equation knapsack / bin packing / knapsack (binary columns).
    if let Some(b) = knapsack_capacity(row, kinds, 1.0, rhs) {
        if is_eq(lhs, rhs) {
            out.push(RowType::EquationKnapsack);
            return finish(out);
        }
        out.push(bin_packing_or_knapsack(row, b));
        if lhs.is_infinite() {
            return finish(out);
        }
        rhs = f64::INFINITY;
    }
    if !opts.miplib_strict && rhs.is_infinite() && lhs.is_finite() {
        // Deviation (module docs): the `≥` side as `Σ −aᵢxᵢ ≤ −lhs`.
        if let Some(b) = knapsack_capacity(row, kinds, -1.0, -lhs) {
            out.push(bin_packing_or_knapsack(row, b));
            return finish(out);
        }
    }

    // Integer knapsack.
    if is_int_knapsack(row, p, 1.0, rhs) {
        out.push(RowType::IntegerKnapsack);
        if lhs.is_infinite() {
            return finish(out);
        }
        rhs = f64::INFINITY;
    }
    if !opts.miplib_strict
        && rhs.is_infinite()
        && lhs.is_finite()
        && is_int_knapsack(row, p, -1.0, -lhs)
    {
        out.push(RowType::IntegerKnapsack);
        return finish(out);
    }

    // Mixed binary (no general integers) or general linear, on what is left.
    let t = if row.cols.iter().all(|&j| kinds[j] != VarKind::Integer) {
        RowType::MixedBinary
    } else {
        RowType::GeneralLinear
    };
    out.push(t);
    if ranged(lhs, rhs) {
        out.push(t);
    }
    finish(out)
}

/// Literal of a binary column: `(column, complemented)`.
type Lit = (usize, bool);

/// The literals of a set-packing / set-partitioning shaped row (≥ 2 of them),
/// written as `Σ literals ≤ 1`, or `None`.
fn clique_of(row: &Row, lo: f64, hi: f64, kinds: &[VarKind]) -> Option<Vec<Lit>> {
    if row.cols.len() < 2 {
        return None;
    }
    let scale = row.vals[0].abs();
    if !row
        .cols
        .iter()
        .zip(&row.vals)
        .all(|(&j, &a)| kinds[j] == VarKind::Binary && eq(a.abs(), scale))
    {
        return None;
    }
    let nneg = row.vals.iter().filter(|&&a| a < 0.0).count() as f64;
    let npos = row.vals.len() as f64 - nneg;
    // ≤ side: Σ_{a>0} x + Σ_{a<0} x̄ ≤ hi/scale + nneg.
    if !is_inf(hi) && eq(1.0, hi / scale + nneg) {
        let lits = row
            .cols
            .iter()
            .zip(&row.vals)
            .map(|(&j, &a)| (j, a < 0.0))
            .collect();
        return Some(lits);
    }
    // ≥ side negated: Σ_{a>0} x̄ + Σ_{a<0} x ≤ npos − lo/scale.
    if !is_inf(lo) && eq(1.0, npos - lo / scale) {
        let lits = row
            .cols
            .iter()
            .zip(&row.vals)
            .map(|(&j, &a)| (j, a > 0.0))
            .collect();
        return Some(lits);
    }
    None
}

/// Size of the clique extracted from a binary `≤` knapsack side
/// (Atamtürk–Nemhauser–Savelsbergh 2000): the largest `k ≥ 2` such that the
/// two smallest of the `k` largest weights sum above the capacity.
fn knapsack_clique_size(row: &Row, kinds: &[VarKind], sign: f64, bound: f64) -> usize {
    if is_inf(bound) || row.cols.iter().any(|&j| kinds[j] != VarKind::Binary) {
        return 0;
    }
    let mut cap = bound;
    let mut w: Vec<f64> = row
        .vals
        .iter()
        .map(|&a| {
            let a = sign * a;
            if a < 0.0 {
                cap -= a;
            }
            a.abs()
        })
        .collect();
    w.sort_unstable_by(|a, b| b.total_cmp(a));
    let tol = EPS * cap.abs().max(1.0);
    let mut k = 0;
    while k + 1 < w.len() && w[k] + w[k + 1] > cap + tol {
        k += 1;
    }
    if k >= 1 { k + 1 } else { 0 }
}

struct UnionFind(Vec<usize>);

impl UnionFind {
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.0[ra] = rb;
        }
    }
}

/// Classifies `p`: row types, column / objective statistics, structure.
pub fn classify(p: &MipProblem, opts: ClassifyOptions) -> Classification {
    let kinds = var_kinds(p);
    let rows = rows_of(p);
    let (m, n) = (p.num_row, p.num_col);

    // --- rows -----------------------------------------------------------
    let mut row_class = Vec::with_capacity(m);
    let mut hist = [0usize; 17];
    for (i, row) in rows.iter().enumerate() {
        let c = classify_row(row, p.row_lower[i], p.row_upper[i], p, &kinds, opts);
        for t in c.sides.iter().flatten() {
            hist[t.index()] += 1;
        }
        row_class.push(c);
    }

    // --- columns --------------------------------------------------------
    let mut cols = ColumnStats::default();
    let mut col_count = vec![0usize; n];
    for r in &rows {
        for &j in &r.cols {
            col_count[j] += 1;
        }
    }
    for j in 0..n {
        match kinds[j] {
            VarKind::Binary => cols.binary += 1,
            VarKind::Integer => cols.integer += 1,
            VarKind::Continuous => cols.continuous += 1,
        }
        let (l, u) = (p.col_lower[j], p.col_upper[j]);
        match (is_inf(l), is_inf(u)) {
            (true, true) => cols.free += 1,
            (false, true) => cols.lower_only += 1,
            (true, false) => cols.upper_only += 1,
            (false, false) if eq(l, u) => cols.fixed += 1,
            (false, false) => cols.boxed += 1,
        }
        if kinds[j] == VarKind::Integer && (is_inf(l) || is_inf(u) || u - l > 1e4) {
            cols.wide_integer += 1;
        }
        if col_count[j] == 0 {
            cols.empty += 1;
        }
    }

    // --- objective --------------------------------------------------------
    let nz: Vec<(usize, f64)> = p
        .cost
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c != 0.0)
        .map(|(j, &c)| (j, c))
        .collect();
    let mut abs: Vec<f64> = nz.iter().map(|&(_, c)| c).collect();
    abs.sort_unstable_by(f64::total_cmp);
    let mut distinct = 0;
    let mut last: Option<f64> = None;
    for &c in &abs {
        if last.is_none_or(|l| !eq(l, c)) {
            distinct += 1;
            last = Some(c);
        }
    }
    let (pos, neg) = (
        nz.iter().filter(|&&(_, c)| c > 0.0).count(),
        nz.iter().filter(|&&(_, c)| c < 0.0).count(),
    );
    let mags = nz.iter().map(|&(_, c)| c.abs());
    let (cmin, cmax) = mags.fold((f64::INFINITY, 0.0f64), |(a, b), v| (a.min(v), b.max(v)));
    let objective = ObjectiveStats {
        nonzeros: nz.len(),
        zero_objective: nz.is_empty(),
        only_binaries: !nz.is_empty() && nz.iter().all(|&(j, _)| kinds[j] == VarKind::Binary),
        distinct_values: distinct,
        sign: match (pos > 0, neg > 0) {
            (false, false) => ObjSign::Zero,
            (true, false) => ObjSign::NonNegative,
            (false, true) => ObjSign::NonPositive,
            (true, true) => ObjSign::Mixed,
        },
        integral: nz.iter().all(|&(j, c)| p.integer[j] && integral(c)),
        dynamism: if nz.len() > 1 { cmax / cmin } else { 1.0 },
    };

    // --- structure --------------------------------------------------------
    let mut s = StructureStats {
        equality_rows: 0,
        ranged_rows: 0,
        vub_pairs: 0,
        vlb_pairs: 0,
        max_vb_ratio: 0.0,
        network_rows: 0,
        clique_rows: 0,
        distinct_cliques: 0,
        max_clique: 0,
        knapsack_cliques: 0,
        components: 0,
        largest_component: 0,
        density: 0.0,
        coef_min: 0.0,
        coef_max: 0.0,
        row_dynamism_max_log10: 0.0,
        row_dynamism_mean_log10: 0.0,
    };
    let mut vub: HashSet<(usize, usize)> = HashSet::new();
    let mut vlb: HashSet<(usize, usize)> = HashSet::new();
    let mut cliques: HashSet<Vec<Lit>> = HashSet::new();
    // Network scan state: per column, the sign (after row reflection) of the
    // entries accepted so far — at most one +1 and one −1.
    let mut net_plus = vec![false; n];
    let mut net_minus = vec![false; n];
    let mut uf = UnionFind((0..n).collect());
    let (mut amin, mut amax) = (f64::INFINITY, 0.0f64);
    let (mut dyn_sum, mut dyn_rows) = (0.0, 0usize);
    let nnz: usize = rows.iter().map(|r| r.cols.len()).sum();

    for (i, row) in rows.iter().enumerate() {
        let (lo, hi) = (p.row_lower[i], p.row_upper[i]);
        if !is_inf(lo) && !is_inf(hi) {
            if eq(lo, hi) {
                s.equality_rows += 1;
            } else {
                s.ranged_rows += 1;
            }
        }
        if row.cols.is_empty() {
            continue;
        }
        // coefficients
        let (rmin, rmax) = row.vals.iter().fold((f64::INFINITY, 0.0f64), |(a, b), v| {
            (a.min(v.abs()), b.max(v.abs()))
        });
        amin = amin.min(rmin);
        amax = amax.max(rmax);
        let d = (rmax / rmin).log10();
        s.row_dynamism_max_log10 = s.row_dynamism_max_log10.max(d);
        dyn_sum += d;
        dyn_rows += 1;
        // components
        for w in row.cols.windows(2) {
            uf.union(w[0], w[1]);
        }
        // variable upper / lower bounds: a·y + b·x with x binary, y not
        if row.cols.len() == 2 {
            let (y, x, a, b) = match (kinds[row.cols[0]], kinds[row.cols[1]]) {
                (VarKind::Binary, k) if k != VarKind::Binary => {
                    (row.cols[1], row.cols[0], row.vals[1], row.vals[0])
                }
                (k, VarKind::Binary) if k != VarKind::Binary => {
                    (row.cols[0], row.cols[1], row.vals[0], row.vals[1])
                }
                _ => (usize::MAX, 0, 0.0, 0.0),
            };
            if y != usize::MAX {
                // `≤` side bounds y above iff a > 0; `≥` side the opposite.
                if !is_inf(hi) {
                    (if a > 0.0 { &mut vub } else { &mut vlb }).insert((y, x));
                }
                if !is_inf(lo) {
                    (if a > 0.0 { &mut vlb } else { &mut vub }).insert((y, x));
                }
                s.max_vb_ratio = s.max_vb_ratio.max(b.abs() / a.abs());
            }
        }
        // cliques
        if let Some(mut lits) = clique_of(row, lo, hi, &kinds) {
            s.clique_rows += 1;
            s.max_clique = s.max_clique.max(lits.len());
            lits.sort_unstable();
            cliques.insert(lits);
        } else {
            let k = knapsack_clique_size(row, &kinds, 1.0, hi)
                .max(knapsack_clique_size(row, &kinds, -1.0, -lo));
            if k >= 2 {
                s.knapsack_cliques += 1;
            }
        }
        // embedded network: row-scanning addition with reflection
        if row.vals.iter().all(|&a| eq(a.abs(), 1.0)) {
            let mut reflect: Option<bool> = None;
            let mut ok = true;
            for (&j, &a) in row.cols.iter().zip(&row.vals) {
                let (fits_as_is, fits_reflected) = if a > 0.0 {
                    (!net_plus[j], !net_minus[j])
                } else {
                    (!net_minus[j], !net_plus[j])
                };
                let need = match (fits_as_is, fits_reflected) {
                    (true, true) => continue,
                    (true, false) => false,
                    (false, true) => true,
                    (false, false) => {
                        ok = false;
                        break;
                    }
                };
                if reflect.is_some_and(|r| r != need) {
                    ok = false;
                    break;
                }
                reflect = Some(need);
            }
            if ok {
                let r = reflect.unwrap_or(false);
                for (&j, &a) in row.cols.iter().zip(&row.vals) {
                    if (a > 0.0) != r {
                        net_plus[j] = true;
                    } else {
                        net_minus[j] = true;
                    }
                }
                s.network_rows += 1;
            }
        }
    }
    s.vub_pairs = vub.len();
    s.vlb_pairs = vlb.len();
    s.distinct_cliques = cliques.len();
    let mut comp_size = vec![0usize; n];
    for (j, &cnt) in col_count.iter().enumerate() {
        if cnt > 0 {
            let r = uf.find(j);
            comp_size[r] += 1;
        }
    }
    s.components = comp_size.iter().filter(|&&c| c > 0).count();
    s.largest_component = comp_size.iter().copied().max().unwrap_or(0);
    s.density = if m * n > 0 {
        nnz as f64 / (m as f64 * n as f64)
    } else {
        0.0
    };
    if nnz > 0 {
        s.coef_min = amin;
        s.coef_max = amax;
        s.row_dynamism_mean_log10 = dyn_sum / dyn_rows as f64;
    }

    Classification {
        num_rows: m,
        num_cols: n,
        nnz,
        row_class,
        row_histogram: hist,
        columns: cols,
        objective,
        structure: s,
    }
}

fn clamp01(v: f64) -> f64 {
    if v.is_nan() { 0.0 } else { v.clamp(0.0, 1.0) }
}

fn frac(a: usize, b: usize) -> f64 {
    if b == 0 { 0.0 } else { a as f64 / b as f64 }
}

impl Classification {
    /// Fixed-length feature vector, index-aligned with [`FEATURE_NAMES`]
    /// (normalisation documented there), version [`FEATURE_VERSION`].
    pub fn features(&self) -> [f64; 44] {
        let (m, n, nnz) = (self.num_rows, self.num_cols, self.nnz);
        let c = &self.columns;
        let o = &self.objective;
        let s = &self.structure;
        let lg = |v: f64| (1.0 + v).log10();
        let mut f = [0.0; 44];
        f[0] = lg(m as f64) / 7.0;
        f[1] = lg(n as f64) / 7.0;
        f[2] = lg(nnz as f64) / 8.0;
        f[3] = if nnz > 0 {
            (s.density.log10() + 6.0) / 6.0
        } else {
            0.0
        };
        f[4] = lg(frac(nnz, m)) / 4.0;
        f[5] = frac(c.binary, n);
        f[6] = frac(c.integer, n);
        f[7] = frac(c.continuous, n);
        f[8] = frac(c.free, n);
        f[9] = frac(c.fixed, n);
        f[10] = frac(c.boxed, n);
        f[11] = f64::from(u8::from(o.zero_objective));
        f[12] = frac(o.nonzeros, n);
        f[13] = f64::from(u8::from(o.only_binaries));
        f[14] = lg(o.distinct_values as f64) / 6.0;
        f[15] = match o.sign {
            ObjSign::Zero | ObjSign::NonPositive => 0.0,
            ObjSign::NonNegative => 1.0,
            ObjSign::Mixed => 0.5,
        };
        f[16] = f64::from(u8::from(o.integral));
        f[17] = o.dynamism.log10() / 12.0;
        let sides: usize = self.row_histogram.iter().sum();
        for (k, &h) in self.row_histogram.iter().enumerate() {
            f[18 + k] = frac(h, sides);
        }
        f[35] = frac(s.equality_rows, m);
        let nonempty = m - self.row_histogram[RowType::Empty.index()];
        f[36] = frac(s.network_rows, nonempty);
        f[37] = frac(s.clique_rows, m);
        f[38] = frac(s.vub_pairs + s.vlb_pairs, n);
        f[39] = (s.components.max(1) as f64).log10() / 7.0;
        f[40] = frac(s.largest_component, n);
        f[41] = if s.coef_min > 0.0 {
            (s.coef_max / s.coef_min).log10() / 12.0
        } else {
            0.0
        };
        f[42] = s.row_dynamism_max_log10 / 12.0;
        f[43] = lg(s.max_vb_ratio) / 12.0;
        f.map(clamp01)
    }

    /// Side count of one row type.
    pub fn count(&self, t: RowType) -> usize {
        self.row_histogram[t.index()]
    }

    /// Human-readable report.
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        let c = &self.columns;
        let o = &self.objective;
        let s = &self.structure;
        let _ = writeln!(
            out,
            "size:        {} rows, {} cols, {} nnz (density {:.3e})",
            self.num_rows, self.num_cols, self.nnz, s.density
        );
        let _ = writeln!(
            out,
            "columns:     {} binary, {} integer ({} wide), {} continuous; {} empty",
            c.binary, c.integer, c.wide_integer, c.continuous, c.empty
        );
        let _ = writeln!(
            out,
            "bounds:      {} free, {} lower-only, {} upper-only, {} boxed, {} fixed",
            c.free, c.lower_only, c.upper_only, c.boxed, c.fixed
        );
        let _ = writeln!(
            out,
            "objective:   {} nonzeros, sign {}, {} distinct values, only-binaries {}, \
             integral {}, dynamism {:.3e}",
            o.nonzeros,
            o.sign.name(),
            o.distinct_values,
            o.only_binaries,
            o.integral,
            o.dynamism
        );
        let _ = writeln!(out, "row types (MIPLIB 2017, side counts):");
        for t in RowType::ALL {
            let h = self.count(t);
            if h > 0 {
                let _ = writeln!(out, "  {:<20} {h}", t.name());
            }
        }
        let _ = writeln!(
            out,
            "rows:        {} equality, {} ranged",
            s.equality_rows, s.ranged_rows
        );
        let _ = writeln!(
            out,
            "var bounds:  {} VUB + {} VLB pairs, max |b/a| {:.3e}",
            s.vub_pairs, s.vlb_pairs, s.max_vb_ratio
        );
        let _ = writeln!(
            out,
            "cliques:     {} rows ({} distinct, max size {}), {} from knapsacks",
            s.clique_rows, s.distinct_cliques, s.max_clique, s.knapsack_cliques
        );
        let _ = writeln!(
            out,
            "network:     {} of {} rows",
            s.network_rows, self.num_rows
        );
        let _ = writeln!(
            out,
            "components:  {} (largest {} cols)",
            s.components, s.largest_component
        );
        let _ = writeln!(
            out,
            "coefficients: |a| in [{:.3e}, {:.3e}], row dynamism max 1e{:.2} mean 1e{:.2}",
            s.coef_min, s.coef_max, s.row_dynamism_max_log10, s.row_dynamism_mean_log10
        );
        out
    }

    /// Machine-readable report (one JSON object, stable key names).
    pub fn to_json(&self) -> String {
        let c = &self.columns;
        let o = &self.objective;
        let s = &self.structure;
        let num = |v: f64| {
            if v.is_finite() {
                format!("{v}")
            } else {
                "null".to_string()
            }
        };
        let mut j = String::from("{");
        let _ = write!(
            j,
            "\"rows\":{},\"cols\":{},\"nnz\":{},",
            self.num_rows, self.num_cols, self.nnz
        );
        let _ = write!(
            j,
            "\"columns\":{{\"binary\":{},\"integer\":{},\"continuous\":{},\"wide_integer\":{},\
             \"free\":{},\"lower_only\":{},\"upper_only\":{},\"boxed\":{},\"fixed\":{},\"empty\":{}}},",
            c.binary,
            c.integer,
            c.continuous,
            c.wide_integer,
            c.free,
            c.lower_only,
            c.upper_only,
            c.boxed,
            c.fixed,
            c.empty
        );
        let _ = write!(
            j,
            "\"objective\":{{\"nonzeros\":{},\"zero\":{},\"only_binaries\":{},\
             \"distinct_values\":{},\"sign\":\"{}\",\"integral\":{},\"dynamism\":{}}},",
            o.nonzeros,
            o.zero_objective,
            o.only_binaries,
            o.distinct_values,
            o.sign.name(),
            o.integral,
            num(o.dynamism)
        );
        j.push_str("\"row_types\":{");
        let items: Vec<String> = RowType::ALL
            .iter()
            .map(|&t| format!("\"{}\":{}", t.name(), self.count(t)))
            .collect();
        j.push_str(&items.join(","));
        j.push_str("},");
        let _ = write!(
            j,
            "\"structure\":{{\"equality_rows\":{},\"ranged_rows\":{},\"vub_pairs\":{},\
             \"vlb_pairs\":{},\"max_vb_ratio\":{},\"network_rows\":{},\"clique_rows\":{},\
             \"distinct_cliques\":{},\"max_clique\":{},\"knapsack_cliques\":{},\
             \"components\":{},\"largest_component\":{},\"density\":{},\"coef_min\":{},\
             \"coef_max\":{},\"row_dynamism_max_log10\":{},\"row_dynamism_mean_log10\":{}}},",
            s.equality_rows,
            s.ranged_rows,
            s.vub_pairs,
            s.vlb_pairs,
            num(s.max_vb_ratio),
            s.network_rows,
            s.clique_rows,
            s.distinct_cliques,
            s.max_clique,
            s.knapsack_cliques,
            s.components,
            s.largest_component,
            num(s.density),
            num(s.coef_min),
            num(s.coef_max),
            num(s.row_dynamism_max_log10),
            num(s.row_dynamism_mean_log10)
        );
        let _ = write!(j, "\"feature_version\":{FEATURE_VERSION},\"features\":{{");
        let items: Vec<String> = FEATURE_NAMES
            .iter()
            .zip(self.features())
            .map(|(k, v)| format!("\"{k}\":{}", num(v)))
            .collect();
        j.push_str(&items.join(","));
        j.push_str("}}");
        j
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_io::read_mps;

    /// Builds an MPS model: `rows` = (sense, rhs, optional range),
    /// `cols` = (name, kind 'B'/'I'/'C', lb, ub, cost, [(row, coef)]).
    fn mps(rows: &[(char, f64, Option<f64>)], cols: &[ColSpec]) -> String {
        let mut s = String::from("NAME T\nROWS\n N  OBJ\n");
        for (i, (sense, _, _)) in rows.iter().enumerate() {
            let _ = writeln!(s, " {sense}  R{i}");
        }
        s.push_str("COLUMNS\n");
        for (k, (name, kind, _, _, cost, entries)) in cols.iter().enumerate() {
            let int = *kind != 'C';
            if int {
                let _ = writeln!(s, "    M{k}  'MARKER'  'INTORG'");
            }
            let _ = writeln!(s, "    {name}  OBJ  {cost}");
            for (r, a) in entries.iter() {
                let _ = writeln!(s, "    {name}  R{r}  {a}");
            }
            if int {
                let _ = writeln!(s, "    E{k}  'MARKER'  'INTEND'");
            }
        }
        s.push_str("RHS\n");
        for (i, (_, b, _)) in rows.iter().enumerate() {
            let _ = writeln!(s, "    RHS  R{i}  {b}");
        }
        if rows.iter().any(|r| r.2.is_some()) {
            s.push_str("RANGES\n");
            for (i, (_, _, r)) in rows.iter().enumerate() {
                if let Some(r) = r {
                    let _ = writeln!(s, "    RNG  R{i}  {r}");
                }
            }
        }
        s.push_str("BOUNDS\n");
        for (name, _, lb, ub, _, _) in cols {
            if lb.is_infinite() && ub.is_infinite() {
                let _ = writeln!(s, " FR BND {name}");
                continue;
            }
            if lb.is_infinite() {
                let _ = writeln!(s, " MI BND {name}");
            } else {
                let _ = writeln!(s, " LO BND {name}  {lb}");
            }
            if ub.is_finite() {
                let _ = writeln!(s, " UP BND {name}  {ub}");
            }
        }
        s.push_str("ENDATA\n");
        s
    }

    fn run(text: &str, strict: bool) -> Classification {
        let p = MipProblem::from_model(&read_mps(text).expect("test MPS parses"));
        classify(
            &p,
            ClassifyOptions {
                miplib_strict: strict,
            },
        )
    }

    const INF: f64 = f64::INFINITY;

    /// `(name, kind 'B'/'I'/'C', lb, ub, cost, [(row, coef)])`.
    type ColSpec = (&'static str, char, f64, f64, f64, &'static [(usize, f64)]);

    /// Four binaries, two general integers (0..10), two continuous (0..10).
    fn cols_for(entries: &[&'static [(usize, f64)]; 8]) -> Vec<ColSpec> {
        let names = ["B0", "B1", "B2", "B3", "I0", "I1", "C0", "C1"];
        let kinds = ['B', 'B', 'B', 'B', 'I', 'I', 'C', 'C'];
        (0..8)
            .map(|k| {
                let ub = if kinds[k] == 'B' { 1.0 } else { 10.0 };
                (names[k], kinds[k], 0.0, ub, 1.0, entries[k])
            })
            .collect()
    }

    /// Every one of the 17 MIPLIB types, one row each, in table order.
    #[test]
    fn each_miplib_type_once() {
        // rows: 0 empty, 1 free, 2 singleton, 3 aggregation, 4 precedence,
        // 5 varbound, 6 set partition, 7 set packing, 8 set cover,
        // 9 cardinality, 10 invariant knapsack, 11 equation knapsack,
        // 12 bin packing, 13 knapsack, 14 integer knapsack, 15 mixed binary,
        // 16 general linear.
        let rows = [
            ('L', 5.0, None),   // 0: no entries
            ('G', -1e30, None), // 1: −∞ ≤ B0 + C0 (1e30 = infinite)
            ('L', 4.0, None),   // 2: C0 ≤ 4
            ('E', 3.0, None),   // 3: C0 + C1 = 3
            ('L', 2.0, None),   // 4: I0 − I1 ≤ 2
            ('L', 0.0, None),   // 5: C0 − 10 B0 ≤ 0
            ('E', 1.0, None),   // 6: B0 + B1 + B2 = 1
            ('L', 1.0, None),   // 7: B0 + B1 + B2 ≤ 1
            ('G', 1.0, None),   // 8: B0 + B1 + B2 ≥ 1
            ('E', 2.0, None),   // 9: B0 + B1 + B2 + B3 = 2
            ('L', 2.0, None),   // 10: B0 + B1 + B2 ≤ 2
            ('E', 5.0, None),   // 11: 2B0 + 3B1 + 4B2 = 5
            ('L', 0.0, None),   // 12: 3B0 + 4B1 − 7B3 ≤ 0 ⇔ 3B0 + 4B1 + 7B̄3 ≤ 7
            ('L', 6.0, None),   // 13: 3B0 + 4B1 + 5B2 ≤ 6
            ('L', 20.0, None),  // 14: 2I0 + 3I1 + B0 ≤ 20
            ('L', 5.0, None),   // 15: B0 + B1 + 2.5 C0 ≤ 5
            ('L', 5.0, None),   // 16: I0 + C0 − B0 ≤ 5
        ];
        static E: [&[(usize, f64)]; 8] = [
            // B0
            &[
                (1, 1.0),
                (5, -10.0),
                (6, 1.0),
                (7, 1.0),
                (8, 1.0),
                (9, 1.0),
                (10, 1.0),
                (11, 2.0),
                (12, 3.0),
                (13, 3.0),
                (14, 1.0),
                (15, 1.0),
                (16, -1.0),
            ],
            // B1
            &[
                (6, 1.0),
                (7, 1.0),
                (8, 1.0),
                (9, 1.0),
                (10, 1.0),
                (11, 3.0),
                (12, 4.0),
                (13, 4.0),
                (15, 1.0),
            ],
            // B2
            &[
                (6, 1.0),
                (7, 1.0),
                (8, 1.0),
                (9, 1.0),
                (10, 1.0),
                (11, 4.0),
                (13, 5.0),
            ],
            // B3
            &[(9, 1.0), (12, -7.0)],
            // I0
            &[(4, 1.0), (14, 2.0), (16, 1.0)],
            // I1
            &[(4, -1.0), (14, 3.0)],
            // C0
            &[(1, 1.0), (2, 1.0), (3, 1.0), (5, 1.0), (15, 2.5), (16, 1.0)],
            // C1
            &[(3, 1.0)],
        ];
        let c = run(&mps(&rows, &cols_for(&E)), true);
        let got: Vec<RowType> = c.row_class.iter().map(|r| r.primary()).collect();
        assert_eq!(got, RowType::ALL.to_vec());
        // Exactly one side each: histogram is all ones.
        assert_eq!(c.row_histogram, [1; 17], "{:?}", c.row_histogram);
    }

    #[test]
    fn ranged_rows_count_both_sides() {
        // 0: 1 ≤ B0+B1+B2 ≤ 1 handled as partition; use ranged 1 ≤ Σ ≤ 2:
        //    ≤ side invariant knapsack (b=2), ≥ side set covering.
        // 1: ranged singleton 1 ≤ C0 ≤ 3 counts 2.
        // 2: ranged precedence −1 ≤ C0 − C1 ≤ 1 counts 2.
        // 3: ranged general 1 ≤ I0 + C0 ≤ 4 counts 2.
        let rows = [
            ('G', 1.0, Some(1.0)),
            ('G', 1.0, Some(2.0)),
            ('G', -1.0, Some(2.0)),
            ('G', 1.0, Some(3.0)),
        ];
        static E: [&[(usize, f64)]; 8] = [
            &[(0, 1.0)],
            &[(0, 1.0)],
            &[(0, 1.0)],
            &[],
            &[(3, 1.0)],
            &[],
            &[(1, 1.0), (2, 1.0), (3, 1.0)],
            &[(2, -1.0)],
        ];
        let c = run(&mps(&rows, &cols_for(&E)), false);
        assert_eq!(
            c.row_class[0].sides,
            [Some(RowType::InvariantKnapsack), Some(RowType::SetCovering)]
        );
        assert_eq!(c.count(RowType::Singleton), 2);
        assert_eq!(c.count(RowType::Precedence), 2);
        assert_eq!(c.count(RowType::GeneralLinear), 2);
        assert_eq!(c.structure.ranged_rows, 4);
    }

    #[test]
    fn complemented_binaries_and_scaling() {
        // 0: 2B0 + 2B1 − 2B2 ≤ 0  ⇔  B0 + B1 + (1−B2) ≤ 1: set packing.
        // 1: B0 + B1 − B2 = 0      ⇔  B0 + B1 + (1−B2) = 1: set partitioning.
        // 2: −B0 − B1 − B2 ≤ −1    ⇔  Σ ≥ 1 seen from the ≤ side: b = 2 → inv. knapsack.
        let rows = [('L', 0.0, None), ('E', 0.0, None), ('L', -1.0, None)];
        static E: [&[(usize, f64)]; 8] = [
            &[(0, 2.0), (1, 1.0), (2, -1.0)],
            &[(0, 2.0), (1, 1.0), (2, -1.0)],
            &[(0, -2.0), (1, -1.0), (2, -1.0)],
            &[],
            &[],
            &[],
            &[],
            &[],
        ];
        let c = run(&mps(&rows, &cols_for(&E)), true);
        let got: Vec<RowType> = c.row_class.iter().map(|r| r.primary()).collect();
        assert_eq!(
            got,
            vec![
                RowType::SetPacking,
                RowType::SetPartitioning,
                RowType::InvariantKnapsack
            ]
        );
        // rows 0/1 are the same clique {B0, B1, B̄2}; row 2 is
        // B̄0 + B̄1 + B̄2 ≤ 2 — not a clique (capacity 2).
        assert_eq!(c.structure.clique_rows, 2);
        assert_eq!(c.structure.distinct_cliques, 1);
    }

    #[test]
    fn geq_knapsack_strict_vs_default() {
        // 3B0 + 4B1 + 5B2 ≥ 6 : MIPLIB-strict → mixed binary; default → knapsack.
        let rows = [('G', 6.0, None)];
        static E: [&[(usize, f64)]; 8] = [
            &[(0, 3.0)],
            &[(0, 4.0)],
            &[(0, 5.0)],
            &[],
            &[],
            &[],
            &[],
            &[],
        ];
        let text = mps(&rows, &cols_for(&E));
        assert_eq!(
            run(&text, true).row_class[0].primary(),
            RowType::MixedBinary
        );
        assert_eq!(run(&text, false).row_class[0].primary(), RowType::Knapsack);
    }

    #[test]
    fn integer_equality_counts_twice() {
        // 2I0 + 3I1 + B0 = 12 (SCIP: int knapsack on the ≤ side, general
        // linear on the remaining ≥ side). Two terms would be an aggregation.
        let rows = [('E', 12.0, None)];
        static E: [&[(usize, f64)]; 8] = [
            &[(0, 1.0)],
            &[],
            &[],
            &[],
            &[(0, 2.0)],
            &[(0, 3.0)],
            &[],
            &[],
        ];
        let c = run(&mps(&rows, &cols_for(&E)), true);
        assert_eq!(
            c.row_class[0].sides,
            [Some(RowType::IntegerKnapsack), Some(RowType::GeneralLinear)]
        );
    }

    #[test]
    fn transportation_is_a_network_with_one_component() {
        // 2×2 transportation: supply rows Σ_j x_ij ≤ s_i, demand rows
        // Σ_i x_ij ≥ d_j — every column has one +1 in a supply and one +1 in
        // a demand row: network after reflecting the demand rows.
        let rows = [
            ('L', 5.0, None),
            ('L', 5.0, None),
            ('G', 4.0, None),
            ('G', 4.0, None),
        ];
        let cols: Vec<ColSpec> = vec![
            ("X00", 'C', 0.0, INF, 1.0, &[(0, 1.0), (2, 1.0)]),
            ("X01", 'C', 0.0, INF, 2.0, &[(0, 1.0), (3, 1.0)]),
            ("X10", 'C', 0.0, INF, 3.0, &[(1, 1.0), (2, 1.0)]),
            ("X11", 'C', 0.0, INF, 1.0, &[(1, 1.0), (3, 1.0)]),
        ];
        let c = run(&mps(&rows, &cols), false);
        assert_eq!(c.structure.network_rows, 4);
        assert_eq!(c.structure.components, 1);
        assert_eq!(c.structure.largest_component, 4);
        assert_eq!(c.columns.continuous, 4);
        assert_eq!(c.columns.lower_only, 4);
        assert_eq!(c.objective.distinct_values, 3);
        assert_eq!(c.objective.sign, ObjSign::NonNegative);
        assert!(!c.objective.integral, "continuous columns");
    }

    #[test]
    fn non_network_row_is_rejected() {
        // Three rows all containing X0 with +1: the third cannot be added.
        let rows = [('L', 1.0, None), ('L', 1.0, None), ('L', 1.0, None)];
        let cols: Vec<ColSpec> = vec![
            ("X0", 'C', 0.0, INF, 0.0, &[(0, 1.0), (1, 1.0), (2, 1.0)]),
            ("X1", 'C', 0.0, INF, 0.0, &[(0, 1.0)]),
            ("X2", 'C', 0.0, INF, 0.0, &[(2, 1.0)]),
        ];
        let c = run(&mps(&rows, &cols), false);
        assert_eq!(c.structure.network_rows, 2);
        assert!(c.objective.zero_objective);
        assert_eq!(c.objective.sign, ObjSign::Zero);
    }

    #[test]
    fn two_term_continuous_row_strict_vs_default() {
        // C0 + 2 C1 ≤ 4: no binary ⇒ the table's variable bound does not
        // apply (default: mixed binary); published MIPLIB 2017 data files it
        // as a variable bound (strict).
        let rows = [('L', 4.0, None)];
        static E: [&[(usize, f64)]; 8] = [&[], &[], &[], &[], &[], &[], &[(0, 1.0)], &[(0, 2.0)]];
        let text = mps(&rows, &cols_for(&E));
        assert_eq!(
            run(&text, false).row_class[0].primary(),
            RowType::MixedBinary
        );
        assert_eq!(
            run(&text, true).row_class[0].primary(),
            RowType::VariableBound
        );
    }

    #[test]
    fn variable_bounds_and_components() {
        // C0 ≤ 10 B0 (VUB, big-M 10), C1 ≥ 2 B1 (VLB), separate blocks.
        let rows = [('L', 0.0, None), ('G', 0.0, None)];
        static E: [&[(usize, f64)]; 8] = [
            &[(0, -10.0)],
            &[(1, -2.0)],
            &[],
            &[],
            &[],
            &[],
            &[(0, 1.0)],
            &[(1, 1.0)],
        ];
        let c = run(&mps(&rows, &cols_for(&E)), false);
        assert_eq!(c.count(RowType::VariableBound), 2);
        assert_eq!(c.structure.vub_pairs, 1);
        assert_eq!(c.structure.vlb_pairs, 1);
        assert_eq!(c.structure.max_vb_ratio, 10.0);
        assert_eq!(c.structure.components, 2);
        assert_eq!(c.columns.empty, 4);
        assert_eq!(c.structure.coef_min, 1.0);
        assert_eq!(c.structure.coef_max, 10.0);
    }

    #[test]
    fn knapsack_clique_extraction() {
        // 5B0 + 4B1 + 4B2 + 1B3 ≤ 7: the three heaviest pairwise exceed 7 → clique of 3.
        let rows = [('L', 7.0, None)];
        static E: [&[(usize, f64)]; 8] = [
            &[(0, 5.0)],
            &[(0, 4.0)],
            &[(0, 4.0)],
            &[(0, 1.0)],
            &[],
            &[],
            &[],
            &[],
        ];
        let c = run(&mps(&rows, &cols_for(&E)), false);
        assert_eq!(c.row_class[0].primary(), RowType::Knapsack);
        assert_eq!(c.structure.knapsack_cliques, 1);
        assert_eq!(c.structure.clique_rows, 0);
        let p = MipProblem::from_model(&read_mps(&mps(&rows, &cols_for(&E))).unwrap());
        let rows = rows_of(&p);
        assert_eq!(knapsack_clique_size(&rows[0], &var_kinds(&p), 1.0, 7.0), 3);
    }

    #[test]
    fn edge_rows_are_variable_bounds_but_cliques() {
        // x_u + x_v ≤ 1 with binaries: the MIPLIB table files 2-term rows
        // with a binary as variable bound; the clique detector still sees it.
        let rows = [('L', 1.0, None), ('L', 1.0, None)];
        static E: [&[(usize, f64)]; 8] = [
            &[(0, 1.0)],
            &[(0, 1.0), (1, 1.0)],
            &[(1, 1.0)],
            &[],
            &[],
            &[],
            &[],
            &[],
        ];
        let c = run(&mps(&rows, &cols_for(&E)), true);
        assert_eq!(c.count(RowType::VariableBound), 2);
        assert_eq!(c.structure.clique_rows, 2);
        assert_eq!(c.structure.max_clique, 2);
        assert_eq!(c.structure.vub_pairs, 0, "binary-binary rows are not VUBs");
    }

    #[test]
    fn objective_stats_on_binaries() {
        // max 3B0 + 3B1 + 5B2 → min-form costs all negative, integral.
        let text = "\
NAME T
OBJSENSE
 MAX
ROWS
 N  OBJ
 L  R0
COLUMNS
    B0  OBJ  3  R0  1
    B1  OBJ  3  R0  1
    B2  OBJ  5  R0  1
RHS
    RHS  R0  1
BOUNDS
 BV BND B0
 BV BND B1
 BV BND B2
ENDATA
";
        let c = run(text, false);
        assert_eq!(c.columns.binary, 3);
        assert_eq!(c.objective.sign, ObjSign::NonPositive);
        assert!(c.objective.only_binaries);
        assert!(c.objective.integral);
        assert_eq!(c.objective.distinct_values, 2);
        assert_eq!(c.objective.dynamism, 5.0 / 3.0);
        assert_eq!(c.row_class[0].primary(), RowType::SetPacking);
    }

    #[test]
    fn feature_vector_is_bounded_and_named() {
        let rows = [('L', 7.0, None), ('G', 1.0, None)];
        static E: [&[(usize, f64)]; 8] = [
            &[(0, 5.0), (1, 1.0)],
            &[(0, 4.0), (1, 1.0)],
            &[(0, 4.0)],
            &[],
            &[(1, 1.0)],
            &[],
            &[(1, 1e7)],
            &[],
        ];
        let c = run(&mps(&rows, &cols_for(&E)), false);
        let f = c.features();
        assert_eq!(f.len(), FEATURE_NAMES.len());
        assert!(f.iter().all(|v| (0.0..=1.0).contains(v)), "{f:?}");
        let rows_sum: f64 = f[18..35].iter().sum();
        assert!((rows_sum - 1.0).abs() < 1e-12);
        // JSON names every feature exactly once.
        let j = c.to_json();
        for name in FEATURE_NAMES {
            assert_eq!(j.matches(&format!("\"{name}\":")).count(), 1, "{name}");
        }
        assert!(j.starts_with('{') && j.ends_with('}'));
    }
}
