//! `rusolve lint <model.mps> [--format text|json] [--rule ID]... [--min-severity S]
//! [--max-items N] [--with-solve]`, the static model lint.
//!
//! One pass over the port's `HighsLp` (loaded by `mip_highs::read_mip`, so every finding
//! names its rows and columns) finds the modelling problems that make solvers slow or
//! numerically fragile. Nothing here reads or changes solver state: no rule feeds the search,
//! so the default search cannot move. The only rule that solves is `condition-estimate`, and
//! only under `--with-solve`: it solves the LP relaxation on a private `Highs` and reads the
//! basis condition estimate the port computes (`HEkk::computeBasisCondition`).
//!
//! ## Rules and their sources
//!
//! | id | what | source |
//! |---|---|---|
//! | `big-m` | integer coefficient above what the row's activity bound implies | Achterberg et al. 2020 (see [`big_m`]); `HPresolve.cpp:3310-3370` |
//! | `coef-range` `rhs-range` `obj-range` `bound-range` | smallest/largest nonzero magnitude and their ratio | `HighsSolve.cpp:518-600`, `HConst.h:42-46`, common solver numerics guidance |
//! | `scaling` | row/column geometric scale far from the median | `HighsLpUtils.cpp:983,997` |
//! | `near-parallel-rows` | rows equal up to a scalar | `HPresolveParallelRowsAndCols.cpp:15-60` |
//! | `empty-row` `empty-col` `free-row` `free-col-no-cost` `fixed-col` | structure | `HPresolve` empty/free/fixed rules |
//! | `fractional-int-bound` | integer column with a non-integral bound | `HighsLpUtils.cpp` assessIntegrality (rounds bounds) |
//! | `condition-estimate` | basis condition number after an LP solve | `HEkk.cpp:3649` `computeBasisCondition`, `Highs.cpp:2699` `getKappa` |
//!
//! ## Output and exit codes
//!
//! Text by default; `--format json` prints one `rusolve.lint/1` object
//! (`docs/schema/rusolve-lint-1.schema.json`). Exit codes follow the exit-code table (see `mip_highs.rs`):
//! 0 the lint ran (findings do not change the code), 2 usage error, 3 input error,
//! 4 internal error.

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::{Map, Value, json};
use std::collections::HashMap;
use std::process::ExitCode;

// ---------------------------------------------------------------------------------------
// Thresholds. Each names its source.
// ---------------------------------------------------------------------------------------

/// HiGHS `kExcessivelySmallObjectiveCoefficient` (`lp_data/HConst.h:42`): HighsSolve.cpp:565
/// warns about costs below it.
const SMALL_OBJ: f64 = 1e-4;
/// HiGHS `kExcessivelyLargeObjectiveCoefficient` (`HConst.h:43`), warned about at
/// `HighsSolve.cpp:570`.
const LARGE_OBJ: f64 = 1e6;
/// HiGHS `kExcessivelySmallBoundValue` (`HConst.h:44`); also used for row bounds
/// (`HighsSolve.cpp:595`).
const SMALL_BOUND: f64 = 1e-4;
/// HiGHS `kExcessivelyLargeBoundValue` (`HConst.h:45`): continuous bounds and row bounds.
const LARGE_BOUND: f64 = 1e6;
/// HiGHS `kExcessivelyLargeIntegerBoundValue` (`HConst.h:46`): integer-column bounds.
const LARGE_INT_BOUND: f64 = 1e5;
/// HiGHS option `large_matrix_value` default (`HighsOptions.h:829`): `assessLp` rejects a
/// matrix entry whose magnitude reaches it.
const LARGE_MATRIX_VALUE: f64 = 1e15;
/// HiGHS option `small_matrix_value` default (`HighsOptions.h:822`): `assessLp` removes
/// entries below it, and parallel-row detection uses it as its value tolerance
/// (`HPresolveParallelRowsAndCols.cpp:50`).
const SMALL_MATRIX_VALUE: f64 = 1e-9;
/// Ratio of largest to smallest nonzero magnitude: common solver numerics guidance advises matrix
/// coefficients in about [1e-3, 1e6], a ratio of 1e9; beyond it is a
/// warning. One thousand times less (1e6) is reported as info.
const RANGE_RATIO_WARN: f64 = 1e9;
const RANGE_RATIO_INFO: f64 = 1e6;
/// Scale of a row/column (`1/sqrt(min*max)`, HiGHS `equilibrationScaleMatrix`,
/// `HighsLpUtils.cpp:983,997`) relative to the median: 4 orders of magnitude is info, 6 a warning
/// (any scale factor F implies a coefficient range of at least F, so these sit below the range levels).
const SCALE_RATIO_WARN: f64 = 1e6;
const SCALE_RATIO_INFO: f64 = 1e4;
/// A coefficient is reported as big-M when it exceeds the tightest valid value by this
/// factor. HiGHS strengthens any excess above `primal_feastol` (`HPresolve.cpp:3331`); 10 keeps
/// the report to the coefficients that visibly loosen the relaxation.
const BIGM_MIN_RATIO: f64 = 10.0;
/// HiGHS `primal_feasibility_tolerance` default, the slack in the strengthening test
/// (`HPresolve.cpp:3331`).
const PRIMAL_FEASTOL: f64 = 1e-7;
/// Smallest implied-slack the tightening divides by; below it the row is treated as redundant.
const BIGM_MIN_SLACK: f64 = 1e-9;
/// Two normalised rows are parallel when every normalised coefficient differs by at most this
/// (HiGHS `small_matrix_value` tolerance, `HPresolveParallelRowsAndCols.cpp:50`).
const PARALLEL_TOL: f64 = SMALL_MATRIX_VALUE;
/// Bound tolerance when deciding that two parallel rows' ranges do not intersect (HiGHS
/// `primal_feasibility_tolerance`, relative to the bound size).
const PARALLEL_FEASTOL: f64 = 1e-6;
/// Basis condition number levels, after common solver numerics guidance: up to about 1e6 is
/// stable, past 1e9 several of the 16
/// digits of double precision are lost; 1/eps = 4.5e15, so 1e14 leaves no usable digit.
const KAPPA_INFO: f64 = 1e6;
const KAPPA_WARN: f64 = 1e9;
const KAPPA_ERROR: f64 = 1e14;
/// Items listed per finding group in the output unless `--max-items` says otherwise.
const DEFAULT_MAX_ITEMS: usize = 10;

// ---------------------------------------------------------------------------------------
// Rules and severities.
// ---------------------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Severity {
    Info,
    Warning,
    Error,
}

impl Severity {
    fn name(self) -> &'static str {
        match self {
            Severity::Info => "info",
            Severity::Warning => "warning",
            Severity::Error => "error",
        }
    }
    fn parse(s: &str) -> Option<Self> {
        match s {
            "info" => Some(Severity::Info),
            "warning" => Some(Severity::Warning),
            "error" => Some(Severity::Error),
            _ => None,
        }
    }
}

/// `(id, one-line meaning, fix)`; the order is the report order.
const RULES: &[(&str, &str, &str)] = &[
    (
        "big-m",
        "integer coefficient far above the row's implied activity bound (weak big-M)",
        "replace the coefficient by the tightest valid value shown and adjust the row bound as shown",
    ),
    (
        "coef-range",
        "constraint-matrix magnitudes span a wide range",
        "rescale rows and columns (change units) so the nonzeros fall in about [1e-3, 1e6]; drop entries below 1e-9",
    ),
    (
        "rhs-range",
        "constraint bound magnitudes span a wide range or are extreme",
        "change units so right-hand sides stay below 1e6 and above 1e-4",
    ),
    (
        "obj-range",
        "objective coefficient magnitudes span a wide range or are extreme",
        "rescale the objective (divide by a common factor, change units) so costs fall in about [1e-4, 1e6]",
    ),
    (
        "bound-range",
        "variable bound magnitudes span a wide range or are extreme",
        "tighten huge bounds to the largest value the variable can really take, or change its units",
    ),
    (
        "scaling",
        "row or column scale far from the median of the model",
        "multiply the row (or substitute x' = s*x for the column) by the scale factor shown",
    ),
    (
        "near-parallel-rows",
        "two rows are the same up to a scalar",
        "merge them into one ranged row (intersection of the bounds shown); if the ranges do not intersect the model is infeasible",
    ),
    (
        "empty-row",
        "row with no nonzeros",
        "delete the row (it constrains nothing); if its bounds exclude zero the model is infeasible",
    ),
    (
        "empty-col",
        "column with no nonzeros in the matrix",
        "delete the column; if it has a cost its value is set by its bounds alone",
    ),
    (
        "free-row",
        "row with neither a lower nor an upper bound",
        "delete the row: it can never bind",
    ),
    (
        "free-col-no-cost",
        "free continuous column with zero cost in a single row: the row is always satisfiable",
        "delete the column and its row (a free column singleton)",
    ),
    (
        "fixed-col",
        "column with equal lower and upper bound",
        "substitute the value into the rows and delete the column",
    ),
    (
        "fractional-int-bound",
        "integer column with a non-integral bound",
        "round the bound as shown (lower up, upper down); if the rounded bounds cross the model is infeasible",
    ),
    (
        "condition-estimate",
        "estimated condition number of the final LP-relaxation basis",
        "rescale (see scaling and coef-range), tighten big-M values, remove near-parallel rows",
    ),
];

fn rule_index(id: &str) -> Option<usize> {
    RULES.iter().position(|r| r.0 == id)
}

/// One located problem: names, evidence numbers, and a rank (larger = report first).
struct Item {
    rows: Vec<String>,
    cols: Vec<String>,
    evidence: Map<String, Value>,
    fix: Option<String>,
    rank: f64,
}

struct Group {
    rule: usize,
    severity: Severity,
    message: String,
    items: Vec<Item>,
}

/// Findings keyed by (rule, severity).
#[derive(Default)]
struct Findings {
    groups: Vec<Group>,
}

impl Findings {
    fn push(&mut self, rule: &str, severity: Severity, message: Option<String>, item: Item) {
        let rule = rule_index(rule).expect("known rule");
        let message = message.unwrap_or_else(|| RULES[rule].1.to_string());
        let g = match self
            .groups
            .iter_mut()
            .find(|g| g.rule == rule && g.severity == severity && g.message == message)
        {
            Some(g) => g,
            None => {
                self.groups.push(Group {
                    rule,
                    severity,
                    message,
                    items: Vec::new(),
                });
                self.groups.last_mut().unwrap()
            }
        };
        g.items.push(item);
    }
}

/// Terms shown by [`View::render_row`] before it elides the rest.
const RENDER_TERMS: usize = 8;

/// A number as a person writes it: `5`, `0.25`, `1e6`, `-inf`.
fn short(x: f64) -> String {
    if !x.is_finite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let m = x.abs();
    if m != 0.0 && !(1e-3..1e5).contains(&m) {
        format!("{x:e}")
    } else {
        format!("{}", x + 0.0)
    }
}

fn num(x: f64) -> Value {
    if x.is_finite() {
        json!(x + 0.0) // + 0.0 turns -0.0 into 0.0
    } else if x > 0.0 {
        json!("inf")
    } else if x < 0.0 {
        json!("-inf")
    } else {
        json!("nan")
    }
}

fn ev(pairs: &[(&str, Value)]) -> Map<String, Value> {
    pairs
        .iter()
        .map(|(k, v)| (k.to_string(), v.clone()))
        .collect()
}

// ---------------------------------------------------------------------------------------
// Model view.
// ---------------------------------------------------------------------------------------

struct View<'a> {
    lp: &'a HighsLp,
    /// Row-wise nonzeros `(col, value)`, sorted by column, explicit zeros dropped.
    rows: Vec<Vec<(usize, f64)>>,
    col_nnz: Vec<usize>,
    nnz: usize,
}

impl<'a> View<'a> {
    #[allow(clippy::needless_range_loop)] // j indexes both the matrix starts and `col_nnz`
    fn new(lp: &'a HighsLp) -> Self {
        let m = &lp.a_matrix_;
        let nr = lp.num_row_ as usize;
        let nc = lp.num_col_ as usize;
        let mut rows = vec![Vec::new(); nr];
        let mut col_nnz = vec![0usize; nc];
        let mut nnz = 0;
        // read_mip hands over a column-wise matrix; walk it column by column so every row
        // vector comes out sorted by column.
        for j in 0..nc {
            for k in m.start_[j] as usize..m.start_[j + 1] as usize {
                let v = m.value_[k];
                if v == 0.0 {
                    continue;
                }
                rows[m.index_[k] as usize].push((j, v));
                col_nnz[j] += 1;
                nnz += 1;
            }
        }
        View {
            lp,
            rows,
            col_nnz,
            nnz,
        }
    }
    fn rname(&self, i: usize) -> String {
        self.lp
            .row_names_
            .get(i)
            .cloned()
            .unwrap_or_else(|| format!("R{i}"))
    }
    fn cname(&self, j: usize) -> String {
        self.lp
            .col_names_
            .get(j)
            .cloned()
            .unwrap_or_else(|| format!("C{j}"))
    }
    /// Row `i` written out as the model states it (`c1: x - 1e6 b <= 0`), at most
    /// [`RENDER_TERMS`] terms, for a diagnostic's location line (docs/design/diagnostics.md).
    fn render_row(&self, i: usize) -> String {
        let mut e = String::new();
        for (k, &(j, a)) in self.rows[i].iter().enumerate() {
            if k == RENDER_TERMS {
                e.push_str(&format!(" + ... ({} more)", self.rows[i].len() - k));
                break;
            }
            let (sign, m) = if a < 0.0 { ("-", -a) } else { ("+", a) };
            let coef = if m == 1.0 {
                String::new()
            } else {
                format!("{} ", short(m))
            };
            if k == 0 {
                e.push_str(&format!(
                    "{}{coef}{}",
                    if a < 0.0 { "-" } else { "" },
                    self.cname(j)
                ));
            } else {
                e.push_str(&format!(" {sign} {coef}{}", self.cname(j)));
            }
        }
        let (lo, up) = (self.lp.row_lower_[i], self.lp.row_upper_[i]);
        let body = if lo == up {
            format!("{e} = {}", short(up))
        } else if lo.is_finite() && up.is_finite() {
            format!("{} <= {e} <= {}", short(lo), short(up))
        } else if up.is_finite() {
            format!("{e} <= {}", short(up))
        } else if lo.is_finite() {
            format!("{e} >= {}", short(lo))
        } else {
            format!("{e} (free)")
        };
        format!("{}: {body}", self.rname(i))
    }
    /// Column `j`'s bounds and type (`0 <= b <= 1, integer`), for a diagnostic's location line.
    fn render_col(&self, j: usize) -> String {
        let (lo, up) = (self.lp.col_lower_[j], self.lp.col_upper_[j]);
        let kind = if self.is_int(j) { ", integer" } else { "" };
        format!("{} <= {} <= {}{kind}", short(lo), self.cname(j), short(up))
    }
    fn is_int(&self, j: usize) -> bool {
        matches!(
            self.lp.integrality_.get(j),
            Some(HighsVarType::kInteger) | Some(HighsVarType::kSemiInteger)
        )
    }
    fn num_int(&self) -> usize {
        (0..self.lp.num_col_ as usize)
            .filter(|&j| self.is_int(j))
            .count()
    }
}

// ---------------------------------------------------------------------------------------
// Rules.
// ---------------------------------------------------------------------------------------

/// `big-m`: coefficient tightening of an integer column in an inequality.
///
/// Achterberg, Bixby, Gu, Rothberg, Weninger, "Presolve reductions in mixed integer
/// programming", INFORMS J. Computing 32(2), 2020, section "coefficient tightening": for a row
/// `a x <= b` with maximal activity `M = sum max(a_k l_k, a_k u_k)` and slack `d = M - b > 0`, an
/// integer column with `a_j > d` can have `a_j` lowered to `d` and `b` lowered by
/// `(a_j - d) u_j` (and symmetrically `a_j < -d` raised to `-d`, `b` raised by `(-d - a_j) l_j`)
/// without cutting off any integer-feasible point, and the LP relaxation only gets tighter. HiGHS
/// does the same in `HPresolve::dominatedColumns`' row loop, lambda `strengthenCoefs`
/// (`HPresolve.cpp:3310-3349`, applied to `<=` rows at :3351-3361 and to `>=` rows by negation at
/// :3363-3373); like it, only rows with a single finite side are used and `M` must be finite.
/// HiGHS takes the other terms from implied row-column bounds; this uses the column bounds,
/// tightened by the bounds that singleton rows imply (`HPresolve::singletonRow`,
/// `HPresolve.cpp:3205-3325`: a row `lo <= a x_j <= up` becomes `lo/a <= x_j <= up/a` for `a > 0`
/// and `up/a <= x_j <= lo/a` for `a < 0`), which is still the weaker (valid) variant. Agents write
/// `x <= 5` as a constraint row as often as a bound, so both must count.
fn big_m(v: &View, out: &mut Findings) {
    let lp = v.lp;
    let mut col_lo = lp.col_lower_.clone();
    let mut col_up = lp.col_upper_.clone();
    for (i, row) in v.rows.iter().enumerate() {
        if let [(j, a)] = row[..] {
            let (r_lo, r_up) = (lp.row_lower_[i], lp.row_upper_[i]);
            let (lo, up) = if a > 0.0 {
                (r_lo / a, r_up / a)
            } else {
                (r_up / a, r_lo / a)
            };
            if lo > col_lo[j] {
                col_lo[j] = lo;
            }
            if up < col_up[j] {
                col_up[j] = up;
            }
        }
    }
    for (i, row) in v.rows.iter().enumerate() {
        if row.len() < 2 {
            continue;
        }
        let (lo, up) = (lp.row_lower_[i], lp.row_upper_[i]);
        let dir = if lo == f64::NEG_INFINITY && up.is_finite() {
            1.0
        } else if up == f64::INFINITY && lo.is_finite() {
            -1.0
        } else {
            continue;
        };
        let b = if dir > 0.0 { up } else { -lo };
        let mut max_act = 0.0;
        let mut finite = true;
        for &(j, a) in row {
            let val = dir * a;
            let bound = if val > 0.0 { col_up[j] } else { col_lo[j] };
            if !bound.is_finite() {
                finite = false;
                break;
            }
            max_act += val * bound;
        }
        if !finite {
            continue;
        }
        let d = max_act - b;
        if d <= BIGM_MIN_SLACK {
            continue;
        }
        for &(j, a) in row {
            let (l, u) = (col_lo[j], col_up[j]);
            if !v.is_int(j) || l >= u {
                continue;
            }
            let val = dir * a;
            let (new_val, new_b) = if val > d + PRIMAL_FEASTOL {
                (d, b - (val - d) * u)
            } else if val < -d - PRIMAL_FEASTOL {
                (-d, b + (-d - val) * l)
            } else {
                continue;
            };
            let ratio = val.abs() / d;
            if ratio < BIGM_MIN_RATIO {
                continue;
            }
            // report in the row's own orientation
            let (tight_coef, tight_rhs) = (dir * new_val, dir * new_b);
            let side = if dir > 0.0 { "upper" } else { "lower" };
            out.push(
                "big-m",
                Severity::Warning,
                None,
                Item {
                    rows: vec![v.rname(i)],
                    cols: vec![v.cname(j)],
                    evidence: ev(&[
                        ("coefficient", num(a)),
                        ("tightest_coefficient", num(tight_coef)),
                        ("ratio", num(ratio)),
                        ("max_activity_slack", num(d)),
                        (
                            if dir > 0.0 { "row_upper" } else { "row_lower" },
                            num(if dir > 0.0 { up } else { lo }),
                        ),
                        (
                            if dir > 0.0 {
                                "tightest_row_upper"
                            } else {
                                "tightest_row_lower"
                            },
                            num(tight_rhs),
                        ),
                        ("column_lower", num(l)),
                        ("column_upper", num(u)),
                    ]),
                    fix: Some(format!(
                        "set the coefficient of {} in {} to {} and the {} bound to {} (was {} and {})",
                        v.cname(j),
                        v.rname(i),
                        tight_coef,
                        side,
                        tight_rhs,
                        a,
                        if dir > 0.0 { up } else { lo }
                    )),
                    rank: ratio,
                },
            );
        }
    }
}

/// min and max of the finite nonzero magnitudes of `vals`, with the argmin/argmax tags.
struct Range {
    min: f64,
    max: f64,
    min_at: usize,
    max_at: usize,
    count: usize,
}

fn range_of(vals: impl Iterator<Item = (usize, f64)>) -> Option<Range> {
    let mut r = Range {
        min: f64::INFINITY,
        max: 0.0,
        min_at: 0,
        max_at: 0,
        count: 0,
    };
    for (tag, x) in vals {
        let a = x.abs();
        if a == 0.0 || !a.is_finite() {
            continue;
        }
        r.count += 1;
        if a < r.min {
            r.min = a;
            r.min_at = tag;
        }
        if a > r.max {
            r.max = a;
            r.max_at = tag;
        }
    }
    (r.count > 0).then_some(r)
}

fn ratio_severity(ratio: f64) -> Option<Severity> {
    if ratio >= RANGE_RATIO_WARN {
        Some(Severity::Warning)
    } else if ratio >= RANGE_RATIO_INFO {
        Some(Severity::Info)
    } else {
        None
    }
}

/// `coef-range`: `Matrix [min, max]` of `HighsSolve.cpp:518-521`, with the thresholds of
/// `assessLp` (`large_matrix_value`/`small_matrix_value`, `HighsOptions.h:818-829`).
fn coef_range(v: &View, out: &mut Findings) {
    let it = v.rows.iter().enumerate().flat_map(|(i, r)| {
        r.iter()
            .map(move |&(j, a)| (i * (v.lp.num_col_ as usize).max(1) + j, a))
    });
    let nc = (v.lp.num_col_ as usize).max(1);
    let Some(r) = range_of(it) else { return };
    let ratio = r.max / r.min;
    let loc = |tag: usize| (v.rname(tag / nc), v.cname(tag % nc));
    let (min_row, min_col) = loc(r.min_at);
    let (max_row, max_col) = loc(r.max_at);
    let sev = if r.max >= LARGE_MATRIX_VALUE {
        Some(Severity::Error)
    } else if r.min < SMALL_MATRIX_VALUE || ratio >= RANGE_RATIO_WARN {
        Some(Severity::Warning)
    } else {
        ratio_severity(ratio)
    };
    let Some(sev) = sev else { return };
    out.push(
        "coef-range",
        sev,
        None,
        Item {
            rows: vec![min_row.clone(), max_row.clone()],
            cols: vec![min_col.clone(), max_col.clone()],
            evidence: ev(&[
                ("min_abs", num(r.min)),
                ("max_abs", num(r.max)),
                ("ratio", num(ratio)),
                ("min_at", json!(format!("{min_row}/{min_col}"))),
                ("max_at", json!(format!("{max_row}/{max_col}"))),
                ("nonzeros", json!(r.count)),
            ]),
            fix: None,
            rank: ratio,
        },
    );
    // the rows with the widest internal ratio carry the problem
    let mut worst: Vec<(f64, usize, f64, f64)> = Vec::new();
    for (i, row) in v.rows.iter().enumerate() {
        if let Some(rr) = range_of(row.iter().map(|&(j, a)| (j, a))) {
            let q = rr.max / rr.min;
            if q >= RANGE_RATIO_INFO {
                worst.push((q, i, rr.min, rr.max));
            }
        }
    }
    worst.sort_by(|a, b| b.0.total_cmp(&a.0));
    for (q, i, mn, mx) in worst.into_iter().take(5) {
        out.push(
            "coef-range",
            ratio_severity(q).unwrap_or(Severity::Info),
            Some("rows with the widest internal coefficient range".into()),
            Item {
                rows: vec![v.rname(i)],
                cols: vec![],
                evidence: ev(&[
                    ("min_abs", num(mn)),
                    ("max_abs", num(mx)),
                    ("ratio", num(q)),
                ]),
                fix: None,
                rank: q,
            },
        );
    }
}

/// Shared by `rhs-range`, `obj-range`, `bound-range`: severity from the HiGHS excessive-value
/// thresholds (`HighsSolve.cpp:565-607`) and the ratio levels.
#[allow(clippy::too_many_arguments)]
fn vec_range(
    out: &mut Findings,
    rule: &str,
    vals: impl Iterator<Item = (usize, f64)>,
    small: f64,
    large: f64,
    name: &dyn Fn(usize) -> String,
    is_row: bool,
) {
    let Some(r) = range_of(vals) else { return };
    let ratio = r.max / r.min;
    let mut sev = ratio_severity(ratio);
    if r.max > large {
        sev = Some(Severity::Warning);
    } else if r.min < small && sev.is_none() {
        sev = Some(Severity::Info);
    }
    let Some(sev) = sev else { return };
    let (a, b) = (name(r.min_at), name(r.max_at));
    let (rows, cols) = if is_row {
        (vec![a.clone(), b.clone()], vec![])
    } else {
        (vec![], vec![a.clone(), b.clone()])
    };
    out.push(
        rule,
        sev,
        None,
        Item {
            rows,
            cols,
            evidence: ev(&[
                ("min_abs", num(r.min)),
                ("max_abs", num(r.max)),
                ("ratio", num(ratio)),
                ("min_at", json!(a)),
                ("max_at", json!(b)),
                ("count", json!(r.count)),
                ("excessively_small_below", num(small)),
                ("excessively_large_above", num(large)),
            ]),
            fix: None,
            rank: ratio,
        },
    );
}

fn rhs_range(v: &View, out: &mut Findings) {
    let lp = v.lp;
    let vals =
        (0..lp.num_row_ as usize).flat_map(|i| [(i, lp.row_lower_[i]), (i, lp.row_upper_[i])]);
    vec_range(
        out,
        "rhs-range",
        vals,
        SMALL_BOUND,
        LARGE_BOUND,
        &|i| v.rname(i),
        true,
    );
}

fn obj_range(v: &View, out: &mut Findings) {
    let lp = v.lp;
    vec_range(
        out,
        "obj-range",
        (0..lp.num_col_ as usize).map(|j| (j, lp.col_cost_[j])),
        SMALL_OBJ,
        LARGE_OBJ,
        &|j| v.cname(j),
        false,
    );
}

fn bound_range(v: &View, out: &mut Findings) {
    let lp = v.lp;
    let nc = lp.num_col_ as usize;
    for (int, large) in [(false, LARGE_BOUND), (true, LARGE_INT_BOUND)] {
        // HiGHS reports continuous and integer column bounds separately (HighsSolve.cpp:588-608)
        let vals = (0..nc)
            .filter(move |&j| v.is_int(j) == int)
            .flat_map(|j| [(j, lp.col_lower_[j]), (j, lp.col_upper_[j])]);
        vec_range(
            out,
            "bound-range",
            vals,
            SMALL_BOUND,
            large,
            &|j| v.cname(j),
            false,
        );
    }
}

fn median(mut x: Vec<f64>) -> f64 {
    x.sort_by(|a, b| a.total_cmp(b));
    x[x.len() / 2]
}

/// `scaling`: each row's and column's geometric scale `sqrt(min|a| * max|a|)` (the inverse of the
/// equilibration factor `1/sqrt(min*max)` that HiGHS applies, `HighsLpUtils.cpp:983,997`; the
/// geometric-mean idea is Tomlin, "On scaling linear programming problems", Math. Programming
/// Study 4, 1975) against the median over all rows (columns). A row or column whose scale differs
/// from the median by orders of magnitude is the one a scaling pass would have to move most.
fn scaling(v: &View, out: &mut Findings) {
    let nr = v.rows.len();
    let nc = v.lp.num_col_ as usize;
    let mut row_scale: Vec<Option<f64>> = vec![None; nr];
    for (i, row) in v.rows.iter().enumerate() {
        if let Some(r) = range_of(row.iter().map(|&(j, a)| (j, a))) {
            row_scale[i] = Some((r.min * r.max).sqrt());
        }
    }
    let mut cmin = vec![f64::INFINITY; nc];
    let mut cmax = vec![0.0f64; nc];
    for row in &v.rows {
        for &(j, a) in row {
            cmin[j] = cmin[j].min(a.abs());
            cmax[j] = cmax[j].max(a.abs());
        }
    }
    let col_scale: Vec<Option<f64>> = (0..nc)
        .map(|j| (cmax[j] > 0.0).then(|| (cmin[j] * cmax[j]).sqrt()))
        .collect();
    for (kind, scales) in [("row", &row_scale), ("column", &col_scale)] {
        let present: Vec<f64> = scales.iter().flatten().copied().collect();
        if present.len() < 2 {
            continue;
        }
        let med = median(present);
        let mut hits: Vec<(f64, usize, f64)> = Vec::new();
        for (k, s) in scales.iter().enumerate() {
            if let Some(s) = *s {
                let f = (s / med).max(med / s);
                if f >= SCALE_RATIO_INFO {
                    hits.push((f, k, s));
                }
            }
        }
        hits.sort_by(|a, b| b.0.total_cmp(&a.0));
        for (f, k, s) in hits {
            let (rows, cols, nm) = if kind == "row" {
                (vec![v.rname(k)], vec![], v.rname(k))
            } else {
                (vec![], vec![v.cname(k)], v.cname(k))
            };
            out.push(
                "scaling",
                ratio_sev_scale(f),
                None,
                Item {
                    rows,
                    cols,
                    evidence: ev(&[
                        ("kind", json!(kind)),
                        ("scale", num(s)),
                        ("median_scale", num(med)),
                        ("factor_from_median", num(f)),
                        ("suggested_multiplier", num(med / s)),
                    ]),
                    fix: Some(format!(
                        "multiply {kind} {nm} by {:e} (its scale {:e} against the median {:e})",
                        med / s,
                        s,
                        med
                    )),
                    rank: f,
                },
            );
        }
    }
}

fn ratio_sev_scale(f: f64) -> Severity {
    if f >= SCALE_RATIO_WARN {
        Severity::Warning
    } else {
        Severity::Info
    }
}

/// `near-parallel-rows`: HiGHS's parallel-row detection (`HPresolve::detectParallelRowsAndCols`,
/// `HPresolveParallelRowsAndCols.cpp:15`) scales each row by its signed largest entry (the
/// first one, by column index, among those within `small_matrix_value` of the maximum, lines
/// 35-50) and hashes the scaled rows. Here the rows are bucketed by sparsity pattern, sorted
/// lexicographically on their scaled values inside a bucket, and adjacent rows compared with
/// the same `small_matrix_value` tolerance; rows of one entry are bounds, not constraints, and
/// are skipped.
fn near_parallel(v: &View, out: &mut Findings) {
    let lp = v.lp;
    let mut norm: Vec<Option<(f64, Vec<f64>)>> = vec![None; v.rows.len()];
    let mut buckets: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, row) in v.rows.iter().enumerate() {
        if row.len() < 2 {
            continue;
        }
        let mx = row.iter().fold(0.0f64, |m, &(_, a)| m.max(a.abs()));
        let s = row
            .iter()
            .find(|&&(_, a)| a.abs() >= mx - SMALL_MATRIX_VALUE)
            .map(|&(_, a)| a)
            .unwrap();
        // FNV-1a over the column pattern
        let mut h: u64 = 0xcbf29ce484222325;
        for &(j, _) in row {
            h = (h ^ (j as u64 + 1)).wrapping_mul(0x100000001b3);
        }
        norm[i] = Some((s, row.iter().map(|&(_, a)| a / s).collect()));
        buckets.entry(h).or_default().push(i);
    }
    for (_, mut ids) in buckets {
        if ids.len() < 2 {
            continue;
        }
        ids.sort_by(|&a, &b| {
            let (va, vb) = (&norm[a].as_ref().unwrap().1, &norm[b].as_ref().unwrap().1);
            va.iter()
                .zip(vb)
                .map(|(x, y)| x.total_cmp(y))
                .find(|o| o.is_ne())
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.cmp(&b))
        });
        for w in ids.windows(2) {
            let (a, b) = (w[0], w[1]);
            let (sa, va) = norm[a].as_ref().unwrap();
            let (sb, vb) = norm[b].as_ref().unwrap();
            // same pattern is guaranteed by the hash only with high probability: check it
            if v.rows[a]
                .iter()
                .map(|e| e.0)
                .ne(v.rows[b].iter().map(|e| e.0))
            {
                continue;
            }
            let diff = va
                .iter()
                .zip(vb)
                .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
            if diff > PARALLEL_TOL {
                continue;
            }
            // each row's bounds on the normalised expression `sum v_k x_k`
            let interval = |i: usize, s: f64| {
                let (lo, up) = (lp.row_lower_[i] / s, lp.row_upper_[i] / s);
                if s > 0.0 { (lo, up) } else { (up, lo) }
            };
            let (la, ua) = interval(a, *sa);
            let (lb, ub) = interval(b, *sb);
            let lo = la.max(lb);
            let up = ua.min(ub);
            let tol = PARALLEL_FEASTOL * (1.0 + lo.abs().min(1e30).max(up.abs().min(1e30)));
            let disjoint = lo > up + tol;
            // a's interval inside b's (or the reverse): one row implies the other
            let within = |l1: f64, u1: f64, l2: f64, u2: f64| {
                l1 >= l2 - PARALLEL_FEASTOL * (1.0 + l2.abs().min(1e30))
                    && u1 <= u2 + PARALLEL_FEASTOL * (1.0 + u2.abs().min(1e30))
            };
            let dominated = within(la, ua, lb, ub) || within(lb, ub, la, ua);
            let (sev, msg) = if disjoint {
                (
                    Severity::Error,
                    Some(
                        "parallel rows whose bounds do not intersect: the model is infeasible"
                            .to_string(),
                    ),
                )
            } else if dominated {
                (
                    Severity::Warning,
                    Some(
                        "parallel rows, one implies the other: the dominated row is redundant"
                            .to_string(),
                    ),
                )
            } else {
                (
                    Severity::Info,
                    Some(
                        "parallel rows with complementary bounds: together one ranged row"
                            .to_string(),
                    ),
                )
            };
            out.push(
                "near-parallel-rows",
                sev,
                msg,
                Item {
                    rows: vec![v.rname(a), v.rname(b)],
                    cols: v.rows[a].iter().take(6).map(|e| v.cname(e.0)).collect(),
                    evidence: ev(&[
                        ("max_normalised_difference", num(diff)),
                        ("scale_ratio", num(sa / sb)),
                        ("normalised_lower_a", num(la)),
                        ("normalised_upper_a", num(ua)),
                        ("normalised_lower_b", num(lb)),
                        ("normalised_upper_b", num(ub)),
                        ("merged_lower", num(lo)),
                        ("merged_upper", num(up)),
                        ("row_entries", json!(v.rows[a].len())),
                    ]),
                    fix: None,
                    rank: if disjoint {
                        1e300
                    } else {
                        v.rows[a].len() as f64
                    },
                },
            );
        }
    }
}

/// `empty-row`, `empty-col`, `free-row`, `free-col-no-cost`, `fixed-col`: the structural cases
/// HPresolve removes (empty rows and columns, free rows, free column singletons, fixed
/// columns), here only reported.
fn structure(v: &View, out: &mut Findings) {
    let lp = v.lp;
    let sense = if lp.sense_ == ObjSense::kMaximize {
        -1.0
    } else {
        1.0
    };
    for (i, row) in v.rows.iter().enumerate() {
        let (lo, up) = (lp.row_lower_[i], lp.row_upper_[i]);
        let free = lo == f64::NEG_INFINITY && up == f64::INFINITY;
        if row.is_empty() {
            let infeasible = lo > PRIMAL_FEASTOL || up < -PRIMAL_FEASTOL;
            out.push(
                "empty-row",
                if infeasible {
                    Severity::Error
                } else {
                    Severity::Warning
                },
                infeasible.then(|| {
                    "empty row whose bounds exclude zero: the model is infeasible".to_string()
                }),
                Item {
                    rows: vec![v.rname(i)],
                    cols: vec![],
                    evidence: ev(&[("lower", num(lo)), ("upper", num(up))]),
                    fix: None,
                    rank: infeasible as u8 as f64,
                },
            );
        } else if free {
            out.push(
                "free-row",
                Severity::Warning,
                None,
                Item {
                    rows: vec![v.rname(i)],
                    cols: vec![],
                    evidence: ev(&[("entries", json!(row.len()))]),
                    fix: None,
                    rank: row.len() as f64,
                },
            );
        }
    }
    for j in 0..lp.num_col_ as usize {
        let (l, u, cost) = (lp.col_lower_[j], lp.col_upper_[j], lp.col_cost_[j]);
        if v.col_nnz[j] == 0 {
            let e = cost * sense;
            let unbounded = (e > 0.0 && l == f64::NEG_INFINITY) || (e < 0.0 && u == f64::INFINITY);
            let (sev, msg) = if unbounded {
                (
                    Severity::Error,
                    Some("column in no row whose cost drives it to an infinite bound: the model is unbounded".to_string()),
                )
            } else if cost == 0.0 {
                (Severity::Warning, None)
            } else {
                (
                    Severity::Info,
                    Some("column in no row: its cost and bounds alone fix its value".to_string()),
                )
            };
            out.push(
                "empty-col",
                sev,
                msg,
                Item {
                    rows: vec![],
                    cols: vec![v.cname(j)],
                    evidence: ev(&[("cost", num(cost)), ("lower", num(l)), ("upper", num(u))]),
                    fix: None,
                    rank: unbounded as u8 as f64,
                },
            );
        } else if v.col_nnz[j] == 1
            && cost == 0.0
            && l == f64::NEG_INFINITY
            && u == f64::INFINITY
            && !v.is_int(j)
        {
            let row = v
                .rows
                .iter()
                .position(|r| r.iter().any(|e| e.0 == j))
                .unwrap();
            out.push(
                "free-col-no-cost",
                Severity::Warning,
                None,
                Item {
                    rows: vec![v.rname(row)],
                    cols: vec![v.cname(j)],
                    evidence: ev(&[
                        ("row_lower", num(lp.row_lower_[row])),
                        ("row_upper", num(lp.row_upper_[row])),
                    ]),
                    fix: None,
                    rank: 0.0,
                },
            );
        }
        if l == u && l.is_finite() {
            out.push(
                "fixed-col",
                Severity::Info,
                None,
                Item {
                    rows: vec![],
                    cols: vec![v.cname(j)],
                    evidence: ev(&[("value", num(l)), ("entries", json!(v.col_nnz[j]))]),
                    fix: None,
                    rank: v.col_nnz[j] as f64,
                },
            );
        }
    }
}

/// `fractional-int-bound`: HiGHS rounds integer bounds inward (`ceil` lower, `floor` upper) when
/// it assesses integrality (`HighsLpUtils.cpp` `assessBounds`/`HighsLp::applyScale` callers, and
/// MIP presolve `HPresolve::changeColLower/Upper` on integer columns).
fn fractional_int_bound(v: &View, out: &mut Findings) {
    let lp = v.lp;
    for j in 0..lp.num_col_ as usize {
        if !v.is_int(j) {
            continue;
        }
        let (l, u) = (lp.col_lower_[j], lp.col_upper_[j]);
        let frac = |x: f64| x.is_finite() && (x - x.round()).abs() > PRIMAL_FEASTOL;
        if !(frac(l) || frac(u)) {
            continue;
        }
        let (nl, nu) = (
            if frac(l) { l.ceil() } else { l },
            if frac(u) { u.floor() } else { u },
        );
        let crossed = nl > nu;
        out.push(
            "fractional-int-bound",
            if crossed {
                Severity::Error
            } else {
                Severity::Warning
            },
            crossed.then(|| {
                "integer bounds cross after rounding: the model is infeasible".to_string()
            }),
            Item {
                rows: vec![],
                cols: vec![v.cname(j)],
                evidence: ev(&[
                    ("lower", num(l)),
                    ("upper", num(u)),
                    ("rounded_lower", num(nl)),
                    ("rounded_upper", num(nu)),
                ]),
                fix: Some(format!("set the bounds of {} to [{nl}, {nu}]", v.cname(j))),
                rank: crossed as u8 as f64,
            },
        );
    }
}

/// `condition-estimate` (opt-in): solve the LP relaxation on a private `Highs` with presolve
/// off and read `HEkk::computeBasisCondition` (`HEkk.cpp:3649`, port:
/// `simplex/ekk/solve_return.rs`), the approximation `Highs::getKappa` returns
/// (`Highs.cpp:2699-2706`, `exact = false`).
fn condition_estimate(lp: &HighsLp, v: &View, out: &mut Findings) -> Result<(), String> {
    let mut relax = lp.clone();
    relax.integrality_.clear();
    let mut highs = Highs::default();
    for (k, val) in [
        ("presolve", "off"),
        ("output_flag", "false"),
        ("threads", "1"),
    ] {
        if highs.set_option_value_str(k, val) == HighsStatus::kError {
            return Err(format!("cannot set option {k}"));
        }
    }
    if highs.pass_model(relax) == HighsStatus::kError {
        return Err("pass_model failed".into());
    }
    if highs.run() == HighsStatus::kError {
        return Err("LP solve failed".into());
    }
    let status = highs.get_model_status();
    if status != HighsModelStatus::kOptimal || !highs.ekk_instance_.status_.has_invert {
        out.push(
            "condition-estimate",
            Severity::Info,
            Some(
                "no condition estimate: the LP relaxation did not end with an optimal basis"
                    .to_string(),
            ),
            Item {
                rows: vec![],
                cols: vec![],
                evidence: ev(&[("model_status", json!(format!("{status:?}")))]),
                fix: None,
                rank: 0.0,
            },
        );
        return Ok(());
    }
    let kappa = highs
        .ekk_instance_
        .compute_basis_condition_lp(&highs.model_lp_, false, false);
    let sev = if kappa >= KAPPA_ERROR {
        Severity::Error
    } else if kappa >= KAPPA_WARN {
        Severity::Warning
    } else if kappa >= KAPPA_INFO {
        Severity::Info
    } else {
        return Ok(());
    };
    out.push(
        "condition-estimate",
        sev,
        None,
        Item {
            rows: vec![],
            cols: vec![],
            evidence: ev(&[
                ("kappa_estimate", num(kappa)),
                ("rows", json!(v.rows.len())),
                (
                    "lp_iterations",
                    json!(highs.get_info().simplex_iteration_count),
                ),
            ]),
            fix: None,
            rank: kappa,
        },
    );
    Ok(())
}

// ---------------------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------------------

pub(crate) struct Args {
    path: String,
    json: bool,
    rules: Vec<usize>,
    min_severity: Severity,
    max_items: usize,
    with_solve: bool,
}

fn usage() -> &'static str {
    "usage: rusolve lint <model.(mps|lp|nl|mof.json)[.gz|.bz2|.xz|.zst]> [--format text|json] [--rule ID]... \
     [--min-severity info|warning|error] [--max-items N] [--with-solve]"
}

pub(crate) fn parse_args(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        path: String::new(),
        json: false,
        rules: Vec::new(),
        min_severity: Severity::Info,
        max_items: DEFAULT_MAX_ITEMS,
        with_solve: false,
    };
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let mut val = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value"))
        };
        match arg.as_str() {
            "--format" => match val("--format")?.as_str() {
                "json" => a.json = true,
                "text" => a.json = false,
                o => return Err(format!("--format expects text or json, got {o}")),
            },
            "--rule" => {
                let id = val("--rule")?;
                let k = rule_index(&id).ok_or_else(|| {
                    format!(
                        "unknown rule {id}; rules: {}",
                        RULES.iter().map(|r| r.0).collect::<Vec<_>>().join(", ")
                    )
                })?;
                a.rules.push(k);
            }
            "--min-severity" => {
                let s = val("--min-severity")?;
                a.min_severity = Severity::parse(&s).ok_or_else(|| {
                    format!("--min-severity expects info, warning or error, got {s}")
                })?;
            }
            "--max-items" => {
                let s = val("--max-items")?;
                a.max_items = s
                    .parse()
                    .map_err(|_| format!("--max-items expects a count, got {s}"))?;
            }
            "--with-solve" => a.with_solve = true,
            s if s.starts_with("--") => return Err(format!("unknown flag {s}")),
            s => {
                if !a.path.is_empty() {
                    return Err(format!("unexpected argument {s}"));
                }
                a.path = s.to_string();
            }
        }
    }
    if a.path.is_empty() {
        return Err("missing model path".into());
    }
    Ok(a)
}

/// HiGHS treats any bound at or beyond `infinite_bound` (default 1e20, `HighsOptions.h`) as infinite.
const INFINITE_BOUND: f64 = 1e20;

fn clamp_infinite(lp: &mut HighsLp) {
    for x in lp
        .col_lower_
        .iter_mut()
        .chain(lp.col_upper_.iter_mut())
        .chain(lp.row_lower_.iter_mut())
        .chain(lp.row_upper_.iter_mut())
    {
        if *x >= INFINITE_BOUND {
            *x = f64::INFINITY;
        } else if *x <= -INFINITE_BOUND {
            *x = f64::NEG_INFINITY;
        }
    }
}

pub fn run(args: &[String]) -> ExitCode {
    let a = match parse_args(args) {
        Ok(a) => a,
        Err(e) => {
            return crate::diag::fail(
                crate::diag::wants_json(args),
                rusolve_core::diagnostic::Diagnostic::error("RS0006", format!("rusolve lint: {e}"))
                    .with_help(usage()),
                2,
            );
        }
    };
    let lp = match crate::mip_highs::read_mip(&a.path) {
        Ok(lp) => lp,
        Err(e) => {
            return crate::diag::fail_read(
                a.json,
                &a.path,
                e.strip_prefix(&format!("read {}: ", a.path)).unwrap_or(&e),
            );
        }
    };
    match lint_core(lp, &a) {
        Ok((out, text)) => {
            if a.json {
                println!("{out}");
            } else {
                print!("{text}");
            }
            ExitCode::SUCCESS
        }
        Err(e) => {
            let code = if e.code == 2 { "RS0003" } else { "RS0090" };
            crate::diag::fail_msg(a.json, code, e.code, format!("rusolve lint: {}", e.msg))
        }
    }
}

/// The body of `rusolve lint` on a loaded model: the `rusolve.lint/1` object
/// and the text report. Shared by the CLI and the MCP `lint` tool.
pub(crate) fn lint_core(
    mut lp: HighsLp,
    a: &Args,
) -> Result<(Value, String), crate::cli_core::CoreError> {
    lp.a_matrix_.ensure_colwise();
    clamp_infinite(&mut lp);
    let v = View::new(&lp);
    let want = |id: &str| a.rules.is_empty() || a.rules.contains(&rule_index(id).unwrap());
    let mut f = Findings::default();
    if want("big-m") {
        big_m(&v, &mut f);
    }
    if want("coef-range") {
        coef_range(&v, &mut f);
    }
    if want("rhs-range") {
        rhs_range(&v, &mut f);
    }
    if want("obj-range") {
        obj_range(&v, &mut f);
    }
    if want("bound-range") {
        bound_range(&v, &mut f);
    }
    if want("scaling") {
        scaling(&v, &mut f);
    }
    if want("near-parallel-rows") {
        near_parallel(&v, &mut f);
    }
    if [
        "empty-row",
        "empty-col",
        "free-row",
        "free-col-no-cost",
        "fixed-col",
    ]
    .iter()
    .any(|r| want(r))
    {
        let mut s = Findings::default();
        structure(&v, &mut s);
        f.groups
            .extend(s.groups.into_iter().filter(|g| want(RULES[g.rule].0)));
    }
    if want("fractional-int-bound") {
        fractional_int_bound(&v, &mut f);
    }
    if a.with_solve
        && want("condition-estimate")
        && let Err(e) = condition_estimate(&lp, &v, &mut f)
    {
        return Err(crate::cli_core::CoreError::internal(format!(
            "condition-estimate: {e}"
        )));
    }
    // order: severity descending, then rule order
    f.groups.retain(|g| g.severity >= a.min_severity);
    f.groups
        .sort_by(|x, y| y.severity.cmp(&x.severity).then(x.rule.cmp(&y.rule)));
    for g in &mut f.groups {
        g.items.sort_by(|x, y| y.rank.total_cmp(&x.rank));
    }
    let count = |s: Severity| {
        f.groups
            .iter()
            .filter(|g| g.severity == s)
            .map(|g| g.items.len())
            .sum::<usize>()
    };
    let (ne, nw, ni) = (
        count(Severity::Error),
        count(Severity::Warning),
        count(Severity::Info),
    );
    let json_out = {
        let findings: Vec<Value> = f
            .groups
            .iter()
            .map(|g| {
                json!({
                    "rule": RULES[g.rule].0,
                    "severity": g.severity.name(),
                    "message": g.message,
                    "fix": RULES[g.rule].2,
                    "count": g.items.len(),
                    "items": g.items.iter().take(a.max_items).map(|it| {
                        let mut o = json!({
                            "rows": it.rows,
                            "columns": it.cols,
                            "evidence": it.evidence,
                        });
                        if let Some(fix) = &it.fix {
                            o["fix"] = json!(fix);
                        }
                        o
                    }).collect::<Vec<_>>(),
                })
            })
            .collect();
        let rules_run: Vec<&str> = RULES
            .iter()
            .enumerate()
            .filter(|(k, r)| {
                (a.rules.is_empty() || a.rules.contains(k))
                    && (r.0 != "condition-estimate" || a.with_solve)
            })
            .map(|(_, r)| r.0)
            .collect();
        let out = json!({
            "schema_version": "rusolve.lint/1",
            "input": {
                "path": a.path,
                "rows": lp.num_row_,
                "columns": lp.num_col_,
                "nonzeros": v.nnz,
                "integer_columns": v.num_int(),
            },
            "rules_run": rules_run,
            "min_severity": a.min_severity.name(),
            "summary": {"errors": ne, "warnings": nw, "infos": ni},
            "findings": findings,
        });
        out
    };
    let mut text = String::new();
    macro_rules! tprintln {
        ($($arg:tt)*) => {{
            text.push_str(&format!($($arg)*));
            text.push('\n');
        }};
    }
    tprintln!(
        "rusolve lint {}: {} rows, {} columns, {} nonzeros, {} integer",
        a.path,
        lp.num_row_,
        lp.num_col_,
        v.nnz,
        v.num_int()
    );
    for g in &f.groups {
        tprintln!(
            "[{}] {} ({}): {}",
            g.severity.name(),
            RULES[g.rule].0,
            g.items.len(),
            g.message
        );
        for it in g.items.iter().take(a.max_items) {
            let mut s = String::from("    ");
            if !it.rows.is_empty() {
                s += &format!("rows {} ", it.rows.join(","));
            }
            if !it.cols.is_empty() {
                s += &format!("cols {} ", it.cols.join(","));
            }
            let e: Vec<String> = it
                .evidence
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect();
            tprintln!("{s}| {}", e.join(" "));
            if let Some(fix) = &it.fix {
                tprintln!("      fix: {fix}");
            }
        }
        if g.items.len() > a.max_items {
            tprintln!("    ... {} more", g.items.len() - a.max_items);
        }
        tprintln!("    fix: {}", RULES[g.rule].2);
    }
    tprintln!("summary: {ne} error(s), {nw} warning(s), {ni} info");
    Ok((json_out, text))
}

/// The cheap lint rules as diagnostics for `solve` and `mip` (warnings on stderr and the
/// envelope's `diagnostics` array). Runs `big-m`, `coef-range`, `bound-range` and
/// `free-col-no-cost` on a copy of the model, so nothing here can touch the search; one
/// diagnostic per finding group, at most [`WARN_ITEMS`] items named. Skipped (empty) when the
/// matrix is above [`WARN_MAX_NNZ`] nonzeros, where even a linear pass is not free.
/// Source of the rules: the `big_m` and `coef_range` docs above; the diagnostic layout is
/// rustc's (<https://rustc-dev-guide.rust-lang.org/diagnostics.html>).
pub(crate) fn warnings(lp: &HighsLp) -> Vec<rusolve_core::diagnostic::Diagnostic> {
    use rusolve_core::diagnostic::{Diagnostic, Location};
    if lp.a_matrix_.value_.len() > WARN_MAX_NNZ {
        return Vec::new();
    }
    let mut lp = lp.clone();
    lp.a_matrix_.ensure_colwise();
    clamp_infinite(&mut lp);
    let v = View::new(&lp);
    let mut f = Findings::default();
    big_m(&v, &mut f);
    coef_range(&v, &mut f);
    bound_range(&v, &mut f);
    let mut s = Findings::default();
    structure(&v, &mut s);
    f.groups.extend(
        s.groups
            .into_iter()
            .filter(|g| RULES[g.rule].0 == "free-col-no-cost"),
    );
    f.groups.retain(|g| g.severity >= Severity::Warning);
    f.groups
        .sort_by(|x, y| y.severity.cmp(&x.severity).then(x.rule.cmp(&y.rule)));
    let mut out = Vec::new();
    for g in &mut f.groups {
        g.items.sort_by(|x, y| y.rank.total_cmp(&x.rank));
        let rule = RULES[g.rule].0;
        let code = match rule {
            "big-m" => "RS2001",
            "coef-range" => "RS2002",
            _ => "RS2003",
        };
        let mut msg = format!("{rule}: {}", g.message);
        if rule == "big-m"
            && let Some(it) = g.items.first()
        {
            let num_of = |k: &str| it.evidence.get(k).and_then(Value::as_f64);
            if let (Some(r), Some(col), Some(a), Some(t)) = (
                it.rows.first(),
                it.cols.first(),
                num_of("coefficient"),
                num_of("tightest_coefficient"),
            ) {
                let binary =
                    num_of("column_lower") == Some(0.0) && num_of("column_upper") == Some(1.0);
                let kind = if binary { "binary" } else { "integer" };
                msg = format!(
                    "{rule}: row '{r}': coefficient {} on {kind} '{col}' is far larger in magnitude than the tightest valid value {}",
                    short(a),
                    short(t)
                );
                if g.items.len() > 1 {
                    msg.push_str(&format!(" (and {} more)", g.items.len() - 1));
                }
            }
        }
        let mut d = Diagnostic::warning(code, msg);
        if let Some(it) = g.items.first() {
            // the location line shows the row or column as written; the numbers behind the
            // finding stay in `rusolve lint --format json`
            if let Some(r) = it.rows.first() {
                let rendered = (0..v.rows.len())
                    .find(|&i| v.rname(i) == *r)
                    .map_or_else(|| r.clone(), |i| v.render_row(i));
                d = d.with_location(Location::Row {
                    name: r.clone(),
                    rendered,
                });
            } else if let Some(cn) = it.cols.first() {
                let rendered = (0..v.lp.num_col_ as usize)
                    .find(|&j| v.cname(j) == *cn)
                    .map_or_else(|| cn.clone(), |j| v.render_col(j));
                d = d.with_location(Location::Column {
                    name: cn.clone(),
                    rendered,
                });
            }
        }
        d = d.with_note(format!("{} finding(s) of rule {rule}", g.items.len()));
        for it in g.items.iter().take(WARN_ITEMS) {
            if let Some(fix) = &it.fix {
                d = d.with_help(fix.clone());
            }
        }
        d = d.with_help(format!(
            "`rusolve lint <model>` lists every finding with {rule}; --no-warnings silences these"
        ));
        out.push(d);
    }
    out
}

/// Items named per warning group.
const WARN_ITEMS: usize = 3;
/// Above this many nonzeros the solve-time warnings are skipped.
const WARN_MAX_NNZ: usize = 2_000_000;
