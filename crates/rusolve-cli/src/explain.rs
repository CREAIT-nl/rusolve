//! `rusolve explain <model.mps> [--format json] [--no-verify] [--max-members N]`
//! (Farkas certificates, unbounded rays and the LP infeasible set).
//!
//! Solves the model as an LP (a MIP's integrality is dropped and the result is
//! labelled `"kind": "relaxation"`; a true MIP IIS is later work) and, when
//! it is infeasible or unbounded, returns an explanation that is checked before
//! it is printed.
//!
//! ## Sources (HiGHS 1.16.0, `highs/`)
//!
//! * IIS: `Highs::getIis` (`lp_data/Highs.cpp:2238-2242`) ->
//!   `Highs::getIisInterface` (`lp_data/HighsInterface.cpp:1929-2153`, port:
//!   `rusolve-ekk/src/lp_data/highs_iis_api.rs`). With the default
//!   `iis_strategy = kIisStrategyLight` (0) it only runs the trivial and
//!   row-value-bound tests (`HighsInterface.cpp:2000`), so a general infeasible
//!   LP gets no set. The elastic-LP route is `kIisStrategyFromLp` (2,
//!   `HighsInterface.cpp:2077`) and `kIisStrategyIrreducible` (4,
//!   `HighsInterface.cpp:2139`) shrinks the set to a true IIS. We set 2+4=6.
//!   (`kIisStrategyFromRay` is commented out in 1.16.0, `:2075`.) The set is
//!   `HighsIis::{col_index_, row_index_, col_bound_, row_bound_}`
//!   (`HighsIis.h`); a bound status says which of the bounds is in the set.
//! * Farkas: `Highs::getDualRay` (`Highs.cpp:2141-2144`) ->
//!   `getDualRayInterface` (`HighsInterface.cpp:1580-1726`): one multiplier per
//!   row, re-solving a zero-cost feasibility problem without presolve when the
//!   ray is not at hand.
//! * Ray: `Highs::getPrimalRay` (`Highs.cpp:2185-2189`) ->
//!   `getPrimalRayInterface` (`HighsInterface.cpp:1735-`): one entry per column.
//!
//! ## Verification, independent of HiGHS's own `lpOk`
//!
//! * the set is infeasible: the model restricted to it (rows outside the set
//!   free, columns outside the set free, each bound outside the set dropped, zero
//!   cost) is re-solved;
//! * the set is irreducible: for each bound in the set, drop just that bound and
//!   re-solve; the rest must be feasible;
//! * Farkas: for multipliers `y` the combination `y'Ax` is confined to
//!   `[ylo, yhi]` by the row bounds and to `[clo, chi]` by the column bounds
//!   (`c = A'y`); `0 >= gap > 0` follows when the intervals are disjoint
//!   (Farkas' lemma for `{l <= Ax <= u, lb <= x <= ub}`; Chvatal, *Linear
//!   Programming*, ch. 6). The sign of HiGHS's ray is not assumed: both `y` and
//!   `-y` are tried;
//! * ray `d`: `A d` must not leave the row-bound directions, `d` must not leave
//!   the column-bound directions, and the objective must strictly improve along
//!   it (Rockafellar, *Convex Analysis*, Thm 8.7: a recession direction).
//!
//! ## Exit codes
//!
//! 0 ran (whatever the status; `verified:false` is reported, not signalled),
//! 2 usage error, 3 input error, 4 internal error.

use crate::mip_highs::read_mip;
use rusolve_ekk::lp_data::highs_iis::{
    HighsIis, kIisBoundStatusBoxed, kIisBoundStatusFree, kIisBoundStatusLower,
    kIisBoundStatusUpper, kIisModelStatusIrreducible, kIisModelStatusReducible,
};
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::{Map, Value, json};
use std::process::ExitCode;

/// `schema_version` of the JSON output; the schema is
/// `docs/schema/rusolve-explain-1.schema.json`.
pub const SCHEMA_VERSION: &str = "rusolve.explain/1";

const USAGE: &str = "usage: rusolve explain <model.mps> [--format text|json] [--no-verify] [--max-members N] [--mip-iis [--per-test-node-limit N] [--time-limit SECONDS]]";

/// `iis_strategy` for `getIis`: kIisStrategyFromLp (2) + kIisStrategyIrreducible (4).
const IIS_STRATEGY: i32 = 6;

pub(crate) struct Args {
    pub(crate) path: String,
    pub(crate) json: bool,
    pub(crate) verify: bool,
    pub(crate) max_members: usize,
    pub(crate) mip_iis: bool,
    pub(crate) per_test_node_limit: i32,
    pub(crate) time_limit: Option<f64>,
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut path = None;
    let mut json = false;
    let mut verify = true;
    let mut max_members = 1000usize;
    let mut mip_iis = false;
    let mut per_test_node_limit = 10000i32;
    let mut time_limit = None;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--format" => match it.next().map(String::as_str) {
                Some("json") => json = true,
                Some("text") => json = false,
                other => return Err(format!("--format wants text or json, got {other:?}")),
            },
            "--no-verify" => verify = false,
            "--mip-iis" => mip_iis = true,
            "--per-test-node-limit" => {
                let v = it.next().ok_or("--per-test-node-limit needs a value")?;
                per_test_node_limit = v
                    .parse()
                    .map_err(|_| format!("--per-test-node-limit: not a number: {v}"))?;
            }
            "--time-limit" => {
                let v = it.next().ok_or("--time-limit needs a value")?;
                let t: f64 = v
                    .parse()
                    .map_err(|_| format!("--time-limit: not a number: {v}"))?;
                if t.is_nan() || t < 0.0 {
                    return Err(format!("--time-limit: must be >= 0: {v}"));
                }
                time_limit = Some(t);
            }
            "--max-members" => {
                let v = it.next().ok_or("--max-members needs a value")?;
                max_members = v
                    .parse()
                    .map_err(|_| format!("--max-members: not a number: {v}"))?;
            }
            s if s.starts_with("--") => return Err(format!("unknown option {s}")),
            s => {
                if path.replace(s.to_string()).is_some() {
                    return Err("more than one model path".into());
                }
            }
        }
    }
    Ok(Args {
        path: path.ok_or("missing model path")?,
        json,
        verify,
        max_members,
        mip_iis,
        per_test_node_limit,
        time_limit,
    })
}

// ---------------------------------------------------------------- model helpers

pub(crate) fn row_name(lp: &HighsLp, i: usize) -> String {
    lp.row_names_
        .get(i)
        .cloned()
        .unwrap_or_else(|| format!("row{i}"))
}

pub(crate) fn col_name(lp: &HighsLp, j: usize) -> String {
    lp.col_names_
        .get(j)
        .cloned()
        .unwrap_or_else(|| format!("col{j}"))
}

/// Nonzeros `(row, col, value)` of the constraint matrix, either format.
pub(crate) fn entries(lp: &HighsLp) -> Vec<(usize, usize, f64)> {
    let m = &lp.a_matrix_;
    let colwise = m.is_colwise();
    let n_major = if colwise { lp.num_col_ } else { lp.num_row_ } as usize;
    let mut out = Vec::with_capacity(m.value_.len());
    for major in 0..n_major {
        for k in m.start_[major] as usize..m.start_[major + 1] as usize {
            let minor = m.index_[k] as usize;
            let (i, j) = if colwise {
                (minor, major)
            } else {
                (major, minor)
            };
            out.push((i, j, m.value_[k]));
        }
    }
    out
}

pub(crate) fn fin(v: f64) -> Value {
    if v.is_finite() { json!(v) } else { Value::Null }
}

pub(crate) fn side_name(lo: bool, hi: bool) -> &'static str {
    match (lo, hi) {
        (true, true) => "both",
        (true, false) => "lower",
        (false, true) => "upper",
        (false, false) => "none",
    }
}

// ---------------------------------------------------------------- feasibility oracle

/// Is the zero-cost LP `lp` (bounds already edited) infeasible?
/// `Ok(true)` infeasible, `Ok(false)` feasible, `Err` the solver could not tell.
/// Inconsistent bounds are decided by inspection (HiGHS `assessBounds` does not
/// reject `lower > upper`, presolve detects it; checking first keeps the oracle
/// independent of that path).
pub(crate) fn is_infeasible(lp: &HighsLp) -> Result<bool, String> {
    for j in 0..lp.num_col_ as usize {
        if lp.col_lower_[j] > lp.col_upper_[j] {
            return Ok(true);
        }
    }
    for i in 0..lp.num_row_ as usize {
        if lp.row_lower_[i] > lp.row_upper_[i] {
            return Ok(true);
        }
    }
    let mut h = Highs::with_options(crate::base_options());
    h.set_option_value_bool("output_flag", false);
    if h.pass_model(lp.clone()) == HighsStatus::kError {
        return Err("pass_model failed".into());
    }
    if h.run() == HighsStatus::kError {
        return Err("run failed".into());
    }
    match h.get_model_status() {
        HighsModelStatus::kInfeasible => Ok(true),
        HighsModelStatus::kOptimal | HighsModelStatus::kModelEmpty => Ok(false),
        s => Err(format!("model status {}", util_model_status_to_string(s))),
    }
}

// ---------------------------------------------------------------- IIS

/// One IIS member: which of its bounds are in the set.
pub(crate) struct Member {
    pub(crate) index: usize,
    pub(crate) lower: bool,
    pub(crate) upper: bool,
    /// a column in the set that carries no bound of its own (kIisBoundStatusFree)
    pub(crate) free_col: bool,
}

pub(crate) fn decode(status: HighsInt) -> (bool, bool) {
    (
        status == kIisBoundStatusLower || status == kIisBoundStatusBoxed,
        status == kIisBoundStatusUpper || status == kIisBoundStatusBoxed,
    )
}

/// The restricted LP of a candidate set: zero cost, rows outside the set free,
/// bounds outside the set dropped, with the bound named by `skip` also dropped.
/// `skip` is `(is_row, index, lower_side)`.
fn restricted(
    lp: &HighsLp,
    rows: &[Member],
    cols: &[Member],
    skip: Option<(bool, usize, bool)>,
) -> HighsLp {
    let mut r = lp.clone();
    r.col_cost_.iter_mut().for_each(|x| *x = 0.0);
    r.offset_ = 0.0;
    r.integrality_.clear();
    r.row_lower_.iter_mut().for_each(|b| *b = -kHighsInf);
    r.row_upper_.iter_mut().for_each(|b| *b = kHighsInf);
    r.col_lower_.iter_mut().for_each(|b| *b = -kHighsInf);
    r.col_upper_.iter_mut().for_each(|b| *b = kHighsInf);
    for m in rows {
        if m.lower && skip != Some((true, m.index, true)) {
            r.row_lower_[m.index] = lp.row_lower_[m.index];
        }
        if m.upper && skip != Some((true, m.index, false)) {
            r.row_upper_[m.index] = lp.row_upper_[m.index];
        }
    }
    for m in cols {
        if m.lower && skip != Some((false, m.index, true)) {
            r.col_lower_[m.index] = lp.col_lower_[m.index];
        }
        if m.upper && skip != Some((false, m.index, false)) {
            r.col_upper_[m.index] = lp.col_upper_[m.index];
        }
    }
    r
}

struct Verification {
    performed: bool,
    set_infeasible: Option<bool>,
    irreducible: Option<bool>,
    checks: Vec<Value>,
    failures: Vec<String>,
}

impl Verification {
    fn verified(&self) -> bool {
        self.performed && self.failures.is_empty()
    }
}

fn verify_iis(lp: &HighsLp, rows: &[Member], cols: &[Member]) -> Verification {
    let mut v = Verification {
        performed: true,
        set_infeasible: None,
        irreducible: None,
        checks: Vec::new(),
        failures: Vec::new(),
    };
    match is_infeasible(&restricted(lp, rows, cols, None)) {
        Ok(true) => v.set_infeasible = Some(true),
        Ok(false) => {
            v.set_infeasible = Some(false);
            v.failures
                .push("the set is feasible when re-solved".to_string());
        }
        Err(e) => {
            v.failures
                .push(format!("could not decide feasibility of the set: {e}"));
        }
    }
    if v.set_infeasible != Some(true) {
        return v;
    }
    let mut irreducible = true;
    let mut check = |is_row: bool, m: &Member, lower: bool, v: &mut Verification| {
        let name = if is_row {
            row_name(lp, m.index)
        } else {
            col_name(lp, m.index)
        };
        let result = is_infeasible(&restricted(lp, rows, cols, Some((is_row, m.index, lower))));
        let (feasible, note) = match result {
            Ok(false) => (true, None),
            Ok(true) => (false, Some("still infeasible without it")),
            Err(_) => (false, Some("solver could not decide")),
        };
        if !feasible {
            irreducible = false;
            v.failures.push(format!(
                "dropping the {} bound of {} {name}: {}",
                if lower { "lower" } else { "upper" },
                if is_row { "row" } else { "column" },
                note.unwrap_or("")
            ));
        }
        v.checks.push(json!({
            "dropped": {"type": if is_row {"row"} else {"column"}, "name": name,
                        "bound": if lower {"lower"} else {"upper"}},
            "rest_feasible": feasible,
        }));
    };
    for m in rows {
        if m.lower {
            check(true, m, true, &mut v);
        }
        if m.upper {
            check(true, m, false, &mut v);
        }
    }
    for m in cols {
        if m.lower {
            check(false, m, true, &mut v);
        }
        if m.upper {
            check(false, m, false, &mut v);
        }
    }
    v.irreducible = Some(irreducible);
    v
}

pub(crate) fn members_json(
    lp: &HighsLp,
    ms: &[Member],
    is_row: bool,
    limit: usize,
    truncated: &mut bool,
) -> Vec<Value> {
    if ms.len() > limit {
        *truncated = true;
    }
    ms.iter()
        .take(limit)
        .map(|m| {
            let (name, lo, hi) = if is_row {
                (
                    row_name(lp, m.index),
                    lp.row_lower_[m.index],
                    lp.row_upper_[m.index],
                )
            } else {
                (
                    col_name(lp, m.index),
                    lp.col_lower_[m.index],
                    lp.col_upper_[m.index],
                )
            };
            let mut o = Map::new();
            o.insert("name".into(), json!(name));
            o.insert("index".into(), json!(m.index));
            o.insert("bound".into(), json!(side_name(m.lower, m.upper)));
            if m.free_col {
                o.insert("note".into(), json!("in the set with none of its bounds"));
            }
            o.insert("lower".into(), if m.lower { fin(lo) } else { Value::Null });
            o.insert("upper".into(), if m.upper { fin(hi) } else { Value::Null });
            Value::Object(o)
        })
        .collect()
}

// ---------------------------------------------------------------- Farkas

struct Farkas {
    multipliers: Vec<(usize, f64)>,
    row_lo: f64,
    row_hi: f64,
    col_lo: f64,
    col_hi: f64,
    /// `> 0` iff the combination is contradictory; `NaN` when an infinite bound is involved
    gap: f64,
    /// `"row_lower_above_col_upper"` or `"col_lower_above_row_upper"`
    direction: &'static str,
    combined: Vec<(usize, f64)>,
    tolerance: f64,
    proves: bool,
}

/// Interval of `sum_i y_i * r_i` over `r in [rl, ru]`.
fn combination_range(y: &[(usize, f64)], lo: &[f64], hi: &[f64]) -> (f64, f64) {
    let (mut a, mut b) = (0.0, 0.0);
    for &(i, yi) in y {
        let (l, u) = if yi > 0.0 {
            (lo[i], hi[i])
        } else {
            (hi[i], lo[i])
        };
        a += yi * l;
        b += yi * u;
    }
    (a, b)
}

/// Farkas check of one sign of the ray; see the module comment.
fn farkas_for(lp: &HighsLp, y: &[(usize, f64)]) -> Farkas {
    let ymax = y.iter().fold(0.0f64, |m, &(_, v)| m.max(v.abs()));
    let mut cvec = vec![0.0f64; lp.num_col_ as usize];
    let mut touched = vec![false; lp.num_col_ as usize];
    let mut yv = vec![0.0f64; lp.num_row_ as usize];
    for &(i, v) in y {
        yv[i] = v;
    }
    for (i, j, a) in entries(lp) {
        if yv[i] != 0.0 {
            cvec[j] += yv[i] * a;
            touched[j] = true;
        }
    }
    // A column whose combined coefficient is rounding noise relative to the ray
    // is zero (an infinite bound times 1e-17 must not void a proof).
    let noise = 1e-9 * ymax;
    let combined: Vec<(usize, f64)> = (0..cvec.len())
        .filter(|&j| touched[j] && cvec[j].abs() > noise)
        .map(|j| (j, cvec[j]))
        .collect();
    let (row_lo, row_hi) = combination_range(y, &lp.row_lower_, &lp.row_upper_);
    let (col_lo, col_hi) = combination_range(&combined, &lp.col_lower_, &lp.col_upper_);
    // Rows give  row_lo <= c'x <= row_hi ; columns give col_lo <= c'x <= col_hi.
    let g1 = row_lo - col_hi;
    let g2 = col_lo - row_hi;
    let (gap, direction) = if g1 >= g2 {
        (g1, "row_lower_above_col_upper")
    } else {
        (g2, "col_lower_above_row_upper")
    };
    let scale = [row_lo, row_hi, col_lo, col_hi]
        .iter()
        .filter(|v| v.is_finite())
        .fold(1.0f64, |m, v| m.max(v.abs()));
    let tolerance = 1e-9 * scale;
    Farkas {
        multipliers: y.to_vec(),
        row_lo,
        row_hi,
        col_lo,
        col_hi,
        gap,
        direction,
        combined,
        tolerance,
        proves: gap.is_finite() && gap > tolerance,
    }
}

fn farkas_json(lp: &HighsLp, f: &Farkas, limit: usize) -> Value {
    let mut ms = f.multipliers.clone();
    ms.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
    let mut comb = f.combined.clone();
    comb.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
    json!({
        "available": true,
        "num_nonzero": f.multipliers.len(),
        "multipliers": ms.iter().take(limit).map(|&(i, v)| json!({
            "row": row_name(lp, i), "index": i, "multiplier": v})).collect::<Vec<_>>(),
        "truncated": f.multipliers.len() > limit,
        "proof": {
            "statement": "the rows combined with these multipliers give a linear form c'x with c = A'y; the row bounds confine it to [row_lo, row_hi], the column bounds to [col_lo, col_hi]; the intervals are disjoint, so 0 >= gap > 0",
            "combined_form": comb.iter().take(limit).map(|&(j, v)| json!({
                "column": col_name(lp, j), "index": j, "coefficient": v})).collect::<Vec<_>>(),
            "combined_form_nonzero": f.combined.len(),
            "row_bounds_interval": [fin(f.row_lo), fin(f.row_hi)],
            "column_bounds_interval": [fin(f.col_lo), fin(f.col_hi)],
            "direction": f.direction,
            "gap": fin(f.gap),
            "tolerance": f.tolerance,
            "verified": f.proves,
        },
    })
}

// ---------------------------------------------------------------- primal ray

struct RayCheck {
    objective_rate: f64,
    max_row_violation: f64,
    max_col_violation: f64,
    tolerance: f64,
    failures: Vec<String>,
}

fn check_ray(lp: &HighsLp, d: &[f64]) -> RayCheck {
    let dmax = d.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let tol = 1e-7 * dmax.max(1.0);
    let mut ad = vec![0.0f64; lp.num_row_ as usize];
    for (i, j, a) in entries(lp) {
        ad[i] += a * d[j];
    }
    let mut failures = Vec::new();
    // Direction `v` is allowed by the bounds `[lo, hi]` iff it does not point out
    // of a finite side that is the only one, or is zero when both are finite.
    let viol = |v: f64, lo: f64, hi: f64| -> f64 {
        match (lo.is_finite(), hi.is_finite()) {
            (true, true) => v.abs(),
            (true, false) => (-v).max(0.0),
            (false, true) => v.max(0.0),
            (false, false) => 0.0,
        }
    };
    let (mut mr, mut mc) = (0.0f64, 0.0f64);
    for (i, &adi) in ad.iter().enumerate() {
        let e = viol(adi, lp.row_lower_[i], lp.row_upper_[i]);
        if e > tol {
            failures.push(format!(
                "A d leaves the bound direction of row {} by {e:.3e}",
                row_name(lp, i)
            ));
        }
        mr = mr.max(e);
    }
    for (j, &dj) in d.iter().enumerate() {
        let e = viol(dj, lp.col_lower_[j], lp.col_upper_[j]);
        if e > tol {
            failures.push(format!(
                "d leaves the bound direction of column {} by {e:.3e}",
                col_name(lp, j)
            ));
        }
        mc = mc.max(e);
    }
    let rate: f64 = lp.col_cost_.iter().zip(d).map(|(cc, x)| cc * x).sum();
    let improving = match lp.sense_ {
        ObjSense::kMinimize => rate < -tol * 1e-2,
        ObjSense::kMaximize => rate > tol * 1e-2,
    };
    if !improving {
        failures.push(format!(
            "the objective does not improve along the ray (rate {rate:.3e})"
        ));
    }
    RayCheck {
        objective_rate: rate,
        max_row_violation: mr,
        max_col_violation: mc,
        tolerance: tol,
        failures,
    }
}

// ---------------------------------------------------------------- driver

/// Solve the LP with the IIS strategy set, resolving "unbounded or infeasible"
/// without presolve. Shared with `diagnose`.
pub(crate) fn solve_for_explanation(lp: HighsLp) -> Result<(Highs, HighsModelStatus), ()> {
    let mut h = Highs::with_options(crate::base_options());
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_int("iis_strategy", IIS_STRATEGY);
    if h.pass_model(lp) == HighsStatus::kError || h.run() == HighsStatus::kError {
        return Err(());
    }
    let mut status = h.get_model_status();
    if status == HighsModelStatus::kUnboundedOrInfeasible {
        // HiGHS may stop at "unbounded or infeasible" when presolve cannot tell
        // (`allow_unbounded_or_infeasible`, Highs.cpp:4855-4872); resolve without presolve.
        h.set_option_value_str("presolve", "off");
        if h.run() == HighsStatus::kError {
            return Err(());
        }
        status = h.get_model_status();
    }
    Ok((h, status))
}

pub fn run(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let a = match parse(args) {
        Ok(a) => a,
        Err(e) => return crate::diag::fail_usage_text(wj, "explain", &format!("{e}\n{USAGE}")),
    };
    let lp = match read_mip(&a.path) {
        Ok(lp) => lp,
        Err(e) => return crate::diag::fail_read_str(a.json, &a.path, &e),
    };
    let num_integer = lp
        .integrality_
        .iter()
        .filter(|t| **t != HighsVarType::kContinuous)
        .count();
    if a.mip_iis && num_integer > 0 {
        return crate::explain_mip::run_mip(&a, lp, num_integer);
    }
    match explain_core(lp, &a) {
        Ok((out, text)) => {
            if a.json {
                println!("{out}");
            } else {
                print!("{text}");
            }
            ExitCode::SUCCESS
        }
        Err(e) => crate::diag::fail_core(a.json, "explain", &e),
    }
}

/// The `Args` of [`explain_core`] for a model that is already loaded (the MCP
/// server). `path` is only the name reported in `model.path`.
pub(crate) fn core_args(path: &str, verify: bool, max_members: usize) -> Args {
    Args {
        path: path.to_string(),
        json: true,
        verify,
        max_members,
        mip_iis: false,
        per_test_node_limit: 10000,
        time_limit: None,
    }
}

/// The body of `rusolve explain` on a loaded model (relaxation IIS or LP
/// explanation): the `rusolve.explain/1` object and the text report. Shared by
/// the CLI and the MCP `explain` tool.
pub(crate) fn explain_core(
    mut lp: HighsLp,
    a: &Args,
) -> Result<(Value, String), crate::cli_core::CoreError> {
    use crate::cli_core::CoreError;
    let num_integer = lp
        .integrality_
        .iter()
        .filter(|t| **t != HighsVarType::kContinuous)
        .count();
    let is_mip = num_integer > 0;
    // The LP relaxation: drop integrality (a true MIP IIS: `--mip-iis`, explain_mip.rs).
    lp.integrality_.clear();
    let kind = if is_mip { "relaxation" } else { "lp" };
    let lp_copy = lp.clone();

    let (mut h, status) = solve_for_explanation(lp)
        .map_err(|()| CoreError::internal(format!("the solver failed on {}", a.path)))?;
    let lp = lp_copy;
    let status_text = crate::envelope::result_status(status);

    let mut out = Map::new();
    out.insert("schema_version".into(), json!(SCHEMA_VERSION));
    out.insert(
        "model".into(),
        json!({"path": a.path, "num_rows": lp.num_row_, "num_cols": lp.num_col_,
               "num_integer": num_integer}),
    );
    out.insert("kind".into(), json!(kind));
    out.insert("status".into(), json!(status_text));
    out.insert(
        "note".into(),
        if is_mip {
            json!("model is a MIP: the explanation is for its LP relaxation (integrality dropped); a true MIP IIS is not computed")
        } else {
            Value::Null
        },
    );
    let mut text = String::new();
    if is_mip {
        text.push_str(
            "note: MIP model; explaining its LP relaxation (a true MIP IIS is later work)\n",
        );
    }
    text.push_str(&format!("status ({kind}): {status_text}\n"));
    out.insert("explanation".into(), Value::Null);
    out.insert("verified".into(), Value::Null);

    match status {
        HighsModelStatus::kInfeasible => {
            let (e, t, ok) = explain_infeasible(&mut h, &lp, a, kind);
            out.insert("explanation".into(), e);
            out.insert("verified".into(), json!(ok));
            text.push_str(&t);
        }
        HighsModelStatus::kUnbounded => {
            let (e, t, ok) = explain_unbounded(&mut h, &lp, a);
            out.insert("explanation".into(), e);
            out.insert("verified".into(), json!(ok));
            text.push_str(&t);
        }
        HighsModelStatus::kOptimal => text.push_str("feasible and bounded: nothing to explain\n"),
        _ => text.push_str("no infeasibility or unboundedness explanation for this status\n"),
    }
    Ok((Value::Object(out), text))
}

/// The IIS members HiGHS returned (rows, columns); both empty when `getIis`
/// failed or returned no usable set. Shared with `diagnose`.
pub(crate) fn iis_members(iis: &HighsIis, st: HighsStatus) -> (Vec<Member>, Vec<Member>) {
    let mut rows: Vec<Member> = Vec::new();
    let mut cols: Vec<Member> = Vec::new();
    if st != HighsStatus::kError
        && iis.valid_
        && !(iis.row_index_.is_empty() && iis.col_index_.is_empty())
    {
        for (k, &i) in iis.row_index_.iter().enumerate() {
            let (lower, upper) = decode(iis.row_bound_[k]);
            rows.push(Member {
                index: i as usize,
                lower,
                upper,
                free_col: false,
            });
        }
        for (k, &j) in iis.col_index_.iter().enumerate() {
            let b = iis.col_bound_[k];
            let (lower, upper) = decode(b);
            cols.push(Member {
                index: j as usize,
                lower,
                upper,
                free_col: b == kIisBoundStatusFree,
            });
        }
    }
    (rows, cols)
}

/// Nonzero Farkas multipliers `(row, y)` of a raw dual ray (entries below
/// `1e-12 * max|y|` are noise). Shared with `diagnose`.
pub(crate) fn ray_multipliers(ray: &[f64]) -> Vec<(usize, f64)> {
    let ymax = ray.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    ray.iter()
        .enumerate()
        .filter(|(_, v)| **v != 0.0 && v.abs() > 1e-12 * ymax)
        .map(|(i, v)| (i, *v))
        .collect()
}

fn explain_infeasible(h: &mut Highs, lp: &HighsLp, a: &Args, kind: &str) -> (Value, String, bool) {
    let mut t = String::new();
    let mut e = Map::new();
    e.insert("type".into(), json!("infeasible"));
    let mut ok = true;

    // (a) the IIS
    let mut iis = HighsIis::default();
    let st = h.get_iis(&mut iis);
    let (rows, cols) = iis_members(&iis, st);
    let highs_claim = if iis.status_ == kIisModelStatusIrreducible {
        "irreducible"
    } else if iis.status_ == kIisModelStatusReducible {
        "reducible"
    } else {
        "unknown"
    };
    if rows.is_empty() && cols.is_empty() {
        ok = false;
        e.insert("iis".into(), Value::Null);
        e.insert(
            "iis_error".into(),
            json!(format!(
                "HiGHS returned no usable infeasible set (get_iis returned {st:?}, valid: {}, iis status code: {}); set verification is skipped",
                iis.valid_, iis.status_
            )),
        );
        t.push_str("IIS: none returned by HiGHS\n");
    }
    let mut verification_json = Value::Null;
    if !rows.is_empty() || !cols.is_empty() {
        let mut truncated = false;
        let rj = members_json(lp, &rows, true, a.max_members, &mut truncated);
        let cj = members_json(lp, &cols, false, a.max_members, &mut truncated);
        e.insert(
            "iis".into(),
            json!({
                "kind": kind,
                "highs_claim": highs_claim,
                "num_rows": rows.len(), "num_cols": cols.len(),
                "rows": rj, "columns": cj, "truncated": truncated,
            }),
        );
        t.push_str(&format!(
            "IIS ({highs_claim} per HiGHS): {} rows, {} columns\n",
            rows.len(),
            cols.len()
        ));
        for m in rows.iter().take(a.max_members) {
            t.push_str(&format!(
                "  row {}: {} bound ({} <= . <= {})\n",
                row_name(lp, m.index),
                side_name(m.lower, m.upper),
                lp.row_lower_[m.index],
                lp.row_upper_[m.index]
            ));
        }
        for m in cols.iter().take(a.max_members) {
            if m.lower || m.upper {
                t.push_str(&format!(
                    "  column {}: {} bound ({} <= . <= {})\n",
                    col_name(lp, m.index),
                    side_name(m.lower, m.upper),
                    lp.col_lower_[m.index],
                    lp.col_upper_[m.index]
                ));
            } else {
                t.push_str(&format!(
                    "  column {}: in the set with none of its bounds\n",
                    col_name(lp, m.index)
                ));
            }
        }
        if truncated {
            t.push_str("  (listing truncated by --max-members)\n");
        }
        if a.verify {
            let v = verify_iis(lp, &rows, &cols);
            if !v.verified() {
                ok = false;
            }
            t.push_str(&format!(
                "verification: set infeasible: {}, irreducible: {}\n",
                tri(v.set_infeasible),
                tri(v.irreducible)
            ));
            for f in &v.failures {
                t.push_str(&format!("  FAILED: {f}\n"));
            }
            verification_json = json!({
                "performed": true,
                "set_infeasible": v.set_infeasible,
                "irreducible": v.irreducible,
                "checks": v.checks,
                "failures": v.failures,
                "verified": v.verified(),
            });
        } else {
            ok = false;
            t.push_str("verification skipped (--no-verify)\n");
            verification_json = json!({
                "performed": false, "set_infeasible": null, "irreducible": null,
                "checks": [], "failures": [], "verified": false,
            });
        }
    }
    e.insert("verification".into(), verification_json);

    // (b) the Farkas certificate
    let mut has_ray = false;
    let mut ray = vec![0.0f64; lp.num_row_ as usize];
    let rst = h.get_dual_ray(&mut has_ray, Some(&mut ray));
    let mut farkas_ok = false;
    if rst != HighsStatus::kError && has_ray {
        let y = ray_multipliers(&ray);
        // HiGHS's sign convention is not assumed: try y and -y, keep the proof.
        let pos = farkas_for(lp, &y);
        let neg_y: Vec<(usize, f64)> = y.iter().map(|&(i, v)| (i, -v)).collect();
        let neg = farkas_for(lp, &neg_y);
        let f = if neg.proves && !pos.proves { neg } else { pos };
        farkas_ok = f.proves;
        e.insert("farkas".into(), farkas_json(lp, &f, a.max_members));
        t.push_str(&format!(
            "Farkas certificate: {} nonzero row multipliers\n",
            f.multipliers.len()
        ));
        let mut ms = f.multipliers.clone();
        ms.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
        for &(i, v) in ms.iter().take(a.max_members) {
            t.push_str(&format!("  {} * row {}\n", v, row_name(lp, i)));
        }
        t.push_str(&format!(
            "  rows confine the combination to [{}, {}], column bounds to [{}, {}]\n",
            f.row_lo, f.row_hi, f.col_lo, f.col_hi
        ));
        t.push_str(&format!(
            "  gap {} > {:.1e}: proves 0 >= gap > 0 -> {}\n",
            f.gap,
            f.tolerance,
            if f.proves { "VERIFIED" } else { "NOT verified" }
        ));
    } else {
        e.insert(
            "farkas".into(),
            json!({"available": false, "reason": "HiGHS returned no dual ray"}),
        );
        t.push_str("Farkas certificate: not available\n");
    }
    if !farkas_ok {
        ok = false;
    }
    t.push_str(&format!("verified: {ok}\n"));
    e.insert("verified".into(), json!(ok));
    (Value::Object(e), t, ok)
}

fn tri(v: Option<bool>) -> &'static str {
    match v {
        Some(true) => "yes",
        Some(false) => "NO",
        None => "unknown",
    }
}

fn explain_unbounded(h: &mut Highs, lp: &HighsLp, a: &Args) -> (Value, String, bool) {
    let mut e = Map::new();
    e.insert("type".into(), json!("unbounded"));
    let mut t = String::new();
    let mut has_ray = false;
    let mut d = vec![0.0f64; lp.num_col_ as usize];
    let st = h.get_primal_ray(&mut has_ray, Some(&mut d));
    if st == HighsStatus::kError || !has_ray {
        e.insert(
            "ray".into(),
            json!({"available": false, "reason": "HiGHS returned no primal ray"}),
        );
        e.insert("verified".into(), json!(false));
        t.push_str("primal ray: not available\nverified: false\n");
        return (Value::Object(e), t, false);
    }
    let rc = check_ray(lp, &d);
    let mut nz: Vec<(usize, f64)> = d
        .iter()
        .enumerate()
        .filter(|(_, v)| **v != 0.0)
        .map(|(j, v)| (j, *v))
        .collect();
    nz.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
    let ok = rc.failures.is_empty();
    e.insert(
        "ray".into(),
        json!({
            "available": true,
            "num_nonzero": nz.len(),
            "direction": nz.iter().take(a.max_members).map(|&(j, v)| json!({
                "column": col_name(lp, j), "index": j, "value": v})).collect::<Vec<_>>(),
            "truncated": nz.len() > a.max_members,
            "objective_sense": if lp.sense_ == ObjSense::kMinimize { "minimize" } else { "maximize" },
            "objective_rate": rc.objective_rate,
            "max_row_violation": rc.max_row_violation,
            "max_col_violation": rc.max_col_violation,
            "tolerance": rc.tolerance,
            "failures": rc.failures,
        }),
    );
    e.insert("verified".into(), json!(ok));
    t.push_str(&format!(
        "primal ray: {} nonzero columns (objective changes by {:.6e} per unit step)\n",
        nz.len(),
        rc.objective_rate
    ));
    for &(j, v) in nz.iter().take(a.max_members) {
        t.push_str(&format!("  {} : {}\n", col_name(lp, j), v));
    }
    for f in &rc.failures {
        t.push_str(&format!("  FAILED: {f}\n"));
    }
    t.push_str(&format!("verified: {ok}\n"));
    (Value::Object(e), t, ok)
}
