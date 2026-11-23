//! `rusolve.Model` (docs/design/python-api.md): variables, constraints, objective, solve and
//! in-place changes, over `rusolve_model::ModelBuilder` and the ported `Highs`.
//!
//! Persistence ("Semantics" in the design): until the first `solve()` the model lives only in
//! the builder. The first `solve()` moves it into a ported `Highs` (`pass_model`). After that,
//! new columns and rows stay in the builder until the next `solve()`, which sends them with one
//! `add_cols` and one `add_rows`; changes to entities the solver already has go straight to its
//! change API, which keeps the basis, so the re-solve is warm-started.

use std::sync::{Arc, OnceLock};
use std::time::Instant;

use pyo3::exceptions::{PyIndexError, PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyInt, PyList, PyRange, PyString, PyTuple};
use rusolve_ekk::HighsInt;
use rusolve_ekk::lp_data::hconst::{HighsModelStatus, HighsVarType, ObjSense};
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_run::Highs as Port;
use rusolve_ekk::lp_data::highs_status::HighsStatus;
use rusolve_ekk::prelude::HighsHessian;
use rusolve_model::{Fill, FillBool, ModelBuilder, Names, Sense};

use crate::expr::{Constraint, LinExpr, Var, Vars, new_model_id};
use crate::ingest::{self, Cols, F64s, Flags, IntSlice, Ints, RowKind, build_err};
use crate::solution::Solution;

const INF: f64 = f64::INFINITY;

/// `HighsOptions::recommended()` (the CLI default profile) with the log off.
fn base_options() -> &'static HighsOptions {
    static BASE: OnceLock<HighsOptions> = OnceLock::new();
    BASE.get_or_init(|| {
        let mut o = HighsOptions::recommended();
        o.output_flag = false;
        o.sync_log_options();
        o
    })
}

/// The result status name. Copy of crates/rusolve-cli/src/envelope.rs::result_status (and of the
/// `NodeLimit` rule of `envelope()`), kept here so the bindings do not depend on the CLI.
fn status_name(
    run: HighsStatus,
    s: HighsModelStatus,
    nodes: i64,
    max_nodes: HighsInt,
) -> &'static str {
    use HighsModelStatus::*;
    if run == HighsStatus::kError {
        return "Error";
    }
    match s {
        kOptimal => "Optimal",
        kInfeasible => "Infeasible",
        kUnbounded => "Unbounded",
        kUnboundedOrInfeasible => "InfeasibleOrUnbounded",
        kTimeLimit => "TimeLimit",
        kIterationLimit => "IterationLimit",
        kSolutionLimit if nodes >= max_nodes as i64 => "NodeLimit",
        kSolutionLimit => "SolutionLimit",
        kInterrupt | kHighsInterrupt => "Interrupted",
        kNotset => "NotSolved",
        kMemoryLimit => "MemoryLimit",
        kModelEmpty => "ModelEmpty",
        kObjectiveBound => "ObjectiveBound",
        kObjectiveTarget => "ObjectiveTarget",
        kUnknown => "Unknown",
        kLoadError | kModelError | kPresolveError | kSolveError | kPostsolveError => "Error",
    }
}

#[derive(Clone, Copy, PartialEq)]
enum RowSense {
    Le,
    Ge,
    Eq,
}

fn row_sense(s: &str) -> PyResult<RowSense> {
    Ok(match s {
        "<=" | "=<" | "L" | "l" => RowSense::Le,
        ">=" | "=>" | "G" | "g" => RowSense::Ge,
        "==" | "=" | "E" | "e" => RowSense::Eq,
        _ => {
            return Err(crate::diag::model_err(format!(
                "unknown constraint sense {s:?}: expected \"<=\", \">=\" or \"==\" (or \"L\", \"G\", \"E\")"
            )));
        }
    })
}

#[inline]
fn sense_bounds(s: RowSense, rhs: f64) -> (f64, f64) {
    match s {
        RowSense::Le => (-INF, rhs),
        RowSense::Ge => (rhs, INF),
        RowSense::Eq => (rhs, rhs),
    }
}

fn obj_sense(s: &str) -> PyResult<Sense> {
    match s {
        "min" | "minimize" | "Minimize" | "MIN" => Ok(Sense::Minimize),
        "max" | "maximize" | "Maximize" | "MAX" => Ok(Sense::Maximize),
        _ => Err(crate::diag::model_err(format!(
            "unknown objective sense {s:?}: expected \"min\" or \"max\""
        ))),
    }
}

fn port_sense(s: Sense) -> ObjSense {
    match s {
        Sense::Minimize => ObjSense::kMinimize,
        Sense::Maximize => ObjSense::kMaximize,
    }
}

/// Row senses: one for all, or one per row.
enum Senses {
    One(RowSense),
    Each(Vec<RowSense>),
}

impl Senses {
    fn read(o: &Bound<'_, PyAny>) -> PyResult<Senses> {
        if let Ok(s) = o.cast::<PyString>() {
            return Ok(Senses::One(row_sense(s.to_str()?)?));
        }
        let mut v = Vec::new();
        for s in o.try_iter()? {
            let s = s?;
            let s = s
                .cast::<PyString>()
                .map_err(|_| PyTypeError::new_err("sense must be a str or a sequence of str"))?;
            v.push(row_sense(s.to_str()?)?);
        }
        Ok(Senses::Each(v))
    }
    #[inline]
    fn get(&self, k: usize) -> PyResult<RowSense> {
        match self {
            Senses::One(s) => Ok(*s),
            Senses::Each(v) => v.get(k).copied().ok_or_else(|| {
                crate::diag::model_err(format!("sense: only {} values for more rows", v.len()))
            }),
        }
    }
}

/// How the bounds of a batch of rows are given.
enum BoundSpec {
    /// `sense` and `rhs`.
    Sense(Senses, F64s),
    /// `lo` and/or `hi`.
    Range(Option<F64s>, Option<F64s>),
    /// Nothing: every row must be a Constraint.
    None,
}

impl BoundSpec {
    fn read(
        sense: Option<&Bound<'_, PyAny>>,
        rhs: Option<&Bound<'_, PyAny>>,
        lo: Option<&Bound<'_, PyAny>>,
        hi: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<BoundSpec> {
        let lo = ingest::opt_f64s(lo, "lo")?;
        let hi = ingest::opt_f64s(hi, "hi")?;
        let sense = sense.filter(|s| !s.is_none());
        let rhs = rhs.filter(|s| !s.is_none());
        match (sense, rhs) {
            (Some(s), Some(r)) => {
                if lo.is_some() || hi.is_some() {
                    return Err(crate::diag::model_err(
                        "give either sense and rhs, or lo and/or hi, not both",
                    ));
                }
                Ok(BoundSpec::Sense(Senses::read(s)?, F64s::read(r, "rhs")?))
            }
            (Some(_), None) => Err(crate::diag::model_err("a constraint sense needs an rhs")),
            (None, Some(_)) => Err(crate::diag::model_err("an rhs needs a constraint sense")),
            (None, None) if lo.is_none() && hi.is_none() => Ok(BoundSpec::None),
            (None, None) => Ok(BoundSpec::Range(lo, hi)),
        }
    }

    /// The bounds of row `k` of the batch, for a plain row with constant `c`.
    #[inline]
    fn row(&self, k: usize, c: f64) -> PyResult<(f64, f64)> {
        let (lo, hi) = match self {
            BoundSpec::Sense(s, r) => {
                let r = per_item_checked(r, k, "rhs")?;
                sense_bounds(s.get(k)?, r)
            }
            BoundSpec::Range(lo, hi) => (
                lo.as_ref()
                    .map_or(Ok(-INF), |v| per_item_checked(v, k, "lo"))?,
                hi.as_ref()
                    .map_or(Ok(INF), |v| per_item_checked(v, k, "hi"))?,
            ),
            BoundSpec::None => {
                return Err(crate::diag::model_err(
                    "this row needs bounds: give sense and rhs (e.g. \"<=\", 3.0), or lo/hi, or \
                     pass a Constraint (e.g. x + y <= 3)",
                ));
            }
        };
        Ok((lo - c, hi - c))
    }

    fn is_none(&self) -> bool {
        matches!(self, BoundSpec::None)
    }

    /// Per-row lower/upper fills for `n` matrix rows.
    fn fills(&self, n: usize) -> PyResult<(FillOwned, FillOwned)> {
        match self {
            BoundSpec::Sense(Senses::One(s), r) => {
                r.per_item(n, "rhs")?;
                Ok(match s {
                    RowSense::Le => (FillOwned::Scalar(-INF), FillOwned::Ref),
                    RowSense::Ge => (FillOwned::Ref, FillOwned::Scalar(INF)),
                    RowSense::Eq => (FillOwned::Ref, FillOwned::Ref),
                })
            }
            BoundSpec::None => Err(crate::diag::model_err(
                "add_matrix needs bounds: give sense and rhs, or lo and/or hi",
            )),
            _ => {
                let mut lo = Vec::with_capacity(n);
                let mut hi = Vec::with_capacity(n);
                for k in 0..n {
                    let (l, h) = self.row(k, 0.0)?;
                    lo.push(l);
                    hi.push(h);
                }
                if let BoundSpec::Sense(Senses::Each(v), _) = self
                    && v.len() != n
                {
                    return Err(crate::diag::model_err(format!(
                        "sense: expected {n} values, got {}",
                        v.len()
                    )));
                }
                Ok((FillOwned::Owned(lo), FillOwned::Owned(hi)))
            }
        }
    }

    /// Check that per-row inputs have exactly `n` values (after a batch of `n` rows).
    fn check_len(&self, n: usize) -> PyResult<()> {
        let check = |v: &F64s, what: &str| match v.len() {
            Some(m) if m != n => Err(crate::diag::model_err(format!(
                "{what}: expected {n} values, got {m}"
            ))),
            _ => Ok(()),
        };
        match self {
            BoundSpec::Sense(s, r) => {
                check(r, "rhs")?;
                if let Senses::Each(v) = s
                    && v.len() != n
                {
                    return Err(crate::diag::model_err(format!(
                        "sense: expected {n} values, got {}",
                        v.len()
                    )));
                }
            }
            BoundSpec::Range(lo, hi) => {
                if let Some(v) = lo {
                    check(v, "lo")?;
                }
                if let Some(v) = hi {
                    check(v, "hi")?;
                }
            }
            BoundSpec::None => {}
        }
        Ok(())
    }
}

/// A bound fill of a matrix batch: a scalar, the `rhs` input itself, or computed values.
enum FillOwned {
    Scalar(f64),
    Ref,
    Owned(Vec<f64>),
}

impl FillOwned {
    fn fill<'a>(&'a self, spec: &'a BoundSpec) -> Fill<'a> {
        match self {
            FillOwned::Scalar(c) => Fill::Scalar(*c),
            FillOwned::Owned(v) => Fill::Slice(v),
            FillOwned::Ref => match spec {
                BoundSpec::Sense(_, r) => r.fill(),
                _ => unreachable!("FillOwned::Ref only for a sense/rhs spec"),
            },
        }
    }
}

#[inline]
fn per_item_checked(v: &F64s, k: usize, what: &str) -> PyResult<f64> {
    match v {
        F64s::Scalar(c) => Ok(*c),
        _ => v.slice().get(k).copied().ok_or_else(|| {
            crate::diag::model_err(format!(
                "{what}: only {} values for more rows",
                v.slice().len()
            ))
        }),
    }
}

fn names_arg(o: Option<&Bound<'_, PyAny>>, n: usize, what: &str) -> PyResult<Option<Names>> {
    let Some(o) = o.filter(|o| !o.is_none()) else {
        return Ok(None);
    };
    if let Ok(s) = o.cast::<PyString>() {
        return Ok(Some(Names::Prefix(s.to_str()?.to_string())));
    }
    let v: Vec<String> = o.extract().map_err(|_| {
        PyTypeError::new_err(format!("{what}: expected a str prefix or a list of str"))
    })?;
    if v.len() != n {
        return Err(crate::diag::model_err(format!(
            "{what}: expected {n} names, got {}",
            v.len()
        )));
    }
    Ok(Some(Names::List(v)))
}

/// Sort `(index, value...)` entries by index and keep the last value of a repeated index (the
/// port's `*_set` changes want strictly increasing sets, as HiGHS's).
fn sorted_unique<T: Copy>(mut v: Vec<(HighsInt, T)>) -> (Vec<HighsInt>, Vec<T>) {
    v.sort_by_key(|e| e.0); // stable: later entries stay after earlier ones
    let mut set: Vec<HighsInt> = Vec::with_capacity(v.len());
    let mut vals: Vec<T> = Vec::with_capacity(v.len());
    for (j, x) in v {
        if set.last() == Some(&j) {
            *vals.last_mut().expect("same length") = x;
        } else {
            set.push(j);
            vals.push(x);
        }
    }
    (set, vals)
}

fn port_err(what: &str, s: HighsStatus) -> PyResult<()> {
    if s == HighsStatus::kError {
        Err(PyRuntimeError::new_err(format!(
            "{what} failed in the solver"
        )))
    } else {
        Ok(())
    }
}

/// Row indices: an int, a `range`, or a sequence/buffer of ints, checked against `0..limit`.
fn row_list(o: &Bound<'_, PyAny>, limit: usize) -> PyResult<Vec<usize>> {
    let v: Vec<i64> = if o.is_instance_of::<PyInt>() {
        vec![o.extract()?]
    } else {
        let ints = Ints::read(o, "rows", None)?;
        let s = ints.slice();
        (0..s.len()).map(|k| s.get(k)).collect()
    };
    v.into_iter()
        .map(|i| {
            if i < 0 || i as u64 >= limit as u64 {
                Err(crate::diag::index_err(format!(
                    "row {i} is out of range (there are {limit} rows)"
                )))
            } else {
                Ok(i as usize)
            }
        })
        .collect()
}

/// An optimisation model. See docs/design/python-api.md.
///
/// >>> import rusolve as rs
/// >>> m = rs.Model()                       # minimise by default; Model(sense="max")
/// >>> x = m.add_vars(3, binary=True)       # a block of 3 binary variables
/// >>> r = m.add_constr(x[0] + x[1] + x[2] <= 2)
/// >>> m.maximize(5 * x[0] + 4 * x[1] + 3 * x[2])
/// >>> sol = m.solve()
/// >>> sol.status, sol.objective
/// ('Optimal', 9.0)
#[pyclass(module = "rusolve")]
pub struct Model {
    id: u64,
    b: ModelBuilder,
    /// The solver, from the first `solve()` on.
    port: Option<Box<Port>>,
    sense: Sense,
    /// Coefficients of (synced row, unsynced column), applied after the column is added.
    pending_coefs: Vec<(HighsInt, HighsInt, f64)>,
    /// The last solve set options; the next one starts again from the base options.
    options_dirty: bool,
    // --- convex QP: the objective's Hessian. -------------------------------------------------
    /// `Some(entries)`: the full symmetric square Hessian of `objective = c'x + x'Qx/2` as
    /// `(row, col, value)` (both triangles stored), set by `set_quadratic_objective`.
    quad: Option<Vec<(usize, usize, f64)>>,
    /// `quad` changed since the last `pass_hessian` into the solver.
    quad_dirty: bool,
    /// The convexity report of the last analysis (`quad_info()`).
    quad_info: Option<crate::quad::QuadInfo>,
    // --- What `Solution.to_json()` and the diagnostics need. -------------------------
    /// The name given to `Model(name=...)`; `input.path` of the result envelope.
    name: Option<String>,
    /// Column names of the solver's model, shared by the Solutions; reset when columns arrive.
    names: Option<Arc<Vec<String>>>,
    /// The settings the last `solve()` call changed, as the envelope's `run.options`.
    last_options: Vec<(String, String)>,
    // --- The Python line that added each block of rows (off by default). -------------
    /// `capture_source`: record the calling line of every `add_constr`/`add_constrs`/`add_matrix`.
    capture_source: bool,
    /// `(first row, end row, path, line)` in row order; looked up by `source_of`.
    row_sources: Vec<(usize, usize, Arc<str>, usize)>,
    /// Counters of the last `solve(callback=...)`, for `callback_stats`.
    callback_stats: Option<Py<PyDict>>,
}

impl Model {
    /// Remember the calling Python line for rows `first..end` (only when `capture_source`).
    fn record_source(&mut self, first: usize, end: usize) {
        if !self.capture_source || end <= first {
            return;
        }
        let at = Python::attach(crate::diag::caller_location);
        if let Some(rusolve_core::diagnostic::Location::Python { path, line, .. }) = at {
            self.row_sources.push((first, end, Arc::from(path), line));
        }
    }

    /// Add `source` (see `source_of`) to each row member of an explain result.
    fn annotate_sources(&self, py: Python<'_>, v: &Bound<'_, PyAny>) -> PyResult<()> {
        let Ok(iis) = v
            .get_item("explanation")
            .and_then(|e| e.get_item("iis"))
            .and_then(|i| i.get_item("rows"))
        else {
            return Ok(());
        };
        for m in iis.try_iter()? {
            let m = m?;
            let row: usize = m.get_item("index")?.extract()?;
            if let Some(d) = self.source_dict(py, row)? {
                m.set_item("source", d)?;
            }
        }
        Ok(())
    }

    fn source_dict<'py>(
        &self,
        py: Python<'py>,
        row: usize,
    ) -> PyResult<Option<Bound<'py, PyDict>>> {
        let Some((_, _, path, line)) = self.row_sources.iter().find(|s| s.0 <= row && row < s.1)
        else {
            return Ok(None);
        };
        let d = PyDict::new(py);
        d.set_item("path", &**path)?;
        d.set_item("line", line)?;
        let src: String = py
            .import("linecache")?
            .call_method1("getline", (&**path, *line))?
            .extract()?;
        d.set_item("source_line", src.trim())?;
        Ok(Some(d))
    }

    fn synced_cols(&self) -> usize {
        self.b.synced_cols()
    }
    fn synced_rows(&self) -> usize {
        self.b.synced_rows()
    }

    fn check_var(&self, v: &Var) -> PyResult<()> {
        if v.model != self.id {
            return Err(crate::diag::model_err(
                "cannot combine variables of two different Models: this Var belongs to another Model",
            ));
        }
        Ok(())
    }

    fn col_index(&self, o: &Bound<'_, PyAny>) -> PyResult<usize> {
        let j = ingest::col_of(o, self.id)?;
        if j < 0 || j as u64 >= self.b.num_cols() as u64 {
            return Err(crate::diag::index_err(format!(
                "variable {j} is out of range (there are {} variables)",
                self.b.num_cols()
            )));
        }
        Ok(j as usize)
    }

    /// Add one row (any row form) with the bounds spec; returns its index.
    fn add_one(
        &mut self,
        row: &Bound<'_, PyAny>,
        spec: &BoundSpec,
        k: usize,
        name: Option<&str>,
    ) -> PyResult<usize> {
        let kind = ingest::push_row(&mut self.b, row, self.id).map_err(|e| with_name(e, name))?;
        let bounds = match kind {
            RowKind::Bounded(lo, hi) => {
                if !spec.is_none() {
                    self.b.abort_row();
                    return Err(crate::diag::model_err(
                        "this row is a Constraint, which has its own bounds: do not give sense, \
                         rhs, lo or hi as well",
                    ));
                }
                (lo, hi)
            }
            RowKind::Plain(c) => match spec.row(k, c) {
                Ok(b) => b,
                Err(e) => {
                    self.b.abort_row();
                    return Err(e);
                }
            },
        };
        self.b
            .finish_row(bounds.0, bounds.1)
            .map_err(|e| build_err(e, name))
    }

    /// Move pending columns and rows to the solver (creating it at the first call).
    fn sync(&mut self) -> PyResult<()> {
        let Some(port) = self.port.as_mut() else {
            let lp = self.b.take().into_lp();
            let mut port = Box::new(Port::with_options(base_options().clone()));
            port_err("pass_model", port.pass_model(lp))?;
            self.port = Some(port);
            self.names = None;
            return Ok(());
        };
        let pending =
            self.b.num_cols() > self.b.synced_cols() || self.b.num_rows() > self.b.synced_rows();
        if pending {
            self.names = None;
            let batch = self.b.take();
            let c = &batch.cols;
            let nc = c.cost.len();
            if nc > 0 {
                let starts = vec![0 as HighsInt; nc];
                port_err(
                    "add_cols",
                    port.add_cols(
                        nc as HighsInt,
                        &c.cost,
                        &c.lower,
                        &c.upper,
                        0,
                        &starts,
                        &[],
                        &[],
                    ),
                )?;
                let (from, to) = (c.first as HighsInt, (c.first + nc - 1) as HighsInt);
                if !c.integrality.is_empty() {
                    port_err(
                        "change_cols_integrality",
                        port.change_cols_integrality_interval(from, to, &c.integrality),
                    )?;
                }
                for (k, name) in c.names.iter().enumerate() {
                    port_err(
                        "pass_col_name",
                        port.pass_col_name(from + k as HighsInt, name),
                    )?;
                }
            }
            let r = &batch.rows;
            let nr = r.lower.len();
            if nr > 0 {
                port_err(
                    "add_rows",
                    port.add_rows(
                        nr as HighsInt,
                        &r.lower,
                        &r.upper,
                        r.index.len() as HighsInt,
                        &r.start[..nr],
                        &r.index,
                        &r.value,
                    ),
                )?;
                for (k, name) in r.names.iter().enumerate() {
                    port_err(
                        "pass_row_name",
                        port.pass_row_name((r.first + k) as HighsInt, name),
                    )?;
                }
            }
        }
        for (i, j, v) in std::mem::take(&mut self.pending_coefs) {
            port_err("change_coeff", port.change_coeff(i, j, v))?;
        }
        Ok(())
    }

    /// Apply the solve options to the solver (from the base options when the previous solve
    /// set some).
    #[allow(clippy::too_many_arguments)]
    fn apply_options(
        &mut self,
        time_limit: Option<f64>,
        mip_gap: Option<f64>,
        threads: Option<i32>,
        presolve: Option<&str>,
        seed: Option<i32>,
        verbose: bool,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<()> {
        let port = self.port.as_mut().expect("synced");
        if self.options_dirty {
            port.pass_options(base_options());
            self.options_dirty = false;
        }
        let any = time_limit.is_some()
            || mip_gap.is_some()
            || threads.is_some()
            || presolve.is_some()
            || seed.is_some()
            || verbose
            || options.is_some_and(|d| !d.is_empty());
        if !any {
            return Ok(());
        }
        self.options_dirty = true;
        let set = |s: HighsStatus, name: &str, v: String| {
            if s == HighsStatus::kError {
                Err(crate::diag::model_err(format!(
                    "invalid value {v} for option {name:?}"
                )))
            } else {
                Ok(())
            }
        };
        if let Some(t) = time_limit {
            set(
                port.set_option_value_double("time_limit", t),
                "time_limit",
                t.to_string(),
            )?;
        }
        if let Some(g) = mip_gap {
            set(
                port.set_option_value_double("mip_rel_gap", g),
                "mip_rel_gap",
                g.to_string(),
            )?;
        }
        if let Some(t) = threads {
            set(
                port.set_option_value_int("threads", t),
                "threads",
                t.to_string(),
            )?;
        }
        if let Some(v) = presolve {
            set(
                port.set_option_value_str("presolve", v),
                "presolve",
                format!("{v:?}"),
            )?;
        }
        if let Some(s) = seed {
            set(
                port.set_option_value_int("random_seed", s),
                "random_seed",
                s.to_string(),
            )?;
        }
        if verbose {
            port.set_option_value_bool("output_flag", true);
        }
        if let Some(d) = options {
            for (k, v) in d.iter() {
                let name: String = k
                    .extract()
                    .map_err(|_| PyTypeError::new_err("option names must be str"))?;
                let s = crate::set_option_any(port, &name, &v)?;
                if s == HighsStatus::kError {
                    return Err(crate::diag::model_err(format!(
                        "unknown option {name:?}, or invalid value {} for it",
                        v.repr().map(|r| r.to_string()).unwrap_or_default()
                    )));
                }
            }
        }
        Ok(())
    }

    fn solution(&mut self, run: HighsStatus, time: f64) -> Solution {
        let port = self.port.as_ref().expect("synced");
        let info = port.get_info();
        let ms = port.get_model_status();
        let sol = port.get_solution();
        let lp = port.get_lp();
        let is_mip = lp
            .integrality_
            .iter()
            .any(|t| *t != HighsVarType::kContinuous);
        let status = status_name(run, ms, info.mip_node_count, port.options_.mip_max_nodes);
        let ok = run != HighsStatus::kError
            && sol.value_valid
            && sol.col_value.len() == lp.num_col_ as usize
            && match ms {
                HighsModelStatus::kOptimal => true,
                HighsModelStatus::kInfeasible
                | HighsModelStatus::kUnboundedOrInfeasible
                | HighsModelStatus::kUnbounded => false,
                _ => {
                    info.primal_solution_status
                        == rusolve_ekk::lp_data::hconst::kSolutionStatusFeasible
                }
            };
        let fin = |v: f64| v.is_finite().then_some(v);
        let objective = if ok {
            Some(info.objective_function_value)
        } else {
            None
        };
        let (bound, gap) = if is_mip {
            (
                fin(info.mip_dual_bound),
                if ok { fin(info.mip_gap) } else { None },
            )
        } else if ms == HighsModelStatus::kOptimal {
            (objective, Some(0.0))
        } else {
            (None, None)
        };
        let duals_ok = ok && !is_mip && sol.dual_valid;
        // The rusolve.result/1 envelope is built here, while the solver still holds this
        // solve; the column names are shared between solves until a column is added.
        let names = self
            .names
            .get_or_insert_with(|| Arc::new(lp.col_names_.clone()))
            .clone();
        let mut env = rusolve_cli::api::env_parts(
            port,
            run,
            time,
            self.name
                .as_deref()
                .unwrap_or(rusolve_cli::api::DEFAULT_NAME),
            &self.last_options,
        );
        if self.quad.is_some() {
            env.set_convex_qp();
        }
        Solution {
            env: Some(env),
            names,
            status,
            objective,
            bound,
            gap,
            x: ok.then(|| sol.col_value.clone()),
            row_activity: ok.then(|| sol.row_value.clone()),
            duals: duals_ok.then(|| sol.row_dual.clone()),
            reduced_costs: duals_ok.then(|| sol.col_dual.clone()),
            nodes: if is_mip {
                info.mip_node_count.max(0)
            } else {
                0
            },
            iterations: (info.simplex_iteration_count.max(0) + info.ipm_iteration_count.max(0))
                as i64,
            time,
            model: self.id,
            lp: matches!(
                ms,
                HighsModelStatus::kInfeasible
                    | HighsModelStatus::kUnbounded
                    | HighsModelStatus::kUnboundedOrInfeasible
            )
            .then(|| Arc::new(lp.clone())),
            quad: self.quad.is_some(),
            cache: [
                PyOnceLock::new(),
                PyOnceLock::new(),
                PyOnceLock::new(),
                PyOnceLock::new(),
            ],
        }
    }

    /// Dense objective (length num_cols) and offset from any objective form.
    fn dense_objective(&self, obj: &Bound<'_, PyAny>) -> PyResult<(Vec<f64>, f64)> {
        let n = self.b.num_cols();
        let mut c = vec![0.0; n];
        let mut offset = 0.0;
        let mut add = |j: i64, v: f64| -> PyResult<()> {
            if j < 0 || j as u64 >= n as u64 {
                return Err(crate::diag::index_err(format!(
                    "objective: variable {j} is out of range (there are {n} variables)"
                )));
            }
            c[j as usize] += v;
            Ok(())
        };
        if obj.is_none() {
        } else if let Ok(d) = obj.cast::<PyDict>() {
            for (k, v) in d.iter() {
                add(ingest::col_of(&k, self.id)?, ingest::coef(&v)?)?;
            }
        } else if let Ok(e) = obj.cast_exact::<LinExpr>() {
            let e = e.try_borrow()?;
            if e.model != 0 && e.model != self.id {
                return Err(crate::diag::model_err(crate::expr::MIXED_MODELS));
            }
            e.terms
                .with(|t| t.iter().try_for_each(|&(j, v)| add(j as i64, v)))?;
            offset = e.constant;
        } else if let Ok(v) = obj.cast_exact::<Var>() {
            let v = v.get();
            self.check_var(v)?;
            add(v.index as i64, 1.0)?;
        } else if obj.is_instance_of::<Constraint>() {
            return Err(PyTypeError::new_err(
                "the objective is an expression, not a Constraint",
            ));
        } else if let Ok(t) = obj.cast_exact::<PyTuple>()
            && t.len() == 2
        {
            let idx = Ints::read(&t.get_item(0)?, "objective indices", Some(self.id))?;
            let val = F64s::read_array(&t.get_item(1)?, "objective values")?;
            let (is, vs) = (idx.slice(), val.slice());
            if is.len() != vs.len() {
                return Err(crate::diag::model_err(format!(
                    "objective (indices, values): {} indices but {} values",
                    is.len(),
                    vs.len()
                )));
            }
            for (k, &v) in vs.iter().enumerate() {
                add(is.get(k), v)?;
            }
        } else {
            let v = F64s::read_array(obj, "objective")?;
            let s = v.slice();
            if s.len() != n {
                return Err(crate::diag::model_err(format!(
                    "objective: a dense objective needs {n} values (one per variable), got {}",
                    s.len()
                )));
            }
            c.copy_from_slice(s);
        }
        for (j, &v) in c.iter().enumerate() {
            if !v.is_finite() {
                return Err(crate::diag::model_err(format!(
                    "objective: the cost of variable {j} is {v}; costs must be finite"
                )));
            }
        }
        if !offset.is_finite() {
            return Err(crate::diag::model_err(
                "objective: the constant must be finite",
            ));
        }
        Ok((c, offset))
    }

    /// Hand a changed Hessian to the solver (`Highs::passHessian`). A model that never had one
    /// does no QP work at all.
    fn sync_quad(&mut self) -> PyResult<()> {
        let port = self.port.as_mut().expect("synced");
        if let Some(entries) = &self.quad {
            // The MIQP refusal on every solve, the convexity analysis when Q or the
            // objective sense changed.
            let maximize = self.sense == Sense::Maximize;
            let checked = crate::quad::assess(port.get_lp(), entries, maximize, self.quad_dirty)?;
            if let Some((shifted, info)) = checked {
                let h = crate::quad::highs_hessian(&shifted, self.b.num_cols());
                port_err("pass_hessian", port.pass_hessian(h))?;
                self.quad_info = Some(info);
                self.quad_dirty = false;
            }
        } else if self.quad_dirty {
            port_err("pass_hessian", port.pass_hessian(HighsHessian::default()))?;
            self.quad_info = None;
            self.quad_dirty = false;
        }
        Ok(())
    }

    fn set_sense_inner(&mut self, s: Sense) -> PyResult<()> {
        // the convexity verdict depends on the sense (PSD for min, NSD for max)
        self.quad_dirty |= self.quad.is_some() && s != self.sense;
        self.sense = s;
        self.b.set_sense(s);
        if let Some(p) = self.port.as_mut() {
            port_err(
                "change_objective_sense",
                p.change_objective_sense(port_sense(s)),
            )?;
        }
        Ok(())
    }
}

impl Model {
    /// `input.path` / `model.path` of the diagnostics: the model's name, or a placeholder.
    fn display_name(&self) -> String {
        self.name
            .clone()
            .unwrap_or_else(|| rusolve_cli::api::DEFAULT_NAME.to_string())
    }

    fn refuse_quad(&self, what: &str) -> PyResult<()> {
        if self.quad.is_some() {
            return Err(crate::diag::model_err(format!(
                "{what}() does not support a model with a quadratic objective yet"
            )));
        }
        Ok(())
    }

    /// One value per column from a Solution, a sequence or array, or a dict keyed by Var, index
    /// or column name.
    fn assignment(&self, o: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
        let n = self.b.num_cols();
        if let Ok(sol) = o.cast::<Solution>() {
            let sol = sol.get();
            if sol.model != self.id {
                return Err(crate::diag::model_err(
                    "this Solution belongs to another Model",
                ));
            }
            return sol.x.clone().ok_or_else(|| {
                crate::diag::model_err(format!(
                    "the Solution has no primal values (status {})",
                    sol.status
                ))
            });
        }
        let values: Vec<f64> = if let Ok(d) = o.cast::<PyDict>() {
            let mut by_name: Vec<(String, f64)> = Vec::new();
            let mut by_col: Vec<Option<f64>> = vec![None; n];
            for (k, v) in d.iter() {
                let v: f64 = v.extract()?;
                if let Ok(name) = k.extract::<String>() {
                    by_name.push((name, v));
                } else {
                    let j = self.col_index(&k)?;
                    if by_col[j].replace(v).is_some() {
                        return Err(crate::diag::model_err(format!(
                            "variable {j} appears twice"
                        )));
                    }
                }
            }
            if !by_name.is_empty() && by_col.iter().any(Option::is_some) {
                return Err(PyTypeError::new_err(
                    "the dict mixes column names with Var or int keys",
                ));
            }
            if !by_name.is_empty() {
                let lp = &self.port.as_ref().expect("synced").model_lp_;
                return rusolve_cli::api::columns_by_name(lp, &by_name)
                    .map_err(crate::diag::model_err);
            }
            let missing = by_col.iter().filter(|v| v.is_none()).count();
            if missing > 0 {
                return Err(crate::diag::model_err(format!(
                    "the assignment misses {missing} of {n} variables"
                )));
            }
            by_col.into_iter().flatten().collect()
        } else {
            F64s::read_array(o, "solution")?.slice().to_vec()
        };
        if values.len() != n {
            return Err(crate::diag::model_err(format!(
                "expected {n} values (one per variable), got {}",
                values.len()
            )));
        }
        Ok(values)
    }
}

/// `json.loads` of a JSON value: the dict a diagnostic returns.
fn json_to_py(py: Python<'_>, v: &serde_json::Value) -> PyResult<Py<PyAny>> {
    let text = serde_json::to_string(v).map_err(|e| crate::diag::model_err(e.to_string()))?;
    Ok(py.import("json")?.call_method1("loads", (text,))?.unbind())
}

/// The `run.options` pairs of a solve call: the named arguments that were given, then the
/// `options=` dict, in order.
fn run_options(
    time_limit: Option<f64>,
    mip_gap: Option<f64>,
    threads: Option<i32>,
    presolve: Option<&str>,
    seed: Option<i32>,
    verbose: bool,
    options: Option<&Bound<'_, PyDict>>,
) -> Vec<(String, String)> {
    let mut v = Vec::new();
    if let Some(t) = time_limit {
        v.push(("time_limit".to_string(), t.to_string()));
    }
    if let Some(g) = mip_gap {
        v.push(("mip_rel_gap".to_string(), g.to_string()));
    }
    if let Some(t) = threads {
        v.push(("threads".to_string(), t.to_string()));
    }
    if let Some(p) = presolve {
        v.push(("presolve".to_string(), p.to_string()));
    }
    if let Some(s) = seed {
        v.push(("random_seed".to_string(), s.to_string()));
    }
    if verbose {
        v.push(("output_flag".to_string(), "true".to_string()));
    }
    if let Some(d) = options {
        for (k, val) in d.iter() {
            v.push((
                k.to_string(),
                val.str().map(|s| s.to_string()).unwrap_or_default(),
            ));
        }
    }
    v
}

fn with_name(e: PyErr, name: Option<&str>) -> PyErr {
    match name {
        None => e,
        Some(n) => Python::attach(|py| {
            let msg = format!("{} (constraint {n:?})", e.value(py));
            if e.is_instance_of::<PyIndexError>(py) {
                crate::diag::index_err(msg)
            } else if e.is_instance_of::<PyValueError>(py) {
                crate::diag::model_err(msg)
            } else {
                PyErr::from_type(e.get_type(py), msg)
            }
        }),
    }
}

#[pymethods]
impl Model {
    #[new]
    #[pyo3(signature = (sense="min", name=None))]
    fn new(sense: &str, name: Option<String>) -> PyResult<Self> {
        let s = obj_sense(sense)?;
        let mut b = ModelBuilder::new();
        b.set_sense(s);
        if let Some(n) = &name {
            b.set_name(n.clone());
        }
        Ok(Model {
            id: new_model_id(),
            b,
            port: None,
            sense: s,
            pending_coefs: Vec::new(),
            options_dirty: false,
            quad: None,
            quad_dirty: false,
            quad_info: None,
            name,
            names: None,
            last_options: Vec::new(),
            capture_source: false,
            row_sources: Vec::new(),
            callback_stats: None,
        })
    }

    // ------------------------------------------------------------ variables

    /// Add `n` variables; returns a `Vars` block. `lb`, `ub`, `cost`: scalar or one per
    /// variable (numpy arrays are read in place); `integer`: bool or one per variable;
    /// `binary=True`: integer with default bounds [0, 1]; `names`: a list of `n` str or a prefix `p` giving
    /// `p[0]`, `p[1]`, ...
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model()
    /// >>> x = m.add_vars(3, lb=0, ub=[1, 2, 3], names="x")
    /// >>> len(x), m.num_vars, m.var_name(x[2])
    /// (3, 3, 'x[2]')
    #[pyo3(signature = (n, lb=None, ub=None, cost=None, integer=None, binary=false, names=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_vars(
        &mut self,
        n: usize,
        lb: Option<&Bound<'_, PyAny>>,
        ub: Option<&Bound<'_, PyAny>>,
        cost: Option<&Bound<'_, PyAny>>,
        integer: Option<&Bound<'_, PyAny>>,
        binary: bool,
        names: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Vars> {
        if self.b.num_cols() + n > i32::MAX as usize {
            return Err(crate::diag::model_err(
                "too many variables (at most 2^31 - 1)",
            ));
        }
        let names = names_arg(names, n, "names")?;
        let lb = ingest::opt_f64s(lb, "lb")?;
        let ub = ingest::opt_f64s(ub, "ub")?;
        let cost = ingest::opt_f64s(cost, "cost")?.unwrap_or(F64s::Scalar(0.0));
        let flags = integer
            .filter(|o| !o.is_none())
            .map(|o| Flags::read(o, "integer"))
            .transpose()?;
        // binary: integer with default bounds [0, 1]; an explicit lb/ub still applies
        let lb = lb.unwrap_or(F64s::Scalar(0.0));
        let ub = ub.unwrap_or(F64s::Scalar(if binary { 1.0 } else { INF }));
        lb.per_item(n, "lb")?;
        ub.per_item(n, "ub")?;
        cost.per_item(n, "cost")?;
        let integer = match (&flags, binary) {
            (_, true) => FillBool::Scalar(true),
            (None, false) => FillBool::Scalar(false),
            (Some(Flags::Scalar(f)), false) => FillBool::Scalar(*f),
            (Some(Flags::Each(v)), false) => {
                if v.len() != n {
                    return Err(crate::diag::model_err(format!(
                        "integer: expected {n} values, got {}",
                        v.len()
                    )));
                }
                FillBool::Slice(v)
            }
        };
        let r = self
            .b
            .add_cols(n, cost.fill(), lb.fill(), ub.fill(), integer)
            .map_err(|e| build_err(e, None))?;
        if let Some(names) = names {
            self.b
                .set_col_names(r.clone(), names)
                .map_err(|e| build_err(e, None))?;
        }
        Ok(Vars {
            start: r.start as i64,
            step: 1,
            len: n,
            model: self.id,
        })
    }

    /// Add one variable; returns a `Var`.
    #[pyo3(signature = (lb=None, ub=None, cost=0.0, integer=false, binary=false, name=None))]
    fn add_var(
        &mut self,
        lb: Option<f64>,
        ub: Option<f64>,
        cost: f64,
        integer: bool,
        binary: bool,
        name: Option<String>,
    ) -> PyResult<Var> {
        let lb = lb.unwrap_or(0.0);
        let ub = ub.unwrap_or(if binary { 1.0 } else { INF });
        if self.b.num_cols() >= i32::MAX as usize {
            return Err(crate::diag::model_err(
                "too many variables (at most 2^31 - 1)",
            ));
        }
        let r = self
            .b
            .add_cols(
                1,
                Fill::Scalar(cost),
                Fill::Scalar(lb),
                Fill::Scalar(ub),
                FillBool::Scalar(integer || binary),
            )
            .map_err(|e| build_err(e, None))?;
        if let Some(name) = name {
            self.b
                .set_col_names(r.clone(), Names::List(vec![name]))
                .map_err(|e| build_err(e, None))?;
        }
        Ok(Var::new(r.start as i32, self.id))
    }

    // ------------------------------------------------------------ constraints

    /// Add one row; returns its index. `row`: a dict `{var: coef}`, an `(indices, values)`
    /// pair, a LinExpr or Var (then give `sense` and `rhs`, or `lo`/`hi`), or a Constraint
    /// (`x + y <= 3`, which carries its bounds).
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_vars(2)
    /// >>> m.add_constr(x[0] + 2 * x[1] <= 10, name="cap")        # a Constraint
    /// 0
    /// >>> m.add_constr({x[0]: 1.0, x[1]: 1.0}, ">=", 1.0)         # a dict row, with sense and rhs
    /// 1
    #[pyo3(signature = (row, sense=None, rhs=None, name=None, *, lo=None, hi=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_constr(
        &mut self,
        row: &Bound<'_, PyAny>,
        sense: Option<&Bound<'_, PyAny>>,
        rhs: Option<&Bound<'_, PyAny>>,
        name: Option<String>,
        lo: Option<&Bound<'_, PyAny>>,
        hi: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<usize> {
        if let Some(s) = sense
            && !s.is_none()
            && !s.is_instance_of::<PyString>()
        {
            return Err(PyTypeError::new_err(
                "sense must be a str: \"<=\", \">=\" or \"==\"",
            ));
        }
        let spec = BoundSpec::read(sense, rhs, lo, hi)?;
        if let Some(v) = [&spec].iter().find_map(|s| match s {
            BoundSpec::Sense(_, r) => r.len(),
            BoundSpec::Range(l, h) => l
                .as_ref()
                .and_then(F64s::len)
                .or(h.as_ref().and_then(F64s::len)),
            BoundSpec::None => None,
        }) {
            return Err(PyTypeError::new_err(format!(
                "add_constr adds one row: rhs/lo/hi must be numbers (got {v} values); use add_constrs for many"
            )));
        }
        let r = self.add_one(row, &spec, 0, name.as_deref())?;
        if let Some(n) = name {
            self.b
                .set_row_names(r..r + 1, Names::List(vec![n]))
                .map_err(|e| build_err(e, None))?;
        }
        self.record_source(r, r + 1);
        Ok(r)
    }

    /// Add many rows; returns their index range. `rows`: an iterable of rows (as for
    /// `add_constr`) or Constraints; `sense`: one str or one per row; `rhs`/`lo`/`hi`: scalar
    /// or one per row; `names`: a list of str or a prefix.
    #[pyo3(signature = (rows, sense=None, rhs=None, names=None, *, lo=None, hi=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_constrs<'py>(
        &mut self,
        py: Python<'py>,
        rows: &Bound<'py, PyAny>,
        sense: Option<&Bound<'py, PyAny>>,
        rhs: Option<&Bound<'py, PyAny>>,
        names: Option<&Bound<'py, PyAny>>,
        lo: Option<&Bound<'py, PyAny>>,
        hi: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyRange>> {
        let spec = BoundSpec::read(sense, rhs, lo, hi)?;
        let first = self.b.num_rows();
        let result = (|| -> PyResult<usize> {
            let mut k = 0;
            if let Ok(l) = rows.cast_exact::<PyList>() {
                for i in 0..l.len() {
                    self.add_one(&l.get_item(i)?, &spec, k, None)?;
                    k += 1;
                }
            } else if let Ok(t) = rows.cast_exact::<PyTuple>() {
                for row in t.iter() {
                    self.add_one(&row, &spec, k, None)?;
                    k += 1;
                }
            } else {
                for row in rows.try_iter()? {
                    self.add_one(&row?, &spec, k, None)?;
                    k += 1;
                }
            }
            spec.check_len(k)?;
            if let Some(n) = names_arg(names, k, "names")? {
                self.b
                    .set_row_names(first..first + k, n)
                    .map_err(|e| build_err(e, None))?;
            }
            Ok(k)
        })();
        match result {
            Ok(k) => {
                self.record_source(first, first + k);
                PyRange::new(py, first as isize, (first + k) as isize)
            }
            Err(e) => {
                self.b.truncate_rows(first);
                Err(e)
            }
        }
    }

    /// Add the rows of a matrix; returns their index range. `A`: scipy.sparse (csr, csc, coo,
    /// or anything with `tocsr()`), a 2-D float buffer (dense numpy), or a list of lists;
    /// `vars`: the variables of `A`'s columns (default: column j is variable j).
    #[pyo3(signature = (a, sense=None, rhs=None, vars=None, names=None, *, lo=None, hi=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_matrix<'py>(
        &mut self,
        py: Python<'py>,
        a: &Bound<'py, PyAny>,
        sense: Option<&Bound<'py, PyAny>>,
        rhs: Option<&Bound<'py, PyAny>>,
        vars: Option<&Bound<'py, PyAny>>,
        names: Option<&Bound<'py, PyAny>>,
        lo: Option<&Bound<'py, PyAny>>,
        hi: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyRange>> {
        let spec = BoundSpec::read(sense, rhs, lo, hi)?;
        let m = read_matrix(a)?;
        let (nrows, ncols) = m.shape;
        let col_map: Option<Vec<i32>> = match vars.filter(|v| !v.is_none()) {
            None => None,
            Some(v) => {
                let cols = Cols::read(v, self.id, self.b.num_cols(), "vars")?;
                if cols.len() != ncols {
                    return Err(crate::diag::model_err(format!(
                        "vars: the matrix has {ncols} columns but {} variables were given",
                        cols.len()
                    )));
                }
                match cols {
                    Cols::Range(r) if r.start == 0 => None,
                    c => Some((0..c.len()).map(|k| c.get(k) as i32).collect()),
                }
            }
        };
        let names = names_arg(names, nrows, "names")?;
        let (lf, uf) = spec.fills(nrows)?;
        let (lower, upper) = (lf.fill(&spec), uf.fill(&spec));
        let cm = col_map.as_deref();
        let data = m.data.slice();
        macro_rules! dispatch {
            ($f:ident, $($pre:expr,)* ; $p:expr, $i:expr) => {
                match ($p, $i) {
                    (IntSlice::I32(p), IntSlice::I32(i)) => self.b.$f($($pre,)* p, i, data, cm, lower, upper),
                    (IntSlice::I32(p), IntSlice::I64(i)) => self.b.$f($($pre,)* p, i, data, cm, lower, upper),
                    (IntSlice::I64(p), IntSlice::I32(i)) => self.b.$f($($pre,)* p, i, data, cm, lower, upper),
                    (IntSlice::I64(p), IntSlice::I64(i)) => self.b.$f($($pre,)* p, i, data, cm, lower, upper),
                }
            };
        }
        let r = match &m.kind {
            MatKind::Csr(p, i) => dispatch!(add_rows_csr, ; p.slice(), i.slice()),
            MatKind::Csc(p, i) => dispatch!(add_rows_csc, nrows, ; p.slice(), i.slice()),
            MatKind::Coo(r, c) => dispatch!(add_rows_coo, nrows, ; r.slice(), c.slice()),
            MatKind::Dense => self.b.add_rows_dense(nrows, ncols, data, cm, lower, upper),
        }
        .map_err(|e| build_err(e, None))?;
        if let Some(n) = names {
            self.b
                .set_row_names(r.clone(), n)
                .map_err(|e| build_err(e, None))?;
        }
        self.record_source(r.start, r.end);
        PyRange::new(py, r.start as isize, r.end as isize)
    }

    // ------------------------------------------------------------ objective

    /// Replace the objective. `obj`: a dict `{var: coef}`, LinExpr (its constant becomes the
    /// offset), Var, `(indices, values)`, or a dense sequence/buffer of length `num_vars`.
    #[pyo3(signature = (obj, sense=None))]
    fn set_objective(&mut self, obj: &Bound<'_, PyAny>, sense: Option<&str>) -> PyResult<()> {
        let (c, offset) = self.dense_objective(obj)?;
        if let Some(s) = sense {
            let s = obj_sense(s)?;
            self.set_sense_inner(s)?;
        }
        let s = self.synced_cols();
        if let Some(p) = self.port.as_mut() {
            if s > 0 {
                port_err(
                    "change_cols_cost",
                    p.change_cols_cost_interval(0, s as HighsInt - 1, &c[..s]),
                )?;
            }
            port_err("change_objective_offset", p.change_objective_offset(offset))?;
        }
        self.b.clear_unsynced_costs();
        for (j, &v) in c.iter().enumerate().skip(s) {
            if v != 0.0 {
                self.b.set_col_cost(j, v).map_err(|e| build_err(e, None))?;
            }
        }
        self.b.set_offset(offset);
        Ok(())
    }

    /// Set the quadratic part of the objective: `objective = c'x + (1/2) x'Qx`.
    ///
    /// `Q`: a square matrix (scipy.sparse, a 2-D array or a list of lists) or a dict
    /// `{(i, j): value}`; or the keywords `rows`, `cols`, `vals` (triplets). Duplicate entries are
    /// summed and an asymmetric Q is replaced by `(Q + Q')/2`, which gives the same objective.
    /// A Q smaller than the number of variables covers the first variables. `None` clears it.
    /// Convexity (PSD for `min`, NSD for `max`) is checked at `solve()`; see
    /// `rusolve_core::qp_convexity`.
    #[pyo3(signature = (q=None, *, rows=None, cols=None, vals=None))]
    fn set_quadratic_objective(
        &mut self,
        q: Option<&Bound<'_, PyAny>>,
        rows: Option<&Bound<'_, PyAny>>,
        cols: Option<&Bound<'_, PyAny>>,
        vals: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let entries = crate::quad::read_hessian(q, rows, cols, vals, self.id, self.b.num_cols())?;
        self.quad = entries;
        self.quad_dirty = true;
        Ok(())
    }

    /// What the convexity check found at the last `solve()` of a model with a quadratic
    /// objective: `{"class", "convexity", "method", "shift", "min_eigenvalue_estimate"}`, or None.
    fn quad_info<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let Some(i) = &self.quad_info else {
            return Ok(None);
        };
        let d = PyDict::new(py);
        d.set_item("class", &i.class)?;
        d.set_item("convexity", i.verdict)?;
        d.set_item("method", i.method)?;
        d.set_item("shift", i.shift)?;
        d.set_item("min_eigenvalue_estimate", i.min_eig)?;
        Ok(Some(d))
    }

    fn minimize(&mut self, obj: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_objective(obj, Some("min"))
    }

    fn maximize(&mut self, obj: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_objective(obj, Some("max"))
    }

    /// "min" or "max".
    #[getter]
    fn sense(&self) -> &'static str {
        match self.sense {
            Sense::Minimize => "min",
            Sense::Maximize => "max",
        }
    }

    // ------------------------------------------------------------ solve

    /// Solve (the GIL is released meanwhile) and return a Solution snapshot.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_var(ub=2.5); m.maximize(x)
    /// >>> m.solve(time_limit=10, threads=1).objective
    /// 2.5
    ///
    /// ``callback`` is any object whose methods are hooks (read the node LP, the
    /// incumbent and the cut pool, add cuts and lazy constraints and solutions, steer
    /// branching and node selection, interrupt); see ``rusolve.CallbackContext``. A hook
    /// that raises stops the solve and ``solve`` raises ``rusolve.CallbackError``.
    /// ``separate``, ``branching_score`` and ``select_node`` need the model's own column
    /// numbers, so they switch presolve off for that solve; ``lazy_rounds`` bounds the
    /// ``lazy_constraints`` re-solves.
    ///
    /// >>> import rusolve as rs
    /// >>> class Stop:
    /// ...     def on_progress(self, ctx: rs.CallbackContext) -> None:
    /// ...         ctx.interrupt()
    /// >>> m = rs.Model(); xs = m.add_vars(3, ub=1, binary=True)
    /// >>> _ = m.add_constr(xs[0] + xs[1] + xs[2] <= 2); m.maximize(xs[0] + xs[1] + xs[2])
    /// >>> m.solve(callback=Stop()).status in ("Optimal", "Interrupted")
    /// True
    #[pyo3(signature = (time_limit=None, mip_gap=None, threads=None, presolve=None, seed=None, verbose=false, options=None, callback=None, lazy_rounds=100))]
    #[allow(clippy::too_many_arguments)]
    fn solve(
        &mut self,
        py: Python<'_>,
        time_limit: Option<f64>,
        mip_gap: Option<f64>,
        threads: Option<i32>,
        presolve: Option<&Bound<'_, PyAny>>,
        seed: Option<i32>,
        verbose: bool,
        options: Option<&Bound<'_, PyDict>>,
        callback: Option<&Bound<'_, PyAny>>,
        lazy_rounds: usize,
    ) -> PyResult<Solution> {
        // presolve: True -> "on", False -> "off", a str passes through raw ("choose")
        let presolve: Option<String> = match presolve.filter(|p| !p.is_none()) {
            None => None,
            Some(p) if p.is_instance_of::<pyo3::types::PyBool>() => {
                Some(if p.is_truthy()? { "on" } else { "off" }.to_string())
            }
            Some(p) => Some(p.extract::<String>().map_err(|_| {
                PyTypeError::new_err(
                    "presolve must be a bool or a str (\"on\", \"off\", \"choose\")",
                )
            })?),
        };
        self.sync()?;
        self.sync_quad()?;
        self.apply_options(
            time_limit,
            mip_gap,
            threads,
            presolve.as_deref(),
            seed,
            verbose,
            options,
        )?;
        self.last_options = run_options(
            time_limit,
            mip_gap,
            threads,
            presolve.as_deref(),
            seed,
            verbose,
            options,
        );
        let callback = callback.filter(|c| !c.is_none());
        if callback.is_some() {
            self.callback_stats = None;
        }
        if let Some(cb) = callback
            && crate::callback::Bridge::needs_presolve_off(cb)
        {
            // hooks run in the presolved column space: use the model's own columns
            self.options_dirty = true;
            let port: &mut Port = self.port.as_mut().expect("synced");
            port.set_option_value_str("presolve", "off");
            // the run record must show what actually ran, so overwrite an earlier entry
            match self.last_options.iter_mut().find(|(k, _)| k == "presolve") {
                Some(e) => e.1 = "off".into(),
                None => self.last_options.push(("presolve".into(), "off".into())),
            }
        }
        let port: &mut Port = self.port.as_mut().expect("synced");
        // an LP solve after a MIP solve must not report the MIP's statistics
        port.mip_stats_ = None;
        let bridge = match callback {
            Some(cb) => Some(crate::callback::Bridge::install(py, cb, port)?),
            None => None,
        };
        let rows_before = port.model_lp_.num_row_;
        let t0 = Instant::now();
        // The port is plain owned Rust data (Send); `&mut self` keeps every other Python
        // thread out of this Model (a concurrent call gets "Already borrowed").
        let (run, lazy) = py.detach(|| {
            port.reset_scheduler_for_threads();
            match &bridge {
                Some(b) => b.run(port, lazy_rounds),
                None => (port.run(), None),
            }
        });
        let time = t0.elapsed().as_secs_f64();
        let mut sol = self.solution(run, time);
        let Some(bridge) = bridge else {
            return Ok(sol);
        };
        // Lazy rows live in the solver only for this solve; the Model stays as built.
        let port: &mut Port = self.port.as_mut().expect("synced");
        let added = port.model_lp_.num_row_ - rows_before;
        if added > 0 {
            if let Some(v) = sol.row_activity.as_mut() {
                v.truncate(rows_before as usize);
            }
            if let Some(v) = sol.duals.as_mut() {
                v.truncate(rows_before as usize);
            }
            port_err(
                "delete_rows_interval",
                port.delete_rows_interval(rows_before, rows_before + added - 1),
            )?;
        }
        if let Some(l) = &lazy
            && !l.accepted
            && (sol.status == "Optimal" || sol.status == "NotSolved")
        {
            // the candidate violates a row the callback would have added: not a solution
            // a hook that raised stopped the solve too: report it as interrupted, not as the
            // round limit (docs/design/hooks.md)
            sol.status = if bridge.interrupted() || bridge.failed() {
                "Interrupted"
            } else {
                "LazyRoundLimit"
            };
            sol.x = None;
            sol.objective = None;
            sol.row_activity = None;
            sol.gap = None;
        }
        let stats = bridge.finish(py, port, lazy.as_ref())?;
        self.callback_stats = Some(stats);
        if let Some(e) = bridge.error(py) {
            if let Ok(v) = Py::new(py, sol) {
                let _ = e.value(py).setattr("solution", v);
            }
            return Err(e);
        }
        Ok(sol)
    }

    /// Counters of the last ``solve(callback=...)``: ``calls`` per hook, ``errors``,
    /// ``interrupted``, the cuts the solver was offered and accepted, node-selection and
    /// branching overrides, user solutions accepted/rejected, and for lazy constraints the
    /// ``lazy`` rounds and rows added. None before such a solve.
    ///
    /// >>> import rusolve as rs
    /// >>> class Watch:
    /// ...     def on_improving_solution(self, ctx: rs.CallbackContext) -> None:
    /// ...         pass
    /// >>> m = rs.Model(); x = m.add_var(ub=3, integer=True); m.maximize(x)
    /// >>> _ = m.solve(callback=Watch())
    /// >>> stats = m.callback_stats
    /// >>> stats is not None and stats["errors"]
    /// 0
    #[getter]
    fn callback_stats(&self, py: Python<'_>) -> Option<Py<PyDict>> {
        self.callback_stats.as_ref().map(|d| d.clone_ref(py))
    }

    /// Record the Python line that adds each constraint (off by default). With it on,
    /// `source_of(row)` names the line and `explain()`/`iis()` put it on every row of the
    /// subsystem as `member["source"]`. Cost: one frame lookup per `add_constr`/`add_constrs`/
    /// `add_matrix` call (see crates/rusolve-py/README.md, "Source lines").
    #[getter]
    fn capture_source(&self) -> bool {
        self.capture_source
    }
    #[setter]
    fn set_capture_source(&mut self, on: bool) {
        self.capture_source = on;
    }

    /// The Python line that added constraint `row` (an int index): `{"path", "line",
    /// "source_line"}`, or None when `capture_source` was off then, or the row is unknown.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); m.capture_source = True; x = m.add_var()
    /// >>> r = m.add_constr(x <= 1)
    /// >>> src = m.source_of(r)
    /// >>> src is not None and src["source_line"]
    /// 'r = m.add_constr(x <= 1)'
    fn source_of(&self, py: Python<'_>, row: usize) -> PyResult<Py<PyAny>> {
        Ok(match self.source_dict(py, row)? {
            Some(d) => d.into_any().unbind(),
            None => py.None(),
        })
    }

    /// Hand everything pending to the solver without solving (internal: lets benchmarks
    /// time the build apart from the solve).
    fn _sync(&mut self) -> PyResult<()> {
        self.sync()
    }

    // ------------------------------------------------------------ changes

    /// Set the bounds of variables (a Var, int, Vars or sequence); `lb`/`ub` scalar or one per
    /// variable.
    fn set_bounds(
        &mut self,
        vars: &Bound<'_, PyAny>,
        lb: &Bound<'_, PyAny>,
        ub: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let cols = Cols::read(vars, self.id, self.b.num_cols(), "vars")?;
        let n = cols.len();
        let (lb, ub) = (F64s::read(lb, "lb")?, F64s::read(ub, "ub")?);
        let (l, u) = (lb.per_item(n, "lb")?, ub.per_item(n, "ub")?);
        // validate everything first: a failed call changes nothing
        for k in 0..n {
            if l.get(k).is_nan() || u.get(k).is_nan() {
                return Err(crate::diag::model_err(format!(
                    "bound of variable {} is NaN",
                    cols.get(k)
                )));
            }
        }
        let s = self.synced_cols();
        let mut synced = Vec::new();
        for k in 0..n {
            let (j, lo, hi) = (cols.get(k), l.get(k), u.get(k));
            if j < s {
                synced.push((j as HighsInt, (lo, hi)));
            } else {
                self.b
                    .set_col_bounds(j, lo, hi)
                    .map_err(|e| build_err(e, None))?;
            }
        }
        if !synced.is_empty() {
            let (set, v) = sorted_unique(synced);
            let (lo, hi): (Vec<f64>, Vec<f64>) = v.into_iter().unzip();
            let p = self.port.as_mut().expect("synced");
            port_err(
                "change_cols_bounds",
                p.change_cols_bounds_set(set.len() as HighsInt, &set, &lo, &hi),
            )?;
        }
        Ok(())
    }

    /// Set the objective coefficients of variables (the others keep theirs).
    fn set_cost(&mut self, vars: &Bound<'_, PyAny>, cost: &Bound<'_, PyAny>) -> PyResult<()> {
        let cols = Cols::read(vars, self.id, self.b.num_cols(), "vars")?;
        let n = cols.len();
        let cost = F64s::read(cost, "cost")?;
        let c = cost.per_item(n, "cost")?;
        for k in 0..n {
            let v = c.get(k);
            if !v.is_finite() {
                return Err(crate::diag::model_err(format!(
                    "cost of variable {} is {v}; costs must be finite",
                    cols.get(k)
                )));
            }
        }
        let s = self.synced_cols();
        let mut synced = Vec::new();
        for k in 0..n {
            let (j, v) = (cols.get(k), c.get(k));
            if j < s {
                synced.push((j as HighsInt, v));
            } else {
                self.b.set_col_cost(j, v).map_err(|e| build_err(e, None))?;
            }
        }
        if !synced.is_empty() {
            let (set, v) = sorted_unique(synced);
            let p = self.port.as_mut().expect("synced");
            port_err(
                "change_cols_cost",
                p.change_cols_cost_set(set.len() as HighsInt, &set, &v),
            )?;
        }
        Ok(())
    }

    /// Make variables semi-continuous (`integer=False`: the value is 0 or within its bounds)
    /// or semi-integer (`integer=True`: 0 or an integer within its bounds).
    #[pyo3(signature = (vars, integer=false))]
    fn set_semi(&mut self, vars: &Bound<'_, PyAny>, integer: bool) -> PyResult<()> {
        let cols = Cols::read(vars, self.id, self.b.num_cols(), "vars")?;
        let s = self.synced_cols();
        let t = if integer {
            HighsVarType::kSemiInteger
        } else {
            HighsVarType::kSemiContinuous
        };
        let mut synced = Vec::new();
        for k in 0..cols.len() {
            let j = cols.get(k);
            if j < s {
                synced.push((j as HighsInt, t));
            } else {
                self.b
                    .set_col_semi(j, integer)
                    .map_err(|e| build_err(e, None))?;
            }
        }
        if !synced.is_empty() {
            let (set, v) = sorted_unique(synced);
            let p = self.port.as_mut().expect("synced");
            port_err(
                "change_cols_integrality",
                p.change_cols_integrality_set(set.len() as HighsInt, &set, &v),
            )?;
        }
        Ok(())
    }

    /// Make variables integer (`flag=True`) or continuous.
    #[pyo3(signature = (vars, flag=None))]
    fn set_integer(
        &mut self,
        vars: &Bound<'_, PyAny>,
        flag: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let cols = Cols::read(vars, self.id, self.b.num_cols(), "vars")?;
        let n = cols.len();
        let flags = match flag.filter(|f| !f.is_none()) {
            None => Flags::Scalar(true),
            Some(f) => Flags::read(f, "flag")?,
        };
        if let Flags::Each(v) = &flags
            && v.len() != n
        {
            return Err(crate::diag::model_err(format!(
                "flag: expected {n} values, got {}",
                v.len()
            )));
        }
        let s = self.synced_cols();
        let mut synced = Vec::new();
        for k in 0..n {
            let j = cols.get(k);
            let f = match &flags {
                Flags::Scalar(f) => *f,
                Flags::Each(v) => v[k],
            };
            if j < s {
                let t = if f {
                    HighsVarType::kInteger
                } else {
                    HighsVarType::kContinuous
                };
                synced.push((j as HighsInt, t));
            } else {
                self.b
                    .set_col_integer(j, f)
                    .map_err(|e| build_err(e, None))?;
            }
        }
        if !synced.is_empty() {
            let (set, v) = sorted_unique(synced);
            let p = self.port.as_mut().expect("synced");
            port_err(
                "change_cols_integrality",
                p.change_cols_integrality_set(set.len() as HighsInt, &set, &v),
            )?;
        }
        Ok(())
    }

    /// Set the bounds `lo <= row <= hi` of rows (an int, range or sequence).
    fn set_row_bounds(
        &mut self,
        rows: &Bound<'_, PyAny>,
        lo: &Bound<'_, PyAny>,
        hi: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let rows = row_list(rows, self.b.num_rows())?;
        let n = rows.len();
        let (lo, hi) = (F64s::read(lo, "lo")?, F64s::read(hi, "hi")?);
        let (l, u) = (lo.per_item(n, "lo")?, hi.per_item(n, "hi")?);
        for (k, &i) in rows.iter().enumerate() {
            if l.get(k).is_nan() || u.get(k).is_nan() {
                return Err(crate::diag::model_err(format!("bound of row {i} is NaN")));
            }
        }
        let s = self.synced_rows();
        let mut synced = Vec::new();
        for (k, &i) in rows.iter().enumerate() {
            let (a, b) = (l.get(k), u.get(k));
            if i < s {
                synced.push((i as HighsInt, (a, b)));
            } else {
                self.b
                    .set_row_bounds(i, a, b)
                    .map_err(|e| build_err(e, None))?;
            }
        }
        if !synced.is_empty() {
            let (set, v) = sorted_unique(synced);
            let (lo, hi): (Vec<f64>, Vec<f64>) = v.into_iter().unzip();
            let p = self.port.as_mut().expect("synced");
            port_err(
                "change_rows_bounds",
                p.change_rows_bounds_set(set.len() as HighsInt, &set, &lo, &hi),
            )?;
        }
        Ok(())
    }

    /// Set the coefficient of `var` in row `row` (0 removes it).
    fn set_coef(&mut self, row: usize, var: &Bound<'_, PyAny>, value: f64) -> PyResult<()> {
        let j = self.col_index(var)?;
        if row >= self.b.num_rows() {
            return Err(crate::diag::index_err(format!(
                "row {row} is out of range (there are {} rows)",
                self.b.num_rows()
            )));
        }
        if !value.is_finite() {
            return Err(crate::diag::model_err(format!(
                "coefficient of variable {j} in row {row} is {value}; coefficients must be finite"
            )));
        }
        if row < self.synced_rows() {
            if j < self.synced_cols() {
                let p = self.port.as_mut().expect("synced");
                port_err(
                    "change_coeff",
                    p.change_coeff(row as HighsInt, j as HighsInt, value),
                )?;
            } else {
                self.pending_coefs
                    .push((row as HighsInt, j as HighsInt, value));
            }
            Ok(())
        } else {
            self.b
                .set_coef(row, j, value)
                .map_err(|e| build_err(e, None))
        }
    }

    // ------------------------------------------------------------ queries

    #[getter]
    fn num_vars(&self) -> usize {
        self.b.num_cols()
    }
    #[getter]
    fn num_constrs(&self) -> usize {
        self.b.num_rows()
    }
    #[getter]
    fn num_nonzeros(&self) -> usize {
        let synced = self
            .port
            .as_ref()
            .map_or(0, |p| p.get_num_nz().max(0) as usize);
        synced + self.b.num_pending_nonzeros()
    }

    /// The name of a variable, or None.
    fn var_name(&self, v: &Bound<'_, PyAny>) -> PyResult<Option<String>> {
        let j = self.col_index(v)?;
        if j < self.synced_cols() {
            let lp = self.port.as_ref().expect("synced").get_lp();
            return Ok(lp.col_names_.get(j).filter(|s| !s.is_empty()).cloned());
        }
        Ok(self.b.col_name(j).map(str::to_string))
    }

    /// The name of a row, or None.
    fn constr_name(&self, r: usize) -> PyResult<Option<String>> {
        if r >= self.b.num_rows() {
            return Err(crate::diag::index_err(format!(
                "row {r} is out of range (there are {} rows)",
                self.b.num_rows()
            )));
        }
        if r < self.synced_rows() {
            let lp = self.port.as_ref().expect("synced").get_lp();
            return Ok(lp.row_names_.get(r).filter(|s| !s.is_empty()).cloned());
        }
        Ok(self.b.row_name(r).map(str::to_string))
    }

    // ------------------------------------------------------------ diagnostics
    //
    // Each method is the Python call of a CLI subcommand and returns the same `rusolve.<cmd>/1`
    // object, parsed (`json.loads` of what `--format json` prints). They read the model as the
    // solver holds it, so they never change it; `ranging` and `feasrelax` solve a copy.

    /// Check an assignment against this model, like `rusolve check`: row and column bounds,
    /// integrality, objective value, and the worst violations first. `solution` is a
    /// `Solution` of this model, a sequence or array of values (one per variable), or a dict
    /// `{Var or name: value}` covering every variable. Returns the `rusolve.check/1` dict;
    /// `verdict` is `"feasible"` or says what failed. A quadratic objective is not supported.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_vars(2, ub=4); _ = m.add_constr(x[0] + x[1] <= 5)
    /// >>> m.maximize(x[0] + x[1])
    /// >>> m.check([4.0, 4.0])["verdict"]
    /// 'infeasible'
    /// >>> m.check(m.solve())["verdict"]
    /// 'feasible'
    #[pyo3(signature = (solution, tol=None, mip_tol=None, max_findings=None))]
    fn check(
        &mut self,
        py: Python<'_>,
        solution: &Bound<'_, PyAny>,
        tol: Option<f64>,
        mip_tol: Option<f64>,
        max_findings: Option<usize>,
    ) -> PyResult<Py<PyAny>> {
        self.refuse_quad("check")?;
        self.sync()?;
        let values = self.assignment(solution)?;
        let name = self.display_name();
        let lp = &self.port.as_ref().expect("synced").model_lp_;
        let v = rusolve_cli::api::check(lp, &values, tol, mip_tol, max_findings, &name)
            .map_err(crate::diag::model_err)?;
        json_to_py(py, &v)
    }

    /// Why is the model infeasible (or unbounded): an irreducible infeasible subsystem of
    /// rows and bounds, or the unbounded ray, like `rusolve explain`. Returns the
    /// `rusolve.explain/1` dict; `verify=True` re-solves the subsystem to prove it infeasible.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_var(); _ = m.add_constr(x >= 3, name="lo"); _ = m.add_constr(x <= 1, name="hi")
    /// >>> e = m.explain()
    /// >>> e["schema_version"], e["status"]
    /// ('rusolve.explain/1', 'Infeasible')
    #[pyo3(signature = (verify=true, max_members=1000))]
    fn explain(&mut self, py: Python<'_>, verify: bool, max_members: usize) -> PyResult<Py<PyAny>> {
        // The port's LP has no Hessian, so every verdict would ignore the quadratic term:
        // refuse before spending a solve on it.
        self.refuse_quad("explain")?;
        self.sync()?;
        let name = self.display_name();
        let lp = &self.port.as_ref().expect("synced").model_lp_;
        let v = py
            .detach(|| rusolve_cli::api::explain(lp, &name, verify, max_members))
            .map_err(crate::diag::model_err)?;
        let out = json_to_py(py, &v)?;
        if self.capture_source {
            self.annotate_sources(py, out.bind(py))?;
        }
        Ok(out)
    }

    /// Alias of [`explain`](Model.explain): the IIS of an infeasible model.
    #[pyo3(signature = (verify=true, max_members=1000))]
    fn iis(&mut self, py: Python<'_>, verify: bool, max_members: usize) -> PyResult<Py<PyAny>> {
        self.explain(py, verify, max_members)
    }

    /// Static model lint, like `rusolve lint`: big-M coefficients, wide coefficient ranges,
    /// free or unbounded columns, ... Returns the `rusolve.lint/1` dict. `rules` limits the
    /// check to rule ids, `min_severity` is `"info"`, `"warning"` or `"error"`.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_var(); y = m.add_var(binary=True)
    /// >>> _ = m.add_constr(x <= 1e9 * y); m.minimize(x)
    /// >>> m.lint()["schema_version"]
    /// 'rusolve.lint/1'
    #[pyo3(signature = (rules=None, min_severity=None, max_items=None, with_solve=false))]
    fn lint(
        &mut self,
        py: Python<'_>,
        rules: Option<Vec<String>>,
        min_severity: Option<String>,
        max_items: Option<usize>,
        with_solve: bool,
    ) -> PyResult<Py<PyAny>> {
        if with_solve {
            self.refuse_quad("lint(with_solve=True)")?;
        }
        self.sync()?;
        let name = self.display_name();
        let rules = rules.unwrap_or_default();
        let lp = &self.port.as_ref().expect("synced").model_lp_;
        let v = py
            .detach(|| {
                rusolve_cli::api::lint(
                    lp,
                    &name,
                    &rules,
                    min_severity.as_deref(),
                    max_items,
                    with_solve,
                )
            })
            .map_err(crate::diag::model_err)?;
        json_to_py(py, &v)
    }

    /// Formulation error or data error, like `rusolve diagnose --format json`: for an
    /// infeasible model, a heuristic verdict (`formulation`, `data`, `unclear`, or `none` when
    /// the model is not infeasible) with its evidence and a suggested fix, from the IIS, the
    /// Farkas ray, the minimum-count relaxation and the name families. Returns the
    /// `rusolve.diagnose/1` dict. A MIP is diagnosed through its LP relaxation; a quadratic
    /// objective is not supported. `time_limit` (seconds) bounds the relaxation.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_var(name="x"); _ = m.add_constr(x >= 3); _ = m.add_constr(x <= 1)
    /// >>> d = m.diagnose()
    /// >>> d["schema_version"], d["status"]
    /// ('rusolve.diagnose/1', 'Infeasible')
    #[pyo3(signature = (time_limit=30.0))]
    fn diagnose(&mut self, py: Python<'_>, time_limit: f64) -> PyResult<Py<PyAny>> {
        self.refuse_quad("diagnose")?;
        self.sync()?;
        let name = self.display_name();
        let lp = &self.port.as_ref().expect("synced").model_lp_;
        let v = py
            .detach(|| rusolve_cli::api::diagnose(lp, &name, time_limit))
            .map_err(crate::diag::model_err)?;
        json_to_py(py, &v)
    }

    /// Sensitivity ranging of the optimal solution, like `rusolve ranging`: how far each
    /// cost, bound and right-hand side can move before the basis changes. Solves a copy of the
    /// model. `rows` and `cols` are name patterns (default all). Returns the
    /// `rusolve.ranging/1` dict. A quadratic objective is not supported.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_vars(2, ub=4); _ = m.add_constr(x[0] + x[1] <= 5)
    /// >>> m.maximize(x[0] + 2 * x[1])
    /// >>> m.ranging()["schema_version"]
    /// 'rusolve.ranging/1'
    #[pyo3(signature = (rows=None, cols=None))]
    fn ranging(
        &mut self,
        py: Python<'_>,
        rows: Option<Vec<String>>,
        cols: Option<Vec<String>>,
    ) -> PyResult<Py<PyAny>> {
        self.refuse_quad("ranging")?;
        self.sync()?;
        let name = self.display_name();
        let (rows, cols) = (rows.unwrap_or_default(), cols.unwrap_or_default());
        let lp = &self.port.as_ref().expect("synced").model_lp_;
        let v = py
            .detach(|| rusolve_cli::api::ranging(lp, &name, &rows, &cols))
            .map_err(crate::diag::model_err)?;
        json_to_py(py, &v)
    }

    /// Feasibility relaxation of an infeasible model, like `rusolve feasrelax`: the cheapest
    /// weighted (`mode="weighted"`, the default) or fewest (`"min-count"`) bound and row
    /// relaxations that make it feasible. `penalties` is the `--penalties` document as a dict.
    /// Returns the `rusolve.feasrelax/1` dict.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_var(); _ = m.add_constr(x >= 3); _ = m.add_constr(x <= 1)
    /// >>> r = m.feasrelax()
    /// >>> r["schema_version"]
    /// 'rusolve.feasrelax/1'
    #[pyo3(signature = (mode=None, global_lower_penalty=None, global_upper_penalty=None, global_rhs_penalty=None, big_m=None, time_limit=None, penalties=None))]
    #[allow(clippy::too_many_arguments)]
    fn feasrelax(
        &mut self,
        py: Python<'_>,
        mode: Option<String>,
        global_lower_penalty: Option<f64>,
        global_upper_penalty: Option<f64>,
        global_rhs_penalty: Option<f64>,
        big_m: Option<f64>,
        time_limit: Option<f64>,
        penalties: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        // verify() and objective_original use the port's HighsLp, which has no Hessian.
        self.refuse_quad("feasrelax")?;
        if let Some(m) = &mode
            && m != "weighted"
            && m != "min-count"
        {
            return Err(crate::diag::model_err(format!(
                "mode: expected \"weighted\" or \"min-count\", got \"{m}\""
            )));
        }
        let mut flags: Vec<String> = Vec::new();
        if let Some(m) = mode {
            flags.extend(["--mode".to_string(), m]);
        }
        for (v, flag) in [
            (global_lower_penalty, "--global-lower-penalty"),
            (global_upper_penalty, "--global-upper-penalty"),
            (global_rhs_penalty, "--global-rhs-penalty"),
            (big_m, "--big-m"),
            (time_limit, "--time-limit"),
        ] {
            if let Some(x) = v {
                flags.extend([flag.to_string(), x.to_string()]);
            }
        }
        let penalties: Option<String> = match penalties {
            Some(p) if !p.is_none() => Some(
                py.import("json")?
                    .call_method1("dumps", (p,))?
                    .extract::<String>()?,
            ),
            _ => None,
        };
        self.sync()?;
        let name = self.display_name();
        let lp = &self.port.as_ref().expect("synced").model_lp_;
        let v = py
            .detach(|| rusolve_cli::api::feasrelax(lp, &name, &flags, penalties.as_deref()))
            .map_err(crate::diag::model_err)?;
        json_to_py(py, &v)
    }

    /// Write the model to a file in the format the extension names: `.mps`, `.lp`, either
    /// with `.gz`. A quadratic objective needs `.lp` (the MPS writer has no quadratic section);
    /// semi-continuous variables are not supported. The file reads back with `rusolve solve`, `rusolve check` and the other
    /// commands.
    ///
    /// >>> import rusolve as rs, tempfile, os
    /// >>> m = rs.Model(); x = m.add_vars(2, ub=1, names="x"); m.maximize(x[0] + x[1])
    /// >>> path = os.path.join(tempfile.mkdtemp(), "m.lp")
    /// >>> m.write(path); os.path.getsize(path) > 0
    /// True
    fn write(&mut self, path: &str) -> PyResult<()> {
        self.sync()?;
        let lp = &self.port.as_ref().expect("synced").model_lp_;
        let hessian = self.quad.as_deref().unwrap_or(&[]);
        rusolve_cli::api::write_model(lp, hessian, path).map_err(crate::diag::model_err)
    }

    /// The Var of column `j`.
    fn var(&self, j: usize) -> PyResult<Var> {
        if j >= self.b.num_cols() {
            return Err(crate::diag::index_err(format!(
                "variable {j} is out of range (there are {} variables)",
                self.b.num_cols()
            )));
        }
        Ok(Var::new(j as i32, self.id))
    }

    fn __repr__(&self) -> String {
        format!(
            "Model({} vars, {} constrs, {} nonzeros, {})",
            self.b.num_cols(),
            self.b.num_rows(),
            self.num_nonzeros(),
            self.sense()
        )
    }
}

// ---------------------------------------------------------------- matrices

pub(crate) enum MatKind {
    Csr(Ints, Ints),
    Csc(Ints, Ints),
    Coo(Ints, Ints),
    Dense,
}

pub(crate) struct Matrix {
    pub(crate) kind: MatKind,
    pub(crate) data: F64s,
    pub(crate) shape: (usize, usize),
}

/// Read a scipy.sparse matrix by its attributes (no scipy import), a 2-D buffer, or a list
/// of lists.
pub(crate) fn read_matrix(a: &Bound<'_, PyAny>) -> PyResult<Matrix> {
    if let Ok(fmt) = a.getattr("format")
        && let Ok(fmt) = fmt.extract::<String>()
        && let Ok(shape) = a.getattr("shape")
    {
        let (nr, nc): (usize, usize) = shape
            .extract()
            .map_err(|_| crate::diag::model_err("the matrix must be 2-D"))?;
        let ints = |name: &str| -> PyResult<Ints> { Ints::read(&a.getattr(name)?, name, None) };
        let data = || F64s::read_array(&a.getattr("data")?, "data");
        let kind = match fmt.as_str() {
            "csr" => MatKind::Csr(ints("indptr")?, ints("indices")?),
            "csc" => MatKind::Csc(ints("indptr")?, ints("indices")?),
            "coo" => MatKind::Coo(ints("row")?, ints("col")?),
            _ => return read_matrix(&a.call_method0("tocsr")?),
        };
        return Ok(Matrix {
            kind,
            data: data()?,
            shape: (nr, nc),
        });
    }
    if let Ok(b) = pyo3::buffer::PyUntypedBuffer::get(a)
        && b.dimensions() == 2
    {
        let (nr, nc) = (b.shape()[0], b.shape()[1]);
        let data = if let Ok(t) = b.as_typed::<f64>()
            && b.is_c_contiguous()
        {
            let _ = t;
            F64s::Buf(b.into_typed::<f64>()?)
        } else {
            // other dtypes or layouts: convert row by row through Python
            let mut v = Vec::with_capacity(nr * nc);
            for row in a.try_iter()? {
                let r = F64s::read_array(&row?, "matrix row")?;
                v.extend_from_slice(r.slice());
            }
            F64s::Owned(v)
        };
        return Ok(Matrix {
            kind: MatKind::Dense,
            data,
            shape: (nr, nc),
        });
    }
    // a list of lists
    let mut v = Vec::new();
    let mut nr = 0;
    let mut nc = None;
    for row in a.try_iter().map_err(|_| {
        PyTypeError::new_err(format!(
            "add_matrix: expected a scipy.sparse matrix, a 2-D array or a list of lists, got {}",
            ingest::type_name(a)
        ))
    })? {
        let r = F64s::read_array(&row?, "matrix row")?;
        let s = r.slice();
        if *nc.get_or_insert(s.len()) != s.len() {
            return Err(crate::diag::model_err(
                "add_matrix: rows of different lengths",
            ));
        }
        v.extend_from_slice(s);
        nr += 1;
    }
    Ok(Matrix {
        kind: MatKind::Dense,
        data: F64s::Owned(v),
        shape: (nr, nc.unwrap_or(0)),
    })
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Model>()?;
    m.add_class::<crate::callback::CallbackContext>()?;
    Ok(())
}
