//! Python callbacks for `Model.solve(callback=...)` (docs/design/hooks.md section "Python").
//!
//! A callback is any Python object; each method it defines is one hook, looked up once per solve
//! (a missing method costs nothing at run time):
//!
//! | method | fires | returns |
//! |---|---|---|
//! | `on_progress(ctx)` | every MIP limit check (HiGHS `kCallbackMipInterrupt`) | nothing; `ctx.interrupt()` stops |
//! | `on_improving_solution(ctx)` | a new incumbent (`kCallbackMipImprovingSolution`) | nothing |
//! | `on_cut_pool(ctx)` | the cut pool is read (`kCallbackMipGetCutPool`) | nothing |
//! | `user_solution(ctx)` | the solver asks for a solution (`kCallbackMipUserSolution`) | a list of values or None |
//! | `separate(ctx)` | after the built-in separators (`MipHooks::separate`) | cuts `a x <= rhs` |
//! | `branching_score(ctx)` | candidate selection (`MipHooks::branching_score`) | one score per candidate |
//! | `select_node(ctx)` | node selection (`MipHooks::select_node`) | an index into `ctx.open_nodes` |
//! | `lazy_constraints(ctx)` | an optimal candidate (`Highs::run_with_lazy_constraints`) | violated rows |
//!
//! Failure policy: an exception inside a hook never unwinds into the solver. The first one is
//! kept, counted, and from then on no hook is called (each falls back to the solver's own rule)
//! and the solve is interrupted at the next limit check; after the solve `Model.solve` raises
//! `rusolve.CallbackError` (code RS4008) with the exception as `__cause__` and the partial
//! `Solution` as `.solution`. Hooks of the `separate`, `branching_score` and `select_node` kind
//! run in the presolved column space, so installing one switches presolve off for that solve
//! (docs/design/hooks.md); `lazy_constraints`, `user_solution` and the `on_*` events use the
//! model's own columns.
//!
//! Provenance: callback and lazy-constraint interfaces as in common MIP solver APIs (SCIP's
//! constraint-handler documentation for lazy constraints; Dantzig, Fulkerson, Johnson, Operations Research 2 (1954) 393-410, for the
//! subtour example); the HiGHS callback system is ported in `lp_data/highs_callback.rs`.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rusolve_core::diagnostic::Diagnostic;
use rusolve_ekk::HighsInt;
use rusolve_ekk::lp_data::hconst::*;
use rusolve_ekk::lp_data::highs_callback::{HighsCallbackFunctionType, HighsCallbackOutput};
use rusolve_ekk::lp_data::highs_run::Highs as Port;
use rusolve_ekk::mip::hooks::*;
use rusolve_ekk::mip::lazy::{LazyConstraints, LazyReport};

use crate::expr::Constraint;

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(|e| e.into_inner())
}

/// `(indices, values, lower, upper)`.
type PoolRow = (Vec<i32>, Vec<f64>, f64, f64);
/// `(lower_bound, estimate, depth, num_domain_changes)`.
type OpenNodeRow = (f64, f64, i32, usize);

/// What a hook of the object given to `Model.solve(callback=...)` is told. Fields that do not
/// apply to the running hook (`event`) are None.
///
/// A cover-cut separator (the hook runs in the model's own columns; the cut says that
/// picking all three items is not allowed):
///
///     >>> import rusolve as rs
///     >>> class Cover:
///     ...     def separate(self, ctx: rs.CallbackContext) -> list[tuple[list[int], list[float], float]]:
///     ...         return [([0, 1, 2], [1.0, 1.0, 1.0], 2.0)]
///     >>> m = rs.Model(); x = m.add_vars(3, ub=1, binary=True)
///     >>> _ = m.add_constr(2 * x[0] + 2 * x[1] + 2 * x[2] <= 5)
///     >>> m.maximize(x[0] + x[1] + x[2])
///     >>> m.solve(callback=Cover(), threads=1).objective
///     2.0
#[pyclass(module = "rusolve", name = "CallbackContext")]
#[derive(Default)]
pub struct CallbackContext {
    /// The hook that is running (`"separate"`, `"on_progress"`, ...).
    #[pyo3(get)]
    event: String,
    /// Column values: the incumbent, the node LP solution or the lazy candidate.
    #[pyo3(get)]
    solution: Option<Vec<f64>>,
    /// Objective of that solution (or the node LP objective in `separate`).
    #[pyo3(get)]
    objective: Option<f64>,
    /// Best dual (lower, when minimising) bound.
    #[pyo3(get)]
    bound: Option<f64>,
    #[pyo3(get)]
    gap: Option<f64>,
    #[pyo3(get)]
    nodes: Option<i64>,
    #[pyo3(get)]
    lp_iterations: Option<i64>,
    #[pyo3(get)]
    running_time: Option<f64>,
    /// Node (local) column lower and upper bounds, in the hook's column space.
    #[pyo3(get)]
    lower: Option<Vec<f64>>,
    #[pyo3(get)]
    upper: Option<Vec<f64>>,
    #[pyo3(get)]
    depth: Option<i64>,
    /// `branching_score`: the fractional integer columns as `(column, value)`.
    #[pyo3(get)]
    candidates: Option<Vec<(i32, f64)>>,
    /// `select_node`: open nodes as `(lower_bound, estimate, depth, num_domain_changes)`.
    #[pyo3(get)]
    open_nodes: Option<Vec<OpenNodeRow>>,
    /// `on_cut_pool`: rows as `(indices, values, lower, upper)`.
    #[pyo3(get)]
    cut_pool: Option<Vec<PoolRow>>,
    interrupt: bool,
    stop: Option<Arc<AtomicBool>>,
}

#[pymethods]
impl CallbackContext {
    /// Ask the solver to stop at its next limit check; the solve ends with status
    /// `"Interrupted"` and the best solution so far.
    fn interrupt(&mut self) {
        self.interrupt = true;
        if let Some(s) = &self.stop {
            s.store(true, Ordering::Release);
        }
    }
}

struct Has {
    on_progress: bool,
    on_improving_solution: bool,
    on_cut_pool: bool,
    user_solution: bool,
    separate: bool,
    branching_score: bool,
    select_node: bool,
    lazy_constraints: bool,
}

impl Has {
    fn of(cb: &Bound<'_, PyAny>) -> Self {
        let h = |n: &str| cb.hasattr(n).unwrap_or(false);
        Has {
            on_progress: h("on_progress"),
            on_improving_solution: h("on_improving_solution"),
            on_cut_pool: h("on_cut_pool"),
            user_solution: h("user_solution"),
            separate: h("separate"),
            branching_score: h("branching_score"),
            select_node: h("select_node"),
            lazy_constraints: h("lazy_constraints"),
        }
    }
    /// A hook that runs in the presolved column space.
    fn presolved_space(&self) -> bool {
        self.separate || self.branching_score || self.select_node
    }
}

struct Failure {
    hook: &'static str,
    err: PyErr,
}

/// State shared by the solver closures and `Model.solve`.
pub(crate) struct Bridge {
    cb: Py<PyAny>,
    has: Has,
    failure: Mutex<Option<Failure>>,
    errors: AtomicU64,
    calls: Mutex<BTreeMap<&'static str, u64>>,
    stop: Arc<AtomicBool>,
    num_col: usize,
}

impl Bridge {
    pub(crate) fn failed(&self) -> bool {
        lock(&self.failure).is_some()
    }

    fn record(&self, hook: &'static str, err: PyErr) {
        self.errors.fetch_add(1, Ordering::Relaxed);
        let mut f = lock(&self.failure);
        if f.is_none() {
            *f = Some(Failure { hook, err });
        }
        self.stop.store(true, Ordering::Release);
    }

    fn ctx(&self, event: &str) -> CallbackContext {
        CallbackContext {
            event: event.to_string(),
            stop: Some(Arc::clone(&self.stop)),
            ..Default::default()
        }
    }

    /// Call hook `hook`; `conv` turns its return value into the Rust answer. `None` when a hook
    /// failed before or in this call (the caller then falls back to the solver's rule).
    fn call<T>(
        &self,
        hook: &'static str,
        ctx: CallbackContext,
        conv: impl FnOnce(&Bound<'_, PyAny>) -> PyResult<T>,
    ) -> Option<T> {
        if self.failed() {
            return None;
        }
        Python::attach(|py| {
            let run = || -> PyResult<T> {
                let ctx = Py::new(py, ctx)?;
                let r = self.cb.bind(py).call_method1(hook, (ctx.clone_ref(py),))?;
                if ctx.borrow(py).interrupt {
                    self.stop.store(true, Ordering::Release);
                }
                conv(&r)
            };
            *lock(&self.calls).entry(hook).or_insert(0) += 1;
            match run() {
                Ok(v) => Some(v),
                Err(e) => {
                    self.record(hook, e);
                    None
                }
            }
        })
    }

    fn user_cuts(&self, hook: &'static str, ctx: CallbackContext) -> Vec<UserCut> {
        self.call(hook, ctx, parse_cuts).unwrap_or_default()
    }
}

/// Cuts from the Python answer: None, one cut, or an iterable of cuts. A cut is a `Constraint`
/// (`x[0] + x[1] <= 1`, `lo <= expr <= hi` gives up to two rows), `(indices, values, rhs)` or
/// `(coefficients, rhs)` with coefficients a dict `{column: value}` or `[(column, value), ...]`.
fn parse_cuts(r: &Bound<'_, PyAny>) -> PyResult<Vec<UserCut>> {
    if r.is_none() {
        return Ok(Vec::new());
    }
    let mut out = Vec::new();
    if r.cast::<Constraint>().is_ok() {
        parse_cut(r, &mut out)?;
        return Ok(out);
    }
    for item in r.try_iter()? {
        parse_cut(&item?, &mut out)?;
    }
    Ok(out)
}

fn merged(terms: impl IntoIterator<Item = (i32, f64)>) -> (Vec<HighsInt>, Vec<f64>) {
    let mut m: BTreeMap<i32, f64> = BTreeMap::new();
    for (j, v) in terms {
        *m.entry(j).or_insert(0.0) += v;
    }
    (m.keys().copied().collect(), m.values().copied().collect())
}

fn parse_cut(item: &Bound<'_, PyAny>, out: &mut Vec<UserCut>) -> PyResult<()> {
    if let Ok(c) = item.cast::<Constraint>() {
        let c = c.get();
        let terms: Vec<(i32, f64)> = c.body.with(|t| t.to_vec());
        if c.hi.is_finite() {
            let (indices, values) = merged(terms.iter().copied());
            out.push(UserCut {
                indices,
                values,
                rhs: c.hi,
            });
        }
        if c.lo.is_finite() {
            let (indices, values) = merged(terms.iter().map(|&(j, v)| (j, -v)));
            out.push(UserCut {
                indices,
                values,
                rhs: -c.lo,
            });
        }
        return Ok(());
    }
    let bad = || {
        PyTypeError::new_err(
            "a cut is a Constraint, (indices, values, rhs) or (coefficients, rhs) \
             with coefficients a dict {column: value} or a list of (column, value)",
        )
    };
    let tup: Vec<Bound<'_, PyAny>> = item.extract().map_err(|_| bad())?;
    match tup.len() {
        3 => {
            let indices: Vec<i32> = tup[0].extract()?;
            let values: Vec<f64> = tup[1].extract()?;
            out.push(UserCut {
                indices,
                values,
                rhs: tup[2].extract()?,
            });
        }
        2 => {
            let pairs: Vec<(i32, f64)> = match tup[0].cast::<PyDict>() {
                Ok(d) => d
                    .iter()
                    .map(|(k, v)| Ok((k.extract()?, v.extract()?)))
                    .collect::<PyResult<_>>()?,
                Err(_) => tup[0].extract().map_err(|_| bad())?,
            };
            let (indices, values) = merged(pairs);
            out.push(UserCut {
                indices,
                values,
                rhs: tup[1].extract()?,
            });
        }
        _ => return Err(bad()),
    }
    Ok(())
}

/// The `MipHooks` side of the bridge.
struct PyHooks(Arc<Bridge>);

impl MipHooks for PyHooks {
    fn name(&self) -> &str {
        "python"
    }

    fn branching_score(&mut self, info: &BranchingInfo<'_>) -> Option<Vec<f64>> {
        let b = &self.0;
        if !b.has.branching_score {
            return None;
        }
        let mut ctx = b.ctx("branching_score");
        ctx.depth = Some(info.depth as i64);
        ctx.bound = Some(info.node_lower_bound);
        ctx.solution = Some(info.solution.to_vec());
        ctx.lower = Some(info.col_lower.to_vec());
        ctx.upper = Some(info.col_upper.to_vec());
        ctx.candidates = Some(info.candidates.iter().map(|c| (c.col, c.value)).collect());
        let n = info.candidates.len();
        b.call("branching_score", ctx, |r| {
            if r.is_none() {
                return Ok(None);
            }
            let v: Vec<f64> = r.extract()?;
            if v.len() != n {
                return Err(PyValueError::new_err(format!(
                    "branching_score returned {} scores for {n} candidates",
                    v.len()
                )));
            }
            Ok(Some(v))
        })
        .flatten()
    }

    fn select_node(&mut self, open: &[OpenNodeInfo]) -> Option<usize> {
        let b = &self.0;
        if !b.has.select_node {
            return None;
        }
        let mut ctx = b.ctx("select_node");
        ctx.open_nodes = Some(
            open.iter()
                .map(|o| (o.lower_bound, o.estimate, o.depth, o.num_domain_changes))
                .collect(),
        );
        let n = open.len();
        b.call("select_node", ctx, |r| {
            if r.is_none() {
                return Ok(None);
            }
            let k: i64 = r.extract()?;
            if k < 0 || k as usize >= n {
                return Err(PyValueError::new_err(format!(
                    "select_node returned {k}, outside 0..{n}"
                )));
            }
            Ok(Some(k as usize))
        })
        .flatten()
    }

    fn separate(&mut self, info: &SeparationInfo<'_>) -> Vec<UserCut> {
        let b = &self.0;
        if !b.has.separate {
            return Vec::new();
        }
        let mut ctx = b.ctx("separate");
        ctx.objective = Some(info.lp_objective);
        ctx.nodes = Some(info.num_nodes);
        ctx.solution = Some(info.solution.to_vec());
        ctx.lower = Some(info.col_lower.to_vec());
        ctx.upper = Some(info.col_upper.to_vec());
        b.user_cuts("separate", ctx)
    }
}

/// The lazy-constraint oracle of `run_with_lazy_constraints`.
struct PyLazy(Arc<Bridge>);

impl LazyConstraints for PyLazy {
    fn lazy_constraints(&mut self, solution: &[f64]) -> Vec<UserCut> {
        let b = &self.0;
        let mut ctx = b.ctx("lazy_constraints");
        ctx.solution = Some(solution.to_vec());
        let cuts = b.user_cuts("lazy_constraints", ctx);
        if b.stop.load(Ordering::Acquire) {
            // An interrupt or a failure must not read as "accepted": a malformed row ends the
            // loop with `accepted = false` (docs/design/hooks.md, Lazy constraints).
            return vec![UserCut {
                indices: Vec::new(),
                values: Vec::new(),
                rhs: 0.0,
            }];
        }
        cuts
    }
}

fn out_ctx(b: &Bridge, event: &str, out: &HighsCallbackOutput) -> CallbackContext {
    let mut c = b.ctx(event);
    c.running_time = Some(out.running_time);
    c.nodes = Some(out.mip_node_count);
    c.lp_iterations = Some(out.mip_total_lp_iterations);
    let fin = |v: f64| v.is_finite().then_some(v);
    c.objective = fin(out.objective_function_value);
    c.bound = fin(out.mip_dual_bound);
    c.gap = fin(out.mip_gap);
    if !out.mip_solution.is_empty() {
        c.solution = Some(out.mip_solution.clone());
    }
    c
}

impl Bridge {
    /// Install the callback on `port`. Returns the bridge; `finish` removes everything again.
    pub(crate) fn install(
        py: Python<'_>,
        cb: &Bound<'_, PyAny>,
        port: &mut Port,
    ) -> PyResult<Arc<Bridge>> {
        let has = Has::of(cb);
        let any = has.on_progress
            || has.on_improving_solution
            || has.on_cut_pool
            || has.user_solution
            || has.presolved_space()
            || has.lazy_constraints;
        if !any {
            return Err(crate::diag::model_err(
                "the callback object defines no hook: define one or more of on_progress, \
                 on_improving_solution, on_cut_pool, user_solution, separate, branching_score, \
                 select_node, lazy_constraints",
            ));
        }
        let _ = py;
        let num_col = port.model_lp_.num_col_ as usize;
        let b = Arc::new(Bridge {
            cb: cb.clone().unbind(),
            has,
            failure: Mutex::new(None),
            errors: AtomicU64::new(0),
            calls: Mutex::new(BTreeMap::new()),
            stop: Arc::new(AtomicBool::new(false)),
            num_col,
        });
        let me = Arc::clone(&b);
        port.set_callback(HighsCallbackFunctionType::new(
            move |ty, _msg, out, data_in| {
                let b = &*me;
                match ty {
                    t if t == kCallbackMipInterrupt => {
                        if b.has.on_progress {
                            let ctx = out_ctx(b, "on_progress", out);
                            b.call("on_progress", ctx, |_| Ok(()));
                        }
                        if b.stop.load(Ordering::Acquire) {
                            data_in.user_interrupt = true;
                        }
                    }
                    t if t == kCallbackMipImprovingSolution && b.has.on_improving_solution => {
                        let ctx = out_ctx(b, "on_improving_solution", out);
                        b.call("on_improving_solution", ctx, |_| Ok(()));
                    }
                    t if t == kCallbackMipGetCutPool && b.has.on_cut_pool => {
                        let mut ctx = out_ctx(b, "on_cut_pool", out);
                        let rows = out.cutpool_num_cut.max(0) as usize;
                        let mut pool = Vec::with_capacity(rows);
                        for r in 0..rows {
                            let (s, e) = (
                                out.cutpool_start[r] as usize,
                                out.cutpool_start
                                    .get(r + 1)
                                    .map_or(out.cutpool_index.len(), |&x| x as usize),
                            );
                            pool.push((
                                out.cutpool_index[s..e].to_vec(),
                                out.cutpool_value[s..e].to_vec(),
                                out.cutpool_lower[r],
                                out.cutpool_upper[r],
                            ));
                        }
                        ctx.cut_pool = Some(pool);
                        b.call("on_cut_pool", ctx, |_| Ok(()));
                    }
                    t if t == kCallbackMipUserSolution && b.has.user_solution => {
                        let ctx = b.ctx("user_solution");
                        let n = b.num_col;
                        let sol = b.call("user_solution", ctx, |r| {
                            if r.is_none() {
                                return Ok(None);
                            }
                            let v: Vec<f64> = r.extract()?;
                            if v.len() != n {
                                return Err(PyValueError::new_err(format!(
                                    "user_solution returned {} values for {n} columns",
                                    v.len()
                                )));
                            }
                            Ok(Some(v))
                        });
                        if let Some(Some(v)) = sol {
                            data_in.user_has_solution = true;
                            data_in.user_solution = v;
                        }
                    }
                    _ => {}
                }
                if b.stop.load(Ordering::Acquire) {
                    data_in.user_interrupt = true;
                }
            },
        ));
        port.start_callback(kCallbackMipInterrupt);
        if b.has.on_improving_solution {
            port.start_callback(kCallbackMipImprovingSolution);
        }
        if b.has.on_cut_pool {
            port.start_callback(kCallbackMipGetCutPool);
        }
        if b.has.user_solution {
            port.start_callback(kCallbackMipUserSolution);
        }
        if b.has.presolved_space() {
            port.set_mip_hooks(Box::new(PyHooks(Arc::clone(&b))));
        }
        Ok(b)
    }

    /// True when a hook needs the model's own column space (presolve off).
    pub(crate) fn needs_presolve_off(cb: &Bound<'_, PyAny>) -> bool {
        Has::of(cb).presolved_space()
    }

    /// Run the lazy-constraint loop (or a plain run when the callback has no lazy hook).
    pub(crate) fn run(
        self: &Arc<Self>,
        port: &mut Port,
        lazy_rounds: usize,
    ) -> (
        rusolve_ekk::lp_data::highs_status::HighsStatus,
        Option<LazyReport>,
    ) {
        if self.has.lazy_constraints {
            let mut oracle = PyLazy(Arc::clone(self));
            let (st, rep) = port.run_with_lazy_constraints(&mut oracle, lazy_rounds);
            (st, Some(rep))
        } else {
            (port.run(), None)
        }
    }

    /// Remove the callback and hooks from `port`; the counters as a dict for `Model.callback_stats`.
    pub(crate) fn finish(
        &self,
        py: Python<'_>,
        port: &mut Port,
        lazy: Option<&LazyReport>,
    ) -> PyResult<Py<PyDict>> {
        let hs = port.mip_hook_stats();
        port.clear_mip_hooks();
        let n = port.model_lp_.num_col_ as usize;
        port.callback_.clear(n);
        port.options_.log_options.user_callback = None;
        let d = PyDict::new(py);
        let calls = PyDict::new(py);
        for (k, v) in lock(&self.calls).iter() {
            calls.set_item(k, v)?;
        }
        d.set_item("calls", calls)?;
        d.set_item("errors", self.errors.load(Ordering::Relaxed))?;
        d.set_item("interrupted", self.stop.load(Ordering::Acquire))?;
        d.set_item("cuts_offered", hs.user_cuts_offered)?;
        d.set_item("cuts_accepted", hs.user_cuts_accepted)?;
        d.set_item("cuts_duplicate", hs.user_cuts_duplicate)?;
        d.set_item("cuts_rejected", hs.user_cuts_rejected_malformed)?;
        d.set_item("branching_overrides", hs.branching_overrides)?;
        d.set_item("node_select_overrides", hs.node_select_overrides)?;
        d.set_item(
            "user_solutions_accepted",
            hs.user_solutions_checked_feasible,
        )?;
        d.set_item("user_solutions_rejected", hs.user_solutions_rejected)?;
        if let Some(l) = lazy {
            let z = PyDict::new(py);
            z.set_item("rounds", l.rounds)?;
            z.set_item("rows_added", l.rows_added)?;
            z.set_item("rows_rejected", l.rows_rejected_malformed)?;
            z.set_item("rows_not_violated", l.rows_not_violated)?;
            z.set_item("accepted", l.accepted)?;
            d.set_item("lazy", z)?;
        }
        Ok(d.unbind())
    }

    pub(crate) fn interrupted(&self) -> bool {
        self.stop.load(Ordering::Acquire) && !self.failed()
    }

    /// The typed error for the first failure (None when no hook failed). The original exception
    /// is the `__cause__`.
    pub(crate) fn error(&self, py: Python<'_>) -> Option<PyErr> {
        let f = lock(&self.failure).take()?;
        let n = self.errors.load(Ordering::Relaxed);
        let ty = f
            .err
            .get_type(py)
            .name()
            .map(|s| s.to_string())
            .unwrap_or_else(|_| "exception".into());
        let msg = format!(
            "the Python callback hook `{}` raised {ty}: {} ({n} hook error{}; the solve was \
             stopped and no later hook was called)",
            f.hook,
            f.err.value(py),
            if n == 1 { "" } else { "s" },
        );
        let d = Diagnostic::error("RS4008", msg)
            .with_help("the original exception is `__cause__`; the partial result is `.solution`")
            .with_help(
                "fix the hook (a hook must return the type its docstring names) and solve again",
            );
        let e = crate::diag::raise(py, "CallbackError", d);
        e.set_cause(py, Some(f.err));
        Some(e)
    }
}
