//! `rusolve.Solution` (docs/design/python-api.md, "solution").
//!
//! A snapshot of one `Model.solve()`: the vectors are copied out of the solver, so later
//! model changes and re-solves never alter it. Arrays become numpy float64 arrays on first
//! access (one `numpy.empty` plus one buffer copy, no rust-numpy dependency) and are cached;
//! without numpy they are lists.

use pyo3::buffer::PyBuffer;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyFloat, PyInt, PyList};

use crate::expr::{Constraint, LinExpr, MIXED_MODELS, Var, Vars};
use crate::ingest;

/// The solve results. Built by `Model::solve`.
#[pyclass(frozen, module = "rusolve")]
pub struct Solution {
    pub(crate) status: &'static str,
    pub(crate) objective: Option<f64>,
    pub(crate) bound: Option<f64>,
    pub(crate) gap: Option<f64>,
    pub(crate) x: Option<Vec<f64>>,
    pub(crate) row_activity: Option<Vec<f64>>,
    pub(crate) duals: Option<Vec<f64>>,
    pub(crate) reduced_costs: Option<Vec<f64>>,
    pub(crate) nodes: i64,
    pub(crate) iterations: i64,
    pub(crate) time: f64,
    pub(crate) model: u64,
    pub(crate) cache: [PyOnceLock<Py<PyAny>>; 4],
    /// The `rusolve.result/1` envelope of this solve and the column names it reports.
    pub(crate) env: Option<rusolve_cli::api::EnvParts>,
    pub(crate) names: std::sync::Arc<Vec<String>>,
    /// The model has a quadratic objective: the LP-only ray is not a ray of the QP, so no RS2011.
    pub(crate) quad: bool,
    /// The `general_constraints` object of the solve (objects, routes, violations, `verified`),
    /// present when the model has indicator, SOS, piecewise-linear or function constraints.
    pub(crate) gc_report: Option<serde_json::Value>,
    /// The `cp` object of the solve (engine, objects, routes, violations, `verified`), present when
    /// the model has scheduling or combinatorial objects.
    pub(crate) cp_report: Option<serde_json::Value>,
    /// The model as solved, kept only for an infeasible or unbounded verdict: `raise_for_status`
    /// computes the infeasible set or the ray from it, at raise time and not on every solve.
    pub(crate) lp: Option<std::sync::Arc<rusolve_ekk::prelude::HighsLp>>,
}

/// `numpy.empty`, or None when numpy cannot be imported.
static NUMPY_EMPTY: PyOnceLock<Option<Py<PyAny>>> = PyOnceLock::new();

/// A float64 numpy array with the values of `v` (a list without numpy).
pub(crate) fn to_array(py: Python<'_>, v: &[f64]) -> PyResult<Py<PyAny>> {
    let empty = NUMPY_EMPTY.get_or_init(py, || {
        py.import("numpy")
            .and_then(|m| m.getattr("empty"))
            .map(Bound::unbind)
            .ok()
    });
    match empty {
        Some(f) => {
            let a = f.bind(py).call1((v.len(), "float64"))?;
            PyBuffer::<f64>::get(&a)?.copy_from_slice(py, v)?;
            Ok(a.unbind())
        }
        None => Ok(PyList::new(py, v)?.into_any().unbind()),
    }
}

impl Solution {
    fn envelope_value(&self, include_solution: bool) -> PyResult<serde_json::Value> {
        let env = self
            .env
            .as_ref()
            .ok_or_else(|| crate::diag::model_err("this solution carries no result envelope"))?;
        let values = if include_solution {
            self.x.as_deref()
        } else {
            None
        };
        let mut v = rusolve_cli::api::envelope_value(&env.envelope(), &self.names, values)
            .map_err(crate::diag::model_err)?;
        if let Some(r) = &self.gc_report {
            v["general_constraints"] = r.clone();
        }
        if let Some(r) = &self.cp_report {
            v["cp"] = r.clone();
        }
        Ok(v)
    }

    fn array(&self, py: Python<'_>, k: usize, v: &Option<Vec<f64>>) -> PyResult<Py<PyAny>> {
        match v {
            None => Ok(py.None()),
            Some(v) => self.cache[k]
                .get_or_try_init(py, || to_array(py, v))
                .map(|a| a.clone_ref(py)),
        }
    }

    fn primal(&self) -> PyResult<&[f64]> {
        self.x
            .as_deref()
            .ok_or_else(|| crate::diag::no_solution(self.status))
    }

    fn col(&self, x: &[f64], j: i64) -> PyResult<f64> {
        if j < 0 || j as usize >= x.len() {
            return Err(crate::diag::index_err(format!(
                "variable index {j} is out of range (the solution has {} variables)",
                x.len()
            )));
        }
        Ok(x[j as usize])
    }

    fn check_model(&self, m: u64) -> PyResult<()> {
        if m != 0 && m != self.model {
            Err(crate::diag::model_err(format!(
                "{MIXED_MODELS}: this belongs to another Model than the Solution"
            )))
        } else {
            Ok(())
        }
    }

    fn terms_value(&self, x: &[f64], t: &[(i32, f64)]) -> PyResult<f64> {
        let mut s = 0.0;
        for &(j, c) in t {
            s += c * self.col(x, j as i64)?;
        }
        Ok(s)
    }

    /// The value of a Var, int index, LinExpr, Constraint body or dict row.
    fn value_of(&self, key: &Bound<'_, PyAny>) -> PyResult<f64> {
        let x = self.primal()?;
        if let Ok(v) = key.cast_exact::<Var>() {
            let v = v.get();
            self.check_model(v.model)?;
            return self.col(x, v.index as i64);
        }
        if key.is_instance_of::<PyInt>() {
            return self.col(x, key.extract()?);
        }
        if let Ok(e) = key.cast_exact::<LinExpr>() {
            let e = e.try_borrow()?;
            self.check_model(e.model)?;
            return Ok(e.terms.with(|t| self.terms_value(x, t))? + e.constant);
        }
        if let Ok(c) = key.cast_exact::<Constraint>() {
            let c = c.get();
            self.check_model(c.model)?;
            return c.body.with(|t| self.terms_value(x, t));
        }
        if let Ok(d) = key.cast::<PyDict>() {
            let mut s = 0.0;
            for (k, v) in d.iter() {
                let (j, m) = ingest::var_key(&k)?;
                self.check_model(m)?;
                s += ingest::coef(&v)? * self.col(x, j as i64)?;
            }
            return Ok(s);
        }
        if key.is_instance_of::<PyFloat>() {
            return Err(PyTypeError::new_err("a variable index must be an int"));
        }
        if let Ok(j) = key.extract::<i64>() {
            return self.col(x, j);
        }
        Err(PyTypeError::new_err(format!(
            "sol[...] takes a Var, an int, a LinExpr, a Constraint, a dict row or a Vars, got {}",
            ingest::type_name(key)
        )))
    }
}

#[pymethods]
impl Solution {
    /// "Optimal", "Infeasible", "Unbounded", "InfeasibleOrUnbounded", "TimeLimit", ...
    #[getter]
    fn status(&self) -> &'static str {
        self.status
    }
    /// Raise the typed error of a non-optimal status (`InfeasibleError`, `UnboundedError`,
    /// `LimitError`, each with `.diagnostics`); return `self` when the status is Optimal.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_var(); _ = m.add_constr(x >= 3); _ = m.add_constr(x <= 1)
    /// >>> try:
    /// ...     m.solve().raise_for_status()
    /// ... except rs.InfeasibleError as e:
    /// ...     e.code
    /// 'RS4004'
    fn raise_for_status(slf: PyRef<'_, Self>) -> PyResult<PyRef<'_, Self>> {
        let py = slf.py();
        let post = match &slf.lp {
            // the CLI's bounded post-solve diagnostics (RS2010 / RS2011), computed now
            Some(lp) => {
                let (lp, status, time, quad) = (lp.clone(), slf.status, slf.time, slf.quad);
                py.detach(move || crate::diag::post_solve(&lp, status, time, quad))
            }
            None => Vec::new(),
        };
        match crate::diag::status_err(slf.status, post) {
            None => Ok(slf),
            Some(e) => Err(e),
        }
    }
    /// True if a primal solution is available (`x` is not None).
    #[getter]
    fn ok(&self) -> bool {
        self.x.is_some()
    }
    #[getter]
    fn objective(&self) -> Option<f64> {
        self.objective
    }
    #[getter]
    fn bound(&self) -> Option<f64> {
        self.bound
    }
    #[getter]
    fn gap(&self) -> Option<f64> {
        self.gap
    }
    #[getter]
    fn x(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.array(py, 0, &self.x)
    }
    #[getter]
    fn row_activity(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.array(py, 1, &self.row_activity)
    }
    /// Row duals (LP only).
    #[getter]
    fn duals(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.array(py, 2, &self.duals)
    }
    /// Column duals (LP only).
    #[getter]
    fn reduced_costs(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.array(py, 3, &self.reduced_costs)
    }
    /// For a model with general constraints: `{"objects", "routes", "added_columns",
    /// "added_rows", "violations", "verified"}`, checked on your objects (not on the generated
    /// rows). None for a model without them.
    #[getter]
    fn general_constraints(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match &self.gc_report {
            None => Ok(py.None()),
            Some(v) => {
                let text =
                    serde_json::to_string(v).map_err(|e| crate::diag::model_err(e.to_string()))?;
                Ok(py.import("json")?.call_method1("loads", (text,))?.unbind())
            }
        }
    }
    /// For a model with scheduling or combinatorial objects: `{"engine", "objects", "checked",
    /// "violations", "verified", "routes", ...}`, checked on your objects (not on the generated
    /// rows). None for a model without them.
    #[getter]
    fn cp(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match &self.cp_report {
            None => Ok(py.None()),
            Some(v) => {
                let text =
                    serde_json::to_string(v).map_err(|e| crate::diag::model_err(e.to_string()))?;
                Ok(py.import("json")?.call_method1("loads", (text,))?.unbind())
            }
        }
    }
    #[getter]
    fn nodes(&self) -> i64 {
        self.nodes
    }
    #[getter]
    fn iterations(&self) -> i64 {
        self.iterations
    }
    /// Wall time of the solve in seconds.
    #[getter]
    fn time(&self) -> f64 {
        self.time
    }

    /// The value of a Var / int index / LinExpr / Constraint body / dict row; a `Vars` gives
    /// an array of its values.
    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(vs) = key.cast_exact::<Vars>() {
            let vs = vs.get();
            self.check_model(vs.model)?;
            let x = self.primal()?;
            let mut v = Vec::with_capacity(vs.len);
            for k in 0..vs.len {
                v.push(self.col(x, vs.col(k) as i64)?);
            }
            return to_array(py, &v);
        }
        Ok(PyFloat::new(py, self.value_of(key)?).into_any().unbind())
    }

    /// `{key: value}` for a `{key: var}` mapping (a list of values for a sequence), keeping
    /// only values `> above` when `above` is given.
    #[pyo3(signature = (mapping, above=None))]
    fn values(
        &self,
        py: Python<'_>,
        mapping: &Bound<'_, PyAny>,
        above: Option<f64>,
    ) -> PyResult<Py<PyAny>> {
        let keep = |v: f64| above.is_none_or(|a| v > a);
        if let Ok(d) = mapping.cast::<PyDict>() {
            let out = PyDict::new(py);
            for (k, var) in d.iter() {
                let v = self.value_of(&var)?;
                if keep(v) {
                    out.set_item(k, v)?;
                }
            }
            return Ok(out.into_any().unbind());
        }
        let out = PyList::empty(py);
        for item in mapping.try_iter()? {
            let v = self.value_of(&item?)?;
            if keep(v) {
                out.append(v)?;
            }
        }
        Ok(out.into_any().unbind())
    }

    /// The keys of a `{key: var}` mapping (the items of a sequence) whose value is
    /// `> threshold`: the picks of binary variables.
    #[pyo3(signature = (mapping, threshold=0.5))]
    fn selected(
        &self,
        py: Python<'_>,
        mapping: &Bound<'_, PyAny>,
        threshold: f64,
    ) -> PyResult<Py<PyList>> {
        let out = PyList::empty(py);
        if let Ok(d) = mapping.cast::<PyDict>() {
            for (k, var) in d.iter() {
                if self.value_of(&var)? > threshold {
                    out.append(k)?;
                }
            }
        } else {
            for item in mapping.try_iter()? {
                let item = item?;
                if self.value_of(&item)? > threshold {
                    out.append(item)?;
                }
            }
        }
        Ok(out.unbind())
    }

    /// The result as the `rusolve.result/1` JSON envelope the CLI prints for `--format json`
    /// (schema `docs/schema/rusolve-result-1.schema.json`), as a string. With
    /// `include_solution=True` (the default here, unlike the CLI) the object carries
    /// `solution.columns`, one `{name, value}` per variable in column order, when a primal
    /// solution exists. `indent` pretty-prints.
    ///
    /// >>> import rusolve as rs, json
    /// >>> m = rs.Model(); x = m.add_vars(2, ub=1); m.maximize(x[0] + x[1])
    /// >>> env = json.loads(m.solve().to_json())
    /// >>> env["schema_version"], env["status"], env["objective"]
    /// ('rusolve.result/1', 'Optimal', 2.0)
    #[pyo3(signature = (include_solution=true, indent=None))]
    fn to_json(&self, include_solution: bool, indent: Option<usize>) -> PyResult<String> {
        let v = self.envelope_value(include_solution)?;
        let text = match indent {
            None => serde_json::to_string(&v),
            Some(n) => {
                let pad = vec![b' '; n];
                let fmt = serde_json::ser::PrettyFormatter::with_indent(&pad);
                let mut out = Vec::new();
                let mut ser = serde_json::Serializer::with_formatter(&mut out, fmt);
                serde::Serialize::serialize(&v, &mut ser)
                    .map(|_| String::from_utf8(out).unwrap_or_default())
            }
        };
        text.map_err(|e| crate::diag::model_err(e.to_string()))
    }

    /// [`to_json`](Solution.to_json) as a dict.
    ///
    /// >>> import rusolve as rs
    /// >>> m = rs.Model(); x = m.add_var(ub=3); m.maximize(x)
    /// >>> m.solve().to_dict(include_solution=False)["objective"]
    /// 3.0
    #[pyo3(signature = (include_solution=true))]
    fn to_dict(&self, py: Python<'_>, include_solution: bool) -> PyResult<Py<PyAny>> {
        let v = self.envelope_value(include_solution)?;
        let text = serde_json::to_string(&v).map_err(|e| crate::diag::model_err(e.to_string()))?;
        Ok(py.import("json")?.call_method1("loads", (text,))?.unbind())
    }

    fn __repr__(&self) -> String {
        let f = |v: Option<f64>| v.map_or("None".to_string(), |v| format!("{v}"));
        let mut s = format!(
            "Solution(status={:?}, objective={}",
            self.status,
            f(self.objective)
        );
        if self.nodes > 0 || self.gap.is_some_and(|g| g != 0.0) {
            s += &format!(
                ", bound={}, gap={}, nodes={}",
                f(self.bound),
                f(self.gap),
                self.nodes
            );
        }
        s += &format!(", iterations={}, time={:.3}s)", self.iterations, self.time);
        s
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Solution>()?;
    Ok(())
}
