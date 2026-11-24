//! Scheduling and combinatorial (CP) objects of `rusolve.Model` (docs/design/cp-model.md):
//! intervals, no-overlap, cumulative, all-different, element and table.
//!
//! The model keeps the *user's* objects in `Model::cp` (the typed objects of `rusolve-io`, the
//! same ones the MCP `build_model` tool fills). A solve of a model that has any of them goes
//! through `rusolve_cli::cp_solve::solve_with`, which picks the engine (today the MIP
//! linearization), checks the answer on the objects and cuts it back to the user's columns. The
//! model's own solver state is never changed, and a model without these objects never reaches
//! this module (`Model::solve` tests `CpModel::is_empty`).

use std::time::Instant;

use pyo3::prelude::*;
use pyo3::types::PyInt;
use rusolve_io::{
    AllDifferent, CpModel, CpTerm, Cumulative, Element, Interval as IoInterval, NoOverlap, Table,
    VariableBound,
};

use super::Model;
use crate::diag::{cp_err, cp_error, cp_solve_error, cp_unsupported};
use crate::expr::{MIXED_MODELS, Var};
use crate::ingest;
use crate::solution::Solution;

/// An interval of a `Model` (`Model.add_interval`): `start + size = end` while it is present.
/// The handle that `add_no_overlap` and `add_cumulative` take.
#[pyclass(frozen, module = "rusolve")]
pub struct Interval {
    pub(crate) index: usize,
    pub(crate) name: String,
    start: i32,
    end: i32,
    size: CpTerm,
    presence: Option<i32>,
    model: u64,
}

#[pymethods]
impl Interval {
    /// The interval's name.
    #[getter]
    fn name(&self) -> &str {
        &self.name
    }

    /// Its position among the model's intervals.
    #[getter]
    fn index(&self) -> usize {
        self.index
    }

    /// The start variable.
    #[getter]
    fn start(&self) -> Var {
        Var::new(self.start, self.model)
    }

    /// The end variable (the one `add_interval` made when `end` was not given).
    #[getter]
    fn end(&self) -> Var {
        Var::new(self.end, self.model)
    }

    /// The size: an int, or the Var that holds it.
    #[getter]
    fn size(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match self.size {
            CpTerm::Const(c) => Ok(c.into_pyobject(py)?.into_any().unbind()),
            CpTerm::Col(j) => Ok(Py::new(py, Var::new(j as i32, self.model))?.into_any()),
        }
    }

    /// The presence binary, or None for an interval that is always present.
    #[getter]
    fn presence(&self) -> Option<Var> {
        self.presence.map(|p| Var::new(p, self.model))
    }

    fn __repr__(&self) -> String {
        format!("Interval({:?})", self.name)
    }
}

impl Model {
    /// `(lower, upper)` of an integer column with finite bounds, else RS4011 naming `obj`.
    fn cp_need_int(&self, obj: &str, what: &str, j: usize) -> PyResult<(f64, f64)> {
        let (lo, hi, integer) = self.gc_col_info(j);
        let name = self.gc_col_label(j);
        if !integer {
            return Err(cp_err(format!(
                "{obj}: {what}: variable '{name}' is not integer (CP objects are defined over integer variables; declare it with `add_var(integer=True)`)"
            )));
        }
        if !(lo.is_finite() && hi.is_finite()) {
            return Err(cp_err(format!(
                "{obj}: {what}: variable '{name}' needs finite bounds (it has [{lo}, {hi}])"
            )));
        }
        Ok((lo, hi))
    }

    /// An integer constant or a Var of this model.
    fn cp_term(&self, o: &Bound<'_, PyAny>, obj: &str, what: &str) -> PyResult<CpTerm> {
        if o.is_instance_of::<Var>() {
            return Ok(CpTerm::Col(self.col_index(o)?));
        }
        match ingest::number(o)? {
            Some(v) if v.is_finite() && v.fract() == 0.0 && v.abs() < 9e15 => {
                Ok(CpTerm::Const(v as i64))
            }
            _ => Err(cp_err(format!(
                "{obj}: {what} must be an integer or an integer Var, got {}",
                ingest::type_name(o)
            ))),
        }
    }

    fn cp_terms(&self, o: &Bound<'_, PyAny>, obj: &str, what: &str) -> PyResult<Vec<CpTerm>> {
        o.try_iter()
            .map_err(|_| {
                cp_err(format!(
                    "{obj}: {what} must be a sequence, got {}",
                    ingest::type_name(o)
                ))
            })?
            .map(|x| self.cp_term(&x?, obj, what))
            .collect()
    }

    /// The intervals of `o` (an iterable of `Interval` handles of this model) as indices.
    fn cp_intervals(&self, o: &Bound<'_, PyAny>, obj: &str) -> PyResult<Vec<usize>> {
        let it = o.try_iter().map_err(|_| {
            cp_err(format!(
                "{obj}: intervals must be a sequence of the Interval objects `add_interval` returns, got {}",
                ingest::type_name(o)
            ))
        })?;
        let mut out = Vec::new();
        for x in it {
            let x = x?;
            let Ok(iv) = x.cast_exact::<Interval>() else {
                return Err(cp_err(format!(
                    "{obj}: expected an Interval from `add_interval`, got {}",
                    ingest::type_name(&x)
                )));
            };
            let iv = iv.get();
            if iv.model != self.id {
                return Err(crate::diag::model_err(format!(
                    "{MIXED_MODELS}: interval '{}' belongs to another Model",
                    iv.name
                )));
            }
            out.push(iv.index);
        }
        Ok(out)
    }

    /// The object name: `requested` if it is free, else `<prefix><k>`.
    fn cp_name(
        requested: Option<String>,
        prefix: &str,
        count: usize,
        used: impl Fn(&str) -> bool,
    ) -> PyResult<String> {
        match requested {
            Some(n) => {
                if n.is_empty() || n.chars().any(char::is_whitespace) {
                    return Err(cp_err(format!(
                        "name {n:?}: a CP object name is not empty and has no whitespace"
                    )));
                }
                if used(&n) {
                    return Err(cp_err(format!(
                        "name {n:?} is already used by another object of the same kind"
                    )));
                }
                Ok(n)
            }
            None => {
                let mut k = count;
                loop {
                    let n = format!("{prefix}{k}");
                    if !used(&n) {
                        return Ok(n);
                    }
                    k += 1;
                }
            }
        }
    }

    /// Structural validation of the objects against the model as it is now.
    fn cp_validate_now(&self) -> PyResult<()> {
        let n = self.b.num_cols();
        let bounds: Vec<VariableBound> = (0..n)
            .map(|j| {
                let (lower, upper, integer) = self.gc_col_info(j);
                VariableBound {
                    lower,
                    upper,
                    integer,
                }
            })
            .collect();
        let names: Vec<String> = (0..n).map(|j| self.gc_col_label(j)).collect();
        self.cp.validate(&bounds, &names).map_err(|e| cp_error(&e))
    }

    /// Run `validate` after `add` pushed an object; `undo` takes it back when it is invalid.
    fn cp_commit(&mut self, undo: impl FnOnce(&mut CpModel)) -> PyResult<()> {
        if let Err(e) = self.cp_validate_now() {
            undo(&mut self.cp);
            return Err(e);
        }
        Ok(())
    }

    // ------------------------------------------------------------ the add_* methods

    pub(super) fn cp_interval(
        &mut self,
        start: &Bound<'_, PyAny>,
        size: &Bound<'_, PyAny>,
        end: Option<&Bound<'_, PyAny>>,
        presence: Option<&Bound<'_, PyAny>>,
        name: Option<String>,
    ) -> PyResult<Interval> {
        let used = |n: &str| self.cp.intervals.iter().any(|o| o.name == n);
        let name = Self::cp_name(name, "interval", self.cp.intervals.len(), used)?;
        let obj = format!("interval '{name}'");
        let s = self.col_index(start)?;
        let size_t = self.cp_term(size, &obj, "size")?;
        let presence = match presence.filter(|p| !p.is_none()) {
            Some(p) => Some(self.col_index(p)?),
            None => None,
        };
        let e = match end.filter(|e| !e.is_none()) {
            Some(e) => self.col_index(e)?,
            None => {
                // the end column the interval needs: integer, between the smallest and the
                // largest start + size
                let (slo, shi) = self.cp_need_int(&obj, "start", s)?;
                let (zlo, zhi) = match size_t {
                    CpTerm::Const(c) => {
                        if c < 0 {
                            return Err(cp_err(format!(
                                "{obj}: size is {c}; it must be nonnegative"
                            )));
                        }
                        (c as f64, c as f64)
                    }
                    CpTerm::Col(j) => self.cp_need_int(&obj, "size", j)?,
                };
                let v = self.add_var(
                    Some(slo + zlo),
                    Some(shi + zhi),
                    0.0,
                    true,
                    false,
                    Some(format!("{name}_end")),
                )?;
                v.index as usize
            }
        };
        self.cp.intervals.push(IoInterval {
            name: name.clone(),
            start: s,
            size: size_t,
            end: e,
            presence,
        });
        let index = self.cp.intervals.len() - 1;
        self.cp_commit(|cp| {
            cp.intervals.pop();
        })?;
        Ok(Interval {
            index,
            name,
            start: s as i32,
            end: e as i32,
            size: size_t,
            presence: presence.map(|p| p as i32),
            model: self.id,
        })
    }

    pub(super) fn cp_no_overlap(
        &mut self,
        intervals: &Bound<'_, PyAny>,
        name: Option<String>,
    ) -> PyResult<String> {
        let used = |n: &str| self.cp.no_overlap.iter().any(|o| o.name == n);
        let name = Self::cp_name(name, "no_overlap", self.cp.no_overlap.len(), used)?;
        let obj = format!("no_overlap '{name}'");
        let list = self.cp_intervals(intervals, &obj)?;
        self.cp.no_overlap.push(NoOverlap {
            name: name.clone(),
            intervals: list,
        });
        self.cp_commit(|cp| {
            cp.no_overlap.pop();
        })?;
        Ok(name)
    }

    pub(super) fn cp_cumulative(
        &mut self,
        intervals: &Bound<'_, PyAny>,
        demands: &Bound<'_, PyAny>,
        capacity: &Bound<'_, PyAny>,
        name: Option<String>,
    ) -> PyResult<String> {
        let used = |n: &str| self.cp.cumulative.iter().any(|o| o.name == n);
        let name = Self::cp_name(name, "cumulative", self.cp.cumulative.len(), used)?;
        let obj = format!("cumulative '{name}'");
        let list = self.cp_intervals(intervals, &obj)?;
        // a single number or Var is the demand of every interval
        let demands = if demands.is_instance_of::<Var>() || demands.is_instance_of::<PyInt>() {
            vec![self.cp_term(demands, &obj, "demand")?; list.len()]
        } else {
            self.cp_terms(demands, &obj, "demand")?
        };
        let capacity = self.cp_term(capacity, &obj, "capacity")?;
        self.cp.cumulative.push(Cumulative {
            name: name.clone(),
            intervals: list,
            demands,
            capacity,
        });
        self.cp_commit(|cp| {
            cp.cumulative.pop();
        })?;
        Ok(name)
    }

    pub(super) fn cp_all_different(
        &mut self,
        vars: &Bound<'_, PyAny>,
        name: Option<String>,
    ) -> PyResult<String> {
        let used = |n: &str| self.cp.all_different.iter().any(|o| o.name == n);
        let name = Self::cp_name(name, "all_different", self.cp.all_different.len(), used)?;
        let cols = self.gc_cols(vars, "vars")?;
        self.cp.all_different.push(AllDifferent {
            name: name.clone(),
            cols,
        });
        self.cp_commit(|cp| {
            cp.all_different.pop();
        })?;
        Ok(name)
    }

    pub(super) fn cp_element(
        &mut self,
        index: &Bound<'_, PyAny>,
        values: &Bound<'_, PyAny>,
        target: &Bound<'_, PyAny>,
        name: Option<String>,
    ) -> PyResult<String> {
        let used = |n: &str| self.cp.element.iter().any(|o| o.name == n);
        let name = Self::cp_name(name, "element", self.cp.element.len(), used)?;
        let obj = format!("element '{name}'");
        let index = self.col_index(index)?;
        let target = self.col_index(target)?;
        let values = self.cp_terms(values, &obj, "value")?;
        self.cp.element.push(Element {
            name: name.clone(),
            index,
            values,
            target,
        });
        self.cp_commit(|cp| {
            cp.element.pop();
        })?;
        Ok(name)
    }

    pub(super) fn cp_table(
        &mut self,
        vars: &Bound<'_, PyAny>,
        tuples: &Bound<'_, PyAny>,
        negated: bool,
        name: Option<String>,
    ) -> PyResult<String> {
        let used = |n: &str| self.cp.table.iter().any(|o| o.name == n);
        let name = Self::cp_name(name, "table", self.cp.table.len(), used)?;
        let obj = format!("table '{name}'");
        let cols = self.gc_cols(vars, "vars")?;
        let rows = tuples.try_iter().map_err(|_| {
            cp_err(format!(
                "{obj}: tuples must be a sequence of sequences of integers, got {}",
                ingest::type_name(tuples)
            ))
        })?;
        let mut list = Vec::new();
        for row in rows {
            let row = row?;
            let mut t = Vec::new();
            for v in self.cp_terms(&row, &obj, "tuple value")? {
                match v {
                    CpTerm::Const(c) => t.push(c),
                    CpTerm::Col(_) => {
                        return Err(cp_err(format!(
                            "{obj}: tuple values are integer constants, not variables"
                        )));
                    }
                }
            }
            list.push(t);
        }
        self.cp.table.push(Table {
            name: name.clone(),
            cols,
            tuples: list,
            negated,
        });
        self.cp_commit(|cp| {
            cp.table.pop();
        })?;
        Ok(name)
    }

    // ------------------------------------------------------------ solve and refusals

    /// RS4012 when the model has CP objects and `what` would ignore them.
    pub(super) fn refuse_cp(&self, what: &str) -> PyResult<()> {
        match self.cp.first_object() {
            Some(o) => Err(cp_unsupported(what, &o)),
            None => Ok(()),
        }
    }

    /// Solve the linearization of the model on a separate solver, check the answer on the user's
    /// objects, and report on the user's columns and rows. Options are those `solve()` just
    /// applied. The engine is picked by `cp_solve::select_engine`: the linearization, or
    /// OR-Tools CP-SAT with `mip_engine=cpsat`.
    pub(super) fn solve_cp(&mut self, py: Python<'_>) -> PyResult<Solution> {
        let t0 = Instant::now();
        let port = self.port.as_ref().expect("synced");
        let user = port.model_lp_.clone();
        let opts = port.get_options().clone();
        let solved = rusolve_cli::cp_solve::solve_with(user, self.cp.clone(), opts, |h| {
            py.detach(|| {
                h.reset_scheduler_for_threads();
                h.run()
            })
        })
        .map_err(|e| cp_solve_error(&e))?;
        let time = t0.elapsed().as_secs_f64();
        Ok(self.solution_from(
            &solved.highs,
            Some(&solved.user_lp),
            None,
            Some((solved.report, solved.engine)),
            solved.run,
            time,
        ))
    }
}
