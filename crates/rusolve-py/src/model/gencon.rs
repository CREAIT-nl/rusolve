//! General constraints of `rusolve.Model` (docs/design/python-api.md, "general constraints"):
//! indicator, SOS1/SOS2, piecewise-linear, abs/min/max and and/or/not.
//!
//! The model keeps the *user's* objects in `Model::gc` (the typed objects of `rusolve-io`, the
//! same ones the MPS and LP readers fill, so there is no second representation). A solve of a
//! model that has any of them linearizes a copy of the model through
//! `rusolve_cli::gc_solve` (the route `rusolve mip` takes for the same objects in a file), solves
//! that copy, checks the answer against the objects and cuts it back to the user's columns and
//! rows. The model's own solver state is never changed, and a model without these objects never
//! reaches this module (`Model::solve` tests `GcModel::has_general`).

use std::time::Instant;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use rusolve_ekk::lp_data::hconst::HighsVarType;
use rusolve_ekk::lp_data::highs_lp::HighsLp;
use rusolve_ekk::lp_data::highs_run::Highs as Port;
use rusolve_io::{
    FuncCon, FuncKind, GcModel, GcTol, Indicator, Pwl, RowSense, SemiVar, SosKind, SosSet,
    VariableBound,
};

use super::{INF, Model, port_err};
use crate::diag::{gc_err, gc_unsupported};
use crate::expr::{Constraint, MIXED_MODELS};
use crate::ingest::{self, Cols, F64s};
use crate::solution::Solution;

/// Linearization limit on a derived big-M, as `rusolve mip --gc-bigm-max`.
const BIGM_MAX: f64 = 1e6;

impl Model {
    /// `(lower, upper, integer)` of column `j` as the model has it now.
    pub(super) fn gc_col_info(&self, j: usize) -> (f64, f64, bool) {
        if let Some(info) = self.b.unsynced_col_info(j) {
            return info;
        }
        let lp = &self
            .port
            .as_ref()
            .expect("a synced column has a solver")
            .model_lp_;
        (
            lp.col_lower_[j],
            lp.col_upper_[j],
            lp.integrality_
                .get(j)
                .is_some_and(|t| *t == HighsVarType::kInteger),
        )
    }

    /// The column's name, or `c<j>` (what `check` and the envelope call an unnamed column).
    pub(super) fn gc_col_label(&self, j: usize) -> String {
        let named = if j < self.synced_cols() {
            self.port
                .as_ref()
                .and_then(|p| p.model_lp_.col_names_.get(j).cloned())
                .filter(|s| !s.is_empty())
        } else {
            self.b.col_name(j).map(str::to_string)
        };
        named.unwrap_or_else(|| format!("c{j}"))
    }

    fn gc_require_binary(&self, what: &str, j: usize) -> PyResult<()> {
        let (lo, hi, integer) = self.gc_col_info(j);
        if integer && lo >= 0.0 && hi <= 1.0 {
            return Ok(());
        }
        Err(gc_err(format!(
            "{what}: variable '{}' is not binary (it needs integer type and bounds within [0, 1]; \
             declare it with `add_var(binary=True)`)",
            self.gc_col_label(j)
        )))
    }

    /// The object name: `requested` if it is free, else `<prefix><k>` for the first free `k`.
    fn gc_name(&self, requested: Option<String>, prefix: &str, first: usize) -> PyResult<String> {
        let g = &self.gc;
        let used = |n: &str| {
            g.indicators.iter().any(|o| o.name == n)
                || g.sos.iter().any(|o| o.name == n)
                || g.pwl.iter().any(|o| o.name == n)
                || g.funcs.iter().any(|o| o.name == n)
        };
        match requested {
            Some(n) => {
                if n.is_empty() || n.chars().any(|c| c.is_whitespace() || c == ':') {
                    return Err(gc_err(format!(
                        "name {n:?}: a general-constraint name is not empty and has no space or ':' (it is written to LP and MPS files)"
                    )));
                }
                if used(&n) {
                    return Err(gc_err(format!(
                        "name {n:?} is already used by another general constraint"
                    )));
                }
                Ok(n)
            }
            None => {
                let mut k = first;
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

    pub(super) fn gc_cols(&self, o: &Bound<'_, PyAny>, what: &str) -> PyResult<Vec<usize>> {
        let cols = Cols::read(o, self.id, self.b.num_cols(), what)?;
        Ok((0..cols.len()).map(|k| cols.get(k)).collect())
    }

    fn gc_floats(o: &Bound<'_, PyAny>, what: &str) -> PyResult<Vec<f64>> {
        Ok(F64s::read_array(o, what)?.slice().to_vec())
    }

    // ------------------------------------------------------------ the add_* methods

    pub(super) fn gc_indicator(
        &mut self,
        var: &Bound<'_, PyAny>,
        value: &Bound<'_, PyAny>,
        constraint: &Bound<'_, PyAny>,
        name: Option<String>,
    ) -> PyResult<String> {
        let name = self.gc_name(name, "ind", self.gc.indicators.len())?;
        let what = format!("indicator '{name}'");
        let z = self.col_index(var)?;
        self.gc_require_binary(&what, z)?;
        let bin_val = match ingest::number(value)? {
            Some(v) if v == 0.0 || v == 1.0 => v == 1.0,
            _ => {
                return Err(gc_err(format!(
                    "{what}: the value must be 0 or 1 (the value of the binary that switches the row on)"
                )));
            }
        };
        let Ok(c) = constraint.cast_exact::<Constraint>() else {
            return Err(PyTypeError::new_err(format!(
                "{what}: the row must be a Constraint such as `x + y <= 5` (not added to the model), got {}",
                ingest::type_name(constraint)
            )));
        };
        let c = c.get();
        if c.model != 0 && c.model != self.id {
            return Err(crate::diag::model_err(format!(
                "{MIXED_MODELS}: the indicator's Constraint belongs to another Model"
            )));
        }
        let (sense, rhs) = if c.lo == c.hi {
            (RowSense::Equal, c.lo)
        } else if c.lo.is_finite() && c.hi == INF {
            (RowSense::GreaterEqual, c.lo)
        } else if c.hi.is_finite() && c.lo == -INF {
            (RowSense::LessEqual, c.hi)
        } else {
            return Err(gc_err(format!(
                "{what}: the row must be one-sided (<= or >=) or an equality, not a range or a free row; \
                 use two indicators for a range"
            )));
        };
        let terms: Vec<(usize, f64)> = c
            .body
            .with(|t| t.iter().map(|&(j, a)| (j as usize, a)).collect());
        self.gc.indicators.push(Indicator {
            name: name.clone(),
            bin: z,
            bin_val,
            terms,
            sense,
            rhs,
        });
        Ok(name)
    }

    pub(super) fn gc_sos(
        &mut self,
        kind: SosKind,
        vars: &Bound<'_, PyAny>,
        weights: Option<&Bound<'_, PyAny>>,
        name: Option<String>,
    ) -> PyResult<String> {
        let prefix = if kind == SosKind::Sos1 {
            "sos1_"
        } else {
            "sos2_"
        };
        let name = self.gc_name(name, prefix, self.gc.sos.len())?;
        let what = format!("sos '{name}'");
        let cols = self.gc_cols(vars, "vars")?;
        if cols.is_empty() {
            return Err(gc_err(format!("{what} has no members")));
        }
        let w = match weights.filter(|w| !w.is_none()) {
            None => (1..=cols.len()).map(|k| k as f64).collect(),
            Some(w) => Self::gc_floats(w, "weights")?,
        };
        if w.len() != cols.len() {
            return Err(gc_err(format!(
                "{what}: {} variables but {} weights (one weight per variable, in the order the set follows)",
                cols.len(),
                w.len()
            )));
        }
        if let Some(bad) = w.iter().find(|v| !v.is_finite()) {
            return Err(gc_err(format!("{what}: weight {bad} is not finite")));
        }
        let mut sorted = cols.clone();
        sorted.sort_unstable();
        if let Some(d) = sorted.windows(2).find(|p| p[0] == p[1]) {
            return Err(gc_err(format!(
                "{what}: variable '{}' appears twice",
                self.gc_col_label(d[0])
            )));
        }
        let mut ws = w.clone();
        ws.sort_by(f64::total_cmp);
        if let Some(d) = ws.windows(2).find(|p| p[0] == p[1]) {
            return Err(gc_err(format!(
                "{what}: weight {} appears twice (weights give the order of the set and must differ)",
                d[0]
            )));
        }
        self.gc.sos.push(SosSet {
            name: name.clone(),
            kind,
            priority: 1.0,
            members: cols.into_iter().zip(w).collect(),
        });
        Ok(name)
    }

    pub(super) fn gc_pwl(
        &mut self,
        x: &Bound<'_, PyAny>,
        y: &Bound<'_, PyAny>,
        xs: &Bound<'_, PyAny>,
        ys: &Bound<'_, PyAny>,
        name: Option<String>,
    ) -> PyResult<String> {
        let name = self.gc_name(name, "pwl", self.gc.pwl.len())?;
        let what = format!("pwl '{name}'");
        let (cx, cy) = (self.col_index(x)?, self.col_index(y)?);
        if cx == cy {
            return Err(gc_err(format!("{what}: x and y are the same variable")));
        }
        let (xs, ys) = (Self::gc_floats(xs, "xs")?, Self::gc_floats(ys, "ys")?);
        if xs.len() != ys.len() {
            return Err(gc_err(format!(
                "{what}: {} x breakpoints but {} y values (one y per breakpoint)",
                xs.len(),
                ys.len()
            )));
        }
        if xs.len() < 2 {
            return Err(gc_err(format!(
                "{what} needs at least two breakpoints, got {}",
                xs.len()
            )));
        }
        if xs.iter().chain(&ys).any(|v| !v.is_finite()) {
            return Err(gc_err(format!("{what}: breakpoints must be finite")));
        }
        if let Some(k) = xs.windows(2).position(|w| w[1] < w[0]) {
            return Err(gc_err(format!(
                "{what}: x breakpoints must be sorted in non-decreasing order (x[{}] = {} follows x[{k}] = {})",
                k + 1,
                xs[k + 1],
                xs[k]
            )));
        }
        self.gc.pwl.push(Pwl {
            name: name.clone(),
            x: cx,
            y: cy,
            xs,
            ys,
        });
        Ok(name)
    }

    pub(super) fn gc_func(
        &mut self,
        kind: FuncKind,
        res: &Bound<'_, PyAny>,
        args: &Bound<'_, PyAny>,
        constant: Option<f64>,
        name: Option<String>,
    ) -> PyResult<String> {
        let name = self.gc_name(name, kind.word(), self.gc.funcs.len())?;
        let what = format!("{} '{name}'", kind.word());
        let r = self.col_index(res)?;
        let args = self.gc_cols(args, "args")?;
        let single = matches!(kind, FuncKind::Abs | FuncKind::Not);
        if args.is_empty() || (single && args.len() != 1) {
            return Err(gc_err(format!(
                "{what}: {} takes {} operand(s), got {}",
                kind.word(),
                if single {
                    "exactly one"
                } else {
                    "at least one"
                },
                args.len()
            )));
        }
        if let Some(c) = constant
            && !c.is_finite()
        {
            return Err(gc_err(format!("{what}: the constant {c} is not finite")));
        }
        if matches!(kind, FuncKind::And | FuncKind::Or | FuncKind::Not) {
            for &j in args.iter().chain(std::iter::once(&r)) {
                self.gc_require_binary(&what, j)?;
            }
        }
        self.gc.funcs.push(FuncCon {
            name: name.clone(),
            kind,
            res: r,
            args,
            constant,
        });
        Ok(name)
    }

    // ------------------------------------------------------------ solve, check, write

    /// The model's objects plus its semi-continuous columns (which live in the solver's
    /// integrality), the form the linearization and the checker take.
    pub(super) fn gc_with_semi(&self, lp: &HighsLp) -> GcModel {
        let mut g = self.gc.clone();
        for (col, t) in lp.integrality_.iter().enumerate() {
            match t {
                HighsVarType::kSemiContinuous => g.semi.push(SemiVar {
                    col,
                    integer: false,
                }),
                HighsVarType::kSemiInteger => g.semi.push(SemiVar { col, integer: true }),
                _ => {}
            }
        }
        g
    }

    /// The structural checks of the file readers, on the model as it is now: a later
    /// `set_bounds` or `set_integer` can have made an object invalid since it was added.
    fn gc_validate(&self, lp: &HighsLp, gc: &GcModel) -> PyResult<()> {
        let bounds: Vec<VariableBound> = (0..lp.num_col_ as usize)
            .map(|j| VariableBound {
                lower: lp.col_lower_[j],
                upper: lp.col_upper_[j],
                integer: lp
                    .integrality_
                    .get(j)
                    .is_some_and(|t| *t == HighsVarType::kInteger),
            })
            .collect();
        gc.validate(&bounds, &lp.col_names_).map_err(gc_err)
    }

    /// RS4010 when the model has general constraints and `what` would ignore them.
    pub(super) fn refuse_gc(&self, what: &str) -> PyResult<()> {
        match self.gc.first_general() {
            Some(o) => Err(gc_unsupported(what, &o)),
            None => Ok(()),
        }
    }

    /// Solve the linearization of the model on a separate solver, check the answer against the
    /// user's objects (rounded binaries), repair trickle flow once, and report on the user's
    /// columns and rows. Options are those `solve()` just applied.
    pub(super) fn solve_gc(&mut self, py: Python<'_>) -> PyResult<Solution> {
        let t0 = Instant::now();
        let port = self.port.as_ref().expect("synced");
        let user = port.model_lp_.clone();
        let opts = port.get_options().clone();
        let gc = self.gc_with_semi(&user);
        self.gc_validate(&user, &gc)?;
        let mut lp = user;
        let plan = rusolve_cli::gc_solve::prepare(&mut lp, gc, BIGM_MAX)
            .map_err(|e| gc_err(format!("cannot linearize {e}")))?;
        let tol = GcTol {
            integrality: opts.mip_feasibility_tolerance,
            feasibility: opts.mip_feasibility_tolerance,
        };
        let mut h = Port::with_options(opts);
        port_err("pass_model", h.pass_model(lp))?;
        let run = py.detach(|| {
            h.reset_scheduler_for_threads();
            h.run()
        });
        let out = rusolve_cli::gc_solve::finish(&plan, &mut h, tol);
        let time = t0.elapsed().as_secs_f64();
        Ok(self.solution_from(&h, Some(&plan.user_lp), Some(out.report), None, run, time))
    }
}
