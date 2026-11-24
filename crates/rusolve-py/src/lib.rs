//! Python bindings, extension module `rusolve._rusolve` (re-exported by `python/rusolve`).
//!
//! Two APIs live here:
//! - `rusolve.Model` (docs/design/python-api.md): the fast native API, in `expr.rs`,
//!   `ingest.rs`, `model.rs` and `solution.rs` over `rusolve_model::ModelBuilder`.
//! - `rusolve.Highs` (this file): the highspy-shaped compatibility class.
//!
//! The rest of this header describes `Highs`.
//!
//! The API mirrors `highspy` (HiGHS 1.16.0, `highs/highs_bindings/highs_bindings.cpp` and
//! `highs/highspy/highs.py`): method names, the `HighsLp` field names with trailing
//! underscores, the `HighsModelStatus`/`HighsStatus`/`ObjSense`/`MatrixFormat`/`HighsVarType`
//! enums with their `k...` names, and `setOptionValue`/`getOptionValue` by option name.
//! Everything is a thin layer over the ported `rusolve_ekk::lp_data::highs_run::Highs`.
//!
//! Defaults equal the CLI's: `Highs()` starts from `HighsOptions::recommended()` (the
//! profile `crates/rusolve-cli/src/main.rs::base_options` applies), and
//! `Highs(profile="highs-1.15.1")` from `HighsOptions::default()`.
//!
//! Not in highspy and not here: callbacks, `HighsModel` (Hessians), basis access, IIS.

mod callback;
mod diag;
mod expr;
mod ingest;
mod model;
mod quad;
mod solution;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyFloat, PyInt, PyString};
use rusolve_ekk::HighsInt;
use rusolve_ekk::lp_data::hconst::kHighsInf;
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_run::Highs as Port;
use rusolve_ekk::lp_data::model_utils::util_model_status_to_string;
use rusolve_ekk::util::c_printf::c_fmt_g;

// ---------------------------------------------------------------- enums

macro_rules! py_enum {
    ($(#[$m:meta])* $name:ident <-> $rust:ident { $($v:ident = $n:expr),+ $(,)? }) => {
        $(#[$m])*
        #[pyclass(eq, eq_int, from_py_object, module = "rusolve")]
        #[derive(Clone, Copy, PartialEq, Eq, Debug)]
        #[allow(non_camel_case_types)]
        pub enum $name { $($v = $n),+ }
        impl From<$rust> for $name {
            fn from(v: $rust) -> Self { match v { $($rust::$v => $name::$v),+ } }
        }
        impl From<$name> for $rust {
            fn from(v: $name) -> Self { match v { $($name::$v => $rust::$v),+ } }
        }
    };
}

py_enum!(
    /// highspy `HighsStatus`.
    HighsStatus <-> RStatus { kError = -1, kOk = 0, kWarning = 1 }
);
py_enum!(
    /// highspy `HighsModelStatus` (HConst.h).
    HighsModelStatus <-> RModelStatus {
        kNotset = 0, kLoadError = 1, kModelError = 2, kPresolveError = 3, kSolveError = 4,
        kPostsolveError = 5, kModelEmpty = 6, kOptimal = 7, kInfeasible = 8,
        kUnboundedOrInfeasible = 9, kUnbounded = 10, kObjectiveBound = 11,
        kObjectiveTarget = 12, kTimeLimit = 13, kIterationLimit = 14, kUnknown = 15,
        kSolutionLimit = 16, kInterrupt = 17, kMemoryLimit = 18, kHighsInterrupt = 19,
    }
);
py_enum!(
    /// highspy `ObjSense`.
    ObjSense <-> RObjSense { kMinimize = 1, kMaximize = -1 }
);
py_enum!(
    /// highspy `MatrixFormat`.
    MatrixFormat <-> RMatrixFormat { kColwise = 1, kRowwise = 2, kRowwisePartitioned = 3 }
);
py_enum!(
    /// highspy `HighsVarType`.
    HighsVarType <-> RVarType {
        kContinuous = 0, kInteger = 1, kSemiContinuous = 2, kSemiInteger = 3, kImplicitInteger = 4,
    }
);

use rusolve_ekk::lp_data::hconst::HighsModelStatus as RModelStatus;
use rusolve_ekk::lp_data::hconst::HighsVarType as RVarType;
use rusolve_ekk::lp_data::hconst::MatrixFormat as RMatrixFormat;
use rusolve_ekk::lp_data::hconst::ObjSense as RObjSense;
use rusolve_ekk::lp_data::highs_status::HighsStatus as RStatus;

// ---------------------------------------------------------------- model objects

/// highspy `HighsSparseMatrix`.
#[pyclass(module = "rusolve", skip_from_py_object)]
#[derive(Clone)]
pub struct HighsSparseMatrix {
    #[pyo3(get, set)]
    pub format_: MatrixFormat,
    #[pyo3(get, set)]
    pub num_col_: i32,
    #[pyo3(get, set)]
    pub num_row_: i32,
    #[pyo3(get, set)]
    pub start_: Vec<i32>,
    #[pyo3(get, set)]
    pub index_: Vec<i32>,
    #[pyo3(get, set)]
    pub value_: Vec<f64>,
}

#[pymethods]
impl HighsSparseMatrix {
    #[new]
    fn new() -> Self {
        Self {
            format_: MatrixFormat::kColwise,
            num_col_: 0,
            num_row_: 0,
            start_: vec![],
            index_: vec![],
            value_: vec![],
        }
    }
}

/// highspy `HighsLp`. `a_matrix_` is a shared `HighsSparseMatrix`, so
/// `lp.a_matrix_.start_ = [...]` works as in highspy.
#[pyclass(module = "rusolve")]
pub struct HighsLp {
    #[pyo3(get, set)]
    pub num_col_: i32,
    #[pyo3(get, set)]
    pub num_row_: i32,
    #[pyo3(get, set)]
    pub col_cost_: Vec<f64>,
    #[pyo3(get, set)]
    pub col_lower_: Vec<f64>,
    #[pyo3(get, set)]
    pub col_upper_: Vec<f64>,
    #[pyo3(get, set)]
    pub row_lower_: Vec<f64>,
    #[pyo3(get, set)]
    pub row_upper_: Vec<f64>,
    #[pyo3(get, set)]
    pub a_matrix_: Py<HighsSparseMatrix>,
    #[pyo3(get, set)]
    pub sense_: ObjSense,
    #[pyo3(get, set)]
    pub offset_: f64,
    #[pyo3(get, set)]
    pub model_name_: String,
    #[pyo3(get, set)]
    pub col_names_: Vec<String>,
    #[pyo3(get, set)]
    pub row_names_: Vec<String>,
    #[pyo3(get, set)]
    pub integrality_: Vec<HighsVarType>,
}

#[pymethods]
impl HighsLp {
    #[new]
    fn new(py: Python<'_>) -> PyResult<Self> {
        Ok(Self {
            num_col_: 0,
            num_row_: 0,
            col_cost_: vec![],
            col_lower_: vec![],
            col_upper_: vec![],
            row_lower_: vec![],
            row_upper_: vec![],
            a_matrix_: Py::new(py, HighsSparseMatrix::new())?,
            sense_: ObjSense::kMinimize,
            offset_: 0.0,
            model_name_: String::new(),
            col_names_: vec![],
            row_names_: vec![],
            integrality_: vec![],
        })
    }
}

impl HighsLp {
    /// The port's `HighsLp` for this object.
    fn to_port(&self, py: Python<'_>) -> PyResult<rusolve_ekk::lp_data::highs_lp::HighsLp> {
        let m = self.a_matrix_.borrow(py);
        let mut lp = rusolve_ekk::lp_data::highs_lp::HighsLp {
            num_col_: self.num_col_,
            num_row_: self.num_row_,
            col_cost_: self.col_cost_.clone(),
            col_lower_: self.col_lower_.clone(),
            col_upper_: self.col_upper_.clone(),
            row_lower_: self.row_lower_.clone(),
            row_upper_: self.row_upper_.clone(),
            sense_: self.sense_.into(),
            offset_: self.offset_,
            model_name_: self.model_name_.clone(),
            col_names_: self.col_names_.clone(),
            row_names_: self.row_names_.clone(),
            integrality_: self.integrality_.iter().map(|&t| t.into()).collect(),
            ..Default::default()
        };
        lp.a_matrix_.format_ = m.format_.into();
        lp.a_matrix_.num_col_ = self.num_col_;
        lp.a_matrix_.num_row_ = self.num_row_;
        lp.a_matrix_.start_ = m.start_.clone();
        lp.a_matrix_.index_ = m.index_.clone();
        lp.a_matrix_.value_ = m.value_.clone();
        Ok(lp)
    }
}

/// highspy `HighsInfo` (the subset the port fills that scripts read).
#[pyclass(module = "rusolve", get_all, skip_from_py_object)]
#[derive(Clone)]
pub struct HighsInfo {
    pub valid: bool,
    pub objective_function_value: f64,
    pub mip_dual_bound: f64,
    pub mip_gap: f64,
    pub mip_node_count: i64,
    pub simplex_iteration_count: i32,
    pub ipm_iteration_count: i32,
}

/// highspy `HighsSolution`.
#[pyclass(module = "rusolve", get_all, set_all, skip_from_py_object)]
#[derive(Clone)]
pub struct HighsSolution {
    pub value_valid: bool,
    pub dual_valid: bool,
    pub col_value: Vec<f64>,
    pub col_dual: Vec<f64>,
    pub row_value: Vec<f64>,
    pub row_dual: Vec<f64>,
}

// ---------------------------------------------------------------- Highs

const DOUBLE_OPTIONS: &[&str] = &[
    "dual_feasibility_tolerance",
    "objective_bound",
    "primal_feasibility_tolerance",
    "primal_simplex_bound_perturbation_multiplier",
    "time_limit",
    "qp_regularization_value",
    "iis_time_limit",
    "icrash_starting_weight",
    "pdlp_optimality_tolerance",
    "mip_rel_gap",
    "mip_abs_gap",
];

fn status(s: RStatus) -> HighsStatus {
    s.into()
}

/// `setOptionValue`'s dispatch by Python type (bool, int, float, str), as highspy; also used
/// by `Model.solve(options=...)`. An int given for a double option is accepted.
pub(crate) fn set_option_any(
    h: &mut Port,
    name: &str,
    value: &Bound<'_, PyAny>,
) -> PyResult<RStatus> {
    Ok(if value.is_instance_of::<PyBool>() {
        h.set_option_value_bool(name, value.extract()?)
    } else if value.is_instance_of::<PyInt>() {
        // an int given for a double option (e.g. time_limit=10) is accepted, as in highspy
        if DOUBLE_OPTIONS.contains(&name) {
            h.set_option_value_double(name, value.extract::<f64>()?)
        } else {
            h.set_option_value_int(name, value.extract()?)
        }
    } else if value.is_instance_of::<PyFloat>() {
        h.set_option_value_double(name, value.extract()?)
    } else if value.is_instance_of::<PyString>() {
        h.set_option_value_str(name, &value.extract::<String>()?)
    } else {
        return Err(PyValueError::new_err(
            "option value must be bool, int, float or str",
        ));
    })
}

/// highspy `Highs`, over the ported `Highs::run`.
#[pyclass(unsendable, module = "rusolve")]
pub struct Highs {
    inner: Port,
}

#[pymethods]
impl Highs {
    /// `profile`: `"recommended"` (the CLI default) or `"highs-1.15.1"`
    /// (crates/rusolve-cli/src/main.rs::base_options).
    #[new]
    #[pyo3(signature = (profile="recommended"))]
    fn new(profile: &str) -> PyResult<Self> {
        let mut options = match profile {
            "recommended" => HighsOptions::recommended(),
            "highs-1.15.1" => HighsOptions::default(),
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown profile {other:?}: expected \"recommended\" or \"highs-1.15.1\""
                )));
            }
        };
        options.sync_log_options();
        Ok(Self {
            inner: Port::with_options(options),
        })
    }

    /// highs.py `version()`.
    fn version(&self) -> String {
        env!("CARGO_PKG_VERSION").to_string()
    }

    /// Read an MPS file (same reader and conversion as `rusolve mip`, with names).
    #[pyo3(name = "readModel")]
    fn read_model(&mut self, path: &str) -> HighsStatus {
        match read_mip(path) {
            Ok(lp) => status(self.inner.pass_model(lp)),
            Err(e) => {
                eprintln!("{e}");
                HighsStatus::kError
            }
        }
    }

    /// `passModel(lp)` or highspy's flat overload
    /// `passModel(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower,
    /// col_upper, row_lower, row_upper, a_start, a_index, a_value[, integrality])`.
    #[pyo3(name = "passModel", signature = (*args))]
    fn pass_model(
        &mut self,
        py: Python<'_>,
        args: &Bound<'_, pyo3::types::PyTuple>,
    ) -> PyResult<HighsStatus> {
        let lp = if args.len() == 1 {
            let lp: PyRef<'_, HighsLp> = args.get_item(0)?.extract()?;
            lp.to_port(py)?
        } else if args.len() == 14 || args.len() == 15 {
            let mut lp = rusolve_ekk::lp_data::highs_lp::HighsLp {
                num_col_: args.get_item(0)?.extract()?,
                num_row_: args.get_item(1)?.extract()?,
                sense_: args.get_item(4)?.extract::<ObjSense>()?.into(),
                offset_: args.get_item(5)?.extract()?,
                col_cost_: args.get_item(6)?.extract()?,
                col_lower_: args.get_item(7)?.extract()?,
                col_upper_: args.get_item(8)?.extract()?,
                row_lower_: args.get_item(9)?.extract()?,
                row_upper_: args.get_item(10)?.extract()?,
                ..Default::default()
            };
            lp.a_matrix_.format_ = args.get_item(3)?.extract::<MatrixFormat>()?.into();
            lp.a_matrix_.num_col_ = lp.num_col_;
            lp.a_matrix_.num_row_ = lp.num_row_;
            lp.a_matrix_.start_ = args.get_item(11)?.extract()?;
            lp.a_matrix_.index_ = args.get_item(12)?.extract()?;
            lp.a_matrix_.value_ = args.get_item(13)?.extract()?;
            if args.len() == 15 {
                let it: Vec<HighsVarType> = args.get_item(14)?.extract()?;
                lp.integrality_ = it.into_iter().map(Into::into).collect();
            }
            lp
        } else {
            return Err(PyValueError::new_err(
                "passModel takes a HighsLp, or the 14/15 flat arguments of highspy",
            ));
        };
        Ok(status(self.inner.pass_model(lp)))
    }

    /// `setOptionValue(name, value)`: bool, int, float or str by Python type, as highspy.
    #[pyo3(name = "setOptionValue")]
    fn set_option_value(&mut self, name: &str, value: &Bound<'_, PyAny>) -> PyResult<HighsStatus> {
        Ok(status(set_option_any(&mut self.inner, name, value)?))
    }

    /// `getOptionValue(name)`: returns `(HighsStatus, value)` as highspy does.
    #[pyo3(name = "getOptionValue")]
    fn get_option_value(
        &mut self,
        py: Python<'_>,
        name: &str,
    ) -> PyResult<(HighsStatus, Py<PyAny>)> {
        let o = &self.inner.options_;
        if DOUBLE_OPTIONS.contains(&name) {
            let v = match name {
                "dual_feasibility_tolerance" => o.dual_feasibility_tolerance,
                "objective_bound" => o.objective_bound,
                "primal_feasibility_tolerance" => o.primal_feasibility_tolerance,
                "primal_simplex_bound_perturbation_multiplier" => {
                    o.primal_simplex_bound_perturbation_multiplier
                }
                "time_limit" => o.time_limit,
                "qp_regularization_value" => o.qp_regularization_value,
                "iis_time_limit" => o.iis_time_limit,
                "icrash_starting_weight" => o.icrash_starting_weight,
                "pdlp_optimality_tolerance" => o.pdlp_optimality_tolerance,
                "mip_rel_gap" => o.mip_rel_gap,
                _ => o.mip_abs_gap,
            };
            return Ok((HighsStatus::kOk, v.into_pyobject(py)?.into_any().unbind()));
        }
        // The typed getters log an error on a kind mismatch; probe quietly.
        let (of, lc) = (
            self.inner.options_.output_flag,
            self.inner.options_.log_to_console,
        );
        self.inner.options_.output_flag = false;
        self.inner.options_.log_to_console = false;
        self.inner.options_.sync_log_options();
        let found: Option<Py<PyAny>> = if let Some(v) = self.inner.get_option_value_bool(name) {
            Some(PyBool::new(py, v).to_owned().into_any().unbind())
        } else if let Some(v) = self.inner.get_option_value_int(name) {
            Some(v.into_pyobject(py)?.into_any().unbind())
        } else {
            self.inner
                .get_option_value_str(name)
                .map(|v| v.into_pyobject(py).map(|o| o.into_any().unbind()))
                .transpose()?
        };
        self.inner.options_.output_flag = of;
        self.inner.options_.log_to_console = lc;
        self.inner.options_.sync_log_options();
        Ok(match found {
            Some(v) => (HighsStatus::kOk, v),
            None => (HighsStatus::kError, py.None()),
        })
    }

    fn run(&mut self) -> HighsStatus {
        status(self.inner.run())
    }

    #[pyo3(name = "getModelStatus")]
    fn get_model_status(&self) -> HighsModelStatus {
        self.inner.get_model_status().into()
    }

    #[pyo3(name = "modelStatusToString")]
    fn model_status_to_string(&self, s: HighsModelStatus) -> String {
        util_model_status_to_string(s.into()).to_string()
    }

    #[pyo3(name = "getInfo")]
    fn get_info(&self) -> HighsInfo {
        let i = self.inner.get_info();
        HighsInfo {
            valid: i.valid,
            objective_function_value: i.objective_function_value,
            mip_dual_bound: i.mip_dual_bound,
            mip_gap: i.mip_gap,
            mip_node_count: i.mip_node_count,
            simplex_iteration_count: i.simplex_iteration_count,
            ipm_iteration_count: i.ipm_iteration_count,
        }
    }

    #[pyo3(name = "getObjectiveValue")]
    fn get_objective_value(&self) -> f64 {
        self.inner.get_info().objective_function_value
    }

    #[pyo3(name = "getSolution")]
    fn get_solution(&self) -> HighsSolution {
        let s = self.inner.get_solution();
        HighsSolution {
            value_valid: s.value_valid,
            dual_valid: s.dual_valid,
            col_value: s.col_value.clone(),
            col_dual: s.col_dual.clone(),
            row_value: s.row_value.clone(),
            row_dual: s.row_dual.clone(),
        }
    }

    #[pyo3(name = "getNumCol")]
    fn get_num_col(&self) -> i32 {
        self.inner.get_num_col()
    }

    #[pyo3(name = "getNumRow")]
    fn get_num_row(&self) -> i32 {
        self.inner.get_num_row()
    }

    /// A `HighsLp` copy of the current model.
    #[pyo3(name = "getLp")]
    fn get_lp(&self, py: Python<'_>) -> PyResult<HighsLp> {
        let p = self.inner.get_lp();
        let mut lp = HighsLp::new(py)?;
        lp.num_col_ = p.num_col_;
        lp.num_row_ = p.num_row_;
        lp.col_cost_ = p.col_cost_.clone();
        lp.col_lower_ = p.col_lower_.clone();
        lp.col_upper_ = p.col_upper_.clone();
        lp.row_lower_ = p.row_lower_.clone();
        lp.row_upper_ = p.row_upper_.clone();
        lp.sense_ = p.sense_.into();
        lp.offset_ = p.offset_;
        lp.model_name_ = p.model_name_.clone();
        lp.col_names_ = p.col_names_.clone();
        lp.row_names_ = p.row_names_.clone();
        lp.integrality_ = p.integrality_.iter().map(|&t| t.into()).collect();
        {
            let mut m = lp.a_matrix_.borrow_mut(py);
            m.format_ = p.a_matrix_.format_.into();
            m.num_col_ = p.num_col_;
            m.num_row_ = p.num_row_;
            m.start_ = p.a_matrix_.start_.clone();
            m.index_ = p.a_matrix_.index_.clone();
            m.value_ = p.a_matrix_.value_.clone();
        }
        Ok(lp)
    }

    #[pyo3(name = "changeColBounds")]
    fn change_col_bounds(&mut self, col: i32, lower: f64, upper: f64) -> HighsStatus {
        status(self.inner.change_col_bounds(col, lower, upper))
    }

    #[pyo3(name = "changeRowBounds")]
    fn change_row_bounds(&mut self, row: i32, lower: f64, upper: f64) -> HighsStatus {
        status(self.inner.change_row_bounds(row, lower, upper))
    }

    #[pyo3(name = "changeColCost")]
    fn change_col_cost(&mut self, col: i32, cost: f64) -> HighsStatus {
        status(self.inner.change_col_cost(col, cost))
    }

    #[pyo3(name = "changeCoeff")]
    fn change_coeff(&mut self, row: i32, col: i32, value: f64) -> HighsStatus {
        status(self.inner.change_coeff(row, col, value))
    }

    #[pyo3(name = "changeObjectiveSense")]
    fn change_objective_sense(&mut self, sense: ObjSense) -> HighsStatus {
        status(self.inner.change_objective_sense(sense.into()))
    }

    /// `addRow(lower, upper, num_nz, indices, values)`.
    #[pyo3(name = "addRow")]
    fn add_row(
        &mut self,
        lower: f64,
        upper: f64,
        num_nz: i32,
        indices: Vec<i32>,
        values: Vec<f64>,
    ) -> PyResult<HighsStatus> {
        check_nz(num_nz, &indices, &values)?;
        Ok(status(
            self.inner.add_row(lower, upper, num_nz, &indices, &values),
        ))
    }

    /// `addCol(cost, lower, upper, num_nz, indices, values)`.
    #[pyo3(name = "addCol")]
    fn add_col(
        &mut self,
        cost: f64,
        lower: f64,
        upper: f64,
        num_nz: i32,
        indices: Vec<i32>,
        values: Vec<f64>,
    ) -> PyResult<HighsStatus> {
        check_nz(num_nz, &indices, &values)?;
        Ok(status(
            self.inner
                .add_col(cost, lower, upper, num_nz, &indices, &values),
        ))
    }

    /// Raw-style solution file after highspy `writeSolution(path, style=0)` (the format of
    /// HighsModelUtils.cpp::writeModelSolution, lines 268-356, without the primal/dual
    /// feasibility status line and the objective/gap lines; names are the model's, or
    /// `C<j>`/`R<i>` as HiGHS's default naming).
    #[pyo3(name = "writeSolution", signature = (path, style=0))]
    fn write_solution(&self, path: &str, style: i32) -> PyResult<HighsStatus> {
        if style != 0 {
            return Err(PyValueError::new_err(
                "only solution style 0 (raw) is implemented",
            ));
        }
        let lp = self.inner.get_lp();
        let sol = self.inner.get_solution();
        let name = |names: &[String], p: &str, i: usize| {
            names.get(i).cloned().unwrap_or_else(|| format!("{p}{i}"))
        };
        let g = |v: f64| c_fmt_g(v, 15);
        let mut out = format!(
            "Model status\n{}\n\n# Primal solution values\n",
            util_model_status_to_string(self.inner.get_model_status())
        );
        if sol.value_valid {
            out += &format!(
                "Objective {}\n# Columns {}\n",
                g(self.get_objective_value()),
                lp.num_col_
            );
            for (j, v) in sol.col_value.iter().enumerate() {
                out += &format!("{} {}\n", name(&lp.col_names_, "C", j), g(*v));
            }
            out += &format!("# Rows {}\n", lp.num_row_);
            for (i, v) in sol.row_value.iter().enumerate() {
                out += &format!("{} {}\n", name(&lp.row_names_, "R", i), g(*v));
            }
        } else {
            out += "None\n";
        }
        out += "\n# Dual solution values\n";
        if sol.dual_valid {
            out += &format!("# Columns {}\n", lp.num_col_);
            for (j, v) in sol.col_dual.iter().enumerate() {
                out += &format!("{} {}\n", name(&lp.col_names_, "C", j), g(*v));
            }
            out += &format!("# Rows {}\n", lp.num_row_);
            for (i, v) in sol.row_dual.iter().enumerate() {
                out += &format!("{} {}\n", name(&lp.row_names_, "R", i), g(*v));
            }
        } else {
            out += "None\n";
        }
        match std::fs::write(path, out) {
            Ok(()) => Ok(HighsStatus::kOk),
            Err(e) => {
                eprintln!("writeSolution {path}: {e}");
                Ok(HighsStatus::kError)
            }
        }
    }

    fn clear(&mut self) -> HighsStatus {
        self.inner.clear_model();
        HighsStatus::kOk
    }
}

fn check_nz(num_nz: i32, indices: &[i32], values: &[f64]) -> PyResult<()> {
    if num_nz < 0 || indices.len() < num_nz as usize || values.len() < num_nz as usize {
        return Err(PyValueError::new_err(
            "num_nz exceeds the length of indices/values",
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------- MPS reading

/// `HighsLp` with integrality and names, read in MPS file order. This is a copy of
/// crates/rusolve-cli/src/mip_highs.rs::read_mip (`pub(crate)` in the CLI, left there) plus
/// `crates/rusolve-cli/src/ekk.rs::highs_lp_from_model` (HMpsFF.cpp rules, see its doc),
/// with the MPS row/column names set on the model.
fn read_mip(path: &str) -> Result<rusolve_ekk::lp_data::highs_lp::HighsLp, String> {
    let (model, order) =
        rusolve_io::read_mps_file_file_order(path).map_err(|e| format!("read {path}: {e}"))?;
    let mut lp = highs_lp_from_model(&model, &order);
    // Highs::readModel (Highs.cpp:814) names the model after the file
    lp.model_name_ = extract_model_name(path);
    if model.bounds.iter().any(|b| b.integer) || !model.gc.semi.is_empty() {
        lp.integrality_ = model
            .bounds
            .iter()
            .map(|b| {
                if b.integer {
                    RVarType::kInteger
                } else {
                    RVarType::kContinuous
                }
            })
            .collect();
        // MPS `SC`/`SI` columns
        for s in &model.gc.semi {
            lp.integrality_[s.col] = if s.integer {
                RVarType::kSemiInteger
            } else {
                RVarType::kSemiContinuous
            };
        }
    }
    if let Some(what) = model.gc.first_general() {
        return Err(format!(
            "read {path}: {what} is not supported by the Python API yet"
        ));
    }
    if model.col_names.len() == model.num_variables {
        lp.col_names_ = model.col_names.clone();
    }
    if model.row_names.len() == model.num_constraints {
        lp.row_names_ = model.row_names.clone();
    }
    Ok(lp)
}

/// HiGHS: io/Filereader.cpp::extractModelName (lines 75-87); copy of the CLI's.
fn extract_model_name(filename: &str) -> String {
    let mut name = filename.to_string();
    if let Some(found) = name.rfind(['/', '\\']) {
        name = name[found + 1..].to_string();
    }
    let mut found = name.rfind('.');
    let ext = match found {
        Some(f) => &name[f + 1..],
        None => name.as_str(),
    };
    if ext == "gz" {
        let f = found.expect("\"gz\" extension has a dot or is the whole name");
        name.truncate(f);
        found = name.rfind('.');
    }
    if let Some(f) = found {
        name.truncate(f);
    }
    name
}

/// Copy of crates/rusolve-cli/src/ekk.rs::highs_lp_from_model (HMpsFF.cpp:136-170, :630-642,
/// :911-956, :1065-1078, :1554-1565): columns in file order, zeros skipped, first of
/// duplicate entries wins, RANGES applied to the row bounds.
fn highs_lp_from_model(
    m: &rusolve_io::Model,
    order: &rusolve_io::MpsFileOrder,
) -> rusolve_ekk::lp_data::highs_lp::HighsLp {
    use rusolve_io::{ObjSense as IoSense, RowSense};
    let n = m.num_variables;
    let nr = m.num_constraints;
    let mut col_cost = vec![0.0; n];
    let mut has_cost = vec![false; n];
    for &(j, v) in &order.cost_entries {
        if v != 0.0 && !has_cost[j] {
            col_cost[j] = v;
            has_cost[j] = true;
        }
    }
    let mut count = vec![0usize; n + 1];
    for &(_, j, v) in &order.entries {
        if v != 0.0 {
            count[j + 1] += 1;
        }
    }
    for j in 0..n {
        count[j + 1] += count[j];
    }
    let mut next = count.clone();
    let total = count[n];
    let mut index = vec![0 as HighsInt; total];
    let mut value = vec![0.0; total];
    for &(i, j, v) in &order.entries {
        if v != 0.0 {
            index[next[j]] = i as HighsInt;
            value[next[j]] = v;
            next[j] += 1;
        }
    }
    let mut start: Vec<HighsInt> = Vec::with_capacity(n + 1);
    let mut kept_index = Vec::with_capacity(total);
    let mut kept_value = Vec::with_capacity(total);
    let mut seen = vec![usize::MAX; nr];
    start.push(0);
    for j in 0..n {
        for k in count[j]..count[j + 1] {
            let i = index[k] as usize;
            if seen[i] != j {
                seen[i] = j;
                kept_index.push(index[k]);
                kept_value.push(value[k]);
            }
        }
        start.push(kept_index.len() as HighsInt);
    }
    let inf = kHighsInf;
    let mut row_lower = vec![0.0; nr];
    let mut row_upper = vec![0.0; nr];
    for i in 0..nr {
        let b = m.rhs[i];
        let (mut lo, mut up) = match m.senses[i] {
            RowSense::GreaterEqual => (b, inf),
            RowSense::Equal => (b, b),
            RowSense::LessEqual => (-inf, b),
        };
        let r = m.ranges[i];
        if r != 0.0 {
            match m.senses[i] {
                RowSense::LessEqual => lo = up - r.abs(),
                RowSense::GreaterEqual => up = lo + r.abs(),
                RowSense::Equal if r < 0.0 => lo = up - r.abs(),
                RowSense::Equal => up = lo + r.abs(),
            }
        }
        row_lower[i] = lo;
        row_upper[i] = up;
    }
    let mut lp = rusolve_ekk::lp_data::highs_lp::HighsLp {
        num_col_: n as HighsInt,
        num_row_: nr as HighsInt,
        col_cost_: col_cost,
        col_lower_: m.bounds.iter().map(|b| b.lower).collect(),
        col_upper_: m.bounds.iter().map(|b| b.upper).collect(),
        row_lower_: row_lower,
        row_upper_: row_upper,
        sense_: match m.obj_sense {
            IoSense::Minimize => RObjSense::kMinimize,
            IoSense::Maximize => RObjSense::kMaximize,
        },
        offset_: m.obj_offset,
        ..Default::default()
    };
    lp.a_matrix_.format_ = RMatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = n as HighsInt;
    lp.a_matrix_.num_row_ = nr as HighsInt;
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = kept_index;
    lp.a_matrix_.value_ = kept_value;
    lp
}

/// Tells the OR-Tools loader where the `ortools` Python package of this environment is (the
/// directory holding `.libs/`), because inside Python the executable is the interpreter and the
/// loader's paths relative to the executable do not apply (`rusolve_cpsat::runtime`). Loads nothing.
#[pyfunction]
fn _set_ortools_dir(dir: std::path::PathBuf) {
    rusolve_cpsat::runtime::set_hint(dir);
}

/// Whether the OR-Tools library for `mip_engine="cpsat"` can be loaded, without solving:
/// `{"available", "library", "version", "searched", "error"}`. Loads the library when it is found.
#[pyfunction]
fn cpsat_info(py: Python<'_>) -> PyResult<Bound<'_, pyo3::types::PyDict>> {
    use pyo3::types::PyDict;
    let a = rusolve_cpsat::runtime::availability();
    let d = PyDict::new(py);
    d.set_item("available", a.available)?;
    d.set_item("library", a.library.map(|p| p.display().to_string()))?;
    d.set_item("version", a.version)?;
    d.set_item("searched", a.searched)?;
    d.set_item("error", a.error)?;
    Ok(d)
}

#[pymodule]
fn _rusolve(m: &Bound<'_, PyModule>) -> PyResult<()> {
    expr::register(m)?;
    model::register(m)?;
    solution::register(m)?;
    m.add_class::<Highs>()?;
    m.add_class::<HighsLp>()?;
    m.add_class::<HighsSparseMatrix>()?;
    m.add_class::<HighsInfo>()?;
    m.add_class::<HighsSolution>()?;
    m.add_class::<HighsStatus>()?;
    m.add_class::<HighsModelStatus>()?;
    m.add_class::<ObjSense>()?;
    m.add_class::<MatrixFormat>()?;
    m.add_class::<HighsVarType>()?;
    m.add_function(wrap_pyfunction!(cpsat_info, m)?)?;
    m.add_function(wrap_pyfunction!(_set_ortools_dir, m)?)?;
    m.add("kHighsInf", f64::INFINITY)?;
    m.add("kHighsIInf", i32::MAX)?;
    Ok(())
}
