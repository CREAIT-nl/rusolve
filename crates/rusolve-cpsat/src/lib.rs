//! Google OR-Tools CP-SAT (Apache-2.0) as an optional engine.
//!
//! Two parts:
//! * [`translate`]: a model (`HighsLp`, whose columns must all be integer, plus the optional CP
//!   objects of `rusolve_io::cp`) to a `CpModelProto`, or a typed [`Refusal`] naming the first
//!   column or row CP-SAT cannot take exactly. No scaling, no rounding of coefficients or rows.
//! * [`solve`]: serialize the model and the `SatParameters`, call the OR-Tools C API
//!   (`ortools/sat/c_api/cp_solver_c.h`, `SolveCpModelWithParameters`), decode the
//!   `CpSolverResponse` into a [`CpOutcome`].
//!
//! The crate does not link OR-Tools: [`runtime`] opens the shared library the first time a solve
//! needs it (search order in that module's docs). The protobuf types in [`pb`] are generated from `proto/cp_model.proto` and
//! `proto/sat_parameters.proto`, copied from OR-Tools 9.15 (Apache-2.0, `proto/LICENSE-ORTOOLS`).
//!
//! Nothing here decides that a point is right: the caller must check the returned point against the
//! original model before reporting it (`rusolve-cli` does, with its own checker).

use prost::Message;
use rusolve_ekk::prelude::HighsLp;
use rusolve_io::cp::CpModel;
use std::ffi::{c_int, c_void};
use std::fmt;

pub mod runtime;

/// The protobuf types of `cp_model.proto` and `sat_parameters.proto` (package
/// `operations_research.sat`).
#[allow(clippy::all, clippy::pedantic, missing_docs)]
pub mod pb {
    include!(concat!(env!("OUT_DIR"), "/operations_research.sat.rs"));
}

mod translate;
pub use translate::{ObjectiveKind, Refusal, Translation, translate};

// SAFETY: `free` is the C library's, which the Rust standard library already links;
// `SolveCpModelWithParameters` (see [`runtime::SolveFn`]) stores a `malloc`ed buffer and its length
// in `*cres` / `*cres_len`; the caller owns that buffer and releases it with `free`. The C API
// aborts (absl CHECK) on a buffer that is not a valid message, so [`solve_proto`] only ever
// passes buffers it has just encoded itself.
unsafe extern "C" {
    fn free(p: *mut c_void);
}

/// Solver parameters, a subset of `SatParameters`.
#[derive(Clone, Debug)]
pub struct CpParams {
    /// Wall-clock limit in seconds (`max_time_in_seconds`); `None` is no limit.
    pub time_limit: Option<f64>,
    /// `num_workers`; 0 lets CP-SAT choose (one per core).
    pub workers: i32,
    pub seed: i32,
    /// `relative_gap_limit`.
    pub rel_gap: f64,
    /// `absolute_gap_limit`.
    pub abs_gap: f64,
    /// `log_search_progress` (CP-SAT writes its log to stdout).
    pub log: bool,
}

impl Default for CpParams {
    /// CP-SAT's own defaults for the gaps (relative 0, absolute 1e-4), one worker, seed 0.
    fn default() -> Self {
        Self {
            time_limit: None,
            workers: 1,
            seed: 0,
            rel_gap: 0.0,
            abs_gap: 1e-4,
            log: false,
        }
    }
}

/// What CP-SAT said, before any interpretation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CpStatus {
    /// `OPTIMAL`: proved optimal (within the gap limits), or, without an objective, feasible.
    Optimal,
    /// `FEASIBLE`: a solution, optimality not proved (a limit stopped the search).
    Feasible,
    /// `INFEASIBLE`: proved.
    Infeasible,
    /// `UNKNOWN`: no solution, nothing proved (a limit stopped the search).
    Unknown,
}

impl CpStatus {
    pub fn name(self) -> &'static str {
        match self {
            CpStatus::Optimal => "OPTIMAL",
            CpStatus::Feasible => "FEASIBLE",
            CpStatus::Infeasible => "INFEASIBLE",
            CpStatus::Unknown => "UNKNOWN",
        }
    }
}

/// The result of one CP-SAT solve.
#[derive(Clone, Debug)]
pub struct CpOutcome {
    pub status: CpStatus,
    /// The point, one value per column of the model, when CP-SAT returned a solution.
    pub x: Option<Vec<f64>>,
    /// CP-SAT's `objective_value` (the model's sense, offset included); `None` without a point.
    pub objective: Option<f64>,
    /// CP-SAT's `best_objective_bound` in the same terms (a lower bound when minimizing, an upper
    /// bound when maximizing); `None` when it is not finite.
    pub bound: Option<f64>,
    pub wall_time_s: f64,
    /// The `num_workers` value that ran.
    pub workers: i32,
    pub conflicts: i64,
    pub branches: i64,
    /// `integer` (the objective went in as `CpObjectiveProto`) or `floating_point`.
    pub objective_kind: ObjectiveKind,
    /// Variables added to the model beyond its columns (for the no-overlap size-0 rule).
    pub aux_variables: usize,
    /// `9.15.6755` or `unknown` ([`runtime::version_near`]).
    pub ortools_version: String,
    /// The library that ran, as [`runtime::Loaded::path`].
    pub ortools_library: std::path::PathBuf,
}

/// Why [`solve`] did not return an outcome.
#[derive(Clone, Debug, PartialEq)]
pub enum CpError {
    /// The model cannot be translated exactly.
    Refused(Refusal),
    /// CP-SAT's validator rejected the model; the message is CP-SAT's own.
    ModelInvalid(String),
    /// The serialized model is larger than the C API's `int` length can carry.
    TooLarge(usize),
    /// The OR-Tools library is not installed (or not loadable); the message lists the paths
    /// searched and how to install it.
    Unavailable(runtime::LoadError),
    /// A response that does not decode (a bug).
    BadResponse(String),
}

impl fmt::Display for CpError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CpError::Refused(r) => write!(f, "{r}"),
            CpError::ModelInvalid(m) => write!(f, "CP-SAT rejected the model: {m}"),
            CpError::TooLarge(n) => write!(
                f,
                "the CP-SAT model is {n} bytes, more than the C API accepts (2 GiB)"
            ),
            CpError::Unavailable(e) => write!(f, "{e}"),
            CpError::BadResponse(m) => {
                write!(f, "CP-SAT returned a response that does not decode: {m}")
            }
        }
    }
}

impl std::error::Error for CpError {}

impl From<Refusal> for CpError {
    fn from(r: Refusal) -> Self {
        CpError::Refused(r)
    }
}

/// Translates and solves. `lp` is the model's integer program (every column integer, finite
/// bounds, integral coefficients and right-hand sides); `cp` its CP objects, if any.
///
/// The returned point is CP-SAT's, unchecked: check it against the original model before use.
pub fn solve(lp: &HighsLp, cp: Option<&CpModel>, params: &CpParams) -> Result<CpOutcome, CpError> {
    let t = translate(lp, cp)?;
    solve_proto(&t, lp.num_col_ as usize, params)
}

/// [`solve`] for a model that is already translated; the first `num_cols` variables are the
/// model's columns.
pub fn solve_proto(
    t: &Translation,
    num_cols: usize,
    params: &CpParams,
) -> Result<CpOutcome, CpError> {
    let lib = runtime::library().map_err(|e| CpError::Unavailable(e.clone()))?;
    let workers = if params.workers > 0 {
        params.workers
    } else {
        std::thread::available_parallelism().map_or(1, |n| n.get() as i32)
    };
    let sat = pb::SatParameters {
        max_time_in_seconds: params.time_limit,
        num_workers: Some(workers),
        random_seed: Some(params.seed),
        relative_gap_limit: Some(params.rel_gap),
        absolute_gap_limit: Some(params.abs_gap),
        log_search_progress: Some(params.log),
        ..Default::default()
    };
    let req = t.model.encode_to_vec();
    let par = sat.encode_to_vec();
    let (Ok(req_len), Ok(par_len)) = (c_int::try_from(req.len()), c_int::try_from(par.len()))
    else {
        return Err(CpError::TooLarge(req.len()));
    };
    let mut res: *mut c_void = std::ptr::null_mut();
    let mut res_len: c_int = 0;
    // SAFETY: `lib.solve` is the symbol looked up in the loaded library, which stays mapped for
    // the life of the process. Both buffers are live for the call and `*_len` is their length;
    // the two out pointers are valid for writes. On return `res` holds `res_len` bytes allocated
    // with `malloc` by OR-Tools, copied out below and then released with `free`, once.
    let bytes = unsafe {
        (lib.solve)(
            req.as_ptr().cast(),
            req_len,
            par.as_ptr().cast(),
            par_len,
            &mut res,
            &mut res_len,
        );
        if res.is_null() || res_len < 0 {
            return Err(CpError::BadResponse("no response buffer".into()));
        }
        let copy = std::slice::from_raw_parts(res.cast::<u8>(), res_len as usize).to_vec();
        free(res);
        copy
    };
    let r = pb::CpSolverResponse::decode(bytes.as_slice())
        .map_err(|e| CpError::BadResponse(e.to_string()))?;
    let status = match pb::CpSolverStatus::try_from(r.status) {
        Ok(pb::CpSolverStatus::Optimal) => CpStatus::Optimal,
        Ok(pb::CpSolverStatus::Feasible) => CpStatus::Feasible,
        Ok(pb::CpSolverStatus::Infeasible) => CpStatus::Infeasible,
        Ok(pb::CpSolverStatus::Unknown) => CpStatus::Unknown,
        Ok(pb::CpSolverStatus::ModelInvalid) => {
            return Err(CpError::ModelInvalid(r.solution_info.clone()));
        }
        Err(_) => return Err(CpError::BadResponse(format!("status {}", r.status))),
    };
    let has_point =
        matches!(status, CpStatus::Optimal | CpStatus::Feasible) && r.solution.len() >= num_cols;
    let x = has_point.then(|| r.solution[..num_cols].iter().map(|&v| v as f64).collect());
    let finite = |v: f64| v.is_finite().then_some(v);
    Ok(CpOutcome {
        status,
        objective: if has_point {
            finite(r.objective_value)
        } else {
            None
        },
        bound: finite(r.best_objective_bound),
        x,
        wall_time_s: r.wall_time,
        workers,
        conflicts: r.num_conflicts,
        branches: r.num_branches,
        objective_kind: t.objective_kind,
        aux_variables: t.aux_variables,
        ortools_version: lib.version.clone(),
        ortools_library: lib.path.clone(),
    })
}
