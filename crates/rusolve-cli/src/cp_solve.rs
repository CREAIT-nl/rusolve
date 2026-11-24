//! Solve path for models with scheduling and combinatorial (CP) objects: intervals, no-overlap,
//! cumulative, all-different, element and table (`rusolve_io::cp`).
//!
//! No file format carries these objects, so the path is reached from the Python API and from the
//! MCP `build_model` tool, never from `rusolve solve <file>`. The router (`route.rs`) still names
//! the class (`cp`) and the engine for the envelope and for a model that has them next to
//! general constraints.
//!
//! ## The engine hook
//!
//! [`solve_with`] is the one place that picks the engine for a CP model: it asks
//! [`select_engine`] and runs the chosen engine.
//!
//! * [`CpEngine::Linearization`] (default): the objects are linearized into the port's `HighsLp`
//!   (`rusolve_io::cp_lin`), solved, checked on the user's objects and cut back to the user's
//!   columns.
//! * [`CpEngine::CpSat`] (`mip_engine=cpsat`): the user's model and objects go to OR-Tools CP-SAT
//!   (`rusolve-cpsat`, loaded at run time, `docs/design/cpsat-engine.md`) with no linearization.
//!   The point is checked on the rows and on the user's objects by [`crate::cpsat_solve::run_cp`]
//!   before it is kept, so a point that fails either is an error (RS0090), never an answer.
//!   The library missing, or a model CP-SAT cannot take exactly, is a typed error (RS0012).
//!
//! Both return a [`CpSolved`] whose `highs` holds the solution in the user's columns and whose
//! `report` comes from [`report`]. Python and MCP both call [`solve_with`].
//!
//! A model without CP objects never reaches this module.

use rusolve_core::diagnostic::Diagnostic;
use rusolve_ekk::lp_data::hconst::{HighsModelStatus, HighsVarType, kSolutionStatusNone};
use rusolve_ekk::lp_data::highs_lp::HighsLp;
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::HighsStatus;
use rusolve_io::cp_lin::{self, CpLinInput, CpLinOptions};
use rusolve_io::gc_lin::Extension;
use rusolve_io::{CpError, CpErrorKind, CpModel, CpTol, VariableBound};
use serde_json::{Value, json};

/// `engine` of the envelope for a CP model solved by linearization.
pub const ENGINE_LINEARIZATION: &str = "cp-linearization";

/// An engine that can solve a CP model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CpEngine {
    /// MIP linearization of the objects, solved by the HiGHS port.
    Linearization,
    /// OR-Tools CP-SAT, selected with `mip_engine=cpsat`.
    CpSat,
}

impl CpEngine {
    /// The name the envelope's `engine` carries.
    pub fn name(self) -> &'static str {
        match self {
            CpEngine::Linearization => ENGINE_LINEARIZATION,
            CpEngine::CpSat => crate::cpsat_solve::ENGINE,
        }
    }
}

/// The engine that solves a CP model under `options`: CP-SAT when `mip_engine=cpsat`, else the
/// linearization. A `mip_engine` value that is neither is refused by [`solve_with`].
pub fn select_engine(_cp: &CpModel, options: &HighsOptions) -> CpEngine {
    if options.mip_engine == "cpsat" {
        CpEngine::CpSat
    } else {
        CpEngine::Linearization
    }
}

/// Why [`solve_with`] did not solve: an object error (`RS4011`, `RS4013`) or an engine error.
#[derive(Debug, Clone, PartialEq)]
pub enum CpSolveError {
    Cp(CpError),
    /// The selected engine could not run: the code (`RS0004` bad `mip_engine`, `RS0012` library
    /// missing or model refused, `RS0090` internal) and the message.
    Engine {
        code: &'static str,
        message: String,
    },
}

impl From<CpError> for CpSolveError {
    fn from(e: CpError) -> Self {
        CpSolveError::Cp(e)
    }
}

impl std::fmt::Display for CpSolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CpSolveError::Cp(e) => write!(f, "{e}"),
            CpSolveError::Engine { message, .. } => write!(f, "{message}"),
        }
    }
}

impl CpSolveError {
    /// The registered diagnostic.
    pub fn diagnostic(&self) -> Diagnostic {
        match self {
            CpSolveError::Cp(e) => diagnostic(e),
            CpSolveError::Engine { code, message } => Diagnostic::error(code, message.clone())
                .with_help(format!("`rusolve --explain {code}` lists the cases")),
        }
    }
}

/// What `prepare` keeps for the solve and the report.
pub struct CpPlan {
    /// The model before any derived column or row.
    pub user_lp: HighsLp,
    pub cp: CpModel,
    pub ext: Extension,
}

/// The column bounds and integrality of `lp`, as [`CpModel::validate`] reads them.
pub fn bounds_of(lp: &HighsLp) -> Vec<VariableBound> {
    (0..lp.num_col_ as usize)
        .map(|j| VariableBound {
            lower: lp.col_lower_[j],
            upper: lp.col_upper_[j],
            integer: lp
                .integrality_
                .get(j)
                .is_some_and(|t| *t == HighsVarType::kInteger),
        })
        .collect()
}

/// Structural checks of `cp` against the model as it is now.
pub fn validate(lp: &HighsLp, cp: &CpModel) -> Result<(), CpError> {
    cp.validate(&bounds_of(lp), &lp.col_names_)
}

/// The registered diagnostic of a CP error: `RS4011` (invalid object) or `RS4013` (the fallback
/// cannot take it).
pub fn diagnostic(e: &CpError) -> Diagnostic {
    match e.kind {
        CpErrorKind::Invalid => Diagnostic::error("RS4011", e.to_string())
            .with_help("`rusolve --explain RS4011` lists the cases"),
        CpErrorKind::Unsupported => Diagnostic::error("RS4013", e.to_string()).with_help(format!(
            "`rusolve --explain RS4013` lists the limits of the MIP linearization; {}",
            cpsat_hint()
        )),
    }
}

/// The sentence that names the engine that takes what the linearization refuses.
fn cpsat_hint() -> String {
    if crate::cpsat_solve::available() {
        "CP-SAT is installed here: set `mip_engine=cpsat` to solve this model with it".into()
    } else {
        "OR-Tools CP-SAT takes it: install it (`install.sh --with-cpsat` or `pip install \
         'rusolve[cpsat]'`) and set `mip_engine=cpsat`"
            .into()
    }
}

/// Linearizes `cp` into `lp` (appending binaries and rows). Returns the plan holding the user's
/// model.
pub fn prepare(lp: &mut HighsLp, cp: CpModel, opt: CpLinOptions) -> Result<CpPlan, CpError> {
    let user_lp = lp.clone();
    let n = lp.num_col_ as usize;
    let integer: Vec<bool> = (0..n)
        .map(|j| {
            lp.integrality_
                .get(j)
                .is_some_and(|t| *t == HighsVarType::kInteger)
        })
        .collect();
    let names: Vec<String> = if lp.col_names_.len() == n {
        lp.col_names_.clone()
    } else {
        Vec::new()
    };
    let ext = cp_lin::linearize(
        &cp,
        &CpLinInput {
            col_lower: &lp.col_lower_,
            col_upper: &lp.col_upper_,
            col_integer: &integer,
            col_names: &names,
        },
        opt,
    )?;
    crate::gc_solve::extend_lp(lp, &ext, "cp_");
    Ok(CpPlan { user_lp, cp, ext })
}

/// The `cp` object of the envelope: objects by kind, how each was linearized, the checker's
/// violations on the user's objects, and the verdict. `checked` is false when there was no
/// solution to check.
pub fn report(
    cp: &CpModel,
    ext: Option<&Extension>,
    engine: &str,
    x: Option<(&[f64], &[String])>,
    tol: CpTol,
) -> Value {
    let violations: Vec<Value> = match x {
        Some((x, names)) => cp
            .check(x, names, tol)
            .into_iter()
            .map(
                |v| json!({"kind": v.kind, "name": v.name, "amount": v.amount, "detail": v.detail}),
            )
            .collect(),
        None => Vec::new(),
    };
    let mut r = json!({
        "engine": engine,
        "objects": cp.counts().iter().filter(|(_, c)| *c > 0).map(|(k, c)| (k.to_string(), json!(c))).collect::<serde_json::Map<_, _>>(),
        "checked": x.is_some(),
        "violations": violations,
        "verified": x.is_some() && violations.is_empty(),
    });
    if let Some(e) = ext {
        r["routes"] = Value::Array(crate::gc_solve::routes_json(e));
        r["added_columns"] = json!(e.cols.len());
        r["added_rows"] = json!(e.rows.len());
    }
    r
}

/// After the solve of the linearization: check the incumbent against the user's objects, repair
/// the trickle flow once, cut the solution to the user's size, and never leave an answer that
/// fails the checker looking feasible. Returns the `cp` report.
///
/// Trickle flow: a binary that is `1e-7` off its value lets a big-M row hold a little slack. The
/// repair fixes every integer column at its rounded value and solves once more (the model is then
/// an LP over the continuous columns); the event is reported either way. With `downgrade` an
/// answer that still fails the checker is turned into "no solution" and the status into Unknown.
pub fn finish(plan: &CpPlan, highs: &mut Highs, tol: CpTol, downgrade: bool) -> Value {
    let n = plan.user_lp.num_col_ as usize;
    let m = plan.user_lp.num_row_ as usize;
    let names = &plan.user_lp.col_names_;
    let usable = |h: &Highs| h.solution_.value_valid && h.solution_.col_value.len() >= n;
    let mut trickle = None;
    if usable(highs) {
        let before = plan.cp.check(&highs.solution_.col_value[..n], names, tol);
        if !before.is_empty() {
            let total = highs.solution_.col_value.len();
            for j in 0..total {
                let integer = highs
                    .get_lp()
                    .integrality_
                    .get(j)
                    .is_some_and(|t| *t == HighsVarType::kInteger);
                if integer {
                    let lp = highs.get_lp();
                    let r = highs.solution_.col_value[j]
                        .round()
                        .clamp(lp.col_lower_[j], lp.col_upper_[j]);
                    highs.change_col_bounds(j as i32, r, r);
                }
            }
            highs.solution_ = Default::default();
            highs.run();
            let after = if usable(highs) {
                plan.cp
                    .check(&highs.solution_.col_value[..n], names, tol)
                    .len()
            } else {
                0
            };
            trickle = Some(json!({
                "violations_before": before.len(),
                "violations_after": after,
                "status": format!("{:?}", highs.get_model_status()),
            }));
        }
    }
    let mut rep = if usable(highs) {
        report(
            &plan.cp,
            Some(&plan.ext),
            ENGINE_LINEARIZATION,
            Some((&highs.solution_.col_value[..n], names)),
            tol,
        )
    } else {
        report(&plan.cp, Some(&plan.ext), ENGINE_LINEARIZATION, None, tol)
    };
    if let Some(t) = trickle {
        rep["trickle_flow"] = t;
    }
    // cut the solution to the user's columns and rows
    let s = &mut highs.solution_;
    s.col_value.truncate(n);
    s.row_value.truncate(m);
    s.col_dual.truncate(n);
    s.row_dual.truncate(m);
    if downgrade && rep["checked"] == true && rep["verified"] == false {
        reject_answer(highs);
    }
    rep
}

/// An answer that fails the checker is never reported as feasible: no solution, status Unknown.
pub fn reject_answer(highs: &mut Highs) {
    highs.solution_.value_valid = false;
    highs.model_status_ = HighsModelStatus::kUnknown;
    highs.info_.primal_solution_status = kSolutionStatusNone;
}

/// What a CP solve hands back.
pub struct CpSolved {
    /// The model as the user built it (counts, names, column count of the report).
    pub user_lp: HighsLp,
    /// The solver after the solve, its solution cut to the user's columns and rows.
    pub highs: Highs,
    pub run: HighsStatus,
    /// The `cp` object of the envelope.
    pub report: Value,
    /// The engine that produced the answer (`cp-linearization` or `cp-sat`).
    pub engine: &'static str,
}

/// THE ENGINE HOOK: solves `lp` with the CP objects `cp` on the engine [`select_engine`] picks.
/// `run` runs the HiGHS port on the linearized model (the Python API passes a closure that
/// releases the GIL; [`solve`] passes `Highs::run`); an engine that does not use the port ignores
/// it.
pub fn solve_with(
    lp: HighsLp,
    cp: CpModel,
    options: HighsOptions,
    run: impl FnOnce(&mut Highs) -> HighsStatus,
) -> Result<CpSolved, CpSolveError> {
    validate(&lp, &cp)?;
    if !matches!(options.mip_engine.as_str(), "highs" | "cpsat") {
        return Err(CpSolveError::Engine {
            code: "RS0004",
            message: format!(
                "mip_engine={}: the engines are \"highs\" (the port, default) and \"cpsat\" (Google OR-Tools CP-SAT)",
                options.mip_engine
            ),
        });
    }
    match select_engine(&cp, &options) {
        CpEngine::CpSat => {
            let tol = CpTol {
                integrality: options.mip_feasibility_tolerance,
            };
            let mut highs = Highs::with_options(options);
            highs.pass_model(lp.clone());
            highs.reset_scheduler_for_threads();
            let cpsat =
                crate::cpsat_solve::run_cp(&mut highs, Some(&cp), "cp model").map_err(|f| {
                    CpSolveError::Engine {
                        code: if f.exit() == 4 { "RS0090" } else { "RS0012" },
                        message: f.message(),
                    }
                })?;
            let n = lp.num_col_ as usize;
            let usable = highs.solution_.value_valid && highs.solution_.col_value.len() >= n;
            let mut report = report(
                &cp,
                None,
                crate::cpsat_solve::ENGINE,
                usable.then(|| (&highs.solution_.col_value[..n], lp.col_names_.as_slice())),
                tol,
            );
            report["cpsat"] = cpsat;
            Ok(CpSolved {
                user_lp: lp,
                highs,
                run: HighsStatus::kOk,
                report,
                engine: crate::cpsat_solve::ENGINE,
            })
        }
        CpEngine::Linearization => {
            let mut lin = lp;
            let plan = prepare(&mut lin, cp, CpLinOptions::default())?;
            let tol = CpTol {
                integrality: options.mip_feasibility_tolerance,
            };
            let mut highs = Highs::with_options(options);
            highs.pass_model(lin);
            let run_status = run(&mut highs);
            let report = finish(&plan, &mut highs, tol, true);
            Ok(CpSolved {
                user_lp: plan.user_lp,
                highs,
                run: run_status,
                report,
                engine: ENGINE_LINEARIZATION,
            })
        }
    }
}

/// [`solve_with`] on the calling thread.
pub fn solve(lp: HighsLp, cp: CpModel, options: HighsOptions) -> Result<CpSolved, CpSolveError> {
    solve_with(lp, cp, options, |h| h.run())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_io::{AllDifferent, CpTerm, Cumulative, Element, Interval, NoOverlap, Table};

    impl CpSolveError {
        /// The object error; any other is a failed test.
        fn cp(self) -> CpError {
            match self {
                CpSolveError::Cp(e) => e,
                other => panic!("expected an object error, got {other}"),
            }
        }
    }

    fn lp_of(text: &str) -> HighsLp {
        let dir = std::env::temp_dir().join(format!("rusolve-cp-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        // One file per call: tests run in parallel, and two models of the same length used to
        // share a file name and overwrite each other.
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = dir.join(format!("m{n}.lp"));
        std::fs::write(&path, text).unwrap();
        crate::mip_highs::read_mip_full(path.to_str().unwrap())
            .unwrap()
            .0
    }

    fn col(lp: &HighsLp, name: &str) -> usize {
        lp.col_names_.iter().position(|n| n == name).unwrap()
    }

    fn opts() -> HighsOptions {
        HighsOptions {
            output_flag: false,
            threads: 1,
            ..HighsOptions::default()
        }
    }

    fn iv(name: &str, s: usize, size: i64, e: usize) -> Interval {
        Interval {
            name: name.into(),
            start: s,
            size: CpTerm::Const(size),
            end: e,
            presence: None,
        }
    }

    /// Minimum objective over every integer point of the `decision` columns' boxes, by
    /// enumeration. `fill` completes a point (the other columns, which the model's own rows
    /// determine) and returns its objective, or `None` when a row is violated; the user's objects
    /// are then the only other constraint test. Keep the enumeration small: the product of the
    /// decision boxes is what runs.
    fn brute(
        lp: &HighsLp,
        cp: &CpModel,
        decision: &[usize],
        fill: &dyn Fn(&mut Vec<f64>) -> Option<f64>,
    ) -> Option<f64> {
        let n = lp.num_col_ as usize;
        let ranges: Vec<(i64, i64)> = decision
            .iter()
            .map(|&j| (lp.col_lower_[j] as i64, lp.col_upper_[j] as i64))
            .collect();
        assert!(
            ranges
                .iter()
                .map(|r| (r.1 - r.0 + 1) as u64)
                .product::<u64>()
                <= 2_000_000,
            "brute force box too large"
        );
        let mut x: Vec<i64> = ranges.iter().map(|r| r.0).collect();
        let mut best: Option<f64> = None;
        loop {
            let mut xf = vec![0.0; n];
            for (k, &j) in decision.iter().enumerate() {
                xf[j] = x[k] as f64;
            }
            if let Some(obj) = fill(&mut xf)
                && cp.check(&xf, &lp.col_names_, CpTol::default()).is_empty()
            {
                best = Some(best.map_or(obj, |b: f64| b.min(obj)));
            }
            let mut k = 0;
            loop {
                if k == x.len() {
                    return best;
                }
                if x[k] < ranges[k].1 {
                    x[k] += 1;
                    break;
                }
                x[k] = ranges[k].0;
                k += 1;
            }
        }
    }

    /// Every column is a decision column and the objective is `cost . x`.
    fn brute_all(lp: &HighsLp, cp: &CpModel, cost: &[f64]) -> Option<f64> {
        let all: Vec<usize> = (0..lp.num_col_ as usize).collect();
        brute(lp, cp, &all, &|x: &mut Vec<f64>| {
            Some(x.iter().zip(cost).map(|(a, b)| a * b).sum())
        })
    }

    /// Solves through the hook and compares with the brute-force optimum `want` on the user's
    /// objects.
    fn agrees(lp: HighsLp, cp: CpModel, want: Option<f64>) -> CpSolved {
        let s = solve(lp, cp, opts()).unwrap();
        match want {
            None => assert!(!s.highs.solution_.value_valid, "{}", s.report),
            Some(w) => {
                assert!(s.highs.solution_.value_valid, "{}", s.report);
                assert_eq!(s.report["verified"], true, "{}", s.report);
                let got = s.highs.get_info().objective_function_value;
                assert!((got - w).abs() < 1e-6, "got {got} want {w}: {}", s.report);
            }
        }
        s
    }

    #[test]
    fn job_shop_two_machines_matches_brute_force() {
        // 3 jobs, 2 machines (the 2-job-2-machine core of a job shop, with a third job on
        // machine 1 only): minimize the makespan
        let lp = lp_of(
            "Minimize\n obj: mk\nSubject To\n m1: mk - e1 >= 0\n m2: mk - e2 >= 0\n m3: mk - e3 >= 0\n m4: mk - e4 >= 0\n m5: mk - e5 >= 0\n pa: s2 - e1 >= 0\n pb: s4 - e3 >= 0\nBounds\n 0 <= s1 <= 8\n 0 <= e1 <= 8\n 0 <= s2 <= 8\n 0 <= e2 <= 8\n 0 <= s3 <= 8\n 0 <= e3 <= 8\n 0 <= s4 <= 8\n 0 <= e4 <= 8\n 0 <= s5 <= 8\n 0 <= e5 <= 8\n 0 <= mk <= 8\nGeneral\n s1 e1 s2 e2 s3 e3 s4 e4 s5 e5 mk\nEnd\n",
        );
        let c = |n: &str| col(&lp, n);
        let cp = CpModel {
            intervals: vec![
                iv("a1", c("s1"), 2, c("e1")),
                iv("a2", c("s2"), 1, c("e2")),
                iv("b1", c("s3"), 2, c("e3")),
                iv("b2", c("s4"), 2, c("e4")),
                iv("c1", c("s5"), 1, c("e5")),
            ],
            no_overlap: vec![
                NoOverlap {
                    name: "machine1".into(),
                    intervals: vec![0, 2, 4],
                },
                NoOverlap {
                    name: "machine2".into(),
                    intervals: vec![1, 3],
                },
            ],
            ..Default::default()
        };
        // decision columns: the five starts; the ends and the makespan follow from the rows
        let starts: Vec<usize> = ["s1", "s2", "s3", "s4", "s5"]
            .iter()
            .map(|n| c(n))
            .collect();
        let ends: Vec<usize> = ["e1", "e2", "e3", "e4", "e5"]
            .iter()
            .map(|n| c(n))
            .collect();
        let sizes = [2.0, 1.0, 2.0, 2.0, 1.0];
        let (mk, n) = (c("mk"), lp.num_col_ as usize);
        let fill = |x: &mut Vec<f64>| {
            x.resize(n, 0.0);
            for k in 0..5 {
                x[ends[k]] = x[starts[k]] + sizes[k];
            }
            let m = ends.iter().map(|&j| x[j]).fold(0.0, f64::max);
            x[mk] = m;
            // rows: a2 after a1, b2 after b1; every column within [0, 8]
            let ok = x[starts[1]] >= x[ends[0]]
                && x[starts[3]] >= x[ends[2]]
                && x.iter().all(|&v| v <= 8.0);
            ok.then_some(m)
        };
        let want = brute(&lp, &cp, &starts, &fill);
        assert!(want.is_some());
        let s = agrees(lp, cp, want);
        assert_eq!(s.engine, "cp-linearization");
        assert_eq!(
            s.highs.solution_.col_value.len(),
            11,
            "cut to the user's columns"
        );
    }

    #[test]
    fn rcpsp_cumulative_matches_brute_force() {
        // 4 tasks of sizes 2, 1, 2, 1 and demands 2, 2, 1, 2 on a resource of capacity 3
        let lp = lp_of(
            "Minimize\n obj: mk\nSubject To\n c1: mk - e1 >= 0\n c2: mk - e2 >= 0\n c3: mk - e3 >= 0\n c4: mk - e4 >= 0\n p: s3 - e1 >= 0\nBounds\n 0 <= s1 <= 4\n 0 <= e1 <= 6\n 0 <= s2 <= 4\n 0 <= e2 <= 6\n 0 <= s3 <= 4\n 0 <= e3 <= 6\n 0 <= s4 <= 4\n 0 <= e4 <= 6\n 0 <= mk <= 6\nGeneral\n s1 e1 s2 e2 s3 e3 s4 e4 mk\nEnd\n",
        );
        let c = |n: &str| col(&lp, n);
        let cp = CpModel {
            intervals: vec![
                iv("t1", c("s1"), 2, c("e1")),
                iv("t2", c("s2"), 1, c("e2")),
                iv("t3", c("s3"), 2, c("e3")),
                iv("t4", c("s4"), 1, c("e4")),
            ],
            cumulative: vec![Cumulative {
                name: "crew".into(),
                intervals: vec![0, 1, 2, 3],
                demands: vec![2, 2, 1, 2].into_iter().map(CpTerm::Const).collect(),
                capacity: CpTerm::Const(3),
            }],
            ..Default::default()
        };
        let starts: Vec<usize> = ["s1", "s2", "s3", "s4"].iter().map(|n| c(n)).collect();
        let ends: Vec<usize> = ["e1", "e2", "e3", "e4"].iter().map(|n| c(n)).collect();
        let sizes = [2.0, 1.0, 2.0, 1.0];
        let (mk, n) = (c("mk"), lp.num_col_ as usize);
        let fill = |x: &mut Vec<f64>| {
            x.resize(n, 0.0);
            for k in 0..4 {
                x[ends[k]] = x[starts[k]] + sizes[k];
            }
            let m = ends.iter().map(|&j| x[j]).fold(0.0, f64::max);
            x[mk] = m;
            // row: task 3 after task 1; every column within its bound
            let ok = x[starts[2]] >= x[ends[0]] && x.iter().all(|&v| v <= 6.0);
            ok.then_some(m)
        };
        let want = brute(&lp, &cp, &starts, &fill);
        assert!(want.is_some());
        let s = agrees(lp, cp, want);
        let routes = s.report["routes"].as_array().unwrap();
        assert!(
            routes
                .iter()
                .any(|r| r["route"] == "cumulative_time_indexed")
        );
    }

    #[test]
    fn optional_intervals_with_presence() {
        // two optional jobs on one machine, each worth 3 (a) and 2 (b) when present; the machine
        // is free for 3 time units: only one fits next to the other's size
        let lp = lp_of(
            "Maximize\n obj: 3 pa + 2 pb\nSubject To\nBounds\n 0 <= sa <= 3\n 0 <= ea <= 3\n 0 <= sb <= 3\n 0 <= eb <= 3\nGeneral\n sa ea sb eb\nBinaries\n pa pb\nEnd\n",
        );
        let c = |n: &str| col(&lp, n);
        let mut a = iv("a", c("sa"), 2, c("ea"));
        a.presence = Some(c("pa"));
        let mut b = iv("b", c("sb"), 2, c("eb"));
        b.presence = Some(c("pb"));
        let cp = CpModel {
            intervals: vec![a, b],
            no_overlap: vec![NoOverlap {
                name: "m".into(),
                intervals: vec![0, 1],
            }],
            ..Default::default()
        };
        let mut cost = vec![0.0; lp.num_col_ as usize];
        cost[c("pa")] = -3.0;
        cost[c("pb")] = -2.0;
        // brute force minimizes; the model maximizes the same thing
        let want = brute_all(&lp, &cp, &cost).unwrap();
        assert_eq!(want, -3.0);
        let s = solve(lp, cp, opts()).unwrap();
        assert_eq!(s.report["verified"], true, "{}", s.report);
        assert!((s.highs.get_info().objective_function_value - 3.0).abs() < 1e-6);
    }

    #[test]
    fn all_different_element_table_match_brute_force() {
        let lp = lp_of(
            "Minimize\n obj: 3 x + 2 y + z + 5 i + t\nSubject To\nBounds\n 0 <= x <= 3\n 0 <= y <= 3\n 0 <= z <= 3\n 0 <= i <= 2\n 0 <= t <= 9\nGeneral\n x y z i t\nEnd\n",
        );
        let c = |n: &str| col(&lp, n);
        let cp = CpModel {
            all_different: vec![AllDifferent {
                name: "d".into(),
                cols: vec![c("x"), c("y"), c("z")],
            }],
            element: vec![Element {
                name: "e".into(),
                index: c("i"),
                values: vec![CpTerm::Const(4), CpTerm::Col(c("x")), CpTerm::Const(7)],
                target: c("t"),
            }],
            table: vec![
                Table {
                    name: "tab".into(),
                    cols: vec![c("x"), c("i")],
                    tuples: vec![vec![1, 0], vec![2, 1], vec![3, 2], vec![0, 2]],
                    negated: false,
                },
                Table {
                    name: "neg".into(),
                    cols: vec![c("y"), c("z")],
                    tuples: vec![vec![0, 1], vec![1, 0], vec![2, 3]],
                    negated: true,
                },
            ],
            ..Default::default()
        };
        let mut cost = vec![0.0; lp.num_col_ as usize];
        for (n, v) in [("x", 3.0), ("y", 2.0), ("z", 1.0), ("i", 5.0), ("t", 1.0)] {
            cost[c(n)] = v;
        }
        let want = brute_all(&lp, &cp, &cost);
        assert!(want.is_some());
        let s = agrees(lp, cp, want);
        let routes: Vec<&str> = s.report["routes"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| r["route"].as_str().unwrap())
            .collect();
        for want in [
            "all_different_value_binaries",
            "element_selection",
            "table_tuple_binaries",
            "table_negated_disjunctive",
        ] {
            assert!(routes.contains(&want), "{routes:?}");
        }
    }

    #[test]
    fn infeasible_cp_is_not_reported_feasible() {
        // two size-3 jobs in a window of 5 on one machine
        let lp = lp_of(
            "Minimize\n obj: s1\nSubject To\nBounds\n 0 <= s1 <= 2\n 0 <= e1 <= 5\n 0 <= s2 <= 2\n 0 <= e2 <= 5\nGeneral\n s1 e1 s2 e2\nEnd\n",
        );
        let c = |n: &str| col(&lp, n);
        let cp = CpModel {
            intervals: vec![iv("a", c("s1"), 3, c("e1")), iv("b", c("s2"), 3, c("e2"))],
            no_overlap: vec![NoOverlap {
                name: "m".into(),
                intervals: vec![0, 1],
            }],
            ..Default::default()
        };
        let mut cost = vec![0.0; 4];
        cost[c("s1")] = 1.0;
        let want = brute_all(&lp, &cp, &cost);
        assert_eq!(want, None);
        let s = agrees(lp, cp, want);
        assert_eq!(s.report["checked"], false);
        assert_eq!(s.report["verified"], false);
    }

    #[test]
    fn refusals_are_typed() {
        let wide = "Minimize\n obj: s1\nSubject To\nBounds\n 0 <= s1 <= 5000\n 0 <= e1 <= 5002\nGeneral\n s1 e1\nEnd\n";
        let lp = lp_of(wide);
        let c = |n: &str| col(&lp, n);
        let cum = |size, demand, cap| CpModel {
            intervals: vec![Interval {
                name: "a".into(),
                start: c("s1"),
                size,
                end: c("e1"),
                presence: None,
            }],
            cumulative: vec![Cumulative {
                name: "r".into(),
                intervals: vec![0],
                demands: vec![demand],
                capacity: cap,
            }],
            ..Default::default()
        };
        // horizon above the limit
        let e = solve(
            lp.clone(),
            cum(CpTerm::Const(2), CpTerm::Const(1), CpTerm::Const(1)),
            opts(),
        )
        .err()
        .unwrap()
        .cp();
        assert_eq!(e.kind, CpErrorKind::Unsupported);
        assert!(
            e.message.contains("horizon") && e.message.contains("CP-SAT"),
            "{e}"
        );
        assert_eq!(diagnostic(&e).code, "RS4013");
        // variable size (the box is narrowed first so the horizon is not what refuses)
        let lp2 = lp_of(
            "Minimize\n obj: s1\nSubject To\nBounds\n 0 <= s1 <= 5\n 0 <= e1 <= 9\n 1 <= z <= 3\nGeneral\n s1 e1 z\nEnd\n",
        );
        let c2 = |n: &str| col(&lp2, n);
        let cp = CpModel {
            intervals: vec![Interval {
                name: "a".into(),
                start: c2("s1"),
                size: CpTerm::Col(c2("z")),
                end: c2("e1"),
                presence: None,
            }],
            cumulative: vec![Cumulative {
                name: "r".into(),
                intervals: vec![0],
                demands: vec![CpTerm::Const(1)],
                capacity: CpTerm::Const(1),
            }],
            ..Default::default()
        };
        let e = solve(lp2.clone(), cp, opts()).err().unwrap().cp();
        assert!(e.message.contains("variable size"), "{e}");
        // a no-overlap size variable that can be 0
        let lp3 = lp_of(
            "Minimize\n obj: s1\nSubject To\nBounds\n 0 <= s1 <= 5\n 0 <= e1 <= 9\n 0 <= s2 <= 5\n 0 <= e2 <= 9\n 0 <= z <= 3\nGeneral\n s1 e1 s2 e2 z\nEnd\n",
        );
        let c3 = |n: &str| col(&lp3, n);
        let cp = CpModel {
            intervals: vec![
                Interval {
                    name: "a".into(),
                    start: c3("s1"),
                    size: CpTerm::Col(c3("z")),
                    end: c3("e1"),
                    presence: None,
                },
                iv("b", c3("s2"), 2, c3("e2")),
            ],
            no_overlap: vec![NoOverlap {
                name: "m".into(),
                intervals: vec![0, 1],
            }],
            ..Default::default()
        };
        let e = solve(lp3, cp, opts()).err().unwrap().cp();
        assert!(e.message.contains("can be 0"), "{e}");
        // an invalid object is RS4011
        let cp = CpModel {
            no_overlap: vec![NoOverlap {
                name: "m".into(),
                intervals: vec![3],
            }],
            ..Default::default()
        };
        let e = solve(lp2, cp, opts()).err().unwrap().cp();
        assert_eq!(e.kind, CpErrorKind::Invalid);
        assert_eq!(diagnostic(&e).code, "RS4011");
    }

    #[test]
    fn a_failing_answer_is_never_feasible() {
        let lp = lp_of(
            "Minimize\n obj: s1\nSubject To\nBounds\n 0 <= s1 <= 2\n 0 <= e1 <= 5\n 0 <= s2 <= 2\n 0 <= e2 <= 5\nGeneral\n s1 e1 s2 e2\nEnd\n",
        );
        let c = |n: &str| col(&lp, n);
        let cp = CpModel {
            intervals: vec![iv("a", c("s1"), 3, c("e1")), iv("b", c("s2"), 1, c("e2"))],
            no_overlap: vec![NoOverlap {
                name: "m".into(),
                intervals: vec![0, 1],
            }],
            ..Default::default()
        };
        let mut plan_lp = lp.clone();
        let plan = prepare(&mut plan_lp, cp, CpLinOptions::default()).unwrap();
        let mut h = Highs::with_options(opts());
        h.pass_model(plan_lp);
        h.run();
        assert!(h.solution_.value_valid);
        // tamper: both jobs start at 0
        for (name, v) in [("s1", 0.0), ("e1", 3.0), ("s2", 0.0), ("e2", 1.0)] {
            h.solution_.col_value[c(name)] = v;
        }
        let rep = finish(&plan, &mut h, CpTol::default(), true);
        // the repair fixes the integer columns at their values; they overlap, so the re-solve is
        // infeasible and the answer is dropped, never reported feasible
        assert_eq!(rep["trickle_flow"]["violations_before"], 1, "{rep}");
        assert_eq!(rep["verified"], false, "{rep}");
        assert!(!h.solution_.value_valid, "{rep}");
        // an incumbent that fails the checker is turned into "no solution": reject_answer
        h.solution_.value_valid = true;
        reject_answer(&mut h);
        assert!(!h.solution_.value_valid);
        assert_eq!(h.get_model_status(), HighsModelStatus::kUnknown);
    }
}
