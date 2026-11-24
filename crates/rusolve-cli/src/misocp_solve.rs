//! Mixed-integer SOCP and convex MIQCQP by LP-based outer approximation: the glue between the
//! conic model files, the Clarabel route (`conic_solve.rs`) and the port's MIP search with
//! `rusolve_ekk::mip::oa` (which holds the cone cuts, the separator and the exactness oracle, and
//! cites the sources). Design: `docs/design/misocp-oa.md`.
//!
//! Only `rusolve solve` / `rusolve mip` with `misocp_oa` on and a model with integer columns and
//! cones or quadratic rows get here (`route.rs`, class `mi_conic`); with the lever off the model
//! is refused as before. The flow of [`run`]:
//!
//! 1. [`build`]: every cone, every convex quadratic row and the Hessian of the objective become
//!    cones `||y(x)|| <= t(x)` with `y`, `t` affine (`Q = LLᵀ`, `½||Lᵀx||² <= u` is
//!    `||(u - 1, √2 Lᵀx)|| <= u + 1`, Lobo et al. 1998). The LP keeps the linear rows and the
//!    integrality; the cones are rows added as cuts;
//! 2. the continuous SOCP relaxation is solved by Clarabel; its optimum is the first cut point
//!    (and, rounded and polished, the first incumbent);
//! 3. the outer-approximation search runs on the LP (tree separator, lazy-constraint loop); every
//!    candidate is polished: the integers are fixed and the continuous SOCP is solved again, which
//!    gives a feasible point that the cuts alone never reach exactly;
//! 4. the answer is the best checker-feasible point (a polished one or the accepted candidate);
//!    it is reported `Optimal` only when the checker confirms it and the bound proves it.

use rusolve_ekk::lp_data::highs_option_records as registry;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::oa::{
    Aff, OaConfig, OaForm, OaStats, OaStatsHandle, Polisher, QuadForm, add_cone_rows, run_oa_with,
};
use rusolve_ekk::prelude::*;
use rusolve_io::conic::ConeKind;
use rusolve_io::{Model, MpsFileOrder, ObjSense as IoSense};
use serde_json::{Value, json};
use std::collections::HashSet;
use std::process::ExitCode;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::check_conic::{DEFAULT_TOL, check_conic_solution, row_interval};
use crate::conic_solve::{ConicError, col_name, refusal_error, row_name, solve_model};
use crate::envelope::{NativeSolved, native_envelope};
use crate::quad_convexity::{Refusal, factor_psd};

/// Cones with at least this many `y` components are written in the extended formulation by
/// `misocp_oa_disaggregate` (a cone of two components is already three columns).
const DISAGG_MIN_Y: usize = 3;
const SQRT2: f64 = std::f64::consts::SQRT_2;

/// The settings of one run.
#[derive(Clone, Debug, Default)]
pub struct Settings {
    pub disaggregate: bool,
    pub no_polish: bool,
    pub lazy_only: bool,
}

/// What [`build`] produces: the OA form and the column of the objective's epigraph, if any.
pub struct Built {
    pub form: OaForm,
    pub obj_eta: Option<usize>,
    pub cones: usize,
}

fn lin_rows(m: &Model) -> Vec<Vec<(usize, f64)>> {
    let (cp, ri, va) = (m.matrix.col_ptr(), m.matrix.row_idx(), m.matrix.values());
    let mut lin: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m.num_constraints];
    for j in 0..m.num_variables {
        for k in cp[j]..cp[j + 1] {
            lin[ri[k] as usize].push((j, va[k]));
        }
    }
    lin
}

/// The cone `½||Lᵀx||² <= u` for `u = bound + sum_j ucoef_j x_j` as `(t, y)`.
fn quad_cone(
    f: &crate::quad_convexity::Factor,
    ucoef: &[(usize, f64)],
    bound: f64,
) -> (Aff, Vec<Aff>) {
    let t = Aff::new(ucoef, bound + 1.0);
    let mut y = vec![Aff::new(ucoef, bound - 1.0)];
    for cn in &f.l {
        let pairs: Vec<(usize, f64)> = f
            .vars
            .iter()
            .zip(cn)
            .filter(|&(_, &v)| v != 0.0)
            .map(|(&j, &v)| (j, SQRT2 * v))
            .collect();
        y.push(Aff::new(&pairs, 0.0));
    }
    (t, y)
}

/// The OA form of the cone constraints of `m`, or the typed refusal.
pub fn build(m: &Model, disaggregate: bool) -> Result<Built, ConicError> {
    let n = m.num_variables;
    for cn in &m.conic.cones {
        if !matches!(cn.kind, ConeKind::Soc | ConeKind::RotatedSoc) {
            return Err(ConicError::UnsupportedCone {
                cone: cn.name.clone(),
                kind: cn.kind.label(),
            });
        }
    }
    let maximize = m.obj_sense == IoSense::Maximize;
    let sgn = if maximize { -1.0 } else { 1.0 };
    let mut form = OaForm::new(&QuadForm { n, entries: vec![] }, maximize, false);
    let mut pending: Vec<(String, Aff, Vec<Aff>)> = Vec::new();
    let mut obj_eta = None;

    // objective Hessian: eta >= ½ xᵀ(sgn Q)x, as a cone on a new column
    if !m.hessian.is_empty() {
        let lower: Vec<_> = m
            .hessian
            .iter()
            .filter(|&&(i, j, _)| i >= j)
            .copied()
            .collect();
        match factor_psd(&lower, sgn) {
            Ok(f) => {
                let eta = form.add_extra_col(sgn, 0.0, kHighsInf, "oa_obj_eta".to_string());
                obj_eta = Some(eta);
                let (t, y) = quad_cone(&f, &[(eta, 1.0)], 0.0);
                pending.push(("obj".to_string(), t, y));
            }
            Err(Refusal::Indefinite { var, .. }) => {
                return Err(ConicError::NonconvexObjective {
                    var: col_name(m, var),
                    sense: if maximize {
                        "maximization"
                    } else {
                        "minimization"
                    },
                });
            }
            Err(Refusal::TooLarge { .. }) => {
                return Err(ConicError::Solver(
                    "objective Hessian is too large to test".into(),
                ));
            }
        }
    }

    // quadratic rows: each finite side with a quadratic term is a cone
    let lin = lin_rows(m);
    let mut quad_of = vec![usize::MAX; m.num_constraints];
    for (k, qr) in m.conic.quad_rows.iter().enumerate() {
        quad_of[qr.row] = k;
    }
    for i in 0..m.num_constraints {
        let k = quad_of[i];
        if k == usize::MAX {
            continue;
        }
        let (lo, hi) = row_interval(m, i);
        let name = row_name(m, i);
        let mut sides: Vec<(&'static str, f64, f64)> = Vec::new();
        if hi.is_finite() {
            sides.push(("upper", 1.0, hi));
        }
        if lo.is_finite() {
            sides.push(("lower", -1.0, lo));
        }
        for (side, sign, bound) in sides {
            let f = match factor_psd(&m.conic.quad_rows[k].q, sign) {
                Ok(f) => f,
                Err(r) => return Err(refusal_error(m, &name, side, r)),
            };
            if f.l.is_empty() {
                continue;
            }
            if lo == hi {
                return Err(ConicError::EqualityQuadRow { row: name });
            }
            // sign (a'x + ½xᵀQx) <= sign bound  <=>  ½||Lᵀx||² <= u,  u = sign bound - sign a'x
            let ucoef: Vec<(usize, f64)> = lin[i].iter().map(|&(j, v)| (j, -sign * v)).collect();
            let (t, y) = quad_cone(&f, &ucoef, sign * bound);
            pending.push((format!("{name}_{side}"), t, y));
        }
    }

    // user cones
    for cn in &m.conic.cones {
        let v = &cn.vars;
        let e = |j: usize, s: f64| Aff::new(&[(j, s)], 0.0);
        let (t, y) = match cn.kind {
            ConeKind::Soc => (e(v[0], 1.0), v[1..].iter().map(|&j| e(j, 1.0)).collect()),
            _ => {
                // 2 x0 x1 >= ||w||²  <=>  ||(x0 - x1, √2 w)|| <= x0 + x1
                let mut y = vec![Aff::new(&[(v[0], 1.0), (v[1], -1.0)], 0.0)];
                y.extend(v[2..].iter().map(|&j| e(j, SQRT2)));
                (Aff::new(&[(v[0], 1.0), (v[1], 1.0)], 0.0), y)
            }
        };
        pending.push((cn.name.clone(), t, y));
    }

    let cones = pending.len();
    for (name, t, y) in pending {
        if disaggregate && y.len() >= DISAGG_MIN_Y {
            form.add_cone_disaggregated(&name, t, y);
        } else {
            form.add_cone(t, y);
        }
    }
    Ok(Built {
        form,
        obj_eta,
        cones,
    })
}

/// The LP the port solves: the user's linear rows, bounds and integrality (the quadratic part of
/// a quadratic row is dropped, which is valid because every finite side is convex), the heads of
/// the cones bounded below by 0, then the extra columns of `form`.
fn build_lp(m: &Model, form: &OaForm) -> HighsLp {
    let order = MpsFileOrder::of_csc(m);
    let mut lp = crate::ekk::highs_lp_from_model(m, &order);
    crate::ekk::fill_names(&mut lp, m);
    lp.integrality_ = m
        .bounds
        .iter()
        .map(|b| {
            if b.integer {
                HighsVarType::kInteger
            } else {
                HighsVarType::kContinuous
            }
        })
        .collect();
    for cn in &m.conic.cones {
        let heads: &[usize] = match cn.kind {
            ConeKind::Soc => &cn.vars[..1],
            _ => &cn.vars[..2],
        };
        for &j in heads {
            if lp.col_lower_[j] < 0.0 {
                lp.col_lower_[j] = 0.0;
            }
        }
    }
    form.extend_lp(&mut lp);
    lp
}

/// What the polish learned.
#[derive(Default)]
struct Best {
    /// Best feasible objective in the model's sense and its point (model columns).
    obj: Option<f64>,
    x: Vec<f64>,
    calls: u64,
    solved: u64,
    improved: u64,
    seconds: f64,
}

type BestHandle = Arc<Mutex<Best>>;

fn lock<T>(m: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    m.lock().unwrap_or_else(|e| e.into_inner())
}

fn better(maximize: bool, a: f64, b: f64) -> bool {
    if maximize { a > b } else { a < b }
}

/// Everything the polish needs, shared with the closure the oracle owns.
struct PolishCtx {
    model: Arc<Model>,
    form: Arc<OaForm>,
    obj_eta: Option<usize>,
    ints: Vec<usize>,
    best: BestHandle,
    deadline: Option<Instant>,
    seen: HashSet<Vec<i64>>,
}

impl PolishCtx {
    /// Fixes the integer columns of `cand` and solves the continuous SOCP. The point is kept when
    /// the checker confirms it. Returns the point in OA columns (for the next solve's start).
    fn polish(&mut self, cand: &[f64]) -> Option<Vec<f64>> {
        let m = &*self.model;
        let key: Vec<i64> = self.ints.iter().map(|&j| cand[j].round() as i64).collect();
        if !self.seen.insert(key) {
            return None;
        }
        let left = self
            .deadline
            .map(|d| d.saturating_duration_since(Instant::now()).as_secs_f64());
        if left == Some(0.0) {
            return None;
        }
        let t0 = Instant::now();
        let mut fixed = m.clone();
        for &j in &self.ints {
            let b = &mut fixed.bounds[j];
            let v = cand[j].round().clamp(b.lower, b.upper);
            b.lower = v;
            b.upper = v;
            b.integer = false;
        }
        for b in &mut fixed.bounds {
            b.integer = false;
        }
        let solved = solve_model(&fixed, left).ok();
        let mut out = None;
        let mut g = lock(&self.best);
        g.calls += 1;
        if let Some(s) = solved
            && s.status == "Optimal"
            && s.x.len() == m.num_variables
        {
            let mut x = s.x.clone();
            for &j in &self.ints {
                x[j] = x[j].round();
            }
            let chk = check_conic_solution(m, &x, DEFAULT_TOL, 0, "", "");
            if chk.verdict == "feasible" {
                g.solved += 1;
                let obj = chk.objective.value;
                let maximize = m.obj_sense == IoSense::Maximize;
                if g.obj.is_none_or(|b| better(maximize, obj, b)) {
                    g.obj = Some(obj);
                    g.x = x.clone();
                    g.improved += 1;
                }
                let mut z = self.form.lift(&x);
                if let Some(e) = self.obj_eta {
                    let q: f64 = 0.5
                        * m.hessian
                            .iter()
                            .map(|&(i, j, v)| v * x[i] * x[j])
                            .sum::<f64>();
                    z[e] = if maximize { -q } else { q };
                }
                out = Some(z);
            }
        }
        g.seconds += t0.elapsed().as_secs_f64();
        out
    }
}

/// The parsed command line of a mixed-integer conic model.
struct Flags {
    json: bool,
    include_solution: bool,
    time_limit: Option<f64>,
    solution_file: Option<String>,
    opts: Vec<(String, String)>,
    settings: Settings,
}

/// `misocp_oa*` and other `key=value` options from `--opt` flags and `RUSOLVE_EKK_OPTS`
/// (the environment first, so a flag wins), in order.
pub fn opt_entries(flags: &[String]) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    if let Ok(list) = std::env::var("RUSOLVE_EKK_OPTS") {
        for kv in list.split(',').filter(|s| !s.is_empty()) {
            if let Some((k, v)) = kv.split_once('=') {
                out.push((k.to_string(), v.to_string()));
            }
        }
    }
    let mut it = flags.iter();
    while let Some(a) = it.next() {
        if a == "--opt"
            && let Some(kv) = it.next()
            && let Some((k, v)) = kv.split_once('=')
        {
            out.push((k.to_string(), v.to_string()));
        }
    }
    out
}

/// Whether the lever `misocp_oa` is on in `flags` / the environment (the last setting wins).
pub fn lever_on(flags: &[String]) -> bool {
    opt_entries(flags)
        .iter()
        .rfind(|(k, _)| k == "misocp_oa")
        .is_some_and(|(_, v)| matches!(v.as_str(), "true" | "1"))
}

fn parse_flags(cmd: &str, flags: &[String]) -> Result<Flags, String> {
    let mut f = Flags {
        json: false,
        include_solution: false,
        time_limit: None,
        solution_file: None,
        opts: Vec::new(),
        settings: Settings::default(),
    };
    let mut it = flags.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--format" => match it.next().map(String::as_str) {
                Some("json") => f.json = true,
                Some("text") => {}
                _ => return Err("--format needs json or text".into()),
            },
            "--solution" | "--include-solution" => f.include_solution = true,
            "--time-limit" => match it.next().and_then(|v| v.parse::<f64>().ok()) {
                Some(t) if t >= 0.0 && t.is_finite() => f.time_limit = Some(t),
                _ => return Err("--time-limit: want a non-negative number of seconds".into()),
            },
            "--solution-file" => match it.next() {
                Some(p) if p.ends_with(".json") => f.solution_file = Some(p.clone()),
                Some(_) => {
                    return Err(
                        "--solution-file: only a .json path is supported for this model class"
                            .into(),
                    );
                }
                None => return Err("--solution-file needs a path".into()),
            },
            "--opt" => {
                if it.next().is_none() {
                    return Err("--opt needs key=value".into());
                }
            }
            "--no-warnings" | "--no-presolve" => {}
            "--threads" | "--seed" => {
                if it.next().is_none() {
                    return Err(format!("{a} needs a value"));
                }
            }
            other => {
                return Err(format!(
                    "{cmd}: argument '{other}' is not supported for a mixed-integer conic model (supported: --format json|text, --time-limit, --include-solution, --solution-file X.json, --opt key=value)"
                ));
            }
        }
    }
    f.opts = opt_entries(flags);
    let on = |key: &str| {
        f.opts
            .iter()
            .rfind(|(k, _)| k == key)
            .is_some_and(|(_, v)| matches!(v.as_str(), "true" | "1"))
    };
    f.settings = Settings {
        disaggregate: on("misocp_oa_disaggregate"),
        no_polish: on("misocp_oa_no_polish"),
        lazy_only: on("misocp_oa_lazy_only"),
    };
    if f.include_solution && !f.json {
        return Err("--solution needs --format json".into());
    }
    Ok(f)
}

/// The outcome of a solve in the model's own sense.
struct Outcome {
    status: &'static str,
    x: Vec<f64>,
    objective: f64,
    dual_bound: f64,
    nodes: i64,
    detail: Value,
}

/// Solves `model` by outer approximation. `opts` are other registry options for the port.
fn solve(
    model: &Model,
    s: &Settings,
    opts: &[(String, String)],
    time_limit: Option<f64>,
) -> Result<Outcome, ConicError> {
    let t_start = Instant::now();
    let deadline = time_limit.map(|t| t_start + Duration::from_secs_f64(t));
    let left =
        |d: Option<Instant>| d.map(|d| d.saturating_duration_since(Instant::now()).as_secs_f64());
    let maximize = model.obj_sense == IoSense::Maximize;
    let n = model.num_variables;
    let built = build(model, s.disaggregate)?;
    let form = built.form;
    let ints: Vec<usize> = (0..n).filter(|&j| model.bounds[j].integer).collect();

    let mut detail = json!({
        "engine": "oa",
        "cones": built.cones,
        "oa_cones": form.cones.len(),
        "disaggregated": s.disaggregate,
        "extra_columns": form.extra.len(),
        "tree_separation": !s.lazy_only,
        "polish_enabled": !s.no_polish,
    });
    let fail = |status: &'static str, mut detail: Value, why: &str| {
        detail["note"] = json!(why);
        Ok(Outcome {
            status,
            x: Vec::new(),
            objective: f64::NAN,
            dual_bound: f64::NAN,
            nodes: 0,
            detail,
        })
    };

    // 1. the continuous SOCP relaxation
    let mut relaxed = model.clone();
    for b in &mut relaxed.bounds {
        b.integer = false;
    }
    let root = solve_model(&relaxed, left(deadline))?;
    detail["root_socp"] = json!({
        "status": root.status,
        "clarabel_status": root.clarabel_status,
        "objective_value": if root.status == "Optimal" { Some(root.objective) } else { None },
        "seconds": root.time_s,
    });
    match root.status {
        "Infeasible" => {
            return fail(
                "Infeasible",
                detail,
                "the continuous relaxation is infeasible (verified certificate)",
            );
        }
        "Unbounded" => return fail("Unknown", detail, "the continuous relaxation is unbounded"),
        _ => {}
    }
    let root_x = (root.status == "Optimal" && root.x.len() == n).then(|| root.x.clone());

    // 2. the model for the port
    let lp = build_lp(model, &form);
    let glo = Arc::new(lp.col_lower_.clone());
    let ghi = Arc::new(lp.col_upper_.clone());
    let mut options = crate::base_options();
    options.output_flag = false;
    options.threads = 1;
    options.time_limit = left(deadline).unwrap_or(kHighsInf);
    // The cuts stop at a violation of 2 * mip_feasibility_tolerance, and near a cone's apex the
    // objective error is about the square root of that (1e-6 gives 7e-4 on the unit disc), which
    // keeps the gap above mip_rel_gap. A tighter value makes the bound close (`--opt` overrides).
    options.mip_feasibility_tolerance = 1e-9;
    for (k, _) in opts {
        if k.starts_with("misocp_oa") {
            continue;
        }
        if registry::option_type_name(k).is_none() {
            return Err(ConicError::Solver(format!("unknown option '{k}'")));
        }
    }
    options.sync_log_options();
    let mut highs = Highs::with_options(options);
    for (k, v) in opts {
        if k.starts_with("misocp_oa") {
            continue;
        }
        if highs.set_option_value_str(k, v) == HighsStatus::kError {
            return Err(ConicError::Solver(format!("invalid option value {k}={v}")));
        }
    }
    highs.pass_model(lp);
    let form = Arc::new(form);
    let stats: OaStatsHandle = Arc::new(Mutex::new(OaStats::default()));

    let best: BestHandle = Arc::new(Mutex::new(Best::default()));
    let ctx = Arc::new(Mutex::new(PolishCtx {
        model: Arc::new(model.clone()),
        form: Arc::clone(&form),
        obj_eta: built.obj_eta,
        ints: ints.clone(),
        best: Arc::clone(&best),
        deadline,
        seen: HashSet::new(),
    }));

    // root cuts: initial tangents, the relaxation's tangent, the disaggregation rows
    let root_z = root_x.as_ref().map(|x| {
        let mut z = form.lift(x);
        if let Some(e) = built.obj_eta {
            let q: f64 = 0.5
                * model
                    .hessian
                    .iter()
                    .map(|&(i, j, v)| v * x[i] * x[j])
                    .sum::<f64>();
            z[e] = if maximize { -q } else { q };
        }
        z
    });
    let root_rows = add_cone_rows(&mut highs, &form, root_z.as_deref(), &glo, &ghi, &stats);
    detail["root_rows"] = json!(root_rows);

    // the rounded relaxation is the first incumbent candidate
    if !s.no_polish
        && let Some(z) = &root_z
        && let Some(point) = lock(&ctx).polish(z)
    {
        let sol = HighsSolution {
            col_value: point,
            ..Default::default()
        };
        let _ = highs.set_solution(&sol);
    }
    if !s.lazy_only {
        // hooks run in the presolved column space: the cut indices are the model's own only
        // with presolve off
        highs.set_option_value_str("presolve", "off");
    }
    let polisher: Option<Polisher> = if s.no_polish {
        None
    } else {
        let shared = Arc::clone(&ctx);
        Some(Box::new(move |cand: &[f64]| lock(&shared).polish(cand)))
    };
    let out = run_oa_with(
        &mut highs,
        Arc::clone(&form),
        glo,
        ghi,
        Arc::clone(&stats),
        OaConfig {
            tree: !s.lazy_only,
            max_rounds: 2000,
            ..OaConfig::default()
        },
        deadline,
        polisher,
    );

    // 3. the answer: the best checker-feasible point and an honest status
    let hs = highs.get_model_status();
    let info_bound = highs.get_info().mip_dual_bound;
    let mut cand: Vec<(f64, Vec<f64>, &'static str)> = Vec::new();
    let mut accepted_checker = None;
    if highs.solution_.value_valid && highs.solution_.col_value.len() >= n {
        let x: Vec<f64> = highs.solution_.col_value[..n].to_vec();
        let mut xr = x.clone();
        for &j in &ints {
            xr[j] = xr[j].round();
        }
        let chk = check_conic_solution(model, &xr, DEFAULT_TOL, 0, "", "");
        accepted_checker = Some(chk.verdict);
        // The cuts leave the cones feasible only within tolerance, so the point's objective can
        // sit a little below the true one (a square-root effect near an apex). The continuous
        // SOCP at its integer part gives the exact point; it replaces the candidate.
        let exact = !s.no_polish && {
            let mut g = lock(&ctx);
            let key: Vec<i64> = ints.iter().map(|&j| xr[j] as i64).collect();
            g.seen.remove(&key);
            let z = form.lift(&xr);
            g.polish(&z).is_some()
        };
        if !exact && chk.verdict == "feasible" {
            cand.push((chk.objective.value, xr, "oa"));
        }
    }
    {
        let g = lock(&best);
        if let Some(o) = g.obj {
            cand.push((o, g.x.clone(), "polish"));
        }
    }
    cand.sort_by(|a, b| {
        if maximize {
            b.0.total_cmp(&a.0)
        } else {
            a.0.total_cmp(&b.0)
        }
    });
    let incumbent = cand.into_iter().next();
    let bound = if info_bound.is_finite() {
        info_bound
    } else if maximize {
        f64::INFINITY
    } else {
        f64::NEG_INFINITY
    };
    let rel_gap = highs.options_.mip_rel_gap;
    let mut status: &'static str = match hs {
        HighsModelStatus::kOptimal if out.lazy.accepted => "Optimal",
        HighsModelStatus::kOptimal => "IterationLimit",
        HighsModelStatus::kInfeasible => "Infeasible",
        HighsModelStatus::kTimeLimit => "TimeLimit",
        HighsModelStatus::kInterrupt => "Interrupted",
        HighsModelStatus::kIterationLimit => "IterationLimit",
        _ => "Unknown",
    };
    let mut bound_inconsistent = false;
    let mut objective = f64::NAN;
    let mut x = Vec::new();
    let mut source = Value::Null;
    if let Some((obj, xi, src)) = &incumbent {
        objective = *obj;
        x = xi.clone();
        source = json!(src);
        if bound.is_finite() {
            // A valid bound never passes the incumbent by more than the checker's own slack: the
            // checker accepts rows and cones violated by DEFAULT_TOL relative to their scale, and
            // a Clarabel point is exact only to its tolerance, so a checked point can sit a
            // little below the true optimum (clay0303m: 1.2e-6 relative). Ten times that
            // slack is an error of the bound or of the cuts.
            let tol = 10.0 * DEFAULT_TOL * (1.0 + obj.abs());
            let passed = if maximize {
                bound < *obj - tol
            } else {
                bound > *obj + tol
            };
            let gap = (obj - bound).abs();
            let closed = gap <= rel_gap.max(1e-9) * obj.abs().max(1e-10) || gap <= 1e-9;
            if passed {
                // a bound that cuts off a checked feasible point is wrong: say so, claim nothing
                bound_inconsistent = true;
                status = "Unknown";
            } else if status == "Optimal" && !closed {
                status = "IterationLimit";
            } else if status == "Infeasible" {
                // a feasible point exists: the relaxation cannot be infeasible
                bound_inconsistent = true;
                status = "Unknown";
            }
        }
    } else if status == "Optimal" {
        // an accepted candidate the checker rejects proves nothing
        status = "Unknown";
    }
    if matches!(status, "TimeLimit" | "IterationLimit" | "Unknown") && incumbent.is_none() {
        x = Vec::new();
    }
    if hs == HighsModelStatus::kUnbounded || hs == HighsModelStatus::kUnboundedOrInfeasible {
        status = "Unknown";
    }
    let g = lock(&best);
    let st = &out.stats;
    detail["rounds"] = json!(out.lazy.rounds);
    detail["accepted"] = json!(out.lazy.accepted);
    detail["nodes_total"] = json!(out.lazy.mip_nodes);
    detail["cuts"] = json!({
        "root_tangent": st.root_tangents,
        "tree": st.tree_cuts,
        "tree_calls": st.tree_calls,
        "lazy_new": st.lazy_cuts_new,
        "lazy_from_tree": st.lazy_cuts_from_pool,
        "dropped_unsafe": st.cuts_unsafe,
        "rows_added_to_model": out.lazy.rows_added,
    });
    detail["polish"] = json!({
        "calls": g.calls,
        "feasible": g.solved,
        "improved": g.improved,
        "seconds": g.seconds,
    });
    detail["incumbent_source"] = source;
    detail["oa_solution_checker"] = json!(accepted_checker);
    detail["bound_inconsistent"] = json!(bound_inconsistent);
    detail["port_status"] = json!(format!("{hs:?}"));
    drop(g);
    Ok(Outcome {
        status,
        x,
        objective,
        dual_bound: bound,
        nodes: out.lazy.mip_nodes,
        detail,
    })
}

/// `rusolve solve|mip <model> [flags]` for a model with integer columns and cones or quadratic
/// rows when `misocp_oa` is on. `args` is what followed the command name.
pub fn run(cmd: &str, args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let Some(path) = args.first() else {
        return crate::diag::fail_msg(
            wj,
            "RS0006",
            2,
            format!("rusolve {cmd}: missing model path"),
        );
    };
    let flags = match parse_flags(cmd, &args[1..]) {
        Ok(f) => f,
        Err(e) => return crate::diag::fail_msg(wj, "RS0003", 2, format!("rusolve {cmd}: {e}")),
    };
    for (k, _) in &flags.opts {
        if !k.starts_with("misocp_oa") && registry::option_type_name(k).is_none() {
            return crate::diag::unknown_option(flags.json, cmd, k);
        }
    }
    let model = match rusolve_io::read_model_file_conic(path) {
        Ok(m) => m,
        Err(e) => return crate::diag::fail_read(flags.json, path, &e.to_string()),
    };
    let t0 = Instant::now();
    let solved = match solve(&model, &flags.settings, &flags.opts, flags.time_limit) {
        Ok(s) => s,
        Err(e) => {
            return crate::diag::fail_msg(flags.json, "RS0012", 3, format!("rusolve {cmd}: {e}"));
        }
    };
    let time_s = t0.elapsed().as_secs_f64();
    if let Some(p) = &flags.solution_file {
        let cols: Vec<Value> = if solved.x.is_empty() {
            Vec::new()
        } else {
            (0..model.num_variables)
                .map(|j| json!({"name": col_name(&model, j), "value": solved.x[j]}))
                .collect()
        };
        if let Err(e) = std::fs::write(
            p,
            serde_json::to_string_pretty(&json!({"columns": cols})).unwrap_or_default(),
        ) {
            return crate::diag::fail_msg(
                false,
                "RS0090",
                4,
                format!("rusolve {cmd}: cannot write {p}: {e}"),
            );
        }
    }
    let mut detail = solved.detail.clone();
    if !solved.x.is_empty() {
        let chk = check_conic_solution(&model, &solved.x, DEFAULT_TOL, 0, "", "");
        detail["checker"] = json!({
            "verdict": chk.verdict,
            "tolerance": chk.tolerance,
            "violations": chk.summary.total_violations,
        });
    }
    if flags.json {
        let nsol = NativeSolved {
            path,
            model: &model,
            status: solved.status,
            objective: solved.objective,
            dual_bound: solved.dual_bound,
            nodes: solved.nodes,
            lp_iterations: 0,
            time_s,
            seed: 0,
            threads: 1,
            time_limit: flags.time_limit,
            node_limit: None,
            x: &solved.x,
            include_solution: flags.include_solution,
        };
        let mut v = serde_json::to_value(native_envelope(&nsol)).expect("envelope serialises");
        v["solver"]["engine"] = json!("highs-port");
        v["misocp_oa"] = detail;
        println!("{v}");
    } else {
        println!("status: {}", solved.status);
        if !solved.x.is_empty() {
            println!("objective: {}", solved.objective);
        }
        if solved.dual_bound.is_finite() {
            println!("dual bound: {}", solved.dual_bound);
        }
        println!(
            "misocp_oa: {} rounds, {} nodes ({:.3}s)",
            detail["rounds"], solved.nodes, time_s
        );
        if let Some(cn) = detail.get("checker") {
            println!(
                "checker: {} ({} violations)",
                cn["verdict"], cn["violations"]
            );
        }
    }
    ExitCode::SUCCESS
}
