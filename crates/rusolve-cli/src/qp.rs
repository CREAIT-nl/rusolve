//! Convex QP front: the convexity check in front of every solve that carries a Hessian.
//!
//! A model with no Hessian never reaches this file's solve-side code: `rusolve solve` calls
//! [`assess`] only when `Model::hessian` is non-empty, so LP and MIP take exactly their old path.
//!
//! What it does for a QP:
//! * [`assess`] runs `rusolve_core::qp_convexity::analyze` on the objective's Q (PSD for a
//!   minimisation, NSD for a maximisation), counts the variable types for the QPLIB-style class
//!   code, and refuses with a typed error (a nonconvex Q with the smallest eigenvalue estimate
//!   and the variables it lives on, or a Hessian together with integer variables, which is
//!   a MIQP) before the engine starts;
//! * a Q that is indefinite by rounding only gets `Q + shift*I` (`-Q + shift*I` for a
//!   maximisation, i.e. `Q - shift*I`), with the shift reported;
//! * [`qp_json`] is the `qp` object of the `rusolve.result/1` envelope: class code, convexity
//!   verdict, method, shift and the engine that ran.
//!
//! The engine in force is HiGHS's choice for a QP: QUASS (the active-set code) unless
//! `solver=ipm|hipo` with `RUSOLVE_EKK_HIPO=1` selects HiPO (`lp_data/highs_qp.rs:60`).

use rusolve_core::qp_convexity::{Convexity, Method, Verdict, analyze, class_code};
use rusolve_ekk::prelude::{HessianFormat, HighsHessian, HighsInt};
use rusolve_io::{Model, ObjSense};
use serde_json::{Value, json};

/// The verdict for a model with a Hessian.
pub struct QpAssessment {
    pub convexity: Convexity,
    /// QPLIB-style three-letter class (Furini et al. 2019).
    pub class: String,
    pub maximize: bool,
}

/// A typed refusal: `message` for stderr, `json` for `--format json`.
pub struct QpRefusal {
    pub message: String,
    pub json: Value,
}

fn var_name(m: &Model, j: usize) -> String {
    match m.col_names.get(j) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("c{j}"),
    }
}

/// Count `(continuous, binary, other integer)` columns and whether any bound is finite.
fn type_counts(m: &Model) -> (usize, usize, usize, bool) {
    let (mut c, mut b, mut i, mut fin) = (0, 0, 0, false);
    for v in &m.bounds {
        fin |= v.lower.is_finite() || v.upper.is_finite();
        if !v.integer {
            c += 1;
        } else if v.lower >= 0.0 && v.upper <= 1.0 {
            b += 1;
        } else {
            i += 1;
        }
    }
    (c, b, i, fin)
}

/// The convexity verdict of a model with integer columns and a Hessian, for the outer-
/// approximation route: `Ok` only when Q is positive semidefinite (negative semidefinite for a
/// maximisation) with no rounding shift, because the tangent cuts are valid for the Q as written
/// and for no shifted one. The message says why not otherwise.
pub fn assess_miqp(model: &Model) -> Result<QpAssessment, String> {
    let maximize = model.obj_sense == ObjSense::Maximize;
    let (cont, bin, other, finite_bound) = type_counts(model);
    let cvx = analyze(model.num_variables, &model.hessian, maximize);
    let class = class_code(&cvx, cont, bin, other, model.num_constraints, finite_bound);
    match cvx.verdict {
        Verdict::Convex => Ok(QpAssessment {
            convexity: cvx,
            class,
            maximize,
        }),
        Verdict::ConvexWithShift => Err(format!(
            "Q is convex only after a rounding shift of {:e}; outer approximation needs the exact convex Q",
            cvx.shift
        )),
        Verdict::Nonconvex => Err(
            "Q is not positive semidefinite (negative semidefinite for a maximisation)".to_string(),
        ),
        Verdict::Unverified => {
            Err("convexity of Q could not be verified within the work cap".to_string())
        }
    }
}

/// Test the objective Hessian of `model` (which must have one). `Err` for a nonconvex Q or a
/// QP with integer variables.
pub fn assess(model: &Model) -> Result<QpAssessment, QpRefusal> {
    let maximize = model.obj_sense == ObjSense::Maximize;
    let (cont, bin, other, finite_bound) = type_counts(model);
    let cvx = analyze(model.num_variables, &model.hessian, maximize);
    let class = class_code(&cvx, cont, bin, other, model.num_constraints, finite_bound);
    if cvx.verdict == Verdict::Nonconvex {
        let vars: Vec<String> = cvx
            .witness
            .iter()
            .map(|&(j, _)| var_name(model, j))
            .collect();
        let lam = cvx.min_eig_estimate.unwrap_or(f64::NAN);
        // For a maximisation the form tested is -Q: report on Q itself.
        let (kind, shown, word) = if maximize {
            (
                "largest",
                -lam,
                "maximisation needs Q negative semidefinite",
            )
        } else {
            (
                "smallest",
                lam,
                "minimisation needs Q positive semidefinite",
            )
        };
        let message = format!(
            "nonconvex quadratic objective: Q has {kind} eigenvalue {shown:.6e} (estimate by {}), \
             on variables {}; a {word}; nonconvex QP is not supported yet (class {class})",
            match cvx.method {
                Method::Diagonal => "a principal minor",
                _ => "Lanczos",
            },
            vars.join(", ")
        );
        let json = json!({
            "schema_version": "rusolve.error/1",
            "error": {
                "kind": "nonconvex_quadratic_objective",
                "message": message,
                "class": class,
                "eigenvalue": shown,
                "eigenvalue_kind": kind,
                "variables": vars,
            }
        });
        return Err(QpRefusal { message, json });
    }
    if cont + bin + other > 0 && bin + other > 0 {
        let message = format!(
            "quadratic objective with {} integer variable(s): MIQP is not supported yet (class {class})",
            bin + other
        );
        let json = json!({
            "schema_version": "rusolve.error/1",
            "error": {"kind": "miqp_unsupported", "message": message, "class": class}
        });
        return Err(QpRefusal { message, json });
    }
    Ok(QpAssessment {
        convexity: cvx,
        class,
        maximize,
    })
}

/// The Hessian handed to the port: unchanged unless the check asked for a shift.
pub fn solver_hessian(model: &Model, a: &QpAssessment) -> HighsHessian {
    if a.convexity.shift == 0.0 {
        return crate::ekk::highs_hessian_from_model(model);
    }
    let n = model.num_variables;
    let sh = if a.maximize {
        -a.convexity.shift
    } else {
        a.convexity.shift
    };
    let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut has_diag = vec![false; n];
    for &(i, j, v) in &model.hessian {
        if i == j {
            has_diag[j] = true;
            cols[j].push((i, v + sh));
        } else {
            cols[j].push((i, v));
        }
    }
    for j in 0..n {
        if !has_diag[j] {
            cols[j].push((j, sh));
        }
    }
    let mut h = HighsHessian {
        dim_: n as HighsInt,
        format_: HessianFormat::kSquare,
        start_: vec![0],
        index_: Vec::new(),
        value_: Vec::new(),
    };
    for c in &mut cols {
        c.sort_by_key(|&(i, _)| i);
        for &(i, v) in c.iter() {
            h.index_.push(i as HighsInt);
            h.value_.push(v);
        }
        h.start_.push(h.index_.len() as HighsInt);
    }
    h
}

/// The `qp` object of the envelope.
pub fn qp_json(a: &QpAssessment, route: &QpRoute) -> Value {
    let c = &a.convexity;
    json!({
        "class": a.class,
        "convexity": match c.verdict {
            Verdict::Convex => "convex",
            Verdict::ConvexWithShift => "convex_after_shift",
            Verdict::Nonconvex => "nonconvex",
            Verdict::Unverified => "unverified",
        },
        "method": match c.method {
            Method::Empty => "empty",
            Method::Diagonal => "diagonal",
            Method::Gershgorin => "gershgorin",
            Method::Ldlt => "ldlt",
            Method::Lanczos => "lanczos",
        },
        "shift": c.shift,
        "min_eigenvalue_estimate": c.min_eig_estimate,
        "hessian_nonzeros": c.nnz,
        "engine": engine_name(),
        "engine_selection": if route.explicit { "user" } else { "default" },
        "fallback": route.fallback_from.is_some(),
        "fallback_from_status": route.fallback_from,
    })
}

/// How the engine was chosen for this solve (HiPO first, QUASS on a non-Optimal HiPO).
pub struct QpRoute {
    /// The user named `solver`; no automatic choice, no fallback.
    pub explicit: bool,
    /// The status HiPO ended with when QUASS was run instead.
    pub fallback_from: Option<String>,
}

impl QpRoute {
    pub fn new(explicit: bool) -> Self {
        Self {
            explicit,
            fallback_from: None,
        }
    }
}

/// The one text-mode line that says what ran (stdout, not stderr).
pub fn text_line(a: &QpAssessment, route: &QpRoute) -> String {
    let c = &a.convexity;
    let how = match (&route.fallback_from, route.explicit) {
        (Some(s), _) => format!(" (HiPO ended {s}, fell back)"),
        (None, true) => " (solver chosen by the user)".to_string(),
        (None, false) => String::new(),
    };
    format!(
        "qp: class {}, {:?} ({:?}), shift {:e}, engine {}{how}",
        a.class,
        c.verdict,
        c.method,
        c.shift,
        engine_name()
    )
}

/// The engine the last QP solve ran.
pub fn engine_name() -> &'static str {
    // `lp_solver_auto` can pick HiPO for a QP from the shape of A alone, so the options after
    // the run are not evidence; the QP dispatch records what it ran.
    if rusolve_ekk::lp_data::highs_qp::LAST_QP_ENGINE_HIPO
        .load(std::sync::atomic::Ordering::Relaxed)
    {
        "hipo"
    } else {
        "quass"
    }
}

/// The message of a command that has no QP support and would otherwise drop the Hessian.
pub fn unsupported(command: &str) -> String {
    format!(
        "the model has a quadratic objective, which `rusolve {command}` ignores; \
         solve it with `rusolve solve` (convex QP)"
    )
}

/// `envelope` (a `rusolve.result/1` JSON text) with the `qp` object added.
pub fn with_qp(envelope: &str, a: &QpAssessment, route: &QpRoute) -> String {
    match serde_json::from_str::<Value>(envelope) {
        Ok(Value::Object(mut o)) => {
            o.insert("qp".to_string(), qp_json(a, route));
            Value::Object(o).to_string()
        }
        _ => envelope.to_string(),
    }
}
