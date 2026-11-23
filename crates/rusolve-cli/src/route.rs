//! Problem-class router: one classifier and one routing table for
//! `rusolve solve` and `rusolve mip`, so the two commands agree on what a model is and who
//! solves it.
//!
//! The model is read once, here, and only when a byte scan ([`rusolve_io::may_have_special_content`])
//! says it may hold a quadratic term, a quadratic row, a cone, a semi-continuous column or a
//! general constraint. A plain LP or MIP is never read here: it falls straight through to
//! `run_legacy`, which is the path it always took (`tools/bench/identity.sh`).
//!
//! | class | engine | `solve` | `mip` |
//! |---|---|---|---|
//! | `lp`, `mip` | the LP / MIP path | as before | as before |
//! | `convex_qp` | `hipo` or `quass` (as the `qp` block) | QP path | QP path |
//! | `gc`, `gc_mip` | `gc-linearization` (`gc_solve.rs` + HiGHS port) | gc route | gc route |
//! | `conic` | `clarabel` (`conic_solve.rs`) | conic route | conic route |
//! | `nonconvex_qp`, `miqp`, `nonconvex_qcqp`, `mi_conic`, `gc_quadratic` | none | typed error, exit 3 | typed error, exit 3 |
//!
//! `rusolve classify --format json` reports the same class name.

use rusolve_io::Model;
use serde_json::{Value, json};
use std::process::ExitCode;
use std::sync::OnceLock;

/// What a model is, by the objects it carries.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    Lp,
    Mip,
    ConvexQp,
    NonconvexQp,
    Miqp,
    /// Continuous convex QCQP or SOCP (quadratic rows and cones, with or without a Q objective).
    Conic,
    NonconvexQcqp,
    MiConic,
    /// Indicator / SOS / PWL / min / max / abs / and / or / semi-continuous, no integers.
    Gc,
    GcMip,
    /// A general constraint together with a quadratic objective, row or cone.
    GcQuadratic,
}

impl Class {
    pub fn name(self) -> &'static str {
        match self {
            Class::Lp => "lp",
            Class::Mip => "mip",
            Class::ConvexQp => "convex_qp",
            Class::NonconvexQp => "nonconvex_qp",
            Class::Miqp => "miqp",
            Class::Conic => "conic",
            Class::NonconvexQcqp => "nonconvex_qcqp",
            Class::MiConic => "mi_conic",
            Class::Gc => "gc",
            Class::GcMip => "gc_mip",
            Class::GcQuadratic => "gc_quadratic",
        }
    }

    /// The engine that solves the class; `None` for a class nothing solves yet.
    pub fn engine(self) -> Option<&'static str> {
        match self {
            Class::Lp | Class::Mip => Some("highs-port"),
            Class::ConvexQp => Some("hipo"),
            Class::Conic => Some("clarabel"),
            Class::Gc | Class::GcMip => Some("gc-linearization"),
            _ => None,
        }
    }

    /// `convex` / `nonconvex` for a class with a quadratic object, else `None`.
    pub fn convexity(self) -> Option<&'static str> {
        match self {
            Class::ConvexQp | Class::Conic => Some("convex"),
            Class::NonconvexQp | Class::NonconvexQcqp => Some("nonconvex"),
            _ => None,
        }
    }
}

/// A classification: the class and, for a class nothing solves, why.
pub struct Verdict {
    pub class: Class,
    pub reason: String,
    /// The QP front's own typed error (`error` object of `rusolve.error/1`: eigenvalue, variables),
    /// kept for a nonconvex Q and a MIQP so the refusal says as much as it did before the router.
    pub detail: Option<Value>,
}

/// The class of a model that has been read. `Lp` / `Mip` mean "the legacy paths".
pub fn classify_model(m: &Model) -> Verdict {
    let q_obj = !m.hessian.is_empty();
    let conic = !m.conic.is_empty();
    let gc = m.gc.has_general() || !m.gc.semi.is_empty();
    let ints = m.bounds.iter().filter(|b| b.integer).count();
    let v = |class, reason: String| Verdict {
        class,
        reason,
        detail: None,
    };
    if gc && (q_obj || conic) {
        let what =
            m.gc.first_general()
                .unwrap_or_else(|| "semi-continuous columns".into());
        return v(
            Class::GcQuadratic,
            format!(
                "the model has {what} together with quadratic terms (objective, rows or cones); \
                 no engine handles general constraints on a quadratic model"
            ),
        );
    }
    if conic {
        if let Some(j) = m.bounds.iter().position(|b| b.integer) {
            let name = m.col_names.get(j).cloned().unwrap_or_default();
            return v(
                Class::MiConic,
                format!(
                    "the model has quadratic rows or cones and {ints} integer variable(s) (column {name} is integer); \
                     mixed-integer conic models are not supported yet"
                ),
            );
        }
        return match crate::conic_solve::build_form(m) {
            Err(
                e @ (crate::conic_solve::ConicError::NonconvexRow { .. }
                | crate::conic_solve::ConicError::NonconvexObjective { .. }
                | crate::conic_solve::ConicError::EqualityQuadRow { .. }),
            ) => v(Class::NonconvexQcqp, e.to_string()),
            // anything else (an unsupported cone, a size limit) is the conic route's own error
            _ => v(Class::Conic, String::new()),
        };
    }
    if q_obj {
        return match crate::qp::assess(m) {
            Ok(_) => v(Class::ConvexQp, String::new()),
            Err(r) => {
                // a Q with integers is a MIQP first; without integers it is a nonconvex QP
                let class = if ints > 0 {
                    Class::Miqp
                } else {
                    Class::NonconvexQp
                };
                Verdict {
                    class,
                    reason: r.message,
                    detail: r.json.get("error").cloned(),
                }
            }
        };
    }
    if gc {
        return v(
            if ints > 0 { Class::GcMip } else { Class::Gc },
            String::new(),
        );
    }
    v(if ints > 0 { Class::Mip } else { Class::Lp }, String::new())
}

/// The class and engine of the run in progress, for the JSON envelope.
static ROUTE: OnceLock<(Class, Option<&'static str>)> = OnceLock::new();

pub fn note(class: Class) {
    let _ = ROUTE.set((class, class.engine()));
}

/// `(class, engine, convexity)` for the envelope. Without a router verdict (a plain LP or MIP,
/// which is never read by the router) the class comes from the integer count and the engine is
/// the solver stack that ran it.
pub fn envelope_fields(
    integer_cols: i64,
    stack_engine: &str,
) -> (String, String, Option<&'static str>) {
    match ROUTE.get() {
        Some((c, e)) => (
            c.name().to_string(),
            if *c == Class::ConvexQp {
                // the engine that produced the answer, as the `qp` block names it (hipo / quass)
                crate::qp::engine_name().to_string()
            } else {
                e.map_or_else(|| stack_engine.to_string(), str::to_string)
            },
            c.convexity(),
        ),
        None => (
            if integer_cols > 0 { "mip" } else { "lp" }.to_string(),
            stack_engine.to_string(),
            None,
        ),
    }
}

fn wants_json(args: &[String]) -> bool {
    args.windows(2)
        .any(|w| w[0] == "--format" && w[1] == "json")
}

/// The typed error of a class nothing solves: stderr text, and under `--format json` a
/// `rusolve.error/1` object on stdout. Exit code 3 (input the solver cannot honour).
pub fn refuse(cmd: &str, args: &[String], verdict: &Verdict) -> ExitCode {
    let c = verdict.class;
    let message = format!("{} models are not supported: {}", c.name(), verdict.reason);
    eprintln!("rusolve {cmd}: {message}");
    if wants_json(args) {
        let mut err = verdict.detail.clone().unwrap_or_else(|| json!({}));
        if let Some(q) = err.get("class").cloned() {
            err["qplib_class"] = q;
        }
        if err.get("kind").is_none() {
            err["kind"] = json!("unsupported_class");
        }
        err["message"] = json!(message);
        err["class"] = json!(c.name());
        err["reason"] = json!(verdict.reason);
        if let Some(cv) = c.convexity() {
            err["convexity"] = json!(cv);
        }
        // The same registered diagnostic the text path would carry, additively
        let d = rusolve_core::diagnostic::Diagnostic::error(
            "RS0012",
            format!("rusolve {cmd}: {message}"),
        )
        .with_help("`rusolve --explain RS0012`");
        let diags: Value =
            serde_json::from_str(&d.to_json()).unwrap_or_else(|_| json!({"code": "RS0012"}));
        println!(
            "{}",
            json!({"schema_version": "rusolve.error/1", "error": err, "diagnostics": [diags]})
        );
    }
    ExitCode::from(3)
}

/// `solve` flags to the `mip` flag set, for a class `solve` hands to the `mip` engine. `--ekk`
/// is the old no-op spelling and is dropped; an engine choice (`--legacy-lp`, `--ipm`) is an error.
fn solve_to_mip_flags(flags: &[String]) -> Result<Vec<String>, String> {
    let mut out = Vec::new();
    let mut it = flags.iter();
    while let Some(f) = it.next() {
        match f.as_str() {
            "--ekk" => {}
            "--csv" | "--solution" | "--include-solution" | "--no-presolve" | "--no-warnings" => {
                out.push(f.clone())
            }
            "--format" | "--time-limit" | "--solution-file" => {
                let v = it.next().ok_or_else(|| format!("{f} needs a value"))?;
                out.push(f.clone());
                out.push(v.clone());
            }
            other => {
                return Err(format!(
                    "argument '{other}' is not supported for this model class \
                     (supported: --format, --csv, --solution, --time-limit, --solution-file, --no-presolve, --no-warnings)"
                ));
            }
        }
    }
    Ok(out)
}

/// `mip` flags to the `solve` flag set, for a class `mip` hands to the QP engine. `--seed` and
/// `--threads` are accepted and ignored, as on the conic route.
fn mip_to_solve_flags(flags: &[String]) -> Result<Vec<String>, String> {
    let mut out = Vec::new();
    let mut it = flags.iter();
    while let Some(f) = it.next() {
        match f.as_str() {
            "--no-warnings" | "--csv" | "--solution" | "--include-solution" | "--no-presolve" => {
                out.push(f.clone())
            }
            "--format" | "--time-limit" | "--solution-file" => {
                let v = it.next().ok_or_else(|| format!("{f} needs a value"))?;
                out.push(f.clone());
                out.push(v.clone());
            }
            "--seed" | "--threads" => {
                it.next().ok_or_else(|| format!("{f} needs a value"))?;
            }
            other => {
                return Err(format!(
                    "argument '{other}' is not supported for this model class \
                     (supported: --format, --csv, --solution, --time-limit, --solution-file, --no-presolve, --no-warnings, --seed, --threads)"
                ));
            }
        }
    }
    Ok(out)
}

/// The router. `Some(code)` when it ran the command; `None` when the model is a plain LP or MIP
/// (or unreadable, or the class keeps its legacy path) and `run_legacy` must run untouched.
pub fn dispatch(cmd: &str, path: &str, args: &[String]) -> Option<ExitCode> {
    if !rusolve_io::may_have_special_content(path) {
        return None;
    }
    let model = rusolve_io::read_model_file_conic(path).ok()?;
    let verdict = classify_model(&model);
    let class = verdict.class;
    if matches!(class, Class::Lp | Class::Mip) {
        return None;
    }
    note(class);
    let rest = &args[1..];
    let usage = |e: String| {
        Some(crate::diag::fail_msg(
            crate::diag::wants_json(args),
            "RS0003",
            2,
            format!("rusolve {cmd}: {e}"),
        ))
    };
    match class {
        Class::Conic => Some(crate::conic_solve::run(cmd, args)),
        Class::ConvexQp => {
            if cmd == "solve" {
                // the QP path of `solve` is the legacy one (ekk.rs)
                return None;
            }
            match mip_to_solve_flags(rest) {
                Ok(mut f) => {
                    let mut a = vec!["solve".to_string(), path.to_string()];
                    a.append(&mut f);
                    Some(crate::run_legacy(&a))
                }
                Err(e) => usage(e),
            }
        }
        Class::Gc | Class::GcMip => {
            if cmd == "mip" {
                return None;
            }
            match solve_to_mip_flags(rest) {
                Ok(mut f) => {
                    let mut a = vec![path.to_string()];
                    a.append(&mut f);
                    Some(crate::mip_highs::run_mip_highs(&a))
                }
                Err(e) => usage(e),
            }
        }
        _ => Some(refuse(cmd, args, &verdict)),
    }
}

/// `class` for `rusolve classify --format json`, from the same classifier.
pub fn class_of_file(path: &str) -> Option<&'static str> {
    let m = rusolve_io::read_model_file_conic(path).ok()?;
    Some(classify_model(&m).class.name())
}
