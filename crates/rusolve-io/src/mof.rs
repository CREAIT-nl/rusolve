//! MathOptFormat (`.mof.json`) reader, linear and mixed-integer subset.
//!
//! Source: the MathOptFormat specification and JSON schema,
//! github.com/jump-dev/MathOptFormat (`README.md`,
//! `schemas/mof.1.9.schema.json`; Legat, Dowson, Garcia, Lubin,
//! "MathOptInterface: a data structure for mathematical optimization
//! problems", INFORMS J. Computing 34(2), 2022, §5). A model is
//! `{"version": {"major": 1, "minor": m}, "variables": [{"name"}],
//! "objective": {"sense", "function"}, "constraints": [{"name"?,
//! "function", "set"}]}`; variables are free unless a constraint bounds
//! them; duplicate terms in a function are summed (schema, ScalarAffineFunction).
//!
//! Supported: schema major version 1 (all minors to date); objective
//! `min` / `max` / `feasibility` with a `Variable` or `ScalarAffineFunction`
//! (its `constant` becomes the objective offset); `Variable`-in-set
//! constraints read as variable bounds (`LessThan`, `GreaterThan`,
//! `EqualTo`, `Interval`, several bounds on one variable intersect) or
//! integrality (`Integer`, `ZeroOne`, the latter also bounding to [0, 1]);
//! `ScalarAffineFunction`-in-`LessThan` / `GreaterThan` / `EqualTo` /
//! `Interval` constraints read as rows, the function's `constant` moved
//! into the bounds.
//!
//! Rejected with an error naming the feature: other functions
//! (quadratic, nonlinear, vector), other sets (semicontinuous,
//! semi-integer, parameters, cones, SOS, complementarity, indicator, ...),
//! an affine function in `Integer` / `ZeroOne`, vector objectives (multi-
//! objective), unknown or duplicate variable names, major versions other
//! than 1, and a bound or coefficient that is not a JSON number (the format
//! has no infinity literal; an absent bound is how a side is left open).
//! Optional fields the solver does not use (`primal_start`, `dual_start`,
//! `description`, `author`) are ignored.

use std::collections::HashMap;

use crate::model::{Model, ObjSense, RowSense, VariableBound};
use rusolve_core::{CscMatrix, Triplet};
use serde_json::Value;

/// A MathOptFormat read failure: the JSON location (e.g.
/// `constraints[3].set`) and why.
#[derive(Debug, Clone, PartialEq)]
pub struct MofError {
    pub location: String,
    pub message: String,
}

impl std::fmt::Display for MofError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "MathOptFormat read error at {}: {}",
            self.location, self.message
        )
    }
}

impl std::error::Error for MofError {}

fn err<T>(location: impl Into<String>, message: impl Into<String>) -> Result<T, MofError> {
    Err(MofError {
        location: location.into(),
        message: message.into(),
    })
}

fn unsupported(what: &str) -> String {
    format!("{what} not supported (rusolve reads linear and mixed-integer linear MathOptFormat)")
}

fn field<'a>(v: &'a Value, key: &str, at: &str) -> Result<&'a Value, MofError> {
    match v.get(key) {
        Some(x) => Ok(x),
        None => err(at, format!("missing field \"{key}\"")),
    }
}

fn string<'a>(v: &'a Value, key: &str, at: &str) -> Result<&'a str, MofError> {
    match field(v, key, at)?.as_str() {
        Some(s) => Ok(s),
        None => err(format!("{at}.{key}"), "expected a string"),
    }
}

fn number(v: &Value, key: &str, at: &str) -> Result<f64, MofError> {
    match field(v, key, at)? {
        Value::Number(n) => match n.as_f64() {
            Some(x) => Ok(x),
            None => err(format!("{at}.{key}"), "number out of range"),
        },
        other => err(
            format!("{at}.{key}"),
            format!("expected a JSON number, found {other}"),
        ),
    }
}

fn array<'a>(v: &'a Value, key: &str, at: &str) -> Result<&'a [Value], MofError> {
    match field(v, key, at)?.as_array() {
        Some(a) => Ok(a),
        None => err(format!("{at}.{key}"), "expected an array"),
    }
}

/// A scalar affine function `Σ coef·x_j + constant`.
struct Affine {
    terms: Vec<(usize, f64)>,
    constant: f64,
    /// Written as `{"type": "Variable"}` (a bound, when in a constraint).
    single_variable: Option<usize>,
}

fn function(f: &Value, at: &str, index: &HashMap<&str, usize>) -> Result<Affine, MofError> {
    let var = |name: &str, at: String| match index.get(name) {
        Some(&j) => Ok(j),
        None => err(at, format!("unknown variable \"{name}\"")),
    };
    match string(f, "type", at)? {
        "Variable" => {
            let j = var(string(f, "name", at)?, format!("{at}.name"))?;
            Ok(Affine {
                terms: vec![(j, 1.0)],
                constant: 0.0,
                single_variable: Some(j),
            })
        }
        "ScalarAffineFunction" => {
            let mut terms = Vec::new();
            for (k, t) in array(f, "terms", at)?.iter().enumerate() {
                let tat = format!("{at}.terms[{k}]");
                let j = var(string(t, "variable", &tat)?, format!("{tat}.variable"))?;
                terms.push((j, number(t, "coefficient", &tat)?));
            }
            Ok(Affine {
                terms,
                constant: number(f, "constant", at)?,
                single_variable: None,
            })
        }
        other => err(
            format!("{at}.type"),
            unsupported(&format!("function type \"{other}\"")),
        ),
    }
}

/// An interval set: (lower, upper), infinite where open.
fn interval(set: &Value, kind: &str, at: &str) -> Result<(f64, f64), MofError> {
    Ok(match kind {
        "LessThan" => (f64::NEG_INFINITY, number(set, "upper", at)?),
        "GreaterThan" => (number(set, "lower", at)?, f64::INFINITY),
        "EqualTo" => {
            let v = number(set, "value", at)?;
            (v, v)
        }
        "Interval" => (number(set, "lower", at)?, number(set, "upper", at)?),
        _ => unreachable!("caller matched the interval kinds"),
    })
}

/// Read a MathOptFormat JSON document into a [`Model`].
pub fn read_mof(bytes: &[u8]) -> Result<Model, MofError> {
    let doc: Value = match serde_json::from_slice(bytes) {
        Ok(v) => v,
        Err(e) => {
            return err(
                format!("line {} column {}", e.line(), e.column()),
                format!("invalid JSON: {e}"),
            );
        }
    };
    if !doc.is_object() {
        return err("top level", "expected a JSON object");
    }
    let version = field(&doc, "version", "top level")?;
    let major = number(version, "major", "version")?;
    if major != 1.0 {
        return err(
            "version.major",
            unsupported(&format!("MathOptFormat major version {major}")) + "; only version 1.x",
        );
    }
    let name = match doc.get("name") {
        Some(Value::String(s)) => s.clone(),
        _ => String::new(),
    };

    // Variables.
    let vars = array(&doc, "variables", "top level")?;
    let mut col_names = Vec::with_capacity(vars.len());
    let mut index: HashMap<&str, usize> = HashMap::with_capacity(vars.len());
    for (j, v) in vars.iter().enumerate() {
        let at = format!("variables[{j}]");
        let n = string(v, "name", &at)?;
        if index.insert(n, j).is_some() {
            return err(at, format!("duplicate variable name \"{n}\""));
        }
        col_names.push(n.to_string());
    }
    let n = vars.len();
    let mut bounds = vec![
        VariableBound {
            lower: f64::NEG_INFINITY,
            upper: f64::INFINITY,
            integer: false,
        };
        n
    ];

    // Objective.
    let obj = field(&doc, "objective", "top level")?;
    let mut objective = vec![0.0; n];
    let mut obj_offset = 0.0;
    let obj_sense = match string(obj, "sense", "objective")? {
        "min" => ObjSense::Minimize,
        "max" => ObjSense::Maximize,
        "feasibility" => ObjSense::Minimize,
        other => return err("objective.sense", format!("unknown sense \"{other}\"")),
    };
    if obj["sense"] != "feasibility" {
        let f = field(obj, "function", "objective")?;
        if let Some(
            t @ ("VectorAffineFunction"
            | "VectorOfVariables"
            | "VectorQuadraticFunction"
            | "VectorNonlinearFunction"),
        ) = f.get("type").and_then(Value::as_str)
        {
            return err(
                "objective.function.type",
                unsupported(&format!("multiple objectives (\"{t}\")")),
            );
        }
        let a = function(f, "objective.function", &index)?;
        for (j, c) in a.terms {
            objective[j] += c;
        }
        obj_offset = a.constant;
    }

    // Constraints.
    let mut senses = Vec::new();
    let mut rhs = Vec::new();
    let mut ranges = Vec::new();
    let mut row_names = Vec::new();
    let mut triplets: Vec<Triplet> = Vec::new();
    for (k, c) in array(&doc, "constraints", "top level")?.iter().enumerate() {
        let at = format!("constraints[{k}]");
        let set = field(c, "set", &at)?;
        let sat = format!("{at}.set");
        let kind = string(set, "type", &sat)?;
        let fat = format!("{at}.function");
        let f = function(field(c, "function", &at)?, &fat, &index)?;
        match (kind, f.single_variable) {
            ("Integer" | "ZeroOne", Some(j)) => {
                let b = &mut bounds[j];
                b.integer = true;
                if kind == "ZeroOne" {
                    b.lower = b.lower.max(0.0);
                    b.upper = b.upper.min(1.0);
                }
            }
            ("Integer" | "ZeroOne", None) => {
                return err(
                    fat,
                    unsupported(&format!(
                        "ScalarAffineFunction in {kind} (integrality of an expression)"
                    )),
                );
            }
            ("LessThan" | "GreaterThan" | "EqualTo" | "Interval", Some(j)) => {
                let (lo, hi) = interval(set, kind, &sat)?;
                let b = &mut bounds[j];
                b.lower = b.lower.max(lo);
                b.upper = b.upper.min(hi);
            }
            ("LessThan" | "GreaterThan" | "EqualTo" | "Interval", None) => {
                let (lo, hi) = interval(set, kind, &sat)?;
                let (lo, hi) = (lo - f.constant, hi - f.constant);
                let (sense, b, r) = if lo == hi {
                    (RowSense::Equal, lo, 0.0)
                } else if lo > hi {
                    return err(sat, format!("empty interval [{lo}, {hi}]"));
                } else if lo.is_finite() && hi.is_finite() {
                    (RowSense::LessEqual, hi, hi - lo)
                } else if hi.is_finite() {
                    (RowSense::LessEqual, hi, 0.0)
                } else if lo.is_finite() {
                    (RowSense::GreaterEqual, lo, 0.0)
                } else {
                    return err(sat, "interval with no finite bound");
                };
                let i = senses.len();
                senses.push(sense);
                rhs.push(b);
                ranges.push(r);
                row_names.push(match c.get("name") {
                    Some(Value::String(s)) => s.clone(),
                    _ => String::new(),
                });
                triplets.extend(f.terms.iter().map(|&(j, v)| (i, j, v)));
            }
            (other, _) => {
                return err(
                    format!("{sat}.type"),
                    unsupported(&format!("set type \"{other}\"")),
                );
            }
        }
    }
    let m = senses.len();
    let matrix = CscMatrix::from_triplets(m, n, &triplets);
    Ok(Model {
        name,
        num_variables: n,
        num_constraints: m,
        obj_sense,
        obj_name: None,
        obj_offset,
        objective,
        senses,
        rhs,
        ranges,
        bounds,
        matrix,
        row_names,
        col_names,
        hessian: Vec::new(),
        gc: Default::default(),
        conic: Default::default(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(s: &str) -> Result<Model, MofError> {
        read_mof(s.as_bytes())
    }

    /// `{version 1.7, variables x y, objective min x, constraints ...}`
    fn doc(objective: &str, constraints: &str) -> String {
        format!(
            r#"{{"version": {{"major": 1, "minor": 7}},
                "variables": [{{"name": "x"}}, {{"name": "y"}}],
                "objective": {objective},
                "constraints": [{constraints}]}}"#
        )
    }

    const MIN_X: &str = r#"{"sense": "min", "function": {"type": "Variable", "name": "x"}}"#;

    #[test]
    fn rows_bounds_and_constants() {
        let m = read(&doc(
            r#"{"sense": "max", "function": {"type": "ScalarAffineFunction",
                "terms": [{"coefficient": 2, "variable": "x"}, {"coefficient": 1, "variable": "x"}],
                "constant": -4.5}}"#,
            r#"{"name": "r", "function": {"type": "ScalarAffineFunction",
                 "terms": [{"coefficient": 1, "variable": "x"}, {"coefficient": -1, "variable": "y"}],
                 "constant": 1}, "set": {"type": "Interval", "lower": 2, "upper": 8}},
               {"function": {"type": "Variable", "name": "x"}, "set": {"type": "GreaterThan", "lower": -1}},
               {"function": {"type": "Variable", "name": "x"}, "set": {"type": "LessThan", "upper": 6}},
               {"function": {"type": "Variable", "name": "y"}, "set": {"type": "ZeroOne"}}"#,
        ))
        .unwrap();
        assert_eq!(m.obj_sense, ObjSense::Maximize);
        assert_eq!(m.objective, [3.0, 0.0], "duplicate terms summed");
        assert_eq!(m.obj_offset, -4.5);
        assert_eq!(m.num_constraints, 1);
        // 2 <= x - y + 1 <= 8  ->  1 <= x - y <= 7
        assert_eq!(
            (m.senses[0], m.rhs[0], m.ranges[0]),
            (RowSense::LessEqual, 7.0, 6.0)
        );
        assert_eq!(m.row_names, ["r"]);
        let b = |j: usize| (m.bounds[j].lower, m.bounds[j].upper, m.bounds[j].integer);
        assert_eq!(b(0), (-1.0, 6.0, false));
        assert_eq!(b(1), (0.0, 1.0, true));
    }

    #[test]
    fn variables_are_free_by_default_and_feasibility_has_no_objective() {
        let m = read(&doc(r#"{"sense": "feasibility"}"#, "")).unwrap();
        assert_eq!(m.bounds[0].lower, f64::NEG_INFINITY);
        assert_eq!(m.bounds[0].upper, f64::INFINITY);
        assert_eq!(m.objective, [0.0, 0.0]);
    }

    #[test]
    fn unsupported_content_is_named() {
        let msg = |o: &str, c: &str| read(&doc(o, c)).unwrap_err().to_string();
        let quad = r#"{"sense": "min", "function": {"type": "ScalarQuadraticFunction",
            "affine_terms": [], "quadratic_terms": [], "constant": 0}}"#;
        assert!(msg(quad, "").contains("function type \"ScalarQuadraticFunction\""));
        let vec = r#"{"sense": "min", "function": {"type": "VectorAffineFunction",
            "terms": [], "constants": [0, 0]}}"#;
        assert!(msg(vec, "").contains("multiple objectives"));
        let semi = r#"{"function": {"type": "Variable", "name": "x"},
            "set": {"type": "Semicontinuous", "lower": 1, "upper": 2}}"#;
        assert!(msg(MIN_X, semi).contains("set type \"Semicontinuous\""));
        let sos = r#"{"function": {"type": "VectorOfVariables", "variables": ["x", "y"]},
            "set": {"type": "SOS1", "weights": [1, 2]}}"#;
        assert!(msg(MIN_X, sos).contains("function type \"VectorOfVariables\""));
        let int_expr = r#"{"function": {"type": "ScalarAffineFunction",
            "terms": [{"coefficient": 1, "variable": "x"}], "constant": 0}, "set": {"type": "Integer"}}"#;
        assert!(msg(MIN_X, int_expr).contains("integrality of an expression"));
        let unknown =
            r#"{"function": {"type": "Variable", "name": "z"}, "set": {"type": "Integer"}}"#;
        assert!(
            msg(MIN_X, unknown).contains("constraints[0].function.name: unknown variable \"z\"")
        );
        let inf = r#"{"function": {"type": "Variable", "name": "x"}, "set": {"type": "LessThan", "upper": "Inf"}}"#;
        assert!(msg(MIN_X, inf).contains("expected a JSON number"));
        let v0 = doc(MIN_X, "").replace("\"major\": 1", "\"major\": 0");
        assert!(
            read(&v0)
                .unwrap_err()
                .to_string()
                .contains("major version 0")
        );
        let dup = doc(MIN_X, "").replace("\"name\": \"y\"", "\"name\": \"x\"");
        assert!(
            read(&dup)
                .unwrap_err()
                .to_string()
                .contains("duplicate variable name")
        );
        assert!(
            read("{\"version\": ")
                .unwrap_err()
                .to_string()
                .contains("invalid JSON")
        );
    }
}
