//! The `build_model` tool. An agent that has no MPS or LP file describes the model as JSON
//! (variables with bounds and type, rows with named terms, an objective) and gets a handle.
//!
//! The JSON is written out as free-format MPS and read back by the same reader every other model
//! goes through, so a built model and a loaded one are the same object and nothing in the solver
//! path is special. MPS layout: IBM MPSX / the free-format conventions the reader takes
//! (<https://lpsolve.sourceforge.net/5.5/mps-format.htm>): `OBJSENSE` section, `ROWS`, `COLUMNS`
//! with `MARKER` lines around integers, `RHS` (the objective row's entry is minus the constant),
//! `RANGES`, `BOUNDS`. Every column gets explicit bound lines, so no default (an integer column
//! without bounds is `[0, 1]` in some readers) is ever relied on.
//!
//! Names are validated here, with the JSON path of the offender, because a name with a space
//! cannot be written to free-format MPS.

use serde_json::{Map, Value};
use std::collections::HashMap;
use std::fmt::Write;

type Args = Map<String, Value>;

fn path_err(at: &str, msg: impl AsRef<str>) -> String {
    format!("{at}: {}", msg.as_ref())
}

fn name_of(v: &Map<String, Value>, at: &str) -> Result<String, String> {
    let Some(n) = v.get("name").and_then(Value::as_str) else {
        return Err(path_err(at, "\"name\" is required (a string)"));
    };
    if n.is_empty() || n.chars().any(char::is_whitespace) {
        return Err(path_err(
            at,
            format!("name {n:?} must be non-empty and contain no whitespace"),
        ));
    }
    Ok(n.to_string())
}

/// A bound: absent is `default`; `null` is the infinity of `inf_sign`.
fn bound(
    v: &Map<String, Value>,
    key: &str,
    default: f64,
    inf: f64,
    at: &str,
) -> Result<f64, String> {
    match v.get(key) {
        None => Ok(default),
        Some(Value::Null) => Ok(inf),
        Some(x) => x
            .as_f64()
            .ok_or_else(|| path_err(&format!("{at}.{key}"), "must be a number or null")),
    }
}

fn number(v: &Map<String, Value>, key: &str, at: &str) -> Result<Option<f64>, String> {
    match v.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(x) => x
            .as_f64()
            .map(Some)
            .ok_or_else(|| path_err(&format!("{at}.{key}"), "must be a number")),
    }
}

fn fmt(x: f64) -> String {
    // `{}` of an f64 is the shortest string that reads back to the same value
    format!("{x}")
}

struct Var {
    name: String,
    lower: f64,
    upper: f64,
    integer: bool,
    cost: f64,
}

/// `terms` as `{"var": coef}` or `[{"var": "x", "coef": 2}]`, in the order given.
fn terms_of(
    v: Option<&Value>,
    vars: &HashMap<String, usize>,
    at: &str,
) -> Result<Vec<(usize, f64)>, String> {
    let at = format!("{at}.terms");
    let mut out: Vec<(usize, f64)> = Vec::new();
    let mut add = |name: &str, c: &Value, here: &str| -> Result<(), String> {
        let Some(&j) = vars.get(name) else {
            return Err(path_err(
                here,
                format!("unknown variable '{name}' (declare it in \"variables\")"),
            ));
        };
        let c = c
            .as_f64()
            .ok_or_else(|| path_err(here, "a coefficient must be a number"))?;
        if out.iter().any(|(k, _)| *k == j) {
            return Err(path_err(here, format!("variable '{name}' appears twice")));
        }
        out.push((j, c));
        Ok(())
    };
    match v {
        None | Some(Value::Null) => {}
        Some(Value::Object(m)) => {
            for (k, c) in m {
                add(k, c, &format!("{at}.{k}"))?;
            }
        }
        Some(Value::Array(a)) => {
            for (i, t) in a.iter().enumerate() {
                let here = format!("{at}[{i}]");
                let Some(o) = t.as_object() else {
                    return Err(path_err(
                        &here,
                        "must be an object {\"var\": name, \"coef\": number}",
                    ));
                };
                let Some(n) = o.get("var").and_then(Value::as_str) else {
                    return Err(path_err(&here, "\"var\" is required (a variable name)"));
                };
                let c = o
                    .get("coef")
                    .ok_or_else(|| path_err(&here, "\"coef\" is required (a number)"))?;
                add(n, c, &here)?;
            }
        }
        Some(_) => {
            return Err(path_err(
                &at,
                "must be an object {variable: coefficient} or an array of {var, coef}",
            ));
        }
    }
    Ok(out)
}

/// The MPS text of a `build_model` argument object.
pub fn to_mps(args: &Args) -> Result<String, String> {
    let sense = match args.get("sense").and_then(Value::as_str) {
        None | Some("minimize") => false,
        Some("maximize") => true,
        Some(other) => {
            return Err(path_err(
                "sense",
                format!("{other:?} is not \"minimize\" or \"maximize\""),
            ));
        }
    };
    let model_name = args
        .get("name")
        .and_then(Value::as_str)
        .filter(|n| !n.is_empty() && !n.chars().any(char::is_whitespace))
        .unwrap_or("model");
    let Some(vars_json) = args.get("variables").and_then(Value::as_array) else {
        return Err(path_err("variables", "is required (an array of variables)"));
    };
    let mut vars: Vec<Var> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    for (i, v) in vars_json.iter().enumerate() {
        let at = format!("variables[{i}]");
        let Some(o) = v.as_object() else {
            return Err(path_err(&at, "must be an object"));
        };
        let name = name_of(o, &at)?;
        if index.insert(name.clone(), i).is_some() {
            return Err(path_err(&at, format!("duplicate variable name '{name}'")));
        }
        let kind = o
            .get("type")
            .and_then(Value::as_str)
            .unwrap_or("continuous");
        let (mut lower, mut upper, integer) = match kind {
            "continuous" => (0.0, f64::INFINITY, false),
            "integer" => (0.0, f64::INFINITY, true),
            "binary" => (0.0, 1.0, true),
            other => {
                return Err(path_err(
                    &format!("{at}.type"),
                    format!("{other:?} is not continuous, integer or binary"),
                ));
            }
        };
        lower = bound(o, "lower", lower, f64::NEG_INFINITY, &at)?;
        upper = bound(o, "upper", upper, f64::INFINITY, &at)?;
        if lower > upper {
            return Err(path_err(
                &at,
                format!("variable '{name}' has lower {lower} > upper {upper}"),
            ));
        }
        let cost = number(o, "objective", &at)?.unwrap_or(0.0);
        vars.push(Var {
            name,
            lower,
            upper,
            integer,
            cost,
        });
    }
    // the objective: per-variable `objective` and/or an `objective` object with terms and offset
    let mut offset = 0.0;
    if let Some(ob) = args.get("objective") {
        let Some(ob) = ob.as_object() else {
            return Err(path_err(
                "objective",
                "must be an object {\"terms\": {variable: coefficient}, \"offset\": number}",
            ));
        };
        for (j, c) in terms_of(ob.get("terms"), &index, "objective")? {
            vars[j].cost += c;
        }
        offset = number(ob, "offset", "objective")?.unwrap_or(0.0);
    }
    struct Row {
        name: String,
        terms: Vec<(usize, f64)>,
        lower: f64,
        upper: f64,
    }
    let mut rows: Vec<Row> = Vec::new();
    let mut row_index: HashMap<String, usize> = HashMap::new();
    if let Some(rj) = args.get("rows") {
        let Some(rj) = rj.as_array() else {
            return Err(path_err("rows", "must be an array of rows"));
        };
        for (i, r) in rj.iter().enumerate() {
            let at = format!("rows[{i}]");
            let Some(o) = r.as_object() else {
                return Err(path_err(&at, "must be an object"));
            };
            let name = name_of(o, &at)?;
            if row_index.insert(name.clone(), i).is_some() {
                return Err(path_err(&at, format!("duplicate row name '{name}'")));
            }
            let terms = terms_of(o.get("terms"), &index, &at)?;
            let (lower, upper) = match o.get("sense").and_then(Value::as_str) {
                Some(s) => {
                    let Some(rhs) = number(o, "rhs", &at)? else {
                        return Err(path_err(&at, "\"sense\" needs \"rhs\" (a number)"));
                    };
                    if o.contains_key("lower") || o.contains_key("upper") {
                        return Err(path_err(
                            &at,
                            "give \"sense\"+\"rhs\" or \"lower\"/\"upper\", not both",
                        ));
                    }
                    match s {
                        "<=" | "L" => (f64::NEG_INFINITY, rhs),
                        ">=" | "G" => (rhs, f64::INFINITY),
                        "=" | "==" | "E" => (rhs, rhs),
                        other => {
                            return Err(path_err(
                                &format!("{at}.sense"),
                                format!("{other:?} is not <=, >= or ="),
                            ));
                        }
                    }
                }
                None => (
                    bound(o, "lower", f64::NEG_INFINITY, f64::NEG_INFINITY, &at)?,
                    bound(o, "upper", f64::INFINITY, f64::INFINITY, &at)?,
                ),
            };
            if lower > upper {
                return Err(path_err(
                    &at,
                    format!("row '{name}' has lower {lower} > upper {upper}"),
                ));
            }
            rows.push(Row {
                name,
                terms,
                lower,
                upper,
            });
        }
    }
    if vars.is_empty() {
        return Err(path_err("variables", "must declare at least one variable"));
    }
    // ---- write it
    let obj = "obj_";
    let obj_name = if row_index.contains_key(obj) {
        "obj__rusolve"
    } else {
        obj
    };
    let mut s = String::new();
    writeln!(s, "NAME {model_name}").unwrap();
    if sense {
        s.push_str("OBJSENSE\n    MAX\n");
    }
    s.push_str("ROWS\n");
    writeln!(s, " N {obj_name}").unwrap();
    // kind of each row: (type letter, rhs, range)
    let mut kinds: Vec<(char, f64, Option<f64>)> = Vec::new();
    for r in &rows {
        let (lo, up) = (r.lower, r.upper);
        let k = match (lo.is_finite(), up.is_finite()) {
            (false, false) => ('N', 0.0, None), // a free row: no constraint
            (false, true) => ('L', up, None),
            (true, false) => ('G', lo, None),
            (true, true) if lo == up => ('E', lo, None),
            (true, true) => ('G', lo, Some(up - lo)),
        };
        writeln!(s, " {} {}", k.0, r.name).unwrap();
        kinds.push(k);
    }
    s.push_str("COLUMNS\n");
    let mut by_col: Vec<Vec<(usize, f64)>> = vec![Vec::new(); vars.len()];
    for (i, r) in rows.iter().enumerate() {
        for &(j, c) in &r.terms {
            by_col[j].push((i, c));
        }
    }
    let mut in_int = false;
    for (j, v) in vars.iter().enumerate() {
        if v.integer && !in_int {
            s.push_str("    MARKER MARKER INTORG\n");
            in_int = true;
        } else if !v.integer && in_int {
            s.push_str("    MARKER MARKER INTEND\n");
            in_int = false;
        }
        let mut wrote = false;
        if v.cost != 0.0 {
            writeln!(s, "    {} {obj_name} {}", v.name, fmt(v.cost)).unwrap();
            wrote = true;
        }
        for &(i, c) in &by_col[j] {
            writeln!(s, "    {} {} {}", v.name, rows[i].name, fmt(c)).unwrap();
            wrote = true;
        }
        if !wrote {
            // a column with no entry still has to exist
            writeln!(s, "    {} {obj_name} 0", v.name).unwrap();
        }
    }
    if in_int {
        s.push_str("    MARKER MARKER INTEND\n");
    }
    s.push_str("RHS\n");
    if offset != 0.0 {
        writeln!(s, "    rhs {obj_name} {}", fmt(-offset)).unwrap();
    }
    for (r, k) in rows.iter().zip(&kinds) {
        if k.1 != 0.0 {
            writeln!(s, "    rhs {} {}", r.name, fmt(k.1)).unwrap();
        }
    }
    if kinds.iter().any(|k| k.2.is_some()) {
        s.push_str("RANGES\n");
        for (r, k) in rows.iter().zip(&kinds) {
            if let Some(range) = k.2 {
                writeln!(s, "    rng {} {}", r.name, fmt(range)).unwrap();
            }
        }
    }
    s.push_str("BOUNDS\n");
    for v in &vars {
        let n = &v.name;
        match (v.lower.is_finite(), v.upper.is_finite()) {
            (false, false) => writeln!(s, " FR bnd {n}").unwrap(),
            (true, true) if v.lower == v.upper => {
                writeln!(s, " FX bnd {n} {}", fmt(v.lower)).unwrap()
            }
            (lo, up) => {
                if lo {
                    writeln!(s, " LO bnd {n} {}", fmt(v.lower)).unwrap();
                } else {
                    writeln!(s, " MI bnd {n}").unwrap();
                }
                if up {
                    writeln!(s, " UP bnd {n} {}", fmt(v.upper)).unwrap();
                }
            }
        }
    }
    s.push_str("ENDATA\n");
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn args(v: Value) -> Args {
        v.as_object().unwrap().clone()
    }

    #[test]
    fn writes_a_small_model() {
        let m = to_mps(&args(json!({
            "sense": "maximize",
            "variables": [
                {"name": "x", "upper": 4, "objective": 3},
                {"name": "y", "type": "binary", "objective": 2}
            ],
            "rows": [{"name": "c", "terms": {"x": 1, "y": 2}, "sense": "<=", "rhs": 6}]
        })))
        .unwrap();
        assert!(m.contains("OBJSENSE\n    MAX"));
        assert!(m.contains(" L c"));
        assert!(m.contains("MARKER MARKER INTORG"));
        assert!(m.contains(" UP bnd x 4"));
    }

    #[test]
    fn errors_name_the_json_path() {
        let e = to_mps(&args(json!({
            "variables": [{"name": "x"}],
            "rows": [{"name": "c", "terms": {"z": 1}, "upper": 1}]
        })))
        .unwrap_err();
        assert!(
            e.starts_with("rows[0].terms.z:") && e.contains("'z'"),
            "{e}"
        );
        let e = to_mps(&args(json!({"variables": [{"name": "x y"}]}))).unwrap_err();
        assert!(e.starts_with("variables[0]:"), "{e}");
    }
}
