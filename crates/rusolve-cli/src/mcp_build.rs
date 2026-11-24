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

/// An entry of a `cp` list with its name.
type NamedEntry<'a> = (String, &'a Map<String, Value>);

/// The `cp` section of a `build_model` argument object: scheduling and combinatorial objects
/// over the declared variables (named by string). No file format carries them, so the handle
/// keeps the parsed [`rusolve_io::CpModel`] next to the model it writes.
///
/// ```json
/// "cp": {
///   "intervals":     [{"name": "a", "start": "s_a", "size": 3, "end": "e_a", "presence": "p_a"}],
///   "no_overlap":    [{"name": "m1", "intervals": ["a", "b"]}],
///   "cumulative":    [{"name": "crew", "intervals": ["a", "b"], "demands": [2, 1], "capacity": 3}],
///   "all_different": [{"name": "d", "variables": ["x", "y"]}],
///   "element":       [{"name": "e", "index": "i", "values": [4, "x", 7], "target": "t"}],
///   "table":         [{"name": "t", "variables": ["x", "y"], "tuples": [[0, 1], [1, 0]], "negated": false}]
/// }
/// ```
/// `size`, `demands`, `values` and `capacity` take a number or a variable name; `presence` is
/// optional. Every error names its JSON path (`cp.no_overlap[0].intervals[2]`).
pub fn parse_cp(args: &Args) -> Result<Option<rusolve_io::CpModel>, String> {
    use rusolve_io::{
        AllDifferent, CpModel, CpTerm, Cumulative, Element, Interval, NoOverlap, Table,
    };
    let Some(cp) = args.get("cp") else {
        return Ok(None);
    };
    let Some(cp) = cp.as_object() else {
        return Err(path_err(
            "cp",
            "must be an object with intervals, no_overlap, cumulative, all_different, element, table",
        ));
    };
    const SECTIONS: [&str; 6] = [
        "intervals",
        "no_overlap",
        "cumulative",
        "all_different",
        "element",
        "table",
    ];
    if let Some(k) = cp.keys().find(|k| !SECTIONS.contains(&k.as_str())) {
        return Err(path_err(
            &format!("cp.{k}"),
            format!("unknown section; the sections are {}", SECTIONS.join(", ")),
        ));
    }
    let mut index: HashMap<String, usize> = HashMap::new();
    for (i, v) in args
        .get("variables")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .enumerate()
    {
        if let Some(n) = v.get("name").and_then(Value::as_str) {
            index.insert(n.to_string(), i);
        }
    }
    let list = |key: &str| -> Result<Vec<NamedEntry<'_>>, String> {
        let Some(v) = cp.get(key) else {
            return Ok(Vec::new());
        };
        let Some(a) = v.as_array() else {
            return Err(path_err(&format!("cp.{key}"), "must be an array"));
        };
        a.iter()
            .enumerate()
            .map(|(i, o)| {
                let at = format!("cp.{key}[{i}]");
                let o = o
                    .as_object()
                    .ok_or_else(|| path_err(&at, "must be an object"))?;
                Ok((at, o))
            })
            .collect()
    };
    let var = |o: &Map<String, Value>, key: &str, at: &str| -> Result<usize, String> {
        let Some(n) = o.get(key).and_then(Value::as_str) else {
            return Err(path_err(
                &format!("{at}.{key}"),
                "is required (a variable name)",
            ));
        };
        index
            .get(n)
            .copied()
            .ok_or_else(|| path_err(&format!("{at}.{key}"), format!("no variable named '{n}'")))
    };
    let var_at = |v: &Value, at: &str| -> Result<usize, String> {
        let n = v
            .as_str()
            .ok_or_else(|| path_err(at, "must be a variable name"))?;
        index
            .get(n)
            .copied()
            .ok_or_else(|| path_err(at, format!("no variable named '{n}'")))
    };
    let term = |v: &Value, at: &str| -> Result<CpTerm, String> {
        match v {
            Value::String(n) => index
                .get(n.as_str())
                .map(|&j| CpTerm::Col(j))
                .ok_or_else(|| path_err(at, format!("no variable named '{n}'"))),
            Value::Number(x) => x
                .as_i64()
                .map(CpTerm::Const)
                .ok_or_else(|| path_err(at, "must be an integer or a variable name")),
            _ => Err(path_err(at, "must be an integer or a variable name")),
        }
    };
    let term_of = |o: &Map<String, Value>, key: &str, at: &str| -> Result<CpTerm, String> {
        match o.get(key) {
            Some(v) => term(v, &format!("{at}.{key}")),
            None => Err(path_err(
                &format!("{at}.{key}"),
                "is required (an integer or a variable name)",
            )),
        }
    };
    let array_of = |o: &Map<String, Value>, key: &str, at: &str| -> Result<Vec<Value>, String> {
        o.get(key)
            .and_then(Value::as_array)
            .cloned()
            .ok_or_else(|| path_err(&format!("{at}.{key}"), "is required (an array)"))
    };
    let mut out = CpModel::default();
    let mut interval_ix: HashMap<String, usize> = HashMap::new();
    for (at, o) in list("intervals")? {
        let name = name_of(o, &at)?;
        if interval_ix
            .insert(name.clone(), out.intervals.len())
            .is_some()
        {
            return Err(path_err(&at, format!("duplicate interval name '{name}'")));
        }
        let presence = match o.get("presence") {
            None | Some(Value::Null) => None,
            Some(_) => Some(var(o, "presence", &at)?),
        };
        out.intervals.push(Interval {
            name,
            start: var(o, "start", &at)?,
            size: term_of(o, "size", &at)?,
            end: var(o, "end", &at)?,
            presence,
        });
    }
    let ivs = |o: &Map<String, Value>, at: &str| -> Result<Vec<usize>, String> {
        array_of(o, "intervals", at)?
            .iter()
            .enumerate()
            .map(|(k, v)| {
                let here = format!("{at}.intervals[{k}]");
                let n = v
                    .as_str()
                    .ok_or_else(|| path_err(&here, "must be an interval name"))?;
                interval_ix
                    .get(n)
                    .copied()
                    .ok_or_else(|| path_err(&here, format!("no interval named '{n}'")))
            })
            .collect()
    };
    for (at, o) in list("no_overlap")? {
        out.no_overlap.push(NoOverlap {
            name: name_of(o, &at)?,
            intervals: ivs(o, &at)?,
        });
    }
    for (at, o) in list("cumulative")? {
        let intervals = ivs(o, &at)?;
        let demands = array_of(o, "demands", &at)?
            .iter()
            .enumerate()
            .map(|(k, v)| term(v, &format!("{at}.demands[{k}]")))
            .collect::<Result<Vec<_>, _>>()?;
        if demands.len() != intervals.len() {
            return Err(path_err(
                &format!("{at}.demands"),
                format!(
                    "has {} entries but intervals has {}; one demand per interval",
                    demands.len(),
                    intervals.len()
                ),
            ));
        }
        out.cumulative.push(Cumulative {
            name: name_of(o, &at)?,
            intervals,
            demands,
            capacity: term_of(o, "capacity", &at)?,
        });
    }
    for (at, o) in list("all_different")? {
        let cols = array_of(o, "variables", &at)?
            .iter()
            .enumerate()
            .map(|(k, v)| var_at(v, &format!("{at}.variables[{k}]")))
            .collect::<Result<Vec<_>, _>>()?;
        out.all_different.push(AllDifferent {
            name: name_of(o, &at)?,
            cols,
        });
    }
    for (at, o) in list("element")? {
        let values = array_of(o, "values", &at)?
            .iter()
            .enumerate()
            .map(|(k, v)| term(v, &format!("{at}.values[{k}]")))
            .collect::<Result<Vec<_>, _>>()?;
        out.element.push(Element {
            name: name_of(o, &at)?,
            index: var(o, "index", &at)?,
            values,
            target: var(o, "target", &at)?,
        });
    }
    for (at, o) in list("table")? {
        let cols = array_of(o, "variables", &at)?
            .iter()
            .enumerate()
            .map(|(k, v)| var_at(v, &format!("{at}.variables[{k}]")))
            .collect::<Result<Vec<_>, _>>()?;
        let mut tuples = Vec::new();
        for (r, row) in array_of(o, "tuples", &at)?.iter().enumerate() {
            let here = format!("{at}.tuples[{r}]");
            let row = row
                .as_array()
                .ok_or_else(|| path_err(&here, "must be an array of integers"))?;
            if row.len() != cols.len() {
                return Err(path_err(
                    &here,
                    format!(
                        "has {} values but the table has {} variables",
                        row.len(),
                        cols.len()
                    ),
                ));
            }
            tuples.push(
                row.iter()
                    .enumerate()
                    .map(|(k, v)| {
                        v.as_i64()
                            .ok_or_else(|| path_err(&format!("{here}[{k}]"), "must be an integer"))
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            );
        }
        let negated = match o.get("negated") {
            None => false,
            Some(v) => v
                .as_bool()
                .ok_or_else(|| path_err(&format!("{at}.negated"), "must be true or false"))?,
        };
        out.table.push(Table {
            name: name_of(o, &at)?,
            cols,
            tuples,
            negated,
        });
    }
    Ok(Some(out))
}

/// The JSON path of the `cp` entry a [`rusolve_io::CpError`] names (`interval 'a'` becomes
/// `cp.intervals[0]`), so a validation error points at its entry like a parse error does.
pub fn cp_error_path(cp: &rusolve_io::CpModel, object: &str) -> Option<String> {
    let (kind, rest) = object.split_once(" '")?;
    let name = rest.strip_suffix('\'')?;
    let (section, names): (&str, Vec<&str>) = match kind {
        "interval" => (
            "intervals",
            cp.intervals.iter().map(|o| o.name.as_str()).collect(),
        ),
        "no_overlap" => (
            "no_overlap",
            cp.no_overlap.iter().map(|o| o.name.as_str()).collect(),
        ),
        "cumulative" => (
            "cumulative",
            cp.cumulative.iter().map(|o| o.name.as_str()).collect(),
        ),
        "all_different" => (
            "all_different",
            cp.all_different.iter().map(|o| o.name.as_str()).collect(),
        ),
        "element" => (
            "element",
            cp.element.iter().map(|o| o.name.as_str()).collect(),
        ),
        "table" => ("table", cp.table.iter().map(|o| o.name.as_str()).collect()),
        _ => return None,
    };
    names
        .iter()
        .position(|n| *n == name)
        .map(|i| format!("cp.{section}[{i}]"))
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

    #[test]
    fn cp_section_parses_and_errors_name_the_path() {
        let a = args(json!({
            "variables": [
                {"name": "s1", "type": "integer", "upper": 9},
                {"name": "e1", "type": "integer", "upper": 9},
                {"name": "s2", "type": "integer", "upper": 9},
                {"name": "e2", "type": "integer", "upper": 9}
            ],
            "cp": {
                "intervals": [
                    {"name": "a", "start": "s1", "size": 2, "end": "e1"},
                    {"name": "b", "start": "s2", "size": "s1", "end": "e2"}
                ],
                "no_overlap": [{"name": "m", "intervals": ["a", "b"]}],
                "cumulative": [{"name": "c", "intervals": ["a", "b"], "demands": [1, 2], "capacity": 2}]
            }
        }));
        let cp = parse_cp(&a).unwrap().unwrap();
        assert_eq!(
            (cp.intervals.len(), cp.no_overlap.len(), cp.cumulative.len()),
            (2, 1, 1)
        );
        assert_eq!(parse_cp(&args(json!({"variables": []}))).unwrap(), None);
        let mut bad = a.clone();
        bad["cp"]["no_overlap"][0]["intervals"][1] = json!("zz");
        let e = parse_cp(&bad).unwrap_err();
        assert!(
            e.starts_with("cp.no_overlap[0].intervals[1]:") && e.contains("zz"),
            "{e}"
        );
        let mut bad = a.clone();
        bad["cp"]["cumulative"][0]["demands"] = json!([1]);
        assert!(
            parse_cp(&bad)
                .unwrap_err()
                .starts_with("cp.cumulative[0].demands:")
        );
        let mut bad = a.clone();
        bad["cp"]["intervals"][0]["end"] = json!("nope");
        assert!(
            parse_cp(&bad)
                .unwrap_err()
                .starts_with("cp.intervals[0].end:")
        );
        let mut bad = a;
        bad["cp"]["bogus"] = json!([]);
        assert!(parse_cp(&bad).unwrap_err().starts_with("cp.bogus:"));
        assert_eq!(
            cp_error_path(&cp, "interval 'b'").as_deref(),
            Some("cp.intervals[1]")
        );
    }
}
