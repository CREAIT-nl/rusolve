//! `rusolve diagnose <model.mps>` (formulation or data?), through the real
//! binary, on generated family-structured models with planted errors.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

/// `docs/schema/rusolve-diagnose-1.schema.json`, the normative output schema.
const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-diagnose-1.schema.json");

/// Fixed-format-free MPS writer for the generators below.
struct Mps {
    rows: Vec<(char, String)>,
    cols: Vec<(String, Vec<(String, f64)>)>,
    rhs: Vec<(String, f64)>,
}

impl Mps {
    fn new() -> Self {
        Mps {
            rows: vec![('N', "obj".into())],
            cols: Vec::new(),
            rhs: Vec::new(),
        }
    }
    fn row(&mut self, sense: char, name: &str, rhs: f64) {
        self.rows.push((sense, name.into()));
        self.rhs.push((name.into(), rhs));
    }
    fn col(&mut self, name: &str, cost: f64, coefs: Vec<(String, f64)>) {
        let mut v = vec![("obj".to_string(), cost)];
        v.extend(coefs);
        self.cols.push((name.into(), v));
    }
    fn text(&self) -> String {
        let mut s = String::from("NAME GEN\nROWS\n");
        for (k, n) in &self.rows {
            s.push_str(&format!(" {k}  {n}\n"));
        }
        s.push_str("COLUMNS\n");
        for (c, v) in &self.cols {
            for (r, x) in v {
                s.push_str(&format!("    {c}  {r}  {x}\n"));
            }
        }
        s.push_str("RHS\n");
        for (r, x) in &self.rhs {
            s.push_str(&format!("    rhs  {r}  {x}\n"));
        }
        s.push_str("ENDATA\n");
        s
    }
}

const NS: usize = 6;
const ND: usize = 8;

/// Transportation model with families `supply[i]` (<= s_i), `demand[j]`
/// (>= d_j or <= d_j when flipped), `minship[i]` (>= m_i) and arc columns
/// `x[i,j]`. Feasible as generated: total supply 660, demand 540, minimum
/// shipment 0.7 of supply.
fn transport(demand_scale_at: Option<(usize, f64)>, flip_demand: bool, min_frac: f64) -> String {
    let mut m = Mps::new();
    let s = |i: usize| 100.0 + 10.0 * (i % 3) as f64;
    let d = |j: usize| {
        let base = 60.0 + 5.0 * (j % 4) as f64;
        match demand_scale_at {
            Some((k, f)) if k == j => base * f,
            _ => base,
        }
    };
    for i in 0..NS {
        m.row('L', &format!("supply[{i}]"), s(i));
    }
    for j in 0..ND {
        m.row(
            if flip_demand { 'L' } else { 'G' },
            &format!("demand[{j}]"),
            d(j),
        );
    }
    for i in 0..NS {
        m.row('G', &format!("minship[{i}]"), (min_frac * s(i)).round());
    }
    for i in 0..NS {
        for j in 0..ND {
            m.col(
                &format!("x[{i},{j}]"),
                1.0 + ((i * 7 + j * 3) % 5) as f64,
                vec![
                    (format!("supply[{i}]"), 1.0),
                    (format!("demand[{j}]"), 1.0),
                    (format!("minship[{i}]"), 1.0),
                ],
            );
        }
    }
    m.text()
}

/// Assignment model `assign_row_i`, `assign_col_j` (names with `_` indices).
/// `row_coef` 2 double-counts the row family; `rhs_at` plants one row rhs.
fn assignment(row_coef: f64, rhs_at: Option<(usize, f64)>) -> String {
    let n = 6;
    let mut m = Mps::new();
    for i in 0..n {
        let rhs = match rhs_at {
            Some((k, v)) if k == i => v,
            _ => 1.0,
        };
        m.row('E', &format!("assign_row_{i}"), rhs);
    }
    for j in 0..n {
        m.row('E', &format!("assign_col_{j}"), 1.0);
    }
    for i in 0..n {
        for j in 0..n {
            m.col(
                &format!("y_{i}_{j}"),
                ((i * 3 + j * 5) % 7) as f64,
                vec![
                    (format!("assign_row_{i}"), row_coef),
                    (format!("assign_col_{j}"), 1.0),
                ],
            );
        }
    }
    m.text()
}

fn write_tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-diagnose-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn diagnose_raw(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("diagnose")
        .args(args)
        .output()
        .unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

fn diagnose(name: &str, mps: &str) -> Value {
    let p = write_tmp(name, mps);
    let (code, stdout, stderr) = diagnose_raw(&[p.to_str().unwrap(), "--format", "json"]);
    assert_eq!(code, 0, "exit code, stderr: {stderr}");
    serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not one JSON object: {e}\n{stdout}"))
}

fn signals(d: &Value) -> Vec<String> {
    d["evidence"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| e["signal"].as_str().unwrap().to_string())
        .collect()
}

#[test]
fn diagnose_feasible_model_reports_nothing() {
    for (n, m) in [
        ("tfeas.mps", transport(None, false, 0.7)),
        ("afeas.mps", assignment(1.0, None)),
    ] {
        let d = diagnose(n, &m);
        assert_eq!(d["verdict"], "none", "{d}");
        assert_eq!(d["status"], "Optimal", "{d}");
        assert!(d["suggested_fix"].is_null());
        assert_eq!(d["evidence"].as_array().unwrap().len(), 0);
    }
}

#[test]
fn diagnose_transport_demand_100x_is_data() {
    let d = diagnose("tdata.mps", &transport(Some((3, 100.0)), false, 0.7));
    assert_eq!(d["status"], "Infeasible");
    assert_eq!(d["verdict"], "data", "{d}");
    assert!(d["confidence"].as_f64().unwrap() >= 0.5, "{d}");
    let fix = &d["suggested_fix"];
    assert_eq!(fix["type"], "data");
    let cells = fix["cells"].as_array().unwrap();
    assert_eq!(cells.len(), 1, "{fix}");
    assert_eq!(cells[0]["name"], "demand[3]");
    assert!(signals(&d).contains(&"b_median_fix_restores_feasibility".to_string()));
}

#[test]
fn diagnose_transport_flipped_demand_direction_is_formulation() {
    // demand[j] <= d_j instead of >= d_j, with minimum shipments above total demand
    let d = diagnose("tflip.mps", &transport(None, true, 0.95));
    assert_eq!(d["status"], "Infeasible");
    assert_eq!(d["verdict"], "formulation", "{d}");
    let fix = &d["suggested_fix"];
    assert_eq!(fix["type"], "formulation");
    let fams: Vec<&str> = fix["families"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f.as_str().unwrap())
        .collect();
    assert!(fams.contains(&"demand[*]"), "{fams:?}");
    assert!(fams.contains(&"minship[*]"), "{fams:?}");
    assert!(signals(&d).contains(&"a_typical_data_still_infeasible".to_string()));
}

#[test]
fn diagnose_assignment_double_counted_balance_is_formulation() {
    let d = diagnose("adouble.mps", &assignment(2.0, None));
    assert_eq!(d["status"], "Infeasible");
    assert_eq!(d["verdict"], "formulation", "{d}");
    assert_eq!(d["suggested_fix"]["type"], "formulation");
    let fams = d["families"].as_array().unwrap();
    let stems: Vec<&str> = fams.iter().map(|f| f["stem"].as_str().unwrap()).collect();
    assert!(
        stems.contains(&"assign_row") && stems.contains(&"assign_col"),
        "{stems:?}"
    );
}

#[test]
fn diagnose_assignment_rhs_100x_is_data() {
    let d = diagnose("adata.mps", &assignment(1.0, Some((3, 100.0))));
    assert_eq!(d["status"], "Infeasible");
    assert_eq!(d["verdict"], "data", "{d}");
    let cells = d["suggested_fix"]["cells"].as_array().unwrap();
    assert_eq!(cells[0]["name"], "assign_row_3");
}

#[test]
fn diagnose_unnamed_structure_is_not_overconfident() {
    // x + y <= 3, x >= 2, y >= 2 under names without families
    let mps = "NAME T\nROWS\n N obj\n L cap\n G nx\n G ny\nCOLUMNS\n x obj 1 cap 1\n x nx 1\n y obj 1 cap 1\n y ny 1\nRHS\n rhs cap 3 nx 2\n rhs ny 2\nENDATA\n";
    let d = diagnose("tiny.mps", mps);
    assert_eq!(d["status"], "Infeasible");
    assert_eq!(d["naming_informative"], false);
    assert!(d["confidence"].as_f64().unwrap() <= 0.6, "{d}");
}

#[test]
fn diagnose_usage_and_input_errors() {
    let (code, _, _) = diagnose_raw(&[]);
    assert_eq!(code, 2);
    let (code, _, _) = diagnose_raw(&["/no/such/model.mps"]);
    assert_eq!(code, 3);
}

fn conforms(v: &Value, s: &Value, root: &Value, path: &str, errs: &mut Vec<String>) -> bool {
    let mut local = Vec::new();
    if let Some(r) = s.get("$ref").and_then(Value::as_str) {
        let name = r.trim_start_matches("#/$defs/");
        return conforms(v, &root["$defs"][name], root, path, errs);
    }
    if let Some(alts) = s.get("oneOf").and_then(Value::as_array) {
        let ok = alts
            .iter()
            .filter(|a| conforms(v, a, root, path, &mut Vec::new()))
            .count();
        if ok != 1 {
            errs.push(format!("{path}: matches {ok} of the oneOf alternatives"));
            return false;
        }
    }
    if let Some(t) = s.get("type") {
        let types: Vec<&str> = match t {
            Value::String(x) => vec![x.as_str()],
            Value::Array(a) => a.iter().filter_map(Value::as_str).collect(),
            _ => vec![],
        };
        let hit = types.iter().any(|t| match *t {
            "null" => v.is_null(),
            "object" => v.is_object(),
            "array" => v.is_array(),
            "string" => v.is_string(),
            "boolean" => v.is_boolean(),
            "integer" => v.is_i64() || v.is_u64(),
            "number" => v.is_number(),
            _ => false,
        });
        if !hit {
            local.push(format!("{path}: {v} is not of type {t}"));
        }
    }
    if let Some(c) = s.get("const").filter(|c| *c != v) {
        local.push(format!("{path}: {v} != const {c}"));
    }
    if let Some(e) = s
        .get("enum")
        .and_then(Value::as_array)
        .filter(|e| !e.contains(v))
    {
        local.push(format!("{path}: {v} not in enum {e:?}"));
    }
    if let Some(obj) = v.as_object() {
        for f in s
            .get("required")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
        {
            let f = f.as_str().unwrap();
            if !obj.contains_key(f) {
                local.push(format!("{path}.{f} is required and missing"));
            }
        }
        if let Some(props) = s.get("properties").and_then(Value::as_object) {
            for (k, sub) in props {
                if let Some(val) = obj.get(k) {
                    conforms(val, sub, root, &format!("{path}.{k}"), &mut local);
                }
            }
            if s.get("additionalProperties") == Some(&Value::Bool(false)) {
                for k in obj.keys() {
                    if !props.contains_key(k) {
                        local.push(format!("{path}.{k} is not in the schema"));
                    }
                }
            }
        }
    }
    if let (Some(arr), Some(items)) = (v.as_array(), s.get("items")) {
        for (i, x) in arr.iter().enumerate() {
            conforms(x, items, root, &format!("{path}[{i}]"), &mut local);
        }
    }
    let ok = local.is_empty();
    errs.extend(local);
    ok
}

#[test]
fn diagnose_json_matches_its_schema() {
    let schema: Value = serde_json::from_str(SCHEMA).unwrap();
    assert_eq!(
        schema["properties"]["schema_version"]["const"],
        "rusolve.diagnose/1"
    );
    for (name, mps) in [
        ("s1.mps", transport(Some((3, 100.0)), false, 0.7)),
        ("s2.mps", transport(None, true, 0.95)),
        ("s3.mps", assignment(1.0, None)),
        ("s4.mps", assignment(2.0, None)),
    ] {
        let d = diagnose(name, &mps);
        let mut errs = Vec::new();
        conforms(&d, &schema, &schema, "$", &mut errs);
        assert!(
            errs.is_empty(),
            "{name} does not match the schema: {errs:?}\n{d}"
        );
        for f in schema["required"].as_array().unwrap() {
            assert!(d.get(f.as_str().unwrap()).is_some(), "{name}: missing {f}");
        }
    }
}
