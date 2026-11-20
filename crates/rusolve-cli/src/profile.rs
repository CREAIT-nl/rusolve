//! `rusolve profile` and the profile file format.
//!
//! A profile is a small versioned JSON file (`"schema_version": "rusolve.profile/1"`):
//! a base (`default` = the HiGHS 1.15.1 settings, `recommended` = what the CLI starts from),
//! the options that differ from that base, an optional feature fingerprint of the instances it
//! was tuned on, and the provenance and evidence of the run that produced it.
//!
//! ```text
//! rusolve profile diff <a> <b>          # a, b: a profile file, `default` or `recommended`
//! rusolve profile show <file>           # the options a file changes, and its provenance
//! rusolve profile validate <file>       # exit 0, or errors with nearest-name hints
//! rusolve profile new [--base B] [--name N] [--opt k=v]... [-o file]
//! rusolve mip model.mps --profile file  # (mip_highs.rs) apply one
//! ```
//!
//! ## Why JSON and not TOML
//!
//! `serde_json` is already a dependency of this crate, the `--format json` outputs, the MCP
//! tools and the Python crate speak JSON, and an agent writes and reads JSON without a second
//! parser. TOML would add a crate and a second syntax for the same data. The price is that
//! JSON has no comments, so a profile carries a free-text `name` and the `provenance` object
//! instead. Closed schema: an unknown top-level key is an error, because the schema version is
//! what says which keys exist (`docs/schema/rusolve-profile-1.schema.json`).
//!
//! ## What a profile may set
//!
//! Only registry options whose `affects` class is `search` (the strategy knobs of
//! `rusolve options`). Tolerances, limits, output and runtime options are refused with the
//! reason: a profile can make the solver
//! slower, never wrong. Every value is checked against the registry type and range
//! (`highs_option_records.rs`), and an unknown name gets the nearest-name hint
//! ([`registry::unknown_option_hint`]), the same one `--opt` and `--options-file` give.
//!
//! ## Exit codes
//!
//! 0 done; 2 the profile (or a flag) is invalid, the same code `--options-file` uses for an
//! invalid value; 3 a file is missing, unreadable or not JSON; 4 internal.
//!
//! ## Fingerprint
//!
//! The per-feature centroid, spread, min and max of the training instances' static features
//! (`rusolve classify --format json`, `FEATURE_NAMES`) and the radius, with the validity check
//! of [`Fingerprint::verdict`]: the distance is
//! `sqrt(mean(((f - centroid) / max(spread, 0.02))^2))`. The verdict is reported by
//! `rusolve mip --profile` on stderr; it never changes what is applied (the design applies in
//! every case and warns).

use rusolve_ekk::lp_data::highs_option_records::{self as registry, Affects, OptionValue};
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_run::Highs;
use serde_json::{Map, Value, json};
use std::collections::BTreeMap;
use std::process::ExitCode;

pub const PROFILE_SCHEMA: &str = "rusolve.profile/1";
pub const DIFF_SCHEMA: &str = "rusolve.profile-diff/1";
const USAGE: &str = "usage: rusolve profile diff <a> <b> [--format json]\n       \
     rusolve profile show <file|default|recommended> [--format json]\n       \
     rusolve profile validate <file>\n       \
     rusolve profile new [--base default|recommended] [--name N] [--opt key=value]... [-o file]\n  \
     a, b: a profile file, `default` (HiGHS 1.15.1 settings) or `recommended`";
const TOP_KEYS: [&str; 7] = [
    "schema_version",
    "name",
    "base",
    "options",
    "fingerprint",
    "provenance",
    "evidence",
];

/// Options the registry files under `search` that nevertheless change what is asked, when to
/// stop, or what is written: a profile may not set them (mip-profiles.md section 5: limits,
/// tolerances and anything that changes the answer or its reporting stay with the caller).
const NOT_FOR_PROFILES: &[&str] = &[
    "time_limit",
    "iis_time_limit",
    "simplex_iteration_limit",
    "ipm_iteration_limit",
    "pdlp_iteration_limit",
    "qp_iteration_limit",
    "objective_bound",
    "objective_target",
    "solve_relaxation",
    "mip_root_presolve_only",
    "mip_max_improving_sols",
    "mip_improving_solution_save",
    "mip_improving_solution_report_sparse",
    "mip_min_logging_interval",
    "mip_report_level",
    "timeless_log",
    "presolve_rule_logging",
    "user_objective_scale",
    "user_bound_scale",
    "allow_unbounded_or_infeasible",
    "infinite_bound",
];

/// The options a profile is a delta on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Base {
    Default,
    Recommended,
}

impl Base {
    pub fn name(self) -> &'static str {
        match self {
            Base::Default => "default",
            Base::Recommended => "recommended",
        }
    }
    pub fn options(self) -> HighsOptions {
        match self {
            Base::Default => HighsOptions::default(),
            Base::Recommended => HighsOptions::recommended(),
        }
    }
    pub fn parse(s: &str) -> Option<Base> {
        match s {
            "default" => Some(Base::Default),
            "recommended" => Some(Base::Recommended),
            _ => None,
        }
    }
}

/// Feature summary of the training instances (mip-profiles.md section 1).
#[derive(Clone, Debug, PartialEq)]
pub struct Fingerprint {
    pub features_version: u32,
    pub features: Vec<String>,
    pub centroid: Vec<f64>,
    pub spread: Vec<f64>,
    pub min: Vec<f64>,
    pub max: Vec<f64>,
    pub radius: f64,
    pub instances: usize,
}

/// How an instance sits relative to a fingerprint (mip-profiles.md section 3).
#[derive(Clone, Debug, PartialEq)]
pub struct Verdict {
    /// `in_region`, `near_region`, `out_of_region` or `incompatible_fingerprint`.
    pub kind: &'static str,
    pub distance: f64,
    /// The three features that contribute most to the distance: (name, contribution).
    pub top: Vec<(String, f64)>,
}

/// Standardised distance with the 0.02 floor on the spread (mip-profiles.md section 3).
fn std_distance(f: &[f64], centroid: &[f64], spread: &[f64]) -> (f64, Vec<f64>) {
    let terms: Vec<f64> = f
        .iter()
        .zip(centroid)
        .zip(spread)
        .map(|((x, c), s)| ((x - c) / s.max(0.02)).powi(2))
        .collect();
    let mean = terms.iter().sum::<f64>() / terms.len().max(1) as f64;
    (mean.sqrt(), terms)
}

impl Fingerprint {
    /// Summarise the feature vectors of the training instances.
    pub fn from_vectors(names: &[&str], vectors: &[Vec<f64>], features_version: u32) -> Self {
        let n = names.len();
        let k = vectors.len().max(1) as f64;
        let mut centroid = vec![0.0; n];
        let mut min = vec![f64::INFINITY; n];
        let mut max = vec![f64::NEG_INFINITY; n];
        for v in vectors {
            for i in 0..n {
                centroid[i] += v[i] / k;
                min[i] = min[i].min(v[i]);
                max[i] = max[i].max(v[i]);
            }
        }
        let spread: Vec<f64> = (0..n)
            .map(|i| {
                (vectors
                    .iter()
                    .map(|v| (v[i] - centroid[i]).powi(2))
                    .sum::<f64>()
                    / k)
                    .sqrt()
            })
            .collect();
        let radius = vectors
            .iter()
            .map(|v| std_distance(v, &centroid, &spread).0)
            .fold(0.0, f64::max);
        Fingerprint {
            features_version,
            features: names.iter().map(|s| s.to_string()).collect(),
            centroid,
            spread,
            min,
            max,
            radius,
            instances: vectors.len(),
        }
    }

    pub fn verdict(&self, f: &[f64], current_version: u32) -> Verdict {
        if self.features_version != current_version || f.len() != self.features.len() {
            return Verdict {
                kind: "incompatible_fingerprint",
                distance: f64::NAN,
                top: vec![],
            };
        }
        let (d, terms) = std_distance(f, &self.centroid, &self.spread);
        let mut order: Vec<usize> = (0..terms.len()).collect();
        order.sort_by(|&a, &b| terms[b].total_cmp(&terms[a]));
        let top = order
            .iter()
            .take(3)
            .map(|&i| (self.features[i].clone(), terms[i]))
            .collect();
        let inside = f
            .iter()
            .enumerate()
            .all(|(i, x)| *x >= self.min[i] - 0.05 && *x <= self.max[i] + 0.05);
        let kind = if d <= self.radius && inside {
            "in_region"
        } else if d <= 2.0 * self.radius {
            "near_region"
        } else {
            "out_of_region"
        };
        Verdict {
            kind,
            distance: d,
            top,
        }
    }

    pub fn to_value(&self) -> Value {
        json!({
            "features_version": self.features_version,
            "features": self.features,
            "centroid": self.centroid, "spread": self.spread,
            "min": self.min, "max": self.max,
            "radius": self.radius, "instances": self.instances,
        })
    }

    fn from_value(v: &Value, errors: &mut Vec<String>) -> Option<Fingerprint> {
        let obj = v.as_object().or_else(|| {
            errors.push("fingerprint: expected an object".into());
            None
        })?;
        let nums = |key: &str, errors: &mut Vec<String>| -> Vec<f64> {
            match obj.get(key).and_then(Value::as_array) {
                Some(a) => a.iter().filter_map(Value::as_f64).collect(),
                None => {
                    errors.push(format!("fingerprint.{key}: expected an array of numbers"));
                    vec![]
                }
            }
        };
        let features: Vec<String> = obj
            .get("features")
            .and_then(Value::as_array)
            .map(|a| {
                a.iter()
                    .filter_map(|s| s.as_str().map(str::to_string))
                    .collect()
            })
            .unwrap_or_else(|| {
                errors.push("fingerprint.features: expected an array of names".into());
                vec![]
            });
        let (centroid, spread) = (nums("centroid", errors), nums("spread", errors));
        let (min, max) = (nums("min", errors), nums("max", errors));
        for (k, len) in [
            ("centroid", centroid.len()),
            ("spread", spread.len()),
            ("min", min.len()),
            ("max", max.len()),
        ] {
            if len != features.len() {
                errors.push(format!(
                    "fingerprint.{k}: {len} values for {} features",
                    features.len()
                ));
            }
        }
        let version = obj.get("features_version").and_then(Value::as_u64);
        if version.is_none() {
            errors.push("fingerprint.features_version: expected an integer".into());
        }
        Some(Fingerprint {
            features_version: version.unwrap_or(0) as u32,
            features,
            centroid,
            spread,
            min,
            max,
            radius: obj.get("radius").and_then(Value::as_f64).unwrap_or(0.0),
            instances: obj.get("instances").and_then(Value::as_u64).unwrap_or(0) as usize,
        })
    }
}

/// A parsed, validated profile.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile {
    pub name: Option<String>,
    pub base: Base,
    /// The options that differ from the base, sorted by name.
    pub options: BTreeMap<String, OptionValue>,
    pub fingerprint: Option<Fingerprint>,
    pub provenance: Option<Value>,
    pub evidence: Option<Value>,
}

pub fn value_json(v: &OptionValue) -> Value {
    match v {
        OptionValue::Bool(b) => json!(b),
        OptionValue::Int(i) => json!(i),
        OptionValue::Double(d) if d.is_finite() => json!(d),
        OptionValue::Double(d) => json!(if *d > 0.0 { "inf" } else { "-inf" }),
        OptionValue::Str(s) => json!(s),
    }
}

fn record(name: &str) -> Option<&'static registry::OptionRecord> {
    registry::OPTION_RECORDS.iter().find(|r| r.name == name)
}

fn affects_word(a: Affects) -> &'static str {
    match a {
        Affects::Search => "search",
        Affects::Tolerance => "tolerance",
        Affects::Semantics => "semantics",
        Affects::Output => "output",
        Affects::Runtime => "runtime",
    }
}

/// Check one JSON value against the registry (name, search class, type, range). The error
/// text starts with the JSON path (`options.<name>`), so an agent can find the line.
pub fn validate_option(name: &str, v: &Value) -> Result<OptionValue, String> {
    let path = format!("options.{name}");
    let Some(rec) = record(name) else {
        return Err(match registry::unknown_option_hint(name) {
            Some(h) => format!("{path}: unknown option; {h}"),
            None => format!("{path}: unknown option (`rusolve options` lists every name)"),
        });
    };
    if NOT_FOR_PROFILES.contains(&name) {
        return Err(format!(
            "{path}: a profile may not set this option: it is a limit, a target or an output \
             setting, not a search strategy (set it with --opt on the command line instead)"
        ));
    }
    if rec.affects != Affects::Search {
        return Err(format!(
            "{path}: a profile may only set search options, and this one affects '{}' \
             (set it with --opt on the command line instead)",
            affects_word(rec.affects)
        ));
    }
    let ty = registry::option_type_name(name).unwrap_or("?");
    match ty {
        "bool" => v
            .as_bool()
            .map(OptionValue::Bool)
            .ok_or_else(|| format!("{path}: expected a JSON boolean (true/false), got {v}")),
        "HighsInt" => {
            let i = v
                .as_i64()
                .ok_or_else(|| format!("{path}: expected an integer, got {v}"))?;
            let (lo, hi) = registry::option_int_bounds(name).unwrap_or((i32::MIN, i32::MAX));
            if i < lo as i64 || i > hi as i64 {
                return Err(format!("{path}: {i} is outside the range [{lo}, {hi}]"));
            }
            Ok(OptionValue::Int(i as i32))
        }
        "double" => {
            let d = match v {
                Value::String(s) if s == "inf" => f64::INFINITY,
                Value::String(s) if s == "-inf" => f64::NEG_INFINITY,
                _ => v
                    .as_f64()
                    .ok_or_else(|| format!("{path}: expected a number, got {v}"))?,
            };
            let (lo, hi) = registry::option_double_bounds(name).unwrap_or((f64::MIN, f64::MAX));
            if d < lo || d > hi || d.is_nan() {
                return Err(format!("{path}: {d} is outside the range [{lo}, {hi}]"));
            }
            Ok(OptionValue::Double(d))
        }
        _ => v
            .as_str()
            .map(|s| OptionValue::Str(s.to_string()))
            .ok_or_else(|| format!("{path}: expected a string, got {v}")),
    }
}

/// Parse a string (`--opt k=v`, a candidate list) into a typed value, then validate it.
pub fn parse_option_text(name: &str, text: &str) -> Result<OptionValue, String> {
    let v = match registry::option_type_name(name) {
        Some("bool") => match text {
            "true" => json!(true),
            "false" => json!(false),
            _ => {
                return Err(format!(
                    "options.{name}: expected true or false, got '{text}'"
                ));
            }
        },
        Some("HighsInt") => text
            .parse::<i64>()
            .map(|i| json!(i))
            .map_err(|_| format!("options.{name}: expected an integer, got '{text}'"))?,
        Some("double") => match text {
            "inf" | "-inf" => json!(text),
            _ => text
                .parse::<f64>()
                .map(|d| json!(d))
                .map_err(|_| format!("options.{name}: expected a number, got '{text}'"))?,
        },
        _ => json!(text),
    };
    validate_option(name, &v)
}

impl Profile {
    pub fn new(base: Base) -> Self {
        Profile {
            name: None,
            base,
            options: BTreeMap::new(),
            fingerprint: None,
            provenance: None,
            evidence: None,
        }
    }

    /// Validate a whole file; every problem is reported, not only the first.
    pub fn from_value(v: &Value) -> Result<Profile, Vec<String>> {
        let mut errors = Vec::new();
        let Some(obj) = v.as_object() else {
            return Err(vec!["the profile must be a JSON object".into()]);
        };
        match obj.get("schema_version").and_then(Value::as_str) {
            Some(PROFILE_SCHEMA) => {}
            Some(other) => errors.push(format!(
                "schema_version: '{other}' is not supported; this binary reads '{PROFILE_SCHEMA}'"
            )),
            None => errors.push(format!(
                "schema_version: missing; expected \"{PROFILE_SCHEMA}\""
            )),
        }
        for key in obj.keys() {
            if !TOP_KEYS.contains(&key.as_str()) {
                errors.push(format!(
                    "{key}: unknown key; the keys of {PROFILE_SCHEMA} are {}",
                    TOP_KEYS.join(", ")
                ));
            }
        }
        let base = match obj.get("base").and_then(Value::as_str) {
            None => Base::Recommended,
            Some(s) => Base::parse(s).unwrap_or_else(|| {
                errors.push(format!("base: '{s}' is not default or recommended"));
                Base::Recommended
            }),
        };
        let mut profile = Profile::new(base);
        profile.name = obj.get("name").and_then(Value::as_str).map(str::to_string);
        match obj.get("options") {
            None => {}
            Some(Value::Object(o)) => {
                for (k, val) in o {
                    match validate_option(k, val) {
                        Ok(ov) => {
                            profile.options.insert(k.clone(), ov);
                        }
                        Err(e) => errors.push(e),
                    }
                }
            }
            Some(_) => errors.push("options: expected an object of name: value".into()),
        }
        if let Some(f) = obj.get("fingerprint").filter(|f| !f.is_null()) {
            profile.fingerprint = Fingerprint::from_value(f, &mut errors);
        }
        profile.provenance = obj.get("provenance").cloned().filter(|p| !p.is_null());
        profile.evidence = obj.get("evidence").cloned().filter(|p| !p.is_null());
        if errors.is_empty() {
            // The typed setter is the last word: a value the registry check passed but the
            // port refuses (for instance an enum-like string) must fail here, not in a solve.
            if let Err(e) = profile.effective() {
                errors.push(e);
            }
        }
        if errors.is_empty() {
            Ok(profile)
        } else {
            Err(errors)
        }
    }

    pub fn to_value(&self) -> Value {
        let mut o = Map::new();
        o.insert("schema_version".into(), json!(PROFILE_SCHEMA));
        if let Some(n) = &self.name {
            o.insert("name".into(), json!(n));
        }
        o.insert("base".into(), json!(self.base.name()));
        let options: Map<String, Value> = self
            .options
            .iter()
            .map(|(k, v)| (k.clone(), value_json(v)))
            .collect();
        o.insert("options".into(), Value::Object(options));
        if let Some(f) = &self.fingerprint {
            o.insert("fingerprint".into(), f.to_value());
        }
        if let Some(p) = &self.provenance {
            o.insert("provenance".into(), p.clone());
        }
        if let Some(e) = &self.evidence {
            o.insert("evidence".into(), e.clone());
        }
        Value::Object(o)
    }

    /// The options a solve starts from: the base with the profile's options applied through
    /// the same typed setter `--opt` uses.
    pub fn effective(&self) -> Result<HighsOptions, String> {
        let mut base = self.base.options();
        let (output_flag, console) = (base.output_flag, base.log_to_console);
        base.output_flag = false;
        let mut highs = Highs::with_options(base);
        for (name, value) in &self.options {
            let status = highs.set_option_value_str(name, &value.to_token());
            if status == rusolve_ekk::prelude::HighsStatus::kError {
                return Err(format!(
                    "options.{name}: the port refused the value {}",
                    value.to_token()
                ));
            }
        }
        let mut out = highs.get_options().clone();
        out.output_flag = output_flag;
        out.log_to_console = console;
        out.sync_log_options();
        Ok(out)
    }

    /// Record `effective` as a delta on `base`: every registry option whose value differs.
    pub fn from_options(base: Base, effective: &HighsOptions) -> Profile {
        let mut p = Profile::new(base);
        for d in diff_options(&base.options(), effective) {
            p.options.insert(d.name.to_string(), d.b);
        }
        p
    }

    pub fn read(path: &str) -> Result<Profile, (u8, String)> {
        let text =
            std::fs::read_to_string(path).map_err(|e| (3, format!("{path}: cannot read: {e}")))?;
        let v: Value =
            serde_json::from_str(&text).map_err(|e| (3, format!("{path}: not JSON: {e}")))?;
        Profile::from_value(&v).map_err(|errs| {
            (
                2,
                crate::diag::INVALID_CONTENT.to_string()
                    + &errs
                        .iter()
                        .map(|e| format!("{path}: {e}"))
                        .collect::<Vec<_>>()
                        .join("\n"),
            )
        })
    }
}

/// One option that differs between two settings.
pub struct DiffEntry {
    pub name: &'static str,
    pub a: OptionValue,
    pub b: OptionValue,
    pub affects: &'static str,
    pub source: &'static str,
}

/// Every registry option whose value differs between `a` and `b`, in registry order.
pub fn diff_options(a: &HighsOptions, b: &HighsOptions) -> Vec<DiffEntry> {
    registry::OPTION_RECORDS
        .iter()
        .filter_map(|r| {
            let (va, vb) = (
                registry::get_option_value(a, r.name)?,
                registry::get_option_value(b, r.name)?,
            );
            (va != vb).then_some(DiffEntry {
                name: r.name,
                a: va,
                b: vb,
                affects: affects_word(r.affects),
                source: r.source,
            })
        })
        .collect()
}

/// `default`, `recommended`, or a profile file.
fn resolve(spec: &str) -> Result<(String, HighsOptions, Option<Profile>), (u8, String)> {
    if let Some(b) = Base::parse(spec) {
        return Ok((spec.to_string(), b.options(), None));
    }
    let p = Profile::read(spec)?;
    let eff = p
        .effective()
        .map_err(|e| (2, format!("{}{spec}: {e}", crate::diag::INVALID_CONTENT)))?;
    Ok((spec.to_string(), eff, Some(p)))
}

/// Static feature vector of a model, in `FEATURE_NAMES` order.
pub fn feature_vector(path: &str) -> Result<Vec<f64>, String> {
    let model = rusolve_io::read_model_file(path).map_err(|e| format!("read {path}: {e}"))?;
    model.gc.reject("rusolve profile")?;
    let problem = rusolve_mip::MipProblem::from_model(&model);
    let cl = rusolve_mip::classify::classify(&problem, Default::default());
    Ok(cl.features().to_vec())
}

/// The verdict line `rusolve mip --profile` prints on stderr. `None` when the profile has no
/// fingerprint or the model cannot be featurised (the solve reports its own read error).
pub fn region_message(profile: &Profile, model_path: &str) -> Option<String> {
    let fp = profile.fingerprint.as_ref()?;
    let f = feature_vector(model_path).ok()?;
    let v = fp.verdict(&f, rusolve_mip::classify::FEATURE_VERSION);
    let top = v
        .top
        .iter()
        .map(|(n, c)| format!("{n} ({c:.2})"))
        .collect::<Vec<_>>()
        .join(", ");
    Some(format!(
        "rusolve mip: profile {}: {} (distance {:.3}, radius {:.3}){}",
        profile.name.as_deref().unwrap_or("(unnamed)"),
        v.kind,
        v.distance,
        fp.radius,
        if top.is_empty() {
            String::new()
        } else {
            format!("; largest contributions: {top}")
        }
    ))
}

fn fail(json: bool, code: u8, msg: &str) -> ExitCode {
    crate::diag::fail_coded(json, "profile", code, msg)
}

fn show_value(v: &OptionValue) -> String {
    v.to_token()
}

pub fn run(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args) || args.iter().any(|a| a == "--json");
    let Some(sub) = args.first() else {
        return fail(wj, 2, &format!("missing subcommand\n{USAGE}"));
    };
    let mut json_out = false;
    let mut positional: Vec<&String> = Vec::new();
    let mut base = Base::Recommended;
    let mut name: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut sets: Vec<String> = Vec::new();
    let mut it = args[1..].iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--format" => match it.next().map(String::as_str) {
                Some("json") => json_out = true,
                Some("text") => {}
                _ => return fail(wj, 2, "--format needs json or text"),
            },
            "--json" => json_out = true,
            "--base" => match it.next().and_then(|s| Base::parse(s)) {
                Some(b) => base = b,
                None => return fail(wj, 2, "--base needs default or recommended"),
            },
            "--name" => match it.next() {
                Some(n) => name = Some(n.clone()),
                None => return fail(wj, 2, "--name needs a value"),
            },
            "-o" | "--output" => match it.next() {
                Some(p) => out_path = Some(p.clone()),
                None => return fail(wj, 2, "-o needs a path"),
            },
            "--opt" => match it.next() {
                Some(kv) => sets.push(kv.clone()),
                None => return fail(wj, 2, "--opt needs key=value"),
            },
            "--help" | "-h" => {
                println!("{USAGE}");
                return ExitCode::SUCCESS;
            }
            s if s.starts_with("--") => {
                return fail(wj, 2, &format!("flag '{s}' is not supported\n{USAGE}"));
            }
            _ => positional.push(a),
        }
    }
    match sub.as_str() {
        "diff" => {
            let [a, b] = positional.as_slice() else {
                return fail(
                    wj,
                    2,
                    &format!("missing or extra argument for `profile diff`\n{USAGE}"),
                );
            };
            let ((an, ao, _), (bn, bo, _)) = match (resolve(a), resolve(b)) {
                (Ok(x), Ok(y)) => (x, y),
                (Err((c, m)), _) | (_, Err((c, m))) => return fail(wj, c, &m),
            };
            let diffs = diff_options(&ao, &bo);
            if json_out {
                let entries: Vec<Value> = diffs
                    .iter()
                    .map(|d| {
                        json!({"name": d.name, "a": value_json(&d.a), "b": value_json(&d.b),
                               "affects": d.affects, "source": d.source})
                    })
                    .collect();
                println!(
                    "{}",
                    json!({"schema_version": DIFF_SCHEMA, "a": an, "b": bn,
                           "identical": diffs.is_empty(), "differences": entries})
                );
            } else if diffs.is_empty() {
                println!("{an} and {bn} set every option alike");
            } else {
                println!("{an} -> {bn}: {} option(s) differ", diffs.len());
                for d in &diffs {
                    println!(
                        "  {} : {} -> {}   [{}; {}]",
                        d.name,
                        show_value(&d.a),
                        show_value(&d.b),
                        d.affects,
                        d.source
                    );
                }
            }
            ExitCode::SUCCESS
        }
        "show" => {
            let [a] = positional.as_slice() else {
                return fail(
                    wj,
                    2,
                    &format!("missing or extra argument for `profile show`\n{USAGE}"),
                );
            };
            let (n, eff, p) = match resolve(a) {
                Ok(x) => x,
                Err((c, m)) => return fail(wj, c, &m),
            };
            let profile = p.unwrap_or_else(|| {
                Profile::from_options(Base::parse(&n).unwrap_or(Base::Recommended), &eff)
            });
            if json_out {
                println!("{}", profile.to_value());
            } else {
                println!(
                    "profile {} (base {}): {} option(s) changed",
                    profile.name.as_deref().unwrap_or(&n),
                    profile.base.name(),
                    profile.options.len()
                );
                for (k, v) in &profile.options {
                    println!("  {k} = {}", show_value(v));
                }
                if let Some(f) = &profile.fingerprint {
                    println!(
                        "fingerprint: {} instances, features_version {}, radius {:.3}",
                        f.instances, f.features_version, f.radius
                    );
                }
                if let Some(e) = &profile.evidence {
                    println!("evidence: {}", e);
                }
            }
            ExitCode::SUCCESS
        }
        "validate" => {
            let [a] = positional.as_slice() else {
                return fail(
                    wj,
                    2,
                    &format!("missing or extra argument for `profile validate`\n{USAGE}"),
                );
            };
            match Profile::read(a) {
                Ok(p) => {
                    println!("{a}: valid {PROFILE_SCHEMA}, {} option(s)", p.options.len());
                    ExitCode::SUCCESS
                }
                Err((c, m)) => fail(wj, c, &m),
            }
        }
        "new" => {
            let mut p = Profile::new(base);
            p.name = name;
            for kv in &sets {
                let Some((k, v)) = kv.split_once('=') else {
                    return fail(wj, 2, &format!("--opt {kv}: expected key=value"));
                };
                match parse_option_text(k, v) {
                    Ok(ov) => {
                        p.options.insert(k.to_string(), ov);
                    }
                    Err(e) => return fail(wj, 2, &e),
                }
            }
            let text = serde_json::to_string_pretty(&p.to_value()).expect("profile serialises");
            match out_path {
                Some(path) => {
                    if let Err(e) = std::fs::write(&path, text + "\n") {
                        return fail(wj, 4, &format!("cannot write {path}: {e}"));
                    }
                }
                None => println!("{text}"),
            }
            ExitCode::SUCCESS
        }
        other => fail(wj, 2, &format!("unknown subcommand '{other}'\n{USAGE}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_through_json_keeps_every_option() {
        let mut p = Profile::new(Base::Recommended);
        p.name = Some("rt".into());
        p.options
            .insert("mip_heuristic_effort".into(), OptionValue::Double(0.15));
        p.options
            .insert("mip_pscost_minreliable".into(), OptionValue::Int(16));
        p.options
            .insert("mip_detect_symmetry".into(), OptionValue::Bool(false));
        let text = serde_json::to_string(&p.to_value()).unwrap();
        let back = Profile::from_value(&serde_json::from_str(&text).unwrap()).unwrap();
        assert_eq!(back, p);
        let eff = back.effective().unwrap();
        assert_eq!(eff.mip_heuristic_effort, 0.15);
        assert_eq!(eff.mip_pscost_minreliable, 16);
        assert!(!eff.mip_detect_symmetry);
        // writing the effective options back as a delta gives the same profile
        let again = Profile::from_options(Base::Recommended, &eff);
        assert_eq!(again.options, p.options);
    }

    #[test]
    fn invalid_values_get_paths_and_hints() {
        let v = json!({"schema_version": PROFILE_SCHEMA, "options": {
            "mip_pscost_minrelaible": 4,
            "mip_heuristic_effort": 7.0,
            "mip_detect_symmetry": "yes",
            "time_limit": 5.0}});
        let errs = Profile::from_value(&v).unwrap_err();
        let all = errs.join("\n");
        assert!(all.contains("options.mip_pscost_minrelaible: unknown option; did you mean mip_pscost_minreliable?"), "{all}");
        assert!(
            all.contains("options.mip_heuristic_effort: 7 is outside the range"),
            "{all}"
        );
        assert!(
            all.contains("options.mip_detect_symmetry: expected a JSON boolean"),
            "{all}"
        );
        assert!(
            all.contains("options.time_limit: a profile may not set this option"),
            "{all}"
        );
    }

    #[test]
    fn wrong_schema_version_and_unknown_key_fail() {
        let errs = Profile::from_value(&json!({"schema_version": "rusolve.profile/2", "opts": {}}))
            .unwrap_err();
        assert!(errs.iter().any(|e| e.starts_with("schema_version")));
        assert!(errs.iter().any(|e| e.starts_with("opts: unknown key")));
    }

    #[test]
    fn fingerprint_verdicts() {
        let names = ["a", "b"];
        let train = vec![vec![0.2, 0.5], vec![0.4, 0.5], vec![0.3, 0.5]];
        let fp = Fingerprint::from_vectors(&names, &train, 1);
        assert_eq!(fp.verdict(&[0.3, 0.5], 1).kind, "in_region");
        assert_eq!(fp.verdict(&[0.95, 0.5], 1).kind, "out_of_region");
        assert_eq!(fp.verdict(&[0.3, 0.5], 2).kind, "incompatible_fingerprint");
    }
}
