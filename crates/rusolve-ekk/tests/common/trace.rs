//! Simplex traces (`simplex/<lp>.<dual|primal>.jsonl`): listing, summaries,
//! and the lockstep first-divergence reporter.
//!
//! A port produces its own trace as a `Vec<serde_json::Value>`, one record
//! per HiGHS hook, built with [`super::events`] (or [`super::emit::Obj`]).
//! [`compare_traces`] walks both streams together and stops at the first
//! record that differs, reporting:
//! - record index, fixture line, and type (expected vs actual);
//! - the last iteration, solve phase, and rebuild seen before it;
//! - every field mismatch in that record (bit patterns + decimals, hashes);
//! - the preceding records for context;
//! - a hint naming the HiGHS function(s) most likely responsible.
//!
//! Iteration counts must match exactly: reinversion is triggered by the
//! synthetic clock (`total_synthetic_tick_ >= build_synthetic_tick_`), so a
//! wrong tick changes the record sequence (a `rebuild` where an `iter` was
//! expected), and that is reported as a type mismatch.

use std::fmt;
use std::fmt::Write as _;

use serde_json::Value;

use super::check::{Mismatch, diff_values, fmt_f64};
use super::fixture::{Fixture, Location, manifest, record_type};
use super::hex64;

/// Which simplex solver produced a trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraceKind {
    /// `simplex_strategy = 1` (`HEkkDual`)
    Dual,
    /// `simplex_strategy = 4` (`HEkkPrimal`)
    Primal,
}

impl TraceKind {
    pub fn suffix(self) -> &'static str {
        match self {
            TraceKind::Dual => ".dual.jsonl",
            TraceKind::Primal => ".primal.jsonl",
        }
    }

    /// The `simplex_strategy` option value the oracle ran with.
    pub fn simplex_strategy(self) -> i32 {
        match self {
            TraceKind::Dual => 1,
            TraceKind::Primal => 4,
        }
    }
}

/// Every trace of `kind` in the manifest (committed first, then
/// target-only), as manifest names for [`super::fixture::load`].
pub fn trace_files(kind: TraceKind) -> Vec<String> {
    let mut v: Vec<_> = manifest().select("simplex/", kind.suffix());
    v.sort_by_key(|e| (e.location != Location::Committed, e.file.clone()));
    v.into_iter().map(|e| e.file.clone()).collect()
}

/// Only the committed traces of `kind` (always present).
pub fn committed_trace_files(kind: TraceKind) -> Vec<String> {
    manifest()
        .select("simplex/", kind.suffix())
        .into_iter()
        .filter(|e| e.location == Location::Committed)
        .map(|e| e.file.clone())
        .collect()
}

/// LP name of a trace file: `simplex/afiro.dual.jsonl` -> `afiro`.
pub fn trace_lp_name(file: &str) -> &str {
    let base = file.rsplit('/').next().unwrap_or(file);
    base.split('.').next().unwrap_or(base)
}

/// Counts per record kind plus the final outcome, for logs and sanity checks.
#[derive(Clone, Debug, Default)]
pub struct TraceSummary {
    pub file: String,
    pub records: usize,
    pub iters: usize,
    pub iter_none: usize,
    pub inverts: usize,
    pub rebuilds: usize,
    /// `(rebuild reason, count)` sorted by reason.
    pub rebuild_reasons: Vec<(i32, usize)>,
    pub cleanups: usize,
    pub phase_ends: usize,
    /// `iteration` of the last `iter`.
    pub last_iteration: Option<i32>,
    pub model_status: Option<i32>,
    pub objective: Option<f64>,
    pub simplex_iteration_count: Option<i32>,
}

impl fmt::Display for TraceSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: {} records, {} iter, {} iter_none, {} invert, {} rebuild {:?}, {} cleanup, {} phase_end, last iteration {:?}, model_status {:?}, objective {}, simplex_iteration_count {:?}",
            self.file,
            self.records,
            self.iters,
            self.iter_none,
            self.inverts,
            self.rebuilds,
            self.rebuild_reasons,
            self.cleanups,
            self.phase_ends,
            self.last_iteration,
            self.model_status,
            self.objective.map(fmt_f64).unwrap_or_else(|| "None".into()),
            self.simplex_iteration_count,
        )
    }
}

pub fn summarize(fx: &Fixture) -> TraceSummary {
    let mut s = TraceSummary {
        file: fx.file.clone(),
        records: fx.len(),
        ..Default::default()
    };
    for r in fx.recs() {
        match r.t() {
            "iter" => {
                s.iters += 1;
                s.last_iteration = Some(r.int("iteration"));
            }
            "iter_none" => s.iter_none += 1,
            "invert" => s.inverts += 1,
            "rebuild" => {
                s.rebuilds += 1;
                let reason = r.int("reason");
                match s.rebuild_reasons.iter_mut().find(|(k, _)| *k == reason) {
                    Some((_, c)) => *c += 1,
                    None => s.rebuild_reasons.push((reason, 1)),
                }
            }
            "cleanup" => s.cleanups += 1,
            "phase_end" => s.phase_ends += 1,
            "result" => {
                s.model_status = Some(r.int("model_status"));
                s.objective = Some(r.f64("objective_function_value"));
                s.simplex_iteration_count = Some(r.int("simplex_iteration_count"));
            }
            _ => {}
        }
    }
    s.rebuild_reasons.sort();
    s
}

/// Options for [`compare_traces`].
#[derive(Clone, Debug)]
pub struct TraceCompare {
    /// Fields removed from both sides before diffing, as `type/dotted.path`
    /// (`*` matches every type), e.g. `ekk_solve_enter/options`,
    /// `*/scalars.cost_perturbation_base`. Use only for fields the component
    /// under test does not produce yet, and say so in the test.
    pub skip: Vec<String>,
    /// When set, only records of these types are compared (both streams are
    /// filtered first). Lets a test pin e.g. just `iter` + `rebuild` while
    /// other hooks are not wired. Positions in reports are still the
    /// original fixture lines.
    pub only_types: Option<Vec<String>>,
    /// Stop after this many compared records (for staged bring-up).
    pub max_records: Option<usize>,
    /// How many preceding records to print as context.
    pub context: usize,
}

impl Default for TraceCompare {
    fn default() -> Self {
        TraceCompare {
            skip: Vec::new(),
            only_types: None,
            max_records: None,
            context: 6,
        }
    }
}

impl TraceCompare {
    pub fn skip(mut self, spec: &str) -> Self {
        assert!(
            spec.contains('/'),
            "skip spec must be \"type/dotted.path\", got {spec:?}"
        );
        self.skip.push(spec.to_string());
        self
    }

    pub fn only(mut self, types: &[&str]) -> Self {
        self.only_types = Some(types.iter().map(|s| s.to_string()).collect());
        self
    }

    pub fn max_records(mut self, n: usize) -> Self {
        self.max_records = Some(n);
        self
    }

    fn keeps(&self, t: &str) -> bool {
        self.only_types
            .as_ref()
            .is_none_or(|v| v.iter().any(|x| x == t))
    }

    fn strip(&self, v: &Value) -> Value {
        let mut v = v.clone();
        let t = record_type(&v).to_string();
        for spec in &self.skip {
            let (ty, path) = spec.split_once('/').expect("checked in skip()");
            if ty != "*" && ty != t {
                continue;
            }
            remove_path(&mut v, path);
        }
        v
    }
}

fn remove_path(v: &mut Value, path: &str) {
    let mut parts = path.split('.').peekable();
    let mut cur = v;
    while let Some(p) = parts.next() {
        let Some(obj) = cur.as_object_mut() else {
            return;
        };
        if parts.peek().is_none() {
            obj.remove(p);
            return;
        }
        match obj.get_mut(p) {
            Some(next) => cur = next,
            None => return,
        }
    }
}

/// What kind of divergence was found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DivergenceKind {
    /// Same record type, different field values.
    Fields,
    /// Different record type at the same position (e.g. the port rebuilt
    /// where HiGHS iterated): the control flow diverged.
    Type,
    /// The port's trace has fewer records.
    ActualEnded,
    /// The port's trace has more records.
    ExpectedEnded,
}

/// The first point where the port's trace leaves the fixture.
#[derive(Clone, Debug)]
pub struct TraceDivergence {
    pub file: String,
    pub kind: DivergenceKind,
    /// 0-based position among compared records.
    pub position: usize,
    /// 1-based fixture line of the expected record (`None` if the fixture ended).
    pub line: Option<usize>,
    /// 0-based index into the actual stream (`None` if it ended).
    pub actual_index: Option<usize>,
    pub expected_type: Option<String>,
    pub actual_type: Option<String>,
    /// `iteration` of the last matching `iter`, if any.
    pub last_iteration: Option<i64>,
    /// Last `solve_phase` seen in matching records.
    pub last_phase: Option<i64>,
    /// `(fixture line, reason)` of the last matching `rebuild`.
    pub last_rebuild: Option<(usize, i64)>,
    pub mismatches: Vec<Mismatch>,
    /// One-line renderings of the preceding expected records.
    pub context: Vec<String>,
    /// HiGHS functions most likely responsible.
    pub hints: Vec<String>,
}

impl fmt::Display for TraceDivergence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let what = match self.kind {
            DivergenceKind::Fields => "field mismatch",
            DivergenceKind::Type => "record type mismatch (control flow diverged)",
            DivergenceKind::ActualEnded => "actual trace ended early",
            DivergenceKind::ExpectedEnded => "actual trace has extra records",
        };
        writeln!(f, "trace {} diverges: {what}", self.file)?;
        writeln!(
            f,
            "  at compared record #{} (fixture line {}, actual index {})",
            self.position,
            self.line.map_or("<end>".into(), |l| l.to_string()),
            self.actual_index.map_or("<end>".into(), |i| i.to_string()),
        )?;
        writeln!(
            f,
            "  expected type {:?}, actual type {:?}",
            self.expected_type.as_deref().unwrap_or("<none>"),
            self.actual_type.as_deref().unwrap_or("<none>")
        )?;
        writeln!(
            f,
            "  last matching: iteration {:?}, solve_phase {:?}, rebuild {}",
            self.last_iteration,
            self.last_phase,
            self.last_rebuild.map_or("<none>".into(), |(l, r)| format!(
                "at line {l} (reason {r})"
            ))
        )?;
        if !self.mismatches.is_empty() {
            writeln!(f, "  {} field mismatch(es):", self.mismatches.len())?;
            for m in &self.mismatches {
                writeln!(f, "    {}", m.to_string().replace('\n', "\n    "))?;
            }
        }
        if !self.hints.is_empty() {
            writeln!(f, "  likely responsible:")?;
            for h in &self.hints {
                writeln!(f, "    - {h}")?;
            }
        }
        if !self.context.is_empty() {
            writeln!(f, "  preceding records:")?;
            for c in &self.context {
                writeln!(f, "    {c}")?;
            }
        }
        Ok(())
    }
}

/// Compare the port's trace with a fixture trace, record by record.
pub fn compare_traces(
    expected: &Fixture,
    actual: &[Value],
    cfg: &TraceCompare,
) -> Result<(), Box<TraceDivergence>> {
    let exp: Vec<(usize, &Value)> = expected
        .records
        .iter()
        .enumerate()
        .filter(|(_, v)| cfg.keeps(record_type(v)))
        .collect();
    let act: Vec<(usize, &Value)> = actual
        .iter()
        .enumerate()
        .filter(|(_, v)| cfg.keeps(record_type(v)))
        .collect();
    let limit = cfg.max_records.unwrap_or(usize::MAX);

    let mut last_iteration = None;
    let mut last_phase = None;
    let mut last_rebuild = None;
    let mut pos = 0;
    loop {
        if pos >= limit {
            return Ok(());
        }
        let e = exp.get(pos);
        let a = act.get(pos);
        let (kind, mismatches) = match (e, a) {
            (None, None) => return Ok(()),
            (Some(_), None) => (DivergenceKind::ActualEnded, Vec::new()),
            (None, Some(_)) => (DivergenceKind::ExpectedEnded, Vec::new()),
            (Some((ei, ev)), Some((_, av))) => {
                let (et, at) = (record_type(ev), record_type(av));
                if et != at {
                    (DivergenceKind::Type, Vec::new())
                } else {
                    let ctx = format!("{}:{} [{et}]", expected.file, ei + 1);
                    let m = diff_values(&ctx, &cfg.strip(ev), &cfg.strip(av));
                    if m.is_empty() {
                        track(
                            ev,
                            *ei,
                            &mut last_iteration,
                            &mut last_phase,
                            &mut last_rebuild,
                        );
                        pos += 1;
                        continue;
                    }
                    (DivergenceKind::Fields, m)
                }
            }
        };
        let expected_type = e.map(|(_, v)| record_type(v).to_string());
        let actual_type = a.map(|(_, v)| record_type(v).to_string());
        let hints = hints_for(
            &kind,
            expected_type.as_deref(),
            actual_type.as_deref(),
            &mismatches,
            e.map(|(i, _)| format!("{}:{} [", expected.file, i + 1)),
        );
        let context = exp[pos.saturating_sub(cfg.context)..pos]
            .iter()
            .map(|(i, v)| one_line(*i, v))
            .chain(
                a.filter(|_| kind != DivergenceKind::Fields)
                    .map(|(i, v)| format!("actual[{i}]: {}", one_line_body(v))),
            )
            .collect();
        return Err(Box::new(TraceDivergence {
            file: expected.file.clone(),
            kind,
            position: pos,
            line: e.map(|(i, _)| i + 1),
            actual_index: a.map(|(i, _)| *i),
            expected_type,
            actual_type,
            last_iteration,
            last_phase,
            last_rebuild,
            mismatches,
            context,
            hints,
        }));
    }
}

/// [`compare_traces`], panicking with the full divergence report.
#[track_caller]
pub fn assert_trace_matches(expected: &Fixture, actual: &[Value], cfg: &TraceCompare) {
    if let Err(d) = compare_traces(expected, actual, cfg) {
        panic!("{d}");
    }
}

fn int_at(v: &Value, path: &[&str]) -> Option<i64> {
    let mut cur = v;
    for p in path {
        cur = cur.get(p)?;
    }
    cur.as_i64()
}

fn track(
    v: &Value,
    i: usize,
    last_iteration: &mut Option<i64>,
    last_phase: &mut Option<i64>,
    last_rebuild: &mut Option<(usize, i64)>,
) {
    let t = record_type(v);
    if t == "iter" {
        *last_iteration = int_at(v, &["iteration"]).or(*last_iteration);
    }
    if let Some(p) = int_at(v, &["solve_phase"])
        .or_else(|| int_at(v, &["p", "solve_phase"]))
        .or_else(|| int_at(v, &["state", "solve_phase"]))
    {
        *last_phase = Some(p);
    }
    if t == "rebuild" {
        *last_rebuild = int_at(v, &["reason"]).map(|r| (i + 1, r));
    }
}

/// Fields shown in one-line record renderings, by path.
const SUMMARY_FIELDS: &[&[&str]] = &[
    &["algorithm"],
    &["iteration"],
    &["reason"],
    &["solve_phase"],
    &["row_out"],
    &["variable_in"],
    &["variable_out"],
    &["p", "variable_in"],
    &["p", "row_out"],
    &["p", "variable_out"],
    &["num_flips"],
    &["update_count"],
    &["model_status"],
];

fn one_line_body(v: &Value) -> String {
    let mut s = format!("[{}]", record_type(v));
    for path in SUMMARY_FIELDS {
        let mut cur = Some(v);
        for p in *path {
            cur = cur.and_then(|c| c.get(p));
        }
        if let Some(x) = cur {
            let rendered = match x.as_str() {
                Some(st) => st.to_string(),
                None => x.to_string(),
            };
            let _ = write!(s, " {}={}", path.join("."), rendered);
        }
    }
    if let Some(h) = v.get("h_basicIndex").and_then(Value::as_str) {
        let _ = write!(s, " h_basicIndex={h}");
    }
    if let Some(x) = v
        .get("total_synthetic_tick")
        .and_then(Value::as_str)
        .and_then(|s| hex64(s).ok())
    {
        let _ = write!(s, " tick={:e}", f64::from_bits(x));
    }
    s
}

fn one_line(i: usize, v: &Value) -> String {
    format!("line {}: {}", i + 1, one_line_body(v))
}

/// `(record type or "*", field path prefix relative to the record, hint)`.
/// The first matching row wins per mismatch; distinct hints are listed.
const HINTS: &[(&str, &str, &str)] = &[
    // dual iteration
    (
        "iter",
        "row_out",
        "HEkkDualRHS::chooseNormal (CHUZR: infeasibility/weight ratio, random start via HighsRandom)",
    ),
    (
        "iter",
        "computed_edge_weight",
        "HEkkDual::chooseRow (recomputed DSE weight) / HEkk::dual_edge_weight_ from the previous iteration",
    ),
    (
        "iter",
        "variable_out",
        "HEkkDual::chooseRow / HEkkDualRHS::chooseNormal",
    ),
    (
        "iter",
        "move_out",
        "HEkkDual::chooseRow (baseValue vs baseLower/baseUpper)",
    ),
    (
        "iter",
        "delta_primal",
        "HEkkDual::chooseRow (primal infeasibility of the leaving row)",
    ),
    (
        "iter",
        "variable_in",
        "HEkkDualRow::chooseMakepack/choosePossible/chooseFinal (CHUZC + BFRT) on HEkk::tableauRowPrice output",
    ),
    (
        "iter",
        "alpha_row",
        "HEkkDualRow::chooseFinal (pivot from row_ap) / HighsSparseMatrix::priceByRow*",
    ),
    (
        "iter",
        "theta_dual",
        "HEkkDualRow::chooseFinal / chooseFinalWorkGroupQuad|Heap (HighsCDouble, maxheapsort)",
    ),
    (
        "iter",
        "num_flips",
        "HEkkDualRow::chooseFinal BFRT groups (workData/workGroup)",
    ),
    (
        "iter",
        "flips",
        "HEkkDualRow::chooseFinal BFRT groups (workData order after pdqsort/maxheapsort)",
    ),
    (
        "iter",
        "alpha_col",
        "HEkkDual::updateFtran -> HSimplexNla::ftran -> HFactor::ftranCall (incl. update etas)",
    ),
    (
        "iter",
        "theta_primal",
        "HEkkDual::updatePrimal (delta_primal / alpha_col)",
    ),
    (
        "iter",
        "numericalTrouble",
        "HEkkDual::updateVerify (alpha_col vs alpha_row)",
    ),
    (
        "iter",
        "edge_weight_mode",
        "HEkkDual::possiblyUseLiDualSteepestEdge / DSE-to-Devex switch",
    ),
    (
        "iter",
        "update_count",
        "HEkk::updatePivots / HFactor::update",
    ),
    (
        "iter",
        "updated_dual_objective_value",
        "HEkkDual::updateDual (objective change) / HEkkDualRow::updateFlip",
    ),
    (
        "iter",
        "total_synthetic_tick",
        "synthetic tick accounting: HVector::synthetic_tick in ftran/btran/PRICE, HEkk::total_synthetic_tick_",
    ),
    (
        "iter",
        "row_ep_density",
        "HEkk::updateOperationResultDensity (row_ep after HEkk::unitBtran)",
    ),
    (
        "iter",
        "row_ap_density",
        "HEkk::updateOperationResultDensity (row_ap after HEkk::tableauRowPrice)",
    ),
    (
        "iter",
        "col_aq_density",
        "HEkk::updateOperationResultDensity (col_aq after HEkk::pivotColumnFtran)",
    ),
    (
        "iter",
        "row_DSE_density",
        "HEkk::updateOperationResultDensity (col_DSE after HEkkDual::updateFtranDSE)",
    ),
    (
        "iter",
        "edge_weight_row_out_after",
        "HEkk::updateDualSteepestEdgeWeights / HEkkDual::updateFtranDSE",
    ),
    (
        "iter",
        "h_dual_edge_weight",
        "HEkk::updateDualSteepestEdgeWeights / HEkkDual::updateFtranDSE / HFactor::ftranCall on row_ep",
    ),
    (
        "iter",
        "h_workDual",
        "HEkkDual::updateDual / HEkkDualRow::updateDual / HEkkDual::shiftCost",
    ),
    (
        "iter",
        "h_baseValue",
        "HEkkDual::updatePrimal / HEkkDualRHS::updatePrimal / HEkkDual::updateFtranBFRT",
    ),
    ("iter", "h_basicIndex", "HEkk::updatePivots"),
    (
        "iter",
        "h_workValue",
        "HEkkPrimal::update / HEkk::flipBound",
    ),
    // primal iteration
    (
        "iter",
        "p.variable_in",
        "HEkkPrimal::chooseColumn / hyperChooseColumn* (CHUZC with edge weights)",
    ),
    (
        "iter",
        "p.move_in",
        "HEkkPrimal::chooseColumn (sign of workDual)",
    ),
    (
        "iter",
        "p.row_out",
        "HEkkPrimal::chooseRow (two-pass Harris ratio test) / assessPivot",
    ),
    ("iter", "p.variable_out", "HEkkPrimal::chooseRow"),
    (
        "iter",
        "p.alpha_col",
        "HEkk::pivotColumnFtran -> HFactor::ftranCall",
    ),
    (
        "iter",
        "p.alpha_row",
        "HEkk::unitBtran + HEkk::tableauRowPrice (HEkkPrimal::update)",
    ),
    (
        "iter",
        "p.theta",
        "HEkkPrimal::chooseRow / HEkkPrimal::update",
    ),
    ("iter", "p.value_in", "HEkkPrimal::update"),
    (
        "iter",
        "p.numericalTrouble",
        "HEkkPrimal::assessPivot / update",
    ),
    (
        "iter",
        "p.h_edge_weight",
        "HEkkPrimal::updateDevex / updatePrimalSteepestEdgeWeights",
    ),
    (
        "iter",
        "p.edge_weight_in",
        "HEkkPrimal::updateDevex / updatePrimalSteepestEdgeWeights",
    ),
    (
        "iter",
        "updated_primal_objective_value",
        "HEkkPrimal::update (objective change)",
    ),
    ("iter", "p", "HEkkPrimal::iterate"),
    (
        "iter_none",
        "",
        "HEkkDual::iterate / HEkkPrimal::iterate (no-pivot exit: optimality or rebuild trigger)",
    ),
    // rebuild / cleanup / starts
    (
        "*",
        "hash.workCost",
        "HEkk::initialiseCost (cost perturbation, random values) / HEkkDual::shiftBack / cleanup",
    ),
    ("*", "hash.workShift", "HEkkDual::shiftCost / shiftBack"),
    (
        "*",
        "hash.workLower",
        "HEkk::initialiseBound (phase-1 boxes / bound perturbation)",
    ),
    (
        "*",
        "hash.workUpper",
        "HEkk::initialiseBound (phase-1 boxes / bound perturbation)",
    ),
    (
        "*",
        "hash.workValue",
        "HEkk::initialiseNonbasicValueAndMove / HEkk::flipBound",
    ),
    (
        "*",
        "hash.nonbasicMove",
        "HEkk::initialiseNonbasicValueAndMove / flipBound / HEkkDual::correctDualInfeasibilities",
    ),
    (
        "*",
        "hash.workDual",
        "HEkk::computeDual (btran of basic costs + PRICE)",
    ),
    (
        "*",
        "hash.baseValue",
        "HEkk::computePrimal (ftran of -N x_N)",
    ),
    (
        "*",
        "hash.baseLower",
        "HEkk::computePrimal / updatePivots (base bounds)",
    ),
    (
        "*",
        "hash.baseUpper",
        "HEkk::computePrimal / updatePivots (base bounds)",
    ),
    (
        "*",
        "hash.dual_edge_weight",
        "HEkk::computeDualSteepestEdgeWeights / HEkkDual::rebuild weight handling",
    ),
    (
        "*",
        "hash.basicIndex",
        "HEkk::computeFactor (HFactor::build permutes basicIndex) / updatePivots",
    ),
    (
        "rebuild",
        "reason",
        "rebuild trigger in HEkkDual::iterate / HEkkPrimal::iterate (update limit, synthetic clock, optimality)",
    ),
    (
        "rebuild",
        "dualInfeasCount",
        "HEkkDual::rebuild -> HEkk::computeSimplexDualInfeasible",
    ),
    ("rebuild", "rebuild_reason_after", "HEkkDual::rebuild"),
    (
        "*",
        "scalars.total_synthetic_tick",
        "synthetic tick accounting (HVector::synthetic_tick, HFactor::build_synthetic_tick)",
    ),
    (
        "*",
        "scalars.build_synthetic_tick",
        "HEkk::computeFactor (HFactor::build_synthetic_tick)",
    ),
    (
        "*",
        "scalars.cost_perturbation",
        "HEkk::initialiseCost (perturbation base / max abs cost)",
    ),
    (
        "*",
        "scalars.updated_dual_objective_value",
        "HEkkDual::updateDual / HEkk::computeDualObjectiveValue",
    ),
    (
        "*",
        "scalars.dual_objective_value",
        "HEkk::computeDualObjectiveValue",
    ),
    (
        "*",
        "scalars.primal_objective_value",
        "HEkk::computePrimalObjectiveValue",
    ),
    (
        "*",
        "scalars.num_primal_infeasibilities",
        "HEkk::computeSimplexPrimalInfeasible",
    ),
    (
        "*",
        "scalars.sum_primal_infeasibilities",
        "HEkk::computeSimplexPrimalInfeasible",
    ),
    (
        "*",
        "scalars.num_dual_infeasibilities",
        "HEkk::computeSimplexDualInfeasible",
    ),
    (
        "*",
        "scalars.sum_dual_infeasibilities",
        "HEkk::computeSimplexDualInfeasible",
    ),
    (
        "*",
        "scalars.update_count",
        "HEkk::updatePivots / computeFactor resets",
    ),
    (
        "*",
        "scalars.iteration_count",
        "HEkkDual::iterate / HEkkPrimal::iterate",
    ),
    ("*", "scalars.", "HEkk scalar bookkeeping (see field name)"),
    (
        "*",
        "state.",
        "HEkkPrimal state (solve_phase / rebuild_reason / edge weights)",
    ),
    (
        "ekk_solve_initialised",
        "numTotRandomValue",
        "HEkk::initialiseSimplexLpRandomVectors (the SECOND draw on HEkk::random_)",
    ),
    (
        "ekk_solve_initialised",
        "numTotPermutation",
        "HEkk::initialiseSimplexLpRandomVectors (HighsRandom::shuffle)",
    ),
    (
        "ekk_solve_initialised",
        "numColPermutation",
        "HEkk::initialiseSimplexLpRandomVectors (HighsRandom::shuffle)",
    ),
    (
        "ekk_solve_initialised",
        "ar_matrix",
        "HEkk::initialisePartitionedRowwiseMatrix / HighsSparseMatrix::createRowwisePartitioned",
    ),
    (
        "ekk_solve_initialised",
        "basis",
        "HEkk::initialiseForSolve / setBasis (logical basis)",
    ),
    ("*", "arrays.workCost", "HEkk::initialiseCost"),
    ("*", "arrays.workLower", "HEkk::initialiseBound"),
    ("*", "arrays.workUpper", "HEkk::initialiseBound"),
    ("*", "arrays.workRange", "HEkk::initialiseBound"),
    (
        "*",
        "arrays.workValue",
        "HEkk::initialiseNonbasicValueAndMove",
    ),
    ("*", "arrays.workDual", "HEkk::computeDual"),
    ("*", "arrays.base", "HEkk::computePrimal"),
    (
        "*",
        "arrays.dual_edge_weight",
        "HEkkDual::initialiseSolve / HEkk::computeDualSteepestEdgeWeights",
    ),
    (
        "*",
        "status",
        "HighsSimplexStatus bookkeeping (HEkk::initialiseForSolve / computeFactor / rebuild)",
    ),
    (
        "invert",
        "",
        "HEkk::computeFactor -> HSimplexNla::invert -> HFactor::build",
    ),
    (
        "phase_end",
        "",
        "HEkkDual::solvePhase1/2 / HEkkPrimal::solvePhase1/2 exit logic",
    ),
    (
        "cleanup",
        "",
        "HEkkDual::cleanup / HEkkPrimal::cleanup (remove perturbations, recompute)",
    ),
    (
        "ekk_solve_exit",
        "arrays.workDual",
        "HEkk::returnFromSolve (zeroBasicDuals)",
    ),
    ("ekk_solve_exit", "", "HEkk::returnFromSolve"),
    (
        "result",
        "",
        "HEkk unscaling / Highs::run solution and basis assembly",
    ),
    (
        "ekk_solve_enter",
        "",
        "LP intake: scaled LP, options, HighsSimplexStatus before HEkk::initialiseForSolve",
    ),
];

fn hint_for_field(t: &str, rel: &str) -> Option<&'static str> {
    HINTS
        .iter()
        .find(|(ty, prefix, _)| (*ty == "*" || *ty == t) && rel.starts_with(prefix))
        .map(|(_, _, h)| *h)
}

fn hints_for(
    kind: &DivergenceKind,
    et: Option<&str>,
    at: Option<&str>,
    mismatches: &[Mismatch],
    ctx_prefix: Option<String>,
) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    let mut add = |h: &str| {
        if !out.iter().any(|x| x == h) {
            out.push(h.to_string());
        }
    };
    match kind {
        DivergenceKind::Type => {
            let pair = (et.unwrap_or(""), at.unwrap_or(""));
            match pair {
                ("rebuild", "iter")
                | ("iter", "rebuild")
                | ("rebuild", "iter_none")
                | ("iter_none", "rebuild") => {
                    add(
                        "rebuild trigger differs: compare total_synthetic_tick vs build_synthetic_tick (synthetic clock), update_count vs update_limit, and the rebuild_reason set in HEkkDual::iterate / HEkkPrimal::iterate",
                    );
                }
                ("invert", _) | (_, "invert") => add(
                    "HEkk::computeFactor called at a different point (rebuild / reinversion trigger)",
                ),
                ("phase_end", _) | (_, "phase_end") | ("cleanup", _) | (_, "cleanup") => {
                    add(
                        "phase / termination decision differs: HEkkDual::solvePhase1/2, HEkkPrimal::solvePhase1/2, and the optimality / infeasibility tests after rebuild",
                    );
                }
                _ => add(
                    "control flow diverged; check the preceding record's rebuild_reason / solve_phase",
                ),
            }
        }
        DivergenceKind::ActualEnded => {
            add("the port stopped early: termination test or missing hook emission")
        }
        DivergenceKind::ExpectedEnded => {
            add("the port continued past HiGHS' end: termination test or extra hook emission")
        }
        DivergenceKind::Fields => {
            let t = et.unwrap_or("");
            let prefix = ctx_prefix.unwrap_or_default();
            for m in mismatches {
                // m.path = "<file>:<line> [<type>].<rel>"
                let rel = m
                    .path
                    .strip_prefix(&prefix)
                    .and_then(|s| s.split_once("].").map(|(_, r)| r))
                    .unwrap_or(&m.path);
                if let Some(h) = hint_for_field(t, rel) {
                    add(h);
                }
            }
        }
    }
    out
}
