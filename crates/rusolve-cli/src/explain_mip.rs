//! `rusolve explain <model.mps> --mip-iis [--per-test-node-limit N] [--time-limit S]`:
//! an irreducible infeasible subsystem (IIS) of a MIP itself, labelled
//! `"kind": "mip"` (MIP part).
//!
//! ## Definition used
//!
//! The elements are the finite bounds of rows and of columns, one element per
//! side (a ranged or equality row has two). Integrality of every column is
//! always kept: a MIP IIS here is a set of bounds that is infeasible as a MIP
//! (all columns still integer where the model says so, all other bounds and rows
//! dropped), and from which no single bound can be removed without making the
//! MIP feasible. Same granularity as the LP IIS of `explain.rs`.
//!
//! ## Sources
//!
//! * Deletion filter: J. W. Chinneck, *Feasibility and Infeasibility in
//!   Optimization: Algorithms and Computational Methods*, Springer 2008, the
//!   deletion filter for finding an IIS (originating with Dravnieks, and Chinneck
//!   and Dravnieks, "Locating minimal infeasible constraint sets in linear
//!   programs", ORSA J. Computing 3(2), 1991). Each element is tentatively removed;
//!   if the rest is still infeasible it stays out, otherwise it is a member.
//!   The chapter/section numbers of the book were not available offline; none is cited.
//! * Applying such filters to integer programs, with the feasibility test being
//!   a MIP solve: O. Guieu and J. W. Chinneck, "Analyzing infeasible
//!   mixed-integer and integer linear programs", INFORMS Journal on Computing
//!   11(1), 63-77, 1999 (title, venue and pages checked; the full text was not
//!   read, so we claim only the idea of filters with MIP feasibility tests, not
//!   their specific algorithms).
//! * Grouping (block removal) to cut the number of MIP solves: the standard
//!   speed-up of deletion filters (testing a block of elements at once, halving
//!   the block when the rest becomes feasible). Our block schedule is our own;
//!   the soundness argument is below.
//! * Seed: when the LP relaxation is infeasible, its IIS (ported HiGHS
//!   `getIis`, as in `explain.rs`) is infeasible as a MIP too (the relaxation has
//!   more solutions), so the filter starts from it instead of from all bounds.
//!   `rusolve feasrelax` is not used: its elastic solution is for the
//!   relaxation and gives no guarantee that the set it points at is infeasible
//!   as a MIP, which a seed needs.
//!
//! ## Soundness
//!
//! A test is a MIP solve of the zero-cost model restricted to a set of
//! elements. It answers Infeasible (HiGHS proved it), Feasible (a solution was
//! found AND it passes the solution checker `check::check_solution` against the
//! restricted model: bounds, rows and integrality) or Undetermined (a node or
//! time limit, any other status, or a witness that fails the check). Dropping
//! elements only relaxes a model, so infeasibility is monotone: an element
//! removed after an Infeasible test stays removable from every subset of the
//! tested set, and an element whose removal gave Feasible is necessary in every
//! subset of the tested set that contains it. Hence the final set is infeasible
//! and each Feasible-decided member is necessary. An element whose test was
//! Undetermined is kept (the set stays infeasible) but reported, and
//! `irreducible` is then `false`.
//!
//! Verification (unless `--no-verify`) re-tests the final set and every member
//! from scratch, independent of the search.

use crate::check::{Tolerances, check_solution};
use crate::explain::{Args, Member, col_name, members_json, row_name, side_name};
use rusolve_ekk::lp_data::highs_iis::HighsIis;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::{Map, Value, json};
use std::process::ExitCode;
use std::time::Instant;

/// `iis_strategy` for the relaxation seed: kIisStrategyFromLp + kIisStrategyIrreducible.
const IIS_STRATEGY: i32 = 6;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Elem {
    is_row: bool,
    index: usize,
    lower: bool,
}

enum Outcome {
    Infeasible,
    Feasible,
    Undetermined(String),
}

struct Oracle<'a> {
    lp: &'a HighsLp,
    node_limit: i32,
    time_limit: Option<f64>,
    start: Instant,
    solves: usize,
}

impl Oracle<'_> {
    /// The zero-cost MIP restricted to `set`: other rows and column bounds free,
    /// integrality kept.
    fn build(&self, set: &[Elem]) -> HighsLp {
        let lp = self.lp;
        let mut r = lp.clone();
        r.col_cost_.iter_mut().for_each(|x| *x = 0.0);
        r.offset_ = 0.0;
        r.row_lower_.iter_mut().for_each(|b| *b = -kHighsInf);
        r.row_upper_.iter_mut().for_each(|b| *b = kHighsInf);
        r.col_lower_.iter_mut().for_each(|b| *b = -kHighsInf);
        r.col_upper_.iter_mut().for_each(|b| *b = kHighsInf);
        for e in set {
            match (e.is_row, e.lower) {
                (true, true) => r.row_lower_[e.index] = lp.row_lower_[e.index],
                (true, false) => r.row_upper_[e.index] = lp.row_upper_[e.index],
                (false, true) => r.col_lower_[e.index] = lp.col_lower_[e.index],
                (false, false) => r.col_upper_[e.index] = lp.col_upper_[e.index],
            }
        }
        r
    }

    fn test(&mut self, set: &[Elem]) -> Outcome {
        let r = self.build(set);
        // Inconsistent bounds are decided by inspection (as explain.rs does).
        if (0..r.num_col_ as usize).any(|j| r.col_lower_[j] > r.col_upper_[j])
            || (0..r.num_row_ as usize).any(|i| r.row_lower_[i] > r.row_upper_[i])
        {
            return Outcome::Infeasible;
        }
        let remaining = self
            .time_limit
            .map(|t| t - self.start.elapsed().as_secs_f64());
        if remaining.is_some_and(|t| t <= 0.0) {
            return Outcome::Undetermined("time limit reached before the test".into());
        }
        self.solves += 1;
        let mut h = Highs::with_options(crate::base_options());
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_int("mip_max_nodes", self.node_limit);
        if let Some(t) = remaining {
            h.set_option_value_double("time_limit", t);
        }
        if h.pass_model(r.clone()) == HighsStatus::kError || h.run() == HighsStatus::kError {
            return Outcome::Undetermined("the solver reported an error".into());
        }
        let status = h.get_model_status();
        if status == HighsModelStatus::kInfeasible {
            return Outcome::Infeasible;
        }
        let sol = h.get_solution();
        let has = sol.value_valid
            && sol.col_value.len() == r.num_col_ as usize
            && (status == HighsModelStatus::kOptimal
                || h.get_info().primal_solution_status == kSolutionStatusFeasible);
        if !has {
            return Outcome::Undetermined(format!(
                "no verdict within the limits (model status {})",
                util_model_status_to_string(status)
            ));
        }
        let tol = Tolerances {
            primal_feasibility_tolerance: 1e-7,
            mip_feasibility_tolerance: 1e-6,
            applied: 1e-6,
        };
        match check_solution(&r, &sol.col_value, tol, 1, "", "") {
            Ok(rep) if rep.verdict == "feasible" => Outcome::Feasible,
            _ => Outcome::Undetermined(
                "a solution was returned but it fails the solution checker".into(),
            ),
        }
    }
}

fn all_elems(lp: &HighsLp) -> Vec<Elem> {
    let mut v = Vec::new();
    let mut push = |is_row, index, lower| {
        v.push(Elem {
            is_row,
            index,
            lower,
        })
    };
    for i in 0..lp.num_row_ as usize {
        if lp.row_lower_[i].is_finite() {
            push(true, i, true);
        }
        if lp.row_upper_[i].is_finite() {
            push(true, i, false);
        }
    }
    for j in 0..lp.num_col_ as usize {
        if lp.col_lower_[j].is_finite() {
            push(false, j, true);
        }
        if lp.col_upper_[j].is_finite() {
            push(false, j, false);
        }
    }
    v
}

/// The ported LP IIS of the zero-cost relaxation as elements; `None` if the
/// relaxation is feasible or HiGHS returns no usable set.
fn relaxation_seed(lp: &HighsLp) -> Option<Vec<Elem>> {
    let mut r = lp.clone();
    r.integrality_.clear();
    r.col_cost_.iter_mut().for_each(|x| *x = 0.0);
    r.offset_ = 0.0;
    let mut h = Highs::with_options(crate::base_options());
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_int("iis_strategy", IIS_STRATEGY);
    if h.pass_model(r) == HighsStatus::kError || h.run() == HighsStatus::kError {
        return None;
    }
    if h.get_model_status() != HighsModelStatus::kInfeasible {
        return None;
    }
    let mut iis = HighsIis::default();
    if h.get_iis(&mut iis) == HighsStatus::kError || !iis.valid_ {
        return None;
    }
    let mut v = Vec::new();
    for (k, &i) in iis.row_index_.iter().enumerate() {
        let (lo, up) = crate::explain::decode(iis.row_bound_[k]);
        if lo {
            v.push(Elem {
                is_row: true,
                index: i as usize,
                lower: true,
            });
        }
        if up {
            v.push(Elem {
                is_row: true,
                index: i as usize,
                lower: false,
            });
        }
    }
    for (k, &j) in iis.col_index_.iter().enumerate() {
        let (lo, up) = crate::explain::decode(iis.col_bound_[k]);
        if lo {
            v.push(Elem {
                is_row: false,
                index: j as usize,
                lower: true,
            });
        }
        if up {
            v.push(Elem {
                is_row: false,
                index: j as usize,
                lower: false,
            });
        }
    }
    (!v.is_empty()).then_some(v)
}

/// Grouped deletion filter on an infeasible `seed`. Returns `(set, undetermined)`
/// where `set` is infeasible and contains `undetermined`.
fn deletion_filter(o: &mut Oracle, seed: Vec<Elem>) -> (Vec<Elem>, Vec<Elem>) {
    let mut fixed: Vec<Elem> = Vec::new(); // removal gave Feasible: necessary
    let mut undet: Vec<Elem> = Vec::new(); // test undetermined: kept, reported
    let mut cand = seed;
    let mut g = cand.len().div_ceil(2).max(1);
    while !cand.is_empty() {
        let mut kept: Vec<Elem> = Vec::new();
        let mut rest = cand; // in order; `kept` holds what was passed over in this pass
        let mut pos = 0usize;
        while pos < rest.len() {
            let end = (pos + g).min(rest.len());
            // set without the block = fixed + undet + kept + rest after the block
            let mut set: Vec<Elem> = Vec::with_capacity(fixed.len() + undet.len() + rest.len());
            set.extend_from_slice(&fixed);
            set.extend_from_slice(&undet);
            set.extend_from_slice(&kept);
            set.extend_from_slice(&rest[end..]);
            match o.test(&set) {
                Outcome::Infeasible => {
                    // the block is not needed: drop it for good
                    rest.drain(pos..end);
                }
                Outcome::Feasible if end - pos == 1 => {
                    fixed.push(rest[pos]);
                    pos = end;
                }
                Outcome::Undetermined(_) if end - pos == 1 => {
                    undet.push(rest[pos]);
                    pos = end;
                }
                _ => {
                    // a block that holds a member (or no verdict): look at it
                    // again with smaller blocks
                    kept.extend_from_slice(&rest[pos..end]);
                    pos = end;
                }
            }
        }
        cand = kept;
        if g == 1 {
            break;
        }
        g = g.div_ceil(2);
    }
    let mut set = fixed;
    set.extend_from_slice(&undet);
    (set, undet)
}

fn to_members(set: &[Elem], rows: bool) -> Vec<Member> {
    let mut idx: Vec<usize> = set
        .iter()
        .filter(|e| e.is_row == rows)
        .map(|e| e.index)
        .collect();
    idx.sort_unstable();
    idx.dedup();
    idx.into_iter()
        .map(|i| Member {
            index: i,
            lower: set
                .iter()
                .any(|e| e.is_row == rows && e.index == i && e.lower),
            upper: set
                .iter()
                .any(|e| e.is_row == rows && e.index == i && !e.lower),
            free_col: false,
        })
        .collect()
}

fn elem_json(lp: &HighsLp, e: &Elem) -> Value {
    json!({
        "type": if e.is_row { "row" } else { "column" },
        "name": if e.is_row { row_name(lp, e.index) } else { col_name(lp, e.index) },
        "bound": if e.lower { "lower" } else { "upper" },
    })
}

fn tri(v: Option<bool>) -> &'static str {
    match v {
        Some(true) => "yes",
        Some(false) => "NO",
        None => "unknown",
    }
}

fn finish(a: &Args, out: Map<String, Value>, text: &str) -> ExitCode {
    if a.json {
        println!("{}", Value::Object(out));
    } else {
        print!("{text}");
    }
    ExitCode::SUCCESS
}

pub(crate) fn run_mip(a: &Args, lp: HighsLp, num_integer: usize) -> ExitCode {
    let mut o = Oracle {
        lp: &lp,
        node_limit: a.per_test_node_limit,
        time_limit: a.time_limit,
        start: Instant::now(),
        solves: 0,
    };
    let all = all_elems(&lp);

    // Establish an infeasible starting set (the relaxation's IIS if it has one).
    let mut seed_kind = "all_bounds";
    let mut start: Option<Vec<Elem>> = None;
    let mut full_outcome = None;
    if let Some(s) = relaxation_seed(&lp)
        && matches!(o.test(&s), Outcome::Infeasible)
    {
        seed_kind = "relaxation_iis";
        start = Some(s);
    }
    if start.is_none() {
        match o.test(&all) {
            Outcome::Infeasible => start = Some(all.clone()),
            other => full_outcome = Some(other),
        }
    }
    let seed_elements = start.as_ref().map_or(all.len(), Vec::len);

    let mut out = Map::new();
    out.insert(
        "schema_version".into(),
        json!(crate::explain::SCHEMA_VERSION),
    );
    out.insert(
        "model".into(),
        json!({"path": a.path, "num_rows": lp.num_row_, "num_cols": lp.num_col_,
               "num_integer": num_integer}),
    );
    out.insert("kind".into(), json!("mip"));
    let mut text = String::new();

    let Some(start) = start else {
        let (status, note) = match full_outcome {
            Some(Outcome::Feasible) => (
                "Feasible",
                "the MIP is feasible (a solution passing the solution check was found): nothing to explain"
                    .to_string(),
            ),
            Some(Outcome::Undetermined(why)) => (
                "Unknown",
                format!(
                    "the MIP could not be decided within the limits ({why}): no infeasible set is claimed"
                ),
            ),
            _ => ("Unknown", "internal: no outcome".to_string()),
        };
        out.insert("status".into(), json!(status));
        out.insert("note".into(), json!(note));
        out.insert("explanation".into(), Value::Null);
        out.insert("verified".into(), Value::Null);
        text.push_str(&format!("status (mip): {status}\n{note}\n"));
        return finish(a, out, &text);
    };

    let (set, undet) = deletion_filter(&mut o, start);
    let search_solves = o.solves;
    let rows = to_members(&set, true);
    let cols = to_members(&set, false);
    let mut truncated = false;
    let rj = members_json(&lp, &rows, true, a.max_members, &mut truncated);
    let cj = members_json(&lp, &cols, false, a.max_members, &mut truncated);

    text.push_str("status (mip): Infeasible\n");
    text.push_str(&format!(
        "IIS (mip): {} rows, {} columns, irreducible: {}\n",
        rows.len(),
        cols.len(),
        if undet.is_empty() {
            "yes"
        } else {
            "NOT proven"
        }
    ));
    for m in rows.iter().take(a.max_members) {
        text.push_str(&format!(
            "  row {}: {} bound ({} <= . <= {})\n",
            row_name(&lp, m.index),
            side_name(m.lower, m.upper),
            lp.row_lower_[m.index],
            lp.row_upper_[m.index]
        ));
    }
    for m in cols.iter().take(a.max_members) {
        text.push_str(&format!(
            "  column {}: {} bound ({} <= . <= {})\n",
            col_name(&lp, m.index),
            side_name(m.lower, m.upper),
            lp.col_lower_[m.index],
            lp.col_upper_[m.index]
        ));
    }
    for u in &undet {
        let d = elem_json(&lp, u);
        text.push_str(&format!(
            "  undetermined: {} bound of {} {} (test hit a limit; kept in the set)\n",
            d["bound"].as_str().unwrap_or(""),
            d["type"].as_str().unwrap_or(""),
            d["name"].as_str().unwrap_or("")
        ));
    }

    // Verification, independent of the search.
    let mut v_json = json!({
        "performed": false, "set_infeasible": null, "irreducible": null,
        "checks": [], "failures": [], "verified": false,
    });
    let mut v_ok = false;
    if a.verify {
        let mut failures: Vec<String> = Vec::new();
        let mut checks: Vec<Value> = Vec::new();
        let set_inf = match o.test(&set) {
            Outcome::Infeasible => Some(true),
            Outcome::Feasible => {
                failures.push("the set is feasible when re-solved".into());
                Some(false)
            }
            Outcome::Undetermined(w) => {
                failures.push(format!("could not decide feasibility of the set: {w}"));
                None
            }
        };
        let mut irreducible = None;
        if set_inf == Some(true) {
            let mut irr = true;
            for (k, e) in set.iter().enumerate() {
                let mut rest = set.clone();
                rest.remove(k);
                let (feasible, why) = match o.test(&rest) {
                    Outcome::Feasible => (true, None),
                    Outcome::Infeasible => (false, Some("still infeasible without it".to_string())),
                    Outcome::Undetermined(w) => (false, Some(w)),
                };
                if !feasible {
                    irr = false;
                    let d = elem_json(&lp, e);
                    failures.push(format!(
                        "dropping the {} bound of {} {}: {}",
                        d["bound"].as_str().unwrap_or(""),
                        d["type"].as_str().unwrap_or(""),
                        d["name"].as_str().unwrap_or(""),
                        why.unwrap_or_default()
                    ));
                }
                checks.push(json!({
                    "dropped": elem_json(&lp, e),
                    "rest_feasible": feasible,
                    "witness_checked": feasible,
                }));
            }
            irreducible = Some(irr);
        }
        v_ok = failures.is_empty();
        text.push_str(&format!(
            "verification: set infeasible: {}, irreducible: {}\n",
            tri(set_inf),
            tri(irreducible)
        ));
        for f in &failures {
            text.push_str(&format!("  FAILED: {f}\n"));
        }
        v_json = json!({
            "performed": true, "set_infeasible": set_inf, "irreducible": irreducible,
            "checks": checks, "failures": failures, "verified": v_ok,
        });
    } else {
        text.push_str("verification skipped (--no-verify)\n");
    }
    let verification_solves = o.solves - search_solves;
    let irreducible = undet.is_empty();
    let ok = v_ok && irreducible;
    text.push_str(&format!(
        "MIP solves: {search_solves} in the search, {verification_solves} in verification (seed: {seed_kind}, {seed_elements} bounds)\n"
    ));
    text.push_str(&format!("verified: {ok}\n"));

    let undet_json: Vec<Value> = undet.iter().map(|x| elem_json(&lp, x)).collect();
    let mut e = Map::new();
    e.insert("type".into(), json!("infeasible"));
    e.insert(
        "iis".into(),
        json!({
            "kind": "mip", "num_rows": rows.len(), "num_cols": cols.len(),
            "rows": rj, "columns": cj, "truncated": truncated,
            "irreducible": irreducible,
            "undetermined": undet_json,
        }),
    );
    e.insert(
        "mip_search".into(),
        json!({
            "algorithm": "grouped deletion filter with MIP feasibility tests",
            "sources": [
                "Chinneck, Feasibility and Infeasibility in Optimization, Springer 2008 (deletion filter)",
                "Guieu and Chinneck, INFORMS J. Computing 11(1), 63-77, 1999 (filters on MIPs)"
            ],
            "seed": seed_kind,
            "seed_elements": seed_elements,
            "num_solves": search_solves,
            "num_verification_solves": verification_solves,
            "per_test_node_limit": a.per_test_node_limit,
            "time_limit": a.time_limit,
            "undetermined": undet_json,
        }),
    );
    e.insert("verification".into(), v_json);
    e.insert("verified".into(), json!(ok));
    out.insert("status".into(), json!("Infeasible"));
    out.insert(
        "note".into(),
        json!("irreducible infeasible subsystem of the MIP itself (integrality of every column kept); the bounds listed are the members"),
    );
    out.insert("explanation".into(), Value::Object(e));
    out.insert("verified".into(), json!(ok));
    finish(a, out, &text)
}
