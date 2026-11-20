//! `rusolve diagnose <model.mps>`
//! says, for an infeasible model, whether the conflict is most likely a
//! *formulation* error (the generating code builds constraint families that
//! contradict each other) or a *data* error (a few numbers out of line), with
//! the evidence and the concrete thing to inspect.
//!
//! ```text
//! rusolve diagnose <model.mps> [--format text|json] [--time-limit S]
//! ```
//!
//! ## Why a new subcommand and not `explain --cause`
//!
//! `explain` reports a checked certificate (IIS, Farkas ray) and its schema is
//! about proof. `diagnose` makes a *heuristic judgement* from several inputs
//! (IIS, Farkas ray, the minimum-count relaxation, family statistics, extra
//! solves) and carries a verdict, a confidence and ranked evidence: a different
//! contract (it can be wrong, and says how sure it is), so it gets its own
//! schema `rusolve.diagnose/1` rather than loosening the proof-shaped one. The
//! building blocks are reused, not copied: `explain::{solve_for_explanation,
//! iis_members, ray_multipliers, is_infeasible}` and
//! `feasrelax::min_count_relaxation`.
//!
//! ## Families (ours)
//!
//! Modelling layers name rows and columns `stem` + index groups. The stem is
//! found by [`family_stem`], which repeatedly strips from the end:
//!
//! * a bracketed group `[..]`, `(..)`, `{..}`, `<..>` (`cap[3,5]`, `cap(3)`);
//! * a run of digits and the separator before it (`cap_3_5`, `flow#12`,
//!   `flow-12`, `flow.12`, `flow12`, `R09`), one `_ . # - : , @ /` or space.
//!
//! Names that carry symbolic indices (`supply_Rotterdam`) have no digits to
//! strip, so a second pass ([`assign_families`]) drops trailing separator
//! tokens until the prefix is shared by at least two names. A name that shares
//! nothing is a family of one. Rows and columns form separate family sets.
//! The families are only informative when some family has at least three
//! members and a stem of two or more characters (`naming_informative`);
//! otherwise (Netlib's `R09`/`X05`) confidence is reduced and a note says so.
//!
//! ## Signals
//!
//! Every signal below is cited; where no source exists the rule is ours and the
//! rationale is in the comment on it. The weights are ours (heuristic scores,
//! not calibrated probabilities).
//!
//! * IIS and Farkas ray: HiGHS `Highs::getIis`/`getDualRay` via `explain`.
//!   Chinneck, *Feasibility and Infeasibility in Optimization* (Springer 2008):
//!   an IIS, and a minimum-cardinality set of constraints whose removal gives
//!   feasibility, localise an infeasibility so the analyst inspects the few
//!   implicated constraints. That book does not split causes into
//!   "modelling vs data"; that split, and everything below, is ours.
//! * (b) data: members of the IIS or of the minimum-count relaxation whose
//!   right-hand side or bound is an outlier within its family by the *modified
//!   z-score* `0.6745 (x - median) / MAD`, flagged at `|z| > 3.5` (Iglewicz &
//!   Hoaglin, *How to Detect and Handle Outliers*, ASQC 1993). When MAD is 0
//!   they use `(x - median) / (1.253314 * meanAD)`; when that is 0 too the
//!   other members are all equal and any different value is an outlier. Ours:
//!   a family needs at least [`MIN_FAMILY`] members with that bound kind, its
//!   median must be nonzero, and the value must differ from the median by at
//!   least the median itself ([`MIN_REL_DEV`]: a z of 4 on a spread of 1e-9 is
//!   rounding, milder spreads are ordinary variation, and we look for gross
//!   errors such as a 100x typo). The decisive test is a
//!   counterfactual (ours): replace the outliers by the family median and
//!   re-solve; if the model becomes feasible, fixing just those cells restores
//!   feasibility.
//! * (a) formulation: the IIS covers at least half (and at least three) of the
//!   members of two or more families (a family-wide contradiction), ideally
//!   with one coefficient pattern shared by the IIS members of a family (ours:
//!   a contradiction produced by generator code repeats). The decisive test is
//!   again a counterfactual (ours): replace every right-hand side/bound of the
//!   families that touch the IIS or the relaxation by its family median
//!   ("typical data" everywhere); if the model is *still* infeasible no few
//!   cells are to blame, the constraint pattern itself conflicts (or the whole
//!   data is scaled wrongly, a unit error in the generator, which is also a
//!   code fix).
//! * minimum-count relaxation (`feasrelax --mode min-count`, Chinneck 2008
//!   ch. 6 elastic programming): relaxed cells that are outliers support data;
//!   relaxed cells that are ordinary family members support formulation (the
//!   solver had to give up a typical constraint).
//! * (c) Farkas concentration (ours): the effective number of rows carrying the
//!   ray, `N = (sum|y|)^2 / sum y^2` (the inverse Simpson index; Hill number of
//!   order 2), is at most 3: a few rows decide, mild support for data. A family
//!   with at least three ray rows covering at least half of it with near-equal
//!   magnitudes (`(sum|y|)^2 / (k sum y^2) >= 0.8`): the family acts as a block,
//!   mild support for formulation. This signal is weak because a plain
//!   supply-versus-demand total mismatch also spreads the ray evenly.
//!
//! The verdict is `data` or `formulation` when the score margin is at least 2,
//! else `unclear` (confidence 0; the leaning stays visible in `scores`). `confidence` is `margin / (data + formulation)` capped at
//! 0.95 (scaled by 0.6 without informative names, and capped at 0.6 when no
//! counterfactual confirmed the verdict): a heuristic rank, not a probability.
//!
//! The model is solved as an LP (integrality dropped, as `explain` does).
//!
//! Exit codes: 0 ran, 2 usage, 3 input, 4 internal.

use crate::explain::{
    col_name, entries, iis_members, is_infeasible, ray_multipliers, row_name, solve_for_explanation,
};
use crate::feasrelax::min_count_relaxation;
use crate::mip_highs::read_mip;
use rusolve_ekk::lp_data::highs_iis::HighsIis;
use rusolve_ekk::prelude::{HighsLp, HighsModelStatus, HighsStatus, HighsVarType, kHighsInf};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::process::ExitCode;

pub const SCHEMA_VERSION: &str = "rusolve.diagnose/1";

const USAGE: &str = "usage: rusolve diagnose <model.mps> [--format text|json] [--time-limit S]";

/// Smallest family (members with the bound kind) the outlier statistic is applied to.
const MIN_FAMILY: usize = 5;
/// Iglewicz & Hoaglin (1993): modified z-scores above 3.5 are potential outliers.
const Z_OUTLIER: f64 = 3.5;
/// Ours: minimum relative distance from the median for an outlier (at least a
/// doubling, or a sign flip: the errors we look for are gross, a 100x or 1000x
/// typo, a unit slip; milder spreads are ordinary data variation).
const MIN_REL_DEV: f64 = 1.0;
/// Ours: at most this many outlier cells count as "a few".
const FEW_CELLS: usize = 5;
/// Ours: score margin needed for a verdict other than `unclear`.
const DECISIVE_MARGIN: f64 = 2.0;
/// Ours: families touched by the conflict beyond which it counts as diffuse.
const MAX_LOCAL_FAMILIES: usize = 4;

// ---------------------------------------------------------------- args

struct Args {
    path: String,
    json: bool,
    time_limit: f64,
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut path = None;
    let mut json = false;
    let mut time_limit = 30.0;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--format" => match it.next().map(String::as_str) {
                Some("json") => json = true,
                Some("text") => json = false,
                other => return Err(format!("--format wants text or json, got {other:?}")),
            },
            "--time-limit" => {
                let v = it.next().ok_or("--time-limit needs a value")?;
                time_limit = v
                    .parse()
                    .map_err(|_| format!("--time-limit: not a number: {v}"))?;
            }
            s if s.starts_with("--") => return Err(format!("unknown option {s}")),
            s => {
                if path.replace(s.to_string()).is_some() {
                    return Err("more than one model path".into());
                }
            }
        }
    }
    Ok(Args {
        path: path.ok_or("missing model path")?,
        json,
        time_limit,
    })
}

// ---------------------------------------------------------------- families

fn is_sep(ch: char) -> bool {
    matches!(ch, '_' | '.' | '#' | '-' | ':' | ',' | '@' | '/' | ' ')
}

/// The name without its trailing index groups: see the module comment.
pub(crate) fn family_stem(name: &str) -> String {
    let mut s: Vec<char> = name.trim().chars().collect();
    while let Some(&last) = s.last() {
        let opener = match last {
            ']' => Some('['),
            ')' => Some('('),
            '}' => Some('{'),
            '>' => Some('<'),
            _ => None,
        };
        if let Some(open) = opener {
            if let Some(p) = s.iter().rposition(|&x| x == open) {
                s.truncate(p);
                continue;
            }
            break;
        }
        if last.is_ascii_digit() {
            while s.last().is_some_and(|c| c.is_ascii_digit()) {
                s.pop();
            }
            if s.last().is_some_and(|&c| is_sep(c)) {
                s.pop();
            }
            continue;
        }
        if is_sep(last) {
            s.pop();
            continue;
        }
        break;
    }
    s.into_iter().collect()
}

/// Stem candidates of one name, longest first: the digit/bracket stem, then
/// that stem with trailing separator-delimited tokens dropped one at a time.
fn stem_candidates(name: &str) -> Vec<String> {
    let mut out = vec![family_stem(name)];
    loop {
        let cur = out.last().unwrap().clone();
        match cur.rfind(is_sep) {
            Some(p) if p > 0 => out.push(cur[..p].trim_end_matches(is_sep).to_string()),
            _ => break,
        }
        if out.last().unwrap().is_empty() {
            out.pop();
            break;
        }
    }
    out
}

/// One stem per name: the longest candidate that at least two names share, else
/// the name's own digit/bracket stem.
pub(crate) fn assign_families(names: &[String]) -> Vec<String> {
    let cands: Vec<Vec<String>> = names.iter().map(|n| stem_candidates(n)).collect();
    let mut count: HashMap<&str, usize> = HashMap::new();
    for c in &cands {
        for s in c {
            *count.entry(s.as_str()).or_insert(0) += 1;
        }
    }
    cands
        .iter()
        .map(|c| {
            c.iter()
                .find(|s| count[s.as_str()] >= 2)
                .unwrap_or(&c[0])
                .clone()
        })
        .collect()
}

struct Family {
    is_row: bool,
    stem: String,
    members: Vec<usize>,
}

struct Families {
    list: Vec<Family>,
    /// family index of every row, then of every column
    of_row: Vec<usize>,
    of_col: Vec<usize>,
}

impl Families {
    fn build(lp: &HighsLp) -> Self {
        let mut list: Vec<Family> = Vec::new();
        let mut index: HashMap<(bool, String), usize> = HashMap::new();
        let mut place = |is_row: bool, stems: Vec<String>| -> Vec<usize> {
            stems
                .into_iter()
                .enumerate()
                .map(|(i, st)| {
                    let k = *index.entry((is_row, st.clone())).or_insert_with(|| {
                        list.push(Family {
                            is_row,
                            stem: st,
                            members: Vec::new(),
                        });
                        list.len() - 1
                    });
                    list[k].members.push(i);
                    k
                })
                .collect()
        };
        let rn: Vec<String> = (0..lp.num_row_ as usize).map(|i| row_name(lp, i)).collect();
        let cn: Vec<String> = (0..lp.num_col_ as usize).map(|j| col_name(lp, j)).collect();
        let of_row = place(true, assign_families(&rn));
        let of_col = place(false, assign_families(&cn));
        Families {
            list,
            of_row,
            of_col,
        }
    }

    fn of(&self, is_row: bool, i: usize) -> usize {
        if is_row {
            self.of_row[i]
        } else {
            self.of_col[i]
        }
    }

    fn informative(&self) -> bool {
        self.list
            .iter()
            .any(|f| f.members.len() >= 3 && f.stem.chars().count() >= 2)
    }

    fn label(&self, k: usize) -> String {
        let f = &self.list[k];
        format!(
            "{}{}[*]",
            if f.is_row { "" } else { "column " },
            if f.stem.is_empty() {
                "(unnamed)"
            } else {
                &f.stem
            }
        )
    }
}

// ---------------------------------------------------------------- cells

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Side {
    L,
    U,
    E,
}

impl Side {
    fn name(self) -> &'static str {
        match self {
            Side::L => "lower",
            Side::U => "upper",
            Side::E => "equality",
        }
    }
}

/// One right-hand side or column bound: the unit the data-error signal looks at.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Cell {
    is_row: bool,
    index: usize,
    side: Side,
}

fn bounds(lp: &HighsLp, is_row: bool, i: usize) -> (f64, f64) {
    if is_row {
        (lp.row_lower_[i], lp.row_upper_[i])
    } else {
        (lp.col_lower_[i], lp.col_upper_[i])
    }
}

/// The value of the `side` bound of a member, when it has one of that kind:
/// equality members `lo == up`, lower `lo` finite and `lo != up`, upper likewise.
fn side_value(lo: f64, up: f64, side: Side) -> Option<f64> {
    match side {
        Side::E => (lo == up && lo.is_finite()).then_some(lo),
        Side::L => (lo > -kHighsInf && lo != up).then_some(lo),
        Side::U => (up < kHighsInf && lo != up).then_some(up),
    }
}

fn cell_value(lp: &HighsLp, c: &Cell) -> f64 {
    let (lo, up) = bounds(lp, c.is_row, c.index);
    match c.side {
        Side::U => up,
        _ => lo,
    }
}

fn set_cell(lp: &mut HighsLp, c: &Cell, v: f64) {
    let (lo, up) = if c.is_row {
        (&mut lp.row_lower_[c.index], &mut lp.row_upper_[c.index])
    } else {
        (&mut lp.col_lower_[c.index], &mut lp.col_upper_[c.index])
    };
    match c.side {
        Side::L => *lo = v,
        Side::U => *up = v,
        Side::E => {
            *lo = v;
            *up = v;
        }
    }
}

fn side_of(lo: f64, up: f64, lower: bool, upper: bool) -> Vec<Side> {
    if lower && upper && lo == up {
        vec![Side::E]
    } else {
        let mut v = Vec::new();
        if lower {
            v.push(Side::L);
        }
        if upper {
            v.push(Side::U);
        }
        v
    }
}

// ---------------------------------------------------------------- robust statistics

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Median and modified z-score of `x` against `pop` (see the module comment).
fn robust_z(pop: &[f64], x: f64) -> (f64, f64) {
    let mut p = pop.to_vec();
    let med = median(&mut p);
    let mut dev: Vec<f64> = pop.iter().map(|v| (v - med).abs()).collect();
    let mad = median(&mut dev);
    let z = if mad > 0.0 {
        0.6745 * (x - med) / mad
    } else {
        let mean_ad = pop.iter().map(|v| (v - med).abs()).sum::<f64>() / pop.len() as f64;
        if mean_ad > 0.0 {
            (x - med) / (1.253314 * mean_ad)
        } else if x == med {
            0.0
        } else {
            f64::INFINITY.copysign(x - med)
        }
    };
    (med, z)
}

struct Outlier {
    cell: Cell,
    value: f64,
    median: f64,
    z: f64,
    family: usize,
    pop: usize,
}

/// Values of the members of family `k` that have a bound of kind `side`.
fn population(lp: &HighsLp, fams: &Families, k: usize, side: Side) -> Vec<f64> {
    let f = &fams.list[k];
    f.members
        .iter()
        .filter_map(|&i| {
            let (lo, up) = bounds(lp, f.is_row, i);
            side_value(lo, up, side)
        })
        .collect()
}

fn outlier_of(lp: &HighsLp, fams: &Families, c: &Cell) -> Option<Outlier> {
    let k = fams.of(c.is_row, c.index);
    let pop = population(lp, fams, k, c.side);
    if pop.len() < MIN_FAMILY {
        return None;
    }
    let x = cell_value(lp, c);
    let (med, z) = robust_z(&pop, x);
    // Ours: a family whose median is 0 (`row_i = 0` balances with a few nonzero
    // members) is a structural pattern, not a scale; a nonzero member there is a
    // different kind of constraint, not an out-of-scale number.
    if med == 0.0 {
        return None;
    }
    let rel = (x - med).abs() / med.abs();
    (z.abs() > Z_OUTLIER && rel >= MIN_REL_DEV).then_some(Outlier {
        cell: *c,
        value: x,
        median: med,
        z,
        family: k,
        pop: pop.len(),
    })
}

// ---------------------------------------------------------------- evidence

struct Evidence {
    id: &'static str,
    supports: &'static str,
    weight: f64,
    text: String,
}

impl Evidence {
    fn json(&self) -> Value {
        json!({"signal": self.id, "supports": self.supports, "weight": self.weight,
               "detail": self.text})
    }
}

fn cell_text(lp: &HighsLp, c: &Cell) -> String {
    format!(
        "{} {} {} bound {}",
        if c.is_row { "row" } else { "column" },
        if c.is_row {
            row_name(lp, c.index)
        } else {
            col_name(lp, c.index)
        },
        c.side.name(),
        cell_value(lp, c)
    )
}

fn fnum(x: f64) -> String {
    if x.is_infinite() {
        if x > 0.0 {
            "+inf".into()
        } else {
            "-inf".into()
        }
    } else {
        format!("{x:.4}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    }
}

fn zjson(z: f64) -> Value {
    if z.is_finite() { json!(z) } else { Value::Null }
}

pub fn run(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let a = match parse(args) {
        Ok(a) => a,
        Err(e) => return crate::diag::fail_usage_text(wj, "diagnose", &format!("{e}\n{USAGE}")),
    };
    let lp = match read_mip(&a.path) {
        Ok(lp) => lp,
        Err(e) => return crate::diag::fail_read_str(a.json, &a.path, &e),
    };
    match diagnose_core(lp, &a.path, a.time_limit) {
        Ok((out, text)) => emit(&a, out, &text),
        Err(()) => crate::diag::fail_msg(
            a.json,
            "RS0090",
            4,
            format!("diagnose: the solver failed on {}", a.path),
        ),
    }
}

/// The diagnosis of a loaded model: the `rusolve.diagnose/1` object and its text rendering.
/// `Err` when the solver fails. Shared by the CLI and the Python `Model.diagnose`.
pub(crate) fn diagnose_core(
    mut lp: HighsLp,
    path: &str,
    time_limit: f64,
) -> Result<(Value, String), ()> {
    let n_int = lp
        .integrality_
        .iter()
        .filter(|t| **t != HighsVarType::kContinuous)
        .count();
    lp.integrality_.clear();
    let kind = if n_int > 0 { "relaxation" } else { "lp" };
    let (mut h, status) = match solve_for_explanation(lp.clone()) {
        Ok(r) => r,
        Err(()) => return Err(()),
    };
    let status_text = crate::envelope::result_status(status);
    let mut out = json!({
        "schema_version": SCHEMA_VERSION,
        "input": {"path": path, "rows": lp.num_row_, "cols": lp.num_col_,
                  "integer_cols": n_int},
        "kind": kind,
        "status": status_text,
        "verdict": "none",
        "confidence": Value::Null,
        "summary": "",
        "naming_informative": Value::Null,
        "scores": {"data": 0.0, "formulation": 0.0},
        "families": [],
        "evidence": [],
        "suggested_fix": Value::Null,
        "notes": [],
    });
    if status != HighsModelStatus::kInfeasible {
        let msg =
            format!("the model is not infeasible (status {status_text}): nothing to diagnose");
        out["summary"] = json!(msg);
        return Ok((out, msg));
    }

    let mut notes: Vec<String> = Vec::new();
    if n_int > 0 {
        notes
            .push("MIP model: the diagnosis is for its LP relaxation (integrality dropped)".into());
    }

    // ---- inputs: IIS, Farkas ray, minimum-count relaxation
    let mut iis = HighsIis::default();
    let st = h.get_iis(&mut iis);
    let (irows, icols) = iis_members(&iis, st);
    if irows.is_empty() && icols.is_empty() {
        notes.push("HiGHS returned no IIS: the IIS-based signals are missing".into());
    }
    let mut has_ray = false;
    let mut ray = vec![0.0f64; lp.num_row_ as usize];
    let rst = h.get_dual_ray(&mut has_ray, Some(&mut ray));
    let y = if rst != HighsStatus::kError && has_ray {
        ray_multipliers(&ray)
    } else {
        notes.push("HiGHS returned no Farkas ray: signal (c) is missing".into());
        Vec::new()
    };
    let relaxed = match min_count_relaxation(&lp, time_limit) {
        Ok(Some((cells, ok))) => {
            if !ok {
                notes.push("the minimum-count relaxation did not verify".into());
            }
            Some(cells)
        }
        Ok(None) => {
            notes.push("the model cannot be relaxed: relaxation signals are missing".into());
            None
        }
        Err(e) => {
            notes.push(format!(
                "minimum-count relaxation failed ({e}): its signals are missing"
            ));
            None
        }
    };

    let fams = Families::build(&lp);
    let informative = fams.informative();
    if !informative {
        notes.push(
            "row/column names carry no usable family structure (no stem with 3+ members): \
             confidence reduced, family signals are weak"
                .into(),
        );
    }

    // ---- cells of the IIS and of the relaxation, per family
    let mut iis_cells: BTreeSet<Cell> = BTreeSet::new();
    let mut iis_members_of: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    // Ours: a column in the IIS only through its default lower bound 0 (x >= 0)
    // says nothing about the model's structure; such members do not count
    // towards a family spanning the IIS.
    let mut trivial_cols: BTreeSet<usize> = BTreeSet::new();
    for (is_row, ms) in [(true, &irows), (false, &icols)] {
        for m in ms {
            if !(m.lower || m.upper) {
                continue;
            }
            let (lo, up) = bounds(&lp, is_row, m.index);
            if !is_row && m.lower && !m.upper && lo == 0.0 {
                trivial_cols.insert(m.index);
            }
            for side in side_of(lo, up, m.lower, m.upper) {
                iis_cells.insert(Cell {
                    is_row,
                    index: m.index,
                    side,
                });
            }
            iis_members_of
                .entry(fams.of(is_row, m.index))
                .or_default()
                .insert(m.index);
        }
    }
    let mut rel_cells: BTreeSet<Cell> = BTreeSet::new();
    let mut rel_members_of: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for r in relaxed.iter().flatten() {
        let (lo, up) = bounds(&lp, r.is_row, r.index);
        let side = if lo == up {
            Side::E
        } else if r.lower_side {
            Side::L
        } else {
            Side::U
        };
        rel_cells.insert(Cell {
            is_row: r.is_row,
            index: r.index,
            side,
        });
        rel_members_of
            .entry(fams.of(r.is_row, r.index))
            .or_default()
            .insert(r.index);
    }

    // Ours: a bug in generator code sits in one or two constraint blocks, so a
    // formulation contradiction is localised in a few families. An IIS or
    // relaxation touching more than MAX_LOCAL_FAMILIES families is a
    // whole-model resource conflict (a staged or networked model running out of
    // capacity), which the family signals cannot tell from a data problem.
    let significant = iis_members_of
        .keys()
        .chain(rel_members_of.keys())
        .copied()
        .collect::<BTreeSet<usize>>()
        .into_iter()
        .filter(|k| {
            fams.list[*k].is_row
                || rel_members_of.contains_key(k)
                || iis_members_of[k].iter().any(|i| !trivial_cols.contains(i))
        })
        .count();
    let diffuse = significant > MAX_LOCAL_FAMILIES;

    // ---- (b) outliers among the implicated cells
    let mut outliers: Vec<Outlier> = Vec::new();
    for c in iis_cells.union(&rel_cells) {
        if let Some(o) = outlier_of(&lp, &fams, c) {
            outliers.push(o);
        }
    }
    outliers.sort_by(|p, q| q.z.abs().total_cmp(&p.z.abs()));
    let is_outlier = |c: &Cell| outliers.iter().any(|o| o.cell == *c);

    let mut ev: Vec<Evidence> = Vec::new();
    let (mut sd, mut sf) = (0.0f64, 0.0f64);
    let mut add = |ev: &mut Vec<Evidence>, id, supports, weight: f64, text: String| {
        match supports {
            "data" => sd += weight,
            "formulation" => sf += weight,
            _ => {}
        }
        ev.push(Evidence {
            id,
            supports,
            weight,
            text,
        });
    };

    // counterfactual 1: replace the outliers by the family median
    // a counterfactual that supports the verdict was run and agreed
    let mut confirmed = false;
    let mut data_fix: Option<bool> = None;
    if !outliers.is_empty() && outliers.len() <= FEW_CELLS {
        let mut fixed = lp.clone();
        for o in &outliers {
            set_cell(&mut fixed, &o.cell, o.median);
        }
        match is_infeasible(&fixed) {
            Ok(inf) => data_fix = Some(!inf),
            Err(e) => notes.push(format!("outlier counterfactual undecided: {e}")),
        }
    }
    if !outliers.is_empty() {
        let list: Vec<String> = outliers
            .iter()
            .take(8)
            .map(|o| {
                format!(
                    "{} (family median {}, robust z {})",
                    cell_text(&lp, &o.cell),
                    fnum(o.median),
                    if o.z.is_finite() {
                        format!("{:.1}", o.z)
                    } else {
                        "inf: every other member is equal".into()
                    }
                )
            })
            .collect();
        ev.push(Evidence {
            id: "b_outlier_cells",
            supports: "neutral",
            weight: 0.0,
            text: format!(
                "{} implicated cell(s) are outliers within their family: {}",
                outliers.len(),
                list.join("; ")
            ),
        });
    }
    match data_fix {
        Some(true) => add(
            &mut ev,
            "b_median_fix_restores_feasibility",
            "data",
            {
                confirmed = true;
                4.0
            },
            format!(
                "replacing the {} outlier cell(s) by their family median makes the model feasible",
                outliers.len()
            ),
        ),
        Some(false) => ev.push(Evidence {
            id: "b_median_fix_restores_feasibility",
            supports: "neutral",
            weight: 0.0,
            text:
                "replacing the outlier cell(s) by their family median leaves the model infeasible"
                    .into(),
        }),
        None => {}
    }
    if outliers.len() > FEW_CELLS {
        ev.push(Evidence {
            id: "b_many_outliers",
            supports: "neutral",
            weight: 0.0,
            text: format!(
                "{} outlier cells: too many to be 'a few numbers out of line'",
                outliers.len()
            ),
        });
    }
    if let Some(cells) = &relaxed
        && !cells.is_empty()
        && cells.len() <= FEW_CELLS
    {
        {
            let all_out = rel_cells.iter().all(&is_outlier);
            let none_out = rel_cells.iter().all(|c| !is_outlier(c));
            let list: Vec<String> = rel_cells.iter().map(|c| cell_text(&lp, c)).collect();
            if all_out {
                add(
                    &mut ev,
                    "b_relaxation_is_outliers",
                    "data",
                    1.5,
                    format!(
                        "the minimum-count relaxation gives up {} cell(s), all outliers in their family: {}",
                        list.len(),
                        list.join("; ")
                    ),
                );
            } else if none_out {
                add(
                    &mut ev,
                    "a_relaxation_is_ordinary_members",
                    "formulation",
                    if diffuse { 0.5 } else { 1.0 },
                    format!(
                        "the minimum-count relaxation gives up {} ordinary (non-outlier) cell(s): {}",
                        list.len(),
                        list.join("; ")
                    ),
                );
            }
        }
    }

    // ---- (a) family span of the IIS, coefficient consistency
    let row_coefs = {
        let mut rc: Vec<Vec<f64>> = vec![Vec::new(); lp.num_row_ as usize];
        for (i, _, v) in entries(&lp) {
            rc[i].push(v);
        }
        rc
    };
    let signature = |is_row: bool, i: usize| -> String {
        let (lo, up) = bounds(&lp, is_row, i);
        let sense = format!("{}/{}", lo > -kHighsInf, up < kHighsInf);
        if !is_row {
            return sense;
        }
        let mut d: Vec<String> = row_coefs[i].iter().map(|v| format!("{v:.9e}")).collect();
        d.sort();
        d.dedup();
        format!("{sense}|{}", d.join(","))
    };
    let mut spanning: Vec<usize> = Vec::new();
    let mut fam_rows_json: Vec<Value> = Vec::new();
    let touched: BTreeSet<usize> = iis_members_of
        .keys()
        .chain(rel_members_of.keys())
        .copied()
        .collect();
    let mut consistent_all = true;
    for &k in &touched {
        let f = &fams.list[k];
        let n_iis = iis_members_of.get(&k).map_or(0, BTreeSet::len);
        let n_rel = rel_members_of.get(&k).map_or(0, BTreeSet::len);
        let cover = n_iis as f64 / f.members.len() as f64;
        let n_span = iis_members_of.get(&k).map_or(0, |s| {
            s.iter()
                .filter(|i| f.is_row || !trivial_cols.contains(i))
                .count()
        });
        let mut consistency = Value::Null;
        if n_span >= 3 && n_span as f64 / f.members.len() as f64 >= 0.5 {
            spanning.push(k);
            let mut cnt: HashMap<String, usize> = HashMap::new();
            for &i in iis_members_of.get(&k).unwrap() {
                *cnt.entry(signature(f.is_row, i)).or_insert(0) += 1;
            }
            let modal = cnt.values().copied().max().unwrap_or(0) as f64 / n_iis as f64;
            if modal < 0.8 {
                consistent_all = false;
            }
            consistency = json!(modal);
        }
        let mut senses: BTreeMap<&str, usize> = BTreeMap::new();
        for &i in &f.members {
            let (lo, up) = bounds(&lp, f.is_row, i);
            let s = if lo == up {
                "="
            } else if lo > -kHighsInf && up < kHighsInf {
                "range"
            } else if lo > -kHighsInf {
                ">="
            } else if up < kHighsInf {
                "<="
            } else {
                "free"
            };
            *senses.entry(s).or_insert(0) += 1;
        }
        fam_rows_json.push(json!({
            "kind": if f.is_row {"row"} else {"column"},
            "stem": f.stem,
            "size": f.members.len(),
            "in_iis": n_iis,
            "in_relaxation": n_rel,
            "iis_coverage": cover,
            "coefficient_consistency": consistency,
            "senses": senses,
            "outlier_cells": outliers.iter().filter(|o| o.family == k).count(),
        }));
    }
    if diffuse {
        ev.push(Evidence {
            id: "diffuse_conflict",
            supports: "neutral",
            weight: 0.0,
            text: format!(
                "the conflict touches {significant} families (more than {MAX_LOCAL_FAMILIES}): a whole-model resource conflict, not a localised family contradiction; the family signals are not used"
            ),
        });
    }
    if spanning.len() >= 2 && !diffuse {
        let names: Vec<String> = spanning.iter().map(|&k| fams.label(k)).collect();
        add(
            &mut ev,
            "a_iis_spans_families",
            "formulation",
            1.0,
            format!(
                "the IIS covers at least half of the members of {} families ({}): a family-wide contradiction",
                names.len(),
                names.join(", ")
            ),
        );
        if consistent_all {
            add(
                &mut ev,
                "a_consistent_coefficients",
                "formulation",
                0.5,
                "within each of those families the IIS members share one coefficient pattern and sense"
                    .into(),
            );
        }
    }

    // counterfactual 2: typical data everywhere in the touched families
    let mut homogenised: Option<bool> = None;
    // Ours: the probe is only as good as "the median is a typical value"; where a
    // family's values spread widely (a time-staged or size-graded family) the
    // replacement distorts the model and the result counts for less.
    let mut heterogeneity = 0.0f64;
    if informative && !touched.is_empty() {
        let mut typical = lp.clone();
        let mut changed = false;
        for &k in &touched {
            let f = &fams.list[k];
            for side in [Side::E, Side::L, Side::U] {
                let pop = population(&lp, &fams, k, side);
                if pop.len() < MIN_FAMILY {
                    continue;
                }
                let mut p = pop.clone();
                let med = median(&mut p);
                let mut dev: Vec<f64> = pop.iter().map(|v| (v - med).abs()).collect();
                let mad = median(&mut dev);
                heterogeneity = heterogeneity.max(if med != 0.0 {
                    mad / med.abs()
                } else if mad > 0.0 {
                    f64::INFINITY
                } else {
                    0.0
                });
                for &i in &f.members {
                    let (lo, up) = bounds(&lp, f.is_row, i);
                    if side_value(lo, up, side).is_some() {
                        set_cell(
                            &mut typical,
                            &Cell {
                                is_row: f.is_row,
                                index: i,
                                side,
                            },
                            med,
                        );
                        changed = true;
                    }
                }
            }
        }
        if changed {
            match is_infeasible(&typical) {
                Ok(inf) => homogenised = Some(inf),
                Err(e) => notes.push(format!("typical-data counterfactual undecided: {e}")),
            }
        }
    }
    match homogenised {
        Some(true) => {
            let homogeneous = heterogeneity <= 0.5 && !diffuse;
            add(
                &mut ev,
                "a_typical_data_still_infeasible",
                "formulation",
                if homogeneous {
                    3.0
                } else if diffuse {
                    0.5
                } else {
                    1.0
                },
                format!(
                    "with every bound of the implicated families replaced by its family median the model is still infeasible: no few cells are to blame, the constraint pattern (or a systematic scaling) conflicts{}",
                    if homogeneous {
                        String::new()
                    } else {
                        format!(" (weakened: {})", if diffuse { "the conflict is diffuse".to_string() } else { format!("the families are heterogeneous, MAD/median up to {}", fnum(heterogeneity)) })
                    }
                ),
            );
            if homogeneous {
                confirmed = true;
            }
        }
        Some(false) => ev.push(Evidence {
            id: "a_typical_data_still_infeasible",
            supports: "neutral",
            weight: 0.0,
            text: "with every bound of the implicated families replaced by its family median the model is feasible: the conflict depends on the values"
                .into(),
        }),
        None => {}
    }

    // ---- (c) Farkas concentration
    let sum_abs: f64 = y.iter().map(|&(_, v)| v.abs()).sum();
    let sum_sq: f64 = y.iter().map(|&(_, v)| v * v).sum();
    let mut n_eff = Value::Null;
    if !y.is_empty() {
        let n = sum_abs * sum_abs / sum_sq;
        n_eff = json!(n);
        if n <= 3.0 + 1e-9 {
            add(
                &mut ev,
                "c_farkas_concentrated",
                "data",
                0.5,
                format!("the Farkas ray is carried by about {n:.1} rows: a few rows decide"),
            );
        }
        let mut by_fam: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
        for &(i, v) in &y {
            by_fam.entry(fams.of_row[i]).or_default().push(v.abs());
        }
        let blocks: Vec<String> = by_fam
            .iter()
            .filter(|(k, v)| {
                let size = fams.list[**k].members.len();
                let s: f64 = v.iter().sum();
                let q: f64 = v.iter().map(|x| x * x).sum();
                v.len() >= 3 && 2 * v.len() >= size && s * s / (v.len() as f64 * q) >= 0.8
            })
            .map(|(k, _)| fams.label(*k))
            .collect();
        if !blocks.is_empty() {
            add(
                &mut ev,
                "c_farkas_spread_over_family",
                "formulation",
                0.5,
                format!(
                    "the Farkas ray spreads near-evenly over family block(s) {}",
                    blocks.join(", ")
                ),
            );
        }
    }

    // ---- verdict
    let margin = (sd - sf).abs();
    let verdict = if margin >= DECISIVE_MARGIN {
        if sd > sf { "data" } else { "formulation" }
    } else {
        "unclear"
    };
    let mut confidence = if sd + sf > 0.0 {
        (margin / (sd + sf)).min(0.95)
    } else {
        0.0
    };
    if !informative {
        confidence *= 0.6;
    }
    // `unclear` carries no confidence: the leaning is in `scores`.
    if verdict == "unclear" {
        confidence = 0.0;
    }
    // Ours: without a confirming counterfactual the verdict rests on structure
    // alone; do not report more than moderate confidence.
    if !confirmed {
        confidence = confidence.min(0.6);
    }

    // ---- suggested fix
    let fix = match verdict {
        "data" => {
            let cells: Vec<Value> = outliers
                .iter()
                .map(|o| {
                    json!({"kind": if o.cell.is_row {"row"} else {"column"},
                           "name": if o.cell.is_row {row_name(&lp, o.cell.index)} else {col_name(&lp, o.cell.index)},
                           "bound": o.cell.side.name(), "value": o.value,
                           "family": fams.label(o.family), "family_median": o.median,
                           "family_population": o.pop,
                           "robust_z": zjson(o.z)})
                })
                .collect();
            let cells = if cells.is_empty() {
                rel_cells
                    .iter()
                    .map(|c| {
                        json!({"kind": if c.is_row {"row"} else {"column"},
                               "name": if c.is_row {row_name(&lp, c.index)} else {col_name(&lp, c.index)},
                               "bound": c.side.name(), "value": cell_value(&lp, c)})
                    })
                    .collect()
            } else {
                cells
            };
            let text = format!(
                "check the input data feeding these cell(s): {}",
                cells
                    .iter()
                    .map(|c| format!(
                        "{} {} (value {})",
                        c["kind"].as_str().unwrap_or(""),
                        c["name"].as_str().unwrap_or(""),
                        c["value"]
                    ))
                    .collect::<Vec<_>>()
                    .join("; ")
            );
            Some(json!({"type": "data", "text": text, "cells": cells, "families": []}))
        }
        "formulation" => {
            let names: Vec<String> = spanning
                .iter()
                .map(|&k| fams.label(k))
                .chain(
                    touched
                        .iter()
                        .filter(|k| !spanning.contains(k))
                        .map(|&k| fams.label(k)),
                )
                .collect();
            let text = format!(
                "check the code that generates the constraint families {} (direction of the inequalities, signs and multiplicity of coefficients, double counting); the data are not the cause",
                names.join(" and ")
            );
            Some(json!({"type": "formulation", "text": text, "cells": [], "families": names}))
        }
        _ => None,
    };

    // ---- summary and output
    let summary = match verdict {
        "data" => "the conflict looks like a data error: a few numbers out of line".to_string(),
        "formulation" => {
            "the conflict looks like a modelling error: constraint families that contradict each other"
                .to_string()
        }
        _ => "unclear: the evidence does not separate a data error from a modelling error".to_string(),
    };
    let mut text = format!(
        "diagnose {}  status {status_text}\nverdict: {verdict}  (confidence {:.2}; scores data {sd}, formulation {sf})\n{summary}\n",
        path, confidence
    );
    text.push_str("evidence:\n");
    for e in &ev {
        text.push_str(&format!(
            "  [{} {}] {}\n",
            e.supports,
            if e.weight > 0.0 {
                format!("+{}", e.weight)
            } else {
                String::new()
            },
            e.text
        ));
    }
    text.push_str("families touched by the IIS or the relaxation:\n");
    for f in &fam_rows_json {
        text.push_str(&format!(
            "  {} {}[*]: {} members, {} in IIS, {} relaxed\n",
            f["kind"].as_str().unwrap_or(""),
            f["stem"].as_str().unwrap_or(""),
            f["size"],
            f["in_iis"],
            f["in_relaxation"]
        ));
    }
    if let Some(f) = &fix {
        text.push_str(&format!(
            "suggested fix: {}\n",
            f["text"].as_str().unwrap_or("")
        ));
    }
    for n in &notes {
        text.push_str(&format!("note: {n}\n"));
    }
    out["verdict"] = json!(verdict);
    out["confidence"] = json!(confidence);
    out["summary"] = json!(summary);
    out["naming_informative"] = json!(informative);
    out["scores"] = json!({"data": sd, "formulation": sf});
    out["families"] = json!(fam_rows_json);
    out["evidence"] = json!(ev.iter().map(Evidence::json).collect::<Vec<_>>());
    out["suggested_fix"] = fix.unwrap_or(Value::Null);
    out["relaxation"] =
        relaxed.as_ref().map_or(Value::Null, |cells| {
            json!(cells
            .iter()
            .map(|r| json!({
                "kind": if r.is_row {"row"} else {"column"},
                "name": if r.is_row {row_name(&lp, r.index)} else {col_name(&lp, r.index)},
                "bound": if r.lower_side {"lower"} else {"upper"},
                "original": r.original, "relaxed_to": r.relaxed_to, "violation": r.violation,
            }))
            .collect::<Vec<_>>())
        });
    out["farkas"] = json!({"num_rows": y.len(), "effective_rows": n_eff});
    out["notes"] = json!(notes);
    Ok((out, text))
}

fn emit(a: &Args, out: Value, text: &str) -> ExitCode {
    if a.json {
        println!("{out}");
    } else {
        println!("{}", text.trim_end());
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stems_of_modelling_layer_names() {
        for (n, s) in [
            ("cap[3,5]", "cap"),
            ("cap_3_5", "cap"),
            ("cap(3)", "cap"),
            ("cap{3}", "cap"),
            ("cap<3,4>", "cap"),
            ("flow#12", "flow"),
            ("flow-12", "flow"),
            ("flow.12.4", "flow"),
            ("flow12", "flow"),
            ("R09", "R"),
            ("supply[a,b][3]", "supply"),
            ("x_1_2[7]", "x"),
            ("balance", "balance"),
            ("12", ""),
        ] {
            assert_eq!(family_stem(n), s, "{n}");
        }
    }

    #[test]
    fn symbolic_indices_group_by_shared_prefix() {
        let names: Vec<String> = [
            "supply_Rotterdam",
            "supply_Antwerp",
            "demand_Paris",
            "demand_Berlin",
            "total",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        assert_eq!(
            assign_families(&names),
            vec!["supply", "supply", "demand", "demand", "total"]
        );
    }

    #[test]
    fn robust_z_flags_a_planted_outlier() {
        let pop = [10.0, 12.0, 11.0, 13.0, 1000.0, 10.5];
        let (med, z) = robust_z(&pop, 1000.0);
        assert!(med > 10.0 && med < 13.0);
        assert!(z > Z_OUTLIER);
        // all other members equal: any other value is an infinite outlier
        let (_, z) = robust_z(&[1.0; 5], 5.0);
        assert!(z.is_infinite());
    }
}
