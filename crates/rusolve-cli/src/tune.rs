//! `rusolve tune` — tune the search options for an instance family.
//!
//! ```text
//! rusolve tune <dir | list.txt> [-o profile.json] [--format json]
//!     [--base recommended|default] [--seeds 3] [--holdout 0.25]
//!     [--budget-work 1e7] [--run-work-cap 2e6] [--time-limit 30] [--jobs 1]
//!     [--known optima.txt] [--candidates file|-] [--candidates-only]
//!     [--max-candidates 16] [--race-seed 1] [--iterations K] [--elite E] [--name NAME]
//! ```
//!
//! `<dir>` holds the family's `.mps`/`.mps.gz` files; `<list.txt>` has one path per line
//! (relative to the list, `#` comments). The target is a family, not a single instance
//! (per-family tuning plus racing, not a per-instance selector).
//!
//! ## What it does
//!
//! 1. **Split.** Instances are ordered by the SHA-256 of their file and the first
//!    `ceil(holdout * n)` become the held-out set (deterministic, mip-profiles.md 2.1). Racing
//!    never sees them.
//! 2. **Candidates.** The baseline (`--base`, default `recommended`) plus the one-exchange
//!    neighbourhood of a small search space [`SPACE`] (each option moved alone to each other
//!    value: the local move of ParamILS, Hutter, Hoos, Leyton-Brown, Stuetzle, JAIR 36, 2009),
//!    plus the candidates in `--candidates` (the LLM-proposal mode, below).
//! 3. **Racing** ([`crate::tune_race`], F-Race / irace defaults): blocks are (instance, seed)
//!    pairs in a seeded random order; every live candidate runs each block; after each block
//!    from the fifth on, a Friedman test with Conover's post-hoc drops the candidates that are
//!    significantly worse than the best at confidence 0.95. The baseline is never dropped
//!    during the race (elitist rule). Racing ends when the work budget is spent, the blocks run
//!    out, or only the baseline is left.
//! 4. **Iterated rounds** (`--iterations K`, `--elite E`; irace, López-Ibáñez, Dubois-Lacoste,
//!    Pérez Cáceres, Birattari, Stütze, "The irace package: Iterated racing for automatic
//!    algorithm configuration", Operations Research Perspectives 3 (2016) 43-58): one race is
//!    round 1. At the end of a round the candidates the final test cannot separate from the best
//!    are its survivors, ranked by mean rank; the top `E` of them (plus the baseline, which is
//!    carried into every round and never dropped) are the elite. Round r+1 races that elite again
//!    plus new candidates, each one an elite with **one more** option of [`SPACE`] changed — so
//!    round 2 searches the pair moves, round 3 the triples, and option interactions become
//!    reachable at all. New candidates are sampled from that neighbourhood with the same
//!    deterministic [`Rng`] seeded from `--race-seed`, capped at `--max-candidates` per round, and
//!    duplicates of an option set already raced are dropped. The shuffled block list is split
//!    into `K` consecutive disjoint windows so every round gets its own blocks and no block is
//!    run twice; the one `--budget-work` budget covers all rounds, and the run stops when it is
//!    spent, the blocks run out, the neighbourhood is exhausted, or K rounds are done.
//!    `--iterations 1` (the default) is exactly the single race of step 3 with no round 2.
//!    After irace (López-Ibáñez et al. 2016, cited at the round loop).
//! 5. **Choice.** After a last test on all blocks run, the surviving candidates are the ones the
//!    data cannot tell apart from the best. Among the tuned ones (not the baseline) the one
//!    that changes the fewest options from the base is the challenger; ties go to the lower
//!    mean cost. It stands against the baseline
//!    only if it is cheaper on the training blocks. The baseline itself is not in this
//!    ranking: it changes nothing, so it would always win a fewest-changes rule and the tool
//!    could never propose anything; the held-out comparison against it is the gate instead.
//!    With several rounds this is the last round's test and the last round's blocks, since a
//!    candidate first raced in round 3 has no data on round 1's blocks.
//! 6. **Hold-out.** The winner and the baseline run on the held-out instances with all seeds.
//!    The winner is written into the profile only if the 95% bootstrap interval of the
//!    geometric-mean cost ratio (winner / baseline) lies below 1.0 and no answer was wrong;
//!    otherwise the profile has no options and the evidence records the negative result
//!    (mip-profiles.md 2.7).
//!
//! ## Metric and budget
//!
//! The cost of a run is one of two measures, chosen by [`Metric::for_jobs`] and recorded as
//! `cost_metric` in the JSON output:
//!
//! * `--jobs 1` (the default): the deterministic work counter, the same `stats.work` the
//!   result envelope reports (`rusolve_ekk::work::WorkSnapshot::work`, machine-independent);
//! * `--jobs` > 1: [`work_units`] (LP iterations plus [`NODE_WORK`] per node), because the work
//!   counters are process-wide atomics and concurrent solves in one process mix into each other.
//!
//! Either way the cost is capped at `--run-work-cap` with a PAR2-style penalty of twice the cap
//! for a run that does not finish (during the solve the cap is enforced by a node limit, the
//! port's only per-run knob); the wall clock `--time-limit` is only a safety net. Ranks, the
//! reported shifted geometric mean (shift [`SHIFT`]) and the `--budget-work` budget are in the
//! chosen unit, so both flags are in that unit too. A tick is far finer than an
//! LP-iteration-equivalent ([`TICKS_PER_WORK_UNIT`]), so when the flag is absent the default of
//! each is scaled by the metric: 1e7 / 2e6 with `iters_nodes`, 2e12 / 4e11 with `c1_work`, which
//! budgets about the same effort on either. Pass `--budget-work` and `--run-work-cap` explicitly
//! when you know the family, and read `cost_metric` to know what the numbers mean.
//!
//! ## Wrong-answer check
//!
//! A run that reports `Optimal` with an objective that differs from the known optimum
//! (relative 1e-6), or `Infeasible`/`Unbounded` when an optimum is known, disqualifies the
//! candidate at once. `--known` is a file of `name value` lines (also `name,value`, and the
//! MIPLIB `.solu` lines `=opt= name value`). An instance without a known value is checked
//! against the first optimal objective the baseline reported for it, and that is recorded in
//! the evidence (`source: baseline_consensus`).
//!
//! ## LLM-proposal mode
//!
//! `--candidates file|-` reads JSON from a file or stdin: an array, or `{"candidates": [...]}`,
//! of `{"name": "...", "options": {"mip_heuristic_effort": 0.15, ...}}` (a profile file works
//! too). An agent writes it after reading `rusolve probe` / `classify --format json` and a
//! solve log. Each entry is validated like a profile file (errors, exit 2) and then raced
//! like any other candidate, so a wrong guess costs budget, never correctness. The binary calls
//! no model API. `--candidates-only` races the baseline and the proposals alone.
//!
//! ## Exit codes
//!
//! 0 done (also when no candidate beat the baseline: that is a result); 2 usage error or an
//! invalid candidate file; 3 an input file is missing, unreadable or does not parse; 4 internal.

use crate::profile::{
    self, Base, Fingerprint, Profile, parse_option_text, validate_option, value_json,
};
use crate::tune_race::{self as race, Rng};
use rusolve_ekk::HighsInt;
use rusolve_ekk::lp_data::hconst::HighsModelStatus;
use rusolve_ekk::lp_data::highs_lp::HighsLp;
use rusolve_ekk::lp_data::highs_option_records::OptionValue;
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::lp_data::model_utils::util_model_status_to_string;
use rusolve_ekk::prelude::HighsStatus;
use serde_json::{Value, json};
use std::collections::{BTreeMap, HashMap};
use std::process::ExitCode;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

pub const TUNE_SCHEMA: &str = "rusolve.tune/1";
const USAGE: &str = "usage: rusolve tune <dir|list.txt> [-o profile.json] [--format json] \
[--base recommended|default] [--seeds N] [--holdout F] [--budget-work W] [--run-work-cap W] \
[--time-limit S] [--jobs N] [--known file] [--candidates file|-] [--candidates-only] \
[--max-candidates N] [--race-seed N] [--iterations K] [--elite E] [--name NAME]\n       \
    rusolve tune --space     (print the searched options as JSON; no instances needed)";

/// Cost weight of one branch-and-bound node, in LP-iteration equivalents. Only the
/// [`Metric::ItersNodes`] metric and the per-run node cap use it.
const NODE_WORK: f64 = 10.0;
/// Shift of the shifted geometric mean, in work units (about one small root solve).
const SHIFT: f64 = 1000.0;
const ALPHA: f64 = 0.05;
const FIRST_TEST: usize = 5;
const WRONG_TOL: f64 = 1e-6;

/// How many work ticks one [`work_units`] costs, roughly. Measured on the `tests/tune.rs`
/// family (18-column knapsacks, base `recommended`, seed 0): 2.9e7 `stats.work` per solve against
/// 145 LP-iteration-equivalents, so about 2e5. The work counter counts touched nonzeros, the placeholder counts
/// iterations, so the two differ by orders of magnitude.
///
/// Used for exactly two things, never inside a cost: the defaults of `--budget-work` and
/// `--run-work-cap` (so a default run budgets about the same effort on either metric), and
/// undoing the metric when the node limit that enforces the run cap inside the solve is derived
/// ([`Cfg::node_cap`], which is therefore the same node count on either metric).
const TICKS_PER_WORK_UNIT: f64 = 2.0e5;

/// The cost of a finished run in the [`Metric::ItersNodes`] unit: LP iterations plus
/// [`NODE_WORK`] per node. Cheap and safe under concurrency, but only a proxy for the work
/// actually done; see [`Metric`] for when the work counter replaces it.
pub fn work_units(lp_iterations: i64, nodes: i64) -> f64 {
    lp_iterations as f64 + NODE_WORK * nodes as f64
}

/// Which number a run's cost is. The work counters (`rusolve_ekk::work`) are process-wide atomics,
/// so they are only attributable to one solve when solves do not overlap: `--jobs 1`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Metric {
    /// `stats.work` of the envelope: simplex synthetic ticks plus the non-simplex units at
    /// `NONSIMPLEX_UNIT_TICKS` each (`rusolve_ekk::work::WorkSnapshot::work`). Deterministic
    /// and machine-independent.
    C1Work,
    /// [`work_units`]: LP iterations plus a constant per node.
    ItersNodes,
}

impl Metric {
    pub fn for_jobs(jobs: usize) -> Self {
        if jobs <= 1 {
            Self::C1Work
        } else {
            Self::ItersNodes
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::C1Work => "c1_work",
            Self::ItersNodes => "iters_nodes",
        }
    }
    pub fn definition(self) -> String {
        match self {
            Self::C1Work => format!(
                "stats.work of rusolve.result/1: simplex synthetic ticks + {} per propagation, \
                 presolve and FeasibilityJump unit (rusolve-ekk/src/work.rs)",
                rusolve_ekk::work::NONSIMPLEX_UNIT_TICKS
            ),
            Self::ItersNodes => {
                format!("lp_iterations + {NODE_WORK} * nodes (work counter unavailable)")
            }
        }
    }
}

/// The one-exchange search space: search options of `rusolve options` and the values to try
/// (the base's own value is skipped). Every entry is checked against the registry by the
/// `tune_space_is_valid` test. Options with a documented effect on the tree search or the
/// primal heuristics (HiGHS 1.16.0 option descriptions), not the rusolve levers, which are
/// evaluated separately.
pub const SPACE: &[(&str, &[&str])] = &[
    ("mip_heuristic_effort", &["0.01", "0.15", "0.3"]),
    ("mip_pscost_minreliable", &["0", "4", "16"]),
    ("mip_lp_age_limit", &["5", "20"]),
    ("mip_pool_age_limit", &["10", "60"]),
    ("mip_detect_symmetry", &["true", "false"]),
    ("mip_allow_restart", &["true", "false"]),
    ("mip_allow_cut_separation_at_nodes", &["true", "false"]),
    ("mip_heuristic_run_feasibility_jump", &["true", "false"]),
    ("mip_heuristic_run_rins", &["true", "false"]),
    ("mip_heuristic_run_rens", &["true", "false"]),
    ("mip_heuristic_run_root_reduced_cost", &["true", "false"]),
    ("mip_heuristic_run_zi_round", &["true", "false"]),
    ("mip_heuristic_run_shifting", &["true", "false"]),
    ("mip_conflict_resolution", &["true", "false"]),
    ("mip_cut_selection_hybrid", &["true", "false"]),
    ("mip_symmetry_dynamic", &["true", "false"]),
];

struct Cfg {
    input: String,
    out: Option<String>,
    json: bool,
    base: Base,
    seeds: usize,
    holdout: f64,
    budget: Option<f64>,
    run_cap: Option<f64>,
    time_limit: f64,
    jobs: usize,
    known: Option<String>,
    candidates: Option<String>,
    candidates_only: bool,
    max_candidates: usize,
    race_seed: u64,
    name: Option<String>,
    /// `--iterations`: racing rounds (irace's `nbIterations`). 1 is the single race.
    iterations: usize,
    /// `--elite`: how many survivors of a round are carried on and mutated (irace's elites).
    elite: usize,
}

impl Cfg {
    /// The metric this run scores with.
    fn metric(&self) -> Metric {
        Metric::for_jobs(self.jobs)
    }
    /// How many ticks one [`work_units`] is, for turning the flag defaults into metric units.
    fn scale(&self) -> f64 {
        if self.metric() == Metric::C1Work {
            TICKS_PER_WORK_UNIT
        } else {
            1.0
        }
    }
    /// `--budget-work` in the metric's unit.
    fn budget(&self) -> f64 {
        self.budget.unwrap_or(1e7 * self.scale())
    }
    /// `--run-work-cap` in the metric's unit.
    fn run_cap(&self) -> f64 {
        self.run_cap.unwrap_or(2e6 * self.scale())
    }
    /// The node limit that enforces `--run-work-cap` inside the solve. The cap is divided back
    /// out of the metric first, so the same nominal cap is the same node count whichever metric
    /// scores the run (the port has no work-unit limit, only `mip_max_nodes`). Never below 100
    /// nodes: that is what the 1000-unit minimum of `--run-work-cap` has always meant, so a run
    /// always searches a little before it is capped, however the cap is worded.
    fn node_cap(&self) -> HighsInt {
        (((self.run_cap() / self.scale()) / NODE_WORK).clamp(100.0, i32::MAX as f64)) as HighsInt
    }
}

/// A named set of options an agent proposed.
type Proposal = (String, BTreeMap<String, OptionValue>);

struct Inst {
    name: String,
    path: String,
    sha256: String,
    lp: HighsLp,
}

#[derive(Clone, Debug)]
struct Run {
    status: HighsModelStatus,
    objective: Option<f64>,
    lp_iterations: i64,
    nodes: i64,
    error: bool,
    /// `stats.work` of this solve, when the metric can attribute it (see [`Metric`]).
    c1_work: Option<u64>,
}

impl Run {
    fn terminal(&self) -> bool {
        !self.error
            && matches!(
                self.status,
                HighsModelStatus::kOptimal
                    | HighsModelStatus::kInfeasible
                    | HighsModelStatus::kUnbounded
                    | HighsModelStatus::kUnboundedOrInfeasible
            )
    }
    fn work(&self, metric: Metric) -> f64 {
        match metric {
            // A run that produced no MIP stats (an LP-only or failed solve) has no count.
            Metric::C1Work => self
                .c1_work
                .map_or_else(|| work_units(self.lp_iterations, self.nodes), |w| w as f64),
            Metric::ItersNodes => work_units(self.lp_iterations, self.nodes),
        }
    }
    /// The racing cost: work, capped; PAR2 penalty when the run did not finish.
    fn cost(&self, cap: f64, metric: Metric) -> f64 {
        if self.terminal() {
            self.work(metric).min(cap)
        } else {
            2.0 * cap
        }
    }
}

#[derive(Clone, PartialEq)]
enum State {
    Alive,
    Eliminated(usize),
    /// Survived round `r` but was not taken into its elite: it is not raced again.
    Held(usize),
    Disqualified(String),
}

impl State {
    /// The state as the JSON output names it.
    fn json_name(&self) -> String {
        match self {
            State::Alive => "survived".to_string(),
            State::Eliminated(b) => format!("eliminated_after_block_{b}"),
            State::Held(r) => format!("survived_round_{r}_not_elite"),
            State::Disqualified(r) => format!("disqualified: {r}"),
        }
    }
    /// The state as the text report names it.
    fn text_name(&self) -> String {
        match self {
            State::Alive => "survived".to_string(),
            State::Eliminated(b) => format!("dropped after block {b}"),
            State::Held(r) => format!("survived round {r}, not elite"),
            State::Disqualified(r) => format!("DISQUALIFIED {r}"),
        }
    }
}

struct Cand {
    name: String,
    origin: &'static str,
    /// 1-based racing round this candidate first appeared in (irace's iterations).
    round: usize,
    delta: BTreeMap<String, OptionValue>,
    effective: HighsOptions,
    state: State,
    /// cost per block index, `None` where the candidate was not run
    costs: Vec<Option<f64>>,
}

/// One racing round: the window of blocks it ran, who raced, who survived its final test,
/// which of those became the elite that seeds the next round, and why it stopped.
struct Round {
    index: usize,
    /// candidate indices that raced in this round
    raced: Vec<usize>,
    /// the round's final noise group: survivors the test cannot separate from the best
    survivors: Vec<usize>,
    /// survivors carried into the next round (the baseline first, then the top `--elite`)
    elites: Vec<usize>,
    /// mean rank of each entry of `survivors` in this round's blocks (1 = cheapest)
    mean_rank: Vec<f64>,
    /// block indices this round owned: `[from, to)` of the shuffled block list
    from: usize,
    to: usize,
    final_p: Option<f64>,
    stopped: &'static str,
    new_candidates: usize,
}

impl Cand {
    fn changes(&self) -> usize {
        self.delta.len()
    }
    fn evaluated(&self) -> usize {
        self.costs.iter().flatten().count()
    }
    /// Shifted geometric mean of the costs over the blocks it ran.
    fn sgm(&self) -> f64 {
        sgm(self.costs.iter().flatten().copied())
    }
}

impl Cand {
    /// Cost relative to the baseline on exactly the blocks this candidate ran (the shifted
    /// geometric means of the two, over the same blocks): below 1.0 is cheaper.
    fn ratio_to(&self, base: &Cand) -> f64 {
        let (mut c, mut b) = (Vec::new(), Vec::new());
        for (x, y) in self.costs.iter().zip(&base.costs) {
            if let (Some(x), Some(y)) = (x, y) {
                c.push(*x);
                b.push(*y);
            }
        }
        (sgm(c.into_iter()) + SHIFT) / (sgm(b.into_iter()) + SHIFT)
    }
}

fn sgm(costs: impl Iterator<Item = f64>) -> f64 {
    let (mut s, mut n) = (0.0, 0usize);
    for c in costs {
        s += (c + SHIFT).ln();
        n += 1;
    }
    if n == 0 {
        f64::NAN
    } else {
        (s / n as f64).exp() - SHIFT
    }
}

struct Wrong {
    candidate: String,
    instance: String,
    seed: u64,
    status: String,
    objective: Option<f64>,
    expected: f64,
    source: &'static str,
}

fn usage_err(json: bool, msg: &str) -> ExitCode {
    crate::diag::fail_usage_text(json, "tune", msg)
}

fn stem(path: &str) -> String {
    let b = path.rsplit('/').next().unwrap_or(path);
    let b = b.strip_suffix(".gz").unwrap_or(b);
    b.strip_suffix(".mps").unwrap_or(b).to_string()
}

fn parse_known(text: &str) -> HashMap<String, f64> {
    let mut m = HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let toks: Vec<&str> = line
            .split(|c: char| c.is_whitespace() || c == ',')
            .filter(|t| !t.is_empty())
            .collect();
        let (name, val) = match toks.as_slice() {
            ["=opt=", name, val, ..] => (*name, *val),
            [tag, ..] if tag.starts_with('=') => continue,
            [name, val, ..] => (*name, *val),
            _ => continue,
        };
        if let Ok(v) = val.parse::<f64>() {
            m.insert(stem(name), v);
        }
    }
    m
}

fn load_instances(input: &str) -> Result<Vec<Inst>, (u8, String)> {
    let mut paths: Vec<String> = Vec::new();
    let meta = std::fs::metadata(input).map_err(|e| (3, format!("{input}: {e}")))?;
    if meta.is_dir() {
        let rd = std::fs::read_dir(input).map_err(|e| (3, format!("{input}: {e}")))?;
        for e in rd.flatten() {
            let p = e.path().to_string_lossy().to_string();
            if p.ends_with(".mps") || p.ends_with(".mps.gz") {
                paths.push(p);
            }
        }
        paths.sort();
    } else {
        let text = std::fs::read_to_string(input).map_err(|e| (3, format!("{input}: {e}")))?;
        let dir = std::path::Path::new(input)
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_default();
        for line in text.lines().map(str::trim) {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let p = std::path::Path::new(line);
            paths.push(if p.is_absolute() {
                line.to_string()
            } else {
                dir.join(p).to_string_lossy().to_string()
            });
        }
    }
    if paths.is_empty() {
        return Err((3, format!("{input}: no .mps or .mps.gz instances found")));
    }
    paths
        .into_iter()
        .map(|p| {
            let lp = crate::mip_highs::read_mip(&p).map_err(|e| (3, e))?;
            let sha256 =
                crate::envelope::sha256_of_file(&p).map_err(|e| (3, format!("{p}: {e}")))?;
            Ok(Inst {
                name: stem(&p),
                path: p,
                sha256,
                lp,
            })
        })
        .collect()
}

fn run_once(inst: &Inst, effective: &HighsOptions, seed: u64, cfg: &Cfg) -> Run {
    let mut o = effective.clone();
    o.output_flag = false;
    o.threads = 1;
    o.random_seed = seed as HighsInt;
    o.time_limit = cfg.time_limit;
    // comparable, checkable answers: the same gap settings for every run (mip_perf.py does the same)
    o.mip_rel_gap = 0.0;
    o.mip_abs_gap = 1e-9;
    o.mip_max_nodes = cfg.node_cap();
    o.sync_log_options();
    let mut h = Highs::with_options(o);
    let failed =
        h.pass_model(inst.lp.clone()) == HighsStatus::kError || h.run() == HighsStatus::kError;
    let info = h.get_info();
    // The work of this solve: `Highs::mip_stats_.work`, the same value `envelope.rs` prints as
    // `stats.work` (`mip_highs.rs` hands `highs.mip_stats_.as_ref()` to `envelope::Solved`). Read
    // only with `--jobs 1`: the counters behind it are process-wide atomics, so solves that
    // overlap in this process cannot be told apart.
    let c1_work = if cfg.metric() == Metric::C1Work {
        h.mip_stats_.as_ref().map(|m| m.work)
    } else {
        None
    };
    Run {
        status: h.get_model_status(),
        objective: h
            .solution_
            .value_valid
            .then_some(info.objective_function_value),
        lp_iterations: info.simplex_iteration_count as i64,
        nodes: info.mip_node_count,
        error: failed,
        c1_work,
    }
}

/// Run `f` over `items` on `jobs` threads, results in input order.
fn par_map<T: Sync, R: Send>(items: &[T], jobs: usize, f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<Option<R>>> = Mutex::new((0..items.len()).map(|_| None).collect());
    std::thread::scope(|s| {
        for _ in 0..jobs.max(1).min(items.len().max(1)) {
            s.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::SeqCst);
                    if i >= items.len() {
                        break;
                    }
                    let r = f(&items[i]);
                    out.lock().unwrap()[i] = Some(r);
                }
            });
        }
    });
    out.into_inner()
        .unwrap()
        .into_iter()
        .map(|r| r.expect("every item ran"))
        .collect()
}

fn make_cand(
    name: String,
    origin: &'static str,
    round: usize,
    delta: BTreeMap<String, OptionValue>,
    base: Base,
    nblocks: usize,
) -> Result<Cand, String> {
    let mut p = Profile::new(base);
    p.options = delta.clone();
    let effective = p.effective()?;
    Ok(Cand {
        name,
        origin,
        round,
        delta,
        effective,
        state: State::Alive,
        costs: vec![None; nblocks],
    })
}

/// Candidates from the proposal file: validated like a profile, every error reported.
fn read_proposals(path: &str, base: Base) -> Result<Vec<Proposal>, (u8, String)> {
    let text = if path == "-" {
        let mut s = String::new();
        std::io::Read::read_to_string(&mut std::io::stdin(), &mut s)
            .map_err(|e| (3, format!("stdin: {e}")))?;
        s
    } else {
        std::fs::read_to_string(path).map_err(|e| (3, format!("{path}: {e}")))?
    };
    let v: Value =
        serde_json::from_str(&text).map_err(|e| (3, format!("{path}: not JSON: {e}")))?;
    let list = match &v {
        Value::Array(a) => a.clone(),
        Value::Object(o) if o.get("candidates").is_some_and(Value::is_array) => {
            o["candidates"].as_array().unwrap().clone()
        }
        _ => {
            return Err((
                2,
                format!("{path}: expected an array or {{\"candidates\": [...]}}"),
            ));
        }
    };
    let mut errors = Vec::new();
    let mut out = Vec::new();
    for (i, c) in list.iter().enumerate() {
        let at = format!("{path}: candidates[{i}]");
        let Some(opts) = c.get("options").and_then(Value::as_object) else {
            errors.push(format!("{at}.options: expected an object of name: value"));
            continue;
        };
        if let Some(b) = c.get("base").and_then(Value::as_str)
            && b != base.name()
        {
            errors.push(format!(
                "{at}.base: '{b}' differs from the tuning base '{}'",
                base.name()
            ));
        }
        let mut delta = BTreeMap::new();
        for (k, val) in opts {
            match validate_option(k, val) {
                Ok(ov) => {
                    delta.insert(k.clone(), ov);
                }
                Err(e) => errors.push(format!("{at}.{e}")),
            }
        }
        let name = c
            .get("name")
            .and_then(Value::as_str)
            .map(str::to_string)
            .unwrap_or_else(|| format!("proposal-{i}"));
        out.push((name, delta));
    }
    if errors.is_empty() {
        Ok(out)
    } else {
        Err((2, errors.join("\n")))
    }
}

fn iso_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0) as i64;
    // civil-from-days (H. Hinnant, "chrono-Compatible Low-Level Date Algorithms")
    let (days, rem) = (secs.div_euclid(86400), secs.rem_euclid(86400));
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

fn parse_args(args: &[String]) -> Result<Cfg, String> {
    let mut cfg = Cfg {
        input: String::new(),
        out: None,
        json: false,
        base: Base::Recommended,
        seeds: 3,
        holdout: 0.25,
        budget: None,
        run_cap: None,
        time_limit: 30.0,
        jobs: 1,
        known: None,
        candidates: None,
        candidates_only: false,
        max_candidates: 16,
        race_seed: 1,
        name: None,
        iterations: 1,
        elite: 3,
    };
    let mut it = args.iter();
    let need = |it: &mut std::slice::Iter<String>, flag: &str| -> Result<String, String> {
        it.next().cloned().ok_or(format!("{flag} needs a value"))
    };
    let num = |s: String, flag: &str| -> Result<f64, String> {
        s.parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v >= 0.0)
            .ok_or(format!("{flag}: '{s}' is not a non-negative number"))
    };
    // a count that must be at least 1: `--iterations 0` and `--elite 0` are usage errors, not
    // silently-clamped no-ops (an error names the bad input and the fix)
    let count = |s: String, flag: &str| -> Result<usize, String> {
        s.parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v >= 1.0)
            .map(|v| v as usize)
            .ok_or(format!("{flag}: '{s}' is not a whole number of at least 1"))
    };
    while let Some(a) = it.next() {
        match a.as_str() {
            "-o" | "--output" => cfg.out = Some(need(&mut it, a)?),
            "--format" => match need(&mut it, a)?.as_str() {
                "json" => cfg.json = true,
                "text" => {}
                _ => return Err("--format needs json or text".into()),
            },
            "--base" => {
                cfg.base =
                    Base::parse(&need(&mut it, a)?).ok_or("--base needs default or recommended")?
            }
            "--seeds" => cfg.seeds = num(need(&mut it, a)?, a)?.max(1.0) as usize,
            "--holdout" => cfg.holdout = num(need(&mut it, a)?, a)?.min(0.5),
            "--budget-work" => cfg.budget = Some(num(need(&mut it, a)?, a)?),
            "--run-work-cap" => cfg.run_cap = Some(num(need(&mut it, a)?, a)?.max(1000.0)),
            "--time-limit" => cfg.time_limit = num(need(&mut it, a)?, a)?,
            "--jobs" => cfg.jobs = num(need(&mut it, a)?, a)?.max(1.0) as usize,
            "--known" => cfg.known = Some(need(&mut it, a)?),
            "--candidates" => cfg.candidates = Some(need(&mut it, a)?),
            "--candidates-only" => cfg.candidates_only = true,
            "--max-candidates" => cfg.max_candidates = num(need(&mut it, a)?, a)? as usize,
            "--race-seed" => cfg.race_seed = num(need(&mut it, a)?, a)? as u64,
            "--iterations" => cfg.iterations = count(need(&mut it, a)?, a)?,
            "--elite" => cfg.elite = count(need(&mut it, a)?, a)?,
            "--name" => cfg.name = Some(need(&mut it, a)?),
            "--help" | "-h" => return Err(String::new()),
            s if s.starts_with('-') && s != "-" => {
                return Err(format!("flag '{s}' is not supported"));
            }
            _ if cfg.input.is_empty() => cfg.input = a.clone(),
            _ => return Err(format!("unexpected argument '{a}'")),
        }
    }
    if cfg.input.is_empty() {
        return Err("missing <dir|list.txt>".into());
    }
    if cfg.candidates_only && cfg.candidates.is_none() {
        return Err("--candidates-only needs --candidates".into());
    }
    Ok(cfg)
}

pub fn run(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    // `--space`: the options a race mutates, for an agent that proposes candidates (needs no
    // instances, runs nothing).
    if args.iter().any(|a| a == "--space") {
        let opts: Vec<Value> = SPACE
            .iter()
            .map(|(n, vals)| json!({"name": n, "values": vals}))
            .collect();
        println!(
            "{}",
            json!({"schema_version": "rusolve.tune-space/1", "iterations_note":
                "one-exchange candidates change one option to one of `values` (the base's own \
                 value is skipped); later rounds change more options at once; options outside \
                 this list reach a race through --candidates",
                "options": opts})
        );
        return ExitCode::SUCCESS;
    }
    let cfg = match parse_args(args) {
        Ok(c) => c,
        Err(m) if m.is_empty() => {
            println!("{USAGE}");
            return ExitCode::SUCCESS;
        }
        Err(m) => return usage_err(wj, &format!("{m}\n{USAGE}")),
    };
    match tune(&cfg) {
        Ok(()) => ExitCode::SUCCESS,
        Err((code, msg)) => crate::diag::fail_coded(cfg.json, "tune", code, &msg),
    }
}

fn tune(cfg: &Cfg) -> Result<(), (u8, String)> {
    let mut insts = load_instances(&cfg.input)?;
    let mut known: HashMap<String, f64> = HashMap::new();
    if let Some(k) = &cfg.known {
        known = parse_known(&std::fs::read_to_string(k).map_err(|e| (3, format!("{k}: {e}")))?);
    }
    // split by hash order (mip-profiles.md 2.1)
    insts.sort_by(|a, b| a.sha256.cmp(&b.sha256));
    let n = insts.len();
    let mut n_hold = ((cfg.holdout * n as f64).ceil() as usize).min(n / 2);
    if n < 4 {
        n_hold = 0;
    }
    let (hold, train) = insts.split_at(n_hold);
    if train.len() < 2 {
        return Err((3, format!("{n} instance(s): at least 2 are needed to tune")));
    }
    if hold.is_empty() {
        eprintln!(
            "rusolve tune: warning: no held-out instance ({n} instances); \
             the result cannot be confirmed and no option will be written"
        );
    }
    let seeds: Vec<u64> = (0..cfg.seeds as u64).collect();

    // ---- candidates
    let metric = cfg.metric();
    let mut blocks: Vec<(usize, u64)> = (0..train.len())
        .flat_map(|i| seeds.iter().map(move |&s| (i, s)))
        .collect();
    Rng(cfg.race_seed).shuffle(&mut blocks);
    let nb = blocks.len();
    let mut cands: Vec<Cand> = vec![
        make_cand(
            format!("base:{}", cfg.base.name()),
            "baseline",
            1,
            BTreeMap::new(),
            cfg.base,
            nb,
        )
        .map_err(|e| (4, e))?,
    ];
    let base_opts = cfg.base.options();
    let mut seen: Vec<BTreeMap<String, OptionValue>> = vec![BTreeMap::new()];
    let mut add = |cands: &mut Vec<Cand>,
                   name: String,
                   origin: &'static str,
                   round: usize,
                   mut delta: BTreeMap<String, OptionValue>|
     -> Result<bool, String> {
        // an option set to the base's own value changes nothing: drop it
        delta.retain(|k, v| {
            profile::value_json(v)
                != profile::value_json(
                    &rusolve_ekk::lp_data::highs_option_records::get_option_value(&base_opts, k)
                        .expect("validated option"),
                )
        });
        if seen.contains(&delta) {
            return Ok(false);
        }
        seen.push(delta.clone());
        cands.push(make_cand(name, origin, round, delta, cfg.base, nb)?);
        Ok(true)
    };
    if let Some(p) = &cfg.candidates {
        for (name, delta) in read_proposals(p, cfg.base)? {
            add(&mut cands, name, "proposal", 1, delta)
                .map_err(|e| (2, e))
                .map(|_| ())?;
        }
    }
    if !cfg.candidates_only {
        let mut moves: Vec<(String, BTreeMap<String, OptionValue>)> = Vec::new();
        for (name, values) in SPACE {
            for v in *values {
                let ov = parse_option_text(name, v).map_err(|e| (4, e))?;
                let cur =
                    rusolve_ekk::lp_data::highs_option_records::get_option_value(&base_opts, name);
                if cur.as_ref() == Some(&ov) {
                    continue;
                }
                moves.push((
                    format!("{name}={v}"),
                    BTreeMap::from([(name.to_string(), ov)]),
                ));
            }
        }
        if moves.len() > cfg.max_candidates {
            Rng(cfg.race_seed ^ 0xA5A5).shuffle(&mut moves);
            moves.truncate(cfg.max_candidates);
        }
        for (name, delta) in moves {
            add(&mut cands, name, "one_exchange", 1, delta)
                .map(|_| ())
                .map_err(|e| (4, e))?;
        }
    }
    if cands.len() < 2 {
        return Err((2, "no candidate besides the baseline to race".into()));
    }

    // ---- iterated racing (irace). Round 1 races the candidates built above; every later round
    // races the elite of the previous round plus candidates that change one more option of SPACE.
    // Citation: M. López-Ibáñez, S. Dubois-Lacoste, L. Pérez Cáceres, T. Birattari, T. Stützle,
    // "The irace package: Iterated racing for automatic algorithm configuration", Operations
    // Research Perspectives 3 (2016) 43-58.
    let mut consensus: HashMap<String, f64> = HashMap::new();
    let mut wrong: Vec<Wrong> = Vec::new();
    let (mut spent, mut evaluations) = (0.0f64, 0usize);
    let mut drop_log: Vec<Value> = Vec::new();
    let mut tests = 0usize;
    let mut rounds: Vec<Round> = Vec::new();
    // The shuffled blocks are split into `iterations` disjoint consecutive windows, so each round
    // owns its blocks and no block is ever run twice. K = 1 makes the window the whole list,
    // which is the single race of the previous behaviour.
    let window = nb.div_ceil(cfg.iterations);
    let mut stop_reason = "iterations_done";
    for round_index in 1..=cfg.iterations {
        let from = (round_index - 1) * window;
        let to = (from + window).min(nb);
        if from >= to {
            stop_reason = "blocks_exhausted";
            break;
        }
        let raced: Vec<usize> = (0..cands.len())
            .filter(|&c| cands[c].state == State::Alive)
            .collect();
        let mut stopped = "round_complete";
        let mut bi = from;
        while bi < to {
            if spent >= cfg.budget() {
                stopped = "work_budget_spent";
                break;
            }
            let live: Vec<usize> = (0..cands.len())
                .filter(|&c| cands[c].state == State::Alive)
                .collect();
            if live.len() == 1 {
                stopped = "only_the_baseline_survives";
                break;
            }
            let (ii, seed) = blocks[bi];
            let inst = &train[ii];
            let results = par_map(&live, cfg.jobs, |&c| {
                run_once(inst, &cands[c].effective, seed, cfg)
            });
            evaluations += live.len();
            // the baseline is `live[0]` (candidate 0 is never dropped): it sets the consensus first
            for (&c, r) in live.iter().zip(&results) {
                spent += r.work(metric);
                if c == 0
                    && r.status == HighsModelStatus::kOptimal
                    && !known.contains_key(&inst.name)
                {
                    consensus
                        .entry(inst.name.clone())
                        .or_insert(r.objective.unwrap_or(f64::NAN));
                }
                if let Some(w) = check_answer(&cands[c].name, inst, seed, r, &known, &consensus) {
                    if c != 0 {
                        cands[c].state = State::Disqualified(format!(
                            "wrong answer on {} (seed {seed}): {} {:?}, expected {}",
                            w.instance, w.status, w.objective, w.expected
                        ));
                    }
                    wrong.push(w);
                }
                cands[c].costs[bi] = Some(r.cost(cfg.run_cap(), metric));
            }
            bi += 1;
            eprintln!(
                "rusolve tune: round {round_index} block {bi}/{nb} ({} seed {seed}): {} live, work spent {:.3e} of {:.3e}",
                inst.name,
                cands.iter().filter(|c| c.state == State::Alive).count(),
                spent,
                cfg.budget()
            );
            // Friedman test on this round's blocks so far (F-Race / irace firstTest, eachTest = 1)
            let alive: Vec<usize> = (0..cands.len())
                .filter(|&c| cands[c].state == State::Alive)
                .collect();
            if bi - from >= FIRST_TEST && alive.len() > 2 {
                tests += 1;
                let table = table_of(&cands, &alive, from, bi);
                for col in race::race_drops(&table, ALPHA) {
                    let c = alive[col];
                    if c != 0 {
                        cands[c].state = State::Eliminated(bi);
                        drop_log.push(json!({"candidate": cands[c].name, "after_blocks": bi,
                                             "round": round_index}));
                    }
                }
            }
        }
        // ---- the round's last test: the survivors are the ones the data cannot separate
        let alive: Vec<usize> = (0..cands.len())
            .filter(|&c| cands[c].state == State::Alive)
            .collect();
        let table = table_of(&cands, &alive, from, bi);
        let final_test = race::friedman(&table, ALPHA);
        let mut group = alive.clone();
        if let Some(f) = &final_test
            && f.p_value < ALPHA
        {
            let best = f.rank_sums.iter().cloned().fold(f64::INFINITY, f64::min);
            group = (0..alive.len())
                .filter(|&j| f.rank_sums[j] - best <= f.critical_difference)
                .map(|j| alive[j])
                .collect();
        }
        let mean_ranks = race::mean_ranks(&table);
        let rank_of = |c: usize| -> f64 {
            alive
                .iter()
                .position(|&a| a == c)
                .map_or(f64::INFINITY, |j| mean_ranks[j])
        };
        rounds.push(Round {
            index: round_index,
            raced,
            survivors: group.clone(),
            elites: vec![],
            mean_rank: group.iter().map(|&c| rank_of(c)).collect(),
            from,
            to: bi,
            final_p: final_test.as_ref().map(|f| f.p_value),
            stopped,
            new_candidates: 0,
        });
        if round_index == cfg.iterations {
            stop_reason = "iterations_done";
            break;
        }
        if stopped == "work_budget_spent" || spent >= cfg.budget() {
            stop_reason = "work_budget_spent";
            break;
        }
        if to >= nb {
            stop_reason = "blocks_exhausted";
            break;
        }
        // ---- the elite: the baseline (always carried, never dropped) plus the best survivors by
        // mean rank. Everything else that was still alive is held out of the later rounds.
        let mut elite: Vec<usize> = group.iter().copied().filter(|&c| c != 0).collect();
        elite.sort_by(|&a, &b| rank_of(a).total_cmp(&rank_of(b)).then(a.cmp(&b)));
        elite.truncate(cfg.elite);
        for &c in &alive {
            if c != 0 && !elite.contains(&c) {
                cands[c].state = State::Held(round_index);
            }
        }
        // ---- the next round's candidates: each elite with one more option of SPACE changed
        let mut moves: Vec<(String, BTreeMap<String, OptionValue>)> = Vec::new();
        for &e in &elite {
            moves.extend(elite_neighbourhood(&cands[e], &base_opts));
        }
        Rng(cfg.race_seed ^ (0xB000 + round_index as u64)).shuffle(&mut moves);
        moves.truncate(cfg.max_candidates);
        let mut added = 0usize;
        for (name, delta) in moves {
            if add(&mut cands, name, "elite_mutation", round_index + 1, delta)
                .map_err(|e| (4, e))?
            {
                added += 1;
            }
        }
        rounds.last_mut().unwrap().elites =
            std::iter::once(0).chain(elite.iter().copied()).collect();
        rounds.last_mut().unwrap().new_candidates = added;
        if added == 0 {
            stop_reason = "neighbourhood_exhausted";
            break;
        }
        stop_reason = "next_round";
    }
    debug_assert!(!rounds.is_empty());
    let last = rounds.last().expect("at least one round ran");
    let (done_from, done_to) = (last.from, last.to);
    let blocks_done: usize = rounds.iter().map(|r| r.to - r.from).sum();
    let group = last.survivors.clone();
    let final_test_p = last.final_p;
    let alive: Vec<usize> = (0..cands.len())
        .filter(|&c| cands[c].state == State::Alive)
        .collect();
    // The challenger: among the tuned candidates the last round's data cannot tell apart from the
    // best, the one that changes the fewest options (ties: lower mean cost). It only stands
    // against the baseline if it is cheaper on the blocks of that round; the hold-out decides.
    let winner = group
        .iter()
        .copied()
        .filter(|&c| c != 0)
        .min_by(|&a, &b| {
            cands[a].changes().cmp(&cands[b].changes()).then(
                cands[a]
                    .sgm_window(done_from, done_to)
                    .total_cmp(&cands[b].sgm_window(done_from, done_to)),
            )
        })
        .filter(|&c| {
            cands[c].sgm_window(done_from, done_to) < cands[0].sgm_window(done_from, done_to)
        })
        .unwrap_or(0);

    // ---- hold-out
    let mut holdout_json = Value::Null;
    let mut confirmed = false;
    let mut hold_wrong = 0usize;
    let mut ratio_json = Value::Null;
    if winner != 0 && !hold.is_empty() {
        let tasks: Vec<(usize, u64, usize)> = (0..hold.len())
            .flat_map(|i| {
                seeds
                    .iter()
                    .flat_map(move |&s| [(i, s, 0usize), (i, s, winner)])
            })
            .collect();
        let results = par_map(&tasks, cfg.jobs, |&(i, s, c)| {
            run_once(&hold[i], &cands[c].effective, s, cfg)
        });
        evaluations += tasks.len();
        for (&(i, _s, c), r) in tasks.iter().zip(&results) {
            spent += r.work(metric);
            if c == 0
                && r.status == HighsModelStatus::kOptimal
                && !known.contains_key(&hold[i].name)
            {
                consensus
                    .entry(hold[i].name.clone())
                    .or_insert(r.objective.unwrap_or(f64::NAN));
            }
        }
        let mut per: Vec<(String, f64, f64)> = Vec::new();
        for (i, inst) in hold.iter().enumerate() {
            let (mut b, mut w, mut k) = (0.0, 0.0, 0.0);
            for (&(ti, s, c), r) in tasks.iter().zip(&results) {
                if ti != i {
                    continue;
                }
                if let Some(wr) = check_answer(&cands[c].name, inst, s, r, &known, &consensus) {
                    if c != 0 {
                        hold_wrong += 1;
                    }
                    wrong.push(wr);
                }
                let cost = r.cost(cfg.run_cap(), metric);
                if c == 0 {
                    b += (cost + SHIFT).ln();
                    k += 1.0;
                } else {
                    w += (cost + SHIFT).ln();
                }
            }
            per.push((
                inst.name.clone(),
                (b / k).exp() - SHIFT,
                (w / k).exp() - SHIFT,
            ));
        }
        let logs: Vec<f64> = per
            .iter()
            .map(|(_, b, w)| ((w + SHIFT) / (b + SHIFT)).ln())
            .collect();
        let (m, lo, hi) = race::bootstrap_mean_ci(&logs, 4000, cfg.race_seed ^ 0x5EED);
        confirmed = hi < 0.0 && hold_wrong == 0;
        ratio_json = json!({"geomean_ratio": m.exp(), "ci95": [lo.exp(), hi.exp()],
                            "resamples": 4000, "unit": "per-instance mean over seeds, shifted"});
        holdout_json = json!({
            "instances": per.iter().map(|(n, b, w)| json!({"name": n, "baseline": b, "tuned": w})).collect::<Vec<_>>(),
            "baseline_sgm": sgm(per.iter().map(|p| p.1)),
            "tuned_sgm": sgm(per.iter().map(|p| p.2)),
            "ratio": ratio_json,
            "wrong_answers": hold_wrong,
        });
    }
    let verdict = if winner == 0 {
        "baseline_kept"
    } else if hold.is_empty() {
        "unconfirmed_no_holdout"
    } else if hold_wrong > 0 {
        "rejected_wrong_answer_on_holdout"
    } else if confirmed {
        "tuned_beats_baseline_on_holdout"
    } else {
        "holdout_not_confirmed"
    };

    // ---- the profile with its evidence
    let mut profile = Profile::new(cfg.base);
    profile.name = Some(cfg.name.clone().unwrap_or_else(|| {
        let dir = std::path::Path::new(&cfg.input);
        format!(
            "tuned-{}",
            dir.file_stem().and_then(|s| s.to_str()).unwrap_or("family")
        )
    }));
    if confirmed {
        profile.options = cands[winner].delta.clone();
    }
    let names = rusolve_mip::classify::FEATURE_NAMES;
    let vectors: Result<Vec<Vec<f64>>, String> = train
        .iter()
        .map(|i| profile::feature_vector(&i.path))
        .collect();
    // What an agent proposing candidates reads (LLM-proposal mode): the static features of
    // each training instance and what every candidate cost on it.
    let mut instance_rows: Vec<Value> = Vec::new();
    if let Ok(v) = &vectors {
        for (ii, inst) in train.iter().enumerate() {
            let mut cost_by_cand = serde_json::Map::new();
            for c in cands.iter() {
                let xs: Vec<f64> = blocks
                    .iter()
                    .enumerate()
                    .filter(|(_, (i, _))| *i == ii)
                    .filter_map(|(b, _)| c.costs.get(b).copied().flatten())
                    .collect();
                let mean = if xs.is_empty() {
                    f64::NAN
                } else {
                    xs.iter().sum::<f64>() / xs.len() as f64
                };
                cost_by_cand.insert(c.name.clone(), nan_null(mean));
            }
            instance_rows.push(json!({
                "name": inst.name, "sha256": inst.sha256,
                "features": names.iter().zip(v[ii].iter())
                    .map(|(n, x)| (n.to_string(), nan_null(*x))).collect::<BTreeMap<_, _>>(),
                "mean_cost_by_candidate": cost_by_cand,
            }));
        }
    }
    match vectors {
        Ok(v) => {
            profile.fingerprint = Some(Fingerprint::from_vectors(
                &names,
                &v,
                rusolve_mip::classify::FEATURE_VERSION,
            ))
        }
        Err(e) => eprintln!("rusolve tune: warning: no fingerprint: {e}"),
    }
    let base_c = &cands[0];
    let win_c = &cands[winner];
    let inst_role =
        |i: &Inst, role: &str| json!({"name": i.name, "sha256": i.sha256, "role": role});
    let cand_rows: Vec<Value> = cands
        .iter()
        .enumerate()
        .map(|(ci, c)| {
            json!({
                "name": c.name, "origin": c.origin, "round": c.round,
                "options": c.delta.iter().map(|(k, v)| (k.clone(), value_json(v))).collect::<BTreeMap<_, _>>(),
                "changed_options": c.changes(),
                "blocks_evaluated": c.evaluated(),
                "shifted_geomean_work": nan_null(c.sgm()),
                "ratio_to_baseline_same_blocks": nan_null(c.ratio_to(&cands[0])),
                "state": c.state.json_name(),
                "in_noise_group": group.contains(&ci),
            })
        })
        .collect();
    // One entry per racing round: who raced, who survived its test, who became the elite.
    let rounds_json: Vec<Value> = rounds
        .iter()
        .map(|r| {
            json!({
                "round": r.index,
                "candidates_raced": r.raced.len(),
                "candidates": r.raced.iter().map(|&c| cands[c].name.clone()).collect::<Vec<_>>(),
                "survivors": r.survivors.iter().map(|&c| cands[c].name.clone()).collect::<Vec<_>>(),
                "survivor_mean_ranks": r
                    .mean_rank
                    .iter()
                    .map(|v| nan_null(*v))
                    .collect::<Vec<_>>(),
                "elites": r.elites.iter().map(|&c| cands[c].name.clone()).collect::<Vec<_>>(),
                "next_round_candidates": r.new_candidates,
                "blocks": [r.from, r.to],
                "blocks_run": r.to - r.from,
                "final_test_p": r.final_p,
                "stopped": r.stopped,
            })
        })
        .collect();
    let wrong_json: Vec<Value> = wrong
        .iter()
        .map(|w| json!({"candidate": w.candidate, "instance": w.instance, "seed": w.seed,
            "status": w.status, "objective": w.objective, "expected": w.expected, "source": w.source}))
        .collect();
    let train_ratio = {
        let (b, w) = (
            base_c.sgm_window(done_from, done_to),
            win_c.sgm_window(done_from, done_to),
        );
        (w + SHIFT) / (b + SHIFT)
    };
    let race_json = json!({
        "test": "friedman_conover", "confidence": 1.0 - ALPHA, "first_test": FIRST_TEST, "each_test": 1,
        "block": "(instance, seed)", "blocks_planned": nb, "blocks_run": blocks_done,
        "tests_run": tests, "eliminations": drop_log,
        "final_test_p": final_test_p,
        "work_spent": spent, "budget_work": cfg.budget(),
        "candidates_raced": cands.len(), "survivors": alive.len(), "noise_group": group.len(),
        "iterations": rounds.len(), "iterations_asked": cfg.iterations, "elite": cfg.elite,
        "stopped": stop_reason, "rounds": rounds_json,
    });
    profile.provenance = Some(json!({
        "created": iso_now(),
        "tool": format!("rusolve tune {}", env!("CARGO_PKG_VERSION")),
        "metric": {"kind": metric.name(), "definition": metric.definition(),
                   "aggregate": "shifted_geomean", "shift": SHIFT, "run_cap": cfg.run_cap(),
                   "unfinished_penalty": "2 x cap"},
        "instances": train.iter().map(|i| inst_role(i, "train")).chain(hold.iter().map(|i| inst_role(i, "holdout"))).collect::<Vec<_>>(),
        "seeds": seeds, "evaluations": evaluations,
        "baseline": {"train": nan_null(base_c.sgm_window(done_from, done_to)), "holdout": holdout_json.get("baseline_sgm").cloned().unwrap_or(Value::Null)},
        "tuned": {"train": nan_null(win_c.sgm_window(done_from, done_to)), "holdout": holdout_json.get("tuned_sgm").cloned().unwrap_or(Value::Null)},
        "race_seed": cfg.race_seed, "base": cfg.base.name(),
        "iterations": rounds.len(), "cost_metric": metric.name(),
    }));
    let mut consensus_names: Vec<&String> = consensus.keys().collect();
    consensus_names.sort();
    profile.evidence = Some(json!({
        "verdict": verdict,
        "winner_candidate": win_c.name,
        "winner_round": win_c.round,
        "winner_options": win_c.delta.iter().map(|(k, v)| (k.clone(), value_json(v))).collect::<BTreeMap<_, _>>(),
        "train_ratio_tuned_over_baseline": nan_null(train_ratio),
        "holdout": holdout_json,
        "race": race_json,
        "candidates": cand_rows,
        "train_instances": instance_rows,
        "wrong_answers": wrong_json,
        "known_optima": {"file": cfg.known, "instances_known": known.len(),
                         "baseline_consensus_used_for": consensus_names},
    }));

    let profile_value = profile.to_value();
    // the profile must pass its own validator before it is written
    Profile::from_value(&profile_value)
        .map_err(|e| (4, format!("internal: own profile invalid: {e:?}")))?;
    if let Some(path) = &cfg.out {
        std::fs::write(
            path,
            serde_json::to_string_pretty(&profile_value).unwrap() + "\n",
        )
        .map_err(|e| (4, format!("cannot write {path}: {e}")))?;
    }
    if cfg.json {
        let out = json!({
            "schema_version": TUNE_SCHEMA,
            "solver": {"name": crate::envelope::SOLVER_NAME, "version": env!("CARGO_PKG_VERSION")},
            "input": cfg.input, "base": cfg.base.name(),
            "split": {"train": train.len(), "holdout": hold.len(), "seeds": cfg.seeds},
            "iterations": rounds.len(), "cost_metric": metric.name(),
            "verdict": verdict,
            "wrong_answers": wrong.len(),
            "profile_path": cfg.out,
            "profile": profile_value,
        });
        println!("{out}");
    } else {
        print_report(
            cfg,
            train.len(),
            hold.len(),
            &cands,
            &group,
            winner,
            verdict,
            &ratio_json,
            train_ratio,
            blocks_done,
            nb,
            spent,
            wrong.len(),
            &profile,
            &rounds,
            metric,
            stop_reason,
        );
    }
    Ok(())
}

fn nan_null(x: f64) -> Value {
    if x.is_finite() { json!(x) } else { Value::Null }
}

impl Cand {
    /// Shifted geometric mean over the blocks `[from, to)` of the race.
    fn sgm_window(&self, from: usize, to: usize) -> f64 {
        sgm(self
            .costs
            .iter()
            .skip(from)
            .take(to - from)
            .flatten()
            .copied())
    }
}

/// The candidates that differ from `elite` by exactly one more option of [`SPACE`]: an option the
/// elite does not change yet, set to a value other than the base's (the same value filter the
/// round-1 one-exchange moves use). This is the irace sampling step restricted to the discrete
/// neighbourhood rusolve has: mutating an elite of r changes gives candidates of r+1 changes, so
/// round 2 searches pairs and round 3 triples.
fn elite_neighbourhood(
    elite: &Cand,
    base_opts: &HighsOptions,
) -> Vec<(String, BTreeMap<String, OptionValue>)> {
    let mut out = Vec::new();
    for (name, values) in SPACE {
        if elite.delta.contains_key(*name) {
            continue;
        }
        let cur = rusolve_ekk::lp_data::highs_option_records::get_option_value(base_opts, name);
        for v in *values {
            let Ok(ov) = parse_option_text(name, v) else {
                continue;
            };
            if cur.as_ref() == Some(&ov) {
                continue;
            }
            let mut delta = elite.delta.clone();
            delta.insert((*name).to_string(), ov);
            out.push((format!("{}+{name}={v}", elite.name), delta));
        }
    }
    out
}

/// blocks x candidates cost table for the Friedman test, over the blocks `[from, to)`. Every
/// listed candidate ran every one of those blocks.
fn table_of(cands: &[Cand], cols: &[usize], from: usize, to: usize) -> Vec<Vec<f64>> {
    (from..to)
        .map(|b| {
            cols.iter()
                .map(|&c| cands[c].costs[b].expect("a live candidate ran every block"))
                .collect()
        })
        .collect()
}

/// The wrong-answer check against the known optimum, or the baseline's consensus.
fn check_answer(
    cand: &str,
    inst: &Inst,
    seed: u64,
    r: &Run,
    known: &HashMap<String, f64>,
    consensus: &HashMap<String, f64>,
) -> Option<Wrong> {
    let (expected, source) = match known.get(&inst.name) {
        Some(v) => (*v, "known_optimum"),
        None => (*consensus.get(&inst.name)?, "baseline_consensus"),
    };
    let bad = match r.status {
        HighsModelStatus::kOptimal => r
            .objective
            .is_none_or(|o| (o - expected).abs() > WRONG_TOL * expected.abs().max(1.0)),
        HighsModelStatus::kInfeasible | HighsModelStatus::kUnbounded => true,
        _ => false,
    };
    bad.then(|| Wrong {
        candidate: cand.to_string(),
        instance: inst.name.clone(),
        seed,
        status: util_model_status_to_string(r.status).to_string(),
        objective: r.objective,
        expected,
        source,
    })
}

#[allow(clippy::too_many_arguments)]
fn print_report(
    cfg: &Cfg,
    n_train: usize,
    n_hold: usize,
    cands: &[Cand],
    group: &[usize],
    winner: usize,
    verdict: &str,
    ratio: &Value,
    train_ratio: f64,
    blocks_done: usize,
    blocks: usize,
    spent: f64,
    wrong: usize,
    profile: &Profile,
    rounds: &[Round],
    metric: Metric,
    stop_reason: &str,
) {
    println!(
        "rusolve tune: {} instances ({n_train} train, {n_hold} held out), {} seeds, base {}",
        n_train + n_hold,
        cfg.seeds,
        cfg.base.name()
    );
    println!(
        "race: {} candidates over {} round(s) (elite {}), {blocks_done} of {blocks} blocks run, \
         work spent {spent:.3e} of {:.3e}, metric {}",
        cands.len(),
        rounds.len(),
        cfg.elite,
        cfg.budget(),
        metric.name()
    );
    for r in rounds {
        let elite: Vec<&str> = r
            .elites
            .iter()
            .filter(|&&c| c != 0)
            .map(|&c| cands[c].name.as_str())
            .collect();
        println!(
            "  round {}: {} candidates, {} blocks, {} survived{}{} ({})",
            r.index,
            r.raced.len(),
            r.to - r.from,
            r.survivors.len(),
            if r.new_candidates > 0 {
                format!(", {} new for round {}", r.new_candidates, r.index + 1)
            } else {
                String::new()
            },
            if elite.is_empty() {
                String::new()
            } else {
                format!("; elite: {}", elite.join(", "))
            },
            if r.stopped == "round_complete" {
                "round complete"
            } else {
                r.stopped
            },
        );
    }
    if stop_reason != "iterations_done" {
        println!("  race stopped: {stop_reason}");
    }
    println!(
        "  {:<46} {:>3} {:>4} {:>6} {:>12}  state",
        "candidate", "chg", "rnd", "blocks", "vs base*"
    );
    for (i, c) in cands.iter().enumerate() {
        println!(
            "  {:<46} {:>3} {:>4} {:>6} {:>12.3}  {}{}{}",
            c.name,
            c.changes(),
            c.round,
            c.evaluated(),
            c.ratio_to(&cands[0]),
            c.state.text_name(),
            if group.contains(&i) {
                " [noise group]"
            } else {
                ""
            },
            if i == winner { " <- winner" } else { "" },
        );
    }
    println!("  (* shifted geomean work relative to the baseline on the blocks the candidate ran)");
    println!("train: tuned/baseline shifted geomean = {train_ratio:.3}");
    if !ratio.is_null() {
        println!(
            "holdout: ratio {:.3}, 95% CI [{:.3}, {:.3}]",
            ratio["geomean_ratio"].as_f64().unwrap_or(f64::NAN),
            ratio["ci95"][0].as_f64().unwrap_or(f64::NAN),
            ratio["ci95"][1].as_f64().unwrap_or(f64::NAN)
        );
    }
    println!("wrong answers: {wrong}");
    println!("verdict: {verdict}");
    println!("winner first appeared in round {}", cands[winner].round);
    if profile.options.is_empty() {
        println!("profile options: (none changed)");
    } else {
        for (k, v) in &profile.options {
            println!("profile option: {k} = {}", v.to_token());
        }
    }
    match &cfg.out {
        Some(p) => println!("profile written to {p}"),
        None => println!("(no -o given: profile not written)"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_ekk::lp_data::highs_option_records;

    fn a(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn tune_space_is_valid() {
        for (name, values) in SPACE {
            for v in *values {
                parse_option_text(name, v).unwrap_or_else(|e| panic!("{e}"));
            }
        }
    }

    #[test]
    fn iterations_and_elite_default_and_zero_is_refused() {
        let c = parse_args(&a(&["fam"])).unwrap();
        assert_eq!((c.iterations, c.elite, c.jobs), (1, 3, 1));
        let c = parse_args(&a(&["fam", "--iterations", "5", "--elite", "7"])).unwrap();
        assert_eq!((c.iterations, c.elite), (5, 7));
        for bad in [
            a(&["--iterations", "0"]),
            a(&["--elite", "0"]),
            a(&["--iterations", "-2"]),
            a(&["--elite", "many"]),
        ] {
            let mut v = a(&["fam"]);
            v.extend(bad);
            assert!(parse_args(&v).is_err(), "{v:?} must be a usage error");
        }
    }

    #[test]
    fn the_c1_metric_needs_serial_solves_and_scales_the_budget_defaults() {
        assert_eq!(Metric::for_jobs(1), Metric::C1Work);
        assert_eq!(Metric::for_jobs(2), Metric::ItersNodes);
        let c1 = parse_args(&a(&["fam"])).unwrap();
        let par = parse_args(&a(&["fam", "--jobs", "4"])).unwrap();
        assert_eq!(
            (c1.metric().name(), par.metric().name()),
            ("c1_work", "iters_nodes")
        );
        // the defaults of the two budget flags are in the metric's unit...
        assert_eq!((c1.budget(), c1.run_cap()), (2e12, 4e11));
        assert_eq!((par.budget(), par.run_cap()), (1e7, 2e6));
        // ...but the node limit that enforces the cap inside the solve is the same either way
        assert_eq!(c1.node_cap(), par.node_cap());
        assert_eq!(par.node_cap(), 200_000);
        // an explicit cap is read in the metric's unit: 8e11 ticks is 4e6 work units, so 4e5 nodes
        let e = parse_args(&a(&["fam", "--run-work-cap", "8e11"])).unwrap();
        assert_eq!(e.node_cap(), 400_000);
        // and the old 1000-unit minimum still leaves 100 nodes of search
        let t = parse_args(&a(&["fam", "--run-work-cap", "1000"])).unwrap();
        assert_eq!(t.node_cap(), 100);
    }

    #[test]
    fn mutating_an_elite_changes_exactly_one_more_option() {
        let base_opts = Base::Recommended.options();
        let mut delta = BTreeMap::new();
        delta.insert(
            "mip_heuristic_run_rens".to_string(),
            OptionValue::Bool(false),
        );
        let elite = make_cand(
            "elite".to_string(),
            "one_exchange",
            1,
            delta,
            Base::Recommended,
            1,
        )
        .unwrap();
        let moves = elite_neighbourhood(&elite, &base_opts);
        assert!(!moves.is_empty(), "there are 15 options left to change");
        for (name, d) in &moves {
            assert_eq!(d.len(), 2, "{name} adds exactly one option to the elite");
            assert_eq!(d["mip_heuristic_run_rens"], OptionValue::Bool(false));
            let (k, v) = d
                .iter()
                .find(|(k, _)| *k != "mip_heuristic_run_rens")
                .unwrap();
            assert!(
                SPACE.iter().any(|(n, vals)| *n == k
                    && vals.iter().any(|t| parse_option_text(n, t).unwrap() == *v)),
                "{k} is not a SPACE option at a SPACE value"
            );
            let base = highs_option_records::get_option_value(&base_opts, k).expect("in registry");
            assert_ne!(value_json(v), value_json(&base), "{name} changes nothing");
        }
        // a candidate that already changes every option of SPACE has no neighbourhood left
        let all: BTreeMap<String, OptionValue> = SPACE
            .iter()
            .map(|(n, vals)| (n.to_string(), parse_option_text(n, vals[0]).unwrap()))
            .collect();
        let full = Cand {
            name: "full".to_string(),
            origin: "test",
            round: 1,
            delta: all,
            effective: base_opts.clone(),
            state: State::Alive,
            costs: vec![None; 1],
        };
        assert!(elite_neighbourhood(&full, &base_opts).is_empty());
    }

    #[test]
    fn known_file_formats() {
        let m = parse_known("=opt= foo 3.5\n# c\nbar.mps.gz,7\nbaz 1e3\n=best= q 9\n");
        assert_eq!(m["foo"], 3.5);
        assert_eq!(m["bar"], 7.0);
        assert_eq!(m["baz"], 1000.0);
        assert!(!m.contains_key("q"));
    }

    #[test]
    fn civil_date_is_correct() {
        assert_eq!(iso_now().len(), 20);
    }
}
