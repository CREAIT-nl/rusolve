//! LP-based branch-and-bound (of).
//!
//! Node selection: plunging into the preferred child with HiGHS's backtrack
//! rules, then hybrid best-estimate / best-bound picks from the queue
//! ([`crate::nodequeue`]). Branching: reliability branching with strong
//! branching on unreliable candidates, warm-started from the node basis
//! ([`crate::branching`]) —,.
//!
//! Pruning uses an *upper limit* below the incumbent value: with an integral
//! objective the next better solution is at least one unit lower, so nodes
//! whose bound exceeds `⌊z − ½⌋ + feastol` cannot improve it (HiGHS
//! `computeNewUpperLimit`, HighsMipSolverData.cpp:628, MIT; reimplemented for
//! scale 1 only — HiGHS also detects rational scales).
// TODO(perf): node copies carry full bound vectors; a domain change stack with
// backtracking (HighsDomain changeBound/backtrack :2023/:2294) replaces them in
// Together with propagation.

use crate::alns::{self, Alns};
use crate::branching::{BranchConfig, Brancher, SbHost, Selection};
use crate::conflict::{self, ConflictPool, ConflictStats, Explain, Revisit};
use crate::cuts::{self, CliqueFixStats, CutConfig, CutStats};
use crate::engine::ekk::EkkEngine;
use crate::engine::{
    CutRow, LpEngine, LpError, LpOutcome, LpStats, LpStatus, RusolveEngine, WarmStart,
};
use crate::engine_select::{MipEngine, ekk_requested};
use crate::heuristics::diving::DiveRule;
use crate::heuristics::rounding::Ctx;
use crate::heuristics::{self, HeuristicConfig, Locks, Source};
use crate::nodequeue::{self, BranchSide, NodeSelector, OpenNode};
use crate::problem::MipProblem;
use crate::propagate::Propagator;
use crate::race::Link;
use crate::redcost::{self, RedcostStats, RootRedcost};
use crate::restart::{self, Carry, RestartRequest};
use crate::submip::{self, FixingRateStats, SubMipConfig, SubMipStats};
use rusolve_lp::SimplexConfig;
use rusolve_lp::rng::XorShift64;
use std::rc::Rc;
use std::time::{Duration, Instant};

/// Absolute row/bound violation above which an LP point is rejected as an
/// incumbent (the gate's independent checker uses 1e-6 scaled).
const CANDIDATE_TOL: f64 = 1e-6;

/// Randomized-rounding trials on the first root LP point while there is no
/// incumbent (`Search::root_incumbent`). On jobshop_bigM s1–s3 one trial
/// fails (the rounded orders form a cycle, so the start-time LP is
/// infeasible) and 20 trials succeed 7–10 times, so 10 trials all miss with
/// probability about 1% (0.65^10).
const ROOT_RR_TRIALS: usize = 10;

#[derive(Debug, Clone, Copy)]
pub struct MipConfig {
    /// `|x − round(x)| ≤ integrality_tol` counts as integral.
    pub integrality_tol: f64,
    /// Margin in the cutoff (HiGHS uses its primal feasibility tolerance).
    pub feastol: f64,
    pub node_limit: Option<u64>,
    /// Wall-clock seconds.
    pub time_limit: Option<f64>,
    pub lp: SimplexConfig,
    /// Run MIP presolve (crate::presolve) once before the tree.
    pub presolve: bool,
    pub heuristics: HeuristicConfig,
    /// Root cutting planes.
    pub cuts: CutConfig,
    /// Re-solve node and cut-round LPs by warm dual simplex from the parent's
    /// (previous round's) basis; off = every LP cold, the A/B lever.
    pub warm_start: bool,
    /// Branching variable and node selection.
    pub branching: BranchConfig,
    /// Reduced-cost fixing at nodes and from the root duals; off = the A/B
    /// lever.
    pub redcost: bool,
    /// RENS at the root after the cut loop (unconditional, HiGHS-faithful)
    /// and RINS (RENS while there is no incumbent) at queue nodes behind
    /// the wall-time budget `submip::SUBMIP_TIME_SHARE`, as sub-MIPs
    /// ([`crate::submip`]); `--submip` opts in — default off: no trigger
    /// shape won on our benchmark set, `--no-submip` is the A/B lever.
    pub submip: bool,
    /// Build the root incumbent chain's fixing point by an
    /// LP-guided dive (coefficient rule, warm node LPs) instead of the raw
    /// post-cut LP point. HiGHS reaches its RENS/RINS neighbourhood by
    /// exactly such a dive (HighsPrimalHeuristics.cpp:394+ — `HighsSearch`
    /// with one-level backtracking; our static fixing was the documented
    /// deviation, submip.rs "Deviations", and its measured failure mode:
    /// the raw point fixes too little, so the restriction stays almost as
    /// big as the problem — air05 16.1 s for nothing, int_production's
    /// 25 ms point steering trees wrong). The dive's integral point also
    /// seeds the incumbent ladder directly (its intermediate simple
    /// roundings run as in the tree), and its incumbents arm — never
    /// disarm — the restart chain (`restart_arm_limit`). Default
    /// off: default-on was measured and blocked by a defect
    /// (mik's restart round 2 spins in a phase-1 LP cycle past the time
    /// limit — mik Optimal → Timeout on the default path);
    /// `--root-dive` opts in (implying `--submip`, the full chain;
    /// `--root-dive --no-submip` = dive only), `--no-root-dive` is a no-op
    /// lever for A/B scripts.
    pub root_dive: bool,
    /// Adaptive large-neighbourhood search: a bandit over RENS, RINS,
    /// mutation and crossover sub-MIPs at queue nodes, under a node and
    /// wall-time budget ([`crate::alns`]). It replaces the in-tree
    /// RINS/RENS of `submip`; off = the A/B lever.
    pub alns: bool,
    /// Local-MIP improver ([`heuristics::localmip`], Lin/Zou/Cai 2024):
    /// polishes every new incumbent within a FJ-relative effort budget
    /// ([`heuristics::localmip::POLISH_EFFORT_MULT`] per call,
    /// [`heuristics::localmip::TOTAL_EFFORT_MULT`] per search), run from
    /// `new_incumbent` so it complements the finders (FJ, rounding, diving,
    /// RINS) instead of competing with them. `--localmip` opts in,
    /// `--no-localmip` is the A/B lever.
    pub localmip: bool,
    /// Conflict analysis: the Farkas proof of an infeasible node LP becomes a
    /// stored conflict that prunes and tightens later nodes
    /// ([`crate::conflict`]); `--no-conflict` is the A/B lever.
    pub conflict: bool,
    /// Root restart: re-run MIP presolve on the root-tightened problem with
    /// the active cuts appended when enough integers got fixed, and redo the
    /// root ([`crate::restart`]);
    /// `--no-restart` is the A/B lever. Needs presolve on; off while racing.
    pub restart: bool,
    /// Root clique-table fixings: cliques from the model rows plus the
    /// objective cutoff row and the root cut pool, propagated into the global
    /// box on every improving incumbent (HiGHS `HighsCliqueTable`, cited in
    /// `cuts::clique`). Measured inert: without probing (HiGHS's
    /// actual fixings seed) no clique member ever reaches its literal
    /// value in the global box, so the cascade never fires — 0
    /// fixings on every probe instance, pg byte-identical base vs cand.
    /// Default off; `--clique-fixings` is the A/B lever.
    pub clique_fixings: bool,
    /// Root probing ([`crate::probing`]): fix free binaries one at a
    /// time, propagate, and apply the implied bound tightenings to the
    /// global box before the root LP (after HiGHS
    /// `HighsImplications::computeImplications` / `HPresolve::runProbing`,
    /// MIT). The seed the clique table and the restart fixings-%
    /// check were measured to lack: probing's fixings put clique members
    /// at their literals (arm the `--clique-fixings` cascade) and count as
    /// inactive integers in `restart::inactive_integers`. Default off;
    /// `--root-probing` is the A/B lever.
    pub root_probing: bool,
    /// Improvement-stalled restart: nodes the tree must process without an
    /// incumbent improvement before the restart may fire
    /// ([`restart::STALL_NODES`]; the config field exists so tests and A/Bs
    /// can move the bar without recompiling).
    pub restart_stall_nodes: u64,
    /// Improvement-stalled restart: seconds the stall window must also span
    /// ([`restart::STALL_SECS`]).
    pub restart_stall_secs: f64,
    /// Root-gap bar of every restart trigger ([`restart::MIN_ROOT_GAP`],
    /// our own design, no HiGHS counterpart). Lever
    /// (`--restart-min-gap`): the bar measures the gap against the
    /// incumbent, so a good incumbent (or a `--cutoff`) disarms the very
    /// reduced-cost-fixing restart it feeds.
    pub restart_min_gap: f64,
    /// Also learn conflicts from nodes cut off by the incumbent (the
    /// dual-bound proof). Default OFF, measured:
    /// with watch lists it wins on gap (−14% time, −25% nodes) but costs
    /// multi_knapsack (+30% time even with the conflict branching score) and
    /// int_production (+13% time; +the score, catastrophic tree growth).
    /// `--conflict-cutoff` is the opt-in.
    pub conflict_cutoff: bool,
    /// Feed the conflict pool into the branching score ([`crate::pseudocost`],
    /// HiGHS HighsPseudocost.h:310–417 with the decaying weight of
    /// HighsDomain.cpp:3833–3841). Default OFF, measured: it wins
    /// lights_out (−19% time), gap and multi_knapsack (with
    /// `--conflict-cutoff`), but hijacks branching on degenerate-LP
    /// instances where pseudocosts are flat — enlight_hard 2.7× nodes to
    /// timeout, market_split +13%.
    /// `--conflict-score` is the opt-in.
    pub conflict_score: bool,
    /// Dual steepest edge on warm node solves, weights carried in the basis
    /// snapshot. `false` = Dantzig warm pricing.
    pub warm_dse: bool,
    /// Lever (`--warm-unfit-last`): a node snapshot that no longer fits
    /// the LP rows warm-starts from the engine's latest basis instead of a
    /// cold primal ([`RusolveEngine::set_unfit_last`]). Default off: inert
    /// without `--tree-cuts`, and the 59-instance A/B of `--tree-cuts
    /// --warm-unfit-last --tree-cut-effort 1` wins cutoff runs but loses
    /// default-run time.
    pub warm_unfit_last: bool,
    /// Concurrent tree searches racing on the presolved problem
    /// ([`crate::race`]); 1 = a single search.
    pub threads: usize,
    /// Seed of the search's random choices (heuristics, Feasibility Jump);
    /// racers get distinct seeds.
    pub seed: u64,
    /// Audit lever (`--cutoff V`, default `None`): a caller-sense
    /// objective cutoff for the problem this config is applied to. The
    /// search starts with `upper_limit` at it and never lets an incumbent
    /// raise the limit above it — HiGHS's `objective_bound`
    /// (HighsMipSolverData.cpp:799–800, MIT), used to give the search a
    /// known optimum and isolate the dual side.
    /// Sub-MIPs never inherit it (their objective space differs).
    pub objective_cutoff: Option<f64>,
    /// Audit lever (`RUSOLVE_AUDIT=1` in the CLI, default off): emit
    /// `audit ...` stderr lines (presolve sizes, root bounds and times,
    /// incumbents) and time the LP solves into `LpStats::lp_nanos`.
    /// Never inherited by sub-MIPs. No effect on the search itself.
    pub audit: bool,
    /// Lever: cold LPs at the root node
    /// (the root LP, a restart round's root, cold fallbacks during the
    /// root cut loop) run on the dual simplex; every other LP keeps
    /// `lp`. HiGHS solves MIP LPs by dual simplex (HighsLpRelaxation.cpp,
    /// default `simplex_strategy`, MIT). `--dual-lp` is the all-LP lever.
    pub root_lp_dual: bool,
}

impl Default for MipConfig {
    fn default() -> Self {
        Self {
            integrality_tol: 1e-6,
            feastol: 1e-6,
            node_limit: None,
            time_limit: None,
            // Node LPs are never presolved: HiGHS switches LP presolve off
            // for the whole tree (HighsMipSolverData.cpp:1004) — per-node
            // reductions cost more than they save on small bound changes, and
            // LP-presolve postsolve bugs became wrong MIP answers here
            // (set_partition[infeasible] s2).
            // Cold root/node LPs stay primal: a cold dual default lost on
            // our benchmark set; `--dual-lp` opts in.
            lp: SimplexConfig {
                presolve: false,
                ..SimplexConfig::default()
            },
            presolve: true,
            heuristics: HeuristicConfig::default(),
            cuts: CutConfig::default(),
            warm_start: true,
            branching: BranchConfig::default(),
            redcost: true,
            // Measured, stays off. The wall-time budget (submip.rs)
            // bounds the old 1.14x/1.20x time loss, but no shape of the
            // trigger wins: with the root call unconditional the
            // static-fixing first incumbent steers whole trees wrong
            // (int_production s1 Optimal->TimeLimit) and one root RENS on a
            // big LP costs 16 s (air05); budget-gating the root keeps cheap
            // parents from ever collecting the wins (gap/gen_int_knapsack
            // 0.4-0.5x time with it).
            submip: false,
            // Prototype lever, default off: measured back to default-off.
            // Dive-only was the default-on
            // candidate (the polish never paid: air05 1.40x / neos-860300
            // 1.36x / sp150x300d 1.58x are sub-MIP setup costs, and mik's
            // winning incumbent -52 266 is dive-delivered, not polish), but
            // the dive's better point arms the restart and
            // mik's restart round 2 hits a phase-1 LP spin past the
            // time limit. Optimal → Timeout on
            // the default path disqualifies default-on until that is fixed;
            // `--root-dive` stays the opt-in lever.
            root_dive: false,
            alns: false,
            localmip: false,
            conflict: true,
            // Default on; `--no-restart` is the lever.
            restart: true,
            // Measured inert (0 fixings on every probe
            // instance; the seed fixings live in HiGHS's probing)
            // — default off; `--clique-fixings` is the lever.
            clique_fixings: false,
            // Default off; `--root-probing` is the lever.
            root_probing: false,
            restart_stall_nodes: restart::STALL_NODES,
            restart_stall_secs: restart::STALL_SECS,
            restart_min_gap: restart::MIN_ROOT_GAP,
            conflict_cutoff: false,
            conflict_score: false,
            warm_dse: false,
            warm_unfit_last: false,
            threads: 1,
            seed: 0,
            objective_cutoff: None,
            audit: false,
            root_lp_dual: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MipStatus {
    Optimal,
    Infeasible,
    Unbounded,
    NodeLimit,
    TimeLimit,
}

#[derive(Debug, Clone)]
pub struct MipSolution {
    pub status: MipStatus,
    /// Caller-sense incumbent value; NaN without an incumbent, ∓∞ when
    /// Unbounded.
    pub objective: f64,
    /// Caller-sense dual bound.
    pub bound: f64,
    /// Best integer-feasible point found (empty when none).
    pub x: Vec<f64>,
    pub nodes: u64,
    pub lp_iterations: u64,
    /// Root cut loop statistics.
    pub cuts: CutStats,
    /// When and where the first incumbent appeared (None without one).
    pub first_incumbent: Option<FirstIncumbent>,
    /// Every accepted incumbent in improving order, first included: the
    /// incumbent ladder triage needed and had to reconstruct from
    /// short-limit reruns.
    pub incumbents: Vec<FirstIncumbent>,
    /// Warm/cold LP solve counters.
    pub lp_stats: LpStats,
    /// Reduced-cost fixing counters.
    pub redcost: RedcostStats,
    /// RINS / RENS sub-MIP counters.
    pub submip: SubMipStats,
    /// Conflict analysis counters.
    pub conflict: ConflictStats,
    /// Root clique-table fixings counters.
    pub cliques: CliqueFixStats,
    /// Root-probing counters.
    pub probing: crate::probing::ProbingStats,
    /// Root restarts performed before this answer.
    pub restarts: u32,
    /// Set by `bnb` when the root asks for a restart; consumed by
    /// `presolve::solve_presolved`'s restart loop and never observed by
    /// external callers (every other entry point forces `restart` off).
    pub(crate) restart: Option<Box<RestartRequest>>,
}

#[derive(Debug, Clone, Copy)]
pub struct FirstIncumbent {
    /// The incumbent's objective (raw `p.objective`, as the tree sees it).
    pub value: f64,
    pub seconds: f64,
    /// Nodes processed before it was found (0 = before the root LP).
    pub nodes: u64,
    pub source: Source,
}

#[derive(Debug, Clone, PartialEq)]
pub enum MipError {
    Lp(LpError),
    /// A node LP below a bounded root came back unbounded: numerical trouble,
    /// since child regions are subsets of the root region.
    UnboundedNode,
    /// A node LP reported Optimal with a point violating the node's rows or
    /// bounds by `violation`; the node bounds make it reproducible.
    InvalidLpSolution {
        violation: f64,
        lower: Vec<f64>,
        upper: Vec<f64>,
    },
}

impl std::fmt::Display for MipError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MipError::Lp(e) => write!(f, "node LP failed: {e}"),
            MipError::UnboundedNode => write!(f, "node LP unbounded below a bounded root"),
            MipError::InvalidLpSolution { violation, .. } => {
                write!(
                    f,
                    "node LP returned an infeasible 'optimal' point (violation {violation:e})"
                )
            }
        }
    }
}

impl std::error::Error for MipError {}

impl From<LpError> for MipError {
    fn from(e: LpError) -> Self {
        MipError::Lp(e)
    }
}

/// Solves `p` to optimality (within `feastol` on the objective) or a limit.
pub fn solve(p: &MipProblem, config: &MipConfig) -> Result<MipSolution, MipError> {
    if config.presolve {
        return crate::presolve::solve_presolved(p, config);
    }
    if config.threads > 1 {
        return crate::race::solve(p, config);
    }
    // A restart needs the presolve chain to re-enter and to post the answer
    // back through every level: without presolve there is nothing to restart.
    solve_tree(
        p,
        &MipConfig {
            restart: false,
            ..*config
        },
        None,
        None,
    )
}

/// One racer's search (`config.presolve` is off: the race shares one
/// presolved problem, so a restart would re-presolve under one racer only —
/// disabled).
pub(crate) fn solve_racing(
    p: &MipProblem,
    config: &MipConfig,
    link: Link<'_>,
) -> Result<MipSolution, MipError> {
    solve_tree(
        p,
        &MipConfig {
            restart: false,
            ..*config
        },
        Some(link),
        None,
    )
}

pub(crate) fn solve_tree(
    p: &MipProblem,
    config: &MipConfig,
    race: Option<Link<'_>>,
    carry: Option<Carry>,
) -> Result<MipSolution, MipError> {
    let start = Instant::now();
    let Some(root) = root_node(p, config.integrality_tol) else {
        return Ok(finish(p, MipStatus::Infeasible, None, f64::INFINITY, 0, 0));
    };
    let mut search = Search::new(p, engine(p, config), config, start);
    search.race = race;
    // Only the top-level, non-racing search may ask for a restart.
    search.restart_armed = config.restart && search.race.is_none();
    // Root probing: tighten the global box before anything else touches
    // it (the root LP, the cut loop and the restart fixings-% check all see
    // the probed box). Skipped inside racers — one probing pass per solve,
    // like HiGHS's presolve-shared pass (the race shares one presolved
    // problem; probing each racer would multiply the cost).
    if config.root_probing && search.race.is_none() {
        // Pass 1 has no incumbent to cut against: cutoff ∞ (pass 2, at the
        // end of the root cut loop, probes with the incumbent's cutoff).
        search.root_probing(f64::INFINITY);
        if search.global_infeasible {
            return Ok(finish(p, MipStatus::Infeasible, None, f64::INFINITY, 0, 0));
        }
    }
    if let Some(carry) = carry {
        search.nodes = carry.nodes;
        search.lp_iterations = carry.lp_iterations;
        search.heur_lp_iterations = carry.heur_lp_iterations;
        search.restarts_done = carry.restarts_done;
        search.first_incumbent = carry.first_incumbent;
        search.incumbents = carry.incumbents;
        if let Some((z, x)) = carry.incumbent {
            search.install_incumbent(z, x);
        }
    }
    search.pre_root_heuristics(&root);
    let end = match search.run(root.clone()) {
        Ok(end) => end,
        // A solve that crossed the search deadline INSIDE the LP (the
        // engine's `set_deadline` bailout; the root, a restart round's root
        // or any other site the depth>0 branch-without-LP arm does not
        // cover) reports the budget end, not an error — the same verdict
        // `limit_hit` produces between nodes. Gated on the time limit
        // actually being hit: a root LP that reaches the pure iteration
        // cap with time left stays an honest engine error. No dual bound
        // is known at this boundary, so −∞ (an unfinished search, never a
        // false bound), except the root LP's own optimum, when it was
        // reached before the bailout (co-100 crossed the deadline inside a
        // cut-loop re-solve and reported −∞ over a solved root LP).
        Err(MipError::Lp(LpError::IterationLimit))
            if config
                .time_limit
                .is_some_and(|t| start.elapsed().as_secs_f64() >= t) =>
        {
            TreeEnd::Limit(MipStatus::TimeLimit, search.root_lp_bound)
        }
        Err(e) => return Err(e),
    };
    if config.audit {
        search.print_tree_profile();
    }
    let first_incumbent = search.first_incumbent;
    let ladder = std::mem::take(&mut search.incumbents);
    let mut sol = conclude(p, config, start, root, search, end)?;
    sol.first_incumbent = sol.first_incumbent.or(first_incumbent);
    if sol.incumbents.is_empty() {
        sol.incumbents = ladder;
    }
    Ok(sol)
}

/// Audit counters of the node LPs.
#[derive(Debug, Default)]
struct TreeProfile {
    solves: u64,
    nanos: u64,
    rows_sum: u64,
    rows_max: usize,
    cut_rows_sum: u64,
    /// Node LPs whose snapshot did not fit the rows (recovered from the
    /// latest basis or solved cold) and their iterations; node LPs of
    /// queue (not plunge) nodes and their iterations.
    unfit_lps: u64,
    unfit_iters: u64,
    queue_lps: u64,
    queue_iters: u64,
    /// The next node-LP count at which a `tree_rows` trajectory line prints.
    next_trace: u64,
}

impl TreeProfile {
    fn node_lp(&mut self, nanos: u64, rows: usize, (lp_cuts, pool): (usize, usize), t: f64) {
        self.solves += 1;
        self.nanos += nanos;
        self.rows_sum += rows as u64;
        self.rows_max = self.rows_max.max(rows);
        self.cut_rows_sum += lp_cuts as u64;
        if self.solves >= self.next_trace {
            eprintln!(
                "audit tree_rows lps={} rows={rows} lp_cuts={lp_cuts} pool={pool} t={t:.6}",
                self.solves
            );
            self.next_trace = (self.next_trace * 4).max(16);
        }
    }
}

/// The node-LP engine: [`RusolveEngine`] by default, the HiGHS-port
/// `EkkEngine` under `RUSOLVE_MIP_EKK=1` (its options come from the HiGHS MIP defaults, not `config.lp`).
fn engine(p: &MipProblem, config: &MipConfig) -> MipEngine {
    if ekk_requested() {
        return MipEngine::Ekk(Box::new(EkkEngine::new(p)));
    }
    MipEngine::Rusolve(rusolve_engine(p, config))
}

fn rusolve_engine(p: &MipProblem, config: &MipConfig) -> RusolveEngine {
    let mut e = RusolveEngine::new(p, config.lp);
    e.set_warm(config.warm_start);
    e.set_warm_dse(config.warm_dse);
    e.set_unfit_last(config.warm_unfit_last);
    e.set_time_lps(config.audit);
    e
}

/// Turns the tree's end state into a solution (runs the feasibility search
/// for an unbounded root without an incumbent).
fn conclude<E: LpEngine>(
    p: &MipProblem,
    config: &MipConfig,
    start: Instant,
    root: Node,
    mut search: Search<'_, E>,
    end: TreeEnd,
) -> Result<MipSolution, MipError> {
    let (nodes, iters) = (search.nodes, search.lp_iterations);
    let cut_stats = std::mem::take(&mut search.cut_stats);
    let mut lp_stats = search.engine.stats();
    lp_stats.no_lp_branches = search.no_lp_branches;
    lp_stats.unsolved_pruned = search.unsolved_pruned;
    let redcost_stats = search.redcost_stats;
    let submip_stats = search.submip_stats;
    let conflict_stats = search.conflict_stats;
    let clique_stats = search.clique_stats;
    let probing_stats = search.probing_stats;
    let mut sol = match end {
        TreeEnd::Exhausted => {
            let status = match search.incumbent {
                Some(_) => MipStatus::Optimal,
                None => MipStatus::Infeasible,
            };
            let bound = search.incumbent.as_ref().map_or(f64::INFINITY, |i| i.0);
            finish(p, status, search.incumbent, bound, nodes, iters)
        }
        TreeEnd::Limit(status, open_bound) => {
            let bound = search
                .incumbent
                .as_ref()
                .map_or(open_bound, |i| i.0.min(open_bound));
            finish(p, status, search.incumbent, bound, nodes, iters)
        }
        TreeEnd::Restart(req) => {
            // Consumed by `presolve::solve_presolved`'s restart loop (the
            // only caller that arms a search). The status is deliberately a
            // non-final limit: should the request ever escape, the answer
            // reads as an unfinished search, never as a false optimality.
            let bound = search
                .incumbent
                .as_ref()
                .map_or(req.bound, |i| i.0.min(req.bound));
            let mut sol = finish(
                p,
                MipStatus::TimeLimit,
                search.incumbent.clone(),
                bound,
                nodes,
                iters,
            );
            sol.restarts = req.carry.restarts_done;
            sol.restart = Some(req);
            sol
        }
        TreeEnd::RootUnbounded => {
            // Meyer (1974): with rational data, unbounded relaxation + any
            // feasible point ⇒ unbounded MIP. A heuristic point settles it;
            // otherwise ask the feasibility question by B&B.
            if let Some((_, x)) = search.incumbent {
                let mut sol = finish(
                    p,
                    MipStatus::Unbounded,
                    None,
                    f64::NEG_INFINITY,
                    nodes,
                    iters,
                );
                sol.objective = p.caller_objective(f64::NEG_INFINITY);
                sol.x = x;
                sol.cuts = cut_stats;
                sol.lp_stats = lp_stats;
                sol.redcost = redcost_stats;
                sol.submip = submip_stats;
                sol.conflict = conflict_stats;
                return Ok(sol);
            }
            let fp = p.feasibility_version();
            let mut feas = Search::new(&fp, engine(&fp, config), config, start);
            feas.race = search.race.as_ref().map(Link::stop_only);
            feas.nodes = nodes;
            feas.lp_iterations = iters;
            let end = feas.run(root)?;
            let (nodes, iters) = (feas.nodes, feas.lp_iterations);
            lp_stats.absorb(&feas.engine.stats());
            // Unbounded needs a verified feasible point, whatever stopped the
            // search; a limit without one stays a limit, never Unbounded.
            let status = match (end, &feas.incumbent) {
                (_, Some(_)) => MipStatus::Unbounded,
                (TreeEnd::Limit(status, _), None) => status,
                (TreeEnd::Exhausted, None) => MipStatus::Infeasible,
                // The feasibility search is never armed for a restart
                // (`restart_armed` is set only by `solve_tree`), so this arm
                // is unreachable; kept for exhaustiveness.
                (TreeEnd::Restart(_), None) => MipStatus::TimeLimit,
                // A zero-cost LP cannot be unbounded: the engine misreported,
                // and mapping it to Infeasible would be a silent wrong answer.
                (TreeEnd::RootUnbounded, None) => return Err(MipError::UnboundedNode),
            };
            let mut sol = finish(p, status, None, f64::NEG_INFINITY, nodes, iters);
            sol.first_incumbent = feas.first_incumbent;
            sol.incumbents = std::mem::take(&mut feas.incumbents);
            if status == MipStatus::Unbounded {
                sol.objective = p.caller_objective(f64::NEG_INFINITY);
                sol.x = feas.incumbent.map(|i| i.1).unwrap_or_default();
            } else if status == MipStatus::Infeasible {
                sol.bound = p.caller_objective(f64::INFINITY);
            }
            sol
        }
    };
    sol.cuts = cut_stats;
    sol.lp_stats = lp_stats;
    sol.redcost = redcost_stats;
    sol.submip = submip_stats;
    sol.conflict = conflict_stats;
    sol.cliques = clique_stats;
    sol.probing = probing_stats;
    Ok(sol)
}

fn finish(
    p: &MipProblem,
    status: MipStatus,
    incumbent: Option<(f64, Vec<f64>)>,
    bound: f64,
    nodes: u64,
    lp_iterations: u64,
) -> MipSolution {
    let (objective, x) = match incumbent {
        Some((z, x)) => (p.caller_objective(z), x),
        None => (f64::NAN, vec![]),
    };
    MipSolution {
        status,
        objective,
        bound: p.caller_objective(bound),
        x,
        nodes,
        lp_iterations,
        first_incumbent: None,
        incumbents: vec![],
        cuts: CutStats::default(),
        lp_stats: LpStats::default(),
        redcost: RedcostStats::default(),
        submip: SubMipStats::default(),
        conflict: ConflictStats::default(),
        cliques: CliqueFixStats::default(),
        probing: crate::probing::ProbingStats::default(),
        restarts: 0,
        restart: None,
    }
}

/// Root bounds with integer columns rounded inward; `None` if that empties a
/// domain.
/// Integer columns whose LP value is more than `tol` from integral, as
/// `(column, value)`. Each value is clamped into the node box first, as
/// HiGHS does when it collects fractional integers (HighsLpRelaxation.cpp:
/// 1467–1472, MIT). Without the clamp, an LP point a hair outside a bound
/// (x = 1.00009 with upper 1) would branch to a child identical to its
/// parent forever. Branching rules take their candidates from here.
fn fractional_integers<'a>(
    integer: &'a [bool],
    x: &'a [f64],
    lower: &'a [f64],
    upper: &'a [f64],
    tol: f64,
) -> impl Iterator<Item = (usize, f64)> + 'a {
    x.iter().enumerate().filter_map(move |(j, &v)| {
        if !integer[j] {
            return None;
        }
        let v = v.clamp(lower[j], upper[j]);
        let f = v - v.floor();
        (f.min(1.0 - f) > tol).then_some((j, v))
    })
}

fn root_node(p: &MipProblem, tol: f64) -> Option<Node> {
    let mut lower = p.col_lower.clone();
    let mut upper = p.col_upper.clone();
    for j in 0..p.num_col {
        if p.integer[j] {
            lower[j] = (lower[j] - tol).ceil();
            upper[j] = (upper[j] + tol).floor();
        }
        if lower[j] > upper[j] {
            return None;
        }
    }
    Some(Node {
        lower,
        upper,
        bound: f64::NEG_INFINITY,
        depth: 0,
        branched: None,
        inherited: Vec::new(),
        warm: None,
        estimate: f64::NEG_INFINITY,
        side: None,
    })
}

/// Appends the columns whose bound differs between the boxes `(lo, hi)` and
/// `(lower, upper)`. Node boxes only tighten, so a difference is a move.
fn moved_cols(lo: &[f64], hi: &[f64], lower: &[f64], upper: &[f64], out: &mut Vec<usize>) {
    for j in 0..lower.len() {
        if lo[j] != lower[j] || hi[j] != upper[j] {
            out.push(j);
        }
    }
}

#[derive(Clone)]
struct Node {
    lower: Vec<f64>,
    upper: Vec<f64>,
    /// Parent LP objective: a valid lower bound for this subtree.
    bound: f64,
    depth: u32,
    /// Column whose bound the branching changed (`None` at the root): the
    /// only propagation seed besides the objective row, since the parent's
    /// bounds were already propagated.
    branched: Option<usize>,
    /// Columns whose bound moved between the parent's box before its
    /// propagation and this child's creation (reduced-cost tightening,
    /// propagation, conflict propagation): with the branch and this node's
    /// own propagation moves, the columns a direct plunge child must revisit
    /// in the conflict pool.
    inherited: Vec<usize>,
    /// The parent's optimal basis, shared by both children (HiGHS keeps it
    /// per plunge node, `NodeData::nodeBasis`, HighsSearch.cpp:960). Queued
    /// nodes keep it too, unlike HiGHS' `HighsNodeQueue::OpenNode`
    /// (HighsNodeQueue.h:142-172): dropping it and starting a popped node from
    /// the engine's current basis measured 1.60x time sgm (mixed) with 29/105
    /// trees and two objectives changed and no peak-RSS movement, because a
    /// foreign basis costs us a cold two-phase re-solve where HiGHS repairs one
    /// (`recoverBasis`, HighsSearch.cpp:788-793).
    warm: Option<Rc<WarmStart>>,
    /// Parent's best-estimate of the subtree's integer optimum (node queue).
    estimate: f64,
    /// The branching that created this node and the parent LP objective,
    /// for pseudocost observations (`None` at the root).
    side: Option<(BranchSide, f64)>,
}

impl OpenNode for Node {
    fn lower_bound(&self) -> f64 {
        self.bound
    }
    fn estimate(&self) -> f64 {
        self.estimate
    }
    fn depth(&self) -> u32 {
        self.depth
    }
    fn side(&self) -> Option<BranchSide> {
        self.side.map(|s| s.0)
    }
}

impl Node {
    /// Child on the `up`/down side of `col` at LP value `value`, warm-started
    /// from `warm` (the parent's basis).
    #[allow(clippy::too_many_arguments)] // one child is its full box plus lineage
    fn child(
        &self,
        col: usize,
        value: f64,
        up: bool,
        bound: f64,
        parent_obj: f64,
        warm: Option<Rc<WarmStart>>,
        inherited: Vec<usize>,
    ) -> Node {
        let mut n = Node {
            lower: self.lower.clone(),
            upper: self.upper.clone(),
            bound,
            depth: self.depth + 1,
            branched: Some(col),
            inherited,
            warm,
            estimate: self.estimate,
            side: Some((
                BranchSide {
                    col,
                    point: value,
                    up,
                },
                parent_obj,
            )),
        };
        if up {
            n.lower[col] = value.ceil();
        } else {
            n.upper[col] = value.floor();
        }
        n
    }
}

/// What became of a tree node whose LP did not solve
/// ([`Search::unsolved_node`]).
enum Unsolved {
    /// Children to explore, first one first (boxed to keep the enum small).
    Branch(Box<(Node, Node)>),
    /// The fallback re-solve of an all-integers-fixed node finished.
    Solved(LpOutcome),
    /// The fallback failed too: the node is dropped.
    Pruned,
}

enum TreeEnd {
    Exhausted,
    RootUnbounded,
    /// Limit status and the best bound among unexplored nodes.
    Limit(MipStatus, f64),
    /// The root asks for a restart: the
    /// request travels out through `MipSolution::restart` to the presolve
    /// loop, which re-enters `solve_tree` on the reduced problem. Boxed to
    /// keep `TreeEnd` small.
    Restart(Box<RestartRequest>),
}

struct Search<'a, E: LpEngine> {
    p: &'a MipProblem,
    engine: E,
    propagator: Propagator,
    config: &'a MipConfig,
    start: Instant,
    /// Every improving solution has an integral objective value.
    integral_objective: bool,
    incumbent: Option<(f64, Vec<f64>)>,
    /// Nodes with a bound above this cannot improve the incumbent.
    upper_limit: f64,
    nodes: u64,
    lp_iterations: u64,
    first_incumbent: Option<FirstIncumbent>,
    /// Every accepted incumbent, in improving order (the ladder).
    incumbents: Vec<FirstIncumbent>,
    locks: Locks,
    rng: XorShift64,
    /// LP iterations spent inside heuristics (part of `lp_iterations`).
    heur_lp_iterations: u64,
    /// Root bounds after root propagation: the domain heuristics fix from.
    global_lower: Vec<f64>,
    global_upper: Vec<f64>,
    cut_stats: CutStats,
    /// Root cut-loop state, reused for node separation.
    separator: Option<cuts::Separator>,
    brancher: Brancher,
    num_integer: usize,
    /// First root LP iterations and solved non-root node LPs (count,
    /// iterations): they scale the tree LP iteration cap.
    root_lp_iterations: u64,
    /// The first root LP's optimal objective (−∞ until solved): the bound a
    /// Deadline bailout reports.
    root_lp_bound: f64,
    node_lp_solves: u64,
    node_lp_iterations: u64,
    /// Audit (`MipConfig::audit` only): node-LP wall time and the LP
    /// row counts (total, cut rows, max) summed over node LP solves.
    prof: TreeProfile,
    /// Basis of the node being branched: where strong-branching probes start.
    sb_start: Option<Rc<WarmStart>>,
    /// Root reduced costs, re-applied to the global box at every new upper
    /// limit (crate::redcost).
    root_redcost: RootRedcost,
    redcost_stats: RedcostStats,
    /// Root clique table (row cliques at build, objective cliques per
    /// incumbent): the incumbent → fixings link (`cuts::clique`).
    cliques: Option<cuts::CliqueTable>,
    clique_stats: cuts::CliqueFixStats,
    /// Root-probing counters (the box tightening runs in `solve_tree`
    /// before the root LP; the stats ride the search out for reporting).
    probing_stats: crate::probing::ProbingStats,
    submip_stats: SubMipStats,
    /// Wall-clock seconds spent in sub-MIPs ([`submip::SUBMIP_TIME_SHARE`]'s
    /// budget; part of the elapsed time, charged by `run_submip`).
    submip_secs: f64,
    /// Fixing-rate observations of this search's sub-MIPs and the stream
    /// that draws their targets and seeds.
    submip_rates: FixingRateStats,
    submip_rng: XorShift64,
    /// `--localmip` accounting: effort spent (charged by `run_localmip`
    /// against [`heuristics::localmip::TOTAL_EFFORT_MULT`]) and the
    /// call/improving counters for `RUSOLVE_LOCALMIP_LOG`.
    localmip_effort: usize,
    localmip_calls: u64,
    localmip_found: u64,
    /// ALNS state (`config.alns`), the node count when the incumbent last
    /// improved, and the root node's wall time: ALNS's budget inputs.
    alns: Option<Alns>,
    incumbent_node: u64,
    root_secs: f64,
    /// Conflicts learnt from infeasible node LPs (`config.conflict`).
    conflicts: ConflictPool,
    conflict_stats: ConflictStats,
    /// Root reduced-cost fixing (with global propagation) proved that no
    /// point beats the incumbent: the search is over.
    global_infeasible: bool,
    /// Incumbent exchange and stop flag when racing.
    race: Option<Link<'a>>,
    /// Root restarts already performed before this search round
    /// (`Carry::restarts_done`).
    restarts_done: u32,
    /// Tree nodes whose LP did not solve and were branched on without
    /// one, and all-integers-fixed nodes pruned after the fallback re-solve
    /// failed too (`LpStats::no_lp_branches` / `unsolved_pruned`).
    no_lp_branches: u64,
    unsolved_pruned: u64,
    /// Only the top-level search of a presolve round may ask for a restart:
    /// armed by `solve_tree`, never by the feasibility sub-search of an
    /// unbounded root.
    restart_armed: bool,
    /// The incumbent limit as it stood when the root cut loop ended,
    /// stashed before the root sub-MIP/dive block can publish a better point.
    /// The restart check's gap bar ([`restart::MIN_ROOT_GAP`], both restart kinds) is
    /// measured against *this*, not the current limit: the bar's population
    /// separation was calibrated on cut-loop-end incumbents (the small
    /// families check in at 0.013–0.047 with the rounding ladder's points,
    /// mik at 0.163), and a point the chain's own first link published
    /// seconds before the check must arm the restart it feeds, not disarm it
    /// (mik: dive incumbent −52 266 closes the check-time gap to 0.0395,
    /// blocking both restart kinds; the cut-loop limit −46 719 gives 0.1629 and the
    /// restart fires). The dive's better incumbent still rides the restart —
    /// `RestartRequest::incumbent` carries the current one. No HiGHS
    /// counterpart needed: HiGHS has no root-gap bar (its restarts ride
    /// presolve fixings, HighsMipSolverData.cpp:2486–2505), so no incumbent
    /// can block them.
    restart_arm_limit: f64,
    /// Stall restart : the root context for an
    /// improvement-stalled restart, stashed when the root check found
    /// the no-fixings shape (wide gap, cut mass) and consumed at a queue node
    /// once the tree has stalled. `None` = not armed.
    stall_restart: Option<StallRestart>,
}

/// What the root stashes for a stall restart: the cut log and final
/// root LP of the round that armed it, plus the stall window. The window
/// re-opens whenever the incumbent improves (or appears) — an improving tree
/// is not stalled — and the restart fires only at a queue node once
/// [`restart::STALL_NODES`] nodes and [`restart::STALL_SECS`] seconds have
/// passed without that happening and the gap over the tree's own dual bound
/// is still [`restart::MIN_ROOT_GAP`] wide.
struct StallRestart {
    /// The root cut loop's accepted cuts (see [`restart::active_cuts`]).
    log: Vec<(cuts::CutClass, CutRow)>,
    /// The final root LP point the cuts are filtered against.
    x: Vec<f64>,
    /// The final root LP objective (internal sense): the fallback dual bound
    /// of the restart request.
    lp_obj: f64,
    /// Incumbent value when the window opened (`None` = no incumbent).
    watch: Option<f64>,
    /// `Search::nodes` / search-elapsed seconds when the window opened.
    since_nodes: u64,
    since_secs: f64,
}

impl StallRestart {
    /// The stall verdict for the current state: fire, or keep waiting (after
    /// re-opening the window if the incumbent significantly improved since
    /// it was last seen — a crawl under [`restart::STALL_IMPROVEMENT`] does
    /// not re-open it). No gap leg in-tree, deliberately: the tree's own
    /// open bound converges while the incumbent crawls (pg's dual side sits
    /// 0.2 % from optimal while its primal takes 20 s), so a gap-over-tree
    /// condition would block the exact instance the stall restart exists for. The
    /// gap bar that matters was checked at arm time, over the root LP.
    fn verdict(
        &mut self,
        incumbent: Option<f64>,
        nodes: u64,
        elapsed_secs: f64,
        stall_nodes: u64,
        stall_secs: f64,
    ) -> StallRestartVerdict {
        let significant = match (self.watch, incumbent) {
            // First sighting (also covers a carried None → Some transition):
            // a new plan restarts the clock.
            (None, _) => incumbent.is_some(),
            // Incumbent vanished: keep the window as it was.
            (Some(_), None) => false,
            (Some(w), Some(z)) => (w - z) / w.abs().max(1.0) >= restart::STALL_IMPROVEMENT,
        };
        if significant {
            self.watch = incumbent;
            self.since_nodes = nodes;
            self.since_secs = elapsed_secs;
            return StallRestartVerdict::Wait;
        }
        if nodes - self.since_nodes >= stall_nodes && elapsed_secs - self.since_secs >= stall_secs {
            StallRestartVerdict::Fire
        } else {
            StallRestartVerdict::Wait
        }
    }
}

enum StallRestartVerdict {
    Fire,
    Wait,
}

impl<'a, E: LpEngine> Search<'a, E> {
    fn new(p: &'a MipProblem, mut engine: E, config: &'a MipConfig, start: Instant) -> Self {
        let integral_objective = p
            .cost
            .iter()
            .zip(&p.integer)
            .all(|(&c, &int)| c == 0.0 || (int && (c - c.round()).abs() <= 1e-9));
        // The search's whole time budget as a solve-internal deadline on
        // the engine. `limit_hit` fires only between nodes; without this a
        // single cycling node LP (mik's restart round 2)
        // spins past the wall-clock limit inside one solve. The deadline
        // is constant for the search, so it is set once here; every engine
        // call — root, node, restart-round root, probes — respects it.
        engine.set_deadline(
            config
                .time_limit
                .map(|t| start + Duration::from_secs_f64(t)),
        );
        Self {
            p,
            engine,
            propagator: Propagator::new(p, config.feastol),
            config,
            start,
            integral_objective,
            incumbent: None,
            upper_limit: config
                .objective_cutoff
                .map_or(f64::INFINITY, |c| p.caller_objective(c)),
            nodes: 0,
            lp_iterations: 0,
            first_incumbent: None,
            incumbents: Vec::new(),
            locks: Locks::new(p),
            rng: XorShift64::new(0x5EED ^ config.seed.wrapping_mul(0x9E37_79B9_7F4A_7C15)),
            heur_lp_iterations: 0,
            global_lower: p.col_lower.clone(),
            global_upper: p.col_upper.clone(),
            cut_stats: CutStats::default(),
            separator: None,
            brancher: Brancher::new(
                &p.integer,
                config.integrality_tol,
                config.feastol,
                config.branching,
            ),
            num_integer: p.integer.iter().filter(|&&i| i).count(),
            root_lp_iterations: 0,
            root_lp_bound: f64::NEG_INFINITY,
            node_lp_solves: 0,
            node_lp_iterations: 0,
            prof: TreeProfile::default(),
            sb_start: None,
            root_redcost: RootRedcost::default(),
            redcost_stats: RedcostStats::default(),
            // HiGHS builds its clique table in presolve (HPresolve.cpp:
            // 1494–1512); ours is per search round, which also gives every
            // restart round a fresh table over the reduced problem.
            cliques: config.clique_fixings.then(|| {
                let t = cuts::root_clique_table(p, &p.col_lower, &p.col_upper, config.feastol);
                if !t.is_empty() && std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
                    eprintln!(
                        "root clique table: {} cliques, {} entries",
                        t.num_cliques(),
                        t.num_entries()
                    );
                }
                t
            }),
            clique_stats: cuts::CliqueFixStats::default(),
            probing_stats: crate::probing::ProbingStats::default(),
            submip_stats: SubMipStats::default(),
            submip_secs: 0.0,
            submip_rates: FixingRateStats::default(),
            submip_rng: XorShift64::new(0x5B3D ^ config.seed.wrapping_mul(0x9E37_79B9_7F4A_7C15)),
            localmip_effort: 0,
            localmip_calls: 0,
            localmip_found: 0,
            alns: config.alns.then(|| Alns::new(config.seed)),
            incumbent_node: 0,
            root_secs: 0.0,
            conflicts: ConflictPool::new(p.num_col, p.integer.iter().filter(|&&i| i).count()),
            conflict_stats: ConflictStats::default(),
            global_infeasible: false,
            race: None,
            restarts_done: 0,
            no_lp_branches: 0,
            unsolved_pruned: 0,
            restart_armed: false,
            restart_arm_limit: f64::INFINITY,
            stall_restart: None,
        }
    }

    /// Rounding heuristics on a fractional node-LP point `x`. Simple
    /// rounding (no LP) runs everywhere; at the root ZI round, randomized
    /// rounding and shifting follow (HiGHS evaluateRootNode order,
    /// HighsMipSolverData.cpp:2110–2115); below the root, randomized rounding
    /// only while there is no incumbent (HiGHS dive start,
    /// HighsMipSolver.cpp:766) and the LP budget allows: 10 000 iterations
    /// plus `effort` × all LP iterations (HighsMipSolverData.cpp:683–687).
    /// Fractional then coefficient diving run within the same budget at
    /// depths 0, 10, 20, … (SCIP's default diving frequency 10, Achterberg
    /// 2007 §9.2).
    // Heuristic LPs do share the node engine, and do move its `last` off the
    // node basis — but that costs nothing today, so no separate workspace was
    // built for them: every tree LP
    // carries its own start basis (children `lp.basis`; probes via `sb_start`,
    // with `LpEngine::probe` restoring `last`; the reduced-cost re-solve from
    // `lp.basis`, and `redcost_node` runs before this call), `duals()` is keyed
    // on `Rc::ptr_eq` so a stale snapshot costs a transpose re-solve at worst,
    // and the shared `StdCache` refresh measures std=4 ms against a ~0.8 s LP
    // census (loop=595, setup=135). The
    // state only becomes load-bearing if queued nodes start from `last`, which
    // was measured to lose 1.60x/2.33x on our benchmark set.
    fn lp_heuristics(&mut self, node: &Node, x: &[f64]) {
        let hc = self.config.heuristics;
        let budget =
            hc.effort * self.lp_iterations as f64 + 10_000.0 - self.heur_lp_iterations as f64;
        let budget_ok = budget > 0.0;
        let root = node.depth == 0;
        let below_root = !root && hc.randomized_rounding && self.incumbent.is_none() && budget_ok;
        let dive = hc.diving && budget_ok && node.depth.is_multiple_of(10);
        if !(hc.simple_rounding || root || below_root || dive) {
            return;
        }
        let mut ctx = Ctx {
            p: self.p,
            locks: &self.locks,
            prop: &mut self.propagator,
            engine: &mut self.engine,
            lower: &self.global_lower,
            upper: &self.global_upper,
            cutoff: self.upper_limit,
            feastol: self.config.feastol,
            rng: &mut self.rng,
            lp_iterations: 0,
            found: Vec::new(),
            deadline: self
                .config
                .time_limit
                .map(|t| self.start + Duration::from_secs_f64(t)),
        };
        if hc.simple_rounding {
            ctx.simple_rounding(x);
        }
        if root {
            if hc.zi_round {
                ctx.zi_round(x);
            }
            if hc.randomized_rounding {
                ctx.randomized_rounding(x);
            }
            if hc.shifting {
                ctx.shifting(x);
            }
        } else if below_root {
            ctx.randomized_rounding(x);
        }
        if dive {
            for rule in [DiveRule::Fractional, DiveRule::Coefficient] {
                let left = budget as u64 - ctx.lp_iterations.min(budget as u64);
                let _ = ctx.dive(rule, &node.lower, &node.upper, x, left);
            }
        }
        let (found, iters) = (ctx.found, ctx.lp_iterations);
        self.heur_lp_iterations += iters;
        self.lp_iterations += iters;
        for (y, source) in found {
            self.try_incumbent(y, source);
        }
    }

    /// The root sequence between the first LP and the cut loop (HiGHS
    /// evaluateRootNode, HighsMipSolverData.cpp:2109–2119, MIT): ZI round,
    /// randomized rounding and shifting on the first LP point, then
    /// `evaluateRootLp` (:1870–1898), which propagates the global domain
    /// against the new cutoff and re-solves the LP if bounds moved. The cut
    /// loop then separates on the tightened domain, so big-M rows over
    /// otherwise unbounded continuous columns (job-shop start times bounded
    /// only through the makespan) give finite cut coefficients. Cuts derived
    /// there are valid for every improving solution, as in HiGHS. Returns the
    /// root LP to separate on, re-solved warm from the first LP's basis; an
    /// infeasible or cut-off re-solve prunes the root.
    ///
    /// Unlike HiGHS (one trial), randomized rounding is repeated up to
    /// `ROOT_RR_TRIALS` times while there is no incumbent, stopping at the
    /// first point found or when the heuristic LP budget is spent: repeated
    /// independent trials are how randomized rounding boosts its success
    /// probability (Raghavan & Thompson 1987, Combinatorica 7(4)).
    fn root_incumbent(&mut self, node: &mut Node, lp: LpOutcome) -> Result<LpOutcome, MipError> {
        let hc = self.config.heuristics;
        let budget =
            hc.effort * self.lp_iterations as f64 + 10_000.0 - self.heur_lp_iterations as f64;
        let trials = if self.incumbent.is_none() {
            ROOT_RR_TRIALS
        } else {
            1
        };
        let mut heuristic_lps = false;
        if self
            .most_fractional(&lp.x, &node.lower, &node.upper)
            .is_some()
        {
            let mut ctx = Ctx {
                p: self.p,
                locks: &self.locks,
                prop: &mut self.propagator,
                engine: &mut self.engine,
                lower: &self.global_lower,
                upper: &self.global_upper,
                cutoff: self.upper_limit,
                feastol: self.config.feastol,
                rng: &mut self.rng,
                lp_iterations: 0,
                found: Vec::new(),
                deadline: self
                    .config
                    .time_limit
                    .map(|t| self.start + Duration::from_secs_f64(t)),
            };
            if hc.zi_round {
                ctx.zi_round(&lp.x);
            }
            if hc.randomized_rounding {
                for _ in 0..trials {
                    ctx.randomized_rounding(&lp.x);
                    if !ctx.found.is_empty() || ctx.lp_iterations as f64 >= budget {
                        break;
                    }
                }
            }
            if hc.shifting {
                ctx.shifting(&lp.x);
            }
            let (found, iters) = (ctx.found, ctx.lp_iterations);
            heuristic_lps = iters > 0;
            self.heur_lp_iterations += iters;
            self.lp_iterations += iters;
            for (y, source) in found {
                self.try_incumbent(y, source);
            }
        }
        let mut moved = false;
        if self.incumbent.is_some() {
            let (lower, upper) = (node.lower.clone(), node.upper.clone());
            // The incumbents just found applied the first root LP's reduced
            // costs to the global box (`new_incumbent`); propagate them too
            // (HiGHS evaluateRootLp, HighsMipSolverData.cpp:1946–1953).
            if self.global_infeasible || self.intersect_global(node).is_none() {
                return Ok(LpOutcome {
                    status: LpStatus::Infeasible,
                    ..lp
                });
            }
            if !self
                .propagator
                .propagate(&mut node.lower, &mut node.upper, None, self.upper_limit)
            {
                return Ok(LpOutcome {
                    status: LpStatus::Infeasible,
                    ..lp
                });
            }
            self.global_lower.clone_from(&node.lower);
            self.global_upper.clone_from(&node.upper);
            moved = node.lower != lower || node.upper != upper;
        }
        // The heuristic LPs left the engine's last basis on a rounded point,
        // and the cut loop's re-solves start from that basis; re-solving warm
        // from the root basis (0 pivots if no bound moved) restores it.
        if !moved && !(heuristic_lps && lp.basis.is_some()) {
            return Ok(lp);
        }
        let lp = self.engine.solve_from(
            &node.lower,
            &node.upper,
            lp.basis.as_deref(),
            self.upper_limit,
        )?;
        self.lp_iterations += lp.iterations;
        Ok(lp)
    }

    /// The root box and LP after a cutoff improved: intersects the global
    /// box into the node, propagates it, and re-solves warm if it moved.
    /// `None` = the root is pruned.
    fn root_after_cutoff(
        &mut self,
        node: &mut Node,
        lp: LpOutcome,
    ) -> Result<Option<LpOutcome>, MipError> {
        let (lower, upper) = (node.lower.clone(), node.upper.clone());
        if self.global_infeasible || self.intersect_global(node).is_none() {
            return Ok(None);
        }
        if !self
            .propagator
            .propagate(&mut node.lower, &mut node.upper, None, self.upper_limit)
        {
            return Ok(None);
        }
        self.global_lower.clone_from(&node.lower);
        self.global_upper.clone_from(&node.upper);
        if node.lower == lower && node.upper == upper {
            return Ok(Some(lp));
        }
        let out = self.engine.solve_from(
            &node.lower,
            &node.upper,
            lp.basis.as_deref(),
            self.upper_limit,
        )?;
        self.lp_iterations += out.iterations;
        match out.status {
            LpStatus::Infeasible | LpStatus::Cutoff => Ok(None),
            _ => Ok(Some(out)),
        }
    }

    /// Heuristics that need no LP, run once on the root box before the root
    /// LP (HiGHS order, HighsMipSolver.cpp:187–209: trivial heuristics, then
    /// Feasibility Jump).
    fn pre_root_heuristics(&mut self, root: &Node) {
        let hc = self.config.heuristics;
        if hc.trivial {
            for x in heuristics::trivial(self.p, &root.lower, &root.upper) {
                self.try_incumbent(x, Source::Trivial);
            }
        }
        if hc.feasibility_jump {
            let start = self.incumbent.as_ref().map(|inc| inc.1.clone());
            let mut cfg = heuristics::fj::FjConfig::new(self.config.feastol);
            cfg.seed = self.config.seed;
            // The co-100 instance spent 30 s past a
            // 60 s limit in FJ restarts before the root LP started.
            cfg.deadline = self.deadline();
            if let Some(r) = heuristics::fj::feasibility_jump(
                self.p,
                &root.lower,
                &root.upper,
                start.as_deref(),
                cfg,
            ) {
                let restart = (hc.fj_restart_effort > 0.0).then(|| r.clone());
                self.try_incumbent(r.x, Source::FeasibilityJump);
                if let Some(first) = restart {
                    let b = heuristics::fj::restart_from_best(
                        self.p,
                        &root.lower,
                        &root.upper,
                        first,
                        cfg,
                        hc.fj_restart_effort,
                    );
                    self.try_incumbent(b.x, Source::FeasibilityJump);
                }
            }
        }
    }

    /// `--localmip`: polish the just-installed incumbent with the
    /// [`heuristics::localmip`] improver
    /// (Lin, Zou & Cai 2024). Called from
    /// `new_incumbent` for every incumbent except the improver's own output
    /// (it continues improving internally), so it runs off incumbent-like
    /// points as the design intends and complements the finders — FJ,
    /// rounding, diving, RINS — instead of competing with them. Budget in
    /// FJ's effort currency ([`heuristics::fj::effort_cap`] multiples):
    /// [`heuristics::localmip::POLISH_EFFORT_MULT`] per call,
    /// [`heuristics::localmip::TOTAL_EFFORT_MULT`] over the search; a call
    /// that returns no point is charged its full budget, so the total caps
    /// the worst case. The found point goes through the same incumbent path
    /// (`try_incumbent`) every heuristic uses.
    fn run_localmip(&mut self) {
        if !self.config.localmip || self.incumbent.is_none() {
            return;
        }
        if self.race.as_ref().is_some_and(Link::stopped) {
            return;
        }
        let cap = heuristics::fj::effort_cap(self.p);
        let left =
            (heuristics::localmip::TOTAL_EFFORT_MULT * cap).saturating_sub(self.localmip_effort);
        let budget = (heuristics::localmip::POLISH_EFFORT_MULT * cap).min(left);
        if budget == 0 {
            return;
        }
        self.localmip_calls += 1;
        let start = self.incumbent.as_ref().unwrap().1.clone();
        let out = heuristics::localmip::local_mip(self.p, &start, budget, self.rng.next_u64());
        self.localmip_effort += out.as_ref().map_or(budget, |r| r.effort.min(budget));
        if std::env::var_os("RUSOLVE_LOCALMIP_LOG").is_some() {
            eprintln!(
                "localmip: calls={} found={} effort={}/{} budget={budget} elapsed={:.3}s",
                self.localmip_calls,
                self.localmip_found,
                self.localmip_effort,
                heuristics::localmip::TOTAL_EFFORT_MULT * cap,
                self.start.elapsed().as_secs_f64(),
            );
        }
        if let Some(r) = out
            && self.try_incumbent(r.x, Source::LocalMip)
        {
            self.localmip_found += 1;
        }
    }

    /// Audit: one `audit tree` line per search round — node-LP solves,
    /// iterations, wall time and mean/max LP rows, strong-branching time,
    /// and the node-separation cost split.
    fn print_tree_profile(&self) {
        let pr = &self.prof;
        let c = &self.cut_stats;
        let l = self.engine.stats();
        let ms = |ns: u64| ns as f64 * 1e-6;
        let n = pr.solves.max(1) as f64;
        let (lp_cuts, pool) = self.separator.as_ref().map_or((0, 0), |s| s.pool_sizes());
        eprintln!(
            "audit tree node_lps={} node_iters={} node_lp_ms={:.1} rows_mean={:.1} \
             rows_max={} cut_rows_mean={:.1} probe_ms={:.1} probe_iters={} \
             tree_calls={} tree_rounds={} tree_applied={} tree_lp_iters={} \
             sep_ms={:.1} pool_ms={:.1} resolve_ms={:.1} age_ms={:.1} aged_rows={} \
             root_sep_ms={:.1} root_pool_ms={:.1} root_resolve_ms={:.1} \
             lp_cuts_end={lp_cuts} pool_end={pool} unfit_lps={} unfit_iters={} \
             queue_lps={} queue_iters={} t={:.6}",
            pr.solves,
            self.node_lp_iterations,
            ms(pr.nanos),
            pr.rows_sum as f64 / n,
            pr.rows_max,
            pr.cut_rows_sum as f64 / n,
            ms(l.probe_nanos),
            l.probe_iterations,
            c.tree_calls,
            c.tree_rounds,
            c.tree_applied,
            c.tree_lp_iterations,
            ms(c.tree_sep_nanos),
            ms(c.tree_pool_nanos),
            ms(c.tree_resolve_nanos),
            ms(c.tree_age_nanos),
            c.tree_aged_rows,
            ms(c.sep_nanos),
            ms(c.pool_nanos),
            ms(c.resolve_nanos),
            pr.unfit_lps,
            pr.unfit_iters,
            pr.queue_lps,
            pr.queue_iters,
            self.start.elapsed().as_secs_f64()
        );
    }

    /// HiGHS `moreHeuristicsAllowed` (HighsMipSolverData.cpp:672–731): the
    /// sub-MIP heuristics share the LP-iteration budget of the other
    /// heuristics, 10 000 plus `effort` × all LP iterations early in the
    /// search (:683–687). HiGHS's later branch extrapolates the tree size
    /// from the pruned tree weight, which this tree does not track; there
    /// the budget is the plain `effort` share (the sub-MIP form, :679).
    fn more_heuristics_allowed(&self) -> bool {
        let effort = self.config.heuristics.effort;
        let (heur, total) = (self.heur_lp_iterations as f64, self.lp_iterations as f64);
        if self.nodes < 1000 {
            heur < total * effort + 10_000.0
        } else {
            heur < total * effort
        }
    }

    /// One RENS (no incumbent, or `root`) or RINS call on the LP point `x`
    /// (HighsMipSolver.cpp:772–797; root RENS HighsMipSolverData.cpp:2396–
    /// 2400), then the sub-MIP's point through the incumbent path every
    /// heuristic uses. Limits follow HiGHS: node limit `200 + nodes / 20`
    /// (HighsPrimalHeuristics.cpp:632, :943), the remaining time, the
    /// incumbent as cutoff, and no recursion (sub-MIPs run with sub-MIP
    /// heuristics off and one thread).
    fn run_submip(&mut self, x: &[f64], root: bool) {
        let remaining = match self.config.time_limit {
            Some(t) => {
                let left = t - self.start.elapsed().as_secs_f64();
                if left <= 0.0 {
                    return;
                }
                Some(left)
            }
            None => None,
        };
        if self.race.as_ref().is_some_and(Link::stopped) {
            return;
        }
        // Cost-aware cap on NODE calls (ours, not HiGHS): a node call
        // would otherwise be free to run its whole `200 + nodes/20` node
        // limit — seconds of cold-engine work on a large parent, as
        // measured on int_production (a 12.6 s solve pushed past a 20 s
        // limit by single late RINS calls on master scheduling). The root
        // call keeps HiGHS's `remaining` limit; its cost is handled by
        // charging, not capping (call site above).
        let cap = if root {
            remaining
        } else {
            let left =
                submip::time_share_left(self.submip_secs, self.start.elapsed().as_secs_f64());
            Some(remaining.unwrap_or(f64::INFINITY).min(left.max(0.0)))
        };
        let rins = self.incumbent.is_some() && !root;
        let cutoff = self.incumbent.as_ref().map(|i| i.0);
        let cfg = SubMipConfig {
            node_limit: 200 + self.nodes / 20,
            time_limit: cap,
            target_fixing_rate: self.submip_rates.target(&mut self.submip_rng),
            cutoff,
            seed: self.submip_rng.next_u64(),
            mip: MipConfig {
                submip: false,
                alns: false,
                // The Local-MIP improver polishes in the outer search; a
                // hard time-boxed sub-MIP spending its seconds on LS polish
                // instead of the restriction solve is the measured anti-pattern.
                localmip: false,
                // Root-only probing — a sub-MIP is hard time-boxed
                // (SUBMIP_TIME_SHARE) and a probe contingent inside it
                // would starve the actual solve (HiGHS probes sub-MIPs
                // with tightened limits; not ported, root scope).
                root_probing: false,
                // HiGHS disables restarts inside sub-MIPs (`!mipsolver.submip`,
                // HighsMipSolverData.cpp:2487).
                restart: false,
                threads: 1,
                time_limit: cap,
                node_limit: None,
                objective_cutoff: None,
                audit: false,
                ..*self.config
            },
            ..SubMipConfig::default()
        };
        let call_start = Instant::now();
        let report = match (&self.incumbent, rins) {
            (Some((_, inc)), true) => submip::rins_report(self.p, inc, x, &cfg),
            _ => submip::rens_report(self.p, x, &cfg),
        };
        // Charge the call's wall seconds (setup included): the iteration
        // budget above cannot see restriction build, presolve and engine
        // setup, which is what made sub-MIPs expensive on cheap parents
        // (submip::SUBMIP_TIME_SHARE).
        let call_secs = call_start.elapsed().as_secs_f64();
        self.submip_secs += call_secs;
        self.submip_stats.seconds += call_secs;
        if std::env::var_os("RUSOLVE_SUBMIP_LOG").is_some() {
            eprintln!(
                "submip root={root} rins={rins} nodes={} outcome={:?} rate={:.2} found={} cutoff={cutoff:?} secs={:.4} tot={:.3}s elapsed={:.3}s",
                self.nodes,
                report.outcome,
                report.fixing_rate,
                report.x.is_some(),
                call_secs,
                self.submip_secs,
                self.start.elapsed().as_secs_f64(),
            );
        }
        self.submip_rates.observe(&report);
        let st = &mut self.submip_stats;
        if matches!(
            report.outcome,
            submip::SubMipOutcome::Skipped | submip::SubMipOutcome::FixingsInfeasible
        ) {
            st.skipped += 1;
        } else if rins {
            st.rins_calls += 1;
        } else {
            st.rens_calls += 1;
        }
        // HiGHS scales the sub-MIP's LP iterations by the share of columns
        // it kept free (HighsPrimalHeuristics.cpp:216–222).
        let ints = self.num_integer as f64;
        let conts = (self.p.num_col as f64) - ints;
        let free = ((1.0 - report.fixing_rate) * ints + conts) / (ints + conts).max(1.0);
        let iters = (free * report.lp_iterations as f64) as u64;
        st.nodes += report.nodes;
        st.lp_iterations += iters;
        self.heur_lp_iterations += iters;
        self.lp_iterations += iters;
        if let Some(y) = report.x
            && self.try_incumbent(y, Source::SubMip)
        {
            if rins {
                self.submip_stats.rins_found += 1;
            } else {
                self.submip_stats.rens_found += 1;
            }
        }
    }

    /// `--root-dive`: the LP-guided dive that builds the root sub-MIP's
    /// fixing point. Dives (coefficient rule, warm node LPs, one-level
    /// backtracking — HiGHS reaches its RENS/RINS neighbourhood by exactly
    /// this dive, `HighsPrimalHeuristics.cpp:394+`; our static fixing was
    /// the documented deviation, submip.rs) from the post-cut root LP
    /// point, installs every point its intermediate simple roundings find,
    /// and returns the dive's final integral point for the sub-MIP to
    /// polish. `None` = the dive failed or hit the heuristic LP budget:
    /// the caller falls back to the raw LP point (the static-fixing shape).
    fn root_dive(&mut self, node: &Node, x: &[f64]) -> Option<Vec<f64>> {
        self.most_fractional(x, &node.lower, &node.upper)?;
        let hc = self.config.heuristics;
        let budget =
            hc.effort * self.lp_iterations as f64 + 10_000.0 - self.heur_lp_iterations as f64;
        if budget <= 0.0 {
            return None;
        }
        let mut ctx = Ctx {
            p: self.p,
            locks: &self.locks,
            prop: &mut self.propagator,
            engine: &mut self.engine,
            lower: &self.global_lower,
            upper: &self.global_upper,
            cutoff: self.upper_limit,
            feastol: self.config.feastol,
            rng: &mut self.rng,
            lp_iterations: 0,
            found: Vec::new(),
            deadline: self
                .config
                .time_limit
                .map(|t| self.start + Duration::from_secs_f64(t)),
        };
        let point = ctx.dive(
            DiveRule::Coefficient,
            &node.lower,
            &node.upper,
            x,
            budget as u64,
        );
        let (found, iters) = (ctx.found, ctx.lp_iterations);
        self.heur_lp_iterations += iters;
        self.lp_iterations += iters;
        for (y, source) in found {
            self.try_incumbent(y, source);
        }
        point
    }

    /// One ALNS call on the node LP point `x` if its budget allows
    /// ([`crate::alns`]); `lower` is the global dual bound, for the reward.
    /// Accounting follows [`Self::run_submip`].
    fn run_alns(&mut self, x: &[f64], lower: f64) {
        if self.race.as_ref().is_some_and(Link::stopped) {
            return;
        }
        let elapsed = self.start.elapsed().as_secs_f64();
        let remaining = self.config.time_limit.map(|t| t - elapsed);
        if remaining.is_some_and(|r| r <= 0.0) {
            return;
        }
        let Some(mut a) = self.alns.take() else {
            return;
        };
        let call = alns::Call {
            lp_x: x,
            incumbent: self.incumbent.as_ref().map(|(z, y)| (*z, y.as_slice())),
            parent_nodes: self.nodes,
            nodes_since_incumbent: self
                .incumbent
                .as_ref()
                .map(|_| self.nodes - self.incumbent_node),
            elapsed,
            remaining,
            root_secs: self.root_secs,
        };
        let run = a.run(self.p, &call, self.config);
        let Some(mut run) = run else {
            self.alns = Some(a);
            return;
        };
        let arm = alns::Arm::ALL.iter().position(|&k| k == run.arm).unwrap();
        let r = &run.report;
        let ints = self.num_integer as f64;
        let conts = (self.p.num_col as f64) - ints;
        let free = ((1.0 - r.fixing_rate) * ints + conts) / (ints + conts).max(1.0);
        let iters = (free * r.lp_iterations as f64) as u64;
        let st = &mut self.submip_stats;
        if matches!(
            r.outcome,
            submip::SubMipOutcome::Skipped | submip::SubMipOutcome::FixingsInfeasible
        ) {
            st.skipped += 1;
        } else {
            st.alns_calls[arm] += 1;
        }
        st.nodes += r.nodes;
        st.lp_iterations += iters;
        st.alns_seconds += run.secs;
        self.heur_lp_iterations += iters;
        self.lp_iterations += iters;
        if std::env::var_os("RUSOLVE_SUBMIP_LOG").is_some() {
            eprintln!(
                "alns node={} arm={} nodes={} outcome={:?} rate={:.2} found={} secs={:.4}",
                self.nodes,
                run.arm.name(),
                r.nodes,
                r.outcome,
                r.fixing_rate,
                r.x.is_some(),
                run.secs
            );
        }
        let improved = match run.report.x.take() {
            Some(y) => self.try_incumbent(y, Source::SubMip),
            None => false,
        };
        if improved {
            self.submip_stats.alns_found[arm] += 1;
            // `new_incumbent` ran with the state taken out: pool it here.
            if let Some((z, y)) = &self.incumbent {
                a.add_solution(*z, y);
            }
        }
        let upper = self.incumbent.as_ref().map_or(f64::INFINITY, |i| i.0);
        a.observe(&run, improved, upper, lower);
        self.alns = Some(a);
    }

    /// The incumbent check for heuristic points: finite, within bounds and
    /// rows by `CANDIDATE_TOL`, integral by `integrality_tol` (integer
    /// entries are then snapped exactly). Rejected points are dropped.
    fn try_incumbent(&mut self, mut x: Vec<f64>, source: Source) -> bool {
        if x.len() != self.p.num_col || x.iter().any(|v| !v.is_finite()) {
            return false;
        }
        let v = self.p.violation(&x);
        if v.bound > CANDIDATE_TOL
            || v.row > CANDIDATE_TOL
            || v.integrality > self.config.integrality_tol
        {
            return false;
        }
        for (xj, &int) in x.iter_mut().zip(&self.p.integer) {
            if int {
                *xj = xj.round();
            }
        }
        self.new_incumbent(x, source)
    }

    fn run(&mut self, root: Node) -> Result<TreeEnd, MipError> {
        let mut tree: NodeSelector<Node> = NodeSelector::new(self.config.branching.node_selection);
        let mut plunge = Some(root);
        let mut pruned_at = self.upper_limit;
        loop {
            let own = self.incumbent.as_ref().map(|i| i.0);
            if let Some(x) = self.race.as_mut().and_then(|r| r.pull(own)) {
                self.try_incumbent(x, Source::Racer);
            }
            if self.global_infeasible {
                return Ok(TreeEnd::Exhausted);
            }
            if self.upper_limit < pruned_at {
                pruned_at = self.upper_limit;
                tree.prune(pruned_at);
            }
            // `lineage`: this node continues the node processed last as its
            // direct child (the plunge path). Only then are the conflict
            // watches guaranteed to have last been set along this box's own
            // lineage, which is what makes `Revisit::Columns` sound; a node
            // from the queue or the plunge stack may come from any lineage
            // and re-checks the whole pool.
            let (mut node, from_queue, lineage) = match plunge.take() {
                Some(n) => {
                    let lin = n.side.is_some();
                    (n, false, lin)
                }
                None => {
                    self.brancher.new_dive();
                    let next = tree.next(
                        &self.brancher.pseudocost,
                        self.nodes,
                        self.upper_limit,
                        self.config.feastol,
                    );
                    match next {
                        Some((n, nodequeue::Source::Queue)) => {
                            self.brancher.new_subtree();
                            // HiGHS ages its conflict pool once per search
                            // pass (HighsMipSolver.cpp:839).
                            self.conflict_stats.aged_out += self.conflicts.age() as u64;
                            (n, true, false)
                        }
                        Some((n, nodequeue::Source::Plunge)) => (n, false, false),
                        None => return Ok(TreeEnd::Exhausted),
                    }
                }
            };
            if node.bound > self.upper_limit {
                continue;
            }
            if self.nodes == 1 {
                self.root_secs = self.start.elapsed().as_secs_f64();
            }
            // Only a tree that outlives its root starts the other racers.
            if self.nodes == 1
                && let Some(race) = &self.race
            {
                race.root_done();
            }
            if let Some(status) = self.limit_hit() {
                let open_bound = tree.open_bound().min(node.bound);
                return Ok(TreeEnd::Limit(status, open_bound));
            }
            self.nodes += 1;
            // Stall restart : the root
            // stashed the no-fixings restart shape (wide gap, cut mass) and
            // the tree has now demonstrated it is not clean — the incumbent
            // has not significantly moved for a full stall window. Checked
            // at every non-root node (a small tree can be one plunge whose
            // queue siblings die unpruned-popped; the restart discards the
            // whole tree either way, so the boundary does not matter). The
            // verdict re-opens the window on every significant improvement,
            // so a tree still making primal progress never restarts
            // (Anderson et al. 2019 restart the same way: on demonstrated
            // size, not on root shape).
            if node.depth > 0
                && self.restart_armed
                && self.restarts_done == 0
                && self.stall_restart.is_some()
            {
                let bound = tree.open_bound().min(node.bound);
                let z = self.incumbent.as_ref().map(|i| i.0);
                let secs = self.start.elapsed().as_secs_f64();
                let verdict = self.stall_restart.as_mut().unwrap().verdict(
                    z,
                    self.nodes,
                    secs,
                    self.config.restart_stall_nodes,
                    self.config.restart_stall_secs,
                );
                if matches!(verdict, StallRestartVerdict::Fire) {
                    let lb = self.stall_restart.take().unwrap();
                    if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
                        eprintln!(
                            "stall restart at node {}: incumbent {:?} stalled since node {}, gap {:.4} over {:.4}",
                            self.nodes,
                            z,
                            lb.since_nodes,
                            self.root_gap(bound),
                            bound
                        );
                    }
                    return Ok(TreeEnd::Restart(Box::new(
                        self.stall_restart_request(&lb, bound),
                    )));
                }
            }
            // Global tightenings since the node was created (root reduced-
            // cost fixing) enter its box and seed its propagation.
            let Some(mut seeds) = self.intersect_global(&mut node) else {
                continue;
            };
            let before = node.side.map(|_| (node.lower.clone(), node.upper.clone()));
            seeds.extend(node.branched);
            let feasible = if seeds.len() <= 1 {
                self.propagator.propagate(
                    &mut node.lower,
                    &mut node.upper,
                    node.branched,
                    self.upper_limit,
                )
            } else {
                self.propagator.propagate_cols(
                    &mut node.lower,
                    &mut node.upper,
                    &seeds,
                    self.upper_limit,
                )
            };
            // Watch-list propagation: a plunge child revisits only conflicts watching a
            // column that moved since its parent's box was propagated —
            // the moves inherited from the parent, the columns the global
            // intersect moved into its box, and its own propagation's. A
            // dense revisit set is a full scan in disguise: traversing and
            // sorting the watch lists then costs more than the guarded
            // slot walk, so it degrades to `Revisit::All`.
            let feasible = feasible
                && if lineage {
                    let mut revisit = std::mem::take(&mut node.inherited);
                    if let Some((lo, hi)) = &before {
                        moved_cols(lo, hi, &node.lower, &node.upper, &mut revisit);
                    }
                    revisit.extend(seeds.iter().copied());
                    revisit.sort_unstable();
                    revisit.dedup();
                    if 2 * revisit.len() >= self.p.num_col {
                        self.conflict_propagate(&mut node, Revisit::All)
                    } else {
                        self.conflict_propagate(&mut node, Revisit::Columns(&revisit))
                    }
                } else {
                    self.conflict_propagate(&mut node, Revisit::All)
                };
            // The box this node's conflict propagation finished on, when
            // there was one: the base for what its children inherit. With
            // it the inherited set is only this node's post-propagation
            // moves (reduced-cost fixing on its LP); without the snapshot
            // the parent's own propagation footprint would enter every
            // child's revisit set, and a wide propagated node would make it
            // dense all the way down a dive.
            let invariant = if lineage && feasible && self.conflicts.len() > 1 {
                Some((node.lower.clone(), node.upper.clone()))
            } else {
                None
            };
            if let (Some((side, _)), Some((lo, hi))) = (node.side, &before) {
                let ps = &mut self.brancher.pseudocost;
                if feasible {
                    let inferences = (0..lo.len())
                        .filter(|&j| {
                            j != side.col && (lo[j] != node.lower[j] || hi[j] != node.upper[j])
                        })
                        .count();
                    ps.add_inference_observation(side.col, inferences, side.up);
                } else {
                    ps.add_cutoff_observation(side.col, side.up);
                }
            }
            // Columns the parent's box moved after its conflict propagation
            // (reduced-cost fixing on this node's LP), for the children's
            // conflict revisit set (see `Node::inherited`).
            let inh_cols = |node: &Node| -> Vec<usize> {
                let (lo, hi) = match invariant.as_ref().or(before.as_ref()) {
                    Some(b) => b,
                    None => return Vec::new(),
                };
                let mut v = Vec::new();
                moved_cols(lo, hi, &node.lower, &node.upper, &mut v);
                v
            };
            if !feasible {
                continue;
            }
            if node.depth == 0 {
                self.global_lower.clone_from(&node.lower);
                self.global_upper.clone_from(&node.upper);
            }
            // Warm dual reopt from the parent basis, stopped once the bound
            // passes the upper limit. Below
            // the root a cold fallback is capped like every tree LP.
            if node.depth > 0 {
                let cap = self.tree_lp_iteration_cap();
                self.engine.set_iteration_limit(Some(cap));
            }
            // The root node's cold LPs may take the dual, and with a
            // time limit a root LP that fails (iteration cap, numerical
            // give-up) re-runs once on the other simplex with no cap (HiGHS:
            // no root iteration limit, HighsMipSolverData.cpp:2018; strategy
            // toggle on failure, HighsLpRelaxation.cpp:1236-1251). A root
            // still unfinished at the deadline ends as a clean TimeLimit.
            self.engine.set_root_lp(
                node.depth == 0 && self.config.root_lp_dual,
                node.depth == 0 && self.config.time_limit.is_some(),
            );
            let prof_t0 = (self.config.audit && node.depth > 0).then(Instant::now);
            let prof_unfit = {
                let l = self.engine.stats();
                l.unfit_recovered + l.fallback_unfit
            };
            let solved = self.engine.solve_from(
                &node.lower,
                &node.upper,
                node.warm.as_deref(),
                self.upper_limit,
            );
            if let Some(t0) = prof_t0 {
                let l = self.engine.stats();
                let iters = solved.as_ref().map_or(0, |o| o.iterations);
                if l.unfit_recovered + l.fallback_unfit > prof_unfit {
                    self.prof.unfit_lps += 1;
                    self.prof.unfit_iters += iters;
                }
                if from_queue {
                    self.prof.queue_lps += 1;
                    self.prof.queue_iters += iters;
                }
                self.prof.node_lp(
                    t0.elapsed().as_nanos() as u64,
                    self.engine.num_lp_rows(),
                    self.separator.as_ref().map_or((0, 0), |s| s.pool_sizes()),
                    self.start.elapsed().as_secs_f64(),
                );
            }
            self.engine.set_iteration_limit(None);
            let mut lp = match solved {
                // A stalled node LP: branch without it instead of aborting
                // the search (HighsSearch.cpp:1347–1418). A numerical
                // give-up too — HiGHS branches on every unsolved status
                // alike — and with every integer column fixed the node gets
                // HiGHS's fallback re-solve instead of an error.
                Err(e @ (LpError::IterationLimit | LpError::Simplex(_))) if node.depth > 0 => {
                    if e == LpError::IterationLimit {
                        self.lp_iterations += self.tree_lp_iteration_cap();
                    }
                    match self.unsolved_node(&node, inh_cols(&node))? {
                        Unsolved::Branch(kids) => {
                            let (first, second) = *kids;
                            tree.push_sibling(second);
                            plunge = Some(first);
                            continue;
                        }
                        Unsolved::Solved(lp) => lp,
                        Unsolved::Pruned => continue,
                    }
                }
                lp => lp?,
            };
            self.lp_iterations += lp.iterations;
            if node.depth == 0 {
                self.root_lp_iterations = lp.iterations;
                if lp.status == LpStatus::Optimal && self.root_lp_bound == f64::NEG_INFINITY {
                    self.root_lp_bound = lp.objective;
                }
                if self.config.audit {
                    eprintln!(
                        "audit root_lp status={:?} obj={:e} iters={} t={:.6}",
                        lp.status,
                        lp.objective,
                        lp.iterations,
                        self.start.elapsed().as_secs_f64()
                    );
                }
            } else {
                self.node_lp_solves += 1;
                self.node_lp_iterations += lp.iterations;
            }
            self.observe_branching(&node, &lp);
            if node.depth > 0
                && self.config.cuts.tree
                && lp.status == LpStatus::Optimal
                && let Some(sep) = self.separator.as_mut()
            {
                // HiGHS ages the LP relaxation's cut rows after every
                // optimal node LP — dive children included
                // (HighsSearch.cpp:958, performAging after resolveLp) — not
                // only where it separates, so rows slack for k node LPs
                // leave even though separation runs at backtracked queue
                // nodes only. Without this the rows of a plunge's first
                // separated node ride along for the whole dive.
                sep.age_node(
                    &mut self.engine,
                    &lp.x,
                    self.config.feastol,
                    self.config.cuts.tree_age_ticks,
                    &mut self.cut_stats,
                );
            }
            if node.depth == 0 {
                self.add_root_redcost(&node, &lp);
            }
            if node.depth == 0 && lp.status == LpStatus::Optimal && lp.objective <= self.upper_limit
            {
                lp = self.root_incumbent(&mut node, lp)?;
            }
            // Cut log for a root restart:
            // every cut the pool accepted; `active_cuts` keeps the tight ones
            // at the final root LP point. Collected only when a restart can
            // fire, to avoid cloning every generated cut otherwise.
            let mut root_cut_log: Vec<(cuts::CutClass, CutRow)> = Vec::new();
            if node.depth == 0 && lp.status == LpStatus::Optimal && lp.objective <= self.upper_limit
            {
                let deadline = self
                    .config
                    .time_limit
                    .and_then(|t| Duration::try_from_secs_f64(t).ok())
                    .and_then(|d| self.start.checked_add(d));
                let race = &self.race;
                let stop = || {
                    deadline.is_some_and(|d| Instant::now() >= d)
                        || race.as_ref().is_some_and(Link::stopped)
                };
                let armed = self.restart_armed;
                // The cut loop's rounding hook (cuts/mod.rs): the
                // root ladder on the loop's own LP points while no cutoff
                // row exists, under the same heuristic LP budget as
                // `root_incumbent`. One RR trial per round: ten-trial rounds
                // were measured and rejected (int_production s1 tree
                // 90,499 -> 139,598 nodes). `cutoff` is the loop's
                // current cutoff-row rhs (INFINITY before the first row),
                // so repair LPs only accept improvements on what the row
                // carries. Points are proposed, verified by the loop (row)
                // and adopted below (incumbent).
                let mut root_heur =
                    |x: &[f64], engine: &mut E, cutoff: f64| -> (Vec<(Vec<f64>, Source)>, u64) {
                        let hc = self.config.heuristics;
                        let budget = hc.effort * self.lp_iterations as f64 + 10_000.0
                            - self.heur_lp_iterations as f64;
                        if budget <= 0.0 || self.incumbent.is_some() {
                            return (Vec::new(), 0);
                        }
                        let mut ctx = Ctx {
                            p: self.p,
                            locks: &self.locks,
                            prop: &mut self.propagator,
                            engine,
                            lower: &self.global_lower,
                            upper: &self.global_upper,
                            cutoff,
                            feastol: self.config.feastol,
                            rng: &mut self.rng,
                            lp_iterations: 0,
                            found: Vec::new(),
                            deadline: self
                                .config
                                .time_limit
                                .map(|t| self.start + Duration::from_secs_f64(t)),
                        };
                        if hc.zi_round {
                            ctx.zi_round(x);
                        }
                        if hc.randomized_rounding {
                            ctx.randomized_rounding(x);
                        }
                        if hc.shifting {
                            ctx.shifting(x);
                        }
                        let (found, iters) = (ctx.found, ctx.lp_iterations);
                        self.heur_lp_iterations += iters;
                        self.lp_iterations += iters;
                        (found, iters)
                    };
                let (out, separator) = if armed {
                    cuts::root_separation_logged(
                        self.p,
                        &mut self.engine,
                        &node.lower,
                        &node.upper,
                        lp,
                        self.config,
                        &stop,
                        &mut self.cut_stats,
                        Some(&mut root_cut_log),
                        self.upper_limit,
                        &mut root_heur,
                    )?
                } else {
                    cuts::root_separation(
                        self.p,
                        &mut self.engine,
                        &node.lower,
                        &node.upper,
                        lp,
                        self.config,
                        &stop,
                        &mut self.cut_stats,
                        self.upper_limit,
                        &mut root_heur,
                    )?
                };
                lp = out;
                self.separator = separator;
                self.lp_iterations += self.cut_stats.lp_iterations;
                if self.config.audit {
                    eprintln!(
                        "audit root_cuts rounds={} status={:?} obj={:e} t={:.6}",
                        self.cut_stats.rounds,
                        lp.status,
                        lp.objective,
                        self.start.elapsed().as_secs_f64()
                    );
                }
                // The cut pool is the root's conflict structure on the
                // root-dive instances (the model rows alone carry none there) —
                // feed the accepted cut rows into the fixings table (HiGHS
                // extractCliquesFromCut) so the incumbent link has cliques
                // to propagate through the rest of the solve.
                if self.config.clique_fixings
                    && let (Some(t), Some(sep)) = (self.cliques.as_mut(), self.separator.as_ref())
                {
                    let before = t.num_cliques();
                    cuts::extend_clique_table_with_cuts(t, sep, self.config.feastol);
                    if t.num_cliques() > before && std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
                        eprintln!(
                            "clique table: +{} from root cuts ({} cliques, {} entries)",
                            t.num_cliques() - before,
                            t.num_cliques(),
                            t.num_entries()
                        );
                    }
                }
                if self.cut_stats.rounds > 0 {
                    self.add_root_redcost(&node, &lp);
                }
                // Points the loop used for its cutoff row become incumbents
                // here (verified once more by `try_incumbent`); an
                // improvement propagates through the root box like the
                // sub-MIP's does.
                let found = std::mem::take(&mut self.cut_stats.root_incumbents);
                let before = self.incumbent.as_ref().map(|i| i.0);
                for (y, source) in found {
                    self.try_incumbent(y, source);
                }
                if self.incumbent.as_ref().map(|i| i.0) != before {
                    lp = match self.root_after_cutoff(&mut node, lp)? {
                        Some(lp) => lp,
                        None => continue,
                    };
                }
                // Post-cut-loop heuristic pass. HiGHS runs RENS/feasibility
                // pump on the root solution after the separation rounds and
                // measures the post-separation fixing rate against that
                // incumbent, gated by `moreHeuristicsAllowed()`
                // (HighsMipSolverData.cpp:2392-2400). Load-bearing for the
                // restart: without it mik restarts off a -534 root point, the
                // restart throws away a working tree, and the base's 42 s solve
                // turns into a 60 s timeout (measured 2026-10-24); with it mik
                // restarts off -46719 and solves in 29 s / 16k nodes.
                //
                // Our design on top of HiGHS's budget gate: the pass runs only
                // when the restart is already live on the *pre-pass* numbers
                // (rate trigger fired, gap over the incumbent above
                // `MIN_ROOT_GAP`). The pass is the restart machinery's
                // incumbent step; where a restart cannot pay for itself, the
                // pass only perturbs the trajectory that branching would have
                // run cleanly (int_production s1: 0 pre-pass fixings, no
                // incumbent, budget-gated pass alone still lost the solve at
                // 20 s vs the no-restart binary's 13.2 s; mik: 59/250 = 23.6%
                // pre-pass fixings and a huge gap).
                if armed
                    && self.cut_stats.rounds > 0
                    && self.more_heuristics_allowed()
                    && !self.brancher.fractional(&lp.x).is_empty()
                    && self.restart_rate_triggered(
                        &self.global_lower,
                        &self.global_upper,
                        &self.p.integer,
                    )
                    && self.root_gap(lp.objective) >= self.config.restart_min_gap
                {
                    lp = self.root_incumbent(&mut node, lp)?;
                }
            } else if from_queue
                && self.config.cuts.tree
                && lp.status == LpStatus::Optimal
                && lp.objective <= self.upper_limit
                && let Some(sep) = self.separator.as_mut()
            {
                // HiGHS separates when a node is taken from the queue, not
                // at plunge children (HighsMipSolver.cpp:684–693, :826–837).
                let before = self.cut_stats.tree_lp_iterations;
                lp = sep.separate_node(
                    self.p,
                    &mut self.engine,
                    &node.lower,
                    &node.upper,
                    lp,
                    self.config,
                    self.upper_limit,
                    self.lp_iterations,
                    &mut self.cut_stats,
                )?;
                self.lp_iterations += self.cut_stats.tree_lp_iterations - before;
            }
            // Stash the incumbent limit at the end of the root cut
            // loop, before the root sub-MIP/dive block below can publish a
            // better point — the restart check's gap bar measures against
            // this (see `restart_arm_limit`).
            if node.depth == 0 {
                self.restart_arm_limit = self.upper_limit;
            }
            if (self.config.submip || self.config.root_dive)
                && node.depth == 0
                && lp.status == LpStatus::Optimal
                && lp.objective <= self.upper_limit
                && (self.incumbent.is_none() || self.more_heuristics_allowed())
                && !self.brancher.fractional(&lp.x).is_empty()
            // The dive (default-on) costs LP iterations only and draws
            // on the same heuristic budget as the rest (the gate above is
            // HiGHS's `moreHeuristicsAllowed` root-RENS gate,
            // HighsMipSolverData.cpp:2392–2400: unconditional while there
            // is no incumbent, budget-gated after). The sub-MIP polish
            // stays on the opt-in `--submip` chain: its seconds pay
            // restriction build + presolve + a cold engine that cheap
            // parents never earn back (sp150x300d 1.59x/air05 1.41x
            // overhead). On that chain the polish
            // runs on the dive's point when the dive reached one; the
            // raw-post-cut-point fallback stays reachable only when the
            // dive is off (`--submip` alone), whose anti-patterns are the
            // documented opt-in risks (int_production s1 steering, air05
            // 16.1 s at equal nodes).
            {
                // HiGHS runs RENS after the root cut loop whatever the
                // incumbent, unless heuristics are out of budget
                // (HighsMipSolverData.cpp:2392–2400). `--root-dive`:
                // build the fixing point by the LP-guided dive first (its
                // integral completions seed the ladder on the way).
                let dived = self
                    .config
                    .root_dive
                    .then(|| self.root_dive(&node, &lp.x))
                    .flatten();
                if self.config.submip {
                    let x = dived.unwrap_or_else(|| lp.x.clone());
                    let before = self.incumbent.as_ref().map(|i| i.0);
                    self.run_submip(&x, true);
                    if self.incumbent.as_ref().map(|i| i.0) != before {
                        // A new cutoff after the cut loop: propagate it into
                        // the root box and re-solve over it, as after the
                        // first-LP heuristics (`root_incumbent`, HiGHS
                        // evaluateRootLp, HighsMipSolverData.cpp:1946–1953).
                        lp = match self.root_after_cutoff(&mut node, lp)? {
                            Some(lp) => lp,
                            None => continue,
                        };
                    }
                }
            }
            match lp.status {
                LpStatus::Infeasible => {
                    self.learn_conflict(&node, false);
                    continue;
                }
                // Cutoff: `lp.objective` is a valid bound above the limit.
                // The dual multipliers that proved it explain the node
                // the way a Farkas ray does (HighsSearch.cpp:189–192 /
                // :222–225, `addBoundExceedingConflict`).
                LpStatus::Cutoff => {
                    if self.config.conflict_cutoff {
                        self.learn_conflict(&node, true);
                    }
                    continue;
                }
                LpStatus::Unbounded if node.depth == 0 => return Ok(TreeEnd::RootUnbounded),
                LpStatus::Unbounded => return Err(MipError::UnboundedNode),
                // Only strong-branching probes stop at an iteration cap.
                LpStatus::IterationLimit => return Err(MipError::Lp(LpError::IterationLimit)),
                LpStatus::Optimal => {}
            }
            if lp.objective > self.upper_limit {
                continue;
            }
            let Some(mut lp) = self.redcost_node(&mut node, lp)? else {
                continue;
            };
            clamp_integers(&mut lp.x, &node.lower, &node.upper, &self.p.integer);
            let mut frac = self.brancher.fractional(&lp.x);
            if frac.is_empty() && self.candidate_violation(&lp.x) > CANDIDATE_TOL {
                // The engine trusts a warm point by a violation relative to
                // max(1, |rhs|); an incumbent needs it absolute, so a big-M
                // row can pass the first and fail the second (jobshop_bigM
                // s35 after the BFRT ratio test, b4a5172). As HiGHS does for
                // an integral point the unscaled LP violates, try it rounded
                // and clamped (HighsLpRelaxation.cpp:1630–1648); failing
                // that, re-solve the node cold before giving up.
                let snapped = snap_integers(lp.x.clone(), &self.p.integer);
                if self.candidate_violation(&snapped) <= CANDIDATE_TOL {
                    lp.x = snapped;
                } else {
                    let cold =
                        match self
                            .engine
                            .solve_from(&node.lower, &node.upper, None, f64::INFINITY)
                        {
                            // The cold re-solve did not finish either — the
                            // node is unsolved, handled like a stalled node LP.
                            Err(LpError::IterationLimit | LpError::Simplex(_))
                                if node.depth > 0 =>
                            {
                                match self.unsolved_node(&node, inh_cols(&node))? {
                                    Unsolved::Branch(kids) => {
                                        let (first, second) = *kids;
                                        tree.push_sibling(second);
                                        plunge = Some(first);
                                        continue;
                                    }
                                    Unsolved::Solved(lp) => lp,
                                    Unsolved::Pruned => continue,
                                }
                            }
                            cold => cold?,
                        };
                    self.lp_iterations += cold.iterations;
                    match cold.status {
                        LpStatus::Infeasible | LpStatus::Cutoff => continue,
                        LpStatus::Optimal if cold.objective > self.upper_limit => continue,
                        LpStatus::Optimal => {}
                        LpStatus::Unbounded => return Err(MipError::UnboundedNode),
                        LpStatus::IterationLimit => {
                            return Err(MipError::Lp(LpError::IterationLimit));
                        }
                    }
                    lp = cold;
                    clamp_integers(&mut lp.x, &node.lower, &node.upper, &self.p.integer);
                    frac = self.brancher.fractional(&lp.x);
                }
            }
            if node.depth == 0 && self.config.audit {
                eprintln!(
                    "audit root_heur t={:.6}",
                    self.start.elapsed().as_secs_f64()
                );
            }
            // Pass 2: re-probe at the end of the root cut loop with the
            // incumbent's cutoff active. Pass 1 (pre-LP, cutoff ∞) found zero
            // fixings on the cover-type instances — their row implications
            // are one-sided, so a direction only becomes infeasible once the
            // objective row `cᵀx ≤ cutoff` propagates against it; this is
            // the pass HiGHS's separation runs
            // (`separateImpliedBounds`, HighsImplications.cpp:583–641), and
            // its fixings land in the global box just in time for the
            // restart check below to measure them. Skipped in racers like
            // pass 1 (one probing pass per solve). An emptied box here means
            // no point beats the cutoff: the incumbent is optimal.
            if node.depth == 0
                && self.config.root_probing
                && self.race.is_none()
                && self.upper_limit.is_finite()
                && !frac.is_empty()
            {
                self.root_probing(self.upper_limit);
                if self.global_infeasible {
                    return Ok(TreeEnd::Exhausted);
                }
            }
            if node.depth == 0 && self.config.audit {
                eprintln!(
                    "audit root_probe t={:.6}",
                    self.start.elapsed().as_secs_f64()
                );
            }
            // Root restart: the root's
            // fixings are final here (propagation + root redcost + node-local
            // redcost all mirrored into the box) and the LP is still
            // fractional — where HiGHS measures `percentageInactiveIntegers`
            // after separation (HighsMipSolverData.cpp:2486–2505).
            // Permanent diagnostic (RUSOLVE_MIP_LOG): why the root restart
            // check does not run when it does not (pruned/integral root).
            if node.depth == 0
                && std::env::var_os("RUSOLVE_MIP_LOG").is_some()
                && !(self.restart_armed
                    && !frac.is_empty()
                    && self.restarts_done < restart::MAX_RESTARTS)
            {
                eprintln!(
                    "restart not checked: armed={} frac={} lp={:.3} restarts_done={}",
                    self.restart_armed,
                    frac.len(),
                    lp.objective,
                    self.restarts_done
                );
            }
            // Never restart past the time limit. HiGHS's root returns
            // on `checkLimits()` between every stage ahead of its restart
            // decision (HighsMipSolverData.cpp:2313–2437 → :2486); ours had
            // no check, so a cut loop stopped by the deadline (co-100) went
            // on to a restart round whose root never solved — reporting −∞
            // over a solved, cut-tightened root. The node loop's
            // `limit_hit` now ends the search on this root's bound.
            if node.depth == 0
                && self.restart_armed
                && !frac.is_empty()
                && self.restarts_done < restart::MAX_RESTARTS
                && !self
                    .config
                    .time_limit
                    .is_some_and(|t| self.start.elapsed().as_secs_f64() >= t)
            {
                match self.restart_request(&lp, &root_cut_log) {
                    Some(req) => return Ok(TreeEnd::Restart(Box::new(req))),
                    // No rate-rule restart — arm the in-tree
                    // stall trigger when the root shows the no-fixings
                    // shape HiGHS serves through clique substitutions.
                    None => self.arm_stall_restart(&lp, &root_cut_log),
                }
            }
            if frac.is_empty() {
                // Trust but verify: an incumbent is only ever an LP point, so
                // a wrong LP answer would become a wrong MIP answer.
                let worst = self.candidate_violation(&lp.x);
                if worst > CANDIDATE_TOL {
                    return Err(MipError::InvalidLpSolution {
                        violation: worst,
                        lower: node.lower,
                        upper: node.upper,
                    });
                }
                let Some(col) = self.unsnappable(&lp.x) else {
                    self.new_incumbent(snap_integers(lp.x, &self.p.integer), Source::Branching);
                    continue;
                };
                // Branch on the near-integral column that a big-M row
                // amplifies, outside strong branching and pseudocosts: its
                // distance to the bound (≤ integrality_tol) would blow up a
                // per-unit cost observation.
                let (z, value) = (lp.objective, lp.x[col]);
                let inh = inh_cols(&node);
                let down = node.child(
                    col,
                    value,
                    false,
                    z,
                    f64::NEG_INFINITY,
                    lp.basis.clone(),
                    inh.clone(),
                );
                let up = node.child(col, value, true, z, f64::NEG_INFINITY, lp.basis, inh);
                tree.push_sibling(up);
                plunge = Some(down);
                continue;
            }
            self.lp_heuristics(&node, &lp.x);
            if self.alns.is_some()
                && from_queue
                && node.depth > 0
                && lp.objective <= self.upper_limit
            {
                let lower = tree.open_bound().min(node.bound);
                self.run_alns(&lp.x, lower);
            } else if self.config.submip
                && from_queue
                && node.depth > 0
                && lp.objective <= self.upper_limit
                && self.more_heuristics_allowed()
                // Cost-aware trigger (ours, not HiGHS): a sub-MIP pays
                // setup the iteration budget cannot see, so a node call
                // needs its wall-time share of the parent tree left to
                // fund it (submip::SUBMIP_TIME_SHARE); the call is then
                // capped at that budget inside `run_submip`.
                && submip::time_share_left(
                    self.submip_secs,
                    self.start.elapsed().as_secs_f64(),
                ) >= submip::SUBMIP_MIN_CALL_SECS
            {
                let x = lp.x.clone();
                self.run_submip(&x, false);
            }
            if lp.objective > self.upper_limit {
                continue;
            }
            if node.depth == 0 {
                self.brancher.set_root_solution(&lp.x);
                if self.config.audit {
                    eprintln!(
                        "audit root_done obj={:e} nodes={} t={:.6}",
                        lp.objective,
                        self.nodes,
                        self.start.elapsed().as_secs_f64()
                    );
                }
            }
            node.estimate = self.brancher.pseudocost.best_estimate(
                lp.objective,
                &frac,
                self.config.feastol,
                self.num_integer,
            );
            self.sb_start = lp.basis.clone();
            let mut brancher = std::mem::take(&mut self.brancher);
            let total = self.lp_iterations;
            let sel = brancher.select(
                self,
                &node.lower,
                &node.upper,
                lp.objective,
                &lp.x,
                &frac,
                total,
            );
            self.brancher = brancher;
            self.sb_start = None;
            let z = lp.objective;
            match sel {
                Selection::Reduce {
                    col,
                    value,
                    up,
                    bound,
                } => {
                    let inh = inh_cols(&node);
                    plunge = Some(node.child(col, value, up, bound.max(z), z, lp.basis, inh));
                }
                Selection::Branch {
                    col,
                    value,
                    down_bound,
                    up_bound,
                    up_first,
                } => {
                    let inh = inh_cols(&node);
                    let down = node.child(
                        col,
                        value,
                        false,
                        down_bound.max(z),
                        z,
                        lp.basis.clone(),
                        inh.clone(),
                    );
                    let up = node.child(col, value, true, up_bound.max(z), z, lp.basis, inh);
                    let (first, second) = if up_first { (up, down) } else { (down, up) };
                    tree.push_sibling(second);
                    plunge = Some(first);
                }
            }
        }
    }

    /// Runs the conflicts `revisit` selects on the node box and
    /// re-propagates what they tightened (HiGHS propagates its conflict pool
    /// inside `HighsDomain::propagate`, HighsDomain.cpp:1692, :1865).
    /// Returns `false` iff the box is infeasible.
    fn conflict_propagate(&mut self, node: &mut Node, revisit: Revisit<'_>) -> bool {
        if !self.config.conflict || self.conflicts.is_empty() {
            return true;
        }
        let full = matches!(revisit, Revisit::All);
        let mut cols: Vec<usize> = match revisit {
            Revisit::All => Vec::new(),
            Revisit::Columns(c) => c.to_vec(),
        };
        for _ in 0..5 {
            let mut changed = Vec::new();
            let r = self.conflicts.propagate(
                &mut node.lower,
                &mut node.upper,
                &self.p.integer,
                self.config.feastol,
                if full {
                    Revisit::All
                } else {
                    Revisit::Columns(&cols)
                },
                &mut changed,
                &mut self.conflict_stats,
            );
            if r.is_err() {
                return false;
            }
            if changed.is_empty() {
                return true;
            }
            changed.sort_unstable();
            changed.dedup();
            let before = if full {
                None
            } else {
                Some((node.lower.clone(), node.upper.clone()))
            };
            if !self.propagator.propagate_cols(
                &mut node.lower,
                &mut node.upper,
                &changed,
                self.upper_limit,
            ) {
                self.conflict_stats.pruned += 1;
                return false;
            }
            // Row propagation can tighten columns beyond the ones the
            // conflicts named; the next round revisits conflicts watching
            // every bound that moved.
            if let Some((lo, hi)) = before {
                cols.clear();
                moved_cols(&lo, &hi, &node.lower, &node.upper, &mut cols);
            }
        }
        true
    }

    /// The node LP is infeasible: turn its Farkas proof into a conflict
    /// (HiGHS `HighsSearch::addInfeasibleConflict`, HighsSearch.cpp:200–235).
    ///
    /// With `cutoff` the proof is the dual-bound one of a node cut off by
    /// the incumbent (`addBoundExceedingConflict`, HighsSearch.cpp:189–192):
    /// `(c − Aᵀy)ᵀx ≤ U − yᵀb`, valid for every solution better than the
    /// cutoff.
    fn learn_conflict(&mut self, node: &Node, cutoff: bool) {
        if !self.config.conflict || node.depth == 0 {
            return;
        }
        let feastol = self.config.feastol;
        let raw = if cutoff {
            self.engine.cutoff_proof()
        } else {
            self.engine.infeasibility_proof()
        };
        let proof = raw.and_then(|(coefs, rhs)| {
            conflict::clean_proof(
                &coefs,
                rhs,
                &self.global_lower,
                &self.global_upper,
                &self.p.integer,
                feastol,
            )
        });
        let Some(proof) = proof else {
            self.conflict_stats.no_proof += 1;
            return;
        };
        let pc = &self.brancher.pseudocost;
        let (up, down) = (
            |j: usize| pc.num_observations_up(j),
            |j: usize| pc.num_observations_down(j),
        );
        let ex = Explain {
            global_lower: &self.global_lower,
            global_upper: &self.global_upper,
            integer: &self.p.integer,
            feastol,
            up_count: &up,
            down_count: &down,
        };
        match conflict::explain_leq(&proof, &node.lower, &node.upper, &ex) {
            Some(lits) if !lits.is_empty() => {
                self.conflict_stats.proofs += 1;
                if cutoff {
                    self.conflict_stats.cutoff_proofs += 1;
                }
                // Credit the conflict to the branching statistics: the
                // decaying weight grows once per explained conflict and is
                // then added to every literal's column on the literal's
                // side — a lower-bound literal credits the up side
                // (HighsDomain.cpp:3833–3841, :3910–3916). Gated by
                // `conflict_score`: measured to hijack branching on
                // degenerate LPs.
                if self.config.conflict_score {
                    let pc = &mut self.brancher.pseudocost;
                    pc.increase_conflict_weight();
                    for l in &lits {
                        pc.add_conflict_observation(l.col as usize, !l.upper);
                    }
                }
                if self.conflicts.add(lits, &self.p.integer, feastol) {
                    self.conflict_stats.added += 1;
                } else {
                    self.conflict_stats.rejected += 1;
                }
            }
            _ => self.conflict_stats.no_proof += 1,
        }
    }

    /// Pseudocost bookkeeping for the node LP just solved (before any node
    /// separation, so the gain is the branching's alone): a cost observation
    /// from the parent objective, or a cutoff observation when the LP is
    /// infeasible or cut off (HiGHS maps the objective limit to infeasible,
    /// HighsLpRelaxation.cpp:1263).
    fn observe_branching(&mut self, node: &Node, lp: &LpOutcome) {
        let Some((side, parent_obj)) = node.side else {
            return;
        };
        let ps = &mut self.brancher.pseudocost;
        match lp.status {
            LpStatus::Infeasible | LpStatus::Cutoff => ps.add_cutoff_observation(side.col, side.up),
            // A parent without an LP objective gives no cost observation.
            LpStatus::Optimal if parent_obj.is_finite() => {
                let target = if side.up {
                    side.point.ceil()
                } else {
                    side.point.floor()
                };
                let objdelta = (lp.objective - parent_obj).max(0.0);
                ps.add_observation(side.col, target - side.point, objdelta);
            }
            _ => {}
        }
    }

    /// HiGHS caps every tree LP at max(10000, 10·avg node LP iterations,
    /// 1.5·first root LP iterations) (HighsMipSolver.cpp:840–843). Without
    /// it a stalled degenerate LP burns the full simplex limit (set_packing,
    /// Seven 100k-iteration cold stalls, 45 s).
    fn tree_lp_iteration_cap(&self) -> u64 {
        let avg = self
            .node_lp_iterations
            .checked_div(self.node_lp_solves)
            .unwrap_or(0);
        (10 * avg).max(3 * self.root_lp_iterations / 2).max(10_000)
    }

    /// A tree node whose LP did not solve (iteration cap, numerical
    /// give-up): branch on an integer column chosen without an LP, bound
    /// inherited (HighsSearch.cpp:1347–1418); with every integer column
    /// fixed, [`Self::resolve_unsolved_fixed`] decides. `inh` is the
    /// children's conflict revisit set.
    fn unsolved_node(&mut self, node: &Node, inh: Vec<usize>) -> Result<Unsolved, MipError> {
        let pick = self.brancher.select_without_lp(
            &node.lower,
            &node.upper,
            &self.p.cost,
            self.incumbent.is_some(),
        );
        let Some((col, point, up_first)) = pick else {
            return Ok(match self.resolve_unsolved_fixed(node)? {
                Some(lp) => Unsolved::Solved(lp),
                None => Unsolved::Pruned,
            });
        };
        self.no_lp_branches += 1;
        let (b, no_lp) = (node.bound, f64::NEG_INFINITY);
        let warm = node.warm.clone();
        let down = node.child(col, point, false, b, no_lp, warm.clone(), inh.clone());
        let up = node.child(col, point, true, b, no_lp, warm, inh);
        Ok(Unsolved::Branch(Box::new(if up_first {
            (up, down)
        } else {
            (down, up)
        })))
    }

    /// An all-integers-fixed node whose LP did not solve, which used to
    /// abort the whole solve. HiGHS re-solves such a node on a fresh LP
    /// (model rows only) with presolve, then primal simplex, then IPM, and
    /// failing all three declares it infeasible with a warning
    /// (HighsSearch.cpp:1420–1476). Ours re-solves it cold and uncapped
    /// with the strategy toggle on failure ([`LpEngine::set_root_lp`]),
    /// cut rows kept (valid inequalities; a cut-free LP would need a second
    /// engine). Still unsolved with time left, it is pruned as in HiGHS and
    /// counted (`LpStats::unsolved_pruned`); past the deadline the error
    /// stands and the search ends as a TimeLimit, so the node's bound
    /// is never dropped from a limit answer.
    fn resolve_unsolved_fixed(&mut self, node: &Node) -> Result<Option<LpOutcome>, MipError> {
        self.engine.set_root_lp(false, true);
        let solved = self
            .engine
            .solve_from(&node.lower, &node.upper, None, self.upper_limit);
        self.engine.set_root_lp(false, false);
        match solved {
            Ok(lp) if lp.status != LpStatus::IterationLimit => Ok(Some(lp)),
            _ if self
                .deadline()
                .is_some_and(|d| std::time::Instant::now() >= d) =>
            {
                Err(MipError::Lp(LpError::IterationLimit))
            }
            _ => {
                self.unsolved_pruned += 1;
                if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
                    eprintln!(
                        "mip: node with all integer columns fixed failed its LP twice; pruned \
                         (HiGHS: declared infeasible)"
                    );
                }
                Ok(None)
            }
        }
    }

    /// The wall-clock end of the search (`start + time_limit`), the same
    /// instant the engine holds.
    fn deadline(&self) -> Option<std::time::Instant> {
        self.config
            .time_limit
            .map(|t| self.start + Duration::from_secs_f64(t))
    }

    fn limit_hit(&self) -> Option<MipStatus> {
        // A stopped racer's answer is discarded (crate::race::merge).
        if self.race.as_ref().is_some_and(Link::stopped) {
            return Some(MipStatus::TimeLimit);
        }
        if self.config.node_limit.is_some_and(|n| self.nodes >= n) {
            return Some(MipStatus::NodeLimit);
        }
        if self
            .config
            .time_limit
            .is_some_and(|t| self.start.elapsed().as_secs_f64() >= t)
        {
            return Some(MipStatus::TimeLimit);
        }
        None
    }

    /// Integer column whose LP value is farthest from integral; ties go to
    /// the lowest index.
    fn most_fractional(&self, x: &[f64], lower: &[f64], upper: &[f64]) -> Option<(usize, f64)> {
        let tol = self.config.integrality_tol;
        let mut best: Option<(usize, f64)> = None;
        let mut best_score = tol;
        for (j, v) in fractional_integers(&self.p.integer, x, lower, upper, tol) {
            let f = v - v.floor();
            let score = f.min(1.0 - f);
            if score > best_score {
                best_score = score;
                best = Some((j, v));
            }
        }
        best
    }

    /// Worst absolute row or bound violation of an incumbent candidate.
    fn candidate_violation(&self, x: &[f64]) -> f64 {
        let v = self.p.violation(x);
        v.row.max(v.bound)
    }

    /// An LP point within `integrality_tol` of integral is not yet a MIP
    /// solution: with big-M rows (M = 1e6) a binary at 1e-6 buys a whole unit
    /// of slack (mip_diff jobshop_bigM_wideM_s34: 214.00002 for a true 215).
    /// Returns the integer column with the largest nonzero rounding error if
    /// rounding the integers breaks feasibility, `None` if the rounded point
    /// is feasible (or nothing rounds).
    fn unsnappable(&self, x: &[f64]) -> Option<usize> {
        let snapped = snap_integers(x.to_vec(), &self.p.integer);
        let v = self.p.violation(&snapped);
        if v.row.max(v.bound) <= CANDIDATE_TOL {
            return None;
        }
        (0..x.len())
            .filter(|&j| self.p.integer[j] && x[j] != snapped[j])
            .max_by(|&a, &b| {
                (x[a] - snapped[a])
                    .abs()
                    .total_cmp(&(x[b] - snapped[b]).abs())
            })
    }

    /// Installs `x` if it improves the incumbent; returns whether it did.
    fn new_incumbent(&mut self, x: Vec<f64>, source: Source) -> bool {
        let z = self.p.objective(&x);
        if self.incumbent.as_ref().is_some_and(|inc| inc.0 <= z) {
            return false;
        }
        let record = FirstIncumbent {
            value: z,
            seconds: self.start.elapsed().as_secs_f64(),
            nodes: self.nodes,
            source,
        };
        if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
            eprintln!(
                "incumbent: {z:.4} at node {}, {:.3}s ({source:?})",
                record.nodes, record.seconds
            );
        }
        if self.config.audit {
            eprintln!(
                "audit incumbent obj={z:e} nodes={} t={:.6} source={source:?}",
                record.nodes, record.seconds
            );
        }
        if self.first_incumbent.is_none() {
            self.first_incumbent = Some(record);
        }
        self.incumbents.push(record);
        if source != Source::Racer
            && let Some(race) = &self.race
        {
            race.publish(z, &x);
        }
        let feastol = self.config.feastol;
        self.upper_limit = if self.integral_objective {
            (z - 0.5).floor() + feastol
        } else {
            (z - feastol).min(z.next_down())
        };
        if let Some(c) = self.config.objective_cutoff {
            self.upper_limit = self.upper_limit.min(self.p.caller_objective(c));
        }
        if let Some(a) = &mut self.alns {
            a.add_solution(z, &x);
        }
        self.incumbent_node = self.nodes;
        self.incumbent = Some((z, x));
        self.apply_root_redcost();
        // The clique-table link in HiGHS's order (extractObjCliques
        // after propagateRootRedcost, HighsMipSolverData.cpp:1576–1591).
        self.apply_clique_fixings();
        if source != Source::LocalMip {
            // `--localmip`: polish the new incumbent (no-op unless the lever
            // is on). The improver's own accepted points do not re-trigger
            // it — its run already continues from each internal improvement.
            self.run_localmip();
        }
        true
    }

    /// Installs the incumbent carried across a restart (its value and point
    /// are already re-referenced into this round's space by
    /// `presolve::solve_presolved`). The upper limit is recomputed from this
    /// problem's objective integrality; the ladder and `first_incumbent` are
    /// restored from the [`Carry`] separately, so nothing is recorded here.
    fn install_incumbent(&mut self, z: f64, x: Vec<f64>) {
        let feastol = self.config.feastol;
        self.upper_limit = if self.integral_objective {
            (z - 0.5).floor() + feastol
        } else {
            (z - feastol).min(z.next_down())
        };
        if let Some(c) = self.config.objective_cutoff {
            self.upper_limit = self.upper_limit.min(self.p.caller_objective(c));
        }
        self.incumbent_node = self.nodes;
        self.incumbent = Some((z, x));
        // The restart carry is this round's first incumbent: run the clique
        // link on it too (HiGHS re-runs extractObjCliques per incumbent, and
        // its presolve runs it once on the carried upper limit,
        // HPresolve.cpp:1494–1512).
        self.apply_clique_fixings();
    }

    /// HiGHS's post-separation fixing-rate trigger: `fixingRate >= 2.5`, or
    /// any fixing at all before the first restart
    /// (HighsMipSolverData.cpp:2486–2505), on the global box (see
    /// `restart_request` for why global).
    fn restart_rate_triggered(
        &self,
        global_lower: &[f64],
        global_upper: &[f64],
        integer: &[bool],
    ) -> bool {
        let fixed = restart::inactive_integers(
            global_lower,
            global_upper,
            &self.p.col_lower,
            &self.p.col_upper,
            integer,
        );
        let rate = 100.0 * fixed as f64 / self.num_integer.max(1) as f64;
        rate >= 2.5 || (rate > 0.0 && self.restarts_done == 0)
    }

    /// Relative distance of the incumbent from the root LP's dual bound —
    /// our proxy for "the search still has real work" (see
    /// `restart::MIN_ROOT_GAP`); infinite while no incumbent exists.
    fn root_gap(&self, lp_objective: f64) -> f64 {
        if self.upper_limit.is_finite() {
            (self.upper_limit - lp_objective).abs().max(0.0) / self.upper_limit.abs().max(1.0)
        } else {
            f64::INFINITY
        }
    }

    /// [`Self::root_gap`] against the stashed cut-loop-end incumbent limit
    /// (see [`Self::restart_arm_limit`]): the bar the restart check and
    /// stall-restart arming measure, so a point the root sub-MIP/dive published
    /// between the cut loop and the check arms the restart instead of
    /// disarming it.
    fn restart_arm_gap(&self, lp_objective: f64) -> f64 {
        if self.restart_arm_limit.is_finite() {
            (self.restart_arm_limit - lp_objective).abs().max(0.0)
                / self.restart_arm_limit.abs().max(1.0)
        } else {
            f64::INFINITY
        }
    }

    /// The restart request if the root's fixing rate passes HiGHS's
    /// post-separation thresholds (`fixingRate >= 2.5`, or any fixing at all
    /// before the first restart, HighsMipSolverData.cpp:2486–2505), counting
    /// integers the root fixed beyond the model's own bounds,
    /// and the root gap gate
    /// (`restart::MIN_ROOT_GAP`, our design) says a restart can pay.
    ///
    /// The box is the *global* one, matching HiGHS measuring
    /// `percentageInactiveIntegers` on the solver's domain, not a node box:
    /// the root's propagation fixings live there via the sync in
    /// `root_incumbent`/`root_after_cutoff`, and the root reduced-cost
    /// fixings applied after the last sync (`apply_root_redcost`, and
    /// `redcost_node`'s push-up at :2125) only ever reach the global box.
    /// At depth 0 the node box never goes below it, so global is the tightest
    /// root domain either way.
    fn restart_request(
        &self,
        lp: &LpOutcome,
        cut_log: &[(cuts::CutClass, CutRow)],
    ) -> Option<RestartRequest> {
        let fixed = restart::inactive_integers(
            &self.global_lower,
            &self.global_upper,
            &self.p.col_lower,
            &self.p.col_upper,
            &self.p.integer,
        );
        let rate = 100.0 * fixed as f64 / self.num_integer.max(1) as f64;
        // Against the cut-loop-end limit (`restart_arm_gap`), so the
        // root sub-MIP/dive's own incumbent cannot close this bar.
        let root_gap = self.restart_arm_gap(lp.objective);
        if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
            eprintln!(
                "restart check: {rate:.1}% inactive integers ({fixed} of {}), {} cuts tight, incumbent {}, root gap {:.4}",
                self.num_integer,
                restart::active_cuts(cut_log, &lp.x, self.config.feastol).len(),
                if self.upper_limit.is_finite() {
                    format!("{:.4}", self.upper_limit)
                } else {
                    "none".to_string()
                },
                root_gap
            );
        }
        if !self.restart_rate_triggered(&self.global_lower, &self.global_upper, &self.p.integer) {
            return None;
        }
        if root_gap < self.config.restart_min_gap {
            return None;
        }
        Some(RestartRequest {
            lower: self.global_lower.clone(),
            upper: self.global_upper.clone(),
            cuts: restart::active_cuts(cut_log, &lp.x, self.config.feastol),
            incumbent: self.incumbent.clone(),
            bound: lp.objective,
            carry: Carry {
                // The presolve loop forward-maps the incumbent through the
                // new presolve before installing it.
                incumbent: None,
                nodes: self.nodes,
                lp_iterations: self.lp_iterations,
                heur_lp_iterations: self.heur_lp_iterations,
                restarts_done: self.restarts_done + 1,
                first_incumbent: self.first_incumbent,
                incumbents: self.incumbents.clone(),
            },
        })
    }

    /// Stall-restart arming: the root check found the shape HiGHS's fixing-rate
    /// rule cannot serve (no fixings — we have no clique table) but a restart
    /// could still pay in: the gap is still [`restart::MIN_ROOT_GAP`] wide
    /// over a fractional root LP and the loop left tight cut rows the
    /// re-presolve can deduce fixings from (pg: gap 0.31, 2 473 tight cuts).
    /// The decision itself moves into the tree
    /// ([`StallRestart::verdict`]): root-only inputs provably cannot separate this
    /// shape from trees that finish clean (int_production s1 / jobshop s1
    /// measure the same 0 % / gap ≥ 0.10 / cuts ≥ 1 on our benchmark set), so the
    /// restart waits for a demonstrated stall. First restart only: a stall
    /// round that produces no fixings must not repeat (our analog of HiGHS's
    /// `restart_presolve_reduction_limit ≥ 1` cycle guard,
    /// HighsMipSolverData.cpp:1390–1402).
    fn arm_stall_restart(&mut self, lp: &LpOutcome, cut_log: &[(cuts::CutClass, CutRow)]) {
        if self.restarts_done != 0 {
            return;
        }
        let tight = restart::active_cuts(cut_log, &lp.x, self.config.feastol).len();
        // Against the cut-loop-end limit (`restart_arm_gap`), same as
        // the rate rule above — the dive's point must arm, not disarm.
        let gap = self.restart_arm_gap(lp.objective);
        // Cut-mass bar (`restart::STALL_RESTART_CUT_DENSITY`): without enough tight
        // rows per model row the re-presolve reduces nothing and the restart
        // is pure overhead (enlight_hard measured 2.1× worse, see the const).
        let mass = (restart::STALL_RESTART_CUT_DENSITY * self.p.num_row as f64).max(1.0);
        if (tight as f64) < mass || gap < self.config.restart_min_gap {
            if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
                eprintln!(
                    "stall restart not armed: {tight} cuts tight (< {mass:.0} wanted), root gap {gap:.4}"
                );
            }
            return;
        }
        if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
            eprintln!(
                "stall restart armed: {tight} cuts tight, root gap {gap:.4}, incumbent {:?}",
                self.incumbent.as_ref().map(|i| i.0)
            );
        }
        self.stall_restart = Some(StallRestart {
            log: cut_log.to_vec(),
            x: lp.x.clone(),
            lp_obj: lp.objective,
            watch: self.incumbent.as_ref().map(|i| i.0),
            since_nodes: self.nodes,
            since_secs: self.start.elapsed().as_secs_f64(),
        });
    }

    /// The stall-restart request, built from the stashed root context when
    /// the tree stalled. Same contract as [`Self::restart_request`]; the box
    /// is the current global one (tighter than the root's, still valid), the
    /// incumbent is whatever the tree has reached, and the carried dual
    /// bound is the tree's own open bound when it passed the root LP's.
    fn stall_restart_request(&self, lb: &StallRestart, bound: f64) -> RestartRequest {
        RestartRequest {
            lower: self.global_lower.clone(),
            upper: self.global_upper.clone(),
            cuts: restart::active_cuts(&lb.log, &lb.x, self.config.feastol),
            incumbent: self.incumbent.clone(),
            bound: lb.lp_obj.max(bound),
            carry: Carry {
                incumbent: None,
                nodes: self.nodes,
                lp_iterations: self.lp_iterations,
                heur_lp_iterations: self.heur_lp_iterations,
                restarts_done: self.restarts_done + 1,
                first_incumbent: self.first_incumbent,
                incumbents: self.incumbents.clone(),
            },
        }
    }

    /// Intersects the node box with the global box; returns the columns
    /// whose bounds moved, or `None` if the box became empty.
    fn intersect_global(&self, node: &mut Node) -> Option<Vec<usize>> {
        let mut moved = Vec::new();
        for j in 0..node.lower.len() {
            let (lo, up) = (self.global_lower[j], self.global_upper[j]);
            let mut hit = false;
            if lo > node.lower[j] {
                node.lower[j] = lo;
                hit = true;
            }
            if up < node.upper[j] {
                node.upper[j] = up;
                hit = true;
            }
            if hit {
                if node.lower[j] > node.upper[j] {
                    return None;
                }
                moved.push(j);
            }
        }
        Some(moved)
    }

    /// Stores the reduced costs of a root LP solved over the global box and,
    /// with an upper limit, applies them at once (HiGHS calls
    /// `addRootRedcost` then `propagateRootRedcost` after every root LP,
    /// HighsMipSolverData.cpp:1946–1953, HighsSeparation.cpp:80–85, :169–175).
    fn add_root_redcost(&mut self, node: &Node, lp: &LpOutcome) {
        if !self.config.redcost || lp.status != LpStatus::Optimal {
            return;
        }
        let Some(duals) = self.engine.duals(lp) else {
            return;
        };
        self.root_redcost.add(
            &duals,
            lp.objective,
            &node.lower,
            &node.upper,
            &self.p.integer,
            self.config.feastol,
        );
        self.apply_root_redcost();
    }

    /// Tightens the global box from the stored root reduced costs for the
    /// current upper limit, then propagates it (HiGHS new-incumbent path,
    /// HighsMipSolverData.cpp:1576–1578). An empty box ends the search.
    fn apply_root_redcost(&mut self) {
        if !self.config.redcost || self.global_infeasible {
            return;
        }
        let mut changed = Vec::new();
        let applied = self.root_redcost.apply(
            self.upper_limit,
            &mut self.global_lower,
            &mut self.global_upper,
            self.config.feastol,
            &mut changed,
        );
        match applied {
            Ok(0) => {}
            Ok(n) => {
                self.redcost_stats.root_tightenings += n;
                if !self.propagator.propagate_cols(
                    &mut self.global_lower,
                    &mut self.global_upper,
                    &changed,
                    self.upper_limit,
                ) {
                    self.global_infeasible = true;
                }
            }
            Err(redcost::Infeasible) => self.global_infeasible = true,
        }
    }

    /// Root probing over the global box, in the two HiGHS root-scope
    /// passes. Pass 1 runs before the root LP with no cutoff (`HPresolve::
    /// runProbing` inside MIP presolve); pass 2 runs at the end of the root
    /// cut loop with the incumbent's cutoff — the probing HiGHS's separation
    /// runs under the domain cutoff (`separateImpliedBounds`,
    /// HighsImplications.cpp:583–641) — where a direction infeasible under
    /// `cᵀx ≤ cutoff` fixes the binary toward the incumbent's side; those
    /// fixings are what the restart fixings-% check measures and what
    /// seeds the clique table. Tightens `global_lower`/`global_upper` in
    /// place (the root node inherits the box through `intersect_global`);
    /// `stats.infeasible` empties the box — with a finite cutoff that proves
    /// the incumbent optimal (the caller ends the search on it), with ∞ the
    /// problem is infeasible. MIT — see `crate::probing` for the full port
    /// inventory.
    fn root_probing(&mut self, cutoff: f64) {
        let stats = crate::probing::run(
            self.p,
            &mut self.propagator,
            &self.locks,
            self.cliques.as_mut(),
            &mut self.global_lower,
            &mut self.global_upper,
            crate::probing::Params {
                feastol: self.config.feastol,
                deadline: self
                    .config
                    .time_limit
                    .map(|t| self.start + std::time::Duration::from_secs_f64(t)),
                contingent: crate::probing::CONTINGENT,
                cutoff,
            },
        );
        if stats.infeasible {
            self.global_infeasible = true;
        }
        if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
            let fixed_pct = 100.0
                * restart::inactive_integers(
                    &self.global_lower,
                    &self.global_upper,
                    &self.p.col_lower,
                    &self.p.col_upper,
                    &self.p.integer,
                ) as f64
                / self.num_integer.max(1) as f64;
            let cutoff = if cutoff.is_finite() {
                format!("{cutoff:.4}")
            } else {
                "none".to_string()
            };
            eprintln!(
                "root probing (cutoff {cutoff}): {} of {} candidates probed, \
                 {} infeasible directions, {} fixings ({fixed_pct:.1}% of integers), \
                 {} tightenings, {} cliques, {:.3}s",
                stats.probed,
                stats.candidates,
                stats.infeasible_probes,
                stats.fixings,
                stats.tightenings,
                stats.cliques,
                stats.seconds
            );
        }
        self.probing_stats = self.probing_stats.merge(stats);
    }

    /// The incumbent → fixings link, in HiGHS's new-incumbent order
    /// (`addSolutionImprovingSolution`, HighsMipSolverData.cpp:1576–1591:
    /// propagate → `propagateRootRedcost` → `extractObjCliques` — this runs
    /// right after `apply_root_redcost`). Re-extracts cliques from the
    /// objective cutoff row and propagates every clique into the global box,
    /// interleaved with row propagation to a fixpoint; an empty box proves
    /// the incumbent optimal and ends the search. HiGHS runs this on every
    /// improving incumbent wherever it was found (its domain is global, and
    /// the tightened bounds reach open nodes through `performBounding`);
    /// the table itself stays root-built (HiGHS also adds cliques
    /// in-tree from cuts and probing).
    fn apply_clique_fixings(&mut self) {
        if self.global_infeasible || self.cliques.is_none() {
            return;
        }
        let Some(table) = self.cliques.as_mut() else {
            return;
        };
        let feastol = self.config.feastol;
        let added = match table.extract_objective_cliques(
            &self.p.cost,
            &self.global_lower,
            &self.global_upper,
            &self.p.integer,
            self.upper_limit,
            feastol,
        ) {
            cuts::ObjCliques::CutoffInfeasible => {
                self.clique_stats.cutoff_pruned = true;
                self.global_infeasible = true;
                if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
                    eprintln!("clique cutoff: box minimum above the upper limit");
                }
                return;
            }
            cuts::ObjCliques::Added(k) => k,
        };
        self.clique_stats.obj_cliques += added as u64;
        let mut round_fixings = 0u64;
        loop {
            let mut changed = Vec::new();
            let table = self.cliques.as_ref().unwrap();
            if !table.propagate_fixings(
                &mut self.global_lower,
                &mut self.global_upper,
                feastol,
                &mut changed,
            ) {
                self.global_infeasible = true;
                self.clique_stats.fixings += round_fixings;
                return;
            }
            if changed.is_empty() {
                break;
            }
            round_fixings += changed.len() as u64;
            if !self.propagator.propagate_cols(
                &mut self.global_lower,
                &mut self.global_upper,
                &changed,
                self.upper_limit,
            ) {
                self.global_infeasible = true;
                break;
            }
        }
        self.clique_stats.fixings += round_fixings;
        if (added > 0 || round_fixings > 0) && std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
            let table = self.cliques.as_ref().unwrap();
            eprintln!(
                "clique: +{added} obj cliques, {round_fixings} fixings ({} cliques, {} entries)",
                table.num_cliques(),
                table.num_entries()
            );
        }
    }

    /// Node-local reduced-cost fixing after the node LP `lp` (HiGHS
    /// `propagateRedCost` + `localdom.propagate()`, HighsSearch.cpp:1000–1029):
    /// tightens integer bounds of the node box, propagates them, and when
    /// the box no longer holds the LP point re-solves the LP warm from its
    /// basis and repeats (HiGHS re-evaluates the node). At the root the
    /// tightenings are global. `None` = the node is pruned.
    fn redcost_node(
        &mut self,
        node: &mut Node,
        lp: LpOutcome,
    ) -> Result<Option<LpOutcome>, MipError> {
        /// Re-solves per node; HiGHS repeats until nothing moves.
        const MAX_RESOLVES: usize = 2;
        if !self.config.redcost
            || self.incumbent.is_none()
            || lp.status != LpStatus::Optimal
            || !self.upper_limit.is_finite()
        {
            return Ok(Some(lp));
        }
        let mut lp = lp;
        let mut counted = false;
        let mut tightened = false;
        for round in 0..=MAX_RESOLVES {
            let Some(duals) = self.engine.duals(&lp) else {
                break;
            };
            if !counted {
                self.redcost_stats.nodes_with_duals += 1;
                counted = true;
            }
            let (lower0, upper0) = (node.lower.clone(), node.upper.clone());
            let mut changed = Vec::new();
            let fixed = redcost::fix_local(
                &duals,
                lp.objective,
                self.upper_limit,
                &self.p.integer,
                &mut node.lower,
                &mut node.upper,
                self.config.feastol,
                &mut changed,
            );
            self.redcost_stats.local_tightenings += changed.len() as u64;
            if !changed.is_empty() && !tightened {
                tightened = true;
                self.redcost_stats.nodes_tightened += 1;
            }
            if fixed.is_err()
                || !self.propagator.propagate_cols(
                    &mut node.lower,
                    &mut node.upper,
                    &changed,
                    self.upper_limit,
                )
            {
                self.redcost_stats.pruned += 1;
                return Ok(None);
            }
            if changed.is_empty() {
                break;
            }
            if node.depth == 0 {
                for j in 0..node.lower.len() {
                    self.global_lower[j] = self.global_lower[j].max(node.lower[j]);
                    self.global_upper[j] = self.global_upper[j].min(node.upper[j]);
                }
            }
            let inside = lp.x.iter().enumerate().all(|(j, &v)| {
                let tol = 1e-9 * v.abs().max(1.0);
                v >= node.lower[j] - tol && v <= node.upper[j] + tol
            });
            if inside {
                break;
            }
            if round == MAX_RESOLVES {
                // The tree code needs an LP point inside the node box (an
                // integral one is verified against the rows).
                node.lower = lower0;
                node.upper = upper0;
                break;
            }
            // Propagation moved the box off the LP point: re-solve.
            self.redcost_stats.resolves += 1;
            if node.depth > 0 {
                let cap = self.tree_lp_iteration_cap();
                self.engine.set_iteration_limit(Some(cap));
            }
            // The root node's cold LPs may take the dual, and with a
            // time limit a root LP that fails (iteration cap, numerical
            // give-up) re-runs once on the other simplex with no cap (HiGHS:
            // no root iteration limit, HighsMipSolverData.cpp:2018; strategy
            // toggle on failure, HighsLpRelaxation.cpp:1236-1251). A root
            // still unfinished at the deadline ends as a clean TimeLimit.
            self.engine.set_root_lp(
                node.depth == 0 && self.config.root_lp_dual,
                node.depth == 0 && self.config.time_limit.is_some(),
            );
            let solved = self.engine.solve_from(
                &node.lower,
                &node.upper,
                lp.basis.as_deref(),
                self.upper_limit,
            );
            self.engine.set_iteration_limit(None);
            match solved {
                Ok(out) => {
                    self.lp_iterations += out.iterations;
                    match out.status {
                        LpStatus::Infeasible | LpStatus::Cutoff => {
                            self.redcost_stats.pruned += 1;
                            return Ok(None);
                        }
                        LpStatus::Optimal if out.objective > self.upper_limit => {
                            self.redcost_stats.pruned += 1;
                            return Ok(None);
                        }
                        LpStatus::Optimal => lp = out,
                        _ => {
                            node.lower = lower0;
                            node.upper = upper0;
                            break;
                        }
                    }
                }
                // The tightenings stay valid; without an LP over them the
                // node keeps its previous box and LP.
                Err(e) => {
                    if e == LpError::IterationLimit {
                        self.lp_iterations += self.tree_lp_iteration_cap();
                    }
                    node.lower = lower0;
                    node.upper = upper0;
                    break;
                }
            }
        }
        Ok(Some(lp))
    }
}

/// Rounds the integer entries of an LP point to exact integers, as HiGHS does
/// before its repair LP (HighsMipSolverData.cpp:1194, MIT).
fn snap_integers(mut x: Vec<f64>, integer: &[bool]) -> Vec<f64> {
    for (xj, &int) in x.iter_mut().zip(integer) {
        if int {
            *xj = xj.round();
        }
    }
    x
}

/// Clamps the integer entries of an LP point into the node box, the same
/// clamp `Search::most_fractional` applies (HighsLpRelaxation.cpp:1467–1472),
/// done in place so that branching candidates, pseudocost fractional parts
/// and strong-branching points all see the clamped values.
fn clamp_integers(x: &mut [f64], lower: &[f64], upper: &[f64], integer: &[bool]) {
    for (j, xj) in x.iter_mut().enumerate() {
        if integer[j] {
            *xj = xj.min(upper[j]).max(lower[j]);
        }
    }
}

/// Strong branching runs on the search's own LP engine and propagator, and
/// hands integral LP points to the incumbent logic.
impl<E: LpEngine> SbHost for Search<'_, E> {
    /// A warm dual probe from the branched node's basis under the cutoff,
    /// capped at `sb_iteration_limit` iterations (HiGHS strong-branches by
    /// dual simplex from the node basis under the objective limit,
    /// HighsSearch.cpp:405, :574). A cold fallback (no basis, rejected
    /// snapshot) is capped like any tree LP.
    fn sb_solve(&mut self, lower: &[f64], upper: &[f64]) -> Result<LpOutcome, LpError> {
        let cap = self.tree_lp_iteration_cap();
        let limit = self.sb_iteration_limit();
        self.engine.set_iteration_limit(Some(cap));
        let start = self.sb_start.clone();
        let lp = self
            .engine
            .probe(lower, upper, start.as_deref(), self.upper_limit, limit);
        self.engine.set_iteration_limit(None);
        let mut lp = lp;
        match &mut lp {
            Ok(out) => {
                self.lp_iterations += out.iterations;
                if out.status == LpStatus::Optimal {
                    clamp_integers(&mut out.x, lower, upper, &self.p.integer);
                }
            }
            Err(LpError::IterationLimit) => self.lp_iterations += cap,
            Err(_) => {}
        }
        lp
    }

    fn sb_iteration_limit(&self) -> u64 {
        let avg = self
            .node_lp_iterations
            .checked_div(self.node_lp_solves)
            .unwrap_or(0);
        self.config.branching.sb_iteration_limit.resolve(
            avg,
            self.root_lp_iterations,
            self.tree_lp_iteration_cap(),
        )
    }

    fn sb_propagate(&mut self, lower: &mut [f64], upper: &mut [f64], col: usize) -> bool {
        self.propagator
            .propagate(lower, upper, Some(col), self.upper_limit)
    }

    fn upper_limit(&self) -> f64 {
        self.upper_limit
    }

    fn sb_integral(&mut self, x: Vec<f64>) {
        let v = self.p.violation(&x);
        if v.row.max(v.bound) <= CANDIDATE_TOL && self.unsnappable(&x).is_none() {
            self.new_incumbent(snap_integers(x, &self.p.integer), Source::Branching);
        }
    }

    fn sb_stop(&self) -> bool {
        self.race.as_ref().is_some_and(Link::stopped)
            || self
                .config
                .time_limit
                .is_some_and(|t| self.start.elapsed().as_secs_f64() >= t)
    }

    fn sb_nodes(&self) -> u64 {
        self.nodes
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_core::CscMatrix;
    use rusolve_io::read_mps;
    use rusolve_lp::rng::XorShift64;

    fn stall_restart(watch: Option<f64>, since_nodes: u64, since_secs: f64) -> StallRestart {
        StallRestart {
            log: Vec::new(),
            x: Vec::new(),
            lp_obj: 0.0,
            watch,
            since_nodes,
            since_secs,
        }
    }

    /// Stall-restart window logic (`StallRestart::verdict`): a significant
    /// incumbent improvement (≥ `restart::STALL_IMPROVEMENT` relative)
    /// re-opens the window; a crawl below it does not; firing needs both
    /// stall bars over the *current* window.
    #[test]
    fn stall_restart_verdict_rearms_only_on_significant_improvement() {
        // No incumbent ever and the bars passed: the cold shape fires too
        // (enlight_hard armed incumbent-less; a tightened re-presolve helps
        // a cold search most — restart.rs MIN_ROOT_GAP docs).
        let mut lb = stall_restart(None, 1, 0.0);
        assert!(matches!(
            lb.verdict(None, 500, 5.0, 300, 2.0),
            StallRestartVerdict::Fire
        ));
        // First sighting: significant → re-arm at (500, 5.0).
        assert!(matches!(
            lb.verdict(Some(100.0), 500, 5.0, 300, 2.0),
            StallRestartVerdict::Wait
        ));
        // Bars unmet: wait.
        assert!(matches!(
            lb.verdict(Some(100.0), 600, 5.5, 300, 2.0),
            StallRestartVerdict::Wait
        ));
        // Both bars met, no movement: fire.
        assert!(matches!(
            lb.verdict(Some(100.0), 800, 7.1, 300, 2.0),
            StallRestartVerdict::Fire
        ));
        // A 0.05 % crawl does not re-arm: same window still fires.
        let mut lb = stall_restart(Some(100.0), 500, 5.0);
        assert!(matches!(
            lb.verdict(Some(99.95), 800, 7.1, 300, 2.0),
            StallRestartVerdict::Fire
        ));
        // A 0.5 % improvement is significant: re-arm at (800, 7.1) — the
        // node bar now needs 1100 again, not 800.
        assert!(matches!(
            lb.verdict(Some(99.5), 800, 7.1, 300, 2.0),
            StallRestartVerdict::Wait
        ));
        assert!(matches!(
            lb.verdict(Some(99.5), 1000, 9.0, 300, 2.0),
            StallRestartVerdict::Wait
        ));
        assert!(matches!(
            lb.verdict(Some(99.5), 1100, 9.2, 300, 2.0),
            StallRestartVerdict::Fire
        ));
        // Incumbent vanished (cannot happen in-tree; carried restarts
        // re-install before the search runs): the window is untouched, so
        // an already-stalled window keeps its bars.
        let mut lb = stall_restart(Some(100.0), 500, 5.0);
        assert!(matches!(
            lb.verdict(None, 800, 7.1, 300, 2.0),
            StallRestartVerdict::Fire
        ));
    }

    /// The restart bar measures the gap against the
    /// cut-loop-end incumbent limit, not the current one — a point the root
    /// sub-MIP/dive published between the cut loop and the check must not
    /// close the bar (`restart_arm_gap` vs `root_gap`, mik's 0.0395-vs-
    /// 0.1629 shape at lp −54 331).
    #[test]
    fn restart_arm_gap_ignores_post_cut_loop_incumbents() {
        let inf = f64::INFINITY;
        let p = MipProblem {
            num_col: 2,
            num_row: 1,
            matrix: CscMatrix::from_triplets(1, 2, &[(0, 0, 2.0), (0, 1, 3.0)]),
            row_lower: vec![-inf],
            row_upper: vec![4.5],
            cost: vec![-5.0, -4.0],
            maximize: false,
            col_lower: vec![0.0; 2],
            col_upper: vec![1.0; 2],
            integer: vec![true; 2],
        };
        let config = MipConfig::default();
        let mut search = Search::new(&p, engine(&p, &config), &config, Instant::now());
        // Cut loop ended at −46 719; the dive then published −52 266; the
        // root LP bound is −54 331 (mik's shape, scaled signs unchanged).
        search.restart_arm_limit = -46_719.0;
        search.upper_limit = -52_266.0;
        let lp = -54_331.0;
        assert!(search.root_gap(lp) < restart::MIN_ROOT_GAP);
        assert!(search.restart_arm_gap(lp) >= restart::MIN_ROOT_GAP);
        // Stash unset (cold root, no incumbent): both infinite, never block.
        search.restart_arm_limit = f64::INFINITY;
        search.upper_limit = f64::INFINITY;
        assert_eq!(search.restart_arm_gap(0.0), f64::INFINITY);
    }

    fn solve_mps(text: &str) -> MipSolution {
        let p = MipProblem::from_model(&read_mps(text).unwrap());
        let sol = solve(&p, &MipConfig::default()).unwrap();
        if !sol.x.is_empty() {
            let v = p.violation(&sol.x);
            assert!(
                v.bound <= 1e-6 && v.row <= 1e-6 && v.integrality <= 1e-6,
                "{v:?}"
            );
        }
        sol
    }

    /// E2e round-trip through the presolve chain with restart armed:
    /// the answer with `restart: true`
    /// must equal the `restart: false` answer. (Whether the trigger FIRES is
    /// a property of root reduction strength, verified by the gate's
    /// per-instance restart counts — mik reaches 23.6% and restarts; every
    /// small hand-crafted instance here is solved to proof by root cutoff
    /// propagation before the check, so no synthetic needle is kept: the
    /// trigger's counting/threshold logic is unit-tested in
    /// `restart::tests::inactive_integers_counts_only_root_fixings`.)
    /// Eight binaries in one knapsack `2a+3b+4c+5d+6e+7f + 8g+9h ≤ 17`,
    /// strictly decreasing profit/weight; distinct weights so presolve keeps
    /// every column (identical ones collapse as parallel, an all-even row
    /// gcd-rounds to an integral root). Root cuts off: lifted covers close
    /// the root gap.
    #[test]
    fn root_restart_roundtrip_postsolves_and_carries() {
        let mut mps = String::from("NAME RSTART\nROWS\n N  OBJ\n L  CAP\nCOLUMNS\n");
        for (j, (cost, w)) in [
            (-10.0, 2),
            (-14.0, 3),
            (-17.0, 4),
            (-20.0, 5),
            (-21.0, 6),
            (-21.0, 7),
            (-4.0, 8),
            (-4.0, 9),
        ]
        .iter()
        .enumerate()
        {
            mps.push_str(&format!(" X{j} OBJ {cost}  CAP {w}\n"));
        }
        mps.push_str("RHS\n    RHS CAP 17\nBOUNDS\n");
        for j in 0..8 {
            mps.push_str(&format!(" BV BND X{j}\n"));
        }
        mps.push_str("ENDATA\n");
        let p = MipProblem::from_model(&read_mps(&mps).unwrap());
        let base = MipConfig {
            cuts: crate::cuts::CutConfig {
                enabled: false,
                ..Default::default()
            },
            ..MipConfig::default()
        };
        let off = solve(
            &p,
            &MipConfig {
                restart: false,
                ..base
            },
        )
        .unwrap();
        let sol = solve(&p, &base).unwrap();
        assert_eq!(sol.status, MipStatus::Optimal);
        assert!(
            (sol.objective - off.objective).abs() < 1e-6,
            "{} vs restart-off {}",
            sol.objective,
            off.objective
        );
    }

    #[test]
    fn binary_knapsack() {
        // max 10a + 13b + 7c + 8d  s.t.  4a + 6b + 3c + 5d ≤ 10 ⇒ a,b (23)
        let sol = solve_mps(
            "\
NAME KNAP
OBJSENSE
 MAX
ROWS
 N  OBJ
 L  CAP
COLUMNS
    M  'MARKER'  'INTORG'
    A  OBJ  10  CAP  4
    B  OBJ  13  CAP  6
    C  OBJ  7  CAP  3
    D  OBJ  8  CAP  5
    M  'MARKER'  'INTEND'
RHS
    RHS  CAP  10
ENDATA
",
        );
        assert_eq!(sol.status, MipStatus::Optimal);
        assert!((sol.objective - 23.0).abs() < 1e-9, "{}", sol.objective);
        assert!((sol.bound - 23.0).abs() < 1e-9);
    }

    #[test]
    fn parity_row_is_infeasible() {
        // 2x − 2y = 1 has no integer solution.
        let sol = solve_mps(
            "\
NAME PAR
ROWS
 N  OBJ
 E  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  OBJ  1  R  2
    Y  OBJ  1  R  -2
    M  'MARKER'  'INTEND'
RHS
    RHS  R  1
BOUNDS
 UP BND X 10
 UP BND Y 10
ENDATA
",
        );
        assert_eq!(sol.status, MipStatus::Infeasible);
    }

    #[test]
    fn unbounded_relaxation_infeasible_mip() {
        // min −z: relaxation unbounded along continuous z; the integer
        // part 2x − 2y = 1 over a finite box has no solution ⇒ Infeasible,
        // not Unbounded (Meyer 1974 path). Unbounded *integer* boxes here
        // need gcd reasoning (presolve): plain B&B never terminates.
        let sol = solve_mps(
            "\
NAME UNB_INF
ROWS
 N  OBJ
 E  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  R  2
    Y  R  -2
    M  'MARKER'  'INTEND'
    Z  OBJ  -1
RHS
    RHS  R  1
BOUNDS
 UP BND X 6
 UP BND Y 6
ENDATA
",
        );
        assert_eq!(sol.status, MipStatus::Infeasible);
    }

    #[test]
    fn unbounded_mip() {
        // max x, x − 2y ≤ 0, x, y ∈ ℤ₊ ⇒ unbounded.
        let sol = solve_mps(
            "\
NAME UNB
OBJSENSE
 MAX
ROWS
 N  OBJ
 L  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  OBJ  1  R  1
    Y  R  -2
    M  'MARKER'  'INTEND'
RHS
    RHS  R  0
BOUNDS
 LI BND X 0
 LI BND Y 0
ENDATA
",
        );
        assert_eq!(sol.status, MipStatus::Unbounded);
        assert_eq!(sol.objective, f64::INFINITY);
    }

    #[test]
    fn limits_inside_unbounded_relaxation_path_never_report_unbounded() {
        // Both instances take the RootUnbounded → feasibility-B&B path; a
        // node or time limit hit there must surface as that limit.
        let unb_inf = "\
NAME UNB_INF
ROWS
 N  OBJ
 E  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  R  2
    Y  R  -2
    M  'MARKER'  'INTEND'
    Z  OBJ  -1
RHS
    RHS  R  1
BOUNDS
 UP BND X 6
 UP BND Y 6
ENDATA
";
        let unb = "\
NAME UNB
OBJSENSE
 MAX
ROWS
 N  OBJ
 L  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  OBJ  1  R  1
    Y  R  -2
    M  'MARKER'  'INTEND'
RHS
    RHS  R  0
BOUNDS
 LI BND X 0
 LI BND Y 0
ENDATA
";
        let limited = |text: &str, node_limit, time_limit, heuristics| {
            let p = MipProblem::from_model(&read_mps(text).unwrap());
            let cfg = MipConfig {
                node_limit,
                time_limit,
                heuristics,
                ..MipConfig::default()
            };
            (p.clone(), solve(&p, &cfg).unwrap())
        };
        for heur in [HeuristicConfig::none(), HeuristicConfig::default()] {
            for (n, t) in [
                (Some(1), None),
                (Some(2), None),
                (Some(3), None),
                (None, Some(0.0)),
            ] {
                let (_, sol) = limited(unb_inf, n, t, heur);
                assert_ne!(sol.status, MipStatus::Unbounded, "UNB_INF {n:?} {t:?}");
                let (p, sol) = limited(unb, n, t, heur);
                if sol.status == MipStatus::Unbounded {
                    let v = p.violation(&sol.x);
                    assert!(v.row.max(v.bound).max(v.integrality) <= 1e-6, "{v:?}");
                } else {
                    assert!(matches!(
                        sol.status,
                        MipStatus::NodeLimit | MipStatus::TimeLimit
                    ));
                    // A limit may carry a pre-root heuristic incumbent
                    // (verified like any other); otherwise no point at all.
                    if sol.x.is_empty() {
                        assert!(sol.objective.is_nan());
                    } else {
                        let v = p.violation(&sol.x);
                        assert!(v.row.max(v.bound).max(v.integrality) <= 1e-6, "{v:?}");
                    }
                }
            }
        }
        // Root LP is the only node allowed: without heuristics the
        // feasibility search gets none; with them the trivial zero point is
        // already feasible, so the unbounded root settles it (Meyer 1974).
        let (_, sol) = limited(unb, Some(1), None, HeuristicConfig::none());
        assert_eq!(sol.status, MipStatus::NodeLimit);
        assert_eq!(sol.bound, f64::INFINITY);
        let (_, sol) = limited(unb, Some(1), None, HeuristicConfig::default());
        assert_eq!(sol.status, MipStatus::Unbounded);
    }

    #[test]
    fn unbounded_with_free_integers_settled_by_heuristic_point() {
        // rand_294 (mip_diff --only-random --cases 400 --seed 3): unbounded
        // relaxation, free integer columns. The zero-objective feasibility
        // B&B never terminates on unbounded integer boxes; the Feasibility
        // Jump point proves feasibility, hence Unbounded (Meyer 1974).
        let text = "\
NAME rand_294
ROWS
 N  OBJ
 L  R0
 L  R1
 L  R2
COLUMNS
    M0  'MARKER'  'INTORG'
    X0  OBJ  3.0  R1  1.5
    X0  R2  1.0
    X1  OBJ  -1.0  R2  1.0
    X2  OBJ  1.0  R0  1.0
    X2  R2  1.0
    X3  OBJ  -0.5  R2  1.0
    X4  OBJ  1.0  R0  1.0
    X4  R1  1.0  R2  1.0
    M5  'MARKER'  'INTEND'
    X5  OBJ  3.0  R1  4.0
    X5  R2  1.0
RHS
    RHS  R0  3.0  R1  12.5
    RHS  R2  7.5
RANGES
    RNG  R2  21.0
BOUNDS
 FR BND X0
 LO BND X1 -3.0
 FR BND X2
 FR BND X3
 FX BND X4 1.0
 UP BND X5 4.5
ENDATA
";
        let p = MipProblem::from_model(&read_mps(text).unwrap());
        let cfg = MipConfig {
            time_limit: Some(5.0),
            ..MipConfig::default()
        };
        let sol = solve(&p, &cfg).unwrap();
        assert_eq!(sol.status, MipStatus::Unbounded);
        let v = p.violation(&sol.x);
        assert!(v.bound <= 1e-6 && v.row <= 1e-6 && v.integrality <= 1e-6);
    }

    #[test]
    fn heuristics_off_and_on_agree_with_enumeration() {
        let mut rng = XorShift64::new(4242);
        let configs = [
            HeuristicConfig::default(),
            HeuristicConfig::none(),
            HeuristicConfig::all(),
        ]
        .map(|heuristics| MipConfig {
            heuristics,
            ..MipConfig::default()
        });
        for case in 0..150 {
            let n = 2 + (rng.next_u64() % 5) as usize;
            let m = 1 + (rng.next_u64() % 4) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let sols = configs.each_ref().map(|c| solve(&p, c).unwrap());
            let z = enumerate(&p);
            for (k, s) in sols.iter().enumerate() {
                assert_eq!(s.status, sols[0].status, "case {case} config {k}");
                if s.status == MipStatus::Optimal {
                    assert!(
                        (s.objective - sols[0].objective).abs() <= 1e-6,
                        "case {case}"
                    );
                }
                if let Some(z) = z {
                    assert!((p.caller_objective(s.objective) - z).abs() <= 1e-6);
                }
            }
        }
    }

    #[test]
    fn mixed_integer_with_continuous_part() {
        // min −x − y, x int, y cont: x + y ≤ 3.5, x − y ≤ 0.6, y ≤ 2.2 ⇒
        // x = 1 (x ≤ y + 0.6 ≤ 2.8 ⇒ x ≤ 2: x=2 ⇒ y ≥ 1.4, y ≤ 1.5 ⇒ −3.5;
        // x=1 ⇒ y ≤ 2.2 ⇒ −3.2) ⇒ optimum −3.5.
        let sol = solve_mps(
            "\
NAME MIX
ROWS
 N  OBJ
 L  R1
 L  R2
COLUMNS
    M  'MARKER'  'INTORG'
    X  OBJ  -1  R1  1
    X  R2  1
    M  'MARKER'  'INTEND'
    Y  OBJ  -1  R1  1
    Y  R2  -1
RHS
    RHS  R1  3.5  R2  0.6
BOUNDS
 UP BND X 10
 UP BND Y 2.2
ENDATA
",
        );
        assert_eq!(sol.status, MipStatus::Optimal);
        assert!((sol.objective + 3.5).abs() < 1e-9, "{}", sol.objective);
    }

    /// `--localmip` wiring: the improver polishes the pre-root trivial
    /// incumbent. `node_limit: Some(0)` stops the tree before any LP, so the only incumbent path is
    /// the pre-root one — flag off must leave the trivial zero point, flag
    /// on must improve it with a verified point.
    #[test]
    fn localmip_lever_polishes_the_pre_root_incumbent() {
        let inf = f64::INFINITY;
        // max 5a + 4b + 3c s.t. 2a + 3b + c ≤ 5 (the localmip.rs lift-move
        // instance): the trivial incumbent is the all-zero point, optimum −9.
        let p = MipProblem {
            num_col: 3,
            num_row: 1,
            matrix: CscMatrix::from_triplets(1, 3, &[(0, 0, 2.0), (0, 1, 3.0), (0, 2, 1.0)]),
            row_lower: vec![-inf],
            row_upper: vec![5.0],
            cost: vec![-5.0, -4.0, -3.0],
            maximize: false,
            col_lower: vec![0.0; 3],
            col_upper: vec![1.0; 3],
            integer: vec![true; 3],
        };
        let base = MipConfig {
            node_limit: Some(0),
            heuristics: HeuristicConfig {
                // FJ's objective continuation finds this optimum itself;
                // with it off, the only improvement path is the lever.
                feasibility_jump: false,
                fj_restart_effort: 0.0,
                ..HeuristicConfig::default()
            },
            ..MipConfig::default()
        };
        let off = solve(&p, &base).unwrap();
        assert_eq!(off.objective, 0.0, "flag off keeps the trivial point");
        let on = solve(
            &p,
            &MipConfig {
                localmip: true,
                ..base
            },
        )
        .unwrap();
        let v = p.violation(&on.x);
        assert!(
            v.bound <= 1e-6 && v.row <= 1e-6 && v.integrality <= 1e-6,
            "{v:?}"
        );
        assert!(
            on.objective < off.objective,
            "flag on polished {} off the trivial point",
            on.objective
        );
    }

    #[test]
    fn node_limit_reports_valid_bound() {
        let p = random_pure_integer(&mut XorShift64::new(7), 6, 4);
        let cfg = MipConfig {
            node_limit: Some(1),
            ..MipConfig::default()
        };
        let sol = solve(&p, &cfg).unwrap();
        if sol.status == MipStatus::NodeLimit {
            let full = solve(&p, &MipConfig::default()).unwrap();
            assert!(sol.bound <= full.objective + 1e-9);
        }
    }

    /// Random pure-integer problem with small finite boxes, for enumeration.
    fn random_pure_integer(rng: &mut XorShift64, n: usize, m: usize) -> MipProblem {
        let pick =
            |rng: &mut XorShift64, xs: &[f64]| xs[(rng.next_u64() % xs.len() as u64) as usize];
        let mut trip = Vec::new();
        for i in 0..m {
            for j in 0..n {
                let a = pick(rng, &[0.0, 0.0, 1.0, 2.0, -1.0, 3.0, -2.0, 0.5]);
                if a != 0.0 {
                    trip.push((i, j, a));
                }
            }
        }
        let inf = f64::INFINITY;
        let (mut row_lower, mut row_upper) = (Vec::new(), Vec::new());
        for _ in 0..m {
            let b = pick(rng, &[-2.0, 0.0, 1.0, 3.0, 4.5, 6.0]);
            let (lo, hi) = match rng.next_u64() % 4 {
                0 => (-inf, b),
                1 => (b, inf),
                2 => (b, b),
                _ => (b - 2.0, b + 1.5),
            };
            row_lower.push(lo);
            row_upper.push(hi);
        }
        let col_lower: Vec<f64> = (0..n).map(|_| pick(rng, &[-2.0, -1.0, 0.0, 0.0])).collect();
        let col_upper = col_lower
            .iter()
            .map(|&l| l + pick(rng, &[1.0, 2.0, 3.0, 4.0]))
            .collect();
        MipProblem {
            num_col: n,
            num_row: m,
            matrix: CscMatrix::from_triplets(m, n, &trip),
            row_lower,
            row_upper,
            cost: (0..n)
                .map(|_| pick(rng, &[-3.0, -1.0, 0.0, 1.0, 2.0, 1.5]))
                .collect(),
            maximize: rng.next_u64().is_multiple_of(4),
            col_lower,
            col_upper,
            integer: vec![true; n],
        }
    }

    /// Exhaustive enumeration over the integer box: an oracle independent
    /// of the LP and the tree.
    fn enumerate(p: &MipProblem) -> Option<f64> {
        let mut x: Vec<f64> = p.col_lower.clone();
        let mut best: Option<f64> = None;
        loop {
            let v = p.violation(&x);
            if v.row <= 1e-9 {
                let z = p.objective(&x);
                best = Some(best.map_or(z, |b: f64| b.min(z)));
            }
            let mut j = 0;
            loop {
                if j == p.num_col {
                    return best;
                }
                if x[j] < p.col_upper[j] {
                    x[j] += 1.0;
                    break;
                }
                x[j] = p.col_lower[j];
                j += 1;
            }
        }
    }

    #[test]
    fn matches_enumeration_on_random_pure_integer() {
        let mut rng = XorShift64::new(20261022);
        let (mut feasible, mut infeasible) = (0, 0);
        for case in 0..300 {
            let n = 2 + (rng.next_u64() % 5) as usize;
            let m = 1 + (rng.next_u64() % 4) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let sol = solve(&p, &MipConfig::default()).unwrap();
            match enumerate(&p) {
                Some(z) => {
                    feasible += 1;
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    let ours = p.caller_objective(sol.objective);
                    assert!((ours - z).abs() <= 1e-6, "case {case}: {ours} vs {z}");
                    let v = p.violation(&sol.x);
                    assert!(v.row <= 1e-6 && v.bound <= 1e-6 && v.integrality <= 1e-6);
                }
                None => {
                    infeasible += 1;
                    assert_eq!(sol.status, MipStatus::Infeasible, "case {case}");
                }
            }
        }
        assert!(feasible > 50 && infeasible > 20, "{feasible}/{infeasible}");
    }

    /// Conflicts learnt from Farkas proofs only ever remove infeasible
    /// points: statuses and optima match enumeration with them on,
    /// presolve off so the tree is deep.
    #[test]
    fn conflict_matches_enumeration_on_random_pure_integer() {
        let mut rng = XorShift64::new(20261024);
        let mut cfg = MipConfig {
            conflict: true,
            presolve: false,
            conflict_cutoff: true,
            ..MipConfig::default()
        };
        // Gomory (on by default) prunes these tiny trees before any
        // conflict arises; this test wants a deep tree to exercise conflicts.
        cfg.cuts.gomory = false;
        let mut stats = ConflictStats::default();
        for case in 0..400 {
            let n = 3 + (rng.next_u64() % 6) as usize;
            let m = 1 + (rng.next_u64() % 5) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let sol = solve(&p, &cfg).unwrap();
            stats.proofs += sol.conflict.proofs;
            stats.cutoff_proofs += sol.conflict.cutoff_proofs;
            stats.added += sol.conflict.added;
            stats.tightenings += sol.conflict.tightenings;
            stats.pruned += sol.conflict.pruned;
            match enumerate(&p) {
                Some(z) => {
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    let ours = p.caller_objective(sol.objective);
                    assert!((ours - z).abs() <= 1e-6, "case {case}: {ours} vs {z}");
                }
                None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
            }
        }
        assert!(stats.added > 0, "conflicts never exercised: {stats:?}");
    }

    /// Sub-MIP heuristics only ever add verified incumbents: statuses and
    /// optima match enumeration with them on.
    #[test]
    fn submip_matches_enumeration_on_random_pure_integer() {
        let mut rng = XorShift64::new(20261023);
        let cfg = MipConfig {
            submip: true,
            ..MipConfig::default()
        };
        for case in 0..150 {
            let n = 2 + (rng.next_u64() % 6) as usize;
            let m = 1 + (rng.next_u64() % 4) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let sol = solve(&p, &cfg).unwrap();
            match enumerate(&p) {
                Some(z) => {
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    let ours = p.caller_objective(sol.objective);
                    assert!((ours - z).abs() <= 1e-6, "case {case}: {ours} vs {z}");
                }
                None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
            }
        }
    }

    /// Warm steepest edge with carried weights (cuts on, presolve off so the
    /// tree is deep) keeps the enumeration optima.
    #[test]
    fn warm_dse_matches_enumeration_on_random_pure_integer() {
        let mut rng = XorShift64::new(20261025);
        let cfg = MipConfig {
            warm_dse: true,
            presolve: false,
            ..MipConfig::default()
        };
        for case in 0..200 {
            let n = 3 + (rng.next_u64() % 6) as usize;
            let m = 1 + (rng.next_u64() % 5) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let sol = solve(&p, &cfg).unwrap();
            match enumerate(&p) {
                Some(z) => {
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    let ours = p.caller_objective(sol.objective);
                    assert!((ours - z).abs() <= 1e-6, "case {case}: {ours} vs {z}");
                }
                None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
            }
        }
    }

    /// The root node's cold LPs on the dual simplex keep the
    /// enumeration optima.
    #[test]
    fn root_lp_dual_matches_enumeration_on_random_pure_integer() {
        let mut rng = XorShift64::new(20261029);
        let cfg = MipConfig {
            root_lp_dual: true,
            presolve: false,
            ..MipConfig::default()
        };
        for case in 0..200 {
            let n = 3 + (rng.next_u64() % 6) as usize;
            let m = 1 + (rng.next_u64() % 5) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let sol = solve(&p, &cfg).unwrap();
            match enumerate(&p) {
                Some(z) => {
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    let ours = p.caller_objective(sol.objective);
                    assert!((ours - z).abs() <= 1e-6, "case {case}: {ours} vs {z}");
                }
                None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
            }
        }
    }

    /// With a time limit a root LP cut off by the configured iteration
    /// cap re-runs uncapped on the dual (HiGHS gives the root LP no cap and
    /// toggles the simplex on failure); without one the cap still ends it
    /// with the iteration-limit error.
    #[test]
    fn root_lp_failure_retries_on_the_other_simplex_under_a_time_limit() {
        let mut rng = XorShift64::new(20261030);
        let p = (0..50)
            .map(|_| random_pure_integer(&mut rng, 8, 5))
            .find(|p| {
                let mut e = RusolveEngine::new(p, MipConfig::default().lp);
                e.solve(&p.col_lower, &p.col_upper)
                    .is_ok_and(|o| o.status == LpStatus::Optimal && o.iterations >= 2)
            })
            .expect("a root LP that needs pivots");
        let base = MipConfig {
            presolve: false,
            lp: SimplexConfig {
                max_iterations: 1,
                ..MipConfig::default().lp
            },
            ..MipConfig::default()
        };
        assert!(matches!(
            solve(&p, &base),
            Err(MipError::Lp(LpError::IterationLimit))
        ));
        let timed = MipConfig {
            time_limit: Some(60.0),
            ..base
        };
        let sol = solve(&p, &timed).unwrap();
        assert!(
            matches!(sol.status, MipStatus::Optimal | MipStatus::Infeasible),
            "{:?}",
            sol.status
        );
        assert!(sol.lp_stats.root_retries >= 1);
    }

    /// Tree node LPs that cannot finish (every solve capped at one
    /// pivot) never abort the solve: the search branches without them and
    /// re-solves all-integers-fixed nodes on the other simplex, uncapped
    /// (HighsSearch.cpp:1347–1476), and still matches enumeration.
    #[test]
    fn unsolved_node_lps_never_abort_the_search() {
        let mut rng = XorShift64::new(20261029);
        let cfg = MipConfig {
            presolve: false,
            time_limit: Some(60.0),
            // Cold node LPs, so the cap reaches them (warm reopts finish);
            // no root heuristics, whose basis-restoring root re-solve is not
            // a tree node LP.
            warm_start: false,
            heuristics: HeuristicConfig::none(),
            lp: SimplexConfig {
                max_iterations: 1,
                ..MipConfig::default().lp
            },
            ..MipConfig::default()
        };
        let mut no_lp = 0;
        for case in 0..60 {
            let n = 2 + (rng.next_u64() % 6) as usize;
            let m = 1 + (rng.next_u64() % 4) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let sol = solve(&p, &cfg).unwrap_or_else(|e| panic!("case {case}: {e}"));
            assert_eq!(sol.lp_stats.unsolved_pruned, 0, "case {case}");
            no_lp += sol.lp_stats.no_lp_branches;
            match enumerate(&p) {
                Some(z) => {
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    let ours = p.caller_objective(sol.objective);
                    assert!((ours - z).abs() <= 1e-6, "case {case}: {ours} vs {z}");
                }
                None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
            }
        }
        assert!(no_lp > 0, "no node LP failed: the test exercises nothing");
    }

    /// ALNS on (crate::alns) keeps the enumeration optima.
    #[test]
    fn alns_matches_enumeration_on_random_pure_integer() {
        let mut rng = XorShift64::new(20261024);
        let cfg = MipConfig {
            alns: true,
            ..MipConfig::default()
        };
        for case in 0..150 {
            let n = 2 + (rng.next_u64() % 6) as usize;
            let m = 1 + (rng.next_u64() % 4) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let sol = solve(&p, &cfg).unwrap();
            match enumerate(&p) {
                Some(z) => {
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    let ours = p.caller_objective(sol.objective);
                    assert!((ours - z).abs() <= 1e-6, "case {case}: {ours} vs {z}");
                }
                None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
            }
        }
    }

    /// Racing (crate::race) gives the single search's answers: statuses and
    /// objectives match enumeration, with and without presolve, and a node
    /// limit reports a bound no better than the optimum.
    #[test]
    fn racing_matches_enumeration_on_random_pure_integer() {
        let mut rng = XorShift64::new(20261023);
        for case in 0..150 {
            let n = 2 + (rng.next_u64() % 5) as usize;
            let m = 1 + (rng.next_u64() % 4) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let config = MipConfig {
                threads: 2 + case % 3,
                presolve: case % 2 == 0,
                ..MipConfig::default()
            };
            let sol = solve(&p, &config).unwrap();
            match enumerate(&p) {
                Some(z) => {
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    let ours = p.caller_objective(sol.objective);
                    assert!((ours - z).abs() <= 1e-6, "case {case}: {ours} vs {z}");
                    assert_eq!(sol.bound, sol.objective, "case {case}");
                    let v = p.violation(&sol.x);
                    assert!(v.row <= 1e-6 && v.bound <= 1e-6 && v.integrality <= 1e-6);
                    let limited = solve(
                        &p,
                        &MipConfig {
                            node_limit: Some(1),
                            ..config
                        },
                    )
                    .unwrap();
                    assert!(
                        p.caller_objective(limited.bound) <= z + 1e-6,
                        "case {case}: bound {} above optimum {z}",
                        limited.bound
                    );
                }
                None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
            }
        }
    }

    /// The unbounded-relaxation paths under racing: a proven Unbounded MIP
    /// and an Infeasible one behind an unbounded relaxation.
    #[test]
    fn racing_settles_unbounded_relaxations() {
        let unbounded = "\
NAME UNB
OBJSENSE
 MAX
ROWS
 N  OBJ
 L  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  OBJ  1  R  1
    Y  R  -2
    M  'MARKER'  'INTEND'
RHS
    RHS  R  0
BOUNDS
 LI BND X 0
 LI BND Y 0
ENDATA
";
        let infeasible = "\
NAME UNB_INF
ROWS
 N  OBJ
 E  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  R  2
    Y  R  -2
    M  'MARKER'  'INTEND'
    Z  OBJ  -1
RHS
    RHS  R  1
BOUNDS
 UP BND X 6
 UP BND Y 6
ENDATA
";
        let config = MipConfig {
            threads: 4,
            ..MipConfig::default()
        };
        let p = MipProblem::from_model(&read_mps(unbounded).unwrap());
        let sol = solve(&p, &config).unwrap();
        assert_eq!(sol.status, MipStatus::Unbounded);
        assert_eq!(sol.objective, f64::INFINITY);
        let p = MipProblem::from_model(&read_mps(infeasible).unwrap());
        assert_eq!(solve(&p, &config).unwrap().status, MipStatus::Infeasible);
    }

    #[test]
    fn root_incumbent_bounds_start_times_before_cuts() {
        // One machine, three jobs (p = 2, 3, 4), Manne's big-M disjunctions
        // with M = 9: z_ij = 1 puts j after i. Start times and the makespan
        // C have no upper bound until a root incumbent's cutoff propagates.
        let p = MipProblem::from_model(
            &read_mps(
                "\
NAME JS
ROWS
 N  OBJ
 G  MK0
 G  MK1
 G  MK2
 G  A01
 G  B01
 G  A02
 G  B02
 G  A12
 G  B12
COLUMNS
    C  OBJ  1  MK0  1
    C  MK1  1  MK2  1
    S0  MK0  -1  A01  1
    S0  B01  -1  A02  1
    S0  B02  -1
    S1  MK1  -1  A01  -1
    S1  B01  1  A12  1
    S1  B12  -1
    S2  MK2  -1  A02  -1
    S2  B02  1  A12  -1
    S2  B12  1
    M  'MARKER'  'INTORG'
    Z01  A01  9  B01  -9
    Z02  A02  9  B02  -9
    Z12  A12  9  B12  -9
    M  'MARKER'  'INTEND'
RHS
    RHS  MK0  2  MK1  3
    RHS  MK2  4  A01  3
    RHS  B01  -7  A02  4
    RHS  B02  -7  A12  4
    RHS  B12  -6
BOUNDS
 UP BND Z01 1
 UP BND Z02 1
 UP BND Z12 1
ENDATA
",
            )
            .unwrap(),
        );
        let config = MipConfig::default();
        let mut root = root_node(&p, config.integrality_tol).unwrap();
        let mut search = Search::new(
            &p,
            RusolveEngine::new(&p, config.lp),
            &config,
            Instant::now(),
        );
        assert!(root.upper.iter().any(|u| u.is_infinite()));
        let lp = search.engine.solve(&root.lower, &root.upper).unwrap();
        assert!(
            search
                .most_fractional(&lp.x, &root.lower, &root.upper)
                .is_some(),
            "{:?}",
            lp.x
        );
        let lp = search.root_incumbent(&mut root, lp).unwrap();
        let z = search
            .incumbent
            .as_ref()
            .expect("root rounding finds a schedule")
            .0;
        assert!(z >= 9.0 - 1e-9, "{z}");
        assert!(root.upper.iter().all(|u| u.is_finite()), "{:?}", root.upper);
        assert_eq!(search.global_upper, root.upper);
        assert!(lp.status == LpStatus::Infeasible || lp.objective <= search.upper_limit);
    }

    /// Wraps the real engine and pushes the last column 2e-6 below its
    /// value on every warm solve: integral points that fail the absolute
    /// incumbent check while passing the engine's relative one.
    struct Drift {
        inner: RusolveEngine,
        drifted: std::rc::Rc<std::cell::Cell<u32>>,
    }

    impl LpEngine for Drift {
        fn solve(&mut self, lower: &[f64], upper: &[f64]) -> Result<LpOutcome, LpError> {
            self.inner.solve(lower, upper)
        }

        fn solve_from(
            &mut self,
            lower: &[f64],
            upper: &[f64],
            start: Option<&WarmStart>,
            cutoff: f64,
        ) -> Result<LpOutcome, LpError> {
            let mut out = self.inner.solve_from(lower, upper, start, cutoff)?;
            if start.is_some() && out.status == LpStatus::Optimal {
                *out.x.last_mut().unwrap() -= 2e-6;
                self.drifted.set(self.drifted.get() + 1);
            }
            Ok(out)
        }
    }

    #[test]
    fn violated_integral_point_is_resolved_cold_not_an_error() {
        // min x + y + z/1000, 1000x + 1000y + z >= 1500: optimum 1.5 at
        // x + y = 1, z = 500. Every warm integral point violates the row
        // (or z >= 0) by 2e-6; snapping cannot fix a continuous column.
        let p = MipProblem::from_model(
            &read_mps(
                "\
NAME DRIFT
ROWS
 N  OBJ
 G  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  OBJ  1  R  1000
    Y  OBJ  1  R  1000
    M  'MARKER'  'INTEND'
    Z  OBJ  0.001  R  1
RHS
    RHS  R  1500
BOUNDS
 UP BND X 3
 UP BND Y 3
ENDATA
",
            )
            .unwrap(),
        );
        let config = MipConfig {
            heuristics: HeuristicConfig::none(),
            ..MipConfig::default()
        };
        let drifted = std::rc::Rc::new(std::cell::Cell::new(0));
        let engine = Drift {
            inner: RusolveEngine::new(&p, config.lp),
            drifted: std::rc::Rc::clone(&drifted),
        };
        let root = root_node(&p, config.integrality_tol).unwrap();
        let mut search = Search::new(&p, engine, &config, Instant::now());
        search.run(root).expect("no InvalidLpSolution");
        assert!(drifted.get() > 0, "the test must exercise warm solves");
        let (z, x) = search.incumbent.as_ref().expect("an incumbent");
        assert!((z - 1.5).abs() < 1e-6, "{z}");
        let v = p.violation(x);
        assert!(v.row.max(v.bound) <= CANDIDATE_TOL, "{v:?}");
    }

    #[test]
    fn most_fractional_ignores_bound_overshoot() {
        // x = 1.00009 on a binary is integral once clamped into [0, 1];
        // branching on it would make the down child equal its parent.
        let p = MipProblem::from_model(
            &read_mps(
                "\
NAME ONE
ROWS
 N  OBJ
 L  R
COLUMNS
    M  'MARKER'  'INTORG'
    X  OBJ  1  R  1
    Y  OBJ  1  R  1
    M  'MARKER'  'INTEND'
RHS
    RHS  R  1
BOUNDS
 UP BND X 1
 UP BND Y 1
ENDATA
",
            )
            .unwrap(),
        );
        let config = MipConfig::default();
        let search = Search::new(
            &p,
            RusolveEngine::new(&p, config.lp),
            &config,
            Instant::now(),
        );
        let (lower, upper) = ([0.0, 0.0], [1.0, 1.0]);
        assert_eq!(
            search.most_fractional(&[1.00009, 0.0], &lower, &upper),
            None
        );
        assert_eq!(
            search.most_fractional(&[1.00009, 0.4], &lower, &upper),
            Some((1, 0.4))
        );
        assert_eq!(
            fractional_integers(&p.integer, &[1.5, 0.4], &lower, &upper, 1e-6).collect::<Vec<_>>(),
            [(1, 0.4)]
        );
    }

    #[test]
    fn big_m_does_not_turn_integrality_tolerance_into_objective() {
        // min t + z, t + 1e6 z >= 1: the LP optimum z = 1e-6, t = 0 is
        // within integrality tolerance, but z = 0 or 1 both cost 1.
        let p = MipProblem::from_model(
            &read_mps(
                "\
NAME BIGM
ROWS
 N  OBJ
 G  R
COLUMNS
    T  OBJ  1  R  1
    M  'MARKER'  'INTORG'
    Z  OBJ  1  R  1000000
    M  'MARKER'  'INTEND'
RHS
    RHS  R  1
BOUNDS
 UP BND Z 1
ENDATA
",
            )
            .unwrap(),
        );
        let config = MipConfig::default();
        let search = Search::new(
            &p,
            RusolveEngine::new(&p, config.lp),
            &config,
            Instant::now(),
        );
        let lp_point = [0.0, 1e-6];
        assert!(search.brancher.fractional(&lp_point).is_empty());
        assert_eq!(search.unsnappable(&lp_point), Some(1));
        assert_eq!(search.unsnappable(&[0.0, 1.0 - 1e-7]), None);
        assert_eq!(snap_integers(vec![0.5, 1.0 - 1e-7], &p.integer), [0.5, 1.0]);
        let sol = solve(&p, &config).unwrap();
        assert!((sol.objective - 1.0).abs() < 1e-9, "{}", sol.objective);
    }

    /// The combinatorial conflict `x0 = 1 ⇒ x1 = x2 = 0 ⇒ COVER
    /// broken` is invisible to per-row reasoning (presolve and the root LP
    /// leave x0 free), but one probe fixes x0 = 0. The flag must deliver
    /// that fixing, keep the optimum, and stay out of the default path.
    #[test]
    fn root_probing_fixes_conflicted_binary_and_keeps_optimum() {
        let mps = "\
NAME PROBE\n\
ROWS\n\
 N  OBJ\n\
 L  PICK\n\
 G  COVER\n\
 L  SLACKR\n\
COLUMNS\n\
    X0  OBJ 1  PICK 1\n\
    X0  SLACKR 1\n\
    X1  OBJ 1  PICK 1  COVER 1\n\
    X2  OBJ 1  PICK 1  COVER 1\n\
    T   SLACKR 1\n\
RHS\n\
    RHS PICK 1  COVER 1  SLACKR 5\n\
BOUNDS\n\
    BV BND X0\n\
    BV BND X1\n\
    BV BND X2\n\
    UP BND T 10\n\
ENDATA\n";
        let p = MipProblem::from_model(&read_mps(mps).unwrap());
        let base = MipConfig {
            presolve: false,
            ..MipConfig::default()
        };
        let off = solve(&p, &base).unwrap();
        assert_eq!(off.status, MipStatus::Optimal);
        assert_eq!(off.probing.probed, 0, "default path never probes");
        let on = solve(
            &p,
            &MipConfig {
                root_probing: true,
                ..base
            },
        )
        .unwrap();
        assert_eq!(on.status, MipStatus::Optimal);
        assert!((on.objective - off.objective).abs() < 1e-9);
        assert!(on.probing.probed > 0);
        assert!(on.probing.fixings >= 1, "{:?}", on.probing);
        assert_eq!(on.probing.infeasible_probes, 1, "{:?}", on.probing);
        // x0 = 0 in the returned point (0; exactly one of x1/x2 = 1).
        assert!(on.x[0].abs() < 1e-9);
        // And with the full default machinery (presolve on) the flag stays
        // compatible: same answer, probing still counted.
        let full = solve(
            &p,
            &MipConfig {
                root_probing: true,
                ..MipConfig::default()
            },
        )
        .unwrap();
        assert_eq!(full.status, MipStatus::Optimal);
        assert!((full.objective - off.objective).abs() < 1e-9);
    }
}
