//! Reliability branching with strong branching on unreliable candidates, and
//! the child-direction rule.
//!
//! Algorithm (Achterberg, Koch & Martin, "Branching rules revisited", Oper.
//! Res. Lett. 33:42–54, 2005) as realised in HiGHS
//! `HighsSearch::selectBranchingCandidate` (reference/HiGHS/highs/mip/
//! HighsSearch.cpp:247–707, MIT licence, Copyright (c) HiGHS developers),
//! read as a reference text and reimplemented:
//!
//! * candidates with reliable pseudocosts get their estimated gains; the
//!   others start at +∞ and are strong-branched one side at a time, best
//!   current score first, until the best candidate is reliable on both sides
//!   or the strong-branching iteration budget runs out (:388–404);
//! * each strong-branch LP also scores every other fractional candidate whose
//!   value in that LP solution already sits on one of its own branch sides
//!   (`analyzeSolution` :414–533);
//! * a side that propagation or the LP proves infeasible / cut off reduces
//!   the node to the other side instead of branching (:556–675);
//! * strong-branch LPs are warm dual probes from the node basis under the
//!   incumbent cutoff (HiGHS `Playground::solveLp`, HighsLpRelaxation.h:129,
//!   objective limit :405); a probe stopped by the cutoff counts as an
//!   infeasible side, as HiGHS maps `kObjectiveBound` to infeasible
//!   (HighsLpRelaxation.cpp:1263);
//! * budget: `base + ½·(non-strong-branching LP iterations)`, with
//!   `min_reliable` shrunk linearly once half of it is used and set to 0
//!   beyond it (`HighsSearch::branch` :1112–1127) — but payoff-aware:
//!   the ceiling only tightens toward SCIP's 50% quota
//!   (`sbiterquot`/`sbiterofs`, branch_relpscost.c) after a window of dry
//!   probes, and never while probes keep cutting sides off or moving the
//!   bound, because such strong branching is load-bearing (neos-860300:
//!   trimming it multiplies the tree 4–8×);
//! * investment depth ([`BranchConfig::sb_invest`]): the
//!   reliability threshold is not a constant while strong branching fits
//!   SCIP's sbiterquot budget `sbiterofs + sbiterquot·(LP iterations)` —
//!   it grows to [`SB_INVEST_MAX_MINREL`] so pseudo-costs are calibrated
//!   on strong-branch probes deeper into the tree, falling back to the
//!   configured base past the quota ([`invest_depth`]). The tree-loss
//!   families (lights_out at 10–31% share) measured 0.2–0.6× nodes; the
//!   high-share ones (fcn[fixed] 73%) stay at the base;
//! * child direction: HiGHS's default `kRootSol` rule — average inferences,
//!   weighted by how far the LP value moved from the root (or sub-root)
//!   solution (:1166–1206) — or nearest rounding, our measured default;
//! * a node whose LP fails is branched on pseudocost scores alone
//!   (:1347–1418, [`Brancher::select_without_lp`]).
//!
//! Beyond HiGHS, from Achterberg, Koch & Martin 2005 §3 and Achterberg's
//! thesis (*Constraint Integer Programming*, TU Berlin 2007, §5.4): an
//! optional lookahead stop (λ), a per-node cap on strong-branched candidates,
//! unscored candidates taken in pseudocost-score order (HiGHS takes them in
//! column order), and a per-probe dual-simplex iteration cap whose dual
//! bound stands in for the exact gain when the cap is hit. HiGHS caps a
//! probe only at its tree-LP limit and scores a capped probe 0. The cap can
//! scale with the node LPs a probe refines — clamp(2 × avg node-LP
//! iterations, 10, 500), [`SbProbeLimit::Adaptive`] — but
//! measured strictly negative in both directions (raising it to ≈150–470
//! cost +3–28% time on the dense MIPLIB slice neos-860300/air05, whose
//! capped probes' dual-bound scores cover more candidates per budget;
//! flooring it at 10 truncated lights_out's tiny warm probes, 2.2× nodes),
//! so the default stays the fixed 100 and `Adaptive` is the
//! `sbiters=adaptive` opt-in.
//!
//! Not ported: the per-column open-node counts tie-break in candidate
//! selection (needs a column index of the queue), the LP degeneracy
//! adjustment of `min_reliable` (needs node-LP duals), and
//! conflict scores.

use crate::engine::{LpError, LpOutcome, LpStatus};
pub use crate::nodequeue::NodeSelection;
use crate::pseudocost::{DEFAULT_MIN_RELIABLE, Pseudocost};

/// Strong-branching budget base of both our default and the
/// [`BranchConfig::highs`] preset (HighsSearch.cpp:1116; the cold variant used 10k).
/// Cutting it to SCIP's `sbiterofs` 1 000 unconditionally measured strictly
/// negative where it binds — neos-860300 45 s Optimal → time limit (61% of
/// its LP iterations are strong branching, and they are load-bearing),
/// air05 +38% time, lights_out 7× nodes — so the 50% quota now applies only
/// through the dry-probe payoff guard ([`SB_DRY_PROBE_LIMIT`]).
pub const HIGHS_SB_BASE_ITERATIONS: u64 = 100_000;
/// Reliability threshold, HiGHS's `mip_pscost_minreliable` 8. Warm sweeps
/// found no better value: 4 and 16
/// move solved counts by ±1 of 54 on the mixed tier, within run-to-run noise.
/// Count-stability halving (once pseudocosts hold `min_reliable`
/// observations per integer column) measured 5.6× nodes on lights_out — an
/// average count says nothing about per-column stability — and was dropped.
pub const SB_MIN_RELIABLE: u32 = DEFAULT_MIN_RELIABLE;
/// Dual-simplex iterations per warm probe; a capped probe scores with its
/// dual bound (Achterberg 2007, §5.4). HiGHS caps probes only at the tree LP
/// cap. Measured optimum of {25, 50, 100, 200, 400} on the mixed tier (43/54
/// solved in 10 s; 1.7% of probes capped) and kept as the default despite
/// the adaptive cap: clamp(2 × avg node-LP iterations, 10, 500)
/// (`sbiters=adaptive`) measured +3.3% LP iterations on neos-860300 (avg 74
/// per node LP → cap 148; its 530 capped probes' dual bounds are as good a
/// score at ⅔ the cost), +26% time on air05 (69% of probes capped at 100,
/// avg 234 per node LP), and 2.2× nodes on lights_out (avg 2.8 per node LP →
/// cap 10 truncates its tiny warm probes mid-dual).
pub const SB_PROBE_FIXED_ITERATIONS: u64 = 100;
/// Clamp bounds of the [`SbProbeLimit::Adaptive`] cap: a
/// probe refines one node LP, so it may cost a couple of them, but never
/// less than a degenerate dual step's worth of escape (10) and never half a
/// tree LP (500).
pub const SB_PROBE_MIN_ITERATIONS: u64 = 10;
pub const SB_PROBE_MAX_ITERATIONS: u64 = 500;
/// Consecutive probes with neither a cut-off side nor an appreciable bound
/// move before the strong-branching budget tightens toward the 50% quota
/// (dry probes are the only ones the quota may trim).
pub const SB_DRY_PROBE_LIMIT: u32 = 16;
/// Always-allowed strong-branching iterations under the tightened quota:
/// SCIP `branch_relpscost`'s `sbiterquot` quota is 50% of all LP iterations
/// (`sbiterofs` 1 000 is its offset).
pub const SB_QUOTA_OFFSET: u64 = 1_000;
/// Reliability depth the [`BranchConfig::sb_invest`] policy grows to while
/// strong branching is cheap relative to the quota: 4× the
/// base 8. Measured on the tree-loss families — lights_out s3
/// 23 411 → 4 560 nodes (0.195×) at 0.43× time, gap s1 0.40× nodes,
/// general_ip 0.58× — because pseudo-costs then keep being calibrated on
/// strong-branch probes through the whole tree instead of going stale after
/// the warmup nodes (HiGHS's 218-node lights_out solve never leaves its
/// warmup: 39 SB-iterations/node; ours left it after ~1 000 probes and ran
/// the remaining 22 000 nodes on those early estimates).
pub const SB_INVEST_MAX_MINREL: u32 = 32;
/// SB-iteration share below which the investment policy holds
/// [`SB_INVEST_MAX_MINREL`] — unused by the binary quota arm itself, kept
/// as documentation of the measured headroom of the winning families
/// (lights_out 10–31%, gap 19% final share at default).
pub const SB_INVEST_SHARE_FULL: f64 = 0.25;
/// Always-allowed strong-branching iterations of the investment quota:
/// SCIP `branch_relpscost`'s `sbiterofs` default, shared with the
/// dry-probe guard.
pub const SB_INVEST_QUOTA_OFFSET: u64 = SB_QUOTA_OFFSET;
/// SB-iteration share of all LP iterations the investment policy may
/// consume, on top of SCIP's offset — the `sbiterquot` 0.5 quota itself
/// (`branch_relpscost.c`; the dry-probe guard already codified its constants). Families
/// already spending more than half their iterations probing
/// (fixed_charge_network[fixed] 73%, block_angular 70% at default)
/// measured node-neutral and time-negative on more depth; the offset
/// keeps the root warmup's structural share spike from disarming the
/// policy on probe-cheap families (lights_out's root reads 72% and still
/// solves 4× smaller with depth).
pub const SB_INVEST_SHARE_BASE: f64 = 0.50;

/// Per-probe strong-branching iteration cap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SbProbeLimit {
    /// The tree LP iteration cap, as in HiGHS: it caps a probe only there
    /// (`Playground::solveLp`, HighsLpRelaxation.h:129).
    Tree,
    /// A fixed iteration count (the A/B lever `sbiters=n`).
    Fixed(u64),
    /// clamp(2 × average node-LP iterations, 100, 500): a probe
    /// refines one node LP, so on instances with expensive node LPs it may
    /// cost a couple of them; on cheap warm nodes it stays at the measured
    /// fixed cap (100). Before the first node LP, the root LP's iterations
    /// stand in.
    Adaptive,
}

impl SbProbeLimit {
    /// Resolves the cap: `avg_node_iters` is the measured average node-LP
    /// iteration count (0 before the first node LP), `root_iters` the root
    /// LP's, `tree_cap` the tree LP iteration cap.
    pub fn resolve(&self, avg_node_iters: u64, root_iters: u64, tree_cap: u64) -> u64 {
        match *self {
            Self::Tree => tree_cap,
            Self::Fixed(n) => n,
            Self::Adaptive => {
                let avg = if avg_node_iters > 0 {
                    avg_node_iters
                } else {
                    root_iters
                };
                avg.saturating_mul(2)
                    .clamp(SB_PROBE_MIN_ITERATIONS, SB_PROBE_MAX_ITERATIONS)
            }
        }
    }
}
/// Candidates strong-branched per node, 0 = no cap. Caps of 8 and 20 gave no
/// consistent gain (mixed: 42 vs 43 solved with `sbiters=100`), so we keep
/// HiGHS's uncapped loop and let the lookahead stop it.
pub const SB_MAX_CANDIDATES: u32 = 0;
/// Lookahead λ = 8 (Achterberg, Koch & Martin 2005, §3), the value
/// that measured best on strong-branching-dominated families (set_cover). λ = 4 lost a solve on mixed and λ = 0 or 16 were
/// neutral.
pub const SB_LOOKAHEAD: u32 = 8;
/// HiGHS `small_matrix_value`, its `epsilon` in the search.
const EPSILON: f64 = 1e-9;

const UP: u8 = 1;
const DOWN: u8 = 2;

/// What strong branching needs from the tree search.
pub trait SbHost {
    /// Solves the relaxation under the given column bounds.
    fn sb_solve(&mut self, lower: &[f64], upper: &[f64]) -> Result<LpOutcome, LpError>;
    /// Iteration cap of one [`SbHost::sb_solve`]; a capped solve is charged
    /// this much against the strong-branching budget.
    fn sb_iteration_limit(&self) -> u64;
    /// Propagates after a bound change on `col`; `false` iff infeasible.
    fn sb_propagate(&mut self, lower: &mut [f64], upper: &mut [f64], col: usize) -> bool;
    /// Current pruning limit: nodes with a larger bound are discarded.
    fn upper_limit(&self) -> f64;
    /// An integer-feasible LP point found while strong branching.
    fn sb_integral(&mut self, x: Vec<f64>);
    /// A global limit (time) was hit: finish with the scores at hand.
    fn sb_stop(&self) -> bool;
    /// Nodes processed by the search so far: the investment policy arms
    /// only on trees large enough to amortize its probes
    /// ([`SB_INVEST_ARM_NODES`]).
    fn sb_nodes(&self) -> u64;
}

/// Which child a dive enters first after branching.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChildRule {
    /// HiGHS's default `kRootSol` (HighsSearch.cpp:1166–1206): average
    /// inferences per side, weighted by how far the LP value moved away from
    /// the root / sub-root LP solution.
    RootSol,
    /// The side the LP value rounds to (ties up).
    Nearest,
}

/// Knobs of branching and node selection. The HiGHS rules are all here;
/// the defaults are what measured best on our generator families (measured on our benchmark set).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BranchConfig {
    /// Pseudocost observations per side before strong branching is skipped
    /// for a candidate; 0 disables strong branching (HiGHS
    /// `mip_pscost_minreliable`, default 8).
    pub min_reliable: u32,
    /// Base of the strong-branching LP-iteration budget, on top of half the
    /// non-strong-branching iterations. SCIP's `sbiterofs` (1 000); HiGHS
    /// uses 100 000 (`HIGHS_SB_BASE_ITERATIONS`).
    pub sb_base_iterations: u64,
    /// Dual-simplex iterations per strong-branch probe; a capped probe
    /// scores with its dual bound (Achterberg 2007 §5.4). `Adaptive` =
    /// clamp(2 × avg node-LP iterations, 100, 500); `Tree` =
    /// the tree LP cap, as in HiGHS.
    pub sb_iteration_limit: SbProbeLimit,
    /// Candidates strong-branched per node at most; 0 = no cap, as in HiGHS.
    pub sb_max_candidates: u32,
    /// Stop strong branching once this many strong-branch LPs in a row left
    /// the best fully scored candidate unchanged (the lookahead λ of
    /// Achterberg, Koch & Martin 2005, §3); 0 = off, as in HiGHS, whose loop
    /// stops only on the budget (HighsSearch.cpp:390–399).
    pub sb_lookahead: u32,
    /// Reliability-depth investment: grow the effective
    /// `min_reliable` to [`SB_INVEST_MAX_MINREL`] while strong branching
    /// fits SCIP's sbiterquot budget ([`invest_depth`]), back to the
    /// configured base past it. The tree-loss families carry low SB
    /// shares (lights_out 10–31%, gap 19%) and won big from depth; the
    /// high-share ones (fcn[fixed] 73%) only paid time.
    pub sb_invest: bool,
    pub child: ChildRule,
    pub node_selection: NodeSelection,
}

impl BranchConfig {
    /// HiGHS 1.15's defaults, as read from its source.
    pub fn highs() -> Self {
        Self {
            min_reliable: DEFAULT_MIN_RELIABLE,
            sb_base_iterations: HIGHS_SB_BASE_ITERATIONS,
            sb_iteration_limit: SbProbeLimit::Tree,
            sb_max_candidates: 0,
            sb_lookahead: 0,
            sb_invest: false,
            child: ChildRule::RootSol,
            node_selection: NodeSelection::Hybrid,
        }
    }

    /// Applies `key=value` overrides separated by commas (the CLI's
    /// `--branch` / `RUSOLVE_BRANCH` lever for A/B runs): `minrel`,
    /// `sbbase`, `sbiters` (`adaptive`/`tree`/a count), `sbcands`,
    /// `lookahead`, `invest` (`on`/`off`), `child` (`rootsol`/`nearest`),
    /// `nodesel` (`hybrid`/`bestbound`), or the preset `highs`.
    pub fn apply(&mut self, spec: &str) -> Result<(), String> {
        for item in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            if item == "highs" {
                *self = Self::highs();
                continue;
            }
            let (k, v) = item
                .split_once('=')
                .ok_or_else(|| format!("branch option '{item}' is not key=value"))?;
            let bad = || format!("bad value for branch option {k}: {v}");
            match k {
                "minrel" => self.min_reliable = v.parse().map_err(|_| bad())?,
                "sbbase" => self.sb_base_iterations = v.parse().map_err(|_| bad())?,
                "sbiters" if v == "tree" => self.sb_iteration_limit = SbProbeLimit::Tree,
                "sbiters" if v == "adaptive" => self.sb_iteration_limit = SbProbeLimit::Adaptive,
                "sbiters" => {
                    self.sb_iteration_limit = SbProbeLimit::Fixed(v.parse().map_err(|_| bad())?)
                }
                "sbcands" => self.sb_max_candidates = v.parse().map_err(|_| bad())?,
                "lookahead" => self.sb_lookahead = v.parse().map_err(|_| bad())?,
                "invest" => {
                    self.sb_invest = match v {
                        "on" | "true" => true,
                        "off" | "false" => false,
                        _ => return Err(bad()),
                    }
                }
                "child" => {
                    self.child = match v {
                        "rootsol" => ChildRule::RootSol,
                        "nearest" => ChildRule::Nearest,
                        _ => return Err(bad()),
                    }
                }
                "nodesel" => {
                    self.node_selection = match v {
                        "hybrid" => NodeSelection::Hybrid,
                        "bestbound" => NodeSelection::BestBound,
                        _ => return Err(bad()),
                    }
                }
                _ => return Err(format!("unknown branch option '{k}'")),
            }
        }
        Ok(())
    }
}

/// Measured defaults (warm strong branching). They differ from [`BranchConfig::highs`] where measured:
/// * `child: Nearest` — `kRootSol` ranks sides by inference counts and
///   root-solution distance; on our families it dives the wrong way: medium
///   set sgm 65 → 50 nodes with `Nearest`, multi_knapsack 4420 → 847 (cold
///   strong branching).
/// * `sb_iteration_limit`, `sb_max_candidates`, `sb_lookahead`,
///   `min_reliable` — see the constants.
impl Default for BranchConfig {
    fn default() -> Self {
        Self {
            min_reliable: SB_MIN_RELIABLE,
            sb_base_iterations: HIGHS_SB_BASE_ITERATIONS,
            sb_iteration_limit: SbProbeLimit::Fixed(SB_PROBE_FIXED_ITERATIONS),
            sb_max_candidates: SB_MAX_CANDIDATES,
            sb_lookahead: SB_LOOKAHEAD,
            sb_invest: false,
            child: ChildRule::Nearest,
            node_selection: NodeSelection::Hybrid,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Selection {
    /// Branch on `col` at LP value `value`; `*_bound` are valid lower bounds
    /// for the children; `up_first` is the side to plunge into.
    Branch {
        col: usize,
        value: f64,
        down_bound: f64,
        up_bound: f64,
        up_first: bool,
    },
    /// One side of `col` is empty or cut off: the node reduces to the `up`
    /// (`col ≥ ⌈value⌉`) or down side, with lower bound `bound`.
    Reduce {
        col: usize,
        value: f64,
        up: bool,
        bound: f64,
    },
}

#[derive(Debug, Clone, Default)]
pub struct Brancher {
    pub pseudocost: Pseudocost,
    config: BranchConfig,
    integer: Vec<bool>,
    integrality_tol: f64,
    feastol: f64,
    /// Per-column UP/DOWN flags: side already strong-branched in this dive.
    reliable_at_node: Vec<u8>,
    touched: Vec<usize>,
    root_sol: Vec<f64>,
    subroot_sol: Vec<f64>,
    /// Probes since the last side cut off or bound appreciably moved: the
    /// payoff memory that guards the budget quota.
    sb_dry_probes: u32,
    /// Last investment depth (diagnostic tracing only).
    last_depth: u32,
    /// LP iterations spent in strong branching so far.
    pub sb_iterations: u64,
}

/// Nodes before the investment policy may arm: the depth pays by
/// amortizing probe cost over the remaining tree, so trees that finish
/// inside the warmup horizon cannot pay it back. Measured on the mixed
/// tier: every family whose base tree stays below ≈1.4k
/// nodes regressed time 1.10–1.47× under unconditional depth (their node
/// gains 0.85–1.02× never amortized), while trees in the 6.9k–23k range
/// that armed either won (lights_out 0.93× time at 0.66× nodes) or stayed
/// near-neutral; 2 000 sits just above the regressed cluster, late enough
/// to spare it and early enough that mid-size winners (gap 7.4k,
/// int_production s3 11.9k) still arm inside their warmup phase.
pub const SB_INVEST_ARM_NODES: u64 = 2_000;
/// Effective reliability depth under the investment policy:
/// [`SB_INVEST_MAX_MINREL`] while the tree has shown itself large enough
/// to amortize the probes ([`SB_INVEST_ARM_NODES`]) and strong branching
/// stays inside SCIP's sbiterquot budget `sbiterofs + sbiterquot·(all LP
/// iterations)` (`branch_relpscost.c`), the configured base otherwise —
/// never below `base`, so a user-set `minrel` only ever adds. The offset
/// is what keeps the root warmup's structural share spike (probes dominate
/// the few root LP iterations) from disarming the policy on probe-cheap
/// families.
fn invest_depth(base: u32, nodes: u64, sb_iterations: u64, total_lp_iterations: u64) -> u32 {
    if nodes >= SB_INVEST_ARM_NODES {
        let quota =
            SB_INVEST_QUOTA_OFFSET + (SB_INVEST_SHARE_BASE * total_lp_iterations as f64) as u64;
        if sb_iterations <= quota {
            return base.max(SB_INVEST_MAX_MINREL);
        }
    }
    base
}

impl Brancher {
    pub fn new(integer: &[bool], integrality_tol: f64, feastol: f64, config: BranchConfig) -> Self {
        let mut pseudocost = Pseudocost::new(integer.len());
        pseudocost.set_min_reliable(config.min_reliable);
        Self {
            pseudocost,
            config,
            integer: integer.to_vec(),
            integrality_tol,
            feastol,
            reliable_at_node: vec![0; integer.len()],
            touched: Vec::new(),
            root_sol: Vec::new(),
            subroot_sol: Vec::new(),
            sb_dry_probes: 0,
            last_depth: 0,
            sb_iterations: 0,
        }
    }

    /// A new dive starts: strong-branching results of the previous one no
    /// longer stand in for reliability (HighsSearch.cpp:1871).
    pub fn new_dive(&mut self) {
        for &j in &self.touched {
            self.reliable_at_node[j] = 0;
        }
        self.touched.clear();
    }

    /// A node was taken from the queue: the sub-root solution resets
    /// (`installNode`, HighsSearch.cpp:870).
    pub fn new_subtree(&mut self) {
        self.subroot_sol.clear();
    }

    pub fn set_root_solution(&mut self, x: &[f64]) {
        self.root_sol = x.to_vec();
    }

    /// Branching decision for a node whose LP failed (HiGHS falls back the
    /// same way, HighsSearch.cpp:1347–1418): among integer columns not yet
    /// fixed, the best pseudocost score at the domain midpoint `point`
    /// (`⌊point⌋`/`⌈point⌉` split the domain). Plunges up first when the
    /// column's cost would worsen the objective and an incumbent exists,
    /// else toward more average inferences, ties down. Returns `(col, point,
    /// up_first)`; `None` when every integer column is fixed.
    ///
    /// Not ported: the 1e6 degeneracy factor HiGHS sets for this pick and
    /// the least-recent-domain-change tie-break.
    pub fn select_without_lp(
        &self,
        lower: &[f64],
        upper: &[f64],
        cost: &[f64],
        have_incumbent: bool,
    ) -> Option<(usize, f64, bool)> {
        let ps = &self.pseudocost;
        let mut best: Option<(f64, usize, f64)> = None;
        for j in (0..self.integer.len()).filter(|&j| self.integer[j]) {
            let (lo, hi) = (lower[j], upper[j]);
            if hi - lo < 0.5 {
                continue;
            }
            let point = match (lo.is_finite(), hi.is_finite()) {
                (true, true) => (0.5 * (lo + hi + 0.5)).floor() + 0.5,
                (true, false) => lo + 0.5,
                (false, true) => hi - 0.5,
                (false, false) => 0.5,
            };
            let score = ps.score_frac(j, point);
            if best.is_none_or(|(s, _, _)| score > s) {
                best = Some((score, j, point));
            }
        }
        let (_, col, point) = best?;
        let up = if cost[col].abs() > self.feastol && have_incumbent {
            cost[col] > 0.0
        } else {
            ps.avg_inferences_up(col) > ps.avg_inferences_down(col) + self.feastol
        };
        Some((col, point, up))
    }

    fn mark(&mut self, col: usize, flag: u8) {
        if self.reliable_at_node[col] == 0 {
            self.touched.push(col);
        }
        self.reliable_at_node[col] |= flag;
    }

    /// Integer columns with a fractional LP value, in column order.
    pub fn fractional(&self, x: &[f64]) -> Vec<(usize, f64)> {
        x.iter()
            .enumerate()
            .filter(|&(j, &v)| self.integer[j] && (v - v.round()).abs() > self.integrality_tol)
            .map(|(j, &v)| (j, v))
            .collect()
    }

    /// Chooses the branching column at a node whose LP (objective
    /// `lp_objective`, point `x`) has the fractional candidates `frac`.
    /// `total_lp_iterations` counts every LP iteration so far, strong
    /// branching included.
    #[allow(clippy::too_many_arguments)]
    pub fn select<H: SbHost>(
        &mut self,
        host: &mut H,
        lower: &[f64],
        upper: &[f64],
        lp_objective: f64,
        x: &[f64],
        frac: &[(usize, f64)],
        total_lp_iterations: u64,
    ) -> Selection {
        debug_assert!(!frac.is_empty());
        // Investment depth: the effective reliability
        // threshold grows with the SB share's remaining headroom under the
        // sbiterquot quota, on trees large enough to amortize the probes.
        // It stays set after the call — the node queue's best-estimate
        // blending weights by it too — and is recomputed from the
        // configured base at every node, so it never compounds. With the
        // policy off, min_reliable is left exactly as it was.
        if self.config.sb_invest {
            let grown = invest_depth(
                self.config.min_reliable,
                host.sb_nodes(),
                self.sb_iterations,
                total_lp_iterations,
            );
            if std::env::var_os("RUSOLVE_MIP_LOG").is_some() && grown != self.last_depth {
                eprintln!(
                    "sb invest: depth -> {grown} at sb={} of {} lp iters",
                    self.sb_iterations, total_lp_iterations
                );
                self.last_depth = grown;
            }
            self.pseudocost.set_min_reliable(grown);
        }
        let min_rel = self.pseudocost.min_reliable();
        let mut max_sb = 0;
        if min_rel > 0 {
            let sb = self.sb_iterations;
            max_sb = self.config.sb_base_iterations + (total_lp_iterations.saturating_sub(sb) >> 1);
            // Payoff guard: probes that keep cutting off sides or moving
            // the bound are load-bearing — trimming them on neos-860300
            // multiplies the tree 4–8× — so the budget never drops below
            // HiGHS's formula while they flow. Only a dry window of probes
            // (no cutoff, no appreciable move) tightens the ceiling toward
            // SCIP's 50% quota (`sbiterquot` 0.5, `sbiterofs` offset,
            // branch_relpscost.c).
            if self.sb_dry_probes >= SB_DRY_PROBE_LIMIT {
                let quota = SB_QUOTA_OFFSET + total_lp_iterations / 2;
                max_sb = max_sb.min(quota);
            }
            if sb > max_sb {
                self.pseudocost.set_min_reliable(0);
            } else if sb > max_sb >> 1 {
                let half = max_sb >> 1;
                let ratio = (sb - half) as f64 / (max_sb - half) as f64;
                let reduced = (f64::from(min_rel) - ratio * f64::from(min_rel - 1)) as u32;
                self.pseudocost.set_min_reliable(min_rel.min(reduced));
            }
        }
        let sel = self.select_inner(host, lower, upper, lp_objective, frac, max_sb);
        self.pseudocost.set_min_reliable(min_rel);
        match sel {
            Selection::Branch {
                col,
                value,
                down_bound,
                up_bound,
                ..
            } => Selection::Branch {
                col,
                value,
                down_bound,
                up_bound,
                up_first: match self.config.child {
                    ChildRule::RootSol => self.up_first(col, value, lp_objective, x),
                    ChildRule::Nearest => value - value.floor() >= 0.5,
                },
            },
            r => r,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn select_inner<H: SbHost>(
        &mut self,
        host: &mut H,
        lower: &[f64],
        upper: &[f64],
        lp_objective: f64,
        frac: &[(usize, f64)],
        max_sb: u64,
    ) -> Selection {
        let nf = frac.len();
        let mut up_score = vec![f64::INFINITY; nf];
        let mut down_score = vec![f64::INFINITY; nf];
        let mut up_rel = vec![false; nf];
        let mut down_rel = vec![false; nf];
        let mut up_bound = vec![lp_objective; nf];
        let mut down_bound = vec![lp_objective; nf];
        let ps = &self.pseudocost;
        for (k, &(col, v)) in frac.iter().enumerate() {
            let flags = self.reliable_at_node[col];
            let reliable = ps.is_reliable(col);
            if reliable || flags & UP != 0 {
                up_score[k] = ps.cost_up(col, v);
                up_rel[k] = true;
            }
            if reliable || flags & DOWN != 0 {
                down_score[k] = ps.cost_down(col, v);
                down_rel[k] = true;
            }
        }
        let mut min_score = self.feastol;
        let lookahead = self.config.sb_lookahead;
        let (mut known_best, mut stale) = (None, 0u32);
        let max_cands = self.config.sb_max_candidates as usize;
        let mut probed = vec![false; nf];
        let mut num_probed = 0usize;
        loop {
            let mut must_stop = self.sb_iterations >= max_sb
                || host.sb_stop()
                || (lookahead > 0 && stale >= lookahead);
            let mut cand = self.best_score(
                frac,
                &up_score,
                &down_score,
                &mut up_rel,
                &mut down_rel,
                &mut min_score,
                must_stop,
            );
            // Candidate cap: a candidate probed on one side may finish the
            // other; a new one ends strong branching.
            if !must_stop && max_cands > 0 && num_probed >= max_cands && !probed[cand] {
                must_stop = true;
                cand = self.best_score(
                    frac,
                    &up_score,
                    &down_score,
                    &mut up_rel,
                    &mut down_rel,
                    &mut min_score,
                    true,
                );
            }
            let (col, v) = frac[cand];
            if lookahead > 0 {
                let ps = &self.pseudocost;
                let best = (0..nf)
                    .filter(|&k| up_score[k].is_finite() && down_score[k].is_finite())
                    .map(|k| (ps.score(frac[k].0, up_score[k], down_score[k]), k))
                    .max_by(|a, b| a.0.total_cmp(&b.0))
                    .map(|(_, k)| k);
                if best.is_some() && best == known_best {
                    stale += 1;
                } else {
                    stale = 0;
                }
                known_best = best;
            }
            if (up_rel[cand] && down_rel[cand]) || must_stop {
                return Selection::Branch {
                    col,
                    value: v,
                    down_bound: down_bound[cand],
                    up_bound: up_bound[cand],
                    up_first: false,
                };
            }
            // Evaluate the side that looks worse first: it is the likelier
            // one to be cut off (HighsSearch.cpp:690).
            let ps = &self.pseudocost;
            let up = !(!down_rel[cand]
                && (up_rel[cand]
                    || (down_score[cand], ps.avg_inferences_down(col))
                        >= (up_score[cand], ps.avg_inferences_up(col))));
            let other_bound = if up { down_bound[cand] } else { up_bound[cand] };
            let reduce = Selection::Reduce {
                col,
                value: v,
                up: !up,
                bound: other_bound,
            };
            let mut cl = lower.to_vec();
            let mut cu = upper.to_vec();
            if up {
                cl[col] = v.ceil();
            } else {
                cu[col] = v.floor();
            }
            if !host.sb_propagate(&mut cl, &mut cu, col) {
                self.pseudocost.add_cutoff_observation(col, up);
                self.sb_dry_probes = 0;
                return reduce;
            }
            let inferences = (0..cl.len())
                .filter(|&j| j != col && (cl[j] != lower[j] || cu[j] != upper[j]))
                .count();
            self.pseudocost
                .add_inference_observation(col, inferences, up);
            let flag = if up { UP } else { DOWN };
            if !probed[cand] {
                probed[cand] = true;
                num_probed += 1;
            }
            // The child's LP bound, and its point when solved to optimality.
            // The probe counts as dry (quota guard) unless it cuts a side
            // off or moves the bound past the EPSILON noise floor.
            self.sb_dry_probes = self.sb_dry_probes.saturating_add(1);
            let solved = match host.sb_solve(&cl, &cu) {
                Ok(lp) => {
                    self.sb_iterations += lp.iterations;
                    match lp.status {
                        LpStatus::Optimal => Some((lp.objective, Some(lp.x))),
                        // A capped probe's dual bound is a valid child bound
                        // and stands in for the gain (Achterberg 2007 §5.4).
                        LpStatus::IterationLimit if lp.objective.is_finite() => {
                            Some((lp.objective, None))
                        }
                        // Cut off by the objective limit counts as
                        // infeasible (HighsLpRelaxation.cpp:1263).
                        LpStatus::Infeasible | LpStatus::Cutoff => {
                            self.pseudocost.add_cutoff_observation(col, up);
                            self.sb_dry_probes = 0;
                            return reduce;
                        }
                        LpStatus::IterationLimit | LpStatus::Unbounded => None,
                    }
                }
                Err(LpError::IterationLimit) => {
                    self.sb_iterations += host.sb_iteration_limit();
                    None
                }
                Err(_) => None,
            };
            let Some((child_obj, x)) = solved else {
                // LP trouble (or an unbounded child of a bounded node): score
                // the candidate 0 so it is avoided if possible
                // (HighsSearch.cpp:676–685).
                up_score[cand] = 0.0;
                down_score[cand] = 0.0;
                up_rel[cand] = true;
                down_rel[cand] = true;
                self.mark(col, UP | DOWN);
                continue;
            };
            let mut objdelta = (child_obj - lp_objective).max(0.0);
            if objdelta <= EPSILON {
                objdelta = 0.0;
            } else {
                // The probe moved the child's bound: payoff (quota guard).
                self.sb_dry_probes = 0;
            }
            if up {
                up_score[cand] = objdelta;
                up_rel[cand] = true;
                up_bound[cand] = child_obj;
            } else {
                down_score[cand] = objdelta;
                down_rel[cand] = true;
                down_bound[cand] = child_obj;
            }
            self.mark(col, flag);
            let delta = if up { v.ceil() - v } else { v.floor() - v };
            self.pseudocost.add_observation(col, delta, objdelta);
            if let Some(x) = x {
                self.analyze_solution(
                    host,
                    &cl,
                    &cu,
                    &x,
                    objdelta,
                    col,
                    frac,
                    &mut up_score,
                    &mut down_score,
                );
                if self.fractional(&x).is_empty() {
                    host.sb_integral(x);
                }
            }
            if child_obj > host.upper_limit() {
                self.sb_dry_probes = 0;
                return reduce;
            }
        }
    }

    /// `selectBestScore` (HighsSearch.cpp:337–383): candidates whose score
    /// is still unknown rank first (∞) until the final selection, where the
    /// pseudocost estimate stands in. Scores at or below `min_score` (which
    /// tracks 1e-3 × the best known two-sided gain) count as reliable. Among
    /// unknown candidates the pseudocost estimate decides (HiGHS: column
    /// order), so a capped strong-branching loop probes the most promising
    /// ones (Achterberg, Koch & Martin 2005 §3).
    #[allow(clippy::too_many_arguments)]
    fn best_score(
        &self,
        frac: &[(usize, f64)],
        up_score: &[f64],
        down_score: &[f64],
        up_rel: &mut [bool],
        down_rel: &mut [bool],
        min_score: &mut f64,
        final_selection: bool,
    ) -> usize {
        let ps = &self.pseudocost;
        let old_min = *min_score;
        let mut best = 0;
        let (mut best_score, mut best_tie) = (-1.0, -1.0);
        for (k, &(col, v)) in frac.iter().enumerate() {
            if up_score[k] <= old_min {
                up_rel[k] = true;
            }
            if down_score[k] <= old_min {
                down_rel[k] = true;
            }
            let s = 1e-3
                * (if up_rel[k] { up_score[k] } else { 0.0 }).min(if down_rel[k] {
                    down_score[k]
                } else {
                    0.0
                });
            *min_score = min_score.max(s);
            let score = if up_score[k] <= old_min || down_score[k] <= old_min {
                ps.score(col, up_score[k].min(old_min), down_score[k].min(old_min))
            } else if up_score[k] == f64::INFINITY || down_score[k] == f64::INFINITY {
                if final_selection {
                    ps.score_frac(col, v)
                } else {
                    f64::INFINITY
                }
            } else {
                ps.score(col, up_score[k], down_score[k])
            };
            let tie = if score == f64::INFINITY {
                ps.score_frac(col, v)
            } else {
                0.0
            };
            if score > best_score || (score == best_score && tie > best_tie) {
                best_score = score;
                best_tie = tie;
                best = k;
            }
        }
        best
    }

    /// Scores other candidates from a strong-branch LP solution that already
    /// lies on one of their branch sides, when propagating that side does
    /// not contradict the solution (`analyzeSolution`, HighsSearch.cpp:414).
    #[allow(clippy::too_many_arguments)]
    fn analyze_solution<H: SbHost>(
        &mut self,
        host: &mut H,
        lower: &[f64],
        upper: &[f64],
        sol: &[f64],
        objdelta: f64,
        branched: usize,
        frac: &[(usize, f64)],
        up_score: &mut [f64],
        down_score: &mut [f64],
    ) {
        let tol = self.feastol;
        for (k, &(col, v)) in frac.iter().enumerate() {
            if col == branched {
                continue;
            }
            let (down_val, up_val) = (v.floor(), v.ceil());
            let up = if sol[col] <= down_val + tol {
                false
            } else if sol[col] >= up_val - tol {
                true
            } else {
                continue;
            };
            let needs_check = if up {
                lower[col] < up_val - tol
            } else {
                upper[col] > down_val + tol
            };
            if needs_check && !self.side_consistent(host, lower, upper, sol, col, up, v) {
                continue;
            }
            if objdelta <= tol {
                let delta = if up { up_val - v } else { down_val - v };
                self.pseudocost.add_observation(col, delta, objdelta);
                self.mark(col, if up { UP } else { DOWN });
            }
            if up {
                up_score[k] = up_score[k].min(objdelta);
            } else {
                down_score[k] = down_score[k].min(objdelta);
            }
        }
    }

    /// Applies the `up`/down branch of `col` to the strong-branch child box,
    /// propagates, and checks that `sol` satisfies every derived bound.
    #[allow(clippy::too_many_arguments)]
    fn side_consistent<H: SbHost>(
        &self,
        host: &mut H,
        lower: &[f64],
        upper: &[f64],
        sol: &[f64],
        col: usize,
        up: bool,
        v: f64,
    ) -> bool {
        let tol = self.feastol;
        let mut l = lower.to_vec();
        let mut u = upper.to_vec();
        if up {
            l[col] = v.ceil();
        } else {
            u[col] = v.floor();
        }
        if !host.sb_propagate(&mut l, &mut u, col) {
            return false;
        }
        (0..l.len()).all(|j| j == col || (l[j] <= sol[j] + tol && u[j] >= sol[j] - tol))
    }

    /// `kRootSol` child selection (HighsSearch.cpp:1166): prefer the side
    /// with more average inferences, scaled by how far the value moved away
    /// from the root or sub-root LP solution in that direction.
    fn up_first(&mut self, col: usize, value: f64, lp_objective: f64, x: &[f64]) -> bool {
        let ps = &self.pseudocost;
        let mut down_prio = ps.avg_inferences_down(col) + EPSILON;
        let mut up_prio = ps.avg_inferences_up(col) + EPSILON;
        let (down_val, up_val) = (value.floor(), value.ceil());
        let reference = if !self.subroot_sol.is_empty() {
            Some(self.subroot_sol[col])
        } else {
            if lp_objective.is_finite() {
                self.subroot_sol = x.to_vec();
            }
            self.root_sol.get(col).copied()
        };
        if let Some(r) = reference {
            let r = r.clamp(down_val, up_val);
            up_prio *= 1.0 + (value - r);
            down_prio *= 1.0 + (r - value);
        }
        up_prio + EPSILON >= down_prio
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Host over a tiny LP: min Σ c_j x_j with the single row Σ x_j ≥ r,
    /// solved greedily (all c_j > 0); counts solves.
    struct Toy {
        cost: Vec<f64>,
        r: f64,
        solves: usize,
        limit: f64,
        found: Vec<Vec<f64>>,
    }

    impl SbHost for Toy {
        fn sb_solve(&mut self, lower: &[f64], upper: &[f64]) -> Result<LpOutcome, LpError> {
            self.solves += 1;
            let mut x = lower.to_vec();
            let mut need = self.r - x.iter().sum::<f64>();
            let mut order: Vec<usize> = (0..x.len()).collect();
            order.sort_by(|&a, &b| self.cost[a].total_cmp(&self.cost[b]));
            for j in order {
                if need <= 0.0 {
                    break;
                }
                let add = (upper[j] - x[j]).min(need);
                x[j] += add;
                need -= add;
            }
            if need > 1e-9 {
                return Ok(LpOutcome {
                    status: LpStatus::Infeasible,
                    objective: f64::INFINITY,
                    x: vec![],
                    iterations: 1,
                    basis: None,
                });
            }
            let objective = x.iter().zip(&self.cost).map(|(a, b)| a * b).sum();
            Ok(LpOutcome {
                status: LpStatus::Optimal,
                objective,
                x,
                iterations: 1,
                basis: None,
            })
        }
        fn sb_iteration_limit(&self) -> u64 {
            1
        }
        fn sb_propagate(&mut self, _: &mut [f64], _: &mut [f64], _: usize) -> bool {
            true
        }
        fn upper_limit(&self) -> f64 {
            self.limit
        }
        fn sb_integral(&mut self, x: Vec<f64>) {
            self.found.push(x);
        }
        fn sb_stop(&self) -> bool {
            false
        }
        fn sb_nodes(&self) -> u64 {
            SB_INVEST_ARM_NODES
        }
    }

    #[test]
    fn strong_branching_scores_unreliable_candidates() {
        // min x0 + 2x1 + 3x2, x0 + x1 + x2 ≥ 1.5, x ∈ [0,1]: LP x = (1, .5, 0).
        let mut host = Toy {
            cost: vec![1.0, 2.0, 3.0],
            r: 1.5,
            solves: 0,
            limit: f64::INFINITY,
            found: vec![],
        };
        let (lo, hi) = (vec![0.0; 3], vec![1.0; 3]);
        let lp = host.sb_solve(&lo, &hi).unwrap();
        host.solves = 0;
        let mut b = Brancher::new(&[true; 3], 1e-6, 1e-6, BranchConfig::highs());
        b.set_root_solution(&lp.x);
        let frac = b.fractional(&lp.x);
        assert_eq!(frac, vec![(1, 0.5)]);
        let sel = b.select(&mut host, &lo, &hi, lp.objective, &lp.x, &frac, 0);
        // Down: x1 = 0 ⇒ x = (1, 0, .5), obj 2.5; up: x1 = 1 ⇒ x = (.5, 1, 0),
        // obj 2.5. Both +0.5 over the LP's 2.
        match sel {
            Selection::Branch {
                col,
                down_bound,
                up_bound,
                ..
            } => {
                assert_eq!(col, 1);
                assert!((down_bound - 2.5).abs() < 1e-12);
                assert!((up_bound - 2.5).abs() < 1e-12);
            }
            other => panic!("{other:?}"),
        }
        assert_eq!(host.solves, 2);
        assert!(host.found.is_empty());
        assert_eq!(b.pseudocost.num_observations_up(1), 1);
        assert_eq!(b.pseudocost.num_observations_down(1), 1);
        // Same dive: both sides are now reliable-at-node, no new solves.
        b.select(&mut host, &lo, &hi, lp.objective, &lp.x, &frac, 0);
        assert_eq!(host.solves, 2);
    }

    #[test]
    fn cut_off_side_reduces_the_node() {
        let mut host = Toy {
            cost: vec![1.0, 2.0, 3.0],
            r: 1.5,
            solves: 0,
            limit: 2.4,
            found: vec![],
        };
        let (lo, hi) = (vec![0.0; 3], vec![1.0; 3]);
        let lp = host.sb_solve(&lo, &hi).unwrap();
        let mut b = Brancher::new(&[true; 3], 1e-6, 1e-6, BranchConfig::highs());
        let frac = b.fractional(&lp.x);
        // The down side (evaluated first) costs 2.5 > 2.4 ⇒ only up survives.
        let sel = b.select(&mut host, &lo, &hi, lp.objective, &lp.x, &frac, 0);
        match sel {
            Selection::Reduce { col, up, .. } => assert_eq!((col, up), (1, true)),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn zero_budget_falls_back_to_pseudocosts() {
        let mut host = Toy {
            cost: vec![1.0, 2.0, 3.0],
            r: 1.5,
            solves: 0,
            limit: f64::INFINITY,
            found: vec![],
        };
        let (lo, hi) = (vec![0.0; 3], vec![1.0; 3]);
        let lp = host.sb_solve(&lo, &hi).unwrap();
        host.solves = 0;
        let mut b = Brancher::new(&[true; 3], 1e-6, 1e-6, BranchConfig::highs());
        b.pseudocost.set_min_reliable(0);
        let frac = b.fractional(&lp.x);
        assert!(matches!(
            b.select(&mut host, &lo, &hi, lp.objective, &lp.x, &frac, 0),
            Selection::Branch { col: 1, .. }
        ));
        assert_eq!(host.solves, 0);
    }

    /// Host where branching on column j moves the objective by `w[j]` on
    /// either side and leaves every other column at ½.
    struct Stub {
        w: Vec<f64>,
        solves: usize,
    }

    impl SbHost for Stub {
        fn sb_solve(&mut self, lower: &[f64], upper: &[f64]) -> Result<LpOutcome, LpError> {
            self.solves += 1;
            let mut x = vec![0.5; self.w.len()];
            let mut objective = 0.0;
            for j in 0..x.len() {
                if lower[j] > 0.0 || upper[j] < 1.0 {
                    x[j] = lower[j];
                    objective += self.w[j];
                }
            }
            Ok(LpOutcome {
                status: LpStatus::Optimal,
                objective,
                x,
                iterations: 1,
                basis: None,
            })
        }
        fn sb_iteration_limit(&self) -> u64 {
            1
        }
        fn sb_propagate(&mut self, _: &mut [f64], _: &mut [f64], _: usize) -> bool {
            true
        }
        fn upper_limit(&self) -> f64 {
            f64::INFINITY
        }
        fn sb_integral(&mut self, _: Vec<f64>) {}
        fn sb_stop(&self) -> bool {
            false
        }
        fn sb_nodes(&self) -> u64 {
            SB_INVEST_ARM_NODES
        }
    }

    #[test]
    fn lookahead_stops_strong_branching_early() {
        let n = 10;
        let w: Vec<f64> = (0..n).map(|j| 10.0 - j as f64).collect();
        let frac: Vec<(usize, f64)> = (0..n).map(|j| (j, 0.5)).collect();
        let (lo, hi, x) = (vec![0.0; n], vec![1.0; n], vec![0.5; n]);
        let run = |lookahead| {
            let config = BranchConfig {
                sb_lookahead: lookahead,
                ..BranchConfig::highs()
            };
            let mut b = Brancher::new(&vec![true; n], 1e-6, 1e-6, config);
            let mut host = Stub {
                w: w.clone(),
                solves: 0,
            };
            let sel = b.select(&mut host, &lo, &hi, 0.0, &x, &frac, 0);
            (sel, host.solves)
        };
        // Without lookahead every unreliable candidate is strong-branched.
        let (sel, solves) = run(0);
        assert!(matches!(sel, Selection::Branch { col: 0, .. }));
        assert_eq!(solves, 2 * n);
        // Column 0 is complete and best after 2 LPs; two more LPs that do
        // not displace it end strong branching.
        let (sel, solves) = run(2);
        assert!(matches!(sel, Selection::Branch { col: 0, .. }));
        assert_eq!(solves, 5);
    }

    #[test]
    fn adaptive_probe_cap_scales_with_node_lps() {
        use SbProbeLimit::*;
        // The cap scales with the node LPs a probe refines, floored at a
        // degenerate dual step's worth, before the first node LP via the
        // root's iterations.
        assert_eq!(Adaptive.resolve(0, 7, 10_000), 14);
        assert_eq!(Adaptive.resolve(0, 0, 10_000), SB_PROBE_MIN_ITERATIONS);
        assert_eq!(Adaptive.resolve(5, 4_000, 10_000), SB_PROBE_MIN_ITERATIONS);
        assert_eq!(Adaptive.resolve(120, 4_000, 10_000), 240);
        assert_eq!(Adaptive.resolve(400, 10, 10_000), SB_PROBE_MAX_ITERATIONS);
        assert_eq!(Fixed(100).resolve(120, 4_000, 10_000), 100);
        assert_eq!(Tree.resolve(120, 4_000, 10_000), 10_000);
    }

    #[test]
    fn sb_budget_limits_the_share_of_tree_iterations() {
        // base 0 + half of the 2n non-SB iterations = n: one iteration per
        // probe stops strong branching after n of the 2n probes it would
        // otherwise run (SCIP sbiterquot pattern).
        let n = 10;
        let w: Vec<f64> = (0..n).map(|j| 10.0 - j as f64).collect();
        let frac: Vec<(usize, f64)> = (0..n).map(|j| (j, 0.5)).collect();
        let (lo, hi, x) = (vec![0.0; n], vec![1.0; n], vec![0.5; n]);
        let config = BranchConfig {
            sb_base_iterations: 0,
            sb_lookahead: 0,
            ..BranchConfig::default()
        };
        let mut b = Brancher::new(&vec![true; n], 1e-6, 1e-6, config);
        let mut host = Stub { w, solves: 0 };
        assert!(matches!(
            b.select(&mut host, &lo, &hi, 0.0, &x, &frac, 2 * n as u64),
            Selection::Branch { .. }
        ));
        assert_eq!(host.solves, n);
    }

    #[test]
    fn dry_probes_tighten_the_budget_only_without_payoff() {
        let n = 20;
        let (lo, hi, x) = (vec![0.0; n], vec![1.0; n], vec![0.5; n]);
        let frac: Vec<(usize, f64)> = (0..n).map(|j| (j, 0.5)).collect();
        let config = BranchConfig {
            sb_lookahead: 0,
            ..BranchConfig::default()
        };
        // A dry round — every probe's bound stays at the node LP's — arms
        // the quota (at least `SB_DRY_PROBE_LIMIT` probes, no payoff); with
        // strong branching past the 50% quota (120k of 200k iterations) but
        // inside HiGHS's formula (140k), the next node probes nothing.
        let mut b = Brancher::new(&vec![true; n], 1e-6, 1e-6, config);
        let mut host = Stub {
            w: vec![0.0; n],
            solves: 0,
        };
        b.select(&mut host, &lo, &hi, 0.0, &x, &frac, 200_000);
        let dry = host.solves;
        assert!(dry >= SB_DRY_PROBE_LIMIT as usize);
        b.sb_iterations = 120_000;
        b.select(&mut host, &lo, &hi, 0.0, &x, &frac, 200_000);
        assert_eq!(host.solves, dry);

        // The same budget state with payoff flowing — a probe whose child
        // the cutoff removes — keeps the loose budget and probes on.
        let (lo3, hi3) = (vec![0.0; 3], vec![1.0; 3]);
        let mut toy = Toy {
            cost: vec![1.0, 2.0, 3.0],
            r: 1.5,
            solves: 0,
            limit: 2.4,
            found: vec![],
        };
        let lp = toy.sb_solve(&lo3, &hi3).unwrap();
        toy.solves = 0;
        let mut b2 = Brancher::new(&[true; 3], 1e-6, 1e-6, config);
        let frac3 = b2.fractional(&lp.x);
        assert!(matches!(
            b2.select(&mut toy, &lo3, &hi3, lp.objective, &lp.x, &frac3, 200_000),
            Selection::Reduce { .. }
        ));
        assert_eq!(toy.solves, 1);
        b2.sb_iterations = 120_000;
        b2.select(&mut toy, &lo3, &hi3, lp.objective, &lp.x, &frac3, 200_000);
        assert_eq!(toy.solves, 2);
    }

    #[test]
    fn invest_depth_arms_inside_the_sbiterquot_budget() {
        // Armed on a large tree while sb iterations fit SCIP's sbiterofs +
        // 50% quota — including a root warmup spike like sb=559 of 778
        // total (quota 1389) — and disarmed past it, e.g. a probe-expensive
        // family at sb=29k of 40k (quota 21k). Below the node threshold the
        // policy never arms: small trees cannot amortize the probes. A
        // user-set larger base is never lowered; base 8 grows to the
        // investment maximum.
        assert_eq!(
            invest_depth(8, SB_INVEST_ARM_NODES, 0, 1),
            SB_INVEST_MAX_MINREL
        );
        assert_eq!(
            invest_depth(8, SB_INVEST_ARM_NODES, 559, 778),
            SB_INVEST_MAX_MINREL
        );
        assert_eq!(
            invest_depth(8, SB_INVEST_ARM_NODES, 1_389, 778),
            SB_INVEST_MAX_MINREL
        );
        assert_eq!(invest_depth(8, SB_INVEST_ARM_NODES, 29_000, 40_000), 8);
        assert_eq!(invest_depth(8, 0, 0, 1), 8);
        assert_eq!(invest_depth(40, SB_INVEST_ARM_NODES, 0, 1), 40);
    }

    #[test]
    fn invest_depth_extends_the_warmup_and_persists_between_nodes() {
        // With the policy on and zero SB share, a column needs the grown
        // depth of observations per side before its pseudo-cost is trusted;
        // the grown depth also persists after the call (the node queue's
        // best-estimate blending weights by it), while the budget decay is
        // still reverted. With the policy off nothing changes.
        let n = 10;
        let w: Vec<f64> = (0..n).map(|j| 10.0 - j as f64).collect();
        let frac: Vec<(usize, f64)> = (0..n).map(|j| (j, 0.5)).collect();
        let (lo, hi, x) = (vec![0.0; n], vec![1.0; n], vec![0.5; n]);
        let run = |invest| {
            let config = BranchConfig {
                sb_invest: invest,
                ..BranchConfig::default()
            };
            let mut b = Brancher::new(&vec![true; n], 1e-6, 1e-6, config);
            let mut host = Stub {
                w: w.clone(),
                solves: 0,
            };
            let sel = b.select(&mut host, &lo, &hi, 0.0, &x, &frac, 10_000 * n as u64);
            (sel, host.solves, b.pseudocost.min_reliable())
        };
        let (sel0, solves0, rel0) = run(false);
        let (sel1, solves1, rel1) = run(true);
        assert!(matches!(sel0, Selection::Branch { .. }));
        assert!(matches!(sel1, Selection::Branch { .. }));
        assert_eq!(rel0, BranchConfig::default().min_reliable);
        assert_eq!(rel1, SB_INVEST_MAX_MINREL);
        // Zero share at entry: the first node invests, so the warmup of the
        // best candidate pair needs more probes than the base depth's.
        assert!(solves1 >= solves0, "{solves1} vs {solves0}");
    }

    #[test]
    fn branching_without_lp_splits_an_unfixed_integer_domain() {
        let mut b = Brancher::new(&[true, true, false], 1e-6, 1e-6, BranchConfig::default());
        let cost = [1.0, 0.0, 0.0];
        // Column 1 is fixed, column 2 continuous: only column 0 qualifies.
        let (lo, hi) = ([0.0, 0.0, 0.0], [3.0, 0.0, 5.0]);
        let (col, point, up) = b.select_without_lp(&lo, &hi, &cost, true).unwrap();
        assert_eq!(col, 0);
        assert_eq!(point, 1.5);
        // Positive cost and an incumbent: the worsening (up) side first.
        assert!(up);
        // No incumbent and no inference data: down first.
        assert!(!b.select_without_lp(&lo, &hi, &cost, false).unwrap().2);
        b.pseudocost.add_inference_observation(0, 3, true);
        assert!(b.select_without_lp(&lo, &hi, &cost, false).unwrap().2);
        assert!(b.select_without_lp(&lo, &lo, &cost, true).is_none());
    }
}
