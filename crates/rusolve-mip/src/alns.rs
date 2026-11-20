//! Adaptive large-neighbourhood search (ALNS): a bandit chooses which
//! sub-MIP neighbourhood to search next and learns from what each one buys.
//!
//! Every call picks one [`Arm`] (RENS, RINS, mutation, crossover), fixes
//! integer columns to that arm's per-arm adaptive target rate, solves the
//! restriction with [`submip::solve_restricted`] under a node and time
//! budget, and rewards the arm by what it found per node spent. The
//! neighbourhood sources:
//!
//! * **RENS**: fix integral LP values; when too few, round the
//!   nearly integral ones; box the rest to their two roundings.
//! * **RINS**: fix where the incumbent and the node LP agree.
//! * **Mutation**: fix a random subset of columns to the incumbent.
//! * **Crossover** (solution merging): fix where two pool solutions,
//!   drawn with a bias towards the better ones, agree.
//!
//! RINS, mutation and crossover top up a short fixing from their reference
//! solution, and every arm trims an oversized one, so that each call sits
//! near the arm's target rate.
//!
//! Provenance:
//!
//! * G. Hendel, "Adaptive large neighborhood search for mixed integer
//!   programming", Math. Prog. Comp. 14:185–221 (2022): the bandit over
//!   neighbourhoods, the reward, and the per-neighbourhood fixing rate.
//! * SCIP `heur_alns.c` and `bandit_exp3ix.c` (Apache-2.0, © 2002–2026
//!   Zuse Institute Berlin; read as a reference text and reimplemented):
//!   - fixing-rate adaptation (`updateFixingRate`): start at the maximum
//!     0.9, minimum 0.3, step 0.2 decaying by 0.75 per call down to 0.01;
//!     the rate goes down after an optimal or infeasible sub-MIP and up
//!     after a node limit without a solution;
//!   - target fixings within ±0.1 of the rate (`neighborhoodFixVariables`,
//!     `DEFAULT_FIXTOL`/`DEFAULT_UNFIXTOL`), and no solve below half the
//!     target (`heurExecAlns`, after the sub-SCIP's presolve);
//!   - the node budget (`determineLimits`): `0.1 · (found + 1)/(calls + 1)
//!     · parent nodes + 500 − used − 100 · calls`, at most 5000, scaled by
//!     `(K − uninitialised + 1)/(K + 1)`, and no call while it is below the
//!     target node count. That count starts at 50 and grows by ×1.05 + 1
//!     after each node limit without a solution (`updateTargetNodeLimit`);
//!   - wait 25 nodes after a new incumbent (`DEFAULT_WAITINGNODES`);
//!   - the reward (`getReward`): with a new incumbent,
//!     `0.5 + 0.5 · (0.8 + 0.2 · closed gap)/(effort + 1)`, with effort
//!     `nodes/100 · (1 − fixing rate)`; without one,
//!     `0.5 · max(0, 1 − nodes · (1 − fixing rate)/target nodes)`;
//!   - crossover draws 2 solutions biased to the better ones and fixes
//!     where they agree (`varFixingsCrossover`); mutation fixes
//!     `⌊rate · n⌋ + 1` random columns to the incumbent
//!     (`varFixingsMutation`); RENS rounds the most nearly integral
//!     fractional columns when the integral ones fall short
//!     (`varFixingsRens`);
//!   - EXP3-IX: `γ_t = sqrt(ln K / (4 t K))`, `η = 2γ`, loss estimate
//!     `(1 − r)/(p + γ)`, weight `·= exp(−η · loss)`.
//!     G. Neu, "Explore no more: improved high-probability regret bounds
//!     for non-stochastic bandits", NeurIPS 2015, introduced EXP3-IX.
//!
//! Deviations (deliberate):
//!
//! * **Time share (ours).** A sub-MIP here pays its own presolve and LP
//!   engine setup. On second-scale instances that made `--submip` a net
//!   loss even when trees shrank (measured on our benchmark set).
//!   So, on top of SCIP's node budget, a call runs only while the ALNS
//!   wall time plus the expected cost of one more call stays within
//!   [`TIME_SHARE`] of the elapsed solve time. The expected cost is the
//!   mean past call, or half the root time before the first call. Each
//!   call's time limit is the share left. A parent solve that finishes
//!   within a few root times never pays for a sub-MIP.
//! * Unavailable arms (no incumbent, fewer than two pool solutions) are
//!   left out of the draw and the distribution is renormalised over the
//!   rest; the loss estimate uses the probability the arm was actually
//!   drawn with. An arm whose fixing cannot run is penalised with reward 0
//!   and the draw repeats, as in SCIP.
//! * Top-up fixings take columns in seeded random order. SCIP orders them
//!   by graph distance, reduced cost and pseudocost (`alnsFixMoreVariables`).
//!   // TODO(perf): the SCIP priorities need the variable graph and
//!   // pseudocosts here; measure random order first.
//! * The sub-MIP stops at an absolute node limit (SCIP: stall nodes) and
//!   has no solution limit; DINS, local branching, proximity and the
//!   zero-objective neighbourhood are not implemented (they need extra rows
//!   or objectives in the restriction).

use crate::bnb::{MipConfig, MipStatus};
use crate::problem::MipProblem;
use crate::submip::{self, Fixing, SubMipConfig, SubMipOutcome, SubMipReport};
use rusolve_lp::rng::XorShift64;
use std::time::Instant;

pub const ARMS: usize = 4;

/// Share of the elapsed solve time that ALNS may spend (ours; the time
/// analogue of SCIP's 0.1 node quotient).
pub const TIME_SHARE: f64 = 0.1;

const MIN_FIXING_RATE: f64 = 0.3;
const MAX_FIXING_RATE: f64 = 0.9;
const FIXING_RATE_START_INC: f64 = 0.2;
const FIXING_RATE_DECAY: f64 = 0.75;
const LRATE_MIN: f64 = 0.01;
const FIX_TOL: f64 = 0.1;
const UNFIX_TOL: f64 = 0.1;
const NODES_QUOT: f64 = 0.1;
const NODES_OFFSET: f64 = 500.0;
const MIN_NODES: u64 = 50;
const MAX_NODES: u64 = 5000;
const TARGET_NODE_FACTOR: f64 = 1.05;
const WAITING_NODES: u64 = 25;
const REWARD_CONTROL: f64 = 0.8;
const REWARD_BASELINE: f64 = 0.5;
const POOL_SIZE: usize = 10;
const CROSSOVER_SOLS: usize = 2;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arm {
    Rens,
    Rins,
    Mutation,
    Crossover,
}

impl Arm {
    pub const ALL: [Arm; ARMS] = [Arm::Rens, Arm::Rins, Arm::Mutation, Arm::Crossover];

    pub fn name(self) -> &'static str {
        match self {
            Arm::Rens => "RENS",
            Arm::Rins => "RINS",
            Arm::Mutation => "mutation",
            Arm::Crossover => "crossover",
        }
    }
}

/// Per-arm target fixing rate (SCIP `NH_FIXINGRATE`).
#[derive(Debug, Clone, Copy)]
struct FixingRate {
    target: f64,
    increment: f64,
}

impl FixingRate {
    fn new() -> Self {
        // SCIP starts at the most conservative (largest) rate.
        Self {
            target: MAX_FIXING_RATE,
            increment: FIXING_RATE_START_INC,
        }
    }

    fn update(&mut self, outcome: &SubMipOutcome, found: bool) {
        match outcome {
            // Solved out or infeasible: the neighbourhood was too small.
            SubMipOutcome::FixingsInfeasible
            | SubMipOutcome::Solved(MipStatus::Optimal | MipStatus::Infeasible) => {
                self.target = (self.target - self.increment).max(MIN_FIXING_RATE);
            }
            // Out of nodes with nothing to show: too large.
            SubMipOutcome::Solved(MipStatus::NodeLimit) if !found => {
                self.target = (self.target + self.increment).min(MAX_FIXING_RATE);
            }
            _ => {}
        }
        self.increment = (self.increment * FIXING_RATE_DECAY).max(LRATE_MIN);
    }
}

/// EXP3-IX over the arms (SCIP `bandit_exp3ix.c`; Neu 2015).
#[derive(Debug, Clone)]
struct Exp3Ix {
    weights: [f64; ARMS],
    iter: u64,
}

impl Exp3Ix {
    fn new() -> Self {
        Self {
            weights: [1.0; ARMS],
            iter: 1,
        }
    }

    /// Draws an arm among `available`; returns it with the probability it
    /// was drawn with.
    fn select(&self, available: &[bool; ARMS], rng: &mut XorShift64) -> Option<(usize, f64)> {
        let sum: f64 = (0..ARMS)
            .filter(|&a| available[a])
            .map(|a| self.weights[a])
            .sum();
        let last = (0..ARMS).rev().find(|&a| available[a])?;
        let r = rng.next_f64() * sum;
        let mut acc = 0.0;
        for a in (0..ARMS).filter(|&a| available[a]) {
            acc += self.weights[a];
            if r <= acc || a == last {
                return Some((a, self.weights[a] / sum));
            }
        }
        unreachable!("`last` is available")
    }

    fn update(&mut self, arm: usize, reward: f64, prob: f64) {
        let k = ARMS as f64;
        let gamma = (k.ln() / (4.0 * self.iter as f64 * k)).sqrt();
        let eta = 2.0 * gamma;
        let loss = (1.0 - reward) / (prob + gamma);
        self.weights[arm] *= (-eta * loss).exp();
        // Renormalise so that repeated losses never underflow every weight.
        let max = self.weights.iter().copied().fold(0.0, f64::max);
        if max > 0.0 && max < 1e-100 {
            self.weights.iter_mut().for_each(|w| *w /= max);
        }
        self.iter += 1;
    }
}

/// What the search knows when it offers ALNS a call.
pub struct Call<'x> {
    /// Node LP solution.
    pub lp_x: &'x [f64],
    /// Incumbent value and point.
    pub incumbent: Option<(f64, &'x [f64])>,
    /// Nodes the parent search has processed.
    pub parent_nodes: u64,
    /// Nodes since the incumbent last improved.
    pub nodes_since_incumbent: Option<u64>,
    /// Seconds since the solve started, and until its time limit.
    pub elapsed: f64,
    pub remaining: Option<f64>,
    /// Seconds the root node took.
    pub root_secs: f64,
}

/// One finished ALNS call, to be [`Alns::observe`]d once the search has
/// tried the point.
pub struct Run {
    pub arm: Arm,
    pub report: SubMipReport,
    pub secs: f64,
    prob: f64,
    old_upper: f64,
}

#[derive(Debug, Clone)]
pub struct Alns {
    bandit: Exp3Ix,
    rates: [FixingRate; ARMS],
    rng: XorShift64,
    /// Best solutions seen, ascending by (internal, minimization) objective.
    pool: Vec<(f64, Vec<f64>)>,
    used_nodes: u64,
    calls: u64,
    found: u64,
    target_nodes: u64,
    runs: [u64; ARMS],
    seconds: f64,
    timed_calls: u64,
}

impl Alns {
    pub fn new(seed: u64) -> Self {
        Self {
            bandit: Exp3Ix::new(),
            rates: [FixingRate::new(); ARMS],
            rng: XorShift64::new(0xA125 ^ seed.wrapping_mul(0x9E37_79B9_7F4A_7C15)),
            pool: Vec::new(),
            used_nodes: 0,
            calls: 0,
            found: 0,
            target_nodes: MIN_NODES,
            runs: [0; ARMS],
            seconds: 0.0,
            timed_calls: 0,
        }
    }

    /// Records a new incumbent for crossover.
    pub fn add_solution(&mut self, z: f64, x: &[f64]) {
        if self.pool.iter().any(|(w, y)| *w == z && y == x) {
            return;
        }
        let at = self.pool.partition_point(|(w, _)| *w <= z);
        self.pool.insert(at, (z, x.to_vec()));
        self.pool.truncate(POOL_SIZE);
    }

    /// Node and time limits for a call now, or `None` when the budget says
    /// wait.
    fn limits(&self, call: &Call) -> Option<(u64, f64)> {
        if call
            .nodes_since_incumbent
            .is_some_and(|k| k < WAITING_NODES)
        {
            return None;
        }
        let quot = NODES_QUOT * (self.found as f64 + 1.0) / (self.calls as f64 + 1.0);
        let uninit = self.runs.iter().filter(|&&r| r == 0).count() as f64;
        let k = ARMS as f64;
        let stall = (quot * call.parent_nodes as f64 + NODES_OFFSET
            - self.used_nodes as f64
            - 100.0 * self.calls as f64)
            .min(MAX_NODES as f64)
            * ((k - uninit + 1.0) / (k + 1.0));
        if stall < self.target_nodes as f64 {
            return None;
        }
        let expected = if self.timed_calls == 0 {
            0.5 * call.root_secs
        } else {
            self.seconds / self.timed_calls as f64
        };
        let share = TIME_SHARE * call.elapsed - self.seconds;
        if share < expected {
            return None;
        }
        let time = call.remaining.map_or(share, |r| share.min(r));
        (time > 0.0).then_some((stall as u64, time))
    }

    /// Runs one ALNS call if the budget allows: draws an arm, builds its
    /// neighbourhood and solves it. `None` when nothing was solved and the
    /// budget is untouched.
    pub fn run(&mut self, p: &MipProblem, call: &Call, mip: &MipConfig) -> Option<Run> {
        let (node_limit, time_limit) = self.limits(call)?;
        let start = Instant::now();
        let mut available = [
            true,
            call.incumbent.is_some(),
            call.incumbent.is_some(),
            self.pool.len() >= CROSSOVER_SOLS,
        ];
        loop {
            let (a, prob) = self.bandit.select(&available, &mut self.rng)?;
            let arm = Arm::ALL[a];
            let cfg = SubMipConfig {
                node_limit,
                time_limit: Some(time_limit),
                target_fixing_rate: self.rates[a].target,
                min_fixing_rate: 0.0,
                cutoff: call.incumbent.map(|(z, _)| z),
                seed: self.rng.next_u64(),
                propagate_fixings: false,
                mip: MipConfig {
                    submip: false,
                    alns: false,
                    // Polish stays in the outer search (bnb.rs `run_submip`
                    // / `run_localmip`), never inside a bandit arm's solve.
                    localmip: false,
                    // Sub-MIPs never restart (bnb.rs `run_submip`, HiGHS
                    // `!mipsolver.submip`).
                    restart: false,
                    threads: 1,
                    time_limit: Some(time_limit),
                    node_limit: None,
                    objective_cutoff: None,
                    audit: false,
                    ..*mip
                },
                ..SubMipConfig::default()
            };
            let Some(fix) = self.fixing(p, arm, call, &cfg) else {
                // SCIP: a neighbourhood that cannot run scores 0 and the
                // bandit draws again among the others.
                self.bandit.update(a, 0.0, prob);
                available[a] = false;
                continue;
            };
            let report = if fix.rate() < 0.5 * cfg.target_fixing_rate {
                SubMipReport {
                    x: None,
                    fixing_rate: fix.rate(),
                    outcome: SubMipOutcome::Skipped,
                    sub_cols: 0,
                    sub_rows: 0,
                    nodes: 0,
                    lp_iterations: 0,
                    infeasible_rates: Vec::new(),
                }
            } else {
                submip::solve_restricted(p, &fix, cfg.cutoff, &cfg)
            };
            return Some(Run {
                arm,
                report,
                secs: start.elapsed().as_secs_f64(),
                prob,
                old_upper: call.incumbent.map_or(f64::INFINITY, |(z, _)| z),
            });
        }
    }

    /// Learns from `run`: `improved` says whether the search accepted its
    /// point, `upper` is the incumbent value after it and `lower` the
    /// global dual bound.
    pub fn observe(&mut self, run: &Run, improved: bool, upper: f64, lower: f64) {
        let a = Arm::ALL.iter().position(|&x| x == run.arm).unwrap();
        let r = &run.report;
        self.seconds += run.secs;
        self.timed_calls += 1;
        self.calls += 1;
        self.runs[a] += 1;
        self.used_nodes += r.nodes;
        if improved {
            self.found += 1;
        }
        self.rates[a].update(&r.outcome, improved);
        if self.calls >= ARMS as u64
            && r.outcome == SubMipOutcome::Solved(MipStatus::NodeLimit)
            && !improved
        {
            self.target_nodes =
                ((self.target_nodes as f64 * TARGET_NODE_FACTOR) as u64 + 1).min(MAX_NODES);
        }
        let reward = reward(r, improved, run.old_upper, upper, lower, self.target_nodes);
        self.bandit.update(a, reward, run.prob);
    }

    /// The fixing of `arm` at its target rate, or `None` when the arm
    /// cannot run (SCIP `SCIP_DIDNOTRUN`/`SCIP_DIDNOTFIND`).
    fn fixing(
        &mut self,
        p: &MipProblem,
        arm: Arm,
        call: &Call,
        cfg: &SubMipConfig,
    ) -> Option<Fixing> {
        let feastol = cfg.mip.feastol;
        let target = cfg.target_fixing_rate;
        let limit = cfg
            .cutoff
            .map_or(f64::INFINITY, |z| submip::upper_limit(p, z, feastol));
        let mut fix = Fixing::new(p, cfg, limit);
        if fix.order.is_empty() {
            return None;
        }
        let order = fix.order.clone();
        let near = |a: f64, b: f64| (a - b).abs() <= feastol;
        let (proposals, refsol): (Vec<(usize, f64)>, Option<Vec<f64>>) = match arm {
            Arm::Rens => {
                let x = call.lp_x;
                let fixed = order
                    .iter()
                    .filter(|&&j| near(x[j], x[j].round()))
                    .map(|&j| (j, x[j].round()))
                    .collect();
                (fixed, None)
            }
            Arm::Rins => {
                let (_, inc) = call.incumbent?;
                let x = call.lp_x;
                let fixed = order
                    .iter()
                    .filter(|&&j| near(x[j], inc[j]))
                    .map(|&j| (j, inc[j].round()))
                    .collect();
                (fixed, Some(inc.to_vec()))
            }
            Arm::Mutation => {
                let (_, inc) = call.incumbent?;
                let n = order.len();
                let want = (target * n as f64) as usize + 1;
                if n <= want {
                    return None;
                }
                let fixed = order[..want].iter().map(|&j| (j, inc[j].round())).collect();
                (fixed, None)
            }
            Arm::Crossover => {
                let sols = self.draw_crossover()?;
                let best = &self.pool[sols[0]].1;
                let fixed = order
                    .iter()
                    .filter(|&&j| sols[1..].iter().all(|&s| near(self.pool[s].1[j], best[j])))
                    .map(|&j| (j, best[j].round()))
                    .collect();
                (fixed, Some(best.clone()))
            }
        };
        // Trim an oversized fixing to the target (SCIP unfixes at random
        // above target + UNFIX_TOL; `order` is already random), or apply
        // it whole.
        let n = order.len() as f64;
        let trim = proposals.len() as f64 / n > target + UNFIX_TOL;
        for (j, v) in proposals {
            if trim && fix.rate() >= target {
                break;
            }
            if fix.can_fix(j, v) {
                fix.tighten(j, v, v);
            }
        }
        if fix.rate() < target - FIX_TOL {
            match (arm, refsol) {
                (Arm::Rens, _) => rens_top_up(&mut fix, call.lp_x, target, feastol),
                (_, Some(r)) => {
                    for &j in &order {
                        if fix.rate() >= target {
                            break;
                        }
                        let v = r[j].round();
                        if fix.can_fix(j, v) {
                            fix.tighten(j, v, v);
                        }
                    }
                    if fix.rate() < target - FIX_TOL {
                        return None;
                    }
                }
                (_, None) => {}
            }
        }
        if arm == Arm::Rens {
            // Box the free fractional columns to their two roundings.
            let x = call.lp_x;
            for &j in &order {
                let (lo, hi) = (fix.lower[j], fix.upper[j]);
                if lo < hi {
                    let down = (x[j] + feastol).floor().clamp(lo, hi);
                    let up = (x[j] - feastol).ceil().clamp(lo, hi);
                    if lo < down || up < hi {
                        fix.tighten(j, down.min(up), up.max(down));
                    }
                }
            }
        }
        Some(fix)
    }

    /// Indices into the pool of the crossover parents, best first: SCIP's
    /// draw from ever-better prefixes of the sorted pool.
    fn draw_crossover(&mut self) -> Option<Vec<usize>> {
        let mut need = CROSSOVER_SOLS;
        let mut last = self.pool.len();
        if last < need {
            return None;
        }
        let mut sols = vec![0; need];
        while need > 0 {
            if last == need {
                sols[..need].copy_from_slice(&(0..need).collect::<Vec<_>>());
                break;
            }
            // Uniform in [need − 1, last − 1].
            let span = (last - need + 1) as u64;
            let next = need - 1 + (self.rng.next_u64() % span) as usize;
            sols[need - 1] = next;
            need -= 1;
            last = next;
        }
        Some(sols)
    }
}

/// RENS top-up (SCIP `varFixingsRens`): round the fractional columns
/// closest to an integer first until the target is met.
fn rens_top_up(fix: &mut Fixing, x: &[f64], target: f64, feastol: f64) {
    let mut frac: Vec<(f64, usize)> = fix
        .order
        .iter()
        .filter(|&&j| fix.lower[j] < fix.upper[j])
        .map(|&j| {
            let f = x[j] - x[j].floor();
            (f.min(1.0 - f), j)
        })
        .filter(|&(d, _)| d > feastol)
        .collect();
    frac.sort_by(|a, b| a.0.total_cmp(&b.0));
    for (_, j) in frac {
        if fix.rate() >= target {
            break;
        }
        let v = x[j].round();
        if fix.can_fix(j, v) {
            fix.tighten(j, v, v);
        }
    }
}

/// SCIP `getReward` (see the module header).
fn reward(
    r: &SubMipReport,
    improved: bool,
    old_upper: f64,
    upper: f64,
    lower: f64,
    target_nodes: u64,
) -> f64 {
    let free = 1.0 - r.fixing_rate;
    if improved {
        let gap = if old_upper.is_infinite() || upper <= lower {
            1.0
        } else {
            ((old_upper - upper) / (old_upper - lower)).clamp(0.0, 1.0)
        };
        let effort = r.nodes as f64 / 100.0 * free;
        let raw = (REWARD_CONTROL + (1.0 - REWARD_CONTROL) * gap) / (effort + 1.0);
        REWARD_BASELINE + (1.0 - REWARD_BASELINE) * raw
    } else {
        let used = r.nodes as f64 * free;
        REWARD_BASELINE * (1.0 - used / target_nodes as f64).max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{LpEngine, LpStatus, RusolveEngine};
    use rusolve_core::CscMatrix;

    #[test]
    fn fixing_rate_moves_and_decays() {
        let mut f = FixingRate::new();
        assert_eq!(f.target, MAX_FIXING_RATE);
        f.update(&SubMipOutcome::Solved(MipStatus::Optimal), false);
        assert!((f.target - 0.7).abs() < 1e-12);
        assert!((f.increment - 0.15).abs() < 1e-12);
        f.update(&SubMipOutcome::Solved(MipStatus::NodeLimit), false);
        assert!((f.target - 0.85).abs() < 1e-12);
        // A node limit that found something keeps the rate.
        f.update(&SubMipOutcome::Solved(MipStatus::NodeLimit), true);
        assert!((f.target - 0.85).abs() < 1e-12);
        for _ in 0..50 {
            f.update(&SubMipOutcome::FixingsInfeasible, false);
        }
        assert_eq!(f.target, MIN_FIXING_RATE);
        assert_eq!(f.increment, LRATE_MIN);
    }

    #[test]
    fn bandit_prefers_the_rewarded_arm_and_skips_unavailable_ones() {
        let mut b = Exp3Ix::new();
        let mut rng = XorShift64::new(7);
        let all = [true; ARMS];
        for _ in 0..200 {
            let (a, p) = b.select(&all, &mut rng).unwrap();
            b.update(a, if a == 2 { 1.0 } else { 0.0 }, p);
        }
        let mut hits = [0; ARMS];
        for _ in 0..1000 {
            hits[b.select(&all, &mut rng).unwrap().0] += 1;
        }
        assert!(hits[2] > 700, "{hits:?}");
        let only = [true, false, false, false];
        for _ in 0..100 {
            assert_eq!(b.select(&only, &mut rng).unwrap(), (0, 1.0));
        }
        assert!(b.select(&[false; ARMS], &mut rng).is_none());
    }

    #[test]
    fn reward_matches_scip_formula() {
        let mut r = SubMipReport {
            x: None,
            fixing_rate: 0.5,
            outcome: SubMipOutcome::Solved(MipStatus::NodeLimit),
            sub_cols: 0,
            sub_rows: 0,
            nodes: 40,
            lp_iterations: 0,
            infeasible_rates: Vec::new(),
        };
        // Failure: 0.5 · (1 − 40 · 0.5 / 50).
        assert!((reward(&r, false, 10.0, 10.0, 0.0, 50) - 0.3).abs() < 1e-12);
        // Success closing half the gap: effort 0.2,
        // 0.5 + 0.5 · (0.8 + 0.2 · 0.5) / 1.2.
        r.nodes = 40;
        let want = 0.5 + 0.5 * 0.9 / 1.2;
        assert!((reward(&r, true, 10.0, 5.0, 0.0, 50) - want).abs() < 1e-12);
        // First solution: full gap credit.
        r.nodes = 0;
        assert!((reward(&r, true, f64::INFINITY, 5.0, 0.0, 50) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pool_keeps_the_best_distinct_solutions_sorted() {
        let mut a = Alns::new(1);
        for k in (0..15).rev() {
            a.add_solution(k as f64, &[k as f64]);
        }
        a.add_solution(0.0, &[0.0]);
        assert_eq!(a.pool.len(), POOL_SIZE);
        assert!(a.pool.windows(2).all(|w| w[0].0 < w[1].0));
        assert_eq!(a.pool[0].0, 0.0);
    }

    #[test]
    fn crossover_draw_is_two_distinct_pool_indices_best_first() {
        let mut a = Alns::new(3);
        for k in 0..6 {
            a.add_solution(k as f64, &[k as f64]);
        }
        for _ in 0..200 {
            let s = a.draw_crossover().unwrap();
            assert_eq!(s.len(), 2);
            assert!(s[0] < s[1] && s[1] < 6, "{s:?}");
        }
        let mut b = Alns::new(3);
        b.add_solution(1.0, &[1.0]);
        assert!(b.draw_crossover().is_none());
        b.add_solution(2.0, &[2.0]);
        assert_eq!(b.draw_crossover().unwrap(), vec![0, 1]);
    }

    #[test]
    fn budget_waits_for_root_time_and_new_incumbents() {
        let a = Alns::new(1);
        let x = [0.5];
        let mut call = Call {
            lp_x: &x,
            incumbent: None,
            parent_nodes: 100,
            nodes_since_incumbent: None,
            elapsed: 1.0,
            remaining: None,
            root_secs: 1.0,
        };
        // 10 % of 1 s < half the root time.
        assert!(a.limits(&call).is_none());
        call.elapsed = 6.0;
        let (nodes, time) = a.limits(&call).unwrap();
        // (0.1 · 100 + 500) · (0 + 1)/(4 + 1) = 102.
        assert_eq!(nodes, 102);
        assert!((time - 0.6).abs() < 1e-12);
        call.nodes_since_incumbent = Some(10);
        assert!(a.limits(&call).is_none());
        call.remaining = Some(0.25);
        call.nodes_since_incumbent = Some(30);
        assert_eq!(a.limits(&call).unwrap().1, 0.25);
    }

    /// Two-row 0/1 knapsack, `max Σ v_j x_j`, with weights drawn from `rng`.
    fn knapsack(rng: &mut XorShift64, n: usize) -> MipProblem {
        let mut trip = Vec::new();
        for i in 0..2 {
            for j in 0..n {
                trip.push((i, j, 1.0 + (rng.next_u64() % 9) as f64));
            }
        }
        MipProblem {
            num_col: n,
            num_row: 2,
            matrix: CscMatrix::from_triplets(2, n, &trip),
            row_lower: vec![f64::NEG_INFINITY; 2],
            row_upper: vec![2.0 * n as f64; 2],
            cost: (0..n)
                .map(|_| -(1.0 + (rng.next_u64() % 20) as f64))
                .collect(),
            maximize: false,
            col_lower: vec![0.0; n],
            col_upper: vec![1.0; n],
            integer: vec![true; n],
        }
    }

    #[test]
    fn calls_return_only_feasible_improving_points() {
        let mut rng = XorShift64::new(20261023);
        let mip = MipConfig::default();
        let mut arms_run = [false; ARMS];
        let mut improved_any = 0;
        for case in 0..10 {
            let p = knapsack(&mut rng, 16 + case);
            let mut e = RusolveEngine::new(&p, mip.lp);
            let lp = e.solve(&p.col_lower, &p.col_upper).unwrap();
            assert_eq!(lp.status, LpStatus::Optimal);
            let mut a = Alns::new(case as u64);
            let mut inc = vec![0.0; p.num_col];
            let mut z = p.objective(&inc);
            let start = z;
            a.add_solution(z, &inc);
            // A second pool point: the first column alone.
            let mut y = inc.clone();
            y[0] = 1.0;
            if p.violation(&y).row <= 1e-9 {
                a.add_solution(p.objective(&y), &y);
            }
            for k in 0..30 {
                let call = Call {
                    lp_x: &lp.x,
                    incumbent: Some((z, &inc)),
                    parent_nodes: 10_000 * (k + 1),
                    nodes_since_incumbent: Some(1000),
                    elapsed: 1e6,
                    remaining: None,
                    root_secs: 0.0,
                };
                let Some(run) = a.run(&p, &call, &mip) else {
                    break;
                };
                arms_run[Arm::ALL.iter().position(|&x| x == run.arm).unwrap()] = true;
                let improved = match &run.report.x {
                    Some(x) => {
                        let v = p.violation(x);
                        assert!(v.row <= 1e-6 && v.bound <= 1e-6 && v.integrality <= 1e-6);
                        let w = p.objective(x);
                        assert!(w < z, "case {case}: {w} does not beat {z}");
                        z = w;
                        inc = x.clone();
                        a.add_solution(z, &inc);
                        true
                    }
                    None => false,
                };
                a.observe(&run, improved, z, lp.objective);
            }
            if z < start {
                improved_any += 1;
            }
        }
        assert!(improved_any >= 5, "{improved_any}");
        assert!(arms_run.iter().filter(|&&r| r).count() >= 3, "{arms_run:?}");
    }
}
