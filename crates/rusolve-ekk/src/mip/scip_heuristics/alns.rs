//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! `mip_alns`: bandit-scheduled adaptive large neighbourhood search on the port stack.
//! Arms: RENS, RINS, mutation,
//! crossover. The code is a translation of the native stack's `crates/rusolve-mip/src/alns.rs`
//! (ours), which derives from SCIP `heur_alns.c` and `bandit_exp3ix.c`:
//!
//! * G. Hendel, "Adaptive large neighborhood search for mixed integer programming", Math. Prog.
//!   Comp. 14:185-221 (2022): the bandit over neighbourhoods, the reward, the fixing rate.
//! * G. Neu, "Explore no more: improved high-probability regret bounds for non-stochastic
//!   bandits", NeurIPS 2015: EXP3-IX (`bandit_exp3ix.c`).
//! * SCIP 11.0.0 (32dd425) heur_alns.c: `updateFixingRate` (start 0.9, min 0.3, step 0.2 decaying
//!   by 0.75 down to 0.01), `neighborhoodFixVariables` (target within +-0.1 of the rate),
//!   `determineLimits` (`0.1 (found+1)/(calls+1) nnodes + 500 - used - 100 calls`, at most 5000,
//!   scaled by `(K - uninitialised + 1)/(K + 1)`), `updateTargetNodeLimit` (50, then x1.05 + 1),
//!   `getReward`, `varFixingsRens`, `varFixingsMutation`, `varFixingsCrossover`;
//!   `DEFAULT_WAITINGNODES` 25.
//!
//! Everything runs on the main thread at the SCIP tree-heuristic hook sites, where no parallel lock
//! is active. The restriction is solved by `s6_solve_sub_mip_src`, the SCIP-heuristic sub-MIP copy; its
//! solution reaches `try_solution` with `kSolutionSourceAlns`, so the tree and the pool see it.

use super::*;
use crate::mip::mip_solver::HighsMipSolver;
use crate::mip::mip_solver_data::kSolutionSourceAlns;
use crate::mip::primal_heuristics::rens_get_fix_val;
use crate::work::WorkSnapshot;

pub(crate) const ALNS_ARMS: usize = 4;
pub(crate) const ALNS_RENS: usize = 0;
pub(crate) const ALNS_RINS: usize = 1;
pub(crate) const ALNS_MUTATION: usize = 2;
pub(crate) const ALNS_CROSSOVER: usize = 3;
pub(crate) const ALNS_ARM_NAMES: [&str; ALNS_ARMS] = ["RENS", "RINS", "mutation", "crossover"];

/// SCIP heur_alns.c: `DEFAULT_MINFIXINGRATE_SUBMIP` side of the adaptive range (0.3) and its start.
pub(crate) const ALNS_MIN_FIXING_RATE: f64 = 0.3;
pub(crate) const ALNS_MAX_FIXING_RATE: f64 = 0.9;
pub(crate) const ALNS_FIXING_RATE_START_INC: f64 = 0.2;
pub(crate) const ALNS_FIXING_RATE_DECAY: f64 = 0.75;
pub(crate) const ALNS_LRATE_MIN: f64 = 0.01;
/// `DEFAULT_FIXTOL` / `DEFAULT_UNFIXTOL`.
pub(crate) const ALNS_FIX_TOL: f64 = 0.1;
pub(crate) const ALNS_UNFIX_TOL: f64 = 0.1;
pub(crate) const ALNS_NODES_QUOT: f64 = 0.1;
pub(crate) const ALNS_NODES_OFFSET: f64 = 500.0;
pub(crate) const ALNS_MIN_NODES: i64 = 50;
pub(crate) const ALNS_MAX_NODES: i64 = 5000;
pub(crate) const ALNS_TARGET_NODE_FACTOR: f64 = 1.05;
/// `DEFAULT_WAITINGNODES`.
pub(crate) const ALNS_WAITING_NODES: i64 = 25;
pub(crate) const ALNS_REWARD_CONTROL: f64 = 0.8;
pub(crate) const ALNS_REWARD_BASELINE: f64 = 0.5;
/// Crossover parents per call (the native arm; SCIP `DEFAULT_NSOLS` 2).
pub(crate) const ALNS_CROSSOVER_SOLS: usize = 2;
/// Work share of the main search that ALNS may spend: SCIP's node quotient 0.1 (Hendel 2022),
/// in work units (`src/work.rs`) instead of nodes.
pub(crate) const ALNS_WORK_SHARE: f64 = 0.1;
/// Work ALNS may spend before the share applies, as a fraction of the root's work (frozen at the
/// first hook). Ours: the native stack's "half a root before the first call" rule in work units.
pub(crate) const ALNS_WORK_OFFSET_ROOT_FRACTION: f64 = 0.5;
/// Seed of the bandit RNG (ours); shifted by `random_seed` like the crossover seed.
pub(crate) const ALNS_RANDSEED: u32 = 0xA125;

/// How a restriction's sub-MIP ended, as far as the fixing rate cares.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum AlnsOutcome {
    /// Proven infeasible (or cut off): the neighbourhood was too small.
    Infeasible,
    /// Solved out: too small.
    Optimal,
    /// Out of nodes: too large unless something was found.
    NodeLimit,
    /// Terminated from outside, or never ran.
    Other,
}

/// Per-arm target fixing rate (SCIP `NH_FIXINGRATE`, `updateFixingRate`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct AlnsFixingRate {
    pub target: f64,
    pub increment: f64,
}

impl AlnsFixingRate {
    pub(crate) fn new() -> Self {
        // SCIP starts at the most conservative (largest) rate.
        Self {
            target: ALNS_MAX_FIXING_RATE,
            increment: ALNS_FIXING_RATE_START_INC,
        }
    }

    pub(crate) fn update(&mut self, outcome: AlnsOutcome, found: bool) {
        match outcome {
            AlnsOutcome::Infeasible | AlnsOutcome::Optimal => {
                self.target = (self.target - self.increment).max(ALNS_MIN_FIXING_RATE);
            }
            AlnsOutcome::NodeLimit if !found => {
                self.target = (self.target + self.increment).min(ALNS_MAX_FIXING_RATE);
            }
            _ => {}
        }
        self.increment = (self.increment * ALNS_FIXING_RATE_DECAY).max(ALNS_LRATE_MIN);
    }
}

/// EXP3-IX over the arms (SCIP `bandit_exp3ix.c`; Neu 2015).
#[derive(Clone, Debug)]
pub(crate) struct AlnsExp3Ix {
    pub weights: [f64; ALNS_ARMS],
    pub iter: u64,
}

impl AlnsExp3Ix {
    pub(crate) fn new() -> Self {
        Self {
            weights: [1.0; ALNS_ARMS],
            iter: 1,
        }
    }

    /// Draws an arm among `available` with the uniform variate `u` in `[0, 1]`; returns it with
    /// the probability it was drawn with, renormalised over the available arms.
    pub(crate) fn select(&self, available: &[bool; ALNS_ARMS], u: f64) -> Option<(usize, f64)> {
        let sum: f64 = (0..ALNS_ARMS)
            .filter(|&a| available[a])
            .map(|a| self.weights[a])
            .sum();
        let last = (0..ALNS_ARMS).rev().find(|&a| available[a])?;
        let r = u * sum;
        let mut acc = 0.0;
        for a in (0..ALNS_ARMS).filter(|&a| available[a]) {
            acc += self.weights[a];
            if r <= acc || a == last {
                return Some((a, self.weights[a] / sum));
            }
        }
        unreachable!("`last` is available")
    }

    /// `gamma_t = sqrt(ln K / (4 t K))`, `eta = 2 gamma`, loss `(1 - r)/(p + gamma)`,
    /// weight `*= exp(-eta loss)`.
    pub(crate) fn update(&mut self, arm: usize, reward: f64, prob: f64) {
        let k = ALNS_ARMS as f64;
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

/// SCIP `getReward` (heur_alns.c): with a new incumbent
/// `0.5 + 0.5 (0.8 + 0.2 closed gap)/(effort + 1)`, effort `nodes/100 (1 - fixing rate)`; without
/// one `0.5 max(0, 1 - nodes (1 - fixing rate)/target nodes)`.
pub(crate) fn alns_reward(
    nodes: i64,
    fixing_rate: f64,
    improved: bool,
    old_upper: f64,
    upper: f64,
    lower: f64,
    target_nodes: i64,
) -> f64 {
    let free = 1.0 - fixing_rate;
    if improved {
        let gap = if old_upper.is_infinite() || upper <= lower || lower.is_infinite() {
            1.0
        } else {
            ((old_upper - upper) / (old_upper - lower)).clamp(0.0, 1.0)
        };
        let effort = nodes as f64 / 100.0 * free;
        let raw = (ALNS_REWARD_CONTROL + (1.0 - ALNS_REWARD_CONTROL) * gap) / (effort + 1.0);
        ALNS_REWARD_BASELINE + (1.0 - ALNS_REWARD_BASELINE) * raw
    } else {
        let used = nodes as f64 * free;
        ALNS_REWARD_BASELINE * (1.0 - used / target_nodes as f64).max(0.0)
    }
}

/// SCIP `determineLimits`, with the work clamp of the design: `work_node_cap` is the number of
/// sub-MIP nodes the remaining work budget pays for. Returns `(max nodes, stall nodes)` or `None`
/// when the budget is below the target node count.
pub(crate) fn alns_determine_limits(
    found: u64,
    calls: u64,
    uninitialised_arms: u64,
    used_nodes: i64,
    parent_nodes: i64,
    target_nodes: i64,
    work_node_cap: f64,
) -> Option<(i64, i64)> {
    let quot = ALNS_NODES_QUOT * (found as f64 + 1.0) / (calls as f64 + 1.0);
    let k = ALNS_ARMS as f64;
    let mut nodes =
        (quot * parent_nodes as f64 + ALNS_NODES_OFFSET - used_nodes as f64 - 100.0 * calls as f64)
            .min(ALNS_MAX_NODES as f64)
            * ((k - uninitialised_arms as f64 + 1.0) / (k + 1.0));
    nodes = nodes.min(work_node_cap);
    if nodes < target_nodes as f64 {
        return None;
    }
    Some((nodes as i64, target_nodes))
}

/// The bandit's state for one solve (`HighsMipSolverData::alns`).
#[derive(Clone, Debug)]
pub(crate) struct AlnsState {
    pub bandit: AlnsExp3Ix,
    pub rates: [AlnsFixingRate; ALNS_ARMS],
    pub rng: S4RandNumGen,
    pub calls: u64,
    pub found: u64,
    /// Sub-MIP nodes spent so far (SCIP `usednodes`).
    pub used_nodes: i64,
    pub target_nodes: i64,
    pub runs: [u64; ALNS_ARMS],
    pub found_by_arm: [u64; ALNS_ARMS],
    /// Work spent in sub-MIP calls (the ALNS work ledger).
    pub alns_work: u64,
    pub initialised: bool,
    /// Work allowed before the share applies, frozen at the first hook.
    pub offset_work: f64,
}

impl AlnsState {
    pub(crate) fn new(random_seed: HighsInt) -> Self {
        Self {
            bandit: AlnsExp3Ix::new(),
            rates: [AlnsFixingRate::new(); ALNS_ARMS],
            rng: S4RandNumGen::new(s4_initialize_random_seed_from(ALNS_RANDSEED, random_seed)),
            calls: 0,
            found: 0,
            used_nodes: 0,
            target_nodes: ALNS_MIN_NODES,
            runs: [0; ALNS_ARMS],
            found_by_arm: [0; ALNS_ARMS],
            alns_work: 0,
            initialised: false,
            offset_work: 0.0,
        }
    }

    fn uniform(&mut self) -> f64 {
        self.rng.get_real(0.0, 1.0)
    }

    /// Work per sub-MIP node: the mean of earlier calls. `None` before the first call, which is
    /// then bounded by the node formula alone (the root's work is no estimate of a node's).
    pub(crate) fn work_per_node(&self) -> Option<f64> {
        (self.used_nodes > 0).then(|| (self.alns_work as f64 / self.used_nodes as f64).max(1.0))
    }

    /// Budget gate and node limits at the hook; `main_work` is the solve's work minus the ledger.
    pub(crate) fn limits(&self, main_work: u64, parent_nodes: i64) -> Option<(i64, i64)> {
        let allowance = ALNS_WORK_SHARE * main_work as f64 + self.offset_work;
        let spent = self.alns_work as f64;
        if spent > allowance {
            return None;
        }
        let work_node_cap = self
            .work_per_node()
            .map_or(f64::INFINITY, |w| (allowance - spent) / w);
        let uninit = self.runs.iter().filter(|&&r| r == 0).count() as u64;
        alns_determine_limits(
            self.found,
            self.calls,
            uninit,
            self.used_nodes,
            parent_nodes,
            self.target_nodes,
            work_node_cap,
        )
    }

    /// Learns from one finished call (SCIP `heurExecAlns` after the sub-SCIP solve).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn observe(
        &mut self,
        arm: usize,
        prob: f64,
        outcome: AlnsOutcome,
        nodes: i64,
        work: u64,
        fixing_rate: f64,
        improved: bool,
        old_upper: f64,
        upper: f64,
        lower: f64,
    ) {
        self.calls += 1;
        self.runs[arm] += 1;
        self.used_nodes += nodes;
        self.alns_work += work;
        if improved {
            self.found += 1;
            self.found_by_arm[arm] += 1;
        }
        self.rates[arm].update(outcome, improved);
        if self.calls >= ALNS_ARMS as u64 && outcome == AlnsOutcome::NodeLimit && !improved {
            self.target_nodes = (((self.target_nodes as f64 * ALNS_TARGET_NODE_FACTOR) as i64) + 1)
                .min(ALNS_MAX_NODES);
        }
        let reward = alns_reward(
            nodes,
            fixing_rate,
            improved,
            old_upper,
            upper,
            lower,
            self.target_nodes,
        );
        self.bandit.update(arm, reward, prob);
    }
}

/// A restriction under construction: the global bounds plus the fixings of one neighbourhood.
struct Restriction {
    lower: Vec<f64>,
    upper: Vec<f64>,
    /// Integer columns currently fixed, and all integer columns.
    nfixed: usize,
    nint: usize,
    /// Free integer columns in seeded random order (SCIP `alnsFixMoreVariables` orders them by
    /// graph distance, reduced cost and pseudocost; phase B).
    order: Vec<usize>,
}

impl Restriction {
    fn new(ms: &HighsMipSolver, rng: &mut S4RandNumGen) -> Option<Self> {
        let d = &ms.mipdata_.domains[0];
        let lower = d.col_lower_.clone();
        let upper = d.col_upper_.clone();
        let nint = ms.mipdata_.integral_cols.len();
        if nint == 0 {
            return None;
        }
        let mut nfixed = 0;
        let mut order: Vec<usize> = Vec::new();
        for &icol in &ms.mipdata_.integral_cols {
            let j = icol as usize;
            if lower[j] == upper[j] {
                nfixed += 1;
            } else {
                order.push(j);
            }
        }
        if order.is_empty() {
            return None;
        }
        // Fisher-Yates with SCIP's randomGetInt.
        for i in (1..order.len()).rev() {
            let k = s6_get_int(rng, 0, i as HighsInt) as usize;
            order.swap(i, k);
        }
        Some(Self {
            lower,
            upper,
            nfixed,
            nint,
            order,
        })
    }

    fn rate(&self) -> f64 {
        self.nfixed as f64 / self.nint as f64
    }

    fn can_fix(&self, j: usize, v: f64) -> bool {
        self.lower[j] < self.upper[j] && self.lower[j] <= v && v <= self.upper[j]
    }

    fn fix(&mut self, j: usize, v: f64) {
        self.lower[j] = v;
        self.upper[j] = v;
        self.nfixed += 1;
    }

    fn try_fix(&mut self, j: usize, v: f64) -> bool {
        let ok = self.can_fix(j, v);
        if ok {
            self.fix(j, v);
        }
        ok
    }
}

/// What the arm proposes: a value per column (in `order`'s priority) or a ready fixing.
struct Proposal {
    /// `(column, value)` in random priority order.
    fixings: Vec<(usize, f64)>,
    /// Reference point to top up a short fixing from.
    refsol: Option<Vec<f64>>,
}

fn alns_near(a: f64, b: f64, feastol: f64) -> bool {
    (a - b).abs() <= feastol
}

/// The restriction of `arm` at its target fixing rate, or `None` when the arm cannot run
/// (SCIP `SCIP_DIDNOTRUN` / `SCIP_DIDNOTFIND`).
fn alns_neighbourhood(
    ms: &mut HighsMipSolver,
    st: &mut AlnsState,
    arm: usize,
) -> Option<Restriction> {
    let feastol = ms.mipdata_.feastol;
    let target = st.rates[arm].target;
    let mut fix = Restriction::new(ms, &mut st.rng)?;
    let lpx: Vec<f64> = ms.mipdata_.rootlpsol.clone();
    let inc: Vec<f64> = ms.mipdata_.incumbent.clone();
    let n_cols = mip_model!(ms).num_col_ as usize;
    let have_lp = lpx.len() == n_cols;
    let have_inc = inc.len() == n_cols;

    let prop: Proposal = match arm {
        // RENS: fix the columns whose relaxation value is integral.
        ALNS_RENS => {
            if !have_lp {
                return None;
            }
            let fixings = fix
                .order
                .iter()
                .filter(|&&j| alns_near(lpx[j], lpx[j].round(), feastol))
                .map(|&j| (j, lpx[j].round()))
                .collect();
            Proposal {
                fixings,
                refsol: None,
            }
        }
        // RINS: fix where the incumbent and the relaxation agree.
        ALNS_RINS => {
            if !have_lp || !have_inc {
                return None;
            }
            let fixings = fix
                .order
                .iter()
                .filter(|&&j| alns_near(lpx[j], inc[j], feastol))
                .map(|&j| (j, inc[j].round()))
                .collect();
            Proposal {
                fixings,
                refsol: Some(inc),
            }
        }
        // Mutation: `floor(rate n) + 1` random integer columns to the incumbent.
        ALNS_MUTATION => {
            if !have_inc {
                return None;
            }
            let want = (target * fix.nint as f64) as usize + 1;
            let extra = want.saturating_sub(fix.nfixed);
            if extra == 0 || fix.order.len() <= extra {
                return None;
            }
            let fixings = fix.order[..extra]
                .iter()
                .map(|&j| (j, inc[j].round()))
                .collect();
            Proposal {
                fixings,
                refsol: None,
            }
        }
        // Crossover: fix where two pool solutions, drawn by SCIP's randomised selection of the
        // crossover, agree.
        ALNS_CROSSOVER => {
            if !have_inc {
                return None;
            }
            let integrality: Vec<HighsVarType> = mip_model!(ms).integrality_.clone();
            let col_lower = fix.lower.clone();
            let col_upper = fix.upper.clone();
            let s6 = ms.mipdata_.s6.as_mut()?;
            if s6.pool.sols.len() < ALNS_CROSSOVER_SOLS || s6.xover_rng.is_none() {
                return None;
            }
            let mut selection: Vec<i32> = vec![0; ALNS_CROSSOVER_SOLS];
            if !crate::mip::scip_heuristics::crossover::s6_select_sols_randomized(
                s6,
                &mut selection,
                ALNS_CROSSOVER_SOLS,
            ) {
                return None;
            }
            let mut agreed: Vec<(usize, f64)> = Vec::new();
            let _ = crate::mip::scip_heuristics::crossover::s6_fix_variables(
                s6,
                &integrality,
                &col_lower,
                &col_upper,
                &selection,
                ALNS_CROSSOVER_SOLS,
                &mut agreed,
            );
            let best = s6.pool.sols[selection[0] as usize].x.clone();
            let mut value: Vec<Option<f64>> = vec![None; n_cols];
            for (j, v) in agreed {
                value[j] = Some(v.round());
            }
            let fixings = fix
                .order
                .iter()
                .filter_map(|&j| value[j].map(|v| (j, v)))
                .collect();
            Proposal {
                fixings,
                refsol: Some(best),
            }
        }
        _ => return None,
    };

    // Trim an oversized fixing to the target (SCIP unfixes at random above target + UNFIX_TOL;
    // `order` is random already), or apply it whole. Mutation is exact by construction.
    let order: Vec<usize> = fix.order.clone();
    let proposed = prop.fixings.len();
    let trim = arm != ALNS_MUTATION
        && (fix.nfixed + proposed) as f64 / fix.nint as f64 > target + ALNS_UNFIX_TOL;
    for (j, v) in prop.fixings {
        if trim && fix.rate() >= target {
            break;
        }
        fix.try_fix(j, v);
    }
    if fix.rate() < target - ALNS_FIX_TOL {
        if arm == ALNS_RENS {
            // Round the fractional columns closest to their HiGHS RENS fixing value first.
            let dom = &ms.mipdata_.domains[0];
            let mut frac: Vec<(f64, usize, f64)> = order
                .iter()
                .filter(|&&j| {
                    fix.lower[j] < fix.upper[j] && !alns_near(lpx[j], lpx[j].round(), feastol)
                })
                .map(|&j| {
                    let v = rens_get_fix_val(ms, dom, j as HighsInt, lpx[j]);
                    ((v - lpx[j]).abs(), j, v)
                })
                .collect();
            frac.sort_by(|a, b| a.0.total_cmp(&b.0));
            for (_, j, v) in frac {
                if fix.rate() >= target {
                    break;
                }
                fix.try_fix(j, v.round());
            }
        } else if let Some(r) = &prop.refsol {
            for &j in &order {
                if fix.rate() >= target {
                    break;
                }
                fix.try_fix(j, r[j].round());
            }
            if fix.rate() < target - ALNS_FIX_TOL {
                return None;
            }
        }
    }
    if arm == ALNS_RENS {
        // Box the free fractional columns to their two roundings.
        for &j in &order {
            let (lo, hi) = (fix.lower[j], fix.upper[j]);
            if lo < hi {
                let down = (lpx[j] + feastol).floor().clamp(lo, hi);
                let up = (lpx[j] - feastol).ceil().clamp(lo, hi);
                if lo < down || up < hi {
                    fix.lower[j] = down.min(up);
                    fix.upper[j] = up.max(down);
                }
            }
        }
    }
    // No call below half the target fixing rate (SCIP `heurExecAlns`); measured before the
    // sub-MIP's own presolve.
    if fix.rate() < 0.5 * target {
        return None;
    }
    Some(fix)
}

fn alns_classify(ok: bool, status: HighsModelStatus) -> AlnsOutcome {
    match status {
        HighsModelStatus::kInfeasible => AlnsOutcome::Infeasible,
        HighsModelStatus::kOptimal if ok => AlnsOutcome::Optimal,
        HighsModelStatus::kSolutionLimit | HighsModelStatus::kIterationLimit => {
            AlnsOutcome::NodeLimit
        }
        _ => AlnsOutcome::Other,
    }
}

/// One scheduling step at a hook site: runs at most one sub-MIP if the work ledger and the
/// node budget allow. `after_root` is the single call right after the root, which is exempt from
/// the 25-node wait for a new incumbent.
pub(crate) fn alns_step(ms: &mut HighsMipSolver, after_root: bool) {
    if ms.submip || !ms.options_mip_.mip_alns || ms.terminate() {
        return;
    }
    if ms.mipdata_.nodequeue.empty() || ms.mipdata_.alns.is_none() {
        return;
    }
    let nnodes = ms.mipdata_.num_nodes;
    let have_inc = !ms.mipdata_.incumbent.is_empty();
    let n_cols = mip_model!(ms).num_col_ as usize;
    let have_lp = ms.mipdata_.rootlpsol.len() == n_cols;
    let pool_len = ms.mipdata_.s6.as_ref().map_or(0, |s| s.pool.sols.len());
    if have_inc && !after_root {
        // SCIP `DEFAULT_WAITINGNODES`: wait after a new incumbent.
        let best_node = ms
            .mipdata_
            .s6
            .as_ref()
            .and_then(|s| s.pool.best().map(|b| b.nodenum))
            .unwrap_or(0);
        if nnodes - best_node < ALNS_WAITING_NODES {
            return;
        }
    }
    let mut st = ms.mipdata_.alns.take().expect("checked above");
    alns_step_inner(ms, &mut st, nnodes, have_inc, have_lp, pool_len);
    ms.mipdata_.alns = Some(st);
}

fn alns_step_inner(
    ms: &mut HighsMipSolver,
    st: &mut AlnsState,
    nnodes: i64,
    have_inc: bool,
    have_lp: bool,
    pool_len: usize,
) {
    let total_work = WorkSnapshot::now().since(&ms.stats_.work_start).work();
    let main_work = total_work.saturating_sub(st.alns_work);
    if !st.initialised {
        st.initialised = true;
        st.offset_work = ALNS_WORK_OFFSET_ROOT_FRACTION * main_work as f64;
    }
    let Some((maxnodes, stallnodes)) = st.limits(main_work, nnodes) else {
        return;
    };
    let mut available = [
        have_lp,
        have_lp && have_inc,
        have_inc,
        have_inc && pool_len >= ALNS_CROSSOVER_SOLS,
    ];
    loop {
        let u = st.uniform();
        let Some((arm, prob)) = st.bandit.select(&available, u) else {
            return;
        };
        let Some(fix) = alns_neighbourhood(ms, st, arm) else {
            // SCIP: a neighbourhood that cannot run scores 0 and the bandit draws again.
            st.bandit.update(arm, 0.0, prob);
            available[arm] = false;
            continue;
        };
        let fixing_rate = fix.rate();
        let lp: HighsLp = mip_model!(ms).clone();
        let basis: HighsBasis = ms.mipdata_.firstrootbasis.clone();
        let old_upper = ms.mipdata_.upper_bound;
        let w0 = WorkSnapshot::now();
        let (ok, subnodes, status) = detach!(ms.mipdata_.heuristics, |h| h.s6_solve_sub_mip_src(
            ms,
            0,
            lp,
            basis,
            fixing_rate,
            fix.lower,
            fix.upper,
            kHighsIInf,
            maxnodes as HighsInt,
            stallnodes as HighsInt,
            f64::INFINITY,
            kSolutionSourceAlns,
        ));
        let work = WorkSnapshot::now().since(&w0).work();
        let upper = ms.mipdata_.upper_bound;
        let improved = upper < old_upper;
        st.observe(
            arm,
            prob,
            alns_classify(ok, status),
            subnodes,
            work,
            fixing_rate,
            improved,
            old_upper,
            upper,
            ms.mipdata_.lower_bound,
        );
        let s = &mut ms.stats_;
        s.alns_calls += 1;
        s.alns_calls_by_arm[arm] += 1;
        s.alns_nodes += subnodes;
        s.alns_work = st.alns_work;
        if improved {
            s.alns_found += 1;
            s.alns_found_by_arm[arm] += 1;
        }
        crate::io::highs_io::highs_log_dev(
            &ms.options_mip_.log_options,
            HighsLogType::kInfo,
            &format!(
                "ALNS {} rate {:.2} nodes {} work {} improved {}\n",
                ALNS_ARM_NAMES[arm], fixing_rate, subnodes, work, improved
            ),
        );
        return;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bandit_update_matches_hand_values() {
        // K = 4, t = 1: gamma = sqrt(ln 4 / 16) = 0.2943524..., eta = 2 gamma.
        let gamma = (4f64.ln() / 16.0).sqrt();
        assert!((gamma - 0.294_352_3).abs() < 1e-6);
        let mut b = AlnsExp3Ix::new();
        // Arm 1 drawn with p = 0.25 and reward 0.5: loss = 0.5 / (0.25 + gamma),
        // weight = exp(-2 gamma loss) = 0.582 (hand: 0.5 / 0.544352 = 0.918531;
        // exp(-0.588705 * 0.918531) = 0.582378).
        b.update(1, 0.5, 0.25);
        assert!((b.weights[1] - 0.582_318).abs() < 1e-5, "{:?}", b.weights);
        assert_eq!(b.weights[0], 1.0);
        assert_eq!(b.iter, 2);
        // Reward 1 leaves the weight (zero loss).
        b.update(0, 1.0, 0.3);
        assert_eq!(b.weights[0], 1.0);
        // Third step (t = 3): gamma = sqrt(ln 4 / 48); reward 0 at p = 0.5 has loss 1/(0.5 + gamma).
        let g3 = (4f64.ln() / 48.0).sqrt();
        b.update(2, 0.0, 0.5);
        let want = (-2.0 * g3 / (0.5 + g3)).exp();
        assert!((b.weights[2] - want).abs() < 1e-12);
    }

    #[test]
    fn bandit_selection_follows_weights_and_availability() {
        let mut b = AlnsExp3Ix::new();
        b.weights = [1.0, 2.0, 3.0, 4.0];
        let all = [true; ALNS_ARMS];
        // Cumulative 1, 3, 6, 10 over sum 10.
        assert_eq!(b.select(&all, 0.05).unwrap(), (0, 0.1));
        assert_eq!(b.select(&all, 0.25).unwrap(), (1, 0.2));
        assert_eq!(b.select(&all, 0.5).unwrap(), (2, 0.3));
        assert_eq!(b.select(&all, 1.0).unwrap(), (3, 0.4));
        // Unavailable arms are left out and the rest renormalised.
        let some = [true, false, true, false];
        assert_eq!(b.select(&some, 0.9).unwrap(), (2, 0.75));
        assert_eq!(b.select(&some, 0.1).unwrap(), (0, 0.25));
        assert!(b.select(&[false; ALNS_ARMS], 0.5).is_none());
    }

    #[test]
    fn bandit_learns_the_rewarded_arm() {
        let mut b = AlnsExp3Ix::new();
        let mut rng = S4RandNumGen::new(7);
        let all = [true; ALNS_ARMS];
        for _ in 0..200 {
            let (a, p) = b.select(&all, rng.get_real(0.0, 1.0)).unwrap();
            b.update(a, if a == 2 { 1.0 } else { 0.0 }, p);
        }
        let mut hits = [0; ALNS_ARMS];
        for _ in 0..1000 {
            hits[b.select(&all, rng.get_real(0.0, 1.0)).unwrap().0] += 1;
        }
        assert!(hits[2] > 700, "{hits:?}");
    }

    #[test]
    fn fixing_rate_moves_and_decays() {
        let mut f = AlnsFixingRate::new();
        assert_eq!(f.target, ALNS_MAX_FIXING_RATE);
        f.update(AlnsOutcome::Optimal, false);
        assert!((f.target - 0.7).abs() < 1e-12);
        assert!((f.increment - 0.15).abs() < 1e-12);
        f.update(AlnsOutcome::NodeLimit, false);
        assert!((f.target - 0.85).abs() < 1e-12);
        // A node limit that found something keeps the rate.
        f.update(AlnsOutcome::NodeLimit, true);
        assert!((f.target - 0.85).abs() < 1e-12);
        f.update(AlnsOutcome::Other, false);
        assert!((f.target - 0.85).abs() < 1e-12);
        for _ in 0..50 {
            f.update(AlnsOutcome::Infeasible, false);
        }
        assert_eq!(f.target, ALNS_MIN_FIXING_RATE);
        assert_eq!(f.increment, ALNS_LRATE_MIN);
    }

    #[test]
    fn reward_matches_scip_formula() {
        // Failure: 0.5 (1 - 40 * 0.5 / 50) = 0.3.
        assert!((alns_reward(40, 0.5, false, 10.0, 10.0, 0.0, 50) - 0.3).abs() < 1e-12);
        // Failure far past the target scores 0.
        assert_eq!(alns_reward(400, 0.0, false, 10.0, 10.0, 0.0, 50), 0.0);
        // Success closing half the gap: effort 0.2, 0.5 + 0.5 (0.8 + 0.2 0.5)/1.2.
        let want = 0.5 + 0.5 * 0.9 / 1.2;
        assert!((alns_reward(40, 0.5, true, 10.0, 5.0, 0.0, 50) - want).abs() < 1e-12);
        // First solution: full gap credit and no effort.
        assert!((alns_reward(0, 0.5, true, f64::INFINITY, 5.0, 0.0, 50) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn limit_arithmetic_with_the_work_clamp() {
        // Nothing run yet, 100 parent nodes: (0.1 * 100 + 500) * 1/5 = 102.
        let (max, stall) =
            alns_determine_limits(0, 0, 4, 0, 100, ALNS_MIN_NODES, f64::INFINITY).unwrap();
        assert_eq!((max, stall), (102, 50));
        // Two calls, one success, all arms initialised, 10000 nodes, 300 used:
        // quot = 0.1 * 2 / 3; (666.67 + 500 - 300 - 200) = 666.67 * 5/5 -> 666.
        let (max, _) =
            alns_determine_limits(1, 2, 0, 300, 10_000, ALNS_MIN_NODES, f64::INFINITY).unwrap();
        assert_eq!(max, 666);
        // The cap at 5000 nodes.
        let (max, _) =
            alns_determine_limits(5, 1, 0, 0, 1_000_000, ALNS_MIN_NODES, f64::INFINITY).unwrap();
        assert_eq!(max, 5000);
        // The work clamp binds and can switch the call off.
        let (max, _) = alns_determine_limits(0, 0, 4, 0, 100, ALNS_MIN_NODES, 80.0).unwrap();
        assert_eq!(max, 80);
        assert!(alns_determine_limits(0, 0, 4, 0, 100, ALNS_MIN_NODES, 49.0).is_none());
        // Budget exhausted by earlier calls.
        assert!(alns_determine_limits(0, 5, 0, 900, 100, ALNS_MIN_NODES, f64::INFINITY).is_none());
    }

    #[test]
    fn work_gate_and_clamp_in_the_state() {
        let mut st = AlnsState::new(0);
        st.initialised = true;
        st.offset_work = 0.0;
        // Before the first call only the node formula bounds the call: (0.1 * 100 + 500)/5.
        assert_eq!(st.work_per_node(), None);
        assert_eq!(st.limits(0, 100), Some((102, 50)));
        // One call of 100 nodes and 1000 work: 10 work per node. Allowance
        // 0.1 * 10000 + 0 = 1000 is spent exactly, so the cap is 0 nodes.
        st.alns_work = 1000;
        st.used_nodes = 100;
        st.calls = 1;
        st.runs[0] = 1;
        assert_eq!(st.work_per_node(), Some(10.0));
        assert_eq!(st.limits(10_000, 100_000), None);
        // Allowance 0.1 * 20000 = 2000: 1000 left = 100 nodes at 10 work per node; the node
        // formula (0.05 * 100000 + 500 - 100 - 100) * 2/5 = 2000 is above, so the clamp binds.
        assert_eq!(st.limits(20_000, 100_000), Some((100, 50)));
        // 0.1 * 10500 = 1050 allowance: 5 nodes, below the target 50.
        assert_eq!(st.limits(10_500, 100_000), None);
        // The offset adds to the allowance.
        st.offset_work = 500.0;
        assert_eq!(st.limits(10_500, 100_000), Some((55, 50)));
        // Spent beyond the allowance: off.
        st.alns_work = 5000;
        assert_eq!(st.limits(10_500, 100_000), None);
    }

    #[test]
    fn observe_updates_ledger_rate_target_and_bandit() {
        let mut st = AlnsState::new(0);
        for k in 0..4 {
            st.observe(
                k,
                0.25,
                AlnsOutcome::NodeLimit,
                40,
                400,
                0.5,
                false,
                10.0,
                10.0,
                0.0,
            );
        }
        assert_eq!(st.calls, 4);
        assert_eq!(st.used_nodes, 160);
        assert_eq!(st.alns_work, 1600);
        assert!(st.runs.iter().all(|&r| r == 1));
        // The fourth call (calls >= K) raised the target: 50 * 1.05 = 52, + 1.
        assert_eq!(st.target_nodes, 53);
        assert_eq!(st.found, 0);
        st.observe(
            2,
            0.25,
            AlnsOutcome::Optimal,
            5,
            50,
            0.9,
            true,
            10.0,
            8.0,
            0.0,
        );
        assert_eq!((st.found, st.found_by_arm[2]), (1, 1));
        assert!(st.rates[2].target < ALNS_MAX_FIXING_RATE);
    }

    // ---- solver-level tests -------------------------------------------------------------

    /// `m`-row 0/1 multi-knapsack, `max v.x`, weights and values from a fixed LCG; hard enough
    /// that the root does not close it.
    fn multi_knapsack(n: usize, m: usize) -> HighsLp {
        use crate::lp_data::hconst::{HighsVarType, MatrixFormat, ObjSense};
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as f64
        };
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = m as HighsInt;
        lp.sense_ = ObjSense::kMaximize;
        let mut weights = vec![vec![0.0; n]; m];
        for j in 0..n {
            let mut sum = 0.0;
            for row in weights.iter_mut() {
                row[j] = 10.0 + (next() % 40.0).floor();
                sum += row[j];
            }
            lp.col_cost_
                .push((sum / m as f64 + (next() % 10.0).floor()).floor());
        }
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf; m];
        lp.row_upper_ = weights
            .iter()
            .map(|r| (r.iter().sum::<f64>() / 3.0).floor())
            .collect();
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = m as HighsInt;
        lp.a_matrix_.start_ = (0..=n).map(|j| (j * m) as HighsInt).collect();
        for j in 0..n {
            for (i, row) in weights.iter().enumerate() {
                lp.a_matrix_.index_.push(i as HighsInt);
                lp.a_matrix_.value_.push(row[j]);
            }
        }
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp
    }

    /// What a solve leaves behind for the assertions.
    struct Run {
        obj: f64,
        nodes: i64,
        lp_iters: i64,
        calls: i64,
        found: i64,
        by_source: [i64; 27],
        incumbent: Vec<f64>,
        alns_state: bool,
        s6_state: bool,
        work: crate::work::WorkSnapshot,
    }

    fn solve(alns: bool, crossover: bool, submip: bool, n: usize) -> Run {
        std::thread::spawn(move || {
            crate::parallel::highs_parallel::initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.mip_alns = alns;
            options.mip_heuristic_run_crossover = crossover;
            let lp = multi_knapsack(n, 3);
            let mut callback = HighsCallback::default();
            let solution = HighsSolution::default();
            let level = HighsInt::from(submip);
            let mut ms =
                HighsMipSolver::new(&mut callback, &options, &lp, &solution, submip, level);
            ms.run();
            let r = Run {
                obj: ms.mipdata_.upper_bound,
                nodes: ms.mipdata_.num_nodes,
                lp_iters: ms.mipdata_.total_lp_iterations,
                calls: ms.stats_.alns_calls,
                found: ms.stats_.alns_found,
                by_source: ms.stats_.improving_by_source,
                incumbent: ms.mipdata_.incumbent.clone(),
                alns_state: ms.mipdata_.alns.is_some(),
                s6_state: ms.mipdata_.s6.is_some(),
                work: ms.collect_stats().work_parts,
            };
            drop(ms);
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            r
        })
        .join()
        .unwrap()
    }

    const N: usize = 48;

    #[test]
    fn off_creates_nothing_and_matches_other_levers_off() {
        let off = solve(false, false, false, N);
        assert!(!off.alns_state && !off.s6_state);
        assert_eq!((off.calls, off.found), (0, 0));
        assert_eq!(off.by_source[kSolutionSourceAlns as usize + 1], 0);
        // ALNS off with crossover on is the path, untouched by this lever.
        let xo = solve(false, true, false, N);
        assert!(xo.s6_state && !xo.alns_state);
        assert_eq!(xo.calls, 0);
        assert!((off.obj - xo.obj).abs() < 1e-9);
    }

    #[test]
    fn on_reaches_the_same_optimum_and_fires() {
        let off = solve(false, false, false, N);
        let on = solve(true, false, false, N);
        assert!(on.alns_state && on.s6_state);
        assert!(on.calls >= 1, "ALNS never ran");
        assert!(
            (off.obj - on.obj).abs() < 1e-9,
            "off {} on {}",
            off.obj,
            on.obj
        );
    }

    // The work counters are per solve (crates/rusolve-ekk/src/work.rs), so other tests running
    // beside this one do not enter its ALNS budget. The separate-process check is
    // crates/rusolve-cli/tests/alns_determinism.rs.
    #[test]
    fn on_is_deterministic_at_one_thread() {
        let a = solve(true, false, false, N);
        let b = solve(true, false, false, N);
        assert_eq!(a.nodes, b.nodes);
        assert_eq!(a.lp_iters, b.lp_iters);
        assert_eq!((a.calls, a.found), (b.calls, b.found));
        assert_eq!(a.by_source, b.by_source);
        assert_eq!(a.obj.to_bits(), b.obj.to_bits());
        assert_eq!(a.incumbent, b.incumbent);
    }

    // Two solves at once in one process each see only their own work (per-solve counters).
    #[test]
    fn concurrent_solves_have_their_own_work() {
        let alone_a = solve(true, false, false, N);
        let alone_b = solve(true, false, false, N - 8);
        assert!(alone_a.work.work() > 0 && alone_b.work.work() > 0);
        assert_ne!(alone_a.work, alone_b.work);
        let ha = std::thread::spawn(|| solve(true, false, false, N));
        let hb = std::thread::spawn(|| solve(true, false, false, N - 8));
        let (a, b) = (ha.join().unwrap(), hb.join().unwrap());
        assert_eq!(a.work, alone_a.work);
        assert_eq!(b.work, alone_b.work);
        assert_eq!(
            (a.nodes, a.calls, a.found),
            (alone_a.nodes, alone_a.calls, alone_a.found)
        );
        assert_eq!(
            (b.nodes, b.calls, b.found),
            (alone_b.nodes, alone_b.calls, alone_b.found)
        );
    }

    #[test]
    fn skipped_inside_a_submip() {
        let sub = solve(true, false, true, N);
        assert!(!sub.alns_state);
        assert_eq!(sub.calls, 0);
    }

    #[test]
    fn arm_availability() {
        std::thread::spawn(|| {
            crate::parallel::highs_parallel::initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.mip_alns = true;
            let lp = multi_knapsack(N, 3);
            let mut callback = HighsCallback::default();
            let solution = HighsSolution::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let mut st = AlnsState::new(0);
            // After the solve the root domain is closed; reopen it to the model's bounds so
            // the neighbourhoods see the situation of a hook.
            let (lo, up) = (
                mip_model!(ms).col_lower_.clone(),
                mip_model!(ms).col_upper_.clone(),
            );
            ms.mipdata_.domains[0].col_lower_ = lo;
            ms.mipdata_.domains[0].col_upper_ = up;
            assert!(!ms.mipdata_.incumbent.is_empty());
            // RENS needs only the root relaxation; its restriction sits near the target rate.
            let rens = alns_neighbourhood(&mut ms, &mut st, ALNS_RENS).expect("RENS");
            assert!(rens.rate() >= 0.5 * st.rates[ALNS_RENS].target);
            // With an incumbent RINS and mutation run.
            let mutation = alns_neighbourhood(&mut ms, &mut st, ALNS_MUTATION).expect("mutation");
            let want = (0.9 * ms.mipdata_.integral_cols.len() as f64) as usize + 1;
            assert_eq!(mutation.nfixed, want);
            // One pool point: crossover cannot run.
            if let Some(s6) = ms.mipdata_.s6.as_mut() {
                s6.pool.sols.truncate(1);
            }
            assert!(alns_neighbourhood(&mut ms, &mut st, ALNS_CROSSOVER).is_none());
            // No incumbent: only RENS is left.
            ms.mipdata_.incumbent.clear();
            assert!(alns_neighbourhood(&mut ms, &mut st, ALNS_RINS).is_none());
            assert!(alns_neighbourhood(&mut ms, &mut st, ALNS_MUTATION).is_none());
            assert!(alns_neighbourhood(&mut ms, &mut st, ALNS_RENS).is_some());
            // No relaxation point: nothing fixes from it.
            ms.mipdata_.rootlpsol.clear();
            assert!(alns_neighbourhood(&mut ms, &mut st, ALNS_RENS).is_none());
            drop(ms);
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }

    #[test]
    fn classify_outcomes() {
        assert_eq!(
            alns_classify(false, HighsModelStatus::kInfeasible),
            AlnsOutcome::Infeasible
        );
        assert_eq!(
            alns_classify(true, HighsModelStatus::kOptimal),
            AlnsOutcome::Optimal
        );
        assert_eq!(
            alns_classify(true, HighsModelStatus::kSolutionLimit),
            AlnsOutcome::NodeLimit
        );
        assert_eq!(
            alns_classify(false, HighsModelStatus::kTimeLimit),
            AlnsOutcome::Other
        );
    }
}
