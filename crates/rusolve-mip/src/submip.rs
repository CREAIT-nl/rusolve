//! Sub-MIP primal heuristics: RINS and RENS.
//!
//! Both restrict the problem to a neighbourhood of a reference point by
//! fixing (RINS) or boxing (RENS) integer columns, then solve the restriction
//! with [`bnb::solve`] under a node limit, a time limit and an objective
//! cutoff, map the result back and re-check it on the original problem.
//!
//! * **RINS** fixes integer columns on which the incumbent and an LP solution
//!   agree; the sub-MIP searches the disagreeing ones below the incumbent.
//! * **RENS** boxes every integer column to `[⌊x̄⌋, ⌈x̄⌉]` of an LP solution
//!   `x̄`, which fixes the integral ones; the sub-MIP searches the roundings
//!   of `x̄`. It needs no incumbent.
//!
//! Fixed columns are eliminated before the sub-MIP is built (their row
//! contributions move into the row bounds), rows left empty or redundant
//! under the remaining bounds are dropped, and the cutoff becomes one more
//! row `cᵀx ≤ upper_limit` over the free columns, since [`MipConfig`] has no
//! cutoff field.
//!
//! Provenance:
//!
//! * E. Danna, E. Rothberg & C. Le Pape, "Exploring relaxation induced
//!   neighborhoods to improve MIP solutions", Math. Prog. 102:71–90 (2005)
//!   — RINS: fix the integer variables on which incumbent and LP relaxation
//!   agree, solve the rest as a truncated sub-MIP with the incumbent cutoff.
//! * T. Berthold, "RENS: the optimal rounding", Math. Prog. Comp. 6:33–54
//!   (2014) — RENS: fix integral LP values, bound the fractional ones to
//!   their two roundings, skip when too few columns get fixed.
//! * HiGHS `HighsPrimalHeuristics.cpp` (reference/HiGHS/highs/mip/, MIT,
//!   © 2016–2025 the HiGHS authors; read as a reference text and
//!   reimplemented, not translated):
//!   - `solveSubMip` (:77): the sub-MIP is the original rows under the new
//!     column bounds, with node limit 200 at the call sites (:369, :627,
//!     :938) and the parent's upper limit as objective bound; a sub-MIP
//!     found infeasible is an *infeasible observation* of its fixing rate,
//!     one that improved the upper limit a *success observation* (:229–:243);
//!   - `determineTargetFixingRate` (:249): fixing stops at a target rate
//!     drawn from `[low, high]`, both 0.6 a priori, `high ← 0.9 ·` mean
//!     infeasible rate, `low ← min(low, 0.9 ·` mean success rate`)`,
//!     `high ← max(high, 1.1 ·` mean success rate`)` ([`FixingRateStats`]);
//!   - `RENS` (:394): the rounding box `[⌊x̄ + feastol⌋, ⌈x̄ − feastol⌉]`,
//!     no sub-MIP below fixing rate 0.1 (:608), and a sub-MIP
//!     infeasible at its root is retried with half the fixing rate (:652);
//!   - `RINS` (:661): agreement means `|x̄ − x_inc| ≤ feastol` and fixes to
//!     the nearest integer; below fixing rate 0.25 the neighbourhood is
//!     extended by fixing further integer columns to `calcFixVal` (:377)
//!     values, closest to `x̄` first (:704, :836–:842); no sub-MIP below 0.1
//!     (:919); the same halving retry (:964);
//!   - the upper limit `⌊z − ½⌋ + feastol` for an integral objective and
//!     `z − feastol` otherwise (`HighsMipSolverData::computeNewUpperLimit`,
//!     HighsMipSolverData.cpp:628).
//! * Dropping rows that the column bounds already satisfy: E. D. Andersen &
//!   K. D. Andersen, "Presolving in linear programming", Math. Prog.
//!   71:221–245 (1995), §3 (forcing and redundant constraints).
//!
//! Deviations (deliberate):
//!
//! * Static fixing. HiGHS reaches the neighbourhood by an LP-guided dive
//!   (fix ≥ 10 % more per round, propagate, resolve the LP, backtrack on
//!   infeasibility); here the fixings are chosen once from the given `x̄`, in
//!   a seeded random column order, until the target rate is met. The dive
//!   needs warm-started node LPs to be affordable. With
//!   [`SubMipConfig::propagate_fixings`] each fixing is propagated as in the
//!   dive (`localdom.propagate()` after every branching, e.g. :496, :799),
//!   implied fixings count towards the rate as in `getFixingRate` (:296),
//!   and a fixing that propagates to infeasibility is undone and skipped
//!   where the dive would backtrack; off by default, as measured.
//!   // TODO(perf): revisit the dive once the dual simplex warm start lands
//!   // (measured static numbers on our benchmark set).
//! * No dive fallback below the minimum fixing rate (HiGHS then runs a
//!   10-node depth-first search on the dive tree); the heuristic just skips.
//! * The cutoff is a row, not an objective bound inside the sub-MIP's
//!   pruning, and no presolve, pseudocost or conflict information is shared
//!   between parent and sub-MIP.
//! * The returned point is re-checked with [`MipProblem::violation`] and
//!   must be strictly better than the cutoff value it was asked to beat.

use crate::bnb::{self, MipConfig, MipStatus};
use crate::problem::MipProblem;
use crate::propagate::Propagator;
use rusolve_core::CscMatrix;
use rusolve_lp::rng::XorShift64;

/// Absolute violation accepted on the mapped-back point (bnb's own
/// candidate tolerance).
const ACCEPT_TOL: f64 = 1e-6;

/// Wall-clock share of the parent *tree* time that sub-MIP calls may spend
/// (charging their seconds, setup included, cumulatively).
///
/// New — not a HiGHS port. HiGHS budgets heuristics in LP iterations only
/// (`moreHeuristicsAllowed`, HighsMipSolverData.cpp:672–731: 10 000 plus
/// `mip_heuristic_effort` × all iterations early, sub-MIP iterations charged
/// back scaled by the free-column share, HighsPrimalHeuristics.cpp:216–222).
/// Iterations cannot see what made our sub-MIPs expensive: each call builds
/// the restriction, re-presolves it and starts a cold engine, so a call that
/// burns no iterations still costs wall time — the iteration-budget-only
/// schedule measured 1.14×/1.20× time at 0.67×/0.76× nodes
/// on our benchmark set, and under a 10 ms cap
/// `int_production` s1 still ran 90 setup-dominated RINS calls (~1 s, none
/// improving). Charging seconds and demanding the same `effort` share
/// (0.05, the HiGHS default of `mip_heuristic_effort`) of the time the
/// parent tree itself has spent keeps the search time near `1 + share` of
/// the no-sub-MIP time while self-scaling: the budget grows only as the
/// parent does. Call *setup* is not capped (the time limit only reaches the
/// inner solve) — it is charged after, so an overshooting call pushes the
/// budget negative and the next one waits for the parent to earn it back;
/// a single big-LP call can still overshoot once (air05: 16 s).
pub const SUBMIP_TIME_SHARE: f64 = 0.05;

/// Wall-clock budget a sub-MIP *node* call must still have ahead of it to
/// start at all: below this the restriction build, presolve and engine
/// setup cost more than the truncated search could return. The root RENS
/// is exempt — HiGHS-faithful unconditional (gating it made the whole
/// schedule a no-op).
pub(crate) const SUBMIP_MIN_CALL_SECS: f64 = 0.01;

/// Remaining sub-MIP wall budget under [`SUBMIP_TIME_SHARE`]:
/// `share × parent-only time − spent`, negative when over budget.
/// `submip_secs` is the search's accumulated sub-MIP time, `elapsed_secs`
/// the search's total (sub-MIP time included). Gates and caps the node
/// calls only; the root RENS is unconditional but its seconds are charged,
/// so a costly root call spends the node calls' budget before they start.
pub(crate) fn time_share_left(submip_secs: f64, elapsed_secs: f64) -> f64 {
    SUBMIP_TIME_SHARE * (elapsed_secs - submip_secs).max(0.0) - submip_secs
}

#[derive(Debug, Clone, Copy)]
pub struct SubMipConfig {
    /// Node limit of the sub-MIP (HiGHS: 200 + parent nodes / 20).
    pub node_limit: u64,
    /// Wall-clock seconds for the sub-MIP.
    pub time_limit: Option<f64>,
    /// Fixing stops once this fraction of the unfixed integer columns is
    /// fixed; the caller draws it from [`FixingRateStats::target`].
    pub target_fixing_rate: f64,
    /// No sub-MIP below this fixing rate: the restriction would be almost
    /// as hard as the problem.
    pub min_fixing_rate: f64,
    /// RINS extends its neighbourhood with rounding fixings up to this rate.
    pub rins_min_fixing_rate: f64,
    /// Internal (minimization) objective the RENS solution must beat; RINS
    /// always uses its incumbent.
    pub cutoff: Option<f64>,
    /// Seeds the column order in which fixings are applied.
    pub seed: u64,
    /// Propagate every fixing and count implied fixings towards the rate,
    /// as HiGHS's dive does. Off by default: on gap it shrank the
    /// neighbourhood to whole jobs and lost closure at equal time.
    pub propagate_fixings: bool,
    /// Tolerances and LP settings of the sub-MIP (its limits are replaced).
    pub mip: MipConfig,
}

impl Default for SubMipConfig {
    fn default() -> Self {
        Self {
            node_limit: 200,
            time_limit: None,
            target_fixing_rate: 0.6,
            min_fixing_rate: 0.1,
            rins_min_fixing_rate: 0.25,
            cutoff: None,
            seed: 1,
            propagate_fixings: false,
            mip: MipConfig::default(),
        }
    }
}

/// How a sub-MIP run ended.
#[derive(Debug, Clone, PartialEq)]
pub enum SubMipOutcome {
    /// Fixing rate below `min_fixing_rate` (or nothing to fix): not solved.
    Skipped,
    /// The fixings alone violate a row: infeasible without a search.
    FixingsInfeasible,
    /// The sub-MIP ran and ended with this status.
    Solved(MipStatus),
    /// The sub-MIP failed (LP error); nothing is reported.
    Failed,
}

#[derive(Debug, Clone)]
pub struct SubMipReport {
    /// Verified, strictly improving point of the original problem.
    pub x: Option<Vec<f64>>,
    /// Fixed integer columns / integer columns not fixed in `p`.
    pub fixing_rate: f64,
    pub outcome: SubMipOutcome,
    /// Columns and rows of the sub-MIP actually solved.
    pub sub_cols: usize,
    pub sub_rows: usize,
    pub nodes: u64,
    pub lp_iterations: u64,
    /// Fixing rates of the attempts in this call whose restriction was
    /// infeasible (the last attempt included), for [`FixingRateStats`].
    pub infeasible_rates: Vec<f64>,
}

/// RINS from a feasible `incumbent` and an LP solution `lp_x`.
pub fn rins(
    p: &MipProblem,
    incumbent: &[f64],
    lp_x: &[f64],
    cfg: &SubMipConfig,
) -> Option<Vec<f64>> {
    rins_report(p, incumbent, lp_x, cfg).x
}

/// RENS from an LP solution `lp_x`, beating `cfg.cutoff` when one is given.
pub fn rens(p: &MipProblem, lp_x: &[f64], cfg: &SubMipConfig) -> Option<Vec<f64>> {
    rens_report(p, lp_x, cfg).x
}

pub fn rins_report(
    p: &MipProblem,
    incumbent: &[f64],
    lp_x: &[f64],
    cfg: &SubMipConfig,
) -> SubMipReport {
    let cutoff = p.objective(incumbent);
    let limit = upper_limit(p, cutoff, cfg.mip.feastol);
    with_retries(p, Some(cutoff), cfg, |target| {
        rins_fixing(p, incumbent, lp_x, cfg, target, limit)
    })
}

pub fn rens_report(p: &MipProblem, lp_x: &[f64], cfg: &SubMipConfig) -> SubMipReport {
    let limit = cfg
        .cutoff
        .map_or(f64::INFINITY, |z| upper_limit(p, z, cfg.mip.feastol));
    with_retries(p, cfg.cutoff, cfg, |target| {
        rens_fixing(p, lp_x, cfg, target, limit)
    })
}

fn rins_fixing(
    p: &MipProblem,
    incumbent: &[f64],
    lp_x: &[f64],
    cfg: &SubMipConfig,
    target: f64,
    limit: f64,
) -> Fixing {
    let mut fix = Fixing::new(p, cfg, limit);
    let order = fix.order.clone();
    for &j in &order {
        if fix.rate() >= target {
            break;
        }
        let v = lp_x[j].round();
        if (lp_x[j] - incumbent[j]).abs() <= cfg.mip.feastol && fix.can_fix(j, v) {
            fix.tighten(j, v, v);
        }
    }
    // Extension (HiGHS RINS :836): round further columns, closest to the
    // LP value first, until the neighbourhood is small enough — never past
    // the target, which HiGHS checks after every fixing.
    let extend_to = cfg.rins_min_fixing_rate.min(target);
    if fix.rate() < extend_to {
        let mut rest: Vec<(f64, usize)> = order
            .iter()
            .filter(|&&j| fix.lower[j] < fix.upper[j])
            .map(|&j| {
                let v = clamp_fix(&fix, j, calc_fix_val(lp_x[j], p.cost[j]));
                ((v - lp_x[j]).abs(), j)
            })
            .collect();
        // Stable sort keeps the random order among ties (HiGHS: hash).
        rest.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (_, j) in rest {
            if fix.rate() >= extend_to {
                break;
            }
            let v = clamp_fix(&fix, j, calc_fix_val(lp_x[j], p.cost[j]));
            if fix.can_fix(j, v) {
                fix.tighten(j, v, v);
            }
        }
    }
    fix
}

fn rens_fixing(
    p: &MipProblem,
    lp_x: &[f64],
    cfg: &SubMipConfig,
    target: f64,
    limit: f64,
) -> Fixing {
    let feastol = cfg.mip.feastol;
    let mut fix = Fixing::new(p, cfg, limit);
    let order = fix.order.clone();
    for &j in &order {
        if fix.rate() >= target {
            break;
        }
        let (lo, hi) = (fix.lower[j], fix.upper[j]);
        let down = (lp_x[j] + feastol).floor().clamp(lo, hi);
        let up = (lp_x[j] - feastol).ceil().clamp(lo, hi);
        if lo < down || up < hi {
            fix.tighten(j, down.min(up), up.max(down));
        }
    }
    fix
}

/// Solves the restriction `build(target)`; when it is infeasible at the
/// root, retries with half the fixing rate (HiGHS RENS :652 and RINS :964
/// halve `maxfixingrate` and backtrack the dive to half its depth; here the
/// same column order makes the smaller fixing a prefix of the larger one).
fn with_retries(
    p: &MipProblem,
    cutoff: Option<f64>,
    cfg: &SubMipConfig,
    build: impl Fn(f64) -> Fixing,
) -> SubMipReport {
    let mut target = cfg.target_fixing_rate;
    let mut infeasible_rates = Vec::new();
    loop {
        let mut report = solve_restricted(p, &build(target), cutoff, cfg);
        let root_infeasible = match report.outcome {
            SubMipOutcome::FixingsInfeasible => true,
            SubMipOutcome::Solved(MipStatus::Infeasible) => report.nodes <= 1,
            _ => false,
        };
        if matches!(
            report.outcome,
            SubMipOutcome::FixingsInfeasible | SubMipOutcome::Solved(MipStatus::Infeasible)
        ) {
            infeasible_rates.push(report.fixing_rate);
        }
        let next = 0.5 * report.fixing_rate;
        // The rate must shrink, or the RINS extension could loop.
        if !root_infeasible
            || report.fixing_rate == 0.0
            || next < cfg.min_fixing_rate
            || next >= target
        {
            report.infeasible_rates = infeasible_rates;
            return report;
        }
        target = next;
    }
}

/// HiGHS `calcFixVal` (:377) with no root-LP direction: round towards the
/// worse objective, to the nearest integer when the cost is zero.
fn calc_fix_val(v: f64, cost: f64) -> f64 {
    if cost > 0.0 {
        v.ceil()
    } else if cost < 0.0 {
        v.floor()
    } else {
        (v + 0.5).floor()
    }
}

fn clamp_fix(fix: &Fixing, j: usize, v: f64) -> f64 {
    v.min(fix.upper[j]).max(fix.lower[j])
}

/// Column bounds of the neighbourhood under construction. With
/// `propagate_fixings` every fixing is propagated (rows and the objective
/// cutoff), as HiGHS propagates its local domain after each dive branching,
/// so implied fixings count towards the rate; a fixing whose propagation
/// proves infeasibility is undone and skipped (HiGHS backtracks the dive).
// TODO(perf): bounds are snapshotted and the rate recounted per fixing
// (O(n) each); a domain change stack would make both incremental.
pub(crate) struct Fixing {
    pub(crate) lower: Vec<f64>,
    pub(crate) upper: Vec<f64>,
    /// Integer columns not fixed after the initial propagation, in seeded
    /// random order: the denominator of the fixing rate.
    pub(crate) order: Vec<usize>,
    fixed: usize,
    prop: Option<Propagator>,
    limit: f64,
    /// Propagating `p` under the cutoff alone is infeasible.
    infeasible: bool,
}

impl Fixing {
    pub(crate) fn new(p: &MipProblem, cfg: &SubMipConfig, limit: f64) -> Self {
        let feastol = cfg.mip.feastol;
        let mut lower = p.col_lower.clone();
        let mut upper = p.col_upper.clone();
        for j in 0..p.num_col {
            if p.integer[j] {
                lower[j] = (lower[j] - feastol).ceil();
                upper[j] = (upper[j] + feastol).floor();
            }
        }
        let mut prop = cfg.propagate_fixings.then(|| Propagator::new(p, feastol));
        let infeasible = lower.iter().zip(&upper).any(|(l, u)| l > u)
            || prop
                .as_mut()
                .is_some_and(|pr| !pr.propagate(&mut lower, &mut upper, None, limit));
        let mut order: Vec<usize> = (0..p.num_col)
            .filter(|&j| !infeasible && p.integer[j] && lower[j] < upper[j])
            .collect();
        let mut rng = XorShift64::new(cfg.seed);
        for i in (1..order.len()).rev() {
            let k = (rng.next_u64() % (i as u64 + 1)) as usize;
            order.swap(i, k);
        }
        Self {
            lower,
            upper,
            order,
            fixed: 0,
            prop,
            limit,
            infeasible,
        }
    }

    pub(crate) fn rate(&self) -> f64 {
        if self.order.is_empty() {
            0.0
        } else {
            self.fixed as f64 / self.order.len() as f64
        }
    }

    /// `j` is still free and `v` lies in its domain.
    pub(crate) fn can_fix(&self, j: usize, v: f64) -> bool {
        self.lower[j] < self.upper[j] && self.lower[j] <= v && v <= self.upper[j]
    }

    /// Tightens `j` to `[lo, hi]` (and propagates); returns `false`, with
    /// the bounds restored, if propagation proves that infeasible.
    pub(crate) fn tighten(&mut self, j: usize, lo: f64, hi: f64) -> bool {
        let Some(prop) = self.prop.as_mut() else {
            if self.lower[j] < self.upper[j] && lo >= hi {
                self.fixed += 1;
            }
            self.lower[j] = lo;
            self.upper[j] = hi;
            return true;
        };
        let (saved_lower, saved_upper) = (self.lower.clone(), self.upper.clone());
        self.lower[j] = lo;
        self.upper[j] = hi;
        if !prop.propagate(&mut self.lower, &mut self.upper, Some(j), self.limit) {
            self.lower = saved_lower;
            self.upper = saved_upper;
            return false;
        }
        self.fixed = self
            .order
            .iter()
            .filter(|&&k| self.lower[k] >= self.upper[k])
            .count();
        true
    }
}

pub(crate) fn solve_restricted(
    p: &MipProblem,
    fix: &Fixing,
    cutoff: Option<f64>,
    cfg: &SubMipConfig,
) -> SubMipReport {
    let mut report = SubMipReport {
        x: None,
        fixing_rate: fix.rate(),
        outcome: SubMipOutcome::Skipped,
        sub_cols: 0,
        sub_rows: 0,
        nodes: 0,
        lp_iterations: 0,
        infeasible_rates: Vec::new(),
    };
    if fix.infeasible {
        report.outcome = SubMipOutcome::FixingsInfeasible;
        return report;
    }
    if fix.order.is_empty() || report.fixing_rate < cfg.min_fixing_rate {
        return report;
    }
    let upper_limit = cutoff.map(|z| upper_limit(p, z, cfg.mip.feastol));
    let Some((sub, free)) = restrict(p, &fix.lower, &fix.upper, upper_limit) else {
        report.outcome = SubMipOutcome::FixingsInfeasible;
        return report;
    };
    report.sub_cols = sub.num_col;
    report.sub_rows = sub.num_row;
    let mip = MipConfig {
        node_limit: Some(cfg.node_limit),
        time_limit: cfg.time_limit,
        ..cfg.mip
    };
    let mut x = fix.lower.clone();
    if sub.num_col == 0 {
        // Everything fixed and every row satisfied: the fixings are the
        // point (an empty bnb incumbent would be indistinguishable from none).
        report.outcome = SubMipOutcome::Solved(MipStatus::Optimal);
    } else {
        let sol = match bnb::solve(&sub, &mip) {
            Ok(sol) => sol,
            Err(_) => {
                report.outcome = SubMipOutcome::Failed;
                return report;
            }
        };
        report.outcome = SubMipOutcome::Solved(sol.status);
        report.nodes = sol.nodes;
        report.lp_iterations = sol.lp_iterations;
        if sol.x.is_empty() {
            return report;
        }
        for (k, &j) in free.iter().enumerate() {
            x[j] = sol.x[k];
        }
    }
    let v = p.violation(&x);
    let feasible = v.bound <= ACCEPT_TOL && v.row <= ACCEPT_TOL && v.integrality <= ACCEPT_TOL;
    let improves = cutoff.is_none_or(|z| p.objective(&x) < z);
    if feasible && improves {
        report.x = Some(x);
    }
    report
}

/// Objective value an improving solution must not exceed (HiGHS
/// `computeNewUpperLimit` without gaps).
pub(crate) fn upper_limit(p: &MipProblem, z: f64, feastol: f64) -> f64 {
    let integral = p
        .cost
        .iter()
        .zip(&p.integer)
        .all(|(&c, &int)| c == 0.0 || (int && (c - c.round()).abs() <= 1e-9));
    if integral {
        (z - 0.5).floor() + feastol
    } else {
        (z - feastol).min(z.next_down())
    }
}

/// The restriction of `p` to `[lower, upper]` with fixed columns eliminated
/// and the objective row `cᵀx ≤ upper_limit` appended; also returns the
/// original index of each sub-MIP column. `None` when the fixings violate a
/// row whatever the free columns do.
fn restrict(
    p: &MipProblem,
    lower: &[f64],
    upper: &[f64],
    upper_limit: Option<f64>,
) -> Option<(MipProblem, Vec<usize>)> {
    let free: Vec<usize> = (0..p.num_col).filter(|&j| lower[j] < upper[j]).collect();
    let mut new_col = vec![usize::MAX; p.num_col];
    for (k, &j) in free.iter().enumerate() {
        new_col[j] = k;
    }
    // Per row: fixed activity, and min/max activity of the free part.
    let m = p.num_row;
    let (mut fixed_act, mut min_act, mut max_act) = (vec![0.0; m], vec![0.0; m], vec![0.0; m]);
    let mut entries: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
    let (ptr, idx, val) = (p.matrix.col_ptr(), p.matrix.row_idx(), p.matrix.values());
    for j in 0..p.num_col {
        for k in ptr[j]..ptr[j + 1] {
            let (i, a) = (idx[k] as usize, val[k]);
            if new_col[j] == usize::MAX {
                fixed_act[i] += a * lower[j];
            } else {
                let (lo, hi) = if a > 0.0 {
                    (a * lower[j], a * upper[j])
                } else {
                    (a * upper[j], a * lower[j])
                };
                min_act[i] += lo;
                max_act[i] += hi;
                entries[i].push((new_col[j], a));
            }
        }
    }
    let tol = ACCEPT_TOL;
    let mut trip = Vec::new();
    let (mut row_lower, mut row_upper) = (Vec::new(), Vec::new());
    for i in 0..m {
        let lo = p.row_lower[i] - fixed_act[i];
        let hi = p.row_upper[i] - fixed_act[i];
        // Infinite bound sums are NaN-free: ±∞ + finite stays ±∞.
        if min_act[i] > hi + tol || max_act[i] < lo - tol {
            return None;
        }
        if min_act[i] >= lo - tol && max_act[i] <= hi + tol {
            continue;
        }
        let r = row_lower.len();
        trip.extend(entries[i].iter().map(|&(k, a)| (r, k, a)));
        row_lower.push(lo);
        row_upper.push(hi);
    }
    let cost: Vec<f64> = free.iter().map(|&j| p.cost[j]).collect();
    if let Some(ul) = upper_limit {
        let fixed_obj: f64 = (0..p.num_col)
            .filter(|&j| new_col[j] == usize::MAX)
            .map(|j| p.cost[j] * lower[j])
            .sum();
        let r = row_lower.len();
        trip.extend(
            cost.iter()
                .enumerate()
                .filter(|&(_, &c)| c != 0.0)
                .map(|(k, &c)| (r, k, c)),
        );
        row_lower.push(f64::NEG_INFINITY);
        row_upper.push(ul - fixed_obj);
    }
    let n = free.len();
    let num_row = row_lower.len();
    let sub = MipProblem {
        num_col: n,
        num_row,
        matrix: CscMatrix::from_triplets(num_row, n, &trip),
        row_lower,
        row_upper,
        cost,
        maximize: false,
        col_lower: free.iter().map(|&j| lower[j]).collect(),
        col_upper: free.iter().map(|&j| upper[j]).collect(),
        integer: free.iter().map(|&j| p.integer[j]).collect(),
    };
    Some((sub, free))
}

/// Per-solve sub-MIP counters, reported on [`crate::MipSolution`].
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SubMipStats {
    /// Calls that built and solved a restriction, per heuristic.
    pub rens_calls: u64,
    pub rins_calls: u64,
    /// Calls that produced a strictly improving incumbent.
    pub rens_found: u64,
    pub rins_found: u64,
    /// Calls the fixing rate or an infeasible fixing turned away.
    pub skipped: u64,
    /// Nodes and (HiGHS-adjusted) LP iterations spent in sub-MIPs.
    pub nodes: u64,
    pub lp_iterations: u64,
    /// Wall-clock seconds spent in sub-MIPs, setup included: the quantity
    /// [`SUBMIP_TIME_SHARE`] budgets. Iteration-based budgets cannot see
    /// this cost (each call pays restriction build, presolve and a cold
    /// engine whatever it iterates), which is why it is tracked separately.
    pub seconds: f64,
    /// ALNS ([`crate::alns`]) calls that solved a restriction and calls
    /// that improved the incumbent, per [`crate::alns::Arm`].
    pub alns_calls: [u64; crate::alns::ARMS],
    pub alns_found: [u64; crate::alns::ARMS],
    /// ALNS wall-clock seconds, setup included.
    pub alns_seconds: f64,
}

/// Fixing-rate observations across calls (HiGHS
/// `determineTargetFixingRate`, :249, and the bookkeeping at :229–:243).
#[derive(Debug, Clone, Default)]
pub struct FixingRateStats {
    infeasible_sum: f64,
    infeasible_count: u32,
    success_sum: f64,
    success_count: u32,
}

impl FixingRateStats {
    /// Target fixing rate for the next call.
    pub fn target(&self, rng: &mut XorShift64) -> f64 {
        let (mut low, mut high) = (0.6_f64, 0.6_f64);
        if self.infeasible_count > 0 {
            high = 0.9 * self.infeasible_sum / self.infeasible_count as f64;
            low = low.min(high);
        }
        if self.success_count > 0 {
            let rate = self.success_sum / self.success_count as f64;
            low = low.min(0.9 * rate);
            high = high.max(1.1 * rate);
        }
        rng.uniform(low, high)
    }

    /// Records a finished call: each infeasible restriction says its rate
    /// was too high, an improving solution that the final rate was good.
    pub fn observe(&mut self, report: &SubMipReport) {
        for &rate in &report.infeasible_rates {
            self.infeasible_sum += rate;
            self.infeasible_count += 1;
        }
        if report.x.is_some() {
            self.success_sum += report.fixing_rate;
            self.success_count += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{LpEngine, LpStatus, RusolveEngine};

    /// `rows`: (coefficients, lower, upper).
    fn mip(rows: &[(&[f64], f64, f64)], lower: &[f64], upper: &[f64], cost: &[f64]) -> MipProblem {
        let n = cost.len();
        let mut trip = Vec::new();
        for (i, (a, _, _)) in rows.iter().enumerate() {
            for (j, &v) in a.iter().enumerate() {
                if v != 0.0 {
                    trip.push((i, j, v));
                }
            }
        }
        MipProblem {
            num_col: n,
            num_row: rows.len(),
            matrix: CscMatrix::from_triplets(rows.len(), n, &trip),
            row_lower: rows.iter().map(|r| r.1).collect(),
            row_upper: rows.iter().map(|r| r.2).collect(),
            cost: cost.to_vec(),
            maximize: false,
            col_lower: lower.to_vec(),
            col_upper: upper.to_vec(),
            integer: vec![true; n],
        }
    }

    fn lp_solution(p: &MipProblem) -> Option<Vec<f64>> {
        let mut e = RusolveEngine::new(p, MipConfig::default().lp);
        let out = e.solve(&p.col_lower, &p.col_upper).ok()?;
        (out.status == LpStatus::Optimal).then_some(out.x)
    }

    fn feasible(p: &MipProblem, x: &[f64]) -> bool {
        let v = p.violation(x);
        v.bound <= 1e-6 && v.row <= 1e-6 && v.integrality <= 1e-6
    }

    fn full_fixing() -> SubMipConfig {
        SubMipConfig {
            target_fixing_rate: 1.0,
            min_fixing_rate: 0.0,
            ..SubMipConfig::default()
        }
    }

    #[test]
    fn rins_improves_a_knapsack_incumbent() {
        // max 10a + 13b + 7c + 8d + e  s.t.  4a + 6b + 3c + 5d + e ≤ 10.
        // Incumbent (0,0,1,1,1) = 16; the LP (1, ½, 1, 0, 0) agrees on c
        // only, the extension fixes one more column, and every completion
        // with c = 1 that the sub-MIP can reach beats 16.
        let p = mip(
            &[(&[4.0, 6.0, 3.0, 5.0, 1.0], f64::NEG_INFINITY, 10.0)],
            &[0.0; 5],
            &[1.0; 5],
            &[-10.0, -13.0, -7.0, -8.0, -1.0],
        );
        let inc = [0.0, 0.0, 1.0, 1.0, 1.0];
        let lp_x = lp_solution(&p).unwrap();
        let cfg = SubMipConfig {
            min_fixing_rate: 0.0,
            ..SubMipConfig::default()
        };
        let report = rins_report(&p, &inc, &lp_x, &cfg);
        let x = report.x.expect("RINS should improve 16");
        assert!(feasible(&p, &x));
        assert!(p.objective(&x) < p.objective(&inc));
    }

    #[test]
    fn rens_with_integral_lp_fixes_everything() {
        // LP optimum is integral: RENS fixes all columns, solves an empty
        // sub-MIP and returns the LP point.
        let p = mip(
            &[(&[1.0, 1.0], 1.0, 1.0)],
            &[0.0, 0.0],
            &[1.0, 1.0],
            &[1.0, 2.0],
        );
        let report = rens_report(&p, &[1.0, 0.0], &full_fixing());
        assert_eq!(report.fixing_rate, 1.0);
        assert_eq!(report.sub_cols, 0);
        assert_eq!(report.x, Some(vec![1.0, 0.0]));
    }

    #[test]
    fn infeasible_cutoff_is_detected_without_a_search() {
        // min x0 + x1, x0 + x1 = 1: nothing beats objective 1, and
        // propagating the cutoff row against the equality proves it.
        let p = mip(
            &[(&[1.0, 1.0], 1.0, 1.0)],
            &[0.0, 0.0],
            &[1.0, 1.0],
            &[1.0, 1.0],
        );
        let cfg = SubMipConfig {
            cutoff: Some(1.0),
            propagate_fixings: true,
            ..full_fixing()
        };
        let report = rens_report(&p, &[0.5, 0.5], &cfg);
        assert_eq!(report.outcome, SubMipOutcome::FixingsInfeasible);
        assert_eq!(report.infeasible_rates, vec![0.0]);
        assert!(report.x.is_none());
    }

    #[test]
    fn propagation_counts_implied_fixings() {
        // Fixing x0 = 1 in x0 + x1 + x2 = 1 fixes the other two: rate 1
        // after one fixing, and the point is complete.
        let p = mip(
            &[(&[1.0, 1.0, 1.0], 1.0, 1.0)],
            &[0.0; 3],
            &[1.0; 3],
            &[3.0, 2.0, 1.0],
        );
        let cfg = SubMipConfig {
            propagate_fixings: true,
            ..full_fixing()
        };
        let report = rens_report(&p, &[1.0, 0.0, 0.0], &cfg);
        assert_eq!(report.fixing_rate, 1.0);
        assert_eq!(report.sub_cols, 0);
        assert!(feasible(&p, report.x.as_ref().unwrap()));
    }

    #[test]
    fn root_infeasible_restriction_is_retried_at_half_rate() {
        // Both columns forced to 1 against x0 + x1 = 1 (bypassing
        // propagation, as an LP-infeasible restriction would look), then
        // the retry fixes one column through the normal path.
        let p = mip(
            &[(&[1.0, 1.0], 1.0, 1.0)],
            &[0.0, 0.0],
            &[1.0, 1.0],
            &[1.0, 2.0],
        );
        let cfg = full_fixing();
        let report = with_retries(&p, None, &cfg, |target| {
            let mut fix = Fixing::new(&p, &cfg, f64::INFINITY);
            if target >= 1.0 {
                fix.lower = vec![1.0, 1.0];
                fix.fixed = 2;
            } else {
                let j = fix.order[0];
                fix.tighten(j, 0.0, 0.0);
            }
            fix
        });
        assert_eq!(report.infeasible_rates, vec![1.0]);
        assert_eq!(report.fixing_rate, 0.5);
        assert!(feasible(&p, report.x.as_ref().unwrap()));
        let mut stats = FixingRateStats::default();
        stats.observe(&report);
        let mut rng = XorShift64::new(3);
        let t = stats.target(&mut rng);
        // Infeasible at 1.0, success at 0.5: [min(0.6, 0.45), max(0.9, 0.55)].
        assert!((0.45..=0.9).contains(&t), "{t}");
    }

    #[test]
    fn time_share_gate() {
        // Nothing spent yet: the full share of 0 elapsed time is 0.
        assert_eq!(time_share_left(0.0, 0.0), 0.0);
        // A root RENS of 5 ms must be amortized by 100 ms of parent tree
        // time at share 0.05 before any node call fits the minimum budget.
        assert!(time_share_left(0.005, 0.01) < SUBMIP_MIN_CALL_SECS);
        assert!(time_share_left(0.005, 0.105) < SUBMIP_MIN_CALL_SECS);
        // The budget grows only with parent (non-sub-MIP) time: 1 s of it
        // funds 50 ms of sub-MIPs.
        assert!((time_share_left(0.005, 1.005) - 0.045).abs() < 1e-12);
        assert!((time_share_left(0.05, 1.05) - 0.0).abs() < 1e-12);
        assert!(time_share_left(0.06, 1.06) < 0.0);
    }

    #[test]
    fn low_fixing_rate_skips() {
        let p = mip(
            &[(&[1.0, 1.0, 1.0, 1.0], 0.0, 2.0)],
            &[0.0; 4],
            &[1.0; 4],
            &[-1.0; 4],
        );
        // Fractional everywhere: RENS fixes nothing.
        let report = rens_report(&p, &[0.5; 4], &SubMipConfig::default());
        assert_eq!(report.outcome, SubMipOutcome::Skipped);
        assert_eq!(report.fixing_rate, 0.0);
    }

    #[test]
    fn target_rate_follows_observations() {
        let mut rng = XorShift64::new(9);
        let stats = FixingRateStats::default();
        assert_eq!(stats.target(&mut rng), 0.6);
        let stats = FixingRateStats {
            success_sum: 0.8,
            success_count: 1,
            ..FixingRateStats::default()
        };
        for _ in 0..100 {
            let t = stats.target(&mut rng);
            assert!((0.6..0.88).contains(&t), "{t}");
        }
    }

    /// Random pure-integer problems (same shape as bnb's tests) with an
    /// enumeration oracle: whatever RINS/RENS return is feasible, never
    /// better than the optimum and strictly better than the cutoff.
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
            let (lo, hi) = match rng.next_u64() % 3 {
                0 => (-inf, b),
                1 => (b, inf),
                _ => (b - 2.0, b + 1.5),
            };
            row_lower.push(lo);
            row_upper.push(hi);
        }
        let col_lower: Vec<f64> = (0..n).map(|_| pick(rng, &[-2.0, -1.0, 0.0, 0.0])).collect();
        let col_upper = col_lower
            .iter()
            .map(|&l| l + pick(rng, &[1.0, 2.0, 3.0]))
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
            maximize: false,
            col_lower,
            col_upper,
            integer: vec![true; n],
        }
    }

    /// All feasible points of a small integer box, by enumeration.
    fn feasible_points(p: &MipProblem) -> Vec<Vec<f64>> {
        let mut x = p.col_lower.clone();
        let mut out = Vec::new();
        loop {
            if p.violation(&x).row <= 1e-9 {
                out.push(x.clone());
            }
            let mut j = 0;
            loop {
                if j == p.num_col {
                    return out;
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
    fn random_problems_against_enumeration() {
        let mut rng = XorShift64::new(20261022);
        let (mut rins_found, mut rens_found, mut checked) = (0, 0, 0);
        for case in 0..300 {
            let n = 3 + (rng.next_u64() % 4) as usize;
            let m = 1 + (rng.next_u64() % 3) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let points = feasible_points(&p);
            let Some(lp_x) = lp_solution(&p) else {
                continue;
            };
            if points.is_empty() {
                continue;
            }
            checked += 1;
            let opt = points
                .iter()
                .map(|x| p.objective(x))
                .fold(f64::INFINITY, f64::min);
            let cfg = SubMipConfig {
                target_fixing_rate: [0.3, 0.6, 1.0][case % 3],
                min_fixing_rate: 0.0,
                seed: case as u64 + 1,
                propagate_fixings: case % 2 == 1,
                ..SubMipConfig::default()
            };
            // Worst feasible point as incumbent: leaves room to improve.
            let worst = points
                .iter()
                .max_by(|a, b| p.objective(a).total_cmp(&p.objective(b)))
                .unwrap();
            if let Some(x) = rins(&p, worst, &lp_x, &cfg) {
                rins_found += 1;
                assert!(feasible(&p, &x), "case {case}");
                let z = p.objective(&x);
                assert!(z < p.objective(worst), "case {case}");
                assert!(z >= opt - 1e-9, "case {case}: {z} < {opt}");
            }
            let report = rens_report(&p, &lp_x, &cfg);
            if let Some(x) = report.x {
                rens_found += 1;
                assert!(feasible(&p, &x), "case {case}");
                let z = p.objective(&x);
                assert!(z >= opt - 1e-9, "case {case}: {z} < {opt}");
            }
            // With a cutoff at the optimum nothing can be returned.
            let tight = SubMipConfig {
                cutoff: Some(opt),
                ..cfg
            };
            assert!(rens(&p, &lp_x, &tight).is_none(), "case {case}");
            let best = points
                .iter()
                .min_by(|a, b| p.objective(a).total_cmp(&p.objective(b)))
                .unwrap();
            assert!(rins(&p, best, &lp_x, &cfg).is_none(), "case {case}");
        }
        assert!(checked > 100, "{checked}");
        assert!(
            rins_found > 20 && rens_found > 20,
            "{rins_found}/{rens_found}"
        );
    }

    #[test]
    fn rens_full_box_covers_every_rounding() {
        // RENS (full box, generous limit) against the best rounding found
        // by enumeration inside the box.
        let mut rng = XorShift64::new(77);
        let mut compared = 0;
        for case in 0..200 {
            let n = 3 + (rng.next_u64() % 4) as usize;
            let m = 1 + (rng.next_u64() % 3) as usize;
            let p = random_pure_integer(&mut rng, n, m);
            let Some(lp_x) = lp_solution(&p) else {
                continue;
            };
            let mut boxed = p.clone();
            for (j, &v) in lp_x.iter().enumerate() {
                boxed.col_lower[j] = (v + 1e-6).floor().max(p.col_lower[j]);
                boxed.col_upper[j] = (v - 1e-6).ceil().min(p.col_upper[j]);
            }
            let best = feasible_points(&boxed)
                .iter()
                .map(|x| p.objective(x))
                .fold(f64::INFINITY, f64::min);
            let cfg = SubMipConfig {
                node_limit: 100_000,
                propagate_fixings: case % 2 == 1,
                ..full_fixing()
            };
            // Tightenings whose propagation fails are skipped, so the
            // restriction contains every feasible point of the box: an exact
            // sub-MIP is at least as good as the best rounding.
            let report = rens_report(&p, &lp_x, &cfg);
            if let Some(x) = &report.x {
                assert!(feasible(&p, x), "case {case}");
            }
            // Skipped only when propagation alone fixes every integer
            // column (rate 0/0, HiGHS skips that too).
            if best.is_finite() && report.outcome != SubMipOutcome::Skipped {
                let z = report.x.as_ref().map(|x| p.objective(x));
                assert!(
                    z.is_some_and(|z| z <= best + 1e-6),
                    "case {case}: {z:?} vs {best}"
                );
            }
            compared += 1;
        }
        assert!(compared > 100, "{compared}");
    }
}
