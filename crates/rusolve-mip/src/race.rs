//! Racing: `threads` complete tree searches of
//! the same presolved problem run concurrently and differ in their random
//! seed, their row/column order and a few strategy switches. They share
//! incumbents, and the first search that proves a result stops the others.
//! Racers other than the default one start when its tree outlives the root
//! node, or after `RAMP_UP`.
//!
//! Sources:
//! * Fischetti & Monaci, "Exploiting erraticism in search", Operations
//!   Research 62(1), 2014. Tree size and solve time vary wildly with
//!   irrelevant changes (row/column order, seeds), so racing K perturbed
//!   copies and keeping the first to finish exploits that variance rather
//!   than suffering it.
//! * Carvajal, Ahmed, Nemhauser, Furman, Goel & Shao, "Using diversification,
//!   communication and parallelism to solve mixed-integer linear programs",
//!   Operations Research Letters 42(2), 2014. Concurrent B&B with different
//!   parameter settings plus incumbent sharing beats both the default and
//!   the best single setting, and sharing incumbents is most of the gain.
//!
//! This is not a HiGHS port: HiGHS 1.15 runs one tree (its threads parallelise
//! LP and presolve work), so racing is the first item of the beyond-HiGHS
//! work.
//!
//! The result depends on thread timing, so a multi-threaded solve is not
//! deterministic. Statuses and objective values are deterministic up to
//! tolerances, but `x`, node counts and limit-time bounds are not. With
//! `threads = 1` the solve does not come through here and is unchanged.

use crate::bnb::{self, MipConfig, MipError, MipSolution, MipStatus};
use crate::heuristics::Source;
use crate::nodequeue::NodeSelection;
use crate::problem::MipProblem;
use crate::pseudocost::DEFAULT_MIN_RELIABLE;
use rusolve_lp::rng::XorShift64;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex, PoisonError};
use std::time::{Duration, Instant};

/// Longest wait of racers 1.. for racer 0's tree to outlive its root
/// (ramp-up). Instances solved at the root then never start racers (they
/// cost ~5 ms per instance when released at the root regardless, on
/// bin_packing 5 → 10 ms). On the 2026-10-22 mixed
/// tier, racing every racer from the start cost 1.2–1.6× time on the
/// 1-node families (p_median 0.18 → 0.25 s, network_design[fixed] s1
/// 1.16 → 1.48 s). The cap keeps big roots (air05: one node in 60 s)
/// diversified: racers' heuristics found 48369 against 50263 there.
// TODO(perf): 1 s is a guess between those two cases; tune on the miplib tier.
const RAMP_UP: Duration = Duration::from_secs(1);

/// State shared by all racers.
pub(crate) struct Shared {
    stop: AtomicBool,
    /// Bumped on every improvement of `best`, under its lock: racers poll it
    /// at each node and only lock when it changed.
    version: AtomicU64,
    /// Best incumbent (internal objective, point in the shared column order).
    best: Mutex<Option<(f64, Vec<f64>)>>,
    /// Ramp-up gate: racers other than racer 0 start once it is set.
    released: Mutex<bool>,
    wake: Condvar,
}

impl Shared {
    fn new() -> Self {
        Self {
            stop: AtomicBool::new(false),
            version: AtomicU64::new(0),
            best: Mutex::new(None),
            released: Mutex::new(false),
            wake: Condvar::new(),
        }
    }

    /// Opens the ramp-up gate.
    fn release(&self) {
        *self.released.lock().unwrap_or_else(PoisonError::into_inner) = true;
        self.wake.notify_all();
    }

    /// Waits until the gate opens or `timeout` passes.
    fn wait_release(&self, timeout: Duration) {
        let open = self.released.lock().unwrap_or_else(PoisonError::into_inner);
        drop(
            self.wake
                .wait_timeout_while(open, timeout, |open| !*open)
                .unwrap_or_else(PoisonError::into_inner),
        );
    }

    fn best(&self) -> std::sync::MutexGuard<'_, Option<(f64, Vec<f64>)>> {
        self.best.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

/// One racer's handle on the race: incumbent exchange through its column
/// permutation, and the stop flag.
pub(crate) struct Link<'a> {
    shared: &'a Shared,
    /// Column `k` of the racer's problem is column `perm[k]` of the shared one.
    perm: Option<&'a [usize]>,
    /// Exchange incumbents (off for the zero-objective feasibility search,
    /// whose objective values do not compare with the race's).
    share: bool,
    /// `Shared::version` at the last pull.
    seen: u64,
    /// Racer 0, whose first node after the root starts the others.
    leader: bool,
}

impl<'a> Link<'a> {
    pub(crate) fn stopped(&self) -> bool {
        self.shared.stop.load(Ordering::Relaxed)
    }

    /// A link to the same race that only observes the stop flag.
    pub(crate) fn stop_only(&self) -> Link<'a> {
        Link {
            shared: self.shared,
            perm: None,
            share: false,
            seen: 0,
            leader: false,
        }
    }

    /// Called when the search takes its first node after the root: racer 0
    /// then starts the others.
    pub(crate) fn root_done(&self) {
        if self.leader {
            self.shared.release();
        }
    }

    /// Offers the racer's new incumbent `x` (objective `z`) to the others.
    pub(crate) fn publish(&self, z: f64, x: &[f64]) {
        if !self.share {
            return;
        }
        let mut best = self.shared.best();
        if best.as_ref().is_some_and(|b| b.0 <= z) {
            return;
        }
        let x = match self.perm {
            Some(perm) => {
                let mut y = vec![0.0; x.len()];
                for (k, &j) in perm.iter().enumerate() {
                    y[j] = x[k];
                }
                y
            }
            None => x.to_vec(),
        };
        *best = Some((z, x));
        self.shared.version.fetch_add(1, Ordering::Release);
    }

    /// The shared incumbent in the racer's column order, if it changed since
    /// the last pull and beats `own` (the racer's incumbent value).
    pub(crate) fn pull(&mut self, own: Option<f64>) -> Option<Vec<f64>> {
        if !self.share || self.shared.version.load(Ordering::Acquire) == self.seen {
            return None;
        }
        let best = self.shared.best();
        self.seen = self.shared.version.load(Ordering::Relaxed);
        let (z, x) = best.as_ref()?;
        if own.is_some_and(|o| o <= *z) {
            return None;
        }
        Some(match self.perm {
            Some(perm) => perm.iter().map(|&j| x[j]).collect(),
            None => x.clone(),
        })
    }
}

/// A racer's problem and settings.
struct Racer {
    config: MipConfig,
    /// Permuted copy of the problem and its column permutation (`None` =
    /// the shared problem as is).
    permuted: Option<(MipProblem, Vec<usize>, Vec<usize>)>,
}

/// One-knob deviations from the base config, in racer order (racer `k ≥ 1`
/// takes `PORTFOLIO[(k - 1) % len]`). Chosen and ordered greedily on the
/// 2026-10-23 knob sweep (medium + mixed families,
/// training seeds 11–15): each entry
/// is the variant that most lowers the shifted geometric mean of the
/// per-instance best time over the racers before it — the greedy portfolio
/// construction of Xu, Hoos & Leyton-Brown, "Hydra: Automatically
/// configuring algorithms for portfolio-based selection", AAAI 2010. Racer
/// counts 2/3/4 cover 0.84/0.74/0.69× the default's time on that sweep
/// (0.91/0.82/0.73× when each generator is held out of the choice), against
/// 0.89/0.83/0.82× for the earlier best-bound / mod-k / child-rule cycle.
const PORTFOLIO: [Switch; 7] = [
    Switch::Gomory,
    Switch::NoStrongBranching,
    Switch::Modk,
    Switch::TreeCuts,
    Switch::NoCuts,
    Switch::NodeSelection,
    Switch::MoreStrongBranching,
];

/// A strategy switch, applied as a move away from the base config so that a
/// base that already has the setting still gets a different racer.
#[derive(Clone, Copy)]
enum Switch {
    Gomory,
    /// `min_reliable = 0`: pseudocost branching without strong branching.
    NoStrongBranching,
    Modk,
    TreeCuts,
    NoCuts,
    NodeSelection,
    /// Twice the base `min_reliable`.
    MoreStrongBranching,
}

impl Switch {
    fn apply(self, config: &mut MipConfig) {
        let b = &mut config.branching;
        match self {
            Switch::Gomory => config.cuts.gomory = !config.cuts.gomory,
            Switch::NoStrongBranching => {
                b.min_reliable = if b.min_reliable == 0 {
                    DEFAULT_MIN_RELIABLE
                } else {
                    0
                }
            }
            Switch::Modk => config.cuts.modk = !config.cuts.modk,
            Switch::TreeCuts => config.cuts.tree = !config.cuts.tree,
            Switch::NoCuts => config.cuts.enabled = !config.cuts.enabled,
            Switch::NodeSelection => {
                b.node_selection = match b.node_selection {
                    NodeSelection::Hybrid => NodeSelection::BestBound,
                    NodeSelection::BestBound => NodeSelection::Hybrid,
                }
            }
            Switch::MoreStrongBranching => {
                b.min_reliable = (2 * b.min_reliable).max(DEFAULT_MIN_RELIABLE)
            }
        }
    }
}

/// Racer `k` of the portfolio. Racer 0 is `base` itself, so a race never
/// does worse than its default search except through lost CPU. The others
/// change seed and row/column order (erraticism, Fischetti & Monaci 2014)
/// and apply one [`PORTFOLIO`] switch each (diversification, Carvajal et
/// al. 2014).
fn racer(p: &MipProblem, base: &MipConfig, k: usize) -> Racer {
    let mut config = *base;
    config.threads = 1;
    if k == 0 {
        return Racer {
            config,
            permuted: None,
        };
    }
    config.seed = base.seed.wrapping_add(k as u64);
    PORTFOLIO[(k - 1) % PORTFOLIO.len()].apply(&mut config);
    let mut rng = XorShift64::new(config.seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));
    let cols = shuffled(p.num_col, &mut rng);
    let rows = shuffled(p.num_row, &mut rng);
    Racer {
        config,
        permuted: Some((p.permuted(&cols, &rows), cols, rows)),
    }
}

/// A uniformly random permutation of `0..n` (Fisher–Yates).
fn shuffled(n: usize, rng: &mut XorShift64) -> Vec<usize> {
    let mut v: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = (rng.next_u64() % (i as u64 + 1)) as usize;
        v.swap(i, j);
    }
    v
}

/// Solves `p` (already presolved, or with presolve off) with
/// `config.threads` racers.
pub(crate) fn solve(p: &MipProblem, config: &MipConfig) -> Result<MipSolution, MipError> {
    let start = Instant::now();
    let shared = Shared::new();
    let winner = AtomicUsize::new(usize::MAX);
    let results: Vec<Option<Result<MipSolution, MipError>>> = std::thread::scope(|s| {
        let (shared, winner) = (&shared, &winner);
        let handles: Vec<_> = (1..config.threads)
            .map(|k| {
                s.spawn(move || {
                    shared.wait_release(RAMP_UP);
                    if shared.stop.load(Ordering::Relaxed) {
                        return None;
                    }
                    let mut r = racer(p, config, k);
                    // Late starters get what is left of the time limit.
                    r.config.time_limit = config
                        .time_limit
                        .map(|t| (t - start.elapsed().as_secs_f64()).max(0.0));
                    Some(run_racer(p, &r, k, shared, winner))
                })
            })
            .collect();
        // Racer 0 runs on the calling thread; the others block on the gate.
        let lead = {
            let _gate = Gate(shared);
            run_racer(p, &racer(p, config, 0), 0, shared, winner)
        };
        let mut results = vec![Some(lead)];
        results.extend(
            handles
                .into_iter()
                .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e))),
        );
        results
    });
    let best = shared
        .best
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner);
    merge(p, results, winner.into_inner(), best)
}

/// Opens the ramp-up gate when racer 0 returns, and stops the race if it
/// panicked (the scope would otherwise wait for full searches to finish
/// before propagating the panic).
struct Gate<'a>(&'a Shared);

impl Drop for Gate<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.stop.store(true, Ordering::Relaxed);
        }
        self.0.release();
    }
}

/// Runs racer `k`, maps its point back to the shared column order and
/// claims the win if it proves a result.
fn run_racer(
    p: &MipProblem,
    r: &Racer,
    k: usize,
    shared: &Shared,
    winner: &AtomicUsize,
) -> Result<MipSolution, MipError> {
    let (q, perm) = match &r.permuted {
        Some((q, cols, _)) => (q, Some(cols.as_slice())),
        None => (p, None),
    };
    let link = Link {
        shared,
        perm,
        share: true,
        seen: 0,
        leader: k == 0,
    };
    let mut out = bnb::solve_racing(q, &r.config, link);
    if let (Ok(sol), Some(perm)) = (&mut out, perm)
        && !sol.x.is_empty()
    {
        let mut x = vec![0.0; sol.x.len()];
        for (k, &j) in perm.iter().enumerate() {
            x[j] = sol.x[k];
        }
        sol.x = x;
    }
    if let Ok(sol) = &out
        && proves(sol, shared)
        && winner
            .compare_exchange(usize::MAX, k, Ordering::AcqRel, Ordering::Acquire)
            .is_ok()
    {
        shared.stop.store(true, Ordering::Relaxed);
        shared.release();
    }
    out
}

/// Whether a racer's finished solve answers for the race. An infeasibility
/// claim while the race holds a verified feasible point is wrong (every
/// racer searches the same problem) and does not stop the others.
fn proves(sol: &MipSolution, shared: &Shared) -> bool {
    match sol.status {
        MipStatus::Optimal | MipStatus::Unbounded => true,
        MipStatus::Infeasible => shared.best().is_none(),
        MipStatus::NodeLimit | MipStatus::TimeLimit => false,
    }
}

/// The race's answer: the winner's status, or the limit status with the
/// strongest racer bound, plus the race's best incumbent and summed work.
fn merge(
    p: &MipProblem,
    results: Vec<Option<Result<MipSolution, MipError>>>,
    winner: usize,
    best: Option<(f64, Vec<f64>)>,
) -> Result<MipSolution, MipError> {
    // Maps caller-sense values to internal ones (and back: it is its own
    // inverse).
    let internal = |v: f64| p.caller_objective(v);
    // Racers that never left the ramp-up gate did no work; the winner's
    // index shifts past them.
    let winner = match winner {
        usize::MAX => usize::MAX,
        w => results[..w].iter().filter(|r| r.is_some()).count(),
    };
    let results: Vec<Result<MipSolution, MipError>> = results.into_iter().flatten().collect();
    let k = if winner != usize::MAX {
        winner
    } else {
        match results.iter().position(Result::is_ok) {
            Some(k) => k,
            None => return results.into_iter().next().expect("threads ≥ 2"),
        }
    };
    let mut nodes = 0;
    let mut iters = 0;
    let mut lp_stats = crate::engine::LpStats::default();
    let mut first = None::<bnb::FirstIncumbent>;
    let mut bound = f64::NEG_INFINITY;
    for r in results.iter().flatten() {
        nodes += r.nodes;
        iters += r.lp_iterations;
        lp_stats.absorb(&r.lp_stats);
        if let Some(f) = r.first_incumbent
            && f.source != Source::Racer
            && first.is_none_or(|g| f.seconds < g.seconds)
        {
            first = Some(f);
        }
        bound = bound.max(internal(r.bound));
    }
    // Unbounded and Infeasible answers stand as the winner reported them.
    let mut sol = results.into_iter().nth(k).expect("k < threads")?;
    sol.nodes = nodes;
    sol.lp_iterations = iters;
    sol.lp_stats = lp_stats;
    sol.first_incumbent = first;
    if matches!(
        sol.status,
        MipStatus::Optimal | MipStatus::NodeLimit | MipStatus::TimeLimit
    ) && let Some((z, x)) = best
        && (sol.x.is_empty() || z < internal(sol.objective))
    {
        sol.objective = p.caller_objective(z);
        sol.x = x;
    }
    match sol.status {
        MipStatus::Optimal => sol.bound = sol.objective,
        MipStatus::NodeLimit | MipStatus::TimeLimit => {
            if !sol.x.is_empty() {
                bound = bound.min(internal(sol.objective));
            }
            sol.bound = p.caller_objective(bound);
        }
        MipStatus::Infeasible | MipStatus::Unbounded => {}
    }
    Ok(sol)
}
