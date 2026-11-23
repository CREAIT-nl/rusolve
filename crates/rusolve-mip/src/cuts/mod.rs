//! Root cutting-plane loop.
//!
//! Separators need the LP point x* and the rows; Gomory cuts (tableau.rs) also
//! take the tableau-row weights from [`LpEngine::tableau_weights`]. Every cut is
//! globally valid (derived from the root domain), lives in the LP for the
//! whole tree, and is checked in tests against enumerated integer points.
//!
//! Round structure, stall rule and round limit follow HiGHS
//! `HighsMipSolverData::evaluateRootNode` (reference/HiGHS/highs/mip/
//! HighsMipSolverData.cpp:2157–2280: progress of the LP point along the
//! smoothed cut direction, α = 1/3, stall after 3 rounds without 1% progress
//! and 0.1% objective movement; `maxSepaRounds = 2·√maxTreeSizeLog2`
//! :1990–1993, :1045) and `HighsSeparation::separationRound`
//! (HighsSeparation.cpp:43–180: implied bounds, then the row separators,
//! then pool selection). LP-row aging follows
//! `HighsLpRelaxation::performAging` (HighsLpRelaxation.cpp:629–677, age
//! limit `mip_lp_age_limit = 10`). HiGHS is MIT licensed; it was read as a
//! reference text and reimplemented.
//!
//! On top of the HiGHS stall rule the loop stops after [`CutConfig::weak_rounds`]
//! rounds in a row that each move the bound by less than a fraction of the
//! bound gain so far: on set-cover-like families cuts close little of the
//! gap after 2–5 rounds, and SCIP runs 1.6–2.1× faster there with cuts off. The rule is relative
//! to the gain, as HiGHS's objective test is (`(obj − first) ≤ (prev − first)
//! · 1.001`), not to `|obj|`: a `|obj|`-relative test stopped `gap` after 4
//! rounds at 812.4 against HiGHS's 830.6, because the objective (≈ 810)
//! dwarfs the root gap (≈ 40). A loop that gains nothing drops its cuts.
//! Before the bound first moves the budget is [`CutConfig::flat_weak_rounds`]
//! (degenerate roots can gain late).

mod clique;
mod generate;
mod implied;
mod modk;
mod path;
mod pool;
mod relax;
mod tableau;
#[cfg(test)]
mod tests;
mod transform;

use crate::bnb::MipConfig;
use crate::engine::{CutRow, LpEngine, LpError, LpOutcome, LpStatus};
use crate::problem::MipProblem;
pub(crate) use clique::{CliqueTable, ObjCliques};
use generate::{CutGen, Method};
use implied::VarBounds;
use pool::CutPool;
use relax::{Aggregator, Relaxation};
use std::time::Instant;
use tableau::TableauState;
use transform::TransformedLp;

/// Coefficients below this are zero (HiGHS small_matrix_value).
const EPSILON: f64 = 1e-9;

/// The root fixings clique table over the model rows of `p`. HiGHS
/// builds its table in presolve (`HPresolve.cpp:1494–1512`); ours is built
/// per search round from the presolved problem, which is the same model the
/// round's nodes branch on.
pub(crate) fn root_clique_table(
    p: &MipProblem,
    lower: &[f64],
    upper: &[f64],
    feastol: f64,
) -> CliqueTable {
    let relax = Relaxation::new(p, lower, upper, EPSILON, feastol);
    CliqueTable::extract(&relax, feastol, EPSILON)
}

/// Extend the root fixings clique table with the cliques implied by the
/// root loop's cut rows (HiGHS `extractCliquesFromCut`,
/// HighsCliqueTable.cpp:1078–1266). Cut rows are valid inequalities, so the
/// same set-packing/prefix extraction over them yields valid cliques — on
/// root-dive instances the model rows carry no conflict structure and the
/// separated knapsack rows are where it lives (pg: 2 473 tight cuts).
pub(crate) fn extend_clique_table_with_cuts(
    table: &mut CliqueTable,
    sep: &Separator,
    feastol: f64,
) {
    let relax = &sep.relax;
    table.add_rows(&relax.rows, relax, feastol, EPSILON);
}

/// Which cut families run and how long the root loop may go.
#[derive(Debug, Clone, Copy)]
pub struct CutConfig {
    pub enabled: bool,
    /// `x ≤ c·y + d` from single-row variable bounds.
    pub implied_bounds: bool,
    /// Clique cuts from the conflict graph of the rows.
    pub cliques: bool,
    /// Lift each separated clique with the zero-weight literals adjacent
    /// to all its members (random greedy, HiGHS `ADD_ZERO_WEIGHT_VARS`,
    /// HighsCliqueTable.cpp:25, :1621–1682). Default on: irp root gap
    /// closed 0.39 → 0.79, air05 1.3–1.5× faster, no solve lost and
    /// sgm 0.99× / 1.00× on our benchmark set; nw04's cutoff run is 2.3×
    /// slower.
    pub clique_extend: bool,
    /// Row aggregation along continuous-column paths feeding lifted covers
    /// and c-MIR (also single rows: knapsack covers).
    pub aggregation: bool,
    /// Mod-k (zero-half for k = 2) cuts from combinations of tight rows.
    /// Off by default: it closes the root gap like HiGHS's does (set_cover
    /// 0 → 78%, set_packing 56 → 77%, max_indep_set 70 → 86%,
    /// gen_int_knapsack 8 → 26%). With cold cut-round re-solves,
    /// our medium benchmark tier measured 1.08x time overall (set_cover
    /// 3.56x, set_packing 1.42x) for 0.96x nodes. With warm cut rounds
    /// it is still 1.05x time for 1.00x nodes
    /// (seeds 1–3): set_cover 1.22x time / 0.61x nodes, max_indep_set
    /// 1.67x / 0.61x, gen_int_knapsack 1.28x / 0.81x, set_packing
    /// 1.46x / 1.18x, set_partition 1.36x / 1.53x.
    // TODO(perf): the re-solve is no longer the cost; the rounds are (more
    // separation, more rows in every node LP). Revisit with per-round
    // efficacy filtering or a mod-k round budget, against these numbers.
    pub modk: bool,
    /// Gomory mixed-integer cuts from tableau rows of the fractional basic
    /// integers, selected by the HiGHS tableau
    /// separator's fractionality/density score under an LP-iteration budget
    /// (tableau.rs) and gated by the pool's efficacy/parallelism filter.
    /// On by default: int_production flips from unsolved-in-60s
    /// to solved with Gomory.
    pub gomory: bool,
    /// Gomory cuts longer than `max(gomory_density · average model row length,
    /// 0.1 · numCol)` are discarded (tableau.rs `max_cut_len`); 0 = no cap.
    /// Our rule, not HiGHS's: its tableau separator accepts any cut that
    /// survives cut-generation cancellation (`100 + 0.15·numCol`,
    /// generate.rs), and a sweep of this cap on the medium/mixed tiers found
    /// crossing thresholds — gen_int_knapsack[eq] regresses unless the cap
    /// is ≤ ~0.6·avg row length while plain gen_int_knapsack wants its
    /// near-full-support winners kept — so 0.6 is the best single value.
    pub gomory_density: f64,
    /// Round limit; `None` = HiGHS `2·√maxTreeSizeLog2`.
    pub max_rounds: Option<u32>,
    /// A round improving the root bound by less than this fraction of the
    /// bound gain since the first root LP counts as weak (a round with no
    /// gain yet is weak unless it moves the bound at all).
    pub min_rel_improvement: f64,
    /// This many weak rounds in a row stop the root loop. 4 lets a few
    /// flat rounds pass before a jump (knapsack_cover_s1 closes at round 4
    /// after three flat ones; graph_coloring jumps at round 3). Measured on
    /// the medium tier, 2 vs 3: time 0.87x vs 0.89x, nodes 0.69x vs
    /// 0.78x, both 50/51 solved (old |obj| rule: 48/51).
    pub weak_rounds: u32,
    /// The weak-round budget while the loop has not moved the bound yet
    /// (every round so far flat). HiGHS has no such rule — only its stall
    /// test, which keeps going while the LP point moves — and on degenerate
    /// roots the first gain can arrive late: neos-911970's bound sits at
    /// the LP value for 4 rounds, then rises 23.26 → 47.26 over rounds 5–11.
    /// `weak_rounds` applies once the bound
    /// has moved. 8 (vs the old shared 4): neos-911970 first-root gap closed
    /// 0.00 → 0.76; synthetic medium/mixed node counts identical, 59-instance
    /// benchmark set 10/59 solved either way.
    pub flat_weak_rounds: u32,
    /// Remove every cut when the root loop did not move the bound: such
    /// cuts only make each node LP dearer (bin_packing: 16+ flat rounds,
    /// 17x slower than stopping early).
    pub drop_if_no_gain: bool,
    /// Drop cuts with slack at the final root LP before branching (keeps x*
    /// optimal; makes every node LP smaller).
    pub drop_slack_at_end: bool,
    /// Efficacy variant of cut-row removal, porting HiGHS
    /// `HighsLpRelaxation::removeObsoleteRows` (HighsLpRelaxation.cpp:536,
    /// MIT; read as a reference text): drop LP cut rows whose slack is basic
    /// in the final basis (HiGHS's exact test; the zero-dual tail of
    /// [`obsolete_row_mask`] is the fallback without basis status). This
    /// flag applies it after every root round's re-solve, where HiGHS calls
    /// it after every `rootSeparationRound` (HighsMipSolverData.cpp:2098);
    /// measured negative on lights_out (round-drop churns the pool: s1
    /// 34,563 nodes vs 7,490), stays probe-only.
    pub efficacy_rounds: bool,
    /// And before branching (HighsMipSolverData.cpp:2481): the final drop
    /// uses the basis-status test instead of the activity-slack test of
    /// [`Self::drop_slack_at_end`] — additionally dropping tight-but-basic
    /// (degenerate) rows the slack test keeps. Measured: lights_out
    /// sgm nodes 0.60× (s1 7,490 → 2,905) but gap 1.22× nodes,
    /// int_production 3→2 solved, fcn[fixed] 1.42× nodes — fails the family
    /// gates, so default off.
    pub efficacy_end: bool,
    /// Ours, not HiGHS's: park the Gomory separator for the rest of the
    /// root loop after a round where it ran and none of its cuts were
    /// selected. ON by default: HiGHS has no such rule and removing it
    /// wins lights_out s2 (7,079 → 4,185) and int_production nodes
    /// (0.86×), but fixed_charge_network[fixed] regresses 1.64× time /
    /// 1.95× nodes — the park stays load-bearing.
    /// `--no-gomory-park` turns it off for probing.
    pub gomory_park: bool,
    /// TEMP(cut-rejection probe): drop every cut row before branching.
    pub drop_all_at_end: bool,
    /// TEMP(cut-rejection probe): drop only Gomory-class rows before branching.
    pub drop_gomory_at_end: bool,
    /// Keep the search's incumbent cutoff inside the root cut loop: add
    /// `cᵀx ≤ upper_limit` as a permanent LP row (never aged, never
    /// efficacy-dropped) whenever the cutoff is finite at loop entry, and
    /// round every round's LP point (the search's rounding ladder, see
    /// `root_separation_logged`'s `heur` hook) while no incumbent exists,
    /// adding the row the moment a point verifies. The objective inequality
    /// cutting off the incumbent's value is standard (Wolsey, *Integer
    /// Programming*, 1998, ch. 8); incumbent-interleaved root rounds follow
    /// HiGHS `evaluateRootNode` (HighsMipSolverData.cpp:2109–2119, 2392–2400,
    /// MIT). Measured: int_production
    /// s1 root bound 40,544 → 42,952 (71 % → 85 % of the LP→incumbent gap)
    /// with a sub-MIP-strength cutoff; the row re-energizes the tableau
    /// separator, so rounds keep producing after the pool saturates.
    pub root_cutoff_row: bool,
    /// Keep the rounding hook running after the cutoff row has landed, so a
    /// better verified point tightens the row in place (the bounded form of
    /// HiGHS's incumbent-driven root restart). Off by default: measured on
    /// int_production s1, the tightenings steer the tree worse (first
    /// incumbent 90,865 → 78,029 moves the tree 90,499 → 152,365 nodes and
    /// past 20 s, where the single landing at the first verified point wins).
    /// The machinery stays live for
    /// a family whose rounding ladder actually reaches a binding incumbent.
    pub root_cutoff_tighten: bool,
    /// Separate (globally valid) cuts at tree nodes taken from the node
    /// queue, like HiGHS `mip_allow_cut_separation_at_nodes`. Off by default,
    /// re-measured with warm node LPs + aging (seeds 1–3 vs `bd7ca03`):
    /// the aging below takes the tax from medium 1.26x / mixed 1.58x to
    /// 1.09x / 1.53x (limit 10; limit 20 is worse: 1.24x / 1.64x), still
    /// failing the ≤ 1.00x gate. Nodes do fall (gap 0.16x, int_production
    /// 0.35x, lights_out 0.53x) but the families that gain nodes lose more
    /// time to rows (lights_out 2.80x, fcn[fixed] 2.42x, elastic_cover
    /// 2.29x, market_split 1.84x, multi_knapsack 1.46x). The remaining lever
    /// is cut *quality*, not budget: efficacy-end (the HiGHS rule that keeps 0 cuts
    /// in lights_out's LP) and a per-round selection cap — 4,410 tree cuts
    /// on lights_out_s1 alone. The knob stays reachable as `--tree-cuts`.
    pub tree: bool,
    /// Node separation may spend at most this fraction of all LP iterations
    /// (plus [`TREE_EFFORT_OFFSET`]); see [`tree_separation_allowed`].
    pub tree_effort: f64,
    /// HiGHS `mip_lp_age_limit` (HighsOptions.h:1176–1180): an LP cut row
    /// slack for this many aging epochs (node LPs once tree cuts are on,
    /// see [`Separator::age_node`]) leaves the LP
    /// (HighsLpRelaxation::performAging, HighsLpRelaxation.cpp:629–677).
    pub lp_age_limit: u32,
    /// Probe every fractional binary of a root round's LP point once
    /// (fix to 1, propagate; fix to 0, propagate) and feed the binary
    /// implications to the separator's clique table as 2-cliques before
    /// clique separation; an infeasible direction becomes a bound cut.
    /// This is the probing half of HiGHS `separateImpliedBounds`
    /// (HighsImplications.cpp:583–641, `runProbing` :298–389,
    /// `computeImplications` :15–140; called first in every
    /// `separationRound`, HighsSeparation.cpp:85–90), MIT, read as a
    /// reference text. Off: in the audit A/B it loses enlight_hard and costs
    /// 2.5–3.4 s of probing on neos-860300 for +0.05 gap closed.
    pub probe_cliques: bool,
    /// HiGHS's two-speed LP-row aging in the tree. Every optimal node
    /// LP only *ticks* (a slack row that has aged before ages one more; a
    /// tight row resets to 0; nothing is removed — HiGHS
    /// `performAging(false)` after `evaluateNode`, HighsSearch.cpp:958,
    /// HighsLpRelaxation.cpp:640–676 `age += (deleteRows || age != 0)`).
    /// Rows leave only at the aging epochs of node separation calls (a
    /// separation round's re-solve, or a skipped separation:
    /// HighsSeparation.cpp:166, :221–226, `performAging(true)`). Off = every node LP is a deleting
    /// epoch. Measured on the tree-cut config: the per-node-LP deletions
    /// strip rows that stored node snapshots still hold nonbasic, and those
    /// snapshots then warm-start from an unrelated basis.
    pub tree_age_ticks: bool,
}

impl Default for CutConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            implied_bounds: true,
            cliques: true,
            clique_extend: true,
            aggregation: true,
            modk: false,
            gomory: true,
            gomory_density: 0.6,
            max_rounds: None,
            min_rel_improvement: 1e-3,
            weak_rounds: 4,
            flat_weak_rounds: 8,
            drop_if_no_gain: true,
            drop_slack_at_end: true,
            efficacy_rounds: false,
            efficacy_end: false,
            gomory_park: true,
            drop_all_at_end: false,
            drop_gomory_at_end: false,
            root_cutoff_row: true,
            root_cutoff_tighten: false,
            tree: false,
            tree_effort: 0.1,
            lp_age_limit: 10,
            probe_cliques: false,
            tree_age_ticks: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CutClass {
    ImpliedBound,
    Clique,
    KnapsackCover,
    MixedBinaryCover,
    MixedIntegerCover,
    Cmir,
    ModK,
    Gomory,
}

impl CutClass {
    pub const ALL: [CutClass; 8] = [
        CutClass::ImpliedBound,
        CutClass::Clique,
        CutClass::KnapsackCover,
        CutClass::MixedBinaryCover,
        CutClass::MixedIntegerCover,
        CutClass::Cmir,
        CutClass::ModK,
        CutClass::Gomory,
    ];

    pub fn name(self) -> &'static str {
        match self {
            CutClass::ImpliedBound => "implied_bound",
            CutClass::Clique => "clique",
            CutClass::KnapsackCover => "knapsack_cover",
            CutClass::MixedBinaryCover => "mixed_bin_cover",
            CutClass::MixedIntegerCover => "mixed_int_cover",
            CutClass::Cmir => "cmir",
            CutClass::ModK => "modk",
            CutClass::Gomory => "gomory",
        }
    }

    fn from_method(m: Method) -> Self {
        match m {
            Method::KnapsackCover => CutClass::KnapsackCover,
            Method::MixedBinaryCover => CutClass::MixedBinaryCover,
            Method::MixedIntegerCover => CutClass::MixedIntegerCover,
            Method::Cmir => CutClass::Cmir,
        }
    }
}

/// What the cut loops did (root unless prefixed `tree_`), per [`CutClass`] (indexed by `class as usize`).
#[derive(Debug, Clone, Default)]
pub struct CutStats {
    pub rounds: u32,
    /// Caller-sense root LP bound before / after cuts (NaN if no loop ran).
    pub root_bound_before: f64,
    pub root_bound_after: f64,
    pub generated: [u64; CutClass::ALL.len()],
    pub applied: [u64; CutClass::ALL.len()],
    /// Cuts still in the LP when branching starts.
    pub in_lp: [u64; CutClass::ALL.len()],
    /// Generated cuts whose base row aggregated two or more LP rows.
    pub aggregated: u64,
    /// Node separation: rounds run and cuts added to the LP in the tree.
    pub tree_rounds: u64,
    pub tree_applied: u64,
    /// Root rounds in which the Gomory separator actually ran (not parked
    /// by the round gate; see `TableauState::parked`).
    pub gomory_calls: u64,
    /// The incumbent cutoff row `cᵀx ≤ upper_limit` was added to the root LP
    /// (`CutConfig::root_cutoff_row`): 1 when added at loop entry from the
    /// search's cutoff, 2 also when a mid-loop rounding point added it.
    pub cutoff_rows: u32,
    /// Feasible points the mid-loop rounding hook returned (verified and
    /// used for the row, and handed back for `try_incumbent`).
    pub root_incumbents: Vec<(Vec<f64>, crate::heuristics::Source)>,
    pub tree_lp_iterations: u64,
    /// In-loop probing (`CutConfig::probe_cliques`): binaries probed,
    /// 2-cliques stored, infeasible directions (bound cuts), seconds.
    pub probe_probed: u64,
    pub probe_cliques: u64,
    pub probe_fixings: u64,
    /// Implied-bound cuts from cached probing implications entering the pool.
    pub probe_ibcuts: u64,
    pub probe_s: f64,
    /// Node separation calls (queue nodes within the effort cap).
    pub tree_calls: u64,
    /// Cut-loop re-solves (root or tree) that failed; the loop stopped
    /// on the last good LP instead of failing the solve.
    pub lp_failures: u64,
    pub lp_iterations: u64,
    pub time_s: f64,
    /// Per-round cost split (root loop here; node rounds in the
    /// `tree_` copies): separator wall time, pool selection, the re-solve
    /// after adding cuts, LP-row aging (with the row removal), and the LP
    /// rows that aging removed.
    pub sep_nanos: u64,
    pub pool_nanos: u64,
    pub resolve_nanos: u64,
    pub age_nanos: u64,
    pub aged_rows: u64,
    pub tree_sep_nanos: u64,
    pub tree_pool_nanos: u64,
    pub tree_resolve_nanos: u64,
    /// Node-LP aging ticks ([`Separator::age_node`]) count here too.
    pub tree_age_nanos: u64,
    pub tree_aged_rows: u64,
}

/// Root clique-table fixings counters.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CliqueFixStats {
    /// Cliques extracted from the objective cutoff row over the solve.
    pub obj_cliques: u64,
    /// Binary bound tightenings the clique propagation produced (post-root
    /// incumbents are excluded by design; see `bnb::Search::root_phase`).
    pub fixings: u64,
    /// An extraction proved the box's minimum objective above the cutoff
    /// (the incumbent is optimal).
    pub cutoff_pruned: bool,
}

/// The root cut loop's rounding hook ([`root_separation`]): proposes
/// candidate points for one round's LP point (with their [`Source`] and the
/// LP iterations its repair step spent), under the loop's current cutoff
/// (`f64::INFINITY` before the first row lands).
pub(crate) type RootHeur<'a, E> =
    dyn FnMut(&[f64], &mut E, f64) -> (Vec<(Vec<f64>, crate::heuristics::Source)>, u64) + 'a;

/// An LP cut row in engine order: its pool index and age, `pool: None` for
/// the incumbent cutoff row `cᵀx ≤ upper_limit`, which no pool bookkeeping,
/// aging or efficacy rule ever touches (it is a valid inequality for every
/// solution the search still wants, not a separator cut).
#[derive(Clone, Copy)]
struct LpCut {
    pool: Option<usize>,
    age: u32,
}

/// Cut-loop state that outlives the root: the pool, the separators' tables
/// and the LP-cut bookkeeping, all over the global (root) domain so every cut
/// stays valid for the whole tree. The tree reuses it for node separation
/// (HighsSeparation::separate, HighsSeparation.cpp:181–228).
pub(crate) struct Separator {
    relax: Relaxation,
    vb: VarBounds,
    cliques: Option<CliqueTable>,
    pool: CutPool,
    cutgen: CutGen,
    /// The tableau separator's own generator (isolated RNG stream).
    tableau_cutgen: CutGen,
    agg: Aggregator,
    tableau: TableauState,
    rng: u64,
    /// Clique-extension shuffles (own stream: toggling the extension must
    /// not perturb the path separator's draws).
    clique_rng: u64,
    /// Pool index and age of each LP cut, in engine order.
    lp_cuts: Vec<LpCut>,
    /// Whether the incumbent cutoff row `cᵀx ≤ upper_limit` is in the LP.
    cutoff_row: bool,
    /// Objective value of the newest cutoff row in the LP (`f64::INFINITY`
    /// before the first one); a new row only lands when it beats this, so
    /// the row tightens with the incumbent like HiGHS's `upper_bound_`.
    cutoff_ub: f64,
    epochs: u32,
    /// `CutConfig::lp_age_limit` (HiGHS `mip_lp_age_limit`).
    lp_age_limit: u32,
    /// Final root LP objective (HiGHS `rootlpsolobj`).
    root_obj: f64,
    /// In-loop probing state (`CutConfig::probe_cliques`).
    probe: Option<ProbeState>,
}

/// In-loop probing (HiGHS `HighsImplications::separateImpliedBounds`'s
/// probing pass): which binaries were probed, the row propagator the probes
/// run on, and scratch boxes.
struct ProbeState {
    prop: crate::propagate::Propagator,
    probed: Vec<bool>,
    lo_s: Vec<f64>,
    hi_s: Vec<f64>,
    /// Cached implications per literal `2·j + val`: `(col, is_upper, bound)`
    /// (HiGHS `implications[2*col+val]`, HighsImplications.cpp:15–140);
    /// empty when not probed or the direction is infeasible.
    implics: Vec<Vec<(usize, bool, f64)>>,
}

impl Separator {
    fn new(p: &MipProblem, lower: &[f64], upper: &[f64], cc: &CutConfig, feastol: f64) -> Self {
        let relax = Relaxation::new(p, lower, upper, EPSILON, feastol);
        let vb = if cc.implied_bounds {
            VarBounds::extract(&relax, feastol, EPSILON)
        } else {
            VarBounds::empty(p.num_col)
        };
        let probing = cc.cliques && cc.probe_cliques;
        let cliques = if cc.cliques {
            // With in-loop probing an empty table still fills up (2-cliques
            // from the probes), so it is kept.
            Some(CliqueTable::extract(&relax, feastol, EPSILON))
                .filter(|t| probing || !t.is_empty())
        } else {
            None
        };
        let probe = probing.then(|| ProbeState {
            prop: crate::propagate::Propagator::new(p, feastol),
            probed: vec![false; p.num_col],
            lo_s: vec![0.0; p.num_col],
            hi_s: vec![0.0; p.num_col],
            implics: vec![Vec::new(); 2 * p.num_col],
        });
        let agg = Aggregator::new(relax.num_ext());
        Self {
            relax,
            vb,
            cliques,
            pool: CutPool::new(),
            cutgen: CutGen::new(p.num_col, feastol, EPSILON, 0x5eed_c0de),
            // The tableau separator generates from its own stream: sharing
            // one CutGen makes every Gomory cover-heuristic call advance the
            // RNG that the path/mod-k covers of later rounds draw from, so
            // turning Gomory on perturbs the other separators' cuts.
            tableau_cutgen: CutGen::new(p.num_col, feastol, EPSILON, 0x05ee_d901),
            agg,
            tableau: TableauState::default(),
            rng: 0x2545_f491_4f6c_dd1d,
            clique_rng: 0x9e37_79b9_7f4a_7c15,
            lp_cuts: Vec::new(),
            cutoff_row: false,
            cutoff_ub: f64::INFINITY,
            epochs: 0,
            lp_age_limit: cc.lp_age_limit,
            root_obj: f64::NEG_INFINITY,
            probe,
        }
    }

    /// Probe the fractional binaries of `x` not probed yet, in both
    /// directions, over the root box (HiGHS `separateImpliedBounds`,
    /// HighsImplications.cpp:593–615: every fractional binary without
    /// cached implications is probed once, before the clique separator
    /// runs; `computeImplications` :97–112 turns each binary implication
    /// `j = v ⇒ k = b` into the 2-clique `{j = v, k = 1 − b}`). An
    /// infeasible direction proves the opposite value; it enters the pool
    /// as the bound cut HiGHS adds there (:656–661, :723–728). Cliques are
    /// derived without the incumbent cutoff, so they hold for every
    /// integer-feasible point, like the rest of the table.
    fn probe_round(
        &mut self,
        feastol: f64,
        stats: &mut CutStats,
        log: &mut Option<&mut Vec<(CutClass, CutRow)>>,
    ) {
        let (Some(ps), Some(table)) = (self.probe.as_mut(), self.cliques.as_mut()) else {
            return;
        };
        let start = Instant::now();
        let relax = &self.relax;
        for j in 0..relax.num_col {
            if ps.probed[j] || !relax.is_binary(j) {
                continue;
            }
            let v = relax.x[j];
            if v <= feastol || v >= 1.0 - feastol {
                continue;
            }
            ps.probed[j] = true;
            stats.probe_probed += 1;
            let mut fix: Option<(f64, f64)> = None; // (coef, rhs) of a bound cut
            for val in [true, false] {
                let dir = crate::probing::probe(
                    j,
                    val,
                    &mut ps.prop,
                    &relax.lower,
                    &relax.upper,
                    &mut ps.lo_s,
                    &mut ps.hi_s,
                    f64::INFINITY,
                );
                if dir.infeasible {
                    // j = 1 infeasible ⇒ x_j ≤ 0; j = 0 infeasible ⇒ −x_j ≤ −1.
                    fix = Some(if val { (1.0, 0.0) } else { (-1.0, -1.0) });
                    continue;
                }
                for &(k, b) in &dir.fixed_binary {
                    if table.add_probing_clique((j, val), (k, b <= 0.5)) {
                        stats.probe_cliques += 1;
                    }
                }
                let lit = &mut ps.implics[2 * j + val as usize];
                lit.extend(dir.upper.iter().map(|&(k, b)| (k, true, b)));
                lit.extend(dir.lower.iter().map(|&(k, b)| (k, false, b)));
            }
            if let Some((a, rhs)) = fix
                && self
                    .pool
                    .add(&[j], &[a], rhs, true, CutClass::ImpliedBound)
                    .is_some()
            {
                stats.probe_fixings += 1;
                stats.generated[CutClass::ImpliedBound as usize] += 1;
                if let Some(log) = log.as_deref_mut() {
                    log.push((
                        CutClass::ImpliedBound,
                        CutRow {
                            inds: vec![j],
                            vals: vec![a],
                            rhs,
                        },
                    ));
                }
            }
        }
        // Implied-bound cuts from the cached implications (HiGHS
        // separateImpliedBounds second pass, HighsImplications.cpp:~650–760,
        // MIT): `x_j = 1 ⇒ x_k ≤ b` gives `x_k + (u_k − b)·x_j ≤ u_k`;
        // `x_j = 1 ⇒ x_k ≥ b` gives `−x_k + (b − l_k)·x_j ≤ −l_k`; the
        // `x_j = 0` side substitutes `1 − x_j`.
        for j in 0..relax.num_col {
            if !ps.probed[j] {
                continue;
            }
            let v = relax.x[j];
            if v <= feastol || v >= 1.0 - feastol {
                continue;
            }
            for val in [true, false] {
                for &(k, is_upper, b) in &ps.implics[2 * j + val as usize] {
                    let (lk, uk) = (relax.lower[k], relax.upper[k]);
                    let inds = [k, j];
                    let mut vals = [0.0; 2];
                    let mut rhs;
                    if is_upper {
                        if b + feastol >= uk {
                            continue;
                        }
                        vals[0] = 1.0;
                        vals[1] = uk - b;
                        rhs = uk;
                    } else {
                        if b - feastol <= lk {
                            continue;
                        }
                        vals[0] = -1.0;
                        vals[1] = b - lk;
                        rhs = -lk;
                    }
                    if !val {
                        // x_j = 0 side: coefficient on (1 − x_j).
                        rhs -= vals[1];
                        vals[1] = -vals[1];
                    }
                    let viol = vals[0] * relax.x[k] + vals[1] * v - rhs;
                    if viol <= feastol {
                        continue;
                    }
                    if inds[0] == inds[1] {
                        continue;
                    }
                    let integral = relax.is_integral(k);
                    if self
                        .pool
                        .add(&inds, &vals, rhs, integral, CutClass::ImpliedBound)
                        .is_some()
                    {
                        stats.probe_ibcuts += 1;
                        stats.generated[CutClass::ImpliedBound as usize] += 1;
                        if let Some(log) = log.as_deref_mut() {
                            log.push((
                                CutClass::ImpliedBound,
                                CutRow {
                                    inds: inds.to_vec(),
                                    vals: vals.to_vec(),
                                    rhs,
                                },
                            ));
                        }
                    }
                }
            }
        }
        stats.probe_s += start.elapsed().as_secs_f64();
    }

    /// One separation round at the optimal point of `lp` (solved under the
    /// node bounds `lower`/`upper`): implied bounds, cliques, path
    /// aggregation, mod-k, then pool selection (HighsSeparation::separationRound,
    /// HighsSeparation.cpp:43–180). Returns `None` when no cut was selected,
    /// otherwise the re-solved LP; when the re-solve fails the round's rows
    /// are removed again and the error is returned for the caller to stop
    /// on `lp`.
    #[allow(clippy::too_many_arguments)]
    fn round<E: LpEngine>(
        &mut self,
        engine: &mut E,
        lower: &[f64],
        upper: &[f64],
        lp: &LpOutcome,
        cc: &CutConfig,
        feastol: f64,
        stats: &mut CutStats,
        mut log: Option<&mut Vec<(CutClass, CutRow)>>,
    ) -> Result<Option<LpOutcome>, LpError> {
        self.relax.set_solution(&lp.x);
        self.probe_round(feastol, stats, &mut log);
        let sep_start = Instant::now();
        let relax = &mut self.relax;
        let pool = &mut self.pool;

        if cc.implied_bounds {
            for (inds, vals, rhs) in self.vb.separate(relax, feastol) {
                let integral = row_integral(relax, &inds, &vals);
                if pool
                    .add(&inds, &vals, rhs, integral, CutClass::ImpliedBound)
                    .is_some()
                {
                    stats.generated[CutClass::ImpliedBound as usize] += 1;
                    if let Some(log) = log.as_deref_mut() {
                        log.push((CutClass::ImpliedBound, CutRow { inds, vals, rhs }));
                    }
                }
            }
        }

        if let Some(table) = self.cliques.as_mut() {
            for (inds, vals, rhs) in
                table.separate(&lp.x, feastol, cc.clique_extend, &mut self.clique_rng)
            {
                if pool
                    .add(&inds, &vals, rhs, true, CutClass::Clique)
                    .is_some()
                {
                    stats.generated[CutClass::Clique as usize] += 1;
                    if let Some(log) = log.as_deref_mut() {
                        log.push((CutClass::Clique, CutRow { inds, vals, rhs }));
                    }
                }
            }
        }

        let mut trans = (cc.aggregation || cc.modk || cc.gomory)
            .then(|| TransformedLp::new(relax, &self.vb, feastol));
        if let (true, Some(trans)) = (cc.aggregation, trans.as_mut()) {
            path::separate(
                relax,
                trans,
                &mut self.agg,
                &mut self.cutgen,
                &mut self.rng,
                feastol,
                |cut, path_len| {
                    let class = CutClass::from_method(cut.method);
                    if pool
                        .add(&cut.inds, &cut.vals, cut.rhs, cut.integral, class)
                        .is_none()
                    {
                        return false;
                    }
                    stats.generated[class as usize] += 1;
                    if path_len > 1 {
                        stats.aggregated += 1;
                    }
                    if let Some(log) = log.as_deref_mut() {
                        log.push((
                            class,
                            CutRow {
                                inds: cut.inds,
                                vals: cut.vals,
                                rhs: cut.rhs,
                            },
                        ));
                    }
                    true
                },
            );
        }
        if let (true, Some(trans)) = (cc.modk, trans.as_mut()) {
            modk::separate(
                relax,
                trans,
                &mut self.agg,
                &mut self.cutgen,
                feastol,
                EPSILON,
                |cut| {
                    let class = CutClass::ModK;
                    if pool
                        .add(&cut.inds, &cut.vals, cut.rhs, cut.integral, class)
                        .is_none()
                    {
                        return false;
                    }
                    stats.generated[class as usize] += 1;
                    if let Some(log) = log.as_deref_mut() {
                        log.push((
                            class,
                            CutRow {
                                inds: cut.inds,
                                vals: cut.vals,
                                rhs: cut.rhs,
                            },
                        ));
                    }
                    true
                },
            );
        }

        let gomory_ran = cc.gomory && !self.tableau.parked;
        if let (true, Some(trans)) = (gomory_ran, trans.as_mut()) {
            tableau::separate(
                engine,
                lower,
                upper,
                lp,
                relax,
                trans,
                &mut self.agg,
                &mut self.tableau_cutgen,
                &mut self.tableau,
                stats.lp_iterations,
                feastol,
                EPSILON,
                cc.gomory_density,
                |cut| {
                    let class = CutClass::Gomory;
                    if pool
                        .add(&cut.inds, &cut.vals, cut.rhs, cut.integral, class)
                        .is_none()
                    {
                        return false;
                    }
                    stats.generated[class as usize] += 1;
                    if let Some(log) = log.as_deref_mut() {
                        log.push((
                            class,
                            CutRow {
                                inds: cut.inds,
                                vals: cut.vals,
                                rhs: cut.rhs,
                            },
                        ));
                    }
                    true
                },
            );
        }
        if gomory_ran {
            stats.gomory_calls += 1;
        }

        let pool_start = Instant::now();
        stats.sep_nanos += (pool_start - sep_start).as_nanos() as u64;
        let selected = pool.separate(&lp.x, lower, upper, feastol);
        stats.pool_nanos += pool_start.elapsed().as_nanos() as u64;
        if std::env::var("RUSOLVE_CUT_DEBUG").is_ok() {
            let mut row = String::new();
            for cls in CutClass::ALL {
                let sel = selected
                    .iter()
                    .filter(|&&c| pool.cuts[c].class == cls)
                    .count();
                if stats.generated[cls as usize] > 0 || sel > 0 {
                    row.push_str(&format!(
                        " {}:{}g+{}s",
                        cls.name(),
                        stats.generated[cls as usize],
                        sel
                    ));
                }
            }
            eprintln!(
                "  sep round {} pool {} avail {}{row}",
                stats.rounds,
                pool.cuts.len(),
                pool.cuts.iter().filter(|c| !c.deleted).count()
            );
        }
        if selected.is_empty() {
            return Ok(None);
        }
        // Round-level Gomory gate (ours, not HiGHS's): a round that ran the
        // tableau separator and had none of its cuts selected parks it for
        // the rest of the root loop (TableauState::parked;
        // generated-but-unselected cuts keep coming back while the LP point
        // barely moves). HiGHS has no such rule and the efficacy-end default
        // interacts badly with it (lights_out s2: efficacy-end 12,001 with
        // the park, 7,599 without), so it is opt-in via CutConfig.
        if cc.gomory_park
            && gomory_ran
            && !selected
                .iter()
                .any(|&c| pool.cuts[c].class == CutClass::Gomory)
        {
            self.tableau.parked = true;
        }
        let rows: Vec<CutRow> = selected
            .iter()
            .map(|&c| {
                let cut = &pool.cuts[c];
                CutRow {
                    inds: cut.inds.clone(),
                    vals: cut.vals.clone(),
                    rhs: cut.rhs,
                }
            })
            .collect();
        let resolve_start = Instant::now();
        engine.add_cuts(&rows);
        for &c in &selected {
            let cut = &pool.cuts[c];
            relax.push_cut(&cut.inds, &cut.vals, cut.rhs, cut.integral);
            stats.applied[cut.class as usize] += 1;
            self.lp_cuts.push(LpCut {
                pool: Some(c),
                age: 0,
            });
        }
        let next = match engine.solve(lower, upper) {
            Ok(next) => next,
            Err(e) => {
                // The re-solve failed. Take this round's rows back out
                // so the engine and relaxation match `lp` again (still
                // optimal for them); the caller stops on it, as HiGHS ends
                // its separation loop on a non-optimal re-solve
                // (HighsSeparation.cpp:78, :210).
                let n = self.lp_cuts.len();
                let first = n - selected.len();
                let keep: Vec<bool> = (0..n).map(|i| i < first).collect();
                engine.retain_cuts(&keep);
                relax.retain_cuts(&keep);
                relax.set_solution(&lp.x);
                for e in self.lp_cuts.drain(first..) {
                    if let Some(c) = e.pool {
                        stats.applied[pool.cuts[c].class as usize] -= 1;
                        pool.lp_cut_removed(c);
                    }
                }
                return Err(e);
            }
        };
        stats.resolve_nanos += resolve_start.elapsed().as_nanos() as u64;
        stats.lp_iterations += next.iterations;
        if next.status == LpStatus::Optimal {
            let age_start = Instant::now();
            self.relax.set_solution(&next.x);
            self.epochs += 1;
            stats.aged_rows += age_lp_cuts(
                &mut self.relax,
                engine,
                &mut self.pool,
                &mut self.lp_cuts,
                self.epochs,
                self.lp_age_limit,
                feastol,
            ) as u64;
            stats.age_nanos += age_start.elapsed().as_nanos() as u64;
        }
        Ok(Some(next))
    }

    /// Separation at a tree node whose LP `lp` is optimal under the node
    /// bounds: rounds continue while the objective keeps moving by more than
    /// 1% of its gain over the root bound (HighsSeparation.cpp:185–218).
    /// Cuts are derived from the global domain, so they enter the LP for the
    /// whole tree like root cuts and age out the same way. Aging of nodes
    /// that run no round is the caller's [`Self::age_node`] tick, not
    /// repeated here.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn separate_node<E: LpEngine>(
        &mut self,
        p: &MipProblem,
        engine: &mut E,
        lower: &[f64],
        upper: &[f64],
        mut lp: LpOutcome,
        config: &MipConfig,
        upper_limit: f64,
        total_lp_iterations: u64,
        stats: &mut CutStats,
    ) -> Result<LpOutcome, LpError> {
        let feastol = config.feastol;
        if lp.status != LpStatus::Optimal {
            return Ok(lp);
        }
        if !tree_separation_allowed(&config.cuts, stats, total_lp_iterations) {
            // HiGHS still ages when it skips separation
            // (HighsSeparation.cpp:221–226); the caller's [`Self::age_node`]
            // tick already counted this node LP, so returning without one
            // here would double-age it — unless that tick does not delete
            // (`tree_age_ticks`), when this is the deleting epoch.
            if config.cuts.tree_age_ticks {
                self.age_epoch(engine, &lp.x, feastol, stats);
            }
            return Ok(lp);
        }
        // Node rounds feed the shared pool but report only into the tree
        // counters: `generated`/`applied` describe the root loop.
        let mut scratch = CutStats::default();
        while lp.status == LpStatus::Optimal
            && lp.objective < upper_limit
            && has_fractional(p, &lp.x, config.integrality_tol)
        {
            let last = lp.objective;
            let before: u64 = scratch.applied.iter().sum();
            // A failed re-solve ends node separation on the node's
            // last good LP (HighsSeparation.cpp:210).
            let Ok(Some(next)) = self
                .round(
                    engine,
                    lower,
                    upper,
                    &lp,
                    &config.cuts,
                    feastol,
                    &mut scratch,
                    None,
                )
                .inspect_err(|_| stats.lp_failures += 1)
            else {
                break;
            };
            stats.tree_rounds += 1;
            stats.tree_applied += scratch.applied.iter().sum::<u64>() - before;
            lp = next;
            if lp.status != LpStatus::Optimal
                || lp.objective - self.root_obj <= (last - self.root_obj).max(feastol) * 1.01
            {
                break;
            }
        }
        stats.tree_lp_iterations += scratch.lp_iterations;
        stats.tree_calls += 1;
        stats.tree_sep_nanos += scratch.sep_nanos;
        stats.tree_pool_nanos += scratch.pool_nanos;
        stats.tree_resolve_nanos += scratch.resolve_nanos;
        stats.tree_age_nanos += scratch.age_nanos;
        stats.tree_aged_rows += scratch.aged_rows;
        Ok(lp)
    }

    /// Ages the LP cuts against one more node-LP point `x`, called for every
    /// optimal node LP (dive children included). HiGHS ticks the row ages
    /// after every optimal node LP without deleting
    /// (`HighsSearch::evaluateNodeSolution`, HighsSearch.cpp:958,
    /// `performAging()` with `deleteRows = false`) and deletes only in the
    /// aging that follows a separation round or a skipped separation
    /// (HighsSeparation.cpp:166, :224, `performAging(true)`). We merge the
    /// two modes into one epoch-with-deletion per node LP: our deletion
    /// keeps the next solve's basis (no immediate re-optimise, unlike
    /// HiGHS's `removeCuts`), and one epoch per node LP is the same age
    /// accounting as HiGHS's tick. Rows slack across k node LPs then leave
    /// even though separation runs only at backtracked queue nodes — aging
    /// per separation call alone lets the rows of a dive's first separated
    /// node ride along for the whole dive.
    pub(crate) fn age_node<E: LpEngine>(
        &mut self,
        engine: &mut E,
        x: &[f64],
        feastol: f64,
        ticks_only: bool,
        stats: &mut CutStats,
    ) {
        if ticks_only {
            // `CutConfig::tree_age_ticks`: HiGHS `performAging(false)`.
            let start = Instant::now();
            self.relax.set_solution(x);
            let m0 = self.relax.num_model_rows;
            for (k, entry) in self.lp_cuts.iter_mut().enumerate() {
                if entry.pool.is_none() {
                    continue;
                }
                if self.relax.activity[m0 + k] < self.relax.rows[m0 + k].upper - feastol {
                    entry.age += u32::from(entry.age != 0);
                } else {
                    entry.age = 0;
                }
            }
            stats.tree_age_nanos += start.elapsed().as_nanos() as u64;
            return;
        }
        self.age_epoch(engine, x, feastol, stats);
    }

    /// One deleting aging epoch at the node LP point `x` (HiGHS
    /// `performAging(true)`).
    fn age_epoch<E: LpEngine>(
        &mut self,
        engine: &mut E,
        x: &[f64],
        feastol: f64,
        stats: &mut CutStats,
    ) {
        let start = Instant::now();
        self.relax.set_solution(x);
        self.epochs += 1;
        stats.tree_aged_rows += age_lp_cuts(
            &mut self.relax,
            engine,
            &mut self.pool,
            &mut self.lp_cuts,
            self.epochs,
            self.lp_age_limit,
            feastol,
        ) as u64;
        stats.tree_age_nanos += start.elapsed().as_nanos() as u64;
    }

    /// Audit: cuts in the LP and live (not deleted) pool cuts.
    pub(crate) fn pool_sizes(&self) -> (usize, usize) {
        (
            self.lp_cuts.len(),
            self.pool.cuts.iter().filter(|c| !c.deleted).count(),
        )
    }
}

/// A re-solve inside the root cut loop failed (iteration limit, the
/// deadline, simplex trouble). The loop keeps its last good LP and stops,
/// as HiGHS breaks out of separation on any non-optimal re-solve
/// (HighsSeparation.cpp:78, :210) rather than failing the solve.
fn lp_failed(stats: &mut CutStats, round: u32, e: &LpError) {
    stats.lp_failures += 1;
    if std::env::var("RUSOLVE_CUT_DEBUG").is_ok() {
        eprintln!(
            "round {round}: cut-loop LP re-solve failed ({e}) — loop ends on the last good LP"
        );
    }
}

/// Initial node-separation allowance in LP iterations.
pub(crate) const TREE_EFFORT_OFFSET: u64 = 10_000;

/// Candidate-point tolerances mirroring `bnb::CANDIDATE_TOL` (a rounding
/// point must violate bounds/rows by at most this and integrality by the
/// configured tolerance before it may set the root cutoff).
const CANDIDATE_TOL: f64 = 1e-6;

/// Verify a proposed rounding point and return its internal objective value
/// (integers snapped, same check as `bnb::Search::try_incumbent`).
fn verified_internal_objective(p: &MipProblem, x: &[f64], int_tol: f64) -> Option<f64> {
    if x.len() != p.num_col || x.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let v = p.violation(x);
    if v.bound > CANDIDATE_TOL || v.row > CANDIDATE_TOL || v.integrality > int_tol {
        return None;
    }
    let mut obj = 0.0;
    for (j, (&xj, &int)) in x.iter().zip(&p.integer).enumerate() {
        obj += p.cost[j] * if int { xj.round() } else { xj };
    }
    obj.is_finite().then_some(obj)
}

/// Add the incumbent cutoff row `cᵀx ≤ ub` to the engine and the relaxation
/// (`ub` internal space, finite). The row is slack at the current LP point
/// (`lp.objective ≤ ub` wherever this is called), so it never needs a
/// re-solve to keep `x*` optimal, and it is tracked as a pool-less
/// [`LpCut`] so the aging and efficacy bookkeeping stays uniform — but it
/// never reaches the tree: `root_separation` strips pool-less rows before
/// branching.
fn add_cutoff_row<E: LpEngine>(
    sep: &mut Separator,
    engine: &mut E,
    p: &MipProblem,
    ub: f64,
) -> bool {
    let strictly_better = ub
        .partial_cmp(&sep.cutoff_ub)
        .is_some_and(|o| o == std::cmp::Ordering::Less);
    if !strictly_better || (sep.cutoff_ub.is_finite() && !moved(ub, sep.cutoff_ub)) {
        return false;
    }
    let (mut inds, mut vals) = (Vec::new(), Vec::new());
    for (j, &c) in p.cost.iter().enumerate() {
        if c != 0.0 {
            inds.push(j);
            vals.push(c);
        }
    }
    if inds.is_empty() {
        return false;
    }
    if !engine.add_cuts(&[CutRow {
        inds: inds.clone(),
        vals: vals.clone(),
        rhs: ub,
    }]) {
        return false;
    }
    sep.relax.push_cut(&inds, &vals, ub, false);
    sep.lp_cuts.push(LpCut { pool: None, age: 0 });
    sep.cutoff_row = true;
    sep.cutoff_ub = ub;
    true
}

/// Effort cap on node separation, shaped like HiGHS
/// `HighsMipSolverData::moreHeuristicsAllowed` (HighsMipSolverData.cpp:
/// 669–700: heuristic LP iterations below `effort · total + 10000`). HiGHS
/// itself separates at every queue node because its re-solves are warm dual
/// simplex; ours spend a full re-solve per round on top of the row tax, so
/// the cap keeps separation to the top of the tree where it pays most.
/// The `--tree-cut-effort` knob is live but only bites once node separation
/// has spent the offset (measured): on int_production_mixed_s1 effort
/// 0 vs 0.1 is 90,263 nodes (time-limit) vs 26,065 solved, while on fcn the
/// cap is never reached and every effort prints identically, so it looks inert there.
/// Re-measured effort 0.02 vs 0.1 at aging limit 10: lights_out improves
/// (1.58 → 1.16 s), int_production and fcn s3 regress (14.8 → 18.0 s,
/// 4.8 → 8.5 s); 0.1 stays the default.
fn tree_separation_allowed(cc: &CutConfig, stats: &CutStats, total_lp_iterations: u64) -> bool {
    cc.tree
        && (stats.tree_lp_iterations as f64)
            < cc.tree_effort * total_lp_iterations as f64 + TREE_EFFORT_OFFSET as f64
}

/// Runs the root cut loop on `lp` (the optimal, fractional root LP under
/// `lower`/`upper`). Adds the selected cuts to `engine` and returns the final
/// root LP plus the separator state for the tree (`None` when the engine
/// cannot hold cuts or cuts are off; `lp` is then returned unchanged).
/// `upper_limit` is the search's current cutoff (internal space,
/// `f64::INFINITY` without an incumbent); with `CutConfig::root_cutoff_row`
/// it enters the loop as the objective row, and `heur` — called on every
/// round's LP point while no row has landed, and (under
/// [`CutConfig::root_cutoff_tighten`) after as well — proposes rounding points
/// (with their [`Source`] and the LP iterations its repair step spent) from
/// which the loop verifies one, adds the row, and hands all of them back in
/// [`CutStats::root_incumbents`] for the search to adopt as incumbents.
#[allow(clippy::too_many_arguments)]
pub(crate) fn root_separation<E: LpEngine>(
    p: &MipProblem,
    engine: &mut E,
    lower: &[f64],
    upper: &[f64],
    lp: LpOutcome,
    config: &MipConfig,
    stop: &dyn Fn() -> bool,
    stats: &mut CutStats,
    upper_limit: f64,
    heur: &mut RootHeur<'_, E>,
) -> Result<(LpOutcome, Option<Separator>), LpError> {
    root_separation_logged(
        p,
        engine,
        lower,
        upper,
        lp,
        config,
        stop,
        stats,
        None,
        upper_limit,
        heur,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn root_separation_logged<E: LpEngine>(
    p: &MipProblem,
    engine: &mut E,
    lower: &[f64],
    upper: &[f64],
    mut lp: LpOutcome,
    config: &MipConfig,
    stop: &dyn Fn() -> bool,
    stats: &mut CutStats,
    mut log: Option<&mut Vec<(CutClass, CutRow)>>,
    upper_limit: f64,
    heur: &mut RootHeur<'_, E>,
) -> Result<(LpOutcome, Option<Separator>), LpError> {
    let start = Instant::now();
    let cc = config.cuts;
    let feastol = config.feastol;
    if !cc.enabled || lp.status != LpStatus::Optimal || !engine.add_cuts(&[]) {
        return Ok((lp, None));
    }
    stats.root_bound_before = p.caller_objective(lp.objective);
    stats.root_bound_after = stats.root_bound_before;

    let mut sep = Separator::new(p, lower, upper, &cc, feastol);
    let max_rounds = cc
        .max_rounds
        .unwrap_or_else(|| default_max_rounds(p, lower, upper));

    // The search's incumbent is already in hand: its cutoff row joins the
    // relaxation before round 1 (slack at `x*`, so no re-solve).
    if cc.root_cutoff_row && add_cutoff_row(&mut sep, engine, p, upper_limit) {
        stats.cutoff_rows += 1;
        if std::env::var("RUSOLVE_CUT_DEBUG").is_ok() {
            eprintln!("cutoff row from entry incumbent {upper_limit:.6}");
        }
    }

    let first_obj = lp.objective;
    let first_x = lp.x.clone();
    let mut prev_obj = first_obj;
    let mut avg_dir = vec![0.0; p.num_col];
    let mut smooth_progress = 0.0;
    let mut stall = 0;
    let mut weak = 0;
    // Landing grants no extra rounds. An earlier probe suggested the row
    // re-energizes the tableau separator, and the first shipped build
    // granted 10 grace rounds + a stall/weak reset on landing; measured
    // (network_design s3), that stretched the root loop 9 → 19 rounds
    // for +0.8 bound and +60 ms wall — the row's job (capture the
    // incumbent and its bound) is done at landing, and HiGHS's incumbent
    // cutoff likewise enters as an objective bound without extending the
    // cut loop.
    while stall < 3 && has_fractional(p, &lp.x, config.integrality_tol) {
        if stats.rounds >= max_rounds || stop() {
            break;
        }
        stats.rounds += 1;
        let next = match sep.round(
            engine,
            lower,
            upper,
            &lp,
            &cc,
            feastol,
            stats,
            log.as_deref_mut(),
        ) {
            Ok(Some(next)) => next,
            Ok(None) => {
                if std::env::var("RUSOLVE_CUT_DEBUG").is_ok() {
                    eprintln!("round {}: no cuts selected — loop ends", stats.rounds);
                }
                break;
            }
            Err(e) => {
                lp_failed(stats, stats.rounds, &e);
                break;
            }
        };
        lp = next;
        if lp.status != LpStatus::Optimal {
            // Infeasible: the cuts prove no integer point exists.
            break;
        }
        stats.root_bound_after = p.caller_objective(lp.objective);
        if std::env::var("RUSOLVE_CUT_DEBUG").is_ok() {
            let ngen: u64 = stats.generated.iter().sum();
            eprintln!(
                "round {}: obj {:.6} gen {} stall {} weak {}",
                stats.rounds, lp.objective, ngen, stall, weak
            );
        }

        // HiGHS removes obsolete (basic, y = 0) cut rows after every root
        // round's re-solve (HighsMipSolverData.cpp:2098), not only before
        // branching: the next round then separates over the smaller LP, at
        // the same point but with a different row mix in the pool. As in
        // HiGHS (`removeCuts` → `optimizeLp`) the drop re-optimises: the
        // bound cannot move (only zero-dual rows leave), the vertex may.
        if cc.efficacy_rounds
            && let Some(mask) = obsolete_row_mask(&mut *engine, &lp, sep.lp_cuts.len(), feastol)
            && drop_obsolete_rows(&mut sep, engine, &mask)
        {
            // Only zero-dual basic rows left, so on a failed re-solve `lp`
            // stays optimal for the smaller LP and its basis still fits.
            match engine.solve(lower, upper) {
                Ok(next) => lp = next,
                Err(e) => {
                    lp_failed(stats, stats.rounds, &e);
                    break;
                }
            }
            stats.lp_iterations += lp.iterations;
            if lp.status == LpStatus::Optimal {
                sep.relax.set_solution(&lp.x);
                stats.root_bound_after = p.caller_objective(lp.objective);
            }
        }

        // Mid-loop incumbent: the
        // loop runs blind while no cutoff exists, and blind it saturates its
        // pool on int_production at 71 % of the gap while HiGHS, whose root
        // has an incumbent from 0.1 s, reaches 92.6 %. Round the LP point
        // through the search's rounding ladder (HiGHS likewise interleaves
        // heuristics through the root loop); a verified point lands the
        // cutoff row — the loop itself is not extended. With
        // `root_cutoff_tighten` the ladder keeps running and every better
        // point tightens the row; by default it stops at the first landing
        // (further tightenings steer int_production's tree worse than
        // the single landing).
        if cc.root_cutoff_row && (!sep.cutoff_row || cc.root_cutoff_tighten) {
            let (points, heur_lps) = heur(&lp.x, engine, sep.cutoff_ub);
            if heur_lps > 0 {
                // The repair LPs left the engine's basis on a rounded box;
                // restore the cut LP's basis before the next round (same
                // box, warm from `lp.basis`: a handful of pivots at most).
                // Same rows as `lp`: a failed restore keeps it.
                match engine.solve_from(lower, upper, lp.basis.as_deref(), f64::INFINITY) {
                    Ok(next) => lp = next,
                    Err(e) => {
                        lp_failed(stats, stats.rounds, &e);
                        break;
                    }
                }
                stats.lp_iterations += lp.iterations;
                if lp.status != LpStatus::Optimal {
                    break;
                }
                sep.relax.set_solution(&lp.x);
                stats.root_bound_after = p.caller_objective(lp.objective);
            }
            let mut best: Option<f64> = None;
            let n_points = points.len();
            for (y, src) in points {
                if let Some(ub) = verified_internal_objective(p, &y, config.integrality_tol) {
                    if best.is_none_or(|b| ub < b) {
                        best = Some(ub);
                    }
                    stats.root_incumbents.push((y, src));
                }
            }
            if let Some(ub) = best
                && add_cutoff_row(&mut sep, engine, p, ub)
            {
                let first = stats.cutoff_rows == 0;
                stats.cutoff_rows += 1;
                if std::env::var("RUSOLVE_CUT_DEBUG").is_ok() {
                    eprintln!(
                        "round {}: cutoff row {} at {ub:.6} ({n_points} heuristic points)",
                        stats.rounds,
                        if first { "added" } else { "tightened" }
                    );
                }
                // The first row enters cold (warm on the current basis
                // leaves the row separator-invisible; the cold landing
                // trajectory, not the incumbent alone, carries
                // int_production's tree win — all-warm was measured and s1
                // went 90,499 → TimeLimit). Cold costs a full root
                // re-solve, so it is only paid while the loop can monetize
                // it: with at least half the round budget left
                // (network_design s3 lands at round 9 of 15 and paid +30 ms
                // on a 54 ms solve for a row that changed nothing; landed
                // warm there, the instance runs 1.02x base and the row is
                // stripped at the end anyway). Later landings and
                // tightenings land warm.
                let cold = first && stats.rounds * 2 <= max_rounds;
                match engine.solve_from(
                    lower,
                    upper,
                    if cold { None } else { lp.basis.as_deref() },
                    f64::INFINITY,
                ) {
                    Ok(next) => lp = next,
                    Err(e) => {
                        // The row is slack at `lp`: take it back out and
                        // stop on `lp`. Its point still reaches the
                        // search through `root_incumbents`.
                        let n = sep.lp_cuts.len();
                        let keep: Vec<bool> = (0..n).map(|i| i + 1 < n).collect();
                        engine.retain_cuts(&keep);
                        sep.relax.retain_cuts(&keep);
                        sep.relax.set_solution(&lp.x);
                        sep.lp_cuts.pop();
                        sep.cutoff_row = sep.lp_cuts.iter().any(|e| e.pool.is_none());
                        lp_failed(stats, stats.rounds, &e);
                        break;
                    }
                }
                stats.lp_iterations += lp.iterations;
                if lp.status != LpStatus::Optimal {
                    break;
                }
                sep.relax.set_solution(&lp.x);
                stats.root_bound_after = p.caller_objective(lp.objective);
                if std::env::var("RUSOLVE_CUT_DEBUG").is_ok() {
                    eprintln!("  re-solve: obj {:.6}", lp.objective);
                }
            } else if std::env::var("RUSOLVE_CUT_DEBUG").is_ok() && n_points > 0 {
                eprintln!(
                    "round {}: heuristic proposed {n_points} points, none verified",
                    stats.rounds
                );
            }
        }

        // Stall detection (HighsMipSolverData.cpp:2228–2266).
        let n = stats.rounds as f64;
        let cur: Vec<f64> = first_x.iter().zip(&lp.x).map(|(a, b)| a - b).collect();
        let sq: f64 = cur.iter().map(|v| v * v).sum();
        let progress = if sq > 0.0 {
            let scale = 1.0 / sq.sqrt();
            let (mut asq, mut dot) = (0.0, 0.0);
            for (a, &c) in avg_dir.iter_mut().zip(&cur) {
                *a = (scale * c - *a) / n;
                asq += *a * *a;
                dot += *a * c;
            }
            if asq > 0.0 { dot / asq.sqrt() } else { 0.0 }
        } else {
            0.0
        };
        if stats.rounds == 1 {
            smooth_progress = progress;
        } else {
            let next = (2.0 / 3.0) * smooth_progress + progress / 3.0;
            if next < smooth_progress * 1.01
                && (lp.objective - first_obj) <= (prev_obj - first_obj) * 1.001
            {
                stall += 1;
            } else {
                stall = 0;
            }
            smooth_progress = next;
        }
        // Low improvement relative to the gain so far.
        if weak_round(first_obj, prev_obj, lp.objective, cc.min_rel_improvement) {
            weak += 1;
            let budget = if moved(first_obj, lp.objective) {
                cc.weak_rounds
            } else {
                cc.flat_weak_rounds
            };
            if weak >= budget {
                break;
            }
        } else {
            weak = 0;
        }
        prev_obj = lp.objective;
    }

    // The cutoff row is a root-loop device only. It is slack at the
    // final root LP (its rhs is an incumbent value, above the relaxation
    // bound) and structurally inert in the tree — node LPs re-solve under
    // an objective cutoff ≤ the row's rhs, so its dual price is zero —
    // yet it taxes every node LP as one dense objective row (medium tier: market_split 1.12× time
    // at 10,990 → 14,346 nodes). The point it carries has already reached
    // the search as an incumbent (`cut_stats.root_incumbents`, adopted in
    // `bnb.rs` after this function returns), so strip the row before
    // branching: keep the steering, drop the tax. Running first also keeps
    // every drop path below looking only at pool cuts.
    let mut stripped = false;
    if cc.root_cutoff_row
        && lp.status == LpStatus::Optimal
        && sep.lp_cuts.iter().any(|e| e.pool.is_none())
    {
        let keep: Vec<bool> = sep.lp_cuts.iter().map(|e| e.pool.is_some()).collect();
        engine.retain_cuts(&keep);
        sep.relax.retain_cuts(&keep);
        sep.lp_cuts.retain(|e| e.pool.is_some());
        sep.cutoff_row = false;
        // Warm: the row is slack with a zero dual, so removing it
        // leaves the current basis dual feasible and the re-solve is a
        // few-pivot formality that keeps `lp.basis` valid for the tree.
        // A failed re-solve keeps `lp`: still optimal without the slack
        // row, and its basis still fits.
        match engine.solve_from(lower, upper, lp.basis.as_deref(), f64::INFINITY) {
            Ok(next) => {
                lp = next;
                stats.lp_iterations += lp.iterations;
                if lp.status == LpStatus::Optimal {
                    sep.relax.set_solution(&lp.x);
                    stats.root_bound_after = p.caller_objective(lp.objective);
                }
            }
            Err(e) => lp_failed(stats, stats.rounds, &e),
        }
        stripped = true;
    }
    if stripped && std::env::var("RUSOLVE_CUT_DEBUG").is_ok() {
        eprintln!("cutoff row stripped before branching");
    }

    if cc.drop_if_no_gain
        && lp.status == LpStatus::Optimal
        && !sep.lp_cuts.is_empty()
        && !moved(first_obj, lp.objective)
    {
        let keep: Vec<bool> = sep.lp_cuts.iter().map(|e| e.pool.is_none()).collect();
        engine.retain_cuts(&keep);
        sep.relax.retain_cuts(&keep);
        for e in &sep.lp_cuts {
            if let Some(c) = e.pool {
                sep.pool.lp_cut_removed(c);
            }
        }
        sep.lp_cuts.retain(|e| e.pool.is_none());
        // On a failed re-solve `lp` stays: the dropped cuts are valid
        // inequalities, so its objective is still a valid dual bound, and
        // its basis carries `cut_ids` to fit the rows that remain.
        match engine.solve(lower, upper) {
            Ok(next) => {
                lp = next;
                stats.lp_iterations += lp.iterations;
                sep.relax.set_solution(&lp.x);
            }
            Err(e) => lp_failed(stats, stats.rounds, &e),
        }
    }

    // TEMP(cut-rejection probe): drop every cut row before branching.
    if cc.drop_all_at_end && lp.status == LpStatus::Optimal && !sep.lp_cuts.is_empty() {
        let keep: Vec<bool> = sep.lp_cuts.iter().map(|e| e.pool.is_none()).collect();
        engine.retain_cuts(&keep);
        sep.relax.retain_cuts(&keep);
        for e in &sep.lp_cuts {
            if let Some(c) = e.pool {
                sep.pool.lp_cut_removed(c);
            }
        }
        sep.lp_cuts.retain(|e| e.pool.is_none());
        match engine.solve(lower, upper) {
            Ok(next) => {
                lp = next;
                stats.lp_iterations += lp.iterations;
                sep.relax.set_solution(&lp.x);
                stats.root_bound_after = p.caller_objective(lp.objective);
            }
            Err(e) => lp_failed(stats, stats.rounds, &e),
        }
    }

    // TEMP(cut-rejection probe): drop only Gomory-class rows before branching.
    if cc.drop_gomory_at_end && lp.status == LpStatus::Optimal && !sep.lp_cuts.is_empty() {
        let keep: Vec<bool> = sep
            .lp_cuts
            .iter()
            .map(|e| {
                e.pool
                    .is_none_or(|c| sep.pool.cuts[c].class != CutClass::Gomory)
            })
            .collect();
        if keep.iter().any(|k| !k) {
            engine.retain_cuts(&keep);
            sep.relax.retain_cuts(&keep);
            let mut k = 0;
            sep.lp_cuts.retain(|e| {
                k += 1;
                if !keep[k - 1] {
                    sep.pool
                        .lp_cut_removed(e.pool.expect("Gomory rows are pool cuts"));
                }
                keep[k - 1]
            });
            match engine.solve(lower, upper) {
                Ok(next) => {
                    lp = next;
                    stats.lp_iterations += lp.iterations;
                    sep.relax.set_solution(&lp.x);
                    stats.root_bound_after = p.caller_objective(lp.objective);
                }
                Err(e) => lp_failed(stats, stats.rounds, &e),
            }
        }
    }

    // Before branching, HiGHS removes obsolete rows once more
    // (HighsMipSolverData.cpp:2481). The efficacy variant drops every row
    // with a zero dual, where the slack test keeps tight-but-inefficacious
    // (degenerate) rows; without engine duals it falls back to slack. With
    // duals the slack rule is subsumed: surviving rows are nonbasic with a
    // nonzero dual, hence tight.
    let mut efficacy_end_applied = false;
    if cc.efficacy_end
        && lp.status == LpStatus::Optimal
        && !sep.lp_cuts.is_empty()
        && let Some(mask) = obsolete_row_mask(&mut *engine, &lp, sep.lp_cuts.len(), feastol)
    {
        efficacy_end_applied = true;
        if drop_obsolete_rows(&mut sep, engine, &mask) {
            // Zero-dual basic rows only: a failed re-solve keeps `lp`.
            match engine.solve(lower, upper) {
                Ok(next) => {
                    lp = next;
                    stats.lp_iterations += lp.iterations;
                    if lp.status == LpStatus::Optimal {
                        sep.relax.set_solution(&lp.x);
                        stats.root_bound_after = p.caller_objective(lp.objective);
                    }
                }
                Err(e) => lp_failed(stats, stats.rounds, &e),
            }
        }
    }

    let relax = &mut sep.relax;
    if cc.drop_slack_at_end
        && !efficacy_end_applied
        && lp.status == LpStatus::Optimal
        && !sep.lp_cuts.is_empty()
    {
        let keep: Vec<bool> = (0..sep.lp_cuts.len())
            .map(|k| {
                sep.lp_cuts[k].pool.is_none()
                    || relax.activity[relax.num_model_rows + k]
                        >= relax.rows[relax.num_model_rows + k].upper - feastol
            })
            .collect();
        engine.retain_cuts(&keep);
        relax.retain_cuts(&keep);
        let x = relax.x.clone();
        relax.set_solution(&x);
        let mut k = 0;
        sep.lp_cuts.retain(|e| {
            k += 1;
            if !keep[k - 1] {
                sep.pool
                    .lp_cut_removed(e.pool.expect("slack-dropped rows are pool cuts"));
            }
            keep[k - 1]
        });
    }
    for e in &sep.lp_cuts {
        if let Some(c) = e.pool {
            stats.in_lp[sep.pool.cuts[c].class as usize] += 1;
        }
    }
    sep.root_obj = lp.objective;
    // The Gomory park is a root-loop rule (TableauState::parked): node
    // separation (`--tree-cuts`) starts with the separator live again.
    sep.tableau.parked = false;
    // In-loop probing derives globally valid cliques from the root box
    // only; node separation would probe under node bounds, so it ends here.
    sep.probe = None;
    stats.time_s = start.elapsed().as_secs_f64();
    Ok((lp, Some(sep)))
}

/// Which LP cut rows are obsolete, aligned to `lp_cuts` order: HiGHS's
/// `removeObsoleteRows` test is `row_status == kBasic`
/// (HighsLpRelaxation.cpp:540, MIT; read as a reference text); where the
/// engine gives no basis status, a zero sign-clipped row dual is the proxy —
/// additionally true for tight-but-degenerate nonbasic rows, which HiGHS
/// keeps. `None` when the engine offers neither.
fn obsolete_row_mask<E: LpEngine>(
    engine: &mut E,
    lp: &LpOutcome,
    n_lp: usize,
    dual_tol: f64,
) -> Option<Vec<bool>> {
    if let Some(basic) = engine.cut_rows_basic(lp)
        && basic.len() == n_lp
    {
        return Some(basic);
    }
    let d = engine.duals(lp)?;
    let y = &d.row_duals;
    // The cut rows are the LAST `n_lp` dual entries: the engine's row-form
    // template can hold more model rows than the relaxation has (a ranged
    // row becomes two half-space rows), so only the tail index is
    // meaningful.
    if y.len() < n_lp {
        return None;
    }
    Some(
        y[y.len() - n_lp..]
            .iter()
            .map(|yi| yi.abs() <= dual_tol)
            .collect(),
    )
}

/// Drops the flagged LP cut rows. Dropping an obsolete row changes neither
/// x* nor the objective (its dual is zero: basic, or degenerate at bound),
/// and the dropped cut re-enters the pool at age 1 (`lpCutRemoved`, :550).
/// Returns whether any row left; HiGHS re-optimises after `removeCuts`
/// (`optimizeLp`, HighsLpRelaxation.cpp:600) to keep the basis valid, and
/// callers must do the same re-solve.
fn drop_obsolete_rows<E: LpEngine>(sep: &mut Separator, engine: &mut E, drop: &[bool]) -> bool {
    if drop.len() != sep.lp_cuts.len() || drop.iter().all(|&d| !d) {
        return false;
    }
    let keep: Vec<bool> = drop
        .iter()
        .zip(&sep.lp_cuts)
        .map(|(&d, e)| e.pool.is_none() || !d)
        .collect();
    if keep.iter().all(|&k| k) {
        return false;
    }
    engine.retain_cuts(&keep);
    sep.relax.retain_cuts(&keep);
    let x = sep.relax.x.clone();
    sep.relax.set_solution(&x);
    let mut k = 0;
    sep.lp_cuts.retain(|e| {
        k += 1;
        if drop[k - 1]
            && let Some(c) = e.pool
        {
            sep.pool.lp_cut_removed(c);
            false
        } else {
            keep[k - 1]
        }
    });
    true
}

/// HighsLpRelaxation::performAging(true): every epoch, LP cuts with slack
/// (our stand-in for "basic"; no basis status yet) age; old ones leave.
/// `agelimit` is HiGHS `mip_lp_age_limit` (default 10); the periodic check
/// runs only each `agelimit/2` epochs (HighsLpRelaxation.cpp:646–653).
fn age_lp_cuts<E: LpEngine>(
    relax: &mut Relaxation,
    engine: &mut E,
    pool: &mut CutPool,
    lp_cuts: &mut Vec<LpCut>,
    epochs: u32,
    agelimit: u32,
    feastol: f64,
) -> usize {
    let agelimit = if !epochs.is_multiple_of((agelimit / 2).max(2)) {
        u32::MAX
    } else {
        epochs.min(agelimit)
    };
    let m0 = relax.num_model_rows;
    let mut keep = vec![true; lp_cuts.len()];
    let mut any = false;
    for (k, entry) in lp_cuts.iter_mut().enumerate() {
        let Some(c) = entry.pool else {
            continue; // the cutoff row never ages (it is the incumbent cutoff)
        };
        let row = &relax.rows[m0 + k];
        if relax.activity[m0 + k] < row.upper - feastol {
            entry.age += 1;
            if entry.age > agelimit {
                keep[k] = false;
                any = true;
                pool.lp_cut_removed(c);
            }
        } else {
            entry.age = 0;
        }
    }
    if any {
        engine.retain_cuts(&keep);
        relax.retain_cuts(&keep);
        let x = relax.x.clone();
        relax.set_solution(&x);
        let mut k = 0;
        lp_cuts.retain(|_| {
            k += 1;
            keep[k - 1]
        });
    }
    keep.iter().filter(|&&k| !k).count()
}

/// The bound moved beyond round-off between `from` and `to` (minimization).
fn moved(from: f64, to: f64) -> bool {
    to - from > 1e-9 * to.abs().max(1.0)
}

/// A root round from `prev` to `cur` is weak when it adds less than `rel` of
/// the gain since `first`, or does not move the bound at all.
fn weak_round(first: f64, prev: f64, cur: f64, rel: f64) -> bool {
    !moved(prev, cur) || cur - prev < rel * (cur - first)
}

fn has_fractional(p: &MipProblem, x: &[f64], tol: f64) -> bool {
    x.iter()
        .zip(&p.integer)
        .any(|(&v, &int)| int && (v - v.round()).abs() > tol)
}

/// All columns integer with integral coefficients ⇒ integral activity.
fn row_integral(relax: &Relaxation, inds: &[usize], vals: &[f64]) -> bool {
    inds.iter()
        .zip(vals)
        .all(|(&j, &v)| relax.integer[j] && (v - v.round()).abs() <= EPSILON)
}

/// `2·√(Σ_int ⌈log2(min(1024, 1 + u − l))⌉)` (HighsMipSolverData.cpp:1045,
/// :1990–1993).
fn default_max_rounds(p: &MipProblem, lower: &[f64], upper: &[f64]) -> u32 {
    let mut log2 = 0.0;
    for j in 0..p.num_col {
        if p.integer[j] {
            log2 += (1024f64.min(1.0 + upper[j] - lower[j])).log2().ceil();
        }
    }
    (2.0 * f64::sqrt(log2)) as u32
}
