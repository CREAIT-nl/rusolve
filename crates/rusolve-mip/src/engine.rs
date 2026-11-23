//! Node-relaxation LP interface.
//!
//! The tree search only ever asks "solve the relaxation under these column
//! bounds". Keeping that behind [`LpEngine`] lets warm starts, dual simplex
//! reoptimization, duals and cut rows arrive from the LP solver without touching
//! the search code.
//!
//! Warm starts: every optimal solve hands back
//! a [`WarmStart`] snapshot; children re-solve from their parent's snapshot
//! and cut rounds from the previous round's, by warm dual simplex with the
//! incumbent cutoff (`rusolve_lp::dual::reoptimize`; usage after
//! the 3-node ladder proof 45516c4). The scheme follows HiGHS (MIT, read as a
//! reference text): plunge children carry the parent basis
//! (`NodeData::nodeBasis`, reference/HiGHS/highs/mip/HighsSearch.cpp:960,
//! :1855–1863, used only while its row count still fits), node LPs stop on an
//! objective bound (`HighsLpRelaxation::setObjectiveLimit`,
//! HighsLpRelaxation.cpp:1114; `kObjectiveBound` prunes, :1263), new cut rows
//! enter with their slack basic (`addCuts`, :509) and only basic-slack cut
//! rows are removed (`removeObsoleteRows`, :536).
//!
//! Since E-row basis fidelity landed in the LP solver, it accepts every
//! snapshot it exported, so the engine no longer stops offering warm starts
//! after repeated rejections. A warm Infeasible that comes with a verified
//! Farkas ray (`ReoptResult::dual_ray`) prunes the node as is; without a
//! ray it is still re-checked by a cold solve.
//!
//! Duals: [`LpEngine::duals`] gives the
//! row duals and reduced costs of an optimal relaxation in the internal
//! minimisation sense, for reduced-cost fixing. Cold solves carry the LP's
//! own duals (presolve-off path) and so do warm solves
//! (`ReoptResult::duals_valid`): either
//! is kept with the snapshot of its solve and served with no factorisation.
//! Only when the LP returned none does the engine fall back to solving
//! `B_Tᵀ y_T = c_B` over the basic columns and tight rows of the exported
//! basis (`basis_row_duals`, one sparse LU). Either way the fixing code only
//! uses the Lagrangian bound
//! [`LpDuals::bound`], which is valid for *any* sign-feasible `y`, so a
//! slightly wrong `y` weakens the fixing but never makes it wrong.

use crate::problem::MipProblem;
use rusolve_lp::lu::SparseLu;
use rusolve_lp::{
    Basis, BasisStatus, DseMode, DseWeights, Lp, ReoptSession, ReoptStatus, SimplexConfig,
    SimplexError, SimplexSolver, SimplexStatus,
};
use std::rc::Rc;

pub mod ekk;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LpStatus {
    Optimal,
    Infeasible,
    Unbounded,
    /// The LP cannot beat the cutoff passed to [`LpEngine::solve_from`]:
    /// `objective` is a valid lower bound on the LP optimum above the cutoff
    /// and `x` is empty. Only returned when a finite cutoff was given.
    Cutoff,
    /// A strong-branching probe ([`LpEngine::probe`]) used up its iteration
    /// cap: `objective` is the dual bound reached (a valid lower bound on the
    /// LP optimum; `-∞` when none is known) and `x` is empty.
    IterationLimit,
}

/// A solved node relaxation. `objective` is in the internal minimization
/// sense; `x` is empty unless `status == Optimal`.
#[derive(Debug, Clone)]
pub struct LpOutcome {
    pub status: LpStatus,
    pub objective: f64,
    pub x: Vec<f64>,
    pub iterations: u64,
    /// Optimal basis to warm-start later solves from (children, the next cut
    /// round). `None` when the engine has none to give.
    pub basis: Option<Rc<WarmStart>>,
}

/// An optimal basis of the engine's LP together with the identity of the cut
/// rows it was taken over, so it stays usable after cuts come and go.
#[derive(Debug, Clone)]
pub struct WarmStart {
    basis: Basis,
    /// Dual steepest-edge weights of the basic variables of `basis`, when
    /// the solve that produced it ran DSE (`RusolveEngine::set_warm_dse`).
    dse: Option<DseWeights>,
    /// Ids of the cut rows present at snapshot time, in row order.
    cut_ids: Rc<[u64]>,
}

/// How node and cut-round LPs were solved (warm dual reopt vs cold primal).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LpStats {
    /// Solves started cold: no snapshot to start from, or warm start off.
    pub cold: u64,
    /// Solves answered by a warm dual reoptimization.
    pub warm: u64,
    /// Snapshots that did not fit the LP (a tight cut row removed, a status
    /// the new bounds cannot hold); the solve ran cold.
    pub fallback_unfit: u64,
    /// Under `RusolveEngine::set_unfit_last` (off by default): unfit
    /// snapshots replaced by the engine's latest basis instead of a cold
    /// solve (HiGHS keeps its current basis when a node's stored one has the
    /// wrong row count, HighsSearch.cpp:790). Not counted in
    /// `fallback_unfit`.
    pub unfit_recovered: u64,
    /// Snapshots the LP's warm entry rejected (singular, not dual feasible)
    /// or warm runs its dual engine abandoned; `reoptimize` solved cold.
    pub fallback_rejected: u64,
    /// Warm runs that hit the iteration limit; re-solved cold.
    pub fallback_iter_limit: u64,
    /// Warm runs that failed (LP error, Unbounded, a point violating the
    /// rows); re-solved cold.
    pub fallback_error: u64,
    /// Warm Infeasible verdicts re-checked by a cold solve, and how many of
    /// those the cold solve overturned.
    pub infeasible_checks: u64,
    pub infeasible_overturned: u64,
    /// Warm Infeasible verdicts that came with a verified Farkas ray
    /// (`ReoptResult::dual_ray`) and pruned without a cold re-check.
    pub infeasible_trusted: u64,
    /// Failed snapshot-free cold solves re-run on the other simplex
    /// ([`LpEngine::set_root_lp`] `retry`).
    pub root_retries: u64,
    /// Counted by the search rather than the engine: tree nodes whose
    /// LP did not solve and were branched on without one, and
    /// all-integers-fixed nodes pruned after their fallback re-solve failed
    /// as well (`bnb::Search::unsolved_node`).
    pub no_lp_branches: u64,
    pub unsolved_pruned: u64,
    /// Warm solves stopped by the objective cutoff.
    pub cutoffs: u64,
    pub warm_iterations: u64,
    pub cold_iterations: u64,
    /// Strong-branching probes ([`LpEngine::probe`]), how many of them stopped
    /// at their iteration cap, and their iterations (also counted in
    /// `warm_iterations` / `cold_iterations`).
    pub probes: u64,
    pub probe_iter_limits: u64,
    pub probe_iterations: u64,
    /// [`LpEngine::duals`] answers: the LP's own duals (cold solves), duals
    /// recomputed from the exported basis (warm solves), and requests that
    /// got none (basis count mismatch, singular or inaccurate system).
    pub duals_lp: u64,
    pub duals_basis: u64,
    pub duals_failed: u64,
    /// Audit (`MipConfig::audit`, off by default): wall nanoseconds
    /// spent inside [`LpEngine::solve_from`] / [`LpEngine::probe`], and the
    /// probe share of it. Zero unless the engine's timing is switched on.
    pub lp_nanos: u64,
    pub probe_nanos: u64,
}

impl LpStats {
    pub fn fallbacks(&self) -> u64 {
        self.fallback_unfit
            + self.fallback_rejected
            + self.fallback_iter_limit
            + self.fallback_error
    }

    /// Adds `other`'s counters to `self`.
    pub fn absorb(&mut self, other: &LpStats) {
        self.cold += other.cold;
        self.warm += other.warm;
        self.fallback_unfit += other.fallback_unfit;
        self.unfit_recovered += other.unfit_recovered;
        self.fallback_rejected += other.fallback_rejected;
        self.fallback_iter_limit += other.fallback_iter_limit;
        self.fallback_error += other.fallback_error;
        self.infeasible_checks += other.infeasible_checks;
        self.infeasible_overturned += other.infeasible_overturned;
        self.infeasible_trusted += other.infeasible_trusted;
        self.root_retries += other.root_retries;
        self.no_lp_branches += other.no_lp_branches;
        self.unsolved_pruned += other.unsolved_pruned;
        self.cutoffs += other.cutoffs;
        self.warm_iterations += other.warm_iterations;
        self.cold_iterations += other.cold_iterations;
        self.probes += other.probes;
        self.probe_iter_limits += other.probe_iter_limits;
        self.probe_iterations += other.probe_iterations;
        self.duals_lp += other.duals_lp;
        self.duals_basis += other.duals_basis;
        self.duals_failed += other.duals_failed;
        self.lp_nanos += other.lp_nanos;
        self.probe_nanos += other.probe_nanos;
    }
}

/// Dual information of an optimal relaxation in the internal minimisation
/// sense: row duals `y` (sign-feasible: `y ≤ 0` on `≤` rows, `y ≥ 0` on `≥`
/// rows) and reduced costs `d = c − Aᵀy` over the relaxation's rows,
/// including the cut rows held when it was solved.
///
/// For every `x` satisfying the rows, `cᵀx = yᵀAx + dᵀx ≥ yᵀb + dᵀx`, so
/// over a box `[l, u]`, `cᵀx ≥ yᵀb + Σ_j min(d_j l_j, d_j u_j)` — the
/// Lagrangian bound [`Self::bound`]. With the optimal `y` it equals the LP
/// objective; the textbook rule follows from it: a column at its lower
/// bound with `d_j > 0` satisfies `x_j ≤ l_j + (cutoff − z)/d_j` in every
/// point with objective at most `cutoff`.
#[derive(Debug, Clone, PartialEq)]
pub struct LpDuals {
    /// `d_j = c_j − yᵀa_j` per column.
    pub reduced_costs: Vec<f64>,
    /// The sign-feasible row duals `y` themselves, in the relaxation's row
    /// order (model rows first, cut rows in engine order after). A cut row's
    /// entry is its efficacy signal: HiGHS reads exactly this to reset row
    /// ages (`HighsLpRelaxation::performAging`, HighsLpRelaxation.cpp:670)
    /// and to drop obsolete rows (`removeObsoleteRows`, :536, via `kBasic`).
    pub row_duals: Vec<f64>,
    /// `|c_j| + Σ_i |y_i a_ij|`: the magnitude `d_j` was computed from.
    scale: Vec<f64>,
    /// `yᵀb`.
    row_term: f64,
    /// `Σ_i |y_i b_i|`.
    row_scale: f64,
}

/// Relative round-off allowance of a computed `d_j` (`1e-9` of the terms it
/// was computed from): below it a reduced cost is treated as zero against
/// an infinite bound, where trusting its sign would make the bound −∞.
const DUAL_NOISE: f64 = 1e-9;

impl LpDuals {
    /// Builds the duals from row duals `y` over the rows of `lp`, clipping
    /// any wrong-signed entry to zero (which keeps [`Self::bound`] valid).
    fn from_rows(lp: &Lp, mut y: Vec<f64>) -> Self {
        for (yi, &sense) in y.iter_mut().zip(&lp.senses) {
            match sense {
                'L' => *yi = yi.min(0.0),
                'G' => *yi = yi.max(0.0),
                _ => {}
            }
        }
        let mut reduced_costs = Vec::with_capacity(lp.num_col);
        let mut scale = Vec::with_capacity(lp.num_col);
        let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
        for j in 0..lp.num_col {
            // πᵀa_j off the CSC column (ascending rows — the dense zip order;
            // stored zeros contribute ±0.0 no-ops to both accumulators).
            let (mut dot, mut mag) = (0.0, 0.0);
            for k in ptr[j]..ptr[j + 1] {
                let t = rval[k] * y[ridx[k] as usize];
                dot += t;
                mag += t.abs();
            }
            reduced_costs.push(lp.objective[j] - dot);
            scale.push(lp.objective[j].abs() + mag);
        }
        let row_term = y.iter().zip(&lp.rhs).map(|(yi, b)| yi * b).sum();
        let row_scale = y.iter().zip(&lp.rhs).map(|(yi, b)| (yi * b).abs()).sum();
        Self {
            reduced_costs,
            row_duals: y,
            scale,
            row_term,
            row_scale,
        }
    }

    /// The Lagrangian lower bound `yᵀb + Σ_j min(d_j l_j, d_j u_j)` on the
    /// objective of every point in the box satisfying the rows, less a
    /// round-off margin. A reduced cost within [`DUAL_NOISE`] of zero against
    /// an infinite bound counts as zero; `None` if a larger one meets an
    /// infinite bound (no finite bound from these duals).
    pub fn bound(&self, lower: &[f64], upper: &[f64]) -> Option<f64> {
        let (mut sum, mut mag) = (self.row_term, self.row_scale);
        for (j, &d) in self.reduced_costs.iter().enumerate() {
            let b = if d > 0.0 {
                lower[j]
            } else if d < 0.0 {
                upper[j]
            } else {
                continue;
            };
            if b.is_finite() {
                sum += d * b;
                mag += (d * b).abs();
            } else if d.abs() > DUAL_NOISE * self.scale[j] {
                return None;
            }
        }
        Some(sum - 1e-12 * (1.0 + mag))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum LpError {
    Simplex(SimplexError),
    IterationLimit,
}

impl std::fmt::Display for LpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LpError::Simplex(e) => write!(f, "{e}"),
            LpError::IterationLimit => write!(f, "LP iteration limit"),
        }
    }
}

impl std::error::Error for LpError {}

/// A globally valid inequality `Σ_k vals[k]·x[inds[k]] ≤ rhs` appended to the
/// relaxation (a cutting plane). Indices are distinct.
#[derive(Debug, Clone, PartialEq)]
pub struct CutRow {
    pub inds: Vec<usize>,
    pub vals: Vec<f64>,
    pub rhs: f64,
}

pub trait LpEngine {
    /// Solves the relaxation with column bounds replaced by `col_lower` /
    /// `col_upper`. Callers guarantee `col_lower[j] ≤ col_upper[j]`. The
    /// engine may start from the basis of its previous optimal solve.
    fn solve(&mut self, col_lower: &[f64], col_upper: &[f64]) -> Result<LpOutcome, LpError>;

    /// [`Self::solve`] starting from `start` (cold when `None`), allowed to
    /// stop with [`LpStatus::Cutoff`] once the LP bound provably exceeds
    /// `cutoff` (`f64::INFINITY` = no cutoff). The default ignores both.
    fn solve_from(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        _start: Option<&WarmStart>,
        _cutoff: f64,
    ) -> Result<LpOutcome, LpError> {
        self.solve(col_lower, col_upper)
    }

    /// A strong-branching probe: [`Self::solve_from`] with the warm run
    /// capped at `iter_limit` iterations and no cold re-solve when the cap is
    /// hit — it reports [`LpStatus::IterationLimit`] with the dual bound
    /// reached instead. HiGHS strong-branches by dual simplex from the node
    /// basis under the objective limit (`Playground::solveLp`,
    /// reference/HiGHS/highs/mip/HighsLpRelaxation.h:129;
    /// HighsSearch.cpp:405, :574). The default is an uncapped `solve_from`.
    fn probe(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        start: Option<&WarmStart>,
        cutoff: f64,
        _iter_limit: u64,
    ) -> Result<LpOutcome, LpError> {
        self.solve_from(col_lower, col_upper, start, cutoff)
    }

    /// Caps the iterations of later cold solves (`None` restores the
    /// configured limit); a cold solve that reaches it returns
    /// [`LpError::IterationLimit`]. Engines without a cap ignore it.
    fn set_iteration_limit(&mut self, _limit: Option<u64>) {}

    /// Root-node cold-solve policy for later solves until the next call:
    /// `dual` runs cold solves on the dual
    /// simplex; `retry` lets a snapshot-free cold solve that fails (iteration
    /// cap or a numerical give-up) re-run once on the other simplex with no
    /// iteration cap, so only the deadline ends it. HiGHS solves MIP LPs by
    /// dual simplex, gives the root LP no iteration limit
    /// (`setIterationLimit()` = kHighsIInf before the first root solve,
    /// HighsMipSolverData.cpp:2018, HighsLpRelaxation.h:370) and toggles the
    /// simplex strategy when an LP fails (HighsLpRelaxation.cpp:1236-1251;
    /// all MIT). Engines without either ignore it.
    fn set_root_lp(&mut self, _dual: bool, _retry: bool) {}

    /// Solve-internal wall-clock deadline for every later solve
    /// (`None` = off). A solve whose pivot loop crosses it returns
    /// [`LpError::IterationLimit`] (or `LpStatus::IterationLimit` for
    /// probes), the same contract the iteration cap has. Engines without a
    /// clock ignore it. HiGHS checks its timer inside the simplex iteration
    /// loop (`HEkk::bailout`, reference/HiGHS/highs/simplex/HEkk.cpp:
    /// 3441–3461, MIT).
    fn set_deadline(&mut self, _deadline: Option<std::time::Instant>) {}

    /// Appends cut rows that every later [`Self::solve`] must respect; they
    /// are valid for the whole tree, so they are never tied to a node.
    /// Returns `false` (and adds nothing) if the engine cannot hold cuts, in
    /// which case callers skip separation. The default holds none.
    fn add_cuts(&mut self, _cuts: &[CutRow]) -> bool {
        false
    }

    /// Drops previously added cut rows: `keep[k]` refers to the k-th cut
    /// currently held, in insertion order (`keep.len()` = number held).
    fn retain_cuts(&mut self, _keep: &[bool]) {}

    /// Audit: rows of the LP as solved (model rows in row form plus
    /// cut rows). The default reports 0.
    fn num_lp_rows(&self) -> usize {
        0
    }

    /// Warm/cold solve counters so far.
    fn stats(&self) -> LpStats {
        LpStats::default()
    }

    /// Row duals and reduced costs of the optimal relaxation `lp` returned
    /// by this engine (its basis must still fit the engine's current rows),
    /// or `None` when the engine cannot provide them. The default has none.
    fn duals(&mut self, _lp: &LpOutcome) -> Option<LpDuals> {
        None
    }

    /// Per cut row currently held, in engine order, whether its slack is
    /// BASIC in the optimal basis of `lp` — HiGHS's exact `removeObsoleteRows`
    /// / `performAging` test (`row_status[i] == kBasic`,
    /// HighsLpRelaxation.cpp:540, :659). `None` when the basis does not
    /// line up with the current rows or the engine has none.
    fn cut_rows_basic(&mut self, _lp: &LpOutcome) -> Option<Vec<bool>> {
        None
    }

    /// The Farkas proof of the latest solve when it returned
    /// [`LpStatus::Infeasible`] with a verified ray: `(coefs, rhs)` with
    /// `Σ_j coefs[j]·x_j ≤ rhs` valid for every point satisfying the rows
    /// (cuts included), and violated by the column bounds of that solve.
    /// `None` without a ray. The default has none.
    fn infeasibility_proof(&self) -> Option<(Vec<f64>, f64)> {
        None
    }

    /// The dual-bound proof of the latest solve when it returned
    /// [`LpStatus::Cutoff`] (HiGHS `computeDualProof` after
    /// `kObjectiveBound`, HighsLpRelaxation.cpp:732–790): `(coefs, rhs)` with
    /// `Σ_j coefs[j]·x_j ≤ rhs` valid for every point that satisfies the rows
    /// (cuts included) with objective at most the cutoff of that solve, and
    /// violated by the column bounds of that solve. `None` without one.
    fn cutoff_proof(&self) -> Option<(Vec<f64>, f64)> {
        None
    }

    /// For each column of `cols` that is basic in the basis of the optimal
    /// relaxation `lp` (solved under `lower`/`upper`), the row weights of its
    /// tableau row ([`TableauAgg`]); nonbasic columns are omitted. `None`
    /// when the engine has no tableau access or the basis does not fit the
    /// current rows. The default has none.
    fn tableau_weights(
        &mut self,
        _lp: &LpOutcome,
        _lower: &[f64],
        _upper: &[f64],
        _cols: &[usize],
    ) -> Option<Vec<TableauAgg>> {
        None
    }
}

/// The row weights `y = e_posᵀB⁻¹` of the tableau row whose basic variable is
/// column `col`, over the separator's rows (model rows by index, then LP cut
/// rows in order): the aggregation `Σ y_r (a_rᵀx − s_r) = 0` has the basic
/// column at coefficient ±1 and, with the other basic variables eliminated,
/// is that tableau row (HighsTableauSeparator.cpp:143–201).
#[derive(Debug, Clone, PartialEq)]
pub struct TableauAgg {
    pub col: usize,
    pub weights: Vec<(usize, f64)>,
}

/// Relative row/bound violation above which a warm "optimal" point is not
/// trusted and the LP is re-solved cold.
const WARM_POINT_TOL: f64 = 1e-6;

/// Largest `rows × (cols + rows)` for which the engine answers tableau rows
/// at all. It gated the old full-position scan (`TABLEAU_SCAN_LIMIT` named
/// it); the per-candidate BTRAN path keeps the same gate until its
/// relaxation has a measured A/B of its own.
const TABLEAU_SCAN_LIMIT: usize = 60_000_000;

/// Floor of the warm iteration limit (see [`RusolveEngine::solve_warm`]).
const MIN_WARM_ITERATIONS: u64 = 200;

/// [`LpEngine`] over `rusolve-lp`.
///
/// The dense row-form [`Lp`] is built once; only bounds change per node.
/// Two-sided rows become an `L` and a `G` row, free rows are dropped. A solve
/// with a [`WarmStart`] runs warm dual simplex (`dual::reoptimize`) with the
/// cutoff; one without runs the cold two-phase primal. Every warm failure
/// falls back to the cold solve and is counted in [`LpStats`].
pub struct RusolveEngine {
    template: Lp,
    /// The configured LP settings.
    config: SimplexConfig,
    /// `config` with `max_iterations` lowered to the cap of
    /// [`LpEngine::set_iteration_limit`]; what `solver` runs (with the
    /// [`LpEngine::set_root_lp`] policy applied on top).
    capped: SimplexConfig,
    /// Root policy ([`LpEngine::set_root_lp`]): cold solves on the dual
    /// simplex, and a failed snapshot-free cold solve retried uncapped on
    /// the other simplex.
    root_dual: bool,
    root_retry: bool,
    solver: SimplexSolver,
    /// Rows of `template` that come from the model; cut rows follow them.
    num_model_rows: usize,
    /// The model row each of those rows came from, and the model's row count
    /// (two-sided rows give two LP rows, free rows none).
    row_origin: Vec<usize>,
    num_orig_rows: usize,
    /// Id of each cut row in `template`, in row order (ids only grow, so the
    /// list is sorted).
    cut_ids: Rc<[u64]>,
    next_cut_id: u64,
    /// Snapshot of the latest optimal solve: where [`LpEngine::solve`] starts.
    last: Option<Rc<WarmStart>>,
    /// Warm starts on (`false` = every solve cold, the pre-warm behaviour).
    warm: bool,
    /// Warm solves price with dual steepest edge, carrying the weights in
    /// the snapshot (`false` = Dantzig, the historic warm behaviour).
    warm_dse: bool,
    /// An unfit snapshot warm-starts from `last` instead of solving
    /// cold (see [`Self::set_unfit_last`]).
    unfit_last: bool,
    /// Most iterations any cold solve took (sets the warm iteration limit).
    max_cold_iterations: u64,
    /// Deadline handed to the LP solver (see [`Self::set_deadline`]);
    /// re-applied to every rebuilt solver so the iteration-cap swaps below
    /// never drop it.
    lp_deadline: Option<std::time::Instant>,
    stats: LpStats,
    /// Row duals of the latest optimal solve (cold, or warm with
    /// `ReoptResult::duals_valid`) and the snapshot of that solve, so
    /// [`LpEngine::duals`] can hand them out without a solve.
    lp_duals: Option<(Rc<WarmStart>, Vec<f64>)>,
    /// Farkas ray of the latest solve when it ended Infeasible with a
    /// verified one (rows of `template` at that time); cleared by every solve.
    last_ray: Option<Vec<f64>>,
    /// Row multipliers `y` and the cutoff `U` of the latest warm solve that
    /// ended `Cutoff` with a proof; cleared by every solve.
    last_cutoff: Option<(Vec<f64>, f64)>,
    /// DualStandard kept warm across node LPs (bounds-only
    /// edits refresh it in place; see `rusolve_lp::StdCache`). Invalidated
    /// by `template_gen` bumps.
    std_cache: rusolve_lp::StdCache,
    /// Bumped by EVERY structural template edit (cut rows added/removed).
    /// Same-gen calls promise `template` differs only in bound values —
    /// the `StdCache` contract.
    template_gen: u64,
    /// Audit: time every solve into `stats.lp_nanos` (off by default).
    time_lps: bool,
}

impl RusolveEngine {
    pub fn new(p: &MipProblem, config: SimplexConfig) -> Self {
        // Each original row maps to 0, 1 or 2 row-form rows.
        let mut senses = Vec::new();
        let mut rhs = Vec::new();
        let mut targets: Vec<Vec<usize>> = Vec::with_capacity(p.num_row);
        let mut row_origin = Vec::new();
        for i in 0..p.num_row {
            let (lo, hi) = (p.row_lower[i], p.row_upper[i]);
            let mut t = Vec::new();
            if lo == hi {
                t.push(senses.len());
                senses.push('E');
                rhs.push(lo);
            } else {
                if hi.is_finite() {
                    t.push(senses.len());
                    senses.push('L');
                    rhs.push(hi);
                }
                if lo.is_finite() {
                    t.push(senses.len());
                    senses.push('G');
                    rhs.push(lo);
                }
            }
            row_origin.resize(row_origin.len() + t.len(), i);
            targets.push(t);
        }
        let num_row = senses.len();
        // CSC build: entries replicate into each mapped target row; per column
        // the source rows ascend and targets[r] = [hi, lo] ascends within and
        // across rows, so the pushed column is strictly ascending like the
        // dense scatter (values are bit-copies; CscBuilder drops zeros, which
        // the MipProblem matrix has none of).
        let mut b = rusolve_core::CscBuilder::new(num_row);
        let mut col: Vec<(usize, f64)> = Vec::new();
        let (ptr, idx, val) = (p.matrix.col_ptr(), p.matrix.row_idx(), p.matrix.values());
        for j in 0..p.num_col {
            col.clear();
            for k in ptr[j]..ptr[j + 1] {
                for &r in &targets[idx[k] as usize] {
                    col.push((r, val[k]));
                }
            }
            b.push_col(&col);
        }
        Self {
            template: Lp {
                num_col: p.num_col,
                num_row,
                a: b.into_matrix(),
                rhs,
                objective: p.cost.clone(),
                maximize: false,
                senses,
                lower: p.col_lower.clone(),
                upper: p.col_upper.clone(),
            },
            config,
            capped: config,
            root_dual: false,
            root_retry: false,
            solver: SimplexSolver::new(config),
            num_model_rows: num_row,
            row_origin,
            num_orig_rows: p.num_row,
            cut_ids: Rc::from(Vec::new()),
            next_cut_id: 0,
            last: None,
            warm: true,
            warm_dse: false,
            unfit_last: false,
            max_cold_iterations: 0,
            lp_deadline: None,
            stats: LpStats::default(),
            time_lps: false,
            lp_duals: None,
            last_ray: None,
            last_cutoff: None,
            std_cache: rusolve_lp::StdCache::new(),
            template_gen: 0,
        }
    }

    /// Turns warm starts off (every solve cold, as before warm starts) or on.
    pub fn set_warm(&mut self, on: bool) {
        self.warm = on;
    }

    /// Dual steepest edge on warm solves, with the weights kept in the
    /// snapshot (HiGHS keeps them with the basis, HEkk::putIterate).
    /// Audit lever: accumulate solve wall time in `LpStats::lp_nanos`.
    pub fn set_time_lps(&mut self, on: bool) {
        self.time_lps = on;
    }

    pub fn set_warm_dse(&mut self, on: bool) {
        self.warm_dse = on;
    }

    /// The cold solver for `capped` under the root policy. Warm solves
    /// keep `capped` itself: the policy is about cold starts only.
    fn rebuild_solver(&mut self) {
        let mut cfg = self.capped;
        cfg.dual |= self.root_dual;
        let mut solver = SimplexSolver::new(cfg);
        // The cap swap rebuilds the solver; the deadline is engine
        // state and must survive it (the tree swaps the cap around
        // every node LP).
        solver.set_deadline(self.lp_deadline);
        self.solver = solver;
    }

    /// Lever: a snapshot that does
    /// not fit the current rows (node cut separation dropped a row that was
    /// tight in it) warm-starts from the engine's latest optimal basis
    /// instead of the cold primal. HiGHS does the same: a node basis whose
    /// row count differs from the LP's is not installed and the LP keeps its
    /// current basis (HighsSearch.cpp:788–793, MIT, ERGO-Code/HiGHS); the
    /// dual simplex then repairs the bound changes from there.
    /// The same lever lets the warm entry flip boxed dual infeasibilities
    /// ([`SimplexConfig::warm_flip_boxed`]): `last` was optimal under other
    /// bounds, and HiGHS's dual simplex repairs exactly that by flipping.
    pub fn set_unfit_last(&mut self, on: bool) {
        self.unfit_last = on;
        self.config.warm_flip_boxed = on;
        self.capped.warm_flip_boxed = on;
    }

    /// Row-free LP: each column independently sits at its cheaper bound.
    fn solve_bounds_only(&self, col_lower: &[f64], col_upper: &[f64]) -> LpOutcome {
        let mut x = Vec::with_capacity(col_lower.len());
        for (j, &c) in self.template.objective.iter().enumerate() {
            let (lo, hi) = (col_lower[j], col_upper[j]);
            let v = if c > 0.0 {
                lo
            } else if c < 0.0 {
                hi
            } else if lo.is_finite() {
                lo
            } else if hi.is_finite() {
                hi
            } else {
                0.0
            };
            if !v.is_finite() {
                return LpOutcome {
                    status: LpStatus::Unbounded,
                    objective: f64::NEG_INFINITY,
                    x: vec![],
                    iterations: 0,
                    basis: None,
                };
            }
            x.push(v);
        }
        let objective = self
            .template
            .objective
            .iter()
            .zip(&x)
            .map(|(c, v)| c * v)
            .sum();
        LpOutcome {
            status: LpStatus::Optimal,
            objective,
            x,
            iterations: 0,
            basis: None,
        }
    }

    /// Wraps an optimal basis of the current LP and remembers it as `last`.
    fn remember(&mut self, basis: Option<Basis>, dse: Option<DseWeights>) -> Option<Rc<WarmStart>> {
        let ws = Rc::new(WarmStart {
            basis: basis?,
            dse,
            cut_ids: Rc::clone(&self.cut_ids),
        });
        self.last = Some(Rc::clone(&ws));
        Some(ws)
    }

    /// The cold solve of the current template (two-phase primal unless the
    /// configuration or the root policy picks the dual). `toggle` runs
    /// the other simplex instead, without the iteration cap (the deadline
    /// still applies).
    fn solve_cold(&mut self, toggle: bool) -> Result<LpOutcome, LpError> {
        let sol = if toggle {
            let mut solver = SimplexSolver::new(SimplexConfig {
                max_iterations: usize::MAX,
                dual: !(self.capped.dual || self.root_dual),
                ..self.capped
            });
            solver.set_deadline(self.lp_deadline);
            solver.solve(&self.template)
        } else {
            self.solver.solve(&self.template)
        }
        .map_err(LpError::Simplex)?;
        let iterations = (sol.phase1_iterations + sol.phase2_iterations) as u64;
        self.stats.cold_iterations += iterations;
        self.max_cold_iterations = self.max_cold_iterations.max(iterations);
        let status = match sol.status {
            SimplexStatus::Optimal => LpStatus::Optimal,
            SimplexStatus::Infeasible => LpStatus::Infeasible,
            SimplexStatus::Unbounded => LpStatus::Unbounded,
            SimplexStatus::IterationLimit => return Err(LpError::IterationLimit),
        };
        let optimal = status == LpStatus::Optimal;
        if status == LpStatus::Infeasible {
            self.last_ray = sol.dual_ray.clone();
        }
        let basis = if optimal {
            self.remember(sol.basis, None)
        } else {
            None
        };
        if let Some(ws) = &basis
            && sol.duals_valid
            && sol.row_duals.len() == self.template.num_row
        {
            self.lp_duals = Some((Rc::clone(ws), sol.row_duals));
        }
        Ok(LpOutcome {
            status,
            objective: sol.objective,
            x: if optimal { sol.x } else { vec![] },
            iterations,
            basis,
        })
    }

    /// `ws` re-expressed over the current rows and bounds, or `None` when it
    /// cannot be. Cut rows added since the snapshot enter with their slack
    /// basic (HiGHS `addCuts`); removed cut rows must have had a basic slack,
    /// as HiGHS `removeObsoleteRows` only removes those (dropping a tight row
    /// leaves one basic variable too many, and which one must leave is a
    /// tableau question); nonbasic columns move to a side their current
    /// bounds have (dual feasibility is then checked by the LP's warm entry).
    fn adapt(&self, ws: &WarmStart) -> Option<Basis> {
        let m0 = self.num_model_rows;
        let snap = &ws.basis;
        if snap.col_status.len() != self.template.num_col
            || snap.row_status.len() != m0 + ws.cut_ids.len()
        {
            return None;
        }
        let row_status = if Rc::ptr_eq(&ws.cut_ids, &self.cut_ids) || ws.cut_ids == self.cut_ids {
            snap.row_status.clone()
        } else {
            let old = &ws.cut_ids;
            let mut rows = Vec::with_capacity(self.template.num_row);
            rows.extend_from_slice(&snap.row_status[..m0]);
            let mut k = 0;
            for &id in self.cut_ids.iter() {
                while k < old.len() && old[k] < id {
                    if snap.row_status[m0 + k] != BasisStatus::Basic {
                        return None;
                    }
                    k += 1;
                }
                if k < old.len() && old[k] == id {
                    rows.push(snap.row_status[m0 + k]);
                    k += 1;
                } else {
                    rows.push(BasisStatus::Basic);
                }
            }
            if snap.row_status[m0 + k..]
                .iter()
                .any(|&s| s != BasisStatus::Basic)
            {
                return None;
            }
            rows
        };
        let t = &self.template;
        let mut col_status = snap.col_status.clone();
        for (j, s) in col_status.iter_mut().enumerate() {
            let (lo, hi) = (t.lower[j], t.upper[j]);
            *s = match *s {
                BasisStatus::Basic => BasisStatus::Basic,
                BasisStatus::Fixed if lo == hi => BasisStatus::Fixed,
                BasisStatus::AtUpper if hi.is_finite() => BasisStatus::AtUpper,
                BasisStatus::AtLower | BasisStatus::Fixed if lo.is_finite() => BasisStatus::AtLower,
                // A free column sits nonbasic at 0, exported as AtLower.
                BasisStatus::AtLower if hi == f64::INFINITY => BasisStatus::AtLower,
                _ if hi.is_finite() => BasisStatus::AtUpper,
                _ if lo.is_finite() => BasisStatus::AtLower,
                _ => return None,
            };
        }
        Some(Basis {
            col_status,
            row_status,
        })
    }

    /// The snapshot's DSE weights over the current rows, in the order
    /// [`DseWeights`] documents (basic columns, then basic rows). Cut rows
    /// added since the snapshot enter with a basic slack and weight 1 (the
    /// old rows' weights are unchanged by a new basic-slack row: the inverse
    /// gains a zero column, so `||e_i^T B^-1||` is the same); removed cut rows
    /// were basic (`adapt` checks) and just drop their entry. `None` when the
    /// snapshot has no weights or they do not fit it.
    fn adapt_dse(&self, ws: &WarmStart) -> Option<DseWeights> {
        let w = ws.dse.as_ref()?;
        let m0 = self.num_model_rows;
        let snap = &ws.basis;
        let nb_cols = snap
            .col_status
            .iter()
            .filter(|&&s| s == BasisStatus::Basic)
            .count();
        let nb_rows = snap
            .row_status
            .iter()
            .filter(|&&s| s == BasisStatus::Basic)
            .count();
        if w.basic.len() != nb_cols + nb_rows || snap.row_status.len() != m0 + ws.cut_ids.len() {
            return None;
        }
        if Rc::ptr_eq(&ws.cut_ids, &self.cut_ids) || ws.cut_ids == self.cut_ids {
            return Some(w.clone());
        }
        // Weight of each old row (None = nonbasic).
        let mut it = w.basic[nb_cols..].iter();
        let old_rows: Vec<Option<f32>> = snap
            .row_status
            .iter()
            .map(|&s| (s == BasisStatus::Basic).then(|| *it.next().unwrap_or(&1.0)))
            .collect();
        let mut basic = w.basic[..nb_cols].to_vec();
        basic.extend(old_rows[..m0].iter().flatten());
        let old = &ws.cut_ids;
        let mut k = 0;
        for &id in self.cut_ids.iter() {
            while k < old.len() && old[k] < id {
                k += 1; // removed cut row (basic slack): entry dropped
            }
            if k < old.len() && old[k] == id {
                basic.extend(old_rows[m0 + k]);
                k += 1;
            } else {
                basic.push(1.0);
            }
        }
        Some(DseWeights { basic })
    }

    /// Fallback for [`LpEngine::duals`] when the LP returned no valid duals
    /// (the LP gives them for every Optimal reopt and cold solve with presolve
    /// off). Row duals of an optimal basis of the current template: `y_i = 0` on
    /// rows whose slack is basic, and on the tight rows `T` (nonbasic slack,
    /// and every `E` row) `B_Tᵀ y_T = c_B`, where `B_T` is the matrix of the
    /// basic columns restricted to `T` — the textbook `Bᵀy = c_B` with the
    /// unit slack columns eliminated (Chvátal, *Linear Programming*, 1983,
    /// ch. 8). Redundant `E` rows whose artificial stays basic (which the
    /// snapshot cannot show) are dropped via [`Self::independent_rows`].
    /// `None` if `|T|` cannot match the number of basic columns, the system
    /// is singular, or its residual is not small.
    fn basis_row_duals(&self, basis: &Basis) -> Option<Vec<f64>> {
        let t = &self.template;
        let mut tight: Vec<usize> = (0..t.num_row)
            .filter(|&i| t.senses[i] == 'E' || basis.row_status[i] != BasisStatus::Basic)
            .collect();
        let basic: Vec<usize> = (0..t.num_col)
            .filter(|&j| basis.col_status[j] == BasisStatus::Basic)
            .collect();
        if tight.len() > basic.len() {
            // E-row status is a placeholder, so the basic artificials of
            // redundant E rows (flow conservation) are unknown: keep a
            // maximal independent set of E rows; the dropped ones get y = 0,
            // which is their dual when their artificial is basic.
            tight = Self::independent_rows(t, &tight, &basic)?;
        }
        if tight.len() != basic.len() {
            return None;
        }
        let mut y = vec![0.0; t.num_row];
        if basic.is_empty() {
            return Some(y);
        }
        let columns: Vec<Vec<f64>> = basic
            .iter()
            .map(|&j| tight.iter().map(|&i| t.a.get(i, j)).collect())
            .collect();
        let cb: Vec<f64> = basic.iter().map(|&j| t.objective[j]).collect();
        // Sparse Markowitz LU (the dual engine's own factor; 1e-11 is its
        // absolute singularity floor); the residual check guards accuracy.
        let lu = SparseLu::factor(&columns, 1e-11).ok()?;
        let mut yt = vec![0.0; tight.len()];
        lu.solve_transpose(&cb, &mut yt);
        // One step of iterative refinement, then a residual check.
        let residual = |yt: &[f64]| -> Vec<f64> {
            columns
                .iter()
                .zip(&cb)
                .map(|(col, &c)| c - col.iter().zip(yt).map(|(a, y)| a * y).sum::<f64>())
                .collect()
        };
        let r = residual(&yt);
        let mut dy = vec![0.0; tight.len()];
        lu.solve_transpose(&r, &mut dy);
        for (y, d) in yt.iter_mut().zip(&dy) {
            *y += d;
        }
        let worst = residual(&yt)
            .iter()
            .zip(&cb)
            .map(|(r, c)| r.abs() / (1.0 + c.abs()))
            .fold(0.0, f64::max);
        if worst.is_nan() || worst > 1e-9 || yt.iter().any(|v| !v.is_finite()) {
            return None;
        }
        for (&i, &v) in tight.iter().zip(&yt) {
            y[i] = v;
        }
        Some(y)
    }

    /// `tight` (non-E rows first in row order, then E rows) restricted to
    /// `basic.len()` rows whose restriction to the basic columns is linearly
    /// independent: every non-E row is kept, E rows enter greedily
    /// (incremental Gaussian elimination). `None` if a non-E row is
    /// dependent or too few rows remain.
    fn independent_rows(t: &Lp, tight: &[usize], basic: &[usize]) -> Option<Vec<usize>> {
        let k = basic.len();
        let order = tight
            .iter()
            .filter(|&&i| t.senses[i] != 'E')
            .chain(tight.iter().filter(|&&i| t.senses[i] == 'E'));
        let mut kept = Vec::with_capacity(k);
        let mut echelon: Vec<(usize, Vec<f64>)> = Vec::with_capacity(k);
        for &i in order {
            if kept.len() == k {
                break;
            }
            let mut row: Vec<f64> = basic.iter().map(|&j| t.a.get(i, j)).collect();
            let scale = row.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            for (p, e) in &echelon {
                let f = row[*p] / e[*p];
                if f != 0.0 {
                    for (r, v) in row.iter_mut().zip(e) {
                        *r -= f * v;
                    }
                }
            }
            let (p, piv) = row.iter().enumerate().fold((0, 0.0f64), |b, (c, v)| {
                if v.abs() > b.1 { (c, v.abs()) } else { b }
            });
            if piv > 1e-9 * scale && scale > 0.0 {
                kept.push(i);
                echelon.push((p, row));
            } else if t.senses[i] != 'E' {
                return None;
            }
        }
        kept.sort_unstable();
        (kept.len() == k).then_some(kept)
    }

    /// Largest row or bound violation of `x` on the current template,
    /// relative to `max(1, |rhs|)` / `max(1, |bound|)`.
    fn violation(&self, x: &[f64]) -> f64 {
        let t = &self.template;
        if x.len() != t.num_col || x.iter().any(|v| !v.is_finite()) {
            return f64::INFINITY;
        }
        let mut worst: f64 = 0.0;
        for (j, &v) in x.iter().enumerate() {
            let below = (t.lower[j] - v) / t.lower[j].abs().max(1.0);
            let above = (v - t.upper[j]) / t.upper[j].abs().max(1.0);
            worst = worst.max(below).max(above);
        }
        let mut act = vec![0.0; t.num_row];
        // A·x straight off the CSC (the core matvec: ascending-j accumulation
        // per row, exactly the dense loop's order — zeros are ±0.0 no-ops).
        t.a.matvec(x, &mut act);
        for ((&a, &r), &sense) in act.iter().zip(&t.rhs).zip(&t.senses) {
            let d = match sense {
                'L' => a - r,
                'G' => r - a,
                _ => (a - r).abs(),
            };
            worst = worst.max(d / r.abs().max(1.0));
        }
        worst
    }

    /// Warm dual reoptimization from `ws`. `Ok(None)` = solve cold instead
    /// (the reason is already counted in `stats`).
    ///
    /// The iteration limit is twice the largest cold solve seen (at least
    /// [`MIN_WARM_ITERATIONS`]): a warm run that needs more is stalling, and
    /// the cold primal is the known-good path. A probe (`probe_limit`) runs
    /// at most that many iterations and reports the cap instead of falling
    /// back.
    // TODO(perf): the limit is a guess; re-tune against the fallback counts
    // once node LPs are sparse.
    fn solve_warm(
        &mut self,
        ws: &WarmStart,
        cutoff: f64,
        probe_limit: Option<u64>,
    ) -> Result<Option<LpOutcome>, LpError> {
        // Holds `last` when it stands in for an unfit `ws` (lever).
        let recovered: Option<Rc<WarmStart>>;
        let (basis, ws) = match self.adapt(ws) {
            Some(b) => (b, ws),
            None => {
                recovered = self.last.clone().filter(|_| self.unfit_last);
                match recovered.as_deref().and_then(|l| Some((self.adapt(l)?, l))) {
                    Some(found) => {
                        self.stats.unfit_recovered += 1;
                        found
                    }
                    None => {
                        self.stats.fallback_unfit += 1;
                        return Ok(None);
                    }
                }
            }
        };
        // `reoptimize` prunes at bound ≥ cutoff; the tree prunes at
        // bound > cutoff, hence the next float up.
        let cut = cutoff.is_finite().then(|| cutoff.next_up());
        let limit = probe_limit
            .unwrap_or((2 * self.max_cold_iterations).max(MIN_WARM_ITERATIONS))
            .try_into()
            .unwrap_or(usize::MAX);
        let carried = if self.warm_dse {
            self.adapt_dse(ws)
        } else {
            None
        };
        let mode = match (&carried, self.warm_dse) {
            (Some(w), _) => DseMode::Carry(w),
            (None, true) => DseMode::Fresh,
            (None, false) => DseMode::Off,
        };
        // The search deadline reaches warm solves too, not only cold
        // solves.
        let Ok(r) = rusolve_lp::reoptimize_dse_cached_until(
            &self.template,
            &self.capped,
            &basis,
            cut,
            Some(limit),
            mode,
            &mut self.std_cache,
            self.template_gen,
            self.lp_deadline,
        ) else {
            self.stats.fallback_error += 1;
            return Ok(None);
        };
        let iterations = r.iters as u64;
        if r.warm {
            self.stats.warm_iterations += iterations;
        } else {
            // The LP's warm entry rejected the snapshot (or its dual engine
            // hit distress) and it solved cold itself: that answer is the
            // cold primal's and is used as is.
            self.stats.fallback_rejected += 1;
            self.stats.cold_iterations += iterations;
        }
        let done = |status, objective, x, basis| {
            Ok(Some(LpOutcome {
                status,
                objective,
                x,
                iterations,
                basis,
            }))
        };
        match r.status {
            ReoptStatus::Optimal => {
                if r.warm && self.violation(&r.x) > WARM_POINT_TOL {
                    self.stats.fallback_error += 1;
                    return Ok(None);
                }
                if r.warm {
                    self.stats.warm += 1;
                }
                let basis = self.remember(r.basis, r.dse_weights);
                // The LP's own duals, in the template's space; kept the
                // way a cold solve's are so `duals()` needs no LU.
                if let Some(ws) = &basis
                    && r.duals_valid
                    && r.row_duals.len() == self.template.num_row
                {
                    self.lp_duals = Some((Rc::clone(ws), r.row_duals));
                }
                done(LpStatus::Optimal, r.objective, r.x, basis)
            }
            ReoptStatus::Pruned => {
                if r.warm {
                    self.stats.warm += 1;
                    self.stats.cutoffs += 1;
                }
                if r.prune_duals.len() == self.template.num_row && cutoff.is_finite() {
                    self.last_cutoff = Some((r.prune_duals, cutoff));
                }
                done(LpStatus::Cutoff, r.obj_bound, vec![], None)
            }
            ReoptStatus::Infeasible if !r.warm => {
                self.last_ray = r.dual_ray;
                done(LpStatus::Infeasible, f64::INFINITY, vec![], None)
            }
            ReoptStatus::Infeasible if r.dual_ray.is_some() => {
                // The LP verified this Farkas ray against the original
                // rows (`yᵀrhs > sup yᵀ(Ax+s)` over the bounds), so the
                // verdict is a proof and needs no cold re-check.
                self.stats.warm += 1;
                self.stats.infeasible_trusted += 1;
                self.last_ray = r.dual_ray;
                done(LpStatus::Infeasible, f64::INFINITY, vec![], None)
            }
            ReoptStatus::Infeasible => {
                // A warm Infeasible without a ray prunes the node only once
                // the cold primal (phase 1, the verdict the tree has always
                // trusted) agrees.
                self.stats.warm += 1;
                self.stats.infeasible_checks += 1;
                Ok(None)
            }
            ReoptStatus::Unbounded if !r.warm => {
                done(LpStatus::Unbounded, f64::NEG_INFINITY, vec![], None)
            }
            ReoptStatus::Unbounded => {
                self.stats.fallback_error += 1;
                Ok(None)
            }
            ReoptStatus::IterationLimit if !r.warm => Err(LpError::IterationLimit),
            ReoptStatus::IterationLimit if probe_limit.is_some() => {
                self.stats.warm += 1;
                self.stats.probe_iter_limits += 1;
                let bound = if r.obj_bound.is_nan() {
                    f64::NEG_INFINITY
                } else {
                    r.obj_bound
                };
                done(LpStatus::IterationLimit, bound, vec![], None)
            }
            ReoptStatus::IterationLimit => {
                self.stats.fallback_iter_limit += 1;
                Ok(None)
            }
        }
    }
}

impl RusolveEngine {
    /// [`LpEngine::solve_from`] and [`LpEngine::probe`] (`probe_limit`).
    fn solve_capped(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        start: Option<&WarmStart>,
        cutoff: f64,
        probe_limit: Option<u64>,
    ) -> Result<LpOutcome, LpError> {
        if self.template.num_row == 0 {
            return Ok(self.solve_bounds_only(col_lower, col_upper));
        }
        self.template.lower.copy_from_slice(col_lower);
        self.template.upper.copy_from_slice(col_upper);
        self.last_ray = None;
        self.last_cutoff = None;
        let checks = self.stats.infeasible_checks;
        let mut spent = 0;
        match start.filter(|_| self.warm) {
            Some(ws) => {
                let before = self.stats.warm_iterations + self.stats.cold_iterations;
                if let Some(out) = self.solve_warm(ws, cutoff, probe_limit)? {
                    return Ok(out);
                }
                spent = self.stats.warm_iterations + self.stats.cold_iterations - before;
            }
            None => self.stats.cold += 1,
        }
        // A start without any snapshot (the root LP, a restart round's
        // root) that fails on the cap or gives up numerically re-runs once on
        // the other simplex, uncapped, while time is left (HiGHS's strategy
        // toggle, HighsLpRelaxation.cpp:1236-1251, and uncapped root LP,
        // HighsMipSolverData.cpp:2018). A warm run's cold fallback and a
        // probe keep the cap and the single attempt.
        let retry = self.root_retry && start.is_none() && probe_limit.is_none();
        let mut out = match self.solve_cold(false) {
            Err(LpError::IterationLimit | LpError::Simplex(_))
                if retry
                    && self
                        .lp_deadline
                        .is_some_and(|d| std::time::Instant::now() < d) =>
            {
                self.stats.root_retries += 1;
                self.solve_cold(true)?
            }
            r => r?,
        };
        if self.stats.infeasible_checks > checks && out.status != LpStatus::Infeasible {
            self.stats.infeasible_overturned += 1;
        }
        out.iterations += spent;
        Ok(out)
    }
}

impl RusolveEngine {
    /// Row weights of the tableau rows of the wanted basic columns, one
    /// `binv_row` BTRAN per candidate — the shape HiGHS uses, one
    /// `getBasisInverseRowSparse` per candidate (`lp_data/Highs.cpp:2272`,
    /// `mip/HighsTableauSeparator.cpp:144`). The session's cached warm state
    /// costs one factorization per call and its column→position map tracks
    /// the basis with it, so no position is ever scanned. `y_i` is the
    /// weight of LP row `i` directly: a row's slack/artificial column is
    /// ±`e_i`, so the old `tableau_row` slack slots held `±y_i` bit-exactly
    /// and the slot sign flips (`rusolve-lp dual.rs` `to_dual_standard`; the
    /// slack slots follow `StdLayout`) were returning `y_i`. A model row
    /// split into two LP rows gets both weights summed (at most one of the
    /// two slacks is nonbasic). Nonbasic candidates are omitted; `cols` is
    /// expected distinct (the separator builds it from distinct columns).
    fn candidate_weights(
        &self,
        session: &mut ReoptSession,
        cols: &[usize],
    ) -> Option<Vec<TableauAgg>> {
        if !session.warm_available() {
            return None;
        }
        let mut out = Vec::with_capacity(cols.len());
        for &j in cols {
            let Some(pos) = session.basic_position(j) else {
                continue;
            };
            let y = session.binv_row(pos)?;
            let mut weights: Vec<(usize, f64)> = Vec::new();
            for (i, &yi) in y.iter().enumerate() {
                if yi == 0.0 {
                    continue;
                }
                let r = if i < self.num_model_rows {
                    self.row_origin[i]
                } else {
                    self.num_orig_rows + (i - self.num_model_rows)
                };
                match weights.last_mut() {
                    Some(last) if last.0 == r => last.1 += yi,
                    _ => weights.push((r, yi)),
                }
            }
            out.push(TableauAgg { col: j, weights });
        }
        Some(out)
    }
}

impl LpEngine for RusolveEngine {
    fn solve(&mut self, col_lower: &[f64], col_upper: &[f64]) -> Result<LpOutcome, LpError> {
        let start = self.last.clone();
        self.solve_from(col_lower, col_upper, start.as_deref(), f64::INFINITY)
    }

    fn solve_from(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        start: Option<&WarmStart>,
        cutoff: f64,
    ) -> Result<LpOutcome, LpError> {
        let t0 = self.time_lps.then(std::time::Instant::now);
        let out = self.solve_capped(col_lower, col_upper, start, cutoff, None);
        if let Some(t0) = t0 {
            self.stats.lp_nanos += t0.elapsed().as_nanos() as u64;
        }
        out
    }

    /// The node's own basis stays the engine's `last` afterwards: probe
    /// bases are thrown away, as HiGHS's `Playground` restores the node
    /// iterate (HighsLpRelaxation.h:146).
    fn probe(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        start: Option<&WarmStart>,
        cutoff: f64,
        iter_limit: u64,
    ) -> Result<LpOutcome, LpError> {
        let last = self.last.clone();
        let t0 = self.time_lps.then(std::time::Instant::now);
        let out = self.solve_capped(col_lower, col_upper, start, cutoff, Some(iter_limit));
        if let Some(t0) = t0 {
            let ns = t0.elapsed().as_nanos() as u64;
            self.stats.lp_nanos += ns;
            self.stats.probe_nanos += ns;
        }
        self.last = last;
        if let Ok(o) = &out {
            self.stats.probes += 1;
            self.stats.probe_iterations += o.iterations;
        }
        out
    }

    fn set_iteration_limit(&mut self, limit: Option<u64>) {
        let max_iterations = limit.map_or(self.config.max_iterations, |l| {
            usize::try_from(l)
                .unwrap_or(usize::MAX)
                .min(self.config.max_iterations)
        });
        if max_iterations != self.capped.max_iterations {
            self.capped = SimplexConfig {
                max_iterations,
                ..self.config
            };
            self.rebuild_solver();
        }
    }

    fn set_root_lp(&mut self, dual: bool, retry: bool) {
        self.root_retry = retry;
        if dual != self.root_dual {
            self.root_dual = dual;
            self.rebuild_solver();
        }
    }

    fn set_deadline(&mut self, deadline: Option<std::time::Instant>) {
        self.lp_deadline = deadline;
        self.solver.set_deadline(deadline);
    }

    fn add_cuts(&mut self, cuts: &[CutRow]) -> bool {
        if cuts.is_empty() {
            return true;
        }
        let t = &mut self.template;
        let mut ids = self.cut_ids.to_vec();
        // Absorb the CSC, append one row per cut at the tail (cut
        // entries sorted by column — the push_row ascending contract; zero
        // entries stay structural exactly as the dense push(0.0) version left
        // them).
        let mut b = rusolve_core::CscBuilder::from_matrix(&t.a);
        for cut in cuts {
            let mut entries: Vec<(usize, f64)> = cut
                .inds
                .iter()
                .copied()
                .zip(cut.vals.iter().copied())
                .collect();
            entries.sort_unstable_by_key(|&(j, _)| j);
            b.push_row(&entries);
            t.senses.push('L');
            t.rhs.push(cut.rhs);
            t.num_row += 1;
            ids.push(self.next_cut_id);
            self.next_cut_id += 1;
        }
        t.a = b.into_matrix();
        self.cut_ids = Rc::from(ids);
        self.template_gen += 1; // rows appended: StdCache contract
        true
    }

    fn retain_cuts(&mut self, keep: &[bool]) {
        let m0 = self.num_model_rows;
        let t = &mut self.template;
        debug_assert_eq!(keep.len(), t.num_row - m0);
        let keep_row = |i: usize| i < m0 || keep[i - m0];
        // Row compaction map: model rows keep their slots, kept cuts re-pack.
        let mut map = vec![rusolve_core::NO_ROW; t.num_row];
        let mut next = 0usize;
        for (i, slot) in map.iter_mut().enumerate() {
            if keep_row(i) {
                *slot = next;
                next += 1;
            }
        }
        // Matrix filter, one CSC pass: surviving entries move bit-exact and
        // the ascending-row layout survives any subsequence (same recipe as
        // ReoptSession::remove_rows).
        let (optr, oridx, oval) = (t.a.col_ptr(), t.a.row_idx(), t.a.values());
        let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(oval.len());
        for j in 0..t.num_col {
            for k in optr[j]..optr[j + 1] {
                let i = oridx[k] as usize;
                if map[i] != rusolve_core::NO_ROW {
                    triplets.push((map[i], j, oval[k]));
                }
            }
        }
        t.a = rusolve_core::CscMatrix::from_triplets(next, t.num_col, &triplets);
        let mut i = 0;
        t.senses.retain(|_| {
            i += 1;
            keep_row(i - 1)
        });
        let mut i = 0;
        t.rhs.retain(|_| {
            i += 1;
            keep_row(i - 1)
        });
        t.num_row = t.rhs.len();
        if keep.iter().any(|&k| !k) {
            let ids: Vec<u64> = self
                .cut_ids
                .iter()
                .zip(keep)
                .filter_map(|(&id, &k)| k.then_some(id))
                .collect();
            self.cut_ids = Rc::from(ids);
            self.template_gen += 1; // rows dropped: StdCache contract
        }
    }

    /// Row weights of the tableau rows of `cols` (see
    /// [`LpEngine::tableau_weights`]). The session owns the template while it
    /// runs; the template goes back afterwards.
    fn tableau_weights(
        &mut self,
        lp: &LpOutcome,
        lower: &[f64],
        upper: &[f64],
        cols: &[usize],
    ) -> Option<Vec<TableauAgg>> {
        if lp.status != LpStatus::Optimal || self.template.num_row == 0 || cols.is_empty() {
            return None;
        }
        let ws = lp.basis.as_ref()?;
        self.template.lower.copy_from_slice(lower);
        self.template.upper.copy_from_slice(upper);
        // Per call the session costs one warm_setup factorization plus one
        // BTRAN per candidate. The old full-position scan refused
        // LPs above `m·(n+m) > 6e7`; that gate stays until relaxing it has
        // its own measured A/B.
        let t = &self.template;
        if t.num_row.saturating_mul(t.num_col + t.num_row) > TABLEAU_SCAN_LIMIT {
            return None;
        }
        let basis = self.adapt(ws)?;
        let empty = Lp {
            num_col: 0,
            num_row: 0,
            a: rusolve_core::CscMatrix::from_triplets(0, 0, &[]),
            rhs: Vec::new(),
            objective: Vec::new(),
            maximize: false,
            senses: Vec::new(),
            lower: Vec::new(),
            upper: Vec::new(),
        };
        let mut session = ReoptSession::new(
            std::mem::replace(&mut self.template, empty),
            &self.config,
            Some(basis),
        );
        let out = self.candidate_weights(&mut session, cols);
        std::mem::swap(session.lp_mut(), &mut self.template);
        out
    }

    fn num_lp_rows(&self) -> usize {
        self.template.num_row
    }

    fn stats(&self) -> LpStats {
        self.stats
    }

    /// `cᵀx = (c − Aᵀy)ᵀx + yᵀAx ≥ (c − Aᵀy)ᵀx + yᵀb` for every
    /// row-feasible `x` once `y` is clipped to `y ≤ 0` on `L` and `y ≥ 0` on
    /// `G` rows (wrong-signed entries to 0, as [`LpDuals::from_rows`]); with
    /// `cᵀx ≤ U` this gives `(c − Aᵀy)ᵀx ≤ U − yᵀb`, which the node box
    /// violates exactly when the Lagrangian bound exceeds `U`.
    fn cutoff_proof(&self) -> Option<(Vec<f64>, f64)> {
        let (y, u) = self.last_cutoff.as_ref()?;
        let t = &self.template;
        if y.len() != t.num_row || t.maximize {
            return None;
        }
        let y: Vec<f64> = y
            .iter()
            .zip(&t.senses)
            .map(|(&v, &s)| match s {
                'L' => v.min(0.0),
                'G' => v.max(0.0),
                _ => v,
            })
            .collect();
        // Aᵀy straight off the CSC (per column ascending rows — the dense
        // zip order; zero multipliers contribute ±0.0 no-ops, design §3.1).
        let mut coefs = vec![0.0; t.num_col];
        t.a.matvec_transpose(&y, &mut coefs);
        for (c, &obj) in coefs.iter_mut().zip(&t.objective) {
            *c = obj - *c;
        }
        let rhs = u - y.iter().zip(&t.rhs).map(|(yi, b)| yi * b).sum::<f64>();
        Some((coefs, rhs))
    }

    /// `yᵀAx ≥ yᵀb` holds for every row-feasible `x` when `y ≤ 0` on `L`
    /// rows and `y ≥ 0` on `G` rows (the slack terms `y_i s_i` are then
    /// non-positive), so `(−Aᵀy)ᵀx ≤ −yᵀb`; the ray's identity
    /// (`yᵀb > sup yᵀ(Ax+s)` over the bounds) says the bounds violate it.
    fn infeasibility_proof(&self) -> Option<(Vec<f64>, f64)> {
        let y = self.last_ray.as_ref()?;
        let t = &self.template;
        if y.len() != t.num_row {
            return None;
        }
        let mut coefs = vec![0.0; t.num_col];
        t.a.matvec_transpose(y, &mut coefs);
        for c in coefs.iter_mut() {
            *c = -*c;
        }
        let rhs = -y.iter().zip(&t.rhs).map(|(yi, b)| yi * b).sum::<f64>();
        Some((coefs, rhs))
    }

    /// The LP's own duals when `lp` was the latest solve over the current rows
    /// (cold, or warm with `duals_valid`); otherwise recomputed from `lp.basis` re-expressed over the current
    /// rows (cut rows added since enter with a basic slack, i.e. `y = 0`).
    fn duals(&mut self, lp: &LpOutcome) -> Option<LpDuals> {
        if lp.status != LpStatus::Optimal {
            return None;
        }
        if self.template.num_row == 0 {
            self.stats.duals_lp += 1;
            return Some(LpDuals::from_rows(&self.template, Vec::new()));
        }
        let Some(ws) = lp.basis.as_ref() else {
            self.stats.duals_failed += 1;
            return None;
        };
        if let Some((cws, y)) = &self.lp_duals
            && Rc::ptr_eq(cws, ws)
            && cws.cut_ids == self.cut_ids
        {
            self.stats.duals_lp += 1;
            return Some(LpDuals::from_rows(&self.template, y.clone()));
        }
        match self.adapt(ws).and_then(|b| self.basis_row_duals(&b)) {
            Some(y) => {
                self.stats.duals_basis += 1;
                Some(LpDuals::from_rows(&self.template, y))
            }
            None => {
                self.stats.duals_failed += 1;
                None
            }
        }
    }

    /// The snapshot basis of `lp`, re-expressed over the current rows by
    /// [`Self::adapt`] (which refuses when a removed row was nonbasic), read
    /// at the cut rows' slack status. No LU: statuses, not duals.
    fn cut_rows_basic(&mut self, lp: &LpOutcome) -> Option<Vec<bool>> {
        if lp.status != LpStatus::Optimal {
            return None;
        }
        let ws = lp.basis.as_ref()?;
        let basis = self.adapt(ws)?;
        let m0 = self.num_model_rows;
        Some(
            basis.row_status[m0..]
                .iter()
                .map(|s| *s == BasisStatus::Basic)
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_io::read_mps;

    /// The LP of a 3-node ladder proof: root optimum
    /// (0.9, 0, 2.4, 6), fractional.
    const LADDER: &str = "\
NAME LADDER
OBJSENSE
 MAX
ROWS
 N  OBJ
 L  R1
 L  R2
 L  R3
COLUMNS
    X0  OBJ  3  R1  1
    X0  R2  4  R3  2
    X1  OBJ  5  R1  1
    X1  R2  3  R3  5
    X2  OBJ  4  R1  2
    X2  R2  1  R3  3
    X3  OBJ  2  R1  1
    X3  R2  2  R3  1
RHS
    RHS  R1  12  R2  18
    RHS  R3  15
BOUNDS
 UP BND X0 6
 UP BND X1 6
 UP BND X2 6
 UP BND X3 6
ENDATA
";

    fn problem() -> MipProblem {
        MipProblem::from_model(&read_mps(LADDER).unwrap())
    }

    fn engine(p: &MipProblem, warm: bool) -> RusolveEngine {
        let mut e = RusolveEngine::new(
            p,
            SimplexConfig {
                presolve: false,
                ..SimplexConfig::default()
            },
        );
        e.set_warm(warm);
        e
    }

    /// Bounds of the two children of the first fractional column of `x`.
    fn children(p: &MipProblem, x: &[f64]) -> [(Vec<f64>, Vec<f64>); 2] {
        let j = (0..p.num_col)
            .find(|&j| (x[j] - x[j].round()).abs() > 1e-6)
            .expect("fractional root");
        let mut down = (p.col_lower.clone(), p.col_upper.clone());
        down.1[j] = x[j].floor();
        let mut up = (p.col_lower.clone(), p.col_upper.clone());
        up.0[j] = x[j].ceil();
        [down, up]
    }

    /// A warm Infeasible child with a Farkas ray prunes without the cold
    /// re-check.
    #[test]
    fn warm_infeasible_with_a_ray_prunes_without_a_cold_check() {
        let p = problem();
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        let ws = root.basis.clone().unwrap();
        // R1: x0+x1+2x2+x3 <= 12; forcing every column to its upper bound 6
        // makes the left side 30.
        let lo = vec![6.0; p.num_col];
        let hi = p.col_upper.clone();
        let out = e.solve_from(&lo, &hi, Some(&ws), f64::INFINITY).unwrap();
        assert_eq!(out.status, LpStatus::Infeasible);
        let s = e.stats();
        assert_eq!(s.infeasible_trusted, 1, "{s:?}");
        assert_eq!((s.infeasible_checks, s.cold), (0, 1), "{s:?}");
    }

    /// A hand-built infeasible node: its Farkas ray gives a proof row that
    /// the node box violates, the explanation keeps only the bounds the proof
    /// needs, those bounds alone are still infeasible for the LP, and the
    /// stored conflict prunes the node and tightens a relaxed sibling.
    #[test]
    fn farkas_proof_explains_and_cuts_off_the_node() {
        use crate::conflict::{
            ConflictPool, ConflictStats, Explain, Revisit, clean_proof, explain_leq,
        };
        let p = problem();
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        let ws = root.basis.clone().unwrap();
        // R1: x0+x1+2x2+x3 <= 12. The node raises x0, x1, x2 to 4 (left
        // side 16) and x3 to 1 (a bound the proof does not need).
        let lo = vec![4.0, 4.0, 4.0, 1.0];
        let hi = p.col_upper.clone();
        let out = e.solve_from(&lo, &hi, Some(&ws), f64::INFINITY).unwrap();
        assert_eq!(out.status, LpStatus::Infeasible);
        let (coefs, rhs) = e.infeasibility_proof().expect("trusted ray");
        // The proof is globally valid: every LP-feasible point satisfies it
        // (the root optimum does) and the node box does not.
        let lhs = |x: &[f64]| coefs.iter().zip(x).map(|(a, v)| a * v).sum::<f64>();
        assert!(lhs(&root.x) <= rhs + 1e-7, "{} > {rhs}", lhs(&root.x));
        let box_min: f64 = coefs
            .iter()
            .enumerate()
            .map(|(j, &a)| if a > 0.0 { a * lo[j] } else { a * hi[j] })
            .sum();
        assert!(box_min > rhs + 1e-9, "{box_min} <= {rhs}");
        let proof = clean_proof(&coefs, rhs, &p.col_lower, &p.col_upper, &p.integer, 1e-6)
            .expect("finite root box");
        let zero = |_: usize| 0;
        let ex = Explain {
            global_lower: &p.col_lower,
            global_upper: &p.col_upper,
            integer: &p.integer,
            feastol: 1e-6,
            up_count: &zero,
            down_count: &zero,
        };
        let lits = explain_leq(&proof, &lo, &hi, &ex).expect("explained");
        assert!(!lits.is_empty() && lits.len() <= 3, "{lits:?}");
        assert!(
            lits.iter().all(|l| l.col != 3),
            "x3 is not needed: {lits:?}"
        );
        // Sufficiency: the literals alone (all else at the root box) make
        // the LP infeasible.
        let (mut llo, lhi) = (p.col_lower.clone(), p.col_upper.clone());
        for l in &lits {
            assert!(!l.upper);
            llo[l.col as usize] = llo[l.col as usize].max(l.val);
        }
        let cold = e.solve_from(&llo, &lhi, None, f64::INFINITY).unwrap();
        assert_eq!(cold.status, LpStatus::Infeasible, "{lits:?}");
        // The pooled conflict prunes the node ...
        let mut pool = ConflictPool::new(p.num_col, p.integer.iter().filter(|&&i| i).count());
        assert!(pool.add(lits.clone(), &p.integer, 1e-6));
        let mut st = ConflictStats::default();
        let (mut nlo, mut nhi, mut changed) = (lo.clone(), hi.clone(), Vec::new());
        let r = pool.propagate(
            &mut nlo,
            &mut nhi,
            &p.integer,
            1e-6,
            Revisit::All,
            &mut changed,
            &mut st,
        );
        assert!(r.is_err() && st.pruned == 1);
        // ... and, with one of its bounds not yet imposed, forbids it.
        let l0 = lits[0];
        let mut slo = lo.clone();
        slo[l0.col as usize] = p.col_lower[l0.col as usize];
        let (mut shi, mut changed) = (hi.clone(), Vec::new());
        pool.propagate(
            &mut slo,
            &mut shi,
            &p.integer,
            1e-6,
            Revisit::All,
            &mut changed,
            &mut st,
        )
        .unwrap();
        assert_eq!(changed, vec![l0.col as usize]);
        // continuous bound: stored shifted by feastol, negated back to l0.val
        assert!(shi[l0.col as usize] <= l0.val + 1e-9);
        assert!(shi[l0.col as usize] < hi[l0.col as usize]);
    }

    #[test]
    fn children_reoptimize_warm_from_the_parent_basis() {
        let p = problem();
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        assert_eq!(root.status, LpStatus::Optimal);
        let ws = root.basis.clone().expect("root basis");
        let (mut warm_iters, mut cold_iters) = (0, 0);
        for (lo, hi) in children(&p, &root.x) {
            let w = e.solve_from(&lo, &hi, Some(&ws), f64::INFINITY).unwrap();
            let c = engine(&p, false).solve(&lo, &hi).unwrap();
            assert_eq!(w.status, c.status);
            assert!((w.objective - c.objective).abs() < 1e-9);
            assert!(w.basis.is_some(), "children hand a basis to grandchildren");
            warm_iters += w.iterations;
            cold_iters += c.iterations;
        }
        let s = e.stats();
        assert_eq!((s.cold, s.warm, s.fallbacks()), (1, 2, 0), "{s:?}");
        assert!(
            warm_iters < cold_iters,
            "warm {warm_iters} vs cold {cold_iters}"
        );
    }

    #[test]
    fn cutoff_prunes_with_a_valid_bound() {
        let p = problem();
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        let ws = root.basis.clone().unwrap();
        let [a, b] = children(&p, &root.x);
        let ta = engine(&p, false).solve(&a.0, &a.1).unwrap().objective;
        let tb = engine(&p, false).solve(&b.0, &b.1).unwrap().objective;
        let ((lo, hi), worse) = if ta > tb { (a, ta) } else { (b, tb) };
        let cutoff = worse - 1e-6;
        let out = e.solve_from(&lo, &hi, Some(&ws), cutoff).unwrap();
        assert_eq!(out.status, LpStatus::Cutoff);
        assert!(out.objective > cutoff && out.objective <= worse + 1e-9);
        // Above the child's optimum the cutoff must not fire.
        let out = e.solve_from(&lo, &hi, Some(&ws), worse + 1e-3).unwrap();
        assert_eq!(out.status, LpStatus::Optimal);
        assert!((out.objective - worse).abs() < 1e-9);
    }

    /// Whenever a warm solve ends `Cutoff` the engine hands out a proof
    /// row that the node box violates and that the root optimum (objective
    /// below every cutoff tried) satisfies, so it is globally valid for the
    /// improving solutions; it is gone after the next solve.
    #[test]
    fn cutoff_proof_is_violated_by_the_box_and_valid_at_the_root() {
        let p = problem();
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        let ws = root.basis.clone().unwrap();
        let mut proofs = 0;
        for (lo, hi) in children(&p, &root.x) {
            let exact = engine(&p, false).solve(&lo, &hi).unwrap().objective;
            for k in 1..20 {
                let cutoff = root.objective + (exact - root.objective) * (k as f64 / 20.0);
                let out = e.solve_from(&lo, &hi, Some(&ws), cutoff).unwrap();
                if out.status != LpStatus::Cutoff {
                    assert!(e.cutoff_proof().is_none());
                    continue;
                }
                let Some((coefs, rhs)) = e.cutoff_proof() else {
                    continue;
                };
                proofs += 1;
                let lhs = |x: &[f64]| coefs.iter().zip(x).map(|(a, v)| a * v).sum::<f64>();
                assert!(root.objective <= cutoff + 1e-9);
                assert!(lhs(&root.x) <= rhs + 1e-7, "{} > {rhs}", lhs(&root.x));
                let box_min: f64 = coefs
                    .iter()
                    .enumerate()
                    .map(|(j, &a)| if a > 0.0 { a * lo[j] } else { a * hi[j] })
                    .sum();
                assert!(box_min > rhs - 1e-7, "{box_min} vs {rhs}");
            }
        }
        assert!(proofs > 0, "no cutoff proof was ever produced");
        let _ = e.solve(&p.col_lower, &p.col_upper).unwrap();
        assert!(e.cutoff_proof().is_none());
    }

    #[test]
    fn capped_probe_reports_a_valid_bound_and_keeps_the_node_basis() {
        let p = problem();
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        let ws = root.basis.clone().unwrap();
        for (lo, hi) in children(&p, &root.x) {
            let exact = engine(&p, false).solve(&lo, &hi).unwrap();
            let full = e.probe(&lo, &hi, Some(&ws), f64::INFINITY, 1000).unwrap();
            assert_eq!(full.status, exact.status);
            assert!((full.objective - exact.objective).abs() < 1e-9);
            let capped = e.probe(&lo, &hi, Some(&ws), f64::INFINITY, 0).unwrap();
            match capped.status {
                LpStatus::IterationLimit => {
                    assert!(capped.x.is_empty());
                    assert!(capped.objective <= exact.objective + 1e-9);
                }
                // A child the parent basis already solves needs no pivot.
                LpStatus::Optimal => assert!((capped.objective - exact.objective).abs() < 1e-9),
                s => panic!("{s:?}"),
            }
            // Probes never replace the node basis `solve` starts from.
            assert!(Rc::ptr_eq(e.last.as_ref().unwrap(), &ws));
        }
        let s = e.stats();
        assert_eq!(s.probes, 4, "{s:?}");
        assert!(s.probe_iter_limits >= 1, "{s:?}");
    }

    /// Checks the textbook rule on a solved relaxation: for a minimisation,
    /// a column at its lower bound with `d_j > 0` costs at least `d_j` per
    /// unit it moves up (so `x_j ≤ l_j + (cutoff − z)/d_j` below `cutoff`),
    /// one at its upper bound with `d_j < 0` at least `−d_j` per unit down;
    /// and the Lagrangian bound equals the LP objective.
    fn check_textbook_rule(p: &MipProblem, lo: &[f64], hi: &[f64], lp: &LpOutcome, d: &LpDuals) {
        let z = lp.objective;
        let bound = d.bound(lo, hi).expect("finite box");
        assert!((bound - z).abs() < 1e-8, "Lagrangian {bound} vs LP {z}");
        let mut checked = 0;
        for j in 0..p.num_col {
            let dj = d.reduced_costs[j];
            for t in [1.0, 2.5] {
                let (mut l2, mut u2) = (lo.to_vec(), hi.to_vec());
                if dj > 1e-9 && lp.x[j] == lo[j] && lo[j] + t <= hi[j] {
                    l2[j] = lo[j] + t;
                } else if dj < -1e-9 && lp.x[j] == hi[j] && hi[j] - t >= lo[j] {
                    u2[j] = hi[j] - t;
                } else {
                    continue;
                }
                let moved = engine(p, false).solve(&l2, &u2).unwrap();
                if moved.status == LpStatus::Optimal {
                    assert!(
                        moved.objective >= z + dj.abs() * t - 1e-7,
                        "col {j} d {dj} t {t}: {} < {z} + |d|·t",
                        moved.objective
                    );
                }
                checked += 1;
            }
        }
        assert!(
            checked > 0,
            "no nonbasic column with a nonzero reduced cost"
        );
    }

    #[test]
    fn reduced_costs_follow_the_textbook_rule_cold_and_warm() {
        // LADDER maximises: internally it is min −c, so this also pins the
        // maximisation sign handling.
        let p = problem();
        assert!(p.maximize);
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        let cold = e.duals(&root).expect("cold duals");
        assert_eq!((e.stats().duals_lp, e.stats().duals_basis), (1, 0));
        // At the root optimum x1 = 0 sits at its lower bound and x3 = 6 at
        // its upper bound: d_1 > 0 and d_3 < 0 in the minimisation sense.
        assert!(root.x[1] == 0.0 && cold.reduced_costs[1] > 0.0, "{cold:?}");
        assert!(root.x[3] == 6.0 && cold.reduced_costs[3] < 0.0, "{cold:?}");
        check_textbook_rule(&p, &p.col_lower, &p.col_upper, &root, &cold);
        // The same LP re-solved warm (0 pivots): duals from the basis agree.
        let ws = root.basis.clone().unwrap();
        let again = e
            .solve_from(&p.col_lower, &p.col_upper, Some(&ws), f64::INFINITY)
            .unwrap();
        assert_eq!(e.stats().warm, 1);
        let warm = e.duals(&again).expect("warm duals");
        // The warm solve's own duals, no LU.
        assert_eq!((e.stats().duals_lp, e.stats().duals_basis), (2, 0));
        for (a, b) in warm.reduced_costs.iter().zip(&cold.reduced_costs) {
            assert!((a - b).abs() < 1e-9, "{warm:?} vs {cold:?}");
        }
        // The LU fallback agrees with the LP's duals on this basis.
        let lu = e
            .basis_row_duals(&e.adapt(&ws).unwrap())
            .expect("basis duals");
        let from_lu = LpDuals::from_rows(&e.template, lu);
        for (a, b) in from_lu.reduced_costs.iter().zip(&warm.reduced_costs) {
            assert!((a - b).abs() < 1e-9, "{from_lu:?} vs {warm:?}");
        }
        // Warm children: duals from their exported bases obey the rule.
        for (lo, hi) in children(&p, &root.x) {
            let w = e.solve_from(&lo, &hi, Some(&ws), f64::INFINITY).unwrap();
            if w.status != LpStatus::Optimal {
                continue;
            }
            let d = e.duals(&w).expect("child duals");
            check_textbook_rule(&p, &lo, &hi, &w, &d);
        }
        assert_eq!(e.stats().duals_failed, 0, "{:?}", e.stats());
        assert_eq!(e.stats().duals_basis, 0, "{:?}", e.stats());
    }

    /// Flow conservation on 3 nodes: the E rows sum to zero, so one
    /// artificial stays basic and the exported basis has fewer basic
    /// columns than E rows. Optimum a12 = a23 = 2, a13 = 0 with d = 3 − 2.
    const FLOW: &str = "\
NAME FLOW
ROWS
 N  OBJ
 E  N1
 E  N2
 E  N3
COLUMNS
    A12  OBJ  1  N1  1
    A12  N2  -1
    A23  OBJ  1  N2  1
    A23  N3  -1
    A13  OBJ  3  N1  1
    A13  N3  -1
RHS
    RHS  N1  2  N3  -2
BOUNDS
 UP BND A12 5
 UP BND A23 5
 UP BND A13 1
ENDATA
";

    #[test]
    fn basis_duals_survive_redundant_equality_rows() {
        let p = MipProblem::from_model(&read_mps(FLOW).unwrap());
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        assert!((root.objective - 4.0).abs() < 1e-9);
        let ws = root.basis.clone().unwrap();
        let again = e
            .solve_from(&p.col_lower, &p.col_upper, Some(&ws), f64::INFINITY)
            .unwrap();
        // The LP's own duals serve a warm solve with no LU ...
        let lp = e.duals(&again).expect("LP duals");
        assert_eq!(e.stats().duals_basis, 0, "{:?}", e.stats());
        // ... the LU fallback (used when `duals_valid` is false) still copes
        // with the basic artificial of a redundant E row.
        e.lp_duals = None;
        let d = e
            .duals(&again)
            .expect("basis duals with a basic artificial");
        for (a, b) in d.reduced_costs.iter().zip(&lp.reduced_costs) {
            assert!((a - b).abs() < 1e-9, "{d:?} vs {lp:?}");
        }
        assert_eq!(e.stats().duals_basis, 1, "{:?}", e.stats());
        assert!((d.reduced_costs[2] - 1.0).abs() < 1e-9, "{d:?}");
        check_textbook_rule(&p, &p.col_lower, &p.col_upper, &again, &d);
    }

    /// `tableau_weights` returns the weights of the candidate's tableau
    /// row, per candidate: aggregating the model rows with them annihilates
    /// every OTHER basic structural column and leaves the candidate's own
    /// column at coefficient 1 — the tableau's identity on its basic columns,
    /// position-order-free and computed here from the raw model matrix.
    /// Nonbasic candidates are omitted (the whole column list is asked for).
    #[test]
    fn tableau_weights_are_the_tableau_row_of_each_candidate() {
        let p = problem();
        let mut e = engine(&p, true);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        let ws = root.basis.clone().unwrap();
        let basic_cols: Vec<usize> = (0..p.num_col)
            .filter(|&j| ws.basis.col_status[j] == BasisStatus::Basic)
            .collect();
        assert!(!basic_cols.is_empty(), "test LP must have a basic column");
        let all: Vec<usize> = (0..p.num_col).collect();
        let aggs = e
            .tableau_weights(&root, &p.col_lower, &p.col_upper, &all)
            .expect("tableau access on the root basis");
        assert_eq!(aggs.len(), basic_cols.len(), "{aggs:?}");
        let acol = |j: usize| -> Vec<f64> {
            let mut c = vec![0.0; p.num_row];
            let (ptr, idx, val) = (p.matrix.col_ptr(), p.matrix.row_idx(), p.matrix.values());
            for k in ptr[j]..ptr[j + 1] {
                c[idx[k] as usize] = val[k];
            }
            c
        };
        for agg in &aggs {
            assert!(
                basic_cols.contains(&agg.col),
                "nonbasic candidate {:?}",
                agg.col
            );
            assert!(!agg.weights.is_empty());
            for &k in &basic_cols {
                let dot: f64 = agg.weights.iter().map(|&(r, w)| w * acol(k)[r]).sum();
                let want = f64::from(k == agg.col);
                assert!(
                    (dot - want).abs() <= 1e-9,
                    "candidate {}: coefficient of basic column {k} is {dot}, want {want}",
                    agg.col
                );
            }
        }
        // Repeated calls answer identically (the session is swapped back in
        // cleanly, the template survives the query).
        let again = e
            .tableau_weights(&root, &p.col_lower, &p.col_upper, &all)
            .expect("second call");
        assert_eq!(again, aggs);
    }

    #[test]
    fn cut_rows_enter_basic_and_leave_when_slack() {
        let p = problem();
        let (lo, hi) = (p.col_lower.clone(), p.col_upper.clone());
        let mut e = engine(&p, true);
        let mut c = engine(&p, false);
        let root = e.solve(&lo, &hi).unwrap();
        assert!((root.x[2] - 2.4).abs() < 1e-9);
        let tight = CutRow {
            inds: vec![2],
            vals: vec![1.0],
            rhs: 2.0,
        };
        let loose = CutRow {
            inds: vec![0, 1, 2, 3],
            vals: vec![1.0; 4],
            rhs: 100.0,
        };
        let same = |e: &mut RusolveEngine, c: &mut RusolveEngine| {
            let (w, k) = (e.solve(&lo, &hi).unwrap(), c.solve(&lo, &hi).unwrap());
            assert_eq!(w.status, LpStatus::Optimal);
            assert!((w.objective - k.objective).abs() < 1e-9);
        };
        for cut in [&tight, &loose] {
            e.add_cuts(std::slice::from_ref(cut));
            c.add_cuts(std::slice::from_ref(cut));
            same(&mut e, &mut c);
        }
        assert_eq!(e.stats().warm, 2, "{:?}", e.stats());
        // Dropping the loose (basic-slack) row keeps the basis.
        e.retain_cuts(&[true, false]);
        c.retain_cuts(&[true, false]);
        same(&mut e, &mut c);
        assert_eq!((e.stats().warm, e.stats().fallbacks()), (3, 0));
        // Dropping the tight row leaves one basic too many: cold, still right.
        e.retain_cuts(&[false]);
        c.retain_cuts(&[false]);
        same(&mut e, &mut c);
        assert_eq!(e.stats().fallback_unfit, 1, "{:?}", e.stats());
    }

    /// Lever: a snapshot left unfit by a dropped tight cut row
    /// warm-starts from the latest basis (HiGHS keeps its current basis,
    /// HighsSearch.cpp:788–793), at the parent bounds and at both children.
    #[test]
    fn unfit_snapshot_recovers_from_latest_basis() {
        let p = problem();
        let (lo, hi) = (p.col_lower.clone(), p.col_upper.clone());
        let mut e = engine(&p, true);
        e.set_unfit_last(true);
        let mut c = engine(&p, false);
        e.solve(&lo, &hi).unwrap();
        let tight = CutRow {
            inds: vec![2],
            vals: vec![1.0],
            rhs: 2.0,
        };
        e.add_cuts(std::slice::from_ref(&tight));
        c.add_cuts(std::slice::from_ref(&tight));
        let cut_lp = e.solve(&lo, &hi).unwrap();
        let stale = cut_lp.basis.clone().expect("optimal snapshot");
        e.retain_cuts(&[false]);
        c.retain_cuts(&[false]);
        // `last` is the stale snapshot itself: no fit anywhere, cold.
        let root = e.solve(&lo, &hi).unwrap();
        assert_eq!(e.stats().fallback_unfit, 1, "{:?}", e.stats());
        let mut bounds = vec![(lo.clone(), hi.clone())];
        bounds.extend(children(&p, &root.x));
        for (l, u) in &bounds {
            let w = e.solve_from(l, u, Some(&stale), f64::INFINITY).unwrap();
            let k = c.solve(l, u).unwrap();
            assert_eq!(w.status, k.status);
            if k.status == LpStatus::Optimal {
                assert!((w.objective - k.objective).abs() < 1e-9);
            }
        }
        let st = e.stats();
        assert_eq!(st.fallback_unfit, 1, "{st:?}");
        assert_eq!(st.unfit_recovered, bounds.len() as u64, "{st:?}");
    }
}
