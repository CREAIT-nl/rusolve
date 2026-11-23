//! Two-phase bounded-variable revised primal simplex over [`ProductFormFactorization`].
//!
//! Scope: continuous LPs. Finite lower bounds are shifted out in
//! preprocessing; upper bounds (finite or infinite) are handled natively via
//! bound flipping (no bound-slack rows). Free variables (`lower = -inf`) are
//! rejected as `Unsupported`.
//!
//! Algorithm:
//! * Preprocessing: shift lower bounds to 0 (rhs and reported objective
//!   adjusted); negate rows with negative rhs (L↔G flip); slack columns
//!   (L:+1, G:−1, E:none); one artificial per row; start on the artificial
//!   basis with values = rhs ≥ 0 (always primal feasible).
//! * Phase 1 minimizes the artificial sum; optimum > tol ⇒ Infeasible.
//!   Artificials are locked out of entering in phase 2 (basic ones stuck at 0
//!   are harmless; the row-dependency cleanup is a later nicety).
//! * Pricing: Dantzig — most negative reduced cost among nonbasic-at-lb, most
//!   positive among nonbasic-at-ub; ties by lowest index. After
//!   `bland_after` consecutive θ=0 pivots, Bland's rule (lowest index)
//!   engages — the anti-cycling guard — and is reported in the solution.
//! * Ratio test: full bounded-variable version. Basic values decreasing leave
//!   at lb; basic values increasing toward a finite ub flip out at ub; the
//!   entering variable itself may saturate at its ub and flip sides without a
//!   basis change. No blocker and no ub cap ⇒ Unbounded.
//! * Updates: single eta via [`ProductFormFactorization::update`], threshold
//!   refactorization inside PFI; on `RefactorizeNeeded` the factorization is
//!   rebuilt and the pivot retried once.
//!
//! Provenance: revised simplex with bounded variables per Chvátal, *Linear
//! Programming* (1983) — pricing, ratio test, two-phase scheme; Bland's rule:
//! Bland (1977); bound-flipping side bookkeeping follows the standard shape
//! described in HiGHS' simplex interface notes (MIT; consulted at boundary
//! level only, no code copied). Eta-file semantics: [`crate::pfi`]
//! (Forrest & Gold BASP; Huangfu & Hall arXiv:1503.01889).
//!
//! No unsafe. f64 throughout.
//!
//! `// TODO(perf)` markers name the planned upgrades: CSC columns, incremental
//! FTRAN/BTRAN pricing, steepest-edge, dense pricing
//! scan today is O(m·total) per iteration.

use crate::pfi::ProductFormFactorization;
use crate::{FactorError, PfiConfig};
use rusolve_core::{Lp, NO_ROW, SparseLp};

// `Lp` IS `rusolve_core::Lp` — CSC matrix + row
// form. The crate-local dense duplicate (`col: Vec<Vec<f64>>`) and its two
// field-copy `From` impls are gone; `rusolve_lp::Lp` re-exports the core type.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimplexStatus {
    Optimal,
    Infeasible,
    Unbounded,
    /// Iteration cap hit without termination (tune [`SimplexConfig`]).
    IterationLimit,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SimplexError {
    Unsupported(String),
    Factor(FactorError),
}

impl std::fmt::Display for SimplexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SimplexError::Unsupported(what) => write!(f, "unsupported: {what}"),
            SimplexError::Factor(e) => write!(f, "factorization failure: {e}"),
        }
    }
}

impl std::error::Error for SimplexError {}

impl From<FactorError> for SimplexError {
    fn from(e: FactorError) -> Self {
        SimplexError::Factor(e)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SimplexConfig {
    /// Pivot / ratio tolerance (factorization-side: blocking-row detection).
    pub pivot_tol: f64,
    /// Reduced-cost pricing tolerance. Deliberately looser than pivot_tol:
    /// on ill-conditioned degenerate bases (scsd1) BTRAN noise reaches ~1e-8,
    /// and noise-level "improvements" on columns with B⁻¹a ≤ 0 were being
    /// declared unbounded rays (false Infeasible in phase 1). Reduced costs
    /// below this are treated as zero ⇒ optimality loss ~tol·Σ|x| change,
    /// far under the 1e-5 bench agreement bar.
    pub pricing_tol: f64,
    /// Consecutive θ=0 pivots before a Bland episode engages.
    pub bland_after: usize,
    /// Iteration budget of ONE Bland episode. The episode releases to
    /// normal pricing (stall perturbation rungs stay active) once it has
    /// made `bland_after` non-degenerate pivots (the degenerate pocket
    /// cleared), this budget is exhausted, or the randomized-deviation rung
    /// fires. Re-entry requires the stall counter to reach `bland_after` ×
    /// 2^episodes (capped ×256) — a growing backoff, so a genuine cycle
    /// keeps getting Bland coverage while a pocket Bland cannot break gives
    /// fast pricing room instead of ever-longer episodes. 0 restores the
    /// legacy TERMINAL behaviour (sticky Bland, never released) as the A/B
    /// lever. Release policy mirrors HiGHS's documented perturbation/
    /// bound-switch release (maxPivot iterations or a non-degenerate
    /// streak, both public).
    pub bland_episode: usize,
    pub max_iterations: usize,
    pub pfi: PfiConfig,
    /// Run the presolve reduction stage before simplex.
    pub presolve: bool,
    /// Fixpoint cap for the presolve loop.
    pub presolve_passes: usize,
    /// Implied-equality conversion + far-side forcing (HiGHS port):
    /// measured with DTE v2 firing on converted rows and still harmful on
    /// Netlib — scrs8 phase 2 312→1827 pivots, e226 +51% iters, geo-mean
    /// 10.3× vs 9.6× off (measured on our benchmark set). The
    /// converted rows are not DTE feeders on this corpus; they are dense
    /// forcing rows our SH/Harris phase 2 grinds on. Default OFF;
    /// `RUSOLVE_PS_IMPEQ=1` remains the A/B lever.
    pub presolve_implied_eq: bool,
    /// Duplicate-col merge: ported + feasibility-verified but net-negative
    /// on Netlib (15,560→16,050 iters, tuff +21%) ⇒ default OFF; `RUSOLVE_PS_PARCOL=1` lever.
    pub presolve_parcol: bool,
    /// Aggregator (implied-free column substitution, HiGHS port): net win
    /// on the 47 MIPLIB relaxations (-6.1% dual iters) and Netlib (-16.9%)
    /// ⇒ default ON;
    /// `RUSOLVE_PRESOLVE_AGG=0` is the CLI lever.
    pub presolve_agg: bool,
    /// Parallel-column dominance fixing (HiGHS :7399–7633 LP arms):
    /// one of a proportional pair is fixed at a bound when the other can
    /// absorb the compensation (opposing bound infinite or row-implied).
    /// Default ON; `RUSOLVE_PS_NO_PARCOLDOM=1` is the A/B lever.
    pub presolve_parcol_dom: bool,
    /// Singleton-cost VUB implied-equality chain (HiGHS
    /// `isImpliedEquation` :263–281 via a column-singleton dual bound):
    /// a cost-carrying column singleton whose improving bound is implied
    /// pins its row's dual away from zero ⇒ row tight at every optimum ⇒
    /// converted to an equation and substituted by DTE. Default ON;
    /// `RUSOLVE_PS_NO_VUBEQ=1` is the A/B lever.
    pub presolve_vub_eq: bool,
    /// Gating experiment: run the (expensive) double-word refinement on
    /// FTRAN b̄ solves too. BTRAN π is ALWAYS refined — pricing is the
    /// correctness surface that made scsd1/scsd6/scrs8 work. If the
    /// bench holds 26/28 with this false, per-pivot cost drops ~40%.
    pub refine_ftran: bool,
    /// Arm the phase-1 exact-tie stability preference
    /// unconditionally (the burst trigger — ≥ 3 sub-1e-7 pivots in the last
    /// 10 — is the default). Set only by the dirty-phase-1-exit re-solve
    /// rung in `solve_dispatch`, whose measured effect is e226 prim+pre
    /// Optimal; never a user lever.
    pub p1_tie_pref_always: bool,
    /// Set on the re-solve rung's config so the rung's own dirty
    /// exits return Limit instead of recursing (the `cleanup_done`
    /// pattern, config-level).
    pub p1_rung_used: bool,
    /// Steepest-edge pricing (rate = d²/‖B⁻¹a‖² with maintained weights;
    /// Goldfarb & Reid 1977, "A practicable steepest-edge simplex algorithm";
    /// Forrest & Goldfarb 1992, "Steepest-edge simplex algorithms for linear
    /// programming").
    pub sh_pricing: bool,
    /// Power-of-2 equilibration scaling (crate::scaling) wrapped around
    /// the presolve+simplex pipeline in solve(). Netlib A/B:
    /// big wins on the degenerate grinders (tuff 2205→518, scrs8 −219,
    /// sctap1 −145, standmps −141, finnis −132) but e226 577→3840 (phase-1
    /// stall interaction; its range improvement is 1837×, so no gate
    /// threshold separates it from the winners) and capri +155 / stair +134.
    /// Not default because the acceptance bar is "must not regress iters
    /// anywhere": DEFAULT OFF behind the `RUSOLVE_SCALE=1` CLI lever until
    /// the scaled-phase-1 stall on e226 is understood.
    pub scaling: bool,
    /// Phase-2 sub-tolerance ENTERING tie perturbation — the phase-2
    /// counterpart of the phase-1 tie perturbation (Beale 1968; Chvátal 1983
    /// §7.3), same safety
    /// invariant: the acceptance test stays the UNPERTURBED `|d| > ptol`
    /// improver set, eps < ptol/2 only reorders the entering ranking.
    /// Default ON (fires only on 4×bland_after Bland stalls, so non-stalling
    /// runs stay bit-identical); `RUSOLVE_NO_P2PERT=1` is the CLI lever back.
    pub p2_tie_perturb: bool,
    /// Phase-1 FROZEN-VERTEX escape. The phase-1 tie perturbation builds its
    /// eps vector ONCE from a fixed hash; on qap10/brandy the resulting entering
    /// ranking + Bland leaving tie-break is a stable fp attractor — w_true
    /// bit-frozen for 98k pivots at ONE vertex.
    /// When w (basic-artificial sum) and the basic-artificial count are
    /// BIT-identical across a 4·bland_after stall window — through a drift
    /// correction (refactorize + `refresh_x_b`), so drift cannot masquerade
    /// as a freeze — the escape engages: (a) eps regenerates from a fresh
    /// deterministic seed (HiGHS releases/re-applies cost perturbation
    /// rather than trusting one fixed vector — HEkkPrimal, Iballe et al.
    /// 2023 public docs; Beale 1968), measured INSUFFICIENT alone (999
    /// regenerations left qap10's vertex closed: every improver pivots at
    /// θ=0), and (b) randomized LEAVING tie-break among exact min-ratio
    /// ties (reservoir draws on the run's deterministic RNG): θ is
    /// unchanged, so feasibility and the improver set are untouched — the
    /// false-Infeasible class that killed entering-randomization in phase 1
    /// (tuff) cannot occur — while each θ=0 pivot moves to a DIFFERENT
    /// adjacent vertex. Engaged until w moves again, then releases.
    /// Non-stalling and w-moving runs never fire either half ⇒ bit-identical
    /// to the one-shot perturbation. `RUSOLVE_NO_P1REPERT=1` is the A/B lever back.
    pub p1_re_perturb: bool,
    /// Run the dual simplex engine (`crate::dual`) instead of the primal
    /// two-phase loop. Library default OFF; the `rusolve solve` CLI turns it
    /// ON, `RUSOLVE_PRIMAL=1` is the lever back. The MIP solver keeps its own
    /// `MipConfig::lp` (primal, `--dual-lp` opt-in).
    /// Lever-off is bit-identical to the primal-only path; a dual run that hits
    /// numerical distress falls back to the primal solver for that LP.
    pub dual: bool,
    /// Dual accuracy regime: the
    /// dual loop runs UNREFINED solves, breaks CHUZR / ratio-group ties on a
    /// 1e-9 relative tolerance class by smallest index, and monitors
    /// |α_col − α_row| per pivot (HiGHS `HEkkDual::updateVerify`), abandoning
    /// the pivot, reinverting and re-arming refinement on trouble.
    /// Component mask `RUSOLVE_DUAL_ACC_PARTS` (default 25). The engine
    /// runs the regime by DEFAULT (once the dano3_3 ladder blow-up that
    /// first held it back was fixed); `RUSOLVE_DUAL_ACC=0` turns it off
    /// for LP and MIP alike. This field only forces it on.
    pub dual_acc: bool,
    /// Cold-start routing policy: on a
    /// COLD primal-engine solve (`dual == false`), check the slack-basis
    /// start for primal infeasibility first and route the LP to the dual
    /// engine when the start violates — the artificial-basis phase 1 that
    /// the primal engine would otherwise pay owns 63–100 % of its
    /// iterations on all 19 blowup instances of our 47 relaxed MIPLIB LPs
    /// while the dual engine solves the same LPs at 0.62–1.71× HiGHS
    /// iterations (measured on our benchmark set). This is HiGHS's
    /// own rule, applied at its unscaled re-solve (HApp.h:424–439, MIT:
    /// "There are unscaled primal infeasibilities, so use dual simplex");
    /// HiGHS picks the dual on 46/46 simplex-solved instances of the same
    /// set. A feasible start keeps the primal path (phase 1 is then
    /// near-free); a routed dual run that hits numerical distress falls
    /// back to the cold primal solve, and the dirty-exit rung stays
    /// primal (`p1_rung_used` skips routing).
    ///
    /// DEFAULT OFF — a default flip needs the full MIP benchmark: the
    /// MIP solver's cold LPs run through this same dispatch, and a dual
    /// default there measured as a loss (tree perturbation swamps the
    /// per-LP win). `RUSOLVE_P1_ROUTE_DUAL=1` is the CLI lever.
    pub p1_route_dual: bool,
    /// Direct CSC-only `Standard` solve path for [`SparseLp`]
    /// input. `RUSOLVE_SPARSE=1` is the
    /// CLI lever. Only `solve_sparse` consumes it: with the lever on AND
    /// presolve/dual off (both consume the [`Lp`] frame) the solve path
    /// builds `Standard` straight from the caller's CSC — no `Lp` hop at
    /// all. Lever off = the (bit-identical, no-fp-arithmetic) `to_lp()` hop
    /// and delegate to `solve`; since `Lp` itself is CSC, the two
    /// positions differ only in that hop.
    pub sparse: bool,
    /// (MIP lever `--warm-unfit-last`, default off): a warm dual entry
    /// whose only dual infeasibilities sit on boxed nonbasic columns flips
    /// them to their other bound instead of rejecting the snapshot, as
    /// HiGHS `HEkkDual::correctDualInfeasibilities` does
    /// (simplex/HEkkDual.cpp:2388–2450, MIT, ERGO-Code/HiGHS). Infeasible
    /// duals on a column with an infinite side still reject (HiGHS shifts
    /// costs there; we keep the honest cold fallback).
    pub warm_flip_boxed: bool,
}

impl Default for SimplexConfig {
    fn default() -> Self {
        Self {
            pivot_tol: 1e-9,
            pricing_tol: 1e-7,
            bland_after: 25,
            // One bland_after-burst episode per engagement. Netlib
            // sweep: beats sticky-Bland on every stalling instance (scsd6
            // p2 277→214, tuff 278→198, scsd1 166→131, standata/standmps
            // −27/−5), bit-identical elsewhere; the re-entry threshold
            // backoff still lets the stall counter reach the 4× rung
            // (engage at ≥100 fires the tie-perturbation rung immediately).
            bland_episode: 25,
            max_iterations: 100_000,
            pfi: PfiConfig::default(),
            presolve: true,
            presolve_passes: 50,
            presolve_implied_eq: false,
            presolve_parcol: true,
            presolve_agg: true,
            presolve_parcol_dom: true,
            presolve_vub_eq: true,
            refine_ftran: true,
            p1_tie_pref_always: false,
            p1_rung_used: false,
            warm_flip_boxed: false,
            // SH default ON: Netlib A/B 25.8x
            // wall (tie vs 25.9x off) but -37% total iters (60,480 vs
            // 95,526), wins on every >0.05s instance except three <=1.10x,
            // scsd6 recovered (0.73x), 27/28, every check passing with SH
            // on; the scsd6-SH 2-cycle that held SH off was cured by the
            // artificial deletion rules. RUSOLVE_NO_SH=1 = lever back.
            sh_pricing: true,
            scaling: false,
            // Phase-2 entering tie perturbation: default ON — it is
            // stall-triggered only (4×bland_after Bland stalls), so non-
            // stalling runs never build eps ⇒ bit-identical without it.
            // RUSOLVE_NO_P2PERT=1 = lever back (CLI).
            p2_tie_perturb: true,
            // Phase-1 re-perturbation: default ON — fires only when a fired
            // phase-1 eps is w-frozen across a 4·bland_after stall window, so
            // every other run stays bit-identical (see field docs).
            p1_re_perturb: true,
            // Dual engine: library default OFF; RUSOLVE_DUAL=1 = lever on (CLI).
            dual: false,
            // Forces the dual accuracy regime on; default OFF (see field docs).
            dual_acc: false,
            // Cold-start routing: default OFF (see field docs).
            // RUSOLVE_P1_ROUTE_DUAL=1 = lever on.
            p1_route_dual: false,
            // Sparse solve path: default OFF — the delegated `Lp` path
            // is the full-featured one (presolve, dual, warm starts) and the
            // two positions are bit-identical anyway; RUSOLVE_SPARSE=1 =
            // lever on (CLI).
            sparse: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SimplexSolution {
    pub status: SimplexStatus,
    /// Caller-space objective value (original sense, original variable space).
    pub objective: f64,
    pub x: Vec<f64>,
    pub phase1_iterations: usize,
    pub phase2_iterations: usize,
    pub bland_triggered: bool,
    /// True when an IterationLimit verdict came from the phase-1
    /// exit gate (pricing said "no improver" but basic values are out of
    /// range — the basis is numerically unusable). `solve_dispatch` runs one
    /// bounded re-solve rung on this marker (stability-armed tie preference);
    /// every other IterationLimit keeps it false.
    pub dirty_p1_exit: bool,
    /// Row duals π (length `num_row`, caller row order/senses), populated
    /// only when [`Self::duals_valid`]. Sign convention and worked example:
    /// [`Self::reduced_costs`] (the pair `(π, d)` is defined by the single
    /// identity `d_j = c_j − πᵀa_j` in the caller's own space).
    pub row_duals: Vec<f64>,
    /// Reduced costs d (length `num_col`, ORIGINAL/pre-presolve column order),
    /// populated only when [`Self::duals_valid`].
    ///
    /// Sign convention: π and d are
    /// reported in the CALLER's space and sense — original rows, original
    /// columns, objective `c` exactly as given in [`Lp::objective`] (NO
    /// minimisation normalisation) — and satisfy for every column j
    ///
    /// ```text
    /// d_j = c_j − πᵀa_j
    /// ```
    ///
    /// At an Optimal solution the sign pattern is
    /// * minimize (`maximize == false`): `d_j ≥ 0` for j at its lower bound,
    ///   `d_j ≤ 0` at its upper bound, `d_j = 0` strictly between bounds;
    ///   `π_i ≤ 0` on 'L' rows, `π_i ≥ 0` on 'G' rows, free on 'E' rows;
    /// * maximize: everything flips — `d_j ≤ 0` at lower bound, `d_j ≥ 0` at
    ///   upper bound; `π_i ≥ 0` on 'L' rows, `π_i ≤ 0` on 'G' rows.
    ///
    /// Worked sign example (maximization, one 'L' row):
    /// `max 3x + 2y  s.t.  x + y ≤ 4,  x, y ≥ 0` ⇒ `x* = 4`, obj 12,
    /// `π = +3` (≥ 0 on 'L' ✓), `d = (3 − 3, 2 − 3) = (0, −1)`: x is basic
    /// with `d = 0`; y sits at its lower bound with `d = −1 ≤ 0` — raising y
    /// from 0 costs one unit of objective per unit, which is exactly why
    /// maximizing pricing looks for `d > 0` at the lower bound. πᵢ doubles
    /// as the sensitivity `∂obj*/∂rhsᵢ` while the optimal basis is unchanged
    /// (here +1 on rhs buys +3 of objective).
    ///
    /// This is HiGHS's own reported-solution convention (`row_dual` /
    /// `col_dual` for the caller's sense, no flip needed) — pinned by
    /// `row_duals_match_highspy`.
    pub reduced_costs: Vec<f64>,
    /// True iff [`Self::row_duals`] and [`Self::reduced_costs`] describe an
    /// Optimal solution of the CALLER's LP (lengths `num_row`/`num_col`).
    ///
    /// LOUD SCOPE NOTE: exact duals are produced ONLY on the presolve-OFF
    /// path — which is what MIP node LPs use — by BOTH engines: the primal
    /// two-phase exit and the dual engine (`SimplexConfig::dual`, and every
    /// `crate::dual::reoptimize` Optimal). On the presolve-ON path dual
    /// postsolve is not implemented, so even at Optimal these fields stay
    /// empty and this flag stays false: duals of a REDUCED LP do not
    /// post-solve back to the caller's LP, and reporting them would be
    /// reporting wrong numbers. `duals_valid == false` means "no
    /// information", never "zero".
    pub duals_valid: bool,
    /// `Unbounded` certificate (presolve-OFF path only): a direction d
    /// (length `num_col`, caller column space) such that from the last
    /// finite iterate x̂ of the pivot loop (not reported), `x̂ + θ·d` stays
    /// feasible for every θ ≥ 0 — equivalently `A·d` is 0 on 'E' rows, ≤ 0
    /// on 'L' rows and ≥ 0 on 'G' rows; `d_j > 0` forces `upper[j] = +inf`
    /// and `d_j < 0` forces `lower[j] = −inf`; and the objective runs away
    /// along it: `cᵀd < 0` for minimize, `cᵀd > 0` for maximize. Assembled
    /// from the FTRAN'd entering column at the unbounded exit, i.e. the
    /// direction that already survived the pivot loop's fresh-factor ray
    /// verification. `None` on every other status/path (honest, not zero).
    pub primal_ray: Option<Vec<f64>>,
    /// `Infeasible` certificate (presolve-OFF path only): Farkas vector
    /// π_f (length `num_row`, caller row order/senses) with
    ///
    /// ```text
    /// π_fᵀ rhs  >  sup { π_fᵀ(A x + s) : lower ≤ x ≤ upper, s in row-slack ranges }
    /// ```
    ///
    /// where the row-slack ranges are `sᵢ ≥ 0` ('L'), `sᵢ ≤ 0` ('G'),
    /// `sᵢ = 0` ('E'). The sup is finite only when `π_f ≤ 0` on 'L' rows,
    /// `π_f ≥ 0` on 'G' rows and, per column j, `(π_fᵀa_j) ≤ 0` whenever
    /// `upper[j] = +inf`, `≥ 0` whenever `lower[j] = −inf`, `= 0` when free;
    /// then no feasible point can satisfy its row system, because any
    /// feasible point would need `π_fᵀrhs` to equal a value strictly below
    /// `π_fᵀrhs`. Taken from the phase-1 exit duals (w > 0), mapped through
    /// the row negations, and returned ONLY after this identity verifies
    /// numerically against the ORIGINAL rows (debug-asserted; release
    /// double-check — a failed verification reports `None`, never a wrong
    /// ray). `None` on the presolve-ON path and on phase-1-ray exits.
    pub dual_ray: Option<Vec<f64>>,
    /// Optimal basis snapshot for warm starts, populated only at Optimal on
    /// the presolve-OFF path — both the primal two-phase engine and the
    /// dual engine export it. LOUD SCOPE NOTE (same pattern as
    /// [`Self::duals_valid`]): on the
    /// presolve-ON path the terminal basis belongs to the REDUCED LP and
    /// basis postsolve is not implemented, so this stays `None` there: a
    /// reduced-LP basis is not a basis of the caller's LP, and reporting it
    /// would be reporting wrong numbers. `None` means "no information",
    /// never "any basis will do". The exact caller-visible indexing
    /// convention is documented on [`Basis`]; feed the snapshot back via
    /// [`SimplexSolver::solve_with`] or [`crate::dual::reoptimize`].
    pub basis: Option<Basis>,
}

/// Status of one variable in a [`Basis`] snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BasisStatus {
    /// The variable is basic in the optimal basis (its value may sit
    /// degenerate exactly on a bound; it is still basic).
    Basic,
    /// Nonbasic at the variable's lower bound (row slack: the row is tight).
    AtLower,
    /// Nonbasic at the variable's upper bound.
    AtUpper,
    /// Nonbasic with `lower == upper` (fixed variable). Import requires
    /// `lower[j] == upper[j]`.
    Fixed,
}

/// Caller-visible basis snapshot. One status per ORIGINAL column and per
/// row slack/surplus variable, in an indexing the MIP solver can round-trip
/// across node bound changes.
///
/// Indexing convention (EXACT — this is the contract):
///
/// * `col_status[j]` — status of ORIGINAL column `j`, i.e. the caller's own
///   [`Lp`] column order (NOT any internal standard-form order), interpreted
///   in the CALLER's bound frame (`Lp::lower[j]`, `Lp::upper[j]`):
///   * `Basic` — variable basic in the optimal basis.
///   * `AtLower` — nonbasic at `lower[j]`; only exported when `lower[j]` is
///     finite. A variable with `lower = -inf`, finite `upper` is internally
///     flip-shifted; its only nonbasic position is `x = upper`, so such
///     variables export `AtUpper`, never `AtLower` (import rejects `AtLower`
///     on them).
///   * `AtUpper` — nonbasic at `upper[j]`; only exported when `upper[j]` is
///     finite. Import rejects `AtUpper` on variables with `upper = +inf`.
///   * `Fixed` — nonbasic with `lower[j] == upper[j]` (export never emits
///     `Fixed` for a basic variable: `Basic` wins, round-trip fidelity first;
///     import requires `lower[j] == upper[j]`).
///   * Fully free variables (`lower = -inf, upper = +inf`, internally
///     ±-split into two columns): export `Basic` iff either half is basic,
///     else `AtLower` (both halves nonbasic, `x = 0`). Import puts the x⁺
///     half basic for `Basic` — if the original basis had the x⁻ half basic,
///     the refactorization / feasibility / dual-feasibility checks reject the
///     warm start and the solver honestly falls back to a cold solve (rare;
///     MIP node LPs with bounded variables never hit this).
/// * `row_status[i]` — status of row `i`'s slack/surplus variable, defined in
///   the ORIGINAL row frame, independent of any internal row negation:
///   `'L'` rows `s_i = rhs_i − aᵢᵀx ≥ 0`, `'G'` rows `s_i = aᵢᵀx − rhs_i ≥ 0`
///   (`'G'` is surplus). Only `Basic` (row loose) and `AtLower` (slack at 0 ⇔
///   row tight) are ever exported; `AtUpper`/`Fixed` are rejected on import.
///   `'E'` rows have no slack variable, so their entry describes the row's
///   FIXED ARTIFICIAL instead:
///   `Basic` iff that artificial is basic (a degenerate redundant equality
///   carried by the artificial at 0), `AtLower` otherwise (the equality is
///   carried by structural columns). Import HONORS the entry — `Basic` puts
///   the row's artificial back into the basis, `AtLower` keeps it out. The
///   old placeholder-plus-positional-fill contract exported a basis that
///   could not say WHICH E-row artificials were basic, so the warm import
///   reconstructed a DIFFERENT (often singular or dual-infeasible) basis and
///   every entry was rejected on E-row-heavy LPs; HiGHS likewise keeps a
///   meaningful status for every row.
/// * Counting: the basic set implied by a snapshot is
///   `{Basic columns} ∪ {Basic row slacks} ∪ {E rows whose entry is Basic}`,
///   which for a basis exported from this LP equals `num_row` minus the
///   number of basic NON-E-row artificials the snapshot cannot see (primal
///   phase-2 artificials pinned at 0 on degenerate rows of `'L'`/`'G'` rows
///   — the dual standard form has no such columns, so a DUAL export is
///   exactly described and its re-import is lossless). The primal warm entry
///   accepts the shortfall and backfills artificial columns of rows that can
///   still take one (its primal-feasibility check rejects a wrong
///   reconstruction via cold fallback); the dual warm entry fills it from its
///   E-row fixed-artificial pool (its only pool) and cold-falls-back when
///   that cannot complete a square basis.
///
/// A `Basis` exported from an LP may only be fed to LPs of the SAME shape
/// (`num_col`, `num_row`, `senses` identical); bound and objective changes are
/// the intended use, row/column insertions are not covered, and any
/// shape mismatch falls back to a cold solve (honest, never a wrong answer).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Basis {
    /// Status per ORIGINAL column (length `num_col`, caller column order).
    pub col_status: Vec<BasisStatus>,
    /// Status per row: the slack/surplus variable of a non-'E' row, the fixed
    /// artificial of an 'E' row (length `num_row`, caller row order —
    /// see the [`Basis`] docs for the exact contract).
    pub row_status: Vec<BasisStatus>,
}

impl Default for SimplexSolution {
    /// Field-completion default for constructors that carry NO dual/ray/
    /// basis information (iteration-limit snapshots, the dual engine, the
    /// presolve-ON postsolve path, presolve terminal verdicts). The
    /// `status`/`objective` values are placeholders — every real site
    /// overrides them explicitly.
    fn default() -> Self {
        Self {
            status: SimplexStatus::IterationLimit,
            objective: f64::NAN,
            x: Vec::new(),
            phase1_iterations: 0,
            phase2_iterations: 0,
            bland_triggered: false,
            dirty_p1_exit: false,
            row_duals: Vec::new(),
            reduced_costs: Vec::new(),
            duals_valid: false,
            primal_ray: None,
            dual_ray: None,
            basis: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Structural,
    Slack,
    Artificial,
}

/// Internal standard form after bound shift + slacks + artificials.
/// All rows are equalities with rhs ≥ 0; all variables lb = 0.
/// Free variables: `lb = -inf, ub` finite ⇒ flip-shift x = ub − x̃
/// (x̃ ≥ 0); fully free ⇒ split x = x⁺ − x⁻ (both ≥ 0). `orig`/`cols_of`
/// map original variables to their (1 or 2) structural columns; `shift` /
/// `flip` describe the reconstruction x_j = shift + x̃ | shift − x̃ | x⁺ − x⁻.
///
/// Columns are stored ONLY in
/// the CSC view — the dense `cols: Vec<Vec<f64>>` array (≈230 MB per MIP
/// node LP at CFL scale) is gone. Every dense-column consumer goes through
/// [`Standard::gather_col`] (or the CSC triple directly); pricing already
/// used the CSC view bit-exactly, so the column index space — and every
/// iteration count — is unchanged.
///
/// `pub` ONLY so [`bench_hooks`] can name it for the materialization bench;
/// every field is private and the type is opaque to outside callers.
pub struct Standard {
    m: usize,
    /// Column-sparse view of the column set: column `j` owns `csc_val[csc_ptr[j]
    /// .. csc_ptr[j + 1]]`, paired row indices in `csc_row`. Entries are stored in
    /// increasing row order — the order the dense column walk visits them —
    /// so a sparse dot product sums the same products in the same sequence as
    /// the dense one. Zeros contribute `x * 0.0 == 0.0`, and `acc + 0.0` leaves
    /// `acc` bit-unchanged (for every `acc` but `-0.0`, which lands on `+0.0`),
    /// so skipping them is bit-exact, not an approximation.
    csc_ptr: Vec<usize>,
    csc_row: Vec<usize>,
    csc_val: Vec<f64>,
    b: Vec<f64>,
    ub: Vec<f64>,
    kind: Vec<Kind>,
    num_structural: usize,
    /// Caller-space objective coefficient per structural column.
    cost: Vec<f64>,
    /// Original variable → its structural columns (1, or 2 for split frees).
    cols_of: Vec<Vec<usize>>,
    /// Per-original bound shift (finite lb, or ub for flip-shifted).
    shift: Vec<f64>,
    /// Per-original: reconstruct as shift − x̃ (flip-shifted MI bounds).
    flip: Vec<bool>,
}

impl Standard {
    /// Total column count (structural + slack + artificial).
    #[inline]
    fn ncol(&self) -> usize {
        self.csc_ptr.len() - 1
    }

    /// Sparse `(row, value)` list of column `v`, ascending by row with exact
    /// zeros dropped (the shape `ProductFormFactorization::build_sparse`
    /// wants; bit-identical values to `gather_col(v, 1.0, ..)`).
    fn sparse_col(&self, v: usize) -> Vec<(usize, f64)> {
        let mut c: Vec<(usize, f64)> = (self.csc_ptr[v]..self.csc_ptr[v + 1])
            .filter(|&k| self.csc_val[k] != 0.0)
            .map(|k| (self.csc_row[k], self.csc_val[k]))
            .collect();
        c.sort_unstable_by_key(|&(r, _)| r);
        c
    }

    /// The only dense-column access API: zero-fill `buf`
    /// (length m) and scatter column `v` scaled by `sgn` (±1.0, exact).
    /// Bit-identical to the old `buf.copy_from_slice(&st.cols[v])` /
    /// `.map(|a| -a)`: `fill(0.0)` + scatter produces the same dense vector
    /// the `Vec<f64>` clone held (structurally, not approximately),
    /// and `sgn * v == ∓v` is exact and sign-symmetric (the
    /// `update_eta_matches_prerecomputed` precedent). The O(m) fill is a
    /// memset ≤1 % of the refined FTRAN it feeds; a touched-list is a
    /// `TODO(perf)` gated on a bench number.
    fn gather_col(&self, v: usize, sgn: f64, buf: &mut [f64]) {
        buf.fill(0.0);
        for k in self.csc_ptr[v]..self.csc_ptr[v + 1] {
            buf[self.csc_row[k]] = sgn * self.csc_val[k];
        }
    }
}

/// Bound used for the per-variable shift: finite lb, else finite ub (flip),
/// else 0 (fully free ⇒ split, no shift).
fn shift_of(lp: &Lp, j: usize) -> f64 {
    if lp.lower[j].is_finite() {
        lp.lower[j]
    } else if lp.upper[j].is_finite() {
        lp.upper[j]
    } else {
        0.0
    }
}

fn to_standard(lp: &Lp) -> Result<Standard, SimplexError> {
    let (n, m) = (lp.num_col, lp.num_row);
    if lp.a.ncol() != n
        || lp.a.nrow() != m
        || lp.rhs.len() != m
        || lp.objective.len() != n
        || lp.lower.len() != n
        || lp.upper.len() != n
        || lp.senses.len() != m
    {
        return Err(SimplexError::Unsupported(
            "inconsistent LP dimensions".into(),
        ));
    }
    for (j, &lb) in lp.lower.iter().enumerate() {
        if lb.is_nan() || lp.upper[j].is_nan() {
            return Err(SimplexError::Unsupported(format!(
                "NaN bound on variable {j}"
            )));
        }
        if lb == f64::NEG_INFINITY && lp.upper[j] == f64::INFINITY {
            // fully free: handled by the ±split below (no shift)
        } else if lb == f64::NEG_INFINITY && lp.upper[j] == f64::NEG_INFINITY {
            return Err(SimplexError::Unsupported(format!(
                "empty bounds on variable {j}: ub = -inf"
            )));
        }
        if lp.upper[j] < lb {
            return Err(SimplexError::Unsupported(format!(
                "empty bounds on variable {j}: ub < lb"
            )));
        }
    }
    for (i, &s) in lp.senses.iter().enumerate() {
        if !matches!(s, 'L' | 'G' | 'E') {
            return Err(SimplexError::Unsupported(format!("row {i} sense '{s}'")));
        }
    }

    // The row-form `Lp` carries its matrix in CSC, so the standard-form
    // construction is the CSC-only one (`to_standard_sparse`). The hop
    // through `SparseLp::from_lp` expands each single-sense row to exactly
    // itself (one finite side), reproducing senses/rhs bit-for-bit; the
    // validation above ran first, so `from_lp` cannot panic. Stored zeros
    // dropped at load were ±0.0 no-ops in the dense-era b-shift
    // — the Netlib A/B spine pins the whole path bit-identical.
    to_standard_sparse(&SparseLp::from_lp(lp))
}

/// CSC-only `Standard` construction from [`SparseLp`]:
/// the same column set, the same
/// expansion order (equality rows once; ranged rows hi-L first then lo-G —
/// I2, via the core's authoritative `row_map`), the same shift/negate/flip/
/// split signs, and the same zero-drop (I3) as [`to_standard`]. The csc
/// arrays match the dense construction's ascending scan entry-for-entry
/// (gated by `to_standard_sparse_matches_dense_bit_identical`), and no dense
/// column is ever materialized: the b-shift scatters over the nnz of the
/// shifted columns (O(nnz) instead of O(n·m)), slacks and artificials land
/// as single-entry appended columns in the same order.
fn to_standard_sparse(lp: &SparseLp) -> Result<Standard, SimplexError> {
    let n = lp.a.ncol();
    if lp.col_lower.len() != n || lp.col_upper.len() != n || lp.objective.len() != n {
        return Err(SimplexError::Unsupported(
            "inconsistent LP dimensions".into(),
        ));
    }
    // `SparseLp::new` already rejected NaN/inverted bounds and stored zeros
    // (I3); the one remaining `to_standard` empty-bound shape is ub = -inf.
    for (j, &lb) in lp.col_lower.iter().enumerate() {
        if lb == f64::NEG_INFINITY && lp.col_upper[j] == f64::NEG_INFINITY {
            return Err(SimplexError::Unsupported(format!(
                "empty bounds on variable {j}: ub = -inf"
            )));
        }
    }
    // Dense-row expansion (I2): `row_map` is authoritative, so senses/rhs are
    // derived to land on exactly its slot order (ascending, hi-L before lo-G).
    let map = lp.row_map();
    let mut sense: Vec<char> = Vec::new();
    let mut b: Vec<f64> = Vec::new();
    for (r, slot) in map.iter().enumerate() {
        let (lo, hi) = (lp.row_lower[r], lp.row_upper[r]);
        if slot[0] != NO_ROW {
            debug_assert_eq!(slot[0], sense.len(), "row_map not ascending");
            if lo == hi {
                sense.push('E');
                b.push(lo);
            } else {
                sense.push('L');
                b.push(hi);
            }
        }
        if slot[1] != NO_ROW {
            debug_assert_eq!(slot[1], sense.len(), "row_map not ascending");
            sense.push('G');
            b.push(lo);
        }
    }
    let m = sense.len();
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());

    // Shift each variable to a nonnegative x̂ (the `shift_of` choice), the
    // b-shift scattered over each shifted column's nnz. Per accumulator b[t]
    // the operation sequence is identical to the dense loop (ascending j,
    // one product v·s per nonzero; skipped structural zeros contribute
    // ±0.0 no-ops).
    let shifts: Vec<f64> = (0..n)
        .map(|j| {
            let (lb, ub) = (lp.col_lower[j], lp.col_upper[j]);
            if lb.is_finite() {
                lb
            } else if ub.is_finite() {
                ub
            } else {
                0.0
            }
        })
        .collect();
    for (j, &s) in shifts.iter().enumerate() {
        if s != 0.0 {
            for k in ptr[j]..ptr[j + 1] {
                let (r, v) = (ridx[k] as usize, rval[k]);
                for &t in map[r].iter() {
                    if t != NO_ROW {
                        b[t] -= v * s;
                    }
                }
            }
        }
    }
    // Row normalization: rhs ≥ 0 by negating the row (L↔G, E stays).
    let negated: Vec<bool> = b.iter().map(|&bi| bi < 0.0).collect();
    for i in 0..m {
        if negated[i] {
            b[i] = -b[i];
            sense[i] = match sense[i] {
                'L' => 'G',
                'G' => 'L',
                s => s,
            };
        }
    }

    // Column set built straight into CSC, ascending rows by construction:
    // input columns are ascending (I3 layout contract) and caller row r maps
    // to dense rows map[r] = [hi, lo] with hi < lo, ordered across r.
    let mut csc_ptr: Vec<usize> = Vec::with_capacity(n + 2 * m + 1);
    let mut csc_row: Vec<usize> = Vec::new();
    let mut csc_val: Vec<f64> = Vec::new();
    let mut ub: Vec<f64> = Vec::new();
    let mut kind: Vec<Kind> = Vec::new();
    let mut cost: Vec<f64> = Vec::new();
    let mut cols_of: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut flip = vec![false; n];
    csc_ptr.push(0);
    macro_rules! push_col {
        ($entries:expr, $ubv:expr, $kindv:expr, $costv:expr) => {{
            let entries: &[(usize, f64)] = $entries;
            debug_assert!(entries.windows(2).all(|w| w[0].0 < w[1].0));
            for &(t, x) in entries {
                debug_assert_ne!(x, 0.0, "zeros must not reach the CSC view");
                csc_row.push(t);
                csc_val.push(x);
            }
            csc_ptr.push(csc_row.len());
            ub.push($ubv);
            kind.push($kindv);
            cost.push($costv);
        }};
    }
    let mut base: Vec<(usize, f64)> = Vec::new();
    let mut flipped: Vec<(usize, f64)> = Vec::new();
    for j in 0..n {
        // base column: expanded entries with the row-negation sign folded in
        // (exact sign flip, order-preserving).
        base.clear();
        for k in ptr[j]..ptr[j + 1] {
            let (r, v) = (ridx[k] as usize, rval[k]);
            for &t in map[r].iter() {
                if t != NO_ROW {
                    base.push((t, if negated[t] { -v } else { v }));
                }
            }
        }
        let (lb, ubj) = (lp.col_lower[j], lp.col_upper[j]);
        if lb.is_finite() {
            let c = csc_ptr.len() - 1;
            push_col!(&base, ubj - lb, Kind::Structural, lp.objective[j]);
            cols_of[j].push(c);
        } else if ubj.is_finite() {
            // flip-shift: x = ub − x̂ ⇒ x̂-column is the negated base
            flipped.clear();
            flipped.extend(base.iter().map(|&(t, x)| (t, -x)));
            let c = csc_ptr.len() - 1;
            push_col!(&flipped, f64::INFINITY, Kind::Structural, -lp.objective[j]);
            cols_of[j].push(c);
            flip[j] = true;
        } else {
            // fully free: x = x⁺ − x⁻, both ≥ 0
            let c1 = csc_ptr.len() - 1;
            push_col!(&base, f64::INFINITY, Kind::Structural, lp.objective[j]);
            flipped.clear();
            flipped.extend(base.iter().map(|&(t, x)| (t, -x)));
            let c2 = csc_ptr.len() - 1;
            push_col!(&flipped, f64::INFINITY, Kind::Structural, -lp.objective[j]);
            cols_of[j].extend([c1, c2]);
        }
    }
    let num_structural = csc_ptr.len() - 1;
    for (i, &s) in sense.iter().enumerate() {
        if s != 'E' {
            push_col!(
                &[(i, if s == 'L' { 1.0 } else { -1.0 })],
                f64::INFINITY,
                Kind::Slack,
                0.0
            );
        }
    }
    for i in 0..m {
        push_col!(&[(i, 1.0)], f64::INFINITY, Kind::Artificial, 0.0);
    }

    Ok(Standard {
        m,
        csc_ptr,
        csc_row,
        csc_val,
        b,
        ub,
        kind,
        num_structural,
        cost,
        cols_of,
        shift: shifts,
        flip,
    })
}

/// Not API — bench seam (`benches/numerics.rs`). Exposes the std-form
/// builder and a size probe so the criterion bench can isolate construction
/// cost without widening the solver's public surface. (The dense
/// `to_standard` front door delegates to the CSC build after validating, so
/// there is exactly one construction to bench.)
#[doc(hidden)]
pub mod bench_hooks {
    use super::{Standard, to_standard_sparse};
    use rusolve_core::SparseLp;

    /// Sparse front door: CSC-only `to_standard_sparse`.
    pub fn build_sparse(lp: &SparseLp) -> Standard {
        to_standard_sparse(lp).expect("bench SparseLp must be well-formed")
    }

    /// Work absorbed by the returned Standard (CSC nnz + column count) —
    /// keeps criterion from optimizing the build away.
    pub fn standard_size(st: &Standard) -> usize {
        st.csc_val.len() + st.csc_ptr.len()
    }
}

/// Minimum iterations between maintained-d refreshes triggered by a
/// refactorization (drift control; the d_q consistency check catches the rest).
/// Priceable-and-improving test for a NONBASIC column (bit of the candidate
/// set): same exclusions the former full scan applied, in the same order.
#[inline]
fn is_cand(
    st: &Standard,
    allow_artificial: bool,
    dead_v: bool,
    at_ub_v: bool,
    d: f64,
    ptol: f64,
    v: usize,
) -> bool {
    if st.kind[v] == Kind::Artificial && (!allow_artificial || dead_v) {
        return false;
    }
    if st.ub[v] == 0.0 {
        return false;
    }
    if at_ub_v { d > ptol } else { d < -ptol }
}

const DJ_MIN_AGE: usize = 40;

/// Row-wise copy of the constraint matrix by counting transpose of the CSC;
/// entries within a row come out in ascending column order.
fn build_row_view(st: &Standard) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
    let m = st.m;
    let nnz = st.csc_row.len();
    let mut ptr = vec![0usize; m + 1];
    for &r in &st.csc_row {
        ptr[r + 1] += 1;
    }
    for i in 0..m {
        ptr[i + 1] += ptr[i];
    }
    let mut next = ptr.clone();
    let mut col = vec![0usize; nnz];
    let mut val = vec![0.0f64; nnz];
    for v in 0..st.ncol() {
        for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
            let r = st.csc_row[k];
            col[next[r]] = v;
            val[next[r]] = st.csc_val[k];
            next[r] += 1;
        }
    }
    (ptr, col, val)
}

/// Sparse dot product πᵀa_v over the CSR column v. Bit-identical to the
/// dense zip-sum (see the csc_* doc on `Standard`): same products, same
/// accumulation order, exact zeros skipped.
#[inline]
fn sparse_dot(x: &[f64], st: &Standard, v: usize) -> f64 {
    let (lo, hi) = (st.csc_ptr[v], st.csc_ptr[v + 1]);
    let mut acc = 0.0f64;
    for k in lo..hi {
        acc += x[st.csc_row[k]] * st.csc_val[k];
    }
    acc
}

/// Burn-rollback snapshot taken when the relaxed phase-1 step engages: basis,
/// basic values, at-ub flags, and the probe tuple (w census, n_art_basic)
/// as recorded at engagement (see `p1_snap` in the phase-1 loop).
type P1EngageSnapshot = (Vec<usize>, Vec<f64>, Vec<bool>, (f64, usize));

/// Re-anchor the basic solution: x_b = B⁻¹(b − Σ_{nonbasic at ub} ub·a).
/// The ub-fixed nonbasic contribution is NOT optional — solving B⁻¹b alone
/// silently shifts every basic value (the two unit tests caught exactly
/// that). Clamps fp-noise negatives (−1e-7,0) to 0 like every other refresh
/// site.
fn refresh_x_b(
    pfi: &mut ProductFormFactorization,
    st: &Standard,
    basis: &[usize],
    at_ub: &[bool],
    x_b: &mut [f64],
) {
    let mut rhs: Vec<f64> = st.b.clone();
    let mut in_b = vec![false; st.ncol()];
    for &v in basis.iter() {
        in_b[v] = true;
    }
    for (v, &ib) in in_b.iter().enumerate() {
        if !ib && at_ub[v] && st.ub[v].is_finite() {
            // Scatter over the column's nnz: O(nnz) instead of
            // the old O(m) dense sweep per ub-fixed column, and per
            // accumulator the same −= ub·a products (skipped zeros were
            // ±0.0 no-ops).
            for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                rhs[st.csc_row[k]] -= st.ub[v] * st.csc_val[k];
            }
        }
    }
    pfi.solve_b(&rhs, x_b);
    for xv in x_b.iter_mut() {
        if *xv < 0.0 && *xv > -1e-7 {
            *xv = 0.0;
        }
    }
}

/// Reconstruct caller-space x from basis values (inverse of the bound
/// transforms recorded in `st`: shift, flip-shift, and free-variable split).
/// Check a reported caller-space solution against the caller's OWN rows and
/// bounds (optimal-verdict defense): returns `(worst bound violation, worst
/// row violation)`, both ≥ 0 — the same measure the differential corpus
/// check applies.
/// O(nnz): one bound sweep plus one column scatter per nonzero.
fn primal_violations(lp: &Lp, x: &[f64]) -> (f64, f64) {
    let mut wb = 0.0f64;
    for (j, &v) in x.iter().enumerate() {
        let over = (lp.lower[j] - v).max(v - lp.upper[j]);
        if over > wb {
            wb = over;
        }
    }
    // Row activities via CSC scatter (column-major, ascending rows): per row
    // the accumulation order is ascending j, exactly the dense loop's order —
    // stored zeros and x_j = 0 contribute ±0.0 no-ops.
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
    let mut act = vec![0.0f64; lp.num_row];
    for (j, &xj) in x.iter().enumerate() {
        if xj != 0.0 {
            for k in ptr[j]..ptr[j + 1] {
                let v = rval[k];
                if v != 0.0 {
                    act[ridx[k] as usize] += v * xj;
                }
            }
        }
    }
    let mut wr = 0.0f64;
    for (i, (&sense_i, &rhs_i)) in lp.senses.iter().zip(&lp.rhs).enumerate() {
        let lhs = act[i];
        let v = match sense_i {
            'L' => lhs - rhs_i,
            'G' => rhs_i - lhs,
            _ => (lhs - rhs_i).abs(),
        };
        if v > wr {
            wr = v;
        }
    }
    (wb, wr)
}

/// Feasibility bar for an Optimal verdict: the HiGHS-class
/// `primal_feasibility_tolerance` scaled by solution magnitude — the PASS
/// bar of the differential corpus check (`1e-6 · max(1, |x|∞)`). Clean
/// presolve-off solves land at 1e-9..1e-13, so this never fires on the
/// corpus; the big-M node-LP repros (row viol 42, bound viol 1.02e-4 at |x|∞ = 287)
/// sit far above it.
fn primal_feasible(lp: &Lp, x: &[f64]) -> bool {
    let scale = x.iter().fold(0.0f64, |a, &v| a.max(v.abs())).max(1.0);
    let tol = 1e-6 * scale;
    let (wb, wr) = primal_violations(lp, x);
    wb <= tol && wr <= tol
}

fn extract_x(st: &Standard, basis: &[usize], x_b: &[f64], at_ub: &[bool]) -> Vec<f64> {
    let mut colval = vec![0.0; st.num_structural];
    for (i, &v) in basis.iter().enumerate() {
        if v < st.num_structural {
            colval[v] = x_b[i];
        }
    }
    for c in 0..st.num_structural {
        if at_ub[c] {
            colval[c] = st.ub[c];
        }
    }
    let mut x = vec![0.0; st.cols_of.len()];
    // j indexes x, cols_of and flip/shift jointly.
    #[allow(clippy::needless_range_loop)]
    for j in 0..st.cols_of.len() {
        x[j] = match st.cols_of[j].as_slice() {
            [c] => {
                if st.flip[j] {
                    st.shift[j] - colval[*c]
                } else {
                    colval[*c] + st.shift[j]
                }
            }
            [c1, c2] => colval[*c1] - colval[*c2],
            _ => 0.0,
        };
    }
    x
}

/// Build a caller-visible [`Basis`] snapshot from per-structural-column basic
/// / at-upper flags, per-row-slack basic flags and per-E-row-artificial basic
/// flags (shared by the primal export here and the dual engine's
/// `export_basis_dual`; both standard forms push structural columns in the
/// same original-variable order and slacks in original non-E row order, so
/// the same assembly applies to both).
///
/// `struct_basic`/`struct_at_ub` are indexed by STANDARD column (length
/// `num_structural`), `slack_basic` by non-E row ordinal, and `art_basic` by
/// E-row ordinal (ascending row order — an 'E' row's entry is `Basic`
/// iff its fixed artificial is basic, see [`Basis`] for the contract). The
/// flip/split/fixed status derivation reads only the caller's own bounds —
/// see the [`Basis`] docs for the exact convention.
pub(crate) fn assemble_basis(
    lp: &Lp,
    cols_of: &[Vec<usize>],
    struct_basic: &[bool],
    struct_at_ub: &[bool],
    slack_basic: &[bool],
    art_basic: &[bool],
) -> Basis {
    let mut col_status = vec![BasisStatus::AtLower; lp.num_col];
    for j in 0..lp.num_col {
        let (lb, ub) = (lp.lower[j], lp.upper[j]);
        match cols_of[j].as_slice() {
            [c] => {
                col_status[j] = if struct_basic[*c] {
                    BasisStatus::Basic
                } else if lb == ub {
                    BasisStatus::Fixed
                } else if !lb.is_finite() {
                    // flip-shifted (lb = -inf, finite ub): the only
                    // nonbasic position x̃ = 0 IS caller x = ub.
                    BasisStatus::AtUpper
                } else if struct_at_ub[*c] {
                    BasisStatus::AtUpper
                } else {
                    BasisStatus::AtLower
                };
            }
            [c1, c2] if struct_basic[*c1] || struct_basic[*c2] => {
                // fully-free ±split: Basic iff either half is basic.
                col_status[j] = BasisStatus::Basic;
            }
            _ => {} // every original variable owns 1 or 2 columns; both
                    // halves of a split free var nonbasic ⇒ AtLower (x = 0).
        }
    }
    let mut next = 0usize;
    let mut next_art = 0usize;
    let row_status = lp
        .senses
        .iter()
        .map(|&sn| {
            let basic = if sn == 'E' {
                // The 'E' entry describes the row's fixed artificial.
                let b = art_basic[next_art];
                next_art += 1;
                b
            } else {
                let b = slack_basic[next];
                next += 1;
                b
            };
            if basic {
                BasisStatus::Basic
            } else {
                BasisStatus::AtLower
            }
        })
        .collect();
    Basis {
        col_status,
        row_status,
    }
}

/// Export the primal engine's terminal state as a [`Basis`] snapshot.
/// An 'E' row's artificial is visible through the row entry; a basic
/// NON-E-row artificial therefore still shrinks the implied basic set below
/// `num_row`, which the warm entries handle per the [`Basis`] counting
/// contract.
fn export_basis(st: &Standard, lp: &Lp, basis: &[usize], at_ub: &[bool]) -> Basis {
    let total = st.ncol();
    let mut in_b = vec![false; total];
    for &v in basis {
        in_b[v] = true;
    }
    let n_slacks = total - st.num_structural - st.m;
    let mut slack_basic = vec![false; n_slacks];
    for k in 0..n_slacks {
        slack_basic[k] = in_b[st.num_structural + k];
    }
    // Primal standard form: one artificial per row, row order preserved,
    // so row i's artificial is standard column (total - m) + i.
    let art_start = total - st.m;
    let art_basic: Vec<bool> = (0..lp.num_row)
        .filter(|&i| lp.senses[i] == 'E')
        .map(|i| in_b[art_start + i])
        .collect();
    assemble_basis(
        lp,
        &st.cols_of,
        &in_b[..st.num_structural],
        &at_ub[..st.num_structural],
        &slack_basic,
        &art_basic,
    )
}

/// Warm-start feasibility tolerance: the same 1e-7 budget the phase-1
/// exit test and the `refresh_x_b` round-off clamp use — a basic value
/// outside the box by more than this is a dual-start case, not noise.
const WARM_FEAS_TOL: f64 = 1e-7;

/// Interpret one [`BasisStatus`] for one ORIGINAL column in the primal warm
/// entry. `Some(plan)` = how the standard form must place this variable;
/// `None` = the status contradicts the caller's bounds (contract violation,
/// see [`Basis`]) and the whole warm start falls back to cold.
pub(crate) enum WarmVar {
    /// Standard column basic (for a ±split free var: the x⁺ half enters).
    Basic,
    /// Nonbasic at x̃ = 0.
    NonLower,
    /// Nonbasic at x̃ = ũ (only reachable with finite lb AND finite ub).
    NonUpper,
}

pub(crate) fn warm_plan(status: BasisStatus, lb: f64, ub: f64) -> Option<WarmVar> {
    match status {
        BasisStatus::Basic => Some(WarmVar::Basic),
        BasisStatus::AtLower => {
            // finite lb: nonbasic at lb. flip-shift (lb=-inf): its only
            // nonbasic position IS x = ub, reported AtUpper — AtLower here
            // means the caller edited the status into an impossible slot.
            // fully free: both halves nonbasic at 0 (x = 0), valid.
            if lb.is_finite() || (lb == f64::NEG_INFINITY && ub == f64::INFINITY) {
                Some(WarmVar::NonLower)
            } else {
                None
            }
        }
        BasisStatus::AtUpper => {
            if ub == f64::INFINITY {
                return None; // no upper side to sit at
            }
            // flip-shift: x = ub ⇔ x̃ = 0, i.e. the LOWER side internally.
            Some(if lb.is_finite() {
                WarmVar::NonUpper
            } else {
                WarmVar::NonLower
            })
        }
        BasisStatus::Fixed => {
            if lb == ub {
                Some(WarmVar::NonLower) // value = lb = ub, side irrelevant
            } else {
                None // "Fixed" with a range: which value? never guess.
            }
        }
    }
}

/// Rebuild the primal standard-form basis (`basis` = column per row slot,
/// `at_ub` = nonbasic upper-side flags) from a caller [`Basis`]. `None` ⇒
/// snapshot unusable (shape mismatch, impossible status, over-full basic
/// set): the caller cold-starts, honestly. 'E'-row entries are HONORED:
/// `Basic` re-enters that row's artificial, `AtLower` keeps it out. A basic-
/// set SHORTFALL is then the invisible basic NON-E-row artificials, backfilled
/// from the rows that can still take one (the caller's feasibility check
/// verifies the reconstruction — a backfilled artificial must sit at 0,
/// exactly the condition for validity).
fn import_basis_primal(lp: &Lp, st: &Standard, snap: &Basis) -> Option<(Vec<usize>, Vec<bool>)> {
    if snap.col_status.len() != lp.num_col || snap.row_status.len() != lp.num_row {
        return None; // row/col churn is not this contract
    }
    let m = st.m;
    let total = st.ncol();
    let art_start = total - m;
    let mut at_ub = vec![false; total];
    let mut basis: Vec<usize> = Vec::with_capacity(m);
    for j in 0..lp.num_col {
        let plan = warm_plan(snap.col_status[j], lp.lower[j], lp.upper[j])?;
        match st.cols_of[j].as_slice() {
            [c] => match plan {
                WarmVar::Basic => basis.push(*c),
                WarmVar::NonLower => {}
                WarmVar::NonUpper => at_ub[*c] = true,
            },
            [c1, _c2] => match plan {
                // Known lossiness (Basis docs): the x⁺ half re-enters
                // regardless of which half was originally basic; the
                // feasibility check rejects a wrong reconstruction.
                WarmVar::Basic => basis.push(*c1),
                WarmVar::NonLower => {}
                WarmVar::NonUpper => return None, // halves have no finite upper side
            },
            _ => return None,
        }
    }
    // Row slacks: k-th non-E row owns standard column num_structural + k.
    // 'E' rows carry their artificial's status: Basic re-enters the
    // row's artificial (standard column art_start + i), AtLower keeps it out.
    let mut k = 0usize;
    let mut art_pushed = vec![false; m];
    let mut slack_basic_row = vec![false; m];
    for i in 0..lp.num_row {
        let e_row = lp.senses[i] == 'E';
        match (e_row, snap.row_status[i]) {
            (true, BasisStatus::Basic) => {
                basis.push(art_start + i);
                art_pushed[i] = true;
            }
            (true, BasisStatus::AtLower) => {}
            (true, _) => return None, // artificials have no upper/fixed side
            (false, BasisStatus::Basic) => {
                basis.push(st.num_structural + k);
                slack_basic_row[i] = true;
            }
            (false, BasisStatus::AtLower) => {}
            (false, _) => return None, // slacks have no upper/fixed side
        }
        if !e_row {
            k += 1;
        }
    }
    if basis.len() > m {
        return None; // over-full: not a basis
    }
    // Shortfall = basic NON-E-row artificials the snapshot cannot see.
    // Backfill in row order, skipping rows that cannot take an artificial
    // (its own artificial already honored in, or its slack basic — e_i twice
    // is singular); the caller's primal-feasibility check verifies the
    // reconstruction (a backfilled artificial must sit at 0, exactly the
    // condition for validity).
    for r in 0..m {
        if basis.len() == m {
            break;
        }
        if !art_pushed[r] && !slack_basic_row[r] {
            basis.push(art_start + r);
        }
    }
    if basis.len() != m {
        return None; // pool exhausted: the snapshot does not describe a basis
    }
    Some((basis, at_ub))
}

pub struct SimplexSolver {
    config: SimplexConfig,
    /// Solve-internal wall-clock deadline (`None` = off, the default).
    /// Checked inside the pivot loop; see [`Self::set_deadline`].
    deadline: Option<std::time::Instant>,
    /// True once the Optimal verdict has already been through the primal
    /// cleanup ladder ([`SimplexSolver::primal_violations`] +
    /// `solve_std`/`solve_cold` cleanup arms). The ladder's rung re-solves and
    /// the presolve-off re-solve run with this set, so the audit→re-solve
    /// chain cannot recurse: every rung returns its audit-best point whatever
    /// the verdict, and only the outer entry gets to run the ladder.
    cleanup_done: bool,
}

#[derive(PartialEq, Eq)]
enum PhaseOutcome {
    Optimal,
    Unbounded,
    Limit,
    /// A phase-1 basic value left the credible range entirely
    /// (|x_b| > 1e12; measured e226 detonation reaches 1e19). The basis has
    /// numerically detonated — pivoting on is a refactor-storm crawl
    /// (measured ~65 ms/pivot on e226) that ends at the phase-1 exit gate
    /// on a garbage point anyway. The caller turns this into an honest
    /// Limit marked `dirty_p1_exit`, arming the `solve_dispatch` re-solve
    /// rung. Emitted only with `allow_artificial` (phase 1).
    Detonated,
    /// SH pricing steered into a numerically rank-deficient basis (fake ray:
    /// refactorize Err(Singular) or π/b̄ dual-identity mismatch). Only the
    /// SH path may emit this; the caller restarts phase 2 Dantzig-priced from
    /// the phase-1 exit snapshot — SH can never be MORE wrong than baseline.
    NumericalDistress,
}

/// What a [`SimplexSolver::run_phase`] run leaves behind at its exit basis,
/// for the caller-space dual/ray exports (read-only capture — nothing here is
/// recomputed after the loop; the vectors are the live pricing quantities at
/// the moment the phase terminated).
#[derive(Default)]
struct PhaseExit {
    /// BTRAN duals π (Bᵀπ = c_B in standard space), captured at the Optimal
    /// exit — the same π the final no-improver pricing pass consumed. Left
    /// empty on other outcomes.
    pi: Vec<f64>,
    /// At the Unbounded exit: `(entering column q, direction sign, b̄ =
    /// B⁻¹(sgn·a_q))` — the FTRAN'd column of the very pivot whose ratio
    /// test found no blocker, i.e. the direction that survived the
    /// fresh-factor ray verification.
    ray: Option<(usize, f64, Vec<f64>)>,
}

impl SimplexSolver {
    pub fn new(config: SimplexConfig) -> Self {
        Self {
            config,
            deadline: None,
            cleanup_done: false,
        }
    }

    /// Solve-internal wall-clock deadline for this solver's runs.
    /// The tree-level limit checks (`Search::limit_hit`) fire only between
    /// node LPs, so a solve-internal cycle — mik's restart round 2 spins in
    /// a phase-1 degenerate oscillation (the dano3_3/brandy degeneracy
    /// class) — is bounded only by the per-phase
    /// `max_iterations` cap, and 100k pivots of a 2 700-row LP outlast any
    /// MIP time limit. HiGHS checks the timer inside the iteration loop and
    /// exits with kTimeLimit/kIterationLimit (`HEkk::bailout`,
    /// reference/HiGHS/highs/simplex/HEkk.cpp:3441–3461, MIT; we port the
    /// mechanism, not the code). A deadline hit returns the existing
    /// IterationLimit status — every caller already treats that as "LP not
    /// solved" — so the MIP layer can bail out cleanly. Scope: the cold
    /// two-phase primal engine (`run_phase`), the home of the deadline-bounded spin; warm dual
    /// reoptimization is already bounded per call by its iteration limit.
    pub fn set_deadline(&mut self, deadline: Option<std::time::Instant>) {
        self.deadline = deadline;
    }

    /// Solve: optionally presolve, then simplex, then postsolve.
    ///
    /// Cold start. To warm-start from an exported [`Basis`] snapshot use
    /// [`Self::solve_with`]; the warm dual reopt entry is
    /// [`crate::dual::reoptimize`].
    pub fn solve(&self, lp: &Lp) -> Result<SimplexSolution, SimplexError> {
        self.solve_with(lp, None)
    }

    /// Solve, optionally warm-starting from a [`Basis`] snapshot.
    ///
    /// `Some(basis)` on the presolve-OFF, primal-lever-OFF surface runs the
    /// warm path: refactorize the given basis, verify PRIMAL feasibility of
    /// the induced basic solution (tol `WARM_FEAS_TOL`); on any rejection —
    /// singular refactorization, infeasible basic values, snapshot/shape
    /// contract violation — fall back to the exact cold solve (honest: the
    /// answer is always right, only the iteration count differs; a warm run
    /// reports `phase1_iterations == 0`). A warm basis that is primal
    /// feasible but not optimal is repaired by phase 2 — this is the
    /// primal-side warm start; the MIP node case (dual feasible, primal
    /// infeasible) belongs to [`crate::dual::reoptimize`].
    ///
    /// Scope (same loud contract as [`SimplexSolution::basis`]): with
    /// `config.presolve` on the snapshot is IGNORED — the solver would have
    /// to pivot on the REDUCED LP, whose basis is not the caller's. With
    /// `config.dual` on the snapshot is used as a warm DUAL start via
    /// `reoptimize`.
    pub fn solve_with(
        &self,
        lp: &Lp,
        basis: Option<&Basis>,
    ) -> Result<SimplexSolution, SimplexError> {
        let Some(snap) = basis else {
            return self.solve_cold(lp);
        };
        if self.config.presolve {
            return self.solve_cold(lp); // snapshot belongs to the wrong LP
        }
        if self.config.dual {
            // Lever-on: the snapshot is a warm DUAL start. Pruned is
            // unreachable without a cutoff; the verdicts map 1:1.
            let r = crate::dual::reoptimize(lp, &self.config, snap, None, None)?;
            let status = match r.status {
                crate::dual::ReoptStatus::Optimal => SimplexStatus::Optimal,
                crate::dual::ReoptStatus::Infeasible => SimplexStatus::Infeasible,
                crate::dual::ReoptStatus::Unbounded => SimplexStatus::Unbounded,
                crate::dual::ReoptStatus::Pruned | crate::dual::ReoptStatus::IterationLimit => {
                    SimplexStatus::IterationLimit
                }
            };
            return Ok(SimplexSolution {
                status,
                objective: r.objective,
                x: r.x,
                phase1_iterations: 0,
                phase2_iterations: r.iters,
                bland_triggered: false,
                basis: r.basis,
                ..Default::default()
            });
        }
        match self.solve_warm(lp, snap)? {
            Some(sol) => Ok(sol),
            None => self.solve_dispatch(lp),
        }
    }

    /// The warm primal start. `Ok(None)` ⇒ snapshot rejected before any
    /// pivot: the caller cold-starts.
    fn solve_warm(&self, lp: &Lp, snap: &Basis) -> Result<Option<SimplexSolution>, SimplexError> {
        let mut st = to_standard(lp)?;
        let m = st.m;
        let total = st.ncol();
        let art_start = total - m;
        let Some((mut basis, mut at_ub)) = import_basis_primal(lp, &st, snap) else {
            return Ok(None);
        };
        // Phase-2 semantics for artificials (the boeing2 pin, see the
        // phase-boundary comment in solve_cold): ub = 0, so a backfilled
        // artificial is valid only at value 0 — exactly what the feasibility
        // check below verifies.
        for v in art_start..total {
            st.ub[v] = 0.0;
        }
        let basis_cols: Vec<Vec<(usize, f64)>> = (0..m).map(|i| st.sparse_col(basis[i])).collect();
        let Ok(mut pfi) = ProductFormFactorization::build_sparse(basis_cols, self.config.pfi)
        else {
            return Ok(None); // singular snapshot basis ⇒ cold
        };
        let mut x_b = vec![0.0; m];
        refresh_x_b(&mut pfi, &st, &basis, &at_ub, &mut x_b);
        // Primal feasibility verification: the warm basis must induce a
        // primal feasible basic solution (tol WARM_FEAS_TOL, the same 1e-7
        // budget the phase-1 exit and refresh clamps use). Anything beyond
        // noise is a DUAL-start case (reoptimize), not a primal one.
        for i in 0..m {
            let v = basis[i];
            if x_b[i] < -WARM_FEAS_TOL
                || (st.ub[v].is_finite() && x_b[i] > st.ub[v] + WARM_FEAS_TOL)
            {
                return Ok(None);
            }
            // Snap the accepted noise into the box (refresh_x_b already
            // clamped (−1e-7, 0); this covers the upper side).
            if x_b[i] < 0.0 {
                x_b[i] = 0.0;
            } else if st.ub[v].is_finite() && x_b[i] > st.ub[v] {
                x_b[i] = st.ub[v];
            }
        }

        let sign = if lp.maximize { -1.0 } else { 1.0 };
        let mut cost2 = vec![0.0; total];
        for (d, &cc) in cost2.iter_mut().zip(&st.cost) {
            *d = cc * sign;
        }
        let mut bland = false;
        let mut randomized = false;
        let mut sh_w: Vec<f64> = Vec::new();
        let mut exit2 = PhaseExit::default();
        let (it2, outcome) = self.run_phase(
            &st,
            &mut pfi,
            &mut basis,
            &mut x_b,
            &mut at_ub,
            &cost2,
            false,
            &mut bland,
            &mut randomized,
            &mut sh_w,
            self.config.sh_pricing,
            0,
            &mut exit2,
        )?;
        if outcome == PhaseOutcome::NumericalDistress {
            return Ok(None); // honest: cold restart from scratch
        }
        if outcome == PhaseOutcome::Unbounded {
            // A ray found from a primal feasible start is real (the same
            // fresh-factor verification the cold path relies on).
            let primal_ray = exit2
                .ray
                .as_ref()
                .and_then(|(q, sgn, bbar)| primal_ray_x(lp, &st, &basis, *q, *sgn, bbar));
            return Ok(Some(SimplexSolution {
                status: SimplexStatus::Unbounded,
                objective: f64::NEG_INFINITY,
                x: vec![],
                phase1_iterations: 0,
                phase2_iterations: it2,
                bland_triggered: bland,
                primal_ray,
                ..Default::default()
            }));
        }
        if outcome == PhaseOutcome::Limit || outcome == PhaseOutcome::Detonated {
            // Detonated cannot fire on this warm phase-2 path
            // (allow_artificial=false); treated as Limit defensively.
            return Ok(Some(limit_solution(
                &st,
                &lp.objective,
                &basis,
                &x_b,
                &at_ub,
                0,
                it2,
                bland,
            )));
        }
        // Optimal tail — mirrors solve_cold: final re-anchor on fresh factors,
        // extract, caller-space duals, basis export.
        if pfi.refactorize().is_ok() {
            refresh_x_b(&mut pfi, &st, &basis, &at_ub, &mut x_b);
        }
        let x = extract_x(&st, &basis, &x_b, &at_ub);
        let objective: f64 = lp.objective.iter().zip(&x).map(|(c, xv)| c * xv).sum();
        let (row_duals, reduced_costs, duals_valid) = if exit2.pi.len() == m {
            let (rd, rc) = caller_duals(lp, sign, &exit2.pi);
            (rd, rc, true)
        } else {
            (Vec::new(), Vec::new(), false)
        };
        Ok(Some(SimplexSolution {
            status: SimplexStatus::Optimal,
            objective,
            x,
            phase1_iterations: 0, // warm: no phase 1 was needed
            phase2_iterations: it2,
            bland_triggered: bland,
            dirty_p1_exit: false,
            row_duals,
            reduced_costs,
            duals_valid,
            primal_ray: None,
            dual_ray: None,
            basis: Some(export_basis(&st, lp, &basis, &at_ub)),
        }))
    }

    /// Cold pipeline body (former `solve`): presolve stage + dispatch.
    fn solve_cold(&self, lp: &Lp) -> Result<SimplexSolution, SimplexError> {
        if !self.config.presolve {
            return self.solve_dispatch(lp);
        }
        // Contract: NaN bounds are Unsupported regardless of stage;
        // ±inf bounds (free vars) are handled via flip/split in to_standard.
        if let Some(j) = lp.lower.iter().position(|l| l.is_nan()) {
            return Err(SimplexError::Unsupported(format!(
                "NaN lower bound on variable {j}"
            )));
        }
        // Power-of-2 equilibration scaling around presolve+simplex
        // (module docs in crate::scaling).
        // Gate-rejected instances take a bit-identical clone and the identity
        // Scaling, so the unscaled path stays bit-identical. SimplexConfig
        // ::scaling is the lever (CLI: RUSOLVE_SCALE=1), off by default.
        let (slp, sc, applied) = if self.config.scaling {
            crate::scaling::scale_lp(lp)
        } else {
            (
                lp.clone(),
                crate::scaling::Scaling::identity(lp.num_row, lp.num_col),
                false,
            )
        };
        match crate::presolve::run(
            &slp,
            self.config.presolve_passes,
            self.config.presolve_implied_eq,
            self.config.presolve_parcol,
            self.config.presolve_agg,
            self.config.presolve_parcol_dom,
            self.config.presolve_vub_eq,
        ) {
            crate::presolve::Presolved::Terminal { status } => Ok(SimplexSolution {
                status,
                objective: match status {
                    SimplexStatus::Unbounded => f64::NEG_INFINITY,
                    _ => f64::NAN,
                },
                x: vec![],
                phase1_iterations: 0,
                phase2_iterations: 0,
                bland_triggered: false,
                // Presolve-ON contract: no duals, no rays — the
                // verdict never touches the caller's full row space.
                ..Default::default()
            }),
            crate::presolve::Presolved::Reduced {
                lp: reduced,
                live,
                postsolve,
            } => {
                // Diagnostics (env-gated `RUSOLVE_LP_LOG`, counters-only,
                // default off): presolve reduction size, original -> reduced
                // rows/cols/nonzeros.
                if std::env::var_os("RUSOLVE_LP_LOG").is_some() {
                    eprintln!(
                        "rusolve presolve: rows {}->{} cols {}->{} nnz {}->{}",
                        lp.num_row,
                        reduced.num_row,
                        lp.num_col,
                        reduced.num_col,
                        lp.a.col_ptr()[lp.num_col],
                        reduced.a.col_ptr()[reduced.num_col],
                    );
                }
                let sol = self.solve_dispatch(&reduced)?;
                if sol.status != SimplexStatus::Optimal {
                    return Ok(SimplexSolution {
                        status: sol.status,
                        objective: sol.objective,
                        x: vec![],
                        phase1_iterations: sol.phase1_iterations,
                        phase2_iterations: sol.phase2_iterations,
                        bland_triggered: sol.bland_triggered,
                        // A ray/dual computed on the REDUCED LP is NOT a
                        // certificate for the caller's LP (rows may be gone):
                        // report none rather than wrong ones.
                        ..Default::default()
                    });
                }
                let (mut x, objective) =
                    crate::presolve::restore(&postsolve, &live, &sol.x, sol.objective);
                if applied {
                    // x is in the SCALED caller space; map back to the
                    // original variable space (x_j = c_j x'_j, exact: powers
                    // of 2) and re-report the objective from the ORIGINAL
                    // costs at the unscaled x — the same convention solve_raw
                    // uses, so the value is invariant under scaling.
                    crate::scaling::unscale_x(&sc, &mut x);
                    let objective: f64 = lp.objective.iter().zip(&x).map(|(c, xv)| c * xv).sum();
                    return Ok(SimplexSolution {
                        status: sol.status,
                        objective,
                        x,
                        phase1_iterations: sol.phase1_iterations,
                        phase2_iterations: sol.phase2_iterations,
                        bland_triggered: sol.bland_triggered,
                        // Presolve-ON contract: sol's duals belong to
                        // the reduced+scaled LP — NOT postsolved, NOT reported.
                        ..Default::default()
                    });
                }
                Ok(SimplexSolution {
                    status: sol.status,
                    objective,
                    x,
                    phase1_iterations: sol.phase1_iterations,
                    phase2_iterations: sol.phase2_iterations,
                    bland_triggered: sol.bland_triggered,
                    // Honest reporting: the reduction stage is not pivots.
                    // Presolve-ON contract: reduced-LP duals are
                    // NOT postsolved — duals_valid stays false.
                    ..Default::default()
                })
            }
        }
    }

    /// Dual-engine lever dispatch: with
    /// `config.dual` on, try the dual engine first; `Ok(None)` is its
    /// NumericalDistress signal (e.g. a Singular refactorization with no
    /// artificial pool to repair) and re-runs the primal solver on the same
    /// LP — primal stays ground truth for all A/Bs. Lever off is the exact
    /// primal-only path (bit-identical).
    fn solve_dispatch(&self, lp: &Lp) -> Result<SimplexSolution, SimplexError> {
        let sol = if !self.config.dual
            && self.config.p1_route_dual
            // The dirty-exit rung is a PRIMAL re-solve (stability-armed tie
            // preference); routing must not hijack it into a dual retry.
            && !self.config.p1_rung_used
            && let Some(sol) = self.route_cold_phase1_dual(lp)?
        {
            sol
        } else if self.config.dual
            && let Some(sol) = crate::dual::solve_dual_raw(lp, &self.config, self.deadline)?
        {
            sol
        } else {
            self.solve_raw(lp)?
        };
        // Dirty-phase-1-exit rung (the optimal-verdict ladder pattern, one level
        // down): the phase-1 exit gate marked the point numerically unusable
        // (pricing "no improver" on an out-of-range basis). Re-solve the
        // SAME LP once with the exact-tie stability preference armed from
        // iteration 0 — the measured fix for e226 prim+pre (Optimal
        // 2,363+155; the burst-armed default reacts too late because the
        // detonation starts on healthy pivots steering INTO the degenerate
        // boom). Only dirty exits pay the rung; a tuff-shaped LP never
        // reaches it (clean phase-1 exits), so its measured preference
        // sensitivity is untouched. The rung's own dirty exit returns Limit
        // (`p1_rung_used`), never recursion; adopt only a clean Optimal.
        if sol.dirty_p1_exit && !self.config.p1_rung_used {
            let rung_cfg = SimplexConfig {
                p1_tie_pref_always: true,
                p1_rung_used: true,
                ..self.config
            };
            let rung = SimplexSolver {
                config: rung_cfg,
                deadline: self.deadline,
                cleanup_done: self.cleanup_done,
            };
            match rung.solve_dispatch(lp) {
                Ok(r) if r.status == SimplexStatus::Optimal && !r.dirty_p1_exit => {
                    if std::env::var("RUSOLVE_DEBUG").is_ok() {
                        eprintln!(
                            "rusolve simplex: dirty phase-1 exit rung adopted (obj {})",
                            r.objective
                        );
                    }
                    return Ok(r);
                }
                Ok(_) | Err(_) => {}
            }
        }
        Ok(sol)
    }

    /// Cold-start routing policy: route a COLD primal-engine start that is primal
    /// infeasible to the dual engine (see [`SimplexConfig::p1_route_dual`]).
    ///
    /// Detector: the dual form's slack basis start (`to_dual_standard` +
    /// `initial_basis` + `initial_x_b` — a ±1 diagonal, so no factorization
    /// is paid) is profiled with [`crate::dual::infeas_profile`]; any
    /// violation above the 1e-7 feasibility tolerance means the primal
    /// engine's cold solve would pay its artificial-basis phase 1, which on
    /// the blowup class of our 47 relaxed MIPLIB LPs grinds 20–90k iterations
    /// on LPs the dual solves in 0.3–3k (measured on our benchmark set). A feasible start returns `Ok(None)` and the
    /// caller keeps the primal path. The built [`crate::dual::DualStandard`]
    /// is handed to [`crate::dual::solve_dual_prepared`] so the routed
    /// solve pays the O(nnz) build exactly once; its NumericalDistress exit
    /// (`None`) propagates as "did not route", landing the caller's
    /// `solve_raw` cold-primal fallback.
    fn route_cold_phase1_dual(&self, lp: &Lp) -> Result<Option<SimplexSolution>, SimplexError> {
        let mut st = crate::dual::to_dual_standard(lp, &self.config)?;
        let (basis, diag) = crate::dual::initial_basis(lp, &st);
        let x_b = crate::dual::initial_x_b(&st, &basis, &diag);
        let (_, violations) = crate::dual::infeas_profile(&st, &basis, &x_b);
        if violations == 0 {
            return Ok(None);
        }
        crate::dual::solve_dual_prepared(lp, &mut st, &self.config, self.deadline)
    }

    /// Sparse front door. With
    /// `config.sparse` AND presolve/dual off, runs the CSC-only `Standard`
    /// path (`to_standard_sparse` + `gather_col` column access): no row-form
    /// LP is materialized. Lever off, or presolve/dual on (dense working
    /// state), = `to_lp()` and delegate to
    /// [`SimplexSolver::solve`] — the conversion performs no fp arithmetic,
    /// which keeps the Netlib A/B spine honest while the lever lands.
    ///
    /// LOUD SCOPE NOTE: on the CSC-only path the
    /// Optimal solution reports empty `row_duals`/`reduced_costs`
    /// (`duals_valid == false`) and `basis: None` — the caller-space dual
    /// mapping and the [`Basis`] row-status contract are defined on the
    /// single-sense/one-slack-per-row [`Lp`] frame and do not yet cover
    /// [`SparseLp`] ranged rows. Warm starts (`solve_with`/`reoptimize`) are
    /// [`Lp`]-frame features; feed a [`Basis`] through [`Self::solve_with`]
    /// on the [`Lp`] until sparse export is implemented.
    pub fn solve_sparse(&self, lp: &SparseLp) -> Result<SimplexSolution, SimplexError> {
        if !self.config.sparse
            || self.config.presolve
            || self.config.dual
            // The routing decision lives in solve_dispatch (Lp frame);
            // with the lever on, take the (bit-identical) `Lp` hop so it
            // cannot be silently bypassed by the CSC-only path.
            || self.config.p1_route_dual
        {
            return self.solve(&lp.to_lp());
        }
        self.solve_sparse_raw(lp)
    }

    /// CSC-only Standard construction + the shared two-phase engine:
    /// `to_standard_sparse` builds the column set straight into the CSC view
    /// and every dense-column consumer in the solve path runs on
    /// `gather_col`-filled workspaces — no dense matrix is ever
    /// materialized. The RUSOLVE_DEBUG residual audit materializes a row-form
    /// [`Lp`] for the caller-space half (env-gated, correctness surface, not
    /// hot).
    fn solve_sparse_raw(&self, lp: &SparseLp) -> Result<SimplexSolution, SimplexError> {
        let st = to_standard_sparse(lp)?;
        let audit = if std::env::var("RUSOLVE_DEBUG").is_ok() {
            Some(lp.to_lp())
        } else {
            None
        };
        self.solve_std(st, &lp.objective, lp.maximize, audit.as_ref())
    }

    fn solve_raw(&self, lp: &Lp) -> Result<SimplexSolution, SimplexError> {
        let st = to_standard(lp)?;
        self.solve_std(st, &lp.objective, lp.maximize, Some(lp))
    }

    /// The shared two-phase engine over a built `Standard`: used by the dense
    /// entry ([`Self::solve_raw`]) and the CSC-only entry
    /// ([`Self::solve_sparse`]). `objective`/`maximize` are the caller-space LP
    /// data needed for solution reporting (identical values whichever front
    /// door built `st` — the sparse-input conversions perform no fp arithmetic on them);
    /// `caller` is the caller-space [`Lp`]: always `Some` on the `Lp` path
    /// (caller-space duals, rays, and the basis export are all defined on the
    /// single-sense/one-slack-per-row `Lp` frame), and the [`Lp`] for the
    /// env-gated RUSOLVE_DEBUG residual audit. The sparse path passes `None`
    /// unless the audit is on: the [`Basis`] row-status contract and the caller-space
    /// dual mapping do not yet cover [`SparseLp`]'s ranged rows (two slacks,
    /// no single sense), so the sparse Optimal solution honestly reports
    /// empty duals and `basis: None` — "no information", never wrong
    /// numbers (same convention as the presolve-ON scope notes). Sparse
    /// dual/basis export is not implemented yet.
    fn solve_std(
        &self,
        mut st: Standard,
        objective: &[f64],
        maximize: bool,
        caller: Option<&Lp>,
    ) -> Result<SimplexSolution, SimplexError> {
        let m = st.m;
        let total = st.ncol();
        let art_start = total - m;

        // Initial basis: artificials, B0 = I.
        let mut basis: Vec<usize> = (art_start..total).collect();
        let mut x_b: Vec<f64> = st.b.clone();
        let mut at_ub = vec![false; total];

        // Initial basis = the artificial identity: each basis column is a
        // unit vector, so build_sparse factors it without any dense
        // materialization.
        let basis_cols: Vec<Vec<(usize, f64)>> = (0..m).map(|i| st.sparse_col(basis[i])).collect();
        let mut pfi = ProductFormFactorization::build_sparse(basis_cols, self.config.pfi)?;

        // Phase 1 costs: artificials 1, everything else 0.
        let mut cost1 = vec![0.0; total];
        cost1[art_start..total].iter_mut().for_each(|c| *c = 1.0);

        let mut bland = false;
        let mut randomized = false;
        // Steepest-edge weights live in phase 2 ONLY (bench-measured:
        // maintaining them through phase 1's thousands of Dantzig pivots
        // costs ~2 BTRANs+sweep per pivot — ~60 Gflops wasted on scrs8 —
        // while a one-time exact FTRAN init at the phase boundary costs
        // <0.5 Gflops). The snapshot/fallback machinery below still guards
        // the phase-2 SH path.
        let mut sh_w: Vec<f64> = Vec::new();
        let mut exit1 = PhaseExit::default();
        let (it1, outcome) = self.run_phase(
            &st,
            &mut pfi,
            &mut basis,
            &mut x_b,
            &mut at_ub,
            &cost1,
            true,
            &mut bland,
            &mut randomized,
            &mut sh_w,
            false, // phase 1: Dantzig-priced, weights neither priced nor kept
            0,
            &mut exit1,
        )?;
        if outcome == PhaseOutcome::Limit {
            return Ok(limit_solution(
                &st, objective, &basis, &x_b, &at_ub, it1, 0, bland,
            ));
        }
        if outcome == PhaseOutcome::Detonated {
            // Honest Limit, marked dirty ⇒ solve_dispatch re-solves
            // once with the rung's unconditional tie preference (the base
            // pass's basis is numerically garbage; the rung starts fresh).
            let mut lim = limit_solution(&st, objective, &basis, &x_b, &at_ub, it1, 0, bland);
            lim.dirty_p1_exit = true;
            return Ok(lim);
        }
        if outcome == PhaseOutcome::NumericalDistress {
            // Defensive: distress is emitted only by the SH phase-2 path.
            return Ok(limit_solution(
                &st, objective, &basis, &x_b, &at_ub, it1, 0, bland,
            ));
        }
        if outcome == PhaseOutcome::Unbounded {
            // Phase 1's artificial sum is bounded below by 0: a ray that
            // survived the fresh-factor verification is a numerical artifact;
            // the honest reading is infeasibility (p1_obj below concurs).
            return Ok(SimplexSolution {
                status: SimplexStatus::Infeasible,
                objective: f64::NAN,
                x: vec![],
                phase1_iterations: it1,
                phase2_iterations: 0,
                bland_triggered: bland,
                // No dual_ray here: phase 1 exited on a RAY, not on its
                // optimum, so exit1.pi carries no valid Farkas certificate.
                ..Default::default()
            });
        }
        let p1_obj: f64 = basis
            .iter()
            .zip(&x_b)
            .map(|(&v, &xb)| {
                if st.kind[v] == Kind::Artificial {
                    xb
                } else {
                    0.0
                }
            })
            .sum();
        // Relaxed phase-1 exit gate: pricing-optimality only proves anything
        // about a basis whose values are IN RANGE. A relaxed escape step may
        // end phase 1 with real structural bound violations — the honest
        // verdict is then Limit (those runs sit at the 100k cap today;
        // never a wrong Optimal/Infeasible off an out-of-range basis).
        let p1_viol: f64 = basis
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                if st.kind[v] == Kind::Artificial {
                    0.0
                } else {
                    (-x_b[i]).max(x_b[i] - st.ub[v]).max(0.0)
                }
            })
            .sum();
        if p1_viol > 1e-7 {
            if std::env::var("RUSOLVE_DEBUG").is_ok() {
                eprintln!("rusolve simplex: phase-1 exit with bound violations {p1_viol} ⇒ Limit");
            }
            let mut lim = limit_solution(&st, objective, &basis, &x_b, &at_ub, it1, 0, bland);
            lim.dirty_p1_exit = true; // arms the solve_dispatch re-solve rung
            return Ok(lim);
        }
        if std::env::var("RUSOLVE_DEBUG").is_ok() {
            eprintln!("rusolve simplex: phase-1 exit p1_obj={p1_obj} iters={it1}");
        }
        if p1_obj > 1e-7 {
            if std::env::var("RUSOLVE_DEBUG").is_ok() {
                eprintln!(
                    "rusolve simplex: INFEASIBLE verdict basis={:?} x_b={:?}",
                    basis, x_b
                );
            }
            return Ok(SimplexSolution {
                status: SimplexStatus::Infeasible,
                objective: f64::NAN,
                x: vec![],
                phase1_iterations: it1,
                phase2_iterations: 0,
                bland_triggered: bland,
                // Phase-1 exit duals (w > 0) → Farkas candidate, kept
                // only if it verifies against the original rows (dense
                // caller LP only — sparse path reports honest `None`).
                dual_ray: caller.and_then(|lp| attempt_farkas(lp, &exit1.pi)),
                ..Default::default()
            });
        }

        // Phase 2: pin every artificial to zero. Phase 1 can leave artificials
        // basic at a degenerate 0, and with ub = inf the phase-2 ratio test
        // caps nothing on their increasing side (`st.ub[basis[i]].is_finite()`
        // at the bbar<0 branch) — such an artificial then INFLATES across
        // phase-2 pivots, "satisfies" its standard row, and the matching
        // caller row is violated by exactly its value (boeing2 false Optimal:
        // obj -320.37 vs true -315.02; row 124 violated by 10.84 = the
        // artificial's basic value; the internal-state audit was clean
        // because the std system genuinely WAS satisfied - by the artificial).
        // ub=0 makes pricing skip it (fixed var) and the ratio test cap its
        // increase at 0 ⇒ a degenerate pivot drives it out, or it stays
        // pinned at 0 as the valid redundant equality it then is. Zero the
        // phase-1 epsilon first so the (ub - x_b) cap cannot come out tiny
        // negative.
        for v in art_start..total {
            st.ub[v] = 0.0;
        }
        for (i, &v) in basis.iter().enumerate() {
            if v >= art_start && x_b[i].abs() < 1e-7 {
                x_b[i] = 0.0;
            }
        }

        // Phase 2 costs: true objective (sign-normalized to minimization),
        // taken per structural column (flip/split columns carry signed costs).
        let sign = if maximize { -1.0 } else { 1.0 };
        let mut cost2 = vec![0.0; total];
        for (d, &cc) in cost2.iter_mut().zip(&st.cost) {
            *d = cc * sign;
        }
        // Fresh anti-cycling state at the phase boundary: the sticky
        // Bland/randomized flags were engaged by PHASE 1's stall history on a
        // different objective — carrying them in silently masked phase-2
        // pricing entirely (Bland is checked before the SH rate: bit-identical
        // SH-vs-Dantzig benches until this reset). If phase 2 degenerates on
        // its own, the stall counter re-engages Bland after bland_after.
        // Gated to sh_pricing (DEFAULT ON with artificial deletion: the bench said
        // otherwise — -37% iters, wall tie, see SimplexConfig::default).
        let bland_boundary = bland;
        let randomized_boundary = randomized;
        if self.config.sh_pricing {
            bland = false;
            randomized = false;
        }
        // Phase-2 exit-state snapshot for the SH numerical-distress
        // fallback: SH pricing is exploratory — if it steers into a
        // rank-deficient basis (fake ray), restart phase 2 Dantzig-priced
        // from this exact state. By construction SH-mode can never be more
        // wrong than the baseline verdict.
        let snap_basis = basis.clone();
        let snap_xb = x_b.clone();
        let snap_ub = at_ub.clone();
        // Exact phase-1-exit factorization (LU + eta chain + basis_cols):
        // scsd1's phase-1-exit basis REFACTORIZES singular (col 49, pivot 0)
        // yet its eta chain solves fine — so the restore must be the pfi
        // itself, not a fresh refactorize of the restored basis.
        let snap_pfi = if self.config.sh_pricing {
            Some(pfi.clone())
        } else {
            None
        };
        let mut exit2 = PhaseExit::default();
        let (mut it2, mut outcome) = self.run_phase(
            &st,
            &mut pfi,
            &mut basis,
            &mut x_b,
            &mut at_ub,
            &cost2,
            false,
            &mut bland,
            &mut randomized,
            &mut sh_w,
            self.config.sh_pricing,
            it1,
            &mut exit2,
        )?;
        if outcome == PhaseOutcome::NumericalDistress {
            if std::env::var("RUSOLVE_DEBUG").is_ok() {
                eprintln!("rusolve simplex: SH distress at it={it2} ⇒ Dantzig phase-2 restart");
            }
            basis.copy_from_slice(&snap_basis);
            x_b.copy_from_slice(&snap_xb);
            at_ub.copy_from_slice(&snap_ub);
            if let Some(s) = snap_pfi {
                pfi = s; // byte-identical to the baseline's phase-2 start
            }
            bland = bland_boundary;
            randomized = randomized_boundary;
            sh_w.clear();
            exit2 = PhaseExit::default(); // restart recomputes its own exit data
            let (i2, o2) = self.run_phase(
                &st,
                &mut pfi,
                &mut basis,
                &mut x_b,
                &mut at_ub,
                &cost2,
                false,
                &mut bland,
                &mut randomized,
                &mut sh_w,
                false,
                it1,
                &mut exit2,
            )?;
            it2 = i2;
            outcome = o2;
        }
        if outcome == PhaseOutcome::Unbounded {
            // The verified direction, mapped to caller column space
            // (dense caller LP only — sparse path reports honest `None`).
            let primal_ray = match (caller, exit2.ray.as_ref()) {
                (Some(lp), Some((q, sgn, bbar))) => primal_ray_x(lp, &st, &basis, *q, *sgn, bbar),
                _ => None,
            };
            return Ok(SimplexSolution {
                status: SimplexStatus::Unbounded,
                objective: f64::NEG_INFINITY,
                x: vec![],
                phase1_iterations: it1,
                phase2_iterations: it2 - it1,
                bland_triggered: bland,
                primal_ray,
                ..Default::default()
            });
        }
        if outcome == PhaseOutcome::Limit || outcome == PhaseOutcome::Detonated {
            // Detonated cannot fire on a phase-2 run (allow_artificial=false);
            // treated as Limit defensively.
            return Ok(limit_solution(
                &st,
                objective,
                &basis,
                &x_b,
                &at_ub,
                it1,
                it2 - it1,
                bland,
            ));
        }

        // Extract caller-space solution (undoes shift/flip/split).
        // Final re-anchor: solve x_b = B⁻¹b on fresh factors so the reported
        // solution satisfies the std system to factor accuracy. The
        // incremental −θ·b̄ updates drift on ill-conditioned bases (scsd6:
        // Optimal carried a 1e-4 row residual and 1.2e-3 obj error while the
        // BASIS was correct — pricing saw the true optimum, the reported
        // values were stale). If the basis refactorizes singular (scsd1
        // lesson: its phase-1-exit basis does), the eta chain is the better
        // solver — keep x_b as-is.
        if pfi.refactorize().is_ok() {
            refresh_x_b(&mut pfi, &st, &basis, &at_ub, &mut x_b);
        }
        // FINAL probe (debug): audit the terminal state three ways — basic
        // bounds, std-system residual (including fixed-at-ub non-basics), and
        // caller-space row/bound violations. This audit is what caught the
        // boeing2 false Optimal (a basic artificial inflating past 0 to
        // satisfy its std row while the caller row stayed violated; see the
        // phase-boundary pin). Env-gated correctness surface (not hot), so
        // the residual sweeps gather dense columns on demand; the caller-space
        // half needs the dense LP (the sparse entry materializes one under
        // RUSOLVE_DEBUG only).
        if std::env::var("RUSOLVE_DEBUG").is_ok() {
            let mut worst_lo = (0.0f64, 0usize);
            let mut worst_up = (0.0f64, 0usize);
            for (i, &v) in basis.iter().enumerate() {
                if -x_b[i] > worst_lo.0 {
                    worst_lo = (-x_b[i], v);
                }
                let upv = x_b[i] - st.ub[v];
                if upv.is_finite() && upv > worst_up.0 {
                    worst_up = (upv, v);
                }
            }
            let mut in_b = vec![false; st.ncol()];
            for &v in basis.iter() {
                in_b[v] = true;
            }
            let mut colbuf = vec![0.0; m];
            let mut r: Vec<f64> = st.b.clone();
            for (i, &v) in basis.iter().enumerate() {
                st.gather_col(v, 1.0, &mut colbuf);
                for (ri, cv) in r.iter_mut().zip(&colbuf) {
                    *ri -= x_b[i] * cv;
                }
            }
            for (v, &ib) in in_b.iter().enumerate() {
                if !ib && at_ub[v] && st.ub[v].is_finite() {
                    st.gather_col(v, 1.0, &mut colbuf);
                    for (ri, cv) in r.iter_mut().zip(&colbuf) {
                        *ri -= st.ub[v] * cv;
                    }
                }
            }
            let (wr, wi) = r
                .into_iter()
                .enumerate()
                .fold((0.0f64, 0usize), |a, (i, x)| {
                    if x.abs() > a.0 { (x.abs(), i) } else { a }
                });
            let xl = extract_x(&st, &basis, &x_b, &at_ub);
            let mut worst = (0.0f64, 0usize, 'L');
            let mut wb = (0.0f64, 0usize);
            if let Some(lp) = caller {
                // Row activities off the CSC (ascending-j accumulation per
                // row — same order the dense loop had; zeros are ±0.0 no-ops).
                let (cptr, cridx, crval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
                let mut act = vec![0.0f64; lp.num_row];
                for j in 0..lp.num_col {
                    let xj = xl[j];
                    if xj != 0.0 {
                        for k in cptr[j]..cptr[j + 1] {
                            let v = crval[k];
                            if v != 0.0 {
                                act[cridx[k] as usize] += v * xj;
                            }
                        }
                    }
                }
                for (i, (&sense_i, &rhs_i)) in lp.senses.iter().zip(&lp.rhs).enumerate() {
                    let lhs = act[i];
                    let v = match sense_i {
                        'L' => lhs - rhs_i,
                        'G' => rhs_i - lhs,
                        _ => (lhs - rhs_i).abs(),
                    };
                    if v > worst.0 {
                        worst = (v, i, sense_i);
                    }
                }
                wb = xl.iter().enumerate().fold((0.0f64, 0usize), |a, (j, &xv)| {
                    let v = (lp.lower[j] - xv).max(xv - lp.upper[j]);
                    if v > a.0 { (v, j) } else { a }
                });
            }
            eprintln!(
                "rusolve FINAL: worst x_b<0 = {worst_lo:?}  worst x_b>ub = {worst_up:?}  \
                 worst |Bx_b-b| = {wr:.3e} at row {wi}  worst row violation {worst:?}  \
                 worst bound {wb:?}"
            );
        }
        let x = extract_x(&st, &basis, &x_b, &at_ub);
        let objective: f64 = objective.iter().zip(&x).map(|(c, xv)| c * xv).sum();
        // Caller-space duals from the exported phase-2 exit π (the same
        // π the final no-improver pricing pass used). The mapping is defined
        // on the caller [`Lp`] (single sense + one slack per row), so it
        // runs whenever `caller` is present: always on the `Lp` path; the
        // sparse path exports duals/basis only under RUSOLVE_DEBUG (honest
        // empty otherwise — see the `caller` doc on `solve_std`, and the
        // sparse dual/basis export not implemented yet).
        // Guard on the length: an unfilled exit (defensive; cannot happen on
        // this arm) reports duals_valid=false rather than guessing.
        let (row_duals, reduced_costs, duals_valid) = match (caller, exit2.pi.len() == m) {
            (Some(lp), true) => {
                let (rd, rc) = caller_duals(lp, sign, &exit2.pi);
                (rd, rc, true)
            }
            _ => (Vec::new(), Vec::new(), false),
        };
        let mut sol = SimplexSolution {
            status: SimplexStatus::Optimal,
            objective,
            x,
            phase1_iterations: it1,
            phase2_iterations: it2 - it1,
            bland_triggered: bland,
            dirty_p1_exit: false,
            row_duals,
            reduced_costs,
            duals_valid,
            primal_ray: None,
            dual_ray: None,
            basis: caller.map(|lp| export_basis(&st, lp, &basis, &at_ub)),
        };
        // ---- Defense: never report Optimal on an infeasible point ----
        // Phase 2's exit test is DUAL (no improver in pricing); its PRIMAL
        // invariant `0 ≤ x_B ≤ ub` is maintained only by the ratio test on
        // b̄ from the PFI factorization, and a numerically bad pivot breaks it
        // silently (repro: big-M job-shop node LPs). HiGHS guards the same exit by recomputing the primal from the
        // final factorization and, on violation, running a cleanup instead of
        // reporting Optimal (HEkkPrimal.cpp:586–601 + cleanup() :647–688,
        // MIT — mechanism only). Its cleanup is a dual-simplex run; the warm
        // dual engine here lives in `crate::dual` (separate module), so the
        // defense is a bounded re-solve ladder on the UNPERTURBED numerics:
        // one rung, SH pricing off (the repro evidence: RUSOLVE_NO_SH
        // removes the row violation), re-audited. On the clean path this
        // costs one O(nnz) row+bound sweep and never fires (Netlib/diff
        // corpus land at 1e-9..1e-13). The rung's own audit is suppressed
        // (`cleanup_done`) so the chain terminates; if the rung is still
        // dirty, return whichever point is closer to feasible at the SAME
        // objective — never a different objective, never worse.
        if !self.cleanup_done
            && let Some(lp) = caller
            && !primal_feasible(lp, &sol.x)
        {
            let (wb0, wr0) = primal_violations(lp, &sol.x);
            let rung_cfg = SimplexConfig {
                sh_pricing: false,
                ..self.config
            };
            let rung = SimplexSolver {
                config: rung_cfg,
                deadline: self.deadline,
                cleanup_done: true,
            };
            match rung.solve_raw(lp) {
                Ok(r) if r.status == SimplexStatus::Optimal => {
                    let (wb1, wr1) = primal_violations(lp, &r.x);
                    let clean = primal_feasible(lp, &r.x);
                    let otol = 1e-6 * sol.objective.abs().max(r.objective.abs()).max(1.0);
                    let obj_diff = (r.objective - sol.objective).abs();
                    let better_obj = (maximize && r.objective > sol.objective + otol)
                        || (!maximize && r.objective < sol.objective - otol);
                    let same_obj = obj_diff <= otol;
                    let less_dirty = wb1.max(wr1) < wb0.max(wr0);
                    let adopt = (clean && (same_obj || better_obj)) || (same_obj && less_dirty);
                    eprintln!(
                        "rusolve simplex: primal cleanup — Optimal x violated bounds \
                                 {wb0:.2e} / rows {wr0:.2e}; SH-off re-solve obj {:.9} (Δ {:.2e}) \
                                 viol {:.2e}/{:.2e} ⇒ {}",
                        r.objective,
                        obj_diff,
                        wb1,
                        wr1,
                        if adopt { "adopted" } else { "kept original" }
                    );
                    if adopt {
                        sol = r;
                    }
                }
                Ok(r) => eprintln!(
                    "rusolve simplex: primal cleanup — SH-off re-solve returned \
                             {:?}, kept original Optimal",
                    r.status
                ),
                Err(e) => eprintln!(
                    "rusolve simplex: primal cleanup — SH-off re-solve errored ({e:?}), \
                     kept original Optimal"
                ),
            }
            // Last line of defense: the ladder could not produce a feasible
            // point. Never hand an infeasible x to callers as Optimal (the
            // MIP solver would branch/prune on it); report a non-optimal status
            // that every caller already treats as "LP not solved".
            if !primal_feasible(lp, &sol.x) {
                eprintln!(
                    "rusolve simplex: primal cleanup — point still infeasible after the \
                     ladder ⇒ reporting IterationLimit instead of Optimal"
                );
                sol.status = SimplexStatus::IterationLimit;
            }
        }
        Ok(sol)
    }

    /// One phase's pivot loop. `allow_artificial` gates artificial entering
    /// (true in phase 1 only). Returns (total iterations so far, outcome).
    // State is threaded explicitly (borrow-split over pfi/basis/values);
    // a struct wrapper is a later refactor when more state accrues.
    #[allow(clippy::too_many_arguments)]
    fn run_phase(
        &self,
        st: &Standard,
        pfi: &mut ProductFormFactorization,
        basis: &mut [usize],
        x_b: &mut Vec<f64>,
        at_ub: &mut [bool],
        cost: &[f64],
        allow_artificial: bool,
        bland_triggered: &mut bool,
        randomized: &mut bool,
        sh_w: &mut Vec<f64>,
        sh_active: bool,
        start_iters: usize,
        exit: &mut PhaseExit,
    ) -> Result<(usize, PhaseOutcome), SimplexError> {
        let tol = self.config.pivot_tol;
        let ptol = self.config.pricing_tol;
        let m = st.m;
        let total = st.ncol();
        // Artificial DELETION convention (Bertsimas-Tsitsiklis phase-I
        // handling; Wolkowitz CO350: drop artificials as soon as nonbasic):
        // an artificial that leaves phase-1 pricing can never re-enter — it
        // left at value 0 (ub=inf ⇒ never leaves at ub), the basis stayed
        // primal feasible, so the restricted problem remains feasible and
        // phase 1 still reaches w=0. Makes the basic-artificial count
        // monotonically non-increasing: a well-founded secondary progress
        // measure phase-1 pricing never had. Basis REPAIR may still swap a
        // dead artificial back in (numerical emergency, not pricing); if it
        // leaves again it is re-marked.
        let mut dead: Vec<bool> = vec![false; total];
        // Phase-1 sub-tolerance ENTERING tie perturbation: built ONCE
        // when the stall escalation first fires (see trigger below); empty ⇒
        // off. Non-stalling phase-1 runs never build it ⇒ bit-identical
        // without it (Beale 1968; Chvátal 1983 §7.3).
        let mut p1_eps: Vec<f64> = Vec::new();
        // Phase-2 ENTERING tie perturbation: the phase-2 counterpart of
        // the phase-1 one, reusing p1_tie_eps verbatim. Built ONCE when phase 2
        // reaches the SAME rung as phase 1 (4×bland_after while Bland stalls); empty ⇒
        // off ⇒ non-stalling phase-2 runs stay bit-identical. Same safety
        // invariant: acceptance stays the UNPERTURBED improver test, eps <
        // ptol/2 only reorders the ranking key (SH d²/w or Dantzig |d|) —
        // never the objective bookkeeping, ratio test, or leaving choice.
        // Beale 1968 / Chvátal §7.3 / HiGHS costPerturbation lineage.
        let mut p2_eps: Vec<f64> = Vec::new();
        // Phase-1 re-perturbation state (see SimplexConfig::p1_re_perturb):
        // `p1_eps_seed` 0 = the original one-shot vector; `p1_recheck` is the
        // stall value at which the next freeze probe runs (usize::MAX until
        // eps engages); `p1_w_probe` is the (w, n_art) snapshot from the
        // last probe (frozen ⇒ bit-equal at the next one).
        let mut p1_eps_seed = 0u64;
        let mut p1_recheck = usize::MAX;
        let mut p1_w_probe: (f64, usize) = (f64::NAN, usize::MAX);
        // Randomized-LEAVING escape engaged by the frozen-w probe
        // (see SimplexConfig::p1_re_perturb); releases when w moves again.
        let mut p1_rand_leave = false;
        // HiGHS-style RELAXED phase-1 row choice (crossing structural bound
        // rows on a provable total-infeasibility descent; HiGHS:
        // simplex/HEkkPrimal.cpp:phase1ChooseRow). Engaged/released with
        // `p1_rand_leave` by the freeze probe; `_ever` pins the w-measure
        // switch (total_infeas instead of artificial_sum) for the probes
        // and the exit gate for the rest of the phase.
        let mut p1_relaxed = false;
        let mut p1_relaxed_ever = false;
        let mut p1_frozen_streak = 0u32;
        // Relaxed phase-1 stagnation bar: consecutive probe windows (100 stalls
        // each) without a 0.1% w drop; 400 ⇒ engage the relaxed row choice
        // (dano3_3 oscillates bit-wise without net descent). `p1_w_ref` is
        // the w at the last genuine-descent reset (NaN before the first
        // probe).
        let mut p1_stag_streak = 0u32;
        let mut p1_w_ref = f64::NAN;
        // Relaxed phase-1 burn: engagement is REVERSIBLE. The relaxed walk
        // exists to leave degenerate plateaus the exact ratio test cannot —
        // but on tuff (measured) it INCREASES total infeasibility (3.4e2 →
        // 8.5e2) and freezes; no escalation bar separates "permanent
        // freeze" (qap10: bit-frozen from probe ~17, walks to feasibility
        // once engaged at streak 200, w moves every probe after) from "long
        // natural freeze" (tuff: streak 746 inside a 748-window oscillation
        // with no genuine descent, which v2 escapes ALONE at 77,328 — a
        // bar of 200 derailed it, and 800 starved qap10 of runway: engaged
        // at stall 81,600 the walk cannot finish inside the 100k cap).
        // Burn discriminator (r4, qap10+tuff traces): NOT frozen-probe
        // counting. The walk walks OUT of a freeze through θ=0 degenerate
        // edges, so the strict w census is BIT-frozen while the walk is
        // alive and winning (qap10: 50 consecutive frozen probes under a
        // walk that won 3,399/3,400 with zero dMaxTheta=0 kills — burning
        // there reverted qap10's finished run at stall 26,500 and
        // re-capped it). What genuinely discriminates is the walk MAKING
        // INFEASIBILITY WORSE: tuff-np's walk pushes the census from 3.4e2
        // to 1.05e4 within 700 stalls of engagement (a 30× rise — burn as
        // EARLY as the census allows, 1.2×, so the revert lands while the
        // created violations are still small enough for the θ=0 repair to
        // undo), while qap10 rides flat at 20.0 and dano3_3 descends.
        // Burn ⇒ revert to the exact walk, lock BOTH engagement paths,
        // until genuine descent re-arms. `p1_w_burn_ref` is the
        // descent min of the SAME w census over the whole run (total_infeas
        // once the relaxed step has ever engaged — comparing against the
        // pre-engagement artificial-sum min would mix measures); it is
        // NOT re-baselined at engagement (that reset handed the walk
        // ~20% w headroom before the burn check activated, and dano3_3
        // spent exactly that window destroying feasibility).
        let mut p1_w_burn_ref = f64::NAN;
        let mut p1_relaxed_burned = false;
        // Relaxed phase-1 burn ROLLBACK. The 100-stall probe granularity means a
        // burn is detected only AFTER the relaxed walk has already moved the
        // vertex (brandy: census 9.58e2 → 3.6e3 in ONE probe window; the
        // revert-to-exact then crawls on the damaged vertex, and the burn
        // probe's own w-move even releases the v2 random-tie walk mid-
        // recovery). Burn therefore also UNDOES: basis/x_b/at_ub and the
        // probe tuple snapshotted at engagement are restored, with fresh
        // factors + drift-corrected values (the repair pass re-scans x_b
        // each iteration, so restoring it clears its state). Sound in phase
        // 1 specifically: artificial costs are CONSTANT, so pricing is a
        // pure function of the restored basis — no objective-side state to
        // rewind. The rng is not rewound: the post-burn trajectory is the
        // v2 walk from the engagement vertex, not a bit-identical replay.
        // Taken at every fresh engagement, consumed at every burn (the
        // burned lock prevents re-engagement until genuine descent).
        let mut p1_snap: Option<P1EngageSnapshot> = None;
        // One fresh-factor verification per phase-1 pricing exit
        // attempt; cleared by every real pivot so a later exit is verified
        // again (see the `best`-is-None arm below).
        let mut p1_exit_verified = false;
        let mut p1_events: Vec<(f64, f64)> = Vec::new();
        // Debug counters (RUSOLVE_P1TRACE): relaxed-step census.
        let (mut rx_tried, mut rx_won, mut rx_s0pos, mut rx_dmax0) = (0u64, 0u64, 0u64, 0u64);
        let (mut rx_unstable, mut rx_pivots) = (0u32, 0u32);
        // Sub-1e-7 leaving-pivot history (last 10 pivots, LSB =
        // newest). 3+ set bits arms the exact-tie stability preference —
        // e226's detonation window (9/50) vs tuff's Optimal path (≤ 2/50).
        let mut stab_hist: u16 = 0;
        // Relaxed phase-1 REPAIR: evict a structural basic row that is out of
        // range (a relaxed step crossed its bound) via a θ=0 pivot whose
        // entering column has support in its tableau row; the recorded
        // bound lie (≤ the violation) is re-anchored immediately, so the
        // residual flows into the artificials — the same absorption the
        // drift correction does. Bounded failures ⇒ give up (exit gate
        // then reports Limit, never a wrong verdict).
        let mut repair_head: Option<usize>;
        let mut repair_fails = 0u32;
        let mut p1_repair_off = false;
        let mut iters = start_iters;
        let mut stall = 0usize;
        // Bounded Bland episode state. `bland_triggered` (threaded out to
        // the caller/trace) now records "has engaged at least once"; pricing
        // follows the per-episode countdown below. The old engage was
        // TERMINAL — once Bland engaged it never released, so degenerate
        // packing node LPs crawled Bland's θ=0 walk to the 100k cap while
        // RUSOLVE_BLAND_AFTER=1e9 solved them in 174-313 p2 iters
        // (benches/data/repros/set_packing_s*_node_stall: the bug was in
        // EXITING Bland, not in pricing). Release after
        // `bland_clear` >= bland_after non-degenerate pivots (the pocket
        // demonstrably broke), on budget exhaustion, or when the randomized
        // deviation fires; rungs (tie-perturbation eps) are sticky-once and keep
        // working after release. Re-entry backoff: the stall threshold
        // doubles per episode. HiGHS releases perturbation/bound-switching
        // the same way (public docs: maxPivot budget or a non-degenerate
        // streak).
        let mut bland_left = 0usize;
        let mut bland_clear = 0usize;
        let mut bland_episodes = 0usize;
        let mut bland_sticky = false;
        // One fresh-factor verification allowed per unbounded verdict (the
        // guard consumes itself on the next real pivot).
        let mut ray_verified = false;
        // Deterministic RNG for the randomized pricing fallback: fixed seed +
        // iteration-driven draws ⇒ reproducible runs (benchmark integrity).
        let mut rng =
            crate::rng::XorShift64::new(0x9E37_79B9_7F4A_7C15 ^ ((m as u64) << 32 | total as u64));
        // TEMPORARY diagnostic (env-gated, phase 1 only): see P1Trace below.
        let p1trace = allow_artificial && std::env::var("RUSOLVE_P1TRACE").is_ok();
        let mut tr = P1Trace {
            w_true: f64::NAN,
            ..Default::default()
        };

        let mut pi = vec![0.0; m];
        let mut bbar = vec![0.0; m];
        let mut cb = vec![0.0; m];
        let mut in_basis = vec![false; total];
        // Column-access workspaces (persistent-workspace reuse convention):
        // `aq` holds the sign-scaled FTRAN input, `qcol` the RAW entering
        // column for PFI updates and steepest-edge weight solves. gather_col
        // fully overwrites (fill + scatter) before every read ⇒ reuse is
        // invisible to the fp operation sequence.
        let mut aq = vec![0.0; m];
        let mut qcol = vec![0.0; m];
        // ---- Steepest-edge weights (SH path only) ----
        // Invariant w_j = 1 + ‖B⁻¹a_j‖² (full-tableau column norm, HiGHS
        // convention ⇒ w ≥ 1, no division blow-up). EXACT init at the
        // phase-2 start: one FTRAN per column on the phase-1-exit basis
        // (<0.5 Gflops on Netlib scale). Phase 1 neither prices with nor
        // maintains weights — bench-measured: per-pivot maintenance across
        // phase 1's Dantzig pivots (~60 Gflops on scrs8) swamps the whole
        // benefit. Algebra: Goldfarb & Reid 1977; Forrest & Goldfarb 1992.
        if sh_active && sh_w.is_empty() {
            let mut t = vec![0.0; m];
            sh_w.reserve(total);
            for v in 0..total {
                // The init is one FTRAN per standard-form column ahead
                // of the pivot loop's own deadline check — tens of seconds on
                // k1mushroom (16 419 rows, 1.7M nonzeros), which overran a
                // 60 s limit by 49 s. Poll
                // the deadline per column (an FTRAN dwarfs the vdso read)
                // and leave the weights empty, so a re-entry re-inits them.
                if let Some(deadline) = self.deadline
                    && std::time::Instant::now() >= deadline
                {
                    sh_w.clear();
                    return Ok((iters, PhaseOutcome::Limit));
                }
                st.gather_col(v, 1.0, &mut qcol);
                pfi.solve_b(&qcol, &mut t);
                sh_w.push(1.0 + t.iter().map(|x| x * x).sum::<f64>());
            }
        }
        let w = sh_w;
        let mut pi_pos = vec![0.0; m];
        let mut e_pos = vec![0.0; m];
        let mut mu = vec![0.0; m];
        let mut alpha_raw = vec![0.0; m];
        let mut weights_stale = false;
        let mut sh_report_left = 8u32;
        // Incremental reduced costs (Dantzig 1963 revised-simplex update;
        // HiGHS HEkkPrimal updates dual values from the pivotal row, Hall &
        // McKinnon 2005): d_j is recomputed from a fresh BTRAN only when
        // invalid/stale, otherwise d_j -= (d_q/alpha_rq) alpha_rj with the
        // pivotal row alpha_r = rho^T A priced through the row-wise copy of A
        // (touching only columns hit by rho's nonzeros).
        let (rptr, rcol, rval) = build_row_view(st);
        let mut dvec = vec![0.0; total];
        let mut d_ok = false;
        let mut trouble_check = true; // pivot column/row cross-check enabled
        let mut d_epoch = 0u64;
        // Wall time spent in DISTRESS refactorizations — rebuilds forced
        // by the FTRAN/BTRAN cross-check mismatch, an update refusal, or a
        // repair-pivot re-anchor below (NOT the healthy eta-limit/tick
        // cadence, which never touches those arms). The mik restart-round-2
        // crawl (measured: 172 of 173 dense rebuilds
        // are repair-cascade re-anchors, `epoch` 29→201 sticky-dense at
        // m=1 040, `etas=1` at every event — a rebuild per committed pivot,
        // each ~0.35 s, 98% of wall in DenseLu::factor) hits a distress arm
        // roughly every 10th pivot: too sparse for any consecutive/window
        // count (interleaved clean pivots reset or cancel every one tried),
        // but its REFACTOR WALL dwarfs a healthy solve's — distress events
        // are ms-scale on any well-conditioned basis, so ≥2 s of them is a
        // dead basis by orders of magnitude on either side.
        let mut distress_secs = 0.0f64;
        let mut distress_events = 0u32;
        let mut d_age = 0usize;
        let mut row_alpha = vec![0.0; total];
        let mut row_mark = vec![false; total];
        let mut touched: Vec<usize> = Vec::new();
        // Nonzero positions of the FTRAN column b̄ (ascending): the ratio test,
        // x_B update and cost dot only visit these.
        let mut nzb: Vec<usize> = Vec::with_capacity(m);
        let mut cand_bits = vec![0u64; total.div_ceil(64)];

        loop {
            iters += 1;
            if iters - start_iters > self.config.max_iterations {
                if p1trace {
                    tr.emit(iters, *bland_triggered, !p1_eps.is_empty());
                }
                return Ok((iters, PhaseOutcome::Limit));
            }
            // In-loop deadline check (HiGHS `HEkk::bailout` shape — it
            // too checks every iteration; see `SimplexSolver::set_deadline`).
            // `Instant::now()` is a vdso read, noise against an FTRAN+BTRAN
            // pivot, and `None` (the default) short-circuits to nothing.
            if let Some(deadline) = self.deadline
                && std::time::Instant::now() >= deadline
            {
                if p1trace {
                    tr.emit(iters, *bland_triggered, !p1_eps.is_empty());
                }
                return Ok((iters, PhaseOutcome::Limit));
            }
            if allow_artificial {
                // Detonation circuit breaker: once a basic value has
                // left the credible range the basis is gone (measured e226
                // prim+pre: x_b → 6e19, then ~500 s of refactor-storm
                // crawling to an exit-gate Limit — the storm (per-pivot
                // refactorizes on near-singular etas) begins long before
                // 1e12, so the pivot-instability burst that precedes every
                // measured detonation (≥ 3 sub-1e-7 leavers in the last 10;
                // tuff's Optimal path never exceeds 2 in 50) drops the bar
                // to 1e8. Read-only check on the previous pivot's values
                // (one pivot of granularity); NaN-proof form (a `>` test
                // would sail past NaN). The rung's unconditional tie
                // preference re-solves e226 to Optimal 2,363+155 fresh.
                let mx = x_b.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                let armed = stab_hist.count_ones() >= 3;
                if !mx.is_finite() || mx > 1e12 || (armed && mx > 1e8) {
                    return Ok((iters, PhaseOutcome::Detonated));
                }
            }
            for b in in_basis.iter_mut() {
                *b = false;
            }
            for &v in basis.iter() {
                in_basis[v] = true;
            }
            if sh_active && weights_stale {
                // Basis swapped outside the pivot math (phase-1 repair): the
                // current basis is NOT logical, so init values would be wrong
                // — recompute weights exactly via fresh FTRANs instead
                // (rare path; m is small whenever repair fires).
                for (v, wv) in w.iter_mut().enumerate() {
                    let mut t = vec![0.0; m];
                    st.gather_col(v, 1.0, &mut qcol);
                    pfi.solve_b(&qcol, &mut t);
                    *wv = 1.0 + t.iter().map(|x| x * x).sum::<f64>();
                }
                weights_stale = false;
            }

            // π: solve Bᵀ π = c_B
            let mut pi_fresh = false;
            if !d_ok || (pfi.factor_epoch() != d_epoch && d_age >= DJ_MIN_AGE) {
                for (i, &v) in basis.iter().enumerate() {
                    cb[i] = cost[v];
                }
                pfi.solve_bt(&cb, &mut pi);
                for v in 0..total {
                    dvec[v] = if in_basis[v] {
                        0.0
                    } else {
                        cost[v] - sparse_dot(&pi, st, v)
                    };
                }
                cand_bits.fill(0);
                for v in 0..total {
                    if !in_basis[v]
                        && is_cand(st, allow_artificial, dead[v], at_ub[v], dvec[v], ptol, v)
                    {
                        cand_bits[v >> 6] |= 1u64 << (v & 63);
                    }
                }
                d_ok = true;
                d_epoch = pfi.factor_epoch();
                d_age = 0;
                pi_fresh = true;
            }
            d_age += 1;

            // ---- Pricing (Dantzig → bounded Bland episode → randomized escape) ----
            // Multi-stage fallback on the *decaying* stall counter (below):
            //  1. bland_after: a Bland EPISODE (smallest-index entering). Its
            //     anti-cycling proof holds in exact arithmetic and mostly
            //     holds here; leaving ties break to lowest variable index.
            //     The episode is bounded and releases (see state comment
            //     above) — an engaged-Bland state is never terminal.
            //  2. 4× bland_after *while a Bland episode is stalling*: Bland can
            //     spin slowly on big degenerate LPs (tuff) and its guarantee
            //     erodes under fp noise (scrs8 cycled under it) ⇒ the tie-perturbation
            //     tie-perturbation rungs and the randomized deviation fire,
            //     are sticky-once, and keep working after the episode releases.
            // Sporadic degeneracy (normal in big LPs) never engages stage 1;
            // a genuine cycle accumulates monotonically and does (kb2).
            if bland_left > 0 {
                bland_left -= 1;
                let cleared = !bland_sticky && bland_clear >= self.config.bland_after.max(1);
                if cleared || bland_left == 0 {
                    let reason = if cleared { "stall cleared" } else { "budget" };
                    // Release to normal pricing; reset the stall counter so
                    // re-engagement needs bland_after FRESH stalled pivots
                    // (the growing per-episode budget is the backoff).
                    stall = 0;
                    eprintln!(
                        "rusolve simplex: Bland episode {bland_episodes} released ({reason})"
                    );
                }
            }
            if bland_left == 0 {
                // Re-entry backoff: each successive episode needs TWICE the
                // stalled pivots (capped ×256) before it engages. A failed
                // episode (budget release, no clear) means Bland could not
                // break the pocket it was called for — the healthy state is
                // fast pricing + sticky rungs, so give it growing room.
                let re_enter = self
                    .config
                    .bland_after
                    .max(1)
                    .saturating_mul(1usize << bland_episodes.min(8));
                if stall >= re_enter {
                    bland_episodes += 1;
                    // Phase 1 keeps the legacy STICKY engagement: it has no
                    // randomized escape (randomization there produced false
                    // Infeasible on tuff) and its cure for persistent stalls
                    // is the sticky phase-1 eps rung — releasing Bland there
                    // strands the run with no recovery path, and the Bland stall is a
                    // phase-2 bug. Phase-1 behavior stays bit-identical.
                    let budget = if allow_artificial {
                        0
                    } else {
                        self.config.bland_episode
                    };
                    bland_sticky = budget == 0; // lever: bland_episode=0 = legacy sticky
                    bland_left = if bland_sticky { usize::MAX } else { budget };
                    bland_clear = 0;
                    *bland_triggered = true;
                    eprintln!(
                        "rusolve simplex: Bland episode {bland_episodes} engaged after {stall} stalled pivots (budget {budget})"
                    );
                }
            }
            let use_bland = bland_left > 0;
            // Randomized escape: phase 2 escalates at 4×bland_after stalls;
            // phase 1 needs a far higher bar (16×) — its artificial-sum
            // pricing is noise-sensitive and early random excursions produced
            // false Infeasible (tuff), while scrs8 genuinely cycles there and
            // Bland's guarantee erodes under fp noise.
            // Degenerate-cycle escape (phase 2 only): one-shot random
            // deviations from Bland every 4×bland_after degenerate pivots,
            // with a first-time drift correction (x_b recomputed = B⁻¹b on
            // fresh factors after long −θ·b̄ accumulation).
            // Phase 1 deliberately has NO escape: persistent/one-shot random
            // there produced false Infeasible/Unbounded (tuff) and repeated
            // drift corrections only crawled (scrs8 survives Bland). scrs8's
            // phase-1 noise-driven cycle needs lexicographic pivoting
            // (not implemented).
            // The phase-2 perturbation trigger MUST be evaluated before use_random resets `stall`
            // below: it fires at the SAME rung (4×bland_after while Bland
            // stalls), sticky-once, and does not disturb the existing
            // randomized+drift escalation (which keeps its one-shot cadence;
            // in the ranking chain use_random still wins when it fires).
            if !allow_artificial
                && self.config.p2_tie_perturb
                && use_bland
                && stall >= 4 * self.config.bland_after
                && p2_eps.is_empty()
            {
                p2_eps = p1_tie_eps(total, ptol, 0);
                eprintln!(
                    "rusolve simplex: phase-2 entering tie perturbation engaged after {stall} stalled pivots"
                );
            }
            let use_random = use_bland && !allow_artificial && stall >= 4 * self.config.bland_after;
            if use_random {
                let rung_stall = stall;
                stall = 0; // keep deviations one-shot
                // The deviation REPLACES the episode — drift-corrected
                // factors + fast pricing get the next stretch (the old
                // "one random pivot, Bland resumes" cadence kept the run in
                // Bland's θ=0 crawl forever). If the pocket persists the
                // stall counter re-enters Bland with backoff. Legacy sticky
                // lever keeps the old "Bland resumes between" cadence.
                if !bland_sticky {
                    bland_left = 0;
                }
                if !*randomized {
                    *randomized = true;
                    if pfi.refactorize().is_ok() {
                        refresh_x_b(pfi, st, basis, at_ub, x_b);
                    }
                    eprintln!(
                        "rusolve simplex: randomized deviation + drift correction after {rung_stall} stalled pivots"
                    );
                }
            }
            // Phase-1 ENTERING tie perturbation: the counterpart of the
            // phase-2 randomized escape, firing at the SAME rung
            // (4×bland_after while Bland stalls) where `use_random` refuses
            // to act (it is gated `!allow_artificial`). Evidence
            // (RUSOLVE_P1TRACE, tuff): the phase-1 stall is a θ=0 crawl at
            // ONE vertex (w_true constant over 13k pivots, degenerate tie
            // sets of 30-51 rows containing NO artificials): Bland's
            // lowest-index entering rule walks the tie graph ~13k pivots
            // before escaping. Perturbing the entering COST key reorders
            // that walk. Safety: the improver SET stays defined by the
            // UNPERTURBED test d < -ptol and eps < ptol/2, so perturbation
            // reorders among true improvers only — the phase-1 "no
            // improver ⇒ Infeasible" exit can never be a perturbation
            // artifact (the failure class that killed p2-style
            // randomization in phase 1: tuff false Infeasibles). Sticky
            // once fired, like Bland. Beale 1968; Chvátal 1983 §7.3.
            if allow_artificial
                && use_bland
                && stall >= 4 * self.config.bland_after
                && p1_eps.is_empty()
            {
                p1_eps = p1_tie_eps(total, ptol, p1_eps_seed);
                p1_recheck = stall + 4 * self.config.bland_after.max(1);
                p1_w_probe = artificial_sum(st, basis, x_b);
                eprintln!(
                    "rusolve simplex: phase-1 entering tie perturbation engaged after {stall} stalled pivots"
                );
            }
            let p1_perturb = allow_artificial && !p1_eps.is_empty();
            let p2_perturb = !allow_artificial && !p2_eps.is_empty();
            // HiGHS-style re-perturbation (see SimplexConfig::
            // p1_re_perturb). Fires ONLY while a fired phase-1 eps is stalling:
            // if w (basic-artificial sum) and the basic-artificial count are
            // BIT-identical to the last probe across the stall window, first
            // re-anchor x_b on fresh factors (drift can mask real progress —
            // a re-anchored w that moves is progress, not a freeze), then
            // regenerate eps from a fresh deterministic seed. A cycle of one
            // fixed entering ordering is almost surely not a cycle of the
            // next. Runs where w moves update the probe and stay bit-
            // identical to the one-shot perturbation (the probe is a pure read).
            if p1_perturb && self.config.p1_re_perturb && stall >= p1_recheck {
                p1_recheck = stall + 4 * self.config.bland_after.max(1);
                let (w_now, n_now) = if p1_relaxed_ever {
                    total_infeas(st, basis, x_b)
                } else {
                    artificial_sum(st, basis, x_b)
                };
                let prev_w_ref = if p1_w_ref.is_nan() {
                    p1_w_ref = w_now;
                    w_now
                } else {
                    p1_w_ref
                };
                // The stagnation bar is a RUNNING MINIMUM (never raised):
                // genuine descent must print a new min 0.1% below the best
                // seen. An oscillation that dips and returns (dano3_3:
                // per-window drops that net to zero) therefore accumulates
                // stagnation instead of resetting it.
                if w_now < p1_w_ref {
                    p1_w_ref = w_now;
                }
                let genuine = w_now < prev_w_ref * (1.0 - 1e-3);
                if (w_now, n_now) == p1_w_probe {
                    let mut rolled = false; // set by a burn ROLLBACK below
                    if pfi.refactorize().is_ok() {
                        refresh_x_b(pfi, st, basis, at_ub, x_b);
                        d_ok = false; // re-priced from re-anchored values
                    }
                    let (w_fix, n_fix) = if p1_relaxed_ever {
                        total_infeas(st, basis, x_b)
                    } else {
                        artificial_sum(st, basis, x_b)
                    };
                    if (w_fix, n_fix) == p1_w_probe {
                        p1_eps_seed = p1_eps_seed
                            .wrapping_add(0x9E37_79B9_7F4A_7C15)
                            .wrapping_add(0x517C_1B72);
                        p1_eps = p1_tie_eps(total, ptol, p1_eps_seed);
                        p1_frozen_streak += 1;
                        p1_stag_streak += 1; // frozen ⇒ no new min either
                        // v2: eps reordering alone is provably not enough
                        // (qap10: 999 regens, w frozen — every improver
                        // pivots θ=0 at the closed vertex). Engage the
                        // vertex-mover: random leaver among exact ties.
                        // BAR 50 (not engage-on-first-freeze): the
                        // random leaver overrides the Bland-index tie rule,
                        // and tuff prim (bland_episode=0, agg off, the
                        // phase-2-perturbation acceptance config) escapes a 400-stall pocket
                        // UNDER Bland — arming at the first frozen probe
                        // (streak 1) steered it into a 99k-probe freeze
                        // (bisect: eps-only PASS, with random leaver FAIL).
                        // 50 consecutive frozen probes (5,000 stalled
                        // pivots) is beyond any measured natural escape
                        // (tuff-test escapes at streak 1-2, w moving by
                        // stall 500) and far below brandy's 153-probe
                        // recovery — brandy engages at 50 and keeps
                        // its escape; qap10 (983) and tuff-np (746) cross
                        // the bar inside their frozen stretches as before.
                        if p1_frozen_streak >= 50 {
                            p1_rand_leave = true;
                        }
                        // Relaxed-step engagement is the SINGLE stag bar below
                        // (400 windows without genuine descent), not a
                        // frozen-streak bar. Measured (v5, RUSOLVE_NO_
                        // RELAXED): brandy-np's freeze is long but SELF-
                        // HEALING — bit-move release at stall 30,000, then
                        // Optimal 54,696+92 (obj 1518.5098964881283 =
                        // HiGHS to 2e-16). Any frozen bar that brandy
                        // crosses (200 AND 250 both tested) engages the
                        // walk mid-freeze and destroys that escape even
                        // with burn-rollback: relaxed steps re-cut the
                        // trajectory the escape depends on. The stag bar
                        // waits for PROVEN stagnation (40,000 stalled
                        // pivots without a 0.1% w drop): brandy escapes at
                        // ~300 windows and never engages; qap10's 983-
                        // probe freeze crosses at stall ~41,600 and its
                        // ~20k-pivot walk still finishes inside the cap.
                        // Burn check on the drift-corrected census: the
                        // walk earns its keep as long as it does not make
                        // infeasibility worse (see p1_w_burn_ref decl).
                        if p1_relaxed {
                            if p1_w_burn_ref.is_nan() || w_fix < p1_w_burn_ref {
                                p1_w_burn_ref = w_fix;
                            }
                            if !p1_relaxed_burned && w_fix > 1.2 * p1_w_burn_ref {
                                p1_relaxed = false;
                                p1_relaxed_burned = true;
                                let ref_w = p1_w_burn_ref;
                                p1_w_burn_ref = f64::NAN;
                                if let Some((sb, sx, su, sp)) = p1_snap.take() {
                                    basis.copy_from_slice(&sb);
                                    x_b.copy_from_slice(&sx);
                                    at_ub.copy_from_slice(&su);
                                    if pfi.refactorize().is_ok() {
                                        refresh_x_b(pfi, st, basis, at_ub, x_b);
                                    }
                                    d_ok = false; // re-price from restored basis
                                    p1_w_probe = sp;
                                    rolled = true;
                                }
                                eprintln!(
                                    "rusolve simplex: phase-1 relaxed row choice BURNED at {stall} stalled pivots (w={w_fix:.6e} rose past 1.2× the post-engagement min {ref_w:.6e}) — {}to the exact ratio test until genuine descent",
                                    if rolled { "rolled back" } else { "reverting" }
                                );
                            }
                        }
                        eprintln!(
                            "rusolve simplex: phase-1 freeze escape after {stall} stalled pivots (w={w_fix:.6e} frozen through drift correction, streak {p1_frozen_streak}): eps re-randomized (seed {p1_eps_seed:#x}){}{}",
                            if p1_rand_leave {
                                " + randomized leaving ties"
                            } else {
                                ""
                            },
                            if p1_relaxed {
                                " + relaxed row choice"
                            } else {
                                ""
                            }
                        );
                    }
                    if !rolled {
                        p1_w_probe = (w_fix, n_fix);
                    }
                } else {
                    // w moved since the last probe (bit-move at minimum).
                    // Same burn check on the raw census (progress ⇒ no
                    // drift-correction needed).
                    let mut rolled = false;
                    if p1_relaxed {
                        if w_now < p1_w_burn_ref || p1_w_burn_ref.is_nan() {
                            p1_w_burn_ref = w_now;
                        } else if !p1_relaxed_burned && w_now > 1.2 * p1_w_burn_ref {
                            p1_relaxed = false;
                            p1_relaxed_burned = true;
                            let ref_w = p1_w_burn_ref;
                            p1_w_burn_ref = f64::NAN;
                            if let Some((sb, sx, su, sp)) = p1_snap.take() {
                                basis.copy_from_slice(&sb);
                                x_b.copy_from_slice(&sx);
                                at_ub.copy_from_slice(&su);
                                if pfi.refactorize().is_ok() {
                                    refresh_x_b(pfi, st, basis, at_ub, x_b);
                                }
                                d_ok = false; // re-price from restored basis
                                p1_w_probe = sp;
                                rolled = true;
                            }
                            eprintln!(
                                "rusolve simplex: phase-1 relaxed row choice BURNED at {stall} stalled pivots (w={w_now:.6e} rose past 1.2× the post-engagement min {ref_w:.6e}) — {}to the exact ratio test until genuine descent",
                                if rolled { "rolled back" } else { "reverting" }
                            );
                        }
                    }
                    // Two progress bars:
                    // bit-movement (releases the v2 random-tie walk —
                    // brandy's measured recovery semantics, unchanged) and
                    // RELATIVE descent against the running min (≥ 0.1% —
                    // `genuine`, computed above: the stagnation bar for
                    // the relaxed step. dano3_3 oscillates at w≈8.2e4 with bit-moves
                    // but no net descent for the whole 100k run; the relaxed step
                    // must engage through stagnation, and once engaged
                    // survive oscillation — released only by genuine
                    // descent).
                    if genuine {
                        p1_stag_streak = 0;
                        if p1_relaxed_burned {
                            // Genuine descent is the only re-arm: v2 (or the
                            // problem itself) is making real progress again,
                            // so a LATER plateau may legitimately re-engage
                            // the relaxed step. dano3_3-shaped oscillations (bit-moves
                            // without descent) keep the lock.
                            p1_relaxed_burned = false;
                            p1_w_burn_ref = f64::NAN;
                            eprintln!(
                                "rusolve simplex: phase-1 relaxed row choice re-armed at {stall} stalled pivots — genuine descent"
                            );
                        }
                        if p1_rand_leave {
                            // w demonstrably descending again — hand the
                            // ratio test back to the deterministic
                            // (Bland-index) tie rule.
                            p1_rand_leave = false;
                            p1_relaxed = false;
                            p1_frozen_streak = 0;
                            eprintln!(
                                "rusolve simplex: phase-1 freeze escape released at {stall} stalled pivots — w descending again"
                            );
                        }
                    } else if !rolled {
                        p1_stag_streak += 1;
                        // Bit-move release (v2 semantics) applies only
                        // while the EXACT walk runs: with the relaxed row
                        // choice engaged, a w bit-move without descent is
                        // the walk's own truncated step, not recovery —
                        // releasing here is what capped brandy-np (the
                        // leaver died one probe after engagement, at stall
                        // 25,900, and Bland crawled the rest of 74k
                        // pivots). Hand back only on genuine descent (the
                        // branch above) or once relaxed is off.
                        if p1_rand_leave && !p1_relaxed {
                            p1_rand_leave = false;
                            eprintln!(
                                "rusolve simplex: phase-1 freeze escape released at {stall} stalled pivots — w moving again"
                            );
                        }
                    }
                    if !rolled {
                        p1_w_probe = (w_now, n_now);
                    }
                }
                // Relaxed-step SOLE engagement latch — hoisted OUT of the
                // frozen/w-moved split: a pure bit-freeze (qap10) never
                // enters the w-moved branch, so a latch inside it can
                // never fire (v6 regression: qap10 Limit with streak 983
                // and no engagement). Bar: 600 CONSECUTIVE BIT-FROZEN
                // probes (60,000 stalled pivots with (w,n) identical
                // through drift-corrected re-anchoring). Freezing alone
                // is not evidence and stag without freezing is not enough
                // either — the bar must be crossed by BOTH:
                //   brandy-np (v6, Optimal 54,696+92): its freeze is
                //   TERMINAL-but-escapable — w=9.58e2 bit-frozen for 535
                //   consecutive probes, broken by the phase-1 exit
                //   itself; any bar ≤535 engages mid-freeze, re-cuts the
                //   trajectory and caps the run (measured at bars 200,
                //   250 and at a 400-window stag-only latch).
                //   qap10: 983 frozen probes — crosses 600 at stall
                //   ~61,600, and the measured ~20.8k-pivot walk finishes
                //   phase 1 near ~82,400, inside the 100k cap (bar 800
                //   starves it: engaged 81,600 ⇒ Limit).
                //   dano3_3: oscillates (max frozen streak 109), never
                //   engages — its honest verdict stays the Limit table.
                // The burned lock (set by every burn, cleared only by
                // genuine descent) keeps post-burn windows (and rollback
                // windows, whose burn set the lock) from re-engaging.
                // The snapshot probe tuple is p1_w_probe: both branches
                // have just written it with their own census read.
                if p1_frozen_streak >= 600
                    && p1_stag_streak >= 400
                    && !p1_relaxed_burned
                    && std::env::var_os("RUSOLVE_NO_RELAXED").is_none()
                {
                    if !p1_relaxed {
                        eprintln!(
                            "rusolve simplex: phase-1 escape engaged after {stall} stalled pivots (frozen streak {p1_frozen_streak}, stagnant windows {p1_stag_streak}): relaxed row choice"
                        );
                        p1_snap = Some((basis.to_vec(), x_b.to_vec(), at_ub.to_vec(), p1_w_probe));
                    }
                    p1_relaxed = true;
                    p1_relaxed_ever = true;
                    // v3: NO burn-ref reset — the walk inherits the
                    // descent min established before the plateau, so
                    // a relaxed step that pushes w back above 1.2× the
                    // best w seen burns immediately.
                }
            }
            let p1_probe = std::env::var("RUSOLVE_P1PROBE").is_ok();
            if p1_probe && allow_artificial && stall >= 20_000 && stall.is_multiple_of(20_000) {
                // TEMPORARY diagnostic: at a frozen vertex, do ANY improver
                // columns have θ>0 (strict min-ratio), and how far could a
                // step go if structural bound-blockers were ignored (θ over
                // artificial-row caps only, the relaxed-step ceiling)?
                let mut n_imp = 0u32;
                let mut n_pos = 0u32;
                let mut top: Vec<(f64, usize)> = Vec::new(); // strict θ
                let mut top_relax: Vec<(f64, usize)> = Vec::new();
                for wi in 0..cand_bits.len() {
                    let mut word = cand_bits[wi];
                    while word != 0 {
                        let v = wi * 64 + word.trailing_zeros() as usize;
                        word &= word - 1;
                        n_imp += 1;
                        let sgn = if at_ub[v] { -1.0 } else { 1.0 };
                        st.gather_col(v, sgn, &mut aq);
                        pfi.solve_b(&aq, &mut bbar);
                        let mut theta = st.ub[v];
                        let mut theta_art = st.ub[v];
                        for i in 0..m {
                            let bi_bar = bbar[i];
                            let (cap, _) = if bi_bar > tol {
                                ((x_b[i] / bi_bar).max(0.0), false)
                            } else if bi_bar < -tol && st.ub[basis[i]].is_finite() {
                                (((st.ub[basis[i]] - x_b[i]) / -bi_bar).max(0.0), true)
                            } else {
                                continue;
                            };
                            if cap < theta {
                                theta = cap;
                            }
                            if st.kind[basis[i]] == Kind::Artificial && cap < theta_art {
                                theta_art = cap;
                            }
                        }
                        if theta > 1e-12 {
                            n_pos += 1;
                            top.push((theta, v));
                        }
                        if theta_art > 1e-9 {
                            top_relax.push((theta_art, v));
                        }
                    }
                }
                top.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
                top_relax.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
                eprintln!(
                    "p1probe it={iters} stall={stall} improvers={n_imp} theta>1e-12: {n_pos} \
                     top_theta={:?} top_theta_artonly={:?}",
                    &top[..top.len().min(5)],
                    &top_relax[..top_relax.len().min(5)]
                );
            }
            // Relaxed phase-1 repair pass (see repair_head decl): while any
            // structural basic row is out of range, drive the worst one out
            // with a θ=0 pivot on any live nonbasic column that has support
            // in its tableau row (rho = Bᵀe_pos through the row view — the
            // pivotal-row machinery, reused).
            repair_head = None;
            let mut repair_q: Option<usize> = None;
            if allow_artificial && !p1_repair_off {
                let mut vpos: Option<usize> = None;
                let mut vw = 1e-7f64;
                for (i, &v) in basis.iter().enumerate() {
                    if st.kind[v] == Kind::Artificial {
                        continue;
                    }
                    let d = (-x_b[i]).max(x_b[i] - st.ub[v]);
                    if d > vw {
                        vw = d;
                        vpos = Some(i);
                    }
                }
                if let Some(pos) = vpos {
                    e_pos.fill(0.0);
                    e_pos[pos] = 1.0;
                    pfi.solve_bt(&e_pos, &mut pi_pos);
                    touched.clear();
                    for (i, &r) in pi_pos.iter().enumerate() {
                        if r == 0.0 {
                            continue;
                        }
                        for k in rptr[i]..rptr[i + 1] {
                            let c = rcol[k];
                            if in_basis[c] || dead[c] {
                                continue;
                            }
                            if !row_mark[c] {
                                row_mark[c] = true;
                                touched.push(c);
                            }
                            row_alpha[c] += r * rval[k];
                        }
                    }
                    for &c in touched.iter() {
                        if row_alpha[c].abs() > 1e-9 && st.ub[c] > 0.0 {
                            repair_q = Some(c);
                            repair_head = Some(pos);
                            break;
                        }
                    }
                    for &c in touched.iter() {
                        row_alpha[c] = 0.0;
                        row_mark[c] = false;
                    }
                    if repair_head.is_none() {
                        repair_fails += 1;
                        if repair_fails > 8 {
                            p1_repair_off = true;
                        }
                    }
                }
            }
            let mut n_cand = 0u32;
            let mut best: Option<(usize, bool, f64)> = None; // (var, from_ub, reduced)
            // Candidate bitset (maintained with dvec): bit v set <=> v is
            // nonbasic, priceable and improving. Scanned in ascending v, so
            // the visit order equals the former full 0..total scan.
            for wi in 0..cand_bits.len() {
                let mut word = cand_bits[wi];
                while word != 0 {
                    let v = wi * 64 + word.trailing_zeros() as usize;
                    word &= word - 1;
                    let d = dvec[v];
                    n_cand += 1;
                    let better = match best {
                        None => true,
                        Some((bv, _, br)) => {
                            if use_random {
                                rng.next_u64().is_multiple_of(u64::from(n_cand))
                            } else if p1_perturb {
                                // Phase 1: argmin of (d + eps[v]) over TRUE improvers
                                // (membership already tested on the UNPERTURBED
                                // d above; eps < ptol/2 ⇒ nothing is reclassified
                                // in either direction).
                                d + p1_eps[v] < br + p1_eps[bv]
                            } else if p2_perturb {
                                // Phase 2: argmin of the EXISTING phase-2 key with the
                                // perturbed cost, over TRUE improvers (membership
                                // already tested on the UNPERTURBED d above; eps <
                                // ptol/2 keeps every improver's key the same sign,
                                // so nothing is reclassified either direction).
                                // SH path: (d+eps)²/w on the carried weights (the
                                // d/√w ranking squared); Dantzig path: |d+eps|.
                                if sh_active {
                                    let (dk, bk) = (d + p2_eps[v], br + p2_eps[bv]);
                                    let (dm, bm) =
                                        (dk * dk / w[v].max(1e-8), bk * bk / w[bv].max(1e-8));
                                    let eps = 1e-12 * (dm + bm + 1.0);
                                    dm > bm + eps || ((dm - bm).abs() <= eps && v < bv)
                                } else {
                                    let (dm, bm) = ((d + p2_eps[v]).abs(), (br + p2_eps[bv]).abs());
                                    dm > bm + 1e-12 || ((dm - bm).abs() <= 1e-12 && v < bv)
                                }
                            } else if use_bland {
                                v < bv
                            } else if sh_active && !allow_artificial {
                                // Steepest-edge rate (PHASE 2 ONLY): objective
                                // drop per unit of edge traveled, d²/w with the
                                // carried exact weights. Phase 1 keeps Dantzig —
                                // SH pricing there amplified its degeneracy
                                // (tuff p1 3,863→58,762); phase 2 is where SE
                                // pays (tuff p2 1,380→116).
                                let (dm, bm) = (d * d / w[v].max(1e-8), br * br / w[bv].max(1e-8));
                                let eps = 1e-12 * (dm + bm + 1.0);
                                dm > bm + eps || ((dm - bm).abs() <= eps && v < bv)
                            } else {
                                let (dm, bm) = (d.abs(), br.abs());
                                dm > bm + 1e-12 || ((dm - bm).abs() <= 1e-12 && v < bv)
                            }
                        }
                    };
                    if better {
                        best = Some((v, at_ub[v], d));
                    }
                }
            }
            if let Some(rq) = repair_q {
                // Repair pivot replaces the priced entering column.
                best = Some((rq, at_ub[rq], dvec[rq]));
            }
            let Some((q, from_ub, _d)) = best else {
                if !pi_fresh {
                    // Maintained d says optimal: confirm on a fresh BTRAN.
                    d_ok = false;
                    iters -= 1;
                    continue;
                }
                if allow_artificial && !p1_exit_verified {
                    // Verify a phase-1 pricing exit on FRESH factors
                    // before trusting it — HiGHS forces a rebuild at
                    // possibly-phase-1-optimal (kRebuildReasonPossibly
                    // Phase1Feasible, HEkkPrimal.cpp:912–915; MIT). The
                    // maintained x_b accumulates per-pivot drift on dirty
                    // trajectories (e226 prim+pre exited with structural
                    // violations ~6e19 while the SAME factors' pricing said
                    // "no improver" — consistent garbage). Re-anchor x_b from
                    // a refactorization and re-price: improvers that appear
                    // on the honest π continue the phase; only a second
                    // no-improver on fresh factors exits (the p1 exit gate
                    // then still audits x for real violations).
                    p1_exit_verified = true;
                    if pfi.refactorize().is_ok() {
                        refresh_x_b(pfi, st, basis, at_ub, x_b);
                    }
                    d_ok = false; // re-price from the re-anchored π
                    iters -= 1; // verification is not a simplex iteration
                    continue;
                }
                if p1trace {
                    tr.emit(iters - 1, use_bland, p1_perturb);
                }
                exit.pi.clone_from(&pi); // Pricing π at the exit basis
                return Ok((iters - 1, PhaseOutcome::Optimal));
            };

            // ---- FTRAN of entering column (sign = direction) ----
            // gather_col ≡ the old dense clone/negate-map:
            // fill(0.0) + scatter writes exactly the vector the clone held,
            // and the ±1 scale is exact. O(m + nnz(a_q)) ≤ 1 % of the solve.
            let sgn = if from_ub { -1.0 } else { 1.0 };
            st.gather_col(q, sgn, &mut aq);
            // FTRAN b̄: refinement gated. BTRAN π (pricing) is always
            // refined; b̄ steers the ratio test, where mild drift is survivable.
            if self.config.refine_ftran {
                pfi.solve_b(&aq, &mut bbar);
            } else {
                pfi.solve_b_raw(&aq, &mut bbar);
            }

            nzb.clear();
            for (i, &b) in bbar.iter().enumerate() {
                if b != 0.0 {
                    nzb.push(i);
                }
            }
            // Maintained-d consistency: d_q = c_q - c_B^T B^-1 a_q is free from
            // the FTRAN just done; a mismatch means drift -> full recompute.
            let dq_true = if pi_fresh || repair_head.is_some() {
                dvec[q]
            } else {
                let cbb: f64 = nzb.iter().map(|&i| cb[i] * bbar[i]).sum();
                let dq = cost[q] - sgn * cbb;
                let still = if from_ub { dq > ptol } else { dq < -ptol };
                if !still || (dq - dvec[q]).abs() > 1e-7 * (1.0 + dq.abs()) {
                    d_ok = false;
                    iters -= 1;
                    continue;
                }
                dq
            };

            // ---- Bounded ratio test, Harris tie-break in phase 2 ----
            // Phase 2 is TWO-PASS: (1) true min-ratio THETA — the step stays
            // exactly min-feasible (no basic pushed below the round-off
            // clamp); (2) among rows TIED at that ratio within a near
            // window, leave the row with the LARGEST |b̄| (Harris, ORSA J.
            // Computing 4(1), 1992; Maros ch. 6): its ratio carries the
            // smallest relative error under b̄-drift, and the pivot
            // maximally perturbs the basis — the documented escape from the
            // θ≈0 stall basin where per-pivot objective progress falls
            // below fp resolution (measured: the scsd6 2-cycle and the
            // tuff 76k-iter stall). Under Bland, near-equal
            // |b̄| ties keep the lowest leaving VARIABLE index (kb2
            // anti-cycling evidence).
            // Phase 1 → the OLD single pass VERBATIM (first-come-kept scan
            // order is load-bearing here: Harris applied to phase 1
            // sent tuff into an artificial p1 cycle; and a phase-2 theta
            // shadowing bug here, once fixed, is what the 15,140→2,900
            // scsd1 / 56,878→27,923 scsd6 iteration wins are measured on).
            let mut theta = st.ub[q]; // entering may saturate at its own ub
            let mut leave: Option<usize> = None; // row position
            let mut leave_to_ub = false;
            if allow_artificial {
                // Tie window: exact-arithmetic ties only. Under Bland ties
                // break to the lowest leaving variable index (anti-cycling
                // guarantee); position order is NOT enough (kb2 evidence).
                let tie = 1e-12;
                // Reservoir counter over same-kind tie candidates
                // (only consumed while p1_rand_leave).
                let mut tie_seen = 0u32;
                for &i in nzb.iter() {
                    let bi_bar = bbar[i];
                    let (cap, to_ub) = if bi_bar > tol {
                        ((x_b[i] / bi_bar).max(0.0), false)
                    } else if bi_bar < -tol && st.ub[basis[i]].is_finite() {
                        (((st.ub[basis[i]] - x_b[i]) / -bi_bar).max(0.0), true)
                    } else {
                        continue;
                    };
                    let better = cap < theta - tie
                        || ((cap - theta).abs() <= tie
                            && match leave {
                                None => true,
                                Some(lpos) => {
                                    // Drive-out preference (same sources
                                    // as the deletion rule): among EXACT
                                    // min-ratio ties, an artificial leaver
                                    // beats a structural one — ejecting it
                                    // strictly decreases the basic-artificial
                                    // count, the secondary progress measure.
                                    // Same-kind ties keep the old semantics
                                    // VERBATIM (first-come-kept; under Bland
                                    // lowest leaving variable index). Step
                                    // length (theta) untouched: unlike the
                                    // Harris window, no non-tied row can
                                    // enter the tie set here.
                                    let art_i = st.kind[basis[i]] == Kind::Artificial;
                                    let art_l = st.kind[basis[lpos]] == Kind::Artificial;
                                    if art_i != art_l {
                                        art_i
                                    } else if p1_rand_leave {
                                        // Reservoir-style random
                                        // leaver among exact ties. θ is
                                        // unchanged ⇒ x is unchanged and the
                                        // improver set is untouched (the
                                        // false-Infeasible class that killed
                                        // entering-randomization in phase 1
                                        // cannot occur); the BASIS changes,
                                        // so each θ=0 pivot lands on a
                                        // DIFFERENT adjacent vertex. The
                                        // drive-out preference above stays
                                        // dominant; draws come from the
                                        // run's deterministic RNG ⇒ the
                                        // escape is reproducible.
                                        tie_seen += 1;
                                        rng.next_u64().is_multiple_of(tie_seen as u64 + 1)
                                    } else if self.config.p1_tie_pref_always
                                        && bbar[i].abs() > bbar[lpos].abs()
                                    {
                                        // Rung arm: UNCONDITIONAL
                                        // largest-|b̄| among exact ties —
                                        // the measured e226 prim+pre Optimal
                                        // 2,363+155 semantics.
                                        // Only the dirty-
                                        // phase-1-exit re-solve rung in
                                        // solve_dispatch sets this flag, so
                                        // clean-exit paths (tuff) never see
                                        // it — the measured tuff prim-np
                                        // cycling was the DEFAULT-path arm.
                                        // θ is the tied cap ⇒ x and the
                                        // improver set are UNCHANGED.
                                        true
                                    } else if stab_hist.count_ones() >= 3
                                        && bbar[lpos].abs() < 1e-7
                                        && bbar[i].abs() > bbar[lpos].abs()
                                    {
                                        // Default-path stability
                                        // preference among EXACT ties, armed
                                        // ONLY by a pivot-instability burst:
                                        // ≥ 3 sub-1e-7 pivots within the last
                                        // 10 (the measured discriminator —
                                        // e226's detonation window runs 9/50
                                        // while tuff's whole 77k-iteration
                                        // Optimal path never exceeds 2/50, so
                                        // this cannot fire there). Swaps an
                                        // unstable tie winner (|b| < 1e-7,
                                        // HiGHS dPivotTol for update_count
                                        // ≥ 20) for a larger-|b| tie partner:
                                        // θ is the tied cap ⇒ x and the
                                        // improver set are UNCHANGED (the
                                        // falsified Harris-window variant
                                        // stepped past NEAR-ties and cycled
                                        // tuff; an unconditional or
                                        // winner-only-gated |b| preference
                                        // measurably cycled tuff prim-np as
                                        // a DEFAULT arm too). e226 prim+pre
                                        // detonated through a θ=0 degenerate
                                        // stretch landing |b|~1e-9 pivots
                                        // (basis conditioning decayed, x_b →
                                        // 6e19, phase-1 exit on a garbage
                                        // point): HiGHS's phase-1 row choice
                                        // filters leavers by |α| ≥ 0.1·max
                                        // (phase1ChooseRow, MIT) — this is
                                        // the exact-tie equivalent, burst-
                                        // armed.
                                        true
                                    } else {
                                        use_bland && basis[i] < basis[lpos]
                                    }
                                }
                            });
                    if better {
                        theta = cap;
                        leave = Some(i);
                        leave_to_ub = to_ub;
                    }
                }
                if p1trace && theta < 1e-11 {
                    // TEMP diagnostic: tie-set size at the final theta (same
                    // 1e-12 window as the scan). Extra O(m) pass, degenerate
                    // pivots only, env-gated.
                    let mut ties = 0usize;
                    let mut any_art = false;
                    for &i in nzb.iter() {
                        let bi_bar = bbar[i];
                        let cap = if bi_bar > tol {
                            (x_b[i] / bi_bar).max(0.0)
                        } else if bi_bar < -tol && st.ub[basis[i]].is_finite() {
                            ((st.ub[basis[i]] - x_b[i]) / -bi_bar).max(0.0)
                        } else {
                            continue;
                        };
                        if (cap - theta).abs() <= 1e-12 {
                            ties += 1;
                            if st.kind[basis[i]] == Kind::Artificial {
                                any_art = true;
                            }
                        }
                    }
                    tr.count_degenerate(ties, any_art);
                }
                // ---- Relaxed phase 1: HiGHS-style TOLERANCE-TRUNCATED
                // phase-1 row choice ----
                // Port of HEkkPrimal::phase1ChooseRow (HiGHS, MIT; algorithm
                // + constants, not code). Walk the entering ray over the
                // TRUNCATED infeasibility f_tol(θ) = Σ max(0, violation −
                // tolf): rows inside tolerance are FREE, so a degenerate
                // vertex (many rows exactly at bounds — the measured qap10
                // freeze: strict |violation| slope ≥ 0 on EVERY edge) still
                // exposes descent. Descent rate starts at the sign-correct
                // out-of-range contribution and loses |b̄| at every
                // RELAXED breakpoint (tolf-shifted, as HiGHS's
                // dRelaxTheta/dFeasTheta); dMaxTheta is the last breakpoint
                // with positive rate after subtracting. Leaving row: among
                // TIGHT caps ≤ dMaxTheta (HiGHS scans descending), take the
                // largest cap with |b̄| ≥ 0.1·max. Rows crossed on the way
                // keep violations ≤ tolf + (dMaxTheta − their relax)·|b̄|
                // (HiGHS likewise: created infeasibilities just re-enter the
                // phase-1 objective) — owned by total_infeas (probe), the
                // θ=0 repair pass, and the phase-1 exit gate (> 1e-7 per
                // row ⇒ honest Limit, never a wrong verdict). When there is
                // no truncated descent (grad0 ≤ 0 ⇒ dMaxTheta = 0) the v2
                // exact-tie path stands ⇒ brandy behavior untouched.
                // Basis: HiGHS simplex/HEkkPrimal.cpp:phase1ChooseRow.
                if p1_relaxed {
                    p1_events.clear();
                    let tolf = 5e-8; // truncation tol (HiGHS 1e-7, halved
                    // for margin against the 1e-7 exit gate / fp clamps)
                    let ubq = st.ub[q];
                    let mut grad = 0.0f64; // descent rate of f_tol at 0+
                    for &i in nzb.iter() {
                        let alpha = bbar[i];
                        let x = x_b[i];
                        // truncated range [lo, hi]: artificials live at 0
                        // (violation = |value|), structural rows in [0, ub].
                        let (lo, hi) = if st.kind[basis[i]] == Kind::Artificial {
                            (0.0, 0.0)
                        } else {
                            (0.0, st.ub[basis[i]])
                        };
                        if alpha > tol {
                            // value falls
                            if x > hi + tolf {
                                grad += alpha; // violation shrinking
                                p1_events.push(((x - hi - tolf) / alpha, alpha));
                            }
                            if x > lo - tolf {
                                p1_events.push(((x - lo + tolf) / alpha, alpha));
                            }
                            if x < lo - tolf {
                                grad -= alpha; // worsening, never re-enters
                            }
                        } else if alpha < -tol {
                            // value rises
                            if x < lo - tolf {
                                grad += -alpha; // violation shrinking
                                p1_events.push(((x - lo + tolf) / alpha, -alpha));
                            }
                            if hi.is_finite() && x < hi + tolf {
                                p1_events.push(((x - hi - tolf) / alpha, -alpha));
                            }
                            if x > hi + tolf {
                                grad -= -alpha; // worsening, never re-enters
                            }
                        }
                    }
                    let mut dmax = 0.0f64;
                    if p1trace {
                        rx_tried += 1;
                    }
                    if grad > 0.0 {
                        p1_events.sort_by(|a, b| {
                            a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
                        });
                        let mut g = grad;
                        for &(te, dec) in p1_events.iter() {
                            if te > ubq {
                                break;
                            }
                            g -= dec;
                            if g <= 0.0 {
                                break;
                            }
                            dmax = te;
                        }
                    } else if p1trace {
                        rx_s0pos += 1;
                    }
                    if dmax <= 0.0 && p1trace {
                        rx_dmax0 += 1;
                    }
                    if p1trace && iters.is_multiple_of(5000) {
                        eprintln!(
                            "p1rx it={iters} tried={rx_tried} won={rx_won} grad<=0:{rx_s0pos} dmax=0:{rx_dmax0}"
                        );
                    }
                    if dmax > theta {
                        // Candidates: TIGHT cap ≤ dmax (violated rows have
                        // cap 0 ⇒ always candidates, they want out anyway).
                        // Artificial tight cap = distance to 0.
                        let cap_of = |i: usize| -> Option<(f64, bool)> {
                            let alpha = bbar[i];
                            let x = x_b[i];
                            if st.kind[basis[i]] == Kind::Artificial {
                                if alpha * x > 0.0 {
                                    Some((x / alpha, false))
                                } else {
                                    None
                                }
                            } else if alpha > tol {
                                Some(((x / alpha).max(0.0), false))
                            } else if alpha < -tol && st.ub[basis[i]].is_finite() {
                                Some((((st.ub[basis[i]] - x) / -alpha).max(0.0), true))
                            } else {
                                None
                            }
                        };
                        let mut dmax_alpha = 0.0f64;
                        for &i in nzb.iter() {
                            let Some((cap, _)) = cap_of(i) else {
                                continue;
                            };
                            if cap <= dmax && bbar[i].abs() > dmax_alpha {
                                dmax_alpha = bbar[i].abs();
                            }
                        }
                        let mut rtheta = 0.0f64;
                        let mut rleave: Option<usize> = None;
                        let mut rtoub = false;
                        for &i in nzb.iter() {
                            let Some((cap, to_ub)) = cap_of(i) else {
                                continue;
                            };
                            if cap > dmax || bbar[i].abs() < 0.1 * dmax_alpha {
                                continue;
                            }
                            if rleave.is_none() || cap > rtheta {
                                rtheta = cap;
                                rleave = Some(i);
                                rtoub = to_ub;
                            }
                        }
                        if let Some(pos) = rleave {
                            theta = rtheta;
                            leave = Some(pos);
                            leave_to_ub = rtoub;
                            if p1trace {
                                rx_won += 1;
                            }
                        }
                    }
                }
            } else {
                // Pass 1: true min-ratio theta (the step itself stays
                // exactly min-feasible regardless of which row leaves).
                for &i in nzb.iter() {
                    let bi_bar = bbar[i];
                    let cap = if bi_bar > tol {
                        // .max(0.0) is a no-op in normal arithmetic (x_b >= 0
                        // is the ratio test's own invariant); it only absorbs
                        // fp noise that would otherwise make theta negative.
                        (x_b[i] / bi_bar).max(0.0)
                    } else if bi_bar < -tol && st.ub[basis[i]].is_finite() {
                        ((st.ub[basis[i]] - x_b[i]) / -bi_bar).max(0.0)
                    } else {
                        continue;
                    };
                    if cap < theta {
                        theta = cap;
                    }
                }
                // Pass 2 (Harris): pick the leaving row among the ties.
                // Window: degenerate families (θ up to ~1e-9 at x-scale ~10)
                // group. BUT a row is only a legal tie candidate if it
                // LANDS within the feasibility tolerance of its bound:
                // selecting cap = θ+δ records the leaver at its bound while
                // it truly lands δ·|b̄| away — an ungated window (1e-9) made
                // that per-pivot lie ~1e-9·|b̄| and accumulated to a 0.625
                // bound violation on scsd6-sparse (false Optimal 50.0 vs
                // 50.5, exact residual — measured). Gate the lie
                // at the same 1e-7 the round-off clamp already accepts.
                let win = 1e-9 * (1.0 + theta);
                let mut best_bar = 0.0f64;
                for &i in nzb.iter() {
                    let bi_bar = bbar[i];
                    let (cap, to_ub) = if bi_bar > tol {
                        ((x_b[i] / bi_bar).max(0.0), false)
                    } else if bi_bar < -tol && st.ub[basis[i]].is_finite() {
                        (((st.ub[basis[i]] - x_b[i]) / -bi_bar).max(0.0), true)
                    } else {
                        continue;
                    };
                    if cap > theta + win {
                        continue;
                    }
                    let abar = bi_bar.abs();
                    if (cap - theta) * abar > 1e-7 {
                        continue;
                    }
                    let better = match leave {
                        None => true,
                        Some(lpos) => {
                            abar > best_bar * (1.0 + 1e-12)
                                || (abar >= best_bar * (1.0 - 1e-12)
                                    && use_bland
                                    && basis[i] < basis[lpos])
                        }
                    };
                    if better {
                        best_bar = abar;
                        leave = Some(i);
                        leave_to_ub = to_ub;
                    }
                }
            }

            if let Some(rpos) = repair_head {
                // Repair pivot: θ=0, forced leaver (the violated row). The
                // support came from rho·a_q; the FTRAN above must see it too
                // (same factors) — if it does not, the fp paths disagree:
                // count and eventually give up rather than spin.
                if bbar[rpos].abs() <= tol {
                    repair_fails += 1;
                    if repair_fails > 8 {
                        p1_repair_off = true;
                    }
                    d_ok = false;
                    iters -= 1;
                    continue;
                }
                theta = 0.0;
                leave = Some(rpos);
                // Record the leaver at the bound it came from: over-ub ⇒ ub,
                // under-lb ⇒ lb (the lie equals the violation, no more).
                leave_to_ub = x_b[rpos] - st.ub[basis[rpos]] > 0.0;
            }

            if p1trace && allow_artificial {
                // Diagnostic: sub-1e-7 pivot density + x_b detonation.
                let mx = x_b.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                if leave.is_some_and(|r| bbar[r].abs() < 1e-7) {
                    rx_unstable += 1;
                }
                rx_pivots += 1;
                if rx_pivots.is_multiple_of(50) {
                    eprintln!("p1unstab it={iters} last50={rx_unstable} maxxb={mx:.3e}");
                    rx_unstable = 0;
                }
                if theta > 1e4 || mx > 1e8 || leave.is_some_and(|r| bbar[r].abs() < 1e-7) {
                    eprintln!(
                        "p1boom it={iters} q={q} theta={theta:.3e} piv={:.3e} maxxb={mx:.3e}",
                        leave.map(|r| bbar[r]).unwrap_or(0.0)
                    );
                }
            }
            if allow_artificial {
                // Arm/disarm the tie stability preference by the
                // recent sub-1e-7 pivot burst (see stab_hist decl). Bound
                // flips (leave = None) shift in a stable bit.
                stab_hist =
                    (stab_hist << 1) | u16::from(leave.is_some_and(|r| bbar[r].abs() < 1e-7));
            }
            match leave {
                None if !theta.is_finite() => {
                    // Numerical guard: an "unbounded ray" computed on drifted
                    // etas can be fake — phase 1's artificial sum is bounded
                    // below by 0, so a phase-1 ray is ALWAYS an artifact, and
                    // phase-2 rays must survive fresh factors to be trusted.
                    // Verify once (refactorize + drift-correct x_b), then
                    // re-price/re-FTRAN; declare only if the ray persists.
                    // (scsd1 false-Infeasible was this verdict misread
                    // downstream: phase 1 Unbounded fell through to the
                    // p1_obj>tol branch.)
                    if !ray_verified {
                        ray_verified = true;
                        let rr = pfi.refactorize();
                        if std::env::var("RUSOLVE_DEBUG").is_ok() {
                            let max_bbar = bbar.iter().fold(0.0f64, |a, &v| a.max(v));
                            eprintln!(
                                "rusolve simplex: RAY-VERIFY q={q} from_ub={from_ub} \
                                 refactorize={rr:?} max_bbar={max_bbar}"
                            );
                        }
                        match rr {
                            Ok(()) => {
                                refresh_x_b(pfi, st, basis, at_ub, x_b);
                            }
                            // A numerically singular basis poisons ALL factor
                            // information (blind pricing, fake rays): repair it
                            // by swapping the dependent column for a non-basic
                            // artificial. Sound in phase 1, where the artificial
                            // basis is ground truth; phase 2 keeps the plain
                            // verdict (finite ubs ⇒ no artificial-basis
                            // fallback semantics). scsd1 false-Infeasible.
                            Err(FactorError::Singular { column, .. }) if allow_artificial => {
                                if repair_basis(
                                    st,
                                    pfi,
                                    basis,
                                    x_b,
                                    at_ub,
                                    column,
                                    self.config.pfi,
                                    self.deadline,
                                ) {
                                    weights_stale = true; // basis swapped: re-init
                                    d_ok = false;
                                    iters -= 1; // repair pass is not a pivot
                                    continue;
                                }
                            }
                            // Phase 2 has no artificial-basis fallback: a basis
                            // whose REFACTORIZATION is singular is numerically
                            // rank-deficient — every solve (π, b̄, d) is a
                            // fabrication and the "ray" is fake (scsd1-SH:
                            // Err(Singular{col:49}) then RAY-DECLARE d=-13.4
                            // with πᵀa_q=15.6 vs c_Bᵀb̄=-2.2). SH pricing is
                            // exploratory ⇒ bail to the Dantzig fallback; the
                            // baseline path keeps its historical behavior.
                            Err(_) => {
                                if sh_active && !allow_artificial {
                                    return Ok((iters, PhaseOutcome::NumericalDistress));
                                }
                            }
                        }
                        iters -= 1; // a verification pass is not a pivot
                        continue;
                    }
                    // Consistency identity: πᵀa_q must equal c_Bᵀ(B⁻¹a_q) for
                    // ANY consistent factor state. Mismatch ⇒ factors/basis
                    // corrupt ⇒ the ray is fabricated, not real. Under SH
                    // (exploratory pricing) bail to the Dantzig fallback;
                    // baseline keeps its historical declare-anyway behavior.
                    // sparse_dot ≡ the old dense zip-sum (same products,
                    // ascending-row order, zeros are no-ops).
                    if !pi_fresh {
                        pfi.solve_bt(&cb, &mut pi);
                    }
                    let dot: f64 = sparse_dot(&pi, st, q);
                    let cbtb: f64 = cb.iter().zip(&bbar).map(|(c, b)| c * b).sum();
                    if sh_active && !allow_artificial {
                        let scale = dot.abs().max(cbtb.abs()).max(1.0);
                        if (dot - cbtb).abs() > 1e-6 * scale {
                            return Ok((iters, PhaseOutcome::NumericalDistress));
                        }
                    }
                    if std::env::var("RUSOLVE_DEBUG").is_ok() {
                        let max_bbar = bbar.iter().fold(0.0f64, |a, &v| a.max(v));
                        eprintln!(
                            "rusolve simplex: RAY-DECLARE q={q} d={_d} pi_dot={dot} \
                             cbT_bbar={cbtb} max_bbar={max_bbar}"
                        );
                    }
                    if p1trace {
                        tr.emit(iters, use_bland, p1_perturb);
                    }
                    // The ray that survived the fresh-factor verification
                    // above — hand its FTRAN'd column out, don't recompute.
                    exit.pi.clone_from(&pi);
                    exit.ray = Some((q, sgn, bbar.clone()));
                    return Ok((iters, PhaseOutcome::Unbounded));
                }
                None => {
                    if std::env::var("RUSOLVE_DEBUG").is_ok() {
                        eprintln!(
                            "[dbg] it={iters} q={q} from_ub={from_ub} d={_d} theta={theta} FLIP (no basis change)"
                        );
                    }
                    // Entering saturates at its own ub: θ is applied to the
                    // basic values as usual (θ is below every basic cap by
                    // construction), but the basis itself is unchanged — the
                    // entering variable just flips to its other bound.
                    for &i in nzb.iter() {
                        x_b[i] -= theta * bbar[i];
                        if x_b[i] < 0.0 && x_b[i] > -1e-7 {
                            x_b[i] = 0.0;
                        }
                    }
                    at_ub[q] = !from_ub;
                    {
                        let on = is_cand(st, allow_artificial, dead[q], at_ub[q], dvec[q], ptol, q);
                        let bit = 1u64 << (q & 63);
                        if on {
                            cand_bits[q >> 6] |= bit;
                        } else {
                            cand_bits[q >> 6] &= !bit;
                        }
                    }
                    // FLIP is progress (θ>0), but decays rather than clears:
                    // phase-1 flips can cycle at zero w-change.
                    stall = stall.saturating_sub(1);
                    if bland_left > 0 {
                        bland_clear += 1; // episode release-progress counter
                    }
                }
                Some(pos) => {
                    if std::env::var("RUSOLVE_DEBUG").is_ok() {
                        eprintln!(
                            "[dbg] it={iters} q={q} from_ub={from_ub} d={_d} theta={theta} leave_pos={pos} leaving={} leave_to_ub={leave_to_ub} x_b={:?}",
                            basis[pos], x_b
                        );
                    }
                    if theta < 1e-11 {
                        // Exact-zero AND near-zero (fp noise on degenerate
                        // caps, e.g. 4e-25) count as degenerate (scrs8).
                        stall += 1;
                    } else {
                        // Decay, not reset: sporadic degeneracy in big LPs
                        // oscillates below the engage threshold, while a true
                        // cycle (all-θ0 stretches) accumulates monotonically.
                        stall = stall.saturating_sub(1);
                        if bland_left > 0 {
                            bland_clear += 1; // episode release-progress counter
                        }
                    }
                    let leaving_var = basis[pos];
                    if allow_artificial && st.kind[leaving_var] == Kind::Artificial {
                        dead[leaving_var] = true; // artificial deletion (see dead decl)
                    }

                    // Primal values: x_B ← x_B − θ·b̄, entering slot gets its own value.
                    for &i in nzb.iter() {
                        x_b[i] -= theta * bbar[i];
                        if x_b[i] < 0.0 && x_b[i] > -1e-7 {
                            x_b[i] = 0.0; // round-off clamp
                        } else {
                            // Relaxed phase-1 mirror: a relaxed step may land a
                            // crossed row just past its ub; absorb the
                            // round-off-sized overshoot exactly like the
                            // negative clamp above. (Strict steps never
                            // trigger this — caps bound the value at ub.)
                            let ub = st.ub[basis[i]];
                            if x_b[i] > ub && x_b[i] < ub + 1e-7 {
                                x_b[i] = ub;
                            }
                        }
                    }

                    // ---- SH weight maintenance ----
                    // EXACT incremental update, algebra transcribed from
                    // HiGHS HEkkPrimal::updatePrimalSteepestEdgeWeights (MIT;
                    // formulas only; Forrest & Goldfarb 1992). Invariant w_j = 1 + ‖â_j‖² (full
                    // tableau column ⇒ w ≥ 1, no division blow-up):
                    //   â_pj = (B⁻¹a_j)[pos]  via row_p = B⁻ᵀe_pos   (BTRAN 1)
                    //   μ    = B⁻ᵀâ_q          (PSE BTRAN)             (BTRAN 2)
                    //   μᵀa_j = â_qᵀ(B⁻¹a_j) = μ·a_j                  (cross term)
                    //   λ_j  = â_pj / α_pq
                    //   w_j ← max(1+λ_j², w_j + λ_j²(s²+1) − 2λ_j·μᵀa_j)
                    //   w_leaving = (1+s²)/α_pq²  (exact; becomes nonbasic)
                    //   w_entering = 0            (basic ⇒ never priced)
                    // in_basis still describes the OLD basis here (leaving is
                    // basic, q nonbasic), so the loop skips the leaving col
                    // (set exactly below) and q; both BTRANs precede
                    // pfi.update, while the factors still describe old B.
                    // Pivotal row rho = B^-T e_pos (BTRAN 1, shared by the
                    // reduced-cost update and the SH weights), priced through
                    // the row-wise A: only columns hit by rho's nonzeros.
                    e_pos.fill(0.0);
                    e_pos[pos] = 1.0;
                    pfi.solve_bt(&e_pos, &mut pi_pos);
                    let pivot = sgn * bbar[pos];
                    touched.clear();
                    for (i, &r) in pi_pos.iter().enumerate() {
                        if r == 0.0 {
                            continue;
                        }
                        for k in rptr[i]..rptr[i + 1] {
                            let c = rcol[k];
                            if in_basis[c] {
                                continue;
                            }
                            if !row_mark[c] {
                                row_mark[c] = true;
                                touched.push(c);
                            }
                            row_alpha[c] += r * rval[k];
                        }
                    }
                    // Numerical-trouble check (HiGHS HEkkPrimal updatePivots
                    // style, MIT, mechanism only): the pivot seen through the
                    // column (FTRAN) and through the row (BTRAN) must agree.
                    // A mismatch means the eta chain has lost accuracy (long
                    // chains: fewer BTRANs per iteration also slow the
                    // synthetic refactor clock). Refactorize, re-anchor x_B
                    // and redo this iteration; with a fresh chain proceed.
                    let arq = row_alpha[q];
                    if trouble_check
                        && pfi.eta_count() > 0
                        && (pivot - arq).abs() > 1e-7 * pivot.abs().max(1.0)
                    {
                        for &c in touched.iter() {
                            row_alpha[c] = 0.0;
                            row_mark[c] = false;
                        }
                        // Wall-time distress breaker: time THIS forced
                        // rebuild; a basis whose distress arms keep paying
                        // real seconds is numerically dead — same verdict
                        // class as the |x_b| detonation above (phase 1 →
                        // dirty-exit rung: fresh re-solve with the
                        // unconditional exact-tie stability preference;
                        // phase 2 has no rung: honest Limit). HiGHS bounds
                        // its simplex the same way on numerical distress
                        // (`HEkk::bailout` + non-fixed status,
                        // reference/HiGHS/highs/simplex/HEkk.cpp:3441–3461,
                        // MIT).
                        distress_events += 1;
                        let t_d = std::time::Instant::now();
                        let rf = pfi.refactorize();
                        distress_secs += t_d.elapsed().as_secs_f64();
                        if distress_events >= 3 && distress_secs >= 2.0 {
                            let outcome = if allow_artificial {
                                PhaseOutcome::Detonated
                            } else {
                                PhaseOutcome::Limit
                            };
                            eprintln!(
                                "rusolve simplex: numerical-trouble storm ({distress_events} distress refactorizations, {distress_secs:.2}s spent rebuilding at it={iters}) — {}",
                                if allow_artificial {
                                    "fresh-factor re-solve rung"
                                } else {
                                    "phase limit"
                                }
                            );
                            return Ok((iters, outcome));
                        }
                        match rf {
                            Ok(()) => {
                                refresh_x_b(pfi, st, basis, at_ub, x_b);
                                d_ok = false;
                                iters -= 1; // a re-anchor pass is not a pivot
                                continue;
                            }
                            Err(FactorError::Singular { column, .. }) if allow_artificial => {
                                if repair_basis(
                                    st,
                                    pfi,
                                    basis,
                                    x_b,
                                    at_ub,
                                    column,
                                    self.config.pfi,
                                    self.deadline,
                                ) {
                                    weights_stale = true;
                                } else {
                                    trouble_check = false;
                                }
                                d_ok = false;
                                iters -= 1;
                                continue;
                            }
                            // Rank-deficient and unrepairable: the eta chain is
                            // all we have (scsd1 lesson); stop second-guessing
                            // it for the rest of this phase and redo the pivot.
                            Err(_) => {
                                trouble_check = false;
                                d_ok = false;
                                iters -= 1;
                                continue;
                            }
                        }
                    }
                    let theta_d = dq_true / pivot;
                    for &c in touched.iter() {
                        if c != q {
                            dvec[c] -= theta_d * row_alpha[c];
                            let on =
                                is_cand(st, allow_artificial, dead[c], at_ub[c], dvec[c], ptol, c);
                            let bit = 1u64 << (c & 63);
                            if on {
                                cand_bits[c >> 6] |= bit;
                            } else {
                                cand_bits[c >> 6] &= !bit;
                            }
                        }
                    }
                    dvec[leaving_var] = -theta_d;
                    dvec[q] = 0.0;
                    cb[pos] = cost[q];
                    if sh_active {
                        let mut s2 = 0.0;
                        alpha_raw.fill(0.0);
                        for &i in nzb.iter() {
                            let a = if from_ub { -bbar[i] } else { bbar[i] };
                            alpha_raw[i] = a;
                            s2 += a * a;
                        }
                        pfi.solve_bt(&alpha_raw, &mut mu); // mu = B^-T a_q (PSE)
                        for &v in touched.iter() {
                            if v == q {
                                continue;
                            }
                            // Columns with alpha_rv == 0 keep w unchanged
                            // (lam = 0 => wn = w >= 1 = floor).
                            let pap = row_alpha[v];
                            let (lo, hi) = (st.csc_ptr[v], st.csc_ptr[v + 1]);
                            let mut muj = 0.0;
                            for k in lo..hi {
                                muj += mu[st.csc_row[k]] * st.csc_val[k];
                            }
                            let lam = pap / pivot;
                            let wn = w[v] + lam * lam * (s2 + 1.0) - 2.0 * lam * muj;
                            let floor = 1.0 + lam * lam;
                            w[v] = if wn.is_finite() && wn > floor {
                                wn
                            } else {
                                floor
                            };
                        }
                        w[leaving_var] = (1.0 + s2) / (pivot * pivot);
                        w[q] = 0.0;
                    }
                    for &c in touched.iter() {
                        row_alpha[c] = 0.0;
                        row_mark[c] = false;
                    }

                    // Factorization update: the eta is built from the
                    // caller's b̄ — the very FTRAN the ratio test just
                    // consumed — instead of a duplicated internal solve.
                    // `bbar` holds B⁻¹(sgn·a_q); passing `sgn` UN-applies the
                    // direction sign so the eta describes the RAW column a_q
                    // that the basis gets (exact ±1 flip ⇒ bit-identical to
                    // the deleted internal FTRAN; see `update`'s docs).
                    // The entering column is the gathered dense vector
                    // Same bits the old `&st.cols[q]` slice held.
                    let q_nz = st.sparse_col(q);
                    if let Err(FactorError::RefactorizeNeeded { .. }) =
                        pfi.update_sparse(pos, q_nz, &bbar, &nzb, sgn)
                    {
                        // Wall-time distress breaker (same condition as
                        // the cross-check arm above): a refused pivot pays a
                        // fresh O(m³) rebuild here; seconds of accumulated
                        // distress rebuilds = dead basis.
                        distress_events += 1;
                        let t_d = std::time::Instant::now();
                        let rf = pfi.refactorize();
                        distress_secs += t_d.elapsed().as_secs_f64();
                        if distress_events >= 3 && distress_secs >= 2.0 {
                            let outcome = if allow_artificial {
                                PhaseOutcome::Detonated
                            } else {
                                PhaseOutcome::Limit
                            };
                            eprintln!(
                                "rusolve simplex: numerical-trouble storm ({distress_events} distress refactorizations, {distress_secs:.2}s spent rebuilding at it={iters}, update-refused) — {}",
                                if allow_artificial {
                                    "fresh-factor re-solve rung"
                                } else {
                                    "phase limit"
                                }
                            );
                            return Ok((iters, outcome));
                        }
                        rf?;
                        st.gather_col(q, 1.0, &mut qcol);
                        // The retry runs on FRESH factors, where the caller's
                        // b̄ (computed on the stale eta chain) no longer
                        // applies: recompute it on the rebuilt factorization,
                        // exactly what the former internal solve did here.
                        pfi.refactorize()?;
                        let mut bbar_fresh = vec![0.0; m];
                        pfi.solve_b(&qcol, &mut bbar_fresh);
                        pfi.update(pos, &qcol, &bbar_fresh, 1.0)?;
                    }
                    basis[pos] = q;
                    x_b[pos] = if from_ub { st.ub[q] - theta } else { theta };
                    at_ub[q] = false;
                    at_ub[leaving_var] = leave_to_ub; // leaves at lb (false) or ub (true)
                    cand_bits[q >> 6] &= !(1u64 << (q & 63));
                    {
                        let l = leaving_var;
                        let on = is_cand(st, allow_artificial, dead[l], at_ub[l], dvec[l], ptol, l);
                        let bit = 1u64 << (l & 63);
                        if on {
                            cand_bits[l >> 6] |= bit;
                        } else {
                            cand_bits[l >> 6] &= !bit;
                        }
                    }
                    ray_verified = false; // real pivot ⇒ guard available again
                    p1_exit_verified = false; // next exit re-verifies
                    if repair_head.take().is_some() {
                        // Repair pivot done: re-anchor NOW — the leaver's
                        // bound lie (≤ its violation) must land on the fresh
                        // factors so the residual flows into the artificials
                        // instead of poisoning the next pricing.
                        //
                        // Wall-time distress breaker (measured storm site:
                        // 172 of the crawl's 173 dense rebuilds): on mik's
                        // restart round 2 the basis is rank-deficient enough
                        // that every rebuild exposes a new dependency —
                        // repair pivot, re-anchor, still singular, next
                        // repair — a cascade of numerically-forced pivots,
                        // each paying a full sticky-dense O(m³) rebuild
                        // (~0.35 s at m=1 040, 98% of wall in DenseLu::
                        // factor) while committing
                        // "progress" no stall ladder can price. A repair is
                        // distress by definition: it feeds the same seconds
                        // accumulator as the cross-check and refusal arms.
                        distress_events += 1;
                        let t_d = std::time::Instant::now();
                        let rf = pfi.refactorize();
                        distress_secs += t_d.elapsed().as_secs_f64();
                        if distress_events >= 3 && distress_secs >= 2.0 {
                            let outcome = if allow_artificial {
                                PhaseOutcome::Detonated
                            } else {
                                PhaseOutcome::Limit
                            };
                            eprintln!(
                                "rusolve simplex: numerical-trouble storm ({distress_events} distress refactorizations, {distress_secs:.2}s spent rebuilding at it={iters}, repair cascade) — {}",
                                if allow_artificial {
                                    "fresh-factor re-solve rung"
                                } else {
                                    "phase limit"
                                }
                            );
                            return Ok((iters, outcome));
                        }
                        if rf.is_ok() {
                            refresh_x_b(pfi, st, basis, at_ub, x_b);
                        }
                        d_ok = false;
                    }
                    // SH-AUTOPSY (debug only): maintained weights vs
                    // ground truth for the NEW basis — run AFTER pfi.update
                    // so both sides describe the same basis (the pre-update
                    // check compared new-basis maintained values against
                    // old-basis FTRANs: a bookkeeping artifact, not damage).
                    // Repair on detection ⇒ every report is fresh damage.
                    if sh_active
                        && std::env::var("RUSOLVE_DEBUG").is_ok()
                        && sh_report_left > 0
                        && iters - start_iters < 6000
                    {
                        let stride = (total / 4).max(1);
                        let mut v = stride.saturating_sub(1);
                        while v < total && sh_report_left > 0 {
                            if !in_basis[v] && v != q {
                                let mut t = vec![0.0; m];
                                st.gather_col(v, 1.0, &mut qcol);
                                pfi.solve_b(&qcol, &mut t);
                                let true_w: f64 = 1.0 + t.iter().map(|x| x * x).sum::<f64>();
                                let rel = (true_w - w[v]).abs() / (true_w + w[v] + 1e-300);
                                if rel > 1e-6 {
                                    eprintln!(
                                        "rusolve SH-AUTOPSY it={iters} q={q} from_ub={from_ub} \
                                         pos={pos} v={v} w={:.6e} true={:.6e} rel={rel:.1e}",
                                        w[v], true_w,
                                    );
                                    w[v] = true_w;
                                    sh_report_left -= 1;
                                }
                            }
                            v += stride;
                        }
                    }
                }
            }
            if p1trace {
                // TEMP diagnostic: count the completed pivot (FLIP arm too,
                // with its theta); refresh w_true every 64th pivot; emit one
                // stderr line every 1000 pivots.
                tr.count_pivot(theta, st.kind[q] == Kind::Artificial);
                if tr.n_piv.is_multiple_of(64) {
                    tr.recompute_w(st, basis, x_b);
                }
                if tr.n_piv - tr.last_emit >= 1000 {
                    tr.emit(iters, use_bland, p1_perturb);
                }
            }
        }
    }
}

// ---- Phase-1 sub-tolerance entering tie perturbation ----
/// Deterministic per-column uniform hash in `[0,1)`: the splitmix64
/// finalizer (Steele/Lea/Vigna reference splitmix64, public domain) applied
/// to the column index with a fixed seed offset. Unlike a sequential RNG
/// stream, `hash01(v)` depends on `v` alone ⇒ the eps vector is
/// draw-order-independent and identical across runs and machines.
pub(crate) fn hash01(v: u64) -> f64 {
    let mut z = v.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    // Top 53 bits ⇒ exact multiple of 2^-53, values in [0,1).
    (z >> 11) as f64 / (1u64 << 53) as f64
}

/// eps[v] = ptol · 0.5 · hash01(v + seed) ∈ [0, ptol/2): SUB-tolerance cost
/// perturbation, so the unperturbed improver test `d < -ptol` can never
/// disagree with the perturbed ranking's notion of who is improving by a
/// real amount (an improver's key stays < -ptol/2; Beale 1968 / Chvátal
/// 1983 §7.3 lineage, entering ties only). `seed` makes the
/// vector regenerable — seed 0 is the original one-shot vector bit-exactly, a
/// fresh seed is a fresh deterministic draw (see SimplexConfig::
/// p1_re_perturb).
fn p1_tie_eps(total: usize, ptol: f64, seed: u64) -> Vec<f64> {
    (0..total)
        .map(|v| ptol * 0.5 * hash01((v as u64).wrapping_add(seed)))
        .collect()
}

/// Phase-1 progress probe — the maintained basic-artificial sum (the
/// phase-1 objective) and the basic-artificial count. O(m), read-only; the
/// re-perturbation freeze test compares these BITWISE across stall windows
/// (a θ=0 crawl leaves them bit-identical; any real θ>0 progress breaks
/// them). Same summation order as the phase-1-exit `p1_obj` check.
fn artificial_sum(st: &Standard, basis: &[usize], x_b: &[f64]) -> (f64, usize) {
    let mut w = 0.0f64;
    let mut n = 0usize;
    for (i, &v) in basis.iter().enumerate() {
        if st.kind[v] == Kind::Artificial {
            w += x_b[i];
            n += 1;
        }
    }
    (w, n)
}

/// Relaxed phase 1: total phase-1 infeasibility = artificial-row |values|
/// (row residuals) PLUS structural basic bound violations. Equal to
/// `artificial_sum` whenever every basic is in range (the invariant of
/// non-relaxed phase 1), but unlike it stays honest when a relaxed step
/// has crossed a structural bound (HiGHS: simplex/HEkkPrimal.cpp:phase1ChooseRow).
fn total_infeas(st: &Standard, basis: &[usize], x_b: &[f64]) -> (f64, usize) {
    let mut w = 0.0;
    let mut n = 0usize;
    for (i, &v) in basis.iter().enumerate() {
        if st.kind[v] == Kind::Artificial {
            w += x_b[i].abs();
            n += 1;
        } else {
            if x_b[i] < 0.0 {
                w -= x_b[i];
            }
            let up = x_b[i] - st.ub[v];
            if up.is_finite() && up > 0.0 {
                w += up;
            }
        }
    }
    (w, n)
}

// ---- TEMPORARY phase-1 diagnostic (env-gated: RUSOLVE_P1TRACE) ----
// Pivot-step histogram + basic-artificial w recheck for phase-1 degeneracy
// analysis. Counters only touched when the env var is set ⇒ no behavior
// change otherwise. A later step decides whether to keep or revert this.
#[derive(Default)]
struct P1Trace {
    n_piv: usize,
    last_emit: usize,
    hist: [usize; 5],
    piv_since: usize,
    art_enter: usize,
    deg_piv: usize,
    deg_ties_sum: usize,
    deg_ties_with_art: usize,
    w_true: f64,
    w_prev: f64,
    dw: f64,
    n_art_basic: usize,
}

impl P1Trace {
    fn bucket(theta: f64) -> usize {
        if theta < 1e-15 {
            0
        } else if theta < 1e-11 {
            1
        } else if theta < 1e-8 {
            2
        } else if theta < 1e-4 {
            3
        } else {
            4
        }
    }
    fn count_pivot(&mut self, theta: f64, art_entering: bool) {
        self.n_piv += 1;
        self.piv_since += 1;
        self.hist[Self::bucket(theta)] += 1;
        if art_entering {
            self.art_enter += 1;
        }
    }
    fn count_degenerate(&mut self, ties: usize, any_art: bool) {
        self.deg_piv += 1;
        self.deg_ties_sum += ties;
        if any_art {
            self.deg_ties_with_art += 1;
        }
    }
    fn recompute_w(&mut self, st: &Standard, basis: &[usize], x_b: &[f64]) {
        let mut w = 0.0f64;
        let mut n = 0usize;
        for (i, &v) in basis.iter().enumerate() {
            if st.kind[v] == Kind::Artificial {
                w += x_b[i];
                n += 1;
            }
        }
        if self.w_true.is_nan() {
            self.w_prev = w; // first sample: no prior w to diff against
        } else {
            self.w_prev = self.w_true;
        }
        self.w_true = w;
        self.n_art_basic = n;
        self.dw = self.w_true - self.w_prev;
    }
    fn emit(&mut self, it: usize, use_bland: bool, p1eps: bool) {
        let h = self.hist;
        let avg = if self.deg_piv > 0 {
            self.deg_ties_sum as f64 / self.deg_piv as f64
        } else {
            0.0
        };
        let frac = if self.deg_piv > 0 {
            self.deg_ties_with_art as f64 / self.deg_piv as f64
        } else {
            0.0
        };
        eprintln!(
            "p1trace it={it} w_true={:.6e} dw_last64={:.3e} n_art_basic={} \
             hist=[<1e-15:{},<1e-11:{},<1e-8:{},<1e-4:{},>=1e-4:{}] (since last line) \
             deg_ties_avg={avg:.2} deg_ties_with_art_frac={frac:.2} \
             art_entering={}/{} bland={} p1eps={}",
            self.w_true,
            self.dw,
            self.n_art_basic,
            h[0],
            h[1],
            h[2],
            h[3],
            h[4],
            self.art_enter,
            self.piv_since,
            use_bland as u8,
            p1eps as u8,
        );
        self.hist = [0; 5];
        self.piv_since = 0;
        self.art_enter = 0;
        self.deg_piv = 0;
        self.deg_ties_sum = 0;
        self.deg_ties_with_art = 0;
        self.last_emit = self.n_piv;
    }
}

#[allow(clippy::too_many_arguments)]
/// Repair a numerically singular basis: replace the dependent basis
/// column at slot `bad_pos` with the first non-basic artificial whose
/// refactorization succeeds AND yields a primal-feasible x_B = B⁻¹b ≥ 0.
/// Artificial columns are never deleted, so the fallback pool always exists;
/// phase 1 re-enters the artificial at its true value and continues from a
/// genuinely feasible BFS (w may rise — correctness over monotonicity).
/// Provenance: artificial-in-basis restoration is standard crash-repair
/// practice (BASP ch. 9 / LUSOL REFCTOR semantics); HiGHS reports the same
/// "restoring artificials" step in its phase-1 log.
/// Column `v` of `st` as the ascending `(row, value)` list that
/// `gather_col(v, 1.0, ·)` followed by dropping exact zeros yields
/// (duplicate CSC rows: last write wins, as in `gather_col`).
fn sparse_col(st: &Standard, v: usize) -> Vec<(usize, f64)> {
    let mut c: Vec<(usize, f64)> = (st.csc_ptr[v]..st.csc_ptr[v + 1])
        .map(|k| (st.csc_row[k], st.csc_val[k]))
        .collect();
    c.sort_by_key(|&(i, _)| i); // stable: duplicate rows keep CSC order
    let mut out: Vec<(usize, f64)> = Vec::with_capacity(c.len());
    for (i, x) in c {
        match out.last_mut() {
            Some(last) if last.0 == i => last.1 = x,
            _ => out.push((i, x)),
        }
    }
    out.retain(|&(_, x)| x != 0.0);
    out
}

#[allow(clippy::too_many_arguments)]
fn repair_basis(
    st: &Standard,
    pfi: &mut ProductFormFactorization,
    basis: &mut [usize],
    x_b: &mut Vec<f64>,
    at_ub: &mut [bool],
    bad_pos: usize,
    cfg: PfiConfig,
    deadline: Option<std::time::Instant>,
) -> bool {
    let m = st.m;
    let art_start = st.ncol() - m;
    // Sparse basis columns, gathered ONCE. The old rebuild gathered
    // m dense length-m columns plus a DenseMatrix copy PER CANDIDATE —
    // O(m²) memory on every attempt (cont1 at m = 160,792 is ~200 GB per
    // copy). `sparse_col` reproduces
    // `gather_col` + `sparsify` exactly (ascending rows, exact zeros
    // dropped), so `build_sparse` sees the columns `build` saw.
    let in_basis_set: std::collections::HashSet<usize> = basis.iter().copied().collect();
    let base_cols: Vec<Vec<(usize, f64)>> = basis.iter().map(|&v| sparse_col(st, v)).collect();
    for r in 0..m {
        let cand = art_start + r;
        if in_basis_set.contains(&cand) {
            continue;
        }
        // Each candidate is a full
        // factorization, so the scan is O(m) of them — ex10's MIP root LP
        // (m ≈ 63k) sat here for minutes past its deadline. One clock read
        // per factorization is noise; an expired deadline gives up the
        // repair and the pivot loop's own deadline check returns Limit.
        if deadline.is_some_and(|d| std::time::Instant::now() >= d) {
            return false;
        }
        let mut cols = base_cols.clone();
        cols[bad_pos] = sparse_col(st, cand);
        if let Ok(new_pfi) = ProductFormFactorization::build_sparse(cols, cfg) {
            let mut xb_try = vec![0.0; m];
            new_pfi.solve_b(&st.b, &mut xb_try);
            if !xb_try.iter().all(|&v| v >= -1e-9) {
                continue; // factorizable but primal-infeasible swap: next
            }
            at_ub[basis[bad_pos]] = false; // dependent column leaves at lb
            basis[bad_pos] = cand;
            *pfi = new_pfi;
            *x_b = xb_try;
            for xv in x_b.iter_mut() {
                if *xv < 0.0 && *xv > -1e-7 {
                    *xv = 0.0;
                }
            }
            eprintln!("rusolve simplex: basis repaired (artificial {cand} into slot {bad_pos})");
            return true;
        }
    }
    false
}

// ---- Caller-space dual & ray exports ----

/// Mirror of the standard-form construction's row-negation predicate: row `i`
/// is negated iff the bound-shifted rhs is negative. Deliberately a standalone
/// mirror (a `Standard` field would sit inside the construction region
/// that is kept separate) — the b-shift below is the same CSC scatter
/// `to_standard_sparse` performs (ascending j, per-column ascending rows,
/// nonzeros only), so the verdicts are bit-identical. If the construction's
/// negation rule ever changes, THIS MUST CHANGE WITH IT.
fn negated_rows(lp: &Lp) -> Vec<bool> {
    let mut b: Vec<f64> = lp.rhs.clone();
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
    for j in 0..lp.num_col {
        let s = shift_of(lp, j);
        if s != 0.0 {
            for k in ptr[j]..ptr[j + 1] {
                let v = rval[k];
                if v != 0.0 {
                    b[ridx[k] as usize] -= v * s;
                }
            }
        }
    }
    b.iter().map(|&bi| bi < 0.0).collect()
}

/// Map the phase-2 exit BTRAN duals π (standard space, costs sign-folded by
/// `sign`) to the caller convention documented on
/// [`SimplexSolution::reduced_costs`]: undo the row negations, multiply by
/// `sign` (= −1 for maximize: phase 2 minimises −c, so π scales with it),
/// then read `d_j = c_j − πᵀa_j` off the ORIGINAL columns via
/// [`reduced_costs_from`] — one dense pass, which is exactly why flip/split
/// columns need no per-column sign surgery.
fn caller_duals(lp: &Lp, sign: f64, pi: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let neg = negated_rows(lp);
    let row_duals: Vec<f64> = pi
        .iter()
        .enumerate()
        .map(|(i, &p)| if neg[i] { -sign * p } else { sign * p })
        .collect();
    let reduced_costs = reduced_costs_from(lp, &row_duals);
    (row_duals, reduced_costs)
}

/// The reduced-cost half of the caller convention: for a CALLER-space row
/// dual vector π (already un-negated and sign-mapped), read
/// `d_j = c_j − πᵀa_j` off the ORIGINAL columns and costs — the single
/// identity that defines the reported pair `(π, d)`. Shared by the primal
/// exit (through [`caller_duals`]) and the dual engine's exit, whose standard
/// row space IS the caller's (no negation to undo; the objective-sign factor
/// is folded into the π it passes here).
pub(crate) fn reduced_costs_from(lp: &Lp, row_duals: &[f64]) -> Vec<f64> {
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
    (0..lp.num_col)
        .map(|j| {
            // πᵀa_j straight off the CSC column (ascending rows, the dense
            // loop's order; stored zeros are ±0.0 no-ops).
            let dot: f64 = (ptr[j]..ptr[j + 1])
                .filter(|&k| rval[k] != 0.0)
                .map(|k| row_duals[ridx[k] as usize] * rval[k])
                .sum();
            lp.objective[j] - dot
        })
        .collect()
}

/// Assemble the caller-space unbounded direction from the verified
/// exit: Δ = sgn on the entering column, −b̄ᵢ on the basic columns (the
/// ratio test found no blocker on THIS direction — the entering ub is +inf
/// and every b̄ < −tol row has ub = +inf), 0 elsewhere; then undo the bound
/// transforms (shifts add nothing to a direction, flip negates, split
/// subtracts). An all-zero result (a pure slack/artificial direction: x
/// itself cannot run away) reports `None`.
fn primal_ray_x(
    lp: &Lp,
    st: &Standard,
    basis: &[usize],
    q: usize,
    sgn: f64,
    bbar: &[f64],
) -> Option<Vec<f64>> {
    let mut dcol = vec![0.0; st.ncol()];
    dcol[q] = sgn;
    for (i, &v) in basis.iter().enumerate() {
        dcol[v] -= bbar[i];
    }
    let mut d = vec![0.0; lp.num_col];
    // j indexes x, cols_of and flip jointly (as in extract_x).
    #[allow(clippy::needless_range_loop)]
    for j in 0..lp.num_col {
        d[j] = match st.cols_of[j].as_slice() {
            [c] => {
                if st.flip[j] {
                    -dcol[*c]
                } else {
                    dcol[*c]
                }
            }
            [c1, c2] => dcol[*c1] - dcol[*c2],
            _ => 0.0,
        };
    }
    if d.iter().all(|&v| v == 0.0) {
        None
    } else {
        Some(d)
    }
}

/// Phase-1 exit duals (w > 0) → caller-space Farkas candidate: undo the row
/// negations (no objective-sign factor — phase 1 always minimises the
/// artificial sum; the bound shifts cancel exactly between πᵀb̃ and the
/// column sup, so the certificate is stated in the ORIGINAL rhs space).
/// Returned only if [`farkas_ok`] verifies it against the original rows;
/// the debug build ALSO asserts the verification (a derivation bug must be
/// loud in tests), the release build reports the honest `None`.
fn attempt_farkas(lp: &Lp, pi: &[f64]) -> Option<Vec<f64>> {
    if pi.len() != lp.num_row {
        return None;
    }
    let neg = negated_rows(lp);
    let y: Vec<f64> = pi
        .iter()
        .zip(&neg)
        .map(|(&p, &n)| if n { -p } else { p })
        .collect();
    debug_assert!(
        farkas_ok(lp, &y),
        "phase-1 exit duals must verify as a Farkas certificate: {y:?}"
    );
    farkas_ok(lp, &y).then_some(y)
}

/// The Farkas identity of [`SimplexSolution::dual_ray`] as a numerical
/// check on the ORIGINAL rows/columns/bounds. The row sup is 0 under the
/// row sign conditions ('L': y ≤ 0, 'G': y ≥ 0, 'E': free — the 'L'/'G'
/// conditions are also FINITENESS conditions: a violated sign makes the
/// row-slack sup +∞), so the whole sup collapses to a per-column term
/// `max { φ·x_j : x_j ∈ [lb, ub] }` with φ = yᵀa_j. THAT max lives at the
/// endpoint the sign of φ selects — `φ > 0` ⇒ `ub`, `φ < 0` ⇒ `lb`, `φ = 0` ⇒
/// any — and it is `+∞` (no certificate) only when THAT side is unbounded.
/// Picking the max over whichever endpoints happen to be finite instead would
/// insert a phantom `0.0` for the infinite side and OVERESTIMATE the sup
/// (e.g. `φ = -1, lb = 2, ub = +∞`: true `-2`, phantom `max(-2, 0) = 0`),
/// rejecting valid certificates on LPs with bounded variables (regression test
/// `farkas_ray_survives_bounded_variables`). Only inside the pricing-noise
/// band `|φ| ≤ tol` does the estimate stay endpoint-clamped: there the error
/// is `tol·|endpoint|`, i.e. below the tolerance the identity is priced with.
/// Tolerances are priced-noise-sized (the certificate is born from BTRAN π
/// with |d| ≤ pricing_tol at exit); a gap that cannot clear them reports
/// `None` — an honest None beats a wrong ray.
pub(crate) fn farkas_ok(lp: &Lp, y: &[f64]) -> bool {
    let y_inf = y.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
    let tol = 1e-6 * (1.0 + y_inf);
    let mut lhs = 0.0f64;
    for (i, (&s, &yi)) in lp.senses.iter().zip(y).enumerate() {
        match s {
            'L' if yi > tol => return false,
            'G' if yi < -tol => return false,
            _ => {}
        }
        lhs += yi * lp.rhs[i];
    }
    let mut sup = 0.0f64;
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
    for j in 0..lp.num_col {
        let phi: f64 = (ptr[j]..ptr[j + 1])
            .filter(|&k| rval[k] != 0.0)
            .map(|k| y[ridx[k] as usize] * rval[k])
            .sum();
        let (lb, ub) = (lp.lower[j], lp.upper[j]);
        if phi > tol {
            // Maximized at the upper endpoint: infinite ⇒ no certificate.
            if ub == f64::INFINITY {
                return false;
            }
            sup += phi * ub;
        } else if phi < -tol {
            // Maximized at the lower endpoint (the phantom-zero case).
            if lb == f64::NEG_INFINITY {
                return false;
            }
            sup += phi * lb;
        } else {
            // |φ| ≤ pricing noise: endpoint-clamped estimate, error ≤
            // tol·max(|lb|, |ub|) per column — under the tolerance below.
            let at_lb = if lb.is_finite() { phi * lb } else { 0.0 };
            let at_ub = if ub.is_finite() { phi * ub } else { 0.0 };
            sup += at_lb.max(at_ub);
        }
    }
    let scale = lhs.abs().max(sup.abs()).max(1.0);
    lhs > sup + 1e-9 * scale
}

#[allow(clippy::too_many_arguments)]
/// Iteration-limit snapshot builder: mirrors SimplexSolution's field set 1:1,
/// so the argument list is the struct's shape, not accidental complexity.
fn limit_solution(
    st: &Standard,
    objective: &[f64],
    basis: &[usize],
    x_b: &[f64],
    at_ub: &[bool],
    p1: usize,
    p2: usize,
    bland: bool,
) -> SimplexSolution {
    let x = extract_x(st, basis, x_b, at_ub);
    let objective: f64 = objective.iter().zip(&x).map(|(c, xv)| c * xv).sum();
    SimplexSolution {
        status: SimplexStatus::IterationLimit,
        objective,
        x,
        phase1_iterations: p1,
        phase2_iterations: p2,
        bland_triggered: bland,
        // An unfinished search has no valid duals and no rays.
        ..Default::default()
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    /// Dense columns of an LP for verification code (dense allowed in
    /// tests) — stored zeros land at their positions.
    fn dense_cols(lp: &Lp) -> Vec<Vec<f64>> {
        let mut cols = vec![vec![0.0; lp.num_row]; lp.num_col];
        let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
        for j in 0..lp.num_col {
            for k in ptr[j]..ptr[j + 1] {
                cols[j][ridx[k] as usize] = rval[k];
            }
        }
        cols
    }

    fn lp_of(rows: &[&[f64]], rhs: Vec<f64>, senses: &str, obj: Vec<f64>, ub: Vec<f64>) -> Lp {
        let n = obj.len();
        let m = rhs.len();
        // Dense test data → CSC triplets, zeros stored like the io reader
        // keeps them (`SparseLp::from_lp` then drops them, I3). `take(n)`
        // keeps the dense cm's silent truncation of over-long rows.
        let triplets: Vec<(usize, usize, f64)> = rows
            .iter()
            .enumerate()
            .flat_map(|(i, row)| row.iter().enumerate().take(n).map(move |(j, &v)| (i, j, v)))
            .collect();
        Lp {
            num_col: n,
            num_row: m,
            a: rusolve_core::CscMatrix::from_triplets(m, n, &triplets),
            rhs,
            objective: obj,
            maximize: false,
            senses: senses.chars().collect(),
            lower: vec![0.0; n],
            upper: ub,
        }
    }

    fn solve(
        rows: &[&[f64]],
        rhs: Vec<f64>,
        senses: &str,
        obj: Vec<f64>,
        ub: Vec<f64>,
    ) -> SimplexSolution {
        SimplexSolver::new(SimplexConfig::default())
            .solve(&lp_of(rows, rhs, senses, obj, ub))
            .unwrap()
    }

    /// Assert two solutions are bit-identical (status, objective bits, full x
    /// bits, per-phase iteration counts) — the A/B spine of the sparse input path.
    fn assert_solutions_bit_identical(got: &SimplexSolution, want: &SimplexSolution, what: &str) {
        assert_eq!(got.status, want.status, "{what}: status");
        assert_eq!(
            got.objective.to_bits(),
            want.objective.to_bits(),
            "{what}: objective"
        );
        assert_eq!(
            got.phase1_iterations, want.phase1_iterations,
            "{what}: phase1"
        );
        assert_eq!(
            got.phase2_iterations, want.phase2_iterations,
            "{what}: phase2"
        );
        assert_eq!(got.x.len(), want.x.len(), "{what}: x len");
        for (i, (a, b)) in got.x.iter().zip(&want.x).enumerate() {
            assert_eq!(a.to_bits(), b.to_bits(), "{what}: x[{i}]");
        }
    }

    /// The mixed-shape LP for the sparse-entry tests: every `to_standard`
    /// column kind (plain, shifted-lb, flip-shift MI, fully-free split),
    /// L/G/E senses, negative rhs (row negation), maximize.
    ///
    /// NOTE, resolved by the `farkas_ok` fix: x3's ub
    /// was temporarily 8 because with 5 the LP is infeasible (row 2 forces
    /// x3 ≥ 8) and the phase-1-exit Farkas export debug-asserted on that
    /// shape — its per-column sup `max(φ·lb, φ·ub)` clamped infinite endpoints
    /// to 0, overestimating the sup and rejecting the (mathematically valid)
    /// certificate. `farkas_ok` now maximizes at the endpoint the sign of φ
    /// selects, so the infeasible-and-unbounded-bounds shape verifies and the
    /// branch's original ub = 5 is restored: these LPs now also cover the
    /// Farkas export on a free/MI-variable infeasible LP.
    fn sparse_entry_lps() -> Vec<Lp> {
        let mut a = lp_of(
            &[
                &[1.5, -2.0, 0.0, 4.0],
                &[0.0, 3.0, -1.0, 0.0],
                &[7.0, 0.0, 2.0, -5.0],
            ],
            vec![6.0, -2.0, 9.0],
            "GLE",
            vec![1.0, -2.0, 0.5, 3.0],
            vec![10.0, f64::INFINITY, 5.0, f64::INFINITY],
        );
        a.lower = vec![0.0, 2.0, f64::NEG_INFINITY, f64::NEG_INFINITY];
        a.maximize = true;
        let b = lp_of(
            &[&[1.0, 1.0], &[1.0, 2.0], &[2.0, 1.0]],
            vec![4.0, 6.0, 6.0],
            "GGG",
            vec![1.0, 1.0],
            vec![10.0, 10.0],
        );
        let mut c = lp_of(&[], vec![], "", vec![-1.0], vec![5.0]);
        c.maximize = true;
        vec![a, b, c]
    }

    #[test]
    fn solve_sparse_matches_solve_bit_identical() {
        // Sparse-path spine: with the lever ON (and presolve off, the only mode the
        // sparse path serves), solve_sparse must reproduce solve() bit for
        // bit: same iterations, same objective bits, same x bits.
        for lp in sparse_entry_lps() {
            let want = SimplexSolver::new(SimplexConfig {
                presolve: false,
                ..SimplexConfig::default()
            })
            .solve(&lp)
            .unwrap();
            let sp = SparseLp::from_lp(&lp);
            let got = SimplexSolver::new(SimplexConfig {
                presolve: false,
                sparse: true,
                ..SimplexConfig::default()
            })
            .solve_sparse(&sp)
            .unwrap();
            assert_solutions_bit_identical(&got, &want, "sparse lever on");
        }
    }

    #[test]
    fn solve_sparse_lever_off_and_presolve_on_delegate() {
        // Lever off, and lever on with presolve on (presolve is a
        // dense consumer), must both equal the plain dense solve.
        for presolve in [false, true] {
            for lp in sparse_entry_lps() {
                let want = SimplexSolver::new(SimplexConfig {
                    presolve,
                    ..SimplexConfig::default()
                })
                .solve(&lp)
                .unwrap();
                let sp = SparseLp::from_lp(&lp);
                let got = SimplexSolver::new(SimplexConfig {
                    presolve,
                    sparse: false,
                    ..SimplexConfig::default()
                })
                .solve_sparse(&sp)
                .unwrap();
                assert_solutions_bit_identical(&got, &want, "sparse lever off");
            }
        }
    }

    /// The GGG covering LP from [`sparse_entry_lps`]: every row's slack
    /// start violates (G rows with positive rhs), so the routing detector must
    /// route it to the dual engine.
    fn p30_ggg_lp() -> Lp {
        lp_of(
            &[&[1.0, 1.0], &[1.0, 2.0], &[2.0, 1.0]],
            vec![4.0, 6.0, 6.0],
            "GGG",
            vec![1.0, 1.0],
            vec![10.0, 10.0],
        )
    }

    /// `min x+y` on one L row with positive rhs: the slack-basis start is
    /// feasible (slack basic at 10, no E rows), so the detector must NOT
    /// route and the primal engine keeps the solve.
    fn p30_feasible_start_lp() -> Lp {
        lp_of(
            &[&[1.0, 1.0]],
            vec![10.0],
            "L",
            vec![1.0, 1.0],
            vec![5.0, 5.0],
        )
    }

    #[test]
    fn p1_route_dual_routed_solve_is_the_g17_dual_solve() {
        // Routing spine: with the lever on, a primal-infeasible cold start is
        // solved by the dual engine — bit-identical to the dual lever path
        // (same to_dual_standard + solve_dual_core inputs; the detector only
        // adds the O(m) infeasibility profile in between) — and agrees with
        // the un-routed primal engine on status/objective/x.
        let lp = p30_ggg_lp();
        let routed = SimplexSolver::new(SimplexConfig {
            presolve: false,
            p1_route_dual: true,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        let dual = SimplexSolver::new(SimplexConfig {
            presolve: false,
            dual: true,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_solutions_bit_identical(&routed, &dual, "routed vs dual lever");
        let primal = SimplexSolver::new(SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_eq!(routed.status, SimplexStatus::Optimal);
        assert_eq!(primal.status, SimplexStatus::Optimal);
        assert!((routed.objective - primal.objective).abs() < 1e-9);
        assert_eq!(routed.objective, 4.0);
        for (a, b) in routed.x.iter().zip(&primal.x) {
            assert!((a - b).abs() < 1e-9);
        }
        // The detector route itself: Some(solution), phase-1 iterations 0
        // (no artificial-basis phase 1 was paid).
        let helper = SimplexSolver::new(SimplexConfig {
            presolve: false,
            p1_route_dual: true,
            ..SimplexConfig::default()
        });
        let Some(sol) = helper.route_cold_phase1_dual(&lp).unwrap() else {
            panic!("infeasible start must route");
        };
        assert_eq!(sol.phase1_iterations, 0);
        assert_eq!(sol.status, SimplexStatus::Optimal);
    }

    #[test]
    fn p1_route_dual_feasible_start_stays_primal() {
        // Feasible slack-basis start: the detector returns Ok(None) and the
        // solve is bit-identical to the lever-off run (the primal engine ran;
        // its phase 1 still pivots the artificial out, which pins the path).
        let lp = p30_feasible_start_lp();
        let helper = SimplexSolver::new(SimplexConfig {
            presolve: false,
            p1_route_dual: true,
            ..SimplexConfig::default()
        });
        assert!(helper.route_cold_phase1_dual(&lp).unwrap().is_none());
        let got = helper.solve(&lp).unwrap();
        let want = SimplexSolver::new(SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_solutions_bit_identical(&got, &want, "feasible start, lever on");
        assert!(want.phase1_iterations > 0, "primal phase 1 ran");
    }

    #[test]
    fn p1_route_dual_infeasible_and_presolve_on_agree() {
        // Non-Optimal verdicts must agree too: a routed dual Infeasible on an
        // LP the primal engine also reports Infeasible (G row no box point
        // satisfies), and with presolve ON (default) the routing fires on the
        // REDUCED LP with the same caller-space answer.
        let lp = lp_of(
            &[&[1.0, 2.0], &[2.0, 1.0]],
            vec![6.0, 6.0],
            "GG",
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        );
        let routed = SimplexSolver::new(SimplexConfig {
            presolve: false,
            p1_route_dual: true,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        let primal = SimplexSolver::new(SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_eq!(routed.status, SimplexStatus::Infeasible);
        assert_eq!(primal.status, SimplexStatus::Infeasible);

        let lp = p30_ggg_lp();
        let got = SimplexSolver::new(SimplexConfig {
            p1_route_dual: true,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        let want = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(got.status, SimplexStatus::Optimal);
        assert_eq!(want.status, SimplexStatus::Optimal);
        assert!((got.objective - want.objective).abs() < 1e-9);
    }

    #[test]
    fn csc_matches_dense_dot_bit_exact() {
        // Mixed structure exercising every column kind to_standard emits:
        // plain structural, shifted-lb, flip-shift (lb=-inf, finite ub),
        // fully-free split (±), negated rows (rhs<0), L/G slacks, artificials.
        let mut lp = lp_of(
            &[
                &[1.5, -2.0, 0.0, 4.0],
                &[0.0, 3.0, -1.0, 0.0],
                &[7.0, 0.0, 2.0, -5.0],
            ],
            vec![6.0, -2.0, 9.0],
            "GLE",
            vec![1.0, -2.0, 0.5, 3.0],
            vec![10.0, f64::INFINITY, 5.0, f64::INFINITY],
        );
        lp.lower = vec![0.0, 2.0, f64::NEG_INFINITY, f64::NEG_INFINITY];
        let st = to_standard(&lp).unwrap();
        let mut rng = crate::rng::XorShift64::new(0xDEAD_BEEF_CAFE_1234);
        let x: Vec<f64> = (0..st.m).map(|_| rng.next_f64() * 4.0 - 2.0).collect();
        let mut dense_col = vec![0.0; st.m];
        for v in 0..st.ncol() {
            // CSR well-formedness: strictly ascending rows, no zeros stored.
            let (lo, hi) = (st.csc_ptr[v], st.csc_ptr[v + 1]);
            assert!(hi - lo <= st.m);
            for k in lo + 1..hi {
                assert!(st.csc_row[k - 1] < st.csc_row[k], "row order col {v}");
                assert_ne!(st.csc_val[k], 0.0);
            }
            // Bit-exactness against the dense zip-sum over the gathered
            // column (gather_col is the dense-column oracle now —
            // `to_standard_sparse_matches_dense_bit_identical` independently
            // pins the CSC arrays against the former dense construction).
            st.gather_col(v, 1.0, &mut dense_col);
            let dense: f64 = x.iter().zip(&dense_col).map(|(a, b)| a * b).sum();
            assert_eq!(sparse_dot(&x, &st, v), dense, "col {v}");
            // gather with sgn = −1 is the exact negation on the nnz (the
            // ±1-scale precedent the FTRAN input relies on); structural
            // zeros stay +0.0 by the fill (zero signs normalize).
            let mut neg = vec![0.0; st.m];
            st.gather_col(v, -1.0, &mut neg);
            for k in 0..st.m {
                if dense_col[k] == 0.0 {
                    assert_eq!(neg[k], 0.0, "col {v} row {k}: zero sign");
                } else {
                    assert_eq!(neg[k], -dense_col[k], "col {v} row {k}");
                }
            }
        }
    }

    /// Field-for-field bit-comparison of two `Standard`s (sparse construction
    /// parity): every f64 compared via to_bits, indices/flags by value.
    fn assert_standards_bit_identical(a: &Standard, b: &Standard, what: &str) {
        assert_eq!(a.m, b.m, "{what}: m");
        assert_eq!(a.csc_ptr, b.csc_ptr, "{what}: csc_ptr");
        assert_eq!(a.csc_row, b.csc_row, "{what}: csc_row");
        assert_eq!(a.csc_val.len(), b.csc_val.len(), "{what}: csc nnz");
        for (i, (x, y)) in a.csc_val.iter().zip(&b.csc_val).enumerate() {
            assert_eq!(x.to_bits(), y.to_bits(), "{what}: csc_val[{i}]");
        }
        for (i, (x, y)) in a.b.iter().zip(&b.b).enumerate() {
            assert_eq!(x.to_bits(), y.to_bits(), "{what}: b[{i}]");
        }
        assert_eq!(a.b.len(), b.b.len(), "{what}: b len");
        for (i, (x, y)) in a.ub.iter().zip(&b.ub).enumerate() {
            assert_eq!(x.to_bits(), y.to_bits(), "{what}: ub[{i}]");
        }
        assert_eq!(a.ub.len(), b.ub.len(), "{what}: ub len");
        assert_eq!(a.kind, b.kind, "{what}: kind");
        assert_eq!(a.num_structural, b.num_structural, "{what}: num_structural");
        for (i, (x, y)) in a.cost.iter().zip(&b.cost).enumerate() {
            assert_eq!(x.to_bits(), y.to_bits(), "{what}: cost[{i}]");
        }
        assert_eq!(a.cost.len(), b.cost.len(), "{what}: cost len");
        assert_eq!(a.cols_of, b.cols_of, "{what}: cols_of");
        for (i, (x, y)) in a.shift.iter().zip(&b.shift).enumerate() {
            assert_eq!(x.to_bits(), y.to_bits(), "{what}: shift[{i}]");
        }
        assert_eq!(a.shift.len(), b.shift.len(), "{what}: shift len");
        assert_eq!(a.flip, b.flip, "{what}: flip");
    }

    #[test]
    fn to_standard_delegates_to_the_csc_construction() {
        // `to_standard`'s dense construction is gone — the row-form `Lp`
        // is CSC and the builder IS `to_standard_sparse`. The old dense-oracle
        // unit test retired with it; the load-bearing bit-identity gate is
        // the Netlib A/B spine against the pre-CSC snapshot (28 instances × 4
        // modes, status/obj/iters/x bits). This pin keeps the delegation
        // exact: if anyone ever re-adds independent logic to `to_standard`,
        // this fires instead of silently forking the two constructions.
        for lp in sparse_entry_lps() {
            let direct = to_standard(&lp).unwrap();
            let via_sp = to_standard_sparse(&SparseLp::from_lp(&lp)).unwrap();
            assert_standards_bit_identical(&via_sp, &direct, "delegated vs direct");
        }
    }

    #[test]
    fn to_standard_rejects_malformed_row_form() {
        // The row-form front door keeps its Unsupported contract (dimensions,
        // NaN bounds, inverted bounds, unknown sense) — `SparseLp::from_lp`
        // would panic on all four, so the validation must run first.
        let ok = lp_of(
            &[&[1.0, 1.0], &[1.0, -1.0]],
            vec![4.0, 1.0],
            "LG",
            vec![1.0, 1.0],
            vec![5.0, 5.0],
        );
        assert!(to_standard(&ok).is_ok());
        let bad_sense = {
            let mut lp = ok.clone();
            lp.senses[0] = 'N';
            lp
        };
        match to_standard(&bad_sense) {
            Err(SimplexError::Unsupported(m)) => assert!(m.contains("sense"), "{m}"),
            Err(e) => panic!("unexpected error: {e}"),
            Ok(_) => panic!("expected Unsupported sense"),
        }
        let nan_bound = {
            let mut lp = ok.clone();
            lp.lower[0] = f64::NAN;
            lp
        };
        match to_standard(&nan_bound) {
            Err(SimplexError::Unsupported(m)) => assert!(m.contains("NaN"), "{m}"),
            Err(e) => panic!("unexpected error: {e}"),
            Ok(_) => panic!("expected Unsupported NaN"),
        }
        let inverted = {
            let mut lp = ok.clone();
            lp.upper[0] = -1.0;
            lp
        };
        match to_standard(&inverted) {
            Err(SimplexError::Unsupported(m)) => assert!(m.contains("ub < lb"), "{m}"),
            Err(e) => panic!("unexpected error: {e}"),
            Ok(_) => panic!("expected Unsupported ub < lb"),
        }
        let bad_dims = {
            let mut lp = ok.clone();
            lp.rhs.pop();
            lp
        };
        match to_standard(&bad_dims) {
            Err(SimplexError::Unsupported(m)) => assert!(m.contains("dimensions"), "{m}"),
            Err(e) => panic!("unexpected error: {e}"),
            Ok(_) => panic!("expected Unsupported dims"),
        }
    }

    #[test]
    fn to_standard_sparse_rejects_empty_ub() {
        // The ub = -inf rejection mirrors to_standard's (SparseLp::new accepts
        // the shape, so the solver-side check must survive at this layer).
        let mut lp = lp_of(
            &[&[1.0]],
            vec![1.0],
            "L",
            vec![1.0],
            vec![f64::NEG_INFINITY],
        );
        lp.lower = vec![f64::NEG_INFINITY];
        let sp = SparseLp::from_lp(&lp);
        match to_standard_sparse(&sp) {
            Err(SimplexError::Unsupported(m)) => assert!(m.contains("ub = -inf"), "{m}"),
            Err(e) => panic!("unexpected error: {e}"),
            Ok(_) => panic!("expected Unsupported for ub = -inf"),
        }
    }

    #[test]
    fn simple_g_min() {
        // min 3x+2y s.t. x+y>=5 -> 10 at (0,5)
        let s = solve(
            &[&[1.0, 1.0]],
            vec![5.0],
            "G",
            vec![3.0, 2.0],
            vec![f64::INFINITY; 2],
        );
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 10.0).abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn negated_rhs_row() {
        // min -x s.t. -x >= -3 (i.e. x <= 3) -> -3; exercises row negation
        let s = solve(
            &[&[-1.0]],
            vec![-3.0],
            "G",
            vec![-1.0],
            vec![f64::INFINITY; 1],
        );
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective + 3.0).abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn equality_row() {
        // min x+y s.t. x+y = 7, x <= 5 -> 7
        let s = solve(
            &[&[1.0, 1.0]],
            vec![7.0],
            "E",
            vec![1.0, 1.0],
            vec![5.0, f64::INFINITY],
        );
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 7.0).abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn infeasible_detected() {
        // x <= 1, x >= 3
        let s = solve(
            &[&[1.0], &[1.0]],
            vec![1.0, 3.0],
            "LG",
            vec![1.0],
            vec![f64::INFINITY; 1],
        );
        assert_eq!(s.status, SimplexStatus::Infeasible);
    }

    #[test]
    fn unbounded_detected() {
        // min -x, no rows
        let s = solve(&[], vec![], "", vec![-1.0], vec![f64::INFINITY; 1]);
        assert_eq!(s.status, SimplexStatus::Unbounded);
    }

    #[test]
    fn ub_saturation_and_row_limit() {
        // min -2x - y s.t. x+y<=12, ub x<=4, y<=10 -> x=4, y=8, obj -16
        let s = solve(
            &[&[1.0, 1.0]],
            vec![12.0],
            "L",
            vec![-2.0, -1.0],
            vec![4.0, 10.0],
        );
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective + 16.0).abs() < 1e-6, "obj {}", s.objective);
        assert!((s.x[0] - 4.0).abs() < 1e-6, "x {:?}", s.x);
        assert!((s.x[1] - 8.0).abs() < 1e-6, "x {:?}", s.x);
    }

    #[test]
    fn pure_ub_flip_no_rows() {
        // min -x s.t. ub x<=5 (no constraint rows) -> -5 via ub flip
        let s = solve(&[], vec![], "", vec![-1.0], vec![5.0]);
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective + 5.0).abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn lower_bound_shift() {
        // min x s.t. x >= 5, lb x>=2 -> 5 at x=5 (in original space)
        let mut lp = lp_of(&[&[1.0]], vec![5.0], "G", vec![1.0], vec![f64::INFINITY]);
        lp.lower = vec![2.0];
        let s = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 5.0).abs() < 1e-6, "obj {}", s.objective);
        assert!((s.x[0] - 5.0).abs() < 1e-6, "x {:?}", s.x);
    }

    #[test]
    fn maximize_via_sign() {
        // max x+y s.t. x+y<=10 -> 10
        let mut lp = lp_of(
            &[&[1.0, 1.0]],
            vec![10.0],
            "L",
            vec![1.0, 1.0],
            vec![f64::INFINITY; 2],
        );
        lp.maximize = true;
        let s = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 10.0).abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn degenerate_ties_terminate() {
        // min x+y s.t. x+y>=4, x+2y>=6, 2x+y>=6 -> 4 at (2,2), triple-tight
        let s = solve(
            &[&[1.0, 1.0], &[1.0, 2.0], &[2.0, 1.0]],
            vec![4.0, 6.0, 6.0],
            "GGG",
            vec![1.0, 1.0],
            vec![10.0, 10.0],
        );
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 4.0).abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn steepest_edge_matches_dantzig() {
        // Same degenerate LP + a bounded one, solved with sh_pricing ON:
        // the (provisional) SH rate must reach the identical optima the
        // Dantzig path finds — the mechanism, not the rate quality, is what
        // this gates. Bench measures iteration-count benefit separately.
        let cfg = SimplexConfig {
            sh_pricing: true,
            ..SimplexConfig::default()
        };
        let solver = SimplexSolver::new(cfg);
        let s1 = solver
            .solve(&lp_of(
                &[&[1.0, 1.0], &[1.0, 2.0], &[2.0, 1.0]],
                vec![4.0, 6.0, 6.0],
                "GGG",
                vec![1.0, 1.0],
                vec![10.0, 10.0],
            ))
            .unwrap();
        assert_eq!(s1.status, SimplexStatus::Optimal);
        assert!((s1.objective - 4.0).abs() < 1e-6, "SH obj {}", s1.objective);

        // min x+y s.t. x+y>=5 with ub split (bounded, exercises from-ub):
        // maximize side — reuse the max pattern: max x+y s.t. x+y<=5 -> 5.
        let mut lp = lp_of(
            &[&[1.0, 1.0]],
            vec![5.0],
            "L",
            vec![1.0, 1.0],
            vec![f64::INFINITY; 2],
        );
        lp.maximize = true;
        let s2 = solver.solve(&lp).unwrap();
        assert_eq!(s2.status, SimplexStatus::Optimal);
        assert!(
            (s2.objective - 5.0).abs() < 1e-6,
            "SH max obj {}",
            s2.objective
        );
    }

    #[test]
    fn zero_artificial_stays_basic() {
        // min x s.t. x - y = 0, x <= 1 -> 0
        let s = solve(
            &[&[1.0, -1.0], &[1.0, 0.0]],
            vec![0.0, 1.0],
            "EL",
            vec![1.0, 0.0],
            vec![f64::INFINITY; 2],
        );
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!(s.objective.abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn free_variable_split() {
        // min x + y s.t. x − y = 2, x ∈ [0,5], y free (FR):
        // y = x−2 ⇒ obj = 2x−2 ⇒ x=0, y=−2, obj=−2
        let mut lp = lp_of(
            &[&[1.0, -1.0]],
            vec![2.0],
            "E",
            vec![1.0, 1.0],
            vec![5.0, f64::INFINITY],
        );
        lp.lower = vec![0.0, f64::NEG_INFINITY];
        let s = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective + 2.0).abs() < 1e-6, "obj {}", s.objective);
        assert!((s.x[1] + 2.0).abs() < 1e-6, "x {:?}", s.x);
    }

    #[test]
    fn mi_bound_flip_shift() {
        // min −x s.t. −x ≤ 0 (i.e. x ≥ 0), x ≤ 5, x MI (−inf..5):
        // optimum x=5, obj −5  (exercises the ub-flip shift x = 5 − x̃)
        let mut lp = lp_of(
            &[&[-1.0, 0.0], &[1.0, 0.0]],
            vec![0.0, 5.0],
            "LL",
            vec![-1.0, 0.0],
            vec![5.0, 1.0],
        );
        lp.lower = vec![f64::NEG_INFINITY, 0.0];
        let s = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective + 5.0).abs() < 1e-6, "obj {}", s.objective);
        assert!((s.x[0] - 5.0).abs() < 1e-6, "x {:?}", s.x);
    }

    #[test]
    fn free_variable_unbounded() {
        // min x with x fully free and no rows ⇒ unbounded below
        let mut lp = lp_of(&[], vec![], "", vec![1.0], vec![f64::INFINITY]);
        lp.lower = vec![f64::NEG_INFINITY];
        for cfg in [
            SimplexConfig::default(),
            SimplexConfig {
                presolve: false,
                ..SimplexConfig::default()
            },
        ] {
            let s = SimplexSolver::new(cfg).solve(&lp).unwrap();
            assert_eq!(s.status, SimplexStatus::Unbounded);
        }
    }

    #[test]
    fn nan_bound_rejected() {
        let mut lp = lp_of(&[], vec![], "", vec![1.0], vec![f64::INFINITY]);
        lp.lower = vec![f64::NAN];
        let r = SimplexSolver::new(SimplexConfig::default()).solve(&lp);
        assert!(matches!(r, Err(SimplexError::Unsupported(_))));
    }

    #[test]
    fn blend_like_structure() {
        // min 2a + 3b s.t. a+b>=100 (demand), 2a+b<=150 (capacity), ub a<=40
        // optimum: a=40,b=60 -> 260; check
        let s = solve(
            &[&[1.0, 1.0], &[2.0, 1.0]],
            vec![100.0, 150.0],
            "GL",
            vec![2.0, 3.0],
            vec![40.0, f64::INFINITY],
        );
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 260.0).abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn p1_tie_eps_is_bounded_and_deterministic() {
        // Phase-1 perturbation invariant: eps[v] ∈ [0, ptol/2) (sub-tolerance ⇒ the
        // unperturbed improver test can never disagree with the perturbed
        // ranking), deterministic per column, and genuinely spread (not a
        // constant, which would degenerate to no perturbation at all).
        let ptol = 1e-7;
        let a = p1_tie_eps(1000, ptol, 0);
        let b = p1_tie_eps(1000, ptol, 0);
        assert_eq!(a, b, "eps must be deterministic per column");
        assert!(
            a.iter()
                .all(|&e| (0.0..ptol * 0.5).contains(&e) && e.is_finite())
        );
        let distinct: std::collections::HashSet<u64> =
            a[..256].iter().map(|e| e.to_bits()).collect();
        assert!(
            distinct.len() > 200,
            "eps must spread, got {}/256",
            distinct.len()
        );
    }

    #[test]
    fn p1_tie_eps_seed_zero_is_historical_and_fresh_seeds_differ() {
        // Seed 0 must reproduce the original one-shot vector bit-exactly
        // (non-re-randomizing runs stay bit-identical), and a re-randomized
        // draw must be a genuinely different vector that keeps the
        // sub-tolerance bound (never reclassifies an improver).
        let ptol = 1e-7;
        let hist: Vec<f64> = (0..1000usize)
            .map(|v| ptol * 0.5 * hash01(v as u64))
            .collect();
        assert_eq!(
            p1_tie_eps(1000, ptol, 0),
            hist,
            "seed 0 = original one-shot bits"
        );
        let seed = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(0x517C_1B72);
        let fresh = p1_tie_eps(1000, ptol, seed);
        assert_ne!(fresh, hist, "a fresh seed must reorder");
        assert!(
            fresh
                .iter()
                .all(|&e| (0.0..ptol * 0.5).contains(&e) && e.is_finite())
        );
        assert_eq!(
            fresh,
            p1_tie_eps(1000, ptol, seed),
            "per-seed deterministic"
        );
    }

    #[test]
    fn p1_eps_never_flips_improver_membership() {
        // Phase-1 perturbation property test: random reduced costs clustered AT the ±ptol
        // improvement boundary. (1) The improver SET is decided on the
        // UNPERTURBED d, so membership cannot flip by construction; what the
        // eps bound must guarantee is that the perturbed key keeps every
        // true improver strictly negative (never demoted out of the race by
        // its own eps) and that argmin(d+eps) over improvers selects an
        // improver. (2) The at-ub side (d > ptol): key stays > ptol/2.
        let ptol = 1e-7;
        let eps = p1_tie_eps(8, ptol, 0);
        let mut rng = crate::rng::XorShift64::new(0x1234_5678_9ABC_DEF0);
        for _ in 0..10_000 {
            // d uniform in ±1.1·ptol: dense on both sides of ±ptol.
            let d: Vec<f64> = (0..8)
                .map(|_| (rng.next_f64() * 2.0 - 1.0) * ptol * 1.1)
                .collect();
            let improvers: Vec<usize> = (0..8).filter(|&v| d[v] < -ptol).collect();
            let at_ub_improvers: Vec<usize> = (0..8).filter(|&v| d[v] > ptol).collect();
            for &v in &improvers {
                assert!(
                    d[v] + eps[v] < -ptol / 2.0,
                    "improver v={v} d={} eps={} key not < -ptol/2",
                    d[v],
                    eps[v]
                );
            }
            for &v in &at_ub_improvers {
                assert!(d[v] + eps[v] > ptol / 2.0, "at-ub improver flipped");
            }
            if !improvers.is_empty() {
                let q = improvers
                    .iter()
                    .copied()
                    .min_by(|x, y| (d[*x] + eps[*x]).total_cmp(&(d[*y] + eps[*y])))
                    .unwrap();
                assert!(d[q] < -ptol, "argmin(d+eps) picked a non-improver");
            }
        }
    }

    #[test]
    fn p1_escalation_preserves_optimum() {
        // Phase-1 perturbation end-to-end wiring: bland_after=1 drops the escalation rung to
        // 4 stalls, so any phase-1 degeneracy trips the perturbation. The
        // optima must be untouched (entering SET + phase-1 exit stay
        // unperturbed). Degenerate trio + blend-like LP, presolve ON and OFF.
        let lps = [
            lp_of(
                &[&[1.0, 1.0], &[1.0, 2.0], &[2.0, 1.0]],
                vec![4.0, 6.0, 6.0],
                "GGG",
                vec![1.0, 1.0],
                vec![10.0, 10.0],
            ),
            lp_of(
                &[&[1.0, 1.0], &[2.0, 1.0]],
                vec![100.0, 150.0],
                "GL",
                vec![2.0, 3.0],
                vec![40.0, f64::INFINITY],
            ),
        ];
        let want = [4.0, 260.0];
        for presolve in [true, false] {
            let cfg = SimplexConfig {
                bland_after: 1,
                presolve,
                ..SimplexConfig::default()
            };
            let solver = SimplexSolver::new(cfg);
            for (lp, &w) in lps.iter().zip(&want) {
                let s = solver.solve(lp).unwrap();
                assert_eq!(s.status, SimplexStatus::Optimal);
                assert!(
                    (s.objective - w).abs() < 1e-6,
                    "presolve={presolve} obj {} want {w}",
                    s.objective
                );
            }
        }
    }

    #[test]
    fn p2_eps_never_flips_improver_membership() {
        // Phase-2 perturbation property test (mirrors p1_eps_never_flips_improver_membership):
        // random reduced costs clustered AT the ±ptol boundary. Membership is
        // decided on the UNPERTURBED d (|d| > ptol), so what eps < ptol/2 must
        // guarantee is (1) every true improver's perturbed cost keeps its sign
        // beyond ptol/2 — never demoted across zero, so both ranking keys
        // (Dantzig |d+eps|, SH (d+eps)²/w) stay well-defined and monotone on
        // the improver set; and (2) argmax of either perturbed key over the
        // TRUE improver set picks a true improver.
        let ptol = 1e-7;
        let eps = p1_tie_eps(8, ptol, 0);
        let w: Vec<f64> = (0..8).map(|v| 1.0 + 50.0 * hash01(v + 77)).collect();
        let mut rng = crate::rng::XorShift64::new(0xDEAD_BEEF_CAFE_1234);
        for _ in 0..10_000 {
            // d uniform in ±1.1·ptol: dense on both sides of ±ptol.
            let d: Vec<f64> = (0..8)
                .map(|_| (rng.next_f64() * 2.0 - 1.0) * ptol * 1.1)
                .collect();
            let improvers: Vec<usize> = (0..8).filter(|&v| d[v] < -ptol || d[v] > ptol).collect();
            for &v in &improvers {
                let k = d[v] + eps[v];
                assert!(
                    (d[v] < -ptol && k < -ptol / 2.0) || (d[v] > ptol && k > ptol / 2.0),
                    "improver v={v} d={} eps={} perturbed cost crossed zero",
                    d[v],
                    eps[v]
                );
            }
            if !improvers.is_empty() {
                // Dantzig path: argmax |d+eps|.
                let q = improvers
                    .iter()
                    .copied()
                    .max_by(|x, y| (d[*x] + eps[*x]).abs().total_cmp(&(d[*y] + eps[*y]).abs()))
                    .unwrap();
                assert!(d[q].abs() > ptol, "argmax|d+eps| picked a non-improver");
                // SH path: argmax (d+eps)²/w.
                let q = improvers
                    .iter()
                    .copied()
                    .max_by(|x, y| {
                        let kx = (d[*x] + eps[*x]).powi(2) / w[*x];
                        let ky = (d[*y] + eps[*y]).powi(2) / w[*y];
                        kx.total_cmp(&ky)
                    })
                    .unwrap();
                assert!(d[q].abs() > ptol, "argmax (d+eps)²/w picked a non-improver");
            }
        }
    }

    /// File-based regression for the optimal-verdict defense: node LPs
    /// of `jobshop_bigM[wideM]` seed 2 where the cold two-phase primal
    /// returned Optimal at the RIGHT objective but an x violating a big-M G
    /// row by 42 / a binary bound by 1e-4. Contract asserted: the returned
    /// x satisfies the caller's rows and bounds at the `feas_check.py` bar
    /// (`1e-6 · max(1, |x|∞)`), at the HiGHS-confirmed objective.
    #[test]
    fn c4_optimal_x_is_feasible_on_big_m_node_lps() {
        use rusolve_io::{ObjSense, RowSense};
        for (name, obj_known) in [
            ("jobshop_wideM_s2_node_rowviol", 287.00331502652784),
            ("jobshop_wideM_s2_node_overshoot", 287.0014429437328),
        ] {
            let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../benches/data/repros/")
                .to_string()
                + name
                + ".mps";
            let text = std::fs::read_to_string(&path).expect("repro file present");
            let m = rusolve_io::read_mps(&text).expect("repro parses");
            assert_eq!(m.ranges.iter().filter(|&&r| r != 0.0).count(), 0);
            assert!(!m.bounds.iter().any(|b| b.integer));
            let senses: Vec<char> = m
                .senses
                .iter()
                .map(|s| match s {
                    RowSense::LessEqual => 'L',
                    RowSense::GreaterEqual => 'G',
                    RowSense::Equal => 'E',
                })
                .collect();
            let lp = Lp {
                num_col: m.num_variables,
                num_row: m.num_constraints,
                a: m.matrix.clone(),
                rhs: m.rhs.clone(),
                objective: m.objective.clone(),
                maximize: m.obj_sense == ObjSense::Maximize,
                senses,
                lower: m.bounds.iter().map(|b| b.lower).collect(),
                upper: m.bounds.iter().map(|b| b.upper).collect(),
            };
            let s = SimplexSolver::new(SimplexConfig {
                presolve: false,
                ..SimplexConfig::default()
            })
            .solve(&lp)
            .unwrap();
            assert_eq!(s.status, SimplexStatus::Optimal, "{name}");
            assert!(
                (s.objective - obj_known).abs() <= 1e-6 * obj_known.abs(),
                "{name} obj {}",
                s.objective
            );
            let scale = s.x.iter().fold(1.0f64, |a, &v| a.max(v.abs()));
            let (wb, wr) = primal_violations(&lp, &s.x);
            assert!(
                wb <= 1e-6 * scale,
                "{name}: bound violation {wb:.3e} > tol {}",
                1e-6 * scale
            );
            assert!(
                wr <= 1e-6 * scale,
                "{name}: row violation {wr:.3e} > tol {}",
                1e-6 * scale
            );
        }
    }

    /// File-based regression: Netlib `brandy` with presolve OFF,
    /// primal two-phase, hit the 100k iteration cap in a phase-1 fp cycle —
    /// w_true bit-frozen at 958.03 for ~99k θ=0 pivots at ONE vertex under
    /// sticky Bland + the fixed (one-shot) phase-1 eps vector. With HiGHS-style
    /// re-perturbation (SimplexConfig::p1_re_perturb, default ON) the run
    /// re-randomizes eps on the frozen-w probe and reaches the optimum
    /// (HiGHS primal strategy=4: Optimal 1518.5099).
    #[test]
    fn p19_brandy_no_presolve_primal_reaches_optimum() {
        use rusolve_io::{ObjSense, RowSense};
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../benches/data/netlib/").to_string()
            + "brandy.mps";
        let text = std::fs::read_to_string(&path).expect("brandy present");
        let m = rusolve_io::read_mps(&text).expect("brandy parses");
        let senses: Vec<char> = m
            .senses
            .iter()
            .map(|s| match s {
                RowSense::LessEqual => 'L',
                RowSense::GreaterEqual => 'G',
                RowSense::Equal => 'E',
            })
            .collect();
        let lp = Lp {
            num_col: m.num_variables,
            num_row: m.num_constraints,
            a: m.matrix.clone(),
            rhs: m.rhs.clone(),
            objective: m.objective.clone(),
            maximize: m.obj_sense == ObjSense::Maximize,
            senses,
            lower: m.bounds.iter().map(|b| b.lower).collect(),
            upper: m.bounds.iter().map(|b| b.upper).collect(),
        };
        let s = SimplexSolver::new(SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        let obj_known = 1518.5098957;
        assert!(
            (s.objective - obj_known).abs() <= 1e-6 * obj_known.abs(),
            "obj {}",
            s.objective
        );
        let scale = s.x.iter().fold(1.0f64, |a, &v| a.max(v.abs()));
        let (wb, wr) = primal_violations(&lp, &s.x);
        assert!(wb <= 1e-6 * scale && wr <= 1e-6 * scale, "viol {wb}/{wr}");
    }

    /// Regression: primal with presolve OFF on Netlib `tuff` returned a silent
    /// wrong "Optimal 0" (phase 1 ended early on a drifted eta chain). With
    /// incremental reduced costs there is one BTRAN less per iteration, the
    /// synthetic refactor clock ran slower (352-eta chains) and the chain lost
    /// accuracy unnoticed. The pivot column/row cross-check (HiGHS
    /// updatePivots style) now forces a refactorization. HiGHS: 0.2921477651.
    #[test]
    fn tuff_primal_no_presolve_reaches_the_true_optimum() {
        use rusolve_io::{ObjSense, RowSense};
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../benches/data/netlib/tuff.mps"
        );
        let Ok(text) = std::fs::read_to_string(path) else {
            eprintln!("tuff.mps not present; skipping");
            return;
        };
        let m = rusolve_io::read_mps(&text).expect("tuff parses");
        let senses: Vec<char> = m
            .senses
            .iter()
            .map(|s| match s {
                RowSense::LessEqual => 'L',
                RowSense::GreaterEqual => 'G',
                RowSense::Equal => 'E',
            })
            .collect();
        let lp = Lp {
            num_col: m.num_variables,
            num_row: m.num_constraints,
            a: m.matrix.clone(),
            rhs: m.rhs.clone(),
            objective: m.objective.clone(),
            maximize: m.obj_sense == ObjSense::Maximize,
            senses,
            lower: m.bounds.iter().map(|b| b.lower).collect(),
            upper: m.bounds.iter().map(|b| b.upper).collect(),
        };
        let s = SimplexSolver::new(SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!(
            (s.objective - 0.2921477651).abs() <= 1e-6,
            "tuff objective {}",
            s.objective
        );
        assert!(primal_feasible(&lp, &s.x));
    }

    #[test]
    fn a2_bounded_episodes_terminate_beale_cycle() {
        // Bounded-episode semantics on a degenerate cycle-prone LP:
        // Beale's example, no perturbation, no SH ⇒ the only escape
        // machinery is the Bland episode. With tiny budgets the episode must
        // release, re-enter at a doubled stall threshold, and still land the
        // optimum — bounded coverage, not abandonment of anti-cycling.
        let lp = lp_of(
            &[
                &[0.25, -60.0, -0.04, 9.0],
                &[0.5, -90.0, -0.02, 3.0],
                &[1.0, 0.0, 1.0, 0.0],
            ],
            vec![0.0, 0.0, 1.0],
            "LLL",
            vec![-0.75, 150.0, -0.02, 6.0],
            vec![f64::INFINITY; 4],
        );
        for episode in [8usize, 32, 200] {
            let s = SimplexSolver::new(SimplexConfig {
                bland_after: 1,
                bland_episode: episode,
                p2_tie_perturb: false,
                sh_pricing: false,
                ..SimplexConfig::default()
            })
            .solve(&lp)
            .unwrap();
            assert_eq!(s.status, SimplexStatus::Optimal, "episode={episode}");
            assert!(
                (s.objective - (-0.0480769230769231)).abs() < 1e-6,
                "episode={episode} obj {}",
                s.objective
            );
            assert!(s.bland_triggered, "episode={episode}: never engaged");
        }
    }

    #[test]
    fn p2_escalation_preserves_optimum() {
        // Phase-2 perturbation end-to-end wiring: bland_after=1 drops the escalation rung to
        // 4 stalls, so phase-2 degeneracy trips the perturbation. Optima must
        // be untouched (improver SET, ratio test, leaving choice, objective
        // bookkeeping all stay unperturbed). Degenerate trio + blend + Beale
        // cycle + tied-cost crawl LPs, perturbation on/off × SH on/off ×
        // presolve on/off.
        let lps = [
            (
                lp_of(
                    &[&[1.0, 1.0], &[1.0, 2.0], &[2.0, 1.0]],
                    vec![4.0, 6.0, 6.0],
                    "GGG",
                    vec![1.0, 1.0],
                    vec![10.0, 10.0],
                ),
                4.0,
            ),
            (
                lp_of(
                    &[&[1.0, 1.0], &[2.0, 1.0]],
                    vec![100.0, 150.0],
                    "GL",
                    vec![2.0, 3.0],
                    vec![40.0, f64::INFINITY],
                ),
                260.0,
            ),
            (
                lp_of(
                    &[
                        &[0.25, -60.0, -0.04, 9.0],
                        &[0.5, -90.0, -0.02, 3.0],
                        &[1.0, 0.0, 1.0, 0.0],
                    ],
                    vec![0.0, 0.0, 1.0],
                    "LLL",
                    vec![-0.75, 150.0, -0.02, 6.0],
                    vec![f64::INFINITY; 4],
                ),
                -0.0480769230769231, // Beale's cycling example
            ),
            (
                lp_of(
                    &[
                        &[1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                        &[2.0, 1.0, 1.0, 1.0, 1.0, 0.0],
                        &[1.0, 2.0, 1.0, 1.0, 0.0, 1.0],
                        &[1.0, 1.0, 2.0, 0.0, 1.0, 1.0],
                        &[1.0, 1.0, 1.0, 3.0, 1.0, 1.0],
                    ],
                    vec![4.0, 4.0, 4.0, 4.0, 5.0],
                    "LLLLL",
                    vec![-1.0, -1.0001, -0.9999, -1.0, -1.0, -1.0],
                    vec![1.0; 6],
                ),
                -4.0001,
            ),
        ];
        for presolve in [true, false] {
            for sh in [true, false] {
                for pert in [true, false] {
                    let cfg = SimplexConfig {
                        bland_after: 1,
                        presolve,
                        sh_pricing: sh,
                        p2_tie_perturb: pert,
                        ..SimplexConfig::default()
                    };
                    let solver = SimplexSolver::new(cfg);
                    for (i, (lp, want)) in lps.iter().enumerate() {
                        let s = solver.solve(lp).unwrap();
                        assert_eq!(
                            s.status,
                            SimplexStatus::Optimal,
                            "lp={i} presolve={presolve} sh={sh} pert={pert}"
                        );
                        assert!(
                            (s.objective - want).abs() < 1e-6,
                            "lp={i} presolve={presolve} sh={sh} pert={pert} obj {} want {want}",
                            s.objective
                        );
                    }
                }
            }
        }
    }

    // ---- Row duals, reduced costs, Farkas ray, unbounded ray ----

    fn solve_np(lp: &Lp) -> SimplexSolution {
        SimplexSolver::new(SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        })
        .solve(lp)
        .unwrap()
    }

    /// Property check on the exported `(π, d)` pair: (1) the defining
    /// identity `d_j = c_j − πᵀa_j` holds for ALL j in the caller's space;
    /// (2) caller-sense sign conditions — min: d ≥ 0 at lb / ≤ 0 at ub,
    /// π ≤ 0 on 'L' / ≥ 0 on 'G'; max: everything flips; (3) complementary
    /// slackness — columns strictly between bounds have d = 0, rows strictly
    /// slack have π = 0.
    fn check_duals(lp: &Lp, s: &SimplexSolution) {
        assert!(s.duals_valid, "presolve-off Optimal must export duals");
        assert_eq!(s.row_duals.len(), lp.num_row);
        assert_eq!(s.reduced_costs.len(), lp.num_col);
        let sign = if lp.maximize { -1.0 } else { 1.0 }; // min-normaliser
        let cols = dense_cols(lp);
        for j in 0..lp.num_col {
            let dot: f64 = (0..lp.num_row).map(|i| s.row_duals[i] * cols[j][i]).sum();
            let want = lp.objective[j] - dot;
            assert!(
                (s.reduced_costs[j] - want).abs() <= 1e-9 * (1.0 + want.abs()),
                "identity d_{j} = c_j - pi'a_j: got {}, want {want}",
                s.reduced_costs[j]
            );
        }
        for i in 0..lp.num_row {
            let p = s.row_duals[i] * sign; // min-normalised
            let ctx = format!(
                "row {i} sense {} pi {} maximize {}",
                lp.senses[i], s.row_duals[i], lp.maximize
            );
            match lp.senses[i] {
                'L' => assert!(p <= 1e-7, "{ctx}"),
                'G' => assert!(p >= -1e-7, "{ctx}"),
                _ => {}
            }
            let lhs: f64 = (0..lp.num_col).map(|j| cols[j][i] * s.x[j]).sum();
            let slack = match lp.senses[i] {
                'L' => lp.rhs[i] - lhs,
                'G' => lhs - lp.rhs[i],
                _ => 0.0,
            };
            if slack > 1e-6 {
                assert!(
                    s.row_duals[i].abs() <= 1e-6 * (1.0 + lhs.abs()),
                    "row {i} slack {slack} but pi {}",
                    s.row_duals[i]
                );
            }
        }
        for j in 0..lp.num_col {
            if lp.lower[j] + 1e-7 >= lp.upper[j] {
                continue; // fixed variable: never priced, no sign condition
            }
            let d = s.reduced_costs[j] * sign; // min-normalised
            if s.x[j] <= lp.lower[j] + 1e-7 {
                assert!(d >= -1e-6, "col {j} at lb: d*sign {d}");
            } else if s.x[j] >= lp.upper[j] - 1e-7 {
                assert!(d <= 1e-6, "col {j} at ub: d*sign {d}");
            } else {
                assert!(d.abs() <= 1e-6, "col {j} interior: d*sign {d}");
            }
        }
    }

    #[test]
    fn dual_identity_and_sign_conditions_hold_for_all_columns() {
        // Mixed-sense LP (L/G/E, at-ub column, positive lb, E row), its
        // maximization twin, a negated-row LP (rhs < 0 ⇒ row flip), a
        // ub-saturation LP, the MI flip-shift LP and a free-variable LP —
        // all presolve OFF, all through the full property check.
        let lps = [
            lp_of(
                &[&[1.0, 1.0, 1.0], &[1.0, -1.0, 0.0], &[0.0, 1.0, 1.0]],
                vec![10.0, 1.0, 2.0],
                "LEG",
                vec![-1.0, -2.0, 1.0],
                vec![3.0, f64::INFINITY, f64::INFINITY],
            ),
            lp_of(&[&[-1.0]], vec![-3.0], "G", vec![-1.0], vec![f64::INFINITY]),
            lp_of(
                &[&[1.0, 1.0]],
                vec![12.0],
                "L",
                vec![-2.0, -1.0],
                vec![4.0, 10.0],
            ),
            lp_of(
                &[&[-1.0, 0.0], &[1.0, 0.0]],
                vec![0.0, 5.0],
                "LL",
                vec![-1.0, 0.0],
                vec![5.0, 1.0],
            ),
            lp_of(
                &[&[1.0, -1.0]],
                vec![2.0],
                "E",
                vec![1.0, 1.0],
                vec![5.0, f64::INFINITY],
            ),
        ];
        for mut lp in lps {
            let s = solve_np(&lp);
            assert_eq!(s.status, SimplexStatus::Optimal, "min {:?}", lp.senses);
            check_duals(&lp, &s);
            lp.maximize = true;
            let s = solve_np(&lp);
            assert_eq!(s.status, SimplexStatus::Optimal, "max {:?}", lp.senses);
            check_duals(&lp, &s);
        }
        // The mixed LP needs its positive lb set explicitly (lp_of gives 0).
        let mut lp = lp_of(
            &[&[1.0, 1.0, 1.0], &[1.0, -1.0, 0.0], &[0.0, 1.0, 1.0]],
            vec![10.0, 1.0, 2.0],
            "LEG",
            vec![-1.0, -2.0, 1.0],
            vec![3.0, f64::INFINITY, f64::INFINITY],
        );
        lp.lower = vec![0.0, 1.0, 0.0];
        let s = solve_np(&lp);
        assert_eq!(s.status, SimplexStatus::Optimal);
        check_duals(&lp, &s);
    }

    #[test]
    fn duals_spot_check_known_values_and_complementary_slackness() {
        // max 3x + 2y s.t. x + y ≤ 4 (the doc example): π = +3, d = (0, −1).
        let mut lp = lp_of(
            &[&[1.0, 1.0]],
            vec![4.0],
            "L",
            vec![3.0, 2.0],
            vec![f64::INFINITY; 2],
        );
        lp.maximize = true;
        let s = solve_np(&lp);
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.row_duals[0] - 3.0).abs() < 1e-9, "pi {:?}", s.row_duals);
        assert!(
            (s.reduced_costs[0].abs() < 1e-9 && (s.reduced_costs[1] + 1.0).abs() < 1e-9),
            "d {:?}",
            s.reduced_costs
        );
        // blend: min 2a+3b, a+b≥100 (G, tight), 2a+b≤150 (L, slack 10),
        // a ≤ 40 ⇒ (40,60): π = (3, 0) — the slack row's dual is EXACTLY 0
        // (row complementary slackness) — b interior ⇒ d_b = 0, a at ub ⇒
        // d_a = 2 − 3 = −1 ≤ 0.
        let s = solve_np(&lp_of(
            &[&[1.0, 1.0], &[2.0, 1.0]],
            vec![100.0, 150.0],
            "GL",
            vec![2.0, 3.0],
            vec![40.0, f64::INFINITY],
        ));
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!(
            (s.row_duals[0] - 3.0).abs() < 1e-9 && s.row_duals[1].abs() < 1e-9,
            "pi {:?}",
            s.row_duals
        );
        assert!(
            (s.reduced_costs[0] + 1.0).abs() < 1e-9 && s.reduced_costs[1].abs() < 1e-9,
            "d {:?}",
            s.reduced_costs
        );
    }

    /// Property check, recomputed from the ORIGINAL rows: row sign
    /// conditions, column finiteness conditions, and the strict gap
    /// πᵀrhs > Σ_j max{φ·x_j : x_j ∈ [lb, ub]} with φ = π_fᵀa_j. The per-column
    /// max is taken at the endpoint the SIGN of φ selects (clamping the
    /// infinite endpoint to 0 and maximizing over the finite ones instead
    /// overestimates the sup and rejects valid certificates); only the
    /// `|φ| ≤ tol` noise band keeps the endpoint clamp.
    fn check_farkas(lp: &Lp, y: &[f64]) {
        assert_eq!(y.len(), lp.num_row);
        for (i, &s) in lp.senses.iter().enumerate() {
            match s {
                'L' => assert!(y[i] <= 1e-6, "L row y {} > 0", y[i]),
                'G' => assert!(y[i] >= -1e-6, "G row y {} < 0", y[i]),
                _ => {}
            }
        }
        let y_inf = y.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
        let tol = 1e-6 * (1.0 + y_inf);
        let cols = dense_cols(lp);
        let mut sup = 0.0f64;
        for j in 0..lp.num_col {
            let phi: f64 = (0..lp.num_row).map(|i| y[i] * cols[j][i]).sum();
            let (lb, ub) = (lp.lower[j], lp.upper[j]);
            if phi > tol {
                assert!(
                    ub < f64::INFINITY,
                    "col {j} phi {phi} > 0 with ub = +inf: sup is +inf"
                );
                sup += phi * ub;
            } else if phi < -tol {
                assert!(
                    lb > f64::NEG_INFINITY,
                    "col {j} phi {phi} < 0 with lb = -inf: sup is +inf"
                );
                sup += phi * lb;
            } else {
                let a = if lb.is_finite() { phi * lb } else { 0.0 };
                let b = if ub.is_finite() { phi * ub } else { 0.0 };
                sup += a.max(b);
            }
        }
        let lhs: f64 = y.iter().zip(&lp.rhs).map(|(p, r)| p * r).sum();
        assert!(
            lhs > sup + 1e-9 * lhs.abs().max(sup.abs()).max(1.0),
            "pi'rhs {lhs} must exceed box sup {sup}"
        );
    }

    #[test]
    fn farkas_certificate_certifies_infeasibility() {
        // x ≤ 1 ∧ x ≥ 3 (L/G pair) and x = 1 ∧ x = 2 (E/E pair), presolve
        // OFF: the Infeasible verdict must carry a verifying Farkas vector.
        let lp1 = lp_of(
            &[&[1.0], &[1.0]],
            vec![1.0, 3.0],
            "LG",
            vec![1.0],
            vec![f64::INFINITY; 1],
        );
        let lp2 = lp_of(
            &[&[1.0], &[1.0]],
            vec![1.0, 2.0],
            "EE",
            vec![1.0],
            vec![f64::INFINITY; 1],
        );
        // Infeasibility that only shows through a NEGATED row (rhs < 0):
        // −x ≤ −3 (x ≥ 3) ∧ x ≤ 1 — pins the row-negation mapping of π_f.
        let lp3 = lp_of(
            &[&[-1.0], &[1.0]],
            vec![-3.0, 1.0],
            "LL",
            vec![1.0],
            vec![f64::INFINITY; 1],
        );
        for lp in [&lp1, &lp2, &lp3] {
            let s = solve_np(lp);
            assert_eq!(s.status, SimplexStatus::Infeasible);
            assert!(!s.duals_valid, "Infeasible carries no optimal duals");
            let y = s
                .dual_ray
                .clone()
                .expect("presolve-off Infeasible must attempt a Farkas ray");
            check_farkas(lp, &y);
        }
    }

    /// Regression (`farkas_ok` sup over bounded variables): the per-column
    /// sup of `φ·x_j` over `x_j ∈ [lb, ub]` must be taken at the endpoint the
    /// SIGN of `φ = π_fᵀa_j` maximizes — `+∞` only when THAT side is unbounded.
    /// The first implementation clamped every infinite endpoint to `0.0` and
    /// took `max(φ·lb, φ·ub)` over the finite ones, which OVERESTIMATES the sup
    /// whenever the maximizing side is the finite one and the other side is
    /// infinite (`φ < 0, ub = +∞` ⇒ phantom `0` instead of the true `φ·lb`):
    /// a mathematically valid certificate was then rejected, the debug build
    /// asserted, and the release build returned an honest-but-needless `None`.
    /// These three LPs are the `farkas_certificate_certifies_infeasibility`
    /// shapes with the variable bounded AWAY from the phantom side
    /// (`lower = 2`, `upper = +∞`), so the maximizing endpoint is the lower one.
    #[test]
    fn farkas_ray_survives_bounded_variables() {
        // `x ≤ 1` with `x ≥ 2`: certificate π_f = (-1), lhs = -1, true
        // sup = φ·lb = -2 (φ = -1) ⇒ valid; the clamped estimate said 0.
        let mut lp1 = lp_of(&[&[1.0]], vec![1.0], "L", vec![1.0], vec![f64::INFINITY]);
        lp1.lower = vec![2.0];
        // The L/G pair of `farkas_certificate_certifies_infeasibility` lp1,
        // same bound change: π_f = (-1, 0) gives φ = -1, sup = -2 < lhs = -1.
        let mut lp2 = lp_of(
            &[&[1.0], &[1.0]],
            vec![1.0, 3.0],
            "LG",
            vec![1.0],
            vec![f64::INFINITY],
        );
        lp2.lower = vec![2.0];
        // `x = 1` with `x ≥ 2` (E row ⇒ π_f free): π_f = (-1) gives φ = -1,
        // lhs = -1 > sup = -2.
        let mut lp3 = lp_of(&[&[1.0]], vec![1.0], "E", vec![1.0], vec![f64::INFINITY]);
        lp3.lower = vec![2.0];
        for lp in [&lp1, &lp2, &lp3] {
            let s = solve_np(lp);
            assert_eq!(s.status, SimplexStatus::Infeasible);
            let y = s
                .dual_ray
                .clone()
                .expect("bounded-variable infeasibility must still yield a Farkas ray");
            check_farkas(lp, &y);
        }
    }

    /// Property check on an unbounded direction: A·d respects the row
    /// senses (0/≤0/≥0 for E/L/G), nonzero components move only into
    /// infinite bounds, and cᵀd runs away in the caller's sense.
    fn check_primal_ray(lp: &Lp, d: &[f64]) {
        assert_eq!(d.len(), lp.num_col);
        let cols = dense_cols(lp);
        for i in 0..lp.num_row {
            let ad: f64 = (0..lp.num_col).map(|j| cols[j][i] * d[j]).sum();
            match lp.senses[i] {
                'L' => assert!(ad <= 1e-9, "L row A·d {ad} > 0"),
                'G' => assert!(ad >= -1e-9, "G row A·d {ad} < 0"),
                _ => assert!(ad.abs() <= 1e-9, "E row A·d {ad} ≠ 0"),
            }
        }
        let scale = d.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
        let eps = 1e-8 * scale.max(1.0);
        for (j, &dj) in d.iter().enumerate() {
            if dj > eps {
                assert_eq!(lp.upper[j], f64::INFINITY, "d_{j} > 0, finite ub");
            }
            if dj < -eps {
                assert_eq!(lp.lower[j], f64::NEG_INFINITY, "d_{j} < 0, finite lb");
            }
        }
        let cd: f64 = lp.objective.iter().zip(d).map(|(c, v)| c * v).sum();
        if lp.maximize {
            assert!(cd > 0.0, "maximize but cᵀd = {cd}");
        } else {
            assert!(cd < 0.0, "minimize but cᵀd = {cd}");
        }
    }

    #[test]
    fn primal_ray_certifies_unboundedness() {
        // L row with a diagonal ray; an E row (ray must keep equality);
        // no rows at all with an MI variable (ray along −inf lb); the same
        // first LP maximized (cᵀd flips sign, not the direction cone).
        let lps = [
            lp_of(
                &[&[1.0, -1.0]],
                vec![5.0],
                "L",
                vec![-1.0, -1.0],
                vec![f64::INFINITY; 2],
            ),
            lp_of(
                &[&[1.0, -1.0]],
                vec![0.0],
                "E",
                vec![0.0, -1.0],
                vec![f64::INFINITY; 2],
            ),
            lp_of(&[], vec![], "", vec![-1.0], vec![f64::INFINITY; 1]),
        ];
        for mut lp in lps {
            let s = solve_np(&lp);
            assert_eq!(s.status, SimplexStatus::Unbounded);
            assert!(!s.duals_valid);
            let d = s
                .primal_ray
                .clone()
                .expect("presolve-off Unbounded must carry a ray");
            check_primal_ray(&lp, &d);
            // Same unboundedness seen as a maximization: negate the costs
            // (the ray cone is one-sided — the direction stays, cᵀd flips).
            lp.maximize = true;
            lp.objective = lp.objective.iter().map(|v| -v).collect();
            let s = solve_np(&lp);
            assert_eq!(s.status, SimplexStatus::Unbounded);
            let d = s.primal_ray.clone().expect("max ray");
            check_primal_ray(&lp, &d);
        }
        // MI variable unbounded below (lb = −inf): the pure bound ray.
        let mut lp = lp_of(&[], vec![], "", vec![1.0], vec![10.0]);
        lp.lower = vec![f64::NEG_INFINITY];
        let s = solve_np(&lp);
        assert_eq!(s.status, SimplexStatus::Unbounded);
        let d = s.primal_ray.clone().expect("bound ray");
        check_primal_ray(&lp, &d);
        // Fully free variable: the ±split reconstruction path of the ray.
        let mut lp = lp_of(&[], vec![], "", vec![1.0], vec![f64::INFINITY]);
        lp.lower = vec![f64::NEG_INFINITY];
        let s = solve_np(&lp);
        assert_eq!(s.status, SimplexStatus::Unbounded);
        let d = s.primal_ray.clone().expect("free split ray");
        check_primal_ray(&lp, &d);
    }

    /// Differential sanity: row_duals/reduced_costs vs highspy's
    /// `row_dual`/`col_dual` on three small LPs (fixture committed at
    /// `benches/data/duals/highs_duals.txt`, generated with highspy 1.15.1
    /// the way `tools/diff/simplex_diff.py` builds models). HiGHS reports duals
    /// in the caller's objective sense — the SAME convention we document —
    /// so this compares values directly, no sign flip.
    #[test]
    fn row_duals_match_highspy() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../benches/data/duals/highs_duals.txt"
        );
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let want = |name: &str, key: &str| -> Vec<f64> {
            let line = text
                .lines()
                .find(|l| l.starts_with(&format!("{name} {key} ")))
                .unwrap_or_else(|| panic!("fixture line {name} {key}"));
            line.split_whitespace()
                .skip(2)
                .map(|t| t.parse().unwrap())
                .collect()
        };
        let cases: [(&str, Lp); 3] = [
            // min 3x+2y s.t. x+y >= 5, x,y >= 0
            (
                "min_G",
                lp_of(
                    &[&[1.0, 1.0]],
                    vec![5.0],
                    "G",
                    vec![3.0, 2.0],
                    vec![f64::INFINITY; 2],
                ),
            ),
            // max 3x+2y s.t. x+y <= 4, x,y >= 0
            (
                "max_L",
                lp_of(
                    &[&[1.0, 1.0]],
                    vec![4.0],
                    "L",
                    vec![3.0, 2.0],
                    vec![f64::INFINITY; 2],
                ),
            ),
            // min x+y s.t. x+y = 7 (E), x <= 5 (L row, slack), x in [0,5]
            (
                "min_EL_ub",
                lp_of(
                    &[&[1.0, 1.0], &[1.0, 0.0]],
                    vec![7.0, 5.0],
                    "EL",
                    vec![1.0, 1.0],
                    vec![5.0, f64::INFINITY],
                ),
            ),
        ];
        for (name, mut lp) in cases {
            if name == "max_L" {
                lp.maximize = true;
            }
            let s = solve_np(&lp);
            assert_eq!(s.status, SimplexStatus::Optimal, "{name}");
            let obj_want = want(name, "objective")[0];
            assert!(
                (s.objective - obj_want).abs() <= 1e-9 * obj_want.abs().max(1.0),
                "{name} obj {} want {obj_want}",
                s.objective
            );
            let rd = want(name, "row_duals");
            assert_eq!(rd.len(), s.row_duals.len(), "{name} row count");
            for (i, (a, b)) in s.row_duals.iter().zip(&rd).enumerate() {
                assert!(
                    (a - b).abs() <= 1e-7 * (1.0 + b.abs()),
                    "{name} row_duals[{i}]: ours {a}, highspy {b}"
                );
            }
            let cd = want(name, "col_duals");
            assert_eq!(cd.len(), s.reduced_costs.len(), "{name} col count");
            for (j, (a, b)) in s.reduced_costs.iter().zip(&cd).enumerate() {
                assert!(
                    (a - b).abs() <= 1e-7 * (1.0 + b.abs()),
                    "{name} reduced_costs[{j}]: ours {a}, highspy col_dual {b}"
                );
            }
        }
    }

    #[test]
    fn presolve_on_reports_no_duals_and_no_rays() {
        // Scope contract: the presolve-ON path NEVER reports
        // duals/rays — Optimal, Unbounded and Infeasible alike.
        let s = solve(
            &[&[1.0, 1.0], &[2.0, 1.0]],
            vec![100.0, 150.0],
            "GL",
            vec![2.0, 3.0],
            vec![40.0, f64::INFINITY],
        );
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!(!s.duals_valid);
        assert!(s.row_duals.is_empty() && s.reduced_costs.is_empty());
        assert!(s.primal_ray.is_none() && s.dual_ray.is_none());
        let s = solve(&[], vec![], "", vec![-1.0], vec![f64::INFINITY; 1]);
        assert_eq!(s.status, SimplexStatus::Unbounded);
        assert!(s.primal_ray.is_none());
        let s = solve(
            &[&[1.0], &[1.0]],
            vec![1.0, 3.0],
            "LG",
            vec![1.0],
            vec![f64::INFINITY; 1],
        );
        assert_eq!(s.status, SimplexStatus::Infeasible);
        assert!(s.dual_ray.is_none());
    }

    // ---- Basis export ----

    /// Property check on an exported [`Basis`] against the solution it came
    /// from: nonbasic statuses sit exactly on the caller bound they name,
    /// `Fixed` requires lb == ub, row `AtLower` ⇔ row tight, row `Basic` ⇒
    /// slack ≥ 0, and the implied basic set never exceeds `num_row` (a
    /// shortfall means basic artificials, per the counting contract).
    fn check_basis<'a>(lp: &Lp, s: &'a SimplexSolution) -> &'a Basis {
        let b = s
            .basis
            .as_ref()
            .expect("presolve-off Optimal must export a basis");
        assert_eq!(b.col_status.len(), lp.num_col);
        assert_eq!(b.row_status.len(), lp.num_row);
        let tol = 1e-6;
        let cols = dense_cols(lp);
        let mut n_basic = 0usize;
        for (j, &st) in b.col_status.iter().enumerate() {
            let (lb, ub, x) = (lp.lower[j], lp.upper[j], s.x[j]);
            match st {
                BasisStatus::Basic => n_basic += 1,
                BasisStatus::AtLower => {
                    assert!(lb.is_finite(), "AtLower on -inf lb col {j}");
                    assert!((x - lb).abs() <= tol, "col {j} AtLower but x={x}");
                }
                BasisStatus::AtUpper => {
                    assert!(ub.is_finite(), "AtUpper on +inf ub col {j}");
                    assert!((x - ub).abs() <= tol, "col {j} AtUpper but x={x}");
                }
                BasisStatus::Fixed => {
                    assert_eq!(lb, ub, "Fixed with a range on col {j}");
                    assert!((x - lb).abs() <= tol, "col {j} Fixed but x={x}");
                }
            }
        }
        for (i, &st) in b.row_status.iter().enumerate() {
            let lhs: f64 = (0..lp.num_col).map(|j| cols[j][i] * s.x[j]).sum();
            let slack = match lp.senses[i] {
                'L' => lp.rhs[i] - lhs,
                'G' => lhs - lp.rhs[i],
                _ => 0.0,
            };
            match st {
                BasisStatus::Basic => {
                    if lp.senses[i] != 'E' {
                        assert!(slack >= -tol, "row {i} Basic but slack {slack}");
                    }
                    // 'E' Basic = the row's artificial is basic: a
                    // degenerate redundant equality, nothing primal to check.
                }
                BasisStatus::AtLower => {
                    // 'E' rows legitimately export AtLower (artificial
                    // nonbasic); for a slack it still means row-tight.
                    if lp.senses[i] != 'E' {
                        assert!(slack.abs() <= tol, "row {i} AtLower but slack {slack}");
                    }
                }
                _ => panic!("row status {st:?} never exported"),
            }
        }
        assert!(
            n_basic <= lp.num_row,
            "implied basic set {n_basic} > num_row {}",
            lp.num_row
        );
        b
    }

    #[test]
    fn basis_export_matches_solution_mixed_structure() {
        // Every column/row class at once: finite-lb structural, positive lb,
        // upper-bound-active column, flip-shift (lb=-inf), fully-free split,
        // fixed column, L/G/E rows including a negated-rhs row.
        let mut lp = lp_of(
            &[
                &[1.0, 1.0, 0.0, 1.0, 0.0, 1.0],
                &[2.0, 0.0, 1.0, 0.0, 1.0, 0.0],
                &[-1.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                &[1.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            ],
            vec![8.0, 10.0, -2.0, 5.0],
            "LLGE",
            vec![1.0, -2.0, 3.0, 1.0, -1.0, 0.5],
            vec![6.0, f64::INFINITY, 4.0, 2.0, f64::INFINITY, 2.0],
        );
        lp.lower = vec![0.0, 1.0, f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY, 2.0];
        let s = solve_np(&lp);
        assert_eq!(s.status, SimplexStatus::Optimal);
        check_basis(&lp, &s);
    }

    #[test]
    fn basis_export_none_on_presolve_on_path() {
        // Row x1 + x2 <= 100 is dominated by the box bounds ⇒ presolve
        // reduces (the Reduced arm) ⇒ basis stays None (scope contract,
        // same loud pattern as the duals).
        let lp = lp_of(
            &[&[1.0, 1.0], &[1.0, 1.0]],
            vec![100.0, 6.0],
            "LL",
            vec![-1.0, -1.0],
            vec![2.0, 4.0],
        );
        let s = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!(
            s.basis.is_none(),
            "presolve-ON must not export a reduced-LP basis"
        );
    }

    #[test]
    fn dual_engine_exports_basis() {
        // The dual lever is the other presolve-OFF path: its Optimal must
        // carry a caller-consistent snapshot too.
        let lp = lp_of(
            &[&[1.0, 1.0], &[2.0, 1.0]],
            vec![5.0, 8.0],
            "GG",
            vec![3.0, 2.0],
            vec![10.0, 10.0],
        );
        let s = SimplexSolver::new(SimplexConfig {
            presolve: false,
            dual: true,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        check_basis(&lp, &s);
    }

    // ---- Warm primal start ----

    fn cfg_np() -> SimplexConfig {
        SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        }
    }

    #[test]
    fn warm_reroot_zero_extra_iterations_identical_x() {
        // The export acceptance test, wired through the real import path: re-solve
        // with the exported basis ⇒ 0 extra pivots, identical x/objective.
        let lp = lp_of(
            &[&[1.0, 1.0], &[2.0, 1.0]],
            vec![10.0, 14.0],
            "LL",
            vec![-3.0, -2.0],
            vec![6.0, f64::INFINITY],
        );
        let solver = SimplexSolver::new(cfg_np());
        let cold = solver.solve(&lp).unwrap();
        assert_eq!(cold.status, SimplexStatus::Optimal);
        let b = cold.basis.clone().expect("export");
        let warm = solver.solve_with(&lp, Some(&b)).unwrap();
        assert_eq!(warm.status, SimplexStatus::Optimal);
        assert_eq!(warm.phase1_iterations, 0, "warm: no phase 1");
        assert_eq!(
            warm.phase2_iterations, 0,
            "re-root on the same basis adds no pivots"
        );
        assert_eq!(warm.x.len(), cold.x.len());
        for (a, x) in warm.x.iter().zip(&cold.x) {
            assert!((a - x).abs() < 1e-12, "x drift {a} vs {x}");
        }
        assert!((warm.objective - cold.objective).abs() < 1e-12);
        // The warm run is itself an Optimal presolve-off solve ⇒ it exports
        // a fresh snapshot, which round-trips again (MIP grandchild case).
        let b2 = warm.basis.clone().expect("warm export");
        assert_eq!(b2, b, "same LP + same basis ⇒ same snapshot");
    }

    #[test]
    fn warm_after_rhs_tweak_fewer_iters_same_objective() {
        // min −x s.t. x ≤ 10 ('L' row), x ≤ 5: optimum x=5 with row slack 5
        // basic. Tightening rhs to 6 keeps the basis primal feasible (slack
        // 5→1): warm re-anchors with 0 pivots and the objective is unchanged;
        // cold pays a fresh phase 1.
        let mut lp = lp_of(&[&[1.0]], vec![10.0], "L", vec![-1.0], vec![5.0]);
        let solver = SimplexSolver::new(cfg_np());
        let root = solver.solve(&lp).unwrap();
        let b = root.basis.clone().expect("export");
        lp.rhs[0] = 6.0;
        let warm = solver.solve_with(&lp, Some(&b)).unwrap();
        let cold = solver.solve(&lp).unwrap();
        assert_eq!(warm.status, SimplexStatus::Optimal);
        assert_eq!(cold.status, SimplexStatus::Optimal);
        assert_eq!(warm.phase1_iterations + warm.phase2_iterations, 0);
        assert!(cold.phase1_iterations + cold.phase2_iterations > 0);
        assert!(
            (warm.objective - cold.objective).abs() < 1e-9,
            "warm {} vs cold truth {}",
            warm.objective,
            cold.objective
        );
        assert!(
            (warm.objective - root.objective).abs() < 1e-9,
            "the rhs move is absorbed by the basic slack"
        );
    }

    #[test]
    fn warm_after_cost_change_repaired_by_phase2() {
        // max x + 2y s.t. x + y ≤ 10, 0 ≤ x ≤ 6, 0 ≤ y ≤ 7: root (3, 7).
        // Costs → max x + 0.5y: the root basis stays primal feasible but
        // stops being optimal — warm phase 2 repairs it in a couple of
        // pivots; truth = cold solve of the changed LP.
        let mut lp = lp_of(
            &[&[1.0, 1.0]],
            vec![10.0],
            "L",
            vec![-1.0, -2.0],
            vec![6.0, 7.0],
        );
        let solver = SimplexSolver::new(cfg_np());
        let root = solver.solve(&lp).unwrap();
        let b = root.basis.clone().expect("export");
        lp.objective = vec![-1.0, -0.5];
        let warm = solver.solve_with(&lp, Some(&b)).unwrap();
        let cold = solver.solve(&lp).unwrap();
        assert_eq!(warm.status, SimplexStatus::Optimal);
        assert_eq!(cold.status, SimplexStatus::Optimal);
        assert!(
            (warm.objective - cold.objective).abs() < 1e-9,
            "warm {} vs cold truth {}",
            warm.objective,
            cold.objective
        );
        let wi = warm.phase1_iterations + warm.phase2_iterations;
        let ci = cold.phase1_iterations + cold.phase2_iterations;
        assert!(
            warm.phase1_iterations == 0 && (1..ci).contains(&wi),
            "warm {wi} (phase1 {}) should beat cold {ci}",
            warm.phase1_iterations
        );
    }

    #[test]
    fn warm_rejects_primal_infeasible_basis_falls_back_cold() {
        // rhs 10→4 makes the root basis primal INFEASIBLE (basic slack at
        // −1): the honest behavior is a cold fallback — right answer, and
        // the re-run phase 1 is visible in the iteration counts.
        let mut lp = lp_of(&[&[1.0]], vec![10.0], "L", vec![-1.0], vec![5.0]);
        let solver = SimplexSolver::new(cfg_np());
        let root = solver.solve(&lp).unwrap();
        let b = root.basis.clone().expect("export");
        lp.rhs[0] = 4.0;
        let warm = solver.solve_with(&lp, Some(&b)).unwrap();
        let cold = solver.solve(&lp).unwrap();
        assert_eq!(warm.status, SimplexStatus::Optimal);
        assert!(
            (warm.objective - cold.objective).abs() < 1e-9,
            "warm {} vs cold truth {}",
            warm.objective,
            cold.objective
        );
        assert!(
            warm.phase1_iterations > 0,
            "fallback evidence: phase 1 re-ran"
        );
    }

    #[test]
    fn warm_ignores_shape_mismatched_snapshot() {
        // A snapshot from a differently-shaped LP violates the contract;
        // the solve must come back correct via the cold path, not crash.
        let lp = lp_of(
            &[&[1.0, 1.0]],
            vec![4.0],
            "L",
            vec![-1.0, -1.0],
            vec![f64::INFINITY; 2],
        );
        let other = lp_of(
            &[&[1.0, 1.0], &[1.0, 0.0]],
            vec![4.0, 2.0],
            "LL",
            vec![-1.0, -1.0],
            vec![f64::INFINITY; 2],
        );
        let solver = SimplexSolver::new(cfg_np());
        let b = solver.solve(&other).unwrap().basis.clone().expect("export");
        let s = solver.solve_with(&lp, Some(&b)).unwrap();
        let truth = solver.solve(&lp).unwrap();
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - truth.objective).abs() < 1e-9);
    }

    // ---- Solve-internal deadline (HiGHS HEkk::bailout shape) ----

    #[test]
    fn deadline_expired_bails_out_of_the_pivot_loop_immediately() {
        // The tree's limit checks fire only between
        // nodes, so a solve-internal cycle must bail on its own clock. An
        // expired deadline ends the solve on its FIRST pivot-loop pass,
        // far short of `max_iterations`.
        let lp = lp_of(
            &[&[1.0, 1.0], &[2.0, 1.0]],
            vec![10.0, 14.0],
            "LL",
            vec![-3.0, -2.0],
            vec![6.0, f64::INFINITY],
        );
        let mut solver = SimplexSolver::new(cfg_np());
        solver.set_deadline(Some(std::time::Instant::now()));
        let sol = solver.solve(&lp).unwrap();
        assert_eq!(sol.status, SimplexStatus::IterationLimit);
        assert!(
            sol.phase1_iterations + sol.phase2_iterations <= 1,
            "bail-out is immediate (iterations = {} + {}), not at the cap",
            sol.phase1_iterations,
            sol.phase2_iterations
        );
    }

    #[test]
    fn deadline_in_the_future_solves_normally() {
        let lp = lp_of(
            &[&[1.0, 1.0], &[2.0, 1.0]],
            vec![10.0, 14.0],
            "LL",
            vec![-3.0, -2.0],
            vec![6.0, f64::INFINITY],
        );
        let mut solver = SimplexSolver::new(cfg_np());
        solver.set_deadline(Some(
            std::time::Instant::now() + std::time::Duration::from_secs(3600),
        ));
        let sol = solver.solve(&lp).unwrap();
        assert_eq!(sol.status, SimplexStatus::Optimal);
        assert!((sol.objective - -24.0).abs() < 1e-9);
    }

    // ---- The dual engine honors the same solve-internal deadline ----

    #[test]
    fn dual_engine_honors_the_deadline() {
        // The LP CLI runs uncapped and bounded by `--time-limit`, and the
        // CLI default engine is the dual — so the dual pivot loop must read
        // the clock too. Expired ⇒ immediate IterationLimit-status
        // exit; a future deadline is a no-op.
        let lp = lp_of(
            &[&[1.0, 1.0], &[2.0, 1.0]],
            vec![10.0, 14.0],
            "LL",
            vec![-3.0, -2.0],
            vec![6.0, f64::INFINITY],
        );
        let cfg = SimplexConfig {
            dual: true,
            ..cfg_np()
        };
        let mut solver = SimplexSolver::new(cfg);
        solver.set_deadline(Some(std::time::Instant::now()));
        let sol = solver.solve(&lp).unwrap();
        assert_eq!(sol.status, SimplexStatus::IterationLimit);
        assert!(sol.phase1_iterations + sol.phase2_iterations <= 1);
        solver.set_deadline(Some(
            std::time::Instant::now() + std::time::Duration::from_secs(3600),
        ));
        let sol = solver.solve(&lp).unwrap();
        assert_eq!(sol.status, SimplexStatus::Optimal);
        assert!((sol.objective - -24.0).abs() < 1e-9);
    }
}
