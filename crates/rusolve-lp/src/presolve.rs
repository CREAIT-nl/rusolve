//! Presolve: reduction rules applied to fixpoint before the simplex, plus
//! postsolve restoring the full-space solution.
//!
//! Why this exists: eliminating whole degenerate families before the simplex
//! is often worth more than any pricing refinement.
//!
//! Rules (all applied per pass; loop to fixpoint, capped):
//! * **infeasible row**: activity range cannot reach rhs ⇒ Infeasible;
//! * **empty row**: no nonzeros and satisfied ⇒ drop;
//! * **redundant row**: its best-serving activity bound already satisfies the
//!   sense ⇒ drop;
//! * **forcing row**: satisfiable only at an activity extreme ⇒ fix every
//!   involved variable at the bound achieving that extreme, drop the row;
//! * **single-nonzero row**: `a_j x_j (sense) rhs` tightens that variable's
//!   bound, or fixes it (E row);
//! * **interval bound tightening** on every row (propagated bounds move only
//!   by a relative gain >= `PROP_MIN_GAIN`, which stops Zeno cascades);
//! * **aggregator** (HiGHS `HPresolve::aggregator`, MIT): an implied-free
//!   column is substituted out through an equality row (fill-in capped,
//!   pivot-thresholded);
//! * **fixed variable** (`lb == ub`, incl. via tightening): substitute —
//!   `rhs -= a·v`, `obj_const += c·v`, column removed;
//! * **empty column** (no constraint entries): optimum sits at a bound
//!   (min: `c < 0 ⇒ ub`, else `lb`; symmetric for max); infinite preferred
//!   bound ⇒ Unbounded;
//! * **implied-equality conversion** (HiGHS port): a ranged row whose
//!   implied activity touches the far side becomes an equality for the rest
//!   of presolve, feeding the doubleton-equality chain;
//! * **pure-bound row fold** (HiGHS port): a row whose nonzeros all
//!   sit on fixed variables except one is a pure bound on that variable —
//!   fold the bound in, delete the row.
//! * **forcing-column removal** (Andersen & Andersen 1995 + HiGHS
//!   port): a column with one infinite bound whose dual cost is pinned at
//!   ≈0 over the sense-only dual box is deleted together with ALL its rows
//!   (never via the fixed-value path — it is value-unknown at presolve
//!   time); postsolve reconstructs it at the far-side-overlap value.
//! * **parallel-row merge** (HiGHS port): rows with
//!   proportional coefficients and identical support are merged by
//!   effective-interval intersection (same direction: the dominating scaled
//!   rhs wins; equalities inside the interval win outright; a collapsed
//!   interval [d, d] converts its lower witness to an equality — the
//!   graph20 class; a true range keeps one witness per side), the losers
//!   deleted; no postsolve record (rows carry no primal value and
//!   a loser is implied by the kept row).
//! * **duplicate-column merge** (HiGHS port): cost-equivalent columns
//!   with proportional coefficients (any finite scale) are merged
//!   through
//!   `merged = x_kept + scale·x_dup` (bounds widen to the merged box, rows
//!   untouched); postsolve splits the merged value back into the two bound
//!   boxes (`DupColRec`).
//! * **parallel-column dominance fixing** (HiGHS
//!   `detectParallelRowsAndCols` :7399–7633 LP arms, MIT): one of a
//!   proportional column pair with unequal costs is fixed at a bound when
//!   the other can absorb the compensation (opposing bound infinite or
//!   row-implied); value-known ⇒ rides the `fixed` substitution.
//! * **singleton-cost VUB implied-equality chain** (HiGHS
//!   `isImpliedEquation` :263–281 via a column-singleton dual bound, MIT):
//!   a cost-carrying column singleton whose improving bound is implied by
//!   the current system pins its row's dual away from zero ⇒ the row is
//!   tight at every optimum ⇒ converted to an equation and eaten by the
//!   doubleton chain in the same pass.
//!
//! Provenance: the rule set and the forcing/redundancy activity arithmetic
//! follow the published presolve paper (Achterberg, "Presolve", ZIB Report
//! 18-01, opus4.kobv.de/opus4-zib/files/6037) and HiGHS' presolve
//! documentation (MIT). Interval tightening is textbook row-wise bound
//! propagation. The implied-equality conversion and the pure-bound row fold
//! are ports of HiGHS `HPresolve.cpp` (MIT): see the attribution headers at
//! each rule. The doubleton-rule choice + stay-col bounds are likewise a
//! port of HiGHS `HPresolve.cpp::doubletonEq` :3036–3203 (MIT), LP branch
//! only.
//!
//! Deferred: equality-chain aggregation, row/col dominance, scaling,
//! dual postsolve (presolve emits no duals at all). General-scale (|s| ≠ 1)
//! duplicate-column merging is follow-up work (detection already
//! reports those groups).
//!
//! Column identity: the working LP holds only LIVE columns; `live[j]` maps a
//! work column to its original variable index.

use super::simplex::SimplexStatus;
use rusolve_core::Lp;

/// Feasibility/forcing tolerance for presolve decisions (absolute floor).
const TOL: f64 = 1e-9;

/// Minimum relative gain for a propagated bound to be applied (see `prop_ub`).
const PROP_MIN_GAIN: f64 = 1e-7;

/// Dual feasibility tolerance for dominated-column detection: the
/// HiGHS `dual_feasibility_tolerance` default (1e-7, HConst.h:344–345).
/// Distinct from the primal `TOL` above, which governs forcing/redundancy
/// decisions and stays untouched.
const DUAL_TOL: f64 = 1e-7;

/// Scale tie-break tolerance for the parallel detection: HiGHS'
/// `small_matrix_value` (HConst.h), used there to decide which signed
/// max-|coef| nonzero defines a row/col's scale (:7298–7320 tie-break).
const SMALL_MATRIX: f64 = 1e-12;

/// Parallelism tolerance on NORMALIZED coefficients: HiGHS re-checks
/// candidate pairs against `small_matrix_value` on RAW values (:7620–7632);
/// ours is the normalized-coefficient form, deliberately looser at 1e-9 but
/// applied to values normalized into [−1, 1] — conservative overall.
const PAR_TOL: f64 = 1e-9;

/// Relative slack for forcing/redundancy verdicts: an extreme activity
/// counts as "at" rhs only if the gap is tiny *relative to the data scale*
/// (an absolute tol false-forces on large-coefficient rows — boeing2 showed
/// a false forcing with gap ~2e-11 under 1e-9). Post-scaling data is O(1),
/// so 1e-12 stays above fp accumulation noise while rejecting near-misses.
fn at_rhs(r: f64, a: f64) -> bool {
    (a - r).abs() <= 1e-12 * r.abs().max(1.0)
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum VarState {
    /// Original variable substituted away at this exact value.
    Fixed(f64),
    /// Not yet resolved (becomes `Live` in the final mapping).
    Pending,
}

/// One eliminated two-term equality: `x_orig = (b − am·x_partner)/ak`,
/// replayed in reverse elimination order by [`restore`].
#[derive(Debug, Clone, Copy)]
pub struct ChainRec {
    orig: usize,
    partner: usize,
    ak: f64,
    am: f64,
    b: f64,
    /// Sense this row carried **before** any implied-equality conversion:
    /// `'E'` if it was always an equation, else `'L'`/`'G'` — the
    /// honest record that presolve tightened a ranged row into an equality
    /// here (HiGHS does the same via `HighsPostsolveStack::RowType` at
    /// HPresolve.cpp:3596–3607). Primal reconstruction needs nothing from
    /// it — x satisfies the equality, hence the original inequality — so
    /// this is dual-postsolve groundwork (we emit no duals).
    #[allow(dead_code)] // read once dual postsolve lands
    orig_sense: char,
}

/// One row of a deleted forcing column, captured in its CURRENT (post-fold)
/// state at removal time: original row index + sense (assertion context),
/// rhs at removal (already includes everything folded in so far), and the
/// (original col, coefficient) entries of ALL live columns of the row —
/// including the forcing column itself (its entry is its coefficient).
/// Columns already substituted are NOT entries: their contribution lives in
/// `rhs` and must not be re-added at reconstruction.
#[derive(Debug, Clone)]
pub struct ForcedRow {
    /// Original row index (assertion/log context only; reconstruction needs
    /// nothing from it).
    row: usize,
    /// Row sense at removal time.
    sense: char,
    /// rhs at removal time (later rhs moves on this row cannot happen: the
    /// row is deleted with the column).
    rhs: f64,
    /// (original col index, coefficient) for every live column of the row,
    /// the forcing column included.
    entries: Vec<(usize, f64)>,
}

/// One eliminated FORCING column (Andersen & Andersen 1995): the column had one infinite bound and dual
/// cost ≈ 0, so its rows only bound it from the finite side. Reconstruction
/// (reverse replay, primal-only) picks x_col as the far-side-overlap value:
/// min(finite_bound, all row upper bounds) for dir +1, symmetric max for
/// dir −1. Not the `fixed` path: the column is never value-known at presolve
/// time.
#[derive(Debug, Clone)]
pub struct ForcingRec {
    /// Original index of the removed column.
    orig: usize,
    /// +1: lower bound was −∞ (rows bound x_col from ABOVE, pick the min);
    /// −1: upper bound was +∞ (rows bound BELOW, pick the max).
    dir: i8,
    /// The finite bound that survived (ub for dir +1, lb for dir −1).
    finite_bound: f64,
    /// The rows deleted together with the column (HiGHS removeRow :4611).
    rows: Vec<ForcedRow>,
}

/// One merged DUPLICATE column (HiGHS `duplicateColumn` +
/// DuplicateColumn::undo, MIT): column `dup` (original index `orig`) left
/// the LP through the identity `merged = x_kept + scale·x_dup` with any
/// finite nonzero scale (HiGHS `duplicateColumn`
/// HighsPostsolveStack.h:513 accepts any real colScale for LP — the ±1
/// restriction was the unported integer-arm machinery; guard: finite scale).
/// Because the duplicate's
/// coefficients are `scale·` the kept column's on EVERY row, the mutation
/// touches no row coefficients at all: it only widens the kept column's
/// bounds to the merged box and drops the duplicate. Reconstruction splits
/// the merged value m back in HiGHS's ANCHOR form
/// (`DuplicateColumn::undo`, HighsPostsolveStack.cpp:700–737):
/// `anchor = kept_lo ≠ −inf ? kept_lo : min(0, kept_hi)` (finite by
/// construction), `x_dup = trim((m − anchor)/scale → dup box)`,
/// `x_kept = trim fired ? m − scale·x_dup : anchor`.
/// Case-checked for any finite scale ≠ 0 and ANY bounds (the ±1 argument
/// generalizes verbatim; the anchor step is REQUIRED for infinite
/// merge-time lower bounds — the naive (m − kept_lo)/scale evaluates to
/// ±inf for a free kept column and the trim never fires, leaking ±inf
/// into both halves; caught by the agg fuzz, seed 90): m lies in the
/// stored merged box (the merged
/// variable's declared bounds, which only tighten after the merge), and the
/// split then lands BOTH parts inside their stored merge-time boxes — hence
/// inside the ORIGINAL bounds (merge-time bounds are tightened originals):
/// an untrimmed split gives x_kept = anchor (a box point by construction);
/// a trim that FIRES implies the
/// inequality its side of x_kept needs (e.g. scale>0, y = dup_hi fires
/// only when (m−anchor)/scale > dup_hi ⇒ x_kept > kept_lo, and the merged
/// box's own upper gives x_kept ≤ kept_hi).
/// Rows see `a·x_kept + a_dup·x_dup = a·(x_kept + scale·x_dup) = a·m` up to
/// the PAR_TOL parallelism slack. Cost drift ≤ DUAL_TOL·|x_dup| (the objDiff
/// gate) — the same accepted price as the forcing-column removal.
#[derive(Debug, Clone)]
pub struct DupColRec {
    /// Original index of the removed (duplicate) column.
    orig: usize,
    /// Original index of the surviving (kept) column, holding the merged
    /// value at replay time.
    kept: usize,
    /// Finite nonzero scale: `merged = x_kept + scale·x_dup`.
    scale: f64,
    /// Kept column's bounds at merge time = merged variable's bounds BEFORE
    /// this merge (needed by the split; later tightening only shrinks m).
    kept_lo: f64,
    /// Upper mate of `kept_lo`. The split formula reads only `kept_lo`
    /// (m ≤ stored box upper holds automatically: the merged column's own
    /// bounds enforce it), so this is dual-postsolve groundwork — same
    /// posture as `ChainRec::orig_sense`.
    #[allow(dead_code)] // read once dual postsolve lands
    kept_hi: f64,
    /// Duplicate column's bounds at merge time.
    dup_lo: f64,
    dup_hi: f64,
}

/// One column eliminated by the AGGREGATOR (free-column substitution through
/// an equality row; HiGHS `HPresolve::aggregator` :6688 + `substituteFreeCol`
/// :3420, MIT). Row `r: a·x_orig + Σ_k a_k·x_k = b` was deleted together with
/// the column; the column is recovered from the row identity
/// `x_orig = (b − Σ_k a_k·x_k)/a` over the row's OTHER live columns at
/// elimination time (columns fixed earlier live in `b`; columns eliminated
/// later are resolved first by the reverse replay). Because the substitution
/// is only done when the implied range of x_orig over the remaining bounds
/// already lies inside its own bounds (`agg_is_free`), the recovered value
/// satisfies them and every other row sees the identical activity as before.
#[derive(Debug, Clone)]
pub struct AggRec {
    orig: usize,
    /// Pivot coefficient of `orig` in the deleted equality.
    a: f64,
    /// rhs of the deleted equality at elimination time.
    b: f64,
    /// (original col, coefficient) of the other live columns of the row.
    rest: Vec<(usize, f64)>,
}

/// A postsolve elimination record, replayed by [`restore`] in reverse
/// creation order. ONE vector, interleaved: a chain record's partner (or a
/// forcing record's stored column) may have been eliminated by a record of
/// the OTHER kind created later, and reverse order over the merged sequence
/// is what guarantees every dependency is resolved before it is read.
#[derive(Debug, Clone)]
enum PostRec {
    Chain(ChainRec),
    Forcing(ForcingRec),
    DupCol(DupColRec),
    Agg(AggRec),
}

/// Maps a reduced-space solution back to the original space.
#[derive(Debug, Clone)]
pub struct Postsolve {
    /// Original index → fixed value (only meaningful for originals absent
    /// from `live`; live originals resolve via the reduced solution).
    fixed: Vec<Option<f64>>,
    /// Constant the substitutions added to the objective (caller space).
    obj_const: f64,
    /// All eliminations (chain + forcing), in elimination order.
    recs: Vec<PostRec>,
}

/// Outcome of presolve.
#[derive(Debug)]
// Reduced carries the whole LP; boxing it only taxes the hot path for size.
#[allow(clippy::large_enum_variant)]
pub enum Presolved {
    Reduced {
        lp: Lp,
        /// Reduced column → original variable index.
        live: Vec<usize>,
        postsolve: Postsolve,
    },
    Terminal {
        status: SimplexStatus,
    },
}

/// Run the reduction loop. `max_passes` guards against limit cycles.
/// `implied_eq` gates the implied-equality conversion + far-side
/// forcing block: ported and unit-tested; measured harmful on Netlib
/// BEFORE (scrs8 p2 312→1830) and AGAIN after the doubletonEq port gave
/// converted rows a chain path (scrs8 p2 312→1827, e226 +51%) ⇒ default
/// off, lever stays.
/// `parcol` gates the duplicate-column merge: ported, unit-tested,
/// feasibility-verified on every firing instance; measured NET-NEGATIVE
/// on Netlib (15,560 → 16,050 iters: tuff +21.0% unexplained, standmps
/// +7.4%, standata +6.5% vs scorpion −16.4%) ⇒ default off,
/// `RUSOLVE_PS_PARCOL=1` is the A/B lever (row merge is separate and
/// default ON, lever `RUSOLVE_PS_NO_PARROW=1`).
/// `agg` gates the aggregator (implied-free column substitution):
/// default ON, net win on the 47 MIPLIB relaxations and Netlib;
/// `RUSOLVE_PRESOLVE_AGG=0` is the CLI lever.
/// `vub_eq` gates the singleton-cost VUB implied-equality chain (see
/// `singleton_cost_eq_conversions`): a cost-carrying column singleton
/// whose improving bound is row-implied pins its row's dual away from
/// zero, so the row is tight at every optimum and is converted to an
/// equation for the DTE chain (mc11/beasleyC3 class);
/// `RUSOLVE_PS_NO_VUBEQ=1` is the
/// CLI lever.
pub fn run(
    lp: &Lp,
    max_passes: usize,
    implied_eq: bool,
    parcol: bool,
    agg: bool,
    parcol_dom: bool,
    vub_eq: bool,
) -> Presolved {
    // Memory guard for the largest LPs (measured on our benchmark set):
    // the reduction
    // loop's working set scales with nnz, and above the budget the
    // reduction is skipped outright. Levers (repo convention):
    // RUSOLVE_PS_NO_MEMGUARD=1 disables; RUSOLVE_PS_MEMGUARD_BYTES=<n>
    // overrides the budget (0 = disabled) for A/B.
    let budget = match std::env::var("RUSOLVE_PS_MEMGUARD_BYTES") {
        Ok(v) => v.trim().parse::<u64>().unwrap_or(PS_GUARD_BUDGET_BYTES),
        Err(_) => PS_GUARD_BUDGET_BYTES,
    };
    let budget = if std::env::var("RUSOLVE_PS_NO_MEMGUARD").is_ok() {
        0
    } else {
        budget
    };
    run_with_budget(
        lp, max_passes, implied_eq, parcol, agg, parcol_dom, vub_eq, budget,
    )
}

/// Bytes-per-nnz slope of the fitted presolve working-set estimate.
/// Measured presolve-phase peak RSS across 17 large LP and MIPLIB instances,
/// 95k–20.1M nnz:
/// 137–1550 B/nnz, worst bdry2 (1.50M nnz, 2.32 GB), highest a2864
/// (20.08M nnz, 3.01 GB, 170 B/nnz). The 800 B/nnz + base estimate sits
/// ABOVE every measured instance (max est 15.9 GiB at a2864 vs 3.01 GB
/// actual) and above the whole corpus (max nnz = a2864; tpl-tub-ws1617's
/// 1.8 GB FILE is name-string volume, nnz only 4.7M → correctly kept).
pub(crate) const PS_GUARD_BYTES_PER_NNZ: u64 = 800;
/// Fixed head: reader, name maps, per-row/col vector headers.
pub(crate) const PS_GUARD_BASE_BYTES: u64 = 256 << 20;
/// Default budget: 16 GiB ⇒ skip threshold ≈ 21.3M nnz — just above the
/// corpus max (a2864, 20.08M nnz), so the guard is a NO-OP on netlib
/// (max nnz 4,316), the 47 relaxed MIPLIB LPs (384k), all of our MIPLIB set
/// (s100 1.78M, bohle-class roots est < 16 GiB) and every large-LP
/// instance we test, INCLUDING the big-nnz
/// points whose presolve is measured cheap AND valuable: scpm1 (6.25M,
/// 1.38 GB), square41 (13.57M, 1.87 GB, rows −96%),
/// a2864 (20.08M, 3.01 GB, cols −93%). Purely preventive insurance for
/// inputs beyond the measured envelope (a bdry2-class outlier at 21M nnz
/// would extrapolate to ~30 GB). The 4 GiB candidate was REJECTED by
/// measurement: it would skip square41/a2864-class presolve whose actual
/// cost is < 4 GB and whose reductions are decisive.
pub(crate) const PS_GUARD_BUDGET_BYTES: u64 = 16 << 30;

/// Estimated presolve working set (bytes) for an LP with `nnz` entries:
/// fitted slope + fixed head, see the constants' docs.
fn memguard_estimate(nnz: usize) -> u64 {
    PS_GUARD_BYTES_PER_NNZ
        .saturating_mul(nnz as u64)
        .saturating_add(PS_GUARD_BASE_BYTES)
}

/// Pure skip decision: true when the estimated presolve working set
/// exceeds `budget` bytes (budget 0 = guard disabled). Split out so
/// tests need no process env (edition-2024 `set_var` is unsafe).
fn memguard_skip(nnz: usize, budget: u64) -> bool {
    budget > 0 && memguard_estimate(nnz) > budget
}

/// [`run`] with an explicit guard budget (bytes; 0 disables). `run`'s
/// env-lever reading is the only thing this adds over the reduction loop.
#[allow(clippy::too_many_arguments)] // mirrors `run`'s gate-by-gate signature
pub(crate) fn run_with_budget(
    lp: &Lp,
    max_passes: usize,
    implied_eq: bool,
    parcol: bool,
    agg: bool,
    parcol_dom: bool,
    vub_eq: bool,
    guard_budget: u64,
) -> Presolved {
    let num_orig = lp.num_col;
    if memguard_skip(lp.a.nnz(), guard_budget) {
        // Identity skip: hand back the input LP with an empty postsolve.
        // Downstream this is bit-identical to the presolve-OFF path (same
        // solve_dispatch on a bit-identical matrix; `restore` degenerates
        // to x[orig] = reduced_x[orig], obj + 0), so the guard can only
        // trade reduction benefit for bounded peak RSS, never correctness.
        eprintln!(
            "[presolve] SKIPPED: nnz {} exceeds memory-guard budget (est {:.2} GiB > {:.2} GiB)",
            lp.a.nnz(),
            memguard_estimate(lp.a.nnz()) as f64 / (1u64 << 30) as f64,
            guard_budget as f64 / (1u64 << 30) as f64,
        );
        return Presolved::Reduced {
            lp: lp.clone(),
            live: (0..num_orig).collect(),
            postsolve: Postsolve {
                fixed: vec![None; num_orig],
                obj_const: 0.0,
                recs: Vec::new(),
            },
        };
    }
    // Sparse working matrix with STABLE ids (see `Sp`): column j / row i keep
    // their original index for the whole run, deletions only clear alive
    // flags. Because deletion preserves relative order, "ascending id" is
    // exactly the "ascending position" order the former dense implementation
    // iterated in, so every sum is accumulated in the identical order.
    let t_start = std::time::Instant::now();
    let mut work = Sp::from_lp(lp);
    let mut states: Vec<VarState> = vec![VarState::Pending; num_orig];
    let mut obj_const = 0.0f64;
    let mut recs: Vec<PostRec> = Vec::new();
    // Per ORIGINAL row: the sense it had in the input, before any
    // implied-equality conversion. Consumed by the DTE postsolve record.
    let osenses: Vec<char> = lp.senses.clone();
    // Effectiveness counters for the RUSOLVE_PRESOLVE_STATS=1 probe lever.
    let stats = std::env::var("RUSOLVE_PRESOLVE_STATS").is_ok();
    // A/B levers (repo convention, cf. RUSOLVE_NO_SH): isolate the
    // pure-bound fold for per-instance regression hunting. The implied-eq
    // block's gate is the `implied_eq` parameter (config-level, default
    // off — see run's doc), not an env read here.
    let no_fold = std::env::var("RUSOLVE_PS_NO_FOLD").is_ok();
    // A/B lever: disables the forcing-column MUTATION only (the
    // detection probe under RUSOLVE_PRESOLVE_STATS keeps counting).
    let no_forcecol = std::env::var("RUSOLVE_PS_NO_FORCECOL").is_ok();
    // A/B lever: disables the parallel-row MERGE mutation only
    // (the detection probe keeps counting under the stats flag).
    let no_parrow = std::env::var("RUSOLVE_PS_NO_PARROW").is_ok();
    // A/B lever: disables ONLY the cross-sense interval-collapse
    // arm (L/G pair whose effective interval reduces to a point ⇒ 'E'); the
    // same-sense dominance merge and per-side loser deletion stay under
    // RUSOLVE_PS_NO_PARROW.
    let parrow_eq = std::env::var("RUSOLVE_PS_NO_PARROWEQ").is_err();
    let debug = std::env::var("RUSOLVE_DEBUG").is_ok();
    let mut n_agg = 0usize;
    let mut n_folded = 0usize;
    let mut n_eq_converted = 0usize;
    let mut n_forced = 0usize;
    // Dominated-col probe counters (filled only under the stats flag).
    let mut n_dom_lb = 0usize;
    let mut n_force_lb = 0usize;
    let mut n_force_ub = 0usize;
    let mut n_dom_ub = 0usize;
    let mut n_near_miss = 0usize;
    let mut n_dual_inf = 0usize;
    // Parallel row/col probe counters (stats-gated). PEAK per-pass,
    // not cumulative: the same static duplicate structure is re-detected
    // every pass and must not be double-counted.
    let mut par_row_groups = 0usize;
    let mut par_row_rows = 0usize;
    let mut par_col_groups = 0usize;
    let mut par_col_cols = 0usize;
    // Merge counter (rows deleted by the dominance merge).
    let mut n_parrow = 0usize;
    // Cross-sense collapse conversions (rows flipped to 'E' by the
    // interval-point arm).
    let mut n_parrow_eq = 0usize;
    // Merge counter (duplicate columns merged away).
    let mut n_parcol = 0usize;
    // Parallel-column dominance counter (gate = `parcol_dom` param;
    // env lever RUSOLVE_PS_NO_PARCOLDOM=1 is read by the CLI).
    let mut n_parcol_dom = 0usize;
    // Singleton-cost VUB implied-equality counter (gate = `vub_eq`
    // param; env lever RUSOLVE_PS_NO_VUBEQ=1 is read by the CLI).
    let mut n_vub_eq = 0usize;

    for _pass in 0..max_passes {
        work.compact();
        let sig_before = work.signature();

        // ---------- singleton-cost VUB implied-equality chain ----------
        // Census: HiGHS removes 1520
        // rows + 1520 cols on mc11 and ~1145/1277 on beasleyC3 that our
        // rules never see — 2-nnz 'L' rows `x_k − 212·b_k ≤ 0` where b_k
        // is a column singleton with nonzero cost. HiGHS's route is its
        // dual propagated `implRowDualLower/Upper` feeding
        // `isImpliedEquationAtLower/Upper` (HPresolve.cpp:263–281, call
        // site :3567–3588, MIT). We have no dual propagation, but for a column SINGLETON the same row-dual bound
        // comes from the reduced-cost sign condition alone — see
        // `singleton_cost_eq_conversions`. Runs before the DTE block so
        // the fresh equations are substituted in the same pass. Lever
        // RUSOLVE_PS_NO_VUBEQ=1 disables.
        if vub_eq {
            n_vub_eq += singleton_cost_eq_conversions(&mut work);
        }

        // ---------- 2-nonzero equality rows: exact substitution (chain killer) ----------
        // Eliminating whole equality CHAINS pays (doubleton-equation
        // substitution, Andersen & Andersen 1995); this rule does exactly
        // that for two-term E rows. x_k = (b − a_m x_m)/a_k is
        // substituted everywhere; x_k's bounds propagate onto x_m so the reduced
        // LP stays exactly equivalent; postsolve replays the records in reverse.
        //
        // The subst/stay choice and the stay-col implied bounds
        // are ports from HiGHS HPresolve.cpp `doubletonEq` (MIT License,
        // HiGHS/ERGO Institute). Anchors: `colAtPos1Better`
        // :3052–3101 and the stay-col bounds :3141–3183. Only the LP branch is
        // ported: every `kInteger` arm (:3053–3071) and the integer
        // divisibility check (:3123–3134) are skipped because rusolve-lp is a
        // pure LP with no integrality domain.
        for i in 0..work.m {
            if !work.ralive[i] || work.senses[i] != 'E' {
                continue;
            }
            let nz: Vec<(usize, f64)> = work.row(i).take(3).collect();
            if nz.len() != 2 {
                continue;
            }
            let ((p, ap), (q, aq)) = (nz[0], nz[1]);
            let b = work.rhs[i];
            // Subst/stay choice, HiGHS colAtPos1Better :3052–3101 (LP branch).
            // `cnt(j)` counts rows containing col j EXCLUDING row i; HiGHS'
            // colsize includes row i, so its `colsize == 1` is our `cnt == 0`
            // and the size comparisons are identical. The subst column is the
            // one eliminated; preference order per the C++:
            //  (a) :3081–3086 singleton first — a column present nowhere else
            //      costs zero fill-in to substitute, regardless of coefficients;
            //  (b) :3090–3092 sizes differ AND |coef| ratio within 2× ⇒ the
            //      sparser column;
            //  (c) :3093–3096 otherwise the LARGER |coef| — the substitution
            //      multiplier is −a_m/a_k, so a bigger |a_k| keeps |mult| ≤ 1
            //      and damps magnitude growth in the surviving rows.
            let cnt = |j: usize| work.col(j).filter(|&(r, _)| r != i).count();
            let (cp, cq) = (cnt(p), cnt(q));
            let (ap_abs, aq_abs) = (ap.abs(), aq.abs());
            let p_subst = if cp == 0 {
                true
            } else if cq == 0 {
                false
            } else if cp != cq && ap_abs.max(aq_abs) <= 2.0 * ap_abs.min(aq_abs) {
                cp < cq
            } else {
                ap_abs > aq_abs
            };
            let (k, m, ak, am) = if p_subst {
                (p, q, ap, aq)
            } else {
                (q, p, aq, ap)
            };
            // Implied bounds on the STAY column x_m = (b − a_k·x_k)/a_m from
            // the subst column's bounds (doubletonEq :3141–3174). Sign split:
            // −a_k/a_m > 0 (coefficients of opposite sign) ⇒ the stay lower
            // comes from the subst LOWER and the stay upper from the subst
            // UPPER; same sign flips the pairing. Quotients hitting an
            // infinite subst bound evaluate to the correct ±∞ under IEEE
            // (HiGHS guards ±kHighsInf explicitly at :3153–3173; the infinite
            // arm of the guard is what `b − a_k·(±∞)` computes here).
            let (imp_lo, imp_hi) = if (ak > 0.0) != (am > 0.0) {
                ((b - ak * work.lower[k]) / am, (b - ak * work.upper[k]) / am)
            } else {
                ((b - ak * work.upper[k]) / am, (b - ak * work.lower[k]) / am)
            };
            if imp_lo > work.upper[m] + TOL || imp_hi < work.lower[m] - TOL {
                traced(
                    "infeasible",
                    format!(
                        "chain row {i}: implied [{imp_lo:?}, {imp_hi:?}] vs x{m} [{:?}, {:?}]",
                        work.lower[m], work.upper[m]
                    ),
                );
                return terminal(SimplexStatus::Infeasible);
            }
            // HiGHS gates the two clamps by an ε-guarded tightened flag and
            // stores lowerTightened/upperTightened in the postsolve record
            // (:3177–3183) ONLY so postsolve can RELAX the stay bounds back
            // to their old values (its bounds move both ways through the
            // MIP search). Ours never relax — bounds only tighten inside
            // `run`, and primal reconstruction restores values, not bounds
            // (we emit no duals, primal-only posture) — so the plain clamps below are
            // safe at our TOL=1e-9 discipline: a clamp past the old bound is
            // justified by the exact row identity, and a ≤TOL bound
            // inversion it can leave is caught by the inversion detection
            // block later in this pass.
            work.lower[m] = work.lower[m].max(imp_lo);
            work.upper[m] = work.upper[m].min(imp_hi);
            // Substitute into objective (only x_m's cost changes) and other rows.
            let ck = work.objective[k];
            if ck != 0.0 {
                // Same cancellation cleanup as the aggregator (seed 836707:
                // a -3x1 / x0 chain pair left c_x0 = 1e-16 and the empty
                // column rule declared Unbounded).
                let (old, delta) = (work.objective[m], ck * am / ak);
                let new = old - delta;
                work.objective[m] = if new.abs() <= SMALL_MATRIX * old.abs().max(delta.abs()) {
                    0.0
                } else {
                    new
                };
                obj_const += ck * b / ak;
            }
            let colk: Vec<(usize, f64)> = work.col(k).collect();
            for (r, akr) in colk {
                if r != i {
                    work.sub_entry(r, m, akr * am / ak);
                    work.rhs[r] -= akr * b / ak;
                }
            }
            recs.push(PostRec::Chain(ChainRec {
                orig: k,
                partner: m,
                ak,
                am,
                b,
                orig_sense: osenses[i],
            }));
            // Remove column k, then row i.
            work.kill_col(k);
            work.kill_row(i);
        }

        // ---------- rows: infeasibility / empties / redundancy / forcing ----------
        let mut keep = vec![true; work.m];
        let mut fixes: Vec<(usize, f64)> = Vec::new(); // col, value
        for i in 0..work.m {
            if !work.ralive[i] {
                continue;
            }
            let (mut amin, mut amax) = (0.0f64, 0.0f64);
            let mut lo_at: Vec<(usize, f64)> = Vec::new();
            let mut hi_at: Vec<(usize, f64)> = Vec::new();
            let mut nnz = 0usize;
            let (mut inf_lo, mut inf_hi) = (false, false);
            for (j, a) in work.row(i) {
                nnz += 1;
                let (lb, ub) = (work.lower[j], work.upper[j]);
                if a > 0.0 {
                    if lb.is_finite() {
                        amin += a * lb;
                        lo_at.push((j, lb));
                    } else {
                        inf_lo = true;
                    }
                    if ub.is_finite() {
                        amax += a * ub;
                        hi_at.push((j, ub));
                    } else {
                        inf_hi = true;
                    }
                } else {
                    if ub.is_finite() {
                        amin += a * ub;
                        lo_at.push((j, ub));
                    } else {
                        inf_lo = true;
                    }
                    if lb.is_finite() {
                        amax += a * lb;
                        hi_at.push((j, lb));
                    } else {
                        inf_hi = true;
                    }
                }
            }
            let s = work.senses[i];
            let r = work.rhs[i];
            let infeasible = match s {
                'L' => !inf_lo && amin > r + TOL * r.abs().max(1.0),
                'G' => !inf_hi && amax < r - TOL * r.abs().max(1.0),
                'E' => {
                    (!inf_lo && amin > r + TOL * r.abs().max(1.0))
                        || (!inf_hi && amax < r - TOL * r.abs().max(1.0))
                }
                _ => false,
            };
            if infeasible {
                traced(
                    "infeasible",
                    format!(
                        "row {i} sense {s} rhs {r:?} activity [{amin:?}, {amax:?}] inf=({inf_lo},{inf_hi})"
                    ),
                );
                return terminal(SimplexStatus::Infeasible);
            }
            if nnz == 0 {
                keep[i] = false; // provably satisfied (the infeasibility check
                // above already caught a violated empty row)
                continue;
            }
            let redundant = match s {
                'L' => !inf_hi && amax <= r + TOL,
                'G' => !inf_lo && amin >= r - TOL,
                _ => false,
            };
            if redundant {
                keep[i] = false;
                continue;
            }
            // Forcing tolerance is relative-to-rhs (see at_rhs): an absolute
            // tol false-forces on large-coefficient rows (boeing2 evidence).
            let (force_hi, force_lo) = match s {
                'L' => (!inf_hi && at_rhs(r, amax), false),
                'G' => (false, !inf_lo && at_rhs(r, amin)),
                'E' => (!inf_hi && at_rhs(r, amax), !inf_lo && at_rhs(r, amin)),
                _ => (false, false),
            };
            if force_hi {
                fixes.append(&mut hi_at);
                keep[i] = false;
            } else if force_lo {
                fixes.append(&mut lo_at);
                keep[i] = false;
            }
        }
        work.drop_rows(&keep);
        for (j, v) in fixes {
            // Forcing-conflict check: two forcing rows scanned in ONE
            // pass pin the same column from the same bounds snapshot — each
            // verdict is a valid deduction, so a disagreement beyond TOL
            // (x pinned to its lower bound by one row and its upper bound by
            // another) PROVES the LP infeasible. The old code clamped the
            // later value into the earlier fix's box, silently swallowing the
            // conflict while BOTH rows had already been dropped — a set-
            // partitioning LP reduced to 0 rows and reported Optimal.
            if v < work.lower[j] - TOL || v > work.upper[j] + TOL {
                traced(
                    "infeasible",
                    format!(
                        "forcing conflict on col {j}: pinned {v:?} vs bounds [{:?}, {:?}]",
                        work.lower[j], work.upper[j]
                    ),
                );
                return terminal(SimplexStatus::Infeasible);
            }
            traced("forcing", format!("row fixes col {j} := {v:?}"));
            let v = v.max(work.lower[j]).min(work.upper[j]);
            work.lower[j] = v;
            work.upper[j] = v;
        }

        // ---------- tightening: single-nonzero rows + interval rows ----------
        for i in 0..work.m {
            if !work.ralive[i] {
                continue;
            }
            let entries: Vec<(usize, f64)> = work.row(i).collect();
            let mut nz: Option<(usize, f64)> = None;
            let nnz = entries.len();
            let mut amin = 0.0;
            let mut amax = 0.0;
            let mut inf = false;
            for &(j, a) in &entries {
                nz = Some((j, a));
                let (p, q) = if a > 0.0 {
                    (work.lower[j], work.upper[j])
                } else {
                    (work.upper[j], work.lower[j])
                };
                if p.is_finite() {
                    amin += a * p;
                } else {
                    inf = true;
                }
                if q.is_finite() {
                    amax += a * q;
                } else {
                    inf = true;
                }
            }
            if nnz == 0 {
                continue;
            }
            let (r, s) = (work.rhs[i], work.senses[i]);
            if nnz == 1 {
                let (j, a) = nz.unwrap();
                match s {
                    'E' => {
                        let v = r / a;
                        if v < work.lower[j] - TOL || v > work.upper[j] + TOL {
                            traced(
                                "infeasible",
                                format!(
                                    "single-nonzero E row {i}: x{j} = {v:?} outside [{:?}, {:?}]",
                                    work.lower[j], work.upper[j]
                                ),
                            );
                            return terminal(SimplexStatus::Infeasible);
                        }
                        let v = v.max(work.lower[j]).min(work.upper[j]); // not clamp(): lb may exceed ub by rounding noise
                        traced(
                            "single-E",
                            format!("row {i} fixes x{j} := {v:?} (rhs {r:?} / a {a:?})"),
                        );
                        work.lower[j] = v;
                        work.upper[j] = v;
                    }
                    'L' if a > 0.0 => work.tighten_ub(j, r / a),
                    'L' => work.tighten_lb(j, r / a),
                    'G' if a > 0.0 => work.tighten_lb(j, r / a),
                    'G' => work.tighten_ub(j, r / a),
                    _ => {}
                }
                continue;
            }
            if inf {
                continue;
            }
            for &(j, a) in &entries {
                let (p, q) = if a > 0.0 {
                    (work.lower[j], work.upper[j])
                } else {
                    (work.upper[j], work.lower[j])
                };
                let others_min = amin - a * p;
                let others_max = amax - a * q;
                match s {
                    'L' => {
                        let cap = (r - others_min) / a;
                        if a > 0.0 {
                            work.prop_ub(j, cap);
                        } else {
                            work.prop_lb(j, cap);
                        }
                    }
                    'G' => {
                        let cap = (r - others_max) / a;
                        if a > 0.0 {
                            work.prop_lb(j, cap);
                        } else {
                            work.prop_ub(j, cap);
                        }
                    }
                    'E' => {
                        let cap_from_min = (r - others_min) / a;
                        let cap_from_max = (r - others_max) / a;
                        if a > 0.0 {
                            work.prop_ub(j, cap_from_min);
                            work.prop_lb(j, cap_from_max);
                        } else {
                            work.prop_lb(j, cap_from_min);
                            work.prop_ub(j, cap_from_max);
                        }
                    }
                    _ => {}
                }
            }
        }

        // ---------- aggregator: implied-free column substitution (HiGHS :6688) ----------
        // Runs on freshly tightened bounds; see `aggregate`.
        if agg {
            n_agg += aggregate(&mut work, &mut obj_const, &mut recs);
        }

        // ---------- dominated-column detection probe (detection only) ----------
        // Ported detection logic from HiGHS HPresolve.cpp detectDominatedCol
        // :4514–4551 (MIT), sense-only row dual seeds = impliedDualRowBounds
        // Orig sums; soundness: r_j strictly positive/negative over the dual
        // box ⇒ some optimum at that bound. Runs ONLY when the stats flag is
        // set and mutates nothing, so default behavior stays byte-identical
        // (off-probe). TODO: propagate row dual bounds via
        // strictly-implied col bounds.
        if stats {
            let (cands, near, dual_inf) = dominated_col_candidates(&work, DUAL_TOL);
            for &(_j, dir) in &cands {
                if dir == 1 {
                    n_dom_lb += 1;
                } else {
                    n_dom_ub += 1;
                }
            }
            n_near_miss += near;
            n_dual_inf += dual_inf;
        }

        // ---------- forcing-column removal (mutation + postsolve) ----------
        // Ported from HiGHS (MIT):
        // detection = weaklyDominatedCol forcing fallback, HPresolve.cpp
        // :4553–4617 (the gate and the ±∞ seed sums are the detection probe's,
        // unchanged); removal semantics = removeRow :4611 (the column's
        // rows leave presolve entirely); reconstruction = HighsPostsolveStack
        // .cpp ForcingColumn::undo :275+ (far-side-overlap value, primal
        // only). Algorithm source: Andersen & Andersen, "Presolving in
        // linear programming", Math. Program. 71, 221–245 (1995).
        // Soundness: the ≈0 gate over SENSE-ONLY
        // dual boxes forces cost_c ≈ 0 AND every row of c to bound x_c from
        // the finite side (L a>0 / G a<0 for dir +1, mirrored for dir −1),
        // so deleting column + rows loses nothing: any reduced solution
        // extends to a feasible full one at x_c = min/max(finite bound, the
        // row bounds), with objective drift ≤ tol/dynamism-scale (HiGHS
        // accepts the same). The propagated-box variant is NOT ported here
        // (needs weak-gate dual propagation, which we do not have).
        // One candidate per pass (the first with a finite counterpart bound
        // — a fully free col's "finite" side is ∞: counted, not mutated);
        // the fixpoint loop re-scans. NOT the `fixed` path: c is never
        // value-known at presolve time, so it rides the record vector.
        if stats || !no_forcecol {
            let cands = forcing_col_candidates(&work, DUAL_TOL);
            if stats {
                n_force_lb += cands.iter().filter(|&&(_, d)| d == 1).count();
                n_force_ub += cands.iter().filter(|&&(_, d)| d == -1).count();
            }
            if !no_forcecol
                && let Some(&(j, dir)) = cands.iter().find(|&&(j, d)| {
                    (if d == 1 { work.upper[j] } else { work.lower[j] }).is_finite()
                })
            {
                let finite_bound = if dir == 1 {
                    work.upper[j]
                } else {
                    work.lower[j]
                };
                let mut rows: Vec<ForcedRow> = Vec::new();
                let mut keep = vec![true; work.m];
                let colj: Vec<(usize, f64)> = work.col(j).collect();
                for (i, _) in colj {
                    rows.push(ForcedRow {
                        row: i,
                        sense: work.senses[i],
                        rhs: work.rhs[i],
                        entries: work.row(i).collect(),
                    });
                    keep[i] = false;
                }
                traced(
                    "forcing-col",
                    format!(
                        "orig col {j} dir {dir} finite bound {finite_bound:?}: {} rows deleted",
                        rows.len()
                    ),
                );
                recs.push(PostRec::Forcing(ForcingRec {
                    orig: j,
                    dir,
                    finite_bound,
                    rows,
                }));
                // Delete the column (no fold: value unknown; its cost is
                // ≈ 0 by the gate, the discarded drift is the A&A price)
                // then ALL its rows — later rules never see them.
                work.kill_col(j);
                work.drop_rows(&keep);
            }
        }

        // ---------- parallel-column dominance fixing (mutation) ----------
        // The LP dominance arms of HiGHS detectParallelRowsAndCols
        // (:7399–7540 + fixCol :7601–7617, MIT; see
        // `parallel_col_dominance_fixes`). Proportional columns with
        // unequal costs (or unmergeable boxes) fix one side at a bound —
        // the compensation gate is the strictly-row-implied bound of the
        // surviving column. Value-known fixes ride the `fixed`
        // substitution at the end of the pass; no postsolve record.
        // Census: HiGHS's parallel
        // rule removes 5,659 cols on neos-1445765 and 543 on fiball that
        // our merge-only parcol never sees. Lever
        // RUSOLVE_PS_NO_PARCOLDOM=1 disables the mutation.
        if parcol_dom {
            let nfix = parallel_col_dominance_fixes(&mut work).len();
            n_parcol_dom += nfix;
        }

        // ---------- parallel-row merge (mutation) + detection probe ----------
        // Ported from HiGHS (MIT):
        // detection = `detectParallelRowsAndCols` :7267–7450 (scale-normalized
        // bucketing + proportionality re-check); dominance = the LP branch of
        // its row reduction (keep the row whose scaled bound dominates, delete
        // the parallel losers). Lever RUSOLVE_PS_NO_PARROW=1 disables the
        // mutation, the stats probe keeps counting.
        // NO POSTSOLVE RECORD — verified from the original: HighsPostsolveStack
        // ::DuplicateRow::undo (HighsPostsolveStack.cpp:505–563) writes ONLY
        // row_dual / row_status, never primal values: rows carry no primal
        // component, and a deleted loser is implied by the kept dominant row
        // (same normalized coefficients, weaker scaled rhs), so any x feasible
        // for the kept row satisfies it. Primal-only postsolve (no duals
        // emitted) therefore has nothing to reconstruct.
        // The cross-sense bound intersection HiGHS performs on its
        // interval rows (row_lower/row_upper per row) IS ported now — see the
        // effective-interval arithmetic below. The COLLAPSE arm (interval
        // reduces to a point ⇒ keep one row converted to 'E') is the graph20
        // class, lever RUSOLVE_PS_NO_PARROWEQ=1; a true range keeps one
        // witness per side and deletes only dominated same-side losers.
        if stats || !no_parrow {
            let groups = parallel_row_groups(&work);
            if stats {
                let (rg, rm) = (groups.len(), groups.iter().map(Vec::len).sum::<usize>());
                if rg > par_row_groups {
                    par_row_groups = rg;
                    par_row_rows = rm;
                }
            }
            if !no_parrow {
                let mut keep = vec![true; work.m];
                for grp in groups {
                    // Effective rhs AND effective direction under the shared
                    // normalized coefficients v: row i reads `s_i·(v·x) (sense)
                    // rhs_i` with s_i = signed max|coef|. Dividing by a
                    // NEGATIVE s_i flips the sense, so an 'L' row can bound
                    // v·x from BELOW (e226 rows 49/50: −0.29(x+y) ≤ −0.2662
                    // ⇔ x+y ≥ 0.9179, which DOMINATES x+y ≥ 0.3467 — the
                    // first cut of this rule deleted the wrong one and made
                    // the e226 postsolve x infeasible by 5.8e-2).
                    // The group is the conjunction  v·x ∈ [max lower eff,
                    // min upper eff] ∩ {equalities' eff} — ported from
                    // HiGHS (MIT): the
                    // interval intersection and rowUpperTightened/
                    // rowLowerTightened bookkeeping of
                    // `detectParallelRowsAndCols`, HPresolve.cpp:7819–7893
                    // (duplicateRow keeps one row with the tightened
                    // interval). Expressible arms in
                    // the sense-only contract:
                    // * same direction: min/max eff wins;
                    // * equalities: all eff agree within TOL and lie inside
                    //   the sided interval (within the TOL sliver, the same
                    //   slack-scale posture as the conversion tier) ⇒ keep
                    //   one 'E', delete the rest; disagreement or an eff
                    //   outside the interval ⇒ SKIP (infeasibility is proven
                    //   elsewhere, never silently cut);
                    // * mixed sides whose interval COLLAPSES to a point
                    //   (|u−d| ≤ TOL rel) ⇒ the conjunction IS the equality
                    //   v·x = d; convert the lower witness to 'E' and delete
                    //   everything else (graph20: 185 'G' rhs 0 ∧ 'L' rhs 0
                    //   pairs ⇒ equality ⇒ DTE cascade);
                    // * a true range (d < u beyond TOL) is inexpressible in
                    //   the sense-only contract: keep ONE witness per side,
                    //   delete dominated same-side losers only.
                    let scales: Vec<f64> =
                        grp.iter().map(|&i| scale_normal(work.row(i)).1).collect();
                    let eff: Vec<f64> = grp
                        .iter()
                        .zip(&scales)
                        .map(|(&i, &s)| work.rhs[i] / s)
                        .collect();
                    // +1 ⇒ upper bound on v·x, −1 ⇒ lower bound, 0 ⇒ equality.
                    let dir: Vec<i8> = grp
                        .iter()
                        .zip(&scales)
                        .map(|(&i, &s)| match work.senses[i] {
                            'E' => 0,
                            'L' => {
                                if s > 0.0 {
                                    1
                                } else {
                                    -1
                                }
                            }
                            'G' => {
                                if s > 0.0 {
                                    -1
                                } else {
                                    1
                                }
                            }
                            _ => i8::MIN, // 'N': never grouped
                        })
                        .collect();
                    if dir.contains(&i8::MIN) {
                        continue;
                    }
                    let tol_at = |v: f64| TOL * v.abs().max(1.0);
                    let ubest = (0..grp.len())
                        .filter(|&k| dir[k] == 1)
                        .min_by(|a, b| eff[*a].partial_cmp(&eff[*b]).unwrap());
                    let dbest = (0..grp.len())
                        .filter(|&k| dir[k] == -1)
                        .max_by(|a, b| eff[*a].partial_cmp(&eff[*b]).unwrap());
                    let uval = ubest.map_or(f64::INFINITY, |k| eff[k]);
                    let dval = dbest.map_or(f64::NEG_INFINITY, |k| eff[k]);
                    // (winner index k, extra: convert winner to 'E')
                    let winner: Option<(usize, bool)> =
                        if let Some(e0) = dir.iter().position(|&d| d == 0) {
                            let e = eff[e0];
                            let agree = dir
                                .iter()
                                .enumerate()
                                .all(|(k, &d)| d != 0 || (eff[k] - e).abs() <= tol_at(e));
                            let inside = dval - tol_at(e) <= e && e <= uval + tol_at(e);
                            if agree && inside {
                                Some((e0, false))
                            } else {
                                None // conflicting E's / inconsistent system
                            }
                        } else {
                            match (ubest, dbest) {
                                (Some(u), Some(d)) => {
                                    if parrow_eq && (uval - dval).abs() <= tol_at(uval.max(dval)) {
                                        // Interval collapses to a point: the
                                        // lower witness becomes the equality.
                                        Some((d, true))
                                    } else {
                                        // True range: keep both witnesses
                                        // (encoded via winner u + keep_d).
                                        let _ = u;
                                        None
                                    }
                                }
                                (Some(u), None) => Some((u, false)),
                                (None, Some(d)) => Some((d, false)),
                                (None, None) => None,
                            }
                        };
                    // True-range groups keep one witness per side; every
                    // other group keeps a single winner (or nothing).
                    let keeps: Vec<usize> = match winner {
                        Some((w, _)) => vec![w],
                        None if ubest.is_some()
                            && dbest.is_some()
                            && dir.iter().all(|&d| d != 0) =>
                        {
                            vec![dbest.unwrap(), ubest.unwrap()]
                        }
                        _ => continue,
                    };
                    let convert = winner.is_some_and(|(_, c)| c);
                    if convert {
                        let i = grp[keeps[0]];
                        traced(
                            "par-row-eq",
                            format!(
                                "orig row {} interval collapse: eff [{dval:?}, {uval:?}] ⇒ 'E' (was {:?})",
                                i, work.senses[i]
                            ),
                        );
                        work.senses[i] = 'E';
                        n_parrow_eq += 1;
                    }
                    for (k, &i) in grp.iter().enumerate() {
                        if keeps.contains(&k) {
                            continue;
                        }
                        keep[i] = false;
                        n_parrow += 1;
                        if debug {
                            let nz_i: Vec<(usize, f64)> = work.row(i).collect();
                            let nz_w: Vec<(usize, f64)> = work.row(grp[keeps[0]]).collect();
                            eprintln!(
                                "[presolve] par-row-full: DEL orig {} sense {:?} rhs {:?} scale {:?} nz {:?} | KEEP orig {} sense {:?} rhs {:?} scale {:?} nz {:?}",
                                i,
                                work.senses[i],
                                work.rhs[i],
                                scale_normal(work.row(i)).1,
                                nz_i,
                                grp[keeps[0]],
                                work.senses[grp[keeps[0]]],
                                work.rhs[grp[keeps[0]]],
                                scale_normal(work.row(grp[keeps[0]])).1,
                                nz_w
                            );
                        }
                        traced(
                            "par-row",
                            format!(
                                "orig row {} (eff rhs {:?}) implied by kept orig row(s) {:?} (effs {:?})",
                                i,
                                eff[k],
                                keeps.iter().map(|&w| grp[w]).collect::<Vec<_>>(),
                                keeps.iter().map(|&w| eff[w]).collect::<Vec<_>>(),
                            ),
                        );
                    }
                }
                work.drop_rows(&keep);
            }
        }

        // ---------- duplicate-column merge (mutation) + probe ----------
        // Ported from HiGHS (MIT):
        // detection = `parallel_col_groups` (detectParallelRowsAndCols cols
        // half, cost class = the kMergeParallelCols objDiff branch); merge =
        // postsolve_stack.duplicateColumn (HighsPostsolveStack.h:513) + the
        // merge-bound arithmetic at HPresolve.cpp:7674-7684, GENERAL scale
        // (the ±1 snap was the unported integer-arm machinery —
        // fractionality/divisibility, which pure LP does not need; the kDominanceCol*
        // fixings are the separate dominance rule). The duplicate's coefficients
        // are scale·kept on
        // EVERY row, so the merge touches no row entries: widen the kept
        // column to the merged box, drop the duplicate, record the split
        // (DupColRec). Default OFF per bench (run's doc); `parcol` is the
        // gate, `RUSOLVE_PS_PARCOL=1` the CLI lever.
        if stats || parcol {
            let groups = parallel_col_groups(&work);
            if stats {
                let (cg, cm) = (groups.len(), groups.iter().map(Vec::len).sum::<usize>());
                if cg > par_col_groups {
                    par_col_groups = cg;
                    par_col_cols = cm;
                }
            }
            if parcol {
                // Collect (rep, dup, scale) ops across ALL groups, then
                // apply in descending dup order (kept from the dense
                // implementation, where index shifting required it).
                let mut merges: Vec<(usize, usize, f64)> = Vec::new();
                for grp in groups {
                    let rep = grp[0];
                    let s_rep = scale_normal(work.col(rep)).1;
                    for &d in &grp[1..] {
                        let s = scale_normal(work.col(d)).1 / s_rep;
                        // General scale: HiGHS `duplicateColumn`
                        // (HighsPostsolveStack.h:513 via HPresolve.cpp:7634)
                        // accepts ANY real colScale in LP — the ±1 snap was the
                        // integer-arm machinery (fractionality/divisibility),
                        // which this pure-LP presolve never had. Guard against
                        // an overflowed ratio poisoning the bounds.
                        if !s.is_finite() {
                            continue;
                        }
                        merges.push((rep, d, s));
                    }
                }
                merges.sort_by_key(|&(_, d, _)| std::cmp::Reverse(d));
                for (rep, d, s) in merges {
                    // Merged box m = x_rep + s·x_d (HiGHS merge-bounds
                    // HPresolve.cpp:7674-7684, MIT: side picked by sign of the
                    // scale — for s<0 the duplicate's UPPER feeds the merged
                    // LOWER). IEEE infinity arithmetic is correct for every
                    // reachable operand (a fixed-at-infinity column is
                    // excluded upstream); a NaN could only come from ∓inf +
                    // ±inf — guard and skip conservatively.
                    let (m_lo, m_hi) = if s > 0.0 {
                        (
                            work.lower[rep] + s * work.lower[d],
                            work.upper[rep] + s * work.upper[d],
                        )
                    } else {
                        (
                            work.lower[rep] + s * work.upper[d],
                            work.upper[rep] + s * work.lower[d],
                        )
                    };
                    if m_lo.is_nan() || m_hi.is_nan() {
                        continue;
                    }
                    traced(
                        "par-col",
                        format!(
                            "orig col {d} merged into {rep} (scale {s:?}): merged box [{m_lo:?}, {m_hi:?}]"
                        ),
                    );
                    recs.push(PostRec::DupCol(DupColRec {
                        orig: d,
                        kept: rep,
                        scale: s,
                        kept_lo: work.lower[rep],
                        kept_hi: work.upper[rep],
                        dup_lo: work.lower[d],
                        dup_hi: work.upper[d],
                    }));
                    work.lower[rep] = m_lo;
                    work.upper[rep] = m_hi;
                    // Objective: the kept column's cost IS the merged cost
                    // (objDiff gate ⇒ dup's cost = s·c_rep within DUAL_TOL;
                    // the dropped drift is DupColRec's doc).
                    work.kill_col(d);
                    n_parcol += 1;
                }
            }
        }

        // ---------- implied-equality conversion + far-side forcing (HiGHS port) ----------
        // Ported from HiGHS HPresolve.cpp (MIT). Anchors: rowPresolve's
        // conversion block :3567–3588 and isImpliedEquationAtLower/Upper
        // :263–281, plus checkForcingRow (:4349–4368), which is what HiGHS
        // runs FIRST on exactly this condition. HiGHS' :3567 trigger is a
        // *dual* bound (complementary slackness: implied row-dual lower >
        // dual_feasibility_tolerance ⇒ row tight at its lower side); the
        // primal-activity instance ported here is the same inference an LP
        // without dual-bound propagation can draw soundly: if the implied
        // activity range (from current bounds) touches the far side of a
        // ranged row, every feasible point has activity == rhs. Two tiers:
        // * PINNED within the strict forcing tolerance at_rhs (1e-12 rel):
        //   the row is a forcing row — every variable is pinned at the
        //   bound achieving the extreme, so fix them and delete the row
        //   (HiGHS checkForcingRow). Deleting matters: converted-but-kept
        //   forcing rows leave their full degeneracy in the simplex (scrs8
        //   phase-2 312→1830 pivots when six zero-pinned rows were only
        //   converted), whereas forcing removes it.
        // * Within the looser TOL·max(|rhs|,1) slack (1e-9 rel): fixing on
        //   that argument is what boeing2 false-forced on (see at_rhs doc),
        //   but converting to equality is safe — it can cut off at most a
        //   slack-scale sliver of feasible activity. The row becomes an
        //   equality for the rest of presolve, feeding the doubleton chain
        //   above; the pre-conversion sense is remembered in ChainRec
        //   (HiGHS does the same via HighsPostsolveStack::RowType,
        //   :3596–3607).
        // Tolerances are strictly tighter than HiGHS' 1e-7 dual gate
        // (HConst.h:344-345). MONOTONE (hence no fixpoint oscillation):
        // only L/G → E plus row deletion, senses never move back, E skipped.
        let mut keep = vec![true; work.m];
        // Forcing completions found by the far-side check: (col, value)
        // pairs, applied after drop_rows exactly like the main loop's
        // forcing fixes.
        let mut eq_fixes: Vec<(usize, f64)> = Vec::new();
        for i in 0..work.m {
            if !work.ralive[i] {
                continue;
            }
            let s = work.senses[i];
            if !implied_eq || s == 'E' {
                continue; // gate off / already converted / always an equation
            }
            let (mut amin, mut amax) = (0.0f64, 0.0f64);
            let mut lo_at: Vec<(usize, f64)> = Vec::new();
            let mut hi_at: Vec<(usize, f64)> = Vec::new();
            let (mut inf_lo, mut inf_hi) = (false, false);
            for (j, a) in work.row(i) {
                let (lb, ub) = (work.lower[j], work.upper[j]);
                if a > 0.0 {
                    if lb.is_finite() {
                        amin += a * lb;
                        lo_at.push((j, lb));
                    } else {
                        inf_lo = true;
                    }
                    if ub.is_finite() {
                        amax += a * ub;
                        hi_at.push((j, ub));
                    } else {
                        inf_hi = true;
                    }
                } else {
                    if ub.is_finite() {
                        amin += a * ub;
                        lo_at.push((j, ub));
                    } else {
                        inf_lo = true;
                    }
                    if lb.is_finite() {
                        amax += a * lb;
                        hi_at.push((j, lb));
                    } else {
                        inf_hi = true;
                    }
                }
            }
            let r = work.rhs[i];
            let slack = TOL * r.abs().max(1.0);
            // The far side of a 'L' row is the activity minimum, of a 'G'
            // row the maximum; touching it pins activity at rhs.
            let pinned = match s {
                'L' => !inf_lo && amin >= r - slack,
                'G' => !inf_hi && amax <= r + slack,
                _ => false,
            };
            if !pinned {
                continue;
            }
            let strict = match s {
                'L' => at_rhs(r, amin),
                _ => at_rhs(r, amax),
            };
            if strict {
                traced(
                    "implied-eq-force",
                    format!(
                        "row {i} '{s}' rhs {r:?}: activity [{amin:?}, {amax:?}] pinned ⇒ fix + drop"
                    ),
                );
                eq_fixes.append(if s == 'L' { &mut lo_at } else { &mut hi_at });
                keep[i] = false;
                n_forced += 1;
            } else {
                traced(
                    "implied-eq",
                    format!(
                        "row {i} '{s}' rhs {r:?}: activity [{amin:?}, {amax:?}] within slack of rhs ⇒ equality"
                    ),
                );
                work.senses[i] = 'E';
                n_eq_converted += 1;
            }
        }

        // ---------- pure-bound row fold: extractVarBounds (HiGHS port) ----------
        // Ported from HiGHS HPresolve.cpp (MIT). Anchors: extractVarBounds
        // :8039–8137 (call site rowPresolve :4375) and its LP-side sibling
        // removeRowSingletons :6945–6958 → rowSingleton (fold the row's
        // bound into the column, delete the row as redundant). The HiGHS
        // extractVarBounds itself is MIP-only (VLB/VUB implications against
        // a lone binary; it returns immediately when mipsolver == nullptr);
        // the LP instantiation ported here is the same inference one level
        // down: if every nonzero but one sits on a FIXED variable
        // (lb == ub), the row constrains exactly one free variable — it is
        // a pure bound. Fold the derived bound into that variable and
        // delete the row: after the fold the row is implied by the bound
        // (bounds only ever tighten from here), so primal postsolve is
        // untouched. This catches rows interval tightening must skip (an
        // infinite bound on the free column poisons its activity sums) and
        // singleton rows the redundancy check would only catch next pass.
        // Tolerances: the fold cap is applied unrelaxed, identical to the
        // interval tightening above (identical arithmetic, identical
        // posture); HiGHS relaxes its implied bounds by
        // boundRelax = max(1000,|b|)·primal_feastol/min(1,|a|) (:8012–8035)
        // only because its sums are incrementally maintained — ours are
        // recomputed each pass from exactly-fixed values, so the noise
        // stays ~1 ulp ≪ TOL=1e-9. Infeasibility of an E fold is judged
        // with the same absolute TOL as the existing single-nonzero E rule.
        // `keep`/`eq_fixes` carry the far-side forcing verdicts from the
        // block above; rows already condemned there must not be re-read.
        for i in 0..work.m {
            if !work.ralive[i] || no_fold || !keep[i] {
                continue; // A/B lever / far-side forcing already dropped this row
            }
            let mut free: Option<(usize, f64)> = None; // the one non-fixed nonzero
            let mut fixed_sum = 0.0f64;
            let mut skip = false;
            for (j, a) in work.row(i) {
                if work.lower[j] == work.upper[j] {
                    if !work.lower[j].is_finite() {
                        skip = true; // fixed at infinity: degenerate, hands off
                        break;
                    }
                    fixed_sum += a * work.lower[j];
                } else if free.is_some() {
                    skip = true; // more than one free column: not a pure bound
                    break;
                } else {
                    free = Some((j, a));
                }
            }
            if skip || free.is_none() || !fixed_sum.is_finite() {
                continue;
            }
            let (j, a) = free.unwrap();
            let r = work.rhs[i] - fixed_sum;
            match work.senses[i] {
                'E' => {
                    let v = r / a;
                    if v < work.lower[j] - TOL || v > work.upper[j] + TOL {
                        traced(
                            "infeasible",
                            format!(
                                "pure-bound E row {i}: x{j} = {v:?} outside [{:?}, {:?}]",
                                work.lower[j], work.upper[j]
                            ),
                        );
                        return terminal(SimplexStatus::Infeasible);
                    }
                    let v = v.max(work.lower[j]).min(work.upper[j]); // not clamp(): lb may exceed ub by rounding noise
                    traced("pure-bound", format!("row {i} fixes x{j} := {v:?}"));
                    work.lower[j] = v;
                    work.upper[j] = v;
                }
                'L' => {
                    if a > 0.0 {
                        work.tighten_ub(j, r / a);
                    } else {
                        work.tighten_lb(j, r / a);
                    }
                    traced("pure-bound", format!("row {i} bounds x{j}, row dropped"));
                }
                'G' => {
                    if a > 0.0 {
                        work.tighten_lb(j, r / a);
                    } else {
                        work.tighten_ub(j, r / a);
                    }
                    traced("pure-bound", format!("row {i} bounds x{j}, row dropped"));
                }
                _ => {}
            }
            keep[i] = false;
            n_folded += 1;
        }
        work.drop_rows(&keep);
        for (j, v) in eq_fixes {
            // Same forcing-conflict check as the main block above: the
            // pure-bound fold may have tightened this column's bounds since
            // the far-side scan; a forcing value outside them beyond TOL is
            // a contradiction of two sound deductions ⇒ infeasible, never a
            // silent clamp.
            if v < work.lower[j] - TOL || v > work.upper[j] + TOL {
                traced(
                    "infeasible",
                    format!(
                        "forcing conflict on col {j}: pinned {v:?} vs bounds [{:?}, {:?}]",
                        work.lower[j], work.upper[j]
                    ),
                );
                return terminal(SimplexStatus::Infeasible);
            }
            let v = v.max(work.lower[j]).min(work.upper[j]);
            work.lower[j] = v;
            work.upper[j] = v;
        }

        // ---------- bound inversion detection (tightening may prove
        // infeasibility by pushing lb past ub) ----------
        for j in 0..work.n {
            if !work.calive[j] {
                continue;
            }
            if work.lower[j] > work.upper[j] + TOL {
                traced(
                    "infeasible",
                    format!(
                        "bound inversion on col {j}: lb {:?} > ub {:?}",
                        work.lower[j], work.upper[j]
                    ),
                );
                return terminal(SimplexStatus::Infeasible);
            }
            if work.lower[j] > work.upper[j] {
                // within tolerance: clamp to a point
                let v = (work.lower[j] + work.upper[j]) / 2.0;
                work.lower[j] = v;
                work.upper[j] = v;
            }
        }

        // ---------- fixed variables: substitute out ----------
        for j in 0..work.n {
            if work.calive[j] && work.lower[j] == work.upper[j] {
                let v = work.lower[j];
                work.substitute(j, v, &mut states, &mut obj_const);
            }
        }

        // ---------- empty columns ----------
        for j in 0..work.n {
            if work.calive[j] && work.col(j).next().is_none() {
                let c = work.objective[j];
                let target = if c == 0.0 {
                    // indifferent objective: any bound-feasible interior value
                    work.lower[j].max(0.0).min(work.upper[j])
                } else {
                    let prefer_ub = if lp.maximize { c > 0.0 } else { c < 0.0 };
                    if prefer_ub {
                        work.upper[j]
                    } else {
                        work.lower[j]
                    }
                };
                if target.is_infinite() {
                    // Improving infinite bound: a CANDIDATE unbounded ray,
                    // but an undetected infeasibility elsewhere outranks it
                    // — fuzz seed 1756 (vubeq generator) had presolve
                    // answer Unbounded on an infeasible LP because this
                    // shortcut won the race against the tightening that
                    // proved infeasibility one block later. The simplex
                    // settles Unbounded vs Infeasible correctly in one
                    // ratio test, so leave the column in the reduced LP
                    // instead of terminating presolve here.
                    continue;
                }
                work.substitute(j, target, &mut states, &mut obj_const);
            }
        }

        if work.signature() == sig_before {
            break;
        }
    }

    if stats {
        // Live-entry nnz of the reduced LP (row() filters dead ids and exact
        // zeros — the entries the reduced solve actually sees).
        let live_nnz: usize = (0..work.m)
            .filter(|&i| work.ralive[i])
            .map(|i| work.row(i).count())
            .sum();
        eprintln!(
            "[presolve] rows {} -> {}, cols {} -> {}, nnz {} -> {} (aggregated cols: {n_agg}, pure-bound folds: {n_folded}, implied equalities: {n_eq_converted}, far-side forcings: {n_forced}, dominated-lb: {n_dom_lb}, dominated-ub: {n_dom_ub}, near-miss: {n_near_miss}, dual-infeasible: {n_dual_inf}, forcing-lb: {n_force_lb}, forcing-ub: {n_force_ub}, par-rows: {par_row_groups} ({par_row_rows} rows, {n_parrow} merged, {n_parrow_eq} eq-collapsed), par-cols: {par_col_groups} ({par_col_cols} cols, {n_parcol} merged), par-col-dom fixes: {n_parcol_dom}, singleton-eq: {n_vub_eq})",
            lp.num_row,
            work.nrow,
            lp.num_col,
            work.ncol,
            lp.a.nnz(),
            live_nnz
        );
    }
    let fixed: Vec<Option<f64>> = states
        .iter()
        .map(|s| match s {
            VarState::Fixed(v) => Some(*v),
            VarState::Pending => None,
        })
        .collect();
    let (reduced, live) = work.to_lp();
    if stats {
        eprintln!("[presolve-time] {:.3}s", t_start.elapsed().as_secs_f64());
    }
    hash_probe(&reduced, &live, obj_const, recs.len());
    Presolved::Reduced {
        lp: reduced,
        live,
        postsolve: Postsolve {
            fixed,
            obj_const,
            recs,
        },
    }
}

/// Restore the full-space solution from the reduced solve.
pub fn restore(
    postsolve: &Postsolve,
    live: &[usize],
    reduced_x: &[f64],
    reduced_obj: f64,
) -> (Vec<f64>, f64) {
    let mut x = vec![0.0; postsolve.fixed.len()];
    for (orig, f) in postsolve.fixed.iter().enumerate() {
        if let Some(v) = f {
            x[orig] = *v;
        }
    }
    for (col, &orig) in live.iter().enumerate() {
        x[orig] = reduced_x.get(col).copied().unwrap_or(0.0);
    }
    // Replay eliminations newest-first over the MERGED record sequence: each
    // record's dependencies (chain partner / forcing row entries) are either
    // already resolved (live/fixed) or eliminated by a record created later,
    // which this reverse order visits first — regardless of record kind.
    for rec in postsolve.recs.iter().rev() {
        match rec {
            PostRec::Chain(c) => {
                x[c.orig] = (c.b - c.am * x[c.partner]) / c.ak;
            }
            PostRec::Agg(g) => {
                let rest: f64 = g.rest.iter().map(|&(k, c)| c * x[k]).sum();
                x[g.orig] = (g.b - rest) / g.a;
            }
            PostRec::DupCol(d) => {
                // Split the merged value exactly as HiGHS
                // DuplicateColumn::undo (HighsPostsolveStack.cpp:700–737,
                // MIT): FIRST anchor the kept column
                // to a FINITE value — its merge-time lower, or min(0,
                // kept_hi) when that lower is −inf — then propose
                // y = (m − anchor)/scale and trim y into the duplicate's
                // merge-time box, recomputing the kept value as
                // m − scale·y only when the trim fires. The anchor step is
                // what makes the arithmetic safe for infinite merge-time
                // bounds: the naive (m − kept_lo)/s evaluates to ±inf for a
                // free kept column and the clamp never fires (inf > inf is
                // false), leaking ±inf into BOTH halves (seed-90 fuzz
                // catch; the ±1 snap had the same hole). Every branch
                // keeps x_kept + scale·y = m and both values inside their
                // boxes (interval proof: forcing to a dup bound moves the
                // kept value STRICTLY into its box, by the merged-box
                // bounds m ∈ [kept_lo + s·dup_lo, kept_hi + s·dup_hi] for
                // s > 0, mirrored for s < 0).
                let m = x[d.kept];
                let anchor = if d.kept_lo != f64::NEG_INFINITY {
                    d.kept_lo
                } else {
                    0.0f64.min(d.kept_hi)
                };
                let mut y = (m - anchor) / d.scale;
                let mut recompute = false;
                if y > d.dup_hi {
                    y = d.dup_hi;
                    recompute = true;
                } else if y < d.dup_lo {
                    y = d.dup_lo;
                    recompute = true;
                }
                x[d.orig] = y;
                x[d.kept] = if recompute { m - d.scale * y } else { anchor };
            }
            PostRec::Forcing(f) => {
                // Far-side-overlap value (HighsPostsolveStack ForcingColumn::
                // undo, :275+): every stored row bounds x_col from the dir
                // side (soundness note §Soundness); f64::min/max ignore a NaN
                // operand, which is the guard against a poisoned entry.
                let mut val = f.finite_bound;
                let mut sums: Vec<(f64, f64)> = Vec::with_capacity(f.rows.len());
                for fr in &f.rows {
                    let (mut rest, mut aic) = (0.0f64, 0.0f64);
                    for &(co, cf) in &fr.entries {
                        if co == f.orig {
                            aic = cf;
                        } else {
                            rest += cf * x[co];
                        }
                    }
                    let bound = (fr.rhs - rest) / aic;
                    val = if f.dir == 1 {
                        val.min(bound)
                    } else {
                        val.max(bound)
                    };
                    sums.push((aic, rest));
                }
                x[f.orig] = val;
                // Invariant check: every deleted row must hold at the
                // reconstructed x (panic in debug, traced-log in release).
                for (fr, &(aic, rest)) in f.rows.iter().zip(&sums) {
                    let act = aic * val + rest;
                    let slack = TOL * (1.0 + fr.rhs.abs());
                    let ok = match fr.sense {
                        'L' => act <= fr.rhs + slack,
                        'G' => act >= fr.rhs - slack,
                        'E' => (act - fr.rhs).abs() <= slack,
                        _ => true,
                    };
                    if !ok {
                        debug_assert!(
                            false,
                            "forcing col {} row {} '{}' violated: act {act:?} vs rhs {:?}",
                            f.orig, fr.row, fr.sense, fr.rhs
                        );
                        traced(
                            "forcing-check",
                            format!(
                                "forcing col {} row {} '{}' violated: act {act:?} vs rhs {:?}",
                                f.orig, fr.row, fr.sense, fr.rhs
                            ),
                        );
                    }
                }
            }
        }
    }
    (x, reduced_obj + postsolve.obj_const)
}

/// RUSOLVE_PRESOLVE_HASH=1: FNV-1a over the exact bits of the reduced LP so
/// two builds can be proven to apply identical reductions.
fn hash_probe(work: &Lp, live: &[usize], obj_const: f64, nrecs: usize) {
    if std::env::var("RUSOLVE_PRESOLVE_HASH").is_err() {
        return;
    }
    let mut h: u64 = 0xcbf29ce484222325;
    let mut put = |x: u64| {
        h ^= x;
        h = h.wrapping_mul(0x100000001b3);
    };
    put(work.num_row as u64);
    put(work.num_col as u64);
    put(obj_const.to_bits());
    put(nrecs as u64);
    for j in 0..work.num_col {
        put(live[j] as u64);
        put(work.objective[j].to_bits());
        put(work.lower[j].to_bits());
        put(work.upper[j].to_bits());
        // Stored entries, column-major ascending rows (the dense scan order;
        // the reduced LP stores no zeros — `Sp::to_lp` builds via CscBuilder).
        let ptr = work.a.col_ptr();
        for k in ptr[j]..ptr[j + 1] {
            let a = work.a.values()[k];
            put(if a == 0.0 { 0 } else { a.to_bits() });
        }
    }
    for i in 0..work.num_row {
        put(work.rhs[i].to_bits());
        put(work.senses[i] as u64);
    }
    eprintln!(
        "[presolve-hash] {h:016x} rows {} cols {}",
        work.num_row, work.num_col
    );
}

fn terminal(status: SimplexStatus) -> Presolved {
    if std::env::var("RUSOLVE_PRESOLVE_HASH").is_ok() {
        eprintln!("[presolve-hash] terminal {status:?}");
    }
    if std::env::var("RUSOLVE_DEBUG").is_ok() {
        eprintln!("[presolve] terminal: {status:?}");
    }
    Presolved::Terminal { status }
}

/// Trace a rule application when RUSOLVE_DEBUG=1 (finds false infeasibles).
fn traced(rule: &str, what: impl std::fmt::Display) {
    if std::env::var("RUSOLVE_DEBUG").is_ok() {
        eprintln!("[presolve] {rule}: {what}");
    }
}

/// Sparse working LP for presolve.
///
/// Design reference: HiGHS `HPresolve` keeps the matrix as linked row/column
/// lists with row/column deletion flags (HPresolve.cpp, `Alink/Anext/Aprev`,
/// `rowDeleted/colDeleted`; MIT) so that each reduction touches only the
/// entries involved. This is the same idea in idiomatic Rust: two
/// synchronized sorted entry lists (by column per row, by row per column),
/// STABLE ids, alive flags, and lazy deletion — a dead row/column is
/// filtered out by the iterators and physically purged by [`Sp::compact`]
/// once per pass. The previous implementation stored dense `col[j][i]`
/// vectors and rebuilt all of them in `drop_rows`/`Vec::remove`
/// (O(m·n) per deletion; ~40 s of a 43 s solve on MIPLIB `fiball`).
///
/// Numerical identity with the dense version: ids preserve the relative
/// order of surviving rows/columns, entries are visited in ascending id (the
/// dense scan order), and an entry whose value is exactly 0.0 is treated as
/// absent exactly as the dense `a == 0.0` skips did.
struct Sp {
    m: usize,
    n: usize,
    /// Alive counts (the dense `num_row` / `num_col`).
    nrow: usize,
    ncol: usize,
    /// rows[i] = (col, value), ascending col; includes dead cols until compact.
    rows: Vec<Vec<(usize, f64)>>,
    /// cols[j] = (row, value), ascending row; includes dead rows until compact.
    cols: Vec<Vec<(usize, f64)>>,
    ralive: Vec<bool>,
    calive: Vec<bool>,
    rhs: Vec<f64>,
    senses: Vec<char>,
    objective: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    maximize: bool,
}

type Signature = (usize, usize, Vec<f64>, Vec<f64>, Vec<f64>, Vec<char>);

impl Sp {
    fn from_lp(lp: &Lp) -> Sp {
        let (m, n) = (lp.num_row, lp.num_col);
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
        for j in 0..n {
            for k in ptr[j]..ptr[j + 1] {
                let (i, a) = (ridx[k] as usize, rval[k]);
                if a != 0.0 {
                    cols[j].push((i, a));
                    rows[i].push((j, a));
                }
            }
        }
        Sp {
            m,
            n,
            nrow: m,
            ncol: n,
            rows,
            cols,
            ralive: vec![true; m],
            calive: vec![true; n],
            rhs: lp.rhs.clone(),
            senses: lp.senses.clone(),
            objective: lp.objective.clone(),
            lower: lp.lower.clone(),
            upper: lp.upper.clone(),
            maximize: lp.maximize,
        }
    }

    /// Live nonzeros of row `i`, ascending column id.
    fn row(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.rows[i]
            .iter()
            .copied()
            .filter(move |&(j, a)| a != 0.0 && self.calive[j])
    }

    /// Live nonzeros of column `j`, ascending row id.
    fn col(&self, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.cols[j]
            .iter()
            .copied()
            .filter(move |&(i, a)| a != 0.0 && self.ralive[i])
    }

    /// Purge dead entries / exact zeros so repeated scans stay O(live nnz).
    fn compact(&mut self) {
        let (ra, ca) = (&self.ralive, &self.calive);
        for (i, r) in self.rows.iter_mut().enumerate() {
            if ra[i] {
                r.retain(|&(j, a)| a != 0.0 && ca[j]);
            } else {
                *r = Vec::new();
            }
        }
        for (j, c) in self.cols.iter_mut().enumerate() {
            if ca[j] {
                c.retain(|&(i, a)| a != 0.0 && ra[i]);
            } else {
                *c = Vec::new();
            }
        }
    }

    /// `A[r][c] -= x` in both views (inserting fill-in at the sorted spot).
    ///
    /// Exact-cancellation cleanup (HiGHS addToMatrix drops entries below
    /// small_matrix_value): a result that is rounding noise relative to the
    /// two terms that produced it is stored as an exact 0.0 (= absent).
    /// Without it a cancelled coefficient ~1e-16 survives, the pure-bound
    /// fold turns `1e-16 x >= rhs-noise` into a bogus huge bound, and the
    /// LP is silently cut (fuzz seed 73584: x0 <= 1 from a 2-chain pair).
    fn sub_entry(&mut self, r: usize, c: usize, x: f64) -> f64 {
        let new = match self.cols[c].binary_search_by_key(&r, |e| e.0) {
            Ok(p) => {
                let old = self.cols[c][p].1;
                let v = old - x;
                let v = if v.abs() <= SMALL_MATRIX * old.abs().max(x.abs()) {
                    0.0
                } else {
                    v
                };
                self.cols[c][p].1 = v;
                v
            }
            Err(p) => {
                let v = 0.0 - x;
                self.cols[c].insert(p, (r, v));
                v
            }
        };
        match self.rows[r].binary_search_by_key(&c, |e| e.0) {
            Ok(p) => self.rows[r][p].1 = new,
            Err(p) => self.rows[r].insert(p, (c, new)),
        }
        new
    }

    /// Coefficient of `(r, c)` if it is a live nonzero.
    fn entry(&self, r: usize, c: usize) -> Option<f64> {
        let p = self.rows[r].binary_search_by_key(&c, |e| e.0).ok()?;
        let a = self.rows[r][p].1;
        (a != 0.0 && self.calive[c]).then_some(a)
    }

    fn kill_col(&mut self, j: usize) {
        debug_assert!(self.calive[j]);
        self.calive[j] = false;
        self.ncol -= 1;
    }

    fn kill_row(&mut self, i: usize) {
        debug_assert!(self.ralive[i]);
        self.ralive[i] = false;
        self.nrow -= 1;
    }

    /// Delete every alive row with `keep[i] == false`.
    fn drop_rows(&mut self, keep: &[bool]) {
        for (i, &k) in keep.iter().enumerate() {
            if !k && self.ralive[i] {
                self.kill_row(i);
            }
        }
    }

    /// Remove fixed column `j` at value `v`: fold into rhs/obj, drop the column.
    fn substitute(&mut self, j: usize, v: f64, states: &mut [VarState], obj_const: &mut f64) {
        if v != 0.0 {
            for k in 0..self.cols[j].len() {
                let (i, a) = self.cols[j][k];
                if self.ralive[i] {
                    self.rhs[i] -= a * v;
                }
            }
        }
        *obj_const += self.objective[j] * v;
        states[j] = VarState::Fixed(v);
        self.kill_col(j);
    }

    /// Interval-propagation variants of `tighten_*`: a finite bound moves only
    /// if the improvement exceeds `PROP_MIN_GAIN` (relative). Near a
    /// degenerate vertex two rows contract each other's bounds geometrically
    /// (Zeno) and every step amplifies rounding error by the row-coefficient
    /// ratio, which ends with bounds crossing the true point by ~1e-8 and a
    /// false infeasibility. Skipping sub-threshold steps leaves the bound
    /// valid (merely looser) and breaks the cascade. HiGHS uses the same idea
    /// for continuous columns (minimum relative improvement before a
    /// propagated bound is applied, HPresolveAnalysis / `changeColUpper`
    /// callers in HPresolve.cpp, `kHighsInf`-guarded 1e3*feastol rule).
    fn prop_ub(&mut self, j: usize, cap: f64) {
        let u = self.upper[j];
        if u.is_infinite() || u - cap > PROP_MIN_GAIN * (1.0 + cap.abs()) {
            self.tighten_ub(j, cap);
        }
    }

    fn prop_lb(&mut self, j: usize, cap: f64) {
        let l = self.lower[j];
        if l.is_infinite() || cap - l > PROP_MIN_GAIN * (1.0 + cap.abs()) {
            self.tighten_lb(j, cap);
        }
    }

    fn tighten_ub(&mut self, j: usize, cap: f64) {
        if cap < self.upper[j] {
            self.upper[j] = cap;
        }
    }

    fn tighten_lb(&mut self, j: usize, cap: f64) {
        if cap > self.lower[j] {
            self.lower[j] = cap;
        }
    }

    /// Cheap change signature for the fixpoint test. Senses are part of it so
    /// the loop keeps spinning after an implied-equality conversion (the
    /// doubleton chain must get a pass to fire on the converted row); the
    /// conversion is monotone L/G → E, so it can never oscillate the loop.
    /// Dead ids never change, so full-length vectors compare exactly like
    /// the former alive-only vectors.
    fn signature(&self) -> Signature {
        (
            self.ncol,
            self.nrow,
            self.lower.clone(),
            self.upper.clone(),
            self.rhs.clone(),
            self.senses.clone(),
        )
    }

    /// Materialize the surviving LP (CSC `Lp` is the downstream contract)
    /// plus the reduced-col → original-col map. Live entries land in
    /// ascending rank order (both row and column), so the CSC layout matches
    /// the dense-era scatter bit for bit.
    fn to_lp(&self) -> (Lp, Vec<usize>) {
        let mut rank = vec![usize::MAX; self.m];
        let mut rows_alive = Vec::with_capacity(self.nrow);
        for i in 0..self.m {
            if self.ralive[i] {
                rank[i] = rows_alive.len();
                rows_alive.push(i);
            }
        }
        let live: Vec<usize> = (0..self.n).filter(|&j| self.calive[j]).collect();
        let mut b = rusolve_core::CscBuilder::new(rows_alive.len());
        let mut col: Vec<(usize, f64)> = Vec::new();
        for &j in &live {
            col.clear();
            col.extend(self.col(j).map(|(i, a)| (rank[i], a)));
            b.push_col(&col);
        }
        let lp = Lp {
            num_col: live.len(),
            num_row: rows_alive.len(),
            a: b.into_matrix(),
            rhs: rows_alive.iter().map(|&i| self.rhs[i]).collect(),
            objective: live.iter().map(|&j| self.objective[j]).collect(),
            maximize: self.maximize,
            senses: rows_alive.iter().map(|&i| self.senses[i]).collect(),
            lower: live.iter().map(|&j| self.lower[j]).collect(),
            upper: live.iter().map(|&j| self.upper[j]).collect(),
        };
        (lp, live)
    }
}

/// Activity of one row over the CURRENT column bounds, with the number of
/// infinite contributions kept separately so "activity without column j" is
/// available in O(1) (the `HighsLinearSumBounds` idea, rebuilt per pass).
#[derive(Clone, Copy, Default)]
struct RowSum {
    min_fin: f64,
    max_fin: f64,
    inf_min: u32,
    inf_max: u32,
}

fn agg_row_sum(work: &Sp, i: usize) -> RowSum {
    let mut s = RowSum::default();
    for (k, a) in work.row(i) {
        let (p, q) = if a > 0.0 {
            (work.lower[k], work.upper[k])
        } else {
            (work.upper[k], work.lower[k])
        };
        if p.is_finite() {
            s.min_fin += a * p;
        } else {
            s.inf_min += 1;
        }
        if q.is_finite() {
            s.max_fin += a * q;
        } else {
            s.inf_max += 1;
        }
    }
    s
}

/// Bounds on x_j implied by row `i` (coefficient `a`) from the bounds of the
/// row's OTHER columns: `(lo, hi)`, ±∞ when the row gives no finite bound.
fn agg_implied(work: &Sp, s: &RowSum, i: usize, j: usize, a: f64) -> (f64, f64) {
    let (p, q) = if a > 0.0 {
        (work.lower[j], work.upper[j])
    } else {
        (work.upper[j], work.lower[j])
    };
    let rest_min = if p.is_finite() {
        (s.inf_min == 0).then_some(s.min_fin - a * p)
    } else {
        (s.inf_min == 1).then_some(s.min_fin)
    };
    let rest_max = if q.is_finite() {
        (s.inf_max == 0).then_some(s.max_fin - a * q)
    } else {
        (s.inf_max == 1).then_some(s.max_fin)
    };
    let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
    let sense = work.senses[i];
    let r = work.rhs[i];
    // a·x_j + rest ≤ r  (sense L or E)  ⇒  a·x_j ≤ r − rest_min
    if matches!(sense, 'L' | 'E')
        && let Some(rm) = rest_min
    {
        let c = (r - rm) / a;
        if a > 0.0 {
            hi = hi.min(c);
        } else {
            lo = lo.max(c);
        }
    }
    // a·x_j + rest ≥ r  (sense G or E)  ⇒  a·x_j ≥ r − rest_max
    if matches!(sense, 'G' | 'E')
        && let Some(rm) = rest_max
    {
        let c = (r - rm) / a;
        if a > 0.0 {
            lo = lo.max(c);
        } else {
            hi = hi.min(c);
        }
    }
    (lo, hi)
}

/// Is column `j` implied free in the CURRENT model once equality row `r` is
/// used to eliminate it? HiGHS `isImpliedFree` (:207) asks whether the bounds
/// implied by rows already lie inside the column's own bounds. Rebuilt for the
/// sparse matrix as: for each finite bound of x_j, some row of the column
/// (the pivot row itself included — its implied range is exactly the
/// substituted expression's range) must imply a bound at least as tight.
/// Soundness (the reason this is stated against the CURRENT bounds, not the
/// model's original ones): after the substitution every row `s ≠ r` of the
/// column survives with x_j replaced by its row-`r` expression, so a bound
/// implied by `s` (or by `r`) over the surviving columns' bounds still holds
/// for the recovered x_j, and the column's own bound — which is dropped with
/// the column — is then redundant. No circularity: only bounds of columns
/// that remain in the model are used. `sums[i]` must be current for every
/// row of the column.
fn agg_is_free(work: &Sp, sums: &[RowSum], j: usize) -> bool {
    let (l, u) = (work.lower[j], work.upper[j]);
    if l > u {
        return false;
    }
    let need_lo = l.is_finite();
    let need_hi = u.is_finite();
    let (mut best_lo, mut best_hi) = (f64::NEG_INFINITY, f64::INFINITY);
    if need_lo || need_hi {
        for (i, a) in work.col(j) {
            let (lo, hi) = agg_implied(work, &sums[i], i, j, a);
            best_lo = best_lo.max(lo);
            best_hi = best_hi.min(hi);
        }
    }
    (!need_lo || best_lo >= l - AGG_FREE_TOL * (1.0 + l.abs()))
        && (!need_hi || best_hi <= u + AGG_FREE_TOL * (1.0 + u.abs()))
}

/// Relative slack of the implied-free test: the tightening block that wrote
/// the column's bound used the same arithmetic, so equality up to rounding is
/// the expected case.
const AGG_FREE_TOL: f64 = 1e-9;
/// HiGHS `presolve_substitution_maxfillin` default (HighsOptions.h:1652).
const AGG_MAX_FILLIN: i64 = 10;
/// HiGHS `presolve_pivot_threshold` default (HighsOptions.h:1612).
const AGG_PIVOT_THRESH: f64 = 0.01;

/// The AGGREGATOR (HiGHS `HPresolve::aggregator`, HPresolve.cpp:6688–6806,
/// MIT).
/// Eliminates implied-free columns through equality rows:
/// `x_j = (b − Σ_k a_k·x_k)/a_j` is substituted into every other row and into
/// the objective, then row and column are deleted. Candidate order, the
/// "length-2 always substitutes" exemption, the two-sided pivot-size test
/// (|a_j| ≥ 0.01·max|row| or ≥ 0.01·max|col|), the fill-in cap (10) and the
/// give-up after three fill-in failures are the C++'s. Differences: only
/// equality rows count as dual-implied-free (HiGHS also accepts inequality
/// rows whose implied dual bound pins their sign — needs dual bound
/// propagation we do not have), and the implied-free test is `agg_is_free`.
/// Returns the number of eliminated columns.
fn aggregate(work: &mut Sp, obj_const: &mut f64, recs: &mut Vec<PostRec>) -> usize {
    let mut sums: Vec<RowSum> = (0..work.m)
        .map(|i| {
            if work.ralive[i] {
                agg_row_sum(work, i)
            } else {
                RowSum::default()
            }
        })
        .collect();
    // Candidates: (equality row, column) pairs whose column is implied free
    // at pass start, with the sizes HiGHS sorts on.
    struct Cand {
        r: usize,
        j: usize,
        rsz: usize,
        csz: usize,
    }
    let mut cands: Vec<Cand> = Vec::new();
    let csz: Vec<usize> = (0..work.n)
        .map(|j| {
            if work.calive[j] {
                work.col(j).count()
            } else {
                0
            }
        })
        .collect();
    let mut free_cache: Vec<Option<bool>> = vec![None; work.n];
    for r in 0..work.m {
        if !work.ralive[r] || work.senses[r] != 'E' {
            continue;
        }
        let rsz = work.row(r).count();
        if rsz < 2 {
            continue;
        }
        for (j, a) in work.row(r) {
            if a.abs() < 1e-9 || work.lower[j] >= work.upper[j] {
                continue; // fixed / inverted columns belong to other rules
            }
            let free = *free_cache[j].get_or_insert_with(|| agg_is_free(work, &sums, j));
            if free {
                cands.push(Cand {
                    r,
                    j,
                    rsz,
                    csz: csz[j],
                });
            }
        }
    }
    cands.sort_by(|x, y| {
        let (m1, m2) = (x.rsz.min(x.csz), y.rsz.min(y.csz));
        // minLen == 2 first (:6706), then smaller row*col product, then
        // smaller min length; HiGHS' last tie-break is a hash — ours is the
        // deterministic (row, col) order.
        (m2 == 2)
            .cmp(&(m1 == 2))
            .then((x.rsz * x.csz).cmp(&(y.rsz * y.csz)))
            .then(m1.cmp(&m2))
            .then((x.r, x.j).cmp(&(y.r, y.j)))
    });

    let mut n_done = 0usize;
    let mut nfail = 0;
    for c in &cands {
        let (r, j) = (c.r, c.j);
        if !work.ralive[r] || !work.calive[j] {
            continue;
        }
        let Some(a) = work.entry(r, j) else { continue };
        let rsz = work.row(r).count();
        let csz = work.col(j).count();
        if rsz < 2 || work.lower[j] >= work.upper[j] || !agg_is_free(work, &sums, j) {
            continue;
        }
        // Fill-in is never a problem when the row or the column has length
        // two (:6747); otherwise gate on pivot size and fill-in.
        if rsz != 2 && csz != 2 {
            let max_row = || work.row(r).fold(0.0f64, |m, (_, v)| m.max(v.abs()));
            let max_col = || work.col(j).fold(0.0f64, |m, (_, v)| m.max(v.abs()));
            let first = if rsz < csz { max_row() } else { max_col() };
            if a.abs() < first * AGG_PIVOT_THRESH {
                let second = if rsz < csz { max_col() } else { max_row() };
                if a.abs() < second * AGG_PIVOT_THRESH {
                    continue;
                }
            }
            let mut fillin: i64 = -((rsz + csz) as i64 - 1);
            'fill: for (s, _) in work.col(j) {
                if s == r {
                    continue;
                }
                for (k, _) in work.row(r) {
                    if work.entry(s, k).is_none() {
                        fillin += 1;
                        if fillin > AGG_MAX_FILLIN {
                            break 'fill;
                        }
                    }
                }
            }
            if fillin > AGG_MAX_FILLIN {
                nfail += 1;
                if nfail == 3 {
                    break; // rows/columns are getting too dense (:6779)
                }
                continue;
            }
            nfail = 0;
        }
        // ---- substitute ----
        let b = work.rhs[r];
        let rest: Vec<(usize, f64)> = work.row(r).filter(|&(k, _)| k != j).collect();
        let colj: Vec<(usize, f64)> = work.col(j).filter(|&(s, _)| s != r).collect();
        let cj = work.objective[j];
        if cj != 0.0 {
            for &(k, ark) in &rest {
                let (old, delta) = (work.objective[k], cj * ark / a);
                let new = old - delta;
                // Cancellation cleanup (HiGHS substitute() zeroes costs below
                // small_matrix_value): an exact-in-theory zero left as ±1e-16
                // turns a harmless free column into a bogus "unbounded
                // direction" in the empty-column rule (fuzz seed 51561).
                work.objective[k] = if new.abs() <= SMALL_MATRIX * old.abs().max(delta.abs()) {
                    0.0
                } else {
                    new
                };
            }
            *obj_const += cj * b / a;
        }
        for &(s, asj) in &colj {
            let f = asj / a;
            for &(k, ark) in &rest {
                let delta = f * ark;
                work.sub_entry(s, k, delta);
            }
            work.rhs[s] -= f * b;
        }
        traced(
            "aggregator",
            format!(
                "col {j} eliminated via equality row {r} (row len {rsz}, col len {csz}, pivot {a:?})"
            ),
        );
        recs.push(PostRec::Agg(AggRec {
            orig: j,
            a,
            b,
            rest,
        }));
        work.kill_col(j);
        work.kill_row(r);
        for &(s, _) in &colj {
            sums[s] = agg_row_sum(work, s);
        }
        n_done += 1;
    }
    n_done
}

/// Row-dual box seeds by sense (HiGHS impliedDualRowBounds init), per row id.
fn dual_seeds(work: &Sp) -> Vec<(f64, f64)> {
    work.senses
        .iter()
        .map(|&s| match s {
            'E' => (f64::NEG_INFINITY, f64::INFINITY), // free dual
            'G' => (0.0, f64::INFINITY),
            'L' => (f64::NEG_INFINITY, 0.0),
            _ => (0.0, 0.0), // 'N' (no constraint): dual pinned at 0
        })
        .collect()
}

/// Dominated-column detection probe (detection only, no mutation).
///
/// Ported detection logic from HiGHS HPresolve.cpp `detectDominatedCol`
/// :4514–4551 (MIT), sense-only row dual seeds = `impliedDualRowBounds`
/// Orig sums; soundness: r_j strictly positive/negative over the dual
/// box ⇒ some optimum at that bound.
///
/// Per column j the reduced-cost interval r_j = c_j − Σ_i a_ij π_i is
/// bounded over the per-row dual box seeded by row sense alone ('E' free,
/// 'G' [0,+∞), 'L' (−∞,0], 'N' [0,0]). A column whose r_j is strictly
/// positive over the whole box has some optimum at its lower bound
/// (candidate +1); strictly negative ⇒ some optimum at its upper bound
/// (candidate −1). Columns already fixed are skipped; a strong bound with
/// an infinite target bound is a dual-infeasibility signal (HiGHS
/// Result::kDualInfeasible path), reported separately and never a candidate.
///
/// Recomputed from scratch each pass (no incremental state). The infinite
/// seeds are safe under plain f64: zero entries are never visited (no
/// 0·∞ = NaN), and each max term is ≥ 0 (possibly +∞) while each min term
/// is ≤ 0 (possibly −∞), so the sums cannot cancel infinities.
///
/// TODO: propagate row dual bounds via strictly-implied column
/// bounds (HiGHS impliedDualRowBounds update); the sense-only seeds here
/// are just the Orig sums.
///
/// Returns (candidates, near-miss count, dual-infeasible-signal count),
/// where a near-miss is a fired candidate whose margin lies within
/// [DUAL_TOL, 10·DUAL_TOL] — the sanity signal that the tolerance choice
/// matters (a 10× looser gate would drop these).
fn dominated_col_candidates(
    work: &Sp,
    dual_tol: f64,
) -> (
    Vec<(usize, /*+1 fix-to-lb, −1 fix-to-ub*/ i8)>,
    usize,
    usize,
) {
    let seeds = dual_seeds(work);
    let mut cands = Vec::new();
    let mut near_miss = 0usize;
    let mut dual_infeasible = 0usize;
    for j in 0..work.n {
        if !work.calive[j] || work.lower[j] == work.upper[j] {
            continue; // dead / already fixed
        }
        // Min-space cost (Lp doc: maximize internally minimizes −cᵀx).
        let c = if work.maximize {
            -work.objective[j]
        } else {
            work.objective[j]
        };
        // r_lower = c − Σ max(a·π), r_upper = c − Σ min(a·π) over each
        // row's dual box [dlo, dhi].
        let (mut m_sum, mut n_sum) = (0.0f64, 0.0f64);
        for (i, a) in work.col(j) {
            let (dlo, dhi) = seeds[i];
            if a > 0.0 {
                m_sum += a * dhi;
                n_sum += a * dlo;
            } else {
                m_sum += a * dlo;
                n_sum += a * dhi;
            }
        }
        let r_lower = c - m_sum;
        let r_upper = c - n_sum;
        if r_lower > dual_tol {
            if work.lower[j] == f64::NEG_INFINITY {
                dual_infeasible += 1; // HiGHS kDualInfeasible path: no fix
            } else {
                cands.push((j, 1));
                if r_lower <= 10.0 * dual_tol {
                    near_miss += 1;
                }
            }
        }
        if r_upper < -dual_tol {
            if work.upper[j] == f64::INFINITY {
                dual_infeasible += 1; // same ±∞ guard on the upper side
            } else {
                cands.push((j, -1));
                if -r_upper <= 10.0 * dual_tol {
                    near_miss += 1;
                }
            }
        }
    }
    (cands, near_miss, dual_infeasible)
}

/// Detection: FORCING columns (Andersen & Andersen 1995; HiGHS
/// weaklyDominatedCol forcing fallback :4553–4617, MIT). A col with an
/// INFINITE bound whose dual is pinned at exactly zero (|orig col-dual
/// bound| ≤ dual_tol/dynamism) can be deleted together with ALL its rows:
/// the dual-zero ⇒ objective indifferent, and the sign structure of the
/// zero condition guarantees the removed rows reconstruct by choosing
/// x_col as the far-side-overlap value (Andersen & Andersen 1995).
/// Returns (col, dir) candidates, dir +1
/// for the lb = −∞ side, −1 for the ub = +∞ side (a fully free col may
/// appear in both); the stats counters derive from this vec.
fn forcing_col_candidates(work: &Sp, dual_tol: f64) -> Vec<(usize, /*+1 / −1*/ i8)> {
    let seeds = dual_seeds(work);
    let mut cands = Vec::new();
    for j in 0..work.n {
        if !work.calive[j] || work.lower[j] == work.upper[j] {
            continue;
        }
        let c = if work.maximize {
            -work.objective[j]
        } else {
            work.objective[j]
        };
        let (mut m_sum, mut n_sum) = (0.0f64, 0.0f64);
        let (mut amax, mut amin) = (0.0f64, f64::INFINITY);
        for (i, a) in work.col(j) {
            let (dlo, dhi) = seeds[i];
            if a > 0.0 {
                m_sum += a * dhi;
                n_sum += a * dlo;
            } else {
                m_sum += a * dlo;
                n_sum += a * dhi;
            }
            let aa = a.abs();
            amax = amax.max(aa);
            amin = amin.min(aa);
        }
        if amax == 0.0 {
            continue; // empty column
        }
        let dynamism = amax / amin; // HiGHS computeDynamism ratio form
        let r_lower = c - m_sum;
        let r_upper = c - n_sum;
        // direction +1: infinite LOWER bound, dual box lower ≈ 0 (:4581)
        if work.lower[j] == f64::NEG_INFINITY && r_lower.abs() <= dual_tol / dynamism {
            cands.push((j, 1));
        }
        // direction −1: infinite UPPER bound, dual box upper ≈ 0 (:4585)
        if work.upper[j] == f64::INFINITY && r_upper.abs() <= dual_tol / dynamism {
            cands.push((j, -1));
        }
    }
    cands
}

/// Signed max-abs scale of a sparse vector + its support — the
/// canonical normalization for parallel detection, ported from HiGHS
/// `rowMax`/`colMax` collection in `detectParallelRowsAndCols`
/// HPresolve.cpp:7298–7320 (MIT). Among |values| equal within
/// `SMALL_MATRIX` the entry with the SMALLEST index wins, so proportional
/// vectors (any sign) normalize to the identical coefficient vector.
/// Takes the nonzeros in ascending index order.
fn scale_normal(entries: impl Iterator<Item = (usize, f64)>) -> (Vec<usize>, f64) {
    let mut sup = Vec::new();
    let (mut sv, mut si) = (0.0f64, usize::MAX);
    for (idx, a) in entries {
        if a == 0.0 {
            continue;
        }
        sup.push(idx);
        let av = a.abs();
        if av > sv.abs() + SMALL_MATRIX || ((av - sv.abs()).abs() <= SMALL_MATRIX && idx < si) {
            sv = a;
            si = idx;
        }
    }
    (sup, sv)
}

/// Exact proportionality re-check inside one support bucket: normalized
/// coefficients agree within `PAR_TOL` (HiGHS raw-value re-check
/// :7620–7632, normalized form). Both value slices follow the shared
/// support order.
fn proportional(av: &[f64], bv: &[f64], sa: f64, sb: f64) -> bool {
    av.iter()
        .zip(bv)
        .all(|(&a, &b)| (a / sa - b / sb).abs() <= PAR_TOL)
}

/// Bucket item: (id, scale-normalized scale, values in support order).
type ScaleItem = (usize, f64, Vec<f64>);

/// Detection probe → interval merge: PARALLEL ROW groups —
/// identical support, proportional coefficients (HiGHS bucket pass
/// :7267–7450). Empty rows are skipped (the empty-row rule owns them).
/// The bucket key dropped the raw sense — HiGHS keeps rows as
/// intervals [row_lower, row_upper] and merges the intervals of parallel
/// rows regardless of how the two sides were written; the sense-mixing is
/// now resolved in the mutation (effective-interval intersection), not in
/// detection. A cross-sense pair whose interval COLLAPSES to a point becomes
/// an equality (graph20 class); a true range keeps one witness per side.
fn parallel_row_groups(work: &Sp) -> Vec<Vec<usize>> {
    use std::collections::BTreeMap;
    // support → candidates; support equality is necessary for
    // proportionality, so it is the bucket key (sense is NOT part of it —
    // see doc above).
    let mut buckets: BTreeMap<Vec<usize>, Vec<ScaleItem>> = BTreeMap::new();
    for i in 0..work.m {
        if !work.ralive[i] {
            continue;
        }
        let (sup, scale) = scale_normal(work.row(i));
        if sup.is_empty() {
            continue; // empty row: the emptiness rule deletes it anyway
        }
        if work.senses[i] == 'N' {
            continue; // free row: nothing to bound
        }
        let vals: Vec<f64> = work.row(i).map(|e| e.1).collect();
        buckets.entry(sup).or_default().push((i, scale, vals));
    }
    let mut out = Vec::new();
    for (_, items) in buckets {
        if items.len() < 2 {
            continue;
        }
        // groups: (representative item index, its scale, members)
        let mut groups: Vec<(usize, f64, Vec<usize>)> = Vec::new();
        for (k, (i, si, vi)) in items.iter().enumerate() {
            if let Some(g) = groups
                .iter_mut()
                .find(|(r, sr, _)| proportional(vi, &items[*r].2, *si, *sr))
            {
                g.2.push(*i);
            } else {
                groups.push((k, *si, vec![*i]));
            }
        }
        for (_, _, members) in groups {
            if members.len() >= 2 {
                out.push(members);
            }
        }
    }
    out
}

/// Detection probe (no mutation): PARALLEL COLUMN groups —
/// identical support, proportional coefficients, and equivalent costs in
/// the merge class (|c_rep·s − c_j| ≤ `DUAL_TOL` in min-space, HiGHS
/// objDiff gate :7527–7557): only cost-equivalent parallels are mergeable,
/// the domination cases belong to the dominance rule and are not counted
/// here. Fixed
/// (substitution owns them) and empty columns are skipped.
///
/// `cost_class = false` drops the objDiff restriction and returns ALL
/// proportional groups (any relative cost) — the candidate set for the
/// parallel-column DOMINANCE arms (HiGHS :7469–7540 continuous-
/// continuous branch pairs any proportional columns and decides merge vs
/// domination per pair).
fn parallel_col_groups_filter(work: &Sp, cost_class: bool) -> Vec<Vec<usize>> {
    use std::collections::BTreeMap;
    let cost = |j: usize| {
        if work.maximize {
            -work.objective[j]
        } else {
            work.objective[j]
        }
    };
    let mut buckets: BTreeMap<Vec<usize>, Vec<ScaleItem>> = BTreeMap::new();
    for j in 0..work.n {
        if !work.calive[j] || work.lower[j] == work.upper[j] {
            continue; // dead / fixed: the substitution rule owns it
        }
        let (sup, s) = scale_normal(work.col(j));
        if sup.is_empty() {
            continue; // empty column: the empty-column rule owns it
        }
        let vals: Vec<f64> = work.col(j).map(|e| e.1).collect();
        buckets.entry(sup).or_default().push((j, s, vals));
    }
    let mut out = Vec::new();
    for (_, items) in buckets {
        if items.len() < 2 {
            continue;
        }
        let mut groups: Vec<(usize, f64, Vec<usize>)> = Vec::new();
        for (k, (j, sj, vj)) in items.iter().enumerate() {
            if let Some(g) = groups.iter_mut().find(|(r, sr, _)| {
                let rc = items[*r].0;
                proportional(vj, &items[*r].2, *sj, *sr)
                    && (!cost_class || (cost(rc) * (*sj / *sr) - cost(*j)).abs() <= DUAL_TOL)
            }) {
                g.2.push(*j);
            } else {
                groups.push((k, *sj, vec![*j]));
            }
        }
        for (_, _, members) in groups {
            if members.len() >= 2 {
                out.push(members);
            }
        }
    }
    out
}

fn parallel_col_groups(work: &Sp) -> Vec<Vec<usize>> {
    parallel_col_groups_filter(work, true)
}

/// Row-implied bound caps on one column from the CURRENT bounds of every
/// other column (the same residual arithmetic as the pass's tightening
/// block, recomputed on demand). Returns (implied_lower, implied_upper);
/// a side no row caps is ±∞. This is our recomputed analogue of HiGHS's
/// incrementally-maintained `implColLower/Upper` (same current-bounds
/// basis, so the dominance gates below see exactly what HiGHS's do).
fn implied_col_bounds(work: &Sp, j: usize) -> (f64, f64) {
    let mut lo = f64::NEG_INFINITY;
    let mut hi = f64::INFINITY;
    for (r, a) in work.col(j) {
        let (mut others_min, mut others_max) = (0.0f64, 0.0f64);
        let (mut inf_min, mut inf_max) = (false, false);
        for (k, b) in work.row(r) {
            if k == j {
                continue;
            }
            // contribution of x_k to the row activity: min = b·(bound the
            // sign picks), max mirrored — identical to the tightening block.
            let (p, q) = if b > 0.0 {
                (work.lower[k], work.upper[k])
            } else {
                (work.upper[k], work.lower[k])
            };
            if p.is_finite() {
                others_min += b * p;
            } else {
                inf_min = true;
            }
            if q.is_finite() {
                others_max += b * q;
            } else {
                inf_max = true;
            }
        }
        let rhs = work.rhs[r];
        // Universal caps: divide the row by a and pick the others-side that
        // bounds the quotient for EVERY achievable others value (E-row
        // arithmetic below; one-sided rows take only the side they imply).
        let (cap_lo, cap_hi) = match work.senses[r] {
            'L' => {
                // Σ a·x ≤ rhs: a>0 ⇒ x_j ≤ (rhs − others_min)/a;
                //           a<0 ⇒ x_j ≥ (rhs − others_min)/a (divide flips).
                (
                    if a < 0.0 && !inf_min {
                        Some((rhs - others_min) / a)
                    } else {
                        None
                    },
                    if a > 0.0 && !inf_min {
                        Some((rhs - others_min) / a)
                    } else {
                        None
                    },
                )
            }
            'G' => (
                if a > 0.0 && !inf_max {
                    Some((rhs - others_max) / a)
                } else {
                    None
                },
                if a < 0.0 && !inf_max {
                    Some((rhs - others_max) / a)
                } else {
                    None
                },
            ),
            'E' => (
                if (!inf_max && a > 0.0) || (!inf_min && a < 0.0) {
                    if a > 0.0 {
                        Some((rhs - others_max) / a)
                    } else {
                        Some((rhs - others_min) / a)
                    }
                } else {
                    None
                },
                if (!inf_min && a > 0.0) || (!inf_max && a < 0.0) {
                    if a > 0.0 {
                        Some((rhs - others_min) / a)
                    } else {
                        Some((rhs - others_max) / a)
                    }
                } else {
                    None
                },
            ),
            _ => (None, None),
        };
        if let Some(c) = cap_lo {
            lo = lo.max(c);
        }
        if let Some(c) = cap_hi {
            hi = hi.min(c);
        }
    }
    (lo, hi)
}

/// Parallel-column DOMINANCE fixing — the LP branch of HiGHS
/// `detectParallelRowsAndCols` (HPresolve.cpp:7399–7540 case table,
/// :7601–7633 fixCol arms, MIT License, HiGHS/ERGO Institute). For two proportional columns
/// `dup = scale·col` (coefficient-wise, any real scale — the LP branch
/// pairs `col` = the earlier candidate, HiGHS :7477–7486), the pair is
/// decided by objDiff = c_col·scale − c_dup (min-space):
///
/// * objDiff < −tol (the scaled col is cheaper): fix the expensive `dup`
///   at its lower bound when `col` can absorb the compensation (its
///   opposing bound is infinite or strictly row-implied — `colUpperInf`,
///   :7399–7411); otherwise fix `col` at a bound when `dup` can absorb
///   (:7503–7511).
/// * objDiff > tol: mirrored (:7514–7522).
/// * |objDiff| ≤ tol: the same arms first (:7525–7539); the merge itself
///   is the separate default-OFF parcol rule.
///
/// Soundness (primal, our no-duals posture — HiGHS needs the STRICT
/// margin only for dual postsolve, :7401–7405): the compensating column's
/// bound must be implied by rows FROM THE CURRENT BOUNDS — then for any
/// feasible point with the fixed column shifted to its bound, the
/// compensated point stays inside every row's valid deduction, hence
/// inside its explicit bound. The fixed value is value-known, so removal
/// rides the existing `fixed` substitution: no new postsolve record
/// (fixColToLower/Upper, :7602–7617).
///
/// Returns (col, value) fixes; HiGHS applies them inside its scan (later
/// pairs see earlier fixings), so callers must too.
fn parallel_col_dominance_fixes(work: &mut Sp) -> Vec<(usize, f64)> {
    // HiGHS primal_feasibility_tolerance (default 1e-7): the tolerance for
    // the compensation gate below.
    const PRIMAL_TOL: f64 = 1e-7;
    let cost = |w: &Sp, j: usize| {
        if w.maximize {
            -w.objective[j]
        } else {
            w.objective[j]
        }
    };
    // Compensation identity: with A_dup = scale·A_col, a shift (Δdup, Δcol)
    // preserving every row satisfies Δdup + scale·Δcol = 0, and the
    // objective moves by Δobj = (Δcol/scale)·objDiff. An arm is sound iff
    // (a) the fixed bound is finite, (b) the objective condition holds in
    // the direction the fix moves the pair (up to the DUAL_TOL drift class
    // HiGHS's objDiff gate accepts), and (c) the COMPENSATING column's
    // movement direction is absorbable: its bound on that side is infinite
    // or implied by rows from the current bounds. Recomputed implied caps
    // are valid deductions, and every current bound is itself a valid
    // deduction (monotone tightening under equivalence-preserving rules),
    // so a finite implied cap ≤ the current bound certifies compensation
    // (the WEAK form; HiGHS's strict margin guards dual postsolve, which
    // our primal-only postsolve does not emit). This re-expresses the LP
    // branch of HiGHS's case table (:7503–7539) with each arm stating its
    // own direction — a literal transcription flipped the col-side arms'
    // absorb directions under a negative scale (fuzz seed 145: an
    // anti-parallel cost-equal pair fixed x1 := 0, −9 → 0).
    let absorbs = |w: &Sp, j: usize, up: bool| {
        if up {
            let ub = w.upper[j];
            !ub.is_finite() || implied_col_bounds(w, j).1 <= ub + PRIMAL_TOL
        } else {
            let lb = w.lower[j];
            !lb.is_finite() || implied_col_bounds(w, j).0 >= lb - PRIMAL_TOL
        }
    };
    let groups = parallel_col_groups_filter(work, false);
    let mut fixes: Vec<(usize, f64)> = Vec::new();
    let fix = |w: &mut Sp, j: usize, v: f64, why: String, out: &mut Vec<(usize, f64)>| {
        traced(
            "par-col-dom",
            format!("orig col {j} fixed at {v:?} ({why})"),
        );
        w.lower[j] = v;
        w.upper[j] = v;
        out.push((j, v));
    };
    for grp in groups {
        // Pair the group's first member (the detection representative) with
        // each other member — the (earlier candidate, current) orientation
        // of the HiGHS hash walk. Non-representative pairs resolve through
        // the fixpoint loop once a member is fixed.
        let rep = grp[0];
        for &other in &grp[1..] {
            if work.lower[rep] == work.upper[rep] || work.lower[other] == work.upper[other] {
                continue; // a fix in this block already pinned one side
            }
            let (col, dup) = (rep, other);
            let s_col = scale_normal(work.col(col)).1;
            let s_dup = scale_normal(work.col(dup)).1;
            let scale = s_dup / s_col;
            let obj_diff = cost(work, col) * scale - cost(work, dup);
            let near = obj_diff.abs() <= DUAL_TOL;
            let s_od = scale * obj_diff;
            // dup@lb (dup moves down ⇒ col moves up for scale>0):
            // objective needs objDiff ≤ DUAL_TOL.
            if work.lower[dup].is_finite()
                && obj_diff <= DUAL_TOL
                && absorbs(work, col, scale > 0.0)
            {
                let why = format!(
                    "{} parallel of col {col}, fixed at lower, scale {scale:.4}",
                    if near { "cost-equal" } else { "dominated" }
                );
                fix(work, dup, work.lower[dup], why, &mut fixes);
            }
            // dup@ub (dup up ⇒ col down for scale>0): objDiff ≥ −DUAL_TOL.
            else if work.upper[dup].is_finite()
                && obj_diff >= -DUAL_TOL
                && absorbs(work, col, scale < 0.0)
            {
                let why = format!(
                    "{} parallel of col {col}, fixed at upper, scale {scale:.4}",
                    if near { "cost-equal" } else { "dominated" }
                );
                fix(work, dup, work.upper[dup], why, &mut fixes);
            }
            // col@lb (col down ⇒ dup up for scale>0): s·objDiff ≥ −DUAL_TOL.
            else if work.lower[col].is_finite()
                && s_od >= -DUAL_TOL
                && absorbs(work, dup, scale > 0.0)
            {
                let why = format!("dominates parallel col {dup}, fixed at lower, scale {scale:.4}");
                fix(work, col, work.lower[col], why, &mut fixes);
            }
            // col@ub (col up ⇒ dup down for scale>0): s·objDiff ≤ DUAL_TOL.
            else if work.upper[col].is_finite()
                && s_od <= DUAL_TOL
                && absorbs(work, dup, scale < 0.0)
            {
                let why = format!("dominates parallel col {dup}, fixed at upper, scale {scale:.4}");
                fix(work, col, work.upper[col], why, &mut fixes);
            }
        }
    }
    fixes
}

/// Singleton-cost VUB implied-equality chain. HiGHS converts a ranged
/// row to an equation when its DUAL bound is bounded away from zero
/// (`isImpliedEquationAtLower/Upper`, HPresolve.cpp:263–281, conversion at
/// :3567–3588, MIT License, HiGHS/ERGO Institute).
/// Its dual bounds come from `impliedDualRowBounds` propagation, which we
/// do not have; this rule derives the SAME row-dual
/// bound from a column singleton's reduced-cost sign condition alone.
///
/// Let b be a column singleton (its only row is r, coefficient a ≠ 0)
/// with min-space cost c_b ≠ 0, r of sense 'L'/'G'. WLOG c_b > 0 (min):
/// the improving side is b's LOWER bound. If that bound is implied by the
/// current system (`implied_col_bounds`, the same valid-deduction basis as
/// the dominance gate), relaxing it to −∞ leaves the feasible set
/// unchanged, and then EVERY optimal solution has reduced cost
/// d_b = c_b − a·y_r ≤ 0 (b can only sit at its upper bound or basic):
/// for sign(a) = −σ(row sense 'L', y_r ≤ 0) this forces
/// y_r ≤ c_b/a < 0 strictly — bounded away from 0 by |c_b/a| — so by
/// complementary slackness r is tight at every optimum ⇒ converting r to
/// 'E' keeps every optimum (it may drop suboptimal slack points only).
/// The c_b < 0 arm relaxes the UPPER bound and mirrors; 'G' rows mirror
/// ('L' needs sign(a) = −sign(c_b), 'G' needs +).
///
/// UNBOUNDED CAVEAT: the argument presumes an optimum exists. On an
/// unbounded LP it is vacuous; the conversion could in principle empty
/// the feasible set and misreport Infeasible. Fuzz arbitrates (the
/// randomized test below includes both statuses); no guard unless it
/// fires in practice.
///
/// Fires once per (b, r): the conversion makes r an equation, which this
/// rule and every other rule treats as terminal (senses never move back).
/// The DTE block later in the same pass substitutes b out (it prefers
/// the col singleton, HiGHS :3081–3086), so no postsolve record is
/// needed beyond DTE's own ChainRec.
fn singleton_cost_eq_conversions(work: &mut Sp) -> usize {
    // HiGHS primal_feasibility_tolerance (1e-7): tolerance for the
    // implied-bound gate (weak form — same as the dominance gate).
    const PRIMAL_TOL: f64 = 1e-7;
    let mut n = 0usize;
    for b in 0..work.n {
        if !work.calive[b] || work.lower[b] == work.upper[b] {
            continue; // dead / already fixed — nothing to bound-shift
        }
        let cb = if work.maximize {
            -work.objective[b]
        } else {
            work.objective[b]
        };
        if cb == 0.0 {
            continue; // no dual leverage
        }
        let col: Vec<(usize, f64)> = work.col(b).take(2).collect();
        if col.len() != 1 {
            continue; // not a column singleton
        }
        let (r, a) = col[0];
        let s = work.senses[r];
        if s != 'L' && s != 'G' {
            continue; // already an equation / free row
        }
        // Row-dual sign condition: 'L' ⇒ y_r ≤ 0 needs sign(a) = −σ,
        // 'G' ⇒ y_r ≥ 0 needs sign(a) = +σ (doc comment above).
        let sigma = cb.signum();
        let sign_ok = match s {
            'L' => a * sigma < 0.0,
            _ => a * sigma > 0.0,
        };
        if !sign_ok {
            continue;
        }
        // Margin: every optimum has |y_r| ≥ |c_b/a|; require it above the
        // dual tolerance (HiGHS gates isImpliedEquation at
        // dual_feasibility_tolerance, HConst defaults).
        if cb.abs() / a.abs() <= DUAL_TOL {
            continue;
        }
        // The improving-side bound must be implied by the current system
        // (weak form, as in the dominance gate): relaxing it then leaves
        // the feasible set unchanged.
        let (ilo, ihi) = implied_col_bounds(work, b);
        if cb > 0.0 {
            let lb = work.lower[b];
            if lb.is_finite() && ilo < lb - PRIMAL_TOL {
                continue;
            }
        } else {
            let ub = work.upper[b];
            if ub.is_finite() && ihi > ub + PRIMAL_TOL {
                continue;
            }
        }
        traced(
            "sing-eq",
            format!(
                "col {b} singleton (cost {cb:?}) in row {r} '{s}': \
                 |y| ≥ {:.3e} ⇒ equation; improving bound relaxed",
                cb.abs() / a.abs()
            ),
        );
        if cb > 0.0 {
            work.lower[b] = f64::NEG_INFINITY;
        } else {
            work.upper[b] = f64::INFINITY;
        }
        work.senses[r] = 'E';
        n += 1;
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lp_of(
        rows: &[&[f64]],
        rhs: Vec<f64>,
        senses: &str,
        obj: Vec<f64>,
        lo: Vec<f64>,
        up: Vec<f64>,
    ) -> Lp {
        let n = obj.len();
        let m = rhs.len();
        let triplets: Vec<(usize, usize, f64)> = rows
            .iter()
            .enumerate()
            .flat_map(|(i, row)| {
                // `take(n)`: the dense cm truncated over-long rows silently;
                // keep that leniency for hand-written test data.
                row.iter().enumerate().take(n).map(move |(j, &v)| (i, j, v))
            })
            .collect();
        Lp {
            num_col: n,
            num_row: m,
            a: rusolve_core::CscMatrix::from_triplets(m, n, &triplets),
            rhs,
            objective: obj,
            maximize: false,
            senses: senses.chars().collect(),
            lower: lo,
            upper: up,
        }
    }

    /// Presolve ON must equal presolve OFF on objective & feasibility.
    fn check_agree(lp: &Lp) -> f64 {
        use crate::simplex::{SimplexConfig, SimplexSolver};
        let off_cfg = SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        };
        let on_cfg = SimplexConfig {
            presolve: true,
            ..SimplexConfig::default()
        };
        let a = SimplexSolver::new(off_cfg).solve(lp).unwrap();
        let b = SimplexSolver::new(on_cfg).solve(lp).unwrap();
        assert_eq!(
            a.status, b.status,
            "status mismatch: {:?} vs {:?}",
            a.status, b.status
        );
        if a.status == SimplexStatus::Optimal {
            assert!(
                (a.objective - b.objective).abs() <= 1e-6 * a.objective.abs().max(1.0),
                "obj mismatch: {} vs {}",
                a.objective,
                b.objective
            );
            // postsolve x must be feasible in the ORIGINAL lp
            let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
            let mut act = vec![0.0f64; lp.num_row];
            for j in 0..lp.num_col {
                let xj = b.x[j];
                if xj != 0.0 {
                    for k in ptr[j]..ptr[j + 1] {
                        let v = rval[k];
                        if v != 0.0 {
                            act[ridx[k] as usize] += v * xj;
                        }
                    }
                }
            }
            for (i, s) in lp.senses.iter().enumerate() {
                let lhs = act[i];
                let r = lp.rhs[i];
                match s {
                    'L' => assert!(lhs <= r + 1e-6, "row {i} violated: {lhs} > {r}"),
                    'G' => assert!(lhs >= r - 1e-6, "row {i} violated: {lhs} < {r}"),
                    _ => assert!((lhs - r).abs() <= 1e-6, "row {i} violated: {lhs} != {r}"),
                }
            }
            for j in 0..lp.num_col {
                assert!(
                    b.x[j] >= lp.lower[j] - 1e-6 && b.x[j] <= lp.upper[j] + 1e-6,
                    "bound {j} violated: {:?}",
                    b.x
                );
            }
        }
        b.objective
    }

    #[test]
    fn fixed_variable_substitution() {
        // min x + y s.t. x + y >= 5, x fixed at 2 -> y=3, obj 5
        let lp = lp_of(
            &[&[1.0, 1.0]],
            vec![5.0],
            "G",
            vec![1.0, 1.0],
            vec![2.0, 0.0],
            vec![2.0, 10.0],
        );
        assert!((check_agree(&lp) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn single_nonzero_e_fixes() {
        // min y s.t. 2x = 6 (fixes x=3), x + y >= 4 -> y=1, obj 1
        let lp = lp_of(
            &[&[2.0, 0.0], &[1.0, 1.0]],
            vec![6.0, 4.0],
            "EG",
            vec![0.0, 1.0],
            vec![0.0, 0.0],
            vec![10.0, 10.0],
        );
        assert!((check_agree(&lp) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn forcing_row_fixes_all() {
        // x <= 2 (G row: x >= 2… construct forcing):
        // rows: x <= 2 and x >= 2 meet at 2; wrap into objective:
        // min x+y s.t. x <= 2, -x <= -2, y >= 5-x -> x=2 forced, y=3
        let lp = lp_of(
            &[&[1.0, 0.0], &[-1.0, 0.0], &[1.0, 1.0]],
            vec![2.0, -2.0, 5.0],
            "LLG",
            vec![1.0, 1.0],
            vec![0.0, 0.0],
            vec![9.0, 9.0],
        );
        assert!((check_agree(&lp) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn redundant_row_dropped() {
        // min x s.t. x <= 100 (redundant given ub 5), x >= 2 -> 2
        let lp = lp_of(&[&[1.0]], vec![100.0], "L", vec![1.0], vec![0.0], vec![5.0]);
        assert!((check_agree(&lp) - 0.0).abs() < 1e-6);
        // with a G row to make it nontrivial
        let lp2 = lp_of(
            &[&[1.0], &[1.0]],
            vec![100.0, 2.0],
            "LG",
            vec![1.0],
            vec![0.0],
            vec![5.0],
        );
        assert!((check_agree(&lp2) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn infeasible_row_detected() {
        // min 0 s.t. x <= 1, x >= 3
        let lp = lp_of(
            &[&[1.0], &[1.0]],
            vec![1.0, 3.0],
            "LG",
            vec![0.0],
            vec![0.0],
            vec![9.0],
        );
        let off = {
            use crate::simplex::{SimplexConfig, SimplexSolver};
            SimplexSolver::new(SimplexConfig {
                presolve: false,
                ..Default::default()
            })
            .solve(&lp)
            .unwrap()
        };
        assert_eq!(off.status, SimplexStatus::Infeasible);
        check_agree(&lp);
    }

    #[test]
    fn empty_column_fixed_at_bound() {
        // min -z s.t. x + y >= 3 (z unconstrained in rows), z <= 7 -> z=7
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0]],
            vec![3.0],
            "G",
            vec![0.0, 0.0, -1.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 7.0],
        );
        assert!((check_agree(&lp) + 7.0).abs() < 1e-6);
    }

    #[test]
    fn empty_column_unbounded() {
        // min -z, z unbounded above, no rows touch z
        let lp = lp_of(&[], vec![], "", vec![-1.0], vec![0.0], vec![f64::INFINITY]);
        use crate::simplex::{SimplexConfig, SimplexSolver};
        let s = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(s.status, SimplexStatus::Unbounded);
    }

    #[test]
    fn bound_tightening_chain() {
        // min x + y s.t. x + y >= 4, x + 2y >= 6, x <= 10, y <= 10:
        // tightening: from rows, x >= ? — mainly exercises propagation, obj=4 at (4,0)?? check: (0,3)=3 < 4 violates row1? 0+3=3<4 violates. (4,0): row2 4>=6 no. (2,2): 4>=4 ok, 6>=6 ok -> 4.
        let lp = lp_of(
            &[&[1.0, 1.0], &[1.0, 2.0]],
            vec![4.0, 6.0],
            "GG",
            vec![1.0, 1.0],
            vec![0.0, 0.0],
            vec![10.0, 10.0],
        );
        assert!((check_agree(&lp) - 4.0).abs() < 1e-6);
    }

    #[test]
    fn equality_chain_collapses() {
        // ladder: x1-x2=0, x2-x3=0, x3=5 (E rows), min x1+x2+x3 -> 15
        let lp = lp_of(
            &[&[1.0, -1.0, 0.0], &[0.0, 1.0, -1.0], &[0.0, 0.0, 1.0]],
            vec![0.0, 0.0, 5.0],
            "EEE",
            vec![1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lp) - 15.0).abs() < 1e-6);
    }

    #[test]
    fn equality_chain_with_costs_and_bounds() {
        // x1 = 2*x2 (E), min x1 + x2 s.t. x2 >= 3 -> x2=3, x1=6, obj 9
        let lp = lp_of(
            &[&[1.0, -2.0], &[0.0, 1.0]],
            vec![0.0, 3.0],
            "EG",
            vec![1.0, 1.0],
            vec![0.0, 0.0],
            vec![100.0, 100.0],
        );
        assert!((check_agree(&lp) - 9.0).abs() < 1e-6);
    }

    #[test]
    fn chain_infeasible_via_bounds() {
        // x1 = x2 (E), x1 <= 2 (ub), x2 >= 5 (G row) -> infeasible
        let lp = lp_of(
            &[&[1.0, -1.0], &[0.0, 1.0]],
            vec![0.0, 5.0],
            "EG",
            vec![0.0, 0.0],
            vec![0.0, 0.0],
            vec![2.0, 10.0],
        );
        check_agree(&lp); // both modes must say Infeasible
        use crate::simplex::{SimplexConfig, SimplexSolver};
        let s = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(s.status, SimplexStatus::Infeasible);
    }

    #[test]
    fn long_chain_postsolve_exactness() {
        // chain of 10: x_i - x_{i+1} = 0 for i<10, x_10 = 2.5; min sum -> 25
        let n = 10;
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for i in 0..n - 1 {
            let mut r = vec![0.0; n];
            r[i] = 1.0;
            r[i + 1] = -1.0;
            rows.push(r);
        }
        let mut last = vec![0.0; n];
        last[n - 1] = 1.0;
        rows.push(last);
        let refs: Vec<&[f64]> = rows.iter().map(|r| r.as_slice()).collect();
        let lp = lp_of(
            &refs,
            vec![0.0; n],
            &"E".repeat(n),
            vec![1.0; n],
            vec![0.0; n],
            vec![100.0; n],
        );
        let mut lp = lp;
        // fix x_10 = 2.5 via rhs tweak: last row rhs should be 2.5
        lp.rhs[n - 1] = 2.5;
        assert!((check_agree(&lp) - 25.0).abs() < 1e-6);
    }

    #[test]
    fn empty_row_dropped() {
        // min x s.t. 0*x >= -5 (satisfied empty row), x >= 2
        let lp = lp_of(
            &[&[0.0], &[1.0]],
            vec![-5.0, 2.0],
            "GG",
            vec![1.0],
            vec![0.0],
            vec![5.0],
        );
        assert!((check_agree(&lp) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn pure_bound_row_folded_and_dropped() {
        // extractVarBounds port: min -x3 s.t. x1 + x2 + x3 <= 7 with
        // x1, x2 fixed at 1 ⇒ the row is a pure bound x3 <= 5 (x3 starts
        // with an infinite upper bound, which interval tightening cannot
        // handle while fixed columns poison its activity sums). The fold
        // sets ub₃ = 5 and the row is deleted as implied-by-bounds.
        let lp = lp_of(
            &[&[1.0, 1.0, 1.0]],
            vec![7.0],
            "L",
            vec![0.0, 0.0, -1.0],
            vec![1.0, 1.0, 0.0],
            vec![1.0, 1.0, f64::INFINITY],
        );
        assert!((check_agree(&lp) + 5.0).abs() < 1e-6);
        // direct: the fold + drop must collapse the LP to 0 rows (x3 then
        // goes through the empty-column rule at its new bound)
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => assert_eq!(w.num_row, 0, "pure-bound row survived"),
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    /// Chain records only (tests ignore forcing records).
    fn chain_recs(p: &Postsolve) -> Vec<&ChainRec> {
        p.recs
            .iter()
            .filter_map(|r| match r {
                PostRec::Chain(c) => Some(c),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn implied_equality_conversion_feeds_doubleton_chain() {
        // Implied-equality port (loose tier): min x1 + x2 s.t.
        // x1 + x2 >= 4, x1 ∈ (-∞, 2+5e-10], x2 ∈ [0, 2]. Activity max
        // 4+5e-10 touches the far side within the 1e-9·rhs slack but NOT
        // within the strict at_rhs forcing tolerance ⇒ the row converts to
        // an equality, becomes a doubleton the chain rule eliminates
        // (x1 = 4 − x2), and the postsolve record keeps the ORIGINAL 'G'
        // sense. True optimum 4; the conversion can move it by ≤ slack.
        // (The x2 − x1 ≤ 3 row de-parallelizes the two columns for the
        // duplicate-column merge — non-redundant by the −inf lower bound, and slack at the
        // optimum, where it tightens x2 ≤ 3.5 and is then deleted.)
        let lp = lp_of(
            &[&[1.0, 1.0], &[-1.0, 1.0]],
            vec![4.0, 3.0],
            "GL",
            vec![1.0, 1.0],
            vec![f64::NEG_INFINITY, 0.0],
            vec![2.0 + 5e-10, 2.0],
        );
        assert!((check_agree(&lp) - 4.0).abs() < 1e-6);
        match run(&lp, 50, true, false, true, false, false) {
            Presolved::Reduced {
                lp: w, postsolve, ..
            } => {
                assert_eq!(w.num_row, 0, "converted row should be chain-eliminated");
                let chains = chain_recs(&postsolve);
                assert_eq!(chains.len(), 1);
                assert_eq!(
                    chains[0].orig_sense, 'G',
                    "postsolve must remember the pre-conversion sense"
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    // ---- doubleton-equation choice (HiGHS doubletonEq port) ----
    //
    // Shared crafted LP. Row 0 is the doubleton: 10·x0 + x1 = 20 with
    // cnt(x0)=2 (rows 1,2) and cnt(x1)=1 (row 3), |coef| ratio 10 ≫ 2.
    // v1 (pure fill heuristic) substituted the SPARSER small-coefficient
    // side x1; HiGHS colAtPos1Better :3093–3096 substitutes the LARGER
    // |coef| side x0 instead, and the stay-col bounds (:3141–3174) then
    // tighten x1's upper bound 100 → 20 via x1 = (20 − 10·x0)/1 with
    // x0 ≥ 0.
    fn dte_v2_lp() -> Lp {
        lp_of(
            &[
                &[10.0, 1.0, 0.0, 0.0],
                &[1.0, 0.0, 1.0, 0.0],
                &[1.0, 0.0, 0.0, 1.0],
                &[0.0, 1.0, 0.0, 1.0],
            ],
            vec![20.0, 5.0, 5.0, 1.0],
            "ELLG",
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![10.0, 100.0, 10.0, 10.0],
        )
    }

    #[test]
    fn dte_v2_coef_aware_choice_beats_sparser_side() {
        // (a) ratio > 2 ⇒ the larger-|coef| column is substituted even
        // though the other side is sparser (v1 would have dropped x1).
        match run(&dte_v2_lp(), 50, false, false, true, false, false) {
            Presolved::Reduced { live, .. } => {
                assert!(
                    !live.contains(&0),
                    "x0 (|coef| 10) must be the subst column, live={live:?}"
                );
                assert!(
                    live.contains(&1),
                    "x1 (sparser, |coef| 1) must stay, live={live:?}"
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn dte_v2_stay_col_bounds_tightened_beyond_v1() {
        // (b) the stay column receives the sign-split implied bound:
        // x1 = 20 − 10·x0 with x0 ≥ 0 ⇒ ub(x1) = 20 (was 100). Under
        // v1's choice x1 was substituted away and x0 (already tight) stayed,
        // so no bound anywhere tightened — this tightening is the strength
        // delta of the HiGHS choice + :3141–3174 propagation.
        match run(&dte_v2_lp(), 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, live, .. } => {
                let pos = live.iter().position(|&o| o == 1).expect("x1 live");
                assert!(
                    (w.upper[pos] - 20.0).abs() <= 1e-9,
                    "stay-col ub must tighten to 20, got {:?}",
                    w.upper[pos]
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn dte_v2_postsolve_reconstructs() {
        // (c) with the substitution direction flipped vs v1, the ChainRec
        // replay (x0 = (20 − x1)/10) must still reconstruct a feasible,
        // optimal x in the ORIGINAL space (check_agree pattern).
        // True optimum: x2 = x3 = 0, x1 = 1 (row 3), x0 = 1.9 ⇒ obj 0.
        assert!(check_agree(&dte_v2_lp()).abs() < 1e-6);
    }

    #[test]
    fn far_side_pinned_row_forces_and_drops() {
        // Implied-equality port (strict tier): min x1 + x2 s.t.
        // x1 + x2 >= 4 with x1, x2 ≤ 2 — activity max pins at rhs within
        // the forcing tolerance, so the row FORCES both variables to their
        // upper bounds and is deleted (no equality, no chain record:
        // HiGHS checkForcingRow behavior, scrs8 evidence).
        let lp = lp_of(
            &[&[1.0, 1.0]],
            vec![4.0],
            "G",
            vec![1.0, 1.0],
            vec![0.0, 0.0],
            vec![2.0, 2.0],
        );
        assert!((check_agree(&lp) - 4.0).abs() < 1e-6);
        match run(&lp, 50, true, false, true, false, false) {
            Presolved::Reduced {
                lp: w, postsolve, ..
            } => {
                assert_eq!(w.num_row, 0, "forcing row must be deleted, not converted");
                assert!(
                    chain_recs(&postsolve).is_empty(),
                    "forcing fixes variables, no chain elimination"
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    // ---- dominated-col detection probe (HiGHS detectDominatedCol) ----

    /// Forcing records only (tests ignore chain records).
    fn forcing_recs(p: &Postsolve) -> Vec<&ForcingRec> {
        p.recs
            .iter()
            .filter_map(|r| match r {
                PostRec::Forcing(f) => Some(f),
                _ => None,
            })
            .collect()
    }

    // ---- forcing-column removal (A&A 1995 / HiGHS port) ----

    /// min y s.t. x0 + y ≤ 5, x0 ∈ (−∞, 10] cost 0, y ∈ [2, 10] cost 1.
    /// x0 is the forcing column: 'L' row with a=+1 ⇒ m_sum = 0 ⇒
    /// r_lower = cost − 0 = 0 ⇒ dir +1 candidate at zero cost.
    fn forcing_zero_cost_lp() -> Lp {
        lp_of(
            &[&[1.0, 1.0]],
            vec![5.0],
            "L",
            vec![0.0, 1.0],
            vec![f64::NEG_INFINITY, 2.0],
            vec![10.0, 10.0],
        )
    }

    #[test]
    fn forcing_col_fires_and_reconstructs() {
        let lp = forcing_zero_cost_lp();
        // Detection: exactly one candidate, dir +1.
        assert_eq!(
            forcing_col_candidates(&Sp::from_lp(&lp), DUAL_TOL),
            vec![(0, 1)]
        );
        // Presolve-on obj == presolve-off obj, and the postsolve x satisfies
        // the DELETED row x0 + y ≤ 5 (check_agree loops all original rows).
        assert!((check_agree(&lp) - 2.0).abs() < 1e-6);
        // Record + reconstruction, directly: the row deletion leaves y an
        // empty column fixed at 2, and restore must choose
        // x0 = min(finite ub 10, (5 − 2)/1) = 3.
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced {
                lp: w,
                live,
                postsolve,
            } => {
                assert!(!live.contains(&0), "forcing column must leave live");
                assert_eq!(w.num_row, 0, "the deleted row must leave work");
                let f = forcing_recs(&postsolve);
                assert_eq!(f.len(), 1);
                assert_eq!((f[0].orig, f[0].dir, f[0].finite_bound), (0, 1, 10.0));
                let (x, obj) = restore(&postsolve, &live, &[], 0.0);
                assert_eq!(x, vec![3.0, 2.0], "far-side-overlap x0 = 3");
                assert!((obj - 2.0).abs() < 1e-12);
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn forcing_near_miss_not_fired() {
        // Same LP, cost −10·DUAL_TOL on x0 (dynamism 1): the gate needs
        // |cost − m_sum| ≤ DUAL_TOL/1 and this is 10× that ⇒ no candidate,
        // no mutation. The tiny negative cost pushes x0 up to the row cap:
        // optimum x0 = 5 − 2 = 3, y = 2, obj = 2 − 3·10·DUAL_TOL.
        let lp = lp_of(
            &[&[1.0, 1.0]],
            vec![5.0],
            "L",
            vec![-10.0 * DUAL_TOL, 1.0],
            vec![f64::NEG_INFINITY, 2.0],
            vec![10.0, 10.0],
        );
        assert!(
            forcing_col_candidates(&Sp::from_lp(&lp), DUAL_TOL).is_empty(),
            "10× gate cost must not fire"
        );
        let obj = check_agree(&lp);
        assert!((obj - (2.0 - 3.0e-6)).abs() < 1e-6);
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced {
                live, postsolve, ..
            } => {
                assert!(live.contains(&0), "near-miss column must stay live");
                assert!(forcing_recs(&postsolve).is_empty());
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn forcing_mixed_sign_rows_not_fired() {
        // x0 ∈ (−∞, 10] cost 0 sits in one L row with a>0 AND one G row
        // with a>0: the G row's dual seed [0,+∞) makes m_sum = +∞ ⇒
        // r_lower = −∞ ⇒ no candidate (soundness note: mixed membership
        // with the wrong signs cannot pass the ≈0 gate). Optimum y=0 at
        // x0 ∈ [1, 5]; true obj 0.
        let lp = lp_of(
            &[&[1.0, 1.0], &[1.0, 1.0]],
            vec![5.0, 1.0],
            "LG",
            vec![0.0, 1.0],
            vec![f64::NEG_INFINITY, 0.0],
            vec![10.0, 10.0],
        );
        assert!(
            forcing_col_candidates(&Sp::from_lp(&lp), DUAL_TOL).is_empty(),
            "mixed-sign column must not fire"
        );
        assert!(check_agree(&lp).abs() < 1e-6);
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced {
                live, postsolve, ..
            } => {
                assert!(live.contains(&0));
                assert!(forcing_recs(&postsolve).is_empty());
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn strong_dominated_col_detected() {
        // min 100·x s.t. x ≤ 1, x ∈ [0, 10]: the 'L' row dual box is
        // (−∞, 0], a = +1 ⇒ max(a·π) = 0 ⇒ r_lower = 100 − 0 = 100 >
        // DUAL_TOL ⇒ fix-to-lb candidate (+1); margin ≫ 10·DUAL_TOL so no
        // near-miss, finite lb so no dual-infeasible signal.
        let lp = lp_of(
            &[&[1.0]],
            vec![1.0],
            "L",
            vec![100.0],
            vec![0.0],
            vec![10.0],
        );
        let (cands, near, dual_inf) = dominated_col_candidates(&Sp::from_lp(&lp), DUAL_TOL);
        assert_eq!(cands, vec![(0, 1)], "x must be a fix-to-lb candidate");
        assert_eq!(near, 0);
        assert_eq!(dual_inf, 0);
        // detection-only: the probe mutates nothing, solve behavior unchanged
        assert!((check_agree(&lp) - 0.0).abs() < 1e-6);
    }

    #[test]
    fn near_miss_not_fired() {
        // Same LP with cost 0 ⇒ r_lower = 0, not > DUAL_TOL: no candidate
        // at all (0 does not even reach the near-miss band [tol, 10·tol]).
        let lp = lp_of(&[&[1.0]], vec![1.0], "L", vec![0.0], vec![0.0], vec![10.0]);
        let (cands, near, dual_inf) = dominated_col_candidates(&Sp::from_lp(&lp), DUAL_TOL);
        assert!(cands.is_empty(), "zero reduced cost must not fire");
        assert_eq!(near, 0);
        assert_eq!(dual_inf, 0);
        check_agree(&lp);
    }

    #[test]
    fn free_row_dual_blocks() {
        // min 100·x s.t. x ≥ 1: 'G' row dual box is [0, +∞), a = +1 ⇒
        // r_lower = 100 − 1·(+∞) = −∞ ⇒ NO candidate despite the cost
        // looking dominant (proves the ±∞ interval guard works).
        let lp = lp_of(
            &[&[1.0]],
            vec![1.0],
            "G",
            vec![100.0],
            vec![0.0],
            vec![10.0],
        );
        let (cands, near, dual_inf) = dominated_col_candidates(&Sp::from_lp(&lp), DUAL_TOL);
        assert!(cands.is_empty(), "unbounded row dual must block detection");
        assert_eq!(near, 0);
        assert_eq!(dual_inf, 0);
        check_agree(&lp);
    }

    // ---- parallel row / duplicate col detection probe ----

    #[test]
    fn parallel_rows_detected_not_mutated() {
        // Two proportional 'L' rows: x0+x1 ≤ 5 and 2x0+2x1 ≤ 8 normalize to
        // the same coefficient vector (scales +1 / +2) ⇒ one group of 2.
        // Detection alone keeps BOTH rows in the reduced LP and the solve
        // is unchanged (optimum x0 = 0 at cost 1·x0 with x ∈ [0,10]).
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0], &[2.0, 2.0, 0.0]],
            vec![5.0, 8.0],
            "LL",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        let groups = parallel_row_groups(&Sp::from_lp(&lp));
        assert_eq!(groups, vec![vec![0, 1]], "proportional L rows must group");
        assert!(
            parallel_col_groups(&Sp::from_lp(&lp)).is_empty(),
            "no parallel cols here"
        );
        assert!(check_agree(&lp).abs() < 1e-6);
        // (The default path MERGES this pair — the detection-only
        // "reduced keeps both rows" pin moved to `parrow_near_miss_untouched`;
        // here detection is asserted to be a pure function.)
        assert_eq!(
            parallel_row_groups(&Sp::from_lp(&lp)),
            groups,
            "detection must be pure"
        );
    }

    #[test]
    fn parallel_rows_near_miss_and_support_mismatch_not_grouped() {
        // Near-miss: 2x0+3x1 ≤ 8 vs x0+x1 ≤ 5 ⇒ NOT proportional (normalized
        // second coefficient 1.5 ≠ 1) ⇒ no group.
        let lp = lp_of(
            &[&[1.0, 1.0, 1.0], &[2.0, 3.0, 0.0]],
            vec![5.0, 8.0],
            "LL",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!(
            parallel_row_groups(&Sp::from_lp(&lp)).is_empty(),
            "non-proportional rows must not group"
        );
        // Cross-sense parallel with the SAME support is now GROUPED (the
        // interval arms need cross-sense witnesses), but a true
        // range like [4, 5] keeps both witnesses, so nothing is deleted.
        // This lp2 pair has DIFFERENT support (row 1 lacks x2) ⇒ not even
        // grouped.
        let lp2 = lp_of(
            &[&[1.0, 1.0, 1.0], &[2.0, 2.0, 0.0]],
            vec![5.0, 8.0],
            "LG",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!(
            parallel_row_groups(&Sp::from_lp(&lp2)).is_empty(),
            "different-support parallels must not group"
        );
        // Same-support true range [4, 5]: grouped, but both witnesses kept.
        let lp3 = lp_of(
            &[&[1.0, 1.0, 0.0], &[2.0, 2.0, 0.0]],
            vec![5.0, 8.0],
            "LG",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert_eq!(
            parallel_row_groups(&Sp::from_lp(&lp3)).len(),
            1,
            "cross-sense same-support parallels group"
        );
        check_agree(&lp2);
        check_agree(&lp3);
    }

    #[test]
    fn duplicate_cols_detected_cost_class() {
        // Columns x0 = x1 = (1, 2)ᵀ with equal costs ⇒ one group of 2
        // (scale +1). x2 = (1, 1)ᵀ stays out.
        let lp = lp_of(
            &[&[1.0, 1.0, 1.0], &[2.0, 2.0, 1.0]],
            vec![6.0, 10.0],
            "LL",
            vec![1.0, 1.0, -1.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert_eq!(
            parallel_col_groups(&Sp::from_lp(&lp)),
            vec![vec![0, 1]],
            "equal cols with equal cost must group"
        );
        // Cost differing by 10·DUAL_TOL (beyond the objDiff gate) ⇒ no group.
        let lp2 = lp_of(
            &[&[1.0, 1.0, 1.0], &[2.0, 2.0, 1.0]],
            vec![6.0, 10.0],
            "LL",
            vec![1.0, 1.0 + 10.0 * DUAL_TOL, -1.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!(
            parallel_col_groups(&Sp::from_lp(&lp2)).is_empty(),
            "cost-different parallel cols must not group"
        );
        // Anti-parallel (x1 = −x0 coefficients, cost negated) normalizes to
        // the same vector and passes the objDiff gate (c_r·s = 1·(−1) = −1 =
        // c_1) ⇒ groups with relative scale −1.
        let lp3 = lp_of(
            &[&[1.0, -1.0, 1.0], &[2.0, -2.0, 1.0]],
            vec![6.0, 10.0],
            "LL",
            vec![1.0, -1.0, -1.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert_eq!(
            parallel_col_groups(&Sp::from_lp(&lp3)),
            vec![vec![0, 1]],
            "anti-parallel cost-mirrored cols must group"
        );
        // Agreement + postsolve-feasibility both ways (optimum not asserted:
        // detection alone leaves the solve path untouched).
        check_agree(&lp);
        check_agree(&lp3);
    }

    // ---- parallel-row merge (mutation + dominance) ----

    #[test]
    fn parrow_dominant_row_survives_loser_deleted() {
        // 'L' pair x0+x1 ≤ 5 (eff 5) vs 2x0+2x1 ≤ 8 (eff 4): the SECOND row
        // dominates ⇒ row 0 is deleted with its own rhs (8), no rewrite.
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0], &[2.0, 2.0, 0.0]],
            vec![5.0, 8.0],
            "LL",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lp) - 0.0).abs() < 1e-6);
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced {
                lp: w, postsolve, ..
            } => {
                assert_eq!(w.num_row, 1, "the dominated row must be deleted");
                assert_eq!(
                    w.rhs[0], 8.0,
                    "the DOMINANT row (eff 4 = 8/2) survives as-is"
                );
                assert!(
                    postsolve.recs.is_empty(),
                    "row merge needs NO postsolve record (DuplicateRow::undo is duals-only)"
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
        // 'E' pure duplicates (eff 5 = 5): one copy deleted, the other kept.
        // Three nonzeros with DISTINCT costs so neither the DTE chain (needs
        // a 2-nz row) nor col merge (needs cost-equivalent cols)
        // can pre-empt the row-count assertion. Optimum: x0 = 5, obj 5.
        let lpe = lp_of(
            &[&[1.0, 1.0, 1.0, 0.0], &[2.0, 2.0, 2.0, 0.0]],
            vec![5.0, 10.0],
            "EE",
            vec![1.0, 2.0, 3.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lpe) - 5.0).abs() < 1e-6);
        match run(&lpe, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                assert_eq!(w.num_row, 1, "duplicate E row must be deleted");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parrow_negative_scale_bounds_below_max_eff_wins() {
        // e226 regression (rows 49/50 there): 'L' rows with NEGATIVE
        // coefficients bound the normalized activity from BELOW —
        // −(x0+x1) ≤ −5 ⇔ x0+x1 ≥ 5 DOMINATES −2(x0+x1) ≤ −8 ⇔ ≥ 4.
        // Direction-blind min-eff (the first cut) deleted the dominant row
        // and broke postsolve feasibility; max-eff must keep row 0 (rhs −5).
        let lp = lp_of(
            &[&[-1.0, -1.0, 0.0], &[-2.0, -2.0, 0.0]],
            vec![-5.0, -8.0],
            "LL",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lp) - 0.0).abs() < 1e-6);
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                assert_eq!(w.num_row, 1, "the dominated row must be deleted");
                assert_eq!(w.rhs[0], -5.0, "max-eff (lower-bound) row survives");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parrow_mixed_effective_direction_skipped() {
        // Same raw sense 'L' + same normalized vector, but OPPOSITE scales:
        // x0+x1 ≤ 5 (upper eff 5) vs −(x0+x1) ≤ −2 (lower eff 2). The pair
        // is a slab 2 ≤ x0+x1 ≤ 5 — HiGHS intersects it into a ranged row;
        // we skip: both rows must survive and the slab solve must agree.
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0], &[-1.0, -1.0, 0.0]],
            vec![5.0, -2.0],
            "LL",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert_eq!(
            parallel_row_groups(&Sp::from_lp(&lp)),
            vec![vec![0, 1]],
            "mixed-direction rows ARE detected as parallel"
        );
        assert!((check_agree(&lp) - 0.0).abs() < 1e-6);
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                assert_eq!(w.num_row, 2, "mixed-direction group must be skipped");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parrow_near_miss_untouched() {
        // 2x0+3x1 ≤ 8 vs x0+x1 ≤ 5: not proportional ⇒ no merge, both rows
        // survive to the simplex.
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0], &[2.0, 3.0, 0.0]],
            vec![5.0, 8.0],
            "LL",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lp) - 0.0).abs() < 1e-6);
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                assert_eq!(w.num_row, 2, "near-miss rows must both survive");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parrow_cross_sense_and_conflicting_e_skipped() {
        // Cross-sense parallel ('L' eff 5 vs 'G' eff 4): feasible slab
        // 4 ≤ x0+x1 ≤ 5 — NOT merged (bound intersection not ported), both
        // rows survive and the simplex solves it.
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0], &[2.0, 2.0, 0.0]],
            vec![5.0, 8.0],
            "LG",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lp) - 0.0).abs() < 1e-6);
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                assert_eq!(w.num_row, 2, "cross-sense parallel must be skipped");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
        // Conflicting scaled 'E' pair (eff 5 vs eff 4): genuinely infeasible
        // (x0+x1 = 5 AND = 4). NOT merged (conflicting E skipped by design);
        // the infeasibility surfaces downstream — interval tightening
        // collapses the box until an activity check proves it (the simplex
        // would otherwise). Presolve on and off must agree.
        let lp2 = lp_of(
            &[&[1.0, 1.0, 0.0], &[2.0, 2.0, 0.0]],
            vec![5.0, 8.0],
            "EE",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        check_agree(&lp2);
        use crate::simplex::{SimplexConfig, SimplexSolver};
        let s = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp2)
            .unwrap();
        assert_eq!(s.status, SimplexStatus::Infeasible);
    }

    #[test]
    fn parrow_cross_sense_interval_collapse_to_equality() {
        // Graph20 class: 2·(x0+x1) ≥ 4 ∧ x0+x1 ≤ 2 ⇒ the
        // effective interval for v = x0+x1 is exactly [2, 2] ⇒ the
        // conjunction IS v = 2: the LOWER witness converts to 'E', the upper
        // side is deleted. min −x0 then forces x0 = 2.
        let lp = lp_of(
            &[&[2.0, 2.0, 0.0], &[1.0, 1.0, 0.0]],
            vec![4.0, 2.0],
            "GL",
            vec![-1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lp) + 2.0).abs() < 1e-6, "obj must be −2");
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                // The collapse feeds the singleton-equality substitution:
                // one 'E' with two nonzeros is solved out entirely.
                assert!(w.num_row <= 1, "collapse must reduce to ≤ 1 row");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parrow_equality_inside_range_wins() {
        // HiGHS duplicateRow E arm: {E eff 4, L eff 5, G eff 3} —
        // the equality agrees with itself and lies inside [3, 5] ⇒ keep the
        // 'E', delete both sided witnesses. min −x0 ⇒ x0 = 4.
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0], &[1.0, 1.0, 0.0], &[1.0, 1.0, 0.0]],
            vec![5.0, 3.0, 4.0],
            "LGE",
            vec![-1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lp) + 4.0).abs() < 1e-6, "obj must be −4");
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                // Same as above: the surviving 'E' is a two-nonzeros
                // singleton and gets substituted out; the point is the two
                // sided witnesses are gone.
                assert!(w.num_row <= 1, "only the equality can survive");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parrow_negative_scale_lower_bound_wins() {
        // e226 class: x0+x1 ≥ 0.3467 ('G', scale +1 ⇒ lower
        // eff 0.3467) vs −x0−x1 ≤ −0.9179 ('L', scale −1 ⇒ dividing flips
        // the sense ⇒ lower eff 0.9179). Same side: the NEGATIVE-scale row
        // dominates and the plain 'G' row goes (the earlier sign-blind winner
        // pick deleted the wrong one and broke e226 postsolve).
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0], &[-1.0, -1.0, 0.0]],
            vec![0.3467, -0.9179],
            "GL",
            vec![-1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        );
        assert!((check_agree(&lp) + 10.0).abs() < 1e-6, "obj must be −10");
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                assert_eq!(w.num_row, 1, "dominated same-side row deleted");
                assert_eq!(w.senses[0], 'L', "the negative-scale row kept");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    // ---- forcing-conflict safety ----

    #[test]
    fn forcing_conflict_same_pass_reports_infeasible() {
        // Two E rows scanned in ONE pass pin x from the same bounds
        // snapshot to DIFFERENT extremes: row 0 (-x+y+p = -1) reaches its
        // activity minimum only at x = ub = 1, row 1 (x+z+w = 0) reaches
        // its minimum only at x = lb = 0. Both forcings are valid
        // deductions, so the LP is infeasible. Before the conflict check
        // the later fix was silently clamped into the earlier fix's box
        // while BOTH rows had already been dropped — a false Optimal.
        let lp = lp_of(
            &[&[-1.0, 1.0, 1.0, 0.0, 0.0], &[1.0, 0.0, 0.0, 1.0, 1.0]],
            vec![-1.0, 0.0],
            "EE",
            vec![1.0, 1.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0],
            vec![1.0, 1.0, 1.0, 1.0, 1.0],
        );
        check_agree(&lp); // status must match the presolve-off solve
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Terminal {
                status: SimplexStatus::Infeasible,
            } => {}
            other => panic!("expected Terminal(Infeasible), got {other:?}"),
        }
    }

    #[test]
    fn forcing_conflict_set_partition_s2_repro() {
        // Inline copy of benches/data/repros/presolve_false_optimal_sp_s2.mps:
        // an infeasible set-partitioning LP (11 E rows, 28 binary-
        // bounded columns) that presolve used to reduce to 0 rows and
        // report Optimal 107 after 0 pivots; --no-presolve and HiGHS both
        // say Infeasible. The culprit: two forcing rows in one pass pinned
        // a column reached by a negative coefficient (pinned to ub = 1)
        // and by a positive one (pinned to lb = 0); the clamp swallowed it.
        let rows: &[&[f64]] = &[
            &[
                0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 1., 1., 0., 1., 1., 0., 0., 1., 0., 0.,
                0., 0., 0., 0., 0., 0., 0.,
            ],
            &[
                0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
                0., 0., 1., 0., 0., 0., 0.,
            ],
            &[
                1., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 1., 1., 0., 0., 0., 0.,
                1., 0., 0., 1., 0., 0., 0.,
            ],
            &[
                0., 0., 0., 1., 0., 1., 1., 1., 0., 1., 0., 0., 0., 0., 0., 1., 1., 0., 0., 1., 1.,
                0., 0., 0., 0., 0., 0., 0.,
            ],
            &[
                1., 1., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
                1., 0., 0., 0., 1., 0., 1.,
            ],
            &[
                1., 0., 0., 1., 1., 0., 0., 1., 1., 1., 0., 0., 1., 0., 0., 1., 1., 0., 1., 0., 0.,
                0., 0., 0., 0., 0., 1., 1.,
            ],
            &[
                0., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 1., 1., 1., 0., 0., 0.,
                0., 1., 0., 0., 1., 1., 0.,
            ],
            &[
                1., 0., 1., 0., 0., 1., 0., 1., 0., 1., 1., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0.,
                0., 0., 0., 0., 0., 0., 1.,
            ],
            &[
                1., 0., 0., 1., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0.,
                0., 0., 1., 0., 1., 1., 0., 1.,
            ],
            &[
                0., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 0., 0.,
                0., 0., 1., 0., 0., 0., 1.,
            ],
            &[
                1., 0., 1., 1., 0., 1., 1., 1., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 1.,
                0., 0., 0., 0., 1., 0., 1.,
            ],
        ];
        let obj = vec![
            23., 10., 9., 60., 17., 6., 40., 14., 15., 6., 20., 25., 8., 18., 9., 24., 36., 5.,
            10., 23., 20., 13., 14., 14., 19., 9., 19., 8.,
        ];
        let lo = vec![
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
            0., 0., 0., 0., 0., 0.,
        ];
        let up = vec![
            0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 1., 1., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1.,
            1., 0., 0., 0., 1., 0.,
        ];
        let lp = lp_of(rows, vec![1.0; 11], &"E".repeat(11), obj, lo, up);
        check_agree(&lp); // presolve on and off must now agree on Infeasible
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Terminal {
                status: SimplexStatus::Infeasible,
            } => {}
            other => panic!("expected Terminal(Infeasible), got {other:?}"),
        }
    }

    // ---- parallel-row merge under tightened/fixed shared columns
    // (MIP node LP shape) ----

    #[test]
    fn parrow_tightened_shared_col_effective_rhs() {
        // MIP node LP shape: two parallel 'L' rows with NEGATIVE scale
        // and a shared column fixed by branching (x = 0.65). Each row's
        // signed max-coef is negative, so dividing flips 'L' into a LOWER
        // bound on the normalized activity: the rows read x+y+w >= 1.0
        // (loser) and x+y+w >= 1.3 (winner). w is free, so interval
        // tightening skips both rows and the merge's effective-rhs (rhs /
        // signed-max-coef) comparison is what decides: max-eff wins ⇒
        // winner kept. A raw/min-eff pick (the earlier version) keeps the loser,
        // and the postsolve x (y=0, w=0.35) violates the deleted winner by
        // 0.3 — caught by check_agree's original-row feasibility pass.
        let lp = lp_of(
            &[&[-1.0, -1.0, -1.0], &[-2.0, -2.0, -2.0]],
            vec![-1.0, -2.6],
            "LL",
            vec![0.0, 1.0, 0.0],
            vec![0.65, 0.0, f64::NEG_INFINITY],
            vec![0.65, 2.0, f64::INFINITY],
        );
        check_agree(&lp);
    }

    #[test]
    fn parrow_parallel_equalities_fixed_shared_col() {
        // Parallel 'E' rows whose RAW rhs differ (1.65 vs 3.3) but whose
        // effective rhs agree exactly after scaling by the signed max-coef
        // (1.65 = 3.3/2): pure duplicates ⇒ mergeable even though the
        // shared columns x, z are fixed (MIP node fixings) and their folds
        // land in the rhs only later in the pass. Postsolve x must satisfy
        // BOTH original rows.
        let lp = lp_of(
            &[&[1.0, 1.0, 1.0], &[2.0, 2.0, 2.0]],
            vec![1.65, 3.3],
            "EE",
            vec![0.0, 1.0, 0.0],
            vec![0.65, -1.0, 1.0],
            vec![0.65, 1.0, 1.0],
        );
        let obj = check_agree(&lp);
        assert!(obj.abs() < 1e-6, "y must be pinned at 0, got {obj}");
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                assert!(w.num_row <= 1, "duplicate E must not survive twice");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parrow_shared_col_fixed_mid_presolve() {
        // The shared column x is fixed DURING presolve (DTE on row 2:
        // x + p = 0 over [0,1] ⇒ x = 0) while rows 0/1 are still live
        // parallel rows; w ∈ [0,+∞) poisons interval tightening (infinite
        // ub ⇒ the row is skipped) so the merge decides. The fixed-column
        // fold subtracts
        // a·v from BOTH parallel rows' rhs proportionally to their scales,
        // so the effective-rhs ordering is fold-invariant: eff 1.2 kept ⇒
        // y* = 1.2, obj −1.2; keeping eff 1.3 instead admits y = 1.3 and
        // postsolve x violates the deleted row by 0.1.
        let lp = lp_of(
            &[
                &[1.0, 1.0, 1.0, 0.0],
                &[2.0, 2.0, 2.0, 0.0],
                &[1.0, 0.0, 0.0, 1.0],
            ],
            vec![1.2, 2.6, 0.0],
            "LLE",
            vec![0.0, -1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![1.0, 2.0, f64::INFINITY, 1.0],
        );
        let obj = check_agree(&lp);
        assert!(
            (obj + 1.2).abs() < 1e-6,
            "kept row must bind: obj {obj} != -1.2"
        );
    }

    // ---- duplicate-column merge (mutation + split reconstruction) ----

    /// DupCol records only (tests ignore other record kinds).
    fn dupcol_recs(p: &Postsolve) -> Vec<&DupColRec> {
        p.recs
            .iter()
            .filter_map(|r| match r {
                PostRec::DupCol(d) => Some(d),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn parcol_duplicate_merged_and_split_reconstructed() {
        // x0 = x1 = (1,2)ᵀ cost 1/1, x0 ∈ [0,3], x1 ∈ [0,4], x2 cost −1.
        // Rows force m = x0+x1 = 5 exactly (G row ≥ 5, and 2m+x2 ≤ 10 with
        // x2 ≥ 0 caps m at 5), x2 = 0 ⇒ obj 5. The split of m = 5 must give
        // x1 = clamp(5−0, 0, 4) = 4 and x0 = 1 — both strictly legal, and
        // check_agree re-feasibilizes every original row.
        let lp = lp_of(
            &[&[1.0, 1.0, 1.0], &[2.0, 2.0, 1.0], &[1.0, 1.0, 0.0]],
            vec![6.0, 10.0, 5.0],
            "LLG",
            vec![1.0, 1.0, -1.0],
            vec![0.0, 0.0, 0.0],
            vec![3.0, 4.0, 10.0],
        );
        assert!((check_agree(&lp) - 5.0).abs() < 1e-6);
        match run(&lp, 50, false, true, true, false, false) {
            Presolved::Reduced {
                lp: w,
                live,
                postsolve,
            } => {
                assert!(!live.contains(&1), "duplicate must leave live");
                let d = dupcol_recs(&postsolve);
                assert_eq!(d.len(), 1);
                assert_eq!((d[0].orig, d[0].kept, d[0].scale), (1, 0, 1.0));
                // kept_lo is 1, not 0 (and dup_lo 2, not 0): the G row
                // (m ≥ 5) tightens BOTH lowers before the merge (x0 ≥ 5−4,
                // x1 ≥ 5−3). The record stores merge-time bounds, as the
                // split needs them: m = 5 ⇒ x1 = clamp(5−1, 2, 4) = 4,
                // x0 = 1.
                assert_eq!((d[0].kept_lo, d[0].kept_hi), (1.0, 3.0));
                assert_eq!((d[0].dup_lo, d[0].dup_hi), (2.0, 4.0));
                // The merged variable is FORCED to m = 5 (m ≥ 5 and
                // 2m + x2 ≤ 10 with x2 ≥ 0), so it leaves through the fixed
                // path and the LP collapses to nothing; restore must split
                // m = 5 into x1 = clamp(5−0, 0, 4) = 4, x0 = 1.
                assert!(w.num_row == 0 && live.is_empty(), "LP must collapse");
                let (x, obj) = restore(&postsolve, &live, &[], 0.0);
                assert_eq!(x, vec![1.0, 4.0, 0.0], "merged value split");
                assert!((obj - 5.0).abs() < 1e-12);
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parcol_antiparallel_scale_minus_one() {
        // x1 = −(1,2)ᵀ with mirrored cost −1 (objDiff gate passes at scale
        // −1): merged m = x0 − x1 ∈ [0−4, 3−0]. Optimum m = −4, x2 = 10,
        // obj = m − x2 = −14 ⇒ the split must return x0 = 0, x1 = 4.
        let lp = lp_of(
            &[&[1.0, -1.0, 1.0], &[2.0, -2.0, 1.0]],
            vec![6.0, 10.0],
            "LL",
            vec![1.0, -1.0, -1.0],
            vec![0.0, 0.0, 0.0],
            vec![3.0, 4.0, 10.0],
        );
        assert!((check_agree(&lp) + 14.0).abs() < 1e-6);
        match run(&lp, 50, false, true, true, false, false) {
            Presolved::Reduced {
                live, postsolve, ..
            } => {
                assert!(!live.contains(&1), "duplicate must leave live");
                let d = dupcol_recs(&postsolve);
                assert_eq!(d.len(), 1);
                assert_eq!(d[0].scale, -1.0);
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parcol_general_scale_merged_and_split() {
        // x1 = (2,4)ᵀ = 2·x0 with cost 2 = 2·c0 (objDiff gate
        // passes at scale 2). HiGHS `duplicateColumn` merges any real scale
        // for LP, so this pair MERGES: m = x0 + 2·x1 ∈ [0, 3 + 2·4], rows
        // read m + x2 ≤ {6,10}, obj = m − x2. With x2 ≤ 10 the optimum is
        // m = 0, x2 = 6 ⇒ obj −6; the split of m = 0 gives x1 = 0, x0 = 0.
        let lp = lp_of(
            &[&[1.0, 2.0, 1.0], &[2.0, 4.0, 1.0]],
            vec![6.0, 10.0],
            "LL",
            vec![1.0, 2.0, -1.0],
            vec![0.0, 0.0, 0.0],
            vec![3.0, 4.0, 10.0],
        );
        assert_eq!(
            parallel_col_groups(&Sp::from_lp(&lp)),
            vec![vec![0, 1]],
            "scale-2 pair must be DETECTED"
        );
        assert!((check_agree(&lp) + 6.0).abs() < 1e-6);
        match run(&lp, 50, false, true, true, false, false) {
            Presolved::Reduced {
                live, postsolve, ..
            } => {
                assert!(!live.contains(&1), "duplicate must leave live");
                let d = dupcol_recs(&postsolve);
                assert_eq!(d.len(), 1);
                assert_eq!((d[0].orig, d[0].kept, d[0].scale), (1, 0, 2.0));
                // dup_hi is 2.5, not 4.0: row 1 (4x1 + … ≤ 10 with the
                // lowers at 0) tightens x1's upper BEFORE the merge.
                assert_eq!((d[0].dup_lo, d[0].dup_hi), (0.0, 2.5));
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parcol_split_free_kept_col_uses_finite_anchor() {
        // Regression (agg fuzz seed 90): BOTH merge-time boxes
        // are (−inf, +inf). The naive split (m − kept_lo)/scale evaluates
        // to ±inf for a free kept column and the clamp never fires
        // (inf > inf is false) — ±inf leaks into BOTH halves and the
        // postsolve x is NaN garbage. HiGHS anchors the kept column to
        // min(0, kept_hi) = 0 first. Here m is pinned to 3:
        // m + x2 ≥ 4 ∧ m − x2 ≤ 2 with x2 ≥ 0 forces x2 = 1, m = 3 ⇒
        // split x1 = (3 − 0)/3 = 1, x0 = 0; obj = 1.
        let lp = lp_of(
            &[&[1.0, 3.0, 1.0], &[1.0, 3.0, -1.0]],
            vec![4.0, 2.0],
            "GL",
            vec![0.0, 0.0, 1.0],
            vec![f64::NEG_INFINITY, f64::NEG_INFINITY, 0.0],
            vec![f64::INFINITY, f64::INFINITY, 10.0],
        );
        assert!((check_agree(&lp) - 1.0).abs() < 1e-6, "obj must be 1");
        match run(&lp, 50, false, true, false, false, false) {
            Presolved::Reduced { postsolve, .. } => {
                let d = dupcol_recs(&postsolve);
                assert_eq!(d.len(), 1, "the free pair must merge");
                assert_eq!((d[0].orig, d[0].kept, d[0].scale), (1, 0, 3.0));
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parcol_split_negative_scale_free_cols_finite() {
        // Same free-pair shape at NEGATIVE scale: x1 = (−2,−2)ᵀ = −2·x0
        // (objDiff gate 0·(−2) − 0 = 0 passes). m = x0 − 2x1: the rows
        // 2 + x2 ≤ m ≤ 6 − x2 squeeze x2 ≤ 2 ⇒ optimum x2 = 2, m = 4,
        // obj −2; the split anchors x0 at 0 and gives x1 = (4 − 0)/(−2)
        // = −2 (finite; the naive form leaks ∓inf). check_agree asserts
        // ON == OFF and postsolve-x feasibility in the ORIGINAL rows.
        let lp = lp_of(
            &[&[1.0, -2.0, -1.0], &[1.0, -2.0, 1.0]],
            vec![2.0, 6.0],
            "GL",
            vec![0.0, 0.0, -1.0],
            vec![f64::NEG_INFINITY, f64::NEG_INFINITY, 0.0],
            vec![f64::INFINITY, f64::INFINITY, 10.0],
        );
        assert!((check_agree(&lp) + 2.0).abs() < 1e-6, "obj must be −2");
        match run(&lp, 50, false, true, false, false, false) {
            Presolved::Reduced { postsolve, .. } => {
                let d = dupcol_recs(&postsolve);
                assert_eq!(d.len(), 1, "the negative-scale free pair merges");
                assert_eq!((d[0].orig, d[0].kept, d[0].scale), (1, 0, -2.0));
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn memguard_decision_matches_measured_ladder() {
        // Threshold fit (measured on our benchmark set): NO-OP at every
        // measured instance — all of Netlib, the 47 relaxed MIPLIB LPs, our
        // MIPLIB set and the big-nnz large-LP points whose presolve is measured cheap AND valuable.
        // netlib max nnz (scsd6):
        assert!(!memguard_skip(4_316, PS_GUARD_BUDGET_BYTES));
        // relaxed-MIPLIB max nnz (neos-860300):
        assert!(!memguard_skip(384_329, PS_GUARD_BUDGET_BYTES));
        // MIPLIB s100 root LP:
        assert!(!memguard_skip(1_777_917, PS_GUARD_BUDGET_BYTES));
        // worst measured B/nnz outlier (bdry2, 2.32 GB actual):
        assert!(!memguard_skip(1_500_003, PS_GUARD_BUDGET_BYTES));
        // measured-cheap + valuable big-nnz points that MUST be kept:
        // scpm1 (1.38 GB, zero reductions), square41 (1.87 GB, rows
        // −96%), a2864 (3.01 GB, cols −93%):
        assert!(!memguard_skip(6_250_000, PS_GUARD_BUDGET_BYTES));
        assert!(!memguard_skip(13_566_426, PS_GUARD_BUDGET_BYTES));
        assert!(!memguard_skip(20_078_717, PS_GUARD_BUDGET_BYTES));
        // beyond the corpus (max measured nnz = a2864's 20.08M) — the
        // guard's preventive arm, no current instance reaches it:
        assert!(memguard_skip(33_000_000, PS_GUARD_BUDGET_BYTES));
        // tpl-tub-ws1617: 1.8 GB FILE but sparse (nnz 4.72M, rows 1.15M)
        // — correctly KEPT (presolve removes 72 % of rows, 2.4 GB):
        assert!(!memguard_skip(4_720_567, PS_GUARD_BUDGET_BYTES));
        // disabled budget never skips; a tiny budget always does:
        assert!(!memguard_skip(33_000_000, 0));
        assert!(memguard_skip(1, 1));
    }

    #[test]
    fn memguard_skip_is_identity_through_run_and_restore() {
        // min x0 + x1 s.t. x0 + x1 >= 2, x0 fixed at 1 (lb == ub) — the
        // fixed-variable rule normally deletes a column, so the
        // not-skipped arm proves the LP is otherwise reducible.
        let lp = lp_of(
            &[&[1.0, 1.0]],
            vec![2.0],
            "G",
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 5.0],
        );
        // Guard armed (budget 1 byte): identity hand-back.
        match run_with_budget(&lp, 50, false, false, true, false, false, 1) {
            Presolved::Reduced {
                lp: w,
                live,
                postsolve,
            } => {
                assert_eq!((w.num_row, w.num_col), (lp.num_row, lp.num_col));
                assert_eq!(w.a.nnz(), lp.a.nnz());
                assert_eq!(live, (0..lp.num_col).collect::<Vec<_>>());
                assert_eq!(postsolve.recs.len(), 0);
                assert_eq!(postsolve.obj_const, 0.0);
                // restore degenerates to identity: x and obj pass through.
                let rx = vec![1.0, 3.0];
                let (x, obj) = restore(&postsolve, &live, &rx, 7.5);
                assert_eq!(x, rx);
                assert_eq!(obj, 7.5);
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
        // Guard off (default budget): the fixed column IS eliminated.
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: w, .. } => {
                assert!(
                    w.num_col < lp.num_col,
                    "fixed x0 must be substituted away (got {} cols)",
                    w.num_col
                )
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn parcol_general_scale_split_clamps_at_upper() {
        // Scale 2 with the objective driving the merged value to the TOP of
        // the box so the split must clamp x1 at its upper bound:
        // min −(x0 + 2x1) s.t. x0 + 2x1 ≤ 10, x0 ∈ [0,3], x1 ∈ [0,4].
        // Merged box [0, 3+2·4] = [0,10]; optimum m = 10 ⇒ split
        // x1 = clamp(10/2, 0, 4) = 4 (upper clamp FIRES), x0 = 10 − 2·4 = 2.
        // check_agree verifies presolve-ON == presolve-OFF AND that the
        // postsolve x is feasible in the ORIGINAL rows — the real test that
        // the general-scale split landed inside both boxes.
        let lp = lp_of(
            &[&[1.0, 2.0]],
            vec![10.0],
            "L",
            vec![-1.0, -2.0],
            vec![0.0, 0.0],
            vec![3.0, 4.0],
        );
        assert!((check_agree(&lp) + 10.0).abs() < 1e-6);
        use crate::simplex::{SimplexConfig, SimplexSolver};
        match run(&lp, 50, false, true, true, false, false) {
            Presolved::Reduced {
                lp: w,
                live,
                postsolve,
            } => {
                assert!(!live.contains(&1), "duplicate must leave live");
                let d = dupcol_recs(&postsolve);
                assert_eq!(d.len(), 1);
                assert_eq!(d[0].scale, 2.0);
                // Solve the reduced LP (m → 10) and replay the split.
                let sol = SimplexSolver::new(SimplexConfig {
                    presolve: false,
                    ..SimplexConfig::default()
                })
                .solve(&w)
                .unwrap();
                let (x, obj) = restore(&postsolve, &live, &sol.x, sol.objective);
                assert_eq!(x, vec![2.0, 4.0], "upper-clamped split");
                assert!((obj + 10.0).abs() < 1e-12);
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    // ---- Aggregator (HiGHS HPresolve::aggregator :6688): implied-free
    // column substitution through an equality row ----

    fn agg_recs(p: &Postsolve) -> Vec<&AggRec> {
        p.recs
            .iter()
            .filter_map(|r| match r {
                PostRec::Agg(a) => Some(a),
                _ => None,
            })
            .collect()
    }

    fn agg_postsolve(lp: &Lp) -> Postsolve {
        match run(lp, 50, false, false, true, false, false) {
            Presolved::Reduced { postsolve, .. } => postsolve,
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    #[test]
    fn agg_free_col_positive_pivot_fill_in_and_objective() {
        // x is a free column of the 3-term equality r0 (x + y + z = 6) and
        // also sits in r1 (x + w >= 3): substituting x = 6 - y - z puts
        // fill-in (y, z) into r1 and moves x's cost onto y, z (obj_const 6).
        // Postsolve branch: x recovered from the row identity.
        let lp = lp_of(
            &[
                &[1.0, 1.0, 1.0, 0.0],
                &[1.0, 0.0, 0.0, 1.0],
                &[0.0, 1.0, -1.0, 1.0],
            ],
            vec![6.0, 3.0, 5.0],
            "EGL",
            vec![1.0, 2.0, 3.0, -1.0],
            vec![f64::NEG_INFINITY, 0.0, 0.0, 0.0],
            vec![f64::INFINITY, 4.0, 4.0, 10.0],
        );
        check_agree(&lp);
        let ps = agg_postsolve(&lp);
        assert!(!agg_recs(&ps).is_empty(), "free column must be aggregated");
    }

    #[test]
    fn agg_negative_pivot_bounded_col_singleton_and_obj_const() {
        // -2x + y + z = -4 with x >= 0 only: the row implies x in [2, 5],
        // so the bound x >= 0 is redundant (implied free). x is a column
        // SINGLETON, so the row and the column both vanish; x's cost 5
        // becomes obj_const 5*(-4)/(-2) = 10 plus 2.5*(y+z).
        let lp = lp_of(
            &[&[-2.0, 1.0, 1.0, 0.0], &[0.0, 1.0, 1.0, -1.0]],
            vec![-4.0, 0.5],
            "EG",
            vec![5.0, 1.0, 2.0, 1.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![f64::INFINITY, 3.0, 3.0, 2.0],
        );
        let obj = check_agree(&lp);
        let ps = agg_postsolve(&lp);
        let recs = agg_recs(&ps);
        assert_eq!(recs.len(), 1, "exactly the singleton x is aggregated");
        assert_eq!(recs[0].orig, 0);
        assert!((recs[0].a + 2.0).abs() < 1e-12 && (recs[0].b + 4.0).abs() < 1e-12);
        // r1 forces y + z >= 0.5 + w; optimum w = 0, y = 0.5, z = 0:
        // x = 2.25, obj = 5*2.25 + 0.5 = 11.75.
        assert!((obj - 11.75).abs() < 1e-6, "obj {obj}");
    }

    #[test]
    fn agg_implied_free_through_another_row() {
        // x >= 0 is NOT implied by r0 alone (x = 10 - y - z reaches -6), but
        // r1 (x - w >= 0, w >= 0) implies x >= 0, so x is implied free and r1
        // survives (rewritten) as the guard: -y - z - w >= -10.
        let lp = lp_of(
            &[&[1.0, 1.0, 1.0, 0.0], &[1.0, 0.0, 0.0, -1.0]],
            vec![10.0, 0.0],
            "EG",
            vec![1.0, -1.0, -2.0, -0.5],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![f64::INFINITY, 8.0, 8.0, 5.0],
        );
        check_agree(&lp);
        let ps = agg_postsolve(&lp);
        assert!(
            agg_recs(&ps).iter().any(|r| r.orig == 0),
            "x implied free through r1 must be aggregated"
        );
    }

    #[test]
    fn agg_not_implied_free_is_not_substituted() {
        // Same as above WITHOUT r1: x >= 0 is a real constraint (the row
        // alone allows x down to -6), so substituting x would silently drop
        // it. No record for x may exist; the LP must still agree.
        let lp = lp_of(
            &[&[1.0, 1.0, 1.0]],
            vec![10.0],
            "E",
            vec![1.0, -1.0, -2.0],
            vec![0.0, 0.0, 0.0],
            vec![f64::INFINITY, 8.0, 8.0],
        );
        let obj = check_agree(&lp);
        // x = 10 - y - z >= 0 binds (y = z = 8 would need x = -6): minimise
        // x - y - 2z = 10 - 2y - 3z at z = 8, y = 2 (x = 0): obj -18.
        assert!((obj + 18.0).abs() < 1e-6, "obj {obj}");
        let ps = agg_postsolve(&lp);
        assert!(
            agg_recs(&ps).iter().all(|r| r.orig != 0),
            "x >= 0 is not implied: x must stay"
        );
    }

    #[test]
    fn agg_fill_in_cap_refuses_dense_substitution() {
        // x (implied free by r0: x + Σy = 20) also sits in three inequality
        // rows with disjoint columns z_i: substituting would add 11 entries
        // to each of them (net fill-in 18 > the cap of 10). Nothing may fire.
        let mut rows: Vec<Vec<f64>> = Vec::new();
        let n = 1 + 11 + 3;
        let mut r0 = vec![0.0; n];
        for v in r0.iter_mut().take(12) {
            *v = 1.0;
        }
        rows.push(r0);
        for i in 0..3 {
            let mut r = vec![0.0; n];
            r[0] = 1.0;
            r[12 + i] = 1.0;
            rows.push(r);
        }
        let rr: Vec<&[f64]> = rows.iter().map(|r| r.as_slice()).collect();
        let mut obj = vec![0.0; n];
        obj[0] = 0.5;
        for v in obj.iter_mut().take(12).skip(1) {
            *v = 1.0;
        }
        for v in obj.iter_mut().skip(12) {
            *v = -1.0;
        }
        let mut lo = vec![0.0; n];
        let mut up = vec![1.0; n];
        up[0] = 100.0;
        for u in up.iter_mut().skip(12) {
            *u = 1.5;
        }
        lo[0] = 0.0;
        let lp = lp_of(&rr, vec![20.0, 21.0, 21.0, 21.0], "ELLL", obj, lo, up);
        check_agree(&lp);
        let ps = agg_postsolve(&lp);
        assert!(agg_recs(&ps).is_empty(), "fill-in 18 > 10 must refuse");
    }

    #[test]
    fn agg_cascade_restores_in_reverse_order() {
        // x and y are both free; r0: x + y + z = 6, r1: y + u + v = 5,
        // r2: x + u >= 1. y can be eliminated AFTER x although x's record
        // lists y among its rest columns; reverse replay must resolve y
        // before x.
        let lp = lp_of(
            &[
                &[1.0, 1.0, 1.0, 0.0, 0.0],
                &[0.0, 1.0, 0.0, 1.0, 1.0],
                &[1.0, 0.0, 0.0, 1.0, 0.0],
            ],
            vec![6.0, 5.0, 1.0],
            "EEG",
            vec![1.0, 0.0, 2.0, 1.0, 3.0],
            vec![f64::NEG_INFINITY, f64::NEG_INFINITY, 0.0, 0.0, 0.0],
            vec![f64::INFINITY, f64::INFINITY, 2.0, 3.0, 3.0],
        );
        check_agree(&lp);
        let ps = agg_postsolve(&lp);
        assert!(
            agg_recs(&ps).len() >= 2,
            "both free columns must be aggregated, got {}",
            agg_recs(&ps).len()
        );
    }

    #[test]
    fn agg_random_equality_lps_agree_with_presolve_off() {
        // Differential fuzz: sparse equality-heavy LPs with free / one-sided /
        // boxed columns; presolve ON (aggregator included) must reproduce the
        // presolve-OFF status/objective and the postsolved x must satisfy the
        // ORIGINAL rows and bounds (check_agree asserts all three).
        use crate::rng::XorShift64;
        let mut fired = 0usize;
        let seeds: u64 = std::env::var("RUSOLVE_AGG_FUZZ_SEEDS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(400);
        let from: u64 = std::env::var("RUSOLVE_AGG_FUZZ_FROM")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(1);
        for seed in from..=seeds {
            let mut rng = XorShift64::new(seed);
            let n = 4 + (rng.next_u64() % 5) as usize;
            let m = 2 + (rng.next_u64() % 4) as usize;
            let mut rows: Vec<Vec<f64>> = vec![vec![0.0; n]; m];
            for r in rows.iter_mut() {
                for v in r.iter_mut() {
                    if rng.next_f64() < 0.5 {
                        *v = (1 + (rng.next_u64() % 3) as i32) as f64
                            * if rng.next_f64() < 0.4 { -1.0 } else { 1.0 };
                    }
                }
            }
            let (mut lo, mut up, mut xs) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
            for j in 0..n {
                match rng.next_u64() % 4 {
                    0 => {
                        lo[j] = f64::NEG_INFINITY;
                        up[j] = f64::INFINITY;
                    }
                    1 => {
                        lo[j] = 0.0;
                        up[j] = f64::INFINITY;
                    }
                    2 => {
                        lo[j] = -((rng.next_u64() % 4) as f64);
                        up[j] = (rng.next_u64() % 5) as f64;
                    }
                    _ => {
                        lo[j] = f64::NEG_INFINITY;
                        up[j] = (rng.next_u64() % 5) as f64;
                    }
                }
                xs[j] = (rng.next_u64() % 3) as f64;
                xs[j] = xs[j].clamp(lo[j], up[j]);
            }
            let mut senses = String::new();
            let mut rhs = Vec::new();
            for r in &rows {
                let act: f64 = r.iter().zip(&xs).map(|(a, x)| a * x).sum();
                match rng.next_u64() % 4 {
                    0 => {
                        senses.push('L');
                        rhs.push(act + (rng.next_u64() % 3) as f64);
                    }
                    1 => {
                        senses.push('G');
                        rhs.push(act - (rng.next_u64() % 3) as f64);
                    }
                    _ => {
                        senses.push('E');
                        rhs.push(act);
                    }
                }
            }
            let obj: Vec<f64> = (0..n).map(|_| rng.uniform(-3.0, 3.0).round()).collect();
            let rr: Vec<&[f64]> = rows.iter().map(|r| r.as_slice()).collect();
            let lp = lp_of(&rr, rhs, &senses, obj, lo, up);
            check_agree(&lp);
            if let Presolved::Reduced { postsolve, .. } =
                run(&lp, 50, false, false, true, false, false)
                && !agg_recs(&postsolve).is_empty()
            {
                fired += 1;
            }
        }
        assert!(fired >= 20, "fuzz must exercise the aggregator ({fired})");
    }

    // ---- parallel-column dominance fixing (HiGHS :7399–7633 LP arms) ----

    /// min x + 5y s.t. x + 2y ≥ 6, x + 2y + s ≤ 14 (s ∈ [5,10]).
    /// x and y are proportional (scale 2) on their shared support; r2
    /// implies x ≤ 9 < 10 strictly, so x can absorb y's compensation ⇒
    /// objDiff = 1·2 − 5 = −3 < −tol fixes y at its lower bound 0.
    #[test]
    fn parcol_dom_fixes_expensive_duplicate_at_lower() {
        let lp = lp_of(
            &[&[1.0, 2.0, 0.0], &[1.0, 2.0, 1.0]],
            vec![6.0, 14.0],
            "GL",
            vec![1.0, 5.0, 0.0],
            vec![0.0, 0.0, 5.0],
            vec![10.0, 20.0, 10.0],
        );
        assert!((check_agree(&lp) - 6.0).abs() < 1e-6);
        match run(&lp, 50, false, false, true, true, false) {
            Presolved::Reduced { postsolve, .. } => {
                assert_eq!(
                    postsolve.fixed[1],
                    Some(0.0),
                    "expensive parallel column must be fixed at its lower bound"
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    /// Mirror arm: objDiff > tol and the surviving column's LOWER bound is
    /// strictly implied ⇒ the duplicate is fixed at its UPPER bound.
    /// min 3x + y s.t. x + 2y ≥ 6, x + 2y − s ≥ 2, s ∈ [4,8], y ∈ [0,2]:
    /// scale 2, objDiff = 3·2 − 1 = 5 > 0; r1/r2 imply x ≥ 2 > 0 strictly,
    /// so y := 2 (which is also the true optimum: y is the cheap side).
    #[test]
    fn parcol_dom_mirror_arm_fixes_duplicate_at_upper() {
        let lp = lp_of(
            &[&[1.0, 2.0, 0.0], &[1.0, 2.0, -1.0]],
            vec![6.0, 2.0],
            "GG",
            vec![3.0, 1.0, 0.0],
            vec![0.0, 0.0, 4.0],
            vec![10.0, 2.0, 8.0],
        );
        assert!((check_agree(&lp) - 8.0).abs() < 1e-6);
        match run(&lp, 50, false, false, true, true, false) {
            Presolved::Reduced { postsolve, .. } => {
                assert_eq!(
                    postsolve.fixed[1],
                    Some(2.0),
                    "cheap parallel duplicate must be fixed at its upper bound"
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    /// Compensation blocked ⇒ no fix: same shape as the first test but r2
    /// does not cap x below its explicit bound (40 − 2·0 − 4 = 36 > 10),
    /// so x cannot be certified to absorb y's shift and nothing fires
    /// (the weak gate requires the implied cap ≤ the current bound; after
    /// the tightening block clamps bounds to caps, only a LOOSE cap — or
    /// none — blocks the arm).
    #[test]
    fn parcol_dom_blocked_without_implied_bound() {
        let lp = lp_of(
            &[&[1.0, 2.0, 0.0], &[1.0, 2.0, 1.0]],
            vec![6.0, 40.0],
            "GL",
            vec![1.0, 5.0, 0.0],
            vec![0.0, 0.0, 4.0],
            vec![10.0, 20.0, 10.0],
        );
        assert!((check_agree(&lp) - 6.0).abs() < 1e-6);
        match run(&lp, 50, false, false, true, true, false) {
            Presolved::Reduced { postsolve, .. } => {
                assert_ne!(
                    postsolve.fixed[1],
                    Some(0.0),
                    "dominance must NOT fire when no row caps the compensating column"
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    /// Fuzz seed 145 regression: after the pure-bound fold of
    /// −5·x1 ≥ −15 (x1 ≤ 3), x1 and x2 become an ANTI-PARALLEL cost-equal
    /// pair (scale −1, objDiff 0) on r1. The literal HiGHS case-table
    /// transcription fixed x1 := 0 (0 − 3·0 = 0 in the scaled objective)
    /// and lost the optimum −9; the direction-aware table must leave x1
    /// alone (x2's lower bound is not implied, so nothing absorbs).
    #[test]
    fn parcol_dom_antiparallel_cost_equal_keeps_optimum() {
        let lp = lp_of(
            &[&[0.0, -5.0, 0.0, 0.0, 0.0], &[0.0, 5.0, -5.0, -4.0, 0.0]],
            vec![-15.0, -5.0],
            "GG",
            vec![0.0, -3.0, 3.0, 4.0, 4.0],
            vec![0.0; 5],
            vec![10.0, 10.0, f64::INFINITY, f64::INFINITY, f64::INFINITY],
        );
        assert!((check_agree(&lp) - (-9.0)).abs() < 1e-6);
        match run(&lp, 50, false, false, true, true, false) {
            Presolved::Reduced { postsolve, .. } => {
                assert_ne!(postsolve.fixed[1], Some(0.0), "x1 must not be fixed at 0");
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    /// Differential fuzz for the dominance arms: random LPs whose matrices
    /// are GUARANTEED to contain proportional columns (each seed duplicates
    /// a column with a random positive scale, random cost and bounds).
    /// Presolve ON must agree with presolve OFF on every seed.
    #[test]
    fn parcol_dom_random_parallel_lps_agree_with_presolve_off() {
        use crate::rng::XorShift64;
        let seeds: u64 = std::env::var("RUSOLVE_PARCOLDOM_FUZZ_SEEDS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(400);
        let mut fired = 0usize;
        for seed in 0..seeds {
            let mut rng = XorShift64::new(seed);
            let n = 3 + (rng.next_u64() % 4) as usize;
            let m = 2 + (rng.next_u64() % 4) as usize;
            // base matrix
            let mut cols: Vec<Vec<f64>> = (0..n)
                .map(|_| {
                    (0..m)
                        .map(|_| {
                            if rng.next_f64() < 0.5 {
                                0.0
                            } else {
                                ((rng.next_u64() % 5) as f64 + 1.0)
                                    * if rng.next_f64() < 0.5 { -1.0 } else { 1.0 }
                            }
                        })
                        .collect()
                })
                .collect();
            // duplicate column 0 with a random positive scale (and random
            // cost/bounds): guarantees a proportional pair per seed
            let scale = (rng.next_u64() % 3) as f64 + 1.0;
            cols.push(cols[0].iter().map(|v| v * scale).collect());
            let n = n + 1;
            let mut lo = vec![0.0; n];
            let mut up = vec![10.0; n];
            for j in 0..n {
                match rng.next_u64() % 3 {
                    0 => up[j] = f64::INFINITY,
                    1 => lo[j] = -((rng.next_u64() % 3) as f64),
                    _ => {}
                }
            }
            let obj: Vec<f64> = (0..n).map(|_| (rng.next_u64() % 9) as f64 - 4.0).collect();
            // a feasible point to calibrate rhs
            let xs: Vec<f64> = (0..n).map(|_| (rng.next_u64() % 5) as f64).collect();
            let mut senses = String::new();
            let mut rhs = Vec::new();
            for i in 0..m {
                let act: f64 = cols.iter().zip(&xs).map(|(c, x)| c[i] * x).sum();
                match rng.next_u64() % 3 {
                    0 => {
                        senses.push('L');
                        rhs.push(act + (rng.next_u64() % 3) as f64);
                    }
                    1 => {
                        senses.push('G');
                        rhs.push(act - (rng.next_u64() % 3) as f64);
                    }
                    _ => {
                        senses.push('E');
                        rhs.push(act);
                    }
                }
            }
            let rows: Vec<Vec<f64>> = (0..m)
                .map(|i| cols.iter().map(|c| c[i]).collect())
                .collect();
            let rr: Vec<&[f64]> = rows.iter().map(|r| r.as_slice()).collect();
            if seed
                == std::env::var("RUSOLVE_PARCOLDOM_FUZZ_DUMP")
                    .ok()
                    .and_then(|s| s.parse::<u64>().ok())
                    .unwrap_or(u64::MAX)
            {
                eprintln!(
                    "DUMP rows={rr:?} rhs={rhs:?} senses={senses} obj={obj:?} lo={lo:?} up={up:?}"
                );
            }
            let lp = lp_of(&rr, rhs, &senses, obj, lo, up);
            check_agree(&lp);
            if let Presolved::Reduced { postsolve, .. } =
                run(&lp, 50, false, false, true, true, false)
            {
                // count columns whose original box was OPEN and got fixed:
                // the merge path is OFF, so these are dominance fixes
                fired += (0..lp.num_col)
                    .filter(|&j| postsolve.fixed[j].is_some() && lp.lower[j] < lp.upper[j])
                    .count();
            }
        }
        assert!(
            fired >= 20,
            "fuzz must exercise the dominance arms ({fired})"
        );
    }

    // ---- singleton-cost VUB implied-equality chain (HiGHS :263–281) ----

    /// mc11 shape: min −x + 2b s.t. x − 212·b ≤ 0, x ∈ [0,1], b ∈ [0,1].
    /// b is a column singleton with min-space cost +2 in an 'L' row with
    /// coefficient −212 < 0 ⇒ every optimum has y₀ ≤ −2/212 < 0 ⇒ the row
    /// is tight at every optimum ⇒ converted to an equation and eaten by
    /// the DTE chain in the same pass.
    #[test]
    fn singeq_mc11_shape_converts_and_chain_eats() {
        let lp = lp_of(
            &[&[1.0, -212.0]],
            vec![0.0],
            "L",
            vec![-1.0, 2.0],
            vec![0.0, 0.0],
            vec![1.0, 1.0],
        );
        // optimum: b = 1/212, x = 1 ⇒ obj = −1 + 2/212
        assert!((check_agree(&lp) - (-1.0 + 2.0 / 212.0)).abs() < 1e-9);
        match run(&lp, 50, false, false, true, false, true) {
            Presolved::Reduced { lp: red, .. } => {
                assert_eq!(red.num_row, 0, "converted row must be eaten by the chain");
                assert!(
                    red.num_col <= 1,
                    "at most one of the pair survives, got {}",
                    red.num_col
                );
            }
            // everything fixed is even stronger — presolve solved it outright
            Presolved::Terminal { .. } => {}
        }
        // rule OFF: nothing else can convert or eat this row
        match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: red, .. } => {
                assert_eq!(red.num_row, 1, "row must survive with the rule off");
                assert_eq!(red.num_col, 2);
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    /// Blocked arm: b's improving bound is NOT implied (the other column
    /// is free below) ⇒ no conversion even though every other gate holds.
    /// min −x + 2b s.t. x − 212·b ≤ 0, x ∈ (−∞,1], b ∈ [0,1].
    #[test]
    fn singeq_blocked_when_bound_not_implied() {
        let lp = lp_of(
            &[&[1.0, -212.0]],
            vec![0.0],
            "L",
            vec![-1.0, 2.0],
            vec![f64::NEG_INFINITY, 0.0],
            vec![1.0, 1.0],
        );
        // optimum: x = 1 at b = 1/212 ⇒ same value as the mc11 shape, but
        // reached WITHOUT the conversion (x's own ub supplies the bound)
        assert!((check_agree(&lp) - (-1.0 + 2.0 / 212.0)).abs() < 1e-9);
        match run(&lp, 50, false, false, true, false, true) {
            Presolved::Reduced { lp: red, .. } => {
                assert_eq!(
                    red.num_row, 1,
                    "un-implied improving bound must block the conversion"
                );
            }
            other => panic!("expected Reduced, got {other:?}"),
        }
    }

    /// Unbounded-LP caveat: the dual argument is vacuous without an
    /// optimum, and the conversion must not corrupt the status. Here the
    /// rule fires on (b, row 0) while the (z, w) ray drives the LP
    /// unbounded — the full solve with the rule ON (default config) must
    /// still report Unbounded, agreeing with presolve OFF.
    #[test]
    fn singeq_unbounded_lp_keeps_status() {
        use crate::simplex::{SimplexConfig, SimplexSolver};
        let lp = lp_of(
            &[&[1.0, -212.0, 0.0, 0.0], &[0.0, 0.0, 1.0, -1.0]],
            vec![0.0, 0.0],
            "LG",
            vec![0.0, 1.0, 1.0, -3.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![1.0, 1.0, f64::INFINITY, f64::INFINITY],
        );
        // asserts status agreement (Unbounded == Unbounded) + would panic
        // on any objective/feasibility drift
        let _ = check_agree(&lp);
        let solved = SimplexSolver::new(SimplexConfig {
            presolve: true,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_eq!(
            solved.status,
            crate::simplex::SimplexStatus::Unbounded,
            "rule firing must not flip an unbounded LP's status"
        );
        // and the rule really fired on row 0 in this instance
        let rows_on = match run(&lp, 50, false, false, true, false, true) {
            Presolved::Reduced { lp: red, .. } => red.num_row,
            // presolve resolving the (correct) unboundedness outright is
            // the strongest form of "status preserved"
            Presolved::Terminal { status } => {
                assert_eq!(
                    status,
                    crate::simplex::SimplexStatus::Unbounded,
                    "terminal must be the correct status"
                );
                0
            }
        };
        let rows_off = match run(&lp, 50, false, false, true, false, false) {
            Presolved::Reduced { lp: red, .. } => red.num_row,
            Presolved::Terminal { .. } => panic!("expected Reduced"),
        };
        assert!(
            rows_on < rows_off,
            "rule must fire ({rows_on} vs {rows_off})"
        );
    }

    /// Randomized differential: every seed builds a random LP plus one
    /// appended singleton-cost row (fresh column b, nonzero cost, 1–2
    /// random co-columns, rhs calibrated at a feasible point); presolve
    /// with the rule ON must agree with presolve OFF on status, objective
    /// and original-feasibility of the reconstructed point, and the rule
    /// must actually fire on a healthy fraction of seeds.
    #[test]
    fn singeq_random_lps_agree_with_presolve_off() {
        use crate::rng::XorShift64;
        let seeds: u64 = std::env::var("RUSOLVE_VUBEQ_FUZZ_SEEDS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(400);
        let mut fired = 0usize;
        for seed in 0..seeds {
            let mut rng = XorShift64::new(1_000_000 + seed);
            let n = 3 + (rng.next_u64() % 4) as usize;
            let m = 2 + (rng.next_u64() % 3) as usize;
            // random base matrix + boxes + costs (parcol fuzz recipe)
            let mut cols: Vec<Vec<f64>> = (0..n)
                .map(|_| {
                    (0..m)
                        .map(|_| {
                            if rng.next_f64() < 0.5 {
                                0.0
                            } else {
                                ((rng.next_u64() % 5) as f64 + 1.0)
                                    * if rng.next_u64().is_multiple_of(2) {
                                        -1.0
                                    } else {
                                        1.0
                                    }
                            }
                        })
                        .collect()
                })
                .collect();
            let mut lo = vec![0.0f64; n];
            let mut up = vec![10.0f64; n];
            for j in 0..n {
                match rng.next_u64() % 3 {
                    0 => up[j] = f64::INFINITY,
                    1 => lo[j] = -((rng.next_u64() % 3) as f64),
                    _ => {}
                }
            }
            let mut obj: Vec<f64> = (0..n).map(|_| (rng.next_u64() % 9) as f64 - 4.0).collect();
            // feasible calibration point
            let mut xs: Vec<f64> = (0..n).map(|_| (rng.next_u64() % 5) as f64).collect();
            // append the singleton-cost row: fresh column b (index n) with
            // nonzero cost, coefficient ±(1..3)·{1, 10, 212} for margin
            // variety, plus 1–2 random existing columns.
            let ab = ((rng.next_u64() % 3) as f64 + 1.0)
                * [1.0, 10.0, 212.0][(rng.next_u64() % 3) as usize]
                * if rng.next_u64().is_multiple_of(2) {
                    -1.0
                } else {
                    1.0
                };
            let mut row = vec![0.0; n];
            for _ in 0..1 + (rng.next_u64() % 2) {
                let j = (rng.next_u64() % n as u64) as usize;
                row[j] = ((rng.next_u64() % 4) as f64 + 1.0)
                    * if rng.next_u64().is_multiple_of(2) {
                        -1.0
                    } else {
                        1.0
                    };
            }
            let sense = if rng.next_u64().is_multiple_of(2) {
                'L'
            } else {
                'G'
            };
            let xb = (rng.next_u64() % 3) as f64; // b's calibration value
            // rhs calibrated at (xs, xb) with slack so the row is not
            // forcing on its own
            let act: f64 = row.iter().zip(&xs).map(|(a, x)| a * x).sum::<f64>() + ab * xb;
            let rhs = match sense {
                'L' => act + (rng.next_u64() % 3) as f64,
                _ => act - (rng.next_u64() % 3) as f64,
            };
            for (i, a) in row.iter().enumerate() {
                // push unconditionally (0.0 included) to keep column
                // lengths aligned; explicit zeros are dropped by Sp.
                cols[i].push(*a);
            }
            cols.push(vec![0.0; m]);
            cols[n].push(ab); // b's single entry lives in the appended row
            let bu = [1.0, 5.0, 10.0][(rng.next_u64() % 3) as usize];
            lo.push(0.0);
            up.push(bu);
            let cb = (rng.next_u64() % 4) as f64 + 1.0;
            obj.push(if rng.next_u64().is_multiple_of(2) {
                cb
            } else {
                -cb
            });
            xs.push(xb);
            let m = m + 1;

            let mut senses = String::new();
            let mut rhs_vec = Vec::new();
            for i in 0..m - 1 {
                let act: f64 = cols.iter().zip(&xs).map(|(c, x)| c[i] * x).sum();
                match rng.next_u64() % 3 {
                    0 => {
                        senses.push('L');
                        rhs_vec.push(act + (rng.next_u64() % 3) as f64);
                    }
                    1 => {
                        senses.push('G');
                        rhs_vec.push(act - (rng.next_u64() % 3) as f64);
                    }
                    _ => {
                        senses.push('E');
                        rhs_vec.push(act);
                    }
                }
            }
            senses.push(sense);
            rhs_vec.push(rhs);
            let rows: Vec<Vec<f64>> = (0..m)
                .map(|i| cols.iter().map(|c| c[i]).collect())
                .collect();
            let rr: Vec<&[f64]> = rows.iter().map(|r| r.as_slice()).collect();
            if seed
                == std::env::var("RUSOLVE_VUBEQ_FUZZ_DUMP")
                    .ok()
                    .and_then(|s| s.parse::<u64>().ok())
                    .unwrap_or(u64::MAX)
            {
                eprintln!(
                    "DUMP rows={rr:?} rhs={rhs_vec:?} senses={senses} obj={obj:?} lo={lo:?} up={up:?}"
                );
            }
            let lp = lp_of(&rr, rhs_vec, &senses, obj, lo, up);
            if std::env::var("RUSOLVE_VUBEQ_FUZZ_VERBOSE").is_ok() {
                eprintln!("seed {seed}");
            }
            check_agree(&lp);
            let rows_on = match run(&lp, 50, false, false, true, false, true) {
                Presolved::Reduced { lp: red, .. } => red.num_row,
                Presolved::Terminal { .. } => 0,
            };
            let rows_off = match run(&lp, 50, false, false, true, false, false) {
                Presolved::Reduced { lp: red, .. } => red.num_row,
                Presolved::Terminal { .. } => 0,
            };
            if rows_on < rows_off {
                fired += 1;
            }
        }
        assert!(
            fired >= 20,
            "fuzz must exercise the singleton-eq conversions ({fired})"
        );
    }
}
