//! Product-form factorization (PFI) of the simplex basis.
//!
//! The basis `B` is represented as the initial factorization times a file of
//! eta matrices, `B = B0 · E1 · E2 · … · Em`, where each `Ei` is the identity
//! except for one column. Solves therefore touch only the dense LU of `B0`
//! plus a sparse sweep per eta; updates are O(nnz of one basis column) plus a
//! solve, with no refactorization until the file ages past a threshold.
//!
//! # The update identity (re-derived, not copied)
//!
//! Let `a_q = Σ_i b̄_i B_i` be the entering column in current basis
//! coordinates, replacing basis position `p`. Then `B' = B·E` where `E` is
//! the identity except column `p`, which holds `b̄` (with `E[p][p] = b̄_p`);
//! `det E = b̄_p`, so the update is valid exactly when `b̄_p` stays away from
//! zero — the pivot guard in [`ProductFormFactorization::update`]. From
//! `E = I + (b̄ − e_p) e_pᵀ`, Sherman–Morrison gives
//! `E⁻¹ = I − (b̄ − e_p) e_pᵀ / b̄_p`, i.e. for a vector `y`:
//!
//! ```text
//! t = y[p] / b̄_p        x[i] = y[i] − b̄_i · t  (i ≠ p)      x[p] = t
//! ```
//!
//! and transposing, `E⁻ᵀ y` changes only component `p`:
//!
//! ```text
//! x[p] = (y[p] − Σ_{i≠p} b̄_i·y[i]) / b̄_p        x[j] = y[j]  (j ≠ p)
//! ```
//!
//! Since `B⁻¹ b = E_m⁻¹ … E_1⁻¹ B0⁻¹ b` (rightmost operator acts first),
//! [`solve_b`] does the dense-LU solve with `B0`, then walks the eta file
//! **oldest-first**. Transposing reverses the product:
//! `B⁻ᵀ b = B0⁻ᵀ E1⁻ᵀ … E_m⁻¹ b`, so [`solve_bt`] walks the eta file
//! **newest-first** and applies the `B0ᵀ` solve *last*. This matches the
//! procedure published as pseudocode in
//! Huangfu & Hall, arXiv:1503.01889 (§product form of inverse) and the BASP
//! eta-file description of Forrest & Gold (ACM TOMS 18(3), 1992).
//!
//! [`solve_b`]: ProductFormFactorization::solve_b
//! [`solve_bt`]: ProductFormFactorization::solve_bt

use crate::dense::{DenseLu, DenseMatrix};
use crate::hvec::HVec;
use crate::lu::{FactorWs, SparseLu};
use crate::{FactorError, PfiConfig, PfiEngine};
use std::cell::RefCell;

/// Outcome of a successful [`ProductFormFactorization::update`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UpdateStatus {
    /// The eta was appended and the eta file is still within limits.
    Updated,
    /// The eta was appended and the file reached [`PfiConfig::eta_limit`],
    /// so the factorization was rebuilt from the explicit basis columns and
    /// the eta file is now empty.
    Refactorized,
}

/// One eta matrix: identity except column `pos`, which holds the basis
/// coordinates `val[idx]` of the entering column (the `pos` entry itself is
/// excluded and carried as `pivot`).
#[derive(Debug, Clone)]
struct Eta {
    pos: usize,
    pivot: f64,
    idx: Vec<i32>,
    val: Vec<f64>,
}

/// Product-form basis factorization: dense LU of the initial basis plus an
/// eta file of single-column updates. See the [module docs](self).
#[derive(Debug, Clone)]
pub struct ProductFormFactorization {
    n: usize,
    lu: Engine,
    /// Sticky engine choice: starts from config.engine, flips to
    /// Dense permanently on a sparse Singular (see [`factor_with`]).
    engine_choice: PfiEngine,
    etas: Vec<Eta>,
    /// Number of (re)factorizations performed on this object (build counts
    /// as 0): lets callers tell when the eta file was reset.
    epoch: u64,
    /// Total stored nonzeros over the eta file (+1 per eta for its pivot):
    /// the per-solve work the file adds on top of the LU solve.
    eta_nnz: usize,
    /// Accumulated eta-attributable solve work since the last
    /// refactorization (a `Cell`: solves take `&self`). Compared against the
    /// factorization's own build work by [`ProductFormFactorization::update`]
    /// — refactor once the etas have cost more than a rebuild would
    /// (HiGHS's synthetic-clock reinversion, HEkk.cpp `updateFactor`).
    eta_ticks: std::cell::Cell<f64>,
    refactor_nanos: std::cell::Cell<f64>,
    /// Nonzeros of the *current* basis columns: each column is an ascending
    /// sorted `(row, value)` list with exact zeros dropped. Source of truth
    /// for the basis, kept sparse so `refine`'s residual sweeps cost O(nnz)
    /// per solve instead of the O(m²) of dense iteration; the rare
    /// factorization sites (build + refactor) materialize dense columns via
    /// [`dense_cols`]. Ascending row order is load-bearing: it keeps the
    /// refinement accumulation sequence bit-identical to the old dense loops
    /// (see [`ProductFormFactorization::refine`]).
    basis_nz: Vec<Vec<(usize, f64)>>,
    /// Row-wise mirror of `basis_nz` (entry `(column, value)`, ASCENDING by
    /// column within a row — kept sorted by [`Self::swap_col`] so the
    /// sparse LU's triangular peel can read it in place of a fresh
    /// transpose, see [`SparseLu::factor_sparse_ws`]): also lets the
    /// transpose residual visit only the nonzeros of `x`.
    basis_rows: Vec<Vec<(u32, f64)>>,
    /// `basis_rows` is in sync with `basis_nz`. Both update paths maintain
    /// it incrementally via [`Self::swap_col`] (the dense `update`
    /// used to drop this flag instead, disabling row-wise residuals).
    rows_ok: bool,
    /// Upper bound on |B[i][j]| over the explicit basis columns: the ||B||
    /// factor of the refinement residual floor. Maintained conservatively
    /// (raised on column swaps, never lowered) — an overestimate only makes
    /// refinement slightly more eager, which is safe.
    col_max: f64,
    config: PfiConfig,
    /// Reusable m-length scratch for the refinement stack (see
    /// [`SolveWorkspace`]); owned here so it lives exactly as long as the
    /// factorization and travels with it through `clone`/restore.
    ws: SolveWorkspace,
    /// Persistent factor-build workspace
    /// plus the recycled previous LU (clones start empty; pure scratch).
    fws: FactorWs,
}

/// Persistent per-solve scratch for the refinement stack. Before this, every refined
/// solve allocated 4 fresh `vec![0.0; m]` in `refine` plus one in
/// `solve_bt_raw` — ×5 refined solves per pivot, on top of the raw LU
/// engine's own per-call temporaries (now persistent too, see
/// [`crate::lu::SparseLu::solve`]): ~25–35 Vec allocations per pivot.
/// Buffers are sized once (first use) and reused across refinement rounds,
/// solves, and pivots.
///
/// # Reentrancy
///
/// Each buffer is borrowed for the duration of the solve that uses it, and
/// no solve nests another solve touching the SAME buffer: `refine` holds
/// `r`/`d`/`hi_part`/`lo_part` while its correction round calls the raw
/// solves, which touch only `y` and the LU engine's own disjoint scratch.
///
/// # Bit-exactness
///
/// Every buffer is fully overwritten before its first read each round —
/// the residual sweeps, permute loops and `copy_from_slice` write all m
/// entries — exactly matching the fresh-allocation behavior, so the reuse
/// is invisible to the fp operation sequence (the old `vec![0.0; m]`
/// zero-fills were dead stores).
#[derive(Debug, Clone, Default)]
struct SolveWorkspace {
    /// `refine`: residual `rhs − B̃·x` (forward) or `rhs − B̃ᵀ·x` (transpose).
    r: RefCell<Vec<f64>>,
    /// `refine`: correction solve `B̃⁻¹·r` / `B̃⁻ᵀ·r`.
    d: RefCell<Vec<f64>>,
    /// `refine`: double-word row accumulators of the forward residual.
    hi_part: RefCell<Vec<f64>>,
    lo_part: RefCell<Vec<f64>>,
    /// `solve_bt_raw`: eta-sweep working vector.
    y: RefCell<Vec<f64>>,
    /// Ascending nonzero positions of `y` during a sparse-rhs BTRAN
    /// eta pass (see [`ProductFormFactorization::solve_bt_raw`]). Dead
    /// entries (position zeroed by a later eta) are filtered at use, never
    /// removed — bounded by nnz(rhs) + eta count.
    bt_pat: RefCell<Vec<u32>>,
    /// Index-carrying twins of `r`,
    /// `d`, `y` for the [`ProductFormFactorization::solve_b_hv`] /
    /// [`ProductFormFactorization::solve_bt_hv`] stack.
    r_hv: RefCell<HVec>,
    d_hv: RefCell<HVec>,
    y_hv: RefCell<HVec>,
    /// Double-word residual accumulators of the sparse refine. Unlike
    /// `hi_part`/`lo_part` these are kept ALL-ZERO between uses (each sparse
    /// residual resets exactly the entries it touched), so a residual costs
    /// O(touched) instead of O(m) fills.
    sp_hi: RefCell<Vec<f64>>,
    sp_lo: RefCell<Vec<f64>>,
    /// Ascending nonzero positions of `x` in the dense-vector
    /// residuals ([`ProductFormFactorization::residual_bwd`] /
    /// [`ProductFormFactorization::residual_fwd`]), gathered by ONE scan.
    res_idx: RefCell<Vec<u32>>,
    /// Rows first touched by the sparse forward residual (in
    /// first-touch order; finalization is order free — one write per row).
    res_touch: RefCell<Vec<u32>>,
    /// Epoch marks (`res_mark[i] == res_stamp` ⇔ touched this call);
    /// bumping the stamp clears all marks in O(1).
    res_mark: RefCell<Vec<u32>>,
    res_stamp: std::cell::Cell<u32>,
}

impl SolveWorkspace {
    /// Borrow one scratch buffer resized to `n` (grown on first use). The
    /// contents are stale garbage on reuse — sound because every consumer
    /// writes all `n` entries before reading any (see bit-exactness note).
    fn buf<'a>(&self, cell: &'a RefCell<Vec<f64>>, n: usize) -> std::cell::RefMut<'a, Vec<f64>> {
        let mut v = cell.borrow_mut();
        if v.len() != n {
            v.resize(n, 0.0);
        }
        v
    }

    /// Start a new mark epoch over `n` slots; returns the mark array
    /// and the fresh stamp (marks are reset only on resize / u32 wrap).
    fn next_stamp(&self, n: usize) -> (std::cell::RefMut<'_, Vec<u32>>, u32) {
        let mut mark = self.res_mark.borrow_mut();
        let mut st = self.res_stamp.get().wrapping_add(1);
        if mark.len() != n || st == 0 {
            mark.clear();
            mark.resize(n, 0);
            st = 1;
        }
        self.res_stamp.set(st);
        (mark, st)
    }
}

/// Gather the ascending nonzero positions of `x` into `idx`; returns
/// `false` (with `idx` truncated garbage) as soon as more than `limit` are
/// found — the density cutoff `count·4 > m` ⇔ `count > ⌊m/4⌋`.
fn collect_nz(x: &[f64], idx: &mut Vec<u32>, limit: usize) -> bool {
    idx.clear();
    for (i, &v) in x.iter().enumerate() {
        if v != 0.0 {
            if idx.len() == limit {
                return false;
            }
            idx.push(i as u32);
        }
    }
    true
}

/// Profiling: per-solve sub-phase counters (printed by the dual engine
/// under RUSOLVE_DUAL_STATS; zero cost when off).
#[derive(Debug, Default)]
pub(crate) struct PfiStats {
    pub eta_ns: std::cell::Cell<f64>,
    pub lu_ns: std::cell::Cell<f64>,
    pub refine_ns: std::cell::Cell<f64>,
    pub n_raw: std::cell::Cell<u64>,
    pub n_corr: std::cell::Cell<u64>,
    /// Cadence census: refine round-1 corrections (fires only when
    /// the round-0 residual cleared the fp floor — the drift signal that
    /// plain-precision eta debt is biting).
    pub n_r1: std::cell::Cell<u64>,
    /// Refactorizations by trigger (indexed by [`RfReason`] as u8).
    pub rf_reason: [std::cell::Cell<u64>; RfReason::COUNT],
    /// Eta-file depth (etas, eta_nnz) and in-epoch raw solves at refactor
    /// time (sums over all refactorizations).
    pub rf_depth_sum: std::cell::Cell<u64>,
    pub rf_depth_max: std::cell::Cell<u64>,
    pub rf_nnz_sum: std::cell::Cell<u64>,
    pub rf_solves_sum: std::cell::Cell<u64>,
    /// Raw solves and round-1 corrections in the epoch that just ended
    /// (deltas tracked against `last_raw`/`last_r1`).
    pub rf_r1_sum: std::cell::Cell<u64>,
    pub last_raw: std::cell::Cell<u64>,
    pub last_r1: std::cell::Cell<u64>,
    /// RUSOLVE_RFDROP=1 probe: max-norm distance from the eta-chain solve
    /// to the fresh-factor solve of the same probe rhs, at refactor time
    /// (sum for a mean, max, sample count).
    pub rf_drift_sum: std::cell::Cell<f64>,
    pub rf_drift_max: std::cell::Cell<f64>,
    pub rf_drift_n: std::cell::Cell<u64>,
    // ---- Sweep-cost census ----------------------------------------------
    // Per-eta-apply sparsity counters, all stats-only: they quantify the
    // hyper-sparsity opportunity the transport census Verdict 3 names
    // (HiGHS's mid-solve sweeps run 2-3x UNDER fresh-sweep cost). Raw
    // solve counts here INCLUDE refine's correction solves; the engine-
    // direct share is `n_raw - n_corr`.
    /// BTRAN raw solves (eta pass newest-first).
    pub bt_raw: std::cell::Cell<u64>,
    /// BTRAN solves where the pattern walk was ATTEMPTED (file gate
    /// open at entry).
    pub bt_gate_on: std::cell::Cell<u64>,
    /// BTRAN solves where the pattern survived the whole pass (never
    /// abandoned at `pat_max`).
    pub bt_kept: std::cell::Cell<u64>,
    /// Etas whose dot ran through the pattern walk vs the plain full-idx
    /// walk.
    pub bt_dot_sparse: std::cell::Cell<u64>,
    pub bt_dot_plain: std::cell::Cell<u64>,
    /// Eta index entries VISITED by BTRAN dots (work done), vs entries
    /// that actually contributed a nonzero term (work that mattered):
    /// eff/touched is the apply-sparsity ratio — on a plain walk it is
    /// exactly the y-density over the eta supports.
    pub bt_touched: std::cell::Cell<u64>,
    pub bt_effective: std::cell::Cell<u64>,
    /// Etas whose dot came out nonzero (the eta actually moved `y`).
    pub bt_hits: std::cell::Cell<u64>,
    /// y nonzeros at the end of the eta pass, sampled from the pattern
    /// when kept; solves whose y left pattern scope counted in `bt_y_dense`.
    pub bt_y_nnz_sum: std::cell::Cell<u64>,
    pub bt_y_dense: std::cell::Cell<u64>,
    /// FTRAN raw solves (eta pass oldest-first): etas that fired (t != 0,
    /// full-idx sweep) vs skipped for free, and total idx entries swept.
    pub ft_raw: std::cell::Cell<u64>,
    pub ft_fired: std::cell::Cell<u64>,
    pub ft_skipped: std::cell::Cell<u64>,
    pub ft_work: std::cell::Cell<u64>,
    /// refine() calls by direction (each pays a residual sweep + round-0
    /// correction solve; `n_corr` already counts the solves).
    pub ref_bt: std::cell::Cell<u64>,
    pub ref_ft: std::cell::Cell<u64>,
    /// refine() calls SKIPPED by the `RUSOLVE_REFINE_AGE` young-file lever.
    pub n_ref_skip: std::cell::Cell<u64>,
    /// Candidate (c′) probe: round-0 residual rnorm-vs-floor stats BEFORE
    /// the unconditional correction solve. `ref0_at_floor` counts solves
    /// whose pre-correction residual already sat at the fp floor — for
    /// those, a residual-gated round 0 would skip the correction solve and
    /// the correction's delta was sub-floor noise anyway.
    pub ref0_total: std::cell::Cell<u64>,
    pub ref0_at_floor: std::cell::Cell<u64>,
    pub ref0_ratio_sum: std::cell::Cell<f64>,
    /// Estimated LIVE eta-file bytes at refactor time (eta_nnz·12 B for
    /// idx+val plus 96 B/eta for the two Vec allocations and struct):
    /// the PFI-accumulation hypothesis for the cont1-class simplex-side
    /// RSS growth (pre13-bucket-close).
    pub rf_bytes_sum: std::cell::Cell<u64>,
    pub rf_bytes_max: std::cell::Cell<u64>,
    /// Wall of the residual evaluations
    /// alone (inside `refine_ns`, which also holds the correction solves),
    /// and residual arms taken: sparse `residual_hv`, dense forward,
    /// column-wise transpose, row-wise transpose.
    pub resid_ns: std::cell::Cell<f64>,
    pub res_arm: [std::cell::Cell<u64>; 4],
}

/// Why a refactorization fired. Distinguishes the cadence triggers
/// (clock/cap, HiGHS `kRebuildReasonSyntheticClockSaysInvert` /
/// `kRebuildReasonUpdateLimitReached`, HEkk.cpp `updateFactor`) from pivot
/// distress and the engine's event rebuilds (verification passes, the
/// Optimal tail — HiGHS rebuilds those too, e.g. `kRebuildReasonPossibly
/// Optimal`, so they must not pollute the cadence number).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RfReason {
    /// Synthetic clock expired (`eta_ticks >= build_work`).
    Clock,
    /// Hard eta-file cap (`eta_limit`).
    Cap,
    /// Pivot below `pivot_tol` on a stale eta file (`RefactorizeNeeded`).
    Pivot,
    /// Engine verification/event rebuild (infeasible/unbounded audit).
    Verify,
    /// The Optimal-exit fresh-factor re-anchor.
    Tail,
    /// Any other direct `refactorize()` call.
    Forced,
}

impl RfReason {
    pub const COUNT: usize = 6;

    pub(crate) fn index(self) -> usize {
        match self {
            RfReason::Clock => 0,
            RfReason::Cap => 1,
            RfReason::Pivot => 2,
            RfReason::Verify => 3,
            RfReason::Tail => 4,
            RfReason::Forced => 5,
        }
    }

    fn label(self) -> &'static str {
        match self {
            RfReason::Clock => "clock",
            RfReason::Cap => "cap",
            RfReason::Pivot => "pivot",
            RfReason::Verify => "verify",
            RfReason::Tail => "tail",
            RfReason::Forced => "forced",
        }
    }
}

thread_local! {
    static PFI_STATS: std::rc::Rc<PfiStats> = std::rc::Rc::default();
}

pub(crate) fn pfi_stats() -> Option<std::rc::Rc<PfiStats>> {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    if !*ON.get_or_init(|| std::env::var("RUSOLVE_DUAL_STATS").is_ok()) {
        return None;
    }
    Some(PFI_STATS.with(|s| s.clone()))
}

/// Cadence lever: `RUSOLVE_RFCAD="K[:M]"` retunes the refactor
/// clock. `K` (>= 1) raises the clock threshold to `K * build_work` —
/// ride staleness harder before paying a rebuild; `M` is a minimum-update
/// floor (HiGHS gates its synthetic clock on
/// `kSyntheticTickReinversionMinUpdateCount = 50`, util/HFactorConst.h,
/// MIT source read as reference). Unset ⇒ `None`: the factor's own
/// [`PfiConfig::clock_scale`] / `min_updates` apply (the dual sets
/// K = 4, everything else keeps `1.0`, bit-identical to the calibrated
/// clock since `x * 1.0 == x` bitwise). Set, it overrides every factor.
pub(crate) fn rfcad() -> Option<(f64, usize)> {
    use std::sync::OnceLock;
    static LEVER: OnceLock<Option<(f64, usize)>> = OnceLock::new();
    *LEVER.get_or_init(|| std::env::var("RUSOLVE_RFCAD").ok().map(|v| parse_rfcad(&v)))
}

/// Parse `K[:M]` (see [`rfcad`]); unparseable halves fall back to the
/// default. Split out from env access so tests need no process env
/// (edition-2024 makes `set_var` unsafe in tests).
fn parse_rfcad(value: &str) -> (f64, usize) {
    let mut it = value.split(':');
    let k = it
        .next()
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(1.0)
        .max(1.0);
    let m = it
        .next()
        .and_then(|s| s.trim().parse::<usize>().ok())
        .unwrap_or(0);
    (k, m)
}

/// Sweep-cost lever: `RUSOLVE_BT_GATE=R` retunes the BTRAN
/// pattern-walk FILE gate — tracking is attempted when
/// `eta_nnz + depth > R·n`. The calibrated `R = 4` asks the file to cost
/// ≥ 4n before paying the O(n) rhs scan, but the pattern-walk cost model itself prices
/// the plain pass at `eta_nnz` work: a file of `R'·n` nnz with a ≤ depth-
/// wide pattern already wastes `(R'−ε)·n` per solve (istanbul eta_nnz 72k
/// vs n = 20.3k sits UNDER the 4n gate and plain-walks 72k entries for a
/// y of ≤ 9 nonzeros). `R = 1` makes the gate compare like-for-like
/// (file work vs scan work). UNSET ⇒ 4.0: the integer comparison
/// `a > 4·n` and the f64 comparison `(a) > 4.0·(n)` agree exactly
/// (integers ≤ 2^53), so the default path is unchanged bit-for-bit.
/// The pattern walk itself is bit-exact against the plain walk (pinned by
/// `bt_eta_pattern_walk_bit_exact`), so this lever moves COST, not the
/// pivot path.
fn bt_gate_ratio() -> f64 {
    use std::sync::OnceLock;
    static LEVER: OnceLock<f64> = OnceLock::new();
    *LEVER.get_or_init(|| {
        std::env::var("RUSOLVE_BT_GATE")
            .map(|v| parse_bt_gate(&v))
            .unwrap_or(4.0)
    })
}

/// Parse helper for [`bt_gate_ratio`] (testable without process env).
fn parse_bt_gate(value: &str) -> f64 {
    value
        .trim()
        .parse::<f64>()
        .ok()
        .filter(|r| (0.0..=1e6).contains(r))
        .unwrap_or(4.0)
}

/// The BTRAN pattern-walk file gate (see [`bt_gate_ratio`]): attempt
/// tracking when the plain eta pass (`eta_nnz + depth` work) costs more
/// than the O(n) rhs scan scaled by `ratio`, and the abandonment
/// threshold `n/16` is meaningful at all. Split out so the default-ratio
/// equivalence to the old integer test `a > 4·n` is pinned by a unit test.
#[inline]
fn bt_gate_open(eta_nnz: usize, depth: usize, n: usize, ratio: f64) -> bool {
    n / 16 > 0 && (eta_nnz + depth) as f64 > ratio * n as f64
}

/// Correction-cost lever: `RUSOLVE_REFINE_AGE=F`
/// skips the refinement stack while the eta file is YOUNG — `eta_ticks <
/// F · build_work`, i.e. younger than the fraction F of the age at which
/// the refactor clock fires. Motivation: we measured the eta-chain drift
/// at refactor points (the OLDEST state the file ever reaches) at
/// ≤ 4.3e-10 with zero round-1 corrections on the whole degenerate trio —
/// the round-0 correction solve (a full raw solve + DW residual per
/// refined solve, ~40% of raw solves) may be insuring against drift the
/// corpus never accumulates. Unset ⇒ `None` ⇒ refine always runs
/// (calibrated behavior). The dense engine has no build-work yardstick ⇒
/// never skips. This lever DOES change solve values (path-visible) —
/// unlike [`bt_gate_ratio`] it must pass the iteration gate, not just the
/// bit-exactness argument.
fn refine_age_frac() -> Option<f64> {
    use std::sync::OnceLock;
    static LEVER: OnceLock<Option<f64>> = OnceLock::new();
    *LEVER.get_or_init(|| {
        std::env::var("RUSOLVE_REFINE_AGE")
            .ok()
            .and_then(|v| parse_refine_age(&v))
    })
}

/// Parse helper for [`refine_age_frac`]: a fraction in (0, 1], else None.
fn parse_refine_age(value: &str) -> Option<f64> {
    let f: f64 = value.trim().parse().ok()?;
    (f > 0.0 && f <= 1.0).then_some(f)
}

/// Candidate (c′) lever: `RUSOLVE_REFINE_GATE=1` makes the round-0
/// correction CONDITIONAL on the double-word residual clearing the fp
/// floor — the very test round 1 already applies. When the raw PFI
/// solve's residual already sits at the floor, the correction's delta is
/// sub-floor noise; skipping it saves a full raw solve (LU + eta walk)
/// and the add pass. The probe counters (`ref0_at_floor_pct`) price
/// the opportunity; path-visible on exactly the at-floor solves ⇒ the
/// iteration gate adjudicates. Unset ⇒ round 0 unconditional (calibrated
/// path, no extra fp work on the default arm).
fn refine_gate_on() -> bool {
    use std::sync::OnceLock;
    static LEVER: OnceLock<bool> = OnceLock::new();
    *LEVER.get_or_init(|| {
        std::env::var("RUSOLVE_REFINE_GATE")
            .map(|v| parse_refine_gate(&v))
            .unwrap_or(false)
    })
}

/// Parse helper for [`refine_gate_on`]: "1" (any casing/trim) ⇒ on,
/// anything else off.
fn parse_refine_gate(value: &str) -> bool {
    value.trim() == "1"
}

/// Pivot-distress counter (a pivot below `pivot_tol` rejects the
/// eta file; the dual's retry path then refactors with reason `Pivot`).
fn rf_count_pivot_error() {
    if let Some(st) = pfi_stats() {
        st.rf_reason[RfReason::Pivot.index()].set(st.rf_reason[RfReason::Pivot.index()].get() + 1);
    }
}

/// `RUSOLVE_RFDROP=1` enables the per-refactor drift probe (needs
/// `RUSOLVE_DUAL_STATS=1`; measurement builds only, never the default path).
fn rf_drift_probe_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("RUSOLVE_RFDROP").is_ok())
}

/// Fixed sparse probe BTRAN rhs (4 unit nonzeros, deterministic): mirrors
/// the dual engine's `e_p` beta-BTRAN shape so the drift is measured on the
/// solve class the engine actually pays.
fn rf_probe_rhs(n: usize) -> Vec<f64> {
    thread_local! {
        static RHS: std::cell::RefCell<(usize, Vec<f64>)> =
            const { std::cell::RefCell::new((0, Vec::new())) };
    }
    RHS.with(|c| {
        let mut g = c.borrow_mut();
        if g.0 != n {
            let mut v = vec![0.0; n];
            let mut rng = crate::rng::XorShift64::new(0x5EED_F00D);
            for k in 1..=4.min(n) {
                v[n * k / 5] = rng.uniform(0.5, 2.0);
            }
            *g = (n, v);
        }
        g.1.clone()
    })
}

/// One-line dump of the PFI sub-phase counters (profiling).
pub(crate) fn pfi_stats_line() -> String {
    if std::env::var("RUSOLVE_DUAL_STATS").is_err() {
        return String::new();
    }
    PFI_STATS.with(|s| {
        let rf_total: u64 = s.rf_reason.iter().map(|c| c.get()).sum();
        let rf_cadence: u64 =
            s.rf_reason[RfReason::Clock.index()].get() + s.rf_reason[RfReason::Cap.index()].get();
        let mut reasons = String::new();
        for r in [
            RfReason::Clock,
            RfReason::Cap,
            RfReason::Pivot,
            RfReason::Verify,
            RfReason::Tail,
            RfReason::Forced,
        ] {
            let n = s.rf_reason[r.index()].get();
            if n > 0 {
                if !reasons.is_empty() {
                    reasons.push(' ');
                }
                reasons.push_str(&format!("{}={}", r.label(), n));
            }
        }
        let depth_mean = rf_total as f64;
        let line = format!(
            "rusolve pfi stats: n_raw={} n_corr={} lu={:.0}ms eta={:.0}ms refine={:.0}ms n_r1={} \
             resid={:.0}ms res_arms[hv,fwd,bwd_col,bwd_row]={:?}",
            s.n_raw.get(),
            s.n_corr.get(),
            s.lu_ns.get() / 1e6,
            s.eta_ns.get() / 1e6,
            s.refine_ns.get() / 1e6,
            s.n_r1.get(),
            s.resid_ns.get() / 1e6,
            [
                s.res_arm[0].get(),
                s.res_arm[1].get(),
                s.res_arm[2].get(),
                s.res_arm[3].get()
            ],
        );
        let mut cad = format!(
            "rusolve rf cadence: rf={} cadence={} reasons[{}] depth_mean={:.0} depth_max={} \
             eta_nnz_mean={:.0} solves_per_epoch={:.1} r1_per_epoch={:.2}",
            rf_total,
            rf_cadence,
            if reasons.is_empty() { "none" } else { &reasons },
            s.rf_depth_sum.get() as f64 / depth_mean.max(1.0),
            s.rf_depth_max.get(),
            s.rf_nnz_sum.get() as f64 / depth_mean.max(1.0),
            s.rf_solves_sum.get() as f64 / depth_mean.max(1.0),
            s.rf_r1_sum.get() as f64 / depth_mean.max(1.0),
        );
        if s.rf_drift_n.get() > 0 {
            cad.push_str(&format!(
                " drift_mean={:.2e} drift_max={:.2e} (n={})",
                s.rf_drift_sum.get() / s.rf_drift_n.get() as f64,
                s.rf_drift_max.get(),
                s.rf_drift_n.get(),
            ));
        }
        cad.push_str(&format!(
            " eta_bytes_mean={:.0} eta_bytes_max={}",
            s.rf_bytes_sum.get() as f64 / depth_mean.max(1.0),
            s.rf_bytes_max.get(),
        ));
        // Sweep-cost census line: per-eta-apply sparsity.
        let bt = s.bt_raw.get().max(1) as f64;
        let bt_etas = (s.bt_dot_sparse.get() + s.bt_dot_plain.get()).max(1) as f64;
        let touched = s.bt_touched.get().max(1) as f64;
        let ft_fired = (s.ft_fired.get() + s.ft_skipped.get()).max(1) as f64;
        let sweep = format!(
            "rusolve sweep census: bt_raw={} gate={:.0}% kept={:.0}% dot_sparse={:.0}% \
             touch={} eff={} eff_pct={:.1} hit_pct={:.0} y_nnz_mean={:.0} y_dense_pct={:.0} \
             ft_raw={} fired={} skip={} fired_pct={:.0} ft_work_mean={:.0} \
             refines_bt={} refines_ft={} refines_skip={} direct_raw={} \
             ref0_total={} ref0_at_floor_pct={:.0} ref0_ratio_mean={:.1e}",
            s.bt_raw.get(),
            100.0 * s.bt_gate_on.get() as f64 / bt,
            100.0 * s.bt_kept.get() as f64 / bt,
            100.0 * s.bt_dot_sparse.get() as f64 / bt_etas,
            s.bt_touched.get(),
            s.bt_effective.get(),
            100.0 * s.bt_effective.get() as f64 / touched,
            100.0 * s.bt_hits.get() as f64 / bt_etas,
            s.bt_y_nnz_sum.get() as f64 / s.bt_kept.get().max(1) as f64,
            100.0 * s.bt_y_dense.get() as f64 / bt,
            s.ft_raw.get(),
            s.ft_fired.get(),
            s.ft_skipped.get(),
            100.0 * s.ft_fired.get() as f64 / ft_fired,
            s.ft_work.get() as f64 / s.ft_fired.get().max(1) as f64,
            s.ref_bt.get(),
            s.ref_ft.get(),
            s.n_ref_skip.get(),
            s.n_raw.get().saturating_sub(s.n_corr.get()),
            s.ref0_total.get(),
            100.0 * s.ref0_at_floor.get() as f64 / s.ref0_total.get().max(1) as f64,
            s.ref0_ratio_sum.get() / s.ref0_total.get().max(1) as f64,
        );
        format!("{line}\n{cad}\n{sweep}")
    })
}

/// Row-wise mirror of sparse basis columns (entry `(column, value)`).
fn rows_of(cols: &[Vec<(usize, f64)>]) -> Vec<Vec<(u32, f64)>> {
    let mut rows: Vec<Vec<(u32, f64)>> = vec![Vec::new(); cols.len()];
    for (j, col) in cols.iter().enumerate() {
        for &(i, a) in col {
            rows[i].push((j as u32, a));
        }
    }
    rows
}

/// Ascending-sorted nonzeros of a dense column: exact zeros dropped. The
/// ascending row order is load-bearing (see [`ProductFormFactorization::refine`]
/// for the bit-exactness argument).
fn sparsify(col: &[f64]) -> Vec<(usize, f64)> {
    col.iter()
        .enumerate()
        .filter_map(|(i, &v)| (v != 0.0).then_some((i, v)))
        .collect()
}

/// Materialize dense columns from the sparse basis storage. Only for the rare
/// factorization sites (build + refactor); the per-solve hot path must stay
/// on the sparse [`ProductFormFactorization::basis_nz`].
fn dense_cols(cols_nz: &[Vec<(usize, f64)>], n: usize) -> Vec<Vec<f64>> {
    cols_nz
        .iter()
        .map(|nz| {
            let mut col = vec![0.0; n];
            for &(i, v) in nz {
                col[i] = v;
            }
            col
        })
        .collect()
}

/// Upper bound on |B[i][j]| over the sparse basis columns. Dropping exact
/// zeros cannot change the fold: |0.0| = 0.0 never raises a max that starts
/// at 0.0.
fn cols_max(cols_nz: &[Vec<(usize, f64)>]) -> f64 {
    cols_nz
        .iter()
        .flat_map(|c| c.iter().map(|&(_, v)| v.abs()))
        .fold(0.0f64, f64::max)
}

/// Basis factorization engine: the `solve`/`solve_transpose` method
/// names match, so PFI's hot call sites dispatch through the enum with zero
/// branching at the call sites themselves.
#[derive(Debug, Clone)]
enum Engine {
    Dense(DenseLu),
    Sparse(Box<SparseLu>),
}

impl Engine {
    fn solve(&self, b: &[f64], x: &mut [f64]) {
        match self {
            Engine::Dense(lu) => lu.solve(b, x),
            Engine::Sparse(lu) => lu.solve(b, x),
        }
    }
    fn solve_transpose(&self, b: &[f64], x: &mut [f64]) {
        match self {
            Engine::Dense(lu) => lu.solve_transpose(b, x),
            Engine::Sparse(lu) => lu.solve_transpose(b, x),
        }
    }
}

/// Largest basis dimension the sparse-Singular → dense fallback may
/// materialize. The dense path holds ~3 n×n f64 copies at once
/// (`dense_cols` + `DenseMatrix::from_columns` + `DenseLu::factor`'s
/// working copy): n = 4096 is ~400 MB and a few seconds of O(n³) work,
/// while the large LPs that reached it (irish-electricity, ex10 at
/// m ≈ 60k) took the process past 20 GB and never finished the dense LU.
/// Above the cap the sparse Singular
/// is returned to the caller, which already owns a Singular path (dual:
/// distress ⇒ primal fallback; primal: artificial repair). HiGHS never
/// densifies a basis; it handles rank deficiency by swapping logicals in
/// (`HFactor::buildHandleRankDeficiency`, util/HFactor.cpp, MIT) and by
/// backtracking (`HEkk::getNonsingularInverse`, simplex/HEkk.cpp).
/// `RUSOLVE_PFI_DENSE_MAX=<n>` overrides (lever; a huge value = legacy).
pub(crate) fn dense_fallback_max_n() -> usize {
    use std::sync::OnceLock;
    static CAP: OnceLock<usize> = OnceLock::new();
    *CAP.get_or_init(|| {
        std::env::var("RUSOLVE_PFI_DENSE_MAX")
            .ok()
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(4096)
    })
}

/// Factor `cols` with the chosen engine; a sparse Singular (Markowitz +
/// two-stage threshold can refuse where dense partial pivoting survives —
/// nasty bases exist, scsd1-class) falls back to dense STICKILY: the
/// caller flips `choice` and the safety net holds for the whole solve —
/// up to [`dense_fallback_max_n`]; above it the Singular propagates.
///
/// `rows` is the ascending row mirror of `cols` (or `None`), and `fws`
/// the persistent build workspace (see [`SparseLu::factor_sparse_ws`]);
/// both only remove allocation/transpose work — the factorization is
/// bit-identical to a cold `SparseLu::factor_sparse(cols, …)`.
fn factor_with(
    cols: &[Vec<(usize, f64)>],
    rows: Option<&[Vec<(u32, f64)>]>,
    config: &PfiConfig,
    choice: &mut PfiEngine,
    fws: &mut FactorWs,
) -> Result<Engine, FactorError> {
    if *choice == PfiEngine::Sparse {
        match SparseLu::factor_sparse_ws(cols, rows, config.pivot_tol, fws) {
            Ok(lu) => return Ok(Engine::Sparse(Box::new(lu))),
            Err(e @ FactorError::Singular { .. }) if cols.len() > dense_fallback_max_n() => {
                if std::env::var("RUSOLVE_PFI_DEBUG").is_ok() {
                    eprintln!(
                        "rusolve pfi: sparse LU {e} (n={}, no dense fallback)",
                        cols.len()
                    );
                }
                return Err(e);
            }
            Err(FactorError::Singular { .. }) => {
                *choice = PfiEngine::Dense;
            }
            Err(e) => return Err(e),
        }
    }
    // Dense path (explicit engine choice or sticky Singular fallback): the
    // only place the O(n²) column materialization is still paid.
    let dense = dense_cols(cols, cols.len());
    DenseLu::factor(&DenseMatrix::from_columns(&dense), config.pivot_tol).map(Engine::Dense)
}

/// Cost of one eta entry in a solve relative to one unit of
/// [`SparseLu::build_work`] (~10 ns, measured): a gathered multiply-add is
/// ~1.5 ns.
const ETA_TICK: f64 = 0.15;

/// Knuth exact two-sum (TAOCP Vol 2, Alg. 3.4; Dekker 1971): returns
/// (hi, lo) with hi + lo == a + b exactly.
#[inline]
pub(crate) fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

impl ProductFormFactorization {
    /// Factorizes an initial basis. Returns [`FactorError::Singular`] if the
    /// basis is numerically singular.
    pub fn build(initial_basis: &DenseMatrix, config: PfiConfig) -> Result<Self, FactorError> {
        let basis_nz: Vec<Vec<(usize, f64)>> = (0..initial_basis.n())
            .map(|j| sparsify(initial_basis.col(j)))
            .collect();
        Self::build_sparse(basis_nz, config)
    }

    /// Like [`build`](Self::build) but from sparse basis columns
    /// (`(row, value)` ascending by row, exact zeros omitted): skips the
    /// O(n²) dense materialization a warm start or MIP node re-solve would
    /// otherwise pay on every basis installation.
    pub fn build_sparse(
        basis_nz: Vec<Vec<(usize, f64)>>,
        config: PfiConfig,
    ) -> Result<Self, FactorError> {
        let n = basis_nz.len();
        let col_max = cols_max(&basis_nz);
        let mut engine_choice = config.engine;
        let basis_rows = rows_of(&basis_nz);
        let mut fws = FactorWs::default();
        let lu = factor_with(
            &basis_nz,
            Some(&basis_rows),
            &config,
            &mut engine_choice,
            &mut fws,
        )?;
        Ok(Self {
            n,
            lu,
            engine_choice,
            etas: Vec::new(),
            epoch: 0,
            eta_nnz: 0,
            eta_ticks: std::cell::Cell::new(0.0),
            refactor_nanos: std::cell::Cell::new(0.0),
            basis_rows,
            rows_ok: true,
            basis_nz,
            col_max,
            config,
            ws: SolveWorkspace::default(),
            fws,
        })
    }

    pub fn n(&self) -> usize {
        self.n
    }

    /// Count of refactorizations since build (see the `epoch` field).
    pub fn factor_epoch(&self) -> u64 {
        self.epoch
    }

    /// Number of eta matrices currently in the file.
    pub fn eta_count(&self) -> usize {
        self.etas.len()
    }

    /// Solves `B x = b` for the *current* basis: dense-LU solve with `B0`,
    /// then the eta forward pass oldest → newest
    /// (`x[i] = y[i] − η_i t`, `x[pos] = t` with `t = y[pos] / pivot`).
    ///
    /// # Panics
    ///
    /// Panics if `b.len()` or `x.len()` differ from `n`.
    pub fn solve_b(&self, b: &[f64], x: &mut [f64]) {
        self.solve_b_raw(b, x);
        self.refine(b, x, false);
    }

    /// Unrefined forward solve (LU + eta passes). Public for the
    /// refinement-gating experiment: pricing correctness comes from
    /// BTRAN; FTRAN b̄ mostly steers ratio-test choices, where mild drift is
    /// survivable. Bench decides whether the 26/28 correctness survives.
    pub fn solve_b_raw(&self, b: &[f64], x: &mut [f64]) {
        let st = pfi_stats();
        let t0 = st.as_ref().map(|_| std::time::Instant::now());
        self.eta_ticks
            .set(self.eta_ticks.get() + ETA_TICK * (self.eta_nnz + self.etas.len()) as f64);
        self.lu.solve(b, x);
        if let (Some(s), Some(t)) = (st.as_ref(), t0) {
            s.lu_ns.set(s.lu_ns.get() + t.elapsed().as_nanos() as f64);
        }
        let t1 = st.as_ref().map(|_| std::time::Instant::now());
        let stat = st.as_ref();
        for eta in self.etas.iter() {
            let t = x[eta.pos] / eta.pivot;
            if t != 0.0 {
                if let Some(s) = stat {
                    s.ft_fired.set(s.ft_fired.get() + 1);
                    s.ft_work.set(s.ft_work.get() + eta.idx.len() as u64);
                }
                for (&i, &v) in eta.idx.iter().zip(&eta.val) {
                    let i = i as usize;
                    x[i] -= v * t;
                }
            } else if let Some(s) = stat {
                s.ft_skipped.set(s.ft_skipped.get() + 1);
            }
            x[eta.pos] = t;
        }
        if let (Some(s), Some(t)) = (st.as_ref(), t1) {
            s.eta_ns.set(s.eta_ns.get() + t.elapsed().as_nanos() as f64);
            s.n_raw.set(s.n_raw.get() + 1);
            s.ft_raw.set(s.ft_raw.get() + 1);
        }
    }

    /// One or two residual-correction rounds on a PFI solve, with the
    /// residual evaluated in DOUBLE-WORD (hi+lo) arithmetic: FMA product
    /// errors (`a.mul_add(-b, a*b)` is exact) plus Knuth two-sum
    /// accumulation (Dekker 1971; Hida–Bailey–Arvanitidis double-word
    /// arithmetic; residual-accurate refinement per Maddick–Ghysels–Hogg
    /// practice — plain-precision residuals cap refinement at κ·eps and
    /// buy nothing, measured on refine_recovers_ill_conditioned_solve).
    /// x ← x + B̃⁻¹(rhs − B·x) with residuals against the EXPLICIT current
    /// basis columns, so drifted etas (B̃ ≠ B) get pulled back toward the
    /// true B. Round 0 corrects unconditionally (for ill-conditioned
    /// systems the plain residual sits at round-off while forward error is
    /// κ times larger); round 1 runs only while the residual clears the fp
    /// floor. scsd6's phase-2 wrong-reduced-cost 2-cycle is the motivating
    /// failure (Optimal on-bench even with the plain residual).
    ///
    /// # Sparse residuals are bit-exact
    ///
    /// The residual sweeps ([`residual_bwd`]/[`residual_fwd`]) run per solve
    /// and iterate only the nonzeros of each basis column instead of all m
    /// entries. The old dense loops skipped `a == 0.0` (and `xi == 0.0`)
    /// entries; such a skipped term contributes p = ±0.0, pe = ∓0.0 and
    /// two_sum(s, ±0.0) = (s, ±0.0) up to zero signs, i.e. it is a no-op of
    /// the double-word accumulation (zero signs cancel in the final
    /// rounding). Iterating the column's nonzeros in ASCENDING ROW ORDER,
    /// columns in ascending column order, and keeping the `xi == 0.0`
    /// column-level skip therefore reproduces the identical fp operation
    /// sequence ⇒ bit-identical residual, identical pivot path. Pinned by
    /// `sparse_residual_bit_exact`.
    ///
    /// [`residual_bwd`]: ProductFormFactorization::residual_bwd
    /// [`residual_fwd`]: ProductFormFactorization::residual_fwd
    fn refine(&self, rhs: &[f64], x: &mut [f64], transpose: bool) {
        let m = self.n;
        let st = pfi_stats();
        if let Some(s) = st.as_ref() {
            if transpose {
                s.ref_bt.set(s.ref_bt.get() + 1);
            } else {
                s.ref_ft.set(s.ref_ft.get() + 1);
            }
        }
        // Candidate (c) young-file skip (see `refine_age_frac`): unset ⇒
        // no branch cost beyond one OnceLock-cached Option read.
        if let Some(f) = refine_age_frac() {
            let young = match &self.lu {
                Engine::Sparse(lu) => self.eta_ticks.get() < f * lu.build_work() as f64,
                Engine::Dense(_) => false,
            };
            if young {
                if let Some(s) = st.as_ref() {
                    s.n_ref_skip.set(s.n_ref_skip.get() + 1);
                }
                return;
            }
        }
        let t0 = st.as_ref().map(|_| std::time::Instant::now());
        'outer: {
            // Persistent scratch (see SolveWorkspace's bit-exactness note —
            // every buffer is fully rewritten by the sweeps below each round, so
            // reuse == the old per-call `vec![0.0; m]`). The four borrows do not
            // collide with the correction round's raw solves, which use only the
            // `y` buffer and the LU engine's disjoint scratch.
            let mut r = self.ws.buf(&self.ws.r, m);
            let mut d = self.ws.buf(&self.ws.d, m);
            let mut hi_part = self.ws.buf(&self.ws.hi_part, m);
            let mut lo_part = self.ws.buf(&self.ws.lo_part, m);
            for round in 0..2 {
                // r ← rhs − B·x (forward) or rhs − Bᵀ·x (transpose), double-word.
                let tr = st.as_ref().map(|_| std::time::Instant::now());
                if transpose {
                    self.residual_bwd(rhs, x, &mut r, &mut hi_part, &mut lo_part);
                } else {
                    self.residual_fwd(rhs, x, &mut r, &mut hi_part, &mut lo_part);
                }
                if let (Some(s), Some(t)) = (st.as_ref(), tr) {
                    s.resid_ns
                        .set(s.resid_ns.get() + t.elapsed().as_nanos() as f64);
                }
                // Floor test. Round 1 always applies it (calibrated
                // convergence); round 0 only under the `RUSOLVE_REFINE_GATE`
                // lever (candidate c′) or the stats probe — unset and
                // stats-off, round 0 runs no extra fp work, bit-identical.
                let gate0 = round == 0 && (refine_gate_on() || st.is_some());
                if round > 0 || gate0 {
                    let rnorm = r.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                    let xnorm = x.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                    let rhs_norm = rhs.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                    let floor = (rhs_norm + xnorm * self.col_max * (m as f64)) * 1e-13;
                    if round == 0 {
                        if let Some(s) = st.as_ref() {
                            s.ref0_total.set(s.ref0_total.get() + 1);
                            if rnorm <= floor {
                                s.ref0_at_floor.set(s.ref0_at_floor.get() + 1);
                            }
                            s.ref0_ratio_sum.set(
                                s.ref0_ratio_sum.get() + (rnorm / floor.max(1e-300)).min(1e30),
                            );
                        }
                        if rnorm <= floor && refine_gate_on() {
                            if let Some(s) = st.as_ref() {
                                s.n_ref_skip.set(s.n_ref_skip.get() + 1);
                            }
                            break 'outer;
                        }
                    } else if rnorm <= floor {
                        break 'outer;
                    }
                }
                if transpose {
                    self.solve_bt_raw(&r, &mut d);
                } else {
                    self.solve_b_raw(&r, &mut d);
                }
                if let Some(s) = st.as_ref() {
                    s.n_corr.set(s.n_corr.get() + 1);
                    if round > 0 {
                        // A round-1 correction means the round-0
                        // residual cleared the fp floor — eta-file drift
                        // the double-word residual could not ignore.
                        s.n_r1.set(s.n_r1.get() + 1);
                    }
                }
                let mut dn = 0.0f64;
                for (xi, &di) in x.iter_mut().zip(d.iter()) {
                    *xi += di;
                    dn = dn.max(di.abs());
                }
                // Convergence shortcut: a correction below 1e-8·‖x‖ is noise —
                // the round-1 DW mat-vec would only re-measure it. Saves the
                // second residual pass on every well-conditioned solve (the
                // common case; scsd-class corrections are orders of magnitude
                // larger and proceed to round 1).
                let xn = x.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                if dn <= 1e-8 * xn.max(1.0) {
                    break 'outer;
                }
            }
        }
        if let (Some(s), Some(t)) = (st.as_ref(), t0) {
            s.refine_ns
                .set(s.refine_ns.get() + t.elapsed().as_nanos() as f64);
        }
    }

    /// Transpose residual `r_j = rhs_j − Σ_i B[i][j]·x[i]`, double-word, over
    /// the sparse basis columns. See the bit-exactness note on [`refine`].
    ///
    /// [`refine`]: ProductFormFactorization::refine
    ///
    /// Row-wise over `basis_rows`: only the nonzeros of `x` are visited, in
    /// ascending `i`, so each column accumulator `(s_j, c_j)` receives the
    /// same terms in the same (ascending-row) order as the column-wise dense
    /// loop ⇒ bit-identical. Columns no nonzero of `x` reaches have
    /// `s = c = 0` ⇒ `r_j = rhs_j` exactly.
    fn residual_bwd(
        &self,
        rhs: &[f64],
        x: &[f64],
        r: &mut [f64],
        _s_part: &mut [f64],
        _c_part: &mut [f64],
    ) {
        let arm = |k: usize| {
            if let Some(s) = pfi_stats() {
                s.res_arm[k].set(s.res_arm[k].get() + 1);
            }
        };
        let m = x.len();
        let mut idx = self.ws.res_idx.borrow_mut();
        // ONE scan of x gathers its ascending nonzeros (early exit at
        // the density cutoff) — the old arm counted, then rescanned x twice.
        if !self.rows_ok || !collect_nz(x, &mut idx, m / 4) {
            arm(2);
            // Dense x (dual simplex, or a mirror out of sync): the plain
            // column-wise loop is cheaper than two row-wise passes plus
            // accumulator fills; same terms in the same order.
            for (j, col) in self.basis_nz.iter().enumerate() {
                let mut s = 0.0;
                let mut c = 0.0;
                for &(i, a) in col {
                    let xi = x[i];
                    if xi == 0.0 {
                        continue;
                    }
                    let p = a * xi;
                    let pe = a.mul_add(-xi, p);
                    let (hi, lo) = two_sum(s, p);
                    s = hi;
                    c += lo - pe;
                }
                let (hi, lo) = two_sum(rhs[j], -s);
                r[j] = hi + (lo - c);
            }
            return;
        }
        arm(3);
        // Accumulators are the all-zero `sp_hi`/`sp_lo` (reset over the
        // touched columns below) instead of two O(m) fills per call.
        let mut s_part = self.ws.buf(&self.ws.sp_hi, m);
        let mut c_part = self.ws.buf(&self.ws.sp_lo, m);
        for &i in idx.iter() {
            let i = i as usize;
            let xi = x[i];
            for &(j, a) in &self.basis_rows[i] {
                let j = j as usize;
                let p = a * xi;
                let pe = a.mul_add(-xi, p); // = p − a·xi exactly
                let (hi, lo) = two_sum(s_part[j], p);
                s_part[j] = hi;
                c_part[j] += lo - pe;
            }
        }
        r.copy_from_slice(rhs);
        let (mut mark, stamp) = self.ws.next_stamp(m);
        for &i in idx.iter() {
            for &(j, _) in &self.basis_rows[i as usize] {
                let j = j as usize;
                if mark[j] == stamp {
                    continue;
                }
                mark[j] = stamp;
                // r_j = rhs_j − (s + c), rounded ONCE at the end so the
                // error scales with |r_j| (tiny), not |dot| (huge).
                let (hi, lo) = two_sum(rhs[j], -s_part[j]);
                r[j] = hi + (lo - c_part[j]);
                s_part[j] = 0.0;
                c_part[j] = 0.0;
            }
        }
    }

    /// Forward residual `r ← rhs − Σ_j B[:,j]·x[j]`, double-word, over the
    /// sparse basis columns: row accumulators touched only where the column
    /// has nonzeros, each row still accumulated in the same
    /// (column-ascending) order as the old dense loop. See the bit-exactness
    /// note on [`refine`].
    ///
    /// One scan of `x` gathers its nonzeros (the old sweep walked all m
    /// column headers, then counted x, then — sparse case — walked the
    /// headers again); for a sparse `x` the row accumulators are the
    /// all-zero `sp_hi`/`sp_lo`, seeded with `rhs_i` on first touch (the
    /// value the old full `hi_part = rhs` copy gave them), so the
    /// accumulation sequence per row is unchanged.
    ///
    /// [`refine`]: ProductFormFactorization::refine
    fn residual_fwd(
        &self,
        rhs: &[f64],
        x: &[f64],
        r: &mut [f64],
        hi_part: &mut [f64],
        lo_part: &mut [f64],
    ) {
        if let Some(s) = pfi_stats() {
            s.res_arm[1].set(s.res_arm[1].get() + 1);
        }
        let m = x.len();
        let mut idx = self.ws.res_idx.borrow_mut();
        collect_nz(x, &mut idx, usize::MAX);
        if !self.rows_ok || idx.len() * 4 > m {
            hi_part.copy_from_slice(rhs);
            lo_part.fill(0.0);
            for &j in idx.iter() {
                let j = j as usize;
                let xj = x[j];
                for &(i, a) in &self.basis_nz[j] {
                    let (si, ci) = (&mut hi_part[i], &mut lo_part[i]);
                    let p = a * xj;
                    let pe = a.mul_add(-xj, p); // = p − a·xj exactly
                    let (hi, lo) = two_sum(*si, -p);
                    *si = hi;
                    *ci += lo + pe; // subtract a·xj = p − pe
                }
            }
            for (ri, (hi, lo)) in r.iter_mut().zip(hi_part.iter().zip(lo_part.iter())) {
                *ri = *hi + *lo;
            }
            return;
        }
        let mut hi_s = self.ws.buf(&self.ws.sp_hi, m);
        let mut lo_s = self.ws.buf(&self.ws.sp_lo, m);
        let mut touched = self.ws.res_touch.borrow_mut();
        touched.clear();
        let (mut mark, stamp) = self.ws.next_stamp(m);
        for &j in idx.iter() {
            let j = j as usize;
            let xj = x[j];
            for &(i, a) in &self.basis_nz[j] {
                if mark[i] != stamp {
                    mark[i] = stamp;
                    hi_s[i] = rhs[i];
                    touched.push(i as u32);
                }
                let p = a * xj;
                let pe = a.mul_add(-xj, p);
                let (hi, lo) = two_sum(hi_s[i], -p);
                hi_s[i] = hi;
                lo_s[i] += lo + pe;
            }
        }
        r.copy_from_slice(rhs);
        for &i in touched.iter() {
            let i = i as usize;
            r[i] = hi_s[i] + lo_s[i];
            hi_s[i] = 0.0;
            lo_s[i] = 0.0;
        }
    }

    /// Solves `Bᵀ x = b` for the *current* basis: since
    /// `B⁻ᵀ = B0⁻ᵀ E1⁻ᵀ … Em⁻ᵀ`, the eta transpose pass runs newest-first
    /// (only `x[pos]` changes per eta: `x[pos] = (x[pos] − Σ η_i x[i]) /
    /// pivot`) and the dense-LU transpose solve with `B0ᵀ` applies last.
    ///
    /// # Panics
    ///
    /// Panics if `b.len()` or `x.len()` differ from `n`.
    pub fn solve_bt(&self, b: &[f64], x: &mut [f64]) {
        self.solve_bt_raw(b, x);
        self.refine(b, x, true);
    }

    /// Unrefined transpose solve (LU + eta passes) — the low-accuracy
    /// regime's BTRAN (`RUSOLVE_DUAL_ACC`), mirroring [`Self::solve_b_raw`].
    pub(crate) fn solve_bt_raw(&self, b: &[f64], x: &mut [f64]) {
        assert_eq!(b.len(), self.n, "solve_bt: wrong rhs length");
        assert_eq!(x.len(), self.n, "solve_bt: wrong solution length");
        let st = pfi_stats();
        let t0 = st.as_ref().map(|_| std::time::Instant::now());
        // Persistent eta-sweep scratch; `copy_from_slice` overwrites
        // all n entries (the old `vec![0.0; n]` zero-fill was a dead store).
        self.eta_ticks
            .set(self.eta_ticks.get() + ETA_TICK * (self.eta_nnz + self.etas.len()) as f64);
        let mut y = self.ws.buf(&self.ws.y, self.n);
        y.copy_from_slice(b);
        // BTRAN of a sparse rhs (the dual engine's `e_p` → β, HiGHS
        // `row_ep = B⁻ᵀe_p`, HEkkDual.cpp:2143 area) walks the eta file
        // NEWEST-FIRST, and each eta's dot product `Σ η_i·y[i]` iterates
        // the eta's full index list although `y` carries at most
        // nnz(rhs) + (processed eta count) nonzeros — the entering column's
        // b̄ can be dense while the BTRAN working vector stays tiny (qap10:
        // 65%-dense etas, ≤ 21-nonzero y). Track y's nonzero positions in
        // an ascending list and compute the dot over the INTERSECTION with
        // the eta's (also ascending) index list. Bit-exact: both walks
        // visit the shared entries in ascending order, and skipped entries
        // have y[i] == 0.0, whose products ±0.0 never change an
        // accumulation (x + ±0.0 = x for x ≠ 0; ±0.0 + ∓0.0 = +0.0).
        // Hall & McKinnon 2005 is the hyper-sparse reference; this is our
        // PFI analogue of their reachable-set bookkeeping. Tracking is
        // abandoned once the pattern grows past `n/16` (the SparseLu
        // hyper threshold) and never ATTEMPTED unless the eta file is
        // dense relative to n (a thin file's plain pass is already cheap;
        // the O(n) rhs scan would be pure overhead — krka regressed 24→64ms
        // before this gate).
        let pat_max = self.n / 16;
        // Sweep census (stats-only): was the walk attempted / did the
        // pattern survive the pass? See `PfiStats` field docs.
        let stat = st.as_ref();
        let mut attempted = false;
        let mut kept = false;
        // Tracking only pays on a DENSE eta file (plain pass ≥ ~4n work to
        // cut) — a thin file (krka: network LP, m=12.5k, LU nnz ≈ 4k, sparse
        // etas) makes the O(n) rhs scan pure loss, and a dense rhs hits the
        // early-exit below after ≤ pat_max+1 steps.
        let mut pat = if bt_gate_open(self.eta_nnz, self.etas.len(), self.n, bt_gate_ratio()) {
            attempted = true;
            let mut p = self.ws.bt_pat.borrow_mut();
            p.clear();
            let mut overflow = false;
            for (i, &v) in y.iter().enumerate() {
                if v != 0.0 {
                    p.push(i as u32);
                    if p.len() > pat_max {
                        overflow = true;
                        break;
                    }
                }
            }
            if overflow {
                None
            } else {
                p.reserve(self.etas.len() + 8);
                kept = true;
                Some(p)
            }
        } else {
            None
        };
        for eta in self.etas.iter().rev() {
            // An eta can carry an EMPTY idx (unit-column pivot, no off-pivot
            // entries; krka exercises this) — its dot is 0.0 with zero work,
            // and `ilog2(0)` would panic, so gate on a non-empty index list.
            let idx_len = eta.idx.len();
            let sparse_walk = pat
                .as_ref()
                .is_some_and(|p| idx_len > 0 && p.len() * (idx_len.ilog2() as usize + 2) < idx_len);
            let dot = if sparse_walk {
                let p = pat.as_ref().unwrap();
                // Pattern walk: iterate y's nonzeros (dead entries skipped),
                // binary-search each in the eta's ascending index list with
                // a monotone window (both lists ascend).
                let mut dot = 0.0;
                let mut matched = 0u64;
                let idx = &eta.idx;
                let mut lo = 0usize;
                for &k in p.iter() {
                    let yk = y[k as usize];
                    if yk == 0.0 {
                        continue;
                    }
                    let seg = &idx[lo..];
                    let pos = seg.partition_point(|&j| (j as u32) < k);
                    if pos < seg.len() && seg[pos] as u32 == k {
                        dot += eta.val[lo + pos] * yk;
                        lo += pos + 1;
                        matched += 1;
                    }
                }
                if let Some(s) = stat {
                    s.bt_dot_sparse.set(s.bt_dot_sparse.get() + 1);
                    s.bt_touched.set(s.bt_touched.get() + p.len() as u64);
                    s.bt_effective.set(s.bt_effective.get() + matched);
                }
                dot
            } else if let Some(s) = stat {
                // Stats arm: identical fp sequence to the plain loop (the
                // hoisted yk is the same load), plus work/effect counters.
                let mut dot = 0.0;
                let mut eff = 0u64;
                for (&i, &v) in eta.idx.iter().zip(&eta.val) {
                    let yk = y[i as usize];
                    dot += v * yk;
                    eff += (yk != 0.0) as u64;
                }
                s.bt_dot_plain.set(s.bt_dot_plain.get() + 1);
                s.bt_touched.set(s.bt_touched.get() + idx_len as u64);
                s.bt_effective.set(s.bt_effective.get() + eff);
                dot
            } else {
                let mut dot = 0.0;
                for (&i, &v) in eta.idx.iter().zip(&eta.val) {
                    dot += v * y[i as usize];
                }
                dot
            };
            y[eta.pos] = (y[eta.pos] - dot) / eta.pivot;
            if let Some(s) = stat {
                s.bt_hits.set(s.bt_hits.get() + (dot != 0.0) as u64);
            }
            let mut abandon = false;
            if let Some(p) = pat.as_mut() {
                // This eta may have introduced a nonzero at its pos; drop
                // tracking entirely if the pattern outgrew the threshold.
                let pos32 = eta.pos as u32;
                if y[eta.pos] != 0.0 {
                    let ip = p.partition_point(|&j| j < pos32);
                    if ip == p.len() || p[ip] != pos32 {
                        p.insert(ip, pos32);
                    }
                }
                abandon = p.len() > pat_max;
            }
            if abandon {
                pat = None;
                kept = false;
            }
        }
        let y_live = pat
            .as_ref()
            .map(|p| p.iter().filter(|&&k| y[k as usize] != 0.0).count() as u64);
        drop(pat);
        if let Some(s) = stat {
            s.bt_raw.set(s.bt_raw.get() + 1);
            if attempted {
                s.bt_gate_on.set(s.bt_gate_on.get() + 1);
            }
            if kept {
                s.bt_kept.set(s.bt_kept.get() + 1);
            }
            match y_live {
                Some(nz) => s.bt_y_nnz_sum.set(s.bt_y_nnz_sum.get() + nz),
                None => s.bt_y_dense.set(s.bt_y_dense.get() + 1),
            }
        }
        if let (Some(s), Some(t)) = (st.as_ref(), t0) {
            s.eta_ns.set(s.eta_ns.get() + t.elapsed().as_nanos() as f64);
            s.n_raw.set(s.n_raw.get() + 1);
        }
        let t1 = st.as_ref().map(|_| std::time::Instant::now());
        self.lu.solve_transpose(&y, x);
        if let (Some(s), Some(t)) = (st.as_ref(), t1) {
            s.lu_ns.set(s.lu_ns.get() + t.elapsed().as_nanos() as f64);
        }
    }

    // ---------------------------------------------------------------
    // Index-carrying solve stack.
    //
    // Twins of solve_b / solve_bt / their raw forms / refine / update that
    // take and return [`HVec`]s, so that on a hyper-sparse iteration no
    // pass of the stack costs O(m): the LU entry seeds from the rhs index
    // instead of scanning it, the result arrives with its pattern instead
    // of `fill(0)` + scatter, the FTRAN eta pass records the positions it
    // writes, the BTRAN eta pass keeps the y pattern for free (no O(m) rhs
    // scan, so the thin-file gate is not needed), and the double-word
    // refinement residuals touch only the rows/columns the pattern
    // reaches. Every fp operation sequence is the one of the slice
    // version (ascending visiting orders preserved, skipped terms are ±0.0
    // no-ops — the same argument as the sparse residual sweeps), pinned by
    // `hv_stack_bit_exact`. Hall & McKinnon (2005) §3 is the algorithmic
    // reference; HiGHS `HVector` + `HFactor::ftran/btran` + `HEkk`'s
    // per-iteration vectors (MIT, read as reference) the engineering one.
    // A dense-mode operand falls back to the slice code.
    // ---------------------------------------------------------------

    /// [`Self::solve_b`] (`refine = true`) or [`Self::solve_b_raw`]
    /// on index-carrying vectors. Bit-identical values.
    pub(crate) fn solve_b_hv(&self, b: &HVec, x: &mut HVec, refine: bool) {
        self.solve_b_raw_hv(b, x);
        if refine {
            self.refine_hv(b, x, false);
        }
    }

    /// [`Self::solve_bt`] (`refine = true`) or [`Self::solve_bt_raw`]
    /// on index-carrying vectors. Bit-identical values.
    pub(crate) fn solve_bt_hv(&self, b: &HVec, x: &mut HVec, refine: bool) {
        self.solve_bt_raw_hv(b, x);
        if refine {
            self.refine_hv(b, x, true);
        }
    }

    fn solve_b_raw_hv(&self, b: &HVec, x: &mut HVec) {
        assert_eq!(b.len(), self.n, "solve_b: wrong rhs length");
        assert_eq!(x.len(), self.n, "solve_b: wrong solution length");
        let st = pfi_stats();
        let t0 = st.as_ref().map(|_| std::time::Instant::now());
        self.eta_ticks
            .set(self.eta_ticks.get() + ETA_TICK * (self.eta_nnz + self.etas.len()) as f64);
        match &self.lu {
            Engine::Sparse(lu) => lu.solve_hv(b, x),
            Engine::Dense(lu) => lu.solve(b.vals(), x.dense_mut()),
        }
        if let (Some(s), Some(t)) = (st.as_ref(), t0) {
            s.lu_ns.set(s.lu_ns.get() + t.elapsed().as_nanos() as f64);
        }
        let t1 = st.as_ref().map(|_| std::time::Instant::now());
        let stat = st.as_ref();
        for eta in self.etas.iter() {
            let t = x[eta.pos] / eta.pivot;
            if t != 0.0 {
                if let Some(s) = stat {
                    s.ft_fired.set(s.ft_fired.get() + 1);
                    s.ft_work.set(s.ft_work.get() + eta.idx.len() as u64);
                }
                if x.is_dense() {
                    let xv = x.dense_mut();
                    for (&i, &v) in eta.idx.iter().zip(&eta.val) {
                        xv[i as usize] -= v * t;
                    }
                } else {
                    for (&i, &v) in eta.idx.iter().zip(&eta.val) {
                        *x.entry(i as usize) -= v * t;
                    }
                }
            } else if let Some(s) = stat {
                s.ft_skipped.set(s.ft_skipped.get() + 1);
            }
            x.assign(eta.pos, t);
        }
        if let (Some(s), Some(t)) = (st.as_ref(), t1) {
            s.eta_ns.set(s.eta_ns.get() + t.elapsed().as_nanos() as f64);
            s.n_raw.set(s.n_raw.get() + 1);
            s.ft_raw.set(s.ft_raw.get() + 1);
        }
    }

    fn solve_bt_raw_hv(&self, b: &HVec, x: &mut HVec) {
        assert_eq!(b.len(), self.n, "solve_bt: wrong rhs length");
        assert_eq!(x.len(), self.n, "solve_bt: wrong solution length");
        if b.is_dense() {
            self.solve_bt_raw(b.vals(), x.dense_mut());
            return;
        }
        let st = pfi_stats();
        let t0 = st.as_ref().map(|_| std::time::Instant::now());
        self.eta_ticks
            .set(self.eta_ticks.get() + ETA_TICK * (self.eta_nnz + self.etas.len()) as f64);
        let mut y = self.ws.y_hv.borrow_mut();
        y.ensure_len(self.n);
        y.copy_from(b);
        y.sort();
        let stat = st.as_ref();
        for eta in self.etas.iter().rev() {
            // Same per-eta walk choice as `solve_bt_raw` (the choice
            // never changes bits — both walks visit the shared nonzeros in
            // ascending order).
            let idx_len = eta.idx.len();
            let walk = match y.pattern() {
                Some(p) => idx_len > 0 && p.len() * (idx_len.ilog2() as usize + 2) < idx_len,
                None => false,
            };
            let dot = if walk {
                let p = y.pattern().unwrap();
                let yv = y.vals();
                let idx = &eta.idx;
                let mut dot = 0.0;
                let mut lo = 0usize;
                for &k in p {
                    let yk = yv[k as usize];
                    if yk == 0.0 {
                        continue;
                    }
                    let seg = &idx[lo..];
                    let pos = seg.partition_point(|&j| (j as u32) < k);
                    if pos < seg.len() && seg[pos] as u32 == k {
                        dot += eta.val[lo + pos] * yk;
                        lo += pos + 1;
                    }
                }
                if let Some(s) = stat {
                    s.bt_dot_sparse.set(s.bt_dot_sparse.get() + 1);
                }
                dot
            } else {
                let yv = y.vals();
                let mut dot = 0.0;
                for (&i, &v) in eta.idx.iter().zip(&eta.val) {
                    dot += v * yv[i as usize];
                }
                if let Some(s) = stat {
                    s.bt_dot_plain.set(s.bt_dot_plain.get() + 1);
                }
                dot
            };
            let nv = (y[eta.pos] - dot) / eta.pivot;
            if nv != 0.0 {
                y.touch_sorted(eta.pos);
            }
            y.assign(eta.pos, nv);
        }
        if let (Some(s), Some(t)) = (st.as_ref(), t0) {
            s.eta_ns.set(s.eta_ns.get() + t.elapsed().as_nanos() as f64);
            s.n_raw.set(s.n_raw.get() + 1);
            s.bt_raw.set(s.bt_raw.get() + 1);
        }
        let t1 = st.as_ref().map(|_| std::time::Instant::now());
        match &self.lu {
            Engine::Sparse(lu) => lu.solve_transpose_hv(&y, x),
            Engine::Dense(lu) => lu.solve_transpose(y.vals(), x.dense_mut()),
        }
        if let (Some(s), Some(t)) = (st.as_ref(), t1) {
            s.lu_ns.set(s.lu_ns.get() + t.elapsed().as_nanos() as f64);
        }
    }

    /// Double-word residual over the reachable rows/columns only.
    /// Returns `false` (caller runs the dense sweep) when `x` or `rhs` is
    /// dense, the row mirror is stale (transpose), or the touched set could
    /// outgrow `r`'s pattern capacity. Bit-identical to
    /// [`Self::residual_fwd`] / [`Self::residual_bwd`]: same ascending
    /// visiting order of `x`'s nonzeros, each accumulator starting from the
    /// value the dense sweep starts it from (`rhs_i` resp. 0.0 — the
    /// zero-invariant scratch plus the rhs pattern give exactly that), and
    /// an untouched rhs-only entry finalizes to `rhs_i` in both.
    fn residual_hv(
        &self,
        rhs: &HVec,
        x: &mut HVec,
        r: &mut HVec,
        hi: &mut [f64],
        lo: &mut [f64],
        transpose: bool,
    ) -> bool {
        if x.is_dense() || rhs.is_dense() || (transpose && !self.rows_ok) {
            return false;
        }
        x.sort();
        let xp = x.pattern().unwrap();
        let xv = x.vals();
        let rp = rhs.pattern().unwrap();
        let rv = rhs.vals();
        let mut work = rp.len();
        for &i in xp {
            let i = i as usize;
            if xv[i] != 0.0 {
                work += if transpose {
                    self.basis_rows[i].len()
                } else {
                    self.basis_nz[i].len()
                };
            }
        }
        if work > r.cap() {
            return false;
        }
        r.clear();
        for &i in rp {
            r.set(i as usize, rv[i as usize]);
        }
        if transpose {
            for &i in xp {
                let i = i as usize;
                let xi = xv[i];
                if xi == 0.0 {
                    continue;
                }
                for &(j, a) in &self.basis_rows[i] {
                    let j = j as usize;
                    r.touch(j);
                    let p = a * xi;
                    let pe = a.mul_add(-xi, p);
                    let (h, l) = two_sum(hi[j], p);
                    hi[j] = h;
                    lo[j] += l - pe;
                }
            }
            for k in 0..r.count() {
                let j = r.pattern().unwrap()[k] as usize;
                let (h, l) = two_sum(rv[j], -hi[j]);
                r.put_listed(j, h + (l - lo[j]));
                hi[j] = 0.0;
                lo[j] = 0.0;
            }
        } else {
            for &i in rp {
                hi[i as usize] = rv[i as usize];
            }
            for &j in xp {
                let j = j as usize;
                let xj = xv[j];
                if xj == 0.0 {
                    continue;
                }
                for &(i, a) in &self.basis_nz[j] {
                    r.touch(i);
                    let p = a * xj;
                    let pe = a.mul_add(-xj, p);
                    let (h, l) = two_sum(hi[i], -p);
                    hi[i] = h;
                    lo[i] += l + pe;
                }
            }
            for k in 0..r.count() {
                let i = r.pattern().unwrap()[k] as usize;
                r.put_listed(i, hi[i] + lo[i]);
                hi[i] = 0.0;
                lo[i] = 0.0;
            }
        }
        debug_assert!(!r.is_dense());
        true
    }

    /// [`Self::refine`] on index-carrying vectors — same rounds, floor
    /// test, convergence shortcut and stats; norms are max-folds (order
    /// free) over the patterns, `x += d` runs over `d`'s pattern.
    fn refine_hv(&self, rhs: &HVec, x: &mut HVec, transpose: bool) {
        if rhs.is_dense() || x.is_dense() {
            self.refine(rhs.vals(), x.dense_mut(), transpose);
            return;
        }
        let m = self.n;
        let st = pfi_stats();
        if let Some(s) = st.as_ref() {
            if transpose {
                s.ref_bt.set(s.ref_bt.get() + 1);
            } else {
                s.ref_ft.set(s.ref_ft.get() + 1);
            }
        }
        if let Some(f) = refine_age_frac() {
            let young = match &self.lu {
                Engine::Sparse(lu) => self.eta_ticks.get() < f * lu.build_work() as f64,
                Engine::Dense(_) => false,
            };
            if young {
                if let Some(s) = st.as_ref() {
                    s.n_ref_skip.set(s.n_ref_skip.get() + 1);
                }
                return;
            }
        }
        let t0 = st.as_ref().map(|_| std::time::Instant::now());
        'outer: {
            let mut r = self.ws.r_hv.borrow_mut();
            r.ensure_len(m);
            let mut d = self.ws.d_hv.borrow_mut();
            d.ensure_len(m);
            let rhs_norm = rhs.max_abs();
            for round in 0..2 {
                let tr = st.as_ref().map(|_| std::time::Instant::now());
                let sparse_ok = {
                    let mut hi = self.ws.buf(&self.ws.sp_hi, m);
                    let mut lo = self.ws.buf(&self.ws.sp_lo, m);
                    self.residual_hv(rhs, x, &mut r, &mut hi, &mut lo, transpose)
                };
                if !sparse_ok {
                    let mut hi = self.ws.buf(&self.ws.hi_part, m);
                    let mut lo = self.ws.buf(&self.ws.lo_part, m);
                    let rv = r.dense_mut();
                    if transpose {
                        self.residual_bwd(rhs.vals(), x.vals(), rv, &mut hi, &mut lo);
                    } else {
                        self.residual_fwd(rhs.vals(), x.vals(), rv, &mut hi, &mut lo);
                    }
                } else if let Some(s) = st.as_ref() {
                    s.res_arm[0].set(s.res_arm[0].get() + 1);
                }
                if let (Some(s), Some(t)) = (st.as_ref(), tr) {
                    s.resid_ns
                        .set(s.resid_ns.get() + t.elapsed().as_nanos() as f64);
                }
                let gate0 = round == 0 && (refine_gate_on() || st.is_some());
                if round > 0 || gate0 {
                    let rnorm = r.max_abs();
                    let xnorm = x.max_abs();
                    let floor = (rhs_norm + xnorm * self.col_max * (m as f64)) * 1e-13;
                    if round == 0 {
                        if let Some(s) = st.as_ref() {
                            s.ref0_total.set(s.ref0_total.get() + 1);
                            if rnorm <= floor {
                                s.ref0_at_floor.set(s.ref0_at_floor.get() + 1);
                            }
                            s.ref0_ratio_sum.set(
                                s.ref0_ratio_sum.get() + (rnorm / floor.max(1e-300)).min(1e30),
                            );
                        }
                        if rnorm <= floor && refine_gate_on() {
                            if let Some(s) = st.as_ref() {
                                s.n_ref_skip.set(s.n_ref_skip.get() + 1);
                            }
                            break 'outer;
                        }
                    } else if rnorm <= floor {
                        break 'outer;
                    }
                }
                if transpose {
                    self.solve_bt_raw_hv(&r, &mut d);
                } else {
                    self.solve_b_raw_hv(&r, &mut d);
                }
                if let Some(s) = st.as_ref() {
                    s.n_corr.set(s.n_corr.get() + 1);
                    if round > 0 {
                        s.n_r1.set(s.n_r1.get() + 1);
                    }
                }
                let mut dn = 0.0f64;
                match d.pattern() {
                    Some(p) => {
                        let dv = d.vals();
                        for &i in p {
                            let di = dv[i as usize];
                            *x.entry(i as usize) += di;
                            dn = dn.max(di.abs());
                        }
                    }
                    None => {
                        let xv = x.dense_mut();
                        for (xi, &di) in xv.iter_mut().zip(d.vals()) {
                            *xi += di;
                            dn = dn.max(di.abs());
                        }
                    }
                }
                let xn = x.max_abs();
                if dn <= 1e-8 * xn.max(1.0) {
                    break 'outer;
                }
            }
        }
        if let (Some(s), Some(t)) = (st.as_ref(), t0) {
            s.refine_ns
                .set(s.refine_ns.get() + t.elapsed().as_nanos() as f64);
        }
    }

    /// [`Self::update`] with the entering column and `bbar` as
    /// index-carrying vectors: the eta and the stored basis column are
    /// built from the (sorted) patterns — bit-identical to `update`'s
    /// ascending, zero-dropped scans — and, unlike [`Self::update_sparse`],
    /// no extra synthetic-clock charge is made (the dual BTRANs every
    /// iteration), so the refactor cadence is `update`'s exactly.
    pub(crate) fn update_hv(
        &mut self,
        leaving_pos: usize,
        entering: &mut HVec,
        bbar: &mut HVec,
        bbar_sign: f64,
    ) -> Result<UpdateStatus, FactorError> {
        assert_eq!(entering.len(), self.n, "update: wrong column length");
        assert_eq!(bbar.len(), self.n, "update: wrong bbar length");
        assert!(leaving_pos < self.n, "update: position out of range");
        let pivot = bbar_sign * bbar[leaving_pos];
        if pivot.abs() < self.config.pivot_tol {
            rf_count_pivot_error();
            return Err(FactorError::RefactorizeNeeded { pivot });
        }
        let nz = bbar.nonzeros_sorted();
        let mut idx = Vec::with_capacity(nz.len());
        let mut val = Vec::with_capacity(nz.len());
        for (i, v) in nz {
            let v = bbar_sign * v;
            if i != leaving_pos && v != 0.0 {
                idx.push(i as i32);
                val.push(v);
            }
        }
        let entering_nz = entering.nonzeros_sorted();
        let cmax = entering_nz.iter().fold(0.0f64, |a, &(_, v)| a.max(v.abs()));
        self.swap_col(leaving_pos, entering_nz);
        self.push_eta(leaving_pos, pivot, idx, val, cmax)
    }

    /// Applies the standard simplex eta update: the basis column at position
    /// `leaving_pos` is replaced by `entering_col` (a column of the original
    /// constraint matrix, in the caller's space), recorded as
    /// `B_new = B · E`. The entering column's coordinates in the *current*
    /// basis (`b̄ = B⁻¹ a_q`) are supplied by the caller as `bbar` with
    /// `bbar_sign` — the simplex ratio test already computes exactly this
    /// solve (on the sign-applied column `sgn·a_q`, so pass `bbar_sign = sgn`);
    /// re-solving it here was a full refined FTRAN per pivot of pure
    /// duplicated work.
    ///
    /// The eta is built from `bbar_sign · bbar` — the caller's direction sign
    /// UN-applied, an exact ±1.0 multiply — so it is `B⁻¹ a_q` against the
    /// raw `entering_col`. IEEE round-to-nearest is sign-symmetric
    /// (`fl(−x) = −fl(x)` per operation), so a solve of `−a` is the bitwise
    /// negation of a solve of `a`: with `refine_ftran` on, this eta is
    /// BIT-IDENTICAL to the internal FTRAN `update` used to run (proved by
    /// `update_eta_matches_prerecomputed`). With it off the caller's b̄ is
    /// raw — the validated regime (eta from raw b̄).
    ///
    /// After appending, the file is checked against
    /// [`PfiConfig::eta_limit`] and refactorized if due. A pivot
    /// `|b̄[leaving_pos]|` below `pivot_tol` returns
    /// [`FactorError::RefactorizeNeeded`] (a stale/ill-conditioned eta file)
    /// rather than panicking; the basis state is left untouched.
    ///
    /// # Panics
    ///
    /// Panics if `entering_col.len()` or `bbar.len()` differ from `n`.
    pub fn update(
        &mut self,
        leaving_pos: usize,
        entering_col: &[f64],
        bbar: &[f64],
        bbar_sign: f64,
    ) -> Result<UpdateStatus, FactorError> {
        assert_eq!(entering_col.len(), self.n, "update: wrong column length");
        assert_eq!(bbar.len(), self.n, "update: wrong bbar length");
        assert!(leaving_pos < self.n, "update: position out of range");
        let pivot = bbar_sign * bbar[leaving_pos];
        if pivot.abs() < self.config.pivot_tol {
            rf_count_pivot_error();
            return Err(FactorError::RefactorizeNeeded { pivot });
        }
        let mut idx = Vec::new();
        let mut val = Vec::new();
        for (i, &v) in bbar.iter().enumerate() {
            let v = bbar_sign * v;
            if i != leaving_pos && v != 0.0 {
                idx.push(i as i32);
                val.push(v);
            }
        }
        let cmax = entering_col.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
        self.swap_col(leaving_pos, sparsify(entering_col));
        self.push_eta(leaving_pos, pivot, idx, val, cmax)
    }

    /// [`update`](Self::update) with the entering column given as its sorted
    /// nonzero list and the nonzero positions `bbar_nz` (ascending) of
    /// `bbar` already known: skips the O(n) scans and the dense column. The
    /// stored eta and basis column are bit-identical to `update`'s.
    pub fn update_sparse(
        &mut self,
        leaving_pos: usize,
        entering_nz: Vec<(usize, f64)>,
        bbar: &[f64],
        bbar_nz: &[usize],
        bbar_sign: f64,
    ) -> Result<UpdateStatus, FactorError> {
        assert_eq!(bbar.len(), self.n, "update: wrong bbar length");
        assert!(leaving_pos < self.n, "update: position out of range");

        // Eta column = b̄ of the RAW entering column = `bbar_sign · bbar`
        // (exact for ±1.0): the caller solved `bbar_sign · a_q` for the
        // pivot direction, the eta must describe `a_q` itself.
        let pivot = bbar_sign * bbar[leaving_pos];
        if pivot.abs() < self.config.pivot_tol {
            rf_count_pivot_error();
            return Err(FactorError::RefactorizeNeeded { pivot });
        }

        // Record eta: identity except column `leaving_pos` holding b̄; the
        // diagonal entry is carried as `pivot`, the rest kept sparse
        // (exact zeros dropped — they cost nothing and keep files tidy).
        // `bbar_sign · v` is exact and zero-preserving (±0.0 fails the
        // `!= 0.0` filter either way), so the stored vector is bit-identical
        // to the earlier one built from the internal solve.
        let mut idx = Vec::with_capacity(bbar_nz.len());
        let mut val = Vec::with_capacity(bbar_nz.len());
        for &i in bbar_nz {
            let v = bbar_sign * bbar[i];
            if i != leaving_pos && v != 0.0 {
                idx.push(i as i32);
                val.push(v);
            }
        }
        // The primal loop no longer BTRANs c_B every iteration (incremental
        // reduced costs), so its eta chains grew ~1.5x longer than the
        // refactor clock was calibrated for and lost accuracy (tuff). Charge
        // that skipped solve to the synthetic clock to keep the cadence.
        self.eta_ticks
            .set(self.eta_ticks.get() + ETA_TICK * (self.eta_nnz + self.etas.len()) as f64);
        let cmax = entering_nz.iter().fold(0.0f64, |a, &(_, v)| a.max(v.abs()));
        self.swap_col(leaving_pos, entering_nz);
        self.push_eta(leaving_pos, pivot, idx, val, cmax)
    }

    /// Swap `basis_nz[pos]` for `entering_nz` and keep the `basis_rows`
    /// mirror in step. The dense [`Self::update`] used to drop
    /// `rows_ok` instead of maintaining rows, permanently disabling the
    /// row-wise residual paths of [`Self::refine`] for every dense-update
    /// caller (the dual engine) — a sparse x then paid the O(nnz B)
    /// column sweep anyway (krka's β refines: 1-nonzero solve results over
    /// a 16.7k-nonzero basis). Both residual paths are bit-identical
    /// (pinned by `sparse_residual_bit_exact`), so enabling them cannot
    /// move a trajectory.
    fn swap_col(&mut self, pos: usize, entering_nz: Vec<(usize, f64)>) {
        if !self.rows_ok {
            self.basis_rows = rows_of(&self.basis_nz);
            self.rows_ok = true;
        }
        // Rows stay ascending by column (binary-searched remove /
        // insert instead of swap_remove + push) — the order the LU peel's
        // transpose would produce, so the mirror can stand in for it.
        // Every residual reading the mirror is order-independent per
        // accumulator (pinned by `sparse_residual_bit_exact`).
        for &(i, _) in &self.basis_nz[pos] {
            let row = &mut self.basis_rows[i];
            let k = row.partition_point(|&(j, _)| (j as usize) < pos);
            if k < row.len() && row[k].0 as usize == pos {
                row.remove(k);
            }
        }
        self.basis_nz[pos] = entering_nz;
        for &(i, a) in &self.basis_nz[pos] {
            let row = &mut self.basis_rows[i];
            let k = row.partition_point(|&(j, _)| (j as usize) < pos);
            row.insert(k, (pos as u32, a));
        }
    }

    /// Shared tail of the two update entry points: record the eta (the basis
    /// column is already replaced), widen `col_max`, refactorize if due.
    fn push_eta(
        &mut self,
        leaving_pos: usize,
        pivot: f64,
        idx: Vec<i32>,
        val: Vec<f64>,
        cmax: f64,
    ) -> Result<UpdateStatus, FactorError> {
        self.eta_nnz += idx.len() + 1;
        self.etas.push(Eta {
            pos: leaving_pos,
            pivot,
            idx,
            val,
        });
        self.col_max = self.col_max.max(cmax);

        if let Some(reason) = self.refactor_due() {
            self.refactorize_reason(reason)?;
            return Ok(UpdateStatus::Refactorized);
        }
        Ok(UpdateStatus::Updated)
    }

    /// Refactorization trigger: the hard cap [`PfiConfig::eta_limit`]
    /// (HiGHS `simplex_update_limit`, default 5000, lp_data/HighsOptions.h),
    /// or — sparse engine only — at least [`PfiConfig::min_updates`] updates
    /// AND the etas have cost more solve work than a rebuild costs
    /// (HiGHS HEkk.cpp `updateFactor`, `kSyntheticTickReinversionMinUpdateCount
    /// = 50` in util/HFactorConst.h; algorithm re-derived, not translated).
    /// Reports WHICH trigger fired for the cadence census.
    fn refactor_due(&self) -> Option<RfReason> {
        let n_eta = self.etas.len();
        if n_eta >= self.config.eta_limit {
            return Some(RfReason::Cap);
        }
        if let Engine::Sparse(lu) = &self.lu {
            let (k, m) = rfcad().unwrap_or((self.config.clock_scale, 0));
            if n_eta >= self.config.min_updates.max(m)
                && self.eta_ticks.get() >= lu.build_work() as f64 * k
            {
                return Some(RfReason::Clock);
            }
        }
        None
    }

    /// Rebuilds the LU factorization from the explicit current basis columns
    /// and empties the eta file. This is the numerical hygiene valve: after
    /// enough updates the eta passes accumulate round-off, and re-expanding
    /// resets the factorization to reference quality.
    pub fn refactorize(&mut self) -> Result<(), FactorError> {
        self.refactorize_reason(RfReason::Forced)
    }

    /// [`Self::refactorize`] recording WHY (cadence census).
    pub(crate) fn refactorize_reason(&mut self, reason: RfReason) -> Result<(), FactorError> {
        if let Some(st) = pfi_stats() {
            st.rf_reason[reason.index()].set(st.rf_reason[reason.index()].get() + 1);
            let depth = self.etas.len() as u64;
            st.rf_depth_sum.set(st.rf_depth_sum.get() + depth);
            st.rf_depth_max.set(st.rf_depth_max.get().max(depth));
            st.rf_nnz_sum.set(st.rf_nnz_sum.get() + self.eta_nnz as u64);
            let raw = st.n_raw.get();
            st.rf_solves_sum
                .set(st.rf_solves_sum.get() + raw - st.last_raw.get());
            st.last_raw.set(raw);
            let r1 = st.n_r1.get();
            st.rf_r1_sum.set(st.rf_r1_sum.get() + r1 - st.last_r1.get());
            st.last_r1.set(r1);
            // Sweep-cost census: estimated LIVE eta-file bytes at the
            // refactor instant (end-of-epoch watermark): 12 B/nnz for the
            // idx+val payloads + 96 B/eta for the two Vec allocations and
            // the Eta struct. The cont1-class simplex-side RSS question
            // (pre13-bucket-close) is answered by this vs the RSS probe.
            let bytes = self.eta_nnz as u64 * 12 + 96 * self.etas.len() as u64;
            st.rf_bytes_sum.set(st.rf_bytes_sum.get() + bytes);
            st.rf_bytes_max.set(st.rf_bytes_max.get().max(bytes));
            // RUSOLVE_RFDROP=1 drift probe: solve a fixed sparse probe rhs
            // through the CURRENT eta file, then through the fresh factors,
            // and measure how far the eta-chain answer had drifted. The
            // pre-solve charges eta ticks, but a due refactor zeroes the
            // clock immediately after, so cadence-relevant triggers are
            // unperturbed; forced event refactors see +1 solve in probe
            // runs only (measurement builds, never the default path).
            if rf_drift_probe_on() {
                let mut pre = vec![0.0; self.n];
                self.solve_bt_raw(&rf_probe_rhs(self.n), &mut pre);
                let r = self.refactorize_timed();
                if r.is_ok() {
                    let mut post = vec![0.0; self.n];
                    self.solve_bt_raw(&rf_probe_rhs(self.n), &mut post);
                    // Fresh factors charge ZERO eta ticks (empty file), so
                    // the post-solve cannot shift the next epoch's clock.
                    let scale = post.iter().fold(0.0f64, |a, &v| a.max(v.abs())).max(1e-30);
                    let d = pre
                        .iter()
                        .zip(&post)
                        .fold(0.0f64, |a, (p, q)| a.max((p - q).abs()))
                        / scale;
                    st.rf_drift_sum.set(st.rf_drift_sum.get() + d);
                    st.rf_drift_max.set(st.rf_drift_max.get().max(d));
                    st.rf_drift_n.set(st.rf_drift_n.get() + 1);
                }
                return r;
            }
        }
        self.refactorize_timed()
    }

    fn refactorize_timed(&mut self) -> Result<(), FactorError> {
        let t = std::time::Instant::now();
        let r = self.do_refactorize();
        self.refactor_nanos
            .set(self.refactor_nanos.get() + t.elapsed().as_nanos() as f64);
        r
    }

    fn do_refactorize(&mut self) -> Result<(), FactorError> {
        let rows = self.rows_ok.then_some(&self.basis_rows[..]);
        let lu = factor_with(
            &self.basis_nz,
            rows,
            &self.config,
            &mut self.engine_choice,
            &mut self.fws,
        )?;
        // The superseded LU's buffers seed the next build.
        if let Engine::Sparse(old) = std::mem::replace(&mut self.lu, lu) {
            self.fws.recycle(old);
        }
        self.epoch += 1;
        self.etas.clear();
        self.eta_nnz = 0;
        self.eta_ticks.set(0.0);
        Ok(())
    }

    /// Cumulative nanoseconds spent inside [`Self::refactorize`] (profiling;
    /// read by the dual engine's RUSOLVE_DUAL_STATS block).
    pub fn refactor_nanos(&self) -> f64 {
        self.refactor_nanos.get()
    }

    /// Profiling: LU factor size (L nnz + U nnz + 3m) of the current
    /// engine, or 0 for the dense engine.
    pub(crate) fn lu_solve_work(&self) -> usize {
        match &self.lu {
            Engine::Sparse(lu) => lu.solve_work(),
            Engine::Dense(_) => 0,
        }
    }

    /// Profiling: total nonzeros of the explicit basis columns.
    pub(crate) fn basis_nonzeros(&self) -> usize {
        self.basis_nz.iter().map(|c| c.len()).sum()
    }

    /// Estimated LIVE eta-file bytes right now (see `rf_bytes_sum` in
    /// `PfiStats` for the estimator's rationale). Read by the dual's
    /// `RUSOLVE_STATS_EVERY` periodic dump for cap-hitting giant-LP
    /// census runs (cont1-class) that never reach the Optimal-exit stats.
    pub(crate) fn eta_live_bytes(&self) -> usize {
        self.eta_nnz * 12 + 96 * self.etas.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::XorShift64;
    use crate::test_support::{max_abs_diff, naive_solve, random_well_conditioned};

    fn config() -> PfiConfig {
        PfiConfig {
            pivot_tol: 1e-10,
            eta_limit: 8,
            min_updates: 50,
            clock_scale: 1.0,
            engine: PfiEngine::Dense,
        }
    }

    fn rhs(n: usize, seed: u64) -> Vec<f64> {
        let mut rng = XorShift64::new(seed);
        (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect()
    }

    #[test]
    fn pfi_build_solve_matches_naive_reference() {
        for seed in [11, 22, 33] {
            let a = random_well_conditioned(30, seed);
            let pfi = ProductFormFactorization::build(&a, config()).expect("buildable");
            let b = rhs(30, seed);
            let mut x = vec![0.0; 30];
            pfi.solve_b(&b, &mut x);
            assert!(
                max_abs_diff(&x, &naive_solve(&a, &b)) < 1e-9,
                "seed {seed}: PFI solve_b diverged"
            );
            assert_eq!(pfi.eta_count(), 0);
        }
    }

    #[test]
    fn pfi_build_solve_bt_matches_naive_reference() {
        for seed in [44, 55] {
            let a = random_well_conditioned(28, seed);
            let pfi = ProductFormFactorization::build(&a, config()).expect("buildable");
            let b = rhs(28, seed);
            let mut x = vec![0.0; 28];
            pfi.solve_bt(&b, &mut x);

            let n = a.n();
            let mut at = DenseMatrix::zeros(n);
            for j in 0..n {
                for i in 0..n {
                    at.set(i, j, a.get(j, i));
                }
            }
            assert!(
                max_abs_diff(&x, &naive_solve(&at, &b)) < 1e-9,
                "seed {seed}: PFI solve_bt diverged"
            );
        }
    }

    #[test]
    fn solve_bt_satisfies_inner_product_identity() {
        // <B⁻ᵀ b, v> == <b, B⁻¹ v> for random v — a structural check of the
        // transpose solve independent of any reference implementation.
        let a = random_well_conditioned(25, 77);
        let pfi = ProductFormFactorization::build(&a, config()).expect("buildable");
        let b = rhs(25, 78);
        let v = rhs(25, 79);

        let mut y = vec![0.0; 25]; // B⁻ᵀ b
        pfi.solve_bt(&b, &mut y);
        let mut z = vec![0.0; 25]; // B⁻¹ v
        pfi.solve_b(&v, &mut z);

        let lhs: f64 = y.iter().zip(&v).map(|(p, q)| p * q).sum();
        let rhs_val: f64 = b.iter().zip(&z).map(|(p, q)| p * q).sum();
        assert!(
            (lhs - rhs_val).abs() <= 1e-8 * lhs.abs().max(rhs_val.abs()),
            "inner product identity broken: {lhs} vs {rhs_val}"
        );
    }

    /// The cadence lever parses `K[:M]`, clamps K to >= 1 (a
    /// sub-1 multiplier would refactor MORE often than calibrated), and
    /// falls back to the neutral default on garbage (unset ⇒ (1.0, 0) ⇒
    /// bit-identical default cadence).
    #[test]
    fn rfcad_lever_parses() {
        assert_eq!(parse_rfcad("2"), (2.0, 0));
        assert_eq!(parse_rfcad("2:50"), (2.0, 50));
        assert_eq!(parse_rfcad("1:50"), (1.0, 50));
        assert_eq!(parse_rfcad(" 1.5 : 8 "), (1.5, 8));
        assert_eq!(parse_rfcad("0.5"), (1.0, 0));
        assert_eq!(parse_rfcad("junk"), (1.0, 0));
        assert_eq!(parse_rfcad("3:junk"), (3.0, 0));
    }

    /// Sweep-cost lever: `RUSOLVE_BT_GATE` parses, and the
    /// default ratio (4.0) is DECISION-EQUIVALENT to the old integer gate
    /// `a > 4·n` at every realistic size (f64 integers are exact ≤ 2^53),
    /// so unset is behavior-identical. Ratio 1.0 opens the pattern walk
    /// on the istanbul-class file (eta_nnz 72k at n 20.3k, ≤ 9-nonzero y
    /// plain-walking 72k entries per solve) while keeping the krka-class
    /// thin file (network LP, sparse etas) off.
    #[test]
    fn bt_gate_lever_parses_and_default_matches_old_gate() {
        assert_eq!(parse_bt_gate("1"), 1.0);
        assert_eq!(parse_bt_gate(" 0.5 "), 0.5);
        assert_eq!(parse_bt_gate("junk"), 4.0);
        assert_eq!(parse_bt_gate("-1"), 4.0);
        assert_eq!(parse_bt_gate("1e9"), 4.0);
        assert!(
            !bt_gate_open(1_000_000, 0, 15, 4.0),
            "pat_max=0 ⇒ never track"
        );
        for n in [17usize, 100, 1000, 20346, 44023, 1_000_000] {
            for a in [
                0usize,
                1,
                4 * n,
                4 * n + 1,
                4 * n.saturating_sub(1),
                9 * n + 7,
            ] {
                assert_eq!(
                    bt_gate_open(a, 0, n, 4.0),
                    n / 16 > 0 && a > 4 * n,
                    "n={n} a={a}"
                );
            }
        }
        assert!(!bt_gate_open(16, 0, 17, 4.0) && !bt_gate_open(16, 0, 16, 4.0));
        // istanbul-class file: old gate closed (plain 72k-entry walks),
        // R=1 opens the pattern walk.
        assert!(!bt_gate_open(72_000, 8, 20_346, 4.0));
        assert!(bt_gate_open(72_000, 8, 20_346, 1.0));
        // krka-class thin file: both arms keep tracking off.
        assert!(!bt_gate_open(500, 40, 12_500, 1.0));
    }

    /// Candidate (c) lever parses a fraction in (0, 1] and refuses
    /// everything else (unset ⇒ None ⇒ refine always runs, the
    /// calibrated path).
    #[test]
    fn refine_age_lever_parses() {
        assert_eq!(parse_refine_age("0.05"), Some(0.05));
        assert_eq!(parse_refine_age(" 1 "), Some(1.0));
        assert_eq!(parse_refine_age("0"), None);
        assert_eq!(parse_refine_age("-0.1"), None);
        assert_eq!(parse_refine_age("1.5"), None);
        assert_eq!(parse_refine_age("junk"), None);
    }

    /// Candidate (c′) lever: only "1" turns the round-0 floor gate on.
    #[test]
    fn refine_gate_lever_parses() {
        assert!(parse_refine_gate("1"));
        assert!(parse_refine_gate(" 1"));
        assert!(!parse_refine_gate("0"));
        assert!(!parse_refine_gate("junk"));
        assert!(!parse_refine_gate("2"));
    }

    /// Wilkinson refinement recovers forward accuracy on an ill-conditioned
    /// basis (col₁ ∥ col₀ at distance 1e-6 ⇒ κ·eps ≪ 1, the convergent
    /// regime): raw error is reduced by orders of magnitude, never regresses.
    /// (scsd6's phase-2 2-cycle — wrong reduced costs from drifted etas — is
    /// the real-world failure this fixes; verified on-bench Optimal.)
    #[test]
    fn refine_recovers_ill_conditioned_solve() {
        let n = 12;
        let mut a = random_well_conditioned(n, 5);
        for i in 0..n {
            let c0 = a.get(i, 0);
            let pert = if i == 0 { 1e-6 } else { 0.0 };
            a.set(i, 1, c0 * (1.0 + 1e-6) + pert);
        }
        let mut rng = XorShift64::new(17);
        let x0: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let mut b = vec![0.0; n];
        for (j, &x0j) in x0.iter().enumerate() {
            for (i, bi) in b.iter_mut().enumerate() {
                *bi += a.get(i, j) * x0j;
            }
        }
        let pfi = ProductFormFactorization::build(&a, config()).expect("buildable");
        let rel = |x: &[f64]| {
            x.iter()
                .zip(&x0)
                .map(|(g, e)| (g - e).abs() / e.abs().max(1.0))
                .fold(0.0f64, f64::max)
        };
        let (mut raw, mut refd) = (vec![0.0; n], vec![0.0; n]);
        pfi.solve_b_raw(&b, &mut raw);
        pfi.solve_b(&b, &mut refd);
        let (e_raw, e_ref) = (rel(&raw), rel(&refd));
        // The refined error floors at ~9e-10 here — the representability
        // limit of the fp-constructed rhs, not the solver (measured: further
        // rounds are inert). The assertion that matters: DW-residual
        // refinement contracts, plain residuals do not (probe-measured).
        assert!(
            e_ref <= e_raw * 0.8 + 1e-14 && e_ref < 1e-8,
            "forward: raw={e_raw:.2e} refined={e_ref:.2e}"
        );

        // Transpose direction: residual ‖Bᵀx − c‖∞ must not regress.
        let c: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let resid = |x: &[f64]| {
            (0..n)
                .map(|j| (0..n).map(|i| a.get(i, j) * x[i]).sum::<f64>() - c[j])
                .map(|v| v.abs())
                .fold(0.0f64, f64::max)
        };
        let (mut praw, mut pref) = (vec![0.0; n], vec![0.0; n]);
        pfi.solve_bt_raw(&c, &mut praw);
        pfi.solve_bt(&c, &mut pref);
        let (r_raw, r_ref) = (resid(&praw), resid(&pref));
        assert!(
            r_ref <= r_raw * 1.5 + 1e-12,
            "transpose: raw={r_raw:.2e} refined={r_ref:.2e}"
        );
    }

    #[test]
    fn pfi_build_rejects_singular_basis() {
        let mut a = random_well_conditioned(10, 9);
        for i in 0..10 {
            a.set(i, 4, a.get(i, 3) * 2.0);
        }
        assert!(matches!(
            ProductFormFactorization::build(&a, config()),
            Err(FactorError::Singular { column: 4, .. })
        ));
    }

    /// A column that is (numerically) in the span of the basis with a zero
    /// component at `pos`: guarantees `b̄[pos] == 0` for the pivot guard test.
    fn column_with_zero_at(basis_cols: &[Vec<f64>], pos: usize, seed: u64) -> Vec<f64> {
        let n = basis_cols.len();
        let mut rng = XorShift64::new(seed);
        let coeffs: Vec<f64> = (0..n)
            .map(|i| {
                if i == pos {
                    0.0
                } else {
                    rng.uniform(-1.0, 1.0)
                }
            })
            .collect();
        let mut a = vec![0.0; n];
        for (j, &cj) in coeffs.iter().enumerate() {
            if cj != 0.0 {
                for i in 0..n {
                    a[i] += cj * basis_cols[j][i];
                }
            }
        }
        a
    }

    /// Drive `steps` max-pivot updates, keeping an explicit copy of the
    /// basis so every solve can be checked against the naive reference.
    /// Returns the max |PFI solve − naive solve| residual over random rhs.
    fn run_update_sequence(initial: &DenseMatrix, steps: usize, cfg: PfiConfig, seed: u64) -> f64 {
        let n = initial.n();
        let mut pfi = ProductFormFactorization::build(initial, cfg).expect("buildable");
        let mut basis: Vec<Vec<f64>> = (0..n).map(|j| initial.col(j).to_vec()).collect();
        let mut rng = XorShift64::new(seed);
        let mut worst = 0.0_f64;

        for step in 0..steps {
            // Random entering column; choose the max-|b̄| position as the
            // leaving slot (max-pivot rule keeps updates well conditioned).
            let entering: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let mut bbar = vec![0.0; n];
            pfi.solve_b(&entering, &mut bbar);
            let mut pos = 0;
            for i in 1..n {
                if bbar[i].abs() > bbar[pos].abs() {
                    pos = i;
                }
            }
            let status = pfi
                .update(pos, &entering, &bbar, 1.0)
                .expect("update accepted");
            basis[pos] = entering.clone();
            if step == 0 {
                assert_eq!(status, UpdateStatus::Updated, "eta_limit too small");
            }

            let b = rhs(n, seed.wrapping_add(step as u64 * 31));
            let mut x = vec![0.0; n];
            pfi.solve_b(&b, &mut x);
            let ref_x = naive_solve(&DenseMatrix::from_columns(&basis), &b);
            worst = worst.max(max_abs_diff(&x, &ref_x));

            let mut xt = vec![0.0; n];
            pfi.solve_bt(&b, &mut xt);
            let mut bt = DenseMatrix::zeros(n);
            // B^T[i][j] = B[j][i] = basis[i][j]
            for (i, col) in basis.iter().enumerate() {
                for (j, &v) in col.iter().enumerate() {
                    bt.set(i, j, v);
                }
            }
            let ref_xt = naive_solve(&bt, &b);
            worst = worst.max(max_abs_diff(&xt, &ref_xt));
        }
        worst
    }

    #[test]
    fn single_update_matches_naive_on_new_basis() {
        let a = random_well_conditioned(20, 101);
        assert!(run_update_sequence(&a, 1, config(), 102) < 1e-9);
    }

    #[test]
    fn long_eta_file_updates_match_naive() {
        // 17 updates with eta_limit high enough that nothing auto-refactorizes:
        // exercises the eta forward/back passes over a long file.
        let a = random_well_conditioned(20, 111);
        let no_refactor = PfiConfig {
            pivot_tol: 1e-10,
            eta_limit: 1000,
            min_updates: 50,
            clock_scale: 1.0,
            engine: PfiEngine::Dense,
        };
        let worst = run_update_sequence(&a, 17, no_refactor, 112);
        assert!(worst < 1e-7, "eta-file drift too large: {worst:e}");
    }

    #[test]
    fn updates_with_auto_refactorize_match_naive() {
        // eta_limit=8 inside `config()` forces rebuilds mid-sequence; results
        // must be reference-grade again after each refactorization.
        let a = random_well_conditioned(16, 121);
        let worst = run_update_sequence(&a, 20, config(), 122);
        assert!(worst < 1e-8, "post-refactorize solve diverged: {worst:e}");
    }

    #[test]
    fn eta_update_then_solve_equals_refactorized_solve() {
        let a = random_well_conditioned(24, 131);
        let mut pfi = ProductFormFactorization::build(&a, config()).expect("buildable");
        let mut rng = XorShift64::new(132);
        for _ in 0..5 {
            let entering: Vec<f64> = (0..24).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let mut bbar = vec![0.0; 24];
            pfi.solve_b(&entering, &mut bbar);
            let pos = (0..24)
                .max_by(|i, j| bbar[*i].abs().total_cmp(&bbar[*j].abs()))
                .unwrap();
            pfi.update(pos, &entering, &bbar, 1.0).unwrap();
        }
        assert_eq!(pfi.eta_count(), 5);

        let b = rhs(24, 133);
        let mut x_eta = vec![0.0; 24];
        pfi.solve_b(&b, &mut x_eta);
        let mut x_bt_eta = vec![0.0; 24];
        pfi.solve_bt(&b, &mut x_bt_eta);

        pfi.refactorize().unwrap();
        assert_eq!(pfi.eta_count(), 0);
        let mut x_ref = vec![0.0; 24];
        pfi.solve_b(&b, &mut x_ref);
        let mut x_bt_ref = vec![0.0; 24];
        pfi.solve_bt(&b, &mut x_bt_ref);

        assert!(max_abs_diff(&x_eta, &x_ref) < 1e-9);
        assert!(max_abs_diff(&x_bt_eta, &x_bt_ref) < 1e-9);
    }

    #[test]
    fn refactorize_threshold_triggers_and_resets_file() {
        let a = random_well_conditioned(12, 141);
        let cfg = PfiConfig {
            pivot_tol: 1e-10,
            eta_limit: 5,
            min_updates: 50,
            clock_scale: 1.0,
            engine: PfiEngine::Dense,
        };
        let mut pfi = ProductFormFactorization::build(&a, cfg).expect("buildable");
        let mut rng = XorShift64::new(142);
        let mut statuses = Vec::new();
        for _ in 0..5 {
            let entering: Vec<f64> = (0..12).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let mut bbar = vec![0.0; 12];
            pfi.solve_b(&entering, &mut bbar);
            let pos = (0..12)
                .max_by(|i, j| bbar[*i].abs().total_cmp(&bbar[*j].abs()))
                .unwrap();
            statuses.push(pfi.update(pos, &entering, &bbar, 1.0).unwrap());
        }
        assert_eq!(
            statuses,
            vec![
                UpdateStatus::Updated,
                UpdateStatus::Updated,
                UpdateStatus::Updated,
                UpdateStatus::Updated,
                UpdateStatus::Refactorized
            ]
        );
        assert_eq!(pfi.eta_count(), 0, "file must be empty after threshold");
    }

    #[test]
    fn tiny_update_pivot_returns_refactorize_needed() {
        let a = random_well_conditioned(15, 151);
        let mut pfi = ProductFormFactorization::build(&a, config()).expect("buildable");
        let cols: Vec<Vec<f64>> = (0..15).map(|j| a.col(j).to_vec()).collect();
        let entering = column_with_zero_at(&cols, 7, 152);
        let mut bbar = vec![0.0; 15];
        pfi.solve_b(&entering, &mut bbar);
        let err = pfi
            .update(7, &entering, &bbar, 1.0)
            .expect_err("zero pivot must not be accepted");
        assert!(
            matches!(err, FactorError::RefactorizeNeeded { .. }),
            "{err:?}"
        );
        assert_eq!(pfi.eta_count(), 0, "failed update must leave state clean");
    }

    #[test]
    fn sparse_engine_matches_dense() {
        // Engine parity: same basis, same update sequence, same eta cadence —
        // solves must track the dense engine to refinement accuracy. The
        // engine_choice assert keeps the test honest (a silent sticky
        // fallback to Dense would make it vacuous).
        let a = random_well_conditioned(30, 77);
        let mut dense = ProductFormFactorization::build(&a, config()).unwrap();
        let mut sparse = ProductFormFactorization::build(
            &a,
            PfiConfig {
                engine: PfiEngine::Sparse,
                ..config()
            },
        )
        .unwrap();
        let mut rng = XorShift64::new(999);
        for round in 0..12 {
            let b: Vec<f64> = (0..30).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let (mut xd, mut xs) = (vec![0.0; 30], vec![0.0; 30]);
            dense.solve_b(&b, &mut xd);
            sparse.solve_b(&b, &mut xs);
            assert!(max_abs_diff(&xd, &xs) < 1e-7, "round {round} fwd");
            let (mut yd, mut ys) = (vec![0.0; 30], vec![0.0; 30]);
            dense.solve_bt(&b, &mut yd);
            sparse.solve_bt(&b, &mut ys);
            assert!(max_abs_diff(&yd, &ys) < 1e-7, "round {round} bwd");
            let col: Vec<f64> = (0..30).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let pos = (round * 7) % 30;
            // Each engine builds the eta from ITS OWN solve — the two
            // engines' b̄ agree only to refinement accuracy, not bitwise.
            let (mut bd, mut bs) = (vec![0.0; 30], vec![0.0; 30]);
            dense.solve_b(&col, &mut bd);
            sparse.solve_b(&col, &mut bs);
            dense.update(pos, &col, &bd, 1.0).unwrap();
            sparse.update(pos, &col, &bs, 1.0).unwrap();
        }
        assert_eq!(
            sparse.engine_choice,
            PfiEngine::Sparse,
            "sparse engine silently fell back to dense"
        );
    }

    /// Verbatim copy of the earlier DENSE forward residual loop (iterate all
    /// m entries per column with the `a == 0.0` / `xj == 0.0` skips, same
    /// double-word arithmetic). `sparse_residual_bit_exact` pins the sparse
    /// production loops ([`ProductFormFactorization::residual_fwd`]) to this.
    fn dense_residual_fwd(cols: &[Vec<f64>], rhs: &[f64], x: &[f64]) -> Vec<f64> {
        let mut hi_part = rhs.to_vec();
        let mut lo_part = vec![0.0; cols.len()];
        for (j, col) in cols.iter().enumerate() {
            let xj = x[j];
            if xj == 0.0 {
                continue;
            }
            for ((si, ci), &a) in hi_part.iter_mut().zip(lo_part.iter_mut()).zip(col.iter()) {
                if a == 0.0 {
                    continue;
                }
                let p = a * xj;
                let pe = a.mul_add(-xj, p);
                let (hi, lo) = two_sum(*si, -p);
                *si = hi;
                *ci += lo + pe;
            }
        }
        hi_part
            .iter()
            .zip(&lo_part)
            .map(|(hi, lo)| hi + lo)
            .collect()
    }

    /// Verbatim copy of the earlier DENSE transpose residual loop; companion
    /// of [`dense_residual_fwd`].
    fn dense_residual_bwd(cols: &[Vec<f64>], rhs: &[f64], x: &[f64]) -> Vec<f64> {
        let mut r = vec![0.0; cols.len()];
        for (j, col) in cols.iter().enumerate() {
            let mut s = 0.0;
            let mut c = 0.0;
            for (&a, &xi) in col.iter().zip(x.iter()) {
                if a == 0.0 || xi == 0.0 {
                    continue;
                }
                let p = a * xi;
                let pe = a.mul_add(-xi, p);
                let (hi, lo) = two_sum(s, p);
                s = hi;
                c += lo - pe;
            }
            let (hi, lo) = two_sum(rhs[j], -s);
            r[j] = hi + (lo - c);
        }
        r
    }

    /// The sparse residual sweeps must reproduce the old dense loops
    /// BIT-FOR-BIT (== on f64, not approx) on a genuinely sparse basis,
    /// through the full lifecycle (factor, updates, auto-refactorizes,
    /// refined forward/transpose solves).
    #[test]
    fn sparse_residual_bit_exact() {
        let m = 120;
        let mut rng = XorShift64::new(4242);
        // Random sparse basis (~3% density) with a column diagonal bump to
        // strict diagonal dominance ⇒ nonsingular (Levy–Dobbins–Gantsevich).
        let mut a = DenseMatrix::zeros(m);
        for j in 0..m {
            let mut offdiag_abs_sum = 0.0;
            for i in 0..m {
                if i != j && rng.next_f64() < 0.03 {
                    let v = rng.uniform(-10.0, 10.0);
                    a.set(i, j, v);
                    offdiag_abs_sum += v.abs();
                }
            }
            a.set(j, j, offdiag_abs_sum + 1.0);
        }
        let nnz: usize = (0..m)
            .map(|j| a.col(j).iter().filter(|&&v| v != 0.0).count())
            .sum();
        assert!(nnz * 2 < m * m / 5, "basis not sparse: {nnz} nnz in {m}²");

        let cfg = PfiConfig {
            pivot_tol: 1e-10,
            eta_limit: 5, // forces auto-refactorizes mid-sequence
            min_updates: 50,
            clock_scale: 1.0,
            engine: PfiEngine::Dense,
        };
        let mut pfi = ProductFormFactorization::build(&a, cfg).expect("buildable");

        let check = |pfi: &ProductFormFactorization, rhs: &[f64], x: &[f64], tag: &str| {
            let cols = dense_cols(&pfi.basis_nz, m);
            let mut r = vec![0.0; m];
            let (mut hi, mut lo) = (vec![0.0; m], vec![0.0; m]);
            pfi.residual_fwd(rhs, x, &mut r, &mut hi, &mut lo);
            assert_eq!(r, dense_residual_fwd(&cols, rhs, x), "fwd {tag}");
            pfi.residual_bwd(rhs, x, &mut r, &mut hi, &mut lo);
            assert_eq!(r, dense_residual_bwd(&cols, rhs, x), "bwd {tag}");
        };

        for round in 0..8 {
            let b = rhs(m, 900 + round);
            let mut x = vec![0.0; m];
            pfi.solve_b_raw(&b, &mut x);
            check(&pfi, &b, &x, "raw solve");
            // Structurally sparse x exercises the xi == 0.0 / xj == 0.0 skips.
            let xs: Vec<f64> = x
                .iter()
                .enumerate()
                .map(|(i, &v)| if i % 2 == 0 { 0.0 } else { v })
                .collect();
            check(&pfi, &b, &xs, "sparse x");
            // Hyper-sparse x: 10%-dense, fires the ROW-WISE residual
            // arms that dense-`update`-maintained `basis_rows` now keep
            // available (pre-slice-2 the dual engine's updates left
            // `rows_ok == false` and this path never ran).
            let xh: Vec<f64> = x
                .iter()
                .enumerate()
                .map(|(i, &v)| if i % 10 == 0 { v } else { 0.0 })
                .collect();
            check(&pfi, &b, &xh, "hyper-sparse x");

            // Lifecycle sanity: refined solves still track the naive
            // reference on the explicit current basis (approx by design;
            // bit-exactness is the residual assert above).
            let mut xb = vec![0.0; m];
            pfi.solve_b(&b, &mut xb);
            let cols = dense_cols(&pfi.basis_nz, m);
            let ref_x = naive_solve(&DenseMatrix::from_columns(&cols), &b);
            assert!(
                max_abs_diff(&xb, &ref_x) < 1e-6,
                "round {round} fwd drifted: {:.2e}",
                max_abs_diff(&xb, &ref_x)
            );
            let mut bt = DenseMatrix::zeros(m);
            for (i, col) in cols.iter().enumerate() {
                for (j, &v) in col.iter().enumerate() {
                    bt.set(i, j, v);
                }
            }
            let mut xbt = vec![0.0; m];
            pfi.solve_bt(&b, &mut xbt);
            let ref_xt = naive_solve(&bt, &b);
            assert!(
                max_abs_diff(&xbt, &ref_xt) < 1e-6,
                "round {round} bwd drifted"
            );

            // Pivot with a SPARSE entering column so basis_nz stays sparse
            // after replacement (max-|b̄| leaving slot keeps it well posed).
            let entering: Vec<f64> = (0..m)
                .map(|_| {
                    if rng.next_f64() < 0.25 {
                        rng.uniform(-5.0, 5.0)
                    } else {
                        0.0
                    }
                })
                .collect();
            let mut bbar = vec![0.0; m];
            pfi.solve_b(&entering, &mut bbar);
            let pos = (0..m)
                .max_by(|i, j| bbar[*i].abs().total_cmp(&bbar[*j].abs()))
                .unwrap();
            pfi.update(pos, &entering, &bbar, 1.0)
                .expect("update accepted");
        }
        let b = rhs(m, 999);
        let mut x = vec![0.0; m];
        pfi.solve_b(&b, &mut x);
        check(&pfi, &b, &x, "final state");
    }

    /// The eta built from the CALLER's pre-computed b̄ is bit-identical
    /// to the eta the deleted internal FTRAN used to build — including the
    /// `from_ub` sign flip. Two PFIs are driven through the same pivot
    /// sequence on a random sparse LP-like basis:
    ///
    /// * **A** reproduces the OLD path: an internal-style refined solve of the
    ///   RAW (unsigned) entering column — exactly what `update` solved before
    ///   the change — fed back with `bbar_sign = 1.0`.
    /// * **B** is the NEW path: the caller solves the SIGN-applied column
    ///   (`sgn·a_q`, the simplex ratio-test b̄) and passes `bbar` + `sgn`.
    ///
    /// After every pivot the newly appended eta vectors must match BITWISE
    /// (pivot bits, index lists, value bits) and subsequent refined forward/
    /// transpose solves must be bit-identical, so the whole pivot path
    /// downstream is unchanged. This is the debug-proof the sign-symmetry
    /// argument in [`ProductFormFactorization::update`] rests on.
    #[test]
    fn update_eta_matches_prerecomputed() {
        let m = 80;
        let mut rng = XorShift64::new(0x611b_0b1a);
        // Sparse diagonally dominant basis (same recipe as
        // sparse_residual_bit_exact): LP-like structure, refined solves on.
        let mut a = DenseMatrix::zeros(m);
        for j in 0..m {
            let mut offdiag_abs_sum = 0.0;
            for i in 0..m {
                if i != j && rng.next_f64() < 0.04 {
                    let v = rng.uniform(-10.0, 10.0);
                    a.set(i, j, v);
                    offdiag_abs_sum += v.abs();
                }
            }
            a.set(j, j, offdiag_abs_sum + 1.0);
        }
        let cfg = PfiConfig {
            pivot_tol: 1e-10,
            eta_limit: 1000, // no auto-refactorize: every eta stays inspectable
            min_updates: usize::MAX,
            clock_scale: 1.0,
            engine: PfiEngine::Sparse,
        };
        let mut old = ProductFormFactorization::build(&a, cfg).expect("buildable A");
        let mut new = ProductFormFactorization::build(&a, cfg).expect("buildable B");

        for step in 0..25 {
            // Sparse LP-like entering column; alternate the direction sign to
            // cover both the from-lb (sgn=+1) and from-ub (sgn=−1) call shapes.
            let entering: Vec<f64> = (0..m)
                .map(|_| {
                    if rng.next_f64() < 0.25 {
                        rng.uniform(-5.0, 5.0)
                    } else {
                        0.0
                    }
                })
                .collect();
            let sgn = if step % 2 == 0 { 1.0 } else { -1.0 };

            // A: OLD path — internal solve of the RAW column, unsigned eta.
            let mut bbar_internal = vec![0.0; m];
            old.solve_b(&entering, &mut bbar_internal);
            let pos = (0..m)
                .max_by(|i, j| bbar_internal[*i].abs().total_cmp(&bbar_internal[*j].abs()))
                .unwrap();
            old.update(pos, &entering, &bbar_internal, 1.0)
                .expect("old-path update");

            // B: NEW path — caller solves the SIGN-applied column (exactly
            // what the simplex ratio test hands to `update`) plus the sign.
            let aq_s: Vec<f64> = entering.iter().map(|&v| sgn * v).collect();
            let mut bbar = vec![0.0; m];
            new.solve_b(&aq_s, &mut bbar);
            new.update(pos, &entering, &bbar, sgn)
                .expect("new-path update");

            // Bitwise eta comparison (the freshly appended eta).
            let (ea, eb) = (&old.etas[step], &new.etas[step]);
            assert_eq!(ea.pos, eb.pos, "step {step}: eta pos");
            assert_eq!(ea.pivot.to_bits(), eb.pivot.to_bits(), "step {step}: pivot");
            assert_eq!(ea.idx, eb.idx, "step {step}: eta index list");
            assert_eq!(ea.val.len(), eb.val.len(), "step {step}: eta nnz");
            for (k, (x, y)) in ea.val.iter().zip(&eb.val).enumerate() {
                assert_eq!(x.to_bits(), y.to_bits(), "step {step}: eta val {k}");
            }

            // Subsequent solves bit-identical on both factorizations.
            let b = rhs(m, 17_000 + step as u64);
            let (mut x1, mut x2) = (vec![0.0; m], vec![0.0; m]);
            old.solve_b(&b, &mut x1);
            new.solve_b(&b, &mut x2);
            assert_eq!(x1, x2, "step {step}: forward solve diverged");
            let (mut t1, mut t2) = (vec![0.0; m], vec![0.0; m]);
            old.solve_bt(&b, &mut t1);
            new.solve_bt(&b, &mut t2);
            assert_eq!(t1, t2, "step {step}: transpose solve diverged");
        }
    }

    /// The pattern-aware BTRAN eta dot (sparse-rhs walk) must equal
    /// the plain newest-first loop BIT-FOR-BIT across the full lifecycle —
    /// unit-rhs (β) BTRANs with a long file of mixed-density etas, plus
    /// dense-rhs BTRANs that must fall back to the plain loop. Compares
    /// against a reference implementation of the OLD pass.
    #[test]
    fn bt_eta_pattern_walk_bit_exact() {
        let m = 300;
        let mut rng = XorShift64::new(2024);
        // Sparse diagonally dominant basis (nonsingular by Levy–Dobbins).
        let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        for j in 0..m {
            let mut off = 0.0f64;
            for i in 0..m {
                if i != j && rng.next_f64() < 0.05 {
                    let v = rng.uniform(-5.0, 5.0);
                    cols[j].push((i, v));
                    off += v.abs();
                }
            }
            cols[j].push((j, off + 1.0));
        }
        let cfg = PfiConfig {
            pivot_tol: 1e-10,
            eta_limit: 10_000, // never auto-refactor: exercise the file
            min_updates: 50,
            clock_scale: 1.0,
            engine: PfiEngine::Sparse,
        };
        let mut pfi = ProductFormFactorization::build_sparse(cols, cfg).unwrap();
        // A long eta file of mixed-density etas (some dense b̄ columns).
        for step in 0..80 {
            let mut entering: Vec<(usize, f64)> = Vec::new();
            if step % 3 == 0 {
                for i in 0..m {
                    entering.push((i, rng.uniform(-1.0, 1.0))); // dense
                }
            } else {
                for _ in 0..4 {
                    entering.push(((rng.next_u64() % m as u64) as usize, rng.uniform(-1.0, 1.0)));
                }
                entering.sort_unstable_by_key(|&(i, _)| i);
            }
            let mut aq = vec![0.0; m];
            for &(i, a) in &entering {
                aq[i] = a;
            }
            let mut bbar = vec![0.0; m];
            pfi.solve_b_raw(&aq, &mut bbar);
            let nz: Vec<usize> = bbar
                .iter()
                .enumerate()
                .filter(|&(_, &v)| v != 0.0)
                .map(|(i, _)| i)
                .collect();
            let pos = nz[(rng.next_u64() % nz.len().max(1) as u64) as usize];
            pfi.update_sparse(pos, entering, &bbar, &nz, 1.0).unwrap();
        }
        assert!(pfi.eta_count() >= 40, "eta file too short to exercise");
        // Pin the krka regression: an eta with an EMPTY idx (unit-column
        // pivot — b̄ carries only the pivot nonzero) is legal, its dot is
        // 0.0, and the pattern gate must not evaluate ilog2(0) (which
        // panicked before the idx_len > 0 guard). The bit-exact checks
        // below then run over a file CONTAINING this eta.
        {
            let pos = 5usize;
            let mut bbar = vec![0.0; m];
            bbar[pos] = 1.25;
            pfi.update_sparse(pos, vec![(pos, 2.5)], &bbar, &[pos], 1.0)
                .unwrap();
            assert!(pfi.etas.last().unwrap().idx.is_empty());
        }

        // Reference: the earlier plain eta pass + LU transpose, hand-rolled.
        let plain_bt = |pfi: &ProductFormFactorization, b: &[f64], x: &mut [f64]| {
            let mut y = b.to_vec();
            for eta in pfi.etas.iter().rev() {
                let mut dot = 0.0;
                for (&i, &v) in eta.idx.iter().zip(&eta.val) {
                    dot += v * y[i as usize];
                }
                y[eta.pos] = (y[eta.pos] - dot) / eta.pivot;
            }
            pfi.lu.solve_transpose(&y, x);
        };
        for seed in [1u64, 2, 3, 77] {
            // Sparse rhs (unit vector = the dual engine's e_p → β)…
            let e = (seed * 41 % m as u64) as usize;
            let mut rhs_v = vec![0.0; m];
            rhs_v[e] = 1.0;
            let (mut got, mut want) = (vec![0.0; m], vec![0.0; m]);
            pfi.solve_bt_raw(&rhs_v, &mut got);
            plain_bt(&pfi, &rhs_v, &mut want);
            for (g, w) in got.iter().zip(&want) {
                assert_eq!(g.to_bits(), w.to_bits(), "unit rhs seed {seed}");
            }
            // …and a dense rhs (fallback path).
            let d: Vec<f64> = (0..m).map(|_| rng.uniform(-1.0, 1.0)).collect();
            pfi.solve_bt_raw(&d, &mut got);
            plain_bt(&pfi, &d, &mut want);
            for (g, w) in got.iter().zip(&want) {
                assert_eq!(g.to_bits(), w.to_bits(), "dense rhs seed {seed}");
            }
        }
    }

    /// The HVec solve stack (`solve_b_hv` / `solve_bt_hv`, raw and
    /// refined, plus `update_hv`) must reproduce the slice stack value for
    /// value across a long update sequence — hyper-sparse solves (m above
    /// the LU hyper threshold), dense blow-ups (dense entering columns)
    /// and the dense-mode fallbacks — and keep the refactor clock in step.
    /// Zero signs are the only allowed difference (an unlisted position is
    /// +0.0 where a dense loop may have written −0.0).
    #[test]
    fn hv_stack_bit_exact() {
        let m = 600;
        let mut rng = XorShift64::new(39);
        let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        for j in 0..m {
            let mut off = 0.0f64;
            // Block-diagonal sparsity (blocks of 6): solves of 1–3-nonzero
            // rhs reach a handful of rows, i.e. stay hyper-sparse.
            for i in (j / 6) * 6..((j / 6) * 6 + 6).min(m) {
                if i != j && rng.next_f64() < 0.4 {
                    let v = rng.uniform(-5.0, 5.0);
                    cols[j].push((i, v));
                    off += v.abs();
                }
            }
            cols[j].push((j, off + 1.0));
            cols[j].sort_unstable_by_key(|&(i, _)| i);
        }
        let cfg = PfiConfig {
            pivot_tol: 1e-10,
            eta_limit: 10_000,
            min_updates: 50,
            clock_scale: 1.0,
            engine: PfiEngine::Sparse,
        };
        let mut a = ProductFormFactorization::build_sparse(cols.clone(), cfg).unwrap();
        let mut b = ProductFormFactorization::build_sparse(cols, cfg).unwrap();
        let same = |g: &[f64], w: &[f64], what: &str| {
            for (i, (g, w)) in g.iter().zip(w).enumerate() {
                assert!(
                    g.to_bits() == w.to_bits() || (*g == 0.0 && *w == 0.0),
                    "{what}: entry {i}: {g:e} vs {w:e}"
                );
            }
        };
        let covers = |h: &HVec, what: &str| {
            if let Some(p) = h.pattern() {
                let mut on = vec![false; h.len()];
                for &i in p {
                    on[i as usize] = true;
                }
                for (i, &v) in h.vals().iter().enumerate() {
                    assert!(v == 0.0 || on[i], "{what}: unlisted nonzero at {i}");
                }
            }
        };
        let mut hv_rhs = HVec::new(m);
        let mut hv_x = HVec::new(m);
        let mut hv_col = HVec::new(m);
        let mut x = vec![0.0; m];
        let mut dense_seen = 0;
        let mut sparse_seen = 0;
        for step in 0..120 {
            // Entering column: mostly 1–3 nonzeros, every 9th dense.
            let mut aq = vec![0.0; m];
            if step % 9 == 4 {
                for v in aq.iter_mut() {
                    *v = rng.uniform(-1.0, 1.0);
                }
            } else {
                for _ in 0..1 + step % 3 {
                    aq[(rng.next_u64() % m as u64) as usize] = rng.uniform(-1.0, 1.0);
                }
            }
            hv_col.clear();
            for (i, &v) in aq.iter().enumerate().rev() {
                if v != 0.0 {
                    hv_col.set(i, v); // descending: exercises the sort
                }
            }
            for refine in [false, true] {
                if refine {
                    a.solve_b(&aq, &mut x);
                } else {
                    a.solve_b_raw(&aq, &mut x);
                }
                b.solve_b_hv(&hv_col, &mut hv_x, refine);
                same(hv_x.vals(), &x, "ftran");
                covers(&hv_x, "ftran");
                if hv_x.is_dense() {
                    dense_seen += 1;
                } else {
                    sparse_seen += 1;
                }
            }
            // Unit and 2-sparse BTRANs (β-shaped rhs), raw and refined.
            for refine in [false, true] {
                let e = (rng.next_u64() % m as u64) as usize;
                let mut rhs = vec![0.0; m];
                rhs[e] = 1.0;
                if step % 2 == 0 {
                    rhs[(e * 7 + 3) % m] = -0.5;
                }
                hv_rhs.clear();
                for (i, &v) in rhs.iter().enumerate() {
                    if v != 0.0 {
                        hv_rhs.set(i, v);
                    }
                }
                let mut want = vec![0.0; m];
                if refine {
                    a.solve_bt(&rhs, &mut want);
                } else {
                    a.solve_bt_raw(&rhs, &mut want);
                }
                b.solve_bt_hv(&hv_rhs, &mut hv_x, refine);
                same(hv_x.vals(), &want, "btran");
                covers(&hv_x, "btran");
            }
            // Pivot on a nonzero of the (refined) b̄; both sides update.
            a.solve_b(&aq, &mut x);
            b.solve_b_hv(&hv_col, &mut hv_x, true);
            let nz: Vec<usize> = (0..m).filter(|&i| x[i].abs() > 1e-3).collect();
            if nz.is_empty() {
                continue;
            }
            let pos = nz[(rng.next_u64() % nz.len() as u64) as usize];
            let sa = a.update(pos, &aq, &x, 1.0).unwrap();
            let sb = b.update_hv(pos, &mut hv_col, &mut hv_x, 1.0).unwrap();
            assert_eq!(sa, sb, "update status step {step}");
            assert_eq!(a.basis_nz, b.basis_nz, "basis column step {step}");
            assert_eq!(a.etas.len(), b.etas.len());
            if let (Some(ea), Some(eb)) = (a.etas.last(), b.etas.last()) {
                assert_eq!((ea.pos, &ea.idx), (eb.pos, &eb.idx), "eta step {step}");
                same(&ea.val, &eb.val, "eta values");
            }
            assert_eq!(a.eta_ticks.get().to_bits(), b.eta_ticks.get().to_bits());
        }
        assert!(
            dense_seen > 0 && sparse_seen > 0,
            "{dense_seen} {sparse_seen}"
        );
    }
}
