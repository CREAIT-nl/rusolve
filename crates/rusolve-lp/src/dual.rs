//! Bounded-variable DUAL simplex engine over the shifted-bounds form.
//!
//! Lever: `SimplexConfig::dual` / CLI
//! `RUSOLVE_DUAL=1`, default OFF — the primal path (`simplex.rs`) is untouched
//! and stays bit-identical when the lever is off.
//!
//! Contents:
//! * [`DualStandard`] (shifted-bounds form + Wächter–Bixby cost-sign dual
//!   start + cost shifting + E-row fixed artificials), `initial_basis` /
//!   `initial_x_b`, the start-infeasibility diagnostic, and lever plumbing.
//!   `solve_dual_raw` emits diagnostics and falls back to the primal solver.
//! * Classical dual engine (leaving rule + dual ratio test with flips).
//! * DuPh1 Σ-infeasibility overlay.
//!
//! Provenance (HiGHS anchors are MIT-licensed source, used for
//! boundaries/formulas only):
//! * Shifted-bounds form, cost-sign dual start, cost shifting: Wächter &
//!   Bixby, Math. Prog. Comp. 121:231 (2006) — dual start / cost shifting.
//! * Bounded dual simplex ratio test with bound flips: Fortuny-Amat & McCarl,
//!   Math. Prog. Study 15 (1981); Gill–Murray–Wright §4.6; HiGHS
//!   `HDualSimplex::dualRatioTest` (flip-then-enter event ordering).
//! * Dual simplex correctness argument: Chvátal, *Linear Programming* (1983)
//!   §4.7.
//! * Σ-infeasibility dual phase 1 (DuPh1): Klein & Gondzio, Math. Prog. Comp.
//!   2023; HiGHS `HDualSimplex` phase-1 flag.
//! * Fixed basic artificials for E rows + delete-once pricing: textbook
//!   bounded simplex (Bertsimas–Tsitsiklis) + our artificial-deletion
//!   convention and the `ub = 0` pin of leftover basic artificials (found on
//!   boeing2, where a basic artificial left at a degenerate 0 by phase 1
//!   was pivoted off zero in phase 2).
//!
//! No `unsafe`. f64 throughout.

use crate::FactorError;
use crate::hvec::HVec;
use crate::lu_ft::BasisFactorization;
use crate::pfi::UpdateStatus;
use crate::simplex::{
    Basis, BasisStatus, SimplexConfig, SimplexError, SimplexSolution, SimplexSolver, SimplexStatus,
    farkas_ok,
};
use rusolve_core::Lp;

/// Column role (dual-side analogue of the private `Kind` in simplex.rs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Structural,
    Slack,
    /// E-row fixed artificial: bounds [0, 0], basic at start, deleted from
    /// pricing on first exit (the delete-once convention).
    Artificial,
}

/// Dual-feasible shifted-bounds standard form:
/// `min cᵀx̃ s.t. Ax̃ = b, 0 ≤ x̃ ≤ ũ` plus slack / fixed-artificial columns.
/// NO row negation and NO rhs-sign preprocessing — the dual engine handles
/// negative rhs natively (negation is primal-phase-1 machinery
/// and would perturb the initial infeasibility profile we measure honestly).
#[derive(Debug, Clone)]
pub struct DualStandard {
    pub m: usize,
    /// Column-sparse matrix (the source of truth; there is no dense
    /// `cols` twin): column `v` owns `csc_val[csc_ptr[v]..csc_ptr[v + 1]]`
    /// with row indices in `csc_row`, ascending; exact zeros are never
    /// stored (same convention as `Standard::csc_*` in simplex.rs).
    pub csc_ptr: Vec<usize>,
    pub csc_row: Vec<usize>,
    pub csc_val: Vec<f64>,
    pub b: Vec<f64>,
    /// Shifted upper bound ũ_j ∈ (0, +∞]; 0 for the fixed E-row artificials.
    pub ub: Vec<f64>,
    pub kind: Vec<Kind>,
    pub num_structural: usize,
    /// TRUE cost sign-normalized to MINIMIZATION (objective·sign). Effective
    /// cost is `cost[j] + s[j]`.
    pub cost: Vec<f64>,
    /// Cost shifting: per-variable nonnegative shift, nonzero
    /// only while the variable is nonbasic on the lower side with ũ = +∞.
    /// Released (zeroed) when the variable enters the basis.
    pub s: Vec<f64>,
    /// Original variable → its structural columns (1, or 2 for split frees).
    pub cols_of: Vec<Vec<usize>>,
    /// Per-original bound shift (finite lb, or ub for flip-shifted).
    pub shift: Vec<f64>,
    /// Per-original: reconstruct as shift − x̃ (flip-shifted MI bounds).
    pub flip: Vec<bool>,
    /// Nonbasic side: true ⇒ variable sits at ũ, false ⇒ at 0. Set at the
    /// start by the cost-sign rule, maintained by flips/pivots.
    pub at_ub: Vec<bool>,
    /// Cost perturbation (the stall rung): signed sub-tolerance
    /// cost offsets added wherever the PIVOT LOOP prices (`cost + s + pert`).
    /// All-zero until the escalation ladder engages it on a plateau/stall
    /// (esc ≥ 2·bland_after); it then lives for the rest of the solve. The
    /// optimality audit prices TRUE costs (`cost + s`, never `pert`), so it is
    /// the audit cleanup pass by construction: Optimal is only ever certified
    /// at true-cost KKT. The cutoff bound prices d at TRUE costs and is
    /// valid for ANY π (weak duality), so the perturbed π only makes it
    /// marginally looser, never unsound.
    pub pert: Vec<f64>,
}

impl DualStandard {
    /// Total column count (structural + slack + artificial).
    #[inline]
    fn ncol(&self) -> usize {
        self.csc_ptr.len() - 1
    }

    /// Sparse `(row, value)` list of column `v`, ascending by row with exact
    /// zeros dropped (the shape `ProductFormFactorization::build_sparse`
    /// wants; bit-identical values to `gather_col(v, 1.0, ..)`).
    fn sparse_col(&self, v: usize) -> Vec<(usize, f64)> {
        (self.csc_ptr[v]..self.csc_ptr[v + 1])
            .filter(|&k| self.csc_val[k] != 0.0)
            .map(|k| (self.csc_row[k], self.csc_val[k]))
            .collect()
    }

    /// The only dense-column access API (twin of `Standard::gather_col` in
    /// simplex.rs): zero-fill `buf` (length m) and scatter column `v` scaled
    /// by `sgn` (±1.0, exact). `fill(0.0)` + scatter reproduces the dense
    /// `Vec<f64>` the deleted `cols` field held, structurally not
    /// approximately; the O(m) fill is a memset next to the
    /// FTRAN it feeds.
    fn gather_col(&self, v: usize, sgn: f64, buf: &mut [f64]) {
        buf.fill(0.0);
        for k in self.csc_ptr[v]..self.csc_ptr[v + 1] {
            buf[self.csc_row[k]] = sgn * self.csc_val[k];
        }
    }
}

/// Cost-shift target ε_s (our tolerance convention): ptol/10, below ptol so
/// cleanup pricing cannot mistake a shifted dual for a real infeasibility.
pub fn eps_s(config: &SimplexConfig) -> f64 {
    config.pricing_tol / 10.0
}

/// Dual refactor cadence (`RUSOLVE_DUAL_RFCAD=0` = off): the dual's
/// factors fire the cost clock at 4 × build work instead of the calibrated
/// 1 ×. A per-phase census against HiGHS found INVERT the
/// largest excess over HiGHS on qap10/ex10 and second on irish: we rebuilt
/// every 11-31 pivots where HiGHS's synthetic clock reinverts every
/// 100-316, because our `build_work` under-counts the real kernel cost
/// (refactor wall 1.9-3.7× eta-solve wall at the clock's break-even point).
/// Same clock as HiGHS (HEkkDual.cpp:2192-2195 / HEkk.cpp `updateFactor`;
/// MIT source read as reference text, algorithm re-derived), rescaled.
/// HiGHS's 50-update floor (`kSyntheticTickReinversionMinUpdateCount`) is
/// NOT adopted: with our PFI's dense etas it sent ex10 into a 142k-iteration
/// phase-1 walk. Measured effect: instructions per
/// iteration stay flat (saved rebuilds are repaid by deeper eta solves);
/// the win is fewer iterations under ACC (47 relaxed MIPLIB LPs Σ −16 %, irish −19 %).
/// The primal keeps the calibrated clock.
const DUAL_CLOCK_SCALE: f64 = 4.0;
const DUAL_RFCAD_DEFAULT: bool = true;

/// The [`PfiConfig`] the dual's factor builds use (see
/// [`DUAL_CLOCK_SCALE`]).
pub(crate) fn dual_pfi_config(config: &SimplexConfig) -> crate::PfiConfig {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    let on = *ON.get_or_init(|| match std::env::var("RUSOLVE_DUAL_RFCAD").as_deref() {
        Ok("0") => false,
        Ok(_) => true,
        Err(_) => DUAL_RFCAD_DEFAULT,
    });
    let mut pfi = config.pfi;
    if on {
        pfi.clock_scale = DUAL_CLOCK_SCALE;
    }
    pfi
}

/// Bound used for the per-variable shift: finite lb, else finite ub (flip),
/// else 0 (fully free ⇒ split, no shift). Same semantics as simplex::shift_of.
fn shift_of(lp: &Lp, j: usize) -> f64 {
    if lp.lower[j].is_finite() {
        lp.lower[j]
    } else if lp.upper[j].is_finite() {
        lp.upper[j]
    } else {
        0.0
    }
}

/// Convert a (presolved) [`Lp`] into the dual standard form with the
/// Wächter–Bixby cost-sign start recorded in `at_ub` / `s`.
pub fn to_dual_standard(lp: &Lp, config: &SimplexConfig) -> Result<DualStandard, SimplexError> {
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
        if lb == f64::NEG_INFINITY && lp.upper[j] == f64::NEG_INFINITY {
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
    for (i, &sn) in lp.senses.iter().enumerate() {
        if !matches!(sn, 'L' | 'G' | 'E') {
            return Err(SimplexError::Unsupported(format!("row {i} sense '{sn}'")));
        }
    }
    let es = eps_s(config);
    let sign = if lp.maximize { -1.0 } else { 1.0 };

    // Bound shifts. NO row negation.
    // The b-shift scatters over each column's nnz (ascending rows — the
    // dense loop's order; stored zeros were ±0.0 no-ops).
    let shifts: Vec<f64> = (0..n).map(|j| shift_of(lp, j)).collect();
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
    let mut b = lp.rhs.clone();
    for (j, &sh) in shifts.iter().enumerate() {
        if sh != 0.0 {
            for k in ptr[j]..ptr[j + 1] {
                let v = rval[k];
                if v != 0.0 {
                    b[ridx[k] as usize] -= v * sh;
                }
            }
        }
    }

    // Columns go straight into the CSC arrays: entry k of column j
    // comes from the lp CSC column at ascending row order, exactly the
    // entries the old zero-fill + scatter + derive pass kept — stored ±0.0
    // are dropped, negated copies negate only the
    // kept nonzeros (`-x` is exact, `-(±0.0)` was dropped either way).
    let mut csc_ptr: Vec<usize> = Vec::with_capacity(n + 2 * m + 1);
    let mut csc_row: Vec<usize> = Vec::new();
    let mut csc_val: Vec<f64> = Vec::new();
    let mut ub: Vec<f64> = Vec::new();
    let mut kind: Vec<Kind> = Vec::new();
    let mut cost: Vec<f64> = Vec::new();
    let mut at_ub: Vec<bool> = Vec::new();
    let mut cols_of: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut flip = vec![false; n];
    csc_ptr.push(0);

    {
        // Start-side rule: cost ≥ 0 ⇒ lower side (x̃ = 0);
        // cost < 0 with finite ũ > 0 ⇒ upper side (x̃ = ũ); cost < 0 with
        // ũ = +∞ ⇒ lower side + cost shift (applied after this loop).
        // The bound-choice rule IS the W&B dual start: with the slack basis
        // (π = 0) the effective reduced cost equals the effective cost.
        let push = |entries: &[(usize, f64)],
                    ubj: f64,
                    c: f64,
                    csc_ptr: &mut Vec<usize>,
                    csc_row: &mut Vec<usize>,
                    csc_val: &mut Vec<f64>,
                    ub: &mut Vec<f64>,
                    kind: &mut Vec<Kind>,
                    cost: &mut Vec<f64>,
                    at_ub: &mut Vec<bool>| {
            let at_upper = c < 0.0 && ubj.is_finite() && ubj > 0.0;
            // Column index BEFORE the closing ptr is appended (csc_ptr
            // starts holding the 0 sentinel — the simplex push_col! idiom).
            let idx = csc_ptr.len() - 1;
            for &(r, x) in entries {
                csc_row.push(r);
                csc_val.push(x);
            }
            csc_ptr.push(csc_row.len());
            ub.push(ubj);
            kind.push(Kind::Structural);
            cost.push(c);
            at_ub.push(at_upper);
            idx
        };
        for j in 0..n {
            let entries: Vec<(usize, f64)> = (ptr[j]..ptr[j + 1])
                .filter(|&k| rval[k] != 0.0)
                .map(|k| (ridx[k] as usize, rval[k]))
                .collect();
            let (lb, ubj) = (lp.lower[j], lp.upper[j]);
            if lb.is_finite() {
                let c = lp.objective[j] * sign;
                let col = push(
                    &entries,
                    ubj - lb,
                    c,
                    &mut csc_ptr,
                    &mut csc_row,
                    &mut csc_val,
                    &mut ub,
                    &mut kind,
                    &mut cost,
                    &mut at_ub,
                );
                cols_of[j].push(col);
            } else if ubj.is_finite() {
                // flip-shift: x = ub − x̃ ⇒ column negated, cost negated.
                let c = -lp.objective[j] * sign;
                let neg: Vec<(usize, f64)> = entries.iter().map(|&(r, x)| (r, -x)).collect();
                let col = push(
                    &neg,
                    f64::INFINITY,
                    c,
                    &mut csc_ptr,
                    &mut csc_row,
                    &mut csc_val,
                    &mut ub,
                    &mut kind,
                    &mut cost,
                    &mut at_ub,
                );
                cols_of[j].push(col);
                flip[j] = true;
            } else {
                // fully free: x = x⁺ − x⁻, both ≥ 0 (existing free-column plumbing).
                let c = lp.objective[j] * sign;
                let col = push(
                    &entries,
                    f64::INFINITY,
                    c,
                    &mut csc_ptr,
                    &mut csc_row,
                    &mut csc_val,
                    &mut ub,
                    &mut kind,
                    &mut cost,
                    &mut at_ub,
                );
                cols_of[j].push(col);
                let neg: Vec<(usize, f64)> = entries.iter().map(|&(r, x)| (r, -x)).collect();
                let col = push(
                    &neg,
                    f64::INFINITY,
                    -c,
                    &mut csc_ptr,
                    &mut csc_row,
                    &mut csc_val,
                    &mut ub,
                    &mut kind,
                    &mut cost,
                    &mut at_ub,
                );
                cols_of[j].push(col);
            }
        }
    }
    let num_structural = csc_ptr.len() - 1;

    // Cost-shift init: every starting-infeasible infinite-side
    // variable (effective cost < 0 at the lower side, no bound to flip to)
    // gets s := ε_s − c so its effective reduced cost is ε_s > 0 at π = 0.
    // Every row contributes exactly one slack/artificial column below, and
    // those are never shifted (cost 0), so one allocation covers the final
    // length (the old double resize landed on the same number).
    let mut s = vec![0.0; num_structural + m];
    for v in 0..num_structural {
        if cost[v] < 0.0 && !at_ub[v] && ub[v] == f64::INFINITY {
            s[v] = es - cost[v];
        }
    }

    // Slacks (L: +e_i, G: −e_i; bounds [0, +∞], cost 0, basic) and E-row
    // fixed artificials (+e_i, bounds [0, 0], cost 0, basic).
    for i in 0..m {
        if lp.senses[i] != 'E' {
            csc_row.push(i);
            csc_val.push(if lp.senses[i] == 'L' { 1.0 } else { -1.0 });
            csc_ptr.push(csc_row.len());
            ub.push(f64::INFINITY);
            kind.push(Kind::Slack);
            cost.push(0.0);
            at_ub.push(false);
        }
    }
    for i in 0..m {
        if lp.senses[i] == 'E' {
            csc_row.push(i);
            csc_val.push(1.0);
            csc_ptr.push(csc_row.len());
            ub.push(0.0); // fixed: the basic value IS the row violation
            kind.push(Kind::Artificial);
            cost.push(0.0);
            at_ub.push(false);
        }
    }

    Ok(DualStandard {
        // Perturbation vector: all-zero until the stall rung
        // engages it (see the `pert` field docs).
        pert: vec![0.0; cost.len()],
        m,
        csc_ptr,
        csc_row,
        csc_val,
        b,
        ub,
        kind,
        num_structural,
        cost,
        s,
        cols_of,
        shift: shifts,
        flip,
        at_ub,
    })
}

/// Bring a cached [`DualStandard`] up to date with bounds-only
/// edits of its `Lp` (the MIP node loop's between-calls change; see
/// [`StdCache`] for the gen contract). Returns `false` when the edit is NOT
/// absorbable — a column changed its finiteness CATEGORY (slot count / flip
/// mapping) or the bounds went invalid (NaN, ub < lb) — and the caller then
/// rebuilds through `to_dual_standard`, whose validation reports properly.
/// The refresh never tightens: sibling nodes can relax bounds too; only the
/// per-column CATEGORY has to match what the cache was built with.
///
/// Incremental updates, from the build equations in `to_dual_standard`:
///   · `shift_j = shift_of(lp, j)`; `b_i −= a_ij·Δshift_j` over the changed
///     column's nonzeros only (the full build scatters every column).
///   · ũ_slot = ub − lb for the lb-finite category; flip/free slots carry ∞.
///
/// The old bounds are re-derived from the built form itself (`shift`, `flip`
/// and `cols_of[j]`), so the cache stores no shadow arrays. Solve state
/// (`pert`, `s`, `at_ub`) is reset to EXACTLY a fresh build's start: the
/// warm import overwrites `at_ub` from the snapshot anyway, but `s`'s cost-shift
/// init and the W&B start flags are read before that.
fn refresh_dual_standard(lp: &Lp, config: &SimplexConfig, st: &mut DualStandard) -> bool {
    let es = eps_s(config);
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
    for j in 0..lp.num_col {
        let (lb, ub) = (lp.lower[j], lp.upper[j]);
        if lb.is_nan() || ub.is_nan() || ub < lb {
            return false; // let to_dual_standard produce the honest error
        }
        let free = st.cols_of[j].len() == 2; // built as an x⁺/x⁻ split
        let flip = st.flip[j]; // built with lb = −∞, ub finite
        let lb_finite_slot = !free && !flip; // built with lb finite
        if free != (!lb.is_finite() && !ub.is_finite()) {
            return false; // free ⇔ both bounds infinite
        }
        if lb_finite_slot != lb.is_finite() {
            return false; // lb finiteness fixes the slot's direction
        }
        let sh_new = shift_of(lp, j);
        let sh_old = st.shift[j];
        if sh_new != sh_old {
            let d = sh_new - sh_old;
            for k in ptr[j]..ptr[j + 1] {
                let v = rval[k];
                if v != 0.0 {
                    st.b[ridx[k] as usize] -= v * d;
                }
            }
            st.shift[j] = sh_new;
        }
        if lb_finite_slot {
            st.ub[st.cols_of[j][0]] = ub - lb;
        }
    }
    // Reset the solve state to the fresh-build start: pert zero, `s` at the
    // cost-shift init, `at_ub` at the W&B start rule — mirroring the build's own
    // loops, which are SLOT-wise: at-upper tests the slot box width ũ
    // (finite and > 0), and the cost shift applies to every slot with
    // cost < 0, ũ = ∞ and lower start — flip slots and both halves of a
    // free split included (a negative-cost free half IS cost-shifted).
    st.pert.iter_mut().for_each(|p| *p = 0.0);
    st.s.iter_mut().for_each(|s| *s = 0.0);
    st.at_ub.iter_mut().for_each(|a| *a = false);
    for v in 0..st.num_structural {
        let c = st.cost[v];
        let at_upper = c < 0.0 && st.ub[v].is_finite() && st.ub[v] > 0.0;
        st.at_ub[v] = at_upper;
        if c < 0.0 && !at_upper && st.ub[v] == f64::INFINITY {
            st.s[v] = es - c;
        }
    }
    true
}

/// Initial basis: every row's slack (L/G) or fixed artificial (E) ⇒ B = ±I.
/// Returns (column per row slot, diagonal sign ±1).
pub fn initial_basis(lp: &Lp, st: &DualStandard) -> (Vec<usize>, Vec<f64>) {
    let m = st.m;
    let mut basis = vec![usize::MAX; m];
    let mut diag = vec![0.0; m];
    // Slack columns were pushed for non-E rows in row order; artificials for
    // E rows after them.
    let n_slacks = lp.senses.iter().filter(|&&sn| sn != 'E').count();
    let art_start = st.num_structural + n_slacks;
    let mut next_slack = st.num_structural;
    let mut next_art = art_start;
    for i in 0..m {
        if lp.senses[i] == 'E' {
            basis[i] = next_art;
            next_art += 1;
            diag[i] = 1.0;
        } else {
            basis[i] = next_slack;
            next_slack += 1;
            diag[i] = if lp.senses[i] == 'L' { 1.0 } else { -1.0 };
        }
    }
    (basis, diag)
}

/// Initial basic values: x_b = B⁻¹(b − Σ_{nb at ũ} ũ_j a_j).
/// The ub-fixed nonbasic contribution is NOT optional (the primal engine's
/// unit tests caught exactly that bug). Do NOT clamp negatives: they are the
/// DuPh input, not noise. Only |v| < 1e-12 noise snaps to 0.
pub fn initial_x_b(st: &DualStandard, basis: &[usize], diag: &[f64]) -> Vec<f64> {
    let m = st.m;
    let total = st.ncol();
    let mut rhs = st.b.clone();
    let mut in_b = vec![false; total];
    for &v in basis.iter() {
        in_b[v] = true;
    }
    for (v, &ib) in in_b.iter().enumerate() {
        if !ib && st.at_ub[v] && st.ub[v].is_finite() {
            for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                rhs[st.csc_row[k]] -= st.ub[v] * st.csc_val[k];
            }
        }
    }
    // B = diag(±1) ⇒ B⁻¹r = diag·r.
    let mut x_b = vec![0.0; m];
    for i in 0..m {
        let mut r = rhs[i] * diag[i];
        if r.abs() < 1e-12 {
            r = 0.0;
        }
        x_b[i] = r;
    }
    x_b
}

/// Sparse dot πᵀa_v over the CSC column v (dual.rs twin of simplex.rs's
/// private `sparse_dot`; same bit-exactness argument).
#[inline]
pub fn csc_dot(x: &[f64], st: &DualStandard, v: usize) -> f64 {
    let (lo, hi) = (st.csc_ptr[v], st.csc_ptr[v + 1]);
    let mut acc = 0.0f64;
    for k in lo..hi {
        acc += x[st.csc_row[k]] * st.csc_val[k];
    }
    acc
}

/// Row-wise (CSR) copy of the dual standard form's matrix, built once per
/// solve for the pivotal-row PRICE `α = ρᵀ A_N` (HiGHS keeps the same
/// `ar_matrix_`: `HEkk::tableauRowPrice`, simplex/HEkk.cpp:2845, row-wise
/// copy built in `HEkk::initialiseSimplexLpBasisAndFactor`/`ar_matrix_`).
struct RowMat {
    ptr: Vec<usize>,
    col: Vec<u32>,
    val: Vec<f64>,
}

impl RowMat {
    /// Transpose of the CSC arrays; entries of a row come out in ascending
    /// column order, entries of a column in ascending row order (so the
    /// per-column accumulation order of [`price_row`] equals [`csc_dot`]).
    fn from_dual(st: &DualStandard) -> Self {
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
        let mut col = vec![0u32; nnz];
        let mut val = vec![0.0f64; nnz];
        for v in 0..st.ncol() {
            for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                let r = st.csc_row[k];
                col[next[r]] = v as u32;
                val[next[r]] = st.csc_val[k];
                next[r] += 1;
            }
        }
        RowMat { ptr, col, val }
    }
}

/// Set of row positions that are (or were since the last rebuild) primal
/// infeasible basics — the hyper-sparse CHUZR candidate set of
/// Hall & McKinnon (2005), *Hyper-sparsity in the revised simplex method and
/// how to exploit it*, Comput. Optim. Appl. 32:259–283, as used by HiGHS
/// `simplex/HEkkDualRHS.cpp` `chooseNormal` :37 / `updateInfeasList` :375 /
/// `createInfeasList` :436 (MIT; `reference/HiGHS`, read as a reference
/// text, reimplemented for this engine, not translated).
///
/// Like HiGHS the list is an OVER-APPROXIMATION: rows are added when an
/// x_b change makes them violate, never removed; the CHUZR scan
/// re-checks each entry, so stale entries cost one comparison and the
/// invariant `{violating rows} ⊆ list` holds between rebuilds. Rebuild =
/// the full 0..m scan (which the caller falls back to when the list grows
/// past a fraction of m, or after anything recomputed x_b wholesale).
/// HiGHS instead prunes the list with a merit cutoff (its hyper mode); we
/// keep it whole because this engine's plateau ladder needs an exact Σv and
/// CHUZR needs full violator coverage, and because the β/b̄ patterns are not
/// exposed by the PFI, maintenance piggybacks on the existing
/// x_b update loops behind a `b̄[i] != 0` test instead of an HVector index
/// list.
struct InfeasList {
    list: Vec<usize>,
    mark: Vec<u8>,
}

impl InfeasList {
    fn new(m: usize) -> Self {
        InfeasList {
            list: Vec::new(),
            mark: vec![0; m],
        }
    }
    /// Drop the current entries (marks released); the caller rebuilds from a
    /// full scan that `push`es the violators back in.
    fn begin_rebuild(&mut self) {
        for &i in self.list.iter() {
            self.mark[i] = 0;
        }
        self.list.clear();
    }
    #[inline]
    fn push(&mut self, i: usize) {
        if self.mark[i] == 0 {
            self.mark[i] = 1;
            self.list.push(i);
        }
    }
}

/// Hyper-sparse passes: the pivot's x_B update clamps
/// round-off negatives `x_b[i] ∈ (−FEAS_TOL, 0)` to 0 on EVERY row, so a
/// pattern-only update must still visit each row that could hold such a
/// value. Only two writers can leave one behind without clamping it: the
/// flip-batch updates (rows recorded here) and the wholesale recomputes
/// (audit, verification refresh, accuracy-trouble computePrimal, backtrack — `all`);
/// every other x_b write lands in `[0, ub]`. With `all` clear, clamping the
/// b̄ pattern plus the recorded rows reproduces the all-row loop exactly.
///
/// Hyper CHUZR: with `track` on, the same writers also log every row
/// whose x_b moved since the last CHUZR (`chg`, deduplicated), or
/// `chg_all` for a dense/wholesale write — the [`HyperChuzr`] dirty set.
struct XbClamp {
    all: bool,
    rows: Vec<usize>,
    cap: usize,
    track: bool,
    chg_all: bool,
    chg: Vec<u32>,
    chg_mark: Vec<bool>,
}

impl XbClamp {
    fn new(m: usize, track: bool) -> Self {
        XbClamp {
            all: true, // the initial x_b was computed wholesale
            rows: Vec::new(),
            cap: (m / 10).max(1),
            track,
            chg_all: true,
            chg: Vec::new(),
            chg_mark: if track { vec![false; m] } else { Vec::new() },
        }
    }

    /// A wholesale x_b write is pending/just happened (audit, refresh,
    /// computePrimal, backtrack).
    fn wholesale(&mut self) {
        self.all = true;
        self.chg_all = true;
    }

    /// Record the rows a flip batch moved (its b̄ pattern), or `all` for a
    /// dense b̄.
    fn note(&mut self, pat: Option<&[u32]>) {
        self.note_chg(pat);
        match pat {
            Some(p) if !self.all && self.rows.len() + p.len() <= self.cap => {
                self.rows.extend(p.iter().map(|&i| i as usize));
            }
            _ => {
                self.all = true;
                self.rows.clear();
            }
        }
    }

    /// Hyper-CHUZR dirty rows: x_b moved on `pat` (None ⇒ every row).
    fn note_chg(&mut self, pat: Option<&[u32]>) {
        if !self.track || self.chg_all {
            return;
        }
        match pat {
            Some(p) if self.chg.len() + p.len() <= self.cap => {
                for &i in p {
                    if !self.chg_mark[i as usize] {
                        self.chg_mark[i as usize] = true;
                        self.chg.push(i);
                    }
                }
            }
            _ => self.chg_all = true,
        }
    }

    /// Forget the dirty rows (a full CHUZR scan just re-read every row).
    fn reset_chg(&mut self) {
        for &i in &self.chg {
            self.chg_mark[i as usize] = false;
        }
        self.chg.clear();
        self.chg_all = false;
    }
}

/// Neumaier-compensated running sum (the hyper CHUZR's incremental Σv: a
/// long add/subtract chain must not drift against the dense re-sum).
#[derive(Clone, Copy, Default)]
struct CompSum {
    hi: f64,
    lo: f64,
}

impl CompSum {
    #[inline]
    fn add(&mut self, x: f64) {
        let t = self.hi + x;
        if self.hi.abs() >= x.abs() {
            self.lo += (self.hi - t) + x;
        } else {
            self.lo += (x - t) + self.hi;
        }
        self.hi = t;
    }

    #[inline]
    fn value(&self) -> f64 {
        self.hi + self.lo
    }
}

/// Hyper-sparse CHUZR: leaving-row choice
/// without the O(m) violation scan. Algorithm after Hall & McKinnon 2005
/// §5 (hyper-sparse CHUZR: a candidate set above a cutoff, refreshed by a
/// full pass only when exhausted) and HiGHS `HEkkDualRHS::chooseNormal` /
/// `updateInfeasList` (MIT; read as a reference, reimplemented here).
///
/// Invariant: every row whose current score exceeds `cutoff` is in `set`
/// (rows outside it scored ≤ `cutoff` when last read and have not moved
/// since — any x_b/ub_b/DSE-weight change on a row puts it in the
/// [`XbClamp`] dirty log or invalidates the structure). With the pivot
/// loop's tie = 0 the CHUZR result is the order-independent lexicographic
/// max of (score, −basis index), so the best set member equals the full
/// scan's choice exactly. Only Σv differs: it is maintained incrementally
/// (compensated) instead of re-summed in row order — a rounding-level
/// change to the plateau ladder's V test, hence the lever.
struct HyperChuzr {
    on: bool,
    valid: bool,
    k: usize,
    cutoff: f64,
    set: Vec<u32>,
    in_set: Vec<bool>,
    vcache: Vec<f64>,
    vtot: CompSum,
    nviol: usize,
    prev: usize,
    since: usize,
    exhausted: bool,
    scratch: Vec<f64>,
    /// Adaptive fallback (see [`HyperChuzr::gate`]): running mean of
    /// hits per build, visits left in the current cool-down, its next
    /// length, and "the next build starts a probe".
    hits_ema: f64,
    cool: usize,
    cool_len: usize,
    probe: bool,
    // stats
    n_hit: usize,
    n_build: usize,
    set_sum: usize,
    n_off: usize,
}

impl HyperChuzr {
    fn new(m: usize, on: bool) -> Self {
        HyperChuzr {
            on,
            valid: false,
            k: (m / HCHUZR_K_DIV).clamp(HCHUZR_K_MIN, HCHUZR_K_MAX),
            cutoff: f64::INFINITY,
            set: Vec::new(),
            in_set: if on { vec![false; m] } else { Vec::new() },
            vcache: if on { vec![0.0; m] } else { Vec::new() },
            vtot: CompSum::default(),
            nviol: 0,
            prev: usize::MAX,
            since: 0,
            exhausted: false,
            scratch: Vec::new(),
            hits_ema: HCHUZR_EMA_INIT,
            cool: 0,
            cool_len: HCHUZR_COOL_MIN,
            probe: false,
            n_hit: 0,
            n_build: 0,
            set_sum: 0,
            n_off: 0,
        }
    }

    /// Adaptive fallback: `false` = this CHUZR visit runs the plain
    /// full scan WITHOUT the hyper path (no `choose`, no `build`). A build
    /// costs one extra O(m) pass plus a selection on top of the full scan
    /// it follows, so a candidate set that only survives a hit or two is a
    /// net loss (supportcase10: dense b̄ / flip batches dirty the set every
    /// few iterations). While the running mean of hits per build is below
    /// [`HCHUZR_MIN_HITS`] the hyper path sits out a cool-down (doubling
    /// on each consecutive failed probe, reset by a probe that pays), then
    /// probes again — it re-enters when the iterations turn sparse.
    /// Same idea as HiGHS dropping `HEkkDualRHS` to DENSE mode (`workCount
    /// = -numRow`) when the candidate list is too large or `col_aq_density`
    /// ≥ 0.05 at `createInfeasList` and retrying at the next rebuild
    /// (simplex/HEkkDualRHS.cpp:436-498, MIT, read as reference); the cost
    /// signal here is the measured hit rate instead of the density.
    /// Exactness is unaffected: the full scan is the reference path.
    fn gate(&mut self) -> bool {
        if self.cool > 0 {
            self.cool -= 1;
            if self.cool == 0 {
                self.probe = true;
            }
            return false;
        }
        true
    }

    /// Fold the life of the structure that just died (`since` hits) into
    /// the running mean; `true` = stop and cool down instead of building.
    fn retire(&mut self) -> bool {
        if self.probe {
            // A probe's verdict is its own first life only.
            self.probe = false;
            self.hits_ema = HCHUZR_MIN_HITS;
            return false;
        }
        if self.n_build == 0 {
            return false;
        }
        self.hits_ema = (1.0 - HCHUZR_EMA_W) * self.hits_ema + HCHUZR_EMA_W * self.since as f64;
        if self.hits_ema < HCHUZR_MIN_HITS {
            self.cool = self.cool_len;
            self.cool_len = (self.cool_len * 2).min(HCHUZR_COOL_MAX);
            self.valid = false;
            self.n_off += 1;
            return true;
        }
        if self.hits_ema >= 2.0 * HCHUZR_MIN_HITS {
            self.cool_len = HCHUZR_COOL_MIN;
        }
        false
    }

    #[inline]
    fn score(v: f64, w: f64, dse_on: bool) -> f64 {
        if dse_on { v * v / w } else { v }
    }

    /// Re-read every row after a full CHUZR scan: violation cache, Σv,
    /// violator count, and the cutoff = the (k+1)-th largest violator score
    /// (−∞ when at most k rows violate, so the set is then every violator).
    fn build(&mut self, x_b: &[f64], ub_b: &[f64], dse_w: &[f64], dse_on: bool) {
        if self.retire() {
            return;
        }
        if self.exhausted && self.since < self.k / 8 {
            // Exhausted fast: the top-k churns faster than it is consumed.
            self.k = (self.k * 2).min(HCHUZR_K_MAX.max(self.vcache.len() / 8));
        }
        self.exhausted = false;
        for &i in &self.set {
            self.in_set[i as usize] = false;
        }
        self.set.clear();
        self.scratch.clear();
        self.vtot = CompSum::default();
        for i in 0..x_b.len() {
            let v = violation(ub_b[i], x_b[i]);
            self.vcache[i] = v;
            self.vtot.add(v);
            if v > FEAS_TOL {
                self.scratch.push(Self::score(v, dse_w[i], dse_on));
            }
        }
        self.nviol = self.scratch.len();
        self.cutoff = if self.scratch.len() > self.k {
            let k = self.k;
            let (_, kth, _) = self
                .scratch
                .select_nth_unstable_by(k, |a, b| b.total_cmp(a));
            *kth
        } else {
            f64::NEG_INFINITY
        };
        for i in 0..x_b.len() {
            let v = self.vcache[i];
            if v > FEAS_TOL && Self::score(v, dse_w[i], dse_on) > self.cutoff {
                self.in_set[i] = true;
                self.set.push(i as u32);
            }
        }
        self.valid = true;
        self.since = 0;
        self.n_build += 1;
    }

    /// `true` when the dirty log forces the structure stale.
    fn stale(&self, xc: &XbClamp) -> bool {
        !self.valid || xc.chg_all
    }

    /// The hyper CHUZR: `Some((ppos, vbest, sbest, Σv))` exactly as the
    /// full scan would choose (`ppos = usize::MAX` ⇒ no violator), or
    /// `None` ⇒ run the full scan (and `build`).
    #[allow(clippy::too_many_arguments)]
    fn choose(
        &mut self,
        x_b: &[f64],
        ub_b: &[f64],
        dse_w: &[f64],
        dse_on: bool,
        basis: &[usize],
        xc: &mut XbClamp,
        list_thr_m: usize,
    ) -> Option<(usize, f64, f64, f64)> {
        if self.stale(xc) {
            self.valid = false;
            return None;
        }
        let cutoff = self.cutoff;
        let reread = |i: usize, me: &mut Self| {
            let v = violation(ub_b[i], x_b[i]);
            let old = me.vcache[i];
            if v != old {
                me.vtot.add(v);
                me.vtot.add(-old);
                me.nviol += usize::from(v > FEAS_TOL);
                me.nviol -= usize::from(old > FEAS_TOL);
                me.vcache[i] = v;
            }
            if v > FEAS_TOL && !me.in_set[i] && Self::score(v, dse_w[i], dse_on) > cutoff {
                me.in_set[i] = true;
                me.set.push(i as u32);
            }
        };
        for k in 0..xc.chg.len() {
            let i = xc.chg[k] as usize;
            xc.chg_mark[i] = false;
            reread(i, self);
        }
        xc.chg.clear();
        if self.prev != usize::MAX {
            reread(self.prev, self); // flip-only exits clamp x_b[p] unlogged
        }
        if list_thr_m > 0 && self.nviol <= list_thr_m {
            // The full scan would switch CHUZR to the infeasibility list
            // here: hand over so the list engages on the same iteration.
            self.valid = false;
            return None;
        }
        let mut ppos = usize::MAX;
        let mut vbest = 0.0f64;
        let mut sbest = 0.0f64;
        let mut k = 0;
        while k < self.set.len() {
            let i = self.set[k] as usize;
            let v = self.vcache[i];
            let sc = Self::score(v, dse_w[i], dse_on);
            if !(v > FEAS_TOL && sc > cutoff) {
                self.in_set[i] = false;
                self.set.swap_remove(k);
                continue;
            }
            // The pivot loop's comparator at tie = 0.
            if ppos == usize::MAX || sc > sbest || (sc >= sbest && basis[i] < basis[ppos]) {
                vbest = v;
                sbest = sc;
                ppos = i;
            }
            k += 1;
        }
        if ppos == usize::MAX && cutoff != f64::NEG_INFINITY {
            self.valid = false; // exhausted: rows below the cutoff may lead
            self.exhausted = true;
            return None;
        }
        self.set_sum += self.set.len();
        self.n_hit += 1;
        self.since += 1;
        if self.set.len() > 8 * self.k {
            self.valid = false; // grew past its purpose: rebuild next visit
        }
        Some((ppos, vbest, sbest, self.vtot.value()))
    }
}

/// `δ -= σθ·β` and `π̂ -= σθ·β` (the edge-travel commit) over β's
/// pattern; `coef = σ·θ` (the old loops' `sigma * theta * beta[i]` parses
/// as `(σ·θ)·β_i`, so hoisting the product is bit-identical, and a zero
/// β_i term is a ±0.0 no-op).
fn dual_edge_step(delta: &mut HVec, piw: &mut [f64], beta: &HVec, coef: f64) {
    match beta.pattern() {
        Some(pat) => {
            for &i in pat {
                let i = i as usize;
                let s = coef * beta[i];
                *delta.entry(i) -= s;
                piw[i] -= s;
            }
        }
        None => {
            let dv = delta.dense_mut();
            for (i, &b) in beta.vals().iter().enumerate() {
                let s = coef * b;
                dv[i] -= s;
                piw[i] -= s;
            }
        }
    }
}

/// `RUSOLVE_DUAL_DEBTSHIFT`: the π half of [`dual_edge_step`] — the
/// debt itself goes to a basic cost shift instead of δ.
fn piw_edge_step(piw: &mut [f64], beta: &HVec, coef: f64) {
    match beta.pattern() {
        Some(pat) => {
            for &i in pat {
                let i = i as usize;
                piw[i] -= coef * beta[i];
            }
        }
        None => {
            for (w, &b) in piw.iter_mut().zip(beta.vals()) {
                *w -= coef * b;
            }
        }
    }
}

/// Flip-batch x_B update `x_b -= b̄` (+ infeasibility-list
/// maintenance) over b̄'s ascending pattern — same per-row operations and
/// the same ascending push order as the all-row loop; rows with b̄_i = 0
/// are no-ops there (`x − ±0.0`, no push).
fn xb_flip_update(
    x_b: &mut [f64],
    bbar: &mut HVec,
    maintain: bool,
    infeas: &mut InfeasList,
    ub_b: &[f64],
    xclamp: &mut XbClamp,
) {
    bbar.sort();
    let mut row = |i: usize| {
        let d = bbar[i];
        x_b[i] -= d;
        if maintain && d != 0.0 && infeas.mark[i] == 0 && violation(ub_b[i], x_b[i]) > FEAS_TOL {
            infeas.push(i);
        }
    };
    match bbar.pattern() {
        Some(pat) => {
            for &i in pat {
                row(i as usize);
            }
        }
        None => {
            for i in 0..bbar.len() {
                row(i);
            }
        }
    }
    xclamp.note(bbar.pattern());
}

/// The pivot's x_B update `x_b += psgn·b̄·step` with the round-off
/// clamp and list maintenance. Pattern mode visits b̄'s ascending pattern
/// and then clamps the [`XbClamp`] rows; the all-row loop runs when b̄ is
/// dense or a wholesale x_b write is pending. Bit-identical to the old
/// all-row loop (see `XbClamp`).
#[allow(clippy::too_many_arguments)]
fn xb_pivot_update(
    x_b: &mut [f64],
    bbar: &mut HVec,
    psgn: f64,
    step: f64,
    maintain: bool,
    infeas: &mut InfeasList,
    ub_b: &[f64],
    xclamp: &mut XbClamp,
) {
    bbar.sort();
    let mut row = |i: usize| {
        let d = psgn * bbar[i] * step;
        x_b[i] += d;
        if x_b[i] < 0.0 && x_b[i] > -FEAS_TOL {
            x_b[i] = 0.0; // round-off clamp (the clamp budget)
        }
        if maintain && d != 0.0 && infeas.mark[i] == 0 && violation(ub_b[i], x_b[i]) > FEAS_TOL {
            infeas.push(i);
        }
    };
    match bbar.pattern() {
        Some(pat) if !xclamp.all => {
            for &i in pat {
                row(i as usize);
            }
            for &i in &xclamp.rows {
                if x_b[i] < 0.0 && x_b[i] > -FEAS_TOL {
                    x_b[i] = 0.0;
                }
            }
        }
        _ => {
            for i in 0..bbar.len() {
                row(i);
            }
        }
    }
    xclamp.note_chg(bbar.pattern());
    xclamp.all = false;
    xclamp.rows.clear();
}

/// Pivotal row `α_v = βᵀ a_v` with its nonzero pattern (the HiGHS `row_ap`
/// HVector: dense array + index list).
struct PriceBuf {
    alpha: Vec<f64>,
    touched: Vec<bool>,
    list: Vec<usize>,
    /// Nonzero pattern of the priced `beta` (the HiGHS `rho` index vector,
    /// `HVector` in `HighsSparseMatrix::priceByRow`): built once per PRICE
    /// by the same pass that decides the technique, then driving the
    /// row-wise copy so hyper-sparse β touches only its pattern's rows
    /// (Hall & McKinnon 2005, hyper-sparse price).
    bpat: Vec<usize>,
    /// The last PRICE ran row-wise, so `list` is the STRUCTURAL
    /// pattern of βᵀA (every column with an entry in a β row), not just the
    /// nonzero results — a column outside it has no entry in any β row.
    structural: bool,
}

impl PriceBuf {
    fn new(total: usize) -> Self {
        PriceBuf {
            alpha: vec![0.0; total],
            touched: vec![false; total],
            list: Vec::new(),
            bpat: Vec::new(),
            structural: false,
        }
    }

    fn clear(&mut self) {
        for &j in &self.list {
            self.alpha[j] = 0.0;
            self.touched[j] = false;
        }
        self.list.clear();
    }

    /// PRICE the row `beta`. Density switch of HiGHS
    /// `HEkk::choosePriceTechnique` (simplex/HEkk.cpp:2830,
    /// `density_for_column_price_switch = 0.75`): a row-wise pass loops only
    /// over the nonzero entries of `beta`; above 75% density the column-wise
    /// pass wins (no scatter, `HighsSparseMatrix::priceByColumn`). Both
    /// accumulate each column's terms in ascending row order from 0.0, so
    /// the values are bit-identical to [`csc_dot`] (zero terms add nothing).
    /// The column pass skips basic/fixed columns; the row pass leaves them
    /// in the pattern (callers filter by `in_basis`).
    ///
    /// Hyper-sparse: the pattern of `beta` is collected in ONE pass
    /// (the technique decision reads its length — the old code paid two
    /// O(m) dense scans per PRICE, a count pass plus the row pass's own
    /// branch over all m rows), and the row-wise copy then walks only the
    /// pattern. The pattern is ascending, so the per-column accumulation
    /// order — and every value — is unchanged.
    // Measured: inlining this loop into the (huge) solve_dual_core
    // makes its cost sensitive to the host function's register allocation;
    // a forced function boundary keeps the hyper-sparse row pass at its
    // out-of-line speed (krka price 25→67 ms when inlined after the
    // row-list code landed).
    #[inline(never)]
    fn price(
        &mut self,
        st: &DualStandard,
        rm: Option<&RowMat>,
        beta: &[f64],
        beta_pat: Option<&[u32]>,
        in_basis: &[bool],
    ) {
        self.clear();
        let m = st.m;
        self.bpat.clear();
        match beta_pat {
            // β's ascending index list (from the HVec BTRAN) replaces
            // the O(m) nonzero scan; same ascending nonzero list.
            Some(pat) => {
                for &i in pat {
                    if beta[i as usize] != 0.0 {
                        self.bpat.push(i as usize);
                    }
                }
            }
            None => {
                for (i, &b) in beta.iter().enumerate() {
                    if b != 0.0 {
                        self.bpat.push(i);
                    }
                }
            }
        }
        self.structural = false;
        match rm {
            Some(rm) if self.bpat.len() as f64 <= 0.75 * m as f64 => {
                self.structural = true;
                for &i in &self.bpat {
                    let b = beta[i];
                    for k in rm.ptr[i]..rm.ptr[i + 1] {
                        let j = rm.col[k] as usize;
                        if !self.touched[j] {
                            self.touched[j] = true;
                            self.list.push(j);
                        }
                        self.alpha[j] += b * rm.val[k];
                    }
                }
            }
            _ => {
                for v in 0..st.ncol() {
                    if in_basis[v] || st.ub[v] == 0.0 {
                        continue;
                    }
                    let a = csc_dot(beta, st, v);
                    if a != 0.0 {
                        self.alpha[v] = a;
                        self.touched[v] = true;
                        self.list.push(v);
                    }
                }
            }
        }
    }
}

/// `d_j += coef * α_j` over the nonzeros of the pivotal row (nonbasic,
/// non-fixed columns, `skip` excluded) with the cost-shift maintenance the
/// old per-pass scan did: a shifted column that fell below ε_s is bumped
/// back to ε_s. HiGHS `HEkkDualRow::updateDual` (HEkkDualRow.cpp:548) with
/// `coef = -θ_d`.
fn update_dj(
    dj: &mut [f64],
    st: &mut DualStandard,
    pbuf: &PriceBuf,
    in_basis: &[bool],
    es: f64,
    coef: f64,
    skip: usize,
) {
    for &v in &pbuf.list {
        if in_basis[v] || st.ub[v] == 0.0 || v == skip {
            continue;
        }
        let mut d = dj[v] + coef * pbuf.alpha[v];
        if st.s[v] > 0.0 && d < es {
            st.s[v] += es - d;
            d = es;
        }
        dj[v] = d;
    }
}

/// The z telescope: the flip-prefix contribution Σ d'_j·Δx̃_j, priced AFTER
/// the event's dual step landed in dj (`update_dj` already ran, so dj[v]
/// IS the post-step effective d; true d = dj − s − pert). `upper` = the
/// candidate SAT at its ũ before the flip (Δx̃ = −ũ); a lower-side
/// candidate flips up (Δx̃ = +ũ). Same slice the flip batch itself walks.
#[inline]
fn z_flip_terms(
    z: &mut f64,
    cands: &[(usize, f64, f64, bool)],
    gq: usize,
    st: &DualStandard,
    dj: &[f64],
) {
    for &(v, _, _, upper) in &cands[..gq] {
        let dt = dj[v] - st.s[v] - st.pert[v];
        *z += if upper { -dt * st.ub[v] } else { dt * st.ub[v] };
    }
}

/// The dual-objective identity at the CURRENT maintained state —
/// caller space, konst ± (πᵀb̃ + Σ_{nonbasic at ũ} d_true·ũ). This is the
/// z seed AND the zcheck reference: it recomputes the identity from the
/// live (piw, dj, at_ub) triple, so comparing it against the incrementally
/// maintained z_val isolates the EVENT BOOKKEEPING (the box-min
/// `caller_dual_bound` at the raw factors π is a different dual point
/// whenever δ ≠ 0 — a legitimate gap the verify-on-crossing absorbs).
fn z_of_state(
    st: &DualStandard,
    in_basis: &[bool],
    at_ub: &[bool],
    dj: &[f64],
    piw: &[f64],
    konst: f64,
    maximize: bool,
) -> f64 {
    let mut z: f64 = piw.iter().zip(&st.b).map(|(p, b)| p * b).sum();
    for v in 0..st.ncol() {
        if !in_basis[v] && at_ub[v] && st.ub[v] != 0.0 {
            z += (dj[v] - st.s[v] - st.pert[v]) * st.ub[v];
        }
    }
    if maximize { konst - z } else { konst + z }
}

/// Reconstruct caller-space x (inverse of the bound transforms; ports the
/// primal `extract_x` logic 1:1).
pub fn extract_x(
    st: &DualStandard,
    lp: &Lp,
    basis: &[usize],
    x_b: &[f64],
    at_ub: &[bool],
) -> Vec<f64> {
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
    let mut x = vec![0.0; lp.num_col];
    #[allow(clippy::needless_range_loop)]
    for j in 0..lp.num_col {
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

/// Primal infeasibility of a basic value: v_i = max(−x_i, x_i − ũ_i, 0).
#[inline]
pub fn violation(ubv: f64, xv: f64) -> f64 {
    (-xv).max(xv - ubv).max(0.0)
}

/// Start primal infeasibility profile: (Σv, count above the
/// 1e-7 feasibility tolerance).
pub fn infeas_profile(st: &DualStandard, basis: &[usize], x_b: &[f64]) -> (f64, usize) {
    let mut sum = 0.0;
    let mut cnt = 0;
    for (i, &v) in basis.iter().enumerate() {
        let viol = violation(st.ub[v], x_b[i]);
        sum += viol;
        if viol > 1e-7 {
            cnt += 1;
        }
    }
    (sum, cnt)
}

/// Export the dual engine's terminal state as a caller-visible [`Basis`]
/// snapshot. The dual standard form pushes structural columns in the
/// same original-variable order and slacks in original non-E row order as the
/// primal form, so [`crate::simplex::assemble_basis`] applies unchanged; the
/// E-row fixed artificials are visible through their rows' entries, and
/// since the dual standard form has no OTHER artificials, a dual export is
/// exactly described by the [`Basis`] contract — its re-import is lossless.
pub(crate) fn export_basis_dual(
    lp: &Lp,
    st: &DualStandard,
    basis: &[usize],
    at_ub: &[bool],
) -> Basis {
    let total = st.ncol();
    let mut in_b = vec![false; total];
    for &v in basis {
        in_b[v] = true;
    }
    let n_slacks = st.kind.iter().filter(|&&k| k == Kind::Slack).count();
    let mut slack_basic = vec![false; n_slacks];
    for k in 0..n_slacks {
        slack_basic[k] = in_b[st.num_structural + k];
    }
    // E-row fixed artificials follow the slacks in ascending E-row order
    // (StdLayout): E-row ordinal eo ↔ standard column num_structural +
    // n_slacks + eo.
    let art_start = st.num_structural + n_slacks;
    let art_basic: Vec<bool> = (0..lp.num_row)
        .filter(|&i| lp.senses[i] == 'E')
        .enumerate()
        .map(|(eo, _)| in_b[art_start + eo])
        .collect();
    crate::simplex::assemble_basis(
        lp,
        &st.cols_of,
        &in_b[..st.num_structural],
        &at_ub[..st.num_structural],
        &slack_basic,
        &art_basic,
    )
}

/// Standard-form variable behind each `Basic` status of a snapshot, in the
/// [`DseWeights`] order (columns ascending, then rows ascending; an E row is
/// its fixed artificial, other rows their slack).
fn basic_vars_in_snapshot_order(lp: &Lp, st: &DualStandard, snap: &Basis) -> Vec<usize> {
    let n_slacks = st.kind.iter().filter(|&&k| k == Kind::Slack).count();
    let art_start = st.num_structural + n_slacks;
    let mut out = Vec::new();
    for j in 0..lp.num_col {
        if snap.col_status[j] == BasisStatus::Basic {
            out.push(st.cols_of[j][0]);
        }
    }
    let (mut k, mut eo) = (0usize, 0usize);
    for i in 0..lp.num_row {
        let v = if lp.senses[i] == 'E' {
            eo += 1;
            art_start + eo - 1
        } else {
            k += 1;
            st.num_structural + k - 1
        };
        if snap.row_status[i] == BasisStatus::Basic {
            out.push(v);
        }
    }
    out
}

/// Carried weights of the snapshot's basic variables, mapped to the slots of
/// the imported `basis` (HiGHS gathers by variable after a re-INVERT,
/// HEkk.cpp:1875–1878). `None` if the counts disagree; slots whose variable
/// the snapshot did not list (artificials filling a shortfall) get 1.
fn import_dse(
    lp: &Lp,
    st: &DualStandard,
    snap: &Basis,
    w: &DseWeights,
    basis: &[usize],
) -> Option<Vec<f64>> {
    if snap.col_status.len() != lp.num_col || snap.row_status.len() != lp.num_row {
        return None;
    }
    let vars = basic_vars_in_snapshot_order(lp, st, snap);
    if vars.len() != w.basic.len() {
        return None;
    }
    let mut by_var = vec![1.0f64; st.ncol()];
    for (&v, &x) in vars.iter().zip(&w.basic) {
        by_var[v] = (x as f64).max(1e-4);
    }
    Some(basis.iter().map(|&v| by_var[v]).collect())
}

/// Weights of the exported snapshot's basic variables (inverse of
/// [`import_dse`]); a basic variable not in `basis` (cannot happen for the
/// dual engine's own export) reads 1.
fn export_dse(
    lp: &Lp,
    st: &DualStandard,
    basis: &[usize],
    dse_w: &[f64],
    snap: &Basis,
) -> DseWeights {
    let mut slot_of = vec![usize::MAX; st.ncol()];
    for (i, &v) in basis.iter().enumerate() {
        slot_of[v] = i;
    }
    let basic = basic_vars_in_snapshot_order(lp, st, snap)
        .into_iter()
        .map(|v| match slot_of[v] {
            usize::MAX => 1.0,
            s => dse_w[s] as f32,
        })
        .collect();
    DseWeights { basic }
}

/// Rebuild the dual standard-form basis + nonbasic-side flags from a caller
/// [`Basis`] (warm entry). Shares the per-variable status interpretation
/// with the primal warm start ([`crate::simplex::warm_plan`]).
///
/// Differences from the primal import: the E-row fixed artificials are the
/// only artificial pool. 'E'-row entries are HONORED: `Basic` re-enters
/// that row's artificial, `AtLower` keeps it out — the old positional fill
/// reconstructed a DIFFERENT basis from an exported snapshot (which artificial
/// was basic was invisible), and the verification then honestly rejected it.
/// A residual shortfall (basic NON-E-row artificials of a primal export,
/// invisible to the snapshot) still fills from the E-artificial pool — the
/// only pool here — and a shortfall that cannot be covered, or an over-full
/// set, rejects the snapshot (cold fallback, per the [`Basis`] counting
/// contract).
fn import_basis_dual(lp: &Lp, st: &DualStandard, snap: &Basis) -> Option<(Vec<usize>, Vec<bool>)> {
    if snap.col_status.len() != lp.num_col || snap.row_status.len() != lp.num_row {
        return None; // row/col churn is not this contract
    }
    let m = st.m;
    let total = st.ncol();
    let mut at_ub = vec![false; total];
    let mut basis: Vec<usize> = Vec::with_capacity(m);
    for j in 0..lp.num_col {
        let plan = crate::simplex::warm_plan(snap.col_status[j], lp.lower[j], lp.upper[j])?;
        match st.cols_of[j].as_slice() {
            [c] => match plan {
                crate::simplex::WarmVar::Basic => basis.push(*c),
                crate::simplex::WarmVar::NonLower => {}
                crate::simplex::WarmVar::NonUpper => at_ub[*c] = true,
            },
            [c1, _c2] => match plan {
                // Lossy free-var handling (see the Basis docs): the x⁺ half
                // re-enters; the dual-feasibility check rejects a wrong
                // reconstruction.
                crate::simplex::WarmVar::Basic => basis.push(*c1),
                crate::simplex::WarmVar::NonLower => {}
                crate::simplex::WarmVar::NonUpper => return None,
            },
            _ => return None,
        }
    }
    // Slacks in ascending non-E row order; 'E' rows carry their fixed
    // artificial's status, artificials in ascending E-row order.
    let n_slacks = st.kind.iter().filter(|&&k| k == Kind::Slack).count();
    let art_start = st.num_structural + n_slacks;
    let mut k = 0usize;
    let mut eo = 0usize;
    for i in 0..lp.num_row {
        if lp.senses[i] == 'E' {
            match snap.row_status[i] {
                BasisStatus::Basic => basis.push(art_start + eo),
                BasisStatus::AtLower => {}
                _ => return None, // fixed artificials have no upper/fixed side
            }
            eo += 1;
        } else {
            match snap.row_status[i] {
                BasisStatus::Basic => basis.push(st.num_structural + k),
                BasisStatus::AtLower => {}
                _ => return None, // slacks have no upper/fixed side
            }
            k += 1;
        }
    }
    let deficit = m.saturating_sub(basis.len());
    if basis.len() + deficit != m {
        return None; // over-full: not a basis
    }
    if deficit > 0 {
        // Residual shortfall = basic NON-E-row artificials of a primal
        // export; fill from the E artificials still nonbasic, in order.
        let mut in_b = vec![false; total];
        for &v in basis.iter() {
            in_b[v] = true;
        }
        let mut filled = 0usize;
        for (v, &was_basic) in in_b.iter().enumerate().skip(art_start) {
            if filled == deficit {
                break;
            }
            if !was_basic {
                basis.push(v);
                filled += 1;
            }
        }
        if filled < deficit {
            return None; // basic non-E artificials beyond the pool: no fit here
        }
    }
    // Sitting at ũ requires a finite ũ (export guarantees it; guard anyway).
    for (v, &at_u) in at_ub.iter().enumerate() {
        if at_u && !st.ub[v].is_finite() {
            return None;
        }
    }
    Some((basis, at_ub))
}

/// Verify the DUAL-feasibility premise of the warm entry: at the basis
/// duals π (effective costs — the invariant the loop maintains, so `s` is
/// included), every nonbasic column must have d̃ ≥ −ptol on the lower side
/// and d̃ ≤ +ptol on the upper side. A violation means the snapshot is not
/// the dual-simplex case (e.g. costs changed since export): the caller
/// falls back honestly rather than running dual simplex on a non-dual-
/// feasible start. Fixed columns (ũ = 0, the E-row artificials among them)
/// are excluded — pricing skips them.
fn warm_dual_feasible(
    st: &DualStandard,
    pfi: &BasisFactorization,
    basis: &[usize],
    at_ub: &[bool],
    config: &SimplexConfig,
) -> bool {
    let m = st.m;
    let total = st.ncol();
    let mut cb = vec![0.0; m];
    for (i, &v) in basis.iter().enumerate() {
        cb[i] = st.cost[v] + st.s[v];
    }
    let mut pi = vec![0.0; m];
    pfi.solve_bt(&cb, &mut pi);
    let mut in_b = vec![false; total];
    for &v in basis {
        in_b[v] = true;
    }
    let tol = config.pricing_tol;
    for v in 0..total {
        if in_b[v] || st.ub[v] == 0.0 {
            continue;
        }
        let d = st.cost[v] + st.s[v] - csc_dot(&pi, st, v);
        if !at_ub[v] {
            if d < -tol {
                return false;
            }
        } else if d > tol {
            return false;
        }
    }
    true
}

/// [`warm_setup`]'s output: (basis, basic values, nonbasic sides, factorization).
type WarmSetup = (Vec<usize>, Vec<f64>, Vec<bool>, BasisFactorization);

/// Build the verified warm-start state of an imported basis (warm entry,
/// shared with [`ReoptSession`] so the two never drift): restore the
/// "basic shifts are 0" invariant (`to_dual_standard` pre-shifts by cost
/// sign alone; a warm-basic shifted column would poison `cb`), factorize the
/// snapshot basis, re-anchor the basic values, and VERIFY the dual-feasibility
/// premise. `None` is the honest rejection the callers fall back cold on:
/// singular basis, or a start that is not the dual-simplex case at all.
fn warm_setup(
    st: &mut DualStandard,
    config: &SimplexConfig,
    wb: Vec<usize>,
    wat: Vec<bool>,
) -> Option<WarmSetup> {
    let m = st.m;
    for &v in wb.iter() {
        st.s[v] = 0.0;
    }
    // Basis straight from the CSC columns (skips the O(m²) dense
    // materialization `DenseMatrix::from_columns` paid per warm entry;
    // `build_sparse` receives exactly what `sparsify` of the dense columns
    // produced — same nonzeros, ascending rows).
    let basis_nz: Vec<Vec<(usize, f64)>> = (0..m).map(|i| st.sparse_col(wb[i])).collect();
    let Ok(mut p) = BasisFactorization::build_sparse(basis_nz, dual_pfi_config(config)) else {
        return None; // singular snapshot basis ⇒ honest fallback
    };
    let mut xb = vec![0.0; m];
    dual_refresh(st, &mut p, &wb, &wat, &mut xb);
    if !warm_dual_feasible(st, &p, &wb, &wat, config) {
        let mut wat = wat;
        if !(config.warm_flip_boxed && warm_flip_boxed(st, &p, &wb, &mut wat, config)) {
            return None; // not the dual-simplex premise ⇒ fallback
        }
        // The flips moved nonbasic values: re-anchor the basic ones.
        dual_refresh(st, &mut p, &wb, &wat, &mut xb);
        return Some((wb, xb, wat, p));
    }
    Some((wb, xb, wat, p))
}

/// Warm-entry repair ([`SimplexConfig::warm_flip_boxed`]): flip every
/// dual-infeasible nonbasic column with a finite upper side `ũ` to the side
/// its reduced cost prefers (HiGHS `correctDualInfeasibilities` flips boxed
/// columns, simplex/HEkkDual.cpp:2388–2450, MIT, ERGO-Code/HiGHS). `false`
/// (nothing changed that matters: the caller rejects) when a violation sits
/// on a column with no finite upper side, which no flip can repair. Same
/// duals and tolerance as [`warm_dual_feasible`], so a `true` return leaves
/// a start that check accepts.
fn warm_flip_boxed(
    st: &DualStandard,
    pfi: &BasisFactorization,
    basis: &[usize],
    at_ub: &mut [bool],
    config: &SimplexConfig,
) -> bool {
    let m = st.m;
    let total = st.ncol();
    let mut cb = vec![0.0; m];
    for (i, &v) in basis.iter().enumerate() {
        cb[i] = st.cost[v] + st.s[v];
    }
    let mut pi = vec![0.0; m];
    pfi.solve_bt(&cb, &mut pi);
    let mut in_b = vec![false; total];
    for &v in basis {
        in_b[v] = true;
    }
    let tol = config.pricing_tol;
    let mut flips = Vec::new();
    for v in 0..total {
        if in_b[v] || st.ub[v] == 0.0 {
            continue;
        }
        let d = st.cost[v] + st.s[v] - csc_dot(&pi, st, v);
        let wrong = if at_ub[v] { d > tol } else { d < -tol };
        if wrong {
            if !st.ub[v].is_finite() {
                return false;
            }
            flips.push(v);
        }
    }
    for v in flips {
        at_ub[v] = !at_ub[v];
    }
    true
}

/// Valid dual bound of the basis duals π in CALLER objective space (the
/// objective cutoff): `π·b + Σ_{nonbasic, d_j < 0} d_j·ũ_j` bounds the shifted-form
/// optimum from below (min sense) because `cᵀx̃ = πᵀb + Σ_j d_j x̃_j` and
/// each nonbasic term is ≥ its best box completion; basic columns are
/// excluded (d_b = 0 exactly at the basis duals). The TRUE costs (no cost-shift
/// shift) are used, so the bound belongs to the caller's LP. Caller-space
/// mapping: `cᵀx = const + Σ_v cost_v x̃_v` for minimize and
/// `cᵀx = const − Σ_v cost_v x̃_v` for maximize, with
/// `const = Σ_j c_j·shift_j` (the bound-shift bookkeeping of
/// `to_dual_standard`). For minimize the result is a LOWER bound on the
/// optimum; for maximize an UPPER bound.
fn caller_dual_bound(lp: &Lp, st: &DualStandard, pi: &[f64], in_basis: &[bool]) -> f64 {
    let mut l: f64 = pi.iter().zip(&st.b).map(|(p, b)| p * b).sum();
    for (v, &basic) in in_basis.iter().enumerate() {
        if basic {
            continue;
        }
        let d = st.cost[v] - csc_dot(pi, st, v);
        if d < 0.0 {
            if !st.ub[v].is_finite() {
                // Unbounded improving side: the bound degenerates to −inf
                // (min) — never prunes. Honest, not clamped.
                l = f64::NEG_INFINITY;
                break;
            }
            l += d * st.ub[v];
        }
    }
    let konst: f64 = lp.objective.iter().zip(&st.shift).map(|(c, s)| c * s).sum();
    if lp.maximize { konst - l } else { konst + l }
}

/// Termination of [`reoptimize`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReoptStatus {
    /// Solved to optimality below the cutoff (or no cutoff).
    Optimal,
    /// The running dual bound proved the LP cannot beat the cutoff: for
    /// minimize `obj_bound ≥ cutoff`, for maximize `obj_bound ≤ cutoff`
    /// (both = the caller-sense dual bound). Fires either as an early stop
    /// mid-loop, or — when the LP solved to an optimum that itself violates
    /// the cutoff — after optimality (exact `obj_bound` = the optimum).
    Pruned,
    /// LP infeasible: prune the node; `obj_bound = +inf` (min) / `−inf` (max).
    Infeasible,
    /// LP unbounded: `obj_bound = −inf` (min) / `+inf` (max).
    Unbounded,
    /// Iteration budget exhausted: `obj_bound` = the last bound evaluated
    /// (NaN if no cutoff was given ⇒ no bound was computed).
    IterationLimit,
}

/// Result of [`reoptimize`] — the MIP search's node-LP answer.
#[derive(Debug, Clone)]
pub struct ReoptResult {
    pub status: ReoptStatus,
    /// Valid bound on the LP optimum in the CALLER's objective space and
    /// sense (min: lower bound; max: upper bound). Exact at Optimal, and at
    /// Pruned when the prune was decided post-optimality.
    pub obj_bound: f64,
    /// Dual-engine outer-loop visits consumed (the node-LP pivot cost the
    /// MIP search budgets around; cold fallbacks report their own solve's
    /// phase-1 + phase-2 total).
    pub iters: usize,
    /// `false` ⇒ the snapshot was rejected (contract violation, singular
    /// refactorization, dual-INFeasible start — the bound-change premise
    /// broken) or the dual engine hit distress: the result comes from an
    /// honest COLD solve, not from a warm reopt.
    pub warm: bool,
    /// Optimal primal solution (empty on every other status).
    pub x: Vec<f64>,
    /// Exact caller objective at Optimal (NaN otherwise).
    pub objective: f64,
    /// Final basis for the grandchildren (Optimal only; the dual engine
    /// exports its terminal basis, see [`export_basis_dual`]).
    pub basis: Option<Basis>,
    /// Row duals π at an Optimal verdict (length `num_row`, caller row
    /// order/senses), in the SAME convention as
    /// [`crate::simplex::SimplexSolution::row_duals`]: the pair
    /// `(π, d)` is defined by the single identity `d_j = c_j − πᵀa_j` in
    /// the caller's own space. Empty unless [`Self::duals_valid`].
    pub row_duals: Vec<f64>,
    /// Reduced costs d at an Optimal verdict (length `num_col`, caller
    /// column order), same convention as
    /// [`crate::simplex::SimplexSolution::reduced_costs`]. Empty unless
    /// [`Self::duals_valid`].
    pub reduced_costs: Vec<f64>,
    /// True iff [`Self::row_duals`]/[`Self::reduced_costs`] describe an
    /// Optimal solution of the CALLER's LP. A warm Optimal always sets it;
    /// the warm-to-cold fallback sets it exactly when that cold solve
    /// produced valid duals (presolve stays off on this contract, so an
    /// Optimal fallback does carry them); every non-Optimal status leaves
    /// it false. `false` means "no information", never "zero" — MIP may
    /// drop its per-node basis refactorization (`basis_row_duals`) exactly
    /// when this is true.
    pub duals_valid: bool,
    /// `Infeasible` certificate: Farkas vector y over
    /// the CALLER's rows (length `num_row`, caller row order/senses) with
    /// `yᵀ rhs > sup { yᵀ(A x + s) : bounds, row-slack ranges }` — the same
    /// identity as [`crate::simplex::SimplexSolution::dual_ray`], verified
    /// numerically against the ORIGINAL rows before it is returned; an
    /// unverified candidate is reported as `None`, never as a wrong ray.
    /// Populated at `Infeasible` from both the warm dual loop's Farkas row
    /// and the cold fallback's phase-1 exit; `None` on every other status.
    pub dual_ray: Option<Vec<f64>>,
    /// At an early `Pruned` stop, the row multipliers π (length
    /// `num_row`, caller row order and sense, the convention of
    /// [`Self::row_duals`]) whose Lagrangian bound
    /// `π·b + Σ_{d_j<0} d_j ũ_j` (the one `obj_bound` reports) crossed the
    /// cutoff. NOT necessarily sign-correct per row sense (the point is only
    /// dual feasible up to the running shifts): a consumer building a proof
    /// clips wrong-signed entries to 0, as `LpDuals::from_rows` does. Empty
    /// on every other status (a post-hoc `Pruned` carries `row_duals`).
    pub prune_duals: Vec<f64>,
    /// Final dual steepest-edge weights of an `Optimal` warm verdict when
    /// [`reoptimize_dse`] ran with `DseMode::Fresh` or `Carry`; `None`
    /// otherwise (the cold fallback never has them).
    pub dse_weights: Option<DseWeights>,
}

/// Warm-reopt census (`RUSOLVE_WARM_STATS=1`): thread-local accumulation
/// of where warm reoptimization time goes across a whole MIP run. The
/// per-solve `RUSOLVE_DUAL_STATS` lines only print at the Optimal exit, so
/// pruned nodes — the majority of a MIP tree — are invisible there; this
/// census counts every call and prints one line per 512. Timers are in ns
/// accumulated; `iters` sums pivot iterations so per-iteration furniture can
/// be priced. Diagnostic lever only: no behavior, off by default.
#[derive(Default)]
struct WarmCensus {
    calls: usize,
    warm_ok: usize,
    fallback: usize,
    pruned: usize,
    optimal: usize,
    infeas: usize,
    iters: u64,
    /// `to_dual_standard` (the single build — `reoptimize_dse`'s copy is
    /// passed through; the earlier duplicate second build is gone).
    t_std: f64,
    /// `import_basis_dual`.
    t_import: f64,
    /// `warm_setup` (LU factorization + `dual_refresh` + premise check).
    t_setup: f64,
    /// `RowMat::from_dual` (the row-wise PRICE copy, nnz ≥ 20k only).
    t_rowmat: f64,
    /// The whole `solve_dual_core` warm call.
    t_loop: f64,
    /// Per-iteration fresh-π BTRAN of the basic costs (`need_pi`).
    t_pi_bt: f64,
    /// Per-iteration `caller_dual_bound` (the O(nnz) cutoff scan).
    t_bound: f64,
}

std::thread_local! {
    static WARM_CENSUS: std::cell::RefCell<WarmCensus> = const {
        std::cell::RefCell::new(WarmCensus {
            calls: 0,
            warm_ok: 0,
            fallback: 0,
            pruned: 0,
            optimal: 0,
            infeas: 0,
            iters: 0,
            t_std: 0.0,
            t_import: 0.0,
            t_setup: 0.0,
            t_rowmat: 0.0,
            t_loop: 0.0,
            t_pi_bt: 0.0,
            t_bound: 0.0,
        })
    };
}

/// Record `nanos` into a census field (caller checks the lever first).
fn census_acc(f: impl FnOnce(&mut WarmCensus)) {
    WARM_CENSUS.with_borrow_mut(f);
}

/// One census line per 512 recorded calls.
fn census_tick() {
    const EVERY: usize = 512;
    WARM_CENSUS.with_borrow(|c| {
        if c.calls % EVERY != 0 {
            return;
        }
        let ms = |ns: f64| ns / 1e6;
        eprintln!(
            "rusolve warm census: calls={} warm={} fb={} pruned={} opt={} infeas={} \
             iters={} | ms: std={:.0} import={:.0} setup={:.0} rowmat={:.0} \
             loop={:.0} pibt={:.0} bound={:.0}",
            c.calls,
            c.warm_ok,
            c.fallback,
            c.pruned,
            c.optimal,
            c.infeas,
            c.iters,
            ms(c.t_std),
            ms(c.t_import),
            ms(c.t_setup),
            ms(c.t_rowmat),
            ms(c.t_loop),
            ms(c.t_pi_bt),
            ms(c.t_bound),
        );
    });
}

/// Warm DUAL reoptimization with an objective cutoff — the MIP node
/// case: a child LP differs from the parent by bound changes that keep the
/// parent basis DUAL feasible (the premise, verified at entry) while making
/// it primal infeasible; dual simplex then reoptimizes in a few pivots and
/// stops as soon as the running dual bound proves the node cannot beat the
/// incumbent.
///
/// * `cutoff`: caller-sense objective cutoff — prune when the dual bound
///   reaches it: `bound ≥ cutoff` for minimize, `bound ≤ cutoff` for
///   maximize. `None` = solve to optimality.
/// * `iter_limit`: `None` = `config.max_iterations`.
/// * Presolve is NOT applied (MIP node LPs run presolve-off; the snapshot
///   contract lives on the caller's LP). A rejected snapshot or dual-engine
///   distress falls back to a cold primal solve of the SAME LP (`warm =
///   false`, verdicts mapped into [`ReoptStatus`], cutoff honored
///   post-hoc: an optimum beyond the cutoff reports `Pruned` with the exact
///   objective as `obj_bound`).
pub fn reoptimize(
    lp: &Lp,
    config: &SimplexConfig,
    basis: &Basis,
    cutoff: Option<f64>,
    iter_limit: Option<usize>,
) -> Result<ReoptResult, SimplexError> {
    reoptimize_dse(lp, config, basis, cutoff, iter_limit, DseMode::Off)
}

/// Dual steepest-edge weights `‖e_iᵀB⁻¹‖²` of the BASIC variables of a
/// [`Basis`] snapshot, so a warm re-solve can keep pricing with them instead
/// of restarting at 1 (task 2; HiGHS keeps `dual_edge_weight_` with the
/// basis it belongs to: HEkk::putIterate/getIterate, HEkk.cpp:3812–3836, and
/// scatters/gathers them by variable across a re-INVERT, HEkk.cpp:1799–1803,
/// :1875–1878). Layout: one entry per `Basic` status of the snapshot, in
/// ascending caller order (columns first, then rows; an E row's entry is its
/// fixed artificial). `f32` keeps a node snapshot small.
#[derive(Debug, Clone, PartialEq)]
pub struct DseWeights {
    pub basic: Vec<f32>,
}

/// How [`reoptimize_dse`] prices the leaving row.
#[derive(Debug, Clone, Copy)]
pub enum DseMode<'a> {
    /// Dantzig (largest infeasibility), the historic warm behaviour.
    Off,
    /// Steepest edge from unit weights (exact only at a slack basis).
    Fresh,
    /// Steepest edge from the carried weights of the offered basis.
    Carry(&'a DseWeights),
}

/// [`reoptimize`] with a leaving-row pricing mode. At an `Optimal` verdict
/// with `dse != Off` the result carries the final weights for the children.
pub fn reoptimize_dse(
    lp: &Lp,
    config: &SimplexConfig,
    basis: &Basis,
    cutoff: Option<f64>,
    iter_limit: Option<usize>,
    dse: DseMode<'_>,
) -> Result<ReoptResult, SimplexError> {
    reoptimize_dse_cached(
        lp,
        config,
        basis,
        cutoff,
        iter_limit,
        dse,
        &mut StdCache::new(),
        0,
    )
}

/// A DualStandard kept alive across the warm reoptimizes of one
/// `Lp`. The MIP engine's node loop calls [`reoptimize_dse_cached`] many
/// times per second on the SAME template with only bound VALUES edited
/// between calls (branching/probing mutate `lower`/`upper` in place;
/// measured on our benchmark set: the rebuild was 34 % of
/// tracked warm time on mik-class instances — `cols_of` alone is one heap
/// Vec per column). The cache skips that rebuild whenever the caller's
/// `gen` promise holds:
///
/// **Contract:** between two calls with the SAME `gen`, only
/// `lp.lower`/`lp.upper` VALUES changed. ANY other edit — matrix rows
/// added/removed (`add_cuts`/`retain_cuts`), rhs, senses, objective,
/// dimensions — must bump `gen`. A same-`gen` call with changed structure
/// is undefined (stale slots ⇒ wrong answers, and no check can catch a
/// same-shape rhs/objective rewrite cheaply).
///
/// Bounds edits are absorbed incrementally: a column whose bound CATEGORY
/// (which of lower/upper is finite — it fixes the slot count and the flip
/// mapping) changed forces one full rebuild; otherwise `b` is updated by
/// the shift delta over just that column's nonzeros, the box width ũ is
/// rewritten, and the solve state (`s`, `pert`, `at_ub`) is reset to
/// exactly what a fresh `to_dual_standard` would carry.
pub struct StdCache {
    generation: u64,
    st: Option<DualStandard>,
}

impl Default for StdCache {
    fn default() -> Self {
        Self::new()
    }
}

impl StdCache {
    pub fn new() -> Self {
        Self {
            generation: 0,
            st: None,
        }
    }
}

/// [`reoptimize_dse`] on a reusable [`StdCache`] (see its contract). A
/// cache miss (first call, gen bump, dimension mismatch) pays the ordinary
/// `to_dual_standard` rebuild and adopts the result.
#[allow(clippy::too_many_arguments)] // the plain signature + cache + gen (precedent: solve_dual_core)
pub fn reoptimize_dse_cached(
    lp: &Lp,
    config: &SimplexConfig,
    basis: &Basis,
    cutoff: Option<f64>,
    iter_limit: Option<usize>,
    dse: DseMode<'_>,
    cache: &mut StdCache,
    generation: u64,
) -> Result<ReoptResult, SimplexError> {
    reoptimize_dse_cached_until(
        lp, config, basis, cutoff, iter_limit, dse, cache, generation, None,
    )
}

/// [`reoptimize_dse_cached`] under a wall-clock `deadline` (the MIP search's
/// time limit). Past it the warm dual stops as at its iteration limit
/// (`ReoptStatus::IterationLimit`, a valid dual bound) and a cold fallback
/// gets the same deadline (`SimplexSolver::set_deadline`). Before this,
/// warm node/cut-loop LPs ran without a
/// clock and carried co-100, fast0507 and bnatt400 4–16 s past a 60 s
/// limit. `None` = the clock-free behaviour above, bit for bit.
#[allow(clippy::too_many_arguments)]
pub fn reoptimize_dse_cached_until(
    lp: &Lp,
    config: &SimplexConfig,
    basis: &Basis,
    cutoff: Option<f64>,
    iter_limit: Option<usize>,
    dse: DseMode<'_>,
    cache: &mut StdCache,
    generation: u64,
    deadline: Option<std::time::Instant>,
) -> Result<ReoptResult, SimplexError> {
    let warm_stats = std::env::var("RUSOLVE_WARM_STATS").is_ok();
    let mark = std::time::Instant::now();
    // Same gen + same shape ⇒ bounds-only edits: refresh in place. The
    // refresh declines (category change / invalid bounds) into a full
    // rebuild through the validating constructor, so every rebuild error
    // path is to_dual_standard's own.
    let reusable = cache.generation == generation
        && cache.st.is_some()
        && cache
            .st
            .as_ref()
            .is_some_and(|st| st.m == lp.num_row && st.cols_of.len() == lp.num_col);
    let st: &mut DualStandard = if reusable {
        let st = cache.st.as_mut().expect("checked some");
        if !refresh_dual_standard(lp, config, st) {
            *st = to_dual_standard(lp, config)?;
        }
        st
    } else {
        // Drop any old pair first: a build error must leave the cache
        // empty (never stale structure under a still-matching gen).
        cache.st = None;
        let built = to_dual_standard(lp, config)?;
        cache.st = Some(built);
        cache.generation = generation;
        cache.st.as_mut().expect("just stored")
    };
    let t_std = mark.elapsed();
    let mark = std::time::Instant::now();
    let warm = import_basis_dual(lp, st, basis);
    let t_import = mark.elapsed();
    if warm_stats {
        census_acc(|c| {
            c.calls += 1;
            c.t_std += t_std.as_nanos() as f64;
            c.t_import += t_import.as_nanos() as f64;
        });
    }
    let Some((wb, wat)) = warm else {
        if warm_stats {
            census_acc(|c| c.fallback += 1);
            census_tick();
        }
        return cold_reopt(lp, config, cutoff, deadline);
    };
    let mut last_bound: Option<f64> = None;
    let mut dse_out: Option<DseWeights> = None;
    let core_mark = std::time::Instant::now();
    let core = solve_dual_core(
        lp,
        st,
        config,
        Some((wb, wat)),
        cutoff,
        iter_limit,
        &mut last_bound,
        match dse {
            DseMode::Carry(w) => Some((Some((basis, w)), &mut dse_out)),
            DseMode::Fresh => Some((None, &mut dse_out)),
            DseMode::Off => None,
        },
        // Warm MIP re-optimizations keep their own caps; the clock is the
        // caller's (the search's deadline, else none).
        deadline,
    );
    let t_loop = core_mark.elapsed();
    if warm_stats {
        census_acc(|c| c.t_loop += t_loop.as_nanos() as f64);
        census_tick();
    }
    match core? {
        None => {
            if warm_stats {
                census_acc(|c| {
                    c.fallback += 1;
                    c.warm_ok += 1;
                });
                census_tick();
            }
            cold_reopt(lp, config, cutoff, deadline)
        }
        Some(DualExit::Pruned { bound, iters, pi }) => {
            if warm_stats {
                census_acc(|c| {
                    c.warm_ok += 1;
                    c.pruned += 1;
                    c.iters += iters as u64;
                });
                census_tick();
            }
            Ok(ReoptResult {
                status: ReoptStatus::Pruned,
                obj_bound: bound,
                iters,
                warm: true,
                x: vec![],
                objective: f64::NAN,
                basis: None,
                // An early prune has no optimal basis to take duals from.
                row_duals: vec![],
                reduced_costs: vec![],
                duals_valid: false,
                dual_ray: None,
                // The π that proved the cutoff, in caller row space/sense.
                prune_duals: {
                    let sign = if lp.maximize { -1.0 } else { 1.0 };
                    pi.iter().map(|&p| sign * p).collect()
                },
                dse_weights: None,
            })
        }
        Some(DualExit::Sol(sol)) => {
            if warm_stats {
                census_acc(|c| {
                    c.warm_ok += 1;
                    c.iters += (sol.phase1_iterations + sol.phase2_iterations) as u64;
                    match sol.status {
                        SimplexStatus::Optimal => c.optimal += 1,
                        SimplexStatus::Infeasible => c.infeas += 1,
                        _ => {}
                    }
                });
                census_tick();
            }
            let mut r = map_reopt_status(lp, &sol, cutoff);
            if sol.status == SimplexStatus::IterationLimit {
                r.obj_bound = last_bound.unwrap_or(f64::NAN);
            }
            r.iters = sol.phase1_iterations + sol.phase2_iterations;
            r.warm = true;
            r.dse_weights = dse_out;
            Ok(r)
        }
    }
}

/// Map a [`SimplexSolution`] verdict into the result shape (shared by
/// the warm and the cold-fallback paths).
fn map_reopt_status(lp: &Lp, sol: &SimplexSolution, cutoff: Option<f64>) -> ReoptResult {
    let beyond = cutoff.is_some_and(|c| {
        if lp.maximize {
            sol.objective <= c
        } else {
            sol.objective >= c
        }
    });
    let (status, obj_bound, objective) = match sol.status {
        SimplexStatus::Optimal => {
            if beyond {
                (ReoptStatus::Pruned, sol.objective, sol.objective)
            } else {
                (ReoptStatus::Optimal, sol.objective, sol.objective)
            }
        }
        SimplexStatus::Infeasible => (
            ReoptStatus::Infeasible,
            if lp.maximize {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            },
            f64::NAN,
        ),
        SimplexStatus::Unbounded => (
            ReoptStatus::Unbounded,
            if lp.maximize {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            },
            f64::NAN,
        ),
        SimplexStatus::IterationLimit => (ReoptStatus::IterationLimit, f64::NAN, f64::NAN),
    };
    ReoptResult {
        status,
        obj_bound,
        iters: sol.phase1_iterations + sol.phase2_iterations,
        warm: false,
        x: sol.x.clone(),
        objective,
        basis: sol.basis.clone(),
        // The verdict's dual/ray evidence rides along on both the
        // warm and the cold-fallback path — the presolve-OFF solvers fill
        // it (Optimal ⇒ duals, Infeasible ⇒ verified Farkas ray), the
        // presolve-ON surface never does.
        row_duals: sol.row_duals.clone(),
        reduced_costs: sol.reduced_costs.clone(),
        duals_valid: sol.duals_valid,
        dual_ray: sol.dual_ray.clone(),
        // A post-hoc prune solved to optimality: its duals are the proof
        // multipliers too.
        prune_duals: if status == ReoptStatus::Pruned && sol.duals_valid {
            sol.row_duals.clone()
        } else {
            Vec::new()
        },
        dse_weights: None,
    }
}

/// Honest cold fallback for [`reoptimize`]: cold primal solve of the same
/// LP (presolve stays OFF — the reoptimize contract is snapshot-shaped and
/// presolve would break it; the config's dual lever is kept, its own
/// distress fallback then runs the primal solver).
fn cold_reopt(
    lp: &Lp,
    config: &SimplexConfig,
    cutoff: Option<f64>,
    deadline: Option<std::time::Instant>,
) -> Result<ReoptResult, SimplexError> {
    let cfg = SimplexConfig {
        presolve: false,
        ..*config
    };
    let mut solver = SimplexSolver::new(cfg);
    solver.set_deadline(deadline);
    let sol = solver.solve(lp)?;
    Ok(map_reopt_status(lp, &sol, cutoff))
}

/// Primal-feasibility tolerance (our tolerance convention): the existing
/// round-off clamp budget in simplex.rs (`x_b ∈ (−1e-7,0) → 0`) — a leave or
/// flip recorded at a bound must lie within the clamp.
const FEAS_TOL: f64 = 1e-7;

/// Standard-form nonzeros below which the row-wise PRICE and incremental d
/// stay off: the per-pass saving is negligible there and the changed
/// rounding only moves ties (MIP root LPs, tiny sweeps stay on the
/// validated fresh-pi path, bit for bit).
const DUAL_FAST_MIN_NNZ: usize = 20_000;

/// Dual-engine entry point (lever dispatch in `SimplexSolver::solve`).
/// `Ok(Some(sol))` ⇒ dual verdict; `Ok(None)` ⇒ NumericalDistress-style
/// fallback: the caller re-runs the primal solver on the same LP (the
/// primal stays ground truth for all A/Bs).
///
/// Classical dual simplex: Dantzig leaving (max violation),
/// event-loop dual ratio test with bound flips, landing check, pivot + PFI
/// retry, optimality check (primal-step arm + one-shift cleanup),
/// Infeasible/Unbounded with one verification pass each. No DuPh1, no
/// perturbation.
/// What [`solve_dual_core`] terminates on. `Sol` = a verdict over the
/// caller's LP; `Pruned` = the running dual bound crossed the objective
/// cutoff without solving to optimality — only reachable with
/// `cutoff == Some`. `Ok(None)` (not an arm) keeps the NumericalDistress /
/// warm-rejection signal: the caller falls back.
enum DualExit {
    Sol(SimplexSolution),
    Pruned {
        bound: f64,
        iters: usize,
        /// The standard-form π whose [`caller_dual_bound`] crossed the cutoff,
        /// before the caller-sense mapping.
        pi: Vec<f64>,
    },
}

/// Cold dual solve. `deadline`: solve-internal wall-clock bound,
/// checked every iteration exactly like the primal engine's deadline check
/// (`SimplexSolver::set_deadline`; HiGHS `HEkk::bailout`, HEkk.cpp
/// 3441–3461, MIT — mechanism, not code); a hit returns the
/// IterationLimit-status point. `None` = no clock read at all.
pub fn solve_dual_raw(
    lp: &Lp,
    config: &SimplexConfig,
    deadline: Option<std::time::Instant>,
) -> Result<Option<SimplexSolution>, SimplexError> {
    let mut st = to_dual_standard(lp, config)?;
    solve_dual_prepared(lp, &mut st, config, deadline)
}

/// The `solve_dual_raw` body over a caller-prepared [`DualStandard`]:
/// the cold-start routing detector (`simplex.rs` `route_cold_phase1_dual`)
/// builds the standard form to test the slack-basis start, then hands it
/// here so the routed solve does not pay a second O(nnz) build.
pub fn solve_dual_prepared(
    lp: &Lp,
    st: &mut DualStandard,
    config: &SimplexConfig,
    deadline: Option<std::time::Instant>,
) -> Result<Option<SimplexSolution>, SimplexError> {
    match solve_dual_core(lp, st, config, None, None, None, &mut None, None, deadline)? {
        None => Ok(None),
        Some(DualExit::Sol(sol)) => Ok(Some(sol)),
        // cutoff == None ⇒ the check never runs.
        Some(DualExit::Pruned { .. }) => unreachable!("prune without cutoff"),
    }
}

/// Farkas ray from the dual engine's empty-event-set Infeasible exit.
/// Row p is basic-violated (σ = +1 below 0, −1 above its upper) and no
/// column can move it back: `β = B⁻ᵀe_p` then certifies infeasibility of
/// the standard form as `y = −σ·β` — `yᵀb = −σ·x_b[p] > 0` by violation,
/// and the empty event set is exactly the per-column sign/bound conditions
/// [`farkas_ok`] checks. The dual standard form never negates rows and the
/// objective sign lives entirely in `cost`, so `y` maps to caller space by
/// the identity (contrast `attempt_farkas`, which undoes row negations).
/// Same ray discipline as the primal path: the debug build ASSERTS the
/// verification (a derivation bug must be loud in tests), the release
/// build reports the honest `None` — an unverified ray is never returned.
fn dual_farkas(lp: &Lp, sigma: f64, beta: &[f64]) -> Option<Vec<f64>> {
    let y: Vec<f64> = beta.iter().map(|&v| -sigma * v).collect();
    debug_assert!(
        farkas_ok(lp, &y),
        "empty-event-set row must verify as a Farkas certificate: {y:?}"
    );
    farkas_ok(lp, &y).then_some(y)
}

/// Warm-DSE plumbing: (basis + weights to import, slot receiving the exported weights).
type WarmDseIo<'a> = (
    Option<(&'a Basis, &'a DseWeights)>,
    &'a mut Option<DseWeights>,
);

/// Re-anchor schedule for incremental-π transport (`RUSOLVE_DUAL_INCP`,
/// see the lever block in [`solve_dual_core`]). Split out from the env and
/// the pivot loop so the schedule is unit-testable without process env
/// (the `parse_ft_flag` precedent). Rules:
///
/// * an EVENT (refactorization, audit step, perturbation, hand-over, shift
///   release — anywhere the maintained representation is invalidated today)
///   forces a re-anchor on the next pass;
/// * on the fresh-d̃ path (`dj_valid_inc == false`: no incremental d̃, or it
///   just went invalid) the cadence forces one every `cad` PIVOTS — the
///   rank-one π̂ steps are where the rounding debt accrues;
/// * on the incremental-d̃ path with d̃ still VALID, the cadence never
///   forces a re-anchor (only events do): that path pays no π-BTRAN today
///   and a cadence anchor would only add a computeDual re-derivation.
#[derive(Clone, Copy, Debug)]
struct PiAnchor {
    pending: bool,
    since: usize,
    cad: usize,
}

impl PiAnchor {
    /// Starts armed: `piw` begins zeroed, so the first pass must anchor.
    fn new(cad: usize) -> Self {
        Self {
            pending: true,
            since: 0,
            cad,
        }
    }
    /// Is a re-anchor due this pass? `dj_valid_inc` = the incremental d̃ is
    /// live (the path skipping the π-BTRAN by construction).
    fn due(&self, dj_valid_inc: bool) -> bool {
        self.pending || (!dj_valid_inc && self.since >= self.cad)
    }
    /// A re-anchor (π-BTRAN) is happening this pass: disarm and restart the
    /// cadence clock.
    fn anchor(&mut self) {
        self.pending = false;
        self.since = 0;
    }
    /// A pivot landed: one rank-one π̂ step of accumulated rounding debt.
    fn pivot(&mut self) {
        self.since += 1;
    }
    /// Representation invalidated: force the next pass to anchor.
    fn event(&mut self) {
        self.pending = true;
    }
}

/// The last basis that factorized (HiGHS `backtracking_basis_`,
/// HEkk::putBacktrackingBasis/getBacktrackingBasis, MIT) and the state
/// that makes it a consistent dual iterate again.
#[derive(Default)]
struct Backtrack {
    on: bool,
    valid: bool,
    count: usize,
    max: usize,
    epoch: u64,
    eta_limit: usize,
    basis: Vec<usize>,
    at_ub: Vec<bool>,
    s: Vec<f64>,
    dse_w: Vec<f64>,
    delta: Vec<f64>,
}

/// The dual pivot loop over the shifted-bounds form, cold (warm-start
/// parameters all `None` ⇒ bit-identical to a plain cold `solve_dual_raw`) or warm:
///
/// * `warm == Some((basis, at_ub))` — start from a verified import
///   (`import_basis_dual` did the shape/count checks; this function still
///   refactorizes, recomputes x_b, and VERIFIES dual feasibility, returning
///   `Ok(None)` — an honest rejection — on any failure).
/// * `cutoff == Some(cut)` — after each iteration's BTRAN, evaluate the
///   valid dual bound of the maintained basis duals in CALLER objective
///   space and stop as soon as it proves the LP cannot beat the cutoff:
///   `bound ≥ cut` for minimize, `bound ≤ cut` for maximize (the bound for
///   max is an UPPER bound on the caller objective). `last_bound` records
///   the most recent evaluation for the iteration-limit report.
/// * `iter_limit` replaces `config.max_iterations` when given.
#[allow(clippy::too_many_arguments)]
fn solve_dual_core(
    lp: &Lp,
    // The caller's standard form, borrowed: `reoptimize_dse` shares its
    // import copy instead of rebuilding (formerly a second full
    // `to_dual_standard` per warm node LP), and the MIP
    // engine hands in a CACHED DualStandard that survives across node LPs
    // (bounds-only edits refresh it in place — see `StdCache`).
    st: &mut DualStandard,
    config: &SimplexConfig,
    warm: Option<(Vec<usize>, Vec<bool>)>,
    cutoff: Option<f64>,
    iter_limit: Option<usize>,
    last_bound: &mut Option<f64>,
    // Warm steepest edge (task 2): `Some((carried, out))` turns DSE on for a
    // warm solve, `carried` = the offered snapshot and its weights (None =
    // unit weights) and `out` receives the final weights at Optimal.
    mut warm_dse: Option<WarmDseIo<'_>>,
    // Solve-internal wall-clock deadline (see `solve_dual_raw`).
    deadline: Option<std::time::Instant>,
) -> Result<Option<DualExit>, SimplexError> {
    let warm_stats = warm.is_some() && std::env::var("RUSOLVE_WARM_STATS").is_ok();
    let st = &mut *st;
    let m = st.m;
    let total = st.ncol();
    let debug = std::env::var("RUSOLVE_DEBUG").is_ok();
    let dtrace = std::env::var("RUSOLVE_DTRACE").is_ok();
    let autopsy_every = std::env::var("RUSOLVE_AUTOPSY_EVERY")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&v| v > 0)
        .unwrap_or(64);
    let dprog = std::env::var("RUSOLVE_DPROG").is_ok();
    let t0 = std::time::Instant::now();
    let cold = warm.is_none();
    // Cold solves keep a cost shift on a column after it enters the basis
    // (HiGHS shifts persist while basic and are removed by cleanup): the
    // old release-on-entry changed the basic cost, hence pi, and planted
    // dual infeasibilities all over the nonbasics (e226/scagr25/share1b
    // bails with perturbation + DSE). The one-shot cleanup drops them all.
    // Lever: RUSOLVE_DUAL_KEEPSHIFT=0 restores release-on-entry.
    let keep_shift = cold
        && std::env::var("RUSOLVE_DUAL_KEEPSHIFT")
            .map(|v| v != "0")
            .unwrap_or(true);
    let (mut basis, mut x_b, mut at_ub, mut pfi) = match warm {
        None => {
            let (b, diag) = initial_basis(lp, st);
            let xb = initial_x_b(st, &b, &diag);
            let at = st.at_ub.clone();
            // Initial basis B = ±I: column i of the basis
            // matrix is diag[i]·e_i — one nonzero per column, so the sparse
            // build skips the O(m²) dense ±I materialization (800 MB at
            // m = 10k) the cold start used to pay.
            let unit: Vec<Vec<(usize, f64)>> = (0..m).map(|i| vec![(i, diag[i])]).collect();
            let p = BasisFactorization::build_sparse(unit, dual_pfi_config(config))?;
            (b, xb, at, p)
        }
        Some((wb, wat)) => {
            // Warm entry (same [`warm_setup`] the session uses).
            let setup_mark = std::time::Instant::now();
            let setup = warm_setup(st, config, wb, wat);
            let t_setup = setup_mark.elapsed();
            if warm_stats {
                census_acc(|c| c.t_setup += t_setup.as_nanos() as f64);
            }
            match setup {
                Some(setup) => setup,
                None => return Ok(None), // singular / dual-infeasible ⇒ fallback
            }
        }
    };
    // Cold solves perturb up front (HiGHS); RUSOLVE_DUAL_PERT=0 = lever off.
    if cold
        && std::env::var("RUSOLVE_DUAL_PERT")
            .map(|v| v != "0")
            .unwrap_or(true)
    {
        highs_cost_perturbation(st);
    }
    if debug {
        let (sumv, cnt) = infeas_profile(st, &basis, &x_b);
        let shifted = st.s.iter().filter(|&&sv| sv > 0.0).count();
        let e_rows = st.kind.iter().filter(|&&k| k == Kind::Artificial).count();
        eprintln!(
            "rusolve dual: start m={m} structural={} e_rows={e_rows} sumv={sumv:.6e} \
             n_infeas={cnt} shifted={shifted}",
            st.num_structural,
        );
    }

    let ptol = config.pricing_tol;
    let pivtol = config.pivot_tol;
    let es = eps_s(config);

    let mut in_basis = vec![false; total];
    // The MAINTAINED dual point is π + δ: π is recomputed from the basis
    // factors each iteration and δ is the accumulated edge travel paid for
    // by flip/saturation batches that did not pivot (δ-commit).
    let mut piw = vec![0.0; m];
    // Hyper-sparse passes: δ, β, b̄, the flip rhs, a_q, e_p
    // and τ are index-carrying vectors (HiGHS `HVector`; Hall & McKinnon
    // 2005): every per-iteration pass over them — solve entries/results,
    // the x_B / DSE / δ updates and their clears — costs O(pattern) on a
    // hyper-sparse iteration instead of O(m). Same values as the dense
    // loops they replace (pinned by `pfi::tests::hv_stack_bit_exact`).
    let mut delta = HVec::new(m);
    let mut beta = HVec::new(m);
    // CHUZR row cache: ub[basis[i]] per row position. The leaving-row scan
    // (and the infeas-list maintenance) otherwise gather the bound through
    // the basis permutation — random access over the whole column space per
    // row; cached, the scan is two streamed arrays. Only the pivot row's
    // entry changes per pivot (O(1)); anything that moves the basis outside
    // the pivot loop (the optimality audit) re-arms a wholesale rebuild.
    let mut ub_b: Vec<f64> = basis.iter().map(|&v| st.ub[v]).collect();
    let mut ub_b_dirty = false;
    let mut cb = vec![0.0; m];
    let mut bbar = HVec::new(m);
    let mut ftmp = HVec::new(m);
    // Raw entering column a_q, gathered from the CSC at every pivot/flip
    // (the primal engine's `qcol` precedent — no per-iteration allocation).
    let mut qcol = HVec::new(m);
    let mut e_p = HVec::new(m); // unit vector: O(1) clear by its pattern
    // x_B round-off clamp bookkeeping (see `XbClamp`).
    // Hyper CHUZR (see `HyperChuzr`): lever + size gate.
    let mut hz = HyperChuzr::new(m, hchuzr_on(m));
    let mut xclamp = XbClamp::new(m, hz.on);
    // `in_basis` is maintained by the pivot's (leaving, entering) swap
    // instead of an O(total) rebuild per iteration; anything that moves the
    // basis wholesale (audit, backtrack) re-arms the rebuild.
    let mut in_basis_dirty = true;
    let mut in_basis_swap: Option<(usize, usize)> = None;
    // Dual steepest-edge weights w_i = ‖e_iᵀB⁻¹‖² per basis slot (Forrest &
    // Goldfarb 1992; HiGHS HEkkDual::updateFtranDSE/updateDual). Exactly 1
    // at the slack/artificial start (B = ±I); warm starts begin at 1 as
    // well (HiGHS "devex-like" reset). RUSOLVE_DUAL_DSE=0 = Dantzig lever.
    // DSE runs on COLD solves only. Warm (MIP node) re-solves restart the
    // weights at 1, which reshapes the tree (mixed tier: lights_out 1.77x
    // nodes) versus the Dantzig CHUZR they were tuned with.
    // Opt-in warm DSE (`reoptimize_dse`) carries
    // the weights in the snapshot as HiGHS does (HEkk::putIterate, HEkk.cpp:3812);
    // off by default: the MIP gate failed on tree reshaping, not on the weights.
    let dse_on = (cold || warm_dse.is_some())
        && std::env::var("RUSOLVE_DUAL_DSE")
            .map(|v| v != "0")
            .unwrap_or(true);
    let harris = std::env::var("RUSOLVE_DUAL_HARRIS")
        .map(|v| v != "0")
        .unwrap_or(true);
    // Probe lever: the best-pivot gate below which the Harris
    // relaxation engages (default 1e-5, the dano3_3-calibrated value).
    let harris_gate: f64 = std::env::var("RUSOLVE_DUAL_HARRIS_GATE")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1e-5);
    // Lever `RUSOLVE_DUAL_LARGE_ALPHA=1` (DEFAULT OFF): HiGHS
    // HEkkDualRow::chooseFinalLargeAlpha (MIT; `reference/HiGHS`, read as a
    // reference text). When the BFRT walk's final group holds no pivot
    // above min(0.1·max|ā| over the walked candidates, 1), back off to the
    // LAST earlier group that does — a shorter dual step (fewer flips)
    // instead of a tiny pivot. Motivated by irish-electricity, whose
    // blow-up starts at θ=0 pivots with |ā| ~ 1e-6.
    let large_alpha = std::env::var("RUSOLVE_DUAL_LARGE_ALPHA").is_ok_and(|v| v != "0");
    let mut gstarts: Vec<usize> = Vec::new();
    let mut n_large_alpha = 0usize;
    // Lever `RUSOLVE_DUAL_TA=1` (DEFAULT OFF): HiGHS's update-count
    // dependent CHUZC pivot tolerance (HEkkDualRow::choosePossible, MIT:
    // Ta = 1e-9 / 3e-8 / 1e-6 for < 10 / < 20 / ≥ 20 updates since the
    // last INVERT) in place of the flat `config.pivot_tol`.
    let ta_ladder = std::env::var("RUSOLVE_DUAL_TA").is_ok_and(|v| v != "0");
    let dinf_trace = std::env::var("RUSOLVE_DUAL_DINF").is_ok();
    // Diagnostic: per-visit state of one column (`RUSOLVE_DUAL_WATCH=j`).
    let watch_col: Option<usize> = std::env::var("RUSOLVE_DUAL_WATCH")
        .ok()
        .and_then(|v| v.trim().parse().ok());
    // `RUSOLVE_DUAL_SHIFT0` (zero-step cost shift, see the pivot site),
    // DEFAULT ON as mode 2 (shift + shiftBack); `=0` off, `=1` no shiftBack.
    // Measured: with FLIPPIV it takes dano3_3 to
    // Optimal on every perturbation seed and ex10 to Optimal at 26k it.
    let shift0_mode = std::env::var("RUSOLVE_DUAL_SHIFT0").unwrap_or_else(|_| "2".into());
    // Only where shifts persist while basic (`keep_shift`, cold solves): a
    // warm re-solve releases the entering column's shift on entry a few
    // lines after the pivot site sets it, so the shift there only
    // invalidates dj (a full recompute) — node-identical, ~2 % MIP time.
    let zero_step_shift = shift0_mode != "0" && keep_shift;
    let mut n_zero_shift = 0usize;
    // `RUSOLVE_DUAL_SHIFT0=2` adds HiGHS's shiftBack (HEkkDual::shiftBack,
    // HEkkDual.cpp:2215-2222, called from updateDual at :2120, MIT): the
    // zero-step shift a column picked up on entering is removed when it
    // LEAVES, so shifts do not accumulate over a long degenerate run (the
    // one-shot cleanup otherwise pays for every one of them).
    let zero_shift_back = zero_step_shift && shift0_mode == "2";
    // `RUSOLVE_DUAL_DEBTSHIFT` (default OFF): a flip-only / saturate /
    // accuracy-trouble step advances the duals to π(θ*) without a pivot, and
    // that was committed as the edge offset δ −= σθβ, β = B⁻ᵀe_p — a debt
    // only a later pivot on the SAME row pays; a pivot on any other row
    // drops it with the fresh BTRAN, leaving every flipped column at its
    // new bound with its old reduced-cost sign (seen on irish-electricity
    // at iteration 7,496). δ = −σθβ is exactly what a cost shift of −σθ on
    // the basic variable of row p does to π = c̃_Bᵀ B⁻¹ (HiGHS keeps every
    // such dual move in its shifted costs, HEkkDual::shiftCost,
    // HEkkDual.cpp:2072, MIT), so the lever books the debt there: π stays
    // exact across pivots on other rows, and the shift is backed out when
    // that variable leaves (the shiftBack slot `zshift`), i.e. the same-row
    // pivot pays it exactly as δ did. Never after the one-shot cleanup:
    // shifts added then would reach the Optimal certificate.
    let debt_shift = std::env::var("RUSOLVE_DUAL_DEBTSHIFT").is_ok_and(|v| v != "0");
    let mut n_debt_shift = 0usize;
    let mut zshift: Vec<f64> = if zero_shift_back || debt_shift {
        vec![0.0; st.cost.len()]
    } else {
        Vec::new()
    };
    let mut n_shift_back = 0usize;
    // `RUSOLVE_DUAL_FLIPPIV` (see the flip-prefix exits), DEFAULT ON
    // for cold solves only; `=0` restores the flip-only δ commit, `=1` also
    // pivots in warm (MIP node) re-solves — there the pivot makes the
    // final-group binary basic (fractional) where the flip left it at a
    // bound, and lights_out_mixed doubled its B&B nodes (seeds 1–3).
    let flip_pivot = match std::env::var("RUSOLVE_DUAL_FLIPPIV").as_deref() {
        Ok("0") => false,
        Ok("1") => true,
        _ => cold,
    };
    let mut n_flip_pivot = 0usize;
    // Levers `RUSOLVE_DUAL_CTIE` / `RUSOLVE_DUAL_RTIE`: break CHUZC
    // (equal |ā| in the pivot group) and CHUZR (equal DSE score) ties by
    // a fixed hashed permutation of the variable index instead of the
    // lowest index. HiGHS does the same: chooseFinalLargeAlpha breaks
    // equal-alpha ties by `workNumTotPermutation` (HEkkDualRow.cpp:512-520,
    // MIT; `reference/HiGHS`, read as a reference text), and its CHUZR
    // scan starts at a random row offset (HEkkDualRHS::chooseNormal).
    // Lowest-index ties are Bland-like and walk a ±1-data degenerate face
    // slowly (ex10).
    let ctie = std::env::var("RUSOLVE_DUAL_CTIE").is_ok_and(|v| v != "0");
    let rtie = std::env::var("RUSOLVE_DUAL_RTIE").is_ok_and(|v| v != "0");
    // Lever `RUSOLVE_DUAL_HGROUP`: HiGHS's Harris-relaxed BFRT
    // groups + backward large-alpha group choice (see the ratio site).
    let hgroup = std::env::var("RUSOLVE_DUAL_HGROUP").is_ok_and(|v| v != "0");
    // Lever `RUSOLVE_DUAL_CORRECT`: dual-infeasibility correction by
    // cost shifting at every computeDual (see the dj recompute site).
    let correct_dual = std::env::var("RUSOLVE_DUAL_CORRECT").is_ok_and(|v| v != "0");
    let mut n_correct = 0usize;
    let mut n_correct_max = 0.0f64;
    // Lagged z-progress: the ladder's
    // z test compares against the z of the last *gain*, so one dip (a
    // cost shift, a noisy refactor) re-arms the whole plateau window even
    // though z keeps climbing. ex10's optimal-face drain climbs z by
    // ~1e-9/visit (89 % of visits dz > 0, only 32 % above the ZFLOOR bar)
    // and walked DuPh1 → hand-over → Bland → distress. With
    // RUSOLVE_DUAL_ZNET=N a visit also counts as progress when z beats
    // the z of an evaluated visit ≥ N visits back by the same noise bar —
    // an exact cycle returns to its own z, so it never passes.
    let znet_lag: usize = std::env::var("RUSOLVE_DUAL_ZNET")
        .ok()
        .and_then(|v| v.trim().parse().ok())
        .unwrap_or(0);
    let mut znet_hist: std::collections::VecDeque<(usize, f64)> = Default::default();
    let mut n_znet = 0usize;
    let tie_key = |v: usize| -> u64 {
        (v as u64 ^ 0x9E37_79B9_7F4A_7C15)
            .wrapping_mul(0xBF58_476D_1CE4_E5B9)
            .rotate_left(31)
    };
    // Backtracking lever (see the update site), DEFAULT OFF:
    // `RUSOLVE_DUAL_BACKTRACK=1` arms it, and only above the dense-fallback
    // cap (below it the dense LU still rescues a sparse Singular). Measured
    // negative on irish-electricity: the restored basis replays the same
    // pivot sequence into the same ~3e-11 singular pivot at update limits
    // 5 → 2 → 1. Off ⇒ singular above
    // the cap goes straight to the primal fallback.
    let mut bt = Backtrack {
        on: m > crate::pfi::dense_fallback_max_n()
            && std::env::var("RUSOLVE_DUAL_BACKTRACK").is_ok_and(|v| v != "0"),
        max: 16,
        epoch: u64::MAX,
        ..Backtrack::default()
    };
    let mut dse_w = vec![1.0f64; m];
    if let Some((Some((snap, w)), _)) = &warm_dse
        && let Some(carried) = import_dse(lp, st, snap, w, &basis)
    {
        dse_w = carried;
    }
    let mut tau = HVec::new(m);

    let mut iters = 0usize;
    let mut stall = 0usize;
    // Degeneracy/plateau escalation (ladder):
    //  · V plateau (Σv drops <10% for `bland_after` visits) or stall ≥
    //    bland_after, with V > 0 ⇒ DuPh1 selection overlay: wider
    //    entering tie window scored by expected V-absorption
    //    min(ũ·|ā|, v_p), sticky until primal feasibility (then the
    //    optimality audit decides). The tuff autopsy forced the plateau trigger: the
    //    max-violation loop ran an EXACT 2-cycle (V frozen at 8.5, pivots
    //    landing l=8.5/2.5) that a degenerate-pivot stall counter can never
    //    see.
    //  · DuPh1 with no V drop for 4·bland_after visits ⇒ sticky-random
    //    leaving among the top-3 violating rows (deterministic XorShift64,
    //    primal-engine seeding).
    //  · esc ≥ 2·bland_after ⇒ cost-perturbation rung: engage
    //    the hash-signed sub-tolerance cost tilt ONCE (one step after
    //    DuPh1); it lives for the rest of the solve and the true-cost
    //    optimality audit is its cleanup by construction. Breaks exact θ=0
    //    cycles, which no selection rule can.
    //  · stall ≥ 8·bland_after ⇒ dual Bland (lowest-variable-index leaving
    //    row, lowest-index min-ratio entering; Chvátal §4.7/§7.3). Sticky.
    //  · stall ≥ 16·bland_after ⇒ NumericalDistress ⇒ primal fallback.
    let mut bland_triggered = false; // escalation engaged (DuPh1 or Bland)
    let mut bland_bland = false; // dual-Bland sticky flag (deterministic rung)
    let mut rng =
        crate::rng::XorShift64::new(0x9E37_79B9_7F4A_7C15 ^ ((m as u64) << 32 | total as u64));
    let mut duph1 = false;
    let mut duph1_iters = 0usize;
    let mut no_v_drop = 0usize;
    let mut v_anchor = f64::INFINITY;
    let mut z_anchor = f64::NEG_INFINITY;
    let mut z_reanchors = 0usize;
    let mut srnd = false; // sticky-random top-3 rung (sticky-on per episode)
    let mut sticky_v: Option<usize> = None;
    // Cost-perturbation rung: hash-signed
    // sub-tolerance cost perturbation, engaged ONCE on a plateau/stall the
    // tie-scoring overlay could not clear (esc ≥ 2·bland_after) and kept
    // for the rest of the solve — the optimality audit prices TRUE costs
    // (`cost + s`, never `pert`), so it is the audit cleanup pass by
    // construction and Optimal is only ever certified at true-cost KKT
    // (the honesty rule). The market_split warm-reopt trace showed
    // exact θ=0 4-visit cycles that DuPh1 and sticky-random pass through
    // UNCHANGED — a frozen dual point makes every selection rule re-index
    // the same transition graph, so only tilting π off the degenerate face
    // can break such a cycle. `RUSOLVE_DUAL_NOPERT=1` is the lever back.
    let perturb_on = std::env::var("RUSOLVE_DUAL_NOPERT").is_err();
    // z-gain noise floor: the
    // plateau ladder's "real dual progress" test is relative
    // (1e-9·(1+|z|)), which at |z|≈6e9 (pds-*) demands +6 per window while
    // honest perturbation-scale progress is 1e-2..1e-1 — a false plateau
    // that walks a progressing solve to Bland/distress. The floor caps the
    // threshold at K·ε·Σ|terms|, the rounding envelope of the z sum itself,
    // so only gains rounding cannot explain count. K=0 restores the old test.
    let zgain_floor_k = zgain_floor_k();
    // Fresh-π ratio-pass rescan skip: `RUSOLVE_DUAL_RESCAN_SKIP=1`.
    // Flag-off: not bit-identical (a skipped column can miss an ulp-level
    // re-bump) and no measured win on s250r10.
    let ratio_rescan_skip = std::env::var("RUSOLVE_DUAL_RESCAN_SKIP")
        .map(|v| v.trim() == "1")
        .unwrap_or(false);
    // Lazy z-pass threshold (see the ladder): 0 = evaluate every visit.
    let zlazy_at = zlazy_at(config.bland_after, m + total >= ZLAZY_MIN_DIM);
    let mut n_zpass = 0usize;
    let mut pert_used = st.pert.iter().any(|&p| p != 0.0); // rung engaged at most once per solve
    // One fresh-factor verification allowed per Infeasible / Unbounded
    // candidate (the ray_verified pattern of the primal engine).
    let mut infeas_verified = false;
    let mut unbounded_verified = false;
    let mut cleanup_done = false;

    // Pivotal-row PRICE (HEkkDual::chooseColumn -> HEkk::tableauRowPrice,
    // simplex/HEkkDual.cpp:1546, HEkk.cpp:2845) from a row-wise copy, and
    // incrementally maintained reduced costs (HEkkDual::updateDual
    // HEkkDual.cpp:2065 / HEkkDualRow::updateDual HEkkDualRow.cpp:548);
    // recomputed from the factors at every refactorization and after
    // anything outside the pivot loop touches the duals (audit,
    // perturbation, verification refresh). Levers: RUSOLVE_DUAL_ROWPRICE=0
    // (column-wise PRICE), RUSOLVE_DUAL_INCD=0 (reduced costs from pi every
    // pass). The incremental d needs no per-pass BTRAN of the basic costs;
    // the warm cutoff bound does, so warm solves keep the fresh-pi path.
    let big = st.csc_row.len() >= DUAL_FAST_MIN_NNZ;
    let rowprice = big
        && std::env::var("RUSOLVE_DUAL_ROWPRICE")
            .map(|v| v != "0")
            .unwrap_or(true);
    // Incremental d + the incremental dual objective now also run on
    // WARM solves with a cutoff — the cutoff test no longer needs fresh π
    // every iteration (see the z bookkeeping below). Threshold decided by
    // A/B on our benchmark set: 10k engages the ≥15k-nnz
    // node LPs (cap_facility class) where the removed π-BTRAN + O(nnz)
    // scan is a pure per-node win (0.77×, bit-identical trees), while the
    // tiny-LP families (gap ≈1k nnz standard form) keep their paths —
    // engaging there sped per-node work but CHANGED TREE SHAPE (gap nodes
    // ×1.5-1.9, one lost solve): the eliminated tax is tiny in absolute
    // terms while the path-sensitivity risk is size-independent.
    let warm_fast_min_nnz = std::env::var("RUSOLVE_WARM_INCD_MIN_NNZ")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(10_000);
    let fast_min_nnz = if cold {
        DUAL_FAST_MIN_NNZ
    } else {
        warm_fast_min_nnz
    };
    let mut incd = st.csc_row.len() >= fast_min_nnz
        && std::env::var("RUSOLVE_DUAL_INCD")
            .map(|v| v != "0")
            .unwrap_or(true);
    // Safety net: skipping the per-pass BTRAN changes the refactor schedule
    // and the rounding, and degenerate LPs are path-sensitive (dano3_3 went
    // from 38k iterations to the 100k limit). The first ladder rung hands
    // the SAME solve over to the fresh-pi path from the current basis
    // (nothing else changes: x_b, at_ub, the shifts and the committed
    // offset δ are kept by both paths). Warm
    // solves get the same rung (their escape used to be only the MIP-level
    // 2×cap cold fallback); cold solves with an iteration limit keep the
    // historic gate so probe caps stay bit-stable.
    // Dual accuracy regime (`RUSOLVE_DUAL_ACC=1` / `SimplexConfig::
    // dual_acc`). The question: are our 6-8 raw solves/iteration (HiGHS: 3-4) the coupled
    // symptom of an inaccurate factorization/update? Census answer: NO.
    // |α_col − α_row|/min stays ≤ 3e-11 on the UNREFINED solves of
    // qap10/dano3_3/istanbul, four orders under HiGHS's trouble tolerance.
    // The refinement corrections are load-bearing for a different reason:
    // on 0/±1 structure the refined solves return EXACT ties, and the
    // CHUZR/ratio-group selections break those ties by index; unrefined
    // ulp noise turns the same choices into noise-driven ones and the
    // ladder grinds (graph20 ×44 iterations, istanbul ×2.5). So the package
    // is:
    //   1  unrefined loop solves (HiGHS runs no refinement at all),
    //   8  the per-pivot accuracy monitor (HEkkDual::updateVerify,
    //      HEkkDual.cpp:2050 → HEkk::reinvertOnNumericalTrouble,
    //      HEkk.cpp:3015; kNumericalTroubleTolerance = 1e-7,
    //      util/HFactorConst.h:57): on trouble the pivot is abandoned,
    //      reinverted, x_B recomputed, and refinement re-armed for the rest
    //      of the solve (see the pivot site),
    //  16  tolerance ties (relative 1e-9) in the CHUZR score and the ratio
    //      group's |ā| choice, broken by smallest index: Bland's rule
    //      (Bland 1977, Math. Oper. Res. 2:103-107) applied to the
    //      tolerance-equality class, so unrefined noise cannot split a tie
    //      the refined solves would have reported exact. tie = 0 is the
    //      historic exact comparison, bit for bit.
    // Default mask 25 = 1|8|16. Measured negative components, kept for
    // retry (`RUSOLVE_DUAL_ACC_PARTS`): 2 computePrimal at every reinversion
    // (HEkkDual::rebuild, HEkkDual.cpp:1086-1091; breaks qap10/dano3_3 with
    // unrefined solves), 4 no fresh-π hand-over (doubles qap10 iterations,
    // no ms/it gain). `RUSOLVE_DUAL_ACC=2` forces the hand-over.
    // HiGHS (MIT, `reference/HiGHS`) is read as a reference text only.
    // Default ON once the robust-ratio levers made dano3_3 robust under it
    // (A/B on our benchmark set: 47 relaxed MIPLIB LPs, dual shifted
    // geometric mean 0.974, 38.2 → 29.2 s, no wrong answer; MIP set
    // medium 0.90 / mixed 1.00; square41 ×0.70, ex10 ×0.59, supportcase10
    // ×0.66, pds-40 ×1.13). `RUSOLVE_DUAL_ACC=0` restores refined solves.
    let acc_env = std::env::var("RUSOLVE_DUAL_ACC").ok();
    let acc = config.dual_acc || acc_env.as_deref().is_none_or(|v| v != "0");
    let acc_parts: u32 = if acc {
        std::env::var("RUSOLVE_DUAL_ACC_PARTS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(25)
    } else {
        0
    };
    let acc_handover = acc_env.as_deref() == Some("2") || acc_parts & 4 == 0;
    let acc_primal = acc_parts & 2 != 0;
    let tie = if acc_parts & 16 != 0 { 1e-9 } else { 0.0 };
    let mut refine_on = acc_parts & 1 == 0;
    let mut n_trouble = 0usize;
    let hand_over = incd
        && warm_dse.is_none()
        && (iter_limit.is_none() || !cold)
        && (!acc || acc_handover)
        && std::env::var("RUSOLVE_DUAL_HANDOVER")
            .map(|v| v != "0")
            .unwrap_or(true);
    // Incremental-π transport (`RUSOLVE_DUAL_INCP=1`, DEFAULT OFF).
    // A census found the fresh-π path
    // re-anchoring π by a full BTRAN of the basic costs on EVERY iteration
    // — 95% of qap10 iterations, 73 μs/it = 18% of its wall (and qap10 is
    // 74% of the corpus gap). The maintained dual point (π̂ = piw, δ folded
    // in) is already carried across the loop: flip/saturate exits advance
    // it by the edge-travel step (piw -= σθβ), and the pivot's re-anchor is
    // in exact arithmetic NOTHING BUT the rank-one dual step π̂ ← π̂ + t·β,
    // t = d̃_q/ā_q — the same unclamped step `update_dj` applies to d̃
    // (d̃_j -= t·ā_j for every priced column ⇔ Δπ̂ = t·β, since d̃_j is
    // c̃_j − π̂ᵀa_j and ā_j = βᵀa_j; β-orthogonality of the other basics
    // makes the leaving column's d̃ = −t fall out too, matching the
    // `dj[leaving] = -t` bookkeeping below). So with the lever on, pivots
    // advance π̂ by that rank-one step and the BTRAN fires only at a
    // re-anchor — events + cadence per [`PiAnchor`], cadence
    // `RUSOLVE_DUAL_INCP_CAD` (default 20 pivots = qap10's refactor
    // period). Lineage: this is HiGHS's own policy — duals maintained
    // incrementally across pivots, computeDual (a fresh BTRAN) only at
    // reinversions/invalidations (HEkkDual::updateDual HEkkDual.cpp:2065,
    // HEkk::computeDual; MIT, reimplemented).
    // Staleness is bounded by the anchor interval; the escalation ladder
    // (DuPh1 → cost perturbation → Bland → distress⇒primal, the
    // repair-cascade guard) is π-agnostic and still catches mispricing
    // cascades. (A ladder-engages-anchor lever was tried and removed:
    // it neither rescued dano3_3's iteration blow-up nor helped qap10.)
    // Flag OFF ⇒
    // bit-identical: `need_pi` keeps its exact original expression and the
    // pivot arm touches nothing.
    let incp = std::env::var("RUSOLVE_DUAL_INCP").is_ok_and(|v| v != "0");
    let incp_cad = std::env::var("RUSOLVE_DUAL_INCP_CAD")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&c| c > 0)
        .unwrap_or(20);
    // Measurement lever: at every re-anchor, compare the maintained π̂
    // (before it is overwritten) against the anchored point the re-anchor
    // produces — the drift the anchor just erased, ‖π̂ − (BTRAN π + δ)‖∞ /
    // (1 + ‖·‖∞). Rounding-level drift ⇒ trajectory moves are degenerate-LP
    // path-sensitivity, not mispricing; growing drift ⇒ structural bug in
    // the rank-one path. Default OFF, zero cost when unset.
    let picheck = incp && std::env::var("RUSOLVE_DUAL_PICHECK").is_ok_and(|v| v != "0");
    let mut pi_ref = if picheck { vec![0.0; m] } else { Vec::new() };
    let (mut pic_max, mut pic_sum, mut pic_n, mut pic_evt_max) = (0.0f64, 0.0f64, 0usize, 0.0f64);
    let mut anchor = PiAnchor::new(incp_cad);
    let rowmat = if rowprice {
        let mark = std::time::Instant::now();
        let rm = RowMat::from_dual(st);
        let t = mark.elapsed();
        if warm_stats {
            census_acc(|c| c.t_rowmat += t.as_nanos() as f64);
        }
        Some(rm)
    } else {
        None
    };
    let mut pbuf = PriceBuf::new(total);
    // Incremental reduced costs across consecutive optimality audit primal
    // steps (the post-cleanup primal tail). Needs the row-wise PRICE copy.
    let mut acache = AuditCache::new(total, m, rowmat.is_some() && audit_inc_on());
    let mut passes = 0u64; // monotone loop-pass counter (iters can step back)
    let mut dj = vec![0.0f64; total];
    let mut dj_valid = false;
    // Incremental dual objective (the spec's "cutoff test without π"):
    // z_val tracks the caller-space TRUE-cost bound
    //   konst ± (πᵀb̃ + Σ_{nonbasic} d_j·x̃_j)
    // of the MAINTAINED dual point (π + δ, dj) — an identity, so it is
    // exactly `caller_dual_bound`'s box-min while the dual-simplex
    // feasibility invariant holds (at-lower d ≥ 0, at-upper d ≤ 0), and
    // both differ only by tolerance-level dual debt otherwise. Events:
    //   · pivot:      Δz = t·x_b[p]⁰ + Σ_flips d'_j·Δx̃_j + d'_r·exit_bound
    //                (+ (s_q + pert_q)·x̃_q⁰ for a shifted/perturbed enterer)
    //   · δ-commit /
    //     saturate:   Δz = −σθ·x_b[p]⁰ + Σ_flips d'_j·Δx̃_j (+ q's own flip)
    // with t = dj[q]/ā_q the unclamped dual step (the same t update_dj
    // uses), x_b[p]⁰ the PRE-event basic value of the leaving row, and
    // d'_j the post-step true-cost reduced cost dj_j + coef·α_j − s_j −
    // pert_j. Derivation: the π-term and the Σ Δd_j x̃_j term of the
    // identity telescope against βᵀb̃ − Σ ā_j x̃_j = x_b[p]. Lineage: HiGHS
    // keeps a dual objective incrementally across pivots
    // (HEkkDual::updateDual, HEkkDual.cpp:2065; HEkkDualRow::updateDual
    // HEkkDualRow.cpp:548; computeDualObjectiveValue HEkk.cpp:1948 — MIT,
    // reimplemented: HiGHS's nonbasic-sum form drops the πᵀb̃ term and the
    // leaving-variable exit-bound term because theirs only feeds REPORTING;
    // ours feeds prune decisions, so both terms are kept and every crossing
    // is re-verified against the exact `caller_dual_bound` before pruning).
    // Seeded wherever dj is recomputed from the factors (refactor cadence,
    // audit, perturbation) — O(m + at-upper count), no extra csc pass.
    // RUSOLVE_DUAL_ZCHECK=1 recomputes-and-compares every iteration (tests).
    let z_on = cutoff.is_some();
    let konst: f64 = if z_on {
        lp.objective.iter().zip(&st.shift).map(|(c, s)| c * s).sum()
    } else {
        0.0
    };
    let mut z_val = 0.0f64;
    let mut z_valid = false;
    let mut z_last_t: Option<f64> = None; // last pivot's dual step (zcheck trace)
    let ztrace = std::env::var("RUSOLVE_DUAL_ZTRACE").is_ok();
    // ztrace prints each event's Δz alongside the state identity, so a bad
    // term is localizable to its event type (pivot / commit / saturate).
    macro_rules! ztrace {
        ($kind:expr, $dz:expr) => {
            if ztrace {
                eprintln!("zevent it={iters} {} dz={:.9e}", $kind, $dz);
            }
        };
    }
    let mut pi_v = if z_on { vec![0.0f64; m] } else { Vec::new() };
    let zcheck = std::env::var("RUSOLVE_DUAL_ZCHECK").is_ok();
    // Which leaving row the committed edge offset δ came from (NO_ROW =
    // δ is zero, MIXED_ROW = several rows or unknown). A pivot on row p
    // keeps d incrementally exact only if δ ∝ β_p (every other basic
    // column is β-orthogonal); otherwise the fresh BTRAN of the new basis
    // DROPS the δ debt on the old row's basic column and d is recomputed.
    const NO_ROW: usize = usize::MAX;
    const MIXED_ROW: usize = usize::MAX - 1;
    let mut delta_row = NO_ROW;
    let mut cands: Vec<(usize, f64, f64, bool)> = Vec::new();
    // Hyper-sparse CHUZR (Hall & McKinnon 2005; HiGHS HEkkDualRHS, MIT):
    // choose the leaving row from the infeasible-basics list when it is
    // valid and small, full 0..m scan (doubling as rebuild) otherwise.
    // Lever: RUSOLVE_DUAL_INF_LIST=<fraction of m> (0 disables).
    // Default 0.05 calibrated on 47 relaxed MIPLIB LPs (interleaved): the list
    // wins on the time-dominant instances (dano3_3 −11%, graph20 −12%,
    // qap10 −4%) and loses above ~0.10 on the 15-25%-violator mid-dense
    // ones (krka, rmatr100-p10, istanbul: random gathers don't beat a
    // streaming scan when the row arrays stay cache-resident); 0.25 put
    // the worst instance at +9% vs +2% at 0.05.
    let list_thr_m = (std::env::var("RUSOLVE_DUAL_INF_LIST")
        .map(|v| v.parse::<f64>().unwrap_or(0.0))
        .unwrap_or(0.05)
        * m as f64) as usize;
    let mut infeas = InfeasList::new(m);
    let mut infeas_ok = false;
    // RUSOLVE_DUAL_STATS=1: per-solve cost counters (per-iteration
    // profiling; off by default, no cost when unset beyond the adds).
    let stats = std::env::var("RUSOLVE_DUAL_STATS").is_ok();
    // Sweep-cost census: periodic mid-run stats dump, only
    // meaningful together with RUSOLVE_DUAL_STATS (see the loop head).
    let stats_dump_every: usize = if stats {
        std::env::var("RUSOLVE_STATS_EVERY")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    } else {
        0
    };
    let (mut cand_sum, mut cand_max, mut n_ratio_sort, mut n_price, mut price_pat_sum) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut n_scan, mut n_scan_fb, mut scan_groups_max) = (0usize, 0usize, 0usize);
    // Degenerate-pivot count (the 3488-classification below), for the
    // LP-gap decomposition. Touched only under `stats`.
    let mut n_degen = 0usize;
    let (mut t_chuzr, mut t_beta_bt, mut t_price, mut t_ratio, mut t_qcol, mut t_dse, mut t_upd) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    // Sub-splits of the chuzr region: the 0..m violation scan versus the
    // O(m+total) z-plateau re-anchor pass (which the timer used to swallow).
    let (mut t_scan, mut t_zpass) = (0.0f64, 0.0f64);
    let (mut n_listscan, mut n_fullscan, mut infeas_len_sum) = (0usize, 0usize, 0usize);
    // Profiling: per-iteration event rates for the transport phases.
    let (mut n_pi_bt, mut n_beta_nz, mut n_beta_solves) = (0usize, 0usize, 0usize);
    // Census: dj invalidations by cause [hand-over, perturbation, audit,
    // verify, shift-release, delta-row, emergency-refactor x2].
    let mut inv_cnt = [0usize; 8];
    // Re-anchor cause split under RUSOLVE_DUAL_INCP (event vs cadence).
    let (mut n_anchor_evt, mut n_anchor_cad) = (0usize, 0usize);
    // HiGHS-style accuracy monitor census (HEkk::reinvertOnNumericalTrouble,
    // HEkk.cpp:3015): |α_col − α_row| / min(|α_col|, |α_row|) per pivot, and
    // the DSE weight error |w̃_p − ‖β‖²| / ‖β‖² — stats-only counters.
    let (mut n_acc, mut n_acc7, mut n_acc9, mut acc_max) = (0usize, 0usize, 0usize, 0.0f64);
    let (mut n_dsew_bad, mut dsew_max) = (0usize, 0.0f64);
    let (mut n_xdrift, mut n_xdrift6, mut xdrift_max) = (0usize, 0usize, 0.0f64);
    let (mut n_revealed, mut revealed_max) = (0usize, 0.0f64);
    let mut t_flip = 0.0f64;
    let mut t_sel = 0.0f64;
    let (mut flip_sum, mut flip_max, mut flip_batches) = (0usize, 0usize, 0usize);
    macro_rules! acc {
        ($acc:expr, $mark:expr) => {
            if stats {
                $acc += $mark.elapsed().as_nanos() as f64;
            }
        };
    }

    // The per-region time line (Optimal exit + `RUSOLVE_STATS_EVERY` dumps:
    // cap-hitting giant LPs never reach the Optimal exit).
    macro_rules! dual_time_line {
        () => {{
            let wall = t0.elapsed().as_secs_f64() * 1e3;
            let rf = pfi.refactor_nanos() / 1e6;
            eprintln!(
                "rusolve dual time: wall={wall:.0}ms chuzr={:.0} (scan={:.0} zpass={:.0}) \
                 sel={:.0} bt_beta={:.0} \
                 price={:.0} ratio={:.0} flip={:.0} qcol={:.0} dse={:.0} \
                 update={:.0} (of which refactor={rf:.0}ms) \
                 flip_batches={flip_batches} flip_mean={:.1} flip_max={flip_max} \
                 pi_bt={n_pi_bt} reanchor_pct={:.1} incp={incp} \
                 anchor_evt={n_anchor_evt} anchor_cad={n_anchor_cad} \
                 pic_drift={:.2e}/{:.2e} evt={:.2e} beta_nz_mean={:.0} inv={inv_cnt:?} \
                 audit_inc={}/{} zpass_n={n_zpass}",
                t_chuzr / 1e6,
                t_scan / 1e6,
                t_zpass / 1e6,
                t_sel / 1e6,
                t_beta_bt / 1e6,
                t_price / 1e6,
                t_ratio / 1e6,
                t_flip / 1e6,
                t_qcol / 1e6,
                t_dse / 1e6,
                t_upd / 1e6,
                flip_sum as f64 / flip_batches.max(1) as f64,
                100.0 * n_pi_bt as f64 / iters.max(1) as f64,
                pic_max,
                pic_sum / pic_n.max(1) as f64,
                pic_evt_max,
                n_beta_nz as f64 / n_beta_solves.max(1) as f64,
                acache.hits,
                acache.fresh,
            );
        }};
    }

    let max_it = iter_limit.unwrap_or(config.max_iterations);
    'outer: loop {
        iters += 1;
        passes += 1;
        // The clock is read on the first and then every 16th iteration
        // (warm MIP node LPs now carry the search deadline too, and
        // run this loop hot; the first read keeps an expired deadline's
        // exit immediate, as a test pins).
        if iters > max_it
            || (iters % 16 == 1 && deadline.is_some_and(|d| std::time::Instant::now() >= d))
        {
            return Ok(Some(DualExit::Sol(limit_dual(
                st,
                lp,
                &basis,
                &x_b,
                &at_ub,
                duph1_iters,
                iters,
                bland_triggered,
            ))));
        }
        if in_basis_dirty {
            for b in in_basis.iter_mut() {
                *b = false;
            }
            for &v in basis.iter() {
                in_basis[v] = true;
            }
            in_basis_dirty = false;
            in_basis_swap = None;
        } else if let Some((leaving, entering)) = in_basis_swap.take() {
            in_basis[leaving] = false;
            in_basis[entering] = true;
        }
        debug_assert!(
            in_basis.iter().filter(|&&b| b).count() == m && basis.iter().all(|&v| in_basis[v]),
            "incremental in_basis out of sync with basis"
        );
        // Backtracking snapshot: the state right after every successful
        // (re)factorization (HiGHS HEkk::putBacktrackingBasis — basis, cost
        // shifts and edge weights). Taken at the loop head, where basis,
        // sides, shifts and the fresh factors agree.
        if bt.on && pfi.eta_count() == 0 && pfi.factor_epoch() != bt.epoch {
            bt.epoch = pfi.factor_epoch();
            bt.valid = true;
            bt.basis.clone_from(&basis);
            bt.at_ub.clone_from(&at_ub);
            bt.s.clone_from(&st.s);
            bt.dse_w.clone_from(&dse_w);
            bt.delta.clear();
            bt.delta.extend_from_slice(delta.vals());
        }
        // Escalation ladder on esc = max(degenerate-pivot stall, V-plateau
        // visits): the tuff cycle mixes θ=0 and θ=3e-5 pivots so `stall`
        // oscillates 0↔1 forever — only the V-plateau counter is monotone
        // under it. (no_v_drop is the previous visit's; off-by-one is
        // harmless at these thresholds.)
        if dprog && iters.is_multiple_of(500) {
            eprintln!(
                "rusolve dual: it={iters} t={:.1}s z_anchor={z_anchor:.6e} stall={stall} no_v_drop={no_v_drop} v_anchor={v_anchor:.4e} debt_drops={} flip_pivot={n_flip_pivot} debt_shift={n_debt_shift}",
                t0.elapsed().as_secs_f64(),
                inv_cnt[5]
            );
        }
        // Sweep-cost census: periodic stats dump for cap-hitting
        // giant-LP runs (cont1-class) that never reach the Optimal-exit
        // stats line. Requires RUSOLVE_DUAL_STATS; unset ⇒ zero cost.
        if stats_dump_every > 0 && iters.is_multiple_of(stats_dump_every) {
            eprintln!(
                "rusolve dual @it={iters} t={:.1}s",
                t0.elapsed().as_secs_f64()
            );
            eprintln!("{}", crate::pfi::pfi_stats_line());
            eprintln!(
                "rusolve pfi live: etas={} live_bytes={} lu_work={}",
                pfi.eta_count(),
                pfi.eta_live_bytes(),
                pfi.lu_solve_work(),
            );
            eprintln!("{}", crate::lu::lu_fact_line());
            dual_time_line!();
            // The per-phase timers too, cumulative to this point — a
            // giant LP that never reaches the Optimal-exit `dual time` line
            // still needs its ms/it attribution (ex10 / irish-electricity).
            eprintln!(
                "rusolve dual time@: duph1={duph1_iters} refactors={} chuzr={:.0} sel={:.0} \
                 bt_beta={:.0} price={:.0} ratio={:.0} flip={:.0} qcol={:.0} dse={:.0} \
                 update={:.0} (refactor={:.0}ms) pi_bt={n_pi_bt} beta_nz_mean={:.0} \
                 chuzr_list_scans={n_listscan} chuzr_full_scans={n_fullscan} \
                 large_alpha={n_large_alpha} hchuzr_hits={} hchuzr_builds={} \
                 hchuzr_set_mean={:.0} hchuzr_k={} hchuzr_off={} zero_shift={n_zero_shift} \
                 correct={n_correct} correct_max={n_correct_max:.2e} znet={n_znet} shift_back={n_shift_back} flip_pivot={n_flip_pivot} debt_shift={n_debt_shift}",
                pfi.factor_epoch(),
                t_chuzr / 1e6,
                t_sel / 1e6,
                t_beta_bt / 1e6,
                t_price / 1e6,
                t_ratio / 1e6,
                t_flip / 1e6,
                t_qcol / 1e6,
                t_dse / 1e6,
                t_upd / 1e6,
                pfi.refactor_nanos() / 1e6,
                n_beta_nz as f64 / n_beta_solves.max(1) as f64,
                hz.n_hit,
                hz.n_build,
                hz.set_sum as f64 / hz.n_hit.max(1) as f64,
                hz.k,
                hz.n_off,
            );
        }
        let esc = stall.max(no_v_drop);
        if esc >= 16 * config.bland_after {
            // Always logged: a silent dual→primal
            // handoff made rescue rows look like dual parity in sweeps.
            let _ = debug;
            eprintln!(
                "rusolve dual: distress after {esc} stalled/plateaued visits ⇒ primal fallback"
            );
            return Ok(None);
        }
        let use_bland = bland_bland;

        // ---- (a) Leaving row: max violation; DuPh1 sticky-random top-3
        // override; Bland: lowest violating variable index ----
        let mut ppos = usize::MAX;
        let mut vbest = 0.0f64;
        let mut sbest = 0.0f64; // CHUZR score: v²/w (DSE) or v (Dantzig)
        let mut vtotal = 0.0f64;
        // Top-3 violating rows by violation (sticky-random's candidate set),
        // as (v, pos) descending.
        let mut top3: [(f64, usize); 3] = [(0.0, usize::MAX); 3];
        let tmark = if stats {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let smark = tmark;
        // Choose from the infeasible-basics list when it is valid and
        // small enough to beat the full scan (Hall & McKinnon 2005; HiGHS
        // HEkkDualRHS::chooseNormal SPARSE mode). The full-scan branch
        // doubles as the list rebuild; rows at or below FEAS_TOL drop out
        // there (their ≤ tol·m mass is also omitted from vtotal on list
        // passes — a documented, tolerance-level change to the plateau
        // ladder's Σv).
        let use_list = list_thr_m > 0 && infeas_ok && infeas.list.len() <= list_thr_m;
        // Hyper CHUZR: exact-argmax shortcut for the full-scan branch
        // (tie = 0 comparator only; Bland / sticky-random keep the scan).
        let hz_try =
            hz.on && !use_list && !use_bland && !srnd && tie == 0.0 && !ub_b_dirty && hz.gate();
        let hz_res = if hz_try {
            hz.choose(&x_b, &ub_b, &dse_w, dse_on, &basis, &mut xclamp, list_thr_m)
        } else {
            None
        };
        let hz_hit = hz_res.is_some();
        macro_rules! chuzr_row {
            ($i:expr) => {{
                let i = $i;
                let v = violation(ub_b[i], x_b[i]);
                vtotal += v;
                if v > FEAS_TOL {
                    if v > top3[2].0 {
                        top3[2] = (v, i);
                        top3.sort_by(|a, b| b.0.total_cmp(&a.0));
                    }
                    let better = if ppos == usize::MAX {
                        true
                    } else if use_bland {
                        basis[i] < basis[ppos]
                    } else {
                        let sc = if dse_on { v * v / dse_w[i] } else { v };
                        sc > sbest * (1.0 + tie)
                            || (sc >= sbest * (1.0 - tie)
                                && if rtie {
                                    tie_key(basis[i]) < tie_key(basis[ppos])
                                } else {
                                    basis[i] < basis[ppos]
                                })
                    };
                    if better {
                        vbest = v;
                        sbest = if dse_on { v * v / dse_w[i] } else { v };
                        ppos = i;
                    }
                }
                v
            }};
        }
        if let Some((hp, hv, _, ht)) = hz_res {
            ppos = hp;
            vbest = hv;
            vtotal = ht;
            hz.prev = hp;
            // The full scan this replaces would have emptied the list
            // (begin_rebuild) and left it disabled (> thr violators).
            infeas_ok = false;
        } else if use_list {
            if stats {
                n_listscan += 1;
                infeas_len_sum += infeas.list.len();
            }
            let n = infeas.list.len();
            for k in 0..n {
                let i = infeas.list[k];
                chuzr_row!(i);
            }
        } else {
            if stats {
                n_fullscan += 1;
            }
            // Rebuild costs a push per violator — a dense set would pay it
            // every single scan for a list it never uses. So dense scans
            // only COUNT (one add on the branch the scan already takes);
            // the iteration the set turns sparse pays one extra O(m) pass
            // to build. Amortized: one build per dense→sparse transition.
            infeas.begin_rebuild(); // stale entries ≤ thr deep: cheap
            if ub_b_dirty {
                for (i, &v) in basis.iter().enumerate() {
                    ub_b[i] = st.ub[v];
                }
                ub_b_dirty = false;
            }
            if list_thr_m > 0 {
                let mut cnt = 0usize;
                for i in 0..m {
                    if chuzr_row!(i) > FEAS_TOL {
                        cnt += 1;
                    }
                }
                if cnt <= list_thr_m {
                    for i in 0..m {
                        if violation(ub_b[i], x_b[i]) > FEAS_TOL {
                            infeas.push(i);
                        }
                    }
                    infeas_ok = true;
                } else {
                    // The list was just emptied by begin_rebuild. Leaving it
                    // enabled would let the next CHUZR scan an EMPTY list and
                    // read "no violators" (phantom audit) while `cnt` rows
                    // violate — the over-approximation invariant only
                    // survives rebuilds that push. Dense until it turns
                    // sparse: disable until the next successful rebuild.
                    infeas_ok = false;
                }
            } else {
                for i in 0..m {
                    chuzr_row!(i);
                }
            }
            if hz_try {
                hz.build(&x_b, &ub_b, &dse_w, dse_on);
                xclamp.reset_chg();
                hz.prev = ppos;
            }
        }
        if !hz_try {
            // A visit the hyper path sat out may clamp its own x_b[p]
            // unlogged (flip-only exit): drop the structure.
            hz.valid = false;
        }
        // List stays exact-as-over-approximation only if every x_b change is
        // seen; maintain it from the update loops while it is in play.
        let maintain = list_thr_m > 0 && infeas_ok && infeas.list.len() <= list_thr_m;
        if let Some(mk) = smark {
            acc!(t_scan, mk);
        }
        // V plateau bookkeeping (dual-objective-plateau proxy: Σv is the
        // dual phase-1 objective up to row scaling). A new low (≥10% under
        // the anchor) RESETS the whole ladder — the degenerate pocket
        // demonstrably broke, so DuPh1/srnd/Bland release (the primal
        // engine's lesson: escalation that never releases is itself a defect).
        if vbest > FEAS_TOL {
            // Dual-objective progress (the real dual progress measure, as
            // Σv is not monotone under steepest edge / flips): z = cᵀx at
            // the current basic solution with the pivot-loop costs. It is
            // nondecreasing along dual pivots, so strict gains rule out
            // cycling and reset the ladder like a V drop does.
            // Lazy z: the pass is
            // O(m + total) per visit — 20-28% of pds-100/stormG2 wall — but
            // its only consumer is the ladder, whose first rung sits at
            // esc = bland_after. Below `zlazy_at` the visit cannot reach a
            // rung, so the pass is skipped and the visit counts as no-gain;
            // the cumulative anchor test (z vs the z of the last gain) then
            // catches up on the first evaluated visit: a dual-progressing
            // solve is reset there, a plateau escalates at the same rungs,
            // at most `zlazy_at` visits later than the per-visit test. Not
            // trajectory-neutral (a skipped gain does not re-anchor V), hence
            // the ZLAZY_MIN_DIM size gate.
            let need_z = zlazy_at == 0
                || dtrace
                || z_anchor == f64::NEG_INFINITY
                || stall.max(no_v_drop) + 1 >= zlazy_at;
            let (z, z_gain, z_thr) = if need_z {
                let zmark = if stats {
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                let (z, zabs) = ladder_z(st, &basis, &x_b, &at_ub, &in_basis, total);
                if let Some(mk) = zmark {
                    acc!(t_zpass, mk);
                }
                n_zpass += 1;
                let z_thr = zgain_threshold(z_anchor, zabs, zgain_floor_k);
                (
                    z,
                    z_anchor == f64::NEG_INFINITY || z > z_anchor + z_thr,
                    z_thr,
                )
            } else {
                (z_anchor, false, 0.0)
            };
            if dtrace {
                eprintln!(
                    "dtrace Z it={iters} z={z:.9e} anchor={z_anchor:.9e} V={vtotal:.4e} dz={:.3e} thr={z_thr:.3e}",
                    z - z_anchor
                );
            }
            if z_gain {
                z_anchor = z;
            }
            let mut z_net = false;
            if znet_lag > 0 && need_z {
                // Keep the newest entry that is ≥ lag visits old at the front.
                while znet_hist.len() >= 2 && znet_hist[1].0 + znet_lag <= iters {
                    znet_hist.pop_front();
                }
                if let Some(&(_, z_old)) = znet_hist.front()
                    && !z_gain
                    && z > z_old + z_thr
                {
                    z_net = true;
                    n_znet += 1;
                }
                znet_hist.push_back((iters, z));
            }
            if vtotal < 0.9 * v_anchor || z_gain || z_net {
                v_anchor = vtotal;
                no_v_drop = 0;
                stall = 0;
                bland_bland = false;
            } else {
                no_v_drop += 1;
            }
            // Stale-anchor release: a numerical spike (dano3_3: V 1e7, z
            // jumping to ~the optimum at it 1250, then settling back to 288)
            // leaves z_anchor above any honest z, so real progress never
            // counts as a gain and the ladder walks to distress. z is
            // nondecreasing only along a dual feasible path; if it sits well
            // below the anchor once the ladder is mid-way, the anchor was a
            // spike. Re-anchor (bounded, so a z-oscillation cycle still
            // escalates).
            if stall.max(no_v_drop) >= 4 * config.bland_after
                && z_reanchors < 64
                && z < z_anchor - 1e-6 * (1.0 + z_anchor.abs())
            {
                z_reanchors += 1;
                z_anchor = z;
                v_anchor = vtotal;
                no_v_drop = 0;
                stall = 0;
                bland_bland = false;
                srnd = false;
            }
            let esc = stall.max(no_v_drop);
            if incd && hand_over && esc >= 2 * config.bland_after {
                // Rung (2x bland_after) of the ladder under the incremental d: switch
                // this solve to the fresh-pi path, whose ladder behaviour
                // is the validated one.
                if debug {
                    eprintln!("rusolve dual: incremental-d hand-over at it={iters} (esc={esc})");
                }
                incd = false;
                dj_valid = false;
                inv_cnt[0] += 1; // Census: dj invalidation cause
                z_valid = false; // dj stops being maintained ⇒ z with it
                if incp {
                    anchor.event(); // hand over on a freshly anchored π
                }
            }
            if !duph1 && esc >= config.bland_after {
                duph1 = true;
                bland_triggered = true;
                if debug {
                    eprintln!(
                        "rusolve dual: DuPh1 engaged at it={iters} (stall={stall}, V={vtotal:.6e}, esc={esc})"
                    );
                }
            }
            // The cost-perturbation rung engages one ladder step AFTER
            // DuPh1 (esc ≥ 2·bland_after): a benign plateau that tie
            // scoring alone can walk out of gets that chance first (scrs8's
            // cold-dual path costs +485 iters when the tilt lands at 25 —
            // measured by a landing-point sweep); a plateau
            // DuPh1 is blind to — an EXACT θ=0 cycle, where every selection
            // rule re-indexes the same frozen dual point's transition graph
            // — gets the tilt while the warm budget still fits the landing
            // (market_split: engaged at esc=50, Optimal by it≈53 < 100).
            if perturb_on && !pert_used && esc >= 2 * config.bland_after {
                engage_cost_perturbation(st, config);
                dj_valid = false;
                inv_cnt[1] += 1; // Census: dj invalidation cause
                z_valid = false; // true costs moved under z: reseed with the new dj
                if incp {
                    anchor.event(); // effective costs moved: re-anchor π
                }
                pert_used = true;
                bland_triggered = true;
                if debug {
                    eprintln!(
                        "rusolve dual: cost perturbation engaged at it={iters} (stall={stall}, V={vtotal:.6e}, esc={esc})"
                    );
                }
            }
            if duph1 {
                duph1_iters += 1;
                if esc >= 4 * config.bland_after {
                    srnd = true;
                }
            }
            // Bland: the deterministic termination rung (sticky through the
            // episode, releases on a V low with the rest of the ladder).
            if esc >= 8 * config.bland_after && !bland_bland {
                bland_bland = true;
                srnd = false;
                sticky_v = None;
                eprintln!("rusolve dual: Bland engaged after {esc} stalled/plateaued visits");
            }
        } else {
            // Primal feasible: back to classical selection (handoff;
            // the optimality audit decides optimal-vs-continue from here).
            duph1 = false;
            srnd = false;
            sticky_v = None;
            no_v_drop = 0;
            v_anchor = f64::INFINITY;
        }
        if let Some(mk) = tmark {
            acc!(t_chuzr, mk);
        }
        // t_sel: sticky/audit/π-BTRAN/computeDual region between chuzr and
        // the β BTRAN — the re-anchor work a refactorization triggers.
        let tmark_a = if stats {
            Some(std::time::Instant::now())
        } else {
            None
        };
        if srnd && !use_bland && hz_hit {
            // srnd engaged on this visit after a hyper CHUZR, which
            // keeps no top-3 — rebuild it exactly as the 0..m scan would.
            for i in 0..m {
                let v = violation(ub_b[i], x_b[i]);
                if v > FEAS_TOL && v > top3[2].0 {
                    top3[2] = (v, i);
                    top3.sort_by(|a, b| b.0.total_cmp(&a.0));
                }
            }
        }
        if srnd && !use_bland {
            // Sticky-random top-3: keep the sticky row while it violates,
            // else draw a fresh one from the current top-3 violators.
            let sticky_pos = sticky_v.and_then(|sv| {
                (0..m).find(|&i| basis[i] == sv && violation(st.ub[basis[i]], x_b[i]) > FEAS_TOL)
            });
            ppos = match sticky_pos {
                Some(pos) => pos,
                None => {
                    let n = top3.iter().filter(|&&(_, i)| i != usize::MAX).count();
                    let (_, pos) = top3[rng.next_u64() as usize % n.max(1)];
                    sticky_v = Some(basis[pos]);
                    pos
                }
            };
            vbest = violation(st.ub[basis[ppos]], x_b[ppos]);
        }
        if vbest <= FEAS_TOL {
            // No bound violation ⇒ optimality check. The audit can
            // move x_b/basis wholesale ⇒ the list and the bound cache are
            // rebuilt next visit.
            infeas_ok = false;
            ub_b_dirty = true;
            in_basis_dirty = true; // The audit may pivot (primal steps)
            xclamp.wholesale(); // ... and rewrite x_b wholesale
            let outcome = audit(
                st,
                config,
                &mut pfi,
                &mut basis,
                &mut x_b,
                &mut at_ub,
                &in_basis,
                delta.dense_mut(),
                &mut cleanup_done,
                &mut unbounded_verified,
                &mut dse_w,
                &mut iters,
                dtrace,
                &mut acache,
                rowmat.as_ref(),
                passes,
            )?;
            match outcome {
                Audit::Optimal => {
                    let x = extract_x(st, lp, &basis, &x_b, &at_ub);
                    let objective: f64 = lp.objective.iter().zip(&x).map(|(c, xv)| c * xv).sum();
                    // The caller-space dual pair at the final basis, in
                    // the SAME convention as SimplexSolution::reduced_costs
                    // (the single identity d_j = c_j − πᵀa_j in the caller's
                    // own rows/columns/sense). The dual standard form has NO
                    // row negation — its row space IS the caller's — so the
                    // only mapping is the objective sign: the engine
                    // minimises sign·c, hence π_caller = sign·π_std. The optimality
                    // audit certifies Optimal only at δ = 0, so this BTRAN at
                    // the true (unshifted) basic costs IS the verified dual
                    // point the audit priced at.
                    let mut cb = vec![0.0; m];
                    for (i, &v) in basis.iter().enumerate() {
                        cb[i] = st.cost[v] + st.s[v]; // basic shifts are 0
                    }
                    let mut pi = vec![0.0; m];
                    pfi.solve_bt(&cb, &mut pi);
                    let sign = if lp.maximize { -1.0 } else { 1.0 };
                    let row_duals: Vec<f64> = pi.iter().map(|&p| sign * p).collect();
                    let reduced_costs = crate::simplex::reduced_costs_from(lp, &row_duals);
                    if debug {
                        final_probe(st, lp, &basis, &x_b, &at_ub);
                    }
                    let snapshot = export_basis_dual(lp, st, &basis, &at_ub);
                    if stats {
                        eprintln!(
                            "rusolve dual stats: m={m} total={total} nnz={} iters={iters} \
                             degen={n_degen} refactors={} eta_file={} price_calls={n_price} \
                             price_pattern_mean={:.0} ratio_sorts={n_ratio_sort} \
                             cand_mean={:.0} cand_max={cand_max} \
                             scan_ok={} scan_fallback={} scan_groups_max={} \
                             chuzr_list_mean={:.0} chuzr_list_scans={} chuzr_full_scans={}",
                            st.csc_row.len(),
                            pfi.factor_epoch(),
                            pfi.eta_count(),
                            price_pat_sum as f64 / n_price.max(1) as f64,
                            cand_sum as f64 / n_ratio_sort.max(1) as f64,
                            n_scan,
                            n_scan_fb,
                            scan_groups_max,
                            infeas_len_sum as f64 / n_listscan.max(1) as f64,
                            n_listscan,
                            n_fullscan,
                        );
                        eprintln!(
                            "rusolve dual accuracy: pivots={n_acc} alpha_gt1e-7={n_acc7} \
                             alpha_gt1e-9={n_acc9} alpha_max={acc_max:.2e} \
                             dsew_gt0.5={n_dsew_bad} dsew_max={dsew_max:.2e} \
                             xdrift_n={n_xdrift} xdrift_gt1e-6={n_xdrift6} \
                             xdrift_max={xdrift_max:.2e} acc={acc} trouble={n_trouble} \
                             revealed={n_revealed} revealed_max={revealed_max:.2e} \
                             refine_on={refine_on}"
                        );
                        eprintln!("{}", crate::pfi::pfi_stats_line());
                        eprintln!("{}", crate::lu::lu_fact_line());
                        eprintln!(
                            "rusolve pfi factor: lu_work={} etas_now={} basis_nnz={} m={m}",
                            pfi.lu_solve_work(),
                            pfi.eta_count(),
                            pfi.basis_nonzeros(),
                        );
                        dual_time_line!();
                    }
                    if dse_on && let Some((_, out)) = warm_dse.as_mut() {
                        **out = Some(export_dse(lp, st, &basis, &dse_w, &snapshot));
                    }
                    return Ok(Some(DualExit::Sol(SimplexSolution {
                        status: SimplexStatus::Optimal,
                        objective,
                        x,
                        // DuPh1-mode visits are the dual engine's phase 1.
                        phase1_iterations: duph1_iters,
                        // The audit round that PROVED optimality is not a
                        // pivot (primal engine convention: the no-improver
                        // pass returns without counting itself).
                        phase2_iterations: iters - 1 - duph1_iters,
                        bland_triggered,
                        // Dual engine has no primal phase 1; the dirty-exit marker is
                        // never set on this path.
                        dirty_p1_exit: false,
                        // Presolve-OFF path ⇒ the basis snapshot is the
                        // caller's LP's own optimal basis (export from the
                        // dual standard form; see `export_basis_dual`).
                        basis: Some(snapshot),
                        // Valid caller-space duals at Optimal (rays
                        // stay honest-None here: an Optimal has none).
                        row_duals,
                        reduced_costs,
                        duals_valid: true,
                        primal_ray: None,
                        dual_ray: None,
                    })));
                }
                Audit::Stepped => {
                    dj_valid = false; // the audit priced fresh and may have moved the basis/shifts
                    inv_cnt[2] += 1; // Census: dj invalidation cause
                    if incp {
                        anchor.event(); // basis/shifts moved under π̂: re-anchor
                    }
                    z_valid = false; // ... so z reseeds with the fresh dj
                    if let Some(mk) = tmark_a {
                        acc!(t_sel, mk);
                    }
                    continue 'outer;
                }
                Audit::Unbounded => {
                    return Ok(Some(DualExit::Sol(SimplexSolution {
                        status: SimplexStatus::Unbounded,
                        objective: f64::NEG_INFINITY,
                        x: vec![],
                        phase1_iterations: duph1_iters,
                        phase2_iterations: iters - duph1_iters,
                        bland_triggered,
                        ..Default::default() // No rays from the dual engine yet
                    })));
                }
                Audit::Distress => {
                    // Always logged: silent audit
                    // distress hand off to primal without any trace.
                    eprintln!(
                        "rusolve dual: final primal audit distress at {iters} iters ⇒ primal fallback"
                    );
                    return Ok(None);
                }
            }
        }

        let p = ppos;
        if dtrace {
            eprintln!(
                "dtrace LEAVING it={iters} p={} vbest={vbest:.3e} vp={:.6e}",
                basis[p], x_b[p]
            );
        }
        let sigma = if x_b[p] < 0.0 { 1.0 } else { -1.0 };
        let exit_bound = if sigma > 0.0 { 0.0 } else { st.ub[basis[p]] };

        // ---- (b) BTRANs: π (effective costs) and the leaving row β ----
        // Two refined BTRANs per iteration (the 1-BTRAN estimate
        // presumes incremental dual maintenance, which we do not do: π
        // recomputed from factors is the audit-grade source of d̃).
        // Effective costs include the perturbation while it is live;
        // caller_dual_bound's estimate is stated at TRUE costs and is
        // valid for ANY π (an identity, not a dual-feasibility appeal),
        // so the perturbed π needs no special handling at the cutoff.
        let need_pi = if incp {
            anchor.due(incd && dj_valid)
        } else {
            !incd || !dj_valid
        };
        // Cause of THIS pass's anchor, captured before `anchor()` disarms
        // (event anchors — audit step, perturbation, refactor — legitimately
        // jump the dual point; only cadence anchors measure rank-one drift).
        let anchor_by_event = anchor.pending;
        if stats && need_pi {
            n_pi_bt += 1;
            if incp {
                if anchor.pending {
                    n_anchor_evt += 1;
                } else {
                    n_anchor_cad += 1;
                }
            }
        }
        if need_pi {
            let mark = std::time::Instant::now();
            for (i, &v) in basis.iter().enumerate() {
                cb[i] = st.cost[v] + st.s[v] + st.pert[v]; // basic shifts are 0 by invariant
            }
            if picheck {
                pi_ref.copy_from_slice(&piw); // the maintained point the anchor overwrites
            }
            if refine_on {
                pfi.solve_bt(&cb, &mut piw);
            } else {
                pfi.solve_bt_raw(&cb, &mut piw);
            }
            if warm_stats {
                census_acc(|c| c.t_pi_bt += mark.elapsed().as_nanos() as f64);
            }
            if incp {
                anchor.anchor();
                // The anchored point can differ from the maintained π̂ by the
                // accumulated rank-one rounding debt: re-derive d̃ from it
                // (computeDual below), exactly as today at every anchor.
                dj_valid = false;
            }
        }
        // Objective cutoff (warm reopt): π·b + Σ_{nonbasic, d<0} d·ũ is
        // a valid lower bound (min) on the LP optimum for any π — evaluate
        // it BEFORE the δ-commit and stop as soon as it proves the node
        // cannot beat the incumbent. Sound per evaluation (no monotonicity
        // assumed); the cold path never pays a cent (cutoff == None).
        // With the incremental d (z bookkeeping above) the per-iteration
        // estimate is the maintained z_val — no π BTRAN, no O(nnz) scan.
        // A crossing is only a HINT: the exact `caller_dual_bound` from a
        // fresh-factors π (one BTRAN + one O(nnz) pass, crossing iterations
        // only) re-verifies before the prune, so the decision is always
        // made on the exact bound and an estimate that drifted ahead of it
        // just reseeds and continues. When z is not maintained (incd off /
        // recompute cadence / post-hand-over) piw is the raw factors π here
        // and the exact bound is computed directly — today's semantics.
        if let Some(cut) = cutoff {
            let mark = std::time::Instant::now();
            let bound = if z_valid {
                z_val
            } else {
                caller_dual_bound(lp, st, &piw, &in_basis)
            };
            if warm_stats {
                census_acc(|c| c.t_bound += mark.elapsed().as_nanos() as f64);
            }
            if z_valid && zcheck {
                // Reference = the identity at the MAINTAINED dual point
                // (fresh factors π_v plus the live δ — piw itself goes
                // stale through pivots while need_pi is false, but dj and
                // delta stay exact): recomputing it from scratch isolates
                // the z event bookkeeping from representation staleness.
                for (i, &v) in basis.iter().enumerate() {
                    cb[i] = st.cost[v] + st.s[v] + st.pert[v];
                }
                pfi.solve_bt(&cb, &mut pi_v);
                // π_m = fresh factors π + live δ; d recomputed FROM π_m per
                // at-upper column (dj itself is the thing the incremental
                // path could drift, so it cannot be its own reference).
                let pi: Vec<f64> = pi_v.iter().zip(delta.vals()).map(|(p, d)| p + d).collect();
                let mut zs: f64 = pi.iter().zip(&st.b).map(|(p, b)| p * b).sum();
                for v in 0..total {
                    if !in_basis[v] && at_ub[v] && st.ub[v] != 0.0 {
                        zs += (st.cost[v] - csc_dot(&pi, st, v)) * st.ub[v];
                    }
                }
                let exact = if lp.maximize { konst - zs } else { konst + zs };
                let drift = (z_val - exact).abs();
                // dj-staleness probe: if dj itself drifted from the fresh
                // π_m duals, z inherits that drift through its terms —
                // pre-existing incremental-d behavior (bounded by the
                // refactor refresh cadence), not a z-event bug.
                let dj_drift = (0..total)
                    .filter(|&v| !in_basis[v] && at_ub[v] && st.ub[v] != 0.0)
                    .map(|v| {
                        (dj[v] - st.s[v] - st.pert[v] - (st.cost[v] - csc_dot(&pi, st, v))).abs()
                    })
                    .fold(0.0f64, f64::max);
                let zthr = std::env::var("RUSOLVE_DUAL_ZTHRESH")
                    .ok()
                    .and_then(|v| v.parse::<f64>().ok())
                    .unwrap_or(1e-6);
                if drift > zthr * (1.0 + exact.abs()) {
                    eprintln!(
                        "rusolve dual: zcheck it={iters} z={z_val:.9} state={exact:.9} \
                         drift={drift:.3e} dj_drift={dj_drift:.3e} t_pivot={:?}",
                        z_last_t
                    );
                }
            }
            *last_bound = Some(bound);
            let crossed = if lp.maximize {
                bound <= cut
            } else {
                bound >= cut
            };
            if crossed {
                let (exact, pi) = if z_valid {
                    for (i, &v) in basis.iter().enumerate() {
                        cb[i] = st.cost[v] + st.s[v] + st.pert[v];
                    }
                    pfi.solve_bt(&cb, &mut pi_v);
                    (caller_dual_bound(lp, st, &pi_v, &in_basis), pi_v.clone())
                } else {
                    (bound, piw.clone())
                };
                let still = if lp.maximize {
                    exact <= cut
                } else {
                    exact >= cut
                };
                if still {
                    return Ok(Some(DualExit::Pruned {
                        bound: exact,
                        iters,
                        pi,
                    }));
                }
                // Estimate crossed, exact did not: absorb the drift into z
                // and keep solving.
                z_val = exact;
            }
        }
        if need_pi {
            for (w, &dl) in piw.iter_mut().zip(delta.vals()) {
                *w += dl; // the maintained dual point is π + δ
            }
            if picheck {
                // Drift the anchor just erased: maintained π̂ vs the point
                // the re-anchor establishes (fresh factors π + committed δ).
                let mut err = 0.0f64;
                let mut refmax = 0.0f64;
                for (o, &n) in pi_ref.iter().zip(&piw) {
                    err = err.max((o - n).abs());
                    refmax = refmax.max(n.abs());
                }
                let rel = err / (1.0 + refmax);
                if anchor_by_event {
                    // Includes the BTRAN solve error of the fresh factors —
                    // the same noise the per-iteration re-anchor path lives
                    // on; recorded separately, it is the reference level.
                    pic_evt_max = pic_evt_max.max(rel);
                } else {
                    pic_max = pic_max.max(rel);
                    pic_sum += rel;
                    pic_n += 1;
                }
            }
        }
        if incd && !dj_valid {
            // computeDual (HEkk::computeDual): d from the factors, with the
            // cost-shift maintenance of the old per-pass scan.
            for v in 0..total {
                if in_basis[v] || st.ub[v] == 0.0 {
                    continue;
                }
                let mut d = st.cost[v] + st.s[v] + st.pert[v] - csc_dot(&piw, st, v);
                if st.s[v] > 0.0 && d < es {
                    st.s[v] += es - d;
                    d = es;
                }
                if correct_dual {
                    // `RUSOLVE_DUAL_CORRECT`: HEkkDual::correctDualInfeasibilities
                    // (HEkkDual.cpp:2388-2500, MIT; `reference/HiGHS`, read
                    // as a reference text), run at every computeDual like
                    // HiGHS's rebuild. A nonbasic dual infeasible by ≥ Td gets
                    // a cost shift that puts it at (1 + rand)·Td on its
                    // feasible side. HiGHS flips boxed columns instead unless
                    // phase 2 is forced; our cold starts are dual feasible by
                    // construction, which is HiGHS's force_phase2
                    // regime, so we always shift and never move the primal.
                    // The optimality audit prices true costs and drops the shifts.
                    let inf = if at_ub[v] { d } else { -d };
                    if inf >= ptol {
                        let r = crate::simplex::hash01((v as u64) ^ ((n_correct as u64) << 24));
                        let nd = if at_ub[v] {
                            -(1.0 + r) * ptol
                        } else {
                            (1.0 + r) * ptol
                        };
                        st.s[v] += nd - d;
                        d = nd;
                        n_correct += 1;
                        n_correct_max = n_correct_max.max(inf);
                    }
                }
                dj[v] = d;
            }
            dj_valid = true;
            delta_row = if delta.vals().iter().any(|&x| x != 0.0) {
                MIXED_ROW
            } else {
                NO_ROW
            };
            // (re)seed the incremental dual objective from the fresh
            // representation — konst ± (πᵀb̃ + Σ d_true·x̃) at the
            // maintained point, O(m + at-upper count).
            if z_on {
                z_val = z_of_state(st, &in_basis, &at_ub, &dj, &piw, konst, lp.maximize);
                z_valid = true;
                ztrace!("seed", 0.0);
            }
        }
        if let Some(mk) = tmark_a {
            acc!(t_sel, mk);
        }
        let tmark = if stats {
            Some(std::time::Instant::now())
        } else {
            None
        };
        // e_p is a unit vector: O(1) swap of the previous 1.0 instead of an
        // O(m) fill per iteration (the hyper-sparse lesson — per-pivot
        // O(m) furniture dominates once the pivots themselves are cheap).
        // An HVec clear is O(pattern) — the previous 1.0 only.
        e_p.clear();
        e_p.set(p, 1.0);
        pfi.solve_bt_hv(&e_p, &mut beta, refine_on);
        // Ascending pattern: PRICE's row-wise accumulation order and the
        // DSE Σβ² below are ascending-index sums.
        beta.sort();
        if let Some(mk) = tmark {
            acc!(t_beta_bt, mk);
        }
        if stats {
            n_beta_solves += 1;
            n_beta_nz += beta.nnz();
        }
        let tmark = if stats {
            Some(std::time::Instant::now())
        } else {
            None
        };
        pbuf.price(st, rowmat.as_ref(), beta.vals(), beta.pattern(), &in_basis);
        if let Some(mk) = tmark {
            acc!(t_price, mk);
        }
        if stats {
            n_price += 1;
            price_pat_sum += pbuf.list.len();
        }

        // ---- (c,d) Dual ratio test for the PINNED leaving row p: the
        // bound-flipping ratio test (BFRT). One scan prices every event
        // candidate; candidates are grouped by ratio; the entering
        // variable is chosen from the FIRST ratio group whose accumulated
        // potential row change Σ ũ·|ā| covers the row violation
        // |x_b[p] − exit_bound| (or from the LAST group when the
        // candidates run out); every candidate at a STRICTLY EARLIER
        // ratio flips — one batch, one FTRAN. Upper-side candidates (at
        // their ũ) are pivot candidates too: entering FROM THE UPPER SIDE
        // (x̃: ũ → l) is the same pivot with a mirrored step.
        //
        // Why: the old rule —
        // flips applied at EVERY ratio below the first entering event,
        // entering = first entering event with the largest |ā| — let a
        // single boxed flip overshoot the pinned row over its far bound
        // (farcross): the batch was kept, the row abandoned, and pg rode
        // an exact 4-visit farcross cycle that cost perturbation could
        // not see (the duals TRAVEL around that cycle; the transitions
        // are flip geometry, not dual ranking). BFRT kills the overshoot
        // by construction: the crossing flip is a pivot candidate, never
        // a flip; the flip prefix moves row p by strictly less than its
        // violation (toward the bound, never past it); and the duals
        // stop at the CHOSEN pivot's ratio, below which every flipped
        // dual has crossed zero — dual feasibility is exact, not widened
        // (our measurements against widened windows still hold).
        // A θ=0 pivot flips NOTHING: ratios clamp at 0, so no group
        // precedes it (HiGHS's `if (workTheta == 0) workCount = 0;`
        // falls out of the grouping). Candidates tied WITH the pivot
        // (inside its ratio window) are left unflipped; their drift
        // ≤ win·|ā| is tolerance-level, the optimality audit's upper-side
        // primal-step arm is the backstop, and the ratio-0 debt-recovery
        // δ-commit path clears it.
        //
        // Why the row must stay pinned: with a collapsed outer re-scan
        // the t27 autopsy showed a flip for row 5 undone by a saturation
        // for row 4 (one column ping-ponging between its bounds, both
        // events θ=0, state exactly repeated forever). Within one pinned
        // row β and the basis are fixed, so each column's sa = σā has a
        // fixed sign and the saturate arm removes each candidate from
        // the event set permanently ⇒ the inner loop terminates.
        //
        // Dual path π(θ) = π − σθβ ⇒ d̃_j(θ) = d̃_j + σθ ā_j.
        // The scan prices at the MAINTAINED dual point π + δ (piw). The
        // pivot re-anchors the duals at π(θ*) and zeroes δ — correct with
        // a flip prefix because flips do not touch the basis: the new
        // basis' BTRAN reproduces π(θ*) (ā_exit = 1, d̃_q(θ*) = 0, and
        // every other basic column is β-orthogonal at position p). A
        // saturate or a non-pivot exit instead COMMITS δ ← δ − σθ*β, so
        // every flipped column stays dual feasible at the next pass.
        // (Without paying the edge travel, flips leave dual debt behind
        // in π and the ratio-0 debt re-entry pongs forever — the
        // case_54 autopsy.) The .max(0.0) ratio clamps are tolerance
        // guards, plus the recovery path after the optimality audit zeroes δ at
        // a clean candidate: debt columns then flip or re-enter at ratio
        // 0 at true π until a pivot re-anchors the duals.
        //
        // Lineage: HiGHS HEkkDualRow::choosePossible/chooseFinal (BFRT,
        // MIT; `reference/HiGHS`, read as a reference text,
        // reimplemented for the pinned-row loop, not translated). The
        // chooseFinalLargeAlpha BACKWARD scan (a ≥10%-of-max-alpha pivot
        // in a LATER group, trading dual travel for a bigger pivot) is
        // NOT ported. TODO(perf): revisit with sweep numbers.
        // First pass of this row's event loop (see the shift
        // maintenance below).
        let mut ratio_first = true;
        'ratio: loop {
            let tmark_r = if stats {
                Some(std::time::Instant::now())
            } else {
                None
            };
            // Candidate: (v, ā, ratio, enters_from_upper). Both event
            // types are dual-zero events; the side only says which bound
            // the column sits at. Entering event: LOWER side with σā < 0,
            // ratio d̃/(−σā). Flip event: UPPER side with σā > 0, ratio
            // −d̃/(σā) — flip it, or pivot it down from its ũ.
            // `st.ub[v] == 0` skips fixed columns: E-row artificials are
            // deleted from pricing by the fixed bound itself.
            cands.clear();
            let hg = hgroup && !use_bland;
            if !incd {
                // cost-shift maintenance (fresh-pi path): every shifted
                // nonbasic column is brought to d >= eps_s before pricing.
                // A later pass of the SAME pinned row only moved π by
                // −σθβ (saturate/flip commits), which leaves csc_dot
                // bit-identical for every column with no entry in a β row;
                // with a structural (row-wise) PRICE those are exactly the
                // columns outside `pbuf.list`, already at d ≥ ε_s from the
                // first pass (only an ulp-level re-bump of a just-bumped
                // column is forgone). s250r10: the per-pass O(nnz) rescan
                // was most of the 72 s ratio region after the hand-over.
                let maintain_one = |v: usize, st: &mut DualStandard| {
                    if in_basis[v] || st.ub[v] == 0.0 || st.s[v] <= 0.0 {
                        return;
                    }
                    let d = st.cost[v] + st.s[v] + st.pert[v] - csc_dot(&piw, st, v);
                    if d < es {
                        st.s[v] += es - d;
                    }
                };
                if ratio_first || !pbuf.structural || !ratio_rescan_skip {
                    for v in 0..total {
                        maintain_one(v, st);
                    }
                } else {
                    for &v in &pbuf.list {
                        maintain_one(v, st);
                    }
                }
            }
            ratio_first = false;
            let pivtol_it = if ta_ladder {
                pivtol.max(match pfi.eta_count() {
                    0..10 => 1e-9,
                    10..20 => 3e-8,
                    _ => 1e-6,
                })
            } else {
                pivtol
            };
            for &v in &pbuf.list {
                if in_basis[v] || st.ub[v] == 0.0 {
                    continue;
                }
                let abar = pbuf.alpha[v];
                let d = if incd {
                    dj[v]
                } else {
                    st.cost[v] + st.s[v] + st.pert[v] - csc_dot(&piw, st, v)
                };
                // (bumped shifted columns sit at d = eps_s up to rounding;
                // the legacy scan clamped to es exactly)
                let d = if st.s[v] > 0.0 && d < es { es } else { d };
                let sa = sigma * abar;
                // HGROUP keeps the signed ratio (dual·move/α, negative for
                // a tolerance-level wrong-signed dual): HiGHS's relaxed
                // groups price `dual + Td` with the sign intact.
                let lo = if hg { f64::NEG_INFINITY } else { 0.0 };
                if !at_ub[v] {
                    if sa < -pivtol_it {
                        cands.push((v, abar, (d / -sa).max(lo), false));
                    }
                } else if sa > pivtol_it {
                    cands.push((v, abar, (-d / sa).max(lo), true));
                }
            }
            if let Some(w) = watch_col {
                let dw = st.cost[w] + st.s[w] + st.pert[w] - csc_dot(&piw, st, w);
                let pos = cands.iter().position(|c| c.0 == w);
                eprintln!(
                    "watch it={iters} col={w} d={dw:.4e} dj={:.4e} at_ub={} basic={} alpha={:.4e} \
                     sigma={sigma} cand={:?} incd={incd} listed={}",
                    dj.get(w).copied().unwrap_or(f64::NAN),
                    at_ub[w],
                    in_basis[w],
                    pbuf.alpha[w],
                    pos.map(|i| cands[i].2),
                    pbuf.list.contains(&w)
                );
            }
            if cands.is_empty() {
                // Event set empty with row p violated: no column can move
                // x_p into its box at the current bounds — Farkas row
                // candidate, valid regardless of the other rows. One
                // fresh-factor verification pass, then the verdict.
                if !infeas_verified {
                    infeas_verified = true;
                    match pfi.refactorize_reason(crate::pfi::RfReason::Verify) {
                        Ok(()) => {
                            xclamp.wholesale(); // x_b recomputed wholesale
                            dual_refresh(st, &mut pfi, &basis, &at_ub, &mut x_b)
                        }
                        Err(_) => return Ok(None), // Singular ⇒ primal fallback
                    }
                    infeas_ok = false; // x_b recomputed wholesale
                    iters -= 1; // a verification pass is not a pivot
                    dj_valid = false;
                    inv_cnt[3] += 1; // Census: dj invalidation cause
                    z_valid = false; // representation refreshed: reseed with it
                    if incp {
                        anchor.event(); // fresh factors: re-anchor π̂
                    }
                    continue 'outer;
                }
                return Ok(Some(DualExit::Sol(SimplexSolution {
                    status: SimplexStatus::Infeasible,
                    objective: f64::NAN,
                    x: vec![],
                    phase1_iterations: duph1_iters,
                    phase2_iterations: iters - duph1_iters,
                    bland_triggered,
                    // Farkas ray from the empty-event-set row, verified
                    // before it leaves (honest None if the
                    // numerical candidate fails verification).
                    dual_ray: dual_farkas(lp, sigma, beta.vals()),
                    ..Default::default() // No duals from the dual engine
                })));
            }
            if stats {
                n_ratio_sort += 1;
                cand_sum += cands.len();
                cand_max = cand_max.max(cands.len());
            }
            // Group walk (HiGHS chooseFinal step 2): the first ratio
            // group whose accumulated change covers the violation, or
            // the last group when the candidates run out. An
            // infinite-ũ candidate makes its group cross (change = inf),
            // so one-sided candidates are never left behind by a flip —
            // the pivot stops at or before their ratio.
            let total_delta = (x_b[p] - exit_bound).abs();
            // Pivot-group SELECTION. The walk only needs the candidates
            // ordered up to the pivot group, but sorting the whole set was
            // the dual engine's top hotspot on dano3_3 (≈18 % of samples
            // inside `slice::sort`, profiled). For large
            // sets pick the group with HiGHS-style work-group passes over
            // the unsorted vector instead (HEkkDualRow::chooseFinalWorkGroupQuad,
            // MIT; `reference/HiGHS`, read as a reference text — our group
            // semantics stay the exact 1e-12-window ones, not their
            // relaxed-ratio groups): each pass is a min-ratio scan of the
            // ungrouped tail plus a partition swapping that group into the
            // vector prefix, so a walk of depth k costs ~2k linear scans
            // instead of an O(n log n) sort of the whole set. A walk deeper
            // than the cap falls back to the old sort+walk untouched
            // (worst case: a few wasted scans); the selected group is the
            // same SET either way — only the flip-prefix summation order
            // can differ.
            const SCAN_MIN_CANDS: usize = 64; // below: a sort is already cheap
            let mut gq = usize::MAX; // pivot-group start
            let mut gend = 0usize; // pivot-group end (exclusive)
            gstarts.clear(); // group starts (LARGE_ALPHA back-off / HGROUP)
            if hg {
                // HGROUP: HiGHS HEkkDualRow::chooseFinalWorkGroupQuad +
                // chooseFinalLargeAlpha (HEkkDualRow.cpp:317-373, 495-525;
                // MIT, `reference/HiGHS`, read as a reference text and
                // reimplemented on our candidate vector). Group k holds
                // every remaining candidate with ratio <= theta_k, where
                // theta_0 = min (d_j + Td)/alpha_j over all candidates and
                // theta_{k+1} is the same Harris bound over the candidates
                // left after group k, so a group is never empty and every
                // member of an earlier group has crossed zero at any later
                // group's ratio (flips stay dual feasible). Groups
                // accumulate u*alpha until the violation is covered. The
                // pivot group is then the LAST group whose largest alpha
                // exceeds min(0.1 * max alpha, 1), walking backwards: a
                // later group with a healthy pivot beats the break group's
                // tiny one. Unlike the exact-window walk, a candidate tied
                // with the pivot inside Td ends dual infeasible by <= Td;
                // the zero-step shift (`RUSOLVE_DUAL_SHIFT0`) absorbs a
                // wrong-signed pivot dual exactly as HiGHS's shiftCost.
                let mut cum = 0.0f64;
                let mut sel = f64::INFINITY;
                for c in cands.iter() {
                    sel = sel.min(c.2 + ptol / c.1.abs());
                }
                let mut cur = 0usize;
                loop {
                    gstarts.push(cur);
                    let mut remain = f64::INFINITY;
                    let mut j = cur;
                    for k in cur..cands.len() {
                        let a = cands[k].1.abs();
                        if cands[k].2 <= sel {
                            cum += st.ub[cands[k].0] * a;
                            cands.swap(j, k);
                            j += 1;
                        } else {
                            remain = remain.min(cands[k].2 + ptol / a);
                        }
                    }
                    cur = j;
                    sel = remain;
                    if cum >= total_delta || cur == cands.len() {
                        break;
                    }
                }
                gstarts.push(cur);
                let amax = cands[..cur].iter().fold(0.0f64, |a, c| a.max(c.1.abs()));
                let final_cmp = (0.1 * amax).min(1.0);
                for g in (0..gstarts.len() - 1).rev() {
                    let (a, b) = (gstarts[g], gstarts[g + 1]);
                    if cands[a..b].iter().any(|c| c.1.abs() > final_cmp) {
                        gq = a;
                        gend = b;
                        if g + 2 < gstarts.len() {
                            n_large_alpha += 1;
                        }
                        break;
                    }
                }
                if gq == usize::MAX {
                    gq = gstarts[gstarts.len() - 2];
                    gend = cur;
                }
            }
            const SCAN_MAX_GROUPS: usize = 4; // depth cap before sort fallback
            if !hg && cands.len() >= SCAN_MIN_CANDS {
                let mut cum = 0.0f64;
                let mut head = 0usize; // end of the settled pre-pivot groups
                let mut cur = 0usize; // end of the current group
                let mut scan_groups = 0usize;
                loop {
                    let mut g0 = f64::INFINITY;
                    for c in &cands[cur..] {
                        if c.2 < g0 {
                            g0 = c.2;
                        }
                    }
                    let win = 1e-12 * (1.0 + g0);
                    if large_alpha {
                        gstarts.push(cur);
                    }
                    let mut j = cur;
                    for k in cur..cands.len() {
                        if cands[k].2 <= g0 + win {
                            cum += st.ub[cands[k].0] * cands[k].1.abs();
                            cands.swap(j, k);
                            j += 1;
                        }
                    }
                    cur = j;
                    scan_groups += 1;
                    if cum >= total_delta || cur == cands.len() {
                        gq = head;
                        gend = cur;
                        break;
                    }
                    head = cur;
                    if scan_groups >= SCAN_MAX_GROUPS {
                        break; // deep walk: hand the tail to the sort path
                    }
                }
                if stats {
                    scan_groups_max = scan_groups_max.max(scan_groups);
                }
            }
            if hg {
                // grouped above
            } else if gq == usize::MAX {
                if cands.len() >= SCAN_MIN_CANDS {
                    n_scan_fb += 1; // the scan gave up on a deep walk
                }
                cands.sort_by(|a, b| a.2.total_cmp(&b.2).then(a.0.cmp(&b.0)));
                gstarts.clear();
                let mut cum = 0.0f64;
                let mut i = 0usize;
                while i < cands.len() {
                    let g0 = cands[i].2;
                    let win = 1e-12 * (1.0 + g0);
                    if large_alpha {
                        gstarts.push(i);
                    }
                    let mut j = i;
                    while j < cands.len() && cands[j].2 <= g0 + win {
                        cum += st.ub[cands[j].0] * cands[j].1.abs();
                        j += 1;
                    }
                    if cum >= total_delta || j == cands.len() {
                        gq = i;
                        gend = j;
                        break;
                    }
                    i = j;
                }
            } else {
                n_scan += 1;
            }
            // Harris relaxation of the pivot group (Maros 2003 ch. 6; HiGHS
            // HEkkDualRow::chooseFinal works with the same relaxed
            // `d + Td` ratios): every candidate whose exact ratio is within
            // the smallest relaxed ratio min (d_j + Td)/|ā_j| joins the
            // group, so the largest-|ā| pick can avoid tiny pivots (the
            // dano3_3 blow-up: |ā| ~ 1e-8 pivots took V to 1e11). Skipped
            // candidates end up with d ≥ −Td, below the audit's tolerance.
            // Gated on a small best pivot: an always-on relaxation leaves
            // the skipped candidates dual infeasible (≤ Td) and, without
            // HiGHS's cost shifting, the debt stalled 5/47 LPs.
            let best_in_group = cands[gq.min(cands.len().saturating_sub(1))..gend.min(cands.len())]
                .iter()
                .fold(0.0f64, |a, c| a.max(c.1.abs()));
            if harris && !hg && !use_bland && gq != usize::MAX && best_in_group < harris_gate {
                let mut tmax = f64::INFINITY;
                for c in &cands[gq..] {
                    tmax = tmax.min(c.2 + ptol / c.1.abs());
                }
                // Partition extension: every tail candidate with ratio
                // ≤ tmax joins the group. On a sorted vector this selects
                // the same set as the old contiguous `while cands[gend]
                // .2 <= tmax` run; within the group the largest-|ā| pick
                // is order-independent.
                let mut j = gend;
                for k in gend..cands.len() {
                    if cands[k].2 <= tmax {
                        cands.swap(j, k);
                        j += 1;
                    }
                }
                gend = j;
            }
            // LARGE_ALPHA back-off (lever above). cands[..gend] is in group
            // order on both selection paths; gstarts holds the starts.
            if large_alpha && !hg && !use_bland && gq != usize::MAX && gstarts.len() > 1 {
                let amax = cands[..gend].iter().fold(0.0f64, |a, c| a.max(c.1.abs()));
                let final_cmp = (0.1 * amax).min(1.0);
                let best_final = cands[gq..gend].iter().fold(0.0f64, |a, c| a.max(c.1.abs()));
                if best_final <= final_cmp {
                    for g in (0..gstarts.len() - 1).rev() {
                        let (a, b) = (gstarts[g], gstarts[g + 1]);
                        if cands[a..b].iter().any(|c| c.1.abs() > final_cmp) {
                            gq = a;
                            gend = b;
                            n_large_alpha += 1;
                            break;
                        }
                    }
                }
            }

            // Pivot choice within the group (HiGHS chooseFinal step 3,
            // largest |alpha|): Harris' dual analogue (Maros ch. 6);
            // under Bland, lowest variable index (anti-cycling pair of
            // the lowest-index leaving rule). Under DuPh1 the group is
            // additionally scored by expected V absorption
            // min(ũ·|ā|, v_p): landing candidates (absorb = v_p) beat
            // saturating ones (absorb = ũ·|ā| < v_p), larger pivots
            // first within a tier — the old entering arm's scoring,
            // carried over to the wider BFRT candidate group.
            let mut q = usize::MAX;
            let mut abar_q = 0.0f64;
            let mut theta = 0.0f64;
            let mut q_upper = false;
            let mut best_absorb = 0.0f64;
            for &(v, ab, r, upper) in &cands[gq..gend] {
                let absorb = (st.ub[v] * ab.abs()).min(vbest);
                let better = if use_bland {
                    q == usize::MAX || v < q
                } else if duph1 {
                    q == usize::MAX
                        || absorb > best_absorb * (1.0 + 1e-12)
                        || (absorb >= best_absorb * (1.0 - 1e-12)
                            && (ab.abs() > abar_q.abs() * (1.0 + tie)
                                || (ab.abs() >= abar_q.abs() * (1.0 - tie) && v < q)))
                } else {
                    q == usize::MAX
                        || ab.abs() > abar_q.abs() * (1.0 + tie)
                        || (ab.abs() >= abar_q.abs() * (1.0 - tie)
                            && if ctie || hg {
                                tie_key(v) < tie_key(q)
                            } else {
                                v < q
                            })
                };
                if better {
                    q = v;
                    abar_q = ab;
                    // HGROUP ratios are signed: theta = max(d_q/alpha_q, 0)
                    // (HEkkDualRow.cpp:242-246); a no-op otherwise.
                    theta = r.max(0.0);
                    q_upper = upper;
                    best_absorb = absorb;
                }
            }
            if q == usize::MAX {
                // Unreachable by construction (the group walk yields a
                // non-empty group); distress-safe.
                return Ok(None);
            }
            if let Some(mk) = tmark_r {
                acc!(t_ratio, mk);
            }

            // Flip prefix: every candidate at a strictly earlier ratio
            // (groups before the pivot's). Each moves row p TOWARD
            // exit_bound by ũ·|σā| (sign-checked both directions) and
            // the group walk guarantees their sum is below the violation
            // — so no flip overshoots the row, and at the dual stop
            // π(θ*) every flipped dual has crossed zero on the side
            // being left. One batch, one FTRAN (the whole-batch form
            // kept from the earlier rule).
            let tmark_f = if stats {
                Some(std::time::Instant::now())
            } else {
                None
            };
            // The leaving row's PRE-event basic value — the dual-step
            // term of the z telescope (t·x_b[p]⁰ / −σθ·x_b[p]⁰) is stated
            // at the state BEFORE this event's flips move x̃_N.
            let xb_p0 = x_b[p];
            ftmp.clear();
            let mut nf = 0usize;
            for &(v, _, _, upper) in &cands[..gq] {
                nf += 1;
                debug_assert!(st.ub[v].is_finite());
                let u = st.ub[v];
                at_ub[v] = !upper;
                if upper {
                    for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                        *ftmp.entry(st.csc_row[k]) -= u * st.csc_val[k]; // f = Σ Δx̃·a with Δ = −ũ
                    }
                } else {
                    for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                        *ftmp.entry(st.csc_row[k]) += u * st.csc_val[k]; // Δ = +ũ
                    }
                }
            }
            if nf > 0 {
                if stats {
                    flip_sum += nf;
                    flip_max = flip_max.max(nf);
                    flip_batches += 1;
                }
                // Sign note: substituting x̃_j: ũ→0 into
                // x_b = B⁻¹(b − Σ_{at ũ} ũ a) makes the basics RISE by
                // B⁻¹(Σ ũ a) — f carries Δ, so "−B⁻¹f" is that change.
                pfi.solve_b_hv(&ftmp, &mut bbar, config.refine_ftran && refine_on);
                xb_flip_update(
                    &mut x_b,
                    &mut bbar,
                    maintain,
                    &mut infeas,
                    &ub_b,
                    &mut xclamp,
                );
                let ubp = st.ub[basis[p]];
                let vrow = violation(ubp, x_b[p]);
                if dtrace {
                    eprintln!(
                        "dtrace flipbatch it={iters} p={} vp={:.3e} vrow={vrow:.3e} th={theta:.3e} nf={nf}",
                        basis[p], x_b[p],
                    );
                }
                // Numerics guards: by construction the prefix can
                // neither fix nor cross row p (cum < total_delta, every
                // contribution toward the bound). If tolerance-level
                // arithmetic nonetheless does, the duals still advanced
                // to π(θ*) — commit δ (no pivot re-anchors them here)
                // and re-pick the row; flip-only activity FEEDS the
                // escalation ladder (a pivot, not a flip, is the
                // progress signal).
                //
                // `RUSOLVE_DUAL_FLIPPIV`: that commit is a δ debt on row p, and the row is no
                // longer violated, so the next pivot lands on another row
                // and the fresh BTRAN drops the debt — every column the
                // prefix flipped is then left at its new bound with its
                // OLD reduced-cost sign, dual infeasible by |d| (irish it
                // 7,496: one flip of α = 384.34 covering a 384.34
                // violation exactly, d = −833 planted at lower). HiGHS
                // never takes a flip-only step: the BFRT always pivots the
                // final-group column, with whatever primal step remains
                // (HEkkDual::updatePrimal after chooseFinal, MIT). With
                // the lever both exits fall through to the pivot, whose
                // basis prices π(θ*) exactly.
                let fall_through = flip_pivot && {
                    let crossed = (if x_b[p] < 0.0 { 1.0 } else { -1.0 }) != sigma;
                    vrow <= FEAS_TOL || crossed
                };
                if fall_through {
                    n_flip_pivot += 1;
                } else if vrow <= FEAS_TOL {
                    x_b[p] = x_b[p].clamp(0.0, ubp);
                    if debt_shift && !cleanup_done {
                        // `RUSOLVE_DUAL_DEBTSHIFT`: the edge debt as a basic
                        // cost shift (see the declaration).
                        let v = basis[p];
                        st.s[v] -= sigma * theta;
                        zshift[v] -= sigma * theta;
                        n_debt_shift += 1;
                        if incp {
                            anchor.event();
                        }
                        piw_edge_step(&mut piw, &beta, sigma * theta);
                    } else {
                        dual_edge_step(&mut delta, &mut piw, &beta, sigma * theta);
                    }
                    if incd {
                        update_dj(&mut dj, st, &pbuf, &in_basis, es, sigma * theta, usize::MAX);
                        if !(debt_shift && !cleanup_done) {
                            delta_row = if delta_row == NO_ROW || delta_row == p {
                                p
                            } else {
                                MIXED_ROW
                            };
                        }
                        if z_valid {
                            let z0 = z_val;
                            z_flip_terms(&mut z_val, &cands, gq, st, &dj);
                            z_val -= sigma * theta * xb_p0;
                            ztrace!("commit", z_val - z0);
                        }
                    }
                    stall += 1;
                    if let Some(mk) = tmark_f {
                        acc!(t_flip, mk);
                    }
                    continue 'outer;
                }
                if !fall_through && (if x_b[p] < 0.0 { 1.0 } else { -1.0 }) != sigma {
                    if dtrace {
                        eprintln!("dtrace farcross it={iters} p={}", basis[p]);
                    }
                    if debt_shift && !cleanup_done {
                        // `RUSOLVE_DUAL_DEBTSHIFT`: the edge debt as a basic
                        // cost shift (see the declaration).
                        let v = basis[p];
                        st.s[v] -= sigma * theta;
                        zshift[v] -= sigma * theta;
                        n_debt_shift += 1;
                        if incp {
                            anchor.event();
                        }
                        piw_edge_step(&mut piw, &beta, sigma * theta);
                    } else {
                        dual_edge_step(&mut delta, &mut piw, &beta, sigma * theta);
                    }
                    if incd {
                        update_dj(&mut dj, st, &pbuf, &in_basis, es, sigma * theta, usize::MAX);
                        if !(debt_shift && !cleanup_done) {
                            delta_row = if delta_row == NO_ROW || delta_row == p {
                                p
                            } else {
                                MIXED_ROW
                            };
                        }
                        if z_valid {
                            let z0 = z_val;
                            z_flip_terms(&mut z_val, &cands, gq, st, &dj);
                            z_val -= sigma * theta * xb_p0;
                            ztrace!("commit", z_val - z0);
                        }
                    }
                    stall += 1;
                    if let Some(mk) = tmark_f {
                        acc!(t_flip, mk);
                    }
                    continue 'outer;
                }
                // Row p still violated with the same side: fall through
                // to the pivot, which pays the edge offset by
                // re-anchoring the duals (δ = 0).
                if let Some(mk) = tmark_f {
                    acc!(t_flip, mk);
                }
            }

            // Landing check: the pivot step must lie in
            // [0, ũ_q + FEAS_TOL]. (The doc writes σ ā_q in the
            // denominator; at σ = −1 that inverts the sign — the pivot
            // equation is authoritative.) Beyond the window q saturates
            // at its own bound: a flip of one by pivot, dual-feasible
            // because its d̃(θ*) = 0 sits on either side.
            let uq = st.ub[q];
            let (mut l, mut step) = if q_upper {
                // Enters from its upper bound: x̃_q: ũ → l. Removing
                // step = ũ − l from the "at ũ" sum moves x_b by
                // +step·b̄; the row equation is x_b[p] + ā_q·step =
                // exit_bound.
                let step_raw = (exit_bound - x_b[p]) / abar_q;
                if !step_raw.is_finite() {
                    return Ok(None);
                }
                (uq - step_raw, step_raw)
            } else {
                // Enters from its lower bound (the classic arm):
                // x̃_q: 0 → l; x_b falls by l·b̄, row equation
                // x_b[p] − ā_q·l = exit_bound.
                let l_raw = (x_b[p] - exit_bound) / abar_q;
                if !l_raw.is_finite() {
                    return Ok(None);
                }
                (l_raw, l_raw)
            };
            if (if q_upper { step } else { l }) > uq + FEAS_TOL {
                // q saturates at its own bound: a flip batch of one
                // (lower→upper for the classic side, upper→lower for
                // the mirrored side), dual-feasible because
                // d̃_q(θ*) = 0 sits on either side.
                at_ub[q] = !q_upper;
                ftmp.clear();
                for k in st.csc_ptr[q]..st.csc_ptr[q + 1] {
                    // f carries Δx̃ = +ũ (join the "at ũ" sum) or
                    // −ũ (leave it).
                    if q_upper {
                        *ftmp.entry(st.csc_row[k]) -= uq * st.csc_val[k];
                    } else {
                        *ftmp.entry(st.csc_row[k]) += uq * st.csc_val[k];
                    }
                }
                pfi.solve_b_hv(&ftmp, &mut bbar, config.refine_ftran && refine_on);
                xb_flip_update(
                    &mut x_b,
                    &mut bbar,
                    maintain,
                    &mut infeas,
                    &ub_b,
                    &mut xclamp,
                );
                // No pivot: pay for the edge travel this pass bought.
                if debt_shift && !cleanup_done {
                    // `RUSOLVE_DUAL_DEBTSHIFT`: the edge debt as a basic
                    // cost shift (see the declaration).
                    let v = basis[p];
                    st.s[v] -= sigma * theta;
                    zshift[v] -= sigma * theta;
                    n_debt_shift += 1;
                    if incp {
                        anchor.event();
                    }
                    piw_edge_step(&mut piw, &beta, sigma * theta);
                } else {
                    dual_edge_step(&mut delta, &mut piw, &beta, sigma * theta);
                }
                if incd {
                    update_dj(&mut dj, st, &pbuf, &in_basis, es, sigma * theta, usize::MAX);
                    if !(debt_shift && !cleanup_done) {
                        delta_row = if delta_row == NO_ROW || delta_row == p {
                            p
                        } else {
                            MIXED_ROW
                        };
                    }
                    if z_valid {
                        let z0 = z_val;
                        z_flip_terms(&mut z_val, &cands, gq, st, &dj);
                        z_val -= sigma * theta * xb_p0;
                        // q's own flip by saturation (skip = MAX ⇒ dj[q] is
                        // post-step): upper-side enters leave ũ (Δx̃ = −ũ).
                        let dt = dj[q] - st.s[q] - st.pert[q];
                        z_val += if q_upper { -dt * uq } else { dt * uq };
                        ztrace!("saturate", z_val - z0);
                    }
                }
                let ubp = st.ub[basis[p]];
                let vrow = violation(ubp, x_b[p]);
                if dtrace {
                    eprintln!(
                        "dtrace saturate it={iters} q={q} p={} vp={:.3e} vrow={vrow:.3e}",
                        basis[p], x_b[p],
                    );
                }
                if vrow <= FEAS_TOL {
                    x_b[p] = x_b[p].clamp(0.0, ubp);
                    stall = stall.saturating_sub(1);
                    continue 'outer;
                }
                if (if x_b[p] < 0.0 { 1.0 } else { -1.0 }) != sigma {
                    continue 'outer; // far-bound crossing: stale ratios
                }
                continue 'ratio; // next event for the SAME row
            }
            // Land inside the box: clamp the within-tolerance overshoot.
            if q_upper {
                step = step.clamp(0.0, uq);
                l = uq - step;
            } else {
                l = l.clamp(0.0, uq);
                step = l;
            }

            // ---- Pivot (mechanics mirror the primal Some(pos) arm) ----
            let tmark = if stats {
                Some(std::time::Instant::now())
            } else {
                None
            };
            // Gather a_q with its pattern (`gather_col`'s values: the
            // CSC entries assigned, everything else 0).
            qcol.clear();
            for k in st.csc_ptr[q]..st.csc_ptr[q + 1] {
                qcol.set(st.csc_row[k], st.csc_val[k]);
            }
            pfi.solve_b_hv(&qcol, &mut bbar, config.refine_ftran && refine_on);
            if let Some(mk) = tmark {
                acc!(t_qcol, mk);
            }
            if stats {
                let ac = bbar[p];
                let rel = (ac - abar_q).abs() / ac.abs().min(abar_q.abs()).max(1e-300);
                n_acc += 1;
                n_acc7 += usize::from(rel > 1e-7);
                n_acc9 += usize::from(rel > 1e-9);
                acc_max = acc_max.max(rel);
            }
            // Accuracy monitor (HEkk::reinvertOnNumericalTrouble,
            // HEkk.cpp:3019-3029: |abs-difference| / min abs, trouble above
            // 1e-7 once updates exist). On trouble the pivot is ABANDONED,
            // as in HiGHS: updateVerify sets rebuild_reason and
            // updateDual/updatePrimal/updatePivots all return early
            // (HEkkDual.cpp:1231-1250, :2055, :2069), so rebuild() runs on
            // the unchanged basis. Committing the pivot first and
            // reinverting after (the first cut of this lever) is what sent
            // dano3_3 to the primal fallback: a row-wise ā_q = −3.7e-9
            // whose column-wise value was −1e-18 got through the update,
            // and the emergency refactor then found the basis singular.
            // The flips already applied stay; the duals advance to π(θ*)
            // exactly as the flip-batch commit guard above does (δ), so the
            // flipped columns stay dual feasible.
            let trouble = acc_parts & 8 != 0 && pfi.eta_count() > 0 && {
                let (ac, ar) = (bbar[p].abs(), abar_q.abs());
                (ac - ar).abs() / ac.min(ar).max(1e-300) > 1e-7
            };
            if trouble {
                if debt_shift && !cleanup_done {
                    // `RUSOLVE_DUAL_DEBTSHIFT`: the edge debt as a basic
                    // cost shift (see the declaration).
                    let v = basis[p];
                    st.s[v] -= sigma * theta;
                    zshift[v] -= sigma * theta;
                    n_debt_shift += 1;
                    if incp {
                        anchor.event();
                    }
                    piw_edge_step(&mut piw, &beta, sigma * theta);
                } else {
                    dual_edge_step(&mut delta, &mut piw, &beta, sigma * theta);
                }
                if pfi.refactorize_reason(crate::pfi::RfReason::Pivot).is_err() {
                    return Ok(None);
                }
                n_trouble += 1;
                // Our analogue of HiGHS raising the Markowitz threshold
                // (HEkk.cpp:3031-3055): refinement re-armed for the rest of
                // this solve.
                refine_on = true;
                dj_valid = false; // computeDual (HEkkDual.cpp:1074) from π + δ
                z_valid = false;
                if incp {
                    anchor.event();
                }
                // computePrimal + infeasibility list (HEkkDual.cpp:1086-1091).
                xclamp.wholesale(); // x_b recomputed wholesale
                let (nr, vr, dz) =
                    acc_compute_primal(st, &mut pfi, &basis, &at_ub, &mut x_b, stats);
                if z_anchor.is_finite() {
                    z_anchor += dz; // same representation as the next z
                }
                n_revealed += nr;
                revealed_max = revealed_max.max(vr);
                infeas_ok = false;
                stall += 1;
                if debug {
                    eprintln!(
                        "rusolve dual: numerical trouble at it={iters} (alpha col {:.3e} row {abar_q:.3e}): pivot abandoned, reinverted, refinement re-armed",
                        bbar[p]
                    );
                }
                continue 'outer;
            }
            // Zero-step cost shift (HiGHS HEkkDualRow::chooseFinal
            // HEkkDualRow.cpp:242–246 + HEkkDual::updateDual
            // HEkkDual.cpp:2072–2076, MIT; read as a reference text): when the
            // entering column's reduced cost at the maintained point has the
            // WRONG sign (or is zero) the dual step is θ_d = 0, and the
            // entering cost is shifted by −d_q so the new basis prices it at
            // exactly zero. Without it the unclamped t = d_q/ā_q walks the
            // duals BACKWARDS along the pivotal row and plants dual
            // infeasibilities |t|·|ā_j| on every row column — the ex10
            // collapse. The shift lives
            // in `s` and is dropped by the audit cleanup like every other —
            // hence none after that one-shot cleanup: a shift added then
            // would survive into the audit's Optimal certificate.
            if zero_step_shift && !cleanup_done {
                let dq = if incd && dj_valid {
                    dj[q]
                } else {
                    st.cost[q] + st.s[q] + st.pert[q]
                        - piw.iter().zip(qcol.vals()).map(|(p, a)| p * a).sum::<f64>()
                };
                let wrong = if q_upper { dq >= 0.0 } else { dq <= 0.0 };
                if wrong && dq != 0.0 {
                    st.s[q] -= dq;
                    n_zero_shift += 1;
                    if zero_shift_back && !cleanup_done {
                        zshift[q] -= dq;
                    }
                    if incd && dj_valid {
                        dj[q] = 0.0;
                    }
                }
            }
            // DSE update (Forrest–Goldfarb 1992 / HiGHS updateDual): τ = B⁻¹ρ
            // with ρ = B⁻ᵀe_p = β on the OLD factors; w_p is recomputed
            // exactly from ρ, the others follow the rank-one recurrence
            // (floored at 1e-4 like HiGHS).
            if dse_on {
                let tmark = if stats {
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                let ap = bbar[p];
                if ap.abs() > 1e-12 {
                    pfi.solve_b_hv(&beta, &mut tau, false);
                    // Σβ² over the ascending pattern (β sorted after its
                    // BTRAN): the dense ascending sum minus +0.0 terms.
                    let wp: f64 = match beta.pattern() {
                        Some(pat) => pat
                            .iter()
                            .map(|&i| beta[i as usize] * beta[i as usize])
                            .sum(),
                        None => beta.vals().iter().map(|b| b * b).sum(),
                    };
                    if stats {
                        let e = (dse_w[p] - wp).abs() / wp.max(1e-300);
                        n_dsew_bad += usize::from(e > 0.5);
                        dsew_max = dsew_max.max(e);
                    }
                    // r = b̄_i/α_p is 0 off b̄'s pattern: the pattern loop
                    // updates exactly the rows the all-row loop did.
                    let mut dse_row = |i: usize| {
                        if i == p {
                            return;
                        }
                        let r = bbar[i] / ap;
                        if r != 0.0 {
                            dse_w[i] = (dse_w[i] - 2.0 * r * tau[i] + r * r * wp).max(1e-4);
                        }
                    };
                    match bbar.pattern() {
                        Some(pat) => pat.iter().for_each(|&i| dse_row(i as usize)),
                        None => (0..m).for_each(&mut dse_row),
                    }
                    dse_w[p] = (wp / (ap * ap)).max(1e-4);
                } else {
                    dse_w.fill(1.0);
                    hz.valid = false; // Every score moved
                }
                if let Some(mk) = tmark {
                    acc!(t_dse, mk);
                }
            }
            if incp {
                // The pivot's rank-one dual step, π̂ ← π̂ + t·β with the
                // SAME unclamped t = d̃_q/ā_q the incremental-d̃ arm applies
                // below (d̃_j -= t·ā_j ⇔ Δπ̂ = t·β). d̃_q is recomputed from
                // the maintained π̂ exactly as the ratio scan priced it
                // (incl. the shifted-column clamp at es); `qcol` is a_q as
                // gathered for the FTRAN, `beta` this iteration's B⁻ᵀe_p.
                // Degenerate pivots (t = 0) leave π̂ untouched, matching the
                // frozen duals on the edge.
                let dq = st.cost[q] + st.s[q] + st.pert[q]
                    - piw.iter().zip(qcol.vals()).map(|(p, a)| p * a).sum::<f64>();
                let dq = if st.s[q] > 0.0 && dq < es { es } else { dq };
                let t_pi = dq / abar_q;
                if t_pi != 0.0 {
                    for (w, &bb) in piw.iter_mut().zip(beta.vals()) {
                        *w += t_pi * bb;
                    }
                }
                anchor.pivot();
            }
            // Lower entry: x̃ rose by step = l, x_b falls by l·b̄.
            // Upper entry: x̃ fell by step, x_b rises by step·b̄.
            let psgn = if q_upper { 1.0 } else { -1.0 };
            xb_pivot_update(
                &mut x_b,
                &mut bbar,
                psgn,
                step,
                maintain,
                &mut infeas,
                &ub_b,
                &mut xclamp,
            );
            let leaving = basis[p];
            basis[p] = q;
            if in_basis_swap.is_some() {
                in_basis_dirty = true; // never expected: two pivots per pass
            }
            in_basis_swap = Some((leaving, q));
            ub_b[p] = st.ub[q]; // the only row whose cached bound changes
            x_b[p] = l; // lands within tolerance of its bound: no violation to add
            if dtrace {
                eprintln!(
                    "dtrace it={iters} leave={leaving}({exit_bound}) q={q} l={l:.6e} th={theta:.3e} ab={abar_q:.3e} drc={:.2e} stall={stall} nc={}/nf={}{}",
                    (abar_q - bbar[p]).abs() / abar_q.abs().min(bbar[p].abs()),
                    cands.len(),
                    nf,
                    if q_upper { "u" } else { "" },
                );
            }
            at_ub[q] = false;
            if !keep_shift {
                if st.s[q] != 0.0 {
                    dj_valid = false; // basic cost changes => pi changes
                    inv_cnt[4] += 1; // Census: dj invalidation cause
                    if incp {
                        anchor.event(); // basic cost moved: π̂'s rank-one path cannot see it
                    }
                }
                st.s[q] = 0.0; // shift released on entering the basis
                if zero_shift_back || debt_shift {
                    zshift[q] = 0.0;
                }
            }
            at_ub[leaving] = sigma < 0.0; // leaves nonbasic at exit_bound
            if delta_row != NO_ROW && delta_row != p {
                dj_valid = false; // δ debt on another row is dropped by the fresh BTRAN
                inv_cnt[5] += 1; // Census: dj invalidation cause
            }
            delta_row = NO_ROW;
            if z_valid && !(incd && dj_valid) {
                z_valid = false; // dj re-derived next pass ⇒ reseed z with it
            }
            if incd && dj_valid {
                // HEkkDual::updateDual (HEkkDual.cpp:2065): theta_d =
                // d_q / alpha_q (UNclamped: the new basis' duals have
                // d_B = 0 exactly), d_j -= theta_d * alpha_j over the
                // pivotal row, d_q = 0, d_leaving = -theta_d.
                let t = dj[q] / abar_q;
                // The z telescope's q-term bookkeeping needs the PRE-step
                // shift/perturbation of the entering column (the effective
                // step zeroes d_eff_q; the true-d residual is exactly
                // (s_q + pert_q)·x̃_q⁰).
                let (sq, pq, xq0) = (st.s[q], st.pert[q], if q_upper { uq } else { 0.0 });
                if dinf_trace {
                    // Diagnostic: sign of d_q before the step (a wrong-
                    // signed d_q makes the unclamped t walk the duals the
                    // wrong way along the row) and the row's infeasibility
                    // after it.
                    let dq_bad = if q_upper { dj[q] > ptol } else { dj[q] < -ptol };
                    let mut nbad = 0usize;
                    let mut worst = 0.0f64;
                    for &v in &pbuf.list {
                        if in_basis[v] || st.ub[v] == 0.0 || v == q {
                            continue;
                        }
                        let d = dj[v] - t * pbuf.alpha[v];
                        let viol = if at_ub[v] { d - ptol } else { -ptol - d };
                        if viol > 0.0 {
                            nbad += 1;
                            worst = worst.max(viol);
                        }
                    }
                    if dq_bad || worst > 1e-6 {
                        eprintln!(
                            "dinf it={iters} q={q} dq={:.3e} ab={abar_q:.3e} t={t:.3e} th={theta:.3e} \
                             bad_after={nbad} worst={worst:.3e} up={q_upper} s_q={:.1e}",
                            dj[q], st.s[q]
                        );
                    }
                }
                update_dj(&mut dj, st, &pbuf, &in_basis, es, -t, q);
                dj[q] = 0.0;
                dj[leaving] = -t;
                if st.s[leaving] > 0.0 && dj[leaving] < es {
                    st.s[leaving] += es - dj[leaving];
                    dj[leaving] = es;
                }
                if z_valid {
                    let z0 = z_val;
                    z_flip_terms(&mut z_val, &cands, gq, st, &dj);
                    z_val += t * xb_p0
                        + (sq + pq) * xq0
                        + (dj[leaving] - st.s[leaving] - st.pert[leaving]) * exit_bound;
                    z_last_t = Some(t);
                    ztrace!("pivot", z_val - z0);
                }
            }
            if zero_shift_back || debt_shift {
                let zs = std::mem::take(&mut zshift[leaving]);
                if zs != 0.0 && !cleanup_done {
                    // Nonbasic now: its cost no longer enters π, only its
                    // own reduced cost and (at upper) its z term move.
                    st.s[leaving] -= zs;
                    if incd && dj_valid {
                        dj[leaving] -= zs;
                    }
                    if z_valid {
                        z_val -= zs * exit_bound;
                    }
                    n_shift_back += 1;
                }
            }
            // The new basis' BTRAN reproduces π(θ*): the edge offset is
            // paid, and with it the whole flip prefix — flips do not
            // touch the basis, d̃_q(θ*) = 0 and every other basic column
            // is β-orthogonal at p, so π(θ*) IS the new basis dual point.
            // (With the row pinned, a pivot can never plant a violated
            // column at the exit side — ā_exit = 1 gives d_exit(θ*) =
            // σθ* with the correct sign — the flip-pong is dead by
            // construction.)
            delta.clear(); // O(touched) instead of an O(m) zero pass
            // Degeneracy signal for the escalation ladder: in DUAL simplex
            // the degenerate pivot is θ*≈0 (duals frozen on the edge), not
            // a small primal landing — the tuff cycle pivots l=8.5/2.5 with
            // θ=0/3e-5 and never tripped an l-only counter, so Bland never
            // engaged and the swap-cycle ran 98k iterations. Upper entry
            // reads its degeneracy from the step (ũ − l), not from l.
            if theta < 1e-11 || (if q_upper { uq - l } else { l }) < 1e-11 {
                stall += 1;
                if stats {
                    n_degen += 1;
                }
            } else {
                stall = stall.saturating_sub(1);
            }
            // The eta is built from the caller's b̄; retry once on
            // fresh factors with a recomputed b̄ (primal engine verbatim).
            let tmark = if stats {
                Some(std::time::Instant::now())
            } else {
                None
            };
            let upd = pfi.update_hv(p, &mut qcol, &mut bbar, 1.0);
            if let Some(mk) = tmark {
                acc!(t_upd, mk);
            }
            let mut singular_fail = false;
            match upd {
                Ok(UpdateStatus::Updated) => {}
                // computeDual at every reinversion bounds the drift of d.
                Ok(UpdateStatus::Refactorized) => {
                    dj_valid = false;
                    inv_cnt[6] += 1; // Census: dj invalidation cause
                    if incp {
                        anchor.event(); // fresh factors: the natural re-anchor point
                    }
                    if stats {
                        // Census: drift of the maintained x_b from the
                        // fresh-factor computePrimal at this reinversion
                        // (fresh factors carry no etas ⇒ the probe solve
                        // charges the refactor clock nothing: path-faithful).
                        let mut xs = vec![0.0; m];
                        dual_refresh(st, &mut pfi, &basis, &at_ub, &mut xs);
                        let mut worst = 0.0f64;
                        for (a, b) in xs.iter().zip(x_b.iter()) {
                            worst = worst.max((a - b).abs() / (1.0 + b.abs()));
                        }
                        n_xdrift += 1;
                        n_xdrift6 += usize::from(worst > 1e-6);
                        xdrift_max = xdrift_max.max(worst);
                    }
                    if acc_primal {
                        // computePrimal (HEkkDual.cpp:1086) + infeasibility
                        // list rebuild (:1090-1091): x_B from the fresh factors.
                        xclamp.wholesale(); // x_b recomputed wholesale
                        let (nr, vr, dz) =
                            acc_compute_primal(st, &mut pfi, &basis, &at_ub, &mut x_b, stats);
                        if z_anchor.is_finite() {
                            z_anchor += dz; // same representation as the next z
                        }
                        n_revealed += nr;
                        revealed_max = revealed_max.max(vr);
                        infeas_ok = false;
                    }
                }
                Err(FactorError::RefactorizeNeeded { .. }) => {
                    if pfi.refactorize_reason(crate::pfi::RfReason::Pivot).is_err() {
                        if !bt.on {
                            return Ok(None);
                        }
                        singular_fail = true;
                    } else {
                        let mut fresh = vec![0.0; m];
                        pfi.solve_b(qcol.vals(), &mut fresh);
                        if pfi.update(p, qcol.vals(), &fresh, 1.0).is_err() {
                            return Ok(None);
                        }
                        dj_valid = false;
                        inv_cnt[7] += 1; // Census: dj invalidation cause
                        if incp {
                            anchor.event(); // emergency refactor: re-anchor π̂
                        }
                    }
                    if acc_primal {
                        xclamp.wholesale(); // x_b recomputed wholesale
                        let (nr, vr, dz) =
                            acc_compute_primal(st, &mut pfi, &basis, &at_ub, &mut x_b, stats);
                        if z_anchor.is_finite() {
                            z_anchor += dz; // same representation as the next z
                        }
                        n_revealed += nr;
                        revealed_max = revealed_max.max(vr);
                        infeas_ok = false;
                    }
                }
                // Above the dense-fallback cap a Singular auto-refactor
                // is no longer rescued by a dense LU (pfi::dense_fallback_max_n);
                // the eta file would keep retrying a failing sparse factor on
                // every update.
                Err(FactorError::Singular { .. }) if m > crate::pfi::dense_fallback_max_n() => {
                    singular_fail = true;
                }
                Err(_) => {}
            }
            if singular_fail {
                // Backtracking (HiGHS HEkk::getNonsingularInverse,
                // simplex/HEkk.cpp, MIT; read as a reference text): return
                // to the last basis that factorized, with the update limit
                // halved to half the updates that reached the singular one.
                let ups = pfi.eta_count();
                if bt.valid && bt.count < bt.max && ups > 1 {
                    bt.count += 1;
                    bt.eta_limit = (ups / 2).max(1);
                    basis.copy_from_slice(&bt.basis);
                    at_ub.copy_from_slice(&bt.at_ub);
                    st.s.copy_from_slice(&bt.s);
                    dse_w.copy_from_slice(&bt.dse_w);
                    delta.load_dense(&bt.delta);
                    in_basis_dirty = true; // Basis restored wholesale
                    xclamp.wholesale(); // ... and x_b refreshed below
                    let mut cfg = dual_pfi_config(config);
                    cfg.eta_limit = bt.eta_limit.min(cfg.eta_limit);
                    let cols: Vec<Vec<(usize, f64)>> =
                        basis.iter().map(|&v| st.sparse_col(v)).collect();
                    let Ok(np) = BasisFactorization::build_sparse(cols, cfg) else {
                        eprintln!("rusolve dual: backtrack basis singular ⇒ primal fallback");
                        return Ok(None);
                    };
                    pfi = np;
                    dual_refresh(st, &mut pfi, &basis, &at_ub, &mut x_b);
                    for (u, &v) in ub_b.iter_mut().zip(basis.iter()) {
                        *u = st.ub[v];
                    }
                    bt.epoch = pfi.factor_epoch();
                    infeas_ok = false;
                    dj_valid = false;
                    z_valid = false;
                    delta_row = NO_ROW;
                    if incp {
                        anchor.event();
                    }
                    eprintln!(
                        "rusolve dual: singular refactor at it={iters} after {ups} updates ⇒ \
                         backtrack #{} (update limit {})",
                        bt.count, bt.eta_limit
                    );
                    continue 'outer;
                }
                eprintln!(
                    "rusolve dual: singular refactor at it={iters} (m={m}) ⇒ primal fallback"
                );
                return Ok(None);
            }
            unbounded_verified = false; // real pivot re-arms the ray guard
            infeas_verified = false; // ... and so does the infeasible guard
            if debug && iters.is_multiple_of(autopsy_every) {
                dual_autopsy(
                    st,
                    &pfi,
                    &basis,
                    &x_b,
                    &at_ub,
                    &in_basis,
                    delta.vals(),
                    ptol,
                );
            }
            continue 'outer;
        }
    }
}

/// `RUSOLVE_DUAL_AUDIT_INC`: incremental audit pricing lever.
fn audit_inc_on() -> bool {
    parse_audit_inc(std::env::var("RUSOLVE_DUAL_AUDIT_INC").ok().as_deref())
}

fn parse_audit_inc(value: Option<&str>) -> bool {
    match value {
        Some(v) => v.trim() != "0",
        None => AUDIT_INC_DEFAULT,
    }
}

/// Default of the incremental audit pricing lever.
const AUDIT_INC_DEFAULT: bool = false;

/// True-cost reduced costs maintained across consecutive optimality audit
/// primal steps — the primal-simplex cleanup tail after the one-shot shift
/// drop (pds-20: 5.4k steps, 45% of iterations, each a full refined BTRAN
/// of the basic costs + an O(nnz) pricing pass). A primal pivot moves the
/// duals by its pivotal row only: `d_j -= θ_d·α_pj` with `θ_d = d_q/α_pq`,
/// `d_q = 0`, `d_leaving = −θ_d` (HiGHS `HEkkPrimal::updateDual`,
/// HEkkPrimal.cpp:1761, MIT; read as a reference text, reimplemented), so
/// a step costs one BTRAN of `e_p` plus a row-wise PRICE over its pattern.
/// The cache is current only for the loop pass right after its own step
/// (`tag`) and within one factor epoch (fresh pricing at every refactor
/// bounds the update drift, HiGHS's rebuild-time `computeDual`); a clean
/// scan on maintained duals re-prices fresh before any Optimal verdict.
struct AuditCache {
    on: bool,
    valid: bool,
    tag: u64,
    epoch: u64,
    d: Vec<f64>,
    e: Vec<f64>,
    rho: Vec<f64>,
    pbuf: PriceBuf,
    hits: usize,
    fresh: usize,
}

impl AuditCache {
    fn new(total: usize, m: usize, on: bool) -> Self {
        let (nt, nm) = if on { (total, m) } else { (0, 0) };
        AuditCache {
            on,
            valid: false,
            tag: 0,
            epoch: 0,
            d: vec![0.0; nt],
            e: vec![0.0; nm],
            rho: vec![0.0; nm],
            pbuf: PriceBuf::new(nt),
            hits: 0,
            fresh: 0,
        }
    }

    fn usable(&self, pass: u64, epoch: u64) -> bool {
        self.on && self.valid && self.tag == pass && self.epoch == epoch
    }

    /// Pivotal-row dual update for the primal step `q` entering at `pos`
    /// (`leaving` = basis[pos]); `alpha_q` = (B⁻¹a_q)[pos] on the OLD
    /// basis. Must run before the factorization is updated.
    #[allow(clippy::too_many_arguments)]
    fn pivot_update(
        &mut self,
        st: &DualStandard,
        pfi: &BasisFactorization,
        rm: Option<&RowMat>,
        in_basis: &[bool],
        pos: usize,
        q: usize,
        leaving: usize,
        alpha_q: f64,
    ) {
        if rm.is_none() || alpha_q == 0.0 {
            self.valid = false;
            return;
        }
        self.e[pos] = 1.0;
        pfi.solve_bt(&self.e, &mut self.rho);
        self.e[pos] = 0.0;
        self.pbuf.price(st, rm, &self.rho, None, in_basis);
        let theta_d = self.d[q] / alpha_q;
        for &j in &self.pbuf.list {
            if in_basis[j] || st.ub[j] == 0.0 || j == q {
                continue;
            }
            self.d[j] -= theta_d * self.pbuf.alpha[j];
        }
        self.d[q] = 0.0;
        self.d[leaving] = -theta_d;
    }
}

/// Outcome of the optimality check.
enum Audit {
    /// Primal feasible + dual clean (shifts dropped): report Optimal.
    Optimal,
    /// Took a primal step / flip / cleanup drop / δ drop / verification
    /// pass — re-run the outer loop.
    Stepped,
    /// Verified unbounded (primal-step ratio test with no blocker).
    Unbounded,
    /// Singular refactorization: fall back to the primal solver.
    Distress,
}

/// The optimality check: recompute π from the effective costs, price
/// at the MAINTAINED dual point π + δ (the edge offset the flip batches
/// committed), full nonbasic scan of the effective duals d̃ (with cost-shift
/// maintenance), then either a primal step on a truly violated
/// column (the state "primal feasible, truly dual infeasible" is
/// primal-simplex territory), the δ drop at an otherwise-clean candidate
/// (Optimal is only ever certified at δ = 0, i.e. at the true basis
/// duals), the one-time shift cleanup, or the Optimal verdict.
#[allow(clippy::too_many_arguments)]
fn audit(
    st: &mut DualStandard,
    config: &SimplexConfig,
    pfi: &mut BasisFactorization,
    basis: &mut [usize],
    x_b: &mut [f64],
    at_ub: &mut [bool],
    in_basis: &[bool],
    delta: &mut [f64],
    cleanup_done: &mut bool,
    unbounded_verified: &mut bool,
    dse_w: &mut [f64],
    iters: &mut usize,
    dtrace: bool,
    cache: &mut AuditCache,
    rm: Option<&RowMat>,
    pass: u64,
) -> Result<Audit, SimplexError> {
    let m = st.m;
    let total = st.ncol();
    let ptol = config.pricing_tol;
    let es = eps_s(config);
    let mut q = usize::MAX;
    let mut worst = 0.0f64;
    let mut any_shifted = false;
    let delta_zero = delta.iter().all(|&d| d == 0.0);
    // Price on the maintained true-cost duals when the previous loop
    // pass was this audit's own primal step (nothing else moved them). A
    // clean scan there is NOT a certificate: it falls through to the fresh
    // pricing below, so Optimal is only ever declared on fresh duals.
    if cache.usable(pass, pfi.factor_epoch()) && delta_zero {
        for v in 0..total {
            if in_basis[v] || st.ub[v] == 0.0 {
                continue;
            }
            let d = cache.d[v];
            let viol = if at_ub[v] { d - ptol } else { -ptol - d };
            if viol > 0.0 && viol > worst {
                worst = viol;
                q = v;
            }
        }
    }
    if q != usize::MAX {
        cache.hits += 1;
    } else {
        cache.valid = false;
        if cache.on {
            cache.fresh += 1;
        }
        let mut pi = vec![0.0; m];
        let mut cb = vec![0.0; m];
        for (i, &v) in basis.iter().enumerate() {
            // TRUE costs (no cost perturbation): the audit IS the cleanup
            // pass — Optimal is certified only at true-cost KKT, and true
            // violations here take the primal-step arm.
            cb[i] = st.cost[v] + st.s[v];
        }
        pfi.solve_bt(&cb, &mut pi);
        // Price at the MAINTAINED dual point π + δ: flip/saturation batches
        // that did not pivot advanced the duals along the leaving-row edge by
        // δ, and those are the duals whose feasibility Optimal must respect.
        for (w, &dl) in pi.iter_mut().zip(delta.iter()) {
            *w += dl;
        }

        for v in 0..total {
            if in_basis[v] || st.ub[v] == 0.0 {
                if cache.on && in_basis[v] {
                    cache.d[v] = 0.0;
                }
                continue;
            }
            let mut d = st.cost[v] + st.s[v] - csc_dot(&pi, st, v);
            // (zero-step / correction shifts can be negative — any
            // nonzero shift needs the cleanup; the ε_s maintenance stays
            // on s > 0 only.)
            any_shifted |= st.s[v] != 0.0;
            if st.s[v] > 0.0 && d < es {
                st.s[v] += es - d;
                d = es;
            }
            if cache.on {
                cache.d[v] = d;
            }
            // Real dual infeasibility (a shifted column sits at ε_s > 0 after
            // the maintenance bump, so it can never be flagged here).
            let viol = if at_ub[v] { d - ptol } else { -ptol - d };
            if viol > 0.0 && viol > worst {
                worst = viol;
                q = v;
            }
        }

        any_shifted |= basis.iter().any(|&v| st.s[v] != 0.0);
        // Shift-free at δ = 0: the priced d̃ are exactly the true-cost
        // reduced costs a primal step updates by its pivotal row.
        cache.valid = cache.on && !any_shifted && delta_zero;
        cache.epoch = pfi.factor_epoch();
    }
    if q != usize::MAX {
        // ---- Primal step on the violating column q ----
        // Direction from its bound: increase from lower / decrease from
        // upper; the existing bounded PRIMAL ratio test, feasibility side
        // only (Harris two-pass pair). Other columns' duals may transiently
        // worsen — each step strictly improves the objective, and the audit
        // re-runs before every Optimal candidate.
        let sgn = if at_ub[q] { -1.0 } else { 1.0 };
        // gather_col ≡ the old dense clone / negated map: the
        // sgn scale is ±1.0, exact and sign-symmetric.
        let mut aq = vec![0.0; m];
        let mut aq_raw = vec![0.0; m];
        st.gather_col(q, sgn, &mut aq);
        st.gather_col(q, 1.0, &mut aq_raw);
        let mut bbar = vec![0.0; m];
        if config.refine_ftran {
            pfi.solve_b(&aq, &mut bbar);
        } else {
            pfi.solve_b_raw(&aq, &mut bbar);
        }
        let tol = config.pivot_tol;
        let mut theta = st.ub[q];
        let mut leave: Option<usize> = None;
        let mut leave_to_ub = false;
        for i in 0..m {
            let bb = bbar[i];
            let cap = if bb > tol {
                (x_b[i] / bb).max(0.0)
            } else if bb < -tol && st.ub[basis[i]].is_finite() {
                ((st.ub[basis[i]] - x_b[i]) / -bb).max(0.0)
            } else {
                continue;
            };
            if cap < theta {
                theta = cap;
            }
        }
        // Harris pass 2 (near window + the scsd6 lie gate at the clamp
        // budget; no Bland here — the ladder owns escalation).
        let win = 1e-9 * (1.0 + theta);
        let mut best_bar = 0.0f64;
        for i in 0..m {
            let bb = bbar[i];
            let (cap, to_ub) = if bb > tol {
                ((x_b[i] / bb).max(0.0), false)
            } else if bb < -tol && st.ub[basis[i]].is_finite() {
                (((st.ub[basis[i]] - x_b[i]) / -bb).max(0.0), true)
            } else {
                continue;
            };
            if cap > theta + win {
                continue;
            }
            let abar = bb.abs();
            if (cap - theta) * abar > FEAS_TOL {
                continue;
            }
            let better = match leave {
                None => true,
                Some(_) => abar > best_bar * (1.0 + 1e-12),
            };
            if better {
                best_bar = abar;
                leave = Some(i);
                leave_to_ub = to_ub;
            }
        }
        if dtrace {
            eprintln!(
                "dtrace audit q={q} sgn={sgn} theta={theta:.3e} leave={:?} to_ub={leave_to_ub} it={}",
                leave, *iters
            );
        }
        match leave {
            None if !theta.is_finite() => {
                // No blocker and no ub cap in the improving direction:
                // Unbounded candidate — verify once on fresh factors.
                if !*unbounded_verified {
                    *unbounded_verified = true;
                    cache.valid = false;
                    match pfi.refactorize_reason(crate::pfi::RfReason::Verify) {
                        Ok(()) => dual_refresh(st, pfi, basis, at_ub, x_b),
                        Err(_) => return Ok(Audit::Distress),
                    }
                    *iters -= 1; // verification is not a pivot
                    return Ok(Audit::Stepped);
                }
                return Ok(Audit::Unbounded);
            }
            None => {
                // q saturates at its own ub: flip, no basis change. It
                // TRAVERSED its full bound span, so it lands on the far
                // side (primal engine convention `at_ub[q] = !from_ub`).
                for i in 0..m {
                    x_b[i] -= theta * bbar[i];
                    if x_b[i] < 0.0 && x_b[i] > -FEAS_TOL {
                        x_b[i] = 0.0;
                    }
                }
                at_ub[q] = sgn > 0.0;
                // δ STAYS: no basis change, and the flip moved q onto its
                // dual-feasible side at π + δ (the violation sides of a
                // column at its lower vs upper bound are opposite).
            }
            Some(pos) => {
                let ep0 = pfi.factor_epoch();
                if cache.valid {
                    cache.pivot_update(st, pfi, rm, in_basis, pos, q, basis[pos], sgn * bbar[pos]);
                }
                for i in 0..m {
                    x_b[i] -= theta * bbar[i];
                    if x_b[i] < 0.0 && x_b[i] > -FEAS_TOL {
                        x_b[i] = 0.0;
                    }
                }
                dse_w.fill(1.0); // basis changed outside the DSE recurrence
                let leaving = basis[pos];
                basis[pos] = q;
                x_b[pos] = if sgn < 0.0 { st.ub[q] - theta } else { theta };
                at_ub[q] = false;
                st.s[q] = 0.0;
                at_ub[leaving] = leave_to_ub;
                // Real pivot ⇒ basis changes ⇒ π recomputed from factors is
                // authoritative; the flip-edge offset δ is subsumed.
                for dv in delta.iter_mut() {
                    *dv = 0.0;
                }
                match pfi.update(pos, &aq_raw, &bbar, sgn) {
                    Err(FactorError::RefactorizeNeeded { .. }) => {
                        if pfi.refactorize_reason(crate::pfi::RfReason::Pivot).is_err() {
                            return Ok(Audit::Distress);
                        }
                        let mut fresh = vec![0.0; m];
                        pfi.solve_b(&aq_raw, &mut fresh);
                        if pfi.update(pos, &aq_raw, &fresh, sgn).is_err() {
                            return Ok(Audit::Distress);
                        }
                    }
                    // See the main loop's update site.
                    Err(FactorError::Singular { .. }) if m > crate::pfi::dense_fallback_max_n() => {
                        return Ok(Audit::Distress);
                    }
                    _ => {}
                }
                if pfi.factor_epoch() != ep0 {
                    cache.valid = false; // fresh factors: re-price fresh (drift reset)
                }
            }
        }
        cache.tag = pass + 1;
        *unbounded_verified = false;
        return Ok(Audit::Stepped);
    }

    if delta.iter().any(|&d| d != 0.0) {
        // Otherwise clean, but flip batches have advanced the maintained
        // duals along an edge NO PIVOT has ever certified. The Optimal
        // certificate (complementary slackness) is stated at the TRUE basis
        // duals π, so drop the offset and re-audit: any real violation at π
        // then takes the primal-step arm (which pivots) or the ratio-0
        // scan path. Optimal is only ever returned at δ = 0.
        for d in delta.iter_mut() {
            *d = 0.0;
        }
        cache.valid = false;
        *iters -= 1; // an offset drop is not a pivot either
        return Ok(Audit::Stepped);
    }

    if any_shifted && !*cleanup_done {
        // Cleanup (at most ONE full drop): drop every remaining
        // shift and re-audit at TRUE costs; real violations then take the
        // primal-step arm above. (Maintenance bumps only fire on s > 0, so
        // after this drop nothing re-shifts: the cleanup is one-shot.)
        *cleanup_done = true;
        if dtrace {
            let ns = st.s.iter().filter(|&&x| x > 0.0).count();
            eprintln!("dtrace cleanup-drop shifted={ns} it={}", *iters);
        }
        for sv in st.s.iter_mut() {
            *sv = 0.0;
        }
        *iters -= 1; // a cleanup drop is not a pivot either
        return Ok(Audit::Stepped);
    }

    // Optimal: final re-anchor on fresh factors (primal tail convention; if
    // the refactorization is singular the eta chain is the better solver).
    if pfi.refactorize_reason(crate::pfi::RfReason::Tail).is_ok() {
        dual_refresh(st, pfi, basis, at_ub, x_b);
    }
    Ok(Audit::Optimal)
}

/// Cost perturbation: the escalation ladder's answer to an EXACT
/// θ=0 cycle. Such a cycle is a property of the FROZEN dual point — every
/// selection rule (max-violation leaving, largest-|ā| entering, DuPh1 tie
/// scoring, sticky-random top-3) re-indexes the same transition graph, so
/// none can break it; only moving π off the degenerate face can. The tilt
/// is hash-signed and strictly SUB-tolerance:
///
/// ```text
/// pert_v = base_v · (2·hash01(v ⊕ 0xD17) − 1)
/// base_v = eps_s · (1 + |cost_v|)     (structural / slack columns)
/// base_v = rowbase = 1e-12            (E-row fixed artificials — the
///                                      HiGHS row constant, HEkk::
///                                      initialiseCost;
///                                      perturbing a basic artificial's
///                                      cost is a perturbed equality-row
///                                      dual, the tie-break for degenerate
///                                      equality pivots)
/// ```
///
/// Lineage: Chvátal §7.3 / Beale 1968 (perturbation against degenerate
/// cycling); the hash01 finalizer is the
/// primal engine's (deterministic, draw-order-independent). A flat
/// `colbase = 1e-7` would equal our `pricing_tol`, so the rung
/// uses `eps_s = ptol/10` — strictly sub-tolerance, the criterion that
/// the true-cost optimality audit (which prices `cost + s`, never `pert`) can
/// never disagree with the perturbed ranking by a real amount. An
/// incidence-calibrated colbase (HiGHS logs a column base of 5.63e-8 on
/// tuff and 1.69e-6 on stair) remains open — TODO(perf): calibrate on the
/// degenerate LPs.
///
/// Honesty: the perturbation NEVER enters the optimality audit, the
/// reported dual pair, or `caller_dual_bound`'s d-vector (all price TRUE
/// costs), so it can only change WHICH feasible-optimal vertex is found,
/// never the verdict. (The cutoff bound consumes the perturbed π, which
/// weak-duality allows: the bound stays valid, only marginally looser.)
fn engage_cost_perturbation(st: &mut DualStandard, config: &SimplexConfig) {
    let colbase = eps_s(config);
    let rowbase = 1e-12; // the HiGHS row constant (HEkk::initialiseCost)
    for (v, p) in st.pert.iter_mut().enumerate() {
        let base = if st.kind[v] == Kind::Artificial {
            rowbase
        } else {
            colbase * (1.0 + st.cost[v].abs())
        };
        *p = base * (2.0 * crate::simplex::hash01((v as u64) ^ 0xD17) - 1.0);
    }
}

/// `RUSOLVE_DUAL_PERT_SEED` (see [`highs_cost_perturbation`]); 0 when unset.
fn pert_seed() -> u64 {
    use std::sync::OnceLock;
    static SEED: OnceLock<u64> = OnceLock::new();
    *SEED.get_or_init(|| {
        std::env::var("RUSOLVE_DUAL_PERT_SEED")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    })
}

/// `RUSOLVE_DUAL_ZFLOOR=K`: the z-gain noise-floor multiplier of
/// the plateau ladder; `0` disables the floor (the plain relative test).
fn zgain_floor_k() -> f64 {
    parse_zgain_floor(std::env::var("RUSOLVE_DUAL_ZFLOOR").ok().as_deref())
}

fn parse_zgain_floor(value: Option<&str>) -> f64 {
    value
        .and_then(|v| v.trim().parse::<f64>().ok())
        .filter(|k| k.is_finite() && *k >= 0.0)
        .unwrap_or(ZGAIN_FLOOR_K)
}

/// `RUSOLVE_DUAL_ZLAZY`: `1` = lazy ladder z-pass on large LPs
/// (m + total ≥ ZLAZY_MIN_DIM), `0` = per-visit pass. Lazy visits skip the
/// pass while esc + 1 < bland_after / 2.
fn zlazy_at(bland_after: usize, big: bool) -> usize {
    let on = match std::env::var("RUSOLVE_DUAL_ZLAZY").ok().as_deref() {
        Some(v) => v.trim() != "0",
        None => ZLAZY_DEFAULT,
    };
    if on && big {
        (bland_after / 2).max(1)
    } else {
        0
    }
}

/// The ladder's dual objective z = Σ c̃_B·x_B + Σ_{nonbasic at ũ} c̃·ũ
/// (c̃ = cost + shift + perturbation) and zabs = Σ|terms|, summed in the
/// order (basic rows ascending, then variables ascending).
///
/// Hyper-sparse passes: bit-identical to the plain
/// two loops while touching less. The at-ũ loop skips 64-variable blocks
/// with no at_ub flag through a branch-free OR fold: on stormG2 that is a
/// few thousand qualifying variables out of 1.5M (0.45 → 0.06 ms per pass).
/// A skipped variable adds no term, so the ascending order of the
/// remaining terms is unchanged.
fn ladder_z(
    st: &DualStandard,
    basis: &[usize],
    x_b: &[f64],
    at_ub: &[bool],
    in_basis: &[bool],
    total: usize,
) -> (f64, f64) {
    let mut z = 0.0f64;
    let mut zabs = 0.0f64;
    for (i, &v) in basis.iter().enumerate() {
        // Unconditional: a `x_b != 0` skip would be bit-identical too, but on
        // stormG2 (~40% nonzero x_b) its mispredicts made the loop 1.5x
        // slower than the plain gathers.
        let t = (st.cost[v] + st.s[v] + st.pert[v]) * x_b[i];
        z += t;
        zabs += t.abs();
    }
    const BLK: usize = 64;
    for (c, blk) in at_ub[..total].chunks(BLK).enumerate() {
        if !blk.iter().fold(false, |a, &b| a | b) {
            continue;
        }
        for (k, &u) in blk.iter().enumerate() {
            let v = c * BLK + k;
            if u && !in_basis[v] && st.ub[v].is_finite() {
                let t = (st.cost[v] + st.s[v] + st.pert[v]) * st.ub[v];
                z += t;
                zabs += t.abs();
            }
        }
    }
    (z, zabs)
}

/// `RUSOLVE_DUAL_HCHUZR`: `1` = hyper CHUZR on LPs with
/// m ≥ HCHUZR_MIN_M, `0` = always the full violation scan, `2` = every
/// size (validation: exercises the path on relaxed MIPLIB LPs / MIP node LPs).
fn hchuzr_on(m: usize) -> bool {
    match std::env::var("RUSOLVE_DUAL_HCHUZR")
        .ok()
        .as_deref()
        .map(str::trim)
    {
        Some("0") => false,
        Some("2") => true,
        Some(_) => m >= HCHUZR_MIN_M,
        None => HCHUZR_DEFAULT && m >= HCHUZR_MIN_M,
    }
}

/// Default of the hyper-CHUZR lever.
/// ON since the adaptive fallback below (measured on our benchmark set): stormG2_1000
/// 3.9x iterations in 120 s; supportcase10/square41/nug08-3rd/pds-40
/// iteration-identical with no wall tax.
const HCHUZR_DEFAULT: bool = true;
/// Size gate: below it the O(m) scan is not worth a candidate set.
const HCHUZR_MIN_M: usize = 20_000;
/// Candidate-set size k = m / DIV clamped to [MIN, MAX] (doubled when the
/// set exhausts within k/8 visits).
const HCHUZR_K_DIV: usize = 1000;
const HCHUZR_K_MIN: usize = 256;
const HCHUZR_K_MAX: usize = 4096;
/// Adaptive fallback ([`HyperChuzr::gate`]): minimum running mean of
/// CHUZR hits per candidate-set build, its EMA weight and starting value,
/// and the cool-down range in CHUZR visits.
const HCHUZR_MIN_HITS: f64 = 3.0;
const HCHUZR_EMA_W: f64 = 0.25;
const HCHUZR_EMA_INIT: f64 = 16.0;
const HCHUZR_COOL_MIN: usize = 32;
const HCHUZR_COOL_MAX: usize = 1024;

/// Default of the lazy z-pass lever.
const ZLAZY_DEFAULT: bool = true;

/// Size gate of the lazy z-pass.
/// The pass is O(m + total) per visit: 8-28% of wall on four large-LP
/// targets (m + total 0.28M-1.8M), noise below. Lazy is NOT
/// trajectory-neutral — skipped visits do not reset the V anchor — and on
/// the degenerate qap10 (m + total 7.8k) that alone took the dual from 34k
/// to 73k iterations, so it is confined to sizes where the pass pays.
const ZLAZY_MIN_DIM: usize = 100_000;

/// Default z-gain floor multiplier. See `zgain_threshold`.
const ZGAIN_FLOOR_K: f64 = 1000.0;

/// The ladder's z-gain threshold: the relative 1e-9·(1+|z_anchor|) test,
/// capped (when `k > 0`) at `k·ε·zabs` — `zabs = Σ|terms|` of the z sum
/// bounds its rounding error to O(n·ε·zabs) in the worst case and
/// O(√n·ε·zabs) typically (Higham, *Accuracy and Stability of Numerical
/// Algorithms*, 2nd ed., §4.2), so a gain above k·ε·zabs with k ≫ √n is
/// not rounding.
fn zgain_threshold(z_anchor: f64, zabs: f64, k: f64) -> f64 {
    let rel = 1e-9 * (1.0 + z_anchor.abs());
    if k > 0.0 {
        rel.min(k * f64::EPSILON * zabs)
    } else {
        rel
    }
}

/// Up-front HiGHS-style cost perturbation for COLD dual solves (the missing
/// piece behind the silent primal fallbacks, measured on our benchmark set:
/// dual-degenerate LPs — set covering/partitioning, assignment — stalled the
/// unperturbed engine at θ = 0 until the 400-visit distress rung).
///
/// Source: HiGHS `HEkk::initialiseCost` (MIT, reference/HiGHS/highs/simplex/
/// HEkk.cpp), read as a reference text and reimplemented over our shifted
/// form; lineage Koberstein 2005 (PhD thesis, Paderborn) §4.2 and Maros 2003.
/// Base `5e-7 · m` with `m = max|c|` (`sqrt(sqrt(m))` above 100; capped at 1
/// when fewer than 1% of columns are boxed); each column moves by
/// `(1 + u_j)(|c_j| + 1) · base`, `u_j ∈ [0,1)` hash-random, in the
/// direction that keeps the cost-sign dual start feasible: up for ũ = ∞
/// (our shifted columns all sit at 0 with d ≥ 0), away from zero for boxed
/// columns. Fixed columns and the two halves of a split free variable are
/// not perturbed (HiGHS: "free - no perturb"). Slack / artificial columns
/// get HiGHS's row constant `1e-12 · (0.5 − u)`.
///
/// The perturbation is above `pricing_tol`, so unlike the sub-tolerance
/// escalation rung it can leave true-cost dual infeasibilities at the
/// perturbed optimum; the optimality audit prices TRUE costs and finishes with
/// primal steps (HiGHS: the primal "cleanup" phase). Optimal is therefore
/// still only certified at true-cost KKT.
fn highs_cost_perturbation(st: &mut DualStandard) {
    let ns = st.num_structural;
    let mut max_abs = 0.0f64;
    for v in 0..ns {
        max_abs = max_abs.max(st.cost[v].abs());
    }
    if max_abs == 0.0 {
        // Pure feasibility LP: HiGHS only logs "perturb as if max cost was
        // 1" and then perturbs by 0, leaving every pivot at θ ≈ 0 (krka:
        // 400-visit false stall). Perturb as logged.
        max_abs = 1.0;
    }
    if max_abs > 100.0 {
        max_abs = max_abs.sqrt().sqrt();
    }
    let boxed = (0..ns)
        .filter(|&v| st.ub[v].is_finite() && st.ub[v] > 0.0)
        .count();
    let total = st.ncol().max(1);
    if (boxed as f64) / (total as f64) < 0.01 {
        max_abs = max_abs.min(1.0);
    }
    let base = 5e-7 * max_abs;
    let mut split = vec![false; st.ncol()];
    for cs in &st.cols_of {
        if cs.len() == 2 {
            split[cs[0]] = true;
            split[cs[1]] = true;
        }
    }
    // Chaos-spread instrument (`RUSOLVE_DUAL_PERT_SEED=k`, default 0 =
    // the historic draw, bit-identical): re-seeds the perturbation's random
    // multipliers so a path-sensitivity A/B can be judged against the
    // spread of statistically equivalent paths (measured: istanbul 1 of 5
    // seeds 9×, dano3_3 3 of 3 to the limit).
    let seed = pert_seed();
    for v in 0..st.ncol() {
        let u =
            crate::simplex::hash01((v as u64) ^ 0xD17 ^ seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        if v >= ns {
            st.pert[v] = (0.5 - u) * 1e-12;
        } else if split[v] || st.ub[v] == 0.0 {
            // free / fixed: no perturbation
        } else {
            let xp = (1.0 + u) * (st.cost[v].abs() + 1.0) * base;
            st.pert[v] = if st.ub[v].is_finite() && st.cost[v] < 0.0 {
                -xp
            } else {
                xp
            };
        }
    }
}

/// Dual-side re-anchor: x_b = B⁻¹(b − Σ_{nonbasic at ũ} ũ a) WITHOUT the
/// primal round-off clamp — x_b negatives are the dual engine's work input;
/// only |v| < 1e-12 noise snaps to 0.
fn dual_refresh(
    st: &DualStandard,
    pfi: &mut BasisFactorization,
    basis: &[usize],
    at_ub: &[bool],
    x_b: &mut [f64],
) {
    let mut rhs = st.b.clone();
    let total = st.ncol();
    let mut in_b = vec![false; total];
    for &v in basis.iter() {
        in_b[v] = true;
    }
    for (v, &ib) in in_b.iter().enumerate() {
        if !ib && at_ub[v] && st.ub[v].is_finite() {
            for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                rhs[st.csc_row[k]] -= st.ub[v] * st.csc_val[k];
            }
        }
    }
    pfi.solve_b(&rhs, x_b);
    for xv in x_b.iter_mut() {
        if xv.abs() < 1e-12 {
            *xv = 0.0;
        }
    }
}

/// computePrimal at a reinversion (HEkkDual::rebuild, HEkkDual.cpp:1086):
/// x_B from the fresh factors via [`dual_refresh`]. With `stats`, returns the
/// rows the recompute turned violated (> FEAS_TOL) and their worst violation
/// — the hidden primal debt the incremental x_B had been carrying — and the
/// jump dz of the ladder's progress measure z = c̃ᵀx̃ the recompute caused:
/// the caller rebases z_anchor by it, so a change of REPRESENTATION never
/// reads as lost (or gained) dual progress.
fn acc_compute_primal(
    st: &DualStandard,
    pfi: &mut BasisFactorization,
    basis: &[usize],
    at_ub: &[bool],
    x_b: &mut [f64],
    stats: bool,
) -> (usize, f64, f64) {
    let old = x_b.to_vec();
    dual_refresh(st, pfi, basis, at_ub, x_b);
    let mut dz = 0.0f64;
    for (i, &v) in basis.iter().enumerate() {
        dz += (st.cost[v] + st.s[v] + st.pert[v]) * (x_b[i] - old[i]);
    }
    let (mut n, mut worst) = (0usize, 0.0f64);
    if stats {
        for (i, &v) in basis.iter().enumerate() {
            let ub = st.ub[v];
            let (vo, vn) = (violation(ub, old[i]), violation(ub, x_b[i]));
            if vn > FEAS_TOL && vo <= FEAS_TOL {
                n += 1;
                worst = worst.max(vn);
            }
        }
    }
    (n, worst, dz)
}

/// Iteration-limit snapshot (mirrors the primal `limit_solution`).
#[allow(clippy::too_many_arguments)]
fn limit_dual(
    st: &DualStandard,
    lp: &Lp,
    basis: &[usize],
    x_b: &[f64],
    at_ub: &[bool],
    p1: usize,
    iters: usize,
    bland_triggered: bool,
) -> SimplexSolution {
    let x = extract_x(st, lp, basis, x_b, at_ub);
    let objective: f64 = lp.objective.iter().zip(&x).map(|(c, xv)| c * xv).sum();
    SimplexSolution {
        status: SimplexStatus::IterationLimit,
        objective,
        x,
        phase1_iterations: p1,
        phase2_iterations: iters - p1,
        bland_triggered,
        ..Default::default() // Unfinished search ⇒ no duals, no rays
    }
}

/// Debug audit (every 64th pivot, RUSOLVE_DEBUG): the std-system
/// residual |x_b − B⁻¹(b − Nx_N)| and the dual sign conditions priced at the
/// maintained point π + δ. Sign errors in a dual ratio test are otherwise
/// invisible until they produce a wrong Optimal (an earlier wrong-answer hunt's lesson) —
/// with the δ-commit the maintained point should stay dual feasible between
/// audits, so a warning here means the edge bookkeeping itself is broken.
#[allow(clippy::too_many_arguments)]
fn dual_autopsy(
    st: &DualStandard,
    pfi: &BasisFactorization,
    basis: &[usize],
    x_b: &[f64],
    at_ub: &[bool],
    in_basis: &[bool],
    delta: &[f64],
    ptol: f64,
) {
    let m = st.m;
    let total = st.ncol();
    // The caller's `in_basis` is the loop-head copy: after a pivot the
    // leaving column is still flagged basic there, so an at-ũ exit would
    // drop its ũ·a term and fake a residual. Rebuild it.
    let _ = in_basis;
    let mut in_now = vec![false; total];
    for &v in basis {
        in_now[v] = true;
    }
    let in_basis = &in_now[..];
    let mut r = st.b.clone();
    for (i, &v) in basis.iter().enumerate() {
        for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
            r[st.csc_row[k]] -= x_b[i] * st.csc_val[k];
        }
    }
    for (v, &ib) in in_basis.iter().enumerate().take(total) {
        if !ib && at_ub[v] && st.ub[v].is_finite() {
            for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                r[st.csc_row[k]] -= st.ub[v] * st.csc_val[k];
            }
        }
    }
    let mut y = vec![0.0; m];
    pfi.solve_b(&r, &mut y);
    // r already subtracts the basic contribution, so y = B⁻¹r is the drift
    // of x_b from the std system itself: a consistent state gives y ≈ 0.
    let mut resid = 0.0f64;
    let mut scale = 1.0f64;
    for i in 0..m {
        resid = resid.max(y[i].abs());
        scale = scale.max(x_b[i].abs());
    }
    let mut cb = vec![0.0; m];
    for (i, &v) in basis.iter().enumerate() {
        cb[i] = st.cost[v] + st.s[v] + st.pert[v];
    }
    let mut pi = vec![0.0; m];
    pfi.solve_bt(&cb, &mut pi);
    for (w, &dl) in pi.iter_mut().zip(delta.iter()) {
        *w += dl;
    }
    let mut worst_dual = 0.0f64;
    let mut worst_v = usize::MAX;
    for v in 0..total {
        if in_basis[v] || st.ub[v] == 0.0 {
            continue;
        }
        let d = st.cost[v] + st.s[v] + st.pert[v] - csc_dot(&pi, st, v);
        let viol = if at_ub[v] { d - ptol } else { -ptol - d };
        if viol > worst_dual {
            worst_dual = viol;
            worst_v = v;
        }
    }
    if resid > 1e-6 * scale || worst_dual > 1e-4 {
        eprintln!(
            "rusolve dual AUTOPSY: residual {resid:.3e} (scale {scale:.3e}) \
             worst dual viol {worst_dual:.3e} at col {worst_v}"
        );
        if worst_v != usize::MAX {
            eprintln!(
                "  worst col {worst_v}: cost={:e} s={:e} pert={:e} ub={:e} at_ub={}",
                st.cost[worst_v], st.s[worst_v], st.pert[worst_v], st.ub[worst_v], at_ub[worst_v]
            );
        }
    }
}

/// Terminal debug probe: basic bounds, std-system
/// residual, caller-space rows/bounds, E-row artificial cleanliness, and the
/// effective-vs-true cost audit (every shift must be gone at Optimal).
fn final_probe(st: &DualStandard, lp: &Lp, basis: &[usize], x_b: &[f64], at_ub: &[bool]) {
    let mut worst_basic = (0.0f64, 0usize); // (violation, column)
    for (i, &v) in basis.iter().enumerate() {
        let viol = violation(st.ub[v], x_b[i]);
        if viol > worst_basic.0 {
            worst_basic = (viol, v);
        }
    }
    let mut in_b = vec![false; st.ncol()];
    for &v in basis.iter() {
        in_b[v] = true;
    }
    let mut r: Vec<f64> = st.b.clone();
    for (i, &v) in basis.iter().enumerate() {
        for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
            r[st.csc_row[k]] -= x_b[i] * st.csc_val[k];
        }
    }
    for (v, &ib) in in_b.iter().enumerate() {
        if !ib && at_ub[v] && st.ub[v].is_finite() {
            for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                r[st.csc_row[k]] -= st.ub[v] * st.csc_val[k];
            }
        }
    }
    let (wr, wi) = r
        .into_iter()
        .enumerate()
        .fold((0.0f64, 0usize), |a, (i, x)| {
            if x.abs() > a.0 { (x.abs(), i) } else { a }
        });
    let xl = extract_x(st, lp, basis, x_b, at_ub);
    // Row activities off the caller CSC (ascending-j accumulation per row —
    // the dense loop's order; stored zeros are ±0.0 no-ops).
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
    let mut worst = (0.0f64, 0usize, 'L');
    for (i, (&sn, &rhs_i)) in lp.senses.iter().zip(&lp.rhs).enumerate() {
        let lhs = act[i];
        let v = match sn {
            'L' => lhs - rhs_i,
            'G' => rhs_i - lhs,
            _ => (lhs - rhs_i).abs(),
        };
        if v > worst.0 {
            worst = (v, i, sn);
        }
    }
    let wb = xl.iter().enumerate().fold((0.0f64, 0usize), |a, (j, &xv)| {
        let v = (lp.lower[j] - xv).max(xv - lp.upper[j]);
        if v > a.0 { (v, j) } else { a }
    });
    let shifted = st.s.iter().filter(|&&sv| sv != 0.0).count();
    let art_bad = basis
        .iter()
        .zip(x_b)
        .filter(|(v, xb)| {
            st.kind[**v] == Kind::Artificial && violation(st.ub[**v], **xb) > FEAS_TOL
        })
        .count();
    eprintln!(
        "rusolve dual FINAL: worst basic bound violation {worst_basic:?}  \
         worst |Bx_b-b| = {wr:.3e} at row {wi}  worst row violation {worst:?}  \
         worst bound {wb:?}  shifted_residue={shifted}  basic_artificials_dirty={art_bad}"
    );
}

// =====================================================================
// The incremental row handle for the MIP cut loop
// =====================================================================

/// Column-block layout of the dual standard-form space, i.e. the index space of
/// a [`TableauRow::coeffs`] and of any other standard-form vector this module
/// exposes. `to_dual_standard` builds the blocks in this order, so the layout
/// is derivable from the caller's LP alone (no internal state needed):
///
/// * `[0, num_structural)` — structural columns, one per ORIGINAL column in
///   the caller's own order (I1: column stability), EXCEPT a fully-free
///   variable (`lower = -inf`, `upper = +inf`) which takes TWO consecutive
///   slots, the x⁺ half then the x⁻ half. So
///   `num_structural = num_col + #fully-free columns`.
/// * `[num_structural, + num_slack)` — the slack/surplus column of each
///   NON-'E' row, in ascending caller-row order (the same subsequence as
///   [`Basis::row_status`] restricted to non-'E' rows).
/// * the remaining `num_artificial` slots — the fixed artificial of each 'E'
///   row, in ascending caller-row order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StdLayout {
    pub num_structural: usize,
    pub num_slack: usize,
    pub num_artificial: usize,
}

impl StdLayout {
    /// Total standard-form column count (= a tableau row's length).
    pub fn total(&self) -> usize {
        self.num_structural + self.num_slack + self.num_artificial
    }
}

/// What sits BASIC at one tableau position, mapped back to caller space so the
/// MIP search can decide whether the row is a cut candidate at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableauBasic {
    /// ORIGINAL column `j` is basic: `TableauRow::basic_value` is its value,
    /// and a fractional integer `j` is exactly the Gomory row to cut on.
    Column(usize),
    /// One half of a fully-free ORIGINAL column `j` is basic (`minus` selects
    /// the x⁻ half). Not a plain caller-variable row — skip for cut generation
    /// (the caller's `x_j = x⁺ - x⁻` is not basic in the tableau sense).
    FreeHalf { col: usize, minus: bool },
    /// The slack/surplus of ORIGINAL row `i` is basic (that row is loose).
    Slack(usize),
    /// The fixed artificial of `'E'` row `i` is basic (degenerate equality
    /// row). No caller variable behind it — skip for cut generation.
    Artificial(usize),
}

/// One row of the simplex tableau at the session basis: the `position`-th
/// row of `B⁻¹A` (`binv_row(basic_position)`), produced by BTRAN on the
/// session's factorization.
#[derive(Debug, Clone)]
pub struct TableauRow {
    /// Basis position this row is the tableau row of (`0 .. num_rows`).
    pub position: usize,
    /// Which caller variable is basic at that position.
    pub basic: TableauBasic,
    /// `x_B[position]` at the session basis: the current value of `basic`.
    /// Gomory's `f0 = basic_value - floor(basic_value)`; `f0` away from 0/1 is
    /// what makes this row a cut candidate.
    pub basic_value: f64,
    /// `ā_j = e_positionᵀ B⁻¹ a_j` over ALL standard-form columns (original
    /// structural columns AND slacks AND 'E'-row artificials), in the
    /// [`StdLayout`] block order. Basic positions carry the identity
    /// (`ā = 1` at this row's own column, `0` at every other basic column);
    /// the Gomory coefficients are the entries at the NONBASIC positions,
    /// which the caller identifies from its [`Basis`] snapshot:
    /// `col_status[j] != Basic` for the structural slot(s) of `j`,
    /// `row_status[i] != Basic` for the slack slot of a non-'E' row `i`.
    pub coeffs: Vec<f64>,
}

/// One caller row to append with [`ReoptSession::add_rows`] — a cut, a
/// no-good row, a clique row: `aᵀx ≤ rhs` ('L'), `aᵀx ≥ rhs` ('G'), `aᵀx = rhs`
/// ('E'), sparse in the CALLER's column space.
#[derive(Debug, Clone)]
pub struct CutRow {
    /// `(column, value)` pairs; columns not listed are 0.0, and zero-valued
    /// entries are dropped on append (the sparse LP input's zero-drop
    /// convention).
    pub entries: Vec<(usize, f64)>,
    /// `'L'` (≤), `'G'` (≥), `'E'` (=).
    pub sense: char,
    /// Finite right-hand side (infinite rows are the caller's ±inf bound
    /// bookkeeping, not a row: the standard form cannot carry them).
    pub rhs: f64,
}

impl CutRow {
    pub fn new(entries: Vec<(usize, f64)>, sense: char, rhs: f64) -> Self {
        Self {
            entries,
            sense,
            rhs,
        }
    }
}

/// The cached warm state of a session basis: [`warm_setup`]'s verified output,
/// kept alive so that the session answers many tableau-row queries from ONE
/// factorization per basis instead of rebuilding one per row.
struct WarmState {
    st: DualStandard,
    basis: Vec<usize>,
    x_b: Vec<f64>,
    at_ub: Vec<bool>,
    pfi: BasisFactorization,
    /// Standard column → caller-space identity ([`TableauBasic`] mapping).
    owner: Vec<TableauBasic>,
    /// Caller column → the basis position whose basic variable is that
    /// column, `usize::MAX` when the column is nonbasic. Built once with the
    /// warm state (an O(m) pass over `basis`/`owner`), so it tracks every
    /// pivot batch exactly as the factorization it lives beside does. A
    /// fully-free column's halves are not the caller's variable and never
    /// appear here, matching [`TableauBasic::Column`].
    col_pos: Vec<usize>,
}

/// Standard column → caller-space identity, the inverse of the block mapping
/// documented on [`StdLayout`].
fn column_owners(lp: &Lp, st: &DualStandard) -> Vec<TableauBasic> {
    let mut owner = vec![TableauBasic::Artificial(usize::MAX); st.ncol()];
    for (j, cols) in st.cols_of.iter().enumerate() {
        match cols.as_slice() {
            [c] => owner[*c] = TableauBasic::Column(j),
            [c1, c2] => {
                owner[*c1] = TableauBasic::FreeHalf {
                    col: j,
                    minus: false,
                };
                owner[*c2] = TableauBasic::FreeHalf {
                    col: j,
                    minus: true,
                };
            }
            _ => {}
        }
    }
    let mut k = 0usize;
    for i in 0..lp.num_row {
        if lp.senses[i] != 'E' {
            owner[st.num_structural + k] = TableauBasic::Slack(i);
            k += 1;
        }
    }
    for (a, i) in (0..lp.num_row).filter(|&i| lp.senses[i] == 'E').enumerate() {
        owner[st.num_structural + k + a] = TableauBasic::Artificial(i);
    }
    owner
}

/// Caller column → basis position, the inverse of the cached warm
/// state's `basis`/`owner` pair (one pass, in position order — a column is
/// basic in at most one position, and free halves map to no caller column).
fn column_positions(num_col: usize, basis: &[usize], owner: &[TableauBasic]) -> Vec<usize> {
    let mut col_pos = vec![usize::MAX; num_col];
    for (pos, &v) in basis.iter().enumerate() {
        if let TableauBasic::Column(j) = owner[v]
            && j < num_col
        {
            col_pos[j] = pos;
        }
    }
    col_pos
}

/// Build (and verify) the warm state of a snapshot on an LP. `None` = the
/// snapshot is unusable for this LP (shape mismatch, singular basis, dual
/// infeasible): the caller keeps the LP and works cold, never a wrong answer.
fn build_warm(lp: &Lp, config: &SimplexConfig, snap: &Basis) -> Option<WarmState> {
    let mut st = to_dual_standard(lp, config).ok()?;
    let (wb, wat) = import_basis_dual(lp, &st, snap)?;
    let (basis, x_b, at_ub, pfi) = warm_setup(&mut st, config, wb, wat)?;
    let owner = column_owners(lp, &st);
    let col_pos = column_positions(lp.num_col, &basis, &owner);
    Some(WarmState {
        st,
        basis,
        x_b,
        at_ub,
        pfi,
        owner,
        col_pos,
    })
}

/// The BTRAN `y = e_positionᵀ B⁻¹` of a cached warm state (shared by
/// [`ReoptSession::tableau_row`] and [`ReoptSession::binv_row`], so both
/// answer from the same solve of the same factorization).
fn bt_row(w: &WarmState, position: usize) -> Option<Vec<f64>> {
    if position >= w.basis.len() {
        return None;
    }
    let mut e = vec![0.0; w.st.m];
    e[position] = 1.0;
    let mut y = vec![0.0; w.st.m];
    w.pfi.solve_bt(&e, &mut y);
    Some(y)
}

/// Which caller bound a standard-side nonbasic position means (Basis contract):
/// the standard lower side is the caller's `lower` when the column has one, and
/// the flip-shift column's UPPER bound otherwise (`x = ub − x̃`, `x̃ = 0`).
fn caller_side(lp: &Lp, j: usize, at_standard_lower: bool) -> BasisStatus {
    if at_standard_lower {
        if lp.lower[j].is_finite() {
            BasisStatus::AtLower
        } else {
            BasisStatus::AtUpper
        }
    } else {
        BasisStatus::AtUpper
    }
}

/// Which side of the [`Basis`] snapshot a repair touches.
enum Slot {
    Col(usize),
    Row(usize),
}

/// Incremental LP handle: the add/remove-rows-keeping-the-basis + tableau-rows
/// surface for the MIP cut loop, on top of the basis snapshot +
/// warm-reoptimize contract. One coherent shape: the session OWNS the caller's
/// [`Lp`] and the current [`Basis`], keeps them in lockstep through row edits,
/// and answers reoptimization and tableau queries from them.
///
/// Cut-loop usage (the shape the MIP search uses):
///
/// ```text
/// let root = solver.solve(&lp)?;                       // cold, presolve OFF
/// let mut s = ReoptSession::new(lp, cfg, root.basis);
/// loop {
///     for pos in 0..s.num_rows() {                     // Gomory rows
///         let r = s.tableau_row(pos);                  // …
///     }
///     s.add_rows(&cuts)?;                              // append
///     let r = s.reoptimize(cutoff, iter_limit)?;       // warm dual pivots
///     if r.status != ReoptStatus::Optimal { break; }
/// }
/// s.remove_rows(&aged)?;                               // cut aging
/// ```
///
/// Why the appended basis stays warm: a new row enters with its slack ('L'/'G')
/// or its 'E'-row artificial BASIC, so the new basis matrix is block
/// triangular with a ±1 diagonal block on the new column and the basis duals
/// are the OLD duals with a 0 appended — every old column keeps its reduced
/// cost, so dual feasibility is preserved EXACTLY and only the new row's basic
/// value is primal infeasible (by the cut violation). `reoptimize` then repairs
/// it in dual pivots, and re-verifies the premise anyway: a snapshot that fails
/// the check falls back to a cold solve (`warm = false`), never a wrong answer.
///
/// Indexing: row indices are APPEND-ONLY and stable across `add_rows` (the
/// dense mirror of the sparse LP input's append-only row invariant);
/// `remove_rows` is the one operation that renumbers, and it returns the
/// old→new map (`rusolve_core::NO_ROW` for deleted rows) so the caller can
/// remap its own cut-index tables. Column indices never move.
pub struct ReoptSession {
    lp: Lp,
    config: SimplexConfig,
    basis: Option<Basis>,
    warm: Option<WarmState>,
}

impl ReoptSession {
    /// Wrap an LP plus the basis currently valid for it — normally
    /// [`crate::simplex::SimplexSolution::basis`] from a presolve-OFF solve,
    /// or the basis a previous [`ReoptSession::reoptimize`] exported.
    /// `None` = no snapshot yet: the first `reoptimize` solves cold and the
    /// session adopts whatever basis that solve exports. A snapshot whose
    /// lengths do not match the LP (another LP's basis, or a shape edited under
    /// [`Self::lp_mut`]) is dropped on entry — the session starts cold rather
    /// than carrying a basis it cannot index.
    pub fn new(lp: Lp, config: &SimplexConfig, basis: Option<Basis>) -> Self {
        let mut s = Self {
            lp,
            config: *config,
            basis,
            warm: None,
        };
        s.drop_inconsistent_basis();
        s
    }

    /// A snapshot must describe the LP it is attached to. A mismatch (a basis
    /// from another LP, or a shape edited under `lp_mut`) is DROPPED rather than
    /// carried: row edits then index safely and the next solve is honestly cold.
    fn drop_inconsistent_basis(&mut self) {
        let ok = self.basis.as_ref().is_some_and(|b| {
            b.col_status.len() == self.lp.num_col && b.row_status.len() == self.lp.num_row
        });
        if !ok {
            self.basis = None;
            self.warm = None;
        }
    }

    /// The LP the session currently holds (row order = caller's, plus appended
    /// rows; column order unchanged).
    pub fn lp(&self) -> &Lp {
        &self.lp
    }

    /// The live basis snapshot, if the session has one.
    pub fn basis(&self) -> Option<&Basis> {
        self.basis.as_ref()
    }

    /// Mutable access for bound and objective edits. Those
    /// keep the snapshot a basis, and [`reoptimize`](Self::reoptimize)
    /// re-verifies the dual-feasibility premise every time, so a premise broken
    /// by an edit here reports `warm = false` with a cold solve rather than a
    /// wrong answer. ROW edits must go through [`Self::add_rows`] /
    /// [`Self::remove_rows`], which maintain the snapshot in lockstep.
    pub fn lp_mut(&mut self) -> &mut Lp {
        self.warm = None;
        &mut self.lp
    }

    /// Current row count (the length of `Basis::row_status`).
    pub fn num_rows(&self) -> usize {
        self.lp.num_row
    }

    /// Standard-form block layout (the tableau-row indexing contract).
    pub fn layout(&self) -> StdLayout {
        let num_structural: usize = (0..self.lp.num_col)
            .map(|j| {
                if self.lp.lower[j].is_finite() || self.lp.upper[j].is_finite() {
                    1
                } else {
                    2 // fully free ⇒ x⁺/x⁻ split
                }
            })
            .sum();
        let num_artificial = self.lp.senses.iter().filter(|&&s| s == 'E').count();
        StdLayout {
            num_structural,
            num_slack: self.lp.num_row - num_artificial,
            num_artificial,
        }
    }

    /// Append rows, keeping the basis warm: each new row's slack ('L'/'G')
    /// or fixed artificial ('E') is recorded BASIC in the snapshot, which is
    /// what preserves dual feasibility (see the type docs). Row indices are
    /// appended at the tail and existing indices do not move.
    ///
    /// A ranged cut `lo ≤ aᵀx ≤ hi` is TWO rows in this dense row form, and the
    /// caller appends them in the I2 order (the `hi` side as 'L' first, then
    /// the `lo` side as 'G') so the standard-form row order matches the
    /// `SparseLp` expansion.
    ///
    /// Validation is all-or-nothing: a rejected call leaves the LP and the
    /// snapshot untouched.
    pub fn add_rows(&mut self, rows: &[CutRow]) -> Result<(), SimplexError> {
        self.drop_inconsistent_basis();
        let n = self.lp.num_col;
        for (k, r) in rows.iter().enumerate() {
            if !matches!(r.sense, 'L' | 'G' | 'E') {
                return Err(SimplexError::Unsupported(format!(
                    "add_rows: row {k} sense '{}'",
                    r.sense
                )));
            }
            if !r.rhs.is_finite() {
                return Err(SimplexError::Unsupported(format!(
                    "add_rows: row {k} rhs {} is not finite",
                    r.rhs
                )));
            }
            let mut seen = vec![false; n];
            for &(j, _) in &r.entries {
                if j >= n {
                    return Err(SimplexError::Unsupported(format!(
                        "add_rows: row {k} references column {j} of {n}"
                    )));
                }
                if seen[j] {
                    return Err(SimplexError::Unsupported(format!(
                        "add_rows: row {k} repeats column {j}"
                    )));
                }
                seen[j] = true;
            }
        }
        // Matrix growth through the lease: absorb the CSC, append one row
        // per CutRow (tail index, existing entries verbatim — the same
        // append-only contract `SparseLp` documents). Entries are sorted by
        // column (the debug-ascending contract of `push_row`; values move
        // bit-exact), zero entries stay structural exactly as the dense
        // push(0.0) version left them.
        let mut bd = rusolve_core::CscBuilder::from_matrix(&self.lp.a);
        for r in rows {
            let mut entries: Vec<(usize, f64)> = r.entries.clone();
            entries.sort_unstable_by_key(|&(j, _)| j);
            bd.push_row(&entries);
            self.lp.rhs.push(r.rhs);
            self.lp.senses.push(r.sense);
            self.lp.num_row += 1;
            if let Some(b) = self.basis.as_mut() {
                // Slack (non-'E') or fixed artificial ('E') enters the basis.
                // With visible 'E'-row artificials the 'E' entry is MEANINGFUL, not a
                // placeholder: `Basic` is exactly the statement "this row's
                // artificial is basic", which is what keeps the appended
                // basis square when the snapshot is re-imported.
                b.row_status.push(BasisStatus::Basic);
            }
        }
        self.lp.a = bd.into_matrix();
        self.warm = None;
        Ok(())
    }

    /// Delete rows (cut aging) and return the old→new row index map
    /// (`rusolve_core::NO_ROW` for deleted rows, the same sentinel convention as
    /// [`rusolve_core::SparseLp::row_map`]).
    ///
    /// What happens to the basis:
    /// * the deleted row's basic row variable — the slack of a non-'E' row,
    ///   or the fixed artificial of an 'E' row whose entry is `Basic` (the
    ///   entry is meaningful for 'E' rows) — leaves WITH the row: the basic set
    ///   loses exactly that column, so it stays square AND feasible — the
    ///   clean case, and `reoptimize` continues warm (test
    ///   `session_remove_loose_row_shrinks_basis_warm`);
    /// * the deleted row was **TIGHT** (slack nonbasic), or an **'E'** row
    ///   whose artificial was nonbasic (entry `AtLower`): one basic variable
    ///   per such row must LEAVE the basis. The repair sheds a basic variable that already sits exactly on a
    ///   bound (`|x_B| ≤ FEAS_TOL` at the lower side, or `|x_B − ũ| ≤ FEAS_TOL`
    ///   at the upper side) and records it nonbasic on that side: deleting a row
    ///   cannot break the feasibility of a point, and the reduced basic solution
    ///   IS that same point with the deleted component dropped, so the snapshot
    ///   stays a square, PRIMAL feasible basis of the reduced LP;
    /// * no such candidate exists (every basic is a fixed artificial, a free
    ///   half, or strictly interior): the snapshot is DROPPED and the next
    ///   `reoptimize` solves cold.
    ///
    /// HONEST LIMIT of the repair (measured, not assumed): deleting a tight row
    /// changes the dual of EVERY column with a coefficient in that row
    /// (`d'_v = d_v + a_{v,r} π_r`), so no basis survives such a deletion
    /// dual-feasibly in general — the shed snapshot is therefore usually
    /// REJECTED by the warm entry's dual-feasibility check and the answer comes
    /// from the cold fallback (correct, `warm = false`: the tight-row tests
    /// observe exactly that). A dual-side repair (HiGHS deletes rows and lets
    /// its phase-1-flavoured dual loop re-establish feasibility) is the
    /// `TODO(perf)` here, to be justified by an aging-heavy MIP run, not by
    /// taste.
    pub fn remove_rows(&mut self, rows: &[usize]) -> Result<Vec<usize>, SimplexError> {
        self.drop_inconsistent_basis();
        let m = self.lp.num_row;
        let mut del: Vec<usize> = rows.to_vec();
        del.sort_unstable();
        del.dedup();
        if let Some(&top) = del.last()
            && top >= m
        {
            return Err(SimplexError::Unsupported(format!(
                "remove_rows: row {top} of {m}"
            )));
        }
        // Removals whose basic row variable (slack, or 'E'-row artificial —
        // the row status is meaningful for 'E' rows) does NOT leave
        // with the row must shed a basic variable each: tight non-'E' rows
        // and 'E' rows whose artificial is nonbasic.
        let mut need_shed = 0usize;
        if self.basis.is_some() {
            for &i in &del {
                let leaves_with_row = self
                    .basis
                    .as_ref()
                    .is_some_and(|b| b.row_status[i] == BasisStatus::Basic);
                if !leaves_with_row {
                    need_shed += 1;
                }
            }
        }
        if need_shed > 0 && !self.shed_basics(need_shed, &del) {
            self.basis = None;
        }
        // Old→new row map up front (the matrix filter needs it for remapping).
        let mut gone = vec![false; m];
        for &i in &del {
            gone[i] = true;
        }
        let mut map = vec![rusolve_core::NO_ROW; m];
        let mut next = 0usize;
        for i in 0..m {
            if !gone[i] {
                map[i] = next;
                next += 1;
            }
        }
        // Matrix filter, one pass over the CSC: surviving entries move
        // bit-exact (stored zeros stay stored — the dense `Vec::remove`
        // kept them too), and the ascending-row layout survives any
        // subsequence, so the remap cannot disorder a column.
        let (optr, oridx, oval) = (self.lp.a.col_ptr(), self.lp.a.row_idx(), self.lp.a.values());
        let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(oval.len());
        for j in 0..self.lp.num_col {
            for k in optr[j]..optr[j + 1] {
                let i = oridx[k] as usize;
                if !gone[i] {
                    triplets.push((map[i], j, oval[k]));
                }
            }
        }
        self.lp.a = rusolve_core::CscMatrix::from_triplets(next, self.lp.num_col, &triplets);
        for &i in del.iter().rev() {
            self.lp.rhs.remove(i);
            self.lp.senses.remove(i);
            self.lp.num_row -= 1;
            if let Some(b) = self.basis.as_mut() {
                b.row_status.remove(i);
            }
        }
        self.warm = None;
        Ok(map)
    }

    /// Move `k` basic variables onto a bound they already sit on (see
    /// [`Self::remove_rows`]). False = no safe repair found.
    fn shed_basics(&mut self, k: usize, del: &[usize]) -> bool {
        if self.basis.is_none() {
            return false;
        }
        self.ensure_warm();
        let mut picks: Vec<(Slot, BasisStatus)> = Vec::with_capacity(k);
        if let Some(w) = self.warm.as_ref() {
            for (pos, &v) in w.basis.iter().enumerate() {
                if picks.len() == k {
                    break;
                }
                let (xb, ub, at_ub) = (w.x_b[pos], w.st.ub[v], w.at_ub[v]);
                let at_lower_side = !at_ub && xb.abs() <= FEAS_TOL;
                let at_upper_side = at_ub && ub.is_finite() && (xb - ub).abs() <= FEAS_TOL;
                if !(at_lower_side || at_upper_side) {
                    continue;
                }
                match w.owner[v] {
                    // A free half or a fixed artificial has no honest
                    // caller-visible status to move to: skip it.
                    TableauBasic::Column(j) => {
                        picks.push((Slot::Col(j), caller_side(&self.lp, j, at_lower_side)));
                    }
                    TableauBasic::Slack(i) if !del.contains(&i) => {
                        // A slack's only bound is its lower side (0 ⇔ row tight).
                        picks.push((Slot::Row(i), BasisStatus::AtLower));
                    }
                    _ => {}
                }
            }
        } else {
            return false;
        }
        if picks.len() < k {
            return false;
        }
        let Some(b) = self.basis.as_mut() else {
            return false;
        };
        for (slot, st) in picks {
            match slot {
                Slot::Col(j) => b.col_status[j] = st,
                Slot::Row(i) => b.row_status[i] = st,
            }
        }
        true
    }

    /// Re-solve the LP the session currently holds, warm from its snapshot:
    /// `cutoff` prunes as soon as the caller-sense dual bound proves the
    /// LP cannot beat it, `iter_limit` defaults to `config.max_iterations`.
    /// An Optimal result adopts the exported basis as the new snapshot (the
    /// grandchildren continue from it); any other status keeps the current
    /// snapshot, which is still a basis of the same LP.
    pub fn reoptimize(
        &mut self,
        cutoff: Option<f64>,
        iter_limit: Option<usize>,
    ) -> Result<ReoptResult, SimplexError> {
        let r = match self.basis.clone() {
            Some(b) => reoptimize(&self.lp, &self.config, &b, cutoff, iter_limit)?,
            None => cold_reopt(&self.lp, &self.config, cutoff, None)?,
        };
        if r.status == ReoptStatus::Optimal {
            self.basis = r.basis.clone();
        }
        // The pivot loop walked away from the basis the cache was factorized
        // for; the next tableau query rebuilds it (one factorization per basis,
        // never one per row).
        self.warm = None;
        Ok(r)
    }

    /// The tableau row of one basic position: BTRAN `e_positionᵀ B⁻¹` on
    /// the session's factorization, dotted with every standard-form column.
    ///
    /// Cost: the FIRST call after a solve / row edit / basis change factorizes
    /// the snapshot once ([`warm_setup`], the same verified entry the dual loop
    /// uses); after that each row is one BTRAN plus one sparse dot per column,
    /// so a full Gomory sweep over the fractional integer rows does not rebuild
    /// anything. `None` = no usable snapshot for this LP (see
    /// [`Self::reoptimize`]'s cold-fallback note) or `position` out of range.
    ///
    /// `coeffs` indexing: [`StdLayout`] block order — structural columns in the
    /// caller's column order (fully-free variables take the x⁺ then x⁻ slot),
    /// then the slacks of the non-'E' rows in ascending row order, then the
    /// 'E'-row artificials in ascending row order. This is the same
    /// caller-space convention as [`Basis`]: the structural slot of caller
    /// column `j` and the slack slot of caller row `i` are exactly the entries
    /// `col_status[j]` / `row_status[i]` describe, so nonbasics are identified
    /// from the snapshot without any internal knowledge.
    pub fn tableau_row(&mut self, position: usize) -> Option<TableauRow> {
        self.ensure_warm();
        let w = self.warm.as_ref()?;
        let y = bt_row(w, position)?;
        let coeffs: Vec<f64> = (0..w.st.ncol()).map(|v| csc_dot(&y, &w.st, v)).collect();
        Some(TableauRow {
            position,
            basic: w.owner[w.basis[position]],
            basic_value: w.x_b[position],
            coeffs,
        })
    }

    /// Does the session hold a usable factorization of its snapshot?
    /// Builds the cached warm state if missing (the same
    /// one-factorization-per-basis contract as [`Self::tableau_row`]).
    /// `false` = no usable snapshot for this LP (see
    /// [`Self::reoptimize`]'s cold-fallback note); every tableau query then
    /// answers `None` until the basis changes.
    pub fn warm_available(&mut self) -> bool {
        self.ensure_warm();
        self.warm.is_some()
    }

    /// The basis position whose basic variable is caller column `j`,
    /// `None` when `j` is nonbasic (or a fully-free column — its halves are
    /// not the caller's variable, and no tableau position holds `x_j`
    /// itself). A cached lookup (see [`WarmState::col_pos`]), so the Gomory
    /// separator locates its candidates without walking the positions.
    pub fn basic_position(&mut self, j: usize) -> Option<usize> {
        self.ensure_warm();
        let pos = *self.warm.as_ref()?.col_pos.get(j)?;
        (pos != usize::MAX).then_some(pos)
    }

    /// `y = e_positionᵀ B⁻¹` over the LP rows (`y.len() = num_rows`,
    /// `y_i` the weight of caller row `i`) — the BTRAN of [`Self::tableau_row`]
    /// WITHOUT the per-column dots. Those dots are only what the full tableau
    /// row needs; a Gomory aggregation needs just `y`, and reads it directly
    /// because a row's slack/artificial column is ±`e_i` (`to_dual_standard`:
    /// `'L'` `+e_i`, `'G'` `−e_i`, `'E'` artificial `+e_i`), so the tableau
    /// row's slack slot holds `±y_i` bit-exactly and the row's weight is
    /// exactly `y_i`. HiGHS answers the same query per candidate with
    /// `Highs::getBasisInverseRowSparse` (`lp_data/Highs.cpp:2272`) and
    /// aggregates with it (`mip/HighsTableauSeparator.cpp:144`).
    pub fn binv_row(&mut self, position: usize) -> Option<Vec<f64>> {
        self.ensure_warm();
        bt_row(self.warm.as_ref()?, position)
    }

    /// Build the cached warm state if it is not there (invalidated by every
    /// mutation the session performs itself).
    fn ensure_warm(&mut self) {
        if self.warm.is_none() {
            self.warm = self
                .basis
                .as_ref()
                .and_then(|b| build_warm(&self.lp, &self.config, b));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::DenseMatrix;
    use crate::rng::XorShift64;

    /// The hyper CHUZR returns exactly the full scan's (ppos, vbest)
    /// under random sparse x_b / DSE-weight moves logged through XbClamp,
    /// and its Σv tracks the dense re-sum to rounding.
    #[test]
    fn p39_hyper_chuzr_matches_full_scan() {
        let m = 3000;
        let mut rng = XorShift64::new(7);
        let mut x_b: Vec<f64> = (0..m).map(|_| rng.uniform(-1.0, 2.0)).collect();
        let ub_b: Vec<f64> = (0..m)
            .map(|i| if i % 3 == 0 { f64::INFINITY } else { 1.0 })
            .collect();
        let mut dse_w: Vec<f64> = (0..m).map(|_| rng.uniform(0.5, 4.0)).collect();
        // Ties in score: duplicate some values so the basis-index
        // tie-break is exercised.
        for i in (0..m).step_by(97) {
            x_b[i] = -0.75;
            dse_w[i] = 1.0;
        }
        let basis: Vec<usize> = (0..m).map(|i| (i * 7919) % m).collect();
        let mut hz = HyperChuzr::new(m, true);
        hz.k = 64; // small set: force many exhaustion rebuilds
        hz.hits_ema = f64::INFINITY; // exercise the set, never the cool-down
        let mut xc = XbClamp::new(m, true);
        let full = |x_b: &[f64], dse_w: &[f64]| {
            let (mut ppos, mut vbest, mut sbest, mut vt) = (usize::MAX, 0.0, 0.0, 0.0);
            for i in 0..m {
                let v = violation(ub_b[i], x_b[i]);
                vt += v;
                if v > FEAS_TOL {
                    let sc = v * v / dse_w[i];
                    if ppos == usize::MAX || sc > sbest || (sc >= sbest && basis[i] < basis[ppos]) {
                        ppos = i;
                        vbest = v;
                        sbest = sc;
                    }
                }
            }
            (ppos, vbest, vt)
        };
        let (mut hits, mut builds) = (0, 0);
        for _ in 0..4000 {
            let (fp, fv, fvt) = full(&x_b, &dse_w);
            match hz.choose(&x_b, &ub_b, &dse_w, true, &basis, &mut xc, 0) {
                Some((hp, hv, _, ht)) => {
                    hits += 1;
                    assert_eq!((hp, hv.to_bits()), (fp, fv.to_bits()));
                    assert!(
                        (ht - fvt).abs() <= 1e-9 * (1.0 + fvt.abs()),
                        "{ht} vs {fvt}"
                    );
                    hz.prev = hp;
                }
                None => {
                    builds += 1;
                    hz.build(&x_b, &ub_b, &dse_w, true);
                    xc.reset_chg();
                    hz.prev = fp;
                }
            }
            if fp == usize::MAX {
                break;
            }
            // Pivot-like move: the chosen row lands feasible (unlogged, as
            // a flip-only exit's clamp), a few rows move (logged, with
            // their DSE weights), sometimes a new large violator appears.
            x_b[fp] = 0.5;
            let mut pat: Vec<u32> = (0..4).map(|_| (rng.next_u64() % m as u64) as u32).collect();
            pat.sort_unstable();
            pat.dedup();
            for &i in &pat {
                let i = i as usize;
                x_b[i] += rng.uniform(-0.6, 0.6);
                dse_w[i] = (dse_w[i] * rng.uniform(0.5, 2.0)).max(1e-4);
            }
            xc.note_chg(Some(&pat));
        }
        assert!(hits > 1000 && builds > 5, "hits={hits} builds={builds}");
    }

    /// The hyper-CHUZR fallback turns off after builds that buy no
    /// hits, sits out a cool-down, probes, and re-enters on a paying probe.
    #[test]
    fn hyper_chuzr_gate_cools_down_and_reenters() {
        let mut hz = HyperChuzr::new(4, true);
        // Builds that die with zero hits drag the mean under the floor.
        let mut offs = 0;
        for _ in 0..16 {
            hz.n_build = 1;
            hz.since = 0;
            if hz.retire() {
                offs += 1;
                break;
            }
        }
        assert_eq!(offs, 1);
        let cool = hz.cool;
        assert_eq!(cool, HCHUZR_COOL_MIN);
        for _ in 0..cool {
            assert!(!hz.gate());
        }
        assert!(hz.probe);
        assert!(hz.gate());
        // The probe build itself always proceeds, at the floor.
        assert!(!hz.retire());
        // A failing probe life cools down again, twice as long.
        hz.since = 0;
        assert!(hz.retire());
        assert_eq!(hz.cool, 2 * HCHUZR_COOL_MIN);
        for _ in 0..hz.cool {
            assert!(!hz.gate());
        }
        assert!(hz.gate() && !hz.retire());
        // A paying probe life keeps it on and resets the cool-down length.
        hz.since = 40;
        assert!(!hz.retire());
        assert_eq!(hz.cool_len, HCHUZR_COOL_MIN);
        assert!(hz.gate());
    }

    #[test]
    fn p36_levers_parse_and_floor_caps_threshold() {
        assert_eq!(parse_zgain_floor(None), ZGAIN_FLOOR_K);
        assert_eq!(parse_zgain_floor(Some("1000")), 1000.0);
        assert_eq!(parse_zgain_floor(Some("0")), 0.0);
        assert_eq!(parse_zgain_floor(Some("-3")), ZGAIN_FLOOR_K);
        assert_eq!(parse_zgain_floor(Some("x")), ZGAIN_FLOOR_K);
        assert_eq!(parse_audit_inc(None), AUDIT_INC_DEFAULT);
        assert!(parse_audit_inc(Some("1")));
        assert!(!parse_audit_inc(Some("0")));
        // k = 0: the plain relative test, untouched.
        let rel = 1e-9 * (1.0 + 6e9);
        assert_eq!(zgain_threshold(-6e9, 6e9, 0.0), rel);
        // pds-scale objective: the floor caps the +6 demand at ~1e-3.
        let t = zgain_threshold(-6e9, 6e9, 1000.0);
        assert!(t < 2e-3 && t > 1e-3, "{t}");
        // heavy cancellation (zabs ≫ |z|): the relative test stays the cap.
        assert_eq!(zgain_threshold(1.0, 1e20, 1000.0), 1e-9 * 2.0);
    }

    /// row-major test data → CSC (dense allowed in tests; zeros stored, the
    /// same convention as the io reader's `Lp`)
    fn cm(rows: &[&[f64]], n: usize) -> rusolve_core::CscMatrix {
        let m = rows.len();
        let triplets: Vec<(usize, usize, f64)> = rows
            .iter()
            .enumerate()
            .flat_map(|(i, row)| {
                // `take(n)`: the dense cm truncated over-long rows silently.
                row.iter().enumerate().take(n).map(move |(j, &v)| (i, j, v))
            })
            .collect();
        rusolve_core::CscMatrix::from_triplets(m, n, &triplets)
    }

    /// Dense form of standard column `v` (dense allowed in tests).
    fn dcol(st: &DualStandard, v: usize) -> Vec<f64> {
        let mut c = vec![0.0; st.m];
        st.gather_col(v, 1.0, &mut c);
        c
    }

    /// Dense columns of a test LP for verification code (stored zeros land
    /// at their positions).
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

    fn cfg() -> SimplexConfig {
        SimplexConfig::default()
    }

    /// Hand LP covering every mapping-table class: lb-shift with finite ub,
    /// lb-shift with inf ub and negative cost (⇒ cost-shifted start),
    /// lb=-inf/finite-ub flip, fully free split, and all three row senses
    /// (incl. a negative-rhs G row kept un-negated on purpose).
    fn kitchen_sink(maximize: bool) -> Lp {
        //      x0   x1   x2   x3   x4
        //  L:   1    2    0   -1    0   <=  4
        //  G:   0    1    1    0   -3   >= -2      (negative rhs)
        //  E:   2    0    0    1   -1   ==  1
        // Bounded in both senses: E forces d3 = d4 on any ray, L then forces
        // d4 >= 2*d1, and G (with d2 <= 0) forces d2 >= 3*d4 - d1 >= 0.
        Lp {
            num_col: 5,
            num_row: 3,
            a: cm(
                &[
                    &[1.0, 2.0, 0.0, -1.0, 0.0],
                    &[0.0, 1.0, 1.0, 0.0, -3.0],
                    &[2.0, 0.0, 0.0, 1.0, -1.0],
                ],
                5,
            ),
            rhs: vec![4.0, -2.0, 1.0],
            objective: vec![1.0, -2.0, 0.0, 3.0, -0.5],
            maximize,
            senses: "LGE".chars().collect(),
            // x0: [0,10]   x1: [2,inf)   x2: (-inf,5] (flip)   x3: free
            // x4: [0,inf)
            lower: vec![0.0, 2.0, f64::NEG_INFINITY, f64::NEG_INFINITY, 0.0],
            upper: vec![10.0, f64::INFINITY, 5.0, f64::INFINITY, f64::INFINITY],
        }
    }

    #[test]
    fn pi_anchor_schedule() {
        // The re-anchor schedule behind RUSOLVE_DUAL_INCP.
        let mut a = PiAnchor::new(20);
        // Armed at start (piw begins zeroed): anchor regardless of d̃ state.
        assert!(a.due(false) && a.due(true));
        a.anchor();
        // Fresh-d̃ path: the cadence gates re-anchors every `cad` PIVOTS…
        for _ in 0..19 {
            a.pivot();
            assert!(!a.due(false), "must not anchor before the cadence is due");
        }
        a.pivot();
        assert!(a.due(false));
        a.anchor();
        // …and the pivot clock restarts at every anchor.
        a.pivot();
        assert!(!a.due(false));
        // Incremental-d̃ path with a LIVE d̃: the cadence never forces an
        // anchor (that path pays no π-BTRAN today; only events do).
        for _ in 0..64 {
            a.pivot();
            assert!(!a.due(true), "live d̃ ⇒ cadence anchors stay off");
        }
        // An event arms the next pass on BOTH paths, and anchoring disarms
        // it and zeroes the cadence clock.
        a.event();
        assert!(a.due(true) && a.due(false));
        a.anchor();
        assert!(!a.due(true) && !a.due(false));
        assert_eq!(a.since, 0);
    }

    #[test]
    fn bound_mapping_matches_design_table() {
        let lp = kitchen_sink(false);
        let st = to_dual_standard(&lp, &cfg()).unwrap();
        // Structural column order mirrors to_standard: x0, x1, x2(flip),
        // x3+/x3- (split), x4 ⇒ 6 columns.
        assert_eq!(st.num_structural, 6);
        // ũ per the bound mapping: ub−lb | inf | inf | inf,inf | inf
        assert_eq!(st.ub[0], 10.0);
        assert!(st.ub[1..6].iter().all(|u| u.is_infinite()));
        // No row negation: b = rhs − Σ A·shift with shifts = (0, 2, 5, 0, 0).
        assert!((st.b[0] - (4.0 - 2.0 * 2.0)).abs() < 1e-12);
        assert!(
            (st.b[1] - (-2.0 - 2.0 - 5.0)).abs() < 1e-12,
            "negative rhs kept: {:?}",
            st.b
        );
        assert!((st.b[2] - 1.0).abs() < 1e-12);
        // Flip column (x2): negated entries, cost −c = 0, shift = ub = 5.
        assert_eq!(dcol(&st, 2), vec![0.0, -1.0, 0.0]);
        assert_eq!(st.cost[2], 0.0);
        assert!(st.flip[2]);
        assert_eq!(st.shift[2], 5.0);
        // Free split (x3): columns a and −a, costs c and −c.
        assert_eq!(dcol(&st, 3), vec![-1.0, 0.0, 1.0]);
        assert_eq!(dcol(&st, 4), vec![1.0, 0.0, -1.0]);
        assert_eq!((st.cost[3], st.cost[4]), (3.0, -3.0));
        // Cost-sign start sides (min problem):
        //  x0 c=+1 ⇒ lower; x1 c=−2 ũ=inf ⇒ lower + shifted;
        //  x2 c̃=0 ⇒ lower; x3+ c=+3 lower, x3− c=−3 ũ=inf ⇒ lower + shifted;
        //  x4 c=−0.5 ũ=inf ⇒ lower + shifted.
        assert!(st.at_ub.iter().take(6).all(|&u| !u));
        assert!(st.s[1] > 0.0 && st.s[4] > 0.0 && st.s[5] > 0.0);
        assert_eq!(st.s[0], 0.0);
        assert_eq!(st.s[2], 0.0);
        assert_eq!(st.s[3], 0.0);
        assert!((st.cost[1] + st.s[1] - eps_s(&cfg())).abs() < 1e-15);
        // Slack/artificial columns: L +e_0, G −e_1, E artificial +e_2 ũ=0.
        assert_eq!(
            &st.kind[6..9],
            &[Kind::Slack, Kind::Slack, Kind::Artificial]
        );
        assert_eq!(dcol(&st, 6), vec![1.0, 0.0, 0.0]);
        assert_eq!(dcol(&st, 7), vec![0.0, -1.0, 0.0]);
        assert_eq!(dcol(&st, 8), vec![0.0, 0.0, 1.0]);
        assert_eq!(st.ub[8], 0.0);
        assert_eq!(st.s.len(), 9);
    }

    #[test]
    fn maximize_negates_costs_and_flips_start_sides() {
        let min = to_dual_standard(&kitchen_sink(false), &cfg()).unwrap();
        let max = to_dual_standard(&kitchen_sink(true), &cfg()).unwrap();
        for v in 0..min.num_structural {
            assert_eq!(max.cost[v], -min.cost[v]);
        }
        for v in min.num_structural..max.cost.len() {
            assert_eq!(max.cost[v], 0.0); // slacks/artificials stay 0
        }
        // x0 cost −1 (< 0) with finite ŵ = 10 ⇒ starts at the upper side.
        assert!(max.at_ub[0]);
        assert_eq!(max.s[0], 0.0); // finite ũ ⇒ flip to ub, no shift needed
    }

    #[test]
    fn initial_basis_and_x_b_match_naive() {
        let lp = kitchen_sink(false);
        let st = to_dual_standard(&lp, &cfg()).unwrap();
        let (basis, diag) = initial_basis(&lp, &st);
        // Row 0 (L) → its slack, row 1 (G) → its slack, row 2 (E) → artificial.
        assert_eq!(basis, vec![6, 7, 8]);
        assert_eq!(diag, vec![1.0, -1.0, 1.0]);
        let x_b = initial_x_b(&st, &basis, &diag);
        // Naive reference: rhs' = b − Σ_{nb at ũ} ũ a; x_b = diag·rhs'.
        let mut rhs = st.b.clone();
        for v in 0..st.ncol() {
            if !basis.contains(&v) && st.at_ub[v] && st.ub[v].is_finite() {
                for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                    rhs[st.csc_row[k]] -= st.ub[v] * st.csc_val[k];
                }
            }
        }
        for i in 0..st.m {
            let want = rhs[i] * diag[i];
            assert!(
                (x_b[i] - if want.abs() < 1e-12 { 0.0 } else { want }).abs() < 1e-12,
                "x_b[{i}] = {} want {want}",
                x_b[i]
            );
        }
        // x_b is NOT clamped: b[1] = −9 with diag −1 gives x_b[1] = +9 here;
        // a negative x_b must stay negative when it occurs (it is DuPh input).
        assert_eq!(x_b[1], 9.0);
    }

    #[test]
    fn extract_x_round_trips_a_guessed_x() {
        let lp = kitchen_sink(false);
        let mut st = to_dual_standard(&lp, &cfg()).unwrap();
        // 3 basis slots. Guess: x0, x1 basic at x̃ = 3, 5 ⇒ x = (3, 7);
        // x4 basic at 4; x2 nonbasic at lower ⇒ x̃ = 0 ⇒ x = shift − 0 = 5;
        // x3 halves nonbasic at lower ⇒ x = 0.
        let basis = vec![0usize, 1, 5];
        let x_b = vec![3.0, 5.0, 4.0];
        for v in 0..st.ncol() {
            st.at_ub[v] = false;
        }
        let x = extract_x(&st, &lp, &basis, &x_b, &st.at_ub);
        for (j, (&g, &w)) in x.iter().zip(&[3.0, 7.0, 5.0, 0.0, 4.0]).enumerate() {
            assert!((g - w).abs() < 1e-12, "x[{j}]: {g} vs {w}");
        }
        // Upper side: x0 nonbasic at ũ (x = 10) and x1 basic at x̃ = 5 with
        // split x3⁺ basic at 2.5, x3⁻ at lower (x = 2.5); x2 basic at x̃ = 2
        // ⇒ x = 5 − 2 = 3 (flip).
        let mut st2 = to_dual_standard(&lp, &cfg()).unwrap();
        for v in 0..st2.ncol() {
            st2.at_ub[v] = false;
        }
        st2.at_ub[0] = true; // x̃0 = ũ0 = 10 ⇒ x0 = 10
        // slots: x1 basic 5 ⇒ x1 = 7; x3+ basic 2.5; x2 basic 2 ⇒ x2 = 3.
        let basis2 = vec![1usize, 3, 2];
        let x_b2 = vec![5.0, 2.5, 2.0];
        let x2 = extract_x(&st2, &lp, &basis2, &x_b2, &st2.at_ub);
        let want2 = [10.0, 7.0, 3.0, 2.5, 0.0];
        for (j, (&g, &w)) in x2.iter().zip(&want2).enumerate() {
            assert!((g - w).abs() < 1e-12, "x[{j}]: {g} vs {w} (x = {x2:?})");
        }
        let _ = &mut st;
    }

    /// Property test: after the cost-sign start, EVERY
    /// nonbasic column is dual feasible w.r.t. the effective costs at π = 0
    /// (slack basis ⇒ d̃ = c̃ + s̃ exactly): lower side ⇒ c + s ≥ −ptol;
    /// upper side ⇒ c ≤ +ptol. Plus: shifted ⇒ ũ = +∞ and nonbasic-lower.
    #[test]
    fn start_is_dual_feasible_property() {
        let mut rng = XorShift64::new(0xD17_2026_1022_0000);
        let c = cfg();
        let ptol = c.pricing_tol;
        for _ in 0..80 {
            let n = 3 + (rng.next_u64() % 6) as usize;
            let m = 2 + (rng.next_u64() % 5) as usize;
            let senses: Vec<char> = (0..m)
                .map(|_| ['L', 'G', 'E'][(rng.next_u64() % 3) as usize])
                .collect();
            let a = cm(
                &(0..m)
                    .map(|_| (0..n).map(|_| rng.uniform(-2.0, 2.0)).collect::<Vec<f64>>())
                    .collect::<Vec<_>>()
                    .iter()
                    .map(|r| r.as_slice())
                    .collect::<Vec<_>>(),
                n,
            );
            let (lower, upper): (Vec<f64>, Vec<f64>) = (0..n)
                .map(|_| match rng.next_u64() % 4 {
                    0 => (0.0, f64::INFINITY),
                    1 => (rng.uniform(0.0, 3.0), rng.uniform(4.0, 10.0)),
                    2 => (f64::NEG_INFINITY, rng.uniform(1.0, 6.0)),
                    _ => (f64::NEG_INFINITY, f64::INFINITY),
                })
                .unzip();
            let lp = Lp {
                num_col: n,
                num_row: m,
                a,
                rhs: (0..m).map(|_| rng.uniform(-8.0, 8.0)).collect(),
                objective: (0..n).map(|_| rng.uniform(-4.0, 4.0)).collect(),
                maximize: rng.next_u64().is_multiple_of(2),
                senses,
                lower,
                upper,
            };
            let st = to_dual_standard(&lp, &c).unwrap();
            let (basis, diag) = initial_basis(&lp, &st);
            for v in 0..st.ncol() {
                if basis.contains(&v) {
                    continue;
                }
                assert!(
                    st.kind[v] == Kind::Structural,
                    "only structurals start nonbasic"
                );
                let d = st.cost[v] + st.s[v]; // π = 0 ⇒ d̃ = effective cost
                if st.at_ub[v] {
                    assert!(d <= ptol, "upper-side dual {d} > ptol");
                    assert!(st.cost[v] < 0.0, "upper start only for negative cost");
                } else {
                    assert!(d >= -ptol, "lower-side dual {d} < -ptol");
                }
                if st.s[v] > 0.0 {
                    assert_eq!(st.ub[v], f64::INFINITY, "only infinite ũ is shifted");
                    assert!(!st.at_ub[v]);
                    assert!((st.cost[v] + st.s[v] - eps_s(&c)).abs() < 1e-12);
                }
            }
            // initial_x_b vs the naive B⁻¹(b − Σ_{at ub} ũa) (independent path).
            let x_b = initial_x_b(&st, &basis, &diag);
            for i in 0..m {
                let mut r = st.b[i];
                for v in 0..st.ncol() {
                    if !basis.contains(&v) && st.at_ub[v] && st.ub[v].is_finite() {
                        r -= st.ub[v] * dcol(&st, v)[i];
                    }
                }
                r *= diag[i];
                if r.abs() < 1e-12 {
                    r = 0.0;
                }
                assert!((x_b[i] - r).abs() <= 1e-9 * r.abs().max(1.0), "x_b[{i}]");
            }
        }
    }

    // ---- Classical dual engine ----

    /// Run the dual engine directly (NOT through the dispatch fallback: a
    /// NumericalDistress `None` must FAIL these tests, not silently pass on
    /// the primal solver).
    fn solve_dual(lp: &Lp) -> SimplexSolution {
        solve_dual_raw(lp, &cfg(), None)
            .unwrap()
            .expect("dual engine returned Distress fallback")
    }

    fn lp6(
        rows: &[&[f64]],
        rhs: Vec<f64>,
        senses: &str,
        obj: Vec<f64>,
        lower: Vec<f64>,
        upper: Vec<f64>,
    ) -> Lp {
        let n = obj.len();
        Lp {
            num_col: n,
            num_row: rhs.len(),
            a: cm(rows, n),
            rhs,
            objective: obj,
            maximize: false,
            senses: senses.chars().collect(),
            lower,
            upper,
        }
    }

    #[test]
    fn dual_min_x_ge_5_is_one_pivot() {
        // min x s.t. x ≥ 5: the textbook dual start enters x at d = +1.
        let s = solve_dual(&lp6(
            &[&[1.0]],
            vec![5.0],
            "G",
            vec![1.0],
            vec![0.0],
            vec![f64::INFINITY],
        ));
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 5.0).abs() < 1e-6, "obj {}", s.objective);
        assert!((s.x[0] - 5.0).abs() < 1e-6, "x {:?}", s.x);
        assert_eq!(s.phase2_iterations, 1, "must be a single dual pivot");
    }

    #[test]
    fn dual_flip_repairs_row_then_primal_step_repairs_debt() {
        // min −x−y s.t. x+y≥4, x−y≥2, 0≤x,y≤6. Start: both at ub ⇒ slack of
        // row 2 is −2 (violated). The ratio test has ONLY the flip of y
        // (upper side, σā>0): the flip alone repairs the row, ending an
        // iteration with no basis change and leaving y on dual DEBT
        // (d̃ = −1 at its new lower side). The optimality audit then takes the
        // primal step y↑ until row 2 goes tight: x=6, y=4, obj −10.
        let s = solve_dual(&lp6(
            &[&[1.0, 1.0], &[1.0, -1.0]],
            vec![4.0, 2.0],
            "GG",
            vec![-1.0, -1.0],
            vec![0.0, 0.0],
            vec![6.0, 6.0],
        ));
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective + 10.0).abs() < 1e-6, "obj {}", s.objective);
        assert!((s.x[0] - 6.0).abs() < 1e-6, "x {:?}", s.x);
        assert!((s.x[1] - 4.0).abs() < 1e-6, "x {:?}", s.x);
    }

    #[test]
    fn dual_landing_saturation_converts_to_flip() {
        // min x+y s.t. x+y≥10, x≤1, y≤100: the entering candidate saturates
        // at ũ_x = 1 (landing L = 10 > 1 ⇒ converts to a flip at its own
        // bound), then y pivots in at 9. Optimum obj 10 at (1, 9).
        let s = solve_dual(&lp6(
            &[&[1.0, 1.0]],
            vec![10.0],
            "G",
            vec![1.0, 1.0],
            vec![0.0, 0.0],
            vec![1.0, 100.0],
        ));
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 10.0).abs() < 1e-6, "obj {}", s.objective);
        assert!((s.x[0] - 1.0).abs() < 1e-6, "x {:?}", s.x);
        assert!((s.x[1] - 9.0).abs() < 1e-6, "x {:?}", s.x);
    }

    #[test]
    fn dual_e_row_only_system() {
        // min x+y s.t. x+y = 7, x ≤ 5: E-row fixed artificial starts basic
        // at +7 (violation above the zero-width bound, σ = −1). x enters,
        // saturates at ũ_x = 5 (flip), y pivots in at 2, the artificial
        // leaves and is deleted. Optimum 7 at (5, 2).
        let s = solve_dual(&lp6(
            &[&[1.0, 1.0]],
            vec![7.0],
            "E",
            vec![1.0, 1.0],
            vec![0.0, 0.0],
            vec![5.0, f64::INFINITY],
        ));
        assert_eq!(s.status, SimplexStatus::Optimal);
        assert!((s.objective - 7.0).abs() < 1e-6, "obj {}", s.objective);
    }

    #[test]
    fn dual_infeasible_verdict() {
        // x ≤ 1, x ≥ 3: after x enters, the L-row slack violation has an
        // empty event set ⇒ verified Infeasible.
        let s = solve_dual(&lp6(
            &[&[1.0], &[1.0]],
            vec![1.0, 3.0],
            "LG",
            vec![1.0],
            vec![0.0],
            vec![f64::INFINITY],
        ));
        assert_eq!(s.status, SimplexStatus::Infeasible);
    }

    #[test]
    fn dual_unbounded_verdict() {
        // min −x, no rows: primal-feasible at the start; cleanup drops the
        // cost shift, the audit sees d = −1 and the primal-step ratio test
        // has no blocker ⇒ verified Unbounded.
        let s = solve_dual(&lp6(
            &[],
            vec![],
            "",
            vec![-1.0],
            vec![0.0],
            vec![f64::INFINITY],
        ));
        assert_eq!(s.status, SimplexStatus::Unbounded);
    }

    /// The primal engine's hand-LP suite re-solved with the dual engine
    /// same optima, every sense/bound class covered.
    #[test]
    fn dual_matches_primal_on_hand_suite() {
        let cases: Vec<(Lp, f64)> = vec![
            (
                lp6(
                    &[&[1.0, 1.0]],
                    vec![5.0],
                    "G",
                    vec![3.0, 2.0],
                    vec![0.0, 0.0],
                    vec![f64::INFINITY; 2],
                ),
                10.0,
            ), // simple_g_min
            (
                lp6(
                    &[&[-1.0]],
                    vec![-3.0],
                    "G",
                    vec![-1.0],
                    vec![0.0],
                    vec![f64::INFINITY],
                ),
                -3.0,
            ), // negated_rhs_row (no negation on the dual side)
            (
                lp6(
                    &[&[1.0, 1.0]],
                    vec![7.0],
                    "E",
                    vec![1.0, 1.0],
                    vec![0.0, 0.0],
                    vec![5.0, f64::INFINITY],
                ),
                7.0,
            ), // equality_row
            (
                lp6(
                    &[&[1.0, 1.0]],
                    vec![12.0],
                    "L",
                    vec![-2.0, -1.0],
                    vec![0.0, 0.0],
                    vec![4.0, 10.0],
                ),
                -16.0,
            ), // ub_saturation_and_row_limit
            (lp6(&[], vec![], "", vec![-1.0], vec![0.0], vec![5.0]), -5.0), // pure_ub_flip_no_rows
            (
                lp6(
                    &[&[1.0]],
                    vec![5.0],
                    "G",
                    vec![1.0],
                    vec![2.0],
                    vec![f64::INFINITY],
                ),
                5.0,
            ), // lower_bound_shift
            (
                lp6(
                    &[&[1.0, 1.0], &[1.0, 2.0], &[2.0, 1.0]],
                    vec![4.0, 6.0, 6.0],
                    "GGG",
                    vec![1.0, 1.0],
                    vec![0.0, 0.0],
                    vec![10.0, 10.0],
                ),
                4.0,
            ), // degenerate_ties_terminate
            (
                lp6(
                    &[&[1.0, -1.0], &[1.0, 0.0]],
                    vec![0.0, 1.0],
                    "EL",
                    vec![1.0, 0.0],
                    vec![0.0, 0.0],
                    vec![f64::INFINITY; 2],
                ),
                0.0,
            ), // zero_artificial_stays_basic
            (
                lp6(
                    &[&[1.0, -1.0]],
                    vec![2.0],
                    "E",
                    vec![1.0, 1.0],
                    vec![0.0, f64::NEG_INFINITY],
                    vec![5.0, f64::INFINITY],
                ),
                -2.0,
            ), // free_variable_split
            (
                lp6(
                    &[&[1.0, 0.0], &[1.0, 1.0]],
                    vec![3.0, 7.0],
                    "LE",
                    vec![0.0, 1.0],
                    vec![0.0, f64::NEG_INFINITY],
                    vec![3.0, f64::INFINITY],
                ),
                4.0,
            ), // x ≤ 3, x + y = 7, y free ⇒ min y = 4 at x = 3
        ];
        for (k, (mut lp, want)) in cases.into_iter().enumerate() {
            let s = solve_dual(&lp);
            assert_eq!(s.status, SimplexStatus::Optimal, "case {k}");
            assert!(
                (s.objective - want).abs() <= 1e-6 * want.abs().max(1.0),
                "case {k}: obj {} want {want}",
                s.objective
            );
            // Maximize mirror: same LP with the sense flipped and sign-flipped
            // costs must hit the sign-flipped optimum (sign-normalization
            // coverage through the whole engine).
            lp.maximize = true;
            lp.objective = lp.objective.iter().map(|c| -c).collect();
            let s2 = solve_dual(&lp);
            assert_eq!(s2.status, SimplexStatus::Optimal, "case {k} max");
            assert!(
                (s2.objective + want).abs() <= 1e-6 * want.abs().max(1.0),
                "case {k} max: obj {} want {}",
                s2.objective,
                -want
            );
        }
    }

    /// Random small LPs, dual engine vs primal engine (no presolve): status
    /// and objective must agree. The real anti-wrong-answer net stays
    /// `tools/diff/simplex_diff.py` (presolve on and off); this is the cheap in-repo cousin.
    #[test]
    fn dual_matches_primal_on_random_lps() {
        use crate::rng::XorShift64;
        let mut rng = XorShift64::new(0xD17_51A1_0000_0001);
        let c = cfg();
        let (mut agree_opt, mut agree_other, mut distress) = (0, 0, 0);
        for t in 0..120 {
            let n = 2 + (t % 5) as usize;
            let m = 1 + (t % 4) as usize;
            let senses: Vec<char> = (0..m)
                .map(|_| ['L', 'G', 'E'][(rng.next_u64() % 3) as usize])
                .collect();
            let a = cm(
                &(0..m)
                    .map(|_| {
                        (0..n)
                            .map(|_| rng.uniform(-2.0, 2.0).round() + rng.uniform(0.0, 0.4))
                            .collect::<Vec<f64>>()
                    })
                    .collect::<Vec<_>>()
                    .iter()
                    .map(|r| r.as_slice())
                    .collect::<Vec<_>>(),
                n,
            );
            let (lower, upper): (Vec<f64>, Vec<f64>) = (0..n)
                .map(|_| match rng.next_u64() % 4 {
                    0 => (0.0, f64::INFINITY),
                    1 => (rng.uniform(0.0, 2.0), rng.uniform(3.0, 8.0)),
                    2 => (f64::NEG_INFINITY, rng.uniform(1.0, 5.0)),
                    _ => (f64::NEG_INFINITY, f64::INFINITY),
                })
                .unzip();
            let lp = Lp {
                num_col: n,
                num_row: m,
                a,
                rhs: (0..m).map(|_| rng.uniform(-6.0, 6.0)).collect(),
                objective: (0..n).map(|_| rng.uniform(-3.0, 3.0)).collect(),
                maximize: rng.next_u64().is_multiple_of(2),
                senses,
                lower,
                upper,
            };
            let dual = solve_dual_raw(&lp, &c, None).unwrap();
            let primal = crate::simplex::SimplexSolver::new(SimplexConfig {
                presolve: false,
                ..SimplexConfig::default()
            })
            .solve(&lp)
            .unwrap();
            let Some(dual) = dual else {
                // Distress fallback is legal; primal is ground truth there.
                distress += 1;
                continue;
            };
            assert_eq!(
                dual.status, primal.status,
                "t{t}: dual {:?} vs primal {:?}",
                dual.status, primal.status
            );
            if dual.status == SimplexStatus::Optimal {
                assert!(
                    (dual.objective - primal.objective).abs()
                        <= 1e-6 * primal.objective.abs().max(1.0),
                    "t{t}: dual obj {} vs primal {}",
                    dual.objective,
                    primal.objective
                );
                // x feasibility of the dual answer at the caller level.
                for (j, &xv) in dual.x.iter().enumerate() {
                    assert!(
                        xv >= lp.lower[j] - 1e-6 && xv <= lp.upper[j] + 1e-6,
                        "t{t} x[{j}]"
                    );
                }
                let cols = dense_cols(&lp);
                for i in 0..m {
                    let lhs: f64 = (0..n).map(|j| cols[j][i] * dual.x[j]).sum();
                    let v = match lp.senses[i] {
                        'L' => lhs - lp.rhs[i],
                        'G' => lp.rhs[i] - lhs,
                        _ => (lhs - lp.rhs[i]).abs(),
                    };
                    assert!(v <= 1e-5, "t{t} row {i} violated by {v}");
                }
                agree_opt += 1;
            } else {
                agree_other += 1;
            }
        }
        // Statuses must agree everywhere the dual engine dared answer, and
        // the mix must cover the space (optima AND infeasible/unbounded).
        assert_eq!(agree_opt + agree_other + distress, 120);
        assert!(agree_opt >= 30, "too few agreed optima: {agree_opt}");
        assert!(
            agree_other >= 10,
            "too few agreed non-optimal verdicts: {agree_other}"
        );
        assert!(distress <= 12, "distress fallback too often: {distress}");
    }

    /// Regression (the t27 autopsy): a max-violation leaving-row pick let a
    /// flip for one row be UNDONE by a saturation chosen for another row
    /// (one column ping-ponging between its bounds, θ=0 both ways, exact
    /// state repetition ⇒ IterationLimit where the answer is Infeasible).
    /// The fix is the pinned leaving row + the δ-commit; this LP is the
    /// 28th random draw of the parity generator (seed 0xD17_51A1_0000_0001)
    /// and must report Infeasible quickly, never IterationLimit.
    #[test]
    fn dual_pinned_row_detects_infeasible_t27_pong_lp() {
        let col = [
            vec![
                -0.7769525734319898,
                1.0818487639434982,
                -1.8843138323560076,
                -1.8230209098226187,
            ],
            vec![
                0.057016557404297875,
                0.03976708135430665,
                1.2984756157878687,
                0.34976961571945453,
            ],
            vec![
                1.1734168947366015,
                -0.7488577424758041,
                2.174035920263932,
                -1.6637712539861995,
            ],
            vec![
                2.2145773090290177,
                -1.9867559886919068,
                1.2092357600546262,
                0.029130925921242758,
            ],
        ];
        // dense column-major draw (kept verbatim from the parity generator)
        // → CSC, zeros stored
        let a = {
            let triplets: Vec<(usize, usize, f64)> = col
                .iter()
                .enumerate()
                .flat_map(|(j, c)| c.iter().enumerate().map(move |(i, &v)| (i, j, v)))
                .collect();
            rusolve_core::CscMatrix::from_triplets(4, 4, &triplets)
        };
        let lp = Lp {
            num_col: 4,
            num_row: 4,
            a,
            rhs: vec![
                5.732734806176438,
                -2.230645294518506,
                5.394278121603875,
                -1.5304691563315407,
            ],
            objective: vec![
                1.2432058132802997,
                -0.9810602934665367,
                1.3131088278129663,
                1.3075121322648986,
            ],
            maximize: true,
            senses: vec!['G', 'E', 'E', 'G'],
            lower: vec![
                0.23173555782463362,
                0.21364316870212186,
                f64::NEG_INFINITY,
                0.0,
            ],
            upper: vec![
                3.8251766192171934,
                6.6222182214010585,
                1.3820360884691305,
                f64::INFINITY,
            ],
        };
        let s = solve_dual(&lp);
        assert_eq!(s.status, SimplexStatus::Infeasible, "pong regression LP");
        assert!(s.phase2_iterations < 1000, "must converge, not crawl");
        // The cold dual engine carries the verified ray too — this LP is
        // the verifier's torture case ('G'/'E' rows, a free column, an
        // infinite upper bound, maximize).
        let y = s
            .dual_ray
            .expect("dual-engine Infeasible carries a Farkas ray");
        assert!(y.iter().any(|&v| v != 0.0), "ray must be nonzero");
        assert!(farkas_ok(&lp, &y), "ray must verify: {y:?}");
    }

    // ---- Warm dual reoptimization with an objective cutoff ----

    fn cfg_np() -> SimplexConfig {
        SimplexConfig {
            presolve: false,
            ..cfg()
        }
    }

    /// MIP-relaxation-flavoured root LP: max x + 2y s.t. 3x + 4y ≤ 12,
    /// 2x + 2y ≤ 7, 0 ≤ x, y ≤ 10. Root optimum (2, 1.5) — fractional y,
    /// both rows tight (slacks nonbasic at 0, x and y basic), so branching
    /// on y hits exactly the MIP node case: parent basis stays dual
    /// feasible, goes primal infeasible.
    fn node_root() -> Lp {
        Lp {
            num_col: 2,
            num_row: 2,
            a: cm(&[&[3.0, 4.0], &[2.0, 2.0]], 2),
            rhs: vec![12.0, 7.0],
            objective: vec![-1.0, -2.0], // min −(x + 2y)
            maximize: false,
            senses: "LL".chars().collect(),
            lower: vec![0.0, 0.0],
            upper: vec![10.0, 10.0],
        }
    }

    fn cold_solve(lp: &Lp) -> SimplexSolution {
        SimplexSolver::new(cfg_np()).solve(lp).unwrap()
    }

    /// Load a captured warm-reopt repro (`.mps` + `.basis`, README in
    /// `benches/data/repros/`) as the caller LP and offered basis snapshot.
    fn load_repro(name: &str) -> (Lp, Basis) {
        use rusolve_io::{ObjSense, RowSense};
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../../benches/data/repros/");
        let text =
            std::fs::read_to_string(dir.to_string() + name + ".mps").expect("repro LP present");
        let m = rusolve_io::read_mps(&text).expect("repro parses");
        assert!(!m.bounds.iter().any(|b| b.integer));
        assert_eq!(m.ranges.iter().filter(|&&r| r != 0.0).count(), 0);
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
            // `.abs()`: the capture's LO lines carry `-0.0`; IEEE-equal to
            // 0.0 anyway, abs only keeps the test's LP byte-clean.
            lower: m.bounds.iter().map(|b| b.lower.abs()).collect(),
            upper: m.bounds.iter().map(|b| b.upper).collect(),
        };
        let bs = std::fs::read_to_string(dir.to_string() + name + ".basis")
            .expect("repro basis present");
        let mut basis = Basis {
            col_status: vec![BasisStatus::AtLower; lp.num_col],
            row_status: vec![BasisStatus::AtLower; lp.num_row],
        };
        for line in bs.lines().skip(1) {
            let mut it = line.split_whitespace();
            let st = |s: &str| match s {
                "Basic" => BasisStatus::Basic,
                "AtLower" => BasisStatus::AtLower,
                "AtUpper" => BasisStatus::AtUpper,
                "Fixed" => BasisStatus::Fixed,
                other => panic!("bad basis status {other:?}"),
            };
            match it.next() {
                Some("C") => {
                    let j: usize = it.next().unwrap().parse().unwrap();
                    basis.col_status[j] = st(it.next().unwrap());
                }
                Some("R") => {
                    let i: usize = it.next().unwrap().parse().unwrap();
                    basis.row_status[i] = st(it.next().unwrap());
                }
                _ => {}
            }
        }
        (lp, basis)
    }

    // ---- E-row basis fidelity (export/import honors the artificial) ----

    /// E-row-heavy LP with a REDUNDANT equality (E2 = E0 + E1): the row rank
    /// is 4 < 5, so one E-row artificial must stay basic at 0 in every basis,
    /// and the exported snapshot must SHOW which one while the import
    /// puts exactly that one back. G4 is a slack row that is never tight
    /// (every feasible point has x0 ≥ 2, x1 ≥ 0 ⇒ x0 + x1 ≥ 2 > 0.5), so its
    /// artificial can never be basic at an optimum and the dual import is
    /// exact. min −(x0+x1+x2+x3):
    ///
    /// ```text
    /// E0:  x0      + x2      = 4
    /// E1:       x1      + x3 = 3
    /// E2:  x0 + x1 + x2 + x3 = 7      (= E0 + E1, redundant)
    /// E3:  x0           + x3 = 5      (independent)
    ///  x0 + x1           ≥ 0.5    (never tight)
    /// ```
    fn e_row_heavy_lp() -> Lp {
        Lp {
            num_col: 4,
            num_row: 5,
            a: cm(
                &[
                    &[1.0, 0.0, 1.0, 0.0], // E0
                    &[0.0, 1.0, 0.0, 1.0], // E1
                    &[1.0, 1.0, 1.0, 1.0], // E2 = E0 + E1 (redundant)
                    &[1.0, 0.0, 0.0, 1.0], // E3 independent
                    &[1.0, 1.0, 0.0, 0.0], // G4
                ],
                4,
            ),
            rhs: vec![4.0, 3.0, 7.0, 5.0, 0.5],
            objective: vec![-1.0, -1.0, -1.0, -1.0],
            maximize: false,
            senses: "EEEEG".chars().collect(),
            lower: vec![0.0, 0.0, 0.0, 0.0],
            upper: vec![4.0, 3.0, 4.0, 3.0],
        }
    }

    #[test]
    fn e_row_status_export_import_export_is_idempotent() {
        // The contract test: on an E-row-heavy LP, exporting a
        // basis, importing it into the dual standard form, and re-exporting
        // must return the SAME snapshot — including WHICH E-row artificial
        // is basic. Under the placeholder contract the E entries all said
        // Basic and the re-export after a positional fill said whatever the
        // fill picked, not what the basis held.
        let lp = e_row_heavy_lp();
        let root = cold_solve(&lp);
        assert_eq!(root.status, SimplexStatus::Optimal, "test LP must solve");
        let snap = root.basis.clone().expect("presolve-off Optimal exports");
        let st = to_dual_standard(&lp, &cfg_np()).unwrap();
        let (nb, nat) = import_basis_dual(&lp, &st, &snap).expect("the LP's own basis imports");
        // Square basis.
        assert_eq!(nb.len(), lp.num_row);
        assert!(
            nb.iter()
                .all(|&v| st.kind[v] != Kind::Artificial || st.ub[v] == 0.0),
            "only fixed E artificials may stand in for artificials"
        );
        let re = export_basis_dual(&lp, &st, &nb, &nat);
        assert_eq!(
            re, snap,
            "export → import → export must be the identity on an E-row-heavy LP"
        );
    }

    #[test]
    fn e_row_heavy_warm_reopt_accepts_the_exported_basis() {
        // The MIP-facing acceptance: a snapshot the cold primal exported for
        // an LP reoptimizes the UNCHANGED LP warm (the block_angular /
        // network_design rejection class). The redundant E row forces a
        // basic artificial, so the old positional fill could not be sure to re-pick it.
        let lp = e_row_heavy_lp();
        let root = cold_solve(&lp);
        let snap = root.basis.clone().expect("basis");
        let r = reoptimize(&lp, &cfg_np(), &snap, None, None).unwrap();
        assert!(r.warm, "the LP's own basis must not be rejected");
        assert_eq!(r.status, ReoptStatus::Optimal);
        assert!(
            (r.objective - root.objective).abs() < 1e-9,
            "warm {} vs cold {}",
            r.objective,
            root.objective
        );
        assert!(
            r.iters < root.phase1_iterations + root.phase2_iterations + 1,
            "reopt of the unchanged LP takes no pivots: iters {}",
            r.iters
        );
    }

    #[test]
    fn e_row_heavy_warm_primal_start_skips_phase_one() {
        // The primal warm entry honors the E entries too: re-solving the
        // unchanged LP from its own export reconstructs the basis exactly
        // (phase-1 iterations 0, Optimal).
        let lp = e_row_heavy_lp();
        let root = cold_solve(&lp);
        let snap = root.basis.clone().expect("basis");
        let warm = SimplexSolver::new(cfg_np())
            .solve_with(&lp, Some(&snap))
            .unwrap();
        assert_eq!(warm.status, SimplexStatus::Optimal);
        assert_eq!(warm.phase1_iterations, 0, "warm start skipped phase 1");
        assert!(
            (warm.objective - root.objective).abs() < 1e-9,
            "warm {} vs cold {}",
            warm.objective,
            root.objective
        );
    }

    // ---- Duals on ReoptResult (same convention as SimplexSolution) ----

    /// The sign-convention audit for a ReoptResult dual pair: lengths, the
    /// reduced-cost identity `d_j = c_j − πᵀa_j` for EVERY column, and the
    /// sign/complementary-slackness pattern against the reported x —
    /// minimize: d ≥ 0 at lower, d ≤ 0 at upper, d = 0 strictly interior,
    /// π ≤ 0 on 'L', π ≥ 0 on 'G', π_i ≠ 0 ⇒ row i tight; maximize flips
    /// every sign (the convention of simplex.rs
    /// `SimplexSolution::reduced_costs`).
    fn check_reopt_duals(lp: &Lp, r: &ReoptResult) {
        assert!(r.duals_valid, "Optimal must carry duals_valid");
        assert_eq!(r.row_duals.len(), lp.num_row);
        assert_eq!(r.reduced_costs.len(), lp.num_col);
        let s = if lp.maximize { -1.0 } else { 1.0 };
        let tol = 1e-6;
        let cols = dense_cols(lp);
        for j in 0..lp.num_col {
            let dot: f64 = (0..lp.num_row).map(|i| r.row_duals[i] * cols[j][i]).sum();
            let want = lp.objective[j] - dot;
            assert!(
                (r.reduced_costs[j] - want).abs() <= 1e-8 * (1.0 + want.abs()),
                "col {j}: d {} want c−πᵀa {want}",
                r.reduced_costs[j]
            );
            let (lb, ub, xj, d) = (lp.lower[j], lp.upper[j], r.x[j], r.reduced_costs[j]);
            let at_lb = (xj - lb).abs() <= 1e-7;
            let at_ub = ub.is_finite() && (xj - ub).abs() <= 1e-7;
            if at_lb && !at_ub {
                assert!(s * d >= -tol, "col {j}: x at lb, s·d {d} < 0");
            } else if at_ub && !at_lb {
                assert!(s * d <= tol, "col {j}: x at ub, s·d {d} > 0");
            } else if lb < xj - 1e-7 && xj + 1e-7 < ub {
                assert!(d.abs() <= tol, "col {j}: x interior, d {d} ≠ 0");
            }
        }
        for i in 0..lp.num_row {
            let p = r.row_duals[i];
            match lp.senses[i] {
                'L' => {
                    assert!(s * p <= tol, "row {i} 'L': s·π {p} > 0");
                    let lhs: f64 = (0..lp.num_col).map(|j| cols[j][i] * r.x[j]).sum();
                    if s * p < -tol {
                        assert!(
                            (lhs - lp.rhs[i]).abs() <= 1e-5,
                            "row {i}: π ≠ 0 but row not tight (lhs {lhs}, rhs {})",
                            lp.rhs[i]
                        );
                    }
                }
                'G' => {
                    assert!(s * p >= -tol, "row {i} 'G': s·π {p} < 0");
                    let lhs: f64 = (0..lp.num_col).map(|j| cols[j][i] * r.x[j]).sum();
                    if s * p > tol {
                        assert!(
                            (lhs - lp.rhs[i]).abs() <= 1e-5,
                            "row {i}: π ≠ 0 but row not tight (lhs {lhs}, rhs {})",
                            lp.rhs[i]
                        );
                    }
                }
                _ => {} // 'E': free sign, no slack
            }
        }
    }

    #[test]
    fn reopt_optimal_duals_match_cold_b4_solve() {
        // Acceptance: the warm reopt of the unchanged LP carries the SAME
        // caller-space dual pair as the cold primal solve of that LP
        // (both presolve off), and it passes the convention audit.
        let lp = node_root();
        let cold = cold_solve(&lp);
        assert!(
            cold.duals_valid,
            "primal presolve-off Optimal carries duals"
        );
        let snap = cold.basis.clone().expect("basis");
        let r = reoptimize(&lp, &cfg_np(), &snap, None, None).unwrap();
        assert!(r.warm, "test LP must answer warm");
        assert_eq!(r.status, ReoptStatus::Optimal);
        check_reopt_duals(&lp, &r);
        for i in 0..lp.num_row {
            assert!(
                (r.row_duals[i] - cold.row_duals[i]).abs()
                    <= 1e-8 * (1.0 + cold.row_duals[i].abs()),
                "row_duals[{i}]: warm {}, cold {}",
                r.row_duals[i],
                cold.row_duals[i]
            );
        }
        for j in 0..lp.num_col {
            assert!(
                (r.reduced_costs[j] - cold.reduced_costs[j]).abs()
                    <= 1e-8 * (1.0 + cold.reduced_costs[j].abs()),
                "reduced_costs[{j}]: warm {}, cold {}",
                r.reduced_costs[j],
                cold.reduced_costs[j]
            );
        }
    }

    #[test]
    fn reopt_child_optimal_duals_valid_and_match_truth() {
        // The MIP node case: a bound-change child reoptimized warm comes back
        // with valid duals matching a cold solve of the child LP — what
        // reduced-cost fixing consumes.
        let mut lp = node_root();
        let root = cold_solve(&lp);
        let b = root.basis.clone().expect("root export");
        lp.upper[1] = 1.0;
        let r = reoptimize(&lp, &cfg_np(), &b, None, None).unwrap();
        assert_eq!(r.status, ReoptStatus::Optimal);
        assert!(r.warm, "premise holds on this crafted child");
        check_reopt_duals(&lp, &r);
        let truth = cold_solve(&lp);
        assert!((r.objective - truth.objective).abs() < 1e-9);
        for i in 0..lp.num_row {
            assert!(
                (r.row_duals[i] - truth.row_duals[i]).abs()
                    <= 1e-8 * (1.0 + truth.row_duals[i].abs()),
                "row_duals[{i}]: warm {}, cold truth {}",
                r.row_duals[i],
                truth.row_duals[i]
            );
        }
    }

    #[test]
    fn dual_lever_cold_optimal_reports_duals_on_every_bound_class() {
        // The shared exit also feeds the cold dual engine (RUSOLVE_DUAL=1).
        // kitchen_sink covers every bound-mapping class (lb-shift, flip,
        // fully-free split) and all three row senses, both objective senses:
        // the dual engine's π must agree with the primal engine's duals.
        for maximize in [false, true] {
            let lp = kitchen_sink(maximize);
            let cfg = SimplexConfig {
                presolve: false,
                dual: true,
                ..cfg()
            };
            let s = SimplexSolver::new(cfg).solve(&lp).unwrap();
            assert_eq!(s.status, SimplexStatus::Optimal, "maximize {maximize}");
            assert!(s.duals_valid, "dual engine must report duals");
            let r = ReoptResult {
                status: ReoptStatus::Optimal,
                obj_bound: s.objective,
                iters: 0,
                warm: true,
                x: s.x.clone(),
                objective: s.objective,
                basis: None,
                row_duals: s.row_duals.clone(),
                reduced_costs: s.reduced_costs.clone(),
                duals_valid: s.duals_valid,
                dual_ray: None,
                prune_duals: Vec::new(),
                dse_weights: None,
            };
            check_reopt_duals(&lp, &r);
            let primal = cold_solve(&lp);
            assert!(primal.duals_valid);
            assert!(
                (s.objective - primal.objective).abs() <= 1e-9 * (1.0 + primal.objective.abs())
            );
            for i in 0..lp.num_row {
                assert!(
                    (s.row_duals[i] - primal.row_duals[i]).abs()
                        <= 1e-7 * (1.0 + primal.row_duals[i].abs()),
                    "max {maximize} row_duals[{i}]: dual engine {}, primal {}",
                    s.row_duals[i],
                    primal.row_duals[i]
                );
            }
        }
    }

    /// `warm_flip_boxed`: the root optimum has x nonbasic at 0 (with y
    /// basic). Moving x to its upper bound 10 leaves π unchanged but makes
    /// d_x > 0 sit on the wrong side: the only dual infeasibility, on a
    /// boxed column. Off, the warm entry rejects (cold fallback); on, the
    /// flip restores the start and the warm answer is the root optimum.
    #[test]
    fn warm_entry_flips_boxed_dual_infeasible_column() {
        let lp = node_root();
        let root = cold_solve(&lp);
        let mut b = root.basis.clone().expect("root export");
        assert_eq!(b.col_status[0], BasisStatus::AtLower);
        b.col_status[0] = BasisStatus::AtUpper;
        let off = reoptimize(&lp, &cfg_np(), &b, None, None).unwrap();
        assert!(
            !off.warm,
            "dual-infeasible start must reject without the flip"
        );
        let cfg = SimplexConfig {
            warm_flip_boxed: true,
            ..cfg_np()
        };
        let on = reoptimize(&lp, &cfg, &b, None, None).unwrap();
        assert!(on.warm, "the boxed flip repairs the start");
        assert_eq!(on.status, ReoptStatus::Optimal);
        assert!(
            (on.objective - root.objective).abs() < 1e-9,
            "obj {}",
            on.objective
        );
        for (a, x) in on.x.iter().zip(&root.x) {
            assert!((a - x).abs() < 1e-9, "x {a} vs {x}");
        }
    }

    #[test]
    fn warm_dual_reopt_upper_bound_child() {
        // Child y ≤ 1: parent basis dual feasible (premise holds), primal
        // infeasible (y basic at 1.5 > 1). Truth: y=1, x=2.5, obj −4.5.
        let mut lp = node_root();
        let root = cold_solve(&lp);
        let b = root.basis.clone().expect("root export");
        lp.upper[1] = 1.0;
        let r = reoptimize(&lp, &cfg_np(), &b, None, None).unwrap();
        assert!(r.warm, "bound change must not break dual feasibility");
        assert_eq!(r.status, ReoptStatus::Optimal);
        assert!((r.objective + 4.5).abs() < 1e-9, "obj {}", r.objective);
        let truth = cold_solve(&lp);
        assert!((r.objective - truth.objective).abs() < 1e-9);
        for (a, x) in r.x.iter().zip(&truth.x) {
            assert!((a - x).abs() < 1e-7, "x drift {a} vs {x}");
        }
        // The grandchildren need this: Optimal reopt exports a basis.
        assert!(r.basis.is_some());
    }

    #[test]
    fn warm_dual_reopt_lower_bound_child_beats_cold() {
        // Child y ≥ 2 (the other branch): shifted form goes primal
        // infeasible (x̃_x basic at −2). Truth: x=4/3, y=2, obj −16/3.
        let mut lp = node_root();
        let root = cold_solve(&lp);
        let b = root.basis.clone().expect("root export");
        lp.lower[1] = 2.0;
        let r = reoptimize(&lp, &cfg_np(), &b, None, None).unwrap();
        assert!(r.warm);
        assert_eq!(r.status, ReoptStatus::Optimal);
        let truth = cold_solve(&lp);
        assert!(
            (r.objective - truth.objective).abs() < 1e-9,
            "warm {} vs cold truth {}",
            r.objective,
            truth.objective
        );
        let cold = truth.phase1_iterations + truth.phase2_iterations;
        assert!(
            r.iters <= cold,
            "warm {r:?} iters {} should not exceed cold {cold}",
            r.iters
        );
    }

    #[test]
    fn cutoff_prunes_when_bound_reaches_it() {
        // Child y ≤ 1 (optimum −4.5). A cutoff of −4.6 (incumbent −4.6 in
        // min sense: prune once the bound says the child cannot go below)
        // must report Pruned with obj_bound ≥ cutoff; a cutoff of −4.0 sits
        // BELOW the dual path (−5 → −4.5, never reaching −4.0) ⇒ Optimal.
        let mut lp = node_root();
        let root = cold_solve(&lp);
        let b = root.basis.clone().expect("root export");
        lp.upper[1] = 1.0;
        let r = reoptimize(&lp, &cfg_np(), &b, Some(-4.6), None).unwrap();
        assert_eq!(r.status, ReoptStatus::Pruned, "bound ≥ −4.6 ⇒ prune");
        assert!(r.obj_bound >= -4.6 - 1e-9, "obj_bound {}", r.obj_bound);
        let r2 = reoptimize(&lp, &cfg_np(), &b, Some(-4.0), None).unwrap();
        assert_eq!(r2.status, ReoptStatus::Optimal);
        assert!((r2.objective + 4.5).abs() < 1e-9);
    }

    /// The multipliers returned at an early prune give a GLOBALLY valid
    /// proof row `(c − Aᵀy)ᵀx ≤ U − yᵀb` (y clipped to the row-sense signs):
    /// it holds on a grid over the ROOT box for every row-feasible point
    /// with `cᵀx ≤ U`, and the child's box violates it (that is the prune).
    #[test]
    fn prune_duals_give_a_globally_valid_proof_row() {
        let lp0 = node_root();
        let root = cold_solve(&lp0);
        let lp0_cols = dense_cols(&lp0);
        let b = root.basis.clone().expect("root export");
        let mut lp = node_root();
        lp.upper[1] = 1.0;
        let u = -4.6;
        let r = reoptimize(&lp, &cfg_np(), &b, Some(u), None).unwrap();
        assert_eq!(r.status, ReoptStatus::Pruned);
        assert_eq!(r.prune_duals.len(), lp.num_row);
        let y: Vec<f64> = r
            .prune_duals
            .iter()
            .zip(&lp.senses)
            .map(|(&v, &s)| match s {
                'L' => v.min(0.0),
                'G' => v.max(0.0),
                _ => v,
            })
            .collect();
        let cols = dense_cols(&lp);
        let coef: Vec<f64> = (0..2)
            .map(|j| lp.objective[j] - (0..2).map(|i| cols[j][i] * y[i]).sum::<f64>())
            .collect();
        let rhs = u - y.iter().zip(&lp.rhs).map(|(a, b)| a * b).sum::<f64>();
        // Child box minimum activity exceeds rhs.
        let child_min: f64 = (0..2)
            .map(|j| {
                if coef[j] > 0.0 {
                    coef[j] * lp.lower[j]
                } else {
                    coef[j] * lp.upper[j]
                }
            })
            .sum();
        assert!(
            child_min > rhs - 1e-9,
            "child box must violate: {child_min} vs {rhs}"
        );
        // Global validity on the root box.
        let n = 200;
        for a in 0..=n {
            for c in 0..=n {
                let x = [10.0 * a as f64 / n as f64, 10.0 * c as f64 / n as f64];
                let rows_ok = (0..2).all(|i| {
                    let act: f64 = (0..2).map(|j| lp0_cols[j][i] * x[j]).sum();
                    act <= lp0.rhs[i] + 1e-12
                });
                let obj: f64 = (0..2).map(|j| lp0.objective[j] * x[j]).sum();
                if rows_ok && obj <= u {
                    let act: f64 = (0..2).map(|j| coef[j] * x[j]).sum();
                    assert!(act <= rhs + 1e-9, "proof cuts off feasible x={x:?}");
                }
            }
        }
    }

    /// Exact `‖e_kᵀB⁻¹‖²` for the basic variables of `b`, in the
    /// [`DseWeights`] order (basic columns ascending, then basic rows).
    fn exact_dse(lp: &Lp, b: &Basis) -> Vec<f64> {
        let m = lp.num_row;
        let dense = dense_cols(lp);
        let mut cols: Vec<Vec<f64>> = Vec::new();
        for j in 0..lp.num_col {
            if b.col_status[j] == BasisStatus::Basic {
                cols.push(dense[j].clone());
            }
        }
        for i in 0..m {
            if b.row_status[i] == BasisStatus::Basic {
                let mut e = vec![0.0; m];
                e[i] = 1.0;
                cols.push(e);
            }
        }
        assert_eq!(cols.len(), m);
        // Gauss-Jordan on [B | I]: B[r][c] = cols[c][r].
        let mut a: Vec<Vec<f64>> = (0..m)
            .map(|r| {
                let mut row: Vec<f64> = (0..m).map(|c| cols[c][r]).collect();
                row.extend((0..m).map(|c| if c == r { 1.0 } else { 0.0 }));
                row
            })
            .collect();
        for c in 0..m {
            let piv = (c..m)
                .max_by(|&x, &y| a[x][c].abs().total_cmp(&a[y][c].abs()))
                .unwrap();
            a.swap(c, piv);
            let d = a[c][c];
            for v in a[c].iter_mut() {
                *v /= d;
            }
            for r in 0..m {
                if r != c {
                    let f = a[r][c];
                    if f != 0.0 {
                        let src = a[c].clone();
                        for (v, s) in a[r].iter_mut().zip(&src) {
                            *v -= f * s;
                        }
                    }
                }
            }
        }
        (0..m)
            .map(|k| a[k][m..].iter().map(|v| v * v).sum())
            .collect()
    }

    /// Task 2: the weights a `Fresh` warm solve exports from a slack start
    /// (where unit weights ARE exact) equal the exact `‖e_kᵀB⁻¹‖²` of the
    /// final basis, and carrying them into a child solve keeps them exact.
    #[test]
    fn warm_dse_weights_stay_exact_and_carry() {
        // min x0 + 2 x1 + 3 x2, three G rows: the all-slack basis is dual
        // feasible (c >= 0) and primal infeasible, so the dual pivots.
        let lp = Lp {
            num_col: 3,
            num_row: 3,
            a: cm(&[&[1.0, 1.0, 1.0], &[2.0, 1.0, 0.0], &[0.0, 1.0, 2.0]], 3),
            rhs: vec![4.0, 3.0, 2.0],
            objective: vec![1.0, 2.0, 3.0],
            maximize: false,
            senses: "GGG".chars().collect(),
            lower: vec![0.0; 3],
            upper: vec![10.0; 3],
        };
        let slack = Basis {
            col_status: vec![BasisStatus::AtLower; 3],
            row_status: vec![BasisStatus::Basic; 3],
        };
        let r = reoptimize_dse(&lp, &cfg_np(), &slack, None, None, DseMode::Fresh).unwrap();
        assert!(r.warm && r.status == ReoptStatus::Optimal);
        let b = r.basis.clone().unwrap();
        let w = r.dse_weights.clone().expect("weights at Optimal");
        let exact = exact_dse(&lp, &b);
        assert_eq!(w.basic.len(), exact.len());
        for (g, e) in w.basic.iter().zip(&exact) {
            assert!(
                (*g as f64 - e).abs() <= 1e-4 * (1.0 + e),
                "weights {:?} vs exact {exact:?}",
                w.basic
            );
        }
        // Child: tighten x0 <= 0.5; carry the parent's weights.
        let mut child = lp.clone();
        child.upper[0] = 0.5;
        let rc = reoptimize_dse(&child, &cfg_np(), &b, None, None, DseMode::Carry(&w)).unwrap();
        assert!(rc.warm && rc.status == ReoptStatus::Optimal);
        let truth = cold_solve(&child);
        assert!((rc.objective - truth.objective).abs() < 1e-9);
        let wc = rc.dse_weights.expect("child weights");
        let exact_c = exact_dse(&child, &rc.basis.unwrap());
        for (g, e) in wc.basic.iter().zip(&exact_c) {
            assert!(
                (*g as f64 - e).abs() <= 1e-4 * (1.0 + e),
                "child weights {:?} vs exact {exact_c:?}",
                wc.basic
            );
        }
    }

    #[test]
    fn dual_infeasible_start_falls_back_cold() {
        // Objective change (not a bound change) breaks the dual-feasibility
        // premise: the honest outcome is a cold fallback — right answer,
        // warm = false.
        let mut lp = node_root();
        let root = cold_solve(&lp);
        let b = root.basis.clone().expect("root export");
        lp.objective = vec![-1.0, -0.5]; // y no longer worth its rows
        let r = reoptimize(&lp, &cfg_np(), &b, None, None).unwrap();
        assert!(!r.warm, "cost change must be rejected as a dual start");
        assert_eq!(r.status, ReoptStatus::Optimal);
        let truth = cold_solve(&lp);
        assert!(
            (r.objective - truth.objective).abs() < 1e-9,
            "fallback {} vs cold truth {}",
            r.objective,
            truth.objective
        );
    }

    #[test]
    fn reopt_proves_infeasible_child() {
        // x ≥ 4 AND y ≤ 1: 3·4 + 4·1 = 16 > 12 — infeasible node. The dual
        // engine proves it from the parent basis; the bound convention for
        // a min problem is +inf (always prune).
        let mut lp = node_root();
        let root = cold_solve(&lp);
        let b = root.basis.clone().expect("root export");
        lp.lower[0] = 4.0;
        lp.upper[1] = 1.0;
        let r = reoptimize(&lp, &cfg_np(), &b, None, None).unwrap();
        assert_eq!(r.status, ReoptStatus::Infeasible);
        assert!(r.warm, "the premise held; infeasibility is the verdict");
        assert_eq!(r.obj_bound, f64::INFINITY);
    }

    #[test]
    fn reopt_infeasible_child_carries_verified_farkas_ray() {
        // An Infeasible verdict from the WARM path must come with a
        // proof — same node as `reopt_proves_infeasible_child`, and the
        // `warm` assertion below pins that the ray is the dual engine's own
        // empty-event-set certificate, not the primal fallback's.
        let mut lp = node_root();
        let root = cold_solve(&lp);
        let b = root.basis.clone().expect("root export");
        lp.lower[0] = 4.0;
        lp.upper[1] = 1.0;
        let r = reoptimize(&lp, &cfg_np(), &b, None, None).unwrap();
        assert_eq!(r.status, ReoptStatus::Infeasible);
        assert!(r.warm, "ray below must be the dual engine's");
        let y = r.dual_ray.expect("infeasible reopt carries a Farkas ray");
        assert_eq!(y.len(), lp.num_row);
        assert!(y.iter().any(|&v| v != 0.0), "ray must be nonzero");
        assert!(farkas_ok(&lp, &y), "ray must verify: {y:?}");
    }

    /// File-based repro: a warm node reopt of a
    /// tiny market_split (medium seed 1) node LP — 26 columns, 6 rows (3 E),
    /// a dual-degenerate basis exported by the cold primal — stalled to the
    /// MIP's 200-iteration warm limit while no cold solve needed > 100.
    /// Contract: the warm reopt of the captured LP+basis must land Optimal
    /// warm, at the cold objective, in < 100 iterations.
    #[test]
    fn c2_warm_reopt_of_market_split_stall_lands_fast() {
        // Cold truth (primal, presolve off): the objective the MIP re-solve
        // used after falling back.
        let (lp, basis) = load_repro("warm_stall_market_split_s1");
        let truth = cold_solve(&lp);
        assert_eq!(truth.status, SimplexStatus::Optimal);
        let r = reoptimize(&lp, &cfg_np(), &basis, None, Some(200)).unwrap();
        assert!(r.warm, "the captured snapshot must import warm");
        assert_eq!(r.status, ReoptStatus::Optimal);
        eprintln!("c2 warm reopt: {} iters", r.iters);
        assert!(
            (r.objective - truth.objective).abs() <= 1e-6 * truth.objective.abs().max(1.0),
            "warm obj {} != cold truth {}",
            r.objective,
            truth.objective
        );
        assert!(
            r.iters < 100,
            "warm reopt took {} iters (must be < 100)",
            r.iters
        );
    }

    /// The three captured root-cut-round reopts of MIPLIB `pg`
    /// (`benches/data/repros/warm_stall_pg_cut{0,1,2}`, README there). At
    /// capture time every one was imported, dropped at the distress rung,
    /// and re-solved cold (4.5k-7.4k iterations vs HiGHS's 113-129 warm).
    /// Mechanism: a boxed flip overshooting the pinned row over its far
    /// bound (see the ratio-test comment). Acceptance: each lands WARM at
    /// the cold objective.
    #[test]
    fn c7_warm_reopt_of_pg_cut_rounds_lands_warm() {
        for name in [
            "warm_stall_pg_cut0",
            "warm_stall_pg_cut1",
            "warm_stall_pg_cut2",
        ] {
            let (lp, basis) = load_repro(name);
            let truth = cold_solve(&lp);
            assert_eq!(truth.status, SimplexStatus::Optimal, "{name} cold");
            let r = reoptimize(&lp, &cfg_np(), &basis, None, None).unwrap();
            eprintln!(
                "c7 {name}: warm={} status={:?} iters={} obj={} truth={}",
                r.warm, r.status, r.iters, r.objective, truth.objective
            );
            assert_eq!(r.status, ReoptStatus::Optimal, "{name}");
            assert!(r.warm, "{name}: must land WARM, not fall back cold");
            assert!(
                (r.objective - truth.objective).abs() <= 1e-6 * truth.objective.abs().max(1.0),
                "{name}: warm obj {} != cold truth {}",
                r.objective,
                truth.objective
            );
            assert!(
                r.iters <= 300,
                "{name}: warm reopt took {} iters (must be <= 300)",
                r.iters
            );
        }
    }

    /// End-to-end proof: three nodes of a tiny MIP's LP relaxations —
    /// root (cold primal) → branch on the first fractional column → both
    /// children reoptimized WARM from the root basis, plus the cutoff
    /// prune. This is the number the MIP search plans node budgets around.
    #[test]
    fn s4_three_node_mip_ladder() {
        // Root LP: max 3x0 + 5x1 + 4x2 + 2x3  s.t.
        //   x0 +  x1 + 2x2 +  x3 ≤ 12
        //  4x0 + 3x1 +  x2 + 2x3 ≤ 18
        //  2x0 + 5x1 + 3x2 +  x3 ≤ 15          0 ≤ x ≤ 6
        // Root optimum (0.9, 0, 2.4, 6) — fractional, branchable.
        let root_lp = Lp {
            num_col: 4,
            num_row: 3,
            a: cm(
                &[
                    &[1.0, 1.0, 2.0, 1.0],
                    &[4.0, 3.0, 1.0, 2.0],
                    &[2.0, 5.0, 3.0, 1.0],
                ],
                4,
            ),
            rhs: vec![12.0, 18.0, 15.0],
            objective: vec![-3.0, -5.0, -4.0, -2.0],
            maximize: false,
            senses: "LLL".chars().collect(),
            lower: vec![0.0; 4],
            upper: vec![6.0; 4],
        };
        let root = cold_solve(&root_lp);
        assert_eq!(root.status, SimplexStatus::Optimal);
        let root_basis = root.basis.clone().expect("root basis");
        let branch = (0..4).find(|&j| {
            let v = root.x[j];
            (v - v.round()).abs() > 1e-6
                && v.round() >= root_lp.lower[j]
                && v.round() <= root_lp.upper[j]
                && v.round() != v
        });
        let j = branch.expect("root relaxation must be fractional");
        let floor_j = root.x[j].floor();
        let ceil_j = root.x[j].ceil();

        let mut child_a = root_lp.clone();
        child_a.upper[j] = floor_j;
        let mut child_b = root_lp.clone();
        child_b.lower[j] = ceil_j;

        // Truth: cold solves of both children.
        let truth_a = cold_solve(&child_a);
        let truth_b = cold_solve(&child_b);
        assert_eq!(truth_a.status, SimplexStatus::Optimal);
        assert_eq!(truth_b.status, SimplexStatus::Optimal);

        // Warm dual reopts from the ROOT basis.
        let ra = reoptimize(&child_a, &cfg_np(), &root_basis, None, None).unwrap();
        let rb = reoptimize(&child_b, &cfg_np(), &root_basis, None, None).unwrap();
        assert!(ra.warm && rb.warm, "bound changes keep dual feasibility");
        assert_eq!(ra.status, ReoptStatus::Optimal);
        assert_eq!(rb.status, ReoptStatus::Optimal);
        assert!(
            (ra.objective - truth_a.objective).abs() < 1e-9,
            "A warm {} vs truth {}",
            ra.objective,
            truth_a.objective
        );
        assert!(
            (rb.objective - truth_b.objective).abs() < 1e-9,
            "B warm {} vs truth {}",
            rb.objective,
            truth_b.objective
        );
        let warm_iters = ra.iters + rb.iters;
        let cold_iters = (truth_a.phase1_iterations + truth_a.phase2_iterations)
            + (truth_b.phase1_iterations + truth_b.phase2_iterations);
        assert!(
            ra.iters < truth_a.phase1_iterations + truth_a.phase2_iterations
                && rb.iters < truth_b.phase1_iterations + truth_b.phase2_iterations,
            "per-child warm ({}, {}) must beat cold ({}, {})",
            ra.iters,
            rb.iters,
            truth_a.phase1_iterations + truth_a.phase2_iterations,
            truth_b.phase1_iterations + truth_b.phase2_iterations
        );
        eprintln!(
            "warm-vs-cold: warm {warm_iters} pivots vs cold {cold_iters} ⇒ {:.1}× cheaper ({} + {} vs {} + {})",
            cold_iters as f64 / warm_iters as f64,
            ra.iters,
            rb.iters,
            truth_a.phase1_iterations + truth_a.phase2_iterations,
            truth_b.phase1_iterations + truth_b.phase2_iterations,
        );

        // Cutoff: incumbent = the better child's optimum minus a hair —
        // the WORSE child must Pruned-fire (min sense: bound ≥ cutoff),
        // the better child must still solve Optimal below it.
        assert!(
            (truth_a.objective - truth_b.objective).abs() > 1e-4,
            "test assumption: children differ"
        );
        let (worse_lp, worse_truth, other_lp, other_truth) =
            if truth_a.objective > truth_b.objective {
                (&child_a, &truth_a, &child_b, &truth_b)
            } else {
                (&child_b, &truth_b, &child_a, &truth_a)
            };
        let cutoff = worse_truth.objective - 1e-6;
        let rw = reoptimize(worse_lp, &cfg_np(), &root_basis, Some(cutoff), None).unwrap();
        assert_eq!(
            rw.status,
            ReoptStatus::Pruned,
            "cutoff below the worse child's optimum"
        );
        assert!(rw.obj_bound >= cutoff - 1e-9, "obj_bound {}", rw.obj_bound);
        let ro = reoptimize(other_lp, &cfg_np(), &root_basis, Some(cutoff), None).unwrap();
        assert_eq!(ro.status, ReoptStatus::Optimal);
        assert!(
            (ro.objective - other_truth.objective).abs() < 1e-9,
            "promising child still solves exactly: {} vs {}",
            ro.objective,
            other_truth.objective
        );

        // Grandchild wiring: the child's own exported basis re-enters with
        // zero pivots on the unchanged LP.
        let g = reoptimize(
            &child_a,
            &cfg_np(),
            ra.basis.as_ref().expect("child basis"),
            None,
            None,
        )
        .unwrap();
        assert!(g.warm);
        assert_eq!(g.status, ReoptStatus::Optimal);
        assert_eq!(g.iters, 0, "re-root on the child's own basis: no pivots");
        assert!((g.objective - ra.objective).abs() < 1e-12);
    }

    // ---- The incremental row handle (ReoptSession) ----

    /// Cut-loop target LP (the s4 ladder root): max 3x0+5x1+4x2+2x3, three 'L'
    /// rows, 0 ≤ x ≤ 6 (min form: costs negated). Root optimum (0.9, 0, 2.4, 6)
    /// at obj −24.3 with row 0 LOOSE (its slack basic), rows 1/2 tight, x1 at
    /// its lower bound, x3 at its upper bound — i.e. the loose-row removal case
    /// and a fractional, cuttable relaxation in one LP.
    fn cut_lp() -> Lp {
        Lp {
            num_col: 4,
            num_row: 3,
            a: cm(
                &[
                    &[1.0, 1.0, 2.0, 1.0],
                    &[4.0, 3.0, 1.0, 2.0],
                    &[2.0, 5.0, 3.0, 1.0],
                ],
                4,
            ),
            rhs: vec![12.0, 18.0, 15.0],
            objective: vec![-3.0, -5.0, -4.0, -2.0],
            maximize: false,
            senses: "LLL".chars().collect(),
            lower: vec![0.0; 4],
            upper: vec![6.0; 4],
        }
    }

    /// Degenerate LP for the TIGHT-row removal case: max x + 2y s.t. x+y ≤ 4,
    /// y ≤ 3, x ≤ 1 ⇒ the optimum (1,3) is where all three rows meet, so the
    /// basic solution carries a slack pinned at 0 — the variable the repair can
    /// shed without moving the point.
    fn degenerate_lp() -> Lp {
        Lp {
            num_col: 2,
            num_row: 3,
            a: cm(&[&[1.0, 1.0], &[0.0, 1.0], &[1.0, 0.0]], 2),
            rhs: vec![4.0, 3.0, 1.0],
            objective: vec![-1.0, -2.0],
            maximize: false,
            senses: "LLL".chars().collect(),
            lower: vec![0.0; 2],
            upper: vec![10.0; 2],
        }
    }

    /// The MIP search's side of the tableau-row indexing contract: caller identity ⇒
    /// standard slot, computed ONLY from the [`Lp`] and [`StdLayout`] (i.e. the
    /// mapping the tableau_row docs promise is derivable without internals —
    /// this helper is the proof of that claim, not a test-only oracle).
    fn std_slot_of(lp: &Lp, layout: StdLayout, b: TableauBasic) -> usize {
        let free_before = |j: usize| {
            (0..j)
                .filter(|&k| lp.lower[k] == f64::NEG_INFINITY && lp.upper[k] == f64::INFINITY)
                .count()
        };
        match b {
            TableauBasic::Column(j) => j + free_before(j),
            TableauBasic::FreeHalf { col, minus } => col + free_before(col) + minus as usize,
            TableauBasic::Slack(i) => {
                layout.num_structural + (0..i).filter(|&k| lp.senses[k] != 'E').count()
            }
            TableauBasic::Artificial(i) => {
                layout.num_structural
                    + layout.num_slack
                    + (0..i).filter(|&k| lp.senses[k] == 'E').count()
            }
        }
    }

    #[test]
    fn tableau_rows_match_scratch_binv_a() {
        // Agreement: every tableau row of the session basis against
        // e_posᵀ B⁻¹ A rebuilt from scratch — the DENSE LU (a different
        // factorization code path than the session's product form) solving
        // Bᵀ y = e_pos, then a dense dot with the rebuilt standard columns.
        let lp = kitchen_sink(false);
        let root = cold_solve(&lp);
        assert_eq!(root.status, SimplexStatus::Optimal, "test LP must solve");
        let snap = root.basis.clone().expect("root basis");
        let mut session = ReoptSession::new(lp.clone(), &cfg_np(), Some(snap.clone()));
        use crate::dense::DenseLu;

        let st = to_dual_standard(&lp, &cfg_np()).unwrap();
        let (nb, nat) = import_basis_dual(&lp, &st, &snap).unwrap();
        let m = st.m;
        let bmat = DenseMatrix::from_columns(&(0..m).map(|i| dcol(&st, nb[i])).collect::<Vec<_>>());
        let lu = DenseLu::factor(&bmat, 1e-12).expect("a basis is nonsingular");
        // Reference basic values: x_b = B⁻¹(b − Σ_{nonbasic at ũ} ũ·a).
        let mut in_b = vec![false; st.ncol()];
        for &v in &nb {
            in_b[v] = true;
        }
        let mut rhs = st.b.clone();
        for (v, &ib) in in_b.iter().enumerate() {
            if !ib && nat[v] && st.ub[v].is_finite() {
                for k in st.csc_ptr[v]..st.csc_ptr[v + 1] {
                    rhs[st.csc_row[k]] -= st.ub[v] * st.csc_val[k];
                }
            }
        }
        let mut x_b_ref = vec![0.0; m];
        lu.solve(&rhs, &mut x_b_ref);
        // The standard point of that basis (for the row-content identity).
        let mut x_tilde = vec![0.0; st.ncol()];
        for (i, &v) in nb.iter().enumerate() {
            x_tilde[v] = x_b_ref[i];
        }
        for (v, &at_u) in nat.iter().enumerate() {
            if at_u && !in_b[v] {
                x_tilde[v] = st.ub[v];
            }
        }
        let mut binv_b = vec![0.0; m];
        lu.solve(&st.b, &mut binv_b);

        for pos in 0..m {
            let row = session.tableau_row(pos).expect("tableau row in range");
            assert_eq!(row.coeffs.len(), st.ncol());
            let mut e = vec![0.0; m];
            e[pos] = 1.0;
            let mut y = vec![0.0; m];
            lu.solve_transpose(&e, &mut y);
            for v in 0..st.ncol() {
                let want: f64 = (st.csc_ptr[v]..st.csc_ptr[v + 1])
                    .map(|k| st.csc_val[k] * y[st.csc_row[k]])
                    .sum();
                let tol = 1e-9 * (1.0 + want.abs());
                assert!(
                    (row.coeffs[v] - want).abs() <= tol,
                    "row {pos} col {v}: {} vs scratch {want}",
                    row.coeffs[v]
                );
            }
            // Row content identity: āᵀx̃ = (B⁻¹b)_pos for the FULL standard
            // point x̃ (basics at x_b, nonbasics at 0 or ũ) — the tableau row
            // reproduces the basic value from the nonbasic contributions.
            let dot: f64 = row.coeffs.iter().zip(&x_tilde).map(|(a, xv)| a * xv).sum();
            assert!(
                (dot - binv_b[pos]).abs() <= 1e-8 * (1.0 + binv_b[pos].abs()),
                "row {pos}: āᵀx̃ {dot} vs (B⁻¹b) {}",
                binv_b[pos]
            );
            assert!(
                (row.basic_value - x_b_ref[pos]).abs() <= 1e-8 * (1.0 + x_b_ref[pos].abs()),
                "row {pos}: basic_value {} vs scratch {}",
                row.basic_value,
                x_b_ref[pos]
            );
        }
    }

    #[test]
    fn tableau_row_basic_identity_and_caller_mapping() {
        // The two properties the Gomory generator leans on: the tableau is the
        // identity on its own basic columns, and every `basic` maps to a slot
        // the caller can find from its own Basis snapshot.
        let lp = kitchen_sink(false);
        let root = cold_solve(&lp);
        let snap = root.basis.clone().expect("root basis");
        let mut session = ReoptSession::new(lp.clone(), &cfg_np(), Some(snap.clone()));
        // kitchen_sink: x0 [0,10], x1 [2,inf), x2 flip, x3 free (split ⇒ 2
        // slots), x4 [0,inf) ⇒ 6 structural; rows L/G/E ⇒ 2 slacks + 1
        // artificial.
        let layout = session.layout();
        assert_eq!(
            layout,
            StdLayout {
                num_structural: 6,
                num_slack: 2,
                num_artificial: 1
            }
        );
        let m = lp.num_row;
        assert!(session.tableau_row(m).is_none(), "position out of range");
        let rows: Vec<TableauRow> = (0..m)
            .map(|pos| session.tableau_row(pos).expect("in range"))
            .collect();
        for (pos, row) in rows.iter().enumerate() {
            assert_eq!(row.position, pos);
            for (q, other) in rows.iter().enumerate() {
                let slot = std_slot_of(&lp, layout, other.basic);
                if pos != q {
                    assert_ne!(
                        slot,
                        std_slot_of(&lp, layout, row.basic),
                        "two basic positions share a standard column"
                    );
                }
                let want = if pos == q { 1.0 } else { 0.0 };
                assert!(
                    (row.coeffs[slot] - want).abs() <= 1e-9,
                    "row {pos} at basic slot {slot} ({} of row {q}): {} want {want}",
                    row.position,
                    row.coeffs[slot]
                );
            }
            // The caller-side identity is consistent with the snapshot.
            match row.basic {
                TableauBasic::Column(j) => assert_eq!(
                    snap.col_status[j],
                    BasisStatus::Basic,
                    "tableau says column {j} is basic"
                ),
                TableauBasic::Slack(i) => assert_eq!(
                    snap.row_status[i],
                    BasisStatus::Basic,
                    "tableau says row {i}'s slack is basic"
                ),
                TableauBasic::Artificial(i) => {
                    assert_eq!(lp.senses[i], 'E');
                    // A basic 'E'-row artificial is now VISIBLE — the
                    // snapshot says Basic for exactly that row.
                    assert_eq!(
                        snap.row_status[i],
                        BasisStatus::Basic,
                        "tableau says row {i}'s artificial is basic"
                    );
                }
                TableauBasic::FreeHalf { col, .. } => {
                    assert!(lp.lower[col] == f64::NEG_INFINITY && lp.upper[col] == f64::INFINITY)
                }
            }
        }
        // Every Basic the snapshot can see appears in exactly one tableau row.
        // (A fully-free column's Basic status maps to a FreeHalf slot, and a
        // basic 'E'-row artificial appears as an Artificial row — the
        // E-row entry above checks it is exactly the snapshot's Basic E rows — so
        // both are excluded from this tally.)
        let seen: Vec<TableauBasic> = rows.iter().map(|r| r.basic).collect();
        let visible = (0..lp.num_col)
            .filter(|&j| {
                snap.col_status[j] == BasisStatus::Basic
                    && !(lp.lower[j] == f64::NEG_INFINITY && lp.upper[j] == f64::INFINITY)
            })
            .count()
            + (0..lp.num_row)
                .filter(|&i| lp.senses[i] != 'E' && snap.row_status[i] == BasisStatus::Basic)
                .count();
        assert_eq!(
            seen.iter()
                .filter(|b| matches!(b, TableauBasic::Column(_) | TableauBasic::Slack(_)))
                .count(),
            visible,
            "snapshot Basic entries vs tableau rows (a basic artificial the snapshot \
             cannot see is the counting contract, not a mismatch)"
        );
    }

    #[test]
    fn b6b_binv_row_and_basic_position_agree_with_the_tableau_scan() {
        // Agreement, bit-exact: `binv_row(pos)` holds the weights the
        // `tableau_row` slack/artificial slots hold (`'L'` slot `+y_i`,
        // `'G'` slot `−y_i`, `'E'` artificial slot `+y_i` — the slack column
        // is a single ±1 entry, so the slot dot returns `±y_i` exactly), and
        // `basic_position` is exactly what a walk of `tableau_row(pos).basic`
        // finds. This is the invariant the MIP engine's per-candidate
        // `tableau_weights` relies on to skip the scan.
        let lp = kitchen_sink(false);
        let root = cold_solve(&lp);
        assert_eq!(root.status, SimplexStatus::Optimal, "test LP must solve");
        let snap = root.basis.clone().expect("root basis");
        let mut session = ReoptSession::new(lp.clone(), &cfg_np(), Some(snap.clone()));
        assert!(session.warm_available());
        let layout = session.layout();
        let m = lp.num_row;
        assert!(session.binv_row(m).is_none(), "position out of range");
        for pos in 0..m {
            let y = session.binv_row(pos).expect("in range");
            assert_eq!(y.len(), m);
            let row = session.tableau_row(pos).expect("in range");
            for (i, &yi) in y.iter().enumerate() {
                let (slot, want) = if lp.senses[i] == 'E' {
                    (std_slot_of(&lp, layout, TableauBasic::Artificial(i)), yi)
                } else {
                    let slot = std_slot_of(&lp, layout, TableauBasic::Slack(i));
                    (slot, if lp.senses[i] == 'L' { yi } else { -yi })
                };
                assert_eq!(row.coeffs[slot], want, "row {pos} slot {slot} (row {i})");
            }
        }
        // basic_position ≡ the scan, both directions, on every column.
        let mut scanned = vec![usize::MAX; lp.num_col];
        for pos in 0..m {
            if let TableauBasic::Column(j) = session.tableau_row(pos).expect("in range").basic {
                scanned[j] = pos;
            }
        }
        for (j, &sp) in scanned.iter().enumerate() {
            assert_eq!(
                session.basic_position(j),
                (sp != usize::MAX).then_some(sp),
                "column {j}"
            );
        }
        // A nonbasic column and an out-of-range query answer None.
        let nonbasic = (0..lp.num_col).find(|&j| scanned[j] == usize::MAX);
        if let Some(j) = nonbasic {
            assert_eq!(session.basic_position(j), None);
        }
        assert_eq!(session.basic_position(lp.num_col), None);
    }

    #[test]
    fn b6b_position_map_tracks_the_basis_after_row_edits_and_reoptimization() {
        // The map lives with the cached warm state, which every add_rows /
        // remove_rows / reoptimize invalidates — so it must answer for the
        // CURRENT basis, never a stale position. Cross-check against a fresh
        // scan of tableau rows after each mutation.
        let lp = kitchen_sink(false);
        let root = cold_solve(&lp);
        let snap = root.basis.clone().expect("root basis");
        let mut session = ReoptSession::new(lp.clone(), &cfg_np(), Some(snap));
        let check = |s: &mut ReoptSession, tag: &str| {
            let mut scanned = vec![usize::MAX; lp.num_col];
            for pos in 0..s.num_rows() {
                if let Some(r) = s.tableau_row(pos)
                    && let TableauBasic::Column(j) = r.basic
                {
                    scanned[j] = pos;
                }
            }
            for (j, &sp) in scanned.iter().enumerate() {
                assert_eq!(
                    s.basic_position(j),
                    (sp != usize::MAX).then_some(sp),
                    "{tag}: column {j}"
                );
            }
        };
        check(&mut session, "root");
        session
            .add_rows(&[CutRow::new(vec![(0, 1.0), (4, 2.0)], 'L', 3.0)])
            .expect("append a cut row");
        check(
            &mut session,
            "after add_rows (stale basis, not yet reoptimized)",
        );
        let r = session.reoptimize(None, None).expect("reoptimize");
        assert_eq!(r.status, ReoptStatus::Optimal);
        check(&mut session, "after reoptimize");
    }

    #[test]
    fn session_add_cut_reopts_warm_and_beats_cold() {
        // (a) one violated cut appended to a solved LP. The slack enters
        // basic, the duals are the old ones with a 0 appended, so the reopt is
        // warm and cheaper than a cold solve of the augmented LP — with the
        // same objective.
        let lp = cut_lp();
        let root = cold_solve(&lp);
        assert_eq!(root.status, SimplexStatus::Optimal);
        let sum: f64 = root.x.iter().sum();
        assert!(
            sum > 9.0 + 1e-6,
            "the crafted cut must be violated by the root point: {sum}"
        );
        let mut s = ReoptSession::new(lp.clone(), &cfg_np(), root.basis.clone());
        s.add_rows(&[CutRow::new(
            vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)],
            'L',
            9.0,
        )])
        .unwrap();
        assert_eq!(s.num_rows(), lp.num_row + 1);
        let snap = s.basis().expect("snapshot");
        assert_eq!(snap.row_status.len(), lp.num_row + 1);
        assert_eq!(
            snap.row_status[lp.num_row],
            BasisStatus::Basic,
            "the new row's slack enters the basis"
        );
        let r = s.reoptimize(None, None).unwrap();
        assert!(r.warm, "slack-basic append preserves dual feasibility");
        assert_eq!(r.status, ReoptStatus::Optimal);
        assert!(
            r.objective > root.objective + 1e-6,
            "the cut must cut the relaxation: {} → {}",
            root.objective,
            r.objective
        );
        let truth = cold_solve(s.lp());
        assert_eq!(truth.status, SimplexStatus::Optimal);
        assert!(
            (r.objective - truth.objective).abs() < 1e-9,
            "warm {} vs cold truth {}",
            r.objective,
            truth.objective
        );
        let cold_iters = truth.phase1_iterations + truth.phase2_iterations;
        assert!(
            r.iters < cold_iters,
            "cut-loop reopt {} iters must beat the cold solve of the augmented LP ({cold_iters})",
            r.iters
        );
    }

    #[test]
    fn session_remove_tight_cut_restores_root_optimum() {
        // (b) the cut is TIGHT at the new optimum (its slack nonbasic), so
        // deleting it is the repair case; either the shed repair holds (warm)
        // or the snapshot is dropped and the answer comes from a cold solve —
        // both must restore the root optimum exactly.
        let lp = cut_lp();
        let root = cold_solve(&lp);
        let mut s = ReoptSession::new(lp.clone(), &cfg_np(), root.basis.clone());
        s.add_rows(&[CutRow::new(
            vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)],
            'L',
            9.0,
        )])
        .unwrap();
        let cut = s.reoptimize(None, None).unwrap();
        assert_eq!(cut.status, ReoptStatus::Optimal);
        assert!(
            cut.objective > root.objective + 1e-6,
            "the cut changed the optimum: {}",
            cut.objective
        );
        assert_eq!(
            s.basis().expect("basis").row_status[3],
            BasisStatus::AtLower,
            "crafted so the cut is tight (slack nonbasic) at the new optimum"
        );
        let map = s.remove_rows(&[3]).unwrap();
        assert_eq!(
            map[3],
            rusolve_core::NO_ROW,
            "the deleted row has no new index"
        );
        assert_eq!(&map[..3], &[0, 1, 2], "surviving rows keep their indices");
        assert_eq!(s.num_rows(), 3);
        if let Some(b) = s.basis() {
            // The repair kept the snapshot: square, and the shed side recorded.
            assert_eq!(b.row_status.len(), 3);
            assert_eq!(b.col_status.len(), 4);
        }
        let back = s.reoptimize(None, None).unwrap();
        assert_eq!(back.status, ReoptStatus::Optimal);
        assert!(
            (back.objective - root.objective).abs() < 1e-9,
            "aging the cut out must restore the root optimum: {} vs {}",
            back.objective,
            root.objective
        );
        let truth = cold_solve(s.lp());
        assert!(
            (back.objective - truth.objective).abs() < 1e-9,
            "and it must match the cold truth of the reduced LP"
        );
        assert_eq!(s.lp().rhs, lp.rhs);
        assert_eq!(s.lp().senses, lp.senses);
        assert_eq!(
            s.lp().a,
            lp.a,
            "matrix content survives remove/add round trip"
        );
    }

    #[test]
    fn session_remove_loose_row_shrinks_basis_warm() {
        // (b), the clean case: the deleted row's slack was BASIC, so the
        // basic set shrinks in step with m and the reopt stays warm.
        let lp = cut_lp();
        let root = cold_solve(&lp);
        let cols0 = dense_cols(&lp);
        let lhs0: f64 = (0..lp.num_col).map(|j| cols0[j][0] * root.x[j]).sum();
        assert!(lhs0 < 12.0 - 1e-6, "row 0 must be loose: Ax = {lhs0}");
        assert_eq!(
            root.basis.as_ref().expect("basis").row_status[0],
            BasisStatus::Basic,
            "loose row ⇒ slack basic"
        );
        let mut s = ReoptSession::new(lp.clone(), &cfg_np(), root.basis.clone());
        let map = s.remove_rows(&[0]).unwrap();
        assert_eq!(map[0], rusolve_core::NO_ROW);
        assert_eq!(&map[1..], &[0, 1], "later rows shift down by one");
        assert_eq!(s.num_rows(), 2);
        let r = s.reoptimize(None, None).unwrap();
        assert!(r.warm, "a basic-slack deletion leaves a square, warm basis");
        assert_eq!(r.status, ReoptStatus::Optimal);
        let truth = cold_solve(s.lp());
        assert!(
            (r.objective - truth.objective).abs() < 1e-9,
            "warm {} vs cold truth {}",
            r.objective,
            truth.objective
        );
        assert!(
            r.objective <= root.objective + 1e-9,
            "dropping a constraint cannot worsen the relaxation"
        );
    }

    #[test]
    fn session_remove_tight_row_repairs_the_snapshot_and_matches_truth() {
        // (b), the repair case on a DEGENERATE optimum: three rows meet at the
        // optimum so a basic slack sits pinned at 0, and that is what the repair
        // sheds — the snapshot survives as a square, primal feasible basis of
        // the reduced LP. Deleting a tight row still shifts the duals of every
        // column that has a coefficient in it, so the warm entry decides whether
        // the start is usable; here it honestly declines and the (correct)
        // answer comes from the cold fallback. Both outcomes must reproduce the
        // cold truth of the reduced LP.
        let lp = degenerate_lp();
        let root = cold_solve(&lp);
        assert_eq!(root.status, SimplexStatus::Optimal);
        let cols = dense_cols(&lp);
        for i in 0..lp.num_row {
            let lhs: f64 = (0..lp.num_col).map(|j| cols[j][i] * root.x[j]).sum();
            assert!(
                (lhs - lp.rhs[i]).abs() < 1e-9,
                "row {i} must be tight at the crafted optimum: {lhs} vs {}",
                lp.rhs[i]
            );
        }
        assert_eq!(
            root.basis.as_ref().expect("basis").row_status[1],
            BasisStatus::AtLower,
            "row 1 is tight ⇒ its slack is nonbasic ⇒ a shed is needed"
        );
        let mut s = ReoptSession::new(lp.clone(), &cfg_np(), root.basis.clone());
        let map = s.remove_rows(&[1]).unwrap();
        assert_eq!(&map, &[0, rusolve_core::NO_ROW, 1]);
        let snap = s.basis().expect("the repair keeps a square snapshot");
        assert_eq!(snap.row_status.len(), 2, "one row status per surviving row");
        assert_eq!(
            snap.row_status,
            vec![BasisStatus::AtLower, BasisStatus::AtLower],
            "the degenerate basic slack (old row 2, value 0) was shed to its lower side"
        );
        assert_eq!(snap.col_status, vec![BasisStatus::Basic; 2]);
        let r = s.reoptimize(None, None).unwrap();
        assert_eq!(r.status, ReoptStatus::Optimal);
        let truth = cold_solve(s.lp());
        assert!(
            (r.objective - truth.objective).abs() < 1e-9,
            "warm/cold (warm = {}) obj {} vs truth {}",
            r.warm,
            r.objective,
            truth.objective
        );
    }

    #[test]
    fn session_add_e_row_reopts_to_truth() {
        // (c) an 'E' row has no slack — its fixed artificial enters basic,
        // the deficit the Basis counting contract already covers, and the
        // artificial's value IS the equality violation the dual loop repairs.
        let lp = cut_lp();
        let root = cold_solve(&lp);
        let mut s = ReoptSession::new(lp.clone(), &cfg_np(), root.basis.clone());
        s.add_rows(&[CutRow::new(vec![(0, 1.0), (3, -1.0)], 'E', 0.0)])
            .unwrap();
        assert_eq!(s.layout().num_artificial, 1, "the new row is an E row");
        assert_eq!(s.basis().expect("basis").row_status[3], BasisStatus::Basic);
        let r = s.reoptimize(None, None).unwrap();
        assert_eq!(r.status, ReoptStatus::Optimal);
        let truth = cold_solve(s.lp());
        assert!(
            (r.objective - truth.objective).abs() < 1e-9,
            "warm {} vs cold truth {}",
            r.objective,
            truth.objective
        );
        assert!(
            (r.x[0] - r.x[3]).abs() < 1e-7,
            "the equality must hold at the warm point: {:?}",
            r.x
        );
    }

    #[test]
    fn a_snapshot_that_cannot_describe_the_lp_is_dropped_not_carried() {
        // The lp_mut escape hatch lets a caller change the shape under the
        // session; the contract is that this degrades to COLD, never to a panic
        // or to a warm start on a basis that cannot be indexed.
        let lp = cut_lp();
        let root = cold_solve(&lp);
        let mut s = ReoptSession::new(lp.clone(), &cfg_np(), root.basis.clone());
        assert!(s.basis().is_some());
        // A legitimate 4th row added UNDER the session (bypassing add_rows): the
        // snapshot no longer has a status for it, so it is dropped, not indexed
        // out of bounds.
        {
            let l = s.lp_mut();
            let n = l.num_col;
            let mut b = rusolve_core::CscBuilder::from_matrix(&l.a);
            b.push_row(&(0..n).map(|j| (j, 1.0)).collect::<Vec<_>>());
            l.a = b.into_matrix();
            l.rhs.push(5.0);
            l.senses.push('L');
            l.num_row += 1;
        }
        let map = s.remove_rows(&[0]).unwrap();
        assert_eq!(map.len(), 4);
        assert!(
            s.basis().is_none(),
            "the snapshot no longer describes the LP: dropped, not indexed out of bounds"
        );
        let r = s.reoptimize(None, None).unwrap();
        assert!(!r.warm, "cold by construction");
        let truth = cold_solve(s.lp());
        assert!(
            (r.objective - truth.objective).abs() < 1e-9,
            "cold truth on the edited LP: {} vs {}",
            r.objective,
            truth.objective
        );
        // Same on the column side, and a consistent snapshot survives.
        let mut narrow = root.basis.clone().expect("root basis");
        narrow.col_status.pop();
        let s2 = ReoptSession::new(lp.clone(), &cfg_np(), Some(narrow));
        assert!(
            s2.basis().is_none(),
            "a short snapshot cannot describe the LP: dropped on entry"
        );
        assert_eq!(s2.lp().a, lp.a, "dropping the basis touches no LP data");
        let mut ok = ReoptSession::new(lp.clone(), &cfg_np(), root.basis.clone());
        ok.lp_mut().objective[0] = -2.0;
        ok.remove_rows(&[0]).unwrap();
        assert!(ok.basis().is_some(), "a consistent snapshot survives");
        let r = ok.reoptimize(None, None).unwrap();
        assert_eq!(r.status, ReoptStatus::Optimal);
        let truth = cold_solve(ok.lp());
        assert!(
            (r.objective - truth.objective).abs() < 1e-9,
            "warm {} vs cold truth {}",
            r.objective,
            truth.objective
        );
    }

    #[test]
    fn add_rows_validates_the_whole_batch_before_touching_the_lp() {
        let lp = cut_lp();
        let mut s = ReoptSession::new(lp.clone(), &cfg_np(), None);
        let good = CutRow::new(vec![(0, 1.0)], 'L', 1.0);
        let bads = [
            CutRow::new(vec![(0, 1.0)], 'X', 1.0), // unknown sense
            CutRow::new(vec![(9, 1.0)], 'L', 1.0), // column out of range
            CutRow::new(vec![(0, 1.0), (0, 2.0)], 'L', 1.0), // duplicate column
            CutRow::new(vec![(0, 1.0)], 'L', f64::NAN), // non-finite rhs
            CutRow::new(vec![(0, 1.0)], 'L', f64::INFINITY), // non-finite rhs
        ];
        for bad in &bads {
            assert!(
                s.add_rows(&[good.clone(), bad.clone()]).is_err(),
                "must reject {:?}",
                bad.sense
            );
            assert_eq!(s.num_rows(), lp.num_row, "a rejected batch applies nothing");
            assert_eq!(s.lp().a, lp.a, "and never touches the matrix");
        }
        s.add_rows(&[good, CutRow::new(vec![(1, -1.0)], 'G', -2.0)])
            .unwrap();
        assert_eq!(s.num_rows(), lp.num_row + 2);
        let cols = dense_cols(s.lp());
        assert_eq!(cols[0], vec![1.0, 4.0, 2.0, 1.0, 0.0]);
        assert_eq!(cols[1], vec![1.0, 3.0, 5.0, 0.0, -1.0]);
        assert_eq!(s.lp().rhs[lp.num_row..], vec![1.0, -2.0]);
        assert_eq!(s.lp().senses[lp.num_row + 1], 'G');
        assert!(
            s.remove_rows(&[lp.num_row + 2]).is_err(),
            "remove_rows validates its indices too"
        );
    }
    /// Cold solves run on perturbed costs (HiGHS-style up-front
    /// perturbation, ~5e-7 relative). The reported objective, primal point
    /// and duals must be the TRUE-cost ones: the cleanup removes the
    /// perturbation. A degenerate 4x4 assignment LP with integer costs and
    /// many ties makes the perturbation active; any leak shows up as an
    /// objective off an integer or a reduced-cost identity error far above
    /// 1e-9. (Differential coverage vs HiGHS: tools/diff/dual_diff.py.)
    #[test]
    fn perturbed_cold_dual_reports_true_cost_optimum_and_duals() {
        let n = 4;
        let nv = n * n;
        let cost: Vec<f64> = (0..nv).map(|k| (((k / n) * (k % n)) % 3) as f64).collect();
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for i in 0..n {
            let mut r = vec![0.0; nv];
            for j in 0..n {
                r[i * n + j] = 1.0;
            }
            rows.push(r);
        }
        for j in 0..n {
            let mut r = vec![0.0; nv];
            for i in 0..n {
                r[i * n + j] = 1.0;
            }
            rows.push(r);
        }
        let refs: Vec<&[f64]> = rows.iter().map(|r| r.as_slice()).collect();
        let lp = lp6(
            &refs,
            vec![1.0; 2 * n],
            &"E".repeat(2 * n),
            cost.clone(),
            vec![0.0; nv],
            vec![f64::INFINITY; nv],
        );
        let s = solve_dual(&lp);
        assert_eq!(s.status, SimplexStatus::Optimal);
        // True optimum: a permutation with i*j mod 3 == 0 everywhere exists
        // (identity on row/col 0 plus pairs hitting 0 or 3), cost 0 or a
        // small integer; the exact value is whatever primal says.
        let p = SimplexSolver::new(SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert!(
            (s.objective - p.objective).abs() < 1e-9,
            "{} vs {}",
            s.objective,
            p.objective
        );
        assert!((s.objective - s.objective.round()).abs() < 1e-9);
        assert!(s.duals_valid);
        for j in 0..nv {
            let aj: f64 = (0..2 * n).map(|r| s.row_duals[r] * rows[r][j]).sum();
            let d = cost[j] - aj;
            assert!((d - s.reduced_costs[j]).abs() < 1e-9, "d identity col {j}");
            assert!(d > -1e-9, "dual feasible at true costs, col {j}: {d}");
            if s.x[j] > 1e-9 {
                assert!(d.abs() < 1e-9, "complementarity col {j}: {d}");
            }
        }
    }

    /// `refresh_dual_standard` is the exact-equal twin of a
    /// fresh `to_dual_standard` for every edit it accepts. Upper-bound-only
    /// edits leave every refreshed field BIT-identical to a rebuild (the
    /// shift vector is untouched, so `b` never re-accumulates); lower-bound
    /// edits re-accumulate `b` by deltas, so those are checked to tolerance
    /// while ũ / shift / the `s` cost-shift init stay bit-exact.
    #[test]
    fn p6_refresh_matches_rebuild_after_bound_edits() {
        let mut lp = Lp {
            num_col: 4,
            num_row: 3,
            a: cm(
                &[
                    &[1.0, 2.0, 1.0, 0.5],
                    &[2.0, 1.0, 3.0, 1.0],
                    &[1.0, 1.0, 1.0, 2.0],
                ],
                4,
            ),
            rhs: vec![10.0, 12.0, 9.0],
            objective: vec![-3.0, -5.0, -4.0, -2.0],
            maximize: false,
            senses: "LLL".chars().collect(),
            lower: vec![0.0, 0.0, 0.5, -1.0],
            upper: vec![6.0, 4.0, f64::INFINITY, 3.0],
        };
        let mut st = to_dual_standard(&lp, &cfg_np()).unwrap();
        let bits = |a: &[f64], b: &[f64]| a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits());

        // Upper-only edits (incl. +inf → finite): shift untouched ⇒ the
        // refreshed cache is bit-identical to a rebuild.
        for (j, v) in [(0, 4.0), (2, 7.5), (0, 6.0), (3, 1.25)] {
            lp.upper[j] = v;
            assert!(
                refresh_dual_standard(&lp, &cfg_np(), &mut st),
                "ub edit {j}"
            );
            let fresh = to_dual_standard(&lp, &cfg_np()).unwrap();
            assert!(bits(&st.b, &fresh.b), "b bit-exact after ub edit {j}={v}");
            assert!(bits(&st.ub, &fresh.ub), "ub after ub edit {j}={v}");
            assert!(bits(&st.shift, &fresh.shift));
            assert!(bits(&st.s, &fresh.s));
            assert!(bits(&st.pert, &fresh.pert));
            assert_eq!(st.at_ub, fresh.at_ub);
        }

        // Lower edits re-accumulate b by per-column deltas: tolerance on b,
        // bit-exactness everywhere else (shift, ũ, s are direct writes).
        for (j, v) in [(0, 1.5), (1, 0.0), (3, -0.5), (0, 0.0), (3, 0.75)] {
            lp.lower[j] = v;
            assert!(
                refresh_dual_standard(&lp, &cfg_np(), &mut st),
                "lb edit {j}"
            );
            let fresh = to_dual_standard(&lp, &cfg_np()).unwrap();
            let bmax = fresh.b.iter().fold(0.0f64, |m, &x| m.max(x.abs()));
            assert!(
                st.b.iter()
                    .zip(&fresh.b)
                    .all(|(x, y)| (x - y).abs() <= 1e-9 * bmax.max(1.0)),
                "b within tolerance after lb edit {j}={v}"
            );
            assert!(bits(&st.ub, &fresh.ub), "ub after lb edit {j}={v}");
            assert!(bits(&st.shift, &fresh.shift));
            assert!(bits(&st.s, &fresh.s));
            assert_eq!(st.at_ub, fresh.at_ub);
        }

        // A finiteness-category change (lb −∞ → finite) is declined: the
        // caller rebuilds. The stale cache must say no, not lie.
        let (old_lb, old_ub) = (lp.lower[0], lp.upper[0]);
        lp.lower[0] = f64::NEG_INFINITY;
        assert!(!refresh_dual_standard(&lp, &cfg_np(), &mut st));
        // ... and invalid bounds are declined too (rebuild reports them).
        lp.lower[0] = old_lb;
        lp.upper[0] = old_lb - 1.0; // ub < lb
        assert!(!refresh_dual_standard(&lp, &cfg_np(), &mut st));
        lp.upper[0] = old_ub;
        assert!(refresh_dual_standard(&lp, &cfg_np(), &mut st));
    }

    /// The reset must mirror the build's SLOT-wise at_ub/cost-shift
    /// `s` init for EVERY slot category. The absorbed WIP tested the caller
    /// `lp.upper[j] > 0` and skipped free splits, so flip slots, fixed
    /// columns and negative-cost free halves drifted from a fresh rebuild.
    #[test]
    fn p6_refresh_matches_rebuild_on_all_slot_categories() {
        let mut lp = Lp {
            num_col: 4,
            num_row: 2,
            a: cm(&[&[1.0, 1.0, 1.0, 1.0], &[1.0, 2.0, 1.0, 3.0]], 4),
            rhs: vec![10.0, 12.0],
            objective: vec![2.0, -3.0, -1.0, -1.0],
            maximize: false,
            senses: "LL".chars().collect(),
            lower: vec![f64::NEG_INFINITY, 2.0, f64::NEG_INFINITY, 0.0],
            upper: vec![5.0, 2.0, f64::INFINITY, 4.0],
        }; // j0 flip (lb −∞), j1 fixed, j2 free, j3 plain lb-finite
        let mut st = to_dual_standard(&lp, &cfg_np()).unwrap();
        let bits = |a: &[f64], b: &[f64]| a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits());

        // Bounds-value edits only, across every category: tighten the flip
        // ub, translate the fixed column (lb and ub), tighten a plain ub.
        lp.upper[0] = 4.0;
        lp.lower[1] = 3.0;
        lp.upper[1] = 3.0;
        lp.upper[3] = 2.5;
        assert!(refresh_dual_standard(&lp, &cfg_np(), &mut st));
        let fresh = to_dual_standard(&lp, &cfg_np()).unwrap();
        assert_eq!(
            st.at_ub, fresh.at_ub,
            "at_ub drifted (flip/fixed/free slots)"
        );
        assert!(
            bits(&st.s, &fresh.s),
            "s drifted (free/flip cost-shift init)"
        );
        assert!(bits(&st.ub, &fresh.ub));
        assert!(bits(&st.shift, &fresh.shift));
        let bmax = fresh.b.iter().fold(0.0f64, |m, &x| m.max(x.abs()));
        assert!(
            st.b.iter()
                .zip(&fresh.b)
                .all(|(x, y)| (x - y).abs() <= 1e-9 * bmax.max(1.0)),
            "b drifted after the fixed-column translation"
        );
    }

    /// End-to-end: one StdCache driven through the MIP node
    /// loop's whole edit grammar — value tighten AND relax (siblings),
    /// a category change (ub +inf → finite), a generation bump with an
    /// appended cut row — and every solve matches both the cold truth and
    /// the plain (rebuild-every-time) `reoptimize_dse`.
    #[test]
    fn p6_std_cache_node_loop_matches_truth() {
        let root_lp = Lp {
            num_col: 4,
            num_row: 3,
            a: cm(
                &[
                    &[1.0, 1.0, 2.0, 1.0],
                    &[4.0, 3.0, 1.0, 2.0],
                    &[2.0, 5.0, 3.0, 1.0],
                ],
                4,
            ),
            rhs: vec![12.0, 18.0, 15.0],
            objective: vec![-3.0, -5.0, -4.0, -2.0],
            maximize: false,
            senses: "LLL".chars().collect(),
            lower: vec![0.0; 4],
            upper: vec![6.0, f64::INFINITY, 6.0, 6.0],
        };
        let root = cold_solve(&root_lp);
        assert_eq!(root.status, SimplexStatus::Optimal);
        let basis = root.basis.clone().expect("root basis");

        let mut cache = StdCache::new();
        let mut lp = root_lp.clone();
        // Each step: (edit, expected Optimal?) — every node LP here is
        // feasible and bounded, so truth is always Optimal.
        type Edit<'a> = (&'a str, Box<dyn Fn(&mut Lp)>, u64);
        let edits: Vec<Edit<'_>> = vec![
            (
                "child A: x0 ≤ 0",
                Box::new(|l: &mut Lp| l.upper[0] = 0.0),
                0,
            ),
            (
                "sibling B: x0 ≥ 1 (relax vs cache)",
                Box::new(|l: &mut Lp| {
                    l.upper[0] = 6.0;
                    l.lower[0] = 1.0;
                }),
                0,
            ),
            (
                "probe: x1 ≤ 4.5 (category change +inf→finite)",
                Box::new(|l: &mut Lp| {
                    l.lower[0] = 0.0;
                    l.upper[1] = 4.5;
                }),
                0,
            ),
            (
                "probe revert: x1 back to +inf",
                Box::new(|l: &mut Lp| {
                    l.upper[1] = f64::INFINITY;
                }),
                0,
            ),
            ("branch: x3 ≥ 2", Box::new(|l: &mut Lp| l.lower[3] = 2.0), 0),
        ];
        for (name, edit, generation) in edits {
            edit(&mut lp);
            let truth = cold_solve(&lp);
            assert_eq!(truth.status, SimplexStatus::Optimal, "{name}");
            let cached = reoptimize_dse_cached(
                &lp,
                &cfg_np(),
                &basis,
                None,
                None,
                DseMode::Off,
                &mut cache,
                generation,
            )
            .unwrap();
            let plain = reoptimize_dse(&lp, &cfg_np(), &basis, None, None, DseMode::Off).unwrap();
            assert_eq!(cached.status, ReoptStatus::Optimal, "{name} cached");
            assert_eq!(plain.status, ReoptStatus::Optimal, "{name} plain");
            assert!(
                (cached.objective - truth.objective).abs() < 1e-9,
                "{name}: cached {} vs truth {}",
                cached.objective,
                truth.objective
            );
            assert!(
                (plain.objective - truth.objective).abs() < 1e-9,
                "{name}: plain {} vs truth {}",
                plain.objective,
                truth.objective
            );
        }

        // Generation bump: append a cut row (structure edit) and re-solve.
        // The cache must rebuild through it, not refresh across it.
        lp.num_row += 1;
        lp.rhs.push(20.0);
        lp.senses.push('L');
        lp.a = {
            let optr = lp.a.col_ptr().to_vec();
            let oridx = lp.a.row_idx().to_vec();
            let oval = lp.a.values().to_vec();
            let mut trip: Vec<(usize, usize, f64)> = Vec::with_capacity(oval.len() + 2);
            for j in 0..lp.num_col {
                for k in optr[j]..optr[j + 1] {
                    trip.push((oridx[k] as usize, j, oval[k]));
                }
            }
            trip.push((3, 0, 1.0));
            trip.push((3, 2, 1.0));
            rusolve_core::CscMatrix::from_triplets(4, lp.num_col, &trip)
        };
        let truth = cold_solve(&lp);
        assert_eq!(truth.status, SimplexStatus::Optimal, "cut LP");
        let cached = reoptimize_dse_cached(
            &lp,
            &cfg_np(),
            &basis,
            None,
            None,
            DseMode::Off,
            &mut cache,
            1,
        )
        .unwrap();
        assert_eq!(cached.status, ReoptStatus::Optimal, "cut LP cached");
        assert!(
            (cached.objective - truth.objective).abs() < 1e-9,
            "cut LP: cached {} vs truth {}",
            cached.objective,
            truth.objective
        );

        // And a bounds edit AFTER the bump refreshes against the new shape.
        lp.upper[2] = 3.0;
        let truth = cold_solve(&lp);
        assert_eq!(truth.status, SimplexStatus::Optimal, "post-cut branch");
        let cached = reoptimize_dse_cached(
            &lp,
            &cfg_np(),
            &basis,
            None,
            None,
            DseMode::Off,
            &mut cache,
            1,
        )
        .unwrap();
        assert_eq!(
            cached.status,
            ReoptStatus::Optimal,
            "post-cut branch cached"
        );
        assert!(
            (cached.objective - truth.objective).abs() < 1e-9,
            "post-cut branch: cached {} vs truth {}",
            cached.objective,
            truth.objective
        );
    }
}
