//! Forrest–Tomlin (FT) LU update.
//!
//! Provenance (reimplemented from these sources, not a line-by-line translation):
//! - Forrest & Tomlin (1972), *Updated triangular factors of the basis to
//!   maintain sparsity in the product form simplex method*, Math. Prog. 2.
//! - Suhl & Suhl (1993), *A fast LU update for linear programming*, Ann. OR 43.
//! - Koberstein (2005), PhD thesis, ch. 5.
//! - HiGHS `util/HFactor.cpp` `updateFT` :2300, `ftranFT` :1848,
//!   `btranFT` :1882 (MIT-licensed reference source).
//!
//! Representation: the base `SparseLu` (`P·B·Q = L·U`) is kept frozen on
//! the L side; the U side is re-housed in an *arena* of logical columns
//! that updates rewrite in place, mirrored row-wise by a UR back-link
//! arena (HiGHS `HFactor`'s UR, `updateFT`). Per basis change
//! replacing position `leaving_pos`:
//!
//! - the leaving position's logical column is vacated (pivot marker `-1`,
//!   column emptied) and its pivot *slot* (permuted row index) is reused;
//! - the entering column's FTRAN result — in slot space — is appended as
//!   the last logical column, pivoting at the vacated slot with value
//!   `old_pivot · alpha` (`alpha` = the tableau pivot);
//! - one row-eta `R` is appended to the file: entries `−ep_k · old_pivot`
//!   at the BTRAN-result (`ep`) slots, unit pivot at the vacated slot.
//!
//! With `B = Pᵀ L R₁···R_k Ũ_k Qᵀ` the solves compose as
//! FTRAN `x = Ũₖ⁻¹ Rₖ⁻¹···R₁⁻¹ L⁻¹ P b` (etas oldest-first, then the
//! arena back-sweep) and BTRAN `x = Pᵀ L⁻ᵀ R₁⁻ᵀ···Rₖ⁻ᵀ Ũₖ⁻ᵀ Q c` (arena
//! forward sweep, etas newest-first, then `L⁻ᵀ`), mirroring the HiGHS
//! pass order. Slot space (permuted rows) is the single internal index
//! space; positions enter/leave only through the delivery maps
//! (`logical_col_pos`) and the solve boundaries.
//!
//! Tested for correctness differentially against `pfi.rs` and dense ground
//! truth. Performance work added the bbar-recovered spike,
//! the PFI-calibrated refactor clock, hyper-sparse FTRAN, the UR row
//! arena, sparse update bookkeeping, and the hyper ep sweep. Interstitial
//! ftran/btran (engine-side solve reuse) remain future work.

use std::cell::RefCell;
use std::collections::BinaryHeap;

use crate::hvec::HVec;
use crate::lu::{HYPER_BACKOFF, HYPER_DIV, HYPER_MIN_N, SparseLu};
use crate::pfi::{ProductFormFactorization, RfReason, UpdateStatus, pfi_stats};
use crate::{FactorError, PfiConfig};

/// Synthetic-work charge per unit of R-file/vacated-logical solve
/// overhead — the same constant and model as PFI's refactor clock
/// (`pfi.rs` `ETA_TICK`, HiGHS HEkk.cpp `updateFactor`; see
/// `ProductFormFactorization::refactor_due` for the trigger policy).
const FT_ETA_TICK: f64 = 0.15;

/// Per-row slack in the UR arena at seed (HiGHS `ur_stuff_size`,
/// `HFactor.cpp` :1443 — 5 for the FT update method).
const UR_STUFF: u32 = 5;

/// Selection lever for the basis-update method.
///
/// `false` (default) keeps the product-form eta file (`pfi.rs`) — the
/// behavior every gate was measured against. `true` selects the FT update
/// (`RUSOLVE_FT=1`).
#[derive(Debug, Clone, Copy, Default)]
pub struct FtConfig {
    pub enabled: bool,
}

impl FtConfig {
    /// Read the lever from the environment: `RUSOLVE_FT=1`/`on` enables FT.
    /// Used by the A/B harness before the config is threaded through the
    /// solve options (slice-3 default decision).
    pub fn from_env() -> Self {
        Self {
            enabled: ft_enabled_from_env(),
        }
    }
}

/// True when `RUSOLVE_FT` is set to `1` or `on` (case-insensitive).
pub fn ft_enabled_from_env() -> bool {
    std::env::var("RUSOLVE_FT").is_ok_and(|v| parse_ft_flag(&v))
}

/// Parse a flag value: `1`/`on` (case-insensitive) enable, anything else
/// disables. Split out from env access so tests need no process env
/// (edition-2024 makes `set_var` unsafe in tests).
fn parse_ft_flag(value: &str) -> bool {
    matches!(value.trim().to_ascii_lowercase().as_str(), "1" | "on")
}

/// Forrest–Tomlin updated factorization: a frozen base `SparseLu` (L
/// side + rebuild anchor) plus the evolving U arena and R-eta file
/// described in the module docs. External interface mirrors `pfi.rs`
/// (`solve_b`/`solve_bt`/`update`) so the two are drop-in differentials.
#[derive(Debug, Clone)]
pub struct FtFactorization {
    n: usize,
    base: SparseLu,
    config: PfiConfig,
    /// Current basis columns by position (refactor input, PFI-parity).
    basis_nz: Vec<Vec<(usize, f64)>>,
    /// U arena, column form: entries of logical column `p` live at
    /// `u_index/u_value[u_start[p]..u_last[p]]` (slot indices).
    u_start: Vec<u32>,
    u_last: Vec<u32>,
    u_index: Vec<u32>,
    u_value: Vec<f64>,
    /// Logical → pivot slot (`-1` = vacated by an update).
    u_pivot_index: Vec<i32>,
    /// Slot → live logical (total: every slot keeps exactly one live
    /// logical; updates reassign the vacated slot to the appended column).
    u_pivot_lookup: Vec<u32>,
    /// Logical → pivot value (`old_pivot · alpha` after an update).
    u_pivot_value: Vec<f64>,
    /// Logical → basis position of the delivered FTRAN component (`-1` =
    /// vacated; live logicals deliver every position exactly once).
    logical_col_pos: Vec<i32>,
    /// Position → live logical (inverse of `logical_col_pos`).
    pos_logical: Vec<i32>,
    /// R-eta file: eta `i` has unit pivot at slot `pf_pivot_slot[i]` and
    /// entries `pf_index/pf_value[pf_start[i]..pf_start[i+1]]`.
    pf_pivot_slot: Vec<u32>,
    pf_start: Vec<u32>,
    pf_index: Vec<u32>,
    pf_value: Vec<f64>,
    /// UR row mirror of the arena (HiGHS `HFactor` UR, `updateFT` in
    /// `HFactor.cpp`): for pivot slot `r`,
    /// the live logicals holding a U entry at slot `r` — logical ids only
    /// (HiGHS also solves from UR row-wise and mirrors the values; ours is
    /// update-only). Makes the pivot-slot deletion sweep O(row nnz)
    /// instead of O(arena nnz). Segments share one buffer with per-row
    /// slack; an overflowing row relocates to the buffer end with
    /// `1.1× + 5` slack (HiGHS's amortized arena).
    ur_start: Vec<u32>,
    ur_lastp: Vec<u32>,
    ur_space: Vec<u32>,
    ur_index: Vec<u32>,
    /// Accumulated solve-work charge of the R file and the vacated
    /// logicals since the last refactorization (a `Cell`: solves take
    /// `&self`). Drives the synthetic-tick refactor clock — the same
    /// policy PFI's eta file is calibrated against (`pfi.rs`
    /// `refactor_due`; HiGHS HEkk.cpp `updateFactor`). Without it the
    /// hard `eta_limit` lets the R file grow ~20× longer than PFI's,
    /// taxing every solve and drifting the degenerate-instance
    /// trajectories (qap10/dano3_3 hit IterationLimit).
    eta_ticks: std::cell::Cell<f64>,
    /// Refactorizations since build.
    epoch: u64,
    /// Hyper-sparse FTRAN attempts left to skip after a reach blow-up
    /// (mirror of `SparseLu::hyper_skip`).
    hyper_skip: std::cell::Cell<u32>,
    /// Same backoff for the hyper ep arm of the update (a shared counter
    /// would let ep blow-ups starve FTRAN hyper attempts and vice versa).
    ep_hyper_skip: std::cell::Cell<u32>,
    /// Hyper-sparse FTRAN solves taken (test visibility: proves the arm
    /// ran rather than every entry gate declining to dense).
    hyper_took: std::cell::Cell<u64>,
    /// Persistent solve scratch.
    ws: RefCell<FtWorkspace>,
    /// Hyper-sparse FTRAN scratch (Gilbert–Peierls reach over the Û
    /// arena, logical-space marks + slot-space pattern; see
    /// [`Self::ftran_hyper`]).
    hyper: RefCell<FtHyper>,
    /// Refinement scratch, disjoint from `ws` (see [`FtRefineWorkspace`]).
    wsr: RefCell<FtRefineWorkspace>,
    /// Max |basis entry| — the conservative refinement-floor factor
    /// (overestimates are safe; PFI keeps the same invariant).
    col_max: f64,
}

#[derive(Debug, Clone, Default)]
struct FtWorkspace {
    /// Slot-space FTRAN vector.
    w: Vec<f64>,
    /// Slot-space BTRAN vector (the fused Û⁻ᵀ sweep lives here).
    v: Vec<f64>,
    /// Slot-space spike column (`u`, update scratch).
    spike: Vec<f64>,
    /// Slot stamps for `spat` dedup (per-update epoch in `stamp`).
    slot_mark: Vec<u32>,
    stamp: u32,
    /// Slots scattered into the spike by [`FtFactorization::update_with_bbar`]
    /// (unique; values may still cancel to zero).
    spat: Vec<u32>,
    /// Slot-space ep sweep vector — ALL ZERO at rest (both ep arms restore
    /// it, which is what lets the hyper arm's dense finish read exact zeros
    /// for unreached slots; the [`Self::btran_impl`] sweeps may leave this
    /// dirty, so it is never shared with them).
    vh: Vec<f64>,
    /// Nonzero slots of the ep sweep `z = Û⁻ᵀQᵀe_leaving` with their
    /// values, sorted ascending by slot — the R-eta source.
    zpat: Vec<(u32, f64)>,
}

impl FtFactorization {
    /// Build from sparse basis columns (`(row, value)` ascending by row),
    /// PFI-parity: seeds the arena from the base `SparseLu` factors.
    pub fn build_sparse(
        basis_nz: Vec<Vec<(usize, f64)>>,
        config: PfiConfig,
    ) -> Result<Self, FactorError> {
        let n = basis_nz.len();
        let base = SparseLu::factor_sparse(&basis_nz, config.pivot_tol)?;
        let col_max = basis_nz
            .iter()
            .flat_map(|c| c.iter())
            .map(|(_, v)| v.abs())
            .fold(0.0f64, f64::max);
        let mut f = Self {
            n,
            base,
            config,
            basis_nz,
            u_start: Vec::new(),
            u_last: Vec::new(),
            u_index: Vec::new(),
            u_value: Vec::new(),
            u_pivot_index: Vec::new(),
            u_pivot_lookup: Vec::new(),
            u_pivot_value: Vec::new(),
            logical_col_pos: Vec::new(),
            pos_logical: Vec::new(),
            pf_pivot_slot: Vec::new(),
            pf_start: vec![0],
            pf_index: Vec::new(),
            pf_value: Vec::new(),
            ur_start: Vec::new(),
            ur_lastp: Vec::new(),
            ur_space: Vec::new(),
            ur_index: Vec::new(),
            eta_ticks: std::cell::Cell::new(0.0),
            epoch: 0,
            hyper_skip: std::cell::Cell::new(0),
            ep_hyper_skip: std::cell::Cell::new(0),
            hyper_took: std::cell::Cell::new(0),
            ws: RefCell::new(FtWorkspace::default()),
            hyper: RefCell::new(FtHyper::default()),
            wsr: RefCell::new(FtRefineWorkspace::default()),
            col_max,
        };
        f.seed_arena();
        Ok(f)
    }

    pub fn n(&self) -> usize {
        self.n
    }

    /// Refactorizations since build (PFI-parity).
    pub fn factor_epoch(&self) -> u64 {
        self.epoch
    }

    /// Number of R-etas currently in the file.
    pub fn eta_count(&self) -> usize {
        self.pf_pivot_slot.len()
    }

    /// Profiling passthrough.
    pub fn lu_solve_work(&self) -> usize {
        self.base.solve_work()
    }

    /// Profiling passthrough.
    pub fn basis_nonzeros(&self) -> usize {
        self.basis_nz.iter().map(|c| c.len()).sum()
    }

    /// Re-seed the arena and clear the R file from the current basis
    /// columns (build + post-refactor state).
    fn seed_arena(&mut self) {
        let seed = self.base.ft_seed();
        self.u_start = seed.u_start;
        self.u_last = seed.u_last;
        self.u_index = seed.u_index;
        self.u_value = seed.u_value;
        self.u_pivot_value = seed.u_pivot_value;
        let n = self.n;
        // Identity maps at build: logical p = permuted column p, pivoting
        // at slot p; delivery follows the base column permutation Q.
        self.u_pivot_index = (0..n as i32).collect();
        self.u_pivot_lookup = (0..n as u32).collect();
        let col_perm = self.base.col_perm().to_vec();
        self.logical_col_pos = col_perm.iter().map(|&p| p as i32).collect();
        self.pos_logical = vec![-1; n];
        for (k, &p) in col_perm.iter().enumerate() {
            self.pos_logical[p] = k as i32;
        }
        self.pf_pivot_slot.clear();
        self.pf_index.clear();
        self.pf_value.clear();
        self.pf_start = vec![0];
        self.eta_ticks.set(0.0);
        self.seed_ur();
    }

    /// Rebuild the UR row mirror from the freshly seeded arena (HiGHS
    /// build, `HFactor.cpp` :1440-1470: exact per-row counts plus
    /// `ur_stuff_size = 5` slack).
    fn seed_ur(&mut self) {
        let n = self.n;
        self.ur_start.clear();
        self.ur_lastp.clear();
        self.ur_space.clear();
        self.ur_index.clear();
        if n == 0 {
            return;
        }
        let mut counts = vec![0u32; n];
        for p in 0..n {
            let (a, bnd) = (self.u_start[p] as usize, self.u_last[p] as usize);
            for k in a..bnd {
                counts[self.u_index[k] as usize] += 1;
            }
        }
        self.ur_start.reserve(n);
        self.ur_lastp.reserve(n);
        self.ur_space.reserve(n);
        let mut start = 0u32;
        for &c in &counts {
            self.ur_start.push(start);
            self.ur_lastp.push(start);
            self.ur_space.push(UR_STUFF);
            start += c + UR_STUFF;
        }
        self.ur_index.resize(start as usize, u32::MAX);
        for p in 0..n {
            let (a, bnd) = (self.u_start[p] as usize, self.u_last[p] as usize);
            for k in a..bnd {
                let r = self.u_index[k] as usize;
                let put = self.ur_lastp[r] as usize;
                self.ur_index[put] = p as u32;
                self.ur_lastp[r] = put as u32 + 1;
            }
        }
    }

    /// Make room in UR row `r` by relocating its segment to the buffer end
    /// with `1.1× + 5` slack (HiGHS `updateFT` relocation, verbatim policy).
    fn relocate_ur_row(&mut self, r: usize) {
        let start = self.ur_start[r] as usize;
        let count = (self.ur_lastp[r] - self.ur_start[r]) as usize;
        let new_start = self.ur_index.len();
        let new_space = count / 10 + 5;
        self.ur_index
            .resize(new_start + count + new_space, u32::MAX);
        self.ur_index.copy_within(start..start + count, new_start);
        self.ur_start[r] = new_start as u32;
        self.ur_lastp[r] = (new_start + count) as u32;
        self.ur_space[r] = new_space as u32;
    }

    /// Append logical `q` to UR row `r` (creating slack when the row is
    /// full). The caller guarantees `q` has a live arena entry at slot `r`.
    fn ur_push(&mut self, r: usize, q: u32) {
        if self.ur_space[r] == 0 {
            self.relocate_ur_row(r);
        }
        self.ur_space[r] -= 1;
        let put = self.ur_lastp[r] as usize;
        self.ur_index[put] = q;
        self.ur_lastp[r] = put as u32 + 1;
    }

    /// Remove logical `q` from UR row `r` (swap-with-last). The back-link
    /// must exist — HiGHS's UR bookkeeping is exact, and a missing link
    /// means the arena and mirror have diverged.
    fn ur_remove(&mut self, r: usize, q: u32) {
        let (mut i, mut last) = (self.ur_start[r] as usize, self.ur_lastp[r] as usize);
        while i < last && self.ur_index[i] != q {
            i += 1;
        }
        debug_assert!(i < last, "UR back-link missing: row {r}, logical {q}");
        if i < last {
            last -= 1;
            self.ur_index[i] = self.ur_index[last];
            self.ur_lastp[r] = last as u32;
            self.ur_space[r] += 1;
        }
    }

    /// Rebuild from the explicit basis columns; called when the synthetic
    /// tick clock trips (R-file/vacated-logical solve work ≥ rebuild
    /// work, PFI's calibrated policy) or at the hard `eta_limit` cap.
    #[allow(dead_code)] // API parity with PFI; callers use refactorize_reason
    pub(crate) fn refactorize(&mut self) -> Result<(), FactorError> {
        self.refactorize_reason(RfReason::Forced)
    }

    /// [`Self::refactorize`] recording WHY (cadence census; the FT
    /// path counts reasons + depth only — the census runs on the PFI path).
    pub(crate) fn refactorize_reason(&mut self, reason: RfReason) -> Result<(), FactorError> {
        if let Some(st) = pfi_stats() {
            st.rf_reason[reason.index()].set(st.rf_reason[reason.index()].get() + 1);
            let depth = self.pf_pivot_slot.len() as u64;
            st.rf_depth_sum.set(st.rf_depth_sum.get() + depth);
            st.rf_depth_max.set(st.rf_depth_max.get().max(depth));
        }
        // Transactional: build from a clone so a failed factorization
        // leaves `self` intact (`mem::take` + `?` gutted the object — the
        // cod105 panic: an empty `basis_nz` indexed by a later update).
        // Refactors are rare; the clone is noise against a full factor.
        let cols = self.basis_nz.clone();
        self.base = SparseLu::factor_sparse(&cols, self.config.pivot_tol)?;
        self.basis_nz = cols;
        self.seed_arena();
        self.epoch += 1;
        Ok(())
    }

    /// FTRAN `B x = b` (positions in, positions out), exposing the
    /// slot-space sweep vector as `aq` for [`Self::update`].
    ///
    /// `w` is kept ALL-ZERO at rest (the `SparseLu::hs.y` discipline —
    /// untouched entries are true zeros, which the hyper-sparse path's
    /// dense finish and entry subtractions rely on): both arms restore
    /// it before returning.
    fn ftran_slots(&self, b: &[f64], x: &mut [f64], w: &mut [f64]) {
        if self.ftran_hyper(b, x, w) {
            return;
        }
        self.ftran_pre_u(b, w);
        self.u_solve(w, x);
        w.fill(0.0);
    }

    /// Hyper-sparse FT FTRAN (`SparseLu::solve_hyper`'s Gilbert–Peierls
    /// structure, Hall & McKinnon 2005, ported onto the FT chain):
    /// reach-limited L sweep → R-eta pass with pattern growth →
    /// reach-limited Û sweep (heap over logicals, values in slot space).
    /// Bit-identical to the dense path (same op sequence per touched
    /// logical; untouched slots hold exact zeros by the `w` rest
    /// discipline). Returns `false` when the entry gate
    /// (size/density/backoff) rejects — the caller runs the dense path.
    fn ftran_hyper(&self, b: &[f64], x: &mut [f64], w: &mut [f64]) -> bool {
        let n = self.n;
        let skip = self.hyper_skip.get();
        if skip > 0 {
            self.hyper_skip.set(skip - 1);
            return false;
        }
        if n < HYPER_MIN_N {
            return false;
        }
        let mut hyb = self.hyper.borrow_mut();
        let hy = &mut *hyb;
        let (took, sparse_l) = self.base.ftran_l_hyper(b, w, &mut hy.lpat);
        if !took {
            debug_assert!(w.iter().all(|v| *v == 0.0));
            return false;
        }
        self.charge_solve_work();
        self.hyper_took.set(self.hyper_took.get() + 1);
        if !sparse_l {
            // Dense L finish: no pattern — dense etas + dense Û sweep.
            self.apply_r_etas_dense(w);
            self.u_solve(w, x);
            w.fill(0.0);
            return true;
        }
        // Slot-space pattern stamps; `seeds` holds the slot pattern (L
        // sweep output, then eta pivots that flip zero → nonzero).
        hy.sstamp = hy.sstamp.wrapping_add(1);
        if hy.sstamp == 0 {
            hy.smark.fill(0);
            hy.sstamp = 1;
        }
        if hy.smark.len() != n {
            hy.smark.clear();
            hy.smark.resize(n, 0);
        }
        if hy.seeds.capacity() < n {
            hy.seeds = Vec::with_capacity(n);
        }
        let stamp = hy.sstamp;
        hy.seeds.clear();
        for &k in &hy.lpat {
            hy.smark[k as usize] = stamp;
            hy.seeds.push(k);
        }
        // R-eta pass, identical arithmetic to the dense path (unmarked
        // slots read exact zeros — the `w` rest invariant), recording
        // pivot slots that flip zero → nonzero into the pattern.
        for i in 0..self.pf_pivot_slot.len() {
            let s = self.pf_pivot_slot[i] as usize;
            let mut v1 = w[s];
            let (a, bnd) = (self.pf_start[i] as usize, self.pf_start[i + 1] as usize);
            for k in a..bnd {
                v1 -= w[self.pf_index[k] as usize] * self.pf_value[k];
            }
            w[s] = v1;
            if v1 != 0.0 && hy.smark[s] != stamp {
                hy.smark[s] = stamp;
                hy.seeds.push(s as u32);
            }
        }
        // Û sweep over the live logicals owning the pattern slots. Max-
        // heap = descending logical order (the dense back-sweep order).
        // `touched` collects every marked logical so the solve can
        // restore the `w` rest invariant afterwards.
        hy.stamp = hy.stamp.wrapping_add(1);
        if hy.stamp == 0 {
            hy.mark.fill(0);
            hy.stamp = 1;
        }
        let nl = self.u_pivot_index.len();
        if hy.mark.len() != nl {
            hy.mark.clear();
            hy.mark.resize(nl, 0);
        }
        hy.heap.clear();
        hy.touched.clear();
        let lstamp = hy.stamp;
        let mark_logical = |hy: &mut FtHyper, q: u32| {
            if hy.mark[q as usize] != lstamp {
                hy.mark[q as usize] = lstamp;
                hy.heap.push(q);
                hy.touched.push(q);
            }
        };
        for si in 0..hy.seeds.len() {
            let q = self.u_pivot_lookup[hy.seeds[si] as usize];
            mark_logical(hy, q);
        }
        // Hyper delivery writes only reached positions.
        x.fill(0.0);
        let thr = n / HYPER_DIV;
        let mut dense = false;
        while let Some(p) = hy.heap.pop() {
            let p = p as usize;
            let s = self.u_pivot_index[p] as usize;
            let z = w[s] / self.u_pivot_value[p];
            w[s] = z;
            x[self.logical_col_pos[p] as usize] = z;
            if z != 0.0 {
                let (a, bnd) = (self.u_start[p] as usize, self.u_last[p] as usize);
                for k in a..bnd {
                    let i = self.u_index[k] as usize;
                    w[i] -= self.u_value[k] * z;
                    mark_logical(hy, self.u_pivot_lookup[i]);
                }
            }
            if hy.touched.len() > thr {
                // Dense finish from the next logical in sweep order
                // (descending; everything above p is already done). Safe
                // because unmarked slots hold true zeros (`w` rest).
                dense = true;
                for pd in (0..p).rev() {
                    let s = self.u_pivot_index[pd];
                    if s < 0 {
                        continue;
                    }
                    let s = s as usize;
                    let z = w[s] / self.u_pivot_value[pd];
                    w[s] = z;
                    x[self.logical_col_pos[pd] as usize] = z;
                    if z != 0.0 {
                        let (a, bnd) = (self.u_start[pd] as usize, self.u_last[pd] as usize);
                        for k in a..bnd {
                            w[self.u_index[k] as usize] -= self.u_value[k] * z;
                        }
                    }
                }
                break;
            }
        }
        if dense {
            self.hyper_skip.set(HYPER_BACKOFF);
            w.fill(0.0);
        } else {
            for &p in &hy.touched {
                w[self.u_pivot_index[p as usize] as usize] = 0.0;
            }
        }
        true
    }

    /// The R-eta application of the dense FTRAN path (oldest first,
    /// unconditional) — also the pre-Û half of [`Self::ftran_pre_u`].
    fn apply_r_etas_dense(&self, w: &mut [f64]) {
        for i in 0..self.pf_pivot_slot.len() {
            let s = self.pf_pivot_slot[i] as usize;
            let mut v1 = w[s];
            let (a, bnd) = (self.pf_start[i] as usize, self.pf_start[i + 1] as usize);
            for k in a..bnd {
                v1 -= w[self.pf_index[k] as usize] * self.pf_value[k];
            }
            w[s] = v1;
        }
    }

    /// FTRAN up to (and excluding) the Û sweep: `w ← R_k⁻¹···R_1⁻¹ L⁻¹ P b`
    /// in slot space. This `u` vector is what the FT spike column is made
    /// of (see [`Self::update`]).
    fn ftran_pre_u(&self, b: &[f64], w: &mut [f64]) {
        let n = self.n;
        debug_assert_eq!(b.len(), n, "ftran: wrong rhs length");
        debug_assert_eq!(w.len(), n, "ftran: wrong work length");
        self.charge_solve_work();
        // L⁻¹ P b — slot space.
        self.base.ftran_l(b, w);
        // R etas, oldest first (unit pivots, in-place row updates).
        for i in 0..self.pf_pivot_slot.len() {
            let s = self.pf_pivot_slot[i] as usize;
            let mut v1 = w[s];
            let (a, bnd) = (self.pf_start[i] as usize, self.pf_start[i + 1] as usize);
            for k in a..bnd {
                v1 -= w[self.pf_index[k] as usize] * self.pf_value[k];
            }
            w[s] = v1;
        }
    }

    /// Charge the R-file/vacated-logical solve overhead to the synthetic
    /// tick clock (PFI's `ETA_TICK` model; see [`Self::eta_ticks`]).
    /// Deliberately NOT charged: arena nnz growth — tried with the UR
    /// slice, it refactored too early and walked dano3_3 back into
    /// IterationLimit (dual A/B on 47 relaxed MIPLIB LPs). The R-file term alone
    /// is the cadence the degenerate instances are calibrated against.
    fn charge_solve_work(&self) {
        let overhead = (self.pf_index.len()
            + self.pf_pivot_slot.len()
            + self.u_pivot_index.len().saturating_sub(self.n)) as f64;
        self.eta_ticks
            .set(self.eta_ticks.get() + FT_ETA_TICK * overhead);
    }

    /// Û⁻¹ back-sweep over logical columns: `z_p` lands at the pivot slot
    /// (written back so `w` holds the delivered values at pivot slots) and
    /// is delivered to its position immediately (later-processed columns
    /// may still write the slot; the delivered value is final).
    fn u_solve(&self, w: &mut [f64], x: &mut [f64]) {
        for p in (0..self.u_pivot_index.len()).rev() {
            let s = self.u_pivot_index[p];
            if s < 0 {
                continue;
            }
            let s = s as usize;
            let z = w[s] / self.u_pivot_value[p];
            w[s] = z;
            x[self.logical_col_pos[p] as usize] = z;
            if z != 0.0 {
                let (a, bnd) = (self.u_start[p] as usize, self.u_last[p] as usize);
                for k in a..bnd {
                    w[self.u_index[k] as usize] -= self.u_value[k] * z;
                }
            }
        }
    }

    /// The update's ep row: `z = Û⁻ᵀQᵀe_leaving` — the pre-eta forward
    /// sweep ONLY. The R-eta identity consumes this mid-chain vector; the
    /// etas and `L⁻ᵀ` tail that [`Self::btran_body`] runs afterwards feed
    /// nothing here (the old update ran the full BTRAN and discarded the
    /// result). Fills `ws.zpat` with the sweep's `(slot, value)` nonzeros.
    ///
    /// Arms mirror [`Self::ftran_hyper`]: a reach-limited forward sweep
    /// driven through the UR back-links in Kahn topological order — the
    /// UR edge slot s (owner q) → column q′ means "`q′` reads my slot", so
    /// a column is evaluated only once every slot it reads is final
    /// (counted-down from its live-entry count; zero-contribution
    /// columns never become ready, which is exactly their value) — dense
    /// finish on reach blow-up, dense fallback on the entry gate. `ws.vh`
    /// is all zero at rest — that is what lets both the hyper sweep and
    /// its dense finish read exact zeros for unreached slots (the values
    /// the dense form computes for them). `z[s_q0] = 1/pivot_q0` at the
    /// seed logical, every other reached value bit-identical to the dense
    /// form.
    fn ep_sweep(&self, leaving_pos: usize, ws: &mut FtWorkspace) {
        self.charge_solve_work();
        let n = self.n;
        let nl = self.u_pivot_index.len();
        let q0 = self.pos_logical[leaving_pos] as usize;
        debug_assert!(q0 < nl && self.u_pivot_index[q0] >= 0);
        ws.zpat.clear();
        let FtWorkspace { vh, zpat, .. } = ws;
        let vh = vh.resize_to(n);
        let skip = self.ep_hyper_skip.get();
        let try_hyper = skip == 0 && n >= HYPER_MIN_N;
        if skip > 0 {
            self.ep_hyper_skip.set(skip - 1);
        }
        if !try_hyper {
            self.ep_dense(q0, vh, zpat);
            return;
        }
        let thr = n / HYPER_DIV;
        let mut hyb = self.hyper.borrow_mut();
        let hy = &mut *hyb;
        hy.stamp = hy.stamp.wrapping_add(1);
        if hy.stamp == 0 {
            hy.mark.fill(0);
            hy.stamp = 1;
        }
        if hy.mark.len() != nl {
            hy.mark.clear();
            hy.mark.resize(nl, 0);
        }
        let lstamp = hy.stamp;
        hy.eheap.clear();
        hy.touched.clear(); // slot-space nonzero list for the zero-rest restore
        hy.eheap.push(std::cmp::Reverse(q0 as u32));
        hy.mark[q0] = lstamp;
        let mut marked = 1usize;
        let mut dense = false;
        while let Some(std::cmp::Reverse(q)) = hy.eheap.pop() {
            let q = q as usize;
            let s = self.u_pivot_index[q] as usize;
            let mut acc = if q == q0 { 1.0 } else { 0.0 };
            let (a, bnd) = (self.u_start[q] as usize, self.u_last[q] as usize);
            for k in a..bnd {
                acc -= self.u_value[k] * vh[self.u_index[k] as usize];
            }
            let z = acc / self.u_pivot_value[q];
            vh[s] = z;
            if z != 0.0 {
                zpat.push((s as u32, z));
                hy.touched.push(s as u32);
                // Reached slot s taxes the columns holding entries at it.
                // Every entry points at a strictly lower logical (arena
                // invariant), so popping the min candidate finalizes all
                // of its contributor slots first.
                let (rs, re) = (self.ur_start[s] as usize, self.ur_lastp[s] as usize);
                for &q2 in &self.ur_index[rs..re] {
                    if hy.mark[q2 as usize] != lstamp {
                        hy.mark[q2 as usize] = lstamp;
                        hy.eheap.push(std::cmp::Reverse(q2));
                        marked += 1;
                    }
                }
            }
            if marked > thr {
                // Dense finish from the next logical in sweep order
                // (ascending; unreached slots below hold exact zeros).
                dense = true;
                for p in (q + 1)..nl {
                    let s = self.u_pivot_index[p];
                    if s < 0 {
                        continue;
                    }
                    let mut acc = if p == q0 { 1.0 } else { 0.0 };
                    let (a, bnd) = (self.u_start[p] as usize, self.u_last[p] as usize);
                    for k in a..bnd {
                        acc -= self.u_value[k] * vh[self.u_index[k] as usize];
                    }
                    let z = acc / self.u_pivot_value[p];
                    vh[s as usize] = z;
                    if z != 0.0 {
                        zpat.push((s as u32, z));
                    }
                }
                break;
            }
        }
        if dense {
            self.ep_hyper_skip.set(HYPER_BACKOFF);
            vh.fill(0.0);
        } else {
            for &s in &hy.touched {
                vh[s as usize] = 0.0;
            }
        }
        zpat.sort_unstable_by_key(|&(s, _)| s);
    }

    /// Dense ep form: full ascending live-logical sweep (the arithmetic of
    /// `btran_body`'s forward half with `b = e_leaving`), then restore the
    /// `vh` rest invariant.
    fn ep_dense(&self, q0: usize, vh: &mut [f64], zpat: &mut Vec<(u32, f64)>) {
        let nl = self.u_pivot_index.len();
        for p in 0..nl {
            let s = self.u_pivot_index[p];
            if s < 0 {
                continue;
            }
            let mut acc = if p == q0 { 1.0 } else { 0.0 };
            let (a, bnd) = (self.u_start[p] as usize, self.u_last[p] as usize);
            for k in a..bnd {
                acc -= self.u_value[k] * vh[self.u_index[k] as usize];
            }
            let z = acc / self.u_pivot_value[p];
            vh[s as usize] = z;
            if z != 0.0 {
                zpat.push((s as u32, z));
            }
        }
        vh.fill(0.0);
        zpat.sort_unstable_by_key(|&(s, _)| s);
    }

    /// Solves `B x = b` for the current basis (no iterative refinement —
    /// the PFI `refine` pass is deliberately out of scope here).
    ///
    /// # Panics
    ///
    /// Panics if `b.len()` or `x.len()` differ from `n`.
    pub fn solve_b(&self, b: &[f64], x: &mut [f64]) {
        assert_eq!(b.len(), self.n, "solve_b: wrong rhs length");
        assert_eq!(x.len(), self.n, "solve_b: wrong solution length");
        {
            let mut ws = self.ws.borrow_mut();
            let w = ws.w.resize_to(self.n);
            self.ftran_slots(b, x, w);
        }
        self.refine(b, x, false);
    }

    /// Unrefined forward solve — the correction kernel and the A/B lever's
    /// raw path (`BasisFactorization::solve_b_raw` maps here).
    ///
    /// # Panics
    ///
    /// Panics if `b.len()` or `x.len()` differ from `n`.
    pub fn solve_b_raw(&self, b: &[f64], x: &mut [f64]) {
        assert_eq!(b.len(), self.n, "solve_b: wrong rhs length");
        assert_eq!(x.len(), self.n, "solve_b: wrong solution length");
        let mut ws = self.ws.borrow_mut();
        let w = ws.w.resize_to(self.n);
        self.ftran_slots(b, x, w);
    }

    /// Solves `Bᵀ x = b` for the current basis, refined (PFI parity — the
    /// degenerate-instance trajectory drift without it shows up as
    /// IterationLimit on dano3_3/qap10 in the A/B on 47 relaxed MIPLIB LPs).
    ///
    /// # Panics
    ///
    /// Panics if `b.len()` or `x.len()` differ from `n`.
    pub fn solve_bt(&self, b: &[f64], x: &mut [f64]) {
        assert_eq!(b.len(), self.n, "solve_bt: wrong rhs length");
        assert_eq!(x.len(), self.n, "solve_bt: wrong solution length");
        self.btran_impl(b, x);
        self.refine(b, x, true);
    }

    /// One or two double-word residual-correction rounds, mirroring
    /// [`ProductFormFactorization::refine`] (same fp sequence in the dense
    /// loops, same floor constants — see that method's provenance: Dekker
    /// two-sum, Maddick–Ghysels–Hogg residual-accurate refinement). FT
    /// keeps no row index, so both residual sweeps use the column-wise
    /// dense paths only.
    fn refine(&self, rhs: &[f64], x: &mut [f64], transpose: bool) {
        let m = self.n;
        let mut wsr = self.wsr.borrow_mut();
        let FtRefineWorkspace { r, d, hi, lo } = &mut *wsr;
        let r = r.resize_to(m);
        let d = d.resize_to(m);
        let hi = hi.resize_to(m);
        let lo = lo.resize_to(m);
        for round in 0..2 {
            if transpose {
                self.residual_bwd(rhs, x, r, hi, lo);
            } else {
                self.residual_fwd(rhs, x, r, hi, lo);
            }
            if round > 0 {
                let rnorm = r.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                let xnorm = x.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                let rhs_norm = rhs.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
                let floor = (rhs_norm + xnorm * self.col_max * (m as f64)) * 1e-13;
                if rnorm <= floor {
                    return;
                }
            }
            // Correction solve into `d`; the raw kernels borrow `ws`, which
            // is disjoint from `wsr`.
            if transpose {
                self.btran_impl(r, d);
            } else {
                let mut ws = self.ws.borrow_mut();
                let w = ws.w.resize_to(m);
                self.ftran_slots(r, d, w);
            }
            for (xi, &di) in x.iter_mut().zip(d.iter()) {
                *xi += di;
            }
            let xn = x.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
            let dn = d.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
            if dn <= xn * 1e-8 {
                return;
            }
        }
    }

    /// `r ← rhs − B·x` in double-word arithmetic (column-wise dense path —
    /// bit-identical term order to PFI's non-row-indexed fallback).
    fn residual_fwd(&self, rhs: &[f64], x: &[f64], r: &mut [f64], hi: &mut [f64], lo: &mut [f64]) {
        hi.copy_from_slice(rhs);
        lo.fill(0.0);
        for (j, col) in self.basis_nz.iter().enumerate() {
            let xj = x[j];
            if xj == 0.0 {
                continue;
            }
            for &(i, a) in col {
                let (si, ci) = (&mut hi[i], &mut lo[i]);
                let pd = a * xj;
                let pe = a.mul_add(-xj, pd);
                let (h, l) = crate::pfi::two_sum(*si, -pd);
                *si = h;
                *ci += l + pe;
            }
        }
        for (ri, (h, l)) in r.iter_mut().zip(hi.iter().zip(lo.iter())) {
            *ri = *h + *l;
        }
    }

    /// `r ← rhs − Bᵀ·x` in double-word arithmetic (column-wise dense path).
    fn residual_bwd(&self, rhs: &[f64], x: &[f64], r: &mut [f64], _s: &mut [f64], _c: &mut [f64]) {
        for (j, col) in self.basis_nz.iter().enumerate() {
            let mut sum = 0.0;
            let mut c = 0.0;
            for &(i, a) in col {
                let xi = x[i];
                if xi == 0.0 {
                    continue;
                }
                let pd = a * xi;
                let pe = a.mul_add(-xi, pd);
                let (h, l) = crate::pfi::two_sum(sum, pd);
                sum = h;
                c += l - pe;
            }
            let (h, l) = crate::pfi::two_sum(rhs[j], -sum);
            r[j] = h + (l - c);
        }
    }

    /// BTRAN body (the engine's `solve_bt` path). The update does not use
    /// this: it needs only the pre-eta Û⁻ᵀ sweep — see
    /// [`Self::ep_sweep`].
    fn btran_impl(&self, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        let mut ws = self.ws.borrow_mut();
        let v = ws.v.resize_to(n);
        self.btran_body(b, x, v);
    }

    /// [`btran_impl`](Self::btran_impl) with caller-owned scratch.
    ///
    /// The Û⁻ᵀ sweep runs directly in slot space: every entry of column
    /// `p` points at a strictly lower logical `q` (the arena invariant —
    /// spike columns append last and the deletion sweep clears vacated
    /// slots), so `v[s_q]` is final when `p` reads it. That fuses the old
    /// position gather, logical-space `y` vector, and slot scatter into
    /// one loop and drops the per-entry `u_pivot_lookup` indirection —
    /// same fp operation sequence, one load per entry instead of two.
    fn btran_body(&self, b: &[f64], x: &mut [f64], v: &mut [f64]) {
        let nl = self.u_pivot_index.len();
        self.charge_solve_work();
        // Û⁻ᵀ forward sweep, slot space: v[s_p] = (b[pos_p] −
        // Σ_{(s,v) ∈ col p} v·v[s]) / pivot_p. Every slot is written by
        // exactly one live logical, so no clearing pass is needed.
        for p in 0..nl {
            let s = self.u_pivot_index[p];
            if s < 0 {
                continue;
            }
            let mut acc = b[self.logical_col_pos[p] as usize];
            let (a, bnd) = (self.u_start[p] as usize, self.u_last[p] as usize);
            for k in a..bnd {
                acc -= self.u_value[k] * v[self.u_index[k] as usize];
            }
            v[s as usize] = acc / self.u_pivot_value[p];
        }
        // R etas, newest first (column updates).
        for i in (0..self.pf_pivot_slot.len()).rev() {
            let s = self.pf_pivot_slot[i] as usize;
            let m = v[s];
            if m != 0.0 {
                let (a, bnd) = (self.pf_start[i] as usize, self.pf_start[i + 1] as usize);
                for k in a..bnd {
                    v[self.pf_index[k] as usize] -= self.pf_value[k] * m;
                }
            }
        }
        self.base.btran_lt(v, x);
    }

    /// FT update: the basis column at `leaving_pos` is replaced by
    /// `entering_nz` (`(row, value)` ascending, exact zeros omitted).
    /// The spike is computed by an internal pre-Û FTRAN of the entering
    /// column — callers that already FTRAN'd the entering column (the dual
    /// engine) should prefer [`Self::update_with_bbar`], which recovers the
    /// spike from the caller's solve result instead. On numerical trouble
    /// returns [`FactorError::RefactorizeNeeded`].
    pub fn update(
        &mut self,
        leaving_pos: usize,
        entering_nz: Vec<(usize, f64)>,
    ) -> Result<UpdateStatus, FactorError> {
        assert!(leaving_pos < self.n, "update: position out of range");
        // u = the entering column FTRAN'd through L and the R etas only
        // (pre-Û sweep) — the spike column source. The workspace is taken
        // out of its `RefCell` for the whole update (`update` holds
        // `&mut self`, so `get_mut` is free; `mem::take` sidesteps the
        // `&self`-method-vs-scratch borrow overlap without copying).
        let n = self.n;
        let mut ws = std::mem::take(self.ws.get_mut());
        {
            let FtWorkspace { v, w, spike, .. } = &mut ws;
            // `v` doubles as the dense rhs here: the update owns `ws`
            // exclusively, and the BTRAN solve restates every slot it uses.
            let b = v.resize_to(n);
            b.fill(0.0);
            for &(i, v) in &entering_nz {
                b[i] = v;
            }
            let w = w.resize_to(n);
            self.ftran_pre_u(b, w);
            let sp = spike.resize_to(n);
            sp.copy_from_slice(w);
            // Restore the all-zero-at-rest invariant of `w` (the hyper
            // FTRAN path depends on it).
            w.fill(0.0);
        }
        // No spike pattern: the append falls back to the dense scan of `u`.
        ws.spat.clear();
        self.update_from_spike(leaving_pos, entering_nz, &mut ws, false)
    }

    /// FT update reusing the caller's FTRAN: `bbar` is the engine's solve
    /// of `bbar_sign·a_q` (`bbar = B⁻¹(bbar_sign·a_q)` in position space,
    /// refined or raw — both fine). The slot-space spike `u` is recovered
    /// by inverting the Û delivery sweep: `u[s_p] = pivot_p·z_p +
    /// Σ_{columns with an entry at s_p} v·z_q` with `z = bbar_sign·bbar`
    /// (entries of a column point at lower logicals, so summing every
    /// column's `v·z_q` contribution into the entry slots, in any order,
    /// plus the pivot terms, reconstructs the pre-Û vector exactly — one
    /// arena pass, no L/R re-solve, and the spike is consistent with the
    /// very `bbar` the engine already consumed). Columns with `z_p = 0`
    /// contribute exact zeros and are skipped; the touched slots are
    /// stamped into `ws.spat` so the spike append never scans `u`.
    pub fn update_with_bbar(
        &mut self,
        leaving_pos: usize,
        entering_nz: Vec<(usize, f64)>,
        bbar: &[f64],
        bbar_sign: f64,
    ) -> Result<UpdateStatus, FactorError> {
        assert!(leaving_pos < self.n, "update: position out of range");
        assert_eq!(bbar.len(), self.n, "update_with_bbar: wrong bbar length");
        let n = self.n;
        let mut ws = std::mem::take(self.ws.get_mut());
        {
            let FtWorkspace {
                spike: u,
                slot_mark,
                stamp,
                spat,
                ..
            } = &mut ws;
            let u = u.resize_to(n);
            u.fill(0.0);
            if slot_mark.len() != n {
                slot_mark.clear();
                slot_mark.resize(n, 0);
            }
            *stamp = stamp.wrapping_add(1);
            if *stamp == 0 {
                slot_mark.fill(0);
                *stamp = 1;
            }
            let st = *stamp;
            spat.clear();
            // Entry scatter (ascending logical order — the fp operation
            // sequence of the dense two-pass form, zero columns skipped).
            for p in 0..self.u_pivot_index.len() {
                if self.u_pivot_index[p] < 0 {
                    continue;
                }
                let z_p = bbar_sign * bbar[self.logical_col_pos[p] as usize];
                if z_p == 0.0 {
                    continue;
                }
                let (a, bnd) = (self.u_start[p] as usize, self.u_last[p] as usize);
                for k in a..bnd {
                    let i = self.u_index[k] as usize;
                    u[i] += self.u_value[k] * z_p;
                    if slot_mark[i] != st {
                        slot_mark[i] = st;
                        spat.push(i as u32);
                    }
                }
            }
            // Pivot terms (second pass, same order as the dense form).
            for p in 0..self.u_pivot_index.len() {
                let s = self.u_pivot_index[p];
                if s < 0 {
                    continue;
                }
                let z_p = bbar_sign * bbar[self.logical_col_pos[p] as usize];
                if z_p == 0.0 {
                    continue;
                }
                let i = s as usize;
                u[i] += self.u_pivot_value[p] * bbar_sign * bbar[self.logical_col_pos[p] as usize];
                if slot_mark[i] != st {
                    slot_mark[i] = st;
                    spat.push(i as u32);
                }
            }
            // Ascending slot order: the spike append must lay its entries
            // exactly as the dense-scan form did (arena entry order fixes
            // the fp order of every later sweep).
            spat.sort_unstable();
        }
        self.update_from_spike(leaving_pos, entering_nz, &mut ws, true)
    }

    /// FT bookkeeping from a computed spike (`ws.spike`, slot space):
    /// BTRAN the R-eta reference `z = Û⁻ᵀQᵀe_leaving`, vacate the leaving
    /// logical, delete the pivot slot's entries from the arena columns,
    /// append the spike as the last logical and the R-eta to the file.
    /// `ws` is the already-taken workspace; it is restored on every path.
    /// `sparse_spike` says `ws.spat` holds the spike's touched slots
    /// (ascending, from [`Self::update_with_bbar`]); otherwise the append
    /// scans `ws.spike` densely.
    fn update_from_spike(
        &mut self,
        leaving_pos: usize,
        entering_nz: Vec<(usize, f64)>,
        ws: &mut FtWorkspace,
        sparse_spike: bool,
    ) -> Result<UpdateStatus, FactorError> {
        // ep row z = Û⁻ᵀQᵀe_leaving: the forward sweep only (hyper-first),
        // nonzeros in `ws.zpat` ascending with their values. The old form
        // ran a full BTRAN per pivot — etas, L⁻ᵀ and all — and threw the
        // position-space result away.
        self.ep_sweep(leaving_pos, ws);
        let outcome = self.update_bookkeep(leaving_pos, entering_nz, ws, sparse_spike);
        // Restore the workspace before any early return keeps it taken.
        *self.ws.get_mut() = std::mem::take(ws);
        outcome
    }

    /// Pure-`&mut self` FT bookkeeping (everything after the ep sweep):
    /// `ws.spike` holds the spike `u` in slot space, `ws.zpat` the ep
    /// nonzeros `(slot, z)` ascending.
    fn update_bookkeep(
        &mut self,
        leaving_pos: usize,
        entering_nz: Vec<(usize, f64)>,
        ws: &mut FtWorkspace,
        sparse_spike: bool,
    ) -> Result<UpdateStatus, FactorError> {
        let n = self.n;
        let FtWorkspace {
            spike, zpat, spat, ..
        } = ws;
        let u: &[f64] = spike.resize_to(n);
        // FT bookkeeping (HiGHS updateFT mirror; see module docs).
        let old_logical = self.pos_logical[leaving_pos];
        debug_assert!(
            old_logical >= 0,
            "update: leaving position has no live logical"
        );
        let old_logical = old_logical as usize;
        let s_p = self.u_pivot_index[old_logical] as usize;
        let pivot = self.u_pivot_value[old_logical];
        // R-eta entries: −z·pivot at the z-nonzero slots (unit pivot at
        // s_p) — straight off the ep sweep's nonzero list.
        let mut pf_idx: Vec<u32> = Vec::with_capacity(zpat.len());
        let mut pf_val: Vec<f64> = Vec::with_capacity(zpat.len());
        for &(a, za) in &**zpat {
            let a = a as usize;
            if a != s_p {
                pf_idx.push(a as u32);
                pf_val.push(-za * pivot);
            }
        }
        // New spike pivot: u[s_p] − Σ pf·u (≡ pivot·alpha — the tableau
        // pivot scaled by the vacated pivot, HiGHS's stored value).
        let mut new_pivot = u[s_p];
        for (&k, &pv) in pf_idx.iter().zip(&pf_val) {
            new_pivot -= pv * u[k as usize];
        }
        if new_pivot.abs() < self.config.pivot_tol {
            return Err(FactorError::RefactorizeNeeded { pivot: new_pivot });
        }
        // Vacate the old logical; its slot is inherited by the spike.
        self.u_pivot_index[old_logical] = -1;
        self.logical_col_pos[old_logical] = -1;
        self.pos_logical[leaving_pos] = -1;
        // Delete the pivot slot's entries from the columns that hold them
        // — over the UR row (HiGHS updateFT "delete pivotal row from U"),
        // O(row nnz) instead of the old O(arena nnz) full sweep. The row
        // is consumed whole: after this no live column references s_p.
        let (rs, re) = (self.ur_start[s_p] as usize, self.ur_lastp[s_p] as usize);
        for &q in &self.ur_index[rs..re] {
            let q = q as usize;
            let mut k = self.u_start[q] as usize;
            let mut last = self.u_last[q] as usize;
            let mut found = false;
            while k < last {
                if self.u_index[k] as usize == s_p {
                    last -= 1;
                    self.u_index[k] = self.u_index[last];
                    self.u_value[k] = self.u_value[last];
                    found = true;
                    break;
                }
                k += 1;
            }
            debug_assert!(
                found,
                "UR listed logical {q} without an entry at slot {s_p}"
            );
            if found {
                self.u_last[q] = last as u32;
            }
        }
        self.ur_lastp[s_p] = rs as u32;
        self.ur_space[s_p] += (re - rs) as u32;
        // Vacate the old logical's column: clear it and drop its UR
        // back-links (HiGHS updateFT "delete pivotal column from UR").
        {
            let mut k = self.u_start[old_logical] as usize;
            let last = self.u_last[old_logical] as usize;
            while k < last {
                let a = self.u_index[k] as usize;
                self.ur_remove(a, old_logical as u32);
                k += 1;
            }
        }
        self.u_last[old_logical] = self.u_start[old_logical];
        // Append the spike (u minus its pivot entry) as the last logical —
        // over the stamped pattern when the spike came from `bbar`, else
        // the dense scan. UR back-links are created alongside (HiGHS
        // updateFT "store column as UR elements").
        let start = self.u_index.len() as u32;
        debug_assert_eq!(self.u_start.last().copied(), Some(start));
        let new_p = self.u_pivot_index.len() as u32;
        if sparse_spike {
            for &a in &*spat {
                let a = a as usize;
                if a != s_p && u[a] != 0.0 {
                    self.u_index.push(a as u32);
                    self.u_value.push(u[a]);
                    self.ur_push(a, new_p);
                }
            }
        } else {
            for (a, &ua) in u.iter().enumerate() {
                if ua != 0.0 && a != s_p {
                    self.u_index.push(a as u32);
                    self.u_value.push(ua);
                    self.ur_push(a, new_p);
                }
            }
        }
        let end = self.u_index.len() as u32;
        self.u_last.push(end);
        self.u_start.push(end);
        self.u_pivot_index.push(s_p as i32);
        self.u_pivot_lookup[s_p] = new_p;
        self.u_pivot_value.push(new_pivot);
        self.logical_col_pos.push(leaving_pos as i32);
        self.pos_logical[leaving_pos] = new_p as i32;
        // R file: entries then pivot slot + segment end.
        self.pf_index.extend_from_slice(&pf_idx);
        self.pf_value.extend_from_slice(&pf_val);
        self.pf_pivot_slot.push(s_p as u32);
        self.pf_start.push(self.pf_index.len() as u32);
        // Basis bookkeeping (refactor input) + conservative col_max growth
        // (PFI keeps the same invariant; overestimates only loosen the
        // refinement floor).
        self.basis_nz[leaving_pos] = entering_nz;
        for &(_, v) in &self.basis_nz[leaving_pos] {
            let a = v.abs();
            if a > self.col_max {
                self.col_max = a;
            }
        }
        // Refactor cadence: the synthetic-tick clock first (PFI's
        // calibrated policy — the R file and the vacated logicals have
        // cost more solve work than a rebuild costs), with the hard
        // `eta_limit` cap retained as the backstop. The clock-less
        // cap-only cadence grew R files ~20× PFI's and cost 2–3.8×
        // wall on 47 relaxed MIPLIB LPs (dual engine, measured).
        let (k, m) = crate::pfi::rfcad().unwrap_or((self.config.clock_scale, 0));
        let due = if self.pf_pivot_slot.len() >= self.config.eta_limit {
            Some(RfReason::Cap)
        } else if self.pf_pivot_slot.len() >= self.config.min_updates.max(m)
            && self.eta_ticks.get() >= self.base.build_work() as f64 * k
        {
            Some(RfReason::Clock)
        } else {
            None
        };
        if let Some(reason) = due {
            self.refactorize_reason(reason)?;
            return Ok(UpdateStatus::Refactorized);
        }
        Ok(UpdateStatus::Updated)
    }
}

/// Engine-facing basis factorization: the product-form eta file (default)
/// or this module's FT update, selected by `RUSOLVE_FT=1` at construction.
/// The dual engine (`dual.rs`) is wired through this adapter so the two
/// update methods are drop-in swappable in the A/B harness; the primal
/// engine keeps `ProductFormFactorization` directly (primal wiring is not
/// done yet).
#[derive(Debug, Clone)]
pub enum BasisFactorization {
    Pfi(Box<ProductFormFactorization>),
    Ft(Box<FtFactorization>),
}

impl BasisFactorization {
    /// Construct per the `RUSOLVE_FT` lever (`1`/`on` → FT).
    pub fn build_sparse(
        basis_nz: Vec<Vec<(usize, f64)>>,
        config: PfiConfig,
    ) -> Result<Self, FactorError> {
        if ft_enabled_from_env() {
            FtFactorization::build_sparse(basis_nz, config).map(|f| Self::Ft(Box::new(f)))
        } else {
            ProductFormFactorization::build_sparse(basis_nz, config).map(|p| Self::Pfi(Box::new(p)))
        }
    }

    pub fn n(&self) -> usize {
        match self {
            Self::Pfi(p) => p.n(),
            Self::Ft(f) => f.n(),
        }
    }

    /// PFI applies iterative refinement; FT has none yet (`TODO(perf)`).
    pub fn solve_b(&self, b: &[f64], x: &mut [f64]) {
        match self {
            Self::Pfi(p) => p.solve_b(b, x),
            Self::Ft(f) => f.solve_b(b, x),
        }
    }

    /// No-refinement variant: for FT this is just [`Self::solve_b`].
    pub fn solve_b_raw(&self, b: &[f64], x: &mut [f64]) {
        match self {
            Self::Pfi(p) => p.solve_b_raw(b, x),
            Self::Ft(f) => f.solve_b(b, x),
        }
    }

    pub fn solve_bt(&self, b: &[f64], x: &mut [f64]) {
        match self {
            Self::Pfi(p) => p.solve_bt(b, x),
            Self::Ft(f) => f.solve_bt(b, x),
        }
    }

    /// No-refinement BTRAN (accuracy regime without refinement): for FT this is just
    /// [`Self::solve_bt`] (FT carries no refinement).
    pub(crate) fn solve_bt_raw(&self, b: &[f64], x: &mut [f64]) {
        match self {
            Self::Pfi(p) => p.solve_bt_raw(b, x),
            Self::Ft(f) => f.solve_bt(b, x),
        }
    }

    /// Dense-column update (the dual engine's call shape). PFI builds its
    /// eta from the caller's `bbar`; FT recovers the slot-space
    /// spike from the same `bbar` — neither side re-solves the entering
    /// column the caller just FTRAN'd.
    pub fn update(
        &mut self,
        leaving_pos: usize,
        entering_col: &[f64],
        bbar: &[f64],
        bbar_sign: f64,
    ) -> Result<UpdateStatus, FactorError> {
        match self {
            Self::Pfi(p) => p.update(leaving_pos, entering_col, bbar, bbar_sign),
            Self::Ft(f) => {
                let entering_nz: Vec<(usize, f64)> = entering_col
                    .iter()
                    .enumerate()
                    .filter(|&(_, &v)| v != 0.0)
                    .map(|(i, &v)| (i, v))
                    .collect();
                f.update_with_bbar(leaving_pos, entering_nz, bbar, bbar_sign)
            }
        }
    }

    /// Index-carrying FTRAN (`refine` selects `solve_b` vs
    /// `solve_b_raw`). FT has no HVec path yet: it runs its slice solve and
    /// returns `x` in dense mode (`TODO(perf)` once FT is a default).
    pub(crate) fn solve_b_hv(&self, b: &HVec, x: &mut HVec, refine: bool) {
        match self {
            Self::Pfi(p) => p.solve_b_hv(b, x, refine),
            Self::Ft(f) => f.solve_b(b.vals(), x.dense_mut()),
        }
    }

    /// Index-carrying BTRAN, see [`Self::solve_b_hv`].
    pub(crate) fn solve_bt_hv(&self, b: &HVec, x: &mut HVec, refine: bool) {
        match self {
            Self::Pfi(p) => p.solve_bt_hv(b, x, refine),
            Self::Ft(f) => f.solve_bt(b.vals(), x.dense_mut()),
        }
    }

    /// [`Self::update`] with index-carrying entering column / `bbar`
    /// (bit-identical eta and basis column for PFI; FT gets the same
    /// sorted nonzero list its dense adapter builds).
    pub(crate) fn update_hv(
        &mut self,
        leaving_pos: usize,
        entering: &mut HVec,
        bbar: &mut HVec,
        bbar_sign: f64,
    ) -> Result<UpdateStatus, FactorError> {
        match self {
            Self::Pfi(p) => p.update_hv(leaving_pos, entering, bbar, bbar_sign),
            Self::Ft(f) => {
                let entering_nz = entering.nonzeros_sorted();
                f.update_with_bbar(leaving_pos, entering_nz, bbar.vals(), bbar_sign)
            }
        }
    }

    pub fn factor_epoch(&self) -> u64 {
        match self {
            Self::Pfi(p) => p.factor_epoch(),
            Self::Ft(f) => f.factor_epoch(),
        }
    }

    pub fn eta_count(&self) -> usize {
        match self {
            Self::Pfi(p) => p.eta_count(),
            Self::Ft(f) => f.eta_count(),
        }
    }

    /// PFI tracks refactor time; FT does not yet (`TODO(perf)` — 0.0 keeps
    /// the engine's log lines honest about that).
    pub fn refactor_nanos(&self) -> f64 {
        match self {
            Self::Pfi(p) => p.refactor_nanos(),
            Self::Ft(_) => 0.0,
        }
    }

    /// Estimated LIVE eta-file bytes (sweep-cost census; see
    /// `ProductFormFactorization::eta_live_bytes`). The FT engine carries
    /// no eta file — 0 keeps the census honest about that.
    pub(crate) fn eta_live_bytes(&self) -> usize {
        match self {
            Self::Pfi(p) => p.eta_live_bytes(),
            Self::Ft(_) => 0,
        }
    }

    pub fn refactorize(&mut self) -> Result<(), FactorError> {
        self.refactorize_reason(RfReason::Forced)
    }

    /// [`Self::refactorize`] with the cadence-census reason tag.
    pub(crate) fn refactorize_reason(&mut self, reason: RfReason) -> Result<(), FactorError> {
        match self {
            Self::Pfi(p) => p.refactorize_reason(reason),
            Self::Ft(f) => f.refactorize_reason(reason),
        }
    }

    /// Profiling: L nnz + U nnz + 3n of the base LU (0 = dense engine).
    pub fn lu_solve_work(&self) -> usize {
        match self {
            Self::Pfi(p) => p.lu_solve_work(),
            Self::Ft(f) => f.lu_solve_work(),
        }
    }

    /// Profiling: total nonzeros of the explicit basis columns.
    pub fn basis_nonzeros(&self) -> usize {
        match self {
            Self::Pfi(p) => p.basis_nonzeros(),
            Self::Ft(f) => f.basis_nonzeros(),
        }
    }
}

/// Scratch for FT refinement rounds — deliberately a separate `RefCell`
/// from [`FtWorkspace`] so `refine` can hold these buffers across the
/// correction solves (which borrow `ws` for their own vectors).
#[derive(Debug, Clone, Default)]
struct FtRefineWorkspace {
    r: Vec<f64>,
    d: Vec<f64>,
    hi: Vec<f64>,
    lo: Vec<f64>,
}

/// Hyper-sparse FTRAN scratch for the Û-arena sweep: `mark` stamps live
/// logicals (heap/`touched` live in logical space, descending order),
/// `smark` stamps slots (the solve's nonzero pattern, grown by the
/// R-eta pass), and `lpat` receives the L-sweep pattern from
/// [`SparseLu::ftran_l_hyper`].
#[derive(Debug, Clone, Default)]
struct FtHyper {
    stamp: u32,
    mark: Vec<u32>,
    heap: BinaryHeap<u32>,
    touched: Vec<u32>,
    sstamp: u32,
    smark: Vec<u32>,
    lpat: Vec<u32>,
    seeds: Vec<u32>,
    /// Min-heap for the ep forward sweep's logical order (see
    /// [`FtFactorization::ep_sweep`]); `mark`/`stamp` serve its marking.
    eheap: BinaryHeap<std::cmp::Reverse<u32>>,
}

/// Resize helper for the persistent scratch buffers.
trait ResizeTo {
    fn resize_to(&mut self, n: usize) -> &mut [f64];
}

impl ResizeTo for Vec<f64> {
    fn resize_to(&mut self, n: usize) -> &mut [f64] {
        if self.len() != n {
            self.resize(n, 0.0);
        }
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pfi::ProductFormFactorization;

    /// Deterministic xorshift64 (self-contained; no env, no global state).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        /// Uniform-ish in [0, 1).
        fn unit(&mut self) -> f64 {
            (self.next() % 10_000) as f64 / 10_000.0
        }
    }

    fn random_col(n: usize, rng: &mut Rng, nnz: usize) -> Vec<(usize, f64)> {
        let mut idx: Vec<usize> = (0..n).collect();
        for i in 0..nnz {
            let j = i + (rng.next() as usize) % (n - i);
            idx.swap(i, j);
        }
        let mut col: Vec<(usize, f64)> = idx[..nnz]
            .iter()
            .map(|&i| {
                (
                    i,
                    if rng.next() & 1 == 0 { -1.0 } else { 1.0 } * (0.5 + rng.unit()),
                )
            })
            .collect();
        col.sort_by_key(|e| e.0);
        col.dedup_by_key(|e| e.0);
        col
    }

    /// Dense Gaussian elimination with partial pivoting (test oracle).
    fn dense_solve(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let mut m = a.to_vec();
        let mut x = b.to_vec();
        for c in 0..n {
            let (mut piv, mut best) = (c, m[c][c].abs());
            for r in c + 1..n {
                if m[r][c].abs() > best {
                    best = m[r][c].abs();
                    piv = r;
                }
            }
            m.swap(c, piv);
            x.swap(c, piv);
            for r in c + 1..n {
                let f = m[r][c] / m[c][c];
                if f != 0.0 {
                    for k in c..n {
                        m[r][k] -= f * m[c][k];
                    }
                    x[r] -= f * x[c];
                }
            }
        }
        for r in (0..n).rev() {
            let mut s = x[r];
            for k in r + 1..n {
                s -= m[r][k] * x[k];
            }
            x[r] = s / m[r][r];
        }
        x
    }

    fn basis_dense(cols: &[Vec<(usize, f64)>], n: usize) -> Vec<Vec<f64>> {
        let mut m = vec![vec![0.0; n]; n];
        for (j, col) in cols.iter().enumerate() {
            for &(i, v) in col {
                m[i][j] = v;
            }
        }
        m
    }

    fn rel_err(got: &[f64], want: &[f64]) -> f64 {
        got.iter()
            .zip(want)
            .map(|(&g, &w)| (g - w).abs() / w.abs().max(1e-9))
            .fold(0.0_f64, f64::max)
    }

    /// Same random pivot sequence through PFI and FT; after every update
    /// both must match the dense rebuild on random FTRAN and BTRAN solves.
    /// `eta_limit` drives both factorizations' refactor cadence. With
    /// `use_bbar`, FT takes the engine's path: the entering column is
    /// FTRAN'd through FT itself (refined) and the update recovers the
    /// spike from that solve — the exact sequence `dual.rs` runs.
    fn differential_run(
        n: usize,
        updates: usize,
        eta_limit: usize,
        seed: u64,
        use_bbar: bool,
    ) -> (u64, usize) {
        let mut rng = Rng(seed);
        let cols: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|j| {
                let mut c = random_col(n, &mut rng, 5);
                c.push((j, 3.0 + rng.unit() * 2.0)); // diagonal dominance
                c.sort_by_key(|e| e.0);
                c.dedup_by_key(|e| e.0);
                c
            })
            .collect();
        let config = PfiConfig {
            eta_limit,
            ..PfiConfig::default()
        };
        let mut pfi = ProductFormFactorization::build_sparse(cols.clone(), config).unwrap();
        let mut ft = FtFactorization::build_sparse(cols.clone(), config).unwrap();
        let mut bd = basis_dense(&cols, n);
        let (mut refactors, mut done, mut attempts) = (0usize, 0usize, 0usize);
        while done < updates {
            attempts += 1;
            assert!(attempts < 10_000, "resampling failed to find pivots");
            let leaving = (rng.next() as usize) % n;
            let entering = random_col(n, &mut rng, 6);
            // Pivot sanity via PFI's pre-update FTRAN (raw, unrefined).
            let mut dense_in = vec![0.0; n];
            for &(i, v) in &entering {
                dense_in[i] = v;
            }
            let mut bbar = vec![0.0; n];
            pfi.solve_b_raw(&dense_in, &mut bbar);
            if bbar[leaving].abs() < 1e-3 {
                continue; // near-singular replacement; resample
            }
            if use_bbar {
                let mut bbar_ft = vec![0.0; n];
                ft.solve_b(&dense_in, &mut bbar_ft);
                ft.update_with_bbar(leaving, entering.clone(), &bbar_ft, 1.0)
                    .unwrap();
            } else {
                ft.update(leaving, entering.clone()).unwrap();
            }
            let nz: Vec<usize> = (0..n).filter(|&i| bbar[i] != 0.0).collect();
            match pfi.update_sparse(leaving, entering.clone(), &bbar, &nz, 1.0) {
                Ok(UpdateStatus::Refactorized) => refactors += 1,
                Ok(_) => {}
                Err(e) => panic!("pfi update failed: {e:?}"),
            }
            // Dense basis bookkeeping.
            for i in 0..n {
                bd[i][leaving] = entering
                    .iter()
                    .find(|&&(j, _)| j == i)
                    .map(|&(_, v)| v)
                    .unwrap_or(0.0);
            }
            done += 1;
            // Differential + ground-truth checks on random vectors.
            for _ in 0..3 {
                let b: Vec<f64> = (0..n).map(|_| rng.unit() * 2.0 - 1.0).collect();
                let truth = dense_solve(&bd, &b);
                let mut xp = vec![0.0; n];
                pfi.solve_b(&b, &mut xp);
                let mut xf = vec![0.0; n];
                ft.solve_b(&b, &mut xf);
                assert!(
                    rel_err(&xf, &truth) < 1e-6,
                    "FT FTRAN vs dense at update {done}: {}",
                    rel_err(&xf, &truth)
                );
                assert!(
                    rel_err(&xp, &truth) < 1e-6,
                    "PFI FTRAN vs dense at update {done}: {}",
                    rel_err(&xp, &truth)
                );
                // BTRAN against the dense transpose.
                let bt: Vec<Vec<f64>> =
                    (0..n).map(|i| (0..n).map(|j| bd[j][i]).collect()).collect();
                let truth_t = dense_solve(&bt, &b);
                let mut xp = vec![0.0; n];
                pfi.solve_bt(&b, &mut xp);
                let mut xf = vec![0.0; n];
                ft.solve_bt(&b, &mut xf);
                assert!(
                    rel_err(&xf, &truth_t) < 1e-6,
                    "FT BTRAN vs dense at update {done}: {}",
                    rel_err(&xf, &truth_t)
                );
                assert!(
                    rel_err(&xp, &truth_t) < 1e-6,
                    "PFI BTRAN vs dense at update {done}: {}",
                    rel_err(&xp, &truth_t)
                );
            }
        }
        (ft.factor_epoch(), refactors)
    }

    #[test]
    fn ft_update_matches_pfi_and_dense() {
        let (ft_epoch, pfi_refactors) = differential_run(40, 60, 12, 0x9E37_79B9_7F4A_7C15, false);
        assert!(ft_epoch >= 3, "expected FT refactors, got {ft_epoch}");
        assert!(pfi_refactors >= 3);
    }

    /// The engine's flow (solve_b → update_with_bbar, refined bbar): the
    /// recovered spike must keep the updated factorization differential
    /// against PFI and the dense rebuild.
    #[test]
    fn ft_update_with_bbar_matches_pfi_and_dense() {
        let (ft_epoch, pfi_refactors) = differential_run(40, 60, 12, 0x9E37_79B9_7F4A_7C15, true);
        assert!(ft_epoch >= 3, "expected FT refactors, got {ft_epoch}");
        assert!(pfi_refactors >= 3);
    }

    /// Hyper-sparse FTRAN arm on a 320-row chain (n ≥ `HYPER_MIN_N`):
    /// after every update, sparse rhs's (≤6 nnz vs the n/16 gate) are
    /// solved by FT — where the entry gate accepts, through the
    /// Gilbert–Peierls L pass, pattern-growth eta pass, and heap Û
    /// sweep — and must match PFI's refined solve. Every 16 updates one
    /// rhs is also checked against the dense rebuild. `hyper_took`
    /// proves the arm actually ran instead of dense fallback throughout.
    #[test]
    fn ft_hyper_ftran_matches_pfi() {
        let n = 320;
        let updates = 48;
        let mut rng = Rng(0x5EED_CAFE_1234_9876);
        let cols: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|j| {
                let mut c = random_col(n, &mut rng, 5);
                c.push((j, 3.0 + rng.unit() * 2.0)); // diagonal dominance
                c.sort_by_key(|e| e.0);
                c.dedup_by_key(|e| e.0);
                c
            })
            .collect();
        let config = PfiConfig {
            eta_limit: 16,
            ..PfiConfig::default()
        };
        let mut pfi = ProductFormFactorization::build_sparse(cols.clone(), config).unwrap();
        let mut ft = FtFactorization::build_sparse(cols.clone(), config).unwrap();
        let mut bd = basis_dense(&cols, n);
        let (mut done, mut attempts) = (0usize, 0usize);
        while done < updates {
            attempts += 1;
            assert!(attempts < 10_000, "resampling failed to find pivots");
            let leaving = (rng.next() as usize) % n;
            let entering = random_col(n, &mut rng, 6);
            let mut dense_in = vec![0.0; n];
            for &(i, v) in &entering {
                dense_in[i] = v;
            }
            let mut bbar = vec![0.0; n];
            pfi.solve_b_raw(&dense_in, &mut bbar);
            if bbar[leaving].abs() < 1e-3 {
                continue; // near-singular replacement; resample
            }
            // Engine flow: FT's own (refined) spike solve feeds the update.
            let mut bbar_ft = vec![0.0; n];
            ft.solve_b(&dense_in, &mut bbar_ft);
            ft.update_with_bbar(leaving, entering.clone(), &bbar_ft, 1.0)
                .unwrap();
            let nz: Vec<usize> = (0..n).filter(|&i| bbar[i] != 0.0).collect();
            pfi.update_sparse(leaving, entering.clone(), &bbar, &nz, 1.0)
                .unwrap();
            for i in 0..n {
                bd[i][leaving] = entering
                    .iter()
                    .find(|&&(j, _)| j == i)
                    .map(|&(_, v)| v)
                    .unwrap_or(0.0);
            }
            done += 1;
            for k in 0..4 {
                let rhs = random_col(n, &mut rng, 6);
                let mut b = vec![0.0; n];
                for &(i, v) in &rhs {
                    b[i] = v;
                }
                let mut xp = vec![0.0; n];
                pfi.solve_b(&b, &mut xp);
                let mut xf = vec![0.0; n];
                ft.solve_b(&b, &mut xf);
                assert!(
                    rel_err(&xf, &xp) < 1e-6,
                    "hyper FT FTRAN vs PFI at update {done}: {}",
                    rel_err(&xf, &xp)
                );
                if done % 16 == 0 && k == 0 {
                    let truth = dense_solve(&bd, &b);
                    assert!(
                        rel_err(&xf, &truth) < 1e-6,
                        "hyper FT FTRAN vs dense at update {done}: {}",
                        rel_err(&xf, &truth)
                    );
                }
            }
        }
        assert!(
            ft.hyper_took.get() > 0,
            "hyper FTRAN arm never ran (entry gate declined every rhs)"
        );
    }

    /// The ep sweep's two arms must agree slot-for-slot: run the hyper
    /// arm (n ≥ `HYPER_MIN_N`, backoff cleared) and the forced-dense arm
    /// for the same leaving position after every update, and cross-check
    /// the shared identity `z[s_p] = 1/pivot_p` at the seed logical.
    #[test]
    fn ft_ep_arms_agree() {
        let n = 320;
        let mut rng = Rng(0xEFEC_DAB8_9674_1235);
        let cols: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|j| {
                let mut c = random_col(n, &mut rng, 5);
                c.push((j, 3.0 + rng.unit() * 2.0));
                c.sort_by_key(|e| e.0);
                c.dedup_by_key(|e| e.0);
                c
            })
            .collect();
        let config = PfiConfig {
            eta_limit: 16,
            ..PfiConfig::default()
        };
        let mut ft = FtFactorization::build_sparse(cols, config).unwrap();
        for it in 0..32 {
            let leaving = (rng.next() as usize) % n;
            // Hyper arm first (clear any ftran backoff; ep counter fresh).
            ft.ep_hyper_skip.set(0);
            let mut ws = std::mem::take(ft.ws.get_mut());
            ft.ep_sweep(leaving, &mut ws);
            let z_hyp = ws.zpat.clone();
            // Forced-dense arm on the same state.
            ft.ep_hyper_skip.set(u32::MAX);
            ft.ep_sweep(leaving, &mut ws);
            let z_dense = ws.zpat.clone();
            *ft.ws.get_mut() = std::mem::take(&mut ws);
            assert_eq!(z_hyp, z_dense, "ep arms differ at iter {it}");
            assert!(!z_hyp.is_empty(), "ep sweep found no nonzeros");
            // Apply a real pivot to move the arena (engine flow).
            let mut dense_in = vec![0.0; n];
            for &(i, v) in &random_col(n, &mut rng, 6) {
                dense_in[i] = v;
            }
            let mut bbar = vec![0.0; n];
            ft.solve_b(&dense_in, &mut bbar);
            let entering: Vec<(usize, f64)> = (0..n)
                .filter(|&i| dense_in[i] != 0.0)
                .map(|i| (i, dense_in[i]))
                .collect();
            let _ = ft.update_with_bbar(leaving, entering, &bbar, 1.0);
        }
    }

    #[test]
    fn ft_aggressive_refactor_path() {
        let (ft_epoch, _) = differential_run(24, 40, 3, 0x0703_C0FF_EE00_1234, true);
        assert!(ft_epoch >= 8, "aggressive refactor cadence, got {ft_epoch}");
    }

    #[test]
    fn ft_lever_flag_parsing() {
        for (value, want) in [
            ("1", true),
            ("on", true),
            ("ON", true),
            (" 1 ", true),
            ("0", false),
            ("off", false),
            ("", false),
        ] {
            assert_eq!(parse_ft_flag(value), want, "RUSOLVE_FT={value:?}");
        }
    }
}
