//! Sparse LU factorization with Markowitz pivot search and threshold
//! partial pivoting.
//!
//! Motivation (Netlib bench): [`crate::dense::DenseLu`] + PFI
//! factors bases densely — O(m²) solves per pivot and no fill control. This
//! module is the sparse engine PFI uses by default:
//!
//! * **Markowitz pivot selection** (H.M. Markowitz, *The Inversion of Sparse
//!   Matrices*, Management Science 17(9), 1971; also Maros, *Computational
//!   Techniques of the Simplex Method*, ch. 7): among numerically admissible
//!   entries, minimize `(row_count−1)·(col_count−1)` fill-in.
//! * **Count-bucket candidate search** (Suhl & Suhl, *Computing sparse LU
//!   factorizations for large-scale linear programming bases*, ORSA J.
//!   Computing 2(4), 1990): active rows and columns live in doubly linked
//!   lists indexed by their nonzero count, so a pivot search visits the
//!   sparsest columns/rows first and stops after a handful of candidates
//!   instead of scanning every row. Count-1 buckets ARE the triangular
//!   (singleton / slack) preprocessing of Suhl & Suhl: column and row
//!   singletons are pivoted first at O(touched) cost, which is where the
//!   bulk of a simplex basis factors. This replaces the former O(m)-per-step
//!   scan that made factorization O(m²) regardless of sparsity (measured
//!   on Netlib).
//! * **Threshold partial pivoting** (Gill/Murray/Wright §4.7; LUSOL uses the
//!   same shape): candidate needs `|a| ≥ pivot_tol · colmax` ⇒ multipliers
//!   bounded by `1/pivot_tol`.
//! * **Sparse factor storage and solves**: L and U are kept as flat
//!   compressed row *and* column arrays; FTRAN/BTRAN run column-oriented
//!   (axpy) sweeps that skip zero entries of the working vector, so a solve
//!   with a sparse right-hand side costs O(entries actually reached) in the
//!   inner loops (the cheap half of Hall & McKinnon 2005 hyper-sparsity; the
//!   symbolic-DFS half is not implemented).
//!
//! Convention: `P·A·Q = L·U` with `row_perm[k]`/`col_perm[k]` the ORIGINAL
//! index at position `k`, so `(P·b)[k] = b[row_perm[k]]` and
//! `(Q·z)[col_perm[k]] = z[k]`.
//!
//! The active submatrix lives as dual sorted-list views — current Schur
//! ROWS and current Schur COLUMNS (LUSOL's row/column-list scheme,
//! Gill/Murray/Wright/Saunders 1985, SOL-85-4) — kept current by
//! right-looking row-outer-product elimination. Pivoting is permutation
//! bookkeeping only: entries never move between rows or columns.
//!
//! Provenance: algorithms from the cited public literature.

use crate::FactorError;
use crate::dense::DenseMatrix;
use crate::hvec::HVec;
use std::cell::RefCell;
use std::collections::BinaryHeap;

/// Profiling: factor_sparse sub-phase counters (appended to the dual
/// engine's RUSOLVE_DUAL_STATS dump; zero cost when off).
#[derive(Debug, Default)]
pub(crate) struct LuFactStats {
    pub fact_ns: std::cell::Cell<f64>,
    pub kernel_ns: std::cell::Cell<f64>,
    pub asm_ns: std::cell::Cell<f64>,
    pub pivot_ns: std::cell::Cell<f64>,
    pub colmax_scan: std::cell::Cell<u64>,
    pub n_recalc: std::cell::Cell<u64>,
    pub n_ins: std::cell::Cell<u64>,
    pub n_rem: std::cell::Cell<u64>,
    pub n_fact: std::cell::Cell<u64>,
    /// Solve census, indexed [FTRAN, BTRAN]:
    /// solves entering the hyper path, skipped by the backoff counter,
    /// refused on rhs count, blown up mid-sweep (dense finish), completed
    /// hyper-sparse, the result-pattern size summed over the hyper
    /// completions, and solves below `HYPER_MIN_N`.
    pub sv_n: [std::cell::Cell<u64>; 2],
    pub sv_skip: [std::cell::Cell<u64>; 2],
    pub sv_cnt: [std::cell::Cell<u64>; 2],
    pub sv_blow: [std::cell::Cell<u64>; 2],
    pub sv_ok: [std::cell::Cell<u64>; 2],
    pub sv_out: [std::cell::Cell<u64>; 2],
    pub sv_dense: [std::cell::Cell<u64>; 2],
    /// L / U strict nonzeros of the most recent factorization.
    pub l_nnz: std::cell::Cell<u64>,
    pub u_nnz: std::cell::Cell<u64>,
    /// Wall of the triangular peel and of the kernel-structure setup
    /// (both inside `fact_ns`; `peel_tr_ns` = the row-wise transpose part
    /// of the peel), the summed kernel size `m - k0` and basis nonzeros.
    pub peel_ns: std::cell::Cell<f64>,
    pub setup_ns: std::cell::Cell<f64>,
    pub kern_rows: std::cell::Cell<u64>,
    pub peel_tr_ns: std::cell::Cell<f64>,
    pub basis_nnz: std::cell::Cell<u64>,
    /// Basis nonzeros of the most recent factorization (fill of the
    /// latest factor = `l_nnz + u_nnz + m - bnnz_last`).
    pub bnnz_last: std::cell::Cell<u64>,
}

#[inline]
fn sv_bump(f: impl Fn(&LuFactStats)) {
    if let Some(s) = lu_fact_stats() {
        f(&s);
    }
}

thread_local! {
    static LU_FACT_STATS: std::rc::Rc<LuFactStats> = std::rc::Rc::default();
}

pub(crate) fn lu_fact_stats() -> Option<std::rc::Rc<LuFactStats>> {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    if !*ON.get_or_init(|| std::env::var("RUSOLVE_DUAL_STATS").is_ok()) {
        return None;
    }
    Some(LU_FACT_STATS.with(|s| s.clone()))
}

/// One-line dump of the factorization counters (profiling).
pub(crate) fn lu_fact_line() -> String {
    if std::env::var("RUSOLVE_DUAL_STATS").is_err() {
        return String::new();
    }
    LU_FACT_STATS.with(|s| {
        let g = |c: &[std::cell::Cell<u64>; 2]| [c[0].get(), c[1].get()];
        let sv = format!(
            " | solves ft/bt: n={:?} skip={:?} cnt={:?} blow={:?} ok={:?} out_mean=[{:.0}, {:.0}] dense={:?} l_nnz={} u_nnz={} peel={:.0}ms (transpose={:.0}ms) setup={:.0}ms kern_mean={:.0} bnnz_mean={:.0} bnnz_last={}",
            g(&s.sv_n),
            g(&s.sv_skip),
            g(&s.sv_cnt),
            g(&s.sv_blow),
            g(&s.sv_ok),
            s.sv_out[0].get() as f64 / s.sv_ok[0].get().max(1) as f64,
            s.sv_out[1].get() as f64 / s.sv_ok[1].get().max(1) as f64,
            g(&s.sv_dense),
            s.l_nnz.get(),
            s.u_nnz.get(),
            s.peel_ns.get() / 1e6,
            s.peel_tr_ns.get() / 1e6,
            s.setup_ns.get() / 1e6,
            s.kern_rows.get() as f64 / s.n_fact.get().max(1) as f64,
            s.basis_nnz.get() as f64 / s.n_fact.get().max(1) as f64,
            s.bnnz_last.get(),
        );
        format!(
            "rusolve lu fact: n_fact={} fact={:.0}ms kernel={:.0}ms pivot={:.0}ms asm={:.0}ms \
             colmax_scan={} recalc={} ins={} rem={}{sv}",
            s.n_fact.get(),
            s.fact_ns.get() / 1e6,
            s.kernel_ns.get() / 1e6,
            s.pivot_ns.get() / 1e6,
            s.asm_ns.get() / 1e6,
            s.colmax_scan.get(),
            s.n_recalc.get(),
            s.n_ins.get(),
            s.n_rem.get(),
        )
    })
}

/// Entry in a sparse factor row/column: (index, value).
type Entry = (usize, f64);

const NONE: usize = usize::MAX;

/// Factor triplet with 4-byte indices: `(orig row, orig col | step,
/// value)` — 16 B instead of 24 B per L/U entry while building.
type Trip = (u32, u32, f64);

/// Persistent workspace of
/// [`SparseLu::factor_sparse_ws`]. Every m-length array of a build lives
/// here between refactorizations; at rest the kernel row/column lists are
/// empty (capacity kept), the count buckets are empty and the `chg` /
/// `dirty` marks are false. `spare` holds the previous factorization, whose
/// buffers the next build reuses. Cloning yields an empty workspace (pure
/// scratch: a cold workspace produces the identical factorization).
///
/// Provenance: build-to-build buffer reuse follows HiGHS (MIT)
/// `util/HFactor.cpp`, where `HFactor::setup` sizes the factor work arrays
/// once and every `build` reuses them; the code here is our own.
#[derive(Debug, Default)]
pub struct FactorWs {
    m: usize,
    rptr: Vec<usize>,
    rent: Vec<(u32, f64)>,
    row_cnt: Vec<u32>,
    col_cnt: Vec<u32>,
    active_col: Vec<bool>,
    active_row: Vec<bool>,
    col_q: Vec<u32>,
    row_q: Vec<u32>,
    u_trips: Vec<Trip>,
    l_trips: Vec<Trip>,
    rows: Vec<Vec<Entry>>,
    cols: Vec<Vec<Entry>>,
    rl: CountLists,
    cl: CountLists,
    colmax: Vec<f64>,
    colmax_row: Vec<usize>,
    chg: Vec<bool>,
    dirty: Vec<bool>,
    cpos: Vec<usize>,
    kr: Vec<usize>,
    kc: Vec<usize>,
    mults: Vec<(usize, f64)>,
    prow: Vec<Entry>,
    chg_list: Vec<usize>,
    fill: Vec<usize>,
    spare: Option<Box<SparseLu>>,
}

impl Clone for FactorWs {
    fn clone(&self) -> Self {
        Self::default()
    }
}

impl FactorWs {
    /// Size the m-length arrays (a no-op at steady state).
    fn ensure(&mut self, m: usize) {
        if self.m == m {
            return;
        }
        let spare = self.spare.take().filter(|s| s.n == m);
        *self = Self {
            m,
            row_cnt: vec![0; m],
            col_cnt: vec![0; m],
            active_col: vec![true; m],
            active_row: vec![true; m],
            rows: vec![Vec::new(); m],
            cols: vec![Vec::new(); m],
            rl: CountLists::new(m),
            cl: CountLists::new(m),
            colmax: vec![0.0; m],
            colmax_row: vec![0; m],
            chg: vec![false; m],
            dirty: vec![false; m],
            cpos: vec![NONE; m],
            spare,
            ..Self::default()
        };
    }

    /// Restore the kernel structures to their rest state over the recorded
    /// kernel ids (also after an error return or the dense-kernel break).
    fn reset_kernel(&mut self) {
        for &i in &self.kr {
            self.rows[i].clear();
            if self.rl.cnt[i] != NONE {
                self.rl.remove(i);
            }
        }
        for &j in &self.kc {
            self.cols[j].clear();
            if self.cl.cnt[j] != NONE {
                self.cl.remove(j);
            }
            self.chg[j] = false;
            self.dirty[j] = false;
        }
        self.kr.clear();
        self.kc.clear();
        self.chg_list.clear();
    }

    /// Hand back a factorization that is no longer needed: its L/U, perm,
    /// inverse and solve-scratch buffers are reused by the next build.
    pub fn recycle(&mut self, lu: Box<SparseLu>) {
        self.spare = Some(lu);
    }
}

/// Flat compressed sparse storage: `ent[ptr[i]..ptr[i+1]]` is line `i`.
#[derive(Debug, Clone, Default)]
struct Flat {
    ptr: Vec<usize>,
    /// Struct-of-arrays entries: 4-byte indices + 8-byte values (12 B per
    /// entry streamed by the solve sweeps instead of a 16 B tuple).
    idx: Vec<u32>,
    val: Vec<f64>,
}

impl Flat {
    /// Counting-sort `(row, col, val)` triplets by row (`by_first`) or by
    /// column (stable) into this Flat's existing buffers. The buffers
    /// are recycled across factorizations and written in place — counts go
    /// to `ptr[line + 2]` so the prefix sum leaves `ptr[line + 1]` as the
    /// scatter cursor of `line`, which ends at the line's end = the next
    /// line's start (no separate cursor array). `idx`/`val` are not
    /// re-zeroed: the scatter writes every slot exactly once.
    fn assign_triplets(&mut self, n: usize, trips: &[Trip], by_first: bool) {
        let ptr = &mut self.ptr;
        ptr.clear();
        ptr.resize(n + 2, 0);
        for &(r, c, _) in trips {
            ptr[(if by_first { r } else { c }) as usize + 2] += 1;
        }
        for i in 2..n + 2 {
            ptr[i] += ptr[i - 1];
        }
        self.idx.resize(trips.len(), 0);
        self.val.resize(trips.len(), 0.0);
        for &(r, c, v) in trips {
            let (line, other) = if by_first { (r, c) } else { (c, r) };
            let f = &mut ptr[line as usize + 1];
            self.idx[*f] = other;
            self.val[*f] = v;
            *f += 1;
        }
        ptr.truncate(n + 1);
    }

    /// [`Self::assign_triplets`] for triplets ALREADY grouped by line
    /// in ascending line order (U by pivot row and L by step are pushed that
    /// way: each elimination step emits its whole U row / L column at once).
    /// A stable counting sort of presorted keys is the identity, so this
    /// linear copy builds the identical Flat without the scatter.
    fn assign_presorted(&mut self, n: usize, trips: &[Trip], by_first: bool) {
        self.ptr.clear();
        self.ptr.push(0);
        self.idx.clear();
        self.val.clear();
        for (p, &(r, c, v)) in trips.iter().enumerate() {
            let (line, other) = if by_first { (r, c) } else { (c, r) };
            debug_assert!(
                line as usize + 1 >= self.ptr.len(),
                "trips not grouped by line"
            );
            while self.ptr.len() <= line as usize {
                self.ptr.push(p);
            }
            self.idx.push(other);
            self.val.push(v);
        }
        self.ptr.resize(n + 1, trips.len());
    }

    #[inline]
    fn line(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (lo, hi) = (self.ptr[i], self.ptr[i + 1]);
        self.idx[lo..hi]
            .iter()
            .zip(&self.val[lo..hi])
            .map(|(&j, &v)| (j as usize, v))
    }

    #[inline]
    fn nnz(&self) -> usize {
        self.val.len()
    }
}

/// Doubly linked lists of ids bucketed by nonzero count (Suhl & Suhl 1990).
#[derive(Debug, Default)]
struct CountLists {
    head: Vec<usize>,
    next: Vec<usize>,
    prev: Vec<usize>,
    cnt: Vec<usize>,
}

impl CountLists {
    fn new(n: usize) -> Self {
        Self {
            head: vec![NONE; n + 2],
            next: vec![NONE; n],
            prev: vec![NONE; n],
            cnt: vec![NONE; n],
        }
    }
    fn insert(&mut self, id: usize, c: usize) {
        self.cnt[id] = c;
        self.prev[id] = NONE;
        let h = self.head[c];
        self.next[id] = h;
        if h != NONE {
            self.prev[h] = id;
        }
        self.head[c] = id;
    }
    fn remove(&mut self, id: usize) {
        let c = self.cnt[id];
        let (p, nx) = (self.prev[id], self.next[id]);
        if p != NONE {
            self.next[p] = nx;
        } else {
            self.head[c] = nx;
        }
        if nx != NONE {
            self.prev[nx] = p;
        }
        self.cnt[id] = NONE;
    }
    fn update(&mut self, id: usize, c: usize) {
        if self.cnt[id] != c {
            self.remove(id);
            self.insert(id, c);
        }
    }
}

/// Pivot candidates examined (columns + rows, once a candidate is in hand)
/// before the Markowitz search stops early. Suhl & Suhl bound the search to
/// a few lists; a small constant keeps the search O(1) in practice.
const SEARCH_LIMIT: usize = 4;

type Cand = (i64, f64, usize, usize); // (score, |a|, orig row, orig col)

/// Bounded Markowitz search over the count buckets. `limit` caps the number
/// of columns+rows examined once a candidate exists (`usize::MAX` = full).
#[allow(clippy::too_many_arguments)]
fn find_pivot(
    rows: &[Vec<Entry>],
    cols: &[Vec<Entry>],
    rl: &CountLists,
    cl: &CountLists,
    colmax: &[f64],
    pivot_tol: f64,
    thresh: f64,
    limit: usize,
    remaining: usize,
) -> Option<Cand> {
    let mut best: Option<Cand> = None;
    let mut examined = 0usize;
    let consider = |best: &mut Option<Cand>, score: i64, v: f64, ri: usize, j: usize| {
        let better = match *best {
            None => true,
            Some((bs, ba, _, _)) => score < bs || (score == bs && v > ba),
        };
        if better {
            *best = Some((score, v, ri, j));
        }
    };
    for c in 1..=remaining {
        // Columns with exactly c active entries.
        let mut j = cl.head[c];
        while j != NONE {
            for &(i, v) in &cols[j] {
                let v = v.abs();
                if v < pivot_tol || v < thresh * colmax[j] {
                    continue;
                }
                let score = (rows[i].len() as i64 - 1) * (c as i64 - 1);
                consider(&mut best, score, v, i, j);
            }
            examined += 1;
            if best.is_some() && examined >= limit {
                return best;
            }
            j = cl.next[j];
        }
        // Rows with exactly c active entries.
        let mut i = rl.head[c];
        while i != NONE {
            for &(j, v) in &rows[i] {
                let v = v.abs();
                if v < pivot_tol || v < thresh * colmax[j] {
                    continue;
                }
                let score = (c as i64 - 1) * (cols[j].len() as i64 - 1);
                consider(&mut best, score, v, i, j);
            }
            examined += 1;
            if best.is_some() && examined >= limit {
                return best;
            }
            i = rl.next[i];
        }
        // Unseen entries sit in rows AND columns of count > c ⇒ score ≥ c².
        if let Some((bs, _, _, _)) = best
            && bs <= (c * c) as i64
        {
            return best;
        }
    }
    best
}

/// Minimum / maximum remaining dimension and minimum density for the dense
/// kernel switch (see the call site).
const DENSE_MIN: usize = 16;
/// Kernel elimination: a row/column whose length exceeds this many
/// times its update count is updated in place (binary search per update)
/// instead of rebuilt by a merge.
const MERGE_RATIO: usize = 8;
const DENSE_MAX: usize = 4000;
const DENSE_DENSITY: f64 = 0.25;

/// Finish the factorization from step `k0` on the remaining active block
/// with dense LU + row partial pivoting. `rows` holds the CURRENT Schur rows
/// (retired rows are empty); the active rows are those not yet in `row_perm`.
#[allow(clippy::too_many_arguments)]
fn dense_kernel_finish(
    k0: usize,
    rows: &[Vec<Entry>],
    (kr, kc): (&[usize], &[usize]),
    (active_row, active_col): (&[bool], &[bool]),
    cpos: &mut [usize],
    pivot_tol: f64,
    row_perm: &mut Vec<usize>,
    col_perm: &mut Vec<usize>,
    u_diag: &mut [f64],
    u_trips: &mut Vec<Trip>,
    l_trips: &mut Vec<Trip>,
    work: &mut usize,
) -> Result<(), FactorError> {
    let m = rows.len();
    let kk = m - k0;
    // The active rows/columns are the still-active kernel ids, in the
    // same ascending order the former O(m) scans produced (every non-kernel
    // id was pivoted by the peel; the kernel loop clears `active_row`).
    let act_rows: Vec<usize> = kr.iter().copied().filter(|&i| active_row[i]).collect();
    let act_cols: Vec<usize> = kc.iter().copied().filter(|&j| active_col[j]).collect();
    debug_assert_eq!(act_rows.len(), kk);
    debug_assert_eq!(act_cols.len(), kk);
    debug_assert!(row_perm.len() == k0);
    for (t, &c) in act_cols.iter().enumerate() {
        cpos[c] = t;
    }
    // Row-major dense block, one row per active row.
    let mut d = vec![0.0f64; kk * kk];
    for (ri, &i) in act_rows.iter().enumerate() {
        for &(j, v) in &rows[i] {
            d[ri * kk + cpos[j]] = v;
        }
    }
    let mut rp: Vec<usize> = (0..kk).collect(); // dense-row order (indices into act_rows)
    for t in 0..kk {
        // Partial pivoting: largest entry in column t among rows t..kk.
        let mut p = t;
        let mut best = d[rp[t] * kk + t].abs();
        for q in t + 1..kk {
            let a = d[rp[q] * kk + t].abs();
            if a > best {
                best = a;
                p = q;
            }
        }
        if best < pivot_tol {
            return Err(FactorError::Singular {
                column: act_cols[t],
                pivot: best,
            });
        }
        rp.swap(t, p);
        let prow_i = rp[t];
        let pivot = d[prow_i * kk + t];
        let k = k0 + t;
        row_perm.push(act_rows[prow_i]);
        col_perm.push(act_cols[t]);
        u_diag[k] = pivot;
        for j in t + 1..kk {
            let v = d[prow_i * kk + j];
            if v != 0.0 {
                u_trips.push((act_rows[prow_i] as u32, act_cols[j] as u32, v));
            }
        }
        // Split the borrow: copy the pivot row tail once.
        let ptail: Vec<f64> = d[prow_i * kk + t + 1..prow_i * kk + kk].to_vec();
        for q in t + 1..kk {
            let ri = rp[q];
            let base = ri * kk;
            let a = d[base + t];
            if a == 0.0 {
                continue;
            }
            let mlt = a / pivot;
            l_trips.push((act_rows[ri] as u32, k as u32, mlt));
            let row = &mut d[base + t + 1..base + kk];
            for (x, &pv) in row.iter_mut().zip(&ptail) {
                *x -= mlt * pv;
            }
        }
    }
    *work += kk * kk * kk / 30 + kk * kk;
    Ok(())
}

/// Workspace for the hyper-sparse solve path (Gilbert & Peierls 1988; Hall &
/// McKinnon 2005 "hyper-sparsity"). `y` is all zero between calls.
#[derive(Debug, Clone, Default)]
struct HyperScratch {
    y: Vec<f64>,
    mark: Vec<u32>,
    stamp: u32,
    heap: BinaryHeap<u32>,
    seeds: Vec<u32>,
    out: Vec<u32>,
}

impl HyperScratch {
    fn new(n: usize) -> Self {
        Self {
            y: vec![0.0; n],
            mark: vec![0; n],
            ..Self::default()
        }
    }
}

/// Hyper-sparse path thresholds: skip below `HYPER_MIN_N` rows; use the
/// sparse sweeps only while the rhs / reach stays under `n / HYPER_DIV`.
pub(crate) const HYPER_MIN_N: usize = 256;
/// Dense fallbacks to take before trying the hyper-sparse path again.
pub(crate) const HYPER_BACKOFF: u32 = 8;
pub(crate) const HYPER_DIV: usize = 16;

/// Hyper-path backoff policy, lever
/// `RUSOLVE_LU_HSKIP`:
/// - `0`: legacy — ONE counter shared by FTRAN and BTRAN, armed for
///   [`HYPER_BACKOFF`] solves by a reach blow-up AND by an rhs-count
///   refusal. In the dual the solve kinds interleave (β BTRAN, a_q /
///   DSE / flip FTRANs, refinement corrections), so one dense DSE FTRAN
///   sent the next eight — mostly hyper-sparse — solves of every kind to
///   the O(m) dense sweeps (census: 30–90 % of all solves on pds-40 /
///   fome13 were backoff skips).
/// - `1`: per-direction counters, armed by blow-ups only (a
///   count refusal costs nothing to detect: O(pattern) on the HVec entry).
///   HiGHS keeps a separate expected density per solve kind for the same
///   reason (`HEkkDual` `row_ep_density` / `col_aq_density` / DSE density
///   feeding `HFactor::ftranL/btranU` `kHyper*` gates, util/HFactor.cpp
///   ~:1545–1776, MIT, read as reference).
/// - `2`: no backoff at all.
/// - `3` (default): per-direction, armed only by the SECOND consecutive blow-up in
///   that direction (an isolated blow-up already finishes densely at about
///   dense cost — `sparse_sweep` continues in place — so skipping the next
///   solves only forfeits their hyper-sparse chance; a run of blow-ups is
///   the consistently-dense regime the backoff exists for).
///
/// Every mode only picks between bit-identical paths. Measured (
/// instructions at equal iterations, mode 0 → 3): pds-40 −12 %, pds-100
/// −32 %, fome13 −23 %, stormG2 −11 %; 47 relaxed MIPLIB LPs geomean 0.995,
/// 0 path diffs.
pub(crate) fn hskip_mode() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        std::env::var("RUSOLVE_LU_HSKIP")
            .ok()
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(HSKIP_DEFAULT)
    })
}

/// Default of [`hskip_mode`].
const HSKIP_DEFAULT: u8 = 3;

/// Dense triangular sweep (same loop as the dense solves), used for the
/// remainder after a hyper-sparse reach blow-up.
fn dense_sweep(flat: &Flat, diag: Option<&[f64]>, asc: bool, y: &mut [f64]) {
    let n = y.len();
    for step in 0..n {
        let k = if asc { step } else { n - 1 - step };
        let yk = y[k];
        if yk != 0.0 {
            let zk = match diag {
                Some(d) => {
                    let z = yk / d[k];
                    y[k] = z;
                    z
                }
                None => yk,
            };
            for (i, v) in flat.line(k) {
                y[i] -= v * zk;
            }
        }
    }
}

/// One triangular sweep restricted to the reachable nonzeros of `hs.seeds`,
/// processing columns in ascending (`asc`) or descending index order with a
/// priority queue. Because the dense sweep also visits columns in index
/// order and skips exact zeros, the floating-point operation sequence is
/// IDENTICAL to the dense loop (bit-exact). Optional division by `diag`
/// (the U sweeps). On reach growth beyond `thr` the sweep finishes densely
/// from the current position; returns `true` when it did so (then
/// `hs.out` is no longer the touched set).
fn sparse_sweep(
    flat: &Flat,
    diag: Option<&[f64]>,
    asc: bool,
    n: usize,
    thr: usize,
    hs: &mut HyperScratch,
) -> bool {
    let key = |k: u32| if asc { !k } else { k };
    hs.stamp = hs.stamp.wrapping_add(1);
    if hs.stamp == 0 {
        hs.mark.fill(0);
        hs.stamp = 1;
    }
    let stamp = hs.stamp;
    hs.heap.clear();
    hs.out.clear();
    for idx in 0..hs.seeds.len() {
        let k = hs.seeds[idx];
        if hs.mark[k as usize] != stamp {
            hs.mark[k as usize] = stamp;
            hs.heap.push(key(k));
            hs.out.push(k);
        }
    }
    while let Some(kk) = hs.heap.pop() {
        let k = key(kk) as usize;
        let yk = hs.y[k];
        if yk != 0.0 {
            let zk = match diag {
                Some(d) => {
                    let z = yk / d[k];
                    hs.y[k] = z;
                    z
                }
                None => yk,
            };
            for (i, v) in flat.line(k) {
                hs.y[i] -= v * zk;
                if hs.mark[i] != stamp {
                    hs.mark[i] = stamp;
                    hs.heap.push(key(i as u32));
                    hs.out.push(i as u32);
                }
            }
        }
        if hs.out.len() > thr {
            // Dense finish from the next index in sweep order.
            let rest: Box<dyn Iterator<Item = usize>> = if asc {
                Box::new(k + 1..n)
            } else {
                Box::new((0..k).rev())
            };
            for kd in rest {
                let ykd = hs.y[kd];
                if ykd != 0.0 {
                    let zk = match diag {
                        Some(d) => {
                            let z = ykd / d[kd];
                            hs.y[kd] = z;
                            z
                        }
                        None => ykd,
                    };
                    for (i, v) in flat.line(kd) {
                        hs.y[i] -= v * zk;
                    }
                }
            }
            return true;
        }
    }
    false
}

/// Column-form seed of the base `U` for the Forrest–Tomlin arena
/// (`lu_ft.rs`). Entries are in permuted slot space; the
/// implicit diagonal is carried in `u_pivot_value` (one pivot per logical
/// position).
#[derive(Debug, Clone)]
pub(crate) struct FtLuSeed {
    pub u_start: Vec<u32>,
    pub u_last: Vec<u32>,
    pub u_index: Vec<u32>,
    pub u_value: Vec<f64>,
    pub u_pivot_value: Vec<f64>,
}

/// LU factors `P·A·Q = L·U` in sparse form.
#[derive(Debug, Clone)]
pub struct SparseLu {
    n: usize,
    row_perm: Vec<usize>,
    col_perm: Vec<usize>,
    /// L rows (strict lower multipliers; unit diagonal implied).
    l_crs: Flat,
    /// L columns (strict lower).
    l_csc: Flat,
    /// U rows, strict upper part only (diagonal in `u_diag`).
    u_crs: Flat,
    /// U columns, strict upper part only.
    u_csc: Flat,
    /// Diagonal of U (explicit; division happens once per row).
    u_diag: Vec<f64>,
    /// Persistent solve temporary.
    /// A solve never nests another solve, so one borrow per call is sound.
    scratch: RefCell<Vec<f64>>,
    /// Inverse permutations (original index -> permuted position) for the
    /// hyper-sparse entry of FTRAN/BTRAN.
    row_inv: Vec<u32>,
    col_inv: Vec<u32>,
    /// Hyper-sparse solve workspace (all-zero `y` at rest).
    hs: RefCell<HyperScratch>,
    /// Solves left to send straight to the dense sweep after a hyper-sparse
    /// attempt turned out dense: the nonzero-count scan and the wasted partial
    /// sweep cost more than they save on workloads (dual simplex) whose
    /// results are dense. Bit-neutral: both paths agree exactly.
    hyper_skip: [std::cell::Cell<u32>; 2],
    /// Mode 3: the previous hyper attempt in this direction blew up
    /// (backoff arms only on the second consecutive blow-up).
    hyper_blew: [std::cell::Cell<bool>; 2],
    /// Work performed by the factorization (entries touched in elimination
    /// plus factor nonzeros): the "build tick" PFI's refactorization trigger
    /// compares accumulated solve cost against (HiGHS's synthetic-clock
    /// idea, HEkk.cpp updateFactor / HFactorConst.h).
    build_work: usize,
}

impl SparseLu {
    /// Factor the matrix whose COLUMNS are `columns` (column-major dense).
    /// Convenience wrapper over [`SparseLu::factor_sparse`].
    pub fn factor(columns: &[Vec<f64>], pivot_tol: f64) -> Result<Self, FactorError> {
        let m = columns.len();
        assert!(
            columns.iter().all(|c| c.len() == m),
            "factor: non-square or ragged basis"
        );
        let sp: Vec<Vec<Entry>> = columns
            .iter()
            .map(|c| {
                c.iter()
                    .enumerate()
                    .filter_map(|(i, &v)| (v != 0.0).then_some((i, v)))
                    .collect()
            })
            .collect();
        Self::factor_sparse(&sp, pivot_tol)
    }

    /// Factor the square matrix whose columns are the sparse lists
    /// `columns[j] = [(row, value), …]`, **ascending by row** with no
    /// duplicate rows (exact zeros are skipped).
    ///
    /// `pivot_tol` is the threshold-pivoting tolerance (admissible entries
    /// need `|a| ≥ pivot_tol · colmax`); multipliers are bounded by
    /// `1/pivot_tol`. Returns [`FactorError::Singular`] when a remaining
    /// column is structurally/numerically zero.
    ///
    /// Cold-workspace wrapper over [`SparseLu::factor_sparse_ws`].
    pub fn factor_sparse(columns: &[Vec<Entry>], pivot_tol: f64) -> Result<Self, FactorError> {
        Self::factor_sparse_ws(columns, None, pivot_tol, &mut FactorWs::default())
    }

    /// [`SparseLu::factor_sparse`] with a
    /// caller-owned persistent [`FactorWs`] and, optionally, the caller's
    /// row-wise copy of the same matrix.
    ///
    /// `rows`, when given, must be the exact row mirror of `columns`:
    /// `rows[i]` lists `(j, a_ij)` **ascending by column**, and `columns`
    /// must then be free of explicit zeros. The triangular peel then reads
    /// the mirror instead of building the transpose — HiGHS keeps the same
    /// basis-matrix row copy across `HFactor::build` calls (`setupMatrix` /
    /// `buildSimple` consume the persistent `MRindex` arrays, util/HFactor.cpp,
    /// MIT, read as reference). The flat transpose the peel builds otherwise
    /// visits each row's entries in the same ascending-column order, so
    /// both sources yield the identical pivot sequence, triplet order and
    /// values (bit-identical factors).
    ///
    /// The workspace keeps every m-length array of the build (peel counts,
    /// queues, kernel row/column lists, count buckets, triplets) between
    /// calls, and recycles the buffers of the previous factorization handed
    /// back via [`FactorWs::recycle`]: at large m the build was dominated by
    /// fresh multi-MB zeroed allocations (page faults), not by arithmetic.
    /// Kernel structures are restored to their all-empty rest state over the
    /// kernel ids only (O(kernel), also on the error paths).
    pub fn factor_sparse_ws(
        columns: &[Vec<Entry>],
        rows: Option<&[Vec<(u32, f64)>]>,
        pivot_tol: f64,
        ws: &mut FactorWs,
    ) -> Result<Self, FactorError> {
        let m = columns.len();
        if m == 0 {
            return Ok(Self::empty());
        }
        debug_assert!(rows.is_none_or(|r| r.len() == m));
        ws.ensure(m);
        let res = Self::factor_core(columns, rows, pivot_tol, ws);
        ws.reset_kernel();
        res
    }

    /// The `m == 0` factorization (also the donor of empty buffers).
    fn empty() -> Self {
        Self {
            n: 0,
            row_perm: vec![],
            col_perm: vec![],
            l_crs: Flat::default(),
            l_csc: Flat::default(),
            u_crs: Flat::default(),
            u_csc: Flat::default(),
            u_diag: vec![],
            scratch: RefCell::new(Vec::new()),
            row_inv: vec![],
            col_inv: vec![],
            hs: RefCell::new(HyperScratch::default()),
            hyper_skip: Default::default(),
            hyper_blew: Default::default(),
            build_work: 0,
        }
    }

    fn factor_core(
        columns: &[Vec<Entry>],
        mirror: Option<&[Vec<(u32, f64)>]>,
        pivot_tol: f64,
        ws: &mut FactorWs,
    ) -> Result<Self, FactorError> {
        let m = columns.len();
        let fst = lu_fact_stats();
        let ft0 = fst.as_ref().map(|_| std::time::Instant::now());

        // Recycled buffers of a previous factorization (or fresh empties).
        let SparseLu {
            mut row_perm,
            mut col_perm,
            mut l_crs,
            mut l_csc,
            mut u_crs,
            mut u_csc,
            u_diag: mut u_diag_step,
            mut scratch,
            mut row_inv,
            mut col_inv,
            mut hs,
            ..
        } = ws.spare.take().map_or_else(Self::empty, |b| *b);
        let FactorWs {
            rptr,
            rent,
            row_cnt,
            col_cnt,
            active_col,
            active_row,
            col_q,
            row_q,
            u_trips,
            l_trips,
            rows,
            cols,
            rl,
            cl,
            colmax,
            colmax_row,
            chg,
            dirty,
            cpos,
            kr,
            kc,
            mults,
            prow,
            chg_list,
            fill,
            ..
        } = ws;

        u_trips.clear(); // (orig row, orig col, v)
        l_trips.clear(); // (orig row, step, mlt)
        // Not re-zeroed: a successful build writes all m pivots (and every
        // entry of the inverse permutations below).
        u_diag_step.resize(m, 0.0);
        row_perm.clear(); // position → original
        col_perm.clear();
        row_perm.reserve(m);
        col_perm.reserve(m);
        let mut work = 0usize;

        // ---- Triangular preprocessing (Suhl & Suhl 1990 "triangularization";
        // HiGHS HFactor::buildSimple): peel column singletons (pivot = the
        // lone entry, U row = the rest of its row, no multipliers) and row
        // singletons (pivot = the lone entry, multipliers = rest of its
        // column / pivot, empty U row). Neither changes any remaining Schur
        // value, so the kernel keeps the original coefficients and no
        // per-row/column lists are built for the (typically dominant)
        // triangular part. Row access is a flat CSR copy (`rptr`/`rent`,
        // each row ascending by column). With the caller's mirror it
        // is a sequential row-by-row copy (independent reads, sequential
        // writes) instead of the column-driven scatter transpose; the peel
        // deliberately reads the flat copy, not the mirror's per-row heap
        // Vecs — basic slacks mostly sit at position == row, so the LIFO
        // peel walks rows near-sequentially and the contiguous copy is
        // cache-friendly (measured: peeling straight off the mirror was
        // ~40 % slower than off the copy, `factor_ws_timing`). ----
        let nnz_total;
        match mirror {
            Some(mr) => {
                let mut nz = 0usize;
                for (j, col) in columns.iter().enumerate() {
                    debug_assert!(col.iter().all(|&(_, v)| v != 0.0));
                    col_cnt[j] = col.len() as u32;
                    nz += col.len();
                }
                rptr.clear();
                rptr.push(0);
                rent.clear();
                for (i, r) in mr.iter().enumerate() {
                    row_cnt[i] = r.len() as u32;
                    rent.extend_from_slice(r);
                    rptr.push(rent.len());
                }
                nnz_total = nz;
                debug_assert_eq!(rent.len(), nz);
            }
            None => {
                rptr.clear();
                rptr.resize(m + 1, 0);
                for col in columns {
                    for &(i, v) in col {
                        if v != 0.0 {
                            rptr[i + 1] += 1;
                        }
                    }
                }
                for i in 0..m {
                    rptr[i + 1] += rptr[i];
                }
                nnz_total = rptr[m];
                fill.clear();
                fill.extend_from_slice(&rptr[..m]);
                rent.clear();
                rent.resize(nnz_total, (0, 0.0));
                for i in 0..m {
                    row_cnt[i] = (rptr[i + 1] - rptr[i]) as u32;
                }
                for (j, col) in columns.iter().enumerate() {
                    let mut c = 0u32;
                    for &(i, v) in col {
                        if v != 0.0 {
                            rent[fill[i]] = (j as u32, v);
                            fill[i] += 1;
                            c += 1;
                        }
                    }
                    col_cnt[j] = c;
                }
            }
        }
        let (rptr, rent) = (&*rptr, &*rent);
        let row_of = |i: usize| &rent[rptr[i]..rptr[i + 1]];
        if let (Some(s), Some(t0)) = (fst.as_ref(), ft0) {
            s.peel_tr_ns
                .set(s.peel_tr_ns.get() + t0.elapsed().as_nanos() as f64);
            s.basis_nnz.set(s.basis_nnz.get() + nnz_total as u64);
            s.bnnz_last.set(nnz_total as u64);
        }
        active_col.fill(true);
        active_row.fill(true);
        let eff_tol0 = pivot_tol.max(1e-4);
        col_q.clear();
        row_q.clear();
        col_q.extend((0..m as u32).filter(|&j| col_cnt[j as usize] == 1));
        row_q.extend((0..m as u32).filter(|&i| row_cnt[i as usize] == 1));
        let mut k0 = 0usize;
        loop {
            if let Some(j) = col_q.pop() {
                let j = j as usize;
                if !active_col[j] || col_cnt[j] != 1 {
                    continue;
                }
                let Some(&(i, v)) = columns[j].iter().find(|&&(i, v)| v != 0.0 && active_row[i])
                else {
                    continue;
                };
                if v.abs() < pivot_tol {
                    continue;
                }
                row_perm.push(i);
                col_perm.push(j);
                u_diag_step[k0] = v;
                let row = row_of(i);
                for &(jj, a) in row {
                    let jj = jj as usize;
                    if jj != j && active_col[jj] {
                        u_trips.push((i as u32, jj as u32, a));
                        col_cnt[jj] -= 1;
                        if col_cnt[jj] == 1 {
                            col_q.push(jj as u32);
                        }
                    }
                }
                work += row.len() + 1;
                active_col[j] = false;
                active_row[i] = false;
                k0 += 1;
            } else if let Some(i) = row_q.pop() {
                let i = i as usize;
                if !active_row[i] || row_cnt[i] != 1 {
                    continue;
                }
                let Some(&(j, v)) = row_of(i).iter().find(|&&(j, _)| active_col[j as usize]) else {
                    continue;
                };
                let j = j as usize;
                if v.abs() < pivot_tol {
                    continue;
                }
                // Multipliers v'/v are bounded by 1/eff_tol only if v is not
                // tiny relative to its column: threshold test.
                let mut cmax = 0.0f64;
                for &(i2, v2) in &columns[j] {
                    if active_row[i2] {
                        cmax = cmax.max(v2.abs());
                    }
                }
                if v.abs() < eff_tol0 * cmax {
                    continue;
                }
                row_perm.push(i);
                col_perm.push(j);
                u_diag_step[k0] = v;
                for &(i2, v2) in &columns[j] {
                    if i2 != i && v2 != 0.0 && active_row[i2] {
                        l_trips.push((i2 as u32, k0 as u32, v2 / v));
                        row_cnt[i2] -= 1;
                        if row_cnt[i2] == 1 {
                            row_q.push(i2 as u32);
                        }
                    }
                }
                work += columns[j].len() + 1;
                active_col[j] = false;
                active_row[i] = false;
                k0 += 1;
            } else {
                break;
            }
        }

        let pe0 = fst.as_ref().map(|_| std::time::Instant::now());
        if let (Some(s), Some(t0), Some(t1)) = (fst.as_ref(), ft0, pe0) {
            s.peel_ns
                .set(s.peel_ns.get() + t1.duration_since(t0).as_nanos() as f64);
            s.kern_rows.set(s.kern_rows.get() + (m - k0) as u64);
        }
        // ---- Kernel: dual sorted-list views of what remains. The kernel
        // ids are recorded so `FactorWs::reset_kernel` restores the rest
        // state in O(kernel). Insertion orders (rows by col ascending, the
        // count buckets by id ascending) are the earlier ones. ----
        kr.clear();
        kc.clear();
        kr.extend((0..m).filter(|&i| active_row[i]));
        kc.extend((0..m).filter(|&j| active_col[j]));
        for &i in kr.iter() {
            debug_assert!(rows[i].is_empty());
            rows[i].reserve(row_cnt[i] as usize);
        }
        for &j in kc.iter() {
            let cj = &mut cols[j];
            debug_assert!(cj.is_empty());
            cj.reserve(col_cnt[j] as usize);
            for &(i, v) in &columns[j] {
                if v != 0.0 && active_row[i] {
                    rows[i].push((j, v)); // j outer ⇒ rows sorted by col
                    cj.push((i, v));
                }
            }
        }
        for &i in kr.iter() {
            rl.insert(i, rows[i].len());
        }
        for &j in kc.iter() {
            cl.insert(j, cols[j].len());
        }
        // Per-column maxima over ACTIVE rows + holder ORIGINAL row (exact:
        // recomputed over the column list whenever the column is touched).
        for &j in kc.iter() {
            recalc_colmax(cols, j, colmax, colmax_row);
        }

        mults.clear(); // (orig row, mlt)
        let mut mbuf: Vec<Entry> = Vec::new(); // merge scratch
        prow.clear(); // pivot row snapshot
        chg_list.clear(); // columns modified this step (`chg` marks)
        // `dirty` marks columns whose holder shrank/vanished this step —
        // the only ones needing an exact `recalc_colmax` (see `touch_colmax`).
        let mut active_nnz: usize = kr.iter().map(|&i| rows[i].len()).sum();

        let kt0 = fst.as_ref().map(|_| std::time::Instant::now());
        if let (Some(s), Some(t0), Some(t1)) = (fst.as_ref(), pe0, kt0) {
            s.setup_ns
                .set(s.setup_ns.get() + t1.duration_since(t0).as_nanos() as f64);
        }
        for k in k0..m {
            // An active row/column with no entries left ⇒ singular.
            if cl.head[0] != NONE {
                return Err(FactorError::Singular {
                    column: cl.head[0],
                    pivot: 0.0,
                });
            }
            if rl.head[0] != NONE {
                // Structurally empty row: report any still-active column.
                let j = (0..m).find(|&j| active_col[j]).unwrap_or(0);
                return Err(FactorError::Singular {
                    column: j,
                    pivot: 0.0,
                });
            }
            let remaining = m - k;

            // ---- Dense-kernel switch: once the Schur complement is dense,
            // sorted-list elimination costs ~10x a dense flop (binary search
            // + Vec shifts per entry). Finish with a dense partial-pivoting
            // LU on the remaining block (LUSOL's dense-bump switch, Gill et
            // al. SOL-85-4 / Suhl & Suhl 1990 "dense window"). ----
            if (DENSE_MIN..=DENSE_MAX).contains(&remaining)
                && (active_nnz as f64) >= DENSE_DENSITY * (remaining as f64) * (remaining as f64)
            {
                dense_kernel_finish(
                    k,
                    rows,
                    (kr, kc),
                    (active_row, active_col),
                    cpos,
                    pivot_tol,
                    &mut row_perm,
                    &mut col_perm,
                    &mut u_diag_step,
                    u_trips,
                    l_trips,
                    &mut work,
                )?;
                break;
            }

            // ---- Markowitz pivot search over the count buckets. TWO-STAGE
            // threshold (LUSOL's pivot threshold, Harris 1990 / Gould et al.
            // 2007 — UT typically 1e-4..1e-1): Stage 1 demands multiplier
            // safety (|a| ≥ max(pivot_tol, 1e-4)·colmax) with a bounded
            // search; if nothing passes, widen to a full search at that
            // threshold, then at the caller's floor; the largest-entry
            // relaxation below still guards true singularity. ----
            let eff_tol = pivot_tol.max(1e-4);
            let pt0 = fst.as_ref().map(|_| std::time::Instant::now());
            let mut best = find_pivot(
                rows,
                cols,
                rl,
                cl,
                colmax,
                pivot_tol,
                eff_tol,
                SEARCH_LIMIT,
                remaining,
            );
            if best.is_none() {
                best = find_pivot(
                    rows,
                    cols,
                    rl,
                    cl,
                    colmax,
                    pivot_tol,
                    pivot_tol,
                    usize::MAX,
                    remaining,
                );
            }
            if let (Some(s), Some(t)) = (fst.as_ref(), pt0) {
                s.pivot_ns
                    .set(s.pivot_ns.get() + t.elapsed().as_nanos() as f64);
            }
            // Fallback: nothing met the threshold — pivot on the largest
            // remaining column entry (documented threshold-pivot relaxation)
            // unless that entry is below the absolute singularity tolerance
            // (same semantics as DenseLu: |pivot| < pivot_tol ⇒ singular).
            let (pr, pc) = match best {
                Some((_, _, i, j)) => (i, j),
                None => {
                    let mut jbest: Option<usize> = None;
                    for j in 0..m {
                        if active_col[j]
                            && !cols[j].is_empty()
                            && colmax[j] > 0.0
                            && jbest.is_none_or(|jb| colmax[j] > colmax[jb])
                        {
                            jbest = Some(j);
                        }
                    }
                    match jbest {
                        Some(j) if colmax[j] >= pivot_tol => (colmax_row[j], j),
                        Some(j) => {
                            return Err(FactorError::Singular {
                                column: j,
                                pivot: colmax[j],
                            });
                        }
                        None => {
                            return Err(FactorError::Singular {
                                column: (0..m).find(|&j| active_col[j]).unwrap_or(0),
                                pivot: 0.0,
                            });
                        }
                    }
                }
            };

            // ---- Pivot value: the exact stored Schur entry. ----
            let pivot = match find_entry(&rows[pr], pc) {
                Some(t) => rows[pr][t].1,
                None => {
                    return Err(FactorError::Singular {
                        column: pc,
                        pivot: 0.0,
                    });
                }
            };
            if pivot == 0.0 {
                return Err(FactorError::Singular {
                    column: pc,
                    pivot: 0.0,
                });
            }

            // ---- Position bookkeeping: step k pivots (pr, pc). ----
            row_perm.push(pr);
            col_perm.push(pc);
            u_diag_step[k] = pivot;

            // ---- Retire the pivot row: its CURRENT row IS the final U row
            // (right-looking: later steps only update multiplier rows). ----
            prow.clear();
            for &(j, v) in &rows[pr] {
                if j != pc {
                    prow.push((j, v));
                    u_trips.push((pr as u32, j as u32, v));
                }
            }
            active_nnz -= rows[pr].len();
            rows[pr].clear();
            rl.remove(pr);
            active_row[pr] = false;

            // ---- Multipliers from the pivot column. ----
            mults.clear();
            for &(i, v) in &cols[pc] {
                if i != pr {
                    mults.push((i, v / pivot));
                }
            }
            for &(i, mlt) in mults.iter() {
                l_trips.push((i as u32, k as u32, mlt));
            }
            work += prow.len() + mults.len() + 1;

            // ---- Right-looking outer-product elimination S −= m·pᵀ, dual
            // views synced. Every touched
            // row and column is rebuilt by ONE sorted merge — O(len +
            // updates); sparse vector addition over sorted packed lists,
            // Duff, Erisman & Reid, Direct Methods for Sparse Matrices, 2nd
            // ed. — instead of a binary search plus a Vec insert/remove
            // shift per updated entry (the memmove hot spot of qap10's
            // kernel). A short update list on a long line keeps the
            // in-place path ([`MERGE_RATIO`]). Values (a − mlt·pv), the
            // exact-zero drops and the per-column `touch_colmax` order
            // (mults order) are those of the former per-entry loop, so the
            // factors are bit-identical. The pc entry is CONSUMED (the L
            // multiplier replaces it). ----
            let (mut n_ins, mut n_rem) = (0u64, 0u64);
            for &(i, mlt) in mults.iter() {
                let ui = &mut rows[i];
                if prow.is_empty() || prow.len() * MERGE_RATIO < ui.len() {
                    if let Some(t) = find_entry(ui, pc) {
                        ui.remove(t);
                        active_nnz -= 1;
                        n_rem += 1;
                    }
                    if prow.is_empty() {
                        continue;
                    }
                    work += prow.len() + ui.len();
                    for &(j, pv) in prow.iter() {
                        let upd = mlt * pv;
                        let p = ui.partition_point(|&(k, _)| k < j);
                        if p < ui.len() && ui[p].0 == j {
                            let nv = ui[p].1 - upd;
                            if nv == 0.0 {
                                active_nnz -= 1;
                                ui.remove(p);
                                n_rem += 2;
                            } else {
                                ui[p].1 = nv;
                            }
                        } else {
                            active_nnz += 1;
                            ui.insert(p, (j, -upd));
                            n_ins += 2;
                        }
                    }
                } else {
                    mbuf.clear();
                    let mut had_pc = false;
                    let mut a = 0usize;
                    for &(j, pv) in prow.iter() {
                        while a < ui.len() && ui[a].0 < j {
                            if ui[a].0 == pc {
                                had_pc = true;
                            } else {
                                mbuf.push(ui[a]);
                            }
                            a += 1;
                        }
                        let upd = mlt * pv;
                        if a < ui.len() && ui[a].0 == j {
                            let nv = ui[a].1 - upd;
                            a += 1;
                            if nv == 0.0 {
                                active_nnz -= 1;
                                n_rem += 2;
                            } else {
                                mbuf.push((j, nv));
                            }
                        } else {
                            active_nnz += 1;
                            n_ins += 2;
                            mbuf.push((j, -upd));
                        }
                    }
                    for &e in &ui[a..] {
                        if e.0 == pc {
                            had_pc = true;
                        } else {
                            mbuf.push(e);
                        }
                    }
                    if had_pc {
                        active_nnz -= 1;
                        n_rem += 1;
                    }
                    work += prow.len() + ui.len() - usize::from(had_pc);
                    std::mem::swap(ui, &mut mbuf);
                }
            }

            // ---- Column views of the same update, then retire the pivot
            // row from them and the pivot column entirely. ----
            for &(j, pv) in prow.iter() {
                let cj = &mut cols[j];
                if mults.len() * MERGE_RATIO < cj.len() {
                    for &(i, mlt) in mults.iter() {
                        let upd = mlt * pv;
                        let p = cj.partition_point(|&(r, _)| r < i);
                        if p < cj.len() && cj[p].0 == i {
                            let nv = cj[p].1 - upd;
                            if nv == 0.0 {
                                cj.remove(p);
                                touch_colmax(j, i, 0.0, colmax, colmax_row, dirty);
                            } else {
                                cj[p].1 = nv;
                                touch_colmax(j, i, nv.abs(), colmax, colmax_row, dirty);
                            }
                        } else {
                            cj.insert(p, (i, -upd));
                            touch_colmax(j, i, upd.abs(), colmax, colmax_row, dirty);
                        }
                    }
                    if let Some(t) = find_entry(cj, pr) {
                        cj.remove(t);
                    }
                } else {
                    mbuf.clear();
                    let mut a = 0usize;
                    for &(i, mlt) in mults.iter() {
                        while a < cj.len() && cj[a].0 < i {
                            if cj[a].0 != pr {
                                mbuf.push(cj[a]);
                            }
                            a += 1;
                        }
                        let upd = mlt * pv;
                        if a < cj.len() && cj[a].0 == i {
                            let nv = cj[a].1 - upd;
                            a += 1;
                            if nv == 0.0 {
                                touch_colmax(j, i, 0.0, colmax, colmax_row, dirty);
                            } else {
                                mbuf.push((i, nv));
                                touch_colmax(j, i, nv.abs(), colmax, colmax_row, dirty);
                            }
                        } else {
                            mbuf.push((i, -upd));
                            touch_colmax(j, i, upd.abs(), colmax, colmax_row, dirty);
                        }
                    }
                    for &e in &cj[a..] {
                        if e.0 != pr {
                            mbuf.push(e);
                        }
                    }
                    std::mem::swap(cj, &mut mbuf);
                }
                if colmax_row[j] == pr {
                    dirty[j] = true;
                }
                if !chg[j] {
                    chg[j] = true;
                    chg_list.push(j);
                }
            }
            if let Some(s) = fst.as_ref() {
                s.n_ins.set(s.n_ins.get() + n_ins);
                s.n_rem.set(s.n_rem.get() + n_rem);
            }
            cols[pc].clear();
            cl.remove(pc);
            active_col[pc] = false;

            // ---- Sync count buckets for exactly the lines modified this
            // step; colmax was maintained incrementally (`touch_colmax`)
            // and only holder-hit columns need the exact rescan. ----
            for &(i, _) in mults.iter() {
                rl.update(i, rows[i].len());
            }
            for &j in chg_list.iter() {
                if active_col[j] {
                    cl.update(j, cols[j].len());
                    if dirty[j] {
                        recalc_colmax(cols, j, colmax, colmax_row);
                        dirty[j] = false;
                    }
                }
            }
            for j in chg_list.iter() {
                chg[*j] = false;
            }
            chg_list.clear();
        }

        // ---- Final remap: ORIGINAL indices → positions (step indices). ----
        let at0 = std::time::Instant::now();
        if let (Some(s), Some(t)) = (fst.as_ref(), kt0) {
            s.kernel_ns
                .set(s.kernel_ns.get() + at0.duration_since(t).as_nanos() as f64);
        }
        row_inv.resize(m, 0);
        col_inv.resize(m, 0);
        for (p, &r) in row_perm.iter().enumerate() {
            row_inv[r] = p as u32;
        }
        for (p, &c) in col_perm.iter().enumerate() {
            col_inv[c] = p as u32;
        }
        for t in l_trips.iter_mut() {
            t.0 = row_inv[t.0 as usize];
        }
        for t in u_trips.iter_mut() {
            t.0 = row_inv[t.0 as usize];
            t.1 = col_inv[t.1 as usize];
        }
        l_crs.assign_triplets(m, l_trips, true);
        l_csc.assign_presorted(m, l_trips, false);
        u_crs.assign_presorted(m, u_trips, true);
        u_csc.assign_triplets(m, u_trips, false);
        // Solve temporaries: `scratch` is fully rewritten before each read;
        // a recycled `HyperScratch` is all-zero `y` at rest with its own
        // consistent mark/stamp pair.
        scratch.get_mut().resize(m, 0.0);
        {
            let h = hs.get_mut();
            if h.y.len() != m || h.mark.len() != m {
                *h = HyperScratch::new(m);
            }
            debug_assert!(h.y.iter().all(|&v| v == 0.0));
        }
        // Fixed per-factorization overhead (~30 ns per row/col measured on
        // m = 7260..19847 near-slack bases: flat peel arrays, remap), in
        // ~10 ns work units.
        work += l_trips.len() + u_trips.len() + 3 * m;
        if let (Some(s), Some(t0)) = (fst.as_ref(), ft0) {
            s.asm_ns
                .set(s.asm_ns.get() + at0.elapsed().as_nanos() as f64);
            s.fact_ns
                .set(s.fact_ns.get() + t0.elapsed().as_nanos() as f64);
            s.n_fact.set(s.n_fact.get() + 1);
            s.l_nnz.set(l_trips.len() as u64);
            s.u_nnz.set(u_trips.len() as u64);
        }

        Ok(Self {
            n: m,
            row_perm,
            col_perm,
            l_crs,
            l_csc,
            u_crs,
            u_csc,
            u_diag: u_diag_step,
            scratch,
            row_inv,
            col_inv,
            hs,
            hyper_skip: Default::default(),
            hyper_blew: Default::default(),
            build_work: work,
        })
    }

    pub fn n(&self) -> usize {
        self.n
    }

    /// Entries touched while factoring (elimination + factor nonzeros).
    pub fn build_work(&self) -> usize {
        self.build_work
    }

    /// Work of one FTRAN/BTRAN pass in the same units as [`build_work`]
    /// (factor nonzeros plus the O(n) sweeps).
    ///
    /// [`build_work`]: SparseLu::build_work
    pub fn solve_work(&self) -> usize {
        self.l_crs.nnz() + self.u_crs.nnz() + 3 * self.n
    }

    /// Hyper-path backoff gate for direction `dir` (0 = FTRAN, 1 =
    /// BTRAN): `true` = this solve goes straight to the dense sweep
    /// (consuming one unit of the backoff). See [`hskip_mode`] for the
    /// policies; every policy only chooses between two bit-identical
    /// paths, so none can move a trajectory.
    #[inline]
    fn hs_skip(&self, dir: usize) -> bool {
        let c = &self.hyper_skip[if hskip_mode() == 0 { 0 } else { dir }];
        let skip = c.get();
        if skip > 0 {
            c.set(skip - 1);
            return true;
        }
        false
    }

    /// Arm the backoff after a dense outcome: `blow` = the reach blew up
    /// mid-sweep (wasted partial work), else the rhs-count refusal (free
    /// to detect, so the backoff policies do not back off on it).
    #[inline]
    fn hs_backoff(&self, dir: usize, blow: bool) {
        match hskip_mode() {
            0 => self.hyper_skip[0].set(HYPER_BACKOFF),
            1 if blow => self.hyper_skip[dir].set(HYPER_BACKOFF),
            3 if blow => {
                if self.hyper_blew[dir].replace(true) {
                    self.hyper_skip[dir].set(HYPER_BACKOFF);
                }
            }
            _ => {}
        }
    }

    /// A hyper attempt in direction `dir` stayed sparse (mode 3 streak reset).
    #[inline]
    fn hs_ok(&self, dir: usize) {
        self.hyper_blew[dir].set(false);
    }

    /// Hyper-sparse entry shared by FTRAN/BTRAN: scatter the nonzeros of `b`
    /// through the inverse permutation, run two reach-limited sweeps and
    /// gather through `out_perm`. Returns `false` (nothing touched) when `b`
    /// is too dense or `n` too small to pay off; the caller then runs the
    /// dense sweeps. Bit-identical to the dense path (see [`sparse_sweep`]).
    #[allow(clippy::too_many_arguments)]
    fn solve_hyper(
        &self,
        b: &[f64],
        x: &mut [f64],
        inv: &[u32],
        first: (&Flat, Option<&[f64]>),
        second: (&Flat, Option<&[f64]>),
        out_perm: &[usize],
    ) -> bool {
        let n = self.n;
        let dir = usize::from(!std::ptr::eq(inv.as_ptr(), self.row_inv.as_ptr()));
        sv_bump(|s| s.sv_n[dir].set(s.sv_n[dir].get() + 1));
        if n < HYPER_MIN_N {
            sv_bump(|s| s.sv_dense[dir].set(s.sv_dense[dir].get() + 1));
            return false;
        }
        if self.hs_skip(dir) {
            sv_bump(|s| s.sv_skip[dir].set(s.sv_skip[dir].get() + 1));
            return false;
        }
        let thr = n / HYPER_DIV;
        if b.iter().filter(|v| **v != 0.0).count() > thr {
            self.hs_backoff(dir, false);
            sv_bump(|s| s.sv_cnt[dir].set(s.sv_cnt[dir].get() + 1));
            return false;
        }
        let mut hs = self.hs.borrow_mut();
        let hs = &mut *hs;
        hs.seeds.clear();
        for (j, &bj) in b.iter().enumerate() {
            if bj != 0.0 {
                let k = inv[j];
                hs.y[k as usize] = bj;
                hs.seeds.push(k);
            }
        }
        let dense1 = sparse_sweep(first.0, first.1, true, n, thr, hs);
        let dense2 = if dense1 {
            dense_sweep(second.0, second.1, false, &mut hs.y);
            true
        } else {
            std::mem::swap(&mut hs.seeds, &mut hs.out);
            sparse_sweep(second.0, second.1, false, n, thr, hs)
        };
        if dense2 {
            self.hs_backoff(dir, true);
            sv_bump(|s| s.sv_blow[dir].set(s.sv_blow[dir].get() + 1));
            // `out_perm` is a permutation: every entry of x is overwritten.
            for k in 0..n {
                x[out_perm[k]] = hs.y[k];
            }
            hs.y.fill(0.0);
        } else {
            self.hs_ok(dir);
            sv_bump(|s| {
                s.sv_ok[dir].set(s.sv_ok[dir].get() + 1);
                s.sv_out[dir].set(s.sv_out[dir].get() + hs.out.len() as u64);
            });
            x.fill(0.0);
            for &k in &hs.out {
                x[out_perm[k as usize]] = hs.y[k as usize];
                hs.y[k as usize] = 0.0;
            }
        }
        true
    }

    /// [`Self::solve_hyper`] driven by
    /// an index-carrying rhs. The entry no longer scans `b` for its nonzeros
    /// (count + scatter were two O(n) passes per solve) and the result is
    /// delivered with its pattern instead of an O(n) `x.fill(0.0)` plus
    /// scatter. Same thresholds and backoff as the dense entry; the sweeps
    /// are the same [`sparse_sweep`] calls, whose fp sequence is fixed by
    /// the heap order (independent of seed order), so values are
    /// bit-identical to [`Self::solve_hyper`] and to the dense solves. HiGHS
    /// `HFactor::ftranL/ftranU` hyper paths take the rhs `HVector`'s index
    /// the same way (util/HFactor.cpp, MIT, read as reference).
    #[allow(clippy::too_many_arguments)]
    fn solve_hyper_hv(
        &self,
        b: &HVec,
        x: &mut HVec,
        inv: &[u32],
        first: (&Flat, Option<&[f64]>),
        second: (&Flat, Option<&[f64]>),
        out_perm: &[usize],
    ) -> bool {
        let n = self.n;
        let dir = usize::from(!std::ptr::eq(inv.as_ptr(), self.row_inv.as_ptr()));
        sv_bump(|s| s.sv_n[dir].set(s.sv_n[dir].get() + 1));
        if n < HYPER_MIN_N {
            sv_bump(|s| s.sv_dense[dir].set(s.sv_dense[dir].get() + 1));
            return false;
        }
        if self.hs_skip(dir) {
            sv_bump(|s| s.sv_skip[dir].set(s.sv_skip[dir].get() + 1));
            return false;
        }
        let thr = n / HYPER_DIV;
        let bv = b.vals();
        let cnt = match b.pattern() {
            Some(p) => p.iter().filter(|&&j| bv[j as usize] != 0.0).count(),
            None => bv.iter().filter(|v| **v != 0.0).count(),
        };
        if cnt > thr {
            self.hs_backoff(dir, false);
            sv_bump(|s| s.sv_cnt[dir].set(s.sv_cnt[dir].get() + 1));
            return false;
        }
        let mut hs = self.hs.borrow_mut();
        let hs = &mut *hs;
        hs.seeds.clear();
        match b.pattern() {
            Some(p) => {
                for &j in p {
                    let bj = bv[j as usize];
                    if bj != 0.0 {
                        let k = inv[j as usize];
                        hs.y[k as usize] = bj;
                        hs.seeds.push(k);
                    }
                }
            }
            None => {
                for (j, &bj) in bv.iter().enumerate() {
                    if bj != 0.0 {
                        let k = inv[j];
                        hs.y[k as usize] = bj;
                        hs.seeds.push(k);
                    }
                }
            }
        }
        let dense1 = sparse_sweep(first.0, first.1, true, n, thr, hs);
        let dense2 = if dense1 {
            dense_sweep(second.0, second.1, false, &mut hs.y);
            true
        } else {
            std::mem::swap(&mut hs.seeds, &mut hs.out);
            sparse_sweep(second.0, second.1, false, n, thr, hs)
        };
        if dense2 {
            self.hs_backoff(dir, true);
            sv_bump(|s| s.sv_blow[dir].set(s.sv_blow[dir].get() + 1));
            let xv = x.dense_mut();
            for k in 0..n {
                xv[out_perm[k]] = hs.y[k];
            }
            hs.y.fill(0.0);
        } else {
            self.hs_ok(dir);
            sv_bump(|s| {
                s.sv_ok[dir].set(s.sv_ok[dir].get() + 1);
                s.sv_out[dir].set(s.sv_out[dir].get() + hs.out.len() as u64);
            });
            x.clear();
            for &k in &hs.out {
                x.set(out_perm[k as usize], hs.y[k as usize]);
                hs.y[k as usize] = 0.0;
            }
        }
        true
    }

    /// FTRAN with an index-carrying rhs/result (see
    /// [`Self::solve_hyper_hv`]); the dense fallback leaves `x` in dense
    /// mode. Bit-identical to [`Self::solve`].
    pub(crate) fn solve_hv(&self, b: &HVec, x: &mut HVec) {
        debug_assert_eq!(b.len(), self.n);
        debug_assert_eq!(x.len(), self.n);
        if self.solve_hyper_hv(
            b,
            x,
            &self.row_inv,
            (&self.l_csc, None),
            (&self.u_csc, Some(&self.u_diag)),
            &self.col_perm,
        ) {
            return;
        }
        self.solve_dense(b.vals(), x.dense_mut());
    }

    /// BTRAN twin of [`Self::solve_hv`]. Bit-identical to
    /// [`Self::solve_transpose`].
    pub(crate) fn solve_transpose_hv(&self, b: &HVec, x: &mut HVec) {
        debug_assert_eq!(b.len(), self.n);
        debug_assert_eq!(x.len(), self.n);
        if self.solve_hyper_hv(
            b,
            x,
            &self.col_inv,
            (&self.u_crs, Some(&self.u_diag)),
            (&self.l_crs, None),
            &self.row_perm,
        ) {
            return;
        }
        self.solve_transpose_dense(b.vals(), x.dense_mut());
    }

    /// FTRAN: solve `A x = b` (A the factored column basis).
    pub fn solve(&self, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        assert_eq!(b.len(), n, "solve: wrong rhs length");
        assert_eq!(x.len(), n, "solve: wrong solution length");
        if self.solve_hyper(
            b,
            x,
            &self.row_inv,
            (&self.l_csc, None),
            (&self.u_csc, Some(&self.u_diag)),
            &self.col_perm,
        ) {
            return;
        }
        self.solve_dense(b, x);
    }

    /// Dense-sweep FTRAN (reference path of [`SparseLu::solve`]).
    fn solve_dense(&self, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        let mut ybuf = self.scratch.borrow_mut();
        if ybuf.len() != n {
            ybuf.resize(n, 0.0);
        }
        let y = &mut *ybuf;
        for (k, yk) in y.iter_mut().enumerate() {
            *yk = b[self.row_perm[k]];
        }
        // L w = y (unit lower): column sweep, skipping zero
        // entries.
        for k in 0..n {
            let yk = y[k];
            if yk != 0.0 {
                for (i, v) in self.l_csc.line(k) {
                    y[i] -= v * yk;
                }
            }
        }
        // U z = w (upper): column sweep from the back.
        for k in (0..n).rev() {
            if y[k] != 0.0 {
                let zk = y[k] / self.u_diag[k];
                y[k] = zk;
                for (i, v) in self.u_csc.line(k) {
                    y[i] -= v * zk;
                }
            }
        }
        // x = Q z
        for k in 0..n {
            x[self.col_perm[k]] = y[k];
        }
    }

    /// BTRAN: solve `Aᵀ x = b`.
    pub fn solve_transpose(&self, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        assert_eq!(b.len(), n, "solve_transpose: wrong rhs length");
        assert_eq!(x.len(), n, "solve_transpose: wrong solution length");
        if self.solve_hyper(
            b,
            x,
            &self.col_inv,
            (&self.u_crs, Some(&self.u_diag)),
            (&self.l_crs, None),
            &self.row_perm,
        ) {
            return;
        }
        self.solve_transpose_dense(b, x);
    }

    /// Dense-sweep BTRAN (reference path of [`SparseLu::solve_transpose`]).
    fn solve_transpose_dense(&self, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        let mut wbuf = self.scratch.borrow_mut();
        if wbuf.len() != n {
            wbuf.resize(n, 0.0);
        }
        let w = &mut *wbuf;
        for (k, wk) in w.iter_mut().enumerate() {
            *wk = b[self.col_perm[k]];
        }
        // Uᵀ p = v (lower, forward): row sweep of U.
        for k in 0..n {
            if w[k] != 0.0 {
                let pk = w[k] / self.u_diag[k];
                w[k] = pk;
                for (j, v) in self.u_crs.line(k) {
                    w[j] -= v * pk;
                }
            }
        }
        // Lᵀ u = p (unit upper, back): row sweep of L.
        for k in (0..n).rev() {
            let wk = w[k];
            if wk != 0.0 {
                for (j, v) in self.l_crs.line(k) {
                    w[j] -= v * wk;
                }
            }
        }
        for k in 0..n {
            x[self.row_perm[k]] = w[k];
        }
    }

    // ------------------------------------------------------------------
    // Forrest–Tomlin support (see `lu_ft.rs`). The L-side
    // solves and the U-arena seed, split out of the fused `solve` /
    // `solve_transpose` so the FT update can own the U side of the
    // factorization while the base `L` stays frozen between refactor
    // rebuilds.
    // ------------------------------------------------------------------

    /// L-pass of FTRAN with the row permutation applied:
    /// `w[a] = b[row_perm[a]]`, then the forward `L⁻¹` column sweep.
    /// `w` is in permuted slot space; the U side is *not* solved here —
    /// the FT arena (`lu_ft.rs`) takes over.
    pub(crate) fn ftran_l(&self, b: &[f64], w: &mut [f64]) {
        let n = self.n;
        debug_assert_eq!(b.len(), n, "ftran_l: wrong rhs length");
        debug_assert_eq!(w.len(), n, "ftran_l: wrong work length");
        for (a, wa) in w.iter_mut().enumerate() {
            *wa = b[self.row_perm[a]];
        }
        for k in 0..n {
            let wk = w[k];
            if wk != 0.0 {
                for (i, v) in self.l_csc.line(k) {
                    w[i] -= v * wk;
                }
            }
        }
    }

    /// Hyper-sparse L-pass of FTRAN for the FT chain (`lu_ft.rs`): if the
    /// rhs qualifies (n ≥ `HYPER_MIN_N`, nonzero count ≤ n/`HYPER_DIV`,
    /// no backoff), scatter the rhs through `row_inv` into the persistent
    /// hyper workspace and run the reach-limited ascending L sweep
    /// ([`sparse_sweep`] over `l_csc` — bit-identical to [`ftran_l`]).
    ///
    /// Returns `(took, stayed_sparse)`:
    /// - `(false, _)`: rhs too dense / on backoff — nothing was written,
    ///   the caller runs the dense [`ftran_l`];
    /// - `(true, false)`: the sweep blew past the threshold and finished
    ///   densely — `w` holds the complete L result (no usable pattern),
    ///   the caller runs the dense R-eta + Û path;
    /// - `(true, true)`: `w` holds the sparse L result and `pattern` the
    ///   reached slot indices (ascending not guaranteed; `w` is zero
    ///   outside the pattern).
    pub(crate) fn ftran_l_hyper(
        &self,
        b: &[f64],
        w: &mut [f64],
        pattern: &mut Vec<u32>,
    ) -> (bool, bool) {
        let n = self.n;
        debug_assert_eq!(b.len(), n, "ftran_l_hyper: wrong rhs length");
        if n < HYPER_MIN_N {
            return (false, false);
        }
        if self.hs_skip(0) {
            return (false, false);
        }
        let thr = n / HYPER_DIV;
        if b.iter().filter(|v| **v != 0.0).count() > thr {
            self.hs_backoff(0, false);
            return (false, false);
        }
        let mut hsb = self.hs.borrow_mut();
        let hs = &mut *hsb;
        hs.seeds.clear();
        for (j, &bj) in b.iter().enumerate() {
            if bj != 0.0 {
                let k = self.row_inv[j];
                hs.y[k as usize] = bj;
                hs.seeds.push(k);
            }
        }
        if sparse_sweep(&self.l_csc, None, true, n, thr, hs) {
            // Dense finish: hand over the full result, no pattern, and
            // back off (the partial sweep cost more than it saved).
            self.hs_backoff(0, true);
            w.copy_from_slice(&hs.y);
            hs.y.fill(0.0);
            (true, false)
        } else {
            self.hs_ok(0);
            pattern.clear();
            pattern.extend_from_slice(&hs.out);
            for &k in &hs.out {
                w[k as usize] = hs.y[k as usize];
                hs.y[k as usize] = 0.0;
            }
            (true, true)
        }
    }

    /// L-pass of BTRAN: the backward `L⁻ᵀ` row sweep on the slot-space
    /// vector `v` (in place), then the scatter `x[row_perm[a]] = v[a]`.
    pub(crate) fn btran_lt(&self, v: &mut [f64], x: &mut [f64]) {
        let n = self.n;
        debug_assert_eq!(v.len(), n, "btran_lt: wrong work length");
        debug_assert_eq!(x.len(), n, "btran_lt: wrong solution length");
        for k in (0..n).rev() {
            let vk = v[k];
            if vk != 0.0 {
                for (j, val) in self.l_crs.line(k) {
                    v[j] -= val * vk;
                }
            }
        }
        for (a, &va) in v.iter().enumerate() {
            x[self.row_perm[a]] = va;
        }
    }

    /// Column permutation `Q`: `col_perm[k]` is the basis position at
    /// permuted column `k` (needed by the FT arena's delivery maps).
    pub(crate) fn col_perm(&self) -> &[usize] {
        &self.col_perm
    }

    /// Seed the FT U-arena (`lu_ft.rs`): the base `U` in column form with
    /// entries left in permuted slot space and the pivots carried in
    /// `u_pivot_value`. At build time logical position `p` *is* permuted
    /// column `p`, pivoting at slot `p` (identity maps); FT updates evolve
    /// the logical ↔ slot ↔ position maps from there.
    pub(crate) fn ft_seed(&self) -> FtLuSeed {
        let n = self.n;
        let mut u_start = Vec::with_capacity(n + 1);
        let mut u_last = Vec::with_capacity(n);
        let mut u_index = Vec::new();
        let mut u_value = Vec::new();
        for p in 0..n {
            u_start.push(u_index.len() as u32);
            for (r, v) in self.u_csc.line(p) {
                u_index.push(r as u32);
                u_value.push(v);
            }
            u_last.push(u_index.len() as u32);
        }
        u_start.push(u_index.len() as u32);
        FtLuSeed {
            u_start,
            u_last,
            u_index,
            u_value,
            u_pivot_value: self.u_diag.clone(),
        }
    }

    /// Expand to `(L, U)` dense matrices in the permuted order so tests can
    /// verify `P·A·Q = L·U` directly.
    pub fn expand(&self) -> (DenseMatrix, DenseMatrix) {
        let n = self.n;
        let mut l = DenseMatrix::identity(n);
        for i in 0..n {
            for (j, v) in self.l_crs.line(i) {
                l.set(i, j, v);
            }
        }
        let mut u = DenseMatrix::zeros(n);
        for i in 0..n {
            u.set(i, i, self.u_diag[i]);
            for (j, v) in self.u_crs.line(i) {
                u.set(i, j, v);
            }
        }
        (l, u)
    }

    /// Total nonzeros in L (strict) + U (incl. diagonal).
    pub fn nnz(&self) -> usize {
        self.l_crs.nnz() + self.u_crs.nnz() + self.n
    }
}

/// Binary search for `idx` in a sorted (idx, value) list.
fn find_entry(list: &[Entry], idx: usize) -> Option<usize> {
    let p = list.partition_point(|&(k, _)| k < idx);
    (p < list.len() && list[p].0 == idx).then_some(p)
}

/// Incremental exact `(colmax, holder)` maintenance for one touched column
/// entry, |value| → `new_abs` (0.0 = removed; an insert passes its |value|).
/// A full `recalc_colmax` rescan per touched column cost 234M scanned
/// entries on qap10 alone (10.2M calls × 23) — but only a shrunken/vanished
/// HOLDER can lower a column max, and raises/equal-arrivals update the max
/// in O(1). The rules are order-independent per step and provably equal a
/// full rescan afterwards (raise ⇒ max is the largest arrival, holder the
/// smallest row among equal-max arrivals — `recalc_colmax`'s strictly-
/// greater walk over ascending rows picks exactly that; no raise ⇒ the
/// incumbent holder, smallest among the old achievers, still achieves the
/// old max), so pivot selection cannot move.
fn touch_colmax(
    j: usize,
    i: usize,
    new_abs: f64,
    colmax: &mut [f64],
    holder: &mut [usize],
    dirty: &mut [bool],
) {
    if new_abs > colmax[j] {
        colmax[j] = new_abs;
        holder[j] = i;
    } else if new_abs == colmax[j] {
        if i < holder[j] {
            holder[j] = i;
        }
    } else if i == holder[j] {
        dirty[j] = true;
    }
}

/// Recompute max |value| and its holder row for column `j` over its
/// (active-row) column list. Exact: retired rows are stripped from every
/// column list at retirement, so whatever remains is the active submatrix.
fn recalc_colmax(cols: &[Vec<Entry>], j: usize, colmax: &mut [f64], holder: &mut [usize]) {
    if let Some(s) = lu_fact_stats() {
        s.n_recalc.set(s.n_recalc.get() + 1);
        s.colmax_scan
            .set(s.colmax_scan.get() + cols[j].len() as u64);
    }
    let mut best = 0.0f64;
    let mut bi = holder[j];
    for &(i, v) in &cols[j] {
        let v = v.abs();
        if v > best {
            best = v;
            bi = i;
        }
    }
    colmax[j] = best;
    holder[j] = bi;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::{DenseLu, DenseMatrix};
    use crate::rng::XorShift64;

    /// Deterministic random matrix with density `d` (values in [-10, 10]).
    fn random_matrix(m: usize, d: f64, rng: &mut XorShift64) -> DenseMatrix {
        let mut a = DenseMatrix::zeros(m);
        for i in 0..m {
            for j in 0..m {
                if rng.next_f64() < d {
                    a.set(i, j, rng.uniform(-10.0, 10.0));
                }
            }
        }
        a
    }

    fn to_columns(a: &DenseMatrix) -> Vec<Vec<f64>> {
        (0..a.n()).map(|j| a.col(j).to_vec()).collect()
    }

    fn max_rel_diff(x: &[f64], y: &[f64]) -> f64 {
        x.iter()
            .zip(y)
            .map(|(a, b)| (a - b).abs() / b.abs().max(1.0))
            .fold(0.0f64, f64::max)
    }

    /// FTRAN/BTRAN agree with DenseLu on sparse random bases.
    #[test]
    fn agrees_with_dense_lu() {
        let mut rng = XorShift64::new(0xC0FFEE);
        for &m in &[4usize, 12, 40] {
            for _attempt in 0..8 {
                let a = random_matrix(m, 0.3, &mut rng);
                let cols = to_columns(&a);
                let Ok(slu) = SparseLu::factor(&cols, 1e-10) else {
                    continue; // rare random singularity: resample
                };
                let dlu = DenseLu::factor(&a, 1e-10).expect("dense factor of same matrix");
                for _ in 0..3 {
                    let b: Vec<f64> = (0..m).map(|_| rng.uniform(-10.0, 10.0)).collect();
                    let (mut xs, mut xd, mut xt) = (vec![0.0; m], vec![0.0; m], vec![0.0; m]);
                    slu.solve(&b, &mut xs);
                    dlu.solve(&b, &mut xd);
                    assert!(
                        max_rel_diff(&xs, &xd) < 1e-8,
                        "m={m} FTRAN {xs:?} vs {xd:?}"
                    );
                    slu.solve_transpose(&b, &mut xt);
                    let mut xd_t = vec![0.0; m];
                    dlu.solve_transpose(&b, &mut xd_t);
                    assert!(max_rel_diff(&xt, &xd_t) < 1e-8, "m={m} BTRAN");
                }
                break;
            }
        }
    }

    /// P·A·Q = L·U reconstruction within fp noise.
    #[test]
    fn reconstruction_permutations() {
        let mut rng = XorShift64::new(0xBEEF);
        let m = 16;
        for _attempt in 0..8 {
            let a = random_matrix(m, 0.35, &mut rng);
            let cols = to_columns(&a);
            let Ok(slu) = SparseLu::factor(&cols, 1e-10) else {
                continue;
            };
            let (l, u) = slu.expand();
            // (LU)[i][j] must equal A[row_perm[i]][col_perm[j]]
            for i in 0..m {
                for j in 0..m {
                    let lu = {
                        let mut s = 0.0;
                        for k in 0..m {
                            s += l.get(i, k) * u.get(k, j);
                        }
                        s
                    };
                    let aij = a.get(slu.row_perm[i], slu.col_perm[j]);
                    assert!(
                        (lu - aij).abs() <= 1e-9 * aij.abs().max(1.0),
                        "recon ({} {}) {lu} vs {aij}",
                        slu.row_perm[i],
                        slu.col_perm[j]
                    );
                }
            }
            break;
        }
    }

    /// Duplicate columns ⇒ structurally singular.
    #[test]
    fn singular_detected() {
        let mut rng = XorShift64::new(42);
        let a = random_matrix(6, 0.5, &mut rng);
        let mut cols = to_columns(&a);
        cols[5] = cols[2].clone();
        let r = SparseLu::factor(&cols, 1e-10);
        assert!(matches!(r, Err(FactorError::Singular { .. })));
    }

    /// Boeing2-class dynamic range: A = D₁·R·D₂ with geometric diagonal
    /// scaling (entries spanning ~1e5, matching boeing2's coefficients).
    /// DenseLu partial pivoting achieves ~6e-13 on this matrix (measured);
    /// threshold+Markowitz must recover x₀ to 1e-6.
    #[test]
    fn badly_scaled_recovery() {
        let m = 12;
        let mut rng = XorShift64::new(7);
        let r = random_matrix(m, 0.6, &mut rng);
        let scales: Vec<f64> = (0..m)
            .map(|i| 10f64.powi(i as i32 * 4 / m as i32 - 2))
            .collect();
        // A[i][j] = s_i · R[i][j] / s_j
        let mut a = DenseMatrix::zeros(m);
        for i in 0..m {
            for j in 0..m {
                a.set(i, j, scales[i] * r.get(i, j) / scales[j]);
            }
        }
        let cols = to_columns(&a);
        let slu = SparseLu::factor(&cols, 1e-10).expect("badly-scaled factor");
        let x0: Vec<f64> = (0..m).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let mut b = vec![0.0; m];
        for (j, &xj) in x0.iter().enumerate() {
            for (i, bi) in b.iter_mut().enumerate() {
                *bi += a.get(i, j) * xj;
            }
        }
        let mut x = vec![0.0; m];
        slu.solve(&b, &mut x);
        let rel = x0
            .iter()
            .zip(&x)
            .map(|(e, g)| (g - e).abs() / e.abs().max(1.0))
            .fold(0.0f64, f64::max);
        assert!(rel < 1e-6, "badly-scaled recovery rel={rel}");
    }

    /// Basis-like structure (slack columns + a coupled kernel, shuffled):
    /// the singleton-first path of the count-bucket search must agree with
    /// DenseLu on FTRAN and BTRAN, and the sparse input path must agree with
    /// the dense wrapper.
    #[test]
    fn triangular_plus_kernel_agrees_with_dense_lu() {
        let mut rng = XorShift64::new(0xFACADE);
        for &(m, kernel) in &[(60usize, 8usize), (150, 20), (200, 0), (120, 50), (90, 90)] {
            let mut a = DenseMatrix::zeros(m);
            // Slack part: identity on rows/cols kernel..m, plus a few
            // sub-diagonal couplings so singleton chains must be peeled.
            for k in kernel..m {
                a.set(k, k, rng.uniform(1.0, 3.0));
                if k + 1 < m && rng.next_f64() < 0.4 {
                    a.set(k + 1, k, rng.uniform(-2.0, 2.0));
                }
            }
            // Dense-ish kernel block with border entries.
            for i in 0..kernel {
                for j in 0..kernel {
                    if rng.next_f64() < 0.5 || i == j {
                        a.set(
                            i,
                            j,
                            rng.uniform(-3.0, 3.0) + if i == j { 4.0 } else { 0.0 },
                        );
                    }
                }
                if kernel < m {
                    a.set(i, kernel + (i % (m - kernel)), rng.uniform(-1.0, 1.0));
                }
            }
            let cols = to_columns(&a);
            let sp: Vec<Vec<(usize, f64)>> = cols
                .iter()
                .map(|c| {
                    c.iter()
                        .enumerate()
                        .filter_map(|(i, &v)| (v != 0.0).then_some((i, v)))
                        .collect()
                })
                .collect();
            let slu = SparseLu::factor_sparse(&sp, 1e-10).expect("sparse factor");
            let dlu = DenseLu::factor(&a, 1e-10).expect("dense factor");
            let b: Vec<f64> = (0..m).map(|_| rng.uniform(-5.0, 5.0)).collect();
            let (mut xs, mut xd) = (vec![0.0; m], vec![0.0; m]);
            slu.solve(&b, &mut xs);
            dlu.solve(&b, &mut xd);
            assert!(max_rel_diff(&xs, &xd) < 1e-8, "m={m} FTRAN");
            slu.solve_transpose(&b, &mut xs);
            dlu.solve_transpose(&b, &mut xd);
            assert!(max_rel_diff(&xs, &xd) < 1e-8, "m={m} BTRAN");
            // Mostly triangular ⇒ factors stay near the matrix nnz.
            let a_nnz: usize = sp.iter().map(|c| c.len()).sum();
            assert!(
                slu.nnz() <= 3 * a_nnz + m,
                "m={m} fill {} vs {a_nnz}",
                slu.nnz()
            );
        }
    }

    #[test]
    fn hyper_sparse_solves_are_bit_identical_to_dense() {
        let mut rng = XorShift64::new(77);
        let n = 600;
        // Sparse, diagonally dominant-ish random matrix.
        let mut cols: Vec<Vec<Entry>> = vec![Vec::new(); n];
        for (j, col) in cols.iter_mut().enumerate() {
            col.push((j, 4.0 + rng.next_f64()));
            for _ in 0..3 {
                let i = (rng.next_u64() as usize) % n;
                if i != j && !col.iter().any(|&(r, _)| r == i) {
                    col.push((i, rng.next_f64() - 0.5));
                }
            }
            col.sort_by_key(|&(r, _)| r);
        }
        let lu = SparseLu::factor_sparse(&cols, 0.1).unwrap();
        for nnz in [1usize, 3, 10, 30, 200] {
            let mut b = vec![0.0; n];
            for _ in 0..nnz {
                b[(rng.next_u64() as usize) % n] = rng.next_f64() - 0.5;
            }
            let (mut x1, mut x2) = (vec![0.0; n], vec![0.0; n]);
            lu.solve(&b, &mut x1);
            lu.solve_dense(&b, &mut x2);
            assert!(x1.iter().zip(&x2).all(|(a, c)| a == c), "ftran nnz={nnz}");
            lu.solve_transpose(&b, &mut x1);
            lu.solve_transpose_dense(&b, &mut x2);
            assert!(x1.iter().zip(&x2).all(|(a, c)| a == c), "btran nnz={nnz}");
        }
    }

    /// Generator: near-slack basis with a Markowitz kernel. `kern` is
    /// the fraction of off-diagonal entries landing on another structural
    /// "home" row (the kernel grows with it; large values reach the dense
    /// finish).
    fn near_slack_cols(m: usize, kern: f64, seed: u64) -> Vec<Vec<Entry>> {
        let mut rng = XorShift64::new(seed);
        let is_home = |i: usize| i % 5 == 2;
        let homes: Vec<usize> = (0..m).filter(|&i| is_home(i)).collect();
        (0..m)
            .map(|i| {
                if !is_home(i) {
                    return vec![(i, rng.uniform(0.5, 2.0))];
                }
                let mut c = vec![(i, rng.uniform(2.0, 3.0))];
                for _ in 0..4 {
                    let r = if rng.next_f64() < kern {
                        homes[(rng.next_u64() % homes.len() as u64) as usize]
                    } else {
                        (rng.next_u64() % m as u64) as usize
                    };
                    if c.iter().all(|&(k, _)| k != r) {
                        c.push((r, rng.uniform(-0.5, 0.5)));
                    }
                }
                c.sort_by_key(|&(k, _)| k);
                c
            })
            .collect()
    }

    fn mirror(cols: &[Vec<Entry>]) -> Vec<Vec<(u32, f64)>> {
        let mut rows = vec![Vec::new(); cols.len()];
        for (j, c) in cols.iter().enumerate() {
            for &(i, v) in c {
                rows[i].push((j as u32, v));
            }
        }
        rows
    }

    fn assert_same_lu(a: &SparseLu, b: &SparseLu, tag: &str) {
        assert_eq!(a.n, b.n, "{tag} n");
        assert_eq!(a.row_perm, b.row_perm, "{tag} row_perm");
        assert_eq!(a.col_perm, b.col_perm, "{tag} col_perm");
        assert_eq!(a.row_inv, b.row_inv, "{tag} row_inv");
        assert_eq!(a.col_inv, b.col_inv, "{tag} col_inv");
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&a.u_diag), bits(&b.u_diag), "{tag} u_diag");
        for (fa, fb, w) in [
            (&a.l_crs, &b.l_crs, "l_crs"),
            (&a.l_csc, &b.l_csc, "l_csc"),
            (&a.u_crs, &b.u_crs, "u_crs"),
            (&a.u_csc, &b.u_csc, "u_csc"),
        ] {
            assert_eq!(fa.ptr, fb.ptr, "{tag} {w} ptr");
            assert_eq!(fa.idx, fb.idx, "{tag} {w} idx");
            assert_eq!(bits(&fa.val), bits(&fb.val), "{tag} {w} val");
        }
        assert_eq!(a.build_work, b.build_work, "{tag} build_work");
    }

    /// The persistent-workspace build (row mirror instead of the
    /// transpose, recycled buffers, O(kernel) reset — including after a
    /// Singular return and after the dense-kernel finish) is bit-identical
    /// to a cold `factor_sparse` on every basis of a sequence.
    #[test]
    fn factor_ws_bit_identical() {
        let m = 3000;
        let mut ws = FactorWs::default();
        let mut ws_t = FactorWs::default();
        for (step, &(kern, seed)) in [
            (0.02, 1u64),
            (0.3, 2),
            (0.9, 3),
            (0.0, 4),
            (0.5, 5),
            (0.05, 6),
        ]
        .iter()
        .enumerate()
        {
            let mut cols = near_slack_cols(m, kern, seed);
            if step == 3 {
                // Singular basis: two identical unit columns. Both builds
                // must fail; the workspace must come back clean.
                cols[1] = cols[0].clone();
                let rows = mirror(&cols);
                assert!(SparseLu::factor_sparse(&cols, 0.1).is_err());
                assert!(SparseLu::factor_sparse_ws(&cols, Some(&rows), 0.1, &mut ws).is_err());
                assert!(SparseLu::factor_sparse_ws(&cols, None, 0.1, &mut ws_t).is_err());
                continue;
            }
            let rows = mirror(&cols);
            let cold = SparseLu::factor_sparse(&cols, 0.1).unwrap();
            let warm = SparseLu::factor_sparse_ws(&cols, Some(&rows), 0.1, &mut ws).unwrap();
            let warm_t = SparseLu::factor_sparse_ws(&cols, None, 0.1, &mut ws_t).unwrap();
            assert_same_lu(&cold, &warm, &format!("mirror step {step}"));
            assert_same_lu(&cold, &warm_t, &format!("transpose step {step}"));
            ws.recycle(Box::new(warm));
            ws_t.recycle(Box::new(warm_t));
        }
    }

    /// Timing probe (not a gate): `cargo test --release -p rusolve-lp
    /// factor_ws_timing -- --ignored --nocapture` with RUSOLVE_DUAL_STATS=1.
    #[test]
    #[ignore]
    fn factor_ws_timing() {
        let m = 200_000;
        let cols = near_slack_cols(m, 0.002, 7);
        let rows = mirror(&cols);
        let mut ws = FactorWs::default();
        let reps: usize = std::env::var("P42_REPS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(50);
        let cold = std::env::var("P42_COLD").is_ok();
        let t = std::time::Instant::now();
        for _ in 0..reps {
            if cold {
                std::hint::black_box(SparseLu::factor_sparse(&cols, 0.1).unwrap());
                continue;
            }
            let lu = SparseLu::factor_sparse_ws(&cols, Some(&rows), 0.1, &mut ws).unwrap();
            ws.recycle(Box::new(lu));
        }
        eprintln!(
            "warm {:.2} ms/build",
            t.elapsed().as_secs_f64() * 1e3 / reps as f64
        );
        eprintln!("{}", lu_fact_line());
    }
}
