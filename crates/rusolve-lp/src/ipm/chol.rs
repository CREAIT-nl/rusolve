//! Supernodal sparse Cholesky of the IPM normal matrix `M = A·Θ·Aᵀ + δI`.
//!
//! Provenance:
//! * fill-reducing ordering: [`super::amd`] (Amestoy, Davis & Duff 1996);
//! * elimination tree with Liu's path compression, postorder, and row-subtree
//!   column counts: Liu, "The role of elimination trees in sparse
//!   factorization", SIMAX 11 (1990); Davis, *Direct Methods for Sparse
//!   Linear Systems* (SIAM 2006), ch. 4;
//! * fundamental supernodes: Liu, Ng & Peyton, "On finding supernodes for
//!   sparse matrix computations", SIMAX 14 (1993);
//! * relaxed supernode amalgamation (merge a child supernode into its parent
//!   when the explicit zeros it adds are few; any topological reordering of
//!   the elimination tree keeps the fill): Ashcraft & Grimes, "The influence
//!   of relaxed supernode partitions on the multifrontal method", ACM TOMS
//!   15(4) (1989); the selection rule (merge the child adding the fewest
//!   artificial nonzeros while that count is ≤ a cap, cap bisected so the
//!   artificial flops are 1–2 % of the total) and the merged-supernode
//!   permutation follow HiGHS HiPO `Analyse::relaxSnMaxNz` /
//!   `Analyse::afterRelaxSn` (MIT, read as reference text; reimplemented);
//! * left-looking supernodal numeric factorization with dense block updates
//!   into a buffer and relative-index scatter: Ng & Peyton, "Block sparse
//!   Cholesky algorithms on advanced uniprocessor computers", SIAM J. Sci.
//!   Comput. 14(5) (1993). We replace the dynamic linked lists of pending
//!   updates by static per-target update lists computed in the symbolic
//!   phase (same updates, fixed order) so that independent subtrees can be
//!   factored concurrently;
//! * dense kernels (packed GEMM, blocked panel Cholesky): [`super::kernels`];
//! * subtree-to-thread scheduling of the elimination tree (independent
//!   subtrees factor concurrently, the remaining top of the tree runs on one
//!   thread with parallel dense kernels): Geist & Ng, "Task scheduling for
//!   parallel sparse Cholesky factorization", Int. J. Parallel Prog. 18(4)
//!   (1989); Pothen & Sun, "A mapping algorithm for parallel sparse Cholesky
//!   factorization", SIAM J. Sci. Comput. 14(5) (1993) (proportional mapping);
//! * tiny-pivot handling for (near-)dependent rows — pivot replaced by a huge
//!   value and its sub-column zeroed so the solve returns 0 in that
//!   coordinate: S. J. Wright, "Modified Cholesky factorizations in
//!   interior-point algorithms for linear programming", SIAM J. Optim. 9(4)
//!   (1999); Y. Zhang, LIPSOL (1996);
//! * pivot modification for the sparse part of a dense-column split (the
//!   modified diagonal entries are corrected by a low-rank term in the Schur
//!   complement): E. D. Andersen, "A modified Schur complement method for
//!   handling dense columns in interior-point methods for linear
//!   programming", ACM TOMS 22(3) (1996).
//!
//! Symbolic analysis runs once per LP (`analyse`); every IPM iteration calls
//! `factor(theta, delta)` and then any number of `solve`s.

use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use rusolve_core::CscMatrix;

use super::amd::amd_order;
use super::kernels::{self, GemmWork, Pivot};

const NONE: u32 = u32::MAX;

/// A subtree's column range with its disjoint slices of `L` and the
/// original diagonal, handed to one factorization thread.
type SubtreeSlice<'a> = (usize, usize, &'a mut [f64], &'a mut [f64]);

/// Relative tiny-pivot threshold: `d ≤ PIVOT_REL · M_jj` means the pivot
/// is cancellation noise (Wright 1999).
const PIVOT_REL: f64 = 1e-14;

fn pivot_rel() -> f64 {
    static V: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("RUSOLVE_IPM_PIVREL")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(PIVOT_REL)
    })
}
/// Experimental lever: when > 0, tiny pivots are perturbed to
/// `pivot_pert() · M_jj` (dynamic regularization) instead of skipped.
fn pivot_pert() -> f64 {
    static V: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("RUSOLVE_IPM_PIVPERT")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0)
    })
}
/// Upper target for the artificial-flop ratio of relaxed amalgamation
/// (HiPO `kUpperRatioRelax` = 0.02, lower = half of it). `RUSOLVE_IPM_RELAX=0`
/// keeps the fundamental supernodes (no amalgamation).
fn relax_target() -> f64 {
    static V: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("RUSOLVE_IPM_RELAX")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.02)
    })
}
/// Factorization threads (`RUSOLVE_IPM_THREADS`, default 1: single-threaded
/// like every other rusolve path and the HiGHS `threads=1` baseline).
pub fn ipm_threads() -> usize {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("RUSOLVE_IPM_THREADS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .map(|t| {
                if t == 0 {
                    std::thread::available_parallelism().map_or(1, |n| n.get())
                } else {
                    t
                }
            })
            .unwrap_or(1)
            .clamp(1, 64)
    })
}
/// Replacement for a skipped pivot (Wright 1999 / LIPSOL use 1e128 on D;
/// on the Cholesky factor L_jj = sqrt(1e128)).
const HUGE_PIVOT: f64 = 1e64;
/// Column block of the dense panel factorization.
const PANEL_NB: usize = 64;
/// Below this many multiply-adds a descendant update uses the plain loop
/// (packing overhead dominates the GEMM kernel).
const GEMM_MIN_WORK: usize = 4096;
/// Minimum flops of one dense kernel call before it is split across threads.
const PAR_MIN_FLOPS: f64 = 2e7;

/// Analysis / factorization statistics (reported by the CLI log and the
/// eval note).
#[derive(Debug, Clone, Default)]
pub struct CholStats {
    pub m: usize,
    pub nnz_a: usize,
    /// Lower-triangle nonzeros of M (incl. diagonal).
    pub nnz_m: usize,
    /// nnz(L) incl. diagonal (true structural count, no artificial zeros).
    pub nnz_l: usize,
    /// Stored lower-trapezoid entries after relaxed amalgamation.
    pub nnz_l_stored: usize,
    pub supernodes: usize,
    /// Fundamental supernodes before amalgamation.
    pub supernodes_fund: usize,
    /// Σ cc_j²: the classic Cholesky flop estimate (true structure).
    pub flops: f64,
    /// The same estimate on the stored (relaxed) structure.
    pub flops_stored: f64,
    /// Columns of A with more than `dense_col_thresh` nonzeros.
    pub dense_cols: usize,
    pub dense_col_thresh: usize,
    /// Share of M's off-diagonal (pattern) contributions Σk_j² owed to the
    /// dense columns.
    pub dense_col_share: f64,
    /// Columns split out of M and recovered by the Schur complement
    /// (`super::normal`), and modified pivots in the last factorization.
    pub dense_split: usize,
    pub dense_mods: usize,
    /// Pivots replaced by `HUGE_PIVOT` in the last factorization.
    pub skipped_pivots: usize,
    /// Of those, pivots not in the structural-dependence set (Θ-induced
    /// cancellation); drives the adaptive regularization in the solver.
    pub cancel_pivots: usize,
    /// Rows flagged structurally dependent by [`NormalChol::factor_detect`].
    pub dependent: usize,
    /// The dense-column split was abandoned for the unsplit factor
    /// (`super::normal`, split health probe).
    pub unsplit: bool,
    pub threads: usize,
    /// Independent subtrees scheduled on worker threads.
    pub subtrees: usize,
    pub t_order: f64,
    pub t_symbolic: f64,
}

/// Error from the symbolic phase.
#[derive(Debug, Clone)]
pub enum CholError {
    /// The normal matrix pattern or the factor exceeds the memory budget.
    TooLarge(String),
}

/// Symbolic cost of factoring `A·Aᵀ` (pattern + AMD + column counts), no
/// numeric storage: used by the cost model that chooses the linear-algebra
/// path before any numeric factorization.
#[derive(Debug, Clone, Copy)]
pub struct SymbolicCost {
    pub nnz_m: usize,
    pub nnz_l: usize,
    pub flops: f64,
}

/// Full symmetric pattern of `A·Aᵀ` (no diagonal) as adjacency lists.
fn normal_pattern(a: &CscMatrix, max_entries: usize) -> Result<(Vec<usize>, Vec<u32>), CholError> {
    let m = a.nrow();
    let n = a.ncol();
    let ptr = a.col_ptr();
    let ri = a.row_idx();
    let mut rcount = vec![0usize; m + 1];
    for &r in ri {
        rcount[r as usize + 1] += 1;
    }
    for i in 0..m {
        rcount[i + 1] += rcount[i];
    }
    let row_ptr = rcount.clone();
    let mut fillp = rcount;
    let mut row_cols = vec![0u32; a.nnz()];
    for j in 0..n {
        for &r in &ri[ptr[j]..ptr[j + 1]] {
            let r = r as usize;
            row_cols[fillp[r]] = j as u32;
            fillp[r] += 1;
        }
    }
    let mut mark = vec![NONE; m];
    let mut adj_ptr = Vec::with_capacity(m + 1);
    adj_ptr.push(0usize);
    let mut adj: Vec<u32> = Vec::new();
    for i in 0..m {
        mark[i] = i as u32;
        for &j in &row_cols[row_ptr[i]..row_ptr[i + 1]] {
            let j = j as usize;
            for &r in &ri[ptr[j]..ptr[j + 1]] {
                if mark[r as usize] != i as u32 {
                    mark[r as usize] = i as u32;
                    adj.push(r as u32);
                }
            }
        }
        if adj.len() > max_entries {
            return Err(CholError::TooLarge(format!(
                "normal-matrix pattern exceeds {max_entries} entries at row {i}/{m}"
            )));
        }
        adj_ptr.push(adj.len());
    }
    Ok((adj_ptr, adj))
}

/// Elimination tree, postorder and column counts of the pattern under
/// the ordering `perm0`. Returns (perm = postordered elimination order,
/// pinv, parent in postorder numbering, column counts incl. diagonal).
fn etree_counts(
    m: usize,
    adj_ptr: &[usize],
    adj: &[u32],
    perm0: &[usize],
) -> (Vec<usize>, Vec<usize>, Vec<u32>, Vec<usize>) {
    let mut pinv0 = vec![0usize; m];
    for (k, &p) in perm0.iter().enumerate() {
        pinv0[p] = k;
    }
    // Elimination tree (Liu) of the permuted pattern.
    let mut parent = vec![NONE; m];
    let mut ancestor = vec![NONE; m];
    for k in 0..m {
        let orig = perm0[k];
        for &r in &adj[adj_ptr[orig]..adj_ptr[orig + 1]] {
            let mut i = pinv0[r as usize];
            if i >= k {
                continue;
            }
            while i != NONE as usize && i < k {
                let nx = ancestor[i];
                ancestor[i] = k as u32;
                if nx == NONE {
                    parent[i] = k as u32;
                    break;
                }
                i = nx as usize;
            }
        }
    }
    // Postorder (iterative DFS over children lists, children ascending).
    let mut child_head = vec![NONE; m];
    let mut sib = vec![NONE; m];
    for v in (0..m).rev() {
        let p = parent[v];
        if p != NONE {
            sib[v] = child_head[p as usize];
            child_head[p as usize] = v as u32;
        }
    }
    let mut post = Vec::with_capacity(m);
    let mut stack: Vec<u32> = Vec::new();
    let mut it_child = child_head;
    for root in 0..m {
        if parent[root] != NONE {
            continue;
        }
        stack.push(root as u32);
        while let Some(&top) = stack.last() {
            let t = top as usize;
            let c = it_child[t];
            if c != NONE {
                it_child[t] = sib[c as usize];
                stack.push(c);
            } else {
                post.push(t);
                stack.pop();
            }
        }
    }
    debug_assert_eq!(post.len(), m);
    let mut postinv = vec![0usize; m];
    for (k, &v) in post.iter().enumerate() {
        postinv[v] = k;
    }
    let perm: Vec<usize> = post.iter().map(|&v| perm0[v]).collect();
    let mut pinv = vec![0usize; m];
    for (k, &p) in perm.iter().enumerate() {
        pinv[p] = k;
    }
    let mut par = vec![NONE; m];
    for v in 0..m {
        if parent[v] != NONE {
            par[postinv[v]] = postinv[parent[v] as usize] as u32;
        }
    }
    // Column counts by row subtrees (cc includes the diagonal).
    let mut cc = vec![1usize; m];
    let mut vmark = vec![NONE; m];
    for k in 0..m {
        vmark[k] = k as u32;
        let orig = perm[k];
        for &r in &adj[adj_ptr[orig]..adj_ptr[orig + 1]] {
            let mut j = pinv[r as usize];
            while j < k && vmark[j] != k as u32 {
                vmark[j] = k as u32;
                cc[j] += 1;
                j = if par[j] == NONE { k } else { par[j] as usize };
            }
        }
    }
    (perm, pinv, par, cc)
}

/// Symbolic cost of the normal matrix of `a` (AMD ordering), or TooLarge.
pub fn symbolic_cost(a: &CscMatrix, max_entries: usize) -> Result<SymbolicCost, CholError> {
    let m = a.nrow();
    let (adj_ptr, adj) = normal_pattern(a, max_entries)?;
    let perm0 = amd_order(m, &adj_ptr, &adj);
    let (_, _, _, cc) = etree_counts(m, &adj_ptr, &adj, &perm0);
    Ok(SymbolicCost {
        nnz_m: adj.len() / 2 + m,
        nnz_l: cc.iter().sum(),
        flops: cc.iter().map(|&c| (c as f64) * (c as f64)).sum(),
    })
}

/// Dense flops of a supernode with `sz` columns and `cl` rows below them:
/// Σ_{i<sz} (sz + cl − i − 1)² (HiPO `Analyse` convention).
fn sn_ops(sz: usize, cl: usize) -> f64 {
    let s = |n: f64| n * (n + 1.0) * (2.0 * n + 1.0) / 6.0;
    let t = (sz + cl) as f64;
    s(t - 1.0) - if cl > 0 { s(cl as f64 - 1.0) } else { 0.0 }
}

/// One relaxation pass with artificial-nonzero cap `max_art`: for each
/// supernode (children first), repeatedly merge the child that creates the
/// fewest artificial nonzeros while that count stays ≤ `max_art`.
/// Returns (merged_into, merged sizes, ratio of artificial flops).
fn relax_pass(
    size0: &[usize],
    clique: &[usize],
    first_child: &[u32],
    next_child: &[u32],
    max_art: u64,
    ops_norelax: f64,
) -> (Vec<u32>, f64) {
    let ns = size0.len();
    let mut size = size0.to_vec();
    let mut fake = vec![0u64; ns];
    let mut merged = vec![NONE; ns];
    let mut fc = first_child.to_vec();
    let mut nx = next_child.to_vec();
    for sn in 0..ns {
        loop {
            let mut best = (u64::MAX, 0usize, NONE);
            let mut c = fc[sn];
            while c != NONE {
                let cu = c as usize;
                let rows_filled = (size[sn] + clique[sn]).saturating_sub(clique[cu]) as u64;
                let tot = rows_filled * size[cu] as u64 + fake[sn] + fake[cu];
                if tot < best.0 || (tot == best.0 && best.1 < size[cu]) {
                    best = (tot, size[cu], c);
                }
                c = nx[cu];
            }
            if best.2 == NONE || best.0 > max_art {
                break;
            }
            let ch = best.2;
            size[sn] += best.1;
            fake[sn] = best.0;
            merged[ch as usize] = sn as u32;
            // unlink ch from sn's child list
            if fc[sn] == ch {
                fc[sn] = nx[ch as usize];
            } else {
                let mut c = fc[sn];
                while nx[c as usize] != ch {
                    c = nx[c as usize];
                }
                nx[c as usize] = nx[ch as usize];
            }
        }
    }
    let mut ops = 0.0;
    for sn in 0..ns {
        if merged[sn] == NONE {
            ops += sn_ops(size[sn], clique[sn]);
        }
    }
    let art = (ops - ops_norelax).max(0.0);
    (merged, art / (art + ops_norelax).max(1.0))
}

/// Relaxed amalgamation (HiPO `relaxSnMaxNz` bisection on the cap so the
/// artificial-flop ratio lands in [target/2, target]). Returns merged_into
/// per fundamental supernode (NONE = root of its merged group).
fn relax_supernodes(sn_first: &[usize], cc: &[usize], sn_parent: &[u32], target: f64) -> Vec<u32> {
    let ns = sn_first.len() - 1;
    if target <= 0.0 || ns < 2 {
        return vec![NONE; ns];
    }
    let size: Vec<usize> = (0..ns).map(|s| sn_first[s + 1] - sn_first[s]).collect();
    let clique: Vec<usize> = (0..ns).map(|s| cc[sn_first[s]] - size[s]).collect();
    let ops_norelax: f64 = (0..ns).map(|s| sn_ops(size[s], clique[s])).sum();
    let mut first_child = vec![NONE; ns];
    let mut next_child = vec![NONE; ns];
    for s in (0..ns).rev() {
        let p = sn_parent[s];
        if p != NONE {
            next_child[s] = first_child[p as usize];
            first_child[p as usize] = s as u32;
        }
    }
    let lower = 0.5 * target;
    let mut cap: u64 = 256;
    let (mut below, mut above): (Option<u64>, Option<u64>) = (None, None);
    let mut best = (f64::INFINITY, cap);
    for _ in 0..20 {
        let (merged, ratio) =
            relax_pass(&size, &clique, &first_child, &next_child, cap, ops_norelax);
        let dist = (ratio - lower).abs().min((ratio - target).abs());
        if dist < best.0 {
            best = (dist, cap);
        }
        if ratio < lower {
            below = Some(cap);
            cap = match above {
                None => cap.saturating_mul(2),
                Some(a) => (cap + a) / 2,
            };
        } else if ratio > target {
            above = Some(cap);
            cap = match below {
                None => cap / 2,
                Some(b) => (b + cap) / 2,
            };
        } else {
            return merged;
        }
        if below.zip(above).is_some_and(|(b, a)| a <= b + 1) {
            break;
        }
    }
    relax_pass(
        &size,
        &clique,
        &first_child,
        &next_child,
        best.1,
        ops_norelax,
    )
    .0
}

/// Immutable symbolic structure (shared by factor threads).
struct Sym {
    m: usize,
    /// perm[k] = original row eliminated k-th; pinv inverse.
    perm: Vec<usize>,
    pinv: Vec<usize>,
    /// A by columns with PERMUTED row indices.
    a_ptr: Vec<usize>,
    a_prow: Vec<u32>,
    a_val: Vec<f64>,
    /// Aᵀ in permuted row order: row k = (column j, value).
    at_ptr: Vec<usize>,
    at_col: Vec<u32>,
    at_val: Vec<f64>,
    /// Supernode partition: columns sn_first[s]..sn_first[s+1].
    sn_first: Vec<usize>,
    /// Row structure of each supernode (sorted, own columns first).
    sn_rows_ptr: Vec<usize>,
    sn_rows: Vec<u32>,
    sn_val_ptr: Vec<usize>,
    snode_of: Vec<u32>,
    /// Static update lists: target s receives updates from (src, first row
    /// position in src) for upd_ptr[s]..upd_ptr[s+1].
    upd_ptr: Vec<usize>,
    upd_src: Vec<u32>,
    upd_pos: Vec<u32>,
    /// Subtree schedule for `threads > 1`: per worker, its subtree ranges
    /// [lo, hi) of supernodes; `top` = the remaining supernodes (ascending).
    plan_threads: usize,
    plan: Vec<Vec<(usize, usize)>>,
    top: Vec<usize>,
}

/// Per-thread numeric workspace.
#[derive(Default)]
struct Ws {
    relmap: Vec<u32>,
    w: Vec<f64>,
    wmark: Vec<u32>,
    touched: Vec<u32>,
    buf: Vec<f64>,
    gw: Vec<GemmWork>,
    skipped: usize,
    cancel: usize,
    modified: Vec<(u32, f64)>,
}

impl Ws {
    fn new(m: usize, gws: usize) -> Self {
        Self {
            relmap: vec![0; m],
            w: vec![0.0; m],
            wmark: vec![NONE; m],
            gw: (0..gws.max(1)).map(|_| GemmWork::default()).collect(),
            ..Default::default()
        }
    }
}

/// Numeric controls shared by factor threads.
struct Ctl<'a> {
    theta: &'a [f64],
    delta: f64,
    /// Reference diagonal (permuted order) enabling pivot modification.
    ref_diag: Option<&'a [f64]>,
    /// Tiny-pivot threshold relative to the assembled diagonal.
    pivrel: f64,
    /// Structurally dependent columns (permuted): always skipped.
    dep: Option<&'a [bool]>,
    mod_rel: f64,
    /// Per-workspace cap on modified pivots (beyond it, tiny pivots skip;
    /// see [`NormalChol::factor_modified`]).
    mod_cap: usize,
    deadline: Option<Instant>,
    abort: AtomicBool,
    /// Threads for the dense kernels of one supernode.
    intra: usize,
}

impl Ctl<'_> {
    fn past_deadline(&self) -> bool {
        if self.abort.load(Ordering::Relaxed) {
            return true;
        }
        if self.deadline.is_some_and(|d| Instant::now() >= d) {
            self.abort.store(true, Ordering::Relaxed);
            return true;
        }
        false
    }
}

/// Normal-equations Cholesky: symbolic structure + numeric storage.
pub struct NormalChol {
    sym: Sym,
    l: Vec<f64>,
    /// Assembled diagonal of M per permuted column (tiny-pivot reference).
    diag_orig: Vec<f64>,
    ws: Vec<Ws>,
    pub stats: CholStats,
    /// Wall-clock deadline checked per supernode, per descendant update and
    /// per panel block by `factor`; when it passes, `factor` returns early
    /// with `aborted = true` (the factor is then unusable).
    pub deadline: Option<Instant>,
    pub aborted: bool,
    /// Pivots modified by `factor_modified` in the last factorization:
    /// (permuted column, amount added to that diagonal entry of M).
    pub modified: Vec<(u32, f64)>,
    /// Structural-dependence flags per permuted column (empty = none).
    dep: Vec<bool>,
}

impl NormalChol {
    /// Symbolic analysis of `A·Aᵀ` for an `m × n` matrix `a`.
    /// `max_entries` caps nnz(pattern of M) and nnz(L) (memory guard).
    pub fn analyse(a: &CscMatrix, max_entries: usize) -> Result<Self, CholError> {
        Self::analyse_opts(a, max_entries, relax_target(), ipm_threads())
    }

    /// [`Self::analyse`] with explicit relaxation target and thread count.
    pub fn analyse_opts(
        a: &CscMatrix,
        max_entries: usize,
        relax: f64,
        threads: usize,
    ) -> Result<Self, CholError> {
        let t0 = Instant::now();
        let m = a.nrow();
        let n = a.ncol();
        let ptr = a.col_ptr();
        let ri = a.row_idx();
        let mut stats = CholStats {
            m,
            nnz_a: a.nnz(),
            threads,
            ..Default::default()
        };
        // Dense-column census (reported; the split is in `super::normal`).
        let thresh = 1000usize.max((10.0 * (m as f64).sqrt()) as usize);
        stats.dense_col_thresh = thresh;
        let (mut tot, mut dense_tot) = (0f64, 0f64);
        for j in 0..n {
            let k = (ptr[j + 1] - ptr[j]) as f64;
            tot += k * k;
            if ptr[j + 1] - ptr[j] > thresh {
                stats.dense_cols += 1;
                dense_tot += k * k;
            }
        }
        stats.dense_col_share = if tot > 0.0 { dense_tot / tot } else { 0.0 };

        let (adj_ptr, adj) = normal_pattern(a, max_entries)?;
        stats.nnz_m = adj.len() / 2 + m;
        let perm0 = amd_order(m, &adj_ptr, &adj);
        stats.t_order = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let (perm_post, _pinv_post, par, cc) = etree_counts(m, &adj_ptr, &adj, &perm0);
        drop(perm0);
        let nnz_l: usize = cc.iter().sum();
        if nnz_l > max_entries {
            return Err(CholError::TooLarge(format!(
                "nnz(L) = {nnz_l} exceeds {max_entries}"
            )));
        }
        stats.nnz_l = nnz_l;
        stats.flops = cc.iter().map(|&c| (c as f64) * (c as f64)).sum();

        // Fundamental supernodes (postorder numbering).
        let mut nchild = vec![0u32; m];
        for v in 0..m {
            if par[v] != NONE {
                nchild[par[v] as usize] += 1;
            }
        }
        let mut fsn_first = vec![0usize];
        for j in 1..m {
            let cont = par[j - 1] == j as u32 && cc[j - 1] == cc[j] + 1 && nchild[j] == 1;
            if !cont {
                fsn_first.push(j);
            }
        }
        if m > 0 {
            fsn_first.push(m);
        } else {
            fsn_first.clear();
            fsn_first.push(0);
        }
        let fns = fsn_first.len() - 1;
        stats.supernodes_fund = fns;
        let mut fsnode_of = vec![0u32; m];
        for s in 0..fns {
            for c in fsn_first[s]..fsn_first[s + 1] {
                fsnode_of[c] = s as u32;
            }
        }
        // Fundamental row structures = structure of each first column.
        let mut frows_ptr = vec![0usize; fns + 1];
        for s in 0..fns {
            frows_ptr[s + 1] = frows_ptr[s] + cc[fsn_first[s]];
        }
        let mut frows = vec![0u32; frows_ptr[fns]];
        let mut fillpos: Vec<usize> = frows_ptr[..fns].to_vec();
        for s in 0..fns {
            frows[fillpos[s]] = fsn_first[s] as u32;
            fillpos[s] += 1;
        }
        {
            let pinv_post = {
                let mut p = vec![0usize; m];
                for (k, &o) in perm_post.iter().enumerate() {
                    p[o] = k;
                }
                p
            };
            let mut vmark = vec![NONE; m];
            for k in 0..m {
                vmark[k] = k as u32;
                let orig = perm_post[k];
                for &r in &adj[adj_ptr[orig]..adj_ptr[orig + 1]] {
                    let mut j = pinv_post[r as usize];
                    while j < k && vmark[j] != k as u32 {
                        vmark[j] = k as u32;
                        let s = fsnode_of[j] as usize;
                        if fsn_first[s] == j {
                            frows[fillpos[s]] = k as u32;
                            fillpos[s] += 1;
                        }
                        j = if par[j] == NONE { k } else { par[j] as usize };
                    }
                }
            }
        }
        drop(adj);
        drop(adj_ptr);
        debug_assert!((0..fns).all(|s| fillpos[s] == frows_ptr[s + 1]));
        let fsn_parent: Vec<u32> = (0..fns)
            .map(|s| {
                let last = fsn_first[s + 1] - 1;
                if par[last] == NONE {
                    NONE
                } else {
                    fsnode_of[par[last] as usize]
                }
            })
            .collect();

        // Relaxed amalgamation and the induced (topological) renumbering.
        let merged = relax_supernodes(&fsn_first, &cc, &fsn_parent, relax);
        let mut received: Vec<Vec<u32>> = vec![Vec::new(); fns];
        let mut newpos = vec![0u32; m];
        let mut sn_first = vec![0usize];
        let mut group_root: Vec<usize> = Vec::new();
        let mut next_col = 0usize;
        let mut stack: Vec<u32> = Vec::new();
        for s in 0..fns {
            if merged[s] != NONE {
                received[merged[s] as usize].push(s as u32);
                continue;
            }
            stack.push(s as u32);
            while let Some(&cur) = stack.last() {
                let cu = cur as usize;
                if !received[cu].is_empty() {
                    let r = std::mem::take(&mut received[cu]);
                    stack.extend(r);
                } else {
                    stack.pop();
                    for c in fsn_first[cu]..fsn_first[cu + 1] {
                        newpos[c] = next_col as u32;
                        next_col += 1;
                    }
                }
            }
            sn_first.push(next_col);
            group_root.push(s);
        }
        debug_assert_eq!(next_col, m);
        let ns = group_root.len();
        stats.supernodes = ns;
        let mut perm = vec![0usize; m];
        for c in 0..m {
            perm[newpos[c] as usize] = perm_post[c];
        }
        drop(perm_post);
        let mut pinv = vec![0usize; m];
        for (k, &p) in perm.iter().enumerate() {
            pinv[p] = k;
        }
        let mut snode_of = vec![0u32; m];
        for s in 0..ns {
            for c in sn_first[s]..sn_first[s + 1] {
                snode_of[c] = s as u32;
            }
        }
        // Rows of a merged supernode: its own columns, then the root's
        // clique (rows of the root below its own columns), renumbered.
        let mut sn_rows_ptr = vec![0usize; ns + 1];
        for s in 0..ns {
            let r = group_root[s];
            let nc_root = fsn_first[r + 1] - fsn_first[r];
            let clique = frows_ptr[r + 1] - frows_ptr[r] - nc_root;
            sn_rows_ptr[s + 1] = sn_rows_ptr[s] + (sn_first[s + 1] - sn_first[s]) + clique;
        }
        let mut sn_rows = vec![0u32; sn_rows_ptr[ns]];
        for s in 0..ns {
            let r = group_root[s];
            let nc_root = fsn_first[r + 1] - fsn_first[r];
            let dst = &mut sn_rows[sn_rows_ptr[s]..sn_rows_ptr[s + 1]];
            let nc = sn_first[s + 1] - sn_first[s];
            for (t, d) in dst[..nc].iter_mut().enumerate() {
                *d = (sn_first[s] + t) as u32;
            }
            for (d, &o) in dst[nc..]
                .iter_mut()
                .zip(&frows[frows_ptr[r] + nc_root..frows_ptr[r + 1]])
            {
                *d = newpos[o as usize];
            }
            dst[nc..].sort_unstable();
            debug_assert!(dst[nc..].iter().all(|&v| v as usize >= sn_first[s + 1]));
        }
        drop(frows);
        let mut sn_val_ptr = vec![0usize; ns + 1];
        let (mut stored, mut fl_stored) = (0usize, 0f64);
        for s in 0..ns {
            let nc = sn_first[s + 1] - sn_first[s];
            let nr = sn_rows_ptr[s + 1] - sn_rows_ptr[s];
            sn_val_ptr[s + 1] = sn_val_ptr[s] + nc * nr;
            stored += nc * nr - nc * (nc - 1) / 2;
            fl_stored += sn_ops(nc, nr - nc) + nr as f64 * 0.0;
        }
        stats.nnz_l_stored = stored;
        stats.flops_stored = fl_stored;
        if sn_val_ptr[ns] > max_entries {
            return Err(CholError::TooLarge(format!(
                "supernodal storage {} exceeds {max_entries}",
                sn_val_ptr[ns]
            )));
        }

        // Static update lists: every supernode k updates each target owning a
        // run of its off-diagonal rows.
        let mut ucount = vec![0usize; ns + 1];
        let mut runs: Vec<(u32, u32, u32)> = Vec::new();
        for k in 0..ns {
            let rl = sn_rows_ptr[k];
            let nr = sn_rows_ptr[k + 1] - rl;
            let nc = sn_first[k + 1] - sn_first[k];
            let mut q = nc;
            while q < nr {
                let t = snode_of[sn_rows[rl + q] as usize];
                runs.push((t, k as u32, q as u32));
                ucount[t as usize + 1] += 1;
                let lt = sn_first[t as usize + 1];
                while q < nr && (sn_rows[rl + q] as usize) < lt {
                    q += 1;
                }
            }
        }
        for s in 0..ns {
            ucount[s + 1] += ucount[s];
        }
        let upd_ptr = ucount.clone();
        let mut fp = ucount;
        let mut upd_src = vec![0u32; runs.len()];
        let mut upd_pos = vec![0u32; runs.len()];
        for &(t, k, q) in &runs {
            let t = t as usize;
            upd_src[fp[t]] = k;
            upd_pos[fp[t]] = q;
            fp[t] += 1;
        }
        drop(runs);

        // A with permuted row indices, and Aᵀ in permuted row order.
        let a_ptr = ptr.to_vec();
        let a_prow: Vec<u32> = ri.iter().map(|&r| pinv[r as usize] as u32).collect();
        let a_val = a.values().to_vec();
        let mut at_ptr = vec![0usize; m + 1];
        for &r in &a_prow {
            at_ptr[r as usize + 1] += 1;
        }
        for i in 0..m {
            at_ptr[i + 1] += at_ptr[i];
        }
        let mut fp = at_ptr.clone();
        let mut at_col = vec![0u32; a.nnz()];
        let mut at_val = vec![0f64; a.nnz()];
        for j in 0..n {
            for k in a_ptr[j]..a_ptr[j + 1] {
                let r = a_prow[k] as usize;
                at_col[fp[r]] = j as u32;
                at_val[fp[r]] = a_val[k];
                fp[r] += 1;
            }
        }
        let nval = sn_val_ptr[ns];
        let mut sym = Sym {
            m,
            perm,
            pinv,
            a_ptr,
            a_prow,
            a_val,
            at_ptr,
            at_col,
            at_val,
            sn_first,
            sn_rows_ptr,
            sn_rows,
            sn_val_ptr,
            snode_of,
            upd_ptr,
            upd_src,
            upd_pos,
            plan_threads: 1,
            plan: Vec::new(),
            top: Vec::new(),
        };
        if threads > 1 {
            sym.build_plan(threads);
            stats.subtrees = sym.plan.iter().map(|p| p.len()).sum();
        }
        stats.t_symbolic = t1.elapsed().as_secs_f64();
        let nws = if threads > 1 { threads } else { 1 };
        Ok(Self {
            ws: (0..nws).map(|_| Ws::new(m, threads)).collect(),
            sym,
            l: vec![0.0; nval],
            diag_orig: vec![0.0; m],
            stats,
            deadline: None,
            aborted: false,
            modified: Vec::new(),
            dep: Vec::new(),
        })
    }

    pub fn dim(&self) -> usize {
        self.sym.m
    }

    /// Numeric factorization of `A·diag(theta)·Aᵀ + δI`.
    pub fn factor(&mut self, theta: &[f64], delta: f64) {
        self.factor_impl(theta, delta, None, 0.0, 0, pivot_rel());
    }

    /// Factor with the looser tiny-pivot threshold `rel` and record the
    /// skipped pivots as *structurally dependent* rows: called once with
    /// Θ = I (the starting-point factor), where no Θ-induced cancellation
    /// exists, so a pivot ≤ rel·M_jj means a (near-)dependent row of A.
    /// Later factorizations always skip those rows (Wright 1999) and count
    /// every other tiny pivot as Θ-cancellation (`stats.cancel_pivots`).
    /// Returns the number of dependent rows.
    pub fn factor_detect(&mut self, theta: &[f64], delta: f64, rel: f64) -> usize {
        self.dep.clear();
        self.factor_impl(theta, delta, None, 0.0, 0, rel.max(pivot_rel()));
        if self.aborted {
            return 0;
        }
        self.dep = (0..self.sym.m)
            .map(|k| self.l_diag(k) == HUGE_PIVOT)
            .collect();
        self.stats.dependent = self.dep.iter().filter(|&&d| d).count();
        self.stats.cancel_pivots = 0;
        self.stats.dependent
    }

    /// Like [`Self::factor`], but a pivot `d ≤ max(PIVREL·M_jj, mod_rel·ref_j)`
    /// (with `ref_diag` in original row order) is replaced by `ref_j` instead
    /// of skipped; the factor is then of `M + Σ c_j e_j e_jᵀ` with the
    /// (permuted column, c_j) pairs in [`Self::modified`] (Andersen 1996).
    /// Every workspace may modify up to `mod_cap + 1` pivots (further tiny
    /// pivots are skipped), so `modified.len() > mod_cap` exactly when the
    /// serial factorization would also exceed the cap; the caller then
    /// refactors without modifications. This keeps the result independent of
    /// the thread count (a per-thread share of the cap did not: dano3_3 at 8
    /// threads took a different path).
    pub fn factor_modified(
        &mut self,
        theta: &[f64],
        delta: f64,
        ref_diag: &[f64],
        mod_rel: f64,
        mod_cap: usize,
    ) {
        let rp: Vec<f64> = self.sym.perm.iter().map(|&o| ref_diag[o]).collect();
        self.factor_impl(theta, delta, Some(&rp), mod_rel, mod_cap + 1, pivot_rel());
    }

    fn factor_impl(
        &mut self,
        theta: &[f64],
        delta: f64,
        ref_diag: Option<&[f64]>,
        mod_rel: f64,
        mod_cap: usize,
        pivrel: f64,
    ) {
        self.stats.skipped_pivots = 0;
        self.stats.cancel_pivots = 0;
        self.aborted = false;
        self.modified.clear();
        for ws in self.ws.iter_mut() {
            ws.wmark.iter_mut().for_each(|v| *v = NONE);
            ws.skipped = 0;
            ws.cancel = 0;
            ws.modified.clear();
        }
        let ctl = Ctl {
            theta,
            delta,
            ref_diag,
            pivrel,
            dep: if self.dep.is_empty() {
                None
            } else {
                Some(&self.dep)
            },
            mod_rel,
            mod_cap,
            deadline: self.deadline,
            abort: AtomicBool::new(false),
            intra: self.stats.threads.max(1),
        };
        let sym = &self.sym;
        let ns = sym.sn_first.len() - 1;
        if sym.plan_threads > 1 && !sym.plan.is_empty() {
            // Independent subtrees on worker threads (disjoint value ranges).
            let mut ranges: Vec<(usize, usize, usize)> = Vec::new();
            for (t, p) in sym.plan.iter().enumerate() {
                for &(lo, hi) in p {
                    ranges.push((lo, hi, t));
                }
            }
            ranges.sort_unstable();
            let mut per_thread: Vec<Vec<SubtreeSlice<'_>>> =
                (0..sym.plan.len()).map(|_| Vec::new()).collect();
            let mut lrest: &mut [f64] = &mut self.l;
            let mut drest: &mut [f64] = &mut self.diag_orig;
            let (mut loff, mut doff) = (0usize, 0usize);
            for &(lo, hi, t) in &ranges {
                let (vlo, vhi) = (sym.sn_val_ptr[lo], sym.sn_val_ptr[hi]);
                let (clo, chi) = (sym.sn_first[lo], sym.sn_first[hi]);
                let tmp = std::mem::take(&mut lrest);
                let (_, tail) = tmp.split_at_mut(vlo - loff);
                let (mine, tail) = tail.split_at_mut(vhi - vlo);
                lrest = tail;
                loff = vhi;
                let tmp = std::mem::take(&mut drest);
                let (_, dtail) = tmp.split_at_mut(clo - doff);
                let (dmine, dtail) = dtail.split_at_mut(chi - clo);
                drest = dtail;
                doff = chi;
                per_thread[t].push((lo, hi, mine, dmine));
            }
            let ctl_ref = &ctl;
            std::thread::scope(|sc| {
                for (jobs, ws) in per_thread.into_iter().zip(self.ws.iter_mut()) {
                    sc.spawn(move || {
                        for (lo, hi, lpart, dpart) in jobs {
                            let base = sym.sn_val_ptr[lo];
                            let cbase = sym.sn_first[lo];
                            for s in lo..hi {
                                if ctl_ref.past_deadline() {
                                    return;
                                }
                                let (desc, rest) = lpart.split_at_mut(sym.sn_val_ptr[s] - base);
                                let sz = sym.sn_val_ptr[s + 1] - sym.sn_val_ptr[s];
                                let c0 = sym.sn_first[s] - cbase;
                                let c1 = sym.sn_first[s + 1] - cbase;
                                if !sym.factor_sn(
                                    s,
                                    desc,
                                    base,
                                    &mut rest[..sz],
                                    &mut dpart[c0..c1],
                                    ctl_ref,
                                    ws,
                                    1,
                                ) {
                                    return;
                                }
                            }
                        }
                    });
                }
            });
            if ctl.abort.load(Ordering::Relaxed) {
                self.aborted = true;
                return;
            }
            let ws0 = &mut self.ws[0];
            for &s in &sym.top {
                if ctl.past_deadline() {
                    self.aborted = true;
                    return;
                }
                let (desc, rest) = self.l.split_at_mut(sym.sn_val_ptr[s]);
                let sz = sym.sn_val_ptr[s + 1] - sym.sn_val_ptr[s];
                let (c0, c1) = (sym.sn_first[s], sym.sn_first[s + 1]);
                if !sym.factor_sn(
                    s,
                    desc,
                    0,
                    &mut rest[..sz],
                    &mut self.diag_orig[c0..c1],
                    &ctl,
                    ws0,
                    ctl.intra,
                ) {
                    self.aborted = true;
                    return;
                }
            }
        } else {
            let ws0 = &mut self.ws[0];
            for s in 0..ns {
                if ctl.past_deadline() {
                    self.aborted = true;
                    return;
                }
                let (desc, rest) = self.l.split_at_mut(sym.sn_val_ptr[s]);
                let sz = sym.sn_val_ptr[s + 1] - sym.sn_val_ptr[s];
                let (c0, c1) = (sym.sn_first[s], sym.sn_first[s + 1]);
                if !sym.factor_sn(
                    s,
                    desc,
                    0,
                    &mut rest[..sz],
                    &mut self.diag_orig[c0..c1],
                    &ctl,
                    ws0,
                    1,
                ) {
                    self.aborted = true;
                    return;
                }
            }
        }
        for ws in self.ws.iter_mut() {
            self.stats.skipped_pivots += ws.skipped;
            self.stats.cancel_pivots += ws.cancel;
            self.modified.append(&mut ws.modified);
        }
    }

    /// Solve `M x = rhs` in place (original row order).
    pub fn solve(&self, rhs: &mut [f64], work: &mut Vec<f64>) {
        let sym = &self.sym;
        let m = sym.m;
        work.resize(m, 0.0);
        for k in 0..m {
            work[k] = rhs[sym.perm[k]];
        }
        self.solve_permuted(work);
        for k in 0..m {
            rhs[sym.perm[k]] = work[k];
        }
    }

    /// Solve `(L·Lᵀ) x = b` in the permuted order, in place.
    pub fn solve_permuted(&self, work: &mut [f64]) {
        self.forward(work);
        self.backward(work);
    }

    /// Forward solve `L z = b` in place (permuted order).
    pub fn forward(&self, work: &mut [f64]) {
        for s in 0..self.sym.sn_first.len() - 1 {
            self.forward_sn(s, work);
        }
    }

    fn forward_sn(&self, s: usize, work: &mut [f64]) {
        let sym = &self.sym;
        let f = sym.sn_first[s];
        let nc = sym.sn_first[s + 1] - f;
        let rl = sym.sn_rows_ptr[s];
        let nr = sym.sn_rows_ptr[s + 1] - rl;
        let rows = &sym.sn_rows[rl..rl + nr];
        let blk = &self.l[sym.sn_val_ptr[s]..sym.sn_val_ptr[s] + nr * nc];
        for j in 0..nc {
            let col = &blk[j * nr..(j + 1) * nr];
            let xj = work[f + j] / col[j];
            work[f + j] = xj;
            if xj != 0.0 {
                for r in j + 1..nr {
                    work[rows[r] as usize] -= col[r] * xj;
                }
            }
        }
    }

    /// Backward solve `Lᵀ x = z` in place (permuted order).
    pub fn backward(&self, work: &mut [f64]) {
        let sym = &self.sym;
        let ns = sym.sn_first.len() - 1;
        for s in (0..ns).rev() {
            let f = sym.sn_first[s];
            let nc = sym.sn_first[s + 1] - f;
            let rl = sym.sn_rows_ptr[s];
            let nr = sym.sn_rows_ptr[s + 1] - rl;
            let rows = &sym.sn_rows[rl..rl + nr];
            let blk = &self.l[sym.sn_val_ptr[s]..sym.sn_val_ptr[s] + nr * nc];
            for j in (0..nc).rev() {
                let col = &blk[j * nr..(j + 1) * nr];
                let mut acc = work[f + j];
                for r in j + 1..nr {
                    acc -= col[r] * work[rows[r] as usize];
                }
                work[f + j] = acc / col[j];
            }
        }
    }

    /// Number of supernodes.
    pub fn n_supernodes(&self) -> usize {
        self.sym.sn_first.len() - 1
    }

    /// Columns (permuted positions) of supernode `s`.
    pub fn sn_cols(&self, s: usize) -> std::ops::Range<usize> {
        self.sym.sn_first[s]..self.sym.sn_first[s + 1]
    }

    /// Sparse forward solve `L z = b` (Gilbert & Peierls, "Sparse partial
    /// pivoting in time proportional to arithmetic operations", SISC 9(5)
    /// (1988); Davis 2006 §3.2, at supernode granularity): `work` holds b
    /// (permuted, zero off `nz`), only the supernodes on the elimination-tree
    /// paths from `nz` to the roots are processed, and their ids are left in
    /// `reach` in ascending (topological) order; z is zero outside their
    /// columns. `mark` (one flag per supernode, all false) is restored.
    pub fn forward_sparse(
        &self,
        work: &mut [f64],
        nz: &[u32],
        mark: &mut [bool],
        reach: &mut Vec<u32>,
    ) {
        let sym = &self.sym;
        reach.clear();
        for &i in nz {
            let mut s = sym.snode_of[i as usize];
            while s != NONE && !mark[s as usize] {
                mark[s as usize] = true;
                reach.push(s);
                let su = s as usize;
                let rl = sym.sn_rows_ptr[su];
                let nr = sym.sn_rows_ptr[su + 1] - rl;
                let nc = sym.sn_first[su + 1] - sym.sn_first[su];
                s = if nr > nc {
                    sym.snode_of[sym.sn_rows[rl + nc] as usize]
                } else {
                    NONE
                };
            }
        }
        reach.sort_unstable();
        for &s in reach.iter() {
            mark[s as usize] = false;
            self.forward_sn(s as usize, work);
        }
    }

    /// `y = (A·diag(theta)·Aᵀ + δI)·x` (original row order), for refinement.
    pub fn matvec(&self, theta: &[f64], delta: f64, x: &[f64], y: &mut [f64]) {
        let sym = &self.sym;
        let n = sym.a_ptr.len() - 1;
        for (k, yv) in y.iter_mut().enumerate() {
            *yv = delta * x[k];
        }
        for j in 0..n {
            let mut s = 0.0;
            for q in sym.a_ptr[j]..sym.a_ptr[j + 1] {
                s += sym.a_val[q] * x[sym.perm[sym.a_prow[q] as usize]];
            }
            if s == 0.0 {
                continue;
            }
            let t = s * theta[j];
            for q in sym.a_ptr[j]..sym.a_ptr[j + 1] {
                y[sym.perm[sym.a_prow[q] as usize]] += sym.a_val[q] * t;
            }
        }
    }

    /// Position of original row `r` in the elimination order.
    pub fn pinv(&self, r: usize) -> usize {
        self.sym.pinv[r]
    }

    /// Original row eliminated at position `k`.
    pub fn perm(&self, k: usize) -> usize {
        self.sym.perm[k]
    }

    /// True when L is diagonal (every column has only its pivot).
    pub fn is_diagonal(&self) -> bool {
        self.stats.nnz_l == self.sym.m
    }

    /// Pivot `L_kk` of permuted column `k`.
    pub fn l_diag(&self, k: usize) -> f64 {
        let s = self.sym.snode_of[k] as usize;
        let nr = self.sym.sn_rows_ptr[s + 1] - self.sym.sn_rows_ptr[s];
        let j = k - self.sym.sn_first[s];
        self.l[self.sym.sn_val_ptr[s] + j * nr + j]
    }
}

impl Sym {
    /// Subtree-to-thread schedule (proportional-mapping style): split the
    /// heaviest subtree into its children until every subtree is at most
    /// `1/(2T)` of the total work, then assign subtrees to threads by
    /// longest-processing-time first. Split roots form the serial `top`.
    fn build_plan(&mut self, threads: usize) {
        let ns = self.sn_first.len() - 1;
        let mut parent = vec![NONE; ns];
        let mut cost = vec![0f64; ns];
        for k in 0..ns {
            let rl = self.sn_rows_ptr[k];
            let nr = self.sn_rows_ptr[k + 1] - rl;
            let nc = self.sn_first[k + 1] - self.sn_first[k];
            if nr > nc {
                parent[k] = self.snode_of[self.sn_rows[rl + nc] as usize];
            }
            cost[k] = sn_ops(nc, nr - nc) + 50.0 * nr as f64;
        }
        let mut sub = cost.clone();
        let mut first_desc: Vec<usize> = (0..ns).collect();
        for k in 0..ns {
            if parent[k] != NONE {
                let p = parent[k] as usize;
                sub[p] += sub[k];
                first_desc[p] = first_desc[p].min(first_desc[k]);
            }
        }
        let mut children: Vec<Vec<usize>> = vec![Vec::new(); ns];
        let mut roots = Vec::new();
        for k in 0..ns {
            if parent[k] == NONE {
                roots.push(k);
            } else {
                children[parent[k] as usize].push(k);
            }
        }
        let total: f64 = roots.iter().map(|&r| sub[r]).sum();
        let cap = total / (2.0 * threads as f64);
        let mut cand: Vec<usize> = roots;
        let mut top = Vec::new();
        for _ in 0..100_000 {
            let Some((i, &r)) = cand
                .iter()
                .enumerate()
                .max_by(|a, b| sub[*a.1].total_cmp(&sub[*b.1]))
            else {
                break;
            };
            if sub[r] <= cap {
                break;
            }
            cand.swap_remove(i);
            top.push(r);
            cand.extend(children[r].iter().copied());
        }
        cand.sort_unstable_by(|a, b| sub[*b].total_cmp(&sub[*a]));
        let mut load = vec![0f64; threads];
        let mut plan: Vec<Vec<(usize, usize)>> = vec![Vec::new(); threads];
        for r in cand {
            let t = (0..threads)
                .min_by(|a, b| load[*a].total_cmp(&load[*b]))
                .unwrap();
            load[t] += sub[r];
            plan[t].push((first_desc[r], r + 1));
        }
        for p in plan.iter_mut() {
            p.sort_unstable();
        }
        top.sort_unstable();
        self.plan_threads = threads;
        self.plan = plan;
        self.top = top;
    }

    /// Assemble, update and factor supernode `s` into `ls` (its value block).
    /// `desc` holds the values of supernodes `< s` starting at value offset
    /// `base`; `diag` the columns of `s`. Returns false on deadline.
    #[allow(clippy::too_many_arguments)]
    fn factor_sn(
        &self,
        s: usize,
        desc: &[f64],
        base: usize,
        ls: &mut [f64],
        diag: &mut [f64],
        ctl: &Ctl,
        ws: &mut Ws,
        intra: usize,
    ) -> bool {
        let f = self.sn_first[s];
        let lend = self.sn_first[s + 1];
        let rows_lo = self.sn_rows_ptr[s];
        let nr = self.sn_rows_ptr[s + 1] - rows_lo;
        ls.iter_mut().for_each(|v| *v = 0.0);
        for t in 0..nr {
            ws.relmap[self.sn_rows[rows_lo + t] as usize] = t as u32;
        }
        let theta = ctl.theta;
        // Assemble columns f..lend of M (lower part).
        for c in f..lend {
            ws.touched.clear();
            for p in self.at_ptr[c]..self.at_ptr[c + 1] {
                let j = self.at_col[p] as usize;
                let t = theta[j] * self.at_val[p];
                if t == 0.0 {
                    continue;
                }
                for q in self.a_ptr[j]..self.a_ptr[j + 1] {
                    let r = self.a_prow[q] as usize;
                    if r >= c {
                        if ws.wmark[r] != c as u32 {
                            ws.wmark[r] = c as u32;
                            ws.touched.push(r as u32);
                        }
                        ws.w[r] += t * self.a_val[q];
                    }
                }
            }
            let col = (c - f) * nr;
            for &r in &ws.touched {
                let r = r as usize;
                ls[col + ws.relmap[r] as usize] += ws.w[r];
                ws.w[r] = 0.0;
            }
            let dpos = col + (c - f);
            ls[dpos] += ctl.delta;
            diag[c - f] = ls[dpos];
        }
        // Descendant updates (static list).
        for u in self.upd_ptr[s]..self.upd_ptr[s + 1] {
            let k = self.upd_src[u] as usize;
            let p = self.upd_pos[u] as usize;
            self.apply_update(k, p, s, desc, base, ls, ws, intra);
            if ctl.past_deadline() {
                return false;
            }
        }
        // Dense panel factorization.
        let nc = lend - f;
        let pr = ctl.pivrel;
        let pp = pivot_pert();
        let dep = ctl.dep;
        let Ws {
            gw,
            skipped,
            cancel,
            modified,
            ..
        } = ws;
        let ref_diag = ctl.ref_diag;
        let mod_rel = ctl.mod_rel;
        let mod_cap = ctl.mod_cap;
        kernels::panel_cholesky(
            ls,
            nr,
            nc,
            PANEL_NB,
            intra,
            gw,
            |j, d| {
                let dorig = diag[j];
                if dep.is_some_and(|dv| dv[f + j]) {
                    *skipped += 1;
                    return Pivot::Skip(HUGE_PIVOT);
                }
                let refj = ref_diag.map(|r| r[f + j]);
                let tiny_ref = refj.is_some_and(|r| d.is_nan() || d <= mod_rel * r);
                if tiny_ref || d.is_nan() || d <= pr * dorig || d <= 1e-300 || !d.is_finite() {
                    if let Some(r) = refj
                        && r.is_finite()
                        && r > 0.0
                        && modified.len() < mod_cap
                    {
                        let dd = if d.is_finite() { d } else { 0.0 };
                        modified.push(((f + j) as u32, r - dd));
                        return Pivot::Replace(r.sqrt());
                    }
                    *skipped += 1;
                    *cancel += 1;
                    if pp > 0.0 && dorig.is_finite() && dorig > 0.0 {
                        return Pivot::Replace((pp * dorig).sqrt());
                    }
                    return Pivot::Skip(HUGE_PIVOT);
                }
                Pivot::Accept
            },
            || !ctl.past_deadline(),
        )
    }

    /// Subtract supernode `k`'s contribution (its rows from position `p`)
    /// from supernode `s`.
    #[allow(clippy::too_many_arguments)]
    fn apply_update(
        &self,
        k: usize,
        p: usize,
        s: usize,
        desc: &[f64],
        base: usize,
        ls: &mut [f64],
        ws: &mut Ws,
        intra: usize,
    ) {
        let fk = self.sn_first[k];
        let nck = self.sn_first[k + 1] - fk;
        let rk_lo = self.sn_rows_ptr[k];
        let nrk = self.sn_rows_ptr[k + 1] - rk_lo;
        let fs = self.sn_first[s];
        let lsn = self.sn_first[s + 1];
        let rows_k = &self.sn_rows[rk_lo..rk_lo + nrk];
        let mut q = p;
        while q < nrk && (rows_k[q] as usize) < lsn {
            q += 1;
        }
        let mr = nrk - p;
        let wc = q - p;
        if ws.buf.len() < mr * wc {
            ws.buf.resize(mr * wc, 0.0);
        }
        let buf = &mut ws.buf[..mr * wc];
        buf.iter_mut().for_each(|v| *v = 0.0);
        let off = self.sn_val_ptr[k] - base;
        let lk = &desc[off..off + nrk * nck];
        if mr * wc * nck >= GEMM_MIN_WORK && wc >= 2 {
            let a = &lk[p..];
            let fl = 2.0 * (mr * wc * nck) as f64;
            if intra > 1 && fl >= PAR_MIN_FLOPS {
                kernels::par_gemm_nt(mr, wc, nck, 1.0, a, nrk, a, nrk, buf, mr, true, &mut ws.gw);
            } else {
                kernels::gemm_nt(
                    mr,
                    wc,
                    nck,
                    1.0,
                    a,
                    nrk,
                    a,
                    nrk,
                    buf,
                    mr,
                    true,
                    &mut ws.gw[0],
                );
            }
        } else {
            for cc in 0..wc {
                let out = &mut buf[cc * mr..(cc + 1) * mr];
                let mut kk = 0;
                // 4-column unroll of the inner product dimension.
                while kk + 4 <= nck {
                    let c0 = &lk[kk * nrk + p..kk * nrk + nrk];
                    let c1 = &lk[(kk + 1) * nrk + p..(kk + 1) * nrk + nrk];
                    let c2 = &lk[(kk + 2) * nrk + p..(kk + 2) * nrk + nrk];
                    let c3 = &lk[(kk + 3) * nrk + p..(kk + 3) * nrk + nrk];
                    let (t0, t1, t2, t3) = (c0[cc], c1[cc], c2[cc], c3[cc]);
                    for rr in cc..mr {
                        out[rr] += c0[rr] * t0 + c1[rr] * t1 + c2[rr] * t2 + c3[rr] * t3;
                    }
                    kk += 4;
                }
                while kk < nck {
                    let c0 = &lk[kk * nrk + p..kk * nrk + nrk];
                    let t0 = c0[cc];
                    if t0 != 0.0 {
                        for rr in cc..mr {
                            out[rr] += c0[rr] * t0;
                        }
                    }
                    kk += 1;
                }
            }
        }
        let nrs = self.sn_rows_ptr[s + 1] - self.sn_rows_ptr[s];
        let relmap = &ws.relmap;
        for cc in 0..wc {
            let jcol = rows_k[p + cc] as usize - fs;
            let colbase = jcol * nrs;
            let src = &ws.buf[cc * mr..(cc + 1) * mr];
            for rr in cc..mr {
                let dst = colbase + relmap[rows_k[p + rr] as usize] as usize;
                ls[dst] -= src[rr];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::XorShift64;

    fn random_a(m: usize, n: usize, per_col: usize, seed: u64) -> CscMatrix {
        let mut rng = XorShift64::new(seed);
        let mut trip = Vec::new();
        for j in 0..n {
            let mut rows: Vec<usize> = (0..per_col)
                .map(|_| (rng.uniform(0.0, m as f64) as usize).min(m - 1))
                .collect();
            rows.sort_unstable();
            rows.dedup();
            for r in rows {
                trip.push((r, j, rng.uniform(-2.0, 2.0)));
            }
        }
        // identity block so A has full row rank
        for i in 0..m {
            trip.push((i, n + i, 1.0));
        }
        CscMatrix::from_triplets(m, n + m, &trip)
    }

    fn dense_m(a: &CscMatrix, theta: &[f64], delta: f64) -> Vec<Vec<f64>> {
        let m = a.nrow();
        let mut d = vec![vec![0.0; m]; m];
        for j in 0..a.ncol() {
            let (p0, p1) = (a.col_ptr()[j], a.col_ptr()[j + 1]);
            for p in p0..p1 {
                for q in p0..p1 {
                    let (r, s) = (a.row_idx()[p] as usize, a.row_idx()[q] as usize);
                    d[r][s] += theta[j] * a.values()[p] * a.values()[q];
                }
            }
        }
        for (i, row) in d.iter_mut().enumerate() {
            row[i] += delta;
        }
        d
    }

    fn check(a: &CscMatrix, seed: u64, relax: f64, threads: usize) {
        let m = a.nrow();
        let mut rng = XorShift64::new(seed + 5);
        let theta: Vec<f64> = (0..a.ncol())
            .map(|_| 10f64.powf(rng.uniform(-4.0, 4.0)))
            .collect();
        let mut ch = NormalChol::analyse_opts(a, usize::MAX, relax, threads).unwrap();
        ch.factor(&theta, 1e-10);
        assert!(!ch.aborted);
        assert_eq!(ch.stats.skipped_pivots, 0);
        let b: Vec<f64> = (0..m).map(|i| (i as f64 * 0.37).sin()).collect();
        let mut x = b.clone();
        let mut work = Vec::new();
        ch.solve(&mut x, &mut work);
        let d = dense_m(a, &theta, 1e-10);
        let mut res = 0f64;
        let mut scale = 0f64;
        for i in 0..m {
            let s: f64 = (0..m).map(|j| d[i][j] * x[j]).sum();
            res = res.max((s - b[i]).abs());
            scale = scale.max(d[i].iter().map(|v| v.abs()).fold(0.0, f64::max));
        }
        assert!(
            res <= 1e-8 * scale.max(1.0),
            "seed {seed} relax {relax} threads {threads}: residual {res}"
        );
        // matvec agrees with the dense product
        let mut y = vec![0.0; m];
        ch.matvec(&theta, 1e-10, &x, &mut y);
        for i in 0..m {
            let s: f64 = (0..m).map(|j| d[i][j] * x[j]).sum();
            assert!((y[i] - s).abs() <= 1e-9 * scale.max(1.0));
        }
    }

    #[test]
    fn factor_solve_matches_dense_residual() {
        for (seed, (m, n, k)) in [
            (1u64, (30, 60, 3)),
            (2, (80, 120, 5)),
            (3, (50, 20, 12)),
            (4, (300, 400, 4)),
            (5, (200, 150, 30)),
        ] {
            let a = random_a(m, n, k, seed + 11);
            for relax in [0.0, 0.02, 0.3] {
                for threads in [1, 3] {
                    check(&a, seed, relax, threads);
                }
            }
        }
    }

    /// Static update lists fix the update order, so the solve is
    /// bit-identical at any thread count (design note §9.1).
    #[test]
    fn threads_give_bit_identical_solves() {
        let a = random_a(300, 400, 4, 17);
        let mut rng = XorShift64::new(23);
        let theta: Vec<f64> = (0..a.ncol())
            .map(|_| 10f64.powf(rng.uniform(-4.0, 4.0)))
            .collect();
        let b: Vec<f64> = (0..a.nrow()).map(|i| (i as f64 * 0.37).sin()).collect();
        let mut sols = Vec::new();
        for threads in [1, 2, 3, 8] {
            let mut ch = NormalChol::analyse_opts(&a, usize::MAX, 0.02, threads).unwrap();
            assert!(
                threads == 1 || ch.stats.subtrees > 1,
                "no parallel subtrees"
            );
            ch.factor(&theta, 1e-10);
            let mut x = b.clone();
            let mut work = Vec::new();
            ch.solve(&mut x, &mut work);
            sols.push(x.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
        }
        for s in &sols[1..] {
            assert_eq!(&sols[0], s);
        }
    }

    #[test]
    fn relaxation_merges_and_keeps_fill() {
        let a = random_a(400, 500, 3, 9);
        let c0 = NormalChol::analyse_opts(&a, usize::MAX, 0.0, 1).unwrap();
        let c1 = NormalChol::analyse_opts(&a, usize::MAX, 0.3, 1).unwrap();
        assert_eq!(c0.stats.nnz_l, c1.stats.nnz_l);
        assert_eq!(c0.stats.nnz_l, c0.stats.nnz_l_stored);
        assert!(c1.stats.supernodes < c0.stats.supernodes);
        assert!(c1.stats.nnz_l_stored >= c1.stats.nnz_l);
    }

    #[test]
    fn dependent_rows_are_skipped_not_nan() {
        // Row 2 = row 0 + row 1 → A·Aᵀ singular.
        let trip = vec![
            (0, 0, 1.0),
            (1, 1, 1.0),
            (2, 0, 1.0),
            (2, 1, 1.0),
            (0, 2, 2.0),
            (2, 2, 2.0),
        ];
        let a = CscMatrix::from_triplets(3, 3, &trip);
        let mut ch = NormalChol::analyse(&a, usize::MAX).unwrap();
        ch.factor(&[1.0, 1.0, 1.0], 0.0);
        assert_eq!(ch.stats.skipped_pivots, 1);
        let mut x = vec![1.0, 2.0, 3.0];
        let mut work = Vec::new();
        ch.solve(&mut x, &mut work);
        assert!(x.iter().all(|v| v.is_finite()));
    }
}
