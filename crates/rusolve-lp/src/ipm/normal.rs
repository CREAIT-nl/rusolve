//! Normal-equations solver with dense-column splitting.
//!
//! `M = A·Θ·Aᵀ + δI` is factored as a sparse part plus a low-rank term:
//! `M = M_s' + U·C·Uᵀ` where `M_s = A_s·Θ_s·A_sᵀ + δI` holds the sparse
//! columns, `M_s'` is `M_s` with its tiny pivots raised to the diagonal of the
//! full `M` (`M_s' = M_s + Σ c_j e_j e_jᵀ`), `U = [A_d, E]` and
//! `C = diag(Θ_d, −c)`. Solves use the Sherman–Morrison–Woodbury identity
//! with the Schur complement `K = C⁻¹ + Uᵀ·M_s'⁻¹·U`.
//!
//! Provenance:
//! * dense-column split with the *modified* Schur complement (tiny pivots of
//!   the sparse part are replaced and corrected through the low-rank term,
//!   so the split stays stable when the dense columns carry a row):
//!   E. D. Andersen, "A modified Schur complement method for handling dense
//!   columns in interior-point methods for linear programming", ACM TOMS
//!   22(3) (1996); E. D. Andersen, C. Roos & T. Terlaky, "On implementing a
//!   primal-dual interior-point method for conic quadratic optimization",
//!   Math. Program. 95 (2003) §dense columns;
//! * Sherman–Morrison–Woodbury: Golub & Van Loan, *Matrix Computations*
//!   (4th ed.) §2.1.4; Wright, *Primal-Dual Interior-Point Methods* (SIAM
//!   1997) ch. 11 (dense columns in normal equations);
//! * candidate dense sets from the sorted column counts and the 40-nonzero
//!   floor: HiGHS IPX `Model::FindDenseColumns` (MIT, read as reference);
//!   the choice among candidates is ours: a symbolic flop/solve cost model
//!   evaluated before any numeric factorization, like HiPO's NE-vs-AS choice by predicted
//!   operation counts (`FactorHiGHSSolver`, MIT);
//! * the block elimination of `K` (Cholesky of the SPD dense-column block,
//!   LU with partial pivoting of the small Schur complement of the modified
//!   pivots): Golub & Van Loan §3.2, §4.2.

use std::time::Instant;

use rusolve_core::CscMatrix;

use super::chol::{CholError, CholStats, NormalChol, symbolic_cost};
use super::kernels::{self, GemmWork, Pivot};

/// Largest dense set considered (`K` is `k × k` dense).
const KMAX: usize = 8000;
/// IPX's floor: a column with ≤ 40 nonzeros is never dense.
const DENSE_MIN_NNZ: usize = 40;
/// Store `V = L⁻¹·U` when it has at most this many nonzeros.
const V_STORE_MAX: usize = 32 << 20;
/// Solves per IPM iteration assumed by the cost model (predictor, correctors,
/// refinement).
const COST_SOLVES: f64 = 8.0;
/// A split must predict at most this share of the unsplit cost.
const COST_GAIN: f64 = 0.7;
/// PCG on the full `M` after a split solve: step cap, relative stop, and
/// the number of steps without halving the best residual that ends it.
const SPLIT_PCG: usize = 30;
const SPLIT_TOL: f64 = 1e-13;
const SPLIT_STALL: usize = 4;
/// Refinement steps when the SMW preconditioner is indefinite.
const SPLIT_REFINE: usize = 4;
/// Split health probe: relative residual that counts as lost accuracy, and
/// when the unsplit factorization is affordable enough to switch to it.
const PROBE_BAD: f64 = 1e-3;
const UNSPLIT_GAIN: f64 = 20.0;
const UNSPLIT_FLOPS: f64 = 1e9;

/// `RUSOLVE_IPM_DENSE`: 0 = never split, 1 = cost model (default),
/// 2 = always take the cheapest non-empty candidate (testing lever).
fn dense_mode() -> u32 {
    static V: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("RUSOLVE_IPM_DENSE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(1)
    })
}
/// `RUSOLVE_IPM_DENSE_MODREL`: a sparse-part pivot `d ≤ modrel · M_jj` is
/// modified (Andersen 1996). TODO(perf): only 1e-8 measured (eval note §3).
fn mod_rel() -> f64 {
    static V: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("RUSOLVE_IPM_DENSE_MODREL")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(1e-8)
    })
}

/// Dense-column choice and its predicted costs (logged).
#[derive(Debug, Clone, Default)]
pub struct DenseChoice {
    pub cols: Vec<usize>,
    /// (k, predicted cost) per evaluated candidate; cost = ∞ if TooLarge.
    pub tried: Vec<(usize, f64)>,
    pub t_select: f64,
}

fn predicted_cost(m: usize, k: usize, nnz_d: usize, c: Option<(f64, usize)>) -> f64 {
    let Some((flops, nnz_l)) = c else {
        return f64::INFINITY;
    };
    let (kf, mf, nl) = (k as f64, m as f64, nnz_l as f64);
    let w_stored = m * k <= V_STORE_MAX;
    let build = if k == 0 {
        0.0
    } else {
        4.0 * nl * kf + 2.0 * kf * nnz_d as f64 + kf * kf * kf / 3.0
    };
    let per_solve = 4.0 * nl
        + if k == 0 {
            0.0
        } else {
            2.0 * mf * kf + 2.0 * kf * kf
        }
        + if w_stored || k == 0 { 0.0 } else { 4.0 * nl };
    flops + build + COST_SOLVES * per_solve
}

fn select_columns(a: &CscMatrix, cols: &[usize]) -> CscMatrix {
    let ptr = a.col_ptr();
    let ri = a.row_idx();
    let v = a.values();
    let mut trip = Vec::new();
    for (jn, &j) in cols.iter().enumerate() {
        for p in ptr[j]..ptr[j + 1] {
            trip.push((ri[p] as usize, jn, v[p]));
        }
    }
    CscMatrix::from_triplets(a.nrow(), cols.len(), &trip)
}

/// Choose the dense columns by the symbolic cost model.
pub fn choose_dense(a: &CscMatrix, max_entries: usize) -> DenseChoice {
    let t0 = Instant::now();
    let mut ch = DenseChoice::default();
    let mode = dense_mode();
    let (m, n) = (a.nrow(), a.ncol());
    if mode == 0 || n == 0 || m == 0 {
        return ch;
    }
    let ptr = a.col_ptr();
    let cnt: Vec<usize> = (0..n).map(|j| ptr[j + 1] - ptr[j]).collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by(|&x, &y| cnt[y].cmp(&cnt[x]).then(x.cmp(&y)));
    let cmax = cnt[order[0]];
    if cmax <= DENSE_MIN_NNZ {
        return ch;
    }
    let tot2: f64 = cnt.iter().map(|&c| (c * c) as f64).sum();
    let mut pre2 = vec![0f64; n + 1];
    for (i, &j) in order.iter().enumerate() {
        pre2[i + 1] = pre2[i] + (cnt[j] * cnt[j]) as f64;
    }
    // Candidates: k(τ) = #{c_j ≥ τ} for τ = cmax/2^i, τ > 40, with the dense
    // set carrying ≥ 25 % of the pattern work Σc².
    let mut cands: Vec<usize> = Vec::new();
    let mut tau = cmax as f64;
    while tau > DENSE_MIN_NNZ as f64 {
        let k = order.partition_point(|&j| cnt[j] as f64 >= tau);
        if k > 0 && k <= KMAX.min(n) && pre2[k] >= 0.25 * tot2 && cands.last() != Some(&k) {
            cands.push(k);
        }
        tau *= 0.5;
    }
    if cands.is_empty() {
        return ch;
    }
    cands.truncate(6);
    let sc = |b: &CscMatrix| {
        symbolic_cost(b, max_entries)
            .ok()
            .filter(|c| c.nnz_l <= max_entries)
            .map(|c| (c.flops, c.nnz_l))
    };
    // The unsplit pattern holds every dense column's clique.
    let c0 = if cmax * (cmax - 1) > max_entries {
        None
    } else {
        sc(a)
    };
    let f0 = predicted_cost(m, 0, 0, c0);
    ch.tried.push((0, f0));
    let mut best = (f0, 0usize);
    let mut best_split = (f64::INFINITY, 0usize);
    for &k in &cands {
        let mut sparse: Vec<usize> = order[k..].to_vec();
        sparse.sort_unstable();
        let nnz_d: usize = order[..k].iter().map(|&j| cnt[j]).sum();
        let f = predicted_cost(m, k, nnz_d, sc(&select_columns(a, &sparse)));
        ch.tried.push((k, f));
        if f < best_split.0 {
            best_split = (f, k);
        }
        if f < COST_GAIN * best.0 || (best.0.is_infinite() && f < best.0) {
            best = (f, k);
        }
    }
    let k = if mode == 2 && best_split.0.is_finite() {
        best_split.1
    } else {
        best.1
    };
    let mut cols = order[..k].to_vec();
    cols.sort_unstable();
    ch.cols = cols;
    ch.t_select = t0.elapsed().as_secs_f64();
    ch
}

/// Normal-equations solver: sparse Cholesky, optionally with dense columns
/// split off and recovered through the modified Schur complement.
pub struct NormalEq {
    pub chol: NormalChol,
    m: usize,
    /// Full A (original rows) for `matvec` and the reference diagonal.
    a_ptr: Vec<usize>,
    a_row: Vec<u32>,
    a_val: Vec<f64>,
    /// Columns of A in the sparse part (A_s column → A column).
    scols: Vec<usize>,
    /// Dense columns (indices into A); their rows in permuted order.
    dcols: Vec<usize>,
    d_ptr: Vec<usize>,
    d_prow: Vec<u32>,
    d_val: Vec<f64>,
    theta_s: Vec<f64>,
    /// Θ and δ of the last factorization (for the refinement residual).
    theta: Vec<f64>,
    delta: f64,
    ref_diag: Vec<f64>,
    /// Current low-rank term: k dense columns then `mods` unit columns.
    mods: Vec<(u32, f64)>,
    /// V = L⁻¹·U (permuted rows) as sparse columns, when stored.
    v_ptr: Vec<usize>,
    v_idx: Vec<u32>,
    v_val: Vec<f64>,
    v_stored: bool,
    /// Block factor of K: L11 (k × k), Y = L11⁻¹K12 (k × q), LU of S (q × q).
    l11: Vec<f64>,
    y: Vec<f64>,
    s_lu: Vec<f64>,
    s_piv: Vec<usize>,
    pub stats: CholStats,
    pub choice: DenseChoice,
    pub deadline: Option<Instant>,
    pub aborted: bool,
    /// Full A and the Cholesky size cap, kept while split so the solver can
    /// fall back to the unsplit factorization (`probe_or_unsplit`).
    full_a: Option<CscMatrix>,
    max_entries: usize,
    /// (split solves, PCG steps), for the log.
    pcg: std::cell::Cell<(usize, usize)>,
}

impl NormalEq {
    pub fn analyse(a: &CscMatrix, max_entries: usize) -> Result<Self, CholError> {
        let choice = choose_dense(a, max_entries);
        Self::analyse_with(a, max_entries, choice)
    }

    pub fn analyse_with(
        a: &CscMatrix,
        max_entries: usize,
        choice: DenseChoice,
    ) -> Result<Self, CholError> {
        let (m, n) = (a.nrow(), a.ncol());
        let dcols = choice.cols.clone();
        let mut is_d = vec![false; n];
        for &j in &dcols {
            is_d[j] = true;
        }
        let scols: Vec<usize> = (0..n).filter(|&j| !is_d[j]).collect();
        let chol = if dcols.is_empty() {
            NormalChol::analyse(a, max_entries)?
        } else {
            NormalChol::analyse(&select_columns(a, &scols), max_entries)?
        };
        let ptr = a.col_ptr();
        let ri = a.row_idx();
        let v = a.values();
        let mut d_ptr = vec![0usize];
        let mut d_prow = Vec::new();
        let mut d_val = Vec::new();
        for &j in &dcols {
            for p in ptr[j]..ptr[j + 1] {
                d_prow.push(chol.pinv(ri[p] as usize) as u32);
                d_val.push(v[p]);
            }
            d_ptr.push(d_prow.len());
        }
        let dcols_empty = dcols.is_empty();
        let mut stats = chol.stats.clone();
        stats.dense_split = dcols.len();
        stats.nnz_a = a.nnz();
        Ok(Self {
            m,
            a_ptr: ptr.to_vec(),
            a_row: ri.iter().map(|&r| r as u32).collect(),
            a_val: v.to_vec(),
            theta_s: vec![0.0; scols.len()],
            theta: Vec::new(),
            delta: 0.0,
            scols,
            dcols,
            d_ptr,
            d_prow,
            d_val,
            ref_diag: vec![0.0; m],
            mods: Vec::new(),
            v_ptr: Vec::new(),
            v_idx: Vec::new(),
            v_val: Vec::new(),
            v_stored: false,
            l11: Vec::new(),
            y: Vec::new(),
            s_lu: Vec::new(),
            s_piv: Vec::new(),
            full_a: if dcols_empty { None } else { Some(a.clone()) },
            max_entries,
            pcg: std::cell::Cell::new((0, 0)),
            chol,
            stats,
            choice,
            deadline: None,
            aborted: false,
        })
    }

    fn split(&self) -> bool {
        !self.dcols.is_empty()
    }

    /// Factor-cost / solve-cost proxies for the corrector count.
    pub fn cost_ratio(&self) -> f64 {
        let k = self.dcols.len() as f64;
        let nl = self.chol.stats.nnz_l as f64;
        let fl = self.chol.stats.flops + 4.0 * nl * k + k * k * k / 3.0;
        let sol = 4.0 * nl + 2.0 * self.m as f64 * k + 2.0 * k * k;
        fl / sol.max(1.0)
    }

    /// Scatter low-rank column `c` (permuted rows) into `v` (zeroed).
    fn u_col(&self, c: usize, v: &mut [f64]) {
        let k = self.dcols.len();
        if c < k {
            for p in self.d_ptr[c]..self.d_ptr[c + 1] {
                v[self.d_prow[p] as usize] = self.d_val[p];
            }
        } else {
            v[self.mods[c - k].0 as usize] = 1.0;
        }
    }

    /// `u_cᵀ·x` (permuted rows).
    fn u_dot(&self, c: usize, x: &[f64]) -> f64 {
        let k = self.dcols.len();
        if c < k {
            (self.d_ptr[c]..self.d_ptr[c + 1])
                .map(|p| self.d_val[p] * x[self.d_prow[p] as usize])
                .sum()
        } else {
            x[self.mods[c - k].0 as usize]
        }
    }

    /// [`NormalChol::factor_detect`] on the unsplit path (the split path's
    /// sparse part is not the row space of A, so it is not probed).
    pub fn factor_detect(&mut self, theta: &[f64], delta: f64, rel: f64) -> usize {
        if self.split() || rel <= 0.0 {
            self.factor(theta, delta);
            return 0;
        }
        self.chol.deadline = self.deadline;
        let d = self.chol.factor_detect(theta, delta, rel);
        self.aborted = self.chol.aborted;
        self.copy_pivot_stats();
        self.stats.dependent = d;
        d
    }

    fn copy_pivot_stats(&mut self) {
        self.stats.skipped_pivots = self.chol.stats.skipped_pivots;
        self.stats.cancel_pivots = self.chol.stats.cancel_pivots;
    }

    pub fn factor(&mut self, theta: &[f64], delta: f64) {
        self.chol.deadline = self.deadline;
        self.aborted = false;
        if !self.split() {
            self.chol.factor(theta, delta);
            self.aborted = self.chol.aborted;
            self.copy_pivot_stats();
            return;
        }
        for (t, &j) in self.theta_s.iter_mut().zip(&self.scols) {
            *t = theta[j];
        }
        self.theta.clear();
        self.theta.extend_from_slice(theta);
        self.delta = delta;
        // Full diagonal of M (original rows) as the modification reference.
        self.ref_diag.iter_mut().for_each(|v| *v = delta);
        for j in 0..self.a_ptr.len() - 1 {
            let t = theta[j];
            for p in self.a_ptr[j]..self.a_ptr[j + 1] {
                self.ref_diag[self.a_row[p] as usize] += t * self.a_val[p] * self.a_val[p];
            }
        }
        let k = self.dcols.len();
        let cap = 200usize.max(k).min(KMAX);
        self.chol
            .factor_modified(&self.theta_s, delta, &self.ref_diag, mod_rel(), cap);
        if self.chol.aborted {
            self.aborted = true;
            return;
        }
        self.mods = std::mem::take(&mut self.chol.modified);
        if self.mods.len() > cap {
            // Too many modified pivots for a dense K: plain factor of M_s'
            // (tiny pivots skip and count as cancellation, which raises the
            // regularization in the solver). Thread-count independent.
            self.mods.clear();
            self.chol.factor(&self.theta_s, delta);
            if self.chol.aborted {
                self.aborted = true;
                return;
            }
        }
        self.mods.sort_unstable_by_key(|v| v.0);
        let q = self.mods.len();
        let kt = k + q;
        self.stats.dense_mods = q;
        let m = self.m;
        // Schur complement K = C⁻¹ + Uᵀ M_s'⁻¹ U (kt × kt, column-major).
        // With M_s' = L·Lᵀ, Uᵀ M_s'⁻¹ U = VᵀV for V = L⁻¹U, whose columns
        // come from sparse forward solves (only the elimination-tree reach of
        // each column of U is touched; neos3, L nearly diagonal with 512k
        // supernodes, spent 17 s per factor in full solves before).
        let mut kmat = vec![0.0; kt * kt];
        if !self.build_v(kt) {
            return;
        }
        if self.v_stored {
            let nnz_v = self.v_idx.len();
            if (nnz_v as f64) > 0.25 * (m * kt) as f64 {
                // Nearly dense V: K = VᵀV by the packed GEMM kernel.
                let mut x = vec![0.0; kt * m];
                for c in 0..kt {
                    for p in self.v_ptr[c]..self.v_ptr[c + 1] {
                        x[c + self.v_idx[p] as usize * kt] = self.v_val[p];
                    }
                }
                let mut gw = GemmWork::default();
                kernels::gemm_nt(kt, kt, m, 1.0, &x, kt, &x, kt, &mut kmat, kt, true, &mut gw);
            } else {
                // Row-wise outer products of the sparse rows of V.
                let mut rp = vec![0usize; m + 1];
                for &i in &self.v_idx {
                    rp[i as usize + 1] += 1;
                }
                for i in 0..m {
                    rp[i + 1] += rp[i];
                }
                let mut fill = rp.clone();
                let mut rc = vec![(0u32, 0.0f64); nnz_v];
                for c in 0..kt {
                    for p in self.v_ptr[c]..self.v_ptr[c + 1] {
                        let i = self.v_idx[p] as usize;
                        rc[fill[i]] = (c as u32, self.v_val[p]);
                        fill[i] += 1;
                    }
                }
                for i in 0..m {
                    let row = &rc[rp[i]..rp[i + 1]];
                    for (t, &(a, va)) in row.iter().enumerate() {
                        for &(b, vb) in &row[t..] {
                            // Columns ascend within a row: b ≥ a (lower part).
                            kmat[b as usize + a as usize * kt] += va * vb;
                        }
                    }
                }
            }
            for c in 0..kt {
                for r in c + 1..kt {
                    kmat[c + r * kt] = kmat[r + c * kt];
                }
            }
        } else {
            let mut v = vec![0.0; m];
            for c in 0..kt {
                if self.deadline.is_some_and(|d| Instant::now() >= d) {
                    self.aborted = true;
                    return;
                }
                v.iter_mut().for_each(|x| *x = 0.0);
                self.u_col(c, &mut v);
                self.chol.solve_permuted(&mut v);
                for r in c..kt {
                    let d = self.u_dot(r, &v);
                    kmat[r + c * kt] = d;
                    kmat[c + r * kt] = d;
                }
            }
        }
        for c in 0..k {
            kmat[c + c * kt] += 1.0 / theta[self.dcols[c]];
        }
        for (t, &(_, cj)) in self.mods.iter().enumerate() {
            let c = k + t;
            kmat[c + c * kt] -= 1.0 / cj;
        }
        self.factor_k(&kmat, k, q);
        self.copy_pivot_stats();
        self.probe_or_unsplit(theta, delta);
    }

    /// Split health check (only with modified pivots): one SMW solve of
    /// `M x = M·u`, u the indicator of the modified rows. When its relative
    /// residual exceeds `PROBE_BAD` the modified Schur complement has lost
    /// all accuracy (nidda late in the IPM: residual 11× the right-hand
    /// side with 13 modified pivots, K22 = −1/c_j + (M_s'⁻¹)_jj cancels;
    /// Andersen 1996 §4 notes the same instability). If the unsplit factor
    /// is affordable (predicted cost finite and ≤ `UNSPLIT_GAIN`× the split
    /// one, or ≤ `UNSPLIT_FLOPS`), switch to it for the rest of the run;
    /// otherwise keep the split and rely on the PCG in `solve`.
    fn probe_or_unsplit(&mut self, theta: &[f64], delta: f64) {
        let m = self.m;
        if self.mods.is_empty() {
            return;
        }
        // x = indicator of the modified rows: the directions where the
        // low-rank correction has to cancel.
        let mut u = vec![0.0; m];
        for &(r, _) in &self.mods {
            u[self.chol.perm(r as usize)] = 1.0;
        }
        let mut b = vec![0.0; m];
        self.matvec(theta, delta, &u, &mut b);
        let bn = b.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        let mut x = b.clone();
        let mut work = Vec::new();
        self.solve_smw(&mut x, &mut work);
        let mut r = vec![0.0; m];
        self.matvec(theta, delta, &x, &mut r);
        let rn = r
            .iter()
            .zip(&b)
            .fold(0.0f64, |a, (ri, bi)| a.max((bi - ri).abs()));
        if rn <= PROBE_BAD * bn {
            return;
        }
        let f0 = self.choice.tried.first().map_or(f64::INFINITY, |t| t.1);
        let k = self.dcols.len();
        let fk = self
            .choice
            .tried
            .iter()
            .find(|t| t.0 == k)
            .map_or(f64::INFINITY, |t| t.1);
        if !(f0.is_finite() && (f0 <= UNSPLIT_GAIN * fk || f0 <= UNSPLIT_FLOPS)) {
            return;
        }
        let Some(a) = self.full_a.take() else {
            return;
        };
        let Ok(chol) = NormalChol::analyse(&a, self.max_entries) else {
            self.full_a = Some(a);
            return;
        };
        if std::env::var_os("RUSOLVE_IPM_LOG").is_some() {
            eprintln!(
                "rusolve ipm: dense split k {k} lost accuracy (probe residual {:.1e} of |b|); unsplit nnz(L) {}",
                rn / bn.max(f64::MIN_POSITIVE),
                chol.stats.nnz_l
            );
        }
        self.chol = chol;
        self.dcols.clear();
        self.mods.clear();
        self.v_ptr.clear();
        self.v_idx = Vec::new();
        self.v_val = Vec::new();
        self.v_stored = false;
        let dense_split = self.stats.dense_split;
        self.stats = self.chol.stats.clone();
        self.stats.dense_split = dense_split;
        self.stats.unsplit = true;
        self.stats.nnz_a = a.nnz();
        self.chol.deadline = self.deadline;
        self.chol.factor(theta, delta);
        self.aborted = self.chol.aborted;
        self.copy_pivot_stats();
    }

    /// V = L⁻¹·U column by column (sparse forward solves). Leaves
    /// `v_stored = false` when V would exceed `V_STORE_MAX` entries (the
    /// caller then forms K and solves with full triangular solves). Returns
    /// false on the deadline.
    fn build_v(&mut self, kt: usize) -> bool {
        let m = self.m;
        let k = self.dcols.len();
        self.v_ptr.clear();
        self.v_ptr.push(0);
        self.v_idx.clear();
        self.v_val.clear();
        self.v_stored = true;
        let mut y = vec![0.0; m];
        let mut mark = vec![false; self.chol.n_supernodes()];
        let mut reach = Vec::new();
        let mut nz: Vec<u32> = Vec::new();
        for c in 0..kt {
            if self.deadline.is_some_and(|d| Instant::now() >= d) {
                self.aborted = true;
                return false;
            }
            nz.clear();
            if c < k {
                for p in self.d_ptr[c]..self.d_ptr[c + 1] {
                    y[self.d_prow[p] as usize] = self.d_val[p];
                    nz.push(self.d_prow[p]);
                }
            } else {
                let r = self.mods[c - k].0;
                y[r as usize] = 1.0;
                nz.push(r);
            }
            self.chol.forward_sparse(&mut y, &nz, &mut mark, &mut reach);
            for &sn in &reach {
                for i in self.chol.sn_cols(sn as usize) {
                    let v = y[i];
                    y[i] = 0.0;
                    if v != 0.0 {
                        self.v_idx.push(i as u32);
                        self.v_val.push(v);
                    }
                }
            }
            if self.v_idx.len() > V_STORE_MAX {
                self.v_stored = false;
                self.v_ptr.clear();
                self.v_idx = Vec::new();
                self.v_val = Vec::new();
                return true;
            }
            self.v_ptr.push(self.v_idx.len());
        }
        true
    }

    /// Block factorization of K: Cholesky of the SPD dense block K11, then
    /// LU (partial pivoting) of S = K22 − K21·K11⁻¹·K12.
    fn factor_k(&mut self, kmat: &[f64], k: usize, q: usize) {
        let kt = k + q;
        self.l11.clear();
        self.l11.resize(k * k, 0.0);
        for c in 0..k {
            self.l11[c * k + c..(c + 1) * k].copy_from_slice(&kmat[c * kt + c..c * kt + k]);
        }
        let diag: Vec<f64> = (0..k).map(|c| kmat[c + c * kt]).collect();
        let mut gw = vec![GemmWork::default()];
        let mut skipped = 0usize;
        kernels::panel_cholesky(
            &mut self.l11,
            k,
            k,
            64,
            1,
            &mut gw,
            |j, d| {
                if d.is_nan() || d <= 1e-14 * diag[j] || d <= 1e-300 {
                    skipped += 1;
                    Pivot::Skip(1e64)
                } else {
                    Pivot::Accept
                }
            },
            || true,
        );
        self.chol.stats.skipped_pivots += skipped;
        // Y = L11⁻¹ K12 (k × q).
        self.y.clear();
        self.y.resize(k * q, 0.0);
        for t in 0..q {
            let col = &mut self.y[t * k..(t + 1) * k];
            for (r, v) in col.iter_mut().enumerate() {
                *v = kmat[r + (k + t) * kt];
            }
            lower_solve(&self.l11, k, col);
        }
        // S = K22 − Yᵀ Y.
        self.s_lu.clear();
        self.s_lu.resize(q * q, 0.0);
        for b in 0..q {
            for a in 0..q {
                let mut s = kmat[(k + a) + (k + b) * kt];
                let (ya, yb) = (&self.y[a * k..(a + 1) * k], &self.y[b * k..(b + 1) * k]);
                for r in 0..k {
                    s -= ya[r] * yb[r];
                }
                self.s_lu[a + b * q] = s;
            }
        }
        self.s_piv = lu_factor(&mut self.s_lu, q);
    }

    /// Solve `M x = rhs` in place (original row order). `work` is scratch.
    /// With a split, the SMW solve is the preconditioner of a conjugate
    /// gradient iteration on the full `M` (at most `SPLIT_PCG` steps, the
    /// iterate with the smallest residual is returned). The SMW solve alone
    /// is exact only in exact arithmetic: when the dense columns carry most
    /// of the diagonal of many rows (cont1: 80k of 121k rows late in the
    /// IPM), `M_s'` is ill-conditioned and the low-rank correction cancels,
    /// leaving a relative error near 1 that iterative refinement (a
    /// Richardson iteration) cannot remove, while PCG removes the few bad
    /// eigendirections in a few steps. Sources: Andersen, Gondzio, Mészáros
    /// & Xu (1996) §5.2 and Adler, Karmarkar, Resende & Veiga, "Data
    /// structures and programming techniques for the implementation of
    /// Karmarkar's algorithm", ORSA J. Comput. 1(2) (1989) (PCG on the
    /// normal equations preconditioned by the sparse-part Cholesky); Golub &
    /// Van Loan §11.5 (PCG).
    pub fn solve(&self, rhs: &mut [f64], work: &mut Vec<f64>) {
        if !self.split() {
            self.chol.solve(rhs, work);
            return;
        }
        let m = self.m;
        let nrm = |v: &[f64]| v.iter().fold(0.0f64, |a, t| a.max(t.abs()));
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let b = rhs.to_vec();
        let bn = nrm(&b);
        if bn == 0.0 {
            return;
        }
        let tol = SPLIT_TOL * bn;
        // x0 = SMW solve; then PCG on the residual equation from x0.
        let mut x = b.clone();
        self.solve_smw(&mut x, work);
        let mut r = vec![0.0; m];
        let mut q = vec![0.0; m];
        self.matvec(&self.theta, self.delta, &x, &mut q);
        for t in 0..m {
            r[t] = b[t] - q[t];
        }
        let mut best = nrm(&r);
        rhs.copy_from_slice(&x);
        let mut steps = 0;
        if best.is_nan() || best <= tol {
            self.pcg.set((self.pcg.get().0 + 1, self.pcg.get().1));
            return;
        }
        let mut z = r.clone();
        self.solve_smw(&mut z, work);
        let mut rz = dot(&r, &z);
        let mut p = z.clone();
        let mut stall = 0;
        let mut indefinite = rz.is_nan() || rz <= 0.0;
        while !indefinite && steps < SPLIT_PCG {
            self.matvec(&self.theta, self.delta, &p, &mut q);
            let pq = dot(&p, &q);
            if pq.is_nan() || pq <= 0.0 {
                indefinite = steps == 0;
                break;
            }
            let alpha = rz / pq;
            for t in 0..m {
                x[t] += alpha * p[t];
                r[t] -= alpha * q[t];
            }
            steps += 1;
            let rn = nrm(&r);
            if rn < best {
                stall = if rn < 0.5 * best { 0 } else { stall + 1 };
                best = rn;
                rhs.copy_from_slice(&x);
            } else {
                stall += 1;
            }
            if rn.is_nan() || rn <= tol || stall >= SPLIT_STALL {
                break;
            }
            z.copy_from_slice(&r);
            self.solve_smw(&mut z, work);
            let rz1 = dot(&r, &z);
            if rz1.is_nan() || rz1 <= 0.0 {
                break;
            }
            let beta = rz1 / rz;
            rz = rz1;
            for t in 0..m {
                p[t] = z[t] + beta * p[t];
            }
        }
        if indefinite {
            // The computed SMW inverse is not positive on r (modified pivots
            // make K indefinite): plain iterative refinement from x0, each
            // step kept only while the residual decreases.
            x.copy_from_slice(rhs);
            self.matvec(&self.theta, self.delta, &x, &mut q);
            for t in 0..m {
                r[t] = b[t] - q[t];
            }
            for _ in 0..SPLIT_REFINE {
                z.copy_from_slice(&r);
                self.solve_smw(&mut z, work);
                for t in 0..m {
                    z[t] += x[t];
                }
                self.matvec(&self.theta, self.delta, &z, &mut q);
                for t in 0..m {
                    q[t] = b[t] - q[t];
                }
                let rn = nrm(&q);
                steps += 1;
                if rn.is_nan() || rn >= best {
                    break;
                }
                best = rn;
                x.copy_from_slice(&z);
                r.copy_from_slice(&q);
                if rn.is_nan() || rn <= tol {
                    break;
                }
            }
            rhs.copy_from_slice(&x);
        }
        self.pcg
            .set((self.pcg.get().0 + 1, self.pcg.get().1 + steps));
    }

    /// (split solves, PCG/refinement steps) since analysis.
    pub fn pcg_counts(&self) -> (usize, usize) {
        self.pcg.get()
    }

    /// One SMW solve `x = M_s'⁻¹b − M_s'⁻¹U·K⁻¹·Uᵀ·M_s'⁻¹b`; with V stored,
    /// `x = L⁻ᵀ(y − V·K⁻¹·Vᵀy)` with `y = L⁻¹b` (two triangular solves).
    fn solve_smw(&self, rhs: &mut [f64], work: &mut Vec<f64>) {
        let m = self.m;
        let k = self.dcols.len();
        let kt = k + self.mods.len();
        work.resize(2 * m, 0.0);
        let (x, v) = work.split_at_mut(m);
        for (kk, xv) in x.iter_mut().enumerate() {
            *xv = rhs[self.chol.perm(kk)];
        }
        let mut t: Vec<f64> = if self.v_stored {
            self.chol.forward(x);
            (0..kt)
                .map(|c| {
                    (self.v_ptr[c]..self.v_ptr[c + 1])
                        .map(|p| self.v_val[p] * x[self.v_idx[p] as usize])
                        .sum()
                })
                .collect()
        } else {
            self.chol.solve_permuted(x);
            (0..kt).map(|c| self.u_dot(c, x)).collect()
        };
        self.solve_k(&mut t);
        if self.v_stored {
            for (c, &zc) in t.iter().enumerate() {
                if zc != 0.0 {
                    for p in self.v_ptr[c]..self.v_ptr[c + 1] {
                        x[self.v_idx[p] as usize] -= self.v_val[p] * zc;
                    }
                }
            }
            self.chol.backward(x);
        } else {
            v.iter_mut().for_each(|e| *e = 0.0);
            for (c, &zc) in t.iter().enumerate() {
                if c < k {
                    for p in self.d_ptr[c]..self.d_ptr[c + 1] {
                        v[self.d_prow[p] as usize] += self.d_val[p] * zc;
                    }
                } else {
                    v[self.mods[c - k].0 as usize] += zc;
                }
            }
            self.chol.solve_permuted(v);
            for (xv, vv) in x.iter_mut().zip(v.iter()) {
                *xv -= vv;
            }
        }
        for (kk, xv) in x.iter().enumerate() {
            rhs[self.chol.perm(kk)] = *xv;
        }
    }

    /// `t ← K⁻¹·t` by the block factorization of K.
    fn solve_k(&self, t: &mut [f64]) {
        let k = self.dcols.len();
        let q = self.mods.len();
        let (t1, t2) = t.split_at_mut(k);
        lower_solve(&self.l11, k, t1);
        for a in 0..q {
            let ya = &self.y[a * k..(a + 1) * k];
            t2[a] -= ya.iter().zip(t1.iter()).map(|(p, r)| p * r).sum::<f64>();
        }
        lu_solve(&self.s_lu, q, &self.s_piv, t2);
        for a in 0..q {
            let ya = &self.y[a * k..(a + 1) * k];
            let za = t2[a];
            for (r, p) in t1.iter_mut().zip(ya) {
                *r -= p * za;
            }
        }
        upper_solve(&self.l11, k, t1);
    }

    /// `y = (A·diag(theta)·Aᵀ + δI)·x` on the full A (original rows).
    pub fn matvec(&self, theta: &[f64], delta: f64, x: &[f64], y: &mut [f64]) {
        if !self.split() {
            self.chol.matvec(theta, delta, x, y);
            return;
        }
        for (yv, xv) in y.iter_mut().zip(x) {
            *yv = delta * xv;
        }
        for j in 0..self.a_ptr.len() - 1 {
            let (p0, p1) = (self.a_ptr[j], self.a_ptr[j + 1]);
            let mut s = 0.0;
            for p in p0..p1 {
                s += self.a_val[p] * x[self.a_row[p] as usize];
            }
            if s == 0.0 {
                continue;
            }
            let t = s * theta[j];
            for p in p0..p1 {
                y[self.a_row[p] as usize] += self.a_val[p] * t;
            }
        }
    }
}

/// Solve `L z = b` in place, L lower (column-major, ld = n).
fn lower_solve(l: &[f64], n: usize, b: &mut [f64]) {
    for j in 0..n {
        let col = &l[j * n..(j + 1) * n];
        let x = b[j] / col[j];
        b[j] = x;
        if x != 0.0 {
            for r in j + 1..n {
                b[r] -= col[r] * x;
            }
        }
    }
}

/// Solve `Lᵀ z = b` in place.
fn upper_solve(l: &[f64], n: usize, b: &mut [f64]) {
    for j in (0..n).rev() {
        let col = &l[j * n..(j + 1) * n];
        let mut acc = b[j];
        for r in j + 1..n {
            acc -= col[r] * b[r];
        }
        b[j] = acc / col[j];
    }
}

/// In-place LU with partial pivoting (column-major, n × n); a zero pivot
/// column is replaced by a huge pivot (that coordinate of the solve → 0).
fn lu_factor(a: &mut [f64], n: usize) -> Vec<usize> {
    let mut piv = vec![0usize; n];
    let scale = a.iter().fold(0.0f64, |s, v| s.max(v.abs())).max(1e-300);
    for j in 0..n {
        let (mut p, mut best) = (j, 0.0f64);
        for r in j..n {
            if a[r + j * n].abs() > best {
                best = a[r + j * n].abs();
                p = r;
            }
        }
        piv[j] = p;
        if p != j {
            for c in 0..n {
                a.swap(j + c * n, p + c * n);
            }
        }
        if best <= 1e-14 * scale {
            a[j + j * n] = 1e64;
            for r in j + 1..n {
                a[r + j * n] = 0.0;
            }
            continue;
        }
        let d = a[j + j * n];
        for r in j + 1..n {
            a[r + j * n] /= d;
        }
        for c in j + 1..n {
            let f = a[j + c * n];
            if f != 0.0 {
                for r in j + 1..n {
                    a[r + c * n] -= a[r + j * n] * f;
                }
            }
        }
    }
    piv
}

fn lu_solve(a: &[f64], n: usize, piv: &[usize], b: &mut [f64]) {
    for j in 0..n {
        b.swap(j, piv[j]);
    }
    for j in 0..n {
        let x = b[j];
        for r in j + 1..n {
            b[r] -= a[r + j * n] * x;
        }
    }
    for j in (0..n).rev() {
        let x = b[j] / a[j + j * n];
        b[j] = x;
        for r in 0..j {
            b[r] -= a[r + j * n] * x;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::XorShift64;

    /// Random sparse A plus `nd` dense columns; `nd` rows only in dense ones.
    fn with_dense(m: usize, n: usize, nd: usize, seed: u64) -> CscMatrix {
        let mut rng = XorShift64::new(seed);
        let mut trip = Vec::new();
        for j in 0..n {
            for _ in 0..3 {
                let r = (rng.uniform(0.0, (m - nd) as f64) as usize).min(m - nd - 1);
                trip.push((r, j, rng.uniform(-2.0, 2.0)));
            }
        }
        for d in 0..nd {
            for i in 0..m {
                if rng.uniform(0.0, 1.0) < 0.8 || i >= m - nd {
                    trip.push((i, n + d, rng.uniform(-1.0, 1.0)));
                }
            }
        }
        for i in 0..m - nd {
            trip.push((i, n + nd + i, 1.0));
        }
        CscMatrix::from_triplets(m, n + nd + m - nd, &trip)
    }

    fn residual(ne: &NormalEq, theta: &[f64], delta: f64, m: usize) -> f64 {
        let b: Vec<f64> = (0..m).map(|i| (i as f64 * 0.71).cos()).collect();
        let mut x = b.clone();
        let mut work = Vec::new();
        ne.solve(&mut x, &mut work);
        let mut y = vec![0.0; m];
        ne.matvec(theta, delta, &x, &mut y);
        let bn = b.iter().fold(0.0f64, |s, v| s.max(v.abs()));
        y.iter()
            .zip(&b)
            .fold(0.0f64, |s, (p, q)| s.max((p - q).abs()))
            / bn
    }

    #[test]
    fn dense_split_solves_full_system() {
        for (seed, (m, n, nd)) in [(1u64, (60, 80, 2)), (2, (150, 200, 5)), (3, (40, 30, 1))] {
            let a = with_dense(m, n, nd, seed);
            let mut rng = XorShift64::new(seed + 99);
            let theta: Vec<f64> = (0..a.ncol())
                .map(|_| 10f64.powf(rng.uniform(-3.0, 3.0)))
                .collect();
            let cols: Vec<usize> = (n..n + nd).collect();
            let choice = DenseChoice {
                cols,
                ..Default::default()
            };
            let mut ne = NormalEq::analyse_with(&a, usize::MAX, choice).unwrap();
            ne.factor(&theta, 1e-10);
            assert!(!ne.aborted);
            // The rows touched only by dense columns force modified pivots.
            assert!(ne.stats.dense_mods >= nd, "mods {}", ne.stats.dense_mods);
            let r = residual(&ne, &theta, 1e-10, m);
            assert!(r <= 1e-7, "seed {seed}: residual {r}");
            let mut full = NormalEq::analyse_with(&a, usize::MAX, DenseChoice::default()).unwrap();
            full.factor(&theta, 1e-10);
            let r0 = residual(&full, &theta, 1e-10, m);
            assert!(r0 <= 1e-7, "unsplit residual {r0}");
        }
    }

    #[test]
    fn cost_model_picks_the_dense_column() {
        let a = with_dense(400, 600, 1, 5);
        let ch = choose_dense(&a, usize::MAX);
        assert_eq!(ch.cols, vec![600], "tried {:?}", ch.tried);
    }

    #[test]
    fn lu_solves_indefinite() {
        let n = 5;
        let mut rng = XorShift64::new(4);
        let a0: Vec<f64> = (0..n * n).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let mut a = a0.clone();
        let piv = lu_factor(&mut a, n);
        let b: Vec<f64> = (0..n).map(|i| i as f64 - 2.0).collect();
        let mut x = b.clone();
        lu_solve(&a, n, &piv, &mut x);
        for i in 0..n {
            let s: f64 = (0..n).map(|j| a0[i + j * n] * x[j]).sum();
            assert!((s - b[i]).abs() < 1e-10);
        }
    }
}
