//! Infeasible-start primal-dual IPM iteration on the scaled internal LP:
//! Mehrotra predictor-corrector + Gondzio multiple centrality correctors,
//! normal equations by [`super::chol::NormalChol`].
//!
//! Provenance:
//! * Newton system and its normal-equations reduction: Wright, *Primal-Dual
//!   Interior-Point Methods* (SIAM 1997) ch. 10–11; residual naming r1..r7
//!   follows HiGHS HiPO `Iterate.cpp` (MIT, read as reference text);
//! * predictor-corrector, σ = (μ_aff/μ)³, step heuristic with γ_f = 0.9
//!   and the starting point: Mehrotra, SIAM J. Optim. 2(4) (1992) §6–7;
//!   Lustig, Marsten & Shanno, ORSA J. Comput. 6 (1994);
//! * multiple centrality correctors and the corrector-count rule from the
//!   factor/solve cost ratio: Gondzio, COAP 6 (1996) §3–4; weighted
//!   acceptance: Colombo & Gondzio, COAP 41 (2008);
//! * static regularization ρ = 1e-12, δ = 1e-10: Altman & Gondzio, OMS 11–12
//!   (1999); HiPO `Parameters.h` (MIT);
//! * infeasibility certificates: Farkas' lemma for `Ax = b, l ≤ x ≤ u` and
//!   the recession-cone ray for dual infeasibility; Andersen et al. (1996)
//!   §5 for the divergence heuristics that trigger the check.

use std::time::Instant;

use super::chol::CholError;
use super::model::IpmLp;
use super::normal::NormalEq;

/// Static primal (ρ) / dual (δ) regularization defaults (HiPO
/// Parameters.h `Regularisation { primal 1e-12, dual 1e-10 }`).
const RHO_DEFAULT: f64 = 1e-12;
const DELTA_DEFAULT: f64 = 1e-10;
/// Fallback regularization for the single retry after NoProgress /
/// Imprecise (needed by milo-v12-6-r2-40-1).
const RETRY_REG: f64 = 1e-6;
/// Structural-dependence probe threshold on the Θ = I starting factor
/// (1e-10 measured as the threshold that rescues qap15).
const DEP_REL: f64 = 1e-10;
/// Adaptive regularization: on Θ-cancellation pivots, raise ρ and δ by this
/// factor and refactor, up to `REG_CAP` (the retry value).
const REG_STEP: f64 = 100.0;
const REG_CAP: f64 = 1e-6;
fn env_f64(key: &str, dflt: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(dflt)
}
const GAMMA_F: f64 = 0.9;
const REFINE_STEPS: usize = 4;
const MCC_BETA: f64 = 0.1;
const MCC_ALPHA_INC: f64 = 0.1;
const MCC_MIN_GAIN: f64 = 0.01;
/// Bounds below which a certificate's infinite-side leakage is negligible:
/// a Farkas/ray certificate is accepted only if it holds for every x with
/// ‖x‖∞ ≤ CERT_RADIUS (unscaled internal space).
const CERT_RADIUS: f64 = 1e9;

#[derive(Debug, Clone)]
pub struct IpmConfig {
    /// pinf, dinf and gap tolerance on the unscaled internal LP.
    pub tol: f64,
    pub max_iter: usize,
    /// Original-LP verification tolerance.
    pub verify_tol: f64,
    pub deadline: Option<Instant>,
    /// Override the Gondzio corrector count (None = cost-ratio rule).
    pub max_correctors: Option<usize>,
    /// Memory guard: max entries of the normal-matrix pattern / factor.
    pub max_chol_entries: usize,
    /// Primal regularization ρ (Θ⁻¹ += ρ) and dual regularization δ
    /// (M += δI).
    pub rho: f64,
    pub delta: f64,
    /// When the first run ends NoProgress / Imprecise, retry once with
    /// ρ = δ = this value (None disables).
    pub retry_reg: Option<f64>,
    /// Tiny-pivot threshold of the structural-dependence probe on the
    /// starting factor (0 disables; `RUSOLVE_IPM_DEPREL`).
    pub dep_rel: f64,
    /// Raise ρ, δ when a factorization hits Θ-cancellation pivots
    /// (`RUSOLVE_IPM_ADAPT`, 0 = off).
    pub adapt_reg: bool,
    pub log: bool,
}

impl Default for IpmConfig {
    fn default() -> Self {
        Self {
            tol: 1e-8,
            max_iter: 200,
            verify_tol: 1e-6,
            deadline: None,
            max_correctors: None,
            max_chol_entries: 300_000_000,
            rho: env_f64("RUSOLVE_IPM_RHO", RHO_DEFAULT),
            delta: env_f64("RUSOLVE_IPM_DELTA", DELTA_DEFAULT),
            retry_reg: match env_f64("RUSOLVE_IPM_RETRY_REG", RETRY_REG) {
                v if v > 0.0 => Some(v),
                _ => None,
            },
            dep_rel: env_f64("RUSOLVE_IPM_DEPREL", DEP_REL),
            adapt_reg: env_f64("RUSOLVE_IPM_ADAPT", 1.0) != 0.0,
            log: std::env::var_os("RUSOLVE_IPM_LOG").is_some(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IpmStatus {
    Optimal,
    PrimalInfeasible,
    DualInfeasible,
    NoProgress,
    /// Internal tolerances met but the original-LP check failed.
    Imprecise,
    IterationLimit,
    TimeLimit,
    /// The normal-matrix factor exceeds the memory guard.
    TooLarge,
}

#[derive(Debug, Clone, Default)]
pub struct IpmStats {
    pub iterations: usize,
    pub correctors: usize,
    pub max_correctors: usize,
    pub m: usize,
    pub n: usize,
    pub nnz_a: usize,
    pub nnz_m: usize,
    pub nnz_l: usize,
    pub supernodes: usize,
    pub dense_cols: usize,
    pub dense_col_share: f64,
    /// Dense columns split out of M (Schur complement, `super::normal`).
    pub dense_split: usize,
    pub skipped_pivots_max: usize,
    pub n_free: usize,
    pub t_order: f64,
    pub t_symbolic: f64,
    pub t_factor: f64,
    pub t_solve: f64,
    pub t_ipm: f64,
    pub pinf: f64,
    pub dinf: f64,
    pub gap: f64,
    /// The run was repeated with `IpmConfig::retry_reg` regularization
    /// (iterations and times include both runs).
    pub retried: bool,
}

/// Interior point of the scaled internal LP.
pub struct IpmPoint {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub zl: Vec<f64>,
    pub zu: Vec<f64>,
}

struct Dir {
    x: Vec<f64>,
    xl: Vec<f64>,
    xu: Vec<f64>,
    y: Vec<f64>,
    zl: Vec<f64>,
    zu: Vec<f64>,
}

impl Dir {
    fn new(n: usize, m: usize) -> Self {
        Self {
            x: vec![0.0; n],
            xl: vec![0.0; n],
            xu: vec![0.0; n],
            y: vec![0.0; m],
            zl: vec![0.0; n],
            zu: vec![0.0; n],
        }
    }
    fn axpy_from(&mut self, base: &Dir, w: f64, c: &Dir) {
        for (d, (b, v)) in [
            (&mut self.x, (&base.x, &c.x)),
            (&mut self.xl, (&base.xl, &c.xl)),
            (&mut self.xu, (&base.xu, &c.xu)),
            (&mut self.y, (&base.y, &c.y)),
            (&mut self.zl, (&base.zl, &c.zl)),
            (&mut self.zu, (&base.zu, &c.zu)),
        ] {
            for i in 0..d.len() {
                d[i] = b[i] + w * v[i];
            }
        }
    }
    fn copy_from(&mut self, o: &Dir) {
        self.x.copy_from_slice(&o.x);
        self.xl.copy_from_slice(&o.xl);
        self.xu.copy_from_slice(&o.xu);
        self.y.copy_from_slice(&o.y);
        self.zl.copy_from_slice(&o.zl);
        self.zu.copy_from_slice(&o.zu);
    }
    fn finite(&self) -> bool {
        [&self.x, &self.y, &self.zl, &self.zu, &self.xl, &self.xu]
            .iter()
            .all(|v| v.iter().all(|t| t.is_finite()))
    }
}

pub struct IpmCore {
    pub status: IpmStatus,
    pub point: IpmPoint,
    pub stats: IpmStats,
}

struct Ctx<'a> {
    lp: &'a IpmLp,
    hl: Vec<bool>,
    hu: Vec<bool>,
    chol: NormalEq,
    theta: Vec<f64>,
    work: Vec<f64>,
    tmp_m: Vec<f64>,
    tmp_m2: Vec<f64>,
    tmp_n: Vec<f64>,
    r7: Vec<f64>,
    t_factor: f64,
    t_solve: f64,
    bnorm: f64,
    refine_steps: usize,
    rho: f64,
    delta: f64,
    adapt_reg: bool,
    /// Per-iteration decay of a raised regularization back toward (ρ0, δ0)
    /// (`RUSOLVE_IPM_REG_DECAY`; 1 = keep the raise).
    reg_decay: f64,
    rho0: f64,
    delta0: f64,
    log: bool,
    reg_raises: usize,
}

impl Ctx<'_> {
    /// Plain solve with the current factor.
    fn solve_plain(&mut self, rhs: &[f64], dy: &mut [f64]) {
        let t0 = Instant::now();
        dy.copy_from_slice(rhs);
        self.chol.solve(dy, &mut self.work);
        self.t_solve += t0.elapsed().as_secs_f64();
    }

    /// Solve M·dy = rhs (theta set, factor current) with refinement.
    fn solve_m(&mut self, rhs: &[f64], dy: &mut [f64]) {
        let t0 = Instant::now();
        let m = self.lp.m;
        dy.copy_from_slice(rhs);
        self.chol.solve(dy, &mut self.work);
        let rnorm = |v: &[f64]| v.iter().fold(0.0f64, |a, t| a.max(t.abs()));
        let bnorm = rnorm(rhs);
        // residual of the unregularized system
        self.chol.matvec(&self.theta, 0.0, dy, &mut self.tmp_m);
        for i in 0..m {
            self.tmp_m[i] = rhs[i] - self.tmp_m[i];
        }
        let mut res = rnorm(&self.tmp_m);
        for _ in 0..3 {
            if res <= 1e-14 * bnorm || !res.is_finite() {
                break;
            }
            let mut corr = std::mem::take(&mut self.tmp_m2);
            corr.clear();
            corr.extend_from_slice(&self.tmp_m);
            self.chol.solve(&mut corr, &mut self.work);
            let mut trial: Vec<f64> = dy.iter().zip(&corr).map(|(a, b)| a + b).collect();
            self.chol.matvec(&self.theta, 0.0, &trial, &mut corr);
            let mut rnew = 0.0f64;
            for i in 0..m {
                corr[i] = rhs[i] - corr[i];
                rnew = rnew.max(corr[i].abs());
            }
            if rnew < res {
                dy.copy_from_slice(&trial);
                self.tmp_m.copy_from_slice(&corr);
                res = rnew;
                trial.clear();
                self.tmp_m2 = corr;
            } else {
                self.tmp_m2 = corr;
                break;
            }
        }
        self.t_solve += t0.elapsed().as_secs_f64();
    }

    /// Newton direction for residuals r1..r6 at iterate `it`.
    #[allow(clippy::too_many_arguments)]
    fn newton(
        &mut self,
        it: &Iter,
        r1: &[f64],
        r2: &[f64],
        r3: &[f64],
        r4: &[f64],
        r5: &[f64],
        r6: &[f64],
        d: &mut Dir,
    ) {
        let lp = self.lp;
        let (m, n) = (lp.m, lp.n);
        for j in 0..n {
            let mut v = r4[j];
            if self.hl[j] {
                v -= (r5[j] + it.zl[j] * r2[j]) / it.xl[j];
            }
            if self.hu[j] {
                v += (r6[j] - it.zu[j] * r3[j]) / it.xu[j];
            }
            self.r7[j] = v;
            self.tmp_n[j] = self.theta[j] * v;
        }
        let mut rhs = vec![0.0; m];
        lp.a.matvec(&self.tmp_n, &mut rhs);
        for i in 0..m {
            rhs[i] += r1[i];
        }
        self.solve_plain(&rhs, &mut d.y);
        lp.a.matvec_transpose(&d.y, &mut self.tmp_n);
        for j in 0..n {
            d.x[j] = self.theta[j] * (self.tmp_n[j] - self.r7[j]);
        }
        // Iterative refinement on the primal Newton equation A·Δx = r1
        // (the block the normal-equations cancellation `r1 + AΘr7` loses):
        // e = r1 − AΔx, M·c = e, Δy += c, Δx += Θ·Aᵀc — keeps
        // Δx = Θ(AᵀΔy − r7), hence the dual block, exact by construction.
        // Accept a step only if ‖e‖∞ decreases (Wright 1997 ch. 11;
        // HiPO `Refine.cpp` refines the full Newton system, MIT).
        let r1n = inf_norm(r1);
        let etol = 1e-10 * r1n + 1e-14 * (1.0 + self.bnorm);
        let mut e = vec![0.0; m];
        lp.a.matvec(&d.x, &mut e);
        for i in 0..m {
            e[i] = r1[i] - e[i];
        }
        let mut en = inf_norm(&e);
        let mut c = vec![0.0; m];
        let mut dxt = vec![0.0; n];
        let mut et = vec![0.0; m];
        for _ in 0..REFINE_STEPS {
            if en <= etol || !en.is_finite() {
                break;
            }
            self.solve_plain(&e, &mut c);
            lp.a.matvec_transpose(&c, &mut self.tmp_n);
            for j in 0..n {
                dxt[j] = d.x[j] + self.theta[j] * self.tmp_n[j];
            }
            lp.a.matvec(&dxt, &mut et);
            for i in 0..m {
                et[i] = r1[i] - et[i];
            }
            let ent = inf_norm(&et);
            if ent.is_nan() || ent >= en {
                break;
            }
            for i in 0..m {
                d.y[i] += c[i];
            }
            std::mem::swap(&mut d.x, &mut dxt);
            std::mem::swap(&mut e, &mut et);
            en = ent;
            self.refine_steps += 1;
        }
        for j in 0..n {
            let dx = d.x[j];
            if self.hl[j] {
                d.xl[j] = dx - r2[j];
                d.zl[j] = (r5[j] - it.zl[j] * d.xl[j]) / it.xl[j];
            } else {
                d.xl[j] = 0.0;
                d.zl[j] = 0.0;
            }
            if self.hu[j] {
                d.xu[j] = r3[j] - dx;
                d.zu[j] = (r6[j] - it.zu[j] * d.xu[j]) / it.xu[j];
            } else {
                d.xu[j] = 0.0;
                d.zu[j] = 0.0;
            }
        }
    }

    /// Θ, factor; with `adapt_reg`, a factor that hit Θ-cancellation pivots
    /// (tiny pivots outside the structural-dependence set) is redone with
    /// ρ = δ = max(ρ, δ)·`REG_STEP` (at most `REG_CAP`) and the raise is kept for the
    /// rest of the run. This replaces skipping such pivots — which makes the
    /// direction inconsistent with r1 — by solving a nearby
    /// regularized system whose error the r1 refinement removes. Sources:
    /// Altman & Gondzio, "Regularized symmetric indefinite systems in
    /// interior point methods for linear and quadratic optimization", OMS
    /// 11–12 (1999) (dynamic regularization); Saunders & Tomlin (1996);
    /// Friedlander & Orban, Math. Prog. Comp. 4 (2012) (primal-dual
    /// regularization keeps the IPM convergent).
    fn factor(&mut self, it: &Iter) {
        let t0 = Instant::now();
        if self.adapt_reg && self.reg_decay > 1.0 && self.rho > self.rho0 {
            self.rho = (self.rho / self.reg_decay).max(self.rho0);
            self.delta = (self.delta / self.reg_decay).max(self.delta0);
        }
        loop {
            for j in 0..self.lp.n {
                let mut d = self.rho;
                if self.hl[j] {
                    d += it.zl[j] / it.xl[j];
                }
                if self.hu[j] {
                    d += it.zu[j] / it.xu[j];
                }
                self.theta[j] = 1.0 / d;
            }
            self.chol.factor(&self.theta, self.delta);
            let cancel = self.chol.stats.cancel_pivots;
            if !self.adapt_reg
                || self.chol.aborted
                || cancel == 0
                || self.rho.min(self.delta) >= REG_CAP
            {
                break;
            }
            let r = (self.rho.max(self.delta) * REG_STEP).min(REG_CAP);
            self.rho = r;
            self.delta = r;
            self.reg_raises += 1;
            if self.log {
                eprintln!(
                    "        {cancel} cancellation pivots: regularization raised to rho {:.0e} delta {:.0e}",
                    self.rho, self.delta
                );
            }
        }
        self.t_factor += t0.elapsed().as_secs_f64();
    }
}

struct Iter {
    x: Vec<f64>,
    xl: Vec<f64>,
    xu: Vec<f64>,
    y: Vec<f64>,
    zl: Vec<f64>,
    zu: Vec<f64>,
}

fn max_step(v: &[f64], dv: &[f64], mask: &[bool]) -> f64 {
    let mut a = f64::INFINITY;
    for j in 0..v.len() {
        if mask[j] && dv[j] < 0.0 {
            a = a.min(-v[j] / dv[j]);
        }
    }
    a
}

fn steps(it: &Iter, d: &Dir, hl: &[bool], hu: &[bool]) -> (f64, f64) {
    let ap = max_step(&it.xl, &d.xl, hl).min(max_step(&it.xu, &d.xu, hu));
    let ad = max_step(&it.zl, &d.zl, hl).min(max_step(&it.zu, &d.zu, hu));
    (ap, ad)
}

fn compl_at(it: &Iter, d: &Dir, ap: f64, ad: f64, hl: &[bool], hu: &[bool]) -> f64 {
    let mut s = 0.0;
    for j in 0..it.x.len() {
        if hl[j] {
            s += (it.xl[j] + ap * d.xl[j]) * (it.zl[j] + ad * d.zl[j]);
        }
        if hu[j] {
            s += (it.xu[j] + ap * d.xu[j]) * (it.zu[j] + ad * d.zu[j]);
        }
    }
    s
}

/// Mehrotra's step-length heuristic (1992 §6): move the blocking variable
/// so its complementarity product equals μ_full/γ_a, never below γ_f·α_max.
fn mehrotra_steps(it: &Iter, d: &Dir, hl: &[bool], hu: &[bool], nb: usize) -> (f64, f64) {
    let (apmax, admax) = steps(it, d, hl, hu);
    let (apm, adm) = (apmax.min(1.0), admax.min(1.0));
    let mu_full = compl_at(it, d, apm, adm, hl, hu) / nb.max(1) as f64;
    let mu_target = mu_full * (1.0 - GAMMA_F);
    let mut ap = 1.0f64;
    if apmax < 1.0 {
        // blocking primal index
        let mut best = (f64::INFINITY, 0usize, false);
        for j in 0..it.x.len() {
            if hl[j] && d.xl[j] < 0.0 {
                let a = -it.xl[j] / d.xl[j];
                if a < best.0 {
                    best = (a, j, false);
                }
            }
            if hu[j] && d.xu[j] < 0.0 {
                let a = -it.xu[j] / d.xu[j];
                if a < best.0 {
                    best = (a, j, true);
                }
            }
        }
        let j = best.1;
        let (v, dv, z, dz) = if best.2 {
            (it.xu[j], d.xu[j], it.zu[j], d.zu[j])
        } else {
            (it.xl[j], d.xl[j], it.zl[j], d.zl[j])
        };
        let zf = z + adm * dz;
        let a = if zf > 0.0 {
            (mu_target / zf - v) / dv
        } else {
            0.0
        };
        ap = a.max(GAMMA_F * apmax).min(apmax * (1.0 - 1e-4));
    }
    let mut ad = 1.0f64;
    if admax < 1.0 {
        let mut best = (f64::INFINITY, 0usize, false);
        for j in 0..it.x.len() {
            if hl[j] && d.zl[j] < 0.0 {
                let a = -it.zl[j] / d.zl[j];
                if a < best.0 {
                    best = (a, j, false);
                }
            }
            if hu[j] && d.zu[j] < 0.0 {
                let a = -it.zu[j] / d.zu[j];
                if a < best.0 {
                    best = (a, j, true);
                }
            }
        }
        let j = best.1;
        let (z, dz, v, dv) = if best.2 {
            (it.zu[j], d.zu[j], it.xu[j], d.xu[j])
        } else {
            (it.zl[j], d.zl[j], it.xl[j], d.xl[j])
        };
        let xf = v + apm * dv;
        let a = if xf > 0.0 {
            (mu_target / xf - z) / dz
        } else {
            0.0
        };
        ad = a.max(GAMMA_F * admax).min(admax * (1.0 - 1e-4));
    }
    (ap.clamp(0.0, 1.0), ad.clamp(0.0, 1.0))
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |a, t| a.max(t.abs()))
}

struct Measures {
    pinf: f64,
    dinf: f64,
    gap: f64,
    pobj: f64,
    dobj: f64,
    mu: f64,
}

/// Run the IPM on `lp` (artificial bounds of free columns may be widened).
pub fn run_ipm(lp: &mut IpmLp, cfg: &IpmConfig) -> IpmCore {
    let t_start = Instant::now();
    let (m, n) = (lp.m, lp.n);
    let mut stats = IpmStats {
        m,
        n,
        nnz_a: lp.a.nnz(),
        n_free: lp.n_free,
        ..Default::default()
    };
    let empty = |status| IpmCore {
        status,
        point: IpmPoint {
            x: vec![0.0; n],
            y: vec![0.0; m],
            zl: vec![0.0; n],
            zu: vec![0.0; n],
        },
        stats: IpmStats::default(),
    };
    let mut chol = match NormalEq::analyse(&lp.a, cfg.max_chol_entries) {
        Ok(c) => c,
        Err(CholError::TooLarge(msg)) => {
            if cfg.log {
                eprintln!("rusolve ipm: {msg}");
            }
            let mut r = empty(IpmStatus::TooLarge);
            r.stats = stats;
            return r;
        }
    };
    chol.deadline = cfg.deadline;
    let cs = &chol.stats;
    stats.nnz_m = cs.nnz_m;
    stats.nnz_l = cs.nnz_l;
    stats.supernodes = cs.supernodes;
    stats.dense_cols = cs.dense_cols;
    stats.dense_col_share = cs.dense_col_share;
    stats.dense_split = cs.dense_split;
    stats.t_order = cs.t_order;
    stats.t_symbolic = cs.t_symbolic;
    // Gondzio (1996) §4: correctors pay off when a factorization costs many
    // solves; K from the factor/solve flop ratio, clamped to [1, 10].
    let ratio = if cs.dense_split > 0 {
        chol.cost_ratio()
    } else {
        cs.flops / (4.0 * cs.nnz_l.max(1) as f64)
    };
    let kmax = cfg
        .max_correctors
        .unwrap_or_else(|| (1.0 + ((ratio - 10.0) / 20.0).floor().max(0.0)).min(10.0) as usize)
        .min(10);
    stats.max_correctors = kmax;
    if cfg.log {
        eprintln!(
            "rusolve ipm: m {} n {} nnz(A) {} nnz(M) {} nnz(L) {} supernodes {} dense cols {} ({:.1}% of M work) order {:.3}s symb {:.3}s correctors<= {}",
            m,
            n,
            lp.a.nnz(),
            cs.nnz_m,
            cs.nnz_l,
            cs.supernodes,
            cs.dense_cols,
            100.0 * cs.dense_col_share,
            cs.t_order,
            cs.t_symbolic,
            kmax
        );
        eprintln!(
            "rusolve ipm: chol fundamental supernodes {} relaxed {} stored nnz(L) {} flops {:.3e} stored {:.3e} threads {} subtrees {}",
            cs.supernodes_fund,
            cs.supernodes,
            cs.nnz_l_stored,
            cs.flops,
            cs.flops_stored,
            cs.threads,
            cs.subtrees
        );
        if !chol.choice.tried.is_empty() {
            eprintln!(
                "rusolve ipm: dense-column split k {} (candidates k:cost {:?}) select {:.3}s",
                cs.dense_split,
                chol.choice
                    .tried
                    .iter()
                    .map(|(k, c)| format!("{k}:{c:.2e}"))
                    .collect::<Vec<_>>(),
                chol.choice.t_select
            );
        }
    }
    let hl: Vec<bool> = lp.l.iter().map(|v| v.is_finite()).collect();
    let hu: Vec<bool> = lp.u.iter().map(|v| v.is_finite()).collect();
    let nb = lp.num_bounds();
    let mut cx = Ctx {
        lp,
        hl: hl.clone(),
        hu: hu.clone(),
        chol,
        theta: vec![1.0; n],
        work: Vec::new(),
        tmp_m: vec![0.0; m],
        tmp_m2: vec![0.0; m],
        tmp_n: vec![0.0; n],
        r7: vec![0.0; n],
        t_factor: 0.0,
        t_solve: 0.0,
        bnorm: inf_norm(&lp.b),
        refine_steps: 0,
        rho: cfg.rho,
        delta: cfg.delta,
        adapt_reg: cfg.adapt_reg,
        reg_decay: env_f64("RUSOLVE_IPM_REG_DECAY", 1.0),
        rho0: cfg.rho,
        delta0: cfg.delta,
        log: cfg.log,
        reg_raises: 0,
    };

    // ---- starting point (Mehrotra 1992 §7, HiPO startingPoint) ----
    let mut it = Iter {
        x: vec![0.0; n],
        xl: vec![0.0; n],
        xu: vec![0.0; n],
        y: vec![0.0; m],
        zl: vec![0.0; n],
        zu: vec![0.0; n],
    };
    {
        let lp = cx.lp;
        let t0 = Instant::now();
        let ndep = cx.chol.factor_detect(&cx.theta, cx.delta, cfg.dep_rel);
        cx.t_factor += t0.elapsed().as_secs_f64();
        if cfg.log && ndep > 0 {
            eprintln!("rusolve ipm: {ndep} structurally dependent rows (always skipped)");
        }
        if cx.chol.aborted {
            let mut r = empty(IpmStatus::TimeLimit);
            stats.t_factor = cx.t_factor;
            stats.t_ipm = t_start.elapsed().as_secs_f64();
            r.stats = stats;
            return r;
        }
        for j in 0..n {
            it.x[j] = 0.0f64.clamp(lp.l[j], lp.u[j]);
        }
        let mut r = vec![0.0; m];
        lp.a.matvec(&it.x, &mut r);
        for i in 0..m {
            r[i] = lp.b[i] - r[i];
        }
        let mut w = vec![0.0; m];
        cx.solve_m(&r, &mut w);
        let mut dx = vec![0.0; n];
        lp.a.matvec_transpose(&w, &mut dx);
        for j in 0..n {
            it.x[j] += dx[j];
        }
        let mut vmin = f64::INFINITY;
        for j in 0..n {
            if hl[j] {
                it.xl[j] = it.x[j] - lp.l[j];
                vmin = vmin.min(it.xl[j]);
            }
            if hu[j] {
                it.xu[j] = lp.u[j] - it.x[j];
                vmin = vmin.min(it.xu[j]);
            }
        }
        let shift = 1.0 + (-1.5 * vmin).max(0.0);
        for j in 0..n {
            if hl[j] {
                it.xl[j] += shift;
            }
            if hu[j] {
                it.xu[j] += shift;
            }
        }
        // y = (AAᵀ)⁻¹ A c, z = c − Aᵀy
        let mut ac = vec![0.0; m];
        lp.a.matvec(&lp.c, &mut ac);
        let mut y = vec![0.0; m];
        cx.solve_m(&ac, &mut y);
        it.y = y;
        let mut aty = vec![0.0; n];
        lp.a.matvec_transpose(&it.y, &mut aty);
        let mut zmin = f64::INFINITY;
        for j in 0..n {
            let z = lp.c[j] - aty[j];
            match (hl[j], hu[j]) {
                (true, true) => {
                    it.zl[j] = z.max(0.0);
                    it.zu[j] = (-z).max(0.0);
                }
                (true, false) => it.zl[j] = z,
                (false, true) => it.zu[j] = -z,
                (false, false) => {}
            }
            if hl[j] {
                zmin = zmin.min(it.zl[j]);
            }
            if hu[j] {
                zmin = zmin.min(it.zu[j]);
            }
        }
        let zshift = 1.0 + (-1.5 * zmin).max(0.0);
        let (mut xz, mut xs, mut zs) = (0.0, 0.0, 0.0);
        for j in 0..n {
            if hl[j] {
                it.zl[j] += zshift;
                xz += it.xl[j] * it.zl[j];
                xs += it.xl[j];
                zs += it.zl[j];
            }
            if hu[j] {
                it.zu[j] += zshift;
                xz += it.xu[j] * it.zu[j];
                xs += it.xu[j];
                zs += it.zu[j];
            }
        }
        if nb > 0 {
            let dxs = 0.5 * xz / zs;
            let dzs = 0.5 * xz / xs;
            for j in 0..n {
                if hl[j] {
                    it.xl[j] += dxs;
                    it.zl[j] += dzs;
                }
                if hu[j] {
                    it.xu[j] += dxs;
                    it.zu[j] += dzs;
                }
            }
        }
    }

    // Unscaled norms for the termination measures (design §6).
    let lp = cx.lp;
    let bnorm = (0..m)
        .map(|i| (lp.b[i] / lp.rscale[i]).abs())
        .fold(0.0f64, f64::max);
    let mut bound_norm = 0.0f64;
    for j in 0..n {
        if lp.art[j] {
            continue;
        }
        if hl[j] {
            bound_norm = bound_norm.max((lp.l[j] * lp.cscale[j]).abs());
        }
        if hu[j] {
            bound_norm = bound_norm.max((lp.u[j] * lp.cscale[j]).abs());
        }
    }
    let cnorm = (0..n)
        .map(|j| (lp.c[j] / lp.cscale[j]).abs())
        .fold(0.0f64, f64::max);
    let pden = 1.0 + bnorm.max(bound_norm);
    let dden = 1.0 + cnorm;

    let mut r1 = vec![0.0; m];
    let mut r2 = vec![0.0; n];
    let mut r3 = vec![0.0; n];
    let mut r4 = vec![0.0; n];
    let mut r5 = vec![0.0; n];
    let mut r6 = vec![0.0; n];
    let zero_m = vec![0.0; m];
    let zero_n = vec![0.0; n];
    let mut aff = Dir::new(n, m);
    let mut dir = Dir::new(n, m);
    let mut cor = Dir::new(n, m);
    let mut trial = Dir::new(n, m);
    let mut best_trial = Dir::new(n, m);
    let mut aty = vec![0.0; n];
    let mut l_cur = lp.l.clone();
    let mut u_cur = lp.u.clone();
    let mut bad_iters = 0usize;
    let mut mu_best = f64::INFINITY;
    let mut status;
    let mut iter = 0usize;
    let mut last;

    loop {
        // ---- residuals and termination ----
        lp.a.matvec(&it.x, &mut r1);
        for i in 0..m {
            r1[i] = lp.b[i] - r1[i];
        }
        lp.a.matvec_transpose(&it.y, &mut aty);
        let (mut p2, mut d4) = (
            inf_norm(
                &r1.iter()
                    .zip(&lp.rscale)
                    .map(|(r, s)| r / s)
                    .collect::<Vec<_>>(),
            ),
            0.0f64,
        );
        let mut compl = 0.0;
        let mut pobj = 0.0;
        let mut dobj: f64 = lp.b.iter().zip(&it.y).map(|(b, y)| b * y).sum();
        for j in 0..n {
            pobj += lp.c[j] * it.x[j];
            r4[j] = lp.c[j] - aty[j] - it.zl[j] + it.zu[j];
            d4 = d4.max((r4[j] / lp.cscale[j]).abs());
            if hl[j] {
                r2[j] = l_cur[j] - it.x[j] + it.xl[j];
                p2 = p2.max((r2[j] * lp.cscale[j]).abs());
                compl += it.xl[j] * it.zl[j];
                dobj += l_cur[j] * it.zl[j];
            } else {
                r2[j] = 0.0;
            }
            if hu[j] {
                r3[j] = u_cur[j] - it.x[j] - it.xu[j];
                p2 = p2.max((r3[j] * lp.cscale[j]).abs());
                compl += it.xu[j] * it.zu[j];
                dobj -= u_cur[j] * it.zu[j];
            } else {
                r3[j] = 0.0;
            }
        }
        let mu = compl / nb.max(1) as f64;
        let meas = Measures {
            pinf: p2 / pden,
            dinf: d4 / dden,
            gap: (pobj - dobj).abs() / (1.0 + 0.5 * (pobj + dobj).abs()),
            pobj,
            dobj,
            mu,
        };
        if cfg.log {
            eprintln!(
                "ipm {:>3} pobj {:+.10e} dobj {:+.10e} pinf {:.2e} dinf {:.2e} gap {:.2e} mu {:.2e}",
                iter, meas.pobj, meas.dobj, meas.pinf, meas.dinf, meas.gap, meas.mu
            );
        }
        let finite = meas.pobj.is_finite() && meas.dobj.is_finite() && meas.mu.is_finite();
        last = meas;
        if !finite {
            status = IpmStatus::NoProgress;
            break;
        }
        if last.pinf <= cfg.tol && last.dinf <= cfg.tol && last.gap <= cfg.tol {
            status = IpmStatus::Optimal;
            break;
        }
        // ---- divergence → certificates (design §6) ----
        if iter >= 3 {
            let ynorm = inf_norm(&it.y);
            let xnorm = inf_norm(&it.x);
            let pdiv = last.dobj > 1e6 * last.pobj.abs().max(1.0) || mu > 1e3 * mu_best;
            let ddiv = last.pobj < -1e6 * last.dobj.abs().max(1.0) || mu > 1e3 * mu_best;
            if (pdiv || ynorm > 1e8) && farkas_ok(lp, &it.y, &hl, &hu) {
                status = IpmStatus::PrimalInfeasible;
                break;
            }
            if (ddiv || xnorm > 1e8)
                && (ray_ok(lp, &it.x, &hl, &hu) || ray_ok(lp, &dir.x, &hl, &hu))
            {
                status = IpmStatus::DualInfeasible;
                break;
            }
        }
        mu_best = mu_best.min(mu);
        if iter >= cfg.max_iter {
            status = IpmStatus::IterationLimit;
            break;
        }
        if cfg.deadline.is_some_and(|d| Instant::now() >= d) {
            status = IpmStatus::TimeLimit;
            break;
        }
        iter += 1;

        // ---- factorization ----
        cx.factor(&it);
        if cx.chol.aborted {
            status = IpmStatus::TimeLimit;
            break;
        }
        stats.skipped_pivots_max = stats.skipped_pivots_max.max(cx.chol.stats.skipped_pivots);

        // ---- predictor ----
        for j in 0..n {
            r5[j] = if hl[j] { -it.xl[j] * it.zl[j] } else { 0.0 };
            r6[j] = if hu[j] { -it.xu[j] * it.zu[j] } else { 0.0 };
        }
        cx.newton(&it, &r1, &r2, &r3, &r4, &r5, &r6, &mut aff);
        if !aff.finite() {
            status = IpmStatus::NoProgress;
            break;
        }
        let (apa, ada) = steps(&it, &aff, &hl, &hu);
        let (apa, ada) = (apa.min(1.0), ada.min(1.0));
        let mu_aff = compl_at(&it, &aff, apa, ada, &hl, &hu) / nb.max(1) as f64;
        let sigma = if mu > 0.0 {
            (mu_aff / mu).powi(3).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let smu = sigma * mu;
        // ---- corrector ----
        for j in 0..n {
            if hl[j] {
                r5[j] = smu - it.xl[j] * it.zl[j] - aff.xl[j] * aff.zl[j];
            }
            if hu[j] {
                r6[j] = smu - it.xu[j] * it.zu[j] - aff.xu[j] * aff.zu[j];
            }
        }
        cx.newton(&it, &r1, &r2, &r3, &r4, &r5, &r6, &mut dir);
        if !dir.finite() {
            status = IpmStatus::NoProgress;
            break;
        }
        if cfg.log {
            // Newton accuracy probe: ‖A·Δx − r1‖∞ relative to ‖r1‖∞ + ‖A‖·‖Δx‖.
            let mut adx = vec![0.0; m];
            lp.a.matvec(&dir.x, &mut adx);
            let mut e = 0.0f64;
            for i in 0..m {
                e = e.max((adx[i] - r1[i]).abs());
            }
            eprintln!(
                "        newton |A dx - r1| {:.2e} |r1| {:.2e} |dx| {:.2e} |dy| {:.2e} skipped {}",
                e,
                inf_norm(&r1),
                inf_norm(&dir.x),
                inf_norm(&dir.y),
                cx.chol.stats.skipped_pivots
            );
        }
        // ---- Gondzio multiple centrality correctors ----
        let (mut ap, mut ad) = steps(&it, &dir, &hl, &hu);
        ap = ap.min(1.0);
        ad = ad.min(1.0);
        let mut k = 0;
        while k < kmax && (ap < 0.95 || ad < 0.95) && nb > 0 {
            let tp = (ap + MCC_ALPHA_INC).min(1.0);
            let td = (ad + MCC_ALPHA_INC).min(1.0);
            let lo = MCC_BETA * smu;
            let hi = smu / MCC_BETA;
            let target = |v: f64| {
                if v < lo {
                    lo - v
                } else if v > hi {
                    (hi - v).max(-hi)
                } else {
                    0.0
                }
            };
            for j in 0..n {
                r5[j] = if hl[j] {
                    target((it.xl[j] + tp * dir.xl[j]) * (it.zl[j] + td * dir.zl[j]))
                } else {
                    0.0
                };
                r6[j] = if hu[j] {
                    target((it.xu[j] + tp * dir.xu[j]) * (it.zu[j] + td * dir.zu[j]))
                } else {
                    0.0
                };
            }
            cx.newton(&it, &zero_m, &zero_n, &zero_n, &zero_n, &r5, &r6, &mut cor);
            if !cor.finite() {
                break;
            }
            // 9-point weight search on [αp·αd, 1] (Colombo & Gondzio 2008).
            let wmin = (ap * ad).clamp(0.0, 1.0);
            let mut best = (ap + ad, -1.0f64, ap, ad);
            for s in 0..9 {
                let w = wmin + (1.0 - wmin) * s as f64 / 8.0;
                trial.axpy_from(&dir, w, &cor);
                let (tp2, td2) = steps(&it, &trial, &hl, &hu);
                let (tp2, td2) = (tp2.min(1.0), td2.min(1.0));
                if tp2 + td2 > best.0 {
                    best = (tp2 + td2, w, tp2, td2);
                    best_trial.copy_from(&trial);
                }
            }
            if best.1 < 0.0 || best.0 < ap + ad + MCC_MIN_GAIN {
                break;
            }
            dir.copy_from(&best_trial);
            ap = best.2;
            ad = best.3;
            k += 1;
            stats.correctors += 1;
        }
        // ---- step ----
        let (ap, ad) = mehrotra_steps(&it, &dir, &hl, &hu, nb);
        if ap.min(ad) < 0.05 {
            bad_iters += 1;
        } else {
            bad_iters = 0;
        }
        for j in 0..n {
            it.x[j] += ap * dir.x[j];
            if hl[j] {
                it.xl[j] = (it.xl[j] + ap * dir.xl[j]).max(f64::MIN_POSITIVE);
                it.zl[j] = (it.zl[j] + ad * dir.zl[j]).max(f64::MIN_POSITIVE);
            }
            if hu[j] {
                it.xu[j] = (it.xu[j] + ap * dir.xu[j]).max(f64::MIN_POSITIVE);
                it.zu[j] = (it.zu[j] + ad * dir.zu[j]).max(f64::MIN_POSITIVE);
            }
        }
        for i in 0..m {
            it.y[i] += ad * dir.y[i];
        }
        // Widen artificial boxes the iterate approaches (keeps r2/r3).
        for j in 0..n {
            if lp.art[j] && it.x[j].abs() > 0.5 * u_cur[j] {
                let nb_ = 2.0 * it.x[j].abs();
                let grow = nb_ - u_cur[j];
                it.xl[j] += grow;
                it.xu[j] += grow;
                l_cur[j] = -nb_;
                u_cur[j] = nb_;
            }
        }
        if cfg.log {
            eprintln!(
                "        sigma {:.2e} ap {:.3} ad {:.3} corr {}",
                sigma, ap, ad, k
            );
        }
        if bad_iters >= 5 {
            // Last chance: the stall may be an infeasibility signature.
            if farkas_ok(lp, &it.y, &hl, &hu) {
                status = IpmStatus::PrimalInfeasible;
            } else if ray_ok(lp, &it.x, &hl, &hu) || ray_ok(lp, &dir.x, &hl, &hu) {
                status = IpmStatus::DualInfeasible;
            } else {
                status = IpmStatus::NoProgress;
            }
            break;
        }
    }
    if matches!(status, IpmStatus::IterationLimit | IpmStatus::NoProgress) {
        if farkas_ok(lp, &it.y, &hl, &hu) {
            status = IpmStatus::PrimalInfeasible;
        } else if ray_ok(lp, &it.x, &hl, &hu) || ray_ok(lp, &dir.x, &hl, &hu) {
            status = IpmStatus::DualInfeasible;
        }
    }
    stats.iterations = iter;
    stats.t_factor = cx.t_factor;
    stats.t_solve = cx.t_solve;
    stats.pinf = last.pinf;
    stats.dinf = last.dinf;
    stats.gap = last.gap;
    stats.t_ipm = t_start.elapsed().as_secs_f64();
    if cfg.log {
        eprintln!(
            "rusolve ipm: {:?} after {} iterations ({} correctors) factor {:.3}s solve {:.3}s total {:.3}s skipped pivots <= {}",
            status,
            iter,
            stats.correctors,
            stats.t_factor,
            stats.t_solve,
            stats.t_ipm,
            stats.skipped_pivots_max
        );
        let (ns, np) = cx.chol.pcg_counts();
        if ns > 0 {
            eprintln!(
                "rusolve ipm: split solves {ns} pcg steps {np} regularization raises {} unsplit {}",
                cx.reg_raises, cx.chol.stats.unsplit
            );
        }
    }
    IpmCore {
        status,
        point: IpmPoint {
            x: it.x,
            y: it.y,
            zl: it.zl,
            zu: it.zu,
        },
        stats,
    }
}

/// Farkas certificate for `Ax = b, l ≤ x ≤ u` (artificial bounds count as
/// infinite): `w = Aᵀy`, margin `bᵀy − sup_{l≤x≤u} wᵀx > 0`, with the
/// infinite-side leakage `Σ|w_j|` small enough that no x with
/// ‖x‖∞ ≤ CERT_RADIUS is feasible. Evaluated in the unscaled internal
/// space, y normalized to ‖y‖∞ = 1.
fn farkas_ok(lp: &IpmLp, yh: &[f64], hl: &[bool], hu: &[bool]) -> bool {
    let (m, n) = (lp.m, lp.n);
    if m == 0 {
        return false;
    }
    let y: Vec<f64> = (0..m).map(|i| yh[i] * lp.rscale[i]).collect();
    let ny = inf_norm(&y);
    if ny.is_nan() || ny <= 0.0 || !ny.is_finite() {
        return false;
    }
    // scaled ŷ normalized so that the unscaled y has unit norm
    let yhn: Vec<f64> = yh.iter().map(|v| v / ny).collect();
    let mut w = vec![0.0; n];
    lp.a.matvec_transpose(&yhn, &mut w); // = C·(Aᵀy) per column
    let mut margin: f64 = lp.b.iter().zip(&yhn).map(|(b, y)| b * y).sum();
    let mut mag: f64 = lp.b.iter().zip(&yhn).map(|(b, y)| (b * y).abs()).sum();
    let mut leak = 0.0f64;
    for j in 0..n {
        let wj = w[j]; // ŵ_j · x̂_j = w_j · x_j
        let fl = hl[j] && !lp.art[j];
        let fu = hu[j] && !lp.art[j];
        if wj > 0.0 {
            if fu {
                margin -= wj * lp.u[j];
                mag += (wj * lp.u[j]).abs();
            } else {
                leak += wj / lp.cscale[j];
            }
        } else if wj < 0.0 {
            if fl {
                margin -= wj * lp.l[j];
                mag += (wj * lp.l[j]).abs();
            } else {
                leak += -wj / lp.cscale[j];
            }
        }
    }
    margin > 1e-7 * (1.0 + mag) && leak * CERT_RADIUS <= 0.5 * margin
}

/// Dual-infeasibility certificate: a direction d (from the diverging x)
/// with A·d ≈ 0, d in the recession cone of the true bounds and cᵀd < 0.
fn ray_ok(lp: &IpmLp, xh: &[f64], hl: &[bool], hu: &[bool]) -> bool {
    let (m, n) = (lp.m, lp.n);
    // unscaled direction, normalized
    let d: Vec<f64> = (0..n).map(|j| xh[j] * lp.cscale[j]).collect();
    let nd = inf_norm(&d);
    if nd.is_nan() || nd <= 0.0 || !nd.is_finite() {
        return false;
    }
    // Drop components that are bounded noise relative to the ray.
    let dh: Vec<f64> = (0..n)
        .map(|j| {
            let v = d[j] / nd;
            if v.abs() < 1e-12 { 0.0 } else { xh[j] / nd }
        })
        .collect();
    let mut ad = vec![0.0; m];
    lp.a.matvec(&dh, &mut ad);
    let mut amag = vec![0.0f64; m];
    {
        let (ptr, ri, va) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
        for j in 0..n {
            for k in ptr[j]..ptr[j + 1] {
                amag[ri[k] as usize] += (va[k] * dh[j]).abs();
            }
        }
    }
    for i in 0..m {
        // unscaled row residual relative to the row's magnitude
        if ad[i].abs() > 1e-9 * amag[i] {
            return false;
        }
    }
    let mut cd = 0.0;
    let mut cmag = 0.0;
    for j in 0..n {
        let v = dh[j] * lp.cscale[j]; // unscaled component
        let fl = hl[j] && !lp.art[j];
        let fu = hu[j] && !lp.art[j];
        if (v > 1e-9 && fu) || (v < -1e-9 && fl) {
            return false;
        }
        cd += lp.c[j] * dh[j];
        cmag += (lp.c[j] * dh[j]).abs();
    }
    cmag > 0.0 && cd < -1e-6 * cmag
}
