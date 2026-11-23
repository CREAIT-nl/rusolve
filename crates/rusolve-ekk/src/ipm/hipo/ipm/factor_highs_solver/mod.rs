//! HiGHS source: ipm/hipo/ipm/FactorHiGHSSolver.h, ipm/hipo/ipm/FactorHiGHSSolver.cpp (lines 1-200)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `FactorHiGHSSolver`: the `LinearSolver` built on FactorHiGHS. The cached
//! C++ references are the [`LsCtx`] argument. `chooseOrdering` borrows
//! `FH_` and the ordering name separately from `self`, so it is an
//! associated function taking them explicitly (its C++ arguments are
//! `kkt_.rowsAS` etc., which live in the context). `chooseOrdering`
//! (choose.rs) and `chooseNla`, `setNla`, `setParallel`
//! (parallel.rs) are the private helpers of `setup` and
//! `analyseAS` / `analyseNE`.

pub mod choose;
pub mod parallel;

use super::info::Info;
use super::kkt_matrix::KktAsPart;
use super::linear_solver::{LinearSolver, LsCtx};
use super::model::Model;
use super::options::Options;
use crate::ipm::hipo::auxiliary::Clock;
use crate::ipm::hipo::auxiliary::logger::{Logger, fix, textline};
use crate::ipm::hipo::factorhighs::factor_highs::FHsolver;
use crate::ipm::hipo::factorhighs::symbolic::Symbolic;
use crate::ipm::hipo::ipm::parameters::kLargeStorageGB;
use crate::ipm::hipo::ipm::status::*;
use crate::ipm::hipo::oracle_trace;
use crate::ipm::hipo::{Int, Int64};
use crate::lp_data::highs_options::kHipoAugmentedString;
use std::sync::atomic::Ordering;

/// HiGHS: FactorHiGHSSolver.h `class FactorHiGHSSolver`
pub struct FactorHiGHSSolver {
    pub valid_: bool,
    pub FH_: FHsolver,
    pub S_: Symbolic,
    pub ordering_AS_: String,
    pub ordering_NE_: String,
    /// rusolve `ipm_hipo_dense_cols`: scaling of the last NE factorisation, for the CG operator
    ne_scaling_: Vec<f64>,
    /// rusolve `ipm_hipo_dense_schur`: Woodbury data of the last NE factorisation
    dw_: Option<DenseWoodbury>,
}

/// rusolve `ipm_hipo_dense_schur`. With M the factored matrix without
/// the dense columns and V = A_D diag(sqrt(d_D)), the full matrix is M + V V^T and
/// (M + V V^T)^-1 = M^-1 - M^-1 V (I + V^T M^-1 V)^-1 V^T M^-1 (Sherman-Morrison-Woodbury; Andersen
/// 1996). `w` holds M^-1 V column by column, `chol` the dense Cholesky factor of the capacitance
/// I + V^T M^-1 V (lower, row major k x k).
struct DenseWoodbury {
    cols: Vec<usize>,
    sqrt_d: Vec<f64>,
    w: Vec<Vec<f64>>,
    chol: Vec<f64>,
}

impl FactorHiGHSSolver {
    /// rusolve `ipm_hipo_dense_cols` (Andersen 1996, "A modified Schur-complement method for
    /// handling dense columns in interior-point methods for linear programming", ACM TOMS 22;
    /// Marxen 1989 for the CG variant): the factor is of M_S = A_S D A_S^T + Rd without the dense
    /// columns, which is the preconditioner. `lhs` holds M_S^-1 rhs on entry; this runs
    /// preconditioned CG on the full M = A D A^T + Rd. Returns false when CG does not converge.
    fn dense_cg(&mut self, ctx: &mut LsCtx<'_>, rhs: &[f64], x: &mut Vec<f64>) -> bool {
        const kTol: f64 = 1e-12;
        const kMaxIter: usize = 300;
        let a = ctx.model.A();
        let n = a.num_col_ as usize;
        let m = rhs.len();
        let rp = ctx.regul.primal;
        let rd = ctx.regul.dual;
        let d: Vec<f64> = (0..n)
            .map(|j| {
                let s = if self.ne_scaling_.is_empty() {
                    1.0
                } else {
                    self.ne_scaling_[j]
                };
                1.0 / (s + rp)
            })
            .collect();
        let apply = |v: &[f64], out: &mut Vec<f64>| {
            out.clear();
            out.extend(v.iter().map(|&t| rd * t));
            for j in 0..n {
                let (b, e) = (a.start_[j] as usize, a.start_[j + 1] as usize);
                let mut t = 0.0;
                for k in b..e {
                    t += a.value_[k] * v[a.index_[k] as usize];
                }
                t *= d[j];
                for k in b..e {
                    out[a.index_[k] as usize] += a.value_[k] * t;
                }
            }
        };
        let dot = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(p, q)| p * q).sum::<f64>();
        let bnorm = dot(rhs, rhs).sqrt();
        if bnorm == 0.0 {
            return true;
        }
        let mut mx = Vec::with_capacity(m);
        apply(x, &mut mx);
        let mut r: Vec<f64> = rhs.iter().zip(&mx).map(|(b, t)| b - t).collect();
        if dot(&r, &r).sqrt() <= kTol * bnorm {
            return true;
        }
        let mut z = r.clone();
        if !self.precond(ctx, &mut z) {
            return false;
        }
        let mut p = z.clone();
        let mut rz = dot(&r, &z);
        let mut ap = Vec::with_capacity(m);
        for it in 0..kMaxIter {
            apply(&p, &mut ap);
            let pap = dot(&p, &ap);
            if !(pap > 0.0) || !pap.is_finite() {
                return false;
            }
            let alpha = rz / pap;
            for i in 0..m {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            if dot(&r, &r).sqrt() <= kTol * bnorm {
                return true;
            }
            z.clone_from(&r);
            if !self.precond(ctx, &mut z) {
                return false;
            }
            let rz_new = dot(&r, &z);
            let beta = rz_new / rz;
            rz = rz_new;
            for i in 0..m {
                p[i] = z[i] + beta * p[i];
            }
        }
        false
    }

    /// Builds the Woodbury data after a NE factorisation. Returns without it when the capacitance is not
    /// positive definite or too large; the plain CG then runs on M^-1 alone.
    fn build_dense_woodbury(&mut self, ctx: &mut LsCtx<'_>) {
        self.dw_ = None;
        if !crate::ipm::hipo::blas::DENSE_SCHUR.load(Ordering::Relaxed) {
            return;
        }
        let a = ctx.model.A();
        let thr = super::kkt_matrix::dense_col_threshold(a, ctx.model.qp());
        if thr == Int::MAX {
            return;
        }
        let cols: Vec<usize> = (0..a.num_col_ as usize)
            .filter(|&j| a.start_[j + 1] - a.start_[j] > thr)
            .collect();
        let k = cols.len();
        let m = a.num_row_ as usize;
        if k == 0 || k > 2000 || k * m > 400_000_000 {
            return;
        }
        let rp = ctx.regul.primal;
        let mut sqrt_d = Vec::with_capacity(k);
        let mut w: Vec<Vec<f64>> = Vec::with_capacity(k);
        let mut v: Vec<Vec<f64>> = Vec::with_capacity(k);
        for &j in &cols {
            let sc = if self.ne_scaling_.is_empty() {
                1.0
            } else {
                self.ne_scaling_[j]
            };
            let sd = (1.0 / (sc + rp)).sqrt();
            sqrt_d.push(sd);
            let mut col = vec![0.0; m];
            for t in a.start_[j] as usize..a.start_[j + 1] as usize {
                col[a.index_[t] as usize] = sd * a.value_[t];
            }
            w.push(col.clone());
            v.push(col);
        }
        // the k solves against one factor are independent: run them on the HiGHS scheduler
        struct WPtr(*mut Vec<f64>);
        // SAFETY: task i touches only w[i]; the index ranges of different tasks are disjoint and
        // the factor is only read (`solve_shared`).
        unsafe impl Sync for WPtr {}
        let wp = WPtr(w.as_mut_ptr());
        let wp = &wp;
        let failed = std::sync::atomic::AtomicBool::new(false);
        let fh = &self.FH_;
        let sym = &self.S_;
        let run = |c0: Int, c1: Int| {
            for i in c0 as usize..c1 as usize {
                // SAFETY: see WPtr
                let x = unsafe { &mut *wp.0.add(i) };
                if fh.solve_shared(sym, x) != 0 {
                    failed.store(true, Ordering::Relaxed);
                }
            }
        };
        crate::parallel::highs_parallel::for_each(0, k as Int, &run, 1);
        if failed.load(Ordering::Relaxed) {
            return;
        }
        let mut c = vec![0.0; k * k];
        for i in 0..k {
            for l in 0..=i {
                let t: f64 = v[i].iter().zip(&w[l]).map(|(p, q)| p * q).sum();
                c[i * k + l] = t + if i == l { 1.0 } else { 0.0 };
            }
        }
        // dense Cholesky C = L L^T (lower, row major)
        for i in 0..k {
            for l in 0..=i {
                let mut t = c[i * k + l];
                for q in 0..l {
                    t -= c[i * k + q] * c[l * k + q];
                }
                if i == l {
                    if !(t > 0.0) || !t.is_finite() {
                        return;
                    }
                    c[i * k + i] = t.sqrt();
                } else {
                    c[i * k + l] = t / c[l * k + l];
                }
            }
        }
        self.dw_ = Some(DenseWoodbury {
            cols,
            sqrt_d,
            w,
            chol: c,
        });
    }

    /// z <- (M + V V^T)^-1 z when the Woodbury data exists, else M^-1 z. False on a solve error.
    fn precond(&mut self, ctx: &LsCtx<'_>, z: &mut Vec<f64>) -> bool {
        if self.FH_.solve(&self.S_, z) != 0 {
            return false;
        }
        let Some(dw) = self.dw_.as_ref() else {
            return true;
        };
        let a = ctx.model.A();
        let k = dw.cols.len();
        let mut t = vec![0.0; k];
        for (i, &j) in dw.cols.iter().enumerate() {
            let mut s = 0.0;
            for q in a.start_[j] as usize..a.start_[j + 1] as usize {
                s += a.value_[q] * z[a.index_[q] as usize];
            }
            t[i] = dw.sqrt_d[i] * s;
        }
        for i in 0..k {
            let mut s = t[i];
            for q in 0..i {
                s -= dw.chol[i * k + q] * t[q];
            }
            t[i] = s / dw.chol[i * k + i];
        }
        for i in (0..k).rev() {
            let mut s = t[i];
            for q in i + 1..k {
                s -= dw.chol[q * k + i] * t[q];
            }
            t[i] = s / dw.chol[i * k + i];
        }
        for i in 0..k {
            for (zi, wi) in z.iter_mut().zip(&dw.w[i]) {
                *zi -= t[i] * wi;
            }
        }
        true
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::FactorHiGHSSolver (lines 14-19)
    /// `FH_(&logger, options.block_size)`; the other initialisers are the
    /// references that became [`LsCtx`].
    pub fn new(options: &Options, logger: &Logger) -> Self {
        FactorHiGHSSolver {
            valid_: false,
            FH_: FHsolver::new(Some(logger), options.block_size),
            S_: Symbolic::new(),
            ordering_AS_: "none".to_string(),
            ordering_NE_: "none".to_string(),
            ne_scaling_: Vec::new(),
            dw_: None,
        }
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::analyseAS (lines 37-58); the body is
    /// [`FactorHiGHSSolver::analyseAS_part`].
    pub fn analyseAS(&mut self, ctx: &mut LsCtx<'_>, S: &mut Symbolic) -> Int {
        let (kkt_AS, _) = ctx.kkt.split_AS_NE();
        self.analyseAS_part(
            &kkt_AS,
            ctx.model,
            ctx.options,
            ctx.logger,
            &mut ctx.info.analyse_AS_time,
            S,
        )
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::analyseAS (lines 37-58), on the AS
    /// members of `kkt_` only, so that `chooseNla`'s AS task runs it beside the NE task
    /// (lines 275-280). C++ `info_.analyse_AS_time` is `analyse_AS_time`.
    pub fn analyseAS_part(
        &mut self,
        kkt_AS: &KktAsPart<'_>,
        model: &Model,
        options: &Options,
        logger: &Logger,
        analyse_AS_time: &mut f64,
        S: &mut Symbolic,
    ) -> Int {
        // Perform analyse phase of augmented system and return symbolic factorisation
        // in object S and the status.

        if kkt_AS.rowsAS.is_empty() || kkt_AS.ptrAS.is_empty() {
            return kStatusErrorAnalyse;
        }

        let m = model.A().num_row_;
        let n = model.A().num_col_;

        // create vector of signs of pivots
        let mut pivot_signs = vec![-1 as Int; (n + m) as usize];
        for i in 0..m {
            pivot_signs[(n + i) as usize] = 1;
        }

        logger.printInfo("Performing AS analyse phase\n");

        let clock = Clock::new();
        let status = FactorHiGHSSolver::chooseOrdering(
            &mut self.FH_,
            options,
            logger,
            kkt_AS.NE_nz_limit.load(Ordering::Relaxed),
            kkt_AS.rowsAS,
            kkt_AS.ptrAS,
            &pivot_signs,
            S,
            &mut self.ordering_AS_,
            "AS",
        );
        *analyse_AS_time += clock.stop();

        status
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::analyseNE (lines 60-78)
    pub fn analyseNE(&mut self, ctx: &mut LsCtx<'_>, S: &mut Symbolic) -> Int {
        // Perform analyse phase of normal equations and return symbolic factorisation
        // in object S and the status. Structure of the matrix must be already
        // computed.

        if ctx.kkt.rowsNE.is_empty() || ctx.kkt.ptrNE.is_empty() {
            return kStatusErrorAnalyse;
        }

        // create vector of signs of pivots
        let pivot_signs = vec![1 as Int; ctx.model.A().num_row_ as usize];

        ctx.logger.printInfo("Performing NE analyse phase\n");

        let clock = Clock::new();
        let status = FactorHiGHSSolver::chooseOrdering(
            &mut self.FH_,
            ctx.options,
            ctx.logger,
            ctx.kkt.NE_nz_limit.load(Ordering::Relaxed),
            &ctx.kkt.rowsNE,
            &ctx.kkt.ptrNE,
            &pivot_signs,
            S,
            &mut self.ordering_NE_,
            "NE",
        );
        ctx.info.analyse_NE_time += clock.stop();

        status
    }
}

impl LinearSolver for FactorHiGHSSolver {
    fn valid(&self) -> bool {
        self.valid_
    }

    fn set_valid(&mut self, valid: bool) {
        self.valid_ = valid;
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::factorAS (lines 84-101)
    fn factorAS(&mut self, ctx: &mut LsCtx<'_>, scaling: &[f64]) -> Int {
        // only execute factorisation if it has not been done yet
        debug_assert!(!self.valid_);

        ctx.kkt.buildASvalues(ctx.model, ctx.info, scaling);

        // set static regularisation, since it may have changed
        self.FH_.setRegularisation(ctx.regul.primal, ctx.regul.dual);

        let clock = Clock::new();
        self.FH_.deadline_ = ctx.options.deadline_;
        if self
            .FH_
            .factorise(&self.S_, &ctx.kkt.rowsAS, &ctx.kkt.ptrAS, &ctx.kkt.valAS)
            != 0
        {
            if self.FH_.timed_out_ {
                return kStatusTimeLimit;
            }
            return kStatusErrorFactorise;
        }
        ctx.info.factor_time += clock.stop();
        ctx.info.factor_number += 1;

        self.valid_ = true;
        kStatusOk
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::factorNE (lines 103-120)
    fn factorNE(&mut self, ctx: &mut LsCtx<'_>, scaling: &[f64]) -> Int {
        // only execute factorisation if it has not been done yet
        debug_assert!(!self.valid_);

        ctx.kkt
            .buildNEvalues(ctx.model, ctx.regul, ctx.info, scaling);
        self.ne_scaling_ = scaling.to_vec();

        // set static regularisation, since it may have changed
        self.FH_.setRegularisation(ctx.regul.primal, ctx.regul.dual);

        let clock = Clock::new();
        self.FH_.deadline_ = ctx.options.deadline_;
        if self
            .FH_
            .factorise(&self.S_, &ctx.kkt.rowsNE, &ctx.kkt.ptrNE, &ctx.kkt.valNE)
            != 0
        {
            if self.FH_.timed_out_ {
                return kStatusTimeLimit;
            }
            return kStatusErrorFactorise;
        }
        ctx.info.factor_time += clock.stop();
        ctx.info.factor_number += 1;

        self.valid_ = true;
        self.build_dense_woodbury(ctx);
        kStatusOk
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::solveAS (lines 126-153)
    fn solveAS(
        &mut self,
        ctx: &mut LsCtx<'_>,
        rhs_x: &[f64],
        rhs_y: &[f64],
        lhs_x: &mut Vec<f64>,
        lhs_y: &mut Vec<f64>,
    ) -> Int {
        // only execute the solve if factorisation is valid
        debug_assert!(self.valid_);

        let n = rhs_x.len();

        // create single rhs
        let mut rhs: Vec<f64> = Vec::with_capacity(rhs_x.len() + rhs_y.len());
        rhs.extend_from_slice(rhs_x);
        rhs.extend_from_slice(rhs_y);

        let clock = Clock::new();
        if self.FH_.solve(&self.S_, &mut rhs) != 0 {
            return kStatusErrorSolve;
        }

        ctx.info.solve_time += clock.stop();
        ctx.info.solve_number += 1;

        ctx.data.back().num_solves += 1;

        // split lhs
        *lhs_x = rhs[..n].to_vec();
        *lhs_y = rhs[n..].to_vec();

        kStatusOk
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::solveNE (lines 155-172)
    fn solveNE(&mut self, ctx: &mut LsCtx<'_>, rhs: &[f64], lhs: &mut Vec<f64>) -> Int {
        // only execute the solve if factorisation is valid
        debug_assert!(self.valid_);

        // initialise lhs with rhs
        *lhs = rhs.to_vec();

        let clock = Clock::new();
        if self.dw_.is_some() {
            if !self.precond(ctx, lhs) {
                return kStatusErrorSolve;
            }
        } else if self.FH_.solve(&self.S_, lhs) != 0 {
            return kStatusErrorSolve;
        }
        let a = ctx.model.A();
        let dense_thr = super::kkt_matrix::dense_col_threshold(a, ctx.model.qp());
        if dense_thr != Int::MAX
            && (0..a.num_col_ as usize).any(|j| a.start_[j + 1] - a.start_[j] > dense_thr)
            && !self.dense_cg(ctx, rhs, lhs)
        {
            return kStatusErrorSolve;
        }

        ctx.info.solve_time += clock.stop();
        ctx.info.solve_number += 1;

        ctx.data.back().num_solves += 1;

        kStatusOk
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::clear (lines 28-31)
    fn clear(&mut self) {
        self.valid_ = false;
        self.FH_.newIter();
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::setup (lines 178-200)
    /// Oracle trace hook: `oracle_trace::hipo_fact` right after `setParallel()`
    /// (`tools/highs-oracle/FORMAT.md`, HiPO records).
    fn setup(&mut self, ctx: &mut LsCtx<'_>) -> Int {
        let clock = Clock::new();

        let status = self.setNla(ctx);
        if status != 0 {
            return status;
        }
        self.setParallel(ctx);

        // oracle_trace::hipo_fact right after setParallel()
        {
            let nla_str = if ctx.options.nla == kHipoAugmentedString {
                "AS"
            } else {
                "NE"
            };
            let ordering = if ctx.options.nla == kHipoAugmentedString {
                &self.ordering_AS_
            } else {
                &self.ordering_NE_
            };
            oracle_trace::hipo_fact(
                nla_str,
                ordering,
                self.S_.flops(),
                self.S_.spops(),
                self.S_.nz() as f64,
            );
        }

        if !ctx.options.timeless_log {
            let mut log_stream = String::new();
            log_stream.push_str(&textline("Analyse time:"));
            log_stream.push_str(&fix(clock.stop(), 0, 2));
            log_stream.push('\n');
            ctx.logger.print(&log_stream);
        }

        self.S_.print(ctx.logger, ctx.logger.debug(1));

        // rusolve: an IPM needs at least ~20 factorisations; at a
        // sustained kFlopRate one factorisation of this size alone overruns the time limit, so
        // hand the LP to IPX at once (no-progress status) instead of hanging in HiPO
        const kFlopRate: f64 = 5.0e10;
        const kMinIter: f64 = 20.0;
        if ctx.options.flops_bailout
            && ctx.options.time_limit > 0.0
            && ctx.options.time_limit.is_finite()
            && kMinIter * self.S_.flops() / kFlopRate > ctx.options.time_limit
        {
            ctx.logger
                .printw("HiPO factorisation too expensive for the time limit\n");
            return kStatusNoProgress;
        }

        // Warn about large memory consumption
        if self.S_.storage() > kLargeStorageGB * 1024.0 * 1024.0 * 1024.0 {
            ctx.logger.printw("Large amount of memory required\n");
        }

        ctx.logger.print("\n");

        kStatusOk
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::flops (lines 665-665)
    fn flops(&self) -> f64 {
        self.S_.flops()
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::spops (lines 666-666)
    fn spops(&self) -> f64 {
        self.S_.spops()
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::nz (lines 667-667)
    fn nz(&self) -> f64 {
        self.S_.nz() as f64
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::getReg (lines 668-670)
    fn getReg(&mut self, reg: &mut Vec<f64>) {
        self.FH_.getRegularisation(&self.S_, reg);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests, hand-derived from FactorHiGHSSolver.cpp.
    use super::*;
    use crate::ipm::hipo::ipm::ipm_data::IpmData;
    use crate::ipm::hipo::ipm::kkt_matrix::KktMatrix;
    use crate::ipm::hipo::ipm::model::Model;
    use crate::ipm::hipo::ipm::parameters::Regularisation;

    #[test]
    fn new_and_clear() {
        // Constructor (lines 14-26) and header defaults; clear (lines 28-31).
        let options = Options::default();
        let logger = Logger::default();
        let mut s = FactorHiGHSSolver::new(&options, &logger);
        assert!(!s.valid_);
        assert_eq!(s.ordering_AS_, "none");
        assert_eq!(s.ordering_NE_, "none");
        s.valid_ = true;
        s.clear();
        assert!(!s.valid());
    }

    #[test]
    fn analyse_without_structure_fails_before_timing() {
        // Lines 41 and 65: empty structure returns kStatusErrorAnalyse before
        // the clock, so the analyse times are untouched.
        let mut kkt = KktMatrix::new();
        let mut options = Options::default();
        let model = Model::default();
        let regul = Regularisation::default();
        let mut info = Info::default();
        let mut data = IpmData::default();
        let logger = Logger::default();
        let mut s = FactorHiGHSSolver::new(&options, &logger);
        let mut ctx = LsCtx {
            kkt: &mut kkt,
            options: &mut options,
            model: &model,
            regul: &regul,
            info: &mut info,
            data: &mut data,
            logger: &logger,
        };
        let mut S = Symbolic::new();
        assert_eq!(s.analyseAS(&mut ctx, &mut S), kStatusErrorAnalyse);
        assert_eq!(s.analyseNE(&mut ctx, &mut S), kStatusErrorAnalyse);
        ctx.kkt.rowsAS = vec![0];
        assert_eq!(s.analyseAS(&mut ctx, &mut S), kStatusErrorAnalyse);
        ctx.kkt.ptrNE = vec![0, 1];
        assert_eq!(s.analyseNE(&mut ctx, &mut S), kStatusErrorAnalyse);
        assert_eq!(info.analyse_AS_time, 0.0);
        assert_eq!(info.analyse_NE_time, 0.0);
        assert_eq!(s.ordering_AS_, "none");
        assert_eq!(s.ordering_NE_, "none");
    }
}
