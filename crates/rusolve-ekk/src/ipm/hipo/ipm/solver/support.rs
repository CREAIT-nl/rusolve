//! rusolve `lp_ipm_support_reduction`. Not in HiGHS.
//!
//! After HiPO converges to the analytic centre of the optimal face, every variable of the
//! optimal-face support is strictly interior, and crossover must push each one to a bound.
//! On datt256_lp crossover from the converged point needs 143784 primal pushes. Here the point is
//! first re-centred to a larger mu, the cost of every lower-bounded variable gets a small positive
//! random perturbation (a weighted L1 objective on the optimal face, which prefers vertices:
//! Megiddo 1991, Andersen and Ye 1996), a few predictor-corrector iterations run, and the true
//! cost is restored. If anything goes wrong the converged point is put back.

use super::Solver;
use crate::ipm::hipo::ipm::status::*;

/// Deterministic pseudo-random in [1, 2) from the column index (splitmix64).
fn rnd(j: usize) -> f64 {
    let mut z = (j as u64).wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    1.0 + (z >> 11) as f64 / (1u64 << 53) as f64
}

impl Solver {
    /// Number of lower-bounded variables whose primal slack exceeds the dual one (the interior
    /// support that crossover will have to push).
    fn supportSize(&self) -> usize {
        let it = self.it_.as_deref().unwrap();
        (0..self.model_.n())
            .filter(|&j| {
                self.model_.hasLb(j) && it.xl[j as usize] > it.zl[j as usize] && {
                    // not at the upper bound either
                    !self.model_.hasUb(j) || it.xu[j as usize] > it.zu[j as usize]
                }
            })
            .count()
    }

    pub fn reduceSupport(&mut self) {
        if !self.Hoptions_.lp_ipm_support_reduction || self.info_.status != kStatusPDFeas {
            return;
        }
        let delta_rel: f64 = 1e-6;
        let mu_target: f64 = 1e-6;
        let sr_tol: f64 = 1e-10;
        let max_extra: i32 = 100;

        let n = self.model_.n() as usize;
        let c0 = self.model_.c_.clone();
        let backup = {
            let it = self.it_.as_deref().unwrap();
            (
                it.x.clone(),
                it.xl.clone(),
                it.xu.clone(),
                it.y.clone(),
                it.zl.clone(),
                it.zu.clone(),
            )
        };
        let support0 = self.supportSize();

        // cost scale: mean absolute nonzero cost
        let (mut s, mut k) = (0.0, 0usize);
        for &c in &c0 {
            if c != 0.0 {
                s += c.abs();
                k += 1;
            }
        }
        let cscale = if k > 0 { s / k as f64 } else { 1.0 };
        let delta = delta_rel * cscale;

        // perturb and re-centre
        {
            let it = self.it_.as_deref_mut().unwrap();
            let model = &mut self.model_;
            for j in 0..n {
                let ji = j as i32;
                if model.hasLb(ji) && it.xl[j] > it.zl[j] {
                    model.c_[j] += delta * rnd(j);
                }
                if model.hasLb(ji) {
                    if it.xl[j] >= it.zl[j] {
                        it.zl[j] = it.zl[j].max(mu_target / it.xl[j].max(1e-300));
                    } else {
                        it.xl[j] = it.xl[j].max(mu_target / it.zl[j].max(1e-300));
                    }
                }
                if model.hasUb(ji) {
                    if it.xu[j] >= it.zu[j] {
                        it.zu[j] = it.zu[j].max(mu_target / it.xu[j].max(1e-300));
                    } else {
                        it.xu[j] = it.xu[j].max(mu_target / it.zu[j].max(1e-300));
                    }
                }
            }
            it.best_mu = 0.0;
            it.best_pinf_ = f64::INFINITY;
            it.best_dinf_ = f64::INFINITY;
            it.bad_iter_ = 0;
            it.residual1234(model);
            it.computeMu(model);
            it.indicators(model);
        }

        let iter0 = self.iter_;
        let maxit0 = self.options_.max_iter;
        let opt0 = self.options_.optimality_tol;
        let hist0 = std::mem::take(&mut self.stall_hist_);
        self.options_.max_iter = iter0 + max_extra;
        self.options_.optimality_tol = sr_tol;
        self.info_.status = kStatusNotRun;

        while self.iter_ < self.options_.max_iter {
            if self.prepareIter() {
                break;
            }
            if self.predictor() {
                break;
            }
            if self.correctors() {
                break;
            }
            self.makeStep();
        }
        let extra = self.iter_ - iter0;
        let ended_ok =
            self.info_.status == kStatusPDFeas || self.it_.as_deref().unwrap().pinf < 1e-6;
        let support_end = self.supportSize();

        // restore the true cost and options
        self.model_.c_ = c0;
        self.options_.max_iter = maxit0;
        self.options_.optimality_tol = opt0;
        self.stall_hist_ = hist0;
        self.info_.status = kStatusPDFeas;

        let finite = {
            let it = self.it_.as_deref().unwrap();
            !it.isNan() && !it.isInf()
        };
        let mut ended_ok = ended_ok;
        if ended_ok && finite {
            // the restored point must still be near-complementary for the true cost, else the
            // perturbation was too large for this instance (woodlands09: gap 0.1)
            let it = self.it_.as_deref_mut().unwrap();
            it.residual1234(&self.model_);
            it.indicators(&self.model_);
            if it.pdgap > 1e-5 {
                ended_ok = false;
            }
        }
        if ended_ok && finite {
            let it = self.it_.as_deref_mut().unwrap();
            let model = &self.model_;
            it.best_mu = 0.0;
            it.residual1234(model);
            it.residual1234(model);
            it.computeMu(model);
            it.indicators(model);
            self.support_reduced_ = true;
            let support1 = self.supportSize();
            self.logger_.print(&format!(
                "Support reduction: {} extra iterations, interior support {} -> {}, pinf {:.2e} dinf {:.2e} gap {:.2e}\n",
                extra, support0, support1, self.it_.as_deref().unwrap().pinf, self.it_.as_deref().unwrap().dinf, self.it_.as_deref().unwrap().pdgap
            ));
        } else {
            let it = self.it_.as_deref_mut().unwrap();
            let model = &self.model_;
            it.x = backup.0;
            it.xl = backup.1;
            it.xu = backup.2;
            it.y = backup.3;
            it.zl = backup.4;
            it.zu = backup.5;
            it.best_mu = 0.0;
            it.residual1234(model);
            it.computeMu(model);
            it.indicators(model);
            self.logger_.print(&format!(
                "Support reduction abandoned after {} iterations (support {} -> {}), converged point restored\n",
                extra, support0, support_end
            ));
        }
    }
}
