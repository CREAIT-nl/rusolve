//! HiGHS source: ipm/hipo/ipm/Iterate.h, ipm/hipo/ipm/Iterate.cpp (lines 1-297)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The IPM iterate, its residuals and the Newton direction. The C++
//! `const Model& model` and `const Regularisation& regul` members are
//! dropped: methods that read them take `model: &Model` /
//! `regul: &Regularisation`. `Rp` / `Rd` (pointers into `total_reg`, or
//! nullptr) are `Option<usize>` offsets into `total_reg`. `pobj`, `dobj`,
//! `pinf`, `dinf`, `pdgap`, `mu` are uninitialised in C++ until
//! `indicators` / `computeMu`; they start at zero.

pub mod part2;

use crate::ipm::hipo::Int;
use crate::ipm::hipo::auxiliary::vector_operations::dotProd;
use crate::ipm::hipo::auxiliary::vector_operations::infNorm;
use crate::ipm::hipo::auxiliary::vector_operations::isInfVector;
use crate::ipm::hipo::auxiliary::vector_operations::isNanVector;
use crate::ipm::hipo::auxiliary::vector_operations::vectorAdd;
use crate::ipm::hipo::ipm::ipm_data::IpmData;
use crate::ipm::hipo::ipm::model::Model;
use crate::ipm::hipo::ipm::parameters::Regularisation;
use crate::ipm::hipo::ipm::parameters::kLargeProduct;
use crate::ipm::hipo::ipm::parameters::kSmallProduct;
use crate::lp_data::hconst::kHighsInf;
use crate::util::highs_utils::cpp_max;
use crate::util::highs_utils::cpp_min;

/// HiGHS: Iterate.h `struct NewtonDir`
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NewtonDir {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub xl: Vec<f64>,
    pub xu: Vec<f64>,
    pub zl: Vec<f64>,
    pub zu: Vec<f64>,
}

impl NewtonDir {
    /// HiGHS: Iterate.cpp::NewtonDir::NewtonDir (lines 10-11)
    pub fn new(m: Int, n: Int) -> Self {
        NewtonDir {
            x: vec![0.0; n as usize],
            y: vec![0.0; m as usize],
            xl: vec![0.0; n as usize],
            xu: vec![0.0; n as usize],
            zl: vec![0.0; n as usize],
            zu: vec![0.0; n as usize],
        }
    }

    /// HiGHS: Iterate.cpp::NewtonDir::clear (lines 13-20)
    pub fn clear(&mut self) {
        for v in self.x.iter_mut() {
            *v = 0.0;
        }
        for v in self.y.iter_mut() {
            *v = 0.0;
        }
        for v in self.xl.iter_mut() {
            *v = 0.0;
        }
        for v in self.xu.iter_mut() {
            *v = 0.0;
        }
        for v in self.zl.iter_mut() {
            *v = 0.0;
        }
        for v in self.zu.iter_mut() {
            *v = 0.0;
        }
    }

    /// HiGHS: Iterate.cpp::NewtonDir::add (lines 22-29)
    pub fn add(&mut self, d: &NewtonDir) {
        vectorAdd(&mut self.x, &d.x, 1.0, 1.0);
        vectorAdd(&mut self.y, &d.y, 1.0, 1.0);
        vectorAdd(&mut self.xl, &d.xl, 1.0, 1.0);
        vectorAdd(&mut self.xu, &d.xu, 1.0, 1.0);
        vectorAdd(&mut self.zl, &d.zl, 1.0, 1.0);
        vectorAdd(&mut self.zu, &d.zu, 1.0, 1.0);
    }
}

/// HiGHS: Iterate.h `struct Residuals`
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Residuals {
    pub r1: Vec<f64>,
    pub r2: Vec<f64>,
    pub r3: Vec<f64>,
    pub r4: Vec<f64>,
    pub r5: Vec<f64>,
    pub r6: Vec<f64>,
}

/// HiGHS: Iterate.h `struct Iterate`
#[derive(Clone, Debug)]
pub struct Iterate {
    pub data: IpmData,
    pub x: Vec<f64>,
    pub xl: Vec<f64>,
    pub xu: Vec<f64>,
    pub y: Vec<f64>,
    pub zl: Vec<f64>,
    pub zu: Vec<f64>,
    pub res: Residuals,
    pub delta: NewtonDir,
    pub pobj: f64,
    pub dobj: f64,
    pub pinf: f64,
    pub dinf: f64,
    pub pdgap: f64,
    pub mu: f64,
    pub scaling: Vec<f64>,
    pub best_mu: f64,
    pub total_reg: Vec<f64>,
    pub Rp: Option<usize>,
    pub Rd: Option<usize>,
    pub ires: Residuals,
    pub bad_iter_: Int,
    pub largest_dx_x_: f64,
    pub largest_dy_y_: f64,
    pub best_pinf_: f64,
    pub best_dinf_: f64,
}

impl Iterate {
    /// HiGHS: Iterate.cpp::Iterate::Iterate (lines 31-38)
    /// `delta(m, n)`, then `clearIter`, `clearRes`, `clearIres`, `best_mu = 0`;
    /// `best_pinf_ = best_dinf_ = kHighsInf`, `Rp = Rd = None`.
    pub fn new(model: &Model) -> Self {
        let mut it = Iterate {
            data: IpmData::default(),
            x: Vec::new(),
            xl: Vec::new(),
            xu: Vec::new(),
            y: Vec::new(),
            zl: Vec::new(),
            zu: Vec::new(),
            res: Residuals::default(),
            delta: NewtonDir::new(model.m(), model.n()),
            pobj: 0.0,
            dobj: 0.0,
            pinf: 0.0,
            dinf: 0.0,
            pdgap: 0.0,
            mu: 0.0,
            scaling: Vec::new(),
            best_mu: 0.0,
            total_reg: Vec::new(),
            Rp: None,
            Rd: None,
            ires: Residuals::default(),
            bad_iter_: 0,
            largest_dx_x_: 0.0,
            largest_dy_y_: 0.0,
            best_pinf_: kHighsInf,
            best_dinf_: kHighsInf,
        };
        it.clearIter(model);
        it.clearRes(model);
        it.clearIres(model);
        it.best_mu = 0.0;
        it
    }

    /// HiGHS: Iterate.cpp::Iterate::isNan (lines 39-42)
    pub fn isNan(&self) -> bool {
        isNanVector(&self.x)
            || isNanVector(&self.xl)
            || isNanVector(&self.xu)
            || isNanVector(&self.y)
            || isNanVector(&self.zl)
            || isNanVector(&self.zu)
    }

    /// HiGHS: Iterate.cpp::Iterate::isInf (lines 43-46)
    pub fn isInf(&self) -> bool {
        isInfVector(&self.x)
            || isInfVector(&self.xl)
            || isInfVector(&self.xu)
            || isInfVector(&self.y)
            || isInfVector(&self.zl)
            || isInfVector(&self.zu)
    }

    /// HiGHS: Iterate.cpp::Iterate::isResNan (lines 47-50)
    pub fn isResNan(&self) -> bool {
        isNanVector(&self.res.r1)
            || isNanVector(&self.res.r2)
            || isNanVector(&self.res.r3)
            || isNanVector(&self.res.r4)
            || isNanVector(&self.res.r5)
            || isNanVector(&self.res.r6)
    }

    /// HiGHS: Iterate.cpp::Iterate::isResInf (lines 51-54)
    pub fn isResInf(&self) -> bool {
        isInfVector(&self.res.r1)
            || isInfVector(&self.res.r2)
            || isInfVector(&self.res.r3)
            || isInfVector(&self.res.r4)
            || isInfVector(&self.res.r5)
            || isInfVector(&self.res.r6)
    }

    /// HiGHS: Iterate.cpp::Iterate::isDirNan (lines 55-58)
    pub fn isDirNan(&self, d: &NewtonDir) -> bool {
        isNanVector(&d.x)
            || isNanVector(&d.xl)
            || isNanVector(&d.xu)
            || isNanVector(&d.y)
            || isNanVector(&d.zl)
            || isNanVector(&d.zu)
    }

    /// HiGHS: Iterate.cpp::Iterate::isDirInf (lines 59-62)
    pub fn isDirInf(&self, d: &NewtonDir) -> bool {
        isInfVector(&d.x)
            || isInfVector(&d.xl)
            || isInfVector(&d.xu)
            || isInfVector(&d.y)
            || isInfVector(&d.zl)
            || isInfVector(&d.zu)
    }

    /// HiGHS: Iterate.cpp::Iterate::computeMu (lines 64-83)
    pub fn computeMu(&mut self, model: &Model) {
        self.mu = 0.0;
        let mut number_finite_bounds: Int = 0;
        for i in 0..model.n() {
            if model.hasLb(i) {
                self.mu += self.xl[i as usize] * self.zl[i as usize];
                number_finite_bounds += 1;
            }
            if model.hasUb(i) {
                self.mu += self.xu[i as usize] * self.zu[i as usize];
                number_finite_bounds += 1;
            }
        }
        self.mu /= number_finite_bounds as f64;

        if self.best_mu > 0.0 {
            self.best_mu = cpp_min(self.best_mu, self.mu);
        } else {
            self.best_mu = self.mu;
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::computeScaling (lines 84-104)
    pub fn computeScaling(&mut self, model: &Model) {
        self.scaling.clear();
        self.scaling.resize(model.n() as usize, 0.0);

        for i in 0..model.n() {
            if model.hasLb(i) {
                self.scaling[i as usize] += self.zl[i as usize] / self.xl[i as usize];
            }
            if model.hasUb(i) {
                self.scaling[i as usize] += self.zu[i as usize] / self.xu[i as usize];
            }

            // slow down the growth of theta
            if self.scaling[i as usize] < 1e-12 {
                self.scaling[i as usize] = (1e-12 * self.scaling[i as usize]).sqrt();
            }
        }

        // compute extremes of theta
        self.data.back().min_theta = kHighsInf;
        self.data.back().max_theta = 0.0;
        for d in self.scaling.iter() {
            if *d != 0.0 {
                self.data.back().min_theta = cpp_min(self.data.back().min_theta, 1.0 / *d);
                self.data.back().max_theta = cpp_max(self.data.back().max_theta, 1.0 / *d);
            }
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::products (lines 105-132)
    pub fn products(&mut self, model: &Model) {
        self.data.back().min_prod = kHighsInf;
        self.data.back().max_prod = 0.0;
        self.data.back().num_small_prod = 0;
        self.data.back().num_large_prod = 0;

        for i in 0..model.n() {
            if model.hasLb(i) {
                let prod = self.xl[i as usize] * self.zl[i as usize] / self.mu;
                self.data.back().min_prod = cpp_min(self.data.back().min_prod, prod);
                self.data.back().max_prod = cpp_max(self.data.back().max_prod, prod);
                if prod < kSmallProduct {
                    self.data.back().num_small_prod += 1;
                }
                if prod > kLargeProduct {
                    self.data.back().num_large_prod += 1;
                }
            }
            if model.hasUb(i) {
                let prod = self.xu[i as usize] * self.zu[i as usize] / self.mu;
                self.data.back().min_prod = cpp_min(self.data.back().min_prod, prod);
                self.data.back().max_prod = cpp_max(self.data.back().max_prod, prod);
                if prod < kSmallProduct {
                    self.data.back().num_small_prod += 1;
                }
                if prod > kLargeProduct {
                    self.data.back().num_large_prod += 1;
                }
            }
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::indicators (lines 134-141)
    pub fn indicators(&mut self, model: &Model) {
        self.primalObj(model);
        self.dualObj(model);
        self.primalInfeasUnscaled(model);
        self.dualInfeasUnscaled(model);
        self.pdGap(model);
        self.products(model);
    }

    /// HiGHS: Iterate.cpp::Iterate::primalObj (lines 143-146)
    pub fn primalObj(&mut self, model: &Model) {
        self.pobj = model.offset() + dotProd(&self.x, model.c());
        if model.qp() {
            self.pobj += model.sense() * model.Q().objective_value(&self.x);
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::dualObj (lines 147-154)
    pub fn dualObj(&mut self, model: &Model) {
        self.dobj = model.offset() + dotProd(&self.y, model.b());
        for i in 0..model.n() {
            if model.hasLb(i) {
                self.dobj += model.lb(i) * self.zl[i as usize];
            }
            if model.hasUb(i) {
                self.dobj -= model.ub(i) * self.zu[i as usize];
            }
        }
        if model.qp() {
            self.dobj -= model.sense() * model.Q().objective_value(&self.x);
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::pdGap (lines 155-158)
    /// HiGHS 1.16 (#3427): ipm/hipo/ipm/Iterate.cpp::Iterate::pdGap (lines 158-168)
    /// (C++ uses the member `model`; passed here as a parameter, as in `indicators`.)
    pub fn pdGap(&mut self, model: &Model) {
        // relative primal-dual gap
        // HiGHS 1.16 #3427 (h16_ipm_gap_offset): ipm/hipo/ipm/Iterate.cpp lines 166-167 (1.16.0)
        self.pdgap = if model.h16_ipm_gap_offset_ {
            (self.pobj - self.dobj).abs()
                / (1.0 + 0.5 * (self.pobj + self.dobj - 2.0 * model.offset()).abs())
        } else {
            (self.pobj - self.dobj).abs() / (1.0 + 0.5 * (self.pobj + self.dobj).abs())
        };
    }

    /// HiGHS: Iterate.cpp::Iterate::primalInfeas (lines 160-166)
    pub fn primalInfeas(&mut self, model: &Model) {
        // relative infinity norm of scaled primal residuals
        self.pinf = infNorm(&self.res.r1);
        self.pinf = cpp_max(self.pinf, infNorm(&self.res.r2));
        self.pinf = cpp_max(self.pinf, infNorm(&self.res.r3));
        self.pinf /= 1.0 + model.normScaledRhs();
    }

    /// HiGHS: Iterate.cpp::Iterate::dualInfeas (lines 167-170)
    pub fn dualInfeas(&mut self, model: &Model) {
        // relative infinity norm of scaled dual residual
        self.dinf = infNorm(&self.res.r4) / (1.0 + model.normScaledObj());
    }

    /// HiGHS: Iterate.cpp::Iterate::primalInfeasUnscaled (lines 171-189)
    pub fn primalInfeasUnscaled(&mut self, model: &Model) {
        // relative infinity norm of unscaled primal residuals
        self.pinf = 0.0;
        for i in 0..model.m() {
            let mut val = self.res.r1[i as usize].abs();
            if model.scaled() {
                val /= model.rowScale(i);
            }
            self.pinf = cpp_max(self.pinf, val);
        }
        for i in 0..model.n() {
            let mut val = self.res.r2[i as usize].abs();
            if model.scaled() {
                val *= model.colScale(i);
            }
            self.pinf = cpp_max(self.pinf, val);

            val = self.res.r3[i as usize].abs();
            if model.scaled() {
                val *= model.colScale(i);
            }
            self.pinf = cpp_max(self.pinf, val);
        }
        self.pinf /= 1.0 + model.normUnscaledRhs();
    }

    /// HiGHS: Iterate.cpp::Iterate::dualInfeasUnscaled (lines 190-199)
    pub fn dualInfeasUnscaled(&mut self, model: &Model) {
        // relative infinity norm of unscaled dual residual
        self.dinf = 0.0;
        for i in 0..model.n() {
            let mut val = self.res.r4[i as usize].abs();
            if model.scaled() {
                val /= model.colScale(i);
            }
            self.dinf = cpp_max(self.dinf, val);
        }
        self.dinf /= 1.0 + model.normUnscaledObj();
    }

    /// HiGHS: Iterate.cpp::Iterate::residual1234 (lines 201-230)
    pub fn residual1234(&mut self, model: &Model) {
        // res1
        self.res.r1 = model.b().to_vec();
        model
            .A()
            .alpha_product_plus_y(-1.0, &self.x, &mut self.res.r1, false);

        // res2
        for i in 0..model.n() {
            if model.hasLb(i) {
                self.res.r2[i as usize] = model.lb(i) - self.x[i as usize] + self.xl[i as usize];
            } else {
                self.res.r2[i as usize] = 0.0;
            }
        }

        // res3
        for i in 0..model.n() {
            if model.hasUb(i) {
                self.res.r3[i as usize] = model.ub(i) - self.x[i as usize] - self.xu[i as usize];
            } else {
                self.res.r3[i as usize] = 0.0;
            }
        }

        // res4
        self.res.r4 = model.c().to_vec();
        model
            .A()
            .alpha_product_plus_y(-1.0, &self.y, &mut self.res.r4, true);
        for i in 0..model.n() {
            if model.hasLb(i) {
                self.res.r4[i as usize] -= self.zl[i as usize];
            }
            if model.hasUb(i) {
                self.res.r4[i as usize] += self.zu[i as usize];
            }
        }
        if model.qp() {
            model
                .Q()
                .alpha_product_plus_y(model.sense(), &self.x, &mut self.res.r4);
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::residual56 (lines 231-245)
    pub fn residual56(&mut self, model: &Model, sigma: f64) {
        for i in 0..model.n() {
            // res5
            if model.hasLb(i) {
                self.res.r5[i as usize] =
                    sigma * self.mu - self.xl[i as usize] * self.zl[i as usize];
            } else {
                self.res.r5[i as usize] = 0.0;
            }

            // res6
            if model.hasUb(i) {
                self.res.r6[i as usize] =
                    sigma * self.mu - self.xu[i as usize] * self.zu[i as usize];
            } else {
                self.res.r6[i as usize] = 0.0;
            }
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::residual7 (lines 247-254)
    pub fn residual7(&self, model: &Model, r: &Residuals) -> Vec<f64> {
        let mut res7 = r.r4.clone();
        for i in 0..model.n() {
            if model.hasLb(i) {
                res7[i as usize] -= (r.r5[i as usize] + self.zl[i as usize] * r.r2[i as usize])
                    / self.xl[i as usize];
            }
            if model.hasUb(i) {
                res7[i as usize] += (r.r6[i as usize] - self.zu[i as usize] * r.r3[i as usize])
                    / self.xu[i as usize];
            }
        }
        res7
    }

    /// HiGHS: Iterate.cpp::Iterate::residual8 (lines 255-271)
    pub fn residual8(
        &self,
        model: &Model,
        regul: &Regularisation,
        r: &Residuals,
        res7: &[f64],
    ) -> Vec<f64> {
        let mut res8 = r.r1.clone();
        let mut temp = res7.to_vec();

        // temp = (Theta^-1+Rp+Q)^-1 * res7
        for i in 0..model.n() {
            let mut denom = self.scaling[i as usize] + regul.primal;
            if model.qp() {
                denom += model.sense() * model.Q().diag(i);
            }
            temp[i as usize] /= denom;
        }

        // res8 += A * temp
        model.A().alpha_product_plus_y(1.0, &temp, &mut res8, false);

        res8
    }

    /// HiGHS: Iterate.cpp::Iterate::clearIter (lines 273-280)
    pub fn clearIter(&mut self, model: &Model) {
        self.x.clear();
        self.x.resize(model.n() as usize, 0.0);
        self.xl.clear();
        self.xl.resize(model.n() as usize, 0.0);
        self.xu.clear();
        self.xu.resize(model.n() as usize, 0.0);
        self.y.clear();
        self.y.resize(model.m() as usize, 0.0);
        self.zl.clear();
        self.zl.resize(model.n() as usize, 0.0);
        self.zu.clear();
        self.zu.resize(model.n() as usize, 0.0);
    }

    /// HiGHS: Iterate.cpp::Iterate::clearRes (lines 281-288)
    pub fn clearRes(&mut self, model: &Model) {
        self.res.r1.clear();
        self.res.r1.resize(model.m() as usize, 0.0);
        self.res.r2.clear();
        self.res.r2.resize(model.n() as usize, 0.0);
        self.res.r3.clear();
        self.res.r3.resize(model.n() as usize, 0.0);
        self.res.r4.clear();
        self.res.r4.resize(model.n() as usize, 0.0);
        self.res.r5.clear();
        self.res.r5.resize(model.n() as usize, 0.0);
        self.res.r6.clear();
        self.res.r6.resize(model.n() as usize, 0.0);
    }

    /// HiGHS: Iterate.cpp::Iterate::clearDir (lines 289-289)
    pub fn clearDir(&mut self) {
        self.delta.clear();
    }

    /// HiGHS: Iterate.cpp::Iterate::clearIres (lines 290-297)
    pub fn clearIres(&mut self, model: &Model) {
        self.ires.r1.clear();
        self.ires.r1.resize(model.m() as usize, 0.0);
        self.ires.r2.clear();
        self.ires.r2.resize(model.n() as usize, 0.0);
        self.ires.r3.clear();
        self.ires.r3.resize(model.n() as usize, 0.0);
        self.ires.r4.clear();
        self.ires.r4.resize(model.n() as usize, 0.0);
        self.ires.r5.clear();
        self.ires.r5.resize(model.n() as usize, 0.0);
        self.ires.r6.clear();
        self.ires.r6.resize(model.n() as usize, 0.0);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn newton_dir_add_and_clear() {
        let mut d = NewtonDir::new(1, 2);
        assert_eq!(d.x, vec![0.0, 0.0]);
        assert_eq!(d.y, vec![0.0]);
        let mut e = NewtonDir::new(1, 2);
        e.x = vec![1.0, 2.0];
        d.add(&e);
        assert_eq!(d.x, vec![1.0, 2.0]);
        d.clear();
        assert_eq!(d.x, vec![0.0, 0.0]);
    }

    #[test]
    fn h16_ipm_gap_offset_hipo_pdgap() {
        // Objective offset 1e6, true objectives near 1, gap 1e-3.
        let mut model = Model::default();
        model.offset_ = 1e6;
        let mut it = Iterate::new(&model);
        it.pobj = 1.0 + 1e6;
        it.dobj = 1.0 - 1e-3 + 1e6;
        let (pobj, dobj) = (it.pobj, it.dobj);

        // Flag off: 1.15.1 formula, denominator ~ 1 + 1e6.
        model.h16_ipm_gap_offset_ = false;
        it.pdGap(&model);
        let off = (pobj - dobj).abs() / (1.0 + 0.5 * (pobj + dobj).abs());
        assert_eq!(it.pdgap, off);
        assert!(it.pdgap < 1e-8);

        // Flag on (HiGHS 1.16 Iterate.cpp lines 166-167): offset-free
        // denominator ~ 2, so the same gap gives a much larger pdgap.
        model.h16_ipm_gap_offset_ = true;
        it.pdGap(&model);
        let on = (pobj - dobj).abs() / (1.0 + 0.5 * (pobj + dobj - 2.0 * model.offset()).abs());
        assert_eq!(it.pdgap, on);
        assert!(it.pdgap > 1e-4);
    }
}
