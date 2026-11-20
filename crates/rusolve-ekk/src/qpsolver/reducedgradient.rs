//! HiGHS source: qpsolver/reducedgradient.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! C++ stores `Gradient& gradient` and `Basis& basis`; here they
//! are parameters (lib.rs rule (a)).

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: reducedgradient.hpp::ReducedGradient (lines 15-93)
#[derive(Clone, Debug, Default)]
pub struct ReducedGradient {
    rg: QpVector,
    uptodate: bool,
}

impl ReducedGradient {
    /// HiGHS: reducedgradient.hpp::ReducedGradient::ReducedGradient (lines 22-23)
    pub fn new(rt: &Runtime) -> Self {
        Self {
            rg: QpVector::new(rt.instance.num_var),
            uptodate: false,
        }
    }

    /// HiGHS: reducedgradient.hpp::ReducedGradient::get (lines 25-30)
    pub fn get(
        &mut self,
        rt: &Runtime,
        basis: &mut Basis,
        gradient: &mut Gradient,
    ) -> &mut QpVector {
        if !self.uptodate {
            self.recompute(rt, basis, gradient);
        }
        &mut self.rg
    }

    /// HiGHS: reducedgradient.hpp::ReducedGradient::recompute (lines 32-36)
    pub fn recompute(&mut self, rt: &Runtime, basis: &mut Basis, gradient: &mut Gradient) {
        self.rg.dim = basis.getinactive().len() as HighsInt;
        let grad = gradient.get_gradient(rt);
        basis.Ztprod(grad, &mut self.rg, false, -1);
        self.uptodate = true;
    }

    /// HiGHS: reducedgradient.hpp::ReducedGradient::reduce (lines 38-65)
    pub fn reduce(&mut self, buffer_d: &QpVector, maxabsd: HighsInt) {
        if !self.uptodate {
            return;
        }
        // C++ lines 42-51: commented-out code, dropped.

        for i in 0..buffer_d.num_nz {
            let idx = buffer_d.index[i as usize];
            if idx == maxabsd {
                continue;
            }
            self.rg.value[idx as usize] -= self.rg.value[maxabsd as usize]
                * buffer_d.value[idx as usize]
                / buffer_d.value[maxabsd as usize];
        }

        self.rg.resparsify();

        self.uptodate = true;
    }

    /// HiGHS: reducedgradient.hpp::ReducedGradient::expand (lines 67-78)
    pub fn expand(&mut self, rt: &Runtime, gradient: &mut Gradient, yp: &QpVector) {
        if !self.uptodate {
            return;
        }

        let grad = gradient.get_gradient(rt);
        let newval = yp.dot(grad);
        self.rg.value.push(newval);
        self.rg.index.push(0);
        // C++17 sequencing: RHS evaluated first (old dim), then num_nz and dim incremented.
        self.rg.index[self.rg.num_nz as usize] = self.rg.dim;
        self.rg.num_nz += 1;
        self.rg.dim += 1;

        self.uptodate = true;
    }

    /// HiGHS: reducedgradient.hpp::ReducedGradient::update (lines 80-92)
    pub fn update(&mut self, alpha: f64, minor: bool) {
        if !self.uptodate {
            return;
        }
        if minor {
            for i in 0..self.rg.num_nz {
                self.rg.value[self.rg.index[i as usize] as usize] *= 1.0 - alpha;
            }
            self.uptodate = true;
        } else {
            self.uptodate = false;
        }
    }
}
