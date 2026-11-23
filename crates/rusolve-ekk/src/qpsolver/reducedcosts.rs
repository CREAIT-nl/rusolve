//! HiGHS source: qpsolver/reducedcosts.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! C++ stores `Basis& basis` and `Gradient& gradient`; here they
//! are parameters (lib.rs rule (a)).

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: reducedcosts.hpp::ReducedCosts (lines 16-44)
#[derive(Clone, Debug, Default)]
pub struct ReducedCosts {
    reducedcosts: QpVector,
    uptodate: bool,
}

impl ReducedCosts {
    /// HiGHS: reducedcosts.hpp::ReducedCosts::ReducedCosts (lines 25-29)
    pub fn new(rt: &Runtime) -> Self {
        Self {
            reducedcosts: QpVector::new(rt.instance.num_var),
            uptodate: false,
        }
    }

    /// HiGHS: reducedcosts.hpp::ReducedCosts::recompute (lines 31-34)
    pub fn recompute(&mut self, rt: &Runtime, basis: &mut Basis, gradient: &mut Gradient) {
        let grad = gradient.get_gradient(rt);
        basis.ftran_into(grad, &mut self.reducedcosts, false, -1);
        self.uptodate = true;
    }

    /// HiGHS: reducedcosts.hpp::ReducedCosts::getReducedCosts (lines 36-41)
    pub fn get_reduced_costs(
        &mut self,
        rt: &Runtime,
        basis: &mut Basis,
        gradient: &mut Gradient,
    ) -> &mut QpVector {
        if !self.uptodate {
            self.recompute(rt, basis, gradient);
        }
        &mut self.reducedcosts
    }

    /// HiGHS: reducedcosts.hpp::ReducedCosts::update (line 43)
    pub fn update(&mut self) {
        self.uptodate = false;
    }
}
