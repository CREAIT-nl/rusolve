//! HiGHS source: qpsolver/pricing.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//!
//!
//! PORT NOTE (ownership, lib.rs rule (a)): every C++ pricing class stores
//! `Runtime& runtime`, `Basis& basis` and `ReducedCosts& redcosts`. Here they
//! are parameters. `price` additionally takes the `Gradient` because
//! `redcosts.getReducedCosts()` may recompute through `gradient.getGradient()`.
//! The C++ `price(x, gradient)` vector argument is unused by every
//! implementation and is dropped; the caller (`Quass::solve`) still
//! evaluates `gradient.getGradient()` BEFORE `price`, exactly as the C++
//! argument evaluation does (it can trigger a recompute).

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: pricing.hpp::Pricing (lines 13-20)
pub trait Pricing {
    /// HiGHS: pricing.hpp::Pricing::price (line 15)
    fn price(
        &mut self,
        rt: &Runtime,
        basis: &mut Basis,
        redcosts: &mut ReducedCosts,
        gradient: &mut Gradient,
        x: &QpVector,
    ) -> HighsInt;

    /// HiGHS: pricing.hpp::Pricing::update_weights (lines 16-17)
    ///
    /// Called from inside `Basis::updatebasis` with that same basis.
    fn update_weights(
        &mut self,
        rt: &Runtime,
        basis: &mut Basis,
        aq: &QpVector,
        ep: &QpVector,
        p: HighsInt,
        q: HighsInt,
    );

    /// HiGHS: pricing.hpp::Pricing::recompute (line 18)
    fn recompute(&mut self, rt: &Runtime, basis: &mut Basis);
}
