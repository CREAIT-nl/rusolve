//! HiGHS source: qpsolver/settings.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete: data only. `std::numeric_limits<HighsInt>::infinity()`
//! is `0` (integers have no infinity), so the three "infinite" integer
//! frequencies/limits default to `0`, exactly as in HiGHS.

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: settings.hpp::RatiotestStrategy (line 15)
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RatiotestStrategy {
    TwoPass = 0,
    Textbook,
}

/// HiGHS: settings.hpp::PricingStrategy (line 17)
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PricingStrategy {
    SteepestEdge = 0,
    DantzigWolfe,
    Devex,
}

/// HiGHS: settings.hpp::Phase1Strategy (line 19)
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Phase1Strategy {
    HIGHS = 0,
    QUASS,
    BOUNDED,
}

/// HiGHS: settings.hpp::Settings (lines 21-83)
#[derive(Clone, Debug)]
pub struct Settings {
    pub ratiotest: RatiotestStrategy,
    pub ratiotest_t: f64,
    pub ratiotest_d: f64,
    /// HiGHS 1.16 (#3153): settings.hpp::Settings (line 25, primal_feasibility_tolerance)
    pub primal_feasibility_tolerance: f64,
    /// HiGHS 1.16 #3153 (h16_qp_bounded_start): flag carrier, set from HighsOptions in
    /// lp_data/highs_qp.rs (no HighsOptions reaches the qpsolver)
    pub h16_qp_bounded_start: bool,
    /// HiGHS 1.16 #3325 (h16_qp_convexity): flag carrier
    pub h16_qp_convexity: bool,
    /// HiGHS 1.16 #3354 (h16_hfactor_pivot): flag carrier for the basis factor's HFactor, set
    /// from HighsOptions in lp_data/highs_qp.rs (C++ `HFactor basisfactor;` runs the
    /// 1.16.0 kernel)
    pub h16_hfactor_pivot: bool,

    pub pricing: PricingStrategy,

    pub pnorm_zero_threshold: f64,
    pub improvement_zero_threshold: f64,
    pub d_zero_threshold: f64,
    pub lambda_zero_threshold: f64,
    pub pQp_zero_threshold: f64,

    pub hessian_regularization_on: bool,
    pub hessian_regularization_value: f64,

    pub allow_hot_start: bool,
    pub phase1strategy: Phase1Strategy,
    pub phase1movefreevarsbasic: bool,
    pub phase1boundfreevars: bool,

    pub reportingfequency: HighsInt,
    pub iteration_log: Eventhandler<Statistics>,
    pub qp_model_status_log: Eventhandler<QpModelStatus>,
    pub nullspace_limit_log: Eventhandler<HighsInt>,
    pub degeneracy_fail_log: Eventhandler<(HighsInt, f64)>,

    pub nullspace_limit: HighsInt,

    pub reinvertfrequency: HighsInt,
    pub gradientrecomputefrequency: HighsInt,
    pub reducedgradientrecomputefrequency: HighsInt,
    pub reducedhessianrecomputefrequency: HighsInt,

    pub iteration_limit: HighsInt,
    pub time_limit: f64,

    pub rowscaling: bool,
    pub varscaling: bool,

    pub perturbation: bool,
}

impl Default for Settings {
    /// HiGHS: settings.hpp::Settings (default member initialisers, lines 22-82)
    fn default() -> Self {
        Self {
            ratiotest: RatiotestStrategy::TwoPass,
            ratiotest_t: 1e-9,
            ratiotest_d: 1e-8,
            primal_feasibility_tolerance: 1e-7,
            h16_qp_bounded_start: false,
            h16_qp_convexity: false,
            h16_hfactor_pivot: false,
            pricing: PricingStrategy::Devex,
            pnorm_zero_threshold: 1e-11,
            improvement_zero_threshold: 1e-4,
            d_zero_threshold: 1e-12,
            lambda_zero_threshold: 1e-9,
            pQp_zero_threshold: 1e-7,
            hessian_regularization_on: false,
            hessian_regularization_value: 1e-7,
            allow_hot_start: false,
            phase1strategy: Phase1Strategy::HIGHS,
            phase1movefreevarsbasic: false,
            phase1boundfreevars: false,
            reportingfequency: 1,
            iteration_log: Eventhandler::default(),
            qp_model_status_log: Eventhandler::default(),
            nullspace_limit_log: Eventhandler::default(),
            degeneracy_fail_log: Eventhandler::default(),
            nullspace_limit: 4000,
            reinvertfrequency: 1000,
            gradientrecomputefrequency: 100,
            // std::numeric_limits<HighsInt>::infinity() == 0
            reducedgradientrecomputefrequency: 0,
            reducedhessianrecomputefrequency: 0,
            iteration_limit: 0,
            time_limit: f64::INFINITY,
            rowscaling: true,
            varscaling: true,
            perturbation: false,
        }
    }
}
