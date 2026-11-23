//! HiGHS source: pdlp/hipdlp/defs.hpp, pdlp/hipdlp/pdhg.cc (lines 1903-1938)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiPDLP enums and parameter structs. The `initialise` methods are
//! defined in `pdhg.cc` but belong to these structs, so they live here.
//! `PdlpIterate`, `pdlp_iterate_ops` and `DetailedTimings` are declared
//! but never used by the solver and are not ported.
//!
//! Defaults: members with a C++ default member initialiser take it; members
//! without one (`StepSizeConfig::primal_step`, `PrimalDualParams::eta`, ...)
//! are indeterminate in C++ and are always assigned before they are read
//! (`PrimalDualParams::initialise`, `initializeStepSizes`); `Default` gives
//! them `0` / the first variant.

use crate::prelude::*;

/// HiGHS: defs.hpp `enum class Device`
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Device {
    #[default]
    CPU,
    GPU,
}

/// HiGHS: defs.hpp `enum class ScalingMethod`
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ScalingMethod {
    #[default]
    NONE,
    RUIZ,
    POCK_CHAMBOLLE,
    L2_NORM,
    COMBINED,
}

/// HiGHS: defs.hpp `enum class RestartStrategy`
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RestartStrategy {
    #[default]
    NO_RESTART,
    FIXED_RESTART,
    ADAPTIVE_RESTART,
    HALPERN_RESTART,
}

/// HiGHS: defs.hpp `enum class StepSizeStrategy`
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StepSizeStrategy {
    #[default]
    FIXED,
    ADAPTIVE,
    MALITSKY_POCK,
    PID,
}

/// HiGHS: defs.hpp `enum class PostSolveRetcode`
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PostSolveRetcode {
    #[default]
    OK = 0,
    INVALID_SOLUTION = 1,
    SCALING_ERROR = 2,
    DIMENSION_MISMATCH = 3,
    NUMERICAL_ERROR = 4,
    CONSTRAINT_VIOLATION = 5,
}

/// HiGHS: defs.hpp `struct StepSizeConfig`
#[derive(Clone, Debug, Default)]
pub struct StepSizeConfig {
    pub primal_step: f64,
    pub dual_step: f64,
    pub beta: f64,
    pub power_method_lambda: f64,
    /// nStepSizeIter
    pub step_size_iter: HighsInt,
}

/// HiGHS: defs.hpp `struct MalitskyPockParams`
#[derive(Clone, Debug)]
pub struct MalitskyPockParams {
    /// Between 0 and 1
    pub step_size_interpolation: f64,
    pub step_size_downscaling_factor: f64,
    pub linesearch_contraction_factor: f64,
}

impl Default for MalitskyPockParams {
    fn default() -> Self {
        Self {
            step_size_interpolation: 0.5,
            step_size_downscaling_factor: 0.7,
            linesearch_contraction_factor: 0.99,
        }
    }
}

impl MalitskyPockParams {
    /// HiGHS: pdhg.cc::MalitskyPockParams::initialise (lines 1929-1933)
    pub fn initialise(&mut self) {
        self.step_size_interpolation = 0.5;
        self.step_size_downscaling_factor = 0.7;
        self.linesearch_contraction_factor = 0.99;
    }
}

/// HiGHS: defs.hpp `struct AdaptiveLinesearchParams`
#[derive(Clone, Debug)]
pub struct AdaptiveLinesearchParams {
    pub step_size_reduction_exponent: f64,
    pub step_size_growth_exponent: f64,
}

impl Default for AdaptiveLinesearchParams {
    fn default() -> Self {
        Self {
            step_size_reduction_exponent: 0.3,
            step_size_growth_exponent: 0.6,
        }
    }
}

impl AdaptiveLinesearchParams {
    /// HiGHS: pdhg.cc::AdaptiveLinesearchParams::initialise (lines 1935-1938)
    pub fn initialise(&mut self) {
        self.step_size_reduction_exponent = 0.3;
        self.step_size_growth_exponent = 0.6;
    }
}

/// HiGHS: defs.hpp `struct PrimalDualParams`
#[derive(Clone, Debug)]
pub struct PrimalDualParams {
    pub eta: f64,
    pub omega: f64,
    pub tolerance: f64,
    pub max_iterations: HighsInt,
    pub device_type: Device,
    pub time_limit: f64,
    // Restart parameters
    pub restart_strategy: RestartStrategy,
    pub fixed_restart_interval: HighsInt,
    pub use_halpern_restart: bool,
    /// 0: standard Halpern, 1: full reflection, in between: over relaxation
    pub halpern_gamma: f64,
    // Scaling parameters
    pub use_ruiz_scaling: bool,
    pub use_pc_scaling: bool,
    pub use_l2_scaling: bool,
    // Ruiz scaling parameters
    pub ruiz_iterations: HighsInt,
    pub ruiz_norm: f64,
    // Pock-Chambolle scaling parameters
    pub pc_alpha: f64,
    // Step sizes strategy
    pub step_size_strategy: StepSizeStrategy,
    pub k_p: f64,
    pub k_i: f64,
    pub k_d: f64,
    pub i_smooth: f64,
    pub malitsky_pock_params: MalitskyPockParams,
    pub adaptive_linesearch_params: AdaptiveLinesearchParams,
    pub log_options_: HighsLogOptions,
    /// rusolve (HiGHS 1.16 backport, default off): the #3469 flag, default false;
    /// set in PDLPSolver::setup from the option.
    pub h16_hipdlp_max: bool,
}

impl Default for PrimalDualParams {
    fn default() -> Self {
        Self {
            eta: 0.0,
            omega: 0.0,
            tolerance: 0.0,
            max_iterations: 0,
            device_type: Device::CPU,
            time_limit: 3600.0,
            restart_strategy: RestartStrategy::NO_RESTART,
            fixed_restart_interval: 0,
            use_halpern_restart: false,
            halpern_gamma: 1.0,
            use_ruiz_scaling: false,
            use_pc_scaling: false,
            use_l2_scaling: false,
            ruiz_iterations: 10,
            ruiz_norm: f64::INFINITY,
            pc_alpha: 1.0,
            step_size_strategy: StepSizeStrategy::FIXED,
            k_p: 0.99,
            k_i: 0.01,
            k_d: 0.0,
            i_smooth: 0.3,
            malitsky_pock_params: MalitskyPockParams::default(),
            adaptive_linesearch_params: AdaptiveLinesearchParams::default(),
            log_options_: HighsLogOptions::default(),
            h16_hipdlp_max: false,
        }
    }
}

impl PrimalDualParams {
    /// HiGHS: pdhg.cc::PrimalDualParams::initialise (lines 1903-1927)
    ///
    /// `i_smooth` is not reset (upstream omission, kept).
    /// `log_options_.clear()` (HighsIO.cpp 357-367 nulls every pointer)
    /// becomes `HighsLogOptions { log_stream: None, output_flag: false,
    /// log_to_console: false, log_dev_level: 0, user_callback: None,
    /// user_callback_active: false }`; `setup` overwrites it before use.
    pub fn initialise(&mut self) {
        self.eta = 0.0;
        self.omega = 0.0;
        self.tolerance = 0.0;
        self.max_iterations = 0;
        self.device_type = Device::CPU;
        self.time_limit = kHighsInf;
        self.restart_strategy = RestartStrategy::NO_RESTART;
        self.fixed_restart_interval = 0;
        self.use_halpern_restart = false;
        self.halpern_gamma = 1.0;
        self.use_ruiz_scaling = false;
        self.use_pc_scaling = false;
        self.use_l2_scaling = false;
        self.ruiz_iterations = 10;
        self.ruiz_norm = f64::INFINITY;
        self.pc_alpha = 1.0;
        self.k_p = 0.99;
        self.k_i = 0.01;
        self.k_d = 0.0;
        self.step_size_strategy = StepSizeStrategy::FIXED;
        self.malitsky_pock_params.initialise();
        self.adaptive_linesearch_params.initialise();
        // log_options_.clear(): HighsIO.cpp 357-367 nulls every pointer.
        self.log_options_ = HighsLogOptions {
            log_stream: None,
            output_flag: false,
            log_to_console: false,
            log_dev_level: 0,
            user_callback: None,
            user_callback_active: false,
        };
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn initialise_resets_to_header_values() {
        let mut p = PrimalDualParams {
            tolerance: 7.0,
            time_limit: 1.0,
            i_smooth: 0.9,
            ruiz_iterations: 3,
            ..Default::default()
        };
        p.malitsky_pock_params.step_size_interpolation = 0.1;
        p.initialise();
        assert_eq!(p.tolerance, 0.0);
        assert_eq!(p.time_limit, kHighsInf);
        assert_eq!(p.ruiz_iterations, 10);
        assert_eq!(p.i_smooth, 0.9);
        assert_eq!(p.malitsky_pock_params.step_size_interpolation, 0.5);
    }
}
