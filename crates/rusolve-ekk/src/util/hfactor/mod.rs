//! HiGHS source: util/HFactor.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Struct definitions. Method bodies are split by HFactor*.cpp line ranges:
//! `setup.rs` (setup*, build, ftranCall/btranCall, update, ...), `build.rs` (buildSimple,
//! rank-deficiency handling, buildFinish), `kernel.rs` (buildKernel + header-inline link/col/row
//! helpers), `solve.rs` (solveMatrixT/solveHyper + ftran/btran L/U/FT/PF/MPF/APF), `update.rs`
//! (updateFT/PF/MPF/APF, get/setInvert), `extend.rs` (HFactorExtend.cpp), `refactor.rs`
//! (HFactorRefactor.cpp), `utils.rs` (HFactorUtils.cpp).
//!
//! Ownership (crate docs, pattern (c)):
//! * C++ `a_start`/`a_index`/`a_value` are `const` pointers into the caller's matrix. Here they
//!   are owned copies (`Vec`s) refreshed by `setup_matrix*`/`setup_general*`: translate
//!   `a_start = a_start_;` as `self.a_start = a_start_.to_vec();`. C++ calls `setupMatrix`
//!   whenever the matrix storage may have moved (HSimplexNla::setPointers), so the copies stay in
//!   sync with the pointer semantics on the simplex path.
//! * C++ `basic_index` (a non-owning `HighsInt*` into `HEkk::basis_.basicIndex_`) is NOT stored:
//!   every method that reads/writes it takes `basic_index: &mut [HighsInt]` as its first
//!   parameter (build, buildSimple, buildKernel, buildHandleRankDeficiency, buildMarkSingC,
//!   buildFinish, update, updateAPF, rebuild, addRows). `getBaseIndex` is dropped.
//! * C++ `build_timer_` points at a local `HighsTimer` in `build()`; here `build_timer_` is an owned
//!   `HighsTimer` that `build` resets (`self.build_timer_ = HighsTimer::new(); self.build_timer_.start();`)
//!   and `buildKernel` reads for the time limit.
//! * C++ `log_data` (unique_ptr backing the pointer fields of `log_options`) is folded into
//!   `log_options` (plain values): `log_data->output_flag = x` -> `self.log_options.output_flag = x`.
//! * The PAMI-only `HVector::next` chain is a raw link; the chained case of C++
//!   `update` (`if (aq->next) updateCFT(...)`) is the `unsafe` [`HFactor::update_multi`],
//!   reached only from `HEkkDual::majorUpdateFactor`; `update` itself is the unchained case.

use crate::prelude::*;

pub mod build;
pub mod debug;
pub mod extend;
pub mod kernel;
pub mod refactor;
pub mod setup;
pub mod solve;
pub mod update;
pub mod utils;

pub use debug::*;

/// HiGHS: HFactor.h kBuildKernelReturnTimeout
pub const kBuildKernelReturnTimeout: HighsInt = -1;

/// HiGHS: HFactor.h::InvertibleRepresentation
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InvertibleRepresentation {
    pub l_pivot_index: Vec<HighsInt>,
    pub l_pivot_lookup: Vec<HighsInt>,
    pub l_start: Vec<HighsInt>,
    pub l_index: Vec<HighsInt>,
    pub l_value: Vec<f64>,
    pub lr_start: Vec<HighsInt>,
    pub lr_index: Vec<HighsInt>,
    pub lr_value: Vec<f64>,
    pub u_pivot_lookup: Vec<HighsInt>,
    pub u_pivot_index: Vec<HighsInt>,
    pub u_pivot_value: Vec<f64>,
    pub u_start: Vec<HighsInt>,
    pub u_last_p: Vec<HighsInt>,
    pub u_index: Vec<HighsInt>,
    pub u_value: Vec<f64>,
    pub ur_start: Vec<HighsInt>,
    pub ur_lastp: Vec<HighsInt>,
    pub ur_space: Vec<HighsInt>,
    pub ur_index: Vec<HighsInt>,
    pub ur_value: Vec<f64>,
    pub pf_start: Vec<HighsInt>,
    pub pf_index: Vec<HighsInt>,
    pub pf_value: Vec<f64>,
    pub pf_pivot_index: Vec<HighsInt>,
    pub pf_pivot_value: Vec<f64>,
}

/// HiGHS: HFactor.h::HFactor
/// All members are `pub` (C++ private members included) so the split impl files can reach them.
/// `Default` = the C++ constructor's initialiser list.
#[derive(Clone, Debug)]
pub struct HFactor {
    pub build_real_tick: f64,
    pub build_synthetic_tick: f64,
    pub rank_deficiency: HighsInt,
    pub row_with_no_pivot: Vec<HighsInt>,
    pub col_with_no_pivot: Vec<HighsInt>,
    pub var_with_no_pivot: Vec<HighsInt>,
    pub refactor_info_: RefactorInfo,
    pub basis_matrix_num_el: HighsInt,
    pub invert_num_el: HighsInt,
    pub kernel_dim: HighsInt,
    pub kernel_num_el: HighsInt,
    pub num_row: HighsInt,
    pub num_col: HighsInt,
    pub num_basic: HighsInt,
    pub inv_num_row: f64,
    pub a_matrix_valid: bool,
    pub pivot_threshold: f64,
    pub pivot_tolerance: f64,
    pub highs_debug_level: HighsInt,
    pub time_limit_: f64,
    pub log_options: HighsLogOptions,
    pub use_original_h_factor_logic: bool,
    pub debug_report_: bool,
    pub basis_matrix_limit_size: HighsInt,
    pub update_method: HighsInt,
    pub nwork: HighsInt,
    pub iwork: Vec<HighsInt>,
    pub dwork: Vec<f64>,
    pub b_var: Vec<HighsInt>,
    pub b_start: Vec<HighsInt>,
    pub b_index: Vec<HighsInt>,
    pub b_value: Vec<f64>,
    pub permute: Vec<HighsInt>,
    pub mc_var: Vec<HighsInt>,
    pub mc_start: Vec<HighsInt>,
    pub mc_count_a: Vec<HighsInt>,
    pub mc_count_n: Vec<HighsInt>,
    pub mc_space: Vec<HighsInt>,
    pub mc_index: Vec<HighsInt>,
    pub mc_value: Vec<f64>,
    pub mc_min_pivot: Vec<f64>,
    pub mr_start: Vec<HighsInt>,
    pub mr_count: Vec<HighsInt>,
    pub mr_space: Vec<HighsInt>,
    pub mr_count_before: Vec<HighsInt>,
    pub mr_index: Vec<HighsInt>,
    pub mwz_column_index: Vec<HighsInt>,
    pub mwz_column_mark: Vec<i8>,
    pub mwz_column_array: Vec<f64>,
    pub col_link_first: Vec<HighsInt>,
    pub col_link_next: Vec<HighsInt>,
    pub col_link_last: Vec<HighsInt>,
    pub row_link_first: Vec<HighsInt>,
    pub row_link_next: Vec<HighsInt>,
    pub row_link_last: Vec<HighsInt>,
    pub l_pivot_lookup: Vec<HighsInt>,
    pub l_pivot_index: Vec<HighsInt>,
    pub l_start: Vec<HighsInt>,
    pub l_index: Vec<HighsInt>,
    pub l_value: Vec<f64>,
    pub lr_start: Vec<HighsInt>,
    pub lr_index: Vec<HighsInt>,
    pub lr_value: Vec<f64>,
    pub u_pivot_lookup: Vec<HighsInt>,
    pub u_pivot_index: Vec<HighsInt>,
    pub u_pivot_value: Vec<f64>,
    pub u_merit_x: HighsInt,
    pub u_total_x: HighsInt,
    pub u_start: Vec<HighsInt>,
    pub u_last_p: Vec<HighsInt>,
    pub u_index: Vec<HighsInt>,
    pub u_value: Vec<f64>,
    pub ur_start: Vec<HighsInt>,
    pub ur_lastp: Vec<HighsInt>,
    pub ur_space: Vec<HighsInt>,
    pub ur_index: Vec<HighsInt>,
    pub ur_value: Vec<f64>,
    pub pf_pivot_value: Vec<f64>,
    pub pf_pivot_index: Vec<HighsInt>,
    pub pf_start: Vec<HighsInt>,
    pub pf_index: Vec<HighsInt>,
    pub pf_value: Vec<f64>,
    pub rhs_: HVector,
    pub a_start: Vec<HighsInt>,
    pub a_index: Vec<HighsInt>,
    pub a_value: Vec<f64>,
    pub build_timer_: HighsTimer,
    /// rusolve (HiGHS 1.16 backport #3352 #3354): pivot flag, default false;
    /// set in HSimplexNla::setup from the option. Copied by Clone.
    pub h16_hfactor_pivot: bool,
}

impl Default for HFactor {
    fn default() -> Self {
        Self {
            build_real_tick: Default::default(),
            build_synthetic_tick: Default::default(),
            rank_deficiency: Default::default(),
            row_with_no_pivot: Default::default(),
            col_with_no_pivot: Default::default(),
            var_with_no_pivot: Default::default(),
            refactor_info_: Default::default(),
            basis_matrix_num_el: Default::default(),
            invert_num_el: Default::default(),
            kernel_dim: Default::default(),
            kernel_num_el: Default::default(),
            num_row: Default::default(),
            num_col: Default::default(),
            num_basic: Default::default(),
            inv_num_row: Default::default(),
            a_matrix_valid: Default::default(),
            pivot_threshold: Default::default(),
            pivot_tolerance: Default::default(),
            highs_debug_level: Default::default(),
            time_limit_: Default::default(),
            log_options: Default::default(),
            use_original_h_factor_logic: Default::default(),
            debug_report_: Default::default(),
            basis_matrix_limit_size: Default::default(),
            update_method: Default::default(),
            nwork: Default::default(),
            iwork: Default::default(),
            dwork: Default::default(),
            b_var: Default::default(),
            b_start: Default::default(),
            b_index: Default::default(),
            b_value: Default::default(),
            permute: Default::default(),
            mc_var: Default::default(),
            mc_start: Default::default(),
            mc_count_a: Default::default(),
            mc_count_n: Default::default(),
            mc_space: Default::default(),
            mc_index: Default::default(),
            mc_value: Default::default(),
            mc_min_pivot: Default::default(),
            mr_start: Default::default(),
            mr_count: Default::default(),
            mr_space: Default::default(),
            mr_count_before: Default::default(),
            mr_index: Default::default(),
            mwz_column_index: Default::default(),
            mwz_column_mark: Default::default(),
            mwz_column_array: Default::default(),
            col_link_first: Default::default(),
            col_link_next: Default::default(),
            col_link_last: Default::default(),
            row_link_first: Default::default(),
            row_link_next: Default::default(),
            row_link_last: Default::default(),
            l_pivot_lookup: Default::default(),
            l_pivot_index: Default::default(),
            l_start: Default::default(),
            l_index: Default::default(),
            l_value: Default::default(),
            lr_start: Default::default(),
            lr_index: Default::default(),
            lr_value: Default::default(),
            u_pivot_lookup: Default::default(),
            u_pivot_index: Default::default(),
            u_pivot_value: Default::default(),
            u_merit_x: Default::default(),
            u_total_x: Default::default(),
            u_start: Default::default(),
            u_last_p: Default::default(),
            u_index: Default::default(),
            u_value: Default::default(),
            ur_start: Default::default(),
            ur_lastp: Default::default(),
            ur_space: Default::default(),
            ur_index: Default::default(),
            ur_value: Default::default(),
            pf_pivot_value: Default::default(),
            pf_pivot_index: Default::default(),
            pf_start: Default::default(),
            pf_index: Default::default(),
            pf_value: Default::default(),
            rhs_: Default::default(),
            a_start: Default::default(),
            a_index: Default::default(),
            a_value: Default::default(),
            build_timer_: Default::default(),
            h16_hfactor_pivot: false,
        }
    }
}

impl HFactor {
    /// HiGHS: HFactor.h::HFactor::getAstart (inline getter)
    pub fn get_astart(&self) -> &[HighsInt] {
        &self.a_start
    }

    /// HiGHS: HFactor.h::HFactor::getAindex (inline getter)
    pub fn get_aindex(&self) -> &[HighsInt] {
        &self.a_index
    }

    /// HiGHS: HFactor.h::HFactor::getAvalue (inline getter)
    pub fn get_avalue(&self) -> &[f64] {
        &self.a_value
    }

    /// HiGHS: HFactor.h::HFactor::setDebugReport (inline setter)
    pub fn set_debug_report(&mut self, debug_report: bool) {
        self.debug_report_ = debug_report;
    }
}
