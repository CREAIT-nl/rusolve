//! HiGHS source: util/HFactor.cpp (lines 160-558)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Setup, build driver, FTRAN/BTRAN entry points, update dispatch.

use super::*;
use crate::prelude::*;

impl HFactor {
    /// HiGHS: HFactor.cpp::HFactor::setup (lines 160-173, 14 lines)
    /// C++ default arguments: pivot_threshold = kDefaultPivotThreshold; pivot_tolerance = kDefaultPivotTolerance; highs_debug_level = kHighsDebugLevelMin; log_options = NULL
    /// Overload `setup(const HighsSparseMatrix&, std::vector<HighsInt>& basic_index, ...)`.
    /// Only `basic_index.size()` is used (as `num_basic`); the pointer is not stored.
    pub fn setup(
        &mut self,
        a_matrix: &HighsSparseMatrix,
        basic_index: &mut Vec<HighsInt>,
        pivot_threshold: f64,
        pivot_tolerance: f64,
        highs_debug_level: HighsInt,
        log_options: Option<&HighsLogOptions>,
    ) {
        let basic_index_size = basic_index.len() as HighsInt;
        // Nothing to do if basic index has no entries, and mustn't try to
        // pass the pointer to entry 0 of a vector of size 0.
        if basic_index_size <= 0 {
            return;
        }
        self.setup_general(
            a_matrix,
            basic_index_size,
            pivot_threshold,
            pivot_tolerance,
            highs_debug_level,
            log_options,
        );
    }

    /// HiGHS: HFactor.cpp::HFactor::setupGeneral (lines 175-186, 12 lines)
    /// C++ default arguments: pivot_threshold = kDefaultPivotThreshold; pivot_tolerance = kDefaultPivotTolerance; highs_debug_level = kHighsDebugLevelMin; log_options = NULL
    /// Overload `setupGeneral(const HighsSparseMatrix*, num_basic, HighsInt* basic_index, ...)`;
    /// the `basic_index` parameter is dropped (only stored as a pointer in C++).
    pub fn setup_general(
        &mut self,
        a_matrix: &HighsSparseMatrix,
        num_basic: HighsInt,
        pivot_threshold: f64,
        pivot_tolerance: f64,
        highs_debug_level: HighsInt,
        log_options: Option<&HighsLogOptions>,
    ) {
        self.setup_general_raw(
            a_matrix.num_col_,
            a_matrix.num_row_,
            num_basic,
            &a_matrix.start_,
            &a_matrix.index_,
            &a_matrix.value_,
            pivot_threshold,
            pivot_tolerance,
            highs_debug_level,
            // C++ null log_options: setupGeneral applies output_flag = false,
            // log_to_console = true, log_dev_level = 0, log_stream = nullptr.
            log_options.unwrap_or(&HighsLogOptions {
                log_stream: None,
                output_flag: false,
                log_to_console: true,
                log_dev_level: 0,
                user_callback: None,
                user_callback_active: false,
            }),
            true,
            kUpdateMethodFt,
        );
    }

    /// HiGHS: HFactor.cpp::HFactor::setup (lines 188-198, 11 lines)
    /// C++ default arguments: pivot_threshold_ = kDefaultPivotThreshold; pivot_tolerance_ = kDefaultPivotTolerance; highs_debug_level_ = kHighsDebugLevelMin; log_options_ = NULL; use_original_HFactor_logic_ = true; update_method_ = kUpdateMethodFt
    /// Overload `setup(num_col, num_row, a_start, a_index, a_value, HighsInt* basic_index, ...)`;
    /// the `basic_index_` parameter is dropped (only forwarded to be stored).
    pub fn setup_raw(
        &mut self,
        num_col_: HighsInt,
        num_row_: HighsInt,
        a_start_: &[HighsInt],
        a_index_: &[HighsInt],
        a_value_: &[f64],
        pivot_threshold_: f64,
        pivot_tolerance_: f64,
        highs_debug_level_: HighsInt,
        log_options_: Option<&HighsLogOptions>,
        use_original_h_factor_logic_: bool,
        update_method_: HighsInt,
    ) {
        self.setup_general_raw(
            num_col_,
            num_row_,
            num_row_,
            a_start_,
            a_index_,
            a_value_,
            pivot_threshold_,
            pivot_tolerance_,
            highs_debug_level_,
            // C++ null log_options: see setup_general.
            log_options_.unwrap_or(&HighsLogOptions {
                log_stream: None,
                output_flag: false,
                log_to_console: true,
                log_dev_level: 0,
                user_callback: None,
                user_callback_active: false,
            }),
            use_original_h_factor_logic_,
            update_method_,
        );
    }

    /// HiGHS: HFactor.cpp::HFactor::setupGeneral (lines 200-341, 142 lines)
    /// Overload `setupGeneral(num_col, num_row, num_basic, a_start, a_index, a_value, HighsInt* basic_index, ...)`;
    /// the `basic_index_` parameter is dropped: `basic_index = basic_index_;` has no Rust counterpart.
    /// `a_start = a_start_;` etc. become owned copies (`to_vec()`); `log_data` folds into `log_options`.
    pub fn setup_general_raw(
        &mut self,
        num_col_: HighsInt,
        num_row_: HighsInt,
        num_basic_: HighsInt,
        a_start_: &[HighsInt],
        a_index_: &[HighsInt],
        a_value_: &[f64],
        pivot_threshold_: f64,
        pivot_tolerance_: f64,
        highs_debug_level_: HighsInt,
        log_options_: &HighsLogOptions,
        use_original_h_factor_logic_: bool,
        update_method_: HighsInt,
    ) {
        // Copy Problem size and (copy of) coefficient matrix
        self.num_row = num_row_;
        self.num_col = num_col_;
        self.num_basic = num_basic_;
        self.inv_num_row = 1.0 / self.num_row as f64;
        self.a_matrix_valid = true;
        self.a_start = a_start_.to_vec();
        self.a_index = a_index_.to_vec();
        self.a_value = a_value_.to_vec();
        self.pivot_threshold = cpp_max(
            kMinPivotThreshold,
            cpp_min(pivot_threshold_, kMaxPivotThreshold),
        );
        self.pivot_tolerance = cpp_max(
            kMinPivotTolerance,
            cpp_min(pivot_tolerance_, kMaxPivotTolerance),
        );
        self.highs_debug_level = highs_debug_level_;
        self.time_limit_ = kHighsInf;
        // log_data folds into log_options. A null C++ log_options_ is passed by
        // the callers as output_flag = false, log_to_console = true,
        // log_dev_level = 0, no stream: exactly what the C++ `!log_options_`
        // branch sets.
        self.log_options.output_flag = log_options_.output_flag;
        self.log_options.log_to_console = log_options_.log_to_console;
        self.log_options.log_dev_level = log_options_.log_dev_level;
        self.log_options.log_stream = log_options_.log_stream;

        self.use_original_h_factor_logic = use_original_h_factor_logic_;
        self.update_method = update_method_;

        // Allocate for working buffer
        self.iwork.reserve((self.num_row * 2) as usize);
        self.dwork.clear();
        self.dwork.resize(self.num_row as usize, 0.0);

        // Find Basis matrix limit size
        self.basis_matrix_limit_size = 0;

        self.iwork.clear();
        self.iwork.resize((self.num_row + 1) as usize, 0);
        for i in 0..self.num_col {
            self.iwork[(self.a_start[(i + 1) as usize] - self.a_start[i as usize]) as usize] += 1;
        }
        let b_max_dim = cpp_max(self.num_row, self.num_basic);
        let mut i = self.num_row;
        let mut counted = 0;
        while i >= 0 && counted < b_max_dim {
            self.basis_matrix_limit_size += i * self.iwork[i as usize];
            counted += self.iwork[i as usize];
            i -= 1;
        }
        self.basis_matrix_limit_size += b_max_dim;

        // Allocate space for basis matrix, L, U factor and Update buffer
        self.b_var.resize(b_max_dim as usize, 0);
        self.b_start.resize((b_max_dim + 1) as usize, 0);
        self.b_index
            .resize(self.basis_matrix_limit_size as usize, 0);
        self.b_value
            .resize(self.basis_matrix_limit_size as usize, 0.0);

        // Allocate space for pivot records
        let permute_max_dim = cpp_max(self.num_row, self.num_basic);
        self.permute.resize(permute_max_dim as usize, 0);

        // Allocate space for Markowitz matrices
        let mc_dim = self.num_basic as usize;
        self.mc_var.resize(mc_dim, 0);
        self.mc_start.resize(mc_dim, 0);
        self.mc_count_a.resize(mc_dim, 0);
        self.mc_count_n.resize(mc_dim, 0);
        self.mc_space.resize(mc_dim, 0);
        self.mc_min_pivot.resize(mc_dim, 0.0);
        self.mc_index.resize(
            (self.basis_matrix_limit_size * kMCExtraEntriesMultiplier) as usize,
            0,
        );
        self.mc_value.resize(
            (self.basis_matrix_limit_size * kMCExtraEntriesMultiplier) as usize,
            0.0,
        );

        let num_row_u = self.num_row as usize;
        self.mr_start.resize(num_row_u, 0);
        self.mr_count.resize(num_row_u, 0);
        self.mr_space.resize(num_row_u, 0);
        self.mr_count_before.resize(num_row_u, 0);
        self.mr_index.resize(
            (self.basis_matrix_limit_size * kMRExtraEntriesMultiplier) as usize,
            0,
        );

        self.mwz_column_mark.clear();
        self.mwz_column_mark.resize(num_row_u, 0);
        self.mwz_column_index.resize(num_row_u, 0);
        self.mwz_column_array.clear();
        self.mwz_column_array.resize(num_row_u, 0.0);

        // Allocate space for count-link-list
        let col_link_max_count = self.num_row;
        let col_link_dim = self.num_basic;
        self.col_link_first.clear();
        self.col_link_first
            .resize((col_link_max_count + 1) as usize, -1);
        self.col_link_next.resize(col_link_dim as usize, 0);
        self.col_link_last.resize(col_link_dim as usize, 0);

        let row_link_max_count = self.num_basic;
        let row_link_dim = self.num_row;
        // Strictly not necessary, but nice if it's the right size
        self.row_link_first
            .resize((row_link_max_count + 1) as usize, 0);
        self.row_link_first.clear();
        self.row_link_first
            .resize((row_link_max_count + 1) as usize, -1);
        self.row_link_next.resize(row_link_dim as usize, 0);
        self.row_link_last.resize(row_link_dim as usize, 0);

        // Allocate space for L factor
        self.l_pivot_lookup.resize(num_row_u, 0);
        self.l_pivot_index.reserve(num_row_u);
        self.l_start.reserve(num_row_u + 1);
        self.l_index
            .reserve((self.basis_matrix_limit_size * kLFactorExtraEntriesMultiplier) as usize);
        self.l_value
            .reserve((self.basis_matrix_limit_size * kLFactorExtraEntriesMultiplier) as usize);

        self.lr_start.reserve(num_row_u + 1);
        self.lr_index
            .reserve((self.basis_matrix_limit_size * kLFactorExtraEntriesMultiplier) as usize);
        self.lr_value
            .reserve((self.basis_matrix_limit_size * kLFactorExtraEntriesMultiplier) as usize);

        // Allocate space for U factor
        self.u_pivot_lookup.resize(num_row_u, 0);
        self.u_pivot_index
            .reserve((self.num_row + kUFactorExtraVectors) as usize);
        self.u_pivot_value
            .reserve((self.num_row + kUFactorExtraVectors) as usize);

        self.u_start
            .reserve((self.num_row + kUFactorExtraVectors + 1) as usize);
        self.u_last_p
            .reserve((self.num_row + kUFactorExtraVectors) as usize);
        self.u_index
            .reserve((self.basis_matrix_limit_size * kUFactorExtraEntriesMultiplier) as usize);
        self.u_value
            .reserve((self.basis_matrix_limit_size * kUFactorExtraEntriesMultiplier) as usize);

        self.ur_start
            .reserve((self.num_row + kUFactorExtraVectors + 1) as usize);
        self.ur_lastp
            .reserve((self.num_row + kUFactorExtraVectors) as usize);
        self.ur_space
            .reserve((self.num_row + kUFactorExtraVectors) as usize);
        self.ur_index
            .reserve((self.basis_matrix_limit_size * kUFactorExtraEntriesMultiplier) as usize);
        self.ur_value
            .reserve((self.basis_matrix_limit_size * kUFactorExtraEntriesMultiplier) as usize);

        // Allocate spaces for Update buffer
        self.pf_pivot_value.reserve(kPFFPivotEntries as usize);
        self.pf_pivot_index.reserve(kPFFPivotEntries as usize);
        self.pf_start.reserve((kPFVectors + 1) as usize);
        self.pf_index
            .reserve((self.basis_matrix_limit_size * kPFEntriesMultiplier) as usize);
        self.pf_value
            .reserve((self.basis_matrix_limit_size * kPFEntriesMultiplier) as usize);

        // Set up the local HVector for use when RHS is
        // std::vector<double>.
        self.rhs_.setup(self.num_row);
        self.rhs_.count = -1;
    }

    /// HiGHS: HFactor.cpp::HFactor::setupMatrix (lines 343-349, 7 lines)
    /// Overload `setupMatrix(const HighsInt* a_start, const HighsInt* a_index, const double* a_value)`: copies.
    pub fn setup_matrix_raw(
        &mut self,
        a_start_: &[HighsInt],
        a_index_: &[HighsInt],
        a_value_: &[f64],
    ) {
        self.a_start = a_start_.to_vec();
        self.a_index = a_index_.to_vec();
        self.a_value = a_value_.to_vec();
        self.a_matrix_valid = true;
    }

    /// HiGHS: HFactor.cpp::HFactor::setupMatrix (lines 351-354, 4 lines)
    /// Overload `setupMatrix(const HighsSparseMatrix*)`.
    pub fn setup_matrix(&mut self, a_matrix: &HighsSparseMatrix) {
        self.setup_matrix_raw(&a_matrix.start_, &a_matrix.index_, &a_matrix.value_);
    }

    /// HiGHS: HFactor.cpp::HFactor::build (lines 356-461, 106 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    /// `HighsTimer build_timer; build_timer_ = &build_timer; build_timer.start();` ->
    /// `self.build_timer_ = HighsTimer::new(); self.build_timer_.start();`
    pub fn build(
        &mut self,
        basic_index: &mut [HighsInt],
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) -> HighsInt {
        // Set up a timer to prevent build running longer than time_limit_,
        // which is kHighsInf by default, and only set to a finite value in
        // HPresolve::removeDependentEquations
        self.build_timer_ = HighsTimer::new();
        self.build_timer_.start();

        let report_lu = false;
        // Ensure that the A matrix is valid for factorization
        debug_assert!(self.a_matrix_valid);
        let factor_timer = FactorTimer;
        // Possibly use the refactorization information!
        if self.refactor_info_.r#use {
            factor_timer.start(FactorReinvert, factor_timer_clock_pointer);
            self.rank_deficiency = self.rebuild(basic_index, factor_timer_clock_pointer);
            factor_timer.stop(FactorReinvert, factor_timer_clock_pointer);
            if self.rank_deficiency == 0 {
                return 0;
            }
        }
        // Refactoring from just the list of basic variables. Initialise the
        // refactorization information.
        self.refactor_info_.clear();
        // Start the timer
        factor_timer.start(FactorInvert, factor_timer_clock_pointer);
        self.build_synthetic_tick = 0.0;
        factor_timer.start(FactorInvertSimple, factor_timer_clock_pointer);
        // Build the L, U factor
        self.build_simple(basic_index);
        if report_lu {
            // printf("\nAfter units and singletons\n"); reportLu(kReportLuBoth, false);
        }
        factor_timer.stop(FactorInvertSimple, factor_timer_clock_pointer);
        factor_timer.start(FactorInvertKernel, factor_timer_clock_pointer);
        let build_kernel_return = self.build_kernel(basic_index);
        factor_timer.stop(FactorInvertKernel, factor_timer_clock_pointer);
        //
        // build_kernel_return of kBuildKernelReturnTimeout (<0) indicates
        // that time limit has been reached, otherwise it's the rank
        // deficiency of the basic variables. If num_basic < num_row, then
        // have to identify the logicals required to complete the basis by
        // continuing as if a full-dimension set of basic variables was rank
        // deficient.
        if build_kernel_return == kBuildKernelReturnTimeout {
            return kBuildKernelReturnTimeout;
        }
        self.rank_deficiency = build_kernel_return;
        let incomplete_basis = self.num_basic < self.num_row;
        if self.rank_deficiency != 0 || incomplete_basis {
            factor_timer.start(FactorInvertDeficient, factor_timer_clock_pointer);
            if self.num_basic == self.num_row {
                crate::highs_log_dev!(
                    &self.log_options,
                    HighsLogType::kWarning,
                    "Rank deficiency of {} identified in basis matrix\n",
                    self.rank_deficiency
                );
            }
            // Singular matrix B: reorder the basic variables so that the
            // singular columns are in the position corresponding to the
            // logical which replaces them
            self.build_handle_rank_deficiency(basic_index);
            self.build_mark_sing_c(basic_index);
            factor_timer.stop(FactorInvertDeficient, factor_timer_clock_pointer);
        }
        if incomplete_basis {
            // Completing the factorization is not relevant if the basis
            // matrix is incomplete, so clear any refactorization information
            // and return
            self.refactor_info_.clear();
            debug_assert!(!self.refactor_info_.r#use);
            let basic_index_rank_deficiency =
                self.rank_deficiency - (self.num_row - self.num_basic);
            return basic_index_rank_deficiency;
        }
        // Complete INVERT
        if report_lu {
            // printf("\nFactored INVERT\n"); reportLu(kReportLuBoth, false);
        }
        factor_timer.start(FactorInvertFinish, factor_timer_clock_pointer);
        self.build_finish(basic_index);
        factor_timer.stop(FactorInvertFinish, factor_timer_clock_pointer);
        //
        // Indicate that the refactorization information is known unless the basis was
        // rank deficient
        if self.rank_deficiency != 0 {
            self.refactor_info_.clear();
        } else {
            // Check that the refactorization information is not (yet) flagged
            // to be used in a future call
            debug_assert!(!self.refactor_info_.r#use);
            // Record build_synthetic_tick to use as the value of
            // build_synthetic_tick if refactorization is performed. Not only
            // is there no build_synthetic_tick measure for refactorization,
            // if there were it would give an unrealistic underestimate of the
            // cost of factorization from scratch
            self.refactor_info_.build_synthetic_tick = self.build_synthetic_tick;
        }

        // Record the number of entries in the INVERT
        self.invert_num_el = self.l_start[self.num_row as usize]
            + self.u_last_p[(self.num_row - 1) as usize]
            + self.num_row;

        self.kernel_dim -= self.rank_deficiency;
        debug_log_rank_deficiency(
            self.highs_debug_level,
            &self.log_options,
            self.rank_deficiency,
            self.basis_matrix_num_el,
            self.invert_num_el,
            self.kernel_dim,
            self.kernel_num_el,
            self.nwork,
        );
        factor_timer.stop(FactorInvert, factor_timer_clock_pointer);
        self.rank_deficiency
    }

    /// HiGHS: HFactor.cpp::HFactor::ftranCall (lines 463-473, 11 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    /// Overload `(HVector& vector, expected_density, clock)`.
    pub fn ftran_call(
        &self,
        vector: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        let use_indices = vector.count >= 0;
        let factor_timer = FactorTimer;
        factor_timer.start(FactorFtran, factor_timer_clock_pointer);
        self.ftran_l(vector, expected_density, factor_timer_clock_pointer);
        self.ftran_u(vector, expected_density, factor_timer_clock_pointer);
        // Possibly find the indices in order
        if use_indices {
            vector.re_index();
        }
        factor_timer.stop(FactorFtran, factor_timer_clock_pointer);
    }

    /// HiGHS: HFactor.cpp::HFactor::ftranCall (lines 475-489, 15 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    /// Overload `(std::vector<double>& vector, clock)`: uses the `rhs_` work vector.
    pub fn ftran_call_vec(
        &mut self,
        vector: &mut Vec<f64>,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        let factor_timer = FactorTimer;
        factor_timer.start(FactorFtran, factor_timer_clock_pointer);
        // Only have to clear the scalars of the local HVector, since the
        // array is moved in. Set the count to -1 to indicate that indices
        // aren't known.
        self.rhs_.clear_scalars();
        self.rhs_.array = std::mem::take(vector);
        self.rhs_.count = -1;
        let expected_density = 1.0;
        // Member passed to a method of the same object (crate root (f)).
        let mut rhs = std::mem::take(&mut self.rhs_);
        self.ftran_call(&mut rhs, expected_density, factor_timer_clock_pointer);
        self.rhs_ = rhs;
        *vector = std::mem::take(&mut self.rhs_.array);
        factor_timer.stop(FactorFtran, factor_timer_clock_pointer);
    }

    /// HiGHS: HFactor.cpp::HFactor::btranCall (lines 491-501, 11 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    /// Overload `(HVector& vector, expected_density, clock)`.
    pub fn btran_call(
        &self,
        vector: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        let use_indices = vector.count >= 0;
        let factor_timer = FactorTimer;
        factor_timer.start(FactorBtran, factor_timer_clock_pointer);
        self.btran_u(vector, expected_density, factor_timer_clock_pointer);
        self.btran_l(vector, expected_density, factor_timer_clock_pointer);
        // Possibly find the indices in order
        if use_indices {
            vector.re_index();
        }
        factor_timer.stop(FactorBtran, factor_timer_clock_pointer);
    }

    /// HiGHS: HFactor.cpp::HFactor::btranCall (lines 503-514, 12 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    /// Overload `(std::vector<double>& vector, clock)`: uses the `rhs_` work vector.
    pub fn btran_call_vec(
        &mut self,
        vector: &mut Vec<f64>,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        // Only have to clear the scalars of the local HVector, since the
        // array is moved in. Set the count to -1 to indicate that indices
        // aren't known.
        self.rhs_.clear_scalars();
        self.rhs_.array = std::mem::take(vector);
        self.rhs_.count = -1;
        let expected_density = 1.0;
        // Member passed to a method of the same object (crate root (f)).
        let mut rhs = std::mem::take(&mut self.rhs_);
        self.btran_call(&mut rhs, expected_density, factor_timer_clock_pointer);
        self.rhs_ = rhs;
        *vector = std::mem::take(&mut self.rhs_.array);
    }

    /// HiGHS: HFactor.cpp::HFactor::update (lines 516-531, 16 lines)
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    /// The C++ `if (aq->next)` -> `updateCFT` branch is PAMI-only: `HVector` has no `next`, so skip it.
    pub fn update(
        &mut self,
        basic_index: &mut [HighsInt],
        aq: &mut HVector,
        ep: &mut HVector,
        i_row: &mut HighsInt,
        hint: &mut HighsInt,
    ) {
        // Updating implies a change of basis. Since the refactorizaion info
        // no longer corresponds to the current basis, it must be
        // invalidated
        self.refactor_info_.clear();
        // Special case: `if (aq->next) { updateCFT(aq, ep, iRow); return; }` is the PAMI
        // chained update, `update_multi`; serial callers never set `next`.

        if self.update_method == kUpdateMethodFt {
            self.update_ft(aq, ep, *i_row);
        }
        if self.update_method == kUpdateMethodPf {
            self.update_pf(aq, *i_row, hint);
        }
        if self.update_method == kUpdateMethodMpf {
            self.update_mpf(aq, ep, *i_row, hint);
        }
        if self.update_method == kUpdateMethodApf {
            self.update_apf(basic_index, aq, ep, *i_row);
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::update (lines 516-528), with the C++ pointer arguments:
    /// `refactor_info_.clear()`, then `updateCFT` when `aq->next` is set, else the same
    /// dispatch as [`HFactor::update`] on `i_row[0]`. The PAMI entry point
    /// (`HEkkDual::majorUpdateFactor` via `HEkk::update_factor_multi`).
    ///
    /// # Safety
    /// As [`HFactor::update_cft`].
    pub unsafe fn update_multi(
        &mut self,
        basic_index: &mut [HighsInt],
        aq: *mut HVector,
        ep: *mut HVector,
        i_row: &mut [HighsInt],
        hint: &mut HighsInt,
    ) {
        // Updating implies a change of basis. Since the refactorizaion info
        // no longer corresponds to the current basis, it must be
        // invalidated
        self.refactor_info_.clear();

        // Special case: if aq->next is non-null, use the chained update
        // SAFETY: caller guarantees aq and ep are valid pointers to live HVectors
        // whose next chains are valid (PAMI protocol: vectors live in a fixed array
        // that does not move while the chain is in use).
        let aq_next = unsafe { (*aq).next };
        if aq_next.is_some() {
            // SAFETY: same protocol as above; update_cft walks the chain and
            // dereferences only within the caller's guarantee.
            unsafe {
                self.update_cft(aq, ep, i_row);
            }
            return;
        }

        if self.update_method == kUpdateMethodFt {
            self.update_ft(unsafe { &mut *aq }, unsafe { &mut *ep }, i_row[0]);
        }
        if self.update_method == kUpdateMethodPf {
            self.update_pf(unsafe { &mut *aq }, i_row[0], hint);
        }
        if self.update_method == kUpdateMethodMpf {
            self.update_mpf(unsafe { &mut *aq }, unsafe { &mut *ep }, i_row[0], hint);
        }
        if self.update_method == kUpdateMethodApf {
            self.update_apf(
                basic_index,
                unsafe { &mut *aq },
                unsafe { &mut *ep },
                i_row[0],
            );
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::setPivotThreshold (lines 533-538, 6 lines)
    /// C++ default arguments: new_pivot_threshold = kDefaultPivotThreshold
    pub fn set_pivot_threshold(&mut self, new_pivot_threshold: f64) -> bool {
        if new_pivot_threshold < kMinPivotThreshold {
            return false;
        }
        if new_pivot_threshold > kMaxPivotThreshold {
            return false;
        }
        self.pivot_threshold = new_pivot_threshold;
        true
    }

    /// HiGHS: HFactor.cpp::HFactor::setTimeLimit (lines 540-544, 5 lines)
    pub fn set_time_limit(&mut self, time_limit: f64) {
        self.time_limit_ = kHighsInf;
        if time_limit < 0.0 {
            return;
        }
        self.time_limit_ = time_limit;
    }

    /// HiGHS: HFactor.cpp::HFactor::luClear (lines 546-558, 13 lines)
    pub fn lu_clear(&mut self) {
        self.l_start.clear();
        self.l_start.push(0);
        self.l_index.clear();
        self.l_value.clear();

        self.u_pivot_index.clear();
        self.u_pivot_value.clear();
        self.u_start.clear();
        self.u_start.push(0);
        self.u_index.clear();
        self.u_value.clear();
    }
}
