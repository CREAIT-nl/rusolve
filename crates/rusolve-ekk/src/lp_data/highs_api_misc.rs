//! HiGHS source: lp_data/Highs.cpp (changeColsIntegrality interval,
//! getBasisTransposeSolve, optimizeModel, passColName, passRowName,
//! passModelName), lp_data/HighsInterface.cpp (changeIntegralityInterface,
//! basisSolveInterface, invertRequirementError), lp_data/HighsLpUtils.cpp
//! (changeLpIntegrality), Highs.h (zeroAllClocks, getPrimalPhase1Dual),
//! util/HighsTimer.h (HighsTimer::zeroAllClocks)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The remaining `Highs` API that `HighsIis::compute`,
//! `Highs::getIisInterface` and `Highs::elasticityFilter` call.
//! `highsVarTypeUserDataNotNull` is dropped (slices are never NULL);
//! `passColName` / `passRowName` keep the error checks and the
//! `col_names_` / `row_names_` writes, and drop the `col_hash_` /
//! `row_hash_` update (no name hash in this port).

use crate::prelude::*;
use std::time::Instant;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::changeColsIntegrality (interval overload,
    /// lines 2934-2955)
    pub fn change_cols_integrality_interval(
        &mut self,
        from_col: HighsInt,
        to_col: HighsInt,
        integrality: &[HighsVarType],
    ) -> HighsStatus {
        self.clear_presolve();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = create_interval(
            &mut index_collection,
            from_col,
            to_col,
            self.model_lp_.num_col_,
        );
        if create_error != 0 {
            highs_log_user(
                &self.options_.log_options,
                HighsLogType::kError,
                &format!(
                    "Interval [{}, {}] supplied to Highs::changeColsIntegrality is out of range [0, {})\n",
                    from_col, to_col, self.model_lp_.num_col_
                ),
            );
            return HighsStatus::kError;
        }
        let call_status = self.change_integrality_interface(&mut index_collection, integrality);
        let mut return_status = HighsStatus::kOk;
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "changeIntegrality",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: HighsInterface.cpp::Highs::changeIntegralityInterface (lines 883-906)
    pub(crate) fn change_integrality_interface(
        &mut self,
        index_collection: &mut HighsIndexCollection,
        integrality: &[HighsVarType],
    ) -> HighsStatus {
        let num_integrality = data_size(index_collection);
        // If a non-positive number of integrality (may) need changing nothing needs
        // to be done
        if num_integrality <= 0 {
            return HighsStatus::kOk;
        }
        // highsVarTypeUserDataNotNull dropped (slices are never NULL)
        // Take a copy of the integrality that can be normalised
        let local_integrality: Vec<HighsVarType> = integrality[..num_integrality as usize].to_vec();
        // If changing the integrality for a set of columns, verify that the
        // set entries are in ascending order
        if index_collection.is_set_ {
            debug_assert!(increasing_set_ok(
                &index_collection.set_,
                0,
                index_collection.dimension_,
                true
            ));
        }
        let mut return_status = change_lp_integrality(
            &mut self.model_lp_,
            index_collection,
            &local_integrality,
            &self.options_,
        );
        debug_assert!(return_status != HighsStatus::kError);
        // Deduce the consequences of new integrality
        self.invalidate_model_status();
        return_status
    }

    /// HiGHS: Highs.cpp::Highs::getBasisTransposeSolve (lines 2370-2395)
    ///
    /// PORT NOTE: `Xrhs` / `solution_vector` are slices (the NULL checks
    /// are dropped); `solution_num_nz` / `solution_indices` are
    /// `Option`s (NULL is the trigger for not identifying indices).
    pub fn get_basis_transpose_solve(
        &mut self,
        xrhs: &[f64],
        solution_vector: &mut [f64],
        solution_num_nz: Option<&mut HighsInt>,
        solution_indices: Option<&mut [HighsInt]>,
    ) -> HighsStatus {
        // Xrhs == NULL check dropped (slice)
        // solution_vector == NULL check dropped (slice)
        // solution_indices can be None - it's the trigger that determines
        // whether they are identified or not
        if !self.ekk_instance_.status_.has_invert {
            return self.invert_requirement_error("getBasisTransposeSolve");
        }
        let num_row = self.model_lp_.num_row_;
        let mut rhs: Vec<f64> = vec![0.0; num_row as usize];
        for row in 0..num_row {
            rhs[row as usize] = xrhs[row as usize];
        }
        self.basis_solve_interface(
            rhs.as_slice(),
            solution_vector,
            solution_num_nz,
            solution_indices,
            true,
        );
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::basisSolveInterface (lines 1480-1561)
    ///
    /// Use `NlaRefs` built from `self.model_lp_` for
    /// `ekk_instance_.btran/ftran`, as `get_basis_inverse_row_sparse`
    /// (highs_api.rs) does.
    pub(crate) fn basis_solve_interface(
        &mut self,
        rhs: &[f64],
        solution_vector: &mut [f64],
        solution_num_nz: Option<&mut HighsInt>,
        solution_indices: Option<&mut [HighsInt]>,
        transpose: bool,
    ) -> HighsStatus {
        let return_status = HighsStatus::kOk;
        let num_row = self.model_lp_.num_row_;
        // For an LP with no rows the solution is vacuous
        if num_row == 0 {
            return return_status;
        }
        // EKK must have an INVERT, but simplex NLA may need the pointer to
        // its LP to be refreshed so that it can use its scale factors
        debug_assert!(self.ekk_instance_.status_.has_invert);
        // Reset the simplex NLA LP and scale pointers for the unscaled LP:
        // setNlaPointersForLpAndScale(lp) -> NlaRefs built from model_lp_
        debug_assert!(!self.model_lp_.is_moved_);
        // Set up solve vector with suitably scaled RHS
        let mut solve_vector = HVector::default();
        solve_vector.setup(num_row);
        solve_vector.clear();
        let mut rhs_num_nz: HighsInt = 0;
        for i_row in 0..num_row {
            if rhs[i_row as usize] != 0.0 {
                solve_vector.index[rhs_num_nz as usize] = i_row;
                rhs_num_nz += 1;
                solve_vector.array[i_row as usize] = rhs[i_row as usize];
            }
        }
        solve_vector.count = rhs_num_nz;
        //
        // Note that solve_vector.count is just used to determine whether
        // hyper-sparse solves should be used. The indices of the nonzeros
        // in the solution are always accumulated. There's no switch (such
        // as setting solve_vector.count = num_row+1) to not do this.
        //
        // Get expected_density from analysis during simplex solve.
        let expected_density: f64 = 1.0;
        // HEkk::btran/ftran with simplex_nla_ pointing at the incumbent LP
        // (HEkk::btran in this port would use ekk_instance_.lp_ instead)
        debug_assert!(self.ekk_instance_.status_.has_nla);
        let refs = NlaRefs::new(
            &self.model_lp_,
            &self.ekk_instance_.basis_.basic_index_,
            &self.ekk_instance_.options_,
        );
        if transpose {
            self.ekk_instance_
                .simplex_nla_
                .btran(&refs, &mut solve_vector, expected_density, None);
        } else {
            self.ekk_instance_
                .simplex_nla_
                .ftran(&refs, &mut solve_vector, expected_density, None);
        }
        // Extract the solution
        let indices_required = solution_indices.is_some();
        if !indices_required {
            // Nonzeros in the solution not required
            if solve_vector.count > num_row {
                // Solution nonzeros not known
                for i_row in 0..num_row {
                    solution_vector[i_row as usize] = solve_vector.array[i_row as usize];
                }
            } else {
                // Solution nonzeros are known
                for i_row in 0..num_row {
                    solution_vector[i_row as usize] = 0.0;
                }
                for i_x in 0..solve_vector.count {
                    let i_row = solve_vector.index[i_x as usize];
                    solution_vector[i_row as usize] = solve_vector.array[i_row as usize];
                }
            }
        } else {
            // Nonzeros in the solution are required
            let solution_indices = solution_indices.expect("checked: solution_indices is Some");
            let solution_num_nz = solution_num_nz.expect("checked: solution_num_nz is Some");
            if solve_vector.count > num_row {
                // Solution nonzeros not known
                // C++ writes `solution_num_nz = 0` (nulls the pointer) and then
                // `solution_indices[*solution_num_nz++]`, which is UB; this
                // is the evident intent.
                *solution_num_nz = 0;
                for i_row in 0..num_row {
                    solution_vector[i_row as usize] = 0.0;
                    if solve_vector.array[i_row as usize] != 0.0 {
                        solution_vector[i_row as usize] = solve_vector.array[i_row as usize];
                        solution_indices[*solution_num_nz as usize] = i_row;
                        *solution_num_nz += 1;
                    }
                }
            } else {
                // Solution nonzeros are known
                for i_row in 0..num_row {
                    solution_vector[i_row as usize] = 0.0;
                }
                for i_x in 0..solve_vector.count {
                    let i_row = solve_vector.index[i_x as usize];
                    solution_vector[i_row as usize] = solve_vector.array[i_row as usize];
                    solution_indices[i_x as usize] = i_row;
                }
                *solution_num_nz = solve_vector.count;
            }
        }
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::invertRequirementError (lines 2953-2958)
    pub(crate) fn invert_requirement_error(&self, method_name: &str) -> HighsStatus {
        debug_assert!(!self.ekk_instance_.status_.has_invert);
        highs_log_user(
            &self.options_.log_options,
            HighsLogType::kError,
            &format!("No invertible representation for {}\n", method_name),
        );
        HighsStatus::kError
    }

    /// HiGHS: Highs.cpp::Highs::optimizeModel (lines 1153-1168)
    ///
    /// PORT NOTE: the profiling scaffolding (1159-1160, 1162-1165) is not ported.
    pub fn optimize_model(&mut self) -> HighsStatus {
        // Level 2a of Highs::run()
        //
        // C++ lines 1156-1157
        let status = self.initialize_multi_threading();
        if status != HighsStatus::kOk {
            return status;
        }

        // C++ lines 1158-1163: profiling scaffolding is not ported
        self.called_optimize_model()
    }

    /// HiGHS: Highs.cpp::Highs::initializeMultiThreading (lines 5111-5135):
    /// `highs::parallel::initialize_scheduler(options_.threads)`, `max_threads_ =
    /// num_threads()`, the "already been initialized" error when `threads != 0` and differs.
    ///
    pub fn initialize_multi_threading(&mut self) -> HighsStatus {
        // C++ line 5112: highs::parallel::initialize_scheduler(this->options_.threads);
        crate::parallel::highs_parallel::initialize_scheduler(self.options_.threads);
        // C++ line 5113: this->max_threads_ = highs::parallel::num_threads();
        self.max_threads_ = crate::parallel::highs_parallel::num_threads();
        // C++ line 5114: HighsLogOptions& log_options = this->options_.log_options;
        // (used inline below as &self.options_.log_options)
        // C++ lines 5115-5124:
        if self.options_.threads != 0 && self.max_threads_ != self.options_.threads {
            highs_log_user(
                &self.options_.log_options,
                HighsLogType::kError,
                &format!(
                    "Option 'threads' is set to {} but global scheduler has already been initialized to use {} threads. The previous scheduler instance can be destroyed by calling Highs::resetGlobalScheduler().\n",
                    self.options_.threads, self.max_threads_
                ),
            );
            return HighsStatus::kError;
        }
        // C++ lines 5125-5130:
        if self.max_threads_ <= 0 {
            highs_log_dev(
                &self.options_.log_options,
                HighsLogType::kError,
                &format!("max_threads() returns {}\n", self.max_threads_),
            );
            debug_assert!(self.max_threads_ > 0);
            return HighsStatus::kError;
        }
        // C++ lines 5131-5132:
        highs_log_dev(
            &self.options_.log_options,
            HighsLogType::kDetailed,
            &format!("Running with {} thread(s)\n", self.max_threads_),
        );
        // C++ line 5134: return HighsStatus::kOk;
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::resetGlobalScheduler (lines 5137-5139)
    /// (`static`, C++ default `blocking = false`): `HighsTaskExecutor::shutdown(blocking)`
    /// on the calling thread's executor.
    pub fn reset_global_scheduler(blocking: bool) {
        // C++ line 5138: HighsTaskExecutor::shutdown(blocking);
        crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(blocking);
    }

    /// rusolve: make the next top-level `run` honour `threads` even
    /// when an earlier run on this OS thread started a scheduler with another worker count.
    /// HiGHS answers that case with the `initialize_multi_threading` error and leaves the
    /// caller to call `resetGlobalScheduler`; a long-lived host (the MCP server, the Python
    /// binding) solves with `threads = 1` and then `threads = 4` on one thread, so it calls
    /// this before `run`. Never call it from inside a solve: it stops the worker threads of
    /// the running executor (sub-solvers inherit the running scheduler, as in C++).
    pub fn reset_scheduler_for_threads(&self) {
        // The reset stops the executor's workers: only its main thread (owner id 0), outside any
        // `run`, may do that.
        debug_assert!(
            !crate::lp_data::highs_run::in_run(),
            "reset_scheduler_for_threads called inside a solve"
        );
        debug_assert!(
            crate::parallel::highs_task_executor::HighsTaskExecutor::this_thread_is_main(),
            "reset_scheduler_for_threads called from an executor worker thread"
        );
        let wanted = self.options_.threads;
        if wanted != 0 && crate::parallel::highs_parallel::num_threads() != wanted {
            Self::reset_global_scheduler(true);
            // `Highs::run` reads the worker deque (LP fold, HiPO auto pick) before it reaches
            // `initialize_multi_threading`; start the new executor now so that read sees it.
            crate::parallel::highs_parallel::initialize_scheduler(wanted);
        }
    }

    /// HiGHS: Highs.h::Highs::zeroAllClocks (line 564)
    pub fn zero_all_clocks(&mut self) {
        self.timer_.zero_all_clocks();
    }

    /// HiGHS: Highs.h::Highs::getPrimalPhase1Dual (lines 1428-1430)
    pub fn get_primal_phase1_dual(&self) -> &[f64] {
        &self.ekk_instance_.primal_phase1_dual_
    }

    /// HiGHS: Highs.cpp::Highs::passColName (lines 743-762)
    pub fn pass_col_name(&mut self, col: HighsInt, name: &str) -> HighsStatus {
        let num_col = self.model_lp_.num_col_;
        if col < 0 || col >= num_col {
            highs_log_user(
                &self.options_.log_options,
                HighsLogType::kError,
                &format!(
                    "Index {} for column name {} is outside the range [0, num_col = {})\n",
                    col, name, num_col
                ),
            );
            return HighsStatus::kError;
        }
        if name.is_empty() {
            highs_log_user(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot define empty column names\n",
            );
            return HighsStatus::kError;
        }
        self.model_lp_
            .col_names_
            .resize(num_col as usize, String::new());
        // col_hash_.update dropped (no name hash in this port)
        self.model_lp_.col_names_[col as usize] = name.to_string();
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::passRowName (lines 764-783)
    pub fn pass_row_name(&mut self, row: HighsInt, name: &str) -> HighsStatus {
        let num_row = self.model_lp_.num_row_;
        if row < 0 || row >= num_row {
            highs_log_user(
                &self.options_.log_options,
                HighsLogType::kError,
                &format!(
                    "Index {} for row name {} is outside the range [0, num_row = {})\n",
                    row, name, num_row
                ),
            );
            return HighsStatus::kError;
        }
        if name.is_empty() {
            highs_log_user(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot define empty row names\n",
            );
            return HighsStatus::kError;
        }
        self.model_lp_
            .row_names_
            .resize(num_row as usize, String::new());
        // row_hash_.update dropped (no name hash in this port)
        self.model_lp_.row_names_[row as usize] = name.to_string();
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::passModelName (lines 785-788)
    pub fn pass_model_name(&mut self, name: &str) -> HighsStatus {
        self.model_lp_.model_name_ = name.to_string();
        HighsStatus::kOk
    }
}

/// HiGHS: HighsSort.cpp::increasingSetOk(const vector<HighsInt>&, HighsInt,
/// HighsInt, bool) (lines 261-287); used only by the debug assert in
/// `change_integrality_interface` (as inlined in util/highs_utils.rs::ok).
fn increasing_set_ok(
    set: &[HighsInt],
    set_entry_lower: HighsInt,
    set_entry_upper: HighsInt,
    strict: bool,
) -> bool {
    let set_num_entries = set.len() as HighsInt;
    let check_bounds = set_entry_lower <= set_entry_upper;
    let mut previous_entry = if check_bounds {
        if strict {
            set_entry_lower - 1
        } else {
            set_entry_lower
        }
    } else {
        -kHighsIInf
    };
    for k in 0..set_num_entries {
        let entry = set[k as usize];
        if strict {
            if entry <= previous_entry {
                return false;
            }
        } else {
            if entry < previous_entry {
                return false;
            }
        }
        if check_bounds && entry > set_entry_upper {
            return false;
        }
        previous_entry = entry;
    }
    true
}

/// HiGHS: HighsLpUtils.cpp::changeLpIntegrality (lines 1715-1762)
pub fn change_lp_integrality(
    lp: &mut HighsLp,
    index_collection: &HighsIndexCollection,
    new_integrality: &[HighsVarType],
    options: &HighsOptions,
) -> HighsStatus {
    let mut return_status = HighsStatus::kOk;
    debug_assert!(ok(index_collection));
    let mut from_k: HighsInt = 0;
    let mut to_k: HighsInt = 0;
    limits(index_collection, &mut from_k, &mut to_k);
    if from_k > to_k {
        return return_status;
    }

    let interval = index_collection.is_interval_;
    let mask = index_collection.is_mask_;
    let col_set = &index_collection.set_;
    let col_mask = &index_collection.mask_;

    // Change the integrality to the user-supplied integrality, according to the
    // technique
    let mut lp_col: HighsInt;
    let mut usr_col: HighsInt = -1;

    // If changing integrality for a problem without an integrality
    // vector (ie an LP), have to create it for the incumbent columns -
    // which are naturally continuous
    if lp.integrality_.len() == 0 {
        lp.integrality_
            .resize(lp.num_col_ as usize, HighsVarType::kContinuous);
    }

    debug_assert!(lp.integrality_.len() as HighsInt == lp.num_col_);
    let mut k = from_k;
    while k < to_k + 1 {
        if interval || mask {
            lp_col = k;
        } else {
            lp_col = col_set[k as usize];
        }
        let col = lp_col;
        if interval {
            usr_col += 1;
        } else {
            usr_col = k;
        }
        if mask && col_mask[col as usize] == 0 {
            k += 1;
            continue;
        }
        lp.integrality_[col as usize] = new_integrality[usr_col as usize];
        k += 1;
    }
    // If integrality_ contains only HighsVarType::kContinuous then
    // clear it
    if !lp.is_mip() {
        lp.integrality_.clear();
    }
    return_status
}

impl HighsTimer {
    /// HiGHS: HighsTimer.h::HighsTimer::zeroAllClocks (lines 90-96)
    ///
    /// PORT NOTE: `HighsTimer` is a single-clock wall stub (util/highs_timer.rs).
    /// Zeroing every clock leaves the run clock stopped with zero time, so
    /// the stub sets `start_time = Instant::now()` and `running = false`.
    pub fn zero_all_clocks(&mut self) {
        self.start_time = Instant::now();
        self.running = false;
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// Column-wise LP builder: `cols[j] = (cost, lower, upper, [(row, value)])`.
    #[allow(clippy::type_complexity)]
    fn lp(cols: &[(f64, f64, f64, &[(HighsInt, f64)])], rows: &[(f64, f64)]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = cols.len() as HighsInt;
        lp.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = vec![0];
        for &(cost, lower, upper, entries) in cols {
            lp.col_cost_.push(cost);
            lp.col_lower_.push(lower);
            lp.col_upper_.push(upper);
            for &(r, v) in entries {
                lp.a_matrix_.index_.push(r);
                lp.a_matrix_.value_.push(v);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        for &(l, u) in rows {
            lp.row_lower_.push(l);
            lp.row_upper_.push(u);
        }
        lp
    }

    fn quiet(lp: HighsLp) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs
    }

    fn iis_lp() -> HighsLp {
        lp(
            &[
                (0.0, 0.0, kHighsInf, &[(0, 1.0), (1, 1.0)]),
                (0.0, 0.0, kHighsInf, &[(0, 1.0), (1, -1.0)]),
            ],
            &[(-kHighsInf, 1.0), (2.0, kHighsInf)],
        )
    }

    #[test]
    fn integrality_and_names() {
        let mut highs = quiet(iis_lp());
        assert_eq!(
            highs.change_cols_integrality_interval(1, 1, &[HighsVarType::kInteger]),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.model_lp_.integrality_,
            vec![HighsVarType::kContinuous, HighsVarType::kInteger]
        );
        assert_eq!(highs.pass_model_name("w1"), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.model_name_, "w1");
        assert_eq!(highs.pass_col_name(0, "x0"), HighsStatus::kOk);
        assert_eq!(highs.pass_row_name(1, "r1"), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.row_names_[1], "r1");
    }

    #[test]
    fn basis_solve_needs_invert() {
        let mut highs = quiet(iis_lp());
        let mut sol = vec![0.0; 2];
        let mut num_nz = 0;
        assert_eq!(
            highs.get_basis_transpose_solve(&[1.0, 0.0], &mut sol, Some(&mut num_nz), None),
            HighsStatus::kError
        );
    }

    // Unit tests

    #[test]
    fn w5_highs_mt_initialize_multi_threading() {
        // threads=0 -> kOk and max_threads_ == 1 (HIGHS_NO_DEFAULT_THREADS)
        let handle = std::thread::spawn(|| {
            let mut highs = Highs::default();
            highs.options_.threads = 0;
            let status = highs.initialize_multi_threading();
            let max_threads = highs.max_threads_;
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            (status, max_threads)
        });
        let (status, max_threads) = handle.join().unwrap();
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(max_threads, 1);

        // threads=1 -> kOk, max_threads_ == 1
        let handle = std::thread::spawn(|| {
            let mut highs = Highs::default();
            highs.options_.threads = 1;
            let status = highs.initialize_multi_threading();
            let max_threads = highs.max_threads_;
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            (status, max_threads)
        });
        let (status, max_threads) = handle.join().unwrap();
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(max_threads, 1);

        // threads=2 on a fresh thread -> kOk, max_threads_ == 2
        let handle = std::thread::spawn(|| {
            let mut highs = Highs::default();
            highs.options_.threads = 2;
            let status = highs.initialize_multi_threading();
            let max_threads = highs.max_threads_;
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            (status, max_threads)
        });
        let (status, max_threads) = handle.join().unwrap();
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(max_threads, 2);

        // threads=1 after a 2-thread init on the same thread -> kError
        let handle = std::thread::spawn(|| {
            let mut highs = Highs::default();
            highs.options_.threads = 2;
            let status1 = highs.initialize_multi_threading();
            assert_eq!(status1, HighsStatus::kOk);
            assert_eq!(highs.max_threads_, 2);

            highs.options_.threads = 1;
            let status2 = highs.initialize_multi_threading();
            let max_threads = highs.max_threads_;
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            (status2, max_threads)
        });
        let (status, max_threads) = handle.join().unwrap();
        assert_eq!(status, HighsStatus::kError);
        assert_eq!(max_threads, 2);

        // rusolve: reset_scheduler_for_threads lets a later run change the worker count
        let handle = std::thread::spawn(|| {
            let mut highs = Highs::default();
            highs.options_.threads = 2;
            assert_eq!(highs.initialize_multi_threading(), HighsStatus::kOk);
            highs.options_.threads = 1;
            highs.reset_scheduler_for_threads();
            // a read between the reset and `run` (LP fold, HiPO pick) must see the new executor
            assert_eq!(crate::parallel::highs_parallel::num_threads(), 1);
            let status = highs.initialize_multi_threading();
            let max_threads = highs.max_threads_;
            highs.options_.threads = 3;
            highs.reset_scheduler_for_threads();
            assert_eq!(crate::parallel::highs_parallel::num_threads(), 3);
            let status3 = highs.initialize_multi_threading();
            let max3 = highs.max_threads_;
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            (status, max_threads, status3, max3)
        });
        let (status, max_threads, status3, max3) = handle.join().unwrap();
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(max_threads, 1);
        assert_eq!(status3, HighsStatus::kOk);
        assert_eq!(max3, 3);
    }
}
