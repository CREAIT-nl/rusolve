//! HiGHS source: mip/HighsMipSolverData.cpp (lines 1-830)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;
use crate::util::highs_utils::{cpp_max, cpp_min, fractionality};

impl HighsMipSolver {
    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::HighsMipSolverData (lines 22-85)
    /// C++: `HighsMipSolverData(HighsMipSolver & mipsolver)`
    ///
    /// PORT NOTE (architecture): the C++ `mipdata_ = decltype(mipdata_)(new
    /// HighsMipSolverData(*this))` is `self.new_data()`: the members are built from
    /// `self` (which still holds the old `mipdata_`, as in C++ while the new object is
    /// constructed; no member constructor reads `mipsolver.mipdata_`), then installed,
    /// then the constructor body runs on the installed data (it registers the pools with
    /// `getDomain()` through `ms`). Members in declaration order (h:70-277); every
    /// initialiser that differs from the Rust `Default` is set explicitly.
    pub fn new_data(&mut self) {
        // lps(1, HighsLpRelaxation(mipsolver)): deque fill ctor = copy of a prototype
        let lp_proto = HighsLpRelaxation::new(self);
        let lp = HighsLpRelaxation::copy_from(self, &lp_proto);
        // domains(1, HighsDomain(mipsolver)): copy of a prototype
        let dom_proto = HighsDomain::new(self);
        let dom = HighsDomain::copy_from(self, &dom_proto);
        let heuristics = HighsPrimalHeuristics::new(self);
        let mut cliquetable = HighsCliqueTable::new(self.num_col());
        // HiGHS 1.16 #3540 (h16_clique_links): the flag is fixed at construction
        cliquetable.set_h16_flags(&self.options_mip_);
        let implications = HighsImplications::new(self);
        let objective_function = HighsObjectiveFunction::new(self);
        let debug_solution = HighsDebugSol::new(self);
        *self.mipdata_ = HighsMipSolverData {
            lps: vec![Slot::new(lp)],
            domains: vec![Slot::new(dom)],
            // pseudocosts(1): one default-constructed HighsPseudocost
            pseudocosts: vec![Slot::new(HighsPseudocost::default())],
            parallel_lock: false,
            heuristics: Slot::new(heuristics),
            cliquetable: Slot::new(cliquetable),
            implications: Slot::new(implications),
            objective_function: Slot::new(objective_function),
            presolve_status: HighsPresolveStatus::kNotSet,
            cliques_extracted: false,
            row_matrix_set: false,
            analytic_center_computed: false,
            analytic_center_status: HighsModelStatus::kNotset,
            detect_symmetries: false,
            firstlpsolobj: -kHighsInf,
            rootlpsolobj: -kHighsInf,
            lower_bound: -kHighsInf,
            upper_bound: kHighsInf,
            upper_limit: kHighsInf,
            optimality_limit: kHighsInf,
            debug_solution,
            ..Default::default()
        };
        // constructor body (cpp 76-84)
        let agelim = self.options_mip_.mip_pool_age_limit;
        let softlimit = self.options_mip_.mip_pool_soft_limit;
        self.mipdata_
            .conflictpools
            .push(Slot::new(HighsConflictPool::new(5 * agelim, softlimit)));
        let ncol = self.num_col();
        // HiGHS 1.16 #3574 (h16_cutpool_age): HighsMipSolverData.cpp lines 81 (1.16.0): the flag is read
        // from options_mip_ at construction and is passed to the new constructor.
        let h16_cutpool_age = self.options_mip_.h16_cutpool_age;
        self.mipdata_
            .cutpools
            .push(Slot::new(HighsCutPool::new_with_age_flag(
                ncol,
                agelim,
                softlimit,
                0,
                h16_cutpool_age,
            )));
        detach!(self.mipdata_.domains[0], |d| d.add_cutpool(self, 0));
        detach!(self.mipdata_.domains[0], |d| d.add_conflict_pool(self, 0));
        let allow = !self.submip;
        self.mipdata_.cliquetable.set_allow_parallel(allow);
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::solutionSourceToString (lines 87-162)
    /// C++: `std::string solutionSourceToString(const int solution_source, const bool code = true) const`
    /// C++ default arguments: `code = true`
    pub fn solution_source_to_string(&self, solution_source: i32, code: bool) -> String {
        // C++ returns from each branch; `(code, text)` pairs keep the same mapping.
        // kSolutionSourceInitial is commented out in the C++ (cpp 91-93) and not listed.
        let (code_str, text): (&str, &str) = if solution_source == kSolutionSourceNone {
            (" ", "None")
        } else if solution_source == kSolutionSourceBranching {
            ("B", "Branching")
        } else if solution_source == kSolutionSourceCentralRounding {
            ("C", "Central rounding")
        } else if solution_source == kSolutionSourceFeasibilityPump {
            ("F", "Feasibility pump")
        } else if solution_source == kSolutionSourceHeuristic {
            ("H", "Heuristic")
        } else if solution_source == kSolutionSourceShifting {
            ("I", "Shifting")
        } else if solution_source == kSolutionSourceFeasibilityJump {
            ("J", "Feasibility jump")
        } else if solution_source == kSolutionSourceSubMip {
            ("L", "Sub-MIP")
        } else if solution_source == kSolutionSourceEmptyMip {
            ("P", "Empty MIP")
        } else if solution_source == kSolutionSourceRandomizedRounding {
            ("R", "Randomized rounding")
        } else if solution_source == kSolutionSourceSolveLp {
            ("S", "Solve LP")
        } else if solution_source == kSolutionSourceEvaluateNode {
            ("T", "Evaluate node")
        } else if solution_source == kSolutionSourceUnbounded {
            ("U", "Unbounded")
        } else if solution_source == kSolutionSourceUserSolution {
            ("X", "User solution")
        } else if solution_source == kSolutionSourceHighsSolution {
            ("Y", "HiGHS solution")
        } else if solution_source == kSolutionSourceZiRound {
            ("Z", "ZI Round")
        } else if solution_source == kSolutionSourceTrivialZ {
            ("z", "Trivial zero")
        } else if solution_source == kSolutionSourceTrivialL {
            ("l", "Trivial lower")
        } else if solution_source == kSolutionSourceTrivialU {
            ("u", "Trivial upper")
        } else if solution_source == kSolutionSourceTrivialP {
            ("p", "Trivial point")
        } else if solution_source == kSolutionSourceCleanup {
            (" ", "")
        // mip_heuristic_run_fpr / mip_heuristic_run_local_mip: mip-heuristics apply_patch.cmake lines 766-767 at 36f7d9d
        } else if solution_source == kSolutionSourceFPR {
            ("A", "FPR")
        } else if solution_source == kSolutionSourceLocalMIP {
            ("M", "Local MIP")
        // rusolve mip_alns (Hendel 2022)
        } else if solution_source == kSolutionSourceAlns {
            ("N", "ALNS")
        // rusolve mip_shift_and_propagate (Berthold, Hendel 2015)
        } else if solution_source == kSolutionSourceShiftProp {
            ("Q", "Shift and propagate")
        } else if solution_source == kSolutionSourceProximity {
            ("O", "Proximity")
        } else if solution_source == kSolutionSourceLocalBranching {
            ("K", "Local branching")
        } else {
            println!(
                "HighsMipSolverData::solutionSourceToString: Unknown source = {}",
                solution_source
            );
            debug_assert_eq!(0, 111, "HighsMipSolverData::solutionSourceToString");
            ("*", "None")
        };
        if code {
            code_str.to_string()
        } else {
            text.to_string()
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::checkSolution (lines 164-187)
    /// C++: `bool checkSolution(const std::vector<double> & solution) const`
    pub fn check_solution(&self, solution: &[f64]) -> bool {
        let feastol = self.mipdata_.feastol;
        for i in 0..self.num_col() {
            let iu = i as usize;
            if solution[iu] < mip_model!(self).col_lower_[iu] - feastol {
                return false;
            }
            if solution[iu] > mip_model!(self).col_upper_[iu] + feastol {
                return false;
            }
            if self.is_col_integer(i) && fractionality(solution[iu], None) > feastol {
                return false;
            }
        }

        for i in 0..self.num_row() {
            let iu = i as usize;
            let mut rowactivity = 0.0;

            let start = self.mipdata_.a_rstart_[iu];
            let end = self.mipdata_.a_rstart_[iu + 1];

            for j in start..end {
                let ju = j as usize;
                rowactivity +=
                    solution[self.mipdata_.a_rindex_[ju] as usize] * self.mipdata_.a_rvalue_[ju];
            }

            if rowactivity > self.row_upper_2(i) + feastol {
                return false;
            }
            if rowactivity < self.row_lower_2(i) - feastol {
                return false;
            }
        }

        true
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::getInfeasibleRows (lines 189-213)
    /// C++: `std::vector<std::tuple<HighsInt, HighsInt, double>> getInfeasibleRows(const std::vector<double> & solution) const`
    pub fn get_infeasible_rows(&self, solution: &[f64]) -> Vec<(HighsInt, HighsInt, f64)> {
        let feastol = self.mipdata_.feastol;
        let mut infeasible_rows: Vec<(HighsInt, HighsInt, f64)> = Vec::new();
        for i in 0..self.num_row() {
            let iu = i as usize;
            let start = self.mipdata_.a_rstart_[iu];
            let end = self.mipdata_.a_rstart_[iu + 1];

            let mut row_activity_quad = HighsCDouble::from(0.0);
            for j in start..end {
                let ju = j as usize;
                row_activity_quad +=
                    HighsCDouble::from(solution[self.mipdata_.a_rindex_[ju] as usize])
                        * self.mipdata_.a_rvalue_[ju];
            }

            let row_activity = f64::from(row_activity_quad);
            if row_activity > self.row_upper_2(i) + feastol {
                let difference = (row_activity - self.row_upper_2(i)).abs();
                infeasible_rows.push((i, 1, difference));
            }
            if row_activity < self.row_lower_2(i) - feastol {
                let difference = (self.row_lower_2(i) - row_activity).abs();
                infeasible_rows.push((i, -1, difference));
            }
        }
        infeasible_rows
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::trySolution (lines 215-244)
    /// C++: `bool trySolution(const std::vector<double> & solution, const int solution_source = kSolutionSourceNone)`
    /// C++ default arguments: `solution_source = kSolutionSourceNone`
    pub fn try_solution(&mut self, solution: &[f64], solution_source: i32) -> bool {
        if solution.len() as HighsInt != self.num_col() {
            return false;
        }

        let feastol = self.mipdata_.feastol;
        let mut obj = HighsCDouble::from(0.0);
        // HiGHS 1.16 #3149 (h16_cdouble): HighsMipSolverData.cpp::trySolution lines 230 (1.16.0)
        let h16_cdouble_obj = self.options_mip_.h16_cdouble;

        for i in 0..self.num_col() {
            let iu = i as usize;
            if solution[iu] < mip_model!(self).col_lower_[iu] - feastol {
                return false;
            }
            if solution[iu] > mip_model!(self).col_upper_[iu] + feastol {
                return false;
            }
            if self.is_col_integer(i) && fractionality(solution[iu], None) > feastol {
                return false;
            }

            if h16_cdouble_obj {
                obj += HighsCDouble::from(self.col_cost_2(i)) * solution[iu];
            } else {
                obj += self.col_cost_2(i) * solution[iu];
            }
        }

        for i in 0..self.num_row() {
            let iu = i as usize;
            let mut rowactivity = 0.0;

            let start = self.mipdata_.a_rstart_[iu];
            let end = self.mipdata_.a_rstart_[iu + 1];

            for j in start..end {
                let ju = j as usize;
                rowactivity +=
                    solution[self.mipdata_.a_rindex_[ju] as usize] * self.mipdata_.a_rvalue_[ju];
            }

            if rowactivity > self.row_upper_2(i) + feastol {
                return false;
            }
            if rowactivity < self.row_lower_2(i) - feastol {
                return false;
            }
        }

        let obj = f64::from(obj);
        self.add_incumbent(solution, obj, solution_source, true, false)
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::solutionRowFeasible (lines 246-263)
    /// C++: `bool solutionRowFeasible(const std::vector<double> & solution) const`
    pub fn solution_row_feasible(&self, solution: &[f64]) -> bool {
        let feastol = self.mipdata_.feastol;
        for i in 0..self.num_row() {
            let iu = i as usize;
            let mut c_double_rowactivity = HighsCDouble::from(0.0);

            let start = self.mipdata_.a_rstart_[iu];
            let end = self.mipdata_.a_rstart_[iu + 1];

            for j in start..end {
                let ju = j as usize;
                c_double_rowactivity +=
                    HighsCDouble::from(solution[self.mipdata_.a_rindex_[ju] as usize])
                        * self.mipdata_.a_rvalue_[ju];
            }

            let rowactivity = f64::from(c_double_rowactivity);
            if rowactivity > self.row_upper_2(i) + feastol {
                return false;
            }
            if rowactivity < self.row_lower_2(i) - feastol {
                return false;
            }
        }
        true
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::trivialHeuristics (lines 265-407)
    /// C++: `HighsModelStatus trivialHeuristics()`
    #[allow(clippy::needless_range_loop)]
    pub fn trivial_heuristics(&mut self) -> HighsModelStatus {
        // HiGHS quirk: the C++ `printf` diagnostics are behind `lc_report = false`, and
        // `save_upper_bound` is only read by them; both are omitted (no behaviour).
        if !self.mipdata_.continuous_cols.is_empty() {
            return HighsModelStatus::kNotset;
        }
        let num_try_heuristic: HighsInt = 4;
        let heuristic_source: [i32; 4] = [
            kSolutionSourceTrivialZ,
            kSolutionSourceTrivialL,
            kSolutionSourceTrivialU,
            kSolutionSourceTrivialP,
        ];

        // C++ copies these (`std::vector<double> col_lower = ...`); `solution` and the
        // incumbent calls below need `&mut self`, so the row / matrix data are copies too.
        let mut col_lower = mip_model!(self).col_lower_.clone();
        let mut col_upper = mip_model!(self).col_upper_.clone();
        let row_lower = mip_model!(self).row_lower_.clone();
        let row_upper = mip_model!(self).row_upper_.clone();
        let matrix_start = mip_model!(self).a_matrix_.start_.clone();
        let matrix_value = mip_model!(self).a_matrix_.value_.clone();

        // Determine the following properties, according to which some
        // trivial heuristics are duplicated or fail immediately
        let mut all_integer_lower_non_positive = true;
        let mut all_integer_lower_zero = true;
        let mut all_integer_lower_finite = true;
        let mut all_integer_upper_finite = true;
        for integer_col in 0..self.mipdata_.numintegercols {
            let i_col = self.mipdata_.integer_cols[integer_col as usize];
            let icu = i_col as usize;
            // Round bounds in to nearest integer
            col_lower[icu] = col_lower[icu].ceil();
            col_upper[icu] = col_upper[icu].floor();
            let legal_bounds = col_lower[icu] <= col_upper[icu]
                && col_lower[icu] < kHighsInf
                && col_upper[icu] > -kHighsInf
                && !col_lower[icu].is_nan()
                && !col_upper[icu].is_nan();
            if !legal_bounds {
                debug_assert!(legal_bounds);
                crate::highs_log_user!(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "HighsMipSolverData::trivialHeuristics() has detected infeasible/illegal bounds [{}, {}] for column {}: MIP is infeasible\n",
                    col_lower[icu],
                    col_upper[icu],
                    i_col
                );
                return HighsModelStatus::kInfeasible;
            }
            // If bounds are inconsistent then MIP is infeasible
            if col_lower[icu] > col_upper[icu] {
                return HighsModelStatus::kInfeasible;
            }

            if col_lower[icu] > 0.0 {
                all_integer_lower_non_positive = false;
            }
            if col_lower[icu] != 0.0 {
                all_integer_lower_zero = false;
            }
            if col_lower[icu] <= -kHighsInf {
                all_integer_lower_finite = false;
            }
            if col_upper[icu] >= kHighsInf {
                all_integer_upper_finite = false;
            }
            // Only continue if one of the properties still holds
            if !(all_integer_lower_non_positive
                || all_integer_lower_zero
                || all_integer_upper_finite)
            {
                break;
            }
        }
        let all_integer_boxed = all_integer_lower_finite && all_integer_upper_finite;

        let feasibility_tolerance = self.options_mip_.mip_feasibility_tolerance;
        // Loop through the trivial heuristics
        let mut solution: Vec<f64> = vec![0.0; self.num_col() as usize];
        for try_heuristic in 0..num_try_heuristic {
            if try_heuristic == 0 {
                // First heuristic is to see whether all-zero for integer
                // variables is feasible
                //
                // If there is a positive lower bound then the heuristic fails
                if !all_integer_lower_non_positive {
                    continue;
                }
                // Determine whether a zero row activity is feasible
                let mut heuristic_failed = false;
                for i_row in 0..self.num_row() {
                    let ir = i_row as usize;
                    if row_lower[ir] > feasibility_tolerance
                        || row_upper[ir] < -feasibility_tolerance
                    {
                        heuristic_failed = true;
                        break;
                    }
                }
                if heuristic_failed {
                    continue;
                }
                solution = vec![0.0; self.num_col() as usize];
            } else if try_heuristic == 1 {
                // Second heuristic is to see whether all-lower for integer
                // variables (if distinct from all-zero) is feasible
                if all_integer_lower_zero {
                    continue;
                }
                // Trivially feasible for columns
                if !self.solution_row_feasible(&col_lower) {
                    continue;
                }
                solution = col_lower.clone();
            } else if try_heuristic == 2 {
                // Third heuristic is to see whether all-upper for integer
                // variables is feasible
                //
                // If there is an infinite upper bound then the heuristic fails
                if !all_integer_upper_finite {
                    continue;
                }
                // Trivially feasible for columns
                if !self.solution_row_feasible(&col_upper) {
                    continue;
                }
                solution = col_upper.clone();
            } else if try_heuristic == 3 {
                // Fourth heuristic is to see whether the "lock point" is feasible
                if !all_integer_boxed {
                    continue;
                }
                for integer_col in 0..self.mipdata_.numintegercols {
                    let i_col = self.mipdata_.integer_cols[integer_col as usize];
                    let icu = i_col as usize;
                    let mut num_positive_values: HighsInt = 0;
                    let mut num_negative_values: HighsInt = 0;
                    for i_el in matrix_start[icu]..matrix_start[icu + 1] {
                        if matrix_value[i_el as usize] > 0.0 {
                            num_positive_values += 1;
                        } else {
                            num_negative_values += 1;
                        }
                    }
                    solution[icu] = if num_positive_values > num_negative_values {
                        col_lower[icu]
                    } else {
                        col_upper[icu]
                    };
                }
                // Trivially feasible for columns
                if !self.solution_row_feasible(&solution) {
                    continue;
                }
            }

            let mut cdouble_obj = HighsCDouble::from(0.0);
            // HiGHS 1.16 #3149 (h16_cdouble): HighsMipSolverData.cpp::trivialHeuristics lines 391-395 (1.16.0)
            let h16_cdouble_obj = self.options_mip_.h16_cdouble;
            for i_col in 0..self.num_col() {
                if h16_cdouble_obj {
                    cdouble_obj +=
                        HighsCDouble::from(self.col_cost_2(i_col)) * solution[i_col as usize];
                } else {
                    cdouble_obj += self.col_cost_2(i_col) * solution[i_col as usize];
                }
            }
            let obj = f64::from(cdouble_obj);
            self.add_incumbent(
                &solution,
                obj,
                heuristic_source[try_heuristic as usize],
                true,
                false,
            );
        }
        HighsModelStatus::kNotset
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::startAnalyticCenterComputation (lines 409-487)
    /// C++: `void startAnalyticCenterComputation(const highs::parallel::TaskGroup & taskGroup)`
    ///
    /// The whole body is `taskGroup.spawn([&]() {..})`. The `[&]` lambda captures `this`; it
    /// is [`AnalyticCenterTask`], which holds field-level pointers to what the lambda reads
    /// and writes (rule (s)), taken here before the spawn.
    pub fn start_analytic_center_computation(&mut self, task_group: &HighsMipTaskGroup<'_, '_>) {
        let task = AnalyticCenterTask {
            options_mip_: std::ptr::addr_of!(self.options_mip_),
            model_: std::ptr::addr_of!(self.model_),
            orig_model_: std::ptr::addr_of!(self.orig_model_),
            presolved_model: std::ptr::addr_of!(self.mipdata_.presolved_model),
            submip_level: self.submip_level,
            // HiGHS 1.16 #3456 / #3522 (h16_ac_centring): the option is read once, at the spawn
            h16_ac_centring: self.options_mip_.h16_ac_centring,
            analytic_center_status: std::ptr::addr_of_mut!(self.mipdata_.analytic_center_status),
            analytic_center: std::ptr::addr_of_mut!(self.mipdata_.analytic_center),
        };
        task_group.spawn(move || {
            // SAFETY: see `AnalyticCenterTask::run`.
            unsafe { task.run() }
        });
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::finishAnalyticCenterComputation (lines 489-541)
    /// HiGHS 1.16 (#3408): HighsMipSolverData.cpp::finishAnalyticCenterComputation (lines 489-545)
    /// C++: `void finishAnalyticCenterComputation(const highs::parallel::TaskGroup & taskGroup)`
    ///
    #[allow(clippy::needless_return)] // C++ `return` at the end of the body kept for control-flow parity
    pub fn finish_analytic_center_computation(&mut self, task_group: &HighsMipTaskGroup<'_, '_>) {
        // PORT NOTE: the C++ `profiling_->mip_` MIP-Timing lines (cpp 491-494, 500-503) are
        // omitted: `HighsProfiling` is a no-op stub whose `mip_` flag is off by default.
        task_group.sync();

        self.mipdata_.analytic_center_computed = true;
        // HiGHS 1.16 #3408 (h16_ac_fix): HighsMipSolverData.cpp lines 505-509 (1.16.0)
        let h16_ac_fix = self.options_mip_.h16_ac_fix;
        let fix_at_centre = if h16_ac_fix {
            self.mipdata_.analytic_center_status == HighsModelStatus::kOptimal
                && self.options_mip_.presolve != kHighsOffString
        } else {
            self.mipdata_.analytic_center_status == HighsModelStatus::kOptimal
        };
        if fix_at_centre {
            let mut nfixed: HighsInt = 0;
            let mut nintfixed: HighsInt = 0;
            for i in 0..self.num_col() {
                let iu = i as usize;
                let bound_range = self.mipdata_.domains[0].col_upper_[iu]
                    - self.mipdata_.domains[0].col_lower_[iu];
                if bound_range == 0.0 {
                    continue;
                }

                // HiGHS 1.16 #3408 (h16_ac_fix): HighsMipSolverData.cpp lines 517-518 (1.16.0)
                let tolerance = if h16_ac_fix {
                    (self.mipdata_.feastol / 10.0) * cpp_min(bound_range, 1.0)
                } else {
                    self.mipdata_.feastol * cpp_min(bound_range, 1.0)
                };

                if self.mipdata_.analytic_center[iu] <= mip_model!(self).col_lower_[iu] + tolerance
                {
                    let lower = mip_model!(self).col_lower_[iu];
                    detach!(self.mipdata_.domains[0], |d| d.change_bound_2(
                        self,
                        HighsBoundType::kUpper,
                        i,
                        lower,
                        HighsDomain_Reason::unspecified()
                    ));
                    if self.mipdata_.domains[0].infeasible() {
                        return;
                    }
                    nfixed += 1;
                    if self.is_col_integer(i) {
                        nintfixed += 1;
                    }
                } else if self.mipdata_.analytic_center[iu]
                    >= mip_model!(self).col_upper_[iu] - tolerance
                {
                    let upper = mip_model!(self).col_upper_[iu];
                    detach!(self.mipdata_.domains[0], |d| d.change_bound_2(
                        self,
                        HighsBoundType::kLower,
                        i,
                        upper,
                        HighsDomain_Reason::unspecified()
                    ));
                    if self.mipdata_.domains[0].infeasible() {
                        return;
                    }
                    nfixed += 1;
                    if self.is_col_integer(i) {
                        nintfixed += 1;
                    }
                }
            }
            if nfixed > 0 {
                crate::highs_log_dev!(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "Fixing {} columns ({} integers) sitting at bound at analytic center\n",
                    nfixed,
                    nintfixed
                );
            }
            detach!(self.mipdata_.domains[0], |d| d.propagate(self));
            if self.mipdata_.domains[0].infeasible() {
                return;
            }
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::startSymmetryDetection (lines 543-560)
    /// C++: `void startSymmetryDetection(const highs::parallel::TaskGroup & taskGroup, std::unique_ptr<SymmetryDetectionData> & symData)`
    ///
    /// The `[&]` lambda passed to `taskGroup.spawn` (captures `this` and the caller's
    /// `symData` by reference) is [`SymmetryDetectionTask`], with field-level pointers
    /// taken before the spawn (rule (s)).
    pub fn start_symmetry_detection(
        &mut self,
        task_group: &HighsMipTaskGroup<'_, '_>,
        sym_data: &mut Option<Box<HighsMipSolverData_SymmetryDetectionData>>,
    ) {
        *sym_data = Some(Box::new(HighsMipSolverData_SymmetryDetectionData::default()));
        let sd = sym_data.as_deref_mut().expect("symData was just created");
        // HiGHS 1.16 #3182 (h16_sym_refine): flag carried into the detection struct before the spawn
        sd.sym_detection.h16_sym_refine = self.options_mip_.h16_sym_refine;
        // HiGHS #3535 (mip_sym_hash64): likewise, with the time limit the 1.16.0 task passes
        sd.sym_detection.mip_sym_hash64 = self.options_mip_.mip_sym_hash64;
        sd.sym_detection.load_model_as_graph(
            &self.mipdata_.presolved_model,
            self.options_mip_.small_matrix_value,
        );
        self.mipdata_.detect_symmetries = sd.sym_detection.initialize_detection();

        if self.mipdata_.detect_symmetries {
            let task = SymmetryDetectionTask {
                timer_: std::ptr::addr_of!(self.timer_),
                time_limit: self.options_mip_.time_limit,
                presolved_model: std::ptr::addr_of!(self.mipdata_.presolved_model),
                sym_data: std::ptr::addr_of_mut!(*sd),
            };
            task_group.spawn(move || {
                // SAFETY: see `SymmetryDetectionTask::run`.
                unsafe { task.run() }
            });
        } else {
            *sym_data = None;
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::finishSymmetryDetection (lines 562-606)
    /// C++: `void finishSymmetryDetection(const highs::parallel::TaskGroup & taskGroup, std::unique_ptr<SymmetryDetectionData> & symData)`
    pub fn finish_symmetry_detection(
        &mut self,
        task_group: &HighsMipTaskGroup<'_, '_>,
        sym_data: &mut Option<Box<HighsMipSolverData_SymmetryDetectionData>>,
    ) {
        task_group.sync();

        let sd = sym_data
            .as_deref_mut()
            .expect("finishSymmetryDetection without symData");
        self.mipdata_.symmetries = std::mem::take(&mut sd.symmetries);
        let detection_time = sd.detection_time;
        let symmetry_time = if self.options_mip_.timeless_log {
            String::new()
        } else {
            format!(" {:.1}s", detection_time)
        };
        crate::highs_log_user!(
            &self.options_mip_.log_options,
            HighsLogType::kInfo,
            "\nSymmetry detection completed in{}\n",
            symmetry_time
        );

        if self.mipdata_.symmetries.num_generators == 0 {
            self.mipdata_.detect_symmetries = false;
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "No symmetry present\n\n"
            );
        } else if self.mipdata_.symmetries.orbitopes.is_empty() {
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "Found {} generator(s)\n\n",
                self.mipdata_.symmetries.num_generators
            );
        } else if self.mipdata_.symmetries.num_perms != 0 {
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "Found {} generator(s) and {} full orbitope(s) acting on {} columns\n\n",
                self.mipdata_.symmetries.num_perms,
                self.mipdata_.symmetries.orbitopes.len(),
                self.mipdata_.symmetries.column_to_orbitope.size()
            );
        } else {
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "Found {} full orbitope(s) acting on {} columns\n\n",
                self.mipdata_.symmetries.orbitopes.len(),
                self.mipdata_.symmetries.column_to_orbitope.size()
            );
        }
        *sym_data = None;

        // `symmetries.orbitopes` and `symmetries.computeStabilizerOrbits(getDomain(), ..)`
        // act on the moved-out `symmetries` (with_symmetries); the domain is detached for
        // the call, and `compute_stabilizer_orbits` only reads the model through `ms`.
        self.with_symmetries(|ms, sym| {
            for orbitope in sym.orbitopes.iter_mut() {
                orbitope.determine_orbitope_type(&mut ms.mipdata_.cliquetable);
            }

            // mip_symmetry_dynamic: SCIP sym_rowcol.c lines 1040-1180 at 32dd425
            if ms.options_mip_.mip_symmetry_dynamic {
                // R7/Q12: sub-MIPs run with mip_detect_symmetry=false
                // (primal_heuristics/part1.rs:244), so data is never built there.
                debug_assert!(!ms.submip);
                // SCIPvarIsIntegral: integer or implied-integer columns
                let is_integral: Vec<bool> = mip_model!(ms)
                    .integrality_
                    .iter()
                    .map(|t| matches!(*t, HighsVarType::kInteger | HighsVarType::kImplicitInteger))
                    .collect();
                ms.mipdata_.symmetry_dynamic =
                    crate::mip::symmetry_dynamic::setup::s3_build(&*sym, &is_integral);
            }

            if sym.num_perms != 0 {
                let mut workspace = StabilizerOrbitWorkspace::default();
                let orbits = detach!(ms.mipdata_.domains[0], |d| sym.compute_stabilizer_orbits(
                    ms,
                    d,
                    &mut workspace
                ));
                ms.mipdata_.global_orbits = orbits;
            }
        });
    }

    /// rusolve: take the finished symmetry data out of `mipdata_` (see `restore_symmetries`).
    pub(crate) fn stash_symmetries(&mut self) -> SymmetryStash {
        SymmetryStash {
            symmetries: std::mem::take(&mut self.mipdata_.symmetries),
            global_orbits: self.mipdata_.global_orbits.take(),
            symmetry_dynamic: self.mipdata_.symmetry_dynamic.take(),
        }
    }

    /// rusolve: put stashed symmetry data back and recompute the stabilizer orbits on the
    /// current root domain, as `finish_symmetry_detection` would have at this point.
    pub(crate) fn restore_symmetries(&mut self, stash: SymmetryStash) {
        self.mipdata_.symmetries = stash.symmetries;
        self.mipdata_.symmetry_dynamic = stash.symmetry_dynamic;
        self.with_symmetries(|ms, sym| {
            if sym.num_perms != 0 {
                let mut workspace = StabilizerOrbitWorkspace::default();
                let orbits = detach!(ms.mipdata_.domains[0], |d| sym.compute_stabilizer_orbits(
                    ms,
                    d,
                    &mut workspace
                ));
                ms.mipdata_.global_orbits = orbits;
            }
        });
    }

    /// rusolve `mip_orbitope_full_fix`: lexicographic orbitopal fixing
    /// (Kaibel, Peinhardt, Pfetsch 2011; SCIP cons_orbitope.c, full case) over all rows of each
    /// orbitope that is not entirely set-packing. It is the same constraint the tree applies to
    /// branched rows only, so the dynamic handling stays valid and nothing is detached.
    pub fn orbitope_full_fix(&mut self) -> HighsInt {
        let mut fixed: HighsInt = 0;
        self.with_symmetries(|ms, sym| {
            for orb in sym.orbitopes.iter() {
                let all_packing = (0..orb.num_rows).all(|r| {
                    orb.row_is_set_packing[r as usize] > HighsOrbitopeMatrix::kRowNotPacking
                });
                if all_packing {
                    continue;
                }
                let rows: Vec<HighsInt> = (0..orb.num_rows).collect();
                fixed += detach!(ms.mipdata_.domains[0], |d| orb
                    .orbital_fixing_for_full_orbitope(ms, &rows, d));
                if ms.mipdata_.domains[0].infeasible() {
                    break;
                }
            }
        });
        if fixed > 0 {
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "Orbitopal fixing: {} columns fixed\n\n",
                fixed
            );
        }
        fixed
    }

    /// rusolve `mip_orbitope_root_fix`: static root fixing from the
    /// packing orbitopes found by symmetry detection, before the first root LP.
    ///
    /// Source: SCIP 10 (Apache-2.0) `cons_orbitope.c` packing/partitioning orbitope propagation
    /// (entries right of the staircase are zero) and the `symmetry` presolver that adds the
    /// orbitopes (Hojny, Pfetsch, "Polytopes associated with symmetry handling", Math. Program.
    /// 175, 2019). The fixing itself reuses the HiGHS frontier routine
    /// (`orbital_fixing_for_packing_orbitope`) with every packing row of the orbitope in matrix
    /// order, which gives that staircase. Dynamic orbital fixing in the tree orders by branching,
    /// so both cannot be combined: the detected symmetries are dropped afterwards.
    /// Returns the number of fixed columns.
    pub fn orbitope_root_fix(&mut self) -> HighsInt {
        let mut fixed: HighsInt = 0;
        let mut keep_dynamic = false;
        let num_col = self.num_col();
        self.with_symmetries(|ms, sym| {
            // Gate: the staircase of a packing orbitope with R packing rows and q columns fixes
            // sum_{i < min(R, q)} (q - 1 - i) entries. Below 2 % of the columns the fixing does
            // not repay giving up the branching-order symmetry handling (measured on comp21-2idx,
            // supportcase18, seymour).
            // TODO(perf): 2 % is a measured split, not a tuned value.
            let mut potential: i64 = 0;
            for orb in sym.orbitopes.iter() {
                let rows = (0..orb.num_rows)
                    .filter(|&r| {
                        orb.row_is_set_packing[r as usize] > HighsOrbitopeMatrix::kRowNotPacking
                    })
                    .count() as i64;
                let q = orb.row_length as i64;
                for i in 0..rows.min(q) {
                    potential += q - 1 - i;
                }
            }
            if (potential as f64) < 0.02 * num_col as f64 {
                return;
            }
            let mut handled = vec![false; sym.orbitopes.len()];
            for (k, orb) in sym.orbitopes.iter().enumerate() {
                let rows: Vec<HighsInt> = (0..orb.num_rows)
                    .filter(|&r| {
                        orb.row_is_set_packing[r as usize] > HighsOrbitopeMatrix::kRowNotPacking
                    })
                    .collect();
                if rows.is_empty() {
                    continue;
                }
                handled[k] = true;
                fixed += detach!(ms.mipdata_.domains[0], |d| orb
                    .orbital_fixing_for_packing_orbitope(ms, &rows, d));
                if ms.mipdata_.domains[0].infeasible() {
                    break;
                }
            }
            if fixed == 0 {
                return;
            }
            // Dynamic symmetry handling stays valid for the columns the static fixing did not
            // touch: when no remaining permutation moves a column of a handled orbitope, the
            // two lexicographic schemes act on disjoint column sets. Otherwise drop everything.
            let mut handled_cols = std::collections::HashSet::<HighsInt>::new();
            for (k, orb) in sym.orbitopes.iter().enumerate() {
                if handled[k] {
                    handled_cols.extend(orb.matrix.iter().copied());
                }
            }
            let disjoint = !sym
                .permutation_columns
                .iter()
                .any(|pc| handled_cols.contains(pc));
            if disjoint {
                let mut k = 0;
                sym.orbitopes.retain(|_| {
                    k += 1;
                    !handled[k - 1]
                });
                sym.column_to_orbitope.clear();
                for (idx, orb) in sym.orbitopes.iter().enumerate() {
                    for &col in orb.matrix.iter() {
                        sym.column_to_orbitope
                            .insert(HighsHashTableEntry::new(col, idx as HighsInt));
                    }
                }
                keep_dynamic = sym.num_perms != 0 || !sym.orbitopes.is_empty();
                if sym.num_perms != 0 {
                    let mut workspace = StabilizerOrbitWorkspace::default();
                    let orbits = detach!(ms.mipdata_.domains[0], |d| sym
                        .compute_stabilizer_orbits(ms, d, &mut workspace));
                    ms.mipdata_.global_orbits = orbits;
                } else {
                    ms.mipdata_.global_orbits = None;
                }
            } else {
                *sym = HighsSymmetries::default();
                ms.mipdata_.global_orbits = None;
            }
        });
        if fixed > 0 {
            if !keep_dynamic {
                self.mipdata_.detect_symmetries = false;
            }
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "Orbitope root fixing: {} columns fixed{}\n\n",
                fixed,
                if keep_dynamic {
                    ", other symmetries kept"
                } else {
                    ""
                }
            );
        }
        fixed
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::limitsToGap (lines 608-626)
    /// C++: `double limitsToGap(const double use_lower_bound, const double use_upper_bound, double & lb, double & ub) const`
    pub fn limits_to_gap(
        &self,
        use_lower_bound: f64,
        use_upper_bound: f64,
        lb: &mut f64,
        ub: &mut f64,
    ) -> f64 {
        let offset = mip_model!(self).offset_;
        *lb = use_lower_bound + offset;
        if lb.abs() <= self.mipdata_.epsilon {
            *lb = 0.0;
        }
        *ub = kHighsInf;
        let mut gap = kHighsInf;
        if use_upper_bound != kHighsInf {
            *ub = use_upper_bound + offset;
            if ub.abs() <= self.mipdata_.epsilon {
                *ub = 0.0;
            }
            *lb = cpp_min(*ub, *lb);
            if *ub == 0.0 {
                gap = if *lb == 0.0 { 0.0 } else { kHighsInf };
            } else {
                gap = (*ub - *lb) / ub.abs();
            }
        }
        gap
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::computeNewUpperLimit (lines 628-667)
    /// C++: `double computeNewUpperLimit(double ub, double mip_abs_gap, double mip_rel_gap) const`
    pub fn compute_new_upper_limit(&self, ub: f64, mip_abs_gap: f64, mip_rel_gap: f64) -> f64 {
        let offset = mip_model!(self).offset_;
        let epsilon = self.mipdata_.epsilon;
        let feastol = self.mipdata_.feastol;
        let mut new_upper_limit: f64;
        if self.mipdata_.objective_function.is_integral() {
            let integral_scale = self.mipdata_.objective_function.integral_scale();
            new_upper_limit = (integral_scale * ub - 0.5).floor() / integral_scale;

            if mip_rel_gap != 0.0 {
                new_upper_limit = cpp_min(
                    new_upper_limit,
                    ub - (mip_rel_gap * (ub + offset).abs() * integral_scale - epsilon).ceil()
                        / integral_scale,
                );
            }

            if mip_abs_gap != 0.0 {
                new_upper_limit = cpp_min(
                    new_upper_limit,
                    ub - (mip_abs_gap * integral_scale - epsilon).ceil() / integral_scale,
                );
            }

            // add feasibility tolerance so that the next best integer feasible solution
            // is definitely included in the remaining search
            new_upper_limit += feastol;
        } else {
            new_upper_limit = cpp_min(ub - feastol, next_down(ub));

            if mip_rel_gap != 0.0 {
                new_upper_limit = cpp_min(new_upper_limit, ub - mip_rel_gap * (ub + offset).abs());
            }

            if mip_abs_gap != 0.0 {
                new_upper_limit = cpp_min(new_upper_limit, ub - mip_abs_gap);
            }
        }

        new_upper_limit
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::moreHeuristicsAllowed (lines 669-731)
    /// C++: `bool moreHeuristicsAllowed() const`
    pub fn more_heuristics_allowed(&self) -> bool {
        let md = &*self.mipdata_;
        // in the beginning of the search and in sub-MIP heuristics we only allow
        // what is proportionally for the currently spent effort plus an initial
        // offset. This is because in a sub-MIP we usually do a truncated search and
        // therefore should not extrapolate the time we spent for heuristics as in
        // the other case. Moreover, since we estimate the total effort for
        // exploring the tree based on the weight of the already pruned nodes, the
        // estimated effort the is not expected to be a good prediction in the
        // beginning.
        if self.submip {
            (md.heuristic_lp_iterations as f64)
                < (md.total_lp_iterations as f64) * md.heuristic_effort
        } else if md.pruned_treeweight < 1e-3
            && md.num_leaves - md.num_leaves_before_run < 10
            && md.num_nodes - md.num_nodes_before_run < 1000
        {
            // in the main MIP solver allow an initial offset of 10000 heuristic LP
            // iterations
            (md.heuristic_lp_iterations as f64)
                < (md.total_lp_iterations as f64) * md.heuristic_effort + 10000.0
        } else if md.heuristic_lp_iterations
            < 100000
                + ((md.total_lp_iterations - md.heuristic_lp_iterations - md.sb_lp_iterations) >> 1)
        {
            // compute the node LP iterations in the current run as only those should be
            // used when estimating the total required LP iterations to complete the
            // search
            let heur_iters_curr_run =
                md.heuristic_lp_iterations - md.heuristic_lp_iterations_before_run;
            let sb_iters_curr_run = md.sb_lp_iterations - md.sb_lp_iterations_before_run;
            let node_iters_curr_run = md.total_lp_iterations
                - md.total_lp_iterations_before_run
                - heur_iters_curr_run
                - sb_iters_curr_run;
            // now estimate the total fraction of LP iterations that we have spent on
            // heuristics by assuming the node iterations of the current run will
            // grow proportional to the pruned weight of the current tree and the
            // iterations spent for anything else are just added as an offset
            let total_heuristic_effort_estim = (md.heuristic_lp_iterations as f64)
                / (((md.total_lp_iterations - node_iters_curr_run) as f64)
                    + (node_iters_curr_run as f64)
                        / cpp_max(0.01, f64::from(md.pruned_treeweight)));
            // since heuristics help most in the beginning of the search, we want to
            // spent the time we have for heuristics in the first 80% of the tree
            // exploration. Additionally we want to spent the proportional effort
            // of heuristics that is allowed in the first 30% of tree exploration as
            // fast as possible, which is why we have the max(0.3/0.8,...).
            // Hence, in the first 30% of the tree exploration we allow to spent all
            // effort available for heuristics in that part of the search as early as
            // possible, whereas after that we allow the part that is proportionally
            // adequate when we want to spent all available time in the first 80%.
            total_heuristic_effort_estim
                < cpp_max(
                    0.3 / 0.8,
                    cpp_min(f64::from(md.pruned_treeweight), 0.8) / 0.8,
                ) * md.heuristic_effort
        } else {
            false
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::removeFixedIndices (lines 733-750)
    /// C++: `void removeFixedIndices()`
    pub fn remove_fixed_indices(&mut self) {
        // `getDomain().isFixed(col)` reads the attached global domain; the four
        // `std::remove_if` + `erase` pairs are `retain(!p)` (order kept).
        let domain = &self.mipdata_.domains[0];
        self.mipdata_
            .integral_cols
            .retain(|&col| !domain.is_fixed(col));
        self.mipdata_
            .integer_cols
            .retain(|&col| !domain.is_fixed(col));
        self.mipdata_
            .implint_cols
            .retain(|&col| !domain.is_fixed(col));
        self.mipdata_
            .continuous_cols
            .retain(|&col| !domain.is_fixed(col));
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::init (lines 752-809)
    /// C++: `void init()`
    pub fn init(&mut self) {
        let num_row = self.num_row();
        let num_col = self.num_col();
        self.mipdata_
            .post_solve_stack
            .initialize_index_maps(num_row, num_col);
        // PORT NOTE: C++ `orig_model_ = model_` points `orig_model_` at the current model.
        // The port's `orig_model_` is owned and `model_` selects it, so the assignment is
        // the identity while `model_` is `Orig`. The call site (`run_mip_presolve`) runs
        // before any `set_model(Presolved)`; checked by the debug_assert below.
        debug_assert!(self.model_ == HighsMipModel::Orig);
        self.mipdata_.feastol = self.options_mip_.mip_feasibility_tolerance;
        self.mipdata_.epsilon = self.options_mip_.small_matrix_value;
        // `cliquetable.buildFrom(mipsolver.orig_model_, *clqtableinit)`: read in place
        // through the parent pointer.
        if let Some(clqtableinit) = self.clqtableinit {
            self.mipdata_
                .cliquetable
                .build_from(Some(&self.orig_model_), &clqtableinit);
        }
        // C++ `highs::parallel::num_threads() > 1 ? option : kHighsIInf`.
        let min_entries_for_parallelism = if crate::parallel::highs_parallel::num_threads() > 1 {
            self.options_mip_
                .mip_min_cliquetable_entries_for_parallelism
        } else {
            kHighsIInf
        };
        self.mipdata_
            .cliquetable
            .set_min_entries_for_parallelism(min_entries_for_parallelism);
        if let Some(implicinit) = self.implicinit {
            detach!(self.mipdata_.implications, |im| im
                .build_from(self, &implicinit));
        }
        self.mipdata_.heuristic_effort = self.options_mip_.mip_heuristic_effort;
        self.mipdata_.detect_symmetries = self.options_mip_.mip_detect_symmetry;

        self.mipdata_.firstlpsolobj = -kHighsInf;
        self.mipdata_.rootlpsolobj = -kHighsInf;
        self.mipdata_.analytic_center_computed = false;
        self.mipdata_.analytic_center_status = HighsModelStatus::kNotset;
        self.mipdata_.max_tree_size_log2 = 0;
        self.mipdata_.num_restarts = 0;
        self.mipdata_.num_restarts_root = 0;
        self.mipdata_.num_improving_sols = 0;
        self.mipdata_.pruned_treeweight = HighsCDouble::from(0.0);
        self.mipdata_.avgrootlpiters = 0.0;
        self.mipdata_.num_nodes = 0;
        self.mipdata_.num_nodes_before_run = 0;
        self.mipdata_.num_leaves = 0;
        self.mipdata_.num_leaves_before_run = 0;
        self.mipdata_.total_repair_lp = 0;
        self.mipdata_.total_repair_lp_feasible = 0;
        self.mipdata_.total_repair_lp_iterations = 0;
        self.mipdata_.total_lp_iterations = 0;
        self.mipdata_.heuristic_lp_iterations = 0;
        self.mipdata_.sepa_lp_iterations = 0;
        self.mipdata_.sb_lp_iterations = 0;
        self.mipdata_.total_lp_iterations_before_run = 0;
        self.mipdata_.heuristic_lp_iterations_before_run = 0;
        self.mipdata_.sepa_lp_iterations_before_run = 0;
        self.mipdata_.sb_lp_iterations_before_run = 0;
        self.mipdata_.num_disp_lines = 0;
        self.mipdata_.num_clique_entries_after_presolve = 0;
        self.mipdata_.num_clique_entries_after_first_presolve = 0;
        self.mipdata_.cliques_extracted = false;
        self.mipdata_.row_matrix_set = false;
        self.mipdata_.lower_bound = -kHighsInf;
        self.mipdata_.upper_bound = kHighsInf;
        self.mipdata_.upper_limit = self.options_mip_.objective_bound;
        self.mipdata_.optimality_limit = self.options_mip_.objective_bound;
        // HiGHS 1.16 #3124 #3521 (h16_mip_workers): HighsMipSolverData.cpp lines 806-809 (1.16.0)
        let h16_mip_workers = self.options_mip_.h16_mip_workers;
        if h16_mip_workers {
            self.mipdata_
                .worker_lp_iterations_stop
                .0
                .store(i64::MAX, std::sync::atomic::Ordering::Relaxed);
        }
        self.mipdata_.primal_dual_integral.initialise();
        if h16_mip_workers {
            // The P-D integral is measured from the start of the MIP solver
            self.mipdata_.primal_dual_integral.prev_time = self.timer_.read();
        }

        let report_level = self.options_mip_.mip_report_level;
        self.mipdata_.dispfreq = if report_level == 0 {
            0
        } else if report_level == 1 {
            2000
        } else {
            100
        };
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::runMipPresolve (lines 811-829)
    /// C++: `void runMipPresolve(const HighsInt presolve_reduction_limit)`
    ///
    /// PORT NOTE (architecture): ownership protocol of `presolve.okSetInput(mipsolver, ..)` /
    /// `presolve.run(postSolveStack)`: the solver
    /// is moved into the `HPresolve` (`Box`, `mem::take(self)`) for the duration and moved
    /// back afterwards; `postSolveStack` is taken out of `mipdata_` first and put back.
    pub fn run_mip_presolve_data(&mut self, presolve_reduction_limit: HighsInt) {
        // PORT NOTE: `timer_.start(presolve_clock)` / `stop(..)` are omitted: the stub
        // HighsTimer has one clock and `start()` restarts it (util/highs_timer.rs), which
        // would reset the MIP run clock that `time_limit` reads.
        let mut post_solve_stack = std::mem::take(&mut self.mipdata_.post_solve_stack);
        let mut presolve = HPresolve::default();
        let ok =
            presolve.ok_set_input_mip(Box::new(std::mem::take(self)), presolve_reduction_limit);
        let (model_status, presolve_status) = if !ok {
            (
                HighsModelStatus::kMemoryLimit,
                HighsPresolveStatus::kOutOfMemory,
            )
        } else {
            let st = presolve.run(&mut post_solve_stack);
            (st, presolve.get_presolve_status())
        };
        *self = *presolve.take_mip_solver();
        self.mipdata_.post_solve_stack = post_solve_stack;
        self.modelstatus_ = model_status;
        self.mipdata_.presolve_status = presolve_status;
        // Report the final presolve reductions unless this is a restart
        if self.options_mip_.presolve != kHighsOffString && self.mipdata_.num_restarts == 0 {
            crate::lp_data::presolve_support::report_presolve_reductions(
                &self.options_mip_.log_options,
                presolve_status,
                &self.orig_model_,
                mip_model!(self),
            );
        }
    }
}

/// C++ `std::nextafter(x, -kHighsInf)`: the next double below `x` (IEEE bit step). The
/// crate has no `nextafter`; this follows the C99 rules: NaN and -inf are returned as is,
/// and both zeros step to `-denorm_min`.
fn next_down(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return x;
    }
    if x == 0.0 {
        return -f64::from_bits(1);
    }
    if x > 0.0 {
        f64::from_bits(x.to_bits() - 1)
    } else {
        f64::from_bits(x.to_bits() + 1)
    }
}

/// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::startAnalyticCenterComputation (lines 411-486)
///
/// The `[&]` lambda that `startAnalyticCenterComputation` spawns. It captures `this`; the
/// port holds field-level pointers to the members the body touches (rule (s), manifest
/// open item 8), taken before the spawn: it reads `mipsolver.options_mip_` and
/// `mipsolver.model_` (`*model_` and `numCol()`), and writes `analyticCenterStatus` and
/// `analyticCenter`.
struct AnalyticCenterTask {
    /// `mipsolver.options_mip_`
    options_mip_: *const HighsOptions,
    /// `mipsolver.model_` (the selector; read when the task runs, as C++ reads the pointer)
    model_: *const HighsMipModel,
    /// `*mipsolver.orig_model_` (what `model_` selects as `Orig`)
    orig_model_: *const HighsLp,
    /// `mipdata_->presolvedModel` (what `model_` selects as `Presolved`)
    presolved_model: *const HighsLp,
    /// Port-only (the oracle `ac` cross-check): `mipsolver.submip_level`, constant for the
    /// life of the solver.
    submip_level: HighsInt,
    /// HiGHS 1.16 #3456 / #3522 (h16_ac_centring): `mipsolver.options_mip_->h16_ac_centring`. Off:
    /// the 1.15.1 HiPO choice and no centring-status check.
    h16_ac_centring: bool,
    /// `analyticCenterStatus`
    analytic_center_status: *mut HighsModelStatus,
    /// `analyticCenter`
    analytic_center: *mut Vec<f64>,
}

// SAFETY: the C++ lambda runs on whichever worker takes it (the owner at `sync`, or a
// thief). The pointees outlive the task: they are members of the `HighsMipSolver` that
// `evaluateRootNode` borrows for the whole life of its `TaskGroup`, whose destructor
// (`cancel` + `taskWait`) runs before that borrow ends (rule (g)). Access protocol: see
// `AnalyticCenterTask::run`.
unsafe impl Send for AnalyticCenterTask {}

impl AnalyticCenterTask {
    /// The lambda body (cpp 411-486).
    ///
    /// PORT NOTE: `ipm.setProfiling` and the profiling clocks are not ported; the
    /// `dump_ipm_lp` / `ipm_logging` blocks are dead (`const bool ... = false`). The HiPO
    /// retry (cpp 470-479) is dead too: `use_hipo` implies `mip_ipm_solver == "hipo"`, so
    /// `mip_ipm_solver == kHighsChooseString` never holds with it.
    ///
    /// # Safety
    /// C++ protocol of `evaluateRootNode`: while this task is pending or running, the main
    /// worker does not write `options_mip_`, `model_`, `orig_model_` or
    /// `mipdata_.presolvedModel` (they change only in `performRestart`, after
    /// `tg.cancel(); tg.taskWait();`), and it reads `analyticCenterStatus` /
    /// `analyticCenter` only after `taskGroup.sync()` (`finishAnalyticCenterComputation`)
    /// or `taskWait()` / `~TaskGroup`.
    unsafe fn run(self) {
        // SAFETY: shared reads of members the main worker does not write while the task
        // is pending (see the function doc).
        let options_mip_ = unsafe { &*self.options_mip_ };
        // SAFETY: as above.
        let model_ = unsafe {
            match *self.model_ {
                HighsMipModel::Orig => &*self.orig_model_,
                HighsMipModel::Presolved => &*self.presolved_model,
            }
        };
        // Highs instantiation
        let mut ipm = Highs::default();
        ipm.set_option_value_bool("output_flag", false);
        // h16_sub_solver_flags: the C++ `Highs ipm` runs 1.16.0 code
        if options_mip_.h16_sub_solver_flags {
            ipm.options_.copy_version_flags_from(options_mip_);
        }
        // Don't use presolve - because this can lead to postsolve putting
        // integer variables onto bounds. This is not just a "less good"
        // AC. It can have implications leading to erroneous fixing of
        // variables and a suboptimal solution declared as optimal.
        ipm.set_option_value_str("presolve", kHighsOffString);
        // Determine the solver
        let mip_ipm_solver = options_mip_.mip_ipm_solver.clone();
        // Currently use IPX by default and take action on failure here if
        // using HiPO.
        // HiGHS 1.16 #3522 (h16_ac_centring): HighsMipSolverData.cpp::startAnalyticCenterComputation
        // lines 431-432 (1.16.0), `use_hipo = false` ("HiPO does not have centring yet").
        let use_hipo = h16_ac_use_hipo(self.h16_ac_centring, &mip_ipm_solver);
        // Later still, pass mip_ipm_solver and take action on failure in
        // solveLp
        let ipm_solver = if use_hipo { kHipoString } else { kIpxString };
        ipm.set_option_value_str("solver", ipm_solver);
        ipm.set_option_value_int("ipm_iteration_limit", 200);
        ipm.set_option_value_str("run_crossover", kHighsOffString);
        ipm.set_option_value_bool("run_centring", true);
        // HiGHS 1.16 #3456 (h16_ac_centring): the IPX centring record (ipx_wrapper.rs
        // `h16_centring_status`) reads `options.h16_ac_centring` of the IPM's own options,
        // which C++ has at its default; the flag is carried to the IPM here.
        if self.h16_ac_centring {
            ipm.set_option_value_bool("h16_ac_centring", true);
        }
        let mut lpmodel: HighsLp = model_.clone();
        lpmodel.col_cost_ = vec![0.0; lpmodel.num_col_ as usize];
        lpmodel.integrality_.clear();
        ipm.pass_model(lpmodel);
        ipm.optimize_lp();
        // const std::vector<double>& sol = ipm.getSolution().col_value;
        let sol = &ipm.get_solution().col_value;
        // `mipsolver.numCol()`
        let num_col = model_.num_col_;
        // Port-only: compare with the oracle `ac` record when RUSOLVE_MIP_AC_FIXTURE is set
        // (this thread's replay state: manifest open items 7 and 10).
        ipm_replay_check(
            kMipIpmRecordAc,
            self.submip_level,
            sol.len() as HighsInt == num_col,
            ipm.get_model_status(),
            ipm.get_solution(),
            ipm.get_basis(),
        );
        if sol.len() as HighsInt != num_col {
            return;
        }
        // SAFETY: the main worker does not touch these two members until `sync` /
        // `taskWait` has returned, i.e. after this write (see the function doc).
        unsafe {
            *self.analytic_center_status = ipm.get_model_status();
            // HiGHS 1.16 #3456 (h16_ac_centring): HighsMipSolverData.cpp::startAnalyticCenterComputation
            // lines 479-484 (1.16.0): the model status alone does not show that the point is
            // centred, so an optimal status without a successful centring becomes kUnknown.
            if self.h16_ac_centring && h16_ac_centring_rejects(ipm.get_info().centring_status) {
                *self.analytic_center_status = HighsModelStatus::kUnknown;
            }
            *self.analytic_center = sol.clone();
        }
    }
}

/// HiGHS 1.16 #3522 (h16_ac_centring): the `use_hipo` of HighsMipSolverData.cpp
/// startAnalyticCenterComputation (1.16.0, lines 431-432) is `false`; off, the 1.15.1
/// expression `mip_ipm_solver == kHipoString`.
fn h16_ac_use_hipo(h16_ac_centring: bool, mip_ipm_solver: &str) -> bool {
    if h16_ac_centring {
        false
    } else {
        mip_ipm_solver == kHipoString
    }
}

/// HiGHS 1.16 #3456 (h16_ac_centring): HighsMipSolverData.cpp::startAnalyticCenterComputation lines
/// 483-484 (1.16.0): `ipm.getInfo().centring_status != kCentringStatusSuccess`.
/// Only called from the `h16_ac_centring` branch of `AnalyticCenterTask::run`.
fn h16_ac_centring_rejects(centring_status: HighsInt) -> bool {
    centring_status != kCentringStatusSuccess
}

/// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::startSymmetryDetection (lines 553-557)
///
/// The `[&]` lambda that `startSymmetryDetection` spawns. It captures `this` (for
/// `mipsolver.timer_`; `symDetection` reads `presolvedModel`, which `loadModelAsGraph`
/// bound it to) and the caller's `symData` by reference. `symData` is not reseated while
/// the task is pending (`startSymmetryDetection` replaces it only after a restart's
/// `taskWait`, `finishSymmetryDetection` resets it after `sync`), so the port points at the
/// boxed `SymmetryDetectionData` itself.
struct SymmetryDetectionTask {
    /// `mipsolver.timer_`
    timer_: *const HighsTimer,
    /// `mipsolver.options_mip_->time_limit` (HighsMipSolverData.cpp lines 559-560 at v1.16.0)
    time_limit: f64,
    /// `mipdata_->presolvedModel` (the model `symDetection` was loaded from)
    presolved_model: *const HighsLp,
    /// `*symData`
    sym_data: *mut HighsMipSolverData_SymmetryDetectionData,
}

// SAFETY: as `AnalyticCenterTask`: the pointees outlive the task (the scoped `TaskGroup`
// of `evaluateRootNode`, rule (g)); access protocol in `SymmetryDetectionTask::run`.
unsafe impl Send for SymmetryDetectionTask {}

impl SymmetryDetectionTask {
    /// The lambda body (cpp 554-556).
    ///
    /// # Safety
    /// C++ protocol of `evaluateRootNode`: the main worker only reads `timer_`, does not
    /// write `mipdata_.presolvedModel` while the task is pending (only `performRestart`
    /// does, after `cancel` + `taskWait`), and touches `*symData` only after
    /// `taskGroup.sync()` (`finishSymmetryDetection`) or `taskWait()` / `~TaskGroup`.
    unsafe fn run(self) {
        // SAFETY: shared reads (see the function doc).
        let (timer_, presolved_model) = unsafe { (&*self.timer_, &*self.presolved_model) };
        // SAFETY: `*symData` is this task's alone until `sync` / `taskWait` returns.
        let sym_data = unsafe { &mut *self.sym_data };
        let start_time = timer_.get_wall_time();
        // HiGHS #3535 (mip_sym_hash64): `run(symmetries, &timer_, time_limit)` in 1.16.0
        let time_limit = if sym_data.sym_detection.mip_sym_hash64 {
            Some((timer_, self.time_limit))
        } else {
            None
        };
        sym_data
            .sym_detection
            .run(presolved_model, &mut sym_data.symmetries, time_limit);
        sym_data.detection_time = timer_.get_wall_time() - start_time;
    }
}

#[cfg(test)]
mod h16_ac_fix_tests {
    //! HiGHS 1.16-ACFIX (#3408): the analytic-centre fixing is off without presolve, and uses
    //! feastol / 10 as its tolerance when `h16_ac_fix` is on.

    use super::*;
    use crate::parallel::highs_parallel::{initialize_scheduler, task_group};
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    /// Three binary columns, one row `x0 + x1 + x2 <= 2`.
    fn three_binary_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![-1.0; 3];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![1.0; 3];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![2.0];
        lp.integrality_ = vec![HighsVarType::kInteger; 3];
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2, 3];
        lp.a_matrix_.index_ = vec![0, 0, 0];
        lp.a_matrix_.value_ = vec![1.0; 3];
        lp
    }

    /// Run `finish_analytic_center_computation` with column 0 half a feastol above its lower
    /// bound (feastol = 1e-6) and the other columns central. Returns the upper and lower
    /// bound of column 0 after the call.
    fn finish_ac_probe(h16_ac_fix: bool, presolve: &'static str) -> (f64, f64) {
        finish_ac_probe_at(h16_ac_fix, presolve, 0.5e-6)
    }

    /// As `finish_ac_probe`, with column 0 of the analytic centre at `x0`.
    fn finish_ac_probe_at(h16_ac_fix: bool, presolve: &'static str, x0: f64) -> (f64, f64) {
        std::thread::spawn(move || {
            initialize_scheduler(1);
            let mut ms = HighsMipSolver::default();
            ms.orig_model_ = three_binary_lp();
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.presolve = presolve.to_string();
            options.h16_ac_fix = h16_ac_fix;
            ms.options_mip_ = options;
            ms.new_data();
            // row-major copy, domain and row activities, as before propagation in the run
            ms.setup_domain_propagation();
            ms.mipdata_.feastol = 1e-6;
            ms.mipdata_.analytic_center_status = HighsModelStatus::kOptimal;
            ms.mipdata_.analytic_center = vec![x0, 0.5, 0.5];
            task_group(|tg| {
                tg.spawn(|| {});
                ms.finish_analytic_center_computation(tg);
            });
            let out = {
                let dom = &ms.mipdata_.domains[0];
                (dom.col_upper_[0], dom.col_lower_[0])
            };
            drop(ms);
            HighsTaskExecutor::shutdown(true);
            out
        })
        .join()
        .unwrap()
    }

    #[test]
    fn h16_ac_fix_on_keeps_column_within_feastol() {
        // 0.5e-6 is above feastol / 10 = 1e-7, so the column is not fixed to its lower bound.
        let (upper, lower) = finish_ac_probe(true, kHighsOnString);
        assert_eq!(upper, 1.0);
        assert_eq!(lower, 0.0);
    }

    #[test]
    fn h16_ac_fix_on_fixes_column_within_tenth_feastol() {
        // HighsMipSolverData.cpp:517-518 (1.16.0): tolerance = (1e-6 / 10) * min(1, 1) = 1e-7,
        // and 0.5e-7 <= 0 + 1e-7, so the upper bound is changed to the lower bound 0.
        let (upper, lower) = finish_ac_probe_at(true, kHighsOnString, 0.5e-7);
        assert_eq!(upper, 0.0);
        assert_eq!(lower, 0.0);
    }

    #[test]
    fn h16_ac_fix_on_presolve_off_skips_fixing() {
        // HighsMipSolverData.cpp:508-509 (1.16.0): no fixing when presolve is off,
        // even for a column inside the 1.16 tolerance.
        let (upper, lower) = finish_ac_probe_at(true, kHighsOffString, 0.5e-7);
        assert_eq!(upper, 1.0);
        assert_eq!(lower, 0.0);
    }

    #[test]
    fn h16_ac_fix_on_presolve_off_skips_fixing_feastol() {
        // HighsMipSolverData.cpp:508-509 (1.16.0): no fixing when presolve is off.
        let (upper, lower) = finish_ac_probe(true, kHighsOffString);
        assert_eq!(upper, 1.0);
        assert_eq!(lower, 0.0);
    }

    #[test]
    fn h16_ac_fix_off_fixes_column_within_feastol() {
        // 1.15.1: tolerance is feastol, so 0.5e-6 <= 0 + 1e-6 fixes the column at its lower bound.
        let (upper, lower) = finish_ac_probe(false, kHighsOnString);
        assert_eq!(upper, 0.0);
        assert_eq!(lower, 0.0);
    }

    #[test]
    fn h16_ac_fix_off_presolve_off_still_fixes() {
        // 1.15.1 has no presolve check, so the fixing still runs with presolve off.
        let (upper, lower) = finish_ac_probe(false, kHighsOffString);
        assert_eq!(upper, 0.0);
        assert_eq!(lower, 0.0);
    }
}

#[cfg(test)]
mod h16_ac_centring_tests {
    //! HiGHS 1.16-ACCENTRE (#3456, #3522): the analytic-centre status and HiPO choice are changed
    //! only when `h16_ac_centring` is on; off, the 1.15.1 values hold.

    use super::*;
    use crate::ipm::ipx_wrapper::h16_centring_status;

    /// Two continuous columns in [0, 1], one row `x0 + x1 <= 1.5`, cost `-x0 - x1`.
    fn two_col_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![-1.0; 2];
        lp.col_lower_ = vec![0.0; 2];
        lp.col_upper_ = vec![1.0; 2];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![1.5];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0; 2];
        lp
    }

    #[test]
    fn h16_ac_centring_on_forces_use_hipo_false() {
        // HighsMipSolverData.cpp:431-432 (1.16.0): use_hipo = false even when mip_ipm_solver is hipo.
        assert!(!h16_ac_use_hipo(true, kHipoString));
        assert!(!h16_ac_use_hipo(true, kIpxString));
    }

    #[test]
    fn h16_ac_centring_off_keeps_use_hipo() {
        // 1.15.1: use_hipo = mip_ipm_solver == kHipoString.
        assert!(h16_ac_use_hipo(false, kHipoString));
        assert!(!h16_ac_use_hipo(false, kIpxString));
    }

    #[test]
    fn h16_ac_centring_on_failed_or_untried_centring_is_rejected() {
        // HighsMipSolverData.cpp:483-484 (1.16.0): centring_tried & failed -> kUnknown, not tried -> kUnknown.
        assert!(h16_ac_centring_rejects(kCentringStatusFailed));
        assert!(h16_ac_centring_rejects(kCentringStatusNone));
    }

    #[test]
    fn h16_ac_centring_on_successful_centring_is_accepted() {
        // Success keeps the model status unchanged.
        assert!(!h16_ac_centring_rejects(kCentringStatusSuccess));
    }

    #[test]
    fn h16_ac_centring_status_mapping() {
        // IpxWrapper.cpp:205-208 (1.16.0): not tried -> None; tried and success -> Success; tried and failed -> Failed.
        assert_eq!(h16_centring_status(0, 0), kCentringStatusNone);
        assert_eq!(h16_centring_status(0, 1), kCentringStatusNone);
        assert_eq!(h16_centring_status(1, 0), kCentringStatusFailed);
        assert_eq!(h16_centring_status(1, 1), kCentringStatusSuccess);
    }

    /// Solve `two_col_lp` with IPX and `run_centring` on, as the analytic-centre task does,
    /// and return the centring status recorded in the IPM's info.
    fn ipx_centring_status(h16_ac_centring: bool) -> HighsInt {
        let mut ipm = Highs::default();
        ipm.set_option_value_bool("output_flag", false);
        ipm.set_option_value_str("presolve", kHighsOffString);
        ipm.set_option_value_str("solver", kIpxString);
        ipm.set_option_value_str("run_crossover", kHighsOffString);
        ipm.set_option_value_bool("run_centring", true);
        if h16_ac_centring {
            ipm.set_option_value_bool("h16_ac_centring", true);
        }
        ipm.pass_model(two_col_lp());
        ipm.optimize_lp();
        ipm.get_info().centring_status
    }

    #[test]
    fn h16_ac_centring_on_ipx_records_centring_status() {
        // Centring is run on this LP, so the record is Success or Failed, never None.
        let status = ipx_centring_status(true);
        assert_ne!(status, kCentringStatusNone);
    }

    #[test]
    fn h16_ac_centring_off_ipx_leaves_centring_status_none() {
        // 1.15.1 has no centring record: the new field stays at its default.
        assert_eq!(ipx_centring_status(false), kCentringStatusNone);
    }
}

/// Symmetry data moved out of `HighsMipSolverData` between the early and the HiGHS-position
/// detection finish.
pub(crate) struct SymmetryStash {
    symmetries: HighsSymmetries,
    global_orbits: Option<std::sync::Arc<StabilizerOrbits>>,
    symmetry_dynamic: Option<Box<crate::mip::symmetry_dynamic::S3SymData>>,
}
