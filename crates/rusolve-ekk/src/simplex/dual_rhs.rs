//! HiGHS source: simplex/HEkkDualRHS.h, simplex/HEkkDualRHS.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Dual simplex primal-infeasibility bookkeeping: CHUZR (`chooseNormal`) with DSE/Devex
//! merit, hyper-sparse infeasibility lists, primal updates.

use crate::prelude::*;

/// HiGHS: HEkkDualRHS.h::HEkkDualRHS
/// # Removed C++ aliases (pattern (a))
///
/// | C++ member      | Rust expression          |
/// |-----------------|--------------------------|
/// | `ekk_instance_` | `ekk` (method parameter) |
/// | `analysis`      | `ekk.analysis_`          |
///
/// `workMark` is C++ `vector<char>` -> `Vec<i8>`.
#[derive(Clone, Debug)]
pub struct HEkkDualRHS {
    pub work_cutoff: f64,
    pub work_count: HighsInt,
    pub work_mark: Vec<i8>,
    pub work_index: Vec<HighsInt>,
    pub work_infeasibility: Vec<f64>,
    pub part_num: HighsInt,
    pub part_num_row: HighsInt,
    pub part_num_col: HighsInt,
    pub part_num_cut: HighsInt,
    pub part_switch: HighsInt,
    pub work_partition: Vec<HighsInt>,
}

impl Default for HEkkDualRHS {
    fn default() -> Self {
        Self {
            work_cutoff: Default::default(),
            work_count: Default::default(),
            work_mark: Default::default(),
            work_index: Default::default(),
            work_infeasibility: Default::default(),
            part_num: Default::default(),
            part_num_row: Default::default(),
            part_num_col: Default::default(),
            part_num_cut: Default::default(),
            part_switch: Default::default(),
            work_partition: Default::default(),
        }
    }
}

impl HEkkDualRHS {
    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::setup (lines 27-35, 9 lines)
    pub fn setup(&mut self, ekk: &mut HEkk) {
        let num_row = ekk.lp_.num_row_;
        self.work_mark.resize(num_row as usize, 0);
        self.work_index.resize(num_row as usize, 0);
        self.work_infeasibility.resize(num_row as usize, 0.0);
        self.part_num = 0;
        self.part_switch = 0;
        // analysis = &ekk_instance_.analysis_ is a removed alias (pattern (a))
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::chooseNormal (lines 37-122, 86 lines)
    /// MUST-PORT (determinism, random stream): CHUZR random start `ekk.random_.integer(n)`; one draw per call branch exactly as in C++.
    pub fn choose_normal(&mut self, ekk: &mut HEkk, ch_index: &mut HighsInt) {
        // NB Must still call HighsInt to maintain sequence of random numbers
        // for code reproducibility!!
        if self.work_count == 0 {
            *ch_index = -1;
            return;
        }

        // Since chooseNormal calls itself, only start the clock if it's not
        // currently running
        let keep_timer_running = ekk.analysis_.simplex_timer_running(ChuzrDualClock, 0);
        if !keep_timer_running {
            ekk.analysis_.simplex_timer_start(ChuzrDualClock, 0);
        }

        if self.work_count < 0 {
            // DENSE mode
            let num_row = -self.work_count;
            // MUST-PORT: random stream draw
            let random_start = ekk.random_.integer_sup(num_row);
            let mut best_merit: f64 = 0.0;
            let mut best_index: HighsInt = -1;
            for section in 0..2 {
                let start = if section == 0 { random_start } else { 0 };
                let end = if section == 0 { num_row } else { random_start };
                for i_row in start..end {
                    if self.work_infeasibility[i_row as usize] > kHighsZero {
                        let my_infeas = self.work_infeasibility[i_row as usize];
                        let my_weight = ekk.dual_edge_weight_[i_row as usize];
                        if best_merit * my_weight < my_infeas {
                            best_merit = my_infeas / my_weight;
                            best_index = i_row;
                        }
                    }
                }
            }
            *ch_index = best_index;
        } else {
            // SPARSE mode
            // MUST-PORT: random stream draw
            let random_start = ekk.random_.integer_sup(self.work_count);
            let mut best_merit: f64 = 0.0;
            let mut best_index: HighsInt = -1;
            for section in 0..2 {
                let start = if section == 0 { random_start } else { 0 };
                let end = if section == 0 {
                    self.work_count
                } else {
                    random_start
                };
                for i in start..end {
                    let i_row = self.work_index[i as usize];
                    if self.work_infeasibility[i_row as usize] > kHighsZero {
                        let my_infeas = self.work_infeasibility[i_row as usize];
                        let my_weight = ekk.dual_edge_weight_[i_row as usize];
                        if best_merit * my_weight < my_infeas {
                            best_merit = my_infeas / my_weight;
                            best_index = i_row;
                        }
                    }
                }
            }

            let mut create_list_again: HighsInt = 0;
            if best_index == -1 {
                create_list_again = (self.work_cutoff > 0.0) as HighsInt;
            } else if best_merit <= self.work_cutoff * 0.99 {
                create_list_again = 1;
            }
            if create_list_again != 0 {
                self.create_infeas_list(ekk, 0.0);
                self.choose_normal(ekk, &mut best_index);
            }
            *ch_index = best_index;
        }
        // Since chooseNormal calls itself, only stop the clock if it's not currently
        // running
        if !keep_timer_running {
            ekk.analysis_.simplex_timer_stop(ChuzrDualClock, 0);
        }
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::chooseMultiGlobal (lines 124-212, 89 lines)
    /// PAMI CHUZR (`HEkkDual::majorChooseRow`).
    ///
    /// C++ `HighsInt* chIndex` is an output array of `chLimit` entries (the caller passes
    /// `choiceIndex.data()`, HEkkDualMulti.cpp:87) -> `ch_index: &mut [HighsInt]`;
    /// C++ `HighsInt* chCount` is a scalar out-parameter -> `ch_count: &mut HighsInt`.
    /// MUST-PORT (determinism, random stream): one `ekk.random_.integer(n)` draw in DENSE
    /// mode, and in SPARSE mode only when `workCount != 0`, exactly as in C++.
    pub fn choose_multi_global(
        &mut self,
        ekk: &mut HEkk,
        ch_index: &mut [HighsInt],
        ch_count: &mut HighsInt,
        ch_limit: HighsInt,
    ) {
        ekk.analysis_.simplex_timer_start(ChuzrDualClock, 0);

        for i in 0..ch_limit {
            ch_index[i as usize] = -1;
        }

        let choose_check: HighsUInt = (ch_limit * 2) as HighsUInt;
        let mut set_p: Vec<(f64, HighsInt)> = Vec::with_capacity(choose_check as usize);

        if self.work_count < 0 {
            // DENSE mode
            let num_row = -self.work_count;
            // MUST-PORT: random stream draw
            let random_start = ekk.random_.integer_sup(num_row);
            let mut cutoff_merit: f64 = 0.0;
            // Now
            for section in 0..2 {
                let start = if section == 0 { random_start } else { 0 };
                let end = if section == 0 { num_row } else { random_start };
                for i_row in start..end {
                    if self.work_infeasibility[i_row as usize] > kHighsZero {
                        let my_infeas = self.work_infeasibility[i_row as usize];
                        let my_weight = ekk.dual_edge_weight_[i_row as usize];
                        if cutoff_merit * my_weight < my_infeas {
                            // Save
                            set_p.push((-my_infeas / my_weight, i_row));
                            // Shrink
                            if set_p.len() >= choose_check as usize {
                                pdqsort(&mut set_p[..]);
                                set_p.resize(ch_limit as usize, (0.0, 0));
                                cutoff_merit = -set_p.last().expect("chLimit > 0").0;
                            }
                        }
                    }
                }
            }
        } else {
            // SPARSE Mode
            let random_start: HighsInt = if self.work_count != 0 {
                // MUST-PORT: random stream draw
                ekk.random_.integer_sup(self.work_count)
            } else {
                // workCount = 0
                0
            };
            let mut cutoff_merit: f64 = 0.0;
            // Now
            for section in 0..2 {
                let start = if section == 0 { random_start } else { 0 };
                let end = if section == 0 {
                    self.work_count
                } else {
                    random_start
                };
                for i in start..end {
                    let i_row = self.work_index[i as usize];
                    if self.work_infeasibility[i_row as usize] > kHighsZero {
                        let my_infeas = self.work_infeasibility[i_row as usize];
                        let my_weight = ekk.dual_edge_weight_[i_row as usize];
                        if cutoff_merit * my_weight < my_infeas {
                            // Save
                            set_p.push((-my_infeas / my_weight, i_row));
                            // Shrink
                            if set_p.len() >= choose_check as usize {
                                pdqsort(&mut set_p[..]);
                                set_p.resize(ch_limit as usize, (0.0, 0));
                                cutoff_merit = -set_p.last().expect("chLimit > 0").0;
                            }
                        }
                    }
                }
            }
        }

        // Store the setP
        pdqsort(&mut set_p[..]);
        if (set_p.len() as HighsInt) > ch_limit {
            set_p.resize(ch_limit as usize, (0.0, 0));
        }
        *ch_count = set_p.len() as HighsInt;
        for i in 0..set_p.len() {
            ch_index[i] = set_p[i].1;
        }
        ekk.analysis_.simplex_timer_stop(ChuzrDualClock, 0);
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::chooseMultiHyperGraphAuto (lines 214-222, 9 lines)
    /// PAMI CHUZR (`HEkkDual::majorChooseRow`).
    pub fn choose_multi_hyper_graph_auto(
        &mut self,
        ekk: &mut HEkk,
        ch_index: &mut [HighsInt],
        ch_count: &mut HighsInt,
        ch_limit: HighsInt,
    ) {
        // Automatically decide to use partition or not
        if self.part_switch != 0 {
            self.choose_multi_hyper_graph_part(ekk, ch_index, ch_count, ch_limit);
        } else {
            self.choose_multi_global(ekk, ch_index, ch_count, ch_limit);
        }
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::chooseMultiHyperGraphPart (lines 224-306, 83 lines)
    /// PAMI CHUZR (`HEkkDual::majorChooseRow`).
    pub fn choose_multi_hyper_graph_part(
        &mut self,
        ekk: &mut HEkk,
        ch_index: &mut [HighsInt],
        ch_count: &mut HighsInt,
        ch_limit: HighsInt,
    ) {
        ekk.analysis_.simplex_timer_start(ChuzrDualClock, 0);

        // Force to use partition method, unless doesn't exist
        if self.part_num != ch_limit {
            self.choose_multi_global(ekk, ch_index, ch_count, ch_limit);
            self.part_switch = 0;
            ekk.analysis_.simplex_timer_stop(ChuzrDualClock, 0);
            return;
        }

        // Initialise
        for i in 0..ch_limit {
            ch_index[i as usize] = -1;
        }
        *ch_count = 0;

        if self.work_count < 0 {
            // DENSE mode
            let num_row = -self.work_count;
            // MUST-PORT: random stream draw
            let random_start = ekk.random_.integer_sup(num_row);
            let mut best_merit: Vec<f64> = vec![0.0; ch_limit as usize];
            let mut best_index: Vec<HighsInt> = vec![-1; ch_limit as usize];
            for section in 0..2 {
                let start = if section == 0 { random_start } else { 0 };
                let end = if section == 0 { num_row } else { random_start };
                for i_row in start..end {
                    if self.work_infeasibility[i_row as usize] > kHighsZero {
                        let i_part = self.work_partition[i_row as usize];
                        let my_infeas = self.work_infeasibility[i_row as usize];
                        let my_weight = ekk.dual_edge_weight_[i_row as usize];
                        if best_merit[i_part as usize] * my_weight < my_infeas {
                            best_merit[i_part as usize] = my_infeas / my_weight;
                            best_index[i_part as usize] = i_row;
                        }
                    }
                }
            }
            let mut count: HighsInt = 0;
            for i in 0..ch_limit {
                if best_index[i as usize] != -1 {
                    ch_index[count as usize] = best_index[i as usize];
                    count += 1;
                }
            }
            *ch_count = count;
        } else {
            // SPARSE mode
            if self.work_count == 0 {
                ekk.analysis_.simplex_timer_stop(ChuzrDualClock, 0);
                return;
            }

            // MUST-PORT: random stream draw
            let random_start = ekk.random_.integer_sup(self.work_count);
            let mut best_merit: Vec<f64> = vec![0.0; ch_limit as usize];
            let mut best_index: Vec<HighsInt> = vec![-1; ch_limit as usize];
            for section in 0..2 {
                let start = if section == 0 { random_start } else { 0 };
                let end = if section == 0 {
                    self.work_count
                } else {
                    random_start
                };
                for i in start..end {
                    let i_row = self.work_index[i as usize];
                    if self.work_infeasibility[i_row as usize] > kHighsZero {
                        let i_part = self.work_partition[i_row as usize];
                        let my_infeas = self.work_infeasibility[i_row as usize];
                        let my_weight = ekk.dual_edge_weight_[i_row as usize];
                        if best_merit[i_part as usize] * my_weight < my_infeas {
                            best_merit[i_part as usize] = my_infeas / my_weight;
                            best_index[i_part as usize] = i_row;
                        }
                    }
                }
            }
            let mut count: HighsInt = 0;
            for i in 0..ch_limit {
                if best_index[i as usize] != -1 {
                    ch_index[count as usize] = best_index[i as usize];
                    count += 1;
                }
            }
            *ch_count = count;
        }

        ekk.analysis_.simplex_timer_stop(ChuzrDualClock, 0);
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::updatePrimal (lines 308-351, 44 lines)
    pub fn update_primal(&mut self, ekk: &mut HEkk, column: &HVector, theta: f64) -> bool {
        ekk.analysis_.simplex_timer_start(UpdatePrimalClock, 0);

        let num_row = ekk.lp_.num_row_;
        let column_count = column.count;

        let tp = ekk.options_.primal_feasibility_tolerance;

        let mut num_excessive_primal: HighsInt = 0;
        let update_primal_in_dense = column_count < 0 || column_count as f64 > 0.4 * num_row as f64;

        let to_entry = if update_primal_in_dense {
            num_row
        } else {
            column_count
        };
        for i_entry in 0..to_entry {
            let i_row = if update_primal_in_dense {
                i_entry
            } else {
                column.index[i_entry as usize]
            };
            ekk.info_.base_value_[i_row as usize] -= theta * column.array[i_row as usize];
            // @primal_infeasibility calculation
            let value = ekk.info_.base_value_[i_row as usize];
            let lower = ekk.info_.base_lower_[i_row as usize];
            let upper = ekk.info_.base_upper_[i_row as usize];
            let mut primal_infeasibility: f64 = 0.0;
            if value < lower - tp {
                primal_infeasibility = lower - value;
            } else if value > upper + tp {
                primal_infeasibility = value - upper;
            }
            if ekk.info_.store_squared_primal_infeasibility {
                self.work_infeasibility[i_row as usize] =
                    primal_infeasibility * primal_infeasibility;
            } else {
                self.work_infeasibility[i_row as usize] = primal_infeasibility.abs();
            }
            if ekk.info_.base_value_[i_row as usize] <= -kExcessivePrimalValue
                || ekk.info_.base_value_[i_row as usize] >= kExcessivePrimalValue
            {
                num_excessive_primal += 1;
            }
        }
        ekk.analysis_.simplex_timer_stop(UpdatePrimalClock, 0);
        // Flag detection of excessive values in return
        if num_excessive_primal != 0 {
            return false;
        }
        true
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::updatePivots (lines 353-373, 21 lines)
    pub fn update_pivots(&mut self, ekk: &mut HEkk, i_row: HighsInt, value: f64) {
        let tp = ekk.options_.primal_feasibility_tolerance;
        ekk.info_.base_value_[i_row as usize] = value;
        // @primal_infeasibility calculation
        let lower = ekk.info_.base_lower_[i_row as usize];
        let upper = ekk.info_.base_upper_[i_row as usize];
        let mut primal_infeasibility: f64 = 0.0;
        if value < lower - tp {
            primal_infeasibility = lower - value;
        } else if value > upper + tp {
            primal_infeasibility = value - upper;
        }
        if ekk.info_.store_squared_primal_infeasibility {
            self.work_infeasibility[i_row as usize] = primal_infeasibility * primal_infeasibility;
        } else {
            self.work_infeasibility[i_row as usize] = primal_infeasibility.abs();
        }
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::updateInfeasList (lines 375-410, 36 lines)
    pub fn update_infeas_list(&mut self, ekk: &mut HEkk, column: &HVector) {
        let column_count = column.count;

        // DENSE mode: disabled
        if self.work_count < 0 {
            return;
        }

        ekk.analysis_.simplex_timer_start(UpdatePrimalClock, 0);

        // work_count is kept in a local across the loop and written back once; the
        // C++ updates the member in place, but the list contents and order are identical.
        let num_entries = column_count.max(0) as usize;
        let mut work_count = self.work_count;
        if self.work_cutoff <= 0.0 {
            // The regular sparse way
            for &i_row in &column.index[..num_entries] {
                let row = i_row as usize;
                if self.work_mark[row] == 0 && self.work_infeasibility[row] != 0.0 {
                    self.work_index[work_count as usize] = i_row;
                    work_count += 1;
                    self.work_mark[row] = 1;
                }
            }
        } else {
            // The hyper sparse way
            for &i_row in &column.index[..num_entries] {
                let row = i_row as usize;
                if self.work_mark[row] == 0
                    && self.work_infeasibility[row] > ekk.dual_edge_weight_[row] * self.work_cutoff
                {
                    self.work_index[work_count as usize] = i_row;
                    work_count += 1;
                    self.work_mark[row] = 1;
                }
            }
        }
        self.work_count = work_count;

        ekk.analysis_.simplex_timer_stop(UpdatePrimalClock, 0);
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::createArrayOfPrimalInfeasibilities (lines 412-434, 23 lines)
    pub fn create_array_of_primal_infeasibilities(&mut self, ekk: &mut HEkk) {
        let num_row = ekk.lp_.num_row_;
        let tp = ekk.options_.primal_feasibility_tolerance;
        for i in 0..num_row {
            // @primal_infeasibility calculation
            let value = ekk.info_.base_value_[i as usize];
            let lower = ekk.info_.base_lower_[i as usize];
            let upper = ekk.info_.base_upper_[i as usize];
            let mut primal_infeasibility: f64 = 0.0;
            if value < lower - tp {
                primal_infeasibility = lower - value;
            } else if value > upper + tp {
                primal_infeasibility = value - upper;
            }
            if ekk.info_.store_squared_primal_infeasibility {
                self.work_infeasibility[i as usize] = primal_infeasibility * primal_infeasibility;
            } else {
                self.work_infeasibility[i as usize] = primal_infeasibility.abs();
            }
        }
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::createInfeasList (lines 436-498, 63 lines)
    pub fn create_infeas_list(&mut self, ekk: &mut HEkk, column_density: f64) {
        let num_row = ekk.lp_.num_row_;

        // 1. Build the full list
        for m in self.work_mark[..num_row as usize].iter_mut() {
            *m = 0;
        }
        self.work_count = 0;
        self.work_cutoff = 0.0;
        for i_row in 0..num_row {
            if self.work_infeasibility[i_row as usize] != 0.0 {
                self.work_mark[i_row as usize] = 1;
                self.work_index[self.work_count as usize] = i_row;
                self.work_count += 1;
            }
        }

        // 2. See if it worth to try to go sparse
        //    (Many candidates, really sparse RHS)
        if self.work_count as f64 > cpp_max(num_row as f64 * 0.01, 500.0) && column_density < 0.05 {
            let icutoff = cpp_max(self.work_count as f64 * 0.001, 500.0) as HighsInt;
            let mut max_merit: f64 = 0.0;
            let mut i_put: usize = 0;
            for i_row in 0..num_row {
                if self.work_mark[i_row as usize] != 0 {
                    let my_merit = self.work_infeasibility[i_row as usize]
                        / ekk.dual_edge_weight_[i_row as usize];
                    if max_merit < my_merit {
                        max_merit = my_merit;
                    }
                    ekk.scattered_dual_edge_weight_[i_put] = -my_merit;
                    i_put += 1;
                }
            }
            // nth_element(dwork, dwork + icutoff, dwork + workCount): only the value at
            // position icutoff is read afterwards, and it is order-independent.
            ekk.scattered_dual_edge_weight_[..self.work_count as usize]
                .select_nth_unstable_by(icutoff as usize, |a, b| {
                    a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
                });
            let cut_merit = -ekk.scattered_dual_edge_weight_[icutoff as usize];
            self.work_cutoff = cpp_min(max_merit * 0.99999, cut_merit * 1.00001);

            // Create again
            for m in self.work_mark[..num_row as usize].iter_mut() {
                *m = 0;
            }
            self.work_count = 0;
            for i_row in 0..num_row {
                if self.work_infeasibility[i_row as usize]
                    >= ekk.dual_edge_weight_[i_row as usize] * self.work_cutoff
                {
                    self.work_index[self.work_count as usize] = i_row;
                    self.work_count += 1;
                    self.work_mark[i_row as usize] = 1;
                }
            }

            // Reduce by drop smaller
            if self.work_count as f64 > icutoff as f64 * 1.5 {
                // Firstly take up "icutoff" number of elements
                let full_count = self.work_count;
                self.work_count = icutoff;
                for i in icutoff..full_count {
                    let i_row = self.work_index[i as usize];
                    if self.work_infeasibility[i_row as usize]
                        > ekk.dual_edge_weight_[i_row as usize] * cut_merit
                    {
                        self.work_index[self.work_count as usize] = i_row;
                        self.work_count += 1;
                    } else {
                        self.work_mark[i_row as usize] = 0;
                    }
                }
            }
        }

        // 3. If there are still too many candidates: disable them
        if self.work_count as f64 > 0.2 * num_row as f64 {
            self.work_count = -num_row;
            self.work_cutoff = 0.0;
        }
    }

    /// HiGHS: HEkkDualRHS.cpp::HEkkDualRHS::assessOptimality (lines 500-535, 36 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn assess_optimality(&mut self, ekk: &mut HEkk) {
        let mut num_work_infeasibilities: HighsInt = 0;
        let mut max_work_infeasibility: f64 = 0.0;
        let num_row = ekk.lp_.num_row_;
        for i_row in 0..num_row {
            if self.work_infeasibility[i_row as usize] > kHighsZero {
                num_work_infeasibilities += 1;
                max_work_infeasibility = cpp_max(
                    self.work_infeasibility[i_row as usize],
                    max_work_infeasibility,
                );
            }
        }
        ekk.compute_simplex_primal_infeasible();
        let num_primal_infeasibilities = ekk.info_.num_primal_infeasibilities;
        let max_primal_infeasibility = ekk.info_.max_primal_infeasibility;
        let regular_report = false;
        if regular_report || (num_work_infeasibilities != 0 && num_primal_infeasibilities == 0) {
            println!(
                "assessOptimality: {:6} rows; workCount = {:4} ({:6.4}) num / max infeasibilities: work = {:4} / {:11.4e}; simplex = {:4} / {:11.4e}: {}",
                num_row,
                self.work_count,
                if self.work_count > 0 {
                    (1.0 * self.work_count as f64) / num_row as f64
                } else {
                    0.0
                },
                num_work_infeasibilities,
                max_work_infeasibility,
                num_primal_infeasibilities,
                max_primal_infeasibility,
                if num_primal_infeasibilities == 0 {
                    "Optimal"
                } else {
                    ""
                }
            );
            if num_work_infeasibilities != 0 && num_primal_infeasibilities == 0 {
                println!(
                    "assessOptimality: call {}; tick {}; iter {}",
                    ekk.debug_solve_call_num_,
                    ekk.debug_initial_build_synthetic_tick_,
                    ekk.iteration_count_
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use super::*;

    /// A dense (`workCount = -numRow`) RHS over `infeas` / `weight`, with a seeded stream.
    fn dense_case(infeas: &[f64], weight: &[f64], seed: HighsUInt) -> (HEkkDualRHS, HEkk) {
        let mut rhs = HEkkDualRHS::default();
        rhs.work_infeasibility = infeas.to_vec();
        rhs.work_index = vec![0; infeas.len()];
        rhs.work_mark = vec![0; infeas.len()];
        rhs.work_count = -(infeas.len() as HighsInt);
        let mut ekk = HEkk::default();
        ekk.dual_edge_weight_ = weight.to_vec();
        ekk.random_ = HighsRandom::new(seed);
        (rhs, ekk)
    }

    /// The stream after the call must equal `expected` (same number of draws).
    fn assert_same_stream(ekk: &mut HEkk, mut expected: HighsRandom) {
        for _ in 0..4 {
            assert_eq!(ekk.random_.integer(), expected.integer());
        }
    }

    #[test]
    fn w5_pami_rhs_choose_multi_row_sets() {
        // `chooseMultiGlobal` (HEkkDualRHS.cpp 124-212), DENSE mode, chLimit = 3
        // (chooseCHECK = 6, so the shrink at 154-158 runs): the 3 rows of largest
        // merit infeas / weight, in decreasing merit order (sort of (-merit, iRow)).
        // Merits: 0:0.5 1:4 2:0 3:3 4:0.25 5:6 6:1 7:0 8:5 9:2 10:0.1
        let infeas = [1.0, 4.0, 0.0, 6.0, 1.0, 6.0, 1.0, 0.0, 5.0, 2.0, 0.1];
        let weight = [2.0, 1.0, 1.0, 2.0, 4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        for seed in 0..8 {
            let (mut rhs, mut ekk) = dense_case(&infeas, &weight, seed);
            let mut expected = ekk.random_.clone();
            expected.integer_sup(infeas.len() as HighsInt);
            let mut ch_index = vec![7; 3];
            let mut ch_count: HighsInt = -5;
            rhs.choose_multi_global(&mut ekk, &mut ch_index, &mut ch_count, 3);
            assert_eq!(ch_count, 3);
            assert_eq!(ch_index, vec![5, 8, 1]);
            assert_same_stream(&mut ekk, expected);
        }

        // Fewer candidates than chLimit: the tail keeps the -1 written at line 128.
        let (mut rhs, mut ekk) = dense_case(&[0.0, 2.0, 0.0, 3.0], &[1.0; 4], 3);
        let mut ch_index = vec![9; 4];
        let mut ch_count: HighsInt = -5;
        rhs.choose_multi_global(&mut ekk, &mut ch_index, &mut ch_count, 4);
        assert_eq!(ch_count, 2);
        assert_eq!(ch_index, vec![3, 1, -1, -1]);

        // Ties in merit: std::pair order breaks them by the smaller row index.
        let (mut rhs, mut ekk) = dense_case(&[1.0, 2.0, 1.0, 2.0, 1.0], &[1.0; 5], 5);
        let mut ch_index = vec![0; 4];
        let mut ch_count: HighsInt = 0;
        rhs.choose_multi_global(&mut ekk, &mut ch_index, &mut ch_count, 4);
        assert_eq!(ch_count, 4);
        assert_eq!(ch_index, vec![1, 3, 0, 2]);
    }

    #[test]
    fn w5_pami_rhs_choose_multi_global_sparse() {
        // SPARSE mode walks workIndex[0..workCount) (lines 170-202).
        let infeas = [3.0, 0.0, 8.0, 1.0, 9.0, 4.0];
        let (mut rhs, mut ekk) = dense_case(&infeas, &[1.0; 6], 11);
        rhs.work_index = vec![3, 0, 5, 2, 0, 0];
        rhs.work_count = 4; // rows 3, 0, 5, 2 (row 4 is not listed)
        let mut expected = ekk.random_.clone();
        expected.integer_sup(4);
        let mut ch_index = vec![0; 2];
        let mut ch_count: HighsInt = 0;
        rhs.choose_multi_global(&mut ekk, &mut ch_index, &mut ch_count, 2);
        assert_eq!(ch_count, 2);
        assert_eq!(ch_index, vec![2, 5]);
        assert_same_stream(&mut ekk, expected);

        // workCount = 0: no random draw (lines 166-169), nothing chosen.
        let (mut rhs, mut ekk) = dense_case(&infeas, &[1.0; 6], 11);
        rhs.work_count = 0;
        let expected = ekk.random_.clone();
        let mut ch_index = vec![4; 3];
        let mut ch_count: HighsInt = 9;
        rhs.choose_multi_global(&mut ekk, &mut ch_index, &mut ch_count, 3);
        assert_eq!(ch_count, 0);
        assert_eq!(ch_index, vec![-1, -1, -1]);
        assert_same_stream(&mut ekk, expected);
    }

    #[test]
    fn w5_pami_rhs_choose_multi_hyper_graph() {
        let infeas = [1.0, 5.0, 2.0, 0.0, 7.0, 3.0, 0.0, 4.0];
        let partition = vec![0, 0, 2, 1, 2, 0, 1, 2];

        // Auto with partSwitch = 0 is chooseMultiGlobal (line 221).
        let (mut rhs, mut ekk) = dense_case(&infeas, &[1.0; 8], 2);
        rhs.work_partition = partition.clone();
        rhs.part_num = 3;
        let mut ch_index = vec![0; 3];
        let mut ch_count: HighsInt = 0;
        rhs.choose_multi_hyper_graph_auto(&mut ekk, &mut ch_index, &mut ch_count, 3);
        assert_eq!((ch_count, ch_index), (3, vec![4, 1, 7]));

        // Part, DENSE: best merit per partition, written in partition order; partition 1
        // has no infeasible row, so it is skipped (lines 260-266).
        let (mut rhs, mut ekk) = dense_case(&infeas, &[1.0; 8], 2);
        rhs.work_partition = partition.clone();
        rhs.part_num = 3;
        rhs.part_switch = 1;
        let mut expected = ekk.random_.clone();
        expected.integer_sup(8);
        let mut ch_index = vec![9; 3];
        let mut ch_count: HighsInt = 0;
        rhs.choose_multi_hyper_graph_auto(&mut ekk, &mut ch_index, &mut ch_count, 3);
        assert_eq!((ch_count, &ch_index[..]), (2, &[1, 4, -1][..]));
        assert_eq!(rhs.part_switch, 1);
        assert_same_stream(&mut ekk, expected);

        // Part, DENSE, a tie inside a partition: the first row met from randomStart wins
        // (strict `<`, line 252).
        let tie = [2.0, 2.0, 2.0, 2.0];
        for seed in 0..8 {
            let (mut rhs, mut ekk) = dense_case(&tie, &[1.0; 4], seed);
            rhs.work_partition = vec![0; 4];
            rhs.part_num = 1;
            let random_start = ekk.random_.clone().integer_sup(4);
            let mut ch_index = vec![9; 1];
            let mut ch_count: HighsInt = 0;
            rhs.choose_multi_hyper_graph_part(&mut ekk, &mut ch_index, &mut ch_count, 1);
            assert_eq!((ch_count, ch_index[0]), (1, random_start));
        }

        // Part, SPARSE over workIndex.
        let (mut rhs, mut ekk) = dense_case(&infeas, &[1.0; 8], 4);
        rhs.work_partition = partition.clone();
        rhs.part_num = 3;
        rhs.work_index = vec![0, 2, 3, 5, 0, 0, 0, 0];
        rhs.work_count = 4;
        let mut ch_index = vec![9; 3];
        let mut ch_count: HighsInt = 0;
        rhs.choose_multi_hyper_graph_part(&mut ekk, &mut ch_index, &mut ch_count, 3);
        assert_eq!((ch_count, &ch_index[..]), (2, &[5, 2, -1][..]));

        // Part, SPARSE, workCount = 0: early return after the initialisation, no draw.
        let (mut rhs, mut ekk) = dense_case(&infeas, &[1.0; 8], 4);
        rhs.work_partition = partition.clone();
        rhs.part_num = 3;
        rhs.work_count = 0;
        let expected = ekk.random_.clone();
        let mut ch_index = vec![9; 3];
        let mut ch_count: HighsInt = 6;
        rhs.choose_multi_hyper_graph_part(&mut ekk, &mut ch_index, &mut ch_count, 3);
        assert_eq!((ch_count, &ch_index[..]), (0, &[-1, -1, -1][..]));
        assert_same_stream(&mut ekk, expected);

        // Part with partNum != chLimit: chooseMultiGlobal, then partSwitch = 0 (229-235).
        let (mut rhs, mut ekk) = dense_case(&infeas, &[1.0; 8], 2);
        rhs.work_partition = partition;
        rhs.part_num = 3;
        rhs.part_switch = 1;
        let mut ch_index = vec![0; 2];
        let mut ch_count: HighsInt = 0;
        rhs.choose_multi_hyper_graph_auto(&mut ekk, &mut ch_index, &mut ch_count, 2);
        assert_eq!((ch_count, ch_index), (2, vec![4, 1]));
        assert_eq!(rhs.part_switch, 0);
    }
}
