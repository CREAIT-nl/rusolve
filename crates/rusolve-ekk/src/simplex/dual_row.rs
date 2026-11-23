//! HiGHS source: simplex/HEkkDualRow.h, simplex/HEkkDualRow.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Dual simplex pivotal row: pack, CHUZC candidate selection with BFRT (Bound Flipping
//! Ratio Test) grouping, the Harris-style final choice, flips and dual updates.

use crate::prelude::*;

// HEkkDualRow.h file-scope constants
const kInitialTotalChange: f64 = 1e-12;
const kInitialRemainTheta: f64 = 1e100;
const kMaxSelectTheta: f64 = 1e18;

/// HiGHS: HEkkDualRow.h::HEkkDualRow
/// # Removed C++ aliases (pattern (a))
///
/// | C++ member              | Rust expression                         |
/// |-------------------------|-----------------------------------------|
/// | `ekk_instance_`         | `ekk` (method parameter)                |
/// | `workMove`              | `ekk.basis_.nonbasic_move_`             |
/// | `workDual`              | `ekk.info_.work_dual_`                  |
/// | `workRange`             | `ekk.info_.work_range_`                 |
/// | `work_devex_index`      | `ekk.info_.devex_index_`                |
/// | `workNumTotPermutation` | `ekk.info_.num_tot_permutation_`        |
/// | `analysis`              | `ekk.analysis_`                         |
///
/// `std::set<HighsInt> freeList` is a `BTreeSet<HighsInt>` (same iteration order).
/// `workData` is `Vec<(HighsInt, f64)>` (`.first` -> `.0`, `.second` -> `.1`).
#[derive(Clone, Debug)]
pub struct HEkkDualRow {
    pub work_size: HighsInt,
    pub free_list: BTreeSet<HighsInt>,
    pub pack_count: HighsInt,
    pub pack_index: Vec<HighsInt>,
    pub pack_value: Vec<f64>,
    pub computed_edge_weight: f64,
    pub work_delta: f64,
    pub work_alpha: f64,
    pub work_theta: f64,
    pub work_pivot: HighsInt,
    pub work_count: HighsInt,
    pub work_data: Vec<(HighsInt, f64)>,
    pub work_group: Vec<HighsInt>,
    pub alt_work_count: HighsInt,
    pub original_work_data: Vec<(HighsInt, f64)>,
    pub sorted_work_data: Vec<(HighsInt, f64)>,
    pub alt_work_group: Vec<HighsInt>,
}

impl Default for HEkkDualRow {
    fn default() -> Self {
        Self {
            work_size: -1,
            free_list: Default::default(),
            pack_count: 0,
            pack_index: Default::default(),
            pack_value: Default::default(),
            computed_edge_weight: 0.0,
            work_delta: 0.0,
            work_alpha: 0.0,
            work_theta: 0.0,
            work_pivot: 0,
            work_count: 0,
            work_data: Default::default(),
            work_group: Default::default(),
            alt_work_count: 0,
            original_work_data: Default::default(),
            sorted_work_data: Default::default(),
            alt_work_group: Default::default(),
        }
    }
}

impl HEkkDualRow {
    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::setupSlice (lines 27-42, 16 lines)
    pub fn setup_slice(&mut self, ekk: &mut HEkk, size: HighsInt) {
        self.work_size = size;

        // Allocate spaces
        self.pack_count = 0;
        self.pack_index.resize(self.work_size as usize, 0);
        self.pack_value.resize(self.work_size as usize, 0.0);

        self.work_count = 0;
        self.work_data.resize(self.work_size as usize, (0, 0.0));
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::setup (lines 44-56, 13 lines)
    /// MUST-PORT (determinism, random stream): `workNumTotPermutation` alias -> `ekk.info_.num_tot_permutation_`.
    pub fn setup(&mut self, ekk: &mut HEkk) {
        // Setup common vectors
        let num_tot = ekk.lp_.num_col_ + ekk.lp_.num_row_;
        self.setup_slice(ekk, num_tot);

        // deleteFreelist() is being called in Phase 1 and Phase 2 since
        // it's in updatePivots(), but create_Freelist() is only called in
        // Phase 2. Hence freeList is not initialised when freeList.empty()
        // is used in deleteFreelist(), clear freeList now.
        self.free_list.clear();
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::clear (lines 58-61, 4 lines)
    pub fn clear(&mut self, ekk: &mut HEkk) {
        self.pack_count = 0;
        self.work_count = 0;
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::chooseMakepack (lines 63-78, 16 lines)
    pub fn choose_makepack(&mut self, ekk: &mut HEkk, row: &HVector, offset: HighsInt) {
        // Pack the indices and values for the row
        //
        // Offset of numCol is used when packing row_ep
        let row_count = row.count;
        for i in 0..row_count {
            let index = row.index[i as usize];
            let value = row.array[index as usize];
            self.pack_index[self.pack_count as usize] = index + offset;
            self.pack_value[self.pack_count as usize] = value;
            self.pack_count += 1;
        }
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::choosePossible (lines 80-103, 24 lines)
    pub fn choose_possible(&mut self, ekk: &mut HEkk) {
        // Determine the possible variables - candidates for CHUZC
        let ta: f64 = if ekk.info_.update_count < 10 {
            1e-9
        } else if ekk.info_.update_count < 20 {
            3e-8
        } else {
            1e-6
        };
        let td: f64 = ekk.options_.dual_feasibility_tolerance;
        let move_out: HighsInt = if self.work_delta < 0.0 { -1 } else { 1 };
        // work_count and work_theta are locals across the scan and written back once;
        // the comparisons and the stored candidates are the same as the C++ member updates.
        let mut work_theta = kHighsInf;
        let mut work_count: HighsInt = 0;
        let pack_count = self.pack_count.max(0) as usize;
        for (&i_col, &pack_value) in self.pack_index[..pack_count]
            .iter()
            .zip(&self.pack_value[..pack_count])
        {
            let mv = ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt;
            let alpha = pack_value * (move_out as f64) * (mv as f64);
            if alpha > ta {
                self.work_data[work_count as usize] = (i_col, alpha);
                work_count += 1;
                let relax = ekk.info_.work_dual_[i_col as usize] * (mv as f64) + td;
                if work_theta * alpha > relax {
                    work_theta = relax / alpha;
                }
            }
        }
        self.work_theta = work_theta;
        self.work_count = work_count;
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::chooseJoinpack (lines 105-115, 11 lines)
    pub fn choose_joinpack(&mut self, ekk: &mut HEkk, other_row: &HEkkDualRow) {
        // Join pack of possible candidates in this row with possible
        // candidates in otherRow
        let other_count = other_row.work_count;
        for j in 0..other_count {
            self.work_data[(self.work_count + j) as usize] = other_row.work_data[j as usize];
        }
        self.work_count += other_count;
        self.work_theta = cpp_min(self.work_theta, other_row.work_theta);
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::chooseFinal (lines 117-315, 199 lines)
    pub fn choose_final(&mut self, ekk: &mut HEkk) -> HighsInt {
        // Chooses the entering variable via BFRT and EXPAND
        //
        // It will
        // (1) reduce the candidates as a small collection
        // (2) choose by BFRT by going over break points
        // (3) choose final by alpha
        // (4) determine final flip variables

        // 1. Reduce by large step BFRT
        ekk.analysis_.simplex_timer_start(Chuzc3Clock, 0);
        let mut report_bfrt = false;
        let debug_bfrt_report_on = false;
        if ekk.debug_iteration_report_ {
            report_bfrt = debug_bfrt_report_on;
            if report_bfrt {
                println!(
                    "HEkkDualRow::chooseFinal Check iter = {}",
                    ekk.iteration_count_
                );
            }
        }
        let mut full_count = self.work_count;
        self.work_count = 0;
        let mut total_change: f64 = 0.0;
        let total_delta = self.work_delta.abs();
        let mut select_theta = 10.0 * self.work_theta + 1e-7;
        loop {
            for i in self.work_count..full_count {
                let i_col = self.work_data[i as usize].0;
                let alpha = self.work_data[i as usize].1;
                let tight = (ekk.basis_.nonbasic_move_[i_col as usize] as f64)
                    * ekk.info_.work_dual_[i_col as usize];
                if alpha * select_theta >= tight {
                    self.work_data.swap(self.work_count as usize, i as usize);
                    self.work_count += 1;
                    total_change += ekk.info_.work_range_[i_col as usize] * alpha;
                }
            }
            select_theta *= 10.0;
            if total_change >= total_delta || self.work_count == full_count {
                break;
            }
        }
        ekk.analysis_.simplex_timer_stop(Chuzc3Clock, 0);
        // 2. Choose by small step BFRT

        // Use the quadratic cost sort for smaller values of workCount,
        // otherwise use the heap-based sort
        let use_quad_sort = true; // workCount < 100;
        let use_heap_sort = !use_quad_sort;
        debug_assert!(use_heap_sort || use_quad_sort);
        if use_quad_sort {
            ekk.analysis_.num_quad_chuzc += 1;
            ekk.analysis_.sum_quad_chuzc_size += self.work_count as f64;
            ekk.analysis_.max_quad_chuzc_size =
                cpp_max(self.work_count, ekk.analysis_.max_quad_chuzc_size);
        } else {
            ekk.analysis_.num_heap_chuzc += 1;
            ekk.analysis_.sum_heap_chuzc_size += self.work_count as f64;
            ekk.analysis_.max_heap_chuzc_size =
                cpp_max(self.work_count, ekk.analysis_.max_heap_chuzc_size);
        }

        if use_heap_sort {
            println!("CHUZC: Using heap sort");
            // Take a copy of workData and workCount for the independent
            // heap-based code
            self.original_work_data = self.work_data.clone();
            self.alt_work_count = self.work_count;
        }
        ekk.analysis_.simplex_timer_start(Chuzc4Clock, 0);
        let mut choose_ok = false;
        if use_quad_sort {
            // Use the O(n^2) quadratic sort for the candidates
            ekk.analysis_.simplex_timer_start(Chuzc4a0Clock, 0);
            choose_ok = self.choose_final_work_group_quad(ekk);
            //    if (!choose_ok) {
            //  choose_ok = quadChooseFinalWorkGroupQuad();
            //    }
            ekk.analysis_.simplex_timer_stop(Chuzc4a0Clock, 0);
        }
        if use_heap_sort {
            // Use the O(n log n) heap sort for the candidates
            ekk.analysis_.simplex_timer_start(Chuzc4a1Clock, 0);
            choose_ok = self.choose_final_work_group_heap(ekk);
            ekk.analysis_.simplex_timer_stop(Chuzc4a1Clock, 0);
        }
        if !choose_ok {
            ekk.analysis_.simplex_timer_stop(Chuzc4Clock, 0);
            return -1;
        }
        // Make sure that there is at least one group according to sorting procedure
        if use_quad_sort {
            debug_assert!(self.work_group.len() > 1);
        }
        if use_heap_sort {
            debug_assert!(self.alt_work_group.len() > 1);
        }

        // 3. Choose large alpha
        ekk.analysis_.simplex_timer_start(Chuzc4bClock, 0);
        let mut break_index: HighsInt = -1;
        let mut break_group: HighsInt = -1;
        let mut alt_break_index: HighsInt = -1;
        let mut alt_break_group: HighsInt = -1;
        // The vectors are moved out of `self` for the call (the callee takes
        // `&mut self` plus slices of them but only reads `ekk`) and moved back.
        if use_quad_sort {
            let pass_work_data = std::mem::take(&mut self.work_data);
            let pass_work_group = std::mem::take(&mut self.work_group);
            self.choose_final_large_alpha(
                ekk,
                &mut break_index,
                &mut break_group,
                self.work_count,
                &pass_work_data,
                &pass_work_group,
            );
            self.work_data = pass_work_data;
            self.work_group = pass_work_group;
        }
        if use_heap_sort {
            let pass_work_data = std::mem::take(&mut self.sorted_work_data);
            let pass_work_group = std::mem::take(&mut self.alt_work_group);
            self.choose_final_large_alpha(
                ekk,
                &mut alt_break_index,
                &mut alt_break_group,
                self.alt_work_count,
                &pass_work_data,
                &pass_work_group,
            );
            self.sorted_work_data = pass_work_data;
            self.alt_work_group = pass_work_group;
        }
        ekk.analysis_.simplex_timer_stop(Chuzc4bClock, 0);

        if !use_quad_sort {
            // If the quadratic sort is not being used, revert to the heap
            // sort results
            break_index = alt_break_index;
            break_group = alt_break_group;
        }
        ekk.analysis_.simplex_timer_start(Chuzc4cClock, 0);

        let move_out: HighsInt = if self.work_delta < 0.0 { -1 } else { 1 };
        debug_assert!(break_index >= 0);
        if use_quad_sort {
            self.work_pivot = self.work_data[break_index as usize].0;
            self.work_alpha = self.work_data[break_index as usize].1
                * (move_out as f64)
                * (ekk.basis_.nonbasic_move_[self.work_pivot as usize] as f64);
        } else {
            self.work_pivot = self.sorted_work_data[break_index as usize].0;
            self.work_alpha = self.sorted_work_data[break_index as usize].1
                * (move_out as f64)
                * (ekk.basis_.nonbasic_move_[self.work_pivot as usize] as f64);
        }
        if ekk.info_.work_dual_[self.work_pivot as usize]
            * (ekk.basis_.nonbasic_move_[self.work_pivot as usize] as f64)
            > 0.0
        {
            self.work_theta = ekk.info_.work_dual_[self.work_pivot as usize] / self.work_alpha;
        } else {
            self.work_theta = 0.0;
        }
        ekk.analysis_.simplex_timer_stop(Chuzc4cClock, 0);

        ekk.analysis_.simplex_timer_start(Chuzc4dClock, 0);

        // 4. Determine BFRT flip index: flip all
        full_count = break_index; // Not used
        let _ = full_count;
        self.work_count = 0;
        if use_quad_sort {
            if report_bfrt {
                // (the C++ printf here is omitted: it indexes workGroup[workGroup.size()])
                self.debug_report_bfrt_var(ekk, -1, &self.work_data);
            }
            for i in 0..self.work_group[break_group as usize] {
                if report_bfrt {
                    self.debug_report_bfrt_var(ekk, i, &self.work_data);
                }
                let i_col = self.work_data[i as usize].0;
                let mv = ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt;
                self.work_data[self.work_count as usize] =
                    (i_col, (mv as f64) * ekk.info_.work_range_[i_col as usize]);
                self.work_count += 1;
            }
            for i in
                self.work_group[break_group as usize]..self.work_group[(break_group + 1) as usize]
            {
                if report_bfrt {
                    self.debug_report_bfrt_var(ekk, i, &self.work_data);
                }
            }
        } else {
            println!(
                "DebugHeapSortCHUZC: Pivot = {:4}; alpha = {:11.4e}; theta = {:11.4e}",
                self.work_pivot, self.work_alpha, self.work_theta
            );
            self.debug_report_bfrt_var(ekk, -1, &self.sorted_work_data);
            for i in 0..self.alt_work_group[break_group as usize] {
                let i_col = self.sorted_work_data[i as usize].0;
                let mv = ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt;
                self.debug_report_bfrt_var(ekk, i, &self.sorted_work_data);
                self.work_data[self.work_count as usize] =
                    (i_col, (mv as f64) * ekk.info_.work_range_[i_col as usize]);
                self.work_count += 1;
            }
            // Look at all entries of final group to see what dual
            // infeasibilities might be created
            debug_assert!(((break_group + 1) as usize) < self.alt_work_group.len());
            let to_i = self.alt_work_group[(break_group + 1) as usize];
            debug_assert!((to_i as usize) <= self.sorted_work_data.len());
            //    HighsInt num_infeasibility = 0;
            let td: f64 = ekk.options_.dual_feasibility_tolerance;
            for i in self.alt_work_group[break_group as usize]..to_i {
                self.debug_report_bfrt_var(ekk, i, &self.sorted_work_data);
                let i_col = self.sorted_work_data[i as usize].0;
                let value = self.sorted_work_data[i as usize].1;
                let mv = ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt;
                let dual = ekk.info_.work_dual_[i_col as usize];
                let new_dual = dual - (move_out as f64) * (mv as f64) * self.work_theta * value;
                let new_dual_infeasibility = (mv as f64) * new_dual;
                let infeasible = new_dual_infeasibility < -td;
                if infeasible {
                    //	num_infeasibility++;
                    self.work_data[self.work_count as usize] =
                        (i_col, (mv as f64) * ekk.info_.work_range_[i_col as usize]);
                    self.work_count += 1;
                    debug_assert!(ekk.info_.work_range_[i_col as usize] < kHighsInf);
                }
            }
        }
        if self.work_theta == 0.0 {
            self.work_count = 0;
        }
        ekk.analysis_.simplex_timer_stop(Chuzc4dClock, 0);

        ekk.analysis_.simplex_timer_start(Chuzc4eClock, 0);
        // Sort workData by .first (iCol) so that columns of A are accessed in order?
        pdqsort(&mut self.work_data[..self.work_count as usize]);
        ekk.analysis_.simplex_timer_stop(Chuzc4eClock, 0);
        ekk.analysis_.simplex_timer_stop(Chuzc4Clock, 0);

        0
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::chooseFinalWorkGroupQuad (lines 317-373, 57 lines)
    pub fn choose_final_work_group_quad(&mut self, ekk: &mut HEkk) -> bool {
        let td = ekk.options_.dual_feasibility_tolerance;
        let full_count = self.work_count;
        self.work_count = 0;
        let mut total_change = kInitialTotalChange;
        let mut select_theta = self.work_theta;
        let total_delta = self.work_delta.abs();
        self.work_group.clear();
        self.work_group.push(0);
        let mut prev_work_count = self.work_count;
        let mut prev_remain_theta = kInitialRemainTheta;
        let mut prev_select_theta = select_theta;

        while select_theta < kMaxSelectTheta {
            let mut remain_theta = kInitialRemainTheta;
            for i in self.work_count..full_count {
                let i_col = self.work_data[i as usize].0;
                let value = self.work_data[i as usize].1;
                let dual = (ekk.basis_.nonbasic_move_[i_col as usize] as f64)
                    * ekk.info_.work_dual_[i_col as usize];
                // Tight satisfy
                if dual <= select_theta * value {
                    self.work_data.swap(self.work_count as usize, i as usize);
                    self.work_count += 1;
                    total_change += value * (ekk.info_.work_range_[i_col as usize]);
                } else if dual + td < remain_theta * value {
                    remain_theta = (dual + td) / value;
                }
            }
            self.work_group.push(self.work_count);

            // Update selectTheta with the value of remainTheta;
            select_theta = remain_theta;
            // Check for no change in this loop - to prevent infinite loop
            if (self.work_count == prev_work_count)
                && (prev_select_theta == select_theta)
                && (prev_remain_theta == remain_theta)
            {
                let num_var = ekk.lp_.num_col_ + ekk.lp_.num_row_;
                debug_dual_chuzc_fail_quad0(
                    &ekk.options_,
                    self.work_count,
                    &self.work_data,
                    num_var,
                    &ekk.info_.work_dual_,
                    select_theta,
                    remain_theta,
                    true,
                );
                return false;
            }
            // Record the initial values of workCount, remainTheta and selectTheta for
            // the next pass through the loop - to check for infinite loop condition
            prev_work_count = self.work_count;
            prev_remain_theta = remain_theta;
            prev_select_theta = select_theta;
            if total_change >= total_delta || self.work_count == full_count {
                break;
            }
        }
        // Check that at least one group has been identified
        if self.work_group.len() <= 1 {
            let num_var = ekk.lp_.num_col_ + ekk.lp_.num_row_;
            debug_dual_chuzc_fail_quad1(
                &ekk.options_,
                self.work_count,
                &self.work_data,
                num_var,
                &ekk.info_.work_dual_,
                select_theta,
                true,
            );
            return false;
        }
        true
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::quadChooseFinalWorkGroupQuad (lines 375-431, 57 lines)
    pub fn quad_choose_final_work_group_quad(&mut self, ekk: &mut HEkk) -> bool {
        let td: HighsCDouble = ekk.options_.dual_feasibility_tolerance.into();
        let full_count = self.work_count;
        self.work_count = 0;
        let mut total_change: HighsCDouble = HighsCDouble::from(kInitialTotalChange);
        let mut select_theta: HighsCDouble = HighsCDouble::from(self.work_theta);
        let total_delta: HighsCDouble = HighsCDouble::from(self.work_delta.abs());
        self.work_group.clear();
        self.work_group.push(0);
        let mut prev_work_count = self.work_count;
        let mut prev_remain_theta: HighsCDouble = HighsCDouble::from(kInitialRemainTheta);
        let mut prev_select_theta: HighsCDouble = select_theta;

        while select_theta < kMaxSelectTheta {
            let mut remain_theta: HighsCDouble = HighsCDouble::from(kInitialRemainTheta);
            for i in self.work_count..full_count {
                let i_col = self.work_data[i as usize].0;
                let value: HighsCDouble = HighsCDouble::from(self.work_data[i as usize].1);
                let dual: HighsCDouble = HighsCDouble::from(
                    (ekk.basis_.nonbasic_move_[i_col as usize] as f64)
                        * ekk.info_.work_dual_[i_col as usize],
                );
                // Tight satisfy
                if dual <= select_theta * value {
                    self.work_data.swap(self.work_count as usize, i as usize);
                    self.work_count += 1;
                    total_change += value * (ekk.info_.work_range_[i_col as usize]);
                } else if dual + td < remain_theta * value {
                    remain_theta = (dual + td) / value;
                }
            }
            self.work_group.push(self.work_count);

            // Update selectTheta with the value of remainTheta;
            select_theta = remain_theta;
            // Check for no change in this loop - to prevent infinite loop
            if (self.work_count == prev_work_count)
                && (prev_select_theta == select_theta)
                && (prev_remain_theta == remain_theta)
            {
                let num_var = ekk.lp_.num_col_ + ekk.lp_.num_row_;
                debug_dual_chuzc_fail_quad0(
                    &ekk.options_,
                    self.work_count,
                    &self.work_data,
                    num_var,
                    &ekk.info_.work_dual_,
                    f64::from(select_theta),
                    f64::from(remain_theta),
                    true,
                );
                return false;
            }
            // Record the initial values of workCount, remainTheta and selectTheta for
            // the next pass through the loop - to check for infinite loop condition
            prev_work_count = self.work_count;
            prev_remain_theta = remain_theta;
            prev_select_theta = select_theta;
            if total_change >= total_delta || self.work_count == full_count {
                break;
            }
        }
        // Check that at least one group has been identified
        if self.work_group.len() <= 1 {
            let num_var = ekk.lp_.num_col_ + ekk.lp_.num_row_;
            debug_dual_chuzc_fail_quad1(
                &ekk.options_,
                self.work_count,
                &self.work_data,
                num_var,
                &ekk.info_.work_dual_,
                f64::from(select_theta),
                true,
            );
            return false;
        }
        true
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::chooseFinalWorkGroupHeap (lines 433-493, 61 lines)
    pub fn choose_final_work_group_heap(&mut self, ekk: &mut HEkk) -> bool {
        let td: f64 = ekk.options_.dual_feasibility_tolerance;
        let full_count = self.alt_work_count;
        let mut total_change: f64 = kInitialTotalChange;
        let mut select_theta: f64 = self.work_theta;
        let total_delta = self.work_delta.abs();
        let mut heap_num_en: HighsInt = 0;
        let mut heap_i: Vec<HighsInt> = vec![0; (full_count + 1) as usize];
        let mut heap_v: Vec<f64> = vec![0.0; (full_count + 1) as usize];
        for i in 0..full_count {
            let i_col = self.original_work_data[i as usize].0;
            let value = self.original_work_data[i as usize].1;
            let dual = (ekk.basis_.nonbasic_move_[i_col as usize] as f64)
                * ekk.info_.work_dual_[i_col as usize];
            let ratio = dual / value;
            if ratio < kMaxSelectTheta {
                heap_num_en += 1;
                heap_i[heap_num_en as usize] = i;
                heap_v[heap_num_en as usize] = ratio;
            }
        }
        maxheapsort(&mut heap_v, &mut heap_i, heap_num_en);

        self.alt_work_count = 0;
        self.alt_work_group.clear();
        self.alt_work_group.push(self.alt_work_count);
        if heap_num_en <= 0 {
            let num_var = ekk.lp_.num_col_ + ekk.lp_.num_row_;
            // No entries in heap = > failure
            debug_dual_chuzc_fail_heap(
                &ekk.options_,
                self.alt_work_count,
                &self.original_work_data,
                num_var,
                &ekk.info_.work_dual_,
                select_theta,
                true,
            );
            return false;
        }
        let mut this_group_first_entry = self.alt_work_count;
        self.sorted_work_data.resize(heap_num_en as usize, (0, 0.0));
        for en in 1..=heap_num_en {
            let i = heap_i[en as usize];
            let i_col = self.original_work_data[i as usize].0;
            let value = self.original_work_data[i as usize].1;
            let dual = (ekk.basis_.nonbasic_move_[i_col as usize] as f64)
                * ekk.info_.work_dual_[i_col as usize];
            if dual > select_theta * value {
                // Breakpoint is in the next group, so record the pointer to its
                // first entry
                self.alt_work_group.push(self.alt_work_count);
                this_group_first_entry = self.alt_work_count;
                select_theta = (dual + td) / value;
                // End loop if all permitted groups have been identified
                if total_change >= total_delta {
                    break;
                }
            }
            // Store the breakpoint
            self.sorted_work_data[self.alt_work_count as usize].0 = i_col;
            self.sorted_work_data[self.alt_work_count as usize].1 = value;
            total_change += value * (ekk.info_.work_range_[i_col as usize]);
            self.alt_work_count += 1;
        }
        if self.alt_work_count > this_group_first_entry {
            self.alt_work_group.push(self.alt_work_count);
        }
        true
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::chooseFinalLargeAlpha (lines 495-529, 35 lines)
    /// MUST-PORT (determinism, random stream): ties broken by `ekk.info_.num_tot_permutation_`.
    pub fn choose_final_large_alpha(
        &mut self,
        ekk: &mut HEkk,
        break_index: &mut HighsInt,
        break_group: &mut HighsInt,
        pass_work_count: HighsInt,
        pass_work_data: &[(HighsInt, f64)],
        pass_work_group: &[HighsInt],
    ) {
        let mut final_compare: f64 = 0.0;
        for i in 0..pass_work_count {
            final_compare = cpp_max(final_compare, pass_work_data[i as usize].1);
        }
        final_compare = cpp_min(0.1 * final_compare, 1.0);
        let count_group = pass_work_group.len() as HighsInt - 1;
        *break_group = -1;
        *break_index = -1;
        let mut i_group = count_group - 1;
        while i_group >= 0 {
            let mut d_max_final: f64 = 0.0;
            let mut i_max_final: HighsInt = -1;
            for i in pass_work_group[i_group as usize]..pass_work_group[(i_group + 1) as usize] {
                if d_max_final < pass_work_data[i as usize].1 {
                    d_max_final = pass_work_data[i as usize].1;
                    i_max_final = i;
                } else if d_max_final == pass_work_data[i as usize].1 {
                    let j_col = pass_work_data[i_max_final as usize].0;
                    let i_col = pass_work_data[i as usize].0;
                    if ekk.info_.num_tot_permutation_[i_col as usize]
                        < ekk.info_.num_tot_permutation_[j_col as usize]
                    {
                        i_max_final = i;
                    }
                }
            }

            if pass_work_data[i_max_final as usize].1 > final_compare {
                *break_index = i_max_final;
                *break_group = i_group;
                break;
            }
            i_group -= 1;
        }
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::updateFlip (lines 531-546, 16 lines)
    pub fn update_flip(&mut self, ekk: &mut HEkk, bfrt_column: &mut HVector) {
        let mut dual_objective_value_change: f64 = 0.0;
        bfrt_column.clear();
        for i in 0..self.work_count {
            let i_col = self.work_data[i as usize].0;
            let change = self.work_data[i as usize].1;
            let mut local_dual_objective_change = change * ekk.info_.work_dual_[i_col as usize];
            local_dual_objective_change *= ekk.cost_scale_;
            dual_objective_value_change += local_dual_objective_change;
            ekk.flip_bound(i_col);
            ekk.lp_.a_matrix_.collect_aj(bfrt_column, i_col, change);
        }
        ekk.info_.updated_dual_objective_value += dual_objective_value_change;
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::updateDual (lines 548-566, 19 lines)
    pub fn update_dual(&mut self, ekk: &mut HEkk, theta: f64) {
        ekk.analysis_.simplex_timer_start(UpdateDualClock, 0);
        let mut dual_objective_value_change: f64 = 0.0;
        let pack_count = self.pack_count.max(0) as usize;
        for (&i_col, &pack_value) in self.pack_index[..pack_count]
            .iter()
            .zip(&self.pack_value[..pack_count])
        {
            // theta * pack_value is formed once; the C++ forms the same product twice.
            let delta_dual = theta * pack_value;
            ekk.info_.work_dual_[i_col as usize] -= delta_dual;
            // Identify the change to the dual objective
            let local_value = ekk.info_.work_value_[i_col as usize];
            let mut local_dual_objective_change =
                (ekk.basis_.nonbasic_flag_[i_col as usize] as f64) * (-local_value * delta_dual);
            local_dual_objective_change *= ekk.cost_scale_;
            dual_objective_value_change += local_dual_objective_change;
        }
        ekk.info_.updated_dual_objective_value += dual_objective_value_change;
        ekk.analysis_.simplex_timer_stop(UpdateDualClock, 0);
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::createFreelist (lines 568-578, 11 lines)
    pub fn create_freelist(&mut self, ekk: &mut HEkk) {
        self.free_list.clear();
        for i in 0..(ekk.lp_.num_col_ + ekk.lp_.num_row_) {
            if ekk.basis_.nonbasic_flag_[i as usize] != 0
                && highs_is_infinity(-ekk.info_.work_lower_[i as usize])
                && highs_is_infinity(ekk.info_.work_upper_[i as usize])
            {
                self.free_list.insert(i);
            }
        }
        //  debugFreeListNumEntries(ekk_instance_, freeList);
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::createFreemove (lines 580-598, 19 lines)
    pub fn create_freemove(&mut self, ekk: &mut HEkk, row_ep: &mut HVector) {
        // TODO: Check with Qi what this is doing and why it's expensive
        if !self.free_list.is_empty() {
            let ta: f64 = if ekk.info_.update_count < 10 {
                1e-9
            } else if ekk.info_.update_count < 20 {
                3e-8
            } else {
                1e-6
            };
            let move_out: HighsInt = if self.work_delta < 0.0 { -1 } else { 1 };
            for &i_var in self.free_list.iter() {
                debug_assert!(i_var < ekk.lp_.num_col_ + ekk.lp_.num_row_);
                let alpha = ekk.lp_.a_matrix_.compute_dot_hvector(row_ep, i_var);
                if alpha.abs() > ta {
                    if alpha * (move_out as f64) > 0.0 {
                        ekk.basis_.nonbasic_move_[i_var as usize] = 1;
                    } else {
                        ekk.basis_.nonbasic_move_[i_var as usize] = -1;
                    }
                }
            }
        }
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::deleteFreemove (lines 599-606, 8 lines)
    pub fn delete_freemove(&mut self, ekk: &mut HEkk) {
        if !self.free_list.is_empty() {
            for &i_var in self.free_list.iter() {
                debug_assert!(i_var < ekk.lp_.num_col_ + ekk.lp_.num_row_);
                ekk.basis_.nonbasic_move_[i_var as usize] = 0;
            }
        }
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::deleteFreelist (lines 608-612, 5 lines)
    pub fn delete_freelist(&mut self, ekk: &mut HEkk, i_var: HighsInt) {
        if !self.free_list.is_empty() {
            if self.free_list.contains(&i_var) {
                self.free_list.remove(&i_var);
            }
        }
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::computeDevexWeight (lines 614-632, 19 lines)
    /// C++ default arguments: slice = -1
    pub fn compute_devex_weight(&mut self, ekk: &mut HEkk, slice: HighsInt) {
        let rp_computed_edge_weight = false;
        self.computed_edge_weight = 0.0;
        for el_n in 0..self.pack_count {
            let vr_n = self.pack_index[el_n as usize];
            if ekk.basis_.nonbasic_flag_[vr_n as usize] == 0 {
                continue;
            }
            let pv =
                (ekk.info_.devex_index_[vr_n as usize] as f64) * self.pack_value[el_n as usize];
            if pv != 0.0 {
                self.computed_edge_weight += pv * pv;
            }
        }
        if rp_computed_edge_weight && slice >= 0 {
            println!(
                "HEkkDualRow::computeDevexWeight: Slice {}; computed_edge_weight = {:11.4e}",
                slice, self.computed_edge_weight
            );
        }
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::debugFindInWorkData (lines 634-640, 7 lines)
    /// Debug check: out of scope; stub returns 0 (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_find_in_work_data(
        &self,
        ekk: &HEkk,
        i_col: HighsInt,
        count: HighsInt,
        work_data_: &[(HighsInt, f64)],
    ) -> HighsInt {
        0
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::debugChooseColumnInfeasibilities (lines 642-673, 32 lines)
    /// Debug check: out of scope; stub returns 0 (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_choose_column_infeasibilities(&self, ekk: &HEkk) -> HighsInt {
        0
    }

    /// HiGHS: HEkkDualRow.cpp::HEkkDualRow::debugReportBfrtVar (lines 674-697, 24 lines)
    /// Debug check: out of scope; stub returns nothing (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_report_bfrt_var(
        &self,
        ekk: &HEkk,
        ix: HighsInt,
        pass_work_data: &[(HighsInt, f64)],
    ) {
    }
}
