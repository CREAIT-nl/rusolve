//! HiGHS source: util/HFactor.cpp (lines 871-1261), util/HFactor.h (inline helpers colInsert..rlinkDel)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Markowitz kernel factorization (pivot search with merit counts, threshold pivoting,
//! active-submatrix update) and its header-inline storage helpers.

use super::*;
use crate::prelude::*;

impl HFactor {
    /// HiGHS: HFactor.cpp::HFactor::buildKernel (lines 871-1261, 391 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    /// Time limit: `build_timer_->read()` -> `self.build_timer_.read()`.
    pub fn build_kernel(&mut self, basic_index: &mut [HighsInt]) -> HighsInt {
        // Deal with the kernel part by 'n-work' pivoting

        let mut fake_search: f64 = 0.0;
        let mut fake_fill: f64 = 0.0;
        let mut fake_eliminate: f64 = 0.0;

        let progress_report = false; // num_basic != num_row;
        let progress_frequency: HighsInt = 10000;
        // Initial timer frequency: may be reduced if iterations get slow
        let mut timer_frequency: HighsInt = 100;
        let mut previous_iteration_time: f64 = 0.0;
        let mut average_iteration_time: f64 = 0.0;
        let check_for_timeout = self.time_limit_ < kHighsInf;
        let mut search_k: HighsInt = 0;

        let check_nwork: HighsInt = -11;
        // `while (nwork-- > 0)`: nwork is decremented even when the test fails.
        loop {
            let nwork_before = self.nwork;
            self.nwork -= 1;
            if !(nwork_before > 0) {
                break;
            }
            if self.nwork == check_nwork {
                self.report_asm();
            }
            // Determine whether to return due to exceeding the time limit
            if check_for_timeout && search_k % timer_frequency == 0 {
                let current_time = self.build_timer_.read();
                let time_difference = current_time - previous_iteration_time;
                previous_iteration_time = current_time;
                let iteration_time = time_difference / (1.0 * timer_frequency as f64);
                average_iteration_time = 0.9 * average_iteration_time + 0.1 * iteration_time;

                if time_difference > self.time_limit_ / 1e3 {
                    timer_frequency = std::cmp::max(1 as HighsInt, timer_frequency / 10);
                }
                let iterations_left: HighsInt = self.kernel_dim - search_k + 1;
                let remaining_time_bound = average_iteration_time * iterations_left as f64;
                let total_time_bound = current_time + remaining_time_bound;
                if current_time > self.time_limit_ || total_time_bound > self.time_limit_ {
                    return kBuildKernelReturnTimeout;
                }
            }

            // 1. Search for the pivot
            let mut j_col_pivot: HighsInt = -1;
            let mut i_row_pivot: HighsInt = -1;
            // 1.1. Setup search merits
            let search_limit: HighsInt = std::cmp::min(self.nwork, 8 as HighsInt);
            let mut search_count: HighsInt = 0;

            let merit_limit: f64 = 1.0 * self.num_basic as f64 * self.num_row as f64;
            let mut merit_pivot: f64 = merit_limit;

            if progress_report && search_k != 0 {
                if search_k % progress_frequency == 0 {
                    let mut min_col_count: HighsInt = kHighsIInf;
                    let mut min_row_count: HighsInt = kHighsIInf;
                    for count in 1..self.num_row {
                        if self.col_link_first[count as usize] >= 0 {
                            min_col_count = count;
                            break;
                        }
                    }
                    for count in 1..self.num_basic {
                        if self.row_link_first[count as usize] >= 0 {
                            min_row_count = count;
                            break;
                        }
                    }
                    println!(
                        "HFactor::buildKernel stage = {:6}: min_col_count = {:3}; min_row_count = {:3}",
                        search_k, min_col_count, min_row_count
                    );
                }
            }
            search_k += 1;
            // 1.2. Search for local singletons
            let mut found_pivot = false;
            if !found_pivot && self.col_link_first[1] != -1 {
                // Not yet found a pivot and there is at least one column
                // singleton
                j_col_pivot = self.col_link_first[1];
                i_row_pivot = self.mc_index[self.mc_start[j_col_pivot as usize] as usize];
                found_pivot = true;
            }
            if !found_pivot && self.row_link_first[1] != -1 {
                i_row_pivot = self.row_link_first[1];
                j_col_pivot = self.mr_index[self.mr_start[i_row_pivot as usize] as usize];
                found_pivot = true;
            }
            let singleton_pivot = found_pivot;
            let mut candidate_pivot_value: f64 = 0.0;
            // 1.3. Major search loop
            //
            // Row count can be more than the number of rows if num_basic >
            // num_row
            let max_count: HighsInt = std::cmp::max(self.num_row, self.num_basic);
            let mut count: HighsInt = 2;
            while !found_pivot && count <= max_count {
                // Column count cannot exceed the number of rows
                if count <= self.num_row {
                    // 1.3.1 Search for columns
                    let mut j = self.col_link_first[count as usize];
                    while j != -1 {
                        let min_pivot = self.mc_min_pivot[j as usize];
                        let start = self.mc_start[j as usize];
                        let end = start + self.mc_count_a[j as usize];
                        for k in start..end {
                            if self.mc_value[k as usize].abs() >= min_pivot {
                                let i = self.mc_index[k as usize];
                                let row_count = self.mr_count[i as usize];
                                let merit_local: f64 =
                                    1.0 * (count - 1) as f64 * (row_count - 1) as f64;
                                if merit_pivot > merit_local {
                                    candidate_pivot_value = self.mc_value[k as usize].abs();
                                    merit_pivot = merit_local;
                                    j_col_pivot = j;
                                    i_row_pivot = i;
                                    found_pivot = found_pivot || (row_count < count);
                                }
                            }
                        }

                        let sc = search_count;
                        search_count += 1;
                        if sc >= search_limit && merit_pivot < merit_limit {
                            found_pivot = true;
                        }
                        if found_pivot {
                            break;
                        }

                        fake_search += count as f64;
                        j = self.col_link_next[j as usize];
                    }
                }

                // Row count cannot exceed the number of basic variables
                if count <= self.num_basic {
                    // 1.3.2 Search for rows
                    let mut i = self.row_link_first[count as usize];
                    while i != -1 {
                        let start = self.mr_start[i as usize];
                        let end = start + self.mr_count[i as usize];
                        for k in start..end {
                            let j = self.mr_index[k as usize];
                            let column_count = self.mc_count_a[j as usize];
                            let merit_local: f64 =
                                1.0 * (count - 1) as f64 * (column_count - 1) as f64;
                            if merit_local < merit_pivot {
                                let mut ifind = self.mc_start[j as usize];
                                while self.mc_index[ifind as usize] != i {
                                    ifind += 1;
                                }
                                if self.mc_value[ifind as usize].abs()
                                    >= self.mc_min_pivot[j as usize]
                                {
                                    candidate_pivot_value = self.mc_value[ifind as usize].abs();
                                    merit_pivot = merit_local;
                                    j_col_pivot = j;
                                    i_row_pivot = i;
                                    found_pivot = found_pivot || (column_count <= count);
                                }
                            }
                        }
                        let sc = search_count;
                        search_count += 1;
                        if sc >= search_limit && merit_pivot < merit_limit {
                            found_pivot = true;
                        }
                        if found_pivot {
                            break;
                        }
                        i = self.row_link_next[i as usize];
                    }

                    fake_search += count as f64;
                }
                count += 1;
            }
            // 1.4. If we found nothing: tell singular
            if i_row_pivot < 0 {
                // To detect the absence of a pivot, it should be sufficient
                // that iRowPivot is (still) -1, but add sanity asserts that
                // jColPivot is (still) -1 and foundPivot is false
                debug_assert!(j_col_pivot < 0);
                debug_assert!(!found_pivot);
                self.rank_deficiency = self.nwork + 1;
                crate::highs_log_dev!(
                    &self.log_options,
                    HighsLogType::kWarning,
                    "Factorization identifies rank deficiency of {}\n",
                    self.rank_deficiency
                );
                return self.rank_deficiency;
            }

            // 2. Elimination other elements by the pivot
            let original_pivotal_row_count = self.mr_count[i_row_pivot as usize];
            let original_pivotal_col_count = self.mc_count_a[j_col_pivot as usize];
            // 2.1. Delete the pivot
            //
            // Remove the pivot row index from the pivotal column of the
            // col-wise matrix. Also decreases the column count
            let pivot_multiplier: f64 = self.col_delete(j_col_pivot, i_row_pivot);
            // Remove the pivot column index from the pivotal row of the
            // row-wise matrix. Also decreases the row count
            self.row_delete(j_col_pivot, i_row_pivot);
            // Remove the pivotal column from the linked list of columns
            // containing it
            self.clink_del(j_col_pivot);
            // Remove the pivotal row from the linked list of rows containing
            // it
            self.rlink_del(i_row_pivot);
            if !singleton_pivot {
                debug_assert!(candidate_pivot_value == pivot_multiplier.abs());
            }
            if pivot_multiplier.abs() < self.pivot_tolerance {
                crate::highs_log_dev!(
                    &self.log_options,
                    HighsLogType::kWarning,
                    "Defer singular pivot = {:11.4}\n",
                    pivot_multiplier
                );
                // Matrix is singular, but defer return since other valid pivots
                // may exist.
                debug_assert!(
                    self.mr_count[i_row_pivot as usize] == original_pivotal_row_count - 1
                );
                if self.mr_count[i_row_pivot as usize] == 0 {
                    // The pivot corresponds to a singleton row. Entry is zeroed,
                    // and do no more since there may be other valid entries in
                    // the pivotal column
                    //
                    // Add the pivotal column to the linked list of columns with
                    // its new count
                    debug_assert!(
                        self.mc_count_a[j_col_pivot as usize] == original_pivotal_col_count - 1
                    );
                    let cnt = self.mc_count_a[j_col_pivot as usize];
                    self.clink_add(j_col_pivot, cnt);
                } else {
                    // Otherwise, other entries in the pivotal column will be
                    // smaller than the pivot, so zero the column
                    // HiGHS 1.16 #3354 (h16_hfactor_pivot): HFactor.cpp line 1194 (1.16.0) passes pivot_multiplier to zeroCol
                    self.zero_col(j_col_pivot, pivot_multiplier);
                    // Add the pivotal row to the linked list of rows with its new
                    // count
                    debug_assert!(
                        self.mr_count[i_row_pivot as usize] == original_pivotal_row_count - 1
                    );
                    let cnt = self.mr_count[i_row_pivot as usize];
                    self.rlink_add(i_row_pivot, cnt);
                }
                // No pivot found, so have to increment nwork
                self.nwork += 1;
                continue;
            }
            self.permute[j_col_pivot as usize] = i_row_pivot;
            debug_assert!(self.mc_var[j_col_pivot as usize] == basic_index[j_col_pivot as usize]);

            self.refactor_info_.pivot_row.push(i_row_pivot);
            self.refactor_info_
                .pivot_var
                .push(basic_index[j_col_pivot as usize]);
            self.refactor_info_.pivot_type.push(kPivotMarkowitz);

            // 2.2. Store active pivot column to L
            let start_a: HighsInt = self.mc_start[j_col_pivot as usize];
            let end_a: HighsInt = start_a + self.mc_count_a[j_col_pivot as usize];
            let mut mwz_column_count: HighsInt = 0;
            for k in start_a..end_a {
                let i_row: HighsInt = self.mc_index[k as usize];
                let value: f64 = self.mc_value[k as usize] / pivot_multiplier;
                self.mwz_column_index[mwz_column_count as usize] = i_row;
                mwz_column_count += 1;
                self.mwz_column_array[i_row as usize] = value;
                self.mwz_column_mark[i_row as usize] = 1;
                self.l_index.push(i_row);
                self.l_value.push(value);
                self.mr_count_before[i_row as usize] = self.mr_count[i_row as usize];
                self.row_delete(j_col_pivot, i_row);
            }
            self.l_start.push(self.l_index.len() as HighsInt);
            fake_fill += (2 * self.mc_count_a[j_col_pivot as usize]) as f64;

            // 2.3. Store non active pivot column to U
            let end_n: HighsInt = start_a + self.mc_space[j_col_pivot as usize];
            let start_n: HighsInt = end_n - self.mc_count_n[j_col_pivot as usize];
            for i in start_n..end_n {
                self.u_index.push(self.mc_index[i as usize]);
                self.u_value.push(self.mc_value[i as usize]);
            }
            self.u_pivot_index.push(i_row_pivot);
            self.u_pivot_value.push(pivot_multiplier);
            self.u_start.push(self.u_index.len() as HighsInt);
            fake_fill += (end_n - start_n) as f64;

            // 2.4. Loop over pivot row to eliminate other column
            let row_start: HighsInt = self.mr_start[i_row_pivot as usize];
            let row_end: HighsInt = row_start + self.mr_count[i_row_pivot as usize];
            for row_k in row_start..row_end {
                // 2.4.1. My pointer
                let i_col: HighsInt = self.mr_index[row_k as usize];
                let my_count: HighsInt = self.mc_count_a[i_col as usize];
                let my_start: HighsInt = self.mc_start[i_col as usize];
                let my_end: HighsInt = my_start + my_count - 1;
                let my_pivot: f64 = self.col_delete(i_col, i_row_pivot);
                self.col_store_n(i_col, i_row_pivot, my_pivot);

                // 2.4.2. Elimination on the overlapping part
                let mut n_fillin: HighsInt = mwz_column_count;
                let mut n_cancel: HighsInt = 0;
                {
                    let range = my_start as usize..my_end as usize;
                    let idx = &self.mc_index[range.clone()];
                    let val = &mut self.mc_value[range];
                    let mark = &mut self.mwz_column_mark[..];
                    let arr = &self.mwz_column_array[..];
                    debug_assert!(mark.len() == arr.len());
                    for (&i_row, v) in idx.iter().zip(val.iter_mut()) {
                        let r = i_row as usize;
                        debug_assert!(r < mark.len());
                        // SAFETY: i_row is a row index stored in the column copy, so
                        // 0 <= i_row < num_row; mwz_column_mark / mwz_column_array both
                        // have length num_row (setup.rs resize), so the access is in bounds.
                        unsafe {
                            if *mark.get_unchecked(r) != 0 {
                                *mark.get_unchecked_mut(r) = 0;
                                n_fillin -= 1;
                                let mut value: f64 = *v;
                                value -= my_pivot * *arr.get_unchecked(r);
                                if value.abs() < kHighsTiny {
                                    value = 0.0;
                                    n_cancel += 1;
                                }
                                *v = value;
                            }
                        }
                    }
                }
                fake_eliminate += mwz_column_count as f64;
                fake_eliminate += (n_fillin * 2) as f64;

                // 2.4.3. Remove cancellation gaps
                if n_cancel > 0 {
                    let mut new_end: HighsInt = my_start;
                    for my_k in my_start..my_end {
                        if self.mc_value[my_k as usize] != 0.0 {
                            self.mc_index[new_end as usize] = self.mc_index[my_k as usize];
                            self.mc_value[new_end as usize] = self.mc_value[my_k as usize];
                            new_end += 1;
                        } else {
                            let r = self.mc_index[my_k as usize];
                            self.row_delete(i_col, r);
                        }
                    }
                    self.mc_count_a[i_col as usize] = new_end - my_start;
                }

                // 2.4.4. Insert fill-in
                if n_fillin > 0 {
                    // 2.4.4.1 Check column size
                    if self.mc_count_a[i_col as usize] + self.mc_count_n[i_col as usize] + n_fillin
                        > self.mc_space[i_col as usize]
                    {
                        // p1&2=active, p3&4=non active, p5=new p1, p7=new p3
                        let p1: HighsInt = self.mc_start[i_col as usize];
                        let p2: HighsInt = p1 + self.mc_count_a[i_col as usize];
                        let p3: HighsInt =
                            p1 + self.mc_space[i_col as usize] - self.mc_count_n[i_col as usize];
                        let p4: HighsInt = p1 + self.mc_space[i_col as usize];
                        self.mc_space[i_col as usize] +=
                            std::cmp::max(self.mc_space[i_col as usize], n_fillin);
                        let p5: HighsInt = self.mc_index.len() as HighsInt;
                        self.mc_start[i_col as usize] = p5;
                        let p7: HighsInt =
                            p5 + self.mc_space[i_col as usize] - self.mc_count_n[i_col as usize];
                        let new_len = (p5 + self.mc_space[i_col as usize]) as usize;
                        self.mc_index.resize(new_len, 0);
                        self.mc_value.resize(new_len, 0.0);
                        self.mc_index
                            .copy_within(p1 as usize..p2 as usize, p5 as usize);
                        self.mc_value
                            .copy_within(p1 as usize..p2 as usize, p5 as usize);
                        self.mc_index
                            .copy_within(p3 as usize..p4 as usize, p7 as usize);
                        self.mc_value
                            .copy_within(p3 as usize..p4 as usize, p7 as usize);
                    }

                    // 2.4.4.2 Fill into column copy
                    {
                        // col_insert inlined on slices: appends at mc_start + mc_count_a.
                        let mut iput = (self.mc_start[i_col as usize]
                            + self.mc_count_a[i_col as usize])
                            as usize;
                        let mut added: HighsInt = 0;
                        let cols = &self.mwz_column_index[..mwz_column_count as usize];
                        let mark = &self.mwz_column_mark[..];
                        let arr = &self.mwz_column_array[..];
                        let mc_index = &mut self.mc_index[..];
                        let mc_value = &mut self.mc_value[..];
                        for &i_row in cols {
                            let r = i_row as usize;
                            debug_assert!(r < mark.len() && r < arr.len());
                            // SAFETY: i_row is a row index of the pivot column (< num_row);
                            // mark and arr have length num_row (setup.rs resize).
                            if unsafe { *mark.get_unchecked(r) } != 0 {
                                let v = -my_pivot * unsafe { *arr.get_unchecked(r) };
                                mc_index[iput] = i_row;
                                mc_value[iput] = v;
                                iput += 1;
                                added += 1;
                            }
                        }
                        self.mc_count_a[i_col as usize] += added;
                    }

                    // 2.4.4.3 Fill into the row copy
                    for i in 0..mwz_column_count {
                        let i_row: HighsInt = self.mwz_column_index[i as usize];
                        if self.mwz_column_mark[i_row as usize] != 0 {
                            // Expand row space
                            if self.mr_count[i_row as usize] == self.mr_space[i_row as usize] {
                                let p1: HighsInt = self.mr_start[i_row as usize];
                                let p2: HighsInt = p1 + self.mr_count[i_row as usize];
                                let p3: HighsInt = self.mr_index.len() as HighsInt;
                                self.mr_start[i_row as usize] = p3;
                                self.mr_space[i_row as usize] *= 2;
                                let new_len = (p3 + self.mr_space[i_row as usize]) as usize;
                                self.mr_index.resize(new_len, 0);
                                self.mr_index
                                    .copy_within(p1 as usize..p2 as usize, p3 as usize);
                            }
                            self.row_insert(i_col, i_row);
                        }
                    }
                }

                // 2.4.5. Reset pivot column mark
                {
                    let mark = &mut self.mwz_column_mark[..];
                    for &r in &self.mwz_column_index[..mwz_column_count as usize] {
                        debug_assert!((r as usize) < mark.len());
                        // SAFETY: r is a row index (< num_row == mark.len()).
                        unsafe { *mark.get_unchecked_mut(r as usize) = 1 };
                    }
                }

                // 2.4.6. Fix max value and link list
                self.col_fix_max(i_col);
                if my_count != self.mc_count_a[i_col as usize] {
                    self.clink_del(i_col);
                    let cnt = self.mc_count_a[i_col as usize];
                    self.clink_add(i_col, cnt);
                }
            }

            // 2.5. Clear pivot column buffer
            for i in 0..mwz_column_count {
                let r = self.mwz_column_index[i as usize];
                self.mwz_column_mark[r as usize] = 0;
            }

            // 2.6. Correct row links for the remain active part
            for i in start_a..end_a {
                let i_row: HighsInt = self.mc_index[i as usize];
                if self.mr_count_before[i_row as usize] != self.mr_count[i_row as usize] {
                    self.rlink_del(i_row);
                    let cnt = self.mr_count[i_row as usize];
                    self.rlink_add(i_row, cnt);
                }
            }
        }
        self.build_synthetic_tick += fake_search * 20.0 + fake_fill * 160.0 + fake_eliminate * 80.0;
        self.rank_deficiency = 0;
        self.rank_deficiency
    }

    /// HiGHS: HFactor.h::HFactor::colInsert (header inline)
    pub fn col_insert(&mut self, i_col: HighsInt, i_row: HighsInt, value: f64) {
        let iput = self.mc_start[i_col as usize] + self.mc_count_a[i_col as usize];
        self.mc_count_a[i_col as usize] += 1;
        self.mc_index[iput as usize] = i_row;
        self.mc_value[iput as usize] = value;
    }

    /// HiGHS: HFactor.h::HFactor::colStoreN (header inline)
    pub fn col_store_n(&mut self, i_col: HighsInt, i_row: HighsInt, value: f64) {
        self.mc_count_n[i_col as usize] += 1;
        let iput = self.mc_start[i_col as usize] + self.mc_space[i_col as usize]
            - self.mc_count_n[i_col as usize];
        self.mc_index[iput as usize] = i_row;
        self.mc_value[iput as usize] = value;
    }

    /// HiGHS: HFactor.h::HFactor::colFixMax (header inline)
    /// HiGHS 1.16 (#3354): HFactor.h::HFactor::colFixMax (lines 522-550)
    pub fn col_fix_max(&mut self, i_col: HighsInt) {
        let mut max_value: f64 = 0.0;
        let k = self.mc_start[i_col as usize] as usize;
        for &v in &self.mc_value[k..k + self.mc_count_a[i_col as usize] as usize] {
            max_value = cpp_max(max_value, v.abs());
        }
        // HiGHS 1.16 #3354 (h16_hfactor_pivot): HFactor.h colFixMax lines 522-550 (1.16.0)
        if self.h16_hfactor_pivot {
            self.mc_min_pivot[i_col as usize] =
                cpp_max(max_value * self.pivot_threshold, self.pivot_tolerance);
        } else {
            self.mc_min_pivot[i_col as usize] = max_value * self.pivot_threshold;
        }
    }

    /// HiGHS: HFactor.h::HFactor::colDelete (header inline)
    pub fn col_delete(&mut self, i_col: HighsInt, i_row: HighsInt) -> f64 {
        let mut idel = self.mc_start[i_col as usize];
        self.mc_count_a[i_col as usize] -= 1;
        let imov = idel + self.mc_count_a[i_col as usize];
        idel += self.mc_index[idel as usize..]
            .iter()
            .position(|&x| x == i_row)
            .expect("col_delete: row not in column") as HighsInt;
        let pivot_multiplier = self.mc_value[idel as usize];
        self.mc_index[idel as usize] = self.mc_index[imov as usize];
        self.mc_value[idel as usize] = self.mc_value[imov as usize];
        pivot_multiplier
    }

    /// HiGHS: HFactor.h::HFactor::rowInsert (header inline)
    pub fn row_insert(&mut self, i_col: HighsInt, i_row: HighsInt) {
        let iput = self.mr_start[i_row as usize] + self.mr_count[i_row as usize];
        self.mr_count[i_row as usize] += 1;
        self.mr_index[iput as usize] = i_col;
    }

    /// HiGHS: HFactor.h::HFactor::rowDelete (header inline)
    pub fn row_delete(&mut self, i_col: HighsInt, i_row: HighsInt) {
        let mut idel = self.mr_start[i_row as usize];
        self.mr_count[i_row as usize] -= 1;
        let imov = idel + self.mr_count[i_row as usize];
        idel += self.mr_index[idel as usize..]
            .iter()
            .position(|&x| x == i_col)
            .expect("row_delete: column not in row") as HighsInt;
        self.mr_index[idel as usize] = self.mr_index[imov as usize];
    }

    /// HiGHS: HFactor.h::HFactor::clinkAdd (header inline)
    pub fn clink_add(&mut self, index: HighsInt, count: HighsInt) {
        let mover = self.col_link_first[count as usize];
        self.col_link_last[index as usize] = -2 - count;
        self.col_link_next[index as usize] = mover;
        self.col_link_first[count as usize] = index;
        if mover >= 0 {
            self.col_link_last[mover as usize] = index;
        }
    }

    /// HiGHS: HFactor.h::HFactor::clinkDel (header inline)
    pub fn clink_del(&mut self, index: HighsInt) {
        let xlast = self.col_link_last[index as usize];
        let xnext = self.col_link_next[index as usize];
        if xlast >= 0 {
            self.col_link_next[xlast as usize] = xnext;
        } else {
            self.col_link_first[(-xlast - 2) as usize] = xnext;
        }
        if xnext >= 0 {
            self.col_link_last[xnext as usize] = xlast;
        }
    }

    /// HiGHS: HFactor.h::HFactor::rlinkAdd (header inline)
    pub fn rlink_add(&mut self, index: HighsInt, count: HighsInt) {
        let mover = self.row_link_first[count as usize];
        self.row_link_last[index as usize] = -2 - count;
        self.row_link_next[index as usize] = mover;
        self.row_link_first[count as usize] = index;
        if mover >= 0 {
            self.row_link_last[mover as usize] = index;
        }
    }

    /// HiGHS: HFactor.h::HFactor::rlinkDel (header inline)
    pub fn rlink_del(&mut self, index: HighsInt) {
        let xlast = self.row_link_last[index as usize];
        let xnext = self.row_link_next[index as usize];
        if xlast >= 0 {
            self.row_link_next[xlast as usize] = xnext;
        } else {
            self.row_link_first[(-xlast - 2) as usize] = xnext;
        }
        if xnext >= 0 {
            self.row_link_last[xnext as usize] = xlast;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// HiGHS 1.16 (#3354, HFactor.h colFixMax lines 522-550): the column of
    /// the TestFactor case with entries -8.1e-11 and -1.6e-10. With the flag on,
    /// mc_min_pivot is max(1.6e-11, pivot_tolerance) = 1e-10. With it off, 1.15.1
    /// keeps max_value * pivot_threshold = 1.6e-11.
    #[test]
    fn h16_hfactor_pivot_col_fix_max_floor() {
        let mut factor = HFactor::default();
        factor.pivot_threshold = 0.1;
        factor.pivot_tolerance = 1e-10;
        factor.mc_start = vec![0];
        factor.mc_count_a = vec![2];
        factor.mc_value = vec![-8.1e-11, -1.6e-10];
        factor.mc_min_pivot = vec![0.0];

        factor.h16_hfactor_pivot = false;
        factor.col_fix_max(0);
        assert!((factor.mc_min_pivot[0] - 1.6e-11).abs() < 1e-20);

        factor.h16_hfactor_pivot = true;
        factor.col_fix_max(0);
        assert_eq!(factor.mc_min_pivot[0], 1e-10);
    }

    /// HiGHS 1.16 (#3354, check/TestFactor.cpp Factor-zero-col-pivot-tolerance):
    /// the near-dependent LP solved with simplex and presolve off.
    fn h16_hfactor_pivot_tolerance_run(flag: bool) -> (HighsStatus, HighsModelStatus) {
        let mut lp = HighsLp {
            num_col_: 4,
            num_row_: 5,
            col_cost_: vec![0.0; 4],
            col_lower_: vec![0.0, 0.0, 0.0, -kHighsInf],
            col_upper_: vec![kHighsInf; 4],
            row_lower_: vec![-kHighsInf; 5],
            row_upper_: vec![0.0, 0.0, 0.0, 0.0, -1.0],
            ..Default::default()
        };
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 5;
        lp.a_matrix_.start_ = vec![0, 4, 8, 11, 14];
        lp.a_matrix_.index_ = vec![0, 2, 3, 4, 0, 2, 3, 4, 0, 1, 4, 1, 2, 4];
        lp.a_matrix_.value_ = vec![
            2.0, -1.0, 3e6, -1.0, -1e13, 1.0, -3000001.0, 1.0, 1.0, 1.0, -2.0, 1e7, -1e7, 1.0,
        ];

        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("presolve", "off"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", "simplex"),
            HighsStatus::kOk
        );
        assert_eq!(highs.set_option_value_int("threads", 1), HighsStatus::kOk);
        assert_eq!(
            highs.set_option_value_bool("h16_hfactor_pivot", flag),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        let run_status = highs.run();
        (run_status, highs.get_model_status())
    }

    /// The 1.16.0 result is kOptimal (TestFactor.cpp REQUIRE).
    #[test]
    fn h16_hfactor_pivot_tolerance_lp() {
        assert_eq!(
            h16_hfactor_pivot_tolerance_run(true),
            (HighsStatus::kOk, HighsModelStatus::kOptimal)
        );
    }

    /// With the flag off the 1.15.1 result is kept: HiGHS 1.15.1 (release highspy) returns
    /// kError with model status kSolveError. A debug 1.15.1 build instead fails
    /// `assert(abs_value < pivot_tolerance)` (1.15.1 HFactor.cpp:1520), which is the bug #3354
    /// fixes, so the release result is only observable with debug assertions off.
    #[test]
    #[cfg(not(debug_assertions))]
    fn h16_hfactor_pivot_tolerance_lp_off() {
        assert_eq!(
            h16_hfactor_pivot_tolerance_run(false),
            (HighsStatus::kError, HighsModelStatus::kSolveError)
        );
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "abs_value < self.pivot_tolerance")]
    fn h16_hfactor_pivot_tolerance_lp_off_trips_cpp_assert() {
        h16_hfactor_pivot_tolerance_run(false);
    }
}
