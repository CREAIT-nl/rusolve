//! HiGHS source: util/HFactor.cpp (lines 2052-2597)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Basis-change updates of the factored form and InvertibleRepresentation copy-in/out.

use super::*;
use crate::prelude::*;

impl HFactor {
    /// HiGHS: HFactor.cpp::HFactor::updateCFT (lines 2052-2298, 247 lines)
    /// Collective Forrest-Tomlin update over the `aq->next` / `ep->next` chains built by
    /// `HEkkDual::majorUpdateFactor` (PAMI). `i_row` is the C++ `HighsInt* iRow` array, one
    /// entry per chained vector.
    ///
    /// # Safety
    /// `aq` and `ep` point to live `HVector`s; every vector reachable through their `next`
    /// links is live, the two chains have the same length, no vector appears twice, and
    /// none of them is otherwise borrowed for the duration of the call.
    pub unsafe fn update_cft(
        &mut self,
        aq: *mut HVector,
        ep: *mut HVector,
        i_row: &mut [HighsInt],
    ) {
        // HiGHS: HFactor.cpp::HFactor::updateCFT (lines 2052-2298)

        // Count the number of updates in the chain
        let mut num_update: HighsInt = 0;
        let mut vec_ptr = aq;
        while !vec_ptr.is_null() {
            num_update += 1;
            // SAFETY: caller guarantees the chain is live and not otherwise borrowed
            let vec_ref = unsafe { &*vec_ptr };
            vec_ptr = match &vec_ref.next {
                Some(next) => next.0.as_ptr(),
                None => std::ptr::null_mut(),
            };
        }

        // Build the aq_work and ep_work arrays
        let mut aq_work: Vec<*mut HVector> = Vec::with_capacity(num_update as usize);
        let mut ep_work: Vec<*mut HVector> = Vec::with_capacity(num_update as usize);

        let mut aq_ptr = aq;
        let mut ep_ptr = ep;
        for _ in 0..num_update {
            aq_work.push(aq_ptr);
            ep_work.push(ep_ptr);
            // SAFETY: caller guarantees the chain is live and not otherwise borrowed
            let aq_ref = unsafe { &*aq_ptr };
            aq_ptr = match &aq_ref.next {
                Some(next) => next.0.as_ptr(),
                None => std::ptr::null_mut(),
            };
            // SAFETY: caller guarantees the chain is live and not otherwise borrowed
            let ep_ref = unsafe { &*ep_ptr };
            ep_ptr = match &ep_ref.next {
                Some(next) => next.0.as_ptr(),
                None => std::ptr::null_mut(),
            };
        }

        // Pivot related buffers
        let pf_np0 = self.pf_pivot_index.len() as HighsInt;
        let mut p_logic: Vec<HighsInt> = vec![0; num_update as usize];
        let mut p_value: Vec<f64> = vec![0.0; num_update as usize];
        let mut p_alpha: Vec<f64> = vec![0.0; num_update as usize];

        for cp in 0..num_update {
            let c_row = i_row[cp as usize];
            let i_logic = self.u_pivot_lookup[c_row as usize];
            p_logic[cp as usize] = i_logic;
            p_value[cp as usize] = self.u_pivot_value[i_logic as usize];
            // SAFETY: aq_work[cp] is a live pointer from the chain
            let aq_ref = unsafe { &*aq_work[cp as usize] };
            p_alpha[cp as usize] = aq_ref.array[c_row as usize];
        }

        // Temporary U pointers
        let mut t_start: Vec<HighsInt> = vec![0; (num_update + 1) as usize];
        let mut t_pivot: Vec<f64> = vec![0.0; num_update as usize];
        t_start[0] = self.u_index.len() as HighsInt;

        // Logically sorted previous row_ep
        let mut sorted_pp: Vec<(HighsInt, HighsInt)> = Vec::new();

        // Major update loop
        for cp in 0..num_update {
            // 1. Expand partial FTRAN result to buffer
            self.iwork.clear();
            // SAFETY: aq_work[cp] is a live pointer from the chain
            let aq_ref = unsafe { &*aq_work[cp as usize] };
            for i in 0..aq_ref.pack_count {
                let index = aq_ref.pack_index[i as usize];
                let value = aq_ref.pack_value[i as usize];
                self.iwork.push(index);
                self.dwork[index as usize] = value;
            }

            // 2. Update partial FTRAN result by recent FT matrix
            for pp in 0..cp {
                let p_row = i_row[pp as usize];
                let mut value = self.dwork[p_row as usize];
                let pf_pp = pp + pf_np0;
                let pf_start_lo = self.pf_start[pf_pp as usize];
                let pf_start_hi = self.pf_start[(pf_pp + 1) as usize];
                for i in pf_start_lo..pf_start_hi {
                    let idx = self.pf_index[i as usize];
                    let val = self.pf_value[i as usize];
                    value -= self.dwork[idx as usize] * val;
                }
                self.iwork.push(p_row); // OK to duplicate
                self.dwork[p_row as usize] = value;
            }

            // 3. Store the partial FTRAN result to matrix U
            let ppaq = self.dwork[i_row[cp as usize] as usize]; // pivot of the partial aq
            self.dwork[i_row[cp as usize] as usize] = 0.0;
            let u_count_x = t_start[cp as usize];
            let u_start_x = u_count_x;
            for &index in &self.iwork {
                let value = self.dwork[index as usize];
                self.dwork[index as usize] = 0.0; // This effectively removes all duplication
                if value.abs() > kHighsTiny {
                    self.u_index.push(index);
                    self.u_value.push(value);
                }
            }
            let u_count_x = self.u_index.len() as HighsInt;
            t_start[(cp + 1) as usize] = u_count_x;
            t_pivot[cp as usize] = p_value[cp as usize] * p_alpha[cp as usize];

            // 4. Expand partial BTRAN result to buffer
            self.iwork.clear();
            // SAFETY: ep_work[cp] is a live pointer from the chain
            let ep_ref = unsafe { &*ep_work[cp as usize] };
            for i in 0..ep_ref.pack_count {
                let index = ep_ref.pack_index[i as usize];
                let value = ep_ref.pack_value[i as usize];
                self.iwork.push(index);
                self.dwork[index as usize] = value;
            }

            // 5. Delete logical later rows (in logical order)
            for isort in 0..cp {
                let pp = sorted_pp[isort as usize].1;
                let p_row = i_row[pp as usize];
                let multiplier = -p_value[pp as usize] * self.dwork[p_row as usize];
                if self.dwork[p_row as usize].abs() > kHighsTiny {
                    // SAFETY: ep_work[pp] is a live pointer from the chain
                    let ep_pp_ref = unsafe { &*ep_work[pp as usize] };
                    for i in 0..ep_pp_ref.pack_count {
                        let index = ep_pp_ref.pack_index[i as usize];
                        let value = ep_pp_ref.pack_value[i as usize];
                        self.iwork.push(index);
                        self.dwork[index as usize] += value * multiplier;
                    }
                }
                self.dwork[p_row as usize] = 0.0; // Force to be 0
            }

            // 6. Update partial BTRAN result by recent U columns
            for pp in 0..cp {
                let kpivot = i_row[pp as usize];
                let mut value = self.dwork[kpivot as usize];
                let t_start_lo = t_start[pp as usize];
                let t_start_hi = t_start[(pp + 1) as usize];
                for k in t_start_lo..t_start_hi {
                    let idx = self.u_index[k as usize];
                    let val = self.u_value[k as usize];
                    value -= self.dwork[idx as usize] * val;
                }
                value /= t_pivot[pp as usize];
                self.iwork.push(kpivot);
                self.dwork[kpivot as usize] = value; // Again OK to duplicate
            }

            // 6.x compute current alpha
            let mut thex: f64 = 0.0;
            for k in u_start_x..u_count_x {
                let index = self.u_index[k as usize];
                let value = self.u_value[k as usize];
                thex += self.dwork[index as usize] * value;
            }
            t_pivot[cp as usize] = ppaq + thex * p_value[cp as usize];

            // 7. Store BTRAN result to FT elimination, update logic helper
            self.dwork[i_row[cp as usize] as usize] = 0.0;
            let pivot_multiplier = -p_value[cp as usize];
            for &index in &self.iwork {
                let value = self.dwork[index as usize];
                self.dwork[index as usize] = 0.0;
                if value.abs() > kHighsTiny {
                    self.pf_index.push(index);
                    self.pf_value.push(value * pivot_multiplier);
                }
            }
            self.pf_pivot_index.push(i_row[cp as usize]);
            self.u_total_x += self.pf_index.len() as HighsInt - *self.pf_start.last().unwrap();
            self.pf_start.push(self.pf_index.len() as HighsInt);

            // 8. Update the sorted ep
            sorted_pp.push((p_logic[cp as usize], cp));
            crate::util::highs_sort::pdqsort(&mut sorted_pp);
        }

        // Now modify the U matrix
        for cp in 0..num_update {
            // 1. Delete pivotal row from U
            let c_index = i_row[cp as usize];
            let c_logic = p_logic[cp as usize];
            self.u_total_x -= self.ur_lastp[c_logic as usize] - self.ur_start[c_logic as usize];
            let mut k = self.ur_start[c_logic as usize];
            while k < self.ur_lastp[c_logic as usize] {
                // Find the pivotal position
                let i_logic = self.u_pivot_lookup[self.ur_index[k as usize] as usize];
                let mut i_find = self.u_start[i_logic as usize];
                self.u_last_p[i_logic as usize] -= 1;
                let i_last = self.u_last_p[i_logic as usize];
                while i_find <= i_last {
                    if self.u_index[i_find as usize] == c_index {
                        break;
                    }
                    i_find += 1;
                }
                // Put last to find, and delete last
                self.u_index[i_find as usize] = self.u_index[i_last as usize];
                self.u_value[i_find as usize] = self.u_value[i_last as usize];
                k += 1;
            }

            // 2. Delete pivotal column from UR
            self.u_total_x -= self.u_last_p[c_logic as usize] - self.u_start[c_logic as usize];
            let mut k = self.u_start[c_logic as usize];
            while k < self.u_last_p[c_logic as usize] {
                // Find the pivotal position
                let i_logic = self.u_pivot_lookup[self.u_index[k as usize] as usize];
                let mut i_find = self.ur_start[i_logic as usize];
                self.ur_lastp[i_logic as usize] -= 1;
                let i_last = self.ur_lastp[i_logic as usize];
                while i_find <= i_last {
                    if self.ur_index[i_find as usize] == c_index {
                        break;
                    }
                    i_find += 1;
                }
                // Put last to find, and delete last
                self.ur_space[i_logic as usize] += 1;
                self.ur_index[i_find as usize] = self.ur_index[i_last as usize];
                self.ur_value[i_find as usize] = self.ur_value[i_last as usize];
                k += 1;
            }

            // 3. Insert the (stored) partial FTRAN to the row matrix
            let u_start_x = t_start[cp as usize];
            let u_end_x = t_start[(cp + 1) as usize];
            self.u_total_x += u_end_x - u_start_x;
            // Store column as UR elements
            for k in u_start_x..u_end_x {
                // Which ETA file
                let i_logic = self.u_pivot_lookup[self.u_index[k as usize] as usize];

                // Move row to the end if necessary
                if self.ur_space[i_logic as usize] == 0 {
                    // Make pointers
                    let row_start = self.ur_start[i_logic as usize];
                    let row_count = self.ur_lastp[i_logic as usize] - row_start;
                    let new_start = self.ur_index.len() as HighsInt;
                    let new_space = (row_count as f64 * 1.1 + 5.0) as HighsInt;

                    // Check matrix UR
                    self.ur_index.resize((new_start + new_space) as usize, 0);
                    self.ur_value.resize((new_start + new_space) as usize, 0.0);

                    // Move elements
                    let i_from = row_start;
                    let i_end = row_start + row_count;
                    let i_to = new_start;
                    self.ur_index
                        .copy_within(i_from as usize..i_end as usize, i_to as usize);
                    self.ur_value
                        .copy_within(i_from as usize..i_end as usize, i_to as usize);

                    // Save new pointers
                    self.ur_start[i_logic as usize] = new_start;
                    self.ur_lastp[i_logic as usize] = new_start + row_count;
                    self.ur_space[i_logic as usize] = new_space - row_count;
                }

                // Put into the next available space
                self.ur_space[i_logic as usize] -= 1;
                let i_put = self.ur_lastp[i_logic as usize];
                self.ur_lastp[i_logic as usize] += 1;
                self.ur_index[i_put as usize] = c_index;
                self.ur_value[i_put as usize] = self.u_value[k as usize];
            }

            // 4. Save pointers
            self.u_start.push(u_start_x);
            self.u_last_p.push(u_end_x);

            self.ur_start.push(self.ur_start[c_logic as usize]);
            self.ur_lastp.push(self.ur_start[c_logic as usize]);
            self.ur_space.push(
                self.ur_space[c_logic as usize] + self.ur_lastp[c_logic as usize]
                    - self.ur_start[c_logic as usize],
            );

            self.u_pivot_lookup[c_index as usize] = self.u_pivot_index.len() as HighsInt;
            self.u_pivot_index[c_logic as usize] = -1;
            self.u_pivot_index.push(c_index);
            self.u_pivot_value.push(t_pivot[cp as usize]);
        }

        // C++ deletes the raw arrays here; Rust Vecs drop automatically.
    }

    /// HiGHS: HFactor.cpp::HFactor::updateFT (lines 2300-2415, 116 lines)
    pub fn update_ft(&mut self, aq: &mut HVector, ep: &mut HVector, i_row: HighsInt) {
        // Store pivot
        let p_logic = self.u_pivot_lookup[i_row as usize];
        let pivot = self.u_pivot_value[p_logic as usize];
        let alpha = aq.array[i_row as usize];
        self.u_pivot_index[p_logic as usize] = -1;

        // Delete pivotal row from U
        let mut k = self.ur_start[p_logic as usize];
        while k < self.ur_lastp[p_logic as usize] {
            // Find the pivotal position
            let i_logic = self.u_pivot_lookup[self.ur_index[k as usize] as usize];
            let mut i_find = self.u_start[i_logic as usize];
            self.u_last_p[i_logic as usize] -= 1;
            let i_last = self.u_last_p[i_logic as usize];
            // PORT NOTE (perf): linear search as a slice scan; i_find ends at the first
            // match, or at i_last + 1 when none matches (as the C++ loop).
            if i_find <= i_last {
                i_find += self.u_index[i_find as usize..=i_last as usize]
                    .iter()
                    .position(|&x| x == i_row)
                    .unwrap_or((i_last - i_find + 1) as usize)
                    as HighsInt;
            }
            // Put last to find, and delete last
            self.u_index[i_find as usize] = self.u_index[i_last as usize];
            self.u_value[i_find as usize] = self.u_value[i_last as usize];
            k += 1;
        }

        // Delete pivotal column from UR
        let mut k = self.u_start[p_logic as usize];
        while k < self.u_last_p[p_logic as usize] {
            // Find the pivotal position
            let i_logic = self.u_pivot_lookup[self.u_index[k as usize] as usize];
            let mut i_find = self.ur_start[i_logic as usize];
            self.ur_lastp[i_logic as usize] -= 1;
            let i_last = self.ur_lastp[i_logic as usize];
            // PORT NOTE (perf): linear search as a slice scan; i_find ends at the first
            // match, or at i_last + 1 when none matches (as the C++ loop).
            if i_find <= i_last {
                i_find += self.ur_index[i_find as usize..=i_last as usize]
                    .iter()
                    .position(|&x| x == i_row)
                    .unwrap_or((i_last - i_find + 1) as usize)
                    as HighsInt;
            }
            // Put last to find, and delete last
            self.ur_space[i_logic as usize] += 1;
            self.ur_index[i_find as usize] = self.ur_index[i_last as usize];
            self.ur_value[i_find as usize] = self.ur_value[i_last as usize];
            k += 1;
        }

        // Store column to U
        self.u_start.push(self.u_index.len() as HighsInt);
        for i in 0..aq.pack_count {
            if aq.pack_index[i as usize] != i_row {
                self.u_index.push(aq.pack_index[i as usize]);
                self.u_value.push(aq.pack_value[i as usize]);
            }
        }
        self.u_last_p.push(self.u_index.len() as HighsInt);
        let u_start_x = *self.u_start.last().unwrap();
        let u_end_x = *self.u_last_p.last().unwrap();
        self.u_total_x += u_end_x - u_start_x + 1;

        // Store column as UR elements
        for k in u_start_x..u_end_x {
            // Which ETA file
            let i_logic = self.u_pivot_lookup[self.u_index[k as usize] as usize];

            // Move row to the end if necessary
            if self.ur_space[i_logic as usize] == 0 {
                // Make pointers
                let row_start = self.ur_start[i_logic as usize];
                let row_count = self.ur_lastp[i_logic as usize] - row_start;
                let new_start = self.ur_index.len() as HighsInt;
                let new_space = (row_count as f64 * 1.1 + 5.0) as HighsInt;

                // Check matrix UR
                self.ur_index.resize((new_start + new_space) as usize, 0);
                self.ur_value.resize((new_start + new_space) as usize, 0.0);

                // Move elements
                let i_from = row_start;
                let i_end = row_start + row_count;
                let i_to = new_start;
                self.ur_index
                    .copy_within(i_from as usize..i_end as usize, i_to as usize);
                self.ur_value
                    .copy_within(i_from as usize..i_end as usize, i_to as usize);

                // Save new pointers
                self.ur_start[i_logic as usize] = new_start;
                self.ur_lastp[i_logic as usize] = new_start + row_count;
                self.ur_space[i_logic as usize] = new_space - row_count;
            }

            // Put into the next available space
            self.ur_space[i_logic as usize] -= 1;
            let i_put = self.ur_lastp[i_logic as usize];
            self.ur_lastp[i_logic as usize] += 1;
            self.ur_index[i_put as usize] = i_row;
            self.ur_value[i_put as usize] = self.u_value[k as usize];
        }

        // Store UR pointers
        self.ur_start.push(self.ur_start[p_logic as usize]);
        self.ur_lastp.push(self.ur_start[p_logic as usize]);
        self.ur_space.push(
            self.ur_space[p_logic as usize] + self.ur_lastp[p_logic as usize]
                - self.ur_start[p_logic as usize],
        );

        // Update pivot count
        self.u_pivot_lookup[i_row as usize] = self.u_pivot_index.len() as HighsInt;
        self.u_pivot_index.push(i_row);
        self.u_pivot_value.push(pivot * alpha);

        // Store row_ep as R matrix
        for i in 0..ep.pack_count {
            if ep.pack_index[i as usize] != i_row {
                self.pf_index.push(ep.pack_index[i as usize]);
                self.pf_value.push(-ep.pack_value[i as usize] * pivot);
            }
        }
        self.u_total_x += (self.pf_index.len() as HighsInt) - *self.pf_start.last().unwrap();

        // Store R matrix pivot
        self.pf_pivot_index.push(i_row);
        self.pf_start.push(self.pf_index.len() as HighsInt);

        // Update total countX
        self.u_total_x -= self.u_last_p[p_logic as usize] - self.u_start[p_logic as usize];
        self.u_total_x -= self.ur_lastp[p_logic as usize] - self.ur_start[p_logic as usize];
    }

    /// HiGHS: HFactor.cpp::HFactor::updatePF (lines 2417-2441, 25 lines)
    pub fn update_pf(&mut self, aq: &mut HVector, i_row: HighsInt, hint: &mut HighsInt) {
        // Check space
        let column_count = aq.pack_count;

        // Copy the pivotal column
        for i in 0..column_count {
            let index = aq.pack_index[i as usize];
            let value = aq.pack_value[i as usize];
            if index != i_row {
                self.pf_index.push(index);
                self.pf_value.push(value);
            }
        }

        // Save pivot
        self.pf_pivot_index.push(i_row);
        self.pf_pivot_value.push(aq.array[i_row as usize]);
        self.pf_start.push(self.pf_index.len() as HighsInt);

        // Check refactor
        self.u_total_x += aq.pack_count;
        if self.u_total_x > self.u_merit_x {
            *hint = 1;
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::updateMPF (lines 2443-2473, 31 lines)
    pub fn update_mpf(
        &mut self,
        aq: &mut HVector,
        ep: &mut HVector,
        i_row: HighsInt,
        hint: &mut HighsInt,
    ) {
        // Store elements
        for i in 0..aq.pack_count {
            self.pf_index.push(aq.pack_index[i as usize]);
            self.pf_value.push(aq.pack_value[i as usize]);
        }
        let p_logic = self.u_pivot_lookup[i_row as usize];
        let u_start_x = self.u_start[p_logic as usize];
        // C++ reads u_start[p_logic + 1] even when p_logic == num_row - 1, past
        // the `u_start.resize(num_row)` in buildFinish (UB). In practice it gets
        // the stale element at index num_row, i.e. the end of the last U column,
        // which equals u_last_p[num_row - 1] because MPF never modifies U.
        let u_end_x = if ((p_logic + 1) as usize) < self.u_start.len() {
            self.u_start[(p_logic + 1) as usize]
        } else {
            self.u_last_p[p_logic as usize]
        };
        for k in u_start_x..u_end_x {
            self.pf_index.push(self.u_index[k as usize]);
            self.pf_value.push(-self.u_value[k as usize]);
        }
        self.pf_index.push(i_row);
        self.pf_value.push(-self.u_pivot_value[p_logic as usize]);
        self.pf_start.push(self.pf_index.len() as HighsInt);

        for i in 0..ep.pack_count {
            self.pf_index.push(ep.pack_index[i as usize]);
            self.pf_value.push(ep.pack_value[i as usize]);
        }
        self.pf_start.push(self.pf_index.len() as HighsInt);

        // Store pivot
        self.pf_pivot_value.push(aq.array[i_row as usize]);

        // Refactor or not
        self.u_total_x += aq.pack_count + ep.pack_count;
        if self.u_total_x > self.u_merit_x {
            *hint = 1;
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::updateAPF (lines 2475-2505, 31 lines)
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    pub fn update_apf(
        &mut self,
        basic_index: &mut [HighsInt],
        aq: &mut HVector,
        ep: &mut HVector,
        i_row: HighsInt,
    ) {
        // Store elements
        for i in 0..aq.pack_count {
            self.pf_index.push(aq.pack_index[i as usize]);
            self.pf_value.push(aq.pack_value[i as usize]);
        }

        let variable_out = basic_index[i_row as usize];
        if variable_out >= self.num_col {
            self.pf_index.push(variable_out - self.num_col);
            self.pf_value.push(-1.0);
        } else {
            let mut k = self.a_start[variable_out as usize];
            while k < self.a_start[(variable_out + 1) as usize] {
                self.pf_index.push(self.a_index[k as usize]);
                self.pf_value.push(-self.a_value[k as usize]);
                k += 1;
            }
        }
        self.pf_start.push(self.pf_index.len() as HighsInt);

        for i in 0..ep.pack_count {
            self.pf_index.push(ep.pack_index[i as usize]);
            self.pf_value.push(ep.pack_value[i as usize]);
        }
        self.pf_start.push(self.pf_index.len() as HighsInt);

        // Store pivot
        self.pf_pivot_value.push(aq.array[i_row as usize]);
    }

    /// HiGHS: HFactor.cpp::HFactor::getInvert (lines 2507-2537, 31 lines)
    pub fn get_invert(&self) -> InvertibleRepresentation {
        let mut invert = InvertibleRepresentation::default();
        invert.l_pivot_index = self.l_pivot_index.clone();
        invert.l_pivot_lookup = self.l_pivot_lookup.clone();
        invert.l_start = self.l_start.clone();
        invert.l_index = self.l_index.clone();
        invert.l_value = self.l_value.clone();
        invert.lr_start = self.lr_start.clone();
        invert.lr_index = self.lr_index.clone();
        invert.lr_value = self.lr_value.clone();
        invert.u_pivot_lookup = self.u_pivot_lookup.clone();
        invert.u_pivot_index = self.u_pivot_index.clone();
        invert.u_pivot_value = self.u_pivot_value.clone();
        invert.u_start = self.u_start.clone();
        invert.u_last_p = self.u_last_p.clone();
        invert.u_index = self.u_index.clone();
        invert.u_value = self.u_value.clone();
        invert.ur_start = self.ur_start.clone();
        invert.ur_lastp = self.ur_lastp.clone();
        invert.ur_space = self.ur_space.clone();
        invert.ur_index = self.ur_index.clone();
        invert.ur_value = self.ur_value.clone();
        invert.pf_start = self.pf_start.clone();
        invert.pf_index = self.pf_index.clone();
        invert.pf_value = self.pf_value.clone();
        invert.pf_pivot_index = self.pf_pivot_index.clone();
        invert.pf_pivot_value = self.pf_pivot_value.clone();
        invert
    }

    /// HiGHS: HFactor.cpp::HFactor::setInvert (lines 2539-2567, 29 lines)
    pub fn set_invert(&mut self, invert: &InvertibleRepresentation) {
        self.l_pivot_index = invert.l_pivot_index.clone();
        self.l_pivot_lookup = invert.l_pivot_lookup.clone();
        self.l_start = invert.l_start.clone();
        self.l_index = invert.l_index.clone();
        self.l_value = invert.l_value.clone();
        self.lr_start = invert.lr_start.clone();
        self.lr_index = invert.lr_index.clone();
        self.lr_value = invert.lr_value.clone();
        self.u_pivot_lookup = invert.u_pivot_lookup.clone();
        self.u_pivot_index = invert.u_pivot_index.clone();
        self.u_pivot_value = invert.u_pivot_value.clone();
        self.u_start = invert.u_start.clone();
        self.u_last_p = invert.u_last_p.clone();
        self.u_index = invert.u_index.clone();
        self.u_value = invert.u_value.clone();
        self.ur_start = invert.ur_start.clone();
        self.ur_lastp = invert.ur_lastp.clone();
        self.ur_space = invert.ur_space.clone();
        self.ur_index = invert.ur_index.clone();
        self.ur_value = invert.ur_value.clone();
        self.pf_start = invert.pf_start.clone();
        self.pf_index = invert.pf_index.clone();
        self.pf_value = invert.pf_value.clone();
        self.pf_pivot_index = invert.pf_pivot_index.clone();
        self.pf_pivot_value = invert.pf_pivot_value.clone();
    }
}

impl InvertibleRepresentation {
    /// HiGHS: HFactor.cpp::InvertibleRepresentation::clear (lines 2569-2597, 29 lines)
    pub fn clear(&mut self) {
        self.l_pivot_index.clear();
        self.l_pivot_lookup.clear();
        self.l_start.clear();
        self.l_index.clear();
        self.l_value.clear();
        self.lr_start.clear();
        self.lr_index.clear();
        self.lr_value.clear();
        self.u_pivot_lookup.clear();
        self.u_pivot_index.clear();
        self.u_pivot_value.clear();
        self.u_start.clear();
        self.u_last_p.clear();
        self.u_index.clear();
        self.u_value.clear();
        self.ur_start.clear();
        self.ur_lastp.clear();
        self.ur_space.clear();
        self.ur_index.clear();
        self.ur_value.clear();
        self.pf_start.clear();
        self.pf_index.clear();
        self.pf_value.clear();
        self.pf_pivot_index.clear();
        self.pf_pivot_value.clear();
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use super::*;

    /// Hand-derived from HFactor.cpp updateCFT (2052-2298) on a 2x2 identity factor with
    /// empty off-diagonal U, a chain of two HVectors (aq1 pivot row 0, alpha 2; aq2 pivot
    /// row 1, alpha 4) and row_ep chain ep1 = [1, 0], ep2 = [0, 1]. Expected values:
    /// t_pivot = [2, 4]; pf column for pivot 1 = {row 0: 0.5}; pivot lookup [2, 3];
    /// u_total_x ends at 1; pf_start = [0, 0, 1].
    #[test]
    fn w5_pami_cft_update_cft_matches_chained_updates() {
        let mut factor = HFactor::default();
        factor.num_row = 2;
        factor.update_method = crate::util::hfactor_const::kUpdateMethodFt;
        factor.u_pivot_lookup = vec![0, 1];
        factor.u_pivot_index = vec![0, 1];
        factor.u_pivot_value = vec![1.0, 1.0];
        factor.u_start = vec![0, 0];
        factor.u_last_p = vec![0, 0];
        factor.ur_start = vec![0, 0];
        factor.ur_lastp = vec![0, 0];
        factor.ur_space = vec![0, 0];
        factor.pf_start = vec![0];
        factor.iwork = vec![0; 2];
        factor.dwork = vec![0.0; 2];

        let mut aq1 = HVector::default();
        aq1.array = vec![2.0, 1.0];
        aq1.pack_count = 2;
        aq1.pack_index = vec![0, 1];
        aq1.pack_value = vec![2.0, 1.0];
        let mut aq2 = HVector::default();
        aq2.array = vec![0.0, 4.0];
        aq2.pack_count = 1;
        aq2.pack_index = vec![1];
        aq2.pack_value = vec![4.0];
        let mut ep1 = HVector::default();
        ep1.array = vec![1.0, 0.0];
        ep1.pack_count = 1;
        ep1.pack_index = vec![0];
        ep1.pack_value = vec![1.0];
        let mut ep2 = HVector::default();
        ep2.array = vec![0.0, 1.0];
        ep2.pack_count = 1;
        ep2.pack_index = vec![1];
        ep2.pack_value = vec![1.0];

        // Chain aq1 -> aq2 and ep1 -> ep2 (C++ `aq->next`, `ep->next`)
        aq1.next = Some(HVectorNext(std::ptr::NonNull::from(&mut aq2)));
        ep1.next = Some(HVectorNext(std::ptr::NonNull::from(&mut ep2)));

        let mut i_row = vec![0, 1];
        // SAFETY: aq1 / ep1 head live chains; aq2 / ep2 outlive the call and are not
        // otherwise borrowed while it runs.
        unsafe {
            factor.update_cft(&mut aq1, &mut ep1, &mut i_row);
        }

        assert_eq!(factor.u_pivot_lookup, vec![2, 3]);
        assert_eq!(factor.u_pivot_index, vec![-1, -1, 0, 1]);
        assert_eq!(factor.u_pivot_value, vec![1.0, 1.0, 2.0, 4.0]);
        assert_eq!(factor.pf_pivot_index, vec![0, 1]);
        assert_eq!(factor.pf_start, vec![0, 0, 1]);
        assert_eq!(factor.pf_index, vec![0]);
        assert_eq!(factor.pf_value, vec![0.5]);
        assert_eq!(factor.u_total_x, 1);
        assert_eq!(factor.u_index, vec![1]);
        assert_eq!(factor.u_value, vec![1.0]);
    }
}
