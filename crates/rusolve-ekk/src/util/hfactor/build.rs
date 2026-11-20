//! HiGHS source: util/HFactor.cpp (lines 560-869, 1263-1527)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! INVERT: singleton passes (buildSimple), rank-deficiency handling, finish.

use super::*;
use crate::prelude::*;

// Unchecked slice access for the INVERT loops; each use sits in an `unsafe`
// block with a SAFETY note naming why the index is in range, and a
// `debug_assert` re-checks it in debug builds.
macro_rules! rd {
    ($s:expr, $i:expr) => {{
        let i = $i as usize;
        debug_assert!(i < $s.len());
        *$s.get_unchecked(i)
    }};
}
macro_rules! wr {
    ($s:expr, $i:expr, $v:expr) => {{
        let i = $i as usize;
        let v = $v;
        debug_assert!(i < $s.len());
        *$s.get_unchecked_mut(i) = v;
    }};
}

impl HFactor {
    /// HiGHS: HFactor.cpp::HFactor::buildSimple (lines 560-869, 310 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    pub fn build_simple(&mut self, basic_index: &mut [HighsInt]) {
        // 0. Clear L and U factor
        self.lu_clear();

        // Set all values of permute to -1 so that unpermuted (rank
        // deficient) columns can be identified
        let permute_dim = self.num_basic;
        self.permute.resize(permute_dim as usize, 0);
        self.permute.fill(-1);

        // 1. Prepare basis matrix and deal with unit columns
        let mut bcount_x: HighsInt = 0;
        let num_row = self.num_row;
        self.mr_count_before[..num_row as usize].fill(0);
        self.nwork = 0;
        let iwork_dim = self.num_basic;
        self.iwork.resize((iwork_dim + 1) as usize, 0);
        self.iwork.fill(0);
        // The loops below index with values read from the basis and the
        // matrix. Everything is sliced from fields sized in `setup`
        // (`b_index`/`b_value` hold `basis_matrix_limit_size`, the largest
        // possible basis matrix; `b_start`/`b_var` `max(num_row, num_basic)`
        // (+1); `permute` `num_basic`; `iwork` `num_basic + 1`;
        // `mr_count_before` `num_row`), and every index is a basis variable,
        // a row of `A` or a running count bounded by those sizes, so each
        // access is in range. The checked version of this code was a large
        // part of INVERT on fome13.
        {
            let num_col = self.num_col;
            let num_basic = self.num_basic;
            let l_len = self.l_index.len() as HighsInt;
            let u_len = self.u_index.len() as HighsInt;
            self.l_start.reserve(num_basic as usize);
            self.u_pivot_index.reserve(num_basic as usize);
            self.u_pivot_value.reserve(num_basic as usize);
            self.u_start.reserve(num_basic as usize);
            self.refactor_info_.pivot_row.reserve(num_basic as usize);
            self.refactor_info_.pivot_var.reserve(num_basic as usize);
            self.refactor_info_.pivot_type.reserve(num_basic as usize);
            let mr = &mut self.mr_count_before[..];
            let b_index = &mut self.b_index[..];
            let b_value = &mut self.b_value[..];
            let b_start = &mut self.b_start[..];
            let b_var = &mut self.b_var[..];
            let iwork = &mut self.iwork[..];
            let permute = &mut self.permute[..];
            let a_start = &self.a_start[..];
            let a_index = &self.a_index[..];
            let a_value = &self.a_value[..];
            let l_start = &mut self.l_start;
            let u_pivot_index = &mut self.u_pivot_index;
            let u_pivot_value = &mut self.u_pivot_value;
            let u_start = &mut self.u_start;
            let ri_row = &mut self.refactor_info_.pivot_row;
            let ri_var = &mut self.refactor_info_.pivot_var;
            let ri_type = &mut self.refactor_info_.pivot_type;
            let mut nwork: HighsInt = 0;
            for i_col in 0..num_basic {
                let i_mat = basic_index[i_col as usize];
                let mut i_row: HighsInt = -1;
                let mut pivot_type: i8 = kPivotIllegal;
                // SAFETY: see the sizes listed above.
                unsafe {
                    if i_mat >= num_col {
                        // 1.1 Logical column
                        let lc_i_row = i_mat - num_col;
                        if rd!(mr, lc_i_row) >= 0 {
                            pivot_type = kPivotLogical;
                            i_row = lc_i_row;
                        } else {
                            wr!(mr, lc_i_row, rd!(mr, lc_i_row) + 1);
                            wr!(b_index, bcount_x, lc_i_row);
                            wr!(b_value, bcount_x, 1.0);
                            bcount_x += 1;
                            wr!(iwork, nwork, i_col);
                            nwork += 1;
                        }
                    } else {
                        // 1.2 Structural column
                        let start = rd!(a_start, i_mat);
                        let count = rd!(a_start, i_mat + 1) - start;
                        let mut ok_unit_col = false;
                        if count == 1 {
                            ok_unit_col =
                                rd!(a_value, start) == 1.0 && rd!(mr, rd!(a_index, start)) >= 0;
                        }
                        if ok_unit_col {
                            pivot_type = kPivotColSingleton;
                            i_row = rd!(a_index, start);
                        } else {
                            for k in start..start + count {
                                let r = rd!(a_index, k);
                                wr!(mr, r, rd!(mr, r) + 1);
                                wr!(b_index, bcount_x, r);
                                wr!(b_value, bcount_x, rd!(a_value, k));
                                bcount_x += 1;
                            }
                            wr!(iwork, nwork, i_col);
                            nwork += 1;
                        }
                    }

                    if i_row >= 0 {
                        // 1.3 Record unit column
                        wr!(permute, i_col, i_row);
                        l_start.push(l_len);
                        u_pivot_index.push(i_row);
                        u_pivot_value.push(1.0);
                        u_start.push(u_len);
                        wr!(mr, i_row, -num_basic);
                        debug_assert!(pivot_type != kPivotIllegal);
                        ri_row.push(i_row);
                        ri_var.push(i_mat);
                        ri_type.push(pivot_type);
                    }
                    wr!(b_start, i_col + 1, bcount_x);
                    wr!(b_var, i_col, i_mat);
                }
            }
            self.nwork = nwork;
        }
        // Record the number of elements in the basis matrix
        self.basis_matrix_num_el = self.num_row - self.nwork + bcount_x;

        self.build_synthetic_tick += (bcount_x * 60 + (self.num_row - self.nwork) * 80) as f64;

        // 2. Search for and deal with singletons
        let mut t2_search: f64 = 0.0;
        let mut t2_store_l: f64 = self.l_index.len() as f64;
        let mut t2_store_u: f64 = self.u_index.len() as f64;
        let mut t2_store_p: f64 = self.nwork as f64;
        {
            let mut nwork = self.nwork;
            let mr = &mut self.mr_count_before[..];
            let b_index = &self.b_index[..];
            let b_value = &self.b_value[..];
            let b_start = &self.b_start[..];
            let b_var = &self.b_var[..];
            let iwork = &mut self.iwork[..];
            let permute = &mut self.permute[..];
            let l_start = &mut self.l_start;
            let l_index = &mut self.l_index;
            let l_value = &mut self.l_value;
            let u_pivot_index = &mut self.u_pivot_index;
            let u_pivot_value = &mut self.u_pivot_value;
            let u_start = &mut self.u_start;
            let u_index = &mut self.u_index;
            let u_value = &mut self.u_value;
            let ri_row = &mut self.refactor_info_.pivot_row;
            let ri_var = &mut self.refactor_info_.pivot_var;
            let ri_type = &mut self.refactor_info_.pivot_type;
            while nwork > 0 {
                let nwork_last = nwork;
                nwork = 0;
                for i in 0..nwork_last {
                    // SAFETY: `i < nwork_last <= num_basic < iwork.len()`;
                    // `i_col < num_basic`; `b_index[k]` for `k` in a column's
                    // range is a row `< num_row`; `b_start[i_col + 1]` is in
                    // range; see the sizes listed in step 1.
                    unsafe {
                        let i_col = rd!(iwork, i);
                        let start = rd!(b_start, i_col);
                        let end = rd!(b_start, i_col + 1);
                        let mut pivot_k: HighsInt = -1;
                        let mut found_row_singleton: HighsInt = 0;
                        let mut count: HighsInt = 0;

                        // 2.1 Search for singleton
                        t2_search += (end - start) as f64;
                        for k in start..end {
                            let mrc = rd!(mr, rd!(b_index, k));
                            if mrc == 1 {
                                pivot_k = k;
                                found_row_singleton = 1;
                                break;
                            }
                            if mrc > 1 {
                                pivot_k = k;
                                count += 1;
                            }
                        }

                        if found_row_singleton != 0 {
                            // 2.2 Deal with row singleton
                            let pivot_multiplier = 1.0 / rd!(b_value, pivot_k);
                            for section in 0..2 {
                                let p0 = if section == 0 { start } else { pivot_k + 1 };
                                let p1 = if section == 0 { pivot_k } else { end };
                                for k in p0..p1 {
                                    let i_row = rd!(b_index, k);
                                    if rd!(mr, i_row) > 0 {
                                        l_index.push(i_row);
                                        l_value.push(rd!(b_value, k) * pivot_multiplier);
                                    } else {
                                        u_index.push(i_row);
                                        u_value.push(rd!(b_value, k));
                                    }
                                    wr!(mr, i_row, rd!(mr, i_row) - 1);
                                }
                            }
                            let i_row = rd!(b_index, pivot_k);
                            wr!(mr, i_row, 0);
                            wr!(permute, i_col, i_row);
                            l_start.push(l_index.len() as HighsInt);

                            u_pivot_index.push(i_row);
                            u_pivot_value.push(rd!(b_value, pivot_k));
                            u_start.push(u_index.len() as HighsInt);
                            debug_assert!(rd!(b_var, i_col) == basic_index[i_col as usize]);

                            ri_row.push(i_row);
                            ri_var.push(basic_index[i_col as usize]);
                            ri_type.push(kPivotRowSingleton);
                        } else if count == 1 {
                            // 2.3 Deal with column singleton
                            for k in start..pivot_k {
                                u_index.push(rd!(b_index, k));
                                u_value.push(rd!(b_value, k));
                            }
                            for k in pivot_k + 1..end {
                                u_index.push(rd!(b_index, k));
                                u_value.push(rd!(b_value, k));
                            }

                            let i_row = rd!(b_index, pivot_k);
                            wr!(mr, i_row, 0);
                            wr!(permute, i_col, i_row);
                            l_start.push(l_index.len() as HighsInt);

                            u_pivot_index.push(i_row);
                            u_pivot_value.push(rd!(b_value, pivot_k));
                            u_start.push(u_index.len() as HighsInt);
                            debug_assert!(rd!(b_var, i_col) == basic_index[i_col as usize]);
                            ri_row.push(i_row);
                            ri_var.push(basic_index[i_col as usize]);
                            ri_type.push(kPivotColSingleton);
                        } else {
                            wr!(iwork, nwork, i_col);
                            nwork += 1;
                        }
                    }
                }

                // No singleton found in the last pass
                if nwork_last == nwork {
                    break;
                }
            }
            self.nwork = nwork;
        }
        t2_store_l = self.l_index.len() as f64 - t2_store_l;
        t2_store_u = self.u_index.len() as f64 - t2_store_u;
        t2_store_p -= self.nwork as f64;

        self.build_synthetic_tick +=
            t2_search * 20.0 + (t2_store_p + t2_store_l + t2_store_u) * 80.0;

        // 3. Prepare the kernel parts
        // 3.1 Prepare row links, row matrix spaces
        self.row_link_first.clear();
        self.row_link_first
            .resize((self.num_basic + 1) as usize, -1);
        self.mr_count.clear();
        self.mr_count.resize(self.num_row as usize, 0);
        let mut mr_count_x: HighsInt = 0;
        // Determine the number of entries in the kernel
        self.kernel_num_el = 0;
        for i_row in 0..self.num_row {
            let count = self.mr_count_before[i_row as usize];
            if count > 0 {
                self.mr_start[i_row as usize] = mr_count_x;
                self.mr_space[i_row as usize] = count * 2;
                mr_count_x += count * 2;
                self.rlink_add(i_row, count);
                self.kernel_num_el += count + 1;
            }
        }
        self.mr_index.resize(mr_count_x as usize, 0);

        // 3.2 Prepare column links, kernel matrix
        self.col_link_first.clear();
        self.col_link_first.resize((self.num_row + 1) as usize, -1);
        let mc_dim = self.num_basic;
        self.mc_index.clear();
        self.mc_value.clear();
        self.mc_count_a.clear();
        self.mc_count_a.resize(mc_dim as usize, 0);
        self.mc_count_n.clear();
        self.mc_count_n.resize(mc_dim as usize, 0);
        let mut mc_count_x: HighsInt = 0;
        for i in 0..self.nwork {
            let i_col = self.iwork[i as usize];
            self.mc_var[i_col as usize] = self.b_var[i_col as usize];
            self.mc_start[i_col as usize] = mc_count_x;
            self.mc_space[i_col as usize] =
                (self.b_start[(i_col + 1) as usize] - self.b_start[i_col as usize]) * 2;
            mc_count_x += self.mc_space[i_col as usize];
            self.mc_index.resize(mc_count_x as usize, 0);
            self.mc_value.resize(mc_count_x as usize, 0.0);
            for k in self.b_start[i_col as usize]..self.b_start[(i_col + 1) as usize] {
                let i_row = self.b_index[k as usize];
                let value = self.b_value[k as usize];
                if self.mr_count_before[i_row as usize] > 0 {
                    self.col_insert(i_col, i_row, value);
                    self.row_insert(i_col, i_row);
                } else {
                    self.col_store_n(i_col, i_row, value);
                }
            }
            self.col_fix_max(i_col);
            self.clink_add(i_col, self.mc_count_a[i_col as usize]);
        }
        self.build_synthetic_tick +=
            ((self.num_row + self.nwork + mc_count_x) * 40 + mr_count_x * 20) as f64;
        // Record the kernel dimension
        self.kernel_dim = self.nwork;
        debug_assert!(
            self.refactor_info_.pivot_row.len() as HighsInt == self.num_basic - self.nwork
        );
    }

    /// HiGHS: HFactor.cpp::HFactor::buildHandleRankDeficiency (lines 1263-1363, 101 lines)
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    pub fn build_handle_rank_deficiency(&mut self, basic_index: &mut [HighsInt]) {
        debug_report_rank_deficiency(
            0,
            self.highs_debug_level,
            &self.log_options,
            self.num_row,
            &self.permute,
            &self.iwork,
            basic_index,
            self.rank_deficiency,
            &self.row_with_no_pivot,
            &self.col_with_no_pivot,
        );
        if self.num_basic < self.num_row {
            self.rank_deficiency += self.num_row - self.num_basic;
        }
        self.row_with_no_pivot
            .resize(self.rank_deficiency as usize, 0);
        self.col_with_no_pivot
            .resize(self.rank_deficiency as usize, 0);
        let mut lc_rank_deficiency: HighsInt = 0;
        if self.num_basic < self.num_row {
            self.iwork.resize(self.num_row as usize, 0);
        } else if self.num_basic > self.num_row {
            self.iwork.resize(self.num_basic as usize, 0);
        }
        for i in 0..self.num_row {
            self.iwork[i as usize] = -1;
        }
        for i in 0..self.num_basic {
            let perm_i = self.permute[i as usize];
            if perm_i >= 0 {
                self.iwork[perm_i as usize] = basic_index[i as usize];
            } else {
                self.col_with_no_pivot[lc_rank_deficiency as usize] = i;
                lc_rank_deficiency += 1;
            }
        }
        if self.num_basic < self.num_row {
            self.permute.resize(self.num_row as usize, 0);
            for i in self.num_basic..self.num_row {
                self.col_with_no_pivot[lc_rank_deficiency as usize] = i;
                lc_rank_deficiency += 1;
                self.permute[i as usize] = -1;
            }
        }
        debug_assert!(lc_rank_deficiency == self.rank_deficiency);
        lc_rank_deficiency = 0;
        for i in 0..self.num_row {
            if self.iwork[i as usize] < 0 {
                self.row_with_no_pivot[lc_rank_deficiency as usize] = i;
                self.iwork[i as usize] = -(lc_rank_deficiency + 1);
                lc_rank_deficiency += 1;
            }
        }
        if self.num_row < self.num_basic {
            for i in self.num_row..self.num_basic {
                self.row_with_no_pivot[lc_rank_deficiency as usize] = i;
                self.iwork[i as usize] = -(lc_rank_deficiency + 1);
                lc_rank_deficiency += 1;
            }
        }
        debug_assert!(lc_rank_deficiency == self.rank_deficiency);
        debug_report_rank_deficiency(
            1,
            self.highs_debug_level,
            &self.log_options,
            self.num_row,
            &self.permute,
            &self.iwork,
            basic_index,
            self.rank_deficiency,
            &self.row_with_no_pivot,
            &self.col_with_no_pivot,
        );
        let row_rank_deficiency = self.rank_deficiency - cpp_max(self.num_basic - self.num_row, 0);
        for k in 0..self.rank_deficiency {
            let i_row = self.row_with_no_pivot[k as usize];
            let i_col = self.col_with_no_pivot[k as usize];
            debug_assert!(self.permute[i_col as usize] == -1);
            self.permute[i_col as usize] = i_row;
            if k < row_rank_deficiency {
                self.l_start.push(self.l_index.len() as HighsInt);
                self.u_pivot_index.push(i_row);
                self.u_pivot_value.push(1.0);
                self.u_start.push(self.u_index.len() as HighsInt);
            }
        }
        debug_report_rank_deficiency(
            2,
            self.highs_debug_level,
            &self.log_options,
            self.num_row,
            &self.permute,
            &self.iwork,
            basic_index,
            self.rank_deficiency,
            &self.row_with_no_pivot,
            &self.col_with_no_pivot,
        );
        debug_report_rank_deficient_asm(
            self.highs_debug_level,
            &self.log_options,
            self.num_row,
            &self.mc_start,
            &self.mc_count_a,
            &self.mc_index,
            &self.mc_value,
            &self.iwork,
            self.rank_deficiency,
            &self.col_with_no_pivot,
            &self.row_with_no_pivot,
        );
    }

    /// HiGHS: HFactor.cpp::HFactor::buildMarkSingC (lines 1365-1398, 34 lines)
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    pub fn build_mark_sing_c(&mut self, basic_index: &mut [HighsInt]) {
        debug_report_mark_sing_c(
            0,
            self.highs_debug_level,
            &self.log_options,
            self.num_row,
            &self.iwork,
            basic_index,
        );

        let basic_index_rank_deficiency =
            self.rank_deficiency - cpp_max(self.num_row - self.num_basic, 0);
        self.var_with_no_pivot
            .resize(self.rank_deficiency as usize, 0);
        for k in 0..self.rank_deficiency {
            let asm_row = self.row_with_no_pivot[k as usize];
            let asm_col = self.col_with_no_pivot[k as usize];
            debug_assert!(asm_row < self.iwork.len() as HighsInt);
            debug_assert!(
                -self.iwork[asm_row as usize] > 0
                    && -self.iwork[asm_row as usize] - 1 < self.rank_deficiency
            );
            // Store negation of 1+ASMcol so that removing column 0 can be
            // identified!
            self.iwork[asm_row as usize] = -(asm_col + 1);
            // Only update basic_index for the true entries
            if asm_col < self.num_basic {
                debug_assert!(k < basic_index_rank_deficiency);
                self.var_with_no_pivot[k as usize] = basic_index[asm_col as usize];
                basic_index[asm_col as usize] = self.num_col + asm_row;
            } else if self.num_basic < self.num_row {
                debug_assert!(asm_col == self.num_basic + k - basic_index_rank_deficiency);
                self.var_with_no_pivot[k as usize] = -1;
            }
        }
        debug_report_mark_sing_c(
            1,
            self.highs_debug_level,
            &self.log_options,
            self.num_row,
            &self.iwork,
            basic_index,
        );
    }

    /// HiGHS: HFactor.cpp::HFactor::buildFinish (lines 1400-1501, 102 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    pub fn build_finish(&mut self, basic_index: &mut [HighsInt]) {
        debug_assert!(self.num_basic >= self.num_row);
        let num_row = self.num_row;
        let nr = num_row as usize;
        // The look up table
        // SAFETY (this function): `u_pivot_index` holds `num_row` rows
        // `< num_row` (a full pivot sequence), the lookups, `iwork`, `lr_*`,
        // `ur_*` and the starts are sized to `num_row` (+1) or to the L / U
        // counts below, and `l_index` / `u_index` entries are rows, so each
        // unchecked access is in range.
        {
            let lookup = &mut self.u_pivot_lookup[..];
            let piv = &self.u_pivot_index[..];
            assert!(piv.len() >= nr && lookup.len() >= nr);
            for i in 0..nr {
                // SAFETY: `i < nr <= piv.len()`; the pivot row is `< num_row <= lookup.len()`.
                unsafe {
                    wr!(lookup, rd!(piv, i), i as HighsInt);
                }
            }
        }
        self.l_pivot_index.clone_from(&self.u_pivot_index);
        self.l_pivot_lookup.clone_from(&self.u_pivot_lookup);

        // LR space
        let l_count_x = self.l_index.len() as HighsInt;
        self.lr_index.resize(l_count_x as usize, 0);
        self.lr_value.resize(l_count_x as usize, 0.0);

        // LR pointer
        self.iwork.clear();
        self.iwork.resize(nr, 0);
        {
            let iwork = &mut self.iwork[..];
            let lookup = &self.l_pivot_lookup[..];
            for &r in &self.l_index[..l_count_x as usize] {
                // SAFETY: `r < num_row`, `lookup[r] < num_row`, `iwork.len() == nr`.
                unsafe {
                    let q = rd!(lookup, r);
                    wr!(iwork, q, rd!(iwork, q) + 1);
                }
            }
        }

        // lr_start[i] is the running sum of the counts before i.
        self.lr_start.clear();
        let mut acc: HighsInt = 0;
        for &cnt in &self.iwork[..nr] {
            self.lr_start.push(acc);
            acc += cnt;
        }
        self.lr_start.push(acc);

        // LR elements
        self.iwork.clear();
        self.iwork.extend_from_slice(&self.lr_start[0..nr]);
        {
            let iwork = &mut self.iwork[..];
            let lookup = &self.l_pivot_lookup[..];
            let l_index = &self.l_index[..];
            let l_value = &self.l_value[..];
            let l_start = &self.l_start[..];
            let lr_index = &mut self.lr_index[..];
            let lr_value = &mut self.lr_value[..];
            for i in 0..nr {
                let index = self.l_pivot_index[i];
                for k in l_start[i]..l_start[i + 1] {
                    // SAFETY: `k` is in `l_start[i]..l_start[i + 1]`, within
                    // the `l_count_x` entries; the target slot is below the
                    // column count `lr_start[row + 1]`.
                    unsafe {
                        let i_row = rd!(lookup, rd!(l_index, k));
                        let i_put = rd!(iwork, i_row);
                        wr!(iwork, i_row, i_put + 1);
                        wr!(lr_index, i_put, index);
                        wr!(lr_value, i_put, rd!(l_value, k));
                    }
                }
            }
        }

        // U pointer
        self.u_start.push(0);
        self.u_last_p.clear();
        self.u_last_p.extend_from_slice(&self.u_start[1..nr + 1]);
        self.u_start.resize(nr, 0);

        // UR space
        let u_count_x = self.u_index.len() as HighsInt;
        let ur_stuff_size: HighsInt = if self.update_method == kUpdateMethodFt {
            5
        } else {
            0
        };
        let ur_count_size = u_count_x + ur_stuff_size * num_row;
        self.ur_index.resize(ur_count_size as usize, 0);
        self.ur_value.resize(ur_count_size as usize, 0.0);

        // UR pointer
        //
        // NB ur_lastp just being used as temporary storage here
        self.ur_lastp.clear();
        self.ur_lastp.resize(nr, 0);
        self.ur_space.clear();
        self.ur_space.resize(nr, ur_stuff_size);
        {
            let ur_lastp = &mut self.ur_lastp[..];
            let lookup = &self.u_pivot_lookup[..];
            for &r in &self.u_index[..u_count_x as usize] {
                // SAFETY: `r < num_row`, `lookup[r] < num_row == ur_lastp.len()`.
                unsafe {
                    let q = rd!(lookup, r);
                    wr!(ur_lastp, q, rd!(ur_lastp, q) + 1);
                }
            }
        }
        // ur_start[i] is the running sum of the counts (plus update space)
        // before i; the C++ fills nr + 1 entries and drops the last.
        self.ur_start.clear();
        let mut acc: HighsInt = 0;
        for &cnt in &self.ur_lastp[..nr] {
            self.ur_start.push(acc);
            acc += cnt + ur_stuff_size;
        }

        // UR element
        //
        // NB ur_lastp initialised here!
        self.ur_lastp.clone_from(&self.ur_start);
        {
            let ur_lastp = &mut self.ur_lastp[..];
            let lookup = &self.u_pivot_lookup[..];
            let u_index = &self.u_index[..];
            let u_value = &self.u_value[..];
            let u_start = &self.u_start[..];
            let u_last_p = &self.u_last_p[..];
            let ur_index = &mut self.ur_index[..];
            let ur_value = &mut self.ur_value[..];
            for i in 0..nr {
                let index = self.u_pivot_index[i];
                for k in u_start[i]..u_last_p[i] {
                    // SAFETY: `k < u_count_x`; the target slot is below the
                    // row's `ur_start` plus its count plus the update space.
                    unsafe {
                        let i_row = rd!(lookup, rd!(u_index, k));
                        let i_put = rd!(ur_lastp, i_row);
                        wr!(ur_lastp, i_row, i_put + 1);
                        wr!(ur_index, i_put, index);
                        wr!(ur_value, i_put, rd!(u_value, k));
                    }
                }
            }
        }

        // Re-factor merit
        self.u_merit_x = (num_row as f64 + (l_count_x + u_count_x) as f64 * 1.5) as HighsInt;
        self.u_total_x = u_count_x;
        if self.update_method == kUpdateMethodPf {
            self.u_merit_x = num_row + u_count_x * 4;
        }
        if self.update_method == kUpdateMethodMpf {
            self.u_merit_x = num_row + u_count_x * 3;
        }

        // Clear update buffer
        self.pf_pivot_value.clear();
        self.pf_pivot_index.clear();
        self.pf_start.clear();
        self.pf_start.push(0);
        self.pf_index.clear();
        self.pf_value.clear();

        if !self.refactor_info_.r#use {
            // Finally, if not calling buildFinish after refactorizing,
            // permute the basic variables
            self.iwork.clear();
            self.iwork
                .extend_from_slice(&basic_index[0..self.num_basic as usize]);
            for i in 0..self.num_basic as usize {
                basic_index[self.permute[i] as usize] = self.iwork[i];
            }
            for i in self.num_row..self.num_basic {
                debug_assert!(basic_index[i as usize] >= self.num_col);
            }
            // Add cost of buildFinish to build_synthetic_tick
            self.build_synthetic_tick += (num_row * 80 + (l_count_x + u_count_x) * 60) as f64;
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::zeroCol (lines 1503-1527, 25 lines)
    /// HiGHS 1.16 (#3354): HFactor.cpp::HFactor::zeroCol (lines 1614-1642)
    pub fn zero_col(&mut self, j_col: HighsInt, pivot_multiplier: f64) {
        let a_count = self.mc_count_a[j_col as usize];
        let a_start = self.mc_start[j_col as usize];
        let a_end = a_start + a_count;
        for i_el in a_start..a_end {
            let abs_value = self.mc_value[i_el as usize].abs();
            let i_row = self.mc_index[i_el as usize];
            let original_row_count = self.mr_count[i_row as usize];
            // Remove the column index from this row of the row-wise
            // matrix. Also decreases the row count
            self.row_delete(j_col, i_row);
            // Remove this row from the linked list of rows containing it
            self.rlink_del(i_row);
            // Add the this row to the linked list of rows with this reduced
            // count
            debug_assert!(self.mr_count[i_row as usize] == original_row_count - 1);
            self.rlink_add(i_row, self.mr_count[i_row as usize]);
            // HiGHS 1.16 #3354 (h16_hfactor_pivot): HFactor.cpp lines 1631-1635 (1.16.0)
            if self.h16_hfactor_pivot {
                debug_assert!(abs_value <= pivot_multiplier.abs());
            } else {
                debug_assert!(abs_value < self.pivot_tolerance);
            }
        }
        // Remove the column from the linked list of columns containing it
        self.clink_del(j_col);
        // Zero the counts of the active and inactive sections of the column
        self.mc_count_a[j_col as usize] = 0;
        self.mc_count_n[j_col as usize] = 0;
    }
}
