//! HiGHS source: util/HFactor.cpp (lines 34-158, 1529-2050)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Triangular solves: the free hyper-sparse kernels `solveMatrixT`/`solveHyper` and
//! FTRAN/BTRAN through L, U and the FT/PF/MPF/APF update forms.

use super::*;
use crate::prelude::*;

/// PORT NOTE (perf): the C++ scatter loop
/// `for (k = start; k < end; k++) rhs_array[index[k]] -= pivot_multiplier * value[k];`
/// over zipped subslices, so only the scatter target is bounds-checked. Same operations
/// in the same order.
#[inline(always)]
fn scatter_sub(
    rhs_array: &mut [f64],
    index: &[HighsInt],
    value: &[f64],
    start: HighsInt,
    end: HighsInt,
    pivot_multiplier: f64,
) {
    if start < end {
        let (start, end) = (start as usize, end as usize);
        for (&i, &v) in index[start..end].iter().zip(&value[start..end]) {
            debug_assert!((i as usize) < rhs_array.len());
            // SAFETY: every index stored in the L/U/LR/UR copies is a row index in
            // 0..num_row, and the RHS array has length >= num_row (HVector::setup(num_row)),
            // so the scatter target is in bounds. Same operation, same order.
            unsafe {
                *rhs_array.get_unchecked_mut(i as usize) -= pivot_multiplier * v;
            }
        }
    }
}

/// HiGHS: HFactor.cpp::solveMatrixT (lines 34-59, 26 lines)
pub fn solve_matrix_t(
    x_start: HighsInt,
    x_end: HighsInt,
    y_start: HighsInt,
    y_end: HighsInt,
    t_index: &[HighsInt],
    t_value: &[f64],
    t_pivot: f64,
    rhs_count: &mut HighsInt,
    rhs_index: &mut [HighsInt],
    rhs_array: &mut [f64],
) {
    // Collect by X
    let mut pivot_multiplier: f64 = 0.0;
    let mut k = x_start;
    while k < x_end {
        pivot_multiplier += t_value[k as usize] * rhs_array[t_index[k as usize] as usize];
        k += 1;
    }

    // Scatter by Y
    if pivot_multiplier.abs() > kHighsTiny {
        let mut work_count = *rhs_count;

        pivot_multiplier /= t_pivot;
        let mut k = y_start;
        while k < y_end {
            let index = t_index[k as usize];
            let value0 = rhs_array[index as usize];
            let value1 = value0 - pivot_multiplier * t_value[k as usize];
            if value0 == 0.0 {
                rhs_index[work_count as usize] = index;
                work_count += 1;
            }
            rhs_array[index as usize] = if value1.abs() < kHighsTiny {
                kHighsZero
            } else {
                value1
            };
            k += 1;
        }

        *rhs_count = work_count;
    }
}

/// HiGHS: HFactor.cpp::solveHyper (lines 61-158, 98 lines)
/// MUST-PORT (determinism, synthetic ticks drive reinversion): updates `rhs.synthetic_tick`; keep every increment and its operand order.
/// `h_pivot_value` is `None` where C++ passes `0` (the L solves); `if (h_pivot_value == 0)` is
/// `h_pivot_value.is_none()`. C++ `h_end = &l_start[1]` (pointer offset) is passed as
/// `&self.l_start[1..]`.
pub fn solve_hyper(
    h_size: HighsInt,
    h_lookup: &[HighsInt],
    h_pivot_index: &[HighsInt],
    h_pivot_value: Option<&[f64]>,
    h_start: &[HighsInt],
    h_end: &[HighsInt],
    h_index: &[HighsInt],
    h_value: &[f64],
    rhs: &mut HVector,
) {
    let mut rhs_count = rhs.count;
    // `work_buffers` sizes `cwork` / `iwork` as C++ `HVectorBase::setup` does.
    rhs.work_buffers();
    let rhs_index = &mut rhs.index;
    let rhs_array = &mut rhs.array;

    // Take count

    // Build list
    let list_mark = &mut rhs.cwork;
    // list_index = iwork[0..], list_stack = iwork[h_size..]
    let iwork = &mut rhs.iwork;
    let stack_off = h_size as usize;
    let mut list_count: HighsInt = 0;

    let mut count_pivot: HighsInt = 0;
    let mut count_entry: HighsInt = 0;

    // PORT NOTE (perf): the depth-first search below indexes eight arrays per step; the
    // accessors skip the bounds checks (debug builds keep them through debug_assert).
    // SAFETY (all `rd!`/`wr!` uses): the C++ reads and writes the same positions unchecked.
    // Every index is a pivot position (< h_start.len() == h_end.len() == h_lookup range),
    // a row index (< h_lookup.len()), or a stack slot (< 2 * h_size + 1 <= iwork.len()); the
    // marks live in `cwork`, which `work_buffers` (above) sizes for all pivot positions. The debug
    // assertions in the macros check this in test builds.
    macro_rules! rd {
        ($s:expr, $i:expr) => {{
            debug_assert!(($i as usize) < $s.len());
            *unsafe { $s.get_unchecked($i as usize) }
        }};
    }
    macro_rules! wr {
        ($s:expr, $i:expr, $v:expr) => {{
            debug_assert!(($i as usize) < $s.len());
            *unsafe { $s.get_unchecked_mut($i as usize) } = $v;
        }};
    }
    for i in 0..rhs_count {
        // Skip touched index
        let rhs_i = rd!(rhs_index, i);
        let i_trans = rd!(h_lookup, rhs_i);
        if rd!(list_mark, i_trans) != 0 {
            continue;
        }

        let mut hi = i_trans; // H matrix pivot index
        let mut hk = rd!(h_start, hi); // H matrix non zero position
        let mut n_stack: HighsInt = -1; // Usage of the stack (-1 not used)

        wr!(list_mark, hi, 1); // Mark this as touched

        loop {
            if hk < rd!(h_end, hi) {
                let h_row = rd!(h_index, hk);
                let hi_sub = rd!(h_lookup, h_row);
                hk += 1;
                if rd!(list_mark, hi_sub) == 0 {
                    // Go to a child
                    wr!(list_mark, hi_sub, 1); // Mark as touched
                    n_stack += 1;
                    wr!(iwork, stack_off + n_stack as usize, hi); // Store current into stack
                    n_stack += 1;
                    wr!(iwork, stack_off + n_stack as usize, hk);
                    hi = hi_sub; // Replace current with child
                    hk = rd!(h_start, hi);
                    if hi >= h_size {
                        count_pivot += 1;
                        count_entry += rd!(h_end, hi) - rd!(h_start, hi);
                    }
                }
            } else {
                wr!(iwork, list_count, hi);
                list_count += 1;
                if n_stack == -1 {
                    // Quit on empty stack
                    break;
                }
                hk = rd!(iwork, stack_off + n_stack as usize); // Back to last in stack
                n_stack -= 1;
                hi = rd!(iwork, stack_off + n_stack as usize);
                n_stack -= 1;
            }
        }
    }

    rhs.synthetic_tick += (count_pivot * 20 + count_entry * 10) as f64;

    // Solve with list
    match h_pivot_value {
        None => {
            rhs_count = 0;
            let mut i_list = list_count - 1;
            while i_list >= 0 {
                let i = iwork[i_list as usize];
                list_mark[i as usize] = 0;
                let pivot_row = h_pivot_index[i as usize];
                let pivot_multiplier = rhs_array[pivot_row as usize];
                if pivot_multiplier.abs() > kHighsTiny {
                    rhs_index[rhs_count as usize] = pivot_row;
                    rhs_count += 1;
                    let start = h_start[i as usize];
                    let end = h_end[i as usize];
                    scatter_sub(rhs_array, h_index, h_value, start, end, pivot_multiplier);
                } else {
                    rhs_array[pivot_row as usize] = 0.0;
                }
                i_list -= 1;
            }
            rhs.count = rhs_count;
        }
        Some(h_pivot_value) => {
            rhs_count = 0;
            let mut i_list = list_count - 1;
            while i_list >= 0 {
                let i = iwork[i_list as usize];
                list_mark[i as usize] = 0;
                let pivot_row = h_pivot_index[i as usize];
                let mut pivot_multiplier = rhs_array[pivot_row as usize];
                if pivot_multiplier.abs() > kHighsTiny {
                    pivot_multiplier /= h_pivot_value[i as usize];
                    rhs_array[pivot_row as usize] = pivot_multiplier;
                    rhs_index[rhs_count as usize] = pivot_row;
                    rhs_count += 1;
                    let start = h_start[i as usize];
                    let end = h_end[i as usize];
                    scatter_sub(rhs_array, h_index, h_value, start, end, pivot_multiplier);
                } else {
                    rhs_array[pivot_row as usize] = 0.0;
                }
                i_list -= 1;
            }
            rhs.count = rhs_count;
        }
    }
}

impl HFactor {
    /// HiGHS: HFactor.cpp::HFactor::ftranL (lines 1529-1584, 56 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    pub fn ftran_l(
        &self,
        rhs: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        let _ = factor_timer_clock_pointer;
        if self.update_method == kUpdateMethodApf {
            debug_assert!(self.update_method != kUpdateMethodApf);
            rhs.tight();
            rhs.pack();
            self.ftran_apf(rhs);
            rhs.tight();
        }

        // Determine style of solve
        let current_density = 1.0 * (rhs.count as f64) * self.inv_num_row;
        let sparse_solve =
            rhs.count < 0 || current_density > kHyperCancel || expected_density > kHyperFtranL;
        if sparse_solve {
            let rhs_index = &mut rhs.index;
            let rhs_array = &mut rhs.array;
            // PORT NOTE (perf): pivot/start arrays sliced to `num_row` so `i` needs no
            // bounds checks.
            let n = self.num_row.max(0) as usize;
            let l_pivot_index = &self.l_pivot_index[..n];
            let l_start = &self.l_start[..=n];
            let l_index = &self.l_index;
            let l_value = &self.l_value;
            // Local accumulation of RHS count
            let mut rhs_count: HighsInt = 0;
            // Transform
            for i in 0..n {
                let pivot_row = l_pivot_index[i];
                let pivot_multiplier = rhs_array[pivot_row as usize];
                if pivot_multiplier.abs() > kHighsTiny {
                    rhs_index[rhs_count as usize] = pivot_row;
                    rhs_count += 1;
                    let start = l_start[i];
                    let end = l_start[i + 1];
                    scatter_sub(rhs_array, l_index, l_value, start, end, pivot_multiplier);
                } else {
                    rhs_array[pivot_row as usize] = 0.0;
                }
            }
            // Save the count
            rhs.count = rhs_count;
        } else {
            // Hyper-sparse solve
            solve_hyper(
                self.num_row,
                &self.l_pivot_lookup,
                &self.l_pivot_index,
                None,
                &self.l_start,
                &self.l_start[1..],
                &self.l_index,
                &self.l_value,
                rhs,
            );
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::btranL (lines 1586-1641, 56 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    pub fn btran_l(
        &self,
        rhs: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        let _ = factor_timer_clock_pointer;
        // Determine style of solve
        let current_density = 1.0 * (rhs.count as f64) * self.inv_num_row;
        let sparse_solve =
            rhs.count < 0 || current_density > kHyperCancel || expected_density > kHyperBtranL;
        if sparse_solve {
            let rhs_index = &mut rhs.index;
            let rhs_array = &mut rhs.array;
            // PORT NOTE (perf): pivot/start arrays sliced to `num_row` so `i` needs no
            // bounds checks.
            let n = self.num_row.max(0) as usize;
            let l_pivot_index = &self.l_pivot_index[..n];
            let lr_start = &self.lr_start[..=n];
            let lr_index = &self.lr_index;
            let lr_value = &self.lr_value;
            // Local accumulation of RHS count
            let mut rhs_count: HighsInt = 0;
            // Transform
            for i in (0..n).rev() {
                let pivot_row = l_pivot_index[i];
                let pivot_multiplier = rhs_array[pivot_row as usize];
                if pivot_multiplier.abs() > kHighsTiny {
                    rhs_index[rhs_count as usize] = pivot_row;
                    rhs_count += 1;
                    rhs_array[pivot_row as usize] = pivot_multiplier;
                    let start = lr_start[i];
                    let end = lr_start[i + 1];
                    scatter_sub(rhs_array, lr_index, lr_value, start, end, pivot_multiplier);
                } else {
                    rhs_array[pivot_row as usize] = 0.0;
                }
            }
            // Save the count
            rhs.count = rhs_count;
        } else {
            // Hyper-sparse solve
            solve_hyper(
                self.num_row,
                &self.l_pivot_lookup,
                &self.l_pivot_index,
                None,
                &self.lr_start,
                &self.lr_start[1..],
                &self.lr_index,
                &self.lr_value,
                rhs,
            );
        }

        if self.update_method == kUpdateMethodApf {
            self.btran_apf(rhs);
            rhs.tight();
            rhs.pack();
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::ftranU (lines 1643-1758, 116 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    pub fn ftran_u(
        &self,
        rhs: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        let _ = factor_timer_clock_pointer;
        debug_assert!(rhs.count >= 0);
        // The update part
        if self.update_method == kUpdateMethodFt {
            self.ftran_ft(rhs);
            rhs.tight();
            rhs.pack();
        } else if self.update_method == kUpdateMethodMpf {
            debug_assert!(self.update_method != kUpdateMethodMpf);
            self.ftran_mpf(rhs);
            rhs.tight();
            rhs.pack();
        }

        // The regular part
        //
        // Determine style of solve
        let current_density = 1.0 * (rhs.count as f64) * self.inv_num_row;
        let sparse_solve =
            rhs.count < 0 || current_density > kHyperCancel || expected_density > kHyperFtranU;
        if sparse_solve {
            // Alias to non constant
            let mut rhs_synthetic_tick: f64 = 0.0;
            // Alias to RHS
            let rhs_index = &mut rhs.index;
            let rhs_array = &mut rhs.array;
            // Alias to factor U
            let u_start = &self.u_start;
            let u_end = &self.u_last_p;
            let u_index = &self.u_index;
            let u_value = &self.u_value;
            // Local accumulation of RHS count
            let mut rhs_count: HighsInt = 0;
            // Transform
            let u_pivot_count = self.u_pivot_index.len() as HighsInt;
            // PORT NOTE (perf): per-pivot arrays sliced to `u_pivot_count` so `i_logic`
            // needs no bounds checks.
            let n = u_pivot_count as usize;
            let num_row = self.num_row as usize;
            let u_pivot_index = &self.u_pivot_index[..n];
            let u_pivot_value = &self.u_pivot_value[..n];
            let u_start = &u_start[..n];
            let u_end = &u_end[..n];
            for i_logic in (0..n).rev() {
                // Skip void
                let pivot_row = u_pivot_index[i_logic];
                if pivot_row != -1 {
                    // Normal part
                    let mut pivot_multiplier = rhs_array[pivot_row as usize];
                    if pivot_multiplier.abs() > kHighsTiny {
                        pivot_multiplier /= u_pivot_value[i_logic];
                        rhs_index[rhs_count as usize] = pivot_row;
                        rhs_count += 1;
                        rhs_array[pivot_row as usize] = pivot_multiplier;
                        let start = u_start[i_logic];
                        let end = u_end[i_logic];
                        if i_logic >= num_row {
                            rhs_synthetic_tick += (end - start) as f64;
                        }
                        scatter_sub(rhs_array, u_index, u_value, start, end, pivot_multiplier);
                    } else {
                        rhs_array[pivot_row as usize] = 0.0;
                    }
                }
            }
            // Save the count
            rhs.count = rhs_count;
            rhs.synthetic_tick +=
                rhs_synthetic_tick * 15.0 + ((u_pivot_count - self.num_row) * 10) as f64;
        } else {
            solve_hyper(
                self.num_row,
                &self.u_pivot_lookup,
                &self.u_pivot_index,
                Some(&self.u_pivot_value),
                &self.u_start,
                &self.u_last_p,
                &self.u_index,
                &self.u_value,
                rhs,
            );
        }
        if self.update_method == kUpdateMethodPf {
            debug_assert!(self.update_method != kUpdateMethodPf);
            self.ftran_pf(rhs);
            rhs.tight();
            rhs.pack();
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::btranU (lines 1760-1846, 87 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    pub fn btran_u(
        &self,
        rhs: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        let _ = factor_timer_clock_pointer;
        if self.update_method == kUpdateMethodPf {
            debug_assert!(self.update_method != kUpdateMethodPf);
            self.btran_pf(rhs);
        }

        // The regular part
        //
        // Determine style of solve
        let current_density = 1.0 * (rhs.count as f64) * self.inv_num_row;
        let sparse_solve =
            rhs.count < 0 || current_density > kHyperCancel || expected_density > kHyperBtranU;
        if sparse_solve {
            // Alias to non constant
            let mut rhs_synthetic_tick: f64 = 0.0;
            // Alias to RHS
            let rhs_array = &mut rhs.array;
            let rhs_index = &mut rhs.index;
            // Alias to factor U
            let ur_start = &self.ur_start;
            let ur_end = &self.ur_lastp;
            let ur_index = &self.ur_index;
            let ur_value = &self.ur_value;
            // Local accumulation of RHS count
            let mut rhs_count: HighsInt = 0;
            // Transform
            let u_pivot_count = self.u_pivot_index.len() as HighsInt;
            // PORT NOTE (perf): per-pivot arrays sliced to `u_pivot_count` so `i_logic`
            // needs no bounds checks.
            let n = u_pivot_count as usize;
            let num_row = self.num_row as usize;
            let u_pivot_index = &self.u_pivot_index[..n];
            let u_pivot_value = &self.u_pivot_value[..n];
            let ur_start = &ur_start[..n];
            let ur_end = &ur_end[..n];
            for i_logic in 0..n {
                // Skip void
                let pivot_row = u_pivot_index[i_logic];
                if pivot_row != -1 {
                    // Normal part
                    let mut pivot_multiplier = rhs_array[pivot_row as usize];
                    if pivot_multiplier.abs() > kHighsTiny {
                        pivot_multiplier /= u_pivot_value[i_logic];
                        rhs_index[rhs_count as usize] = pivot_row;
                        rhs_count += 1;
                        rhs_array[pivot_row as usize] = pivot_multiplier;
                        let start = ur_start[i_logic];
                        let end = ur_end[i_logic];
                        if i_logic >= num_row {
                            rhs_synthetic_tick += (end - start) as f64;
                        }
                        scatter_sub(rhs_array, ur_index, ur_value, start, end, pivot_multiplier);
                    } else {
                        rhs_array[pivot_row as usize] = 0.0;
                    }
                }
            }
            // Save the count
            rhs.count = rhs_count;
            rhs.synthetic_tick +=
                rhs_synthetic_tick * 15.0 + ((u_pivot_count - self.num_row) * 10) as f64;
        } else {
            solve_hyper(
                self.num_row,
                &self.u_pivot_lookup,
                &self.u_pivot_index,
                Some(&self.u_pivot_value),
                &self.ur_start,
                &self.ur_lastp,
                &self.ur_index,
                &self.ur_value,
                rhs,
            );
        }

        // The update part
        debug_assert!(rhs.count >= 0);
        if self.update_method == kUpdateMethodFt {
            rhs.tight();
            rhs.pack();
            self.btran_ft(rhs);
            rhs.tight();
        }
        if self.update_method == kUpdateMethodMpf {
            debug_assert!(self.update_method != kUpdateMethodMpf);
            rhs.tight();
            rhs.pack();
            self.btran_mpf(rhs);
            rhs.tight();
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::ftranFT (lines 1848-1880, 33 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn ftran_ft(&self, vector: &mut HVector) {
        // Alias to non constant
        debug_assert!(vector.count >= 0);
        let mut rhs_count = vector.count;
        let rhs_index = &mut vector.index;
        let rhs_array = &mut vector.array;
        // Alias to PF buffer
        let pf_pivot_count = self.pf_pivot_index.len() as HighsInt;
        let pf_pivot_index = &self.pf_pivot_index;
        let pf_start = &self.pf_start;
        let pf_index = &self.pf_index;
        let pf_value = &self.pf_value;
        for i in 0..pf_pivot_count {
            let i_row = pf_pivot_index[i as usize];
            let value0 = rhs_array[i_row as usize];
            let mut value1 = value0;
            let start = pf_start[i as usize];
            let end = pf_start[(i + 1) as usize];
            if start < end {
                let (start, end) = (start as usize, end as usize);
                for (&j, &v) in pf_index[start..end].iter().zip(&pf_value[start..end]) {
                    value1 -= rhs_array[j as usize] * v;
                }
            }
            // This would skip the situation where they are both zeros
            if value0 != 0.0 || value1 != 0.0 {
                if value0 == 0.0 {
                    rhs_index[rhs_count as usize] = i_row;
                    rhs_count += 1;
                }
                rhs_array[i_row as usize] = if value1.abs() < kHighsTiny {
                    kHighsZero
                } else {
                    value1
                };
            }
        }
        // Save count back
        vector.count = rhs_count;
        vector.synthetic_tick +=
            (pf_pivot_count * 20 + pf_start[pf_pivot_count as usize] * 5) as f64;
        if pf_start[pf_pivot_count as usize] / (pf_pivot_count + 1) < 5 {
            vector.synthetic_tick += (pf_start[pf_pivot_count as usize] * 5) as f64;
        }
    }

    /// HiGHS: HFactor.cpp::HFactor::btranFT (lines 1882-1915, 34 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn btran_ft(&self, vector: &mut HVector) {
        // Alias to non constant
        debug_assert!(vector.count >= 0);
        let mut rhs_count = vector.count;
        let rhs_index = &mut vector.index;
        let rhs_array = &mut vector.array;
        // Alias to PF buffer
        let pf_pivot_count = self.pf_pivot_index.len() as HighsInt;
        let pf_pivot_index = &self.pf_pivot_index;
        let pf_start = &self.pf_start;
        let pf_index = &self.pf_index;
        let pf_value = &self.pf_value;
        // Apply row ETA backward
        let mut rhs_synthetic_tick: f64 = 0.0;
        let mut i = pf_pivot_count - 1;
        while i >= 0 {
            let pivot_row = pf_pivot_index[i as usize];
            let pivot_multiplier = rhs_array[pivot_row as usize];
            if pivot_multiplier != 0.0 {
                let start = pf_start[i as usize];
                let end = pf_start[(i + 1) as usize];
                rhs_synthetic_tick += (end - start) as f64;
                if start < end {
                    let (start, end) = (start as usize, end as usize);
                    for (&i_row, &v) in pf_index[start..end].iter().zip(&pf_value[start..end]) {
                        let value0 = rhs_array[i_row as usize];
                        let value1 = value0 - pivot_multiplier * v;
                        if value0 == 0.0 {
                            rhs_index[rhs_count as usize] = i_row;
                            rhs_count += 1;
                        }
                        rhs_array[i_row as usize] = if value1.abs() < kHighsTiny {
                            kHighsZero
                        } else {
                            value1
                        };
                    }
                }
            }
            i -= 1;
        }
        vector.synthetic_tick += rhs_synthetic_tick * 15.0 + (pf_pivot_count * 10) as f64;
        // Save count back
        vector.count = rhs_count;
    }

    /// HiGHS: HFactor.cpp::HFactor::ftranPF (lines 1917-1950, 34 lines)
    pub fn ftran_pf(&self, vector: &mut HVector) {
        // Alias to PF buffer
        let pf_pivot_count = self.pf_pivot_index.len() as HighsInt;
        let pf_pivot_index = &self.pf_pivot_index;
        let pf_pivot_value = &self.pf_pivot_value;
        let pf_start = &self.pf_start;
        let pf_index = &self.pf_index;
        let pf_value = &self.pf_value;

        // Alias to non constant
        let mut rhs_count = vector.count;
        let rhs_index = &mut vector.index;
        let rhs_array = &mut vector.array;

        // Forwardly
        for i in 0..pf_pivot_count {
            let pivot_row = pf_pivot_index[i as usize];
            let mut pivot_multiplier = rhs_array[pivot_row as usize];
            if pivot_multiplier.abs() > kHighsTiny {
                pivot_multiplier /= pf_pivot_value[i as usize];
                rhs_array[pivot_row as usize] = pivot_multiplier;
                let mut k = pf_start[i as usize];
                while k < pf_start[(i + 1) as usize] {
                    let index = pf_index[k as usize];
                    let value0 = rhs_array[index as usize];
                    let value1 = value0 - pivot_multiplier * pf_value[k as usize];
                    if value0 == 0.0 {
                        rhs_index[rhs_count as usize] = index;
                        rhs_count += 1;
                    }
                    rhs_array[index as usize] = if value1.abs() < kHighsTiny {
                        kHighsZero
                    } else {
                        value1
                    };
                    k += 1;
                }
            }
        }

        // Save count
        vector.count = rhs_count;
    }

    /// HiGHS: HFactor.cpp::HFactor::btranPF (lines 1952-1981, 30 lines)
    pub fn btran_pf(&self, vector: &mut HVector) {
        // Alias to PF buffer
        let pf_pivot_count = self.pf_pivot_index.len() as HighsInt;
        let pf_pivot_index = &self.pf_pivot_index;
        let pf_pivot_value = &self.pf_pivot_value;
        let pf_start = &self.pf_start;
        let pf_index = &self.pf_index;
        let pf_value = &self.pf_value;

        // Alias to non constant
        let mut rhs_count = vector.count;
        let rhs_index = &mut vector.index;
        let rhs_array = &mut vector.array;

        // Backwardly
        let mut i = pf_pivot_count - 1;
        while i >= 0 {
            let pivot_row = pf_pivot_index[i as usize];
            let mut pivot_multiplier = rhs_array[pivot_row as usize];
            let mut k = pf_start[i as usize];
            while k < pf_start[(i + 1) as usize] {
                pivot_multiplier -= pf_value[k as usize] * rhs_array[pf_index[k as usize] as usize];
                k += 1;
            }
            pivot_multiplier /= pf_pivot_value[i as usize];

            if rhs_array[pivot_row as usize] == 0.0 {
                rhs_index[rhs_count as usize] = pivot_row;
                rhs_count += 1;
            }
            rhs_array[pivot_row as usize] = if pivot_multiplier.abs() < kHighsTiny {
                1e-100
            } else {
                pivot_multiplier
            };
            i -= 1;
        }

        // Save count
        vector.count = rhs_count;
    }

    /// HiGHS: HFactor.cpp::HFactor::ftranMPF (lines 1983-1999, 17 lines)
    pub fn ftran_mpf(&self, vector: &mut HVector) {
        // Alias to non constant
        let mut rhs_count = vector.count;
        let rhs_index = &mut vector.index;
        let rhs_array = &mut vector.array;

        // Forwardly
        let pf_pivot_count = self.pf_pivot_value.len() as HighsInt;
        for i in 0..pf_pivot_count {
            solve_matrix_t(
                self.pf_start[(i * 2 + 1) as usize],
                self.pf_start[(i * 2 + 2) as usize],
                self.pf_start[(i * 2) as usize],
                self.pf_start[(i * 2 + 1) as usize],
                &self.pf_index,
                &self.pf_value,
                self.pf_pivot_value[i as usize],
                &mut rhs_count,
                rhs_index,
                rhs_array,
            );
        }

        // Remove cancellation
        vector.count = rhs_count;
    }

    /// HiGHS: HFactor.cpp::HFactor::btranMPF (lines 2001-2016, 16 lines)
    pub fn btran_mpf(&self, vector: &mut HVector) {
        // Alias to non constant
        let mut rhs_count = vector.count;
        let rhs_index = &mut vector.index;
        let rhs_array = &mut vector.array;

        // Backwardly
        let mut i = self.pf_pivot_value.len() as HighsInt - 1;
        while i >= 0 {
            solve_matrix_t(
                self.pf_start[(i * 2) as usize],
                self.pf_start[(i * 2 + 1) as usize],
                self.pf_start[(i * 2 + 1) as usize],
                self.pf_start[(i * 2 + 2) as usize],
                &self.pf_index,
                &self.pf_value,
                self.pf_pivot_value[i as usize],
                &mut rhs_count,
                rhs_index,
                rhs_array,
            );
            i -= 1;
        }

        // Remove cancellation
        vector.count = rhs_count;
    }

    /// HiGHS: HFactor.cpp::HFactor::ftranAPF (lines 2018-2034, 17 lines)
    pub fn ftran_apf(&self, vector: &mut HVector) {
        // Alias to non constant
        let mut rhs_count = vector.count;
        let rhs_index = &mut vector.index;
        let rhs_array = &mut vector.array;

        // Backwardly
        let pf_pivot_count = self.pf_pivot_value.len() as HighsInt;
        let mut i = pf_pivot_count - 1;
        while i >= 0 {
            solve_matrix_t(
                self.pf_start[(i * 2 + 1) as usize],
                self.pf_start[(i * 2 + 2) as usize],
                self.pf_start[(i * 2) as usize],
                self.pf_start[(i * 2 + 1) as usize],
                &self.pf_index,
                &self.pf_value,
                self.pf_pivot_value[i as usize],
                &mut rhs_count,
                rhs_index,
                rhs_array,
            );
            i -= 1;
        }

        // Remove cancellation
        vector.count = rhs_count;
    }

    /// HiGHS: HFactor.cpp::HFactor::btranAPF (lines 2036-2050, 15 lines)
    pub fn btran_apf(&self, vector: &mut HVector) {
        // Alias to non constant
        let mut rhs_count = vector.count;
        let rhs_index = &mut vector.index;
        let rhs_array = &mut vector.array;

        // Forwardly
        let pf_pivot_count = self.pf_pivot_value.len() as HighsInt;
        for i in 0..pf_pivot_count {
            solve_matrix_t(
                self.pf_start[(i * 2) as usize],
                self.pf_start[(i * 2 + 1) as usize],
                self.pf_start[(i * 2 + 1) as usize],
                self.pf_start[(i * 2 + 2) as usize],
                &self.pf_index,
                &self.pf_value,
                self.pf_pivot_value[i as usize],
                &mut rhs_count,
                rhs_index,
                rhs_array,
            );
        }
        vector.count = rhs_count;
    }
}
