//! HiGHS source: presolve/HPresolve.cpp (lines 4987-5221), presolve/HPresolveTest.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Singleton-column stuffing and the presolve rule test. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::singletonColStuffing (lines 4987-5221, 235 lines)
    ///
    /// MUST-PORT (determinism): `pdqsort(vec.begin(), vec.end(), comp)` (C++ 5008) ->
    /// `pdqsort_by(&mut vec, comp)`; ties in the ratio comparator are ordered by pdqsort itself.
    pub fn singleton_col_stuffing(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        // singleton column stuffing
        // see Gamrath, G., Koch, T., Martin, A. et al., Progress in presolving
        // for mixed integer programming, Math. Prog. Comp. 7, 367–398 (2015).

        // count number of fixed columns
        let mut num_fixed_cols: HighsInt = 0;

        #[derive(Clone, Copy)]
        struct Candidate {
            col: HighsInt,
            val: f64,
            multiplier: HighsInt,
        }

        // C++ lambda `isSingleton`
        fn is_singleton(this: &HPresolve, col: HighsInt) -> bool {
            this.col_deleted[col as usize] == 0
                && this.colsize[col as usize] == 1
                && this.model.col_lower_[col as usize] != this.model.col_upper_[col as usize]
        }

        // C++ lambda `sortCols`
        fn sort_cols(this: &HPresolve, vec: &mut Vec<Candidate>) {
            pdqsort_by(vec.as_mut_slice(), |c1, c2| {
                this.model.col_cost_[c1.col as usize] / c1.val
                    < this.model.col_cost_[c2.col as usize] / c2.val
            });
        }

        // C++ lambda `updateActivityBounds`
        fn update_activity_bounds(
            sum_lower: &mut HighsCDouble,
            sum_upper: &mut HighsCDouble,
            sum_lower_finite: &mut bool,
            sum_upper_finite: &mut bool,
            aj: f64,
            lower_sum_bound: f64,
            upper_sum_bound: f64,
        ) {
            *sum_lower_finite = *sum_lower_finite && lower_sum_bound.abs() != kHighsInf;
            *sum_upper_finite = *sum_upper_finite && upper_sum_bound.abs() != kHighsInf;
            if *sum_lower_finite {
                *sum_lower += aj * HighsCDouble::from(lower_sum_bound);
            }
            if *sum_upper_finite {
                *sum_upper += aj * HighsCDouble::from(upper_sum_bound);
            }
        }

        // C++ lambda `addCandidate`
        fn add_candidate(
            this: &HPresolve,
            candidates: &mut Vec<Candidate>,
            col: HighsInt,
            val: f64,
            direction: HighsInt,
            min_weight: &mut f64,
            max_weight: &mut f64,
            num_integer_candidates: &mut usize,
        ) {
            if this.model.integrality_[col as usize] == HighsVarType::kInteger {
                *num_integer_candidates += 1;
            }

            *min_weight = cpp_min(*min_weight, direction as f64 * val);
            *max_weight = cpp_max(*max_weight, direction as f64 * val);
            candidates.push(Candidate {
                col,
                val,
                multiplier: direction,
            });
        }

        // C++ lambda `fixCol`
        fn fix_col(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
            col: HighsInt,
            direction: HighsInt,
        ) -> HPresolveResult {
            if direction < 0 {
                hpresolve_checked_call!(this.fix_col_to_lower(postsolve_stack, col));
            } else {
                hpresolve_checked_call!(this.fix_col_to_upper(postsolve_stack, col));
            }
            HPresolveResult::kOk
        }

        // C++ lambda `computeCandidates`
        #[allow(clippy::too_many_arguments)]
        fn compute_candidates(
            this: &HPresolve,
            row: HighsInt,
            direction: HighsInt,
            candidates: &mut Vec<Candidate>,
            sum_lower: &mut HighsCDouble,
            sum_upper: &mut HighsCDouble,
            sum_lower_finite: &mut bool,
            sum_upper_finite: &mut bool,
            num_integer_candidates: &mut usize,
            min_weight: &mut f64,
            max_weight: &mut f64,
            allow_integer_candidates: bool,
        ) -> bool {
            // vectors for candidates and activity bounds
            candidates.clear();
            candidates.reserve(this.rowsize[row as usize] as usize);
            *sum_lower = HighsCDouble::from(0.0);
            *sum_upper = HighsCDouble::from(0.0);
            *sum_lower_finite = true;
            *sum_upper_finite = true;
            *num_integer_candidates = 0;
            *min_weight = kHighsInf;
            *max_weight = -kHighsInf;

            for nz in this.get_row_vector(row).iter(this) {
                // get column index, coefficient, cost and bounds
                let j = nz.index();
                let aj = direction as f64 * nz.value();
                let cj = this.model.col_cost_[j as usize];
                let mut sum_lower_bound = this.model.col_lower_[j as usize];
                let mut sum_upper_bound = this.model.col_upper_[j as usize];
                let is_candidate = allow_integer_candidates
                    || this.model.integrality_[j as usize] != HighsVarType::kInteger;
                if is_singleton(this, j) {
                    // check singleton
                    if aj > 0.0 {
                        if cj >= 0.0 {
                            // dual fixing: fix to lower bound
                            sum_upper_bound = sum_lower_bound;
                        } else if is_candidate {
                            // candidate for stuffing
                            sum_upper_bound = sum_lower_bound;
                            add_candidate(
                                this,
                                candidates,
                                j,
                                aj,
                                1,
                                min_weight,
                                max_weight,
                                num_integer_candidates,
                            );
                        }
                    } else {
                        // aj < 0
                        debug_assert!(aj < 0.0);
                        if cj <= 0.0 {
                            // dual fixing: fix to upper bound
                            sum_lower_bound = sum_upper_bound;
                        } else if is_candidate {
                            // candidate for stuffing; multiply column with -1
                            sum_lower_bound = sum_upper_bound;
                            add_candidate(
                                this,
                                candidates,
                                j,
                                aj,
                                -1,
                                min_weight,
                                max_weight,
                                num_integer_candidates,
                            );
                        }
                    }
                }
                // update activities
                if aj < 0.0 {
                    std::mem::swap(&mut sum_lower_bound, &mut sum_upper_bound);
                }
                update_activity_bounds(
                    sum_lower,
                    sum_upper,
                    sum_lower_finite,
                    sum_upper_finite,
                    aj,
                    sum_lower_bound,
                    sum_upper_bound,
                );
                if !*sum_lower_finite && !*sum_upper_finite {
                    return false;
                }
            }
            true
        }

        // C++ lambda `checkCandidates`
        #[allow(clippy::too_many_arguments)]
        fn check_candidates(
            this: &HPresolve,
            row: HighsInt,
            direction: HighsInt,
            candidates: &mut Vec<Candidate>,
            sum_lower: &mut HighsCDouble,
            sum_upper: &mut HighsCDouble,
            sum_lower_finite: &mut bool,
            sum_upper_finite: &mut bool,
        ) -> bool {
            // number of integer candidates and weights
            let mut num_integer_candidates: usize = 0;
            let mut min_weight: f64 = 0.0;
            let mut max_weight: f64 = 0.0;

            // compute candidates
            if compute_candidates(
                this,
                row,
                direction,
                candidates,
                sum_lower,
                sum_upper,
                sum_lower_finite,
                sum_upper_finite,
                &mut num_integer_candidates,
                &mut min_weight,
                &mut max_weight,
                true,
            ) {
                // return if there are no integer columns
                if num_integer_candidates == 0 {
                    return true;
                }
                // all columns need to have same weights if we only have integer
                // columns
                if num_integer_candidates == candidates.len() {
                    return min_weight == max_weight;
                }
            }

            // recompute candidates without integer columns
            if num_integer_candidates > 0
                && compute_candidates(
                    this,
                    row,
                    direction,
                    candidates,
                    sum_lower,
                    sum_upper,
                    sum_lower_finite,
                    sum_upper_finite,
                    &mut num_integer_candidates,
                    &mut min_weight,
                    &mut max_weight,
                    false,
                )
            {
                return true;
            }

            false
        }

        // C++ lambda `checkRow`
        fn check_row(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
            num_fixed_cols: &mut HighsInt,
            row: HighsInt,
            rhs: f64,
            direction: HighsInt,
        ) -> HPresolveResult {
            // HiGHS 1.16 #3554 (h16_col_stuffing): HPresolveColStuffing.cpp lines 204-209 (1.16.0)
            let h16_col_stuffing = this.options.h16_col_stuffing;

            // skip row if rhs is not finite
            if direction as f64 * rhs == kHighsInf {
                return HPresolveResult::kOk;
            }

            // vectors for candidates and activity bounds
            let mut candidates: Vec<Candidate> = Vec::new();
            let mut sum_lower = HighsCDouble::default();
            let mut sum_upper = HighsCDouble::default();
            let mut sum_lower_finite = false;
            let mut sum_upper_finite = false;

            // compute candidates
            if !check_candidates(
                this,
                row,
                direction,
                &mut candidates,
                &mut sum_lower,
                &mut sum_upper,
                &mut sum_lower_finite,
                &mut sum_upper_finite,
            ) {
                return HPresolveResult::kOk;
            }

            // sort candidates
            sort_cols(this, &mut candidates);

            // check candidates
            for &t in candidates.iter() {
                // both bounds have to be finite
                if this.model.col_lower_[t.col as usize] == -kHighsInf
                    || this.model.col_upper_[t.col as usize] == kHighsInf
                {
                    break;
                }
                // compute delta (bound difference)
                let delta = (t.multiplier as f64 * t.val)
                    * (HighsCDouble::from(this.model.col_upper_[t.col as usize])
                        - HighsCDouble::from(this.model.col_lower_[t.col as usize]));
                // check if variable can be fixed
                if sum_upper_finite
                    && if h16_col_stuffing {
                        delta - (direction as f64 * rhs - sum_upper) <= this.primal_feastol
                    } else {
                        delta <= direction as f64 * rhs - sum_upper + this.primal_feastol
                    }
                {
                    *num_fixed_cols += 1;
                    hpresolve_checked_call!(fix_col(this, postsolve_stack, t.col, t.multiplier));
                } else if sum_lower_finite
                    && if h16_col_stuffing {
                        delta - (sum_lower - direction as f64 * rhs) <= this.primal_feastol
                    } else {
                        delta <= sum_lower - direction as f64 * rhs + this.primal_feastol
                    }
                {
                    // Previously
                    //
                    // direction * rhs <= sumLower + primal_feastol
                    //
                    // But this led to fixing at -t.multiplier until row activity
                    // was at its bound, which is primal optmal but degenerate,
                    // and does not yield an optimal basis if column costs are
                    // positive
                    *num_fixed_cols += 1;
                    hpresolve_checked_call!(fix_col(this, postsolve_stack, t.col, -t.multiplier));
                }
                // update row activities
                if sum_lower_finite {
                    sum_lower += delta;
                }
                if sum_upper_finite {
                    sum_upper += delta;
                }
            }

            HPresolveResult::kOk
        }

        // consider only non-fixed singleton columns
        if !is_singleton(self, col) {
            return HPresolveResult::kOk;
        }

        // get row index
        let row = self.a_row[self.colhead[col as usize] as usize];

        // return if we have an empty or singleton row, a large row, or the row is ranged
        // HiGHS #3404 (mip_dual_subst_limits): HPresolveColStuffing.cpp lines 235-237 at 3fa5161b1
        let x1_dsub = self.options.mip_dual_subst_limits;
        if self.rowsize[row as usize] <= 1
            || (x1_dsub && self.rowsize[row as usize] > 1000)
            || self.is_ranged(row)
        {
            return HPresolveResult::kOk;
        }

        // check row
        let row_upper = self.model.row_upper_[row as usize];
        hpresolve_checked_call!(check_row(
            self,
            postsolve_stack,
            &mut num_fixed_cols,
            row,
            row_upper,
            1
        ));
        let row_lower = self.model.row_lower_[row as usize];
        hpresolve_checked_call!(check_row(
            self,
            postsolve_stack,
            &mut num_fixed_cols,
            row,
            row_lower,
            -1
        ));

        if num_fixed_cols > 0 {
            crate::highs_log_dev!(
                &self.options.log_options,
                HighsLogType::kDetailed,
                "Singleton column stuffing fixed {} columns\n",
                num_fixed_cols
            );
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolveTest.cpp::HPresolve::presolveRuleTest (lines 12-19, 8 lines)
    pub fn presolve_rule_test(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        debug_assert!(self.options.presolve_rule_test != 0);
        if self.options.presolve_rule_test == kPresolveRuleColStuffing {
            return self.presolve_rule_test_col_stuffing(postsolve_stack);
        }
        // HiGHS #3255 (mip_fm_presolve): HPresolveTest.cpp lines 21-22 (1.16.0)
        if self.options.presolve_rule_test == kH16PresolveRuleFourierMotzkin
            && self.x1b_fm_presolve_enabled()
        {
            return self.presolve_rule_test_fourier_motzkin(postsolve_stack);
        }
        HPresolveResult::kOk
    }

    /// HiGHS #3255: Fourier-Motzkin presolve runs only with both `mip_fm_presolve` and
    /// the row-append infrastructure it needs (`mip_presolve_row_append`).
    pub(crate) fn x1b_fm_presolve_enabled(&self) -> bool {
        self.options.mip_fm_presolve && self.options.mip_presolve_row_append
    }

    /// HiGHS 1.16: HPresolveTest.cpp::HPresolve::presolveRuleTestFourierMotzkin (lines 63-70)
    pub fn presolve_rule_test_fourier_motzkin(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        debug_assert!(self.options.presolve_rule_test == kH16PresolveRuleFourierMotzkin);
        crate::highs_log_user!(
            &self.options.log_options,
            HighsLogType::kInfo,
            "HPresolve::presolveRuleTestFourierMotzkin\n"
        );
        let mut num_cols_eliminated: HighsInt = 0;
        self.fourier_motzkin(postsolve_stack, &mut num_cols_eliminated)
    }

    /// HiGHS: HPresolveTest.cpp::HPresolve::presolveRuleTestColStuffing (lines 20-38, 19 lines)
    pub fn presolve_rule_test_col_stuffing(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        debug_assert!(self.options.presolve_rule_test == kPresolveRuleColStuffing);
        debug_assert!(self.model.num_row_ == 1);
        crate::highs_log_user!(
            &self.options.log_options,
            HighsLogType::kInfo,
            "HPresolve::presolveRuleTestColStuffing\n"
        );
        let mut col: HighsInt = 0;
        while col < self.model.num_col_ {
            if self.col_deleted[col as usize] == 0 {
                let result = self.singleton_col_stuffing(postsolve_stack, col);
                if result != HPresolveResult::kOk {
                    return result;
                }
            }
            col += 1;
        }
        crate::highs_log_user!(
            &self.options.log_options,
            HighsLogType::kInfo,
            "HPresolve::presolveRuleTestColStuffing: Stuffing removed {} rows and {} columns\n",
            self.num_deleted_rows,
            self.num_deleted_cols
        );
        // Possibly remove the row
        self.row_presolve(postsolve_stack, 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// One row `sum_j x_j + s <= row_upper` with `num_bin` binaries `x_j` and one
    /// continuous singleton `s` (column `num_bin`), built by hand as in the other lib tests.
    fn stuffing_lp(num_bin: HighsInt, s_cost: f64, row_upper: f64) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = num_bin + 1;
        lp.num_row_ = 1;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0];
        for _ in 0..num_bin {
            lp.col_cost_.push(-1.0);
            lp.col_lower_.push(0.0);
            lp.col_upper_.push(1.0);
            lp.integrality_.push(HighsVarType::kInteger);
            lp.a_matrix_.index_.push(0);
            lp.a_matrix_.value_.push(1.0);
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        lp.col_cost_.push(s_cost);
        lp.col_lower_.push(0.0);
        lp.col_upper_.push(kHighsInf);
        lp.integrality_.push(HighsVarType::kContinuous);
        lp.a_matrix_.index_.push(0);
        lp.a_matrix_.value_.push(1.0);
        lp.a_matrix_
            .start_
            .push(lp.a_matrix_.index_.len() as HighsInt);
        lp.row_lower_.push(-kHighsInf);
        lp.row_upper_.push(row_upper);
        lp
    }

    /// Runs singleton column stuffing on the singleton column with the option set to `on`
    /// and returns the number of columns it deleted.
    fn stuff_and_count(num_bin: HighsInt, s_cost: f64, row_upper: f64, on: bool) -> HighsInt {
        let lp = stuffing_lp(num_bin, s_cost, row_upper);
        let mut options = HighsOptions::default();
        options.mip_dual_subst_limits = on;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, None));
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(1, num_bin + 1);
        let _ = presolve.singleton_col_stuffing(&mut stack, num_bin);
        presolve.col_deleted.iter().filter(|&&d| d != 0).count() as HighsInt
    }

    #[test]
    fn x1_dsub_row_too_large_skips() {
        // Row of 1001 entries (> 1000): with the option on, column stuffing skips the row
        // and deletes nothing; the same row with the option off is stuffed.
        let num_bin: HighsInt = 1000;
        assert_eq!(stuff_and_count(num_bin, 0.0, 1.0, true), 0);
        assert!(stuff_and_count(num_bin, 0.0, 1.0, false) > 0);
    }

    /// Row 0 is `0.1 x + s <= 2e10` with `x` in `[0, 1e11]` (not a singleton, so it
    /// contributes `fl(0.1) * 1e11` to the row's upper activity sum, which is not exact in
    /// double) and `s` a singleton in `[0, 1e10]` with cost -1 (a stuffing candidate).
    /// Row 1 is `x <= 1e11` only, so that `x` has two entries.
    fn h16_stuffing_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![0.0, -1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![1e11, 1e10];
        lp.integrality_ = vec![HighsVarType::kContinuous; 2];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![2e10, 1e11];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 3];
        lp.a_matrix_.index_ = vec![0, 1, 0];
        lp.a_matrix_.value_ = vec![0.1, 1.0, 1.0];
        lp
    }

    /// Runs singleton column stuffing on `s` (column 1) of `h16_stuffing_lp` with
    /// `h16_col_stuffing` set to `on` and returns the number of columns it deleted.
    fn h16_stuff_count(on: bool) -> HighsInt {
        let lp = h16_stuffing_lp();
        let mut options = HighsOptions::default();
        options.h16_col_stuffing = on;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, None));
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(2, 2);
        let _ = presolve.singleton_col_stuffing(&mut stack, 1);
        presolve.col_deleted.iter().filter(|&&d| d != 0).count() as HighsInt
    }

    #[test]
    fn h16_col_stuffing_on_forms_differ() {
        // Exact activity: sumUpper = fl(0.1) * 1e11 = 1e10 + 5.55e-7 and rhs - sumUpper =
        // 1e10 - 5.55e-7, so delta - (rhs - sumUpper) = 5.55e-7 > primal_feastol (1e-7).
        // The 1.16.0 form `delta - (rhs - sumUpper) <= feastol` rejects the fix to upper.
        // The 1.15.1 form `delta <= (rhs - sumUpper) + feastol` rounds 1e10 - 4.55e-7 to
        // 1e10 and accepts it, so the two forms disagree here.
        assert_eq!(h16_stuff_count(true), 0);
    }

    #[test]
    fn h16_col_stuffing_off_is_1_15_1() {
        // Same input with the option off: the 1.15.1 comparison fixes s to its upper bound.
        assert_eq!(h16_stuff_count(false), 1);
    }
}
