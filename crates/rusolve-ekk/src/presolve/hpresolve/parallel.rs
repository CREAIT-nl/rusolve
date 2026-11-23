//! HiGHS source: presolve/HPresolve.cpp (lines 7267-7960)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! detectParallelRowsAndCols. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::detectParallelRowsAndCols (lines 7267-7960, 694 lines)
    ///
    /// MUST-PORT (determinism): the row/column hashes (`sparse_combine_value` with
    /// `double_hash_code(x) as u64`), the `HighsHashTable<HighsInt, HighsInt> numRowSingletons`
    /// (`index_or_insert`, `find`) and the bucket order of the `unordered_multimap` decide which
    /// parallel pair is found first. `buckets` -> `HighsMultiBuckets` (see its doc for the
    /// `find` / `last = it++` / `erase(last)` / `emplace_hint(last, ...)` pattern).
    /// `colMax`/`rowMax` are `std::vector<std::pair<double, HighsInt>>` -> `Vec<(f64, HighsInt)>`.
    pub fn detect_parallel_rows_and_cols(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::PARALLEL);
        debug_assert!(self.analysis_.allow_rule_[kPresolveRuleParallelRowsAndCols as usize]);
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleParallelRowsAndCols);
        }

        // Rust-only: the C++ `if (debug_report)` printf paths are compiled out
        // (`debug_report` is the constant `false` in HiGHS).
        let debug_report = false;

        // HiGHS 1.16 (h16_parallel_rows): the v1.16.0 gates of #3261, #3520 and #3358 below read
        // this flag once. Off, every site runs the 1.15.1 statements.
        let h16_parallel_rows = self.options.h16_parallel_rows;
        let h16_row_duals = self.options.h16_row_duals;
        // HiGHS 1.16 #3367 (h16_row_duals): 1.16 reads the member may_require_basis_postsolve_
        // (set in presolve under h16_presolve_rules); the 1.15.1 test is the
        // `#if !ENABLE_SPARSIFY_FOR_LP` condition.
        let may_require_basis_postsolve = if h16_row_duals && self.options.h16_presolve_rules {
            self.h16_may_require_basis_postsolve_
        } else {
            self.mipsolver.is_none() && self.options.lp_presolve_requires_basis_postsolve
        };

        // kMergeParallelCols etc. are C++ `constexpr HighsInt` locals.
        const kMergeParallelCols: HighsInt = 0;
        const kDominanceColToUpper: HighsInt = 1;
        const kDominanceColToLower: HighsInt = 2;
        const kDominanceDuplicateColToUpper: HighsInt = 3;
        const kDominanceDuplicateColToLower: HighsInt = 4;

        let mut row_hashes: Vec<u64> = self.rowsize.iter().map(|&x| x as u64).collect();
        let mut col_hashes: Vec<u64> = self.colsize.iter().map(|&x| x as u64).collect();
        let mut row_max: Vec<(f64, HighsInt)> = vec![(0.0, 0); self.rowsize.len()];
        let mut col_max: Vec<(f64, HighsInt)> = vec![(0.0, 0); self.colsize.len()];

        let mut num_row_singletons: HighsHashTable<HighsInt, HighsInt> = HighsHashTable::new();

        let nnz = self.a_value.len();

        // Step 1: Determine scales for rows and columns and remove column singletons
        // from the initial row hashes which are initialized with the row sizes
        for i in 0..nnz {
            if self.a_value[i] == 0.0 {
                continue;
            }
            let ar = self.a_row[i] as usize;
            let ac = self.a_col[i] as usize;
            debug_assert!(self.col_deleted[ac] == 0);
            if self.colsize[ac] == 1 {
                col_max[ac].0 = self.a_value[i];
                row_hashes[ar] = row_hashes[ar].wrapping_sub(1);
                *num_row_singletons.index_or_insert(&self.a_row[i]) += 1;
                continue;
            }
            let abs_val = self.a_value[i].abs();
            let abs_row_max = row_max[ar].0.abs();

            // among the largest values which are equal in tolerance
            // we use the nonzero with the smallest row/column index for the column/row
            // scale so that we ensure that duplicate rows/columns are scaled to have
            // the same sign
            if abs_val >= abs_row_max - self.options.small_matrix_value {
                // we are greater or equal with tolerances, check if we are either
                // strictly larger or equal with a smaller index and remember the signed
                // nonzero if one of those things is the case
                if abs_val > abs_row_max + self.options.small_matrix_value
                    || self.a_col[i] < row_max[ar].1
                {
                    row_max[ar].0 = self.a_value[i];
                    row_max[ar].1 = self.a_col[i];
                }
            }

            let abs_col_max = col_max[ac].0.abs();
            if abs_val >= abs_col_max - self.options.small_matrix_value {
                if abs_val > abs_col_max + self.options.small_matrix_value
                    || self.a_row[i] < col_max[ac].1
                {
                    col_max[ac].0 = self.a_value[i];
                    col_max[ac].1 = self.a_row[i];
                }
            }
        }

        // Step 2: Compute hash values for rows and columns excluding singleton
        // columns
        for i in 0..nnz {
            if self.a_value[i] == 0.0 {
                continue;
            }
            let ar = self.a_row[i] as usize;
            let ac = self.a_col[i] as usize;
            if self.colsize[ac] == 1 {
                col_hashes[ac] = self.a_row[i] as u64;
            } else {
                sparse_combine_value(
                    &mut row_hashes[ar],
                    self.a_col[i],
                    double_hash_code(self.a_value[i] / row_max[ar].0) as u64,
                );
                sparse_combine_value(
                    &mut col_hashes[ac],
                    self.a_row[i],
                    double_hash_code(self.a_value[i] / col_max[ac].0) as u64,
                );
            }
        }

        // Step 3: Loop over the rows and columns and put them into buckets using the
        // computed hash values. Whenever a bucket already contains a row/column,
        // check if we can apply a (nearly) parallel row reduction or a
        // parallel/dominated column reduction.
        let mut buckets = HighsMultiBuckets::default();

        for i in 0..self.model.num_col_ as usize {
            if self.col_deleted[i] != 0 {
                continue;
            }
            let ii = i as HighsInt;
            if self.colsize[i] == 0 {
                hpresolve_checked_call!(self.col_presolve(postsolve_stack, ii));
                continue;
            }
            let key = col_hashes[i];
            let count = buckets.count(key);
            // `auto it = buckets.find(key); decltype(it) last = it;`
            let mut it: usize = 0;
            let mut last: Option<usize> = if count > 0 { Some(0) } else { None };

            let mut del_col: HighsInt = -1;

            while it < count {
                let parallel_col_candidate = buckets.get(key, it);
                last = Some(it);
                it += 1;

                let pc = parallel_col_candidate as usize;

                // we want to check if the columns are parallel, first rule out
                // hash collisions with different size columns
                if self.colsize[i] != self.colsize[pc] {
                    continue;
                }
                // The columns have the same length. Next we determine whether domination
                // is possible in one of the directions, and if it is we designate the
                // dominating column as column 2. The first thing we check is whether the
                // the objective value of one of the (scaled) columns is strictly better
                // then the objective value of the other column which rules out domination
                // in one direction.

                // helpers for checking dominance between parallel columns which is
                // possible for different cases of the variable types: if col can be
                // increased infinitely in which case duplicateCol can be fixed to its
                // lower bound. duplicateCol can be decreased infinitely in which case col
                // can be fixed to its upper bound. for both cases we exploit that the
                // column that remains unfixed can always compensate for the fixed column.
                // This only holds if the compensating column can compensate exactly for
                // feasible value of the fixed column. In the continuous case this
                // trivially holds. In the case where both variables are integer and the
                // scale is +- 1 this also holds trivially. If the scale is > 1 and both
                // variables are integer, this only holds in one direction. We can apply
                // the reduction due to the following reasoning: Applying the scale to
                // col, means we change its meaning and it is not an integer variable
                // anymore, but a variable that moves on multiples of 1/scale. As we have
                // taken care that the scale is >=1 and integral for two integer
                // variables, the scaled column can always exactly compensate for the
                // other column as it can move by 1/k with k being integer. Hence every
                // kth allowed value is integral and no integral value is skipped. If the
                // compensating column is integral
                let mut check_col_impl_bounds = true;
                let mut check_duplicate_col_impl_bounds = true;

                let col: HighsInt;
                let duplicate_col: HighsInt;
                let mut col_scale: f64;
                if self.model.integrality_[i] == HighsVarType::kInteger
                    && self.model.integrality_[pc] == HighsVarType::kInteger
                {
                    // both variables are integral, hence the scale must be integral
                    // therefore first choose the smaller colMax value for col2, then check
                    // integrality of colMax[col1] / colMax[col2].
                    if col_max[i].0.abs() < col_max[pc].0.abs() {
                        col = ii;
                        duplicate_col = parallel_col_candidate;
                    } else {
                        col = parallel_col_candidate;
                        duplicate_col = ii;
                    }

                    let scale_cand = col_max[duplicate_col as usize].0 / col_max[col as usize].0;
                    // HighsUtils.h::fractionality(scale_cand, &colScale), inlined:
                    // `val = round(input); *intval = val; return abs(input - val)`.
                    col_scale = scale_cand.round();
                    let fractionality = (scale_cand - col_scale).abs();
                    if fractionality > self.options.small_matrix_value {
                        continue;
                    }
                    debug_assert!(col_scale.abs() >= 1.0);

                    // if the scale is larger than 1, duplicate column cannot compensate for
                    // all values of scaled col due to integrality as the scaled column
                    // moves on a grid of 1/scale.
                    //
                    // ToDo: Check whether this is too restrictive
                    if col_scale != 1.0 {
                        check_duplicate_col_impl_bounds = false;
                    }
                } else if self.model.integrality_[i] == HighsVarType::kInteger {
                    col = ii;
                    duplicate_col = parallel_col_candidate;
                    col_scale = col_max[duplicate_col as usize].0 / col_max[col as usize].0;

                    // as col is integral and duplicateCol is not col cannot compensate for
                    // duplicate col
                    check_col_impl_bounds = false;
                } else {
                    col = parallel_col_candidate;
                    duplicate_col = ii;
                    col_scale = col_max[duplicate_col as usize].0 / col_max[col as usize].0;

                    // as col might be integral and duplicateCol is not integral. In that
                    // case col cannot compensate for duplicate col
                    check_col_impl_bounds = self.model.integrality_[pc] != HighsVarType::kInteger;
                }

                // C++ lambdas capturing [&]: defined after the assignments above, which
                // is equivalent because they are only called below and read the
                // variables at call time.
                let col_upper_inf = || -> bool {
                    if !check_col_impl_bounds {
                        return false;
                    }
                    if self.mipsolver.is_none() {
                        // for LP we check strict redundancy of the bounds as otherwise dual
                        // postsolve might fail when the bound is used in the optimal solution
                        if col_scale > 0.0 {
                            self.is_upper_strictly_implied(col, None)
                        } else {
                            self.is_lower_strictly_implied(col, None)
                        }
                    } else {
                        // for MIP we do not need dual postsolve so the reduction is valid if
                        // the bound is weakly redundant
                        if col_scale > 0.0 {
                            self.is_upper_implied(col)
                        } else {
                            self.is_lower_implied(col)
                        }
                    }
                };

                let col_lower_inf = || -> bool {
                    if !check_col_impl_bounds {
                        return false;
                    }
                    if self.mipsolver.is_none() {
                        if col_scale > 0.0 {
                            self.is_lower_strictly_implied(col, None)
                        } else {
                            self.is_upper_strictly_implied(col, None)
                        }
                    } else {
                        // for MIP we do not need dual postsolve so the reduction is valid if
                        // the bound is weakly redundant
                        if col_scale > 0.0 {
                            self.is_lower_implied(col)
                        } else {
                            self.is_upper_implied(col)
                        }
                    }
                };

                let duplicate_col_upper_inf = || -> bool {
                    if !check_duplicate_col_impl_bounds {
                        return false;
                    }
                    if self.mipsolver.is_none() {
                        self.is_upper_strictly_implied(duplicate_col, None)
                    } else {
                        self.is_upper_implied(duplicate_col)
                    }
                };

                let duplicate_col_lower_inf = || -> bool {
                    if !check_duplicate_col_impl_bounds {
                        return false;
                    }
                    if self.mipsolver.is_none() {
                        self.is_lower_strictly_implied(duplicate_col, None)
                    } else {
                        self.is_lower_implied(duplicate_col)
                    }
                };

                // HiGHS: `double objDiff = static_cast<double>(col_cost[col] *
                // static_cast<HighsCDouble>(colScale) - col_cost[duplicateCol]);`
                let obj_diff: f64 = f64::from(
                    self.model.col_cost_[col as usize] * HighsCDouble::from(col_scale)
                        - self.model.col_cost_[duplicate_col as usize],
                );

                // if (std::abs(objDiff) > options->small_matrix_value) continue;
                let mut reduction_case = kMergeParallelCols;
                // now do the case distinctions for dominated columns
                // the cases are a lot simpler due to the helper functions
                // for checking the infinite bounds which automatically
                // incorporate the check for the variable types that allow domination.
                if obj_diff < -self.options.dual_feasibility_tolerance {
                    // scaled col is better than duplicate col
                    if col_upper_inf()
                        && self.model.col_lower_[duplicate_col as usize] != -kHighsInf
                    {
                        reduction_case = kDominanceDuplicateColToLower;
                    } else if duplicate_col_lower_inf()
                        && (col_scale < 0.0 || self.model.col_upper_[col as usize] != kHighsInf)
                        && (col_scale > 0.0 || self.model.col_lower_[col as usize] != -kHighsInf)
                    {
                        reduction_case = if col_scale > 0.0 {
                            kDominanceColToUpper
                        } else {
                            kDominanceColToLower
                        };
                    } else {
                        continue;
                    }
                } else if obj_diff > self.options.dual_feasibility_tolerance {
                    // duplicate col is better than scaled col
                    if col_lower_inf() && self.model.col_upper_[duplicate_col as usize] != kHighsInf
                    {
                        reduction_case = kDominanceDuplicateColToUpper;
                    } else if duplicate_col_upper_inf()
                        && (col_scale < 0.0 || self.model.col_lower_[col as usize] != -kHighsInf)
                        && (col_scale > 0.0 || self.model.col_upper_[col as usize] != kHighsInf)
                    {
                        reduction_case = if col_scale > 0.0 {
                            kDominanceColToLower
                        } else {
                            kDominanceColToUpper
                        };
                    } else {
                        continue;
                    }
                } else {
                    if col_upper_inf()
                        && self.model.col_lower_[duplicate_col as usize] != -kHighsInf
                    {
                        reduction_case = kDominanceDuplicateColToLower;
                    } else if col_lower_inf()
                        && self.model.col_upper_[duplicate_col as usize] != kHighsInf
                    {
                        reduction_case = kDominanceDuplicateColToUpper;
                    } else if duplicate_col_upper_inf()
                        && (col_scale < 0.0 || self.model.col_lower_[col as usize] != -kHighsInf)
                        && (col_scale > 0.0 || self.model.col_upper_[col as usize] != kHighsInf)
                    {
                        reduction_case = if col_scale > 0.0 {
                            kDominanceColToLower
                        } else {
                            kDominanceColToUpper
                        };
                    } else if duplicate_col_lower_inf()
                        && (col_scale < 0.0 || self.model.col_upper_[col as usize] != kHighsInf)
                        && (col_scale > 0.0 || self.model.col_lower_[col as usize] != -kHighsInf)
                    {
                        reduction_case = if col_scale > 0.0 {
                            kDominanceColToUpper
                        } else {
                            kDominanceColToLower
                        };
                    }
                }

                if reduction_case == kMergeParallelCols {
                    let x_int = self.model.integrality_[col as usize] == HighsVarType::kInteger;
                    let y_int =
                        self.model.integrality_[duplicate_col as usize] == HighsVarType::kInteger;
                    if x_int {
                        // The only possible reduction if the column parallelism check
                        // succeeds is to merge the two columns into one. If one column is
                        // integral this means we have restrictions on integers and need to
                        // check additional conditions to allow the merging of two integer
                        // columns, or a continuous column and an integer.
                        let illegal_scale: bool;
                        if self.model.integrality_[duplicate_col as usize] != HighsVarType::kInteger
                        {
                            debug_assert!(!y_int);
                            // only one column is integral which cannot be duplicateCol due to
                            // the way we assign the columns above
                            //
                            // Scale must not exceed 1/(y_u-y_l) in magnitude
                            illegal_scale = (col_scale
                                * (self.model.col_upper_[duplicate_col as usize]
                                    - self.model.col_lower_[duplicate_col as usize]))
                                .abs()
                                < 1.0 - self.primal_feastol;
                            if !illegal_scale && debug_report {
                                println!(
                                    "kMergeParallelCols: T-F is {} legal with scale {:.4} and duplicateCol = [{:.4}, {:.4}]",
                                    if illegal_scale { "not" } else { "   " },
                                    col_scale,
                                    self.model.col_lower_[duplicate_col as usize],
                                    self.model.col_upper_[duplicate_col as usize]
                                );
                            }
                        } else {
                            // Both columns integer
                            debug_assert!(x_int && y_int);
                            // Scale must be integer and not exceed (x_u-x_l)+1 in magnitude
                            let scale_limit = self.model.col_upper_[col as usize]
                                - self.model.col_lower_[col as usize]
                                + 1.0
                                + self.primal_feastol;
                            illegal_scale = col_scale.abs() > scale_limit;
                        }
                        if illegal_scale {
                            continue;
                        }
                    } else {
                        // Neither column integer: no problem with
                        debug_assert!(!x_int && !y_int);
                    }
                }

                let mut parallel = true;
                // now check whether the coefficients are actually parallel
                // (`for (colNz : getColumnVector(col))` with a body that calls findNonzero)
                let mut it_col = self.get_column_vector(col);
                while let Some(col_nz) = it_col.next_nz(self) {
                    let duplicate_col_row_pos = self.find_nonzero(col_nz.index(), duplicate_col);
                    parallel = duplicate_col_row_pos != -1;
                    if !parallel {
                        break;
                    }

                    let diff = f64::from(
                        self.a_value[duplicate_col_row_pos as usize]
                            - HighsCDouble::from(col_scale) * col_nz.value(),
                    );
                    parallel = diff.abs() <= self.options.small_matrix_value;
                    if !parallel {
                        break;
                    }
                }

                if !parallel {
                    continue;
                }

                // switch (reductionCase): C++ `break` inside the switch exits the switch,
                // and the unconditional `break` after it exits the candidate loop.
                // Both exit the loop here.
                match reduction_case {
                    kDominanceDuplicateColToLower => {
                        del_col = duplicate_col;
                        if self.colsize[duplicate_col as usize] == 1 {
                            let row = self.a_row[self.colhead[duplicate_col as usize] as usize];
                            *num_row_singletons.index_or_insert(&row) -= 1;
                        }
                        hpresolve_checked_call!(
                            self.fix_col_to_lower(postsolve_stack, duplicate_col)
                        );
                    }
                    kDominanceDuplicateColToUpper => {
                        del_col = duplicate_col;
                        if self.colsize[duplicate_col as usize] == 1 {
                            let row = self.a_row[self.colhead[duplicate_col as usize] as usize];
                            *num_row_singletons.index_or_insert(&row) -= 1;
                        }
                        hpresolve_checked_call!(
                            self.fix_col_to_upper(postsolve_stack, duplicate_col)
                        );
                    }
                    kDominanceColToLower => {
                        del_col = col;
                        if self.colsize[col as usize] == 1 {
                            let row = self.a_row[self.colhead[col as usize] as usize];
                            *num_row_singletons.index_or_insert(&row) -= 1;
                        }
                        hpresolve_checked_call!(self.fix_col_to_lower(postsolve_stack, col));
                    }
                    kDominanceColToUpper => {
                        del_col = col;
                        if self.colsize[col as usize] == 1 {
                            let row = self.a_row[self.colhead[col as usize] as usize];
                            *num_row_singletons.index_or_insert(&row) -= 1;
                        }
                        hpresolve_checked_call!(self.fix_col_to_upper(postsolve_stack, col));
                    }
                    kMergeParallelCols => {
                        let ok_merge = postsolve_stack.duplicate_column(
                            col_scale,
                            self.model.col_lower_[col as usize],
                            self.model.col_upper_[col as usize],
                            self.model.col_lower_[duplicate_col as usize],
                            self.model.col_upper_[duplicate_col as usize],
                            col,
                            duplicate_col,
                            self.model.integrality_[col as usize] == HighsVarType::kInteger,
                            self.model.integrality_[duplicate_col as usize]
                                == HighsVarType::kInteger,
                            self.options.mip_feasibility_tolerance,
                        );
                        // HiGHS 1.16 #3261 (h16_parallel_rows): HPresolveParallelRowsAndCols.cpp lines 364-378 (1.16.0)
                        if h16_parallel_rows {
                            if !ok_merge {
                                if debug_report {
                                    println!(
                                        "HPresolve::detectParallelRowsAndCols Illegal merge prevented"
                                    );
                                }
                                continue;
                            }
                        } else if !ok_merge && debug_report {
                            println!(
                                "HPresolve::detectParallelRowsAndCols Illegal merge prevented"
                            );
                            break;
                        }
                        // When merging a continuous variable into an integer
                        // variable, the integer will become continuous - since any
                        // value in its range can be mapped back to an integer and a
                        // continuous variable. Hence the number of integer
                        // variables in the rows corresponding to the former integer
                        // variable reduces.
                        //
                        // With the opposite - merging an integer variable into a
                        // continuous variable - the retained variable is
                        // continuous, so no action is required
                        let rowsize_int_reduction = self.model.integrality_[duplicate_col as usize]
                            != HighsVarType::kInteger
                            && self.model.integrality_[col as usize] == HighsVarType::kInteger;
                        if rowsize_int_reduction {
                            self.model.integrality_[col as usize] = HighsVarType::kContinuous;
                        }

                        self.mark_changed_col(col);
                        if self.colsize[duplicate_col as usize] == 1 {
                            let row = self.a_row[self.colhead[duplicate_col as usize] as usize];
                            *num_row_singletons.index_or_insert(&row) -= 1;
                        }

                        // compute bounds of merged variable
                        let (merge_lower, merge_upper) = if col_scale > 0.0 {
                            (
                                self.model.col_lower_[col as usize]
                                    + col_scale * self.model.col_lower_[duplicate_col as usize],
                                self.model.col_upper_[col as usize]
                                    + col_scale * self.model.col_upper_[duplicate_col as usize],
                            )
                        } else {
                            (
                                self.model.col_lower_[col as usize]
                                    + col_scale * self.model.col_upper_[duplicate_col as usize],
                                self.model.col_upper_[col as usize]
                                    + col_scale * self.model.col_lower_[duplicate_col as usize],
                            )
                        };

                        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 9407-9409 at 7b890c255
                        if self.options.mip_presolve_clique_table {
                            self.x1b_eliminate_clique_col(col);
                        }

                        // change bounds
                        hpresolve_checked_call!(self.change_col_bounds(
                            col,
                            merge_lower,
                            merge_upper
                        ));

                        // mark duplicate column as deleted
                        self.mark_col_deleted(duplicate_col);
                        // remove all nonzeros of duplicateCol
                        let mut coliter = self.colhead[duplicate_col as usize];
                        while coliter != -1 {
                            debug_assert!(self.a_col[coliter as usize] == duplicate_col);

                            let colpos = coliter;
                            let colrow = self.a_row[coliter as usize];

                            // if an integer column was merged into a continuous one make
                            // sure to update the integral rowsize
                            if rowsize_int_reduction {
                                self.rowsize_integer[colrow as usize] -= 1;
                            }

                            coliter = self.a_next[coliter as usize];

                            self.unlink(colpos);

                            self.reinsert_equation(colrow);
                        }
                        // set cost to zero
                        self.model.col_cost_[duplicate_col as usize] = 0.0;
                        del_col = duplicate_col;

                        // remove implied bounds, since they might in general not be valid
                        // anymore
                        self.reset_col_implied_bounds(col, -1);

                        // if an implicit integer and an integer column were merged, check if
                        // merged continuous column is implicit integer after merge
                        if rowsize_int_reduction
                            && self.model.integrality_[duplicate_col as usize]
                                == HighsVarType::kImplicitInteger
                        {
                            let implied_integer = self.is_implied_integer(col);
                            hpresolve_checked_call!(implied_integer.result());
                            if implied_integer.flag() {
                                self.convert_implied_integer(col, -1, true);
                            }
                        }
                        break;
                    }
                    _ => {}
                }

                break;
            }

            if del_col != -1 {
                if del_col != ii {
                    buckets.erase(key, last.expect("erase of a candidate"));
                }

                // we could have new row singletons since a column was removed. Remove
                // those rows immediately
                hpresolve_checked_call!(self.check_limits(postsolve_stack));
                hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));
            } else {
                buckets.emplace_hint(last, key, ii);
            }
        }

        buckets.clear();

        // HiGHS #3250 (mip_presolve_row_append): HPresolveParallelRowsAndCols.cpp lines
        // 474-488 (1.16.0). Iterate non-cut rows before cut rows so that when a cut is
        // parallel to a non-cut, the non-cut is already in the bucket as the surviving
        // row and the cut is removed. This prevents a cut from absorbing a non-cut
        // constraint which would cause the constraint to be lost when moveCutsToPool
        // removes it. The comparator orders by (isCutRow, index), a strict total order,
        // so the sort is the same as the C++ pdqsort.
        let mut row_order: Vec<HighsInt> = (0..self.model.num_row_).collect();
        if self.options.mip_presolve_row_append {
            row_order.sort_by_key(|&a| (postsolve_stack.is_cut_row(a), a));
        }

        for row_index in 0..self.model.num_row_ as usize {
            let i = row_order[row_index] as usize;
            if self.row_deleted[i] != 0 {
                continue;
            }
            let ii = i as HighsInt;
            // HiGHS 1.16 #3367 (h16_row_duals): HPresolveParallelRowsAndCols.cpp lines 490-494 (1.16.0)
            if self.rowsize[i] <= 1
                || (self.rowsize[i] == 2 && self.is_equation(ii))
                || (h16_row_duals && self.rowsize[i] == self.rowsize_singleton[i])
            {
                hpresolve_checked_call!(self.row_presolve(postsolve_stack, ii));
                continue;
            }
            let key = row_hashes[i];
            let count = buckets.count(key);
            // `auto it = buckets.find(rowHashes[i]); decltype(it) last = it;`
            let mut it: usize = 0;
            let mut last: Option<usize> = if count > 0 { Some(0) } else { None };

            // HiGHS 1.16 #3367 (h16_row_duals): the live per-row singleton count `rowsizeSingleton`
            // replaces the counts taken before the scan (HPresolveParallelRowsAndCols.cpp
            // lines 521-522, 1.16.0)
            let get_num_singletons = |this: &HPresolve, row: HighsInt| -> HighsInt {
                if h16_row_duals {
                    this.rowsize_singleton[row as usize]
                } else {
                    num_row_singletons.find(&row).copied().unwrap_or(0)
                }
            };

            let num_singleton = get_num_singletons(self, ii);

            // #if !ENABLE_SPARSIFY_FOR_LP (ENABLE_SPARSIFY_FOR_LP is 0: this check is active)
            // HiGHS 1.16 #3367 (h16_row_duals): HPresolveParallelRowsAndCols.cpp line 503 (1.16.0)
            if may_require_basis_postsolve && num_singleton != 0 {
                continue;
            }

            let mut del_row: HighsInt = -1;
            if count > 0 {
                self.store_row(ii);
            }
            // `while (it != buckets.end() && it->first == rowHashes[i])`: the range can only
            // shrink inside the loop under h16_row_duals, so reading its length each pass is
            // the 1.15.1 loop when the flag is off.
            while it < buckets.count(key) {
                let parallel_row_cand = buckets.get(key, it);
                last = Some(it);
                it += 1;
                let pc = parallel_row_cand as usize;

                // HiGHS 1.16 #3367 (h16_row_duals): HPresolveParallelRowsAndCols.cpp line 519 (1.16.0)
                if h16_row_duals && self.row_deleted[pc] != 0 {
                    continue;
                }

                let num_singleton = get_num_singletons(self, ii);
                let num_singleton_candidate = get_num_singletons(self, parallel_row_cand);

                // #if !ENABLE_SPARSIFY_FOR_LP (ENABLE_SPARSIFY_FOR_LP is 0: this check is active)
                if may_require_basis_postsolve && num_singleton_candidate != 0 {
                    continue;
                }

                if self.rowsize[i] - num_singleton != self.rowsize[pc] - num_singleton_candidate {
                    continue;
                }

                if num_singleton_candidate > 1 || num_singleton > 1 {
                    // we only handle the case where the rows have at most one extra
                    // singleton except when one row has no extra singleton and is an
                    // equation. In that case we sparsify the other row by adding the
                    // equation and can subsequently solve it as an individual component as
                    // it is a row which only contains singletons
                    if (num_singleton != 0 || !self.is_equation(ii))
                        && (num_singleton_candidate != 0 || !self.is_equation(parallel_row_cand))
                    {
                        continue;
                    }
                } else if num_singleton_candidate != num_singleton {
                    // if only one of the two constraints has an extra singleton,
                    // we require at least one of the constraints to be an equation
                    // if that is the case we can add that equation to the other row
                    // and will make it into either a row singleton or a doubleton equation
                    // which is removed afterwards
                    if !self.is_equation(ii) && !self.is_equation(parallel_row_cand) {
                        continue;
                    }
                }

                let row_scale = row_max[pc].0 / row_max[i].0;
                // check parallel case
                // (`for (rowNz : getStoredRow())` with a body that calls findNonzero)
                let mut parallel = true;
                let mut it_row = self.get_stored_row();
                while let Some(row_nz) = it_row.next_nz(self) {
                    if self.colsize[row_nz.index() as usize] == 1 {
                        // skip singletons
                        continue;
                    }
                    let nz_pos = self.find_nonzero(parallel_row_cand, row_nz.index());
                    parallel = nz_pos != -1;
                    if !parallel {
                        break;
                    }

                    let diff = f64::from(
                        self.a_value[nz_pos as usize]
                            - HighsCDouble::from(row_scale) * row_nz.value(),
                    );
                    parallel = diff.abs() <= self.options.small_matrix_value;
                    if !parallel {
                        break;
                    }
                }
                if !parallel {
                    continue;
                }

                if num_singleton == 0 && num_singleton_candidate == 0 {
                    let mut row_lower_tightened = false;
                    let mut row_upper_tightened = false;
                    let mut new_upper: f64;
                    let mut new_lower: f64;
                    if row_scale > 0.0 {
                        new_upper = self.model.row_upper_[i] * row_scale;
                        new_lower = self.model.row_lower_[i] * row_scale;
                    } else {
                        new_lower = self.model.row_upper_[i] * row_scale;
                        new_upper = self.model.row_lower_[i] * row_scale;
                    }

                    if new_upper < self.model.row_upper_[pc] {
                        if new_upper < self.model.row_lower_[pc] - self.primal_feastol {
                            return HPresolveResult::kPrimalInfeasible;
                        }

                        // HiGHS 1.16 #3520 (h16_parallel_rows): HPresolveParallelRowsAndCols.cpp line 585 (1.16.0)
                        let snap_upper = if h16_parallel_rows {
                            new_upper < self.model.row_lower_[pc]
                        } else {
                            new_upper <= self.model.row_lower_[pc] + self.primal_feastol
                        };
                        if snap_upper {
                            new_upper = self.model.row_lower_[pc];
                        }

                        if new_upper < self.model.row_upper_[pc] {
                            row_upper_tightened = true;
                            if row_scale > 0.0 {
                                let tmp = self.row_dual_lower[i] / row_scale;
                                self.row_dual_lower[i] = self.row_dual_lower[pc] * row_scale;
                                self.row_dual_lower[pc] = tmp;
                            } else {
                                let tmp = self.row_dual_upper[i] / row_scale;
                                self.row_dual_upper[i] = self.row_dual_lower[pc] * row_scale;
                                self.row_dual_lower[pc] = tmp;
                            }

                            if h16_row_duals {
                                // HiGHS 1.16 #3270 (h16_row_duals): HPresolveParallelRowsAndCols.cpp line 600 (1.16.0)
                                self.h16_change_row_upper(parallel_row_cand, new_upper, true);
                            } else {
                                self.model.row_upper_[pc] = new_upper;
                                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 9534 at f2ae2b320 (1.15.1 7850)
                                self.x1_row_modified(pc as HighsInt);
                            }
                        }
                    }

                    if new_lower > self.model.row_lower_[pc] {
                        if new_lower > self.model.row_upper_[pc] + self.primal_feastol {
                            return HPresolveResult::kPrimalInfeasible;
                        }

                        // HiGHS 1.16 #3520 (h16_parallel_rows): HPresolveParallelRowsAndCols.cpp line 608 (1.16.0)
                        let snap_lower = if h16_parallel_rows {
                            new_lower > self.model.row_upper_[pc]
                        } else {
                            new_lower >= self.model.row_upper_[pc] - self.primal_feastol
                        };
                        if snap_lower {
                            new_lower = self.model.row_upper_[pc];
                        }

                        if new_lower > self.model.row_lower_[pc] {
                            // the rows lower bound is tightened
                            // instead of updating the activities of dual constraints, we
                            // can simply swap the bounds on the row duals. If the old
                            // lower bound on the row dual was finite, the new row dual
                            // lower bound is infinite as the new row lower bound must be
                            // a finite value. This infinite contribution, was, however,
                            // already counted from the parallel row. Therefore by
                            // swapping the bounds unlinking the other row will not
                            // decrease the infinity counter, but simply remove a bound
                            // with zero contribution. For a negative scale we need to
                            // swap with the negated upper bound of the row dual of row i.
                            row_lower_tightened = true;
                            if row_scale > 0.0 {
                                let tmp = self.row_dual_upper[i] / row_scale;
                                self.row_dual_upper[i] = self.row_dual_upper[pc] * row_scale;
                                self.row_dual_upper[pc] = tmp;
                            } else {
                                let tmp = self.row_dual_lower[i] / row_scale;
                                self.row_dual_lower[i] = self.row_dual_upper[pc] * row_scale;
                                self.row_dual_upper[pc] = tmp;
                            }

                            if h16_row_duals {
                                // HiGHS 1.16 #3270 (h16_row_duals): HPresolveParallelRowsAndCols.cpp line 634 (1.16.0)
                                self.h16_change_row_lower(parallel_row_cand, new_lower, true);
                            } else {
                                self.model.row_lower_[pc] = new_lower;
                                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 9568 at f2ae2b320 (1.15.1 7884)
                                self.x1_row_modified(pc as HighsInt);
                            }
                        }
                    }
                    // remove implied bounds, since they might in general not be valid
                    // anymore
                    self.reset_row_dual_implied_bounds(parallel_row_cand, -1);

                    postsolve_stack.duplicate_row(
                        parallel_row_cand,
                        row_upper_tightened,
                        row_lower_tightened,
                        ii,
                        row_scale,
                    );
                    del_row = ii;
                    self.mark_row_deleted(ii);
                    for k in 0..self.rowpositions.len() {
                        let rowiter = self.rowpositions[k];
                        self.unlink(rowiter);
                    }
                    break;
                } else if self.is_equation(ii) {
                    // row i is equation and parallel (except for singletons)
                    // add to the row parallelRowCand
                    let stored = self.get_stored_row();
                    hpresolve_checked_call!(self.equality_row_addition(
                        postsolve_stack,
                        ii,
                        parallel_row_cand,
                        -row_scale,
                        stored
                    ));
                    del_row = parallel_row_cand;
                    if h16_row_duals {
                        // HiGHS 1.16 #3367 (h16_row_duals): HPresolveParallelRowsAndCols.cpp lines 661-666 (1.16.0)
                        // the candidate row has changed, remove it from the buckets
                        // (`last = buckets.erase(last)`: `last` now names the next element)
                        let erased = last.expect("erase of a candidate");
                        buckets.erase(key, erased);
                        it = erased;
                        last = Some(erased);
                        // equalityRowAddition calls rowPresolve, which can delete or change
                        // row i and overwrite rowpositions. stop if row i is deleted or only
                        // contains singleton columns
                        if self.row_deleted[i] != 0 || self.rowsize[i] == self.rowsize_singleton[i]
                        {
                            break;
                        }
                        self.store_row(ii);
                    } else if h16_parallel_rows {
                        // HiGHS 1.16 #3358 (h16_parallel_rows): the `break` added in #3358 (3358.diff lines 5-12), as merged
                        break;
                    }
                } else if self.is_equation(parallel_row_cand) {
                    // the row parallelRowCand is an equation; add it to the other row
                    let row_vector = self.get_row_vector(parallel_row_cand);
                    hpresolve_checked_call!(self.equality_row_addition(
                        postsolve_stack,
                        parallel_row_cand,
                        ii,
                        -row_max[i].0 / row_max[pc].0,
                        row_vector
                    ));
                    del_row = ii;
                    // HiGHS 1.16 #3358 (h16_parallel_rows): the `break` added in #3358 (3358.diff lines 13-20), as merged
                    if h16_parallel_rows || h16_row_duals {
                        break;
                    }
                } else {
                    debug_assert_eq!(num_singleton, 1);
                    debug_assert_eq!(num_singleton_candidate, 1);

                    let (row_lower, row_upper) = if row_scale > 0.0 {
                        (
                            self.model.row_lower_[i] * row_scale,
                            self.model.row_upper_[i] * row_scale,
                        )
                    } else {
                        (
                            self.model.row_upper_[i] * row_scale,
                            self.model.row_lower_[i] * row_scale,
                        )
                    };
                    // todo: two inequalities with one singleton. check whether the rows can
                    // be converted to equations by introducing a shared slack variable
                    // which is the case if the singletons have similar properties
                    // (objective sign, bounds, scaled coefficient) and the scaled right
                    // hand sides match. Then the case reduces to adding one equation to the
                    // other and substituting one of the singletons due to the resulting
                    // doubleton equation.
                    let _ = (row_lower, row_upper);
                }
            }

            if h16_row_duals {
                // HiGHS 1.16 #3367 (h16_row_duals): HPresolveParallelRowsAndCols.cpp lines 708-714 (1.16.0)
                // reinsert row i only if it is unchanged and contains non-singleton
                // columns. the rowPresolve called by equalityRowAddition may also have
                // deleted row i
                if del_row != ii
                    && self.row_deleted[i] == 0
                    && self.rowsize[i] != self.rowsize_singleton[i]
                {
                    buckets.emplace_hint(last, key, ii);
                }
                if del_row != -1 {
                    hpresolve_checked_call!(self.check_limits(postsolve_stack));
                }
            } else if del_row != -1 {
                if del_row != ii {
                    buckets.erase(key, last.expect("erase of a candidate"));
                }

                hpresolve_checked_call!(self.check_limits(postsolve_stack));
            } else {
                buckets.emplace_hint(last, key, ii);
            }
        }

        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleParallelRowsAndCols);
        }

        HPresolveResult::kOk
    }
}

#[cfg(test)]
mod h16_tests {
    use super::*;
    use crate::prelude::*;

    fn set_colwise(lp: &mut HighsLp, start: Vec<HighsInt>, index: Vec<HighsInt>, value: Vec<f64>) {
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = start;
        lp.a_matrix_.index_ = index;
        lp.a_matrix_.value_ = value;
    }

    /// Runs `detect_parallel_rows_and_cols` on `lp` and returns the presolve state.
    fn detect(lp: HighsLp, options: &HighsOptions) -> HPresolve {
        detect_with_stack(lp, options).0
    }

    /// As `detect`, also returning the postsolve stack.
    fn detect_with_stack(lp: HighsLp, options: &HighsOptions) -> (HPresolve, HighsPostsolveStack) {
        let num_row = lp.num_row_;
        let num_col = lp.num_col_;
        let mut presolve = HPresolve::default();
        // The rule gates (`allow_rule_`) are set up by the presolve loop; the equality-row
        // path reaches them, so set them up here as matrix.rs does.
        presolve.analysis_.setup(&lp, options, true);
        assert!(presolve.ok_set_input(lp, options, options.presolve_reduction_limit, None));
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(num_row, num_col);
        assert!(matches!(
            presolve.detect_parallel_rows_and_cols(&mut stack),
            HPresolveResult::kOk
        ));
        (presolve, stack)
    }

    /// #3358: equation row 2 (a + 2b + 3c = 5) is parallel, up to one extra singleton each,
    /// to rows 0 and 1, which sit in the same bucket. With #3358 the candidate loop breaks
    /// after the first equality-row addition, so the other row is left as it was. 1.15.1
    /// goes on and adds the equation into both rows, which row presolve then removes. Row 3 keeps
    /// columns a, b and c from being parallel.
    #[test]
    fn h16_parallel_nearly_parallel_break() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 5;
        lp.num_row_ = 4;
        lp.col_cost_ = vec![0.0; 5];
        lp.col_lower_ = vec![0.0; 5];
        lp.col_upper_ = vec![10.0; 5];
        lp.integrality_ = vec![HighsVarType::kContinuous; 5];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf, 5.0, -kHighsInf];
        lp.row_upper_ = vec![10.0, 20.0, 5.0, 100.0];
        // a: rows 0-3 = 1,1,1,1; b: 2,2,2,1; c: 3,3,3,1; s0: row 0 = 1; s1: row 1 = 1
        set_colwise(
            &mut lp,
            vec![0, 4, 8, 12, 13, 14],
            vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1],
            vec![
                1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 1.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0,
            ],
        );
        let mut options = HighsOptions::default();
        // Rows with extra singletons are skipped in LP mode unless this is off.
        options.lp_presolve_requires_basis_postsolve = false;
        options.h16_parallel_rows = false;
        let (off, off_stack) = detect_with_stack(lp.clone(), &options);
        options.h16_parallel_rows = true;
        let lp_row_upper = lp.row_upper_.clone();
        let (on, on_stack) = detect_with_stack(lp, &options);

        // On: exactly one of rows 0 and 1 (the first candidate in the bucket) receives the
        // equation and is removed; the other keeps its four nonzeros and its bound.
        assert_eq!(on.row_deleted[2], 0);
        let on_untouched: Vec<usize> = (0..2).filter(|&r| on.row_deleted[r] == 0).collect();
        assert_eq!(on_untouched.len(), 1);
        let r = on_untouched[0];
        assert_eq!(on.rowsize[r], 4);
        assert_eq!(on.model.row_upper_[r], lp_row_upper[r]);
        // Off: the equation is added into both rows, and row presolve removes both.
        assert!(off.row_deleted[0] != 0);
        assert!(off.row_deleted[1] != 0);
        assert!(on_stack.num_reductions() < off_stack.num_reductions());
    }

    /// Rows 0 and 1 are parallel with row 1 = 2 * row 0 (row scale 0.5). Row 2 and row 3 keep
    /// the columns from being parallel, so only the row pass merges anything.
    /// Row 0 is [4, 10]; row 1 has upper bound `row1_upper`.
    #[test]
    fn h16_parallel_row_bound_tolerance() {
        // Scaled upper bound of row 1 is 4 + 5e-8: within primal feastol of row 0's lower
        // bound 4. 1.15.1 snaps it to 4; 1.16.0 (#3520) keeps it, since it is not below 4.
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 3;
        lp.col_cost_ = vec![0.0; 3];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![10.0; 3];
        lp.integrality_ = vec![HighsVarType::kContinuous; 3];
        lp.row_lower_ = vec![4.0, -kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![10.0, 8.0 + 1e-7, 100.0];
        // col x: rows 0,1,2 = 1,2,1; col y: rows 0,1,2 = 1,2,3; col z: rows 0,1 = 1,2
        set_colwise(
            &mut lp,
            vec![0, 3, 6, 8],
            vec![0, 1, 2, 0, 1, 2, 0, 1],
            vec![1.0, 2.0, 1.0, 1.0, 2.0, 3.0, 1.0, 2.0],
        );
        let mut options = HighsOptions::default();
        options.h16_parallel_rows = false;
        let off = detect(lp.clone(), &options);
        options.h16_parallel_rows = true;
        let on = detect(lp, &options);

        // Off: the 1.15.1 snap gives the bound exactly 4.
        assert_eq!(off.model.row_upper_[0], 4.0);
        assert!(off.row_deleted[1] != 0);
        // On: the bound stays at the scaled value, above 4 and within feastol.
        assert!(on.model.row_upper_[0] > 4.0);
        assert!(on.model.row_upper_[0] <= 4.0 + 1e-7);
        assert!(on.row_deleted[1] != 0);
    }

    /// #3261 MIP example: maximise 2x + 2y, 2x + 2y <= 7, x integer in [0, 1], y semicontinuous
    /// in [2, 3]. With the port on the optimum is 7 (x = 1, y = 2.5). With it off the result is
    /// 6, which is what 1.15.1 returns (checked against highspy 1.15.1); the off path must keep
    /// the 1.15.1 behaviour under rule F.
    #[test]
    fn h16_parallel_pr3261_mip_optimum() {
        for (on, want) in [(true, 7.0), (false, 6.0)] {
            let mut lp = HighsLp::default();
            lp.num_col_ = 2;
            lp.num_row_ = 1;
            lp.sense_ = ObjSense::kMaximize;
            lp.col_cost_ = vec![2.0, 2.0];
            lp.col_lower_ = vec![0.0, 2.0];
            lp.col_upper_ = vec![1.0, 3.0];
            lp.row_lower_ = vec![0.0];
            lp.row_upper_ = vec![7.0];
            set_colwise(&mut lp, vec![0, 1, 2], vec![0, 0], vec![2.0, 2.0]);
            lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kSemiContinuous];
            let mut highs = Highs::default();
            assert_eq!(
                highs.set_option_value_bool("output_flag", false),
                HighsStatus::kOk
            );
            assert_eq!(
                highs.set_option_value_bool("h16_parallel_rows", on),
                HighsStatus::kOk
            );
            assert_eq!(
                highs.set_option_value_str("presolve", "on"),
                HighsStatus::kOk
            );
            assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
            assert_eq!(highs.run(), HighsStatus::kOk);
            assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
            let objective = highs.get_info().objective_function_value;
            assert!(
                (objective - want).abs() <= 1e-6,
                "on={on} objective={objective} want={want}"
            );
        }
    }

    /// Port of TestPresolveRules.cpp "test-parallel-rows-cut-ordering" (1.16.0
    /// lines 552-598). Rows 0 and 1 are parallel ([1, 1]) and row 0 is typed as a cut;
    /// with `mip_presolve_row_append` the cut is removed and the non-cut row 1 kept
    /// (the C++ test checks this after `run`, whose rule test calls
    /// detectParallelRowsAndCols directly). Off, the 1.15.1 order keeps row 0.
    #[test]
    fn x1b_psr_d_parallel_rows_cut_ordering() {
        for on in [true, false] {
            let mut lp = HighsLp::default();
            lp.num_col_ = 2;
            lp.num_row_ = 3;
            lp.col_cost_ = vec![1.0, 2.0];
            lp.col_lower_ = vec![0.0, 0.0];
            lp.col_upper_ = vec![10.0, 10.0];
            lp.row_lower_ = vec![-kHighsInf; 3];
            lp.row_upper_ = vec![5.0, 5.0, 3.0];
            set_colwise(&mut lp, vec![0, 3, 5], vec![0, 1, 2, 0, 1], vec![1.0; 5]);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.mip_presolve_row_append = on;
            let mut presolve = HPresolve::default();
            presolve.analysis_.setup(&lp, &options, true);
            assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, None));
            let mut stack = HighsPostsolveStack::default();
            stack.initialize_index_maps(3, 2);
            // Mark parallel row 0 as a cut
            stack.set_row_type(0, OrigRowType::kCut);
            assert!(matches!(
                presolve.detect_parallel_rows_and_cols(&mut stack),
                HPresolveResult::kOk
            ));
            let deleted: Vec<bool> = (0..2).map(|i| presolve.row_deleted[i] != 0).collect();
            if on {
                assert_eq!(deleted, vec![true, false], "the cut row 0 must be removed");
            } else {
                assert_eq!(deleted, vec![false, true], "1.15.1 keeps the first row");
            }
            assert!(stack.is_cut_row(0) && !stack.is_cut_row(1));
        }
    }
}
