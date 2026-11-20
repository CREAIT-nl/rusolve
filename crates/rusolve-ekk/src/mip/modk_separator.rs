//! HiGHS source: mip/HighsModkSeparator.h, mip/HighsModkSeparator.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use std::collections::HashSet;

use crate::prelude::*;
use crate::util::highs_pdqsort::pdqsort_by;

/// HiGHS: HighsModkSeparator.h::HighsModkSeparator (lines 49-58)
#[derive(Clone, Default, Debug)]
pub struct HighsModkSeparator {
    /// C++ base class `HighsSeparator`.
    pub base: HighsSeparator,
}

/// HiGHS: HighsModkSeparator.cpp::separateModKCuts (lines 25-41)
/// C++: `template <HighsInt k, typename FoundModKCut> static bool separateModKCuts(const std::vector<int64_t> & intSystemValue, const std::vector<HighsInt> & intSystemIndex, const std::vector<HighsInt> & intSystemStart, const HighsCutPool & cutpool, HighsInt numCol, FoundModKCut && foundModKCut)`
// PORT NOTE: `cutpool` is an index into `ms.mipdata_.cutpools`. The C++ lambda
// `foundModKCut` captures `cutGen`, `lpAggregator`, `transLp`, `cutpool` (all by reference)
// and `k` (by reference, reassigned between the calls in `separateLpSolution`); here it
// receives `ms` as its first argument so that it can add cuts while this function holds
// `ms` (`HighsGFkSolve::solve` reports through `&mut dyn FnMut`). `k` is the const generic.
// C++ `GFkSolve.solve<k>` uses the default `kNumRhs = 1`.
pub fn separate_mod_k_cuts<const K: u32>(
    ms: &mut HighsMipSolver,
    int_system_value: &[i64],
    int_system_index: &[HighsInt],
    int_system_start: &[HighsInt],
    cutpool: usize,
    num_col: HighsInt,
    found_mod_k_cut: &mut dyn FnMut(&mut HighsMipSolver, &[HighsGFkSolve_SolutionEntry], i32),
) -> bool {
    let mut gfk_solve = HighsGFkSolve::default();

    let num_cuts = ms.mipdata_.cutpools[cutpool].get_num_cuts();

    gfk_solve.from_csc::<K, 1>(
        int_system_value,
        int_system_index,
        int_system_start,
        num_col + 1,
    );
    gfk_solve.set_rhs::<K, 1>(num_col, 1, 0);
    gfk_solve.solve::<K, 1>(&mut |sol: &[HighsGFkSolve_SolutionEntry], w: i32| {
        found_mod_k_cut(ms, sol, w)
    });

    ms.mipdata_.cutpools[cutpool].get_num_cuts() != num_cuts
}

/// C++ `foundCut` lambda of `HighsModkSeparator::separateLpSolution` (cpp 206-243), with the
/// captured `k` passed in. C++ takes `weights` by non-const reference and sorts it in place; the
/// Rust caller hands a slice, so a sorted local copy is taken (same pdqsort, same `operator<`).
#[allow(clippy::too_many_arguments)]
fn found_mod_k_cut(
    ms: &mut HighsMipSolver,
    lp_relaxation: &HighsLpRelaxation,
    lp_aggregator: &mut HighsLpAggregator,
    cut_gen: &mut HighsCutGeneration,
    trans_lp: &mut HighsTransformedLp,
    integral_scales: &[(HighsInt, f64)],
    used_weights: &mut HashSet<Vec<(HighsInt, HighsUInt)>>,
    inds: &mut Vec<HighsInt>,
    vals: &mut Vec<f64>,
    rhs: &mut f64,
    k: HighsInt,
    weights: &[HighsGFkSolve_SolutionEntry],
) {
    // cuts which come from a single row can already be found with the
    // aggregation heuristic
    if weights.is_empty() {
        return;
    }

    let mut weights = weights.to_vec();
    pdqsort_by(&mut weights, |a, b| a.op_lt(b));
    let key: Vec<(HighsInt, HighsUInt)> = weights.iter().map(|w| (w.index, w.weight)).collect();
    if !used_weights.insert(key) {
        return;
    }

    debug_assert!(lp_aggregator.is_empty());
    for w in weights.iter() {
        // HiGHS quirk: `(w.weight * (k - 1)) % k` is HighsUInt (unsigned) arithmetic.
        let scaled = w.weight.wrapping_mul((k - 1) as HighsUInt) % (k as HighsUInt);
        let weight = integral_scales[w.index as usize].1 * ((scaled as f64) / (k as f64));
        let row = integral_scales[w.index as usize].0;
        lp_aggregator.add_row(ms, lp_relaxation, row, weight);
    }

    lp_aggregator.get_current_aggregation(ms, lp_relaxation, inds, vals, false);

    *rhs = 0.0;
    cut_gen.generate_cut(ms, lp_relaxation, trans_lp, inds, vals, rhs, true);

    if k != 2 {
        lp_aggregator.clear();
        for w in weights.iter() {
            let weight = integral_scales[w.index as usize].1 * ((w.weight as f64) / (k as f64));
            let row = integral_scales[w.index as usize].0;
            lp_aggregator.add_row(ms, lp_relaxation, row, weight);
        }
    }

    lp_aggregator.get_current_aggregation(ms, lp_relaxation, inds, vals, true);

    *rhs = 0.0;
    cut_gen.generate_cut(ms, lp_relaxation, trans_lp, inds, vals, rhs, true);

    lp_aggregator.clear();
}

impl HighsModkSeparator {
    /// HiGHS: HighsModkSeparator.cpp::HighsModkSeparator::separateLpSolution (lines 43-267)
    /// C++: `void separateLpSolution(HighsLpRelaxation & lpRelaxation, HighsLpAggregator & lpAggregator, HighsTransformedLp & transLp, HighsCutPool & cutpool)`
    // PORT NOTE: `ms` is `lpRelaxation.getMipSolver()`; `cutpool` is an index into
    // `ms.mipdata_.cutpools` (see `HighsSeparatorBox` in separator.rs).
    pub fn separate_lp_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &mut HighsLpRelaxation,
        lp_aggregator: &mut HighsLpAggregator,
        trans_lp: &mut HighsTransformedLp,
        cutpool: usize,
    ) {
        // HiGHS never changes the LP in this separator: every use below is shared.
        let lp_relaxation: &HighsLpRelaxation = &*lp_relaxation;
        let lp = lp_relaxation.get_lp();

        let mut skip_row: Vec<bool> = vec![false; lp.num_row_ as usize];

        // mark all rows that have continuous variables with a nonzero solution value
        // in the transformed LP to be skipped
        for &col in ms.mipdata_.continuous_cols.iter() {
            if trans_lp.bound_distance(col) == 0.0 {
                continue;
            }

            let start = lp.a_matrix_.start_[col as usize];
            let end = lp.a_matrix_.start_[col as usize + 1];

            for i in start..end {
                skip_row[lp.a_matrix_.index_[i as usize] as usize] = true;
            }
        }

        let mut cut_gen = HighsCutGeneration::new(ms, lp_relaxation, cutpool);

        let mut integral_scales: Vec<(HighsInt, f64)> = Vec::new();
        let mut int_system_value: Vec<i64> = Vec::new();
        let mut int_system_index: Vec<HighsInt> = Vec::new();
        let mut int_system_start: Vec<HighsInt> = Vec::new();

        int_system_start.push(0);

        let mut inds: Vec<HighsInt> = Vec::new();
        let mut vals: Vec<f64> = Vec::new();
        let mut scale_vals: Vec<f64> = Vec::new();

        let mut upper: Vec<f64> = Vec::new();
        let mut solval: Vec<f64> = Vec::new();
        let mut rhs: f64 = 0.0;

        let lp_solution = lp_relaxation.get_solution();
        let mut num_nonzero_rhs: HighsInt = 0;
        let max_int_row_len: HighsInt = (1000.0 + 0.1 * (lp.num_col_ as f64)) as HighsInt;

        for row in 0..lp.num_row_ {
            if skip_row[row as usize] {
                continue;
            }

            let leq_row: bool = if lp.row_upper_[row as usize] - lp_solution.row_value[row as usize]
                <= ms.mipdata_.feastol
            {
                true
            } else if lp_solution.row_value[row as usize] - lp.row_lower_[row as usize]
                <= ms.mipdata_.feastol
            {
                false
            } else {
                continue;
            };

            let mut rowlen: HighsInt = 0;
            let mut rowinds: &[HighsInt] = &[];
            let mut rowvals: &[f64] = &[];

            lp_relaxation.get_row(ms, row, &mut rowlen, &mut rowinds, &mut rowvals);

            if leq_row {
                rhs = lp.row_upper_[row as usize];
                inds.clear();
                inds.extend_from_slice(&rowinds[..rowlen as usize]);
                vals.clear();
                vals.extend_from_slice(&rowvals[..rowlen as usize]);
            } else {
                debug_assert!(
                    lp_solution.row_value[row as usize] - lp.row_lower_[row as usize]
                        <= ms.mipdata_.feastol
                );

                rhs = -lp.row_lower_[row as usize];
                inds.clear();
                inds.extend_from_slice(&rowinds[..rowlen as usize]);
                vals.clear();
                vals.extend(rowvals[..rowlen as usize].iter().map(|&x| -x));
            }

            let mut integral_positive = false;
            if !trans_lp.transform(
                ms,
                lp_relaxation,
                &mut vals,
                &mut upper,
                &mut solval,
                &mut inds,
                &mut rhs,
                &mut integral_positive,
                true,
                false,
            ) {
                continue;
            }

            rowlen = inds.len() as HighsInt;
            if rowlen > max_int_row_len {
                let mut int_row_len: HighsInt = 0;
                for i in 0..rowlen as usize {
                    if solval[i] <= ms.mipdata_.feastol {
                        continue;
                    }
                    if ms.is_col_continuous(inds[i]) {
                        continue;
                    }
                    int_row_len += 1;
                }

                // skip row if either too long or 0 = 0 row
                if int_row_len > max_int_row_len
                    || (int_row_len == 0 && rhs.abs() <= ms.mipdata_.epsilon)
                {
                    continue;
                }
            }

            let intscale: f64;
            let intrhs: i64;

            if !lp_relaxation.is_row_integral(ms, row) {
                scale_vals.clear();
                for i in 0..rowlen as usize {
                    if ms.is_col_continuous(inds[i]) {
                        continue;
                    }
                    if solval[i] > ms.mipdata_.feastol {
                        scale_vals.push(vals[i]);
                    }
                }

                if rhs.abs() > ms.mipdata_.epsilon {
                    scale_vals.push(-rhs);
                }
                if scale_vals.is_empty() {
                    continue;
                }

                intscale = HighsIntegers::integral_scale_vec(
                    &scale_vals,
                    ms.mipdata_.feastol,
                    ms.mipdata_.epsilon,
                );
                if intscale == 0.0 || intscale > 1e6 {
                    continue;
                }

                intrhs = HighsIntegers::nearest_integer(intscale * rhs);

                for i in 0..rowlen as usize {
                    if ms.is_col_continuous(inds[i]) {
                        continue;
                    }
                    if solval[i] > ms.mipdata_.feastol {
                        int_system_index.push(inds[i]);
                        int_system_value.push(HighsIntegers::nearest_integer(intscale * vals[i]));
                    }
                }
            } else {
                intscale = 1.0;
                intrhs = HighsIntegers::nearest_integer(rhs);

                for i in 0..rowlen as usize {
                    if solval[i] > ms.mipdata_.feastol {
                        int_system_index.push(inds[i]);
                        int_system_value.push(HighsIntegers::nearest_integer(vals[i]));
                    }
                }
            }

            num_nonzero_rhs += (intrhs != 0) as HighsInt;

            int_system_index.push(lp.num_col_);
            int_system_value.push(intrhs);
            int_system_start.push(int_system_value.len() as HighsInt);
            integral_scales.push((row, intscale));
        }

        if integral_scales.is_empty() || num_nonzero_rhs == 0 {
            return;
        }

        // PORT NOTE: C++ `HighsHashTable<std::vector<SolutionEntry>> usedWeights` has no Rust
        // counterpart (HighsHashTable keys are HighsInt or (u32, u32) only). It is used only as a
        // set (insert reports whether the key is new; cleared between the k rounds), so a HashSet
        // of the sorted (index, weight) pairs has the same behaviour.
        let mut used_weights: HashSet<Vec<(HighsInt, HighsUInt)>> = HashSet::new();
        let num_col = lp.num_col_;

        // C++ `foundCut` captures `k` by reference and `k` is reassigned between rounds; each round
        // passes its own constant (2, 3, 5, 7) into `found_mod_k_cut`.
        if separate_mod_k_cuts::<2>(
            ms,
            &int_system_value,
            &int_system_index,
            &int_system_start,
            cutpool,
            num_col,
            &mut |ms_, weights, _rhs_index| {
                found_mod_k_cut(
                    ms_,
                    lp_relaxation,
                    lp_aggregator,
                    &mut cut_gen,
                    trans_lp,
                    &integral_scales,
                    &mut used_weights,
                    &mut inds,
                    &mut vals,
                    &mut rhs,
                    2,
                    weights,
                )
            },
        ) {
            return;
        }

        used_weights.clear();
        if separate_mod_k_cuts::<3>(
            ms,
            &int_system_value,
            &int_system_index,
            &int_system_start,
            cutpool,
            num_col,
            &mut |ms_, weights, _rhs_index| {
                found_mod_k_cut(
                    ms_,
                    lp_relaxation,
                    lp_aggregator,
                    &mut cut_gen,
                    trans_lp,
                    &integral_scales,
                    &mut used_weights,
                    &mut inds,
                    &mut vals,
                    &mut rhs,
                    3,
                    weights,
                )
            },
        ) {
            return;
        }

        used_weights.clear();
        if separate_mod_k_cuts::<5>(
            ms,
            &int_system_value,
            &int_system_index,
            &int_system_start,
            cutpool,
            num_col,
            &mut |ms_, weights, _rhs_index| {
                found_mod_k_cut(
                    ms_,
                    lp_relaxation,
                    lp_aggregator,
                    &mut cut_gen,
                    trans_lp,
                    &integral_scales,
                    &mut used_weights,
                    &mut inds,
                    &mut vals,
                    &mut rhs,
                    5,
                    weights,
                )
            },
        ) {
            return;
        }

        used_weights.clear();
        separate_mod_k_cuts::<7>(
            ms,
            &int_system_value,
            &int_system_index,
            &int_system_start,
            cutpool,
            num_col,
            &mut |ms_, weights, _rhs_index| {
                found_mod_k_cut(
                    ms_,
                    lp_relaxation,
                    lp_aggregator,
                    &mut cut_gen,
                    trans_lp,
                    &integral_scales,
                    &mut used_weights,
                    &mut inds,
                    &mut vals,
                    &mut rhs,
                    7,
                    weights,
                )
            },
        );
    }

    /// HiGHS: HighsModkSeparator.h::HighsModkSeparator::HighsModkSeparator (lines 56-57)
    /// C++: `HighsModkSeparator(const HighsMipSolver & mipsolver)`
    pub fn new(ms: &HighsMipSolver) -> Self {
        HighsModkSeparator {
            base: HighsSeparator::new(ms, kModKSepaString),
        }
    }
}
