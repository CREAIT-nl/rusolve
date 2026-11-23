//! HiGHS source: mip/HighsTableauSeparator.h, mip/HighsTableauSeparator.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;
use crate::util::highs_hash::hash_u64_bytes;
use crate::util::highs_pdqsort::{pdqsort_branchless_by, pdqsort_by};

/// HiGHS: HighsTableauSeparator.h::HighsTableauSeparator (lines 20-32)
#[derive(Clone, Default, Debug)]
pub struct HighsTableauSeparator {
    /// C++ base class `HighsSeparator`.
    pub base: HighsSeparator,
    /// C++: `int64_t numTries`
    pub num_tries: i64,
    /// rusolve gate for mip_tableau_loosen: last LP objective, consecutive calls without gain
    pub loosen_last_obj: f64,
    pub loosen_stall: i64,
    /// rusolve mip_tableau_loosen_auto: `None` until a call tried aggregations, then whether HiGHS's
    /// density cap rejected nearly all of them
    pub loosen_starved: Option<bool>,
    pub loosen_latched_off: bool,
}

/// rusolve gate for mip_tableau_loosen: relative LP gain per call below
/// which a call counts as stalled, and the number of consecutive stalled calls that restores HiGHS's caps.
const TABLEAU_LOOSEN_MIN_GAIN: f64 = 1e-3;
const TABLEAU_LOOSEN_MAX_STALL: i64 = 2;
/// mip_tableau_loosen_auto: share of the first round's aggregations rejected by the density cap
/// (with no cut produced) that counts as starved
const TABLEAU_AUTO_STARVED_SHARE: f64 = 0.9;

/// HiGHS: HighsTableauSeparator.cpp::FractionalInteger (lines 22-37)
#[derive(Clone, Default, Debug)]
pub struct FractionalInteger {
    /// C++: `double fractionality`
    pub fractionality: f64,
    /// C++: `double row_ep_norm2`
    pub row_ep_norm2: f64,
    /// C++: `double score`
    pub score: f64,
    /// C++: `HighsInt basisIndex`
    pub basis_index: HighsInt,
    /// C++: `std::vector<std::pair<HighsInt, double>> row_ep`
    pub row_ep: Vec<(HighsInt, f64)>,
}

/// `std::set`/`std::sort` ordering = C++ `operator<` (equivalence = neither is less).
impl PartialEq for FractionalInteger {
    fn eq(&self, other: &Self) -> bool {
        !self.op_lt(other) && !other.op_lt(self)
    }
}
impl Eq for FractionalInteger {}
impl PartialOrd for FractionalInteger {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for FractionalInteger {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.op_lt(other) {
            std::cmp::Ordering::Less
        } else if other.op_lt(self) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

impl FractionalInteger {
    /// HiGHS: HighsTableauSeparator.cpp::FractionalInteger::operator< (lines 29-31)
    /// C++: `bool operator<(const FractionalInteger & other) const`
    pub fn op_lt(&self, other: &FractionalInteger) -> bool {
        self.score > other.score
    }

    /// HiGHS: HighsTableauSeparator.cpp::FractionalInteger::FractionalInteger (lines 33-33)
    /// C++: `FractionalInteger()`
    // PORT NOTE: `= default` leaves the scalars uninitialised in C++; every use assigns them
    // first, so zero is as good as any value.
    pub fn new() -> Self {
        FractionalInteger::default()
    }

    /// HiGHS: HighsTableauSeparator.cpp::FractionalInteger::FractionalInteger (lines 35-36)
    /// C++: `FractionalInteger(HighsInt basisIndex, double fractionality)`
    pub fn new_2(basis_index: HighsInt, fractionality: f64) -> Self {
        FractionalInteger {
            fractionality,
            row_ep_norm2: 0.0,
            score: -1.0,
            basis_index,
            row_ep: Vec::new(),
        }
    }
}

/// C++ `std::make_pair(a, ha) > std::make_pair(b, hb)`, which is `pair(b) < pair(a)` under
/// `std::pair::operator<` (so a NaN first element falls through to the second, as in C++).
#[allow(clippy::nonminimal_bool)]
fn pair_gt(a: f64, ha: u64, b: f64, hb: u64) -> bool {
    b < a || (!(a < b) && hb < ha)
}

/// C++ `HighsHashHelpers::hash(int64_t)`: an 8-byte key, little-endian `u32` halves.
fn hash_i64(x: i64) -> u64 {
    hash_u64_bytes(x as u32, (x >> 32) as u32)
}

impl HighsTableauSeparator {
    /// HiGHS: HighsTableauSeparator.cpp::HighsTableauSeparator::separateLpSolution (lines 39-244)
    /// C++: `void separateLpSolution(HighsLpRelaxation & lpRelaxation, HighsLpAggregator & lpAggregator, HighsTransformedLp & transLp, HighsCutPool & cutpool)`
    // PORT NOTE: `ms` is `lpRelaxation.getMipSolver()`; `cutpool` is an index into
    // `ms.mipdata_.cutpools` (see `HighsSeparatorBox` in separator.rs).
    // HiGHS quirk: the `score` filter is `!(score <= feastol)` (remove_if), not `score > feastol`,
    // so NaN scores are kept; clippy's `nonminimal_bool` is allowed for that reason.
    #[allow(clippy::nonminimal_bool)]
    pub fn separate_lp_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &mut HighsLpRelaxation,
        lp_aggregator: &mut HighsLpAggregator,
        trans_lp: &mut HighsTransformedLp,
        cutpool: usize,
    ) {
        // C++ `Highs& lpSolver = lpRelaxation.getLpSolver()` is re-fetched at each use: the Rust
        // accessor borrows `lp_relaxation` mutably, and this function changes no LP state.
        if !lp_relaxation.get_lp_solver().has_invert() {
            return;
        }

        if ms.mipdata_.cutpools[cutpool].get_num_available_cuts()
            > ms.options_mip_.mip_pool_soft_limit
        {
            return;
        }

        // PORT NOTE: C++ keeps `basisinds` as a pointer into the solver. The basis is not changed in
        // this function, so a copy holds the same values and frees the borrow of `lp_relaxation`.
        let basisinds: Vec<HighsInt> = lp_relaxation
            .get_lp_solver()
            .get_basic_variables_array()
            .to_vec();
        let num_row = lp_relaxation.num_rows();
        let num_col = lp_relaxation.num_cols();

        let mut cut_gen = HighsCutGeneration::new(ms, lp_relaxation, cutpool);
        // gate: loosened caps are used only while the LP objective still moves between calls
        // the auto lever loosens once the first round showed HiGHS's
        // density cap rejecting (nearly) every aggregated row
        let want_loosen = ms.options_mip_.mip_tableau_loosen
            || (ms.options_mip_.mip_tableau_loosen_auto && self.loosen_starved == Some(true));
        let mut loosen = want_loosen;
        if want_loosen {
            let obj = lp_relaxation.get_objective();
            let gain = obj - self.loosen_last_obj;
            let tol = TABLEAU_LOOSEN_MIN_GAIN * obj.abs().max(1.0);
            if self.loosen_last_obj.is_finite() && gain <= tol {
                self.loosen_stall += 1;
            } else {
                self.loosen_stall = 0;
            }
            self.loosen_last_obj = obj;
            if self.loosen_stall >= TABLEAU_LOOSEN_MAX_STALL {
                self.loosen_latched_off = true;
            }
            if self.loosen_latched_off {
                loosen = false;
            }
        }
        cut_gen.loosen_len = loosen;

        let mut base_row_inds: Vec<HighsInt> = Vec::new();
        let mut base_row_vals: Vec<f64> = Vec::new();

        let lp_solution = lp_relaxation.get_solution();

        let mut fractional_basisvars: Vec<FractionalInteger> = Vec::with_capacity(num_row as usize);
        for i in 0..num_row {
            let my_fractionality: f64 = if basisinds[i as usize] >= num_col {
                let row = basisinds[i as usize] - num_col;

                if !lp_relaxation.is_row_integral(ms, row) {
                    continue;
                }

                fractionality(lp_solution.row_value[row as usize], None)
            } else {
                let col = basisinds[i as usize];
                if ms.is_col_continuous(col) {
                    continue;
                }

                fractionality(lp_solution.col_value[col as usize], None)
            };

            if my_fractionality < 1000.0 * ms.mipdata_.feastol {
                continue;
            }

            fractional_basisvars.push(FractionalInteger::new_2(i, my_fractionality));
        }

        if fractional_basisvars.is_empty() {
            return;
        }
        let mut max_tries: i64 = 5000
            + (self.base.get_num_calls() as i64) * 50
            + (ms.mipdata_.total_lp_iterations - ms.mipdata_.heuristic_lp_iterations) / 10;
        if self.num_tries >= max_tries {
            return;
        }

        max_tries -= self.num_tries;

        let integral_cols_len = ms.mipdata_.integral_cols.len() as HighsInt;
        max_tries = cpp_min(
            max_tries,
            200 + ((0.1 * (cpp_min(num_row, integral_cols_len) as f64)) as i64),
        );

        if fractional_basisvars.len() as i64 > max_tries {
            let edge_wt: Option<&[f64]> = lp_relaxation.get_lp_solver_2().get_dual_edge_weights();
            let num_tries = self.num_tries;
            if let Some(edge_wt) = edge_wt {
                pdqsort_by(&mut fractional_basisvars, |fracint1, fracint2| {
                    let score1 = fracint1.fractionality * (1.0 - fracint1.fractionality)
                        / edge_wt[fracint1.basis_index as usize];
                    let score2 = fracint2.fractionality * (1.0 - fracint2.fractionality)
                        / edge_wt[fracint2.basis_index as usize];
                    pair_gt(
                        score1,
                        hash_i64(num_tries + fracint1.basis_index as i64),
                        score2,
                        hash_i64(num_tries + fracint2.basis_index as i64),
                    )
                });
            } else {
                pdqsort_by(&mut fractional_basisvars, |fracint1, fracint2| {
                    pair_gt(
                        fracint1.fractionality,
                        hash_i64(num_tries + fracint1.basis_index as i64),
                        fracint2.fractionality,
                        hash_i64(num_tries + fracint2.basis_index as i64),
                    )
                });
            }

            fractional_basisvars.truncate(max_tries as usize);
        }

        let mut row_ep_buffer = HVector::default();
        row_ep_buffer.setup(num_row);

        self.num_tries += fractional_basisvars.len() as i64;

        for fracvar in fractional_basisvars.iter_mut() {
            if lp_relaxation
                .get_lp_solver()
                .get_basis_inverse_row_sparse(fracvar.basis_index, &mut row_ep_buffer)
                != HighsStatus::kOk
            {
                continue;
            }

            // handled by other separator
            if row_ep_buffer.count == 1 {
                continue;
            }

            fracvar.row_ep_norm2 = 0.0;
            let mut min_weight = kHighsInf;
            let mut max_weight = 0.0;
            for j in 0..row_ep_buffer.count as usize {
                let row = row_ep_buffer.index[j];
                let weight = row_ep_buffer.array[row as usize];
                let max_abs_row_val = lp_relaxation.get_max_abs_row_val(ms, row);

                let scaled_weight = max_abs_row_val * weight.abs();
                if scaled_weight <= ms.mipdata_.feastol {
                    continue;
                }

                min_weight = cpp_min(min_weight, scaled_weight);
                max_weight = cpp_max(max_weight, scaled_weight);
                fracvar.row_ep_norm2 += scaled_weight * scaled_weight;
                fracvar.row_ep.push((row, weight));
            }

            if fracvar.row_ep.len() <= 1 {
                continue;
            }

            // rusolve lever mip_tableau_loosen: HiGHS caps this at 1e4
            let max_weight_range = if loosen { 1e6 } else { 1e4 };
            if max_weight / min_weight <= max_weight_range {
                fracvar.score =
                    fracvar.fractionality * (1.0 - fracvar.fractionality) / fracvar.row_ep_norm2;
            }
        }

        // erase(remove_if(score <= feastol)) == retain(!(score <= feastol))
        fractional_basisvars.retain(|f| !(f.score <= ms.mipdata_.feastol));

        if fractional_basisvars.is_empty() {
            return;
        }

        pdqsort_branchless_by(&mut fractional_basisvars, |a, b| a.op_lt(b));
        let mut best_score: f64 = -1.0;

        let num_cuts = ms.mipdata_.cutpools[cutpool].get_num_cuts();
        let best_score_fac: [f64; 2] = [0.0025, 0.01];
        let (mut tv_try, mut tv_dens) = (0usize, 0usize);

        for fracvar in fractional_basisvars.iter() {
            if ms.mipdata_.cutpools[cutpool].get_num_cuts() - num_cuts >= 1000 {
                break;
            }

            if fracvar.score
                < best_score_fac
                    [(ms.mipdata_.cutpools[cutpool].get_num_cuts() - num_cuts >= 50) as usize]
                    * best_score
            {
                break;
            }

            debug_assert!(lp_aggregator.is_empty());
            tv_try += 1;
            for row_weight in fracvar.row_ep.iter() {
                lp_aggregator.add_row(ms, lp_relaxation, row_weight.0, row_weight.1);
            }

            lp_aggregator.get_current_aggregation(
                ms,
                lp_relaxation,
                &mut base_row_inds,
                &mut base_row_vals,
                false,
            );

            // HiGHS quirk: `10 * (baseRowInds.size() - fracvar.row_ep.size())` is `size_t`
            // arithmetic, so a shorter aggregation wraps around and the test fires.
            // rusolve lever mip_tableau_loosen: no density cap (HiGHS: 10000 + numCol)
            let max_density = if loosen { usize::MAX / 20 } else { 10000usize };
            if 10usize.wrapping_mul(base_row_inds.len().wrapping_sub(fracvar.row_ep.len()))
                > max_density.wrapping_add(ms.num_col() as usize)
            {
                tv_dens += 1;
                lp_aggregator.clear();
                continue;
            }

            let len = base_row_inds.len() as HighsInt;
            if len > fracvar.row_ep.len() as HighsInt {
                let mut max_abs_val = 0.0;
                let mut min_abs_val = kHighsInf;
                for i in 0..len as usize {
                    if base_row_inds[i] < ms.num_col() {
                        max_abs_val = cpp_max(base_row_vals[i].abs(), max_abs_val);
                        min_abs_val = cpp_min(base_row_vals[i].abs(), min_abs_val);
                    }
                }
                if max_abs_val / min_abs_val > 1e6 {
                    lp_aggregator.clear();
                    continue;
                }
            }

            ms.mipdata_.debug_solution.check_row_aggregation(
                lp_relaxation.get_lp_solver_2().get_lp(),
                &base_row_inds,
                &base_row_vals,
                base_row_inds.len() as HighsInt,
            );

            let mut rhs: f64 = 0.0;
            cut_gen.generate_cut(
                ms,
                lp_relaxation,
                trans_lp,
                &mut base_row_inds,
                &mut base_row_vals,
                &mut rhs,
                false,
            );
            if trans_lp.get_globaldom(ms).infeasible() {
                break;
            }

            lp_aggregator.get_current_aggregation(
                ms,
                lp_relaxation,
                &mut base_row_inds,
                &mut base_row_vals,
                true,
            );
            rhs = 0.0;
            cut_gen.generate_cut(
                ms,
                lp_relaxation,
                trans_lp,
                &mut base_row_inds,
                &mut base_row_vals,
                &mut rhs,
                false,
            );
            if trans_lp.get_globaldom(ms).infeasible() {
                break;
            }

            lp_aggregator.clear();
            if best_score == -1.0 && ms.mipdata_.cutpools[cutpool].get_num_cuts() != num_cuts {
                best_score = fracvar.score;
            }
        }
        if ms.options_mip_.mip_tableau_loosen_auto
            && self.loosen_starved.is_none()
            && tv_try > 0
            && !want_loosen
        {
            let made = ms.mipdata_.cutpools[cutpool].get_num_cuts() - num_cuts;
            self.loosen_starved =
                Some(made == 0 && (tv_dens as f64) >= TABLEAU_AUTO_STARVED_SHARE * (tv_try as f64));
        }
    }

    /// HiGHS: HighsTableauSeparator.h::HighsTableauSeparator::HighsTableauSeparator (lines 30-31)
    /// C++: `HighsTableauSeparator(const HighsMipSolver & mipsolver)`
    pub fn new(ms: &HighsMipSolver) -> Self {
        HighsTableauSeparator {
            base: HighsSeparator::new(ms, kTableauSepaString),
            num_tries: 0,
            loosen_last_obj: f64::NEG_INFINITY,
            loosen_stall: 0,
            loosen_starved: None,
            loosen_latched_off: false,
        }
    }
}
