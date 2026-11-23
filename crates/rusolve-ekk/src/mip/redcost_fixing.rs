//! HiGHS source: mip/HighsRedcostFixing.h, mip/HighsRedcostFixing.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS latest (#3474, #3481, 3771ccf96): HighsRedcostFixing.cpp::HighsRedcostFixing::addRootRedcost,
/// lambda `findLurkingBounds` (lines 236-304).
///
/// `lurking_bounds` and `other_lurking_bounds` are the flattened multimaps (see
/// [`multimap_entries`]). The C++ lambda reads `mipdata_->feastol` and `lower_bound`, which
/// are passed as values.
// `lurkingBounds.emplace_hint(pos, ..)` with `pos = lower_bound(key)`: libc++
// `__tree::__find_leaf(hint, ..)` (SDK usr/include/c++/v1/__tree) sees `key <= *hint` and
// `*prev(hint) < key`, so it links the node immediately before `hint`, i.e. before every entry
// with an equal key. The flattened `insert(pos, ..)` below is that placement.
#[allow(clippy::too_many_arguments)]
fn find_lurking_bounds_x1(
    feastol: f64,
    lower_bound: f64,
    direction: HighsInt,
    input_bound: f64,
    input_other_bound: f64,
    lp_objective: f64,
    red_cost: f64,
    max_num_steps: HighsInt,
    max_num_steps_exp: HighsInt,
    lurking_bounds: &mut Vec<(f64, HighsInt)>,
    other_lurking_bounds: &mut Vec<(f64, HighsInt)>,
) {
    // bounds whose absolute value is at least this are not handled: they
    // cannot be cast to HighsInt (kHighsIInf) or, for 64-bit HighsInt, leave
    // no room for the integer arithmetic in findLurkingBounds and are not
    // exact integers in double precision (kHighsInvMacheps)
    // kHighsInvMacheps = 1/kHighsMacheps = 1/ldexp(1,-52) = 2^52
    const K_HIGHS_INV_MACHEPS: f64 = 4503599627370496.0;
    let max_bound = cpp_min(kHighsIInf as f64, K_HIGHS_INV_MACHEPS);

    // bounds whose absolute value is at least maxBound cannot be handled
    if input_bound.abs() >= max_bound {
        return;
    }
    let bound: HighsInt = input_bound as HighsInt;

    if (direction as f64) * red_cost == kHighsInf {
        lurking_bounds.clear();
        other_lurking_bounds.clear();
        lurking_bounds.push((-kHighsInf, bound));
        return;
    }

    let last_bound: i64 = if (direction as f64) * input_other_bound >= max_bound {
        // other bound cannot be handled, treat it as infinite
        let mut lb = i64::from(bound) + (direction as i64) * (max_num_steps as i64);
        // keep the lurking bounds in the range of HighsInt
        lb = (direction as i64) * cpp_min((direction as i64) * lb, i64::from(kHighsIInf));
        lb
    } else {
        input_other_bound as i64 - (direction as i64)
    };

    let mut step: i64 = 1;
    let range: i64 = (direction as i64) * (last_bound - i64::from(bound));
    if range > (max_num_steps as i64) {
        step = (range + (max_num_steps as i64) - 1) >> (max_num_steps_exp as u32);
    }
    let shift: f64 = (direction as f64) * (1.0 - 10.0 * feastol);
    step *= direction as i64;

    // C++ `for (int64_t lurkingBound = bound; direction * lurkingBound <= direction * lastBound;
    // lurkingBound += step)`. `continue` in C++ still runs the increment, which the labelled
    // block below reproduces by leaving the block and then stepping.
    let mut lurking_bound: i64 = i64::from(bound);
    while (direction as i64) * lurking_bound <= (direction as i64) * last_bound {
        'body: {
            let frac_bound: f64 = (lurking_bound - i64::from(bound)) as f64 + shift;
            let required_cutoff_bound: f64 = frac_bound * red_cost + lp_objective;
            if required_cutoff_bound < lower_bound + feastol {
                break 'body;
            }

            // check if we already have a better lurking bound stored
            let mut useful = true;
            let pos = lurking_bounds.partition_point(|entry| entry.0 < required_cutoff_bound);
            for entry in &lurking_bounds[pos..] {
                useful = (direction as i64) * i64::from(entry.1)
                    >= (direction as i64) * (lurking_bound + step);
                if !useful {
                    break;
                }
            }
            if !useful {
                break 'body;
            }

            // we have no better lurking bound stored; store this lurking bound and check if it
            // dominates one that is already stored
            lurking_bounds.insert(pos, (required_cutoff_bound, lurking_bound as HighsInt));
            let mut it = pos;
            let mut i: usize = 0;
            while i != it {
                if (direction as i64) * i64::from(lurking_bounds[i].1)
                    >= (direction as i64) * lurking_bound
                {
                    lurking_bounds.remove(i);
                    it -= 1;
                } else {
                    i += 1;
                }
            }
        }
        lurking_bound += step;
    }
}

/// HiGHS: HighsRedcostFixing.h::HighsRedcostFixing (lines 25-43)
#[derive(Clone, Default, Debug)]
pub struct HighsRedcostFixing {
    /// C++: `std::vector<std::multimap<double, HighsInt>> lurkingColUpper`
    pub lurking_col_upper: Vec<CppMultimapF64>,
    /// C++: `std::vector<std::multimap<double, HighsInt>> lurkingColLower`
    pub lurking_col_lower: Vec<CppMultimapF64>,
}

/// C++ `static_cast<HighsInt>(double)` as the oracle build (Apple clang, arm64) executes it:
/// `fcvtzs` saturates, so `+inf` and values above `INT_MAX` give `INT_MAX`, `-inf` and values
/// below `INT_MIN` give `INT_MIN`, and NaN gives 0. That is exactly Rust's `as i32`. Used for
/// the column bounds of `addRootRedcost`, which are `+/-inf` for unbounded columns.
// HiGHS quirk: C++ leaves an out-of-range conversion undefined; the port reproduces arm64
// (checked with a clang -O2 probe), not the x86 `cvttsd2si` INT_MIN result.
fn cpp_double_to_int(x: f64) -> HighsInt {
    x as HighsInt
}

/// Flattens a `std::multimap<double, HighsInt>` in iteration order (key ascending, then the
/// values of one key in their stored order).
fn multimap_entries(map: &CppMultimapF64) -> Vec<(f64, HighsInt)> {
    map.iter()
        .flat_map(|(key, values)| values.iter().map(move |value| (key.0, *value)))
        .collect()
}

/// Rebuilds a multimap from entries in iteration order; equal keys keep their order.
fn multimap_from_entries(entries: &[(f64, HighsInt)]) -> CppMultimapF64 {
    let mut map = CppMultimapF64::new();
    for &(key, value) in entries {
        map.entry(OrderedF64(key)).or_default().push(value);
    }
    map
}

/// HiGHS: HighsRedcostFixing.cpp::HighsRedcostFixing::addRootRedcost, lambda
/// `findLurkingBounds` (lines 233-302).
///
/// `lurking_bounds` and `other_lurking_bounds` are the flattened multimaps (see
/// [`multimap_entries`]). The C++ lambda reads `mipdata_->feastol` and `lower_bound`, which
/// are passed as values.
// `lurkingBounds.emplace_hint(pos, ..)` with `pos = lower_bound(key)`: libc++
// `__tree::__find_leaf(hint, ..)` (SDK usr/include/c++/v1/__tree) sees `key <= *hint` and
// `*prev(hint) < key`, so it links the node immediately before `hint`, i.e. before every entry
// with an equal key. The flattened `insert(pos, ..)` below is that placement (pinned by the
// unit test `emplace_hint_places_before_equal_keys`, checked with a clang probe).
#[allow(clippy::too_many_arguments)]
fn find_lurking_bounds(
    feastol: f64,
    lower_bound: f64,
    direction: HighsInt,
    bound: HighsInt,
    other_bound: HighsInt,
    is_other_bound_finite: bool,
    lp_objective: f64,
    red_cost: f64,
    max_num_steps: HighsInt,
    max_num_steps_exp: HighsInt,
    lurking_bounds: &mut Vec<(f64, HighsInt)>,
    other_lurking_bounds: &mut Vec<(f64, HighsInt)>,
) {
    // HiGHS quirk: the integer arithmetic below can overflow for INT_MIN / INT_MAX bounds
    // (C++ undefined behaviour); wrapping ops give the two's-complement result of the
    // arm64 oracle build. `max_num_steps_exp` is in [5, 10], so the shift never reaches 32.
    if (direction as f64) * red_cost == kHighsInf {
        lurking_bounds.clear();
        other_lurking_bounds.clear();
        lurking_bounds.push((-kHighsInf, bound));
        return;
    }

    let last_bound: HighsInt = if !is_other_bound_finite {
        bound.wrapping_add(direction.wrapping_mul(max_num_steps))
    } else {
        other_bound.wrapping_sub(direction)
    };

    let mut step: HighsInt = 1;
    let range: HighsInt = direction.wrapping_mul(last_bound.wrapping_sub(bound));
    if range > max_num_steps {
        step = range.wrapping_add(max_num_steps).wrapping_sub(1) >> max_num_steps_exp;
    }
    let shift: f64 = (direction as f64) * (1.0 - 10.0 * feastol);
    step = step.wrapping_mul(direction);

    // C++ `for (lurkingBound = bound; direction * lurkingBound <= direction * lastBound;
    // lurkingBound += step)`. `continue` in C++ still runs the increment, which the labelled
    // block below reproduces by leaving the block and then stepping.
    let mut lurking_bound: HighsInt = bound;
    while direction.wrapping_mul(lurking_bound) <= direction.wrapping_mul(last_bound) {
        'body: {
            let frac_bound: f64 = f64::from(lurking_bound.wrapping_sub(bound)) + shift;
            let required_cutoff_bound: f64 = frac_bound * red_cost + lp_objective;
            if required_cutoff_bound < lower_bound + feastol {
                break 'body;
            }

            // check if we already have a better lurking bound stored
            let mut useful = true;
            let pos = lurking_bounds.partition_point(|entry| entry.0 < required_cutoff_bound);
            for entry in &lurking_bounds[pos..] {
                useful = direction.wrapping_mul(entry.1)
                    >= direction.wrapping_mul(lurking_bound.wrapping_add(step));
                if !useful {
                    break;
                }
            }
            if !useful {
                break 'body;
            }

            // we have no better lurking bound stored; store this lurking bound and check if it
            // dominates one that is already stored
            lurking_bounds.insert(pos, (required_cutoff_bound, lurking_bound));
            let mut it = pos;
            let mut i: usize = 0;
            while i != it {
                if direction.wrapping_mul(lurking_bounds[i].1)
                    >= direction.wrapping_mul(lurking_bound)
                {
                    lurking_bounds.remove(i);
                    it -= 1;
                } else {
                    i += 1;
                }
            }
        }
        lurking_bound = lurking_bound.wrapping_add(step);
    }
}

impl HighsRedcostFixing {
    /// HiGHS: HighsRedcostFixing.cpp::HighsRedcostFixing::getLurkingBounds (lines 12-35)
    /// C++: `std::vector<std::pair<double, HighsDomainChange>> getLurkingBounds(const HighsMipSolver & mipsolver, const HighsDomain & globaldom) const`
    pub fn get_lurking_bounds(
        &self,
        ms: &HighsMipSolver,
        globaldom: &HighsDomain,
    ) -> Vec<(f64, HighsDomainChange)> {
        let mut domchgs: Vec<(f64, HighsDomainChange)> = Vec::new();
        if self.lurking_col_lower.is_empty() {
            return domchgs;
        }

        for &col in ms.mipdata_.integral_cols.iter() {
            let ci = col as usize;
            for (key, values) in self.lurking_col_lower[ci].iter() {
                for &value in values {
                    if (value as f64) > globaldom.col_lower_[ci] {
                        domchgs.push((
                            key.0,
                            HighsDomainChange {
                                boundval: value as f64,
                                column: col,
                                boundtype: HighsBoundType::kLower,
                            },
                        ));
                    }
                }
            }

            for (key, values) in self.lurking_col_upper[ci].iter() {
                for &value in values {
                    if (value as f64) < globaldom.col_upper_[ci] {
                        domchgs.push((
                            key.0,
                            HighsDomainChange {
                                boundval: value as f64,
                                column: col,
                                boundtype: HighsBoundType::kUpper,
                            },
                        ));
                    }
                }
            }
        }

        domchgs
    }

    /// HiGHS: HighsRedcostFixing.cpp::HighsRedcostFixing::propagateRootRedcost (lines 37-72)
    /// C++: `void propagateRootRedcost(const HighsMipSolver & mipsolver)`
    // PORT NOTE: `mipsolver.mipdata_->getDomain()` is `domains[0]` (attached; the body
    // detaches it around `changeBound` / `propagate`, which take `ms`). `self` is detached from
    // `ms.mipdata_.redcostfixing` by the caller.
    pub fn propagate_root_redcost(&mut self, ms: &mut HighsMipSolver) {
        if self.lurking_col_lower.is_empty() {
            return;
        }
        // Reporting only: the root domain's stack growth across this call.
        let n0 = ms.mipdata_.domains[0].domchgstack_.len();
        self.propagate_root_redcost_body(ms);
        crate::work::add_event(crate::work::Event::RedcostRootRuns, 1);
        crate::work::add_event(
            crate::work::Event::RedcostRootBoundChanges,
            ms.mipdata_.domains[0].domchgstack_.len().saturating_sub(n0) as u64,
        );
    }

    fn propagate_root_redcost_body(&mut self, ms: &mut HighsMipSolver) {
        if self.lurking_col_lower.is_empty() {
            return;
        }

        let num_integral = ms.mipdata_.integral_cols.len();
        for idx in 0..num_integral {
            let col = ms.mipdata_.integral_cols[idx];
            let ci = col as usize;
            let lower_bound = ms.mipdata_.lower_bound;
            let upper_limit = ms.mipdata_.upper_limit;

            // HiGHS quirk: erase(begin, upper_bound(lower_bound)) drops every key <= lower_bound.
            self.lurking_col_lower[ci].retain(|key, _| lower_bound < key.0);
            self.lurking_col_upper[ci].retain(|key, _| lower_bound < key.0);

            // The C++ iterators are not invalidated by `changeBound` (the multimaps are only
            // touched by this worker's redcost fixing, which is detached here), so the values
            // are collected first.
            let lower_values: Vec<HighsInt> = self.lurking_col_lower[ci]
                .range(OrderedF64(upper_limit)..)
                .flat_map(|(_, values)| values.iter().copied())
                .collect();
            for value in lower_values {
                if (value as f64) > ms.mipdata_.domains[0].col_lower_[ci] {
                    detach!(ms.mipdata_.domains[0], |d| d.change_bound_2(
                        ms,
                        HighsBoundType::kLower,
                        col,
                        value as f64,
                        HighsDomain_Reason::unspecified(),
                    ));
                    if ms.mipdata_.domains[0].infeasible() {
                        return;
                    }
                }
            }

            let upper_values: Vec<HighsInt> = self.lurking_col_upper[ci]
                .range(OrderedF64(upper_limit)..)
                .flat_map(|(_, values)| values.iter().copied())
                .collect();
            for value in upper_values {
                if (value as f64) < ms.mipdata_.domains[0].col_upper_[ci] {
                    detach!(ms.mipdata_.domains[0], |d| d.change_bound_2(
                        ms,
                        HighsBoundType::kUpper,
                        col,
                        value as f64,
                        HighsDomain_Reason::unspecified(),
                    ));
                    if ms.mipdata_.domains[0].infeasible() {
                        return;
                    }
                }
            }
        }

        detach!(ms.mipdata_.domains[0], |d| d.propagate(ms));
    }

    /// HiGHS: HighsRedcostFixing.cpp::HighsRedcostFixing::propagateRedCost (lines 74-193)
    /// C++: `static void propagateRedCost(const HighsMipSolver & mipsolver, HighsDomain & localdomain, HighsDomain & globaldom, const HighsLpRelaxation & lp, HighsConflictPool & conflictpool, HighsPseudocost & pseudocost, double upper_limit)`
    // PORT NOTE: `globaldom` / `conflictpool` are indices into `ms.mipdata_.domains` /
    // `conflictpools` (the HighsDomain::conflictAnalyzeReconvergence convention); `localdomain`,
    // `lp` and `pseudocost` are detached by the caller (`HighsSearch`).
    pub fn propagate_red_cost(
        ms: &mut HighsMipSolver,
        localdomain: &mut HighsDomain,
        globaldom: usize,
        lp: &HighsLpRelaxation,
        conflictpool: usize,
        pseudocost: &mut HighsPseudocost,
        upper_limit: f64,
    ) {
        // Reporting only: tally calls, bound changes and cutoffs.
        let n0 = localdomain.domchgstack_.len();
        Self::propagate_red_cost_body(
            ms,
            localdomain,
            globaldom,
            lp,
            conflictpool,
            pseudocost,
            upper_limit,
        );
        crate::work::add_event(crate::work::Event::RedcostLocalRuns, 1);
        crate::work::add_event(
            crate::work::Event::RedcostLocalBoundChanges,
            localdomain.domchgstack_.len().saturating_sub(n0) as u64,
        );
        if localdomain.infeasible() {
            crate::work::add_event(crate::work::Event::RedcostLocalCutoffs, 1);
        }
    }

    fn propagate_red_cost_body(
        ms: &mut HighsMipSolver,
        localdomain: &mut HighsDomain,
        globaldom: usize,
        lp: &HighsLpRelaxation,
        conflictpool: usize,
        pseudocost: &mut HighsPseudocost,
        upper_limit: f64,
    ) {
        let lpredcost: &[f64] = &lp.get_solution().col_dual;
        let lpobjective: f64 = lp.get_objective();
        let gap: HighsCDouble = HighsCDouble::from(upper_limit) - lpobjective;

        let tolerance: f64 = cpp_max(
            10.0 * ms.mipdata_.feastol,
            ms.mipdata_.epsilon * f64::from(gap),
        );

        debug_assert!(!localdomain.infeasible());
        let mut bound_changes: Vec<HighsDomainChange> =
            Vec::with_capacity(ms.mipdata_.integral_cols.len());
        for idx in 0..ms.mipdata_.integral_cols.len() {
            let col = ms.mipdata_.integral_cols[idx];
            let ci = col as usize;
            // lpobj + (col - bnd) * redcost <= cutoffbound
            // (col - bnd) * redcost <= gap
            // redcost * col <= gap + redcost * bnd
            //   case bnd is upper bound  => redcost < 0 :
            //      col >= gap/redcost + ub
            //      <=> redcost < gap / (lb - ub)
            //   case bnd is lower bound  => redcost > 0 :
            //      col <= gap/redcost + lb
            //      <=> redcost > gap / (ub - lb)
            if localdomain.col_upper_[ci] == localdomain.col_lower_[ci] {
                continue;
            }
            if lpredcost[ci].abs() <= tolerance {
                continue;
            }

            let max_increase: f64 =
                lpredcost[ci] * (localdomain.col_upper_[ci] - localdomain.col_lower_[ci]);
            if max_increase > gap {
                debug_assert!(localdomain.col_lower_[ci] != -kHighsInf);
                debug_assert!(lpredcost[ci] > tolerance);
                // HiGHS quirk: `floor` is applied to the compensated sum before `double(..)`.
                let newub: f64 = f64::from(
                    (gap / lpredcost[ci] + localdomain.col_lower_[ci] + ms.mipdata_.feastol)
                        .floor(),
                );
                if newub >= localdomain.col_upper_[ci] {
                    continue;
                }
                debug_assert!(newub < localdomain.col_upper_[ci]);

                if ms.mipdata_.domains[globaldom].is_binary(ms, col) {
                    bound_changes.push(HighsDomainChange {
                        boundval: newub,
                        column: col,
                        boundtype: HighsBoundType::kUpper,
                    });
                } else {
                    localdomain.change_bound_2(
                        ms,
                        HighsBoundType::kUpper,
                        col,
                        newub,
                        HighsDomain_Reason::unspecified(),
                    );
                    if localdomain.infeasible() {
                        return;
                    }
                }
            } else if max_increase < -gap {
                debug_assert!(localdomain.col_upper_[ci] != kHighsInf);
                debug_assert!(lpredcost[ci] < -tolerance);
                // HiGHS quirk: `ceil` is applied to the compensated sum before `double(..)`.
                let newlb: f64 = f64::from(
                    (gap / lpredcost[ci] + localdomain.col_upper_[ci] - ms.mipdata_.feastol).ceil(),
                );

                if newlb <= localdomain.col_lower_[ci] {
                    continue;
                }
                debug_assert!(newlb > localdomain.col_lower_[ci]);

                if ms.mipdata_.domains[globaldom].is_binary(ms, col) {
                    bound_changes.push(HighsDomainChange {
                        boundval: newlb,
                        column: col,
                        boundtype: HighsBoundType::kLower,
                    });
                } else {
                    localdomain.change_bound_2(
                        ms,
                        HighsBoundType::kLower,
                        col,
                        newlb,
                        HighsDomain_Reason::unspecified(),
                    );
                    if localdomain.infeasible() {
                        return;
                    }
                }
            }
        }

        if !bound_changes.is_empty() {
            let mut inds: Vec<HighsInt> = Vec::new();
            let mut vals: Vec<f64> = Vec::new();
            let mut rhs: f64 = 0.0;

            if (ms.options_mip_.mip_redcost_proof_uncapped || bound_changes.len() <= 100)
                && lp.compute_dual_proof(
                    ms,
                    globaldom,
                    upper_limit,
                    &mut inds,
                    &mut vals,
                    &mut rhs,
                    false,
                )
            {
                let old_num_conflicts = ms.mipdata_.conflictpools[conflictpool].get_num_conflicts();
                for domchg in bound_changes.iter() {
                    if localdomain.is_active(domchg) {
                        continue;
                    }
                    let len = inds.len() as HighsInt;
                    localdomain.conflict_analyze_reconvergence(
                        ms,
                        domchg,
                        &inds,
                        &vals,
                        len,
                        rhs,
                        conflictpool,
                        globaldom,
                        pseudocost,
                    );
                }
                let added_constraints = ms.mipdata_.conflictpools[conflictpool].get_num_conflicts()
                    != old_num_conflicts;

                if added_constraints {
                    localdomain.propagate(ms);
                    if localdomain.infeasible() {
                        return;
                    }

                    bound_changes.retain(|domchg| !localdomain.is_active(domchg));
                }

                if !bound_changes.is_empty() {
                    for domchg in bound_changes.iter() {
                        localdomain.change_bound(
                            ms,
                            domchg.clone(),
                            HighsDomain_Reason::unspecified(),
                        );
                        if localdomain.infeasible() {
                            break;
                        }
                    }

                    if !localdomain.infeasible() {
                        localdomain.propagate(ms);
                    }
                }
            } else {
                for domchg in bound_changes.iter() {
                    localdomain.change_bound(ms, domchg.clone(), HighsDomain_Reason::unspecified());
                    if localdomain.infeasible() {
                        break;
                    }
                }

                if !localdomain.infeasible() {
                    localdomain.propagate(ms);
                }
            }
        }
    }

    /// HiGHS: HighsRedcostFixing.cpp::HighsRedcostFixing::addRootRedcost (lines 195-324)
    /// C++: `void addRootRedcost(const HighsMipSolver & mipsolver, const std::vector<double> & lpredcost, double lpobjective)`
    // PORT NOTE: `lpredcost` aliases `getLp().getSolution().col_dual` at every C++ call site;
    // the caller clones it before detaching anything. The body runs
    // `getLp().computeBasicDegenerateDuals(..., getDomain(), getDomain(), getConflictPool(),
    // pseudocost, false)` by detaching `lps[0]`, `domains[0]` (passed as `localdom`; the
    // `globaldom` index 0 is then unused because `getdualproof == false`) and
    // `pseudocosts[0]`. `self` is detached from `ms.mipdata_.redcostfixing` by the caller.
    // HiGHS quirk: that call overwrites `col_dual` (the aliased `lpredcost`) and the C++ loops
    // below read the overwritten values. The body therefore re-reads `lps[0]`'s `col_dual`
    // after the call, and the caller's copy is not used.
    pub fn add_root_redcost(
        &mut self,
        ms: &mut HighsMipSolver,
        lpredcost: &[f64],
        lpobjective: f64,
    ) {
        let _ = lpredcost;
        let num_col = ms.num_col() as usize;
        self.lurking_col_lower
            .resize_with(num_col, CppMultimapF64::new);
        self.lurking_col_upper
            .resize_with(num_col, CppMultimapF64::new);

        // Provided domains won't be used (only used for dual proof)
        let feastol = ms.mipdata_.feastol;
        detach!(ms.mipdata_.lps[0], |lp| detach!(
            ms.mipdata_.domains[0],
            |d| detach!(ms.mipdata_.pseudocosts[0], |pc| lp
                .compute_basic_degenerate_duals(ms, feastol, d, 0, 0, pc, false))
        ));

        let lpredcost: Vec<f64> = ms.mipdata_.lps[0].get_solution().col_dual.clone();

        // Compute maximum number of steps per column with large domain
        // max_steps = 2 ** k, k = max(5, min(10 ,round(log(|D| / 10)))),
        // D = {col : integral_cols | (ub - lb) >= 512}
        // This is to avoid doing 2**10 steps when there's many unbounded columns
        let mut num_redcost_large_domain_cols: HighsInt = 0;
        for idx in 0..ms.mipdata_.integral_cols.len() {
            let ci = ms.mipdata_.integral_cols[idx] as usize;
            let domain = &ms.mipdata_.domains[0];
            if (domain.col_upper_[ci] - domain.col_lower_[ci]) >= 512.0
                && lpredcost[ci].abs() > ms.mipdata_.feastol
            {
                num_redcost_large_domain_cols += 1;
            }
        }
        let mut max_num_steps_exp: HighsInt = 10;
        let (_, mut expshift) = frexp(f64::from(num_redcost_large_domain_cols / 10));
        if expshift > 5 {
            expshift = cpp_min(expshift, max_num_steps_exp);
            max_num_steps_exp = max_num_steps_exp - expshift + 5;
        }
        let max_num_steps: HighsInt = (1u64 << max_num_steps_exp) as HighsInt;

        let feastol = ms.mipdata_.feastol;
        let lower_bound = ms.mipdata_.lower_bound;
        for idx in 0..ms.mipdata_.integral_cols.len() {
            let col = ms.mipdata_.integral_cols[idx];
            let ci = col as usize;
            if lpredcost[ci] > feastol {
                // col <= (cutoffbound - lpobj)/redcost + lb
                // so for lurkub = lb to ub - 1 we can compute the necessary cutoff
                // bound to reach this bound which is:
                //  lurkub = (cutoffbound - lpobj)/redcost + lb
                //  cutoffbound = (lurkub - lb) * redcost + lpobj
                let lb = ms.mipdata_.domains[0].col_lower_[ci];
                let ub = ms.mipdata_.domains[0].col_upper_[ci];
                let mut upper_entries = multimap_entries(&self.lurking_col_upper[ci]);
                let mut lower_entries = multimap_entries(&self.lurking_col_lower[ci]);
                // HiGHS #3474, #3481 (mip_lurking_bounds_robust): HighsRedcostFixing.cpp lines 310-318 at 3771ccf96
                let use_x1 = ms.options_mip_.mip_lurking_bounds_robust;
                if use_x1 {
                    find_lurking_bounds_x1(
                        feastol,
                        lower_bound,
                        1,
                        lb,
                        ub,
                        lpobjective,
                        lpredcost[ci],
                        max_num_steps,
                        max_num_steps_exp,
                        &mut upper_entries,
                        &mut lower_entries,
                    );
                } else {
                    find_lurking_bounds(
                        feastol,
                        lower_bound,
                        1,
                        cpp_double_to_int(lb),
                        cpp_double_to_int(ub),
                        ub != kHighsInf,
                        lpobjective,
                        lpredcost[ci],
                        max_num_steps,
                        max_num_steps_exp,
                        &mut upper_entries,
                        &mut lower_entries,
                    );
                }
                self.lurking_col_upper[ci] = multimap_from_entries(&upper_entries);
                self.lurking_col_lower[ci] = multimap_from_entries(&lower_entries);
            } else if lpredcost[ci] < -feastol {
                // col >= (cutoffbound - lpobj)/redcost + ub
                // so for lurklb = lb + 1 to ub we can compute the necessary cutoff
                // bound to reach this bound which is:
                //  lurklb = (cutoffbound - lpobj)/redcost + ub
                //  cutoffbound = (lurklb - ub) * redcost + lpobj
                let lb = ms.mipdata_.domains[0].col_lower_[ci];
                let ub = ms.mipdata_.domains[0].col_upper_[ci];
                let mut upper_entries = multimap_entries(&self.lurking_col_upper[ci]);
                let mut lower_entries = multimap_entries(&self.lurking_col_lower[ci]);
                // HiGHS #3474, #3481 (mip_lurking_bounds_robust): HighsRedcostFixing.cpp lines 319-327 at 3771ccf96
                let use_x1 = ms.options_mip_.mip_lurking_bounds_robust;
                if use_x1 {
                    find_lurking_bounds_x1(
                        feastol,
                        lower_bound,
                        -1,
                        ub,
                        lb,
                        lpobjective,
                        lpredcost[ci],
                        max_num_steps,
                        max_num_steps_exp,
                        &mut lower_entries,
                        &mut upper_entries,
                    );
                } else {
                    find_lurking_bounds(
                        feastol,
                        lower_bound,
                        -1,
                        cpp_double_to_int(ub),
                        cpp_double_to_int(lb),
                        lb != -kHighsInf,
                        lpobjective,
                        lpredcost[ci],
                        max_num_steps,
                        max_num_steps_exp,
                        &mut lower_entries,
                        &mut upper_entries,
                    );
                }
                self.lurking_col_upper[ci] = multimap_from_entries(&upper_entries);
                self.lurking_col_lower[ci] = multimap_from_entries(&lower_entries);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x1_rcf_infinite_bound_returns() {
        // input_bound = +inf, direction = 1: early return, both vecs stay empty
        let mut lurking: Vec<(f64, HighsInt)> = Vec::new();
        let mut other: Vec<(f64, HighsInt)> = Vec::new();
        find_lurking_bounds_x1(
            1e-6,
            -kHighsInf,
            1,
            kHighsInf,
            100.0,
            0.0,
            1.0,
            1024,
            10,
            &mut lurking,
            &mut other,
        );
        assert!(lurking.is_empty());
        assert!(other.is_empty());

        // input_bound = -inf, direction = -1: early return, both vecs stay empty
        let mut lurking: Vec<(f64, HighsInt)> = Vec::new();
        let mut other: Vec<(f64, HighsInt)> = Vec::new();
        find_lurking_bounds_x1(
            1e-6,
            -kHighsInf,
            -1,
            -kHighsInf,
            -100.0,
            0.0,
            -1.0,
            1024,
            10,
            &mut lurking,
            &mut other,
        );
        assert!(lurking.is_empty());
        assert!(other.is_empty());
    }

    #[test]
    fn x1_rcf_matches_old_finite() {
        // Simple LCG for deterministic pseudo-random bounds
        let mut state: u64 = 123456789;
        let mut next = || -> u64 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state
        };

        let red_costs: [f64; 4] = [0.5, 1.0, -2.5, 7.0];
        let lp_objectives: [f64; 2] = [0.0, 3.0];

        for _case in 0..200 {
            let input_bound: f64 = ((next() % 10001) as f64) - 5000.0;
            let input_other_bound: f64 = ((next() % 10001) as f64) - 5000.0;
            let direction: HighsInt = if next() % 2 == 0 { 1 } else { -1 };
            let red_cost = red_costs[(next() % 4) as usize];
            let lp_objective = lp_objectives[(next() % 2) as usize];

            let feastol = 1e-6;
            let lower_bound = -kHighsInf;
            let max_num_steps: HighsInt = 1024;
            let max_num_steps_exp: HighsInt = 10;

            // Old function: is_other_bound_finite = true, other as HighsInt
            let mut old_lurking: Vec<(f64, HighsInt)> = Vec::new();
            let mut old_other: Vec<(f64, HighsInt)> = Vec::new();
            find_lurking_bounds(
                feastol,
                lower_bound,
                direction,
                input_bound as HighsInt,
                input_other_bound as HighsInt,
                true,
                lp_objective,
                red_cost,
                max_num_steps,
                max_num_steps_exp,
                &mut old_lurking,
                &mut old_other,
            );

            // New function: f64 inputs
            let mut new_lurking: Vec<(f64, HighsInt)> = Vec::new();
            let mut new_other: Vec<(f64, HighsInt)> = Vec::new();
            find_lurking_bounds_x1(
                feastol,
                lower_bound,
                direction,
                input_bound,
                input_other_bound,
                lp_objective,
                red_cost,
                max_num_steps,
                max_num_steps_exp,
                &mut new_lurking,
                &mut new_other,
            );

            assert_eq!(
                old_lurking, new_lurking,
                "case {_case}: input_bound={input_bound}, input_other_bound={input_other_bound}, \
                 direction={direction}, red_cost={red_cost}, lp_objective={lp_objective}"
            );
            assert_eq!(
                old_other, new_other,
                "case {_case}: input_bound={input_bound}, input_other_bound={input_other_bound}, \
                 direction={direction}, red_cost={red_cost}, lp_objective={lp_objective}"
            );
        }
    }

    /// arm64 `fcvtzs` (the oracle build): saturating, NaN -> 0 (clang -O2 probe).
    #[test]
    fn double_to_int_saturates_like_arm64() {
        assert_eq!(cpp_double_to_int(kHighsInf), HighsInt::MAX);
        assert_eq!(cpp_double_to_int(-kHighsInf), HighsInt::MIN);
        assert_eq!(cpp_double_to_int(f64::NAN), 0);
        assert_eq!(cpp_double_to_int(1e10), HighsInt::MAX);
        assert_eq!(cpp_double_to_int(-1e10), HighsInt::MIN);
        assert_eq!(cpp_double_to_int(-3.7), -3);
    }

    /// libc++ `multimap::emplace_hint(lower_bound(k), k, v)` links the new node before the
    /// entries with key `k` (probe: {2:20, 2:21} + hint insert of 2:99 iterates 99, 20, 21).
    #[test]
    fn emplace_hint_places_before_equal_keys() {
        let feastol = 1e-6;
        let key = (0.0 + (1.0 - 10.0 * feastol)) * 1.0 + 0.0;
        let mut lurking = vec![(key, 5)];
        let mut other = Vec::new();
        // direction 1, bound 0, other bound 1 (finite): one candidate, lurking bound 0 at `key`.
        find_lurking_bounds(
            feastol,
            -kHighsInf,
            1,
            0,
            1,
            true,
            0.0,
            1.0,
            1024,
            10,
            &mut lurking,
            &mut other,
        );
        assert_eq!(lurking, vec![(key, 0), (key, 5)]);
        let map = multimap_from_entries(&lurking);
        assert_eq!(multimap_entries(&map), lurking);
    }
}
